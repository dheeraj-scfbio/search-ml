import pandas as pd
import numpy as np
import joblib
import re
import os
from tqdm import tqdm

def clean_feature_name(name):
    """
    Cleans feature names to be compatible with models like XGBoost.
    """
    name = re.sub(r'[,\(\)\-\s]', '_', name)
    name = re.sub(r'[\{\}\[\]:"]', '', name)
    return name

def parse_parameter_file(file_path):
    """
    Parses parameter file (Key: Min,Max).
    """
    params = {}
    if not os.path.exists(file_path):
        print(f"Warning: Parameter file '{file_path}' not found.")
        return None

    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or ':' not in line:
                    continue
                try:
                    key, value = line.split(':')
                    key = key.strip()
                    if ',' in value:
                        min_val, max_val = value.split(',')
                        params[key] = (float(min_val.strip()), float(max_val.strip()))
                    else:
                        val = float(value.strip())
                        params[key] = (val, val)
                except ValueError:
                    pass
    except Exception as e:
        print(f"Error reading parameter file: {e}")
        return None
    return params

def perform_virtual_screening(molecule_csv_path, protein_csv_path, models_dir, output_csv_path, parameter_file_path):
    print("--- Starting Virtual Screening ---")

    # --- 1. Load Models and Scaler ---
    print(f"Loading models from '{models_dir}'...")
    try:
        scaler = joblib.load(os.path.join(models_dir, 'scaler_final.joblib'))
        final_xgb = joblib.load(os.path.join(models_dir, 'final_xgb_model.joblib'))
        final_cat = joblib.load(os.path.join(models_dir, 'final_cat_model.joblib'))
        final_rf = joblib.load(os.path.join(models_dir, 'final_rf_model.joblib'))
        final_lgbm = joblib.load(os.path.join(models_dir, 'final_lgbm_model.joblib'))
        final_meta_model = joblib.load(os.path.join(models_dir, 'final_meta_model.joblib'))
        
        order_path = os.path.join(models_dir, 'final_ensemble_model_order.txt')
        if os.path.exists(order_path):
            with open(order_path, 'r') as f:
                model_order = f.read().splitlines()
        else:
            model_order = ['xgb', 'cat', 'rf', 'lgbm']

        print("All models loaded successfully.")
    except FileNotFoundError as e:
        print(f"Error loading models: {e}")
        return

    # --- 2. Load Protein Features ---
    print(f"Loading protein features from '{protein_csv_path}'...")
    try:
        protein_features_df = pd.read_csv(protein_csv_path)
        protein_features_df = protein_features_df.drop(columns=[col for col in protein_features_df.columns if 'Unnamed' in col], errors='ignore')
        if len(protein_features_df) != 1:
            print(f"Warning: Protein CSV has {len(protein_features_df)} rows. Using the first one.")
        protein_features = protein_features_df.iloc[[0]]
    except FileNotFoundError:
        print(f"Error: Protein file not found.")
        return

    # --- 3. Stream Molecules and Predict ---
    print(f"Streaming molecules from '{molecule_csv_path}'...")
    
    try:
        if hasattr(scaler, 'feature_names_in_'):
            feature_columns_ordered = list(scaler.feature_names_in_)
        else:
            feature_columns_ordered = None
    except AttributeError:
        feature_columns_ordered = None

    all_results = []
    chunksize = 50000 
    
    try:
        print("Processing chunks...")
        for molecule_chunk in pd.read_csv(molecule_csv_path, chunksize=chunksize):
            molecule_chunk = molecule_chunk.drop(columns=[col for col in molecule_chunk.columns if 'Unnamed' in col], errors='ignore')

            # Combine Protein + Ligand
            protein_repeated = pd.concat([protein_features] * len(molecule_chunk), ignore_index=True)
            combined_chunk = pd.concat([molecule_chunk.reset_index(drop=True), protein_repeated], axis=1)

            # Feature Prep
            X_screening = combined_chunk.copy()
            X_screening.columns = [clean_feature_name(col) for col in X_screening.columns]

            if feature_columns_ordered:
                missing = set(feature_columns_ordered) - set(X_screening.columns)
                if missing:
                    print(f"Error: Missing columns: {list(missing)[:5]}")
                    return
                X_screening = X_screening[feature_columns_ordered]
            else:
                X_screening = X_screening.select_dtypes(include=[np.number])

            # Scale
            X_screening_scaled = pd.DataFrame(scaler.transform(X_screening), columns=X_screening.columns, index=X_screening.index)

            # --- PREDICTIONS ---
            # 1. Base Model Predictions (pKi units)
            base_preds = {}
            base_preds['xgb'] = final_xgb.predict(X_screening_scaled)
            base_preds['cat'] = final_cat.predict(X_screening_scaled)
            base_preds['rf'] = final_rf.predict(X_screening_scaled)
            base_preds['lgbm'] = final_lgbm.predict(X_screening_scaled)

            # 2. Stack for Meta-Model
            meta_features = []
            for model_name in model_order:
                if model_name in base_preds:
                    meta_features.append(base_preds[model_name])
            
            meta_X_screening = np.column_stack(meta_features)

            # 3. Final Meta-Model Prediction (pKi units)
            final_pred_pKi = final_meta_model.predict(meta_X_screening)

            # --- STORE RESULTS ---
            chunk_results = molecule_chunk.copy()
            
            # A. Meta-Model Results (All Units)
            chunk_results['PBE (pKi/pKd/pIC50)'] = np.round(final_pred_pKi, 2)
            chunk_results['PBE (kcal/mol)'] = np.round(-1.364 * final_pred_pKi, 2)
            chunk_results['PBE (nM)'] = np.round(10**(9 - final_pred_pKi), 2)
            
            # B. Individual Models (Energy Units ONLY)
            chunk_results['XGB (kcal/mol)'] = np.round(-1.364 * base_preds['xgb'], 2)
            chunk_results['CAT (kcal/mol)'] = np.round(-1.364 * base_preds['cat'], 2)
            chunk_results['RF (kcal/mol)'] = np.round(-1.364 * base_preds['rf'], 2)
            chunk_results['LGBM (kcal/mol)'] = np.round(-1.364 * base_preds['lgbm'], 2)

            all_results.append(chunk_results)
            print(f"Processed {len(molecule_chunk)} molecules...")

    except FileNotFoundError:
        print(f"Error: Molecule file '{molecule_csv_path}' not found.")
        return

    # --- 4. Consolidate and Filter ---
    print("Consolidating results...")
    if not all_results:
        print("No results generated.")
        return

    results_with_predictions = pd.concat(all_results, ignore_index=True)
    
    # Sort by Meta-Model Energy (More negative is better)
    results_with_predictions.sort_values(by='PBE (kcal/mol)', ascending=True, inplace=True)

    # Filter (Using Parameters logic)
    params = parse_parameter_file(parameter_file_path)
    filtered_results = results_with_predictions
    
    if params:
        print(f"Applying filters...")
        initial_count = len(filtered_results)
        for param, (min_val, max_val) in params.items():
            if param in filtered_results.columns:
                filtered_results = filtered_results[
                    (filtered_results[param] >= min_val) & 
                    (filtered_results[param] <= max_val)
                ]
        print(f"Filtered: {len(filtered_results)} / {initial_count} molecules remaining.")

    # --- 5. Final Column Organization (STRICT) ---
    cols = list(filtered_results.columns)
    
    # Identify ID column (PDB_File or similar)
    id_candidates = ['PDB_File', 'Molecule', 'Name', 'ID', 'ZINC_ID', 'SMILES']
    id_col = next((c for c in id_candidates if c in cols), cols[0])
    
    # Define EXACT requested order
    meta_cols = ['PBE (pKi/pKd/pIC50)', 'PBE (kcal/mol)', 'PBE (nM)']
    individual_cols = ['XGB (kcal/mol)', 'CAT (kcal/mol)', 'RF (kcal/mol)', 'LGBM (kcal/mol)']
    
    # Construct final list - EXCLUDING all features
    final_order = [id_col] + meta_cols + individual_cols
    
    # Validation to ensure columns exist
    final_order = [c for c in final_order if c in filtered_results.columns]
    
    results_to_save = filtered_results[final_order].copy()
    
    results_to_save.to_csv(output_csv_path, index=False)
    print(f"\n--- Completed ---")
    print(f"Results saved to '{output_csv_path}'")

if __name__ == '__main__':
    # --- Configuration ---
    MODELS_DIRECTORY = './'
    MOLECULE_CSV = 'dataset.csv' 
    PROTEIN_CSV = 'protein.csv'   
    OUTPUT_CSV = 'results.csv'
    PARAMETER_FILE = 'parameter.txt'

    perform_virtual_screening(
        molecule_csv_path=MOLECULE_CSV,
        protein_csv_path=PROTEIN_CSV,
        models_dir=MODELS_DIRECTORY,
        output_csv_path=OUTPUT_CSV,
        parameter_file_path=PARAMETER_FILE
    )
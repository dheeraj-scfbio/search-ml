import pandas as pd
import numpy as np
import joblib
import re
import os
from tqdm import tqdm

def clean_feature_name(name):
    """Cleans feature names to be compatible with models like XGBoost."""
    name = re.sub(r'[,\(\)\-\s]', '_', name)
    name = re.sub(r'[\{\}\[\]:"]', '', name)
    return name

def perform_reverse_screening(molecule_csv_path, protein_csv_path, models_dir, output_csv_path):
    print("--- Starting Reverse Virtual Screening (1 Ligand vs Many Proteins) ---")

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
        model_order = open(order_path, 'r').read().splitlines() if os.path.exists(order_path) else ['xgb', 'cat', 'rf', 'lgbm']
        
        print("All models loaded successfully.")
    except Exception as e:
        print(f"Error loading models: {e}")
        return

    # --- 2. Load Single Molecule Features ---
    print(f"Loading ligand features from '{molecule_csv_path}'...")
    try:
        mol_df = pd.read_csv(molecule_csv_path)
        mol_df = mol_df.drop(columns=[col for col in mol_df.columns if 'Unnamed' in col], errors='ignore')
        if len(mol_df) != 1:
            print(f"Warning: Molecule CSV has {len(mol_df)} rows. Using the first one.")
        molecule_features = mol_df.iloc[[0]]
    except Exception as e:
        print(f"Error loading ligand file: {e}")
        return

    # --- 3. Stream Proteins and Predict ---
    print(f"Streaming proteins from '{protein_csv_path}'...")
    
    # Get expected feature order from scaler
    feature_columns_ordered = list(scaler.feature_names_in_) if hasattr(scaler, 'feature_names_in_') else None

    all_results = []
    chunksize = 10000 
    
    try:
        for protein_chunk in pd.read_csv(protein_csv_path, chunksize=chunksize):
            protein_chunk = protein_chunk.drop(columns=[col for col in protein_chunk.columns if 'Unnamed' in col], errors='ignore')

            # REVERSE COMBINATION: Combine 1 Molecule + Many Proteins
            molecule_repeated = pd.concat([molecule_features] * len(protein_chunk), ignore_index=True)
            # Ensure indexes align for horizontal concat
            combined_chunk = pd.concat([molecule_repeated.reset_index(drop=True), protein_chunk.reset_index(drop=True)], axis=1)

            # Feature Prep
            X_screening = combined_chunk.copy()
            X_screening.columns = [clean_feature_name(col) for col in X_screening.columns]

            if feature_columns_ordered:
                X_screening = X_screening[feature_columns_ordered]
            else:
                X_screening = X_screening.select_dtypes(include=[np.number])

            # Scale
            X_screening_scaled = pd.DataFrame(scaler.transform(X_screening), columns=X_screening.columns)

            # --- PREDICTIONS ---
            base_preds = {
                'xgb': final_xgb.predict(X_screening_scaled),
                'cat': final_cat.predict(X_screening_scaled),
                'rf': final_rf.predict(X_screening_scaled),
                'lgbm': final_lgbm.predict(X_screening_scaled)
            }

            # Stack for Meta-Model
            meta_X_screening = np.column_stack([base_preds[m] for m in model_order if m in base_preds])
            final_pred_pKi = final_meta_model.predict(meta_X_screening)

            # --- STORE RESULTS ---
            chunk_results = protein_chunk.copy()
            
            # Binding Energy Calculations
            chunk_results['PBE (pKi/pKd/pIC50)'] = np.round(final_pred_pKi, 2)
            chunk_results['PBE (kcal/mol)'] = np.round(-1.364 * final_pred_pKi, 2)
            chunk_results['PBE (nM)'] = np.round(10**(9 - final_pred_pKi), 2)
            
            # Individual Model contributions
            chunk_results['XGB (kcal/mol)'] = np.round(-1.364 * base_preds['xgb'], 2)
            chunk_results['CAT (kcal/mol)'] = np.round(-1.364 * base_preds['cat'], 2)
            chunk_results['RF (kcal/mol)'] = np.round(-1.364 * base_preds['rf'], 2)
            chunk_results['LGBM (kcal/mol)'] = np.round(-1.364 * base_preds['lgbm'], 2)

            all_results.append(chunk_results)
            print(f"Processed {len(protein_chunk)} protein entries...")

    except Exception as e:
        print(f"Error during processing: {e}")
        return

    # --- 4. Consolidate and Save ---
    if not all_results:
        print("No results generated.")
        return

    final_df = pd.concat(all_results, ignore_index=True)
    
    # Sort by Meta-Model Energy (Best binders first)
    final_df.sort_values(by='PBE (kcal/mol)', ascending=True, inplace=True)

    # Column Selection: Identify Protein ID column
    cols = list(final_df.columns)
    id_candidates = ['Protein_ID', 'Entry', 'Accession', 'Protein_Name', 'PDB_ID']
    id_col = next((c for c in id_candidates if c in cols), cols[0])
    
    meta_cols = ['PBE (pKi/pKd/pIC50)', 'PBE (kcal/mol)', 'PBE (nM)']
    indiv_cols = ['XGB (kcal/mol)', 'CAT (kcal/mol)', 'RF (kcal/mol)', 'LGBM (kcal/mol)']
    
    output_order = [id_col] + meta_cols + indiv_cols
    output_order = [c for c in output_order if c in final_df.columns]
    
    final_df[output_order].to_csv(output_csv_path, index=False)
    print(f"\n--- Completed ---\nResults saved to '{output_csv_path}'")

if __name__ == '__main__':
    MODELS_DIRECTORY = './'
    MOLECULE_CSV = 'dataset.csv' # One ligand
    PROTEIN_CSV = 'protein.csv'  # Multiple proteins
    OUTPUT_CSV = 'reverse_screening_results.csv'

    perform_reverse_screening(
        molecule_csv_path=MOLECULE_CSV,
        protein_csv_path=PROTEIN_CSV,
        models_dir=MODELS_DIRECTORY,
        output_csv_path=OUTPUT_CSV
    )

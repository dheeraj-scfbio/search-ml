import pandas as pd
import os

def merge_feature_csvs(protein_csv, pocket_csv, output_csv):
    """
    Merges ligand, protein, and volume feature CSV files into a single file.
    Replicates single-row protein and volume data to match ligand rows.
    
    Args:
        ligand_csv (str): Path to the ligand features CSV file (multiple rows).
        protein_csv (str): Path to the protein features CSV file (single row).
        volume_csv (str): Path to the pocket volume CSV file (single row).
        output_csv (str): Path for the final merged output CSV file.
    """

    try:
        df_protein = pd.read_csv(protein_csv)
        print(f"Loaded '{protein_csv}' with shape: {df_protein.shape}")
    except FileNotFoundError:
        print(f"Error: The file '{protein_csv}' was not found.")
        return

    try:
        df_pocket = pd.read_csv(pocket_csv)
        print(f"Loaded '{pocket_csv}' with shape: {df_pocket.shape}")
    except FileNotFoundError:
        print(f"Error: The file '{pocket_csv}' was not found.")
        return

    # --- Check and replicate protein and volume rows ---
    if len(df_protein) != 1:
        print("Error: Protein feature CSV should contain exactly one row.")
        return
    if len(df_pocket) != 1:
        print("Error: Pocket feature CSV should contain exactly one row.")
        return


    print(f"Merging Protein features on 'PDB_File'...")
    df_merged = pd.merge(df_protein, df_pocket, on='PDB_File', how='inner')
    print(f"Shape after merge: {df_merged.shape}")

    # --- 4. Save the final merged file ---
    df_merged.to_csv(output_csv, index=False)
    print(f"\n--- Success! ---")
    print(f"Successfully merged all features into '{output_csv}'")


if __name__ == "__main__":
    protein_features_file = "protein_properties.csv"
    pocket_features_file = "pocket_properties.csv"
    final_output_file = "protein.csv"

    merge_feature_csvs(
        protein_csv=protein_features_file,
        pocket_csv=pocket_features_file,
        output_csv=final_output_file
    )


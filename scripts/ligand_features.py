import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Descriptors3D, Crippen, AllChem, rdMolDescriptors
import os
import argparse
import csv # Use csv module for memory-efficient writing
from tqdm import tqdm # Import tqdm for a progress bar

# Define the exact column order
DESIRED_COLUMN_ORDER = [
    "PDB_File", "LigVol", "MolWt", "MolLogP", "TPSA", "NumHDonors",
    "NumHAcceptors", "NumRotatableBonds", "NumAromaticRings",
    "FractionCSP3", "NumHeteroatoms", "HeavyAtomCount", "NumAtoms",
    "MolMR", "BertzCT", "NumCarbons", "Asphericity", "Eccentricity",
    "InertialShapeFactor", "NPR1", "NPR2", "PMI1", "PMI2", "PMI3",
    "RadiusOfGyration", "SpherocityIndex"
]

# --- Pre-define descriptor functions ---
DESCRIPTOR_FUNCS_2D = {
    "MolWt": Descriptors.MolWt,
    "MolLogP": Descriptors.MolLogP,
    "TPSA": Descriptors.TPSA,
    "NumHDonors": Descriptors.NumHDonors,
    "NumHAcceptors": Descriptors.NumHAcceptors,
    "NumRotatableBonds": Descriptors.NumRotatableBonds,
    "NumAromaticRings": Descriptors.NumAromaticRings,
    "FractionCSP3": Descriptors.FractionCSP3,
    "NumHeteroatoms": Descriptors.NumHeteroatoms,
    "BertzCT": Descriptors.BertzCT
}

DESCRIPTOR_FUNCS_3D = {
    "LigVol": AllChem.ComputeMolVolume,
    "Asphericity": Descriptors3D.Asphericity,
    "Eccentricity": Descriptors3D.Eccentricity,
    "InertialShapeFactor": Descriptors3D.InertialShapeFactor,
    "NPR1": Descriptors3D.NPR1,
    "NPR2": Descriptors3D.NPR2,
    "PMI1": Descriptors3D.PMI1,
    "PMI2": Descriptors3D.PMI2,
    "PMI3": Descriptors3D.PMI3,
    "RadiusOfGyration": Descriptors3D.RadiusOfGyration,
    "SpherocityIndex": Descriptors3D.SpherocityIndex
}

def count_carbons(mol):
    """Counts the number of carbon atoms in a molecule."""
    if mol is None:
        return 0
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

def calculate_features_sequentially(sdf_file_path, output_csv_file, batch_size=5000):
    """
    Calculates features sequentially and writes to CSV in batches.
    
    Args:
        sdf_file_path (str): Path to the input SDF file.
        output_csv_file (str): Path to the output CSV file.
        batch_size (int): Number of molecules to process before writing to disk.
    """
    if not os.path.exists(sdf_file_path):
        print(f"Error: Input SDF file not found: {sdf_file_path}")
        return

    print(f"Opening SDF file: {sdf_file_path}")
    # We must disable sanitization in the supplier to handle errors manually
    supplier = Chem.SDMolSupplier(sdf_file_path, sanitize=False, removeHs=False)
    if not supplier:
        print("Error: Could not open SDMolSupplier. Check file format.")
        return

    # --- Setup Batch Writing ---
    batch = []
    total_processed = 0
    total_failed = 0
    
    print(f"Processing molecules sequentially, writing in batches of {batch_size}...")

    try:
        with open(output_csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=DESIRED_COLUMN_ORDER, extrasaction='ignore')
            writer.writeheader()

            # Wrap the supplier with tqdm for a progress bar
            for i, mol in enumerate(tqdm(supplier, desc="Processing Molecules")):
                
                # --- Robust Molecule Processing ---
                try:
                    if mol is None:
                        print(f"Warning: Molecule at index {i} is None. Skipping.")
                        total_failed += 1
                        continue
                        
                    # --- Manual Sanitization ---
                    # This is where errors often happen.
                    # We try to sanitize, and if it fails, we skip the molecule.
                    Chem.SanitizeMol(mol)

                    # --- Get Molecule Name ---
                    if mol.HasProp('IDENTIFIER'):
                        mol_name = mol.GetProp('IDENTIFIER')
                    elif mol.HasProp('_Name'):
                        mol_name = mol.GetProp('_Name')
                    else:
                        mol_name = f"Molecule_{i+1}"
                    
                    current_features = {"PDB_File": mol_name}

                    # --- Calculate 2D Descriptors ---
                    for name, func in DESCRIPTOR_FUNCS_2D.items():
                        current_features[name] = func(mol)
                        
                    # --- Calculate 3D Descriptors ---
                    for name, func in DESCRIPTOR_FUNCS_3D.items():
                        current_features[name] = func(mol)
                        
                    # --- Calculate other properties ---
                    current_features["HeavyAtomCount"] = mol.GetNumHeavyAtoms()
                    current_features["NumAtoms"] = mol.GetNumAtoms()
                    current_features["MolMR"] = Crippen.MolMR(mol)
                    current_features["NumCarbons"] = count_carbons(mol)
                    
                    # Add the valid features to our batch
                    batch.append(current_features)
                    total_processed += 1

                except Exception as e:
                    # Catch any error during sanitization or feature calculation
                    print(f"Warning: Failed to process molecule at index {i}. Error: {e}. Skipping.")
                    total_failed += 1
                    continue
                # --- End of Robust Processing Block ---

                # --- Write batch to disk ---
                if len(batch) >= batch_size:
                    writer.writerows(batch)
                    batch = [] # Clear the batch to save RAM
            
            # Write any remaining molecules in the last batch
            if batch:
                writer.writerows(batch)
                
    except Exception as e:
        print(f"\nAn critical error occurred: {e}")
        print("The process may be incomplete.")
    finally:
        print("\n--- Feature Calculation Complete ---")
        print(f"Total valid molecules processed and saved: {total_processed}")
        print(f"Total molecules that failed to process: {total_failed}")
        print(f"Features saved to '{output_csv_file}'")


def main():
    parser = argparse.ArgumentParser(description="Calculate 2D and 3D molecular descriptors from an SDF file sequentially.")
    
    parser.add_argument("-i", "--input", 
                        required=True, 
                        help="Path to the input multi-ligand SDF file.")
    
    parser.add_argument("-o", "--output", 
                        required=True, 
                        help="Path to the output CSV file to save features.")
    
    parser.add_argument("-b", "--batchsize", 
                        type=int, 
                        default=5000, 
                        help="Number of molecules to write per batch. (Default: 5000)")

    args = parser.parse_args()

    # Pass arguments to the main function
    calculate_features_sequentially(args.input, args.output, args.batchsize)

# --- Main execution block ---
if __name__ == "__main__":
    main()

import os
import glob
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem

def find_input_file(base_name):
    """Searches for a file with the given base name and supported extensions."""
    extensions = ['.sdf', '.pdb', '.smi', '.smiles']
    for ext in extensions:
        # Search for exact name + extension
        files = glob.glob(f"{base_name}{ext}")
        if files:
            return files[0], ext
    return None, None

def convert_to_3d_sdf(base_input_name, output_file):
    # 1. Detect the file automatically
    input_path, ext = find_input_file(base_input_name)
    
    if not input_path:
        print(f"[!] Error: No file found matching '{base_input_name}' with supported extensions (.sdf, .pdb, .smi)")
        return

    print(f"[*] Detected file: {input_path} (Format: {ext})")
    mol = None

    # 2. Load the molecule based on detected format
    try:
        if ext in ['.smi', '.smiles']:
            with open(input_path, 'r') as f:
                line = f.readline().split()
                if line:
                    mol = Chem.MolFromSmiles(line[0])
        elif ext == '.pdb':
            mol = Chem.MolFromPDBFile(input_path, removeHs=False)
        elif ext == '.sdf':
            suppl = Chem.SDMolSupplier(input_path, removeHs=False)
            mol = next(suppl) if suppl else None
    except Exception as e:
        print(f"[!] parsing error: {e}")
        return

    if mol is None:
        print("[!] Error: Could not parse molecule. Check input file validity.")
        return

    # 3. Add Hydrogens (Necessary for proper 3D geometry)
    mol = Chem.AddHs(mol)

    # 4. Generate 3D Coordinates using ETKDGv3
    print("[*] Generating 3D coordinates and embedding...")
    params = AllChem.ETKDGv3()
    if AllChem.EmbedMolecule(mol, params) == -1:
        print("[!] Standard 3D embedding failed. Trying with random coordinates...")
        AllChem.EmbedMolecule(mol, useRandomCoords=True)

    # 5. Energy Minimization (MMFF94 Force Field)
    print("[*] Optimizing geometry with MMFF94 force field...")
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except:
        print("[!] MMFF failed. Falling back to UFF optimization...")
        AllChem.UFFOptimizeMolecule(mol)

    # 6. Save output as SDF (Always ensures .sdf extension)
    if not output_file.lower().endswith('.sdf'):
        output_file += ".sdf"

    writer = Chem.SDWriter(output_file)
    writer.write(mol)
    writer.close()

    print(f"[+] Success! 3D structure saved as: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Auto-detect format and convert to 3D SDF.")
    parser.add_argument("input_name", help="The name of the file (without extension)")
    parser.add_argument("output_name", help="The name of the output file")

    args = parser.parse_args()
    convert_to_3d_sdf(args.input_name, args.output_name)
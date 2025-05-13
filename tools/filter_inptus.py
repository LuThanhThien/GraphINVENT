import argparse
from pathlib import Path
from rdkit import Chem

total_found = 0
total = 0
num_safe_backup = 10_000

motif_list = [
    "C1=CC=C(C=C1)C2=CC(=O)C3=CC=CC=C3O2",
    # "C1=CC(=C(C=C1O)C2=CC(=O)C3=C(O2)C=CC=C3O)",
    # "C1=CC=C(C=C1)C2=CC(=O)C3=C(O2)C=CC=C3O",
    # "C1=CC=C(C=C1)C2=CC(=O)C3=C(C2=O)C=CC=C3",
    # "C1=CC=C2C(=C1)C(=O)C3=C(O2)C=CC=C3O",
]

def passes_filters(mol, motif_list):
    # 1) does it match any of the motifs?
    has_motif = any(
        mol.HasSubstructMatch(Chem.MolFromSmiles(smiles))
        for smiles in motif_list
    )
    if not has_motif:
        return False
    if mol.GetNumHeavyAtoms() >= 50:
        return False
    return True

    # 2) does it have a pIC50 property?
    if not mol.HasProp('pIC50'):
        return False

    # 3) pull it out as a float, safely
    try:
        pIC50 = float(mol.GetProp('pIC50'))
        print("pIC50:", pIC50)
    except ValueError:
        # in case the prop is present but not parseable
        return False

    return pIC50 > 5

def save_safe_backup(output_smi, list_valid):
    global total
    global num_safe_backup
    if total % num_safe_backup == 0 and len(list_valid) > 0:
        # append the current list to the output file
        print(f"Saving backup to {output_smi}...")
        with open(output_smi, "a") as fout:
            fout.writelines(list_valid)
        list_valid.clear()  # clear the list after saving
    return list_valid

def filter_smi(input_smi, output_smi, motif_list, get_id=0):
    global total_found
    global total
    num_passed = 0
    if not input_smi.exists():
        print(f"File {input_smi} does not exist.")
        return
    if not output_smi.parent.exists():
        print(f"Creating output directory {output_smi.parent}.")
        output_smi.parent.mkdir(parents=True, exist_ok=True)

    list_valid = []
    with open(input_smi) as fin:
        for i, line in enumerate(fin):
            if i % 1000 == 0:
                print(f"Reading line: {i}, found: {num_passed} molecules.")

            if i == 0:
                continue
            
            list_valid = save_safe_backup(output_smi, list_valid)
            
            total += 1
            parts = line.strip().split()
            if len(parts) < get_id + 1:
                print(f"Skipping line {i}: {line.strip()}")
                continue
            
            smi = parts[get_id]
            # print("SMILES:", smi)

            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue

            if passes_filters(mol, motif_list):
                list_valid.append(f"{smi} {num_passed}\n")
                num_passed += 1

    if num_passed > 0:
        with open(output_smi, "a") as fout:
            fout.writelines(list_valid)
        total_found += num_passed
        print(f"Filtering {input_smi} -> {output_smi}")
        print(f"Filtered {input_smi}: {num_passed} molecules found.")

def main():
    p = argparse.ArgumentParser(
        description="Filter a .smi file by QED, SA, size, benzothiazole motif, etc."
    )
    p.add_argument("input",  help="Folder containing .smi files")
    p.add_argument("--is_file", "-f", action="store_true", help="Input is a file, not a folder")
    p.add_argument("--extension", "-e", default=".msmi", help="File extension to filter")
    p.add_argument("--get_id", "-g", type=int, default=0, help="Column index for SMILES in the input file (default: 0)")
    args = p.parse_args()

    # benzothiazole SMARTS (aromatic form):
    if args.is_file:
        input_file = Path(args.input)
        output_file = input_file.parent / "filtered" / input_file.name.replace(args.extension, ".filtered.smi")
        filter_smi(input_file, output_file, motif_list, get_id=args.get_id)
    else:
        input_folder = Path(args.input)
        output_folder = input_folder.parent / "filtered"
        for smi_file in input_folder.glob(f"*{args.extension}"):
            output_file = output_folder / smi_file.name.replace(args.extension, ".filtered.smi")
            filter_smi(smi_file, output_file, motif_list, get_id=args.get_id)
    print(f"Done. Total found: {total_found}/{total} molecules.")

if __name__=="__main__":
    main()

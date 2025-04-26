import argparse
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import QED

# You’ll need the “sascorer.py” from the RDKit contribs (or from:
# https://gist.github.com/delokal/ff4e60f85b375ddfc1b),
# placed alongside this script:
from sascorer import calculateScore

total_found = 0
total = 0

def passes_filters(mol, qed_thresh, sa_thresh, max_atoms, motif):
    # 1. size
    if mol.GetNumAtoms() > max_atoms:
        return False
    # 2. drug‐likeness
    # if QED.qed(mol) < qed_thresh:
    #     return False
    # if calculateScore(mol) > sa_thresh:
    #     return False
    # 3. substructure
    if not mol.HasSubstructMatch(motif):
        return False
    return True

def filter_smi(input_smi, output_smi, motif, qed_thresh=0.5, sa_thresh=5.0, max_atoms=50):
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
            if i == 0:
                continue
            
            total += 1
            # remove your “Xe” placeholders and empty‐ID lines up front
            if "[Xe]" in line or line.strip().isdigit():
                continue

            parts = line.strip().split()
            if len(parts) < 2:
                continue
            smi, name = parts[0], parts[1]

            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue

            if passes_filters(mol, qed_thresh, sa_thresh, max_atoms, motif):
                # write out SMILES <tab > Name/ID
                list_valid.append(f"{smi} {name}\n")
                # fout.write(f"{smi}\t{name}\n")
                num_passed += 1
                
    if num_passed > 0:
        with open(output_smi, "w") as fout:
            fout.writelines(list_valid)
        total_found += num_passed
        print(f"Filtering {input_smi} -> {output_smi}")
        print(f"Filtered {input_smi}: {num_passed} molecules found.")

def main():
    p = argparse.ArgumentParser(
        description="Filter a .smi file by QED, SA, size, benzothiazole motif, etc."
    )
    p.add_argument("input",  help="Folder containing .smi files")
    p.add_argument("--qed",       type=float, default=0.5, help="minimum QED")
    p.add_argument("--sa",        type=float, default=5.0, help="maximum SA score")
    p.add_argument("--max-atoms", type=int,   default=50,  help="max # of atoms")
    args = p.parse_args()

    # benzothiazole SMARTS (aromatic form):
    motif = Chem.MolFromSmarts("c1ccc2ncsc2c1")
    input_folder = Path(args.input)
    output_folder = input_folder.parent / "filtered"
    for smi_file in input_folder.glob("*.smi"):
        output_file = output_folder / smi_file.name
        filter_smi(smi_file, output_file, motif, args.qed, args.sa, args.max_atoms)
    print(f"Done. Total found: {total_found}/{total} molecules.")
    
if __name__ == "__main__":
    main()
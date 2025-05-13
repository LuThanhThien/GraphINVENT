import math
import random
import rdkit
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.Chem.rdmolfiles import SmilesMolSupplier
from argparse import ArgumentParser

parser = ArgumentParser(description="Visualize molecules from a SMILES file.")
parser.add_argument("smi_file", type=str, help="Path to the SMILES file.")
parser.add_argument("-n", "--n_samples", type=int, default=30, help="Number of samples to visualize.")

args = parser.parse_args()

if __name__ == "__main__":
    # Example usage:
    # smi_file = "path/to/file.smi"
    # smi_file = args.smi_file
    # n_samples = 100
    # n_samples = args.n_samples

    smi_file = args.smi_file
    # Check if the file exists
    try:
        with open(smi_file, 'r') as file:
            pass
    except FileNotFoundError:
        print(f"File {smi_file} not found.")
        exit(1)

    # load molecules from file
    mols = SmilesMolSupplier(smi_file, sanitize=True, nameColumn=-1)

    mols_list = [mol for mol in mols]
    n_samples = min(args.n_samples, len(mols_list))  # sample size should not exceed number of molecules
    mols_sampled = random.sample(mols_list, n_samples)  # sample 100 random molecules to visualize

    mols_per_row = int(math.sqrt(n_samples))            # make a square grid

    smi_file_names = smi_file.split(".")[:-1]  # remove file extension
    smi_file_name = ".".join(smi_file_names)  # name of file without extension
    png_filename = smi_file_name + ".png"  # name of PNG file to create
    labels=list(range(n_samples))       # label structures with a number

    # draw the molecules (creates a PIL image)
    img = MolsToGridImage(mols=mols_sampled,
                        molsPerRow=mols_per_row,
                        legends=[str(i) for i in labels])

    img.save(png_filename)
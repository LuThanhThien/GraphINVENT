from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_derivative(deriv_smiles: str, core_smiles: str) -> bool:
    # parse SMILES
    mol_deriv = Chem.MolFromSmiles(deriv_smiles)
    mol_core  = Chem.MolFromSmiles(core_smiles)
    if mol_deriv is None or mol_core is None:
        raise ValueError("Invalid SMILES")

    # simple substructure match test
    if mol_deriv.HasSubstructMatch(mol_core):
        return True

    # optionally: check scaffold equality instead of full-substructure
    # scaffold_deriv = MurckoScaffold.GetScaffoldForMol(mol_deriv)
    # scaffold_core  = MurckoScaffold.GetScaffoldForMol(mol_core)
    # return scaffold_deriv.ToSmiles() == scaffold_core.ToSmiles()

    return False

# example usage
# core   = "C1=CC=C(C=C1)C2=CC(=O)C3=C(O2)C=CC=C3O"       # e.g. a flavonoid core
core   = "c1ccccc1"
cand   = "C1=CC(=C(C=C1O)C2=CC(=O)C3=C(O2)C=CC=C3O)O"  # candidate derivative
print(is_derivative(cand, core))  # True

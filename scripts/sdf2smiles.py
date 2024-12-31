from rdkit import Chem

suppl = Chem.SDMolSupplier("lig.sdf")

for mol in suppl:
    if mol is not None:
        print(Chem.MolToSmiles(mol))

from rdkit import Chem

mol=Chem.MolFromSmiles('C')
w=Chem.SDWriter('C.sdf')
w.write(mol)
w.close()

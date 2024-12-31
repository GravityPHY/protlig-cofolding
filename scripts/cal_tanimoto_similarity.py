# This script will 
# 1. Process all the ligand in PDBBind to mol2
# 2. Calculate a list of Tanimoto similarities, 
#    each similarity stands for test molecule 
#    compare to each ligand appeared in PDBBind
# 3. Get the highest Tanimoto similarity in the list for each test molecule
# 4. Add the column with highest Tanimoto similarity to the initial dataframe
# 5. Save the updated dataframe in a .csv file and name as "with_tanimoto.csv"
import os
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

def get_pdbbind_ligands(pdbbind_data_path):
    ms=[]
    pdb_lists=os.listdir(pdbbind_data_path)
    for pdb_id in pdb_lists:
        pdb_path=os.path.join(pdbbind_data_path,pdb_id)
        mol2_file=os.path.join(pdb_path,f"{pdb_id}_ligand.mol2")
        try:
            mol = Chem.MolFromMol2File(mol2_file)
            if mol is not None:
                ms.append(mol)
        except Exception as e:
            pass
    return ms


def cal_similarity(ms,SMILES):
    fps = [FingerprintMols.FingerprintMol(x) for x in ms]
    mol = Chem.MolFromSmiles(SMILES)
    bins= [DataStructs.FingerprintSimilarity(FingerprintMols.FingerprintMol(mol),fp) for fp in fps]
    return bins

def get_max_similarity(ms,SMILES):
    fps = [FingerprintMols.FingerprintMol(x) for x in ms]
    mol = Chem.MolFromSmiles(SMILES)
    bins= [DataStructs.FingerprintSimilarity(FingerprintMols.FingerprintMol(mol),fp) for fp in fps]
    return max(bins)



tsv_path="~/PBI/single_chain/umol_benchmark.tsv"
df=pd.read_csv(tsv_path,sep="\t")

# find out sequence contain gap(s) '/' in the sequence
df_filter=df[~df['seqfromatoms'].str.contains('/')]
df_filter=df_filter.dropna(subset=['smiles'])

path="~/yuha1k/data/PDBBind_processed"
ms=get_pdbbind_ligands(path)

tanimoto_scores=[]
for row in df_filter.iterrows():
    tanimoto_scores.append(get_max_similarity(ms,row[1]['smiles']))

df_filter['Tanimoto']=tanimoto_scores
df_filter.to_csv("with_tanimoto.csv",index=False)
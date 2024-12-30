import os
import json
import argparse

parser = argparse.ArgumentParser(description = """Accept .json contains neccesary input information to make prediction with umol""")
parser.add_argument("--json_path",type=str)

#Parse args
args = parser.parse_args()
#Data
json_path = args.json_path

input_data=json.load(open(json_path))
base_path=os.path.dirname(json_path)

OUTDIR=input_data['OUTPUTDIR']
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

fasta_path=os.path.join(OUTDIR,f"{input_data['STRUCID']}.fasta")
with open(fasta_path,'w') as file:
    file.write(f">{input_data['STRUCID']}\n")
    file.write(input_data['EXPRSEQUENCE']) 

ID=input_data['STRUCID']
FASTA=fasta_path
POCKET_INDICES=f"{OUTDIR}{ID}_pocket_indices.npy" #Zero indexed numpy array of what residues are in the pocket (all CBs within 10Å from the ligand)
LIGAND_SMILES=input_data['SMILES'] #Make sure these are canonical as in RDKit. If you do not have SMILES - you can input a .sdf file to 'make_ligand_feats.py'
UNICLUST="/da/dmp/uniclust/uniclust30_2018_08/uniclust30_2018_08"


## Search Uniclust30 with HHblits to generate an MSA (a few minutes)
HHBLITS="/home/yuha1k/umol_conda/hh-suite/build/bin/hhblits"
os.system(f"{HHBLITS} -i {FASTA} -d {UNICLUST} -E 0.001 -all -oa3m {OUTDIR}{ID}.a3m")

os.system(f"python3 /home/yuha1k/umol_conda/src/make_msa_seq_feats.py \
            --input_fasta_path {FASTA} \
            --input_msas {OUTDIR}{ID}.a3m \
            --outdir {OUTDIR}")

os.system(f"python3 /home/yuha1k/umol_conda/src/make_ligand_feats.py \
            --input_smiles '{LIGAND_SMILES}' \
            --outdir {OUTDIR}")

os.system(f"python3 /home/yuha1k/umol_conda/src/make_targetpost_npy.py \
            --outfile '{POCKET_INDICES}' \
            --target_pos {input_data['POCKETPOS']}")

## Predict (a few minutes)
MSA_FEATS=os.path.join(OUTDIR,"msa_features.pkl")
LIGAND_FEATS=os.path.join(OUTDIR,"ligand_inp_features.pkl")
POCKET_PARAMS="/da/dmp/PBI/umol_params/params/params_pocket.npy" #Umol-pocket params in pixie
NO_POCKET_PARAMS="/da/dmp/PBI/umol_params/params/params_no_pocket.npy" #Umol no-pocket params in pixie
NUM_RECYCLES=3


#Change to no-pocket params if no pocket
#Then also leave out the target pos
os.system(f"python3 /home/yuha1k/umol_conda/src/predict.py  \
    --msa_features {MSA_FEATS} \
    --ligand_features {LIGAND_FEATS} \
    --id {ID} \
    --ckpt_params {POCKET_PARAMS} \
    --target_pos {POCKET_INDICES} \
    --num_recycles {NUM_RECYCLES} \
    --outdir {OUTDIR}")
print("prediction complete")
RAW_PDB=os.path.join(OUTDIR,f"{ID}_pred_raw.pdb")
os.system(f"python3 /home/yuha1k/umol_conda/src/relax/align_ligand_conformer.py \
            --pred_pdb {RAW_PDB} \
            --ligand_smiles '{LIGAND_SMILES}'\
            --outdir {OUTDIR}")
os.system(f"grep ATOM {OUTDIR}{ID}_pred_raw.pdb > {OUTDIR}{ID}_pred_protein.pdb")



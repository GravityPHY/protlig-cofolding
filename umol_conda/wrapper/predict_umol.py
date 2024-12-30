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

# make sure the saving directory exists
OUTDIR=input_data['OUTPUTDIR']
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

# make a fasta file in the saving directory
fasta_path=os.path.join(OUTDIR,f"{input_data['STRUCID']}.fasta")
with open(fasta_path,'w') as file:
    file.write(f">{input_data['STRUCID']}\n")
    file.write(input_data['EXPRSEQUENCE']) 

ID=input_data['STRUCID']
FASTA=fasta_path
POCKET_INDICES=f"{OUTDIR}{ID}_pocket_indices.npy" #Zero indexed numpy array of what residues are in the pocket (all CBs within 10Å from the ligand)
LIGAND_SMILES=input_data['SMILES'] #Make sure these are canonical as in RDKit. If you do not have SMILES - you can input a .sdf file to 'make_ligand_feats.py'
UNICLUST="/da/dmp/uniclust/uniclust30_2018_08/uniclust30_2018_08"
RANDOMSEED=input_data['random']
## Search Uniclust30 with HHblits to generate an MSA (a few minutes)
HHBLITS="/home/yuha1k/umol_conda/hh-suite/build/bin/hhblits"
os.system(f"{HHBLITS} -i {FASTA} -d {UNICLUST} -E 0.001 -all -oa3m {OUTDIR}{ID}.a3m -cpu 10")

# prepare MSA feature
os.system(f"python3 /home/yuha1k/umol_conda/src/make_msa_seq_feats.py \
            --input_fasta_path {FASTA} \
            --input_msas {OUTDIR}{ID}.a3m \
            --outdir {OUTDIR}")

# prepare ligand feature
os.system(f"python3 /home/yuha1k/umol_conda/src/make_ligand_feats.py \
            --input_smiles '{LIGAND_SMILES}' \
            --outdir {OUTDIR}")

#os.system(f"python3 /home/yuha1k/umol_conda/src/make_targetpost_npy.py \
#            --outfile '{POCKET_INDICES}' \
#            --target_pos {input_data['POCKETPOS']}")

## Predict (a few minutes)
MSA_FEATS=os.path.join(OUTDIR,"msa_features.pkl")
LIGAND_FEATS=os.path.join(OUTDIR,"ligand_inp_features.pkl")
POCKET_PARAMS="/da/dmp/PBI/umol_params/params/params_pocket.npy" #Umol-pocket params in pixie
NO_POCKET_PARAMS="/da/dmp/PBI/umol_params/params/params_no_pocket.npy" #Umol no-pocket params in pixie
NUM_RECYCLES=3


#Change to no-pocket params if no pocket
#Then also leave out the target pos
if not input_data["POCKETPOS"]:
    PARAMS=NO_POCKET_PARAMS
    POCKET_INDICES=None
else:
    PARAMS=POCKET_PARAMS
    os.system(f"python3 /home/yuha1k/umol_conda/src/make_targetpost_npy.py \
            --outfile '{POCKET_INDICES}' \
            --target_pos {input_data['POCKETPOS']}")

os.system(f"python3 /home/yuha1k/umol_conda/src/predict.py  \
    --msa_features {MSA_FEATS} \
    --ligand_features {LIGAND_FEATS} \
    --id {ID} \
    --ckpt_params {PARAMS} \
    --target_pos {POCKET_INDICES} \
    --num_recycles {NUM_RECYCLES} \
    --outdir {OUTDIR}\
    --randomseed {RANDOMSEED}")

RAW_PDB=os.path.join(OUTDIR,f"{ID}_pred_raw.pdb")
os.system(f"python3 /home/yuha1k/umol_conda/src/relax/align_ligand_conformer.py \
            --pred_pdb {RAW_PDB} \
            --ligand_smiles '{LIGAND_SMILES}'\
            --outdir {OUTDIR}")
os.system(f"grep ATOM {OUTDIR}{ID}_pred_raw.pdb > {OUTDIR}{ID}_pred_protein.pdb")


## Relax the protein (a few minutes)
#This fixes clashes mainly in the protein, but also in the protein-ligand interface.

PRED_PROTEIN=os.path.join(OUTDIR,f"{ID}_pred_protein.pdb")
PRED_LIGAND=os.path.join(OUTDIR,f"{ID}_pred_ligand.sdf")
RESTRAINTS="CA+ligand" # or "protein"
os.system(f"python3 /home/yuha1k/umol_conda/src/relax/openmm_relax.py --input_pdb {PRED_PROTEIN} \
                        --ligand_sdf {PRED_LIGAND} \
                        --file_name {ID} \
                        --restraint_type {RESTRAINTS} \
                        --outdir {OUTDIR}")

#Write plDDT to Bfac column
RAW_COMPLEX=os.path.join(OUTDIR,f"{ID}_pred_raw.pdb")
RELAXED_COMPLEX=os.path.join(OUTDIR,f"{ID}_relaxed_complex.pdb")
os.system(f"python3 /home/yuha1k/umol_conda/src/relax/add_plddt_to_relaxed.py  \
            --raw_complex {RAW_COMPLEX} \
            --relaxed_complex {RELAXED_COMPLEX}  \
            --outdir {OUTDIR}")
print(f"The final relaxed structure can be found at {OUTDIR}{ID}_relaxed_plddt.pdb")

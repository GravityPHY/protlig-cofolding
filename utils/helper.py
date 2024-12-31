import os
import re
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFMCS
from pymol import cmd
from Bio.PDB import PDBParser, PDBIO

def centroid(selection="all"):
    try:
        model = cmd.get_model(selection)
        x = sum([a.coord[0] for a in model.atom])/len(model.atom)
        y = sum([a.coord[1] for a in model.atom])/len(model.atom)
        z = sum([a.coord[2] for a in model.atom])/len(model.atom)
        return np.array([x,y,z])
    except Exception as e:
        return np.array([float('nan'),float('nan'),float('nan')])

def protein_rmsd(crystal_path,predict_path,CA=True):
    """
    Calculate the rmsd of of two provided models
    Parameters:
    crystal_path (str):
    predict_path (str):
    CA (bool): If True, return the CA rmsd; else, return the all-atom rmsd
    Return:
    rmsd (float)
    """
    cmd.reinitialize()
    cmd.load(crystal_path,"crystal")
    cmd.load(predict_path,"prediction")
    if CA:
        rmsd=cmd.align("prediction and name CA","crystal and name CA")[0]
    else:
        rmsd=cmd.align("prediction","crystal")[0]
    return rmsd

def get_pocket_resn(pdb_path,lig_chain,protein_chain,radius=10.0):
    cmd.load(pdb_path,"structure")
    cmd.select("ligand",f"chain {lig_chain} and resn LI1")
    if cmd.count_atoms("ligand")==0:
        cmd.select("ligand",f"chain {protein_chain} and resn LI1")
    if cmd.count_atoms("ligand")==0:
        cmd.select("ligand",f"byres (chain {protein_chain} around {str(radius)}) and resn LI1")
    if cmd.count_atoms("ligand")==0:
        cmd.select("ligand",f"chain {lig_chain} and resn INH")
    if cmd.count_atoms("ligand")==0:
        cmd.select("ligand",f"chain {lig_chain} and resn LG1")
    if cmd.count_atoms("ligand")==0:
        cmd.select("ligand","organic and not polymer")
    elif cmd.count_atoms("ligand")==0:
        raise Exception(f"Please check this file {pdb_path}")
    
    
    #cmd.select("ligand","organic and not polymer")
    #ligand_sele=f"chain {lig_chain}"
    protein_sele=f"chain {protein_chain}"
    fasta_sequence=cmd.get_fastastr(protein_sele)
    sequence=''.join(fasta_sequence.split('\n')[1:])
    seq_len=len(sequence)
    #res_sele=f"({protein_sele}) within 10.0 of ({ligand_sele})"
    res_sele=f"({protein_sele}) within {str(radius)} of ligand"
    cmd.select("res_sele",res_sele)
    #resn=cmd.get_model("res_sele")
    #res_nums=sorted(map(int, resn))
    residue_numbers=[]
    cmd.iterate("res_sele",
                "residue_numbers.append(resi)",
               space={'residue_numbers':residue_numbers})
    residue_numbers=list(set(residue_numbers))
    residue_numbers=[re.sub(r'\D', '', i) for i in residue_numbers]
    residue_numbers=sorted([int(i) for i in residue_numbers if int(i)<seq_len])
    residue_numbers=[str(i) for i in residue_numbers]
    cmd.delete("structure")
    cmd.delete("res_sele")
    residue_numbers_str=",".join(residue_numbers)
    return residue_numbers_str

def pocket_rmsd_delta(crystal_path,predict_path,pocketpos,CA=False,):
    cmd.reinitialize()
    cmd.load(crystal_path,"crystal")
    cmd.load(predict_path,"prediction")
    pocket_str="+".join(pocketpos.split(","))
    cmd.select("crystal_pocket",f"crystal and resi {pocket_str}")
    cmd.select("predict_pocket",f"prediction and resi {pocket_str}")
    if CA:
        pocket_rmsd=cmd.align("predict_pocket and name CA","crystal_pocket and name CA")[0]
    else:
        pocket_rmsd=cmd.align("predict_pocket","crystal_pocket")[0]
    crystal_centroid=centroid("crystal_pocket")
    predict_centroid=centroid("predict_pocket")
    dist = np.linalg.norm(predict_centroid-crystal_centroid)
    return pocket_rmsd, dist

def ligand_centroid_diff(crystal_path,predict_path,cry_lig_chain="A",pred_lig_chain="B",CA=False,):
    cmd.reinitialize()
    cmd.load(crystal_path,"crystal")
    cmd.load(predict_path,"prediction")
    #cmd.select("crystal_pocket",f"crystal and resi {pocket_str}")
    #cmd.select("predict_pocket",f"prediction and resi {pocket_str}")
    if CA:
        cmd.align("prediction and name CA","crystal and name CA")[0]
    else:
        cmd.align("prediction","crystal")[0]
    
    cmd.select("pred_lig",f"prediction and chain {pred_lig_chain}")
    
    cmd.select("crys_lig",f"crystal and chain {cry_lig_chain}")
        
    pred_centroid=centroid("pred_lig")
    crys_centroid=centroid("crys_lig")
    dist = np.linalg.norm(pred_centroid-crys_centroid)
    return dist

def get_coord(mol, indices=None):
    if indices is None:
        indices = tuple(range(mol.GetNumAtoms()))
    output = []
    for atom_id in indices:
        pos = mol.GetConformer().GetAtomPosition(atom_id)
        output.append((pos.x, pos.y, pos.z))
    return tuple(output)

def rmsd_calc(r_coord, m_coord):
    s = 0
    for r, m in zip(r_coord, m_coord):
        s += (r[0] - m[0]) ** 2 + (r[1] - m[1]) ** 2 + (r[2] - m[2]) ** 2
    s = (s / len(r_coord)) ** 0.5
    return s

def ligand_rmsd_0(crystal_path,predict_path,cry_lig_chain="A",pred_lig_chain="B",CA=False,save_dir='./',save=False):
    """
    Parameters:
    crystal_path (str):
    predict_path (str):
    CA (bool): If True, return the CA rmsd; else, return the all-atom rmsd
    save_dir (Path): a path that save ligand in .sdf format 
    save (bool): If True, keep the sdf under save_dir; else, remove the sdf files that were procued
    """
    cmd.reinitialize()
    cmd.load(crystal_path,"crystal")
    cmd.load(predict_path,"prediction")
    cmd.align("prediction","crystal")
    cmd.select("pred_lig",f"prediction and chain {pred_lig_chain}")
    cmd.select("crys_lig",f"crystal and chain {cry_lig_chain}")
    if cmd.count_atoms("crys_lig")==0:
        raise Exception(f"Please check this file {crystal_path}")

    if not os.path.exists(save_dir):
        print(f"{save_dir} does not exist... creating this directory now")
        os.makedirs(save_dir)
        
    pred_file=os.path.join(save_dir,"pred_lig.sdf")
    cmd.save(pred_file,"pred_lig")
    true_file=os.path.join(save_dir,"crys_lig.sdf")
    cmd.save(true_file,"crys_lig")

    pred = Chem.SDMolSupplier(pred_file,removeHs=False, sanitize=False)
    pred_mol = next(pred)

    true = Chem.SDMolSupplier(true_file,removeHs=False, sanitize=False)
    true_mol = next(true)
    
    match_indices = pred_mol.GetSubstructMatches(true_mol, uniquify=False, useChirality=True, maxMatches=10000)
    

    ref_coord = get_coord(true_mol)
    mol_coord = get_coord(pred_mol)
    min_rmsd = float('inf')
    #if not match_indices:
        # posebuster rmsd 
    #    try:
    #        rmsd=modules.rmsd.check_rmsd(pred_mol, true_mol)['results']['rmsd']
    #        return rmsd
    #    except Exception as e:
    #        return e
    for ids in match_indices:
        mol_coord = get_coord(pred_mol, ids)
        s = rmsd_calc(ref_coord, mol_coord)
        if s<min_rmsd:
            min_rmsd=s
    if min_rmsd==float('inf'):
        return None
    return min_rmsd

def CalcLigRMSD(lig1, lig2, rename_lig2=True, output_filename="tmp.pdb"):
    """
    Calculate the Root-mean-square deviation (RMSD) between two prealigned ligands, 
    even when atom names between the two ligands are not matching.
    The symmetry of the molecules is taken into consideration (e.g. tri-methyl groups). 
    Moreover, if one ligand structure has missing atoms (e.g. undefined electron density in the crystal structure), 
    the RMSD is calculated for the maximum common substructure (MCS).

    Parameters
    ----------
    lig1 : RDKit molecule
    lig2 : RDKit molecule
    rename_lig2 : bool, optional
        True to rename the atoms of lig2 according to the atom names of lig1
    output_filename : str, optional
        If rename_lig2 is set to True, a PDB file with the renamed lig2 atoms is written as output.
        This may be useful to check that the RMSD has been "properly" calculated, 
        i.e. that the atoms have been properly matched for the calculation of the RMSD.
    
    Returns
    -------
    rmsd : float
        Root-mean-square deviation between the two input molecules
    """
    # Exclude hydrogen atoms from the RMSD calculation
    #lig1 = Chem.RemoveHs(lig1)
    #lig2 = Chem.RemoveHs(lig2)
    # Extract coordinates
    coordinates_lig2 = lig2.GetConformer().GetPositions()
    coordinates_lig1 = lig1.GetConformer().GetPositions()
  # Calculate the RMSD between the MCS of lig1 and lig2 (useful if e.g. the crystal structures has missing atoms)
    res = rdFMCS.FindMCS([lig1, lig2])
    ref_mol = Chem.MolFromSmarts(res.smartsString)
  # Match the ligands to the MCS
  # For lig2, the molecular symmetry is considered:
  # If 2 atoms are symmetric (3 and 4), two indeces combinations are printed out
  # ((0,1,2,3,4), (0,1,2,4,3)) and stored in mas2_list
    mas1 = list(lig1.GetSubstructMatch(ref_mol))  # match lig1 to MCS
    mas2_list = lig2.GetSubstructMatches(ref_mol, uniquify=False)
  # Reorder the coordinates of the ligands and calculate the RMSD between all possible symmetrical atom matches
    coordinates_lig1 = coordinates_lig1[mas1]
    list_rmsd = []
    for match1 in mas2_list:
        coordinates_lig2_tmp = coordinates_lig2[list(match1)]
        diff = coordinates_lig2_tmp - coordinates_lig1
        list_rmsd.append(np.sqrt((diff * diff).sum() / len(coordinates_lig2_tmp)))  # rmsd
  # Return the minimum RMSD
    lig_rmsd = min(list_rmsd)
  # Write out a PDB file with matched atom names
    if rename_lig2:
        mas2 = mas2_list[np.argmin(list_rmsd)]
        correspondence_key2_item1 = dict(zip(mas2, mas1))
        atom_names_lig1 = [atom1.GetPDBResidueInfo().GetName() for atom1 in lig1.GetAtoms()]
        lig1_ResName = lig1.GetAtoms()[0].GetPDBResidueInfo().GetResidueName()
        for i, atom1 in enumerate(lig2.GetAtoms()):
            atom1.GetPDBResidueInfo().SetResidueName(lig1_ResName)
            if i in correspondence_key2_item1.keys():
                atom1.GetPDBResidueInfo().SetName(atom_names_lig1[correspondence_key2_item1[i]])
        Chem.MolToPDBFile(lig2, output_filename)
    return lig_rmsd

def ligand_rmsd(crystal_path,predict_path,cry_lig_chain="A",pred_lig_chain="B",CA=False,save_dir='./',save=False):
    """
    Parameters:
    crystal_path (str):
    predict_path (str):
    CA (bool): If True, return the CA rmsd; else, return the all-atom rmsd
    save_dir (Path): a path that save ligand in .sdf format 
    save (bool): If True, keep the sdf under save_dir; else, remove the sdf files that were procued
    """
    cmd.reinitialize()
    cmd.load(crystal_path,"crystal")
    cmd.load(predict_path,"prediction")
    cmd.align("prediction","crystal")
    cmd.select("pred_lig",f"prediction and chain {pred_lig_chain}")
    cmd.select("crys_lig",f"crystal and chain {cry_lig_chain}")
    if cmd.count_atoms("crys_lig")==0:
        raise Exception(f"Please check this file {crystal_path}")

    if not os.path.exists(save_dir):
        print(f"{save_dir} does not exist... creating this directory now")
        os.makedirs(save_dir)
        
    pred_file=os.path.join(save_dir,"pred_lig.sdf")
    cmd.save(pred_file,"pred_lig")
    true_file=os.path.join(save_dir,"crys_lig.sdf")
    cmd.save(true_file,"crys_lig")

    pred = Chem.SDMolSupplier(pred_file,removeHs=False, sanitize=False)
    pred_mol = next(pred)

    true = Chem.SDMolSupplier(true_file,removeHs=False, sanitize=False)
    true_mol = next(true)    
    return CalcLigRMSD(true_mol,pred_mol,rename_lig2=False)

def get_plddt(pdb_path,chain_name=None,calc_type="mean"):
    parser = PDBParser()
    structure = parser.get_structure("PDB", pdb_path)
    b_factors = []
    for chain in structure.get_chains():
        if chain.get_id() in chain_name:
            for residue in chain:
                for atom in residue:
                    b_factor = atom.get_bfactor()
                    if b_factor is not None:
                        b_factors.append(b_factor)
    mean_b_factor = np.mean(b_factors)
    return mean_b_factor

def get_pocket_plddt(predict_path,protein_chain="A", pred_lig_chain="B",radius=5,CA=False):
    cmd.reinitialize()
    cmd.load(predict_path,"prediction")
    cmd.select("pred_lig",f"prediction and chain {pred_lig_chain}")
    protein_sele=f"chain {protein_chain}"
    res_sele=f"({protein_sele}) within {str(radius)} of pred_lig"
    cmd.select("res_sele",res_sele)
    residue_numbers=[]
    cmd.iterate("res_sele",
                "residue_numbers.append((resi, resn, b))",
               space=locals())
    bs=[]
    for resi, resn, b in residue_numbers:
        bs.append(b)
    
    return np.mean(bs)
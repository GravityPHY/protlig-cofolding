import argparse
import os
import re
import sys

from rdkit import Chem
from rdkit.Chem import rdFMCS

from read_input import read_pdbqt, read_input
def get_coord(mol, indices=None):
    if indices is None:
        indices = tuple(range(mol.GetNumAtoms()))
    output = []
    for atom_id in indices:
        pos = mol.GetConformer().GetAtomPosition(atom_id)
        output.append((pos.x, pos.y, pos.z))
    return tuple(output)


def rmsd(mol, ref, chirality):

    def rmsd_calc(r_coord, m_coord):
        s = 0
        for r, m in zip(r_coord, m_coord):
            s += (r[0] - m[0]) ** 2 + (r[1] - m[1]) ** 2 + (r[2] - m[2]) ** 2
        s = (s / len(r_coord)) ** 0.5
        return s

    match_indices = mol.GetSubstructMatches(ref, uniquify=False, useChirality=chirality, maxMatches=10000)
    min_rmsd = float('inf')
    if not match_indices:
        mcs = rdFMCS.FindMCS([mol, ref], threshold=1.0,
                             ringMatchesRingOnly=False, completeRingsOnly=False,
                             matchChiralTag=chirality)
        if not mcs:
            return None
        patt = Chem.MolFromSmarts(mcs.smartsString)
        refMatch, molMatch = ref.GetSubstructMatches(patt, uniquify=False), \
                             mol.GetSubstructMatches(patt, uniquify=False)

        for ids_ref in refMatch:
            for ids_mol in molMatch:
                ref_coord = get_coord(ref, ids_ref)
                mol_coord = get_coord(mol, ids_mol)
                s = rmsd_calc(ref_coord, mol_coord)
                if s < min_rmsd:
                    min_rmsd = s
    else:
        ref_coord = get_coord(ref)
        for ids in match_indices:
            mol_coord = get_coord(mol, ids)
            s = rmsd_calc(ref_coord, mol_coord)
            if s < min_rmsd:
                min_rmsd = s

def rmsd_calc(r_coord, m_coord):
    s = 0
    for r, m in zip(r_coord, m_coord):
        s += (r[0] - m[0]) ** 2 + (r[1] - m[1]) ** 2 + (r[2] - m[2]) ** 2
    s = (s / len(r_coord)) ** 0.5
    return s

def get_coord(mol, indices=None):
    if indices is None:
        indices = tuple(range(mol.GetNumAtoms()))
    output = []
    for atom_id in indices:
        pos = mol.GetConformer().GetAtomPosition(atom_id)
        output.append((pos.x, pos.y, pos.z))
    return tuple(output)
chirality=True
#pred={m.GetProp('_Name'): m for m in Chem.SDMolSupplier("pose.sdf") if m}
mols=[mol for mol,mol_name in read_input("pose.sdf")]
mol=next(Chem.SDMolSupplier("pose.sdf"))
ref=next(Chem.SDMolSupplier("ligand.sdf"))#{m.GetProp('_Name'): m for m in Chem.SDMolSupplier("ligand.sdf") if m}
match_indices = mol.GetSubstructMatches(ref, uniquify=False, useChirality=chirality, maxMatches=10000)
print(match_indices)
ref_coord = get_coord(ref)
mol_coord = get_coord(mol)
min_rmsd = float('inf')
for ids in match_indices:
    mol_coord = get_coord(mol, ids)
    s = rmsd_calc(ref_coord, mol_coord)
    if s<min_rmsd:
        min_rmsd=s
print(min_rmsd)
for i,mol in enumerate(mols,1):
    mol_name = mol.GetProp('_Name')
    if isinstance(ref, dict):
        try:
            refmol = ref[mol.GetProp('_Name')]
        except KeyError:
            sys.stderr.write(f'Molecule with name {mol.GetProp("_Name")} is not available '
                                         f'in the reference SDF file\n')
            print(f'{mol_name}\t{i}\tNo matches')
            refmol = None
    else:
        refmol = ref
    rmsd(mol,refmol,chirality)

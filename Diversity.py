#!/usr/bin/env python
# coding: utf-8

from glob import iglob
from sys import argv

# Getting started with filtering rotatable bonds
# Import functions
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker


def main():
    master_folder = argv[1]
    for folder in iglob(f'{master_folder}/*'):
        ligand_preparation(folder)


def ligand_preparation(folder_path):
    # load database
    suppl = Chem.SmilesMolSupplier(f'{folder_path}/actives_final.smi')

    # remove all molecules that have more than 10 rot. bonds
    filtered_rotbonds = []
    for mol in suppl:
        if Descriptors.NumRotatableBonds(mol) <= 10:
            filtered_rotbonds.append(mol)

    # write filtered SMILES to new file
    with Chem.SmilesWriter(f'{folder_path}/RotBonds_filtered.smi') as w:
        for m in filtered_rotbonds:
            w.write(m)

    if len(filtered_rotbonds) == 0:
        print(f'{folder_path}: No molecules left after RotatableBonds filtering')
        return

    # take the prepared file and continue with the TPSA
    # grab Database, read and print TPSA values
    # filter all TPSA over 140

    filtered_TPSA = []
    for mol in filtered_rotbonds:
        if Descriptors.TPSA(mol) <= 140:
            filtered_TPSA.append(mol)

    # write filtered SMILES to new file
    with Chem.SmilesWriter(f'{folder_path}/RotTPSA.smi') as w:
        for m in filtered_TPSA:
            w.write(m)

    if len(filtered_TPSA) == 0:
        print(f'{folder_path}: No molecules left after TPSA filtering')
        return

    # remove everything inducing PAIN
    # import databases and PAINs catalog
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    catalog = FilterCatalog(params)

    # load data, compare with catalog and append non-PAIN compounds
    filtered_PAINs = []

    # separate into PAINful and 'good' compounds
    for mol in filtered_TPSA:
        if not catalog.HasMatch(mol):
            filtered_PAINs.append(mol)

    # write filtered SMILES to new file
    with Chem.SmilesWriter(f'{folder_path}/PAINs_filtered.smi') as w:
        for m in filtered_PAINs:
            w.write(m)

    if len(filtered_PAINs) == 0:
        print(f'{folder_path}: No molecules left after PAINs filtering')
        return

    # last but not least, look for lipinski violators
    # apply lipinski filter, append remaining compounds to mols-file
    filtered_Ro5 = []
    for mol in filtered_PAINs:
        molecular_weight = Descriptors.ExactMolWt(mol)
        logp = Descriptors.MolLogP(mol)
        h_bond_donor = Descriptors.NumHDonors(mol)
        h_bond_acceptors = Descriptors.NumHAcceptors(mol)

        if molecular_weight <= 500 and logp <= 5 and h_bond_donor <= 5 and h_bond_acceptors <= 10:
            filtered_Ro5.append(mol)

    # write filtered SMILES to new file
    with Chem.SmilesWriter(f'{folder_path}/actives_filtered.smi') as w:
        for m in filtered_Ro5:
            w.write(m)

    if len(filtered_Ro5) == 0:
        print(f'{folder_path}: No molecules left after Ro5 filtering')
        return

    # diversify by picking a randomized dataset

    # import data, start rdkit Chem and create Morgan Fingerprint
    with Chem.SmilesMolSupplier(f'{folder_path}/actives_filtered.smi') as suppl:
        ms = [a for a in suppl if a is not None]
    fps = [GetMorganFingerprint(b, 3) for b in ms]
    nfps = len(fps)

    def distij(i, j, fps=fps):
        return 1-DataStructs.DiceSimilarity(fps[i], fps[j])

    # picker for 30 molecules
    picker = MaxMinPicker()
    pickIndices = picker.LazyPick(distij, nfps, min(30, nfps), seed=23)

    # generate molecules instead of indices
    picks = [ms[x] for x in pickIndices]

    # write diversified SMILES to new file
    with Chem.SmilesWriter(f'{folder_path}/ligs_diversified.smi') as w:
        for m in picks:
            w.write(m)

    # Done


if __name__ == '__main__':
    main()

#!/usr/bin/env python

"""
Provides helper functions for converting between AMBER files and OpenEye OEMols.
"""

import numpy as np
from openeye.oechem import *
from openforcefield.typing.engines.smirnoff import *
from networkx.algorithms import isomorphism


def process_and_convert(smiles, reference_file, name):
    """
    Wrapper to map between a SMILES string and externally prepared `mol2` file.

    Parameters
    ----------
    smiles : str
        A SMILES string
    reference_file : str
        mol2 file containing the SMILES string
    name: str
        Title for the molecule

    """
    smiles_mol = process_smiles(smiles, name=name, add_hydrogens=True)
    reference_mol = load_mol2(reference_file, name=name, add_tripos=True)
    reference_to_target_mapping = atom_mapping(reference_mol, smiles_mol)
    smiles_mol = remap_charges(reference_to_target_mapping, reference_mol, smiles_mol)
    smiles_mol = remap_type(reference_to_target_mapping, reference_mol, smiles_mol)
    smiles_mol = remap_coordinates(reference_to_target_mapping, reference_mol, smiles_mol)
    print('Either use the returned `OEMol()` or save as `mol2` to continue...')
    return smiles_mol


def load_mol2(filename, name=None, add_tripos=True):
    ifs = oemolistream()
    if not ifs.open(filename):
        print(f'Unable to open {filename} for reading...')
    for mol in ifs.GetOEGraphMols():
        if add_tripos:
            OETriposAtomNames(mol)
        if name:
            mol.SetTitle(name)
        return mol


def save_mol2(filename, mol, standard=True):
    ofs = oemolostream()
    if not ofs.open(filename):
        print(f'Unable to open {filename} for writing...')
    if standard:
        OEWriteMolecule(ofs, mol)
    ofs.close()


def add_waters_and_ions():
    smiles = ['[Na+]', '[Cl-]', 'O']
    molecules = []
    for molecule in smiles:
        molecules.append(process_smiles(molecule))
    return molecules


def process_smiles(string, name=None, add_hydrogens=True, add_tripos=True):
    mol = OEMol()
    OESmilesToMol(mol, string)
    if add_hydrogens:
        OEAddExplicitHydrogens(mol)
    if add_tripos:
        OETriposAtomNames(mol)
    if name:
        mol.SetTitle(name)
    return mol


def atom_mapping(reference, target):
    reference_topology = create_topology(reference)
    target_topology = create_topology(target)

    reference_graph = create_graph(reference_topology)
    target_graph = create_graph(target_topology)

    reference_to_target_mapping = dict()
    graph_matcher = isomorphism.GraphMatcher(reference_graph, target_graph)
    if graph_matcher.is_isomorphic():
        print('Reference → Target')
        for (reference_atom, target_atom) in graph_matcher.mapping.items():
            reference_to_target_mapping[reference_atom] = target_atom
            reference_name = reference.GetAtom(OEHasAtomIdx(reference_atom)).GetName()
            reference_type = reference.GetAtom(OEHasAtomIdx(reference_atom)).GetType()
            target_name = target.GetAtom(OEHasAtomIdx(target_atom)).GetName()
            target_type = target.GetAtom(OEHasAtomIdx(target_atom)).GetType()

            print(f'({reference_name:4} [{reference_type:2}] ) {reference_atom:3d} → '
                  f'{target_atom:3d} ({target_name:4})')
    else:
        print('Graph is not isomorphic.')
    return reference_to_target_mapping


def create_topology(mol):
    return generateTopologyFromOEMol(mol)


def create_graph(topology):
    return generateGraphFromTopology(topology)


def remap_charges(reference_to_target_mapping, reference_mol, target_mol):
    print('Remap charges...')
    print('Existing → New')
    assert reference_mol.GetMaxAtomIdx() == target_mol.GetMaxAtomIdx()
    for (reference_atom, target_atom) in reference_to_target_mapping.items():
        reference = reference_mol.GetAtom(OEHasAtomIdx(reference_atom))
        reference_chg = reference.GetPartialCharge()

        target = target_mol.GetAtom(OEHasAtomIdx(target_atom))
        # This will be zero, if the target molecule was built from SMILES
        target_chg = target.GetPartialCharge()
        target_name = target.GetName()

        target.SetPartialCharge(reference_chg)
        print(f'({target_name:4}) {target_chg:+04f} → {target.GetPartialCharge():+04f}')
    return target_mol


def remap_type(reference_to_target_mapping, reference_mol, target_mol):
    print('Remap atom types...')
    print('Existing → New')
    assert reference_mol.GetMaxAtomIdx() == target_mol.GetMaxAtomIdx()
    for (reference_atom, target_atom) in reference_to_target_mapping.items():
        reference = reference_mol.GetAtom(OEHasAtomIdx(reference_atom))
        reference_type = reference.GetType()

        target = target_mol.GetAtom(OEHasAtomIdx(target_atom))
        # This will be zero, if the target molecule was built from SMILES
        target_type = target.GetType()
        target_name = target.GetName()

        target.SetType(reference_type)
        print(f'({target_name:4}) {target_type:4} → {target.GetType():4}')
    return target_mol


def remap_coordinates(reference_to_target_mapping, reference_mol, target_mol):
    print('Remap coordinates...')
    print('Existing → New')
    assert reference_mol.GetMaxAtomIdx() == target_mol.GetMaxAtomIdx()
    mapped_coordinates = np.zeros((reference_mol.GetMaxAtomIdx(), 3))
    reference_coordinates = reference_mol.GetCoords()
    for (reference_atom, target_atom) in reference_to_target_mapping.items():
        # print(reference_coordinates[reference_atom])
        # print(mapped_coordinates[target_atom])

        mapped_coordinates[target_atom] = reference_coordinates[reference_atom]
        current_coordinates = target_mol.GetCoords()[target_atom]
        print(f'{current_coordinates} → {mapped_coordinates[target_atom]}')
    target_mol.SetCoords(mapped_coordinates.flatten())
    return target_mol

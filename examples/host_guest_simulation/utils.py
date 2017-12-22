#!/usr/bin/env python

"""
Provides helper functions for converting between AMBER files and OpenEye OEMols.
"""

import numpy as np
import subprocess as sp
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
    smiles_mol = remap_names(reference_to_target_mapping, reference_mol, smiles_mol)
    smiles_mol = remap_charges(reference_to_target_mapping, reference_mol, smiles_mol)
    resname = parse_residue_name(reference_file)
    smiles_mol = remap_residues(reference_to_target_mapping, reference_mol, smiles_mol,
                                resname=resname)
    # smiles_mol = remap_type(reference_to_target_mapping, reference_mol, smiles_mol)
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


def save_mol2(filename, mol, standard=True, add_tripos=True):
    # I am not sure why, but I cannot get this to write the residue name.
    ofs = oemolostream()
    if not ofs.open(filename):
        print(f'Unable to open {filename} for writing...')
    if add_tripos:
        OETriposAtomNames(mol)
    if standard:
        OEWriteMolecule(ofs, mol)
    else:
        # ofs.SetFormat(OEFormat_MOL2H)
        flavor = OEOFlavor_Generic_Default | OEOFlavor_MOL2H_Default | OEOFlavor_MOL2H_AtomTypeNames
        ofs.SetFlavor(OEFormat_MOL2H, flavor)
        OEWriteMol2File(ofs, mol)
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
        print('Determining mapping...')
        print('Reference → Target')
        for (reference_atom, target_atom) in graph_matcher.mapping.items():
            reference_to_target_mapping[reference_atom] = target_atom
            reference_name = reference.GetAtom(OEHasAtomIdx(reference_atom)).GetName()
            reference_type = reference.GetAtom(OEHasAtomIdx(reference_atom)).GetType()
            target_name = target.GetAtom(OEHasAtomIdx(target_atom)).GetName()

            print(f'({reference_name:4} {reference_type:5}) {reference_atom:3d} → '
                  f'{target_atom:3d} ({target_name:4})')
    else:
        print('Graph is not isomorphic.')
    return reference_to_target_mapping


def create_topology(mol):
    return generateTopologyFromOEMol(mol)


def create_graph(topology):
    return generateGraphFromTopology(topology)


def remap_charges(reference_to_target_mapping, reference_mol, target_mol):
    print('Remapping charges...')
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


def remap_names(reference_to_target_mapping, reference_mol, target_mol):
    print('Remapping atom names...')
    print('Reference → Target')
    assert reference_mol.GetMaxAtomIdx() == target_mol.GetMaxAtomIdx()
    for (reference_atom, target_atom) in reference_to_target_mapping.items():
        reference = reference_mol.GetAtom(OEHasAtomIdx(reference_atom))
        reference_name = reference.GetName()

        target = target_mol.GetAtom(OEHasAtomIdx(target_atom))
        target_name = target.GetName()
        target.SetName(reference_name)
        print(f'{reference_name:4} → {target_name:4}')
    return target_mol


def parse_residue_name(input_mol2, path='./'):
    p = sp.Popen(['awk', '{print $8}', input_mol2], cwd=path, stdout=sp.PIPE)
    for line in p.stdout:
        if line.decode("utf-8").split() != []:
            name = line.decode("utf-8").split()[0]
            print(f'Found residue name = {name}')
            return name
        else:
            pass
    return None


def remap_residues(reference_to_target_mapping, reference_mol, target_mol, resname=None):
    print('Remapping residue names and numbers...')
    print('Existing → New')
    assert reference_mol.GetMaxAtomIdx() == target_mol.GetMaxAtomIdx()
    # It's not clear to me that we have to loop over all the atoms in the molecule,
    # if OpenEye knows they are connected properly, then setting the residue name
    # and residue number for one atom should be enough, but keeping this loop
    # seems safe and won't be a bottleneck.
    for (reference_atom, target_atom) in reference_to_target_mapping.items():
        reference = reference_mol.GetAtom(OEHasAtomIdx(reference_atom))
        reference_residue = OEAtomGetResidue(reference)
        reference_resname = reference_residue.GetName()
        reference_resnum = reference_residue.GetResidueNumber()
        # I believe this gets set to 'UNL' if OpenEye can't recognize the residue name.
        # Thus, I'm adding an override to manually set the residue name, and simply.

        target = target_mol.GetAtom(OEHasAtomIdx(target_atom))
        target_name = target.GetName()
        target_residue = OEAtomGetResidue(target)
        target_resname = target_residue.GetName()
        target_resnum = target_residue.GetResidueNumber()

        if resname is not None:
            target_residue.SetName(resname)
        else:
            target_residue.SetName(reference_resname)
        target_residue.SetResidueNumber(reference_resnum)

        print(f'({target_name:4}) {target_resname:4} {target_resnum:4} → '
              f'{target_residue.GetName():4} {target_residue.GetResidueNumber():4}')
    return target_mol


def remap_type(reference_to_target_mapping, reference_mol, target_mol):
    print('Remapping atom types...')
    print('Existing → New')
    assert reference_mol.GetMaxAtomIdx() == target_mol.GetMaxAtomIdx()
    for (reference_atom, target_atom) in reference_to_target_mapping.items():
        reference = reference_mol.GetAtom(OEHasAtomIdx(reference_atom))
        reference_type = reference.GetType()

        target = target_mol.GetAtom(OEHasAtomIdx(target_atom))
        # This will be None, if the target molecule was built from SMILES
        target_type = target.GetType()
        target_name = target.GetName()

        target.SetType(reference_type)
        print(f'({target_name:4}) {target_type:4} → {target.GetType():4}')
    return target_mol


def remap_coordinates(reference_to_target_mapping, reference_mol, target_mol):
    print('Remapping coordinates...')
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

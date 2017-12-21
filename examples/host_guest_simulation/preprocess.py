#!/usr/bin/env python

"""
Create an OpenMM compatible PDB file using AMBER inputs.
"""

import subprocess as sp


def create_pdb_with_conect(solvated_pdb, amber_prmtop, output_pdb, path='./'):
    """
    Create a PDB file containing CONECT records.
    This is not very robust, please manually check the `cpptraj` output.
    `cpptraj` must be in your PATH.

    Parameters
    ----------
    solvated_pdb : str
        Existing solvated structure from e.g., Mobley's Benchmark Sets repository
    amber_prmtop : str
        AMBER (or other) parameters for the residues in the solvated PDB file
    output_pdb : str
        Output PDB file name
    path : str
        Directory for input and output files

    """
    cpptraj = \
    f'''
    parm {amber_prmtop}
    trajin {solvated_pdb}
    trajout {output_pdb} conect
    '''

    cpptraj_input = output_pdb + '.in'
    cpptraj_output = output_pdb + '.out'

    with open(path + cpptraj_input, 'w') as file:
        file.write(cpptraj)
    with open(path + cpptraj_output, 'w') as file:
        p = sp.Popen(['cpptraj', '-i', cpptraj_input], cwd=path,
                     stdout=file, stderr=file)
        output, error = p.communicate()
    if p.returncode == 0:
        print('PDB file written by cpptraj.')
    elif p.returncode == 1:
        print('Error returned by cpptraj.')
        print(f'Output: {output}')
        print(f'Error: {error}')
    else:
        print(f'Output: {output}')
        print(f'Error: {error}')


def prune_conect(input_pdb, output_pdb, path='./'):
    """
    Delete CONECT records that correspond only to water molecules.
    This is necessary to be standards-compliant.
    This is not very robust.

    Parameters
    ----------
    input_pdb : str
        Input PDB file name
    output_pdb : str
        Output PDB file name
    path : str
        Directory for input and output files

    """
    p = sp.Popen(['grep', '-m 1', 'WAT', input_pdb], cwd=path, stdout=sp.PIPE)
    for line in p.stdout:
        first_water_residue = int(float(line.decode("utf-8").split()[1]))
        print(f'First water residue = {first_water_residue}')

    p = sp.Popen(['egrep', '-n', f'CONECT [ ]* {first_water_residue}', input_pdb],
                 cwd=path, stdout=sp.PIPE)
    for line in p.stdout:
        line_to_delete_from = int(float(line.decode("utf-8").split(':')[0]))
        print(f'Found first water CONECT entry at line = {line_to_delete_from}')

    with open(path + output_pdb, 'w') as file:
        sp.Popen(
         ['awk', f'NR < {line_to_delete_from}', input_pdb], cwd=path, stdout=file)

        sp.Popen(['echo', 'END'], cwd=path, stdout=file)
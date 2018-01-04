#!/usr/bin/env python

"""
Postprocess the `mol2` file to correct residue name and numbering.
"""

import subprocess as sp


def extract_coordinates(filename, path='./'):

    p = sp.Popen(['awk', '/@<TRIPOS>ATOM/{flag=1;next}/@<TRIPOS>BOND/{flag=0}flag {print $3,$4,$5}',
                  filename], cwd=path, stdout=sp.PIPE)
    coordinates = []
    for line in p.stdout:
        coordinate = [float(i) for i in line.decode("utf-8").split()]
        coordinates.append(coordinate)

    return coordinates


def mutate_coordinates(filename, coordinates, path='./'):
    p = sp.call(['cp', filename, 'tmp.mol2'], cwd=path, stdout=sp.PIPE)

    # Get index and atom, then write coordinates?
    # Then worry about extra columns? 

    return
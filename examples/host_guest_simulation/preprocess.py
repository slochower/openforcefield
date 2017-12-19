def create_pdb_with_conect(solvated_pdb, amber_prmtop, pdb_output):
    cpptraj = \
    '''
    parm {}
    trajin {}
    trajout {} conect
    '''.format(amber_prmtop, solvated_pdb, pdb_with_conect)

    with open(path + cpptraj_input, 'w') as file:
        file.write(cpptraj)
    cpptraj_output = sp.check_output(['cpptraj', '-i', cpptraj_input], cwd=path)

def prune_conect():
    pass

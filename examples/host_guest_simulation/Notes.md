# Some things that work...

# Some things that don't...
Many of these things might be able to work, given a little bit more ingenuity and tenacity. I don't want to claim they are broken (though some of them certainly seem so), only that I haven't found these paths to give a viable outcome.

- Write `mol2` via ParmEd (e.g., `structure.save(file.mol2)`) after using `oemmutils` to convert an `OEMol` to an OpenMM topology. The `mol2` file won't have any atom types. Didn't investigate why.
- Start from `mol2` files with GAFF atom types.
- Start from `mol2` files with water or ions. These are more easily generated directly as `OEMol`s.
- Build cyclodextrin (or other host macrocycle) from SMILES and expect reasonable coordinates.
- Use OpenEye to generate charges for an entire host molecule. It will take too long and the conformations generated are nonphysical and may impact the quality of the charges.
- Read in a PDB file without `CONECT` records.
- Read in a PDB file with `CONECT` records between water atoms, as can happen with TIP3P.
- Read in a PDB file with water residue named `WAT` instead of `HOH` (**check**)
- Save a cyclic molecule `prmtop`


# To do...

- Create a `frcmod` file via `parmchk`...
- Create aligned PDBs
- See if current `mol2` and `prmtop` files can actually run a simulation
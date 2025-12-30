## Example 2: Converting HSX into Tight-binding Parameters

### Copy-paste

```bash
./hsx2tb
```

### Explanation

To run this example, please follow the instructions found in the `example/README.md` and then in `example/ex_1_siesta/README.md`.  Here we just need to make sure we have the following files
- HSX_out: human-readable version of `*.HSX` file for the Hamiltonian (H), Overlap (S), and position (X) data.
- ORB_INDX: indices labeling the atomic orbitals we consider.
- input.fdf: input file for `SIESTA` calculation, we read the lattice vectors from here.
- kpoints.dat: hand-made file for defining our k-point path.  Takes the same form as the `SIESTA` k-point path block.

When we have these, type
```bash
./hsx2tb
```

This will output the following files 
- `TB_PARAMS`: contains all hopping/overlap parameters and the Hamiltonian/overlap onsite matrices.
- `TB_EIGS`: the resultant bands from the resultant tight-binding parameters.

To make sure we have the proper output files here, compare with the ones found in `/expected_output`.
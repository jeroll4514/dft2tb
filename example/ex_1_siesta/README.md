## Example 1: Working with SIESTA

### Copy-paste

```bash
siesta < input.fdf > output
mv *HSX HSX
mv *ORB_INDX ORB_INDX
./make_hsx_readable
cp {HSX_out,ORB_INDX} ../ex_2_hsx2tb
```

### Explanation

To run this example, please follow the instructions found in the `example/README.md`.  We now need to ensure we are within an environment to run a `SIESTA` DFT calculation.  We will do so here with the provided `input.fdf` file and `*.psf` pseudopotentials (obtained [here](https://nninc.cnf.cornell.edu/) for `SIESTA`).  Note how it includes the `SaveHS = T` flag, which will be an input file later

```bash
siesta < input.fdf > output
```

There will be many output files, but we only care about the following:
- SnS.HSX
- SnS.ORB_IDX

If we look within `dft2tb/src/hsx_m.f90`, we actually need to rename `SnS.HSX` to `HSX`.  Then we can run our executable to turn it into human-readable code.

```bash
mv *HSX HSX
./make_hsx_readable
```

Finally, if we look within `dft2tb/src/hsx2tb.f90`, we actually need to rename `SnS.ORB_INDX` to `ORB_INDX` to properly input the data.

```bash
mv *ORB_INDX ORB_INDX
```

We are now done!  Within this directory should be our `HSX_out` file, which is the human-readable version and also is the input for our next example.  Then we need to just move these files into the next directory:

```bash
cp {HSX_out,ORB_INDX} ../ex_2_hsx2tb
```

To make sure we have the proper output files here, compare with the ones found in `/expected_output`.
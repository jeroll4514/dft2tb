## Example

This directory walks through how one would go about compiling and using this code to determine the tight-binding parameters for monolayer SnS.  First we need to compile the project, first make sure that you are within an environment that has a `gfortran` compiler for FORTRAN and LAPACK.  From the current `example` directory, type

```bash
cd ..
make
```

This will generate two folders:
1. build: this contains all .o and .mod files during compilation process
2. executables: this contains all executables we will need to call.

We need to move the executables into our running directory to ensure the proper files are inputted.

```bash
cp executables/make_hsx_readable example/ex_1_siesta
cp executables/hsx2tb example/ex_2_hsx2tb
cp executables/tb2bands example/ex_3_tb2bands
cd example
```

To go through the whole process, we must do things in the following order
1. Run a siesta calculation to get our Hamiltonian (H) and Overlap (S) matrices.  Then turn then from machine code into human-readable code.
2. Turn these H/S matrices into their associated tight-binding parameters.
3. Use these tight-binding parameters to generate a bandstructure.

This corresponds with three subdirectories here.
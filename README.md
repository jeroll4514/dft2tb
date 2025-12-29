# dft2tb

This repository contains the FORTRAN code developed for my honors undergraduate senior thesis in physics (Spring 2023).  The code inputs the standard .HSX output from a [SIESTA](https://siesta-project.org/siesta/index.html) DFT calculation and computes tight-binding parameters.  These parameters are determined using the [LCAO method](https://doi.org/10.1103/PhysRev.94.1498) for our two-center integrals: given the Hamiltonian/Overlap matrix elements, we solve for the parameters that yield the closest results.

## How to compile
Make sure you are within an environment (e.g. conda) that includes the gfortran compiler for FORTRAN and includes the LAPACK library.  When in the main directory, type
- make $\leftarrow$ compiles the project
- make clean $\leftarrow$ removes all executables and object files created within the compilation process.

All executables are stored within the `executables/` directory upon compilation.

## Directories
- `build`: upon compilation will store all .o and .mod files
- `executables`: upon compilation will store all executables
- `src`: source codes for all processes

## `SIESTA` version notice
In any version of Siesta newer than 5.0 the HSX file structure has changed.

It no longer natively supports changing the file from .HSX machine code into human code.  To make this change, I took the 'human' `write_hsx` code from a prior version (hsx.f90 instead of hsx_m.f90) and copied it into the current `hsx_m.f90` code.  I then had to include `hsx%has_xij = .True.` after reading it from the .HSX file or else the code will not complete.

### Notes
This code reflects my early training in scientific computing and predates my current C++-based research framework.  I did not use Git at the time of its development, so I cleaned up the code for this current repository.

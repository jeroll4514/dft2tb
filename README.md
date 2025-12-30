# dft2tb

This repository contains the FORTRAN code developed for my honors undergraduate senior thesis in physics (Graduated Spring 2023).  The code inputs the standard .HSX output from a [SIESTA](https://siesta-project.org/siesta/index.html) DFT calculation and computes tight-binding parameters.  These parameters are determined using the [LCAO method](https://doi.org/10.1103/PhysRev.94.1498) for our two-center integrals: given the Hamiltonian/Overlap matrix elements, we solve for the parameters that yield the closest results.  If resultant bands are not sufficiently accurate, this allows for a great starting point for optimization processes.

## Quickstart
```bash
git clone https://github.com/jeroll4514/dft2tb.git
cd dft2tb
make
```

## How to compile
Ensure you are within an environment (e.g., conda) that includes the gfortran compiler for FORTRAN and the LAPACK library.  When in the main directory, type
- make $\leftarrow$ compiles the project
- make clean $\leftarrow$ removes all executables and object files created within the compilation process.

All executables are stored within the `executables` directory upon compilation.

## Directories
- `build`: upon compilation, will store all .o and .mod files
- `executables`: upon compilation, will store all executables
- `src`: source codes for all processes

## Example workflow (end-to-end)

A complete worked example is provided in the `example/` directory, demonstrating
the full workflow:

1. Run a SIESTA DFT calculation and extract Hamiltonian/overlap matrices.
2. Convert binary HSX output into a human-readable format.
3. Generate tight-binding parameters and compute the resulting band structure.

Each step is documented with copy-paste commands and explanations in the
corresponding subdirectories:

- `example/ex_1_siesta/`
- `example/ex_2_hsx2tb/`

This example uses monolayer SnS as a reference system.

## Citation

If you use this code in academic work, please cite:

> J. Roll, *Automating the determination of tight-binding representations for two-dimensional ferroelectrics*,
> Physics Undergraduate Honors Thesis, University of Arkansas (2023).  
> https://scholarworks.uark.edu/physuht/22

## `SIESTA` version notice
In any version of Siesta newer than 5.0 the HSX file structure has changed.

It no longer natively supports changing the file from .HSX machine code into human code.  To make this change, I took the 'human' `write_hsx` code from a prior version (hsx.f90 instead of hsx_m.f90) and copied it into the current `hsx_m.f90` code.  I then had to include `hsx%has_xij = .True.` after reading it from the .HSX file or else the code will not complete.

### Notes
This code reflects my early training in scientific computing and predates my current C++-based research framework.  Originally developed during my 2022–2023 honors thesis, I later added a Makefile and reorganized the repository for reproducibility and archival.  This work is preserved here as a record of my development as a computational scientist.

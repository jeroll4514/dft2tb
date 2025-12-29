# DFT2TB

This repository contains the FORTRAN code developed for my honors undergraduate senior thesis in physics (Spring 2023).  The code inputs the standard .HSX output from a [SIESTA](https://siesta-project.org/siesta/index.html) DFT calculation and computes tight-binding parameters.  These parameters are determined using the [LCAO method](https://doi.org/10.1103/PhysRev.94.1498) for our two-center integrals: given the Hamiltonian/Overlap matrix elements, we solve for the parameters that yield the closest results.

### Notes
This code reflects my early training in scientific computing and predates my current C++-based research framework.  I did not use Git at the time of its development, so I cleaned up the code for this current repository.

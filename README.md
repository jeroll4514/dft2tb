# DFT2TB

This repository contains the FORTRAN code developed for my honors undergraduate senior thesis in physics (Spring 2023).  The code inputs the standard .HSX output from a `SIESTA` DFT calculation and computes tight-binding parameters.  These parameters are determined from the LCAO method for our two-center integrals: given the Hamiltonian/Overlap matrix element, we solve for what parameters yield the closest results.

### Notes
This code reflects my early training in scientific computing and predates my current C++-based research framework.  I did not use Git at the time of its development, so I cleaned up the code for this current repository.

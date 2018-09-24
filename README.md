# AMR-LBM-OpenMP-2D
Adaptive Mesh Refinement (AMR) + Lattice-Boltzmann Method (LBM) + OpenMP

Based on the following papers:

A. Fakhari, M. Geier, T. Lee, “A mass-conserving lattice Boltzmann method with dynamic grid refinement for immiscible two-phase flows”, Journal of Computational Physics 315:434–457 (2016) [https://doi.org/10.1016/j.jcp.2016.03.058].

A. Fakhari, D. Bolster, L.-S. Luo, “A weighted multiple-relaxation-time lattice Boltzmann method for multiphase flows and its application in partial coalescence cascades”, Journal of Computational Physics 341:22–43 (2017) [https://doi.org/10.1016/j.jcp.2017.03.062].

This package is for numerical simulation of multiphase flows (rising bubble).
It employs OpenMP parallelization for share-memory machines.

The computational setup is a 2D bubble rising due to buoyancy force.
The computational domain is periodic in the x-direction and no-slip (bounce-back) at the bottom and top boundaries.

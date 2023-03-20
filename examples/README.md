# Workflow of Normal Dynamics simulations with PINDOL

The main steps to run and analyse a Normal Dynamics (ND) simulation are the following:

1. [Dynamical matrix diagonalization](#dynamical-matrix-diagonalization): generate *qmatrix.nd* and *freq.nd* with [QPOINTS](https://github.com/acammarat/phtools/tree/main/qpoints)
2. [Fourier transform of force constants](#fourier-transform-of-force-constants): generate *phi.nd* with [**phind**](https://github.com/acammarat/pindol/tree/main/phind)
3. Run the ND simulation with [**pindol**](https://github.com/acammarat/pindol/tree/main/pindol)
4. Analyse the results
   1. Convert the trajectory, velocities and accelerations in XYZ format
   2. Obtain the effective force constants which include the effect of the temperature

## Dynamical matrix diagonalization

## Fourier transform of force constants


...

## Contributions, bug reports and feature requests

We are happy to accept contributions. To report bugs or request new features, please use the [Issue Tracker](https://github.com/acammarat/pindol/issues). If you use the programs in this repository in your work, please send an email to cammaant [at] fel.cvut.cz - we will collect them and put up a list of outputs.

# Workflow of Normal Dynamics simulations with PINDOL

The main steps to run and analyse a Normal Dynamics (ND) simulation are the following:

1. [Dynamical matrix diagonalization](#dynamical-matrix-diagonalization): generate *qmatrix.nd* and *freq.nd* with [QPOINTS](https://github.com/acammarat/phtools/tree/main/qpoints)
2. [Fourier transform of anharmonic force constants](#fourier-transform-of-anharmonic-force-constants): generate *phi.nd* with [**phind**](https://github.com/acammarat/pindol/tree/main/phind)
3. [Normal Dynamics simulation](#normal-dynamics-simulation): generate the ND trajectory with [**pindol**](https://github.com/acammarat/pindol/tree/main/pindol)
4. [Analyse the results](#analyse-the-results)
   1. Convert the trajectory, velocities and accelerations in XYZ format
   2. Obtain the effective force constants which include the effect of the temperature

At the moment, the geometry reference file must be created as [VASP](https://www.vasp.at) POSCAR format for all the pre- and postprocessing steps; for the ND run, also the [LAMMPS](https://www.lammps.org) format is accepted. In the following example, we will assume that the geometry is in POSCAR format.
   
## Dynamical matrix diagonalization

The [qpoints](https://github.com/acammarat/phtools/tree/main/qpoints) code is used in this step. The harmonic force constants, eigenvectors and frequencies can be obtained, in principle, with any calculator. At the moment, only [PHONOPY]( https://phonopy.github.io/phonopy) is supported; for this reason, the following example will refer to its usage.

First, we need to choose a list of q-points which constitutes the q-set; for example

```
   0.000000000000000    0.000000000000000    0.000000000000000
   0.500000000000000    0.000000000000000    0.000000000000000
   0.500000000000000    0.166666666666667    0.000000000000000
   0.500000000000000    0.333333333333333    0.000000000000000
   0.500000000000000    0.500000000000000    0.000000000000000
   0.500000000000000    0.333333333333333    0.166666666666667
  -0.333333333333333    0.500000000000000    0.166666666666667
   0.166666666666667    0.000000000000000    0.000000000000000
   0.333333333333333    0.000000000000000    0.000000000000000
   0.166666666666667    0.166666666666667    0.000000000000000
   0.333333333333333    0.166666666666667    0.000000000000000
  -0.333333333333333    0.166666666666667    0.000000000000000
  -0.166666666666667    0.166666666666667    0.000000000000000
   0.333333333333333    0.333333333333333    0.000000000000000
  -0.333333333333333    0.333333333333333    0.000000000000000
  -0.333333333333333    0.333333333333333    0.166666666666667
  
  ```
  which is the list of unique points contained in a 5x5x5 Mokhorst-Pack mesh obtained with the aid of the [SPGLIB](https://spglib.readthedocs.io/en/latest) library.



## Fourier transform of anharmonic force constants

## Normal Dynamics simulation

## Analyse the results


...

## Contributions, bug reports and feature requests

We are happy to accept contributions. To report bugs or request new features, please use the [Issue Tracker](https://github.com/acammarat/pindol/issues). If you use the programs in this repository in your work, please send an email to cammaant [at] fel.cvut.cz - we will collect them and put up a list of outputs.

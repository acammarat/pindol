# Workflow of Normal Dynamics simulations with PINDOL

The main steps to run and analyse a Normal Dynamics (ND) simulation are the following:

1. [Dynamical matrix diagonalization](#dynamical-matrix-diagonalization): generate *qmatrix.nd* and *freq.nd* with [QPOINTS](https://github.com/acammarat/phtools/tree/main/qpoints)
2. [Fourier transform of anharmonic force constants](#fourier-transform-of-anharmonic-force-constants): generate *phi.nd* with [**phind**](https://github.com/acammarat/pindol/tree/main/phind)
3. [Normal Dynamics simulation](#normal-dynamics-simulation): generate the ND trajectory with [**pindol**](https://github.com/acammarat/pindol/tree/main/pindol)
4. [Analyse the results](#analyse-the-results)
   1. Convert the trajectory, velocities and accelerations in XYZ format
   2. Obtain the effective force constants which include the effect of the temperature

At the moment, the geometry reference file must be created as [VASP](https://www.vasp.at) POSCAR format for all the pre- and postprocessing steps; for the ND run, also the [LAMMPS](https://www.lammps.org) format is accepted. To show the workflow, we consider the crystalline Si system in the [si](https://github.com/acammarat/pindol/tree/main/examples/si) folder as an example.
   
## Dynamical matrix diagonalization

The [qpoints](https://github.com/acammarat/phtools/tree/main/qpoints) code is used in this step. The harmonic force constants, eigenvectors and frequencies can be obtained, in principle, with any calculator. At the moment, only [PHONOPY](https://phonopy.github.io/phonopy) is supported; for this reason, the following example will refer to its usage.

First, we have to use phonopy to obtain the FORCE_CONSTANTS file containing the force constants calculated at small displacements, for which the harmonic expansion of the potential energy is supposed to be a good approximation. For this step, the reader is referred to the [PHONOPY](https://phonopy.github.io/phonopy) web site.

Then we need to choose a list of q-points which constitutes the q-set; for example

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

which is the list of unique points contained in a 6x6x6 Mokhorst-Pack mesh obtained with the aid of the [SPGLIB](https://spglib.readthedocs.io/en/latest) library. Due to implementation reasons, the set must be ordered as G-H-S and must not contain complex-conjugated couples. To this aim, we can use the [q4phind](https://github.com/acammarat/pindol/tree/main/q4phind) code. If we image to write the above list into the file *qlist.dat*, we can then execute

```
$ q4phind qlist.dat
```
which produces the *qordered.dat* file. Since we include the &Gamma; point, before proceding we need to identify the acoustic modes to be excluded from the phonon-phonon scattering tensor to be calculated in the [next step](#fourier-transform-of-anharmonic-force-constants). To this aim, for example, we can visualize the eigenvectors pattern at &Gamma; by generating an .ascii file and visualize it with v_sim (see ["How to watch animation"](https://phonopy.github.io/phonopy/animation.html#how-to-watch-animation) in the phonopy website.). Let us imagine that the modes labeled as 1, 2 and 3 at &Gamma; are the acoustic ones; then, the list of (q,j) modes to skip is "1 1 1 2 1 3", where the first "1" of each couple refers to the first q-point in the *qordered.dat* list (i.e. &Gamma;) and the second number is the *j* label of the mode.

We are then ready to prepare the *qpt.inp* file for **qpoints**:

```
phonopy         # phonopy executable
POSCAR          # reference geometry file
5 5 5           # mesh used to generate the FORCE_CONSTANTS file
3               # nskip number of modes to skip
1 1  1 2  1 3   # (q,j) list of modes to skip
16
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

where from the 6th line on we pasted the content of the *qordered.dat* file. The first line must contain the path to the phonopy executable including all the necessary options to run it; in this example, phonopy is found in one of the folders listed in the PATH variable. Be executing

```
$ qpoints qpt.inp
```
**qpoints** calls **phonopy** to generate the *qpoints.yaml* file and generates the *qmatrix.nd* and *freq.nd* files, which contain the eigenvectors and eigenfrequencies at each q-point specified in the set.

## Fourier transform of anharmonic force constants

Now we need to calculate the Fourier transform of the anharmonic force constants. The [**phind**](https://github.com/acammarat/pindol/tree/main/phind) code is used in this step. The Normal Dynamics formalism accounts for the anharmonic interactions at any order of the Taylor expansion of the potential energy. At the moment, we implemented only the third-order in **pindol**; accordingly, **phind** is written to calculate only the Fourier transform of the third-order Cartesian force constants.

Several calculators are available to obtain the third-order Cartesian force constants; at the moment, **phind** is interfaced only with [PHONO3PY](https://phonopy.github.io/phono3py). **phono3py** is used to obtain the fc3.hdf5 file containing the Cartesian third-order force constants calculated at small displacements, for which the first-order anharmonic term of the potential energy is supposed to be a good approximation. For this step, the reader is referred to the [PHONO3PY](https://phonopy.github.io/phono3py) web site.

## Normal Dynamics simulation

## Analyse the results


...

## Contributions, bug reports and feature requests

We are happy to accept contributions. To report bugs or request new features, please use the [Issue Tracker](https://github.com/acammarat/pindol/issues). If you use the programs in this repository in your work, please send an email to cammaant [at] fel.cvut.cz - we will collect them and put up a list of outputs.

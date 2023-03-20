# Workflow of Normal Dynamics simulations with PINDOL

The main steps to run and analyse a Normal Dynamics (ND) simulation with **pindol** are the following:

1. [Dynamical matrix diagonalization](#dynamical-matrix-diagonalization): generate *qmatrix.nd* and *freq.nd* with [QPOINTS](https://github.com/acammarat/phtools/tree/main/qpoints)
2. [Fourier transform of anharmonic force constants](#fourier-transform-of-anharmonic-force-constants): generate *phi.nd* with [**phind**](https://github.com/acammarat/pindol/tree/main/phind)
3. [Normal Dynamics simulation](#normal-dynamics-simulation): generate the ND trajectory with [**pindol**](https://github.com/acammarat/pindol/tree/main/pindol)
4. [Analyse the results](#analyse-the-results)
   1. Convert the trajectory, velocities and accelerations in XYZ format with [nd2xyz](https://github.com/acammarat/pindol/tree/main/nd2xyz)
   2. Obtain the effective force constants which include the effect of the temperature with [nd2hiphive](https://github.com/acammarat/pindol/tree/main/nd2hiphive) and [hiPhive](https://hiphive.materialsmodeling.org)

At the moment, the geometry reference file must be created as [VASP](https://www.vasp.at) POSCAR format for all the pre- and postprocessing steps; for the ND run, also the [LAMMPS](https://www.lammps.org) format is accepted. To show the workflow, we consider the crystalline Si system in the [si](https://github.com/acammarat/pindol/tree/main/examples/si) folder as an example.

## Dynamical matrix diagonalization

The [qpoints](https://github.com/acammarat/phtools/tree/main/qpoints) code is used in this step. The harmonic force constants, eigenvectors and frequencies can be obtained, in principle, with any calculator. At the moment, only [PHONOPY](https://phonopy.github.io/phonopy) is supported; for this reason, the following example will refer to its usage.

First, we have to use phonopy to obtain the FORCE_CONSTANTS file containing the force constants calculated at small displacements, for which the harmonic expansion of the potential energy is supposed to be a good approximation. For this step, the reader is referred to the [PHONOPY](https://phonopy.github.io/phonopy) web site. The units of the FORCE_CONSTANTS must be eV/Ang<sup>2</sup>.

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

where from the 6th line on we pasted the content of the *qordered.dat* file. The first line must contain the path to the phonopy executable including all the necessary options to run it; in this example, phonopy is found in one of the folders listed in the PATH variable. By executing

```
$ qpoints qpt.inp
```
**qpoints** calls **phonopy** to generate the *qpoints.yaml* file and generates the *qmatrix.nd* and *freq.nd* files, which contain the eigenvectors and eigenfrequencies at each q-point specified in the set. The eigenvectors are dimensionless while the frequencies are in THz, as in the **phonopy** output.

## Fourier transform of anharmonic force constants

Now we need to calculate the Fourier transform of the anharmonic force constants. The [**phind**](https://github.com/acammarat/pindol/tree/main/phind) code is used in this step. The Normal Dynamics formalism accounts for the anharmonic interactions at any order of the Taylor expansion of the potential energy. At the moment, we implemented only the third-order in **pindol**; accordingly, **phind** is written to calculate the Fourier transform of only the third-order Cartesian force constants.

Several calculators are available to obtain the third-order Cartesian force constants; at the moment, **phind** is interfaced only with [PHONO3PY](https://phonopy.github.io/phono3py). **phono3py** is used to obtain the *fc3.hdf5* file containing the Cartesian third-order force constants calculated at small displacements, for which the first-order anharmonic term of the potential energy is supposed to be a good approximation. For this step, the reader is referred to the [PHONO3PY](https://phonopy.github.io/phono3py) web site. The units of the force constants in *fc3.hdf5* must be eV/Ang<sup>3</sup>.

Then, we use the *fc3_extract.py* python script in the [examples](https://github.com/acammarat/pindol/tree/main/examples) folder, which reads *fc3.hdf5* and writes the *fc3.dat* input file for **phind**:

```
$ python fc3_extract.py
```
The units of the force constants in *fc3.dat* are the same as in *fc3.hdf5*.
We now create the input file *phind.inp* for **phind**:

```
POSCAR         # reference geometry
5 5 5          # supercell dimensions used to generate the fc3.hdf5 file
1              # natom_types
Si 28.0855     # atomic symbol, mass [amu]
0              # write (q,j -> l) map: 0=no, 1=yes
1              # calc phi: 0=no, 1=yes
```
Since **phind** is able to exploit the OpenMP parallelization, it is highly recommended to set the OMP_NUM_THREADS variable:
```
$ export OMP_NUM_THREADS=64
```
where 64 must be substituted with a suitable number according to the number of cores available on your machine. By running
```
$ phind phind.inp
```
**phind** produces the file *phi.nd* which contains the Fourier transform of the Cartesian tensor of the third-order force constants.

## Normal Dynamics simulation

Once we have *qmatrix.nd*, *freq.nd* and *phi.nd*, we are almost ready to run **pindol** to perform our Normal Dynamics simulation. The reader is referred to the [pindol](https://github.com/acammarat/pindol/tree/main/pindol) README.md file for the possible options. Let us imagine that we want to run an NVT simulation at 10K, for 4 ns with 1 fs as time step; the *pindol.inp* looks like

```
ATOMS
  types 1
  atom Si 28.0855
END_ATOMS

REFCONF poscar POSCAR

UNITS real

INITVEL 10.0 666 # temperature, seed

ND
  nvt 10.0 10.0 100.0 # initial temperature, final temperature, characteristic time scale of the thermostat 
  integrator 1.0 22 both 1.d-15 # dt, method, accuracy
  runsteps 4000000
  printcoord 1 normcoord.dat
  printvel 1 normvel.dat
  printacc 1 normacc.dat
END_ND

FINALCONF poscar final.vasp

WRITERESTART final.restart

END
```
A very basic OpenMP parallelization is implemented in **pindol**, which can be enabled by exporting the OMP_NUM_THREADS variable as done above for **phind**. We recommend to test the execution time over different values of the OMP_NUM_THREADS variable on a short run, in order to estimate in advance the optimal number of cores. The ND simulation is then performed by executing
```
$ pindol | tee out
```
where we choose to write the standard output also to the *out* file for convenience. The beginning of the output looks like this:

```
        _           _       _  
  _ __ (_)_ __   __| | ___ | | 
 | '_ \| | '_ \ / _` |/ _ \| | 
 | |_) | | | | | (_| | (_) | | 
 | .__/|_|_| |_|\__,_|\___/|_| 
 |_|                           v1.0  

Running on 1 threads
Reading "pindol.inp"...
  1 atomic type(s) selected: Si
  'POSCAR' poscar reference configuration file will be read
  real units selected
  atomic velocities will be initialized at      10.0000 K
  NVT dynamics selected at      10.0000 K with a Nose-Hoover thermostat and a characteristic timescale of     100.0000 fs
  'final.vasp' poscar final configuration file will be written
  'final.restart'  restart file will be written
Reading "pindol.inp"...DONE.
Reading reference configuration...
  2 atoms found in the reference cell
Reading reference configuration...DONE.
Reading polarization vectors...
  1 q-point found: (0,0,0)
  2 q-point found: (1/2,0,0)
  3 q-point found: (1/2,1/6,0)
  4 q-point found: (1/2,1/3,0)
  5 q-point found: (1/2,1/2,0)
  6 q-point found: (1/2,1/3,1/6)
  7 q-point found: (-1/3,1/2,1/6)
  8 q-point found: (1/6,0,0)
  9 q-point found: (1/3,0,0)
  10 q-point found: (1/6,1/6,0)
  11 q-point found: (1/3,1/6,0)
  12 q-point found: (-1/3,1/6,0)
  13 q-point found: (-1/6,1/6,0)
  14 q-point found: (1/3,1/3,0)
  15 q-point found: (-1/3,1/3,0)
  16 q-point found: (-1/3,1/3,1/6)
  Supercell size: 6 x 6 x 6
Reading polarization vectors...DONE.
Reading matrix of the phonon interaction strength...
  5389 phi elements read
Reading matrix of the phonon interaction strength...DONE.
Setting up the supercell...
  432 atoms in the supercell
Setting up the supercell...DONE.
Atomic velocities initialized at      10.0000 K.
Performing normal dynamics...
#             time    harm. potential energy anh. potential energy      potential_energy        kinetic_energy          total_energy    conserved_quantity   current_temperature
   0.0000000000000       0.0000000000000       0.0000000000000       0.0000000000000       1.5216794479993       1.5216794479993       1.5216794479993       10.418215264719    
   1.0000000000000      0.79488050631480E-02  0.17421750022140E-08  0.79488068053230E-02   1.5137244355497       1.5216732423551       1.5216794517577       10.363751078953    

...
```
At the end of the execution, the output looks like

```
...
   4000000.0000000      0.89127593451890     -0.27484485814697E-05  0.89127318607031      0.60773675481256       1.4990099408829       1.5212621755371       4.1608844387324    
Performing normal dynamics...DONE.
Writing the final configuration...DONE.
Writing the restart file...DONE.
Deallocating variables...DONE.
Total run time is: 4 hours, 25 minutes, 11 seconds and 505 milliseconds.
```
The files *normcoord.dat*, *normvel.dat* and *normacc.dat* specified in the input are written incrementally during the run at the frequency specified in the corresponding keyword. Other files (e.g., final configuration, restart) are written as per user request.


## Analyse the results

### Conversion of normal trajectory, velocities and accelerations in standard Cartesian .xyz format

We assume that in the *pindol.inp* input file, we asked to create the files containing the normal coordinates, velocities and accelerations as *normcoord.dat*, *normvel.dat* and *normacc.dat*, respectively. In order to convert such files into Cartesian .xyz format, we will use [**nd2xyz**](https://github.com/acammarat/pindol/tree/main/nd2xyz). We then create the *nd2xyz.inp* input file which looks like

```
POSCAR                     # reference geometry
normcoord.dat              # ND trajectory ; "0" = no file to convert
normvel.dat                # qdot file     ; "0" = no file to convert
normacc.dat                # qddot file    ; "0" = no file to convert
1                          # natom_types
Si 28.0855
3000000 4000000            # time window: initial and final time
0                          # "1" = use the last step of the input files
                           # regardless the time window specified above
                           # and write a POSCAR_last.vasp configuration

In this example, we want to extract the last nanosecond of trajectory. The Cartesian trajectory (Ang), velocities (Ang/fs) and accelerations (Ang/fs<sup>2</sup>) are written in *ndtrj.xyz*, *ndvel.xyz* and *ndacc.xyz*, respectively. The file *ndtraj_ave.xyz* contains the average atomic position within the specified time window. Such files can then be used for postprocessing to extract any physical quantity as in usual standard Molecular Dynamics (MD) simulations.

### Conversion of normal trajectory, velocities and accelerations in standard Cartesian .xyz format

## Contributions, bug reports and feature requests

We are happy to accept contributions. To report bugs or request new features, please use the [Issue Tracker](https://github.com/acammarat/pindol/issues). If you use the programs in this repository in your work, please send an email to cammaant [at] fel.cvut.cz - we will collect them and put up a list of outputs.

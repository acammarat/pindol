# pindol
Phonon-Inspired Normal Dynamics of Lattices

**pindol** is a package to performs atom dynamics in the NVE and NVT ensembles. It exploits the Normal Dynamics formalism, allowing to perform Ab Initio Molecular Dynamics at the cost of Molecular Dynamics.

The **pindol** package contains the following codes:

- [**q4phind**](https://github.com/acammarat/pindol/tree/main/q4phind) (Preprocessing) Orders a q-point set as GM, H, S, used to generate qmatrix.nd and freq.nd with [QPOINTS](https://github.com/acammarat/phtools/tree/main/qpoints) compatible with [PHIND](https://github.com/acammarat/pindol/tree/main/phind)
- [**phind**](https://github.com/acammarat/pindol/tree/main/phind) (Preprocessing) Calculates the first order anharmonic interaction strength &Phi; to be used for the [PINDOL](https://github.com/acammarat/pindol/tree/main/pindol) code
- [**pindol**](https://github.com/acammarat/pindol/tree/main/pindol) Main code to perform Normal Dynamics simulations
- [**nd2xyz**](https://github.com/acammarat/pindol/tree/main/nd2xyz) (Postprocessing) Converts normal trajectories into XYZ format
- [**nd2hiphive**](https://github.com/acammarat/pindol/tree/main/nd2hiphive) (Postprocessing) Creates input files for [hiPhive](https://hiphive.materialsmodeling.org/) to extract the effective force constants

The preprocessing tool [QPOINTS](https://github.com/acammarat/phtools/tree/main/qpoints) must be used to prepare the qmatrix.nd and freq.nd files.

## Installation

The codes require a fortran compiler (tested with gfortran 11.3.0). After cloning, enter each folder and type

`make`

If the compilation ends successfully, the executable with the name of the code is created.

## Usage

For the usage, please refer to the related README.md file.

## Citation

The users of PINDOL have little formal obligations specified in the [GNU General Public License](http://www.gnu.org/copyleft/gpl.txt).
However, it is common practice in the scientific literature, to acknowledge the efforts of people that have made the research possible.
In this spirit, please cite

A. Cammarata, M. Dasic and P. Nicolini, *Normal Dynamics: solving Newton’s equations in the reciprocal space*, Phys. Rev. Lett **XX**, XXXXX (XXXX) DOI: [xxx](https://doi.org/10.1103/xxx)

A. Cammarata, M. Dasic and P. Nicolini, *Sampling dynamical trajectories in the reciprocal space*, Phys. Rev. B **XX**, XXXXX (XXXX) DOI: [xxx](https://doi.org/10.1103/xxx)

where the formulation used to perform normal dynamics is reported.

## Contributions, bug reports and feature requests

We are happy to accept contributions. To report bugs or request new features, please use the [Issue Tracker](https://github.com/acammarat/pindol/issues). If you use the programs in this repository in your work, please send an email to cammaant [at] fel.cvut.cz - we will collect them and put up a list of outputs.

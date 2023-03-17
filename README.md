# pindol
Phonon-Inspired Normal Dynamics of Lattices

**pindol** is a package to performs atom dynamics in the NVE and NVT ensembles. It exploits the Normal Dynamics formalism, allowing to perform Ab Initio Molecular Dynamics at the cost of Molecular Dynamics.

The **pindol** package contains the following codes:

- [**q4phind**](https://github.com/acammarat/pindol/tree/main/q4phind) Order a q-point set as GM, H, S, used to generate qmatrix.nd and freq.nd with [QPOINTS](https://github.com/acammarat/phtools/tree/main/qpoints) compatible with [PHIND](https://github.com/acammarat/pindol/tree/main/phind)

## Installation

The code requires a fortran compiler (tested with gfortran 11.3.0). After cloning, from this folder it is enough to type

`make`

or you can type

`make pindol`

from the parent directory.

For the integration of the equation of motion, the code makes use of the [DVODE](https://computing.llnl.gov/sites/default/files/dvode.f) subroutine.

If the compilation ends successfully, the executable pindol is created.

## Usage

The program requires a few inputs in order to be run:

- the reference configuration (either in POSCAR or LAMMPS data format);

- two files containing the phonon eigenvector and eigenvalues (frequencies), which can be generated using the auxiliary code PHXXX; 

- the file containing the third-order coupling constants, which can be generated using the auxiliary code PHXXX;

- the input scipt (see the `template.inp` file in this folder for more information about the syntax and the available keywords). 

To execute the code, it is sufficient to run

`pindol`

## Example

XXX {Magari possiamo metterci lo zirconio? (il template.inp e' per Zr...)}

## Citation

 The users of PINDOL have little formal obligations specified in the [GNU General Public License](https://www.gnu.org/licenses/old-licenses/gpl-2.0.txt).
 However, it is common practice in the scientific literature, to acknowledge the efforts of people that have made the research possible.
 In this spirit, please cite

 A. Cammarata, P. Nicolini, M. Dasic, T. Polcar, *XXX*, Phys. Rev. XXX **XXX**, XXX (XXX), DOI: [XXX](XXX)

 where the formulation used to perform normal dynamics is reported.


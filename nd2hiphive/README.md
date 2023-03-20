# nd2hiphive

Processes the normal coordinate and acceleration files to produce the necessary input for [hiPhive](https://hiphive.materialsmodeling.org)

## Installation

The code requires a fortran compiler and the [lapack](https://netlib.org/lapack/) libraries. After cloning, enter the folder and compile it with

`make`

If the compilation ends successfully, the executable nd2hiphive is created.

## Usage

The qmatrix.nd and freq.nd files of the reference geometry must be present in the folder where the code is executed.

The format of the input file is


```

string                  name of the POSCAR file
string                  name of the ND trajectory file
string                  name of the ND Qddot file
int                     number of atomic types (natom_types)
string double           atomic symbol and mass (at_pertype(i), mass_pertype(i) [uma] of atom type 1
...  ...
string double           atomic symbol and mass (at_pertype(i), mass_pertype(i) [uma] of atom type natom_types
double double double    initial time, final time, max skip time


```

where `int` is an integer number, `double` is a real number of type double, and `string` a string of type char. The command line syntax can be shown by using the `-h` option:

```

$ nd2hiphive
            _ ____  _     _       _     _           
  _ __   __| |___ \| |__ (_)_ __ | |__ (_)_   _____ 
 | '_ \ / _` | __) | '_ \| | '_ \| '_ \| \ \ / / _ \
 | | | | (_| |/ __/| | | | | |_) | | | | |\ V /  __/
 |_| |_|\__,_|_____|_| |_|_| .__/|_| |_|_| \_/ \___|
                           |_|                      
                                            0.7

 Syntax: nd2hiPhive <setting file>

```

After the execution, the following files are created:

- ndhiPhive_prim.xyz the reference unit cell
- ndhiPhive_prim_dir.vasp the reference unit cell in POSCAR format
- ndhiPhive_superc.xyz the supercell commensurate with the q-point set
- ndhiPhive_phonopy.py the python code to be added to the example fcs.py script, the latter used to call hiPhive and calculate the effective force constants

## Example

The *nd2hiphive.inp* file is an example of input file. The *fcs.py* file is an example of python script to execute hiPhive

## Citation

The users of ND2XYZ have little formal obligations specified in the [GNU General Public License](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html).
However, it is common practice in the scientific literature, to acknowledge the efforts of people that have made the research possible.
In this spirit, please cite

A. Cammarata, M. Dasic and P. Nicolini, *Normal Dynamics: solving Newton’s equations in the reciprocal space*, Phys. Rev. Lett **XX**, XXXXX (XXXX) DOI: [xxx](https://doi.org/10.1103/xxx)

A. Cammarata, M. Dasic and P. Nicolini, *Sampling dynamical trajectories in the reciprocal space*, Phys. Rev. B **XX**, XXXXX (XXXX) DOI: [xxx](https://doi.org/10.1103/xxx)


# nd2hiphive

Processes the normal coordinate and acceleration files to produce the necessary input for [hiPhive](https://hiphive.materialsmodeling.org).

## Installation

The code requires a fortran compiler and the [lapack](https://netlib.org/lapack/) libraries. After cloning, enter the folder and compile it with

`make`

If the compilation ends successfully, the executable nd2hiphive is created.

## Usage

The *qmatrix.nd* and *freq.nd* files of the reference geometry must be present in the folder where the code is executed.

The format of the input file is


```

string                  name of the POSCAR file
string                  name of the normal coordinates file
string                  name of the normal accelerations file
int                     number of atomic types 
string double           atomic symbol and mass (amu) of the first atom type
...  ...
string double           atomic symbol and mass (amu) of the last atom type
double double double    initial time, final time, max skip time


```

where `int` is an integer number, `double` is a real number of type double, and `string` a string of type char. The *skip time* is meant to skip configurations falling within a time window with maximum random width equal to skip_time; this is to avoid to collect very similar configurations and forces which do not improve the force constants fit in a significant way. The command line syntax can be shown by using the `-h` option:

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

- *ndhiPhive_prim.xyz* the reference unit cell
- *ndhiPhive_prim_dir.vasp* the reference unit cell in POSCAR format
- *ndhiPhive_superc.xyz* the supercell commensurate with the q-point set
- *ndhiPhive_phonopy.py* the python code to be added to the example fcs.py script, the latter used to call **hiPhive** and calculate the effective force constants
- *ndhiPhive_dispfor.xyz* containing displacements (Ang) and forces (eV/Ang) at the extracted trajectory snapshots.

## Example

The *nd2hiphive.inp* file is an example of input file. The *fcs.py* file is an example of python script to execute hiPhive

## Citation

The users of **nd2hiphive** have little formal obligations specified in the [GNU General Public License](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html).
However, it is common practice in the scientific literature, to acknowledge the efforts of people that have made the research possible.
In this spirit, please cite

A. Cammarata, M. Dasic and P. Nicolini, *Normal Dynamics: solving Newton’s equations of motion in the reciprocal space*, DOI:

<!--- [https://dx.doi.org/10.2139/ssrn.4550608](https://dx.doi.org/10.2139/ssrn.4550608) --->


# nd2xyz

Converts normal trajectories into XYZ format.

## Installation

The code requires a fortran compiler. After cloning, enter the folder and compile it with

`make`

If the compilation ends successfully, the executable nd2xyz is created.

## Usage

The *qmatrix.nd* and *freq.nd* files of the reference geometry must be present in the folder where the code is executed.

The format of the input file is


```

string                  name of the POSCAR file
string                  name of the normal coordinate file; 0 = no file to convert
string                  name of the normal velocities file; 0 = no file to convert
string                  name of the normal accelerations file; 0 = no file to convert
int                     number of atomic types
string double           atomic symbol and mass (amu) of the first atom type
...  ...
string double           atomic symbol and mass (amu) of the last atom type
double double           time window: initial and final time
flag                    1 = write POSCAR restart from last configuration

```

where `int` is an integer number, `double` is a real number of type double, and `string` a string of type char. The command line syntax can be shown by using the `-h` option:

```

$ nd2xyz 
            _ ____                 
  _ __   __| |___ \__  ___   _ ____
 | '_ \ / _` | __) \ \/ / | | |_  /
 | | | | (_| |/ __/ >  <| |_| |/ / 
 |_| |_|\__,_|_____/_/\_\\__, /___|
                         |___/     
                           0.7

  Syntax: nd2xyz <setting file>

```

After the execution the following files are created if requested by the user:

- *ndtrj.xyz* the Cartesian trajectory (Ang)
- *ndvel.xyz* the Cartesian velocities (Ang/fs)
- *ndacc.xyz* the Cartesian accelerations (Ang/fs<sup>2</sup>)
- *ndtraj_ave.xyz* the average atomic position within the specified time window (Ang)

## Example

The *nd2xyz.inp* file is an example of input file. 

## Citation

The users of **nd2xyz** have little formal obligations specified in the [GNU General Public License](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html).
However, it is common practice in the scientific literature, to acknowledge the efforts of people that have made the research possible.
In this spirit, please cite

A. Cammarata, M. Dasic and P. Nicolini, *Normal Dynamics: solving Newton’s equations in the reciprocal space*, Phys. Rev. Lett **XX**, XXXXX (XXXX) DOI: [xxx](https://doi.org/10.1103/xxx)

A. Cammarata, M. Dasic and P. Nicolini, *Sampling dynamical trajectories in the reciprocal space*, Phys. Rev. B **XX**, XXXXX (XXXX) DOI: [xxx](https://doi.org/10.1103/xxx)

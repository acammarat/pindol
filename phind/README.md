# phind

Calculates the first order anharmonic interaction strength &Phi; to be used for the [PINDOL](https://github.com/acammarat/pindol/tree/main/pindol) code

## Installation

The code requires a fortran compiler. After cloning, enter the folder and compile it with

`make`

If the compilation ends successfully, the executable phind is created.

## Usage

The format of the input file is


```
string                  name of the POSCAR file
int int int           dimensions of the supercell used to create fc3.dat (ncells(i))
int                   number of atomic types (natom_types)
string double           atomic symbol and mass (at_pertype(i), mass_pertype(i) [uma] of atom type 1
...  ...
string double           atomic symbol and mass (at_pertype(i), mass_pertype(i) [uma] of atom type natom_types
int                   write (q,j -> l) map: 0=no, 1=yes
int                   calc phi: 0=no, 1=yes


```

where `int` and `double` are an integer and a real number of type double, respectively, while `string` is a string of type char. The command line syntax can be shown by using the `-h` option:

```

$ phind
      _                   
    _| |_              _  
   /     \   _ __   __| | 
  ( (| |) ) | '_ \ / _` | 
   \_   _/  | | | | (_| | 
     |_|    |_| |_|\__,_|  
                  4.3

 Syntax: phind <setting file>
 Export the number of threads before executing phind:
 export OMP_NUM_THREADS=

```

After the execution, the file phi.nd is created, to be used with the main [pindol](https://github.com/acammarat/pindol/tree/main/pindol) code. The output file phi.qallowed.nd contains the (q,q',q'') triplets which satisfy the &Delta;(q+q'+q'') selection rule.


## Example

The *phind.inp* file is an example of input file. 

## Citation

The users of PHIND have little formal obligations specified in the [GNU General Public License](http://www.gnu.org/copyleft/gpl.txt).
However, it is common practice in the scientific literature, to acknowledge the efforts of people that have made the research possible.
In this spirit, please cite

A. Cammarata, M. Dasic and P. Nicolini, *Normal Dynamics: solving Newton’s equations in the reciprocal space*, Phys. Rev. Lett **XX**, XXXXX (XXXX) DOI: [xxx](https://doi.org/10.1103/xxx)

A. Cammarata, M. Dasic and P. Nicolini, *Sampling dynamical trajectories in the reciprocal space*, Phys. Rev. B **XX**, XXXXX (XXXX) DOI: [xxx](https://doi.org/10.1103/xxx)

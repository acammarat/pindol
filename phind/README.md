# phind

Calculates the first order anharmonic interaction strength &Phi; to be used by the [PINDOL](https://github.com/acammarat/pindol/tree/main/pindol) code.

## Installation

The code requires a fortran compiler. After cloning, enter the folder and compile it with

`make`

If the compilation ends successfully, the executable **phind** is created.

## Usage

The format of the input file is


```
string                name of the POSCAR file
int int int           dimensions of the supercell used to create fc3.dat
int                   number of atomic types
string double         atomic symbol and mass (amu) of the first atom type
...  ...
string double         atomic symbol and mass (amu) of the last atom type
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

After the execution, the file *phi.nd* is created, to be used with the main [pindol](https://github.com/acammarat/pindol/tree/main/pindol) code. The output file phi.qallowed.nd contains the (q,q',q'') triplets which satisfy the &Delta;(q+q'+q'') selection rule.


## Example

The *phind.inp* file is an example of input file. 

## Citation

The users of **phind** have little formal obligations specified in the [GNU General Public License](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html).
However, it is common practice in the scientific literature, to acknowledge the efforts of people that have made the research possible.
In this spirit, please cite

A. Cammarata, M. Dasic and P. Nicolini, *Normal Dynamics: solving Newton’s equations of motion in the reciprocal space*, DOI:
<!--- [https://dx.doi.org/10.2139/ssrn.4550608](https://dx.doi.org/10.2139/ssrn.4550608) --->

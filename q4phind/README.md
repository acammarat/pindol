# q4phind

Orders a q-point set as GM, H, S, used to generate *qmatrix.nd* and *freq.nd* with [QPOINTS](https://github.com/acammarat/phtools/tree/main/qpoints) and compatible with [PHIND](https://github.com/acammarat/pindol/tree/main/phind)

## Installation

The code requires a fortran compiler. After cloning, enter the folder and compile it with

`make`

If the compilation ends successfully, the executable q4phind is created.

## Usage

The format of the input file is


```

double double double  reduced components of the first q-point
...
double double double  reduced components of the last q-point

```

where `double` is a real number of type double. The command line syntax can be shown by using the `-h` option:

```

$ q4phind 
       _  _       _                   
  __ _| || |    _| |_              _  
 / _` | || |_  /     \   _ __   __| | 
| (_| |__   _|( (| |) ) | '_ \ / _` | 
 \__, |  |_|   \_   _/  | | | | (_| | 
    |_|          |_|    |_| |_|\__,_| 
                              0.1

 Syntax: q4phind <q-point set file>

```

After the execution the file qordered.nd is created. The content can be pasted into the input file of [QPOINTS](https://github.com/acammarat/phtools/tree/main/qpoints) at the point where the list of q-points is specified.

## Citation

The users of **q4phind** have little formal obligations specified in the [GNU General Public License](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html).
However, it is common practice in the scientific literature, to acknowledge the efforts of people that have made the research possible.
In this spirit, please cite

A. Cammarata, M. Dasic and P. Nicolini, *Normal Dynamics: solving Newton’s equations of motion in the reciprocal space*, DOI:
<!--- [https://dx.doi.org/10.2139/ssrn.4550608](https://dx.doi.org/10.2139/ssrn.4550608) --->

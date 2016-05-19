# Parthess User Guide


## Overview

*Parthess* is a program to determine force constants in Molecular Mechanics (MM) by *Partial Hessian Fitting (PHF)* scheme. The details of *PHF* scheme is published on XXX.


## Quick Start
1. Install [Anaconda(Python3 version)](https://www.continuum.io/downloads) 
2. Install [cclib](https://cclib.github.io/). Just run ```pip install cclib```.
3. Download [rxcclib](https://github.com/ruixingw/rxcclib).
4. Setting up the $PYTHONPATH environment and put *rxcclib* in it.
5. Download [Tsubasa](https://github.com/ruixingw/tsubasa).
6. Prepare molecule geometry, and run *Tsubasa* program.
7. From the output of *Tsubasa* program, run *Parthess*.


## HowTo

This program is written in Python3 and used many third-party packages.

Dependencies:
- Python3
- numpy > 1.10
- [cclib](https://cclib.github.io/) > 1.3.2
- [rxcclib](https://github.com/ruixingw/rxcclib)


If you are not an experienced Python programmer, it is **strongly recommended** to use [Anaconda](https://www.continuum.io/downloads) rather than install Python and packages separately. 

*numpy* is already included in the newest version of *Anaconda* with [Intel MKL Optimizations](https://www.continuum.io/blog/developer-blog/anaconda-25-release-now-mkl-optimizations), which gives a much better performance.

*Parthess* requires another code package written by us, temporarily named [*rxcclib*](https://github.com/ruixingw/rxcclib). It defines molecules, MM functions and methods for file parsing. You will need to install it manually.

[*cclib*](https://cclib.github.io/) is not included in Anaconda, you will need to install it manually. Installation of cclib is very easy:

```pip install cclib```

Although, only a small unit conversion utility is used in *cclib*. We plan to contribute some file parsing codes of *rxcclib* to *cclib* in the future.

## Prepare Input File

It is recommended to generate all inputs file by [*Tsubasa* program](https://github.com/ruixingw/tsubasa).

Parthess scheme needs an MM input file and the QM Hessian input file.

- *Parthess* accepts the Gaussian Formatted Checkpoint File (FCHK) as the QM Hessian input. A frequency calculation can yield the Hessian.
- The MM input file follows Gaussian Format, which includes the optimized geometry, atom types, atomic charges, and MM functions. The undetermined force constants are denoted by "XXXXXX". An example for H<sub>2</sub>O<sub>2</sub> is attached.


```
%mem=12gb         
#p amber=softonly geom=connectivity nosymm
iop(4/33=3,7/33=1)
freq

MM-H2O2

0 1
O-oh--0.410780   -3.932366558742   -2.706592019388   -0.283087648364
H-ho-0.410780   -4.104122594899   -1.885283961559    0.212183800621
O-oh--0.410780   -2.475774075762   -2.710437185030   -0.283043419203
H-ho-0.410780   -2.304047969868   -3.531709623942    0.212297456772

 1 2 1.0 3 1.0
 2
 3 4 1.0
 4

AmbTrs ho oh oh ho 0 0 0 0 0.0 XXXXXX 0.0 0.0 1.0
HrmBnd1 ho oh oh XXXXXX 100.2817
HrmStr1 ho oh XXXXXX 0.97434
HrmStr1 oh oh XXXXXX 1.45660
Nonbon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 -1.2
VDW oh  1.7210  0.2104
VDW ho  0.0000  0.0000

```




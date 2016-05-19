# Parthess


## Overview
Parthess is a program to determine force constants in Molecular Mechanics (MM) by *Partial Hessian Fitting* (PHF) scheme. The details of PHF scheme is published on XXX.


## HowTo
This program is written in Python3 and used many third-party packages.
Dependencies:
- Python3 , while only Python 3.5.1 (Anaconda 4.0.0) is tested.
- numpy >1.10
- cclib >1.3.2
- rxcclib (See this page)

It is strongly recommended to use Anaconda rather than install Python and packages separately, if you are not an experienced Python programmer.

We also wrote another library temporarily named *rxcclib*, which defines molecules, MM functions and methods for file parsing. You will need to install it manually.

*cclib* is not included in Anaconda, you will need to install it manually. We plan to contribute some codes of rxcclib to cclib in the future.
Installation of cclib is very easy:

```pip install cclib```


## Prepare Input File
PHF scheme will need a MM input file and the QM Hessian. 

- QM Hessian is usually obtained by a frequency calculation. *Parthess* program accepts the Gaussian Formatted Checkpoint File (FCHK) as the QM Hessian input.
- MM input file includes the equilibrium geometry (same with FCHK file), atom types, atomic charges, and MM functions. The MM input file follows Gaussian Format. An example of H_2O_2 is attached.

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




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

It is strongly recommended to use Anaconda rather than install Python and packages separately.

We also wrote another library temporarily named *rxcclib*, which defines molecules, MM functions and methods for file parsing. 
cclib is not included in Anaconda, you will need to install it first. We plan to contribute some codes of rxcclib to cclib in the future.
Installation of cclib is very easy:
```pip install cclib```




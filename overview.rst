Overview
========
**Parmhess** is a program to determine force field constants for AMBER force field, written in Python_, based on Hessian Fitting methods. It can:

.. _Python: http://www.python.org

- Use *Cartesian* Hessian matrices to determine force constants.
- Derive force constants for quadratic terms: bond-stretching and angle-bending.
- Derive force constants for one-Fourier-term dihedal-torsion and AMBER-type improper-torsion.


*Contact*

- Ruixing Wang: rwang013[%]e.ntu.edu.sg

- Dr.Hajime Hirao: hirao[%]ntu.edu.sg


Or, open an issue in the github page:

- Issue Tracker: https://github.com/ruixingw/parmhess/issues
- Source Code: https://github.com/ruixingw/parmhess


About Parmhess
==============

The code behind was developed in Dr.Hirao's group and is licensed under the LGPL_. The method, tests, and benchmarks are published on:

|    Wang, R., Ozhgibesov, M. and Hirao, H., Partial hessian fitting for determining force constant parameters in molecular mechanics. *J. Comput. Chem.*, 37(26), pp.2349-2359, **2016** (DOI: `10.1002/jcc.24457`_)  (Featured as Insider Cover, DOI: `10.1002/jcc.24485`_)
|

.. _LGPL: http://www.gnu.org/copyleft/lgpl.html
.. _`10.1002/jcc.24457`: http://dx.doi.org/10.1002/jcc.24457
.. _`10.1002/jcc.24485`: http://dx.doi.org/10.1002/jcc.24485

If you use Parmhess in your scientific work, please support our work by citing this article.


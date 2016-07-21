============
Installation
============

Prerequisites
-------------

1. Python3 (Anaconda)

Parmhess is written in Python3_ and runs on Linux.

The easiest way to setup Python3 is to install **Anaconda**, "the leading open data science platform powered by Python." It is freely available at:

    http://www.continuum.io/downloads

Be sure to install Anaconda for Linux (Python3.5), but not Python2.7 nor for other operating system.

If you are an experienced Python3 user and prefer to use your own Python build, make sure that **numpy** and **pyyaml** be installed properly.

2. Gaussian 09 (G09)

Gaussian 09 is used to perform QM and MM calculation and should be properly installed on your machine. Any revision of G09 should work, although, due to `a bug of Revision B.01`_ (search "Gaussian 09 fix" in this page), some attention should be paid to RESP charge calculation in **Tsubasa**. To use G09 B01, See *Job Control Arguments* :doc:`Tsubasa Manual<tsubasa>`.

For MM calculation (**Parmhess**), :code:`g09` is directly called as it is usually very fast. For QM calculation (**Tsubasa**), you could use your job submit command or script (see Tsubasa manual).

3. AmberTools

AmberTools is used to identify atom types and calculate RESP charge from Gaussian output. It is freely available at `Amber website`_ and the installation guide is in `Amber Manual`_. Make sure it is installed properly as **Tsubasa** will call :code:`antechamber` directly from :code:`$PATH`.

Download Parmhess
-----------------
Download the `the newest release of **Parmhess**`__.

.. __ : https://github.com/ruixingw/parmhess/releases

You will get:

1. File "parmhess.py", "classdef.py". The **Parmhess** program.
2. File "tsubasa.py", "config.yml", "vdw.dat". A copy of the **Tsubasa** program, which is used to prepare the input files.
3. Folder "rxcclib".  A copy of the library code package **rxcclib**, the header file of **Tsubasa** and **Parmhess**.


Make it easier to use
---------------------

You may wish to create soft-links for **parmhess.py** and **tsubasa.py** to your :code:`$PATH`. Suppose :code:`~/bin` is in your :code:`$PATH`, then run:

.. code-block:: bash

    ln -s /PathToParmhess/parmhess.py ~/bin
    ln -s /PathToParmhess/tsubasa/tsubasa.py ~/bin

which allows you to directly run :code:`parmhess.py` and :code:`tsubasa.py` at any directory.


.. _`Amber website` : http://ambermd.org/#AmberTools
.. _`Amber Manual` : http://ambermd.org/doc12/
.. _anaconda : https://www.continuum.io/downloads
.. _Python3: https://www.python.org/
.. _`a bug of Revision B.01` : http://ambermd.org/bugfixesat.html
.. _rxcclib: https://github.com/ruixingw/rxcclib
.. _Tsubasa: https://github.com/ruixingw/tsubasa

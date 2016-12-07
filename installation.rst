============
Installation
============

Prerequisites
-------------

1. Python3 (Anaconda)

Parmhess is written in Python3_ and runs on Linux.

The easiest way to setup Python3 is to install **Anaconda**. It is freely available at:

    http://www.continuum.io/downloads

Be sure to install Anaconda for Linux (Python3.X), not Python2.7 nor for other operating system.

If you are an experienced Python3 user and prefer to use your own Python build, make sure that **numpy** and **pyyaml** be installed properly.

2. Gaussian 09 (G09)

Gaussian 09 is used to perform QM and MM calculation and should be properly installed on your machine.

Any revision of G09 should work, although, there is a bug in G09 Rev. B.01, which cause special attention in RESP charge calculation:


..
  
  Gaussian 09 fix (cited from `AmberTools website`_)
  In Gaussian09 rev B.01, the facility to write out the electrostatic potential on a grid of points was inadvertently deleted. This means that antechamber and resp jobs won't work as they should. Fernando Clemente of Gaussian has kindly provided a script to work around the problem. Download the `fixreadinesp.sh`_ file, and follow the instructions there. (Note: you will have to make the script executable by typing chmod +x fixreadinesp.sh.)

After fixing the Gaussian Output, it can be provided to do next steps. See **Job Control Arguments** section in :doc:`Prepare inputs by Tsubasa<tsubasa>`.

For MM calculation (**Parmhess**), :code:`g09` is directly called as it is usually very fast. For QM calculation (**Tsubasa**), you could use your job submitting command or script for queuing system. See **Config file** in :doc:`Tsubasa manual<tsubasa>`.

3. AmberTools

AmberTools is used to identify atom types and calculate RESP charge from Gaussian output. It is freely available at `Amber website`_ and the installation guide is in `Amber Manual`_. **Tsubasa** will directly call :code:`antechamber` to determine atom types and RESP charges.


Download Parmhess
-----------------
You can download the the newest release of **Parmhess** here_.

.. _here : https://github.com/ruixingw/parmhess/releases

You will get:

1. File "parmhess.py", "classdef.py". The **Parmhess** program.
2. File "tsubasa.py", "config.yml", "vdw.dat". A copy of the **Tsubasa** program, which is used to prepare the input files.
3. Folder "rxcclib".  A copy of the library code package **rxcclib**, which is used by **Tsubasa** and **Parmhess**.


Make it easier to use
---------------------

You may wish to create soft-links for **parmhess.py** and **tsubasa.py** to your :code:`$PATH`. Suppose :code:`~/bin` is in your :code:`$PATH`, then run:

.. code-block:: bash

    ln -s /PathToParmhess/parmhess.py ~/bin
    ln -s /PathToParmhess/tsubasa.py ~/bin

Alternatively, you may add Parmhess folder to your :code:`$PATH` by adding the following line to your :code:`~/.bashrc`:

.. code-block:: bash

   export PATH = "/PathToParmhess":$PATH

Either of them enables you to directly run :code:`parmhess.py` and :code:`tsubasa.py` at any directory.


.. _`Amber website` : http://ambermd.org/#AmberTools
.. _`Amber Manual` : http://ambermd.org/doc12/
.. _anaconda : https://www.continuum.io/downloads
.. _Python3: https://www.python.org/
.. _`AmberTools page` : http://ambermd.org/bugfixesat.html
.. _rxcclib: https://github.com/ruixingw/rxcclib
.. _Tsubasa: https://github.com/ruixingw/tsubasa
.. _`fixreadinesp.sh`_: http://ambermd.org/fixreadinesp.sh

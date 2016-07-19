============
Installation
============

Prerequisites
-------------

1. Python3 (Anaconda)

Parmhess is written in Python3_ and runs on Linux.

The easiest way to setup Python3 is to install Anaconda_, which is "the leading open data science platform powered by Python." It is freely available at:

.. _anaconda : https://www.continuum.io/downloads
.. _Python3: https://www.python.org/

    http://www.continuum.io/downloads

Be sure to install Anaconda for Linux (Python3.5), but not Python2.7 nor for other operating system.

If you are an experienced Python3 user and you prefer to use your own Python build, make sure that **numpy** and **pyyaml** be installed properly.

2. Gaussian 09 (G09)

Gaussian 09 is used to perform QM and MM calculation and should be properly installed on your machine. Any revision of G09 should work. However, due to `a bug of Revision B.01`_ (search "Gaussian 09 fix" in this page), some attention should be paid to RESP charge calculation in **Tsubasa**. To use G09 B01, See the chapter: `how to use G09 B01 in Tsubasa`__.

For MM calculation (**Parmhess**), :code:`g09` is directly called, while for QM calculation (**Tsubasa**), you could use your job submit command or script (see Tsubasa manual).

.. _`a bug of Revision B.01` : http://ambermd.org/bugfixesat.html
.. __ : ../latest/tsubasa.html

3. AmberTools

AmberTools is used to identify atom types and calculate RESP charge from Gaussian output. It is freely available at `Amber website`_ and the installation guide is in `Amber Manual`_. Make sure it is installed properly as **Tsubasa** will call :code:`antechamber` directly from :code:`$PATH`.

.. _`Amber website` : http://ambermd.org/#AmberTools
.. _`Amber Manual` : http://ambermd.org/doc12/



Download Parmhess
-----------------
Next, download Parmhess using any one of these commands:

.. code-block:: bash

    git clone https://github.com/ruixingw/parmhess  # Using git
    svn co https://github.com/ruixingw/parmhess  # or using svn
    wget https://github.com/ruixingw/parmhess/archive/master.zip -O parmhess.zip  # or download zip


You will get three things:

1. File  "parmhess.py". The main program.
2. Folder "tsubasa". A copy of the Tsubasa_ program, which is used to prepare the input files.
3. Folder "rxcclib".  A copy of the library code package rxcclib_, which is used by Tsubasa_ and Parmhess.

.. _rxcclib: https://github.com/ruixingw/rxcclib
.. _Tsubasa: https://github.com/ruixingw/tsubasa

Make it easier to use
---------------------

You may wish to create soft-links for **parmhess.py** and **tsubasa/tsubasa.py** to your $PATH, so that you could directly use them at any directory. Suppose :code:`~/bin` is in your :code:`$PATH`, then:

.. code-block:: bash

    ln -s ~/bin /PathToParmhess/parmhess.py
    ln -s ~/bin /PathToParmhess/tsubasa/tsubasa.py

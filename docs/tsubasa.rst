=========================
Prepare inputs by Tsubasa
=========================

**Parmhess** needs the following inputs.

1. Gaussian Formatted Checkpoint file (FCHK) and Gaussian Output (LOG) as QM input, which should be created from a QM frequency calculation. Hessian, forces, and internal coordinates are read from here.
2. An MM input file in Gaussian format. Coordinates, charges, atom types, and MM functions are read from here. Force constants waiting to be determined should be written as :code:`XXXXXX`.
3. An :code:`input.inp` file that includes the name of the aforementioned two files.

All these files can be prepared by **Tsubasa** program automatically.

If you prefer to prepare these files by hand, please refer to :doc:`whatsgoingon`.

.. _`Tsubasa` : http://github.com/ruixingw/tsubasa/


Tsubasa Manual
--------------


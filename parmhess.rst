===============================
Parameterization using Parmhess
===============================

Once the inputs are prepared, simply run:

.. code-block::  bash

   parmhess.py input.inp

Parmhess will read inputs and invoke **g09** command to do necessary MM calculations. These files will be generated:

1. Folder **hffiles**, includes all temporary files. They are used to calculate lowercase :math:`h` (see the paper).
2. Result MM input file including generated parameters. Named as PHF/FHF/IHF_result_mmXXX.com.



Acceptable arguments:



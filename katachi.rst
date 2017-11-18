Run Katachi Amendment
=====================

After parameterization, run:

.. code-block:: bash

  katachi.py ihf_result_mmH2O2.com 0 calcall 100

The result will be named katachi_ihf_result_mmH2O2.com.

Here are the last lines of the example:

::

  AmbTrs   ho  oh  oh  ho   0   0   0   0    0.000  1.552  0.000  0.000   1.0
  HrmBnd1 ho oh oh   68.779  98.2714
  HrmStr1 ho oh   550.931  0.96915
  HrmStr1 oh oh   346.974  1.44835
  Nonbon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 -1.2
  VDW    ho  0.0000  0.0000
  VDW    oh  1.7210  0.210


|

The detailed usage is as follows.

::

  > katachi.py --help
  usage: katachi.py [-h] mmresult loopid opt convthreshold
  positional arguments:
    mmresult       parameterized result MM file.
    loopid         loopid
    opt            opt or calcall
    convthreshold  convergence threshold


*loopid* is the number of starting cycle. It should be 0 if it is a fresh run. If Katachi was interrupted in the middle, this number can be set to that of the latest cycle to restore the process.

*opt* or *calcall* controls the commands used for MM optimization. If *opt* is specified, *opt=(nomicro,cartesian)* will be used for MM optimization, and after the convergence, *opt=(nomicro,cartesian,tight,calcall)* will then be used and start again from the previous result. If calcall is specified, all steps will use the latter keywords. In practice, we found that using calcall is usually better.

*convthreshold* is the cycle number threshold. If the result does not improve any more in the given cycles, the program will stop. In our tests, we used 100 for this threshold; that means, if the result does not improve in 100 cycles, the program will stop and use the best result before.


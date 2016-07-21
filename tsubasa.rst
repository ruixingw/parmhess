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

Acceptable arguments:

::

  [ruixingw@NTU test~]$ tsubasa.py -h
  usage: tsubasa.py [-h] [-i INPUTGEOM] [-c CONFIGFILE]
                  [--readvdw EXTERNALVDWFILE]
                  [--startfrom {freq,resp,antechamber,readmol2}]
                  [--stopafter {opt,freq,resp,antechamber}]

  optional arguments:
    -h, --help            show this help message and exit
    -i INPUTGEOM          Inputfile including molecular specs and connectivity
    -c CONFIGFILE         Tsubasa config file.
    --readvdw EXTERNALVDWFILE
                          If provided, read external vdW parameters from file
    --startfrom {freq,resp,antechamber,readmol2}
                          Start from a certain step.
                          Choices=['opt','freq','resp','antechamber','readmol2']
    --stopafter {opt,freq,resp,antechamber}
                          Stop after a certain step.
                          Choices=['freq','resp','antechamber']

Input Arguments
^^^^^^^^^^^^^^^

1. :code:`-i INPUTGEOM`: INPUTGEOM should be a file with :code:`.gau` extension. If not specified, program will search the current working directory(CWD) for a :code:`.gau` file , then use it if found, or exit if not.

The content of this file should include the geometry and connectivity in Gaussian format. An example is shown below:

::

  [ruixingw@NTU test~]$ cat H2O2.gau
  0 1
   O     -3.932366558742     -2.706592019388     -0.283087648364
   H     -4.104122594899     -1.885283961559      0.212183800621
   O     -2.475774075762     -2.710437185030     -0.283043419203
   H     -2.304047969868     -3.531709623942      0.212297456772

   1 2 1.0 3 1.0
   2
   3 4 1.0
   4
                                #(here is a blank line at last)
  [ruixingw@NTU test~]$ ls
   H2O2.gau
  [ruixingw@NTU test~]$ tsubasa.py   # H2O2.gau will be read as INPUTGEOM

2. :code:`-c CONFIGFILE`: CONFIGFILE should be a file with :code:`.yml` extension. If not specified, program will search the CWD for a :code:`.yml` file; if failed, a template :code:`.yml` file will be copied to CWD, named as same as :code:`.gau` file; if it is found, then check if the name of :code:`.yml` file and :code:`.gau` are the same and if they are not same, raise an error.

::

  [ruixingw@NTU test~]$ tsubasa.py   # H2O2.gau is readed as INPUTGEOM
   WARNING:root:Config file is not found. A template is copied to current directory. Program will now quit.
  [ruixingw@NTU test~]$ ls
  H2O2.gau    H2O2.yml
  [ruixingw@NTU test~]$ tsubasa.py # INPUTGEOM: H2O2.gau ; CONFIGFILE: H2O2.yml


3. :code:`--readvdw EXTERNALVDWFILE`: Tsubasa attached a :code:`vdw.dat` file including the vdW parameters in GAFF. If some certain vdW parameter is missing, you may provide a file to specify the vdW parameters. The format is same to :code:`vdw.dat`, which can be found in the source directory.

::

  [ruixingw@NTU test~]$ cat externalvdw.dat
   Zn          1.395   0.01491700         IOD Zn2+: Li,PF; Merz,Jr; JCTC2013 DOI:10.1021/ct400146w
  [ruixingw@NTU test~]$ tsubasa.py --readvdw externalvdw.dat # vdW for Zn is read


Job Control Arguments
^^^^^^^^^^^^^^^^^^^^^

**Tsubasa** have these steps:

- Optimization (opt)
- Frequency calculation (freq)
- MK charge calculation by Gaussian (resp)
- RESP charge calculation by antechamber (antechamber)
- Read all outputs and generate the MM input file (readmol2).


Normally, for small molecules, just running these steps in sequence should work. However, if the system is larger or complicated, these steps do not always succeed and may result in "Error termination". For such cases, you may check Gaussian options carefully and run them manually. After that, the success calculation can be provided to **Tsubasa** to do the following steps. **Note that the provided file must be named exactly same as the Tsubasa generated one.** Job flow can be controled by these arguments:

-  --startfrom {freq,resp,antechamber,readmol2}

   Choose one from the list. If not specified, the program starts from (opt) as normal.

-   --stopafter {opt,freq,resp,antechamber}

   Choose one from the list. If not specified, the program stops after (readmol2) as normal.



An example of the whole process is:

::

  [rwang013@boonlay-h00 test]$ ls
  H2O2.gau  H2O2.yml
  [rwang013@boonlay-h00 test]$ tsubasa.py
  INFO     Read config from H2O2.yml
  INFO     Runing optimization...                  # Start Optimization
  INFO     Run g09 : myg09boon optH2O2.com         # Submit optH2O2.com to PBS queue system
  INFO     Checking g09 termination for optH2O2.com...
  WARNING  No log file detected. Wait 2s..         # No log file detected due to the delay of queue system
  WARNING  Log file detected: optH2O2.log waiting for termination..
  INFO         ..normal termination
  INFO     Running frequency calculation...        # Start Frequency Calculation
  INFO     Run g09 : myg09boon freqH2O2.com
  INFO     Checking g09 termination for freqH2O2.com...
  182312.boonlay-h00
  WARNING  No log file detected. Wait 2s..
  WARNING  Log file detected: freqH2O2.log waiting for termination..
  INFO         ..normal termination
  INFO     Running RESP calculation...            # Start MK charge calculation
  INFO     Run g09 : myg09a2boon respH2O2.com
  INFO     Checking g09 termination for respH2O2.com...
  182313.boonlay-h00
  WARNING  No log file detected. Wait 2s..
  WARNING  Log file detected: respH2O2.log waiting for termination..
  INFO         ..normal termination
  INFO     Run antechamber:                      # Run antechamber for RESP
  INFO     Runing antechamber:                   # command is as below
  INFO     antechamber -c resp -i respH2O2.log -fi gout -o respH2O2.mol2 -fo mol2 -pf y

  INFO     Format CHK file by:                   # Prepare QM Fchk
  INFO       formchk -3 freqH2O2.chk freqH2O2.fchk
  INFO     Read fchk:freqH2O2.fchk

  [rwang013@boonlay-h00 test]$ ls  # inputs for Parmhess
  freqH2O2.fchk  freqH2O2.log  input.inp  mmH2O2.com  tsubasa/
  [rwang013@boonlay-h00 test]$ cd tsubasa/
  [rwang013@boonlay-h00 tsubasa]$ ls     # Temp files of Tsubasa
  freqH2O2.chk  freqH2O2.fchk  H2O2.gau      H2O2.yml   mmH2O2.com   optH2O2.com  respH2O2.chk  respH2O2.log  freqH2O2.com  freqH2O2.log   H2O2.tsubasa  input.inp  optH2O2.chk  optH2O2.log  respH2O2.com  respH2O2.mol2
  [rwang013@boonlay-h00 tsubasa]$


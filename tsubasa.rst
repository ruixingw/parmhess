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

Files of **Tsubasa**:

1. :code:`tsubasa.py` the main progam.
2. :code:`config.yml` the template config file
3. :code:`vdw.dat` includes vdW parameters from GAFF

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

1. :code:`-i INPUTGEOM`: INPUTGEOM should be a file with :code:`.gau` extension. If not specified, program will search the current working directory(CWD) for a :code:`.gau` file , then use it if found, or exit if failed.

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
                  #note here is a blank line at the end
  [ruixingw@NTU test~]$ ls
   H2O2.gau
  [ruixingw@NTU test~]$ tsubasa.py   # H2O2.gau will be read as INPUTGEOM

2. :code:`-c CONFIGFILE`: CONFIGFILE should be a file with :code:`.yml` extension. If not specified, program will search the CWD for a :code:`.yml` file; if failed, a template :code:`.yml` file will be copied to CWD, named as same as :code:`.gau` file; if it is found, then check if the name of :code:`.yml` file and :code:`.gau` are the same. If they are not same, program will raise an error and exit.

::

  [ruixingw@NTU test~]$ tsubasa.py   # H2O2.gau is readed as INPUTGEOM
   WARNING:root:Config file is not found.
   A template is copied to current directory. Program will now quit.
  [ruixingw@NTU test~]$ ls
  H2O2.gau    H2O2.yml
  [ruixingw@NTU test~]$ tsubasa.py   # INPUTGEOM: H2O2.gau ; CONFIGFILE: H2O2.yml


3. :code:`--readvdw EXTERNALVDWFILE`: Tsubasa attached a :code:`vdw.dat` file which includes the vdW parameters used in GAFF. If some certain vdW parameter is missing, you may provide a file to specify the vdW parameters. The format is same to :code:`vdw.dat`, which can be found in the source directory.

::

  [ruixingw@NTU test~]$ cat externalvdw.dat
   Zn          1.395   0.01491700         JCTC2013 DOI:10.1021/ct400146w
  [ruixingw@NTU test~]$ tsubasa.py --readvdw externalvdw.dat # vdW for Zn is read


Job Control Arguments
^^^^^^^^^^^^^^^^^^^^^

**Tsubasa** do these steps in sequence:

- Optimization (opt)
- Frequency calculation (freq)
- MK charge calculation by Gaussian (resp)
- RESP charge calculation by antechamber (antechamber)
- Read all outputs and generate the MM input file (readmol2).


Normally, just running these steps in sequence should work. However, if the system is larger or complicated, these steps do not always succeed and may result in "Error termination". For such cases, you may check Gaussian options and run calculation manually. The successful results can be provided to **Tsubasa** to do the following steps. **Note that the provided files must be named exactly same as those of Tsubasa generated.** Job flow can be controled by these arguments:

4. :code:`--startfrom {freq,resp,antechamber,readmol2}`  (choose one from the list)

   Read the existing files and restart from the specified step. If not specified, the program starts from the beginning (opt) as normal.

5. :code:`--stopafter {opt,freq,resp,antechamber}`  (choose one from the list)

   Stop after the specified step. If not specified, the program ends as normal.



An example of the whole process is:

::

  [ruixingw@NTU test]$ ls
  H2O2.gau  H2O2.yml
  [ruixingw@NTU test]$ tsubasa.py
  INFO     Read config from H2O2.yml
  INFO     Runing optimization...                  # Start Optimization
  INFO     Run g09 : myg09boon optH2O2.com         # Submit optH2O2.com to PBS
  INFO     Checking g09 termination for optH2O2.com...
  WARNING  No log file detected. Wait 2s..         # (lag of queue system...)
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

  [ruixingw@NTU test]$ ls  # inputs for Parmhess
  freqH2O2.fchk  freqH2O2.log  input.inp  mmH2O2.com  tsubasa/
  [ruixingw@NTU test]$ cd tsubasa/
  [ruixingw@NTU tsubasa]$ ls     # Temporary files of Tsubasa
  freqH2O2.chk  freqH2O2.fchk  H2O2.gau      H2O2.yml   mmH2O2.com   optH2O2.com
  respH2O2.chk  respH2O2.log  freqH2O2.com  freqH2O2.log   H2O2.tsubasa  input.inp
  optH2O2.chk  optH2O2.log  respH2O2.com  respH2O2.mol2
  [ruixingw@NTU tsubasa]$ cd ..
  [ruixingw@NTU test]$ cat mmH2O2.com       # MM input file is ready
  %mem=12gb
  #p amber=softonly geom=connectivity nosymm
  iop(4/33=3,7/33=1)
  freq=intmodes

  MM

  0 1
  O-oh--0.410452   -0.718633164030   -0.118472295063   -0.054617618503
  H-ho-0.410452   -1.023538245240    0.665457463989    0.436707052760
  O-oh--0.410452    0.718637032315    0.118468958072   -0.054573365001
  H-ho-0.410452    1.023507272500   -0.665430766999    0.436820814218

   1 2 1.0 3 1.0
   2
   3 4 1.0
   4

  AmbTrs ho oh oh ho 0 180 0 0 0.0 XXXXXX 0.0 0.0 1.0
  HrmBnd1 ho oh oh XXXXXX 100.2486
  HrmStr1 ho oh XXXXXX 0.97412
  HrmStr1 oh oh XXXXXX 1.45667
  Nonbon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 -1.2
  VDW ho  0.0000  0.0000
  VDW oh  1.7210  0.2104

  [rwang013@boonlay-h00 test]$
  
Config file
^^^^^^^^^^^

The config file includes the commands to run Gaussian and antechamber etc. The format is YAML_. An example is shown below.

.. _YAML: http://yaml.org/


::

  g09rt: myg09boon
  g09a2rt: myg09a2boon

  antechamber: antechamber -c resp
  clean: rm *gaussian*


  opthead: |
    %mem=16gb
    %nproc=12
    #p b3lyp/6-31+g* geom=connectivity
    int=ultrafine symm=(loose,follow)
    opt=(verytight,maxstep=7,notrust)

    opt-title


  opttail: |


  freqhead: |
    %mem=16gb
    %nproc=12
    #p b3lyp/chkbas int=ultrafine symm=loose geom=allcheck guess=tcheck freq=intmodes iop(7/33=1)


  resphead: |
    %mem=16gb
    %nproc=12
    #p b3lyp/chkbas
    iop(6/33=2,6/42=17,6/41=10)
    int=ultrafine symm=loose
    pop=mk
    geom=allcheck guess=tcheck


  resptail: |


  mmhead: |
    %mem=12gb
    #p amber=softonly geom=connectivity nosymm
    iop(4/33=3,7/33=1)
    freq=intmodes

    MM


Keywords:

1. :code:`g09rt`:  command to run Gaussian 09 for (opt, freq). Here, running "myg09boon test.com" should yield "test.log" in the same folder.

2. :code:`g09a2rt`:  command to run Gaussian 09 for resp (to avoid G09 B01 bug). Here, running "myg09a2boon test.com" should yield "test.log" in the same folder.

3. :code:`antechamber`:  command to run antechamber. Charge type may be modified. For large molecule(>100 atoms), "-pl 30" may be added (see :code:`antechamber -h` for details).

4. :code:`clean`: this command will be run at the end of all steps for clean purpose.

5. :code:`opthead`, :code:`opttail`, :code:`freqhead`, :code:`resphead`, :code:`resptail`, :code:`mmhead`: Keywords that are used to run Gaussian. The content of :code:`.gau` file will be pasted betwwen *head* and *tail* section. Blank lines will be inserted between them to form a Gaussian Input file. You may change them in your own need. Normally, keywords in :code:`freqhead` and :code:`mmhead` should not be changed. 


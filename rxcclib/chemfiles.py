# From fchk read Charge, Multiplicity, Coordinates
from __future__ import print_function
import os
import subprocess
import time
import logging
import shutil
import numpy as np
import rxcclib.cclibutils as cclibutils


class rxFileError(Exception):
    def __init(self, value):
        self.value = value

    def __repr__(self):
        return repr(self.value)

    __str__ = __repr__


class File(object):
    def __init__(self, name):
        self.pwd = os.path.abspath(os.path.join('.', name))  # filepath
        self.name = name   # filename with relative path
        self.pwd = os.path.split(self.name)[0]
        self._com = GauCOM(self)
        self._log = GauLOG(self)
        self._fchk = GauFCHK(self)
        self._mol2 = Mol2(self)
        self.default = 'fchk'
    # subfile objects

    @property
    def comname(self):
        return self.name + '.com'

    @property
    def logname(self):
        return self.name + '.log'

    @property
    def chkname(self):
        return self.name + '.chk'

    @property
    def fchkname(self):
        return self.name + '.fchk'

    @property
    def mol2name(self):
        return self.name + '.mol2'

    @property
    def com(self):
        return self._com

    @property
    def log(self):
        return self._log

    @property
    def fchk(self):
        return self._fchk

    @property
    def mol2(self):
        return self._mol2

    # molecule properties; setter should only be called by subfile methods
    @property
    def natoms(self):
        return self._natoms

    @natoms.setter
    def natoms(self, value):
        if not hasattr(self, '_natoms'):
            self._natoms = value
        else:
            if self._natoms != value:
                raise rxFileError(
                    'Error: natoms is already read,' +
                    " and not consistent with new value: Current value is "
                    + str(self._natoms) + ", New value is " + str(value))

    @property
    def multiplicity(self):
        return self._multiplicity

    @multiplicity.setter
    def multiplicity(self, value):
        if not hasattr(self, '_multiplicity'):
            self._multiplicity = value
        else:
            if self._multiplicity != value:
                raise rxFileError(
                    "Error: multiplicity is already read," +
                    " and not consistent with new value: Current value is "
                    + str(self._multiplicity) + ", New value is " + str(value))

    @property
    def totalcharge(self):
        return self._totalcharge

    @totalcharge.setter
    def totalcharge(self, value):
        if not hasattr(self, '_totalcharge'):
            self._totalcharge = value
        else:
            if self._totalcharge != value:
                raise rxFileError(
                    "Error: totalcharge is already read," +
                    " and not consistent with new value: Current value is "
                    + str(self._totalcharge) + ", New value is " + str(value))

    @property
    def atomtypelist(self):
        return self.mol2.atomtypelist

    @property
    def atomchargelist(self):
        return self.mol2.atomchargelist

    def runformchk(self):
        string = 'formchk -3 ' + self.chkname + ' ' + self.fchkname
        logging.info('  ' + string)
        iferror = os.popen(string)
        if iferror.read().find('Error') >= 0:
            raise rxFileError('   Error in formatting' + self.chkname)
            iferror.close()
            return False
        iferror.close()
        return True


class GauFCHK(object):
    def __init__(self, parent):
        self._parent = parent
        self.filename = self._parent.fchkname
        self.readstate = False

    def read(self):
        if self.readstate:
            logging.warning("fchk.read(): fchk is already read")
            return True

        self.coordslist = []
        self.atomlist = [None]
        self.xyz = ''
        self.readstate = True

        logging.info('Read fchk:' + self._parent.fchkname)
        # FCHK parser
        tmpfchkread = 'In ' + self.filename + '.read(), '
        with open(self._parent.fchkname, 'r') as f:
            for string in f:

                if string.find('Charge') == 0:
                    self.totalcharge = int(string.split('I')[1])
                    self._parent.totalcharge = self.totalcharge
                    logging.debug(tmpfchkread + 'Read Charge at line: '
                                  + string)

                if string.find('Multiplicity') == 0:
                    self.multiplicity = int(string.split('I')[1])
                    self._parent.multiplicity = self.multiplicity
                    logging.debug(tmpfchkread + 'Read Multiplicity at line: '
                                  + string)

                if string.find('Atomic numbers') == 0:
                    self.natoms = int(string.split('=')[1])
                    self._parent.natoms = self.natoms
                    logging.debug(tmpfchkread + 'Read Natoms at line: '
                                  + string)
                    for string in f:
                        try:
                            self.atomlist.extend(
                                [int(x) for x in string.split()])
                        except ValueError:
                            logging.debug(tmpfchkread + ' end read atomlist'
                                          ' at line: ' + string)
                            assert len(self.atomlist) == self.natoms+1, (
                                "Error: len(atomlist) != natoms ! ")
                            break

                if string.find('Current cartesian coordinates') == 0:
                    for string in f:
                        try:
                            self.coordslist.extend(
                                [float(x) for x in string.split()])
                        except ValueError:
                            tmp = (tmpfchkread + ' end read Cartesian'
                                   ' Coordinates at line: ' + string)
                            logging.debug(tmp)
                            assert (len(self.coordslist) == self.natoms * 3)
                            break

                # Read Hessian
                if string.find('Cartesian Force Constants') == 0:
                    self.hessian = []
                    for string in f:
                        try:
                            self.hessian.extend(
                                [float(x) for x in string.split()])
                        except ValueError:
                            tmp = ('In ' + self.filename +
                                   '. read(), end read'
                                   ' Hessian at line: ' + string)
                            logging.debug(tmp)
                            assert (len(self.hessian) ==
                                    4.5 * self.natoms ** 2 + 1.5 * self.natoms)
                            break
                # Read Internal Hessian
                if string.find('Internal Force Constants') == 0:
                    self.inthessian = []
                    for string in f:
                        try:
                            self.inthessian.extend(
                                [float(x) for x in string.split()])
                        except ValueError:
                            tmp = ('In ' + self.filename +
                                   '. read(), end read internal'
                                   ' Hessian at line: ' + string)
                            logging.debug(tmp)
                            break

        self.coordslist = [cclibutils.convertor(x, "bohr", "Angstrom")
                           for x in self.coordslist]
        self.coordslist = np.array(self.coordslist)
        self.atomlist = np.array(self.atomlist)
        self.hessian = np.array(self.hessian)
        for i in range(0, len(self.atomlist) - 1):
            tmp = str(self.atomlist[i + 1]) + '   ' + str(self.coordslist[
                3 * i]) + '   ' + str(self.coordslist[
                    3 * i + 1]) + '   ' + str(self.coordslist[3 * i +
                                                              2]) + '\n'
            self.xyz += tmp
        return True

    def findHessianElement(self, i, j):  # i, j: coordinate number
        if i < j:
            i, j = j, i
        num = i * (i - 1) / 2 + j
        num = int(num)
        return self.hessian[num - 1]

    def find33Hessian(self, i, j):  # i, j: atom number
        if i < j:
            i, j = j, i
        tthess = []
        i1 = 3 * (i - 1) + 1
        i2 = 3 * (i - 1) + 2
        i3 = 3 * i
        j1 = 3 * (j - 1) + 1
        j2 = 3 * (j - 1) + 2
        j3 = 3 * j
        tthess.append([self.findHessianElement(i1, j1),
                       self.findHessianElement(i1, j2),
                       self.findHessianElement(i1, j3)])
        tthess.append([self.findHessianElement(i2, j1),
                       self.findHessianElement(i2, j2),
                       self.findHessianElement(i2, j3)])
        tthess.append([self.findHessianElement(i3, j1),
                       self.findHessianElement(i3, j2),
                       self.findHessianElement(i3, j3)])
        tthess = np.array(tthess)
        return tthess

    def findintHessianElement(self, i, j):  # i, j: coordinate number
        if i < j:
            i, j = j, i
        num = i * (i - 1) / 2 + j
        num = int(num)
        return self.inthessian[num - 1]


class Mol2(object):
    def __init__(self, parent):
        self._parent = parent

    def read(self):
        self.atomtypelist = [None]
        self.atomchargelist = [None]
        with open(self._parent.mol2name, 'r') as f:
            for line in f:
                if line.find('ATOM') >= 0 and line.find('@') >= 0:
                    break
            for line in f:
                if line.find('@') >= 0:
                    break
                mol2 = line.split()
                self.atomtypelist.append(mol2[-4])
                self.atomchargelist.append(float(mol2[-1]))
        assert len(self.atomtypelist) == len(self.atomchargelist)
        return True


class GauCOM(object):
    g09rt = 'g09'
    g09a2rt = 'g09'

    def __init__(self, parent):
        self._parent = parent
        self.xyzfile = ''

    @property
    def parent(self):
        return self._parent

    # File operation
    def rung09(self):
        ifchk = 1  # if no chk, add.
        content = ''
        with open(self._parent.comname, 'r') as f:
            for line in f.readlines():
                if line.find('%chk') >= 0:
                    if line[line.find('=') + 1:] != self._parent.chkname:
                        oldname = line[line.find('=') + 1:].strip('\n')
                        newname = self._parent.chkname
                        if os.path.isfile(oldname):
                            try:
                                shutil.copyfile(oldname, newname)
                            except shutil.SameFileError:
                                pass
                    line = '%chk=' + newname + '\n'
                    ifchk = 0
                content += line
        if ifchk == 1:
            content = '%chk=' + self._parent.chkname + '\n' + content
        with open(self._parent.comname, 'w') as f:
            f.write(content)
        logging.info('Run g09 : ' + GauCOM.g09rt + ' ' +
                     self._parent.comname)
        self.running = subprocess.Popen([GauCOM.g09rt, self._parent.comname])

    def rung09a2(self):
        bak = GauCOM.g09rt
        GauCOM.g09rt = GauCOM.g09a2rt
        self.rung09()
        GauCOM.g09rt = bak

    def isover(self, wait=True):
        # Return None if not ended, True if normal term, False if error.
        logging.info('Checking g09 termination for ' + self._parent.comname +
                     '...')
        if wait is False:
            if self.running.poll() is None:
                return None
        elif wait is True:
            self.running.wait()
        else:
            raise
        while True:
            output = ''
            if not os.path.isfile(self._parent.logname):
                logging.warning('No log file detected. Wait 2s..')
                time.sleep(2)
                if os.path.isfile(self._parent.logname):
                    logging.warning('Log file detected: ' +
                                    self._parent.logname,
                                    'waiting for termination..')
                continue

            with open(self._parent.logname, 'r') as f:
                for x in f.readlines()[:-10:-1]:
                    output += x
            if output.find('Normal termination') >= 0:
                logging.debug('    ..normal termination')
                return True
            if output.find('Error termination') >= 0:
                logging.error('Error termination in ' +
                              self._parent.comname)
                raise rxFileError('G09 Error termination')
                return False
            if wait is False:
                return None
            time.sleep(0.5)


class GauLOG(object):
    antecommand = 'antechamber -c resp'

    def __init__(self, parent):
        self._parent = parent
        self.freq = 0

    def getnatoms(self):
        with open(self._parent.logname, 'r') as f:
            for x in f.readlines():
                if x.find('NAtoms') >= 0:
                    self.natoms = int(x.split()[1])
                    return self.natoms
                    break

    def coordslast(self):  # TODO
        return
        with open(self._parent.logname, 'r') as f:
            for x in list(reversed(f.readlines())):
                if x.find('orientation') >= 0:
                    self.orn = x
                    break
        return self.orn

    def getfreq(self):
        freq = []
        p = False
        with open(self._parent.logname, 'r') as f:
            for line in f.readlines():
                if line.find('Frequencies -- ') >= 0:
                    freq.extend(line.split()[2:])
                if line.find('#P') >= 0:
                    if p is False:
                        freq = []
                        p = True
            if freq != []:
                freq = list(map(float, freq))
                # if (len(freq)>3*self.natoms-5):
                #     freq=freq[:len(freq)/2]
            else:
                logging.error("No Freq found")
                raise rxFileError('No Freq found')
        self.freq = freq
        return self.freq

    def runantecham(self):
        logging.info('Runing antechamber: \n')
        command = (GauLOG.antecommand + ' -i ' +
                   self._parent.logname + ' -fi gout -o ' +
                   self._parent.mol2name + ' -fo mol2' +
                   ' -pf y')
        logging.info(command)
        os.system(command)

#!/usr/bin/env python3
from io import StringIO
import numpy as np
import rxcclib.File.chemfiles as rxfile


class DihdFunction(object):
    def __init__(self, periodicity, phase, forceconst, npaths, dihd):
        self.periodicity = periodicity
        self.forceconst = forceconst
        self.phase = phase
        self.npaths = npaths
        self.mydihd = dihd
        self.repr = dihd.repr

    def __str__(self):
        return str(self.forceconst)

    __repr__ = __str__


class MMFunction(object):
    unknownsign = 'XXXXXX'

    def __init__(self, line):
        fun = line.split()
        self.type = None
        self.repr = None

        def newfloat(value):
            if value == MMFunction.unknownsign:
                return MMFunction.unknownsign
            else:
                return float(value)

        if fun[0] == 'AmbTrs':
            self.type = 'dihd'
            self.a = fun[1]
            self.b = fun[2]
            self.c = fun[3]
            self.d = fun[4]
            self.forceconst = None
            self.dihdfunctions = []
            phase = []
            self.npaths = float(fun[13])
            for x in fun[5:9]:
                phase.append(int(x))
            for i, paras in enumerate(fun[9:13]):
                self.dihdfunctions.append(DihdFunction(i + 1,
                                                       phase[i],
                                                       newfloat(paras),
                                                       self.npaths,
                                                       self))
            self.repr = self.a + ' ' + self.b + ' ' + self.c + ' ' + self.d
        elif fun[0] == 'HrmBnd1':
            self.type = 'angle'
            self.a = fun[1]
            self.b = fun[2]
            self.c = fun[3]
            self.forceconst = newfloat(fun[4])
            self.eqvalue = newfloat(fun[5])
            self.repr = self.a + ' ' + self.b + ' ' + self.c
        elif fun[0] == 'HrmStr1':
            self.type = 'bond'
            self.a = fun[1]
            self.b = fun[2]
            self.forceconst = newfloat(fun[3])
            self.eqvalue = newfloat(fun[4])
            self.repr = self.a + ' ' + self.b
        elif fun[0] == 'ImpTrs':
            self.type = 'improper'
            self.a = fun[1]
            self.b = fun[2]
            self.c = fun[3]
            self.d = fun[4]
            self.forceconst = newfloat(fun[5])
            self.phase = newfloat(fun[6])
            self.periodicity = newfloat(fun[7])
            self.repr = self.a + ' ' + self.b + ' ' + self.c + ' ' + self.d
        elif fun[0] == 'VDW':
            self.type = 'vdw'
            self.content = line
            self.atomtype = fun[1]
            self.radius = fun[2]
            self.welldepth = fun[3]
        else:
            self.type = 'else'
            self.content = line


class GauAmberCOM(rxfile.GauCOM):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent._com = self

    def read(self):
        self.atomlist = [None]
        self.atomtypelist = [None]
        self.atomchargelist = [None]
        self.coordslist = []
        self.connectivity = ''
        self.dihdfunc = []
        self.anglefunc = []
        self.bondfunc = []
        self.improperfunc = []
        self.additionfunc = []
        self.vdw = []
        self.xyz = ''
        self.vdwdict = {}
        with open(self.parent.comname, 'r') as f:
            content = f.read()
            tmp = content.split('\n')
            block = ''
            for item in tmp:
                if item.isspace():
                    item = ''
                block += item + '\n'
            block = block.split('\n\n')
            block = [x + '\n' for x in block]

        # route, title, molespecs, connectivity, mmfunctions
        blockindex = [0, 1, 2, 3, 4]
        if block[0].find('allcheck') >= 0:
            blockindex[1] = -1
            blockindex[2] = -1
            blockindex[1:] = [x - 2 for x in blockindex[1:]]
        if block[0].find('connectivity') < 0:
            blockindex[3] = -1
            blockindex[3:] = [x - 1 for x in blockindex[3:]]

        def molespecs(line):
            tmp = line.split()[0]
            if tmp.find('-') >= 0:
                self.atomlist.append(tmp.split('-')[0])
                if tmp.count('-') == 2:
                    tmp = tmp.split('-')
                    self.atomtypelist.append(tmp[1])
                    self.atomchargelist.append(tmp[2])
                elif tmp.count('-') == 3:
                    tmp = tmp.split('-')
                    self.atomtypelist.append(tmp[1])
                    self.atomchargelist.append(-float(tmp[3]))
            else:
                self.atomlist.append(tmp)
            self.coordslist.extend(line.split()[1:4])

        for index, item in enumerate(block):
            if index == blockindex[0]:
                self.route = item
            if index == blockindex[1]:
                self.title = item
            if index == blockindex[2]:
                f = StringIO(item)
                line = next(f)
                self.totalcharge = line.split()[0]
                self.multiplicity = line.split()[1]
                for line in f:
                    molespecs(line)
            if index == blockindex[3]:
                f = StringIO(item)
                for line in f:
                    self.connectivity += line
            if index == blockindex[4]:
                f = StringIO(item)
                for line in f:
                    thisline = MMFunction(line)
                    if thisline.type == 'dihd':
                        self.dihdfunc.append(thisline)
                    elif thisline.type == 'angle':
                        self.anglefunc.append(thisline)
                    elif thisline.type == 'bond':
                        self.bondfunc.append(thisline)
                    elif thisline.type == 'else':
                        self.additionfunc.append(thisline)
                    elif thisline.type == 'vdw':
                        self.vdw.append(thisline)
                        self.vdwdict.update({thisline.atomtype: (
                            thisline.radius, thisline.welldepth)})
                    elif thisline.type == 'improper':
                        self.improperfunc.append(thisline)

        self.coordslist = np.array(self.coordslist)
        for i in range(0, len(self.atomlist) - 1):
            tmp = str(self.atomlist[i + 1]) + '   ' + str(self.coordslist[
                3 * i]) + '   ' + str(self.coordslist[
                    3 * i + 1]) + '   ' + str(self.coordslist[3 * i +
                                                              2]) + '\n'
            self.xyz += tmp
        return True

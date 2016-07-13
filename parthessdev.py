#!/usr/bin/env python3
import argparse
import subprocess
import logging
import itertools
import copy
import os
import shutil
import pdb
from io import StringIO

import numpy as np

import rxcclib.molecules as rxmol
import rxcclib.chemfiles as rxfile


class DihdForceConst(object):
    def __init__(self, value, dihd):
        self.forceconst = value
        self.dihd = dihd
        self.repr = dihd.repr

    # def __repr__(self):
    #     return ("Dihedral forceconst: " + str(self.forceconst) +
    #             " of dihedral " + self.repr)

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
            self.forceconst = []
            self.phase = []
            self.npaths = float(fun[13])
            for paras in fun[9:13]:
                self.forceconst.append(DihdForceConst(newfloat(paras), self))
            for phase in fun[5:9]:
                self.phase.append(int(phase))
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
            self.npaths = newfloat(fun[7])
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
        self.xyzfile = ''
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
                self.xyzfile += line
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


def matchdihd(dihd, func):
    a = (dihd[1].atomtype == func.a or func.a == '*')
    b = (dihd[2].atomtype == func.b or func.b == '*')
    c = (dihd[3].atomtype == func.c or func.c == '*')
    d = (dihd[4].atomtype == func.d or func.d == '*')
    forward = a and b and c and d
    a = (dihd[1].atomtype == func.d or func.d == '*')
    b = (dihd[2].atomtype == func.c or func.c == '*')
    c = (dihd[3].atomtype == func.b or func.b == '*')
    d = (dihd[4].atomtype == func.a or func.a == '*')
    backward = a and b and c and d
    if forward or backward:
        return True
    else:
        return False


def matchbond(bond, func):
    a = (bond[1].atomtype == func.a or func.a == '*')
    b = (bond[2].atomtype == func.b or func.b == '*')
    forward = a and b
    a = (bond[1].atomtype == func.b or func.b == '*')
    b = (bond[2].atomtype == func.a or func.a == '*')
    backward = a and b
    if forward or backward:
        return True
    else:
        return False


def matchangle(angle, func):
    a = (angle[1].atomtype == func.a or func.a == '*')
    b = (angle[2].atomtype == func.b or func.b == '*')
    c = (angle[3].atomtype == func.c or func.c == '*')
    forward = a and b and c
    a = (angle[1].atomtype == func.c or func.c == '*')
    b = (angle[2].atomtype == func.b or func.b == '*')
    c = (angle[3].atomtype == func.a or func.a == '*')
    backward = a and b and c
    if forward or backward:
        return True
    else:
        return False


def matchimproper(improper, func):
    a = (improper[1].atomtype == func.a or func.a == '*')
    b = (improper[2].atomtype == func.b or func.b == '*')
    c = (improper[3].atomtype == func.c or func.c == '*')
    d = (improper[4].atomtype == func.d or func.d == '*')
    forward = a and b and c and d
    if forward:
        return True
    else:
        return False


# Unit Hessian Component Tail
def hesstail(obj, itnlcordL, hessvdwtail, i=0):
    global mmcom
    tailstring = ''
    for item in itnlcordL:
        if type(item) == rxmol.Dihd:
            parm = ['0.000', '0.000', '0.000', '0.000']
            if obj is item:
                parm[i] = '1.000'
            tailstring += 'AmbTrs  ' + ' '.join([
                x.center(3, ' ') for x in item.repr.split()
            ]) + '  ' + ' '.join(
                [str(x).center(3, ' ') for x in item.phase]) + '  ' + ' '.join(
                    [str(x) for x in parm]) + '   ' + str(item.npaths) + '\n'
        else:
            parm = '0.000'
            if obj is item:
                parm = '1.000'
            if type(item) == rxmol.Angle:
                tailstring += 'HrmBnd1  ' + ' '.join([
                    x.center(3, ' ') for x in item.repr.split()
                ]) + '  ' + parm + '  ' + '{:>9.5f}'.format(
                    item.anglevalue) + '\n'
            elif type(item) == rxmol.Bond:
                tailstring += 'HrmStr1  ' + ' '.join([
                    x.center(3, ' ') for x in item.repr.split()
                ]) + '  ' + parm + '  ' + '{:>7.5f}'.format(item.length) + '\n'
            elif type(item) == rxmol.Improper:
                tailstring += 'ImpTrs  ' + ' '.join([
                    x.center(3, ' ') for x in item.repr.split()
                ]) + '  ' + parm + '  ' + '{:6.2f}'.format(
                    item.phase) + '  ' + str(item.npaths) + '\n'
    for x in mmcom.additionfunc:
        tailstring += x.content
    tailstring += hessvdwtail
    tailstring += '\n\n'
    return tailstring


# Hprime(known) Hessian Component Tail
def hprimetail(itnlcordL, hprimevdwtail):
    global mole
    tailstring = ''
    for item in itnlcordL:
        if type(item) == rxmol.Dihd:
            parm = []
            dihd = item
            for item in dihd.forceconst:
                if str(item) == MMFunction.unknownsign:
                    parm.append('0.000')
                else:
                    parm.append(item)
            tailstring += (
                'AmbTrs  ' + ' '.join(
                    [x.center(3, ' ') for x in dihd.repr.split()]) + '  ' +
                ' '.join([str(x).center(3, ' ') for x in dihd.phase]) + '  ' +
                ' '.join([str(x) for x in parm]) + '   ' + str(dihd.npaths) +
                '\n')
        else:
            if item.forceconst == MMFunction.unknownsign:
                parm = '0.000'
            else:
                parm = str(item.forceconst)

            if type(item) == rxmol.Angle:

                tailstring += ('HrmBnd1  ' + ' '.join([
                    x.center(3, ' ') for x in item.repr.split()
                ]) + '  ' + parm + '  ' + '{:>9.5f}'.format(item.anglevalue) +
                               '\n')
            if type(item) == rxmol.Bond:
                tailstring += ('HrmStr1  ' + ' '.join([
                    x.center(3, ' ') for x in item.repr.split()
                ]) + '  ' + parm + '  ' + '{:>7.5f}'.format(item.length) + '\n'
                               )
            if type(item) == rxmol.Improper:
                tailstring += ('ImpTrs  ' + ' '.join(
                    [x.center(3, ' ') for x in item.repr.split()]) + '  ' +
                               parm + '  ' + '{:6.2f}'.format(item.phase) +
                               '  ' + str(item.npaths) + '\n')
    for x in mmcom.additionfunc:
        tailstring += x.content
    tailstring += hprimevdwtail
    tailstring += '\n\n'
    return tailstring


def calcphfgroup(adict, thishprime, qmfchk):
    for key, value in adict.items():
        leftL, rightL = [], []
        a, b = [int(x) for x in key.split('-')]
        hq = qmfchk.fchk.find33Hessian(a, b)
        hp = thishprime.fchk.find33Hessian(a, b)
        hideal = hq - hp
        hideal = [x for item in hideal for x in item]
        hk = []
        for item in value:
            tmp = item.hessfile.fchk.find33Hessian(a, b)
            hk.append([x for row in tmp for x in row])
        rightL = hideal
        leftL = list(zip(*hk))
        leftL = np.array(leftL)
        rightL = np.array(rightL)
        res = np.linalg.lstsq(leftL, rightL)[0]
        for i, item in enumerate(value):
            item.forceconst = res[i]


def summarize(unkL, itnlcordL, originalname, finalhead, method):
    # Summarize
    logging.info('Start Summarizing')
    for func in unkL:
        if func.forceconst == MMFunction.unknownsign:
            res = 0
            i = 0
            for item in itnlcordL:
                if item.func == func:
                    res += float(str(item.forceconst))
                    i += 1
            func.forceconst = res / i
            print(func.type, func.repr, func.forceconst)
        if func.type == 'dihd':
            res = [0.0, 0.0, 0.0, 0.0]
            i = 0
            for item in itnlcordL:
                if item.func == func:
                    res = [float(str(x)) + oldx
                           for x, oldx in zip(item.forceconst, res)]
                    i += 1
            for i, item in enumerate(func.forceconst):
                item.forceconst = res[i]
            print(func.type, func.repr, func.forceconst)

    # Build tailstring
    tailstring = ''
    for dihd in mmcom.dihdfunc:
        parm = []
        for item in dihd.forceconst:
            if str(item) == MMFunction.unknownsign:
                parm.append('0.000')
                logging.critical('Force constant is not'
                                 ' determined for dihedral ' +
                                 dihd.repr)
                raise
            else:
                parm.append(float(str(item)))
        tailstring += 'AmbTrs  ' + ' '.join(
            [x.center(3, ' ') for x in dihd.repr.split()]) + '  ' + ' '.join(
                [str(x).center(3, ' ') for x in dihd.phase]) + '  ' + ' '.join(
                    ['{:>6.3f}'.format(x)
                     for x in parm]) + '   ' + str(dihd.npaths) + '\n'
    for angle in mmcom.anglefunc:
        if angle.forceconst == MMFunction.unknownsign:
            parm = '0.000'
            logging.critical('Force constant is not determined for angle ' +
                             angle.repr)
            raise
        else:
            parm = angle.forceconst
        tailstring += 'HrmBnd1  ' + ' '.join([x.center(
            3, ' ') for x in angle.repr.split()]) + '  ' + '{:>7.3f}'.format(
                parm) + '  ' + '{:>9.5f}'.format(angle.eqvalue) + '\n'
    for bond in mmcom.bondfunc:
        if bond.forceconst == MMFunction.unknownsign:
            parm = '0.000'
            logging.critical('Force constant is not determined for bond ' +
                             bond.repr)
            raise
        else:
            parm = bond.forceconst
        tailstring += 'HrmStr1  ' + ' '.join([x.center(
            3, ' ') for x in bond.repr.split()]) + '  ' + '{:>8.3f}'.format(
                parm) + '  ' + '{:>7.5f}'.format(bond.eqvalue) + '\n'
    for improper in mmcom.improperfunc:
        if improper.forceconst == MMFunction.unknownsign:
            logging.critical('Force constant is not determined for improper ' +
                             improper.repr)
            raise
        else:
            parm = improper.forceconst
        tailstring += 'ImpTrs  ' + ' '.join([
            x.center(3, ' ') for x in improper.repr.split()
        ]) + '  ' + '{:>7.3f}'.format(parm) + '  ' + '{:6.2f}'.format(
            improper.phase) + '  ' + str(improper.npaths) + '\n'
    for x in mmcom.additionfunc:
        tailstring += x.content
    for vdw in mmcom.vdw:
        tailstring += 'VDW  ' + '  ' + vdw.atomtype + \
            '  ' + vdw.radius + '  ' + vdw.welldepth + '\n'
    tailstring += '\n\n'

    logging.info('\n\nresult:\n')
    logging.info(tailstring)

    finalname = method + '_result_' + originalname
    with open(finalname, 'w') as f:
        logging.info('Write result to file' + finalname)
        f.write(finalhead + tailstring)
    shutil.copy(finalname, os.path.join('..', finalname))


def main(args):
    global mmcom
    inputinp = args.inputinp
    quiet = args.quiet
    nocalc = args.nocalc
    with open(inputinp, 'r') as f:
        for line in f:
            if line.find('mmfile') >= 0:
                mmfile = line.split('=')[1].strip('\n')
            elif line.find('qmfchk') >= 0:
                qmfchk = line.split('=')[1].strip('\n')

    # Create work directory
    if nocalc:
        os.chdir('hffiles')
    else:
        if os.path.lexists('hffiles'):
            shutil.rmtree('hffiles')
        os.mkdir('hffiles')
        os.chdir('hffiles')

    originalname = os.path.basename(mmfile)

    # Logging module setting. Print INFO on screen and DEBUG INFO in file======
    logging.basicConfig(filename=originalname + '.phflog',
                        level=logging.DEBUG,
                        filemode='w')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(levelname)-8s %(message)s')
    console.setFormatter(formatter)
    if not quiet:
        logging.getLogger('').addHandler(console)
    #  ==================================================================

    logging.info('Start Hessian Fitting for MM parameterization\n\n')
    logging.info('Release under GNU LGPL License')
    logging.info('If you use this program and/or the Hessian fitting method' +
                 ', cite DOI:XXX \n\n')

    logging.info("Provided mmfile input: " + mmfile + ' and qmfchk input: ' +
                 qmfchk)

    shutil.copy(os.path.join('..', mmfile), '.')
    shutil.copy(os.path.join('..', qmfchk), '.')

    mmfile = os.path.splitext(mmfile)[0]
    qmfchk = os.path.splitext(qmfchk)[0]

    #   Instantialize Molecule and file
    mole = rxmol.Molecule('mmcom')
    mmcom = GauAmberCOM(rxfile.File(mmfile))
    qmfchk = rxfile.File(qmfchk)
    qmfchk.fchk.read()
    mmcom.read()
    mole.readfromxyz(mmcom.xyz)
    mole.readchargefromlist(mmcom.atomchargelist)
    mole.readtypefromlist(mmcom.atomtypelist)
    mole.readconnectivity(mmcom.connectivity)

    #   Count unknown parameters
    nunk = 0
    unkL = []
    unkL.extend(sorted(mmcom.bondfunc, key=lambda x: x.repr))
    unkL.extend(sorted(mmcom.anglefunc, key=lambda x: x.repr))
    unkL.extend(sorted(mmcom.dihdfunc, key=lambda x: x.repr))
    unkL.extend(sorted(mmcom.improperfunc, key=lambda x: x.repr))

    # Count unknowns: special treat to dihedrals;
    # Match improper here by the way
    for item in unkL:
        if item.type == 'dihd':
            for paras in item.forceconst:
                if paras == MMFunction.unknownsign:
                    nunk += 1
                    paras.known = False
                else:
                    paras.known = True
        elif item.type == 'improper':
            # /* match improper
            for atom3 in mole:
                if atom3.atomtype == item.c:
                    permu = list(itertools.permutations(atom3.neighbor, 3))
                    for tu in permu:
                        a = tu[0].atomtype == item.a or item.a == '*'
                        b = tu[1].atomtype == item.b or item.b == '*'
                        c = tu[2].atomtype == item.d or item.d == '*'
                        if a and b and c:
                            mole.addimproper(tu[0].atomnum, tu[1].atomnum,
                                             atom3.atomnum, tu[2].atomnum)
                            break
            # */

            if item.forceconst == MMFunction.unknownsign:
                nunk += 1
                item.known = False
            else:
                item.known = True
        else:
            if item.forceconst == MMFunction.unknownsign:
                nunk += 1
                item.known = False
            else:
                item.known = True

    # Prepare Head section
    hessxyz = ''
    hprimexyz = ''
    finalxyz = ''
    for atom in mole:
        hprimexyz += (
            atom.elementsym + '-' + atom.name + '-' +
            '{:8.6f}'.format(float(atom.atomcharge)) + '   ' + '   '.join(
                ["{: .12f}".format(x) for x in atom.coords]) + '\n')
        hessxyz += (
            atom.elementsym + '-' + atom.name + '-0.000000' + '   ' +
            '    '.join(["{: .12f}".format(x) for x in atom.coords]) + '\n')
        finalxyz += (
            atom.elementsym + '-' + atom.atomtype + '-' +
            '{:8.6f}'.format(float(atom.atomcharge)) + '   ' + '    '.join(
                ["{: .12f}".format(x) for x in atom.coords]) + '\n')

    hesshead = (mmcom.route + '\nhess\n\n' + str(qmfchk.totalcharge) + ' ' +
                str(qmfchk.multiplicity) + '\n' + hessxyz + '\n' +
                mmcom.connectivity + '\n')
    hprimehead = (mmcom.route + '\nhprime\n\n' + str(qmfchk.totalcharge) + ' '
                  + str(qmfchk.multiplicity) + '\n' + hprimexyz + '\n' +
                  mmcom.connectivity + '\n')
    finalhead = (mmcom.route + '\nfinal\n\n' + str(qmfchk.totalcharge) + ' ' +
                 str(qmfchk.multiplicity) + '\n' + finalxyz + '\n' +
                 mmcom.connectivity + '\n')

    # vdW parameters
    vdwdict = mmcom.vdwdict
    hessvdwtail = ''
    hprimevdwtail = ''

    for atom in mole:
        atom.vdwradius = vdwdict[atom.atomtype][0]
        atom.vdwwelldepth = vdwdict[atom.atomtype][1]
        hessvdwtail += (
            'VDW ' + atom.name + '  ' + atom.vdwradius + '  0.0000\n')
        hprimevdwtail += 'VDW ' + atom.name + '  ' + \
            atom.vdwradius + ' ' + atom.vdwwelldepth + '\n'

    # Match internal coordinate and MMFunction
    for dihd in mole.dihdlist.values():
        for dihdfunc in mmcom.dihdfunc:
            if matchdihd(dihd, dihdfunc):
                dihd.func = dihdfunc
                dihd.forceconst = [DihdForceConst(x, dihd)
                                   for x in dihdfunc.forceconst]
                dihd.phase = copy.copy(dihdfunc.phase)
                dihd.npaths = dihdfunc.npaths
                for dihdfc in dihd.forceconst:
                    if dihdfc == MMFunction.unknownsign:
                        dihdfc.known = False
                    else:
                        dihdfc.known = True
    for angle in mole.anglelist.values():
        for anglefunc in mmcom.anglefunc:
            if matchangle(angle, anglefunc):
                angle.func = anglefunc
                angle.forceconst = anglefunc.forceconst
                if angle.forceconst == MMFunction.unknownsign:
                    angle.known = False
                else:
                    angle.known = True
    for bond in mole.bondlist.values():
        for bondfunc in mmcom.bondfunc:
            if matchbond(bond, bondfunc):
                bond.func = bondfunc
                bond.forceconst = bondfunc.forceconst
                if bond.forceconst == MMFunction.unknownsign:
                    bond.known = False
                else:
                    bond.known = True

    for improper in mole.improperlist.values():
        for improperfunc in mmcom.improperfunc:
            if matchimproper(improper, improperfunc):
                improper.func = improperfunc
                improper.forceconst = improperfunc.forceconst
                if improper.forceconst == MMFunction.unknownsign:
                    improper.known = False
                else:
                    improper.known = True

    # Count real force constants of all unknown internal coordinates
    itnlcordL = []
    itnlcordL.extend(sorted(mole.dihdlist.values(), key=lambda x: x.repr))
    itnlcordL.extend(sorted(mole.anglelist.values(), key=lambda x: x.repr))
    itnlcordL.extend(sorted(mole.bondlist.values(), key=lambda x: x.repr))
    itnlcordL.extend(sorted(mole.improperlist.values(), key=lambda x: x.repr))
    realnunk = 0

    for item in itnlcordL:
        if type(item) == rxmol.Dihd:
            for index, parms in enumerate(item.forceconst):
                if str(parms) == MMFunction.unknownsign:
                    realnunk += 1
        elif type(item) == rxmol.Angle:
            if item.forceconst == MMFunction.unknownsign:
                realnunk += 1
        elif type(item) == rxmol.Improper:
            if item.forceconst == MMFunction.unknownsign:
                realnunk += 1
        elif type(item) == rxmol.Bond:
            if item.forceconst == MMFunction.unknownsign:
                realnunk += 1

    # Prepare Unit Hessian File
    hess = []
    num = realnunk
    unkparmL = []
    for obj in itnlcordL:
        if type(obj) == rxmol.Dihd:
            for i, parms in enumerate(obj.forceconst):
                if str(parms) == MMFunction.unknownsign:
                    with open('hess' + str(len(hess)) + '.com', 'w') as f:
                        f.write(hesshead + hesstail(obj, itnlcordL,
                                                    hessvdwtail, i))
                    this = rxfile.File('hess' + str(len(hess)))
                    obj.forceconst[i].hessfile = this
                    this.orig = obj.forceconst[i]
                    unkparmL.append(this.orig)
                    hess.append(this)
                    if not nocalc:
                        this.com.rung09()
                        try:
                            this.com.isover()
                        except rxfile.rxFileError:
                            this.com.rung09()
                            this.com.isover()
                    this.runformchk()
                    this.fchk.read()
                    num -= 1
                    logging.info(str(num + 3) + ' left')
        else:
            if obj.forceconst == MMFunction.unknownsign:
                with open('hess' + str(len(hess)) + '.com', 'w') as f:
                    f.write(hesshead + hesstail(obj, itnlcordL, hessvdwtail))
                this = rxfile.File('hess' + str(len(hess)))
                obj.hessfile = this
                this.orig = obj
                unkparmL.append(this.orig)
                hess.append(this)
                if not nocalc:
                    this.com.rung09()
                    try:
                        this.com.isover()
                    except:
                        this.com.rung09()
                        this.com.isover()
                this.runformchk()
                num -= 1
                this.fchk.read()
                logging.info(str(num + 3) + ' left')
    # Identify Coupled Terms

    links = {}
    for obj in unkparmL:
        if type(obj) == DihdForceConst:
            a = obj.dihd[1].atomnum
            b = obj.dihd[4].atomnum
            if a > b:
                a, b = b, a
            if str(a) + '-' + str(b) in links.keys():
                links[str(a) + '-' + str(b)].append(obj)
            else:
                links.update({str(a) + '-' + str(b): [obj]})
        if type(obj) == rxmol.Angle:
            a = obj[1].atomnum
            b = obj[3].atomnum
            if a > b:
                a, b = b, a
            if str(a) + '-' + str(b) in links.keys():
                links[str(a) + '-' + str(b)].append(obj)
            else:
                links.update({str(a) + '-' + str(b): [obj]})
        if type(obj) == rxmol.Bond:
            a = obj[1].atomnum
            b = obj[2].atomnum
            if a > b:
                a, b = b, a
            if str(a) + '-' + str(b) in links.keys():
                links[str(a) + '-' + str(b)].append(obj)
            else:
                links.update({str(a) + '-' + str(b): [obj]})
        if type(obj) == rxmol.Improper:
            a = obj[1].atomnum
            b = obj[4].atomnum
            if a > b:
                a, b = b, a
            if str(a) + '-' + str(b) in links.keys():
                links[str(a) + '-' + str(b)].append(obj)
            else:
                links.update({str(a) + '-' + str(b): [obj]})

    # store links
    onetwoL = {}
    onetricL = {}
    onetriucL = {}
    onefourL = {}
    onetrifourL = {}
    for key, value in links.items():
        types = []
        for item in value:
            types.append(type(item))
        dbool = DihdForceConst in types
        abool = rxmol.Angle in types
        bbool = rxmol.Bond in types
        ibool = rxmol.Improper in types
        if dbool and abool:
            onetrifourL.update({key: value})
        elif dbool:
            onefourL.update({key: value})
        elif abool and ibool:
            onetricL.update({key: value})
        elif abool:
            onetriucL.update({key: value})
        elif bbool:
            onetwoL.update({key: value})

    # Start parameterization.
    # sequence: onefour-->onetrifour-->onetri(coupled)
    # -->onetri(uncoupled)-->onetwo
    if onefourL:
        onefourhprime = hprimehead + hprimetail(itnlcordL, hprimevdwtail)
        with open('onefourhprime.com', 'w') as f:
            f.write(onefourhprime)
        onefourhprime = rxfile.File('onefourhprime')
        if not nocalc:
            onefourhprime.com.rung09()
            onefourhprime.com.isover()
        onefourhprime.runformchk()
        onefourhprime.fchk.read()
        calcphfgroup(onefourL, onefourhprime, qmfchk)
    if onetrifourL:
        onetrifourhprime = hprimehead + hprimetail(itnlcordL, hprimevdwtail)
        with open('onetrifourhprime.com', 'w') as f:
            f.write(onetrifourhprime)
        onetrifourhprime = rxfile.File('onetrifourhprime')
        if not nocalc:
            onetrifourhprime.com.rung09()
            onetrifourhprime.com.isover()
        onetrifourhprime.runformchk()
        onetrifourhprime.fchk.read()
        calcphfgroup(onetrifourL, onetrifourhprime, qmfchk)

    if onetricL:
        onetrichprime = hprimehead + hprimetail(itnlcordL, hprimevdwtail)
        with open('onetrichprime.com', 'w') as f:
            f.write(onetrichprime)
        onetrichprime = rxfile.File('onetrichprime')
        if not nocalc:
            onetrichprime.com.rung09()
            onetrichprime.com.isover()
        onetrichprime.runformchk()
        onetrichprime.fchk.read()
        calcphfgroup(onetricL, onetrichprime, qmfchk)

    if onetriucL:
        onetriucLhprime = hprimehead + hprimetail(itnlcordL, hprimevdwtail)
        with open('onetriucLhprime.com', 'w') as f:
            f.write(onetriucLhprime)
        onetriucLhprime = rxfile.File('onetriucLhprime')
        if not nocalc:
            onetriucLhprime.com.rung09()
            onetriucLhprime.com.isover()
        onetriucLhprime.runformchk()
        onetriucLhprime.fchk.read()
        calcphfgroup(onetriucL, onetriucLhprime, qmfchk)

    if onetwoL:
        onetwoLhprime = hprimehead + hprimetail(itnlcordL, hprimevdwtail)
        with open('onetwoLhprime.com', 'w') as f:
            f.write(onetwoLhprime)
        onetwoLhprime = rxfile.File('onetwoLhprime')
        if not nocalc:
            onetwoLhprime.com.rung09()
            onetwoLhprime.com.isover()
        onetwoLhprime.runformchk()
        onetwoLhprime.fchk.read()
        calcphfgroup(onetwoL, onetwoLhprime, qmfchk)

    # Summarize PHF and write to file
    summarize(unkL, itnlcordL, originalname, finalhead, 'phf')

    # End of PHF
    # Clean up:

    for item in itnlcordL:
        if type(item) != rxmol.Dihd:
            if item.known is False:
                item.forceconst = MMFunction.unknownsign
        else:
            for parms in item.forceconst:
                if parms.known is False:
                    parms = MMFunction.unknownsign
    for item in unkL:
        if item.type != 'dihd':
            if item.known is False:
                item.forceconst = MMFunction.unknownsign
        else:
            for parms in item.forceconst:
                if parms.known is False:
                    parms.forceconst = MMFunction.unknownsign

    # Start of FHF
    leftL = []
    hideal = []
    hprime = onefourhprime
    # for i in range(0, len(qmfchk.fchk.hessian)):
    #     leftL.append([])
    #     for item in unkparmL:
    #         leftL[-1].append(item.hessfile.fchk.hessian[i])
    hideal = qmfchk.fchk.hessian - hprime.fchk.hessian
    for item in unkparmL:
        tmp = item.hessfile.fchk.hessian
        leftL.append([x for x in tmp])

    leftL = list(zip(*leftL))
    print(len(leftL))
    print(len(leftL[1]))
    print(len(hideal))

    results = np.linalg.lstsq(leftL, hideal)[0]
    for i, item in enumerate(unkparmL):
        item.forceconst = results[i]
    summarize(unkL, itnlcordL, originalname, finalhead, 'fhf')


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('inputinp',
                        default='input.inp',
                        help='input.inp prepared by tsubasa.')
    parser.add_argument('--quiet',
                        '-q',
                        action='store_true',
                        help=('Do not show info on screen'))
    parser.add_argument('--nocalc', '-nc', action='store_true',
                        help=('Use already calculated file.'))

    args = parser.parse_args()
    main(args)

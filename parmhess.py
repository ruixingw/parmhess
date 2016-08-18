#!/usr/bin/env python3
import argparse
import logging
import itertools
import copy
import os
import sys
import shutil
import numpy as np
from io import StringIO
import rxcclib.molecules as rxmol
import rxcclib.chemfiles as rxfile
from classdef import InternalCoordinates
from classdef import DihdForceConst
from classdef import MMFunction
from classdef import GauAmberCOM


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
        print(leftL, rightL)
        for i, item in enumerate(value):
            item.forceconst = res[i]
            print(item.forceconst)
            print(item.repr)


def exactsum(geom, itnlcordL, originalname, hprimehead):
    exacttail = ''
    for item in itnlcordL:
        if type(item) == rxmol.Dihd:
            parm = [float(str(x)) for x in item.forceconst]
            print(parm)
            exacttail += ('AmbTrs  ' + ' '.join(
                [x.center(3, ' ') for x in item.repr.split()]) + '  '
                          + ' '.join(
                              [str(x).center(3, ' ') for x in item.phase]) +
                          '  ' + ' '.join(
                              ['{:>.10f}'.format(x) for x in parm])
                          + '   ' + str(item.npaths) + '\n')
        elif type(item) == rxmol.Angle:
            exacttail += ('HrmBnd1  ' + ' '.join([x.center(
                3, ' ') for x in item.repr.split()])
                          + '  ' + '{:>.10f}'.format(
                              item.forceconst) +
                          '  ' + '{:>9.5f}'.format(item.anglevalue) + '\n')
        elif type(item) == rxmol.Bond:
            exacttail += ('HrmStr1  ' + ' '.join([x.center(
                3, ' ') for x in item.repr.split()])
                          + '  ' + '{:>.10f}'.format(
                              item.forceconst) + '  ' +
                          '{:>7.5f}'.format(item.length) + '\n')
        elif type(item) == rxmol.Improper:
            exacttail += ('ImpTrs  ' + ' '.join([
                x.center(3, ' ') for x in item.repr.split()]) +
                          '  ' + '{:>.10f}'.format(item.forceconst) +
                          '  ' + '{:6.2f}'.format(
                              item.phase) + '  ' +
                          str(item.npaths) + '\n')
    for x in mmcom.additionfunc:
        exacttail += x.content
    for atom in geom:
        exacttail += 'VDW  ' + '  ' + atom.name + \
            '  ' + atom.vdwradius + '  ' + atom.vdwwelldepth + '\n'
    exacttail += '\n\n'

    finalname = 'exactihf_result_' + originalname
    logging.info('Write result to file ' + finalname)
    logging.info(exacttail)
    with open(finalname, 'w') as f:
        f.write(hprimehead + exacttail)
    shutil.copy(finalname, os.path.join('..', finalname))
    return rxfile.File(finalname.strip('.com'))


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
        if func.type == 'dihd':
            res = [0.0, 0.0, 0.0, 0.0]
            i = 0
            for item in itnlcordL:
                if item.func == func:
                    res = [float(str(x)) + oldx
                           for x, oldx in zip(item.forceconst, res)]
                    i += 1
            for index, item in enumerate(func.forceconst):
                item.forceconst = res[index] / i


    # Build tailstring
    tailstring = ''
    for dihd in mmcom.dihdfunc:
        parm = []
        for i, item in enumerate(dihd.forceconst):
            if str(item) == MMFunction.unknownsign:
                parm.append('0.000')
                logging.critical('Force constant is not'
                                 ' determined for dihedral ' + dihd.repr)
                raise
            else:
                tmp = float(str(item))
                if tmp < 0:
                    dihd.phase[i] = 180
                parm.append(abs(tmp))

        tailstring += 'AmbTrs  ' + ' '.join(
            [x.center(3, ' ') for x in dihd.repr.split()]) + '  ' + ' '.join(
                [str(x).center(3, ' ') for x in dihd.phase]) + '  ' + ' '.join(
                    ['{:>.10f}'.format(x)
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
            3, ' ') for x in angle.repr.split()]) + '  ' + '{:>.10f}'.format(
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
            3, ' ') for x in bond.repr.split()]) + '  ' + '{:>.10f}'.format(
                parm) + '  ' + '{:>7.5f}'.format(bond.eqvalue) + '\n'
    for improper in mmcom.improperfunc:
        if improper.forceconst == MMFunction.unknownsign:
            logging.critical('Force constant is not determined for improper ' +
                             improper.repr)
            raise
        else:
            parm = improper.forceconst
            if parm < 0:
                improper.phase = 180
                parm = -parm
        tailstring += 'ImpTrs  ' + ' '.join([
            x.center(3, ' ') for x in improper.repr.split()
        ]) + '  ' + '{:>.10f}'.format(parm) + '  ' + '{:6.2f}'.format(
            improper.phase) + '  ' + str(improper.npaths) + '\n'
    for x in mmcom.additionfunc:
        tailstring += x.content
    for vdw in mmcom.vdw:
        tailstring += 'VDW  ' + '  ' + vdw.atomtype + \
            '  ' + vdw.radius + '  ' + vdw.welldepth + '\n'
    tailstring += '\n\n'

    finalname = method + '_result_' + originalname
    logging.info('Write result to file ' + finalname)
    logging.info(tailstring)
    with open(finalname, 'w') as f:

        f.write(finalhead + tailstring)
    shutil.copy(finalname, os.path.join('..', finalname))


def addlink1(mmfile, itnlcordL):
    fakecontent = '--link1--\n'
    fakecontent += ('%chk=' + mmfile.chkname + '\n')
    fakecontent += ('#p geom=allcheck ')
    fakecontent += ('freq=(readfc,modredundant,intmodes) '
                    'iop(4/33=3,7/33=1,99/5=5)\n\n')
    fakecontent += ('* * K\n* * * K\n* * * * K\n')
    for item in itnlcordL:
        if type(item) == rxmol.Bond:
            fakecontent += (str(item[1].atomnum) + ' ' +
                            str(item[2].atomnum) + ' A\n')
        if type(item) == rxmol.Angle:
            fakecontent += (str(item[1].atomnum) + ' ' +
                            str(item[2].atomnum) + ' ' +
                            str(item[3].atomnum) + ' A\n')
        if type(item) == rxmol.Dihd:
            fakecontent += (str(item[1].atomnum) + ' ' +
                            str(item[2].atomnum) + ' ' +
                            str(item[3].atomnum) + ' ' +
                            str(item[4].atomnum) + ' A\n')
        if type(item) == rxmol.Improper:
            fakecontent += (str(item[1].atomnum) + ' ' +
                            str(item[2].atomnum) + ' ' +
                            str(item[3].atomnum) + ' ' +
                            str(item[4].atomnum) + ' A\n')
    fakecontent += '\n'
    return fakecontent


def main(args):
    global mmcom
    np.set_printoptions(edgeitems=200)
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
                 ', cite DOI:10.1002/jcc.24457 \n\n')

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

    shutil.copy(os.path.join('..', qmfchk.logname), '.')
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
                    success = []
                    for tu in permu:
                        a = tu[0].atomtype == item.a or item.a == '*'
                        b = tu[1].atomtype == item.b or item.b == '*'
                        c = tu[2].atomtype == item.d or item.d == '*'
                        if a and b and c:
                            success.append([tu[0].atomnum, tu[1].atomnum,
                                            atom3.atomnum, tu[2].atomnum])
                    success = sorted(success, key=lambda x: (str(x[1]) +
                                                             str(x[3])))
                    success = success[0]
                    mole.addimproper(*success)
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
                improper.phase = improperfunc.phase
                improper.npaths = improperfunc.npaths
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
            i = 0
            for index, parms in enumerate(item.forceconst):
                if str(parms) == MMFunction.unknownsign:
                    realnunk += 1
                    i += 1
            if i > 1:
                item.multipleterms = True
            else:
                item.multipleterms = False
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
                    this = rxfile.File('hess' + str(len(hess)))
                    with open('hess' + str(len(hess)) + '.com', 'w') as f:
                        f.write('%chk=' + this.chkname + '\n' +
                                hesshead +
                                hesstail(obj, itnlcordL, hessvdwtail, i) +
                                addlink1(this, itnlcordL))
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
                this = rxfile.File('hess' + str(len(hess)))
                with open('hess' + str(len(hess)) + '.com', 'w') as f:
                    f.write('%chk=' + this.chkname + '\n' +
                            hesshead +
                            hesstail(obj, itnlcordL, hessvdwtail)+
                            addlink1(this, itnlcordL))
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
        onefourhprime = rxfile.File('onefourhprime')
        tmp = ('%chk=' + onefourhprime.chkname + '\n' +
               hprimehead + hprimetail(itnlcordL, hprimevdwtail))
        tmp += addlink1(onefourhprime, itnlcordL)
        with open(onefourhprime.comname, 'w') as f:
            f.write(tmp)
        if not nocalc:
            onefourhprime.com.rung09()
            onefourhprime.com.isover()
            onefourhprime.runformchk()
        onefourhprime.fchk.read()
        calcphfgroup(onefourL, onefourhprime, qmfchk)
    if onetrifourL:
        onetrifourhprime = rxfile.File('onetrifourhprime')
        tmp =('%chk=' + onetrifourhprime.chkname + '\n' +
              hprimehead + hprimetail(itnlcordL, hprimevdwtail))
        tmp += addlink1(onetrifourhprime, itnlcordL)
        with open(onetrifourhprime.comname, 'w') as f:
            f.write(tmp)
        if not nocalc:
            onetrifourhprime.com.rung09()
            onetrifourhprime.com.isover()
            onetrifourhprime.runformchk()
        onetrifourhprime.fchk.read()
        calcphfgroup(onetrifourL, onetrifourhprime, qmfchk)

    if onetricL:
        onetrichprime = rxfile.File('onetrichprime')
        tmp = ('%chk=' + onetrichprime.chkname + '\n' +
               hprimehead + hprimetail(itnlcordL, hprimevdwtail))
        tmp += addlink1(onetrichprime, itnlcordL)
        with open(onetrichprime.comname, 'w') as f:
            f.write(tmp)
        if not nocalc:
            onetrichprime.com.rung09()
            onetrichprime.com.isover()
            onetrichprime.runformchk()
        onetrichprime.fchk.read()
        calcphfgroup(onetricL, onetrichprime, qmfchk)

    if onetriucL:
        onetriuchprime = rxfile.File('onetriuchprime')
        tmp = ('%chk=' + onetriuchprime.chkname + '\n' +
               hprimehead + hprimetail(itnlcordL, hprimevdwtail))
        tmp += addlink1(onetriuchprime, itnlcordL)
        with open(onetriuchprime.comname, 'w') as f:
            f.write(tmp)

        if not nocalc:
            onetriuchprime.com.rung09()
            onetriuchprime.com.isover()
            onetriuchprime.runformchk()
        onetriuchprime.fchk.read()
        calcphfgroup(onetriucL, onetriuchprime, qmfchk)


    if onetwoL:
        onetwohprime = rxfile.File('onetwohprime')
        tmp = ('%chk=' + onetwohprime.chkname + '\n' +
               hprimehead + hprimetail(itnlcordL, hprimevdwtail))
        tmp += addlink1(onetwohprime, itnlcordL)
        with open(onetwohprime.comname, 'w') as f:
            f.write(tmp)

        if not nocalc:
            onetwohprime.com.rung09()
            onetwohprime.com.isover()
            onetwohprime.runformchk()
        onetwohprime.fchk.read()
        calcphfgroup(onetwoL, onetwohprime, qmfchk)

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
    try:
        hprime = onefourhprime
    except UnboundLocalError:
        try:
            hprime = onetrifourhprime
        except UnboundLocalError:
            try:
                hprime = onetrichprime
            except UnboundLocalError:
                try:
                    hprime = onetriuchprime
                except UnboundLocalError:
                    hprime = onetwohprime

    hideal = qmfchk.fchk.hessian - hprime.fchk.hessian
    for item in unkparmL:
        tmp = item.hessfile.fchk.hessian
        leftL.append([x for x in tmp])

    leftL = list(zip(*leftL))

    results = np.linalg.lstsq(leftL, hideal)[0]
    for i, item in enumerate(unkparmL):
        item.forceconst = results[i]
    summarize(unkL, itnlcordL, originalname, finalhead, 'fhf')

    # End of FHF
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
    # Start of IHF

    def readintcoords(fileobj):
        intcords = []
        with open(fileobj.logname) as f:
            tmp = f.read()
        tmp = tmp.split('Initial command')
        # assert len(tmp) == 3, fileobj.logname + str(len(tmp))
        tmp = tmp[-1]
        with StringIO(tmp) as f:
            for line in f:
                if line.find('Initial Parameters') < 0:
                    continue
                break
            [next(f) for x in range(0, 4)]
            for line in f:
                if line.find('---') >= 0:
                    break
                line = line.split()
                line = line[2]
                line = line[2:-1]
                line = [int(x) for x in line.split(',')]
                intcords.append(InternalCoordinates(fileobj,
                                                    line))
        return intcords
    # Check internal coordinates order matching
    qmfchk.itnl = readintcoords(qmfchk)
    assert len(qmfchk.itnl) == len(itnlcordL)
    for item in hess:
        item.itnl = readintcoords(item)
        assert len(item.itnl) == len(itnlcordL), item.comname
        for x, y in zip(item.itnl, qmfchk.itnl):
            assert x.atomset == y.atomset
    hprime.itnl = readintcoords(hprime)
    assert len(hprime.itnl) == len(itnlcordL)
    for x, y in zip(hprime.itnl, qmfchk.itnl):
        assert x.atomset == y.atomset

    # Read internal Hessian
    for i in range(1, len(qmfchk.itnl)+1):
        qmfchk.itnl[i-1].hessian = qmfchk.fchk.findintHessianElement(i, i)
        qmfchk.itnl[i-1].force = qmfchk.fchk.intforces[i-1]
        hprime.itnl[i-1].hessian = hprime.fchk.findintHessianElement(i, i)
        hprime.itnl[i-1].force = hprime.fchk.intforces[i-1]
        for item in hess:
            item.itnl[i-1].hessian = item.fchk.findintHessianElement(i, i)
            item.itnl[i-1].force = item.fchk.intforces[i-1]

    # Match itnlcordL and qmfchk.itnl:
    for item in itnlcordL:
        atomset = []
        if type(item) == rxmol.Improper:
            a = item[1].atomnum
            b = item[2].atomnum
            c = item[3].atomnum
            d = item[4].atomnum
            if b > c:
                atomset = [d, c, b, a]
            else:
                atomset = [a, b, c, d]
        else:
            for atom in item:
                atomset.append(atom.atomnum)
        matchres = filter(lambda x: x.atomset == atomset, qmfchk.itnl)
        matchres = list(matchres)
        assert len(matchres) == 1, ('Multiple matching of itnlcordL '
                                    'with qmfchk.itnl in IHF.' + item.repr)
        matchres = matchres[0]
        matchres.myselfitnl = item
        item.gauitnl = matchres

    # Build leftL in itnl order:
    # Hessian Left
    leftL = []
    for item in unkparmL:
        leftL.append([])
        tmp = item.hessfile.itnl
        leftL[-1].extend(x.hessian for x in tmp)

    leftL = list(zip(*leftL))
    # Force Left
    dihdL = list(sorted(mole.dihdlist.values(), key=lambda x: x.repr))
    forcesL = []
    forceleft = []
    for item in dihdL:
        if item.multipleterms is True:
            forcesL.append(item)
    for dihd in forcesL:
        tmp = []
        index = qmfchk.itnl.index(dihd.gauitnl)
        for item in hess:
            if type(item.orig) != DihdForceConst:
                tmp.append(0.0)
                print('out:', item.orig.repr)
                continue
            print('in:', item.orig.repr)
            value = (item.fchk.intforces[index])
            tmp.append(value)
        leftL.append(tmp)
        forceleft.append(tmp)

    # Build b  vector in itnl order:
    # Hessian Right
    rightL = []
    for qm, hp in zip(qmfchk.itnl, hprime.itnl):
        rightL.append(qm.hessian - hp.hessian)
    # Force Right
    forceright = []
    for dihd in forcesL:
        index = qmfchk.itnl.index(dihd.gauitnl)
        tmp = qmfchk.fchk.intforces[index] - hprime.fchk.intforces[index]
        rightL.append(tmp)
        forceright.append(tmp)

    print('forceleft')
    print(np.array(forceleft))
    print('forceright')
    print(np.array(forceright))
    # print(qmfchk.fchk.intforces[index])
    # print(hprime.fchk.intforces[index])
    # print(rightL[-1])
    # print(dihd.gauitnl.atomset)

    leftL = np.array(leftL)
    rightL = np.array(rightL)
#    print(leftL[index])
    res = np.linalg.lstsq(leftL, rightL)[0]

    for i, item in enumerate(unkparmL):
        item.forceconst = res[i]
        print(item.repr, res[i])
    summarize(unkL, itnlcordL, originalname, finalhead, 'ihf')

    exactsum(mole, itnlcordL, originalname, hprimehead)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('inputinp',
                        default='input.inp',
                        help='input.inp prepared by tsubasa.')
    parser.add_argument('--quiet',
                        '-q',
                        action='store_true',
                        help=('Do not show info on screen'))
    parser.add_argument('--nocalc',
                        '-nc',
                        action='store_true',
                        help=('Use already calculated file.'))

    args = parser.parse_args()
    main(args)

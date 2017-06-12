#!/usr/bin/env python3
import argparse
import logging
import itertools
import copy
import os
import shutil
import numpy as np
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
            tailstring += 'AmbTrs  ' + ' '.join(
                [x.center(3, ' ')
                 for x in item.repr.split()]) + '  ' + ' '.join(
                     [str(x).center(3, ' ')
                      for x in item.phase]) + '  ' + ' '.join(
                          [str(x)
                           for x in parm]) + '   ' + str(item.npaths) + '\n'
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
            tailstring += ('AmbTrs  ' + ' '.join(
                [x.center(3, ' ')
                 for x in dihd.repr.split()]) + '  ' + ' '.join(
                     [str(x).center(3, ' ')
                      for x in dihd.phase]) + '  ' + ' '.join(
                          [str(x)
                           for x in parm]) + '   ' + str(dihd.npaths) + '\n')
        else:
            if item.forceconst == MMFunction.unknownsign:
                parm = '0.000'
            else:
                parm = str(item.forceconst)

            if type(item) == rxmol.Angle:

                tailstring += ('HrmBnd1  ' + ' '.join(
                    [x.center(3, ' ')
                     for x in item.repr.split()]) + '  ' + parm + '  ' +
                               '{:>9.5f}'.format(item.anglevalue) + '\n')
            if type(item) == rxmol.Bond:
                tailstring += ('HrmStr1  ' + ' '.join(
                    [x.center(3, ' ')
                     for x in item.repr.split()]) + '  ' + parm + '  ' +
                               '{:>7.5f}'.format(item.length) + '\n')
            if type(item) == rxmol.Improper:
                tailstring += ('ImpTrs  ' + ' '.join([
                    x.center(3, ' ') for x in item.repr.split()
                ]) + '  ' + parm + '  ' + '{:6.2f}'.format(item.phase) + '  ' +
                               str(item.npaths) + '\n')
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
        if func.type == 'dihd':
            res = [0.0, 0.0, 0.0, 0.0]
            i = 0
            for item in itnlcordL:
                if item.func == func:
                    res = [
                        float(str(x)) + oldx
                        for x, oldx in zip(item.forceconst, res)
                    ]
                    i += 1
            for index, item in enumerate(func.forceconst):
                item.forceconst = res[index] / i

    # Build tailstring
    tailstring = ''
    for dihd in mmcom.dihdfunc:
        parm = []
        for i,item in enumerate(dihd.forceconst):
            if str(item) == MMFunction.unknownsign:
                parm.append('0.000')
                logging.critical('Force constant is not'
                                 ' determined for dihedral ' + dihd.repr)
                raise
            else:
                value = float(str(item))
                if value < 0:
                    dihd.phase[i] = 180
                    parm.append(abs(value))
                else:
                    parm.append(value)
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
        tailstring += 'HrmBnd1  ' + ' '.join(
            [x.center(3, ' ')
             for x in angle.repr.split()]) + '  ' + '{:>7.3f}'.format(
                 parm) + '  ' + '{:>9.5f}'.format(angle.eqvalue) + '\n'
    for bond in mmcom.bondfunc:
        if bond.forceconst == MMFunction.unknownsign:
            parm = '0.000'
            logging.critical('Force constant is not determined for bond ' +
                             bond.repr)
            raise
        else:
            parm = bond.forceconst
        tailstring += 'HrmStr1  ' + ' '.join(
            [x.center(3, ' ')
             for x in bond.repr.split()]) + '  ' + '{:>8.3f}'.format(
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

    finalname = method + '_result_' + originalname
    logging.info('Write result to file ' + finalname)
    logging.info(tailstring)
    with open(finalname, 'w') as f:

        f.write(finalhead + tailstring)
    shutil.copy(finalname, os.path.join('..', finalname))
    return finalname


def addlink1(mmfile, itnlcordL):
    #    return ''
    content = '\n--link1--\n'
    content += ('%chk=' + mmfile.chkname + '\n')
    content += ('#p geom=allcheck ')
    content += ('freq=(readfc,modredundant,intmodes) '
                'iop(4/33=3,7/33=1,99/5=5)\n\n')
    content += ('* * K\n* * * K\n* * * * K\n')
    for item in itnlcordL:
        if type(item) == rxmol.Bond:
            content += (
                str(item[1].atomnum) + ' ' + str(item[2].atomnum) + ' A\n')
        if type(item) == rxmol.Angle:
            content += (str(item[1].atomnum) + ' ' + str(item[2].atomnum) + ' '
                        + str(item[3].atomnum) + ' A\n')
        if type(item) == rxmol.Dihd:
            content += (
                str(item[1].atomnum) + ' ' + str(item[2].atomnum) + ' ' +
                str(item[3].atomnum) + ' ' + str(item[4].atomnum) + ' A\n')
        if type(item) == rxmol.Improper:
            content += (
                str(item[1].atomnum) + ' ' + str(item[2].atomnum) + ' ' +
                str(item[3].atomnum) + ' ' + str(item[4].atomnum) + ' A\n')
    content += '\n'
    return content


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'inputinp', default='input.inp', help='input.inp prepared by tsubasa.')
    parser.add_argument(
        '--quiet',
        '-q',
        action='store_true',
        help=('Do not show info on screen'))
    parser.add_argument(
        '--nocalc',
        '-nc',
        action='store_true',
        help=('Use already calculated file.'))

    args = parser.parse_args()
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
    logging.basicConfig(
        filename=originalname + '.phflog', level=logging.DEBUG, filemode='w')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(levelname)-8s %(message)s')
    console.setFormatter(formatter)
    if not quiet:
        logging.getLogger('').addHandler(console)
    #  ==================================================================

    logging.info('Start Hessian Fitting for MM parameterization\n')
    logging.info('This program is released under BSD 3-Clause License')
    logging.info('See details at https://github.com/ruixingw/parmhess\n\n')

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
    bakroute = mmcom.route
    mmcom.route = '%mem=1gb\n#p opt=(nomicro,redundant,modredundant,maxcyc=1)\namber=softonly geom=connectivity nosymm\niop(4/33=3)\n'

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
                    res = []
                    for tu in permu:
                        a = tu[0].atomtype == item.a or item.a == '*'
                        b = tu[1].atomtype == item.b or item.b == '*'
                        c = tu[2].atomtype == item.d or item.d == '*'
                        if a and b and c:
                            res.append([
                                tu[0].atomnum, tu[1].atomnum, atom3.atomnum,
                                tu[2].atomnum
                            ])
                    res = sorted(res, key=lambda x: (str(x[1]) + str(x[3])))
                    res = res[0]
                    mole.addimproper(*res)
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
        hprimexyz += (atom.elementsym + '-' + atom.name + '-' +
                      '{:8.6f}'.format(float(atom.atomcharge)) + '   ' +
                      '   '.join(["{: .12f}".format(x)
                                  for x in atom.coords]) + '\n')
        hessxyz += (atom.elementsym + '-' + atom.name + '-0.000000' + '   ' +
                    '    '.join(["{: .12f}".format(x)
                                 for x in atom.coords]) + '\n')
        finalxyz += (atom.elementsym + '-' + atom.atomtype + '-' +
                     '{:8.6f}'.format(float(atom.atomcharge)) + '   ' +
                     '    '.join(["{: .12f}".format(x)
                                  for x in atom.coords]) + '\n')

    hesshead = (mmcom.route + '\nhess\n\n' + str(qmfchk.totalcharge) + ' ' +
                str(qmfchk.multiplicity) + '\n' + hessxyz + '\n' +
                mmcom.connectivity + '\n')
    hprimehead = (mmcom.route + '\nhprime\n\n' + str(qmfchk.totalcharge) + ' '
                  + str(qmfchk.multiplicity) + '\n' + hprimexyz + '\n' +
                  mmcom.connectivity + '\n')
    finalhead = (bakroute + '\nfinal\n\n' + str(qmfchk.totalcharge) + ' ' +
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
                dihd.forceconst = [
                    DihdForceConst(x, dihd) for x in dihdfunc.forceconst
                ]
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

    # MM freq file
    
    content = ('* * K\n* * * K\n* * * * K\n')
    for item in itnlcordL:
        if type(item) == rxmol.Bond:
            content += (
                str(item[1].atomnum) + ' ' + str(item[2].atomnum) + ' A\n')
        if type(item) == rxmol.Angle:
            content += (str(item[1].atomnum) + ' ' + str(item[2].atomnum) + ' '
                        + str(item[3].atomnum) + ' A\n')
        if type(item) == rxmol.Dihd:
            content += (
                str(item[1].atomnum) + ' ' + str(item[2].atomnum) + ' ' +
                str(item[3].atomnum) + ' ' + str(item[4].atomnum) + ' A\n')
        if type(item) == rxmol.Improper:
            content += (
                str(item[1].atomnum) + ' ' + str(item[2].atomnum) + ' ' +
                str(item[3].atomnum) + ' ' + str(item[4].atomnum) + ' A\n')
    content += '\n'
    modred = content
    content = ('#p geom=allcheck nosymm ')
    content += ('amber=softonly freq=intmodes iop(4/33=3,7/33=1,99/5=5)\n\n')
    mmfreqhead = content

    def newrun(this,tail):
        this.com.rung09()
        try:
            this.com.isover()
        except rxfile.rxFileError:
            pass
        with open("freq_"+this.comname,'w') as f:
            f.write("%chk=" + this.chkname+'\n' + mmfreqhead + '\n' + tail)
        logging.info('  ..freq')
        os.system('g09 freq_'+this.comname)
        logging.info('  ..ok')
        this.runformchk()
        pass

    def newhesstail(*args):
        return modred + hesstail(*args)

    def newhprimetail(*args):
        return modred + hprimetail(*args)
    # Prepare Unit Hessian File
    hess = []
    num = realnunk
    unkparmL = []
    for obj in itnlcordL:
        if type(obj) == rxmol.Dihd:
            for i, parms in enumerate(obj.forceconst):
                if str(parms) == MMFunction.unknownsign:
                    with open('hess' + str(len(hess)) + '.com', 'w') as f:
                        f.write(hesshead + newhesstail(obj, itnlcordL,
                                                    hessvdwtail, i))
                    this = rxfile.File('hess' + str(len(hess)))
                    obj.forceconst[i].hessfile = this
                    this.orig = obj.forceconst[i]
                    unkparmL.append(this.orig)
                    hess.append(this)
                    if not nocalc:
                        newrun(this, hesstail(obj, itnlcordL,
                                              hessvdwtail, i))
                    this.fchk.read()
                    num -= 1
                    logging.info(str(num + 3) + ' left')
        else:
            if obj.forceconst == MMFunction.unknownsign:
                with open('hess' + str(len(hess)) + '.com', 'w') as f:
                    f.write(hesshead + newhesstail(obj, itnlcordL, hessvdwtail))
                this = rxfile.File('hess' + str(len(hess)))
                obj.hessfile = this
                this.orig = obj
                unkparmL.append(this.orig)
                hess.append(this)
                if not nocalc:
                    newrun(this, hesstail(obj, itnlcordL, hessvdwtail))
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
        onefourhprime = hprimehead + newhprimetail(itnlcordL, hprimevdwtail)
        with open('onefourhprime.com', 'w') as f:
            f.write(onefourhprime)
        onefourhprime = rxfile.File('onefourhprime')
        if not nocalc:
            newrun(onefourhprime,hprimetail(itnlcordL, hprimevdwtail))
            ##############################
            # onefourhprime.com.rung09() #
            # onefourhprime.com.isover() #
            # onefourhprime.runformchk() #
            ##############################
        onefourhprime.fchk.read()
        calcphfgroup(onefourL, onefourhprime, qmfchk)
    if onetrifourL:
        onetrifourhprime = hprimehead + newhprimetail(itnlcordL, hprimevdwtail)
        with open('onetrifourhprime.com', 'w') as f:
            f.write(onetrifourhprime)
        onetrifourhprime = rxfile.File('onetrifourhprime')
        if not nocalc:
            newrun(onetrifourhprime,hprimetail(itnlcordL, hprimevdwtail))
            #onetrifourhprime.com.rung09()
            #onetrifourhprime.com.isover()
            #onetrifourhprime.runformchk()
        onetrifourhprime.fchk.read()
        calcphfgroup(onetrifourL, onetrifourhprime, qmfchk)

    if onetricL:
        onetrichprime = hprimehead + newhprimetail(itnlcordL, hprimevdwtail)
        with open('onetrichprime.com', 'w') as f:
            f.write(onetrichprime)
        onetrichprime = rxfile.File('onetrichprime')
        if not nocalc:
            newrun(onetrichprime,hprimetail(itnlcordL, hprimevdwtail))
            #onetrichprime.com.rung09()
            #onetrichprime.com.isover()
            #onetrichprime.runformchk()
        onetrichprime.fchk.read()
        calcphfgroup(onetricL, onetrichprime, qmfchk)

    if onetriucL:
        onetriucLhprime = hprimehead + newhprimetail(itnlcordL, hprimevdwtail)
        with open('onetriucLhprime.com', 'w') as f:
            f.write(onetriucLhprime)
        onetriucLhprime = rxfile.File('onetriucLhprime')
        if not nocalc:
            newrun(onetriucLhprime,hprimetail(itnlcordL, hprimevdwtail))
            #onetriucLhprime.com.rung09()
            #onetriucLhprime.com.isover()
            #onetriucLhprime.runformchk()
        onetriucLhprime.fchk.read()
        calcphfgroup(onetriucL, onetriucLhprime, qmfchk)

    if onetwoL:
        onetwoLhprime = hprimehead + newhprimetail(itnlcordL, hprimevdwtail)
        with open('onetwoLhprime.com', 'w') as f:
            f.write(onetwoLhprime)
        onetwoLhprime = rxfile.File('onetwoLhprime')
        if not nocalc:
            newrun(onetwoLhprime,hprimetail(itnlcordL, hprimevdwtail))
            #onetwoLhprime.com.rung09()
            #onetwoLhprime.com.isover()
            #onetwoLhprime.runformchk()
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
    try:
        hprime = onefourhprime
    except:
        try:
            hprime = onetrifourhprime
        except:
            try:
                hprime = onetrichprime
            except:
                try:
                    hprime = onetriucLhprime
                except:
                    hprime = onetwoLhprime

    hideal = qmfchk.fchk.hessian - hprime.fchk.hessian
    for item in unkparmL:
        tmp = item.hessfile.fchk.hessian
        leftL.append([x for x in tmp])

    leftL = list(zip(*leftL))

    results = np.linalg.lstsq(leftL, hideal)[0]
    for i, item in enumerate(unkparmL):
        item.forceconst = results[i]
    fhffilename = summarize(unkL, itnlcordL, originalname, finalhead, 'fhf')


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

# BHF:
    pairs = []
    pairs.extend(onetwoL.keys())
    pairs.extend(onetricL.keys())
    pairs.extend(onetriucL.keys())
    pairs.extend(onefourL.keys())
    pairs.extend(onetrifourL.keys())
    for i in range(1,mole.natoms+1):
        pairs.append(str(i)+'-'+str(i))

    def getBHessian(fileobj):
        #return fileobj.fchk.hessian
        hessian = []
        for item in pairs:
            a,b = item.split('-')
            a = int(a)
            b = int(b)
            h = fileobj.fchk.find33Hessian(a,b)
            if a!=b:
            #if True:
                for item in h.tolist():
                    for i in item:
                        hessian.append(i)
            else:
                 pass
                 hessian.append(0.5*h[0][0])
                 hessian.append(h[1][0])
                 hessian.append(0.5*h[1][1])
                 hessian.append(h[2][0])
                 hessian.append(h[2][1])
                 hessian.append(0.5*h[2][2])
        return hessian

    def wfhf(fhffile,deltafunc,omegafunc,name):

        # with open(fhffile,'r') as f:
        #     content = f.read()
        # with open(fhffile,'w') as f:
        #     chkname = os.path.splitext(fhffile)[0]+'.chk\n'
        #     f.write(chkname+content)
        # fhffile = rxfile.File(os.path.splitext(fhffile)[0])
        # fhffile.com.rung09()
        # fhffile.com.isover()
        # fhffile.runformchk()
        # fhffile.fchk.read()
        # fhfhessian = np.array(getBHessian(fhffile))

        leftL = []
        hideal = []
        try:
            hprime = onefourhprime
        except:
            try:
                hprime = onetrifourhprime
            except:
                try:
                    hprime = onetrichprime
                except:
                    try:
                        hprime = onetriucLhprime
                    except:
                        hprime = onetwoLhprime
        qmhess = np.array(getBHessian(qmfchk))
        hprimehess = np.array(getBHessian(hprime))

        hideal = qmhess - hprimehess
#       deltaij = deltafunc(fhfhessian,qmhess)
#        omegaij = omegafunc(deltaij)

        for item in unkparmL:
            tmp = np.array(getBHessian(item.hessfile))
            leftL.append([x for x in tmp])

        leftL = np.array(list(zip(*leftL)))
#        hideal = np.array([x*y for x,y in zip(hideal, omegaij)])
#       leftL = np.array([x*y for x,y in zip(omegaij, leftL)])

        results = np.linalg.lstsq(leftL, hideal)[0]

        for i, item in enumerate(unkparmL):
            item.forceconst = results[i]
        summarize(unkL, itnlcordL, originalname, finalhead, name)
        # End 
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

        return results


    # abs diff
    def delta456(mmhess,qmhess):
        return np.array([abs(x-y) for x,y in zip(mmhess,qmhess)])
    # squared diff
    def delta7(mmhess,qmhess):
        return np.array([(x-y)**2 for x,y in zip(mmhess,qmhess)])
    # 1 - x/max
    def omega1(deltaij):
        ma = max(deltaij)
        return np.array([1-x/ma for x in deltaij])
    # 1- x/summ
    def omega47(deltaij):
        summ = sum(deltaij)
        return np.array([(1-x/summ) for x in deltaij])
    # 1- (x/summ)^2
    def omega5(deltaij):
        summ = sum(deltaij)
        return np.array([(1-(x/summ)**2) for x in deltaij])
    # (1-x/summ)^2
    def omega6(deltaij):
        summ = sum(deltaij)
        return np.array([(1-x/summ)**2 for x in deltaij])
    # exp(-lambda*delta_ij/Delta_av)
    def omega8(deltaij):
        av = sum(deltaij)/len(deltaij)
        return np.array([np.exp(-x/av) for x in deltaij])
    def omega0(deltaij):
        return np.array([1 for x in deltaij])


    # method 4:
    #wfhf(fhffilename,delta456,omega0,'bhf')
    #wfhf(fhffilename,delta456,omega1,'wbhf1')
    #wfhf(fhffilename,delta456,omega47,'wfhf4')
    #wfhf(fhffilename,delta456,omega5,'wfhf5')
    #wfhf(fhffilename,delta456,omega6,'wfhf6')
    #wfhf(fhffilename,delta7,omega47,'wfhf7')
    #wfhf(fhffilename,delta7,omega8,'wbhf8')
    wfhf(fhffilename,delta456,omega0,'rfhf')

    # # do iter of wFHF
    # i = 0
    # while True:
    #     i+=1
    #     print(i)
    #     results = np.array(wfhfiter(fhffilename))
    #     if i==1 :
    #         oldresults = results
    #         continue
    #     else:
    #         diff = oldresults - results
    #     oldresults = results
    #     length = np.linalg.norm(diff)
    #     print(length)
    #     if length == 0.000:
    #         print(results)
    #         break


    # Start of IHF

    def readintcoords(fileobj):
        intcords = []
        with open(fileobj.logname) as f:
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
                intcords.append(InternalCoordinates(fileobj, line))
        return intcords

    # Check internal coordinates order matching
    qmfchk.itnl = readintcoords(qmfchk)
    assert len(qmfchk.itnl) == len(itnlcordL)
    for item in hess:
        item.itnl = readintcoords(item)
        assert len(item.itnl) == len(itnlcordL)
        for index, this in enumerate(item.itnl):
            assert this.atomset == qmfchk.itnl[index].atomset
    hprime.itnl = readintcoords(hprime)
    assert len(hprime.itnl) == len(itnlcordL)

    for x, y in zip(hprime.itnl, qmfchk.itnl):
        assert x.atomset == y.atomset

    # Read internal Hessian
    for i in range(1, len(qmfchk.itnl) + 1):
        qmfchk.itnl[i - 1].value = qmfchk.fchk.findintHessianElement(i, i)
        hprime.itnl[i - 1].value = (hprime.fchk.findintHessianElement(i, i))
        for item in hess:
            item.itnl[i - 1].value = item.fchk.findintHessianElement(i, i)
    # Build leftL in itnl order:

    leftL = []
    for item in unkparmL:
        leftL.append([])
        tmp = item.hessfile.itnl
        leftL[-1].extend(x.value for x in tmp)

    leftL = list(zip(*leftL))
    # Build b  vector in itnl order:
    rightL = []
    for qm, hp in zip(qmfchk.itnl, hprime.itnl):
        rightL.append(qm.value - hp.value)
    rightL = np.array(rightL)

    leftL = np.array(leftL)
    rightL = np.array(rightL)
    res = np.linalg.lstsq(leftL, rightL)[0]

    for i, item in enumerate(unkparmL):
        item.forceconst = res[i]
    summarize(unkL, itnlcordL, originalname, finalhead, 'ihf')

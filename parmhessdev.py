#!/usr/bin/env python3
import argparse
import logging
import itertools
import copy
import os
import time
import shutil
import pdb
from multiprocessing import Pool
import numpy as np
from io import StringIO
import rxcclib.Geometry.molecules as rxmol
import rxcclib.File.chemfiles as rxfile
import rxcclib.utils.utils as utils
from rxcclib.Geometry.InternalCoordinates import InternalCoordinates
from rxcclib.File.GauAmberCOM import DihdFunction
from rxcclib.File.GauAmberCOM import MMFunction
from rxcclib.File.GauAmberCOM import GauAmberCOM


class HessFile(rxfile.File):
    def __init__(self, filename, itnl):
        super().__init__(filename)
        self.itnl = itnl
        itnl.hess = self

    def run(self):
        try:
            self.com.rung09()
            self.com.isover()
            self.runformchk()
        except Exception as e:
            logging.warning('The following error occurred during'
                            'HessFile.run() of {}'.format(self.comname))
            logging.warning(e)
            logging.warning('Try once more:')
            self.com.rung09()
            self.com.isover()
            self.runformchk()
        return True

    def getcoorddivs(self):

        ty = type(self.itnl)
        # if bond or angle
        if ty is rxmol.Bond or ty is rxmol.Angle:
            # Force = 2k(x-x0)*firstdiv
            # firstdiv = Force/2
            self.firstdiv = [x / 2 for x in self.fchk.force]
            self.intfirstdiv = [x / 2 for x in self.fchk.intforce]


            # Hessian = 2k*firstdiv1*firstdiv2+2k(x-x0)*seconddiv
            # seconddiv = (Hessian - 2*firstdiv1*firstdiv2) / 2
            self.seconddiv = []
            for i, x in enumerate(self.fchk.hessian):
                row, column = utils.findrowcolumn(i)
                self.seconddiv.append(x / 2 - self.firstdiv[row] *
                                      self.firstdiv[column])
            self.intseconddiv = []
            for i, x in enumerate(self.fchk.inthessian):
                row, column = utils.findrowcolumn(i)
                self.intseconddiv.append(x / 2 - self.intfirstdiv[row] *
                                         self.intfirstdiv[column])

        elif ty is rxmol.Dihd or ty is rxmol.Improper:
            # Force = -nV*sinX*firstdiv2
            # Firstdiv = Force / -sinX
            self.firstdiv = [x / self.itnl.sinX for x in self.fchk.force]
            self.intfirstdiv = [x / self.itnl.sinX
                                for x in self.fchk.intforce]
            # Hessian = -n^2V*cosX*firstdiv1*firstdiv2 - nV*sinX*seconddiv
            # seconddiv = (Hessian + cosX*firstdiv1*firstdiv2) / (-sinX)
            self.seconddiv = []
            for i, x in enumerate(self.fchk.hessian):
                row, column = utils.findrowcolumn(i)
                self.seconddiv.append((x + self.itnl.cosX * self.firstdiv[row]
                                       * self.firstdiv[column])
                                      / (- self.itnl.sinX))


            self.intseconddiv = []
            for i, x in enumerate(self.fchk.inthessian):
                row, column = utils.findrowcolumn(i)
                self.intseconddiv.append((x + self.itnl.cosX *
                                          self.intfirstdiv[row] *
                                          self.intfirstdiv[column])
                                         / (- self.itnl.sinX))

        # def removesmall(L):
        #     fL=[]
        #     for item in L:
        #         if abs(item) < 1.0e-10:
        #             item = 0.0
        #         fL.append(item)
        #     return fL
        # self.firstdiv = removesmall(self.firstdiv)
        # self.intfirstdiv = removesmall(self.intfirstdiv)
        # self.seconddiv = removesmall(self.seconddiv)
        # self.intseconddiv = removesmall(self.intseconddiv)

    def recover(self, n=0):
        # Calc Forces: 0 for harmonic, -n*sin(n*\phi)*firstdiv (phase = 0)
        # Calc Hessian: 2*1stdvs1*1stdvs2
        #      -n^2*cos(n*\phi)*1stdvs1*1stdvs2 -nsin(n*\phi)*2nddvs12
        ty = type(self.itnl)
        if ty is rxmol.Bond or ty is rxmol.Angle:
            self.realforce = [0 for x in self.firstdiv]
            self.realintforce = [0 for x in self.intfirstdiv]

            self.realhessian = []
            self.realinthessian = []
            for i, x in enumerate(self.fchk.hessian):
                row, column = utils.findrowcolumn(i)
                self.realhessian.append(
                    2 * self.firstdiv[row] *
                    self.firstdiv[column]
                )
            for i, x in enumerate(self.fchk.inthessian):
                row, column = utils.findrowcolumn(i)
                self.realinthessian.append(
                    2 * self.intfirstdiv[row] *
                    self.intfirstdiv[column]
                )

        elif ty is rxmol.Dihd or ty is rxmol.Improper:
            sinX = np.sin((n*self.itnl.anglevalue - 180)*np.pi/180)
            cosX = np.cos((n*self.itnl.anglevalue - 180)*np.pi/180)
            self.realforce = [-n * sinX * x for x in self.firstdiv]
            self.realintforce = [-n * sinX * x
                                 for x in self.intfirstdiv]

            self.realhessian = []
            self.realinthessian = []
            for i, x in enumerate(self.fchk.hessian):
                row, column = utils.findrowcolumn(i)
                self.realhessian.append(
                    (-n**2 *
                     cosX * self.firstdiv[row] *
                     self.firstdiv[column] -
                     n * sinX * self.seconddiv[i])
                )
            for i, x in enumerate(self.fchk.inthessian):
                row, column = utils.findrowcolumn(i)
                self.realinthessian.append(
                    (-n**2 *
                     cosX * self.intfirstdiv[row] *
                     self.intfirstdiv[column] -
                     n * sinX * self.intseconddiv[i])
                )
        else:
            raise


class GlobalSetting(object):
    pass


def initset():
    # parse input
    parser = argparse.ArgumentParser()
    parser.add_argument('inputinp',
                        default='input.inp',
                        help='input.inp prepared by tsubasa')
    parser.add_argument('--quiet',
                        '-q',
                        action='store_true',
                        help='remove screen message')
    parser.add_argument('--nocalc',
                        '-nc',
                        action='store_true',
                        help='Use alreday calculated file')
    parser.add_argument('--debug',
                        action='store_true',
                        help='Debug mode')
    args = parser.parse_args()

    with open(args.inputinp, 'r') as f:
        for line in f:
            if line.find('mmfile') >= 0:
                mmfile = line.split('=')[1].strip('\n')
                mmfile = os.path.splitext(mmfile)[0]
            elif line.find('qmfchk') >= 0:
                qmfile = line.split('=')[1].strip('\n')
                qmfile = os.path.splitext(qmfile)[0]

    GlobalSetting.debug = args.debug
    GlobalSetting.nocalc = args.nocalc

    # logging module
    # # to file
    logging.basicConfig(filename=mmfile + '.info',
                        level=logging.DEBUG,
                        filemode='w')
    # # to screen
    if not args.quiet:
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        formatter = logging.Formatter('%(levelname)-8s %(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)

    # Welcome info
    logging.info('Parmhess for Amber Parameterization')
    logging.info('   An implementation of PHF/FHF/IHF method.')
    logging.info('   If you use this program and/or Hessian Fitting methods,'
                 ' cite DOI: 10.1002/jcc.24457\n\n')
    logging.info('Input details:')
    logging.info('   MMFile: ' + mmfile)
    logging.info('   QMFile: ' + qmfile)
    logging.info('   NOCalc: ' + str(GlobalSetting.nocalc))
    logging.info('   Quiet : ' + str(args.quiet) + '\n\n')

    # change working directory
    if GlobalSetting.nocalc is True:
        logging.info('NOCalc is True; So use old files without recalc.')
    if args.quiet is True:
        logging.info('Quiet  is True; So no screen message.')

    existance = os.path.lexists('hffiles')
    if existance and GlobalSetting.nocalc:
        os.chdir('hffiles')
    elif existance and not GlobalSetting.nocalc:
        shutil.rmtree('hffiles')
        os.mkdir('hffiles')
        os.chdir('hffiles')
        shutil.copy(os.path.join('..', mmfile + '.com'), '.')
        shutil.copy(os.path.join('..', qmfile + '.fchk'), '.')
        shutil.copy(os.path.join('..', qmfile + '.log'), '.')
    elif not existance and not GlobalSetting.nocalc:
        os.mkdir('hffiles')
        os.chdir('hffiles')
        shutil.copy(os.path.join('..', mmfile + '.com'), '.')
        shutil.copy(os.path.join('..', qmfile + '.fchk'), '.')
        shutil.copy(os.path.join('..', qmfile + '.log'), '.')
    elif not existance and GlobalSetting.nocalc:
        logging.critical('NOCalc is requested but no "hffiles" folder exists.')
    else:
        raise

    GlobalSetting.mmfile = rxfile.File(mmfile)
    GlobalSetting.qmfile = rxfile.File(qmfile)
    GlobalSetting.qmfile.fchk.read()


def readgeom():
    # Read info from com
    mole = rxmol.Molecule('thisgeometry')
    GlobalSetting.mole = mole
    mmcom = GauAmberCOM(GlobalSetting.mmfile)
    mmcom.read()
    mole.readfromxyz(mmcom.xyz)
    mole.readchargefromlist(mmcom.atomchargelist)
    mole.readtypefromlist(mmcom.atomtypelist)
    mole.readconnectivity(mmcom.connectivity)
    for atom in mole:
        atom.vdwradius = mmcom.vdwdict[atom.atomtype][0]
        atom.vdwwelldepth = mmcom.vdwdict[atom.atomtype][1]

    # Store and count finalfunc
    finalfuncL = []
    finalfuncL.extend(sorted(mmcom.bondfunc, key=lambda x: x.repr))
    finalfuncL.extend(sorted(mmcom.anglefunc, key=lambda x: x.repr))
    finalfuncL.extend(sorted(mmcom.dihdfunc, key=lambda x: x.repr))
    finalfuncL.extend(sorted(mmcom.improperfunc, key=lambda x: x.repr))

    for item in finalfuncL:
        if item.type == 'dihd':
            item.dihdunkterms = []
            item.known = True
            for func in item.dihdfunctions:
                if func.forceconst == MMFunction.unknownsign:
                    item.known = False
                    item.dihdunkterms.append(func.periodicity)

        elif item.type == 'improper':
            # # if improper, find and do mole.addimproper  here
            for atom3 in mole:
                if atom3.atomtype == item.c:
                    permu = itertools.permutations(atom3.neighbor, 3)
                    res = []
                    for tu in permu:
                        a = tu[0].atomtype == item.a or item.a == '*'
                        b = tu[1].atomtype == item.b or item.b == '*'
                        c = tu[2].atomtype == item.d or item.d == '*'
                        if a and b and c:
                            res.append([tu[0].atomnum, tu[1].atomnum,
                                        atom3.atomnum, tu[2].atomnum])
                    res = sorted(res, key=lambda x: (str(x[1]) + str(x[3])))
                    res = res[0]
                    mole.addimproper(*res)
            if item.forceconst == MMFunction.unknownsign:
                item.known = False
            else:
                item.known = True
        # # Bond and Angle
        else:
            if item.forceconst == MMFunction.unknownsign:
                item.known = False
            else:
                item.known = True

    # Match itnl and finalfunc
    itnlcordL = []
    itnlcordL.extend(sorted(mole.dihdlist.values(), key=lambda x: x.repr))
    itnlcordL.extend(sorted(mole.anglelist.values(), key=lambda x: x.repr))
    itnlcordL.extend(sorted(mole.bondlist.values(), key=lambda x: x.repr))
    itnlcordL.extend(sorted(mole.improperlist.values(), key=lambda x: x.repr))

    unkitnlL = []
    knownitnlL = []
    for item in itnlcordL:
        for func in finalfuncL:
            if matchitnlwithfinalfunc(item, func):
                item.func = func
                if type(item) is rxmol.Dihd:
                    item.dihdfunctions = copy.deepcopy(func.dihdfunctions)
                    item.npaths = func.npaths
                elif type(item) is rxmol.Improper:
                    item.forceconst = func.forceconst
                    item.phase = func.phase
                    item.periodicity = func.periodicity
                else:
                    item.forceconst = func.forceconst
                item.known = func.known
                break
        try:
            if item.known is False:
                unkitnlL.append(item)
            elif item.known is True:
                knownitnlL.append(item)
        except AttributeError as e:
            logging.critical(e)
            logging.critical('Matching finalfunc was unsuccessful {} '.format(item.repr))

    return finalfuncL, itnlcordL, unkitnlL, knownitnlL


def matchitnlwithfinalfunc(item, finalfunc):
    if type(item) is rxmol.Dihd and finalfunc.type == 'dihd':
        a = (item[1].atomtype == finalfunc.a or finalfunc.a == '*')
        b = (item[2].atomtype == finalfunc.b or finalfunc.b == '*')
        c = (item[3].atomtype == finalfunc.c or finalfunc.c == '*')
        d = (item[4].atomtype == finalfunc.d or finalfunc.d == '*')
        forward = a and b and c and d
        a = (item[1].atomtype == finalfunc.d or finalfunc.d == '*')
        b = (item[2].atomtype == finalfunc.c or finalfunc.c == '*')
        c = (item[3].atomtype == finalfunc.b or finalfunc.b == '*')
        d = (item[4].atomtype == finalfunc.a or finalfunc.a == '*')
        backward = a and b and c and d
        if forward or backward:
            return True
        else:
            return False
    elif type(item) is rxmol.Angle and finalfunc.type == 'angle':
        a = (item[1].atomtype == finalfunc.a or finalfunc.a == '*')
        b = (item[2].atomtype == finalfunc.b or finalfunc.b == '*')
        c = (item[3].atomtype == finalfunc.c or finalfunc.c == '*')
        forward = a and b and c
        a = (item[1].atomtype == finalfunc.c or finalfunc.c == '*')
        b = (item[2].atomtype == finalfunc.b or finalfunc.b == '*')
        c = (item[3].atomtype == finalfunc.a or finalfunc.a == '*')
        backward = a and b and c
        if forward or backward:
            return True
        else:
            return False
    elif type(item) is rxmol.Bond and finalfunc.type == 'bond':
        a = (item[1].atomtype == finalfunc.a or finalfunc.a == '*')
        b = (item[2].atomtype == finalfunc.b or finalfunc.b == '*')
        forward = a and b
        a = (item[1].atomtype == finalfunc.b or finalfunc.b == '*')
        b = (item[2].atomtype == finalfunc.a or finalfunc.a == '*')
        backward = a and b
        if forward or backward:
            return True
        else:
            return False
    elif type(item) is rxmol.Improper and finalfunc.type == 'improper':
        a = (item[1].atomtype == finalfunc.a or finalfunc.a == '*')
        b = (item[2].atomtype == finalfunc.b or finalfunc.b == '*')
        c = (item[3].atomtype == finalfunc.c or finalfunc.c == '*')
        d = (item[4].atomtype == finalfunc.d or finalfunc.d == '*')
        forward = a and b and c and d
        if forward:
            return True
        else:
            return False
    else:
        return False


def prepHeadAndVDW():
    mole = GlobalSetting.mole
    mmcom = GlobalSetting.mmfile.com
    qmfile = GlobalSetting.qmfile
    hessxyz = ''
    hprimexyz = ''
    finalxyz = ''
    hessvdwtail = 'Nonbon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 -1.2\n'
    hprimevdwtail = 'Nonbon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 -1.2\n'
    finalvdwtail = 'Nonbon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 -1.2\n'
    for atom in mole:
        hessxyz += '{}-{}-{}     {}\n'.format(
            atom.elementsym,
            atom.name,
            '0.000000',
            '   '.join(["{: .12f}".format(x) for x in atom.coords])
        )
        hprimexyz += '{}-{}-{:8.6f}     {}\n'.format(
            atom.elementsym,
            atom.name,
            float(atom.atomcharge),
            '   '.join(["{: .12f}".format(x) for x in atom.coords])
        )
        finalxyz += '{}-{}-{:8.6f}     {}\n'.format(
            atom.elementsym,
            atom.atomtype,
            float(atom.atomcharge),
            '   '.join(["{: .12f}".format(x) for x in atom.coords])
        )

        hessvdwtail += 'VDW  {}  {}  0.0000\n'.format(atom.name,
                                                      atom.vdwradius)
        hprimevdwtail += 'VDW  {}  {}  {}\n'.format(atom.name,
                                                    atom.vdwradius,
                                                    atom.vdwwelldepth)
    for key, value in mmcom.vdwdict.items():
        finalvdwtail += 'VDW  {}  {}  {}\n'.format(key,
                                                   value[0],
                                                   value[1])

    headL = ['', '', '']
    xyzL = [hessxyz, hprimexyz, finalxyz]
    for i in range(3):
        headL[i] = ('{}\n'
                    'hess\n\n'
                    '{} {} \n'
                    '{}\n'
                    '{}\n'.format(
                        mmcom.route,
                        qmfile.totalcharge,
                        qmfile.multiplicity,
                        xyzL[i],
                        mmcom.connectivity
                    ))
    return headL, hessvdwtail, hprimevdwtail, finalvdwtail


def preparehesstail(thisitnl, unkitnlL, hessaddtail):
    resstr = ''
    for item in unkitnlL:
        if type(item) is rxmol.Dihd:
            if item is thisitnl:
                parm = ['627.5095', '0.000', '0.000', '0.000']
            else:
                parm = ['0.000', '0.000', '0.000', '0.000']
                # Gaussian can ignore missing dihedral
                # So just ignore all zero terms.
                continue
            # cosine is made close to zero for better
            # precision of 2nd coordinate derivatives
            # set cos(n*\phi - \gamma) = 0 and sin = 1
            # \gamma = n*\phi + 90
            # 3600: to get a minimum positive number
            # Gaussian do not accept float number for phase
            # so use a near integer.
            thisphase = int((item.anglevalue + 3510) % 360)
            phase = [str(thisphase), '0', '0', '0']
            # and record the actual cos/sin value.
            X = (item.anglevalue - thisphase)
            item.cosX = np.cos(X * np.pi / 180)
            item.sinX = np.sin(X * np.pi / 180)

            resstr += 'AmbTrs  {}  {}  {}  1.0\n'.format(
                ' '.join([x.center(3, ' ') for x in item.repr.split()]),
                ' '.join([x.center(3, ' ') for x in phase]),
                ' '.join([x.center(3, ' ') for x in parm]))
        else:
            if item is thisitnl:
                if type(item) is rxmol.Bond:
                    parm = '2240.877122994058'  # 627.5095 / bohr^2
                else:
                    parm = '627.5095'
            else:
                parm = '0.000'
            if type(item) is rxmol.Angle:
                # set (\theta - \theta_0) = 1
                resstr += 'HrmBnd1  {}  {}  {:>9.5f}\n'.format(
                    ' '.join([x.center(3, ' ') for x in item.repr.split()]),
                    parm, (item.anglevalue + 180 / np.pi))
            elif type(item) is rxmol.Bond:
                resstr += 'HrmStr1  {}  {}  {:>7.5f}\n'.format(
                    ' '.join([x.center(3, ' ') for x in item.repr.split()]),
                    parm, (item.length + 0.5291772086))
            elif type(item) is rxmol.Improper:
                phase = '90.0'
                X = (item.periodicity * item.anglevalue - 90)
                item.cosX = np.cos(X * np.pi / 180)
                item.sinX = np.sin(X * np.pi / 180)

                resstr += 'ImpTrs  {}  {}  {}  {}\n'.format(
                    ' '.join([x.center(3, ' ') for x in item.repr.split()]),
                    parm, phase, item.periodicity)
    return resstr + hessaddtail + '\n\n'


def preorihesstail(thisitnl, itnlcordL, hesstail):
    # prepare original HessTail for comparasion
    resstr = ''
    for item in itnlcordL:
        if type(item) is rxmol.Dihd:
            if item is thisitnl:
                parm = ['0.000', '1.000', '0.000', '0.000']
            else:
                # Gaussian can ignore missing dihedral
                # So just ignore all zero terms.
                continue

            phase = ['0', '180', '0', '0']

            resstr += 'AmbTrs  {}  {}  {}  1.0\n'.format(
                ' '.join([x.center(3, ' ') for x in item.repr.split()]),
                ' '.join([x.center(3, ' ') for x in phase]),
                ' '.join([x.center(3, ' ') for x in parm]))
        else:
            if item is thisitnl:
                parm = '1.000'
            else:
                parm = '0.000'
            if type(item) is rxmol.Angle:
                # set (\theta - \theta_0) = 1
                resstr += 'HrmBnd1  {}  {}  {:>9.5f}\n'.format(
                    ' '.join([x.center(3, ' ') for x in item.repr.split()]),
                    parm, (item.anglevalue))
            elif type(item) is rxmol.Bond:
                resstr += 'HrmStr1  {}  {}  {:>7.5f}\n'.format(
                    ' '.join([x.center(3, ' ') for x in item.repr.split()]),
                    parm, (item.length))
            elif type(item) is rxmol.Improper:
                phase = '180.0'
                resstr += 'ImpTrs  {}  {}  {}  {}\n'.format(
                    ' '.join([x.center(3, ' ') for x in item.repr.split()]),
                    parm, phase, item.periodicity)

    return resstr + hesstail + '\n\n'


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
            content += (str(item[1].atomnum) + ' ' +
                        str(item[2].atomnum) + ' A\n')
        if type(item) == rxmol.Angle:
            content += (str(item[1].atomnum) + ' ' +
                        str(item[2].atomnum) + ' ' +
                        str(item[3].atomnum) + ' A\n')
        if type(item) == rxmol.Dihd:
            content += (str(item[1].atomnum) + ' ' +
                        str(item[2].atomnum) + ' ' +
                        str(item[3].atomnum) + ' ' +
                        str(item[4].atomnum) + ' A\n')
        if type(item) == rxmol.Improper:
            content += (str(item[1].atomnum) + ' ' +
                        str(item[2].atomnum) + ' ' +
                        str(item[3].atomnum) + ' ' +
                        str(item[4].atomnum) + ' A\n')
    content += '\n'
    return content


def makehess(head, tail, unkitnlL, itnlcordL):
    fileL = []
    if GlobalSetting.nocalc is False:
        for item in unkitnlL:
            fileL.append(HessFile('hess' + str(len(fileL)), item))
            with open(item.hess.comname, 'w') as f:
                f.write('%chk= {}\n{}{}{}'.format(
                    item.hess.chkname,
                    head,
                    preorihesstail(item, itnlcordL, tail),
                    addlink1(item.hess, itnlcordL)
                ))
    else:
        for item in unkitnlL:
            fileL.append(HessFile('hess' + str(len(fileL)), item))
            preorihesstail(item, itnlcordL, tail)

    return fileL


def makehprime(head, tail, itnlcordL):
    for item in itnlcordL:
        if type(item) is rxmol.Dihd:
            if item.forceconst == MMfunction.unknownsign:
                parm = '0.0'
        else:


    pass


def cleanup(itnlcordL, finalfuncL):
    for item in itnlcordL:
        if type(item) != rxmol.Dihd:
            if item.known is False:
                item.forceconst = MMFunction.unknownsign
        else:
            for parms in item.dihdfunctions:
                if parms.known is False:
                    parms.forceconst = MMFunction.unknownsign

    for item in finalfuncL:
        if item.type != 'dihd':
            if item.known is False:
                item.forceconst = MMFunction.unknownsign
        else:
            for parms in item.dihdfunctions:
                if parms.known is False:
                    parms.forceconst = MMFunction.unknownsign


def runhess(hessL, para=4):
    if GlobalSetting.nocalc is False:
        p = Pool(para)
        p.map(parallel_runhess, hessL)
    for x in hessL:
        x.fchk.read()
        # x.getcoorddivs()
        # x.recover(n=2)
    return


def parallel_runhess(x):
    x.run()
    return


if __name__ == '__main__':
    # read inputs, set logging module
    initset()

    # read info, identify itnls, funcs
    finalfuncL, itnlcordL, unkitnlL, knownitnlL = readgeom()

    # prepare head for hess ,hprime and finalresult
    heads, hessvdwtail, hprimevdwtail, finalvdwtail = prepHeadAndVDW()

    # prepare hess, inithprime:
    hessL = makehess(heads[0], hessvdwtail, unkitnlL, itnlcordL)

    # run hess
    runhess(hessL, para=4)

    thismole = GlobalSetting.mole


    # for item in knownitnlL:
    #     if type(item) is rxmol.Dihd:
    #         phase = [x.phase for x in item.dihdfunctions]
    #         parm = [x.forceconst for x in item.dihdfunctions]
    #         npaths = item.npaths
    #         resstr += 'AmbTrs  {}  {}  {}  {}\n'.format(
    #             ' '.join([x.center(3, ' ') for x in item.repr.split()]),
    #             ' '.join([x.center(3, ' ') for x in phase]),
    #             ' '.join([x.center(3, ' ') for x in parm]),
    #             npaths)
    #     else:
    #         parm = item.forceconst
    #         if type(item) is rxmol.Angle:
    #             # set (\theta - \theta_0) = 1
    #             resstr += 'HrmBnd1  {}  {}  {:>9.5f}\n'.format(
    #                 ' '.join([x.center(3, ' ') for x in item.repr.split()]),
    #                 parm, (item.anglevalue))
    #         elif type(item) is rxmol.Bond:
    #             resstr += 'HrmStr1  {}  {}  {:>7.5f}\n'.format(
    #                 ' '.join([x.center(3, ' ') for x in item.repr.split()]),
    #                 parm, (item.length))
    #         elif type(item) is rxmol.Improper:
    #             phase = '180.0'
    #             resstr += 'ImpTrs  {}  {}  {}  {}\n'.format(
    #                 ' '.join([x.center(3, ' ') for x in item.repr.split()]),
    #                 parm, phase, item.periodicity)





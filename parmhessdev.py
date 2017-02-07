#!/usr/bin/env python3
import argparse
import logging
import itertools
import copy
import os
import shutil
from multiprocessing import Pool
from io import StringIO
import numpy as np
import rxcclib.Geometry.molecules as rxmol
import rxcclib.File.chemfiles as rxfile
import rxcclib.utils.utils as utils
from rxcclib.File.GauAmberCOM import MMFunction
from rxcclib.File.GauAmberCOM import GauAmberCOM


class HessFile(rxfile.File):
    def __init__(self, filename, itnl):
        super().__init__(filename)
        self.itnl = itnl
        itnl.hess = self
        self.mymolecule = self.itnl[1].mymolecule

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
            self.itnl.firstdiv = [x / 2 for x in self.fchk.force]
            self.itnl.intfirstdiv = [x / 2 for x in self.fchk.intforce]

            # Hessian = 2k*firstdiv1*firstdiv2+2k(x-x0)*seconddiv
            # seconddiv = (Hessian - 2*firstdiv1*firstdiv2) / 2
            self.itnl.seconddiv = []
            for i, x in enumerate(self.fchk.hessian):
                row, column = utils.find_rowcolumn_of_LTri(i)
                self.itnl.seconddiv.append(x / 2 - self.itnl.firstdiv[row] *
                                           self.itnl.firstdiv[column])
            self.itnl.intseconddiv = []
            for i, x in enumerate(self.fchk.inthessian):
                row, column = utils.find_rowcolumn_of_LTri(i)
                self.itnl.intseconddiv.append(x / 2 - self.itnl.intfirstdiv[
                    row] * self.itnl.intfirstdiv[column])

        elif ty is rxmol.Dihd or ty is rxmol.Improper:
            # Force = -nV*sinX*firstdiv2
            # Firstdiv = Force / -sinX
            self.itnl.firstdiv = [x / self.itnl.sinX for x in self.fchk.force]
            self.itnl.intfirstdiv = [
                x / self.itnl.sinX for x in self.fchk.intforce
            ]
            # Hessian = -n^2V*cosX*firstdiv1*firstdiv2 - nV*sinX*seconddiv
            # seconddiv = (Hessian + cosX*firstdiv1*firstdiv2) / (-sinX)
            self.itnl.seconddiv = []
            for i, x in enumerate(self.fchk.hessian):
                row, column = utils.find_rowcolumn_of_LTri(i)
                self.itnl.seconddiv.append(
                    (x + self.itnl.cosX * self.itnl.firstdiv[row] *
                     self.itnl.firstdiv[column]) / (-self.itnl.sinX))

            self.itnl.intseconddiv = []
            for i, x in enumerate(self.fchk.inthessian):
                row, column = utils.find_rowcolumn_of_LTri(i)
                self.itnl.intseconddiv.append(
                    (x + self.itnl.cosX * self.itnl.intfirstdiv[row] *
                     self.itnl.intfirstdiv[column]) / (-self.itnl.sinX))

        # def removesmall(L):
        #     fL=[]
        #     for item in L:
        #         if abs(item) < 1.0e-10:
        #             item = 0.0
        #         fL.append(item)
        #     return fL
        # self.itnl.firstdiv = removesmall(self.itnl.firstdiv)
        # self.itnl.intfirstdiv = removesmall(self.itnl.intfirstdiv)
        # self.itnl.seconddiv = removesmall(self.itnl.seconddiv)
        # self.itnl.intseconddiv = removesmall(self.itnl.intseconddiv)

    def recover(self, n=0):
        # Calc Forces: 0 for harmonic, -n*sin(n*\phi)*firstdiv (phase = 0)
        # Calc Hessian: 2*1stdvs1*1stdvs2
        #      -n^2*cos(n*\phi)*1stdvs1*1stdvs2 -nsin(n*\phi)*2nddvs12
        ty = type(self.itnl)
        if ty is rxmol.Bond or ty is rxmol.Angle:
            self.realforce = [0 for x in self.itnl.firstdiv]
            self.realintforce = [0 for x in self.itnl.intfirstdiv]

            self.realhessian = []
            self.realinthessian = []
            for i, x in enumerate(self.fchk.hessian):
                row, column = utils.find_rowcolumn_of_LTri(i)
                self.realhessian.append(2 * self.itnl.firstdiv[row] *
                                        self.itnl.firstdiv[column])
            for i, x in enumerate(self.fchk.inthessian):
                row, column = utils.find_rowcolumn_of_LTri(i)
                self.realinthessian.append(2 * self.itnl.intfirstdiv[row] *
                                           self.itnl.intfirstdiv[column])

        elif ty is rxmol.Dihd or ty is rxmol.Improper:
            sinX = np.sin((n * self.itnl.anglevalue - 180) * np.pi / 180)

            cosX = np.cos((n * self.itnl.anglevalue - 180) * np.pi / 180)
            self.realforce = [-n * sinX * x for x in self.itnl.firstdiv]
            self.realintforce = [-n * sinX * x for x in self.itnl.intfirstdiv]

            self.realhessian = []
            self.realinthessian = []
            for i, x in enumerate(self.fchk.hessian):
                row, column = utils.find_rowcolumn_of_LTri(i)
                self.realhessian.append((-n**2 * cosX * self.itnl.firstdiv[row]
                                         * self.itnl.firstdiv[column] - n *
                                         sinX * self.itnl.seconddiv[i]))
            for i, x in enumerate(self.fchk.inthessian):
                row, column = utils.find_rowcolumn_of_LTri(i)
                self.realinthessian.append(
                    (-n**2 * cosX * self.itnl.intfirstdiv[row] *
                     self.itnl.intfirstdiv[column] - n * sinX *
                     self.itnl.intseconddiv[i]))
        else:
            raise


class HprimeFile(rxfile.File):
    def __init__(self, filename):
        super().__init__(filename)

    def run(self):
        try:
            self.com.rung09()
            self.com.isover()
            self.runformchk()
        except Exception as e:
            logging.warning('The following error occurred during'
                            'HprimeFile.run() of {}'.format(self.comname))
            logging.warning(e)
            logging.warning('Try once more:')
            self.com.rung09()
            self.com.isover()
            self.runformchk()
        return True


class GlobalSetting(object):
    pass


def initset():
    # parse input
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'inputinp', default='input.inp', help='input.inp prepared by tsubasa')
    parser.add_argument(
        '--quiet', '-q', action='store_true', help='remove screen message')
    parser.add_argument(
        '--nocalc',
        '-nc',
        action='store_true',
        help='Use alreday calculated file')
    parser.add_argument('--debug', action='store_true', help='Debug mode')
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
    GlobalSetting.mmfile = rxfile.File(mmfile)
    GlobalSetting.qmfile = rxfile.File(qmfile)

    # logging module
    # # to file
    logging.basicConfig(
        filename=mmfile + '.info', level=logging.DEBUG, filemode='w')
    # # to screen
    if not args.quiet:
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        formatter = logging.Formatter('%(levelname)-8s %(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)

    if GlobalSetting.nocalc is True:
        logging.info('NOCalc is True; so use old files without recalc.')
    if args.quiet is True:
        logging.info('Quiet  is True; so no screen message.')

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

    existance = os.path.lexists('hffiles')
    if existance and GlobalSetting.nocalc:
        os.chdir('hffiles')
    elif not GlobalSetting.nocalc:
        if existance:
            shutil.rmtree('hffiles')
        os.mkdir('hffiles')
        os.chdir('hffiles')
        shutil.copy(os.path.join('..', mmfile + '.com'), '.')
        shutil.copy(os.path.join('..', qmfile + '.fchk'), '.')
        shutil.copy(os.path.join('..', qmfile + '.log'), '.')
    elif not existance and GlobalSetting.nocalc:
        logging.critical(
            'NOCalc is requested but "hffiles" folder does not exist.')
    else:
        raise Exception("Unexpected condition at Line initset() - chdir()")


def readgeom():
    # Read info from com
    GlobalSetting.qmfile.fchk.read()
    mole = rxmol.Molecule('thisgeometry')
    GlobalSetting.mole = mole
    mmcom = GauAmberCOM(GlobalSetting.mmfile)
    mmcom.read()
    mole.readfromxyz(mmcom.xyz)
    mole.readchargefromlist(mmcom.atomchargelist)
    mole.readtypefromlist(mmcom.atomtypelist)
    mole.readconnectivity(mmcom.connectivity)
    for atom in mole:
        atom.vdwradius = float(mmcom.vdwdict[atom.atomtype][0])
        atom.vdwwelldepth = float(mmcom.vdwdict[atom.atomtype][1])

    # Store and count finalfunc
    finalfuncL = []
    finalfuncL.extend(sorted(mmcom.bondfunc, key=lambda x: x.repr))
    finalfuncL.extend(sorted(mmcom.anglefunc, key=lambda x: x.repr))
    finalfuncL.extend(sorted(mmcom.dihdfunc, key=lambda x: x.repr))
    finalfuncL.extend(sorted(mmcom.improperfunc, key=lambda x: x.repr))

    for item in finalfuncL:
        if item.type == 'dihd':
            item.known = True
            for func in item.dihdfunctions:
                if func.forceconst == MMFunction.unknownsign:
                    func.known = False
                    item.known = False
                else:
                    func.known = True

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
                            res.append([
                                tu[0].atomnum, tu[1].atomnum, atom3.atomnum,
                                tu[2].atomnum
                            ])
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
            logging.critical('Matching finalfunc was unsuccessful {} '.format(
                item.repr))

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
            atom.elementsym, atom.name, '0.000000',
            '   '.join(["{: .12f}".format(x) for x in atom.coords]))
        hprimexyz += '{}-{}-{:8.6f}     {}\n'.format(
            atom.elementsym, atom.name,
            float(atom.atomcharge),
            '   '.join(["{: .12f}".format(x) for x in atom.coords]))
        finalxyz += '{}-{}-{:8.6f}     {}\n'.format(
            atom.elementsym, atom.atomtype,
            float(atom.atomcharge),
            '   '.join(["{: .12f}".format(x) for x in atom.coords]))

        hessvdwtail += 'VDW  {}  {}  0.0000\n'.format(atom.name,
                                                      atom.vdwradius)
        hprimevdwtail += 'VDW  {}  {:6.4f}  {:6.4f}\n'.format(
            atom.name, atom.vdwradius, atom.vdwwelldepth)
    for key, value in mmcom.vdwdict.items():
        finalvdwtail += 'VDW  {}  {}  {}\n'.format(key, value[0], value[1])

    headL = ['', '', '']
    xyzL = [hessxyz, hprimexyz, finalxyz]
    for i in range(3):
        headL[i] = ('{}\n'
                    'hess\n\n'
                    '{} {} \n'
                    '{}\n'
                    '{}\n'.format(mmcom.route, qmfile.totalcharge,
                                  qmfile.multiplicity, xyzL[i],
                                  mmcom.connectivity))
    return headL, [hessvdwtail, hprimevdwtail, finalvdwtail]


def preparehesstail(thisitnl, itnlcordL, hessvdwtail):
    resstr = ''
    for item in itnlcordL:
        if type(item) is rxmol.Dihd:
            if item is thisitnl:
                parm = ['627.5095', '0.000', '0.000', '0.000']
            else:
                continue
                parm = ['0.000', '0.000', '0.000', '0.000']
                # Gaussian can ignore missing dihedral
                # So just ignore all zero terms.

            # cosine is made close to zero for better
            # precision of 2nd coordinate derivatives
            # set cos(n*\phi - \gamma) = 0 and sin = 1
            # \gamma = n*\phi + 90
            # 3600: to get a minimum positive number
            # Gaussian do not accept float number for phase
            # so use a near integer.
            thisphase = int((item.anglevalue + 3510) % 360)
            if thisphase > 180:
                thisphase -= 180
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

    return resstr + hessvdwtail + '\n\n'


# def preorihesstail(thisitnl, itnlcordL, hesstail):
#     # prepare original HessTail for comparasion
#     resstr = ''
#     for item in itnlcordL:
#         if type(item) is rxmol.Dihd:
#             if item is thisitnl:
#                 parm = ['0.000', '1.000', '0.000', '0.000']
#             else:
#                 # Gaussian can ignore missing dihedral
#                 # So just ignore all zero terms.
#                 continue

#             phase = ['0', '180', '0', '0']

#             resstr += 'AmbTrs  {}  {}  {}  1.0\n'.format(
#                 ' '.join([x.center(3, ' ') for x in item.repr.split()]),
#                 ' '.join([x.center(3, ' ') for x in phase]),
#                 ' '.join([x.center(3, ' ') for x in parm]))
#         else:
#             if item is thisitnl:
#                 parm = '1.000'
#             else:
#                 parm = '0.000'
#             if type(item) is rxmol.Angle:
#                 # set (\theta - \theta_0) = 1
#                 resstr += 'HrmBnd1  {}  {}  {:>9.5f}\n'.format(
#                     ' '.join([x.center(3, ' ') for x in item.repr.split()]),
#                     parm, (item.anglevalue))
#             elif type(item) is rxmol.Bond:
#                 resstr += 'HrmStr1  {}  {}  {:>7.5f}\n'.format(
#                     ' '.join([x.center(3, ' ') for x in item.repr.split()]),
#                     parm, (item.length))
#             elif type(item) is rxmol.Improper:
#                 phase = '180.0'
#                 resstr += 'ImpTrs  {}  {}  {}  {}\n'.format(
#                     ' '.join([x.center(3, ' ') for x in item.repr.split()]),
#                     parm, phase, item.periodicity)

#     return resstr + hesstail + '\n\n'


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


def makehess(head, tail, unkitnlL, itnlcordL):
    hessL = []
    if GlobalSetting.nocalc is False:
        for i, item in enumerate(unkitnlL):
            hessL.append(HessFile('hess' + str(i), item))
            with open(item.hess.comname, 'w') as f:
                f.write('%chk= {}\n{}{}{}'.format(
                    item.hess.chkname, head,
                    preparehesstail(item, itnlcordL, tail),
                    addlink1(item.hess, itnlcordL)))
    else:
        for i, item in enumerate(unkitnlL):
            hessL.append(HessFile('hess' + str(i), item))
            preparehesstail(item, itnlcordL, tail)

    return hessL


def preparehprimetail(itnlcordL, hprimevdwtail):
    resstr = ''
    for item in itnlcordL:
        if type(item) is rxmol.Dihd:
            parm = []
            phase = []
            for n in item.dihdfunctions:
                if n.known:
                    parm.append("{:.3f}".format(n.forceconst))
                    phase.append(str(n.phase))
                else:
                    parm.append('0.000')
                    phase.append('0')
            resstr += 'AmbTrs  {}  {}  {}  1.0\n'.format(
                ' '.join([x.center(3, ' ') for x in item.repr.split()]),
                ' '.join([x.center(3, ' ') for x in phase]),
                ' '.join([x.center(3, ' ') for x in parm]))
        else:
            if item.known:
                parm = "{:.3f}".format(item.forceconst)
            else:
                parm = '0.000'
            if type(item) is rxmol.Angle:
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

    return resstr + hprimevdwtail + '\n\n'


def makehprime(head, tail, itnlcordL):
    if GlobalSetting.nocalc is False:
        hprime = HprimeFile("hprime")
        with open(hprime.comname, 'w') as f:
            f.write('%chk={}\n{}{}{}'.format(hprime.chkname, head,
                                             preparehprimetail(itnlcordL,
                                                               tail),
                                             addlink1(hprime, itnlcordL)))
    else:
        hprime = HprimeFile("hprime")
        preparehprimetail(itnlcordL, tail)

    return hprime


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


def runhess(hessL, para=2):
    if GlobalSetting.nocalc is False:
        #        for x in hessL:
        #           x.run()
        p = Pool(para)
        p.map(parallel_runhess, hessL)
    for x in hessL:
        x.fchk.read()
        x.getcoorddivs()
        x.recover(n=2)
    return


def parallel_runhess(x):
    x.run()
    return


def findmatch(itnlcordL):
    def readitnl(fileobj):
        intcords = []
        with open(fileobj.logname) as f:
            tmp = f.read()
        tmp = tmp.split('Initial command')[-1]

        with StringIO(tmp) as f:
            for line in f:
                if line.find('Initial Parameters') < 0:
                    continue
                break
            [next(f) for x in range(0, 4)]
            for line in f:
                if line.find('---') >= 0:
                    break
                line = line.split()[2][2:-1]
                line = [int(x) for x in line.split(',')]
                intcords.append(" ".join([str(x) for x in line]))
        return intcords

    dict = {}
    qmfile = GlobalSetting.qmfile
    gauseq = readitnl(qmfile)

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
        atomset = " ".join([str(x) for x in atomset])
        item.gauseq = gauseq.index(atomset) + 1
    return


def eqforceL(unkitnlL):
    E = []
    for item in unkitnlL:

        if type(
                item
        ) is rxmol.Dihd:  # check unknown numbers; should do this separately to avoid repeated process in eqrightL
            unknum = list(
                filter(lambda x: True if x.forceconst == 'XXXXXX' else False,
                       item.dihdfunctions))
            unknum = len(unknum)
            if unknum < 2:
                continue

        elementfun = lambda x: item.intfirstdiv[x.gauseq - 1]
        eq = []
        for itnl in unkitnlL:
            if type(itnl) is rxmol.Dihd:
                for func in itnl.dihdfunctions:
                    if func.forceconst != 'XXXXXX':
                        continue
                    eq.append(-func.periodicity * np.sin(
                        (func.periodicity * itnl.anglevalue - func.phase) *
                        np.pi / 180) * elementfun(itnl))
            else:
                eq.append(0.0)
                eq.append(2 * elementfun(itnl))
        E.append(eq)
    return E


def eqhessianL(unkitnlL):

    E = []
    for item in unkitnlL:
        firstdiv = lambda x: item.intfirstdiv[x.gauseq - 1]
        elementfun = lambda x: item.intseconddiv[utils.find_num_of_LTri(x.gauseq, x.gauseq)]
        eq = []
        for itnl in unkitnlL:
            if type(itnl) is rxmol.Dihd:
                for func in itnl.dihdfunctions:
                    if func.forceconst != 'XXXXXX':
                        continue
                    eq.append(-func.periodicity**2 * np.cos(
                        (func.periodicity * itnl.anglevalue - func.phase) * np.
                        pi / 180) * firstdiv(itnl)**2 - func.periodicity * np.
                              sin((func.periodicity * itnl.anglevalue - func.
                                   phase) * np.pi / 180) * elementfun(itnl))
            else:
                eq.append(2 * firstdiv(itnl)**2)
                eq.append(2 * elementfun(itnl))
        E.append(eq)
    return E


def eqrightL(hprime, unkitnlL):
    qmfile = GlobalSetting.qmfile
    rightL = []
    for item in unkitnlL:
        if type(item) is rxmol.Dihd:  # check unknown numbers
            unknum = list(
                filter(lambda x: True if x.forceconst == 'XXXXXX' else False,
                       item.dihdfunctions))
            unknum = len(unknum)
            if unknum < 2:
                continue
        rightL.append(-qmfile.fchk.intforce[item.gauseq - 1] +
                      hprime.fchk.intforce[item.gauseq - 1])
    for item in unkitnlL:
        hessnum = utils.find_num_of_LTri(item.gauseq, item.gauseq)
        rightL.append(qmfile.fchk.inthessian[hessnum] - hprime.fchk.inthessian[
            hessnum])

    return np.array(rightL)


def solve(leftL, rightL, unkitnlL):
    result = np.linalg.solve(leftL, rightL)
    res = result
    result = iter(result)
    for item in unkitnlL:
        if type(item) is rxmol.Dihd:
            for func in item.dihdfunctions:
                if func.forceconst == 'XXXXXX':
                    func.forceconst = next(result) * 627.5095
        elif type(item) is rxmol.Angle:
            item.forceconst = next(result) * 627.5095
            item.eqvalue = item.anglevalue - next(result) * 180.0 / np.pi
        elif type(item) is rxmol.Bond:
            item.forceconst = next(result) * 2240.899
            item.eqvalue = item.length - next(result) / 0.5291772086
        elif type(item) is rxmol.Improper:
            raise
    return res


def summarize(finalfuncL, itnlcordL, originalname, finalhead, method):
    # Summarize
    mmcom = GlobalSetting.mmfile.com
    logging.info('Start Summarizing')
    for func in finalfuncL:
        if func.type == 'dihd':
            if func.known is False:
                parm = [0.0,0.0,0.0,0.0]
                phase = [0,0,0,0]
                i = 0
                for item in itnlcordL:
                    if item.func == func:
                        fcs = [x.forceconst for x in item.dihdfunctions]
                        parm = [x+y for x,y in zip(parm,fcs)]
                        i+=1
                func.forceconst = [func.npaths*x/i for x in parm]
                print('haha,', func.forceconst)
                for i, item in enumerate(func.forceconst):
                    if item<0:
                        phase[i]=180
                        func.forceconst[i]*=-1
                func.phase = phase

        elif func.type == 'improper':
                res = 0
                i = 0
                for item in itnlcordL:
                    if item.func == func:
                        res += item.forceconst
                        i += 1
                func.forceconst = res / i
        else:
            if func.forceconst == 'XXXXXX':
                res = 0
                eqres = 0
                i = 0
                for item in itnlcordL:
                    if item.func == func:
                        res += item.forceconst
                        eqres += item.eqvalue
                        i += 1
                func.forceconst = res / i
                func.eqvalue = eqres / i

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
            parm = dihd.forceconst
            phase = dihd.phase

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
        tailstring += 'HrmBnd1  ' + ' '.join(
            [x.center(3, ' ')
             for x in angle.repr.split()]) + '  ' + '{:>.10f}'.format(
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
             for x in bond.repr.split()]) + '  ' + '{:>.10f}'.format(
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
    return finalname


def unavgsummarize(finalfuncL, itnlcordL, originalname, hprimehead, hprimetail,
                   method):
    # Summarize
    mmcom = GlobalSetting.mmfile.com
    logging.info('Start Summarizing')

    # Build tailstring
    tailstring = ''
    for item in itnlcordL:
        if type(item) is rxmol.Dihd:
            phase = [x.phase for x in item.dihdfunctions]
            parm = [item.npaths * x.forceconst for x in item.dihdfunctions]
            tailstring += 'AmbTrs  ' + ' '.join(
                [x.center(3, ' ')
                 for x in item.repr.split()]) + '  ' + ' '.join(
                     [str(x).center(3, ' ') for x in phase]) + '  ' + ' '.join(
                         ['{:>.10f}'.format(x)
                          for x in parm]) + '   ' + str(item.npaths) + '\n'
        elif type(item) is rxmol.Angle:
            if item.forceconst == MMFunction.unknownsign:
                parm = '0.000'
                logging.critical('Force constant is not determined for angle ' +
                                 item.repr)
                raise
            else:
                parm = item.forceconst
                tailstring += 'HrmBnd1  ' + ' '.join(
                    [x.center(3, ' ')
                     for x in item.repr.split()]) + '  ' + '{:>.10f}'.format(
                             parm) + '  ' + '{:>9.5f}'.format(item.eqvalue) + '\n'
        elif type(item) is rxmol.Bond:
            if item.forceconst == MMFunction.unknownsign:
                parm = '0.000'
                logging.critical('Force constant is not determined for bond ' +
                                 item.repr)
                raise
            else:
                parm = item.forceconst
            tailstring += 'HrmStr1  ' + ' '.join(
                [x.center(3, ' ')
                 for x in item.repr.split()]) + '  ' + '{:>.10f}'.format(
                     parm) + '  ' + '{:>7.5f}'.format(item.eqvalue) + '\n'
        elif type(item) is rxmol.Improper:
            if item.forceconst == MMFunction.unknownsign:
                logging.critical('Force constant is not determined for improper ' +
                                 item.repr)
                raise
            else:
                parm = item.forceconst
                if parm < 0:
                    item.phase = 180

                    parm = -parm
            tailstring += 'ImpTrs  ' + ' '.join([
                x.center(3, ' ') for x in item.repr.split()
            ]) + '  ' + '{:>.10f}'.format(parm) + '  ' + '{:6.2f}'.format(
                item.phase) + '  ' + str(item.
                                         npaths) + '\n'

    tailstring += hprimetail + '\n\n'
    finalname = method + '_result_' + originalname 
    logging.info('Write result to file ' + finalname)
    logging.info(tailstring)
    with open(finalname, 'w') as f:
        f.write(hprimehead + tailstring)
    shutil.copy(finalname, os.path.join('..', finalname))
    return finalname


if __name__ == '__main__':
    # read inputs, set logging module
    initset()

    # read info, identify itnls, funcs
    finalfuncL, itnlcordL, unkitnlL, knownitnlL = readgeom()

    # prepare head for hess ,hprime and finalresult
    heads, vdwtails = prepHeadAndVDW()

    # prepare hess, inithprime:
    hessL = makehess(heads[0], vdwtails[0], unkitnlL, itnlcordL)
    hprime = makehprime(heads[1], vdwtails[1], itnlcordL)

    # run hess
    runhess(hessL, para=4)

    if GlobalSetting.nocalc is False:
        hprime.run()
    hprime.fchk.read()

    thismole = GlobalSetting.mole
    np.set_printoptions(linewidth=200)
    findmatch(itnlcordL)
    leftL = eqforceL(unkitnlL)
    leftL.extend(eqhessianL(unkitnlL))
    leftL = np.array(leftL)
    rightL = eqrightL(hprime, unkitnlL)
    res = solve(leftL, rightL, unkitnlL)

    summarize(finalfuncL, itnlcordL, GlobalSetting.mmfile.comname, heads[2], 'ifhf')
    unavgsummarize(finalfuncL, itnlcordL, GlobalSetting.mmfile.comname, heads[1], vdwtails[1], 'unavgifhf')


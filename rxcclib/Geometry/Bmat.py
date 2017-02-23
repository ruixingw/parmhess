#!/usr/bin/env python3
# Construct B-matrix converting Cartesian Coordinates to Internal Coordinates
import os
import shutil
import logging
import argparse
import itertools
import copy
import numpy as np
import sympy as sp
import sympy.vector as spvec
from rxcclib.File.GauAmberCOM import GauAmberCOM
from rxcclib.File.GauAmberCOM import MMFunction
import rxcclib.File.chemfiles as rxfile
import rxcclib.Geometry.molecules as rxmol
import rxcclib.utils.cclibutils as utils


class GlobalSetting(object):
    pass


class BmatrixError(Exception):
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
            for func in item.dihdfunctions:
                item.known = True
                item.dihdunkterms = []
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
                        c = tu[2].atomtype == item.c or item.c == '*'
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
    for item in itnlcordL:
        for func in finalfuncL:
            if matchitnlwithfinalfunc(item, func):
                item.func = func
                if type(item) is rxmol.Dihd:
                    item.dihdfunctions = copy.deepcopy(func.dihdfunctions)
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
        except AttributeError as e:
            logging.critical(e)
            logging.critical('Matching finalfunc was unsuccessful {} '.format(
                item.repr))

    return finalfuncL, itnlcordL, unkitnlL


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


class Bmatrix(object):
    def __init__(self, itnlcordL):

        # self.Bmatrix = np.matrix([[r1], [r2], [r3], ..., [rn]  ])
        self.mymolecule = itnlcordL[0][1].mymolecule
        self.itnlcordL = itnlcordL
        deltax = []
        for atom in self.mymolecule:
            deltax.append(np.array(atom.coords))
        self._deltax = np.array(deltax)
        self.krondelta = np.identity(self.mymolecule.natoms + 1, dtype=int)
        pass

    def itnlsympy(self):
        sp.init_printing(use_unicode=True)
        self.symBmat = []
        natoms = self.mymolecule.natoms
        string = ''
        for i in range(1, natoms + 1):
            for j in ('x', 'y', 'z'):
                string += 'X' + str(i) + j + ' '
        coordsymbols = sp.symbols(string)

        # Symbol <--> Number dict
        self.symboltonum = {}
        for x, y in zip(self.mymolecule.coordslist, coordsymbols):
            self.symboltonum[y] = x

        # connect atom coords and symbolcoords
        for i in range(1, natoms + 1):
            atom = self.mymolecule[i]
            atom.symbolcoords = coordsymbols[3 * i - 3:3 * i]

        TS = spvec.CoordSysCartesian('TS')
        i, j, k = (TS.i, TS.j, TS.k)

        for item in self.itnlcordL:
            Bt = []
            ty = type(item)
            if ty is rxmol.Dihd or ty is rxmol.Improper:
                x, y, z = item[1].symbolcoords
                avec = x * i + y * j + z * k
                x, y, z = item[2].symbolcoords
                bvec = x * i + y * j + z * k
                x, y, z = item[3].symbolcoords
                cvec = x * i + y * j + z * k
                x, y, z = item[4].symbolcoords
                dvec = x * i + y * j + z * k

                v12 = avec - bvec
                v23 = bvec - cvec
                v34 = cvec - dvec

                v12u = v12.normalize()
                v23u = v23.normalize()
                v34u = v34.normalize()

                n1 = v12u.cross(v23u).normalize()
                n2 = v23u.cross(v34u).normalize()

                m1 = n1.cross(v23u)

                x = n1.dot(n2)
                y = m1.dot(n2)

                n1pp = sp.diff(n1, item[1].symbolcoords[0])
                print("n1ppsp: ", n1pp.evalf(subs=self.symboltonum))
                # print('prepare dihd:')
                # dihd = sp.atan2(y, x)
                # print(dihd)
                # print('diff x:')
                # diff1 = sp.diff(dihd, item[1].symbolcoords[0])
                # print(diff1)
                # print(diff1.evalf(subs=self.symboltonum))
                # print('diff y:')
                # diff2 = sp.diff(dihd, item[1].symbolcoords[1])
                # print('diff z:')
                # diff3 = sp.diff(dihd, item[1].symbolcoords[2])
                # print('prep res')
                #res = sp.simplify(diff1 * i + diff2 * j + diff3 * k)
                #print(res.evalf(subs=self.symboltonum))

                with open('wocao.txt', 'w') as f:
                    f.write(str(n1pp))

                quit()
            elif ty is rxmol.Angle:
                continue
                x, y, z = item[1].symbolcoords
                avec = x * i + y * j + z * k
                x, y, z = item[2].symbolcoords
                bvec = x * i + y * j + z * k
                x, y, z = item[3].symbolcoords
                cvec = x * i + y * j + z * k

                v12 = bvec - avec
                v23 = cvec - bvec

                v12u = v12.normalize()
                v23u = v23.normalize()

                angle = sp.pi - sp.acos(v12u.dot(v23u))
                print(angle.evalf(subs=symboltonum) * 180.0 / np.pi)
                print(item.repr)
                for atom in self.mymolecule:
                    if atom not in (item[1], item[2], item[3]):
                        Bt.append([0.0, 0.0, 0.0])
                    else:
                        if atom is item[1]:
                            tmp1 = sp.diff(angle, item[1].symbolcoords[0])
                            tmp2 = sp.diff(angle, item[1].symbolcoords[1])
                            tmp3 = sp.diff(angle, item[1].symbolcoords[2])
                            Bt.append((tmp1 * i + tmp2 * j + tmp3 * k).evalf(
                                subs=symboltonum))
                        elif atom is item[2]:
                            tmp1 = sp.diff(angle, item[2].symbolcoords[0])
                            tmp2 = sp.diff(angle, item[2].symbolcoords[1])
                            tmp3 = sp.diff(angle, item[2].symbolcoords[2])
                            Bt.append((tmp1 * i + tmp2 * j + tmp3 * k).evalf(
                                subs=symboltonum))
                        elif atom is item[3]:
                            tmp1 = sp.diff(angle, item[3].symbolcoords[0])
                            tmp2 = sp.diff(angle, item[3].symbolcoords[1])
                            tmp3 = sp.diff(angle, item[3].symbolcoords[2])
                            Bt.append((tmp1 * i + tmp2 * j + tmp3 * k).evalf(
                                subs=symboltonum))
                self.symBmat.append(Bt)

            elif ty is rxmol.Bond:
                continue
                x, y, z = item[1].symbolcoords
                avec = x * i + y * j + z * k
                x, y, z = item[2].symbolcoords
                bvec = x * i + y * j + z * k

                v12 = bvec - avec
                length = v12.magnitude()

                for atom in self.mymolecule:
                    if atom not in (item[1], item[2]):
                        Bt.extend([0.0, 0.0, 0.0])
                    elif atom is item[1]:
                        tmp1 = sp.diff(length, item[1].symbolcoords[0])
                        tmp2 = sp.diff(length, item[1].symbolcoords[1])
                        tmp3 = sp.diff(length, item[1].symbolcoords[2])
                        Bt.extend([
                            x.evalf(subs=self.symboltonum)
                            for x in [tmp1, tmp2, tmp3]
                        ])
                        tmp11 = sp.diff(tmp1, item[1].symbolcoords[0])
                        print(tmp1, tmp1.evalf(subs=self.symboltonum))
                        print(tmp11, tmp11.evalf(subs=self.symboltonum))
                    elif atom is item[2]:
                        tmp1 = sp.diff(length, item[2].symbolcoords[0])
                        tmp2 = sp.diff(length, item[2].symbolcoords[1])
                        tmp3 = sp.diff(length, item[2].symbolcoords[2])
                        Bt.extend([
                            x.evalf(subs=self.symboltonum)
                            for x in [tmp1, tmp2, tmp3]
                        ])
                self.symBmat.append(Bt)
            else:
                # should never happen
                raise BmatrixError('Unexpected itnl type: ', ty)
                pass
        self.symBmat = np.array(self.symBmat)
        pass

    def firstdiv(self):
        Bmat = []
        for itnl in self.itnlcordL:
            Bt = []
            if type(itnl) is rxmol.Bond:
                this = itnl
                thisatom = [this[1], this[2]]
                v21 = this[1].coords - this[2].coords
                e21 = v21 / np.linalg.norm(v21)

                for atom in self.mymolecule:
                    if atom not in thisatom:
                        Bt.extend([0.0, 0.0, 0.0])
                    else:
                        if atom is this[1]:
                            Bt.extend(e21)
                        elif atom is this[2]:
                            Bt.extend(-e21)

                Bmat.append(Bt)
            elif type(itnl) is rxmol.Angle:
                this = itnl
                thisatom = [*this]
                v31 = this[1].coords - this[2].coords
                v32 = this[3].coords - this[2].coords
                e31 = v31 / np.linalg.norm(v31)
                e32 = v32 / np.linalg.norm(v32)

                s1 = (this.anglecos * e31 - e32) / (np.linalg.norm(v31) *
                                                    this.anglesin)
                s2 = (this.anglecos * e32 - e31) / (np.linalg.norm(v32) *
                                                    this.anglesin)
                for atom in self.mymolecule:
                    if atom not in thisatom:
                        Bt.extend([0.0, 0.0, 0.0])
                    else:
                        if atom is this[1]:
                            s = s1
                        elif atom is this[3]:
                            s = s2
                        elif atom is this[2]:
                            s = -s1 - s2

                        Bt.extend(s)
                Bmat.append(Bt)
            elif type(itnl) is rxmol.Dihd or type(itnl) is rxmol.Improper:
                this = itnl
                thisatom = [*this]
                phi1 = self.mymolecule.angle(
                    *[x.atomnum for x in thisatom[:-1]]).anglevalue
                phi2 = self.mymolecule.angle(
                    *[x.atomnum for x in thisatom[1:]]).anglevalue

                e12 = this[2].coords - this[1].coords
                e23 = this[3].coords - this[2].coords
                e32 = -e23
                e34 = this[4].coords - this[3].coords
                e43 = -e34

                n1 = np.cross(e12, e23)
                n2 = np.cross(e23, e34)
                n1 = n1/np.linalg.norm(n1)
                n2 = n2/np.linalg.norm(n2)

                l12 = np.linalg.norm(e12)
                l23 = np.linalg.norm(e23)
                l34 = np.linalg.norm(e43)

                e12u = e12 / np.linalg.norm(e12)
                e23u = e23 / np.linalg.norm(e23)
                e32u = -e23u
                e34u = e34 / np.linalg.norm(e34)
                e43u = -e34u

                n1223 = np.cross(e12u, e23u)
                n1223 = n1223 / np.linalg.norm(n1223)
                n4332 = np.cross(e43u, e32u)
                n4332 = n4332 / np.linalg.norm(n4332)

                sin123 = np.sin(phi1 * np.pi / 180)
                cos123 = np.cos(phi1 * np.pi / 180)
                sin234 = np.sin(phi2 * np.pi / 180)
                cos234 = np.cos(phi2 * np.pi / 180)

                for atom in self.mymolecule:
                    if atom not in thisatom:
                        Bt.extend([0.0, 0.0, 0.0])
                    else:
                        if atom is this[1]:
                            s1 = -n1223 / (l12 * sin123)
                            #print(s1)
                            one = np.identity(3, dtype=int)
                            minusone = -one

                            n1pp = np.cross(minusone, e23)
                            print("n1pp:", n1pp)
                            cos = this.anglecos*(
                                np.dot(
                                    np.cross(n1pp, e23u), n2)
                            )
                            sin = this.anglesin*(
                                np.dot(n1pp, n2)
                            )
                            s1 = cos - sin
                            #print(s1)
                            s = s1
                        elif atom is this[2]:
                            s2 = (n1223 / (l12 * sin123) - cos123 * n1223 /
                                  (l23 * sin123) + cos234 * n4332 /
                                  (l23 * sin234))
                            s = s2
                        elif atom is this[3]:
                            s3 = (n4332 / (l34 * sin234) - cos234 * n4332 /
                                  (l23 * sin234) + cos123 * n1223 /
                                  (l23 * sin123))
                            s = s3
                        elif atom is this[4]:
                            s4 = -n4332 / (l34 * sin234)
                            s = s4
                        else:
                            # Should never happen
                            print(atom.name, atom.atomnum)
                            raise BmatrixError(
                                'Unexpected atom in dihd of Bmatrix.construct')
                        Bt.extend(s)
                Bmat.append(Bt)
            else:
                # Should never happen
                print(itnl)
                raise BmatrixError("Unexpected itnl type: " + str(type(itnl)))

        self.Bmat = np.array(Bmat)
        return self.Bmat

    def renew(self):
        del self.Bmat
        return self.construct()


if __name__ == '__main__':
    initset()
    finalfuncL, itnlcordL, unkitnlL = readgeom()

    a = Bmatrix(itnlcordL)
    a.firstdiv()
    np.set_printoptions(edgeitems=3, linewidth=45)
    #a.itnlsympy()
    res = 0
    #for i in range(0,12):
    #    res += np.dot(a.Bmat[0][i], a._deltax[i])

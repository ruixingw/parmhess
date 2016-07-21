# Molecules, atoms definition
from __future__ import print_function
import inspect
import logging
from io import StringIO
import numpy as np
import rxcclib.cclibutils as cclibutils


class rxMolError(Exception):
    def __init__(self, value):
        self.value = value

    def __repr__(self):
        return repr(self.value)

    __str__ = __repr__


class DihdForceConst(object):
    def __init__(self, value, dihd):
        self.forceconst = value
        self.dihd = dihd
        self.repr = dihd.repr

    def __repr__(self):
        return repr(self.forceconst)

    __str__ = __repr__

    def _call__(self, value):
        self.value = value


class Molecule(object):
    periotable = cclibutils.PeriodicTable()

    def __init__(self, moleculename):
        if not isinstance(moleculename, str):
            logging.error('Error: Molecule name must be str')
            raise rxMolError('Error: Molecule name must be str')

        self.name = moleculename
        self._atomlist = [None]
        self._bondlist = {}
        self._anglelist = {}
        self._dihdlist = {}
        self._improperlist = {}
        self._elementlist = {}
        for element in Molecule.periotable.element[1:]:
            self._elementlist.update({element: [None]})

    @property
    def natoms(self):
        return len(self._atomlist) - 1

    @property
    def bondlist(self):
        return self._bondlist

    @property
    def anglelist(self):
        return self._anglelist

    @property
    def dihdlist(self):
        return self._dihdlist

    @property
    def improperlist(self):
        return self._improperlist

    # Add atom and internal coordinates

    def addatom(self, idorsym, coords, unit='angstrom'):
        tmpatom = Atom(self, idorsym, coords, unit)
        self._atomlist.append(tmpatom)

    def addbond(self, atomnum1, atomnum2):
        if atomnum1 > atomnum2:
            atomnum1, atomnum2 = atomnum2, atomnum1
        self._bondlist.update({str(atomnum1) + '-' + str(atomnum2):
                               Bond(self, atomnum1, atomnum2)})

    def addangle(self, atomnum1, atomnum2, atomnum3):
        if atomnum1 > atomnum3:
            atomnum1, atomnum3 = atomnum3, atomnum1
        self._anglelist.update(
            {str(atomnum1) + '-' + str(atomnum2) + '-' + str(atomnum3):
             Angle(self, atomnum1, atomnum2, atomnum3)})

    def adddihd(self, atomnum1, atomnum2, atomnum3, atomnum4):
        if atomnum2 > atomnum3:
            atomnum2, atomnum3 = atomnum3, atomnum2
            atomnum1, atomnum4 = atomnum4, atomnum1
        elif atomnum2 == atomnum3:
            if atomnum1 > atomnum4:
                atomnum1, atomnum4 = atomnum4, atomnum1
        self._dihdlist.update({
            str(atomnum1) + '-' + str(atomnum2) + '-' + str(atomnum3) + '-' +
            str(atomnum4): Dihd(self, atomnum1, atomnum2, atomnum3, atomnum4)
        })

    def addimproper(self, atomnum1, atomnum2, atomnum3, atomnum4):
        self._improperlist.update({(str(atomnum1) + '-' + str(atomnum2) + '-' +
                                    str(atomnum3) + '-' + str(atomnum4)):
                                   Improper(self, atomnum1, atomnum2,
                                            atomnum3, atomnum4)})

    # Get atom and internal coordiantes
    def atom(self, atomnum):
        return self._atomlist[atomnum]

    def bond(self, atomnum1, atomnum2):
        if atomnum1 > atomnum2:
            atomnum1, atomnum2 = atomnum2, atomnum1
        return self._bondlist[str(atomnum1) + '-' + str(atomnum2)]

    def angle(self, atomnum1, atomnum2, atomnum3):
        if atomnum1 > atomnum3:
            atomnum1, atomnum3 = atomnum3, atomnum1
        return self._anglelist[str(atomnum1) + '-' + str(atomnum2) + '-' +
                               str(atomnum3)]

    def dihd(self, atomnum1, atomnum2, atomnum3, atomnum4):
        if atomnum2 > atomnum3:
            atomnum2, atomnum3 = atomnum3, atomnum2
            atomnum1, atomnum4 = atomnum4, atomnum1
        return self._dihdlist[str(atomnum1) + '-' + str(atomnum2) + '-' + str(
            atomnum3) + '-' + str(atomnum4)]

    def improper(self, atomnum1, atomnum2, atomnum3, atomnum4):
        return self._improperlist[str(atomnum1) + '-' + str(atomnum2) + '-' +
                                  str(atomnum3) + '-' + str(atomnum4)]

    # mole[5] as mole.atom(5) and slice
    def __getitem__(self, key):
        if isinstance(key, int):
            return self._atomlist[key]
        if isinstance(key, slice):
            start = key.start
            stop = key.stop
            if start is None:
                start = 1
            L = []
            for x in range(start, stop):
                L.append(self._atomlist[x])
            return L

    # iterate atoms
    def __iter__(self):
        for i in range(self.natoms):
            yield self[i + 1]

    # Read structure from xyz
    def readfromxyz(self, string):
        f = StringIO(string)
        for line in f:
            tmp = line.split()
            atom = tmp[0]
            if atom.isdigit():
                atom = int(atom)
            coords = np.array([tmp[1], tmp[2], tmp[3]], dtype=float)
            self.addatom(atom, coords, unit='angstrom')

    # Read connectivity to complete neighbor info
    def readconnectivity(self, conntystring):
        f = StringIO(conntystring)
        for line in f:
            tmp = line.split()
            if not tmp:
                logging.debug('End of connectivity, return.')
                return
            ite = iter(tmp)
            item0 = next(ite)
            if not item0.isdigit():
                logging.debug('End of connectivity, return.')
                return
            a = int(item0)
            try:
                b = int(next(ite))
            except StopIteration:
                continue
            self.addbond(a, b)
            while True:
                try:
                    b = next(ite)
                    b = next(ite)
                    b = int(b)
                    self.addbond(a, b)
                except StopIteration:
                    break
        angles = []
        dihds = []
        for atom1 in self:  # atom1: atom obj
            for atom2 in atom1.neighbor:  # atom2: obj
                if atom2 is atom1:
                    continue
                for atom3 in atom2.neighbor:  # atom3: obj
                    if (atom3 is atom2 or
                            atom3 is atom1):
                        continue
                    a = atom1.atomnum
                    b = atom2.atomnum
                    c = atom3.atomnum
                    if a > c:
                        a, c = c, a
                    angles.append(str(a) + '-' + str(b) + '-' + str(c))
                    for atom4 in atom3.neighbor:
                        if (atom4 is atom3 or atom4 is atom2 or
                                atom4 is atom1):
                            continue
                        a = atom1.atomnum
                        b = atom2.atomnum
                        c = atom3.atomnum
                        d = atom4.atomnum
                        if b > c:
                            b, c = c, b
                            a, d = d, a
                        dihds.append(str(a) + '-' + str(b) + '-' + str(c) + '-'
                                     + str(d))
        angles = list(set(angles))
        dihds = list(set(dihds))
        for item in angles:
            tmp = [int(x) for x in item.split('-')]
            self.addangle(*tmp)
        for item in dihds:
            tmp = [int(x) for x in item.split('-')]
            self.adddihd(*tmp)

    def readtypefromlist(self, L):
        if len(L) != len(self._atomlist):
            logging.error(
                "Error when reading atomtype from list:"
                " length is not consistent with natoms")
            raise rxMolError(
                "Error when reading atomtype from list:"
                " length is not consistent with natoms")
        ite = iter(L)
        next(ite)
        for atom in self:
            atom.atomtype = next(ite)

    def readchargefromlist(self, L):
        if len(L) != len(self._atomlist):
            logging.error(
                "Error when reading atomcharge from list:"
                " length is not consistent with natoms")
            raise rxMolError(
                "Error when reading atomcharge from list:"
                " length is not consistent with natoms")
        ite = iter(L)
        next(ite)
        for atom in self:
            atom.atomcharge = next(ite)


class Atom(object):

    periotable = cclibutils.PeriodicTable()

    def __init__(self, mole, idorsym,
                 coords,
                 unit='angstrom'):  # molecule object,int,[float,float,float]
        # Assertion
        callername = inspect.stack()[1][3]
        assert callername == 'addatom', ("Atom must be added"
                                         " via Molecule.addatom method")
        assert isinstance(mole, Molecule), ("First argument must be"
                                            " a molecule object!. Use "
                                            "Molecule.addatom method to avoid"
                                            " this problem.")

        assert unit != 'bohr' or unit != 'angstrom', ("Coordinate unit"
                                                      " must be bohr or"
                                                      " angstrom")

        self._mymolecule = mole
        if isinstance(idorsym, int):
            self._elementid = idorsym
            try:
                self._elementsym = Atom.periotable.element[self.elementid]
            except KeyError:
                logging.critical(
                    "Error when adding atom:"
                    " Idtosym not defined for atomic no:"
                    + str(self.elementid))
                raise rxMolError(
                    "Error when adding atom: Idtosym"
                    " not defined for atomic no:"
                    + str(self.elementid))
        elif isinstance(idorsym, str):
            self._elementsym = idorsym
            try:
                self._elementid = Atom.periotable.number[self.elementsym]
            except KeyError:
                logging.critical(
                    "Error when adding atom: Idtosym not"
                    " defined for atomic symbol:"
                    + self.elementsym)
                raise rxMolError(
                    "Error when adding atom: Idtosym not"
                    " defined for atomic symbol:"
                    + self.elementsym)
        else:
            logging.critical(
                "Error when adding atom: Expected atomic"
                " NO(int) or symbol(str) for input, received a"
                + str(type(idorsym)))
            raise rxMolError(
                "Error when adding atom: Expected atomic"
                " NO(int) or symbol(str) for input, received a"
                + str(type(idorsym)))

        if unit == 'bohr':
            self.coords = cclibutils.convertor(coords, "bohr", "Angstrom")
        elif unit == 'angstrom':
            self.coords = coords

        self.atomnum = mole.natoms + 1
        self._mymolecule._elementlist[self.elementsym].append(self)
        self.atomtype = self.name
        self._neighbor = []
        self.atomcharge = None

    @property
    def elementid(self):
        return self._elementid

    @property
    def elementsym(self):
        return self._elementsym

    @property
    def neighbor(self):
        return self._neighbor

    def addneighbor(self, atomobj):
        if isinstance(atomobj, Atom):
            self._neighbor.append(atomobj)
            self._neighbor = list(set(self._neighbor))
        else:
            logging.error(
                "Error when adding neighbor: atomobj must be an atom object.")
            raise rxMolError(
                "Error when adding neighbor: atomobj must be an atom object.")

    def delneighbor(self, atomobj):
        if isinstance(atomobj, Atom):
            self._neighbor.remove(atomobj)
        else:
            logging.error(
                "Error when deleting neighbor: "
                "atomobj must be an atom object.")
            raise rxMolError(
                "Error when deleting neighbor: "
                "atomobj must be an atom object.")

    @property
    def name(self):
        index = self.mymolecule._elementlist[self.elementsym].index(self)
        return (str(self.elementsym) + str(index))

    @property
    def mymolecule(self):
        return self._mymolecule

    def __str__(self):
        return "Atom object for atom " + self.name


class Bond(object):
    def __init__(self, mole, a, b):  # self, atomnum a, atomnum b
        if a > b:
            a, b = b, a
        self._a = mole[a]
        self._b = mole[b]
        self.vec = self._a.coords - self._b.coords
        self.repr = self._a.name + ' ' + self._b.name
        self._a.addneighbor(mole[b])
        self._b.addneighbor(mole[a])

    def __getitem__(self, value):
        if value == 1:
            return self._a
        elif value == 2:
            return self._b
        else:
            raise rxMolError("Index for bond object must be 1 or 2.")

    @property
    def length(self):
        return np.linalg.norm(self.vec)

    def __str__(self):
        return "Bond object of bond " + self._a.name + '-' + self._b.name


class Angle(object):
    def __init__(self, mole, a, b, c):
        if a > c:
            a, c = c, a
        self._a = mole[a]
        self._b = mole[b]
        self._c = mole[c]
        self._ab = mole[a].coords - mole[b].coords
        self._bc = mole[b].coords - mole[c].coords
        self.repr = self._a.name + ' ' + self._b.name + ' ' + self._c.name

    def __getitem__(self, value):
        if value == 1:
            return self._a
        elif value == 2:
            return self._b
        elif value == 3:
            return self._c
        else:
            raise rxMolError("Index for angle object must be 1, 2 or 3.")

    @property
    def anglevalue(self):
        v1u = self._ab / np.linalg.norm(self._ab)
        v2u = self._bc / np.linalg.norm(self._bc)
        angle = 180.0 - np.arccos(np.dot(v1u, v2u)) * 180.0 / np.pi
        if np.isnan(angle):
            if (v1u == v2u).all():
                return 0.0
            else:
                return 180.0
        return angle

    def __str__(self):
        return ("Angle object of angle " + self._a.name +
                '-' + self._b.name + '-' + self._c.name)


class Dihd(object):
    def __init__(self, mole, a, b, c, d):
        if b > c:
            b, c = c, b
            a, d = d, a
        self._a = mole[a]
        self._b = mole[b]
        self._c = mole[c]
        self._d = mole[d]
        self.repr = (self._a.name + ' ' + self._b.name + ' ' +
                     self._c.name + ' ' + self._d.name)

    def __getitem__(self, value):
        if value == 1:
            return self._a
        elif value == 2:
            return self._b
        elif value == 3:
            return self._c
        elif value == 4:
            return self._d
        else:
            raise rxMolError(
                "Index for dihedral object must be 1, 2 ,3 or 4.")

    @property
    def dihdvalue(self):
        v1 = self._a.coords - self._b.coords
        v2 = self._b.coords - self._c.coords
        v3 = self._c.coords - self._d.coords
        v1u = v1 / np.linalg.norm(v1)
        v2u = v1 / np.linalg.norm(v2)
        v3u = v1 / np.linalg.norm(v3)
        n1 = np.cross(v1u, v2u)
        n2 = np.cross(v2u, v3u)
        dihd = np.arccos(np.dot(n1, n2)) * 180.0 / np.pi
        if np.isnan(dihd):
            if (n1 == n2).all():
                return 0.0
            else:
                return 180.0
        return dihd

    def __str__(self):
        return ("Dihedral object of dihedral " + self._a.name +
                '-' + self._b.name + '-' + self._c.name + '-' + self._d.name)


class Improper(object):
    def __init__(self, mole, a, b, c, d):
        self._a = mole[a]
        self._b = mole[b]
        self._c = mole[c]
        self._d = mole[d]
        self.repr = (self._a.name + ' ' + self._b.name +
                     ' ' + self._c.name + ' ' + self._d.name)

    def __getitem__(self, value):
        if value == 1:
            return self._a
        elif value == 2:
            return self._b
        elif value == 3:
            return self._c
        elif value == 4:
            return self._d
        else:
            raise rxMolError(
                "Index for improper object must be 1, 2 ,3 or 4.")

    @property
    def impropervalue(self):
        v1 = self._a.coords - self._b.coords
        v2 = self._b.coords - self._c.coords
        v3 = self._c.coords - self._d.coords
        v1u = v1 / np.linalg.norm(v1)
        v2u = v1 / np.linalg.norm(v2)
        v3u = v1 / np.linalg.norm(v3)
        n1 = np.cross(v1u, v2u)
        n2 = np.cross(v2u, v3u)
        improper = np.arccos(np.dot(n1, n2)) * 180.0 / np.pi
        if np.isnan(improper):
            if (n1 == n2).all():
                return 0.0
            else:
                return 180.0
        return improper

    def __str__(self):
        return ("Improper object of improper " + self._a.name +
                '-' + self._b.name + '-' + self._c.name + '-' + self._d.name)


if __name__ == '__main__':
    import doctest
    doctest.testmod()

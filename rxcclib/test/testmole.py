#!/usr/bin/env python3
import rxcclib.molecules as rxmol
import unittest, os


class TestMole(unittest.TestCase):
    def setUp(self):
        self.mole = rxmol.Molecule("H2O2")
        self.assertTrue(isinstance(self.mole, rxmol.Molecule))

        self.mole.addatom('O', np.array([0.0, 0.0, 0.0]), unit='angstrom')
        self.mole.addatom('H', np.array([0.0, 1.0, 0.0]), unit='angstrom')
        self.mole.addatom('O', np.array([0.0, 0.0, 1.0]))
        self.mole.addatom('H', np.array([1.0, 0.0, 0.0]), unit='angstrom')
        self.assertTrue(isinstance(self.mole[1], rxmol.Atom))

        self.mole.addbond(1, 2)
        self.mole.addbond(1, 3)
        self.mole.addbond(3, 4)
        self.assertTrue(isinstance(self.mole.bond(1, 3), rxmol.Bond))

        self.mole.addangle(2, 1, 3)
        self.mole.addangle(4, 3, 1)
        self.assertTrue(isinstance(self.mole.angle(2, 1, 3), rxmol.Angle))

        self.mole.adddihd(2, 1, 3, 4)
        self.assertTrue(isinstance(self.mole.dihd(4, 3, 1, 2), rxmol.Dihd))

        # @property natoms
        self.assertEqual(self.mole.natoms, 4)

        self.assertEqual(self.mole.atom(1).elementid, 8)
        # def __iter__
        ite = iter(self.mole)
        for atom in self.mole:
            print(atom.name)
        # def __next__
        self.assertEqual(next(ite).name, 'O1')
        self.assertEqual(next(ite).name, 'H1')
        ite2 = iter(self.mole)
        self.assertEqual(next(ite2).name, 'O1')
        self.assertEqual(next(ite2).name, 'H1')
        self.assertEqual(next(ite).name, 'O2')
        self.assertEqual(next(ite).name, 'H2')
        self.assertEqual(next(ite2).name, 'O2')
        self.assertEqual(next(ite2).name, 'H2')

        # def __getitem__
        self.assertEqual(str(self.mole[1].name), 'O1')

    def test_mole(self):
        # @property xyz
        np.testing.assert_array_equal(self.mole[1].coords,
                                      np.array([0.0, 0.0, 0.0]))
        np.testing.assert_array_equal(self.mole[2].coords,
                                      np.array([0.0, 1.0, 0.0]))
        np.testing.assert_array_equal(self.mole[3].coords,
                                      np.array([0.0, 0.0, 1.0]))
        np.testing.assert_array_equal(self.mole[4].coords,
                                      np.array([1.0, 0.0, 0.0]))
        # @property atomnum
        self.assertEqual(self.mole[1].atomnum, 1)
        # @property elementid
        self.assertEqual(self.mole[1].elementid, 8)
        # @property elementsym
        self.assertEqual(self.mole[1].elementsym, 'O')
        # @property name
        self.assertEqual(self.mole[1].name, 'O1')
        self.assertEqual(self.mole[1].atomtype, 'O1')
        # @property mymolecule
        self.assertEqual(str(self.mole.atom(1).mymolecule.name), 'H2O2')
        # test_bond(self):
        # Bondlength
        self.assertEqual(self.mole.bond(1, 3).length, 1.0)
        # test_angle(self):
        self.assertEqual(self.mole.angle(2, 1, 3).anglevalue, 90.0)
        # test dihd
        self.assertEqual(self.mole.dihd(4, 3, 1, 2).dihdvalue, 90.0)
        # test_bondfunc(self):
        self.mole[1].atomtype = 'oh'
        self.mole[2].atomtype = 'ho'
        self.mole[3].atomtype = 'oh'
        self.mole[4].atomtype = 'ho'
    # test neighbor

    def test_readfile(self):
        benz = rxmol.Molecule("Benzene")
        with open('samples/ben.xyz', 'r') as f:
            f = f.read()
            benz.readfromxyz(f)
        self.assertEqual(benz[1].name, 'C1')
        self.assertEqual(benz[12].name, 'H6')
        #test_connty
        with open('samples/cnnty.com', 'r') as f:
            cn = f.read()
        benz.readconnectivity(cn)
        self.assertEqual(len(benz.bondlist.values()), 12)
        self.assertEqual(len(benz.anglelist.values()), 18)
        self.assertEqual(len(benz.dihdlist.values()), 24)
        self.assertAlmostEqual(
            benz.angle(1, 2, 3).anglevalue,
            120.00, delta=0.1)
        os.system(
            'rm A* q* Q* p* esout *gaussian* samples/bencom.fchk samples/bencom.chk samples/bencom.log')


if __name__ == '__main__':
    import numpy as np
    unittest.main()

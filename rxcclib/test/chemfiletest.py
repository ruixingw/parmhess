#!/usr/bin/env python3
import rxcclib.File.chemfiles as rxccfile
import rxcclib.Geometry.molecules as rxmol
import unittest, os, logging
import numpy as np
from io import StringIO

rxccfile.GauCOM.g09rt = 'g09'
rxccfile.GauCOM.g09a2rt = 'g09'
os.system(
    'rm A* q* Q* p* esout *Gaussian* samples/bencom.fchk samples/bencom.chk samples/bencom.log')


class TestFile(unittest.TestCase):
    def test_comfchk(self):
        file = rxccfile.File('samples/bencom')
        self.assertIsInstance(file, rxccfile.File)
        self.assertIsInstance(file.com, rxccfile.GauCOM)
        self.assertIsInstance(file.log, rxccfile.GauLOG)
        self.assertIsInstance(file.fchk, rxccfile.GauFCHK)
        self.assertIsInstance(file.mol2, rxccfile.Mol2)
        file.com.rung09()
        file.com.isover()
        file.runformchk()
       # self.assertEqual(file.com.read(), True)
        self.assertEqual(file.fchk.read(), True)
        self.assertEqual(file.natoms, 12)
        self.assertEqual(file.multiplicity, 1)
        self.assertEqual(file.totalcharge, 0)
        self.assertIsInstance(file.fchk.xyz, str)
        tmp = file.fchk.xyz[::-1]
        tmp = float(tmp.split()[0][::-1])
        hess = file.fchk.find33Hessian(3, 5)
        self.assertAlmostEqual(hess[0][0], -2.62909045e-2)
        self.assertAlmostEqual(hess[1][1], 3.38743754e-2)
        self.assertAlmostEqual(hess[2][2], 7.19580040e-3)

    def test_logmol2(self):
        file = rxccfile.File('samples/benresp')
        file.log.runantecham()
        file.mol2.read()
        self.assertEqual(file.atomtypelist[0], None)
        self.assertEqual(file.atomtypelist[1], 'ca')
        self.assertEqual(file.atomtypelist[12], 'ha')
        self.assertEqual(file.atomchargelist[0], None)
        self.assertEqual(file.atomchargelist[1], -0.117738)
        self.assertEqual(file.atomchargelist[12], 0.117738)
    # def test_MMcom(self):
    #     mmfile=rxccfile.File('samples/mmfile')
    #     mmfile.com.read()
    #     xyz=StringIO(mmfile.com.xyz)

    #        self.assertEqual()


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    unittest.main()

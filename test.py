import rxcclib.molecules as rxmol
import parmhessdev as phd
import pdb

a = phd.HessFile('hess0',1)
mole = rxmol.Molecule('123')
a.run()
a.com.read()
mole.readfromxyz(a.fchk.xyz)
mole.readconnectivity(a.com.connectivity)

dihd = mole.dihdlist['1-2-3-4']
a.itnl = dihd
dihd.sinX = 1
dihd.cosX = 0 
a.getcoorddivs()


#!/usr/bin/env python3
#################################################################
# Tsubasa -- Get ready for Hessian Fitting Parameterization
# Automated Optimization, Freq&Hessian, RESP Charge Calculation
# First version by Ruixing at 2-Feb-2016
#################################################################
from __future__ import print_function
import rxcclib.molecules as rxmol
import rxcclib.chemfiles as rxccfile
import os
import sys
import yaml
import argparse
import logging
import shutil

# Parse input
parser = argparse.ArgumentParser()
parser.add_argument('-i',
                    dest='inputgeom',
                    default=False,
                    help=('Inputfile including mo'
                          'lecular specs and connectivity'))
parser.add_argument('-c',
                    dest='configfile',
                    default=False,
                    help='Tsubasa config file.')
parser.add_argument('--readvdw',
                    dest='externalvdwfile',
                    default=False,
                    help='If provided, read external vdW parameters from file',
                    nargs=1)
parser.add_argument(
    '--startfrom',
    default='opt',
    choices=['freq', 'resp', 'antechamber', 'readmol2'],
    help=("Start from a certain step. Choices"
          "=['opt','freq','resp','antechamber','readmol2']"))
parser.add_argument(
    '--stopafter',
    default='readmol2',
    choices=['opt', 'freq', 'resp', 'antechamber'],
    help="Stop after a certain step. Choices=['freq','resp','antechamber']")
args = parser.parse_args()
inputgeom = args.inputgeom
ymlfile = args.configfile
externalvdw = args.externalvdwfile
print(externalvdw)
if externalvdw:
    externalvdw = externalvdw[0]

startfrom = args.startfrom
if startfrom == 'opt':
    startfrom = 0
elif startfrom == 'freq':
    startfrom = 1
    logging.warning('Start from freq')
elif startfrom == 'resp':
    startfrom = 2
    logging.warning('Start from resp')
elif startfrom == 'antechamber':
    startfrom = 3
    logging.warning('Start from antechamber')
elif startfrom == 'readmol2':
    startfrom = 4
    logging.warning('Start from readmol2')

stopafter = args.stopafter
if stopafter == 'opt':
    stopafter = 1
    logging.warning('stop after opt')
elif stopafter == 'freq':
    stopafter = 2
    logging.warning('stop after freq')
elif stopafter == 'resp':
    stopafter = 3
    logging.warning('stop after resp')
elif stopafter == 'antechamber':
    stopafter = 4
    logging.warning('stop after antechamber')
elif stopafter == 'readmol2':
    stopafter = 5

##############
# copy config file to current pathp if no argument is specified
pwd = os.path.split(os.path.realpath(__file__))[0]
if inputgeom is False and ymlfile is False and externalvdw is False:
    gauname = False
    ymlname = False
    vdwname = False
    for filename in os.listdir():
        if filename.find('.gau') >= 0:
            gauname = filename[0:filename.find('.gau')]
        if filename.find('.yml') >= 0:
            ymlname = filename[0:filename.find('.yml')]
        if filename.find('vdw') >= 0:
            externalvdw = filename
    if gauname and ymlname and ymlname != gauname:
        print(ymlname, gauname)
        logging.critical('Inconsistent name of inputgeom and config file. ')
        sys.exit()
    elif gauname and not ymlname:
        shutil.copyfile(
            os.path.join(pwd, 'config.yml'), os.path.join(os.getcwd(),
                                                          gauname + '.yml'))
        logging.warning(
            'Config file is not found. A template'
            ' is copied to current directory. Program will now quit.')
        sys.exit()
    elif gauname and ymlname and ymlname == gauname:
        inputgeom = gauname
        ymlfile = ymlname
    else:
        logging.critical('No inputgeom file found.')
        sys.exit()
###########################

# Logging module setting. Print INFO on screen and DEBUG INFO in file##
logging.basicConfig(filename=inputgeom + '.tsubasa',
                    level=logging.DEBUG,
                    filemode='w')
console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)-8s %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)
##########################################################################

logging.info('Read config from ' + ymlfile + '.yml\n')
if not os.path.isfile(ymlfile + '.yml'):
    logging.critical(ymlname + '.yml does not exist! ')
    sys.exit()

with open(ymlfile + '.yml', 'r') as f:
    yml = f.read()

config = yaml.load(yml)
rxccfile.GauCOM.g09rt = config['g09rt']
rxccfile.GauCOM.g09a2rt = config['g09a2rt']
rxccfile.GauLOG.antecommand = config['antechamber']
clean = config['clean']
opthead = config['opthead'] + '\n'
opttail = config['opttail'] + '\n'
freqhead = config['freqhead'] + '\n'
resphead = config['resphead'] + '\n'
resptail = config['resptail'] + '\n'
mmhead = config['mmhead'] + '\n'

optfile = rxccfile.File('opt' + inputgeom)
freqfile = rxccfile.File('freq' + inputgeom)
respfile = rxccfile.File('resp' + inputgeom)

# Run Calculations
if startfrom < 1:
    with open(inputgeom + '.gau', 'r') as initxyz:
        with open(optfile.comname, 'w') as f:
            f.write(opthead)
            f.write(initxyz.read())
            f.write(opttail)
    logging.info('Runing optimization...')
    optfile.com.rung09()
    try:
        optfile.com.isover()
    except:
        logging.critical('Optimization failed')
        sys.exit()

    if stopafter == 1:
        logging.warning('User request stopping after optimization')
        sys.exit()

if startfrom < 2:
    with open(freqfile.comname, 'w') as f:
        freqhead = ('%chk=' + os.path.split(optfile.chkname)[1]
                    + '\n' + freqhead)
        f.write(freqhead)
    logging.info('Running frequency calculation...')
    freqfile.com.rung09()

    try:
        freqfile.com.isover()
    except:
        logging.critical('Frequency calculation failed')
        sys.exit()

    if stopafter == 2:
        logging.warning('User request stopping after frequency calculation')
        sys.exit()
if startfrom < 3:
    with open(respfile.comname, 'w') as f:
        resphead = '%chk=' + os.path.split(freqfile.chkname)[
            1] + '\n' + resphead
        f.write(resphead)
        f.write(resptail)
    logging.info('Running RESP calculation...')
    respfile.com.rung09a2()
    try:
        respfile.com.isover()
    except:
        logging.critical('MK calculation failed.')
        sys.exit()

    if stopafter == 3:
        logging.warning('User request stop after resp')
        sys.exit()

if startfrom < 4:
    logging.info('Run antechamber:')
    respfile.log.runantecham()
    if stopafter == 4:
        logging.warning('User request stopping after antechamber')
        sys.exit()

logging.info('Format CHK file by: ')

# Read fchk : coordinates, charge, spin, natoms
qmfile = freqfile
qmfile.runformchk()
qmfile.fchk.read()

logging.debug('...done\n\nBuild mmxyz...')
thisgeom = rxmol.Molecule('this')
xyz = ['']
thisgeom.readfromxyz(qmfile.fchk.xyz)

# Read mol2 file to get Charge&Atomtype and assign to atoms
if startfrom < 5:
    respfile.mol2.read()
    thisgeom.readchargefromlist(respfile.mol2.atomchargelist)
    thisgeom.readtypefromlist(respfile.mol2.atomtypelist)
if stopafter == 3:
    logging.warning('User request stop after readmol2')
    sys.exit()

connectivity = ''
logging.debug('Read internal coordinates from connectivity...')
with open(inputgeom + '.gau', 'r') as f:
    for line in f:
        if line == '\n':
            break
    for line in f:
        if line == '\n':
            break
        connectivity += line

thisgeom.readconnectivity(connectivity)

# Build MM input
mmxyz = str(qmfile.totalcharge) + ' ' + str(qmfile.multiplicity) + '\n'


def f2s(fl):
    return "{: .12f}".format(fl)


for i in range(1, qmfile.natoms + 1):
    mmxyz = mmxyz + thisgeom[i].elementsym + '-' + thisgeom[
        i].atomtype + '-' + '{:.6f}'.format(thisgeom[
            i].atomcharge) + '   ' + '   '.join(
                [f2s(x) for x in thisgeom[i].coords]) + '\n'

mmxyz = mmxyz + '\n'


class Bondfunc(object):

    def __init__(self, mole, bondobj):
        a = bondobj[1].atomtype
        b = bondobj[2].atomtype
        if a > b:
            a, b = b, a
        self.link = a + ' ' + b


class Anglefunc(object):

    def __init__(self, molecule, angleobj):
        a = angleobj[1].atomtype
        b = angleobj[2].atomtype
        c = angleobj[3].atomtype
        if a > c:
            a, c = c, a
        self.link = a + ' ' + b + ' ' + c


class Dihdfunc(object):

    def __init__(self, molecule, dihdobj):
        a = dihdobj[1].atomtype
        b = dihdobj[2].atomtype
        c = dihdobj[3].atomtype
        d = dihdobj[4].atomtype
        self.periodicity = 2
        self.phase = 180.0
        self.npaths = 1.0
        if b > c:
            a, d = d, a
            b, c = c, b
        elif b == c:
            if a > d:
                a, d = d, a
        self.link = a + ' ' + b + ' ' + c + ' ' + d

# add funcs
thisgeom.dihdfunc = {}
thisgeom.anglefunc = {}
thisgeom.bondfunc = {}
for bond in thisgeom.bondlist.values():
    thisfunc = Bondfunc(thisgeom, bond)
    thisgeom.bondfunc.update({thisfunc.link: bond})
    bond.func = thisfunc
for angle in thisgeom.anglelist.values():
    thisfunc = Anglefunc(thisgeom, angle)
    thisgeom.anglefunc.update({thisfunc.link: angle})
    angle.func = thisfunc
for dihd in thisgeom.dihdlist.values():
    thisfunc = Dihdfunc(thisgeom, dihd)
    thisgeom.dihdfunc.update({thisfunc.link: dihd})
    dihd.func = thisfunc

# Read internal coordinates, add bond, angle, dihedral and define MM functions
sorteddihd = sorted(thisgeom.dihdfunc.keys(),
                    key=lambda item: item.split()[1] + ' ' + item.split()[2])
sortedangle = sorted(thisgeom.anglefunc.keys(),
                     key=lambda item: item.split()[1] + ' ' + item.split()[0])
sortedbond = sorted(thisgeom.bondfunc.keys(), key=lambda item: item.split()[0])

# Build input file and MMtail(functions)
mmfile = rxccfile.File('mm' + gauname)
mmtail = ''
input = 'natoms=' + str(qmfile.natoms) + '\nmmfile=' + os.path.split(
    mmfile.comname)[1] + '\n' + 'qmfchk=' + os.path.split(freqfile.fchkname)[
        1] + '\n' + 'qmlog=' + os.path.split(freqfile.logname)[1] + '\n'
input = input + '\n\nLink start\n'

# dihedral is assigned n=2 ,phase=180 and Npaths=1 temporarily

for key in sorteddihd:
    # key is sorted dihdfunc
    mmtail = mmtail + 'AmbTrs ' + key + ' 0 180 0 0 0.0 XXXXXX 0.0 0.0 1.0\n'
    # For each key of dihdfunc.key, filter x in
    # dihedral.list.values(obj) to find out whose x.func(obj) match this key
    # 'this' will be the dihd obj who satisfy the condition aforementioned
    this = filter(lambda x: x.func.link == key, thisgeom.dihdlist.values())
    this = list(set(this))
    for x in this:
        input += str(x[1].atomnum) + '-' + str(x[2].atomnum) + '-' + str(x[
            3].atomnum) + '-' + str(x[4].atomnum) + '\n'
    input += 'next  # ' + key + '\n'

for key in sortedangle:
    this = filter(lambda x: x.func.link == key, thisgeom.anglelist.values())
    this = list(set(this))
    total = 0
    for num, x in enumerate(this):
        total += x.anglevalue
        now = total / (num + 1)
        if abs(x.anglevalue - now) > 3 and total != 0:
            logging.warning('Angle ' + x.repr +
                            ' has very different angle value of  {:.4f}'.format(
                                x.anglevalue) + ' compared to ' + x.func.link +
                            ' {:.4f}'.format(now))
        input += str(x[1].atomnum) + '-' + str(x[2].atomnum) + '-' + str(x[
            3].atomnum) + '\n'
    input += 'next  # ' + key + '\n'
    total = total / len(this)
    total = "{:.4f}".format(total)
    mmtail = mmtail + 'HrmBnd1 ' + key + ' XXXXXX ' + total + '\n'
for key in sortedbond:
    this = filter(lambda x: x.func.link == key, thisgeom.bondlist.values())
    this = list(set(this))
    total = 0
    for num, x in enumerate(this):
        total += x.length
        now = total / (num + 1)
        if abs(x.length - now) > 0.1 and total != 0:
            logging.warning('bond ' + x.repr +
                            ' has very different length of {:.4f}'.format(
                                x.length) + ' compared to ' + x.func.link +
                            ' {:.5f}'.format(now))
        input += str(x[1].atomnum) + '-' + str(x[2].atomnum) + '\n'
    input += 'next  # ' + key + '\n'
    total = total / len(this)
    total = "{:.5f}".format(total)
    mmtail = mmtail + 'HrmStr1 ' + key + ' XXXXXX ' + total + '\n'

# Add Nonbon function and vdW parameters
with open('input.inp', 'w') as f:
    f.write(input)
mmtail = mmtail + 'Nonbon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 -1.2\n'

radii = {}
welldepth = {}
# read all vdw in files
with open(os.path.join(pwd, 'vdw.dat'), 'r') as f:
    for string in f:
        item = string.split()
        radii.update({item[0].strip(' '): item[1].strip(' ')})
        welldepth.update({item[0].strip(' '): item[2].strip(' ')})
if externalvdw:
    logging.info('Read user provided vdW parameters from ' + externalvdw)
    with open(externalvdw, 'r') as f:
        for string in f:
            item = string.split()
            radii.update({item[0].strip(' '): item[1].strip(' ')})
            welldepth.update({item[0].strip(' '): item[2].strip(' ')})
item = []

# find existing atomtypes
for atom in thisgeom:
    item.append(atom.atomtype)
item = list(set(item))
for i in range(0, len(item)):
    mmtail += 'VDW ' + item[i] + '  ' + radii[item[i]] + '  ' + welldepth[item[
        i]] + '\n'
mmtail = mmtail + '\n'

# write mmfile
with open(mmfile.comname, 'w') as f:
    f.write(mmhead)
    f.write(mmxyz)
    f.write(connectivity)
    f.write('\n')
    f.write(mmtail)

logging.debug('...done\n\nEND')
os.system(clean)

os.system('mkdir tsubasa')
os.system('mv * tsubasa')
os.system('cp tsubasa/mm* .')
os.system('cp tsubasa/' + os.path.split(freqfile.fchkname)[1] + ' .')
os.system('cp tsubasa/' + os.path.split(freqfile.fchkname)[1] + ' .')
os.system('cp tsubasa/' + os.path.split(freqfile.logname)[1] + ' .')
os.system('cp tsubasa/input.inp .')
os.system(clean)

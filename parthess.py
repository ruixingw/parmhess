#!/usr/bin/env python3
from __future__ import print_function
import rxcclib.molecules as rxmol
import rxcclib.chemfiles as rxfile
import numpy as np
import argparse,os,logging,itertools,copy

############Parse input
parser=argparse.ArgumentParser()
parser.add_argument('mmfile',help="mmfile prepared by tsubasa.")
parser.add_argument('qmfile',help="FCHK file from a frequency calculation.")
parser.add_argument('--delete',help="Delete all temp files. Default option will save files to ./phfdebug",action='store_true')

args=parser.parse_args()
mmfile=args.mmfile
originname=mmfile
qmfile=args.qmfile
delete=args.delete


#### Logging module setting. Print INFO on screen and DEBUG INFO in file###########
logging.basicConfig(filename=mmfile+'.phfout',level=logging.DEBUG,filemode='w')
console=logging.StreamHandler()
console.setLevel(logging.INFO)
formatter=logging.Formatter('%(levelname)-8s %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)
###################################################################################

print('\n\n                         Partial Hessian Fitting for parameterization\n\n')

logging.info("Provided mmfile input: "+mmfile+' and qmfchk input: '+qmfile)

mmfile=mmfile[0:mmfile.find('.')]
qmfile=qmfile[0:qmfile.find('.')]


##prepare del.sh, del old files, print title
with open('del.sh','w') as f:
    f.write('rm -rf *hprime* hess* tmpphffile del.sh\n')
os.system('chmod +x del.sh')


## Instantialize Molecule and file
mole=rxmol.Molecule('mmcom')
mmcom=rxfile.File(mmfile).com
qmfile=rxfile.File(qmfile)
qmfile.fchk.read()
mmcom.read()
mole.readfromxyz(mmcom.xyz)
mole.readchargefromlist(mmcom.atomchargelist)
mole.readtypefromlist(mmcom.atomtypelist)
mole.readconnectivity(mmcom.connectivity)

## Count unknown parameters
nunk=0
nozomuL=[]
nozomuL.extend(mmcom.nozomubondfunc)
nozomuL.extend(mmcom.nozomuanglefunc)
nozomuL.extend(mmcom.nozomudihdfunc)
nozomuL.extend(mmcom.nozomuimproperfunc)

for item in nozomuL:
    if item.type=='dihd':
        for paras in item.forceconst:
            if paras==rxfile.mmfunction.unknownsign:
                nunk+=1
    elif item.type=='improper':
        ##### match improper #######
        for atom3 in mole:
            if atom3.atomtype==item.c:
                permu=list(itertools.permutations(atom3.neighbor,3))
                for tu in permu:
                    a= tu[0].atomtype==item.a or item.a=='*'
                    b= tu[1].atomtype==item.b or item.b=='*'
                    c= tu[2].atomtype==item.d or item.d=='*'
                    if a and b and c:
                        mole.addimproper(tu[0].atomnum,tu[1].atomnum,atom3.atomnum,tu[2].atomnum)
                        break
       #############################

        if item.forceconst==rxfile.mmfunction.unknownsign:
            nunk+=1
    else:
        if item.forceconst==rxfile.mmfunction.unknownsign:
            nunk+=1


# Prepare Head section
hessxyz=''
hprimexyz=''
finalxyz=''
for atom in mole:
    hprimexyz+=atom.atomsym+'-'+atom.name+'-'+'{:8.6f}'.format(float(atom.atomcharge))+'   '+'    '.join(["{: .12f}".format(x) for x in atom.coords])+'\n'
    hessxyz+=atom.atomsym+'-'+atom.name+'-0.000000'+'   '+'    '.join(["{: .12f}".format(x) for x in atom.coords])+'\n'
    finalxyz+=atom.atomsym+'-'+atom.atomtype+'-'+'{:8.6f}'.format(float(atom.atomcharge))+'   '+'    '.join(["{: .12f}".format(x) for x in atom.coords])+'\n'

hesshead=mmcom.commandline+'\nhess\n\n'+str(qmfile.totalcharge)+' '+str(qmfile.multiplicity)+'\n'+hessxyz+'\n'+mmcom.connectivity+'\n'
hprimehead=mmcom.commandline+'\nhprime\n\n'+str(qmfile.totalcharge)+' '+str(qmfile.multiplicity)+'\n'+hprimexyz+'\n'+mmcom.connectivity+'\n'
finalhead=mmcom.commandline+'\nfinal\n\n'+str(qmfile.totalcharge)+' '+str(qmfile.multiplicity)+'\n'+finalxyz+'\n'+mmcom.connectivity+'\n'

##  match MM function and internal coordinates

def matchdihd(dihd,func):
    a=(dihd[1].atomtype==func.a or func.a=='*')
    b=(dihd[2].atomtype==func.b or func.b=='*')
    c=(dihd[3].atomtype==func.c or func.c=='*')
    d=(dihd[4].atomtype==func.d or func.d=='*')
    forward=a and b and c and d
    a=(dihd[1].atomtype==func.d or func.d=='*')
    b=(dihd[2].atomtype==func.c or func.c=='*')
    c=(dihd[3].atomtype==func.b or func.b=='*')
    d=(dihd[4].atomtype==func.a or func.a=='*')
    backward=a and b and c and d
    if forward or backward:
        return True
    else:
        return False
def matchbond(bond,func):
    a=(bond[1].atomtype==func.a or func.a=='*')
    b=(bond[2].atomtype==func.b or func.b=='*')
    forward=a and b
    a=(bond[1].atomtype==func.b or func.b=='*')
    b=(bond[2].atomtype==func.a or func.a=='*')
    backward=a and b
    if forward or backward:
        return True
    else:
        return False
def matchangle(angle,func):
    a=(angle[1].atomtype==func.a or func.a=='*')
    b=(angle[2].atomtype==func.b or func.b=='*')
    c=(angle[3].atomtype==func.c or func.c=='*')
    forward=a and b and c
    a=(angle[1].atomtype==func.c or func.c=='*')
    b=(angle[2].atomtype==func.b or func.b=='*')
    c=(angle[3].atomtype==func.a or func.a=='*')
    backward=a and b and c
    if forward or backward:
        return True
    else:
        return False
def matchimproper(improper,func):
    a=(improper[1].atomtype==func.a or func.a=='*')
    b=(improper[2].atomtype==func.b or func.b=='*')
    c=(improper[3].atomtype==func.c or func.c=='*')
    d=(improper[4].atomtype==func.d or func.d=='*')
    forward=a and b and c and d
    if forward:
        return True
    else:
        return False

# Read vdw data
vdwdict=mmcom.vdwdict
hessvdwtail=''
hprimevdwtail=''

for atom in mole:
    atom.vdwradius=vdwdict[atom.atomtype][0]
    atom.vdwwelldepth=vdwdict[atom.atomtype][1]
    hessvdwtail+='VDW '+atom.name+'  '+atom.vdwradius+'  0.0000\n'
    hprimevdwtail+='VDW '+atom.name+'  '+atom.vdwradius+' '+atom.vdwwelldepth+'\n'



# Match internal coordinate and mmfunction
for dihd in mole.dihdlist.values():
        for dihdfunc in mmcom.nozomudihdfunc:
            if matchdihd(dihd,dihdfunc):
                dihd.nozomufunc=dihdfunc
                dihd.forceconst=[rxmol.dihdforceconst(x,dihd) for x in dihdfunc.forceconst]
                dihd.phase=copy.copy(dihdfunc.phase)
                dihd.npaths=dihdfunc.npaths
for angle in mole.anglelist.values():
        for anglefunc in mmcom.nozomuanglefunc:
            if matchangle(angle,anglefunc):
                angle.nozomufunc=anglefunc
                angle.forceconst=anglefunc.forceconst
for bond in mole.bondlist.values():
        for bondfunc in mmcom.nozomubondfunc:
            if matchbond(bond,bondfunc):
                bond.nozomufunc=bondfunc
                bond.forceconst=bondfunc.forceconst
for improper in mole.improperlist.values():
        for improperfunc in mmcom.nozomuimproperfunc:
            if matchimproper(improper,improperfunc):
                improper.nozomufunc=improperfunc
                improper.forceconst=improperfunc.forceconst


# Count unknown
itnlcordL=[]
itnlcordL.extend(mole.dihdlist.values())
itnlcordL.extend(mole.anglelist.values())
itnlcordL.extend(mole.bondlist.values())
itnlcordL.extend(mole.improperlist.values())
realnunk=0

for item in itnlcordL:
    if type(item)==rxmol.Dihd:
        for parms in item.forceconst:
            if str(parms)=='XXXXXX':
                realnunk+=1
    elif type(item)==rxmol.Angle:
        if item.forceconst=='XXXXXX':
            realnunk+=1
    elif type(item)==rxmol.Improper:
        if item.forceconst=='XXXXXX':
            realnunk+=1
    elif type(item)==rxmol.Bond:
        if item.forceconst=='XXXXXX':
            realnunk+=1

# Unit Hessian Component Tail
def hesstail(obj,i=0):
    tailstring=''
    for item in itnlcordL:
        if type(item)==rxmol.Dihd:
            parm=['0.000','0.000','0.000','0.000']
            if obj is item:
                parm[i]='1.000'
            tailstring+='AmbTrs  '+' '.join([x.center(3,' ') for x in item.repr.split()])+'  '+' '.join([str(x).center(3,' ') for x in item.phase])+'  '+' '.join([str(x) for x in parm])+'   '+str(item.npaths)+'\n'
        else:
            parm='0.000'
            if obj is item:
                parm='1.000'
            if type(item)==rxmol.Angle:
                tailstring+='HrmBnd1  '+' '.join([x.center(3,' ')for x in item.repr.split()])+'  '+parm+'  '+'{:>9.5f}'.format(item.anglevalue)+'\n'
            elif type(item)==rxmol.Bond:
                tailstring+='HrmStr1  '+' '.join([x.center(3,' ')for x in item.repr.split()])+'  '+parm+'  '+'{:>7.5f}'.format(item.length)+'\n'
            elif type(item)==rxmol.Improper:
                tailstring+='ImpTrs  '+' '.join([x.center(3,' ') for x in item.repr.split()])+'  '+parm+'  '+'{:6.2f}'.format(item.phase)+'  '+str(item.npaths)+'\n'

    for x in mmcom.additionfunc:
        tailstring+=x.content
    tailstring+=hessvdwtail
    tailstring+='\n\n'
    return tailstring


# Hprime(known) Hessian Component Tail
def hprimetail():
    tailstring=''
    for dihd in mole.dihdlist.values():
        parm=[]
        for item in dihd.forceconst:
            if str(item)=='XXXXXX':
                parm.append('0.000')
            else:
                parm.append(item)
        tailstring+='AmbTrs  '+' '.join([x.center(3,' ') for x in dihd.repr.split()])+'  '+' '.join([str(x).center(3,' ') for x in dihd.phase])+'  '+' '.join([str(x) for x in parm])+'   '+str(dihd.npaths)+'\n'
    for angle in mole.anglelist.values():
        if angle.forceconst=='XXXXXX':
            parm='0.000'
        else:
            parm=str(angle.forceconst)
        tailstring+='HrmBnd1  '+' '.join([x.center(3,' ')for x in angle.repr.split()])+'  '+parm+'  '+'{:>9.5f}'.format(angle.anglevalue)+'\n'
    for bond in mole.bondlist.values():
        if bond.forceconst=='XXXXXX':
            parm='0.000'
        else:
            parm=str(bond.forceconst)
        tailstring+='HrmStr1  '+' '.join([x.center(3,' ')for x in bond.repr.split()])+'  '+parm+'  '+'{:>7.5f}'.format(bond.length)+'\n'
    for improper in mole.improperlist.values():
        if improper.forceconst=='XXXXXX':
            parm='0.000'
        else:
            parm=str(improper.forceconst)
        tailstring+='ImpTrs  '+' '.join([x.center(3,' ') for x in improper.repr.split()])+'  '+parm+'  '+'{:6.2f}'.format(improper.phase)+'  '+str(improper.npaths)+'\n'
    for x in mmcom.additionfunc:
        tailstring+=x.content
    tailstring+=hprimevdwtail
    tailstring+='\n\n'
    return tailstring

# Prepare Unit Hessian File
hess=[]
num=realnunk
unkparmL=[]
for obj in itnlcordL:
    if type(obj)==rxmol.Dihd:
        obj.hessfile=[None,None,None,None]
        for i,parms in enumerate(obj.forceconst):
            if str(parms)=='XXXXXX':
                with open('hess'+str(len(hess))+'.com','w') as f:
                    f.write(hesshead+hesstail(obj,i))
                this=rxfile.File('hess'+str(len(hess)))
                obj.forceconst[i].hessfile=this
                this.orig=obj.forceconst[i]
                unkparmL.append(this.orig)
                hess.append(this)
                this.com.rung09()
                try:
                    this.com.isover()
                except:
                    this.com.rung09()
                    this.com.isover()
                this.runformchk()
                this.fchk.read()
                num-=1
                logging.info(str(num)+'+3 left')
    else:
         if obj.forceconst=='XXXXXX':
             with open('hess'+str(len(hess))+'.com','w') as f:
                 f.write(hesshead+hesstail(obj))
             this=rxfile.File('hess'+str(len(hess)))
             obj.hessfile=this
             this.orig=obj
             unkparmL.append(this.orig)
             hess.append(this)
             this.com.rung09()
             try:
                 this.com.isover()
             except:
                 this.com.rung09()
                 this.com.isover()
             this.runformchk()
             num-=1
             this.fchk.read()
             logging.info(str(num)+'+3 left')


#Identify Coupled Terms

links={}
for obj in unkparmL:
    if type(obj)==rxmol.dihdforceconst:
        a=obj.dihd[1].atomnum
        b=obj.dihd[4].atomnum
        if a>b: a,b=b,a
        if str(a)+'-'+str(b) in links.keys():
            links[str(a)+'-'+str(b)].append(obj)
        else:
            links.update({str(a)+'-'+str(b):[obj]})
    if type(obj)==rxmol.Angle:
        a=obj[1].atomnum
        b=obj[3].atomnum
        if a>b: a,b=b,a
        if str(a)+'-'+str(b) in links.keys():
            links[str(a)+'-'+str(b)].append(obj)
        else:
            links.update({str(a)+'-'+str(b):[obj]})
    if type(obj)==rxmol.Bond:
        a=obj[1].atomnum
        b=obj[2].atomnum
        if a>b: a,b=b,a
        if str(a)+'-'+str(b) in links.keys():
            links[str(a)+'-'+str(b)].append(obj)
        else:
            links.update({str(a)+'-'+str(b):[obj]})
    if type(obj)==rxmol.Improper:
        a=obj[1].atomnum
        b=obj[4].atomnum
        if a>b: a,b=b,a
        if str(a)+'-'+str(b) in links.keys():
            links[str(a)+'-'+str(b)].append(obj)
        else:
            links.update({str(a)+'-'+str(b):[obj]})

## store links
onetwoL={}
onetricL={}
onetriucL={}
onefourL={}
onetrifourL={}
for key,value in links.items():
    mytype=None
    types=[]
    for item in value:
        types.append(type(item))
    dboll=rxmol.dihdforceconst in types
    aboll=rxmol.Angle in types
    bboll=rxmol.Bond in types
    iboll=rxmol.Improper in types
    if dboll and aboll:
        onetrifourL.update({key:value})
    elif dboll:
        onefourL.update({key:value})
    elif aboll and iboll:
        onetricL.update({key:value})
    elif aboll:
        onetriucL.update({key:value})
    elif bboll:
        onetwoL.update({key:value})
        pass
    pass


# Start parameterization.


def calcgroup(adict,thishprime):
    for key,value in adict.items():
        leftL,rightL=[],[]
        a,b=[int(x) for x in key.split('-')]
        hq=qmfile.fchk.find33Hessian(a,b)
        hp=thishprime.fchk.find33Hessian(a,b)
        hideal=hq-hp
        hideal=[x for item in hideal for x in item]
        hk=[]
        for item in value:
            tmp=item.hessfile.fchk.find33Hessian(a,b)
            hk.append([x for row in tmp for x in row])
        for i in range(0,len(hideal)):
            rightL.append(hideal[i])
            line=[]
            for hks in hk:
                line.append(hks[i])
            leftL.append(line)
        leftL=np.array(leftL)
        rightL=np.array(rightL)
        res=np.linalg.lstsq(leftL,rightL)[0]
        for i,item in enumerate(value):
            item.forceconst=res[i]


# sequence: onefour-->onetrifour-->onetri(coupled)-->onetri(uncoupled)-->onetwo
if len(onefourL)!=0:
    onefourhprime=hprimehead+hprimetail()
    with open('onefourhprime.com','w') as f:
        f.write(onefourhprime)
    onefourhprime=rxfile.File('onefourhprime')
    onefourhprime.com.rung09()
    onefourhprime.com.isover()
    onefourhprime.runformchk()
    onefourhprime.fchk.read()
    calcgroup(onefourL,onefourhprime)
if len(onetrifourL)!=0:
    onetrifourhprime=hprimehead+hprimetail()
    with open('onetrifourhprime.com','w') as f:
        f.write(onetrifourhprime)
    onetrifourhprime=rxfile.File('onetrifourhprime')
    onetrifourhprime.com.rung09()
    onetrifourhprime.com.isover()
    onetrifourhprime.runformchk()
    onetrifourhprime.fchk.read()
    calcgroup(onetrifourL,onetrifourhprime)

if len(onetricL)!=0:
    onetrichprime=hprimehead+hprimetail()
    with open('onetrichprime.com','w') as f:
        f.write(onetrichprime)
    onetrichprime=rxfile.File('onetrichprime')
    onetrichprime.com.rung09()
    onetrichprime.com.isover()
    onetrichprime.runformchk()
    onetrichprime.fchk.read()
    calcgroup(onetricL,onetrichprime)


if len(onetriucL)!=0:
    onetriucLhprime=hprimehead+hprimetail()
    with open('onetriucLhprime.com','w') as f:
         f.write(onetriucLhprime)
    onetriucLhprime=rxfile.File('onetriucLhprime')
    onetriucLhprime.com.rung09()
    onetriucLhprime.com.isover()
    onetriucLhprime.runformchk()
    onetriucLhprime.fchk.read()
    calcgroup(onetriucL,onetriucLhprime)

if len(onetwoL)!=0:
    onetwoLhprime=hprimehead+hprimetail()
    with open('onetwoLhprime.com','w') as f:
        f.write(onetwoLhprime)
    onetwoLhprime=rxfile.File('onetwoLhprime')
    onetwoLhprime.com.rung09()
    onetwoLhprime.com.isover()
    onetwoLhprime.runformchk()
    onetwoLhprime.fchk.read()
    calcgroup(onetwoL,onetwoLhprime)



logging.info('Start Summarizing')
for nozomufunc in nozomuL:
    if nozomufunc.forceconst=='XXXXXX':
        res=0
        i=0
        for item in itnlcordL:
            if item.nozomufunc==nozomufunc:
                res+=item.forceconst
                i+=1
        nozomufunc.forceconst=res/i
        print(nozomufunc.type,nozomufunc.repr,nozomufunc.forceconst)
    if nozomufunc.type=='dihd':
        res=[0.0,0.0,0.0,0.0]
        i=0
        for item in itnlcordL:
            if item.nozomufunc==nozomufunc:
                res=[float(str(x))+oldx for x,oldx in zip(item.forceconst,res)]
                i+=1
        nozomufunc.forceconst=[x/i for x in res]
        print(nozomufunc.type,nozomufunc.repr,nozomufunc.forceconst)





tailstring=''
for dihd in mmcom.nozomudihdfunc:
    parm=[]
    for item in dihd.forceconst:
        if str(item)=='XXXXXX':
            parm.append('0.000')
            logging.critical('Force constant is not determined for dihedral '+dihd.repr)
            raise
        else:
            parm.append(float(str(item)))
    tailstring+='AmbTrs  '+' '.join([x.center(3,' ') for x in dihd.repr.split()])+'  '+' '.join([str(x).center(3,' ') for x in dihd.phase])+'  '+' '.join(['{:>6.3f}'.format(x) for x in parm])+'   '+str(dihd.npaths)+'\n'
for angle in mmcom.nozomuanglefunc:
    if angle.forceconst=='XXXXXX':
        parm='0.000'
        logging.critical('Force constant is not determined for angle '+angle.repr)
        raise
    else:
        parm=angle.forceconst
    tailstring+='HrmBnd1  '+' '.join([x.center(3,' ')for x in angle.repr.split()])+'  '+'{:>7.3f}'.format(parm)+'  '+'{:>9.5f}'.format(angle.eqvalue)+'\n'
for bond in mmcom.nozomubondfunc:
    if bond.forceconst=='XXXXXX':
        parm='0.000'
        logging.critical('Force constant is not determined for bond '+bond.repr)
        raise
    else:
        parm=bond.forceconst
    tailstring+='HrmStr1  '+' '.join([x.center(3,' ')for x in bond.repr.split()])+'  '+'{:>8.3f}'.format(parm)+'  '+'{:>7.5f}'.format(bond.eqvalue)+'\n'
for improper in mmcom.nozomuimproperfunc:
    if improper.forceconst=='XXXXXX':
        logging.critical('Force constant is not determined for improper '+improper.repr)
        raise
    else:
        parm=improper.forceconst
    tailstring+='ImpTrs  '+' '.join([x.center(3,' ') for x in improper.repr.split()])+'  '+'{:>7.3f}'.format(parm)+'  '+'{:6.2f}'.format(improper.phase)+'  '+str(improper.npaths)+'\n'
for x in mmcom.additionfunc:
    tailstring+=x.content
for vdw in mmcom.nozomuvdw:
    tailstring+='VDW  '+'  '+vdw.atomtype+'  '+vdw.radius+'  '+vdw.welldepth+'\n'
tailstring+='\n\n'

logging.info('\n\nresult:\n')
logging.info(tailstring)


with open('phf_result_'+originname,'w') as f:
    logging.info('Write result to file phf_result_'+originname)
    f.write(finalhead+tailstring)


if not delete:
    os.popen('mkdir tmpphffile')
    os.popen('mv *hprime* hess* del.sh tmpphffile')
else:
    os.popen('rm -rf *hprime* hess* gau* Hess* del.sh\n')



# Write Final result file

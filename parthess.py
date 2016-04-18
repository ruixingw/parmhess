#!/usr/bin/env python
from __future__ import print_function
import rx.molecules as rxmol
import rx.chemfiles as rxfile
import numpy as np
from io import StringIO
import argparse,os

parser=argparse.ArgumentParser()
parser.add_argument('mmfile',help="mmfile prepared by tsubasa.")
parser.add_argument('qmfile',help="FCHK file from a frequency calculation.")
parser.add_argument('-d','--delete',help="Delete all temp files. Default option will save files to ./phfdebug",action='store_true')
args=parser.parse_args()
mmfile=args.mmfile
originname=mmfile
qmfile=args.qmfile
delete=args.delete

print(mmfile,qmfile)
mmfile=mmfile[0:mmfile.find('.')]
qmfile=qmfile[0:qmfile.find('.')]
##prepare del.sh, del old files, print title
with open('del.sh','w') as f:
    f.write('rm -rf phfpy* tmp* *hprime* hess* gau* Hess* phfdebug del.sh\n')
os.system('chmod +x del.sh')

print('\n\n                         Partial Hessian Fitting for parameterization\n\n')
if delete:
    print("\n\nWARNING: Delete mode is on, all temperory files would be deleted.\n\n")
#os.system('rm -rf phfpy* gau* hprim* hess* Hess* ')
mole=rxmol.Molecule('mmfile')
mmfile=rxfile.File(mmfile).com
qmfile=rxfile.File(qmfile)
qmfile.readfchk()
mmfile.read()
nunk=0
for item in mmfile.nozomubondfunc:
    if item.value==rxfile.mmfunction.magicnum:
        nunk+=1
for item in mmfile.nozomuanglefunc:
    if item.value==rxfile.mmfunction.magicnum:
        nunk+=1
for item in mmfile.nozomudihdfunc:
    if item.value==rxfile.mmfunction.magicnum:
        nunk+=1

xyz=StringIO(mmfile.xyz)
mole.readfromxyz(xyz)
mole.readchargefromlist(mmfile.atomchargelist)
mole.readtypefromlist(mmfile.atomtypelist)
mole.readconnectivity(mmfile.connectivity)
hessxyz=''
hprimexyz=''
finalxyz=''
for atom in mole:
    hprimexyz+=atom.atomsym+'-'+atom.name+'-'+str(atom.atomcharge)+'   '+'    '.join(["{: .12f}".format(x) for x in atom.coords])+'\n'
    hessxyz+=atom.atomsym+'-'+atom.name+'-0.000000'+'   '+'    '.join(["{: .12f}".format(x) for x in atom.coords])+'\n'
    finalxyz+=atom.atomsym+'-'+atom.atomtype+'-'+str(atom.atomcharge)+'   '+'    '.join(["{: .12f}".format(x) for x in atom.coords])+'\n'

hesshead=mmfile.commandline+'\nhess\n\n'+str(qmfile.totalcharge)+' '+str(qmfile.multiplicity)+'\n'+hessxyz+'\n'+mmfile.connectivity+'\n'
hprimehead=mmfile.commandline+'\nhess\n\n'+str(qmfile.totalcharge)+' '+str(qmfile.multiplicity)+'\n'+hprimexyz+'\n'+mmfile.connectivity+'\n'
finalhead=mmfile.commandline+'\nfinal\n\n'+str(qmfile.totalcharge)+' '+str(qmfile.multiplicity)+'\n'+finalxyz+'\n'+mmfile.connectivity+'\n'
mole.readtypefromlist(mmfile.atomtypelist)

def matchdihd(dihd,func):
    a=(dihd.a.atomtype==func.a or func.a=='*')
    b=(dihd.b.atomtype==func.b or func.b=='*')
    c=(dihd.c.atomtype==func.c or func.c=='*')
    d=(dihd.d.atomtype==func.d or func.d=='*')
    forward=a and b and c and d
    a=(dihd.a.atomtype==func.d or func.d=='*')
    b=(dihd.b.atomtype==func.c or func.c=='*')
    c=(dihd.c.atomtype==func.b or func.b=='*')
    d=(dihd.d.atomtype==func.a or func.a=='*')
    reverse=a and b and c and d
    if forward or reverse:
        return True
    else:
        return False
def matchbond(bond,func):
    a=(bond.a.atomtype==func.a or func.a=='*')
    b=(bond.b.atomtype==func.b or func.b=='*')
    forward=a and b
    a=(bond.a.atomtype==func.b or func.b=='*')
    b=(bond.b.atomtype==func.a or func.a=='*')
    reverse=a and b
    if forward or reverse:
        return True
    else:
        return False
def matchangle(angle,func):
    a=(angle.a.atomtype==func.a or func.a=='*')
    b=(angle.b.atomtype==func.b or func.b=='*')
    c=(angle.c.atomtype==func.c or func.c=='*')
    forward=a and b and c
    a=(angle.a.atomtype==func.c or func.c=='*')
    b=(angle.b.atomtype==func.b or func.b=='*')
    c=(angle.c.atomtype==func.a or func.a=='*')
    reverse=a and b and c
    if forward or reverse:
        return True
    else:
        return False
pwd=os.path.split(os.path.realpath(__file__))[0]
vdwdata={}
with open(os.path.join(pwd,'vdw.dat'),'r') as f:
    for line in f:
        vdwdata.update({line.split()[0]:line.split()[1]+' '+line.split()[2]})
hessvdwtail=''
hprimevdwtail=''
for atom in mole:
    atom.vdwradius=vdwdata[atom.atomtype].split()[0]
    atom.vdwepsilon=vdwdata[atom.atomtype].split()[1]
    hessvdwtail+='VDW '+atom.name+'  '+atom.vdwradius+'  0.0000\n'
    hprimevdwtail+='VDW '+atom.name+'  '+atom.vdwradius+' '+atom.vdwepsilon+'\n'

for dihd in mole.dihdlist.values():
        for dihdfunc in mmfile.nozomudihdfunc:
            if matchdihd(dihd,dihdfunc):
                dihd.nozomufunc=dihdfunc
                dihd.parm=dihd.nozomufunc.value
for angle in mole.anglelist.values():
        for anglefunc in mmfile.nozomuanglefunc:
            if matchangle(angle,anglefunc):
                angle.nozomufunc=anglefunc
                angle.parm=anglefunc.value
for bond in mole.bondlist.values():
        for bondfunc in mmfile.nozomubondfunc:
            if matchbond(bond,bondfunc):
                bond.nozomufunc=bondfunc
                bond.parm=bondfunc.value
def tail(one,vdwtail,heorhp):
    hesstail=''
    for dihd in mole.dihdlist.values():
        if heorhp=='hess':
            parm='0.00 '
            if one is dihd:
                parm='1.00 '
        elif heorhp=='hprime':
            if dihd.parm=='XXXXXX':
                parm='0.00 '
            else:
                parm=str(dihd.parm)+' '
            if one is dihd:
                parm='0.00 '
        hesstail+='AmbTrs '+dihd.mydihdtype+' '
        if dihd.nozomufunc.periodicity==1:
            hesstail+=' '+str(dihd.nozomufunc.phase)+'0 0 0 '+parm+' 0.0 0.0 0.0 '+str(dihd.nozomufunc.npaths)+'\n'
        elif dihd.nozomufunc.periodicity==2:
            hesstail+='0 '+str(dihd.nozomufunc.phase)+' 0 0 0.0 '+parm+' 0.0 0.0 '+str(dihd.nozomufunc.npaths)+'\n'
        elif dihd.nozomufunc.periodicity==3:
            hesstail+='0 0 '+str(dihd.nozomufunc.phase)+' 0 0.0 0.0 '+parm+' 0.0 '+str(dihd.nozomufunc.npaths)+'\n'
        elif dihd.nozomufunc.periodicity==4:
            hesstail+='0 0 0 '+str(dihd.nozomufunc.phase)+' 0.0 0.0 0.0 '+parm+' '+str(dihd.nozomufunc.npaths)+'\n'
        else:
            print("No perio found for ",dihd.mydihdtype)
            quit()
    for angle in mole.anglelist.values():
        if heorhp=='hess':
            parm='0.00 '
            if one is angle:
                parm='1.00 '
        elif heorhp=='hprime':
            if angle.parm=='XXXXXX':
                parm='0.00 '
            else:
                parm=str(angle.parm)+' '
            if one is angle:
                parm='0.00 '

        hesstail+='HrmBnd1 '+angle.myangletype+' '+parm+"{: .3f}".format(angle.anglevalue)+'\n'

    for bond in mole.bondlist.values():
        if heorhp=='hess':
            parm='0.00 '
            if one is bond:
                parm='1.00 '
        elif heorhp=='hprime':
            if bond.parm=='XXXXXX':
                parm='0.00 '
            else:
                parm=str(bond.parm)+' '
            if one is bond:
                parm='0.00 '

        hesstail+='HrmStr1 '+bond.mybondtype+' '+parm+"{: .5f}".format(bond.length)+'\n'

    for addfunc in mmfile.additionfunc:
        if addfunc.content.find('ImpTrs')>=0 and heorhp=='hess':
            continue
        hesstail+=addfunc.content
    hesstail+=vdwtail
    hesstail+='\n\n'
    return hesstail


i=0
hess=[]
mmfunctions=[]
for obj in mole.dihdlist.values():
    mmfunctions.append(obj)
for obj in mole.anglelist.values():
    mmfunctions.append(obj)
for obj in mole.bondlist.values():
    mmfunctions.append(obj)

num=len(mmfunctions)
for obj in mmfunctions:
    with open('hess'+str(len(hess))+'.com','w') as f:
        f.write(hesshead+tail(obj,hessvdwtail,'hess'))
    this=rxfile.File('hess'+str(len(hess)))
    obj.hessfile=this
    hess.append(this)
    hess[-1].rung09()
    hess[-1].isover()
    num-=1
    print(num," left")
    hess[-1].runformchk()
    hess[-1].readfchk()
    i+=1

#Identify coupled terms

links={}
for obj in mmfunctions:
    if type(obj)==rxmol.Dihd:
        a=obj[1].atomnum
        b=obj[4].atomnum
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

for key,value in links.items():
    if len(value)>1:
        for item in links[key]:
            item.ifcouple=key
    else:
        for item in links[key]:
            item.ifcouple=False




# Dihds Hprime
#print(mole.dihd(1,2,3,4).parm)
dihdhprimetail=tail(mole.bond(1,2),hprimevdwtail,'hprime')
dihdhprime=hprimehead+dihdhprimetail
with open('dihdhprime.com','w') as f:
    f.write(dihdhprime)
dihdhprime=rxfile.File('dihdhprime')
dihdhprime.rung09()
dihdhprime.isover()
dihdhprime.runformchk()
dihdhprime.readfchk()

for dihd in mole.dihdlist.values():
    if dihd.parm=='XXXXXX' and dihd.ifcouple==False:
        a=dihd[1].atomnum
        b=dihd[4].atomnum
        hk=dihd.hessfile.find33Hessian(a,b)
        hp=dihdhprime.find33Hessian(a,b)
        hq=qmfile.find33Hessian(a,b)
        up=np.sum((hq-hp)*hk)
        down=np.sum(hk*hk)
        k=up/down
        left=[]
        right=[]
        hideal=hq-hp
        for i in range(0,3):
            for j in range(0,3):
                left.append([hk[i][j]])
                right.append(hideal[i][j])
        left=np.array(left)
        right=np.array(right)
        linres=np.linalg.lstsq(left,right)[0]
        dihd.parmlls=linres

        dihd.parm=k
    if dihd.parm=='XXXXXX' and dihd.ifcouple!=False:

        hk=[]
        a=int(dihd.ifcouple.split('-')[0])
        b=int(dihd.ifcouple.split('-')[1])
        for item in links[dihd.ifcouple]:
            hk.append(item.hessfile.find33Hessian(a,b))

        left=[]
        for ha in hk:
            left.append([np.sum(ha*hb) for hb in hk])
        left=np.array(left)

        hq=qmfile.find33Hessian(a,b)
        hp=dihdhprime.find33Hessian(a,b)
        right=[]
        for ha in hk:
            this=np.sum((hq-hp)*ha)
            right.append(this)
        right=np.array(right)
        result=np.linalg.solve(left,right)
        result=iter(result)
        llsleft=[]
        llsright=[]
        hideal=hq-hp
        for i in range(0,3):
            for j in range(0,3):
                llsleft.append([x[i][j] for x in hk])
                llsright.append(hideal[i][j])
        llsleft=np.array(llsleft)
        llsright=np.array(llsright)
        llsres=np.linalg.lstsq(llsleft,llsright)[0]
        llsres=iter(llsres)
        for item in links[dihd.ifcouple]:
            item.parm=next(result)
            item.parmlls=next(llsres)
#            item.parmlls=

for nozomufunc in mmfile.nozomudihdfunc:
    if nozomufunc.value=='XXXXXX':
        res=0
        i=0
        reslls=0
        for dihd in mole.dihdlist.values():
            if dihd.nozomufunc==nozomufunc:
                res+=dihd.parm
                reslls+=dihd.parmlls
                i+=1
        nozomufunc.value=res/i
        nozomufunc.valuells=reslls/i
        print("ori:",nozomufunc.value)
        print("lls:",nozomufunc.valuells)
# Angles Hprime
#print(mole.dihd(1,2,3,4).parm)
anglehprimetail=tail('123',hprimevdwtail,'hprime')
anglehprime=hprimehead+anglehprimetail
with open('anglehprime.com','w') as f:
    f.write(anglehprime)
anglehprime=rxfile.File('anglehprime')
anglehprime.rung09()
anglehprime.isover()
anglehprime.runformchk()
anglehprime.readfchk()

for angle in mole.anglelist.values():
    if angle.parm=='XXXXXX' and angle.ifcouple==False:
        a=angle[1].atomnum
        b=angle[3].atomnum
        hk=angle.hessfile.find33Hessian(a,b)
        hp=anglehprime.find33Hessian(a,b)
        hq=qmfile.find33Hessian(a,b)
        up=np.sum((hq-hp)*hk)
        down=np.sum(hk*hk)
        k=up/down
        hideal=hq-hp
        left=[]
        right=[]
        for i in range(0,3):
            for j in range(0,3):
                left.append([hk[i][j]])
                right.append(hideal[i][j])
        left=np.array(left)
        right=np.array(right)
        linres=np.linalg.lstsq(left,right)[0]
        angle.parmlls=linres

        angle.parm=k
    if angle.parm=='XXXXXX' and angle.ifcouple!=False:

        hk=[]
        a=int(angle.ifcouple.split('-')[0])
        b=int(angle.ifcouple.split('-')[1])
        for item in links[angle.ifcouple]:
            hk.append(item.hessfile.find33Hessian(a,b))

        left=[]
        for ha in hk:
            left.append([np.sum(ha*hb) for hb in hk])
        left=np.array(left)

        hq=qmfile.find33Hessian(a,b)
        hp=anglehprime.find33Hessian(a,b)
        right=[]
        for ha in hk:
            this=np.sum((hq-hp)*ha)
            right.append(this)
        right=np.array(right)
        result=np.linalg.solve(left,right)
        result=iter(result)

        llsleft=[]
        llsright=[]
        hideal=hq-hp
        for i in range(0,3):
            for j in range(0,3):
                llsleft.append([x[i][j] for x in hk])
                llsright.append(hideal[i][j])
        llsleft=np.array(llsleft)
        llsright=np.array(llsright)
        llsres=np.linalg.lstsq(llsleft,llsright)[0]
        llsres=iter(llsres)

        for item in links[angle.ifcouple]:
            item.parm=next(result)
            item.parmlls=next(llsres)

for nozomufunc in mmfile.nozomuanglefunc:
    if nozomufunc.value=='XXXXXX':
        res=0
        reslls=0
        i=0
        for angle in mole.anglelist.values():
            if angle.nozomufunc==nozomufunc:
                res+=angle.parm
                reslls+=angle.parmlls
                i+=1
        nozomufunc.value=res/i
        nozomufunc.valuells=reslls/i
        print("ori:",nozomufunc.value)
        print("lls:",nozomufunc.valuells)

# Bonds Hprime
#print(mole.dihd(1,2,3,4).parm)
bondhprimetail=tail('123',hprimevdwtail,'hprime')
bondhprime=hprimehead+bondhprimetail
with open('bondhprime.com','w') as f:
    f.write(bondhprime)
bondhprime=rxfile.File('bondhprime')
bondhprime.rung09()
bondhprime.isover()
bondhprime.runformchk()
bondhprime.readfchk()

for bond in mole.bondlist.values():
    if bond.parm=='XXXXXX' and bond.ifcouple==False:
        a=bond[1].atomnum
        b=bond[2].atomnum
        hk=bond.hessfile.find33Hessian(a,b)
        hp=bondhprime.find33Hessian(a,b)
        hq=qmfile.find33Hessian(a,b)
        up=np.sum((hq-hp)*hk)
        down=np.sum(hk*hk)
        k=up/down
        bond.eigenvalue=np.linalg.eig(hq)[0]
        bond.eigenvector=np.linalg.eig(hq)[1]
        vecu=bond.vec/np.linalg.norm(bond.vec)
        bond.proj=[abs(np.dot(vecu,x)) for x in bond.eigenvector]
        kab=[x*y for x,y in zip(bond.eigenvalue,bond.proj)]
        kab=np.sum(kab)
        kab=abs(kab*1185.82157)
        hideal=hq-hp
        left=[]
        right=[]
        middle=[]
        for i in range(0,3):
            for j in range(0,3):
                left.append([hk[i][j]])
                right.append(hideal[i][j])
                if hk[i][j]!=0:
                    middle.append(hideal[i][j]/hk[i][j])
        middle=np.array(middle)
        print(middle)
        print(hideal)
        print(hk)
        middle=np.sum(middle)/len(middle)
        left=np.array(left)
        right=np.array(right)
        linres=np.linalg.lstsq(left,right)[0]

        bond.parmlls=linres
        bond.parmeach=middle
        bond.parm=k
        bond.parmsemi=kab
    if bond.parm=='XXXXXX' and bond.ifcouple!=False:

        hk=[]
        a=int(bond.ifcouple.split('-')[0])
        b=int(bond.ifcouple.split('-')[1])
        for item in links[bond.ifcouple]:
            hk.append(item.hessfile.find33Hessian(a,b))

        left=[]
        for ha in hk:
            left.append([np.sum(ha*hb) for hb in hk])
        left=np.array(left)

        hq=qmfile.find33Hessian(a,b)
        hp=bondhprime.find33Hessian(a,b)
        right=[]
        for ha in hk:
            this=np.sum((hq-hp)*ha)
            right.append(this)
        right=np.array(right)
        result=np.linalg.solve(left,right)
        result=iter(result)

        llsleft=[]
        llsright=[]
        hideal=hq-hp
        for i in range(0,3):
            for j in range(0,3):
                llsleft.append([x[i][j] for x in hk])
                llsright.append(hideal[i][j])
        llsleft=np.array(llsleft)
        llsright=np.array(llsright)
        llsres=np.linalg.lstsq(llsleft,llsright)[0]
        llsres=iter(llsres)

        for item in links[bond.ifcouple]:
            item.parm=next(result)
            item.parmlls=next(llsres)

for nozomufunc in mmfile.nozomubondfunc:
    if nozomufunc.value=='XXXXXX':
        res=0
        ressemi=0
        reslls=0
        resmid=0
        i=0
        for bond in mole.bondlist.values():
            if bond.nozomufunc==nozomufunc:
                ressemi+=bond.parmsemi
                reslls+=bond.parmlls
                res+=bond.parm
                resmid=bond.parmeach
                i+=1
        nozomufunc.value=res/i
        nozomufunc.valuesemi=ressemi/i
        nozomufunc.valuells=reslls/i
        nozomufunc.valuemid=resmid/i
        print("ori:",nozomufunc.value)
        print("lls:",nozomufunc.valuells)
        print("middle:",nozomufunc.valuemid)


# Write final result file
finaltail=''
for item in mmfile.nozomudihdfunc:
    finaltail+='AmbTrs '+item.repr+' '
    parm="{: .3f}".format(item.value)
    if item.periodicity==1:
            finaltail+=' '+str(item.phase)+' 0 0 0 '+parm+' 0.0 0.0 0.0 '+str(item.npaths)+'\n'
    if item.periodicity==2:
            finaltail+='0 '+str(item.phase)+' 0 0 0.0 '+parm+' 0.0 0.0 '+str(item.npaths)+'\n'
    if item.periodicity==3:
            finaltail+='0 0 '+str(item.phase)+' 0 0.0 0.0 '+parm+' 0.0 '+str(item.npaths)+'\n'
    if item.periodicity==4:
            finaltail+='0 0 0 '+str(item.phase)+' 0.0 0.0 0.0 '+parm+' '+str(item.npaths)+'\n'


for item in mmfile.nozomuanglefunc:
    finaltail+='HrmBnd1 '+item.repr+' '
    parm="{: .3f}".format(item.value)
    finaltail+=' '+parm+' {: .4f}'.format(item.eqvalue)+'\n'
for item in mmfile.nozomubondfunc:
    finaltail+='HrmStr1 '+item.repr+' '
    parm="{: .3f}".format(item.value)
    finaltail+=' '+parm+' {: .5f}'.format(item.eqvalue)+'\n'

print(finaltail)

for addfunc in mmfile.additionfunc:
    finaltail+=addfunc.content
for nozovdw in mmfile.nozomuvdw:
    finaltail+=nozovdw.content

finaltail+='\n\n'
with open('phf_result_'+originname,'w') as f:
    f.write(finalhead+finaltail)


if not delete:
    os.system('mkdir phfdebug')
    os.system('mv phfpy* tmp* *hprime* hess* gau* Hess* del.sh phfdebug')
else:
    os.system('rm -rf phfpy* tmp* *hprime* hess* gau* Hess* del.sh\n')



# Write Final result file

#!/usr/bin/env python
from __future__ import print_function
from geomdef import *
from readfiles import *
import numpy as np
import argparse,os

parser=argparse.ArgumentParser()
parser.add_argument('inputfile',help="allinone input file")
parser.add_argument('-d','--debug',help="save all temp file for debug",action='store_true')
args=parser.parse_args()
input=args.inputfile
debug=args.debug

##prepare del.sh, del old files, print title
with open('del.sh','w') as f:
    f.write('rm -rf phfpy* hprim* tmp* hess* gau* Hess* del.sh\n')
os.system('chmod +x del.sh')

print('\n\n                         Partial Hessian Fitting for parameterization\n\n')
if debug:
    print("\n\nWARNING: Debug mode is on, all temperory files would be saved.\n\n")
os.system('rm -rf phfpy* gau* hprim* hess* Hess* ')

print('Read config from:',input)

##Read input. Read MMinput, QMfchk(freq) input, natoms, and Links to a list.
with open(input,'r') as f:
    for i,line in enumerate(f.readlines()):
        if line.find('#')==0:
            continue
        if line.find('#')>=0:
            line=line[:line.find('#')]
        if line.find('natoms')>=0:
            natoms=int(line.split('=')[1])
        if line.find('mmfile')>=0:
            mmfile=line.split('=')[1].strip('\n')
        if line.find('qmfchk')>=0:
            qmfchk=line.split('=')[1].strip('\n')
        if line.find('Link start')>=0:
            linkstart=i
    links=[]
    links.append([])
    p=0
    f.seek(0)
    for line in f.readlines()[linkstart+1:]:
        if line.find('#')==0:
            continue
        if line.find('#')>=0:
            line=line[:line.find('#')]
        if line.find('next')<0:
            links[p].append(line.strip('\n'))
        else:
            p=p+1
            links.append([])

    ifgroup=[]
    ngroup=0
    with open(mmfile,'r') as f:
        nunk=0
        for line in f.readlines():
            if line.find('XXXXXX')>=0:
                ifgroup.append('')
                if line.find('&&&')>=0:
                    ifgroup[nunk]=int(line.split('&&&')[1])
                    if ifgroup[nunk]>ngroup:
                        ngroup=ifgroup[nunk]
                else:
                    ifgroup[nunk]=0
                nunk+=1


print('natoms:',natoms)
print('mmfile:',mmfile)
print('qmfchk:',qmfchk)
print('Number of link blocks:',len(links)-1)
print('Number of unknown parameters:',nunk,'\n\n\n\n')
if len(links)-1!=nunk:
    print("Error in input.inp: Number of link blocks is not equal to unknowns")
    quit()

## Instantialize qmfchk
qmfchk=File(qmfchk.split('.')[0])



##Change g09 command
gauCOM.g09rt='g09'
## Prepare mmInput h(k) as hess[iparm]
print("Prepare hess:")
hcontent=[]
with open(mmfile,'r') as f:

    for iparm in range(0,nunk):
        f.seek(0)
        hcontent.append('')
        nX=-1
        for line in f.readlines():
            # Set charge,vdW to zero
            if line.find('VDW')>=0:
                line=line.split()[0]+'  '+line.split()[1]+'  0.0000\n'
            if line.find('-')>=0 and line.count('.')==4:
                line=line[:line.find('.')-1]+'0.000000'+line[line.find('.')+7:]
            if line.find('XXXXXX')>=0:
                nX+=1
                if iparm==nX:
                    line=line.split('XXXXXX')[0]+'1.000'+line.split('XXXXXX')[1]
                else:
                    line=line.split('XXXXXX')[0]+'0.000'+line.split('XXXXXX')[1]
            else:
                if line.find('AmbTrs')>=0:
                    continue
                elif line.find('HrmBnd1')>=0:
                    line=line.split()[0]+' '+line.split()[1]+' '+line.split()[2]+' '+line.split()[3]+' 0.000 '+line.split()[5]+'\n'
                elif line.find('HrmStr1')>=0:
                    line=line.split()[0]+' '+line.split()[1]+' '+line.split()[2]+' 0.000 '+line.split()[4]+'\n'
                elif line.find('ImpTrs')>=0:
                    line=line.split()[0]+' '+line.split()[1]+' '+line.split()[2]+' '+line.split()[3]+' '+line.split()[4]+' 0.000 '+line.split()[6]+' '+line.split()[7]+'\n'
            hcontent[iparm]+=line.strip('&&&').strip('\n')+'\n'
        with open('hess'+str(iparm)+'.com','w') as w:
            w.write(hcontent[iparm])

## Instantialize hess[iparm], Run g09
hess=[]
for iparm in range(0,nunk):
    hess.append(File('hess'+str(iparm)))
    hess[iparm].rung09()
for iparm in range(0,nunk):
    if not hess[iparm].isover():
        quit()
    else:
        if not hess[iparm].formchk():
            quit()

print('hess prepared.\n\n')




# Start Calculate Parameters:

# Read QMHessian:
if not qmfchk.readHessianFromFchk():
    print("Failure to read QMHessian from:",qmfchk.name)
    quit()


# Read h(k) Hessian:
for iparm in range(0,nunk):
    if not hess[iparm].readHessianFromFchk():
        print("Failure to read MMHessian",iparm,"from:",hess[iparm].fchkname())
        quit()


# Start Loop: Calculate every hprime and parameter

hprime=[]
hprimegrp=['']
# An Hprime is saved to prepare next hprime
Hprime=''
with open(mmfile,'r') as f:
    f.seek(0)
    Hprime=[x.strip('&&&').strip('\n')+'\n' for x in f.readlines()]


maxgroup=ngroup
# igroup : group number
# class linearequation(object):
#     def __init__(self):
#         self.matA=[]
#         self.vecB=[]
#         self.current=0
#         self.high=0
#     def print(self):
#         print("Coeffiecnt matrix is,",self.matA)
#         print("Corresponding vecB is,",self.vecB)
#     def __iter__(self):
#         return self
#     def next(self):
#         if self.current>=self.high:
#             raise StopIteration
#         else:
#             self.current+=1
#             return [self.matA[self.current-1],self.vecB[self.current-1]]



def recreqs(A):  # A[i][0,1] ; 0:matA, 1:vecB; i: i blocks
    matCalc=[]
    vecCalc=[]
    item=[]
    ig=0
    result=0
    nblocks=len(A[0])
    matCalc=list(zip(*[A[i][0] for i in range(0,len(A))]))
    vecCalc=list(zip(*[A[i][1] for i in range(0,len(A))]))

    if debug:
        print("matCalc is:\n",matCalc)
        print("vecCalc is:\n",vecCalc)
        print("\nstart")
    for i in range(0,len(matCalc)):
        if debug:
            print("Selected matCalc is:\n",matCalc[i])
            print("Selected vecCalc is:\n",vecCalc[i])
        rsprint=np.linalg.solve(matCalc[i],vecCalc[i])
        result+=rsprint
        print(rsprint)
        if debug:
            print('stop\n')
        ig+=1
    print("total result=",result)
    result=result/ig
    print("averaged result=",result)
    print("\nthis group had used",ig,"linear equation systems\n")
    return result




for igroup in range(1,maxgroup+1):
    print('Group NO.',igroup)
    # Prepare hprim: set all unknowns to 0.000
    content=''

    for line in Hprime:
        if line.find('XXXXXX')>=0:
            line=line.split('XXXXXX')[0]+'0.000'+line.split('XXXXXX')[1]
        content+=line

    with open('hprimegrp'+str(igroup)+'.com','w') as f:
        f.write(content)

    # Instantialize hprimegrp, rung09, and formchk
    hprimegrp.append(File('hprimegrp'+str(igroup)))
    hprimegrp[igroup].rung09()
    if not hprimegrp[igroup].isover():
        print("Error in hprim:",hprimegrp[igroup].comname())
    else:
        if not hprimegrp[igroup].formchk():
            quit()

    # Read hprimegrp Hessian
    if not hprimegrp[igroup].readHessianFromFchk():
        print("Failure to read Hessian from:",hprimegrp[igroup].fchkname())
        quit()
    # Filter  parameters of this group
    nlink=0
    matA=[]
    vecB=[]
    iparmingroup=-1
    for iparm in range(0,nunk):
        if ifgroup[iparm]!=igroup:
            continue
        matA.append([])
        vecB.append([])
        iparmingroup+=1
        for num,x in enumerate(links[iparm]): ## each atom pair
            i=int(x.split('-')[0])
            j=int(x.split('-')[-1])
            if debug:
                print("Atom i:",i)
                print("Atom j:",j)
     #       print(i,j)
            # read 33 Hessian, instantialize to numpy array
            Hqm=np.array(qmfchk.find33Hessian(i,j),dtype=float)
            Hpm=np.array(hprimegrp[igroup].find33Hessian(i,j),dtype=float)
            hk=np.array(hess[iparm].find33Hessian(i,j),dtype=float)
            if debug:
                print("Hqm=\n",np.array(Hqm))
                print("Hpm=\n",np.array(Hpm))
                print("hk=\n",np.array(hk))
            vecB[iparmingroup].append(np.sum((Hqm-Hpm)*hk))

            matA[iparmingroup].append([])
            for ipa in range(0,nunk):
                    if ifgroup[ipa]==igroup:
                         i=int(x.split('-')[0])
                         j=int(x.split('-')[-1])
                         hb=np.array(hess[ipa].find33Hessian(i,j),dtype=float)
                         if debug:
                             print("hb=\n",np.array(hb))
                             print("hb*hk=\n",np.sum(hb*hk))
                         matA[iparmingroup][num].append(np.sum(hb*hk))

    A=[]
    for i in range(0,iparmingroup+1):    # A[igroup][0]:mat  [1]:vec
        A.append([])
        A[i].append(matA[i])
        A[i].append(vecB[i])
        if debug:
            print("Coeff:\n")
            print(np.array(A[i][0]))
            print("\nvecB:\n")
            print(np.array(A[i][1]))
            print("Another\n")
    result=0
#    ig=0
    result=recreqs(A)
    # print("Total:",result)
    # result=result/ig
    # print("Average:",result)



    # Update Hprime
    newHprime=[]
    i=0
    for line in Hprime:
        if line.find('&&&')>=0 and line.find('XXXXXX')>=0:
            if int(line.split('&&&')[1])==igroup:
                line=line.split('XXXXXX')[0]+"{:.3f}".format(result[i])+line.split('XXXXXX')[1]
                i+=1
        newHprime.append(line)
    Hprime=newHprime





for iparm in range(0,nunk):

    print('NO.',iparm)

    # Prepare hprim: set all unknowns to 0.000
    content=''

    for line in Hprime:
        if line.find('XXXXXX')>=0:
            line=line.split('XXXXXX')[0]+'0.000'+line.split('XXXXXX')[1]
        content+=line

    with open('hprime'+str(iparm)+'.com','w') as f:
        f.write(content)

    # Instantialize hprime, rung09, and formchk
    hprime.append(File('hprime'+str(iparm)))
    if ifgroup[iparm]!=0:
        continue
    hprime[iparm].rung09()
    if not hprime[iparm].isover():
        print("Error in hprim:",hprime[iparm].comname())
    else:
        if not hprime[iparm].formchk():
            quit()

    # Read hprime Hessian
    if not hprime[iparm].readHessianFromFchk():
        print("Failure to read Hessian from:",hprime[iparm].fchkname())
        quit()

    # Start calculate parameter: group 0
    result=0
    for x in links[iparm]: ## each atom pair
        up=0
        down=0
        i=int(x.split('-')[0])
        j=int(x.split('-')[-1])
        # read 3 by 3 Hessian, instantialize to numpy array
        Hqm=np.array(qmfchk.find33Hessian(i,j),dtype=float)
        Hpm=np.array(hprime[iparm].find33Hessian(i,j),dtype=float)
        hk=np.array(hess[iparm].find33Hessian(i,j),dtype=float)

        up=np.sum((Hqm-Hpm)*hk)
        down=np.sum(hk*hk)
        result+=up/down
        if debug:
            print('Hqm:\n',np.array(Hqm))
            print('Hprime:\n',np.array(Hpm))
            print('hk:\n',np.array(hk))
            print('parameter:',up/down,'\n')

    ## Average every atom pairs
    result=result/(len(links[iparm]))
    print('{:3f}'.format(result))
    print("Averaged from",len(links[iparm]),"results")

    ## Update Hprime for next cycle
    newHprime=[]
    ifchange=True
    for line in Hprime:
        if line.find('XXXXXX')>=0 and ifchange:
            line=line[:line.find('XXXXXX')]+'{: .3f}'.format(result)+line[line.find('XXXXXX')+6:]
            ifchange=False
        newHprime.append(line)
    Hprime=newHprime

## Write final result
with open('phf_result_'+mmfile,'w') as f:
    for line in Hprime:
        f.write(line)


if debug:
    print("WARNING: Debug mode is on, all temperory files were saved.")
else:
    os.system('rm -rf gau* hprim* hess* Hess* del.sh')

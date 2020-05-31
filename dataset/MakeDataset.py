#!/usr/bin/env python
''' make a dataset of RNA chains '''
import os
import sys

Lmin=30      # L>=Lmin
Lmax=750     # L<=Lmax
resolu="4.0" # 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 20.0, all
PairMin=10   # Minimal base pair number
c=0.8        # cd-hit-est sequence identity cutoff     

datasetdir=os.path.dirname(os.path.abspath(__file__))
rootdir   =os.path.dirname(datasetdir)
culldir   =os.path.join(rootdir,"cull")
cdhitest  =os.path.join(rootdir,"script/cd-hit-est")

#### list of targets with >=20 atoms per residue ####
atomic_target_list=[]
fp=open(os.path.join(culldir,"all_c1.0_s1.0/list.atomic"),'r')
for line in fp.read().splitlines():
    atomic_target_list.append(line.split()[0])
fp.close()
atomic_target_set=set(atomic_target_list)

#### number of base pairs per target ####
filename=os.path.join(datasetdir,"ListOfBasePairs")
dssrnames=os.path.join(culldir,"all_c1.0_s1.0/DSSR/*.dssr")
cmd='grep "List of" %s > %s'%(dssrnames,filename)
print(cmd)
os.system(cmd)
basepair_dict=dict()
fp=open(filename,'r')
for line in fp.read().splitlines():
    target,basepairs=line.split(':')
    target=os.path.basename(target).split('.')[0]
    basepairs=int(basepairs.split()[2])
    basepair_dict[target]=basepairs
fp.close()

#### list of targets with required length, resolution ####
target_list=[]
sequence_dict=dict()
fp=open(os.path.join(culldir,"pdb_atom.sort.%s_c1.0_s1.0"%resolu),'r')
for block in ('\n'+fp.read()).split('\n>'):
    lines=block.splitlines()
    if len(lines)<2:
        continue
    target=lines[0].split()[0]
    if not target in atomic_target_set:
        print("%s do not have full atoms"%(target))
        continue
    if not target in basepair_dict:
        print("%s do not have dssr"%(target))
        continue
    if basepair_dict[target]<PairMin:
        print("%s base pairs=%d"%(target,basepair_dict[target]))
        continue
    sequence=''.join(lines[1:])
    Lch=len(sequence)
    if Lch<Lmin or Lch>Lmax:
        print("%s Lch=%d"%(target,Lch))
        continue
    target_list.append(target)
    sequence_dict[target]=sequence

#### perform cd-hit-est ####
txt=''
for target in target_list:
    txt+=">%s\n%s\n"%(target,sequence_dict[target])
cdhitinfile=os.path.join(datasetdir,"%s_c1.0.fasta"%resolu)
cdhitoutfile=os.path.join(datasetdir,"%s_c%.1f"%(resolu,c))
fp=open(cdhitinfile,'w')
fp.write(txt)
fp.close()
cmd="%s -r 0 -M 5000 -n 5 -c %f -i %s -o %s"%(
    cdhitest,c,cdhitinfile,cdhitoutfile)
print(cmd)
os.system(cmd)

os.system("rm %s"%cdhitinfile)
os.system("rm %s.clstr"%cdhitoutfile)

#### print list ####
target_list=[]
len_dict=dict()
fp=open(cdhitoutfile,'r')
for block in ('\n'+fp.read()).split('\n>'):
    lines=block.splitlines()
    if len(lines)<2:
        continue
    target=lines[0].split()[0]
    sequence=''.join(lines[1:])
    len_dict[target]=len(sequence)
    target_list.append(target)
fp.close()

txt=''
for target in target_list:
    txt+="%s\t%d\n"%(target,len_dict[target])

fp=open(cdhitoutfile+".len",'w')
fp.write(txt)
fp.close()

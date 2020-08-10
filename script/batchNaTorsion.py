#!/usr/bin/env python
docstring='''
batchNaTorsion.py ../cull/all_c1.0_s1.0/list.atomic ../cull/pdb_atom.sort.4.0_c0.8_s0.8 ../cull/all_c1.0_s1.0/
    Read target list ../cull/all_c1.0_s1.0/list.atomic
    fasta file ../cull/pdb_atom.sort.4.0_c0.8_s0.8, and
    pdb folder ../cull/all_c1.0_s1.0/
    get all torsions at torsion.txt
'''

import sys, os
import textwrap
import subprocess

bindir=os.path.dirname(os.path.abspath(__file__))
NaTorsion_exe=os.path.join(bindir,"NaTorsion")

def batchNaTorsion(inputfasta,target_set,pdbfolder):
    len_dict=dict()
    target_list=[]
    txt=''
    fp=open(inputfasta,'r')
    for block in ('\n'+fp.read()).split('\n>'):
        lines=block.splitlines()
        if len(lines)<2:
            continue
        target=lines[0].split()[0]
        if not target in target_set:
            print("skip unlisted target %s"%target)
        filename=os.path.join(pdbfolder,target+".pdb")
        if not os.path.isfile(filename):
            print("Warning! %s missing"%filename)
            continue
        sequence=''.join(lines[1:])
        len_dict[target]=len(sequence)
        target_list.append(target)
        cmd="%s %s/%s.pdb 7"%(NaTorsion_exe,pdbfolder,target)
        p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        stdout,stderr=p.communicate()
        txt+="#### %s ####\n"%target
        txt+=stdout
    fp.close()
    fp=open("NaTorsion.raw",'w')
    fp.write(txt)
    fp.close()
    os.system("gzip -f NaTorsion.raw")
    return target_list,len_dict

if __name__=="__main__":
    if len(sys.argv)<=3:
        sys.stderr.write(docstring)
        exit()

    listfile  =sys.argv[1]
    inputfasta=sys.argv[2]
    pdbfolder=sys.argv[3]
    
    fp=open(listfile,'r')
    target_set=set([line.split()[0] for line in fp.read().splitlines()])
    fp.close()

    target_list,len_dict=batchNaTorsion(inputfasta,target_set,pdbfolder)
    print("%d targets"%len(target_list))
    print("%d nucleotides"%sum(len_dict.values()))

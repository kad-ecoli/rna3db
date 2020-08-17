#!/usr/bin/env python
docstring='''
makeDimerList.py pdb_atom.fasta list.multichain pdb/data/structures/all/pdb list.nodimer dimer.fasta
    Extract all pairs with chain with L>=10 residues from pdb_atom.fasta.
    Output the chains from the same pdb to list.multichain.
    Calculate inter-chain secondary structures for all listed multichain
    pdb using x3dna-dssr and cssr. Output list of pdb without dimer chain
    to list.nodimer. output the sequence of dimer chain pair with
    at least 5 residue pairs by either dssr or cssr to dimer.fasta, whose
    header has the following fields:
    >pdbID:chainID1:chainID2	length	dssr_base_pairs	cssr_base_pairs
'''

import os, sys
import subprocess
from clean_pdb import readPDB, parseHET, parseNoBB

dssr_exe=os.path.join(os.path.dirname(os.path.abspath(__file__)),"x3dna-dssr")
cssr_exe=os.path.join(os.path.dirname(os.path.abspath(__file__)),"cssr")

def makeMultiChainList(fastaFile,multiChainFile):
    fp=open(fastaFile,'r')
    blocks=fp.read().lstrip('>').split('\n>')
    fp.close()

    pdb_dict=dict()
    seq_dict=dict()
    for block in blocks:
        chainID,sequence=block.strip().splitlines()
        if len(sequence)<10:
            continue
        pdbID=chainID[:4]
        if not pdbID in pdb_dict:
            pdb_dict[pdbID]=[]
        pdb_dict[pdbID].append(chainID)
        seq_dict[chainID]=sequence
    
    for pdbID in pdb_dict.keys():
        if len(pdb_dict[pdbID])<=1:
            for chainID in pdb_dict[pdbID]:
                del seq_dict[chainID]
            del pdb_dict[pdbID]

    txt=''
    for pdbID in pdb_dict:
        chainList=','.join([chain[4:] for chain in pdb_dict[pdbID]])
        txt+="%s\t%s\n"%(pdbID,chainList)
    fp=open(multiChainFile,'w')
    fp.write(txt)
    fp.close()
    return pdb_dict,seq_dict

def ss_merge_dimer(chainFile1,chainFile2):    
    KeepChainNum=1
    removeH=True
    fixHET=2
    delNoBB=3

    PDBtxt1=readPDB(chainFile1,KeepChainNum,removeH,'A')
    PDBtxt1=parseHET(PDBtxt1,fixHET)
    PDBtxt1=parseNoBB(PDBtxt1,delNoBB)

    PDBtxt2=readPDB(chainFile2,KeepChainNum,removeH,'B')
    PDBtxt2=parseHET(PDBtxt2,fixHET)
    PDBtxt2=parseNoBB(PDBtxt2,delNoBB)
    
    p=subprocess.Popen(cssr_exe+" - ",
        shell=True,stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    stdout,stderr=p.communicate(input=PDBtxt1+PDBtxt2)
    bp_cssr=0
    for line in stdout.splitlines():
        if not "WC         " in line and \
           not "Wobble     " in line:
            continue
        bp_cssr+=(line[5:7]=="A." and line[20:22]=="B.")
        bp_cssr+=(line[5:7]=="B." and line[20:22]=="A.")
    if bp_cssr==0:
        return 0,0
    
    p=subprocess.Popen(dssr_exe+" -i=stdin --format=pdb --pair-only ",
        shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate(input=PDBtxt1+PDBtxt2)
    bp_dssr=0
    for line in stdout.splitlines():
        if not "WC         " in line and \
           not "Wobble     " in line:
            continue
        bp_dssr+=(line[5:7]=="A." and line[20:22]=="B.")
        bp_dssr+=(line[5:7]=="B." and line[20:22]=="A.")
    return bp_dssr,bp_cssr

def makeDimerFasta(pdb_dict,seq_dict,prefix,noDimerFile,dimerChainFile):
    parsed_pdb=[]
    if os.path.isfile(noDimerFile) and os.path.isfile(dimerChainFile):
        fp=open(noDimerFile,'r')
        parsed_pdb=fp.read().splitlines()
        fp.close()
        fp=open(dimerChainFile,'r')
        for line in fp.read().splitlines():
            if line.startswith('>'):
                pdbID=line[1:5]
                if not pdbID in parsed_pdb:
                    parsed_pdb.append(pdbID)
        fp.close()
        print("%d previously parsed pdb entries"%len(parsed_pdb))
    
    noDimerTxt   =''
    dimerChainTxt=''
    for pdbID,chainList in pdb_dict.items():
        if pdbID in parsed_pdb:
            continue
        print(pdbID+'\t'+','.join([chain[4:] for chain in pdb_dict[pdbID]]))
        txt=''
        for c1 in range(len(chainList)):
            chain1=chainList[c1]
            chainFile1=os.path.join(prefix, pdbID, chain1+".pdb.gz")
            for c2 in range(c1+1,len(chainList)):
                chain2=chainList[c2]
                chainFile2=os.path.join(prefix, pdbID, chain2+".pdb.gz")
                bp_dssr,bp_cssr=ss_merge_dimer(chainFile1,chainFile2)
                if bp_dssr<5 and bp_cssr<5:
                    continue
                print("Found dimer: %s\t%s"%(chain1,chain2))
                sequence1=seq_dict[chain1]+'\n'+seq_dict[chain2]
                sequence2=seq_dict[chain2]+'\n'+seq_dict[chain1]
                L=len(sequence1)-1
                txt+=">%s:%s:%s\t%d\t%d\t%d\n%s\n"%(
                    pdbID,chain1[4:],chain2[4:],L,bp_dssr,bp_cssr,sequence1)
                if sequence1==sequence2:
                    continue
                txt+=">%s:%s:%s\t%d\t%d\t%d\n%s\n"%(
                    pdbID,chain2[4:],chain1[4:],L,bp_dssr,bp_cssr,sequence2)
        if len(txt):
            dimerChainTxt+=txt
        else:
            noDimerTxt+=pdbID+'\n'

    fp=open(noDimerFile,   'a' if len(parsed_pdb) else 'w')
    fp.write(noDimerTxt)
    fp.close()
    fp=open(dimerChainFile,'a' if len(parsed_pdb) else 'w')
    fp.write(dimerChainTxt)
    fp.close()
    return

if __name__=="__main__":
    if len(sys.argv)!=6:
        sys.stderr.write(docstring)
        exit()

    fastaFile     =sys.argv[1]
    multiChainFile=sys.argv[2]
    prefix        =sys.argv[3]
    noDimerFile   =sys.argv[4]
    dimerChainFile=sys.argv[5]

    pdb_dict,seq_dict=makeMultiChainList(fastaFile,multiChainFile)
    makeDimerFasta(pdb_dict,seq_dict,prefix,noDimerFile,dimerChainFile)

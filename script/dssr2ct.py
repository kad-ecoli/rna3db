#!/usr/bin/env python
docstring='''
dssr2ct.py input.fasta input.dssr output.ct
    convert DSSR format input to CT format output
'''
import sys
import re

'''
CT File Format
A CT (Connectivity Table) file contains secondary structure information for a sequence. These files are saved with a CT extension. When entering a structure to calculate the free energy, the following format must be followed.

Start of first line: number of bases in the sequence
End of first line: title of the structure
Each of the following lines provides information about a given base in the sequence. Each base has its own line, with these elements in order:
Base number: index n
Base (A, C, G, T, U, X)
Index n-1
Index n+1
Number of the base to which n is paired. No pairing is indicated by 0 (zero).
Natural numbering. RNAstructure ignores the actual value given in natural numbering, so it is easiest to repeat n here.
'''

def read_fasta(fastafile):
    sequence=''
    name='dssr2ct'
    txt=''
    if fastafile=='-':
        txt=sys.stdin.read()
    else:
        fp=open(fastafile,'r')
        txt=fp.read()
        fp.close()
    for line in txt.splitlines():
        if line.startswith('>'):
            name=line[1:].split()[0]
        else:
            sequence+=line.upper().replace('T','U')
    return name,sequence

def read_dssr(dssrfile):
    pat=re.compile("^\s*\d+\s[ACGTU](\d+)\s+[ACGTU](\d+)\s+[ACGTU][-][ACGTU]\s+((Wobble )|(WC ))")
    pair_dict=dict()
    fp=open(dssrfile,'r')
    for line in fp.read().splitlines():
        if pat.match(line):
            items=line.strip().split()
            nt1=int(items[1][1:])
            nt2=int(items[2][1:])
            pair_dict[nt1]=nt2
            pair_dict[nt2]=nt1
    fp.close()
    return pair_dict

def dssr2ct(name,sequence,pair_dict,ctfile):
    txt='%5d %s\n'%(len(pair_dict)/2,name)
    for i,nt in enumerate(sequence):
        nt1=i+1
        nt2=0
        nextnt=nt1+1
        if nextnt>len(sequence):
            nextnt=0
        if nt1 in pair_dict:
            nt2=pair_dict[nt1]
        txt+="%5d %s %5d %5d %5d %5d\n"%(nt1,nt,nt1-1,nextnt,nt2,nt1)
    if ctfile=='-':
        sys.stdout.write(txt)
    else:
        fp=open(ctfile,'w')
        fp.write(txt)
        fp.close()
    return

if __name__=="__main__":
    if len(sys.argv)<3:
        sys.stderr.write(docstring)
        exit()

    fastafile=sys.argv[1]
    dssrfile =sys.argv[2]
    ctfile   ="-"
    if len(sys.argv)>3:
        ctfile=sys.argv[3]

    name,sequence=read_fasta(fastafile)
    pair_dict=read_dssr(dssrfile)
    dssr2ct(name,sequence,pair_dict,ctfile)


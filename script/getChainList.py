#!/usr/bin/env python
docstring='''
getChainList na_seqres.txt na_chain.list
    From input fasta na_seqres.txt, get the list of PDB chains and output it to
    chain.list.
'''
import sys

if len(sys.argv)<3:
    sys.stderr.write(docstring)
    exit()

infile=sys.argv[1]
outfile=sys.argv[2]

fp=open(infile,'r')
lines=fp.read().splitlines()
fp.close()

txt=''
for line in lines:
    if line.startswith('>'):
        chain=line.split()[0].lstrip('>').replace('_','')
        txt+=chain+'\n'

fp=open(outfile,'w')
fp.write(txt)
fp.close()

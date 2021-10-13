#!/usr/bin/env python
docstring='''
extractNA pdb_seqres.txt.gz na_seqres.txt
    From input fasta pdb_seqres.txt.gz, extract all nucleotide sequences,
    and save the sequences to output fasta na_seqres.txt
'''
import sys
import gzip

if len(sys.argv)!=3:
    sys.stderr.write(docstring)
    exit()

infile=sys.argv[1]
outfile=sys.argv[2]

if infile.endswith(".gz"):
    fp=gzip.open(infile,'r')
else:
    fp=open(infile,'r')
blocks=fp.read().lstrip('>').split('\n>')
fp.close()

txt=''
for block in blocks:
    lines=block.splitlines()
    if len(lines)<2:
        continue
    if lines[0].split()[1]!="mol:na":
        continue
    sequence=''.join(lines[1:])
    if len(sequence)<10:
        continue
    txt+='>'+lines[0]+'\n'+sequence+'\n'

fp=open(outfile,'w')
fp.write(txt)
fp.close()

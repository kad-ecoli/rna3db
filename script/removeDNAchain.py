#!/usr/bin/env python
docstring='''
removeDNAchain.py na_type.list
    Read molecule type definition na_type.list, and
    remove files for DNA-RNA, DNA/RNA and DNA.
'''
import sys, os

if len(sys.argv)<2:
    sys.stderr.write(docstring)
    exit()

infile=sys.argv[1]

dna_chain_list=[]
fp=open(infile,'r')
for line in fp.read().splitlines():
    chain,mol=line.split()[:2]
    if not mol.startswith("RNA"):
        dna_chain_list.append(chain)
fp.close()

dna_file_list=[]
for chain in dna_chain_list:
    filename=chain+".pdb"
    if os.path.isfile(filename):
        dna_file_list.append(filename)

print("removing %d DNA, DNA/RNA and DNA-RNA chains"%len(dna_file_list))
for filename in dna_file_list:
    os.system("rm "+filename)

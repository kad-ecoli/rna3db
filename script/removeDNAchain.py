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
rna_pdb_list=[]
fp=open(infile,'r')
for line in fp.read().splitlines():
    chain,mol=line.split()[:2]
    if mol.startswith("RNA"):
        rna_pdb_list.append(chain[:4])
    else:
        dna_chain_list.append(chain)
fp.close()
rna_pdb_set=set(rna_pdb_list)

dna_file_list=[]
dna_pdb_list=[]
for chain in dna_chain_list:
    pdb=chain[:4]
    for filename in [chain+".pdb",
                     chain+".pdb.gz",
             pdb+'/'+chain+".pdb",
             pdb+'/'+chain+".pdb.gz"]:
        if os.path.isfile(filename):
            dna_file_list.append(filename)
    if not pdb in rna_pdb_set and not pdb in dna_pdb_list and os.path.isdir(pdb):
            dna_pdb_list.append(pdb)

print("removing %d DNA, DNA/RNA and DNA-RNA chains"%len(dna_file_list))
for filename in dna_file_list:
    os.system("rm "+filename)

print("removing %d DNA only folders"%len(dna_pdb_list))
for pdb in dna_pdb_list:
    os.system("rm -r "+pdb)

#!/usr/bin/env python
docstring='''
getRNAlist na_chain.list na_type.list na_chain.rna
    Filter na_chain.list by selecting only chains for RNA and RNA/DNA hybrid.
    Output the filtered list to na_chain.rna.
'''
import sys

if len(sys.argv)!=4:
    sys.stderr.write(docstring)
    exit()

infile=sys.argv[1]
molfile=sys.argv[2]
outfile=sys.argv[3]

fp=open(molfile,'r')
mol_dict=dict()
for line in fp.read().splitlines():
    chain,mol=line.split('\t')[:2]
    mol_dict[chain]=mol
fp.close()

fp=open(infile,'r')
txt=''
remove_list=[]
for chain in fp.read().splitlines():
    if chain in mol_dict and not mol_dict[chain].startswith("RNA"):
        remove_list.append(chain)
        continue
    txt+=chain+'\n'
fp.close()

fp=open(outfile,'w')
fp.write(txt)
fp.close()

print("remove %d DNA, DNA/RNA or DNA-RNA chains"%len(remove_list))

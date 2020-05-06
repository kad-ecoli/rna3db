#!/usr/bin/env python
docstring='''
getNAtype.py na_chain.list na_type.list
    For each chain listed in na_chain.list, find the chain in current folder,
    and output the molecule type (RNA/DNA) at na_type.list in the format:

    chain mol rna_count dna_count
      |    |     |         |------ number of DNA atoms
      |    |     |---------------- number of RNA atoms
      |    |---------------------- molecule (RNA, RNA/DNA, DNA-RNA, DNA/RNA, DNA)
      |--------------------------- chain ID

    Here, RNA and DNA refers to pure RNA and pure DNA, respectively.
    RNA/DNA are molecules with >= 90% RNA atoms
    DNA-RNA are molecules with <90% but >10% RNA atoms
    DNA/RNA are molecules with <=10% RNA atoms
    If ouput file na_type.list, new entries will be appended to the output.
'''
import sys, os
import gzip

if len(sys.argv)<3:
    sys.stderr.write(docstring)
    exit()

infile=sys.argv[1]
outfile=sys.argv[2]

parsed_chain_list=[]
if os.path.isfile(outfile):
    fp=open(outfile,'r')
    for line in fp.read().splitlines():
        parsed_chain_list.append(line.split()[0])
    fp.close()

fp=open(infile,'r')
chain_list=[]
for chain in fp.read().splitlines():
    if not chain in parsed_chain_list:
        chain_list.append(chain)
fp.close()

type_list=[]
txt=''
for chain in chain_list:
    dna_count=0
    rna_count=0
    filename=chain+".pdb"
    if not os.path.isfile(filename):
        filename=chain+".pdb.gz"
        if not os.path.isfile(filename):
            sys.stderr.write("ERROR! No such file %s/%s\n"%(
                os.getcwd(),filename))
            continue
        fp=gzip.open(filename,'r')
    else:
        fp=open(filename,'r')
    lines=fp.read().splitlines()
    fp.close()
    for line in lines:
        if line.startswith('END'):
            break
        if line.startswith('ATOM  '):
            resn=line[17:20]
            if resn[:2]=="  ":
                rna_count+=1
            elif resn[:2]==" D":
                dna_count+=1
    mol='NA'
    if   rna_count >0 and dna_count==0:
        mol="RNA"
    elif rna_count >0 and dna_count >0:
        if   rna_count>=9*dna_count:
            mol="RNA/DNA"
        elif dna_count>=9*rna_count:
            mol="DNA/RNA"
        else:
            mol="DNA-RNA"
    elif rna_count==0 and dna_count >0:
        mol="DNA"
    txt+="%s\t%s\t%d\t%d\n"%(chain,mol,rna_count,dna_count)
    if dna_count and rna_count:
        print("%s\t%s\t%d\t%d"%(chain,mol,rna_count,dna_count))

fp=open(outfile,'a')
fp.write(txt)
fp.close()

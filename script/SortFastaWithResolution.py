#!/usr/bin/env python
docstring=''' 
SortFastaWithResolution.py resolu.idx pdb_atom.fasta pdb_atom.sorted
    Sort fasta by length. If two sequence has the same length, take the one
    with better (lower) resolution. Output the files to 
    pdb_atom.sorted.{1.5,2.0,2.5,3.0,3.5,4.0,20.0,all}.fasta
'''

import sys
if len(sys.argv)!=4:
    sys.stderr.write(docstring)
    exit()

resolufile=sys.argv[1]
fastafile =sys.argv[2]
outfile   =sys.argv[3]

resolu_dict=dict()
fp=open(resolufile,'r')
lines=fp.read().splitlines()
fp.close()
nmr_dict=dict()
for line in lines:
    if not "\t;\t" in line:
        continue
    items=line.split("\t;\t")
    if len(items)!=2:
        continue
    idcode=items[0].lower()
    resolu=10
    nmr_dict[idcode]=True
    if len(items[1]):
        if float(items[1])>0:
            resolu=float(items[1])
            nmr_dict[idcode]=False
    resolu_dict[idcode]=resolu

fasta_list=[]
fp=open(fastafile,'r')
blocks=fp.read().strip().lstrip('>').split('\n>')
fp.close()
for block in blocks:
    lines=block.splitlines()
    header=lines[0]
    sequence=''.join(lines[1:])
    idcode=header[:4].lower()
    resolu=10
    if idcode in resolu_dict:
        resolu=resolu_dict[idcode]
    L=len(sequence)
    if L<10:
        continue
    fasta_list.append((resolu,L,header,sequence))

fasta_list=sorted(fasta_list)
fasta_list=[(-L,resolu,header,sequence) for resolu,L,header,sequence in fasta_list]
fasta_list=sorted(fasta_list)
fasta_list=[(-L,resolu,header,sequence) for L,resolu,header,sequence in fasta_list]
sequence_list=set()

txt=''
for L,resolu,header,sequence in fasta_list:
    if nmr_dict[header[:4]]:
        resolu="NA"
    txt+=">%s\t%s\t%d\n%s\n"%(header,str(resolu),L,sequence)
fp=open(outfile+".all.fasta",'w')
fp.write(txt)
fp.close()

for cutoff in [1.5,2.0,2.5,3.0,3.5,4.0,20]: # 20 includes NMR
    txt=''
    for L,resolu,header,sequence in fasta_list:
        if nmr_dict[header[:4]]:
            continue
        if resolu>cutoff:
            continue
        txt+=">%s\t%s\t%d\n%s\n"%(header,str(resolu),L,sequence)
    fp=open("%s.%.1f.fasta"%(outfile,cutoff),'w')
    fp.write(txt)
    fp.close()

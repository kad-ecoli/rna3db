#!/usr/bin/python
# fiber -rna -seq=ACGU ACGU.rna.pdb
# fiber      -seq=ACGT ACGT.dna.pdb 
rna="ACGU"
seq=""
for a in rna:
    for b in rna:
        seq+=a+b

for a in rna:
    for b in rna:
        for c in rna:
            seq+=a+b+c

print("fiber -rna -seq=%s %d.rna.pdb"%(seq,len(seq)))
seq=seq.replace('U','T')
print("fiber      -seq=%s %d.dna.pdb"%(seq,len(seq)))

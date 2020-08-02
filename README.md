# RNA3DB: maintaining a local copy of RNA 3D structure #

This package synchronize and curate the updated set of RNA structures from the Protein Data Bank (PDB) ftp site.

### Content ###
[1] ``script/`` scripts and executables for update.
- ``pdb/derived_data/pdb_seqres.txt.gz`` FASTA file for all chains (protein and nucleic acid)
- ``pdb/derived_data/na_seqres.txt`` FASTA file for all nucleic acid chains
- ``pdb/derived_data/na_type.list`` Table of molecule types for nucleic acid chains. This file is in the following format: chain, mol, rna_count, dna_count. rna_count and dna_count are the number of RNA and DNA atoms, respectively. mol is the molecule type decided based on %rna=(rna_count/(rna_count+dna_count)); mol has five values: RNA (%rna=100), RNA/DNA (90%<=%rna<100%), DNA-RNA (10%<%rna<90%), DNA/RNA (0%<%rna<=10%), DNA (%rna=0).
- ``pdb/derived_data/na_chain.list`` List RNA chains
- ``pdb/derived_data/index/resolu.idx`` Table of resolution

[2] ``pdb/derived_data/`` list and sequence of nucleic acid structures.

[3] ``pdb/data/structures/all/pdb/`` raw PDB files.

### Update ###

Structure database:
```bash
script/update.sh
```

Gene Ontology annotation database (optional):
```bash
script/GOA.sh
```

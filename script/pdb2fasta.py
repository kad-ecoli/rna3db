#!/usr/bin/env python
# 2016-04-27 Chengxin Zhang
docstring='''
pdb2fasta.py pdb.pdb > seq.fasta
    convert PDB file pdb.pdb to sequence FASTA file seq.fasta

options:
-PERMISSIVE={TER,ATOM,HETATM} how to treat nonstandatd amino acids
    TER   - (default) only HETATM in the middle of the chain and ATOM;
            discard HETATM at either end of the chain
    ATOM  - only allow ATOM residues
    HETATM- All all ATOM & HETATM residues with CA atom, including ligands

-allowX={false,true} whether to allow residue "X"
    true  - allow any amino acid type, including X
    false - (default) only allow amino acid that can mapped to standard
            residue types

-outfmt={PDB,COFACTOR} how to treat multichain PDB
    PDB      - convert multiple chains PDB to multiple sequence FASTA
    COFACTOR - convert PDB to single sequence FASTA. If there are 
               multiple chains in pdb, they will be treated as one chain.

-SEQRES={false,true} whether to convert from "SEQRES" entries
    false - (default) always convert from "ATOM" entries
    true  - convert from "SEQRES" entries if present, 
            otherwise convert from "ATOM" entries.
            (For PDB format input only)

-mol={rna,protein,rna,dna,all} which macromolecule type to use
    rna     - (default) use " C3'" atoms
    dna     - use " C3'" atoms, the same as rna
    all     - use " CA " or " C3'" atoms
    protein - use " CA " atoms

-dir={folder} -suffix={suffix}
    batch convert a folder of pdb files at {folder} with file name extension
    {suffix}. Example: 
    $ echo 157dA 157dB > list
    $ pdb2fasta -dir ./ list -suffix .pdb.gz
    This will convert ./157dA.pdb.gz and ./157dB.pdb.gz
'''
import sys,os
import shutil
import textwrap
import gzip,tarfile
import random

from MODRES_dicts import dna2rna,modres2na

def code_with_modified_residues():
    aa3to1=dict()
    for resName in modres2na:
        aa3to1[resName]=dna2rna[modres2na[resName]][-1].lower()
    return aa3to1

def pdbbundle2seq(tarball_name="pdb-bundle.tar.gz",PERMISSIVE="MSE",
    outfmt="PDB",allowX=True,mol="all"):
    '''convert best effort/minimum PDB bundle to sequence
    '''
    chain_id_mapping=dict()
    tarball_prefix=os.path.basename(tarball_name).split('.')[0]
    tar=tarfile.open(tarball_name,'r:gz')
    names=tar.getnames()
    PDBid=names[-1].split('-')[0]

    # parse chain id mapping
    fp=tar.extractfile(PDBid+"-chain-id-mapping.txt")
    map_txt=fp.read()
    fp.close()
    names=[] # list of PDB files in tarball
    for section in map_txt.split('\n'+PDBid+"-pdb-bundle"):
        if not ':' in section:
            continue
        idx,section=section.split(".pdb:\n")
        pdb_bundle_name=PDBid+"-pdb-bundle"+idx+".pdb"
        names.append(pdb_bundle_name)
        for line in section.splitlines():
            New_chain_ID,Original_chain_ID=line.split()
            key=pdb_bundle_name.split('.')[0]+':'+New_chain_ID
            value=tarball_prefix+':'+Original_chain_ID
            chain_id_mapping[key]=value

    header_list=[]
    sequence_list=[]
    for pdb_bundle_name in names:
        # parse text in *-pdb-bundle*.pdb
        fp=tar.extractfile(pdb_bundle_name)
        txt=fp.read()
        fp.close()
        header_list_tmp,sequence_list_tmp=pdbtxt2seq(
            txt,pdb_bundle_name,PERMISSIVE,outfmt,allowX,False,mol)
        if outfmt=="PDB":
            header_list+=[chain_id_mapping[h] for h \
                in header_list_tmp]
        sequence_list+=sequence_list_tmp
    if outfmt=="COFACTOR":
        sequence=''.join([''.join(s.splitlines()) for s in sequence_list])
        sequence=textwrap.fill(''.join(sequence),60)
        header=tarball_prefix+'\t'+str(len(sequence))
        header_list=[header]
        sequence_list=[sequence]
    return header_list,sequence_list

def pdb2seq(infile="pdb.pdb", PERMISSIVE="MSE", outfmt="PDB",
    allowX=True, SEQRES=False,mol="all"):
    '''Convert PDB to sequence.
    Return two lists, one for headers and the other for sequence.

    PERMISSIVE - whether allow non-standard residues
        ATOM:   Only allow ATOM residues
        HETATM: Allow all ATOM & HETATM residues, even if they are ligands
        MSE:   (default) Disallow any non-standard amino acid apart from MSE
    '''
    if infile.endswith(".tar.gz"): # best effort/minimum PDB bundle
        return pdbbundle2seq(infile,PERMISSIVE,outfmt,allowX,mol)
    elif infile.endswith(".gz"):
        fp=gzip.open(infile,'rU')
    else:
        fp=open(infile,'rU')
    txt=fp.read()
    fp.close()
    return pdbtxt2seq(txt,infile,PERMISSIVE,outfmt,allowX,SEQRES,mol)

def pdbtxt2seq(txt='',infile='pdb.pdb',PERMISSIVE="TER",outfmt="PDB",
    allowX=True,SEQRES=False,mol="all"):
    '''Convert PDB text "txt" to sequence read from PDB file "infile"
    Return two lists, one for headers and the other for sequence.

    PERMISSIVE - whether allow non-standard residues
        ATOM:   Only allow ATOM residues
        HETATM: Allow all ATOM & HETATM residues, even if they are ligands
        TER: (default) only HETATM in the middle of the chain and ATOM;
            discard HETATM at either end of the chain
    '''
    txt=txt.split("\nENDMDL")[0] # Only the first model
    if not "SEQRES" in txt:
        SEQRES=False # use "ATOM" is "SEQRES" is absent
    mol=mol.lower()

    aa3to1=code_with_modified_residues()
    
    chain_list=[]
    chain_dict=dict() # Each chain will be one key
    het_count_dict=dict()
    for line in txt.splitlines():
        line=line+' '*(80-len(line)) # Each line contains at least 80 char

        if SEQRES:
            if not line[:6]=="SEQRES":
                continue
            chain_id=line[11].replace(' ','_')
            tmp_seq=[]
            for residue in line[19:].split():
                if len(residue)!=3:
                    continue # only convert amino acid
                aa=aa3to1[residue] if residue in aa3to1 else 'X'
                if allowX or aa!='X':
                    tmp_seq.append((len(tmp_seq),aa))
            if tmp_seq:
                if not chain_id in chain_dict:
                    chain_dict[chain_id]=[]
                    chain_list.append(chain_id)
                chain_dict[chain_id]+=tmp_seq
            continue

        if (mol=="protein" and line[12:16]!=" CA ") or \
           (mol in {"rna","dna"} and line[12:16]!=" C3'") or \
           (mol=="all" and line[12:16]!=" CA " and line[12:16]!=" C3'"):
            continue
        if not line[16] in [' ','A']: # remove alternative location
            continue

        residue=line[17:20] # residue name

        if   PERMISSIVE == "ATOM"   and line[0:6]!="ATOM  ":
            continue
        elif PERMISSIVE in {"TER","HETATM"} and not line[0:6] in {"ATOM  ","HETATM"}:
            continue
        
        # underscore for empty chain identifier
        chain_id=line[21].replace(' ','_')
        res_num=int(line[22:26]) # residue sequence number
        aa='X'
        if residue in aa3to1:
            aa=aa3to1[residue]
        elif line[12:16]!=" CA ":
            aa='n'
        if not allowX and aa in ['X','n']:
            continue
        residue_tuple=(res_num,aa)

        if not chain_id in chain_dict:
            if line[0:6]=="HETATM" and PERMISSIVE=="TER":
                continue # remove C5' end HETATM
            chain_dict[chain_id]=[]
            chain_list.append(chain_id)

        if not residue_tuple in chain_dict[chain_id]:
            chain_dict[chain_id].append(residue_tuple)
            if PERMISSIVE=="TER":
                if not chain_id in het_count_dict:
                    het_count_dict[chain_id]=0
                if line[0:6]=="ATOM  ":
                    het_count_dict[chain_id]=0
                elif line[0:6]=="HETATM":
                    het_count_dict[chain_id]+=1
    
    if PERMISSIVE=="TER":
        for chain_id in het_count_dict:
            if het_count_dict[chain_id]>0:
                chain_dict[chain_id]=chain_dict[chain_id][:-het_count_dict[chain_id]]

    header_list=[]
    sequence_list=[]
    PDBID=os.path.basename(infile).split('.')[0]
    #PDBID=os.path.basename(infile).upper().split('.')[0]
    for chain_id in chain_list:
        res_num_list,sequence=zip(*chain_dict[chain_id])
        header_list.append(PDBID+':'+chain_id)
        sequence_list.append(''.join(sequence))

    if outfmt=="COFACTOR":
        if len(sequence_list)>1:
            sys.stderr.write("WARNING! Multichain PDB %s\n"%infile)
        sequence=''.join(sequence_list)
        header=PDBID+'\t'+str(len(sequence))
        sequence=textwrap.fill(''.join(sequence),60)
        header_list=[header]
        sequence_list=[sequence]
    return header_list,sequence_list


def pdb2fasta(infile="pdb.pdb", PERMISSIVE="MSE", outfmt="PDB",
    allowX=True,SEQRES=False, mol="all"):
    '''Convert PDB to FASTA'''
    header_list,sequence_list=pdb2seq(infile,PERMISSIVE,outfmt,allowX,SEQRES,mol)
    fasta_list=['>'+header_list[i]+'\n'+ \
                  sequence_list[i] for i in range(len(sequence_list))]
    return '\n'.join(fasta_list)+'\n'

if __name__=="__main__":
    PERMISSIVE="TER"
    outfmt="PDB"
    allowX=False
    SEQRES=False
    mol="rna"
    prefix=''
    suffix=''
    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-outfmt="):
            outfmt=arg[len("-outfmt="):].upper()
        elif arg.startswith("-PERMISSIVE="):
            PERMISSIVE=arg[len("-PERMISSIVE="):].upper()
        elif arg.startswith("-dir="):
            prefix=arg[len("-dir="):]
        elif arg.startswith("-suffix="):
            suffix=arg[len("-suffix="):]
        elif arg.startswith("-allowX="):
            allowX=(arg[len("-allowX="):].lower()=="true")
        elif arg.startswith("-SEQRES="):
            SEQRES=(arg[len("-SEQRES="):].lower()=="true")
        elif arg.startswith("-mol="):
            mol=arg[len("-mol="):].lower()
        elif arg.startswith("-"):
            sys.stderr.write("ERROR! Unknown argument %s\n"%arg)
            exit()
        else:
            argv.append(arg)

    if len(argv)<1:
        sys.stderr.write(docstring)
    
    for filename in argv:
        if prefix:
            fp=open(filename,'r')
            target_list=fp.read().splitlines()
            fp.close()
            for f in target_list:
                sys.stdout.write(pdb2fasta(
                    prefix+f+suffix, PERMISSIVE=PERMISSIVE, outfmt=outfmt,
                    allowX=allowX, SEQRES=SEQRES, mol=mol))
        else:
            sys.stdout.write(pdb2fasta(
                filename, PERMISSIVE=PERMISSIVE, outfmt=outfmt,
                allowX=allowX, SEQRES=SEQRES, mol=mol))

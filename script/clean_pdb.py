#!/usr/bin/env python
docstring='''
clean_pdb.py input.pdb clean.pdb
    Only keep one alternative location.
    Only keep first MODEL.

Option:
    -StartIndex={1,2,...}
        reindex residue number so that the first residue starts with
        StartIndex. default is no change.
    -NewChainID={_,A,...,Z,a,...,z,0,...,9, }
        Assign a new chain ID to output PDB. Default is no change.
    -KeepChainNum=1
        Number of chains to keep. If 0, keep all chains.
    -fixHET={0,1,2,3}
        How to deal with heteroatoms?
        0 - Keep all residues unchanged
        1 - Convert non-standard nucleotides to standard RNA nucleotides.
            Delete unconvertible nucleotides.
        2 - (Default) Only delete heteroatoms at 5' or 3' termini.
            Convert remaining non-standard nucleotides to standard
            RNA nucleotides. Delete unconvertible nucleotides.
        3 - delete all heteroatomas.
    -delNoBB={0,1,2,3}
        0 - Keep all partial residues
        1 - Delete all residues without P
        2 - Delete all residues without C3'
        3 - (Default) Delete all residues without P and C3'
    -removeH={0,1}
        0 - Keep all hydrogen atoms
        1 - (Default) delete all hydrogen atoms
    -delCol={0,1,2,3}
        whether to remove columns after xyz coordinates
        0 - keep all columns
        1 - remove element and charge
        2 - remove element, charge and Bfactor
        3 - (default) remove element, charge, Bfactor and occupancy

    -dir={folder} -suffix={suffix}
        batch convert a folder of pdb files at {folder} with file name extension
        {suffix}. Example: 
        $ echo 157dA 157dB > list
        $ clean_pdb.py -dir=./ list -suffix=.pdb.gz .pdb
        This will convert ./157dA.pdb.gz and ./157dB.pdb.gz to
        157dA.pdb and 157dB.pdb
'''
import sys,os
import re
import gzip
from MODRES_dicts import dna2rna,modres2na
standard_residues=set(dna2rna.values())

def parseNoBB(PDBtxt,delNoBB):
    '''
    delNoBB={0,1,2,3}
        0 - Keep all partial residues
        1 - Delete all residues without P
        2 - Delete all residues without C3'
        3 - (Default) Delete all residues without P and C3'
    '''
    blocks=[]
    hasP=False
    hasC3=False
    prev_resi=''
    block=''
    for line in PDBtxt.splitlines():
        if not line.startswith('ATOM  ') and not line.startswith('HETATM'):
            if len(block):
                if delNoBB==0 or \
                  (delNoBB==1 and hasP) or \
                  (delNoBB==2 and hasC3) or \
                  (delNoBB==3 and (hasP or hasC3)):
                    blocks.append(block)
            block=''
            hasP=False
            hasC3=False
            blocks.append(line+'\n')
            continue
        
        resi=line[21:26]
        if resi!=prev_resi and len(block):
            if delNoBB==0 or \
              (delNoBB==1 and hasP) or \
              (delNoBB==2 and hasC3) or \
              (delNoBB==3 and (hasP or hasC3)):
              blocks.append(block)
            block=''
            hasP=False
            hasC3=False
        prev_resi=resi
        atom=line[12:16]
        if atom==" P  ":
            hasP=True
        elif atom==" C3'":
            hasC3=True
        block+=line+'\n'

    return ''.join(blocks)


def reindexPDB(PDBtxt,startindex=1):
    '''
    reindex residue number of PDB format text
    
    options:
        startindex - index of first residue
        PDBtxt     - text of input PDB to be reindexed
    '''
    PDBtxt_reindex=''
    current_old_index='' # residue number in origin PDB
    prev_chainID=''

    for line in PDBtxt.splitlines():
        if len(line)<27 or (not line.startswith("ATOM  ") and \
            not line.startswith("HETATM") and not line.startswith("TER")):
            PDBtxt_reindex+=line+'\n'
            continue
        resSeq=line[22:27] # residue sequence number
        current_chainID=line[21] # chain identifier

        if prev_chainID!=current_chainID: # first residue encountered
            current_old_index=resSeq # residue number in origin PDB
            current_new_index=int(startindex)
            prev_chainID=current_chainID
            resSeq_new=str(current_new_index)
            resSeq_new=' '*(4-len(resSeq_new))+resSeq_new+' '
        elif resSeq!=current_old_index:
            current_new_index+=1
            current_old_index=resSeq
            resSeq_new=str(current_new_index)
            resSeq_new=' '*(4-len(resSeq_new))+resSeq_new+' '
        PDBtxt_reindex+=line[:16]+' '+line[17:22]+resSeq_new+line[27:]+'\n'
    return PDBtxt_reindex

def readPDB(infile,KeepChainNum,removeH,NewChainID=''):
    '''
    KeepChainNum=1
       Number of chains to keep. If 0, keep all chains.
    removeH={0,1}
       0 - Keep all hydrogen atoms
       1 - (Default) delete all hydrogen atoms
    NewChainID={A,...,Z,a,...,z,0,...,9, }
       Assign a new chain ID to output PDB. Default is no change.
    '''
    PDBtxt=''
    chainNum=0
    prev_chainID=''
    if infile.endswith(".gz"):
        fp=gzip.open(infile,'r')
    else:
        fp=open(infile,'r')
    for line in fp.read().splitlines():
        if line.startswith("END"):
            break
        elif line.startswith("TER"):
            PDBtxt+="TER\n"
        elif line.startswith("ATOM  ") or line.startswith("HETATM"):
            if len(line)<54:
                continue
            altLoc=line[16]
            if not altLoc in (' ','A'):
                continue
            elif altLoc=='A':
                line=line[:16]+' '+line[17:]
            if removeH:
                if len(line)>=78:
                    element=line[76:78].strip()
                    if element=='H':
                        continue
                else:
                    atomName=line[12:16].strip()
                    if atomName[0] in "1234567890":
                        atomName=atomName[1:]
                    if len(atomName) and atomName[0]=='H':
                        continue
            chainID=line[21]
            if prev_chainID!=chainID:
                chainNum+=1
                prev_chainID=chainID
                if KeepChainNum>0 and chainNum>KeepChainNum:
                    break
            if len(NewChainID):
                line=line[:21]+NewChainID[-1]+line[22:]
            PDBtxt+=line+'\n'
    fp.close()
    return PDBtxt

def parseHET(PDBtxt,fixHET=2):
    '''
    fixHET: How to deal with heteroatoms?
        0 - Keep all residues unchanged
        1 - Convert non-standard nucleotides to standard RNA nucleotides.
            Delete unconvertible nucleotides.
        2 - (Default) Only delete heteroatoms at 5' or 3' termini.
            Convert remaining non-standard nucleotides to standard
            RNA nucleotides. Delete unconvertible nucleotides.
        3 - delete all heteroatomas.
    '''
    PDBlines=[]
    prev_chainID=''
    resCount=0
    hetCount=0
    for line in PDBtxt.splitlines():
        if not line.startswith('ATOM  ') and not line.startswith('HETATM'):
            if fixHET==2 and hetCount:
                PDBlines=PDBlines[:-hetCount]
                hetCount=0
            PDBlines.append(line+'\n')
            continue
        chainID=line[21]
        resName=line[17:20]
        if fixHET==1:
            if not resName in modres2na:
                continue
            resName=dna2rna[modres2na[resName]]
            line="ATOM  "+line[6:17]+resName+line[20:]
        elif fixHET==2:
            if not resName in modres2na:
                continue
            resName=dna2rna[modres2na[resName]]
            if chainID!=prev_chainID:
                resCount=0
                if hetCount:
                    PDBlines=PDBlines[:-hetCount]
                    hetCount=0
            prev_chainID=chainID
            if resCount==0 and line.startswith('HETATM'):
                continue
            resCount+=1
            if line.startswith("ATOM  "):
                hetCount=0
            else:
                hetCount+=1
            line="ATOM  "+line[6:17]+resName+line[20:]
        elif fixHET==3:
            if not resName in standard_residues or not line.startswith('ATOM  '):
                continue
        PDBlines.append(line+'\n')


    if fixHET==2 and hetCount:
        PDBlines=PDBlines[:-hetCount]
    return ''.join(PDBlines)

def parsePDBcols(PDBtxt,delCol):
    '''
    -delCol={0,1,2,3}
        whether to remove columns after xyz coordinates
        0 - keep all columns
        1 - remove element and charge
        2 - remove element, charge and Bfactor
        3 - (default) remove element, charge, Bfactor and occupancy
    '''
    PDB_lines=[]
    for line in PDBtxt.splitlines():
        if line.startswith("ATOM  ") or line.startswith("HETATM"):
            if delCol==1:
                line=line[:66]
            elif delCol==2:
                line=line[:60]
            elif delCol==3:
                line=line[:54]
        PDB_lines.append(line+'\n')
    return ''.join(PDB_lines)

if __name__=="__main__":
    #### parse commandline arguments ####
    StartIndex="NA"
    NewChainID=''
    KeepChainNum=1
    fixHET=2
    delNoBB=3
    removeH=True
    prefix=''
    suffix=''
    delCol=3
    
    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-StartIndex="):
            StartIndex=int(arg[len("-StartIndex="):])
        elif arg.startswith("-NewChainID="):
            NewChainID=arg[len("-NewChainID="):][-1:].replace('_',' ')
        elif arg.startswith("-KeepChainNum="):
            KeepChainNum=int(arg[len("-KeepChainNum="):])
        elif arg.startswith("-fixHET="):
            fixHET=int(arg[len("-fixHET="):])
        elif arg.startswith("-delNoBB="):
            delNoBB=int(arg[len("-delNoBB="):])
        elif arg.startswith("-removeH="):
            removeH=int(arg[len("-removeH="):])
        elif arg.startswith("-delCol="):
            delCol=int(arg[len("-delCol="):])
        elif arg.startswith("-dir="):
            prefix=arg[len("-dir="):]
        elif arg.startswith("-suffix="):
            suffix=arg[len("-suffix="):]
        elif arg.startswith('-'):
            sys.stderr.write("ERROR! Unknown option %s\n"%arg)
            exit()
        else:
            argv.append(arg)

    if len(argv)!=2:
        sys.stderr.write(docstring)
        exit()
    infile=argv[0]
    outfile=argv[1]

    if prefix:
        fp=open(infile,'r')
        infile_list=[]
        outfile_list=[]
        for line in fp.read().splitlines():
            target=line.split()[0]
            filename=os.path.join(prefix,target+suffix)
            if not os.path.isfile(filename):
                filename=os.path.join(prefix,target[:4],target+suffix)
            if not os.path.isfile(filename):
                sys.stderr.write("ERROR! No such file %s/{%s/}%s%s\n"%(
                    prefix,target[:4],target,suffix))
                continue
            infile_list.append(filename)
            outfile_list.append(os.path.join(os.getcwd(),target+outfile))
        fp.close()
    else:
        infile_list=[infile]
        outfile_list=[outfile]

    for infile,outfile in zip(infile_list,outfile_list):
        print("%s => %s"%(infile,outfile))
        #### parse PDB file ####
        PDBtxt=readPDB(infile,KeepChainNum,removeH,NewChainID)
        if fixHET>0:
            PDBtxt=parseHET(PDBtxt,fixHET)
        if delNoBB:
            PDBtxt=parseNoBB(PDBtxt,delNoBB)
        if StartIndex!="NA":
            PDBtxt=reindexPDB(PDBtxt,StartIndex)
        if delCol:
            PDBtxt=parsePDBcols(PDBtxt,delCol)

        #### write PDB file ####
        if outfile=='-':
            sys.stdout.write(PDBtxt)
        else:
            if outfile.endswith(".gz"):
                fp=gzip.open(outfile,'w')
            else:
                fp=open(outfile,'w')
            fp.write(PDBtxt)
            fp.close()

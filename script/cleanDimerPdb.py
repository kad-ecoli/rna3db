#!/usr/bin/env python
docstring='''
cleanDimerPdb.py input.pdb clean.pdb
    Only keep one alternative location.
    Only keep first MODEL.
    Sequentially assign new chain ID of A,...,Z,a...z,0...9 to output PDB
    Sequentially renumber residue index from 1 to L

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
from clean_pdb import readPDB, parseHET, parseNoBB, reindexPDB, parsePDBcols

def getLastResidueIndex(PDBtxt):
    for line in PDBtxt.splitlines()[::-1]:
        if line.startswith("ATOM  "):
            return int(line[22:26])
    return 0

if __name__=="__main__":
    #### parse commandline arguments ####
    KeepChainNum=0
    fixHET=2
    delNoBB=3
    removeH=True
    prefix=''
    suffix=''
    delCol=3
    
    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-fixHET="):
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
    inListFile=argv[0]
    outSuffix =argv[1]

    fp=open(inListFile,'r')
    infile_list=[]
    outfile_list=[]
    for line in fp.read().splitlines():
        line=line.split()[0]
        pdbID=line.split(':')[0]
        chainID_list=line.split(':')[1:]
        filename_list=[]
        for chainID in chainID_list:
            filename=os.path.join(prefix,pdbID+chainID+suffix)
            if not os.path.isfile(filename):
                filename=os.path.join(prefix,pdbID,pdbID+chainID+suffix)
            if not os.path.isfile(filename):
                sys.stderr.write("ERROR! No such file %s/{%s/}%s%s\n"%(
                    prefix,pdbID,pdbID+chainID,suffix))
                filename_list=[]
                break
            filename_list.append(filename)
        if len(filename_list)==0:
            sys.stderr.write("ERROR! Cannot find input file for %s\n"%(line))
            continue
        infile_list.append(list(filename_list))
        outfile_list.append(os.path.join(os.getcwd(),line+outSuffix))
    fp.close()

    NewChainIDlist="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890"
    for filename_list,outfile in zip(infile_list,outfile_list):
        print("%s => %s"%(' + '.join(filename_list),outfile))
        #### parse PDB file ####
        StartIndex=1
        txt=''
        for c,filename in enumerate(filename_list):
            PDBtxt=readPDB(filename,KeepChainNum,removeH,NewChainIDlist[c])
            if fixHET>0:
                PDBtxt=parseHET(PDBtxt,fixHET)
            if delNoBB:
                PDBtxt=parseNoBB(PDBtxt,delNoBB)
            PDBtxt=reindexPDB(PDBtxt,StartIndex)
            if delCol:
                PDBtxt=parsePDBcols(PDBtxt,delCol)
            StartIndex=getLastResidueIndex(PDBtxt)+1
            txt+=PDBtxt

        #### write PDB file ####
        if outfile.endswith(".gz"):
            fp=gzip.open(outfile,'w')
        else:
            fp=open(outfile,'w')
        fp.write(txt)
        fp.close()

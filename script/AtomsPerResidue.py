#!/usr/bin/env python
docstring='''
AtomsPerResidue.py input.pdb
    Report average number of atoms per residue in input.pdb

Options
    -dir={folder} -suffix={suffix}
        batch calculate a folder of pdb files at {folder} with file name
        extension {suffix}. Example: 
        $ echo 157dA 157dB > list
        $ AtomsPerResidue.py -dir=./ list -suffix=.pdb.gz
'''
import sys, os

def AtomsPerResidue(infile):
    AtomNum=0
    ResidueList=[]
    fp=open(infile,'r')
    for line in fp.read().splitlines():
        if not line.startswith("ATOM  ") and not line.startswith("HETATM"):
            continue
        AtomNum+=1
        ResidueList.append(line[21:27])
    fp.close()
    return AtomNum,len(set(ResidueList))

if __name__=="__main__":
    prefix=''
    suffix=''
    
    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-dir="):
            prefix=arg[len("-dir="):]
        elif arg.startswith("-suffix="):
            suffix=arg[len("-suffix="):]
        elif arg.startswith('-'):
            sys.stderr.write("ERROR! Unknown option %s\n"%arg)
            exit()
        else:
            argv.append(arg)

    if len(argv)!=1:
        sys.stderr.write(docstring)
        exit()

    infile=argv[0]
    if prefix:
        fp=open(infile,'r')
        infile_list=[]
        for line in fp.read().splitlines():
            target=line.split()[0]
            infile_list.append(os.path.join(prefix,target+suffix))
        fp.close()
    else:
        infile_list=[infile]

    for infile in infile_list:
        AtomNum,ResidueNum=AtomsPerResidue(infile)
        if len(prefix):
            infile=infile[(len(prefix)+(prefix[-1]!='/')):]
        if len(suffix):
            infile=infile[:-len(suffix)]
        print("%s\t%.2f"%(infile,1.*AtomNum/ResidueNum))

#!/usr/bin/env python
docstring='''
makeCullTable.py cull/
    read data from cull/pdb_atom.*fasta
    draw tables at cull/README.md
'''


from datetime import datetime
import sys, os

def CountSeqFasta(filename):
    ''' count number of sequence in fasta file '''
    fp=open(filename,'r')
    seqnum=0
    for line in fp.read().splitlines():
        seqnum+=line.startswith('>')
    fp.close()
    return seqnum

if __name__=="__main__":
    if len(sys.argv)!=2:
        sys.stderr.write(docstring)
        exit()
    datadir=os.path.abspath(sys.argv[1])

    README_md='''# Representative Sets of RNA 3D Structures #
Last update: %s
RNA chains with at least 10 nucleotides, including RNA chains with up to 10%% DNA nucleotides.
Representative sets:
| Resolution |     all     | -c 1. -s 1. | -c .9 -s .9 | -c .8 -s .8 |   -c 1.0    |   -c 0.9    |   -c 0.8    |
|    :---:   |    :---:    |    :---:    |     :---:   |     :---:   |    :---:    |    :---:    |    :---:    |
'''%(datetime.now()) # $DATE $TOTAL


    for resolu in ["1.5","2.0","2.5","3.0","3.5","4.0","20.0","all"]:
        README_md+=("|     %s    |"%resolu).replace(" 20.0 ","20.0 ")
        seqnum=str(CountSeqFasta("%s/pdb_atom.sort.%s.fasta"%(datadir,resolu)))
        README_md+=' '*(12-len(seqnum))+seqnum+' |'
        for s in [True,False]:
            for c in ["1.0","0.9","0.8"]:
                seqnum=str(CountSeqFasta("%s/pdb_atom.sort.%s_c%s_s%s"%(
                    datadir,resolu,c,c if s else "0.0")))
                README_md+=' '*(12-len(seqnum))+seqnum+' |'
        README_md+='\n'

fp=open(datadir+"/README.md",'w')
fp.write(README_md)
fp.close()

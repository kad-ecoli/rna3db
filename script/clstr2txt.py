#!/usr/bin/env python
docstring='''
clstr2txt input.clstr output.txt
    convert CD-HIT format cluster file to text file, where each line is one cluster;
    members in the same cluster are separated by tab
'''
import sys

def clstr2txt(infile,outfile):
    txt=''
    fp=open(infile,'r')
    for block in ('\n'+fp.read()).split('\n>Cluster ')[1:]:
        member_list=[]
        representative=''
        for line in block.splitlines()[1:]:
            items=line.split(', >')[1].split()
            target=items[0][:-3]
            marker=items[1]
            if marker=='*':
                representative=target
            else:
                member_list.append(target)
        txt+='\t'.join([representative]+member_list)+'\n'
    fp.close()
    fp=open(outfile,'w')
    fp.write(txt)
    fp.close()
    return

if __name__=="__main__":
    if len(sys.argv)!=3:
        sys.stderr.write(docstring)
        exit()

    infile=sys.argv[1]
    outfile=sys.argv[2]
    clstr2txt(infile,outfile)

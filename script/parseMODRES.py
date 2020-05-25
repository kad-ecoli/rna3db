#!/usr/bin/env python
docstring='''
parseMODRES.py MODRES MODRES_dicts.py
    parse MODRES and output it as python dictionary to MODRES_dicts.py
'''

import sys

def parseMODRES(infile):
    txt='''#!/usr/bin/env python
dna2rna={        ' DT':'  U',
    ' DA':'  A', ' DU':'  U', ' DC':'  C', ' DG':'  G', ' DI':'  A',
    '  A':'  A', '  U':'  U', '  C':'  C', '  G':'  G', '  I':'  A',
}

modres2na = {    ' DT':' DT',
    ' DA':' DA', ' DU':' DU', ' DC':' DC', ' DG':' DG', ' DI':' DI',
    '  A':'  A', '  U':'  U', '  C':'  C', '  G':'  G', '  I':'  I',
'''
    
    na_set={   ' DT', 
        ' DA', ' DU', ' DC', ' DG', ' DI',
        '  A', '  U', '  C', '  G', '  I',
    }

    fp=open(infile,'r')
    lines=fp.read().splitlines()
    fp.close()

    modres_dict=dict()
    for line in lines:
        if not line.startswith("MODRES"):
            continue
        resName=line[12:15]
        stdRes =line[24:27]
        if resName in na_set:
            continue
        if not stdRes in na_set:
            continue
        if not resName in modres_dict:
            modres_dict[resName]=dict()
        if not stdRes in modres_dict[resName]:
            modres_dict[resName][stdRes]=0
        modres_dict[resName][stdRes]+=1

    count=-1
    line="    "
    for resName in modres_dict:
        repeat,stdRes=sorted([(modres_dict[resName][stdRes],stdRes
            ) for stdRes in modres_dict[resName]],reverse=True)[0]
        if len(modres_dict[resName])>1:
            print(resName,modres_dict[resName])
        count+=1
        if count==5:
            txt+=line.rstrip()+'\n'
            line="    "
            count=0
        line+="'%s':'%s', "%(resName,stdRes)
    
    if line!="    ":
        txt+=line.rstrip()
    txt+="\n}\n"
    return txt

if __name__=="__main__":
    if len(sys.argv)!=3:
        sys.stderr.write(docstring)
        exit()

    infile =sys.argv[1]
    outfile=sys.argv[2]
    fp=open(outfile,'w')
    fp.write(parseMODRES(infile))
    fp.close()

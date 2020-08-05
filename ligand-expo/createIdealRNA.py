#!/usr/bin/python
''' create IdealRNA.hpp '''
from string import Template

atom_template=Template('        ${RESN}["$ATOM"]=tmp; ${RESN}["$ATOM"][0]=$X; ${RESN}["$ATOM"][1]=$Y; ${RESN}["$ATOM"][2]=$Z;\n')

hpp_txt='''/* IdealRNA.hpp - Idealized RNA nucleotide conformations */
#include <string>
#include <map>
#include <vector>

using namespace std;

/* option: 0 - ideal pdb; 1 - model pdb; 2 - 3dna fiber */
map<string, map<string,vector<float> > >parse_ideal_pdb(int option=2)
{
    vector<float> tmp(3,0);
    map<string, map<string,vector<float> > >ideal_pdb;
    if (option==0)
    {
'''

for resn in ['  A','  U','  C','  G',
             ' DA',' DT',' DC',' DG']:
    fp=open("%s/%s_ideal.pdb"%(resn.strip(),resn.strip()),'r')
    lines=fp.read().splitlines()
    fp.close()

    hpp_txt+="        map<string,vector<float> >%s;\n"%(resn.strip())
    for line in lines:
        if not line.startswith("ATOM  ") or line[77]=='H' or line[12:16]==" OP3":
            continue
        atom=line[12:16]
        x   =line[30:38]
        y   =line[38:46]
        z   =line[46:54]
        hpp_txt+=atom_template.substitute(
            RESN=resn.strip(),
            ATOM=atom,
            X=x,
            Y=y,
            Z=z,
        )
    hpp_txt+='''        ideal_pdb["%s"]=%s;
        map<string,vector<float> >().swap(%s);

'''%(resn,resn.strip(),resn.strip())

hpp_txt+='''    }
    else if (option==1)
    {
'''

for resn in ['  A','  U','  C','  G',
             ' DA',' DT',' DC',' DG']:
    fp=open("%s/%s_model.pdb"%(resn.strip(),resn.strip()),'r')
    lines=fp.read().splitlines()
    fp.close()

    hpp_txt+="        map<string,vector<float> >%s;\n"%(resn.strip())
    for line in lines:
        if not line.startswith("ATOM  ") or line[77]=='H' or line[12:16]==" OP3":
            continue
        atom=line[12:16]
        x   =line[30:38]
        y   =line[38:46]
        z   =line[46:54]
        hpp_txt+=atom_template.substitute(
            RESN=resn.strip(),
            ATOM=atom,
            X=x,
            Y=y,
            Z=z,
        )
    hpp_txt+='''        ideal_pdb["%s"]=%s;
        map<string,vector<float> >().swap(%s);

'''%(resn,resn.strip(),resn.strip())

hpp_txt+='''    }
    else if (option==2)
    {
'''
# generated using 
# fiber -rna -seq=ACGU ACGU.rna.pdb
# fiber      -seq=ACGT ACGT.dna.pdb 
fp=open("ACGU.rna.pdb",'r')
lines=fp.read().splitlines()
fp.close()
fp=open("ACGT.dna.pdb",'r')
lines+=fp.read().splitlines()
fp.close()

for resn in ['  A','  U','  C','  G',
             ' DA',' DT',' DC',' DG']:
    hpp_txt+="        map<string,vector<float> >%s;\n"%(resn.strip())
    for line in lines:
        if not line.startswith("ATOM  ") or line[77]=='H' or \
            line[12:16]==" OP3" or line[21]!='A' or line[17:20]!=resn:
            continue
        atom=line[12:16]
        x   =line[30:38]
        y   =line[38:46]
        z   =line[46:54]
        hpp_txt+=atom_template.substitute(
            RESN=resn.strip(),
            ATOM=atom,
            X=x,
            Y=y,
            Z=z,
        )
    hpp_txt+='''        ideal_pdb["%s"]=%s;
        map<string,vector<float> >().swap(%s);

'''%(resn,resn.strip(),resn.strip())

hpp_txt+='''    }
    vector<float> ().swap(tmp);
    return ideal_pdb;
}
'''

fp=open("IdealRNA.hpp",'w')
fp.write(hpp_txt)
fp.close()

#!/usr/bin/python
''' create IdealDNA.hpp '''
from string import Template

atom_template=Template('        ${RESN}["$ATOM"]=tmp; ${RESN}["$ATOM"][0]=$X; ${RESN}["$ATOM"][1]=$Y; ${RESN}["$ATOM"][2]=$Z;\n')
atom2_template=Template('    ${RESN}["$ATOM"]=tmp; ${RESN}["$ATOM"][0]=$X; ${RESN}["$ATOM"][1]=$Y; ${RESN}["$ATOM"][2]=$Z;\n')

hpp_txt='''/* IdealDNA.hpp - Idealized DNA nucleotide conformations */
#include <string>
#include <map>
#include <vector>

using namespace std;

/* option: 0 - ideal pdb; 1 - model pdb; 2 - 3dna fiber */
map<string, map<string,vector<float> > >parse_ideal_dna(int option=2)
{
    vector<float> tmp(3,0);
    map<string, map<string,vector<float> > >ideal_dna;
    if (option==0)
    {
'''

for resn in [' DA',' DC',' DG',' DT']:
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
    hpp_txt+='''        ideal_dna["%s"]=%s;
        map<string,vector<float> >().swap(%s);

'''%(resn,resn.strip(),resn.strip())

hpp_txt+='''    }
    else if (option==1)
    {
'''

for resn in [' DA',' DC',' DG',' DT']:
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
    hpp_txt+='''        ideal_dna["%s"]=%s;
        map<string,vector<float> >().swap(%s);

'''%(resn,resn.strip(),resn.strip())

hpp_txt+='''    }
    else if (option==2)
    {
'''
# generated using 
# fiber      -seq=ACGT ACGT.dna.pdb 
fp=open("ACGT.dna.pdb",'r')
lines+=fp.read().splitlines()
fp.close()

for resn in [' DA',' DC',' DG',' DT']:
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
    hpp_txt+='''        ideal_dna["%s"]=%s;
        map<string,vector<float> >().swap(%s);

'''%(resn,resn.strip(),resn.strip())

hpp_txt+='''    }
'''

# generated using
# fiber      -seq=AAACAGATCACCCGCTGAGCGGGTTATCTGTTAAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATGATTCAACACCAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTCCTGCTTGAAGACGAGGATGCAGCCGCGGCTGGAGGCGGGGGTGTAGTCGTGGTTTAATACTAGTATTCATCCTCGTCTTGATGCTGGTGTTTATTCTTGTTT 224.dna.pdb

fp=open("224.dna.pdb",'r')
lines=fp.read().splitlines()
fp.close()
r=0
for resn1 in [' DA',' DC',' DG',' DT']:
    for resn2 in [' DA',' DC',' DG',' DT']:
        key=resn1.strip()+resn2.strip()+'0'
        hpp_txt+="    map<string,vector<float> >%s;\n"%key
        r+=1
        resi=str(r)
        resi=' '*(4-len(resi))+resi
        for line in lines:
            if not line.startswith("ATOM  ") or line[77]=='H' or \
                line[12:16]==" OP3" or line[21]!='A' or line[22:26]!=resi or \
                line[17:20]!=resn1:
                continue
            atom=line[12:16]
            x   =line[30:38]
            y   =line[38:46]
            z   =line[46:54]
            hpp_txt+=atom2_template.substitute(
                RESN=key,
                ATOM=atom,
                X=x,
                Y=y,
                Z=z,
            )
        hpp_txt+='''    ideal_dna["%s"]=%s;
    map<string,vector<float> >().swap(%s);

'''%(resn1+resn2+'0',key,key)

        key=resn1.strip()+resn2.strip()+'1'
        hpp_txt+="    map<string,vector<float> >%s;\n"%key
        r+=1
        resi=str(r)
        resi=' '*(4-len(resi))+resi
        for line in lines:
            if not line.startswith("ATOM  ") or line[77]=='H' or \
                line[12:16]==" OP3" or line[21]!='A' or line[22:26]!=resi or \
                line[17:20]!=resn2:
                continue
            atom=line[12:16]
            x   =line[30:38]
            y   =line[38:46]
            z   =line[46:54]
            hpp_txt+=atom2_template.substitute(
                RESN=key,
                ATOM=atom,
                X=x,
                Y=y,
                Z=z,
            )
        hpp_txt+='''    ideal_dna["%s"]=%s;
    map<string,vector<float> >().swap(%s);

'''%(resn1+resn2+'1',key,key)

for resn1 in [' DA',' DC',' DG',' DT']:
    for resn2 in [' DA',' DC',' DG',' DT']:
        for resn3 in [' DA',' DC',' DG',' DT']:
            key=resn1.strip()+resn2.strip()+resn3.strip()+'0'
            hpp_txt+="    map<string,vector<float> >%s;\n"%key
            r+=1
            resi=str(r)
            resi=' '*(4-len(resi))+resi
            for line in lines:
                if not line.startswith("ATOM  ") or line[77]=='H' or \
                    line[12:16]==" OP3" or line[21]!='A' or line[22:26]!=resi or \
                    line[17:20]!=resn1:
                    continue
                atom=line[12:16]
                x   =line[30:38]
                y   =line[38:46]
                z   =line[46:54]
                hpp_txt+=atom2_template.substitute(
                    RESN=key,
                    ATOM=atom,
                    X=x,
                    Y=y,
                    Z=z,
                )
            hpp_txt+='''    ideal_dna["%s"]=%s;
    map<string,vector<float> >().swap(%s);

'''%(resn1+resn2+resn3+'0',key,key)

            key=resn1.strip()+resn2.strip()+resn3.strip()+'1'
            hpp_txt+="    map<string,vector<float> >%s;\n"%key
            r+=1
            resi=str(r)
            resi=' '*(4-len(resi))+resi
            for line in lines:
                if not line.startswith("ATOM  ") or line[77]=='H' or \
                    line[12:16]==" OP3" or line[21]!='A' or line[22:26]!=resi or \
                    line[17:20]!=resn2:
                    continue
                atom=line[12:16]
                x   =line[30:38]
                y   =line[38:46]
                z   =line[46:54]
                hpp_txt+=atom2_template.substitute(
                    RESN=key,
                    ATOM=atom,
                    X=x,
                    Y=y,
                    Z=z,
                )
            hpp_txt+='''    ideal_dna["%s"]=%s;
    map<string,vector<float> >().swap(%s);

'''%(resn1+resn2+resn3+'1',key,key)

            key=resn1.strip()+resn2.strip()+resn3.strip()+'2'
            hpp_txt+="    map<string,vector<float> >%s;\n"%key
            r+=1
            resi=str(r)
            resi=' '*(4-len(resi))+resi
            for line in lines:
                if not line.startswith("ATOM  ") or line[77]=='H' or \
                    line[12:16]==" OP3" or line[21]!='A' or line[22:26]!=resi or \
                    line[17:20]!=resn3:
                    continue
                atom=line[12:16]
                x   =line[30:38]
                y   =line[38:46]
                z   =line[46:54]
                hpp_txt+=atom2_template.substitute(
                    RESN=key,
                    ATOM=atom,
                    X=x,
                    Y=y,
                    Z=z,
                )
            hpp_txt+='''    ideal_dna["%s"]=%s;
    map<string,vector<float> >().swap(%s);

'''%(resn1+resn2+resn3+'2',key,key)

hpp_txt+='''    vector<float> ().swap(tmp);
    return ideal_dna;
}
'''

fp=open("IdealDNA.hpp",'w')
fp.write(hpp_txt)
fp.close()

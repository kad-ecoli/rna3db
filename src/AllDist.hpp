/* calculate backbone torsion angles (omega, psi, phi, kappa, alpha) */
#ifndef AllDist_HPP
#define AllDist_HPP 1
#include <cstring>
#include "PDBParser.hpp"
#include "GeometryTools.hpp"
#include "StringTools.hpp"

using namespace std;

size_t AllDist(vector<char>&chain1_list, vector<char>&chain2_list,
    vector<string>&resn1_list, vector<string>&resn2_list,
    vector<int>&resi1_list,    vector<int>&resi2_list,
    vector<char>&icode1_list,  vector<char>&icode2_list,
    const ModelUnit &pdb_entry, vector<vector<float> >&BPlenMat)
{
    vector<float> tmp_len(10,-1.);

    /* coordinates of current residue */
    vector<float>Pi(3,0);  bool has_Pi;  vector<float>Pj(3,0);  bool has_Pj; 
    vector<float>O5i(3,0); bool has_O5i; vector<float>O5j(3,0); bool has_O5j;
    vector<float>C5i(3,0); bool has_C5i; vector<float>C5j(3,0); bool has_C5j;
    vector<float>C4i(3,0); bool has_C4i; vector<float>C4j(3,0); bool has_C4j;
    vector<float>C3i(3,0); bool has_C3i; vector<float>C3j(3,0); bool has_C3j;
    vector<float>C2i(3,0); bool has_C2i; vector<float>C2j(3,0); bool has_C2j;
    vector<float>C1i(3,0); bool has_C1i; vector<float>C1j(3,0); bool has_C1j;
    vector<float>O4i(3,0); bool has_O4i; vector<float>O4j(3,0); bool has_O4j;
    vector<float>O3i(3,0); bool has_O3i; vector<float>O3j(3,0); bool has_O3j;
    vector<float>Nxi(3,0); bool has_Nxi; vector<float>Nxj(3,0); bool has_Nxj;

    size_t bp, c1, c2, r1, r2, a1, a2;
    char base;
    bp=0;
    for (c1=0;c1<pdb_entry.chains.size(); c1++)
    {
        for (r1=0; r1<pdb_entry.chains[c1].residues.size(); r1++)
        {
            for (c2=c1; c2<pdb_entry.chains.size(); c2++)
            {
                for (r2=(c1==c2)*(r1+1); r2<pdb_entry.chains[c2].residues.size(); r2++)
                {
                    chain1_list.push_back(pdb_entry.chains[c1].chainID);
                    chain2_list.push_back(pdb_entry.chains[c2].chainID);
                    resn1_list.push_back(pdb_entry.chains[c1].residues[r1].resn);
                    resn2_list.push_back(pdb_entry.chains[c2].residues[r2].resn);
                    resi1_list.push_back(pdb_entry.chains[c1].residues[r1].resi);
                    resi2_list.push_back(pdb_entry.chains[c2].residues[r2].resi);
                    icode1_list.push_back(pdb_entry.chains[c1].residues[r1].icode);
                    icode2_list.push_back(pdb_entry.chains[c2].residues[r2].icode);
                    BPlenMat.push_back(tmp_len);
                    
                    has_Pi  = has_Pj  = false;
                    has_O5i = has_O5j = false;
                    has_C5i = has_C5j = false;
                    has_C4i = has_C4j = false;
                    has_C3i = has_C3j = false;
                    has_C2i = has_C2j = false;
                    has_C1i = has_C1j = false;
                    has_O4i = has_O4j = false;
                    has_O3i = has_O3j = false;
                    has_Nxi = has_Nxj = false;

                    base=tolower(pdb_entry.chains[c1].residues[r1].resn[2]);
                    for (a1=0;a1<pdb_entry.chains[c1].residues[r1].atoms.size();a1++)
                    {
                        if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" P  ")
                        {
                            has_Pi=true;
                            Pi=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                        }
                        else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" O5'")
                        {
                            has_O5i=true;
                            O5i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                        }
                        else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" C5'")
                        {
                            has_C5i=true;
                            C5i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                        }
                        else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" C4'")
                        {
                            has_C4i=true;
                            C4i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                        }
                        else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" C3'")
                        {
                            has_C3i=true;
                            C3i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                        }
                        else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" C2'")
                        {
                            has_C2i=true;
                            C2i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                        }
                        else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" C1'")
                        {
                            has_C1i=true;
                            C1i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                        }
                        else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" O4'")
                        {
                            has_O4i=true;
                            O4i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                        }
                        else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" O3'")
                        {
                            has_O3i=true;
                            O3i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                        }
                        else if ((pdb_entry.chains[c1].residues[r1].atoms[a1].name==" N9 "
                                 && (base=='a' || base=='g')) ||
                                (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" N1 "
                                 && (base=='c' || base=='t' || base=='u')))
                        {
                            has_Nxi=true;
                            Nxi=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                        }
                    }
                    
                    base=tolower(pdb_entry.chains[c2].residues[r2].resn[2]);
                    for (a2=0;a2<pdb_entry.chains[c2].residues[r2].atoms.size();a2++)
                    {
                        if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" P  ")
                        {
                            has_Pj=true;
                            Pj=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" O5'")
                        {
                            has_O5j=true;
                            O5j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" C5'")
                        {
                            has_C5j=true;
                            C5j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" C4'")
                        {
                            has_C4j=true;
                            C4j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" C3'")
                        {
                            has_C3j=true;
                            C3j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" C2'")
                        {
                            has_C2j=true;
                            C2j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" C1'")
                        {
                            has_C1j=true;
                            C1j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" O4'")
                        {
                            has_O4j=true;
                            O4j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" O3'")
                        {
                            has_O3j=true;
                            O3j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if ((pdb_entry.chains[c2].residues[r2].atoms[a2].name==" N9 "
                                    && (base=='a' || base=='g')) ||
                                (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" N1 "
                                 && (base=='c' || base=='t' || base=='u')))
                        {
                            has_Nxj=true;
                            Nxj=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                    }


                    if (has_Pi  && has_Pj)  BPlenMat[bp][0]=Points2Distance( Pi, Pj); //   PP: P[i]-P[j]
                    if (has_O5i && has_O5j) BPlenMat[bp][1]=Points2Distance(O5i,O5j); // O5O5: O5'[i]-O5'[j]
                    if (has_C5i && has_C5j) BPlenMat[bp][2]=Points2Distance(C5i,C5j); // C5C5: C5'[i]-C5'[j]
                    if (has_C4i && has_C4j) BPlenMat[bp][3]=Points2Distance(C4i,C4j); // C4C4: C4'[i]-C4'[j]
                    if (has_C3i && has_C3j) BPlenMat[bp][4]=Points2Distance(C3i,C3j); // C3C3: C3'[i]-C3'[j]
                    if (has_C2i && has_C2j) BPlenMat[bp][5]=Points2Distance(C2i,C2j); // C2C2: C2'[i]-C2'[j]
                    if (has_C1i && has_C1j) BPlenMat[bp][6]=Points2Distance(C1i,C1j); // C1C1: C1'[i]-C1'[j]
                    if (has_O4i && has_O4j) BPlenMat[bp][7]=Points2Distance(O4i,O4j); // O4O4: O4'[i]-O4'[j]
                    if (has_O3i && has_O3j) BPlenMat[bp][8]=Points2Distance(O3i,O3j); // O3O3: O3'[i]-O3'[j]
                    if (has_Nxi && has_Nxj) BPlenMat[bp][9]=Points2Distance(Nxi,Nxj); //   NN: N[i]-N[j]
                    bp++;
                }
            }
        }
    }

    /* clean up */
    tmp_len.clear();

    Pi.clear();  Pj.clear();
    O5i.clear(); O5j.clear();
    C5i.clear(); C5j.clear();
    C4i.clear(); C4j.clear();
    C3i.clear(); C3j.clear();
    C2i.clear(); C2j.clear();
    C1i.clear(); C1j.clear();
    O4i.clear(); O4j.clear();
    O3i.clear(); O3j.clear();
    Nxi.clear(); Nxj.clear();
    return bp;
}

#endif

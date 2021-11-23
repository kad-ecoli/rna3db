/* calculate backbone torsion angles (omega, psi, phi, kappa, alpha) */
#ifndef BPtorsion2_HPP
#define BPtorsion2_HPP 1
#include <cstring>
#include "PDBParser.hpp"
#include "GeometryTools.hpp"
#include "StringTools.hpp"

using namespace std;

size_t dummy_base_pair(const ModelUnit &pdb_entry, 
    vector<char>&chain1_list, vector<char>&chain2_list,
    vector<string>&resn1_list, vector<string>&resn2_list,
    vector<int>&resi1_list,   vector<int>&resi2_list,
    vector<char>&icode1_list, vector<char>&icode2_list,
    vector<string>&name_list)
{
    size_t BPcount=0;
    size_t c1,c2,r1,r2;
    for (c1=0;c1<pdb_entry.chains.size();c1++)
    {
        for (r1=0;r1<pdb_entry.chains[c1].residues.size();r1++)
        {
            for (c2=c1;c2<pdb_entry.chains.size();c2++)
            {
                for (r2=(r1+1)*(c1==c2);r2<pdb_entry.chains[c2].residues.size();r2++)
                {
                    chain1_list.push_back(pdb_entry.chains[c1].chainID);
                    chain2_list.push_back(pdb_entry.chains[c2].chainID);
                    resn1_list.push_back(pdb_entry.chains[c1].residues[r1].resn);
                    resn2_list.push_back(pdb_entry.chains[c2].residues[r2].resn);
                    resi1_list.push_back(pdb_entry.chains[c1].residues[r1].resi);
                    resi2_list.push_back(pdb_entry.chains[c2].residues[r2].resi);
                    icode1_list.push_back(pdb_entry.chains[c1].residues[r1].icode);
                    icode2_list.push_back(pdb_entry.chains[c2].residues[r2].icode);
                    name_list.push_back("all         ");
                    BPcount++;
                }
            }
        }
    }
    return BPcount;
}

size_t read_dssr_base_pair(vector<char>&chain1_list, vector<char>&chain2_list,
    vector<string>&resn1_list, vector<string>&resn2_list,
    vector<int>&resi1_list,   vector<int>&resi2_list,
    vector<char>&icode1_list, vector<char>&icode2_list,
    vector<string>&name_list, const string dssrfile, const string bptypes)
{
    size_t BPcount=0;
    string line;
    vector<string> line_vec;
    string nt1,nt2,name;

    vector<string> bp_type_vec;
    split(bptypes, bp_type_vec, ',');
    size_t i;
    for (i=0;i<bp_type_vec.size();i++)
        bp_type_vec[i]+=((string)("           ")).substr(bp_type_vec[i].size());

    ifstream fp;
    bool use_stdin=(dssrfile=="-");
    if (!use_stdin) fp.open(dssrfile.c_str(),ios::in);
    while((use_stdin)?cin.good():fp.good())
    {
        if (use_stdin) getline(cin,line);
        else           getline(fp,line);
        
        if (line.size()==0) continue;
        else if (line[0]=='*')   continue;
        else if (line.substr(0,4)=="List") continue;
        else if (line.substr(0,8)=="     nt1") continue;

        name=line.substr(39,11);
        if (find(bp_type_vec.begin(),bp_type_vec.end(),name)==bp_type_vec.end()) continue;
        nt1=line.substr(5,14);
        nt2=line.substr(20,14);
        BPcount++;

        /* parse nt1 */ 
        split(nt1,line_vec,'.');
        if (line_vec.size()==2)
        {
            chain1_list.push_back(line_vec[0][0]);
            nt1=line_vec[1];
        }
        else chain1_list.push_back(' ');
        for (i=0;i<line_vec.size();i++) line_vec[i].clear();
        line_vec.clear();

        split(nt1,line_vec,'^');
        if (line_vec.size()==2)
        {
            icode1_list.push_back(line_vec[1][0]);
            nt1=line_vec[0];
        }
        else icode1_list.push_back(' ');
        for (i=0;i<line_vec.size();i++) line_vec[i].clear();
        line_vec.clear();

        if (nt1[0]=='D')
        {
            resn1_list.push_back(" "+nt1.substr(0,2));
            resi1_list.push_back(atoi(nt1.substr(2).c_str()));
        }
        else
        {
            resn1_list.push_back("  "+nt1.substr(0,1));
            resi1_list.push_back(atoi(nt1.substr(1).c_str()));
        }

        /* parse nt2 */ 
        split(nt2,line_vec,'.');
        if (line_vec.size()==2)
        {
            chain2_list.push_back(line_vec[0][0]);
            nt2=line_vec[1];
        }
        else chain2_list.push_back(' ');
        for (i=0;i<line_vec.size();i++) line_vec[i].clear();
        line_vec.clear();

        split(nt2,line_vec,'^');
        if (line_vec.size()==2)
        {
            icode2_list.push_back(line_vec[1][0]);
            nt2=line_vec[0];
        }
        else icode2_list.push_back(' ');
        for (i=0;i<line_vec.size();i++) line_vec[i].clear();
        line_vec.clear();

        if (nt2[0]=='D')
        {
            resn2_list.push_back(" "+nt2.substr(0,2));
            resi2_list.push_back(atoi(nt2.substr(2).c_str()));
        }
        else
        {
            resn2_list.push_back("  "+nt2.substr(0,1));
            resi2_list.push_back(atoi(nt2.substr(1).c_str()));
        }
        name_list.push_back(name);
    }
    if (!use_stdin) fp.close();

    /* clean up */
    line.clear();
    nt1.clear();
    nt2.clear();
    name.clear();
    vector<string>().swap(bp_type_vec);
    return BPcount;
}

void BPtorsion2(const vector<char>&chain1_list, const vector<char>&chain2_list,
    const vector<string>&resn1_list, const vector<string>&resn2_list,
    const vector<int>&resi1_list, const vector<int>&resi2_list,
    const vector<char>&icode1_list, const vector<char>&icode2_list,
    const ModelUnit &pdb_entry, const bool show_tor, const bool show_len,
    const bool show_ang, vector<vector<float> >&BPtorMat,
    vector<vector<float> >&BPlenMat, vector<vector<float> >&BPangMat)
{
    vector<float> tmp_tor(3,-360.);
    vector<float> tmp_len(3,-1.);
    vector<float> tmp_ang(3,-360.);

    size_t BPcount=chain1_list.size();

    if (show_tor) BPtorMat.assign(BPcount,tmp_tor);
    if (show_len) BPlenMat.assign(BPcount,tmp_len);
    if (show_ang) BPangMat.assign(BPcount,tmp_ang);

    /* coordinates of current residue */
    vector<float>Pi(3,0);  bool has_Pi;  vector<float>Pj(3,0);  bool has_Pj; 
    vector<float>C4i(3,0); bool has_C4i; vector<float>C4j(3,0); bool has_C4j;
    vector<float>Nxi(3,0); bool has_Nxi; vector<float>Nxj(3,0); bool has_Nxj;

    size_t bp, c1, c2, r1, r2, a1, a2;
    bool found_nt1, found_nt2;
    char base;
    for (bp=0;bp<BPcount;bp++)
    {
        found_nt1=found_nt2=false;
        has_Pi   =has_Pj   =false;
        has_C4i  =has_C4j  =false;
        has_Nxi  =has_Nxj  =false;

        for (c1=0; c1<pdb_entry.chains.size(); c1++)
        {
            if (pdb_entry.chains[c1].chainID!=chain1_list[bp]) continue;
            for (r1=0; r1<pdb_entry.chains[c1].residues.size(); r1++)
            {
                if (pdb_entry.chains[c1].residues[r1].resi !=resi1_list[bp]  ||
                    pdb_entry.chains[c1].residues[r1].resn !=resn1_list[bp]  ||
                    pdb_entry.chains[c1].residues[r1].icode!=icode1_list[bp]) continue;
                found_nt1=true;
                break;
            }
            if (found_nt1) break;
        }
        for (c2=0; c2<pdb_entry.chains.size(); c2++)
        {
            if (pdb_entry.chains[c2].chainID!=chain2_list[bp]) continue;
            for (r2=0; r2<pdb_entry.chains[c2].residues.size(); r2++)
            {
                if (pdb_entry.chains[c2].residues[r2].resi !=resi2_list[bp]  ||
                    pdb_entry.chains[c2].residues[r2].resn !=resn2_list[bp]  ||
                    pdb_entry.chains[c2].residues[r2].icode!=icode2_list[bp]) continue;
                found_nt2=true;
                break;
            }
            if (found_nt2) break;
        }

        if (!found_nt1 || !found_nt2)
        {
            cerr<<"ERROR! Cannot find base pair "
                <<chain1_list[bp]<<'.'<<resi1_list[bp]<<resn1_list[bp]<<icode1_list[bp]<<'\t'
                <<chain2_list[bp]<<'.'<<resi2_list[bp]<<resn2_list[bp]<<icode2_list[bp]<<endl;
            continue;
        }

        base=tolower(pdb_entry.chains[c1].residues[r1].resn[2]);
        for (a1=0;a1<pdb_entry.chains[c1].residues[r1].atoms.size();a1++)
        {
            if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" P  ")
            {
                has_Pi=true;
                Pi=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
            }
            else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" C4'")
            {
                has_C4i=true;
                C4i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
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
            else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" C4'")
            {
                has_C4j=true;
                C4j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
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

        if (show_tor)
        {
            if (has_Pi  && has_C4i && has_C4j && has_Pj ) BPtorMat[bp][0]=rad2deg(Points2Dihedral( Pi,C4i,C4j, Pj)); // PCCP: P[i]-C4'[i]-C4'[j]-P[j]
            if (has_C4i && has_Nxi && has_Nxj && has_C4j) BPtorMat[bp][1]=rad2deg(Points2Dihedral(C4i,Nxi,Nxj,C4j)); // CNNC: C4'[i]-N[i]-N[j]-C4'[j]
            if (has_Pi  && has_Nxi && has_Nxj && has_Pj ) BPtorMat[bp][2]=rad2deg(Points2Dihedral( Pi,Nxi,Nxj, Pj)); // PNNP: P[i]-N[i]-N[j]-P[j]
        }
        if (show_len)
        {
            if (has_Pi  && has_Pj)  BPlenMat[bp][0]=Points2Distance( Pi, Pj); //   PP: P[i]-P[j]
            if (has_C4i && has_C4j) BPlenMat[bp][1]=Points2Distance(C4i,C4j); // C4C4: C4'[i]-C4'[j]
            if (has_Nxi && has_Nxj) BPlenMat[bp][2]=Points2Distance(Nxi,Nxj); //   NN: N[i]-N[j]
        }
        if (show_ang)
        {
            if (has_Pi  && has_C4i && has_Pj  && has_C4j) BPangMat[bp][0]=rad2deg(Points4Angle( Pi,C4i, Pj,C4j)); // aPC: <P[i]C4'[i],P[j]C4'[j]>
            if (has_C4i && has_Nxi && has_C4j && has_Nxj) BPangMat[bp][1]=rad2deg(Points4Angle(C4i,Nxi,C4j,Nxj)); // aCN: <C4'[i]N[i],C4'[j]N[j]>
            if (has_Pi  && has_Nxi && has_Pj  && has_Nxj) BPangMat[bp][2]=rad2deg(Points4Angle( Pi,Nxi, Pj,Nxj)); // aPN: <P[i]N[i],P[j]N[j]>
        }
    }

    /* clean up */
    tmp_tor.clear();
    tmp_len.clear();
    tmp_ang.clear();
    Pi.clear();  Pj.clear();
    C4i.clear(); C4j.clear();
    Nxi.clear(); Nxj.clear();
    return;
}

#endif

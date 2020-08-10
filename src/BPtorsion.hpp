/* calculate backbone torsion angles (omega, psi, phi, kappa, alpha) */
#ifndef BPtorsion_HPP
#define BPtorsion_HPP 1
#include <cstring>
#include "PDBParser.hpp"
#include "GeometryTools.hpp"
#include "StringTools.hpp"

using namespace std;

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

        //cout<<"nt1="<<nt1<<";nt2="<<nt2<<";name="<<name<<'.'<<endl;
        //cout<<chain1_list.back()<<'.'<<resn1_list.back()<<resi1_list.back()<<icode1_list.back()<<';'
            //<<chain2_list.back()<<'.'<<resn2_list.back()<<resi2_list.back()<<icode2_list.back()<<';'<<name<<endl;
    }
    if (!use_stdin) fp.close();

    /* clean up */
    line.clear();
    nt1.clear();
    nt2.clear();
    name.clear();
    //resn1.clear();
    //resn2.clear();
    vector<string>().swap(bp_type_vec);
    return BPcount;
}

void BPtorsion(const vector<char>&chain1_list, const vector<char>&chain2_list,
    const vector<string>&resn1_list, const vector<string>&resn2_list,
    const vector<int>&resi1_list, const vector<int>&resi2_list,
    const vector<char>&icode1_list, const vector<char>&icode2_list,
    const ModelUnit &pdb_entry, const bool show_tor, const bool show_len,
    const bool show_ang, vector<vector<float> >&BPtorMat,
    vector<vector<float> >&BPlenMat, vector<vector<float> >&BPangMat)
{
    vector<float> tmp_tor(21,-360.);
    vector<float> tmp_len(10,-1.);
    vector<float> tmp_ang(21,-360.);

    size_t BPcount=chain1_list.size();

    if (show_tor) BPtorMat.assign(BPcount,tmp_tor);
    if (show_len) BPlenMat.assign(BPcount,tmp_len);
    if (show_ang) BPangMat.assign(BPcount,tmp_ang);

    /* coordinates of previous residue */
    vector<float>prev_Pi(3,0);  bool has_prev_Pi;  vector<float>prev_Pj(3,0);  bool has_prev_Pj; 
    vector<float>prev_O5i(3,0); bool has_prev_O5i; vector<float>prev_O5j(3,0); bool has_prev_O5j;
    vector<float>prev_C5i(3,0); bool has_prev_C5i; vector<float>prev_C5j(3,0); bool has_prev_C5j;
    vector<float>prev_C4i(3,0); bool has_prev_C4i; vector<float>prev_C4j(3,0); bool has_prev_C4j;
    vector<float>prev_C3i(3,0); bool has_prev_C3i; vector<float>prev_C3j(3,0); bool has_prev_C3j;
    vector<float>prev_C2i(3,0); bool has_prev_C2i; vector<float>prev_C2j(3,0); bool has_prev_C2j;
    vector<float>prev_C1i(3,0); bool has_prev_C1i; vector<float>prev_C1j(3,0); bool has_prev_C1j;
    vector<float>prev_O4i(3,0); bool has_prev_O4i; vector<float>prev_O4j(3,0); bool has_prev_O4j;
    vector<float>prev_O3i(3,0); bool has_prev_O3i; vector<float>prev_O3j(3,0); bool has_prev_O3j;
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
    /* coordinates of next residue */
    vector<float>next_Pi(3,0);  bool has_next_Pi;  vector<float>next_Pj(3,0);  bool has_next_Pj; 
    vector<float>next_O5i(3,0); bool has_next_O5i; vector<float>next_O5j(3,0); bool has_next_O5j;
    vector<float>next_C5i(3,0); bool has_next_C5i; vector<float>next_C5j(3,0); bool has_next_C5j;
    vector<float>next_C4i(3,0); bool has_next_C4i; vector<float>next_C4j(3,0); bool has_next_C4j;
    vector<float>next_C3i(3,0); bool has_next_C3i; vector<float>next_C3j(3,0); bool has_next_C3j;
    vector<float>next_C2i(3,0); bool has_next_C2i; vector<float>next_C2j(3,0); bool has_next_C2j;
    vector<float>next_C1i(3,0); bool has_next_C1i; vector<float>next_C1j(3,0); bool has_next_C1j;
    vector<float>next_O4i(3,0); bool has_next_O4i; vector<float>next_O4j(3,0); bool has_next_O4j;
    vector<float>next_O3i(3,0); bool has_next_O3i; vector<float>next_O3j(3,0); bool has_next_O3j;

    size_t bp, c1, c2, r1, r2, a1, a2;
    bool found_nt1, found_nt2;
    char base;
    for (bp=0;bp<BPcount;bp++)
    {
        found_nt1    = found_nt2    = false;
        has_prev_Pi  = has_prev_Pj  = false;
        has_prev_O5i = has_prev_O5j = false;
        has_prev_C5i = has_prev_C5j = false;
        has_prev_C4i = has_prev_C4j = false;
        has_prev_C3i = has_prev_C3j = false;
        has_prev_C2i = has_prev_C2j = false;
        has_prev_C1i = has_prev_C1j = false;
        has_prev_O4i = has_prev_O4j = false;
        has_prev_O3i = has_prev_O3j = false;
        has_Pi       = has_Pj       = false;
        has_O5i      = has_O5j      = false;
        has_C5i      = has_C5j      = false;
        has_C4i      = has_C4j      = false;
        has_C3i      = has_C3j      = false;
        has_C2i      = has_C2j      = false;
        has_C1i      = has_C1j      = false;
        has_O4i      = has_O4j      = false;
        has_O3i      = has_O3j      = false;
        has_Nxi      = has_Nxj      = false;
        has_next_Pi  = has_next_Pj  = false;
        has_next_O5i = has_next_O5j = false;
        has_next_C5i = has_next_C5j = false;
        has_next_C4i = has_next_C4j = false;
        has_next_C3i = has_next_C3j = false;
        has_next_C2i = has_next_C2j = false;
        has_next_C1i = has_next_C1j = false;
        has_next_O4i = has_next_O4j = false;
        has_next_O3i = has_next_O3j = false;

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

        if (r1>0 && (show_tor || show_ang))
        {
            for (a1=0;a1<pdb_entry.chains[c1].residues[r1-1].atoms.size();a1++)
            {
                if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" P  ")
                {
                    has_prev_Pi=true;
                    prev_Pi=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" O5'")
                {
                    has_prev_O5i=true;
                    prev_O5i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" C5'")
                {
                    has_prev_C5i=true;
                    prev_C5i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" C4'")
                {
                    has_prev_C4i=true;
                    prev_C4i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" C3'")
                {
                    has_prev_C3i=true;
                    prev_C3i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" C2'")
                {
                    has_prev_C2i=true;
                    prev_C2i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" C1'")
                {
                    has_prev_C1i=true;
                    prev_C1i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" O4'")
                {
                    has_prev_O4i=true;
                    prev_O4i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" O3'")
                {
                    has_prev_O3i=true;
                    prev_O3i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                }
            }
        }

        if (r2>0 && (show_tor || show_ang))
        {
            for (a2=0;a2<pdb_entry.chains[c2].residues[r2-1].atoms.size();a2++)
            {
                if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" P  ")
                {
                    has_prev_Pj=true;
                    prev_Pj=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                }
                else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" O5'")
                {
                    has_prev_O5j=true;
                    prev_O5j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                }
                else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" C5'")
                {
                    has_prev_C5j=true;
                    prev_C5j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                }
                else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" C4'")
                {
                    has_prev_C4j=true;
                    prev_C4j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                }
                else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" C3'")
                {
                    has_prev_C3j=true;
                    prev_C3j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                }
                else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" C2'")
                {
                    has_prev_C2j=true;
                    prev_C2j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                }
                else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" C1'")
                {
                    has_prev_C1j=true;
                    prev_C1j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                }
                else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" O4'")
                {
                    has_prev_O4j=true;
                    prev_O4j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                }
                else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" O3'")
                {
                    has_prev_O3j=true;
                    prev_O3j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                }
            }
        }

        if (r1<pdb_entry.chains[c1].residues.size()-1 && (show_tor || show_ang))
        {
            for (a1=0;a1<pdb_entry.chains[c1].residues[r1+1].atoms.size();a1++)
            {
                if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" P  ")
                {
                    has_next_Pi=true;
                    next_Pi=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" O5'")
                {
                    has_next_O5i=true;
                    next_O5i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" C5'")
                {
                    has_next_C5i=true;
                    next_C5i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" C4'")
                {
                    has_next_C4i=true;
                    next_C4i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" C3'")
                {
                    has_next_C3i=true;
                    next_C3i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" C2'")
                {
                    has_next_C2i=true;
                    next_C2i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" C1'")
                {
                    has_next_C1i=true;
                    next_C1i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" O4'")
                {
                    has_next_O4i=true;
                    next_O4i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" O3'")
                {
                    has_next_O3i=true;
                    next_O3i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                }
            }
        }

        if (r2<pdb_entry.chains[c2].residues.size()-1 && (show_tor || show_ang))
        {
            for (a2=0;a2<pdb_entry.chains[c2].residues[r2+1].atoms.size();a2++)
            {
                if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" P  ")
                {
                    has_next_Pj=true;
                    next_Pj=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                }
                else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" O5'")
                {
                    has_next_O5j=true;
                    next_O5j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                }
                else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" C5'")
                {
                    has_next_C5j=true;
                    next_C5j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                }
                else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" C4'")
                {
                    has_next_C4j=true;
                    next_C4j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                }
                else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" C3'")
                {
                    has_next_C3j=true;
                    next_C3j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                }
                else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" C2'")
                {
                    has_next_C2j=true;
                    next_C2j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                }
                else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" C1'")
                {
                    has_next_C1j=true;
                    next_C1j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                }
                else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" O4'")
                {
                    has_next_O4j=true;
                    next_O4j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                }
                else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" O3'")
                {
                    has_next_O3j=true;
                    next_O3j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                }
            }
        }

        if (show_tor)
        {
            if (has_prev_Pi  && has_Pi  && has_Pj  && has_next_Pj ) BPtorMat[bp][0] =rad2deg(Points2Dihedral( prev_Pi, Pi, Pj, next_Pj)); // Pm:  P[i-1]-P[i]-P[j]-P[j+1]
            if (has_next_Pi  && has_Pi  && has_Pj  && has_prev_Pj ) BPtorMat[bp][1] =rad2deg(Points2Dihedral( next_Pi, Pi, Pj, prev_Pj)); // Pp:  P[i+1]-P[i]-P[j]-P[j-1]
            if (has_prev_O5i && has_O5i && has_O5j && has_next_O5j) BPtorMat[bp][2] =rad2deg(Points2Dihedral(prev_O5i,O5i,O5j,next_O5j)); // O5m: O5'[i-1]-O5'[i]-O5'[j]-O5'[j+1]
            if (has_next_O5i && has_O5i && has_O5j && has_prev_O5j) BPtorMat[bp][3] =rad2deg(Points2Dihedral(next_O5i,O5i,O5j,prev_O5j)); // O5p: O5'[i+1]-O5'[i]-O5'[j]-O5'[j-1]
            if (has_prev_C5i && has_C5i && has_C5j && has_next_C5j) BPtorMat[bp][4] =rad2deg(Points2Dihedral(prev_C5i,C5i,C5j,next_C5j)); // C5m: C5'[i-1]-C5'[i]-C5'[j]-C5'[j+1]
            if (has_next_C5i && has_C5i && has_C5j && has_prev_C5j) BPtorMat[bp][5] =rad2deg(Points2Dihedral(next_C5i,C5i,C5j,prev_C5j)); // C5p: C5'[i+1]-C5'[i]-C5'[j]-C5'[j-1]
            if (has_prev_C4i && has_C4i && has_C4j && has_next_C4j) BPtorMat[bp][6] =rad2deg(Points2Dihedral(prev_C4i,C4i,C4j,next_C4j)); // C4m: C4'[i-1]-C4'[i]-C4'[j]-C4'[j+1]
            if (has_next_C4i && has_C4i && has_C4j && has_prev_C4j) BPtorMat[bp][7] =rad2deg(Points2Dihedral(next_C4i,C4i,C4j,prev_C4j)); // C4p: C4'[i+1]-C4'[i]-C4'[j]-C4'[j-1]
            if (has_prev_C3i && has_C3i && has_C3j && has_next_C3j) BPtorMat[bp][8] =rad2deg(Points2Dihedral(prev_C3i,C3i,C3j,next_C3j)); // C3m: C3'[i-1]-C3'[i]-C3'[j]-C3'[j+1]
            if (has_next_C3i && has_C3i && has_C3j && has_prev_C3j) BPtorMat[bp][9] =rad2deg(Points2Dihedral(next_C3i,C3i,C3j,prev_C3j)); // C3p: C3'[i+1]-C3'[i]-C3'[j]-C3'[j-1]
            if (has_prev_C2i && has_C2i && has_C2j && has_next_C2j) BPtorMat[bp][10]=rad2deg(Points2Dihedral(prev_C2i,C2i,C2j,next_C2j)); // C2m: C2'[i-1]-C2'[i]-C2'[j]-C2'[j+1]
            if (has_next_C2i && has_C2i && has_C2j && has_prev_C2j) BPtorMat[bp][11]=rad2deg(Points2Dihedral(next_C2i,C2i,C2j,prev_C2j)); // C2p: C2'[i+1]-C2'[i]-C2'[j]-C2'[j-1]
            if (has_prev_C1i && has_C1i && has_C1j && has_next_C1j) BPtorMat[bp][12]=rad2deg(Points2Dihedral(prev_C1i,C1i,C1j,next_C1j)); // C1m: C1'[i-1]-C1'[i]-C1'[j]-C1'[j+1]
            if (has_next_C1i && has_C1i && has_C1j && has_prev_C1j) BPtorMat[bp][13]=rad2deg(Points2Dihedral(next_C1i,C1i,C1j,prev_C1j)); // C1p: C1'[i+1]-C1'[i]-C1'[j]-C1'[j-1]
            if (has_prev_O4i && has_O4i && has_O4j && has_next_O4j) BPtorMat[bp][14]=rad2deg(Points2Dihedral(prev_O4i,O4i,O4j,next_O4j)); // O4m: O4'[i-1]-O4'[i]-O4'[j]-O4'[j+1]
            if (has_next_O4i && has_O4i && has_O4j && has_prev_O4j) BPtorMat[bp][15]=rad2deg(Points2Dihedral(next_O4i,O4i,O4j,prev_O4j)); // O4p: O4'[i+1]-O4'[i]-O4'[j]-O4'[j-1]
            if (has_prev_O3i && has_O3i && has_O3j && has_next_O3j) BPtorMat[bp][16]=rad2deg(Points2Dihedral(prev_O3i,O3i,O3j,next_O3j)); // O3m: O3'[i-1]-O3'[i]-O3'[j]-O3'[j+1]
            if (has_next_O3i && has_O3i && has_O3j && has_prev_O3j) BPtorMat[bp][17]=rad2deg(Points2Dihedral(next_O3i,O3i,O3j,prev_O3j)); // O3p: O3'[i+1]-O3'[i]-O3'[j]-O3'[j-1]
            if (has_Pi       && has_C4i && has_C4j && has_Pj      ) BPtorMat[bp][18]=rad2deg(Points2Dihedral(      Pi,C4i,C4j,      Pj)); // P44P: P[i]-C4'[i]-C4'[j]-P[j]
            if (has_C4i      && has_C1i && has_C1j && has_C4j     ) BPtorMat[bp][19]=rad2deg(Points2Dihedral(     C4i,C1i,C1j,     C4j)); // C4114C: C4'[i]-C1'[i]-C1'[j]-C4'[j]
            if (has_next_Pi  && has_Pi  && has_Pj  && has_next_Pj ) BPtorMat[bp][20]=rad2deg(Points2Dihedral( next_Pi, Pi, Pj, next_Pj)); // PppP: P[i+1]-P[i]-P[j]-P[j+1]
        }
        if (show_len)
        {
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
        }
        if (show_ang)
        {
            if (has_prev_Pi  && has_Pi  && has_next_Pj  && has_Pj ) BPangMat[bp][0] =rad2deg(Points4Angle( prev_Pi, Pi, next_Pj, Pj)); // aPm:      <P[i-1]P[i],P[j+1]P[j]>
            if (has_next_Pi  && has_Pi  && has_prev_Pj  && has_Pj ) BPangMat[bp][1] =rad2deg(Points4Angle( next_Pi, Pi, prev_Pj, Pj)); // aPp:      <P[i+1]P[i],P[j-1]P[j]>
            if (has_prev_O5i && has_O5i && has_next_O5j && has_O5j) BPangMat[bp][2] =rad2deg(Points4Angle(prev_O5i,O5i,next_O5j,O5j)); // aO5m: <O5'[i-1]O5'[i],O5'[j+1]O5'[j]> 
            if (has_next_O5i && has_O5i && has_prev_O5j && has_O5j) BPangMat[bp][3] =rad2deg(Points4Angle(next_O5i,O5i,prev_O5j,O5j)); // aO5p: <O5'[i+1]O5'[i],O5'[j-1]O5'[j]> 
            if (has_prev_C5i && has_C5i && has_next_C5j && has_C5j) BPangMat[bp][4] =rad2deg(Points4Angle(prev_C5i,C5i,next_C5j,C5j)); // aC5m: <C5'[i-1]C5'[i],C5'[j+1]C5'[j]> 
            if (has_next_C5i && has_C5i && has_prev_C5j && has_C5j) BPangMat[bp][5] =rad2deg(Points4Angle(next_C5i,C5i,prev_C5j,C5j)); // aC5p: <C5'[i+1]C5'[i],C5'[j-1]C5'[j]> 
            if (has_prev_C4i && has_C4i && has_next_C4j && has_C4j) BPangMat[bp][6] =rad2deg(Points4Angle(prev_C4i,C4i,next_C4j,C4j)); // aC4m: <C4'[i-1]C4'[i],C4'[j+1]C4'[j]> 
            if (has_next_C4i && has_C4i && has_prev_C4j && has_C4j) BPangMat[bp][7] =rad2deg(Points4Angle(next_C4i,C4i,prev_C4j,C4j)); // aC4p: <C4'[i+1]C4'[i],C4'[j-1]C4'[j]> 
            if (has_prev_C3i && has_C3i && has_next_C3j && has_C3j) BPangMat[bp][8] =rad2deg(Points4Angle(prev_C3i,C3i,next_C3j,C3j)); // aC3m: <C3'[i-1]C3'[i],C3'[j+1]C3'[j]> 
            if (has_next_C3i && has_C3i && has_prev_C3j && has_C3j) BPangMat[bp][9] =rad2deg(Points4Angle(next_C3i,C3i,prev_C3j,C3j)); // aC3p: <C3'[i+1]C3'[i],C3'[j-1]C3'[j]> 
            if (has_prev_C2i && has_C2i && has_next_C2j && has_C2j) BPangMat[bp][10]=rad2deg(Points4Angle(prev_C2i,C2i,next_C2j,C2j)); // aC2m: <C2'[i-1]C2'[i],C2'[j+1]C2'[j]> 
            if (has_next_C2i && has_C2i && has_prev_C2j && has_C2j) BPangMat[bp][11]=rad2deg(Points4Angle(next_C2i,C2i,prev_C2j,C2j)); // aC2p: <C2'[i+1]C2'[i],C2'[j-1]C2'[j]> 
            if (has_prev_C1i && has_C1i && has_next_C1j && has_C1j) BPangMat[bp][12]=rad2deg(Points4Angle(prev_C1i,C1i,next_C1j,C1j)); // aC1m: <C1'[i-1]C1'[i],C1'[j+1]C1'[j]> 
            if (has_next_C1i && has_C1i && has_prev_C1j && has_C1j) BPangMat[bp][13]=rad2deg(Points4Angle(next_C1i,C1i,prev_C1j,C1j)); // aC1p: <C1'[i+1]C1'[i],C1'[j-1]C1'[j]> 
            if (has_prev_O4i && has_O4i && has_next_O4j && has_O4j) BPangMat[bp][14]=rad2deg(Points4Angle(prev_O4i,O4i,next_O4j,O4j)); // aO4m: <O4'[i-1]O4'[i],O4'[j+1]O4'[j]> 
            if (has_next_O4i && has_O4i && has_prev_O4j && has_O4j) BPangMat[bp][15]=rad2deg(Points4Angle(next_O4i,O4i,prev_O4j,O4j)); // aO4p: <O4'[i+1]O4'[i],O4'[j-1]O4'[j]> 
            if (has_prev_O3i && has_O3i && has_next_O3j && has_O3j) BPangMat[bp][16]=rad2deg(Points4Angle(prev_O3i,O3i,next_O3j,O3j)); // aO3m: <O3'[i-1]O3'[i],O3'[j+1]O3'[j]> 
            if (has_next_O3i && has_O3i && has_prev_O3j && has_O3j) BPangMat[bp][17]=rad2deg(Points4Angle(next_O3i,O3i,prev_O3j,O3j)); // aO3p: <O3'[i+1]O3'[i],O3'[j-1]O3'[j]> 
            if (has_Pi       && has_C4i && has_Pj       && has_C4j) BPangMat[bp][18]=rad2deg(Points4Angle(      Pi,C4i,      Pj,C4j)); // aPC:      <P[i]C4'[i],P[j]C4'[j]>
            if (has_C4i      && has_C1i && has_C4j      && has_C1j) BPangMat[bp][19]=rad2deg(Points4Angle(     C4i,C1i,     C4j,C1j)); // aCC:    <C4'[i]C1'[i],C4'[j]C1'[j]>
            if (has_next_Pi  && has_Pi  && has_next_Pj  && has_Pj ) BPangMat[bp][20]=rad2deg(Points4Angle( next_Pi, Pi, next_Pj, Pj)); //aPppP:     <P[i+1]P[i],P[j+1]P[j]>
        }
    }

    /* clean up */
    tmp_tor.clear();
    tmp_len.clear();
    tmp_ang.clear();

    prev_Pi.clear();  prev_Pj.clear();
    prev_O5i.clear(); prev_O5j.clear();
    prev_C5i.clear(); prev_C5j.clear();
    prev_C4i.clear(); prev_C4j.clear();
    prev_C3i.clear(); prev_C3j.clear();
    prev_C2i.clear(); prev_C2j.clear();
    prev_C1i.clear(); prev_C1j.clear();
    prev_O4i.clear(); prev_O4j.clear();
    prev_O3i.clear(); prev_O3j.clear();

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

    next_Pi.clear();  next_Pj.clear();  
    next_O5i.clear(); next_O5j.clear();
    next_C5i.clear(); next_C5j.clear();
    next_C4i.clear(); next_C4j.clear();
    next_C3i.clear(); next_C3j.clear();
    next_C2i.clear(); next_C2j.clear();
    next_C1i.clear(); next_C1j.clear();
    next_O4i.clear(); next_O4j.clear();
    next_O3i.clear(); next_O3j.clear();
    return;
}

#endif

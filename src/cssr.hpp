/* calculate backbone torsion angles (omega, psi, phi, kappa, alpha) */
#ifndef cssr_HPP
#define cssr_HPP 1
#include <cstring>
#include <set>
#include <map>
#include "PDBParser.hpp"
#include "GeometryTools.hpp"
#include "BPstat.hpp"

using namespace std;

/* check if number of left and right bracket between residue r1 and r2 are the same */
bool check_balance_bracket(const vector<char>&dot_bracket, const size_t r1, const size_t r2)
{
    size_t left_count=0;
    size_t right_count=0;
    size_t r_left, r_right;
    if (r1<r2)
    {
        r_left=r1;
        r_right=r2;
    }
    else if (r2<r1)
    {
        r_left=r2;
        r_right=r1;
    }
    else return true; // r1==r2
    size_t r;
    for (r=r_left+1;r<r_right;r++)
    {
        left_count +=(dot_bracket[r]=='<' || dot_bracket[r]=='(' ||
                      dot_bracket[r]=='[' || dot_bracket[r]=='{');
        right_count+=(dot_bracket[r]=='>' || dot_bracket[r]==')' ||
                      dot_bracket[r]==']' || dot_bracket[r]=='}');
    }
    return left_count==right_count;
}

/* check type bracket
 * 0 - . ; 1 - <> ; 2 - () ; 3 - []; 4 - {}*/
size_t check_type_bracket(const vector<char>&dot_bracket, const size_t r1, const size_t r2)
{
    size_t max_type=0;
    size_t cur_type=0;
    size_t r_left, r_right;
    if (r1<r2)
    {
        r_left=r1;
        r_right=r2;
    }
    else if (r2<r1)
    {
        r_left=r2;
        r_right=r1;
    }
    else return 0; // r1==r2
    size_t r;
    for (r=r_left+1;r<r_right;r++)
    {
        switch (dot_bracket[r])
        {
            case '<': cur_type=1; break;
            case '>': cur_type=1; break;
            case '(': cur_type=2; break;
            case ')': cur_type=2; break;
            case '[': cur_type=3; break;
            case ']': cur_type=3; break;
            case '{': cur_type=4; break;
            case '}': cur_type=4; break;
            default:  cur_type=0; break;
        }
        if (cur_type>max_type) max_type=cur_type;
    }
    return max_type;
}

void cssr_bp2dot(const vector<string>&res_str_vec, 
    const vector<size_t> filtered_bp_vec,
    const vector<pair<float,vector<string> > >&bp_vec, vector<char>&dot_bracket)
{
    /* initialize dot_bracket */
    size_t Lch=res_str_vec.size();
    dot_bracket.assign(Lch,'.');

    /* map residue to string */
    vector<pair<string,size_t> > res2int_vec;
    size_t r;
    for (r=0;r<Lch;r++)
        res2int_vec.push_back(pair<string,size_t>(res_str_vec[r],r));
    map<string, size_t> res2int_map(res2int_vec.begin(),res2int_vec.end());
    vector<pair<string, size_t> >().swap(res2int_vec);

    /* convert bp_vec to int */
    vector<size_t> tmp_size_t(4,0); // r2 pairs_in_helix upstream_pairs downstream_pairs
    vector<vector<size_t> >bpint_vec(Lch,tmp_size_t);
    for (r=0;r<Lch;r++) bpint_vec[r][0]=r;
    size_t r1,r2,bp;
    for (bp=0;bp<bp_vec.size();bp++)
    {
        if (find(filtered_bp_vec.begin(), filtered_bp_vec.end(),bp
            )==filtered_bp_vec.end()) continue;
        r1=res2int_map[bp_vec[bp].second[0]];
        r2=res2int_map[bp_vec[bp].second[1]];
        bpint_vec[r1][0]=r2;
        bpint_vec[r2][0]=r1;
    }
    map<string, size_t>().swap(res2int_map);

    /* count the number of neighboring base pairs */
    size_t max_bp=0;
    for (r=0;r<Lch;r++)
    {
        r2=bpint_vec[r][0];
        if (r2==r) continue;
        if (r>0)
        {
            for (r1=r-1; r1>=0; r1--) // find upstream pairs
            {
                if (bpint_vec[r1][0]!=r2+1) break;
                if (r>bpint_vec[r][0] && r1<bpint_vec[r1][0]) break;
                r2++;
                bpint_vec[r][2]++;
                if (r1==0) break;
            }
        }
        r2=bpint_vec[r][0];
        if (r<Lch-1)
        {
            for (r1=r+1; r1<Lch; r1++) // find downstream pairs
            {
                if (bpint_vec[r1][0]!=r2-1) break;
                if (r<bpint_vec[r][0] && r1>bpint_vec[r1][0]) break;
                r2--;
                bpint_vec[r][3]++;
            }
        }
        bpint_vec[r][1]=bpint_vec[r][2]+1+bpint_vec[r][3];
        if (bpint_vec[r][1]>max_bp) max_bp=bpint_vec[r][1];
        //cout<<r<<'\t'<<bpint_vec[r][0]<<'\t'<<bpint_vec[r][1]<<'\t'<<bpint_vec[r][2]<<'\t'<<bpint_vec[r][3]<<endl;
    }

    /* assign nested ss */
    cerr<<"assign nested ss with <="<<max_bp<<" base pairs"<<endl;
    for (;max_bp>0; max_bp--)
    {
        //cout<<"max_bp="<<max_bp<<endl;
        for (r=0;r<Lch;r++)
        {
            if (bpint_vec[r][1]<max_bp) continue;
            r1=r;
            r2=bpint_vec[r][0];
            //cout<<r1<<'\t'<<r2<<'\t'<<bpint_vec[r][1]<<'\t'<<bpint_vec[r][2]<<'\t'<<bpint_vec[r][3]<<endl;
            if (!check_balance_bracket(dot_bracket, r1, r2)) continue;
            for (r1=r-bpint_vec[r][2]; r1<=r+bpint_vec[r][3]; r1++)
            {
                r2=bpint_vec[r1][0];
                if (r1<r2)
                {
                    dot_bracket[r1]='<';
                    dot_bracket[r2]='>';
                }
                else
                {
                    dot_bracket[r1]='>';
                    dot_bracket[r2]='<';
                }
                bpint_vec[r1][0]=r1;
                bpint_vec[r1][1]=0;
                bpint_vec[r2][0]=r2;
                bpint_vec[r2][1]=0;
            }
        }
    }
    
    /* assign knotted ss */
    max_bp=0;
    for (r=0;r<Lch;r++) if (bpint_vec[r][1]>max_bp) max_bp=bpint_vec[r][1];
    if (max_bp)
    {
        cerr<<"assign knotted ss with <="<<max_bp<<" base pairs"<<endl;
        char left_bracket_vec [4]={'<', '(', '[', '{'};
        char right_bracket_vec[4]={'>', ')', ']', '}'};
        size_t bracket_type=0;
        bool too_many_knot=false;
        for (;max_bp>0; max_bp--)
        {
            for (r=0;r<Lch;r++)
            {
                if (bpint_vec[r][1]<max_bp) continue;
                r1=r;
                r2=bpint_vec[r][0];
                bracket_type=check_type_bracket(dot_bracket, r1, r2);
                if (bracket_type>=4)
                {
                    too_many_knot=true;
                    continue;
                }
                for (r1=r-bpint_vec[r][2]; r1<=r+bpint_vec[r][3]; r1++)
                {
                    r2=bpint_vec[r1][0];
                    if (r1<r2)
                    {
                        dot_bracket[r1]=left_bracket_vec[bracket_type];
                        dot_bracket[r2]=right_bracket_vec[bracket_type];
                    }
                    else
                    {
                        dot_bracket[r1]=right_bracket_vec[bracket_type];
                        dot_bracket[r2]=left_bracket_vec[bracket_type];
                    }
                    bpint_vec[r1][0]=r1;
                    bpint_vec[r1][1]=0;
                    bpint_vec[r2][0]=r2;
                    bpint_vec[r2][1]=0;
                }
            }
        }
        if (too_many_knot) cerr<<"warning! some pseudoknots not shown"<<endl;
    }


    /* clean up */
    vector<vector<size_t> >().swap(bpint_vec);
    return;
}

/* if a residue is paired with multiple other residues, just keep the
 * highest scoring residue pair */
void filter_bp(const vector<pair<float,vector<string> > >&bp_vec,
    vector<size_t> &filtered_bp_vec)
{
    vector<string> paired_res_vec;
    vector<pair<float,size_t> > bp_idx_vec;
    size_t bp;
    for (bp=0;bp<bp_vec.size();bp++)
        bp_idx_vec.push_back(pair<float,size_t>(-bp_vec[bp].first,bp));
    sort(bp_idx_vec.begin(), bp_idx_vec.end());
    
    for (bp=0;bp<bp_idx_vec.size();bp++)
    {
        if (-bp_idx_vec[bp].first<0.5) break;
        if (find(paired_res_vec.begin(), paired_res_vec.end(), 
            bp_vec[bp_idx_vec[bp].second].second[0])!=paired_res_vec.end())
            continue;
        if (find(paired_res_vec.begin(), paired_res_vec.end(),
            bp_vec[bp_idx_vec[bp].second].second[1])!=paired_res_vec.end())
            continue;
        filtered_bp_vec.push_back(bp_idx_vec[bp].second);
        paired_res_vec.push_back(bp_vec[bp_idx_vec[bp].second].second[0]);
        paired_res_vec.push_back(bp_vec[bp_idx_vec[bp].second].second[1]);
    }
    vector<string>().swap(paired_res_vec);
    vector<pair<float,size_t> >().swap(bp_idx_vec);
    return;
}

void removeUnspecifiedAtom(ModelUnit &pdb_entry, const string atom)
{
    size_t c,r,a;
    for (c=0;c<pdb_entry.chains.size();c++)
    {
        for (r=0;r<pdb_entry.chains[c].residues.size();r++)
        {
            for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++)
            {
                if (pdb_entry.chains[c].residues[r].atoms[a].name!=atom)
                    continue;
                if (a!=0)
                {
                    pdb_entry.chains[c].residues[r].atoms[0].name=
                    pdb_entry.chains[c].residues[r].atoms[a].name;
                    pdb_entry.chains[c].residues[r].atoms[0].xyz[0]=
                    pdb_entry.chains[c].residues[r].atoms[a].xyz[0];
                    pdb_entry.chains[c].residues[r].atoms[0].xyz[1]=
                    pdb_entry.chains[c].residues[r].atoms[a].xyz[1];
                    pdb_entry.chains[c].residues[r].atoms[0].xyz[2]=
                    pdb_entry.chains[c].residues[r].atoms[a].xyz[2];
                }
                break;
            }
            if (pdb_entry.chains[c].residues[r].atoms.size()<=1) continue;
            for (a=pdb_entry.chains[c].residues[r].atoms.size()-1;a>=1;a--)
            {
                pdb_entry.chains[c].residues[r].atoms[a].name.clear();
                pdb_entry.chains[c].residues[r].atoms[a].xyz.clear();
                pdb_entry.chains[c].residues[r].atoms.pop_back();
            }
        }
    }
    return;
}

inline bool bp_tor_score(
    const vector<float> &c1, const vector<float> &c2,
    const vector<float> &c3, const vector<float> &c4,
    const float mu, const float sd, const float tol, const float weight,
    float &nominator, float &denominator, vector<float>&sd_tor)
{
    float tor=rad2deg(Points2Dihedral(c1,c2,c3,c4));
    if(tor<-180) return false;
    denominator+=weight;
    float diff=fabs(tor-mu); 
    if (diff>180) diff-=180;
    nominator+=weight*(1-diff/(tol* sd));
    sd_tor.push_back(sd);
    return true;
}

inline bool bp_ang_score(
    const vector<float> &c1, const vector<float> &c2,
    const vector<float> &c3, const vector<float> &c4,
    const float mu, const float sd, const float tol, const float weight,
    float &nominator, float &denominator, vector<float>&sd_ang)
{
    float ang=rad2deg(Points4Angle(c1,c2,c3,c4));
    if(ang<-180) return false;
    denominator+=weight;
    nominator+=weight*(1-fabs(ang-mu)/(tol* sd));
    sd_ang.push_back(sd);
    return true;
}

inline bool bp_len_score(const vector<float>&c1, const vector<float>&c2,
    const float mu, const float sd, const float tol, const float weight,
    float &nominator, float &denominator, vector<float>&sd_len)
{
    denominator+=weight;
    nominator+=weight*(1-fabs(Points2Distance(c1,c2)-mu)/(tol*sd));
    sd_len.push_back(sd);
    return true;
}

inline bool bp_nn_score(const bool previnextj, const bool nextiprevj,
    const bool has_prev_ci,       const bool has_next_ci,
    const bool has_prev_cj,       const bool has_next_cj,
    const vector<float> &prev_ci, const vector<float> &next_ci,
    const vector<float> &prev_cj, const vector<float> &next_cj,
    const float mu, const float sd, const float tol, const float weight,
    float &nominator, float &denominator, vector<float>&sd_nn)
{
    if (!(has_prev_ci && has_next_cj) && !(has_next_ci && has_prev_cj))
        return false;
    denominator+=weight;
    float nominator_previnextj=0;
    float nominator_nextiprevj=0;
    if (has_prev_ci  && has_next_cj )
    {
        nominator_previnextj=1-fabs(Points2Distance(prev_ci,next_cj)-mu)/(tol*sd);
        if (!previnextj && nominator_previnextj>0) nominator_previnextj=0;
    }
    if (has_next_ci  && has_prev_cj )
    {
        nominator_nextiprevj=1-fabs(Points2Distance(next_ci,prev_cj)-mu)/(tol*sd);
        if (!nextiprevj && nominator_nextiprevj>0) nominator_nextiprevj=0;
    }

    if (nominator_previnextj>nominator_nextiprevj)
         nominator+=weight*nominator_previnextj;
    else nominator+=weight*nominator_nextiprevj;
    sd_nn.push_back(sd);
    return true;
}

inline float mean_vec(const vector<float>&vec)
{
    float sum=0;
    int i;
    for (i=0;i<vec.size();i++) sum+=vec[i];
    return sum/vec.size();
}

float adjuststd(const vector<float>&sd_len,const vector<float>&sd_nn,
    const vector<float>&sd_tor,const vector<float>&sd_ang,
    const float sd_len_mean, const float sd_nn_mean,
    const float sd_tor_mean, const float sd_ang_mean,
    const float weight_len, const float weight_nn ,
    const float weight_tor, const float weight_ang)
{
    return (weight_len*((sd_len.size()>0)?(sd_len_mean-mean_vec(sd_len)):0)+
            weight_nn *((sd_nn.size() >0)?(sd_nn_mean -mean_vec(sd_nn )):0)+
            weight_tor*((sd_tor.size()>0)?(sd_tor_mean-mean_vec(sd_tor)):0)+
            weight_ang*((sd_ang.size()>0)?(sd_ang_mean-mean_vec(sd_ang)):0))
          /(weight_len*(sd_len.size()>0)+
            weight_nn *(sd_nn.size() >0)+
            weight_tor*(sd_tor.size()>0)+
            weight_ang*(sd_ang.size()>0));
}

void cssr(const ModelUnit &pdb_entry, vector<string>&res_str_vec,
    vector<pair<float,vector<string> > >&bp_vec, const bool interchain)
{
    /* pre-trained parameters */
    float weight_len=1;
    float weight_nn =1;
    float weight_tor=1;
    float weight_ang=1;
    float tol=2; // tolerance, in the unit of standard deviation
    float adjust1=0.2;  // adjust for varying number of tests
    float adjust2=0.05; // adjust for varying std
    float adjust3=0.1;  // baseline value
    float totaltest=(10-1)*(weight_len>0)+
                    (10-1)*(weight_nn >0)+
                    ( 9-1)*(weight_tor>0)+
                    (10-1)*(weight_ang>0);
    
    /* other variables */
    vector<string> tmp_bp(3,"");
    float ang;
    size_t c1, c2, r1, r2, a1, a2;
    char base1,base2;
    char base1prev,base1next;
    char base2prev,base2next;
    char icode,chainID;
    stringstream ss;
    float nominator, denominator;
    float nominator_previnextj;
    float nominator_nextiprevj;
    int   successtest;
    bool previnextj,nextiprevj;
    vector<float> sd_len;
    vector<float> sd_nn;
    vector<float> sd_tor;
    vector<float> sd_ang;
    float sd_len_mean=(  PP_sd+ O5O5_sd+ C5C5_sd+ C4C4_sd+ C3C3_sd+
                       C2C2_sd+ C1C1_sd+ O4O4_sd+ O3O3_sd+ NN_sd)/10;
    float sd_nn_mean =(  PP_sd+ O5O5_sd+ C5C5_sd+ C4C4_sd+ C3C3_sd+
                       C2C2_sd+ C1C1_sd+ O4O4_sd+ O3O3_sd+ NN_sd)/10;
    float sd_tor_mean=(  Pp_sd+  O5p_sd+  C5p_sd+  C4p_sd+ C3p_sd+
                        C2p_sd+  C1p_sd+  O4p_sd+  O3m_sd)/9;
    float sd_ang_mean=( aPp_sd+ aO5p_sd+ aC5p_sd+ aC4p_sd+ aC3p_sd+
                       aC2p_sd+ aC1p_sd+ aO4p_sd+ aO3p_sd+ aCC_sd)/10;

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
    vector<float>prev_Nxi(3,0); bool has_prev_Nxi; vector<float>prev_Nxj(3,0); bool has_prev_Nxj;
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
    vector<float>next_Nxi(3,0); bool has_next_Nxi; vector<float>next_Nxj(3,0); bool has_next_Nxj;

    /* loop over base pairs */
    for (c1=0; c1<pdb_entry.chains.size(); c1++)
    {
        //cerr<<">chain_"<<c1<<':'<<pdb_entry.chains[c1].chainID<<'\t'
            //<<pdb_entry.chains[c1].residues.size()<<endl;
        for (r1=0; r1<pdb_entry.chains[c1].residues.size(); r1++)
        {
            base1=pdb_entry.chains[c1].residues[r1].resn[2];
            base1prev=base1next=0;
            //cerr<<base1<<endl;
            if (!pdb_entry.chains[c1].residues[r1].atoms.size()) continue;
            if ((chainID=pdb_entry.chains[c1].chainID)!=' ') 
                ss<<chainID<<".";
            ss<<base1<<pdb_entry.chains[c1].residues[r1].resi;
            if ((icode=pdb_entry.chains[c1].residues[r1].icode)!=' ')
                ss<<'^'<<icode;
            tmp_bp[0]=ss.str();
            res_str_vec.push_back(tmp_bp[0]);
            ss.str(string());

            has_prev_Pi  = false;
            has_prev_O5i = false;
            has_prev_C5i = false;
            has_prev_C4i = false;
            has_prev_C3i = false;
            has_prev_C2i = false;
            has_prev_C1i = false;
            has_prev_O4i = false;
            has_prev_O3i = false;
            has_prev_Nxi = false;
            has_Pi       = false;
            has_O5i      = false;
            has_C5i      = false;
            has_C4i      = false;
            has_C3i      = false;
            has_C2i      = false;
            has_C1i      = false;
            has_O4i      = false;
            has_O3i      = false;
            has_Nxi      = false;
            has_next_Pi  = false;
            has_next_O5i = false;
            has_next_C5i = false;
            has_next_C4i = false;
            has_next_C3i = false;
            has_next_C2i = false;
            has_next_C1i = false;
            has_next_O4i = false;
            has_next_O3i = false;
            has_next_Nxi = false;

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
                        && (base1=='A' || base1=='G')) ||
                        (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" N1 "
                        && (base1=='C' || base1=='T' || base1=='U')))
                {
                    has_Nxi=true;
                    Nxi=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
            }

            if (r1>0 && (weight_tor>0 || weight_ang>0 || weight_nn>0))
            {
                base1prev=pdb_entry.chains[c1].residues[r1-1].resn[2];
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
                    else if ((pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" N9 "
                            && (base1prev=='A' || base1prev=='G')) ||
                             (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" N1 "
                            && (base1prev=='C' || base1prev=='T' || base1prev=='U')))
                    {
                        has_prev_Nxi=true;
                        prev_Nxi=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                }
            }

            if (r1<pdb_entry.chains[c1].residues.size()-1 && (weight_tor>0 || weight_ang>0 || weight_nn>0))
            {
                base1next=pdb_entry.chains[c1].residues[r1+1].resn[2];
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
                    else if ((pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" N9 "
                            && (base1next=='A' || base1next=='G')) ||
                             (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" N1 "
                            && (base1next=='C' || base1next=='T' || base1next=='U')))
                    {
                        has_next_Nxi=true;
                        next_Nxi=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                }
            }


            for (c2=c1+interchain; c2<pdb_entry.chains.size(); c2++)
            {
                for (r2=(c1==c2)*(r1+1); r2<pdb_entry.chains[c2].residues.size(); r2++)
                {
                    if (!pdb_entry.chains[c2].residues[r2].atoms.size() ||
                        Points2Distance2(
                        pdb_entry.chains[c1].residues[r1].atoms[0].xyz,
                        pdb_entry.chains[c2].residues[r2].atoms[0].xyz)>530)
                        continue;
                    base2=pdb_entry.chains[c2].residues[r2].resn[2];
                    base2prev=base2next=0;
                    if ((base1=='A' &&(base2=='U' || base2=='T')) ||
                        (base1=='C' && base2=='G')||
                       ((base1=='U' || base1=='T')&& base2=='A')  ||
                        (base1=='G' && base2=='C'))
                        tmp_bp[2]=base1+string("-")+base2+" WC";
                    else if ((base1=='G' && base2=='U') ||
                             (base1=='U' && base2=='G')) 
                        tmp_bp[2]=base1+string("-")+base2+" Wobble";
                    else continue;

                    if ((chainID=pdb_entry.chains[c2].chainID)!=' ') 
                        ss<<chainID<<".";
                    ss<<base2<<pdb_entry.chains[c2].residues[r2].resi;
                    if ((icode=pdb_entry.chains[c2].residues[r2].icode)!=' ')
                        ss<<'^'<<icode;
                    tmp_bp[1]=ss.str();
                    ss.str(string());

                    has_prev_Pj  = false;
                    has_prev_O5j = false;
                    has_prev_C5j = false;
                    has_prev_C4j = false;
                    has_prev_C3j = false;
                    has_prev_C2j = false;
                    has_prev_C1j = false;
                    has_prev_O4j = false;
                    has_prev_O3j = false;
                    has_Pj       = false;
                    has_O5j      = false;
                    has_C5j      = false;
                    has_C4j      = false;
                    has_C3j      = false;
                    has_C2j      = false;
                    has_C1j      = false;
                    has_O4j      = false;
                    has_O3j      = false;
                    has_Nxj      = false;
                    has_next_Pj  = false;
                    has_next_O5j = false;
                    has_next_C5j = false;
                    has_next_C4j = false;
                    has_next_C3j = false;
                    has_next_C2j = false;
                    has_next_C1j = false;
                    has_next_O4j = false;
                    has_next_O3j = false;

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
                                && (base2=='A' || base2=='G')) ||
                                (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" N1 "
                                && (base2=='C' || base2=='T' || base2=='U')))
                        {
                            has_Nxj=true;
                            Nxj=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                    }

                    if (r2>0 && (weight_tor>0 || weight_ang>0 || weight_nn>0) && c2-c1+r2-r1>1)
                    {
                        base2prev=pdb_entry.chains[c2].residues[r2-1].resn[2];
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
                            else if ((pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" N9 "
                                  && (base2prev=='A' || base2prev=='G')) ||
                                     (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" N1 "
                                  && (base2prev=='C' || base2prev=='T' || base2prev=='U')))
                            {
                                has_prev_Nxj=true;
                                prev_Nxj=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                        }
                    }

                    if (r2<pdb_entry.chains[c2].residues.size()-1 && 
                       (weight_tor>0 || weight_ang>0 || weight_nn>0) && c2-c1+r2-r1>1)
                    {
                        base2next=pdb_entry.chains[c2].residues[r2+1].resn[2];
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
                            else if ((pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" N9 "
                                  && (base2next=='A' || base2next=='G')) ||
                                     (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" N1 "
                                  && (base2next=='C' || base2next=='T' || base2next=='U')))
                            {
                                has_next_Nxj=true;
                                next_Nxj=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                        }
                    }

                    nominator=denominator=successtest=0;
                    if (weight_len>0)
                    {
                        successtest--;
                        if (has_Pi  && has_Pj ) successtest+=bp_len_score( Pi, Pj,  PP_mu,  PP_sd,tol,weight_len,nominator,denominator,sd_len); //   PP: P[i]-P[j]
                        if (has_O5i && has_O5j) successtest+=bp_len_score(O5i,O5j,O5O5_mu,O5O5_sd,tol,weight_len,nominator,denominator,sd_len); // O5O5: O5'[i]-O5'[j]
                        if (has_C5i && has_C5j) successtest+=bp_len_score(C5i,C5j,C5C5_mu,C5C5_sd,tol,weight_len,nominator,denominator,sd_len); // C5C5: C5'[i]-C5'[j]
                        if (has_C4i && has_C4j) successtest+=bp_len_score(C4i,C4j,C4C4_mu,C4C4_sd,tol,weight_len,nominator,denominator,sd_len); // C4C4: C4'[i]-C4'[j]
                        if (has_C3i && has_C3j) successtest+=bp_len_score(C3i,C3j,C3C3_mu,C3C3_sd,tol,weight_len,nominator,denominator,sd_len); // C3C3: C3'[i]-C3'[j]
                        if (has_C2i && has_C2j) successtest+=bp_len_score(C2i,C2j,C2C2_mu,C2C2_sd,tol,weight_len,nominator,denominator,sd_len); // C2C2: C2'[i]-C2'[j]
                        if (has_C1i && has_C1j) successtest+=bp_len_score(C1i,C1j,C1C1_mu,C1C1_sd,tol,weight_len,nominator,denominator,sd_len); // C1C1: C1'[i]-C1'[j]
                        if (has_O4i && has_O4j) successtest+=bp_len_score(O4i,O4j,O4O4_mu,O4O4_sd,tol,weight_len,nominator,denominator,sd_len); // O4O4: O4'[i]-O4'[j]
                        if (has_O3i && has_O3j) successtest+=bp_len_score(O3i,O3j,O3O3_mu,O3O3_sd,tol,weight_len,nominator,denominator,sd_len); // O3O3: O3'[i]-O3'[j]
                        if (has_Nxi && has_Nxj) successtest+=bp_len_score(Nxi,Nxj,  NN_mu,  NN_sd,tol,weight_len,nominator,denominator,sd_len); //   NN: N[i]-N[j]
                    }
                    if (weight_nn>0)
                    {
                        successtest--;
                        previnextj=((base1prev=='A' &&(base2next=='U' || base2next=='T')) ||
                                    (base1prev=='C' && base2next=='G')||
                                   ((base1prev=='U' || base1prev=='T')&& base2next=='A')  ||
                                    (base1prev=='G' && base2next=='C'));
                        nextiprevj=((base1next=='A' &&(base2prev=='U' || base2prev=='T')) ||
                                    (base1next=='C' && base2prev=='G')||
                                   ((base1next=='U' || base1next=='T')&& base2prev=='A')  ||
                                    (base1next=='G' && base2prev=='C'));
                        if (has_Pi  && has_Pj) //   PP: P[i]-P[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_Pi, has_next_Pi, has_prev_Pj, has_next_Pj,
                                prev_Pi, next_Pi, prev_Pj, next_Pj,   PP_mu,  PP_sd,tol,weight_nn,nominator,denominator,sd_len);
                        if (has_O5i && has_O5j) // O5O5: O5'[i]-O5'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_O5i,has_next_O5i,has_prev_O5j,has_next_O5j,
                                prev_O5i,next_O5i,prev_O5j,next_O5j,O5O5_mu,O5O5_sd,tol,weight_nn,nominator,denominator,sd_len);
                        if (has_C5i && has_C5j) // C5C5: C5'[i]-C5'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_C5i,has_next_C5i,has_prev_C5j,has_next_C5j,
                                prev_C5i,next_C5i,prev_C5j,next_C5j,C5C5_mu,C5C5_sd,tol,weight_nn,nominator,denominator,sd_len);
                        if (has_C4i && has_C4j) // C4C4: C4'[i]-C4'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_C4i,has_next_C4i,has_prev_C4j,has_next_C4j,
                                prev_C4i,next_C4i,prev_C4j,next_C4j,C4C4_mu,C4C4_sd,tol,weight_nn,nominator,denominator,sd_len);
                        if (has_C3i && has_C3j) // C3C3: C3'[i]-C3'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_C3i,has_next_C3i,has_prev_C3j,has_next_C3j,
                                prev_C3i,next_C3i,prev_C3j,next_C3j,C3C3_mu,C3C3_sd,tol,weight_nn,nominator,denominator,sd_len);
                        if (has_C2i && has_C2j) // C2C2: C2'[i]-C2'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_C2i,has_next_C2i,has_prev_C2j,has_next_C2j,
                                prev_C2i,next_C2i,prev_C2j,next_C2j,C2C2_mu,C2C2_sd,tol,weight_nn,nominator,denominator,sd_len);
                        if (has_C1i && has_C1j) // C1C1: C1'[i]-C1'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_C1i,has_next_C1i,has_prev_C1j,has_next_C1j,
                                prev_C1i,next_C1i,prev_C1j,next_C1j,C1C1_mu,C1C1_sd,tol,weight_nn,nominator,denominator,sd_len);
                        if (has_O4i && has_O4j) // O4O4: O4'[i]-O4'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_O4i,has_next_O4i,has_prev_O4j,has_next_O4j,
                                prev_O4i,next_O4i,prev_O4j,next_O4j,O4O4_mu,O4O4_sd,tol,weight_nn,nominator,denominator,sd_len);
                        if (has_O3i && has_O3j) // O3O3: O3'[i]-O3'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_O3i,has_next_O3i,has_prev_O3j,has_next_O3j,
                                prev_O3i,next_O3i,prev_O3j,next_O3j,O3O3_mu,O3O3_sd,tol,weight_nn,nominator,denominator,sd_len);
                        if (has_Nxi && has_Nxj) //   NN: N[i]-N[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_Nxi,has_next_Nxi,has_prev_Nxj,has_next_Nxj,
                                prev_Nxi,next_Nxi,prev_Nxj,next_Nxj,  NN_mu,  NN_sd,tol,weight_nn,nominator,denominator,sd_len);
                    }
                    if (weight_tor>0)
                    {
                        successtest--;
                        if (has_next_Pi  && has_Pi  && has_Pj  && has_next_Pj ) successtest+=bp_tor_score( next_Pi, Pi, Pj, next_Pj, Pp_mu, Pp_sd, tol,weight_tor,nominator,denominator,sd_tor); //  Pp:  P[i+1]-P[i]-P[j]-P[j-1]
                        if (has_next_O5i && has_O5i && has_O5j && has_next_O5j) successtest+=bp_tor_score(next_O5i,O5i,O5j,next_O5j,O5p_mu,O5p_sd, tol,weight_tor,nominator,denominator,sd_tor); // O5p: O5'[i+1]-O5'[i]-O5'[j]-O5'[j-1]
                        if (has_next_C5i && has_C5i && has_C5j && has_next_C5j) successtest+=bp_tor_score(next_C5i,C5i,C5j,next_C5j,C5p_mu,C5p_sd, tol,weight_tor,nominator,denominator,sd_tor); // C5p: C5'[i+1]-C5'[i]-C5'[j]-C5'[j-1]
                        if (has_next_C4i && has_C4i && has_C4j && has_next_C4j) successtest+=bp_tor_score(next_C4i,C4i,C4j,next_C4j,C4p_mu,C4p_sd, tol,weight_tor,nominator,denominator,sd_tor); // C4p: C4'[i+1]-C4'[i]-C4'[j]-C4'[j-1]
                        if (has_next_C3i && has_C3i && has_C3j && has_next_C3j) successtest+=bp_tor_score(next_C3i,C3i,C3j,next_C3j,C3p_mu,C3p_sd, tol,weight_tor,nominator,denominator,sd_tor); // C3p: C3'[i+1]-C3'[i]-C3'[j]-C3'[j-1]
                        if (has_next_C2i && has_C2i && has_C2j && has_next_C2j) successtest+=bp_tor_score(next_C2i,C2i,C2j,next_C2j,C2p_mu,C2p_sd, tol,weight_tor,nominator,denominator,sd_tor); // C2p: C2'[i+1]-C2'[i]-C2'[j]-C2'[j-1]
                        if (has_next_C1i && has_C1i && has_C1j && has_next_C1j) successtest+=bp_tor_score(next_C1i,C1i,C1j,next_C1j,C1p_mu,C1p_sd, tol,weight_tor,nominator,denominator,sd_tor); // C1p: C1'[i+1]-C1'[i]-C1'[j]-C1'[j-1]
                        if (has_next_O4i && has_O4i && has_O4j && has_next_O4j) successtest+=bp_tor_score(next_O4i,O4i,O4j,next_O4j,O4p_mu,O4p_sd, tol,weight_tor,nominator,denominator,sd_tor); // O4p: O4'[i+1]-O4'[i]-O4'[j]-O4'[j-1]
                        if (has_prev_O3i && has_O3i && has_O3j && has_prev_O3j) successtest+=bp_tor_score(prev_O3i,O3i,O3j,prev_O3j,O3m_mu,O3m_sd, tol,weight_tor,nominator,denominator,sd_tor); // O3m: O3'[i-1]-O3'[i]-O3'[j]-O3'[j+1]
                    }
                    if (weight_ang>0)
                    {
                        successtest--;
                        if (has_next_Pi  && has_Pi  && has_prev_Pj  && has_Pj ) successtest+=bp_ang_score( next_Pi, Pi, prev_Pj, Pj, aPp_mu, aPp_sd, tol,weight_ang,nominator,denominator,sd_ang); //  aPp:      <P[i+1]P[i],P[j-1]P[j]>
                        if (has_next_O5i && has_O5i && has_prev_O5j && has_O5j) successtest+=bp_ang_score(next_O5i,O5i,prev_O5j,O5j,aO5p_mu,aO5p_sd, tol,weight_ang,nominator,denominator,sd_ang); // aO5p: <O5'[i+1]O5'[i],O5'[j-1]O5'[j]> 
                        if (has_next_C5i && has_C5i && has_prev_C5j && has_C5j) successtest+=bp_ang_score(next_C5i,C5i,prev_C5j,C5j,aC5p_mu,aC5p_sd, tol,weight_ang,nominator,denominator,sd_ang); // aC5p: <C5'[i+1]C5'[i],C5'[j-1]C5'[j]> 
                        if (has_next_C4i && has_C4i && has_prev_C4j && has_C4j) successtest+=bp_ang_score(next_C4i,C4i,prev_C4j,C4j,aC4p_mu,aC4p_sd, tol,weight_ang,nominator,denominator,sd_ang); // aC4p: <C4'[i+1]C4'[i],C4'[j-1]C4'[j]> 
                        if (has_next_C3i && has_C3i && has_prev_C3j && has_C3j) successtest+=bp_ang_score(next_C3i,C3i,prev_C3j,C3j,aC3p_mu,aC3p_sd, tol,weight_ang,nominator,denominator,sd_ang); // aC3p: <C3'[i+1]C3'[i],C3'[j-1]C3'[j]> 
                        if (has_next_C2i && has_C2i && has_prev_C2j && has_C2j) successtest+=bp_ang_score(next_C2i,C2i,prev_C2j,C2j,aC2p_mu,aC2p_sd, tol,weight_ang,nominator,denominator,sd_ang); // aC2p: <C2'[i+1]C2'[i],C2'[j-1]C2'[j]> 
                        if (has_next_C1i && has_C1i && has_prev_C1j && has_C1j) successtest+=bp_ang_score(next_C1i,C1i,prev_C1j,C1j,aC1p_mu,aC1p_sd, tol,weight_ang,nominator,denominator,sd_ang); // aC1p: <C1'[i+1]C1'[i],C1'[j-1]C1'[j]> 
                        if (has_next_O4i && has_O4i && has_prev_O4j && has_O4j) successtest+=bp_ang_score(next_O4i,O4i,prev_O4j,O4j,aO4p_mu,aO4p_sd, tol,weight_ang,nominator,denominator,sd_ang); // aO4p: <O4'[i+1]O4'[i],O4'[j-1]O4'[j]> 
                        if (has_next_O3i && has_O3i && has_prev_O3j && has_O3j) successtest+=bp_ang_score(next_O3i,O3i,prev_O3j,O3j,aO3p_mu,aO3p_sd, tol,weight_ang,nominator,denominator,sd_ang); // aO3p: <O3'[i+1]O3'[i],O3'[j-1]O3'[j]> 
                        if (has_C4i      && has_C1i && has_C4j      && has_C1j)//successtest+=bp_ang_score(    C4i,C1i,     C4j,C1j, aCC_mu, aCC_sd, tol,weight_ang,nominator,denominator,sd_ang); //  aCC:    <C4'[i]C1'[i],C4'[j]C1'[j]>
                        {
                            ang=rad2deg(Points4Angle(C4i,C1i,C4j,C1j));
                            if(ang>=-180)
                            {
                                successtest++;
                                denominator+=weight_ang; 
                                if (ang>aCC_mu) nominator+=weight_ang;
                                else nominator+=weight_ang*(1-fabs(ang-aCC_mu)/(tol*aCC_sd));
                            }
                        }
                    }
                    
                    if (denominator>0) bp_vec.push_back(pair<float,vector<string> >(
                        nominator/denominator+adjust1*successtest/totaltest+
                        adjust2*adjuststd(sd_len,sd_nn,sd_tor,sd_ang,
                            sd_len_mean,sd_nn_mean,sd_tor_mean,sd_ang_mean,
                            weight_len,weight_nn,weight_tor,weight_ang)+adjust3, tmp_bp));

                    sd_len.clear();
                    sd_nn.clear();
                    sd_tor.clear();
                    sd_ang.clear();
                }
            }
        }
        //cerr<<endl;
    }

    /* clean up */
    vector<string>().swap(tmp_bp);

    prev_Pi.clear();  prev_Pj.clear();
    prev_O5i.clear(); prev_O5j.clear();
    prev_C5i.clear(); prev_C5j.clear();
    prev_C4i.clear(); prev_C4j.clear();
    prev_C3i.clear(); prev_C3j.clear();
    prev_C2i.clear(); prev_C2j.clear();
    prev_C1i.clear(); prev_C1j.clear();
    prev_O4i.clear(); prev_O4j.clear();
    prev_O3i.clear(); prev_O3j.clear();
    prev_Nxi.clear(); prev_Nxj.clear();

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
    next_Nxi.clear(); next_Nxj.clear();
    return;
}

#endif

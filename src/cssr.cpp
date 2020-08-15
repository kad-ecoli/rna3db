const char* docstring=""
"cssr cg.pdb [outfmt] [atom]\n"
"    secondary structure assignment for coarse grain model\n"
"\n"
"outfmt: output format\n"
"    1 - dot bracket format\n"
"    2 - (default) DSSR format\n"
"    4 - confidence score\n"
"atom: atom used for assignment. default is all the following atoms:\n"
"    \" P  \", \" O5'\", \" O4'\", \" O3'\", \" C5'\", \" C4'\", \" C3'\", \" C2'\", \" C1'\"\n"
"    and N (N9 for a/g; N1 for c/t/u)\n"
;

#include <iostream>
#include "PDBParser.hpp"
#include "cssr.hpp"

int main(int argc,char **argv)
{
    if (argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    string pdbfile   =argv[1];
    int outfmt       =(argc>2)?atoi(argv[2]):2;
    string atom      =(argc>3)?argv[3]:"auto";
    
    bool show_dot    =(outfmt%2==1); outfmt/=2;
    bool show_dssr   =(outfmt%2==1); outfmt/=2;
    bool show_conf   =(outfmt%2==1);
    
    int atomic_detail=2;
    int allowX       =1; // only allow ATOM and MSE
    ModelUnit pdb_entry=read_pdb_structure(
        pdbfile.c_str(),atomic_detail,allowX);

    if (atom.size()!=4)
    {
        cerr<<"ERROR! atom name must have 4 characters, including spaces."<<endl;
        return 1;
    }
    else if (atom!="auto") removeUnspecifiedAtom(pdb_entry, atom);

    vector<string> res_str_vec;
    vector<pair<float,vector<string> > > bp_vec;
    cssr(pdb_entry, res_str_vec, bp_vec);

    /* output */
    size_t bp;
    vector<size_t> filtered_bp_vec;
    if (show_dot || show_dssr) filter_bp(bp_vec, filtered_bp_vec);
    if (show_dot)
    {
        vector<char> dot_bracket;
        cssr_bp2dot(res_str_vec, filtered_bp_vec, bp_vec, dot_bracket);
        for (size_t r=0;r<dot_bracket.size();r++) cout<<dot_bracket[r];
        cout<<endl;
        dot_bracket.clear();
    }
    if (show_dssr)
    {
        size_t count=0;
        cout<<"\n"
            <<"****************************************************************************"
            <<"\nList of "<<filtered_bp_vec.size()<<" base pairs\n"
            <<"     nt1            nt2            bp  name        Saenger   LW   DSSR"<<endl;
        for (bp=0;bp<bp_vec.size();bp++)
        {
            if (find(filtered_bp_vec.begin(), filtered_bp_vec.end(),bp
                )==filtered_bp_vec.end()) continue;
            cout<<setw(4)<<right<<++count<<' '
                <<setw(14)<<left<<bp_vec[bp].second[0]<<' '
                <<setw(14)<<left<<bp_vec[bp].second[1]<<' '
                <<bp_vec[bp].second[2];
            if (bp_vec[bp].second[2].substr(4)=="WC")
                cout<<"          19-XIX    cWW  cW-W";
            else if (bp_vec[bp].second[2].substr(4)=="Wobble")
                cout<<"      28-XXVIII cWW  cW-W";
            cout<<endl;
        }
    }
    vector<size_t>().swap(filtered_bp_vec);
    if (show_conf && bp_vec.size())
    {
        sort(bp_vec.begin(), bp_vec.end());
        for (bp=bp_vec.size()-1;bp>0;bp--)
            cout<<bp_vec[bp].second[0]<<'\t'<<bp_vec[bp].second[1]<<'\t'
                <<bp_vec[bp].first<<endl;
    }

    /* clean up */
    vector<string>().swap(res_str_vec);
    vector<pair<float,vector<string> > >().swap(bp_vec);
    vector<ChainUnit>().swap(pdb_entry.chains);
    return 0;
}

const char* docstring=""
"AllDist full.pdb\n"
"    calculate inter-atomic distances in pdb file\n"
"Output: missing distance are marked -1.000\n"
"        [0]   PP: P[i]-P[j]\n"
"        [1] O5O5: O5'[i]-O5'[j]\n"
"        [2] C5C5: C5'[i]-C5'[j]\n"
"        [3] C4C4: C4'[i]-C4'[j]\n"
"        [4] C3C3: C3'[i]-C3'[j]\n"
"        [5] C2C2: C2'[i]-C2'[j]\n"
"        [6] C1C1: C1'[i]-C1'[j]\n"
"        [7] O4O4: O4'[i]-O4'[j]\n"
"        [8] O3O3: O3'[i]-O3'[j]\n"
"        [9]   NN: N[i]-N[j]. N=N9 for a/g; N=N1 for c/t/u\n"
;

#include <iostream>
#include "PDBParser.hpp"
#include "AllDist.hpp"

int main(int argc,char **argv)
{
    if (argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    string pdbfile     =argv[1];
    
    int atomic_detail  =2;
    int allowX         =1; // only allow ATOM and MSE
    ModelUnit pdb_entry=read_pdb_structure(
        pdbfile.c_str(),atomic_detail,allowX);

    vector<char>   chain1_list; // Chain.ResnResi^Icode
    vector<char>   chain2_list;
    vector<string> resn1_list;
    vector<string> resn2_list;
    vector<int>    resi1_list;
    vector<int>    resi2_list;
    vector<char>   icode1_list;
    vector<char>   icode2_list;
    
    vector<vector<float> >BPlenMat;
    size_t BPcount=AllDist(chain1_list, chain2_list, resn1_list, resn2_list,
        resi1_list, resi2_list, icode1_list, icode2_list, pdb_entry,
        BPlenMat);

    cout<<"N c resi  N c resi ";
    cout<<"      PP    O5O5    C5C5    C4C4    C3C3    C2C2"
        <<"    C1C1    O4O4    O3O3      NN";
    cout<<endl;

    for (size_t bp=0;bp<BPcount;bp++)
    {
        cout<<char(tolower(resn1_list[bp][2]))<<' '<<chain1_list[bp]<<' '
            <<setw(4)<<resi1_list[bp]<<icode1_list[bp]<<' '
            <<char(tolower(resn2_list[bp][2]))<<' '<<chain2_list[bp]<<' '
            <<setw(4)<<resi2_list[bp]<<icode2_list[bp]
            <<setiosflags(ios::fixed)<<setprecision(3)
            <<' '<<setw(7)<<BPlenMat[bp][0]
            <<' '<<setw(7)<<BPlenMat[bp][1]
            <<' '<<setw(7)<<BPlenMat[bp][2]
            <<' '<<setw(7)<<BPlenMat[bp][3]
            <<' '<<setw(7)<<BPlenMat[bp][4]
            <<' '<<setw(7)<<BPlenMat[bp][5]
            <<' '<<setw(7)<<BPlenMat[bp][6]
            <<' '<<setw(7)<<BPlenMat[bp][7]
            <<' '<<setw(7)<<BPlenMat[bp][8]
            <<' '<<setw(7)<<BPlenMat[bp][9]
            <<endl;
            BPlenMat[bp].clear();
    }

    /* clean up */
    vector<char>   ().swap(chain1_list);
    vector<char>   ().swap(chain2_list);
    vector<string> ().swap(resn1_list);
    vector<string> ().swap(resn2_list);
    vector<int>    ().swap(resi1_list);
    vector<int>    ().swap(resi2_list);
    vector<char>   ().swap(icode1_list);
    vector<char>   ().swap(icode2_list);
    vector<vector<float> >().swap(BPlenMat);
    vector<ChainUnit>().swap(pdb_entry.chains);
    return 0;
}

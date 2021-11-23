const char* docstring=""
"BPtorsion2 full.pdb full.dssr [option] [bptypes]\n"
"    calculate diheral angles, distance, and angles in base pairs assigned by DSSR\n"
"option: (default is 7=1+2+4)\n"
"    1 - diheral angle. missing diherals are marked -360.00.\n"
"        dihedral *p are marked -360.00 if i-1 and j+1 unpaired\n"
"        dihedral *m are marked -360.00 if i+1 and j-1 unpaired\n"
"        [0] PCCP: P[i]-C4'[i]-C4'[j]-P[j]\n"
"        [1] CNNC: C4'[i]-N[i]-N[j]-C4'[j]\n"
"    2 - distance. missing distance are marked -1.000\n"
"        [0]   PP: P[i]-P[j]\n"
"        [1]   CC: C4'[i]-C4'[j]\n"
"        [2]   NN: N[i]-N[j]. N=N9 for a/g; N=N1 for c/t/u\n"
"    4 - angle. missing angle are marked -360.00\n"
"        angles *p are marked -360.00 if i-1 and j+1 unpaired\n"
"        angles *m are marked -360.00 if i+1 and j-1 unpaired\n"
"        [0]  aPC: <P[i]C4'[i],P[j]C4'[j]>\n"
"        [1]  aCN: <C4'[i]N[i],N[j]C4'[j]>\n"
"bptypes: (default: WC,Wobble) base pair type to consider\n"
"        it must be a comma separated list from the following base pair types\n"
"        WC,--,Wobble,Sheared,Linker,rHoogsteen,Platform,~Sheared,Imino,~rHoogsteen,rWC,Hoogsteen,~Wobble,Calcutta,~Hoogsteen,~rWobble,rWobble\n"
;

#include <iostream>
#include "PDBParser.hpp"
#include "BPtorsion2.hpp"

int main(int argc,char **argv)
{
    if (argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    string pdbfile     =argv[1];
    string dssrfile    =(argc>2)?argv[2]:"";
    int option         =(argc>3)?atoi(argv[3]):7;
    bool show_tor      =(option%2==1); option/=2;
    bool show_len      =(option%2==1); option/=2;
    bool show_ang      =(option%2==1);
    string bptypes     =(argc>4)?argv[4]:"WC,Wobble";
    
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
    vector<string> name_list;
    size_t BPcount=0;
    if (dssrfile.size()) BPcount=read_dssr_base_pair(chain1_list,chain2_list,
        resn1_list,resn2_list,resi1_list,resi2_list,icode1_list,
        icode2_list,name_list,dssrfile,bptypes);
    else BPcount=dummy_base_pair(pdb_entry,chain1_list,chain2_list,resn1_list,resn2_list,
        resi1_list,resi2_list,icode1_list,icode2_list,name_list);
    //for (size_t bp=0;bp<chain1_list.size();bp++) cout<<chain1_list[bp]<<';'
    //    <<resn1_list[bp]<<';'<<resi1_list[bp]<<';'<<icode1_list[bp]<<endl;
    
    vector<vector<float> >BPtorMat;
    vector<vector<float> >BPlenMat;
    vector<vector<float> >BPangMat;
    BPtorsion2(chain1_list, chain2_list, resn1_list, resn2_list,
        resi1_list, resi2_list, icode1_list, icode2_list, pdb_entry,
        show_tor, show_len, show_ang, BPtorMat, BPlenMat, BPangMat);

    cout<<"N c resi  N c resi  name       ";
    if (show_tor)   cout<<"    PCCP    CNNC    PNNP";
    if (show_len)   cout<<"      PP      CC      NN";
    if (show_ang)   cout<<"     aPC     aCN     aPN";
    cout<<endl;

    size_t bp;
    for (bp=0;bp<BPcount;bp++)
    {
        cout<<char(tolower(resn1_list[bp][2]))<<' '<<chain1_list[bp]<<' '
            <<setw(4)<<resi1_list[bp]<<icode1_list[bp]<<' '
            <<char(tolower(resn2_list[bp][2]))<<' '<<chain2_list[bp]<<' '
            <<setw(4)<<resi2_list[bp]<<icode2_list[bp]<<' '
            <<name_list[bp];
        if (show_tor)
        {
            cout<<setiosflags(ios::fixed)<<setprecision(2)
                <<' '<<setw(7)<<BPtorMat[bp][0]
                <<' '<<setw(7)<<BPtorMat[bp][1]
                <<' '<<setw(7)<<BPtorMat[bp][2];
            BPtorMat[bp].clear();
        }
        if (show_len)
        {
            cout<<setiosflags(ios::fixed)<<setprecision(3)
                <<' '<<setw(7)<<BPlenMat[bp][0]
                <<' '<<setw(7)<<BPlenMat[bp][1]
                <<' '<<setw(7)<<BPlenMat[bp][2];
                BPlenMat[bp].clear();
        }
        if (show_ang)
        {
            cout<<setiosflags(ios::fixed)<<setprecision(2)
                <<' '<<setw(7)<<BPangMat[bp][0]
                <<' '<<setw(7)<<BPangMat[bp][1]
                <<' '<<setw(7)<<BPangMat[bp][2];
                BPangMat[bp].clear();
        }
        cout<<endl;
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
    vector<string> ().swap(name_list);
    vector<vector<float> >().swap(BPtorMat);
    vector<vector<float> >().swap(BPlenMat);
    vector<vector<float> >().swap(BPangMat);
    vector<ChainUnit>().swap(pdb_entry.chains);
    return 0;
}

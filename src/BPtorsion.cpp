const char* docstring=""
"BPtorsion full.pdb full.dssr [option] [bptypes]\n"
"    calculate diheral angles, distance, and angles in base pairs assigned by DSSR\n"
"option:\n"
"    1 - diheral angle. missing diherals are marked -360.00.\n"
"        dihedral *p are marked -360.00 if i-1 and j+1 unpaired\n"
"        dihedral *m are marked -360.00 if i+1 and j-1 unpaired\n"
"        [0]  Pm:  P[i-1]-P[i]-P[j]-P[j-1]\n"
"        [1]  Pp:  P[i+1]-P[i]-P[j]-P[j+1]\n"
"        [2]  O5m: O5'[i-1]-O5'[i]-O5'[j]-O5'[j-1]\n"
"        [3]  O5p: O5'[i+1]-O5'[i]-O5'[j]-O5'[j+1]\n"
"        [4]  C5m: C5'[i-1]-C5'[i]-C5'[j]-C5'[j-1]\n"
"        [5]  C5p: C5'[i+1]-C5'[i]-C5'[j]-C5'[j+1]\n"
"        [6]  C4m: C4'[i-1]-C4'[i]-C4'[j]-C4'[j-1]\n"
"        [7]  C4p: C4'[i+1]-C4'[i]-C4'[j]-C4'[j+1]\n"
"        [8]  C3m: C3'[i-1]-C3'[i]-C3'[j]-C3'[j-1]\n"
"        [9]  C3p: C3'[i+1]-C3'[i]-C3'[j]-C3'[j+1]\n"
"        [10] C2m: C2'[i-1]-C2'[i]-C2'[j]-C2'[j-1]\n"
"        [11] C2p: C2'[i+1]-C2'[i]-C2'[j]-C2'[j+1]\n"
"        [12] C1m: C1'[i-1]-C1'[i]-C1'[j]-C1'[j-1]\n"
"        [13] C1p: C1'[i+1]-C1'[i]-C1'[j]-C1'[j+1]\n"
"        [14] O4m: O4'[i-1]-O4'[i]-O4'[j]-O4'[j-1]\n"
"        [15] O4p: O4'[i+1]-O4'[i]-O4'[j]-O4'[j+1]\n"
"        [16] O3m: O3'[i-1]-O3'[i]-O3'[j]-O3'[j-1]\n"
"        [17] O3p: O3'[i+1]-O3'[i]-O3'[j]-O3'[j+1]\n"
"        [18] P44P: P[i]-C4'[i]-C4'[j]-P[j]\n"
"        [19] C4114C: C4'[i]-C1'[i]-C1'[j]-C4'[j]\n"
"    2 - (default) distance. missing distance are marked -1.000\n"
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
"    4 - angle. missing angle are marked -360.00\n"
"        angles *p are marked -360.00 if i-1 and j+1 unpaired\n"
"        angles *m are marked -360.00 if i+1 and j-1 unpaired\n"
"        [0]  aPm:      <P[i-1]P[i],P[j+1]P[j]>\n"
"        [1]  aPp:      <P[i+1]P[i],P[j-1]P[j]>\n"
"        [2]  aO5m: <O5'[i-1]O5'[i],O5'[j+1]O5'[j]>\n"
"        [3]  aO5p: <O5'[i+1]O5'[i],O5'[j-1]O5'[j]>\n"
"        [4]  aC5m: <C5'[i-1]C5'[i],C5'[j+1]C5'[j]>\n"
"        [5]  aC5p: <C5'[i+1]C5'[i],C5'[j-1]C5'[j]>\n"
"        [6]  aC4m: <C4'[i-1]C4'[i],C4'[j+1]C4'[j]>\n"
"        [7]  aC4p: <C4'[i+1]C4'[i],C4'[j-1]C4'[j]>\n"
"        [8]  aC3m: <C3'[i-1]C3'[i],C3'[j+1]C3'[j]>\n"
"        [9]  aC3p: <C3'[i+1]C3'[i],C3'[j-1]C3'[j]>\n"
"        [10] aC2m: <C2'[i-1]C2'[i],C2'[j+1]C2'[j]>\n"
"        [11] aC2p: <C2'[i+1]C2'[i],C2'[j-1]C2'[j]>\n"
"        [12] aC1m: <C1'[i-1]C1'[i],C1'[j+1]C1'[j]>\n"
"        [13] aC1p: <C1'[i+1]C1'[i],C1'[j-1]C1'[j]>\n"
"        [14] aO4m: <O4'[i-1]O4'[i],O4'[j+1]O4'[j]>\n"
"        [15] aO4p: <O4'[i+1]O4'[i],O4'[j-1]O4'[j]>\n"
"        [16] aO3m: <O3'[i-1]O3'[i],O3'[j+1]O3'[j]>\n"
"        [17] aO3p: <O3'[i+1]O3'[i],O3'[j-1]O3'[j]>\n"
"        [18] aPC:      <P[i]C4'[i],P[j]C4'[j]>\n"
"        [19] aCC:    <C4'[i]C1'[i],C4'[j]C1'[j]>\n"
"bptypes: (default: WC,Wobble) base pair type to consider\n"
"        it must be a comma separated list from the following base pair types\n"
"        WC,--,Wobble,Sheared,Linker,rHoogsteen,Platform,~Sheared,Imino,~rHoogsteen,rWC,Hoogsteen,~Wobble,Calcutta,~Hoogsteen,~rWobble,rWobble\n"
;

#include <iostream>
#include "PDBParser.hpp"
#include "BPtorsion.hpp"

int main(int argc,char **argv)
{
    if (argc<3)
    {
        cerr<<docstring;
        return 0;
    }
    string pdbfile     =argv[1];
    string dssrfile    =argv[2];
    int option         =(argc>3)?atoi(argv[3]):2;
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
    size_t BPcount=read_dssr_base_pair(chain1_list,chain2_list,
        resn1_list,resn2_list,resi1_list,resi2_list,icode1_list,
        icode2_list,name_list,dssrfile,bptypes);
    
    vector<vector<float> >BPtorMat;
    vector<vector<float> >BPlenMat;
    vector<vector<float> >BPangMat;
    BPtorsion(chain1_list, chain2_list, resn1_list, resn2_list,
        resi1_list, resi2_list, icode1_list, icode2_list, pdb_entry,
        show_tor, show_len, show_ang, BPtorMat, BPlenMat, BPangMat);

    cout<<"N c resi  N c resi  name       ";
    if (show_tor)   cout<<"      Pm      Pp     O5m     O5p     C5m     C5p"
        <<"     C4m     C4p     C3m     C3p     C2m     C2p     C1m     C1p"
        <<"     O4m     O4p     O3m     O3p    P44P  C4114C";
    if (show_len)   cout<<"      PP    O5O5    C5C5    C4C4    C3C3    C2C2"
        <<"    C1C1    O4O4    O3O3      NN";
    if (show_ang)   cout<<"     aPm     aPp    aO5m    aO5p    aC5m    aC5p"
        <<"    aC4m    aC4p    aC3m    aC3p    aC2m    aC2p    aC1m    aC1p"
        <<"    aO4m    aO4p    aO3m    aO3p     aPC     aCC";
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
                <<' '<<setw(7)<<BPtorMat[bp][2]
                <<' '<<setw(7)<<BPtorMat[bp][3]
                <<' '<<setw(7)<<BPtorMat[bp][4]
                <<' '<<setw(7)<<BPtorMat[bp][5]
                <<' '<<setw(7)<<BPtorMat[bp][6]
                <<' '<<setw(7)<<BPtorMat[bp][7]
                <<' '<<setw(7)<<BPtorMat[bp][8]
                <<' '<<setw(7)<<BPtorMat[bp][9]
                <<' '<<setw(7)<<BPtorMat[bp][10]
                <<' '<<setw(7)<<BPtorMat[bp][11]
                <<' '<<setw(7)<<BPtorMat[bp][12]
                <<' '<<setw(7)<<BPtorMat[bp][13]
                <<' '<<setw(7)<<BPtorMat[bp][14]
                <<' '<<setw(7)<<BPtorMat[bp][15]
                <<' '<<setw(7)<<BPtorMat[bp][16]
                <<' '<<setw(7)<<BPtorMat[bp][17]
                <<' '<<setw(7)<<BPtorMat[bp][18]
                <<' '<<setw(7)<<BPtorMat[bp][19];
            BPtorMat[bp].clear();
        }
        if (show_len)
        {
            cout<<setiosflags(ios::fixed)<<setprecision(3)
                <<' '<<setw(7)<<BPlenMat[bp][0]
                <<' '<<setw(7)<<BPlenMat[bp][1]
                <<' '<<setw(7)<<BPlenMat[bp][2]
                <<' '<<setw(7)<<BPlenMat[bp][3]
                <<' '<<setw(7)<<BPlenMat[bp][4]
                <<' '<<setw(7)<<BPlenMat[bp][5]
                <<' '<<setw(7)<<BPlenMat[bp][6]
                <<' '<<setw(7)<<BPlenMat[bp][7]
                <<' '<<setw(7)<<BPlenMat[bp][8]
                <<' '<<setw(7)<<BPlenMat[bp][9];
                BPlenMat[bp].clear();
        }
        if (show_ang)
        {
            cout<<setiosflags(ios::fixed)<<setprecision(2)
                <<' '<<setw(7)<<BPangMat[bp][0]
                <<' '<<setw(7)<<BPangMat[bp][1]
                <<' '<<setw(7)<<BPangMat[bp][2]
                <<' '<<setw(7)<<BPangMat[bp][3]
                <<' '<<setw(7)<<BPangMat[bp][4]
                <<' '<<setw(7)<<BPangMat[bp][5]
                <<' '<<setw(7)<<BPangMat[bp][6]
                <<' '<<setw(7)<<BPangMat[bp][7]
                <<' '<<setw(7)<<BPangMat[bp][8]
                <<' '<<setw(7)<<BPangMat[bp][9]
                <<' '<<setw(7)<<BPangMat[bp][10]
                <<' '<<setw(7)<<BPangMat[bp][11]
                <<' '<<setw(7)<<BPangMat[bp][12]
                <<' '<<setw(7)<<BPangMat[bp][13]
                <<' '<<setw(7)<<BPangMat[bp][14]
                <<' '<<setw(7)<<BPangMat[bp][15]
                <<' '<<setw(7)<<BPangMat[bp][16]
                <<' '<<setw(7)<<BPangMat[bp][17]
                <<' '<<setw(7)<<BPangMat[bp][18]
                <<' '<<setw(7)<<BPangMat[bp][19];
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

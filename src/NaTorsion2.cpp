const char* docstring=""
"NaTorsion2 full.pdb [option]\n"
"    calculate torsion angles, bond length and bond angles\n"
"option:\n"
"    1 - (default) torsion angle. missing torsions are marked -360.00.\n"
"        [0] pseudo torsion eta:     C4'[-1]-P-C4'-P[+1]\n"
"        [1] pseudo torsion theta:   P-C4'-P[+1]-C4'[+1]\n"
"        [2] pseudo torsion PPC4C1:  P[+1]-P-C4'-C1'\n"
"        [3] pseudo torsion PPC4N:   P[+1]-P-C4'-N\n"
"        [4] pseudo torsion C4PC4C1: C4'[-1]-P-C4'-C1'\n"
"        [5] pseudo torsion C4PC4N:  C4'[-1]-P-C4'-N\n"
"                                    N=N9,C=C4 for a/g; N=N1,C=C2 for c/t/u\n"
"    2 - bond length. missing bond are marked -1.000\n"
"        [0] pseudo bond: P-C4'\n"
"        [1] pseudo bond: C4'-P[+1]\n"
"        [2] pseudo bond: C4'-C1'\n"
"        [3] pseudo bond: C4'-N\n"
"        [4] pseudo bond: P-P[+1]\n"
"        [5] pseudo bond: C4'[-1]-C4'\n"
"    4 - bond angle. missing angle are marked -360.00\n"
"        [0] pseudo angle: C4'[-1]-P-C4'\n"
"        [1] pseudo angle: P-C4'-P[+1]\n"
"        [2] pseudo angle: P-C4'-C1'\n"
"        [3] pseudo angle: P-C4'-N\n"
;

#include <iostream>
#include "PDBParser.hpp"
#include "NaTorsion2.hpp"

int main(int argc,char **argv)
{
    if (argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    string filename    =argv[1];
    int option         =(argc>2)?atoi(argv[2]):1;
    bool show_tor      =(option%2==1); option/=2;
    bool show_len      =(option%2==1); option/=2;
    bool show_ang      =(option%2==1);
    
    int atomic_detail  =2;
    int allowX         =1; // only allow ATOM and MSE
    ModelUnit pdb_entry=read_pdb_structure(
        filename.c_str(),atomic_detail,allowX);

    int c,r; // chain index, residue index;
    vector<vector<float> >NaTorMat;
    vector<vector<float> >NaLenMat;
    vector<vector<float> >NaAngMat;
    cout<<"N c resi ";
    if (show_tor)   cout<<"     eta   theta  PPC4C1   PPC4N C4PC4C1  C4PC4N";
    if (show_len)   cout<<"    PC4'   C4'Pp  C4'C1'    C4'N     PPp C4'mC4'";
    if (show_ang)   cout<<" C4'mPC4'   PC4'Pp  PC4'C1'    PC4'N";
    cout<<endl;
    for (c=0;c<pdb_entry.chains.size();c++)
    {
        NaTorsion2(pdb_entry.chains[c],NaTorMat,NaLenMat,NaAngMat,
            show_tor,show_len,show_ang);
        for (r=0;r<pdb_entry.chains[c].residues.size();r++)
        {
            cout<<char(tolower(pdb_entry.chains[c].residues[r].resn[2]))<<' '
                <<pdb_entry.chains[c].chainID<<' '
                <<setw(4)<<pdb_entry.chains[c].residues[r].resi
                <<pdb_entry.chains[c].residues[r].icode;
            if (show_tor)
            {
                cout<<setiosflags(ios::fixed)<<setprecision(2)
                    <<' '<<setw(7)<<NaTorMat[r][0]
                    <<' '<<setw(7)<<NaTorMat[r][1]
                    <<' '<<setw(7)<<NaTorMat[r][2]
                    <<' '<<setw(7)<<NaTorMat[r][3]
                    <<' '<<setw(7)<<NaTorMat[r][4]
                    <<' '<<setw(7)<<NaTorMat[r][5];
                NaTorMat[r].clear();
            }
            if (show_len)
            {
                cout<<setiosflags(ios::fixed)<<setprecision(3)
                    <<' '<<setw(7)<<NaLenMat[r][0]
                    <<' '<<setw(7)<<NaLenMat[r][1]
                    <<' '<<setw(7)<<NaLenMat[r][2]
                    <<' '<<setw(7)<<NaLenMat[r][3]
                    <<' '<<setw(7)<<NaLenMat[r][4]
                    <<' '<<setw(7)<<NaLenMat[r][5];
                    NaLenMat[r].clear();
            }
            if (show_ang)
            {
                cout<<setiosflags(ios::fixed)<<setprecision(2)
                    <<' '<<setw(8)<<NaAngMat[r][0]
                    <<' '<<setw(8)<<NaAngMat[r][1]
                    <<' '<<setw(8)<<NaAngMat[r][2]
                    <<' '<<setw(8)<<NaAngMat[r][3];
                    NaAngMat[r].clear();
            }
            cout<<endl;
        }
        if (show_tor) NaTorMat.clear();
        if (show_len) NaLenMat.clear();
        if (show_ang) NaAngMat.clear();
    }
    return 0;
}

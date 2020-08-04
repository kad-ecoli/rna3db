const char* docstring="NaTorsion full.pdb\n"
"    Calculate pucker (C1'C2'C3'C4': >0 for C3'-endo; <0 for C2'-endo),\n"
"    pseudo torsion (eta, theta, eta', theta'),\n"
"    backbone torsion (alpha, beta, gamma, delta , epsilon, zeta) and\n"
"    sidechain torsion (chi) from full.pdb.\n"
"    Torsions for missing atoms are marked with '-360.00'.\n";

#include <iostream>
#include "PDBParser.hpp"
#include "NaTorsion.hpp"

int main(int argc,char **argv)
{
    if (argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    
    int atomic_detail=2;
    int allowX=1;        // only allow ATOM and MSE
    ModelUnit pdb_entry=read_pdb_structure(argv[1],atomic_detail,allowX);

    int c,r; // chain index, residue index;
    cout<<"N c resi   pucker     eta   theta    eta'  theta'   alpha    beta   gamma   delta epsilon    zeta     chi"<<endl;
    for (c=0;c<pdb_entry.chains.size();c++)
    {
        vector<vector<float> >NaTorMat=NaTorsion(pdb_entry.chains[c]);
        for (r=0;r<pdb_entry.chains[c].residues.size();r++)
        {
            cout<<char(tolower(pdb_entry.chains[c].residues[r].resn[2]))<<' '
                <<pdb_entry.chains[c].chainID<<' '
                <<setw(4)<<pdb_entry.chains[c].residues[r].resi
                <<pdb_entry.chains[c].residues[r].icode<<' '
                <<setiosflags(ios::fixed)<<setprecision(2)
                <<setw(7)<<NaTorMat[r][0]<<' '
                <<setw(7)<<NaTorMat[r][1]<<' '
                <<setw(7)<<NaTorMat[r][2]<<' '
                <<setw(7)<<NaTorMat[r][3]<<' '
                <<setw(7)<<NaTorMat[r][4]<<' '
                <<setw(7)<<NaTorMat[r][5]<<' '
                <<setw(7)<<NaTorMat[r][6]<<' '
                <<setw(7)<<NaTorMat[r][7]<<' '
                <<setw(7)<<NaTorMat[r][8]<<' '
                <<setw(7)<<NaTorMat[r][9]<<' '
                <<setw(7)<<NaTorMat[r][10]<<' '
                <<setw(7)<<NaTorMat[r][11]<<endl;
        }
        NaTorMat.clear();
    }
    return 0;
}

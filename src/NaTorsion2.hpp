/* calculate backbone torsion angles (omega, psi, phi, kappa, alpha) */
#ifndef NaTorsion2_HPP
#define NaTorsion2_HPP 1
#include <cstring>
#include "PDBParser.hpp"
#include "GeometryTools.hpp"

using namespace std;


void NaTorsion2(ChainUnit& chain, vector<vector<float> >&NaTorMat,
    vector<vector<float> >&NaLenMat, vector<vector<float> >&NaAngMat, 
    const bool show_tor, const bool show_len, const bool show_ang)
{
    int L=chain.residues.size();
    // default torsion angles: omega=phi=psi=360
    vector<float> tmp_tor(6,-360.);
    vector<float> tmp_len(6,-1.);
    vector<float> tmp_ang(4,-360.);
    if (show_tor) NaTorMat.assign(L,tmp_tor);
    if (show_len) NaLenMat.assign(L,tmp_len);
    if (show_ang) NaAngMat.assign(L,tmp_ang);

    // coordinates of previous residue
    vector<float> prev_P(3,0);  bool has_prev_P=false;
    vector<float> prev_C4(3,0); bool has_prev_C4=false;
    // coordinates of current residue
    vector<float> P(3,0.);      bool has_P=false;
    vector<float> C4(3,0.);     bool has_C4=false;
    vector<float> C1(3,0.);     bool has_C1=false;
    vector<float> Nx(3,0.);     bool has_Nx=false;
    // coordinates of next residue
    vector<float> next_P(3,0);  bool has_next_P=false;
    vector<float> next_C4(3,0); bool has_next_C4=false;
    char base=' ';

    int r,a; // residue index, atom index
    
    for (r=0;r<L;r++) 
    {
        // whether the required atom exists
        has_prev_C4=false;
        has_P=false;
        has_C4=false;
        has_C1=false;
        has_Nx=false;
        has_next_P=false;
        has_next_C4=false;
        base=tolower(chain.residues[r].resn[2]);
        
        if (r>0) // find previous residue atoms
        {
            for (a=0;a<chain.residues[r-1].atoms.size();a++)
            {
                if (chain.residues[r-1].atoms[a].name==" P  ")
                {
                    has_prev_P=true;
                    prev_P=chain.residues[r-1].atoms[a].xyz;
                }
                else if (chain.residues[r-1].atoms[a].name==" C4'")
                {
                    has_prev_C4=true;
                    prev_C4=chain.residues[r-1].atoms[a].xyz;
                }
            }
        }

        // find current residue atoms
        for (a=0;a<chain.residues[r].atoms.size();a++)
        {
            if (chain.residues[r].atoms[a].name==" P  ")
            {
                has_P=true;
                P=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" C4'")
            {
                has_C4=true;
                C4=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" C1'")
            {
                has_C1=true;
                C1=chain.residues[r].atoms[a].xyz;
            }
            else if((chain.residues[r].atoms[a].name==" N9 " &&
                    (base=='a' || base=='g')) ||
                    (chain.residues[r].atoms[a].name==" N1 " && 
                    (base=='c' || base=='t' || base=='u')))
            {
                has_Nx=true;
                Nx=chain.residues[r].atoms[a].xyz;
            }
        }

        if (r+1<L) // find next residue atoms
        {
            for (a=0;a<chain.residues[r+1].atoms.size();a++)
            {
                if (chain.residues[r+1].atoms[a].name==" P  ")
                {
                    has_next_P=true;
                    next_P=chain.residues[r+1].atoms[a].xyz;
                }
                else if (chain.residues[r+1].atoms[a].name==" C4'")
                {
                    has_next_C4=true;
                    next_C4=chain.residues[r+1].atoms[a].xyz;
                }
            }
        }
 
        if (show_tor)
        {
            if (has_prev_C4 && has_P && has_C4 && has_next_P) 
                NaTorMat[r][0]=rad2deg(Points2Dihedral(prev_C4,P,C4,next_P));// eta
            if (has_P && has_C4 && has_next_P && has_next_C4) 
                NaTorMat[r][1]=rad2deg(Points2Dihedral(P,C4,next_P,next_C4));// theta
            if (has_next_P && has_P && has_C4 && has_C1) 
                NaTorMat[r][2]=rad2deg(Points2Dihedral(next_P,P,C4,C1));     // P[+1]-P-C4'-C1'
            if (has_next_P && has_P && has_C4 && has_Nx) 
                NaTorMat[r][3]=rad2deg(Points2Dihedral(next_P,P,C4,Nx));     // P[+1]-P-C4'-N
            if (has_prev_C4 && has_P && has_C4 && has_C1)     
                NaTorMat[r][4]=rad2deg(Points2Dihedral(prev_C4,P,C4,C1));    // C4'[-1]-P-C4'-C1'
            if (has_prev_C4 && has_P && has_C4 && has_Nx)    
                NaTorMat[r][5]=rad2deg(Points2Dihedral(prev_C4,P,C4,Nx));    // C4'[-1]-P-C4'-N
        }
        if (show_len)
        {
            if (has_P       && has_C4    ) NaLenMat[r][0] =Points2Distance(P,C4);      // P-C4'
            if (has_C4      && has_next_P) NaLenMat[r][1] =Points2Distance(C4,next_P); // C4'-P[+1]
            if (has_C4      && has_C1    ) NaLenMat[r][2] =Points2Distance(C4,C1);     // C4'-C1'
            if (has_C4      && has_Nx    ) NaLenMat[r][3] =Points2Distance(C4,Nx);     // C4'-N
            if (has_P       && has_next_P) NaLenMat[r][4] =Points2Distance(P,next_P);  // P-P[+1]
            if (has_prev_C4 && has_C4    ) NaLenMat[r][5] =Points2Distance(prev_C4,C4);// C4'[-1]-C4'
        }
        if (show_ang)
        {
            if (has_prev_C4 && has_P  && has_C4    ) 
                NaAngMat[r][0] =rad2deg(Points2Angle(prev_C4,P,C4));// C4'[-1]-P-C4'
            if (has_P       && has_C4 && has_next_P) 
                NaAngMat[r][1] =rad2deg(Points2Angle(P,C4,next_P)); // P-C4'-P[+1]
            if (has_P       && has_C4 && has_C1    )
                NaAngMat[r][2] =rad2deg(Points2Angle(P,C4,C1));     // P-C4'-C1'
            if (has_P       && has_C4 && has_Nx    )
                NaAngMat[r][3] =rad2deg(Points2Angle(P,C4,Nx));     // P-C4'-N
        }
    }
    

    tmp_tor.clear();
    tmp_len.clear();
    tmp_ang.clear();
    prev_C4.clear();
    P.clear();
    C4.clear();
    C1.clear();
    Nx.clear();
    next_P.clear(); 
    next_C4.clear();
    return;
}

#endif

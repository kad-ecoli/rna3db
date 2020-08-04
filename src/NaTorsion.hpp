/* calculate backbone torsion angles (omega, psi, phi, kappa, alpha) */
#ifndef BackboneTorsion_HPP
#define BackboneTorsion_HPP 1
#include <cstring>
#include "PDBParser.hpp"
#include "GeometryTools.hpp"

using namespace std;


/* calculate pucker (C1'C2'C3'C4'), pseudo torsion (eta, theta, eta', theta')
 * backbone torsion (alpha, beta, gamma, delta , epsilon, zeta) and
 * sidechain torsion (chi) from backbone.pdb */
vector<vector<float> > NaTorsion(ChainUnit& chain)
{
    int L=chain.residues.size();
    // default torsion angles: omega=phi=psi=360
    vector<float> tmp_array(12,-360.);
    vector<vector<float> > NaTorMat(L,tmp_array);

    // coordinates of previous residue
    vector<float> prev_C4(3,0); bool has_prev_C4=false;
    vector<float> prev_C1(3,0); bool has_prev_C1=false;
    vector<float> prev_O3(3,0); bool has_prev_O3=false;
    // coordinates of current residue
    vector<float> P(3,0.);      bool has_P=false;
    vector<float> O5(3,0.);     bool has_O5=false;
    vector<float> C5(3,0.);     bool has_C5=false;
    vector<float> C4(3,0.);     bool has_C4=false;
    vector<float> C3(3,0.);     bool has_C3=false;
    vector<float> C2(3,0.);     bool has_C2=false;
    vector<float> C1(3,0.);     bool has_C1=false;
    vector<float> O4(3,0.);     bool has_O4=false;
    vector<float> O3(3,0.);     bool has_O3=false;
    vector<float> Nx(3,0.);     bool has_Nx=false;
    vector<float> Cx(3,0.);     bool has_Cx=false;
    // coordinates of next residue
    vector<float> next_P(3,0);  bool has_next_P=false;
    vector<float> next_O5(3,0); bool has_next_O5=false;
    vector<float> next_C4(3,0); bool has_next_C4=false;
    vector<float> next_C1(3,0); bool has_next_C1=false;
    char base=' ';

    int r,a; // residue index, atom index
    
    for (r=0;r<L;r++) 
    {
        // whether the required atom exists
        has_prev_C4=false;
        has_prev_C1=false;
        has_prev_O3=false;
        has_P=false;
        has_O5=false;
        has_C5=false;
        has_C4=false;
        has_C3=false;
        has_C2=false;
        has_C1=false;
        has_O4=false;
        has_O3=false;
        has_Nx=false;
        has_Cx=false;
        has_next_P=false;
        has_next_O5=false;
        has_next_C4=false;
        has_next_C1=false;
        base=tolower(chain.residues[r].resn[2]);
        
        if (r>0) // find previous residue atoms
        {
            for (a=0;a<chain.residues[r-1].atoms.size();a++)
            {
                if (chain.residues[r-1].atoms[a].name==" C4'")
                {
                    has_prev_C4=true;
                    prev_C4=chain.residues[r-1].atoms[a].xyz;
                }
                else if (chain.residues[r-1].atoms[a].name==" C1'")
                {
                    has_prev_C1=true;
                    prev_C1=chain.residues[r-1].atoms[a].xyz;
                }
                else if (chain.residues[r-1].atoms[a].name==" O3'")
                {
                    has_prev_O3=true;
                    prev_O3=chain.residues[r-1].atoms[a].xyz;
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
            else if (chain.residues[r].atoms[a].name==" O5'")
            {
                has_O5=true;
                O5=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" C5'")
            {
                has_C5=true;
                C5=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" C4'")
            {
                has_C4=true;
                C4=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" C3'")
            {
                has_C3=true;
                C3=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" C2'")
            {
                has_C2=true;
                C2=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" C1'")
            {
                has_C1=true;
                C1=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" O4'")
            {
                has_O4=true;
                O4=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" O3'")
            {
                has_O3=true;
                O3=chain.residues[r].atoms[a].xyz;
            }
            else if((chain.residues[r].atoms[a].name==" N9 " &&
                    (base=='a' || base=='g')) ||
                    (chain.residues[r].atoms[a].name==" N1 " && 
                    (base=='c' || base=='t' || base=='u')))
            {
                has_Nx=true;
                Nx=chain.residues[r].atoms[a].xyz;
            }
            else if((chain.residues[r].atoms[a].name==" C4 " &&
                    (base=='a' || base=='g')) ||
                    (chain.residues[r].atoms[a].name==" C2 " && 
                    (base=='c' || base=='t' || base=='u')))
            {
                has_Cx=true;
                Cx=chain.residues[r].atoms[a].xyz;
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
                else if (chain.residues[r+1].atoms[a].name==" O5'")
                {
                    has_next_O5=true;
                    next_O5=chain.residues[r+1].atoms[a].xyz;
                }
                else if (chain.residues[r+1].atoms[a].name==" C4'")
                {
                    has_next_C4=true;
                    next_C4=chain.residues[r+1].atoms[a].xyz;
                }
                else if (chain.residues[r+1].atoms[a].name==" C1'")
                {
                    has_next_C1=true;
                    next_C1=chain.residues[r+1].atoms[a].xyz;
                }
            }
        }
 
        if (has_C1 && has_C2 && has_C3 && has_C4)
            NaTorMat[r][0]=rad2deg(Points2Dihedral(C1,C2,C3,C4));          // pucker
        
        if (has_prev_C4 && has_P && has_C4 && has_next_P)
            NaTorMat[r][1]=rad2deg(Points2Dihedral(prev_C4,P,C4,next_P));  // eta
        if (has_P && has_C4 && has_next_P && has_next_C4)
            NaTorMat[r][2]=rad2deg(Points2Dihedral(P,C4,next_P,next_C4));  // theta

        if (has_prev_C1 && has_P && has_C1 && has_next_P)
            NaTorMat[r][3]=rad2deg(Points2Dihedral(prev_C1,P,C1,next_P));  // eta'
        if (has_P && has_C1 && has_next_P && has_next_C1)
            NaTorMat[r][4]=rad2deg(Points2Dihedral(P,C1,next_P,next_C1));  // theta'
        
        if (has_prev_O3 && has_P && has_O5 && has_C5)
            NaTorMat[r][5]=rad2deg(Points2Dihedral(prev_O3,P,O5,C5));      // alpha
        if (has_P && has_O5 && has_C5 && has_C4)
            NaTorMat[r][6]=rad2deg(Points2Dihedral(P,O5,C5,C4));           // beta
        if (has_O5 && has_C5 && has_C4 && has_C3)
            NaTorMat[r][7]=rad2deg(Points2Dihedral(O5,C5,C4,C3));          // gamma
        if (has_C5 && has_C4 && has_C3 && has_O3)
            NaTorMat[r][8]=rad2deg(Points2Dihedral(C5,C4,C3,O3));          // delta
        if (has_C4 && has_C3 && has_O3 && has_next_P)
            NaTorMat[r][9]=rad2deg(Points2Dihedral(C4,C3,O3,next_P));      // epsilon
        if (has_C3 && has_O3 && has_next_P && has_next_O5)
            NaTorMat[r][10]=rad2deg(Points2Dihedral(C3,O3,next_P,next_O5));// zeta
        
        if (has_O4 && has_C1 && has_Nx && has_Cx)
            NaTorMat[r][11]=rad2deg(Points2Dihedral(O4,C1,Nx,Cx));         // chi
    }
    

    tmp_array.clear();
    prev_C4.clear();
    prev_C1.clear();
    prev_O3.clear();
    P.clear();
    O5.clear();
    C5.clear();
    C4.clear();
    C3.clear();
    C2.clear();
    C1.clear();
    O4.clear();
    O3.clear();
    Nx.clear();
    Cx.clear();
    next_P.clear(); 
    next_O5.clear();
    next_C4.clear();
    next_C1.clear();
    return NaTorMat;
}

#endif

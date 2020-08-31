#include <string>
#include <vector>
#include <cmath>
#include <cstdlib> 
#include <ctime> 

#include "PDBParser.hpp"
#include "Superpose.hpp"
#include "GeometryTools.hpp"
#include "pseudo_torsion_prob.hpp"

using namespace std;

void parse_simple_rna(map<string, map<string,vector<float> > >&simple_rna,
    map<string,vector<string> >&missing_rna)
{
    vector<string> missing_atom_vec;
    
    //missing_atom_vec.push_back(" P  ");
    missing_atom_vec.push_back(" OP1");
    missing_atom_vec.push_back(" OP2");
    missing_atom_vec.push_back(" O5'");
    missing_atom_vec.push_back(" C5'");
    //missing_atom_vec.push_back(" C4'");
    missing_atom_vec.push_back(" O4'");
    missing_atom_vec.push_back(" C3'");
    missing_atom_vec.push_back(" O3'");
    missing_atom_vec.push_back(" C2'");
    missing_atom_vec.push_back(" O2'");
    //missing_atom_vec.push_back(" C1'");
    missing_atom_vec.push_back(" N9 ");
    missing_atom_vec.push_back(" C8 ");
    missing_atom_vec.push_back(" N7 ");
    missing_atom_vec.push_back(" C5 ");
    missing_atom_vec.push_back(" C6 ");
    missing_atom_vec.push_back(" N6 ");
    missing_atom_vec.push_back(" N1 ");
    missing_atom_vec.push_back(" C2 ");
    missing_atom_vec.push_back(" N3 ");
    missing_atom_vec.push_back(" C4 ");
    missing_rna["  A"]=missing_atom_vec;
    for (size_t a=0;a<missing_atom_vec.size();a++) missing_atom_vec[a].clear();
    missing_atom_vec.clear();

    //missing_atom_vec.push_back(" P  ");
    missing_atom_vec.push_back(" OP1");
    missing_atom_vec.push_back(" OP2");
    missing_atom_vec.push_back(" O5'");
    missing_atom_vec.push_back(" C5'");
    //missing_atom_vec.push_back(" C4'");
    missing_atom_vec.push_back(" O4'");
    missing_atom_vec.push_back(" C3'");
    missing_atom_vec.push_back(" O3'");
    missing_atom_vec.push_back(" C2'");
    missing_atom_vec.push_back(" O2'");
    //missing_atom_vec.push_back(" C1'");
    missing_atom_vec.push_back(" N1 ");
    missing_atom_vec.push_back(" C2 ");
    missing_atom_vec.push_back(" O2 ");
    missing_atom_vec.push_back(" N3 ");
    missing_atom_vec.push_back(" C4 ");
    missing_atom_vec.push_back(" N4 ");
    missing_atom_vec.push_back(" C5 ");
    missing_atom_vec.push_back(" C6 ");
    missing_rna["  C"]=missing_atom_vec;
    for (size_t a=0;a<missing_atom_vec.size();a++) missing_atom_vec[a].clear();
    missing_atom_vec.clear();

    //missing_atom_vec.push_back(" P  ");
    missing_atom_vec.push_back(" OP1");
    missing_atom_vec.push_back(" OP2");
    missing_atom_vec.push_back(" O5'");
    missing_atom_vec.push_back(" C5'");
    //missing_atom_vec.push_back(" C4'");
    missing_atom_vec.push_back(" O4'");
    missing_atom_vec.push_back(" C3'");
    missing_atom_vec.push_back(" O3'");
    missing_atom_vec.push_back(" C2'");
    missing_atom_vec.push_back(" O2'");
    //missing_atom_vec.push_back(" C1'");
    missing_atom_vec.push_back(" N9 ");
    missing_atom_vec.push_back(" C8 ");
    missing_atom_vec.push_back(" N7 ");
    missing_atom_vec.push_back(" C5 ");
    missing_atom_vec.push_back(" C6 ");
    missing_atom_vec.push_back(" O6 ");
    missing_atom_vec.push_back(" N1 ");
    missing_atom_vec.push_back(" C2 ");
    missing_atom_vec.push_back(" N2 ");
    missing_atom_vec.push_back(" N3 ");
    missing_atom_vec.push_back(" C4 ");
    missing_rna["  G"]=missing_atom_vec;
    for (size_t a=0;a<missing_atom_vec.size();a++) missing_atom_vec[a].clear();
    missing_atom_vec.clear();

    //missing_atom_vec.push_back(" P  ");
    missing_atom_vec.push_back(" OP1");
    missing_atom_vec.push_back(" OP2");
    missing_atom_vec.push_back(" O5'");
    missing_atom_vec.push_back(" C5'");
    //missing_atom_vec.push_back(" C4'");
    missing_atom_vec.push_back(" O4'");
    missing_atom_vec.push_back(" C3'");
    missing_atom_vec.push_back(" O3'");
    missing_atom_vec.push_back(" C2'");
    missing_atom_vec.push_back(" O2'");
    //missing_atom_vec.push_back(" C1'");
    missing_atom_vec.push_back(" N1 ");
    missing_atom_vec.push_back(" C2 ");
    missing_atom_vec.push_back(" O2 ");
    missing_atom_vec.push_back(" N3 ");
    missing_atom_vec.push_back(" C4 ");
    missing_atom_vec.push_back(" O4 ");
    missing_atom_vec.push_back(" C5 ");
    missing_atom_vec.push_back(" C6 ");
    missing_rna["  U"]=missing_atom_vec;
    vector<string>().swap(missing_atom_vec);


    vector<float> tmp(3,0);
    
    map<string,vector<float> >A;
    A[" P  "]=tmp; A[" P  "][0]=   3.063; A[" P  "][1]=   8.025; A[" P  "][2]=  -4.135;
    A[" OP1"]=tmp; A[" OP1"][0]=   3.223; A[" OP1"][1]=   8.856; A[" OP1"][2]=  -5.350;
    A[" OP2"]=tmp; A[" OP2"][0]=   1.891; A[" OP2"][1]=   7.121; A[" OP2"][2]=  -4.118;
    A[" O5'"]=tmp; A[" O5'"][0]=   4.396; A[" O5'"][1]=   7.181; A[" O5'"][2]=  -3.881;
    A[" C5'"]=tmp; A[" C5'"][0]=   5.621; A[" C5'"][1]=   7.881; A[" C5'"][2]=  -3.587;
    A[" C4'"]=tmp; A[" C4'"][0]=   6.719; A[" C4'"][1]=   6.889; A[" C4'"][2]=  -3.258;
    A[" O4'"]=tmp; A[" O4'"][0]=   6.486; A[" O4'"][1]=   6.364; A[" O4'"][2]=  -1.919;
    A[" C3'"]=tmp; A[" C3'"][0]=   6.776; A[" C3'"][1]=   5.642; A[" C3'"][2]=  -4.140;
    A[" O3'"]=tmp; A[" O3'"][0]=   7.397; A[" O3'"][1]=   5.908; A[" O3'"][2]=  -5.390;
    A[" C2'"]=tmp; A[" C2'"][0]=   7.567; A[" C2'"][1]=   4.725; A[" C2'"][2]=  -3.208;
    A[" O2'"]=tmp; A[" O2'"][0]=   8.962; A[" O2'"][1]=   4.911; A[" O2'"][2]=  -3.068;
    A[" C1'"]=tmp; A[" C1'"][0]=   6.874; A[" C1'"][1]=   4.999; A[" C1'"][2]=  -1.877;
    A[" N9 "]=tmp; A[" N9 "][0]=   5.666; A[" N9 "][1]=   4.158; A[" N9 "][2]=  -1.647;
    A[" C8 "]=tmp; A[" C8 "][0]=   4.345; A[" C8 "][1]=   4.449; A[" C8 "][2]=  -1.900;
    A[" N7 "]=tmp; A[" N7 "][0]=   3.524; A[" N7 "][1]=   3.496; A[" N7 "][2]=  -1.584;
    A[" C5 "]=tmp; A[" C5 "][0]=   4.348; A[" C5 "][1]=   2.496; A[" C5 "][2]=  -1.086;
    A[" C6 "]=tmp; A[" C6 "][0]=   4.082; A[" C6 "][1]=   1.215; A[" C6 "][2]=  -0.578;
    A[" N6 "]=tmp; A[" N6 "][0]=   2.849; A[" N6 "][1]=   0.697; A[" N6 "][2]=  -0.485;
    A[" N1 "]=tmp; A[" N1 "][0]=   5.134; A[" N1 "][1]=   0.481; A[" N1 "][2]=  -0.167;
    A[" C2 "]=tmp; A[" C2 "][0]=   6.356; A[" C2 "][1]=   1.001; A[" C2 "][2]=  -0.262;
    A[" N3 "]=tmp; A[" N3 "][0]=   6.725; A[" N3 "][1]=   2.179; A[" N3 "][2]=  -0.716;
    A[" C4 "]=tmp; A[" C4 "][0]=   5.655; A[" C4 "][1]=   2.893; A[" C4 "][2]=  -1.121;
    simple_rna["  A"]=A;
    map<string,vector<float> >().swap(A);

    map<string,vector<float> >C;
    C[" P  "]=tmp; C[" P  "][0]=   6.913; C[" P  "][1]=   5.099; C[" P  "][2]=  -6.683;
    C[" OP1"]=tmp; C[" OP1"][0]=   7.496; C[" OP1"][1]=   5.711; C[" OP1"][2]=  -7.898;
    C[" OP2"]=tmp; C[" OP2"][0]=   5.439; C[" OP2"][1]=   4.971; C[" OP2"][2]=  -6.666;
    C[" O5'"]=tmp; C[" O5'"][0]=   7.579; C[" O5'"][1]=   3.667; C[" O5'"][2]=  -6.429;
    C[" C5'"]=tmp; C[" C5'"][0]=   8.987; C[" C5'"][1]=   3.596; C[" C5'"][2]=  -6.135;
    C[" C4'"]=tmp; C[" C4'"][0]=   9.376; C[" C4'"][1]=   2.167; C[" C4'"][2]=  -5.806;
    C[" O4'"]=tmp; C[" O4'"][0]=   8.896; C[" O4'"][1]=   1.851; C[" O4'"][2]=  -4.467;
    C[" C3'"]=tmp; C[" C3'"][0]=   8.750; C[" C3'"][1]=   1.087; C[" C3'"][2]=  -6.688;
    C[" O3'"]=tmp; C[" O3'"][0]=   9.418; C[" O3'"][1]=   0.975; C[" O3'"][2]=  -7.938;
    C[" C2'"]=tmp; C[" C2'"][0]=   8.920; C[" C2'"][1]=  -0.112; C[" C2'"][2]=  -5.756;
    C[" O2'"]=tmp; C[" O2'"][0]=  10.194; C[" O2'"][1]=  -0.710; C[" O2'"][2]=  -5.617;
    C[" C1'"]=tmp; C[" C1'"][0]=   8.486; C[" C1'"][1]=   0.494; C[" C1'"][2]=  -4.425;
    C[" N1 "]=tmp; C[" N1 "][0]=   7.014; C[" N1 "][1]=   0.438; C[" N1 "][2]=  -4.195;
    C[" C2 "]=tmp; C[" C2 "][0]=   6.487; C[" C2 "][1]=  -0.729; C[" C2 "][2]=  -3.648;
    C[" O2 "]=tmp; C[" O2 "][0]=   7.253; C[" O2 "][1]=  -1.661; C[" O2 "][2]=  -3.379;
    C[" N3 "]=tmp; C[" N3 "][0]=   5.148; C[" N3 "][1]=  -0.800; C[" N3 "][2]=  -3.431;
    C[" C4 "]=tmp; C[" C4 "][0]=   4.352; C[" C4 "][1]=   0.233; C[" C4 "][2]=  -3.736;
    C[" N4 "]=tmp; C[" N4 "][0]=   3.054; C[" N4 "][1]=   0.114; C[" N4 "][2]=  -3.504;
    C[" C5 "]=tmp; C[" C5 "][0]=   4.874; C[" C5 "][1]=   1.442; C[" C5 "][2]=  -4.299;
    C[" C6 "]=tmp; C[" C6 "][0]=   6.214; C[" C6 "][1]=   1.492; C[" C6 "][2]=  -4.509;
    simple_rna["  C"]=C;
    map<string,vector<float> >().swap(C);

    map<string,vector<float> >G;
    G[" P  "]=tmp; G[" P  "][0]=   8.572; G[" P  "][1]=   0.556; G[" P  "][2]=  -9.231;
    G[" OP1"]=tmp; G[" OP1"][0]=   9.394; G[" OP1"][1]=   0.756; G[" OP1"][2]= -10.446;
    G[" OP2"]=tmp; G[" OP2"][0]=   7.262; G[" OP2"][1]=   1.245; G[" OP2"][2]=  -9.214;
    G[" O5'"]=tmp; G[" O5'"][0]=   8.359; G[" O5'"][1]=  -1.009; G[" O5'"][2]=  -8.977;
    G[" C5'"]=tmp; G[" C5'"][0]=   9.505; G[" C5'"][1]=  -1.830; G[" C5'"][2]=  -8.683;
    G[" C4'"]=tmp; G[" C4'"][0]=   9.061; G[" C4'"][1]=  -3.241; G[" C4'"][2]=  -8.354;
    G[" O4'"]=tmp; G[" O4'"][0]=   8.485; G[" O4'"][1]=  -3.248; G[" O4'"][2]=  -7.015;
    G[" C3'"]=tmp; G[" C3'"][0]=   7.951; G[" C3'"][1]=  -3.812; G[" C3'"][2]=  -9.236;
    G[" O3'"]=tmp; G[" O3'"][0]=   8.452; G[" O3'"][1]=  -4.267; G[" O3'"][2]= -10.486;
    G[" C2'"]=tmp; G[" C2'"][0]=   7.447; G[" C2'"][1]=  -4.913; G[" C2'"][2]=  -8.304;
    G[" O2'"]=tmp; G[" O2'"][0]=   8.195; G[" O2'"][1]=  -6.103; G[" O2'"][2]=  -8.164;
    G[" C1'"]=tmp; G[" C1'"][0]=   7.407; G[" C1'"][1]=  -4.169; G[" C1'"][2]=  -6.973;
    G[" N9 "]=tmp; G[" N9 "][0]=   6.139; G[" N9 "][1]=  -3.421; G[" N9 "][2]=  -6.743;
    G[" C8 "]=tmp; G[" C8 "][0]=   5.851; G[" C8 "][1]=  -2.097; G[" C8 "][2]=  -6.995;
    G[" N7 "]=tmp; G[" N7 "][0]=   4.628; G[" N7 "][1]=  -1.748; G[" N7 "][2]=  -6.675;
    G[" C5 "]=tmp; G[" C5 "][0]=   4.069; G[" C5 "][1]=  -2.923; G[" C5 "][2]=  -6.175;
    G[" C6 "]=tmp; G[" C6 "][0]=   2.765; G[" C6 "][1]=  -3.169; G[" C6 "][2]=  -5.670;
    G[" O6 "]=tmp; G[" O6 "][0]=   1.822; G[" O6 "][1]=  -2.391; G[" O6 "][2]=  -5.558;
    G[" N1 "]=tmp; G[" N1 "][0]=   2.618; G[" N1 "][1]=  -4.505; G[" N1 "][2]=  -5.268;
    G[" C2 "]=tmp; G[" C2 "][0]=   3.600; G[" C2 "][1]=  -5.472; G[" C2 "][2]=  -5.344;
    G[" N2 "]=tmp; G[" N2 "][0]=   3.260; G[" N2 "][1]=  -6.688; G[" N2 "][2]=  -4.908;
    G[" N3 "]=tmp; G[" N3 "][0]=   4.821; G[" N3 "][1]=  -5.239; G[" N3 "][2]=  -5.818;
    G[" C4 "]=tmp; G[" C4 "][0]=   4.982; G[" C4 "][1]=  -3.949; G[" C4 "][2]=  -6.213;
    simple_rna["  G"]=G;
    map<string,vector<float> >().swap(G);

    map<string,vector<float> >U;
    U[" P  "]=tmp; U[" P  "][0]=   7.514; U[" P  "][1]=  -4.163; U[" P  "][2]= -11.779;
    U[" OP1"]=tmp; U[" OP1"][0]=   8.313; U[" OP1"][1]=  -4.439; U[" OP1"][2]= -12.994;
    U[" OP2"]=tmp; U[" OP2"][0]=   6.784; U[" OP2"][1]=  -2.876; U[" OP2"][2]= -11.762;
    U[" O5'"]=tmp; U[" O5'"][0]=   6.490; U[" O5'"][1]=  -5.365; U[" O5'"][2]= -11.525;
    U[" C5'"]=tmp; U[" C5'"][0]=   7.010; U[" C5'"][1]=  -6.675; U[" C5'"][2]= -11.231;
    U[" C4'"]=tmp; U[" C4'"][0]=   5.874; U[" C4'"][1]=  -7.623; U[" C4'"][2]= -10.902;
    U[" O4'"]=tmp; U[" O4'"][0]=   5.387; U[" O4'"][1]=  -7.318; U[" O4'"][2]=  -9.563;
    U[" C3'"]=tmp; U[" C3'"][0]=   4.631; U[" C3'"][1]=  -7.503; U[" C3'"][2]= -11.784;
    U[" O3'"]=tmp; U[" O3'"][0]=   4.807; U[" O3'"][1]=  -8.156; U[" O3'"][2]= -13.034;
    U[" C2'"]=tmp; U[" C2'"][0]=   3.612; U[" C2'"][1]=  -8.157; U[" C2'"][2]= -10.852;
    U[" O2'"]=tmp; U[" O2'"][0]=   3.598; U[" O2'"][1]=  -9.565; U[" O2'"][2]= -10.712;
    U[" C1'"]=tmp; U[" C1'"][0]=   3.981; U[" C1'"][1]=  -7.510; U[" C1'"][2]=  -9.521;
    U[" N1 "]=tmp; U[" N1 "][0]=   3.318; U[" N1 "][1]=  -6.195; U[" N1 "][2]=  -9.291;
    U[" C2 "]=tmp; U[" C2 "][0]=   2.055; U[" C2 "][1]=  -6.218; U[" C2 "][2]=  -8.751;
    U[" O2 "]=tmp; U[" O2 "][0]=   1.475; U[" O2 "][1]=  -7.251; U[" O2 "][2]=  -8.463;
    U[" N3 "]=tmp; U[" N3 "][0]=   1.473; U[" N3 "][1]=  -4.982; U[" N3 "][2]=  -8.552;
    U[" C4 "]=tmp; U[" C4 "][0]=   2.034; U[" C4 "][1]=  -3.755; U[" C4 "][2]=  -8.841;
    U[" O4 "]=tmp; U[" O4 "][0]=   1.415; U[" O4 "][1]=  -2.713; U[" O4 "][2]=  -8.618;
    U[" C5 "]=tmp; U[" C5 "][0]=   3.362; U[" C5 "][1]=  -3.834; U[" C5 "][2]=  -9.404;
    U[" C6 "]=tmp; U[" C6 "][0]=   3.953; U[" C6 "][1]=  -5.023; U[" C6 "][2]=  -9.608;
    simple_rna["  U"]=U;
    map<string,vector<float> >().swap(U);

    vector<float> ().swap(tmp);
    return;
}

bool fillSimpleMissingRNAatom(ResidueUnit &residue,
    map<string, map<string,vector<float> > >&simple_rna,
    map<string,vector<string> >&missing_rna)
{
    if (missing_rna.count(residue.resn)==0)
    {
        cerr<<"cannot fill known residue "<<residue.resn<<endl;
        return false;
    }

    vector<float> tmp(3,0);
    vector<vector<float> > xyz_list1(residue.atoms.size(),tmp);
    vector<vector<float> > xyz_list2(residue.atoms.size(),tmp);
    vector<vector<float> > RotMatix;  // U
    vector<float> TranVect;  // t

    size_t a;
    for (a=0;a<residue.atoms.size();a++)
    {
        xyz_list1[a][0]=simple_rna[residue.resn][residue.atoms[a].name][0];
        xyz_list1[a][1]=simple_rna[residue.resn][residue.atoms[a].name][1];
        xyz_list1[a][2]=simple_rna[residue.resn][residue.atoms[a].name][2];

        xyz_list2[a][0]=residue.atoms[a].xyz[0];
        xyz_list2[a][1]=residue.atoms[a].xyz[1];
        xyz_list2[a][2]=residue.atoms[a].xyz[2];
    }

    RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);
    
    for (a=0;a<residue.atoms.size();a++)
    {
        xyz_list1[a].clear();
        xyz_list2[a].clear();
    }
    xyz_list1.clear();
    xyz_list2.clear();
                
    AtomUnit atom;
    atom.xyz.assign(3,0);
    for (a=0;a<missing_rna[residue.resn].size();a++)
    {
        atom.name=missing_rna[residue.resn][a];
        ChangeCoor(simple_rna[residue.resn][atom.name],
            RotMatix, TranVect, atom.xyz);
        residue.atoms.push_back(atom);
    }

    atom.name.clear();
    atom.xyz.clear();
    TranVect.clear();
    for (a=0;a<3;a++) RotMatix[a].clear();
    RotMatix.clear();
    return true;
}


void calculateCoordinates(vector<float>&D, 
    const vector<float>&refA, const vector<float>&refB, const vector<float>&refC, 
    const float l, const float ang, const float di)
{
    vector<float>CA(3,0);
    vector<float>CB(3,0);
    subtract(refA,refB,CA);
    subtract(refB,refC,CB);

    /* Plane Parameters */
    float A=(CA[1]*CB[2]) - (CA[2]*CB[1]);
    float B=(CA[2]*CB[0]) - (CA[0]*CB[2]);
    float G=(CA[0]*CB[1]) - (CA[1]*CB[0]);

    /* Dot Product Constant */
    float F=sqrt(CB[0]*CB[0] + CB[1]*CB[1] + CB[2]*CB[2]
        ) * l * cos(deg2rad(ang));

    float cons=B*CB[2]-CB[1]*G;
    cons=sqrt(cons*cons*(-(F*F)*(A*A+B*B+G*G)+(B*B*(CB[0]*CB[0]+CB[2]*CB[2]
        ) + A*A*(CB[1]*CB[1]+CB[2]*CB[2])- (2*A*CB[0]*CB[2]*G) + (
        CB[0]*CB[0]+ CB[1]*CB[1])*G*G- (2*B*CB[1])*(A*CB[0]+CB[2]*G))*l*l));
    float denom=B*B*(CB[0]*CB[0]+CB[2]*CB[2]) + A*A*(CB[1]*CB[1]+CB[2]*CB[2]
        ) - (2*A*CB[0]*CB[2]*G) + G*G*(CB[0]*CB[0]+CB[1]*CB[1]
        ) - (2*B*CB[1])*(A*CB[0]+CB[2]*G);

    float X=((B*B*CB[0]*F)-(A*B*CB[1]*F)+(F*G)*(-A*CB[2]+CB[0]*G)+cons
        )/denom;

    float Y,Z;
    if ((B==0 or CB[2]==0) && (CB[1]==0 or G==0))
    {
        float const1=sqrt( G*G*(-A*A*X*X+(B*B+G*G)*(l-X)*(l+X)));
        Y= ((-A*B*X)+const1)/(B*B+G*G);
        Z= -(A*G*G*X+B*const1)/(G*(B*B+G*G));
    }
    else
    {
        Y= ((A*A*CB[1]*F)*(B*CB[2]-CB[1]*G)+ G*(-F*pow(B*CB[2]-CB[1]*G,2
            ) + CB[0]*cons) - A*( B*B*CB[0]*CB[2]*F- B*CB[0]*CB[1]*F*G + \
            CB[2]*cons)) / ((B*CB[2]-CB[1]*G)*denom);
        Z= ((A*A*CB[2]*F)*(B*CB[2]-CB[1]*G) + (B*F)*pow(B*CB[2]-CB[1]*G,2
            ) + (A*CB[0]*F*G)*(-B*CB[2]+CB[1]*G) - B*CB[0]*cons + \
            A*CB[1]*cons) / ((B*CB[2]-CB[1]*G)*denom);
    }

    
    /* GET THE NEW VECTOR from the orgin */
    vector<float>tmp_array(3,0);
    tmp_array[0]=X;
    tmp_array[1]=Y;
    tmp_array[2]=Z;
    vectorsum(tmp_array, refC, D);

    float angle=di-rad2deg(Points2Dihedral(refA, refB, refC, D));
    CoordinateRotation(D,refC,refB, angle, tmp_array);
    for (int i=0;i<3;i++) D[i]=tmp_array[i];

    /* clean up */
    tmp_array.clear();
    CA.clear();
    CB.clear();
}

void make_chain(vector<vector<float> >&P_coord, vector<vector<float> >&C4_coord,
    vector<vector<float> >&C1_coord, const vector<vector<float> >&rama_table)
{
    int L=C4_coord.size();
    if (L==0) return;

    /* define statistical values for backbone geometry of nucleotide */
    /*
    float PC4   =  3.849; // P-C4'
    float C4Pp  =  3.785; // C4'-P[+1]
    float C4C1  =  2.341; // C4'-C1'

    float C4mPC4= 107.14; // C4'[-1]-P-C4'
    float PC4Pp = 105.29; // P-C4'-P[+1]
    float PC4C1 = 115.69; // P-C4'-C1'

    float i1    =-112.85; // P[+1]-P-C4'-C1'
    */
    
    /* define ideal helical values for backbone geometry of nucleotide */
    float PC4   =  3.927; // P-C4'
    float C4Pp  =  3.870; // C4'-P[+1]
    float C4C1  =  2.346; // C4'-C1'

    float C4mPC4= 100.32; // C4'[-1]-P-C4'
    float PC4Pp =  89.02; // P-C4'-P[+1]
    float PC4C1 = 115.23; // P-C4'-C1'

    float i1    = -98.77; // P[+1]-P-C4'-C1'

    /* initialize first residue */
    P_coord[1][0]=C4Pp;
    P_coord[0][0]=PC4 * cos(deg2rad(PC4Pp));
    P_coord[0][1]=PC4 * sin(deg2rad(PC4Pp));

    calculateCoordinates(C1_coord[0],
        P_coord[1], P_coord[0], C4_coord[0], C4C1, PC4C1, i1);

    /* add subsequent residues */
    int r;
    for (r=1;r<L;r++)
    {
        calculateCoordinates(C4_coord[r],
            P_coord[r-1], C4_coord[r-1], P_coord[r], PC4, C4mPC4, rama_table[r-1][1]);
        calculateCoordinates(P_coord[r+1],
            C4_coord[r-1], P_coord[r], C4_coord[r], C4Pp, PC4Pp, rama_table[r][0]);
        calculateCoordinates(C1_coord[r],
            P_coord[r+1], P_coord[r], C4_coord[r], C4C1, PC4C1, i1);
    }

    /* center coordinate */
    vector<float>mean_coor(3,0);
    for (r=0;r<L;r++)
    {
        mean_coor[0]+=P_coord[r][0];
        mean_coor[1]+=P_coord[r][1];
        mean_coor[2]+=P_coord[r][2];
    }
    mean_coor[0]/=L;
    mean_coor[1]/=L;
    mean_coor[2]/=L;
    for (r=0;r<L;r++)
    {
        P_coord[r][0] -=mean_coor[0];
        P_coord[r][1] -=mean_coor[1];
        P_coord[r][2] -=mean_coor[2];

        C4_coord[r][0]-=mean_coor[0];
        C4_coord[r][1]-=mean_coor[1];
        C4_coord[r][2]-=mean_coor[2];

        C1_coord[r][0]-=mean_coor[0];
        C1_coord[r][1]-=mean_coor[1];
        C1_coord[r][2]-=mean_coor[2];
    }
    P_coord[L][0]-=mean_coor[0];
    P_coord[L][1]-=mean_coor[1];
    P_coord[L][2]-=mean_coor[2];
}

void randomRNAseq(string& sequence, const int L)
{
    int r;
    float p;
    for (r=0;r<L;r++)
    {
        p=1.*rand()/RAND_MAX;
        if      (p<0.25879407) sequence+='a'; // 0.25879407
        else if (p<0.47917813) sequence+='u'; // 0.22038406
        else if (p<0.70690686) sequence+='c'; // 0.22772873
        else                   sequence+='g'; // 0.29309314
    }
    return;
}

void convert_chain(ChainUnit & chain, const string & sequence,
    const vector<vector<float> >&P_coord, const vector<vector<float> >&C4_coord,
    const vector<vector<float> >&C1_coord)
{
    chain.chainID='A';
    int L=C4_coord.size();

    ResidueUnit residue;
    AtomUnit atom;
    atom.xyz.assign(3,0);
    residue.resn="  ";
    residue.atoms.assign(3,atom);
    atom.xyz.clear();
    chain.residues.assign(L,residue);
    residue.atoms.clear();
    residue.resn.clear();

    int r;
    for (r=0;r<L;r++)
    {
        chain.residues[r].resi =r+1;
        chain.residues[r].het  =false;
        chain.residues[r].icode=' ';
        chain.residues[r].resn+=toupper(sequence[r]);

        chain.residues[r].atoms[0].name=" P  ";
        chain.residues[r].atoms[0].xyz[0]=P_coord[r][0];
        chain.residues[r].atoms[0].xyz[1]=P_coord[r][1];
        chain.residues[r].atoms[0].xyz[2]=P_coord[r][2];

        chain.residues[r].atoms[1].name=" C4'";
        chain.residues[r].atoms[1].xyz[0]=C4_coord[r][0];
        chain.residues[r].atoms[1].xyz[1]=C4_coord[r][1];
        chain.residues[r].atoms[1].xyz[2]=C4_coord[r][2];

        chain.residues[r].atoms[2].name=" C1'";
        chain.residues[r].atoms[2].xyz[0]=C1_coord[r][0];
        chain.residues[r].atoms[2].xyz[1]=C1_coord[r][1];
        chain.residues[r].atoms[2].xyz[2]=C1_coord[r][2];
    }
}

void randomRamaTable(vector<vector<float> >&rama_table)
{
    int L=rama_table.size();
    int cdf_sum=pseudo_torsion_cdf[360-1][360-1]; // 248502
    int rnd;
    int i,j,r;
    float eta,theta;
    for (r=0;r<L;r++)
    {
        rnd=(rand() % cdf_sum); // max rand is INT_MAX=2147483647
        for (i=0;i<360;i++) // theta
        {
            if (rnd>=pseudo_torsion_cdf[i][360-1]) continue;
            theta=i-180+0.5;
            for (j=0;j<360;j++) // eta
            {
                if (rnd>=pseudo_torsion_cdf[i][j]) continue;
                eta=j-180+0.5;
                rama_table[r][0]=eta;
                rama_table[r][1]=theta;
                break;
            }
            break;
        }
    }
}

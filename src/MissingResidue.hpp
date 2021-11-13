#ifndef MissingResidue_HPP
#define MissingResidue_HPP 1

#include <iostream>
#include "PDBParser.hpp"
#include "GeometryTools.hpp"

#define RUNIF (float(rand())/RAND_MAX)

string AA1to3(const char A)
{
    if (A=='A') return "ALA";
    if (A=='B') return "ASX";
    if (A=='C') return "CYS";
    if (A=='D') return "ASP";
    if (A=='E') return "GLU";
    if (A=='F') return "PHE";
    if (A=='G') return "GLY";
    if (A=='H') return "HIS";
    if (A=='I') return "ILE";
    if (A=='K') return "LYS";
    if (A=='L') return "LEU";
    if (A=='M') return "MET";
    if (A=='N') return "ASN";
    if (A=='O') return "PYL";
    if (A=='P') return "PRO";
    if (A=='Q') return "GLN";
    if (A=='R') return "ARG";
    if (A=='S') return "SER";
    if (A=='T') return "THR";
    if (A=='U') return "SEC";
    if (A=='V') return "VAL";
    if (A=='W') return "TRP";    
    if (A=='Y') return "TYR";
    if (A=='Z') return "GLX";
    if ('a'<=A && A<='z') return "  "+string(1,char(toupper(A)));
    return "UNK";
}

void readOneSequenceFromFasta(const string & infasta, string &sequence,
    const int mol_opt=0)
{
    ifstream fp;
    bool use_stdin=(infasta=="-");
    if (!use_stdin) fp.open(infasta.c_str(),ios::in);
    string line;
    while(use_stdin?cin.good():fp.good())
    {
        if (use_stdin) getline(cin,line);
        else           getline(fp,line);
        if (line.size()==0 || line[0]=='>' || line[0]=='#') continue;
        sequence+=line;
    }
    line.clear();
    if (mol_opt>0) transform(sequence.begin(), sequence.end(), sequence.begin(),
        [](unsigned char c){ return std::tolower(c); }); // RNA
    else if (mol_opt<0) transform(sequence.begin(), sequence.end(), sequence.begin(),
        [](unsigned char c){ return std::toupper(c); }); // protein
    return;
}

/* >0 for RNA, <0 for protein */
int check_mol(const ChainUnit & chain)
{
    int r;
    int mol_opt=0;
    for (r=0;r<chain.residues.size();r++)
    {
        if (chain.residues[r].resn[0]==' ') mol_opt++;
        else mol_opt--;
    }
    return mol_opt;
}

void existCoor2vec(const string & sequence, const ChainUnit & chain,
    vector<vector<float> >&xyz_full, vector<bool>&miss_list, const int mol_opt)
{
    int L=sequence.size();
    if (mol_opt>0) L*=2;
    miss_list.assign(L,true);
    vector<float>xyz_tmp(3,0);
    xyz_full.assign(L,xyz_tmp);
    int r,resi,a;

    for (r=0;r<chain.residues.size();r++)
    {
        resi=chain.residues[r].resi;
        if (resi<=0 || resi>sequence.size())
        {
            cerr<<"ERROR! residue index "<<resi<<endl;
            continue;
        }
        for (a=0;a<chain.residues[r].atoms.size();a++)
        {
            if (mol_opt<0 && chain.residues[r].atoms[a].name==" CA ")
            {
                miss_list[resi-1]=false;
                xyz_full[resi-1][0]=chain.residues[r].atoms[a].xyz[0];
                xyz_full[resi-1][1]=chain.residues[r].atoms[a].xyz[1];
                xyz_full[resi-1][2]=chain.residues[r].atoms[a].xyz[2];
            }
            else if (mol_opt>0 && chain.residues[r].atoms[a].name==" P  ")
            {
                miss_list[resi*2-2]=false;
                xyz_full[resi*2-2][0]=chain.residues[r].atoms[a].xyz[0];
                xyz_full[resi*2-2][1]=chain.residues[r].atoms[a].xyz[1];
                xyz_full[resi*2-2][2]=chain.residues[r].atoms[a].xyz[2];
            }
            else if (mol_opt>0 && chain.residues[r].atoms[a].name==" C4'")
            {
                miss_list[resi*2-1]=false;
                xyz_full[resi*2-1][0]=chain.residues[r].atoms[a].xyz[0];
                xyz_full[resi*2-1][1]=chain.residues[r].atoms[a].xyz[1];
                xyz_full[resi*2-1][2]=chain.residues[r].atoms[a].xyz[2];
            }
        }
    }

    xyz_tmp.clear();
    return;
}

void appendNewResidue(const string & sequence,ChainUnit & chain,
    const vector<vector<float> > &xyz_full,const vector<bool>&miss_list,
    const int mol_opt)
{
    vector<ResidueUnit> residues_new;
    ResidueUnit res_tmp;
    AtomUnit   atom_tmp;
    res_tmp.het=false;
    res_tmp.icode=' '; 
    atom_tmp.xyz.assign(3,0);
    if (mol_opt<0)
    {
        atom_tmp.name=" CA ";
        res_tmp.atoms.push_back(atom_tmp);
    }
    else if (mol_opt>0)
    {
        atom_tmp.name=" P  ";
        res_tmp.atoms.push_back(atom_tmp);
        atom_tmp.name=" C4'";
        res_tmp.atoms.push_back(atom_tmp);
    }

    map<int,int> resi_map; // key is resi-1, value is position in chain
    int r,resi;
    for (r=0;r<chain.residues.size();r++)
        resi_map[chain.residues[r].resi-1]=r;

    int L=sequence.size();
    for (resi=0;resi<L;resi++)
    {
        if (resi_map.count(resi))
        {
            r=resi_map[resi];
            residues_new.push_back(chain.residues[r]);
        }
        else
        {
            res_tmp.resi=resi+1;
            res_tmp.resn=AA1to3(sequence[resi]);
            if (mol_opt<0) // protein
            {
                res_tmp.atoms[0].xyz[0]=xyz_full[resi][0];
                res_tmp.atoms[0].xyz[1]=xyz_full[resi][1];
                res_tmp.atoms[0].xyz[2]=xyz_full[resi][2];
            }
            else if (mol_opt>0)
            {
                res_tmp.atoms[0].xyz[0]=xyz_full[resi*2][0];
                res_tmp.atoms[0].xyz[1]=xyz_full[resi*2][1];
                res_tmp.atoms[0].xyz[2]=xyz_full[resi*2][2];
                res_tmp.atoms[1].xyz[0]=xyz_full[resi*2+1][0];
                res_tmp.atoms[1].xyz[1]=xyz_full[resi*2+1][1];
                res_tmp.atoms[1].xyz[2]=xyz_full[resi*2+1][2];
            }
            residues_new.push_back(res_tmp);
        }
    }

    chain.residues.clear();
    for (resi=0;resi<L;resi++) chain.residues.push_back(residues_new[resi]);

    string ().swap(res_tmp.resn);
    vector<AtomUnit>().swap(res_tmp.atoms);
    string ().swap(atom_tmp.name);
    vector<float> ().swap(atom_tmp.xyz);
    vector<ResidueUnit>().swap(residues_new);
    return;
}

void RandomSphereSampling(vector<vector<float> >&cart_coor_array,float d=3.8)
{
    float z,theta;
    int i;
    for (i=0;i<cart_coor_array.size();i++)
    {
        z=2*RUNIF-1;
        theta=RUNIF*2*PI;
        cart_coor_array[i][0]=d*sqrt(1-z*z)*cos(theta); // x
        cart_coor_array[i][1]=d*sqrt(1-z*z)*sin(theta); // y
        cart_coor_array[i][2]=d*z;
    }
}

void fillLoop(int rstart, int rend, vector<vector<float> >&xyz_full,
    vector<bool>&miss_list, float d=3.8)
{
    float l=0; // distance between anchors
    int M=rend-rstart+1; // number of missing atoms
    int toggle=-1; // toggle<0: fill left atom; toggle>0: fill right atom

    vector<float> anchor_vec(3,0);
    vector<float> rotate_axis(3,0);
    vector<float> xyz(3,0);
    vector<float> xyz_old(3,0);
    float theta, theta_max;
    float phi;
    float z;
    float angle;

    while (rstart<=rend)
    {
        M=rend-rstart+1; // number of missing atoms
        l=Points2Distance(xyz_full[rstart-1],xyz_full[rend+1]);
        if (M==1 && M*d<=l-d) 
        {
            xyz_full[rstart][0]=(xyz_full[rstart-1][0]+xyz_full[rend+1][0])/2.;
            xyz_full[rstart][1]=(xyz_full[rstart-1][1]+xyz_full[rend+1][1])/2.;
            xyz_full[rstart][2]=(xyz_full[rstart-1][2]+xyz_full[rend+1][2])/2.;
            miss_list[rstart]=false;
            rstart++;
            rend--;
            break;
        }

        /* biased sphere sampling */
        if (M*d<=l-d) theta=0;
        else
        {
            if (M*d>=l+d) theta_max=PI;
            else theta_max=acos((l*l+d*d+M*d*M*d)/(2*l*d));
            theta=theta_max*RUNIF;
        }
        phi=(2*RUNIF-1)*PI;
        if (theta<PI/2)
        {
            xyz[0]=d*cos(phi)*sin(theta);
            xyz[1]=d*sin(phi)*sin(theta);
            xyz[2]=d*cos(theta);
        }
        else
        {
            z=-RUNIF;
            xyz[0]=d*sqrt(1-z*z)*cos(phi);
            xyz[1]=d*sqrt(1-z*z)*sin(phi);
            xyz[2]=d*z;
        }

        /* rotate xyz according to anchor_vec */
        if (l>Extra)
        {
            if (toggle<0) subtract(xyz_full[rend+1],xyz_full[rstart-1],anchor_vec);
            else          subtract(xyz_full[rstart-1],xyz_full[rend+1],anchor_vec);
            if (l>Extra)
            {
                angle=acos(anchor_vec[2]/l);
                /* rotation axis is cross product between anchor_vec and z */
                if (angle>=Extra)
                {
                    rotate_axis[0]= anchor_vec[1];
                    rotate_axis[1]=-anchor_vec[0];
                    rotate_axis[2]=0;
                    xyz_old[0]=xyz[0];
                    xyz_old[1]=xyz[1];
                    xyz_old[2]=xyz[2];
                    CoordinateRotation(xyz_old, rotate_axis, angle, xyz);
                }
            }
        }
        if (toggle<0)
        {
            xyz_full[rstart][0]=xyz_full[rstart-1][0]+xyz[0];
            xyz_full[rstart][1]=xyz_full[rstart-1][1]+xyz[1];
            xyz_full[rstart][2]=xyz_full[rstart-1][2]+xyz[2];
            miss_list[rstart]=false;
            rstart++;
        }
        else
        {
            xyz_full[rend][0]  =xyz_full[rend+1][0]  +xyz[0];
            xyz_full[rend][1]  =xyz_full[rend+1][1]  +xyz[1];
            xyz_full[rend][2]  =xyz_full[rend+1][2]  +xyz[2];
            miss_list[rend]=false;
            rend--;
        }
        toggle=-toggle;
    }
    /* clean up */
    vector<float>().swap(anchor_vec);
    vector<float>().swap(rotate_axis);
    vector<float>().swap(xyz);
    vector<float>().swap(xyz_old);
    rstart=-1;
    rend=-1;
    return;
}

/* main function for atom filling */
void MissingResidue(vector<vector<float> >&xyz_full, vector<bool>miss_list)
{
    int L=miss_list.size();
    int i,j;
    float d=3.8; // distance between adjacent atoms
    
    /* fill C terminl tail */
    int Ltail;
    vector<vector<float> > cart_coor_array;
    Ltail=0;
    for (i=L-1;i>=0;i--)
    {
        if (miss_list[i]==false) break;
        Ltail++;
    }
    if (Ltail)
    {
        cart_coor_array.assign(Ltail,xyz_full[L-1]);
        RandomSphereSampling(cart_coor_array,d);
        j=0;
        for (i=L-Ltail;i<L;i++)
        {
            xyz_full[i][0]=xyz_full[i-1][0]+cart_coor_array[j][0];
            xyz_full[i][1]=xyz_full[i-1][1]+cart_coor_array[j][1];
            xyz_full[i][2]=xyz_full[i-1][2]+cart_coor_array[j][2];
            cart_coor_array[j].clear();
            j++;
        }
        cart_coor_array.clear();
    }
   
    /* fill N terminl tail */
    Ltail=0;
    for (i=0;i<L;i++)
    {
        if (miss_list[i]==false) break;
        Ltail++;
    }
    if (Ltail)
    {
        cart_coor_array.assign(Ltail,xyz_full[0]);
        RandomSphereSampling(cart_coor_array,d);
        j=0;
        for (i=Ltail-1;i>=0;i--)
        {
            xyz_full[i][0]=xyz_full[i+1][0]+cart_coor_array[j][0];
            xyz_full[i][1]=xyz_full[i+1][1]+cart_coor_array[j][1];
            xyz_full[i][2]=xyz_full[i+1][2]+cart_coor_array[j][2];
            cart_coor_array[j].clear();
            j++;
        }
        cart_coor_array.clear();
    }

    /* fill loop in the middle */
    int rstart=-1;
    int rend=-1;
    for (i=Ltail;i<L;i++)
    {
        if (rstart==-1)
        {
            if (miss_list[i]==false) continue;
            rstart=i;
            rend=i;
            continue;
        }
        if (miss_list[i]==true)
        {
            rend=i;
            continue;
        }
        fillLoop(rstart,rend,xyz_full,miss_list,d);
    }
   
    /* clean up */
    vector<vector<float> >().swap(cart_coor_array);
    return;
}
#endif

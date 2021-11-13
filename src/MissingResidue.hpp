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
    string sequence_tmp="";
    int r;
    while(use_stdin?cin.good():fp.good())
    {
        if (use_stdin) getline(cin,line);
        else           getline(fp,line);
        if (line.size()==0 || line[0]=='>' || line[0]=='#') continue;
        sequence_tmp+=line;
    }
    line.clear();
    if (mol_opt>0) for (r=0;r<sequence_tmp.size();r++) 
        sequence+=tolower(sequence_tmp[r]);
    else           for (r=0;r<sequence_tmp.size();r++) 
        sequence+=toupper(sequence_tmp[r]);
    sequence_tmp.clear();
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

void RandomSphereSampling(vector<float> &xyz,float d=3.8)
{
    float z,theta;
    z=2*RUNIF-1;
    theta=RUNIF*2*PI;
    xyz[0]=d*sqrt(1-z*z)*cos(theta); // x
    xyz[1]=d*sqrt(1-z*z)*sin(theta); // y
    xyz[2]=d*z;
}

/* clash_cut=3.4 - van der Waals diameter of carbon + 0.1 */
float getClash(const vector<vector<float> >&xyz_full, 
    const vector<vector<int> >&int_coor_array,
    const vector<bool>&miss_list, const int resi)
{
    int L=xyz_full.size();
    float clash_cut2=4.0*4.0; // non-adjacent atoms
    float clash_nei2=3.5*3.5; // adjacent atoms
    float clash_score=0;
    int i;
    float dx,dy,dz,d2;
    vector<int>   resi_xyz_int(3,0);
    vector<float> resi_xyz(3,0);
    resi_xyz_int[0]=int_coor_array[resi][0];
    resi_xyz_int[1]=int_coor_array[resi][1];
    resi_xyz_int[2]=int_coor_array[resi][2];
    resi_xyz[0]    =      xyz_full[resi][0];
    resi_xyz[1]    =      xyz_full[resi][1];
    resi_xyz[2]    =      xyz_full[resi][2];

    for (i=0;i<xyz_full.size();i++)
    {
        if (miss_list[i]) continue;
        if (resi_xyz_int[0]!=int_coor_array[i][0] ||
            resi_xyz_int[1]!=int_coor_array[i][1] ||
            resi_xyz_int[2]!=int_coor_array[i][2]) continue;
        dx=resi_xyz[0]-xyz_full[i][0];
        dy=resi_xyz[1]-xyz_full[i][1];
        dz=resi_xyz[2]-xyz_full[i][2];
        if (resi==i)             d2=0;
        else if (abs(resi-i)==1) d2=clash_nei2-(dx*dx+dy*dy+dz*dz);
        else                     d2=clash_cut2-(dx*dx+dy*dy+dz*dz);
        if (d2>clash_score) clash_score=d2;
    }
    vector<int>   ().swap(resi_xyz_int);
    vector<float> ().swap(resi_xyz);
    return clash_score;
}

void fillLoop(int &rstart, int &rend, vector<vector<float> >&xyz_full,
    vector<vector<int> >&int_coor_array, vector<bool>&miss_list, float d=3.8)
{
    float l=0; // distance between anchors
    int M=rend-rstart+1; // number of missing atoms
    int toggle=-1; // toggle<0: fill left atom; toggle>0: fill right atom

    vector<float> anchor_vec(3,0);
    vector<float> rotate_axis(3,0);
    vector<float> xyz(3,0);
    vector<float> xyz_old(3,0);
    
    vector<float> coor_prev(3,0);
    vector<int> int_coor_prev(3,0);
    float theta, theta_max;
    float phi;
    float z,zmin;
    float angle;
    float Md;
    int resi,resi_anchor;
    int total_trial=5;
    int Ntrial=0;
    float clash_score=0;
    float clash_score_best=0;
    float shrink=1;

    //cerr<<"filling residues "<<rstart<<" to "<<rend<<endl;

    while (rstart<=rend)
    {
        M=rend-rstart+1; // number of missing atoms
        l=Points2Distance(xyz_full[rstart-1],xyz_full[rend+1]);
        Md=M*d;
        if (M==1 && Md<=l-d) 
        {
            xyz_full[rstart][0]=(xyz_full[rstart-1][0]+xyz_full[rend+1][0])/2.;
            xyz_full[rstart][1]=(xyz_full[rstart-1][1]+xyz_full[rend+1][1])/2.;
            xyz_full[rstart][2]=(xyz_full[rstart-1][2]+xyz_full[rend+1][2])/2.;
            int_coor_array[rstart][0]=xyz_full[rstart][0]/4;
            int_coor_array[rstart][1]=xyz_full[rstart][1]/4;
            int_coor_array[rstart][2]=xyz_full[rstart][2]/4;
            miss_list[rstart]=false;
            rstart++;
            rend--;
            break;
        }
        
        shrink=pow(0.95,M-1);
        if (shrink<0.8) shrink=0.8;
        Md*=shrink;
        //Md*=0.9;

        if (toggle<0) subtract(xyz_full[rend+1],xyz_full[rstart-1],anchor_vec);
        else          subtract(xyz_full[rstart-1],xyz_full[rend+1],anchor_vec);
        angle=0;
        if (l>Extra) angle=acos(anchor_vec[2]/l);
        rotate_axis[0]= anchor_vec[1];
        rotate_axis[1]=-anchor_vec[0];
        rotate_axis[2]=0;

        /* uniform cylinder sampling */
        zmin=(d*d+l*l-Md*Md)/(2*l*d);
        if      (zmin> 1) zmin= 1;
        else if (zmin<-1) zmin=-1;
        if (toggle<0)
        {
            resi=rstart;
            resi_anchor=rstart-1;
        }
        else
        {
            resi=rend;
            resi_anchor=rend+1;
        }
        //cerr<<"miss["<<resi<<"]="<<miss_list[resi]<<endl;

        for (Ntrial=0;Ntrial<total_trial;Ntrial++)
        {        
            phi=(2*RUNIF-1)*PI;
            z=RUNIF*(1-zmin)+zmin;
            xyz[0]=d*sqrt(1-z*z)*cos(phi);
            xyz[1]=d*sqrt(1-z*z)*sin(phi);
            xyz[2]=d*z;

            /* rotate xyz according to anchor_vec */
            if (l>Extra && angle>=Extra)
            {
                xyz_old[0]=xyz[0];
                xyz_old[1]=xyz[1];
                xyz_old[2]=xyz[2];
                CoordinateRotation(xyz_old, rotate_axis, rad2deg(angle), xyz);
            }
            xyz_full[resi][0]=xyz_full[resi_anchor][0]+xyz[0];
            xyz_full[resi][1]=xyz_full[resi_anchor][1]+xyz[1];
            xyz_full[resi][2]=xyz_full[resi_anchor][2]+xyz[2];
            int_coor_array[resi][0]=xyz_full[resi][0]/4;
            int_coor_array[resi][1]=xyz_full[resi][1]/4;
            int_coor_array[resi][2]=xyz_full[resi][2]/4;
            if (zmin>=0.99) break;
            clash_score=getClash(xyz_full,int_coor_array,miss_list,resi);
            if (clash_score<=0) break;
            if (Ntrial==0) clash_score_best=clash_score;
            cerr<<"loop resi="<<resi<<";trial="<<Ntrial<<";clash="<<clash_score<<endl;
            if (clash_score<=clash_score_best)
            {
                coor_prev[0]    =xyz_full[resi][0];
                coor_prev[1]    =xyz_full[resi][1];
                coor_prev[2]    =xyz_full[resi][2];
                int_coor_prev[0]=int_coor_array[resi][0];
                int_coor_prev[1]=int_coor_array[resi][1];
                int_coor_prev[2]=int_coor_array[resi][2];
            }
            else if (clash_score>clash_score_best)
            {
                xyz_full[resi][0]=coor_prev[0];
                xyz_full[resi][1]=coor_prev[1];
                xyz_full[resi][2]=coor_prev[2];
                int_coor_array[resi][0]=int_coor_prev[0];
                int_coor_array[resi][1]=int_coor_prev[1];
                int_coor_array[resi][2]=int_coor_prev[2];
            }
        }
        if (toggle<0)
        {
            miss_list[rstart]=false;
            rstart++;
        }
        else
        {
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
    vector<float>().swap(coor_prev);
    vector<int>  ().swap(int_coor_prev);
    rstart=-1;
    rend=-1;
    return;
}

/* main function for atom filling 
 * d - distance between adjacent atoms */
void MissingResidue(vector<vector<float> >&xyz_full, vector<bool>miss_list,
    const float d=3.8)
{
    int L=miss_list.size();
    int i;
    int Ltail;
    vector<float> coor_prev(3,0);
    vector<float> xyz_rnd(3,0);
    vector<int> int_coor_prev(3,0);
    vector<vector<int> > int_coor_array(L,int_coor_prev);
    for (i=0;i<L;i++) // lattice constant =4
    {
        if (miss_list[i]) continue;
        int_coor_array[i][0]=xyz_full[i][0]/4;
        int_coor_array[i][1]=xyz_full[i][1]/4;
        int_coor_array[i][2]=xyz_full[i][2]/4;
    }

    int total_trial=5;
    int Ntrial=0;
    float clash_score=0;
    float clash_score_best=0;
    
    /* fill C terminl tail */
    Ltail=0;
    for (i=L-1;i>=0;i--)
    {
        if (miss_list[i]==false) break;
        Ltail++;
    }
    if (Ltail)
    {
        for (i=L-Ltail;i<L;i++)
        {
            for (Ntrial=0;Ntrial<total_trial;Ntrial++)
            {
                RandomSphereSampling(xyz_rnd,d);
                xyz_full[i][0]      =xyz_full[i-1][0]+xyz_rnd[0];
                xyz_full[i][1]      =xyz_full[i-1][1]+xyz_rnd[1];
                xyz_full[i][2]      =xyz_full[i-1][2]+xyz_rnd[2];
                int_coor_array[i][0]=xyz_full[i][0]/4;
                int_coor_array[i][1]=xyz_full[i][1]/4;
                int_coor_array[i][2]=xyz_full[i][2]/4;
                clash_score=getClash(xyz_full,int_coor_array,miss_list,i);
                if (clash_score<=0) break;
                if (Ntrial==0) clash_score_best=clash_score;
                cerr<<"C resi="<<i<<";trial="<<Ntrial<<";clash="<<clash_score<<endl;
                if (clash_score<=clash_score_best)
                {
                    coor_prev[0]    =xyz_full[i][0];
                    coor_prev[1]    =xyz_full[i][1];
                    coor_prev[2]    =xyz_full[i][2];
                    int_coor_prev[0]=int_coor_array[i][0];
                    int_coor_prev[1]=int_coor_array[i][1];
                    int_coor_prev[2]=int_coor_array[i][2];
                }
                else if (clash_score>clash_score_best)
                {
                    xyz_full[i][0]=coor_prev[0];
                    xyz_full[i][1]=coor_prev[1];
                    xyz_full[i][2]=coor_prev[2];
                    int_coor_array[i][0]=int_coor_prev[0];
                    int_coor_array[i][1]=int_coor_prev[1];
                    int_coor_array[i][2]=int_coor_prev[2];
                }
            }
            miss_list[i]=false;
        }
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
        for (i=Ltail-1;i>=0;i--)
        {
            for (Ntrial=0;Ntrial<total_trial;Ntrial++)
            {
                RandomSphereSampling(xyz_rnd,d);
                xyz_full[i][0]      =xyz_full[i+1][0]+xyz_rnd[0];
                xyz_full[i][1]      =xyz_full[i+1][1]+xyz_rnd[1];
                xyz_full[i][2]      =xyz_full[i+1][2]+xyz_rnd[2];
                int_coor_array[i][0]=xyz_full[i][0]/4;
                int_coor_array[i][1]=xyz_full[i][1]/4;
                int_coor_array[i][2]=xyz_full[i][2]/4;
                clash_score=getClash(xyz_full,int_coor_array,miss_list,i);
                if (clash_score<=0) break;
                if (Ntrial==0) clash_score_best=clash_score;
                cerr<<"N resi="<<i<<";trial="<<Ntrial<<";clash="<<clash_score<<endl;
                if (clash_score<=clash_score_best)
                {
                    coor_prev[0]    =xyz_full[i][0];
                    coor_prev[1]    =xyz_full[i][1];
                    coor_prev[2]    =xyz_full[i][2];
                    int_coor_prev[0]=int_coor_array[i][0];
                    int_coor_prev[1]=int_coor_array[i][1];
                    int_coor_prev[2]=int_coor_array[i][2];
                }
                else if (clash_score>clash_score_best)
                {
                    xyz_full[i][0]=coor_prev[0];
                    xyz_full[i][1]=coor_prev[1];
                    xyz_full[i][2]=coor_prev[2];
                    int_coor_array[i][0]=int_coor_prev[0];
                    int_coor_array[i][1]=int_coor_prev[1];
                    int_coor_array[i][2]=int_coor_prev[2];
                }
            }
            miss_list[i]=false;
        }
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
        fillLoop(rstart,rend,xyz_full,int_coor_array,miss_list,d);
    }
   
    /* clean up */
    vector<float>().swap(coor_prev);
    vector<float>().swap(xyz_rnd);
    vector<vector<int> >().swap(int_coor_array);
    vector<int>().swap(int_coor_prev);
    return;
}
#endif

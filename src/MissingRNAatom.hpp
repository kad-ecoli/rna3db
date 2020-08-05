#ifndef MissingRNAatom_HPP
#define MissingRNAatom_HPP 1

#include <iostream>
#include "PDBParser.hpp"
#include "Superpose.hpp"
#include "IdealRNA.hpp"

/* return true if residue is standard nucleotide
 * return false otherwise */
bool checkMissingRNAatom(ResidueUnit &residue, 
    map<string, map<string,vector<float> > >&ideal_pdb,
    vector<string> &missing_atom_vec,
    vector<size_t> &nonstd_atom_vec, const int option=0)
{
    missing_atom_vec.clear();
    nonstd_atom_vec.clear();
    if (ideal_pdb.count(residue.resn)==0)
    {
        cout<<residue.resn<<residue.resi<<residue.icode
            <<"\tnon-standard residue"<<endl;
        return false;
    }
    size_t a,i,j;
    for (a=0;a<residue.atoms.size();a++)
        if (ideal_pdb[residue.resn].count(residue.atoms[a].name)==0)
            nonstd_atom_vec.push_back(a);
    if (nonstd_atom_vec.size())
    {
        cout<<residue.resn<<residue.resi<<residue.icode<<"\t"
            <<nonstd_atom_vec.size()<<" non-standard atoms\t";
        for (i=0;i<nonstd_atom_vec.size();i++)
            cout<<residue.atoms[nonstd_atom_vec[i]].name<<' ';
        cout<<endl;
        if (option)
        {
            for (i=0;i<nonstd_atom_vec.size();i++)
                residue.atoms[nonstd_atom_vec[i]].name.clear();
            for (i=0;i<residue.atoms.size();i++)
            {
                if (residue.atoms[i].name.size()) continue;
                for (j=i+1;j<residue.atoms.size();j++)
                {
                    if (residue.atoms[j].name.size()==0) continue;
                    residue.atoms[i].name  =residue.atoms[j].name;
                    residue.atoms[i].xyz[0]=residue.atoms[j].xyz[0];
                    residue.atoms[i].xyz[1]=residue.atoms[j].xyz[1];
                    residue.atoms[i].xyz[2]=residue.atoms[j].xyz[2];
                    residue.atoms[j].name.clear();
                    break;
                }
            }
            for (i=0;i<nonstd_atom_vec.size();i++) residue.atoms.pop_back();
        }
    }
    
    map<string,vector<float> >::iterator it;
    bool found;
    for (it=ideal_pdb[residue.resn].begin(); it!=ideal_pdb[residue.resn].end(); it++)
    {
        found=false;
        for (a=0;a<residue.atoms.size();a++)
        {
            if (residue.atoms[a].name==it->first)
            {
                found=true;
                break;
            }
        }
        if (found==false) missing_atom_vec.push_back(it->first);
    }
    if (missing_atom_vec.size())
    {
        cout<<residue.resn<<residue.resi<<residue.icode<<"\t"
            <<missing_atom_vec.size()<<" missing atoms\t";
        for (a=0;a<missing_atom_vec.size();a++)
            cout<<missing_atom_vec[a]<<' ';
        cout<<endl;
    }
    return true;
}

void fillMissingRNAatom(ResidueUnit &residue,
    map<string, map<string,vector<float> > >&ideal_pdb,
    vector<string>missing_atom_vec)
{
    vector<float> tmp(3,0);
    vector<vector<float> > xyz_list1(residue.atoms.size(),tmp);
    vector<vector<float> > xyz_list2(residue.atoms.size(),tmp);
    vector<vector<float> > RotMatix;  // U
    vector<float> TranVect;  // t

    size_t a;
    for (a=0;a<residue.atoms.size();a++)
    {
        xyz_list1[a][0]=ideal_pdb[residue.resn][residue.atoms[a].name][0];
        xyz_list1[a][1]=ideal_pdb[residue.resn][residue.atoms[a].name][1];
        xyz_list1[a][2]=ideal_pdb[residue.resn][residue.atoms[a].name][2];

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
    for (a=0;a<missing_atom_vec.size();a++)
    {
        atom.name=missing_atom_vec[a];
        ChangeCoor(ideal_pdb[residue.resn][atom.name],
            RotMatix, TranVect, atom.xyz);
        residue.atoms.push_back(atom);
    }

    atom.name.clear();
    atom.xyz.clear();
    TranVect.clear();
    for (a=0;a<3;a++) RotMatix[a].clear();
    RotMatix.clear();
    return;
}

bool MissingRNAatom(ChainUnit &chain, 
    map<string, map<string,vector<float> > >&ideal_pdb,
    const int option=0)
{
    vector<string> missing_atom_vec;
    vector<size_t> nonstd_atom_vec;
    size_t r,i,j;
    vector<size_t> remove_residues_vec;
    for (r=0;r<chain.residues.size();r++)
    {
        checkMissingRNAatom(chain.residues[r], ideal_pdb,
            missing_atom_vec, nonstd_atom_vec, option);
        if (option==2 && missing_atom_vec.size())
        {
            cout<<"remove residue "
                <<chain.residues[r].resn
                <<chain.residues[r].resi
                <<chain.residues[r].icode
                <<" with "
                <<missing_atom_vec.size()
                <<" missing atoms"<<endl;
            chain.residues[r].resn.clear();
            remove_residues_vec.push_back(r);
        }
        if (option>=3 && chain.residues[r].atoms.size()>=3 &&
            missing_atom_vec.size())
        {
            cout<<"fill residue "
                <<chain.residues[r].resn
                <<chain.residues[r].resi
                <<chain.residues[r].icode
                <<" with "
                <<missing_atom_vec.size()
                <<" missing atoms"<<endl;
            fillMissingRNAatom(chain.residues[r], ideal_pdb, missing_atom_vec);
            missing_atom_vec.clear();
        }
        if (option==3 && chain.residues[r].atoms.size()<3)
        {
            cout<<"failed to fill residue "
                <<chain.residues[r].resn
                <<chain.residues[r].resi
                <<chain.residues[r].icode
                <<" with "
                <<chain.residues[r].atoms.size()
                <<" atoms and "
                <<missing_atom_vec.size()
                <<" missing atoms"<<endl;
        }
        if (option==4 && chain.residues[r].atoms.size()<3)
        {
            cout<<"remove residue "
                <<chain.residues[r].resn
                <<chain.residues[r].resi
                <<chain.residues[r].icode
                <<" with "
                <<chain.residues[r].atoms.size()
                <<" atoms"<<endl;
            chain.residues[r].resn.clear();
            remove_residues_vec.push_back(r);
        }
    }
    if (remove_residues_vec.size())
    {
        for (i=0;i<chain.residues.size();i++)
        {
            if (chain.residues[i].resn.size()) continue;
            for (j=i+1;j<chain.residues.size();j++)
            {
                if (chain.residues[j].resn.size()==0) continue;
                chain.residues[i].het  =chain.residues[j].het;
                chain.residues[i].resi =chain.residues[j].resi;
                chain.residues[i].icode=chain.residues[j].icode;
                chain.residues[i].resn =chain.residues[j].resn;
                chain.residues[i].atoms=chain.residues[j].atoms;
                chain.residues[j].resn.clear();
                break;
            }
        }
        for (i=0;i<remove_residues_vec.size();i++) chain.residues.pop_back();
    }
    missing_atom_vec.clear();
    nonstd_atom_vec.clear();
    remove_residues_vec.clear();
    return true;
}

bool MissingRNAatom(ModelUnit &pdb_entry, 
    map<string, map<string,vector<float> > >&ideal_pdb,
    const int option=0)
{
    size_t c;
    for (c=0;c<pdb_entry.chains.size();c++)
        MissingRNAatom(pdb_entry.chains[c], ideal_pdb, option);
    return true;
}

#endif

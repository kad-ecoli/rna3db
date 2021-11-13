const char* docstring=""
"MissingResidue input.fasta input.pdb output.pdb\n"
"    Find and fill missing atoms in input.pdb based on input fasta\n"
"    Output filled model to output.pdb\n"
"\n"
"Input:"
"    input.fasta - input sequence.\n"
"    input.pdb   - input PDB. Should include CA for a protein.\n"
"                  Should include P and C4' for RNA.\n"
"                  Residue index must be consistent with that in input.fasta\n"
"                  Only the first chain is considered\n"
;

#include "MissingResidue.hpp"

int main(int argc,char **argv)
{
    if (argc<4)
    {
        cerr<<docstring<<endl;
        return 1;
    }
    string infasta=argv[1];
    string infile =argv[2];
    string outfile=argv[3];

    srand(111111);

    int atomic_detail=2;
    int allowX=0;
    ModelUnit pdb_entry=read_pdb_structure(infile.c_str(),atomic_detail,allowX);
    
    int mol_opt=check_mol(pdb_entry.chains[0]);
    string sequence="";
    readOneSequenceFromFasta(infasta,sequence,mol_opt);

    vector<vector<float> > xyz_full;
    vector<bool> miss_list;
    existCoor2vec(sequence,pdb_entry.chains[0],xyz_full,miss_list,mol_opt);
    MissingResidue(xyz_full,miss_list);
    appendNewResidue(sequence,pdb_entry.chains[0],xyz_full,miss_list,mol_opt);
    write_pdb_structure(outfile.c_str(),pdb_entry.chains[0]);
    vector<ChainUnit>().swap(pdb_entry.chains);
    return 0;
}

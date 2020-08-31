const char* docstring=
"NAbuilder sequence rama.txt seed out.pdb\n"
"    reconstruct 3D structure from RNA sequence\n"
"    and pseudo torsion angle table 'pseudo.txt'.\n"
"    and output structure to out.pdb\n"
"\n"
"NAbuilder 10 random 147842874 out.pdb\n"
"    generate a 10-nucleotide random structure with random sequence\n"
"    using 147842874 as seed\n"
"\n"
"Input\n"
"    'sequence' - lower case input sequence with nucleotide type a c g u,\n"
"                 such as 'cgggcauccg'. it can also be an integer, in which\n"
"                 case the sequence is randomly generated with given length\n"
"    'rama.txt' - two-column table for eta and theta. if 'random',\n"
"                 use random-walk to generate the structure\n"
"    'seed'     - random seed. this is only useful if either sequence or\n"
"                 torsion angle is random\n"
;

#include <iostream>
#include <vector>
#include <string>

#include "NAbuilder.hpp"

using namespace std;

void readRamaTable(const char* filename, vector<vector<float> >&rama_table)
{
    ifstream fp(filename, ios::in);
    int use_stdin=(strcmp(filename,"-")==0);
    if (!fp && !use_stdin)
    {
        cerr<<"ERROR! Cannot read file "<<filename<<endl;
        return;
    }

    int r=0;
    float theta, eta;
    while (use_stdin?cin.good():fp.good())
    {
        if (use_stdin) cin>>eta>>theta;
        else            fp>>eta>>theta;
        if (theta>=-180) rama_table[r][0]=eta;
        if (eta  >=-180) rama_table[r][1]=theta;
        //cout<<setiosflags(ios::fixed)<<setprecision(2)<<setw(7)
            //<<rama_table[r][0]<<' '<<setw(7)<<rama_table[r][1]<<endl;
        if (++r >= rama_table.size()) break;
    }
}

int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc<2)
    {
        cerr<<docstring;
        return 0;
    }

    if (argc>3) srand(atol(argv[3]));
    else        srand(time(0));

    /* read query sequence */
    string sequence=argv[1];
    int L=sequence.size();
    if ('0'<=sequence[0] && sequence[0]<='9')
    {
        L=atoi(argv[1]);
        sequence.clear();
        randomRNAseq(sequence,L);
        cout<<sequence<<endl;
    }

    /* read backbone torsion table */
    vector<float> tmp_array(2,180);
    tmp_array[0]= 168.46; //   eta: C4'[-1]-P-C4'-P[+1]
    tmp_array[1]=-147.97; // theta: P-C4'-P[+1]-C4'[+1]
    vector<vector<float> >rama_table(L,tmp_array);
    if (argc>2)
    {
        string ramaFileName=argv[2];
        if (ramaFileName!="random") 
            readRamaTable(ramaFileName.c_str(),rama_table);
        else randomRamaTable(rama_table);
        ramaFileName.clear();
    }

    //for (int r=0;r<L;r++) cout<<setiosflags(ios::fixed)<<setprecision(2)
            //<<rama_table[r][0]<<'\t'<<rama_table[r][1]<<endl;
    
    /* make structure */
    tmp_array[0]=tmp_array[1]=0;
    tmp_array.push_back(0);
    vector<vector<float> > P_coord(L+1,tmp_array);
    vector<vector<float> > C4_coord(L,tmp_array);
    vector<vector<float> > C1_coord(L,tmp_array);
    make_chain(P_coord, C4_coord, C1_coord, rama_table);
    ChainUnit chain;
    convert_chain(chain, sequence, P_coord, C4_coord, C1_coord);

    /* fill atoms */
    map<string, map<string,vector<float> > >simple_rna;
    map<string,vector<string> >missing_rna;
    parse_simple_rna(simple_rna, missing_rna);
    for (int r=0;r<L;r++)
        fillSimpleMissingRNAatom(chain.residues[r], simple_rna, missing_rna);
    write_pdb_structure((argc>4?argv[4]:"-"),chain);

    /* clean up */
    vector<float>().swap(tmp_array);
    vector<vector<float> >().swap(P_coord);
    vector<vector<float> >().swap(C4_coord);
    vector<vector<float> >().swap(C1_coord);
    vector<vector<float> >().swap(rama_table);
    vector<ResidueUnit>().swap(chain.residues);
    string().swap(sequence);
    map<string, map<string,vector<float> > >().swap(simple_rna);
    map<string,vector<string> >().swap(missing_rna);
    return 0;
}

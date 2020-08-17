const char* docstring=""
"newChainID input.pdb output.pdb chainID\n"
"    replace chain ID for ATOM and HETATM line of input PDB file to new\n"
"    chainID in output PDB file. default new chainID is '_' (empty)\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>

using namespace std;

void newChainID(const string infile="-", const string outfile="-", char ChainID='_')
{
    if (ChainID=='_') ChainID=' ';
    ifstream fp_in;
    ofstream fp_out;
    if (infile!="-") fp_in.open(infile.c_str(),ios::in);
    if (outfile!="-") fp_out.open(outfile.c_str(),ofstream::out);
    string txt,line;
    while ((infile!="-")?fp_in.good():cin.good())
    {
        if (infile!="-") getline(fp_in,line);
        else getline(cin,line);

        if      (line.substr(0,6)=="DBREF ")
            if (outfile!="-") fp_out<<line.substr(0,12)<<ChainID<<line.substr(13)<<endl;
            else                cout<<line.substr(0,12)<<ChainID<<line.substr(13)<<endl;
        else if (line.substr(0,6)=="SEQADV")
            if (outfile!="-") fp_out<<line.substr(0,16)<<ChainID<<line.substr(17)<<endl;
            else                cout<<line.substr(0,16)<<ChainID<<line.substr(17)<<endl;
        else if (line.substr(0,6)=="SEQRES")
            if (outfile!="-") fp_out<<line.substr(0,11)<<ChainID<<line.substr(12)<<endl;
            else                cout<<line.substr(0,11)<<ChainID<<line.substr(12)<<endl;
        else if (line.substr(0,6)=="MODRES")
            if (outfile!="-") fp_out<<line.substr(0,16)<<ChainID<<line.substr(17)<<endl;
            else                cout<<line.substr(0,16)<<ChainID<<line.substr(17)<<endl;
        else if (line.substr(0,6)=="HET   ")
            if (outfile!="-") fp_out<<line.substr(0,12)<<ChainID<<line.substr(13)<<endl;
            else                cout<<line.substr(0,12)<<ChainID<<line.substr(13)<<endl;
        else if (line.substr(0,6)=="HET   ")
            if (outfile!="-") fp_out<<line.substr(0,12)<<ChainID<<line.substr(13)<<endl;
            else                cout<<line.substr(0,12)<<ChainID<<line.substr(13)<<endl;
        else if (line.substr(0,6)=="HELIX ")
            if (outfile!="-") fp_out<<line.substr(0,19)<<ChainID<<line.substr(20)<<endl;
            else                cout<<line.substr(0,19)<<ChainID<<line.substr(20)<<endl;
        else if (line.substr(0,6)=="SHEET ")
            if (outfile!="-") fp_out<<line.substr(0,21)<<ChainID<<line.substr(22)<<endl;
            else                cout<<line.substr(0,21)<<ChainID<<line.substr(22)<<endl;
        else if (line.substr(0,6)=="SSBOND")
            if (outfile!="-") fp_out<<line.substr(0,15)<<ChainID<<line.substr(16)<<endl;
            else                cout<<line.substr(0,15)<<ChainID<<line.substr(16)<<endl;
        else if (line.substr(0,6)=="SITE  ")
            if (outfile!="-") fp_out<<line.substr(0,22)<<ChainID<<line.substr(23)<<endl;
            else                cout<<line.substr(0,22)<<ChainID<<line.substr(23)<<endl;
        else if (line.substr(0,6)=="ATOM  ")
            if (outfile!="-") fp_out<<line.substr(0,21)<<ChainID<<line.substr(22)<<endl;
            else                cout<<line.substr(0,21)<<ChainID<<line.substr(22)<<endl;
        else if (line.substr(0,6)=="HETATM")
            if (outfile!="-") fp_out<<line.substr(0,21)<<ChainID<<line.substr(22)<<endl;
            else                cout<<line.substr(0,21)<<ChainID<<line.substr(22)<<endl;
        else if (line.substr(0,6)=="TER   " && line.size()>=21)
            if (outfile!="-") fp_out<<line.substr(0,21)<<ChainID<<line.substr(22)<<endl;
            else                cout<<line.substr(0,21)<<ChainID<<line.substr(22)<<endl;
        else
            if (outfile!="-") fp_out<<line<<endl;
            else                cout<<line<<endl;
    }
    fp_in.close();
    fp_out.close();
    return;
}

int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    string infile=argv[1];
    string outfile=(argc<=2)?"-":argv[2];
    char   chainID=(argc<=3)?'_':argv[3][0];
    newChainID(infile,outfile,chainID);
    return 0;
}

#!/bin/bash
# Download standard conformation of nucleotides
# create IdealRNA.hpp

__file__=`readlink -e $0`
ligandexpodir=`dirname $__file__`
cd $ligandexpodir
ligandlist="A
U
C
G
DA
DT
DC
DG"
formatlist=".cif
.xml
_model.sdf
_ideal.sdf
_model.pdb
_ideal.pdb
_D3L1.gif
_D3L3.gif
"

for ligand in $ligandlist;do
    echo $ligand
    mkdir -p $ligand
    folder=`echo $ligand|cut -c1`
    for format in $formatlist;do
        wget -q http://ligand-expo.rcsb.org/reports/$folder/$ligand/$ligand$format -O $ligand/$ligand$format
    done
done

./createIdealRNA.py

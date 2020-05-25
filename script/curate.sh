#!/bin/bash
bindir=`dirname $(readlink -e $0)`
rootdir=`dirname $bindir`

echo "database path $rootdir"
cd $rootdir

#echo "Make MODRES"
#zcat $rootdir/pdb/data/structures/all/pdb/*gz |grep "^MODRES" > \
     #$rootdir/pdb/derived_data/MODRES
#$bindir/parseMODRES.py $rootdir/pdb/derived_data/MODRES $bindir/MODRES_dicts.py

echo "Convert PDB entries to fasta"
cd       $rootdir/pdb/data/structures/all/pdb/
mkdir -p $rootdir/cull

$bindir/pdb2fasta.py -dir=$rootdir/pdb/data/structures/all/pdb/ \
    $rootdir/pdb/derived_data/na_chain.list -suffix=.pdb.gz \
    -allowX=false -mol=rna -PERMISSIVE=TER | cut -f1 -d: \
    > $rootdir/cull/pdb_atom.fasta

echo "Cull sequence"
$bindir/SortFastaWithResolution.py $rootdir/pdb/derived_data/index/resolu.idx $rootdir/cull/pdb_atom.fasta $rootdir/cull/pdb_atom.sort
resolu_list="1.5
2.0
2.5
3.0
3.5
4.0
20.0
all"

cdhitest="$bindir/cd-hit-est -l 9 -r 0 -M 5000 -g 1"
for resolu in $resolu_list;do
    $cdhitest -i $rootdir/cull/pdb_atom.sort.${resolu}.fasta     -n 10 -o $rootdir/cull/pdb_atom.sort.${resolu}_c1.0_s1.0 -c 1.0 -s 1.0
    $cdhitest -i $rootdir/cull/pdb_atom.sort.${resolu}_c1.0_s1.0 -n  9 -o $rootdir/cull/pdb_atom.sort.${resolu}_c0.9_s0.9 -c 0.9 -s 0.9
    $cdhitest -i $rootdir/cull/pdb_atom.sort.${resolu}_c0.9_s0.9 -n  5 -o $rootdir/cull/pdb_atom.sort.${resolu}_c0.8_s0.8 -c 0.8 -s 0.8

    $cdhitest -i $rootdir/cull/pdb_atom.sort.${resolu}_c1.0_s1.0 -n 10 -o $rootdir/cull/pdb_atom.sort.${resolu}_c1.0_s0.0 -c 1.0
    $cdhitest -i $rootdir/cull/pdb_atom.sort.${resolu}_c0.9_s0.9 -n  9 -o $rootdir/cull/pdb_atom.sort.${resolu}_c0.9_s0.0 -c 0.9
    $cdhitest -i $rootdir/cull/pdb_atom.sort.${resolu}_c0.8_s0.8 -n  5 -o $rootdir/cull/pdb_atom.sort.${resolu}_c0.8_s0.0 -c 0.8
done

$bindir/makeCullTable.py $rootdir/cull

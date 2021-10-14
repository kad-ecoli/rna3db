#!/bin/bash
bindir=`dirname $(readlink -e $0)`
rootdir=`dirname $bindir`
dimerdir=$rootdir/dimer

echo "dimer database path $dimerdir"
if [ ! -d "$dimerdir" ];then
    mkdir -p $dimerdir
fi
cd $dimerdir
$bindir/makeDimerList.py $rootdir/cull/pdb_atom.fasta $dimerdir/list.multichain $rootdir/pdb/data/structures/all/pdb $dimerdir/list.nodimer $dimerdir/pdb_atom.fasta

echo "Cull sequence"
$bindir/SortDimerFastaWithResolution.py $rootdir/pdb/derived_data/index/resolu.idx $dimerdir/pdb_atom.fasta $dimerdir/pdb_atom.sort

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
    $cdhitest -i $dimerdir/pdb_atom.sort.${resolu}.fasta     -n 10 -o $dimerdir/pdb_atom.sort.${resolu}_c1.0_s1.0 -c 1.0 -s 1.0
    $cdhitest -i $dimerdir/pdb_atom.sort.${resolu}_c1.0_s1.0 -n  9 -o $dimerdir/pdb_atom.sort.${resolu}_c0.9_s0.9 -c 0.9 -s 0.9
    $cdhitest -i $dimerdir/pdb_atom.sort.${resolu}_c0.9_s0.9 -n  5 -o $dimerdir/pdb_atom.sort.${resolu}_c0.8_s0.8 -c 0.8 -s 0.8

    $cdhitest -i $dimerdir/pdb_atom.sort.${resolu}_c1.0_s1.0 -n 10 -o $dimerdir/pdb_atom.sort.${resolu}_c1.0_s0.0 -c 1.0
    $cdhitest -i $dimerdir/pdb_atom.sort.${resolu}_c0.9_s0.9 -n  9 -o $dimerdir/pdb_atom.sort.${resolu}_c0.9_s0.0 -c 0.9
    $cdhitest -i $dimerdir/pdb_atom.sort.${resolu}_c0.8_s0.8 -n  5 -o $dimerdir/pdb_atom.sort.${resolu}_c0.8_s0.0 -c 0.8
done

$bindir/makeCullTable.py $dimerdir


echo "Clean PDB"
resolu_cs_list="all_c1.0_s1.0"

for resolu_cs in $resolu_cs_list;do
    echo $resolu_cs
    mkdir -p $dimerdir/$resolu_cs/DSSR
    grep '>' $dimerdir/pdb_atom.sort.$resolu_cs|cut -f1|sed 's/>//g' > $dimerdir/$resolu_cs/list
    cd $dimerdir/$resolu_cs

    echo "remove files not listed by $PWD/list"
    for target in `ls *.pdb|cut -f1 -d.`;do
        GREP=`grep $target list`
	if [ -z "$GREP" ];then
	    echo rm $target.pdb
	    rm $target.pdb
	    rm DSSR/$target.dssr
	fi
    done
    
    echo "update $PWD/$list.new"
    for target in `cat list`;do
        if [ ! -s "${target}.pdb" ];then
	    echo $target
	fi
    done > list.new
    
    echo "generate $PWD/*.pdb"
    $bindir/cleanDimerPdb.py -dir=$rootdir/pdb/data/structures/all/pdb/ list.new -suffix=.pdb.gz .pdb
    rm list.new

    echo "generate $PWD/list.atomic"
    $bindir/AtomsPerResidue.py -dir=./ list -suffix=.pdb|grep -P "^[:\w]+\s[2]\d\.\d{2}$"|cut -f1 > list.atomic

    echo "generate $PWD/DSSR/*.dssr"
    for target in `cat list`;do
	    if [ ! -s "$dimerdir/$resolu_cs/DSSR/$target.dssr" ];then
	    echo    x3dna-dssr -i=$target.pdb --pair-only -o=DSSR/$target.dssr
	    $bindir/x3dna-dssr -i=$target.pdb --pair-only -o=DSSR/$target.dssr
	fi
    done
    
    echo "fill missing atoms"
    cd $rootdir/dimer/$resolu_cs
    for target in `cat list`;do
	$bindir/MissingRNAatom ${target}.pdb ${target}.pdb 5
    done
    $bindir/clstr2txt.py $rootdir/dimer/pdb_atom.sort.${resolu_cs}.clstr $rootdir/dimer/pdb_atom.sort.${resolu_cs}.txt
    $bindir/pdb2xyz -dir ./ list -suffix .pdb -atom " C3'" > C3.xyz
    $bindir/pdb2xyz -dir ./ list -suffix .pdb -atom " C4'" > C4.xyz
    $bindir/pdb2xyz -dir ./ list -suffix .pdb -atom " C5'" > C5.xyz
    $bindir/pdb2xyz -dir ./ list -suffix .pdb -atom " O5'" > O5.xyz
    $bindir/pdb2xyz -dir ./ list -suffix .pdb -atom " O3'" > O3.xyz
    $bindir/pdb2xyz -dir ./ list -suffix .pdb -atom " P  " >  P.xyz
    cd $rootdir/dimer
    tar -cjvf ${resolu_cs}.tar.bz2 ${resolu_cs}/list* ${resolu_cs}/*.pdb ${resolu_cs}/*.xyz pdb_atom.sort.${resolu_cs}.txt
done

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
    $cdhitest -i $dimerdir/pdb_atom.sort.${resolu}.fasta     -n  8 -o $dimerdir/pdb_atom.sort.${resolu}_c1.0_s1.0 -c 1.0 -s 1.0
    $cdhitest -i $dimerdir/pdb_atom.sort.${resolu}_c1.0_s1.0 -n  7 -o $dimerdir/pdb_atom.sort.${resolu}_c0.9_s0.9 -c 0.9 -s 0.9
    $cdhitest -i $dimerdir/pdb_atom.sort.${resolu}_c0.9_s0.9 -n  4 -o $dimerdir/pdb_atom.sort.${resolu}_c0.8_s0.8 -c 0.8 -s 0.8

    $cdhitest -i $dimerdir/pdb_atom.sort.${resolu}_c1.0_s1.0 -n  8 -o $dimerdir/pdb_atom.sort.${resolu}_c1.0_s0.0 -c 1.0
    $cdhitest -i $dimerdir/pdb_atom.sort.${resolu}_c0.9_s0.9 -n  7 -o $dimerdir/pdb_atom.sort.${resolu}_c0.9_s0.0 -c 0.9
    $cdhitest -i $dimerdir/pdb_atom.sort.${resolu}_c0.8_s0.8 -n  4 -o $dimerdir/pdb_atom.sort.${resolu}_c0.8_s0.0 -c 0.8
done

$bindir/makeCullTable.py $dimerdir

echo "Clean PDB"
resolu_cs_list="all_c1.0_s1.0"

for resolu_cs in $resolu_cs_list;do
    echo $resolu_cs
    mkdir -p $dimerdir/$resolu_cs/DSSR
    mkdir -p $dimerdir/$resolu_cs/RNAfold
    grep '>' $dimerdir/pdb_atom.sort.$resolu_cs|cut -f1|sed 's/>//g' > $dimerdir/$resolu_cs/list
    cd $dimerdir/$resolu_cs

    echo "remove files not listed by $PWD/list"
    for target in `ls *.pdb|cut -f1 -d.`;do
        GREP=`grep $target list`
        if [ -z "$GREP" ];then
            echo rm $target.pdb
            rm $target.pdb
            rm DSSR/$target.dssr
            rm DSSR/$target.dbn
            rm DSSR/$target.sto
            rm DSSR/$target.cm
	    rm RNAfold/${target}_dp.ps.gz
        fi
    done
    
    echo "generate $PWD/*.pdb"
    $bindir/cleanDimerPdb.py -dir=$rootdir/pdb/data/structures/all/pdb/ list -suffix=.pdb.gz .pdb

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
    cut -f1 $rootdir/dimer/pdb_atom.sort.${resolu_cs} | $bindir/fastaOneLine - $rootdir/dimer/$resolu_cs/list.fasta
    tar -cjvf ${resolu_cs}.tar.bz2 ${resolu_cs}/list* ${resolu_cs}/*.pdb ${resolu_cs}/*.xyz pdb_atom.sort.${resolu_cs}.txt ${resolu_cs}/list.fasta
    
    echo "generate dot bracket"
    cd $rootdir/dimer/$resolu_cs/DSSR
    for target in `cat $rootdir/dimer/$resolu_cs/list.atomic`;do
        if [ -s "$target.dbn" ];then
            continue
        fi
        $bindir/x3dna-dssr -i=$rootdir/dimer/$resolu_cs/$target.pdb
        Ldbn=`tail -1 dssr-2ndstrs.dbn|sed 's/&//g'|wc -c`
        Lseq=`grep -A1 ">$target$" $rootdir/dimer/$resolu_cs/list.fasta|tail -1 |wc -c`
        if [ $Lseq -eq $Ldbn ];then
            tail -1 dssr-2ndstrs.dbn | sed 's/&//g' | sed 's/[A-Za-z{}\[<>]/./g' | sed 's/]/./g' > $target.dbn
        fi
        rm dssr-*
    done

    echo "generate _dp.ps.gz"
    cd $rootdir/dimer/$resolu_cs/RNAfold
    for target in `cat $rootdir/dimer/$resolu_cs/list`;do
        if [ -s "${target}_dp.ps.gz" ];then
            continue
        fi
	echo ">$target" > $target.fasta
	grep -PA2 "^>$target\b" $rootdir/dimer/pdb_atom.sort.$resolu_cs|grep -v '>'|tr "[:lower:]" "[:upper:]" >> $target.fasta
        $bindir/RNAfold -p  $target.fasta
	mv `echo ${target}_dp.ps|sed 's/:/_/g'` ${target}_dp.ps
	gzip ${target}_dp.ps
	rm $target.fasta `echo ${target}_ss.ps|sed "s/:/_/g"`
    done

    echo "re-copy $PWD/*.pdb"
    cd $rootdir/dimer/$resolu_cs
    $bindir/cleanDimerPdb.py -dir=$rootdir/pdb/data/structures/all/pdb/ list -suffix=.pdb.gz .pdb
    cd $rootdir/dimer/$resolu_cs/DSSR
    for target in `cat $rootdir/dimer/$resolu_cs/list`;do
        if [ -s "$target.dbn" ];then
            continue
        fi
        $bindir/cssr $rootdir/dimer/$resolu_cs/$target.pdb 1 | sed 's/[\[{}()]/./g'|sed 's/]/./g'  > $target.dbn
    done
done

echo "make smaller subsets"
resolu_cs_list="all_c1.0_s1.0 all_c0.9_s0.9 all_c0.8_s0.8 all_c0.8_s0.0"

for resolu_cs in $resolu_cs_list;do
    echo $resolu_cs
    mkdir -p $rootdir/dimer/$resolu_cs
    grep '>' $rootdir/dimer/pdb_atom.sort.$resolu_cs|cut -f1|sed 's/>//g' > $rootdir/dimer/$resolu_cs/list
    $bindir/clstr2txt.py $rootdir/dimer/pdb_atom.sort.${resolu_cs}.clstr $rootdir/dimer/pdb_atom.sort.${resolu_cs}.txt
    cut -f1 $rootdir/dimer/pdb_atom.sort.${resolu_cs} > $rootdir/dimer/$resolu_cs/list.fasta
done

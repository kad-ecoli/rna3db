#!/bin/bash
bindir=`dirname $(readlink -e $0)`
rootdir=`dirname $bindir`

echo "database path $rootdir"
cd $rootdir

echo "Make MODRES"
#zcat $rootdir/pdb/data/structures/all/pdb/*/*.pdb.gz |grep "^MODRES" > \
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
    $cdhitest -i $rootdir/cull/pdb_atom.sort.${resolu}.fasta     -n  8 -o $rootdir/cull/pdb_atom.sort.${resolu}_c1.0_s1.0 -c 1.0 -s 1.0
    $cdhitest -i $rootdir/cull/pdb_atom.sort.${resolu}_c1.0_s1.0 -n  7 -o $rootdir/cull/pdb_atom.sort.${resolu}_c0.9_s0.9 -c 0.9 -s 0.9
    $cdhitest -i $rootdir/cull/pdb_atom.sort.${resolu}_c0.9_s0.9 -n  4 -o $rootdir/cull/pdb_atom.sort.${resolu}_c0.8_s0.8 -c 0.8 -s 0.8

    $cdhitest -i $rootdir/cull/pdb_atom.sort.${resolu}_c1.0_s1.0 -n  8 -o $rootdir/cull/pdb_atom.sort.${resolu}_c1.0_s0.0 -c 1.0
    $cdhitest -i $rootdir/cull/pdb_atom.sort.${resolu}_c0.9_s0.9 -n  7 -o $rootdir/cull/pdb_atom.sort.${resolu}_c0.9_s0.0 -c 0.9
    $cdhitest -i $rootdir/cull/pdb_atom.sort.${resolu}_c0.8_s0.8 -n  4 -o $rootdir/cull/pdb_atom.sort.${resolu}_c0.8_s0.0 -c 0.8
done

$bindir/makeCullTable.py $rootdir/cull

echo "Clean PDB"
resolu_cs_list="all_c1.0_s1.0"

for resolu_cs in $resolu_cs_list;do
    echo $resolu_cs
    mkdir -p $rootdir/cull/$resolu_cs/DSSR
    mkdir -p $rootdir/cull/$resolu_cs/RNAfold
    grep '>' $rootdir/cull/pdb_atom.sort.$resolu_cs|cut -f1|sed 's/>//g' > $rootdir/cull/$resolu_cs/list
    cd $rootdir/cull/$resolu_cs

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
    $bindir/clean_pdb.py -dir=$rootdir/pdb/data/structures/all/pdb/ list -suffix=.pdb.gz .pdb -StartIndex=1 -NewChainID=_

    echo "generate $PWD/list.atomic"
    $bindir/AtomsPerResidue.py -dir=./ list -suffix=.pdb|grep -P "^\w+\s[2]\d\.\d{2}$"|cut -f1 > list.atomic

    echo "generate $PWD/DSSR/*.dssr"
    for target in `cat list`;do
            if [ ! -s "$rootdir/cull/$resolu_cs/DSSR/$target.dssr" ];then
            echo    x3dna-dssr -i=$target.pdb --pair-only -o=DSSR/$target.dssr
            $bindir/x3dna-dssr -i=$target.pdb --pair-only -o=DSSR/$target.dssr
        fi
    done

    echo "fill missing atoms"
    cd $rootdir/cull/$resolu_cs
    for target in `cat list`;do
        $bindir/MissingRNAatom ${target}.pdb ${target}.pdb 5
    done
    $bindir/clstr2txt.py $rootdir/cull/pdb_atom.sort.${resolu_cs}.clstr $rootdir/cull/pdb_atom.sort.${resolu_cs}.txt
    $bindir/pdb2xyz -dir ./ list -suffix .pdb -atom " C3'" > C3.xyz
    $bindir/pdb2xyz -dir ./ list -suffix .pdb -atom " C4'" > C4.xyz
    $bindir/pdb2xyz -dir ./ list -suffix .pdb -atom " C5'" > C5.xyz
    $bindir/pdb2xyz -dir ./ list -suffix .pdb -atom " O5'" > O5.xyz
    $bindir/pdb2xyz -dir ./ list -suffix .pdb -atom " O3'" > O3.xyz
    $bindir/pdb2xyz -dir ./ list -suffix .pdb -atom " P  " >  P.xyz
    cd $rootdir/cull
    cut -f1 $rootdir/cull/pdb_atom.sort.${resolu_cs} > $rootdir/cull/$resolu_cs/list.fasta
    tar -cjvf ${resolu_cs}.tar.bz2 ${resolu_cs}/list* ${resolu_cs}/*.pdb ${resolu_cs}/*.xyz pdb_atom.sort.${resolu_cs}.txt ${resolu_cs}/list.fasta
    
    echo "generate dot bracket"
    cd $rootdir/cull/$resolu_cs/DSSR
    for target in `cat $rootdir/cull/$resolu_cs/list.atomic`;do
        if [ -s "$target.dbn" ];then
            continue
        fi
        $bindir/x3dna-dssr -i=$rootdir/cull/$resolu_cs/$target.pdb
        Ldbn=`tail -1 dssr-2ndstrs.dbn|sed 's/&//g'|wc -c`
        Lseq=`grep -A1 ">$target$" $rootdir/cull/$resolu_cs/list.fasta|tail -1 |wc -c`
        if [ $Lseq -eq $Ldbn ];then
            tail -1 dssr-2ndstrs.dbn | sed 's/&//g' | sed 's/[A-Za-z{}\[<>]/./g' | sed 's/]/./g' > $target.dbn
        fi
        rm dssr-*
    done

    echo "generate _dp.ps.gz"
    cd $rootdir/cull/$resolu_cs/RNAfold
    for target in `cat $rootdir/cull/$resolu_cs/list`;do
        if [ -s "${target}_dp.ps.gz" ];then
            continue
        fi
	echo ">$target" > $target.fasta
	grep -PA1 "^>$target\b" $rootdir/cull/pdb_atom.sort.$resolu_cs|grep -v '>'|tr "[:lower:]" "[:upper:]" >> $target.fasta
        $bindir/RNAfold -p  $target.fasta
	gzip ${target}_dp.ps
	rm $target.fasta ${target}_ss.ps
    done

    echo "re-copy $PWD/*.pdb"
    cd $rootdir/cull/$resolu_cs/
    $bindir/clean_pdb.py -dir=$rootdir/pdb/data/structures/all/pdb/ list -suffix=.pdb.gz .pdb -StartIndex=1 -NewChainID=_
    cd $rootdir/cull/$resolu_cs/DSSR
    for target in `cat $rootdir/cull/$resolu_cs/list`;do
        if [ -s "$target.dbn" ];then
            continue
        fi
        $bindir/cssr $rootdir/cull/$resolu_cs/$target.pdb 1 | sed 's/[\[{}()]/./g'|sed 's/]/./g'  > $target.dbn
    done
done

echo "make smaller subsets"
resolu_cs_list="all_c1.0_s1.0 all_c0.9_s0.9 all_c0.8_s0.8 all_c0.8_s0.0"

for resolu_cs in $resolu_cs_list;do
    echo $resolu_cs
    mkdir -p $rootdir/cull/$resolu_cs
    grep '>' $rootdir/cull/pdb_atom.sort.$resolu_cs|cut -f1|sed 's/>//g' > $rootdir/cull/$resolu_cs/list
    $bindir/clstr2txt.py $rootdir/cull/pdb_atom.sort.${resolu_cs}.clstr $rootdir/cull/pdb_atom.sort.${resolu_cs}.txt
    cut -f1 $rootdir/cull/pdb_atom.sort.${resolu_cs} > $rootdir/cull/$resolu_cs/list.fasta
done

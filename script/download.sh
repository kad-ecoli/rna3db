#!/bin/bash
bindir=`dirname $(readlink -e $0)`
rootdir=`dirname $bindir`

echo "database path $rootdir"

echo "Download sequence to $rootdir/pdb/derived_data"
mkdir -p $rootdir/pdb/derived_data/index
wget -q ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz -O $rootdir/pdb/derived_data/pdb_seqres.txt.gz
$bindir/extractNA.py $rootdir/pdb/derived_data/pdb_seqres.txt.gz $rootdir/pdb/derived_data/na_seqres.txt
$bindir/getChainList.py $rootdir/pdb/derived_data/na_seqres.txt $rootdir/pdb/derived_data/na_chain.list
if [ -s "$rootdir/pdb/derived_data/na_type.list" ];then
    $bindir/getRNAlist.py $rootdir/pdb/derived_data/na_chain.list $rootdir/pdb/derived_data/na_type.list $rootdir/pdb/derived_data/na_chain.list
fi

echo "Download PDB entries"
mkdir -p $rootdir/pdb/data/structures/all/pdb/
cd       $rootdir/pdb/data/structures/all/pdb/
for chain in `cat $rootdir/pdb/derived_data/na_chain.list`;do
    if [ ! -s "${chain}.pdb.gz" ];then
        if [ ! -s "$chain.pdb" ];then
            $bindir/fetch.py $chain
        fi
        gzip ${chain}.pdb
    fi
done
rm `ls|grep -P "^[\da-z]{4}\.pdb"` `ls|grep pdb-bundle.tar.gz`
$bindir/getNAtype.py $rootdir/pdb/derived_data/na_chain.list $rootdir/pdb/derived_data/na_type.list
$bindir/getRNAlist.py $rootdir/pdb/derived_data/na_chain.list $rootdir/pdb/derived_data/na_type.list $rootdir/pdb/derived_data/na_chain.list
$bindir/removeDNAchain.py $rootdir/pdb/derived_data/na_type.list

echo "Download resolution"
wget -q ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/resolu.idx -O $rootdir/pdb/derived_data/index/resolu.idx

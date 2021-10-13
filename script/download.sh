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
    pdb=`echo $chain|cut -c1-4`
    if [ ! -s "$pdb" ];then
        mkdir $pdb
    fi
    if [ ! -s "${pdb}/${chain}.pdb.gz" ];then
        cd $rootdir/pdb/data/structures/all/pdb/$pdb
        $bindir/fetch.py $chain
	if [ ! -s "${chain}.pdb" ];then
	    if [ ! -s "${pdb}.cif" ];then
                $bindir/fetch.py -outfmt=cif $pdb
	    fi
	    $bindir/cif2pdb ${pdb}.cif ${chain}.pdb -chain `echo $chain|cut -c5-`
	fi
	if [ -s "${chain}.pdb" ];then
            gzip ${chain}.pdb
	fi
        cd $rootdir/pdb/data/structures/all/pdb/
    fi
done
rm `find -type f|grep -v ".pdb.gz$"`
$bindir/getNAtype.py $rootdir/pdb/derived_data/na_chain.list $rootdir/pdb/derived_data/na_type.list
$bindir/getRNAlist.py $rootdir/pdb/derived_data/na_chain.list $rootdir/pdb/derived_data/na_type.list $rootdir/pdb/derived_data/na_chain.list
$bindir/removeDNAchain.py $rootdir/pdb/derived_data/na_type.list

echo "Download resolution"
wget -q ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/resolu.idx -O $rootdir/pdb/derived_data/index/resolu.idx

echo "Download release date"
for pdb in `cut -c1-4 $rootdir/pdb/derived_data/na_chain.list|uniq`;do 
    echo $pdb `$bindir/get_pdb_release_date.sh $pdb`|sed 's/ /\t/g'
done > $rootdir/pdb/derived_data/na_chain.date 

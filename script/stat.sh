#!/bin/bash
bindir=`dirname $(readlink -e $0)`
rootdir=`dirname $bindir`
culldir="$rootdir/cull"

echo "L vs DSSR"
mkdir -p $rootdir/stat
cd $rootdir/stat
$bindir/statLvsDSSR.py $culldir/all_c1.0_s1.0/list.atomic $culldir/pdb_atom.sort.4.0_c0.8_s0.8 $culldir/all_c1.0_s1.0/DSSR/

echo "NaTorsion"
$bindir/batchNaTorsion.py $culldir/all_c1.0_s1.0/list.atomic $culldir/pdb_atom.sort.4.0_c0.8_s0.8 $culldir/all_c1.0_s1.0/
$bindir/statNaTorsion.py NaTorsion.raw.gz

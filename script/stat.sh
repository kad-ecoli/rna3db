#!/bin/bash
bindir=`dirname $(readlink -e $0)`
rootdir=`dirname $bindir`
culldir="$rootdir/cull"

mkdir -p $rootdir/stat
cd $rootdir/stat

echo "L vs DSSR"
$bindir/statLvsDSSR.py $culldir/all_c1.0_s1.0/list.atomic $culldir/pdb_atom.sort.4.0_c0.8_s0.8 $culldir/all_c1.0_s1.0/DSSR/

echo "NaTorsion"
$bindir/batchNaTorsion.py $culldir/all_c1.0_s1.0/list.atomic $culldir/pdb_atom.sort.4.0_c0.8_s0.8 $culldir/all_c1.0_s1.0/
$bindir/batchNaTorsion2.py $culldir/all_c1.0_s1.0/list.atomic $culldir/pdb_atom.sort.4.0_c0.8_s0.8 $culldir/all_c1.0_s1.0/
$bindir/statNaTorsion.py NaTorsion.raw.gz
$bindir/statNaTorsion2.py NaTorsion.raw.gz
$bindir/statNaTorsion3.py NaTorsion2.raw.gz

echo "BPtorsion"
$bindir/batchBPtorsion.py $culldir/all_c1.0_s1.0/list.atomic $culldir/pdb_atom.sort.4.0_c0.8_s0.8 $culldir/all_c1.0_s1.0/ $culldir/all_c1.0_s1.0/DSSR/
$bindir/batchBPtorsion2.py $culldir/all_c1.0_s1.0/list.atomic $culldir/pdb_atom.sort.4.0_c0.8_s0.8 $culldir/all_c1.0_s1.0/ $culldir/all_c1.0_s1.0/DSSR/
$bindir/statBPtorsion.py BPtorsion.raw.gz
$bindir/statBPtorsion2.py BPtorsion2.raw.gz

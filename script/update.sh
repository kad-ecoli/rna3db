#!/bin/bash
bindir=`dirname $(readlink -e $0)`
rootdir=`dirname $bindir`

echo $bindir/download.sh
$bindir/download.sh

echo $bindir/cull.sh
$bindir/cull.sh

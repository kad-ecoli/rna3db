#!/bin/bash
FILE=`readlink -e $0`
bindir=`dirname $FILE`
rootdir=`dirname $bindir`

echo "Download goa"
if [ ! -d "$rootdir/GO/goa/UNIPROT" ];then
    mkdir -p $rootdir/GO/goa/UNIPROT
fi

wget -q ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/README -O $rootdir/GO/goa/UNIPROT/README
curl -s ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz | zgrep -P '(^RNAcentral)|(^!)' > $rootdir/GO/goa/UNIPROT/goa_uniprot_all.gaf
gzip -f $rootdir/GO/goa/UNIPROT/goa_uniprot_all.gaf

echo "Download obo"
if [ ! -d "$rootdir/obo/go/extensions" ];then
    mkdir -p $rootdir/obo/go/extensions
fi
wget -q http://purl.obolibrary.org/obo/go/go-basic.obo -O $rootdir/obo/go/go-basic.obo
wget -q http://purl.obolibrary.org/obo/go.obo -O $rootdir/obo/go.obo
wget -q http://purl.obolibrary.org/obo/go.owl -O $rootdir/obo/go.owl
wget -q http://purl.obolibrary.org/obo/go/extensions/go-plus.owl -O $rootdir/obo/go/extensions/go-plus.owl

echo "curate goa"
if [ ! -d "$rootdir/stat" ];then
    mkdir -p "$rootdir/stat"
fi

$bindir/GOstat.py $rootdir/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz $rootdir/obo/go/go-basic.obo $rootdir/stat/GOstat.txt

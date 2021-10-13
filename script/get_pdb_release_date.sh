#!/bin/bash
if [ -z "$1" ];then
    echo get_pdb_release_date.sh pdbid
fi
PDBID="$1"
curl -s https://www.rcsb.org/structure/$PDBID  | grep -ohP "Released:\S+?\d{4}-\d{2}-\d{2}"  |grep -ohP "\d{4}-\d{2}-\d{2}"

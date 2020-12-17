#! /usr/bin/bash

# USAGE:
# % train.sh CONF,
# where CONF contains both .lst file and command line arguments
# for mxfold2'

IMG=/home/mah258/Pyle/mxfold2/mxfold2.simg
SEQS=/home/mah258/Pyle/data/pdb-cleaned/VL/VL-sequences
LABS=/home/mah258/Pyle/data/pdb-cleaned/VL/VL-labels
OUTPUT=output.csv
BEST=best.pt

dir=$(dirname $0)
out_conf=$(grep -A1 save-config $1 | tail -n 1)
epochs=$(grep -A1 epochs $1 | tail -n 1)
list=$(grep -A1 list $1 | tail -n 1)
param=$(grep -A1 param $1 | tail -n 1)
in_conf=$(sed -e 'N; /--list\n.*$/d' -e 'N; s/epochs\n.*$/epochs\n1/' $1)

echo "Running epoch 1..."
min_mcc=$(singularity exec $IMG mxfold2 train $list $in_conf; \
        python3 $dir/eval.py $IMG $SEQS $LABS $out_conf $OUTPUT; \
        python3 $dir/get_mcc.py $OUTPUT)

in_conf=

for epoch in $(seq 2 $(($epochs - 1))); do
    echo "Running epoch $epoch..."
    cur_mcc=$(singularity exec $IMG mxfold2 train $list $in_conf; \
        python3 $dir/eval.py $IMG $SEQS $LABS $out_conf $OUTPUT; \
        python3 $dir/get_mcc.py $OUTPUT)
    if [[ $cur_mcc < $min_mcc ]]; then
        min_mcc=$cur_mcc
        cp $param $BEST 
    fi
done

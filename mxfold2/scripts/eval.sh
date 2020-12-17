IMG=/home/mah258/Pyle/mxfold2/mxfold2.simg

singularity exec $IMG python3 $(dirname $0)/eval.py $1 $2 $3 $4

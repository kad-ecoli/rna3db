import os
import sys
import argparse
from glob import glob
import subprocess
from parse import *

def main():
    img, seqs, labs, conf, output = sys.argv[1:]
    labs = sorted(glob(labs+'/*'))
    seqs = sorted(glob(seqs+'/*'))
    f = open(output, 'w')
    f.write('id,f1,mcc\n')
    for lab, seq in zip(labs, seqs):
        name, _ = os.path.splitext(os.path.basename(seq))
        print(f'Predicting {name}...')
        commands = ['singularity', 'exec', img, 'mxfold2', 'predict', f'@{conf}', seq]
        output_bytes = subprocess.check_output(commands)
        output = output_bytes.decode('utf-8')
        seq, pred_lab, BPpred = parse_dot_bracket(output)
        ground_lab, BPpred = read_label(lab)
        F1, MCC = calcF1MCC(seq, ground_lab, pred_lab)
        f.write(f'{name},{F1},{MCC}\n')
    f.close()

if __name__ == '__main__':
    main()

import os
import pandas as pd
import sys
import numpy as np
import argparse 
import subprocess
from pyteomics.fasta import read

parser = argparse.ArgumentParser(description='Find necessary CDSs')
parser.add_argument('-orgs', 
    type=str, 
    help='Organisms payth', 
    required=True
    )
parser.add_argument('-save_way', 
    type=str, 
    help='output', 
    required=True
    )
args = parser.parse_args()

organisms = os.listdir(args.orgs)
save_way = args.save_way

for org in organisms:
    print('Analyzing of {} ...'.format(org))
    stains = os.listdir(args.orgs + org + '/All/')
    for stain in stains:
        if '.DS_Store' in stain:
            continue
        print('Analyzing of {} ...'.format(stain[: -12]))
        desc_stain = read(args.orgs + org + '/All/' + stain)
        description = ''
        for line in desc_stain:
            description += line.description.split()[0] + ', '
        description = description[: -2]
        comand = 'ncbi-acc-download --out {}{}/{}.gbk {}'.format(save_way, org, stain[: -12], description)
        subprocess.call(comand, shell=True)

print('Done!')    
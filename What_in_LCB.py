import argparse 
import pandas as pd 
import os
import sys
from os.path import split
from pyteomics.fasta import read

parser = argparse.ArgumentParser(description='What in SBs')
parser.add_argument('-sbs', 
    type=str, 
    help='SB', 
    required=True
    )
parser.add_argument('-antn', 
    type=str, 
    help='Annotation', 
    required=True
    )
parser.add_argument('-orgs', 
    type=str, 
    help='References of organisms', 
    default=None
    )
parser.add_argument('-save_way', 
    type=str, 
    help='Save way of annotation', 
    required=True
    )
args = parser.parse_args()

SBs = os.listdir(args.sbs)
coord_prot = os.listdir(args.antn)
orgsnisms = os.listdir(args.orgs)
part = ''
condition = ''
coordinate_of_split = ''
check_list = {}

#Maiking spike 
for org in orgsnisms:
    if '.DS_Store' in org:
        continue
    strains = os.listdir('{}{}/All/'.format(args.orgs, org))
    for strain in strains:
        if '.DS_Store' in strain:
            continue
        check_list[strain[: -12]] = {}
        fna = read('{}{}/All/{}'.format(args.orgs, org, strain))
        counter = 0
        index = 0
        for line in fna:
            counter += len(line.sequence)
            if index != 0:
                check_list[strain[: -12]][line.description.split()[0]] = counter
            else:
                check_list[strain[: -12]][line.description.split()[0]] = 0
            index += 1

for organism in SBs:
    if '.DS_store' in organism:
        continue
    for an_org in coord_prot:
        if '.DS_Store' in an_org:
            continue
        if '.txt' not in an_org:
            continue
        if an_org[11: -4] in organism:
            print('Analyzing of {} ...'.format(an_org[11: -4]))
            save_file = open(args.save_way + '{}.txt'.format(an_org[11: -4]), 'w')
            save_file.write('Strain\tSB\tCoord_of_SB\tProtein\tCoord_of_prot\tPart\tCondition\n')
            SB_table = pd.read_csv(args.sbs + organism, 
                header=None,
                names=('SB', 'Strain', 'Coordinates'), 
                sep='\t'
                )
            coord = pd.read_csv(args.antn + an_org, 
                header=None, 
                names=('Strain', 'Contig', 'Coordinates', 'Protein', 'Part'), 
                sep='\t'
                )
            for SB_ind in range(len(SB_table)):
                for ind in range(len(coord)):
                    if SB_table.Strain[SB_ind] == coord.Strain[ind]:
                        for strain in check_list:
                            if strain == SB_table.Strain[SB_ind]:
                                for key in check_list[strain]:
                                    if coord.Contig[ind] in key:
                                        coef = check_list[strain][key]
                                        #print('Analyzing of {} ...'.format(line.description.split()[1]))
                                        start_SB = int(SB_table.Coordinates[SB_ind].split(':')[0][1:])
                                        end_SB = int(SB_table.Coordinates[SB_ind].split(':')[1][: -1])
                                        start_prot = int(coord.Coordinates[ind].split(',')[0][1:])
                                        end_prot = int(coord.Coordinates[ind].split(',')[1][: -1])
                                        if coord.Part[ind] == 'partI':
                                            part = 'partI'
                                        elif coord.Part[ind] == 'partII':
                                            part = 'partII'
                                        else:
                                            part = 'full'
                                        if start_SB <  start_prot and end_SB > end_prot:
                                            condition = 'full protein'
                                            save_file.write('{Strain}\t{SB}\t{Start_SB}:{End_SB}\t{Protein}\t{Start_Prot}:{End_Prot}\t{Part}\t{Condition}\n'.format(
                                                Strain =  SB_table.Strain[SB_ind],
                                                SB =  SB_table.SB[SB_ind],
                                                Start_SB = start_SB,
                                                End_SB = end_SB,
                                                Protein = coord.Protein[ind],
                                                Start_Prot = start_prot,
                                                End_Prot = end_prot,
                                                Part = part,
                                                Condition = condition
                                                ))

                                        elif start_SB > start_prot and  start_SB < end_prot and end_SB > end_prot:
                                            condition = 'front end'
                                            save_file.write('{Strain}\t{SB}\t{Start_SB}:{End_SB}\t{Protein}\t{Start_Prot}:{End_Prot}\t{Part}\t{Condition}\n'.format(
                                                Strain =  SB_table.Strain[SB_ind],
                                                SB =  SB_table.SB[SB_ind],
                                                Start_SB = start_SB,
                                                End_SB = end_SB,
                                                Protein = coord.Protein[ind],
                                                Start_Prot = start_prot,
                                                End_Prot = end_prot,
                                                Part = part,
                                                Condition = condition
                                                ))

                                        elif start_SB < start_prot and end_SB < end_prot and start_prot < end_SB:
                                            condition = 'back end'
                                            save_file.write('{Strain}\t{SB}\t{Start_SB}:{End_SB}\t{Protein}\t{Start_Prot}:{End_Prot}\t{Part}\t{Condition}\n'.format(
                                                Strain =  SB_table.Strain[SB_ind],
                                                SB =  SB_table.SB[SB_ind],
                                                Start_SB = start_SB,
                                                End_SB = end_SB,
                                                Protein = coord.Protein[ind],
                                                Start_Prot = start_prot,
                                                End_Prot = end_prot,
                                                Part = part,
                                                Condition = condition
                                                ))
                        
            save_file.close()
            print('Analizing of {} was done!'.format(an_org[11: -4]))
            
print('Done!')
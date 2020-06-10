import argparse
import sys
import os
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description='Parse CDS')
parser.add_argument('-gbks', help='Organisms path', required=True)
parser.add_argument('-save_way', help='Output derectory', required=True)
args = parser.parse_args()

organisms = os.listdir(args.gbks)

for org in organisms:
    if '.DS_Store' in org:
        continue
    print('Analyze of {}'.format(org))
    savefile = open('{}Annotation_{}.txt'.format(args.save_way, org), 'w')
    strains = os.listdir(args.gbks + org)
    for strain in strains:
        if '.DS_Store' in strain:
            continue
        print('Analyze of {}'.format(strain[: -4]))
        for record in SeqIO.parse(args.gbks + org + '/{}'.format(strain), "genbank"):
            for i in record.features:
                if i.qualifiers.get('product') == None:
                    continue
                if i.qualifiers.get('translation') == None:
                    continue
                coord = str(i.location).replace('>', '')
                coord = coord.replace('<', '')
                if 'join' in coord:
                    savefile.write('{}\t{}\t{}\t{}\tpartI\n'.format(
                        strain[: -4], record.name, coord[5: -1].split(',')[0].replace(':', ',')[: -3].replace(' ', ''), i.qualifiers.get('product')[0].replace(' ', '_')))
                    savefile.write('{}\t{}\t{}\t{}\tpartII\n'.format(
                        strain[: -4], record.name, coord[5: -1].split(',')[1].replace(':', ',')[: -3].replace(' ', ''), i.qualifiers.get('product')[0].replace(' ', '_')))
                else:
                    savefile.write('{}\t{}\t{}\t{}\tfull\n'.format(
                        strain[: -4], record.name, coord.replace(':', ',')[: -3].replace(' ', ''), i.qualifiers.get('product')[0].replace(' ', '_')))
                #savefile.write(i.qualifiers.get('translation')[0] + '\n')

    savefile.close()
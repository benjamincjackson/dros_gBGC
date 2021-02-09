"""
Take Henry/Padraig's custom bed format file for representing multiple alignments
and write new versions of the reference sequence that the file is based on, with
alleles from the all of the aligned species.
"""

import pysam, argparse, csv, sys, gzip, copy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

parser = argparse.ArgumentParser(description="""convert wgabed alignment to fasta files""")

parser.add_argument('reference_fai_file', help='the index for the reference file used to align to, in fasta format - to get the chromosome lengths')
parser.add_argument('wga_bed_file', help='Henry and Padraig\'s custom bed file representation of multple alignments')
args = parser.parse_args()

tbx = pysam.TabixFile(args.wga_bed_file)

# initiate an empty dict:
temp_dict = {}

# read in the index file and fill the dictionary above with the contig lengths
with open(args.reference_fai_file, 'r') as f:
    textreader = csv.reader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
    for line in textreader:
        temp_dict[line[0]] = line[1]

# template_reference_dict = {x: int(temp_dict[x]) * ['.'] for x in temp_dict.keys()}

# then make a dictionary of dictionaries contig sequences,
# by first, getting the list of the species in the multiple alignment from the
# first line of the bed file
with gzip.open(args.wga_bed_file, mode = 'rt') as f:
    textreader = csv.reader(f, delimiter='\t')
    row = next(textreader)
    species_names = row[4].split(',')

# and then make some a dictionary of dictionaries - one for each species, to
# populate with the alleles from the multiple alignment:
reference_dicts = {x: {y: int(temp_dict[y]) * ['.'] for y in temp_dict.keys()} for x in species_names}

for chrom in temp_dict.keys():
# for chrom in ['2L']:
    if tbx.fetch(chrom):
        for row in tbx.fetch(chrom, parser=pysam.asBed()):
            # print(row)
            contig = row[0]
            start = int(row[1])
            stop = int(row[2])

            # make a dictionary of the alleles for the species at this position
            allele_dict_raw = {row[4].split(',')[x]: row[7].split(',')[x] for x in range(len(species_names))}

            # check if there are any gaps in the sim alleles, in which case we want
            # to disregard the equivalent alleles in the mel and yak alignments
            ref_alleles = list(allele_dict_raw[species_names[0]])

            if '-' in ref_alleles:

                indx = [x for x,y in enumerate(ref_alleles) if y != '-']

                # as a check, the length of the index should equal the difference
                # bewteen the start and stop coordinates
                if len(indx) is not (stop - start):
                    print('something went wrong at pos ' + contig + ':' + str(start) + '-' + str(stop))

                # fill a new dictionary with the appropriate length alleles
                for species in species_names:
                    oldalleles = allele_dict_raw[species]
                    newalleles = ''
                    for x in indx:
                        newalleles = newalleles + oldalleles[x]

                    reference_dicts[species][contig][start:stop] = copy.copy(newalleles)

            else:
                for species in species_names:
                    reference_dicts[species][contig][start:stop] = allele_dict_raw[species]

            # print('after the loop: ' + str(reference_dicts['dsim'][contig][start:stop]))
            # print('after the loop: ' + str(reference_dicts['dmel'][contig][start:stop]))
            # print('after the loop: ' + str(reference_dicts['dyak'][contig][start:stop]))
            # print()

# NB change this so that the code is independent of manually entered species names
if reference_dicts['dsim'] == reference_dicts['dmel']:
    if reference_dicts['dsim'] == reference_dicts['dyak']:
        print('oh fuck')

# write some fasta files
for species in species_names:
        filename = species_names[0] + '_' + species + '_alleles.fasta'
        records = [SeqRecord(Seq(''.join(ref_seq), generic_dna), id = contig, description = '') for contig, ref_seq in reference_dicts[species].items()]
        SeqIO.write(records, filename, "fasta")














#

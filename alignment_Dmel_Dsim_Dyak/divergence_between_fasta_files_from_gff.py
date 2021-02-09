# check divergence between mouse and rat reference sequences.

import gzip, argparse, sys, csv, Bio, glob
from pyfaidx import Fasta

parser = argparse.ArgumentParser(description="""check divergence between two fasta files using
                                                every line in a gff file""")

parser.add_argument('input_fasta_file_1', help='first fasta file')
parser.add_argument('input_fasta_file_2', help='second (aligned) fasta file')
parser.add_argument('gff_file', help='a gzipped gff-format file')

args = parser.parse_args()


# get at the reference sequences with pyfaidx
fasta1 = Fasta(args.input_fasta_file_1)
fasta2 = Fasta(args.input_fasta_file_2)

# for key in reference_dict_2:
#     print(reference_dict_2[key])
#
# for key in reference_dict_2['seq2']:
#     print(key)
#
# print(reference_dict_2['seq2'][4])

# initiate a count of differences
diff=0
# initiate a count of total number of tests
denom=0

with gzip.open(args.gff_file, 'rt') as f:
    reader = csv.reader(f, delimiter = '\t')
    for line in reader:
        chrom = line[0]
        start = int(line[3])
        end = int(line[4])

        for pos in range(start, end + 1):

            if chrom not in fasta1 or chrom not in fasta2:
                continue
                
            allele1 = fasta1[chrom][pos - 1].seq.upper()
            allele2 = fasta2[chrom][pos - 1].seq.upper()
            #allele1 = 'N'
            #allele2 = '.'
            #print(line[0] + ' ' + line[1] + ' ' + line[9] + ' first allele is ' + allele1 + ' second allele is ' + allele2)

            # check that both alleles are nucleotides
            LUT = ['A', 'T', 'C', 'G']
            # if allele1 in LUT:
            #     print('it\'s a nucleotide')
            if allele1 in LUT and allele2 in LUT:
                #print('it\'s a nucleotide')
                denom+=1

                if allele1 != allele2:
                    diff+=1


print('differences = ' + str(diff) +', total comparisons = ' + str(denom) + ', proportion = ' + str(diff/denom))

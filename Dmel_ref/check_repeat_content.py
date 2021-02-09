"""
Read in a soft-masked fasta file and get the GC/AT content of upper case sites
(not repeats) and lower case sites (soft masked as repeats).
"""


from pyfaidx import Fasta
import sys, argparse

parser = argparse.ArgumentParser(description="""report GC content stats on a soft-masked fasta file""")
parser.add_argument('fasta_file', help='the fasta file to analyse')
parser.add_argument('--contigs',
                    nargs = '+',
                    help='an optional list of contigs you want to restrict the analysis to')

args = parser.parse_args()

fasta = Fasta(args.fasta_file)

if args.contigs:
    contigs = [x for x in args.contigs]
else:
    contigs = fasta.keys()

lookup = {'A': 0, 'C': 0, 'G': 0, 'T': 0,
          'a': 0, 'c': 0, 'g': 0, 't': 0,
          'N': 0, 'n': 0}

for contig in contigs:
    for x in range(len(fasta[contig])):
        nucleotide = str(fasta[contig][x])
        lookup[nucleotide]+=1


GC_unmasked = lookup['G'] + lookup['C']
AT_unmasked = lookup['A'] + lookup['T']

GC_masked = lookup['g'] + lookup['c']
AT_masked = lookup['a'] + lookup['t']

print('GC unmasked count is : ' + str(GC_unmasked))
print('AT unmasked count is : ' + str(AT_unmasked))
print('GC masked count is : ' + str(GC_masked))
print('AT masked count is : ' + str(AT_masked))

print('GC content of unmasked sequence is: ' + str(GC_unmasked / (GC_unmasked + AT_unmasked)))
print('GC content of masked sequence is: ' + str(GC_masked / (GC_masked + AT_masked)))
print('total GC content is: ' + str((GC_unmasked + GC_masked) / (GC_unmasked + AT_unmasked + GC_masked + AT_masked)))

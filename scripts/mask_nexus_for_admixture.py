"""
Use the coordinates of admixture tracts provided by the Dmel nexus (Lack et al
2015) to rewrite individual fly line's fasta files with admixed regions masked
by Ns
"""

import sys, argparse, glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

parser = argparse.ArgumentParser(description="""mask admixed regions of the D mel Nexus data""")
parser.add_argument('input_directory', help='the folder with the fasta files to mask in it, include trailing \'/\'')
parser.add_argument('output_directory', help='the folder to write the masked data to, include trailing \'/\'')
args = parser.parse_args()

# arguments specified on the command line:
indir = args.input_directory
outdir = args.output_directory

# here's the lengths of the chromosomes for Dmel v5.57
Dmel_ref = {'2L': 23011544,
            'X': 22422827,
            '3L': 24543557,
            '4': 1351857,
            '2R': 21146708,
            '3R': 27905053}

# read in the regions to be masked, and store in a dictionary
# initiate the dictionary:
counter=0
tracts = {}
with open('/Users/ben/data/dros_X_A_mut/Dmel_nexus/dpgp3_admixture_tracts.txt', 'rt') as f:
    for line in f:
        l = line.strip().split()
        sample = l[0]
        chr = l[1].strip('Chr')
        strt = l[2]
        stp = l[3]

        # if the sample isn't in the dictionary yet:
        if not sample in tracts:
            # initiate a dictionary of chromosome arms for the sample
            tracts[sample] = {'2L': [],
                              '2R': [],
                              '3R': [],
                              '3L': [],
                              '4':[],
                              'X': []}

        # add the start and stop coordinates of the admixture mask to the list for
        # this chromosome arm for this sample:
        tracts[sample][chr].append(strt + ':' + stp)

# # write a file with the proportion of genome masked for admixture:
# total_length = sum([int(i) for i in Dmel_ref.values()])
# with open('/Users/ben/data/dros_X_A_mut/Dmel_nexus/proportion_admixture_masked.txt', 'wt') as f:
#     f.write('sample\tproportion_masked\n')
#     for sample in tracts:
#         if sample[0:2] == 'ZI' or sample[0:2] == 'RG':
#             chunks = []
#             for chr in tracts[sample]:
#                 if len(tracts[sample][chr]) > 0:
#                     for pair in tracts[sample][chr]:
#                         strt = int(pair.split(':')[0]) - 1
#                         stp = int(pair.split(':')[1])
#                         gapsize = stp - strt
#                         chunks.append(gapsize)
#             f.write(sample + '\t' + str(round(sum(chunks) / total_length, 5)) + '\n')

for FASTA in glob.glob(indir + '*.fa'):
    fasta_dict = SeqIO.to_dict(SeqIO.parse(FASTA, "fasta"))

    newfilename = FASTA.split('/')[-1]
    sample = newfilename.split('.')[0]

    # print(FASTA)
    # print(newfilename)
    # print(sample)
    # print(fasta_dict)
    # print()

    if sample in tracts:
        for chr in tracts[sample]:
            if len(tracts[sample][chr]) > 0:
                for pair in tracts[sample][chr]:
                    strt = int(pair.split(':')[0]) - 1
                    stp = int(pair.split(':')[1])
                    gapsize = stp - strt

                    tempseq = str(fasta_dict[chr].seq)

                    newrecord = SeqRecord(Seq(tempseq[0:strt] + 'N' * gapsize + tempseq[stp:]), id = chr, name = chr, description = chr)

                    fasta_dict[chr] = newrecord

                    # print(len(tempseq))
                    # print(len(newseq))
                    # print()
                    # print(tracts[sample])
                    # print(chr)
                    # print(strt)
                    # print(stp)
        # print(fasta_dict)

    with open(outdir + newfilename, 'wt') as f_out:
        for i in fasta_dict.keys():
            SeqIO.write(fasta_dict[i], f_out, "fasta")


        # counter+=1
        # if counter > 0:
        #     sys.exit()
























#

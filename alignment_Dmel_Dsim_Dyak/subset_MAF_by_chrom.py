"""
Split the complete Dmel Dsim Dyak MAF file by Dmel chromosome,
(so that I can index it later)
"""

import sys, gzip
from Bio import AlignIO
from Bio.AlignIO import MafIO

counter=0

twoL = []
threeL = []
twoR = []
threeR = []
header = []
# stop = False

with gzip.open('dmel.dsim.dyak.maf.gz', 'rt') as maf:
    for line in maf:
        if line[0] == '#':
            header.append(line.strip())
        else:
            break

    for block in AlignIO.parse(maf, 'maf'):
        # if stop:
        #     break
        for seqrec in block:
            if seqrec.id.startswith('dmel'):
                # print(block)
                chr = seqrec.id.split('.')[1]
                if chr == '2L':
                    twoL.append(block)
                if chr == '2R':
                    twoR.append(block)
                if chr == '3L':
                    threeL.append(block)
                if chr == '3R':
                    threeR.append(block)

                # counter +=1
                # if len(twoL) > 0:
                #     print(twoL[0])
                #     stop = True

                continue

with open('dmel.dsim.dyak.2L.maf', 'wt') as f_2L:
    f_2L.write('\n'.join(header))
    AlignIO.write(twoL, f_2L, 'maf')

with open('dmel.dsim.dyak.2R.maf', 'wt') as f_2R:
    f_2R.write('\n'.join(header))
    AlignIO.write(twoR, f_2R, 'maf')

with open('dmel.dsim.dyak.3L.maf', 'wt') as f_3L:
    f_3L.write('\n'.join(header))
    AlignIO.write(threeL, f_3L, 'maf')

with open('dmel.dsim.dyak.3R.maf', 'wt') as f_3R:
    f_3R.write('\n'.join(header))
    AlignIO.write(threeR, f_3R, 'maf')

idx_2L = MafIO.MafIndex("dmel.dsim.dyak.2L.maf.idx", 'dmel.dsim.dyak.2L.maf', "dmel.2L")
idx_2R = MafIO.MafIndex("dmel.dsim.dyak.2R.maf.idx", 'dmel.dsim.dyak.2R.maf', "dmel.2R")
idx_3L = MafIO.MafIndex("dmel.dsim.dyak.3L.maf.idx", 'dmel.dsim.dyak.3L.maf', "dmel.3L")
idx_3R = MafIO.MafIndex("dmel.dsim.dyak.3R.maf.idx", 'dmel.dsim.dyak.3R.maf', "dmel.3R")














#

"""
Get the alignment score for each Dmel short intron from the wgabed file
"""

import sys, pysam, statistics

isheader=0
counter=0

mydict = {}

wgabedtbx = pysam.TabixFile("../alignment_Dmel_Dsim_Dyak/dmel.dsim.dyak.wga.bed.gz")

# then get a score per intron
with open('../data/table.txt', 'rt') as f:
    for line in f:
        l = line.strip().split()
        if isheader == 0:
            isheader+=1
            continue
        mel_contig = l[0]

        if mel_contig == "X":
            continue

        mel_site = int(l[1])
        mel_name = l[119]

        score = float([x for x in wgabedtbx.fetch(mel_contig, mel_site - 1, mel_site)][0].split()[9])

        if not mel_name in mydict:
            mydict[mel_name] = [score]
        else:
            mydict[mel_name].append(score)

        # counter+=1
        # if counter > 1000:
        #     break

with open('../data/scores.txt', 'wt') as f:
    for intron in mydict:
        f.write(intron + '\t' + str(statistics.mean(mydict[intron])) + '\n')















#

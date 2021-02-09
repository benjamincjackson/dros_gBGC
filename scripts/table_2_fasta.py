"""
write per-intron fasta files from the tabular data I produced
"""

import sys
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

sitecounter=0
isheader=0
isfirst=True
with open('../data/table.txt', 'rt') as f:
    for line in f:
        l = line.strip().split()
        if isheader == 0:
            isheader+=1
            RGnames = l[8:29]
            ZInames = l[29:98]
            MDnames = l[98:119]
            continue

        melcontig = l[0]
        simcontig = l[2]

        mel = l[5]
        sim = l[6]
        yak = l[7]
        RG = l[8:29]
        ZI = l[29:98]
        MD = l[98:119]

        # if any([x not in ['A', 'C', 'G', 'T'] for x in [mel] + [sim] + [yak] + RG + ZI + MD]):
        #     continue

        melname = l[119]
        simname = l[120]
        melGC = l[121]
        simGC = l[122]

        if isfirst:
            filenameold = melcontig + '&' + simcontig + '&' + melname + '&' + simname + '&' + melGC + '&' + simGC

            melstring = mel
            simstring = sim
            yakstring = yak
            RGstrings = [i for i in RG]
            ZIstrings = [i for i in ZI]
            MDstrings = [i for i in MD]

            isfirst = False
            continue

        else:

            filenamenew = melcontig + '&' + simcontig + '&' + melname + '&' + simname + '&' + melGC + '&' + simGC

            if filenamenew == filenameold:
                melstring = melstring + mel
                simstring = simstring + sim
                yakstring = yakstring + yak
                RGstrings = [x + y for x,y in zip(RGstrings,RG)]
                ZIstrings = [x + y for x,y in zip(ZIstrings,ZI)]
                MDstrings = [x + y for x,y in zip(MDstrings,MD)]
                continue

            else:
                alignment = MultipleSeqAlignment([SeqRecord(Seq(melstring, generic_dna), id='dmel', description='')] +
                                                 [SeqRecord(Seq(simstring, generic_dna), id='dsim', description='')] +
                                                 [SeqRecord(Seq(yakstring, generic_dna), id='dyak', description='')] +
                                                 [SeqRecord(Seq(x, generic_dna), id=y, description='') for x,y in zip(RGstrings,RGnames)] +
                                                 [SeqRecord(Seq(x, generic_dna), id=y, description='') for x,y in zip(ZIstrings,ZInames)] +
                                                 [SeqRecord(Seq(x, generic_dna), id=y, description='') for x,y in zip(MDstrings,MDnames)])


                AlignIO.write(alignment, '/Users/ben/data/dros_X_A_mut/fasta_SIs/intersection/' + filenameold + '.fa', "fasta")
                sitecounter=sitecounter+len(melstring)
                filenameold =  melcontig + '&' + simcontig + '&' + melname + '&' + simname + '&' + melGC + '&' + simGC

                melstring = mel
                simstring = sim
                yakstring = yak
                RGstrings = [i for i in RG]
                ZIstrings = [i for i in ZI]
                MDstrings = [i for i in MD]

print(sitecounter)

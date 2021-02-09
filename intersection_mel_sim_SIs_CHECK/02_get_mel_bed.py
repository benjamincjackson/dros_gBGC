import sys, glob, gzip

with gzip.open('mel_SI.bed.gz', 'wt') as f_out:
    for file in glob.glob('/Users/ben/data/dros_X_A_mut/fasta_SIs/ZI/*.pos'):
        chr = file.split('/')[-1].split('.')[0]
        with open(file, 'rt') as f:
            positions = []
            for line in f:
                positions.append(int(line.strip().split()[0]))

            length = len(positions)

            start = positions[0]
            check = positions[0]
            counter = 0
            NEWSTART = False
            for pos in positions:
                if counter == 0:
                    counter+=1
                    continue
                if counter < length:
                    if pos - 1 == check:
                        check = pos
                        counter+=1
                    else:
                        f_out.write(chr + '\t' + str(start - 1) + '\t' + str(check) + '\n')
                        start = pos
                        check = pos
                        counter+=1

                if counter == length:
                    f_out.write(chr + '\t' + str(start - 1) + '\t' + str(pos) + '\n')

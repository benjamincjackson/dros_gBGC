import gzip

with gzip.open('../Dsim_annotation/SIs/SIs.sorted.uniq.gff.gz', 'rt') as f:
    for l in f:
        line = l.strip().split()
        gffchr, gffstart, gffstop = line[0], int(line[3]), int(line[4])
        bedstart = gffstart - 1
        bedstop = gffstop
        print(gffchr + '\t' + str(bedstart) + '\t' + str(bedstop) + '\t' + line[8])

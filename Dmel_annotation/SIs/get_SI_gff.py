"""
Get short introns from the header lines of a precomputed flybase fasta file of
introns
"""

import gffutils, gzip, csv, sys, random, re, pysam

# the gff info:
gff = pysam.TabixFile("../gff/dmel-all-r5.57.uniq.gff.cleaned.gz")

# a dict to look up the non-recombining regions according to Campos et al (2012),
# supplementary table 1:
NRdict = {'X': ['1:380103', '21947405:2277501'],
          '2L': ['21946551:2301154'],
          '2R': ['1:128568'],
          '3L': ['22754668:2454355'],
          '3R': ['-224037:378656']}


i=0
j=0
with gzip.open("../flybase_other/dmel-all-intron-r5.57.fasta.gz", 'rt') as f_in:
    for l in f_in:
        # only read header lines
        if l[0] is not '>':
            continue
        i+=1
        line = l.strip('\n')

        # the name of the intron (which includes FBtr and FBgn info) is the first item:
        name = line.split(' ', 1)[0]

        # followed by all the other info
        items = line.strip('; ').split(' ', 1)[1].split('; ')

        # can make a dict of all the info in the header line
        stuff = {}
        # and then populate it
        for item in items:
            k = item.split('=')[0]
            v = item.split('=')[1]
            stuff[k] = v

        # the SI- relevant bit is the length key:
        if int(stuff['length']) > 65:
            continue

        # here's the contig:
        contig = stuff['loc'].split(':')[0]

        # you can also tell whether the intron is on the forward or reverse strand
        # by the presence of "complement" next to the position information (complement = reverse strand).
        if not re.search('complement', stuff['loc']):
            # here are the coordinates of the intron:
            coordlist = stuff['loc'].split(':')[1].split('..')
            strand = '+'

        if re.search('complement', stuff['loc']):
            coordlist = stuff['loc'].split(':')[1].strip('complement(').strip(')').split('..')
            strand = '-'

        # print(coords)
        # and now test whether any of them lie in the NR region(s) for the appropriate contig:
        skip = False
        if not contig in NRdict:
            continue

        for coord in coordlist:
            for pair in NRdict[contig]:
                start = int(pair.split(':')[0])
                end = int(pair.split(':')[1])
                if int(coord) >= start and int(coord) <= end:
                    skip = True
        if skip:
            continue

        # look up the start/end coordinates in the gff file for overlapping exons/introns etc:

        # if the intron is on the positive strand, then take the start coordinate + 8--30
        if strand == '+':
            x0 = int(coordlist[0]) + 7
            x1 = int(coordlist[0]) + 29
        # if the intron is on the negative strand, then take the end coordinate - 8--30
        if strand == '-':
            x1 = int(coordlist[1]) - 7
            x0 = int(coordlist[1]) - 29

        safeintron = True
        # some weird errors arise without this patch:
        try:
            gff.fetch(contig, x0, x1)
        except UnicodeDecodeError:
            continue
        else:
            # get all the features that overlap 8-30bp of introns <66bp
            for row in gff.fetch(contig, x0, x1, parser=pysam.asTuple()):
                try:
                    info = [x for x in row]
                except UnicodeDecodeError:
                    continue
                else:
                    # we can explicitly check for two overlapping features: introns and exons.
                    # if an exon overlaps with this SI, then let's skip this SI:
                    info = [x for x in row]
                    if info[2] == "exon":
                        safeintron = False
                        break

                    # for introns it's more complicated. E.g. maybe there are multiple short introns at exactly the same position,
                    # i.e. from different transcripts?
                    if info[2] == 'intron':
                        # There will always be a match to the original short intron record.
                        # if the start, stop and strand features all match, then we can ignore it:
                        if coordlist[0] == info[3] and coordlist[1] == info[4] and strand == info[6]:
                            continue

                        # If the intron is long, then we want to skip the focal SI:
                        if (int(info[4]) - int(info[3])) > 65:
                            safeintron = False
                            break
                        else:
                            # do something where there is a short intron but its start/end coordinates don't exactly match
                            # the focal short intron
                            if (int(info[4]) - int(info[3])) < 66:
                                focal_strand = strand
                                test_strand = info[6]
                                # if the introns are on the +ve strand and start in the same place, this is ok:
                                if focal_strand == test_strand and focal_strand == '+':
                                    if int(info[3]) == int(coordlist[0]):
                                        continue
                                    # if they don't start in the same place, it's not ok
                                    else:
                                        safeintron = False
                                        break
                                # same as above, for the -ve strand
                                elif focal_strand == test_strand and focal_strand == '-':
                                    if int(info[4]) == int(coordlist[1]):
                                        continue
                                    else:
                                        safeintron = False
                                        break
                                # if theyre on different strands, it's not ok.
                                else:
                                    safeintron = False
                                    break

            # Print the fasta header line if safeintron = True?
            if safeintron == True:

                # Write a gff file line with SI coordinates?
                print(contig + "\tBen\tshort_intron_8-30\t" + str(x0) + "\t" + str(x1) + "\t.\t" + strand + "\t.\tparent=" + stuff['parent'] +  ";name=" + name.strip('>'))

















#

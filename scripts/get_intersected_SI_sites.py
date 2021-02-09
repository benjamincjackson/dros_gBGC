"""
I think that I want to make a table of sites that are SIs in both Dmel and Dsim,
one row per site, columns for each outgroup (allele), columns for the polymorphism
data, columns for the name of the intron in Dmel, and in Dsim, and GC content
of the column in Dmel and Dsim?

1) read overlap.bed line by line
2) to start with, just deal with rows that have length=1 in the wgabed info
3) use the info on the row, which is a site in dmel coords which is both a Dmel
   and a Dsim short intron to:
        * get the Dmel, Dsim and Dyak alleles
        * get the Dmel polymorphism data
        * get the Dmel intron name (think I can look this up from the original indexed
                                    bed file of short introns?)
        * ??? get the Dmel intron GC content ???
        * get the Dsim coordinate, and use this to:
            * get the Dsim polymorphism data
            * get the Dsim intron name (think I can look this up from the original indexed
                                        bed file of short introns?)
            * ??? get the Dsim intron GC content ???


"""

import argparse, sys, gzip, glob, os, pysam
from pyfaidx import Fasta
from pysam import VariantFile
sys.path.insert(0, '/Users/ben/programs/WGAbed/')
from get_out_of_bed import *
from Bio.SeqUtils import GC
from numpy.random import choice

# A function to extract a portion of a faidxed fasta file when there's a
# possibility that there's no sequence at those coordinates in the fasta file
# (because some of the dros nexus data have missing chromosomes)
def get_from_fasta(sample, contig, start, end):
    try:
        x = mel_fastas[sample.split('/')[-1].strip('.fa')][contig][start - 1: end].seq

    except (ValueError, KeyError):
        return 'N' * (end - (start - 1))

    return x

# a function to return the (haploid) genotype from a diploid tuple for these dros
# inbred lines. If there is residual heterozygosity, then sample an allele randomly
# according to depth (as in Jackson et al 2017)
def get_haploid_genotype(sample_info, site_alleles, strand = '+'):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    # here's the diploid genotype for one individual at one site:
    GT = sample_info['GT']
	# if there's more than one allele here, we want to sample one of them
	# proportional to the depth of coverage of the two
    if len(set(GT)) > 1:
        AD = sample_info['AD']
        draw = choice(GT, 1, AD)[0]
        al = site_alleles[draw]
    else:
        al = site_alleles[GT[0]]

    if strand == '+':
        return(al)
    if strand == '-':
        return(complement[al])

# a function to take a chunk of VCF file and turn it into a dict with the format
# sample_name: sequence (in this case, applying some filters and converting diploid
# base calls to haploid genotypes for these Dros inbred lines). Also need to give it
# the sample size if you want to filter out sites with missing data
def vcfchunk_to_dict(n, variantfile, contig, start, end, strand):
    samples = variantfile.header.samples
    records = variantfile.fetch(contig, start - 1, end)
    seqdict = {sample: '' for sample in samples}

    for record in records:
        # need to apply some filters at this stage:
		# the "n != record.info['AN']" clause means that sites with any missing data will be skipped
        if not 'AN' in record.info or n != record.info['AN'] or len(record.ref) > 1 or record.qual < 30:
			## Skip sites that don't pass the tests entirely by commenting out the
			## two lines below (don't add Ns)
            # for sample in record.samples:
            #     seqdict[sample] = seqdict[sample] + 'N'
            continue

        if record.alts:
            if record.qual < 30 or len(max(record.alts, key = len)) > 1:
				## Skip sites that don't pass the tests entirely by commenting out the
				## two lines below (don't add Ns)
                # for sample in record.samples:
                #     seqdict[sample] = seqdict[sample] + 'N'
                continue

        for sample in record.samples:
            seqdict[sample] = seqdict[sample] + get_haploid_genotype(record.samples[sample],
                                                                     record.alleles,
                                                                     strand = strand)

    return(seqdict)

# A dictionary of non-recombining regions in Dmel, to exclude them from
# downstream analysis
NRdict = {'X': ['1:380103', '21947405:2277501'],
          '2L': ['21946551:2301154'],
          '2R': ['1:128568'],
          '3L': ['22754668:2454355'],
          '3R': ['-224037:378656']}

NRdictsim = {'Scf_X': ['1:282985', '20452015:20829647'],
            'Scf_2L': ['21188738:23539531'],
            'Scf_2R': ['1:2053631'],
            'Scf_3L': ['22157209:24153973'],
            'Scf_3R': []}

# a dict of where to find each Dmel population's alleles
pop_dict = {'ZI': '../Dmel_nexus/dpgp3_sequences/fasta_masked/',
			'RG': '../dpgp2_sequences/fasta_masked/'}

# lists of alleles to retain after downsampling
RG21 = ['RG34', 'RG19', 'RG7', 'RG24', 'RG36', 'RG37N', 'RG39', 'RG4N', 'RG2',
        'RG32N', 'RG9', 'RG13N', 'RG18N', 'RG22', 'RG25', 'RG28', 'RG3', 'RG33',
        'RG38N', 'RG5', 'RG6N']
ZI69 =  ['ZI114N', 'ZI117',  'ZI161',  'ZI181',  'ZI184',  'ZI194',  'ZI207',
         'ZI210',  'ZI211',  'ZI213',  'ZI214',  'ZI219',  'ZI232',  'ZI233',
         'ZI235',  'ZI239',  'ZI250',  'ZI252',  'ZI253',  'ZI254N', 'ZI255',
         'ZI264',  'ZI265',  'ZI267', 'ZI268',  'ZI27',   'ZI271',  'ZI292',
         'ZI296',  'ZI303',  'ZI311N', 'ZI320',  'ZI321',  'ZI324',  'ZI332',
         'ZI333',  'ZI339',  'ZI341',  'ZI344',  'ZI348',  'ZI358',  'ZI364',
         'ZI365',  'ZI368',  'ZI378',  'ZI379',  'ZI384',  'ZI386',  'ZI388',
         'ZI398',  'ZI400',  'ZI402',  'ZI418N', 'ZI420',  'ZI437',  'ZI443',
         'ZI447',  'ZI455N', 'ZI456',  'ZI457',  'ZI460',  'ZI476',  'ZI477',
         'ZI486',  'ZI517',  'ZI523',  'ZI527',  'ZI85',  'ZI90']


sample_list = []

for sample in RG21:
    sample_list.append(pop_dict['RG'] + sample + '.fa')

for sample in ZI69:
    sample_list.append(pop_dict['ZI'] + sample + '.fa')

MD21 = ['MD03', 'MD06', 'MD105', 'MD106', 'MD146', 'MD15', 'MD197', 'MD199',
        'MD201', 'MD221', 'MD224', 'MD225', 'MD233', 'MD235', 'MD238',
        'MD243', 'MD251', 'MD255', 'MD63', 'MD72', 'MD73']

# the bedfile of dsim SIs, for getting the sim intron name
tbx_dsim_SI = pysam.TabixFile("../intersection_mel_sim_SIs/sim_SI.bed.gz")

# Get the GC content of all the short introns in mel and sim, using the bed files
# of SIs and the reference sequences
simref = Fasta('../Dsim_ref/v2.02/dsim-all-chromosome-r2.02.fasta')
melref = Fasta('../Dmel_ref/dmel-all-chromosome-r5.57.fasta')

melGCdict = {}
with open('../intersection_mel_sim_SIs/mel_SI.bed') as mel:
    for line in mel:
        l = line.strip().split()
        melSIcontig = l[0]
        melSIbedstart = int(l[1])
        melSIbedend = int(l[2])
        melSIname = l[3].split(';name=')[1]
        seq = melref[melSIcontig][melSIbedstart:melSIbedend].seq
        GCcontent = round(GC(seq) / 100, 4)
        melGCdict[melSIname] = GCcontent

simGCdict = {}
with open('../intersection_mel_sim_SIs/sim_SI.bed') as sim:
    for line in sim:
        l = line.strip().split()
        simSIcontig = l[0]
        simSIbedstart = int(l[1])
        simSIbedend = int(l[2])
        simSIname = l[3].split(';name=')[1]
        seq = simref[simSIcontig][simSIbedstart:simSIbedend].seq
        GCcontent = round(GC(seq) / 100, 4)
        simGCdict[simSIname] = GCcontent


wgabed = pysam.TabixFile('../alignment_Dmel_Dsim_Dyak/dmel.dsim.dyak.wga.bed.gz')

mel_fastas = {}
for x in sample_list:
    mel_fastas[x.split('/')[-1].strip('.fa')] = Fasta(x)

VCFs = {}
for x in ['Scf_2L', 'Scf_2R', 'Scf_3L', 'Scf_3R', 'Scf_X']:
    VCFs[x] = VariantFile('../Dsim_variant_calling/VCF/' + x + '.MD.g.vcf.gz')

counter=0

with open('../data/table.txt', 'wt') as f_out:
    f_out.write('melcontig\tmelsite\tsimcontig\tsimsite\tsimstrand\tmelref\tsimref\tyakref\t' + \
          '\t'.join(RG21) + '\t' + '\t'.join(ZI69) + '\t' + '\t'.join(MD21) + \
          '\tdmelname\tdsimname\tdmelGC\tdsimGC' + '\n')

    with open('../intersection_mel_sim_SIs/overlap.bed', 'rt') as f:
        for line in f:
            l = line.strip().split()

            # if there is no overlap, then skip this site:
            if l[4] == '.':
                continue

            # populate the things I'm interested in:
            contig = l[4]

# # # # # # # take the overlap of the two bed coordinates in the intersected file
            A_start = int(l[1])
            A_end = int(l[2])

            B_start = int(l[5])
            B_end = int(l[6])

            x = range(A_start, A_end)
            y = range(B_start, B_end)

            xy = sorted(set(x).intersection(y))

            bedstart = xy[0]
            if len(xy) == 1:
                bedend = bedstart + 1
            else:
                bedend = xy[-1] + 1
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

            Dmel_intron_name = l[3].split(';name=')[1]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


            # if the site is length > 1, NB. SKIP FOR NOW
            if (bedend - bedstart) != 1:
                 continue
                # print(line.strip())

    # # # # # test whether this site lies in the NR region(s) for the appropriate contig:
            skip = False
            if not contig in NRdict:
                continue

            for pair in NRdict[contig]:
                NRstart = int(pair.split(':')[0])
                NRend = int(pair.split(':')[1])
                if (bedstart + 1) >= NRstart and bedend <= NRend:
                    skip = True

            if skip:
                continue

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # # # # #  get the reference alleles:

    		# the outgroups Dmel, Dsim, DsimV2, Dyak are in the alignment, can discard DsimV2 for this):
            threeway = intersect2align(chromo = contig, start = bedstart, end = bedend,
                                      wga_bed = wgabed,
                                      ins_rel_ref=True)

    		# check that all three relevant outgroups are present (if not, skip this):
            if not all([x in threeway[0] for x in ['dsim', 'dmel', 'dyak']]):
                continue

    		# the order of the outgroups is: Dmel, Dsim, DsimV2, Dyak, but this doesn't
    		# matter if we make a dictionary with those indices in:
            which_spp = {y: x for x,y in enumerate(threeway[0])}
            threeway_dict = {threeway[0][which_spp['dsim']]: ''.join([x for x in threeway[1][which_spp['dsim']] if x != '-']),
    						 threeway[0][which_spp['dmel']]: ''.join([x for x, y in zip(threeway[1][which_spp['dmel']], threeway[1][which_spp['dsim']]) if y != '-']),
    						 threeway[0][which_spp['dyak']]: ''.join([x for x, y in zip(threeway[1][which_spp['dyak']], threeway[1][which_spp['dsim']]) if y != '-'])}

            # # exclude soft masked sites by checking all alleles are upper case:
            if not all([x in ['A', 'C', 'G', 'T'] for x in threeway_dict.values()]):
                continue

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # # # # #  get the mel sample alleles:

            mel_samples = [x.split('/')[-1].strip('.fa') for x in sample_list]
            mel_sequences = [get_from_fasta(x, contig, bedstart + 1, bedend) for x in sample_list]

            if not all([len(x) == 1 for x in mel_sequences]):
                sys.exit('there is a mel polymorphism allele entry with length > 1')

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # # # # #  get the sim sample alleles:
            # First, need to get the Dsim contig and Dsim coordinates somehow
            sim_contig = l[9].split(',')[which_spp['dsim']]
            sim_start = int(l[10].split(',')[which_spp['dsim']])
            sim_allele = l[11].split(',')[which_spp['dsim']]
            sim_strand = l[12].split(',')[which_spp['dsim']]
            if len(sim_allele) > 1:
                counter+=1
                continue
                # sys.exit('simulans allele is longer than 1')


            # test whether this site lies in the NR region(s) for the appropriate contig:
            skip_sim = False
            if not sim_contig in NRdictsim:
                continue

            if sim_contig != 'Scf_3R':
                for pair in NRdictsim[sim_contig]:
                    NRstartsim = int(pair.split(':')[0])
                    NRendsim = int(pair.split(':')[1])
                    if (sim_start + 1) >= NRstartsim and sim_start + 1 <= NRendsim:
                        skip_sim = True

            if skip_sim:
                continue

            ### NEED TO COMPLEMENT THE SIM POLY ALLELES IF THE STRAND IS NEGATIVE I THINK

            sim_dict = vcfchunk_to_dict(n = 42,
                                        variantfile = VCFs[sim_contig],
                                        contig = sim_contig,
                                        start = sim_start + 1, end = sim_start + 1,
                                        strand = sim_strand)

            if not all([len(x) == 1 for x in sim_dict.values()]):
                continue

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # # # # #  look up the Dsim intron name
            Dsim_intron_line = [str(x) for x in tbx_dsim_SI.fetch(sim_contig, sim_start, sim_start + 1)][0]
            Dsim_intron_name = Dsim_intron_line.split()[3].split(';name=')[1]
            # print(sim_start)
            # print(Dsim_intron_line)
            # print(Dmel_intron_name)
            # print()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # # # # #  get the GC content

            melGC = melGCdict[Dmel_intron_name]
            simGC = simGCdict[Dsim_intron_name]

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # # # # #  write the line

            f_out.write(contig + '\t' + str(bedstart + 1) + '\t' + sim_contig + '\t' + str(sim_start + 1) + '\t' + sim_strand + '\t' + \
                        threeway_dict['dmel'] + '\t' + threeway_dict['dsim'] + '\t' + threeway_dict['dyak'] + '\t' + \
                        '\t'.join(mel_sequences) + '\t' + '\t'.join([sim_dict[x] for x in MD21]) + '\t' + \
                        Dmel_intron_name + '\t' + Dsim_intron_name + '\t' + str(melGC) + '\t' + str(simGC) + '\n')




print('total sites i skipped because the sim allele was >length 1 was: ' + str(counter))






#

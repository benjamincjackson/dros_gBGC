This respository contains all the scripts necessary to replicate the analyses in the paper:

Jackson BC, Charlesworth B. Evidence for a force favouring GC over AT at short intronic sites in *Drosophila simulans* and *D. melanogaster*. *Submitted*

The final, filtered, dataset of genotypes at short intronic sites is available: `data/table.txt`

The pipeline for calling variants in the *D. sim* lines is available [here](https://github.com/benjamincjackson/dsim_variant_pipeline_ref_v2.02), and the resulting VCF files are available from the corresponding author on request.

The genomes for the *D. melanogaster* lines from the Drosophila Population Genomics Project/Drosophila Genome Nexus are publicly available [Lack et al (2015)](https://academic.oup.com/genetics/article/199/4/1229/5935880), [Lack et al (2016)](https://academic.oup.com/mbe/article/33/12/3308/2450097?login).
Scripts that were used to download and process the *D. mel* data are available under `Dmel_nexus/`




----

The rest of the repository is organised (alphabetically) as follows:

`Dmel_annotation/`

Scripts for downloading annotation for the *D. mel* genome from flybase, and for extracting the coordinates of short intronic (SI) sites.

`Dmel_nexus/`

Scripts for downloading Zambian (ZI) (and Rwandan, RG) lines of *D. mel* from the Drosophila Genome Nexus, for converting the data to FASTA format, and information about admixture tracts.

`Dmel_ref/`

Scripts for downloading the *D. mel* reference genome from flybase and for checking its repeat content.

`Dsim_annotation/`

Scripts for downloading annotation for the *D. sim* genome from flybase, and for extracting the coordinates of short intronic (SI) sites.

`Dsim_ref/`

Script for downloading the *D. sim* reference genome from flybase.

`Dyak_ref/`

Script for downloading the *D. yak* reference genome from flybase.

`PAML/`

Scripts and file-structure for running the PAML analyses described in the manuscript

`alignment_Dmel_Dsim_Dyak/`

Scripts for doing things with the multiple alignment between *mel*, *sim* and *yak*, from [Zeng et al (2019)](https://academic.oup.com/mbe/article/36/2/423/5182503?login=true). The multiple alignment was created by Henry Juho, using the pipeline described [here](https://github.com/henryjuho/threeway_fly_alignment).

`data/`

* `data/table.txt`

	Contains the final genotypes (including X-chromosome data) at homologous short intronic sites for ZI, RG and MD populations of *mel* and *sim*, the corresponding alleles at the *sim*, *mel* and *yak* reference sequences, the GC content in either species' 8-30bp region of each intron, and the names of the introns in either species.
	
* `data/<files>*.RData`

	Intermediate data produced by, and then used for the final analysis in, the scripts under `scripts/`
	
`est-sfs/`

Scripts and file-structure for running est-sfs as described in the manuscript

`glemin_gamma/`

R code from the supplement of [Glémin et al (2015)](https://genome.cshlp.org/content/25/8/1215.short), imported for some analyses under `scripts/`

`intersection_mel_sim_SIs/`

Scripts for finding the intersection of *mel* and *sim* short intronic sites, using the multiple alignment and the annotation of the two species' genomes.

`intersection_mel_sim_SIs_CHECK/`

checking the above

`manuscript/plots/`

R code to plot the figures in the main text of the manuscript, using *.RData files in `data/`

`scripts/`

R and python scripts to run the majority of the analyses.

 

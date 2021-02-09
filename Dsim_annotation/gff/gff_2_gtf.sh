gunzip dsim-all-r2.02.gff.gz

awk '!seen[$0]++' dsim-all-r2.02.gff > dsim-all-r2.02.uniq.gff

python2 ../../scripts/gffutils-flybase-convert.py --gff dsim-all-r2.02.uniq.gff \
						 --gtf dsim-all-r2.02.uniq.gtf \
						 --db dsim-all-r2.02.uniq.gtf.db

bgzip dsim-all-r2.02.uniq.gff.cleaned
tabix -pgff dsim-all-r2.02.uniq.gff.cleaned.gz

bgzip dsim-all-r2.02.gff
bgzip dsim-all-r2.02.uniq.gff
bgzip dsim-all-r2.02.uniq.gtf

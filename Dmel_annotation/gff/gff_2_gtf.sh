gunzip dmel-all-r5.57.gff.gz

awk '!seen[$0]++' dmel-all-r5.57.gff > dmel-all-r5.57.uniq.gff

python ../../scripts/gffutils-flybase-convert.py --gff dmel-all-r5.57.uniq.gff \
						 --gtf dmel-all-r5.57.uniq.gtf \
						 --db dmel-all-r5.57.uniq.gtf.db

bgzip dmel-all-r5.57.uniq.gff.cleaned
tabix -pgff dmel-all-r5.57.uniq.gff.cleaned.gz

bgzip dmel-all-r5.57.gff
bgzip dmel-all-r5.57.uniq.gff

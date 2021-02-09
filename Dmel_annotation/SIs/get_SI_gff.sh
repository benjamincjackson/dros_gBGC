python3 get_SI_gff.py > SIs.gff
sort -u -k1,1 -k4,4n -k5,5n SIs.gff > SIs.sorted.uniq.gff
bgzip SIs.sorted.uniq.gff
tabix -pgff SIs.sorted.uniq.gff.gz

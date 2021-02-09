# Get a wgabed file ordered by dsim v2.02 from the maf file using wgaBED utilities
# Need to do it by chromosome and then stick them all together

for CHR in Scf_2L Scf_2R Scf_3L Scf_3R Scf_4 Scf_X
do
  python /Users/ben/programs/WGAbed/maf_to_bed.py -i dmel.dsim.dyak.maf.gz -r dsim -c $CHR >> dsim.dmel.dyak.wga.bed
done

sort -k1,1 -k2,2n dsim.dmel.dyak.wga.bed | bgzip -c > dsim.dmel.dyak.wga.bed.gz
tabix dsim.dmel.dyak.wga.bed.gz

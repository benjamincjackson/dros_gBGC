# Get a wgabed file ordered by dsim v2.02 from the maf file using wgaBED utilities
# Need to do it by chromosome and then stick them all together

for CHR in 2L 2R 3L 3R 4 X
do
  python /Users/ben/programs/WGAbed/maf_to_bed.py -i dmel.dsim.dyak.maf.gz -r dmel -c $CHR >> dmel.dsim.dyak.wga.bed
done

sort -k1,1 -k2,2n dmel.dsim.dyak.wga.bed | bgzip -c > dmel.dsim.dyak.wga.bed.gz
tabix dmel.dsim.dyak.wga.bed.gz

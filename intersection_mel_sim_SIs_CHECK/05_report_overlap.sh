echo "total ZI SI sites: " `zcat < mel_SI.bed.gz | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'`
echo "total MD SI sites: " `zcat < sim_SI.bed.gz | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'`
echo "overlaps is: " `awk '{sum += $14} END {print sum}' overlap.bed`

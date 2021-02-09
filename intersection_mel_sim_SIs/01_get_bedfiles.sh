python3 get_sim_SI_bedfile.py > sim_SI.bed
bgzip -c sim_SI.bed > sim_SI.bed.gz
tabix -pbed sim_SI.bed.gz

python3 get_mel_SI_bedfile.py > mel_SI.bed
bgzip -c mel_SI.bed > mel_SI.bed.gz
tabix -pbed mel_SI.bed.gz

zcat < ../alignment_Dmel_Dsim_Dyak/dmel.dsim.dyak.wga.bed.gz |\
  python ~/programs/WGAbed/non_ref_intersect.py \
         -b sim_SI.bed.gz \
         -q dsim \
         -c Scf_2L > Scf_2L_alignment.wgabed

zcat < ../alignment_Dmel_Dsim_Dyak/dmel.dsim.dyak.wga.bed.gz |\
    python ~/programs/WGAbed/non_ref_intersect.py \
            -b sim_SI.bed.gz \
            -q dsim \
            -c Scf_2R > Scf_2R_alignment.wgabed

zcat < ../alignment_Dmel_Dsim_Dyak/dmel.dsim.dyak.wga.bed.gz |\
    python ~/programs/WGAbed/non_ref_intersect.py \
            -b sim_SI.bed.gz \
            -q dsim \
            -c Scf_3L > Scf_3L_alignment.wgabed

zcat < ../alignment_Dmel_Dsim_Dyak/dmel.dsim.dyak.wga.bed.gz |\
    python ~/programs/WGAbed/non_ref_intersect.py \
            -b sim_SI.bed.gz \
            -q dsim \
            -c Scf_3R > Scf_3R_alignment.wgabed

zcat < ../alignment_Dmel_Dsim_Dyak/dmel.dsim.dyak.wga.bed.gz |\
    python ~/programs/WGAbed/non_ref_intersect.py \
           -b sim_SI.bed.gz \
           -q dsim \
           -c Scf_X > Scf_X_alignment.wgabed

cat Scf_2L_alignment.wgabed Scf_2R_alignment.wgabed Scf_3L_alignment.wgabed Scf_3R_alignment.wgabed Scf_X_alignment.wgabed | sort -k1,1 -k2n,2 -k3n,3 > mel_sites_sim_coords.wga.bed
rm Scf_2L_alignment.wgabed Scf_2R_alignment.wgabed Scf_3L_alignment.wgabed Scf_3R_alignment.wgabed Scf_X_alignment.wgabed

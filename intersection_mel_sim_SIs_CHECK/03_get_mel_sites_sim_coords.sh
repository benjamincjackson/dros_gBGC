cat /Users/ben/Dropbox/Biology_Projects/dros_X_A_mut/alignment_Dmel_Dsim_Dsim_Dyak/dmel.dsim.dyak.dsimV2.wga.bed |\
  python /Users/ben/programs/WGAbed/non_ref_intersect.py \
         -b /Users/ben/Dropbox/Biology_Projects/dros_X_A_mut/intersection_mel_sim_SIs_CHECK/sim_SI.bed.gz \
         -q dsim \
         -c Scf_2L > Scf_2L_alignment.wgabed

 cat /Users/ben/Dropbox/Biology_Projects/dros_X_A_mut/alignment_Dmel_Dsim_Dsim_Dyak/dmel.dsim.dyak.dsimV2.wga.bed |\
    python /Users/ben/programs/WGAbed/non_ref_intersect.py \
            -b /Users/ben/Dropbox/Biology_Projects/dros_X_A_mut/intersection_mel_sim_SIs_CHECK/sim_SI.bed.gz \
            -q dsim \
            -c Scf_2R > Scf_2R_alignment.wgabed

cat /Users/ben/Dropbox/Biology_Projects/dros_X_A_mut/alignment_Dmel_Dsim_Dsim_Dyak/dmel.dsim.dyak.dsimV2.wga.bed |\
    python /Users/ben/programs/WGAbed/non_ref_intersect.py \
            -b /Users/ben/Dropbox/Biology_Projectsata/dros_X_A_mut/intersection_mel_sim_SIs_CHECK/sim_SI.bed.gz \
            -q dsim \
            -c Scf_3L > Scf_3L_alignment.wgabed

cat /Users/ben/Dropbox/Biology_Projects/dros_X_A_mut/alignment_Dmel_Dsim_Dsim_Dyak/dmel.dsim.dyak.dsimV2.wga.bed |\
    python /Users/ben/programs/WGAbed/non_ref_intersect.py \
            -b /Users/ben/Dropbox/Biology_Projects/dros_X_A_mut/intersection_mel_sim_SIs_CHECK/sim_SI.bed.gz \
            -q dsim \
            -c Scf_3R > Scf_3R_alignment.wgabed

cat /Users/ben/Dropbox/Biology_Projects/dros_X_A_mut/alignment_Dmel_Dsim_Dsim_Dyak/dmel.dsim.dyak.dsimV2.wga.bed |\
    python /Users/ben/programs/WGAbed/non_ref_intersect.py \
           -b /Users/ben/Dropbox/Biology_Projects/dros_X_A_mut/intersection_mel_sim_SIs_CHECK/sim_SI.bed.gz \
           -q dsim \
           -c Scf_X > Scf_X_alignment.wgabed

cat Scf_2L_alignment.wgabed Scf_2R_alignment.wgabed Scf_3L_alignment.wgabed Scf_3R_alignment.wgabed Scf_X_alignment.wgabed | sort -k1,1 -k2n,2 -k3n,3 > mel_sites_sim_coords.wga.bed

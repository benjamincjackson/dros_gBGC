mkdir -p output_files_GC_bins

for i in 01 02 03 04 05
do
  cp bin_${i}/run_2/rst output_files_GC_bins/rst_GC_bin_${i}
  cp bin_${i}/run_2/output output_files_GC_bins/output_GC_bin_${i}
done

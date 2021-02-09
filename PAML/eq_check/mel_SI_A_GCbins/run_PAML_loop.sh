for i in {1..5}
do
  mkdir bin_0${i}
  cp input_files_GC_bins/GC_bin_0${i}.PHYLIP.seq baseml.ctl dros_unrooted.tree loop.sh bin_0${i}/

  cd bin_0${i}
	sed -i '.bak' 's/seqfile =/seqfile = GC_bin_0'${i}'.PHYLIP.seq/' baseml.ctl
  ./loop.sh

  cd ..
done

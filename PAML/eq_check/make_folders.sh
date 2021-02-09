for FOL in ../mean_SI_A_GCbins/ ../mel_SI_A_GCbins/ ../sim_SI_A_GCbins/
do
  DIR=`echo $FOL | cut -d'/' -f2`
  mkdir $DIR
  mkdir $DIR/input_files_GC_bins/
  cp $FOL/input_files_GC_bins/* $DIR/input_files_GC_bins/
done

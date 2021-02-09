#!/bin/bash

echo $RANDOM > seed.txt

if [ ! -d output_files/ ]; then
  mkdir output_files
fi

# for the number of bins:
for i in {1..5}
do
	# then run the program
	~/programs/est-sfs-release-2.03/est-sfs config-kimura.txt input_files/GC_bin_0${i}.txt seed.txt output_files/GC_bin_0${i}.sfs output_files/GC_bin_0${i}.pvalues
done

# ~/programs/est-sfs-release-2.03/est-sfs config-kimura.txt input_files/ALL.txt  seed.txt output_files/ALL.sfs output_files/ALL.pvalues

#!/bin/bash


mkdir TMP/
mkdir -p outputfiles/
for i in {1..5}
do
  mkdir outputfiles/GCbin0${i}
  for INPUT in `find inputfiles/GCbin0${i}/*seq`
  do
    FILE=`echo $INPUT | rev | cut -d'/' -f1 | rev`
    NAME=`echo $FILE | cut -d'.' -f1`
    cp $INPUT baseml.ctl dros.tree TMP/
    cd TMP/
    # first, edit the control file to be specific to one bin's input file
    sed -i '.bak' "s/seqfile = /seqfile = $FILE/g" baseml.ctl
    # # then run the program
    ~/programs/paml4.8/bin/baseml baseml.ctl
    cp rst ../outputfiles/GCbin0${i}/${NAME}.rst
    cp output ../outputfiles/GCbin0${i}/${NAME}.output
    rm 2base.t baseml.ctl.bak dros.tree output rst1 baseml.ctl $FILE lnf rst rub
    cd ..
  done
done

rm -r TMP

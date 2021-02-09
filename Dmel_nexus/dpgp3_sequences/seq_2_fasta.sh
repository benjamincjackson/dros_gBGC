if ! $( ls | grep -q samples.txt )
then
  for FILE in `ls dpgp3_Chr2L/*seq`
  do
    SAMPLE=`echo $FILE | cut -d'/' -f2 | cut -d'_' -f1`
    echo $SAMPLE >> samples.txt
  done
fi

if ! $( ls -d */ | grep -q "fasta" )
then
  mkdir fasta
fi

while read SAMPLE
do
  echo -e ">2L" > fasta/${SAMPLE}.fa
  fold -w60 dpgp3_Chr2L/${SAMPLE}_Chr2L.seq >> fasta/${SAMPLE}.fa
  echo -e "\n>2R" >> fasta/${SAMPLE}.fa
  fold -w60 dpgp3_Chr2R/${SAMPLE}_Chr2R.seq >> fasta/${SAMPLE}.fa
  echo -e "\n>3L" >> fasta/${SAMPLE}.fa
  fold -w60 dpgp3_Chr3L/${SAMPLE}_Chr3L.seq >> fasta/${SAMPLE}.fa
  echo -e "\n>3R" >> fasta/${SAMPLE}.fa
  fold -w60 dpgp3_Chr3R/${SAMPLE}_Chr3R.seq >> fasta/${SAMPLE}.fa
  echo -e "\n>X" >> fasta/${SAMPLE}.fa
  fold -w60 dpgp3_ChrX/${SAMPLE}_ChrX.seq >> fasta/${SAMPLE}.fa
  samtools faidx fasta/${SAMPLE}.fa
done < samples.txt

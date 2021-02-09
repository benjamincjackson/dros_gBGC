for i in {1..10}
do
  mkdir run_${i}
  ~/programs/paml4.8/bin/baseml baseml.ctl
  mv 2base.t		output			rst1		lnf			rst			rub run_${i}/
  grep "lnL" run_${i}/output |  cut -d' ' -f7 >> lnlks.txt
done

#!/bin/bash
for k in *.csv; do
  newname=$k.sim;
  python3 ../simsetgen.py 4 $k $newname ;
  #echo $newname;
  #echo "Hello!";
#do echo "((echo 1; echo `pwd`/$k; echo `pwd`/tree.nwk; echo 1; echo d) | /usr/local/bin/HYPHYMP /usr/local/lib/hyphy/TemplateBatchFiles/BUSTED.bf )" | qsub -l walltime=48:00:00 ;
done

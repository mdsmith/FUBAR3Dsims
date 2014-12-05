#!/bin/bash
for k in $1/*.csv.sim.*; do
  simname=${PWD}/$k;
  faname=$simname.fasta;
  (echo $simname; echo $faname) | /usr/local/bin/HYPHYMP ${PWD}/f3dsims.bf;
  #echo "((echo $simname; echo $faname) | /usr/local/bin/HYPHYMP ${PWD}/../f3dsims.bf )" | qsub -l walltime=48:00:00 -q eternity ;
done

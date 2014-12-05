#!/bin/bash
for k in *.csv.sim.*; do
  simname=${PWD}/$k;
  faname=$simname.fasta;
  (echo $simname; echo $faname) | /usr/local/bin/HYPHYMP ${PWD}/../f3dsims.bf;
done

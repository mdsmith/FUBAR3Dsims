#!/bin/bash

#------------------------------------------------------------------------------#
#                             SIMULATIONS FOR FUBAR                            #
#                                   simpar.sh                                  #
#------------------------------------------------------------------------------#
#
# Usage: bash simpar.sh input.csv
#
# This script should take a space-delimited .csv file with no header as its
# argument. It runs a simulation for each line of this input file. input.csv
# should have the following 11 columns, in the order given:
#
# 1. csvfile: path to the binary input .csv for the simulation
#
# 2. simsdir: path to the directory containing FUBAR3Dsims files
#
# 3. ntaxa: number of taxa
#
# 4. mu1n: nonsyn mu_1
#
# 5. sigma1n: nonsyn sigma_1
#
# 6. mu2n: nonsyn mu_2
#
# 7. sigma2n: nonsyn sigma_2
#
# 8. mu1s: syn mu_1
#
# 9. sigma1s: syn sigma_1
#
# 10. mu2s: syn mu_2
#
# 11. sigma2s: syn sigma_2

#------------------------------------------------------------------------------#

while read csvfile simsdir ntaxa mu1n sigma1n mu2n sigma2n mu1s sigma1s mu2s \
sigma2s; 
  
  do newname=$csvfile.sim;
  python3 $simsdir/simsetgen.py $ntaxa $csvfile $newname \
    --nonsyn $mu1n $sigma1n $mu2n $sigma2n \
    --syn $mu1s $sigma1s $mu2s $sigma2s;

  for k in $csvfile.sim.*; do
    simname=$k;
    faname=$simname.fasta;
    (echo $simname; echo $faname) | /usr/local/bin/HYPHYMP $simsdir/f3dsims.bf;
    #echo "((echo $simname; echo $faname) | /usr/local/bin/HYPHYMP ${PWD}/../f3dsims.bf )" | qsub -l walltime=48:00:00 -q eternity ;
  done

  # The resulting FASTA files can be found at $csvfile.sim.*.fasta

done < $1

#------------------------------------------------------------------------------#

exit

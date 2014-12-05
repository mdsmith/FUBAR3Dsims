Simulations for FUBAR3D
=======================

HyPhy branch language sequence simulation scripts for the FUBAR3D model.

To generate settings files from csvs in subdirectory "data" (within this
directory):

    ./genSets.sh data

To generate fasta files from these settings files:

    ./runSets.sh data

To change the number of taxa in a tree, edit the genSets.sh file and change
the integer passed to simsetgen.py (it can be any number, but will be
unrooted).

To change the distributions from which syns and nonsyns are drawn, change the
parameters to simsetgen.py in genSets.sh to include

    --nonsyn mu1 sigma1 mu2 sigma2
and

    --syn mu1 sigma1 mu2 sigma2

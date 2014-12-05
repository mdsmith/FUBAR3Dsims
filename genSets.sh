#!/bin/bash
for k in *.csv; do
  newname=$k.sim;
  python3 ../simsetgen.py 4 $k $newname ;
done

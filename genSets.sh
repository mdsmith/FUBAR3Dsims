#!/bin/bash
for k in $1/*.csv; do
  newname=$k.sim;
  python3 simsetgen.py 16 $k $newname ;
done

#!/bin/bash


iteration=0
lastfile=2496
while [ $iteration -lt $lastfile ]; 
do
  iteration=$(( iteration + 1 ))
  echo $iteration
  infilename=(`head -n $iteration inputfile_8TeV | tail -1`)
  echo $infilename
  outfilename=lhehisto8TeV_$iteration.root
  echo $outfilename
  bsub -q 1nd $PWD/runJob.csh $infilename $outfilename
done

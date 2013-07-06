#!/bin/bash


iteration=0
lastfile=176
while [ $iteration -lt $lastfile ]; 
do
  iteration=$(( iteration + 1 ))
#  echo $iteration
  filename=(`head -n $iteration PDFFilesDY | tail -1`)
  echo $filename
  bsub -q 2nd $PWD/runPDF.csh $filename 0 
done

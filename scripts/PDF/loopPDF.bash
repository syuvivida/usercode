#!/bin/bash


iteration=0
lastfile=2
while [ $iteration -lt $lastfile ]; 
do
  iteration=$(( iteration + 1 ))
#  echo $iteration
  filename=(`head -n $iteration PDFFiles | tail -1`)
  echo $filename
  ./runPDF $filename 13 >& $iteration.Log &

done

#!/bin/bash


scriptname=`basename $0`
EXPECTED_ARGS=3

if [ $# -ne $EXPECTED_ARGS ]
then
    echo "Usage: $scriptname filelist logdirectory numberofjobs"
    echo "Example: ./$scriptname ttbar_input TTTo2L2Nu2B_8TeV-powheg 2160"
    exit 1
fi




iteration=0
lastfile=$3
while [ $iteration -lt $lastfile ]; 
do
  iteration=$(( iteration + 1 ))
  number=`grep -a "patTuple_muon_"$iteration"_" $1 | wc -l`
  rm -rf temp
  if [ $((number)) -gt 1 ]; then 
      echo $iteration "has problem"; 
      grep -a "patTuple_muon_"$iteration"_" $1 > temp
      tempiteration=0
      lasttemp=$((number))
#      echo $lasttemp
      while [ $tempiteration -lt $lasttemp ]; 
 	do
    	  tempiteration=$(( tempiteration + 1 ))
	  infilename=(`head -n $tempiteration temp  | tail -1`)
#	  echo $infilename
	  numbertemp=`grep -a $infilename $2/res/CMSSW_$iteration.stdout | wc -l`
	  if [ $((numbertemp)) -gt 0 ]; then 
	      echo $infilename " is really produced"
	  else
	      echo $infilename " should be removed" 
	      lcg-ls -l "srm://grid71.phy.ncu.edu.tw:8446/srm/managerv2?SFN=${infilename}"
# when you are sure which files to remove, uncomment the following line
#              lcg-del -l "srm://grid71.phy.ncu.edu.tw:8446/srm/managerv2?SFN=${infilename}"
	  fi
      done
  else if [ $((number)) -lt 1 ]; then    
      echo "Missing " $iteration ;
  fi
  fi
done

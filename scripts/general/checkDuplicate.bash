#!/bin/bash


iteration=0
lastfile=395
while [ $iteration -lt $lastfile ]; 
do
  iteration=$(( iteration + 1 ))
  number=`grep -a "/dpm/phy.ncu.edu.tw/home/cms/store/user/syu/DoubleMu/DoubleMu_A/1adea058a76608841d0b8b239b6c699b/patTuple_muon_"$iteration"_" muonoutput | wc -l`
  if [ $((number)) -gt 1 ]; then echo $iteration "has problem"; fi

done

#!/bin/csh

ls /mc/Laurent | grep -a Template.root | awk '{print "./runLoop.csh produceHisto /mc/Laurent/"$9 " "leading" "matching}' leading=$1 matching=$2 | bash
cd /mc/Laurent
mkdir $3
ls /mc/Laurent | grep -a Template_histo.root | awk '{print "mv /mc/Laurent/"$9 " /mc/Laurent/"dirname"/."}' dirname=$3 | bash
cd /mc/Laurent/$3
ls | grep -a histo.root | awk '{print "/mc/Laurent/scripts/runScan.csh "$9}' | bash

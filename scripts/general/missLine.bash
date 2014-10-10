#!/bin/bash


scriptname=`basename $0`
EXPECTED_ARGS=2

if [ $# -ne $EXPECTED_ARGS ]
then
    echo "Usage: $scriptname "
    echo "Example: ./$scriptname file A fileB"
    echo "Check if there is any lines in fileB but not in fileA"
    echo "the files don't have to be sorted" 
    exit 1
fi
awk 'NR==FNR{a[$0]=1;next}!a[$0]' $1 $2


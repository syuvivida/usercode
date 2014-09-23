#!/bin/bash

scriptname=`basename $0`
EXPECTED_ARGS=2

if [ $# -ne $EXPECTED_ARGS ]
then
    echo "Usage: ./$scriptname remote_directory string"
    echo "Example: ./$scriptname AbelianZPrime_ZH_lljj_M800-MADGRAPH flattuple"
    exit 1
fi

echo $1
string=$2
if [[ $1 == */* ]]
then
cmsLs /store/user/syu/$1 | grep -a syu | awk '{print "cmsLs "$5}' | bash | grep -a $string | awk -v my_var=$1 '{print "cmsRm "$5}' | bash
#cmsLs /store/user/syu/$1 | grep -a syu | awk '{print "cmsRmdir "$5}' | bash

#cmsRmdir /store/user/syu/$1

else
cmsLs /store/user/syu/$1 | grep -a syu | awk '{print "cmsLs "$5}' | bash | grep -a syu | awk '{print "cmsLs "$5}' | bash | grep -a $string | awk -v my_var=$1 '{print "cmsRm "$5}' | bash
#cmsLs /store/user/syu/$1 | grep -a syu | awk '{print "cmsLs "$5}' | bash | grep -a syu | awk '{print "cmsRmdir "$5}' | bash

#cmsLs /store/user/syu/$1 | grep -a syu | awk '{print "cmsRmdir "$5}' | bash

#cmsRmdir /store/user/syu/$1

fi

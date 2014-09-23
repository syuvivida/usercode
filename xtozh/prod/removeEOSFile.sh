#!/bin/bash

scriptname=`basename $0`
EXPECTED_ARGS=1

if [ $# -ne $EXPECTED_ARGS ]
then
    echo "Usage: ./$scriptname directory"
    echo "Example: ./$scriptname AbelianZPrime_ZH_lljj_M800-MADGRAPH"
    exit 1
fi

echo $1

cmsLs /store/user/syu/$1 | grep -a syu | awk '{print "cmsLs "$5}' | bash | grep -a syu | awk '{print "cmsLs "$5}' | bash | grep -a root | awk -v my_var=$1 '{print "cmsRm "$5}' | bash

cmsLs /store/user/syu/$1 | grep -a syu | awk '{print "cmsLs "$5}' | bash | grep -a syu | awk '{print "cmsRmdir "$5}' | bash

cmsLs /store/user/syu/$1 | grep -a syu | awk '{print "cmsRmdir "$5}' | bash

cmsRmdir /store/user/syu/$1

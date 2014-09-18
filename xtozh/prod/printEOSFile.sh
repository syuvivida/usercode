#!/bin/bash

scriptname=`basename $0`
EXPECTED_ARGS=2

if [ $# -ne $EXPECTED_ARGS ]
then
    echo "Usage: ./$scriptname directory outputFileName"
    echo "Example: ./$scriptname AbelianZPrime_ZH_lljj_M800-MADGRAPH filelist"
    exit 1
fi

echo $1


cmsLs /store/user/syu/$1 | grep -a syu | awk '{print "cmsLs "$5}' | bash | grep -a syu | awk '{print "cmsLs "$5}' | bash | grep -a root | awk '{print $5}' > $2



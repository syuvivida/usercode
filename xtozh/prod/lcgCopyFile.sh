#!/bin/bash

scriptname=`basename $0`
EXPECTED_ARGS=2

if [ $# -ne $EXPECTED_ARGS ]
then
    echo "Usage: ./$scriptname inutFile directoryname"
    echo "Example: ./$scriptname filelist \"//dpm/phy.ncu.edu.tw/home/cms/store/user/syu/DoublePhotonHighPt/DoubleEl_B/1c27adee06a9333ec3a821914a825c47/\""
    exit 1
fi

awk -F "${2}" '{print "lcg-cp \"srm://grid71.phy.ncu.edu.tw/'${2}'"$2"\" "$2}' ${1} | bash


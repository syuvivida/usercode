#!/bin/bash

scriptname=`basename $0`
EXPECTED_ARGS=3


if [ $# -ne $EXPECTED_ARGS ]
then
    echo "Usage: ./$scriptname remoteDirectory string userid"
    echo "Example: ./$scriptname TTJets_MassiveBinDECAY_TuneZ2star_8TeV flattuple syu"
    exit 1
fi

userid=$3

topdir="temp"
string=$2

if [[ ! -e $topdir ]]; then
    echo "creating directory "$topdir
    mkdir $topdir
fi


prefix="srm://f-dpm001.grid.sinica.edu.tw:8446"
echo $prefix
dpmprefix=$prefix"/dpm/grid.sinica.edu.tw/home/cms/store/user/"$userid"/"$1
echo $dpmprefix



if [[ $1 == */* ]];
then
    removeprefix=`lcg-ls $dpmprefix`
    echo $removeprefix
    dpmprefix=$prefix$removeprefix
    echo $dpmprefix
    
    lcg-ls $dpmprefix | grep -a $string |  awk -F $removeprefix -v my_var=$dpmprefix -v my_dir=$topdir"/"$1 '{print "lcg-del -l -v "my_var"/"$2}' | bash

else
    removeprefix=`lcg-ls $dpmprefix | grep -a $userid |  awk -v my_var=$prefix '{print "lcg-ls "my_var""$1}' | bash`
    echo $removeprefix
    dpmprefix=$prefix$removeprefix
    echo $dpmprefix
    
    lcg-ls $dpmprefix | grep -a $string |  awk -F $removeprefix -v my_var=$dpmprefix -v my_dir=$topdir"/"$1 '{print "lcg-del -l -v "my_var"/"$2}' | bash

fi




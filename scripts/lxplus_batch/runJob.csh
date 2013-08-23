#!/bin/tcsh

setenv WORKDIR /afs/cern.ch/work/s/syu/Markus/CMSSW_5_3_9/src
cd $WORKDIR

setenv SCRAM_ARCH slc5_amd64_gcc462; eval `scramv1 runtime -csh`

cd $WORKDIR/runJob

cmsRun dumpLHE_cfg.py inputFiles="${1}" histoFile="$2"





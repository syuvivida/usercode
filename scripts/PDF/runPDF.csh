#!/bin/tcsh

cd /afs/cern.ch/work/s/syu/PDFReweigh/CMSSW_5_2_5/src/PDF;
setenv SCRAM_ARCH slc5_amd64_gcc462; eval `scramv1 runtime -csh`
setenv LHAPATH /afs/cern.ch/work/s/syu/jetsub/LHAPDF/share/lhapdf/PDFsets/
echo $LHAPATH
root -q -b runPDF.C\(\"myLHAPDF_reweighing\",\"$1\",$2\)

echo "Job finished"

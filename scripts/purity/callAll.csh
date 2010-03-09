#!/bin/csh
rm -rf event.dat
./exec.csh /home/syu/gammajet_pt20_30.root
./exec.csh /home/syu/gammajet_pt30_50.root
./exec.csh /home/syu/gammajet_pt50_80.root
#./exec.csh /home/syu/qcd_pt20_30.root
#./exec.csh /home/syu/qcd_pt30_80.root


#./runLoop.csh loopEntries /home/syu/gammajet_pt20_30.root
#./runLoop.csh loopEntries /home/syu/gammajet_pt30_50.root
#./runLoop.csh loopEntries /home/syu/gammajet_pt50_80.root
#./runLoop.csh loopEntries /home/syu/qcd_pt20_30.root
#./runLoop.csh loopEntries /home/syu/qcd_pt30_80.root

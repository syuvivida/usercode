cd $CMSSW_BASE/src

cvs co -d RandomConeAna UserCode/EYKim/photonStudy/RandomConeAna
cvs co -r $1 RecoEgamma/PhotonIdentification
cvs co -r $1 RecoEgamma/EgammaIsolationAlgos

# patch
mv RandomConeAna/patches/PhotonIsolationCalculator.cc  $CMSSW_BASE/src/RecoEgamma/PhotonIdentification/src
mv RandomConeAna/patches/PhotonIsolationCalculator.h   $CMSSW_BASE/src/RecoEgamma/PhotonIdentification/interface
mv RandomConeAna/patches/EgammaRecHitIsolation.cc      $CMSSW_BASE/src/RecoEgamma/EgammaIsolationAlgos/src
mv RandomConeAna/patches/EgammaRecHitIsolation.h       $CMSSW_BASE/src/RecoEgamma/EgammaIsolationAlgos/interface
mv RandomConeAna/patches/isolationCalculator_cfi.py    $CMSSW_BASE/src/RecoEgamma/PhotonIdentification/python
mv RandomConeAna/patches/EgammaEcalRecHitIsolationProducer.h   $CMSSW_BASE/src/RecoEgamma/EgammaIsolationAlgos/plugins
mv RandomConeAna/patches/EgammaEcalRecHitIsolationProducer.cc  $CMSSW_BASE/src/RecoEgamma/EgammaIsolationAlgos/plugins


rm RandomConAna/patches -rf



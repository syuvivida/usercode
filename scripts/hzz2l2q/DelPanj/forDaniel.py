import FWCore.ParameterSet.Config as cms

isMergedJet = True

## Below are default values for the names of higgs collections
## and whether the analysis is a merged-jet analysis or not
higgsCollectionEE = "hzzeejj:h"
higgsCollectionMM = "hzzmmjj:h"
MERGED            = False


if isMergedJet:
  higgsCollectionEE = "hzzee1j:h"
  higgsCollectionMM = "hzzmm1j:h"
  MERGED            = True


process = cms.Process("COMBINE")
#runOnMC = False # to run on DATA
runOnMC = True

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring(
     'file:/afs/cern.ch/work/s/syu/h2l2qSkimData.root'
#    'file:/data4/syu/patsample/53X_MC/h2l2qSkimData_10_1_CL3.root'
#    'file:/scratch/syu/update_PU/CMSSW_5_3_3_patch3/src/runJob/testData.root'
    )
 )

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1))

process.load("RecoBTag.PerformanceDB.PoolBTagPerformanceDB062012")
process.load("RecoBTag.PerformanceDB.BTagPerformanceDB062012")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

baseJetSel = cms.PSet(
    Jets=cms.InputTag("cleanPatJetsNoPUIsoLept")
)


from DelPanj.TreeMaker.eSel2012HZZ_cff import *
from DelPanj.TreeMaker.muSel2012HZZ_cff import *

process.tree = cms.EDAnalyzer(
	'TreeMaker',
	fillPUweightInfo_ = cms.bool(True),
	fillEventInfo_ = cms.bool(True),
	fillGenInfo_   = cms.bool(True),
	fillMuonInfo_  = cms.bool(True),
	fillElecInfo_  = cms.bool(True),
	fillElecIsoInfo_ = cms.bool(False),
	fillJetInfo_   = cms.bool(True),
	fillMetInfo_   = cms.bool(False),
	fillTrigInfo_  = cms.bool(False),
	fillPhotInfo_  = cms.bool(False),
	fillZZInfo_    = cms.bool(True),
#        hzzeejjTag     = cms.InputTag("hzzeejj:h"),
#        hzzmmjjTag     = cms.InputTag("hzzmmjj:h"),
        hzzeejjTag     = cms.InputTag(higgsCollectionEE),
        hzzmmjjTag     = cms.InputTag(higgsCollectionMM),
        merged_        = cms.bool(MERGED),
	genPartLabel   = cms.InputTag("genParticles"),
	patMuons       = cms.InputTag("userDataSelectedMuons"),
	patElectrons   = cms.InputTag("userDataSelectedElectrons"),
        e2012IDSet     = eSel2012HZZ,
        mu2012IDSet    = muSel2012HZZ,
        eleRhoIso      = cms.InputTag("kt6PFJets","rho"),
        muoRhoIso      = cms.InputTag("kt6PFJetsCentralNeutral", "rho"),
	patMet         = cms.InputTag("patMETs"),
	beamSpotLabel  = cms.InputTag("offlineBeamSpot"),
	patJetPfAk05   = baseJetSel,
	outFileName    = cms.string('outputFileName.root')      
)


process.counter_original = cms.EDAnalyzer('MyEventCounter',
   instance = cms.int32(1) 
)


process.TFileService = cms.Service("TFileService",
      fileName = cms.string("hzz2l2q.root"),
      closeFileFast = cms.untracked.bool(True)
  )


process.p = cms.Path(
	process.counter_original*
        process.tree##Trigger Applied.
	)






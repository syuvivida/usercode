import FWCore.ParameterSet.Config as cms

process = cms.Process("COMBINE")
#runOnMC = False # to run on DATA
runOnMC = True

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring(
    'file:/afs/cern.ch/work/s/syu/GluGluToHToZZTo2L2Q_M-300_8TeV-powheg-pythia6/h2l2qSkimData_2_1_Taz.root',
    'file:/afs/cern.ch/work/s/syu/GluGluToHToZZTo2L2Q_M-300_8TeV-powheg-pythia6/h2l2qSkimData_4_1_1I2.root',
    'file:/afs/cern.ch/work/s/syu/GluGluToHToZZTo2L2Q_M-300_8TeV-powheg-pythia6/h2l2qSkimData_20_1_h3C.root'   
    )
 )

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1))

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("RecoBTag.PerformanceDB.PoolBTagPerformanceDB062012")
process.load("RecoBTag.PerformanceDB.BTagPerformanceDB062012")

process.GlobalTag.globaltag = cms.string('START50_V15::All')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

baseJetSel = cms.PSet(
  Jets=cms.InputTag("cleanPatJetsNoPUIsoLept")
)

from DelPanj.TreeMaker.eSelLvdp2011_cff import *
from DelPanj.TreeMaker.muSelLvdp2011_cff import *

from DelPanj.TreeMaker.eSel2012HZZ_cff import *
from DelPanj.TreeMaker.muSel2012HZZ_cff import *

from DelPanj.TreeMaker.eSel2012Tight_cff import *
from DelPanj.TreeMaker.muSel2012NoIso_cff import *

process.tree = cms.EDAnalyzer(
	'TreeMaker',
	fillPUweightInfo_ = cms.bool(True),
	fillEventInfo_ = cms.bool(True),
	fillGenInfo_   = cms.bool(True),
	fillMuonInfo_  = cms.bool(False),
	fillElecInfo_  = cms.bool(False),
	fillElecIsoInfo_ = cms.bool(False),
	fillJetInfo_   = cms.bool(False),
	fillMetInfo_   = cms.bool(False),
	fillTrigInfo_  = cms.bool(False),
	fillPhotInfo_  = cms.bool(False),
	fillZJetPlant_ = cms.bool(False),
	fillZZInfo_    = cms.bool(True),
        hzzeejjTag = cms.InputTag("hzzeejj:h"),
        hzzmmjjTag = cms.InputTag("hzzmmjj:h"),
	genPartLabel=cms.InputTag("genParticles"),
	patMuons=cms.InputTag("userDataSelectedMuons"),
	patElectrons = cms.InputTag("userDataSelectedElectrons"),
	leadElecPset_ = eSelLvdp2011,
	subLeadElecPset_ = eSelLvdp2011,
        e2012IDSet  =  eSel2012HZZ,
        mu2012IDSet = muSel2012HZZ,
##### debug
#        e2012TagSet = eSel2012Tight,
#        mu2012NoIsoSet = muSel2012NoIso,
#####
        eleRhoIso = cms.InputTag("kt6PFJetsForIso","rho"),
        muoRhoIso = cms.InputTag("kt6PFJetsCentralNeutral", "rho"),
	patMet=cms.InputTag("patMETs"),
	beamSpotLabel=cms.InputTag("offlineBeamSpot"),
	patJetPfAk05 = baseJetSel,
	outFileName=cms.string('outputFileName.root')
	)



process.counter_original = cms.EDAnalyzer('EventCounter',
   instance = cms.int32(1) 
)


process.TFileService = cms.Service("TFileService",
      fileName = cms.string("hzz2l2q_M300.root"),
      closeFileFast = cms.untracked.bool(True)
  )


process.p = cms.Path(
	process.counter_original*
	process.tree##Trigger Applied.
	)
  
 




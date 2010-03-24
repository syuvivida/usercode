import FWCore.ParameterSet.Config as cms


process = cms.Process("ANA")

process.load("FWCore.MessageLogger.MessageLogger_cfi")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
       '/store/relval/CMSSW_3_5_0_pre1/RelValTTbar/GEN-SIM-RECO/STARTUP3X_V14-v1/0006/14920B0A-0DE8-DE11-B138-002618943926.root'
    
    )
                            )

process.GenTrig = cms.EDAnalyzer("GenTrig",
        genLabel = cms.InputTag( "genParticles"),
        HLTLabel = cms.InputTag( "TriggerResults::HLT")
 )


process.GenTrig.dumpHEP = cms.untracked.bool(False);
process.GenTrig.pdgCode = cms.untracked.int32(22);


process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.TFileService = cms.Service("TFileService", fileName = cms.string('debug.root'))

process.p = cms.Path(process.GenTrig)



import FWCore.ParameterSet.Config as cms


process = cms.Process("ANA")

process.load("FWCore.MessageLogger.MessageLogger_cfi")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
       '/store/relval/CMSSW_3_5_0_pre2/RelValGammaJets_Pt_80_120/GEN-SIM-RECO/STARTUP3X_V14-v1/0010/D2ECC336-22EE-DE11-8A55-002618943920.root'
    )
                            )

process.RECOTrigger = cms.EDAnalyzer("RECOTrigger",
        genLabel = cms.InputTag( "genParticles"),
        phoLabel = cms.InputTag( "photons"),
        HLTLabel = cms.InputTag( "TriggerResults::HLT")
                                     )


process.RECOTrigger.dumpHEP = cms.untracked.bool(False);
process.RECOTrigger.pdgCode = cms.untracked.int32(22);

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('START3X_V26::All')
# process.GlobalTag.globaltag = cms.string('GR10_P_V4::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.TFileService = cms.Service("TFileService", fileName = cms.string('mpamatch_rereco.root'))

                                        
# let it run
process.p = cms.Path(
            process.RECOTrigger
            )



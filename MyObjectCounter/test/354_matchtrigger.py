import FWCore.ParameterSet.Config as cms


process = cms.Process("ANA")

process.load("FWCore.MessageLogger.MessageLogger_cfi")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(

    '/store/relval/CMSSW_3_3_6/RelValGammaJets_Pt_80_120/GEN-SIM-RECO/STARTUP3X_V8H-v1/0009/BE03C8D0-9EE4-DE11-80C8-002618943964.root',

#   '/store/relval/CMSSW_3_5_0_pre2/RelValGammaJets_Pt_80_120/GEN-SIM-RECO/STARTUP3X_V14-v1/0010/D2ECC336-22EE-DE11-8A55-002618943920.root',
#         Below is 336patch4 900 GeV data
#    'rfio:/castor/cern.ch/user/h/hbrun/BSCSkim_MinBiasFilter/data/BSCFilter_Jan29_v8_MinBiasFilter_Run123596_88.root'
    )
                            )


process.DataRECOTrigger = cms.EDAnalyzer("DataRECOTrigger",
        genLabel = cms.InputTag( "genParticles"),
        phoLabel = cms.InputTag( "photons"),
        HLTLabel = cms.InputTag( "TriggerResults::HLT")
 )


process.DataRECOTrigger.dumpHEP = cms.untracked.bool(False);
process.DataRECOTrigger.pdgCode = cms.untracked.int32(22);

process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.TFileService = cms.Service("TFileService", fileName = cms.string('trig.root'))

process.p = cms.Path(
                     process.DataRECOTrigger)



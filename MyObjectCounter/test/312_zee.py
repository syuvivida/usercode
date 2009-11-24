import FWCore.ParameterSet.Config as cms


process = cms.Process("Zee")

process.load("FWCore.MessageLogger.MessageLogger_cfi")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_1_2/RelValZEE/GEN-SIM-RECO/MC_31X_V3-v1/0007/0C3128C6-A878-DE11-9CEE-001D09F25208.root',
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_1_2/RelValZEE/GEN-SIM-RECO/MC_31X_V3-v1/0007/4666B35F-9278-DE11-B918-000423D8F63C.root',
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_1_2/RelValZEE/GEN-SIM-RECO/MC_31X_V3-v1/0007/900FF494-9278-DE11-BA53-001D09F28F1B.root',
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_1_2/RelValZEE/GEN-SIM-RECO/MC_31X_V3-v1/0007/A838D597-9078-DE11-AEDD-000423D99896.root',
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_1_2/RelValZEE/GEN-SIM-RECO/MC_31X_V3-v1/0007/C6B66CD2-A978-DE11-A5A2-000423D98BC4.root',
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_1_2/RelValZEE/GEN-SIM-RECO/MC_31X_V3-v1/0007/F0303A91-9278-DE11-AADC-001D09F25456.root'
    )
                            )


process.Zee = cms.EDAnalyzer("Zee",
        eleLabel = cms.InputTag( "gsfElectrons"),
        phoLabel = cms.InputTag( "photons"),
        genLabel = cms.InputTag( "genParticles" )
 )
process.Zee.dumpHEP = cms.untracked.bool(False);


process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.TFileService = cms.Service("TFileService", fileName = cms.string('zee_new.root'))

process.p = cms.Path(#process.patDefaultSequence*
                     process.Zee)



import FWCore.ParameterSet.Config as cms


process = cms.Process("PAT")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('MC_31X_V3::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#    'file:/mc/fakeev/SingleElectronPt100.root'
#     'file:/mc/fakeev/PATuple/PAT_nocheckcharge_pt20.root'
#    'file:/mc/Photon-Jet/PATuple/QCD_photon_test.root'
#     'rfio:/castor/cern.ch/user/g/gaultney/SkimsTest/test.root'
    'file:/home/yunju/CMSSW/CMSSW_3_1_1/src/Photon_Pargun_0906/NoAnti_SinglePhotonFlatPt90_eta06_906.root'
    )
                            )

process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.MyObjectCounter = cms.EDAnalyzer("MyObjectCounter",
        muoLabel = cms.InputTag( "cleanLayer1Muons" ),
        eleLabel = cms.InputTag( "cleanLayer1Electrons" ),
        phoLabel = cms.InputTag( "cleanLayer1Photons" ),
        triggerEvent = cms.InputTag( "patTriggerEvent" ),
        trigger      = cms.InputTag( "patTrigger"),                  matcherName   = cms.string( "myphotonTriggerL1Match")                                 
 )
process.MyObjectCounter.dumpHEP = cms.untracked.bool(False);




process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.TFileService = cms.Service("TFileService", fileName = cms.string('results_eiko.root'))

process.p = cms.Path(process.patDefaultSequence*
                     process.MyObjectCounter)



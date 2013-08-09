import FWCore.ParameterSet.Config as cms

process = cms.Process("dumpLHE")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:MCDBtoEDM_NONE.root')
)

process.dummy = cms.EDAnalyzer("DummyLHEAnalyzer",
    src = cms.InputTag("source"),
    histoutputFile= cms.untracked.string('gen.root')          
)

process.p = cms.Path(process.dummy)



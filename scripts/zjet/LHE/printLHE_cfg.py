import FWCore.ParameterSet.Config as cms

process = cms.Process("printLHE")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:bcmass_test_NONE.root')
)

process.dummy = cms.EDAnalyzer("PrintLHEAnalyzer",
    src = cms.InputTag("source"),
    histoutputFile= cms.untracked.string('print.root')          
)

process.p = cms.Path(process.dummy)



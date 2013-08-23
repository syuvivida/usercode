import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register ('histoFile',
                  '',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Name of root file that stores histograms")

options.parseArguments()

process = cms.Process("dumpLHE")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("LHESource",
			    fileNames  = cms.untracked.vstring(options.inputFiles)
                            )

process.dummy = cms.EDAnalyzer("DummyLHEAnalyzer",
    src = cms.InputTag("source"),
    histoutputFile= cms.untracked.string(options.histoFile)          
)

process.p = cms.Path(process.dummy)



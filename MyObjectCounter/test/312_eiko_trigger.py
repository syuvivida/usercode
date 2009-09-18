import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

process.load("FWCore.MessageLogger.MessageLogger_cfi")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
     'rfio:/castor/cern.ch/user/g/gaultney/SkimsTest/test.root'
    )
                            )


process.myphotonTriggerHLTMatch = cms.EDFilter( "PATTriggerMatcherDRDPtLessByR",
    src     = cms.InputTag( "cleanLayer1Photons" ),
    matched = cms.InputTag( "patTrigger" ),
    andOr          = cms.bool( False ),
    filterIdsEnum  = cms.vstring( 'TriggerPhoton' ),
    filterIds      = cms.vint32( 0 ),
    filterLabels   = cms.vstring( '*' ),
    pathNames      = cms.vstring( '*' ),
    collectionTags = cms.vstring( '*' ),
    maxDPtRel = cms.double( 1.0 ),
    maxDeltaR = cms.double( 0.2 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( True )
)


process.myphotonTriggerL1Match = cms.EDFilter( "PATTriggerMatcherDRDPtLessByR",
    src     = cms.InputTag( "cleanLayer1Photons" ),
    matched = cms.InputTag( "patTrigger" ),
    andOr          = cms.bool( False ),
    filterIdsEnum  = cms.vstring( 'TriggerL1NoIsoEG' ),
    filterIds      = cms.vint32( 0 ),
    filterLabels   = cms.vstring( '*' ),
    pathNames      = cms.vstring( '*' ),
    collectionTags = cms.vstring( '*' ),
    maxDPtRel = cms.double( 1.0 ),
    maxDeltaR = cms.double( 0.2 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( True )
)

process.patTriggerPhotonMatcher = cms.Sequence(
   process.myphotonTriggerL1Match +
   process.myphotonTriggerHLTMatch
)

process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff" )
#process.patTriggerMatcher += process.patTriggerPhotonMatcher
process.patTriggerMatcher += process.myphotonTriggerL1Match
process.patTriggerMatcher.remove( process.patTriggerElectronMatcher )
process.patTriggerMatcher.remove( process.patTriggerMuonMatcher )
process.patTriggerMatcher.remove( process.patTriggerTauMatcher )


process.patTriggerEvent.patTriggerMatches = [ "myphotonTriggerL1Match"]
## configure pat trigger
process.patTrigger.onlyStandAlone = False

## add trigger specific event content to PAT event content
#process.out.outputCommands += patTriggerEventContent
#for matchLabel in process.patTriggerEvent.patTriggerMatches:
#        process.out.outputCommands += [ 'keep patTriggerObjectsedmAssociation_patTriggerEvent_' + matchLabel + '_*' ]



process.MyObjectCounter = cms.EDAnalyzer("MyObjectCounter",
        muoLabel = cms.InputTag( "cleanLayer1Muons" ),
        eleLabel = cms.InputTag( "cleanLayer1Electrons" ),
        phoLabel = cms.InputTag( "cleanLayer1Photons" ),
        triggerEvent = cms.InputTag( "patTriggerEvent" ),
        trigger      = cms.InputTag( "patTrigger"),
#     matcherName = cms.string( "patTriggerPhotonMatcher" )
      matcherName = cms.string( "myphotonTriggerL1Match")                                    
 )
process.MyObjectCounter.dumpHEP = cms.untracked.bool(False);

process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.TFileService = cms.Service("TFileService", fileName = cms.string('results_eiko.root'))

process.p = cms.Path(
                     process.patTriggerSequence*
                     process.MyObjectCounter
                     )



# Auto generated configuration file
# using: 
# Revision: 1.303.2.7 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: PhotonJet_Pt_80_120_7TeV_cfi.py -s GEN --conditions MC_42_V17::All --datatier GEN-SIM --eventcontent RAWSIM --no_exec -n 10
import FWCore.ParameterSet.Config as cms

process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic7TeV2011Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(500000)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.303.2.7 $'),
    annotation = cms.untracked.string('zJet.py nevts:10'),
    name = cms.untracked.string('PyReleaseValidation')
)



# Other statements
process.GlobalTag.globaltag = 'MC_44_V13::All'

process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    displayPythiaCards   = cms.untracked.bool(True),				 
    comEnergy = cms.double(7000.0),				 
    maxEventsToPrint = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(
	pythiaUESettings = cms.vstring(
	'MSTU(21)=1     ! Check on possible errors during program execution', 
	'MSTJ(22)=2     ! Decay those unstable particles', 
	'PARJ(71)=10 .  ! for which ctau  10 mm', 
	'MSTP(33)=0     ! no K factors in hard cross sections', 
	'MSTP(2)=1      ! which order running alphaS', 
	'MSTP(51)=10042 ! structure function chosen (external PDF CTEQ6L1)',
	'MSTP(52)=2     ! work with LHAPDF',

	'PARP(82)=1.832 ! pt cutoff for multiparton interactions', 
	'PARP(89)=1800. ! sqrts for which PARP82 is set', 
	'PARP(90)=0.275 ! Multiple interactions: rescaling power', 

        'MSTP(95)=6     ! CR (color reconnection parameters)',
        'PARP(77)=1.016 ! CR',
        'PARP(78)=0.538 ! CR',

	'PARP(80)=0.1   ! Prob. colored parton from BBR',
	
	'PARP(83)=0.356 ! Multiple interactions: matter distribution parameter', 
	'PARP(84)=0.651 ! Multiple interactions: matter distribution parameter', 

	'PARP(62)=1.025 ! ISR cutoff', 
	
	'MSTP(91)=1     ! Gaussian primordial kT', 
	'PARP(93)=10.0  ! primordial kT-max', 
	'MSTP(81)=0     ! multiple parton interactions 1 is Pythia default',
	'MSTP(111)=0    ! fragmentation',
	'MSTP(82)=4     ! Defines the multi-parton model', 
	),
        processParameters = cms.vstring(
	'MSEL = 0        ! user defined processes',
	'MSUB(15) = 1    ! ff -> Z0 f',
	'MSUB(30) = 1    ! ff -> Z0 g',
	'MDME(174,1) = 0 ! Z decay into d dbar',
	'MDME(175,1) = 0 ! Z decay into u ubar',
	'MDME(176,1) = 0 ! Z decay into s sbar',
	'MDME(177,1) = 0 ! Z decay into c cbar',
	'MDME(178,1) = 0 ! Z decay into b bbar',
	'MDME(179,1) = 0 ! Z decay into t tbar',
	'MDME(182,1) = 1 ! Z decay into e- e+',
	'MDME(183,1) = 0 ! Z decay into nu_e nu_ebar',
	'MDME(184,1) = 0 ! Z decay into mu- mu+',
	'MDME(185,1) = 0 ! Z decay into nu_mu nu_mubar',
	'MDME(186,1) = 0 ! Z decay into tau- tau+',
	'MDME(187,1) = 0 ! Z decay into nu_tau nu_taubar' ,
	'CKIN(3) = 20    ! minimum pt hat for hard interactions',
	'CKIN(4) = -1    ! maximum pt hat for hard interactions',
	),
		
        parameterSets = cms.vstring('pythiaUESettings', 
				    'processParameters')
	)
				 )



process.TFileService = cms.Service("TFileService", fileName =
				   cms.string('7TeV_pythia_zJet_noMPIHad.root'))

baseJetSel = cms.PSet(
  Jets=cms.InputTag("selectedPatJetsPFlow")
)

from DelPanj.TreeMaker.eSelLvdp2011_cff import *
process.tree = cms.EDAnalyzer(
	'TreeMaker',
	fillPUweightInfo_ = cms.bool(False),
	DontDoPUReweight_ = cms.bool(True),#MC
	fillEventInfo_ = cms.bool(False),
	fillGenInfo_   = cms.bool(True),
	fillMuonInfo_  = cms.bool(False),
	fillElecInfo_  = cms.bool(False),
	fillElecIsoInfo_ = cms.bool(False),
	fillJetInfo_   = cms.bool(False),
	fillMetInfo_   = cms.bool(False),
	fillTrigInfo_  = cms.bool(False),
	fillPhotInfo_  = cms.bool(False),
	fillZJetPlant_ = cms.bool(False),
	genPartLabel=cms.InputTag("genParticles"),
	patMuons=cms.InputTag("selectedPatMuonsPFlow"),
	patElectrons = cms.InputTag("selectedPatElectronsPFlow"),
	leadElecPset_ = eSelLvdp2011,
	subLeadElecPset_ = eSelLvdp2011,
	patJetLabel =cms.InputTag("selectedPatJetsPFlow"),
	patMet=cms.InputTag("patMETs"),
	beamSpotLabel=cms.InputTag("offlineBeamSpot"),
	patJetPfAk05 = baseJetSel,
	outFileName=cms.string('outputFileName.root')
	)



# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.end_ana     = cms.Path(process.tree)


# Schedule definition
#process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.RAWSIMoutput_step)

process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,
				process.end_ana)

# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

import FWCore.ParameterSet.Config as cms

process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
readFiles = cms.untracked.vstring()

process.source = cms.Source("LHESource",firstEvent = cms.untracked.uint32(0),skipEvents = cms.untracked.uint32(0),
    fileNames =readFiles

)

readFiles.extend( [
'file:/data4/syu/madgraph_lhe/lhe_files/100_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/101_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/102_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/103_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/104_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/105_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/106_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/107_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/108_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/109_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/10_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/110_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/111_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/112_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/113_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/114_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/115_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/116_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/117_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/118_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/119_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/11_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/120_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/121_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/122_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/123_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/124_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/125_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/126_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/127_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/128_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/129_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/12_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/130_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/131_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/132_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/133_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/134_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/135_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/136_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/137_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/138_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/139_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/13_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/140_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/141_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/142_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/143_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/144_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/145_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/146_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/147_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/148_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/149_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/14_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/150_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/151_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/152_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/153_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/154_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/155_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/156_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/157_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/158_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/159_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/15_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/160_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/161_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/162_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/163_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/164_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/165_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/166_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/167_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/168_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/169_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/170_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/171_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/172_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/173_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/174_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/175_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/176_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/177_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/178_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/179_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/17_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/180_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/181_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/182_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/183_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/184_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/185_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/186_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/187_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/188_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/189_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/18_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/190_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/191_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/192_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/193_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/194_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/195_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/196_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/197_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/198_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/199_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/1_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/200_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/201_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/202_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/203_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/205_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/206_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/207_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/208_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/209_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/20_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/210_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/211_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/212_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/213_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/214_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/215_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/216_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/217_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/218_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/219_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/21_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/220_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/221_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/222_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/223_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/224_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/225_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/226_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/227_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/228_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/229_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/22_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/230_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/231_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/232_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/233_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/234_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/235_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/236_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/237_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/238_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/239_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/23_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/240_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/241_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/242_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/243_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/244_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/245_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/246_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/247_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/248_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/249_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/24_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/250_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/251_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/252_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/253_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/254_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/255_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/256_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/257_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/258_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/259_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/25_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/260_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/261_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/262_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/263_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/264_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/265_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/266_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/267_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/268_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/269_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/26_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/270_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/271_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/272_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/273_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/274_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/275_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/276_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/277_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/278_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/279_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/27_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/280_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/281_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/282_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/283_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/284_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/285_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/286_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/287_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/288_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/289_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/290_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/291_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/292_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/293_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/294_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/295_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/296_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/297_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/298_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/299_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/2_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/300_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/301_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/302_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/304_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/309_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/30_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/311_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/316_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/318_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/319_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/31_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/320_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/324_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/325_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/328_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/329_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/32_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/330_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/332_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/333_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/337_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/339_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/33_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/341_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/342_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/343_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/347_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/349_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/351_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/354_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/358_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/359_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/35_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/360_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/362_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/365_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/36_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/371_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/376_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/377_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/378_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/379_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/37_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/380_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/381_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/382_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/383_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/387_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/388_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/389_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/38_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/390_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/392_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/393_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/395_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/397_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/398_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/39_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/3_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/40_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/41_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/42_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/43_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/44_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/45_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/46_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/47_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/48_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/49_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/4_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/50_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/51_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/52_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/53_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/54_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/55_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/56_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/57_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/58_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/5_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/60_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/61_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/62_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/63_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/64_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/65_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/66_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/67_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/68_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/69_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/6_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/70_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/71_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/72_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/73_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/74_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/75_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/76_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/77_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/78_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/79_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/7_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/80_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/82_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/83_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/84_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/85_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/86_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/87_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/88_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/89_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/8_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/90_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/91_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/92_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/93_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/94_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/95_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/96_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/97_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/98_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/99_unweighted_events.lhe',
'file:/data4/syu/madgraph_lhe/lhe_files/9_unweighted_events.lhe',
	

	]
		  );

process.options = cms.untracked.PSet(
	)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
	version = cms.untracked.string('$Revision: 1.381.2.2 $'),
	annotation = cms.untracked.string('Configuration/GenProduction/python/EightTeV/file:/tmp/ymtzeng/ttbar.root'),
	name = cms.untracked.string('PyReleaseValidation')
	)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
					splitLevel = cms.untracked.int32(0),
					eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
					outputCommands = process.RAWSIMEventContent.outputCommands,
					fileName = cms.untracked.string('sim.root'),
					dataset = cms.untracked.PSet(
	filterName = cms.untracked.string(''),
	dataTier = cms.untracked.string('GEN-SIM')
	),
					SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
	)
					)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'START53_V7A::All'

process.generator = cms.EDFilter("Pythia6HadronizerFilter",
    jetMatching = cms.untracked.PSet(
        MEMAIN_showerkt = cms.double(0),
        MEMAIN_nqmatch = cms.int32(5),
        MEMAIN_minjets = cms.int32(0),
        MEMAIN_maxjets = cms.int32(2),
        MEMAIN_qcut = cms.double(20.0),
        MEMAIN_excres = cms.string(''),
        MEMAIN_etaclmax = cms.double(5.0),
        outTree_flag = cms.int32(0),
        scheme = cms.string('Madgraph'),
        mode = cms.string('auto')
    ),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    displayPythiaCards   = cms.untracked.bool(True),				    comEnergy = cms.double(7000.0),
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

	'MSTP(81)=21    ! multiple parton interactions 1 is Pythia default', 
	'MSTP(82)=4     ! Defines the multi-parton model', 
	),
	processParameters = cms.vstring(
            'MSTJ(1)=1       ! Fragmentation/hadronization on or off',
	    ),
        parameterSets = cms.vstring('pythiaUESettings', 'processParameters')                                )

)

baseJetSel = cms.PSet(
  Jets=cms.InputTag("selectedPatJetsPFlow")
)

from DelPanj.TreeMaker.eSel2012HZZ_cff import *
from DelPanj.TreeMaker.muSel2012HZZ_cff import *

process.tree = cms.EDAnalyzer(
	'TreeMaker',
	fillPUweightInfo_ = cms.bool(False),
	fillEventInfo_ = cms.bool(True),
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
	patMuons=cms.InputTag("userDataSelectedMuons"),
	patElectrons = cms.InputTag("userDataSelectedElectrons"),
        e2012IDSet  =  eSel2012HZZ,
        mu2012IDSet = muSel2012HZZ,
        eleRhoIso = cms.InputTag("kt6PFJetsForIso","rho"),
        muoRhoIso = cms.InputTag("kt6PFJetsCentralNeutral", "rho"),
	patMet=cms.InputTag("patMETs"),
	beamSpotLabel=cms.InputTag("offlineBeamSpot"),
	patJetPfAk05 = baseJetSel,
	outFileName=cms.string('outputFileName.root')
	)



process.TFileService = cms.Service("TFileService", fileName =
				   cms.string('/data4/syu/madgraph_lhe/zJettest_test.root'))

process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
#process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)
process.end_ana     = cms.Path(process.tree)


# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.end_ana)
#process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 


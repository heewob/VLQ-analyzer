import FWCore.ParameterSet.Config as cms
process = cms.Process("analysis")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '106X_upgrade2018_realistic_v16_L1v1'

process.skeletonAnalyzer = cms.EDAnalyzer("skeletonAnalyzer",
fatJetCollection = cms.InputTag("slimmedJetsAK8"),
   genPartCollection = cms.InputTag("prunedGenParticles"),
   jetCollection = cms.InputTag("slimmedJets"),
   bits = cms.InputTag("TriggerResults", "", "HLT"),
   muonCollection= cms.InputTag("slimmedMuons"),
   electronCollection = cms.InputTag("slimmedElectrons"),
tauCollection = cms.InputTag("slimmedTaus"),
 metCollection = cms.InputTag("slimmedMETs"),
 pileupCollection = cms.InputTag("slimmedAddPileupInfo")
)

process.leptonVeto = cms.EDFilter("leptonVeto",
 muonCollection= cms.InputTag("slimmedMuons"),
   electronCollection = cms.InputTag("slimmedElectrons"),
tauCollection = cms.InputTag("slimmedTaus"))
process.source = cms.Source("PoolSource",
fileNames = cms.untracked.vstring("/store/mc/RunIIAutumn18MiniAOD/TprimeTprime_M-1300_TuneCP5_PSweights_13TeV-madgraph-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/270000/054432FA-430C-2C4F-9E15-2019196FB8CC.root")
) # the input file here is stored at a data site. You can use the CMS Data Aggregation Service to find a dataset and get the storage location of the files you are interested in. You can also download the file locally to your src folder (using xrdcp, which you can learn about here https://uscms.org/uscms_at_work/computing/LPC/usingEOSAtLPC.shtml), but the above line needs to be formatted this way: cms.untracked.vstring("file:<file name>')  
process.TFileService = cms.Service("TFileService",fileName = cms.string("TprimeTprime_output_nolep.root")
)
process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(5000000)
)
process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(True),
)
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.p = cms.Path(process.leptonVeto* process.skeletonAnalyzer)

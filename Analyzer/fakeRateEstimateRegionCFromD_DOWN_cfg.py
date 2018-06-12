### PDG IDs ###

A_PDGID = 36
MU_PDGID = 13
TAU_PDGID = 15
ANY_PDGID = 0

### Tau decay types ###

TAU_HAD = 0
TAU_MU = 1
TAU_E = 2
TAU_ALL = 3

### Tau hadronic decay types ###

TAU_ALL_HAD = -1
TAU_1PRONG_0NEUTRAL = 0
TAU_1PRONG_1NEUTRAL = 1
TAU_1PRONG_2NEUTRAL = 2
TAU_1PRONG_3NEUTRAL = 3
TAU_1PRONG_NNEUTRAL = 4
TAU_2PRONG_0NEUTRAL = 5
TAU_2PRONG_1NEUTRAL = 6
TAU_2PRONG_2NEUTRAL = 7
TAU_2PRONG_3NEUTRAL = 8
TAU_2PRONG_NNEUTRAL = 9
TAU_3PRONG_0NEUTRAL = 10
TAU_3PRONG_1NEUTRAL = 11
TAU_3PRONG_2NEUTRAL = 12
TAU_3PRONG_3NEUTRAL = 13
TAU_3PRONG_NNEUTRAL = 14
TAU_RARE = 15

### No consideration of pT rank ###

ANY_PT_RANK = -1

#################
# Initialization
#################
import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
process = cms.Process("CleanJetsAnalyzer")

###################################################
# initialize MessageLogger and output report
###################################################
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
process.options   = cms.untracked.PSet( 
		wantSummary = cms.untracked.bool(True), 
		SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

####################
# Input File List
####################
process.source = cms.Source("PoolSource",
         fileNames = cms.untracked.vstring(
'root://eoscms/FILE_PATHRegionB_selection_NUM.root')
)

process.ggh = cms.EDAnalyzer("FakeRateMiniAODEstimateRegionA",
   outFileName = cms.string('/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/DIRNAME/DIRNAME_Plots_NUM.root'),
   tauTag = cms.InputTag("muHadTauDMIsoSelector"),
   checkBTag = cms.bool(True),
   jetTag = cms.InputTag("slimmedJets"),
   csvBTag = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
   mu3dRMin = cms.double(0.0),
   mu3dRMax = cms.double(0.8),
   tauPtCut = cms.double(10.0),
   diMudRCut = cms.double(1.0),
   tauHadOverlapdRCut = cms.double(.8),
   mu3dROverlapCut = cms.double(.4),
   mu3Tag = cms.InputTag('Mu3'),
   mu12Tag = cms.InputTag('Mu1Mu2'),
   requireRemovedMuon = cms.bool(True),
   rooDataset = cms.bool(True),
#   CanvasName = cms.string("FinalFakeRateDMtoMedIsoOnlyEtavsPtCanvas"),
   HistName = cms.string("FinalFakeRateDMtoMedIsoOnlyEtavsPt"),
   TH2FileName = cms.string('/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/TAUFAKERATES_DOWN.root')
)

process.p2 = cms.Path(
	process.ggh
)

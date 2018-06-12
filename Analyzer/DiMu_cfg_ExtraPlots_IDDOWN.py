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

process.ggh = cms.EDAnalyzer("DiMu_ExtraPlots",
   outFileName = cms.string('/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/DIRNAME/DIRNAME_Plots_NUM.root'),
   mu12Tag = cms.InputTag('Mu1Mu2'),
   tauTag = cms.InputTag('muHadTauDMIsoSelector'),
   checkBTag = cms.bool(True),
   csvBTag = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
   bTagSFShift = cms.string("central"),   
   mu3Tag = cms.InputTag('Mu3'),
   metTag = cms.InputTag('slimmedMETs'),
   jetTag = cms.InputTag("slimmedJets"),
   pileupSummaryInfo = cms.InputTag("slimmedAddPileupInfo", "", "PAT"),
   genEventInfoToken = cms.InputTag("generator", "", "SIM"),
   tauPtCut = cms.double(10.0),
   diMudRCut = cms.double(1.0),
   tauHadOverlapdRCut = cms.double(.8),
   mu3dROverlapCut = cms.double(.4),
   xsec = cms.double(XSEC),
   lumi = cms.double(LUMI_DATA),
   summedWeights = cms.double(SUMMED_WEIGHTS),
   MC = cms.bool(True),
   ApplyFR = cms.bool(False),
   rooDataset = cms.bool(True),
   PileupFileName = cms.string('/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/PileupWeights_69200.root'),
   fpIDs_BToF = cms.string("/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/Zmumu_RunBCDEF_ID_DOWN.root"),
   fpIDs_GH = cms.string("/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/Zmumu_RunGH_ID_DOWN.root"),
   fpISOs_BToF = cms.string("/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/Zmumu_RunBCDEF_ISO.root"),
   fpISOs_GH = cms.string("/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/Zmumu_RunGH_ISO.root"),
   fpTrack = cms.string("/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/Tracking_EfficienciesAndSF_BCDEFGH.root"),
   fpTrigger_BToF = cms.string("/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/EfficienciesAndSF_RunBtoF.root"),
   fpTrigger_GH = cms.string("/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/EfficienciesAndSF_Period4.root"),
   fp_LowPt = cms.string("/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/TNP_LowMuPt_EFFICIENCIES_IDDOWN.root")
)

#########################################################
# this will produce a ref to the original muon collection
#########################################################
process.p2 = cms.Path(
	process.ggh
)

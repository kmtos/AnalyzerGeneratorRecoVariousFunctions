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
'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYHighMass_MedIsoMu2_TauDMAntiMedIso_JAN9/180111_230446/0000/RegionB_selection_101.root')
)

process.ggh = cms.EDAnalyzer("DiMu_ExtraPlots",
   outFileName = cms.string('/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/IRNAME_Plots_NUM.root'),
   mu12Tag = cms.InputTag('Mu1Mu2'),
   tauTag = cms.InputTag('muHadTauDMIsoSelector'),
   mu3Tag = cms.InputTag('Mu3ID'),
   metTag = cms.InputTag('slimmedMETs'),
   jetTag = cms.InputTag("slimmedJets"),
   pileupSummaryInfo = cms.InputTag("slimmedAddPileupInfo", "", "PAT"),
   genEventInfoToken = cms.InputTag("generator", "", "SIM"),
   tauPtCut = cms.double(20.0),
   xsec = cms.double(1),
   lumi = cms.double(1),
   summedWeights = cms.double(1),
   PileupFileName = cms.string('/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/PileupWeights.root')    
)

#########################################################
# this will produce a ref to the original muon collection
#########################################################
process.p2 = cms.Path(
	process.ggh
)
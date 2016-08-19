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
mylist = FileUtils.loadListFromFile('/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_6/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/inFileList_ggH300_a9_CleanJets.txt')

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
readFiles = cms.untracked.vstring(*mylist)
process.source = cms.Source("PoolSource",
    fileNames = readFiles,
    skipEvents = cms.untracked.uint32(0)
    )

process.ggh = cms.EDAnalyzer("GGHAnalyzer_OLD",
   outFileName = cms.string('/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_6/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/DIRNAME/DIRNAME_Plots.root'),
   genParticleTag = cms.InputTag("genParticles", "", ""),
   akJetTag = cms.InputTag("ak4PFJets"),
   muonsTag = cms.InputTag("muons"),
   vtxTag = cms.InputTag('offlinePrimaryVertices'),
   muonMapTag = cms.InputTag("CleanJets", "muonValMap"),
   jetValMapTag = cms.InputTag("CleanJets", "jetValMap", "CLEANJETS"),
   tauRECOTag = cms.InputTag("hpsPFTauProducer", "", "RECO"),
   tauCJTag = cms.InputTag("hpsPFTauProducer", "", "CLEANJETS"),
   pizerosTag = cms.InputTag("hpsPFTauProducer", "pizeros" ),
   looseIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits", "", "CLEANJETS"),
   medIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits", "", "CLEANJETS"),
   tightIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits", "", "CLEANJETS"),
   decayModeFindingTagCJ = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingOldDMs", "", "CLEANJETS"),
   looseIsoTagRECO = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits", "", "RECO"),
   medIsoTagRECO = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits", "", "RECO"),
   tightIsoTagRECO = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits", "", "RECO"),
   decayModeFindingTagRECO = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingOldDMs", "", "RECO"),
   looseIsoTagCJMVA = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVA3oldDMwLT", "", "CLEANJETS"),
   medIsoTagCJMVA = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVA3oldDMwLT", "", "CLEANJETS"),
   tightIsoTagCJMVA = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVA3oldDMwLT", "", "CLEANJETS"),
   looseIsoTagRECOMVA = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVA3oldDMwLT", "", "RECO"),
   medIsoTagRECOMVA = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVA3oldDMwLT", "", "RECO"),
   tightIsoTagRECOMVA = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVA3oldDMwLT", "", "RECO")

)


process.p2 = cms.Path(
	process.ggh
)

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
<<<<<<< HEAD
mylist = FileUtils.loadListFromFile('/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_6/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/inFileList_TT_CleanJets.txt')
=======
mylist = FileUtils.loadListFromFile('/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/inFileList_TT_CleanJets.txt')
>>>>>>> 3ec953b4f9b58e1126afc333ff858db8e9e35a36

process = cms.Process("FAKERATEANALYZER")

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

#########################################
# Rerunning bTaggin on CleanJets Sample
#########################################
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")  #Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff") #Configuration.StandardSequences.FrontierConditions_CMS.GlobalTag_cff")
<<<<<<< HEAD
process.GlobalTag.globaltag = cms.string('80X_mcRun2_asymptotic_v14') #'IDEAL_V9::All'
=======
process.GlobalTag.globaltag = cms.string('76X_dataRun2_v15') #'IDEAL_V9::All'
>>>>>>> 3ec953b4f9b58e1126afc333ff858db8e9e35a36
process.impactParameterTagInfos.jetTracks = cms.InputTag("ak4JetTracksAssociatorAtVertex")
process.ak4JetTracksAssociatorAtVertex.jets = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')
process.ak4JetTracksAssociatorAtVertex.tracks = cms.InputTag("generalTracks")

process.btagging = cms.Sequence(
    process.ak4JetTracksAssociatorAtVertex*
    # impact parameters and IP-only algorithms
    process.impactParameterTagInfos*
    (process.trackCountingHighEffBJetTags +
     process.trackCountingHighPurBJetTags +
     process.jetProbabilityBJetTags +
     process.jetBProbabilityBJetTags +
     # SV tag infos depending on IP tag infos, and SV (+IP) based algos
     process.secondaryVertexTagInfos*
     (process.simpleSecondaryVertexHighEffBJetTags +
      process.simpleSecondaryVertexHighPurBJetTags +
      process.combinedSecondaryVertexBJetTags) +
     process.ghostTrackVertexTagInfos*
     process.ghostTrackBJetTags)
)

####################
# Input File List
####################
readFiles = cms.untracked.vstring(*mylist)
process.source = cms.Source("PoolSource",
    fileNames = readFiles,
    skipEvents = cms.untracked.uint32(0)
    )

process.ggh = cms.EDAnalyzer("FakeRateAnalyzer",
<<<<<<< HEAD
   outFileName = cms.string('/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_6/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/DIRNAME/DIRNAME_Plots.root'),
=======
   outFileName = cms.string('/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/DIRNAME/DIRNAME_Plots.root'),
>>>>>>> 3ec953b4f9b58e1126afc333ff858db8e9e35a36
   genParticleTag = cms.InputTag("genParticles", "", ""),
   genJetTag = cms.InputTag("ak4GenJets", "", ""),
   akJetTag = cms.InputTag("ak4PFJets"),
   muonsTag = cms.InputTag("muons"),
  vtxTag = cms.InputTag('offlinePrimaryVertices'),
   muonMapTag = cms.InputTag("CleanJets", "muonValMap"),
   jetValMapTag = cms.InputTag("CleanJets", "jetValMap", "CLEANJETS"),
   tauCJTag = cms.InputTag("hpsPFTauProducer", "", "CLEANJETS"),
   pizerosTag = cms.InputTag("hpsPFTauProducer", "pizeros" ),
<<<<<<< HEAD
   #looseIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVA3oldDMwLT", "", "CLEANJETS"),
   looseIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits", "", "CLEANJETS"),
   #medIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVA3oldDMwLT", "", "CLEANJETS"),
   medIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits", "", "CLEANJETS"),
   #tightIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVA3oldDMwLT", "", "CLEANJETS"),
   tightIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits", "", "CLEANJETS"),
   decayModeFindingTagCJ = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs", "", "CLEANJETS"),
   genMatchPDGIDTag = cms.int32(6),
   oldJetTag = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS'),
   csvBTag = cms.InputTag("combinedSecondaryVertexBJetTags", "", "FAKERATEANALYZER")
=======
   looseIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVA3oldDMwLT", "", "CLEANJETS"),
   #looseIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits", "", "CLEANJETS"),
   medIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVA3oldDMwLT", "", "CLEANJETS"),
   #medIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits", "", "CLEANJETS"),
   tightIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVA3oldDMwLT", "", "CLEANJETS"),
   #tightIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits", "", "CLEANJETS"),
   decayModeFindingTagCJ = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingOldDMs", "", "CLEANJETS"),
   genMatchPDGIDTag = cms.int32(6),
   oldJetTag = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS'),
   csvBTag = cms.InputTag("pfCombinedMVAV2BJetTags", "", "FAKERATEANALYZER") #combinedSecondaryVertexBJetTags 
>>>>>>> 3ec953b4f9b58e1126afc333ff858db8e9e35a36
)


process.p2 = cms.Path(
	process.btagging*
	process.ggh
)

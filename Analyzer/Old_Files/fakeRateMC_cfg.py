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

process.ggh = cms.EDAnalyzer("FakeRateMCAnalyzer",
   outFileName = cms.string('/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/DIRNAME/DIRNAME_Plots_NUM.root'),
   akJetTag = cms.InputTag("ak4PFJets"),
   tauTag = cms.InputTag("hpsPFTauProducer", "", "SKIM"),
   #looseIsoTag = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits", "", "SKIM"),
   looseIsoTag = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT", "", "SKIM"),
   #medIsoTag = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits", "", "SKIM"),
   medIsoTag = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT", "", "SKIM"),
   #tightIsoTag = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits", "", "SKIM"),
   tightIsoTag = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT", "", "SKIM"),
   decayModeFindingTag = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs", "", "SKIM"),
   isoRawTag = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits", "", "SKIM"),
   oldJetTag = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'SKIM'),
   csvBTag = cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags", "", "SKIM"),
   mu12Tag = cms.InputTag('Isolate'),
   muonsTag = cms.InputTag("muons"),
   muonMapTag = cms.InputTag("CleanJets", "muonValMap"),
   requireRemovedMuon = cms.bool(True),
   checkInvMass = cms.bool(True),
   checkInvMassValue = cms.double(4.0),
   pileupSummaryInfo = cms.InputTag("addPileupInfo", "", "HLT"),
   PileupFileName = cms.string('/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/PileupWeights.root'),
   xsec = cms.double(XSEC),
   lumiData = cms.double(LUMI_DATA),
   summedWeights = cms.double(SUMMED_WEIGHTS),
   genEventInfoToken = cms.InputTag("generator", "", "SIM")

)

#########################################################
# this will produce a ref to the original muon collection
#########################################################
process.muonsRef = cms.EDFilter('MuonRefSelector',
                   src = cms.InputTag('muons'),
                   cut = cms.string('pt > 0.0'),
                   filter = cms.bool(True)
)

process.p2 = cms.Path(
        process.muonsRef*
	process.ggh
)

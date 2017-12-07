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
#mylist = FileUtils.loadListFromFile('/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/inFileList_ggH125_a9_Trig_CJ_NewMVA.txt')
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
#        fileNames = cms.untracked.vstring('root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-19_TuneCUETP8M1_13TeV_madgraph_pythia8/CRAB3_H125_a19_Mu2Loose_DEC2/161213_165941/0000/Background_edm_output_NUM.root')
#        fileNames = cms.untracked.vstring('root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-750_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/CRAB3_H750_a9_Mu2Loose_DEC2/161213_165902/0000/Background_edm_output_NUM.root')
#       fileNames = cms.untracked.vstring('root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-5_TuneCUETP8M1_13TeV_madgraph_pythia8/CRAB3_H125_a5_Mu2Loose_DEC2/161213_165854/0000/Background_edm_output_NUM.root')
#        fileNames = cms.untracked.vstring('root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/CRAB3_H125_a9_Mu2Loose_DEC2/161213_165839/0000/Background_edm_output_NUM.root')

## These are just for the Tau CJ efficiencies
#        fileNames = cms.untracked.vstring('root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-19_TuneCUETP8M1_13TeV_madgraph_pythia8/CRAB3_H125_a19_CJOnly_DEC2/170106_044801/0000/Background_edm_output_NUM.root')
        fileNames = cms.untracked.vstring('root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-750_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/CRAB3_H750_a9_CJOnly_DEC2/170106_174616/0000/Background_edm_output_NUM.root')
#       fileNames = cms.untracked.vstring('root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-5_TuneCUETP8M1_13TeV_madgraph_pythia8/CRAB3_H125_a5_CJOnly_DEC2/170106_174616/0000/Background_edm_output_NUM.root')
#        fileNames = cms.untracked.vstring('root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/CRAB3_H125_a9_CJOnly_DEC2/170106_053407/0000/Background_edm_output_NUM.root')
)

#readFiles = cms.untracked.vstring(*mylist)
#process.source = cms.Source("PoolSource",
#    fileNames = readFiles,
#    skipEvents = cms.untracked.uint32(0)
#    )

process.ggh = cms.EDAnalyzer("GGHAnalyzer_OLD",
   outFileName = cms.string('/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/DIRNAME/DIRNAME_Plots_NUM.root'),
   genParticleTag = cms.InputTag("genParticles", "", ""),
   akJetTag = cms.InputTag("ak4PFJets"),
   muonsTag = cms.InputTag("muons"),
   vtxTag = cms.InputTag('offlinePrimaryVertices'),
   muonMapTag = cms.InputTag("CleanJets", "muonValMap"),
   jetValMapTag = cms.InputTag("CleanJets", "jetValMap", "SKIM"),
   tauRECOTag = cms.InputTag("hpsPFTauProducer", "", "RECO"),
   tauCJTag = cms.InputTag("hpsPFTauProducer", "", "SKIM"),
   pizerosTag = cms.InputTag("hpsPFTauProducer", "pizeros" ),
   #looseIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits", "", "SKIM"),
   looseIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT", "", "SKIM"),
   #medIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits", "", "SKIM"),
   medIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT", "", "SKIM"),
   #tightIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits", "", "SKIM"),
   tightIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT", "", "SKIM"),
   decayModeFindingTagCJ = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs", "", "SKIM"),
   #looseIsoTagRECO = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits", "", "RECO"),
   looseIsoTagRECO = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT", "", "RECO"),   
   #medIsoTagRECO = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits", "", "RECO"),
   medIsoTagRECO = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT", "", "RECO"),   
   #tightIsoTagRECO = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits", "", "RECO"),
   tightIsoTagRECO = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT", "", "RECO"),
   decayModeFindingTagRECO = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs", "", "RECO"),
   isoRawTag = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits", "", "SKIM"),
   oldJetTag = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'SKIM'),
   mu12Tag = cms.InputTag('Isolation'),
   mu3Tag = cms.InputTag('Mu3ID'),
   csvBTag = cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags", "", "SKIM")
)


process.p2 = cms.Path(
	process.ggh
)

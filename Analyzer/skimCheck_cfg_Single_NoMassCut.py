### PTTG IDs ###

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

####################
# Input File List
####################
process.source = cms.Source("PoolSource",
         fileNames = cms.untracked.vstring(
####   SeparatedDiMu ####
#'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYHighMass_SeparatedDiMu_NoMassCut_FEB9/170227_192345/0000/RegionB_selection_NUM.root')
#'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/DY1JetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYLowMass_SeparatedDiMu_NoMassCut_FEB9/170227_192357/0000/RegionB_selection_NUM.root')
#'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/QCD_1000toInf_SeparatedDiMu_NoMassCut_FEB9/170227_192409/0000/RegionB_selection_NUM.root')
#'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/QCD_120to170_SeparatedDiMu_NoMassCut_FEB9/170227_192422/0000/RegionB_selection_NUM.root')
#'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/QCD_15to20_SeparatedDiMu_NoMassCut_FEB9/170227_192435/0000/RegionB_selection_NUM.root')
#'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/QCD_170to300_SeparatedDiMu_NoMassCut_FEB9/170227_192452/0000/RegionB_selection_NUM.root')
#'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/QCD_20to30_SeparatedDiMu_NoMassCut_FEB9/170227_192504/0000/RegionB_selection_NUM.root')
#'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/QCD_300to470_SeparatedDiMu_NoMassCut_FEB9/170227_192516/0000/RegionB_selection_NUM.root')
#'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/QCD_30to50_SeparatedDiMu_NoMassCut_FEB9/170227_192529/0000/RegionB_selection_NUM.root')
#'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/QCD_470to600_SeparatedDiMu_NoMassCut_FEB9/170228_165746/0000/RegionB_selection_NUM.root')
#'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/QCD_50to80_SeparatedDiMu_NoMassCut_FEB9/170227_192557/0000/RegionB_selection_NUM.root')
#'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/QCD_600to800_SeparatedDiMu_NoMassCut_FEB9/170227_192610/0000/RegionB_selection_NUM.root')
#'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/QCD_800to1000_SeparatedDiMu_NoMassCut_FEB9/170227_192623/0000/RegionB_selection_NUM.root')
#'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/QCD_80to120_SeparatedDiMu_NoMassCut_FEB9/170227_192636/0000/RegionB_selection_NUM.root')
#'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-19_TuneCUETP8M1_13TeV_madgraph_pythia8/SignalH125a19_SeparatedDiMu_NoMassCut_FEB9/170227_192044/0000/RegionB_selection_NUM.root')
#'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-5_TuneCUETP8M1_13TeV_madgraph_pythia8/SignalH125a5_SeparatedDiMu_NoMassCut_FEB9/170227_192056/0000/RegionB_selection_NUM.root')
#'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/SignalH125a9_SeparatedDiMu_NoMassCut_FEB9/170227_192109/0000/RegionB_selection_NUM.root')
#'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-750_M-9_TuneCUETP8M1_13TeV_madgraph_pythia8/SignalH750a9_SeparatedDiMu_NoMassCut_FEB9/170227_192121/0000/RegionB_selection_NUM.root')
'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TTBar_SeparatedDiMu_NoMassCut_FEB9/170227_192649/0000/RegionB_selection_NUM.root')


)

process.ggh = cms.EDAnalyzer("SkimCheck",
   outFileName = cms.string('/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/DIRNAME/DIRNAME_NUM.root'),
   akJetTag = cms.InputTag("ak4PFJets"),
   muonsTag = cms.InputTag("muons"),
   muonMapTag = cms.InputTag("CleanJets", "muonValMap"),
   jetValMapTag = cms.InputTag("CleanJets", "jetValMap", "SKIM"),
   tauTag = cms.InputTag("muHadIsoTauSelector", "", "SKIM"),
   tightIsoTag = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT", "", "SKIM"),
   veryTightIsoTag = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT", "", "SKIM"),
   decayModeFindingTag = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs", "", "SKIM"),
   isoRawTag = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits", "", "SKIM"),
   oldJetTag = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'SKIM'),
   mu12Tag = cms.InputTag('Isolate'),
   mu3Tag = cms.InputTag('Mu3ID'),
   csvBTag = cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags", "", "SKIM") 
)


process.p2 = cms.Path(
	process.ggh
)

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
'root://eoscms//eos/cms/store/group/phys_higgs/HiggsExo/ktos/TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/TTBar_SignalRegion_SEP4/170906_205116/0000/RegionB_selection_50.root')
)

process.ggh = cms.EDAnalyzer("FakeRateMCAnalyzer",
   outFileName = cms.string('/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/RNAME_Plots_NUM.root'),
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
   dataPileupFileName = cms.string('/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/pileupHistogram_AllDataNorm.root'),
   probFunctionVariable = cms.vint32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74),
   probValue = cms.vdouble(1.78653e-05 ,2.56602e-05 ,5.27857e-05 ,8.88954e-05 ,0.000109362 ,0.000140973 ,0.000240998 ,0.00071209 ,0.00130121 ,0.00245255 ,0.00502589 ,0.00919534 ,0.0146697 ,0.0204126 ,0.0267586 ,0.0337697 ,0.0401478 ,0.0450159 ,0.0490577 ,0.0524855 ,0.0548159 ,0.0559937 ,0.0554468 ,0.0537687 ,0.0512055 ,0.0476713 ,0.0435312 ,0.0393107 ,0.0349812 ,0.0307413 ,0.0272425 ,0.0237115 ,0.0208329 ,0.0182459 ,0.0160712 ,0.0142498 ,0.012804 ,0.011571 ,0.010547 ,0.00959489 ,0.00891718 ,0.00829292 ,0.0076195 ,0.0069806 ,0.0062025 ,0.00546581 ,0.00484127 ,0.00407168 ,0.00337681 ,0.00269893 ,0.00212473 ,0.00160208 ,0.00117884 ,0.000859662 ,0.000569085 ,0.000365431 ,0.000243565 ,0.00015688 ,9.88128e-05 ,6.53783e-05 ,3.73924e-05 ,2.61382e-05 ,2.0307e-05 ,1.73032e-05 ,1.435e-05 ,1.36486e-05 ,1.35555e-05 ,1.37491e-05 ,1.34255e-05 ,1.33987e-05 ,1.34061e-05 ,1.34211e-05 ,1.34177e-05 ,1.32959e-05 ,1.33287e-05),
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

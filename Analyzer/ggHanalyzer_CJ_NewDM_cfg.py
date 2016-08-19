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
mylist = FileUtils.loadListFromFile('/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_6/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/inFileList_ggH125_a9_CleanJets.txt')

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

#########################################
# Rerunning bTaggin on CleanJets Sample
#########################################
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")  #Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff") #Configuration.StandardSequences.FrontierConditions_CMS.GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('76X_dataRun2_v15') #'IDEAL_V9::All'
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

process.ggh = cms.EDAnalyzer("GGHAnalyzer_IndivCJ",
   outFileName = cms.string('/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_6/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/DIRNAME/DIRNAME_Plots.root'),
   genParticleTag = cms.InputTag("genParticles", "", ""),
   akJetTag = cms.InputTag("ak4PFJets"),
   muonsTag = cms.InputTag("muons"),
   muonMapTag = cms.InputTag("CleanJets", "muonValMap"),
   jetValMapTag = cms.InputTag("CleanJets", "jetValMap", "CLEANJETS"),
   tauCJTag = cms.InputTag("hpsPFTauProducer", "", "CLEANJETS"),
   pizerosTag = cms.InputTag("hpsPFTauProducer", "pizeros" ),
   looseIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits", "", "CLEANJETS"),
   medIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits", "", "CLEANJETS"),
   tightIsoTagCJ = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits", "", "CLEANJETS"),
   decayModeFindingTagCJ = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs", "", "CLEANJETS"),
   isoPtSumTagCJ = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr", "", "CLEANJETS"),
   genMatchedTauVisiblePtMapTagCJ = cms.InputTag("genATauMuMatchedRecoTauSelectorCJ", "TauVisiblePtMap", "CleanJetsAnalyzer"),
   genMatchedTauDecayModeMapTagCJ = cms.InputTag("genATauMuMatchedRecoTauSelectorCJ", "TauDecayModeMap", "CleanJetsAnalyzer"),
   genMatchedTauMatchedMapTagCJ = cms.InputTag("genATauMuMatchedRecoTauSelectorCJ", "TauMatchedMap", "CleanJetsAnalyzer"),
   genMatchedRecoTausCJ = cms.InputTag("genATauMuMatchedRecoTauSelectorCJ", "valMapAccessers", "CleanJetsAnalyzer"),
   looseIsoTagCJMVA = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVA3newDMwLT", "", "CLEANJETS"),
   medIsoTagCJMVA = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVA3newDMwLT", "", "CLEANJETS"),
   tightIsoTagCJMVA = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVA3newDMwLT", "", "CLEANJETS"),
   oldJetTag = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS'),
   csvBTag = cms.InputTag("combinedSecondaryVertexBJetTags", "", "FAKERATEANALYZER")

)

###########################
# GenTauDecayID definition
###########################
ATauTauPSet = cms.PSet(momPDGID = cms.vint32(A_PDGID),
                       chargedHadronPTMin = cms.double(0.0), #should always be 0.0
                       neutralHadronPTMin = cms.double(0.0), #should always be 0.0
                       chargedLeptonPTMin = cms.double(0.0), #should always be 0.0
                       totalPTMin = cms.double(0.0)) #should always be 0.0

process.genATauHadSelectorCJ = cms.EDFilter(
    'GenObjectProducer',
    genParticleTag = cms.InputTag('genParticles'),
    absMatchPDGIDs = cms.vuint32(TAU_PDGID),     #choose a gen tau...
    sisterAbsMatchPDGID = cms.uint32(TAU_PDGID), #...whose sister is another gen tau...
    genTauDecayIDPSet = ATauTauPSet,             #...and whose mother is a pseudoscalar a
    primaryTauDecayType = cms.uint32(TAU_HAD),    #primary tau decay mode is mu...
    sisterTauDecayType = cms.uint32(TAU_MU),    #...sister tau decay mode is hadronic
    primaryTauPTRank = cms.int32(ANY_PT_RANK),  #should always be ANY_PT_RANK
    primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD), #choose TAU_ALL_HAD when the tau decay type is non-hadronic
    sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),     #choose TAU_ALL_HAD when the tau decay type is hadronic and you want any hadronic mode
    primaryTauAbsEtaMax = cms.double(-1),  #|eta| < 2.1 on muon from tau-->mu
    primaryTauPTMin = cms.double(-1),      #pT > 5 GeV on muon from tau-->mu
    countSister = cms.bool(False),          #only put the muon from tau-->mu in the output collection (i.e. object has  |PDG ID| = 13 and status = 1 that is decayed from the tau)
    applyPTCuts = cms.bool(False),          #should always be False
    countKShort = cms.bool(False),          #should always be False
    minNumGenObjectsToPassFilter = cms.uint32(1), #EDFilter only returns true if >=1 tau_Had is found satisfying pT, |eta|, and decay mode cuts
    makeAllCollections = cms.bool(False) #should always be False
    )

process.recoTauSelectorCJ = cms.EDFilter('PFTauRefSelector',
                                        src = cms.InputTag("hpsPFTauProducer", "", "CLEANJETS"),
                                        cut = cms.string('pt > 0.0'),
                                        filter = cms.bool(True)
                                        )

process.genATauMuMatchedRecoTauSelectorCJ = cms.EDFilter(
    'TauMatchedRecoObjectProducer',
    genParticleTag = cms.InputTag('genParticles'),
    selectedGenParticleTag = cms.InputTag('genATauHadSelectorCJ'), #must be a reco::GenParticleRefVector
    recoObjTag = cms.InputTag('recoTauSelectorCJ'),      
    baseRecoObjTag = cms.InputTag("hpsPFTauProducer", "", "CLEANJETS"),
    genTauDecayIDPSet = ATauTauPSet,      #need to know the pseudoscalar a mother
    applyPTCuts = cms.bool(False),        #should always be false
    countKShort = cms.bool(False),        #should always be false
    pTRank = cms.int32(ANY_PT_RANK),      #should always be ANY_PT_RANK
    makeAllCollections = cms.bool(False), #should always be False
    useGenObjPTRank = cms.bool(True),     #should always be True
    nOutputColls = cms.uint32(1),         #should always be 1
    dR = cms.double(0.1),                 #dR criteria for matching
    ifTauColl = cms.bool(True),		  #This Creates a map of GenTaus and their Decay Mode and Visible Pt
    minNumGenObjectsToPassFilter = cms.uint32(1) #EDFilter returns true if >=1 gen-matched reco muon is found
    )

process.p2 = cms.Path(
	process.genATauHadSelectorCJ*
	process.recoTauSelectorCJ*
	process.genATauMuMatchedRecoTauSelectorCJ*	
	process.btagging*
	process.ggh
)

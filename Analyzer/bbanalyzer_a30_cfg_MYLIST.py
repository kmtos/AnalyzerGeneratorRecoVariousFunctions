import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
mylist = FileUtils.loadListFromFile('/afs/cern.ch/user/k/ktos/BBA/CMSSW_7_4_1_patch1/src/BBA/Analyzer/inFileList_a30.txt')
process = cms.Process("BBA")

#process.load("BBA/Analyzer/bbaanalyzer_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_71_V1::All', '')
process.load("Configuration.StandardSequences.MagneticField_cff")

####################
# Message Logger
####################
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
## switch to uncheduled mode
process.options.allowUnscheduled = cms.untracked.bool(True)

####################
# Input File List
####################
readFiles = cms.untracked.vstring(*mylist)
process.source = cms.Source("PoolSource",
    fileNames = readFiles,
    skipEvents = cms.untracked.uint32(0)
    )

############################################################
# Defining matching in DeltaR, sorting by best DeltaR
############################################################
process.bTriggerMatchHLTMu10CentralPFJet30BTagCSV = cms.EDProducer("PATTriggerMatcherDRLessByR",
        src     = cms.InputTag( 'slimmedJets' ), # matcher input collections
        matched = cms.InputTag( 'patTrigger' ),  # selections of trigger objects
        matchedCuts = cms.string( 'type( "TriggerBJet" ) && path( "HLT_Mu10_CentralPFJet30_BTagCSV0p*")' ), # input does not yet have the 'saveTags' parameter in HLT
        maxDPtRel   = cms.double( 0.5 ), # no effect here
        maxDeltaR   = cms.double( 0.5 ), ####### selection of matches
        maxDeltaEta = cms.double( 0.2 ), # no effect here
        resolveAmbiguities    = cms.bool( True ),  # definition of matcher output
        resolveByMatchQuality = cms.bool( True )  # definition of matcher output
)       

process.mTriggerMatchHLTMu10CentralPFJet30BTagCSV = cms.EDProducer("PATTriggerMatcherDRLessByR",
        src     = cms.InputTag( 'slimmedMuons' ),
        matched = cms.InputTag( 'patTrigger' ), # selections of trigger objects
        matchedCuts = cms.string( 'type( "TriggerMuon" ) && path( "HLT_Mu10_CentralPFJet30_BTagCSV0p*")' ), # input does not yet have the 'saveTags' parameter in HLT
        maxDPtRel   = cms.double( 0.5 ), # no effect here
        maxDeltaR   = cms.double( 0.5 ), #### selection of matches
        maxDeltaEta = cms.double( 0.2 ), # no effect here
        resolveAmbiguities    = cms.bool( True ),# definition of matcher output
        resolveByMatchQuality = cms.bool( True )# definition of matcher output
)       

process.bTriggerMatchHLTIsoMu24eta2p1CentralPFJet30BTagCSV = cms.EDProducer("PATTriggerMatcherDRLessByR",
        src     = cms.InputTag( 'slimmedJets' ), # matcher input collections
        matched = cms.InputTag( 'patTrigger' ),  # selections of trigger objects
        matchedCuts = cms.string( 'type( "TriggerBJet" ) && path( "HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV*")' ), # input does not yet have the 'saveTags' parameter in HLT
        maxDPtRel   = cms.double( 0.5 ), # no effect here
        maxDeltaR   = cms.double( 0.5 ), ####### selection of matches
        maxDeltaEta = cms.double( 0.2 ), # no effect here
        resolveAmbiguities    = cms.bool( True ),  # definition of matcher output
        resolveByMatchQuality = cms.bool( True )  # definition of matcher output
)

process.mTriggerMatchHLTIsoMu24eta2p1CentralPFJet30BTagCSV = cms.EDProducer("PATTriggerMatcherDRLessByR",
        src     = cms.InputTag( 'slimmedMuons' ),
        matched = cms.InputTag( 'patTrigger' ), # selections of trigger objects
        matchedCuts = cms.string( 'type( "TriggerMuon" ) && path( "HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV*")' ), # input does not yet have the 'saveTags' parameter in HLT
        maxDPtRel   = cms.double( 0.5 ), # no effect here
        maxDeltaR   = cms.double( 0.5 ), #### selection of matches
        maxDeltaEta = cms.double( 0.2 ), # no effect here
        resolveAmbiguities    = cms.bool( True ),# definition of matcher output
        resolveByMatchQuality = cms.bool( True )# definition of matcher output
)

# talk to output module
process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string("/afs/cern.ch/user/k/ktos/BBA/CMSSW_7_4_1_patch1/src/BBA/Analyzer/BSUB/DIRNAME/edmTriggerOutput.root"),
	outputCommands = process.MINIAODSIMEventContent.outputCommands
)


################################################################################
# Running the matching and setting the the trigger on
################################################################################
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process ) # This is optional and can be omitted.
switchOnTriggerMatching( process, triggerMatchers = [ 'mTriggerMatchHLTIsoMu24eta2p1CentralPFJet30BTagCSV',
						      'bTriggerMatchHLTIsoMu24eta2p1CentralPFJet30BTagCSV',
						      'mTriggerMatchHLTMu10CentralPFJet30BTagCSV',
						      'bTriggerMatchHLTMu10CentralPFJet30BTagCSV'
						    ])
################################################################################
# Analyzer setting
################################################################################
process.bba = cms.EDAnalyzer("BBAAnalyzer",
   outFileName = cms.string('/afs/cern.ch/user/k/ktos/BBA/CMSSW_7_4_1_patch1/src/BBA/Analyzer/BSUB/DIRNAME/DIRNAME.root'),
   prunedGenParticleTag = cms.InputTag("prunedGenParticles"),
   slimmedJetTag      = cms.InputTag("slimmedJets"),
   slimmedTauTag      = cms.InputTag("slimmedTaus"),
   slimmedMuonTag     = cms.InputTag("slimmedMuons"),
   hpsTauTag	      = cms.InputTag("hpsPFTauProducer","","RECO"),
   muonMatchIsoMu     = cms.string( 'mTriggerMatchHLTIsoMu24eta2p1CentralPFJet30BTagCSV' ),
   muonMatchNonIsoMu  = cms.string( 'mTriggerMatchHLTMu10CentralPFJet30BTagCSV' ),
   bMatchNonIsoMu     = cms.string( 'bTriggerMatchHLTMu10CentralPFJet30BTagCSV' ),
   bMatchIsoMu        = cms.string( 'bTriggerMatchHLTIsoMu24eta2p1CentralPFJet30BTagCSV' ),
   trigger            = cms.InputTag( "patTrigger" ),
   triggerEvent       = cms.InputTag( "patTriggerEvent" )
)

process.p = cms.Path(process.bba)

// -*- C++ -*-
//
// Package:    FakeRateAnalyzer
// Class:      FakeRateAnalyzer
// 
/**\class FakeRateAnalyzer FakeRateAnalyzer.cc Analyzer/src/FakeRateAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Kyle Martin Tos
//         Created:  Thu Feb 26 09:51:23 CET 2015
// $Id$
//
//


// system include files
#include <memory>
#include <string>
#include <sstream>
#include <typeinfo>

// user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/JetFloatAssociation.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "HLTrigger/HLTcore/interface/HLTConfigData.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTanalyzers/interface/HLTInfo.h"

#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "AnalyzerGeneratorRecoVariousFunctions/VariousFunctions/interface/VariousFunctions.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"


using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;


//
// class declaration
//

class FakeRateAnalyzer : public edm::EDAnalyzer {
   public:
      typedef reco::JetFloatAssociation::Container JetBCEnergyRatioCollection;
      explicit FakeRateAnalyzer(const edm::ParameterSet&);
      ~FakeRateAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      //delete memory
      void reset(const bool);
    
      // ----------member data ---------------------------
      //pointer to output file object
      TFile* out_;

      //name of output root file
      std::string outFileName_;
      edm::EDGetTokenT<vector<reco::GenParticle> > genParticleTag_;
      edm::EDGetTokenT<vector<reco::GenJet> > genJetTag_;
      edm::EDGetTokenT<vector<reco::PFJet> > akJetTag_;
      edm::EDGetTokenT<vector<reco::Muon> > muonsTag_;
      edm::EDGetTokenT<edm::ValueMap<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > > muonMapTag_;
      edm::EDGetTokenT<edm::ValueMap<reco::PFJetRef> >  jetValMapTag_;
      edm::EDGetTokenT<vector<reco::PFTau> > tauCJTag_;
      edm::EDGetTokenT<vector<reco::RecoTauPiZero> > pizerosTag_;      
      edm::EDGetTokenT<reco::PFTauDiscriminator>  looseIsoTagCJ_;
      edm::EDGetTokenT<reco::PFTauDiscriminator>  medIsoTagCJ_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> tightIsoTagCJ_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> decayModeFindingTagCJ_;
      int genMatchPDGIDTag_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> isoRawTag_;
      edm::EDGetTokenT<reco::PFJetCollection>  oldJetTag_;
      edm::EDGetTokenT<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > mu3Tag_;
      edm::EDGetTokenT<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > mu12Tag_;
      edm::EDGetTokenT<reco::JetTagCollection> csvBTag_;
      //edm::EDGetTokenT<edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference> > csvBTag_;

      //Histograms
      TH1F* NEvents_;   
      TH1F* GenMatchType_; 
      TH1F* GenMatchPDGID_; 
      TH1F* GenMatchPDGIDwithLargeEta_;
      TH1F* NMuRemoved_;
      TH1F* TauMuTauHaddR_;
      TH1F* MassDiLepGEN_;
      TH1F* MassDiLepRECO_;
      TH1F* Mu3LargestPtMatch_;
      TH1F* Mu3PtMatch_;
      TH1F* MatchedTauMuPtLepton_;
      TH1F* MatchedTauHadPtLepton_;
      TH1F* MatchedTauHadEtaLepton_;
      TH1F* MatchedHTLepton_;
      TH1F* MatchedBDiscCSVLepton_;
      TH1F* MatchedTauMuPtJet_;
      TH1F* MatchedTauHadPtJet_;
      TH1F* MatchedTauHadEtaJet_;
      TH1F* MatchedHTJet_;
      TH1F* MatchedBDiscCSVJet_;
      TH1F* NEventsCuts_;
      TH1F* NConstituentsCJ_;
      TH1F* PtOverMassMu1Mu2_;
      TH1F* PtMu1Mu2_;

      TH1F* NTauDecayModeGEN_;
      TH1F* NTauDecayModeRECOCJ_;
      TH1F* NTauDecayModeRECONoCJ_;
      TH1F* NTauDecayModeRECOCJIfGen3_;
      TH1F* NTauDecayModeRECONoCJIfGen3_;
      TH2F* NTausRECOvsCLEANJETS_;
      TH2F* GenDiTaudRvsCJDiTaudR_;
      TH2F* TauHadnConstvsPt_;
      TH2F* RECOPtvsGENVisiblePtCJ_;
      TH2F* RECOPtvsGENVisiblePtRECO_;
      TH2F* RECOCJTauDMvsGENTauDM_;
      TH2F* RECONoCJTauDMvsGENTauDM_;
      TH2F* RECOCJTauDMvsGENTauDMMedIso_;
      TH2F* RECONoCJTauDMvsGENTauDMMedIso_;
      TH2F* GenVisiblePtvsRecoCJPtIfCorrDM_;
      TH2F* GenVisiblePtvsRecoNoCJPtIfCorrDM_;
      TH2F* RECOPtvsGENdRCJ_;
      TH2F* RECOPtvsGENdRRECO_;
      TH2F* RECOPtvsGENdRCJMedIso_;
      TH2F* RECOPtvsGENdRRECOMedIso_;
      TH2F* MatchedTauHadEtavsTauHadPt_;
      TH2F* BDiscCSVvsRawIso_;


      TH1F* MatchedLooseIsoRECOPt_;
      TH1F* MatchedMedIsoRECOPt_;
      TH1F* MatchedTightIsoRECOPt_;
      TH1F* MatchedDMFindRECOPt_;
      TH1F* MatchedRECOPt_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOPt_;
      TGraphAsymmErrors* FinalEffMedIsoRECOPt_;
      TGraphAsymmErrors* FinalEffTightIsoRECOPt_;
      TGraphAsymmErrors* FinalEffDMFindRECOPt_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOPtOverDM_;
      TGraphAsymmErrors* FinalEffMedIsoRECOPtOverDM_;
      TGraphAsymmErrors* FinalEffTightIsoRECOPtOverDM_;

      TH1F* MatchedLooseIsoCJPt_;
      TH1F* MatchedMedIsoCJPt_;
      TH1F* MatchedTightIsoCJPt_;
      TH1F* MatchedDMFindCJPt_;
      TH1F* MatchedCJPt_;
      TGraphAsymmErrors* FinalEffLooseIsoCJPt_;
      TGraphAsymmErrors* FinalEffMedIsoCJPt_;
      TGraphAsymmErrors* FinalEffTightIsoCJPt_;
      TGraphAsymmErrors* FinalEffDMFindCJPt_;
      TGraphAsymmErrors* FinalEffLooseIsoCJPtOverDM_;
      TGraphAsymmErrors* FinalEffMedIsoCJPtOverDM_;
      TGraphAsymmErrors* FinalEffTightIsoCJPtOverDM_;


      TH1F* MatchedLooseIsoRECOdR_;
      TH1F* MatchedMedIsoRECOdR_;
      TH1F* MatchedTightIsoRECOdR_;
      TH1F* MatchedDMFindRECOdR_;
      TH1F* MatchedRECOdR_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOdR_;
      TGraphAsymmErrors* FinalEffMedIsoRECOdR_;
      TGraphAsymmErrors* FinalEffTightIsoRECOdR_;
      TGraphAsymmErrors* FinalEffDMFindRECOdR_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOdROverDM_;
      TGraphAsymmErrors* FinalEffMedIsoRECOdROverDM_;
      TGraphAsymmErrors* FinalEffTightIsoRECOdROverDM_;

      TH1F* MatchedLooseIsoCJdR_;
      TH1F* MatchedMedIsoCJdR_;
      TH1F* MatchedTightIsoCJdR_;
      TH1F* MatchedDMFindCJdR_;
      TH1F* MatchedCJdR_;
      TGraphAsymmErrors* FinalEffLooseIsoCJdR_;
      TGraphAsymmErrors* FinalEffMedIsoCJdR_;
      TGraphAsymmErrors* FinalEffTightIsoCJdR_;
      TGraphAsymmErrors* FinalEffDMFindCJdR_;
      TGraphAsymmErrors* FinalEffLooseIsoCJdROverDM_;
      TGraphAsymmErrors* FinalEffMedIsoCJdROverDM_;
      TGraphAsymmErrors* FinalEffTightIsoCJdROverDM_;


      TH1F* MatchedLooseIsoCJdRPtCut_;
      TH1F* MatchedMedIsoCJdRPtCut_;
      TH1F* MatchedTightIsoCJdRPtCut_;
      TH1F* MatchedDMFindCJdRPtCut_;
      TH1F* MatchedCJdRPtCut_;
      TGraphAsymmErrors* FinalEffLooseIsoCJdRPtCut_;
      TGraphAsymmErrors* FinalEffMedIsoCJdRPtCut_;
      TGraphAsymmErrors* FinalEffTightIsoCJdRPtCut_;
      TGraphAsymmErrors* FinalEffDMFindCJdRPtCut_;

      TH1F* MatchedLooseIsoRECOdRPtCut_;
      TH1F* MatchedMedIsoRECOdRPtCut_;
      TH1F* MatchedTightIsoRECOdRPtCut_;
      TH1F* MatchedDMFindRECOdRPtCut_;
      TH1F* MatchedRECOdRPtCut_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOdRPtCut_;
      TGraphAsymmErrors* FinalEffMedIsoRECOdRPtCut_;
      TGraphAsymmErrors* FinalEffTightIsoRECOdRPtCut_;
      TGraphAsymmErrors* FinalEffDMFindRECOdRPtCut_;

      TH1F* MatchedLooseIsoRECOPtGen_;
      TH1F* MatchedMedIsoRECOPtGen_;
      TH1F* MatchedTightIsoRECOPtGen_;
      TH1F* MatchedDMFindRECOPtGen_;
      TH1F* MatchedRECOPtGen_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOPtGen_;
      TGraphAsymmErrors* FinalEffMedIsoRECOPtGen_;
      TGraphAsymmErrors* FinalEffTightIsoRECOPtGen_;
      TGraphAsymmErrors* FinalEffDMFindRECOPtGen_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOPtGenOverDM_;
      TGraphAsymmErrors* FinalEffMedIsoRECOPtGenOverDM_;
      TGraphAsymmErrors* FinalEffTightIsoRECOPtGenOverDM_;


      TH1F* MatchedLooseIsoCJPtGen_;
      TH1F* MatchedMedIsoCJPtGen_;
      TH1F* MatchedTightIsoCJPtGen_;
      TH1F* MatchedDMFindCJPtGen_;
      TH1F* MatchedCJPtGen_;
      TGraphAsymmErrors* FinalEffLooseIsoCJPtGen_;
      TGraphAsymmErrors* FinalEffMedIsoCJPtGen_;
      TGraphAsymmErrors* FinalEffTightIsoCJPtGen_;
      TGraphAsymmErrors* FinalEffDMFindCJPtGen_;
      TGraphAsymmErrors* FinalEffLooseIsoCJPtGenOverDM_;
      TGraphAsymmErrors* FinalEffMedIsoCJPtGenOverDM_;
      TGraphAsymmErrors* FinalEffTightIsoCJPtGenOverDM_;


      TH1F* OneProngDMCJPt_;
      TH1F* OneProngOnePizDMCJPt_;
      TH1F* OneProngTwoPizDMCJPt_;
      TH1F* ThreeProngDMCJPt_;
      TH1F* OneProngDMRECOPt_;
      TH1F* OneProngOnePizDMRECOPt_;
      TH1F* OneProngTwoPizDMRECOPt_;
      TH1F* ThreeProngDMRECOPt_;
      TH1F* MatchedOneProngCJPt_;
      TH1F* MatchedOneProngOnePizCJPt_;
      TH1F* MatchedOneProngTwoPizCJPt_;
      TH1F* MatchedThreeProngCJPt_;
      TH1F* MatchedOneProngRECOPt_;
      TH1F* MatchedOneProngOnePizRECOPt_;
      TH1F* MatchedOneProngTwoPizRECOPt_;
      TH1F* MatchedThreeProngRECOPt_;
      TGraphAsymmErrors* FinalOneProngDMCJPt_;
      TGraphAsymmErrors* FinalOneProngOnePizDMCJPt_;
      TGraphAsymmErrors* FinalOneProngTwoPizDMCJPt_;
      TGraphAsymmErrors* FinalThreeProngDMCJPt_;
      TGraphAsymmErrors* FinalOneProngDMRECOPt_;
      TGraphAsymmErrors* FinalOneProngOnePizDMRECOPt_;
      TGraphAsymmErrors* FinalOneProngTwoPizDMRECOPt_;
      TGraphAsymmErrors* FinalThreeProngDMRECOPt_;
      TMultiGraph* FinalEffOneProngDMSAMEPt_;
      TMultiGraph* FinalEffOneProngOnePizDMSAMEPt_;
      TMultiGraph* FinalEffOneProngTwoPizDMSAMEPt_;
      TMultiGraph* FinalEffThreeProngDMSAMEPt_;



      TMultiGraph* FinalEffLooseIsoSAMEPt_;
      TMultiGraph* FinalEffMedIsoSAMEPt_;
      TMultiGraph* FinalEffTightIsoSAMEPt_;
      TMultiGraph* FinalEffDMFindSAMEPt_;
      TMultiGraph* FinalEffLooseIsoSAMEdR_;
      TMultiGraph* FinalEffMedIsoSAMEdR_;
      TMultiGraph* FinalEffTightIsoSAMEdR_;
      TMultiGraph* FinalEffDMFindSAMEdR_;
      TMultiGraph* FinalEffLooseIsoSAMEdRPtCut_;
      TMultiGraph* FinalEffMedIsoSAMEdRPtCut_;
      TMultiGraph* FinalEffTightIsoSAMEdRPtCut_;
      TMultiGraph* FinalEffDMFindSAMEdRPtCut_;
      TMultiGraph* FinalEffLooseIsoSAMEPtGen_;
      TMultiGraph* FinalEffMedIsoSAMEPtGen_;
      TMultiGraph* FinalEffTightIsoSAMEPtGen_;
      TMultiGraph* FinalEffDMFindSAMEPtGen_;

      TMultiGraph* FinalEffLooseIsoSAMEPtOverDM_;
      TMultiGraph* FinalEffMedIsoSAMEPtOverDM_;
      TMultiGraph* FinalEffTightIsoSAMEPtOverDM_;
      TMultiGraph* FinalEffLooseIsoSAMEdROverDM_;
      TMultiGraph* FinalEffMedIsoSAMEdROverDM_;
      TMultiGraph* FinalEffTightIsoSAMEdROverDM_;
      TMultiGraph* FinalEffLooseIsoSAMEPtGenOverDM_;
      TMultiGraph* FinalEffMedIsoSAMEPtGenOverDM_;
      TMultiGraph* FinalEffTightIsoSAMEPtGenOverDM_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FakeRateAnalyzer::FakeRateAnalyzer(const edm::ParameterSet& iConfig):
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  genParticleTag_(consumes<vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleTag"))),
  genJetTag_(consumes<vector<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("genJetTag"))),
  akJetTag_(consumes<vector<reco::PFJet> >(iConfig.getParameter<edm::InputTag>("akJetTag"))),
  muonsTag_(consumes<vector<reco::Muon> >(iConfig.getParameter<edm::InputTag>("muonsTag"))),
  muonMapTag_(consumes<edm::ValueMap<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > >(iConfig.getParameter<edm::InputTag>("muonMapTag"))),
  jetValMapTag_(consumes<edm::ValueMap<reco::PFJetRef> >(iConfig.getParameter<edm::InputTag>("jetValMapTag"))),
  tauCJTag_(consumes<vector<reco::PFTau> >(iConfig.getParameter<edm::InputTag>("tauCJTag"))),
  pizerosTag_(consumes<vector<reco::RecoTauPiZero> >(iConfig.getParameter<edm::InputTag>("pizerosTag"))),
  looseIsoTagCJ_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("looseIsoTagCJ"))),
  medIsoTagCJ_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("medIsoTagCJ"))),
  tightIsoTagCJ_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("tightIsoTagCJ"))),
  decayModeFindingTagCJ_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("decayModeFindingTagCJ"))),
  genMatchPDGIDTag_(iConfig.getParameter<int>("genMatchPDGIDTag")),
  isoRawTag_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("isoRawTag"))),
  oldJetTag_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("oldJetTag"))),
  mu3Tag_(consumes<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > >(iConfig.getParameter<edm::InputTag>("mu3Tag"))),
  mu12Tag_(consumes<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > >(iConfig.getParameter<edm::InputTag>("mu12Tag"))),
  csvBTag_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("csvBTag")))
  //csvBTag_(consumes<edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference> >(iConfig.getParameter<edm::InputTag>("csvBTag")))
{
  std::cout << "genMatchPDGIDTag_= " << genMatchPDGIDTag_ << std::endl;
  reset(false);    
}//FakeRateAnalyzer



FakeRateAnalyzer::~FakeRateAnalyzer()
{
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void FakeRateAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "\n<------------THIS IS A NEW EVENT------------>" << std::endl;
  NEvents_->Fill(0);

  //Get gen particle collection
  edm::Handle<reco::GenParticleCollection> pGenParts;
  iEvent.getByToken(genParticleTag_, pGenParts);

  //Get ak4Jets particle collection
  edm::Handle<std::vector<reco::GenJet> > pGenJets;
  iEvent.getByToken(genJetTag_, pGenJets);

  //Get ak4Jets particle collection
  edm::Handle<std::vector<reco::PFJet> > pAkJets;
  iEvent.getByToken(akJetTag_, pAkJets);

  //Get RECO Muons particle collection
  edm::Handle<std::vector<reco::Muon> > pMuons;
  iEvent.getByToken(muonsTag_, pMuons);

  //get jet-muon map
  edm::Handle<edm::ValueMap<reco::MuonRefVector> > pMuonMap;
  iEvent.getByToken(muonMapTag_, pMuonMap);

  //get old-new jet map
  edm::Handle<edm::ValueMap<reco::PFJetRef> > pJetValMap;
  iEvent.getByToken(jetValMapTag_, pJetValMap);

  //Get CleanJets Tau particle collection
  edm::Handle<std::vector<reco::PFTau> > pTausCJ;
  iEvent.getByToken(tauCJTag_, pTausCJ);

  //Get pi-Zero particle collection
  edm::Handle<std::vector<reco::RecoTauPiZero> > pPiZero;
  iEvent.getByToken(pizerosTag_, pPiZero);

  //Get Loose Iso Collection
  Handle<PFTauDiscriminator> pLooseIsoDiscCJ; 
  iEvent.getByToken(looseIsoTagCJ_, pLooseIsoDiscCJ); 

  //Get Medium Iso Collection
  Handle<PFTauDiscriminator> pMedIsoDiscCJ; 
  iEvent.getByToken(medIsoTagCJ_, pMedIsoDiscCJ);

  //Get Tight Iso Collection
  Handle<PFTauDiscriminator> pTightIsoDiscCJ; 
  iEvent.getByToken(tightIsoTagCJ_, pTightIsoDiscCJ);

  //Get Decay Mode Finding Collection
  Handle<PFTauDiscriminator> pDMFindingCJ; 
  iEvent.getByToken(decayModeFindingTagCJ_, pDMFindingCJ);

  //Get IsoRaw  Collection
  Handle<PFTauDiscriminator> pIsoRaw;
  iEvent.getByToken(isoRawTag_, pIsoRaw);

  //Old Jet collection for bTagging
  edm::Handle<reco::PFJetCollection> pOldJets;
  iEvent.getByToken(oldJetTag_, pOldJets);

  //Old Jet collection for bTagging
  edm::Handle<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > pMu3;
  iEvent.getByToken(mu3Tag_, pMu3);

  //Old Jet collection for bTagging
  edm::Handle<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > pMu12;
  iEvent.getByToken(mu12Tag_, pMu12);

  //Get combVertMVA JetTagCollection
  edm::Handle<reco::JetTagCollection> pCSV;
  iEvent.getByToken(csvBTag_, pCSV);

//////////////////////////////
// Begin Analyzer
//////////////////////////////
  ///////////////////////////
  // iterating over the taus
  /////////////////////////// 
  unsigned int DMCJ = 0, MedIsoCJ = 0, LooseIsoCJ = 0, TightIsoCJ = 0;
  double nConstituents = 1000, rawIsoValue = -1;
  for (std::vector<reco::PFTau>::const_iterator iTauCJ = pTausCJ->begin(); iTauCJ != pTausCJ->end(); ++iTauCJ)
  {
    const reco::PFJetRef& tauJetRef = (*iTauCJ).jetRef();
    const reco::MuonRefVector& removedMuons = (*pMuonMap)[tauJetRef];
    std::vector<reco::PFCandidatePtr> JetPFCands = tauJetRef->getPFConstituents();
    nConstituents = JetPFCands.size();
    reco::PFTauRef PFTauRef(pTausCJ, iTauCJ - pTausCJ->begin());

    //find the highest pT associated muon
    std::vector<reco::MuonRef> removedMuonRefs;
    for (reco::MuonRefVector::const_iterator iMuon = removedMuons.begin(); iMuon != removedMuons.end(); ++iMuon) {removedMuonRefs.push_back(*iMuon);}//for iMuon
    std::cout << "\tremovedMuonRefs.size= " << removedMuonRefs.size() << std::endl;
    for (unsigned int iter = 0; iter < removedMuonRefs.size(); iter++)
    {
      for (unsigned int jter = iter + 1; jter < removedMuonRefs.size(); jter++)
      {
        if (removedMuonRefs[jter]->pt() > removedMuonRefs[iter]->pt())
	{
 	  reco::MuonRef TEMPRef = removedMuonRefs[iter];
	  removedMuonRefs[iter] = removedMuonRefs[jter];  
	  removedMuonRefs[jter] = TEMPRef;
        }//if jter > iter
      }//for jter
    }//for iter

    //Calculate Tau_mu tau_H dR
    double Pt_GenMatchCJ = -1, dR_tauRemovedMu = 10000;
    reco::MuonRef removedMuonRef;
    bool removedMu = false;
    int tauDecayModeRECO = -1;
    for (unsigned int iter = 0; iter < removedMuonRefs.size(); iter++)
    {
      double dPhi = reco::deltaPhi(removedMuonRefs[iter]->phi(), iTauCJ->phi() );
      double dR_tauMu = sqrt( (removedMuonRefs[iter]->eta() - iTauCJ->eta() ) * (removedMuonRefs[iter]->eta() - iTauCJ->eta() )  +  dPhi * dPhi );
      std::cout << "\t\t\tMuRef->pt()= " << removedMuonRefs[iter]->pt() << "  \tdR_tauMu= " << dR_tauMu << std::endl;
      if (dR_tauMu < .5)
      {
        NEvents_->Fill(1);
	NConstituentsCJ_->Fill(nConstituents );
        TauHadnConstvsPt_->Fill(nConstituents, iTauCJ->pt() );
        TauMuTauHaddR_->Fill(dR_tauMu );
 	GenDiTaudRvsCJDiTaudR_->Fill(dR_tauMu, dR_tauMu);       
        DMCJ = (*pDMFindingCJ)[PFTauRef];
        MedIsoCJ = (*pMedIsoDiscCJ)[PFTauRef];
        LooseIsoCJ = (*pLooseIsoDiscCJ)[PFTauRef];
        TightIsoCJ = (*pTightIsoDiscCJ)[PFTauRef];
        rawIsoValue = (*pIsoRaw)[PFTauRef];
        tauDecayModeRECO = iTauCJ->decayMode();
	removedMu = true;
	dR_tauRemovedMu = dR_tauMu;
	Pt_GenMatchCJ = iTauCJ->pt();
	removedMuonRef = removedMuonRefs[iter];
	break;
      }//if
    }//for itr

std::cout << "removedMu= " << removedMu << "  abs(iTauCJ->eta() )= " << abs(iTauCJ->eta() ) << " Pt_GenMatchCJ= " << Pt_GenMatchCJ << std::endl;


    if (removedMu && abs(iTauCJ->eta() ) < 2.3 && Pt_GenMatchCJ > 20)
    {
      //Getting Mu1 and Mu2
      reco::MuonRef mu1Ref = reco::MuonRef((*pMu12)[0] );
      reco::MuonRef mu2Ref = reco::MuonRef((*pMu12)[1] );
      reco::LeafCandidate::LorentzVector diMuP4;
      diMuP4 = mu1Ref->p4();
      diMuP4 += mu2Ref->p4();
      PtOverMassMu1Mu2_->Fill(diMuP4.Pt() / diMuP4.M() );
      PtMu1Mu2_->Fill(diMuP4.Pt() );  //mu1Ref->pt() + mu2Ref->pt() );

      ////////////////////////////////////////////////////////////////////////
      // Looping through the Signal Particles searching for the di-tau GM 
      ////////////////////////////////////////////////////////////////////////
      bool leptonGenMatched = false, jetGenMatched = false;
      double dRGMlep = 1000000000000000, dRGMjet = 1000000000000000;
      int leptonPdgIdGM = genMatchPDGIDTag_, jetPdgIdGM = genMatchPDGIDTag_;
      for (reco::GenParticleCollection::const_iterator iGenParticle = pGenParts->begin(); iGenParticle != pGenParts->end(); ++iGenParticle)
      {
        if (abs(iGenParticle->pdgId() ) == genMatchPDGIDTag_ && iGenParticle->status() == 62)
        {
	  // Loop through children
	  reco::LeafCandidate::LorentzVector daughterP4;
          for (unsigned int iChild = 0; iChild < iGenParticle->numberOfDaughters(); ++iChild)
	  {
	    reco::GenParticleRef childRef = iGenParticle->daughterRef(iChild);
            if (iChild == 0)
              daughterP4 = childRef->p4();
	    else
	      daughterP4 += childRef->p4();
            double dPhiGM = reco::deltaPhi(iTauCJ->phi(), childRef->phi() ), dEtaGM = iTauCJ->eta() - childRef->eta(), dRGMTemp = sqrt( dEtaGM*dEtaGM + dPhiGM*dPhiGM );
            std::cout << "\tMom= " << iGenParticle->pdgId() << "\tchild= " << childRef->pdgId() << "\tdRGMTemp= " << dRGMTemp << std::endl;
            if (dRGMTemp < .05 && abs(childRef->pdgId() ) <= 16 && abs(childRef->pdgId() ) >= 11 && dRGMlep > dRGMTemp) 
            {
              leptonPdgIdGM = abs(childRef->pdgId() );
              leptonGenMatched = true;
              dRGMlep = dRGMTemp;
            }//if
            else if (dRGMTemp < .05 && dRGMjet > dRGMTemp)
            {
              jetGenMatched = true;
              jetPdgIdGM = abs(childRef->pdgId() );
              dRGMjet = dRGMTemp;
            }//else

	    //Loop Through grandchildren
            for (unsigned int iChild2 = 0; iChild2 < childRef->numberOfDaughters(); ++iChild2)
            {
              reco::GenParticleRef grandChildRef = childRef->daughterRef(iChild2 );
	      double dPhiGM2 = reco::deltaPhi(iTauCJ->phi(), grandChildRef->phi() ), dEtaGM2 = iTauCJ->eta() - grandChildRef->eta(), dRGM2Temp = sqrt( dEtaGM2*dEtaGM2 + dPhiGM2*dPhiGM2 );
              std::cout << "\tMom= " << iGenParticle->pdgId() << "\tchild= " << childRef->pdgId() << "\tgrandChild= " << grandChildRef->pdgId() << "\tdRGM2Temp= " << dRGM2Temp << std::endl;
              if (dRGM2Temp < .05)
              {
	 	if (abs(grandChildRef->pdgId() ) <= 16 && abs(grandChildRef->pdgId() ) >= 11 && dRGMlep < dRGM2Temp)
		{
		  leptonPdgIdGM = abs(grandChildRef->pdgId() );
                  leptonGenMatched = true;
		  dRGMlep = dRGM2Temp;
		}//if
		else if (dRGMjet < dRGM2Temp)
		{
		  jetGenMatched = true;
		  jetPdgIdGM = abs(grandChildRef->pdgId() );
		  dRGMjet = dRGM2Temp;
		}//else
              }//if drGM < .05
	    }//for iChild2
	  }//for iChild
	  MassDiLepGEN_->Fill(daughterP4.M() );
        }// if leptonGenMatched
      }//for iGen Particle

      /////////////////////////////////////////////////////////////////
      // Looping through the ak4PFJets searching for di-tau GM
      /////////////////////////////////////////////////////////////////
      for (reco::GenJetCollection::const_iterator iGenJet = pGenJets->begin(); iGenJet != pGenJets->end(); ++iGenJet)
      {
        double dPhiGM = reco::deltaPhi(iTauCJ->phi(), iGenJet->phi() ), dEtaGM = iTauCJ->eta() - iGenJet->eta(), dRGMTemp = sqrt( dEtaGM*dEtaGM + dPhiGM*dPhiGM );
        if (dRGMTemp < .05 && dRGMjet > dRGMTemp)
        {
	  dRGMjet = dRGMTemp;
	  jetGenMatched = true;
	  jetPdgIdGM = abs(iTauCJ->pdgId() );
	}//if
      }//for iGenJet

      if (leptonGenMatched  || jetGenMatched)
      {
	std::cout << "\tdRGMlep= " << dRGMlep << " \tdRGMjet= " << dRGMjet << std::endl;
        ///////////////////////////////////////
        // Filling Cut Option Histograms
        ///////////////////////////////////////
        std::cout << "\t\tRemoved Mu" << std::endl;
        NEventsCuts_->Fill(0);
        double sumHT = 0;
        size_t numJets = pAkJets->size();
 	for ( size_t iJet = 0; iJet < numJets; ++iJet ) 
	{
	  reco::PFJetRef jetRef(pAkJets, iJet);
          if (jetRef->pt() > 20.0 && tauJetRef.key() != jetRef.key() )
            sumHT += jetRef->pt();
        }//for pAkJets      
 
        double BDisc = -1;
        for (unsigned int iBTagInfo = 0; iBTagInfo != pCSV->size(); ++iBTagInfo)
        {
          if (reco::PFJetRef(pOldJets, iBTagInfo).key() == ((*iTauCJ).jetRef()).key() )
          {
            const reco::JetTagCollection& bTags = *(pCSV.product());
            BDisc = bTags[iBTagInfo].second;
          }//if
        }//for iBTagInfo
        BDiscCSVvsRawIso_->Fill(BDisc, rawIsoValue );

        if (leptonGenMatched && (!jetGenMatched || dRGMlep < dRGMjet) )
        {
          GenMatchPDGID_->Fill(leptonPdgIdGM );
	  if (abs(iTauCJ->eta() ) > 1.5 )
            GenMatchPDGIDwithLargeEta_->Fill(leptonPdgIdGM );
          MatchedTauHadEtaLepton_->Fill(iTauCJ->eta() );
          MatchedTauHadPtLepton_->Fill(iTauCJ->pt() );
          MatchedTauMuPtLepton_->Fill(removedMuonRefs[0]->pt() );
 	  sumHT = sumHT - mu1Ref->pt() - mu2Ref->pt() - removedMuonRef->pt(); ;
          MatchedHTLepton_->Fill(sumHT);
          MatchedBDiscCSVLepton_->Fill(BDisc );
	  MatchedTauHadEtavsTauHadPt_->Fill(iTauCJ->eta(), iTauCJ->pt() );
          if (leptonPdgIdGM <= 6)
            GenMatchType_->Fill(3);
          if (leptonPdgIdGM == 15 || leptonPdgIdGM == 16)
            GenMatchType_->Fill(0);
          if (leptonPdgIdGM == 13 || leptonPdgIdGM == 14)
            GenMatchType_->Fill(1);
          if (leptonPdgIdGM == 11 || leptonPdgIdGM == 12)
            GenMatchType_->Fill(2);
  
	  std::cout << "Lepton particel is closer" << std::endl;
          NEventsCuts_->Fill(18);
          if (BDisc < .825) //BDisc < .9 && sumHT < 300)
            NEventsCuts_->Fill(1);
          if (BDisc < .5) //BDisc < .9 && sumHT < 300 && MedIsoCJ == 1 && DMCJ == 1)
            NEventsCuts_->Fill(2);
          if (BDisc < .9 )
            NEventsCuts_->Fill(3);
          if (BDisc < .9 && MedIsoCJ == 1 && DMCJ == 1)
            NEventsCuts_->Fill(4);
          if (sumHT < 300)
            NEventsCuts_->Fill(5);
          if (BDisc < .95) //sumHT < 300 && MedIsoCJ == 1 && DMCJ == 1)
            NEventsCuts_->Fill(6);
          if (sumHT < 400 && BDisc < .9)
           NEventsCuts_->Fill(7);
          if (sumHT < 400 && BDisc < .9 && MedIsoCJ == 1 && DMCJ == 1)
            NEventsCuts_->Fill(8);
          if (MedIsoCJ == 1 && DMCJ == 1)
            NEventsCuts_->Fill(20);
          if (LooseIsoCJ == 1 && DMCJ == 1)
            NEventsCuts_->Fill(22);
          if (TightIsoCJ == 1 && DMCJ == 1)
            NEventsCuts_->Fill(23);
          NEventsCuts_->Fill(21);
        }// Lepton Counting

        else if (jetGenMatched && (!leptonGenMatched || dRGMlep > dRGMjet) )    
        {
          GenMatchPDGID_->Fill(jetPdgIdGM );
          if (abs(iTauCJ->eta() ) > 1.5 )
            GenMatchPDGIDwithLargeEta_->Fill(jetPdgIdGM );
          MatchedTauHadEtaJet_->Fill(iTauCJ->eta() );
          MatchedTauHadPtJet_->Fill(iTauCJ->pt() );
          MatchedTauMuPtJet_->Fill(removedMuonRefs[0]->pt() );
          sumHT = sumHT - mu1Ref->pt() - mu2Ref->pt() - removedMuonRef->pt();
          MatchedHTJet_->Fill(sumHT);
          MatchedBDiscCSVJet_->Fill(BDisc );
          MatchedTauHadEtavsTauHadPt_->Fill(iTauCJ->eta(), iTauCJ->pt() );
          if (jetPdgIdGM <= 6)
            GenMatchType_->Fill(3);
          if (jetPdgIdGM == 15 || jetPdgIdGM == 16)
            GenMatchType_->Fill(0);
          if (jetPdgIdGM == 13 || jetPdgIdGM == 14)
            GenMatchType_->Fill(1);
          if (jetPdgIdGM == 11 || jetPdgIdGM == 12)
            GenMatchType_->Fill(2);

          std::cout << "Jet particel is closer" << std::endl;
          NEventsCuts_->Fill(17);
          if (BDisc < .825) //BDisc < .9 && sumHT < 300)
            NEventsCuts_->Fill(9);
          if (BDisc < .5 ) // BDisc < .9 && sumHT < 300 && MedIsoCJ == 1 && DMCJ == 1)
            NEventsCuts_->Fill(10);
          if (BDisc < .9 )
            NEventsCuts_->Fill(11);
          if (BDisc < .9 && MedIsoCJ == 1 && DMCJ == 1)
            NEventsCuts_->Fill(12);
          if (sumHT < 300)
            NEventsCuts_->Fill(13);
          if (BDisc < .95) //sumHT < 300 && MedIsoCJ == 1 && DMCJ == 1)
            NEventsCuts_->Fill(14);
          if (sumHT < 400 && BDisc < .9)
            NEventsCuts_->Fill(15);
          if (sumHT < 400 && BDisc < .9 && MedIsoCJ == 1 && DMCJ == 1)
            NEventsCuts_->Fill(16);
	  if (MedIsoCJ == 1 && DMCJ == 1)
	    NEventsCuts_->Fill(20);
          if (LooseIsoCJ == 1 && DMCJ == 1)
            NEventsCuts_->Fill(22);
          if (TightIsoCJ == 1 && DMCJ == 1)
            NEventsCuts_->Fill(23);
          NEventsCuts_->Fill(21);
        }//if leptonGenMatched  

        if (jetGenMatched && leptonGenMatched)
          NEventsCuts_->Fill(19);
        
	/////////////////////////////////
 	// Checking if removed mu in mu3
	/////////////////////////////////
	double largestPt = -1;
	for (reco::MuonRefVector::const_iterator iMuon = pMu3->begin(); iMuon!=pMu3->end(); ++iMuon)
	{
	  if ((**iMuon).pt() > largestPt)
	    largestPt = (**iMuon).pt();
	}//for iMuon

        for (reco::MuonRefVector::const_iterator iMuon = pMu3->begin(); iMuon!=pMu3->end(); ++iMuon)
        {
	  if (deltaR(**iMuon, *removedMuonRef) < .0001 && (fabs((**iMuon).pt() - removedMuonRef->pt() ) / (**iMuon).pt() ) < .0001 )
	  {
	    Mu3PtMatch_->Fill((**iMuon).pt() );
	    if ( (fabs((**iMuon).pt() - largestPt) / largestPt ) < .0001)
	      Mu3LargestPtMatch_->Fill((**iMuon).pt() );
	  }// if removed muon and mu3 match
	}// for Mu3

        ///////////////////////////
        //Check Di-Muon resonanice
        ///////////////////////////
        MassDiLepRECO_->Fill(diMuP4.M() );
 
	/////////////////////////////////////////////////////////////////
	//Getting the Event Counts for the Reconstruction Efficiencies
	/////////////////////////////////////////////////////////////////
        NTauDecayModeRECOCJ_->Fill(tauDecayModeRECO );
        NEvents_->Fill(3);
        MatchedCJPt_->Fill(Pt_GenMatchCJ );
 	std::cout << "DMCJ= " << DMCJ << "\tTightIsoCJ= " << TightIsoCJ << "\tMedIsoCJ= " << MedIsoCJ << "\tLooseIsoCJ= " << LooseIsoCJ << std::endl;
        if (DMCJ == 1)
          MatchedDMFindCJPt_->Fill(Pt_GenMatchCJ );
        if (TightIsoCJ == 1 && DMCJ == 1)
          MatchedTightIsoCJPt_->Fill(Pt_GenMatchCJ );
        if (MedIsoCJ == 1 && DMCJ == 1)
          MatchedMedIsoCJPt_->Fill(Pt_GenMatchCJ );
        if (LooseIsoCJ == 1 && DMCJ == 1)
          MatchedLooseIsoCJPt_->Fill(Pt_GenMatchCJ );
      
        MatchedCJdR_->Fill(dR_tauRemovedMu );
        if (DMCJ == 1)
          MatchedDMFindCJdR_->Fill(dR_tauRemovedMu );
        if (TightIsoCJ == 1 && DMCJ == 1)
          MatchedTightIsoCJdR_->Fill(dR_tauRemovedMu );
        if (MedIsoCJ == 1 && DMCJ == 1)
          MatchedMedIsoCJdR_->Fill(dR_tauRemovedMu );
        if (LooseIsoCJ == 1 && DMCJ == 1)
          MatchedLooseIsoCJdR_->Fill(dR_tauRemovedMu );
    
        if (Pt_GenMatchCJ > 100)
        {
          MatchedCJdRPtCut_->Fill(dR_tauRemovedMu );
          if (DMCJ == 1)
            MatchedDMFindCJdRPtCut_->Fill(dR_tauRemovedMu );
          if (TightIsoCJ == 1 && DMCJ == 1)
            MatchedTightIsoCJdRPtCut_->Fill(dR_tauRemovedMu );
          if (MedIsoCJ == 1 && DMCJ == 1)
            MatchedMedIsoCJdRPtCut_->Fill(dR_tauRemovedMu );
          if (LooseIsoCJ == 1 && DMCJ == 1)
            MatchedLooseIsoCJdRPtCut_->Fill(dR_tauRemovedMu );
        }//if Pt_GenMatchCJ > 100
      }//if jetGenMatched
    }//if GEN tau1Ref is the had in mu+had and it is matched to CleanJets Jet
  }//iTauCJ
}//End FakeRateAnalyzer::analyze


// ------------ method called once each job just before starting event loop  ------------
void FakeRateAnalyzer::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  //Open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //Book histograms
  NEvents_     = new TH1F("NEvents"    , "", 9, -.5, 8.5);
      NEvents_->GetXaxis()->SetBinLabel(1, "TotalEvents"); 
      NEvents_->GetXaxis()->SetBinLabel(2, "#tau_{#mu} + #tau_{had} Match");
      NEvents_->GetXaxis()->SetBinLabel(3, "Gen #tau_{#mu} + #tau_{had}");
      NEvents_->GetXaxis()->SetBinLabel(4, "CJ Gen Match #tau_{had}");
      NEvents_->GetXaxis()->SetBinLabel(5, "RECO Gen Match #tau_{had}");
      NEvents_->GetXaxis()->SetBinLabel(6, "Event with #tau_{#mu} Removed");
      NEvents_->GetXaxis()->SetBinLabel(7, "Event with no #tau_{#mu} Removed ");
  GenMatchType_     = new TH1F("GenMatchType"    , "", 4, -.5, 3.5);
      GenMatchType_->GetXaxis()->SetBinLabel(1, "#tau GenMatch");
      GenMatchType_->GetXaxis()->SetBinLabel(2, "#mu GenMatch");
      GenMatchType_->GetXaxis()->SetBinLabel(3, "elec GenMatch");
      GenMatchType_->GetXaxis()->SetBinLabel(4, "Had GenMatch");
  GenMatchPDGID_     = new TH1F("GenMatchPDGID"    , "", 1001, -.5, 1000.5);
  GenMatchPDGIDwithLargeEta_     = new TH1F("GenMatchPDGIDwithLargeEta"    , "", 1001, -.5, 1000.5);
  NMuRemoved_            = new TH1F("NMuRemoved"    , "", 2, -.5, 1.5);
      NMuRemoved_->GetXaxis()->SetBinLabel(1, "#tau_{#mu} + #tau_{Had} Event");
      NMuRemoved_->GetXaxis()->SetBinLabel(2, "and a Removed #mu");
  TauMuTauHaddR_       = new TH1F("TauMuTauHaddR"    , "", 25, 0, 8);
  MassDiLepGEN_       = new TH1F("MassDiLepGEN"    , "", 25, 0, 130);
  MassDiLepRECO_       = new TH1F("MassDiLepRECO"    , "", 25, 70, 120);
  Mu3LargestPtMatch_       = new TH1F("Mu3LargestPtMatch"    , "", 25, 0, 130);
  Mu3PtMatch_       = new TH1F("Mu3PtMatch"    , "", 25, 0, 130);
  MatchedTauMuPtLepton_       = new TH1F("MatchedTauMuPtLepton", "", 25, 0, 500);
  MatchedTauHadPtLepton_       = new TH1F("MatchedTauHadPtLepton" , "", 25, 0, 500);
  MatchedTauHadEtaLepton_       = new TH1F("MatchedTauHadEtaLepton", "", 25, -2.5, 2.5);
  MatchedHTLepton_       = new TH1F("MatchedHTLepton", "", 25, 0, 1000);
  MatchedBDiscCSVLepton_       = new TH1F("MatchedBDiscCSVLepton", "", 25, 0, 1);
  MatchedTauMuPtJet_       = new TH1F("MatchedTauMuPtJet", "", 25, 0, 500);
  MatchedTauHadPtJet_       = new TH1F("MatchedTauHadPtJet", "", 25, 0, 500);
  MatchedTauHadEtaJet_       = new TH1F("MatchedTauHadEtaJet", "", 25, -2.5, 2.5);
  MatchedHTJet_       = new TH1F("MatchedHTJet", "", 25, 0, 1000);
  MatchedBDiscCSVJet_       = new TH1F("MatchedBDiscCSVJet", "", 25, 0, 1);
  NEventsCuts_       = new TH1F("NEventsCuts", "", 24, -.5, 23.5);
      NEventsCuts_->GetXaxis()->SetBinLabel(1, "TotalEvents");
      NEventsCuts_->GetXaxis()->SetBinLabel(2, "L: BD < .825"); //BD < . 9 , HT < 300");
      NEventsCuts_->GetXaxis()->SetBinLabel(3, "L: BD < .5"); //BD < . 9 , HT < 300 + MedIso");
      NEventsCuts_->GetXaxis()->SetBinLabel(4, "L: BD < .9");
      NEventsCuts_->GetXaxis()->SetBinLabel(5, "L: BD < .9 + Med Iso");
      NEventsCuts_->GetXaxis()->SetBinLabel(6, "L: HT < 300");
      NEventsCuts_->GetXaxis()->SetBinLabel(7, "L: BD < .95");  //HT < 300 + MedIso");
      NEventsCuts_->GetXaxis()->SetBinLabel(8, "L: HT < 400 , BD < .9");
      NEventsCuts_->GetXaxis()->SetBinLabel(9, "L: HT < 400 , BD < .9 MedIso"); 
      NEventsCuts_->GetXaxis()->SetBinLabel(10, "J: BD < .825"); //BD < . 9 , HT < 300");
      NEventsCuts_->GetXaxis()->SetBinLabel(11, "J: BD < .5"); //BD < . 9 , HT < 300 + MedIso");
      NEventsCuts_->GetXaxis()->SetBinLabel(12, "J: BD < .9");
      NEventsCuts_->GetXaxis()->SetBinLabel(13, "J: BD < .9 + Med Iso");
      NEventsCuts_->GetXaxis()->SetBinLabel(14, "J: HT < 300");
      NEventsCuts_->GetXaxis()->SetBinLabel(15, "J: BD < .95");  // HT < 300 + MedIso");
      NEventsCuts_->GetXaxis()->SetBinLabel(16, "J: HT < 400 , BD < .9");
      NEventsCuts_->GetXaxis()->SetBinLabel(17, "J: HT < 400 , BD < .9 MedIso");
      NEventsCuts_->GetXaxis()->SetBinLabel(18, "Jet GM");
      NEventsCuts_->GetXaxis()->SetBinLabel(19, "Lep GM");
      NEventsCuts_->GetXaxis()->SetBinLabel(20, "Jet and Lep");
      NEventsCuts_->GetXaxis()->SetBinLabel(21, "Removed Mu + Med Iso");
      NEventsCuts_->GetXaxis()->SetBinLabel(22, "Removed Muon");
      NEventsCuts_->GetXaxis()->SetBinLabel(23, "Removed Mu + Loose Iso");
      NEventsCuts_->GetXaxis()->SetBinLabel(24, "Removed Mu + Tight Iso");
  NConstituentsCJ_        = new TH1F("NConstituentsCJ"    , "", 25, 0, 50);
  PtOverMassMu1Mu2_        = new TH1F("PtOverMassMu1Mu2"    , "", 25, 0, 50);
  PtMu1Mu2_        = new TH1F("PtMu1Mu2"    , "", 25, 0, 300);

  NTauDecayModeGEN_        = new TH1F("NTauDecayModeGEN"    , "", 12, -.5, 11.5);
  NTauDecayModeRECOCJ_        = new TH1F("NTauDecayModeRECOCJ"    , "", 12, -.5, 11.5);
  NTauDecayModeRECONoCJ_        = new TH1F("NTauDecayModeRECONoCJ"    , "", 12, -.5, 11.5);
  NTauDecayModeRECOCJIfGen3_        = new TH1F("NTauDecayModeRECOCJIfGen3"    , "", 12, -.5, 11.5);
  NTauDecayModeRECONoCJIfGen3_        = new TH1F("NTauDecayModeRECONoCJIfGen3"    , "", 12, -.5, 11.5);
  NTausRECOvsCLEANJETS_  = new TH2F("NTausRECOvsCLEANJETS" , "", 11, -.5, 10.5, 11, -.5, 10.5);
  GenDiTaudRvsCJDiTaudR_  = new TH2F("GenDiTaudRvsCJDiTaudR" , "", 50, 0, 8, 50, 0, 10);
  TauHadnConstvsPt_  = new TH2F("TauHadnConstvsPt" , "", 50, 0, 50, 30, 0, 160);
  RECOPtvsGENVisiblePtCJ_  = new TH2F("RECOPtvsGENVisiblePtCJ" , "", 50, 0, 160, 50, 0, 160);
  RECOPtvsGENVisiblePtRECO_  = new TH2F("RECOPtvsGENVisiblePtRECO" , "", 50, 0, 160, 50, 0, 160);
  RECOCJTauDMvsGENTauDM_  = new TH2F("RECOCJTauDMvsGENTauDM" , "", 4, .5, 4.5, 4, 0.5, 4.5);
  RECONoCJTauDMvsGENTauDM_  = new TH2F("RECONoCJTauDMvsGENTauDM" , "", 4, .5, 4.5, 4, 0.5, 4.5);
  RECOCJTauDMvsGENTauDMMedIso_  = new TH2F("RECOCJTauDMvsGENTauDMMedIso" , "", 4, .5, 4.5, 4, 0.5, 4.5);
  RECONoCJTauDMvsGENTauDMMedIso_  = new TH2F("RECONoCJTauDMvsGENTauDMMedIso" , "", 4, .5, 4.5, 4, 0.5, 4.5);
  GenVisiblePtvsRecoCJPtIfCorrDM_   = new TH2F("GenVisiblePtvsRecoCJPtIfCorrDM"  , "", 40, 0, 160, 40, 0, 160);
  GenVisiblePtvsRecoNoCJPtIfCorrDM_ = new TH2F("GenVisiblePtvsRecoNoCJPtIfCorrDM", "", 40, 0, 160, 40, 0, 160);
  RECOPtvsGENdRCJ_ = new TH2F("RECOPtvsGENdRCJ", "", 40, 0, 160, 40, 0, 4);
  RECOPtvsGENdRRECO_ = new TH2F("RECOPtvsGENdRRECO", "", 40, 0, 160, 40, 0, 4);
  RECOPtvsGENdRCJMedIso_ = new TH2F("RECOPtvsGENdRCJMedIso", "", 40, 0, 160, 40, 0, 4);
  RECOPtvsGENdRRECOMedIso_ = new TH2F("RECOPtvsGENdRRECOMedIso", "", 40, 0, 160, 40, 0, 4);
  MatchedTauHadEtavsTauHadPt_ = new TH2F("MatchedTauHadEtavsTauHadPt", "", 100, -2.5, 2.5, 100, 0, 300);
  BDiscCSVvsRawIso_ = new TH2F("BDiscCSVvsRawIso", "", 10, 0, 1, 20, 0, 20);

  RECOCJTauDMvsGENTauDM_->GetXaxis()->SetBinLabel(1, "1prong");
  RECOCJTauDMvsGENTauDM_->GetXaxis()->SetBinLabel(2, "1prong+1#pi^{0}");
  RECOCJTauDMvsGENTauDM_->GetXaxis()->SetBinLabel(3, "1prong+2#pi^{0}");
  RECOCJTauDMvsGENTauDM_->GetXaxis()->SetBinLabel(4, "3prong");
  RECOCJTauDMvsGENTauDM_->GetYaxis()->SetBinLabel(1, "1prong");
  RECOCJTauDMvsGENTauDM_->GetYaxis()->SetBinLabel(2, "1prong+1#pi^{0}");
  RECOCJTauDMvsGENTauDM_->GetYaxis()->SetBinLabel(3, "1prong+2#pi^{0}");
  RECOCJTauDMvsGENTauDM_->GetYaxis()->SetBinLabel(4, "3prong");
  
  RECONoCJTauDMvsGENTauDM_->GetXaxis()->SetBinLabel(1, "1prong");
  RECONoCJTauDMvsGENTauDM_->GetXaxis()->SetBinLabel(2, "1prong+1#pi^{0}");
  RECONoCJTauDMvsGENTauDM_->GetXaxis()->SetBinLabel(3, "1prong+2#pi^{0}");
  RECONoCJTauDMvsGENTauDM_->GetXaxis()->SetBinLabel(4, "3prong");
  RECONoCJTauDMvsGENTauDM_->GetYaxis()->SetBinLabel(1, "1prong");
  RECONoCJTauDMvsGENTauDM_->GetYaxis()->SetBinLabel(2, "1prong+1#pi^{0}");
  RECONoCJTauDMvsGENTauDM_->GetYaxis()->SetBinLabel(3, "1prong+2#pi^{0}");
  RECONoCJTauDMvsGENTauDM_->GetYaxis()->SetBinLabel(4, "3prong");

  RECOCJTauDMvsGENTauDMMedIso_->GetXaxis()->SetBinLabel(1, "1prong");
  RECOCJTauDMvsGENTauDMMedIso_->GetXaxis()->SetBinLabel(2, "1prong+1#pi^{0}");
  RECOCJTauDMvsGENTauDMMedIso_->GetXaxis()->SetBinLabel(3, "1prong+2#pi^{0}");
  RECOCJTauDMvsGENTauDMMedIso_->GetXaxis()->SetBinLabel(4, "3prong");
  RECOCJTauDMvsGENTauDMMedIso_->GetYaxis()->SetBinLabel(1, "1prong");
  RECOCJTauDMvsGENTauDMMedIso_->GetYaxis()->SetBinLabel(2, "1prong+1#pi^{0}");
  RECOCJTauDMvsGENTauDMMedIso_->GetYaxis()->SetBinLabel(3, "1prong+2#pi^{0}");
  RECOCJTauDMvsGENTauDMMedIso_->GetYaxis()->SetBinLabel(4, "3prong");

  RECONoCJTauDMvsGENTauDMMedIso_->GetXaxis()->SetBinLabel(1, "1prong");
  RECONoCJTauDMvsGENTauDMMedIso_->GetXaxis()->SetBinLabel(2, "1prong+1#pi^{0}");
  RECONoCJTauDMvsGENTauDMMedIso_->GetXaxis()->SetBinLabel(3, "1prong+2#pi^{0}");
  RECONoCJTauDMvsGENTauDMMedIso_->GetXaxis()->SetBinLabel(4, "3prong");
  RECONoCJTauDMvsGENTauDMMedIso_->GetYaxis()->SetBinLabel(1, "1prong");
  RECONoCJTauDMvsGENTauDMMedIso_->GetYaxis()->SetBinLabel(2, "1prong+1#pi^{0}");
  RECONoCJTauDMvsGENTauDMMedIso_->GetYaxis()->SetBinLabel(3, "1prong+2#pi^{0}");
  RECONoCJTauDMvsGENTauDMMedIso_->GetYaxis()->SetBinLabel(4, "3prong");


  MatchedLooseIsoRECOPt_  = new TH1F("MatchedLooseIsoRECOPt"    , "", 11, 0, 220.0);
  MatchedMedIsoRECOPt_    = new TH1F("MatchedMedIsoRECOPt", "", 11, 0, 220);
  MatchedTightIsoRECOPt_    = new TH1F("MatchedTightIsoRECOPt", "", 11, 0, 220);
  MatchedDMFindRECOPt_    = new TH1F("MatchedDMFindRECOPt"    , "", 11, 0, 220);
  MatchedRECOPt_          = new TH1F("MatchedRECOPt"    , "", 11, 0, 220);

  MatchedLooseIsoCJPt_  = new TH1F("MatchedLooseIsoCJPt"    , "", 11, 0, 220.0);
  MatchedMedIsoCJPt_    = new TH1F("MatchedMedIsoCJPt", "", 11, 0, 220);
  MatchedTightIsoCJPt_    = new TH1F("MatchedTightIsoCJPt", "", 11, 0, 220);
  MatchedDMFindCJPt_    = new TH1F("MatchedDMFindCJPt"    , "", 11, 0, 220);
  MatchedCJPt_          = new TH1F("MatchedCJPt"    , "", 11, 0, 220);

  FinalEffLooseIsoRECOPt_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOPt_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOPt_ = new TGraphAsymmErrors(30);
  FinalEffDMFindRECOPt_ = new TGraphAsymmErrors(30);
  FinalEffLooseIsoRECOPtOverDM_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOPtOverDM_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOPtOverDM_ = new TGraphAsymmErrors(30);

  FinalEffLooseIsoCJPt_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoCJPt_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoCJPt_ = new TGraphAsymmErrors(30);
  FinalEffDMFindCJPt_ = new TGraphAsymmErrors(30);
  FinalEffLooseIsoCJPtOverDM_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoCJPtOverDM_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoCJPtOverDM_ = new TGraphAsymmErrors(30);

  MatchedLooseIsoRECOdR_  = new TH1F("MatchedLooseIsoRECOdR"    , "", 8, 0, 1.5);
  MatchedMedIsoRECOdR_    = new TH1F("MatchedMedIsoRECOdR", "", 8, 0, 1.5);
  MatchedTightIsoRECOdR_    = new TH1F("MatchedTightIsoRECOdR", "", 8, 0, 1.5);
  MatchedDMFindRECOdR_    = new TH1F("MatchedDMFindRECOdR"    , "", 8, 0, 1.5);
  MatchedRECOdR_          = new TH1F("MatchedRECOdR"    , "", 8, 0, 1.5);

  MatchedLooseIsoCJdR_  = new TH1F("MatchedLooseIsoCJdR"    , "", 8, 0, 1.5);
  MatchedMedIsoCJdR_    = new TH1F("MatchedMedIsoCJdR", "", 8, 0, 1.5);
  MatchedTightIsoCJdR_    = new TH1F("MatchedTightIsoCJdR", "", 8, 0, 1.5);
  MatchedDMFindCJdR_    = new TH1F("MatchedDMFindCJdR"    , "", 8, 0, 1.5);
  MatchedCJdR_          = new TH1F("MatchedCJdR"    , "", 8, 0, 1.5);

  FinalEffLooseIsoRECOdR_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOdR_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOdR_ = new TGraphAsymmErrors(30);
  FinalEffDMFindRECOdR_ = new TGraphAsymmErrors(30);
  FinalEffLooseIsoRECOdROverDM_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOdROverDM_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOdROverDM_ = new TGraphAsymmErrors(30);

  FinalEffLooseIsoCJdR_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoCJdR_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoCJdR_ = new TGraphAsymmErrors(30);
  FinalEffDMFindCJdR_ = new TGraphAsymmErrors(30);
  FinalEffLooseIsoCJdROverDM_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoCJdROverDM_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoCJdROverDM_ = new TGraphAsymmErrors(30);

  MatchedLooseIsoRECOdRPtCut_  = new TH1F("MatchedLooseIsoRECOdRPtCut"    , "", 8, 0, 1.5);
  MatchedMedIsoRECOdRPtCut_    = new TH1F("MatchedMedIsoRECOdRPtCut", "", 8, 0, 1.5);
  MatchedTightIsoRECOdRPtCut_    = new TH1F("MatchedTightIsoRECOdRPtCut", "", 8, 0, 1.5);
  MatchedDMFindRECOdRPtCut_    = new TH1F("MatchedDMFindRECOdRPtCut"    , "", 8, 0, 1.5);
  MatchedRECOdRPtCut_          = new TH1F("MatchedRECOdRPtCut"    , "", 8, 0, 1.5);

  MatchedLooseIsoCJdRPtCut_  = new TH1F("MatchedLooseIsoCJdRPtCut"    , "", 8, 0, 1.5);
  MatchedMedIsoCJdRPtCut_    = new TH1F("MatchedMedIsoCJdRPtCut", "", 8, 0, 1.5);
  MatchedTightIsoCJdRPtCut_    = new TH1F("MatchedTightIsoCJdRPtCut", "", 8, 0, 1.5);
  MatchedDMFindCJdRPtCut_    = new TH1F("MatchedDMFindCJdRPtCut"    , "", 8, 0, 1.5);
  MatchedCJdRPtCut_          = new TH1F("MatchedCJdRPtCut"    , "", 8, 0, 1.5);

  FinalEffLooseIsoRECOdRPtCut_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOdRPtCut_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOdRPtCut_ = new TGraphAsymmErrors(30);
  FinalEffDMFindRECOdRPtCut_ = new TGraphAsymmErrors(30);

  FinalEffLooseIsoCJdRPtCut_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoCJdRPtCut_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoCJdRPtCut_ = new TGraphAsymmErrors(30);
  FinalEffDMFindCJdRPtCut_ = new TGraphAsymmErrors(30);


  MatchedLooseIsoRECOPtGen_  = new TH1F("MatchedLooseIsoRECOPtGen"    , "", 11, 0, 220.0);
  MatchedMedIsoRECOPtGen_    = new TH1F("MatchedMedIsoRECOPtGen", "", 11, 0, 220);
  MatchedTightIsoRECOPtGen_    = new TH1F("MatchedTightIsoRECOPtGen", "", 11, 0, 220);
  MatchedDMFindRECOPtGen_    = new TH1F("MatchedDMFindRECOPtGen"    , "", 11, 0, 220);
  MatchedRECOPtGen_          = new TH1F("MatchedRECOPtGen"    , "", 11, 0, 220);

  MatchedLooseIsoCJPtGen_  = new TH1F("MatchedLooseIsoCJPtGen"    , "", 11, 0, 220.0);
  MatchedMedIsoCJPtGen_    = new TH1F("MatchedMedIsoCJPtGen", "", 11, 0, 220);
  MatchedTightIsoCJPtGen_    = new TH1F("MatchedTightIsoCJPtGen", "", 11, 0, 220);
  MatchedDMFindCJPtGen_    = new TH1F("MatchedDMFindCJPtGen"    , "", 11, 0, 220);
  MatchedCJPtGen_          = new TH1F("MatchedCJPtGen"    , "", 11, 0, 220);

  FinalEffLooseIsoRECOPtGen_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOPtGen_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOPtGen_ = new TGraphAsymmErrors(30);
  FinalEffDMFindRECOPtGen_ = new TGraphAsymmErrors(30);
  FinalEffLooseIsoRECOPtGenOverDM_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOPtGenOverDM_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOPtGenOverDM_ = new TGraphAsymmErrors(30);

  FinalEffLooseIsoCJPtGen_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoCJPtGen_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoCJPtGen_ = new TGraphAsymmErrors(30);
  FinalEffDMFindCJPtGen_ = new TGraphAsymmErrors(30);
  FinalEffLooseIsoCJPtGenOverDM_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoCJPtGenOverDM_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoCJPtGenOverDM_ = new TGraphAsymmErrors(30);

  OneProngDMCJPt_ = new TH1F("OneProngDMCJPt"    , "", 11, 0, 220.0);
  OneProngOnePizDMCJPt_ = new TH1F("OneProngOnePizDMCJPt"    , "", 11, 0, 220.0);
  OneProngTwoPizDMCJPt_ = new TH1F("OneProngTwoPizDMCJPt"    , "", 11, 0, 220.0);
  ThreeProngDMCJPt_ = new TH1F("ThreeProngDMCJPt"    , "", 11, 0, 220.0);
  OneProngDMRECOPt_ = new TH1F("OneProngDMRECOPt"    , "", 11, 0, 220.0);
  OneProngOnePizDMRECOPt_ = new TH1F("OneProngOnePizDMRECOPt"    , "", 11, 0, 220.0);
  OneProngTwoPizDMRECOPt_ = new TH1F("OneProngTwoPizDMRECOPt"    , "", 11, 0, 220.0);
  ThreeProngDMRECOPt_ = new TH1F("ThreeProngDMRECOPt"    , "", 11, 0, 220.0);
  MatchedOneProngCJPt_ = new TH1F("MatchedOneProngCJPt", "", 11, 0, 220.0);
  MatchedOneProngOnePizCJPt_ = new TH1F("MatchedOneProngOnePizCJPt", "", 11, 0, 220.0);
  MatchedOneProngTwoPizCJPt_ = new TH1F("MatchedOneProngTwoPizCJPt", "", 11, 0, 220.0);
  MatchedThreeProngCJPt_ = new TH1F("MatchedThreeProngCJPt", "", 11, 0, 220.0);
  MatchedOneProngRECOPt_ = new TH1F("MatchedOneProngRECOPt", "", 11, 0, 220.0);
  MatchedOneProngOnePizRECOPt_ = new TH1F("MatchedOneProngOnePizRECOPt", "", 11, 0, 220.0);
  MatchedOneProngTwoPizRECOPt_ = new TH1F("MatchedOneProngTwoPizRECOPt", "", 11, 0, 220.0);
  MatchedThreeProngRECOPt_ = new TH1F("MatchedThreeProngRECOPt", "", 11, 0, 220.0);
  FinalOneProngDMCJPt_ = new TGraphAsymmErrors(30);
  FinalOneProngOnePizDMCJPt_ = new TGraphAsymmErrors(30);
  FinalOneProngTwoPizDMCJPt_ = new TGraphAsymmErrors(30);
  FinalThreeProngDMCJPt_ = new TGraphAsymmErrors(30);
  FinalOneProngDMRECOPt_ = new TGraphAsymmErrors(30);
  FinalOneProngOnePizDMRECOPt_ = new TGraphAsymmErrors(30);
  FinalOneProngTwoPizDMRECOPt_ = new TGraphAsymmErrors(30);
  FinalThreeProngDMRECOPt_ = new TGraphAsymmErrors(30);
  FinalEffOneProngDMSAMEPt_ = new TMultiGraph();
  FinalEffOneProngOnePizDMSAMEPt_ = new TMultiGraph();
  FinalEffOneProngTwoPizDMSAMEPt_ = new TMultiGraph();
  FinalEffThreeProngDMSAMEPt_ = new TMultiGraph();



  FinalEffLooseIsoSAMEPt_ = new TMultiGraph();
  FinalEffMedIsoSAMEPt_ = new TMultiGraph();
  FinalEffTightIsoSAMEPt_ = new TMultiGraph();
  FinalEffDMFindSAMEPt_ = new TMultiGraph();
  FinalEffLooseIsoSAMEdR_ = new TMultiGraph();
  FinalEffMedIsoSAMEdR_ = new TMultiGraph();
  FinalEffTightIsoSAMEdR_ = new TMultiGraph();
  FinalEffDMFindSAMEdR_ = new TMultiGraph();
  FinalEffLooseIsoSAMEdRPtCut_ = new TMultiGraph();
  FinalEffMedIsoSAMEdRPtCut_ = new TMultiGraph();
  FinalEffTightIsoSAMEdRPtCut_ = new TMultiGraph();
  FinalEffDMFindSAMEdRPtCut_ = new TMultiGraph();
  FinalEffLooseIsoSAMEPtGen_ = new TMultiGraph();
  FinalEffMedIsoSAMEPtGen_ = new TMultiGraph();
  FinalEffTightIsoSAMEPtGen_ = new TMultiGraph();
  FinalEffDMFindSAMEPtGen_ = new TMultiGraph();

  FinalEffLooseIsoSAMEPtOverDM_ = new TMultiGraph();
  FinalEffMedIsoSAMEPtOverDM_ = new TMultiGraph();
  FinalEffTightIsoSAMEPtOverDM_ = new TMultiGraph();
  FinalEffLooseIsoSAMEdROverDM_ = new TMultiGraph();
  FinalEffMedIsoSAMEdROverDM_ = new TMultiGraph();
  FinalEffTightIsoSAMEdROverDM_ = new TMultiGraph();
  FinalEffLooseIsoSAMEPtGenOverDM_ = new TMultiGraph();
  FinalEffMedIsoSAMEPtGenOverDM_ = new TMultiGraph();
  FinalEffTightIsoSAMEPtGenOverDM_ = new TMultiGraph();

}

// ------------ method called once each job just after ending the event loop  ------------
void FakeRateAnalyzer::endJob()
{
  //Make the Canvases
  TCanvas NEventsCanvas("NEvents","",600,600);
  TCanvas GenMatchTypeCanvas("GenMatchType","",600,600);
  TCanvas GenMatchPDGIDCanvas("GenMatchPDGID","",600,600);
  TCanvas GenMatchPDGIDwithLargeEtaCanvas("GenMatchPDGIDwithLargeEta","",600,600);
  TCanvas NMuRemovedCanvas("NMuRemoved","",600,600);
  TCanvas TauMuTauHaddRCanvas("TauMuTauHaddR","",600,600);
  TCanvas MassDiLepGENCanvas("MassDiLepGEN","",600,600);
  TCanvas MassDiLepRECOCanvas("MassDiLepRECO","",600,600);
  TCanvas Mu3LargestPtMatchCanvas("Mu3LargestPtMatch","",600,600);
  TCanvas Mu3PtMatchCanvas("Mu3PtMatch","",600,600);
  TCanvas MatchedTauMuPtLeptonCanvas("MatchedTauMuPtLepton","",600,600);
  TCanvas MatchedTauHadPtLeptonCanvas("MatchedTauHadPtLepton","",600,600);
  TCanvas MatchedTauHadEtaLeptonCanvas("MatchedTauHadEtaLepton","",600,600);
  TCanvas MatchedHTLeptonCanvas("MatchedHTLepton","",600,600); 
  TCanvas MatchedBDiscCSVLeptonCanvas("MatchedBDiscCSVLepton","",600,600);
  TCanvas MatchedTauMuPtJetCanvas("MatchedTauMuPtJet","",600,600);
  TCanvas MatchedTauHadPtJetCanvas("MatchedTauHadPtJet","",600,600);
  TCanvas MatchedTauHadEtaJetCanvas("MatchedTauHadEtaJet","",600,600);
  TCanvas MatchedHTJetCanvas("MatchedHTJet","",600,600);
  TCanvas MatchedBDiscCSVJetCanvas("MatchedBDiscCSVJet","",600,600);
  TCanvas NEventsCutsCanvas("NEventsCuts","",600,600);
  TCanvas NConstituentsCJCanvas("NConstituentsCJ","",600,600);
  TCanvas PtOverMassMu1Mu2Canvas("PtOverMassMu1Mu2","",600,600);
  TCanvas PtMu1Mu2Canvas("PtMu1Mu2","",600,600);

  TCanvas NTauDecayModeGENCanvas("NTauDecayModeGEN","",600,600);
  TCanvas NTauDecayModeRECOCJCanvas("NTauDecayModeRECOCJ","",600,600);
  TCanvas NTauDecayModeRECONoCJCanvas("NTauDecayModeRECONoCJ","",600,600);
  TCanvas NTauDecayModeRECOCJIfGen3Canvas("NTauDecayModeRECOCJIfGen3","",600,600);
  TCanvas NTauDecayModeRECONoCJIfGen3Canvas("NTauDecayModeRECONoCJIfGen3","",600,600);
  TCanvas NTausRECOvsCLEANJETSCanvas("NTausRECOvsCLEANJETS","",600,600);
  TCanvas GenDiTaudRvsCJDiTaudRCanvas("GenDiTaudRvsCJDiTaudR","",600,600);
  TCanvas TauHadnConstvsPtCanvas("TauHadnConstvsPt","",600,600);
  TCanvas RECOPtvsGENVisiblePtCJCanvas("RECOPtvsGENVisiblePtCJ","",600,600);
  TCanvas RECOPtvsGENVisiblePtRECOCanvas("RECOPtvsGENVisiblePtRECO","",600,600);
  TCanvas RECOCJTauDMvsGENTauDMCanvas("RECOCJTauDMvsGENTauDM","",600,600);
  TCanvas RECONoCJTauDMvsGENTauDMCanvas("RECONoCJTauDMvsGENTauDM","",600,600);
  TCanvas RECOCJTauDMvsGENTauDMMedIsoCanvas("RECOCJTauDMvsGENTauDMMedIso","",600,600);
  TCanvas RECONoCJTauDMvsGENTauDMMedIsoCanvas("RECONoCJTauDMvsGENTauDMMedIso","",600,600);
  TCanvas GenVisiblePtvsRecoCJPtIfCorrDMCanvas("GenVisiblePtvsRecoCJPtIfCorrDM","",600,600);
  TCanvas GenVisiblePtvsRecoNoCJPtIfCorrDMCanvas("GenVisiblePtvsRecoNoCJPtIfCorrDM","",600,600);
  TCanvas RECOPtvsGENdRCJCanvas("RECOPtvsGENdRCJ","",600,600);
  TCanvas RECOPtvsGENdRRECOCanvas("RECOPtvsGENdRRECO","",600,600);
  TCanvas RECOPtvsGENdRCJMedIsoCanvas("RECOPtvsGENdRCJMedIso","",600,600);
  TCanvas RECOPtvsGENdRRECOMedIsoCanvas("RECOPtvsGENdRRECOMedIso","",600,600);
  TCanvas MatchedTauHadEtavsTauHadPtCanvas("MatchedTauHadEtavsTauHadPt","",600,600);
  TCanvas BDiscCSVvsRawIsoCanvas("BDiscCSVvsRawIso","",600,600);

  std::cout << "\nRECOCJTauDMvsGENTauDM by column" << std::endl;
  for (int x = 1; x < 5; x++)
  {
    int sum = 0;
    for (int y = 1; y < 5; y++)
      sum += RECOCJTauDMvsGENTauDM_->GetBinContent(x,y);
    for (int y = 1; y < 5; y++)
      std::cout << "(" << x << "," << y << ")= " << RECOCJTauDMvsGENTauDM_->GetBinContent(x,y) << " Column average= " << RECOCJTauDMvsGENTauDM_->GetBinContent(x,y) / sum << std::endl;
  }//for x
  std::cout << "\n\nRECOCJTauDMvsGENTauDM by row" << std::endl;
  for (int y = 1; y < 5; y++)
  {
    int sum = 0;
    for (int x = 1; x < 5; x++)
      sum += RECOCJTauDMvsGENTauDM_->GetBinContent(x,y);
    for (int x = 1; x < 5; x++)
      std::cout << "(" << x << "," << y << ")= " << RECOCJTauDMvsGENTauDM_->GetBinContent(x,y) << "  Row Average= " << RECOCJTauDMvsGENTauDM_->GetBinContent(x,y) / sum << std::endl;
  }//for y
  std::cout << "\n\nRECOCJTauDMvsGENTauDM by total" << std::endl;
  int sumTotal = 0;
  for (int y = 1; y < 5; y++)
  {
    for (int x = 1; x < 5; x++)
      sumTotal += RECOCJTauDMvsGENTauDM_->GetBinContent(x,y);
  }//for y
  for (int y = 1; y < 5; y++)
  { 
    for (int x = 1; x < 5; x++)
      std::cout << "(" << x << "," << y << ")= " << RECOCJTauDMvsGENTauDM_->GetBinContent(x,y) << "  Total Average= " << RECOCJTauDMvsGENTauDM_->GetBinContent(x,y) / sumTotal << std::endl;
  }//for y


  std::cout << "\n\n\n\nRECONoCJTauDMvsGENTauDM by column" << std::endl;
  for (int x = 1; x < 5; x++)
  {
    int sum = 0;
    for (int y = 1; y < 5; y++)
      sum += RECONoCJTauDMvsGENTauDM_->GetBinContent(x,y);
    for (int y = 1; y < 5; y++)
      std::cout << "(" << x << "," << y << ")= " << RECONoCJTauDMvsGENTauDM_->GetBinContent(x,y) << "  Column Average= " << RECONoCJTauDMvsGENTauDM_->GetBinContent(x,y) / sum << std::endl;
  }//for x
  std::cout << "\n\nRECONoCJTauDMvsGENTauDM by row" << std::endl;
  for (int y = 1; y < 5; y++)
  {
    int sum = 0;
    for (int x = 1; x < 5; x++)
      sum += RECONoCJTauDMvsGENTauDM_->GetBinContent(x,y);
    for (int x = 1; x < 5; x++)
      std::cout << "(" << x << "," << y << ")= " << RECONoCJTauDMvsGENTauDM_->GetBinContent(x,y) << "  Row Average= " << RECONoCJTauDMvsGENTauDM_->GetBinContent(x,y) / sum << std::endl;
  }//for y
  std::cout << "\n\nRECONoCJTauDMvsGENTauDM by total" << std::endl;
  sumTotal = 0;
  for (int y = 1; y < 5; y++)
  { 
    for (int x = 1; x < 5; x++)
      sumTotal += RECONoCJTauDMvsGENTauDM_->GetBinContent(x,y);
  }//for y
  for (int y = 1; y < 5; y++)
  { 
    for (int x = 1; x < 5; x++)
      std::cout << "(" << x << "," << y << ")= " << RECONoCJTauDMvsGENTauDM_->GetBinContent(x,y) << "  Total Average= " << RECONoCJTauDMvsGENTauDM_->GetBinContent(x,y) / sumTotal << std::endl;
  }//for y


  std::cout << "\nRECOCJTauDMvsGENTauDMMedIso by column" << std::endl;
  for (int x = 1; x < 5; x++)
  {
    int sum = 0;
    for (int y = 1; y < 5; y++)
      sum += RECOCJTauDMvsGENTauDMMedIso_->GetBinContent(x,y);
    for (int y = 1; y < 5; y++)
      std::cout << "(" << x << "," << y << ")= " << RECOCJTauDMvsGENTauDMMedIso_->GetBinContent(x,y) << "  Column Average= " << RECOCJTauDMvsGENTauDMMedIso_->GetBinContent(x,y) / sum << std::endl;
  }//for x
  std::cout << "\n\nRECOCJTauDMvsGENTauDMMedIso by row" << std::endl;
  for (int y = 1; y < 5; y++)
  {
    int sum = 0;
    for (int x = 1; x < 5; x++)
      sum += RECOCJTauDMvsGENTauDMMedIso_->GetBinContent(x,y);
    for (int x = 1; x < 5; x++)
      std::cout << "(" << x << "," << y << ")= " << RECOCJTauDMvsGENTauDMMedIso_->GetBinContent(x,y) << "  Row Average= " << RECOCJTauDMvsGENTauDMMedIso_->GetBinContent(x,y) / sum << std::endl;
  }//for y
  std::cout << "\n\nRECOCJTauDMvsGENTauDMMedIso by total" << std::endl;
  sumTotal = 0;
  for (int y = 1; y < 5; y++)
  {
    for (int x = 1; x < 5; x++)
      sumTotal += RECOCJTauDMvsGENTauDMMedIso_->GetBinContent(x,y);
  }//for y
  for (int y = 1; y < 5; y++)
  { 
    for (int x = 1; x < 5; x++)
      std::cout << "(" << x << "," << y << ")= " << RECOCJTauDMvsGENTauDMMedIso_->GetBinContent(x,y) << "  Total Average= " << RECOCJTauDMvsGENTauDMMedIso_->GetBinContent(x,y) / sumTotal << std::endl;
  }//for y


  std::cout << "\n\n\n\nRECONoCJTauDMvsGENTauDMMedIso by column" << std::endl;
  for (int x = 1; x < 5; x++)
  {
    int sum = 0;
    for (int y = 1; y < 5; y++)
      sum += RECONoCJTauDMvsGENTauDMMedIso_->GetBinContent(x,y);
    for (int y = 1; y < 5; y++)
      std::cout << "(" << x << "," << y << ")= " << RECONoCJTauDMvsGENTauDMMedIso_->GetBinContent(x,y) << "  Column Average= " << RECONoCJTauDMvsGENTauDMMedIso_->GetBinContent(x,y) / sum << std::endl;
  }//for x
  std::cout << "\n\nRECONoCJTauDMvsGENTauDMMedIso by row" << std::endl;
  for (int y = 1; y < 5; y++)
  {
    int sum = 0;
    for (int x = 1; x < 5; x++)
      sum += RECONoCJTauDMvsGENTauDMMedIso_->GetBinContent(x,y);
    for (int x = 1; x < 5; x++)
      std::cout << "(" << x << "," << y << ")= " << RECONoCJTauDMvsGENTauDMMedIso_->GetBinContent(x,y) << "  Row Average= " << RECONoCJTauDMvsGENTauDMMedIso_->GetBinContent(x,y) / sum << std::endl;
  }//for y
  std::cout << "\n\nRECONoCJTauDMvsGENTauDMMedIso by total" << std::endl;
  sumTotal = 0;
  for (int y = 1; y < 5; y++)
  { 
    for (int x = 1; x < 5; x++)
      sumTotal += RECONoCJTauDMvsGENTauDMMedIso_->GetBinContent(x,y);
  }//for y
  for (int y = 1; y < 5; y++)
  { 
    for (int x = 1; x < 5; x++)
      std::cout << "(" << x << "," << y << ")= " << RECONoCJTauDMvsGENTauDMMedIso_->GetBinContent(x,y) << "  Total Average= " << RECONoCJTauDMvsGENTauDMMedIso_->GetBinContent(x,y) / sumTotal << std::endl;
  }//for y


  TCanvas MatchedLooseIsoRECOPtCanvas("MatchedLooseIsoRECOPt","",600,600);
  TCanvas MatchedMedIsoRECOPtCanvas("MatchedMedIsoRECOPt","",600,600);
  TCanvas MatchedTightIsoRECOPtCanvas("MatchedTightIsoRECOPt","",600,600);
  TCanvas MatchedDMFindRECOPtCanvas("MatchedDMFindRECOPt","",600,600);
  TCanvas MatchedRECOPtCanvas("MatchedRECOPt","",600,600);
  TCanvas FinalEffLooseIsoRECOPtCanvas("FinalEffLooseIsoRECOPt","",600,600);
  TCanvas FinalEffMedIsoRECOPtCanvas("FinalEffMedIsoRECOPt","",600,600);
  TCanvas FinalEffTightIsoRECOPtCanvas("FinalEffTightIsoRECOPt","",600,600);
  TCanvas FinalEffDMFindRECOPtCanvas("FinalEffDMFindRECOPt","",600,600);
  TCanvas FinalEffLooseIsoRECOPtOverDMCanvas("FinalEffLooseIsoRECOPtOverDM","",600,600);
  TCanvas FinalEffMedIsoRECOPtOverDMCanvas("FinalEffMedIsoRECOPtOverDM","",600,600);
  TCanvas FinalEffTightIsoRECOPtOverDMCanvas("FinalEffTightIsoRECOPtOverDM","",600,600);

  TCanvas MatchedLooseIsoCJPtCanvas("MatchedLooseIsoCJPt","",600,600);
  TCanvas MatchedMedIsoCJPtCanvas("MatchedMedIsoCJPt","",600,600);
  TCanvas MatchedTightIsoCJPtCanvas("MatchedTightIsoCJPt","",600,600);
  TCanvas MatchedDMFindCJPtCanvas("MatchedDMFindCJPt","",600,600);
  TCanvas MatchedCJPtCanvas("MatchedCJPt","",600,600);
  TCanvas FinalEffLooseIsoCJPtCanvas("FinalEffLooseIsoCJPt","",600,600);
  TCanvas FinalEffMedIsoCJPtCanvas("FinalEffMedIsoCJPt","",600,600);
  TCanvas FinalEffTightIsoCJPtCanvas("FinalEffTightIsoCJPt","",600,600);
  TCanvas FinalEffDMFindCJPtCanvas("FinalEffDMFindCJPt","",600,600);
  TCanvas FinalEffLooseIsoCJPtOverDMCanvas("FinalEffLooseIsoCJPtOverDM","",600,600);
  TCanvas FinalEffMedIsoCJPtOverDMCanvas("FinalEffMedIsoCJPtOverDM","",600,600);
  TCanvas FinalEffTightIsoCJPtOverDMCanvas("FinalEffTightIsoCJPtOverDM","",600,600);

  TCanvas MatchedLooseIsoRECOdRCanvas("MatchedLooseIsoRECOdR","",600,600);
  TCanvas MatchedMedIsoRECOdRCanvas("MatchedMedIsoRECOdR","",600,600);
  TCanvas MatchedTightIsoRECOdRCanvas("MatchedTightIsoRECOdR","",600,600);
  TCanvas MatchedDMFindRECOdRCanvas("MatchedDMFindRECOdR","",600,600);
  TCanvas MatchedRECOdRCanvas("MatchedRECOdR","",600,600);
  TCanvas FinalEffLooseIsoRECOdRCanvas("FinalEffLooseIsoRECOdR","",600,600);
  TCanvas FinalEffMedIsoRECOdRCanvas("FinalEffMedIsoRECOdR","",600,600);
  TCanvas FinalEffTightIsoRECOdRCanvas("FinalEffTightIsoRECOdR","",600,600);
  TCanvas FinalEffDMFindRECOdRCanvas("FinalEffDMFindRECOdR","",600,600);
  TCanvas FinalEffLooseIsoRECOdROverDMCanvas("FinalEffLooseIsoRECOdROverDM","",600,600);
  TCanvas FinalEffMedIsoRECOdROverDMCanvas("FinalEffMedIsoRECOdROverDM","",600,600);
  TCanvas FinalEffTightIsoRECOdROverDMCanvas("FinalEffTightIsoRECOdROverDM","",600,600);

  TCanvas MatchedLooseIsoCJdRCanvas("MatchedLooseIsoCJdR","",600,600);
  TCanvas MatchedMedIsoCJdRCanvas("MatchedMedIsoCJdR","",600,600);
  TCanvas MatchedTightIsoCJdRCanvas("MatchedTightIsoCJdR","",600,600);
  TCanvas MatchedDMFindCJdRCanvas("MatchedDMFindCJdR","",600,600);
  TCanvas MatchedCJdRCanvas("MatchedCJdR","",600,600);
  TCanvas FinalEffLooseIsoCJdRCanvas("FinalEffLooseIsoCJdR","",600,600);
  TCanvas FinalEffMedIsoCJdRCanvas("FinalEffMedIsoCJdR","",600,600);
  TCanvas FinalEffTightIsoCJdRCanvas("FinalEffTightIsoCJdR","",600,600);
  TCanvas FinalEffDMFindCJdRCanvas("FinalEffDMFindCJdR","",600,600);
  TCanvas FinalEffLooseIsoCJdROverDMCanvas("FinalEffLooseIsoCJdROverDM","",600,600);
  TCanvas FinalEffMedIsoCJdROverDMCanvas("FinalEffMedIsoCJdROverDM","",600,600);
  TCanvas FinalEffTightIsoCJdROverDMCanvas("FinalEffTightIsoCJdROverDM","",600,600);

  TCanvas MatchedLooseIsoRECOdRPtCutCanvas("MatchedLooseIsoRECOdRPtCut","",600,600);
  TCanvas MatchedMedIsoRECOdRPtCutCanvas("MatchedMedIsoRECOdRPtCut","",600,600);
  TCanvas MatchedTightIsoRECOdRPtCutCanvas("MatchedTightIsoRECOdRPtCut","",600,600);
  TCanvas MatchedDMFindRECOdRPtCutCanvas("MatchedDMFindRECOdRPtCut","",600,600);
  TCanvas MatchedRECOdRPtCutCanvas("MatchedRECOdRPtCut","",600,600);
  TCanvas FinalEffLooseIsoRECOdRPtCutCanvas("FinalEffLooseIsoRECOdRPtCut","",600,600);
  TCanvas FinalEffMedIsoRECOdRPtCutCanvas("FinalEffMedIsoRECOdRPtCut","",600,600);
  TCanvas FinalEffTightIsoRECOdRPtCutCanvas("FinalEffTightIsoRECOdRPtCut","",600,600);
  TCanvas FinalEffDMFindRECOdRPtCutCanvas("FinalEffDMFindRECOdRPtCut","",600,600);

  TCanvas MatchedLooseIsoCJdRPtCutCanvas("MatchedLooseIsoCJdRPtCut","",600,600);
  TCanvas MatchedMedIsoCJdRPtCutCanvas("MatchedMedIsoCJdRPtCut","",600,600);
  TCanvas MatchedTightIsoCJdRPtCutCanvas("MatchedTightIsoCJdRPtCut","",600,600);
  TCanvas MatchedDMFindCJdRPtCutCanvas("MatchedDMFindCJdRPtCut","",600,600);
  TCanvas MatchedCJdRPtCutCanvas("MatchedCJdRPtCut","",600,600);
  TCanvas FinalEffLooseIsoCJdRPtCutCanvas("FinalEffLooseIsoCJdRPtCut","",600,600);
  TCanvas FinalEffMedIsoCJdRPtCutCanvas("FinalEffMedIsoCJdRPtCut","",600,600);
  TCanvas FinalEffTightIsoCJdRPtCutCanvas("FinalEffTightIsoCJdRPtCut","",600,600);
  TCanvas FinalEffDMFindCJdRPtCutCanvas("FinalEffDMFindCJdRPtCut","",600,600);

  TCanvas MatchedLooseIsoRECOPtGenCanvas("MatchedLooseIsoRECOPtGen","",600,600);
  TCanvas MatchedMedIsoRECOPtGenCanvas("MatchedMedIsoRECOPtGen","",600,600);
  TCanvas MatchedTightIsoRECOPtGenCanvas("MatchedTightIsoRECOPtGen","",600,600);
  TCanvas MatchedDMFindRECOPtGenCanvas("MatchedDMFindRECOPtGen","",600,600);
  TCanvas MatchedRECOPtGenCanvas("MatchedRECOPtGen","",600,600);
  TCanvas FinalEffLooseIsoRECOPtGenCanvas("FinalEffLooseIsoRECOPtGen","",600,600);
  TCanvas FinalEffMedIsoRECOPtGenCanvas("FinalEffMedIsoRECOPtGen","",600,600);
  TCanvas FinalEffTightIsoRECOPtGenCanvas("FinalEffTightIsoRECOPtGen","",600,600);
  TCanvas FinalEffDMFindRECOPtGenCanvas("FinalEffDMFindRECOPtGen","",600,600);
  TCanvas FinalEffLooseIsoRECOPtGenOverDMCanvas("FinalEffLooseIsoRECOPtGenOverDM","",600,600);
  TCanvas FinalEffMedIsoRECOPtGenOverDMCanvas("FinalEffMedIsoRECOPtGenOverDM","",600,600);
  TCanvas FinalEffTightIsoRECOPtGenOverDMCanvas("FinalEffTightIsoRECOPtGenOverDM","",600,600);

  TCanvas MatchedLooseIsoCJPtGenCanvas("MatchedLooseIsoCJPtGen","",600,600);
  TCanvas MatchedMedIsoCJPtGenCanvas("MatchedMedIsoCJPtGen","",600,600);
  TCanvas MatchedTightIsoCJPtGenCanvas("MatchedTightIsoCJPtGen","",600,600);
  TCanvas MatchedDMFindCJPtGenCanvas("MatchedDMFindCJPtGen","",600,600);
  TCanvas MatchedCJPtGenCanvas("MatchedCJPtGen","",600,600);
  TCanvas FinalEffLooseIsoCJPtGenCanvas("FinalEffLooseIsoCJPtGen","",600,600);
  TCanvas FinalEffMedIsoCJPtGenCanvas("FinalEffMedIsoCJPtGen","",600,600);
  TCanvas FinalEffTightIsoCJPtGenCanvas("FinalEffTightIsoCJPtGen","",600,600);
  TCanvas FinalEffDMFindCJPtGenCanvas("FinalEffDMFindCJPtGen","",600,600);
  TCanvas FinalEffLooseIsoCJPtGenOverDMCanvas("FinalEffLooseIsoCJPtGenOverDM","",600,600);
  TCanvas FinalEffMedIsoCJPtGenOverDMCanvas("FinalEffMedIsoCJPtGenOverDM","",600,600);
  TCanvas FinalEffTightIsoCJPtGenOverDMCanvas("FinalEffTightIsoCJPtGenOverDM","",600,600);

  TCanvas OneProngDMCJPtCanvas("OneProngDMCJPt","",600,600);
  TCanvas OneProngOnePizDMCJPtCanvas("OneProngOnePizDMCJPt","",600,600);
  TCanvas OneProngTwoPizDMCJPtCanvas("OneProngTwoPizDMCJPt","",600,600);
  TCanvas ThreeProngDMCJPtCanvas("ThreeProngDMCJPt","",600,600);
  TCanvas OneProngDMRECOPtCanvas("OneProngDMRECOPt","",600,600);
  TCanvas OneProngOnePizDMRECOPtCanvas("OneProngOnePizDMRECOPt","",600,600);
  TCanvas OneProngTwoPizDMRECOPtCanvas("OneProngTwoPizDMRECOPt","",600,600);
  TCanvas ThreeProngDMRECOPtCanvas("ThreeProngDMRECOPt","",600,600);
  TCanvas MatchedOneProngCJPtCanvas("MatchedOneProngCJPtCanvas","",600,600);
  TCanvas MatchedOneProngOnePizCJPtCanvas("MatchedOneProngOnePizCJPtCanvas","",600,600);
  TCanvas MatchedOneProngTwoPizCJPtCanvas("MatchedOneProngTwoPizCJPtCanvas","",600,600);
  TCanvas MatchedThreeProngCJPtCanvas("MatchedThreeProngCJPtCanvas","",600,600);
  TCanvas MatchedOneProngRECOPtCanvas("MatchedOneProngRECOPtCanvas","",600,600);
  TCanvas MatchedOneProngOnePizRECOPtCanvas("MatchedOneProngOnePizRECOPtCanvas","",600,600);
  TCanvas MatchedOneProngTwoPizRECOPtCanvas("MatchedOneProngTwoPizRECOPtCanvas","",600,600);
  TCanvas MatchedThreeProngRECOPtCanvas("MatchedThreeProngRECOPtCanvas","",600,600);
  TCanvas FinalOneProngDMCJPtCanvas("FinalOneProngDMCJPt","",600,600);
  TCanvas FinalOneProngOnePizDMCJPtCanvas("FinalOneProngOnePizDMCJPt","",600,600);
  TCanvas FinalOneProngTwoPizDMCJPtCanvas("FinalOneProngTwoPizDMCJPt","",600,600);
  TCanvas FinalThreeProngDMCJPtCanvas("FinalThreeProngDMCJPt","",600,600);
  TCanvas FinalOneProngDMRECOPtCanvas("FinalOneProngDMRECOPt","",600,600);
  TCanvas FinalOneProngOnePizDMRECOPtCanvas("FinalOneProngOnePizDMRECOPt","",600,600);
  TCanvas FinalOneProngTwoPizDMRECOPtCanvas("FinalOneProngTwoPizDMRECOPt","",600,600);
  TCanvas FinalThreeProngDMRECOPtCanvas("FinalThreeProngDMRECOPt","",600,600);
  TCanvas FinalEffOneProngDMSAMEPtCanvas("FinalEffOneProngDMSAMEPt","",600,600);
  TCanvas FinalEffOneProngOnePizDMSAMEPtCanvas("FinalEffOneProngOnePizDMSAMEPt","",600,600);
  TCanvas FinalEffOneProngTwoPizDMSAMEPtCanvas("FinalEffOneProngTwoPizDMSAMEPt","",600,600);
  TCanvas FinalEffThreeProngDMSAMEPtCanvas("FinalEffThreeProngDMSAMEPt","",600,600);




  TCanvas FinalEffLooseIsoSAMEPtCanvas("FinalEffLooseIsoSAMEPt","",600,600);
  TCanvas FinalEffMedIsoSAMEPtCanvas("FinalEffMedIsoSAMEPt","",600,600);
  TCanvas FinalEffTightIsoSAMEPtCanvas("FinalEffTightIsoSAMEPt","",600,600);
  TCanvas FinalEffDMFindSAMEPtCanvas("FinalEffDMFindSAMEPt","",600,600);
  TCanvas FinalEffLooseIsoSAMEdRCanvas("FinalEffLooseIsoSAMEdR","",600,600);
  TCanvas FinalEffMedIsoSAMEdRCanvas("FinalEffMedIsoSAMEdR","",600,600);
  TCanvas FinalEffTightIsoSAMEdRCanvas("FinalEffTightIsoSAMEdR","",600,600);
  TCanvas FinalEffDMFindSAMEdRCanvas("FinalEffDMFindSAMEdR","",600,600);
  TCanvas FinalEffLooseIsoSAMEdRPtCutCanvas("FinalEffLooseIsoSAMEdRPtCut","",600,600);
  TCanvas FinalEffMedIsoSAMEdRPtCutCanvas("FinalEffMedIsoSAMEdRPtCut","",600,600);
  TCanvas FinalEffTightIsoSAMEdRPtCutCanvas("FinalEffTightIsoSAMEdRPtCut","",600,600);
  TCanvas FinalEffDMFindSAMEdRPtCutCanvas("FinalEffDMFindSAMEdRPtCut","",600,600);
  TCanvas FinalEffLooseIsoSAMEPtGenCanvas("FinalEffLooseIsoSAMEPtGen","",600,600);
  TCanvas FinalEffMedIsoSAMEPtGenCanvas("FinalEffMedIsoSAMEPtGen","",600,600);
  TCanvas FinalEffTightIsoSAMEPtGenCanvas("FinalEffTightIsoSAMEPtGen","",600,600);
  TCanvas FinalEffDMFindSAMEPtGenCanvas("FinalEffDMFindSAMEPtGen","",600,600);

  TCanvas FinalEffLooseIsoSAMEPtOverDMCanvas("FinalEffLooseIsoSAMEPtOverDM","",600,600);
  TCanvas FinalEffMedIsoSAMEPtOverDMCanvas("FinalEffMedIsoSAMEPtOverDM","",600,600);
  TCanvas FinalEffTightIsoSAMEPtOverDMCanvas("FinalEffTightIsoSAMEPtOverDM","",600,600);
  TCanvas FinalEffLooseIsoSAMEdROverDMCanvas("FinalEffLooseIsoSAMEdROverDM","",600,600);
  TCanvas FinalEffMedIsoSAMEdROverDMCanvas("FinalEffMedIsoSAMEdROverDM","",600,600);
  TCanvas FinalEffTightIsoSAMEdROverDMCanvas("FinalEffTightIsoSAMEdROverDM","",600,600);
  TCanvas FinalEffLooseIsoSAMEPtGenOverDMCanvas("FinalEffLooseIsoSAMEPtGenOverDM","",600,600);
  TCanvas FinalEffMedIsoSAMEPtGenOverDMCanvas("FinalEffMedIsoSAMEPtGenOverDM","",600,600);
  TCanvas FinalEffTightIsoSAMEPtGenOverDMCanvas("FinalEffTightIsoSAMEPtGenOverDM","",600,600);
/*
  FinalEffLooseIsoRECOPt_->Divide(MatchedLooseIsoRECOPt_, MatchedRECOPt_);
  FinalEffMedIsoRECOPt_->Divide(MatchedMedIsoRECOPt_,     MatchedRECOPt_);
  FinalEffTightIsoRECOPt_->Divide(MatchedTightIsoRECOPt_, MatchedRECOPt_);
  FinalEffDMFindRECOPt_->Divide(MatchedDMFindRECOPt_,     MatchedRECOPt_);
  FinalEffLooseIsoCJPt_->Divide(MatchedLooseIsoCJPt_, MatchedCJPt_);
  FinalEffMedIsoCJPt_->Divide(MatchedMedIsoCJPt_,     MatchedCJPt_);
  FinalEffTightIsoCJPt_->Divide(MatchedTightIsoCJPt_, MatchedCJPt_);
  FinalEffDMFindCJPt_->Divide(MatchedDMFindCJPt_,     MatchedCJPt_);
  FinalEffLooseIsoRECOPtOverDM_->Divide(MatchedLooseIsoRECOPt_, MatchedDMFindRECOPt_);
  FinalEffMedIsoRECOPtOverDM_->Divide(MatchedMedIsoRECOPt_,     MatchedDMFindRECOPt_);
  FinalEffTightIsoRECOPtOverDM_->Divide(MatchedTightIsoRECOPt_, MatchedDMFindRECOPt_);
  FinalEffLooseIsoCJPtOverDM_->Divide(MatchedLooseIsoCJPt_, MatchedDMFindCJPt_);
  FinalEffMedIsoCJPtOverDM_->Divide(MatchedMedIsoCJPt_,     MatchedDMFindCJPt_);
  FinalEffTightIsoCJPtOverDM_->Divide(MatchedTightIsoCJPt_, MatchedDMFindCJPt_);

  FinalEffLooseIsoRECOdR_->Divide(MatchedLooseIsoRECOdR_, MatchedRECOdR_);
  FinalEffMedIsoRECOdR_->Divide(MatchedMedIsoRECOdR_,     MatchedRECOdR_);
  FinalEffTightIsoRECOdR_->Divide(MatchedTightIsoRECOdR_, MatchedRECOdR_);
  FinalEffDMFindRECOdR_->Divide(MatchedDMFindRECOdR_,     MatchedRECOdR_);
  FinalEffLooseIsoCJdR_->Divide(MatchedLooseIsoCJdR_, MatchedCJdR_);
  FinalEffMedIsoCJdR_->Divide(MatchedMedIsoCJdR_,     MatchedCJdR_);
  FinalEffTightIsoCJdR_->Divide(MatchedTightIsoCJdR_, MatchedCJdR_);
  FinalEffDMFindCJdR_->Divide(MatchedDMFindCJdR_,     MatchedCJdR_);
  FinalEffLooseIsoRECOdROverDM_->Divide(MatchedLooseIsoRECOdR_, MatchedDMFindRECOdR_);
  FinalEffMedIsoRECOdROverDM_->Divide(MatchedMedIsoRECOdR_,     MatchedDMFindRECOdR_);
  FinalEffTightIsoRECOdROverDM_->Divide(MatchedTightIsoRECOdR_, MatchedDMFindRECOdR_);
  FinalEffLooseIsoCJdROverDM_->Divide(MatchedLooseIsoCJdR_, MatchedDMFindCJdR_);
  FinalEffMedIsoCJdROverDM_->Divide(MatchedMedIsoCJdR_,     MatchedDMFindCJdR_);
  FinalEffTightIsoCJdROverDM_->Divide(MatchedTightIsoCJdR_, MatchedDMFindCJdR_);

  FinalEffLooseIsoRECOdRPtCut_->Divide(MatchedLooseIsoRECOdRPtCut_, MatchedRECOdRPtCut_);
  FinalEffMedIsoRECOdRPtCut_->Divide(MatchedMedIsoRECOdRPtCut_,     MatchedRECOdRPtCut_);
  FinalEffTightIsoRECOdRPtCut_->Divide(MatchedTightIsoRECOdRPtCut_, MatchedRECOdRPtCut_);
  FinalEffDMFindRECOdRPtCut_->Divide(MatchedDMFindRECOdRPtCut_,     MatchedRECOdRPtCut_);
  FinalEffLooseIsoCJdRPtCut_->Divide(MatchedLooseIsoCJdRPtCut_, MatchedCJdRPtCut_);
  FinalEffMedIsoCJdRPtCut_->Divide(MatchedMedIsoCJdRPtCut_,     MatchedCJdRPtCut_);
  FinalEffTightIsoCJdRPtCut_->Divide(MatchedTightIsoCJdRPtCut_, MatchedCJdRPtCut_);
  FinalEffDMFindCJdRPtCut_->Divide(MatchedDMFindCJdRPtCut_,     MatchedCJdRPtCut_);

  FinalEffLooseIsoRECOPtGen_->Divide(MatchedLooseIsoRECOPtGen_, MatchedRECOPtGen_);
  FinalEffMedIsoRECOPtGen_->Divide(MatchedMedIsoRECOPtGen_,     MatchedRECOPtGen_);
  FinalEffTightIsoRECOPtGen_->Divide(MatchedTightIsoRECOPtGen_, MatchedRECOPtGen_);
  FinalEffDMFindRECOPtGen_->Divide(MatchedDMFindRECOPtGen_,     MatchedRECOPtGen_);
  FinalEffLooseIsoCJPtGen_->Divide(MatchedLooseIsoCJPtGen_, MatchedCJPtGen_);
  FinalEffMedIsoCJPtGen_->Divide(MatchedMedIsoCJPtGen_,     MatchedCJPtGen_);
  FinalEffTightIsoCJPtGen_->Divide(MatchedTightIsoCJPtGen_, MatchedCJPtGen_);
  FinalEffDMFindCJPtGen_->Divide(MatchedDMFindCJPtGen_,     MatchedCJPtGen_);
  FinalEffLooseIsoRECOPtGenOverDM_->Divide(MatchedLooseIsoRECOPtGen_, MatchedDMFindRECOPtGen_);
  FinalEffMedIsoRECOPtGenOverDM_->Divide(MatchedMedIsoRECOPtGen_,     MatchedDMFindRECOPtGen_);
  FinalEffTightIsoRECOPtGenOverDM_->Divide(MatchedTightIsoRECOPtGen_, MatchedDMFindRECOPtGen_);
  FinalEffLooseIsoCJPtGenOverDM_->Divide(MatchedLooseIsoCJPtGen_, MatchedDMFindCJPtGen_);
  FinalEffMedIsoCJPtGenOverDM_->Divide(MatchedMedIsoCJPtGen_,     MatchedDMFindCJPtGen_);
  FinalEffTightIsoCJPtGenOverDM_->Divide(MatchedTightIsoCJPtGen_, MatchedDMFindCJPtGen_);

  FinalOneProngDMCJPt_->Divide(OneProngDMCJPt_, MatchedOneProngCJPt_);
  FinalOneProngOnePizDMCJPt_->Divide(OneProngOnePizDMCJPt_, MatchedOneProngOnePizCJPt_);
  FinalOneProngTwoPizDMCJPt_->Divide(OneProngTwoPizDMCJPt_, MatchedOneProngTwoPizCJPt_);
  FinalThreeProngDMCJPt_->Divide(ThreeProngDMCJPt_, MatchedThreeProngCJPt_);
  FinalOneProngDMRECOPt_->Divide(OneProngDMRECOPt_, MatchedOneProngRECOPt_);
  FinalOneProngOnePizDMRECOPt_->Divide(OneProngOnePizDMRECOPt_, MatchedOneProngOnePizRECOPt_);
  FinalOneProngTwoPizDMRECOPt_->Divide(OneProngTwoPizDMRECOPt_, MatchedOneProngTwoPizRECOPt_);
  FinalThreeProngDMRECOPt_->Divide(ThreeProngDMRECOPt_, MatchedThreeProngRECOPt_);
*/
std::cout << "<----------------Declared Canvases-------------->" << std::endl;

  //Format the 1D plots and Draw (canvas, hist, grid, log y, log z, color, size, style, xAxisTitle, xTitleSize, xLabelSize, xTitleOffSet, yAxisTitle, yTitleSize, yLabelSize, yTitleOffset)
  VariousFunctions::formatAndDrawCanvasAndHist1D(NEventsCanvas, NEvents_,
	 1, 0, 0, kBlack, 7, 20, "", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(GenMatchTypeCanvas, GenMatchType_,
         1, 0, 0, kBlack, 7, 20, "", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(GenMatchPDGIDCanvas, GenMatchPDGID_,
         1, 0, 0, kBlack, 7, 20, "PDG ID #", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(GenMatchPDGIDwithLargeEtaCanvas, GenMatchPDGIDwithLargeEta_,
         1, 0, 0, kBlack, 7, 20, "PDG ID # if |#eta| > 1.5", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NMuRemovedCanvas, NMuRemoved_,
	 1, 0, 0, kBlack, 7, 20, "N #mu Removed", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauMuTauHaddRCanvas, TauMuTauHaddR_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(#tau_{mu} + #tau_{H})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MassDiLepGENCanvas, MassDiLepGEN_,
         1, 0, 0, kBlack, 7, 20, "Gen m_{di-lept}", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MassDiLepRECOCanvas, MassDiLepRECO_,
         1, 0, 0, kBlack, 7, 20, "RECO m_{di-lep}", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(Mu3LargestPtMatchCanvas, Mu3LargestPtMatch_,
         1, 0, 0, kBlack, 7, 20, "Largest p_{T} (Mu3) matched to Tau", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(Mu3PtMatchCanvas, Mu3PtMatch_,
         1, 0, 0, kBlack, 7, 20, "p_{T} (Mu3) matched to Tau", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTauMuPtLeptonCanvas, MatchedTauMuPtLepton_,
         1, 0, 0, kBlack, 7, 20, "Lepton: p_{T}(#tau_{mu})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTauHadPtLeptonCanvas, MatchedTauHadPtLepton_,
         1, 0, 0, kBlack, 7, 20, "Lepton: p_{T}(#tau_{H})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTauHadEtaLeptonCanvas, MatchedTauHadEtaLepton_,
         1, 0, 0, kBlack, 7, 20, "Lepton: #eta(#tau_{H})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedHTLeptonCanvas, MatchedHTLepton_,
         1, 0, 0, kBlack, 7, 20, "Lepton: H_{T} of ak4PFJet with p_{T} > 20", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedBDiscCSVLeptonCanvas, MatchedBDiscCSVLepton_,
         1, 0, 0, kBlack, 7, 20, "Lepton: CSV B Discriminant Value", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTauMuPtJetCanvas, MatchedTauMuPtJet_,
         1, 0, 0, kBlack, 7, 20, "Jet: p_{T}(#tau_{mu})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTauHadPtJetCanvas, MatchedTauHadPtJet_,
         1, 0, 0, kBlack, 7, 20, "Jet: p_{T}(#tau_{H})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTauHadEtaJetCanvas, MatchedTauHadEtaJet_,
         1, 0, 0, kBlack, 7, 20, "Jet: #eta(#tau_{H})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedHTJetCanvas, MatchedHTJet_,
         1, 0, 0, kBlack, 7, 20, "Jet: H_{T} of ak4PFJet with p_{T} > 20", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedBDiscCSVJetCanvas, MatchedBDiscCSVJet_,
         1, 0, 0, kBlack, 7, 20, "Jet: CSV B Discriminant Value", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NEventsCutsCanvas, NEventsCuts_,
         1, 0, 0, kBlack, 7, 20, "", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NConstituentsCJCanvas, NConstituentsCJ_,
         1, 0, 0, kBlack, 7, 20, "p_{T}(#mu_{1}#mu_{2}) / M(#mu_{1}#mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(PtOverMassMu1Mu2Canvas, PtOverMassMu1Mu2_,
         1, 0, 0, kBlack, 7, 20, "p_{T}(#mu_{1}#mu_{2}) / Mass(#mu_{1}#mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(PtMu1Mu2Canvas, PtMu1Mu2_,
         1, 0, 0, kBlack, 7, 20, "p_{T}(#mu_{1}#mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(NTauDecayModeGENCanvas, NTauDecayModeGEN_,
         1, 0, 0, kBlack, 7, 20, "TauDecayMode GEN ", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NTauDecayModeRECOCJCanvas, NTauDecayModeRECOCJ_,
	 1, 0, 0, kBlack, 7, 20, "TauDecayMode RECOCJ", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NTauDecayModeRECONoCJCanvas, NTauDecayModeRECONoCJ_,
         1, 0, 0, kBlack, 7, 20, "TauDecayMode RECONoCJ", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NTauDecayModeRECOCJIfGen3Canvas, NTauDecayModeRECOCJIfGen3_,
         1, 0, 0, kBlack, 7, 20, "TauDecayMode RECO CJ if GEN=3", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NTauDecayModeRECONoCJIfGen3Canvas, NTauDecayModeRECONoCJIfGen3_,
         1, 0, 0, kBlack, 7, 20, "TauDecayMode RECO No CJ if GEN=3", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist2D(NTausRECOvsCLEANJETSCanvas, NTausRECOvsCLEANJETS_,
	 1, 0, 0, kBlack, 7, 20, "nRECO #tau's", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(GenDiTaudRvsCJDiTaudRCanvas, GenDiTaudRvsCJDiTaudR_,
	 1, 0, 0, kBlack, 7, 20, "gen #DeltaR(#tau#tau)", .04, .04, 1.1, "CJ #DeltaR(#tau#tau)", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(TauHadnConstvsPtCanvas, TauHadnConstvsPt_,
         1, 0, 0, kBlack, 7, 20, "# Constituents(#tau_{H})", .04, .04, 1.1, "Pt(#tau_{H})", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(RECOPtvsGENVisiblePtCJCanvas, RECOPtvsGENVisiblePtCJ_,
         1, 0, 0, kBlack, 7, 20, "RECO Pt(#tau)", .04, .04, 1.1, "GEN Visible Pt(#tau)", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(RECOPtvsGENVisiblePtRECOCanvas, RECOPtvsGENVisiblePtRECO_,
         1, 0, 0, kBlack, 7, 20, "RECO Pt(#tau)", .04, .04, 1.1, "GEN Visible Pt(#tau)", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(RECOCJTauDMvsGENTauDMCanvas, RECOCJTauDMvsGENTauDM_,
         1, 0, 0, kBlack, 7, 20, "RECO CJ Tau DM", .04, .04, 1.1, "GEN Tau DM", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(RECONoCJTauDMvsGENTauDMCanvas, RECONoCJTauDMvsGENTauDM_,
	 1, 0, 0, kBlack, 7, 20, "RECO No CJ Tau DM", .04, .04, 1.1, "GEN Tau DM", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(RECOCJTauDMvsGENTauDMMedIsoCanvas, RECOCJTauDMvsGENTauDMMedIso_,
         1, 0, 0, kBlack, 7, 20, "RECO CJ Tau DM", .04, .04, 1.1, "GEN Tau DM", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(RECONoCJTauDMvsGENTauDMMedIsoCanvas, RECONoCJTauDMvsGENTauDMMedIso_,
         1, 0, 0, kBlack, 7, 20, "RECO No CJ Tau DM", .04, .04, 1.1, "GEN Tau DM", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(GenVisiblePtvsRecoCJPtIfCorrDMCanvas, GenVisiblePtvsRecoCJPtIfCorrDM_,
         1, 0, 0, kBlack, 7, 20, "RECO CJ Tau p_{T} if RECO DM is correct", .04, .04, 1.1, "GEN Visible p_{T}(#tau)", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(GenVisiblePtvsRecoNoCJPtIfCorrDMCanvas, GenVisiblePtvsRecoNoCJPtIfCorrDM_,
         1, 0, 0, kBlack, 7, 20, "RECO No CJ Tau p_{T} if RECO DM is correct", .04, .04, 1.1, "GEN Visible p_{T}(#tau)", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(RECOPtvsGENdRCJCanvas, RECOPtvsGENdRCJ_,
         1, 0, 0, kBlack, 7, 20, "RECO CJ p_{T}", .04, .04, 1.1, "GEN #DeltaR(#tau#tau)", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(RECOPtvsGENdRRECOCanvas, RECOPtvsGENdRRECO_,
         1, 0, 0, kBlack, 7, 20, "RECO No CJ p_{T}", .04, .04, 1.1, "GEN #DeltaR(#tau#tau)", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(RECOPtvsGENdRCJMedIsoCanvas, RECOPtvsGENdRCJMedIso_,
         1, 0, 0, kBlack, 7, 20, "RECO CJ p_{T} + MedIso", .04, .04, 1.1, "GEN #DeltaR(#tau#tau)", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(RECOPtvsGENdRRECOMedIsoCanvas, RECOPtvsGENdRRECOMedIso_,
         1, 0, 0, kBlack, 7, 20, "RECO No CJ p_{T} + MedIso", .04, .04, 1.1, "GEN #DeltaR(#tau#tau)", .04, .04, 1.6, "", .04, .04, 1.0);
 VariousFunctions::formatAndDrawCanvasAndHist2D(MatchedTauHadEtavsTauHadPtCanvas, MatchedTauHadEtavsTauHadPt_,
        1, 0, 0, kBlack, 7, 20, "GM #eta(#tau_{H})", .04, .04, 1.1, "GM p_{T}(#tau_{H})", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(BDiscCSVvsRawIsoCanvas, BDiscCSVvsRawIso_,
         1, 0, 0, kBlack, 7, 20, "pfBTag CSV Value", .04, .04, 1.1, "Raw Iso Value", .04, .04, 1.6, "", .04, .04, 1.0);



  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoRECOPtCanvas, MatchedLooseIsoRECOPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(RECO Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoRECOPtCanvas, MatchedMedIsoRECOPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(RECO Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoRECOPtCanvas, MatchedTightIsoRECOPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(RECO Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedDMFindRECOPtCanvas, MatchedDMFindRECOPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(RECO Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedRECOPtCanvas, MatchedRECOPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);


  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoRECOPtCanvas, FinalEffLooseIsoRECOPt_,
	 1, 1, 1, kBlack, 1, 20, "Pt(RECO Matched + Loose Iso + DM) / Pt(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoRECOPtCanvas, FinalEffMedIsoRECOPt_,
	 1, 1, 1, kBlack, 1, 20, "Pt(RECO Matched + Med Iso + DM) / Pt(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoRECOPtCanvas, FinalEffTightIsoRECOPt_,
	 1, 1, 1, kBlack, 1, 20, "Pt(RECO Matched + Tight Iso + DM) / Pt(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffDMFindRECOPtCanvas, FinalEffDMFindRECOPt_,
	 1, 1, 1, kBlack, 1, 20, "Pt(RECO Matched + DecayModeFinding) / Pt(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoRECOPtOverDMCanvas, FinalEffLooseIsoRECOPtOverDM_,
         1, 1, 1, kBlack, 1, 20, "Pt(RECO Matched + Loose Iso + DM) / Pt(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoRECOPtOverDMCanvas, FinalEffMedIsoRECOPtOverDM_,
         1, 1, 1, kBlack, 1, 20, "Pt(RECO Matched + Med Iso + DM) / Pt(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoRECOPtOverDMCanvas, FinalEffTightIsoRECOPtOverDM_,
         1, 1, 1, kBlack, 1, 20, "Pt(RECO Matched + Tight Iso + DM) / Pt(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);


  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoCJPtCanvas, MatchedLooseIsoCJPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoCJPtCanvas, MatchedMedIsoCJPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoCJPtCanvas, MatchedTightIsoCJPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedDMFindCJPtCanvas, MatchedDMFindCJPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedCJPtCanvas, MatchedCJPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoCJPtCanvas, FinalEffLooseIsoCJPt_,
	 1, 1, 1, kBlack, 1, 20, "Pt(CleanJets Matched + Loose Iso + DM) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoCJPtCanvas, FinalEffMedIsoCJPt_,
	 1, 1, 1, kBlack, 1, 20, "Pt(CleanJets Matched + Med Iso + DM) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoCJPtCanvas, FinalEffTightIsoCJPt_,
	 1, 1, 1, kBlack, 1, 20, "Pt(CleanJets Matched + Tight Iso + DM) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffDMFindCJPtCanvas, FinalEffDMFindCJPt_,
	 1, 1, 1, kBlack, 1, 20, "Pt(CleanJets Matched + DecayModeFinding) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoCJPtOverDMCanvas, FinalEffLooseIsoCJPtOverDM_,
         1, 1, 1, kBlack, 1, 20, "Pt(CleanJets Matched + Loose Iso + DM) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoCJPtOverDMCanvas, FinalEffMedIsoCJPtOverDM_,
         1, 1, 1, kBlack, 1, 20, "Pt(CleanJets Matched + Med Iso + DM) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoCJPtOverDMCanvas, FinalEffTightIsoCJPtOverDM_,
         1, 1, 1, kBlack, 1, 20, "Pt(CleanJets Matched + Tight Iso + DM) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoRECOdRPtCutCanvas, MatchedLooseIsoRECOdRPtCut_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(RECO Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoRECOdRPtCutCanvas, MatchedMedIsoRECOdRPtCut_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(RECO Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoRECOdRPtCutCanvas, MatchedTightIsoRECOdRPtCut_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(RECO Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedDMFindRECOdRPtCutCanvas, MatchedDMFindRECOdRPtCut_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(RECO Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedRECOdRPtCutCanvas, MatchedRECOdRPtCut_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoRECOdRPtCutCanvas, FinalEffLooseIsoRECOdRPtCut_,
	 1, 1, 1, kBlack, 1, 20, "#DeltaR(RECO Matched + Loose Iso + DM) / #DeltaR(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoRECOdRPtCutCanvas, FinalEffMedIsoRECOdRPtCut_,
	 1, 1, 1, kBlack, 1, 20, "#DeltaR(RECO Matched + Med Iso + DM) / #DeltaR(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoRECOdRPtCutCanvas, FinalEffTightIsoRECOdRPtCut_,
	 1, 1, 1, kBlack, 1, 20, "#DeltaR(RECO Matched + Tight Iso + DM) / #DeltaR(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffDMFindRECOdRPtCutCanvas, FinalEffDMFindRECOdRPtCut_,
	 1, 1, 1, kBlack, 1, 20, "#DeltaR(RECO Matched + DecayModeFinding) / #DeltaR(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoRECOdRCanvas, MatchedLooseIsoRECOdR_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(RECO Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoRECOdRCanvas, MatchedMedIsoRECOdR_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(RECO Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoRECOdRCanvas, MatchedTightIsoRECOdR_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(RECO Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedDMFindRECOdRCanvas, MatchedDMFindRECOdR_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(RECO Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedRECOdRCanvas, MatchedRECOdR_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoRECOdRCanvas, FinalEffLooseIsoRECOdR_,
         1, 1, 1, kBlack, 1, 20, "#DeltaR(RECO Matched + Loose Iso + DM) / #DeltaR(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoRECOdRCanvas, FinalEffMedIsoRECOdR_,
         1, 1, 1, kBlack, 1, 20, "#DeltaR(RECO Matched + Med Iso + DM) / #DeltaR(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoRECOdRCanvas, FinalEffTightIsoRECOdR_,
         1, 1, 1, kBlack, 1, 20, "#DeltaR(RECO Matched + Tight Iso + DM) / #DeltaR(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffDMFindRECOdRCanvas, FinalEffDMFindRECOdR_,
         1, 1, 1, kBlack, 1, 20, "#DeltaR(RECO Matched + DecayModeFinding) / #DeltaR(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoRECOdROverDMCanvas, FinalEffLooseIsoRECOdROverDM_,
         1, 1, 1, kBlack, 1, 20, "#DeltaR(RECO Matched + Loose Iso + DM) / #DeltaR(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoRECOdROverDMCanvas, FinalEffMedIsoRECOdROverDM_,
         1, 1, 1, kBlack, 1, 20, "#DeltaR(RECO Matched + Med Iso + DM) / #DeltaR(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoRECOdROverDMCanvas, FinalEffTightIsoRECOdROverDM_,
         1, 1, 1, kBlack, 1, 20, "#DeltaR(RECO Matched + Tight Iso + DM) / #DeltaR(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoCJdRCanvas, MatchedLooseIsoCJdR_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoCJdRCanvas, MatchedMedIsoCJdR_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoCJdRCanvas, MatchedTightIsoCJdR_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedDMFindCJdRCanvas, MatchedDMFindCJdR_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedCJdRCanvas, MatchedCJdR_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoCJdRCanvas, FinalEffLooseIsoCJdR_,
         1, 1, 1, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Loose Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoCJdRCanvas, FinalEffMedIsoCJdR_,
         1, 1, 1, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Med Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoCJdRCanvas, FinalEffTightIsoCJdR_,
         1, 1, 1, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Tight Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffDMFindCJdRCanvas, FinalEffDMFindCJdR_,
         1, 1, 1, kBlack, 1, 20, "#DeltaR(CleanJets Matched + DecayModeFinding) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoCJdROverDMCanvas, FinalEffLooseIsoCJdROverDM_,
         1, 1, 1, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Loose Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoCJdROverDMCanvas, FinalEffMedIsoCJdROverDM_,
         1, 1, 1, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Med Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoCJdROverDMCanvas, FinalEffTightIsoCJdROverDM_,
         1, 1, 1, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Tight Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoCJdRPtCutCanvas, MatchedLooseIsoCJdRPtCut_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoCJdRPtCutCanvas, MatchedMedIsoCJdRPtCut_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoCJdRPtCutCanvas, MatchedTightIsoCJdRPtCut_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedDMFindCJdRPtCutCanvas, MatchedDMFindCJdRPtCut_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedCJdRPtCutCanvas, MatchedCJdRPtCut_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoCJdRPtCutCanvas, FinalEffLooseIsoCJdRPtCut_,
	 1, 1, 1, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Loose Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoCJdRPtCutCanvas, FinalEffMedIsoCJdRPtCut_,
	 1, 1, 1, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Med Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoCJdRPtCutCanvas, FinalEffTightIsoCJdRPtCut_,
	 1, 1, 1, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Tight Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffDMFindCJdRPtCutCanvas, FinalEffDMFindCJdRPtCut_,
	 1, 1, 1, kBlack, 1, 20, "#DeltaR(CleanJets Matched + DecayModeFinding) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoRECOPtGenCanvas, MatchedLooseIsoRECOPtGen_,
	 1, 0, 0, kBlack, 7, 20, "PtGen(RECO Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoRECOPtGenCanvas, MatchedMedIsoRECOPtGen_,
	 1, 0, 0, kBlack, 7, 20, "PtGen(RECO Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoRECOPtGenCanvas, MatchedTightIsoRECOPtGen_,
	 1, 0, 0, kBlack, 7, 20, "PtGen(RECO Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedDMFindRECOPtGenCanvas, MatchedDMFindRECOPtGen_,
	 1, 0, 0, kBlack, 7, 20, "PtGen(RECO Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedRECOPtGenCanvas, MatchedRECOPtGen_,
	 1, 0, 0, kBlack, 7, 20, "PtGen(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoRECOPtGenCanvas, FinalEffLooseIsoRECOPtGen_,
	 1, 1, 1, kBlack, 1, 20, "PtGen(RECO Matched + Loose Iso + DM) / PtGen(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoRECOPtGenCanvas, FinalEffMedIsoRECOPtGen_,
	 1, 1, 1, kBlack, 1, 20, "PtGen(RECO Matched + Med Iso + DM) / PtGen(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoRECOPtGenCanvas, FinalEffTightIsoRECOPtGen_,
	 1, 1, 1, kBlack, 1, 20, "PtGen(RECO Matched + Tight Iso + DM) / PtGen(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffDMFindRECOPtGenCanvas, FinalEffDMFindRECOPtGen_,
	 1, 1, 1, kBlack, 1, 20, "PtGen(RECO Matched + DecayModeFinding) / PtGen(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoRECOPtGenOverDMCanvas, FinalEffLooseIsoRECOPtGenOverDM_,
         1, 1, 1, kBlack, 1, 20, "PtGen(RECO Matched + Loose Iso + DM) / PtGen(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoRECOPtGenOverDMCanvas, FinalEffMedIsoRECOPtGenOverDM_,
         1, 1, 1, kBlack, 1, 20, "PtGen(RECO Matched + Med Iso + DM) / PtGen(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoRECOPtGenOverDMCanvas, FinalEffTightIsoRECOPtGenOverDM_,
         1, 1, 1, kBlack, 1, 20, "PtGen(RECO Matched + Tight Iso + DM) / PtGen(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoCJPtGenCanvas, MatchedLooseIsoCJPtGen_,
         1, 0, 0, kBlack, 7, 20, "PtGen(CleanJets Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoCJPtGenCanvas, MatchedMedIsoCJPtGen_,
         1, 0, 0, kBlack, 7, 20, "PtGen(CleanJets Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoCJPtGenCanvas, MatchedTightIsoCJPtGen_,
         1, 0, 0, kBlack, 7, 20, "PtGen(CleanJets Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedDMFindCJPtGenCanvas, MatchedDMFindCJPtGen_,
         1, 0, 0, kBlack, 7, 20, "PtGen(CleanJets Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedCJPtGenCanvas, MatchedCJPtGen_,
         1, 0, 0, kBlack, 7, 20, "PtGen(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoCJPtGenCanvas, FinalEffLooseIsoCJPtGen_,
         1, 1, 1, kBlack, 1, 20, "PtGen(CleanJets Matched + Loose Iso + DM) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoCJPtGenCanvas, FinalEffMedIsoCJPtGen_,
         1, 1, 1, kBlack, 1, 20, "PtGen(CleanJets Matched + Med Iso + DM) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoCJPtGenCanvas, FinalEffTightIsoCJPtGen_,
         1, 1, 1, kBlack, 1, 20, "PtGen(CleanJets Matched + Tight Iso + DM) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffDMFindCJPtGenCanvas, FinalEffDMFindCJPtGen_,
         1, 1, 1, kBlack, 1, 20, "PtGen(CleanJets Matched + DecayModeFinding) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoCJPtGenOverDMCanvas, FinalEffLooseIsoCJPtGenOverDM_,
         1, 1, 1, kBlack, 1, 20, "PtGen(CleanJets Matched + Loose Iso + DM) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoCJPtGenOverDMCanvas, FinalEffMedIsoCJPtGenOverDM_,
         1, 1, 1, kBlack, 1, 20, "PtGen(CleanJets Matched + Med Iso + DM) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoCJPtGenOverDMCanvas, FinalEffTightIsoCJPtGenOverDM_,
         1, 1, 1, kBlack, 1, 20, "PtGen(CleanJets Matched + Tight Iso + DM) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngDMCJPtCanvas, OneProngDMCJPt_, 
         1, 0, 0, kBlack, 7, 20, "1 Prong p_{T}(CleanJets Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngOnePizDMCJPtCanvas, OneProngOnePizDMCJPt_, 
         1, 0, 0, kBlack, 7, 20, "1 Prong + 1 #pi^{0} p_{T}(CleanJets Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngTwoPizDMCJPtCanvas, OneProngTwoPizDMCJPt_, 
         1, 0, 0, kBlack, 7, 20, "1 Prong + 2 #pi^{0} p_{T}(CleanJets Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(ThreeProngDMCJPtCanvas, ThreeProngDMCJPt_, 
         1, 0, 0, kBlack, 7, 20, "3 Prong p_{T}(CleanJets Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngDMRECOPtCanvas, OneProngDMRECOPt_, 
         1, 0, 0, kBlack, 7, 20, "1 Prong p_{T}(RECO Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngOnePizDMRECOPtCanvas, OneProngOnePizDMRECOPt_, 
         1, 0, 0, kBlack, 7, 20, "1 Prong + 1 #pi^{0} p_{T}(RECO Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngTwoPizDMRECOPtCanvas, OneProngTwoPizDMRECOPt_, 
         1, 0, 0, kBlack, 7, 20, "1 Prong + 2 #pi^{0} p_{T}(RECO Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(ThreeProngDMRECOPtCanvas, ThreeProngDMRECOPt_, 
         1, 0, 0, kBlack, 7, 20, "3 Prong p_{T}(RECO Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedOneProngCJPtCanvas, MatchedOneProngCJPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong p_{T}(CleanJets Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedOneProngOnePizCJPtCanvas, MatchedOneProngOnePizCJPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 1 #pi^{0} p_{T}(CleanJets Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedOneProngTwoPizCJPtCanvas, MatchedOneProngTwoPizCJPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 2 #pi^{0} p_{T}(CleanJets Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedThreeProngCJPtCanvas, MatchedThreeProngCJPt_,
         1, 0, 0, kBlack, 7, 20, "3 Prong p_{T}(CleanJets Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedOneProngRECOPtCanvas, MatchedOneProngRECOPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong p_{T}(RECO Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedOneProngOnePizRECOPtCanvas, MatchedOneProngOnePizRECOPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 1 #pi^{0} p_{T}(RECO Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedOneProngTwoPizRECOPtCanvas, MatchedOneProngTwoPizRECOPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 2 #pi^{0} p_{T}(RECO Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedThreeProngRECOPtCanvas, MatchedThreeProngRECOPt_,
         1, 0, 0, kBlack, 7, 20, "3 Prong p_{T}(RECO Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngDMCJPtCanvas, FinalOneProngDMCJPt_, 
         1, 1, 1, kBlack, 1, 20, "1 Prong p_{T}(CleanJets Matched + DM) / 1 Prong p_{T}(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngOnePizDMCJPtCanvas, FinalOneProngOnePizDMCJPt_, 
         1, 1, 1, kBlack, 1, 20, "1 Prong + 1 #pi^{0} p_{T}(CleanJets Matched + DM) / 1 Prong + 1 #pi^{0} p_{T}(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngTwoPizDMCJPtCanvas, FinalOneProngTwoPizDMCJPt_, 
         1, 1, 1, kBlack, 1, 20, "1 Prong + 2 #pi^{0} p_{T}(CleanJets Matched + DM) / 1 Prong + 2 #pi^{0} p_{T}(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalThreeProngDMCJPtCanvas, FinalThreeProngDMCJPt_, 
         1, 1, 1, kBlack, 1, 20, "3 Prong p_{T}(CleanJets Matched + DM) / 3 Prong p_{T}(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngDMRECOPtCanvas, FinalOneProngDMRECOPt_, 
         1, 1, 1, kBlack, 1, 20, "1 Prong p_{T}(RECO Matched + DM) / 1 Prong p_{T}(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngOnePizDMRECOPtCanvas, FinalOneProngOnePizDMRECOPt_, 
         1, 1, 1, kBlack, 1, 20, "1 Prong + 1 #pi^{0} p_{T}(RECO Matched + DM) / 1 Prong + 1 #pi^{0} p_{T}(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngTwoPizDMRECOPtCanvas, FinalOneProngTwoPizDMRECOPt_, 
         1, 1, 1, kBlack, 1, 20, "1 Prong + 2 #pi^{0} p_{T}(RECO Matched + DM) / 1 Prong + 2 #pi^{0} p_{T}(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalThreeProngDMRECOPtCanvas, FinalThreeProngDMRECOPt_, 
         1, 1, 1, kBlack, 1, 20, "3 Prong p_{T}(RECO Matched + DM) / 3 Prong p_{T}(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

std::cout << "after formatting" << std::endl;
  
  FinalEffLooseIsoCJPt_->SetLineColor(kRed);
  FinalEffMedIsoCJPt_->SetLineColor(kRed);
  FinalEffTightIsoCJPt_->SetLineColor(kRed);
  FinalEffDMFindCJPt_->SetLineColor(kRed);
  FinalEffLooseIsoCJPtOverDM_->SetLineColor(kRed);
  FinalEffMedIsoCJPtOverDM_->SetLineColor(kRed);
  FinalEffTightIsoCJPtOverDM_->SetLineColor(kRed);

  FinalEffLooseIsoCJdR_->SetLineColor(kRed);
  FinalEffMedIsoCJdR_->SetLineColor(kRed);
  FinalEffTightIsoCJdR_->SetLineColor(kRed);
  FinalEffDMFindCJdR_->SetLineColor(kRed);
  FinalEffLooseIsoCJdROverDM_->SetLineColor(kRed);
  FinalEffMedIsoCJdROverDM_->SetLineColor(kRed);
  FinalEffTightIsoCJdROverDM_->SetLineColor(kRed);

  FinalEffLooseIsoCJdRPtCut_->SetLineColor(kRed);
  FinalEffMedIsoCJdRPtCut_->SetLineColor(kRed);
  FinalEffTightIsoCJdRPtCut_->SetLineColor(kRed);
  FinalEffDMFindCJdRPtCut_->SetLineColor(kRed);

  FinalEffLooseIsoCJPtGen_->SetLineColor(kRed);
  FinalEffMedIsoCJPtGen_->SetLineColor(kRed);
  FinalEffTightIsoCJPtGen_->SetLineColor(kRed);
  FinalEffDMFindCJPtGen_->SetLineColor(kRed);
  FinalEffLooseIsoCJPtGenOverDM_->SetLineColor(kRed);
  FinalEffMedIsoCJPtGenOverDM_->SetLineColor(kRed);
  FinalEffTightIsoCJPtGenOverDM_->SetLineColor(kRed);

  FinalOneProngDMCJPt_->SetLineColor(kRed);
  FinalOneProngOnePizDMCJPt_->SetLineColor(kRed);
  FinalOneProngTwoPizDMCJPt_->SetLineColor(kRed);
  FinalThreeProngDMCJPt_->SetLineColor(kRed);
  
  TLegend *leg = new TLegend(0.1,0.7,0.25,0.9);
  leg->AddEntry(FinalEffDMFindRECOPtGen_, "No CleanJets","L");
  leg->AddEntry(FinalEffDMFindCJPtGen_, "CleanJets","L");

  FinalEffLooseIsoSAMEPtCanvas.cd();
  FinalEffLooseIsoSAMEPt_->Add(FinalEffLooseIsoCJPt_);
  FinalEffLooseIsoSAMEPt_->Add(FinalEffLooseIsoRECOPt_);
  FinalEffLooseIsoSAMEPt_->Draw("ap");
  FinalEffLooseIsoSAMEPt_->GetXaxis()->SetTitle("p_{T}");
  FinalEffLooseIsoSAMEPt_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Loose Isolation / Matched)");
  gPad->Modified();
  leg->Draw();

  FinalEffMedIsoSAMEPtCanvas.cd();
  FinalEffMedIsoSAMEPt_->Add(FinalEffMedIsoCJPt_);
  FinalEffMedIsoSAMEPt_->Add(FinalEffMedIsoRECOPt_);
  FinalEffMedIsoSAMEPt_->Draw("ap");
  FinalEffMedIsoSAMEPt_->GetXaxis()->SetTitle("p_{T}");
  FinalEffMedIsoSAMEPt_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Med Isolation / Matched)");
  gPad->Modified();
  leg->Draw();

  FinalEffTightIsoSAMEPtCanvas.cd();
  FinalEffTightIsoSAMEPt_->Add(FinalEffTightIsoCJPt_);
  FinalEffTightIsoSAMEPt_->Add(FinalEffTightIsoRECOPt_);
  FinalEffTightIsoSAMEPt_->Draw("ap");
  FinalEffTightIsoSAMEPt_->GetXaxis()->SetTitle("p_{T}");
  FinalEffTightIsoSAMEPt_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Tight Isolation / Matched)");
  gPad->Modified();
  leg->Draw();

  FinalEffDMFindSAMEPtCanvas.cd();
  FinalEffDMFindSAMEPt_->Add(FinalEffDMFindCJPt_);
  FinalEffDMFindSAMEPt_->Add(FinalEffDMFindRECOPt_);
  FinalEffDMFindSAMEPt_->Draw("ap");
  FinalEffDMFindSAMEPt_->GetXaxis()->SetTitle("p_{T}");
  FinalEffDMFindSAMEPt_->GetYaxis()->SetTitle("#epsilon (Matched + DM / Matched)");
  gPad->Modified();
  leg->Draw();

  FinalEffLooseIsoSAMEPtOverDMCanvas.cd();
  FinalEffLooseIsoSAMEPtOverDM_->Add(FinalEffLooseIsoCJPtOverDM_);
  FinalEffLooseIsoSAMEPtOverDM_->Add(FinalEffLooseIsoRECOPtOverDM_);
  FinalEffLooseIsoSAMEPtOverDM_->Draw("ap");
  FinalEffLooseIsoSAMEPtOverDM_->GetXaxis()->SetTitle("p_{T}");
  FinalEffLooseIsoSAMEPtOverDM_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Loose Isolation / Matched + DM)");
  gPad->Modified();
  leg->Draw();

  FinalEffMedIsoSAMEPtOverDMCanvas.cd();
  FinalEffMedIsoSAMEPtOverDM_->Add(FinalEffMedIsoCJPtOverDM_);
  FinalEffMedIsoSAMEPtOverDM_->Add(FinalEffMedIsoRECOPtOverDM_);
  FinalEffMedIsoSAMEPtOverDM_->Draw("ap");
  FinalEffMedIsoSAMEPtOverDM_->GetXaxis()->SetTitle("p_{T}");
  FinalEffMedIsoSAMEPtOverDM_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Med Isolation / Matched + DM)");
  gPad->Modified();
  leg->Draw();

  FinalEffTightIsoSAMEPtOverDMCanvas.cd();
  FinalEffTightIsoSAMEPtOverDM_->Add(FinalEffTightIsoCJPtOverDM_);
  FinalEffTightIsoSAMEPtOverDM_->Add(FinalEffTightIsoRECOPtOverDM_);
  FinalEffTightIsoSAMEPtOverDM_->Draw("ap");
  FinalEffTightIsoSAMEPtOverDM_->GetXaxis()->SetTitle("p_{T}");
  FinalEffTightIsoSAMEPtOverDM_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Tight Isolation / Matched + DM)");
  gPad->Modified();
  leg->Draw();


////////////////////////
// For dR
////////////////////////
  FinalEffLooseIsoSAMEdRCanvas.cd();
  FinalEffLooseIsoSAMEdR_->Add(FinalEffLooseIsoCJdR_);
  FinalEffLooseIsoSAMEdR_->Add(FinalEffLooseIsoRECOdR_);
  FinalEffLooseIsoSAMEdR_->Draw("ap");
  FinalEffLooseIsoSAMEdR_->GetXaxis()->SetTitle("#DeltaR (#tau#tau)");
  FinalEffLooseIsoSAMEdR_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Loose Isolation / Matched)");
  gPad->Modified();
  leg->Draw();

  FinalEffMedIsoSAMEdRCanvas.cd();
  FinalEffMedIsoSAMEdR_->Add(FinalEffMedIsoCJdR_);
  FinalEffMedIsoSAMEdR_->Add(FinalEffMedIsoRECOdR_);
  FinalEffMedIsoSAMEdR_->Draw("ap");
  FinalEffMedIsoSAMEdR_->GetXaxis()->SetTitle("#DeltaR (#tau#tau)");
  FinalEffMedIsoSAMEdR_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Med Isolation / Matched)");
  gPad->Modified();
  leg->Draw();

  FinalEffTightIsoSAMEdRCanvas.cd();
  FinalEffTightIsoSAMEdR_->Add(FinalEffTightIsoCJdR_);
  FinalEffTightIsoSAMEdR_->Add(FinalEffTightIsoRECOdR_);
  FinalEffTightIsoSAMEdR_->Draw("ap");
  FinalEffTightIsoSAMEdR_->GetXaxis()->SetTitle("#DeltaR (#tau#tau)");
  FinalEffTightIsoSAMEdR_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Tight Isolation / Matched)");
  gPad->Modified();
  leg->Draw();

  FinalEffDMFindSAMEdRCanvas.cd();
  FinalEffDMFindSAMEdR_->Add(FinalEffDMFindCJdR_);
  FinalEffDMFindSAMEdR_->Add(FinalEffDMFindRECOdR_);
  FinalEffDMFindSAMEdR_->Draw("ap");
  FinalEffDMFindSAMEdR_->GetXaxis()->SetTitle("#DeltaR (#tau#tau)");
  FinalEffDMFindSAMEdR_->GetYaxis()->SetTitle("#epsilon (Matched + DM / Matched)");
  gPad->Modified();
  leg->Draw();

  FinalEffLooseIsoSAMEdROverDMCanvas.cd();
  FinalEffLooseIsoSAMEdROverDM_->Add(FinalEffLooseIsoCJdROverDM_);
  FinalEffLooseIsoSAMEdROverDM_->Add(FinalEffLooseIsoRECOdROverDM_);
  FinalEffLooseIsoSAMEdROverDM_->Draw("ap");
  FinalEffLooseIsoSAMEdROverDM_->GetXaxis()->SetTitle("#DeltaR (#tau#tau)");
  FinalEffLooseIsoSAMEdROverDM_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Loose Isolation / Matched + DM)");
  gPad->Modified();
  leg->Draw();

  FinalEffMedIsoSAMEdROverDMCanvas.cd();
  FinalEffMedIsoSAMEdROverDM_->Add(FinalEffMedIsoCJdROverDM_);
  FinalEffMedIsoSAMEdROverDM_->Add(FinalEffMedIsoRECOdROverDM_);
  FinalEffMedIsoSAMEdROverDM_->Draw("ap");
  FinalEffMedIsoSAMEdROverDM_->GetXaxis()->SetTitle("#DeltaR (#tau#tau)");
  FinalEffMedIsoSAMEdROverDM_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Med Isolation / Matched + DM)");
  gPad->Modified();
  leg->Draw();

  FinalEffTightIsoSAMEdROverDMCanvas.cd();
  FinalEffTightIsoSAMEdROverDM_->Add(FinalEffTightIsoCJdROverDM_);
  FinalEffTightIsoSAMEdROverDM_->Add(FinalEffTightIsoRECOdROverDM_);
  FinalEffTightIsoSAMEdROverDM_->Draw("ap");
  FinalEffTightIsoSAMEdROverDM_->GetXaxis()->SetTitle("#DeltaR (#tau#tau)");
  FinalEffTightIsoSAMEdROverDM_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Tight Isolation / Matched + DM)");
  gPad->Modified();
  leg->Draw();

  FinalEffLooseIsoSAMEdRPtCutCanvas.cd();
  FinalEffLooseIsoSAMEdRPtCut_->Add(FinalEffLooseIsoCJdRPtCut_);
  FinalEffLooseIsoSAMEdRPtCut_->Add(FinalEffLooseIsoRECOdRPtCut_);
  FinalEffLooseIsoSAMEdRPtCut_->Draw("ap");
  FinalEffLooseIsoSAMEdRPtCut_->GetXaxis()->SetTitle("#DeltaR (#tau#tau) with RECO p_{T} > 20");
  FinalEffLooseIsoSAMEdRPtCut_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Loose Isolation / Matched)");
  gPad->Modified();
  leg->Draw();

  FinalEffMedIsoSAMEdRPtCutCanvas.cd();
  FinalEffMedIsoSAMEdRPtCut_->Add(FinalEffMedIsoCJdRPtCut_);
  FinalEffMedIsoSAMEdRPtCut_->Add(FinalEffMedIsoRECOdRPtCut_);
  FinalEffMedIsoSAMEdRPtCut_->Draw("ap");
  FinalEffMedIsoSAMEdRPtCut_->GetXaxis()->SetTitle("#DeltaR (#tau#tau) with RECO p_{T} > 20");
  FinalEffMedIsoSAMEdRPtCut_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Med Isolation / Matched)");
  gPad->Modified();
  leg->Draw();

  FinalEffTightIsoSAMEdRPtCutCanvas.cd();
  FinalEffTightIsoSAMEdRPtCut_->Add(FinalEffTightIsoCJdRPtCut_);
  FinalEffTightIsoSAMEdRPtCut_->Add(FinalEffTightIsoRECOdRPtCut_);
  FinalEffTightIsoSAMEdRPtCut_->Draw("ap");
  FinalEffTightIsoSAMEdRPtCut_->GetXaxis()->SetTitle("#DeltaR (#tau#tau) with RECO p_{T} > 20");
  FinalEffTightIsoSAMEdRPtCut_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Tight Isolation / Matched)");
  gPad->Modified();
  leg->Draw();

  FinalEffDMFindSAMEdRPtCutCanvas.cd();
  FinalEffDMFindSAMEdRPtCut_->Add(FinalEffDMFindCJdRPtCut_);
  FinalEffDMFindSAMEdRPtCut_->Add(FinalEffDMFindRECOdRPtCut_);
  FinalEffDMFindSAMEdRPtCut_->Draw("ap");
  FinalEffDMFindSAMEdRPtCut_->GetXaxis()->SetTitle("#DeltaR (#tau#tau) with RECO p_{T} > 20");
  FinalEffDMFindSAMEdRPtCut_->GetYaxis()->SetTitle("#epsilon (Matched + DM / Matched)");
  gPad->Modified();
  leg->Draw();


////////////////////////
// For Pt Gen
////////////////////////
  FinalEffLooseIsoSAMEPtGenCanvas.cd();
  FinalEffLooseIsoSAMEPtGen_->Add(FinalEffLooseIsoCJPtGen_);
  FinalEffLooseIsoSAMEPtGen_->Add(FinalEffLooseIsoRECOPtGen_);
  FinalEffLooseIsoSAMEPtGen_->Draw("ap");
  FinalEffLooseIsoSAMEPtGen_->GetXaxis()->SetTitle("p_{T} Gen-Visible");
  FinalEffLooseIsoSAMEPtGen_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Loose Isolation / Matched)");
  gPad->Modified();
  leg->Draw();

  FinalEffMedIsoSAMEPtGenCanvas.cd();
  FinalEffMedIsoSAMEPtGen_->Add(FinalEffMedIsoCJPtGen_);
  FinalEffMedIsoSAMEPtGen_->Add(FinalEffMedIsoRECOPtGen_);
  FinalEffMedIsoSAMEPtGen_->Draw("ap");
  FinalEffMedIsoSAMEPtGen_->GetXaxis()->SetTitle("p_{T} Gen-Visible");
  FinalEffMedIsoSAMEPtGen_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Med Isolation / Matched)");
  gPad->Modified();
  leg->Draw();

  FinalEffTightIsoSAMEPtGenCanvas.cd();
  FinalEffTightIsoSAMEPtGen_->Add(FinalEffTightIsoCJPtGen_);
  FinalEffTightIsoSAMEPtGen_->Add(FinalEffTightIsoRECOPtGen_);
  FinalEffTightIsoSAMEPtGen_->Draw("ap");
  FinalEffTightIsoSAMEPtGen_->GetXaxis()->SetTitle("p_{T} Gen-Visible");
  FinalEffTightIsoSAMEPtGen_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Tight Isolation / Matched)");
  gPad->Modified();
  leg->Draw();

  FinalEffDMFindSAMEPtGenCanvas.cd();
  FinalEffDMFindSAMEPtGen_->Add(FinalEffDMFindCJPtGen_);
  FinalEffDMFindSAMEPtGen_->Add(FinalEffDMFindRECOPtGen_);
  FinalEffDMFindSAMEPtGen_->Draw("ap");
  FinalEffDMFindSAMEPtGen_->GetXaxis()->SetTitle("p_{T} Gen-Visible");
  FinalEffDMFindSAMEPtGen_->GetYaxis()->SetTitle("#epsilon (Matched + DM / Matched)");
  gPad->Modified();
  leg->Draw();

  FinalEffLooseIsoSAMEPtGenOverDMCanvas.cd();
  FinalEffLooseIsoSAMEPtGenOverDM_->Add(FinalEffLooseIsoCJPtGenOverDM_);
  FinalEffLooseIsoSAMEPtGenOverDM_->Add(FinalEffLooseIsoRECOPtGenOverDM_);
  FinalEffLooseIsoSAMEPtGenOverDM_->Draw("ap");
  FinalEffLooseIsoSAMEPtGenOverDM_->GetXaxis()->SetTitle("p_{T} Gen-Visible");
  FinalEffLooseIsoSAMEPtGenOverDM_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Loose Isolation / Matched + DM)");
  gPad->Modified();
  leg->Draw();
  
  FinalEffMedIsoSAMEPtGenOverDMCanvas.cd();
  FinalEffMedIsoSAMEPtGenOverDM_->Add(FinalEffMedIsoCJPtGenOverDM_);
  FinalEffMedIsoSAMEPtGenOverDM_->Add(FinalEffMedIsoRECOPtGenOverDM_);
  FinalEffMedIsoSAMEPtGenOverDM_->Draw("ap");
  FinalEffMedIsoSAMEPtGenOverDM_->GetXaxis()->SetTitle("p_{T} Gen-Visible");
  FinalEffMedIsoSAMEPtGenOverDM_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Med Isolation / Matched + DM)");
  gPad->Modified();
  leg->Draw();
  
  FinalEffTightIsoSAMEPtGenOverDMCanvas.cd();
  FinalEffTightIsoSAMEPtGenOverDM_->Add(FinalEffTightIsoCJPtGenOverDM_);
  FinalEffTightIsoSAMEPtGenOverDM_->Add(FinalEffTightIsoRECOPtGenOverDM_);
  FinalEffTightIsoSAMEPtGenOverDM_->Draw("ap");
  FinalEffTightIsoSAMEPtGenOverDM_->GetXaxis()->SetTitle("p_{T} Gen-Visible");
  FinalEffTightIsoSAMEPtGenOverDM_->GetYaxis()->SetTitle("#epsilon (Matched + DM + Tight Isolation / Matched + DM )");
  gPad->Modified();
  leg->Draw();

////////////////////////
// For By Decay Mode
////////////////////////

  FinalEffOneProngDMSAMEPtCanvas.cd();
  FinalEffOneProngDMSAMEPt_->Add(FinalOneProngDMCJPt_);
  FinalEffOneProngDMSAMEPt_->Add(FinalOneProngDMRECOPt_);
  FinalEffOneProngDMSAMEPt_->Draw("ap");
  FinalEffOneProngDMSAMEPt_->GetXaxis()->SetTitle("1 Prong Decay: p_{T}");
  FinalEffOneProngDMSAMEPt_->GetYaxis()->SetTitle("#epsilon (Matched + DM / Matched)");
  gPad->Modified();
  leg->Draw();

  FinalEffOneProngOnePizDMSAMEPtCanvas.cd();
  FinalEffOneProngOnePizDMSAMEPt_->Add(FinalOneProngOnePizDMCJPt_);
  FinalEffOneProngOnePizDMSAMEPt_->Add(FinalOneProngOnePizDMRECOPt_);
  FinalEffOneProngOnePizDMSAMEPt_->Draw("ap");
  FinalEffOneProngOnePizDMSAMEPt_->GetXaxis()->SetTitle("1 Prong + 1 #pi^{0} Decay: p_{T}");
  FinalEffOneProngOnePizDMSAMEPt_->GetYaxis()->SetTitle("#epsilon (Matched + DM / Matched)");
  gPad->Modified();
  leg->Draw();

  FinalEffOneProngTwoPizDMSAMEPtCanvas.cd();
  FinalEffOneProngTwoPizDMSAMEPt_->Add(FinalOneProngTwoPizDMCJPt_);
  FinalEffOneProngTwoPizDMSAMEPt_->Add(FinalOneProngTwoPizDMRECOPt_);
  FinalEffOneProngTwoPizDMSAMEPt_->Draw("ap");
  FinalEffOneProngTwoPizDMSAMEPt_->GetXaxis()->SetTitle("1 Prong + 2 #pi^{0} Decay: p_{T}");
  FinalEffOneProngTwoPizDMSAMEPt_->GetYaxis()->SetTitle("#epsilon (Matched + DM / Matched)");
  gPad->Modified();
  leg->Draw();
  
  FinalEffThreeProngDMSAMEPtCanvas.cd();
  FinalEffThreeProngDMSAMEPt_->Add(FinalThreeProngDMCJPt_);
  FinalEffThreeProngDMSAMEPt_->Add(FinalThreeProngDMRECOPt_);
  FinalEffThreeProngDMSAMEPt_->Draw("ap");
  FinalEffThreeProngDMSAMEPt_->GetXaxis()->SetTitle("3 Prong Decay: p_{T}");
  FinalEffThreeProngDMSAMEPt_->GetYaxis()->SetTitle("#epsilon (Matched + DM / Matched)");
  gPad->Modified();
  leg->Draw();


std::cout << "check before setting axis" << std::endl;

/*  FinalEffLooseIsoSAMEPt_->GetXaxis()->SetTitle("p_{T} (Matched + DM + Loose Isolation) / p_{T} (Matched)");  
  FinalEffMedIsoSAMEPt_->GetXaxis()->SetTitle("p_{T} (Matched + DM + Medium Isolation) / p_{T} (Matched)");
  FinalEffTightIsoSAMEPt_->GetXaxis()->SetTitle("p_{T} (Matched + DM + Tight Isolation) / p_{T} (Matched)");  
  FinalEffDMFindSAMEPt_->GetXaxis()->SetTitle("p_{T} (Matched + DM) / p_{T} (Matched)");
  FinalEffLooseIsoSAMEdR_->GetXaxis()->SetTitle("#DeltaR (Matched + DM + Loose Isolation) / #DeltaR (Matched)");  
  FinalEffMedIsoSAMEdR_->GetXaxis()->SetTitle("#DeltaR (Matched + DM + Medium Isolation) / #DeltaR (Matched)");
  FinalEffTightIsoSAMEdR_->GetXaxis()->SetTitle("#DeltaR (Matched + DM + Tight Isolation) / #DeltaR (Matched)");  
  FinalEffDMFindSAMEdR_->GetXaxis()->SetTitle("#DeltaR (Matched + DM) / #DeltaR (Matched)");
  FinalEffLooseIsoSAMEPtGen_->GetXaxis()->SetTitle("p_{T} Gen-Visible (Matched + DM + Loose Isolation) / p_{T} Gen-Visible (Matched)");
  FinalEffMedIsoSAMEPtGen_->GetXaxis()->SetTitle("p_{T} Gen-Visible (Matched + DM + Medium Isolation) / p_{T} Gen-Visible (Matched)"); 
  FinalEffTightIsoSAMEPtGen_->GetXaxis()->SetTitle("p_{T} Gen-Visible (Matched + DM + Tight Isolation) / p_{T} Gen-Visible (Matched)");
  FinalEffDMFindSAMEPtGen_->GetXaxis()->SetTitle("p_{T} Gen-Visible (Matched + DM) / p_{T} Gen-Visible (Matched)"); 
*/
std::cout << "<----------------Formatted Canvases and Histos-------------->" << std::endl;

  //Write output file
  out_->cd();

  NEventsCanvas.Write();
  GenMatchTypeCanvas.Write();
  GenMatchPDGIDCanvas.Write();
  GenMatchPDGIDwithLargeEtaCanvas.Write();
  NMuRemovedCanvas.Write();
  TauMuTauHaddRCanvas.Write();
  MassDiLepGENCanvas.Write();
  MassDiLepRECOCanvas.Write();
  Mu3LargestPtMatchCanvas.Write();
  Mu3PtMatchCanvas.Write();
  MatchedTauMuPtLeptonCanvas.Write();
  MatchedTauHadPtLeptonCanvas.Write();
  MatchedTauHadEtaLeptonCanvas.Write();
  MatchedHTLeptonCanvas.Write();
  MatchedBDiscCSVLeptonCanvas.Write();
  MatchedTauMuPtJetCanvas.Write();
  MatchedTauHadPtJetCanvas.Write();
  MatchedTauHadEtaJetCanvas.Write();
  MatchedHTJetCanvas.Write();
  MatchedBDiscCSVJetCanvas.Write();
  NEventsCutsCanvas.Write();
  NConstituentsCJCanvas.Write();
  PtOverMassMu1Mu2Canvas.Write();
  PtMu1Mu2Canvas.Write();

  NTauDecayModeGENCanvas.Write();
  NTauDecayModeRECOCJCanvas.Write();
  NTauDecayModeRECONoCJCanvas.Write();
  NTauDecayModeRECOCJIfGen3Canvas.Write();
  NTauDecayModeRECONoCJIfGen3Canvas.Write();
  NTausRECOvsCLEANJETSCanvas.Write();
  GenDiTaudRvsCJDiTaudRCanvas.Write();
  TauHadnConstvsPtCanvas.Write();
  RECOPtvsGENVisiblePtCJCanvas.Write();
  RECOPtvsGENVisiblePtRECOCanvas.Write();
  RECOCJTauDMvsGENTauDMCanvas.Write();
  RECONoCJTauDMvsGENTauDMCanvas.Write();
  RECOCJTauDMvsGENTauDMMedIsoCanvas.Write();
  RECONoCJTauDMvsGENTauDMMedIsoCanvas.Write();
  GenVisiblePtvsRecoCJPtIfCorrDMCanvas.Write();
  GenVisiblePtvsRecoNoCJPtIfCorrDMCanvas.Write();
  RECOPtvsGENdRCJCanvas.Write();
  RECOPtvsGENdRRECOCanvas.Write();
  RECOPtvsGENdRCJMedIsoCanvas.Write();
  RECOPtvsGENdRRECOMedIsoCanvas.Write();
  MatchedTauHadEtavsTauHadPtCanvas.Write();
  BDiscCSVvsRawIsoCanvas.Write();

  MatchedLooseIsoRECOPtCanvas.Write();
  MatchedMedIsoRECOPtCanvas.Write();
  MatchedTightIsoRECOPtCanvas.Write();
  MatchedDMFindRECOPtCanvas.Write();
  MatchedRECOPtCanvas.Write();
  FinalEffLooseIsoRECOPtCanvas.Write();
  FinalEffMedIsoRECOPtCanvas.Write();
  FinalEffTightIsoRECOPtCanvas.Write();
  FinalEffDMFindRECOPtCanvas.Write();
  FinalEffLooseIsoRECOPtOverDMCanvas.Write();
  FinalEffMedIsoRECOPtOverDMCanvas.Write();
  FinalEffTightIsoRECOPtOverDMCanvas.Write();


  MatchedLooseIsoCJPtCanvas.Write();
  MatchedMedIsoCJPtCanvas.Write();
  MatchedTightIsoCJPtCanvas.Write();
  MatchedDMFindCJPtCanvas.Write();
  MatchedCJPtCanvas.Write();
  FinalEffLooseIsoCJPtCanvas.Write();
  FinalEffMedIsoCJPtCanvas.Write();
  FinalEffTightIsoCJPtCanvas.Write();
  FinalEffDMFindCJPtCanvas.Write();
  FinalEffLooseIsoCJPtOverDMCanvas.Write();
  FinalEffMedIsoCJPtOverDMCanvas.Write();
  FinalEffTightIsoCJPtOverDMCanvas.Write();

  MatchedLooseIsoRECOdRCanvas.Write();
  MatchedMedIsoRECOdRCanvas.Write();
  MatchedTightIsoRECOdRCanvas.Write();
  MatchedDMFindRECOdRCanvas.Write();
  MatchedRECOdRCanvas.Write();
  FinalEffLooseIsoRECOdRCanvas.Write();
  FinalEffMedIsoRECOdRCanvas.Write();
  FinalEffTightIsoRECOdRCanvas.Write();
  FinalEffDMFindRECOdRCanvas.Write();
  FinalEffLooseIsoRECOdROverDMCanvas.Write();
  FinalEffMedIsoRECOdROverDMCanvas.Write();
  FinalEffTightIsoRECOdROverDMCanvas.Write();

  MatchedLooseIsoCJdRCanvas.Write();
  MatchedMedIsoCJdRCanvas.Write();
  MatchedTightIsoCJdRCanvas.Write();
  MatchedDMFindCJdRCanvas.Write();
  MatchedCJdRCanvas.Write();
  FinalEffLooseIsoCJdRCanvas.Write();
  FinalEffMedIsoCJdRCanvas.Write();
  FinalEffTightIsoCJdRCanvas.Write();
  FinalEffDMFindCJdRCanvas.Write();
  FinalEffLooseIsoCJdROverDMCanvas.Write();
  FinalEffMedIsoCJdROverDMCanvas.Write();
  FinalEffTightIsoCJdROverDMCanvas.Write();
  
  MatchedLooseIsoRECOdRPtCutCanvas.Write();
  MatchedMedIsoRECOdRPtCutCanvas.Write();
  MatchedTightIsoRECOdRPtCutCanvas.Write();
  MatchedDMFindRECOdRPtCutCanvas.Write();
  MatchedRECOdRPtCutCanvas.Write();
  FinalEffLooseIsoRECOdRPtCutCanvas.Write();
  FinalEffMedIsoRECOdRPtCutCanvas.Write();
  FinalEffTightIsoRECOdRPtCutCanvas.Write();
  FinalEffDMFindRECOdRPtCutCanvas.Write();

  MatchedLooseIsoCJdRPtCutCanvas.Write();
  MatchedMedIsoCJdRPtCutCanvas.Write();
  MatchedTightIsoCJdRPtCutCanvas.Write();
  MatchedDMFindCJdRPtCutCanvas.Write();
  MatchedCJdRPtCutCanvas.Write();
  FinalEffLooseIsoCJdRPtCutCanvas.Write();
  FinalEffMedIsoCJdRPtCutCanvas.Write();
  FinalEffTightIsoCJdRPtCutCanvas.Write();
  FinalEffDMFindCJdRPtCutCanvas.Write();

  MatchedLooseIsoRECOPtGenCanvas.Write();
  MatchedMedIsoRECOPtGenCanvas.Write();
  MatchedTightIsoRECOPtGenCanvas.Write();
  MatchedDMFindRECOPtGenCanvas.Write();
  MatchedRECOPtGenCanvas.Write();
  FinalEffLooseIsoRECOPtGenCanvas.Write();
  FinalEffMedIsoRECOPtGenCanvas.Write();
  FinalEffTightIsoRECOPtGenCanvas.Write();
  FinalEffDMFindRECOPtGenCanvas.Write();
  FinalEffLooseIsoRECOPtGenOverDMCanvas.Write();
  FinalEffMedIsoRECOPtGenOverDMCanvas.Write();
  FinalEffTightIsoRECOPtGenOverDMCanvas.Write();

  MatchedLooseIsoCJPtGenCanvas.Write();
  MatchedMedIsoCJPtGenCanvas.Write();
  MatchedTightIsoCJPtGenCanvas.Write();
  MatchedDMFindCJPtGenCanvas.Write();
  MatchedCJPtGenCanvas.Write();
  FinalEffLooseIsoCJPtGenCanvas.Write();
  FinalEffMedIsoCJPtGenCanvas.Write();
  FinalEffTightIsoCJPtGenCanvas.Write();
  FinalEffDMFindCJPtGenCanvas.Write();
  FinalEffLooseIsoCJPtGenOverDMCanvas.Write();
  FinalEffMedIsoCJPtGenOverDMCanvas.Write();
  FinalEffTightIsoCJPtGenOverDMCanvas.Write();

  OneProngDMCJPtCanvas.Write();
  OneProngOnePizDMCJPtCanvas.Write();
  OneProngTwoPizDMCJPtCanvas.Write();
  ThreeProngDMCJPtCanvas.Write();
  OneProngDMRECOPtCanvas.Write();
  OneProngOnePizDMRECOPtCanvas.Write();
  OneProngTwoPizDMRECOPtCanvas.Write();
  ThreeProngDMRECOPtCanvas.Write();
  MatchedOneProngCJPtCanvas.Write();
  MatchedOneProngOnePizCJPtCanvas.Write();
  MatchedOneProngTwoPizCJPtCanvas.Write();
  MatchedThreeProngCJPtCanvas.Write();
  MatchedOneProngRECOPtCanvas.Write();
  MatchedOneProngOnePizRECOPtCanvas.Write();
  MatchedOneProngTwoPizRECOPtCanvas.Write();
  MatchedThreeProngRECOPtCanvas.Write();
  FinalOneProngDMCJPtCanvas.Write();
  FinalOneProngOnePizDMCJPtCanvas.Write();
  FinalOneProngTwoPizDMCJPtCanvas.Write();
  FinalThreeProngDMCJPtCanvas.Write();
  FinalOneProngDMRECOPtCanvas.Write();
  FinalOneProngOnePizDMRECOPtCanvas.Write();
  FinalOneProngTwoPizDMRECOPtCanvas.Write();
  FinalThreeProngDMRECOPtCanvas.Write();
  FinalEffOneProngDMSAMEPtCanvas.Write();
  FinalEffOneProngOnePizDMSAMEPtCanvas.Write();
  FinalEffOneProngTwoPizDMSAMEPtCanvas.Write();
  FinalEffThreeProngDMSAMEPtCanvas.Write();



  FinalEffLooseIsoSAMEPtCanvas.Write();
  FinalEffMedIsoSAMEPtCanvas.Write();
  FinalEffTightIsoSAMEPtCanvas.Write();
  FinalEffDMFindSAMEPtCanvas.Write();
  FinalEffLooseIsoSAMEdRCanvas.Write();
  FinalEffMedIsoSAMEdRCanvas.Write();
  FinalEffTightIsoSAMEdRCanvas.Write();
  FinalEffDMFindSAMEdRCanvas.Write();
  FinalEffLooseIsoSAMEdRPtCutCanvas.Write();
  FinalEffMedIsoSAMEdRPtCutCanvas.Write();
  FinalEffTightIsoSAMEdRPtCutCanvas.Write();
  FinalEffDMFindSAMEdRPtCutCanvas.Write();
  FinalEffLooseIsoSAMEPtGenCanvas.Write();
  FinalEffMedIsoSAMEPtGenCanvas.Write();
  FinalEffTightIsoSAMEPtGenCanvas.Write();
  FinalEffDMFindSAMEPtGenCanvas.Write();

  FinalEffLooseIsoSAMEPtOverDMCanvas.Write();
  FinalEffMedIsoSAMEPtOverDMCanvas.Write();
  FinalEffTightIsoSAMEPtOverDMCanvas.Write();
  FinalEffLooseIsoSAMEdROverDMCanvas.Write();
  FinalEffMedIsoSAMEdROverDMCanvas.Write();
  FinalEffTightIsoSAMEdROverDMCanvas.Write();
  FinalEffLooseIsoSAMEPtGenOverDMCanvas.Write();
  FinalEffMedIsoSAMEPtGenOverDMCanvas.Write();
  FinalEffTightIsoSAMEPtGenOverDMCanvas.Write();

  out_->Write();
  out_->Close();
std::cout << "DONE" << std::endl;
}//EndJob

// ------------ method called when starting to processes a run  ------------
void FakeRateAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void FakeRateAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void FakeRateAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void FakeRateAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

//Delete Memory
void FakeRateAnalyzer::reset(const bool doDelete)
{
  if ((doDelete) && (NEvents_ != NULL)) delete NEvents_;
  NEvents_ = NULL;
  if ((doDelete) && (GenMatchType_ != NULL)) delete GenMatchType_;
  GenMatchType_ = NULL;
  if ((doDelete) && (GenMatchPDGID_ != NULL)) delete GenMatchPDGID_;
  GenMatchPDGID_ = NULL;
  if ((doDelete) && (GenMatchPDGIDwithLargeEta_ != NULL)) delete GenMatchPDGIDwithLargeEta_;
  GenMatchPDGIDwithLargeEta_ = NULL;
  if ((doDelete) && (NMuRemoved_ != NULL)) delete NMuRemoved_;
  NMuRemoved_ = NULL;
  if ((doDelete) && (TauMuTauHaddR_ != NULL)) delete TauMuTauHaddR_;
  TauMuTauHaddR_ = NULL;
  if ((doDelete) && (MassDiLepGEN_ != NULL)) delete MassDiLepGEN_;
  MassDiLepGEN_ = NULL;
  if ((doDelete) && (MassDiLepRECO_ != NULL)) delete MassDiLepRECO_;
  MassDiLepRECO_ = NULL;
  if ((doDelete) && (Mu3LargestPtMatch_ != NULL)) delete Mu3LargestPtMatch_;
  Mu3LargestPtMatch_ = NULL;
  if ((doDelete) && (Mu3PtMatch_ != NULL)) delete Mu3PtMatch_;
  Mu3PtMatch_ = NULL;
  if ((doDelete) && (MatchedTauMuPtLepton_ != NULL)) delete MatchedTauMuPtLepton_;
  MatchedTauMuPtLepton_ = NULL;
  if ((doDelete) && (MatchedTauHadPtLepton_ != NULL)) delete MatchedTauHadPtLepton_;
  MatchedTauHadPtLepton_ = NULL;
  if ((doDelete) && (MatchedTauHadEtaLepton_ != NULL)) delete MatchedTauHadEtaLepton_;
  MatchedTauHadEtaLepton_ = NULL;
  if ((doDelete) && (MatchedHTLepton_ != NULL)) delete MatchedHTLepton_;
  MatchedHTLepton_ = NULL;
  if ((doDelete) && (MatchedBDiscCSVLepton_ != NULL)) delete MatchedBDiscCSVLepton_;
  MatchedBDiscCSVLepton_ = NULL;
  if ((doDelete) && (MatchedTauMuPtJet_ != NULL)) delete MatchedTauMuPtJet_;
  MatchedTauMuPtJet_ = NULL;
  if ((doDelete) && (MatchedTauHadPtJet_ != NULL)) delete MatchedTauHadPtJet_;
  MatchedTauHadPtJet_ = NULL;
  if ((doDelete) && (MatchedTauHadEtaJet_ != NULL)) delete MatchedTauHadEtaJet_;
  MatchedTauHadEtaJet_ = NULL;
  if ((doDelete) && (MatchedHTJet_ != NULL)) delete MatchedHTJet_;
  MatchedHTJet_ = NULL;
  if ((doDelete) && (MatchedBDiscCSVJet_ != NULL)) delete MatchedBDiscCSVJet_;
  MatchedBDiscCSVJet_ = NULL;
  if ((doDelete) && (NEventsCuts_ != NULL)) delete NEventsCuts_;
  NEventsCuts_ = NULL;
  if ((doDelete) && (NConstituentsCJ_ != NULL)) delete NConstituentsCJ_;
  NConstituentsCJ_ = NULL;
  if ((doDelete) && (PtOverMassMu1Mu2_ != NULL)) delete PtOverMassMu1Mu2_;
  PtOverMassMu1Mu2_ = NULL;
  if ((doDelete) && (PtMu1Mu2_ != NULL)) delete PtMu1Mu2_;
  PtMu1Mu2_ = NULL;

  if ((doDelete) && (NTauDecayModeGEN_ != NULL)) delete NTauDecayModeGEN_;
  NTauDecayModeGEN_ = NULL;
  if ((doDelete) && (NTauDecayModeRECOCJ_ != NULL)) delete NTauDecayModeRECOCJ_;
  NTauDecayModeRECOCJ_ = NULL;
  if ((doDelete) && (NTauDecayModeRECONoCJ_ != NULL)) delete NTauDecayModeRECONoCJ_;
  NTauDecayModeRECONoCJ_ = NULL;
  if ((doDelete) && (NTauDecayModeRECOCJIfGen3_ != NULL)) delete NTauDecayModeRECOCJIfGen3_;
  NTauDecayModeRECOCJIfGen3_ = NULL;
  if ((doDelete) && (NTauDecayModeRECONoCJIfGen3_ != NULL)) delete NTauDecayModeRECONoCJIfGen3_;
  NTauDecayModeRECONoCJIfGen3_ = NULL;
  if ((doDelete) && (GenDiTaudRvsCJDiTaudR_ != NULL)) delete GenDiTaudRvsCJDiTaudR_;
  GenDiTaudRvsCJDiTaudR_ = NULL;
  if ((doDelete) && (TauHadnConstvsPt_ != NULL)) delete TauHadnConstvsPt_;
  TauHadnConstvsPt_ = NULL;
  if ((doDelete) && (RECOPtvsGENVisiblePtCJ_ != NULL)) delete RECOPtvsGENVisiblePtCJ_;
  RECOPtvsGENVisiblePtCJ_ = NULL;
  if ((doDelete) && (RECOPtvsGENVisiblePtRECO_ != NULL)) delete RECOPtvsGENVisiblePtRECO_;
  RECOPtvsGENVisiblePtRECO_ = NULL;
  if ((doDelete) && (RECOCJTauDMvsGENTauDM_ != NULL)) delete RECOCJTauDMvsGENTauDM_;
  RECOCJTauDMvsGENTauDM_ = NULL;
  if ((doDelete) && (RECONoCJTauDMvsGENTauDM_ != NULL)) delete RECONoCJTauDMvsGENTauDM_;
  RECONoCJTauDMvsGENTauDM_ = NULL;
  if ((doDelete) && (RECOCJTauDMvsGENTauDMMedIso_ != NULL)) delete RECOCJTauDMvsGENTauDMMedIso_;
  RECOCJTauDMvsGENTauDMMedIso_ = NULL;
  if ((doDelete) && (RECONoCJTauDMvsGENTauDMMedIso_ != NULL)) delete RECONoCJTauDMvsGENTauDMMedIso_;
  RECONoCJTauDMvsGENTauDMMedIso_ = NULL;
  if ((doDelete) && (GenVisiblePtvsRecoCJPtIfCorrDM_ != NULL)) delete GenVisiblePtvsRecoCJPtIfCorrDM_;
  GenVisiblePtvsRecoCJPtIfCorrDM_ = NULL;
  if ((doDelete) && (GenVisiblePtvsRecoNoCJPtIfCorrDM_ != NULL)) delete GenVisiblePtvsRecoNoCJPtIfCorrDM_;
  GenVisiblePtvsRecoNoCJPtIfCorrDM_ = NULL;
  if ((doDelete) && (RECOPtvsGENdRCJ_ != NULL)) delete RECOPtvsGENdRCJ_;
  RECOPtvsGENdRCJ_ = NULL;
  if ((doDelete) && (RECOPtvsGENdRRECO_ != NULL)) delete RECOPtvsGENdRRECO_;
  RECOPtvsGENdRRECO_ = NULL;
  if ((doDelete) && (RECOPtvsGENdRCJMedIso_ != NULL)) delete RECOPtvsGENdRCJMedIso_;
  RECOPtvsGENdRCJMedIso_ = NULL;
  if ((doDelete) && (RECOPtvsGENdRRECOMedIso_ != NULL)) delete RECOPtvsGENdRRECOMedIso_;
  RECOPtvsGENdRRECOMedIso_ = NULL;
  if ((doDelete) && (MatchedTauHadEtavsTauHadPt_ != NULL)) delete MatchedTauHadEtavsTauHadPt_;
  MatchedTauHadEtavsTauHadPt_ = NULL;
  if ((doDelete) && (BDiscCSVvsRawIso_ != NULL)) delete BDiscCSVvsRawIso_;
  BDiscCSVvsRawIso_ = NULL;

  if ((doDelete) && (MatchedLooseIsoRECOPt_ != NULL)) delete MatchedLooseIsoRECOPt_;
  MatchedLooseIsoRECOPt_ = NULL;
  if ((doDelete) && (MatchedMedIsoRECOPt_ != NULL)) delete MatchedMedIsoRECOPt_;
  MatchedMedIsoRECOPt_ = NULL;
  if ((doDelete) && (MatchedTightIsoRECOPt_ != NULL)) delete MatchedTightIsoRECOPt_;
  MatchedTightIsoRECOPt_ = NULL;
  if ((doDelete) && (MatchedDMFindRECOPt_ != NULL)) delete MatchedDMFindRECOPt_;
  MatchedDMFindRECOPt_ = NULL;
  if ((doDelete) && (MatchedRECOPt_ != NULL)) delete MatchedRECOPt_;
  MatchedRECOPt_ = NULL;
  if ((doDelete) && (FinalEffLooseIsoRECOPt_ != NULL)) delete FinalEffLooseIsoRECOPt_;
  FinalEffLooseIsoRECOPt_ = NULL;
  if ((doDelete) && (FinalEffMedIsoRECOPt_ != NULL)) delete FinalEffMedIsoRECOPt_;
  FinalEffMedIsoRECOPt_ = NULL;
  if ((doDelete) && (FinalEffTightIsoRECOPt_ != NULL)) delete FinalEffTightIsoRECOPt_;
  FinalEffTightIsoRECOPt_ = NULL;
  if ((doDelete) && (FinalEffDMFindRECOPt_ != NULL)) delete FinalEffDMFindRECOPt_;
  FinalEffDMFindRECOPt_ = NULL;

  if ((doDelete) && (MatchedLooseIsoCJPt_ != NULL)) delete MatchedLooseIsoCJPt_;
  MatchedLooseIsoCJPt_ = NULL;
  if ((doDelete) && (MatchedMedIsoCJPt_ != NULL)) delete MatchedMedIsoCJPt_;
  MatchedMedIsoCJPt_ = NULL;
  if ((doDelete) && (MatchedTightIsoCJPt_ != NULL)) delete MatchedTightIsoCJPt_;
  MatchedTightIsoCJPt_ = NULL;
  if ((doDelete) && (MatchedDMFindCJPt_ != NULL)) delete MatchedDMFindCJPt_;
  MatchedDMFindCJPt_ = NULL;
  if ((doDelete) && (MatchedCJPt_ != NULL)) delete MatchedCJPt_;
  MatchedCJPt_ = NULL;
  if ((doDelete) && (FinalEffLooseIsoCJPt_ != NULL)) delete FinalEffLooseIsoCJPt_;
  FinalEffLooseIsoCJPt_ = NULL;
  if ((doDelete) && (FinalEffMedIsoCJPt_ != NULL)) delete FinalEffMedIsoCJPt_;
  FinalEffMedIsoCJPt_ = NULL;
  if ((doDelete) && (FinalEffTightIsoCJPt_ != NULL)) delete FinalEffTightIsoCJPt_;
  FinalEffTightIsoCJPt_ = NULL;
  if ((doDelete) && (FinalEffDMFindCJPt_ != NULL)) delete FinalEffDMFindCJPt_;
  FinalEffDMFindCJPt_ = NULL;


  if ((doDelete) && (MatchedLooseIsoRECOdR_ != NULL)) delete MatchedLooseIsoRECOdR_;
  MatchedLooseIsoRECOdR_ = NULL;
  if ((doDelete) && (MatchedMedIsoRECOdR_ != NULL)) delete MatchedMedIsoRECOdR_;
  MatchedMedIsoRECOdR_ = NULL;
  if ((doDelete) && (MatchedTightIsoRECOdR_ != NULL)) delete MatchedTightIsoRECOdR_;
  MatchedTightIsoRECOdR_ = NULL;
  if ((doDelete) && (MatchedDMFindRECOdR_ != NULL)) delete MatchedDMFindRECOdR_;
  MatchedDMFindRECOdR_ = NULL;
  if ((doDelete) && (MatchedRECOdR_ != NULL)) delete MatchedRECOdR_;
  MatchedRECOdR_ = NULL;
  if ((doDelete) && (FinalEffLooseIsoRECOdR_ != NULL)) delete FinalEffLooseIsoRECOdR_;
  FinalEffLooseIsoRECOdR_ = NULL;
  if ((doDelete) && (FinalEffMedIsoRECOdR_ != NULL)) delete FinalEffMedIsoRECOdR_;
  FinalEffMedIsoRECOdR_ = NULL;
  if ((doDelete) && (FinalEffTightIsoRECOdR_ != NULL)) delete FinalEffTightIsoRECOdR_;
  FinalEffTightIsoRECOdR_ = NULL;
  if ((doDelete) && (FinalEffDMFindRECOdR_ != NULL)) delete FinalEffDMFindRECOdR_;
  FinalEffDMFindRECOdR_ = NULL;

  if ((doDelete) && (MatchedLooseIsoCJdR_ != NULL)) delete MatchedLooseIsoCJdR_;
  MatchedLooseIsoCJdR_ = NULL;
  if ((doDelete) && (MatchedMedIsoCJdR_ != NULL)) delete MatchedMedIsoCJdR_;
  MatchedMedIsoCJdR_ = NULL;
  if ((doDelete) && (MatchedTightIsoCJdR_ != NULL)) delete MatchedTightIsoCJdR_;
  MatchedTightIsoCJdR_ = NULL;
  if ((doDelete) && (MatchedDMFindCJdR_ != NULL)) delete MatchedDMFindCJdR_;
  MatchedDMFindCJdR_ = NULL;
  if ((doDelete) && (MatchedCJdR_ != NULL)) delete MatchedCJdR_;
  MatchedCJdR_ = NULL;
  if ((doDelete) && (FinalEffLooseIsoCJdR_ != NULL)) delete FinalEffLooseIsoCJdR_;
  FinalEffLooseIsoCJdR_ = NULL;
  if ((doDelete) && (FinalEffMedIsoCJdR_ != NULL)) delete FinalEffMedIsoCJdR_;
  FinalEffMedIsoCJdR_ = NULL;
  if ((doDelete) && (FinalEffTightIsoCJdR_ != NULL)) delete FinalEffTightIsoCJdR_;
  FinalEffTightIsoCJdR_ = NULL;
  if ((doDelete) && (FinalEffDMFindCJdR_ != NULL)) delete FinalEffDMFindCJdR_;
  FinalEffDMFindCJdR_ = NULL;

  if ((doDelete) && (MatchedLooseIsoRECOdRPtCut_ != NULL)) delete MatchedLooseIsoRECOdRPtCut_;
  MatchedLooseIsoRECOdRPtCut_ = NULL;
  if ((doDelete) && (MatchedMedIsoRECOdRPtCut_ != NULL)) delete MatchedMedIsoRECOdRPtCut_;
  MatchedMedIsoRECOdRPtCut_ = NULL;
  if ((doDelete) && (MatchedTightIsoRECOdRPtCut_ != NULL)) delete MatchedTightIsoRECOdRPtCut_;
  MatchedTightIsoRECOdRPtCut_ = NULL;
  if ((doDelete) && (MatchedDMFindRECOdRPtCut_ != NULL)) delete MatchedDMFindRECOdRPtCut_;
  MatchedDMFindRECOdRPtCut_ = NULL;
  if ((doDelete) && (MatchedRECOdRPtCut_ != NULL)) delete MatchedRECOdRPtCut_;
  MatchedRECOdRPtCut_ = NULL;
  if ((doDelete) && (FinalEffLooseIsoRECOdRPtCut_ != NULL)) delete FinalEffLooseIsoRECOdRPtCut_;
  FinalEffLooseIsoRECOdRPtCut_ = NULL;
  if ((doDelete) && (FinalEffMedIsoRECOdRPtCut_ != NULL)) delete FinalEffMedIsoRECOdRPtCut_;
  FinalEffMedIsoRECOdRPtCut_ = NULL;
  if ((doDelete) && (FinalEffTightIsoRECOdRPtCut_ != NULL)) delete FinalEffTightIsoRECOdRPtCut_;
  FinalEffTightIsoRECOdRPtCut_ = NULL;
  if ((doDelete) && (FinalEffDMFindRECOdRPtCut_ != NULL)) delete FinalEffDMFindRECOdRPtCut_;
  FinalEffDMFindRECOdRPtCut_ = NULL;

  if ((doDelete) && (MatchedLooseIsoCJdRPtCut_ != NULL)) delete MatchedLooseIsoCJdRPtCut_;
  MatchedLooseIsoCJdRPtCut_ = NULL;
  if ((doDelete) && (MatchedMedIsoCJdRPtCut_ != NULL)) delete MatchedMedIsoCJdRPtCut_;
  MatchedMedIsoCJdRPtCut_ = NULL;
  if ((doDelete) && (MatchedTightIsoCJdRPtCut_ != NULL)) delete MatchedTightIsoCJdRPtCut_;
  MatchedTightIsoCJdRPtCut_ = NULL;
  if ((doDelete) && (MatchedDMFindCJdRPtCut_ != NULL)) delete MatchedDMFindCJdRPtCut_;
  MatchedDMFindCJdRPtCut_ = NULL;
  if ((doDelete) && (MatchedCJdRPtCut_ != NULL)) delete MatchedCJdRPtCut_;
  MatchedCJdRPtCut_ = NULL;
  if ((doDelete) && (FinalEffLooseIsoCJdRPtCut_ != NULL)) delete FinalEffLooseIsoCJdRPtCut_;
  FinalEffLooseIsoCJdRPtCut_ = NULL;
  if ((doDelete) && (FinalEffMedIsoCJdRPtCut_ != NULL)) delete FinalEffMedIsoCJdRPtCut_;
  FinalEffMedIsoCJdRPtCut_ = NULL;
  if ((doDelete) && (FinalEffTightIsoCJdRPtCut_ != NULL)) delete FinalEffTightIsoCJdRPtCut_;
  FinalEffTightIsoCJdRPtCut_ = NULL;
  if ((doDelete) && (FinalEffDMFindCJdRPtCut_ != NULL)) delete FinalEffDMFindCJdRPtCut_;
  FinalEffDMFindCJdRPtCut_ = NULL;

  if ((doDelete) && (MatchedLooseIsoRECOPtGen_ != NULL)) delete MatchedLooseIsoRECOPtGen_;
  MatchedLooseIsoRECOPtGen_ = NULL;
  if ((doDelete) && (MatchedMedIsoRECOPtGen_ != NULL)) delete MatchedMedIsoRECOPtGen_;
  MatchedMedIsoRECOPtGen_ = NULL;
  if ((doDelete) && (MatchedTightIsoRECOPtGen_ != NULL)) delete MatchedTightIsoRECOPtGen_;
  MatchedTightIsoRECOPtGen_ = NULL;
  if ((doDelete) && (MatchedDMFindRECOPtGen_ != NULL)) delete MatchedDMFindRECOPtGen_;
  MatchedDMFindRECOPtGen_ = NULL;
  if ((doDelete) && (MatchedRECOPtGen_ != NULL)) delete MatchedRECOPtGen_;
  MatchedRECOPtGen_ = NULL;
  if ((doDelete) && (FinalEffLooseIsoRECOPtGen_ != NULL)) delete FinalEffLooseIsoRECOPtGen_;
  FinalEffLooseIsoRECOPtGen_ = NULL;
  if ((doDelete) && (FinalEffMedIsoRECOPtGen_ != NULL)) delete FinalEffMedIsoRECOPtGen_;
  FinalEffMedIsoRECOPtGen_ = NULL;
  if ((doDelete) && (FinalEffTightIsoRECOPtGen_ != NULL)) delete FinalEffTightIsoRECOPtGen_;
  FinalEffTightIsoRECOPtGen_ = NULL;
  if ((doDelete) && (FinalEffDMFindRECOPtGen_ != NULL)) delete FinalEffDMFindRECOPtGen_;
  FinalEffDMFindRECOPtGen_ = NULL;

  if ((doDelete) && (OneProngDMCJPt_ != NULL)) delete OneProngDMCJPt_;
  OneProngDMCJPt_ = NULL;
  if ((doDelete) && (OneProngOnePizDMCJPt_ != NULL)) delete OneProngOnePizDMCJPt_;
  OneProngOnePizDMCJPt_ = NULL;
  if ((doDelete) && (OneProngTwoPizDMCJPt_ != NULL)) delete OneProngTwoPizDMCJPt_;
  OneProngOnePizDMCJPt_ = NULL;
  if ((doDelete) && (ThreeProngDMCJPt_ != NULL)) delete ThreeProngDMCJPt_;
  ThreeProngDMCJPt_ = NULL;
  if ((doDelete) && (OneProngDMRECOPt_ != NULL)) delete OneProngDMRECOPt_;
  OneProngDMRECOPt_ = NULL;
  if ((doDelete) && (OneProngOnePizDMRECOPt_ != NULL)) delete OneProngOnePizDMRECOPt_;
  OneProngOnePizDMRECOPt_ = NULL;
  if ((doDelete) && (OneProngTwoPizDMRECOPt_ != NULL)) delete OneProngTwoPizDMRECOPt_;
  OneProngTwoPizDMRECOPt_ = NULL;
  if ((doDelete) && (ThreeProngDMRECOPt_ != NULL)) delete ThreeProngDMRECOPt_;
  ThreeProngDMRECOPt_ = NULL;
  if ((doDelete) && (MatchedOneProngCJPt_ != NULL)) delete MatchedOneProngCJPt_;
  MatchedOneProngCJPt_ = NULL;
  if ((doDelete) && (MatchedOneProngOnePizCJPt_ != NULL)) delete MatchedOneProngOnePizCJPt_;
  MatchedOneProngOnePizCJPt_ = NULL;
  if ((doDelete) && (MatchedOneProngTwoPizCJPt_ != NULL)) delete MatchedOneProngTwoPizCJPt_;
  MatchedOneProngTwoPizCJPt_ = NULL;
  if ((doDelete) && (MatchedThreeProngCJPt_ != NULL)) delete MatchedThreeProngCJPt_;
  MatchedThreeProngCJPt_ = NULL;
  if ((doDelete) && (MatchedOneProngRECOPt_ != NULL)) delete MatchedOneProngRECOPt_;
  MatchedOneProngRECOPt_ = NULL;
  if ((doDelete) && (MatchedOneProngOnePizRECOPt_ != NULL)) delete MatchedOneProngOnePizRECOPt_;
  MatchedOneProngOnePizRECOPt_ = NULL;
  if ((doDelete) && (MatchedOneProngTwoPizRECOPt_ != NULL)) delete MatchedOneProngTwoPizRECOPt_;
  MatchedOneProngTwoPizRECOPt_ = NULL;
  if ((doDelete) && (MatchedThreeProngRECOPt_ != NULL)) delete MatchedThreeProngRECOPt_;
  MatchedThreeProngRECOPt_ = NULL;
  if ((doDelete) && (FinalOneProngDMCJPt_ != NULL)) delete FinalOneProngDMCJPt_;
  FinalOneProngDMCJPt_ = NULL;
  if ((doDelete) && (FinalOneProngOnePizDMCJPt_ != NULL)) delete FinalOneProngOnePizDMCJPt_;
  FinalOneProngOnePizDMCJPt_ = NULL;
  if ((doDelete) && (FinalOneProngTwoPizDMCJPt_ != NULL)) delete FinalOneProngTwoPizDMCJPt_;
  FinalOneProngTwoPizDMCJPt_ = NULL;
  if ((doDelete) && (FinalThreeProngDMCJPt_ != NULL)) delete FinalThreeProngDMCJPt_;
  FinalThreeProngDMCJPt_ = NULL;
  if ((doDelete) && (FinalOneProngDMRECOPt_ != NULL)) delete FinalOneProngDMRECOPt_;
  FinalOneProngDMRECOPt_ = NULL;
  if ((doDelete) && (FinalOneProngOnePizDMRECOPt_ != NULL)) delete FinalOneProngOnePizDMRECOPt_;
  FinalOneProngOnePizDMRECOPt_ = NULL;
  if ((doDelete) && (FinalOneProngTwoPizDMRECOPt_ != NULL)) delete FinalOneProngTwoPizDMRECOPt_;
  FinalOneProngTwoPizDMRECOPt_ = NULL;
  if ((doDelete) && (FinalThreeProngDMRECOPt_ != NULL)) delete FinalThreeProngDMRECOPt_;
  FinalThreeProngDMRECOPt_ = NULL;
  if ((doDelete) && (FinalEffOneProngDMSAMEPt_ != NULL)) delete FinalEffOneProngDMSAMEPt_;
  FinalEffOneProngDMSAMEPt_ = NULL;
  if ((doDelete) && (FinalEffOneProngOnePizDMSAMEPt_ != NULL)) delete FinalEffOneProngOnePizDMSAMEPt_;
  FinalEffOneProngOnePizDMSAMEPt_ = NULL;
  if ((doDelete) && (FinalEffOneProngTwoPizDMSAMEPt_ != NULL)) delete FinalEffOneProngTwoPizDMSAMEPt_;
  FinalEffOneProngTwoPizDMSAMEPt_ = NULL;
  if ((doDelete) && (FinalEffThreeProngDMSAMEPt_ != NULL)) delete FinalEffThreeProngDMSAMEPt_;
  FinalEffThreeProngDMSAMEPt_ = NULL;


  if ((doDelete) && (MatchedLooseIsoCJPtGen_ != NULL)) delete MatchedLooseIsoCJPtGen_;
  MatchedLooseIsoCJPtGen_ = NULL;
  if ((doDelete) && (MatchedMedIsoCJPtGen_ != NULL)) delete MatchedMedIsoCJPtGen_;
  MatchedMedIsoCJPtGen_ = NULL;
  if ((doDelete) && (MatchedTightIsoCJPtGen_ != NULL)) delete MatchedTightIsoCJPtGen_;
  MatchedTightIsoCJPtGen_ = NULL;
  if ((doDelete) && (MatchedDMFindCJPtGen_ != NULL)) delete MatchedDMFindCJPtGen_;
  MatchedDMFindCJPtGen_ = NULL;
  if ((doDelete) && (MatchedCJPtGen_ != NULL)) delete MatchedCJPtGen_;
  MatchedCJPtGen_ = NULL;
  if ((doDelete) && (FinalEffLooseIsoCJPtGen_ != NULL)) delete FinalEffLooseIsoCJPtGen_;
  FinalEffLooseIsoCJPtGen_ = NULL;
  if ((doDelete) && (FinalEffMedIsoCJPtGen_ != NULL)) delete FinalEffMedIsoCJPtGen_;
  FinalEffMedIsoCJPtGen_ = NULL;
  if ((doDelete) && (FinalEffTightIsoCJPtGen_ != NULL)) delete FinalEffTightIsoCJPtGen_;
  FinalEffTightIsoCJPtGen_ = NULL;
  if ((doDelete) && (FinalEffDMFindCJPtGen_ != NULL)) delete FinalEffDMFindCJPtGen_;
  FinalEffDMFindCJPtGen_ = NULL;


}//void FakeRateAnalyzer

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void FakeRateAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FakeRateAnalyzer);

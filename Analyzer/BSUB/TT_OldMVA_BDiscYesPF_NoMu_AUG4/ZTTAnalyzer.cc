// -*- C++ -*-
//
// Package:    ZTTAnalyzer
// Class:      ZTTAnalyzer
// 
/**\class ZTTAnalyzer ZTTAnalyzer.cc Analyzer/src/ZTTAnalyzer.cc

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

class ZTTAnalyzer : public edm::EDAnalyzer {
   public:
      typedef reco::JetFloatAssociation::Container JetBCEnergyRatioCollection;
      explicit ZTTAnalyzer(const edm::ParameterSet&);
      ~ZTTAnalyzer();

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
      edm::EDGetTokenT<vector<reco::PFJet> > akJetTag_;
      edm::EDGetTokenT<vector<reco::Muon> > muonsTag_;
      edm::EDGetTokenT<vector<reco::PFTau> > tauRECOTag_;
      edm::EDGetTokenT<vector<reco::RecoTauPiZero> > pizerosTag_;      
      edm::EDGetTokenT<reco::PFTauDiscriminator>  looseIsoTagRECO_;
      edm::EDGetTokenT<reco::PFTauDiscriminator>  medIsoTagRECO_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> tightIsoTagRECO_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> decayModeFindingTagRECO_;


      //Histograms
      TH1F* NEvents_;   
      TH1F* NMuRemoved_;
      TH1F* TauMuTauHaddR_;
      TH1F* NTauDecayModeGEN_;
      TH2F* NTausRECOvsCLEANJETS_;
      TH2F* TauHadnConstvsPt_;
      TH2F* RECOPtvsGENVisiblePtRECO_;
      TH2F* RECOPtvsGENdRRECO_;
      TH2F* RECOPtvsGENdRRECOMedIso_;


      TH1F* MatchedLooseIsoRECOPt_;
      TH1F* MatchedMedIsoRECOPt_;
      TH1F* MatchedTightIsoRECOPt_;
      TH1F* MatchedDMFindRECOPt_;
      TH1F* MatchedRECOPt_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOPt_;
      TGraphAsymmErrors* FinalEffMedIsoRECOPt_;
      TGraphAsymmErrors* FinalEffTightIsoRECOPt_;
      TGraphAsymmErrors* FinalEffDMFindRECOPt_;

      TH1F* MatchedLooseIsoRECOdR_;
      TH1F* MatchedMedIsoRECOdR_;
      TH1F* MatchedTightIsoRECOdR_;
      TH1F* MatchedDMFindRECOdR_;
      TH1F* MatchedRECOdR_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOdR_;
      TGraphAsymmErrors* FinalEffMedIsoRECOdR_;
      TGraphAsymmErrors* FinalEffTightIsoRECOdR_;
      TGraphAsymmErrors* FinalEffDMFindRECOdR_;

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
      TH1F* OneProngDMRECOPt_;
      TH1F* OneProngOnePizDMRECOPt_;
      TH1F* OneProngTwoPizDMRECOPt_;
      TH1F* ThreeProngDMRECOPt_;
      TH1F* MatchedOneProngRECOPt_;
      TH1F* MatchedOneProngOnePizRECOPt_;
      TH1F* MatchedOneProngTwoPizRECOPt_;
      TH1F* MatchedThreeProngRECOPt_;
      TGraphAsymmErrors* FinalOneProngDMRECOPt_;
      TGraphAsymmErrors* FinalOneProngOnePizDMRECOPt_;
      TGraphAsymmErrors* FinalOneProngTwoPizDMRECOPt_;
      TGraphAsymmErrors* FinalThreeProngDMRECOPt_;
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
ZTTAnalyzer::ZTTAnalyzer(const edm::ParameterSet& iConfig):
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  genParticleTag_(consumes<vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleTag"))),
  akJetTag_(consumes<vector<reco::PFJet> >(iConfig.getParameter<edm::InputTag>("akJetTag"))),
  muonsTag_(consumes<vector<reco::Muon> >(iConfig.getParameter<edm::InputTag>("muonsTag"))),
  tauRECOTag_(consumes<vector<reco::PFTau> >(iConfig.getParameter<edm::InputTag>("tauRECOTag"))),
  pizerosTag_(consumes<vector<reco::RecoTauPiZero> >(iConfig.getParameter<edm::InputTag>("pizerosTag"))),
  looseIsoTagRECO_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("looseIsoTagRECO"))),
  medIsoTagRECO_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("medIsoTagRECO"))),
  tightIsoTagRECO_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("tightIsoTagRECO"))),
  decayModeFindingTagRECO_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("decayModeFindingTagRECO")))
{
  reset(false);    
}//ZTTAnalyzer



ZTTAnalyzer::~ZTTAnalyzer()
{
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void ZTTAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "\n<------------THIS IS A NEW EVENT------------>" << std::endl;
  NEvents_->Fill(0);

  //Get gen particle collection
  edm::Handle<reco::GenParticleCollection> pGenParts;
  iEvent.getByToken(genParticleTag_, pGenParts);

  //Get ak4Jets particle collection
  edm::Handle<std::vector<reco::PFJet> > pAkJets;
  iEvent.getByToken(akJetTag_, pAkJets);

  //Get RECO Muons particle collection
  edm::Handle<std::vector<reco::Muon> > pMuons;
  iEvent.getByToken(muonsTag_, pMuons);

  //Get RECO Taus particle collection
  edm::Handle<std::vector<reco::PFTau> > pTausRECO;
  iEvent.getByToken(tauRECOTag_, pTausRECO);

  //Get pi-Zero particle collection
  edm::Handle<std::vector<reco::RecoTauPiZero> > pPiZero;
  iEvent.getByToken(pizerosTag_, pPiZero);

  //Get Loose Iso Collection
  Handle<PFTauDiscriminator> pLooseIsoDiscRECO;
  iEvent.getByToken(looseIsoTagRECO_, pLooseIsoDiscRECO);

  //Get Medium Iso Collection
  Handle<PFTauDiscriminator> pMedIsoDiscRECO;
  iEvent.getByToken(medIsoTagRECO_, pMedIsoDiscRECO);

  //Get Tight Iso Collection
  Handle<PFTauDiscriminator> pTightIsoDiscRECO;
  iEvent.getByToken(tightIsoTagRECO_, pTightIsoDiscRECO);

  //Get Decay Mode Finding Collection
  Handle<PFTauDiscriminator> pDMFindingRECO;
  iEvent.getByToken(decayModeFindingTagRECO_, pDMFindingRECO);

//////////////////////////////
// Begin Analyzer
//////////////////////////////
  reco::GenParticleRef tau1Ref, tau2Ref;
  bool checkDY = false;
  for (reco::GenParticleCollection::const_iterator iGenParticle = pGenParts->begin(); iGenParticle != pGenParts->end(); ++iGenParticle)
  {
    if((iGenParticle->pdgId() == 23 || iGenParticle->pdgId() == 22) && iGenParticle->numberOfDaughters() == 2 && fabs(iGenParticle->daughter(0)->pdgId() ) == 15 )
    {
      checkDY = true;
      tau1Ref = iGenParticle->daughterRef(0);
      tau2Ref = iGenParticle->daughterRef(1);
     
      //This finds the status 2 taus by seeing if a tau is the daughter of the current tau ref
      while(VariousFunctions::findIfInDaughters(tau1Ref, 15, true) || VariousFunctions::findIfInDaughters(tau2Ref, 15, true))
      {
        if(VariousFunctions::findIfInDaughters(tau1Ref, 15, true))
          tau1Ref = VariousFunctions::findDaughterInDaughters(tau1Ref, 15, true);
        if(VariousFunctions::findIfInDaughters(tau2Ref, 15, true))
          tau2Ref = VariousFunctions::findDaughterInDaughters(tau2Ref, 15, true);
      }//while
    }//ir found a1 and it decayed to taus
  }//for iGen Particle

  if (!checkDY)
    return;
  std::cout << "chekcDY= " << checkDY << std::endl;
  double dPhi_gen = reco::deltaPhi(tau1Ref->phi(), tau2Ref->phi() ), dR_tauMu_gen = sqrt( (tau1Ref->eta() - tau2Ref->eta())*(tau1Ref->eta() - tau2Ref->eta())  +  dPhi_gen * dPhi_gen );
  reco::LeafCandidate::LorentzVector GenTau1Visible = VariousFunctions::sumTauP4(tau1Ref, VariousFunctions::tauDecayMode(tau1Ref), false); 
  reco::LeafCandidate::LorentzVector GenTau2Visible = VariousFunctions::sumTauP4(tau2Ref, VariousFunctions::tauDecayMode(tau2Ref), false);  

//////////////////////////
// Matching to RECO Taus
//////////////////////////
  double dR_GenMatchRECO1 = 1000, dR_GenMatchRECO2 = 1000, dR_GenMatchRECOTEMP1 = 1000, dR_GenMatchRECOTEMP2 = 1000, Pt_GenMatchRECO1 = -1, Pt_GenMatchRECO2 = -1;
  unsigned int DMRECO1 = 0, DMRECO2 = 0, MedIsoRECO1 = 0, MedIsoRECO2 = 0, LooseIsoRECO1 = 0, LooseIsoRECO2 = 0, TightIsoRECO1 = 0, TightIsoRECO2 = 0;
  for (std::vector<reco::PFTau>::const_iterator iTauRECO = pTausRECO->begin(); iTauRECO != pTausRECO->end(); ++iTauRECO)
  {
    //Gen Match the CleanJets taus
    double dPhi1 = reco::deltaPhi(tau1Ref->phi(), iTauRECO->phi() ), dPhi2 = reco::deltaPhi(tau2Ref->phi(), iTauRECO->phi() );
    dR_GenMatchRECOTEMP1 = sqrt( (tau1Ref->eta() - iTauRECO->eta() ) * (tau1Ref->eta() - iTauRECO->eta() )  +  dPhi1 * dPhi1);
    dR_GenMatchRECOTEMP2 = sqrt( (tau2Ref->eta() - iTauRECO->eta() ) * (tau2Ref->eta() - iTauRECO->eta() )  +  dPhi2 * dPhi2);
    std::cout << "dR_GenMatchRECOTEMP1= " << dR_GenMatchRECOTEMP1 << "\tdR_GenMatchRECOTEMP2= " << dR_GenMatchRECOTEMP2 << std::endl;
    reco::PFTauRef PFTauRef(pTausRECO, iTauRECO - pTausRECO->begin());
    //first 2 checks are basic requirements, the 2 in the "or" are making sure it doesn't matche better with the second tau
    if (dR_GenMatchRECOTEMP1 < .1 && dR_GenMatchRECOTEMP1 < dR_GenMatchRECO1 && (dR_GenMatchRECOTEMP1 < dR_GenMatchRECOTEMP2 || dR_GenMatchRECOTEMP2 > dR_GenMatchRECO2) )
    {
      dR_GenMatchRECO1 = dR_GenMatchRECOTEMP1;
      Pt_GenMatchRECO1 = iTauRECO->pt();
      DMRECO1 = (*pDMFindingRECO)[PFTauRef];
      MedIsoRECO1 = (*pMedIsoDiscRECO)[PFTauRef];
      LooseIsoRECO1 = (*pLooseIsoDiscRECO)[PFTauRef];
      TightIsoRECO1 = (*pTightIsoDiscRECO)[PFTauRef];
    }//if iTauRECO matches tau1Ref 
    else if (dR_GenMatchRECOTEMP2 < .1 && dR_GenMatchRECOTEMP2 < dR_GenMatchRECO2 && (dR_GenMatchRECOTEMP2 < dR_GenMatchRECOTEMP1 || dR_GenMatchRECOTEMP1 > dR_GenMatchRECO1) )
    {
      dR_GenMatchRECO2 = dR_GenMatchRECOTEMP2;
      Pt_GenMatchRECO2 = iTauRECO->pt();
      DMRECO2 = (*pDMFindingRECO)[PFTauRef];
      MedIsoRECO2 = (*pMedIsoDiscRECO)[PFTauRef];
      LooseIsoRECO2 = (*pLooseIsoDiscRECO)[PFTauRef]; 
      TightIsoRECO2 = (*pTightIsoDiscRECO)[PFTauRef];
    }//else

  }//iTauRECO

  if (Pt_GenMatchRECO2 > 20)
  {
std::cout << "\t\tPASSED GEN MATCH and MU REMOVAL" << std::endl;
    RECOPtvsGENVisiblePtRECO_->Fill(Pt_GenMatchRECO2, GenTau2Visible.Pt() );
    NEvents_->Fill(4);
    RECOPtvsGENdRRECO_->Fill(Pt_GenMatchRECO2, dR_tauMu_gen);
    if (MedIsoRECO2 == 1 && DMRECO2 == 1)
      RECOPtvsGENdRRECOMedIso_->Fill(Pt_GenMatchRECO2, dR_tauMu_gen);
    MatchedRECOPt_->Fill(Pt_GenMatchRECO2 );
    if (DMRECO2 == 1)
      MatchedDMFindRECOPt_->Fill(Pt_GenMatchRECO2 );
    if (TightIsoRECO2 == 1 && DMRECO2 == 1)
      MatchedTightIsoRECOPt_->Fill(Pt_GenMatchRECO2 );
    if (MedIsoRECO2 == 1 && DMRECO2 == 1)
      MatchedMedIsoRECOPt_->Fill(Pt_GenMatchRECO2 );
    if (LooseIsoRECO2 == 1 && DMRECO2 == 1)
      MatchedLooseIsoRECOPt_->Fill(Pt_GenMatchRECO2 );

    MatchedRECOdR_->Fill(dR_tauMu_gen );
    if (DMRECO2 == 1)
      MatchedDMFindRECOdR_->Fill(dR_tauMu_gen );
    if (TightIsoRECO2 == 1 && DMRECO2 == 1)
      MatchedTightIsoRECOdR_->Fill(dR_tauMu_gen );
    if (MedIsoRECO2 == 1 && DMRECO2 == 1)
      MatchedMedIsoRECOdR_->Fill(dR_tauMu_gen );
    if (LooseIsoRECO2 == 1 && DMRECO2 == 1)
      MatchedLooseIsoRECOdR_->Fill(dR_tauMu_gen );

    if (Pt_GenMatchRECO2 > 20) 
    {
      MatchedRECOdRPtCut_->Fill(dR_tauMu_gen );
      if (DMRECO2 == 1)
        MatchedDMFindRECOdRPtCut_->Fill(dR_tauMu_gen );
      if (TightIsoRECO2 == 1 && DMRECO2 == 1)
        MatchedTightIsoRECOdRPtCut_->Fill(dR_tauMu_gen );
      if (MedIsoRECO2 == 1 && DMRECO2 == 1)
        MatchedMedIsoRECOdRPtCut_->Fill(dR_tauMu_gen );
      if (LooseIsoRECO2 == 1 && DMRECO2 == 1)
        MatchedLooseIsoRECOdRPtCut_->Fill(dR_tauMu_gen );
    }//Pt_GenMatchRECO2 > 20

    MatchedRECOPtGen_->Fill(GenTau2Visible.Pt() );
    if (DMRECO2 == 1)
      MatchedDMFindRECOPtGen_->Fill(GenTau2Visible.Pt() );
    if (TightIsoRECO2 == 1 && DMRECO2 == 1)
      MatchedTightIsoRECOPtGen_->Fill(GenTau2Visible.Pt() );
    if (MedIsoRECO2 == 1 && DMRECO2 == 1)
      MatchedMedIsoRECOPtGen_->Fill(GenTau2Visible.Pt() );
    if (LooseIsoRECO2 == 1 && DMRECO2 == 1)
      MatchedLooseIsoRECOPtGen_->Fill(GenTau2Visible.Pt() );

  }//if GEN tau2Ref is the had in mu+had and it is matched to CleanJets Jet

  else if (Pt_GenMatchRECO1 > 20)
  {
std::cout << "\t\tPASSED GEN MATCH and MU REMOVAL" << std::endl;
   RECOPtvsGENVisiblePtRECO_->Fill(Pt_GenMatchRECO1, GenTau1Visible.Pt() );
   NEvents_->Fill(4);
   RECOPtvsGENdRRECO_->Fill(Pt_GenMatchRECO1, dR_tauMu_gen);
   if (MedIsoRECO1 == 1 && DMRECO1 == 1)
     RECOPtvsGENdRRECOMedIso_->Fill(Pt_GenMatchRECO1, dR_tauMu_gen);
   MatchedRECOPt_->Fill(Pt_GenMatchRECO1 );
   if (DMRECO1 == 1)
     MatchedDMFindRECOPt_->Fill(Pt_GenMatchRECO1 );
   if (TightIsoRECO1 == 1 && DMRECO1 == 1)
     MatchedTightIsoRECOPt_->Fill(Pt_GenMatchRECO1 );
   if (MedIsoRECO1 == 1 && DMRECO1 == 1)
     MatchedMedIsoRECOPt_->Fill(Pt_GenMatchRECO1 );
   if (LooseIsoRECO1 == 1 && DMRECO1 == 1)
     MatchedLooseIsoRECOPt_->Fill(Pt_GenMatchRECO1 );
   
   MatchedRECOdR_->Fill(dR_tauMu_gen );
   if (DMRECO1 == 1)
     MatchedDMFindRECOdR_->Fill(dR_tauMu_gen );
   if (TightIsoRECO1 == 1 && DMRECO1 == 1)
     MatchedTightIsoRECOdR_->Fill(dR_tauMu_gen );
   if (MedIsoRECO1 == 1 && DMRECO1 == 1)
     MatchedMedIsoRECOdR_->Fill(dR_tauMu_gen );
   if (LooseIsoRECO1 == 1 && DMRECO1 == 1)
     MatchedLooseIsoRECOdR_->Fill(dR_tauMu_gen );

    if (Pt_GenMatchRECO1 > 20)
    {
      MatchedRECOdRPtCut_->Fill(dR_tauMu_gen );
      if (DMRECO1 == 1)
        MatchedDMFindRECOdRPtCut_->Fill(dR_tauMu_gen );
      if (TightIsoRECO1 == 1 && DMRECO1 == 1)
        MatchedTightIsoRECOdRPtCut_->Fill(dR_tauMu_gen );
      if (MedIsoRECO1 == 1 && DMRECO1 == 1)
        MatchedMedIsoRECOdRPtCut_->Fill(dR_tauMu_gen );
      if (LooseIsoRECO1 == 1 && DMRECO1 == 1)
        MatchedLooseIsoRECOdRPtCut_->Fill(dR_tauMu_gen );
    }//Pt_GenMatchRECO1 > 20

   
   MatchedRECOPtGen_->Fill(GenTau1Visible.Pt() );
   if (DMRECO1 == 1)
     MatchedDMFindRECOPtGen_->Fill(GenTau1Visible.Pt() );
   if (TightIsoRECO1 == 1 && DMRECO1 == 1)
     MatchedTightIsoRECOPtGen_->Fill(GenTau1Visible.Pt() );
   if (MedIsoRECO1 == 1 && DMRECO1 == 1)
     MatchedMedIsoRECOPtGen_->Fill(GenTau1Visible.Pt() );
   if (LooseIsoRECO1 == 1 && DMRECO1 == 1)
     MatchedLooseIsoRECOPtGen_->Fill(GenTau1Visible.Pt() );
  }//if GEN tau1Ref is the had in mu+had and it is matched to CleanJets Jet

 
}//End ZTTAnalyzer::analyze


// ------------ method called once each job just before starting event loop  ------------
void ZTTAnalyzer::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  //Open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //Book histograms
  NEvents_     = new TH1F("NEvents"    , "", 9, -.5, 8.5);
      NEvents_->GetXaxis()->SetBinLabel(1, "TotalEvents"); 
      NEvents_->GetXaxis()->SetBinLabel(2, "#tau_{#mu} + #tau_{had} Match");
      NEvents_->GetXaxis()->SetBinLabel(3, "Gen #tau_{#mu} + #tau_{had}");
      NEvents_->GetXaxis()->SetBinLabel(5, "RECO Gen Match #tau_{had}");
      NEvents_->GetXaxis()->SetBinLabel(6, "Event with #tau_{#mu} Removed");
      NEvents_->GetXaxis()->SetBinLabel(7, "Event with no #tau_{#mu} Removed ");
  NMuRemoved_            = new TH1F("NMuRemoved"    , "", 2, -.5, 1.5);
      NMuRemoved_->GetXaxis()->SetBinLabel(1, "#tau_{#mu} + #tau_{Had} Event");
      NMuRemoved_->GetXaxis()->SetBinLabel(2, "and a Removed #mu");
  TauMuTauHaddR_       = new TH1F("TauMuTauHaddR"    , "", 100, 0, 100);
  NTauDecayModeGEN_        = new TH1F("NTauDecayModeGEN"    , "", 12, -.5, 11.5);
  NTausRECOvsCLEANJETS_  = new TH2F("NTausRECOvsCLEANJETS" , "", 11, -.5, 10.5, 11, -.5, 10.5);
  TauHadnConstvsPt_  = new TH2F("TauHadnConstvsPt" , "", 50, 0, 50, 30, 0, 220);
  RECOPtvsGENVisiblePtRECO_  = new TH2F("RECOPtvsGENVisiblePtRECO" , "", 50, 0, 220, 50, 0, 220);
  RECOPtvsGENdRRECO_ = new TH2F("RECOPtvsGENdRRECO", "", 40, 0, 220, 40, 0, 4);
  RECOPtvsGENdRRECOMedIso_ = new TH2F("RECOPtvsGENdRRECOMedIso", "", 40, 0, 220, 40, 0, 4);


  MatchedLooseIsoRECOPt_  = new TH1F("MatchedLooseIsoRECOPt"    , "", 11, 0, 220.0);
  MatchedMedIsoRECOPt_    = new TH1F("MatchedMedIsoRECOPt", "", 11, 0, 220);
  MatchedTightIsoRECOPt_    = new TH1F("MatchedTightIsoRECOPt", "", 11, 0, 220);
  MatchedDMFindRECOPt_    = new TH1F("MatchedDMFindRECOPt"    , "", 11, 0, 220);
  MatchedRECOPt_          = new TH1F("MatchedRECOPt"    , "", 11, 0, 220);

  FinalEffLooseIsoRECOPt_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOPt_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOPt_ = new TGraphAsymmErrors(30);
  FinalEffDMFindRECOPt_ = new TGraphAsymmErrors(30);

  MatchedLooseIsoRECOdR_  = new TH1F("MatchedLooseIsoRECOdR"    , "", 25, 0, 5);
  MatchedMedIsoRECOdR_    = new TH1F("MatchedMedIsoRECOdR", "", 25, 0, 5);
  MatchedTightIsoRECOdR_    = new TH1F("MatchedTightIsoRECOdR", "", 25, 0, 5);
  MatchedDMFindRECOdR_    = new TH1F("MatchedDMFindRECOdR"    , "", 25, 0, 5);
  MatchedRECOdR_          = new TH1F("MatchedRECOdR"    , "", 25, 0, 5);

  FinalEffLooseIsoRECOdR_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOdR_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOdR_ = new TGraphAsymmErrors(30);
  FinalEffDMFindRECOdR_ = new TGraphAsymmErrors(30);

  MatchedLooseIsoRECOdRPtCut_  = new TH1F("MatchedLooseIsoRECOdRPtCut"    , "", 25, 0, 5);
  MatchedMedIsoRECOdRPtCut_    = new TH1F("MatchedMedIsoRECOdRPtCut", "", 25, 0, 5);
  MatchedTightIsoRECOdRPtCut_    = new TH1F("MatchedTightIsoRECOdRPtCut", "", 25, 0, 5);
  MatchedDMFindRECOdRPtCut_    = new TH1F("MatchedDMFindRECOdRPtCut"    , "", 25, 0, 5);
  MatchedRECOdRPtCut_          = new TH1F("MatchedRECOdRPtCut"    , "", 25, 0, 5);

  FinalEffLooseIsoRECOdRPtCut_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOdRPtCut_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOdRPtCut_ = new TGraphAsymmErrors(30);
  FinalEffDMFindRECOdRPtCut_ = new TGraphAsymmErrors(30);

  MatchedLooseIsoRECOPtGen_  = new TH1F("MatchedLooseIsoRECOPtGen"    , "", 11, 0, 220.0);
  MatchedMedIsoRECOPtGen_    = new TH1F("MatchedMedIsoRECOPtGen", "", 11, 0, 220);
  MatchedTightIsoRECOPtGen_    = new TH1F("MatchedTightIsoRECOPtGen", "", 11, 0, 220);
  MatchedDMFindRECOPtGen_    = new TH1F("MatchedDMFindRECOPtGen"    , "", 11, 0, 220);
  MatchedRECOPtGen_          = new TH1F("MatchedRECOPtGen"    , "", 11, 0, 220);

  FinalEffLooseIsoRECOPtGen_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOPtGen_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOPtGen_ = new TGraphAsymmErrors(30);
  FinalEffDMFindRECOPtGen_ = new TGraphAsymmErrors(30);
  OneProngDMRECOPt_ = new TH1F("OneProngDMRECOPt"    , "", 11, 0, 220.0);
  OneProngOnePizDMRECOPt_ = new TH1F("OneProngOnePizDMRECOPt"    , "", 11, 0, 220.0);
  OneProngTwoPizDMRECOPt_ = new TH1F("OneProngTwoPizDMRECOPt"    , "", 11, 0, 220.0);
  ThreeProngDMRECOPt_ = new TH1F("ThreeProngDMRECOPt"    , "", 11, 0, 220.0);
  MatchedOneProngRECOPt_ = new TH1F("MatchedOneProngRECOPt", "", 11, 0, 220.0);
  MatchedOneProngOnePizRECOPt_ = new TH1F("MatchedOneProngOnePizRECOPt", "", 11, 0, 220.0);
  MatchedOneProngTwoPizRECOPt_ = new TH1F("MatchedOneProngTwoPizRECOPt", "", 11, 0, 220.0);
  MatchedThreeProngRECOPt_ = new TH1F("MatchedThreeProngRECOPt", "", 11, 0, 220.0);
  FinalOneProngDMRECOPt_ = new TGraphAsymmErrors(30);
  FinalOneProngOnePizDMRECOPt_ = new TGraphAsymmErrors(30);
  FinalOneProngTwoPizDMRECOPt_ = new TGraphAsymmErrors(30);
  FinalThreeProngDMRECOPt_ = new TGraphAsymmErrors(30);
}

// ------------ method called once each job just after ending the event loop  ------------
void ZTTAnalyzer::endJob()
{
  //Make the Canvases
  TCanvas NEventsCanvas("NEvents","",600,600);
  TCanvas NMuRemovedCanvas("NMuRemoved","",600,600);
  TCanvas TauMuTauHaddRCanvas("TauMuTauHaddR","",600,600);
  TCanvas NTauDecayModeGENCanvas("NTauDecayModeGEN","",600,600);
  TCanvas NTausRECOvsCLEANJETSCanvas("NTausRECOvsCLEANJETS","",600,600);
  TCanvas TauHadnConstvsPtCanvas("TauHadnConstvsPt","",600,600);
  TCanvas RECOPtvsGENVisiblePtRECOCanvas("RECOPtvsGENVisiblePtRECO","",600,600);
  TCanvas RECOPtvsGENdRRECOCanvas("RECOPtvsGENdRRECO","",600,600);
  TCanvas RECOPtvsGENdRRECOMedIsoCanvas("RECOPtvsGENdRRECOMedIso","",600,600);

  TCanvas MatchedLooseIsoRECOPtCanvas("MatchedLooseIsoRECOPt","",600,600);
  TCanvas MatchedMedIsoRECOPtCanvas("MatchedMedIsoRECOPt","",600,600);
  TCanvas MatchedTightIsoRECOPtCanvas("MatchedTightIsoRECOPt","",600,600);
  TCanvas MatchedDMFindRECOPtCanvas("MatchedDMFindRECOPt","",600,600);
  TCanvas MatchedRECOPtCanvas("MatchedRECOPt","",600,600);
  TCanvas FinalEffLooseIsoRECOPtCanvas("FinalEffLooseIsoRECOPt","",600,600);
  TCanvas FinalEffMedIsoRECOPtCanvas("FinalEffMedIsoRECOPt","",600,600);
  TCanvas FinalEffTightIsoRECOPtCanvas("FinalEffTightIsoRECOPt","",600,600);
  TCanvas FinalEffDMFindRECOPtCanvas("FinalEffDMFindRECOPt","",600,600);

  TCanvas MatchedLooseIsoRECOdRCanvas("MatchedLooseIsoRECOdR","",600,600);
  TCanvas MatchedMedIsoRECOdRCanvas("MatchedMedIsoRECOdR","",600,600);
  TCanvas MatchedTightIsoRECOdRCanvas("MatchedTightIsoRECOdR","",600,600);
  TCanvas MatchedDMFindRECOdRCanvas("MatchedDMFindRECOdR","",600,600);
  TCanvas MatchedRECOdRCanvas("MatchedRECOdR","",600,600);
  TCanvas FinalEffLooseIsoRECOdRCanvas("FinalEffLooseIsoRECOdR","",600,600);
  TCanvas FinalEffMedIsoRECOdRCanvas("FinalEffMedIsoRECOdR","",600,600);
  TCanvas FinalEffTightIsoRECOdRCanvas("FinalEffTightIsoRECOdR","",600,600);
  TCanvas FinalEffDMFindRECOdRCanvas("FinalEffDMFindRECOdR","",600,600);

  TCanvas MatchedLooseIsoRECOdRPtCutCanvas("MatchedLooseIsoRECOdRPtCut","",600,600);
  TCanvas MatchedMedIsoRECOdRPtCutCanvas("MatchedMedIsoRECOdRPtCut","",600,600);
  TCanvas MatchedTightIsoRECOdRPtCutCanvas("MatchedTightIsoRECOdRPtCut","",600,600);
  TCanvas MatchedDMFindRECOdRPtCutCanvas("MatchedDMFindRECOdRPtCut","",600,600);
  TCanvas MatchedRECOdRPtCutCanvas("MatchedRECOdRPtCut","",600,600);
  TCanvas FinalEffLooseIsoRECOdRPtCutCanvas("FinalEffLooseIsoRECOdRPtCut","",600,600);
  TCanvas FinalEffMedIsoRECOdRPtCutCanvas("FinalEffMedIsoRECOdRPtCut","",600,600);
  TCanvas FinalEffTightIsoRECOdRPtCutCanvas("FinalEffTightIsoRECOdRPtCut","",600,600);
  TCanvas FinalEffDMFindRECOdRPtCutCanvas("FinalEffDMFindRECOdRPtCut","",600,600);

  TCanvas MatchedLooseIsoRECOPtGenCanvas("MatchedLooseIsoRECOPtGen","",600,600);
  TCanvas MatchedMedIsoRECOPtGenCanvas("MatchedMedIsoRECOPtGen","",600,600);
  TCanvas MatchedTightIsoRECOPtGenCanvas("MatchedTightIsoRECOPtGen","",600,600);
  TCanvas MatchedDMFindRECOPtGenCanvas("MatchedDMFindRECOPtGen","",600,600);
  TCanvas MatchedRECOPtGenCanvas("MatchedRECOPtGen","",600,600);
  TCanvas FinalEffLooseIsoRECOPtGenCanvas("FinalEffLooseIsoRECOPtGen","",600,600);
  TCanvas FinalEffMedIsoRECOPtGenCanvas("FinalEffMedIsoRECOPtGen","",600,600);
  TCanvas FinalEffTightIsoRECOPtGenCanvas("FinalEffTightIsoRECOPtGen","",600,600);
  TCanvas FinalEffDMFindRECOPtGenCanvas("FinalEffDMFindRECOPtGen","",600,600);
  TCanvas OneProngDMRECOPtCanvas("OneProngDMRECOPt","",600,600);
  TCanvas OneProngOnePizDMRECOPtCanvas("OneProngOnePizDMRECOPt","",600,600);
  TCanvas OneProngTwoPizDMRECOPtCanvas("OneProngTwoPizDMRECOPt","",600,600);
  TCanvas ThreeProngDMRECOPtCanvas("ThreeProngDMRECOPt","",600,600);
  TCanvas MatchedOneProngRECOPtCanvas("MatchedOneProngRECOPtCanvas","",600,600);
  TCanvas MatchedOneProngOnePizRECOPtCanvas("MatchedOneProngOnePizRECOPtCanvas","",600,600);
  TCanvas MatchedOneProngTwoPizRECOPtCanvas("MatchedOneProngTwoPizRECOPtCanvas","",600,600);
  TCanvas MatchedThreeProngRECOPtCanvas("MatchedThreeProngRECOPtCanvas","",600,600);
  TCanvas FinalOneProngDMRECOPtCanvas("FinalOneProngDMRECOPt","",600,600);
  TCanvas FinalOneProngOnePizDMRECOPtCanvas("FinalOneProngOnePizDMRECOPt","",600,600);
  TCanvas FinalOneProngTwoPizDMRECOPtCanvas("FinalOneProngTwoPizDMRECOPt","",600,600);
  TCanvas FinalThreeProngDMRECOPtCanvas("FinalThreeProngDMRECOPt","",600,600);


  //I CHANGED THIS SO THAT IT IS OVER GEM-MATCH + DM. i ALSO CHANGED THE LABELS
  FinalEffLooseIsoRECOPt_->Divide(MatchedLooseIsoRECOPt_, MatchedDMFindRECOPt_);
  FinalEffMedIsoRECOPt_->Divide(MatchedMedIsoRECOPt_,     MatchedDMFindRECOPt_);
  FinalEffTightIsoRECOPt_->Divide(MatchedTightIsoRECOPt_, MatchedDMFindRECOPt_);
  FinalEffDMFindRECOPt_->Divide(MatchedDMFindRECOPt_,     MatchedDMFindRECOPt_);

  FinalEffLooseIsoRECOdR_->Divide(MatchedLooseIsoRECOdR_, MatchedDMFindRECOdR_);
  FinalEffMedIsoRECOdR_->Divide(MatchedMedIsoRECOdR_,     MatchedDMFindRECOdR_);
  FinalEffTightIsoRECOdR_->Divide(MatchedTightIsoRECOdR_, MatchedDMFindRECOdR_);
  FinalEffDMFindRECOdR_->Divide(MatchedDMFindRECOdR_,     MatchedDMFindRECOdR_);

  FinalEffLooseIsoRECOdRPtCut_->Divide(MatchedLooseIsoRECOdRPtCut_, MatchedDMFindRECOdRPtCut_);
  FinalEffMedIsoRECOdRPtCut_->Divide(MatchedMedIsoRECOdRPtCut_,     MatchedDMFindRECOdRPtCut_);
  FinalEffTightIsoRECOdRPtCut_->Divide(MatchedTightIsoRECOdRPtCut_, MatchedDMFindRECOdRPtCut_);
  FinalEffDMFindRECOdRPtCut_->Divide(MatchedDMFindRECOdRPtCut_,     MatchedDMFindRECOdRPtCut_);

  FinalEffLooseIsoRECOPtGen_->Divide(MatchedLooseIsoRECOPtGen_, MatchedDMFindRECOPtGen_);
  FinalEffMedIsoRECOPtGen_->Divide(MatchedMedIsoRECOPtGen_,     MatchedDMFindRECOPtGen_);
  FinalEffTightIsoRECOPtGen_->Divide(MatchedTightIsoRECOPtGen_, MatchedDMFindRECOPtGen_);
  FinalEffDMFindRECOPtGen_->Divide(MatchedDMFindRECOPtGen_,     MatchedDMFindRECOPtGen_);

  FinalOneProngDMRECOPt_->Divide(OneProngDMRECOPt_, MatchedOneProngRECOPt_);
  FinalOneProngOnePizDMRECOPt_->Divide(OneProngOnePizDMRECOPt_, MatchedOneProngOnePizRECOPt_);
  FinalOneProngTwoPizDMRECOPt_->Divide(OneProngTwoPizDMRECOPt_, MatchedOneProngTwoPizRECOPt_);
  FinalThreeProngDMRECOPt_->Divide(ThreeProngDMRECOPt_, MatchedThreeProngRECOPt_);

std::cout << "<----------------Declared Canvases-------------->" << std::endl;

  //Format the 1D plots and Draw (canvas, hist, grid, log y, log z, color, size, style, xAxisTitle, xTitleSize, xLabelSize, xTitleOffSet, yAxisTitle, yTitleSize, yLabelSize, yTitleOffset)
  VariousFunctions::formatAndDrawCanvasAndHist1D(NEventsCanvas, NEvents_,
	 1, 0, 0, kBlack, 7, 20, "", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NMuRemovedCanvas, NMuRemoved_,
	 1, 0, 0, kBlack, 7, 20, "N #mu Removed", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauMuTauHaddRCanvas, TauMuTauHaddR_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(#tau_{mu} + #tau_{H})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NTauDecayModeGENCanvas, NTauDecayModeGEN_,
         1, 0, 0, kBlack, 7, 20, "TauDecayMode GEN ", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist2D(NTausRECOvsCLEANJETSCanvas, NTausRECOvsCLEANJETS_,
	 1, 0, 0, kBlack, 7, 20, "nRECO #tau's", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(TauHadnConstvsPtCanvas, TauHadnConstvsPt_,
         1, 0, 0, kBlack, 7, 20, "# Constituents(#tau_{H})", .04, .04, 1.1, "Pt(#tau_{H})", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(RECOPtvsGENVisiblePtRECOCanvas, RECOPtvsGENVisiblePtRECO_,
         1, 0, 0, kBlack, 7, 20, "RECO Pt(#tau)", .04, .04, 1.1, "GEN Visible Pt(#tau)", .04, .04, 1.6, "", .04, .04, 1.0);



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
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngDMRECOPtCanvas, OneProngDMRECOPt_, 
         1, 0, 0, kBlack, 7, 20, "1 Prong p_{T}(RECO Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngOnePizDMRECOPtCanvas, OneProngOnePizDMRECOPt_, 
         1, 0, 0, kBlack, 7, 20, "1 Prong + 1 #pi^{0} p_{T}(RECO Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngTwoPizDMRECOPtCanvas, OneProngTwoPizDMRECOPt_, 
         1, 0, 0, kBlack, 7, 20, "1 Prong + 2 #pi^{0} p_{T}(RECO Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(ThreeProngDMRECOPtCanvas, ThreeProngDMRECOPt_, 
         1, 0, 0, kBlack, 7, 20, "3 Prong p_{T}(RECO Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedOneProngRECOPtCanvas, MatchedOneProngRECOPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong p_{T}(RECO Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedOneProngOnePizRECOPtCanvas, MatchedOneProngOnePizRECOPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 1 #pi^{0} p_{T}(RECO Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedOneProngTwoPizRECOPtCanvas, MatchedOneProngTwoPizRECOPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 2 #pi^{0} p_{T}(RECO Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedThreeProngRECOPtCanvas, MatchedThreeProngRECOPt_,
         1, 0, 0, kBlack, 7, 20, "3 Prong p_{T}(RECO Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngDMRECOPtCanvas, FinalOneProngDMRECOPt_, 
         1, 1, 1, kBlack, 1, 20, "1 Prong p_{T}(RECO Matched + DM) / 1 Prong p_{T}(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngOnePizDMRECOPtCanvas, FinalOneProngOnePizDMRECOPt_, 
         1, 1, 1, kBlack, 1, 20, "1 Prong + 1 #pi^{0} p_{T}(RECO Matched + DM) / 1 Prong + 1 #pi^{0} p_{T}(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngTwoPizDMRECOPtCanvas, FinalOneProngTwoPizDMRECOPt_, 
         1, 1, 1, kBlack, 1, 20, "1 Prong + 2 #pi^{0} p_{T}(RECO Matched + DM) / 1 Prong + 2 #pi^{0} p_{T}(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalThreeProngDMRECOPtCanvas, FinalThreeProngDMRECOPt_, 
         1, 1, 1, kBlack, 1, 20, "3 Prong p_{T}(RECO Matched + DM) / 3 Prong p_{T}(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

std::cout << "after formatting" << std::endl;
std::cout << "check before setting axis" << std::endl;

std::cout << "<----------------Formatted Canvases and Histos-------------->" << std::endl;

  //Write output file
  out_->cd();

  NEventsCanvas.Write();
  NMuRemovedCanvas.Write();
  TauMuTauHaddRCanvas.Write();
  NTauDecayModeGENCanvas.Write();
  NTausRECOvsCLEANJETSCanvas.Write();
  TauHadnConstvsPtCanvas.Write();
  RECOPtvsGENVisiblePtRECOCanvas.Write();
  RECOPtvsGENdRRECOCanvas.Write();
  RECOPtvsGENdRRECOMedIsoCanvas.Write();

  MatchedLooseIsoRECOPtCanvas.Write();
  MatchedMedIsoRECOPtCanvas.Write();
  MatchedTightIsoRECOPtCanvas.Write();
  MatchedDMFindRECOPtCanvas.Write();
  MatchedRECOPtCanvas.Write();
  FinalEffLooseIsoRECOPtCanvas.Write();
  FinalEffMedIsoRECOPtCanvas.Write();
  FinalEffTightIsoRECOPtCanvas.Write();
  FinalEffDMFindRECOPtCanvas.Write();

  MatchedLooseIsoRECOdRCanvas.Write();
  MatchedMedIsoRECOdRCanvas.Write();
  MatchedTightIsoRECOdRCanvas.Write();
  MatchedDMFindRECOdRCanvas.Write();
  MatchedRECOdRCanvas.Write();
  FinalEffLooseIsoRECOdRCanvas.Write();
  FinalEffMedIsoRECOdRCanvas.Write();
  FinalEffTightIsoRECOdRCanvas.Write();
  FinalEffDMFindRECOdRCanvas.Write();

  MatchedLooseIsoRECOdRPtCutCanvas.Write();
  MatchedMedIsoRECOdRPtCutCanvas.Write();
  MatchedTightIsoRECOdRPtCutCanvas.Write();
  MatchedDMFindRECOdRPtCutCanvas.Write();
  MatchedRECOdRPtCutCanvas.Write();
  FinalEffLooseIsoRECOdRPtCutCanvas.Write();
  FinalEffMedIsoRECOdRPtCutCanvas.Write();
  FinalEffTightIsoRECOdRPtCutCanvas.Write();
  FinalEffDMFindRECOdRPtCutCanvas.Write();

  MatchedLooseIsoRECOPtGenCanvas.Write();
  MatchedMedIsoRECOPtGenCanvas.Write();
  MatchedTightIsoRECOPtGenCanvas.Write();
  MatchedDMFindRECOPtGenCanvas.Write();
  MatchedRECOPtGenCanvas.Write();
  FinalEffLooseIsoRECOPtGenCanvas.Write();
  FinalEffMedIsoRECOPtGenCanvas.Write();
  FinalEffTightIsoRECOPtGenCanvas.Write();
  FinalEffDMFindRECOPtGenCanvas.Write();
  OneProngDMRECOPtCanvas.Write();
  OneProngOnePizDMRECOPtCanvas.Write();
  OneProngTwoPizDMRECOPtCanvas.Write();
  ThreeProngDMRECOPtCanvas.Write();
  MatchedOneProngRECOPtCanvas.Write();
  MatchedOneProngOnePizRECOPtCanvas.Write();
  MatchedOneProngTwoPizRECOPtCanvas.Write();
  MatchedThreeProngRECOPtCanvas.Write();
  FinalOneProngDMRECOPtCanvas.Write();
  FinalOneProngOnePizDMRECOPtCanvas.Write();
  FinalOneProngTwoPizDMRECOPtCanvas.Write();
  FinalThreeProngDMRECOPtCanvas.Write();

  out_->Write();
  out_->Close();
std::cout << "DONE" << std::endl;
}//EndJob

// ------------ method called when starting to processes a run  ------------
void ZTTAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void ZTTAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void ZTTAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void ZTTAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

//Delete Memory
void ZTTAnalyzer::reset(const bool doDelete)
{
  if ((doDelete) && (NEvents_ != NULL)) delete NEvents_;
  NEvents_ = NULL;
  if ((doDelete) && (NMuRemoved_ != NULL)) delete NMuRemoved_;
  NMuRemoved_ = NULL;
  if ((doDelete) && (TauMuTauHaddR_ != NULL)) delete TauMuTauHaddR_;
  TauMuTauHaddR_ = NULL;
  if ((doDelete) && (NTauDecayModeGEN_ != NULL)) delete NTauDecayModeGEN_;
  NTauDecayModeGEN_ = NULL;
  if ((doDelete) && (NTausRECOvsCLEANJETS_ != NULL)) delete NTausRECOvsCLEANJETS_;
  NTausRECOvsCLEANJETS_ = NULL;
  if ((doDelete) && (TauHadnConstvsPt_ != NULL)) delete TauHadnConstvsPt_;
  TauHadnConstvsPt_ = NULL;
  if ((doDelete) && (RECOPtvsGENVisiblePtRECO_ != NULL)) delete RECOPtvsGENVisiblePtRECO_;
  RECOPtvsGENVisiblePtRECO_ = NULL;
  if ((doDelete) && (RECOPtvsGENdRRECO_ != NULL)) delete RECOPtvsGENdRRECO_;
  RECOPtvsGENdRRECO_ = NULL;
  if ((doDelete) && (RECOPtvsGENdRRECOMedIso_ != NULL)) delete RECOPtvsGENdRRECOMedIso_;
  RECOPtvsGENdRRECOMedIso_ = NULL;

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
  if ((doDelete) && (OneProngDMRECOPt_ != NULL)) delete OneProngDMRECOPt_;
  OneProngDMRECOPt_ = NULL;
  if ((doDelete) && (OneProngOnePizDMRECOPt_ != NULL)) delete OneProngOnePizDMRECOPt_;
  OneProngOnePizDMRECOPt_ = NULL;
  if ((doDelete) && (OneProngTwoPizDMRECOPt_ != NULL)) delete OneProngTwoPizDMRECOPt_;
  OneProngTwoPizDMRECOPt_ = NULL;
  if ((doDelete) && (ThreeProngDMRECOPt_ != NULL)) delete ThreeProngDMRECOPt_;
  ThreeProngDMRECOPt_ = NULL;
  if ((doDelete) && (MatchedOneProngRECOPt_ != NULL)) delete MatchedOneProngRECOPt_;
  MatchedOneProngRECOPt_ = NULL;
  if ((doDelete) && (MatchedOneProngOnePizRECOPt_ != NULL)) delete MatchedOneProngOnePizRECOPt_;
  MatchedOneProngOnePizRECOPt_ = NULL;
  if ((doDelete) && (MatchedOneProngTwoPizRECOPt_ != NULL)) delete MatchedOneProngTwoPizRECOPt_;
  MatchedOneProngTwoPizRECOPt_ = NULL;
  if ((doDelete) && (MatchedThreeProngRECOPt_ != NULL)) delete MatchedThreeProngRECOPt_;
  MatchedThreeProngRECOPt_ = NULL;
  if ((doDelete) && (FinalOneProngDMRECOPt_ != NULL)) delete FinalOneProngDMRECOPt_;
  FinalOneProngDMRECOPt_ = NULL;
  if ((doDelete) && (FinalOneProngOnePizDMRECOPt_ != NULL)) delete FinalOneProngOnePizDMRECOPt_;
  FinalOneProngOnePizDMRECOPt_ = NULL;
  if ((doDelete) && (FinalOneProngTwoPizDMRECOPt_ != NULL)) delete FinalOneProngTwoPizDMRECOPt_;
  FinalOneProngTwoPizDMRECOPt_ = NULL;
  if ((doDelete) && (FinalThreeProngDMRECOPt_ != NULL)) delete FinalThreeProngDMRECOPt_;
  FinalThreeProngDMRECOPt_ = NULL;


}//void ZTTAnalyzer

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ZTTAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZTTAnalyzer);

// -*- C++ -*-
//
// Package:    BBAAnalyzer
// Class:      BBAAnalyzer
// 
/**\class BBAAnalyzer BBAAnalyzer.cc Analyzer/src/BBAAnalyzer.cc

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
#include <map>
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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "HLTrigger/HLTcore/interface/HLTConfigData.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTanalyzers/interface/HLTInfo.h"

#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"

#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "AnalyzerGeneratorRecoVariousFunctions/VariousFunctions/interface/VariousFunctions.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"


using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;


//
// class declaration
//

class BBAAnalyzer : public edm::EDAnalyzer {
   public:
      typedef reco::JetFloatAssociation::Container JetBCEnergyRatioCollection;
      explicit BBAAnalyzer(const edm::ParameterSet&);
      ~BBAAnalyzer();

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

      //gen particle Tag
      edm::InputTag prunedGenParticleTag_;

      //PAT jet Tag
      edm::InputTag slimmedJetTag_;

      //PAT TAU Tag
      edm::InputTag slimmedTauTag_;

      //PAT Muon Tag
      edm::EDGetTokenT< pat::MuonCollection > slimmedMuonTag_;

      //PAT Step 2 Tau Tag
      edm::InputTag hpsTauTag_;

      // Mu BTagCSV Muon Tag
      std::string  muonMatchMu_;
      
      //Non  Mu BTagCSV Muon Tag     
      std::string  muonMatchNonMu_;
      
      //Non  Mu BTagCSV bJet Tag
      std::string  bMatchNonMu_;
      
      // Mu BTagCSV bJet Tag
      std::string  bMatchMu_;

      // input for patTrigger
      edm::InputTag trigger_;

      // input for patTriggerEvent
      edm::EDGetTokenT< pat::TriggerEvent > triggerEventToken_;

      //Histograms
      TH1F* HLTPass_;
      TH1F* NonIsoMuPt_;
      TH1F* IsoMuPt_;
      TH1F* NonIsoBPt_;
      TH1F* IsoBPt_;
      TH1F* DiTaudR_Ptgt40_;
      TH1F* DiTaudR_Ptgt60_;
      TH1F* DiTaudR_Ptgt80_;
      TH1F* DiTaudR_Ptgt100_;

      TH1F* NMu_;
      TH2F* NSlimTauvsNhpsTau_;
      TH2F* NonIsoTrigPtvsRecoPt_;
      TH2F* IsoTrigPtvsRecoPt_;
      TH2F* NonIsoTrigEtavsRecoEta_;

      TH1F* NonIsoMatchDiTauDr_;
      TH1F* IsoMatchDiTauDr_;

      TH1F* NonIsoBCSVTurnOn_;
      TH1F* IsoBCSVTurnOn_;
      TH1F* NonIsoBPtTurnOn_;
      TH1F* IsoBPtTurnOn_;

      TH1F* NonIsoBCSVPassAll_;
      TH1F* IsoBCSVPassAll_;
      TH1F* NonIsoBPtPassAll_;
      TH1F* IsoBPtPassAll_;

      TH1F* NonIsoBCSVPassTrig_;
      TH1F* IsoBCSVPassTrig_;
      TH1F* NonIsoBPtPassTrig_;
      TH1F* IsoBPtPassTrig_;

      TH1F* NonIsoBCSVTurnOnMu_;
      TH1F* IsoBCSVTurnOnMu_;
      TH1F* NonIsoBPtTurnOnMu_;
      TH1F* IsoBPtTurnOnMu_;

      TH1F* NonIsoBCSVPassAllMu_;
      TH1F* IsoBCSVPassAllMu_;
      TH1F* NonIsoBPtPassAllMu_;
      TH1F* IsoBPtPassAllMu_;

      TH1F* NonIsoBCSVPassTrigMu_;
      TH1F* IsoBCSVPassTrigMu_;
      TH1F* NonIsoBPtPassTrigMu_;
      TH1F* IsoBPtPassTrigMu_;

      TH1F* TrigMatchdR_;
      TH1F* CSVlt54Pt_;
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
BBAAnalyzer::BBAAnalyzer(const edm::ParameterSet& iConfig):
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  prunedGenParticleTag_(iConfig.getParameter<edm::InputTag>("prunedGenParticleTag")),
  slimmedJetTag_(iConfig.getParameter<edm::InputTag>("slimmedJetTag")),
  slimmedTauTag_(iConfig.getParameter<edm::InputTag>("slimmedTauTag")),
  slimmedMuonTag_( consumes< pat::MuonCollection >( iConfig.getParameter< edm::InputTag >( "slimmedMuonTag" ) ) ),
  hpsTauTag_(iConfig.getParameter<edm::InputTag>("hpsTauTag")),
  muonMatchMu_( iConfig.getParameter< std::string >( "muonMatchIsoMu" ) ),
  muonMatchNonMu_( iConfig.getParameter< std::string >( "muonMatchNonIsoMu" ) ),
  bMatchNonMu_( iConfig.getParameter< std::string >( "bMatchNonIsoMu" ) ),
  bMatchMu_( iConfig.getParameter< std::string >( "bMatchIsoMu" ) ),
  trigger_( iConfig.getParameter< edm::InputTag >( "trigger" ) ),
  triggerEventToken_( consumes< pat::TriggerEvent >( iConfig.getParameter< edm::InputTag >( "triggerEvent" ) ) )
{
  reset(false);    
}//BBAAnalyzer



BBAAnalyzer::~BBAAnalyzer()
{
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void BBAAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "\n<------------THIS IS A NEW EVENT------------>" << std::endl;
  HLTPass_->Fill(0);
///////////////////////////////////////////////////////
//  Declare Event Handles
///////////////////////////////////////////////////////

  // slimmedJet handle
  edm::Handle<std::vector<pat::Jet> > pSlimmedJetTag;
  iEvent.getByLabel(slimmedJetTag_, pSlimmedJetTag);

  // slimmedTau handle
  edm::Handle<std::vector<pat::Tau> > pSlimmedTauTag;
  iEvent.getByLabel(slimmedTauTag_, pSlimmedTauTag);

  // slimmedTau handle
  edm::Handle<std::vector<reco::PFTau> > pHpsTauTag;
  iEvent.getByLabel(hpsTauTag_, pHpsTauTag);

  // PAT trigger event
  edm::Handle< pat::TriggerEvent > triggerEvent;
  iEvent.getByToken( triggerEventToken_, triggerEvent );

  // PAT object collection
  edm::Handle< pat::MuonCollection > pSlimmedMuonTag;
  iEvent.getByToken( slimmedMuonTag_, pSlimmedMuonTag );

  //Get gen particle collection
  edm::Handle<reco::GenParticleCollection> pPrunGenParts;
  iEvent.getByLabel(prunedGenParticleTag_, pPrunGenParts);

///////////////////////////////////////////////////////
// Get Gen Taus
///////////////////////////////////////////////////////
  reco::GenParticleCollection::const_iterator firstGen = pPrunGenParts->begin();
  reco::GenParticleRef tau1Ref = firstGen->daughterRef(0), tau2Ref = firstGen->daughterRef(0);
  for (reco::GenParticleCollection::const_iterator iGenParticle = pPrunGenParts->begin(); iGenParticle != pPrunGenParts->end(); ++iGenParticle)
  {
    if(iGenParticle->pdgId() == 36 && iGenParticle->numberOfDaughters() == 2)
    {
      tau1Ref = iGenParticle->daughterRef(0);
      tau2Ref = iGenParticle->daughterRef(1);
      break;
    }//if pdgid == 36
  }//for iGen Particle
/*
///////////////////////////////////////////////////////
// Slim Taus to Gen Taus 
///////////////////////////////////////////////////////
  double dR1_slimTau = 1000, dR2_slimTau = 1000, dR1_slimTauTemp = 0, dR2_slimTauTemp = 0;
  double slimTau1Pt = -1, slimTau2Pt = -1, slimTau1Eta = -10, slimTau2Eta = -10, slimTau1Phi = -10, slimTau2Phi = -10;
//  std::cout << "\tNumber of SlimTaus= " << pSlimmedTauTag->size() << std::endl;
  for (auto iTau = pSlimmedTauTag->begin(); iTau != pSlimmedTauTag->end(); ++iTau)
  {
    dR1_slimTauTemp = sqrt( (tau1Ref->eta() - iTau->eta() ) * (tau1Ref->eta() - iTau->eta() ) + (tau1Ref->phi() - iTau->phi() ) * (tau1Ref->phi() - iTau->phi() ) );
    dR2_slimTauTemp = sqrt( (tau2Ref->eta() - iTau->eta() ) * (tau2Ref->eta() - iTau->eta() ) + (tau2Ref->phi() - iTau->phi() ) * (tau2Ref->phi() - iTau->phi() ) );

//    std::cout << "\t\ttau1Ref: Eta= " << tau1Ref->eta() << "  Phi= " << tau1Ref->phi() << "  tau2Ref: Eta= " << tau2Ref->eta() << "  Phi= " << tau2Ref->phi() << std::endl;
//    std::cout << "\t\tiTau:    Eta= " << iTau->eta() << "  Phi= " << iTau->phi() << std::endl;
//    std::cout << "\t\tdR1_slimTau= " << dR1_slimTau << "  dR1_slimTauTemp= " << dR1_slimTauTemp << "  PtDiff1= " << abs(slimTau1Pt - iTau->pt() ) << std::endl;
//    std::cout << "\t\tdR2_slimTau= " << dR2_slimTau << "  dR2_slimTauTemp= " << dR2_slimTauTemp << "  PtDiff2= " << abs(slimTau2Pt - iTau->pt() ) << std::endl;

    //first 2 checks are basic requirements, the 2 in the "or" are making sure it doesn't matche better with the second tau
    if (dR1_slimTauTemp < .4 && dR1_slimTauTemp < dR1_slimTau && (dR1_slimTauTemp < dR2_slimTauTemp || dR2_slimTauTemp > dR2_slimTau) )
    {
      dR1_slimTau = dR1_slimTauTemp;
      slimTau1Pt = iTau->pt();
      slimTau1Eta = iTau->eta();
      slimTau1Phi = iTau->phi();
    }//if iTau matches tau1Ref
    else if (dR2_slimTauTemp < .4 && dR2_slimTauTemp < dR2_slimTau && (dR2_slimTauTemp < dR1_slimTauTemp || dR1_slimTauTemp > dR1_slimTau) )
    {
      dR2_slimTau = dR2_slimTauTemp;
      slimTau2Pt = iTau->pt();
      slimTau2Eta = iTau->eta();
      slimTau2Phi = iTau->phi();
    }//else
  }//for auto iTau
*/
///////////////////////////////////////////////////////
// Match PF Taus to Gen Taus 
///////////////////////////////////////////////////////
  double dR1_PFTau = 1000, dR2_PFTau = 1000, dR1_PFTauTemp = 0, dR2_PFTauTemp = 0;
  double PFTau1Pt = -1, PFTau2Pt = -1, PFTau1Eta = -10, PFTau2Eta = -10, PFTau1Phi = -10, PFTau2Phi = -10;
//  std::cout << "\tNumber of HPSTaus= " << pHpsTauTag->size() << std::endl;
  for (auto iTau = pHpsTauTag->begin(); iTau != pHpsTauTag->end(); ++iTau)
  {
    dR1_PFTauTemp = sqrt( (tau1Ref->eta() - iTau->eta() ) * (tau1Ref->eta() - iTau->eta() ) + (tau1Ref->phi() - iTau->phi() ) * (tau1Ref->phi() - iTau->phi() ) );
    dR2_PFTauTemp = sqrt( (tau2Ref->eta() - iTau->eta() ) * (tau2Ref->eta() - iTau->eta() ) + (tau2Ref->phi() - iTau->phi() ) * (tau2Ref->phi() - iTau->phi() ) );
    
    std::cout << "\t\ttau1Ref: Eta= " << tau1Ref->eta() << "  Phi= " << tau1Ref->phi() << "  tau2Ref: Eta= " << tau2Ref->eta() << "  Phi= " << tau2Ref->phi() << std::endl;
    std::cout << "\t\tiTau:    Eta= " << iTau->eta() << "  Phi= " << iTau->phi() << std::endl;
    std::cout << "\t\tdR1_PFTau= " << dR1_PFTau << "  dR1_PFTauTemp= " << dR1_PFTauTemp << "  PtDiff1= " << abs(PFTau1Pt - iTau->pt() ) << std::endl;
    std::cout << "\t\tdR2_PFTau= " << dR2_PFTau << "  dR2_PFTauTemp= " << dR2_PFTauTemp << "  PtDiff2= " << abs(PFTau2Pt - iTau->pt() ) << std::endl;
    
  //first 2 checks are basic requirements, the 2 in the "or" are making sure it doesn't matche better with the second tau
    if (dR1_PFTauTemp < .4 && dR1_PFTauTemp < dR1_PFTau && (dR1_PFTauTemp < dR2_PFTauTemp || dR2_PFTauTemp > dR2_PFTau) )
    { 
      dR1_PFTau = dR1_PFTauTemp;
      PFTau1Pt = iTau->pt();
      PFTau1Eta = iTau->eta();
      PFTau1Phi = iTau->phi();
    }//if iTau matches tau1Ref 
    else if (dR2_PFTauTemp < .4 && dR2_PFTauTemp < dR2_PFTau && (dR2_PFTauTemp < dR1_PFTauTemp || dR1_PFTauTemp > dR1_PFTau) )
    { 
      dR2_PFTau = dR2_PFTauTemp;
      PFTau2Pt = iTau->pt();
      PFTau2Eta = iTau->eta();
      PFTau2Phi = iTau->phi();
    }//else
  }//for auto iTau
  

  NSlimTauvsNhpsTau_->Fill(pSlimmedTauTag->size(), pHpsTauTag->size() );
  NMu_->Fill( pSlimmedMuonTag->size() );
std::cout << "\tPF:   dR1_PFTau= " << dR1_PFTau <<   "    dR1_PFTauTemp= " << dR1_PFTauTemp <<   "  Pt= " << PFTau1Pt <<   "  Eta= " << PFTau1Eta <<   "  Phi= " << PFTau1Phi << std::endl;
std::cout << "\tPF:   dR2_PFTau= " << dR2_PFTau <<   "    dR2_PFTauTemp= " << dR2_PFTauTemp <<   "  Pt= " << PFTau2Pt <<   "  Eta= " << PFTau2Eta <<   "  Phi= " << PFTau2Phi << std::endl;
  double diTauDr = sqrt( (PFTau1Eta - PFTau2Eta ) * (PFTau1Eta - PFTau2Eta ) + (PFTau1Phi - PFTau2Phi ) * (PFTau1Phi - PFTau2Phi ) );
std::cout << "\tPF:   DiTauDr= " << diTauDr << std::endl;

/*
///////////////////////////////////////////////////////////////////////////////////////////////
// Filling Base Histograms of all jets to Divide with for Trig Effeciency Study 
///////////////////////////////////////////////////////////////////////////////////////////////
std::cout << "\t!!!!!!!!!!!!! TRIGGER MATCHING !!!!!!!!!!!!!!" << std::endl;
  const pat::helper::TriggerMatchHelper matchHelper;
///////////////////////////////////////////////////////
// Loop for Non- Muon Trigger Matching
///////////////////////////////////////////////////////
  double dR_slimMuTau1 = 1000, dR_slimMuTau2 = 1000,  dR_PFMuTau1 = 1000, dR_PFMuTau2 = 1000;
  bool foundMatchAlready = false;
  for( size_t iMuon = 0; iMuon < pSlimmedMuonTag->size(); ++iMuon )
  {
    const pat::TriggerObjectRef trigRefM( matchHelper.triggerMatchObject( pSlimmedMuonTag, iMuon, muonMatchNonMu_, iEvent, *triggerEvent ) );
    if ( trigRefM.isAvailable() && trigRefM.isNonnull() )
    {
      cout << "<---------Non- Muon Trigger Match----------->" << endl; 
      if (!foundMatchAlready)
      {
        HLTPass_->Fill(1);
	foundMatchAlready = true;
      }// if !foundMatchAlready
      NonIsoMuPt_->Fill(trigRefM->pt() );
      dR_slimMuTau1 = sqrt( (trigRefM->eta() - slimTau1Eta)*(trigRefM->eta() - slimTau1Eta)  +  (trigRefM->phi() - slimTau1Phi)*(trigRefM->phi() - slimTau1Phi) );
      dR_slimMuTau2 = sqrt( (trigRefM->eta() - slimTau2Eta)*(trigRefM->eta() - slimTau2Eta)  +  (trigRefM->phi() - slimTau2Phi)*(trigRefM->phi() - slimTau2Phi) );
      dR_PFMuTau1 = sqrt( (trigRefM->eta() - PFTau1Eta)*(trigRefM->eta() - PFTau1Eta)  +  (trigRefM->phi() - PFTau1Phi)*(trigRefM->phi() - PFTau1Phi) );
      dR_PFMuTau2 = sqrt( (trigRefM->eta() - PFTau2Eta)*(trigRefM->eta() - PFTau2Eta)  +  (trigRefM->phi() - PFTau2Phi)*(trigRefM->phi() - PFTau2Phi) );
      std::cout << "\t\tdR_slimMuTau1= " << dR_slimMuTau1 << "  dR_slimMuTau2= " << dR_slimMuTau2 << " dR_PFMuTau1= " << dR_PFMuTau1 << "  dR_PFMuTau2= " << dR_PFMuTau2 << std::endl;
      //check for Non-B
      bool checkDiTau = false;
      for( size_t iJet = 0; iJet < pSlimmedJetTag->size(); ++iJet )    
      {
//        if (pSlimmedJetTag->at(iJet).bDiscriminator("pfCombinedSecondaryVertexV2BJetTags") > .7 && abs(pSlimmedJetTag->at(iJet).eta() ) < 2.6)
//          NonIsoBPtPassAllMu_->Fill(pSlimmedJetTag->at(iJet).pt() );
//        if (pSlimmedJetTag->at(iJet).pt() > 50 && abs(pSlimmedJetTag->at(iJet).eta() ) < 2.6)
//          NonIsoBCSVPassAllMu_->Fill(pSlimmedJetTag->at(iJet).bDiscriminator("pfCombinedSecondaryVertexV2BJetTags"));
 
        const pat::TriggerObjectRef trigRefB( matchHelper.triggerMatchObject( pSlimmedJetTag, iJet, bMatchNonMu_, iEvent, *triggerEvent ) );
        if ( trigRefB.isAvailable() && trigRefB.isNonnull() )
        { 
          cout << "\t<---------Non- B Trigger Match = FULL TRIGGER MATCH----------->" << endl;
//          if (pSlimmedJetTag->at(iJet).pt() > 50)
//            NonIsoBCSVPassTrigMu_->Fill( pSlimmedJetTag->at(iJet).bDiscriminator("pfCombinedSecondaryVertexV2BJetTags") );
//          if ( pSlimmedJetTag->at(iJet).bDiscriminator("pfCombinedSecondaryVertexV2BJetTags") > .7 )
//            NonIsoBPtPassTrigMu_->Fill( pSlimmedJetTag->at(iJet).pt() );

	  if (!checkDiTau && dR1_PFTau != 1000 && dR2_PFTau != 1000)
	  {
	    NonIsoMatchDiTauDr_->Fill(diTauDr );
	    checkDiTau = true;
	  }//if !checkDiTau
	}//iftrigRefB.isAvailable
        
      }//for iJet
    }//if
  }//for


///////////////////////////////////////////////////////
// Loop for  Muon Trigger Matching
///////////////////////////////////////////////////////
  dR_slimMuTau1 = 1000;
  dR_slimMuTau2 = 1000;
  dR_PFMuTau1 = 1000;
  dR_PFMuTau2 = 1000;
  foundMatchAlready = false;
  for( size_t iMuon = 0; iMuon < pSlimmedMuonTag->size(); ++iMuon ) 
  {
    const pat::TriggerObjectRef trigRefM( matchHelper.triggerMatchObject( pSlimmedMuonTag, iMuon, muonMatchMu_, iEvent, *triggerEvent ) );
    if ( trigRefM.isAvailable() && trigRefM.isNonnull() )
    {
      cout << "<-------- Muon Trigger Match----------->" << endl;

      if (!foundMatchAlready)
      {
        HLTPass_->Fill(3);
        foundMatchAlready = true;
      }// if !foundMatchAlready
      IsoMuPt_->Fill(trigRefM->pt() );
      dR_slimMuTau1 = sqrt( (trigRefM->eta() - slimTau1Eta)*(trigRefM->eta() - slimTau1Eta)  +  (trigRefM->phi() - slimTau1Phi)*(trigRefM->phi() - slimTau1Phi) );
      dR_slimMuTau2 = sqrt( (trigRefM->eta() - slimTau2Eta)*(trigRefM->eta() - slimTau2Eta)  +  (trigRefM->phi() - slimTau2Phi)*(trigRefM->phi() - slimTau2Phi) );
      dR_PFMuTau1 = sqrt( (trigRefM->eta() - PFTau1Eta)*(trigRefM->eta() - PFTau1Eta)  +  (trigRefM->phi() - PFTau1Phi)*(trigRefM->phi() - PFTau1Phi) );
      dR_PFMuTau2 = sqrt( (trigRefM->eta() - PFTau2Eta)*(trigRefM->eta() - PFTau2Eta)  +  (trigRefM->phi() - PFTau2Phi)*(trigRefM->phi() - PFTau2Phi) );
      std::cout << "\t\tdR_slimMuTau1= " << dR_slimMuTau1 << "  dR_slimMuTau2= " << dR_slimMuTau2 << " dR_PFMuTau1= " << dR_PFMuTau1 << "  dR_PFMuTau2= " << dR_PFMuTau2 << std::endl;
      bool checkDiTau = false;
      for( size_t iJet = 0; iJet < pSlimmedJetTag->size(); ++iJet )
      { 
//        if (pSlimmedJetTag->at(iJet).bDiscriminator("pfCombinedSecondaryVertexV2BJetTags") > .7 && abs(pSlimmedJetTag->at(iJet).eta() ) < 2.6)
//          IsoBPtPassAllMu_->Fill(pSlimmedJetTag->at(iJet).pt() );
//        if (pSlimmedJetTag->at(iJet).pt() > 50 && abs(pSlimmedJetTag->at(iJet).eta() ) < 2.6)
//          IsoBCSVPassAllMu_->Fill(pSlimmedJetTag->at(iJet).bDiscriminator("pfCombinedSecondaryVertexV2BJetTags"));

        const pat::TriggerObjectRef trigRefB( matchHelper.triggerMatchObject( pSlimmedJetTag, iJet, bMatchMu_, iEvent, *triggerEvent ) );
        if ( trigRefB.isAvailable() && trigRefB.isNonnull() )
        { 
          cout << "\t<--------- B Trigger Match = FULL TRIGGER MATCH----------->" << endl;
//          if (pSlimmedJetTag->at(iJet).pt() > 50)
//            IsoBCSVPassTrigMu_->Fill( pSlimmedJetTag->at(iJet).bDiscriminator("pfCombinedSecondaryVertexV2BJetTags") );
//          if ( pSlimmedJetTag->at(iJet).bDiscriminator("pfCombinedSecondaryVertexV2BJetTags") > .7 )
//            IsoBPtPassTrigMu_->Fill( pSlimmedJetTag->at(iJet).pt() );

          if (!checkDiTau && dR1_PFTau != 1000 && dR2_PFTau != 1000)
          { 
            IsoMatchDiTauDr_->Fill(diTauDr );
            checkDiTau = true;
          }//if !checkDiTau
        }//iftrigRefB.isAvailable
        
      }//for iJet
    }//if
  }//for

///////////////////////////////////////////////////////
// Loop for Non- BJet Trigger Matching
///////////////////////////////////////////////////////
  dR_slimMuTau1 = 1000;
  dR_slimMuTau2 = 1000;
  dR_PFMuTau1 = 1000;
  dR_PFMuTau2 = 1000;
  foundMatchAlready = false;
  bool checkMu = false;
  for (auto iMu = pSlimmedMuonTag->begin(); iMu != pSlimmedMuonTag->end(); ++iMu)
  {
    if (iMu->pt() > 10 && !checkMu)
    {
      checkMu = true;
      for( size_t iJet = 0; iJet < pSlimmedJetTag->size(); ++iJet ) 
      {
        if (pSlimmedJetTag->at(iJet).bDiscriminator("pfCombinedSecondaryVertexV2BJetTags") > .7 && abs(pSlimmedJetTag->at(iJet).eta() ) < 2.6)
          NonIsoBPtPassAllMu_->Fill(pSlimmedJetTag->at(iJet).pt() );
        if (pSlimmedJetTag->at(iJet).pt() > 50 && abs(pSlimmedJetTag->at(iJet).eta() ) < 2.6)
          NonIsoBCSVPassAllMu_->Fill(pSlimmedJetTag->at(iJet).bDiscriminator("pfCombinedSecondaryVertexV2BJetTags"));
    
        const pat::TriggerObjectRef trigRefB( matchHelper.triggerMatchObject( pSlimmedJetTag, iJet, bMatchNonMu_, iEvent, *triggerEvent ) );
        if ( trigRefB.isAvailable() && trigRefB.isNonnull() )
        {
          cout << "<---------Non- BJet Trigger Match----------->" << endl;
          if (!foundMatchAlready)
          {
            HLTPass_->Fill(2);
            foundMatchAlready = true;
          }// if !foundMatchAlready
          NonIsoBPt_->Fill(trigRefB->pt() );
          dR_slimMuTau1 = sqrt( (trigRefB->eta() - slimTau1Eta)*(trigRefB->eta() - slimTau1Eta)  +  (trigRefB->phi() - slimTau1Phi)*(trigRefB->phi() - slimTau1Phi) );
          dR_slimMuTau2 = sqrt( (trigRefB->eta() - slimTau2Eta)*(trigRefB->eta() - slimTau2Eta)  +  (trigRefB->phi() - slimTau2Phi)*(trigRefB->phi() - slimTau2Phi) );
          dR_PFMuTau1 = sqrt( (trigRefB->eta() - PFTau1Eta)*(trigRefB->eta() - PFTau1Eta)  +  (trigRefB->phi() - PFTau1Phi)*(trigRefB->phi() - PFTau1Phi) );
          dR_PFMuTau2 = sqrt( (trigRefB->eta() - PFTau2Eta)*(trigRefB->eta() - PFTau2Eta)  +  (trigRefB->phi() - PFTau2Phi)*(trigRefB->phi() - PFTau2Phi) );
          std::cout << "\t\tdR_slimMuTau1= " << dR_slimMuTau1 << "  dR_slimMuTau2= " << dR_slimMuTau2 << " dR_PFMuTau1= " << dR_PFMuTau1 << "  dR_PFMuTau2= " << dR_PFMuTau2 << std::endl;
          if (pSlimmedJetTag->at(iJet).pt() > 50)
          {
      	    NonIsoBCSVPassTrigMu_->Fill( pSlimmedJetTag->at(iJet).bDiscriminator("pfCombinedSecondaryVertexV2BJetTags") );
    	    double dR=sqrt((trigRefB->eta()-pSlimmedJetTag->at(iJet).eta())*(trigRefB->eta()-pSlimmedJetTag->at(iJet).eta())+(trigRefB->phi()-pSlimmedJetTag->at(iJet).phi())*(trigRefB->phi()-pSlimmedJetTag->at(iJet).phi()));
            std::cout << "\t\tNon pt > 50 dR(trig_b matched_b)= " << dR << std::endl;
    	    TrigMatchdR_->Fill(dR);
          }
          if ( pSlimmedJetTag->at(iJet).bDiscriminator("pfCombinedSecondaryVertexV2BJetTags") > .7 )
          {
    	    NonIsoBPtPassTrigMu_->Fill( pSlimmedJetTag->at(iJet).pt() );
            double dR=sqrt((trigRefB->eta()-pSlimmedJetTag->at(iJet).eta())*(trigRefB->eta()-pSlimmedJetTag->at(iJet).eta())+(trigRefB->phi()-pSlimmedJetTag->at(iJet).phi())*(trigRefB->phi()-pSlimmedJetTag->at(iJet).phi()));
            std::cout << "\t\tNon CSV > .7 dR(trig_b matched_b)= " << dR << std::endl;
    	    TrigMatchdR_->Fill(dR);
          }
          NonIsoTrigPtvsRecoPt_->Fill(trigRefB->pt(), pSlimmedJetTag->at(iJet).pt() );
          NonIsoTrigEtavsRecoEta_->Fill(trigRefB->eta(), pSlimmedJetTag->at(iJet).eta() );
        }//if
      }//for iJet
    }//if iMu->pt() > 10
  }//for iMu

///////////////////////////////////////////////////////
// Loop for  BJet Trigger Matching
///////////////////////////////////////////////////////
  dR_slimMuTau1 = 1000;
  dR_slimMuTau2 = 1000;
  dR_PFMuTau1 = 1000;
  dR_PFMuTau2 = 1000;
  foundMatchAlready = false;
  checkMu = false;
  for (auto iMu = pSlimmedMuonTag->begin(); iMu != pSlimmedMuonTag->end(); ++iMu)
  {
    if (iMu->pt() > 24 && iMu->isLooseMuon() && !checkMu)
    {
      checkMu = true;
      for( size_t iJet = 0; iJet < pSlimmedJetTag->size(); ++iJet ) 
      {
        if (pSlimmedJetTag->at(iJet).bDiscriminator("pfCombinedSecondaryVertexV2BJetTags") > .7 && abs(pSlimmedJetTag->at(iJet).eta() ) < 2.6 )
          IsoBPtPassAllMu_->Fill(pSlimmedJetTag->at(iJet).pt() );
        if (pSlimmedJetTag->at(iJet).pt() > 50 && abs(pSlimmedJetTag->at(iJet).eta() ) < 2.6)
          IsoBCSVPassAllMu_->Fill(pSlimmedJetTag->at(iJet).bDiscriminator("pfCombinedSecondaryVertexV2BJetTags"));
    
        const pat::TriggerObjectRef trigRefB( matchHelper.triggerMatchObject( pSlimmedJetTag, iJet, bMatchMu_, iEvent, *triggerEvent ) );
        if ( trigRefB.isAvailable() && trigRefB.isNonnull() )
        {
          cout << "<--------- BJet Trigger Match----------->" << endl;
          if (!foundMatchAlready)
          {
            HLTPass_->Fill(4);
            foundMatchAlready = true;
          }// if !foundMatchAlready
          if (pSlimmedJetTag->at(iJet).bDiscriminator("pfCombinedSecondaryVertexV2BJetTags") < .54)
    	    CSVlt54Pt_->Fill(pSlimmedJetTag->at(iJet).pt() );
          IsoBPt_->Fill(trigRefB->pt() );
          dR_slimMuTau1 = sqrt( (trigRefB->eta() - slimTau1Eta)*(trigRefB->eta() - slimTau1Eta)  +  (trigRefB->phi() - slimTau1Phi)*(trigRefB->phi() - slimTau1Phi) );
          dR_slimMuTau2 = sqrt( (trigRefB->eta() - slimTau2Eta)*(trigRefB->eta() - slimTau2Eta)  +  (trigRefB->phi() - slimTau2Phi)*(trigRefB->phi() - slimTau2Phi) );
          dR_PFMuTau1 = sqrt( (trigRefB->eta() - PFTau1Eta)*(trigRefB->eta() - PFTau1Eta)  +  (trigRefB->phi() - PFTau1Phi)*(trigRefB->phi() - PFTau1Phi) );
          dR_PFMuTau2 = sqrt( (trigRefB->eta() - PFTau2Eta)*(trigRefB->eta() - PFTau2Eta)  +  (trigRefB->phi() - PFTau2Phi)*(trigRefB->phi() - PFTau2Phi) );
          std::cout << "\t\tdR_slimMuTau1= " << dR_slimMuTau1 << "  dR_slimMuTau2= " << dR_slimMuTau2 << " dR_PFMuTau1= " << dR_PFMuTau1 << "  dR_PFMuTau2= " << dR_PFMuTau2 << std::endl;
          if (pSlimmedJetTag->at(iJet).pt() > 50)
          {
            IsoBCSVPassTrigMu_->Fill(pSlimmedJetTag->at(iJet).bDiscriminator("pfCombinedSecondaryVertexV2BJetTags") );
            double dR=sqrt((trigRefB->eta()-pSlimmedJetTag->at(iJet).eta())*(trigRefB->eta()-pSlimmedJetTag->at(iJet).eta())+(trigRefB->phi()-pSlimmedJetTag->at(iJet).phi())*(trigRefB->phi()-pSlimmedJetTag->at(iJet).phi()));
            std::cout << "\t\t pt > 50 dR(trig_b matched_b)= " << dR << std::endl;
    	    TrigMatchdR_->Fill(dR);
          }
          if (pSlimmedJetTag->at(iJet).bDiscriminator("pfCombinedSecondaryVertexV2BJetTags") > .7 )
          {
            IsoBPtPassTrigMu_->Fill(pSlimmedJetTag->at(iJet).pt() );
            double dR=sqrt((trigRefB->eta()-pSlimmedJetTag->at(iJet).eta())*(trigRefB->eta()-pSlimmedJetTag->at(iJet).eta())+(trigRefB->phi()-pSlimmedJetTag->at(iJet).phi())*(trigRefB->phi()-pSlimmedJetTag->at(iJet).phi()));
            std::cout << "\t\t CSV > .7 dR(trig_b matched_b)= " << dR << std::endl;
    	    TrigMatchdR_->Fill(dR);
          }
          IsoTrigPtvsRecoPt_->Fill(trigRefB->pt(), pSlimmedJetTag->at(iJet).pt() );
        }//if
      }//for
    }//if iMu
  }//for iMu
*/
  double dPhi = reco::deltaPhi(tau1Ref->phi(), tau2Ref->phi() );
  double genDiTauDR = sqrt( (tau1Ref->eta() - tau2Ref->eta())*(tau1Ref->eta() - tau2Ref->eta())  +  (dPhi)*(dPhi) );
  double bestBPt = 0;
  for (auto iJet = pSlimmedJetTag->begin(); iJet < pSlimmedJetTag->end(); ++iJet )
  {
    if (iJet->bDiscriminator("pfCombinedSecondaryVertexV2BJetTags") > .54 && bestBPt < iJet->pt() )
      bestBPt = iJet->pt();
  }//for

  if (bestBPt > 40) 
    DiTaudR_Ptgt40_->Fill(genDiTauDR );
  if (bestBPt > 60)
    DiTaudR_Ptgt60_->Fill(genDiTauDR );
  if (bestBPt > 80)
    DiTaudR_Ptgt80_->Fill(genDiTauDR );
  if (bestBPt > 100)
    DiTaudR_Ptgt100_->Fill(genDiTauDR );


}//End BBAAnalyzer::analyze


// ------------ method called once each job just before starting event loop  ------------
void BBAAnalyzer::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  //Open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //Book histograms
  HLTPass_ = new TH1F("HLTPass", "", 7, -.5, 6.5);
      HLTPass_->GetXaxis()->SetBinLabel(1, "TotalEvents"); 
      HLTPass_->GetXaxis()->SetBinLabel(2, "Muon NonMu Match");
      HLTPass_->GetXaxis()->SetBinLabel(3, "bJet NonMu Match");
      HLTPass_->GetXaxis()->SetBinLabel(4, "Muon Mu Match");
      HLTPass_->GetXaxis()->SetBinLabel(5, "bJet Mu Match");
      HLTPass_->GetXaxis()->SetBinLabel(6, "Muon bTagMu Match");
      HLTPass_->GetXaxis()->SetBinLabel(7, "bJet bTagMu Match");
  NonIsoMuPt_ 	= new TH1F("NonIsoMuPt", "", 50, 0, 100);
  IsoMuPt_ 	= new TH1F("IsoMuPt",    "", 50, 0, 100);
  NonIsoBPt_ 	= new TH1F("NonIsoBPt",  "", 50, 0, 150);
  IsoBPt_ 	= new TH1F("IsoBPt",     "", 50, 0, 150);
  DiTaudR_Ptgt40_       = new TH1F("DiTaudR_Ptgt40",     "", 50, 0, 10);
  DiTaudR_Ptgt60_       = new TH1F("DiTaudR_Ptgt60",     "", 50, 0, 10);
  DiTaudR_Ptgt80_       = new TH1F("DiTaudR_Ptgt80",     "", 50, 0, 10);
  DiTaudR_Ptgt100_       = new TH1F("DiTaudR_Ptgt100",     "", 50, 0, 10);

  NMu_      = new TH1F("NMu", "", 12, -.5, 11.5);
  NSlimTauvsNhpsTau_  = new TH2F("NSlimTauvsNhpsTau", "", 20, -.5, 19.5, 12, -.5, 11.5);
  NonIsoTrigPtvsRecoPt_  = new TH2F("NonIsoTrigPtvsRecoPt", "", 30, 0, 150, 30, 0, 150);
  IsoTrigPtvsRecoPt_  = new TH2F("IsoTrigPtvsRecoPt", "", 30, 0, 150, 30, 0, 150);
  NonIsoTrigEtavsRecoEta_  = new TH2F("NonIsoTrigEtavsRecoEta", "", 30, 0, 3, 30, 0, 3);

  NonIsoMatchDiTauDr_ = new TH1F("NonIsoMatchDiTauDr", "", 50, 0, 8);
  IsoMatchDiTauDr_    = new TH1F("IsoMatchDiTauDr",    "", 50, 0, 8);

  NonIsoBCSVTurnOn_ = new TH1F("NonIsoBCSVTurnOn", "", 33, 0, 1);
  IsoBCSVTurnOn_    = new TH1F("IsoBCSVTurnOn",    "", 33, 0, 1);
  NonIsoBPtTurnOn_  = new TH1F("NonIsoBPtTurnOn",  "", 33, 0, 132);
  IsoBPtTurnOn_     = new TH1F("IsoBPtTurnOn",     "", 33, 0, 132);

  NonIsoBCSVPassAll_ = new TH1F("NonIsoBCSVPassAll", "", 33, 0, 1);
  IsoBCSVPassAll_    = new TH1F("IsoBCSVPassAll",    "", 33, 0, 1);
  NonIsoBPtPassAll_  = new TH1F("NonIsoBPtPassAll",  "", 33, 0, 132);
  IsoBPtPassAll_     = new TH1F("IsoBPtPassAll",     "", 33, 0, 132);

  NonIsoBCSVPassTrig_ = new TH1F("NonIsoBCSVPassTrig", "", 33, 0, 1);
  IsoBCSVPassTrig_    = new TH1F("IsoBCSVPassTrig",    "", 33, 0, 1);
  NonIsoBPtPassTrig_  = new TH1F("NonIsoBPtPassTrig",  "", 33, 0, 132);
  IsoBPtPassTrig_     = new TH1F("IsoBPtPassTrig",     "", 33, 0, 132);

  NonIsoBCSVTurnOnMu_ = new TH1F("NonIsoBCSVTurnOnMu", "", 33, 0, 1);
  IsoBCSVTurnOnMu_    = new TH1F("IsoBCSVTurnOnMu",    "", 33, 0, 1);
  NonIsoBPtTurnOnMu_  = new TH1F("NonIsoBPtTurnOnMu",  "", 33, 0, 132);
  IsoBPtTurnOnMu_     = new TH1F("IsoBPtTurnOnMu",     "", 33, 0, 132);

  NonIsoBCSVPassAllMu_ = new TH1F("NonIsoBCSVPassAllMu", "", 33, 0, 1);
  IsoBCSVPassAllMu_    = new TH1F("IsoBCSVPassAllMu",    "", 33, 0, 1);
  NonIsoBPtPassAllMu_  = new TH1F("NonIsoBPtPassAllMu",  "", 33, 0, 132);
  IsoBPtPassAllMu_     = new TH1F("IsoBPtPassAllMu",     "", 33, 0, 132);

  NonIsoBCSVPassTrigMu_ = new TH1F("NonIsoBCSVPassTrigMu", "", 33, 0, 1);
  IsoBCSVPassTrigMu_    = new TH1F("IsoBCSVPassTrigMu",    "", 33, 0, 1);
  NonIsoBPtPassTrigMu_  = new TH1F("NonIsoBPtPassTrigMu",  "", 33, 0, 132);
  IsoBPtPassTrigMu_     = new TH1F("IsoBPtPassTrigMu",     "", 33, 0, 132);

  TrigMatchdR_     = new TH1F("TrigMatchdR",     "", 33, 0, 2);
  CSVlt54Pt_     = new TH1F("CSVlt54Pt",     "", 33, 0, 120);

}

// ------------ method called once each job just after ending the event loop  ------------
void BBAAnalyzer::endJob()
{



  //Make the Canvases
  TCanvas HLTPassCanvas_("HLTPassCanvas","",600,600);
  TCanvas NonIsoMuPtCanvas_("NonIsoMuPtCanvas","",600,600);
  TCanvas IsoMuPtCanvas_("IsoMuPtCanvas","",600,600);
  TCanvas NonIsoBPtCanvas_("NonIsoBPtCanvas","",600,600);
  TCanvas IsoBPtCanvas_("IsoBPtCanvas","",600,600);
  TCanvas DiTaudR_Ptgt40Canvas_("DiTaudR_Ptgt40Canvas","",600,600);
  TCanvas DiTaudR_Ptgt60Canvas_("DiTaudR_Ptgt60Canvas","",600,600);
  TCanvas DiTaudR_Ptgt80Canvas_("DiTaudR_Ptgt80Canvas","",600,600);
  TCanvas DiTaudR_Ptgt100Canvas_("DiTaudR_Ptgt100Canvas","",600,600);

  TCanvas NMuCanvas_("NMuCanvas","",600,600);
  TCanvas NSlimTauvsNhpsTauCanvas("NSlimTauvsNhpsTau","",600,600);
  TCanvas NonIsoTrigPtvsRecoPtCanvas("NonIsoTrigPtvsRecoPt","",600,600);
  TCanvas IsoTrigPtvsRecoPtCanvas("IsoTrigPtvsRecoPt","",600,600);
  TCanvas NonIsoTrigEtavsRecoEtaCanvas("NonIsoTrigEtavsRecoEta","",600,600);

  TCanvas NonIsoMatchDiTauDrCanvas_("NonIsoMatchDiTauDrCanvas","",600,600);
  TCanvas IsoMatchDiTauDrCanvas_("IsoMatchDiTauDrCanvas","",600,600);   

  TCanvas NonIsoBCSVTurnOnCanvas_("NonIsoBCSVTurnOnCanvas","",600,600); 
  TCanvas IsoBCSVTurnOnCanvas_("IsoBCSVTurnOnCanvas","",600,600);    
  TCanvas NonIsoBPtTurnOnCanvas_("NonIsoBPtTurnOnCanvas","",600,600);  
  TCanvas IsoBPtTurnOnCanvas_("IsoBPtTurnOnCanvas","",600,600);    

  TCanvas NonIsoBCSVPassAllCanvas_("NonIsoBCSVPassAllCanvas","",600,600);
  TCanvas IsoBCSVPassAllCanvas_("IsoBCSVPassAllCanvas","",600,600);  
  TCanvas NonIsoBPtPassAllCanvas_("NonIsoBPtPassAllCanvas","",600,600);
  TCanvas IsoBPtPassAllCanvas_("IsoBPtPassAllCanvas","",600,600);  

  TCanvas NonIsoBCSVPassTrigCanvas_("NonIsoBCSVPassTrigCanvas","",600,600);
  TCanvas IsoBCSVPassTrigCanvas_("IsoBCSVPassTrigCanvas","",600,600);  
  TCanvas NonIsoBPtPassTrigCanvas_("NonIsoBPtPassTrigCanvas","",600,600);
  TCanvas IsoBPtPassTrigCanvas_("IsoBPtPassTrigCanvas","",600,600);  

  NonIsoBCSVTurnOn_->Divide(NonIsoBCSVPassTrig_, NonIsoBCSVPassAll_);
  IsoBCSVTurnOn_->Divide(IsoBCSVPassTrig_, IsoBCSVPassAll_);
  NonIsoBPtTurnOn_->Divide(NonIsoBPtPassTrig_, NonIsoBPtPassAll_);
  IsoBPtTurnOn_->Divide(IsoBPtPassTrig_, IsoBPtPassAll_);

  TCanvas NonIsoBCSVTurnOnMuCanvas_("NonIsoBCSVTurnOnMuCanvas","",600,600);
  TCanvas IsoBCSVTurnOnMuCanvas_("IsoBCSVTurnOnMuCanvas","",600,600);
  TCanvas NonIsoBPtTurnOnMuCanvas_("NonIsoBPtTurnOnMuCanvas","",600,600);
  TCanvas IsoBPtTurnOnMuCanvas_("IsoBPtTurnOnMuCanvas","",600,600);

  TCanvas NonIsoBCSVPassAllMuCanvas_("NonIsoBCSVPassAllMuCanvas","",600,600);
  TCanvas IsoBCSVPassAllMuCanvas_("IsoBCSVPassAllMuCanvas","",600,600);
  TCanvas NonIsoBPtPassAllMuCanvas_("NonIsoBPtPassAllMuCanvas","",600,600);
  TCanvas IsoBPtPassAllMuCanvas_("IsoBPtPassAllMuCanvas","",600,600);

  TCanvas NonIsoBCSVPassTrigMuCanvas_("NonIsoBCSVPassTrigMuCanvas","",600,600);
  TCanvas IsoBCSVPassTrigMuCanvas_("IsoBCSVPassTrigMuCanvas","",600,600);
  TCanvas NonIsoBPtPassTrigMuCanvas_("NonIsoBPtPassTrigMuCanvas","",600,600);
  TCanvas IsoBPtPassTrigMuCanvas_("IsoBPtPassTrigMuCanvas","",600,600);

  TCanvas TrigMatchdRCanvas_("TrigMatchdRCanvas","",600,600);
  TCanvas CSVlt54PtCanvas_("CSVlt54PtCanvas","",600,600);

  NonIsoBCSVTurnOnMu_->Divide(NonIsoBCSVPassTrigMu_, NonIsoBCSVPassAllMu_);
  IsoBCSVTurnOnMu_->Divide(IsoBCSVPassTrigMu_, IsoBCSVPassAllMu_);
  NonIsoBPtTurnOnMu_->Divide(NonIsoBPtPassTrigMu_, NonIsoBPtPassAllMu_);
  IsoBPtTurnOnMu_->Divide(IsoBPtPassTrigMu_, IsoBPtPassAllMu_);

 
std::cout << "<----------------Declared Canvases-------------->" << std::endl;

  //Format the 1D plots and Draw (canvas, hist, grid, log y, log z, color, size, style, xAxisTitle, xTitleSize, xLabelSize, xTitleOffSet, yAxisTitle, yTitleSize, yLabelSize, yTitleOffset)
  VariousFunctions::formatAndDrawCanvasAndHist1D(HLTPassCanvas_, HLTPass_, 1, 0, 0, kBlack, 7, 20, "", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NonIsoMuPtCanvas_, NonIsoMuPt_, 1, 0, 0, kBlack, 7, 20, "pt(#mu) of NonIsoTrig", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoMuPtCanvas_, IsoMuPt_, 1, 0, 0, kBlack, 7, 20, "pt(#mu) of IsoTrig", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NonIsoBPtCanvas_, NonIsoBPt_, 1, 0, 0, kBlack, 7, 20, "pt(b) of NonIsoTrig", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoBPtCanvas_, IsoBPt_, 1, 0, 0, kBlack, 7, 20, "pt(b) of IsoTrig", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiTaudR_Ptgt40Canvas_, DiTaudR_Ptgt40_, 1, 0, 0, kBlack, 7, 20, "dR(#tau#tau) pt(b) > 40", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiTaudR_Ptgt60Canvas_, DiTaudR_Ptgt60_, 1, 0, 0, kBlack, 7, 20, "dR(#tau#tau) pt(b) > 60", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiTaudR_Ptgt80Canvas_, DiTaudR_Ptgt80_, 1, 0, 0, kBlack, 7, 20, "dR(#tau#tau) pt(b) > 80", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiTaudR_Ptgt100Canvas_, DiTaudR_Ptgt100_, 1, 0, 0, kBlack, 7, 20, "dR(#tau#tau) pt(b) > 100", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(NMuCanvas_, NMu_, 1, 0, 0, kBlack, 7, 20, "Number of slimmed muons", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist2D(NSlimTauvsNhpsTauCanvas, NSlimTauvsNhpsTau_, 0, 0, 0, kBlack, 7, 20, "Trig Pt(b)", .04, .04, 1.1, "Reco Pt(b)", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(NonIsoTrigPtvsRecoPtCanvas, NonIsoTrigPtvsRecoPt_, 0, 0, 0, kBlack, 7, 20, "NonIso Trig Pt", .04, .04, 1.1, "RECO Trig Match Pt", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(IsoTrigPtvsRecoPtCanvas, IsoTrigPtvsRecoPt_, 0, 0, 0, kBlack, 7, 20, "Iso Trig Pt", .04, .04, 1.1, "RECO Match Trig Pt", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(NonIsoTrigEtavsRecoEtaCanvas, NonIsoTrigEtavsRecoEta_, 0, 0, 0, kBlack, 7, 20, "NonIso Trig Eta", .04, .04, 1.1, "RECO Trig Match Eta", .04, .04, 1.6, "", .04, .04, 1.0);

  VariousFunctions::formatAndDrawCanvasAndHist1D(NonIsoMatchDiTauDrCanvas_, NonIsoMatchDiTauDr_, 1, 0, 0, kBlack, 7, 20, "Non-Iso Trig #DeltaR(#tau#tau)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoMatchDiTauDrCanvas_, IsoMatchDiTauDr_, 1, 0, 0, kBlack, 7, 20, "Iso Trig #DeltaR(#tau#tau)", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(NonIsoBCSVTurnOnCanvas_, NonIsoBCSVTurnOn_, 1, 0, 0, kBlack, 7, 20, "Non-Iso Trig CSV(bJet) TurnOn pt > 50", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoBCSVTurnOnCanvas_, IsoBCSVTurnOn_, 1, 0, 0, kBlack, 7, 20, "Iso Trig CSV(bJet) TurnOn pt > 50", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NonIsoBPtTurnOnCanvas_, NonIsoBPtTurnOn_, 1, 0, 0, kBlack, 7, 20, "Non-Iso Trig Pt(bJet) TurnOn CSV > .7", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoBPtTurnOnCanvas_, IsoBPtTurnOn_, 1, 0, 0, kBlack, 7, 20, "Iso Trig Pt(bJet) TurnOn CSV > .7", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(NonIsoBCSVPassAllCanvas_, NonIsoBCSVPassAll_, 1, 0, 0, kBlack, 7, 20, "Non-Iso Trig CSV(bJet) PassAll pt > 50", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoBCSVPassAllCanvas_, IsoBCSVPassAll_, 1, 0, 0, kBlack, 7, 20, "Iso Trig CSV(bJet) PassAll pt > 50", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NonIsoBPtPassAllCanvas_, NonIsoBPtPassAll_, 1, 0, 0, kBlack, 7, 20, "Non-Iso Trig Pt(bJet) PassAll CSV > .7", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoBPtPassAllCanvas_, IsoBPtPassAll_, 1, 0, 0, kBlack, 7, 20, "Iso Trig Pt(bJet) PassAll CSV > .7", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(NonIsoBCSVPassTrigCanvas_, NonIsoBCSVPassTrig_, 1, 0, 0, kBlack, 7, 20, "Non-Iso Trig CSV(bJet) PassTrig pt > 50", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoBCSVPassTrigCanvas_, IsoBCSVPassTrig_, 1, 0, 0, kBlack, 7, 20, "Iso Trig CSV(bJet) PassTrig pt > 50", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NonIsoBPtPassTrigCanvas_, NonIsoBPtPassTrig_, 1, 0, 0, kBlack, 7, 20, "Non-Iso Trig Pt(bJet) PassTrig CSV > .7", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoBPtPassTrigCanvas_, IsoBPtPassTrig_, 1, 0, 0, kBlack, 7, 20, "Iso Trig Pt(bJet) PassTrig CSV > .7", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(NonIsoBCSVTurnOnMuCanvas_, NonIsoBCSVTurnOnMu_, 1, 0, 0, kBlack, 7, 20, "Non-Iso Trig CSV(bJet) TurnOn pt > 50", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoBCSVTurnOnMuCanvas_, IsoBCSVTurnOnMu_, 1, 0, 0, kBlack, 7, 20, "Iso Trig CSV(bJet) TurnOn pt > 50", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NonIsoBPtTurnOnMuCanvas_, NonIsoBPtTurnOnMu_, 1, 0, 0, kBlack, 7, 20, "Non-Iso Trig Pt(bJet) TurnOn CSV > .7", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoBPtTurnOnMuCanvas_, IsoBPtTurnOnMu_, 1, 0, 0, kBlack, 7, 20, "Iso Trig Pt(bJet) TurnOn CSV > .7", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(NonIsoBCSVPassAllMuCanvas_, NonIsoBCSVPassAllMu_, 1, 0, 0, kBlack, 7, 20, "Non-Iso Trig CSV(bJet) PassAll pt > 50", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoBCSVPassAllMuCanvas_, IsoBCSVPassAllMu_, 1, 0, 0, kBlack, 7, 20, "Iso Trig CSV(bJet) PassAll pt > 50", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NonIsoBPtPassAllMuCanvas_, NonIsoBPtPassAllMu_, 1, 0, 0, kBlack, 7, 20, "Non-Iso Trig Pt(bJet) PassAll CSV > .7", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoBPtPassAllMuCanvas_, IsoBPtPassAllMu_, 1, 0, 0, kBlack, 7, 20, "Iso Trig Pt(bJet) PassAll CSV > .7", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(NonIsoBCSVPassTrigMuCanvas_, NonIsoBCSVPassTrigMu_, 1, 0, 0, kBlack, 7, 20, "Non-Iso Trig CSV(bJet) PassTrig pt > 50", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoBCSVPassTrigMuCanvas_, IsoBCSVPassTrigMu_, 1, 0, 0, kBlack, 7, 20, "Iso Trig CSV(bJet) PassTrig pt > 50", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NonIsoBPtPassTrigMuCanvas_, NonIsoBPtPassTrigMu_, 1, 0, 0, kBlack, 7, 20, "Non-Iso Trig Pt(bJet) PassTrig CSV > .7", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoBPtPassTrigMuCanvas_, IsoBPtPassTrigMu_, 1, 0, 0, kBlack, 7, 20, "Iso Trig Pt(bJet) PassTrig CSV > .7", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(TrigMatchdRCanvas_, TrigMatchdR_, 1, 0, 0, kBlack, 7, 20, "#DeltaR(trig_b, matched_b)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(CSVlt54PtCanvas_, CSVlt54Pt_, 1, 0, 0, kBlack, 7, 20, "Trig Match CSV < .54 Pt(bJet)", .04, .04, 1.1,  "", .04, .04, 1.0, false);

std::cout << "<----------------Formatted Canvases and Histos-------------->" << std::endl;

  //Write output file
  out_->cd();

  HLTPassCanvas_.Write();
  NonIsoMuPtCanvas_.Write();
  IsoMuPtCanvas_.Write();
  NonIsoBPtCanvas_.Write();
  IsoBPtCanvas_.Write();
  DiTaudR_Ptgt40Canvas_.Write();
  DiTaudR_Ptgt60Canvas_.Write();
  DiTaudR_Ptgt80Canvas_.Write();
  DiTaudR_Ptgt100Canvas_.Write();

  NMuCanvas_.Write();
  NSlimTauvsNhpsTauCanvas.Write();
  NonIsoTrigPtvsRecoPtCanvas.Write();
  IsoTrigPtvsRecoPtCanvas.Write();
  NonIsoTrigEtavsRecoEtaCanvas.Write();

  NonIsoMatchDiTauDrCanvas_.Write();
  IsoMatchDiTauDrCanvas_.Write();
  
  NonIsoBCSVTurnOnCanvas_.Write();
  IsoBCSVTurnOnCanvas_.Write();
  NonIsoBPtTurnOnCanvas_.Write();
  IsoBPtTurnOnCanvas_.Write();

  NonIsoBCSVPassAllCanvas_.Write();
  IsoBCSVPassAllCanvas_.Write();
  NonIsoBPtPassAllCanvas_.Write();
  IsoBPtPassAllCanvas_.Write();
  
  NonIsoBCSVPassTrigCanvas_.Write();
  IsoBCSVPassTrigCanvas_.Write();
  NonIsoBPtPassTrigCanvas_.Write();
  IsoBPtPassTrigCanvas_.Write();
  
  NonIsoBCSVTurnOnMuCanvas_.Write();
  IsoBCSVTurnOnMuCanvas_.Write();
  NonIsoBPtTurnOnMuCanvas_.Write();
  IsoBPtTurnOnMuCanvas_.Write();

  NonIsoBCSVPassAllMuCanvas_.Write();
  IsoBCSVPassAllMuCanvas_.Write();
  NonIsoBPtPassAllMuCanvas_.Write();
  IsoBPtPassAllMuCanvas_.Write();

  NonIsoBCSVPassTrigMuCanvas_.Write();
  IsoBCSVPassTrigMuCanvas_.Write();
  NonIsoBPtPassTrigMuCanvas_.Write();
  IsoBPtPassTrigMuCanvas_.Write();

  TrigMatchdRCanvas_.Write();
  CSVlt54PtCanvas_.Write();

  out_->Write();
  out_->Close();
std::cout << "DONE" << std::endl;
}//EndJob

// ------------ method called when starting to processes a run  ------------
void BBAAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void BBAAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void BBAAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void BBAAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

//Delete Memory
void BBAAnalyzer::reset(const bool doDelete)
{
  if ((doDelete) && (HLTPass_ != NULL)) delete HLTPass_;
  HLTPass_ = NULL;
  if ((doDelete) && (NonIsoMuPt_ != NULL)) delete NonIsoMuPt_;
  NonIsoMuPt_ = NULL;
  if ((doDelete) && (IsoMuPt_ != NULL)) delete IsoMuPt_;
  IsoMuPt_ = NULL;
  if ((doDelete) && (NonIsoBPt_ != NULL)) delete NonIsoBPt_;
  NonIsoBPt_ = NULL;
  if ((doDelete) && (IsoBPt_ != NULL)) delete IsoBPt_;
  IsoBPt_ = NULL;
  if ((doDelete) && (DiTaudR_Ptgt40_ != NULL)) delete DiTaudR_Ptgt40_;
  DiTaudR_Ptgt40_ = NULL;
  if ((doDelete) && (DiTaudR_Ptgt60_ != NULL)) delete DiTaudR_Ptgt60_;
  DiTaudR_Ptgt60_ = NULL;
  if ((doDelete) && (DiTaudR_Ptgt80_ != NULL)) delete DiTaudR_Ptgt80_;
  DiTaudR_Ptgt80_ = NULL;
  if ((doDelete) && (DiTaudR_Ptgt100_ != NULL)) delete DiTaudR_Ptgt100_;
  DiTaudR_Ptgt100_ = NULL;

  if ((doDelete) && (NMu_ != NULL)) delete NMu_;
  NMu_ = NULL;
  if ((doDelete) && (NSlimTauvsNhpsTau_ != NULL)) delete NSlimTauvsNhpsTau_;
  NSlimTauvsNhpsTau_ = NULL;
  if ((doDelete) && (NonIsoTrigPtvsRecoPt_ != NULL)) delete NonIsoTrigPtvsRecoPt_;
  NonIsoTrigPtvsRecoPt_ = NULL;
  if ((doDelete) && (IsoTrigPtvsRecoPt_ != NULL)) delete IsoTrigPtvsRecoPt_;
  IsoTrigPtvsRecoPt_ = NULL;
  if ((doDelete) && (NonIsoTrigEtavsRecoEta_ != NULL)) delete NonIsoTrigEtavsRecoEta_;
  NonIsoTrigEtavsRecoEta_ = NULL;

  if ((doDelete) && (NonIsoMatchDiTauDr_ != NULL)) delete NonIsoMatchDiTauDr_;
  NonIsoMatchDiTauDr_ = NULL;
  if ((doDelete) && (IsoMatchDiTauDr_ != NULL)) delete IsoMatchDiTauDr_;
  IsoMatchDiTauDr_ = NULL;

  if ((doDelete) && (NonIsoBCSVTurnOn_ != NULL)) delete NonIsoBCSVTurnOn_;
  NonIsoBCSVTurnOn_ = NULL;
  if ((doDelete) && (IsoBCSVTurnOn_ != NULL)) delete IsoBCSVTurnOn_;
  IsoBCSVTurnOn_ = NULL;
  if ((doDelete) && (NonIsoBPtTurnOn_ != NULL)) delete NonIsoBPtTurnOn_;
  NonIsoBPtTurnOn_ = NULL;
  if ((doDelete) && (IsoBPtTurnOn_ != NULL)) delete IsoBPtTurnOn_;
  IsoBPtTurnOn_ = NULL;
  
  if ((doDelete) && (NonIsoBCSVPassAll_ != NULL)) delete NonIsoBCSVPassAll_;
  NonIsoBCSVPassAll_ = NULL;
  if ((doDelete) && (IsoBCSVPassAll_ != NULL)) delete IsoBCSVPassAll_;
  IsoBCSVPassAll_ = NULL;
  if ((doDelete) && (NonIsoBPtPassAll_ != NULL)) delete NonIsoBPtPassAll_;
  NonIsoBPtPassAll_ = NULL;
  if ((doDelete) && (IsoBPtPassAll_ != NULL)) delete IsoBPtPassAll_;
  IsoBPtPassAll_ = NULL;
  
  if ((doDelete) && (NonIsoBCSVPassTrig_ != NULL)) delete NonIsoBCSVPassTrig_;
  NonIsoBCSVPassTrig_ = NULL;
  if ((doDelete) && (IsoBCSVPassTrig_ != NULL)) delete IsoBCSVPassTrig_;
  IsoBCSVPassTrig_ = NULL;
  if ((doDelete) && (NonIsoBPtPassTrig_ != NULL)) delete NonIsoBPtPassTrig_;
  NonIsoBPtPassTrig_ = NULL;
  if ((doDelete) && (IsoBPtPassTrig_ != NULL)) delete IsoBPtPassTrig_;
  IsoBPtPassTrig_ = NULL;

  if ((doDelete) && (NonIsoBCSVTurnOnMu_ != NULL)) delete NonIsoBCSVTurnOnMu_;
  NonIsoBCSVTurnOnMu_ = NULL;
  if ((doDelete) && (IsoBCSVTurnOnMu_ != NULL)) delete IsoBCSVTurnOnMu_;
  IsoBCSVTurnOnMu_ = NULL;
  if ((doDelete) && (NonIsoBPtTurnOnMu_ != NULL)) delete NonIsoBPtTurnOnMu_;
  NonIsoBPtTurnOnMu_ = NULL;
  if ((doDelete) && (IsoBPtTurnOnMu_ != NULL)) delete IsoBPtTurnOnMu_;
  IsoBPtTurnOnMu_ = NULL;

  if ((doDelete) && (NonIsoBCSVPassAllMu_ != NULL)) delete NonIsoBCSVPassAllMu_;
  NonIsoBCSVPassAllMu_ = NULL;
  if ((doDelete) && (IsoBCSVPassAllMu_ != NULL)) delete IsoBCSVPassAllMu_;
  IsoBCSVPassAllMu_ = NULL;
  if ((doDelete) && (NonIsoBPtPassAllMu_ != NULL)) delete NonIsoBPtPassAllMu_;
  NonIsoBPtPassAllMu_ = NULL;
  if ((doDelete) && (IsoBPtPassAllMu_ != NULL)) delete IsoBPtPassAllMu_;
  IsoBPtPassAllMu_ = NULL;

  if ((doDelete) && (NonIsoBCSVPassTrigMu_ != NULL)) delete NonIsoBCSVPassTrigMu_;
  NonIsoBCSVPassTrigMu_ = NULL;
  if ((doDelete) && (IsoBCSVPassTrigMu_ != NULL)) delete IsoBCSVPassTrigMu_;
  IsoBCSVPassTrigMu_ = NULL;
  if ((doDelete) && (NonIsoBPtPassTrigMu_ != NULL)) delete NonIsoBPtPassTrigMu_;
  NonIsoBPtPassTrigMu_ = NULL;
  if ((doDelete) && (IsoBPtPassTrigMu_ != NULL)) delete IsoBPtPassTrigMu_;
  IsoBPtPassTrigMu_ = NULL;

  if ((doDelete) && (TrigMatchdR_ != NULL)) delete TrigMatchdR_;
  TrigMatchdR_ = NULL;
  if ((doDelete) && (CSVlt54Pt_ != NULL)) delete CSVlt54Pt_;
  CSVlt54Pt_ = NULL;

}//void BBAAnalyzer

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void BBAAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BBAAnalyzer);

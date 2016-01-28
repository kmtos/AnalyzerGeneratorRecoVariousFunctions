// -*- C++ -*-
//
// Package:    GGHAnalyzer
// Class:      GGHAnalyzer
// 
/**\class GGHAnalyzer GGHAnalyzer.cc Analyzer/src/GGHAnalyzer.cc

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

#include "BBA/VariousFunctions/interface/VariousFunctions.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;


//
// class declaration
//

class GGHAnalyzer : public edm::EDAnalyzer {
   public:
      typedef reco::JetFloatAssociation::Container JetBCEnergyRatioCollection;
      explicit GGHAnalyzer(const edm::ParameterSet&);
      ~GGHAnalyzer();

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
      edm::EDGetTokenT<edm::ValueMap<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > > muonMapTag_;
      edm::EDGetTokenT<edm::ValueMap<reco::PFJetRef> >  jetValMapTag_;
      edm::EDGetTokenT<vector<reco::PFTau> > tauRECOTag_;
      edm::EDGetTokenT<vector<reco::PFTau> > tauCJTag_;
      edm::EDGetTokenT<vector<reco::RecoTauPiZero> > pizerosTag_;      
      edm::EDGetTokenT<reco::PFTauDiscriminator>  looseIsoTagCJ_;
      edm::EDGetTokenT<reco::PFTauDiscriminator>  medIsoTagCJ_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> tightIsoTagCJ_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> decayModeFindingTagCJ_;
      edm::EDGetTokenT<reco::PFTauDiscriminator>  looseIsoTagRECO_;
      edm::EDGetTokenT<reco::PFTauDiscriminator>  medIsoTagRECO_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> tightIsoTagRECO_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> decayModeFindingTagRECO_;


      //Histograms
      TH1F* NEvents_;   
      TH1F* NMuRemoved_;
      TH1F* TauMuTauHaddR_;
      TH1F* NConstituentsCJ_;
      TH2F* NTausRECOvsCLEANJETS_;
      TH2F* GenDiTaudRvsCJDiTaudR_;
      TH2F* TauHadnConstvsPt_;

      TH1F* MatchedLooseIsoRECOPt_;
      TH1F* MatchedMedIsoRECOPt_;
      TH1F* MatchedTightIsoRECOPt_;
      TH1F* MatchedDMFindRECOPt_;
      TH1F* MatchedRECOPt_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOPt_;
      TGraphAsymmErrors* FinalEffMedIsoRECOPt_;
      TGraphAsymmErrors* FinalEffTightIsoRECOPt_;
      TGraphAsymmErrors* FinalEffDMFindRECOPt_;

      TH1F* MatchedLooseIsoCJPt_;
      TH1F* MatchedMedIsoCJPt_;
      TH1F* MatchedTightIsoCJPt_;
      TH1F* MatchedDMFindCJPt_;
      TH1F* MatchedCJPt_;
      TGraphAsymmErrors* FinalEffLooseIsoCJPt_;
      TGraphAsymmErrors* FinalEffMedIsoCJPt_;
      TGraphAsymmErrors* FinalEffTightIsoCJPt_;
      TGraphAsymmErrors* FinalEffDMFindCJPt_;

      TH1F* MatchedLooseIsoRECOdR_;
      TH1F* MatchedMedIsoRECOdR_;
      TH1F* MatchedTightIsoRECOdR_;
      TH1F* MatchedDMFindRECOdR_;
      TH1F* MatchedRECOdR_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOdR_;
      TGraphAsymmErrors* FinalEffMedIsoRECOdR_;
      TGraphAsymmErrors* FinalEffTightIsoRECOdR_;
      TGraphAsymmErrors* FinalEffDMFindRECOdR_;

      TH1F* MatchedLooseIsoCJdR_;
      TH1F* MatchedMedIsoCJdR_;
      TH1F* MatchedTightIsoCJdR_;
      TH1F* MatchedDMFindCJdR_;
      TH1F* MatchedCJdR_;
      TGraphAsymmErrors* FinalEffLooseIsoCJdR_;
      TGraphAsymmErrors* FinalEffMedIsoCJdR_;
      TGraphAsymmErrors* FinalEffTightIsoCJdR_;
      TGraphAsymmErrors* FinalEffDMFindCJdR_;

      TH1F* MatchedLooseIsoRECOPtGen_;
      TH1F* MatchedMedIsoRECOPtGen_;
      TH1F* MatchedTightIsoRECOPtGen_;
      TH1F* MatchedDMFindRECOPtGen_;
      TH1F* MatchedRECOPtGen_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOPtGen_;
      TGraphAsymmErrors* FinalEffMedIsoRECOPtGen_;
      TGraphAsymmErrors* FinalEffTightIsoRECOPtGen_;
      TGraphAsymmErrors* FinalEffDMFindRECOPtGen_;

      TH1F* MatchedLooseIsoCJPtGen_;
      TH1F* MatchedMedIsoCJPtGen_;
      TH1F* MatchedTightIsoCJPtGen_;
      TH1F* MatchedDMFindCJPtGen_;
      TH1F* MatchedCJPtGen_;
      TGraphAsymmErrors* FinalEffLooseIsoCJPtGen_;
      TGraphAsymmErrors* FinalEffMedIsoCJPtGen_;
      TGraphAsymmErrors* FinalEffTightIsoCJPtGen_;
      TGraphAsymmErrors* FinalEffDMFindCJPtGen_;

      TGraphAsymmErrors* FinalEffLooseIsoSAMEPt_;
      TGraphAsymmErrors* FinalEffMedIsoSAMEPt_;
      TGraphAsymmErrors* FinalEffTightIsoSAMEPt_;
      TGraphAsymmErrors* FinalEffDMFindSAMEPt_;
      TGraphAsymmErrors* FinalEffLooseIsoSAMEdR_;
      TGraphAsymmErrors* FinalEffMedIsoSAMEdR_;
      TGraphAsymmErrors* FinalEffTightIsoSAMEdR_;
      TGraphAsymmErrors* FinalEffDMFindSAMEdR_;
      TGraphAsymmErrors* FinalEffLooseIsoSAMEPtGen_;
      TGraphAsymmErrors* FinalEffMedIsoSAMEPtGen_;
      TGraphAsymmErrors* FinalEffTightIsoSAMEPtGen_;
      TGraphAsymmErrors* FinalEffDMFindSAMEPtGen_;


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
GGHAnalyzer::GGHAnalyzer(const edm::ParameterSet& iConfig):
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  genParticleTag_(consumes<vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleTag"))),
  akJetTag_(consumes<vector<reco::PFJet> >(iConfig.getParameter<edm::InputTag>("akJetTag"))),
  muonsTag_(consumes<vector<reco::Muon> >(iConfig.getParameter<edm::InputTag>("muonsTag"))),
  muonMapTag_(consumes<edm::ValueMap<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > >(iConfig.getParameter<edm::InputTag>("muonMapTag"))),
  jetValMapTag_(consumes<edm::ValueMap<reco::PFJetRef> >(iConfig.getParameter<edm::InputTag>("jetValMapTag"))),
  tauRECOTag_(consumes<vector<reco::PFTau> >(iConfig.getParameter<edm::InputTag>("tauRECOTag"))),
  tauCJTag_(consumes<vector<reco::PFTau> >(iConfig.getParameter<edm::InputTag>("tauCJTag"))),
  pizerosTag_(consumes<vector<reco::RecoTauPiZero> >(iConfig.getParameter<edm::InputTag>("pizerosTag"))),
  looseIsoTagCJ_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("looseIsoTagCJ"))),
  medIsoTagCJ_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("medIsoTagCJ"))),
  tightIsoTagCJ_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("tightIsoTagCJ"))),
  decayModeFindingTagCJ_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("decayModeFindingTagCJ"))),
  looseIsoTagRECO_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("looseIsoTagRECO"))),
  medIsoTagRECO_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("medIsoTagRECO"))),
  tightIsoTagRECO_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("tightIsoTagRECO"))),
  decayModeFindingTagRECO_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("decayModeFindingTagRECO")))
{
  reset(false);    
}//GGHAnalyzer



GGHAnalyzer::~GGHAnalyzer()
{
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void GGHAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  //get jet-muon map
  edm::Handle<edm::ValueMap<reco::MuonRefVector> > pMuonMap;
  iEvent.getByToken(muonMapTag_, pMuonMap);

  //get old-new jet map
  edm::Handle<edm::ValueMap<reco::PFJetRef> > pJetValMap;
  iEvent.getByToken(jetValMapTag_, pJetValMap);

  //Get RECO Taus particle collection
  edm::Handle<std::vector<reco::PFTau> > pTausRECO;
  iEvent.getByToken(tauRECOTag_, pTausRECO);

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

  NTausRECOvsCLEANJETS_->Fill(pTausRECO->size(), pTausCJ->size() );
  reco::GenParticleRef tau1Ref, tau2Ref;
  int checkGenMuHad = false;
  for (reco::GenParticleCollection::const_iterator iGenParticle = pGenParts->begin(); iGenParticle != pGenParts->end(); ++iGenParticle)
  {
    if(iGenParticle->pdgId() == 36 && iGenParticle->numberOfDaughters() == 2 && fabs(iGenParticle->daughter(0)->pdgId() ) == 15 )
    {
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

      if (VariousFunctions::tauDecayMode(tau1Ref) == 7  && VariousFunctions::tauDecayMode(tau2Ref) < 5 )
      {
        NEvents_->Fill(2);
        checkGenMuHad = 2;
      }//if tau2 Ref
      else if (VariousFunctions::tauDecayMode(tau2Ref) == 7  && VariousFunctions::tauDecayMode(tau1Ref) < 5 ) 
      {
        NEvents_->Fill(2);
        checkGenMuHad = 1;
      }//if
    }//ir found a1 and it decayed to taus
  }//for iGen Particle

  double dPhi_gen = reco::deltaPhi(tau1Ref->phi(), tau2Ref->phi() ), dR_tauMu_gen = sqrt( (tau1Ref->eta() - tau2Ref->eta())*(tau1Ref->eta() - tau2Ref->eta())  +  dPhi_gen * dPhi_gen );
  reco::LeafCandidate::LorentzVector GenTau1Visible = VariousFunctions::sumTauP4(tau1Ref, VariousFunctions::tauDecayMode(tau1Ref), false); 
  reco::LeafCandidate::LorentzVector GenTau2Visible = VariousFunctions::sumTauP4(tau2Ref, VariousFunctions::tauDecayMode(tau2Ref), false);  
  std::cout << "Tauhad is tauRef #" << checkGenMuHad << "   hpsPFTauProducer CleanJets size()= " << pTausCJ->size() << "  hpsPFTauProducer RECO size()= " << pTausRECO->size() << std::endl;

/////////////////////// 
// Analyze
/////////////////////// 
  double dR_GenMatchCJ1 = 1000, dR_GenMatchCJ2 = 1000, dR_GenMatchCJTEMP1 = 1000, dR_GenMatchCJTEMP2 = 1000, Pt_GenMatchCJ1 = -1, Pt_GenMatchCJ2 = -1;
  unsigned int DMCJ1 = 0, DMCJ2 = 0, MedIsoCJ1 = 0, MedIsoCJ2 = 0, LooseIsoCJ1 = 0, LooseIsoCJ2 = 0, TightIsoCJ1 = 0, TightIsoCJ2 = 0;
  double nConstituents = 1000;
  for (std::vector<reco::PFTau>::const_iterator iTauCJ = pTausCJ->begin(); iTauCJ != pTausCJ->end(); ++iTauCJ)
  {
    const reco::PFJetRef& tauJetRef = (*iTauCJ).jetRef();
    const reco::PFJetRef& tauRECOJetRef = (*pJetValMap)[tauJetRef];
    const reco::MuonRefVector& removedMuons = (*pMuonMap)[tauJetRef];
    std::vector<reco::PFCandidatePtr> JetPFCands = tauJetRef->getPFConstituents();
    nConstituents = JetPFCands.size();
    NMuRemoved_->Fill(removedMuons.size() );

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
 	GenDiTaudRvsCJDiTaudR_->Fill(dR_tauMu_gen, dR_tauMu);       
	break;
      }//if
    }//for itre

    //Gen Match the CleanJets taus
    double dPhi1 = reco::deltaPhi(tau1Ref->phi(), iTauCJ->phi() ), dPhi2 = reco::deltaPhi(tau2Ref->phi(), iTauCJ->phi() );
    dR_GenMatchCJTEMP1 = sqrt( (tau1Ref->eta() - iTauCJ->eta() ) * (tau1Ref->eta() - iTauCJ->eta() )  +  dPhi1 * dPhi1);
    dR_GenMatchCJTEMP2 = sqrt( (tau2Ref->eta() - iTauCJ->eta() ) * (tau2Ref->eta() - iTauCJ->eta() )  +  dPhi2 * dPhi2); 
    reco::PFTauRef PFTauRef(pTausCJ, iTauCJ - pTausCJ->begin()); 
    //first 2 checks are basic requirements, the 2 in the "or" are making sure it doesn't matche better with the second tau
    if (dR_GenMatchCJTEMP1 < .4 && dR_GenMatchCJTEMP1 < dR_GenMatchCJ1 && checkGenMuHad == 1) // && (dR_GenMatchCJTEMP1 < dR_GenMatchCJTEMP2 || dR_GenMatchCJTEMP2 > dR_GenMatchCJ2) )
    {
      dR_GenMatchCJ1 = dR_GenMatchCJTEMP1;
      Pt_GenMatchCJ1 = iTauCJ->pt();
      DMCJ1 = (*pDMFindingCJ)[PFTauRef];
      MedIsoCJ1 = (*pMedIsoDiscCJ)[PFTauRef];
      LooseIsoCJ1 = (*pLooseIsoDiscCJ)[PFTauRef];
      TightIsoCJ1 = (*pTightIsoDiscCJ)[PFTauRef];
    }//if iTauCJ matches tau1Ref 
    else if (dR_GenMatchCJTEMP2 < .4 && dR_GenMatchCJTEMP2 < dR_GenMatchCJ2  && checkGenMuHad == 2)// &&(dR_GenMatchCJTEMP2 < dR_GenMatchCJTEMP1 || dR_GenMatchCJTEMP1 > dR_GenMatchCJ1) )
    {
      dR_GenMatchCJ2 = dR_GenMatchCJTEMP2;
      Pt_GenMatchCJ2 = iTauCJ->pt();
      DMCJ2 = (*pDMFindingCJ)[PFTauRef];
      MedIsoCJ2 = (*pMedIsoDiscCJ)[PFTauRef];
      LooseIsoCJ2 = (*pLooseIsoDiscCJ)[PFTauRef];      
      TightIsoCJ2 = (*pTightIsoDiscCJ)[PFTauRef];
    }//else

    std::cout << "\t\tdR_GenMatchCJ1= " << dR_GenMatchCJ1 << "  \tdR_GenMatchCJTEMP1= " << dR_GenMatchCJTEMP1 << "  \tPtDiff1= " << abs(Pt_GenMatchCJ1 - iTauCJ->pt() ) << std::endl;
    std::cout << "\t\tdR_GenMatchCJ2= " << dR_GenMatchCJ2 << "  \tdR_GenMatchCJTEMP2= " << dR_GenMatchCJTEMP2 << "  \tPtDiff2= " << abs(Pt_GenMatchCJ2 - iTauCJ->pt() ) << std::endl;
  }//for auto iTauCJ

  if (dR_GenMatchCJ2 != 1000 && checkGenMuHad == 2)
  {
std::cout << "dR_GenMatchCJ2= " << dR_GenMatchCJ2 << "  DMCJ2= " << DMCJ2 << "  TightIsoCJ2= " << TightIsoCJ2 << "  MedIsoCJ2= " << MedIsoCJ2 << "  LooseIsoCJ2= " << LooseIsoCJ2 << std::endl;
    NEvents_->Fill(3);
    MatchedCJPt_->Fill(Pt_GenMatchCJ2 );
    if (DMCJ2 == 1)
      MatchedDMFindCJPt_->Fill(Pt_GenMatchCJ2 );
    if (TightIsoCJ2 == 1 && DMCJ2 == 1)
      MatchedTightIsoCJPt_->Fill(Pt_GenMatchCJ2 );
    if (MedIsoCJ2 == 1 && DMCJ2 == 1)
      MatchedMedIsoCJPt_->Fill(Pt_GenMatchCJ2 );
    if (LooseIsoCJ2 == 1 && DMCJ2 == 1)
      MatchedLooseIsoCJPt_->Fill(Pt_GenMatchCJ2 );
  
    MatchedCJdR_->Fill(dR_tauMu_gen );
    if (DMCJ2 == 1)
      MatchedDMFindCJdR_->Fill(dR_tauMu_gen );
    if (TightIsoCJ2 == 1 && DMCJ2 == 1)
      MatchedTightIsoCJdR_->Fill(dR_tauMu_gen );
    if (MedIsoCJ2 == 1 && DMCJ2 == 1)
      MatchedMedIsoCJdR_->Fill(dR_tauMu_gen );
    if (LooseIsoCJ2 == 1 && DMCJ2 == 1)
      MatchedLooseIsoCJdR_->Fill(dR_tauMu_gen );

    MatchedCJPtGen_->Fill(GenTau2Visible.Pt() );
    if (DMCJ2 == 1)
      MatchedDMFindCJPtGen_->Fill(GenTau2Visible.Pt() );
    if (TightIsoCJ2 == 1 && DMCJ2 == 1)
      MatchedTightIsoCJPtGen_->Fill(GenTau2Visible.Pt() );
    if (MedIsoCJ2 == 1 && DMCJ2 == 1)
      MatchedMedIsoCJPtGen_->Fill(GenTau2Visible.Pt() );
    if (LooseIsoCJ2 == 1 && DMCJ2 == 1)
      MatchedLooseIsoCJPtGen_->Fill(GenTau2Visible.Pt() );
  }//if GEN tau2Ref is the had in mu+had and it is matched to CleanJets Jet
  else if (dR_GenMatchCJ1 != 1000 && checkGenMuHad == 1)
  {
std::cout << "dR_GenMatchCJ1= " << dR_GenMatchCJ1 << "  DMCJ1= " << DMCJ1 << "  TightIsoCJ1= " << TightIsoCJ1 << "  MedIsoCJ1= " << MedIsoCJ1 << "  LooseIsoCJ1= " << LooseIsoCJ1 << std::endl;
    NEvents_->Fill(3);
    MatchedCJPt_->Fill(Pt_GenMatchCJ1 );
    if (DMCJ1 == 1)
      MatchedDMFindCJPt_->Fill(Pt_GenMatchCJ1 );
    if (TightIsoCJ1 == 1 && DMCJ1 == 1)
      MatchedTightIsoCJPt_->Fill(Pt_GenMatchCJ1 );
    if (MedIsoCJ1 == 1 && DMCJ1 == 1)
      MatchedMedIsoCJPt_->Fill(Pt_GenMatchCJ1 );
    if (LooseIsoCJ1 == 1 && DMCJ1 == 1)
      MatchedLooseIsoCJPt_->Fill(Pt_GenMatchCJ1 );

    MatchedCJdR_->Fill(dR_tauMu_gen );
    if (DMCJ1 == 1)
      MatchedDMFindCJdR_->Fill(dR_tauMu_gen );
    if (TightIsoCJ1 == 1 && DMCJ1 == 1)
      MatchedTightIsoCJdR_->Fill(dR_tauMu_gen );
    if (MedIsoCJ1 == 1 && DMCJ1 == 1)
      MatchedMedIsoCJdR_->Fill(dR_tauMu_gen );
    if (LooseIsoCJ1 == 1 && DMCJ1 == 1)
      MatchedLooseIsoCJdR_->Fill(dR_tauMu_gen );
    
    MatchedCJPtGen_->Fill(GenTau1Visible.Pt() );
    if (DMCJ1 == 1)
      MatchedDMFindCJPtGen_->Fill(GenTau1Visible.Pt() );
    if (TightIsoCJ1 == 1 && DMCJ1 == 1)
      MatchedTightIsoCJPtGen_->Fill(GenTau1Visible.Pt() );
    if (MedIsoCJ1 == 1 && DMCJ1 == 1)
      MatchedMedIsoCJPtGen_->Fill(GenTau1Visible.Pt() );
    if (LooseIsoCJ1 == 1 && DMCJ1 == 1)
      MatchedLooseIsoCJPtGen_->Fill(GenTau1Visible.Pt() );
  }//if GEN tau1Ref is the had in mu+had and it is matched to CleanJets Jet

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
    reco::PFTauRef PFTauRef(pTausRECO, iTauRECO - pTausRECO->begin());
    //first 2 checks are basic requirements, the 2 in the "or" are making sure it doesn't matche better with the second tau
    if (dR_GenMatchRECOTEMP1 < .4 && dR_GenMatchRECOTEMP1 < dR_GenMatchRECO1 && checkGenMuHad == 1) //&& (dR_GenMatchRECOTEMP1 < dR_GenMatchRECOTEMP2 || dR_GenMatchRECOTEMP2 > dR_GenMatchRECO2) )
    {
      dR_GenMatchRECO1 = dR_GenMatchRECOTEMP1;
      Pt_GenMatchRECO1 = iTauRECO->pt();
      DMRECO1 = (*pDMFindingRECO)[PFTauRef];
      MedIsoRECO1 = (*pMedIsoDiscRECO)[PFTauRef];
      LooseIsoRECO1 = (*pLooseIsoDiscRECO)[PFTauRef];
      TightIsoRECO1 = (*pTightIsoDiscRECO)[PFTauRef];
    }//if iTauRECO matches tau1Ref 
    else if (dR_GenMatchRECOTEMP2 < .4 && dR_GenMatchRECOTEMP2 < dR_GenMatchRECO2 && checkGenMuHad == 2) // && (dR_GenMatchRECOTEMP2 < dR_GenMatchRECOTEMP1 || dR_GenMatchRECOTEMP1 > dR_GenMatchRECO1) )
    {
      dR_GenMatchRECO2 = dR_GenMatchRECOTEMP2;
      Pt_GenMatchRECO2 = iTauRECO->pt();
      DMRECO2 = (*pDMFindingRECO)[PFTauRef];
      MedIsoRECO2 = (*pMedIsoDiscRECO)[PFTauRef];
      LooseIsoRECO2 = (*pLooseIsoDiscRECO)[PFTauRef]; 
      TightIsoRECO2 = (*pTightIsoDiscRECO)[PFTauRef];
    }//else

    std::cout << "\t\tdR_GenMatchRECO1= " << dR_GenMatchRECO1 << "  \tdR_GenMatchRECOTEMP1= " << dR_GenMatchRECOTEMP1 << "  \tPtDiff1= " << abs(Pt_GenMatchRECO1 - iTauRECO->pt() ) << std::endl;
    std::cout << "\t\tdR_GenMatchRECO2= " << dR_GenMatchRECO2 << "  \tdR_GenMatchRECOTEMP2= " << dR_GenMatchRECOTEMP2 << "  \tPtDiff2= " << abs(Pt_GenMatchRECO2 - iTauRECO->pt() ) << std::endl;
  }//iTauRECO

  if (dR_GenMatchRECO2 != 1000 && checkGenMuHad == 2)
  {
std::cout << "dR_GenMatchRECO2= " << dR_GenMatchRECO2 << "  DMRECO2= " << DMRECO2 << "  TightIsoRECO2= " << TightIsoRECO2 << "  MedIsoRECO2= " << MedIsoRECO2 << "  LooseIsoRECO2= " << LooseIsoRECO2 << std::endl;
    NEvents_->Fill(4);
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
  else if (dR_GenMatchRECO1 != 1000 && checkGenMuHad == 1)
  {
std::cout << "dR_GenMatchRECO1= " << dR_GenMatchRECO1 << "  DMRECO1= " << DMRECO1 << "  TightIsoRECO1= " << TightIsoRECO1 << "  MedIsoRECO1= " << MedIsoRECO1 << "  LooseIsoRECO1= " << LooseIsoRECO1 << std::endl;
   NEvents_->Fill(4);
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

 
}//End GGHAnalyzer::analyze


// ------------ method called once each job just before starting event loop  ------------
void GGHAnalyzer::beginJob()
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
      NEvents_->GetXaxis()->SetBinLabel(6, "");
      NEvents_->GetXaxis()->SetBinLabel(7, "");
  NMuRemoved_            = new TH1F("NMuRemoved"    , "", 9, -.5, 8.5);
  TauMuTauHaddR_       = new TH1F("TauMuTauHaddR"    , "", 100, 0, 100);
  NConstituentsCJ_        = new TH1F("NConstituentsCJ"    , "", 50, 0, 50);
  NTausRECOvsCLEANJETS_  = new TH2F("NTausRECOvsCLEANJETS" , "", 11, -.5, 10.5, 11, -.5, 10.5);
  GenDiTaudRvsCJDiTaudR_  = new TH2F("GenDiTaudRvsCJDiTaudR" , "", 50, 0, 10, 50, 0, 10);
  TauHadnConstvsPt_  = new TH2F("TauHadnConstvsPt" , "", 50, 0, 50, 30, 0, 100);



  MatchedLooseIsoRECOPt_  = new TH1F("MatchedLooseIsoRECOPt"    , "", 30, 0, 90.0);
  MatchedMedIsoRECOPt_    = new TH1F("MatchedMedIsoRECOPt", "", 30, 0, 90);
  MatchedTightIsoRECOPt_    = new TH1F("MatchedTightIsoRECOPt", "", 30, 0, 90);
  MatchedDMFindRECOPt_    = new TH1F("MatchedDMFindRECOPt"    , "", 30, 0, 90);
  MatchedRECOPt_          = new TH1F("MatchedRECOPt"    , "", 30, 0, 90);

  MatchedLooseIsoCJPt_  = new TH1F("MatchedLooseIsoCJPt"    , "", 30, 0, 90.0);
  MatchedMedIsoCJPt_    = new TH1F("MatchedMedIsoCJPt", "", 30, 0, 90);
  MatchedTightIsoCJPt_    = new TH1F("MatchedTightIsoCJPt", "", 30, 0, 90);
  MatchedDMFindCJPt_    = new TH1F("MatchedDMFindCJPt"    , "", 30, 0, 90);
  MatchedCJPt_          = new TH1F("MatchedCJPt"    , "", 30, 0, 90);

  FinalEffLooseIsoRECOPt_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOPt_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOPt_ = new TGraphAsymmErrors(30);
  FinalEffDMFindRECOPt_ = new TGraphAsymmErrors(30);

  FinalEffLooseIsoCJPt_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoCJPt_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoCJPt_ = new TGraphAsymmErrors(30);
  FinalEffDMFindCJPt_ = new TGraphAsymmErrors(30);



  MatchedLooseIsoRECOdR_  = new TH1F("MatchedLooseIsoRECOdR"    , "", 30, 0, 4.0);
  MatchedMedIsoRECOdR_    = new TH1F("MatchedMedIsoRECOdR", "", 30, 0, 4);
  MatchedTightIsoRECOdR_    = new TH1F("MatchedTightIsoRECOdR", "", 30, 0, 4);
  MatchedDMFindRECOdR_    = new TH1F("MatchedDMFindRECOdR"    , "", 30, 0, 4);
  MatchedRECOdR_          = new TH1F("MatchedRECOdR"    , "", 30, 0, 4);

  MatchedLooseIsoCJdR_  = new TH1F("MatchedLooseIsoCJdR"    , "", 30, 0, 4.0);
  MatchedMedIsoCJdR_    = new TH1F("MatchedMedIsoCJdR", "", 30, 0, 4);
  MatchedTightIsoCJdR_    = new TH1F("MatchedTightIsoCJdR", "", 30, 0, 4);
  MatchedDMFindCJdR_    = new TH1F("MatchedDMFindCJdR"    , "", 30, 0, 4);
  MatchedCJdR_          = new TH1F("MatchedCJdR"    , "", 30, 0, 4);

  FinalEffLooseIsoRECOdR_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOdR_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOdR_ = new TGraphAsymmErrors(30);
  FinalEffDMFindRECOdR_ = new TGraphAsymmErrors(30);

  FinalEffLooseIsoCJdR_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoCJdR_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoCJdR_ = new TGraphAsymmErrors(30);
  FinalEffDMFindCJdR_ = new TGraphAsymmErrors(30);



  MatchedLooseIsoRECOPtGen_  = new TH1F("MatchedLooseIsoRECOPtGen"    , "", 30, 0, 90.0);
  MatchedMedIsoRECOPtGen_    = new TH1F("MatchedMedIsoRECOPtGen", "", 30, 0, 90);
  MatchedTightIsoRECOPtGen_    = new TH1F("MatchedTightIsoRECOPtGen", "", 30, 0, 90);
  MatchedDMFindRECOPtGen_    = new TH1F("MatchedDMFindRECOPtGen"    , "", 30, 0, 90);
  MatchedRECOPtGen_          = new TH1F("MatchedRECOPtGen"    , "", 30, 0, 90);

  MatchedLooseIsoCJPtGen_  = new TH1F("MatchedLooseIsoCJPtGen"    , "", 30, 0, 90.0);
  MatchedMedIsoCJPtGen_    = new TH1F("MatchedMedIsoCJPtGen", "", 30, 0, 90);
  MatchedTightIsoCJPtGen_    = new TH1F("MatchedTightIsoCJPtGen", "", 30, 0, 90);
  MatchedDMFindCJPtGen_    = new TH1F("MatchedDMFindCJPtGen"    , "", 30, 0, 90);
  MatchedCJPtGen_          = new TH1F("MatchedCJPtGen"    , "", 30, 0, 90);

  FinalEffLooseIsoRECOPtGen_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOPtGen_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOPtGen_ = new TGraphAsymmErrors(30);
  FinalEffDMFindRECOPtGen_ = new TGraphAsymmErrors(30);

  FinalEffLooseIsoCJPtGen_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoCJPtGen_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoCJPtGen_ = new TGraphAsymmErrors(30);
  FinalEffDMFindCJPtGen_ = new TGraphAsymmErrors(30);


  FinalEffLooseIsoSAMEPt_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoSAMEPt_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoSAMEPt_ = new TGraphAsymmErrors(30);
  FinalEffDMFindSAMEPt_ = new TGraphAsymmErrors(30);
  FinalEffLooseIsoSAMEdR_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoSAMEdR_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoSAMEdR_ = new TGraphAsymmErrors(30);
  FinalEffDMFindSAMEdR_ = new TGraphAsymmErrors(30);
  FinalEffLooseIsoSAMEPtGen_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoSAMEPtGen_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoSAMEPtGen_ = new TGraphAsymmErrors(30);
  FinalEffDMFindSAMEPtGen_ = new TGraphAsymmErrors(30);

}

// ------------ method called once each job just after ending the event loop  ------------
void GGHAnalyzer::endJob()
{
  //Make the Canvases
  TCanvas NEventsCanvas("NEvents","",600,600);
  TCanvas NMuRemovedCanvas("NMuRemoved","",600,600);
  TCanvas TauMuTauHaddRCanvas("TauMuTauHaddR","",600,600);
  TCanvas NConstituentsCJCanvas("NConstituentsCJ","",600,600);
  TCanvas NTausRECOvsCLEANJETSCanvas("NTausRECOvsCLEANJETS","",600,600);
  TCanvas GenDiTaudRvsCJDiTaudRCanvas("GenDiTaudRvsCJDiTaudR","",600,600);
  TCanvas TauHadnConstvsPtCanvas("TauHadnConstvsPt","",600,600);

  TCanvas MatchedLooseIsoRECOPtCanvas("MatchedLooseIsoRECOPt","",600,600);
  TCanvas MatchedMedIsoRECOPtCanvas("MatchedMedIsoRECOPt","",600,600);
  TCanvas MatchedTightIsoRECOPtCanvas("MatchedTightIsoRECOPt","",600,600);
  TCanvas MatchedDMFindRECOPtCanvas("MatchedDMFindRECOPt","",600,600);
  TCanvas MatchedRECOPtCanvas("MatchedRECOPt","",600,600);
  TCanvas FinalEffLooseIsoRECOPtCanvas("FinalEffLooseIsoRECOPt","",600,600);
  TCanvas FinalEffMedIsoRECOPtCanvas("FinalEffMedIsoRECOPt","",600,600);
  TCanvas FinalEffTightIsoRECOPtCanvas("FinalEffTightIsoRECOPt","",600,600);
  TCanvas FinalEffDMFindRECOPtCanvas("FinalEffDMFindRECOPt","",600,600);

  TCanvas MatchedLooseIsoCJPtCanvas("MatchedLooseIsoCJPt","",600,600);
  TCanvas MatchedMedIsoCJPtCanvas("MatchedMedIsoCJPt","",600,600);
  TCanvas MatchedTightIsoCJPtCanvas("MatchedTightIsoCJPt","",600,600);
  TCanvas MatchedDMFindCJPtCanvas("MatchedDMFindCJPt","",600,600);
  TCanvas MatchedCJPtCanvas("MatchedCJPt","",600,600);
  TCanvas FinalEffLooseIsoCJPtCanvas("FinalEffLooseIsoCJPt","",600,600);
  TCanvas FinalEffMedIsoCJPtCanvas("FinalEffMedIsoCJPt","",600,600);
  TCanvas FinalEffTightIsoCJPtCanvas("FinalEffTightIsoCJPt","",600,600);
  TCanvas FinalEffDMFindCJPtCanvas("FinalEffDMFindCJPt","",600,600);



  TCanvas MatchedLooseIsoRECOdRCanvas("MatchedLooseIsoRECOdR","",600,600);
  TCanvas MatchedMedIsoRECOdRCanvas("MatchedMedIsoRECOdR","",600,600);
  TCanvas MatchedTightIsoRECOdRCanvas("MatchedTightIsoRECOdR","",600,600);
  TCanvas MatchedDMFindRECOdRCanvas("MatchedDMFindRECOdR","",600,600);
  TCanvas MatchedRECOdRCanvas("MatchedRECOdR","",600,600);
  TCanvas FinalEffLooseIsoRECOdRCanvas("FinalEffLooseIsoRECOdR","",600,600);
  TCanvas FinalEffMedIsoRECOdRCanvas("FinalEffMedIsoRECOdR","",600,600);
  TCanvas FinalEffTightIsoRECOdRCanvas("FinalEffTightIsoRECOdR","",600,600);
  TCanvas FinalEffDMFindRECOdRCanvas("FinalEffDMFindRECOdR","",600,600);

  TCanvas MatchedLooseIsoCJdRCanvas("MatchedLooseIsoCJdR","",600,600);
  TCanvas MatchedMedIsoCJdRCanvas("MatchedMedIsoCJdR","",600,600);
  TCanvas MatchedTightIsoCJdRCanvas("MatchedTightIsoCJdR","",600,600);
  TCanvas MatchedDMFindCJdRCanvas("MatchedDMFindCJdR","",600,600);
  TCanvas MatchedCJdRCanvas("MatchedCJdR","",600,600);
  TCanvas FinalEffLooseIsoCJdRCanvas("FinalEffLooseIsoCJdR","",600,600);
  TCanvas FinalEffMedIsoCJdRCanvas("FinalEffMedIsoCJdR","",600,600);
  TCanvas FinalEffTightIsoCJdRCanvas("FinalEffTightIsoCJdR","",600,600);
  TCanvas FinalEffDMFindCJdRCanvas("FinalEffDMFindCJdR","",600,600);


  TCanvas MatchedLooseIsoRECOPtGenCanvas("MatchedLooseIsoRECOPtGen","",600,600);
  TCanvas MatchedMedIsoRECOPtGenCanvas("MatchedMedIsoRECOPtGen","",600,600);
  TCanvas MatchedTightIsoRECOPtGenCanvas("MatchedTightIsoRECOPtGen","",600,600);
  TCanvas MatchedDMFindRECOPtGenCanvas("MatchedDMFindRECOPtGen","",600,600);
  TCanvas MatchedRECOPtGenCanvas("MatchedRECOPtGen","",600,600);
  TCanvas FinalEffLooseIsoRECOPtGenCanvas("FinalEffLooseIsoRECOPtGen","",600,600);
  TCanvas FinalEffMedIsoRECOPtGenCanvas("FinalEffMedIsoRECOPtGen","",600,600);
  TCanvas FinalEffTightIsoRECOPtGenCanvas("FinalEffTightIsoRECOPtGen","",600,600);
  TCanvas FinalEffDMFindRECOPtGenCanvas("FinalEffDMFindRECOPtGen","",600,600);

  TCanvas MatchedLooseIsoCJPtGenCanvas("MatchedLooseIsoCJPtGen","",600,600);
  TCanvas MatchedMedIsoCJPtGenCanvas("MatchedMedIsoCJPtGen","",600,600);
  TCanvas MatchedTightIsoCJPtGenCanvas("MatchedTightIsoCJPtGen","",600,600);
  TCanvas MatchedDMFindCJPtGenCanvas("MatchedDMFindCJPtGen","",600,600);
  TCanvas MatchedCJPtGenCanvas("MatchedCJPtGen","",600,600);
  TCanvas FinalEffLooseIsoCJPtGenCanvas("FinalEffLooseIsoCJPtGen","",600,600);
  TCanvas FinalEffMedIsoCJPtGenCanvas("FinalEffMedIsoCJPtGen","",600,600);
  TCanvas FinalEffTightIsoCJPtGenCanvas("FinalEffTightIsoCJPtGen","",600,600);
  TCanvas FinalEffDMFindCJPtGenCanvas("FinalEffDMFindCJPtGen","",600,600);

  TCanvas FinalEffLooseIsoSAMEPtCanvas("FinalEffLooseIsoSAMEPt","",600,600);
  TCanvas FinalEffMedIsoSAMEPtCanvas("FinalEffMedIsoSAMEPt","",600,600);
  TCanvas FinalEffTightIsoSAMEPtCanvas("FinalEffTightIsoSAMEPt","",600,600);
  TCanvas FinalEffDMFindSAMEPtCanvas("FinalEffDMFindSAMEPt","",600,600);
  TCanvas FinalEffLooseIsoSAMEdRCanvas("FinalEffLooseIsoSAMEdR","",600,600);
  TCanvas FinalEffMedIsoSAMEdRCanvas("FinalEffMedIsoSAMEdR","",600,600);
  TCanvas FinalEffTightIsoSAMEdRCanvas("FinalEffTightIsoSAMEdR","",600,600);
  TCanvas FinalEffDMFindSAMEdRCanvas("FinalEffDMFindSAMEdR","",600,600);
  TCanvas FinalEffLooseIsoSAMEPtGenCanvas("FinalEffLooseIsoSAMEPtGen","",600,600);
  TCanvas FinalEffMedIsoSAMEPtGenCanvas("FinalEffMedIsoSAMEPtGen","",600,600);
  TCanvas FinalEffTightIsoSAMEPtGenCanvas("FinalEffTightIsoSAMEPtGen","",600,600);
  TCanvas FinalEffDMFindSAMEPtGenCanvas("FinalEffDMFindSAMEPtGen","",600,600);


  FinalEffLooseIsoRECOPt_->Divide(MatchedLooseIsoRECOPt_, MatchedRECOPt_);
  FinalEffMedIsoRECOPt_->Divide(MatchedMedIsoRECOPt_,     MatchedRECOPt_);
  FinalEffTightIsoRECOPt_->Divide(MatchedTightIsoRECOPt_, MatchedRECOPt_);
  FinalEffDMFindRECOPt_->Divide(MatchedDMFindRECOPt_,     MatchedRECOPt_);
  FinalEffLooseIsoCJPt_->Divide(MatchedLooseIsoCJPt_, MatchedCJPt_);
  FinalEffMedIsoCJPt_->Divide(MatchedMedIsoCJPt_,     MatchedCJPt_);
  FinalEffTightIsoCJPt_->Divide(MatchedTightIsoCJPt_, MatchedCJPt_);
  FinalEffDMFindCJPt_->Divide(MatchedDMFindCJPt_,     MatchedCJPt_);



  FinalEffLooseIsoRECOdR_->Divide(MatchedLooseIsoRECOdR_, MatchedRECOdR_);
  FinalEffMedIsoRECOdR_->Divide(MatchedMedIsoRECOdR_,     MatchedRECOdR_);
  FinalEffTightIsoRECOdR_->Divide(MatchedTightIsoRECOdR_, MatchedRECOdR_);
  FinalEffDMFindRECOdR_->Divide(MatchedDMFindRECOdR_,     MatchedRECOdR_);
  FinalEffLooseIsoCJdR_->Divide(MatchedLooseIsoCJdR_, MatchedCJdR_);
  FinalEffMedIsoCJdR_->Divide(MatchedMedIsoCJdR_,     MatchedCJdR_);
  FinalEffTightIsoCJdR_->Divide(MatchedTightIsoCJdR_, MatchedCJdR_);
  FinalEffDMFindCJdR_->Divide(MatchedDMFindCJdR_,     MatchedCJdR_);


  FinalEffLooseIsoRECOPtGen_->Divide(MatchedLooseIsoRECOPtGen_, MatchedRECOPtGen_);
  FinalEffMedIsoRECOPtGen_->Divide(MatchedMedIsoRECOPtGen_,     MatchedRECOPtGen_);
  FinalEffTightIsoRECOPtGen_->Divide(MatchedTightIsoRECOPtGen_, MatchedRECOPtGen_);
  FinalEffDMFindRECOPtGen_->Divide(MatchedDMFindRECOPtGen_,     MatchedRECOPtGen_);
  FinalEffLooseIsoCJPtGen_->Divide(MatchedLooseIsoCJPtGen_, MatchedCJPtGen_);
  FinalEffMedIsoCJPtGen_->Divide(MatchedMedIsoCJPtGen_,     MatchedCJPtGen_);
  FinalEffTightIsoCJPtGen_->Divide(MatchedTightIsoCJPtGen_, MatchedCJPtGen_);
  FinalEffDMFindCJPtGen_->Divide(MatchedDMFindCJPtGen_,     MatchedCJPtGen_);

std::cout << "<----------------Declared Canvases-------------->" << std::endl;

  //Format the 1D plots and Draw (canvas, hist, grid, log y, log z, color, size, style, xAxisTitle, xTitleSize, xLabelSize, xTitleOffSet, yAxisTitle, yTitleSize, yLabelSize, yTitleOffset)
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		NEventsCanvas, NEvents_, 1, 0, 0, kBlack, 7, 20, "", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		NMuRemovedCanvas, NMuRemoved_, 1, 0, 0, kBlack, 7, 20, "N #mu Removed", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		TauMuTauHaddRCanvas, TauMuTauHaddR_, 1, 0, 0, kBlack, 7, 20, "#DeltaR(#tau_{mu} + #tau_{H})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		NConstituentsCJCanvas, NConstituentsCJ_, 1, 0, 0, kBlack, 7, 20, "Tau Mom PDGID", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist2D(
		NTausRECOvsCLEANJETSCanvas, NTausRECOvsCLEANJETS_, 1, 0, 0, kBlack, 7, 20, "nRECO #tau's", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(
		GenDiTaudRvsCJDiTaudRCanvas, GenDiTaudRvsCJDiTaudR_, 1, 0, 0, kBlack, 7, 20, "gen #DeltaR(#tau#tau)", .04, .04, 1.1, "CJ #DeltaR(#tau#tau)", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(
                TauHadnConstvsPtCanvas, TauHadnConstvsPt_, 1, 0, 0, kBlack, 7, 20, "# Constituents(#tau_{H})", .04, .04, 1.1, "Pt(#tau_{H})", .04, .04, 1.6, "", .04, .04, 1.0);

  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedLooseIsoRECOPtCanvas, MatchedLooseIsoRECOPt_, 1, 0, 0, kBlack, 7, 20, "Pt(RECO Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedMedIsoRECOPtCanvas, MatchedMedIsoRECOPt_, 1, 0, 0, kBlack, 7, 20, "Pt(RECO Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedTightIsoRECOPtCanvas, MatchedTightIsoRECOPt_, 1, 0, 0, kBlack, 7, 20, "Pt(RECO Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedDMFindRECOPtCanvas, MatchedDMFindRECOPt_, 1, 0, 0, kBlack, 7, 20, "Pt(RECO Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedRECOPtCanvas, MatchedRECOPt_, 1, 0, 0, kBlack, 7, 20, "Pt(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffLooseIsoRECOPtCanvas, FinalEffLooseIsoRECOPt_, 1, 0, 0, kBlack, 1, 20, "Pt(RECO Matched + Loose Iso + DM) / Pt(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffMedIsoRECOPtCanvas, FinalEffMedIsoRECOPt_, 1, 0, 0, kBlack, 1, 20, "Pt(RECO Matched + Med Iso + DM) / Pt(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffTightIsoRECOPtCanvas, FinalEffTightIsoRECOPt_, 1, 0, 0, kBlack, 1, 20, "Pt(RECO Matched + Tight Iso + DM) / Pt(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffDMFindRECOPtCanvas, FinalEffDMFindRECOPt_, 1, 0, 0, kBlack, 1, 20, "Pt(RECO Matched + DecayModeFinding) / Pt(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedLooseIsoCJPtCanvas, MatchedLooseIsoCJPt_, 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedMedIsoCJPtCanvas, MatchedMedIsoCJPt_, 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedTightIsoCJPtCanvas, MatchedTightIsoCJPt_, 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedDMFindCJPtCanvas, MatchedDMFindCJPt_, 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedCJPtCanvas, MatchedCJPt_, 1, 0, 0, kBlack, 7, 20, "Pt(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffLooseIsoCJPtCanvas, FinalEffLooseIsoCJPt_, 1, 0, 0, kBlack, 1, 20, "Pt(CleanJets Matched + Loose Iso + DM) / Pt(CLeanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffMedIsoCJPtCanvas, FinalEffMedIsoCJPt_, 1, 0, 0, kBlack, 1, 20, "Pt(CleanJets Matched + Med Iso + DM) / Pt(CLeanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffTightIsoCJPtCanvas, FinalEffTightIsoCJPt_, 1, 0, 0, kBlack, 1, 20, "Pt(CleanJets Matched + Tight Iso + DM) / Pt(CLeanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffDMFindCJPtCanvas, FinalEffDMFindCJPt_, 1, 0, 0, kBlack, 1, 20, "Pt(CleanJets Matched + DecayModeFinding) / Pt(CLeanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);



  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedLooseIsoRECOdRCanvas, MatchedLooseIsoRECOdR_, 1, 0, 0, kBlack, 7, 20, "#DeltaR(RECO Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedMedIsoRECOdRCanvas, MatchedMedIsoRECOdR_, 1, 0, 0, kBlack, 7, 20, "#DeltaR(RECO Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedTightIsoRECOdRCanvas, MatchedTightIsoRECOdR_, 1, 0, 0, kBlack, 7, 20, "#DeltaR(RECO Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedDMFindRECOdRCanvas, MatchedDMFindRECOdR_, 1, 0, 0, kBlack, 7, 20, "#DeltaR(RECO Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedRECOdRCanvas, MatchedRECOdR_, 1, 0, 0, kBlack, 7, 20, "#DeltaR(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffLooseIsoRECOdRCanvas, FinalEffLooseIsoRECOdR_, 1, 0, 0, kBlack, 1, 20, "#DeltaR(RECO Matched + Loose Iso + DM) / #DeltaR(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffMedIsoRECOdRCanvas, FinalEffMedIsoRECOdR_, 1, 0, 0, kBlack, 1, 20, "#DeltaR(RECO Matched + Med Iso + DM) / #DeltaR(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffTightIsoRECOdRCanvas, FinalEffTightIsoRECOdR_, 1, 0, 0, kBlack, 1, 20, "#DeltaR(RECO Matched + Tight Iso + DM) / #DeltaR(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffDMFindRECOdRCanvas, FinalEffDMFindRECOdR_, 1, 0, 0, kBlack, 1, 20, "#DeltaR(RECO Matched + DecayModeFinding) / #DeltaR(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedLooseIsoCJdRCanvas, MatchedLooseIsoCJdR_, 1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedMedIsoCJdRCanvas, MatchedMedIsoCJdR_, 1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedTightIsoCJdRCanvas, MatchedTightIsoCJdR_, 1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedDMFindCJdRCanvas, MatchedDMFindCJdR_, 1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedCJdRCanvas, MatchedCJdR_, 1, 0, 0, kBlack, 7, 20, "#DeltaR(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffLooseIsoCJdRCanvas, FinalEffLooseIsoCJdR_, 1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Loose Iso + DM) / #DeltaR(CLeanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffMedIsoCJdRCanvas, FinalEffMedIsoCJdR_, 1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Med Iso + DM) / #DeltaR(CLeanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffTightIsoCJdRCanvas, FinalEffTightIsoCJdR_, 1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Tight Iso + DM) / #DeltaR(CLeanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffDMFindCJdRCanvas, FinalEffDMFindCJdR_, 1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + DecayModeFinding) / #DeltaR(CLeanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);



  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedLooseIsoRECOPtGenCanvas, MatchedLooseIsoRECOPtGen_, 1, 0, 0, kBlack, 7, 20, "PtGen(RECO Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedMedIsoRECOPtGenCanvas, MatchedMedIsoRECOPtGen_, 1, 0, 0, kBlack, 7, 20, "PtGen(RECO Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedTightIsoRECOPtGenCanvas, MatchedTightIsoRECOPtGen_, 1, 0, 0, kBlack, 7, 20, "PtGen(RECO Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedDMFindRECOPtGenCanvas, MatchedDMFindRECOPtGen_, 1, 0, 0, kBlack, 7, 20, "PtGen(RECO Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedRECOPtGenCanvas, MatchedRECOPtGen_, 1, 0, 0, kBlack, 7, 20, "PtGen(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffLooseIsoRECOPtGenCanvas, FinalEffLooseIsoRECOPtGen_, 1, 0, 0, kBlack, 1, 20, "PtGen(RECO Matched + Loose Iso + DM) / PtGen(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffMedIsoRECOPtGenCanvas, FinalEffMedIsoRECOPtGen_, 1, 0, 0, kBlack, 1, 20, "PtGen(RECO Matched + Med Iso + DM) / PtGen(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffTightIsoRECOPtGenCanvas, FinalEffTightIsoRECOPtGen_, 1, 0, 0, kBlack, 1, 20, "PtGen(RECO Matched + Tight Iso + DM) / PtGen(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffDMFindRECOPtGenCanvas, FinalEffDMFindRECOPtGen_, 1, 0, 0, kBlack, 1, 20, "PtGen(RECO Matched + DecayModeFinding) / PtGen(RECO Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedLooseIsoCJPtGenCanvas, MatchedLooseIsoCJPtGen_, 1, 0, 0, kBlack, 7, 20, "PtGen(CleanJets Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedMedIsoCJPtGenCanvas, MatchedMedIsoCJPtGen_, 1, 0, 0, kBlack, 7, 20, "PtGen(CleanJets Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedTightIsoCJPtGenCanvas, MatchedTightIsoCJPtGen_, 1, 0, 0, kBlack, 7, 20, "PtGen(CleanJets Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedDMFindCJPtGenCanvas, MatchedDMFindCJPtGen_, 1, 0, 0, kBlack, 7, 20, "PtGen(CleanJets Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(
		MatchedCJPtGenCanvas, MatchedCJPtGen_, 1, 0, 0, kBlack, 7, 20, "PtGen(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffLooseIsoCJPtGenCanvas, FinalEffLooseIsoCJPtGen_, 1, 0, 0, kBlack, 1, 20, "PtGen(CleanJets Matched + Loose Iso + DM) / PtGen(CLeanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffMedIsoCJPtGenCanvas, FinalEffMedIsoCJPtGen_, 1, 0, 0, kBlack, 1, 20, "PtGen(CleanJets Matched + Med Iso + DM) / PtGen(CLeanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffTightIsoCJPtGenCanvas, FinalEffTightIsoCJPtGen_, 1, 0, 0, kBlack, 1, 20, "PtGen(CleanJets Matched + Tight Iso + DM) / PtGen(CLeanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(
		FinalEffDMFindCJPtGenCanvas, FinalEffDMFindCJPtGen_, 1, 0, 0, kBlack, 1, 20, "PtGen(CleanJets Matched + DecayModeFinding) / PtGen(CLeanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

  FinalEffLooseIsoCJPt_->SetLineColor(kRed);
  FinalEffMedIsoCJPt_->SetLineColor(kRed);
  FinalEffTightIsoCJPt_->SetLineColor(kRed);
  FinalEffDMFindCJPt_->SetLineColor(kRed);

  FinalEffLooseIsoCJdR_->SetLineColor(kRed);
  FinalEffMedIsoCJdR_->SetLineColor(kRed);
  FinalEffTightIsoCJdR_->SetLineColor(kRed);
  FinalEffDMFindCJdR_->SetLineColor(kRed);

  FinalEffLooseIsoCJPtGen_->SetLineColor(kRed);
  FinalEffMedIsoCJPtGen_->SetLineColor(kRed);
  FinalEffTightIsoCJPtGen_->SetLineColor(kRed);
  FinalEffDMFindCJPtGen_->SetLineColor(kRed);

  TLegend *leg = new TLegend(0.1,0.7,0.25,0.9);
  leg->AddEntry(FinalEffDMFindRECOPtGen_, "No CleanJets","L");
  leg->AddEntry(FinalEffDMFindCJPtGen_, "CleanJets","L");

  FinalEffLooseIsoSAMEPtCanvas.cd();
  FinalEffLooseIsoCJPt_->Draw();
  FinalEffLooseIsoRECOPt_->Draw("SAME");
  leg->Draw();

  FinalEffMedIsoSAMEPtCanvas.cd();
  FinalEffMedIsoCJPt_->Draw();
  FinalEffMedIsoRECOPt_->Draw("SAME");
  leg->Draw();

  FinalEffTightIsoSAMEPtCanvas.cd();
  FinalEffTightIsoCJPt_->Draw();
  FinalEffTightIsoRECOPt_->Draw("SAME");
  leg->Draw();

  FinalEffDMFindSAMEPtCanvas.cd();
  FinalEffDMFindCJPt_->Draw();
  FinalEffDMFindRECOPt_->Draw("SAME");
  leg->Draw();

  FinalEffLooseIsoSAMEdRCanvas.cd();
  FinalEffLooseIsoCJdR_->Draw();
  FinalEffLooseIsoRECOdR_->Draw("SAME");
  leg->Draw();

  FinalEffMedIsoSAMEdRCanvas.cd();
  FinalEffMedIsoCJdR_->Draw();
  FinalEffMedIsoRECOdR_->Draw("SAME");
  leg->Draw();

  FinalEffTightIsoSAMEdRCanvas.cd();
  FinalEffTightIsoCJdR_->Draw();
  FinalEffTightIsoRECOdR_->Draw("SAME");
  leg->Draw();

  FinalEffDMFindSAMEdRCanvas.cd();
  FinalEffDMFindCJdR_->Draw();
  FinalEffDMFindRECOdR_->Draw("SAME");
  leg->Draw();

  FinalEffLooseIsoSAMEPtGenCanvas.cd();
  FinalEffLooseIsoCJPtGen_->Draw();
  FinalEffLooseIsoRECOPtGen_->Draw("SAME");
  leg->Draw();

  FinalEffMedIsoSAMEPtGenCanvas.cd();
  FinalEffMedIsoCJPtGen_->Draw();
  FinalEffMedIsoRECOPtGen_->Draw("SAME");
  leg->Draw();

  FinalEffTightIsoSAMEPtGenCanvas.cd();
  FinalEffTightIsoCJPtGen_->Draw();
  FinalEffTightIsoRECOPtGen_->Draw("SAME");
  leg->Draw();

  FinalEffDMFindSAMEPtGenCanvas.cd();
  FinalEffDMFindCJPtGen_->Draw();
  FinalEffDMFindRECOPtGen_->Draw("SAME");
  leg->Draw();

std::cout << "<----------------Formatted Canvases and Histos-------------->" << std::endl;

  //Write output file
  out_->cd();

  NEventsCanvas.Write();
  NMuRemovedCanvas.Write();
  TauMuTauHaddRCanvas.Write();
  NConstituentsCJCanvas.Write();
  NTausRECOvsCLEANJETSCanvas.Write();
  GenDiTaudRvsCJDiTaudRCanvas.Write();
  TauHadnConstvsPtCanvas.Write();

  MatchedLooseIsoRECOPtCanvas.Write();
  MatchedMedIsoRECOPtCanvas.Write();
  MatchedTightIsoRECOPtCanvas.Write();
  MatchedDMFindRECOPtCanvas.Write();
  MatchedRECOPtCanvas.Write();
  FinalEffLooseIsoRECOPtCanvas.Write();
  FinalEffMedIsoRECOPtCanvas.Write();
  FinalEffTightIsoRECOPtCanvas.Write();
  FinalEffDMFindRECOPtCanvas.Write();

  MatchedLooseIsoCJPtCanvas.Write();
  MatchedMedIsoCJPtCanvas.Write();
  MatchedTightIsoCJPtCanvas.Write();
  MatchedDMFindCJPtCanvas.Write();
  MatchedCJPtCanvas.Write();
  FinalEffLooseIsoCJPtCanvas.Write();
  FinalEffMedIsoCJPtCanvas.Write();
  FinalEffTightIsoCJPtCanvas.Write();
  FinalEffDMFindCJPtCanvas.Write();


  MatchedLooseIsoRECOdRCanvas.Write();
  MatchedMedIsoRECOdRCanvas.Write();
  MatchedTightIsoRECOdRCanvas.Write();
  MatchedDMFindRECOdRCanvas.Write();
  MatchedRECOdRCanvas.Write();
  FinalEffLooseIsoRECOdRCanvas.Write();
  FinalEffMedIsoRECOdRCanvas.Write();
  FinalEffTightIsoRECOdRCanvas.Write();
  FinalEffDMFindRECOdRCanvas.Write();

  MatchedLooseIsoCJdRCanvas.Write();
  MatchedMedIsoCJdRCanvas.Write();
  MatchedTightIsoCJdRCanvas.Write();
  MatchedDMFindCJdRCanvas.Write();
  MatchedCJdRCanvas.Write();
  FinalEffLooseIsoCJdRCanvas.Write();
  FinalEffMedIsoCJdRCanvas.Write();
  FinalEffTightIsoCJdRCanvas.Write();
  FinalEffDMFindCJdRCanvas.Write();


  MatchedLooseIsoRECOPtGenCanvas.Write();
  MatchedMedIsoRECOPtGenCanvas.Write();
  MatchedTightIsoRECOPtGenCanvas.Write();
  MatchedDMFindRECOPtGenCanvas.Write();
  MatchedRECOPtGenCanvas.Write();
  FinalEffLooseIsoRECOPtGenCanvas.Write();
  FinalEffMedIsoRECOPtGenCanvas.Write();
  FinalEffTightIsoRECOPtGenCanvas.Write();
  FinalEffDMFindRECOPtGenCanvas.Write();

  MatchedLooseIsoCJPtGenCanvas.Write();
  MatchedMedIsoCJPtGenCanvas.Write();
  MatchedTightIsoCJPtGenCanvas.Write();
  MatchedDMFindCJPtGenCanvas.Write();
  MatchedCJPtGenCanvas.Write();
  FinalEffLooseIsoCJPtGenCanvas.Write();
  FinalEffMedIsoCJPtGenCanvas.Write();
  FinalEffTightIsoCJPtGenCanvas.Write();
  FinalEffDMFindCJPtGenCanvas.Write();

  FinalEffLooseIsoSAMEPtCanvas.Write();
  FinalEffMedIsoSAMEPtCanvas.Write();
  FinalEffTightIsoSAMEPtCanvas.Write();
  FinalEffDMFindSAMEPtCanvas.Write();
  FinalEffLooseIsoSAMEdRCanvas.Write();
  FinalEffMedIsoSAMEdRCanvas.Write();
  FinalEffTightIsoSAMEdRCanvas.Write();
  FinalEffDMFindSAMEdRCanvas.Write();
  FinalEffLooseIsoSAMEPtGenCanvas.Write();
  FinalEffMedIsoSAMEPtGenCanvas.Write();
  FinalEffTightIsoSAMEPtGenCanvas.Write();
  FinalEffDMFindSAMEPtGenCanvas.Write();

  out_->Write();
  out_->Close();
std::cout << "DONE" << std::endl;
}//EndJob

// ------------ method called when starting to processes a run  ------------
void GGHAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void GGHAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void GGHAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void GGHAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

//Delete Memory
void GGHAnalyzer::reset(const bool doDelete)
{
  if ((doDelete) && (NEvents_ != NULL)) delete NEvents_;
  NEvents_ = NULL;
  if ((doDelete) && (NMuRemoved_ != NULL)) delete NMuRemoved_;
  NMuRemoved_ = NULL;
  if ((doDelete) && (TauMuTauHaddR_ != NULL)) delete TauMuTauHaddR_;
  TauMuTauHaddR_ = NULL;
  if ((doDelete) && (NConstituentsCJ_ != NULL)) delete NConstituentsCJ_;
  NConstituentsCJ_ = NULL;
  if ((doDelete) && (NTausRECOvsCLEANJETS_ != NULL)) delete NTausRECOvsCLEANJETS_;
  NTausRECOvsCLEANJETS_ = NULL;
  if ((doDelete) && (GenDiTaudRvsCJDiTaudR_ != NULL)) delete GenDiTaudRvsCJDiTaudR_;
  GenDiTaudRvsCJDiTaudR_ = NULL;
  if ((doDelete) && (TauHadnConstvsPt_ != NULL)) delete TauHadnConstvsPt_;
  TauHadnConstvsPt_ = NULL;

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


}//void GGHAnalyzer

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GGHAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GGHAnalyzer);

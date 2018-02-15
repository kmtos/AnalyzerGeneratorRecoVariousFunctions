// -*- C++ -*-
//
// Package:    FakeRateMiniAODGetRates
// Class:      FakeRateMiniAODGetRates
// 
/**\class FakeRateMiniAODGetRates FakeRateMiniAODGetRates.cc Analyzer/src/FakeRateMiniAODGetRates.cc

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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

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

class FakeRateMiniAODGetRates : public edm::EDAnalyzer {
   public:
      typedef reco::JetFloatAssociation::Container JetBCEnergyRatioCollection;
      explicit FakeRateMiniAODGetRates(const edm::ParameterSet&);
      ~FakeRateMiniAODGetRates();

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
      edm::EDGetTokenT<edm::View<pat::Jet> > jetTag_;
      edm::EDGetTokenT<edm::View<pat::Muon> > muonsTag_;
      edm::EDGetTokenT<edm::View<pat::Tau> > tauTag_;
      std::string vLooseIsoTag_;
      std::string looseIsoTag_;
      std::string medIsoTag_;
      std::string tightIsoTag_;
      std::string vTightIsoTag_;
      std::string decayModeFindingTag_;
      std::string isoRawTag_;
      double mu3dRMin_;
      double mu3dRMax_;
      double tauPtCut_;
//      edm::EDGetTokenT<edm::View<pat::Muon> > mu3Tag_;
      edm::EDGetTokenT<edm::View<pat::Muon> > mu12Tag_;
      bool requireRemovedMuon_;
      bool checkInvMass_;
      double checkInvMassMin_;
      double checkInvMassMax_;


      //Histograms
      TH1F* NEvents_;   
      TH1F* InvMassTauMuMu1_;
      TH1F* InvMassMu1Mu2_;
      TH1F* InvMassTauMuMu2_;

      TH2F* EtavsPtTauLooseIso_;
      TH2F* EtavsPtTauMedIso_;
      TH2F* EtavsPtTauTightIso_;
      TH2F* EtavsPtTauDMFind_;
      TH2F* EtavsPtJet_;
      TH2F* EtavsPtJetSoftMuon_;
      TH2F* EtavsPtJetSoftMuon_noMu_;

      TH1F* TauLooseIsoEta_;
      TH1F* TauMedIsoEta_;
      TH1F* TauTightIsoEta_;
      TH1F* TauDMFindEta_;
      TH1F* JetEta_;
      TH1F* JetEtaWithSoftMuon_;
      TH1F* JetEtaWithSoftMuon_noMu_;
      TGraphAsymmErrors* FinalEffLooseIsoEta_;
      TGraphAsymmErrors* FinalEffMedIsoEta_;
      TGraphAsymmErrors* FinalEffTightIsoEta_;
      TGraphAsymmErrors* FinalEffDMFindEta_;

      TH1F* TauLooseIsoPt_;
      TH1F* TauMedIsoPt_;
      TH1F* TauTightIsoPt_;
      TH1F* TauDMFindPt_;
      TH1F* JetPt_;
      TH1F* JetPtWithSoftMuon_;
      TH1F* JetPtWithSoftMuon_noMu_;
      TGraphAsymmErrors* FinalEffLooseIsoPt_;
      TGraphAsymmErrors* FinalEffMedIsoPt_;
      TGraphAsymmErrors* FinalEffTightIsoPt_;
      TGraphAsymmErrors* FinalEffDMFindPt_;


     TH1F* OneProngDMEta_;
     TH1F* OneProngOnePizDMEta_;
     TH1F* OneProngTwoPizDMEta_;
     TH1F* ThreeProngDMEta_;
     TH1F* OneProngDMPt_;
     TH1F* OneProngOnePizDMPt_;
     TH1F* OneProngTwoPizDMPt_;
     TH1F* ThreeProngDMPt_;
     TGraphAsymmErrors* FinalOneProngDMEta_;
     TGraphAsymmErrors* FinalOneProngOnePizDMEta_;
     TGraphAsymmErrors* FinalOneProngTwoPizDMEta_;
     TGraphAsymmErrors* FinalThreeProngDMEta_;
     TGraphAsymmErrors* FinalOneProngDMPt_;
     TGraphAsymmErrors* FinalOneProngOnePizDMPt_;
     TGraphAsymmErrors* FinalOneProngTwoPizDMPt_;
     TGraphAsymmErrors* FinalThreeProngDMPt_;

     TH1F* OneProngMedIsoEta_;
     TH1F* OneProngOnePizMedIsoEta_;
     TH1F* OneProngTwoPizMedIsoEta_;
     TH1F* ThreeProngMedIsoEta_;
     TH1F* OneProngMedIsoPt_;
     TH1F* OneProngOnePizMedIsoPt_;
     TH1F* OneProngTwoPizMedIsoPt_;
     TH1F* ThreeProngMedIsoPt_;
     TGraphAsymmErrors* FinalOneProngMedIsoEta_;
     TGraphAsymmErrors* FinalOneProngOnePizMedIsoEta_;
     TGraphAsymmErrors* FinalOneProngTwoPizMedIsoEta_;
     TGraphAsymmErrors* FinalThreeProngMedIsoEta_;
     TGraphAsymmErrors* FinalOneProngMedIsoPt_;
     TGraphAsymmErrors* FinalOneProngOnePizMedIsoPt_;
     TGraphAsymmErrors* FinalOneProngTwoPizMedIsoPt_;
     TGraphAsymmErrors* FinalThreeProngMedIsoPt_;


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
FakeRateMiniAODGetRates::FakeRateMiniAODGetRates(const edm::ParameterSet& iConfig):
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  jetTag_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetTag"))),
  muonsTag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonsTag"))),
  tauTag_(consumes<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("tauTag"))),
  vLooseIsoTag_(iConfig.getParameter<std::string>("vLooseIsoTag")),
  looseIsoTag_(iConfig.getParameter<std::string>("looseIsoTag")),
  medIsoTag_(iConfig.getParameter<std::string>("medIsoTag")),
  tightIsoTag_(iConfig.getParameter<std::string>("tightIsoTag")),
  vTightIsoTag_(iConfig.getParameter<std::string>("vTightIsoTag")),
  decayModeFindingTag_(iConfig.getParameter<std::string>("decayModeFindingTag")),
  isoRawTag_(iConfig.getParameter<std::string>("isoRawTag")),
  mu3dRMin_(iConfig.getParameter<double>("mu3dRMin")),
  mu3dRMax_(iConfig.getParameter<double>("mu3dRMax")),
  tauPtCut_(iConfig.getParameter<double>("tauPtCut")),
//  mu3Tag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("mu3Tag"))),
  mu12Tag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("mu12Tag"))),
  requireRemovedMuon_(iConfig.getParameter<bool>("requireRemovedMuon")),
  checkInvMass_(iConfig.getParameter<bool>("checkInvMass")),
  checkInvMassMin_(iConfig.getParameter<double>("checkInvMassMin")),
  checkInvMassMax_(iConfig.getParameter<double>("checkInvMassMax"))
{
  reset(false);    
}//FakeRateMiniAODGetRates



FakeRateMiniAODGetRates::~FakeRateMiniAODGetRates()
{
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void FakeRateMiniAODGetRates::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "\n<------------THIS IS A NEW EVENT------------>" << std::endl;
  NEvents_->Fill(0);

  //Get ak4Jets particle collection
  edm::Handle<edm::View<pat::Jet> > pJets;
  iEvent.getByToken(jetTag_, pJets);

  //Get RECO Muons particle collection
  edm::Handle<edm::View<pat::Muon> > pMuons;
  iEvent.getByToken(muonsTag_, pMuons);

  //Get CleanJets Tau particle collection
  edm::Handle<edm::View<pat::Tau> > pTaus;
  iEvent.getByToken(tauTag_, pTaus);

//  //Old Jet collection for bTagging
//  edm::Handle<edm::View<pat::Muon> > pMu3;
//  iEvent.getByToken(mu3Tag_, pMu3);

  //Old Jet collection for bTagging
  edm::Handle<edm::View<pat::Muon> > pMu12;
  iEvent.getByToken(mu12Tag_, pMu12);

  pat::Muon mu1, mu2;
  double highestMuPt = -1;

  if (pMu12->size() > 2)
  {
    std::cout << "DAMMIT SOMETHIGN MESSED UP" << std::endl;
    return;
  }
  for(edm::View<pat::Muon>::const_iterator iMuon=pMu12->begin(); iMuon!=pMu12->end();++iMuon)
  {
    if (highestMuPt == -1)
    {
      highestMuPt = iMuon->pt();
      mu1 = *iMuon;
    }//if  
    else if (iMuon->pt() > highestMuPt)
    {
      mu2 = mu1;
      mu1 = *iMuon;
    }//else if
    else 
      mu2 = *iMuon;
  }

  reco::LeafCandidate::LorentzVector diMuP4 = mu1.p4() + mu2.p4();
std::cout << "mu1: pt=" << mu1.pt() << "  eta=" << mu1.eta() << "  pdgID=" << mu1.pdgId() << "\nmu2: pt=" << mu2.pt() << "  eta=" << mu2.eta() << "  pdgID=" << mu2.pdgId() << std::endl;
std::cout << "DiMu: pt= " << diMuP4.Pt() << "  eta=" << diMuP4.Eta() << " Mass=" << diMuP4.M() << std::endl;

  if (checkInvMass_ && (diMuP4.M() < checkInvMassMin_ || diMuP4.M() > checkInvMassMax_) )
    return;
  InvMassMu1Mu2_->Fill(diMuP4.M() );


//////////////////////////////
// Begin Analyzer
//////////////////////////////
  ///////////////////////////
  // iterating over the taus
  /////////////////////////// 
  reco::LeafCandidate::LorentzVector  diTauP4;
  double tauDecayMode = -1, DMFind = -1, LooseIso = -1, MedIso = -1, TightIso = -1;
  //double VLooseIso = -1, VTightIso = -1;
  for (edm::View<pat::Tau>::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); ++iTau)
  {
    NEvents_->Fill(1);
    double bestMu3dR = 10000000;
    bool checkSubMu = false;
    pat::Muon mu3;
    for (edm::View<pat::Muon>::const_iterator iMu = pMuons->begin(); iMu != pMuons->end(); ++iMu)
    {
      double currdR = deltaR(*iTau, *iMu), mu1dR = deltaR(mu1, *iMu), mu2dR = deltaR(mu2, *iMu);
      if (mu1dR < .5 || mu2dR < .5)
       continue;
      if (currdR < mu3dRMax_ && currdR > mu3dRMin_ && currdR < bestMu3dR)
      {
        diTauP4 = iTau->p4() + iMu->p4();
        checkSubMu = true;
        mu3 = *iMu;
        tauDecayMode = iTau->decayMode();
        DMFind = iTau->tauID(decayModeFindingTag_);
//        VLooseIso = iTau->tauID(vLooseIsoTag_);
        LooseIso = iTau->tauID(looseIsoTag_);
        MedIso = iTau->tauID(medIsoTag_);
        TightIso = iTau->tauID(tightIsoTag_);
//        VTightIso = iTau->tauID(vTightIsoTag_);
      }//if
    }//for iMu


    if (checkSubMu)
      NEvents_->Fill(2);
    else if (iTau->pt() > tauPtCut_ )
      NEvents_->Fill(3);
    else if (DMFind >= .5)
      NEvents_->Fill(4);
    else if (LooseIso >= .5)
      NEvents_->Fill(5);
    else if (MedIso >= .5)
      NEvents_->Fill(6);
    else if (TightIso >= .5)
      NEvents_->Fill(7);
 
    if ( (!checkSubMu && requireRemovedMuon_) || iTau->pt() < tauPtCut_ || fabs(iTau->eta() ) > 2.4)
      continue;

    if (checkSubMu)
    {
      reco::LeafCandidate::LorentzVector diMuP4_Mu1TauMu, diMuP4_Mu2TauMu;    
      diMuP4_Mu1TauMu = mu1.p4();
      diMuP4_Mu1TauMu += mu3.p4();
 
      diMuP4_Mu2TauMu = mu2.p4();
      diMuP4_Mu2TauMu += mu3.p4();

      InvMassTauMuMu1_->Fill(diMuP4_Mu1TauMu.M() );
      InvMassTauMuMu2_->Fill(diMuP4_Mu2TauMu.M() );
    }//if removed Mu

    if (DMFind >= .5)
    {
      TauDMFindPt_->Fill(iTau->pt() );
      TauDMFindEta_->Fill(iTau->eta() );
      EtavsPtTauDMFind_->Fill(iTau->pt(), iTau->eta() );
      if (tauDecayMode == 0)
      {
        OneProngDMEta_->Fill(iTau->eta() );
        OneProngDMPt_->Fill(iTau->pt() );  
      }//else if tauDecayMode >= .5
      else if (tauDecayMode == 1)
      {
        OneProngOnePizDMEta_->Fill(iTau->eta() );
        OneProngOnePizDMPt_->Fill(iTau->pt() );  
      }//else if tauDecayMode >= .5
      else if (tauDecayMode == 2)
      {
        OneProngTwoPizDMEta_->Fill(iTau->eta() );
        OneProngTwoPizDMPt_->Fill(iTau->pt() );  
      }//else if tauDecayMode >= .5
      else if (tauDecayMode == 10)
      {
        ThreeProngDMEta_->Fill(iTau->eta() );
        ThreeProngDMPt_->Fill(iTau->pt() );  
      }//else if tauDecayMode >= .5
    }//if DMFind >= .5

    if (TightIso >= .5 && DMFind >= .5)
    {
      TauTightIsoPt_->Fill(iTau->pt() );
      TauTightIsoEta_->Fill(iTau->eta() );
      EtavsPtTauTightIso_->Fill(iTau->pt(), iTau->eta() );
    }//if TightIso >= .5 &&DMFind >= .5

    if (MedIso >= .5 && DMFind >= .5)
    {
      TauMedIsoPt_->Fill(iTau->pt() );
      TauMedIsoEta_->Fill(iTau->eta() );
      EtavsPtTauMedIso_->Fill(iTau->pt(), iTau->eta() );
      if (tauDecayMode == 0)
      {
        OneProngMedIsoEta_->Fill(iTau->eta() );
        OneProngMedIsoPt_->Fill(iTau->pt() );  
      }//else if tauDecayMode >= .5
      else if (tauDecayMode == 1)
      {
        OneProngOnePizMedIsoEta_->Fill(iTau->eta() );
        OneProngOnePizMedIsoPt_->Fill(iTau->pt() );
      }//else if tauDecayMode >= .5
      else if (tauDecayMode == 2)
      {
        OneProngTwoPizMedIsoEta_->Fill(iTau->eta() );
        OneProngTwoPizMedIsoPt_->Fill(iTau->pt() );
      }//else if tauDecayMode >= .5
      else if (tauDecayMode == 10)
      {
        ThreeProngMedIsoEta_->Fill(iTau->eta() );
        ThreeProngMedIsoPt_->Fill(iTau->pt() );
      }//else if tauDecayMode == 1
    }//if MedIso == 1 && DMFind >= .5

    if (LooseIso >= .5 && DMFind >= .5)
    {
      TauLooseIsoPt_->Fill(iTau->pt() );
      TauLooseIsoEta_->Fill(iTau->eta() );
      EtavsPtTauLooseIso_->Fill(iTau->pt(), iTau->eta() );
    }//if Loose DMFind == 1
  }//iTau

  for (edm::View<pat::Jet>::const_iterator iJet = pJets->begin(); iJet != pJets->end(); ++iJet)
  {
    if (iJet->pt() > tauPtCut_ && fabs(iJet->eta() ) < 2.4 )
    {
      JetEta_->Fill(iJet->eta() ); 
      JetPt_->Fill(iJet->pt() ); 
      EtavsPtJet_->Fill(iJet->pt(), iJet->eta() );
      for (edm::View<pat::Muon>::const_iterator iMu = pMuons->begin(); iMu != pMuons->end(); ++iMu)
       {  
       double currdR = deltaR(*iJet, *iMu), mu1dR = deltaR(mu1, *iMu), mu2dR = deltaR(mu2, *iMu);
       if (mu1dR < .5 || mu2dR < .5)
         continue;
       if (currdR < mu3dRMax_ && currdR > mu3dRMin_)
        {
          JetEtaWithSoftMuon_->Fill(iJet->eta() );
          JetPtWithSoftMuon_->Fill(iJet->pt() );
          EtavsPtJetSoftMuon_->Fill(iJet->pt(), iJet->eta() );          
          reco::LeafCandidate::LorentzVector jetP4 = iJet->p4() - iMu->p4();
          JetEtaWithSoftMuon_noMu_->Fill(jetP4.Eta() );
          JetPtWithSoftMuon_noMu_->Fill(jetP4.Pt() );
          EtavsPtJetSoftMuon_noMu_->Fill(jetP4.Pt(), jetP4.Eta() );
        }//if iSoftMuon
      }//for i
    }//if 
  }//iJEt
}//End FakeRateMiniAODGetRates::analyze


// ------------ method called once each job just before starting event loop  ------------
void FakeRateMiniAODGetRates::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  //Open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");


  //Book histograms
  NEvents_     = new TH1F("NEvents"    , "", 9, -.5, 8.5);
      NEvents_->GetXaxis()->SetBinLabel(1, "TotalEvents"); 
      NEvents_->GetXaxis()->SetBinLabel(2, "Total #tau's");
      NEvents_->GetXaxis()->SetBinLabel(3, "#tau's With #mu");
      NEvents_->GetXaxis()->SetBinLabel(4, "#tau's p_{T} > 20");
      NEvents_->GetXaxis()->SetBinLabel(5, "Pass DMFind");
      NEvents_->GetXaxis()->SetBinLabel(6, "Pass Loose Iso");
      NEvents_->GetXaxis()->SetBinLabel(7, "Pass Med Iso");
      NEvents_->GetXaxis()->SetBinLabel(8, "Pass Tight Iso");
  InvMassTauMuMu1_     = new TH1F("InvMassTauMuMu1"    , "", 75, 0, 150);
  InvMassMu1Mu2_     = new TH1F("InvMassMu1Mu2"    , "", 75, 0, 150);
  InvMassTauMuMu2_     = new TH1F("InvMassTauMuMu2"    , "", 75, 0, 150);

  Float_t binsx[] = {0, 20, 25, 35, 60, 300};
  Float_t binsy[] = {0, .9, 1.5, 2.4};
  EtavsPtTauLooseIso_  = new TH2F("EtavsPtTauLooseIso" , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx, sizeof(binsy)/sizeof(Float_t) - 1, binsy);
  EtavsPtTauMedIso_  = new TH2F("EtavsPtTauMedIso"     , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx, sizeof(binsy)/sizeof(Float_t) - 1, binsy);
  EtavsPtTauTightIso_  = new TH2F("EtavsPtTauTightIso" , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx, sizeof(binsy)/sizeof(Float_t) - 1, binsy);
  EtavsPtTauDMFind_  = new TH2F("EtavsPtTauDMFind"     , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx, sizeof(binsy)/sizeof(Float_t) - 1, binsy);
  EtavsPtJet_  = new TH2F("EtavsPtJet"                 , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx, sizeof(binsy)/sizeof(Float_t) - 1, binsy);
  EtavsPtJetSoftMuon_  = new TH2F("EtavsPtJetSoftMuon" , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx, sizeof(binsy)/sizeof(Float_t) - 1, binsy);
  EtavsPtJetSoftMuon_noMu_  = new TH2F("EtavsPtJetSoftMuon_noMu" , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx, sizeof(binsy)/sizeof(Float_t) - 1, binsy);

  TauLooseIsoEta_  = new TH1F("TauLooseIsoEta"    , "", 11, -2.4, 2.4);
  TauMedIsoEta_    = new TH1F("TauMedIsoEta", "", 11, -2.4, 2.4);
  TauTightIsoEta_    = new TH1F("TauTightIsoEta", "", 11, -2.4, 2.4);
  TauDMFindEta_    = new TH1F("TauDMFindEta"    , "", 11, -2.4, 2.4);
  JetEta_          = new TH1F("JetEta"    , "", 11, -2.4, 2.4);
  JetEtaWithSoftMuon_          = new TH1F("JetEtaWithSoftMuon"    , "", 11, -2.4, 2.4);
  JetEtaWithSoftMuon_noMu_          = new TH1F("JetEtaWithSoftMuon_noMu"    , "", 11, -2.4, 2.4);

  TauLooseIsoPt_  = new TH1F("TauLooseIsoPt"    , "", 50, 20, 220.0);
  TauMedIsoPt_    = new TH1F("TauMedIsoPt", "", 50, 20, 220.0);
  TauTightIsoPt_    = new TH1F("TauTightIsoPt", "", 50, 20, 220.0);
  TauDMFindPt_    = new TH1F("TauDMFindPt"    , "", 50, 20, 220.0);
  JetPt_          = new TH1F("JetPt"    , "", 50, 20, 220.0);
  JetPtWithSoftMuon_          = new TH1F("JetPtWithSoftMuon"    , "", 50, 20, 220.0);
  JetPtWithSoftMuon_noMu_          = new TH1F("JetPtWithSoftMuon_noMu"    , "", 50, 20, 220.0);

  FinalEffLooseIsoEta_ = new TGraphAsymmErrors(11);
  FinalEffMedIsoEta_ = new TGraphAsymmErrors(11);
  FinalEffTightIsoEta_ = new TGraphAsymmErrors(11);
  FinalEffDMFindEta_ = new TGraphAsymmErrors(11);

  FinalEffLooseIsoPt_ = new TGraphAsymmErrors(11);
  FinalEffMedIsoPt_ = new TGraphAsymmErrors(11);
  FinalEffTightIsoPt_ = new TGraphAsymmErrors(11);
  FinalEffDMFindPt_ = new TGraphAsymmErrors(11);

  OneProngDMEta_ = new TH1F("OneProngDMEta"    , "", 11, -2.4, 2.4);
  OneProngOnePizDMEta_ = new TH1F("OneProngOnePizDMEta"    , "", 11, -2.4, 2.4);
  OneProngTwoPizDMEta_ = new TH1F("OneProngTwoPizDMEta"    , "", 11, -2.4, 2.4);
  ThreeProngDMEta_ = new TH1F("ThreeProngDMEta"    , "", 11, -2.4, 2.4);
  OneProngDMPt_ = new TH1F("OneProngDMPt"    , "", 11, 20, 220.0);
  OneProngOnePizDMPt_ = new TH1F("OneProngOnePizDMPt"    , "", 11, 20, 220.0);
  OneProngTwoPizDMPt_ = new TH1F("OneProngTwoPizDMPt"    , "", 11, 20, 220.0);
  ThreeProngDMPt_ = new TH1F("ThreeProngDMPt"    , "", 11, 20, 220.0);
  FinalOneProngDMEta_ = new TGraphAsymmErrors(11);
  FinalOneProngOnePizDMEta_ = new TGraphAsymmErrors(11);
  FinalOneProngTwoPizDMEta_ = new TGraphAsymmErrors(11);
  FinalThreeProngDMEta_ = new TGraphAsymmErrors(11);
  FinalOneProngDMPt_ = new TGraphAsymmErrors(11);
  FinalOneProngOnePizDMPt_ = new TGraphAsymmErrors(11);
  FinalOneProngTwoPizDMPt_ = new TGraphAsymmErrors(11);
  FinalThreeProngDMPt_ = new TGraphAsymmErrors(11);

  OneProngMedIsoEta_ = new TH1F("OneProngMedIsoEta"    , "", 11, -2.4, 2.4);
  OneProngOnePizMedIsoEta_ = new TH1F("OneProngOnePizMedIsoEta"    , "", 11, -2.4, 2.4);
  OneProngTwoPizMedIsoEta_ = new TH1F("OneProngTwoPizMedIsoEta"    , "", 11, -2.4, 2.4);
  ThreeProngMedIsoEta_ = new TH1F("ThreeProngMedIsoEta"    , "", 11, -2.4, 2.4);
  OneProngMedIsoPt_ = new TH1F("OneProngMedIsoPt"    , "", 11, 20, 220.0);
  OneProngOnePizMedIsoPt_ = new TH1F("OneProngOnePizMedIsoPt"    , "", 11, 20, 220.0);
  OneProngTwoPizMedIsoPt_ = new TH1F("OneProngTwoPizMedIsoPt"    , "", 11, 20, 220.0);
  ThreeProngMedIsoPt_ = new TH1F("ThreeProngMedIsoPt"    , "", 11, 20, 220.0);
  FinalOneProngMedIsoEta_ = new TGraphAsymmErrors(11);
  FinalOneProngOnePizMedIsoEta_ = new TGraphAsymmErrors(11);
  FinalOneProngTwoPizMedIsoEta_ = new TGraphAsymmErrors(11);
  FinalThreeProngMedIsoEta_ = new TGraphAsymmErrors(11);
  FinalOneProngMedIsoPt_ = new TGraphAsymmErrors(11);
  FinalOneProngOnePizMedIsoPt_ = new TGraphAsymmErrors(11);
  FinalOneProngTwoPizMedIsoPt_ = new TGraphAsymmErrors(11);
  FinalThreeProngMedIsoPt_ = new TGraphAsymmErrors(11);

  InvMassTauMuMu1_->Sumw2();
  InvMassMu1Mu2_->Sumw2();
  InvMassTauMuMu2_->Sumw2();

  EtavsPtTauLooseIso_->Sumw2();
  EtavsPtTauMedIso_->Sumw2();
  EtavsPtTauTightIso_->Sumw2();
  EtavsPtTauDMFind_->Sumw2();
  EtavsPtJet_->Sumw2();
  EtavsPtJetSoftMuon_->Sumw2();
  EtavsPtJetSoftMuon_noMu_->Sumw2();

  TauLooseIsoEta_->Sumw2();
  TauMedIsoEta_->Sumw2();
  TauTightIsoEta_->Sumw2();
  TauDMFindEta_->Sumw2();
  JetEta_->Sumw2();
  JetEtaWithSoftMuon_->Sumw2();
  JetEtaWithSoftMuon_noMu_->Sumw2();

  TauLooseIsoPt_->Sumw2();
  TauMedIsoPt_->Sumw2();
  TauTightIsoPt_->Sumw2();
  TauDMFindPt_->Sumw2();
  JetPt_->Sumw2();
  JetPtWithSoftMuon_->Sumw2();
  JetPtWithSoftMuon_noMu_->Sumw2();

  OneProngDMEta_->Sumw2();
  OneProngOnePizDMEta_->Sumw2();
  OneProngTwoPizDMEta_->Sumw2();
  ThreeProngDMEta_->Sumw2();
  OneProngDMPt_->Sumw2();
  OneProngOnePizDMPt_->Sumw2();
  OneProngTwoPizDMPt_->Sumw2();
  ThreeProngDMPt_->Sumw2();

  OneProngMedIsoEta_->Sumw2();
  OneProngOnePizMedIsoEta_->Sumw2();
  OneProngTwoPizMedIsoEta_->Sumw2();
  ThreeProngMedIsoEta_->Sumw2();
  OneProngMedIsoPt_->Sumw2();
  OneProngOnePizMedIsoPt_->Sumw2();
  OneProngTwoPizMedIsoPt_->Sumw2();
  ThreeProngMedIsoPt_->Sumw2();
}

// ------------ method called once each job just after ending the event loop  ------------
void FakeRateMiniAODGetRates::endJob()
{
  //Make the Canvases
  TCanvas NEventsCanvas("NEvents","",600,600);
  TCanvas InvMassTauMuMu1Canvas("InvMassTauMuMu1","",600,600);
  TCanvas InvMassMu1Mu2Canvas("InvMassMu1Mu2","",600,600);
  TCanvas InvMassTauMuMu2Canvas("InvMassTauMuMu2","",600,600);

  TCanvas EtavsPtTauLooseIsoCanvas("EtavsPtTauLooseIso","",600,600);
  TCanvas EtavsPtTauMedIsoCanvas("EtavsPtTauMedIso","",600,600);
  TCanvas EtavsPtTauTightIsoCanvas("EtavsPtTauTightIso","",600,600);
  TCanvas EtavsPtTauDMFindCanvas("EtavsPtTauDMFind","",600,600);
  TCanvas EtavsPtJetCanvas("EtavsPtJet","",600,600);
  TCanvas EtavsPtJetSoftMuonCanvas("EtavsPtJetSoftMuon","",600,600);
  TCanvas EtavsPtJetSoftMuon_noMuCanvas("EtavsPtJetSoftMuon_noMu","",600,600);

  TCanvas TauLooseIsoEtaCanvas("TauLooseIsoEta","",600,600);
  TCanvas TauMedIsoEtaCanvas("TauMedIsoEta","",600,600);
  TCanvas TauTightIsoEtaCanvas("TauTightIsoEta","",600,600);
  TCanvas TauDMFindEtaCanvas("TauDMFindEta","",600,600);
  TCanvas JetEtaCanvas("JetEta","",600,600);
  TCanvas JetEtaWithSoftMuonCanvas("JetEtaWithSoftMuon","",600,600);
  TCanvas JetEtaWithSoftMuon_noMuCanvas("JetEtaWithSoftMuon_noMu","",600,600);
  TCanvas FinalEffLooseIsoEtaCanvas("FinalEffLooseIsoEta","",600,600);
  TCanvas FinalEffMedIsoEtaCanvas("FinalEffMedIsoEta","",600,600);
  TCanvas FinalEffTightIsoEtaCanvas("FinalEffTightIsoEta","",600,600);
  TCanvas FinalEffDMFindEtaCanvas("FinalEffDMFindEta","",600,600);

  TCanvas TauLooseIsoPtCanvas("TauLooseIsoPt","",600,600);
  TCanvas TauMedIsoPtCanvas("TauMedIsoPt","",600,600);
  TCanvas TauTightIsoPtCanvas("TauTightIsoPt","",600,600);
  TCanvas TauDMFindPtCanvas("TauDMFindPt","",600,600);
  TCanvas JetPtCanvas("JetPt","",600,600);
  TCanvas JetPtWithSoftMuonCanvas("JetPtWithSoftMuon","",600,600);
  TCanvas JetPtWithSoftMuon_noMuCanvas("JetPtWithSoftMuon_noMu","",600,600);
  TCanvas FinalEffLooseIsoPtCanvas("FinalEffLooseIsoPt","",600,600);
  TCanvas FinalEffMedIsoPtCanvas("FinalEffMedIsoPt","",600,600);
  TCanvas FinalEffTightIsoPtCanvas("FinalEffTightIsoPt","",600,600);
  TCanvas FinalEffDMFindPtCanvas("FinalEffDMFindPt","",600,600);

  TCanvas OneProngDMEtaCanvas("OneProngDMEta","",600,600);
  TCanvas OneProngOnePizDMEtaCanvas("OneProngOnePizDMEta","",600,600);
  TCanvas OneProngTwoPizDMEtaCanvas("OneProngTwoPizDMEta","",600,600);
  TCanvas ThreeProngDMEtaCanvas("ThreeProngDMEta","",600,600);
  TCanvas OneProngDMPtCanvas("OneProngDMPt","",600,600);
  TCanvas OneProngOnePizDMPtCanvas("OneProngOnePizDMPt","",600,600);
  TCanvas OneProngTwoPizDMPtCanvas("OneProngTwoPizDMPt","",600,600);
  TCanvas ThreeProngDMPtCanvas("ThreeProngDMPt","",600,600);
  TCanvas FinalOneProngDMEtaCanvas("FinalOneProngDMEta","",600,600);
  TCanvas FinalOneProngOnePizDMEtaCanvas("FinalOneProngOnePizDMEta","",600,600);
  TCanvas FinalOneProngTwoPizDMEtaCanvas("FinalOneProngTwoPizDMEta","",600,600);
  TCanvas FinalThreeProngDMEtaCanvas("FinalThreeProngDMEta","",600,600);
  TCanvas FinalOneProngDMPtCanvas("FinalOneProngDMPt","",600,600);
  TCanvas FinalOneProngOnePizDMPtCanvas("FinalOneProngOnePizDMPt","",600,600);
  TCanvas FinalOneProngTwoPizDMPtCanvas("FinalOneProngTwoPizDMPt","",600,600);
  TCanvas FinalThreeProngDMPtCanvas("FinalThreeProngDMPt","",600,600);

  TCanvas OneProngMedIsoEtaCanvas("OneProngMedIsoEta","",600,600);
  TCanvas OneProngOnePizMedIsoEtaCanvas("OneProngOnePizMedIsoEta","",600,600);
  TCanvas OneProngTwoPizMedIsoEtaCanvas("OneProngTwoPizMedIsoEta","",600,600);
  TCanvas ThreeProngMedIsoEtaCanvas("ThreeProngMedIsoEta","",600,600);
  TCanvas OneProngMedIsoPtCanvas("OneProngMedIsoPt","",600,600);
  TCanvas OneProngOnePizMedIsoPtCanvas("OneProngOnePizMedIsoPt","",600,600);
  TCanvas OneProngTwoPizMedIsoPtCanvas("OneProngTwoPizMedIsoPt","",600,600);
  TCanvas ThreeProngMedIsoPtCanvas("ThreeProngMedIsoPt","",600,600);
  TCanvas FinalOneProngMedIsoEtaCanvas("FinalOneProngMedIsoEta","",600,600);
  TCanvas FinalOneProngOnePizMedIsoEtaCanvas("FinalOneProngOnePizMedIsoEta","",600,600);
  TCanvas FinalOneProngTwoPizMedIsoEtaCanvas("FinalOneProngTwoPizMedIsoEta","",600,600);
  TCanvas FinalThreeProngMedIsoEtaCanvas("FinalThreeProngMedIsoEta","",600,600);
  TCanvas FinalOneProngMedIsoPtCanvas("FinalOneProngMedIsoPt","",600,600);
  TCanvas FinalOneProngOnePizMedIsoPtCanvas("FinalOneProngOnePizMedIsoPt","",600,600);
  TCanvas FinalOneProngTwoPizMedIsoPtCanvas("FinalOneProngTwoPizMedIsoPt","",600,600);
  TCanvas FinalThreeProngMedIsoPtCanvas("FinalThreeProngMedIsoPt","",600,600);

std::cout << "<----------------Declared Canvases-------------->" << std::endl;

  //Format the 1D plots and Draw (canvas, hist, grid, log y, log z, color, size, style, xAxisTitle, xTitleSize, xLabelSize, xTitleOffSet, yAxisTitle, yTitleSize, yLabelSize, yTitleOffset)
  VariousFunctions::formatAndDrawCanvasAndHist1D(NEventsCanvas, NEvents_,
	 1, 0, 0, kBlack, 7, 20, "", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassTauMuMu1Canvas, InvMassTauMuMu1_,
	 1, 0, 0, kBlack, 7, 20, "Mass(#tau_{#mu} #mu_{1})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassMu1Mu2Canvas, InvMassMu1Mu2_,
	 1, 0, 0, kBlack, 7, 20, "Mass(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassTauMuMu2Canvas, InvMassTauMuMu2_,
	 1, 0, 0, kBlack, 7, 20, "Mass(#tau_{#mu} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);


  VariousFunctions::formatAndDrawCanvasAndHist2D(EtavsPtTauLooseIsoCanvas, EtavsPtTauLooseIso_,
         1, 0, 0, kBlack, 7, 20, "p_{T} #tau Loose Iso", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(EtavsPtTauMedIsoCanvas, EtavsPtTauMedIso_,
         1, 0, 0, kBlack, 7, 20, "p_{T} #tau Med Iso", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(EtavsPtTauTightIsoCanvas, EtavsPtTauTightIso_,
         1, 0, 0, kBlack, 7, 20, "p_{T} #tau Tight Iso", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(EtavsPtTauDMFindCanvas, EtavsPtTauDMFind_,
         1, 0, 0, kBlack, 7, 20, "p_{T} #tau DMFind", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(EtavsPtJetCanvas, EtavsPtJet_,
         1, 0, 0, kBlack, 7, 20, "p_{T} Jet", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(EtavsPtJetSoftMuonCanvas, EtavsPtJetSoftMuon_,
         1, 0, 0, kBlack, 7, 20, "p_{T}(Jet with Soft Mu)", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(EtavsPtJetSoftMuon_noMuCanvas, EtavsPtJetSoftMuon_noMu_,
         1, 0, 0, kBlack, 7, 20, "p_{T}(Jet with Soft Mu Removed)", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);

  VariousFunctions::formatAndDrawCanvasAndHist1D(TauLooseIsoEtaCanvas, TauLooseIsoEta_,
	 1, 0, 0, kBlack, 7, 20, "Eta(PFTau + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauMedIsoEtaCanvas, TauMedIsoEta_,
	 1, 0, 0, kBlack, 7, 20, "Eta(PFTau + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauTightIsoEtaCanvas, TauTightIsoEta_,
	 1, 0, 0, kBlack, 7, 20, "Eta(PFTau + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauDMFindEtaCanvas, TauDMFindEta_,
	 1, 0, 0, kBlack, 7, 20, "Eta(PFTau + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(JetEtaCanvas, JetEta_,
	 1, 0, 0, kBlack, 7, 20, "Eta(ak4Jets)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(JetEtaWithSoftMuonCanvas, JetEtaWithSoftMuon_,
	 1, 0, 0, kBlack, 7, 20, "Eta(ak4Jets) with Soft Muon", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(JetEtaWithSoftMuon_noMuCanvas, JetEtaWithSoftMuon_noMu_,
	 1, 0, 0, kBlack, 7, 20, "Eta(ak4Jets noMu) with Soft Muon", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoEtaCanvas, FinalEffLooseIsoEta_,
	 1, 1, 1, kBlack, 1, 20, "Eta(PFTau + Loose Iso + DM / nJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoEtaCanvas, FinalEffMedIsoEta_,
	 1, 1, 1, kBlack, 1, 20, "Eta(PFTau + Med Iso + DM / nJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoEtaCanvas, FinalEffTightIsoEta_,
	 1, 1, 1, kBlack, 1, 20, "Eta(PFTau + Tight Iso + DM / nJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffDMFindEtaCanvas, FinalEffDMFindEta_,
	 1, 1, 1, kBlack, 1, 20, "Eta(PFTau + DecayModeFinding / nJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);


  VariousFunctions::formatAndDrawCanvasAndHist1D(TauLooseIsoPtCanvas, TauLooseIsoPt_,
         1, 0, 0, kBlack, 7, 20, "Pt(PFTau + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauMedIsoPtCanvas, TauMedIsoPt_,
         1, 0, 0, kBlack, 7, 20, "Pt(PFTau + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauTightIsoPtCanvas, TauTightIsoPt_,
         1, 0, 0, kBlack, 7, 20, "Pt(PFTau + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauDMFindPtCanvas, TauDMFindPt_,
         1, 0, 0, kBlack, 7, 20, "Pt(PFTau + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(JetPtCanvas, JetPt_,
         1, 0, 0, kBlack, 7, 20, "Pt(ak4Jets)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(JetPtWithSoftMuonCanvas, JetPtWithSoftMuon_,
         1, 0, 0, kBlack, 7, 20, "Pt(ak4Jets) with Soft Muon", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(JetPtWithSoftMuon_noMuCanvas, JetPtWithSoftMuon_noMu_,
         1, 0, 0, kBlack, 7, 20, "Pt(ak4Jets noMu) with Soft Muon", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoPtCanvas, FinalEffLooseIsoPt_,
         1, 1, 1, kBlack, 1, 20, "Pt(PFTau + Loose Iso + DM / nJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoPtCanvas, FinalEffMedIsoPt_,
         1, 1, 1, kBlack, 1, 20, "Pt(PFTau + Med Iso + DM / nJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoPtCanvas, FinalEffTightIsoPt_,
         1, 1, 1, kBlack, 1, 20, "Pt(PFTau + Tight Iso + DM / nJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffDMFindPtCanvas, FinalEffDMFindPt_,
         1, 1, 1, kBlack, 1, 20, "Pt(PFTau + DecayModeFinding / nJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngDMEtaCanvas, OneProngDMEta_,
         1, 0, 0, kBlack, 7, 20, "1 Prong #eta(PFTau + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngOnePizDMEtaCanvas, OneProngOnePizDMEta_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 1 #pi^{0} #eta(PFTau + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngTwoPizDMEtaCanvas, OneProngTwoPizDMEta_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 2 #pi^{0} #eta(PFTau + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(ThreeProngDMEtaCanvas, ThreeProngDMEta_,
         1, 0, 0, kBlack, 7, 20, "3 Prong #eta(PFTau + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngDMPtCanvas, OneProngDMPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong p_{T}(PFTau + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngOnePizDMPtCanvas, OneProngOnePizDMPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 1 #pi^{0} p_{T}(PFTau + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngTwoPizDMPtCanvas, OneProngTwoPizDMPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 2 #pi^{0} p_{T}(PFTau + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(ThreeProngDMPtCanvas, ThreeProngDMPt_,
         1, 0, 0, kBlack, 7, 20, "3 Prong #eta(PFTau + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngDMEtaCanvas, FinalOneProngDMEta_,
         1, 1, 1, kBlack, 1, 20, "1 Prong #eta(PFTau + DM / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngOnePizDMEtaCanvas, FinalOneProngOnePizDMEta_,
         1, 1, 1, kBlack, 1, 20, "1 Prong + 1 #pi^{0} #eta(PFTau + DM / ak4PfJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngTwoPizDMEtaCanvas, FinalOneProngTwoPizDMEta_,
         1, 1, 1, kBlack, 1, 20, "1 Prong + 2 #pi^{0} #eta(PFTau + DM / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalThreeProngDMEtaCanvas, FinalThreeProngDMEta_,
         1, 1, 1, kBlack, 1, 20, "3 Prong #eta(PFTau + DM / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngDMPtCanvas, FinalOneProngDMPt_,
         1, 1, 1, kBlack, 1, 20, "1 Prong p_{T}(PFTau + DM / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngOnePizDMPtCanvas, FinalOneProngOnePizDMPt_,
         1, 1, 1, kBlack, 1, 20, "1 Prong + 1 #pi^{0} p_{T}(PFTau + DM / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngTwoPizDMPtCanvas, FinalOneProngTwoPizDMPt_,
         1, 1, 1, kBlack, 1, 20, "1 Prong + 2 #pi^{0} p_{T}(PFTau + DM / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalThreeProngDMPtCanvas, FinalThreeProngDMPt_,
         1, 1, 1, kBlack, 1, 20, "3 Prong p_{T}(PFTau + DM / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngMedIsoEtaCanvas, OneProngMedIsoEta_,
         1, 0, 0, kBlack, 7, 20, "1 Prong #eta(PFTau + DM + MedIso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngOnePizMedIsoEtaCanvas, OneProngOnePizMedIsoEta_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 1 #pi^{0} #eta(PFTau + DM + MedIso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngTwoPizMedIsoEtaCanvas, OneProngTwoPizMedIsoEta_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 2 #pi^{0} #eta(PFTau + DM + MedIso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(ThreeProngMedIsoEtaCanvas, ThreeProngMedIsoEta_,
         1, 0, 0, kBlack, 7, 20, "3 Prong #eta(PFTau + DM + MedIso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngMedIsoPtCanvas, OneProngMedIsoPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong p_{T}(PFTau + DM + MedIso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngOnePizMedIsoPtCanvas, OneProngOnePizMedIsoPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 1 #pi^{0} p_{T}(PFTau + DM + MedIso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngTwoPizMedIsoPtCanvas, OneProngTwoPizMedIsoPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 2 #pi^{0} p_{T}(PFTau + DM + MedIso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(ThreeProngMedIsoPtCanvas, ThreeProngMedIsoPt_,
         1, 0, 0, kBlack, 7, 20, "3 Prong #eta(PFTau + DM + MedIso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngMedIsoEtaCanvas, FinalOneProngMedIsoEta_,
         1, 1, 1, kBlack, 1, 20, "1 Prong #eta(PFTau + DM + MedIso / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngOnePizMedIsoEtaCanvas, FinalOneProngOnePizMedIsoEta_,
         1, 1, 1, kBlack, 1, 20, "1 Prong + 1 #pi^{0} #eta(PFTau + DM + MedIso / ak4PfJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngTwoPizMedIsoEtaCanvas, FinalOneProngTwoPizMedIsoEta_,
         1, 1, 1, kBlack, 1, 20, "1 Prong + 2 #pi^{0} #eta(PFTau + DM + MedIso / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalThreeProngMedIsoEtaCanvas, FinalThreeProngMedIsoEta_,
         1, 1, 1, kBlack, 1, 20, "3 Prong #eta(PFTau + DM + MedIso / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngMedIsoPtCanvas, FinalOneProngMedIsoPt_,
         1, 1, 1, kBlack, 1, 20, "1 Prong p_{T}(PFTau + DM + MedIso / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngOnePizMedIsoPtCanvas, FinalOneProngOnePizMedIsoPt_,
         1, 1, 1, kBlack, 1, 20, "1 Prong + 1 #pi^{0} p_{T}(PFTau + DM + MedIso / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngTwoPizMedIsoPtCanvas, FinalOneProngTwoPizMedIsoPt_,
         1, 1, 1, kBlack, 1, 20, "1 Prong + 2 #pi^{0} p_{T}(PFTau + DM + MedIso / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalThreeProngMedIsoPtCanvas, FinalThreeProngMedIsoPt_,
         1, 1, 1, kBlack, 1, 20, "3 Prong p_{T}(PFTau + DM + MedIso / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

std::cout << "after formatting" << std::endl;
  
////////////////////////
// For Pt Gen
////////////////////////
std::cout << "<----------------Formatted Canvases and Histos-------------->" << std::endl;

  //Write output file
  out_->cd();

  NEventsCanvas.Write();
  InvMassTauMuMu1Canvas.Write();
  InvMassMu1Mu2Canvas.Write();
  InvMassTauMuMu2Canvas.Write();
 
  EtavsPtTauLooseIsoCanvas.Write();
  EtavsPtTauMedIsoCanvas.Write();
  EtavsPtTauTightIsoCanvas.Write();
  EtavsPtTauDMFindCanvas.Write();
  EtavsPtJetCanvas.Write();
  EtavsPtJetSoftMuonCanvas.Write();
  EtavsPtJetSoftMuon_noMuCanvas.Write();

  TauLooseIsoEtaCanvas.Write();
  TauMedIsoEtaCanvas.Write();
  TauTightIsoEtaCanvas.Write();
  TauDMFindEtaCanvas.Write();
  JetEtaCanvas.Write();
  JetEtaWithSoftMuonCanvas.Write();
  JetEtaWithSoftMuon_noMuCanvas.Write();
  FinalEffLooseIsoEtaCanvas.Write();
  FinalEffMedIsoEtaCanvas.Write();
  FinalEffTightIsoEtaCanvas.Write();
  FinalEffDMFindEtaCanvas.Write();

  TauLooseIsoPtCanvas.Write();
  TauMedIsoPtCanvas.Write();
  TauTightIsoPtCanvas.Write();
  TauDMFindPtCanvas.Write();
  JetPtCanvas.Write();
  JetPtWithSoftMuonCanvas.Write();
  JetPtWithSoftMuon_noMuCanvas.Write();
  FinalEffLooseIsoPtCanvas.Write();
  FinalEffMedIsoPtCanvas.Write();
  FinalEffTightIsoPtCanvas.Write();
  FinalEffDMFindPtCanvas.Write();

  OneProngDMEtaCanvas.Write();
  OneProngOnePizDMEtaCanvas.Write();
  OneProngTwoPizDMEtaCanvas.Write();
  ThreeProngDMEtaCanvas.Write();
  OneProngDMPtCanvas.Write();
  OneProngOnePizDMPtCanvas.Write();
  OneProngTwoPizDMPtCanvas.Write();
  ThreeProngDMPtCanvas.Write();
  FinalOneProngDMEtaCanvas.Write();
  FinalOneProngOnePizDMEtaCanvas.Write();
  FinalOneProngTwoPizDMEtaCanvas.Write();
  FinalThreeProngDMEtaCanvas.Write();
  FinalOneProngDMPtCanvas.Write();
  FinalOneProngOnePizDMPtCanvas.Write();
  FinalOneProngTwoPizDMPtCanvas.Write();
  FinalThreeProngDMPtCanvas.Write();

  OneProngMedIsoEtaCanvas.Write();
  OneProngOnePizMedIsoEtaCanvas.Write();
  OneProngTwoPizMedIsoEtaCanvas.Write();
  ThreeProngMedIsoEtaCanvas.Write();
  OneProngMedIsoPtCanvas.Write();
  OneProngOnePizMedIsoPtCanvas.Write();
  OneProngTwoPizMedIsoPtCanvas.Write();
  ThreeProngMedIsoPtCanvas.Write();
  FinalOneProngMedIsoEtaCanvas.Write();
  FinalOneProngOnePizMedIsoEtaCanvas.Write();
  FinalOneProngTwoPizMedIsoEtaCanvas.Write();
  FinalThreeProngMedIsoEtaCanvas.Write();
  FinalOneProngMedIsoPtCanvas.Write();
  FinalOneProngOnePizMedIsoPtCanvas.Write();
  FinalOneProngTwoPizMedIsoPtCanvas.Write();
  FinalThreeProngMedIsoPtCanvas.Write();

  out_->Write();
  out_->Close();
std::cout << "DONE" << std::endl;
}//EndJob

// ------------ method called when starting to processes a run  ------------
void FakeRateMiniAODGetRates::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void FakeRateMiniAODGetRates::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void FakeRateMiniAODGetRates::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void FakeRateMiniAODGetRates::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

//Delete Memory
void FakeRateMiniAODGetRates::reset(const bool doDelete)
{
  if ((doDelete) && (NEvents_ != NULL)) delete NEvents_;
  NEvents_ = NULL;
  if ((doDelete) && (InvMassTauMuMu1_ != NULL)) delete InvMassTauMuMu1_;
  InvMassTauMuMu1_ = NULL;
  if ((doDelete) && (InvMassMu1Mu2_ != NULL)) delete InvMassMu1Mu2_;
  InvMassMu1Mu2_ = NULL;
  if ((doDelete) && (InvMassTauMuMu2_ != NULL)) delete InvMassTauMuMu2_;
  InvMassTauMuMu2_ = NULL;

  if ((doDelete) && (EtavsPtTauLooseIso_ != NULL)) delete EtavsPtTauLooseIso_;
  EtavsPtTauLooseIso_ = NULL;
  if ((doDelete) && (EtavsPtTauMedIso_ != NULL)) delete EtavsPtTauMedIso_;
  EtavsPtTauMedIso_ = NULL;
  if ((doDelete) && (EtavsPtTauTightIso_ != NULL)) delete EtavsPtTauTightIso_;
  EtavsPtTauTightIso_ = NULL;
  if ((doDelete) && (EtavsPtTauDMFind_ != NULL)) delete EtavsPtTauDMFind_;
  EtavsPtTauDMFind_ = NULL;
  if ((doDelete) && (EtavsPtJet_ != NULL)) delete EtavsPtJet_;
  EtavsPtJet_ = NULL;
  if ((doDelete) && (EtavsPtJetSoftMuon_ != NULL)) delete EtavsPtJetSoftMuon_;
  EtavsPtJetSoftMuon_ = NULL;
  if ((doDelete) && (EtavsPtJetSoftMuon_noMu_ != NULL)) delete EtavsPtJetSoftMuon_noMu_;
  EtavsPtJetSoftMuon_noMu_ = NULL;

  if ((doDelete) && (TauLooseIsoEta_ != NULL)) delete TauLooseIsoEta_;
  TauLooseIsoEta_ = NULL;
  if ((doDelete) && (TauMedIsoEta_ != NULL)) delete TauMedIsoEta_;
  TauMedIsoEta_ = NULL;
  if ((doDelete) && (TauTightIsoEta_ != NULL)) delete TauTightIsoEta_;
  TauTightIsoEta_ = NULL;
  if ((doDelete) && (TauDMFindEta_ != NULL)) delete TauDMFindEta_;
  TauDMFindEta_ = NULL;
  if ((doDelete) && (JetEta_ != NULL)) delete JetEta_;
  JetEta_ = NULL;
  if ((doDelete) && (JetEtaWithSoftMuon_ != NULL)) delete JetEtaWithSoftMuon_;
  JetEtaWithSoftMuon_ = NULL;
  if ((doDelete) && (JetEtaWithSoftMuon_noMu_ != NULL)) delete JetEtaWithSoftMuon_noMu_;
  JetEtaWithSoftMuon_noMu_ = NULL;
  if ((doDelete) && (FinalEffLooseIsoEta_ != NULL)) delete FinalEffLooseIsoEta_;
  FinalEffLooseIsoEta_ = NULL;
  if ((doDelete) && (FinalEffMedIsoEta_ != NULL)) delete FinalEffMedIsoEta_;
  FinalEffMedIsoEta_ = NULL;
  if ((doDelete) && (FinalEffTightIsoEta_ != NULL)) delete FinalEffTightIsoEta_;
  FinalEffTightIsoEta_ = NULL;
  if ((doDelete) && (FinalEffDMFindEta_ != NULL)) delete FinalEffDMFindEta_;
  FinalEffDMFindEta_ = NULL;

  if ((doDelete) && (TauLooseIsoPt_ != NULL)) delete TauLooseIsoPt_;
  TauLooseIsoPt_ = NULL;
  if ((doDelete) && (TauMedIsoPt_ != NULL)) delete TauMedIsoPt_;
  TauMedIsoPt_ = NULL;
  if ((doDelete) && (TauTightIsoPt_ != NULL)) delete TauTightIsoPt_;
  TauTightIsoPt_ = NULL;
  if ((doDelete) && (TauDMFindPt_ != NULL)) delete TauDMFindPt_;
  TauDMFindPt_ = NULL;
  if ((doDelete) && (JetPt_ != NULL)) delete JetPt_;
  JetPt_ = NULL;
  if ((doDelete) && (JetPtWithSoftMuon_ != NULL)) delete JetPtWithSoftMuon_;
  JetPtWithSoftMuon_ = NULL;
  if ((doDelete) && (JetPtWithSoftMuon_noMu_ != NULL)) delete JetPtWithSoftMuon_noMu_;
  JetPtWithSoftMuon_noMu_ = NULL;
  if ((doDelete) && (FinalEffLooseIsoPt_ != NULL)) delete FinalEffLooseIsoPt_;
  FinalEffLooseIsoPt_ = NULL;
  if ((doDelete) && (FinalEffMedIsoPt_ != NULL)) delete FinalEffMedIsoPt_;
  FinalEffMedIsoPt_ = NULL;
  if ((doDelete) && (FinalEffTightIsoPt_ != NULL)) delete FinalEffTightIsoPt_;
  FinalEffTightIsoPt_ = NULL;
  if ((doDelete) && (FinalEffDMFindPt_ != NULL)) delete FinalEffDMFindPt_;
  FinalEffDMFindPt_ = NULL;

  if ((doDelete) && (OneProngDMEta_ != NULL)) delete OneProngDMEta_;
  OneProngDMEta_ = NULL;
  if ((doDelete) && (OneProngOnePizDMEta_ != NULL)) delete OneProngOnePizDMEta_;
  OneProngOnePizDMEta_ = NULL;
  if ((doDelete) && (OneProngTwoPizDMEta_ != NULL)) delete OneProngTwoPizDMEta_;
  OneProngOnePizDMEta_ = NULL;
  if ((doDelete) && (ThreeProngDMEta_ != NULL)) delete ThreeProngDMEta_;
  ThreeProngDMEta_ = NULL;
  if ((doDelete) && (OneProngDMPt_ != NULL)) delete OneProngDMPt_;
  OneProngDMPt_ = NULL;
  if ((doDelete) && (OneProngOnePizDMPt_ != NULL)) delete OneProngOnePizDMPt_;
  OneProngOnePizDMPt_ = NULL;
  if ((doDelete) && (OneProngTwoPizDMPt_ != NULL)) delete OneProngTwoPizDMPt_;
  OneProngTwoPizDMPt_ = NULL;
  if ((doDelete) && (ThreeProngDMPt_ != NULL)) delete ThreeProngDMPt_;
  ThreeProngDMPt_ = NULL;
  if ((doDelete) && (FinalOneProngDMEta_ != NULL)) delete FinalOneProngDMEta_;
  FinalOneProngDMEta_ = NULL;
  if ((doDelete) && (FinalOneProngOnePizDMEta_ != NULL)) delete FinalOneProngOnePizDMEta_;
  FinalOneProngOnePizDMEta_ = NULL;
  if ((doDelete) && (FinalOneProngTwoPizDMEta_ != NULL)) delete FinalOneProngTwoPizDMEta_;
  FinalOneProngTwoPizDMEta_ = NULL;
  if ((doDelete) && (FinalThreeProngDMEta_ != NULL)) delete FinalThreeProngDMEta_;
  FinalThreeProngDMEta_ = NULL;
  if ((doDelete) && (FinalOneProngDMPt_ != NULL)) delete FinalOneProngDMPt_;
  FinalOneProngDMPt_ = NULL;
  if ((doDelete) && (FinalOneProngOnePizDMPt_ != NULL)) delete FinalOneProngOnePizDMPt_;
  FinalOneProngOnePizDMPt_ = NULL;
  if ((doDelete) && (FinalOneProngTwoPizDMPt_ != NULL)) delete FinalOneProngTwoPizDMPt_;
  FinalOneProngTwoPizDMPt_ = NULL;
  if ((doDelete) && (FinalThreeProngDMPt_ != NULL)) delete FinalThreeProngDMPt_;
  FinalThreeProngDMPt_ = NULL;

  if ((doDelete) && (OneProngMedIsoEta_ != NULL)) delete OneProngMedIsoEta_;
  OneProngMedIsoEta_ = NULL;
  if ((doDelete) && (OneProngOnePizMedIsoEta_ != NULL)) delete OneProngOnePizMedIsoEta_;
  OneProngOnePizMedIsoEta_ = NULL;
  if ((doDelete) && (OneProngTwoPizMedIsoEta_ != NULL)) delete OneProngTwoPizMedIsoEta_;
  OneProngOnePizMedIsoEta_ = NULL;
  if ((doDelete) && (ThreeProngMedIsoEta_ != NULL)) delete ThreeProngMedIsoEta_;
  ThreeProngMedIsoEta_ = NULL;
  if ((doDelete) && (OneProngMedIsoPt_ != NULL)) delete OneProngMedIsoPt_;
  OneProngMedIsoPt_ = NULL;
  if ((doDelete) && (OneProngOnePizMedIsoPt_ != NULL)) delete OneProngOnePizMedIsoPt_;
  OneProngOnePizMedIsoPt_ = NULL;
  if ((doDelete) && (OneProngTwoPizMedIsoPt_ != NULL)) delete OneProngTwoPizMedIsoPt_;
  OneProngTwoPizMedIsoPt_ = NULL;
  if ((doDelete) && (ThreeProngMedIsoPt_ != NULL)) delete ThreeProngMedIsoPt_;
  ThreeProngMedIsoPt_ = NULL;
  if ((doDelete) && (FinalOneProngMedIsoEta_ != NULL)) delete FinalOneProngMedIsoEta_;
  FinalOneProngMedIsoEta_ = NULL;
  if ((doDelete) && (FinalOneProngOnePizMedIsoEta_ != NULL)) delete FinalOneProngOnePizMedIsoEta_;
  FinalOneProngOnePizMedIsoEta_ = NULL;
  if ((doDelete) && (FinalOneProngTwoPizMedIsoEta_ != NULL)) delete FinalOneProngTwoPizMedIsoEta_;
  FinalOneProngTwoPizMedIsoEta_ = NULL;
  if ((doDelete) && (FinalThreeProngMedIsoEta_ != NULL)) delete FinalThreeProngMedIsoEta_;
  FinalThreeProngMedIsoEta_ = NULL;
  if ((doDelete) && (FinalOneProngMedIsoPt_ != NULL)) delete FinalOneProngMedIsoPt_;
  FinalOneProngMedIsoPt_ = NULL;
  if ((doDelete) && (FinalOneProngOnePizMedIsoPt_ != NULL)) delete FinalOneProngOnePizMedIsoPt_;
  FinalOneProngOnePizMedIsoPt_ = NULL;
  if ((doDelete) && (FinalOneProngTwoPizMedIsoPt_ != NULL)) delete FinalOneProngTwoPizMedIsoPt_;
  FinalOneProngTwoPizMedIsoPt_ = NULL;
  if ((doDelete) && (FinalThreeProngMedIsoPt_ != NULL)) delete FinalThreeProngMedIsoPt_;
  FinalThreeProngMedIsoPt_ = NULL;

}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void FakeRateMiniAODGetRates::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FakeRateMiniAODGetRates);

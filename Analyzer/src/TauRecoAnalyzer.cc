// -*- C++ -*-
//
// Package:    TauRecoAnalyzer
// Class:      TauRecoAnalyzer
// 
/**\class TauRecoAnalyzer TauRecoAnalyzer.cc Analyzer/src/TauRecoAnalyzer.cc

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


using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;


//
// class declaration
//

class TauRecoAnalyzer : public edm::EDAnalyzer {
   public:
      typedef reco::JetFloatAssociation::Container JetBCEnergyRatioCollection;
      explicit TauRecoAnalyzer(const edm::ParameterSet&);
      ~TauRecoAnalyzer();

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
      TH1F* TauMuHadMass_;
      TH2F* NTausRECOvsCLEANJETS_;
      TH2F* GenDiTaudRvsCJDiTaudR_;

      TH1F* TauMuPt_;
      TH1F* TauHadPt_;
      TH1F* MatchedTightIsoRECOPt_;
      TH1F* MatchedDMFindRECOPt_;
      TH1F* MatchedRECOPt_;
      TH1F* FinalEffLooseIsoRECOPt_;
      TH1F* FinalEffMedIsoRECOPt_;
      TH1F* FinalEffTightIsoRECOPt_;
      TH1F* FinalEffDMFindRECOPt_;

      TH1F* MatchedLooseIsoCJPt_;
      TH1F* MatchedMedIsoCJPt_;
      TH1F* MatchedTightIsoCJPt_;
      TH1F* MatchedDMFindCJPt_;
      TH1F* MatchedCJPt_;
      TH1F* FinalEffLooseIsoCJPt_;
      TH1F* FinalEffMedIsoCJPt_;
      TH1F* FinalEffTightIsoCJPt_;
      TH1F* FinalEffDMFindCJPt_;

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
TauRecoAnalyzer::TauRecoAnalyzer(const edm::ParameterSet& iConfig):
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
}//TauRecoAnalyzer



TauRecoAnalyzer::~TauRecoAnalyzer()
{
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void TauRecoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
/////////////////////// 
// Analyze
/////////////////////// 
  for (std::vector<reco::PFTau>::const_iterator iTauCJ = pTausCJ->begin(); iTauCJ != pTausCJ->end(); ++iTauCJ)
  {
    const reco::PFJetRef& tauJetRef = (*iTauCJ).jetRef();
    const reco::PFJetRef& tauRECOJetRef = (*pJetValMap)[tauJetRef];
    const reco::MuonRefVector& removedMuons = (*pMuonMap)[tauJetRef];
    NMuRemoved_->Fill(removedMuons.size() );
    if (removedMuons.size() > 1)
      NEvents_->Fill(2);

    //find the highest pT associated muon
    std::vector<reco::MuonRef> removedMuonRefs;
    for (reco::MuonRefVector::const_iterator iMuon = removedMuons.begin(); iMuon != removedMuons.end(); ++iMuon) {removedMuonRefs.push_back(*iMuon);}//for iMuon
    std::cout << "\tremovedMuonRefs.size= " << removedMuonRefs.size() << std::endl;
    for (unsigned int iMuon = 0; iMuon < removedMuonRefs.size(); iMuon++)
    {
      for (unsigned int jter = iMuon + 1; jter < removedMuonRefs.size(); jter++)
      {
        if (removedMuonRefs[jter]->pt() > removedMuonRefs[iMuon]->pt())
	{
 	  reco::MuonRef TEMPRef = removedMuonRefs[iMuon];
	  removedMuonRefs[iMuon] = removedMuonRefs[jter];  
	  removedMuonRefs[jter] = TEMPRef;
        }//if jter > iMuon
      }//for jter
    }//for iMuon

    for (unsigned int iMuon = 0; iMuon < removedMuonRefs.size(); iMuon++) {std::cout << "iMuon= " << iMuon << "  removedMuonRefs[iMuon]->pt()= " << removedMuonRefs[iMuon]->pt() << std::endl;}

    //Calculate Tau_mu tau_H dR
//    bool diTauCheck = false; 
    for (unsigned int iMuon = 0; iMuon < removedMuonRefs.size(); iMuon++)
    {
      double dPhi = reco::deltaPhi(removedMuonRefs[iMuon]->phi(), iTauCJ->phi() );
      double dR_tauMu = sqrt( (removedMuonRefs[iMuon]->eta() - iTauCJ->eta() ) * (removedMuonRefs[iMuon]->eta() - iTauCJ->eta() )  +  dPhi * dPhi );
      std::cout << "\t\t\tMuRef->pt()= " << removedMuonRefs[iMuon]->pt() << "  \tdR_tauMu= " << dR_tauMu << std::endl;
      if (dR_tauMu < .5)
      {
//	diTauCheck = true;
        NEvents_->Fill(1);
        TauMuTauHaddR_->Fill(dR_tauMu );
	TauMuPt_->Fill(removedMuonRefs[iMuon]->pt() );
        double muHadMass = (removedMuonRefs[iMuon]->p4() + iTauCJ->p4()).M();
	TauHadPt_->Fill(iTauCJ->pt() );
        TauMuHadMass_->Fill(muHadMass );
	break;
      }//if
    }//for iMuon


  }//for iTauCJ
}//End TauRecoAnalyzer::analyze


// ------------ method called once each job just before starting event loop  ------------
void TauRecoAnalyzer::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  //Open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //Book histograms
  NEvents_     = new TH1F("NEvents"    , "", 9, -.5, 8.5);
      NEvents_->GetXaxis()->SetBinLabel(1, "TotalEvents"); 
      NEvents_->GetXaxis()->SetBinLabel(2, "#tau_{#mu} + #tau_{had} Match");
      NEvents_->GetXaxis()->SetBinLabel(3, "##mu > 1");
      NEvents_->GetXaxis()->SetBinLabel(4, "");
      NEvents_->GetXaxis()->SetBinLabel(5, "");
      NEvents_->GetXaxis()->SetBinLabel(6, "");
      NEvents_->GetXaxis()->SetBinLabel(7, "");
  NMuRemoved_            = new TH1F("NMuRemoved"    , "", 9, -.5, 8.5);
  TauMuTauHaddR_       = new TH1F("TauMuTauHaddR"    , "", 100, 0, 100);
  TauMuHadMass_        = new TH1F("TauMuHadMass"    , "", 100, 0, 100);
  NTausRECOvsCLEANJETS_  = new TH2F("NTausRECOvsCLEANJETS" , "", 11, -.5, 10.5, 11, -.5, 10.5);
  GenDiTaudRvsCJDiTaudR_  = new TH2F("GenDiTaudRvsCJDiTaudR" , "", 50, 0, 10, 50, 0, 10);

  TauMuPt_  = new TH1F("TauMuPt"    , "", 30, 0, 90.0);
  TauHadPt_    = new TH1F("TauHadPt", "", 30, 0, 90);
  MatchedTightIsoRECOPt_    = new TH1F("MatchedTightIsoRECOPt", "", 30, 0, 90);
  MatchedDMFindRECOPt_    = new TH1F("MatchedDMFindRECOPt"    , "", 30, 0, 90);
  MatchedRECOPt_          = new TH1F("MatchedRECOPt"    , "", 30, 0, 90);
  FinalEffLooseIsoRECOPt_ = new TH1F("FinalEffLooseIsoRECOPt"    , "", 30, 0, 90.0);
      FinalEffLooseIsoRECOPt_->GetYaxis()->SetTitle("#epsilon");
  FinalEffMedIsoRECOPt_   = new TH1F("FinalEffMedIsoRECOPt", "", 30, 0, 90);
      FinalEffMedIsoRECOPt_->GetYaxis()->SetTitle("#epsilon");
  FinalEffTightIsoRECOPt_   = new TH1F("FinalEffTightIsoRECOPt", "", 30, 0, 90);
      FinalEffTightIsoRECOPt_->GetYaxis()->SetTitle("#epsilon");
  FinalEffDMFindRECOPt_   = new TH1F("FinalEffDMFindRECOPt"    , "", 30, 0, 90);
      FinalEffDMFindRECOPt_->GetYaxis()->SetTitle("#epsilon");

  MatchedLooseIsoCJPt_  = new TH1F("MatchedLooseIsoCJPt"    , "", 30, 0, 90.0);
  MatchedMedIsoCJPt_    = new TH1F("MatchedMedIsoCJPt", "", 30, 0, 90);
  MatchedTightIsoCJPt_    = new TH1F("MatchedTightIsoCJPt", "", 30, 0, 90);
  MatchedDMFindCJPt_    = new TH1F("MatchedDMFindCJPt"    , "", 30, 0, 90);
  MatchedCJPt_          = new TH1F("MatchedCJPt"    , "", 30, 0, 90);
  FinalEffLooseIsoCJPt_ = new TH1F("FinalEffLooseIsoCJPt"    , "", 30, 0, 90.0);
      FinalEffLooseIsoCJPt_->GetYaxis()->SetTitle("#epsilon");
  FinalEffMedIsoCJPt_   = new TH1F("FinalEffMedIsoCJPt", "", 30, 0, 90);
      FinalEffMedIsoCJPt_->GetYaxis()->SetTitle("#epsilon");
  FinalEffTightIsoCJPt_   = new TH1F("FinalEffTightIsoCJPt", "", 30, 0, 90);
      FinalEffTightIsoCJPt_->GetYaxis()->SetTitle("#epsilon");
  FinalEffDMFindCJPt_   = new TH1F("FinalEffDMFindCJPt"    , "", 30, 0, 90);
      FinalEffDMFindCJPt_->GetYaxis()->SetTitle("#epsilon");

}

// ------------ method called once each job just after ending the event loop  ------------
void TauRecoAnalyzer::endJob()
{
  //Make the Canvases
  TCanvas NEventsCanvas("NEvents","",600,600);
  TCanvas NMuRemovedCanvas("NMuRemoved","",600,600);
  TCanvas TauMuTauHaddRCanvas("TauMuTauHaddR","",600,600);
  TCanvas TauMuHadMassCanvas("TauMuHadMass","",600,600);
  TCanvas NTausRECOvsCLEANJETSCanvas("NTausRECOvsCLEANJETS","",600,600);
  TCanvas GenDiTaudRvsCJDiTaudRCanvas("GenDiTaudRvsCJDiTaudR","",600,600);

  TCanvas TauMuPtCanvas("TauMuPt","",600,600);
  TCanvas TauHadPtCanvas("TauHadPt","",600,600);
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


std::cout << "<----------------Declared Canvases-------------->" << std::endl;

  //Format the 1D plots and Draw (canvas, hist, grid, log y, log z, color, size, style, xAxisTitle, xTitleSize, xLabelSize, xTitleOffSet, yAxisTitle, yTitleSize, yLabelSize, yTitleOffset)
  VariousFunctions::formatAndDrawCanvasAndHist1D(NEventsCanvas, NEvents_, 1, 0, 0, kBlack, 7, 20, "", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NMuRemovedCanvas, NMuRemoved_, 1, 0, 0, kBlack, 7, 20, "N #mu Removed", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauMuTauHaddRCanvas, TauMuTauHaddR_, 1, 0, 0, kBlack, 7, 20, "#DeltaR(#tau_{mu} + #tau_{H})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauMuHadMassCanvas, TauMuHadMass_, 1, 0, 0, kBlack, 7, 20, "M(#tau_{#mu}#tau_{H})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist2D(NTausRECOvsCLEANJETSCanvas, NTausRECOvsCLEANJETS_, 1, 0, 0, kBlack, 7, 20, "nRECO #tau's", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(GenDiTaudRvsCJDiTaudRCanvas, GenDiTaudRvsCJDiTaudR_, 1, 0, 0, kBlack, 7, 20, "gen #DeltaR(#tau#tau)", .04, .04, 1.1, "CJ #DeltaR(#tau#tau)", .04, .04, 1.6, "", .04, .04, 1.0);

  VariousFunctions::formatAndDrawCanvasAndHist1D(TauMuPtCanvas, TauMuPt_, 1, 0, 0, kBlack, 7, 20, "Pt(#mu from found di-#tau)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauHadPtCanvas, TauHadPt_, 1, 0, 0, kBlack, 7, 20, "Pt(#tau_{H} from found di-#tau)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoRECOPtCanvas, MatchedTightIsoRECOPt_, 1, 0, 0, kBlack, 7, 20, "Pt(RECO Matched + Tight Iso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedDMFindRECOPtCanvas, MatchedDMFindRECOPt_, 1, 0, 0, kBlack, 7, 20, "Pt(RECO Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedRECOPtCanvas, MatchedRECOPt_, 1, 0, 0, kBlack, 7, 20, "Pt(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(FinalEffLooseIsoRECOPtCanvas, FinalEffLooseIsoRECOPt_, 1, 0, 0, kBlack, 7, 20, "Pt(RECO Matched + Loose Iso) / Pt(RECO Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(FinalEffMedIsoRECOPtCanvas, FinalEffMedIsoRECOPt_, 1, 0, 0, kBlack, 7, 20, "Pt(RECO Matched + Med Iso) / Pt(RECO Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(FinalEffTightIsoRECOPtCanvas, FinalEffTightIsoRECOPt_, 1, 0, 0, kBlack, 7, 20, "Pt(RECO Matched + Tight Iso) / Pt(RECO Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(FinalEffDMFindRECOPtCanvas, FinalEffDMFindRECOPt_, 1, 0, 0, kBlack, 7, 20, "Pt(RECO Matched + DecayModeFinding) / Pt(RECO Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoCJPtCanvas, MatchedLooseIsoCJPt_, 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + Loose Iso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoCJPtCanvas, MatchedMedIsoCJPt_, 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + Med Iso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoCJPtCanvas, MatchedTightIsoCJPt_, 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + Tight Iso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedDMFindCJPtCanvas, MatchedDMFindCJPt_, 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedCJPtCanvas, MatchedCJPt_, 1, 0, 0, kBlack, 7, 20, "Pt(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(FinalEffLooseIsoCJPtCanvas, FinalEffLooseIsoCJPt_, 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + Loose Iso) / Pt(CLeanJets Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(FinalEffMedIsoCJPtCanvas, FinalEffMedIsoCJPt_, 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + Med Iso) / Pt(CLeanJets Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(FinalEffTightIsoCJPtCanvas, FinalEffTightIsoCJPt_, 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + Tight Iso) / Pt(CLeanJets Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(FinalEffDMFindCJPtCanvas, FinalEffDMFindCJPt_, 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + DecayModeFinding) / Pt(CLeanJets Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);


std::cout << "<----------------Formatted Canvases and Histos-------------->" << std::endl;

  //Write output file
  out_->cd();

  NEventsCanvas.Write();
  NMuRemovedCanvas.Write();
  TauMuTauHaddRCanvas.Write();
  TauMuHadMassCanvas.Write();
  NTausRECOvsCLEANJETSCanvas.Write();
  GenDiTaudRvsCJDiTaudRCanvas.Write();

  TauMuPtCanvas.Write();
  TauHadPtCanvas.Write();
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

  out_->Write();
  out_->Close();
std::cout << "DONE" << std::endl;
}//EndJob

// ------------ method called when starting to processes a run  ------------
void TauRecoAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void TauRecoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void TauRecoAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void TauRecoAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

//Delete Memory
void TauRecoAnalyzer::reset(const bool doDelete)
{
  if ((doDelete) && (NEvents_ != NULL)) delete NEvents_;
  NEvents_ = NULL;
  if ((doDelete) && (NMuRemoved_ != NULL)) delete NMuRemoved_;
  NMuRemoved_ = NULL;
  if ((doDelete) && (TauMuTauHaddR_ != NULL)) delete TauMuTauHaddR_;
  TauMuTauHaddR_ = NULL;
  if ((doDelete) && (TauMuHadMass_ != NULL)) delete TauMuHadMass_;
  TauMuHadMass_ = NULL;
  if ((doDelete) && (NTausRECOvsCLEANJETS_ != NULL)) delete NTausRECOvsCLEANJETS_;
  NTausRECOvsCLEANJETS_ = NULL;
  if ((doDelete) && (GenDiTaudRvsCJDiTaudR_ != NULL)) delete GenDiTaudRvsCJDiTaudR_;
  GenDiTaudRvsCJDiTaudR_ = NULL;

  if ((doDelete) && (TauMuPt_ != NULL)) delete TauMuPt_;
  TauMuPt_ = NULL;
  if ((doDelete) && (TauHadPt_ != NULL)) delete TauHadPt_;
  TauHadPt_ = NULL;
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

}//void TauRecoAnalyzer

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TauRecoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauRecoAnalyzer);

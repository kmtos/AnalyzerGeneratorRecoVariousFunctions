// -*- C++ -*-
//
// Package:    TTbarAnalyzer
// Class:      TTbarAnalyzer
// 
/**\class TTbarAnalyzer TTbarAnalyzer.cc Analyzer/src/TTbarAnalyzer.cc

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
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"

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

class TTbarAnalyzer : public edm::EDAnalyzer {
   public:
      typedef reco::JetFloatAssociation::Container JetBCEnergyRatioCollection;
      explicit TTbarAnalyzer(const edm::ParameterSet&);
      ~TTbarAnalyzer();

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

      //openning the jet output file
      ofstream jetOutput_;

      //name of output root file
      std::string outFileName_;

      //name of output jetOutput file
      std::string jetOutputFileName_;

      //gen particle tag
      edm::InputTag prunedGenParticleTag_;

      //PAT jet tag
      edm::InputTag slimmedJetTag_;

      //PAT TAU tag
      edm::InputTag slimmedTauTag_;

      //Histograms
      TH1F* NEvents_;   
      TH1F* GENDiTauDR_;
      TH1F* nGenTauPart_;
      TH1F* nGenBPart_;
      TH1F* GenTauMomPDGID_;
      TH1F* GenBMomPDGID_;   
      TH2F* GENTauMomPDGIDvsPt_;
      TH2F* GENBMomPDGIDvsPt_;

      TH1F* TauPt_;
      TH1F* TauPtDiff_;
      TH1F* TauDR_;
      TH1F* DiTauDR_;
      TH1F* BMatchFlavor_;

      TH1F* JetPt_;
      TH1F* JetPtDiff_;
      TH1F* JetDR_;

      TH1F* TauIsoDisc_;
      TH1F* TauMediumIsoEff_;
 
      TH2F* TauIsovsPt_;
      TH2F* DiTaudRvsBPt_;
      TH2F* TauNearestBdRvsBPt_;
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
TTbarAnalyzer::TTbarAnalyzer(const edm::ParameterSet& iConfig):
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  jetOutputFileName_(iConfig.getParameter<std::string>("jetOutputFileName")), 
  prunedGenParticleTag_(iConfig.getParameter<edm::InputTag>("prunedGenParticleTag")),
  slimmedJetTag_(iConfig.getParameter<edm::InputTag>("slimmedJetTag")),
  slimmedTauTag_(iConfig.getParameter<edm::InputTag>("slimmedTauTag"))
{
  reset(false);    
}//TTbarAnalyzer



TTbarAnalyzer::~TTbarAnalyzer()
{
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void TTbarAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "\n<------------THIS IS A NEW EVENT------------>" << std::endl;
  NEvents_->Fill(-1);

  // define a jet handle
  edm::Handle<std::vector<pat::Jet> > pSlimmedJetTag;
  iEvent.getByLabel(slimmedJetTag_, pSlimmedJetTag);

  // define a jet handle
  edm::Handle<std::vector<pat::Tau> > pSlimmedTauTag;
  iEvent.getByLabel(slimmedTauTag_, pSlimmedTauTag);

  //Get gen particle collection
  edm::Handle<reco::GenParticleCollection> pPrunGenParts;
  iEvent.getByLabel(prunedGenParticleTag_, pPrunGenParts);


//CHECK THAT THERE ARE AT LEAST 2 taus IN SLIMMED TAUS 
  std::cout << "\tSize of pSlimmedTauTag= " << pSlimmedTauTag->size() << std::endl;
  int slimTauSize = pSlimmedTauTag->size(); 
  NEvents_->Fill(slimTauSize);
  if (slimTauSize < 2)
    return;

//  
//GEN PARTICLES LOOP
//
  reco::GenParticleCollection::const_iterator firstGen = pPrunGenParts->begin();
  reco::GenParticleRef b1Ref = firstGen->daughterRef(0), b2Ref = firstGen->daughterRef(0), tau1Ref = firstGen->daughterRef(0), tau2Ref = firstGen->daughterRef(0);
  int nTau = 0, nB = 0;
  for (reco::GenParticleCollection::const_iterator iGenParticle = pPrunGenParts->begin(); iGenParticle != pPrunGenParts->end(); ++iGenParticle)
  {
    if (abs(iGenParticle->pdgId() ) == 15)
    {
      std::cout << "\t\t\ttau Particle pdgId= " << iGenParticle->pdgId() << "\t\t\t\tmother pdgId= " << iGenParticle->mother(0)->pdgId() << std::endl;
      if (abs(iGenParticle->mother(0)->pdgId() ) != 15)
        nTau++;
      GenTauMomPDGID_->Fill(iGenParticle->mother(0)->pdgId() );
      GENTauMomPDGIDvsPt_->Fill(iGenParticle->mother(0)->pdgId(), iGenParticle->pt() );
    }//if
    if (abs(iGenParticle->pdgId() ) == 5)
    {
      std::cout << "\t\t\tb Particle pdgId= " << iGenParticle->pdgId() << "\t\t\t\tmother pdgId= " << iGenParticle->mother(0)->pdgId() << std::endl;
      if (abs(iGenParticle->mother(0)->pdgId() ) != 5)
        nB++;
      GenBMomPDGID_->Fill(iGenParticle->mother(0)->pdgId() );
      GENBMomPDGIDvsPt_->Fill(iGenParticle->mother(0)->pdgId(), iGenParticle->pt() );
    }//if
  }//for Gen Particle 

std::cout << "\t\t\t\tnB= " << nB << "\tnTau= " << nTau << std::endl;
nGenTauPart_->Fill(nTau);
nGenBPart_->Fill(nB);
//
//END GEN PARTICLES LOOP
//

 
//
//START GEN MATCHING SLIMMED B JETS
//
  double b1pt = -1, b2pt = -1, b1eta = -10, b2eta = -10, b1phi = -10, b2phi = -10;
  double CSV1 = -1, CSV2 = -1, partFlavor1 = -1, partFlavor2 = -1;
  for (auto jet = pSlimmedJetTag->begin(); jet != pSlimmedJetTag->end(); ++jet)
  {
std::cout << "\t\tbjet1:  CSV= " << CSV1 << "  Pt= " << b1pt << std::endl;
std::cout << "\t\tbjet2:  CSV= " << CSV2 << "  Pt= " << b2pt << std::endl;
    if (jet->pt() > b1pt )
    {
      b2pt = b1pt;
      b2eta = b1eta;
      b2phi = b1phi;
      CSV2 = CSV1;
      partFlavor2 = partFlavor1;
      b1pt = jet->pt();
      b1eta = jet->eta();
      b1phi = jet->phi();
      CSV1 = jet->bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
      partFlavor1 = jet->partonFlavour();
    }//if jet->pt > b1pt
    if (jet->pt() > b2pt && jet->pt() < b1pt)
    {
      b2pt = jet->pt();
      b2eta = jet->eta();
      b2phi = jet->phi();
      CSV2 = jet->bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
      partFlavor2 = jet->partonFlavour();
    }//if jet->pt() > b2pt
  }//for auto jet


std::cout << "\tbjet1:  CSV= " << CSV1 << "  Pt= " << b1pt << "  Eta= " << b1eta << "  Phi= " << b1phi << "  Flavor= " << partFlavor1 << std::endl;
std::cout << "\tbjet2:  CSV= " << CSV2 << "  Pt= " << b2pt << "  Eta= " << b2eta << "  Phi= " << b2phi << "  Flavor= " << partFlavor2 << std::endl;

//
//END GEN MATCHING SLIMMED B JETS
//
 
//
//START SLIMMED TAUS GEN MATCHING
//
   double tau1Pt = -1, tau2Pt = -1, tau1Eta = -10, tau2Eta = -10, tau1Phi = -10, tau2Phi = -10;
   int loose1 = 0, loose2 = 0, medium1 = 0, medium2 = 0, tight1 = 0, tight2 = 0;
  for (auto iTau = pSlimmedTauTag->begin(); iTau != pSlimmedTauTag->end(); ++iTau)
  {
std::cout << "\tTau1: Pt= " << tau1Pt << "  Eta= " << tau1Eta << "  Phi= " << tau1Phi << std::endl;
std::cout << "\tTau2: Pt= " << tau2Pt << "  Eta= " << tau2Eta << "  Phi= " << tau2Phi << std::endl;
    if (iTau->pt() > tau1Pt )
    {
      tau2Pt = tau1Pt;
      tau2Eta = tau1Eta;
      tau2Phi = tau1Phi;
      loose2 = loose1;
      medium2 = medium2;
      tight2 = tight2;
      tau1Pt = iTau->pt();
      tau1Eta = iTau->eta();
      tau1Phi = iTau->phi();
      loose1 = iTau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"); 
      medium1 = iTau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
      tight1 = iTau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
      if (iTau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits") == 1)
        TauMediumIsoEff_->Fill(1);
    }//if iTau->pt > tau1Pt 
    if (iTau->pt() > tau2Pt && iTau->pt() < tau1Pt)
    { 
      tau2Pt = iTau->pt();
      tau2Eta = iTau->eta();
      tau2Phi = iTau->phi();
      loose2 = iTau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"); 
      medium2 = iTau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
      tight2 = iTau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
      if (iTau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits") == 1)
        TauMediumIsoEff_->Fill(1);
    }//if iTau->pt() > tau2Pt
  }//for auto iTau

std::cout << "\tTau1: Pt= " << tau1Pt << "  Eta= " << tau1Eta << "  Phi= " << tau1Phi << std::endl;
std::cout << "\tTau2: Pt= " << tau2Pt << "  Eta= " << tau2Eta << "  Phi= " << tau2Phi << std::endl;

//
//END SLIMMED TAUS GEN MATCHING
//




//
//START FILLING HISTOGRAMS
//
  TauPt_->Fill(tau1Pt);
  TauPt_->Fill(tau2Pt);
  DiTauDR_->Fill( sqrt( (tau2Eta - tau1Eta) * (tau2Eta - tau1Eta)  +  (tau2Phi - tau1Phi) * (tau2Phi - tau1Phi) ) );
  JetPt_->Fill(b1pt);
  JetPt_->Fill(b2pt);
  BMatchFlavor_->Fill(partFlavor1);
  BMatchFlavor_->Fill(partFlavor2);
 
  TauIsoDisc_->Fill(0);
  TauIsovsPt_->Fill(0.01, tau1Pt);
  if (loose1 == 1)
  {
    TauIsoDisc_->Fill(1);
    TauIsovsPt_->Fill(1, tau1Pt);
  }//if loose 1
  if (medium1 == 1)
  {
    TauIsoDisc_->Fill(2);
    TauIsovsPt_->Fill(2, tau1Pt);
  }//if medium1
  if (tight1 == 1)
  {
    TauIsoDisc_->Fill(3);
    TauIsovsPt_->Fill(3, tau1Pt);
  }//if tight 1
  TauIsoDisc_->Fill(0);
  TauIsovsPt_->Fill(0.01, tau2Pt);
  if (loose2 == 1)
  {
    TauIsoDisc_->Fill(1);
    TauIsovsPt_->Fill(1, tau2Pt);
  }//if loose 2
  if (medium2 == 1)
  {
    TauIsoDisc_->Fill(2);
    TauIsovsPt_->Fill(2, tau2Pt);
  }//fi medium2
  if (tight2 == 1)
  {
    TauIsoDisc_->Fill(3);
    TauIsovsPt_->Fill(3, tau2Pt);
  }//if tight 2
  double dRtau1b1 = sqrt( (b1eta - tau1Eta) * (b1eta - tau1Eta)  +  (b1phi - tau1Phi) * (b1phi - tau1Phi) ), dRtau1b2 = sqrt( (b2eta - tau1Eta) * (b2eta - tau1Eta)  +  (b2phi - tau1Phi) * (b2phi - tau1Phi) );
  double dRtau2b1 = sqrt( (b1eta - tau2Eta) * (b1eta - tau2Eta)  +  (b1phi - tau2Phi) * (b1phi - tau2Phi) ), dRtau2b2 = sqrt( (b2eta - tau2Eta) * (b2eta - tau2Eta)  +  (b2phi - tau2Phi) * (b2phi - tau2Phi) );
  if (dRtau1b1 < dRtau1b2 && dRtau1b1 < dRtau2b1 && dRtau1b1 < dRtau2b2)
    TauNearestBdRvsBPt_->Fill(dRtau1b1, b1pt);
  if (dRtau2b1 < dRtau1b1 && dRtau2b1 < dRtau1b2 && dRtau2b1 < dRtau2b2)
    TauNearestBdRvsBPt_->Fill(dRtau2b1, b1pt);
  if (dRtau1b2 < dRtau1b1 && dRtau1b2 < dRtau2b1 && dRtau1b2 < dRtau2b2)
    TauNearestBdRvsBPt_->Fill(dRtau1b2, b2pt);
  if (dRtau2b2 < dRtau1b1 && dRtau2b2 < dRtau1b2 && dRtau2b2 < dRtau2b1)
    TauNearestBdRvsBPt_->Fill(dRtau2b2, b2pt);

  if (b1pt > b2pt)
    DiTaudRvsBPt_->Fill( sqrt( (tau2Eta - tau1Eta) * (tau2Eta - tau1Eta)  +  (tau2Phi - tau1Phi) * (tau2Phi - tau1Phi) ),  b1pt);
  if (b1pt < b2pt)
    DiTaudRvsBPt_->Fill( sqrt( (tau2Eta - tau1Eta) * (tau2Eta - tau1Eta)  +  (tau2Phi - tau1Phi) * (tau2Phi - tau1Phi) ),  b2pt);
 
//
//END FILLING HISTOGRAMS
//
}//End TTbarAnalyzer::analyze


// ------------ method called once each job just before starting event loop  ------------
void TTbarAnalyzer::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  //Open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");
  jetOutput_.open (jetOutputFileName_.c_str());

  //Book histograms
  NEvents_     = new TH1F("NEvents"    , "", 9, -1.5, 8.5);
  GENDiTauDR_     = new TH1F("GENDiTauDR"    , "", 100, 0, 8);
  nGenTauPart_     = new TH1F("nGenTauPart"    , "", 21, -.5, 20.5);
  nGenBPart_     = new TH1F("nGenBPart"    , "", 26, -.5, 25.5);
  GenTauMomPDGID_     = new TH1F("GenTauMomPDGID"    , "", 1001, -.5, 1000.5);
  GenBMomPDGID_     = new TH1F("GenBMomPDGID"    , "", 1001, -.5, 1000.5);
  GENTauMomPDGIDvsPt_  = new TH2F("GENTauMomPDGIDvsPt" , "", 601, -.5, 600.5, 80, 0.0, 80);
  GENBMomPDGIDvsPt_  = new TH2F("GENBMomPDGIDvsPt" , "", 601, -.5, 600.5, 150, 0.0, 150);

  TauPt_     = new TH1F("TauPt"    , "", 201, -1, 200.0);
  TauPtDiff_ = new TH1F("TauPtDiff", "", 101, -1, 100.0);
  TauDR_     = new TH1F("TauDR"    , "", 100, 0, .4);
  DiTauDR_   = new TH1F("DiTauDR"  , "", 100, 0, 8);

  JetPt_     = new TH1F("JetPt",     "", 300, -1, 300.0);
  JetPtDiff_ = new TH1F("JetPtDiff", "", 100, -1.1, 100.0);
  JetDR_     = new TH1F("JetDR",     "", 100, 0, .4);
  BMatchFlavor_= new TH1F("BMatchFlavor", "", 7, -.5, 6.5);

  TauIsoDisc_ = new TH1F("TauIsoDisc", "", 4, -.5, 3.5);
      TauIsoDisc_->GetXaxis()->SetBinLabel(1, "Total Taus");
      TauIsoDisc_->GetXaxis()->SetBinLabel(2, "Loose");
      TauIsoDisc_->GetXaxis()->SetBinLabel(3, "Medium");
      TauIsoDisc_->GetXaxis()->SetBinLabel(4, "Tight");
  TauMediumIsoEff_     = new TH1F("TauMediumIsoEff"    , "", 2, -.5, 1.5);
      TauMediumIsoEff_->GetXaxis()->SetBinLabel(1, "Total");
      TauMediumIsoEff_->GetXaxis()->SetBinLabel(2, "Matched");

  TauIsovsPt_  = new TH2F("TauIsovsPt" , "", 4, 0, 4, 200, 0.0, 200);
      TauIsovsPt_->GetXaxis()->SetBinLabel(1, "Total Taus");
      TauIsovsPt_->GetXaxis()->SetBinLabel(2, "Loose");
      TauIsovsPt_->GetXaxis()->SetBinLabel(3, "Medium");
      TauIsovsPt_->GetXaxis()->SetBinLabel(4, "Tight");
  DiTaudRvsBPt_  = new TH2F("DiTaudRvsBPt" , "", 30, 0, 8, 300, 0.0, 300);
  TauNearestBdRvsBPt_  = new TH2F("TauNearestBdRvsBPt" , "", 30, 0, 8, 300, 0.0, 300);
}

// ------------ method called once each job just after ending the event loop  ------------
void TTbarAnalyzer::endJob()
{


  jetOutput_.close();

  //Make the Canvases
  TCanvas NEventsCanvas("NEvents","",600,600);
  TCanvas GENDiTauDRCanvas("GENDiTauDR","",600,600);
  TCanvas nGenTauPartCanvas("nGenTauPart","",600,600);
  TCanvas nGenBPartCanvas("nGenBPart","",600,600);
  TCanvas GenTauMomPDGIDCanvas("GenTauMomPDGID","",600,600);
  TCanvas GenBMomPDGIDCanvas("GenBMomPDGID","",600,600);
  TCanvas GENTauMomPDGIDvsPtCanvas("GENTauMomPDGIDvsPt","",600,600);
  TCanvas GENBMomPDGIDvsPtCanvas("GENBMomPDGIDvsPt","",600,600);

  TCanvas TauPtCanvas("TauPt","",600,600);
  TCanvas TauPtDiffCanvas("TauPtDiff","",600,600);
  TCanvas TauDRCanvas("TauDR","",600,600);
  TCanvas DiTauDRCanvas("DiTauDR","",600,600);

  TCanvas JetPtCanvas("JetPt","",600,600);
  TCanvas JetPtDiffCanvas("JetPtDiff","",600,600);
  TCanvas JetDRCanvas("JetDR","",600,600);
  TCanvas BMatchFlavorCanvas("BMatchFlavor","",600,600);

  TCanvas TauIsoDiscCanvas("TauIsoDisc","",600,600);
  TCanvas TauMediumIsoEffCanvas("TauMediumIsoEff","",600,600);

  TCanvas TauIsovsPtCanvas("TauIsovsPt","",600,600);
  TCanvas DiTaudRvsBPtCanvas("DiTaudRvsBPt","",600,600);
  TCanvas TauNearestBdRvsBPtCanvas("TauNearestBdRvsBPt","",600,600);

std::cout << "<----------------Declared Canvases-------------->" << std::endl;

  //Format the 1D plots and Draw (canvas, hist, grid, log y, log z, color, size, style, xAxisTitle, xTitleSize, xLabelSize, xTitleOffSet, yAxisTitle, yTitleSize, yLabelSize, yTitleOffset)
  VariousFunctions::formatAndDrawCanvasAndHist1D(NEventsCanvas, NEvents_, 1, 0, 0, kBlack, 7, 20, "Number of Taus", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(GENDiTauDRCanvas, GENDiTauDR_, 1, 0, 0, kBlack, 7, 20, "GEN #Delta R(#tau#tau)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(nGenTauPartCanvas, nGenTauPart_, 1, 0, 0, kBlack, 7, 20, "number of Tau Part", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(nGenBPartCanvas, nGenBPart_, 1, 0, 0, kBlack, 7, 20, "number of b Part", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(GenTauMomPDGIDCanvas, GenTauMomPDGID_, 1, 0, 0, kBlack, 7, 20, "Tau Mom PDGID", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(GenBMomPDGIDCanvas, GenBMomPDGID_, 1, 0, 0, kBlack, 7, 20, "b Mom PDGID", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist2D(GENTauMomPDGIDvsPtCanvas, GENTauMomPDGIDvsPt_, 1, 0, 0, kBlack, 7, 20, "GEN Tau PDG ID", .04, .04, 1.1, "Pt(#tau)", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(GENBMomPDGIDvsPtCanvas, GENBMomPDGIDvsPt_, 1, 0, 0, kBlack, 7, 20, "GEN b PDG ID", .04, .04, 1.1, "GEN Pt(b)", .04, .04, 1.6, "", .04, .04, 1.0);

  VariousFunctions::formatAndDrawCanvasAndHist1D(TauPtCanvas, TauPt_, 1, 0, 0, kBlack, 7, 20, "pt(Slimmmed Taus)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauPtDiffCanvas, TauPtDiff_, 1, 0, 0, kBlack, 7, 20, "#Deltapt(Slimmmed Taus and Gen Tau)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauDRCanvas, TauDR_, 1, 0, 0, kBlack, 7, 20, "#Delta R(Slimmmed Taus and Gen Tau)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiTauDRCanvas, DiTauDR_, 1, 0, 0, kBlack, 7, 20, "#Delta R(#tau #tau)", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(JetPtCanvas, JetPt_, 1, 0, 0, kBlack, 7, 20, "pt(Slimmmed Jets)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(JetPtDiffCanvas, JetPtDiff_, 1, 0, 0, kBlack, 7, 20, "#Delta pt(Slimmmed Jets and Gen B)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(JetDRCanvas, JetDR_, 1, 0, 0, kBlack, 7, 20, "#Delta R(Slimmmed Jets and GenB)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(BMatchFlavorCanvas, BMatchFlavor_, 1, 0, 0, kBlack, 7, 20, "Flavor of Matched bJets", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(TauIsoDiscCanvas, TauIsoDisc_, 1, 0, 0, kBlack, 7, 20, "Iso(slimmmed Taus)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauMediumIsoEffCanvas, TauMediumIsoEff_, 1, 0, 0, kBlack, 7, 20, "Taus Passing Medium Isolation", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist2D(TauIsovsPtCanvas, TauIsovsPt_, 1, 0, 0, kBlack, 7, 20, "Tau Iso", .04, .04, 1.1, "Pt(#tau)", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(DiTaudRvsBPtCanvas, DiTaudRvsBPt_, 1, 0, 0, kBlack, 7, 20, "#Delta R(#tau #tau)", .04, .04, 1.1, "Pt(b)", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(TauNearestBdRvsBPtCanvas, TauNearestBdRvsBPt_, 1, 0, 0, kBlack, 7, 20, "#Delta R(#tau nearest b)", .04, .04, 1.1, "Pt(nearest b)", .04, .04, 1.6, "", .04, .04, 1.0);

std::cout << "<----------------Formatted Canvases and Histos-------------->" << std::endl;

  //Write output file
  out_->cd();

  NEventsCanvas.Write();
  GENDiTauDRCanvas.Write();
  nGenTauPartCanvas.Write();
  nGenBPartCanvas.Write();
  GenTauMomPDGIDCanvas.Write();
  GenBMomPDGIDCanvas.Write();
  GENTauMomPDGIDvsPtCanvas.Write();
  GENBMomPDGIDvsPtCanvas.Write();

  TauPtCanvas.Write();
  TauPtDiffCanvas.Write();
  TauDRCanvas.Write();
  DiTauDRCanvas.Write();

  JetPtCanvas.Write();
  JetPtDiffCanvas.Write();
  JetDRCanvas.Write();
  BMatchFlavorCanvas.Write();

  TauIsoDiscCanvas.Write();
  TauMediumIsoEffCanvas.Write();

  TauIsovsPtCanvas.Write();
  DiTaudRvsBPtCanvas.Write();
  TauNearestBdRvsBPtCanvas.Write();

  out_->Write();
  out_->Close();
std::cout << "DONE" << std::endl;
}//EndJob

// ------------ method called when starting to processes a run  ------------
void TTbarAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void TTbarAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void TTbarAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void TTbarAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

//Delete Memory
void TTbarAnalyzer::reset(const bool doDelete)
{
  if ((doDelete) && (NEvents_ != NULL)) delete NEvents_;
  NEvents_ = NULL;
  if ((doDelete) && (GENDiTauDR_ != NULL)) delete GENDiTauDR_;
  GENDiTauDR_ = NULL;
  if ((doDelete) && (nGenTauPart_ != NULL)) delete nGenTauPart_;
  nGenTauPart_ = NULL;
  if ((doDelete) && (nGenBPart_ != NULL)) delete nGenBPart_;
  nGenBPart_ = NULL;
  if ((doDelete) && (GenTauMomPDGID_ != NULL)) delete GenTauMomPDGID_;
  GenTauMomPDGID_ = NULL;
  if ((doDelete) && (GenBMomPDGID_ != NULL)) delete GenBMomPDGID_;
  GenBMomPDGID_ = NULL;
  if ((doDelete) && (GENTauMomPDGIDvsPt_ != NULL)) delete GENTauMomPDGIDvsPt_;
  GENTauMomPDGIDvsPt_ = NULL;
  if ((doDelete) && (GENBMomPDGIDvsPt_ != NULL)) delete GENBMomPDGIDvsPt_;
  GENBMomPDGIDvsPt_ = NULL;

  if ((doDelete) && (TauPt_ != NULL)) delete TauPt_;
  TauPt_ = NULL;
  if ((doDelete) && (TauPtDiff_ != NULL)) delete TauPtDiff_;
  TauPtDiff_ = NULL;
  if ((doDelete) && (TauDR_ != NULL)) delete TauDR_;
  TauDR_ = NULL;
  if ((doDelete) && (DiTauDR_ != NULL)) delete DiTauDR_;
  DiTauDR_ = NULL;

  if ((doDelete) && (JetPt_ != NULL)) delete JetPt_;
  JetPt_ = NULL;
  if ((doDelete) && (JetPtDiff_ != NULL)) delete JetPtDiff_;
  JetPtDiff_ = NULL;
  if ((doDelete) && (JetDR_ != NULL)) delete JetDR_;
  JetDR_ = NULL;
  if ((doDelete) && (BMatchFlavor_ != NULL)) delete BMatchFlavor_;
  BMatchFlavor_ = NULL;

  if ((doDelete) && (TauIsoDisc_ != NULL)) delete TauIsoDisc_;
  TauIsoDisc_ = NULL;
  if ((doDelete) && (TauMediumIsoEff_ != NULL)) delete TauMediumIsoEff_;
  TauMediumIsoEff_ = NULL;

  if ((doDelete) && (TauIsovsPt_ != NULL)) delete TauIsovsPt_;
  TauIsovsPt_ = NULL;
  if ((doDelete) && (DiTaudRvsBPt_ != NULL)) delete DiTaudRvsBPt_;
  DiTaudRvsBPt_ = NULL;
  if ((doDelete) && (TauNearestBdRvsBPt_ != NULL)) delete TauNearestBdRvsBPt_;
  TauNearestBdRvsBPt_ = NULL;

}//void TTbarAnalyzer

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TTbarAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TTbarAnalyzer);

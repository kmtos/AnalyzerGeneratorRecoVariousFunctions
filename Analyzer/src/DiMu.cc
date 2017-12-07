// -*- C++ -*-
//
// Package:    DiMu
// Class:      DiMu
// 
/**\class DiMu DiMu.cc Analyzer/src/DiMu.cc

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
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

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

class DiMu : public edm::EDAnalyzer {
   public:
      typedef reco::JetFloatAssociation::Container JetBCEnergyRatioCollection;
      explicit DiMu(const edm::ParameterSet&);
      ~DiMu();

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
      edm::EDGetTokenT<edm::View<pat::Muon> > mu12Tag_;

      //Histograms
      TH1F* NEvents_;   
      TH1F* InvMassUpsilonRange_;
      TH1F* InvMassZPeakRange_;
      TH1F* InvMassFullRange_;
      TH1F* Mu1Pt_;
      TH1F* Mu2Pt_;
      TH1F* DiMuPt_;
      TH1F* Mu1Eta_;
      TH1F* Mu2Eta_;
      TH1F* DiMuEta_;
      TH1F* DiMudR_;
      TH1F* DiMuPhi_;
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
DiMu::DiMu(const edm::ParameterSet& iConfig):
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  mu12Tag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("mu12Tag")))
{
  reset(false);    
}//DiMu



DiMu::~DiMu()
{
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void DiMu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "\n<------------THIS IS A NEW EVENT------------>" << std::endl;

  //Old Jet collection for bTagging
  edm::Handle<edm::View<pat::Muon> > pMu12;
  iEvent.getByToken(mu12Tag_, pMu12);

  pat::Muon mu1 = pat::Muon( (*pMu12)[0] );
  pat::Muon mu2 = pat::Muon( (*pMu12)[1] );
  reco::LeafCandidate::LorentzVector diMuP4 = mu1.p4() + mu2.p4();
  InvMassZPeakRange_->Fill(diMuP4.M() );
  InvMassUpsilonRange_->Fill(diMuP4.M() );
  InvMassFullRange_->Fill(diMuP4.M() );
  Mu1Pt_->Fill(mu1.pt() );
  Mu2Pt_->Fill(mu2.pt() );
  DiMuPt_->Fill(diMuP4.Pt() );
  Mu1Eta_->Fill(mu1.eta() );
  Mu2Eta_->Fill(mu2.eta() );
  DiMuEta_->Fill(diMuP4.Eta() );
  DiMuPhi_->Fill(diMuP4.Phi() );
  DiMudR_->Fill(reco::deltaR(mu1, mu2) );
   
}//End InvMass::analyze


// ------------ method called once each job just before starting event loop  ------------
void DiMu::beginJob()
{
  std::cout << "Begin Job" << std::endl;
  VariousFunctions::setTDRStyle(true); //NEW

  //Open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");


  //Book histograms
  NEvents_     = new TH1F("NEvents"    , "", 9, -.5, 8.5);
      NEvents_->GetXaxis()->SetBinLabel(1, "TotalEvents"); 
      NEvents_->GetXaxis()->SetBinLabel(2, "#tau_{#mu} + #tau_{had} Match");
      NEvents_->GetXaxis()->SetBinLabel(3, "Gen #tau_{#mu} + #tau_{had}");
      NEvents_->GetXaxis()->SetBinLabel(4, "Gen Match #tau_{had}");
      NEvents_->GetXaxis()->SetBinLabel(5, "Gen Match #tau_{had}");
      NEvents_->GetXaxis()->SetBinLabel(6, "Event with #tau_{#mu} Removed");
      NEvents_->GetXaxis()->SetBinLabel(7, "Event with no #tau_{#mu} Removed ");
  InvMassUpsilonRange_     = new TH1F("InvMassUpsilonRange"    , "", 50, 5, 15);
  InvMassZPeakRange_     = new TH1F("InvMassZPeakRange"    , "", 200, 50, 150);
  InvMassFullRange_     = new TH1F("InvMassFullRange"    , "", 450, 0, 150);
  Mu1Pt_     = new TH1F("Mu1Pt"    , "", 75, 0, 500);
  Mu2Pt_     = new TH1F("Mu2Pt"    , "", 75, 0, 500);
  DiMuPt_     = new TH1F("DiMuPt"    , "", 75, 0, 500);
  Mu1Eta_     = new TH1F("Mu1Eta"    , "", 75, -2.5, 2.5);
  Mu2Eta_     = new TH1F("Mu2Eta"    , "", 75, -2.5, 2.5);
  DiMuEta_     = new TH1F("DiMuEta"    , "", 75, -2.5, 2.5);
  DiMudR_     = new TH1F("DiMudR"    , "", 75, 0, 3.5);
  DiMuPhi_     = new TH1F("DiMuPhi"    , "", 75, 0, 6.5);

}

// ------------ method called once each job just after ending the event loop  ------------
void DiMu::endJob()
{
  //Make the Canvases
  TCanvas NEventsCanvas("NEvents","",600,600);
  TCanvas InvMassUpsilonRangeCanvas("InvMassUpsilonRange","",600,600);
  TCanvas InvMassZPeakRangeCanvas("InvMassZPeakRange","",600,600);
  TCanvas *InvMassFullRangeCanvas = new TCanvas("InvMassFullRange","",600,600);
  TCanvas Mu1PtCanvas("Mu1Pt","",600,600);
  TCanvas Mu2PtCanvas("Mu2Pt","",600,600);
  TCanvas DiMuPtCanvas("DiMuPt","",600,600);
  TCanvas Mu1EtaCanvas("Mu1Eta","",600,600);
  TCanvas Mu2EtaCanvas("Mu2Eta","",600,600);
  TCanvas DiMuEtaCanvas("DiMuEta","",600,600);
  TCanvas DiMudRCanvas("DiMudR","",600,600);
  TCanvas DiMuPhiCanvas("DiMuPhi","",600,600);

std::cout << "<----------------Declared Canvases-------------->" << std::endl;

  //Format the 1D plots and Draw (canvas, hist, grid, log y, log z, color, size, style, xAxisTitle, xTitleSize, xLabelSize, xTitleOffSet, yAxisTitle, yTitleSize, yLabelSize, yTitleOffset)
  VariousFunctions::formatAndDrawCanvasAndHist1D(NEventsCanvas, NEvents_,
	 1, 0, 0, kBlack, 7, 20, "", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassUpsilonRangeCanvas, InvMassUpsilonRange_,
	 1, 0, 0, kBlack, 7, 20, "Mass(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassZPeakRangeCanvas, InvMassZPeakRange_,
	 1, 0, 0, kBlack, 7, 20, "Mass(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
//  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassFullRangeCanvas, InvMassFullRange_,
//	 1, 0, 0, kBlack, 7, 20, "Mass(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(Mu1PtCanvas, Mu1Pt_,
	 1, 0, 0, kBlack, 7, 20, "p_{T}(#mu_{1})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(Mu2PtCanvas, Mu2Pt_,
	 1, 0, 0, kBlack, 7, 20, "p_{T}(#mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiMuPtCanvas, DiMuPt_,
	 1, 0, 0, kBlack, 7, 20, "p_{T}(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(Mu1EtaCanvas, Mu1Eta_,
	 1, 0, 0, kBlack, 7, 20, "#eta(#mu_{1})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(Mu2EtaCanvas, Mu2Eta_,
	 1, 0, 0, kBlack, 7, 20, "#eta(#mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiMuEtaCanvas, DiMuEta_,
	 1, 0, 0, kBlack, 7, 20, "#eta(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiMudRCanvas, DiMudR_,
	 1, 1, 0, kBlack, 7, 20, "#DeltaR(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiMuPhiCanvas, DiMuPhi_,
	 1, 1, 0, kBlack, 7, 20, "#Phi(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::SetCanvas(InvMassFullRangeCanvas, InvMassFullRangeCanvas->GetWw(), InvMassFullRangeCanvas->GetWh());
  VariousFunctions::SetHist(InvMassFullRange_, kBlack, kBlack, kWhite, 1.0, "InvMass(#mu_{1}#mu}_{2})", "# of Fakes");
  VariousFunctions::CMS_lumi(InvMassFullRangeCanvas, 4, 33, true, true, true);

std::cout << "after formatting" << std::endl;
  
////////////////////////
// For Pt Gen
////////////////////////
std::cout << "<----------------Formatted Canvases and Histos-------------->" << std::endl;

  //Write output file
  out_->cd();

  NEventsCanvas.Write();
  InvMassUpsilonRangeCanvas.Write();
  InvMassZPeakRangeCanvas.Write();
  InvMassFullRangeCanvas->Write();
  Mu1PtCanvas.Write();
  Mu2PtCanvas.Write();
  DiMuPtCanvas.Write();
  Mu1EtaCanvas.Write();
  Mu2EtaCanvas.Write();
  DiMuEtaCanvas.Write();
  DiMudRCanvas.Write();
  DiMuPhiCanvas.Write();

  out_->Write();
  out_->Close();
std::cout << "DONE" << std::endl;
}//EndJob

// ------------ method called when starting to processes a run  ------------
void DiMu::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void DiMu::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DiMu::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DiMu::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

//Delete Memory
void DiMu::reset(const bool doDelete)
{
  if ((doDelete) && (NEvents_ != NULL)) delete NEvents_;
  NEvents_ = NULL;
  if ((doDelete) && (InvMassUpsilonRange_ != NULL)) delete InvMassUpsilonRange_;
  InvMassUpsilonRange_ = NULL;
  if ((doDelete) && (InvMassZPeakRange_ != NULL)) delete InvMassZPeakRange_;
  InvMassZPeakRange_ = NULL;
  if ((doDelete) && (InvMassFullRange_ != NULL)) delete InvMassFullRange_;
  InvMassFullRange_ = NULL;
  if ((doDelete) && (Mu1Pt_ != NULL)) delete Mu1Pt_;
  Mu1Pt_ = NULL;
  if ((doDelete) && (Mu2Pt_ != NULL)) delete Mu2Pt_;
  Mu2Pt_ = NULL;
  if ((doDelete) && (DiMuPt_ != NULL)) delete DiMuPt_;
  DiMuPt_ = NULL;
  if ((doDelete) && (Mu1Eta_ != NULL)) delete Mu1Eta_;
  Mu1Eta_ = NULL;
  if ((doDelete) && (Mu2Eta_ != NULL)) delete Mu2Eta_;
  Mu2Eta_ = NULL;
  if ((doDelete) && (DiMuEta_ != NULL)) delete DiMuEta_;
  DiMuEta_ = NULL;
  if ((doDelete) && (DiMudR_ != NULL)) delete DiMudR_;
  DiMudR_ = NULL;
  if ((doDelete) && (DiMuPhi_ != NULL)) delete DiMuPhi_;
  DiMuPhi_ = NULL;

}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiMu::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMu);

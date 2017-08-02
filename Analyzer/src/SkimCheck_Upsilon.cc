// -*- C++ -*-
//
// Package:    SkimCheck_Upsilon
// Class:      SkimCheck_Upsilon
// 
/**\class SkimCheck_Upsilon SkimCheck_Upsilon.cc Analyzer/src/SkimCheck_Upsilon.cc

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

class SkimCheck_Upsilon : public edm::EDAnalyzer {
   public:
      typedef reco::JetFloatAssociation::Container JetBCEnergyRatioCollection;
      explicit SkimCheck_Upsilon(const edm::ParameterSet&);
      ~SkimCheck_Upsilon();

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
      edm::EDGetTokenT<vector<reco::PFJet> > akJetTag_;
      edm::EDGetTokenT<vector<reco::Muon> > muonsTag_;
      edm::EDGetTokenT<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > mu12Tag_;

      //Histograms
      TH1F* MassDiMuRECO_;
      TH1F* NEventsCuts_;
      TH1F* PtOverMassMu1Mu2_;
      TH1F* PtMu1Mu2_;
      TH1F* PtMu1_;
      TH1F* PTMu2_;
      TH1F* dRMu1Mu2_;
      TH1F* RelIsoRaw_;
      TH1F* IsoRawCutBased_;
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
SkimCheck_Upsilon::SkimCheck_Upsilon(const edm::ParameterSet& iConfig):
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  akJetTag_(consumes<vector<reco::PFJet> >(iConfig.getParameter<edm::InputTag>("akJetTag"))),
  muonsTag_(consumes<vector<reco::Muon> >(iConfig.getParameter<edm::InputTag>("muonsTag"))),
  mu12Tag_(consumes<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > >(iConfig.getParameter<edm::InputTag>("mu12Tag")))
{
  reset(false);    
}//SkimCheck_Upsilon



SkimCheck_Upsilon::~SkimCheck_Upsilon()
{
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void SkimCheck_Upsilon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "\n<------------THIS IS A NEW EVENT------------>" << std::endl;

  //Get ak4Jets particle collection
  edm::Handle<std::vector<reco::PFJet> > pAkJets;
  iEvent.getByToken(akJetTag_, pAkJets);

  //Get RECO Muons particle collection
  edm::Handle<std::vector<reco::Muon> > pMuons;
  iEvent.getByToken(muonsTag_, pMuons);

  //Old Jet collection for bTagging
  edm::Handle<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > pMu12;
  iEvent.getByToken(mu12Tag_, pMu12);

  std::cout << "Declared objects" << std::endl;
//////////////////////////////
// Begin Analyzer
//////////////////////////////
  //Getting Mu1 and Mu2
  reco::MuonRef mu1Ref = reco::MuonRef((*pMu12)[0] );
  reco::MuonRef mu2Ref = reco::MuonRef((*pMu12)[1] );
  PtMu1_->Fill(mu1Ref->pt() );
  PTMu2_->Fill(mu2Ref->pt() );
  double dPhi = reco::deltaPhi(mu1Ref->phi(), mu2Ref->phi() ), dEta = mu1Ref->eta() - mu2Ref->eta();
  dRMu1Mu2_->Fill( sqrt( dPhi*dPhi + dEta*dEta) );
  reco::LeafCandidate::LorentzVector diMuP4;
  diMuP4 = mu1Ref->p4();
  diMuP4 += mu2Ref->p4();
  PtOverMassMu1Mu2_->Fill(diMuP4.Pt() / diMuP4.M() );
  MassDiMuRECO_->Fill(diMuP4.M() );
  PtMu1Mu2_->Fill(diMuP4.Pt() );  //mu1Ref->pt() + mu2Ref->pt() );
}//End SkimCheck_Upsilon::analyze


// ------------ method called once each job just before starting event loop  ------------
void SkimCheck_Upsilon::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  //Open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //Book histograms
  MassDiMuRECO_       = new TH1F("MassDiMuRECO"    , "", 150, 5, 15);
  NEventsCuts_       = new TH1F("NEventsCuts", "", 13, -.5, 12.5);
      NEventsCuts_->GetXaxis()->SetBinLabel(1, "TotalEvents");
      NEventsCuts_->GetXaxis()->SetBinLabel(2, "BD < 5"); 
      NEventsCuts_->GetXaxis()->SetBinLabel(3, "BD < .75"); 
      NEventsCuts_->GetXaxis()->SetBinLabel(4, "BD < .9");
      NEventsCuts_->GetXaxis()->SetBinLabel(5, "BD < .9 + Med Iso");
      NEventsCuts_->GetXaxis()->SetBinLabel(6, "HT < 300");
      NEventsCuts_->GetXaxis()->SetBinLabel(7, "HT < 400");
      NEventsCuts_->GetXaxis()->SetBinLabel(8, "HT < 300 , BD < .75"); 
      NEventsCuts_->GetXaxis()->SetBinLabel(9, "HT < 400 , BD < .5");
      NEventsCuts_->GetXaxis()->SetBinLabel(10, "IsoRaw < .5");
      NEventsCuts_->GetXaxis()->SetBinLabel(11, "IsoRaw < 1"); 
      NEventsCuts_->GetXaxis()->SetBinLabel(12, "Tight Iso"); 
      NEventsCuts_->GetXaxis()->SetBinLabel(13, "VTight Iso");
  PtOverMassMu1Mu2_        = new TH1F("PtOverMassMu1Mu2"    , "", 50, 0, 50);
  PtMu1Mu2_        = new TH1F("PtMu1Mu2"    , "", 50, 0, 700);
  PtMu1_        = new TH1F("PtMu1"    , "", 100, 0, 500);
  PTMu2_        = new TH1F("PTMu2"    , "", 100, 0, 500);
  dRMu1Mu2_        = new TH1F("dRMu1Mu2"    , "", 100, 0, 5);
  RelIsoRaw_        = new TH1F("RelIsoRaw"    , "", 100, 0, 5);
  IsoRawCutBased_        = new TH1F("IsoRawCutBased"    , "", 100, 0, 10);

}

// ------------ method called once each job just after ending the event loop  ------------
void SkimCheck_Upsilon::endJob()
{
  //Make the Canvases
  TCanvas MassDiMuRECOCanvas("MassDiMuRECO","",600,600);
  TCanvas NEventsCutsCanvas("NEventsCuts","",600,600);
  TCanvas PtOverMassMu1Mu2Canvas("PtOverMassMu1Mu2","",600,600);
  TCanvas PtMu1Mu2Canvas("PtMu1Mu2","",600,600);
  TCanvas PtMu1Canvas("PtMu1","",600,600);
  TCanvas PTMu2Canvas("PTMu2","",600,600);
  TCanvas dRMu1Mu2Canvas("dRMu1Mu2","",600,600);
  TCanvas RelIsoRawCanvas("RelIsoRaw","",600,600);
  TCanvas IsoRawCutBasedCanvas("IsoRawCutBased","",600,600);

std::cout << "<----------------Declared Canvases-------------->" << std::endl;

  //Format the 1D plots and Draw (canvas, hist, grid, log y, log z, color, size, style, xAxisTitle, xTitleSize, xLabelSize, xTitleOffSet, yAxisTitle, yTitleSize, yLabelSize, yTitleOffset)
  VariousFunctions::formatAndDrawCanvasAndHist1D(MassDiMuRECOCanvas, MassDiMuRECO_,
         1, 0, 0, kBlack, 7, 20, "RECO m_{#mu #mu}", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NEventsCutsCanvas, NEventsCuts_,
         1, 0, 0, kBlack, 7, 20, "", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(PtOverMassMu1Mu2Canvas, PtOverMassMu1Mu2_,
         1, 0, 0, kBlack, 7, 20, "p_{T}(#mu_{1}#mu_{2}) / Mass(#mu_{1}#mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(PtMu1Mu2Canvas, PtMu1Mu2_,
         1, 0, 0, kBlack, 7, 20, "p_{T}(#mu_{1}#mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(PtMu1Canvas, PtMu1_,
         1, 0, 0, kBlack, 7, 20, "p_{T}(#mu_{1}) ", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(PTMu2Canvas, PTMu2_,
         1, 0, 0, kBlack, 7, 20, "p_{T}(#mu_{2}) ", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(dRMu1Mu2Canvas, dRMu1Mu2_,
         1, 0, 0, kBlack, 7, 20, "#Delta R(#mu_{1}#mu_{1})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(RelIsoRawCanvas, RelIsoRaw_,
         1, 0, 0, kBlack, 7, 20, "Raw Isolation Value / p_{T}(#tau_{Had})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoRawCutBasedCanvas, IsoRawCutBased_,
         1, 0, 0, kBlack, 7, 20, "Raw Isolation Value passing Med CutBased", .04, .04, 1.1,  "", .04, .04, 1.0, false);

std::cout << "after formatting" << std::endl;
  
std::cout << "<----------------Formatted Canvases and Histos-------------->" << std::endl;

  //Write output file
  out_->cd();

  MassDiMuRECOCanvas.Write();
  NEventsCutsCanvas.Write();
  PtOverMassMu1Mu2Canvas.Write();
  PtMu1Mu2Canvas.Write();
  PtMu1Canvas.Write();
  PTMu2Canvas.Write();
  dRMu1Mu2Canvas.Write();
  RelIsoRawCanvas.Write();
  IsoRawCutBasedCanvas.Write();

  out_->Write();
  out_->Close();
std::cout << "DONE" << std::endl;
}//EndJob

// ------------ method called when starting to processes a run  ------------
void SkimCheck_Upsilon::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void SkimCheck_Upsilon::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void SkimCheck_Upsilon::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void SkimCheck_Upsilon::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

//Delete Memory
void SkimCheck_Upsilon::reset(const bool doDelete)
{
  if ((doDelete) && (MassDiMuRECO_ != NULL)) delete MassDiMuRECO_;
  MassDiMuRECO_ = NULL;
  if ((doDelete) && (NEventsCuts_ != NULL)) delete NEventsCuts_;
  NEventsCuts_ = NULL;
  if ((doDelete) && (PtOverMassMu1Mu2_ != NULL)) delete PtOverMassMu1Mu2_;
  PtOverMassMu1Mu2_ = NULL;
  if ((doDelete) && (PtMu1Mu2_ != NULL)) delete PtMu1Mu2_;
  PtMu1Mu2_ = NULL;
  if ((doDelete) && (PtMu1_ != NULL)) delete PtMu1_;
  PtMu1_ = NULL;
  if ((doDelete) && (PTMu2_ != NULL)) delete PTMu2_;
  PTMu2_ = NULL;
  if ((doDelete) && (dRMu1Mu2_ != NULL)) delete dRMu1Mu2_;
  dRMu1Mu2_ = NULL;
  if ((doDelete) && (RelIsoRaw_ != NULL)) delete RelIsoRaw_;
  RelIsoRaw_ = NULL;
  if ((doDelete) && (IsoRawCutBased_ != NULL)) delete IsoRawCutBased_;
  IsoRawCutBased_ = NULL;

}//void SkimCheck_Upsilon

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SkimCheck_Upsilon::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SkimCheck_Upsilon);

// -*- C++ -*-
//
// Package:    FakeRateMiniAODMCGetRegionA
// Class:      FakeRateMiniAODMCGetRegionA
// 
/**\class FakeRateMiniAODMCGetRegionA FakeRateMiniAODMCGetRegionA.cc Analyzer/src/FakeRateMiniAODMCGetRegionA.cc

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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

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

class FakeRateMiniAODMCGetRegionA : public edm::EDAnalyzer {
   public:
      typedef reco::JetFloatAssociation::Container JetBCEnergyRatioCollection;
      explicit FakeRateMiniAODMCGetRegionA(const edm::ParameterSet&);
      ~FakeRateMiniAODMCGetRegionA();

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
      TFile* PileupFile;

      //name of output root file
      std::string outFileName_;
      edm::EDGetTokenT<edm::View<pat::Tau> > tauTag_;
      double mu3dRMin_;
      double mu3dRMax_;
      double tauPtCut_;
      edm::EDGetTokenT<edm::View<pat::Muon> > mu3Tag_;
      edm::EDGetTokenT<edm::View<pat::Muon> > mu12Tag_;
      bool requireRemovedMuon_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryInfo_;
      edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
      double xsec_;
      double lumi_;
      double summedWeights_;
      std::string PileupFileName_;


      //Histograms
      TH1F* NEvents_;   
      TH1F* InvMassTauMuMu1_;
      TH1F* InvMassTauMuMu2_;
      TH1F* DiMuInvMassFakeWeight_;
      TH1F* DiMuDiTauInvMassFakeWeight_;
      TH1F* DiTauInvMassFakeWeight_;
      TH1F* PtMu1FakeWeight_;
      TH1F* PtMu2FakeWeight_;
      TH1F* EtaFakeWeight_;
      TH1F* DRFakeWeight_;
      TH1F* DRNoWeighting_;
      TH1F* DiMuInvMassFakeWeightZoom_;
      TH1F* TauVisMass_;
      TH1F* TauVisMassZoom_;

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
FakeRateMiniAODMCGetRegionA::FakeRateMiniAODMCGetRegionA(const edm::ParameterSet& iConfig):

  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  tauTag_(consumes<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("tauTag"))),
  mu3dRMin_(iConfig.getParameter<double>("mu3dRMin")),
  mu3dRMax_(iConfig.getParameter<double>("mu3dRMax")),
  tauPtCut_(iConfig.getParameter<double>("tauPtCut")),
  mu3Tag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("mu3Tag"))),
  mu12Tag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("mu12Tag"))),
  requireRemovedMuon_(iConfig.getParameter<bool>("requireRemovedMuon")),
  pileupSummaryInfo_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummaryInfo"))),
  genEventInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfoToken"))),
  xsec_(iConfig.getParameter<double>("xsec")),
  lumi_(iConfig.getParameter<double>("lumi")),
  summedWeights_(iConfig.getParameter<double>("summedWeights")),
  PileupFileName_(iConfig.getParameter<std::string>("PileupFileName"))
{
  reset(false);    
}//FakeRateMiniAODMCGetRegionA



FakeRateMiniAODMCGetRegionA::~FakeRateMiniAODMCGetRegionA()
{
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void FakeRateMiniAODMCGetRegionA::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "\n<------------THIS IS A NEW EVENT------------>" << std::endl;
  NEvents_->Fill(0);

  //Get CleanJets Tau particle collection
  edm::Handle<edm::View<pat::Tau> > pTaus;
  iEvent.getByToken(tauTag_, pTaus);

  //Old Jet collection for bTagging
  edm::Handle<edm::View<pat::Muon> > pMu3;
  iEvent.getByToken(mu3Tag_, pMu3);

  //Old Jet collection for bTagging
  edm::Handle<edm::View<pat::Muon> > pMu12;
  iEvent.getByToken(mu12Tag_, pMu12);

  double highestMuPt = -1;
  pat::Muon mu1, mu2;
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

  std::cout << "mu1.pt()=" << mu1.pt() << "\tmu2.pt()=" << mu2.pt() << std::endl;

  reco::LeafCandidate::LorentzVector diMuP4 = mu1.p4() + mu2.p4();
  if (diMuP4.M() > 30.0)
    return;

//////////////////////////////
// Begin Analyzer
//////////////////////////////
  ///////////////////////////
  // iterating over the taus
  /////////////////////////// 
  reco::LeafCandidate::LorentzVector  diTauP4, diMuDiTauP4;
  double bestMu3dR = 10000000;
  bool checkEventDiTau = false;
  pat::Muon mu3;
std::cout << "pTaus->size()= " << pTaus->size() << "\tpMu3->size()=" << pMu3->size() << std::endl;
  for (edm::View<pat::Tau>::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); ++iTau)
  {
    if (iTau->pt() < tauPtCut_ || (fabs(iTau->eta() )) > 2.4)
      continue;
    std::cout <<"PASS" << std::endl;
    for (edm::View<pat::Muon>::const_iterator iMu = pMu3->begin(); iMu != pMu3->end(); ++iMu)
    {
      double currdR = deltaR(*iTau, *iMu);
      std::cout << "currdR=" << currdR << "\tmu3dRMax_=" << mu3dRMax_ << "\tmu3dRMin_=" << mu3dRMin_ << std::endl;
      if (currdR < mu3dRMax_  && currdR > mu3dRMin_ && currdR < bestMu3dR)
      {
        diTauP4 = iTau->p4() + iMu->p4();
        diMuDiTauP4 = mu1.p4() + mu2.p4() + iTau->p4() + iMu->p4();
        bestMu3dR = currdR;
        checkEventDiTau = true;
        mu3 = *iMu;
      }//if
    }//for iMu
  }//for iTau

  std::cout << "checkEventDiTau=" << checkEventDiTau << std::endl;
  //Checking that there was a tau with a removed muon
  if (!checkEventDiTau && requireRemovedMuon_)
    return;

  reco::LeafCandidate::LorentzVector Mu1Mu3P4, Mu2Mu3P4;    
  Mu1Mu3P4 = mu1.p4() + mu3.p4();
  Mu2Mu3P4 = mu2.p4() + mu3.p4();

  InvMassTauMuMu1_->Fill(Mu1Mu3P4.M() );
  InvMassTauMuMu2_->Fill(Mu2Mu3P4.M() );


  edm::Handle<std::vector<PileupSummaryInfo> > pPileupSummaryInfo;
  iEvent.getByToken(pileupSummaryInfo_, pPileupSummaryInfo);

  int nTrueVertices = 0;
  if (pPileupSummaryInfo.isValid() && pPileupSummaryInfo->size()>0)
    nTrueVertices = pPileupSummaryInfo->at(1).getTrueNumInteractions();

  PileupFile = new TFile(PileupFileName_.c_str());
  TH1F* Pileup_ = (TH1F*)PileupFile->Get("PileupWeights");
  double pileupWeight = Pileup_->GetBinContent(nTrueVertices);
  PileupFile->Close();

  edm::Handle<GenEventInfoProduct> genEventInfo;
  iEvent.getByToken(genEventInfoToken_, genEventInfo);
  double eventGenWeight = genEventInfo->weight();
  double genWeight = eventGenWeight * lumi_ * xsec_ / summedWeights_;
  std::cout << "genWeight=" << genWeight << "\teventGenWeight=" << eventGenWeight << "\tlumi_=" << lumi_  << "\tsummedWeights_=" << summedWeights_ << "\txsec_=" << xsec_ << "\tpileupWeight=" <<pileupWeight <<  "\tnTrueVertices=" << nTrueVertices << std::endl;


  DiMuInvMassFakeWeight_->Fill(diMuP4.M(),  pileupWeight*genWeight );
  DiMuDiTauInvMassFakeWeight_->Fill(diMuDiTauP4.M(),  pileupWeight*genWeight );
  DiTauInvMassFakeWeight_->Fill(diTauP4.M(),  pileupWeight*genWeight );
  DiMuInvMassFakeWeightZoom_->Fill(diMuP4.M(),  pileupWeight*genWeight );
  TauVisMass_->Fill(diTauP4.M(),  pileupWeight*genWeight );
  TauVisMassZoom_->Fill(diTauP4.M(),  pileupWeight*genWeight );
  PtMu1FakeWeight_->Fill(mu1.pt(), pileupWeight*genWeight );
  PtMu2FakeWeight_->Fill(mu2.pt(), pileupWeight*genWeight );
  EtaFakeWeight_->Fill(mu1.eta(), pileupWeight*genWeight );
  double dR_Mu1Mu2 = deltaR(mu1, mu2);
  DRFakeWeight_->Fill(dR_Mu1Mu2, pileupWeight*genWeight );
  DRNoWeighting_->Fill(dR_Mu1Mu2);
}//End FakeRateMiniAODMCGetRegionA::analyze


// ------------ method called once each job just before starting event loop  ------------
void FakeRateMiniAODMCGetRegionA::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  //Open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");
  
//  TCanvas* FinalFakeRateLooseIsoEtavsPtCanvas = (TCanvas*)TH2File.Get("FinalFakeRateLooseIsoEtavsPtCanvas");
//  TCanvas* FinalFakeRateMedIsoEtavsPtCanvas = (TCanvas*)TH2File.Get("FinalFakeRateMedIsoEtavsPtCanvas");
//  TCanvas* FinalFakeRateTightIsoEtavsPtCanvas = (TCanvas*)TH2File.Get("FinalFakeRateTightIsoEtavsPtCanvas");
//  TCanvas* FinalFakeRateEtavsPtCanvas = (TCanvas*)TH2File.Get("FinalFakeRateEtavsPtCanvas");

//  TH2F* FinalFakeRateLooseIsoEtavsPt_ = (TH2F*)FinalFakeRateLooseIsoEtavsPtCanvas->GetPrimitive("FakeRateLooseIsoEtavsPt");
//  TH2F* FinalFakeRateMedIsoEtavsPt_ = (TH2F*)FinalFakeRateMedIsoEtavsPtCanvas->GetPrimitive("FakeRateMedIsoEtavsPt");
//  TH2F* FinalFakeRateTightIsoEtavsPt_ = (TH2F*)FinalFakeRateTightIsoEtavsPtCanvas->GetPrimitive("FakeRateTightIsoEtavsPt");
//  TH2F* FinalFakeRateEtavsPt_ = (TH2F*)FinalFakeRateEtavsPtCanvas->GetPrimitive("FakeRateEtavsPt");

//  TCanvas* FinalFakeRateLooseIsoEtavsPtSoftMuonCanvas = (TCanvas*)TH2File.Get("FinalFakeRateLooseIsoEtavsPtSoftMuonCanvas");
//  TCanvas* FinalFakeRateMedIsoEtavsPtSoftMuonCanvas = (TCanvas*)TH2File.Get("FinalFakeRateMedIsoEtavsPtSoftMuonCanvas");
//  TCanvas* FinalFakeRateTightIsoEtavsPtSoftMuonCanvas = (TCanvas*)TH2File.Get("FinalFakeRateTightIsoEtavsPtSoftMuonCanvas");
//  TCanvas* FinalFakeRateEtavsPtSoftMuonCanvas = (TCanvas*)TH2File.Get("FinalFakeRateEtavsPtSoftMuonCanvas");

//  TH2F* FinalFakeRateLooseIsoEtavsPtSoftMuon_ = (TH2F*)FinalFakeRateLooseIsoEtavsPtSoftMuonCanvas->GetPrimitive("FakeRateLooseIsoEtavsPtSoftMuon");
//  TH2F* FinalFakeRateMedIsoEtavsPtSoftMuon_ = (TH2F*)FinalFakeRateMedIsoEtavsPtSoftMuonCanvas->GetPrimitive("FakeRateMedIsoEtavsPtSoftMuon");
//  TH2F* FinalFakeRateTightIsoEtavsPtSoftMuon_ = (TH2F*)FinalFakeRateTightIsoEtavsPtSoftMuonCanvas->GetPrimitive("FakeRateTightIsoEtavsPtSoftMuon");
//  TH2F* FinalFakeRateEtavsPtSoftMuon_ = (TH2F*)FinalFakeRateEtavsPtSoftMuonCanvas->GetPrimitive("FakeRateEtavsPtSoftMuon");

//  TCanvas* FinalFakeRateLooseIsoEtavsPtSoftMuon_noMuCanvas = (TCanvas*)TH2File.Get("FinalFakeRateLooseIsoEtavsPtSoftMuon_noMuCanvas");
//  TCanvas* FinalFakeRateMedIsoEtavsPtSoftMuon_noMuCanvas = (TCanvas*)TH2File.Get("FinalFakeRateMedIsoEtavsPtSoftMuon_noMuCanvas");
//  TCanvas* FinalFakeRateTightIsoEtavsPtSoftMuon_noMuCanvas = (TCanvas*)TH2File.Get("FinalFakeRateTightIsoEtavsPtSoftMuon_noMuCanvas");
//  TCanvas* FinalFakeRateEtavsPtSoftMuon_noMuCanvas = (TCanvas*)TH2File.Get("FinalFakeRateEtavsPtSoftMuon_noMuCanvas");

//  TH2F* FinalFakeRateLooseIsoEtavsPtSoftMuon_noMu_ = (TH2F*)FinalFakeRateLooseIsoEtavsPtSoftMuon_noMuCanvas->GetPrimitive("FakeRateLooseIsoEtavsPtSoftMuon_noMu");
//  TH2F* FinalFakeRateMedIsoEtavsPtSoftMuon_noMu_ = (TH2F*)FinalFakeRateMedIsoEtavsPtSoftMuon_noMuCanvas->GetPrimitive("FakeRateMedIsoEtavsPtSoftMuon_noMu");
//  TH2F* FinalFakeRateTightIsoEtavsPtSoftMuon_noMu_ = (TH2F*)FinalFakeRateTightIsoEtavsPtSoftMuon_noMuCanvas->GetPrimitive("FakeRateTightIsoEtavsPtSoftMuon_noMu");
//  TH2F* FinalFakeRateEtavsPtSoftMuon_noMu_ = (TH2F*)FinalFakeRateEtavsPtSoftMuon_noMuCanvas->GetPrimitive("FakeRateEtavsPtSoftMuon_noMu");


  //Book histograms
  NEvents_     = new TH1F("NEvents"    , "", 9, -.5, 8.5);
      NEvents_->GetXaxis()->SetBinLabel(1, "TotalEvents"); 
      NEvents_->GetXaxis()->SetBinLabel(2, "#tau_{#mu} + #tau_{had} Match");
      NEvents_->GetXaxis()->SetBinLabel(3, "Gen #tau_{#mu} + #tau_{had}");
      NEvents_->GetXaxis()->SetBinLabel(4, "Gen Match #tau_{had}");
      NEvents_->GetXaxis()->SetBinLabel(5, "Gen Match #tau_{had}");
      NEvents_->GetXaxis()->SetBinLabel(6, "Event with #tau_{#mu} Removed");
      NEvents_->GetXaxis()->SetBinLabel(7, "Event with no #tau_{#mu} Removed ");
  InvMassTauMuMu1_     = new TH1F("InvMassTauMuMu1"    , "", 10, 0, 150);
  InvMassTauMuMu2_     = new TH1F("InvMassTauMuMu2"    , "", 10, 0, 150);
  DiMuInvMassFakeWeight_     = new TH1F("DiMuInvMassFakeWeight"    , "", 10, 0, 150);
  DiMuDiTauInvMassFakeWeight_     = new TH1F("DiMuDiTauInvMassFakeWeight"    , "", 333, 0, 1000);
  DiTauInvMassFakeWeight_     = new TH1F("DiTauInvMassFakeWeight"    , "", 150, 0, 150);
  PtMu1FakeWeight_     = new TH1F("PtMu1FakeWeight"    , "", 10, 0, 300);
  PtMu2FakeWeight_     = new TH1F("PtMu2FakeWeight"    , "", 10, 0, 300);
  EtaFakeWeight_     = new TH1F("EtaFakeWeight"    , "", 10, -2.5, 2.5);
  DRFakeWeight_     = new TH1F("DRFakeWeight"    , "", 10, 0, 5);
  DRNoWeighting_     = new TH1F("DRNoWeighting"    , "", 10, 0, 5);
  DiMuInvMassFakeWeightZoom_     = new TH1F("DiMuInvMassFakeWeightZoom"    , "", 600, 0, 30);
  TauVisMass_     = new TH1F("TauVisMass"    , "", 10, 0, 150);
  TauVisMassZoom_     = new TH1F("TauVisMassZoom"    , "", 10, 0, 30);
}

// ------------ method called once each job just after ending the event loop  ------------
void FakeRateMiniAODMCGetRegionA::endJob()
{
  //Make the Canvases
  TCanvas NEventsCanvas_("NEventsCanvas","",600,600);
  TCanvas InvMassTauMuMu1Canvas_("InvMassTauMuMu1Canvas","",600,600);
  TCanvas InvMassTauMuMu2Canvas_("InvMassTauMuMu2Canvas","",600,600);
  TCanvas DiMuInvMassFakeWeightCanvas_("DiMuInvMassFakeWeightCanvas","",600,600);
  TCanvas DiMuDiTauInvMassFakeWeightCanvas_("DiMuDiTauInvMassFakeWeightCanvas","",600,600);
  TCanvas DiTauInvMassFakeWeightCanvas_("DiTauInvMassFakeWeightCanvas","",600,600);
  TCanvas PtMu1FakeWeightCanvas_("PtMu1FakeWeightCanvas","",600,600);
  TCanvas PtMu2FakeWeightCanvas_("PtMu2FakeWeightCanvas","",600,600);
  TCanvas EtaFakeWeightCanvas_("EtaFakeWeightCanvas","",600,600);
  TCanvas DRFakeWeightCanvas_("DRFakeWeightCanvas","",600,600);
  TCanvas DRNoWeightingCanvas_("DRNoWeightingCanvas","",600,600);
  TCanvas DiMuInvMassFakeWeightZoomCanvas_("DiMuInvMassFakeWeightZoomCanvas","",600,600);
  TCanvas TauVisMassCanvas_("TauVisMassCanvas","",600,600);
  TCanvas TauVisMassZoomCanvas_("TauVisMassZoomCanvas","",600,600);


std::cout << "<----------------Declared Canvases-------------->" << std::endl;

  //Format the 1D plots and Draw (canvas, hist, grid, log y, log z, color, size, style, xAxisTitle, xTitleSize, xLabelSize, xTitleOffSet, yAxisTitle, yTitleSize, yLabelSize, yTitleOffset)
  VariousFunctions::formatAndDrawCanvasAndHist1D(NEventsCanvas_, NEvents_,
	 1, 0, 0, kBlack, 7, 20, "", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassTauMuMu1Canvas_, InvMassTauMuMu1_,
	 1, 0, 0, kBlack, 7, 20, "M(#tau_{#mu} #mu_{1})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassTauMuMu2Canvas_, InvMassTauMuMu2_,
	 1, 0, 0, kBlack, 7, 20, "M(#tau_{#mu} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiMuInvMassFakeWeightCanvas_, DiMuInvMassFakeWeight_,
         1, 0, 0, kBlack, 1, 20, "M(#mu_{2} #mu_{1})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiMuDiTauInvMassFakeWeightCanvas_, DiMuDiTauInvMassFakeWeight_,
         1, 0, 0, kBlack, 1, 20, "M(#mu_{2}#mu_{1}#mu_{3}#tau_{H})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiTauInvMassFakeWeightCanvas_, DiTauInvMassFakeWeight_,
         1, 0, 0, kBlack, 1, 20, "M(#mu_{3} #tau_{H})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(PtMu1FakeWeightCanvas_, PtMu1FakeWeight_,
         1, 0, 0, kBlack, 1, 20, "p_{T}(#mu_{1})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(PtMu2FakeWeightCanvas_, PtMu2FakeWeight_,
         1, 0, 0, kBlack, 1, 20, "p_{T}(#mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(EtaFakeWeightCanvas_, EtaFakeWeight_,
         1, 0, 0, kBlack, 1, 20, "#eta(#mu_{1})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DRFakeWeightCanvas_, DRFakeWeight_,
         1, 0, 0, kBlack, 1, 20, "#DeltaR(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DRNoWeightingCanvas_, DRNoWeighting_,
         1, 0, 0, kBlack, 1, 20, "#DeltaR(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiMuInvMassFakeWeightZoomCanvas_, DiMuInvMassFakeWeightZoom_,
         1, 0, 0, kBlack, 1, 20, "M(#mu_{2} #mu_{1})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauVisMassCanvas_, TauVisMass_,
         1, 0, 0, kBlack, 1, 20, "M(#tau_{2} #tau_{1})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauVisMassZoomCanvas_, TauVisMassZoom_,
         1, 0, 0, kBlack, 1, 20, "M(#tau_{2} #tau_{1})", .04, .04, 1.1,  "", .04, .04, 1.0, false);

std::cout << "after formatting" << std::endl;
  
////////////////////////
// For Pt Gen
////////////////////////
std::cout << "<----------------Formatted Canvases and Histos-------------->" << std::endl;

  //Write output file
  out_->cd();

  NEvents_->Write();
  InvMassTauMuMu1_->Write();
  InvMassTauMuMu2_->Write();
  DiMuInvMassFakeWeight_->Write();
  DiMuDiTauInvMassFakeWeight_->Write();
  DiTauInvMassFakeWeight_->Write();
  PtMu1FakeWeight_->Write();
  PtMu2FakeWeight_->Write();
  EtaFakeWeight_->Write();
  DRFakeWeight_->Write();
  DRNoWeighting_->Write();
  DiMuInvMassFakeWeightZoom_->Write();
  TauVisMass_->Write();
  TauVisMassZoom_->Write();

  out_->Write();
  out_->Close();
std::cout << "DONE" << std::endl;
}//EndJob

// ------------ method called when starting to processes a run  ------------
void FakeRateMiniAODMCGetRegionA::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void FakeRateMiniAODMCGetRegionA::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void FakeRateMiniAODMCGetRegionA::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void FakeRateMiniAODMCGetRegionA::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

//Delete Memory
void FakeRateMiniAODMCGetRegionA::reset(const bool doDelete)
{
  if ((doDelete) && (NEvents_ != NULL)) delete NEvents_;
  NEvents_ = NULL;
  if ((doDelete) && (InvMassTauMuMu1_ != NULL)) delete InvMassTauMuMu1_;
  InvMassTauMuMu1_ = NULL;
  if ((doDelete) && (InvMassTauMuMu2_ != NULL)) delete InvMassTauMuMu2_;
  InvMassTauMuMu2_ = NULL;
  if ((doDelete) && (DiMuInvMassFakeWeight_ != NULL)) delete DiMuInvMassFakeWeight_;
  DiMuInvMassFakeWeight_ = NULL;
  if ((doDelete) && (DiMuDiTauInvMassFakeWeight_ != NULL)) delete DiMuDiTauInvMassFakeWeight_;
  DiMuDiTauInvMassFakeWeight_ = NULL;
  if ((doDelete) && (DiTauInvMassFakeWeight_ != NULL)) delete DiTauInvMassFakeWeight_;
  DiTauInvMassFakeWeight_ = NULL;
  if ((doDelete) && (PtMu1FakeWeight_ != NULL)) delete PtMu1FakeWeight_;
  PtMu1FakeWeight_ = NULL;
  if ((doDelete) && (PtMu2FakeWeight_ != NULL)) delete PtMu2FakeWeight_;
  PtMu2FakeWeight_ = NULL;
  if ((doDelete) && (EtaFakeWeight_ != NULL)) delete EtaFakeWeight_;
  EtaFakeWeight_ = NULL;
  if ((doDelete) && (DRFakeWeight_ != NULL)) delete DRFakeWeight_;
  DRFakeWeight_ = NULL;
  if ((doDelete) && (DRNoWeighting_ != NULL)) delete DRNoWeighting_;
  DRNoWeighting_ = NULL;
  if ((doDelete) && (DiMuInvMassFakeWeightZoom_ != NULL)) delete DiMuInvMassFakeWeightZoom_;
  DiMuInvMassFakeWeightZoom_ = NULL;
  if ((doDelete) && (TauVisMass_ != NULL)) delete TauVisMass_;
  TauVisMass_ = NULL;
  if ((doDelete) && (TauVisMassZoom_ != NULL)) delete TauVisMassZoom_;
  TauVisMassZoom_ = NULL;


}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void FakeRateMiniAODMCGetRegionA::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FakeRateMiniAODMCGetRegionA);

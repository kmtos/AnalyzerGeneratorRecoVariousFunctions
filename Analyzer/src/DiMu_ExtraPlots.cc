// -*- C++ -*-
//
// Package:    DiMu_ExtraPlots
// Class:      DiMu_ExtraPlots
// 
/**\class DiMu_ExtraPlots DiMu_ExtraPlots.cc Analyzer/src/DiMu_ExtraPlots.cc

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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

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

class DiMu_ExtraPlots : public edm::EDAnalyzer {
   public:
      typedef reco::JetFloatAssociation::Container JetBCEnergyRatioCollection;
      explicit DiMu_ExtraPlots(const edm::ParameterSet&);
      ~DiMu_ExtraPlots();

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
      edm::EDGetTokenT<edm::View<pat::Muon> > mu12Tag_;
      edm::EDGetTokenT<edm::View<pat::Tau> > tauTag_;
      edm::EDGetTokenT<edm::View<pat::Muon> > mu3Tag_;
      edm::EDGetTokenT<edm::View<pat::MET> > metTag_;
      edm::EDGetTokenT<edm::View<pat::Jet> > jetTag_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryInfo_;
      edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
      double tauPtCut_;
      double xsec_; 
      double lumi_; 
      double summedWeights_;
      std::string PileupFileName_;

      //Histograms
      TH1F* NMedIsoTausWithMu3_;   
      TH1F* NEvents_;   
      TH1F* InvMassTauHadMu3_;
      TH1F* PtTauHadMu3_;
      TH1F* InvMassUpsilonRange_;
      TH1F* InvMassZPeakRange_;
      TH1F* InvMassFullRange_;
      TH1F* InvMassDiMuBarrBarr_;
      TH1F* InvMassDiMuBarrEndc_;
      TH1F* InvMassDiMuEndcEndc_;
      TH1F* Mu1Pt_;
      TH1F* Mu2Pt_;
      TH1F* DiMuPt_;
      TH1F* Mu1Eta_;
      TH1F* Mu2Eta_;
      TH1F* DiTauEta_;
      TH1F* DiTauPhi_;
      TH1F* DiMuEta_;
      TH1F* DiMudR_;
      TH1F* DiMuPhi_;
      TH1F* EtMET_;
      TH1F* DiTauDiMuMass_;
      TH1F* DiMuDiTauDeltaPhi_;
      TH1F* METDiTauDeltaPhi_;
      TH1F* METDiMuDeltaPhi_;
      TH2F* PtMu1vsPtMu2_;
      TH1F* PileupWeights_;
      TH1F* GenWeights_;
      TH1F* HighestCSVInclJet_;
      TH1F* HighestCombMVAJet_;
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
DiMu_ExtraPlots::DiMu_ExtraPlots(const edm::ParameterSet& iConfig):
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  mu12Tag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("mu12Tag"))),
  tauTag_(consumes<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("tauTag"))),
  mu3Tag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("mu3Tag"))),
  metTag_(consumes<edm::View<pat::MET> >(iConfig.getParameter<edm::InputTag>("metTag"))),
  jetTag_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetTag"))),
  pileupSummaryInfo_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummaryInfo"))),
  genEventInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfoToken"))),
  tauPtCut_(iConfig.getParameter<double>("tauPtCut")),
  xsec_(iConfig.getParameter<double>("xsec")),
  lumi_(iConfig.getParameter<double>("lumi")),
  summedWeights_(iConfig.getParameter<double>("summedWeights")),
  PileupFileName_(iConfig.getParameter<std::string>("PileupFileName"))
{
  reset(false);    
}//DiMu_ExtraPlots



DiMu_ExtraPlots::~DiMu_ExtraPlots()
{
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void DiMu_ExtraPlots::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "\n<------------THIS IS A NEW EVENT------------>" << std::endl;

  //Old Jet collection for bTagging
  edm::Handle<edm::View<pat::Muon> > pMu12;
  iEvent.getByToken(mu12Tag_, pMu12);

  edm::Handle<edm::View<pat::Tau> > pTaus;
  iEvent.getByToken(tauTag_, pTaus);

  edm::Handle<edm::View<pat::Muon> > pMu3;
  iEvent.getByToken(mu3Tag_, pMu3);

  edm::Handle<edm::View<pat::Jet> > pJets;
  iEvent.getByToken(jetTag_, pJets);

  double highestCSVInc = -1, highestCMVA = -1;
  for (edm::View<pat::Jet>::const_iterator iJet = pJets->begin(); iJet != pJets->end(); ++iJet)
  {
    if (iJet->pt() > tauPtCut_ && fabs(iJet->eta() ) < 2.4  && iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > highestCSVInc) 
      highestCSVInc = iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    if (iJet->pt() > tauPtCut_ && fabs(iJet->eta() ) < 2.4  && iJet->bDiscriminator("pfCombinedMVABJetTags") > highestCMVA) 
      highestCMVA = iJet->bDiscriminator("pfCombinedMVABJetTags");
  }//for iJet
  HighestCSVInclJet_->Fill(highestCSVInc );
  HighestCombMVAJet_->Fill(highestCMVA );
  edm::Handle<edm::View<pat::MET> > pMET;
  iEvent.getByToken(metTag_, pMET);
  pat::MET MET = pMET->at(0);

  edm::Handle<std::vector<PileupSummaryInfo> > pPileupSummaryInfo;
  iEvent.getByToken(pileupSummaryInfo_, pPileupSummaryInfo);

  int nTrueVertices = 0;
  if (pPileupSummaryInfo.isValid() && pPileupSummaryInfo->size()>0)
    nTrueVertices = pPileupSummaryInfo->at(1).getTrueNumInteractions();

  PileupFile = new TFile(PileupFileName_.c_str());
  TH1F* Pileup_ = (TH1F*)PileupFile->Get("PileupWeights");
  double pileupWeight = Pileup_->GetBinContent(nTrueVertices);
  PileupWeights_->Fill(pileupWeight);

  edm::Handle<GenEventInfoProduct> genEventInfo;
  iEvent.getByToken(genEventInfoToken_, genEventInfo);
  double eventGenWeight = genEventInfo->weight();
  double genWeight = eventGenWeight * lumi_ * xsec_ / summedWeights_;
  GenWeights_->Fill(genWeight);
  std::cout << "genWeight=" << genWeight << "\teventGenWeight=" << eventGenWeight << "\tlumi_=" << lumi_  << "\tsummedWeights_=" << summedWeights_ << "\txsec_=" << xsec_ << "\tpileupWeight=" <<pileupWeight <<  "\tnTrueVertices=" << nTrueVertices << std::endl;

  if (pPileupSummaryInfo.isValid()) std::cout << "VALID" << std::endl;
  if (pPileupSummaryInfo->size()>0) std::cout << "SIZE>0" << std::endl;
	
  pat::Muon mu1 = pat::Muon( (*pMu12)[0] );
  pat::Muon mu2 = pat::Muon( (*pMu12)[1] );
  std::cout << "mu1.pt()=" << mu1.pt() << "\tmu1.eta()=" << mu1.eta() << std::endl;

  reco::LeafCandidate::LorentzVector diMuP4 = mu1.p4() + mu2.p4();
  InvMassZPeakRange_->Fill(diMuP4.M(), pileupWeight*genWeight );
  InvMassUpsilonRange_->Fill(diMuP4.M(), pileupWeight*genWeight );
  InvMassFullRange_->Fill(diMuP4.M(), pileupWeight*genWeight );
  if (mu1.eta() < 1.6 && mu2.eta() < 1.6)
    InvMassDiMuBarrBarr_->Fill(diMuP4.M(), pileupWeight*genWeight );
  else if ( (mu1.eta() < 1.6 && mu2.eta() > 1.6) || (mu1.eta()  > 1.6 && mu2.eta() < 1.6) )
    InvMassDiMuBarrEndc_->Fill(diMuP4.M(), pileupWeight*genWeight );
  else if (mu1.eta() > 1.6 && mu2.eta() > 1.6)
    InvMassDiMuEndcEndc_->Fill(diMuP4.M(), pileupWeight*genWeight );
  Mu1Pt_->Fill(mu1.pt(), pileupWeight*genWeight );
  Mu2Pt_->Fill(mu2.pt(), pileupWeight*genWeight );
  PtMu1vsPtMu2_->Fill(mu1.pt(), mu2.pt() , pileupWeight*genWeight );
  DiMuPt_->Fill(diMuP4.Pt(), pileupWeight*genWeight );
  Mu1Eta_->Fill(mu1.eta(), pileupWeight*genWeight );
  Mu2Eta_->Fill(mu2.eta(), pileupWeight*genWeight );
  DiMuEta_->Fill(diMuP4.Eta(), pileupWeight*genWeight );
  DiMuPhi_->Fill(diMuP4.Phi(), pileupWeight*genWeight );
  DiMudR_->Fill(reco::deltaR(mu1, mu2), pileupWeight*genWeight );

  EtMET_->Fill(MET.pt(), pileupWeight*genWeight );
  METDiMuDeltaPhi_->Fill(fabs(diMuP4.Phi() - MET.phi() ), pileupWeight*genWeight );

  reco::LeafCandidate::LorentzVector diTauP4, fourBody;
  double bestdR = 100000000000000000;
  bool checkEventDiTau = false;
  unsigned int TauRemovedMuCount = 0;
  for (edm::View<pat::Tau>::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); ++iTau)
  {
    std::cout << "iTau=>pt=" << iTau->pt() << std::endl;
    bool checkMu3Removed = false;
    for (edm::View<pat::Muon>::const_iterator iMu = pMu3->begin(); iMu != pMu3->end(); ++iMu)
    {
      double currdR = deltaR(*iTau, *iMu);
      std::cout  << "\tdR=" << currdR << std::endl;
      if (currdR < .8 && currdR < bestdR)   
      {
        bestdR = currdR;
        diTauP4 = iTau->p4() + iMu->p4();
        fourBody = iTau->p4() + iMu->p4() + mu1.p4() + mu2.p4();
        checkMu3Removed = true;
        checkEventDiTau = true;
      }//if
    }//for iMu
    if (checkMu3Removed)
      TauRemovedMuCount++;
  }//for iTau
  if (checkEventDiTau)
  {
    DiTauDiMuMass_->Fill(fourBody.M(), pileupWeight*genWeight );
    DiMuDiTauDeltaPhi_->Fill(fabs(diTauP4.Phi() - diMuP4.Phi() ), pileupWeight*genWeight );
    InvMassTauHadMu3_->Fill(diTauP4.M(), pileupWeight*genWeight );
    PtTauHadMu3_->Fill(diTauP4.Pt(), pileupWeight*genWeight );
    DiMuDiTauDeltaPhi_->Fill(fabs(diTauP4.Phi() - MET.phi()), pileupWeight*genWeight );
    DiTauEta_->Fill(diTauP4.Eta(), pileupWeight*genWeight );
    DiTauPhi_->Fill(diTauP4.Phi(), pileupWeight*genWeight );
    METDiTauDeltaPhi_->Fill(fabs(diTauP4.Phi() - MET.phi() ), pileupWeight*genWeight );
  }//for checkMu3Removed 
  NMedIsoTausWithMu3_->Fill(TauRemovedMuCount ); 
}//End InvMass::analyze


// ------------ method called once each job just before starting event loop  ------------
void DiMu_ExtraPlots::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  //Open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");


  //Book histograms
  NMedIsoTausWithMu3_     = new TH1F("NMedIsoTausWithMu3"    , "", 9, -.5, 8.5);
  NEvents_     = new TH1F("NEvents"    , "", 9, -.5, 8.5);
      NEvents_->GetXaxis()->SetBinLabel(1, "TotalEvents"); 
      NEvents_->GetXaxis()->SetBinLabel(2, "#tau_{#mu} + #tau_{had} Match");
      NEvents_->GetXaxis()->SetBinLabel(3, "Gen #tau_{#mu} + #tau_{had}");
      NEvents_->GetXaxis()->SetBinLabel(4, "Gen Match #tau_{had}");
      NEvents_->GetXaxis()->SetBinLabel(5, "Gen Match #tau_{had}");
      NEvents_->GetXaxis()->SetBinLabel(6, "Event with #tau_{#mu} Removed");
      NEvents_->GetXaxis()->SetBinLabel(7, "Event with no #tau_{#mu} Removed ");
  InvMassTauHadMu3_     = new TH1F("InvMassTauHadMu3"    , "", 450, 0, 30);
  PtTauHadMu3_     = new TH1F("PtTauHadMu3"    , "", 500, 0, 500);
  InvMassUpsilonRange_     = new TH1F("InvMassUpsilonRange"    , "", 1000, 5, 15);
  InvMassZPeakRange_     = new TH1F("InvMassZPeakRange"    , "", 10000, 50, 150);
  InvMassFullRange_     = new TH1F("InvMassFullRange"    , "", 15000, 0, 150);
  InvMassDiMuBarrBarr_     = new TH1F("InvMassDiMuBarrBarr"    , "", 3000, 0, 30);
  InvMassDiMuBarrEndc_     = new TH1F("InvMassDiMuBarrEndc"    , "", 3000, 0, 30);
  InvMassDiMuEndcEndc_     = new TH1F("InvMassDiMuEndcEndc"    , "", 3000, 0, 30);
  Mu1Pt_     = new TH1F("Mu1Pt"    , "", 100, 0, 500);
  Mu2Pt_     = new TH1F("Mu2Pt"    , "", 100, 0, 500);
  DiMuPt_     = new TH1F("DiMuPt"    , "", 100, 0, 500);
  Mu1Eta_     = new TH1F("Mu1Eta"    , "", 75, -2.5, 2.5);
  Mu2Eta_     = new TH1F("Mu2Eta"    , "", 75, -2.5, 2.5);
  DiTauEta_     = new TH1F("DiTauEta"    , "", 75, -2.5, 2.5);
  DiTauPhi_     = new TH1F("DiTauPhi"    , "", 75, -3.2, 3.2);
  DiMuEta_     = new TH1F("DiMuEta"    , "", 75, -2.5, 2.5);
  DiMudR_     = new TH1F("DiMudR"    , "", 75, 0, 3.5);
  DiMuPhi_     = new TH1F("DiMuPhi"    , "", 75, -3.2, 3.2);
  EtMET_     = new TH1F("EtMET"    , "", 500, 0, 500);
  DiTauDiMuMass_     = new TH1F("DiTauDiMuMass"    , "", 100000, 0, 1000);
  DiMuDiTauDeltaPhi_     = new TH1F("DiMuDiTauDeltaPhi"    , "", 75, 0, 6.5);
  METDiTauDeltaPhi_     = new TH1F("METDiTauDeltaPhi"    , "", 75, 0, 6.5);
  METDiMuDeltaPhi_     = new TH1F("METDiMuDeltaPhi"    , "", 75, 0, 6.5);
  PtMu1vsPtMu2_  = new TH2F("PtMu1vsPtMu2" , "", 100, 0, 500, 100, 0, 500);
  PileupWeights_     = new TH1F("PileupWeights"    , "", 200, 0, 2);
  GenWeights_     = new TH1F("GenWeights"    , "", 20000, -10000, 10000);
  HighestCSVInclJet_     = new TH1F("HighestCSVInclJet"    , "", 100, 0, 1);
  HighestCombMVAJet_     = new TH1F("HighestCombMVAJet"    , "", 100, 0, 1);
}

// ------------ method called once each job just after ending the event loop  ------------
void DiMu_ExtraPlots::endJob()
{
  //Make the Canvases
  TCanvas NMedIsoTausWithMu3Canvas("NMedIsoTausWithMu3","",600,600);
  TCanvas NEventsCanvas("NEvents","",600,600);
  TCanvas InvMassTauHadMu3Canvas("InvMassTauHadMu3","",600,600);
  TCanvas PtTauHadMu3Canvas("PtTauHadMu3","",600,600);
  TCanvas InvMassUpsilonRangeCanvas("InvMassUpsilonRange","",600,600);
  TCanvas InvMassZPeakRangeCanvas("InvMassZPeakRange","",600,600);
  TCanvas InvMassFullRangeCanvas("InvMassFullRange","",600,600);
  TCanvas InvMassDiMuBarrBarrCanvas("InvMassDiMuBarrBarr","",600,600);
  TCanvas InvMassDiMuBarrEndcCanvas("InvMassDiMuBarrEndc","",600,600);
  TCanvas InvMassDiMuEndcEndcCanvas("InvMassDiMuEndcEndc","",600,600);
  TCanvas Mu1PtCanvas("Mu1Pt","",600,600);
  TCanvas Mu2PtCanvas("Mu2Pt","",600,600);
  TCanvas DiMuPtCanvas("DiMuPt","",600,600);
  TCanvas Mu1EtaCanvas("Mu1Eta","",600,600);
  TCanvas Mu2EtaCanvas("Mu2Eta","",600,600);
  TCanvas DiTauEtaCanvas("DiTauEta","",600,600);
  TCanvas DiTauPhiCanvas("DiTauPhi","",600,600);
  TCanvas DiMuEtaCanvas("DiMuEta","",600,600);
  TCanvas DiMudRCanvas("DiMudR","",600,600);
  TCanvas DiMuPhiCanvas("DiMuPhi","",600,600);
  TCanvas EtMETCanvas("EtMET","",600,600);
  TCanvas DiTauDiMuMassCanvas("DiTauDiMuMass","",600,600);
  TCanvas DiMuDiTauDeltaPhiCanvas("DiMuDiTauDeltaPhi","",600,600);
  TCanvas METDiTauDeltaPhiCanvas("METDiTauDeltaPhi","",600,600);
  TCanvas METDiMuDeltaPhiCanvas("METDiMuDeltaPhi","",600,600);
  TCanvas PtMu1vsPtMu2Canvas("PtMu1vsPtMu2","",600,600);
  TCanvas PileupWeightsCanvas("PileupWeights","",600,600);
  TCanvas GenWeightsCanvas("GenWeights","",600,600);
  TCanvas HighestCSVInclJetCanvas("HighestCSVInclJet","",600,600);
  TCanvas HighestCombMVAJetCanvas("HighestCombMVAJet","",600,600);

std::cout << "<----------------Declared Canvases-------------->" << std::endl;

  //Format the 1D plots and Draw (canvas, hist, grid, log y, log z, color, size, style, xAxisTitle, xTitleSize, xLabelSize, xTitleOffSet, yAxisTitle, yTitleSize, yLabelSize, yTitleOffset)
  VariousFunctions::formatAndDrawCanvasAndHist1D(NMedIsoTausWithMu3Canvas, NMedIsoTausWithMu3_,
	 1, 0, 0, kBlack, .1, 20, "# of MedIso #Taus with #tau_{#mu}", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NEventsCanvas, NEvents_,
	 1, 0, 0, kBlack, .1, 20, "", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassTauHadMu3Canvas, InvMassTauHadMu3_,
	 1, 0, 0, kBlack, .1, 20, "Mass(#tau_{H} #tau_{#mu})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(PtTauHadMu3Canvas, PtTauHadMu3_,
	 1, 0, 0, kBlack, .1, 20, "Pt(#tau_{H} #tau_{#mu})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassUpsilonRangeCanvas, InvMassUpsilonRange_,
	 1, 0, 0, kBlack, .1, 20, "Mass(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassZPeakRangeCanvas, InvMassZPeakRange_,
	 1, 0, 0, kBlack, .1, 20, "Mass(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassFullRangeCanvas, InvMassFullRange_,
	 1, 0, 0, kBlack, .1, 20, "Mass(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassDiMuBarrBarrCanvas, InvMassDiMuBarrBarr_,
         1, 0, 0, kBlack, .1, 20, "Barrel, Barrel Mass(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassDiMuBarrEndcCanvas, InvMassDiMuBarrEndc_,
         1, 0, 0, kBlack, .1, 20, "Barrel, Endcap Mass(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassDiMuEndcEndcCanvas, InvMassDiMuEndcEndc_,
         1, 0, 0, kBlack, .1, 20, "Endcap, Endcap Mass(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(Mu1PtCanvas, Mu1Pt_,
	 1, 0, 0, kBlack, .1, 20, "p_{T}(#mu_{1})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(Mu2PtCanvas, Mu2Pt_,
	 1, 0, 0, kBlack, .1, 20, "p_{T}(#mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiMuPtCanvas, DiMuPt_,
	 1, 0, 0, kBlack, .1, 20, "p_{T}(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(Mu1EtaCanvas, Mu1Eta_,
	 1, 0, 0, kBlack, .1, 20, "#eta(#mu_{1})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(Mu2EtaCanvas, Mu2Eta_,
	 1, 0, 0, kBlack, .1, 20, "#eta(#mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiTauEtaCanvas, DiTauEta_,
	 1, 0, 0, kBlack, .1, 20, "#eta(#tau#tau)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiTauPhiCanvas, DiTauPhi_,
	 1, 0, 0, kBlack, .1, 20, "#phi(#tau#tau)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiMuEtaCanvas, DiMuEta_,
	 1, 0, 0, kBlack, .1, 20, "#eta(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiMudRCanvas, DiMudR_,
	 1, 1, 0, kBlack, .1, 20, "#DeltaR(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiMuPhiCanvas, DiMuPhi_,
	 1, 1, 0, kBlack, .1, 20, "#Phi(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(EtMETCanvas, EtMET_,
	 1, 1, 0, kBlack, .1, 20, "MET E_{T}", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiTauDiMuMassCanvas, DiTauDiMuMass_,
         1, 0, 0, kBlack, .1, 20, "Mass(#mu#mu#tau#tau)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiMuDiTauDeltaPhiCanvas, DiMuDiTauDeltaPhi_,
         1, 1, 0, kBlack, .1, 20, "#DeltaPhi(#mu#mu, #tau#tau)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(METDiTauDeltaPhiCanvas, METDiTauDeltaPhi_,
         1, 1, 0, kBlack, .1, 20, "#DeltaPhi(MET, #tau#tau)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(METDiMuDeltaPhiCanvas, METDiMuDeltaPhi_,
         1, 1, 0, kBlack, .1, 20, "#Delta#Phi(MET, DiMu)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist2D(PtMu1vsPtMu2Canvas, PtMu1vsPtMu2_,
         1, 0, 0, kBlack, .1, 20, "p_{T}(#mu_{1})", .04, .04, 1.1, "p_{T}(#mu_{2})", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist1D(PileupWeightsCanvas, PileupWeights_,
         1, 0, 0, kBlack, .1, 20, "Pileup Weight (data/MC)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(GenWeightsCanvas, GenWeights_,
         1, 0, 0, kBlack, .1, 20, "Gen Weight (EventWeight * LumiData / LumiMC)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(HighestCSVInclJetCanvas, HighestCSVInclJet_,
         1, 0, 0, kBlack, .1, 20, "highest CSV Inclusive value jet in Event", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(HighestCombMVAJetCanvas, HighestCombMVAJet_,
         1, 0, 0, kBlack, .1, 20, "highest Combined MVA value jet in Event", .04, .04, 1.1,  "", .04, .04, 1.0, false);


std::cout << "after formatting" << std::endl;
  
////////////////////////
// For Pt Gen
////////////////////////
std::cout << "<----------------Formatted Canvases and Histos-------------->" << std::endl;

  //Write output file
  out_->cd();

  NMedIsoTausWithMu3Canvas.Write();
  NEventsCanvas.Write();
  InvMassTauHadMu3Canvas.Write();
  PtTauHadMu3Canvas.Write();
  InvMassUpsilonRangeCanvas.Write();
  InvMassZPeakRangeCanvas.Write();
  InvMassFullRangeCanvas.Write();
  InvMassDiMuBarrBarrCanvas.Write();
  InvMassDiMuBarrEndcCanvas.Write();
  InvMassDiMuEndcEndcCanvas.Write();
  Mu1PtCanvas.Write();
  Mu2PtCanvas.Write();
  DiMuPtCanvas.Write();
  Mu1EtaCanvas.Write();
  Mu2EtaCanvas.Write();
  DiTauEtaCanvas.Write();
  DiTauPhiCanvas.Write();
  DiMuEtaCanvas.Write();
  DiMudRCanvas.Write();
  DiMuPhiCanvas.Write();
  EtMETCanvas.Write();
  DiTauDiMuMassCanvas.Write();
  DiMuDiTauDeltaPhiCanvas.Write();
  METDiTauDeltaPhiCanvas.Write();
  METDiMuDeltaPhiCanvas.Write();
  PtMu1vsPtMu2Canvas.Write();
  PileupWeightsCanvas.Write();
  GenWeightsCanvas.Write();
  HighestCSVInclJetCanvas.Write();
  HighestCombMVAJetCanvas.Write();

  out_->Write();
  out_->Close();
std::cout << "DONE" << std::endl;
}//EndJob

// ------------ method called when starting to processes a run  ------------
void DiMu_ExtraPlots::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void DiMu_ExtraPlots::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DiMu_ExtraPlots::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DiMu_ExtraPlots::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

//Delete Memory
void DiMu_ExtraPlots::reset(const bool doDelete)
{
  if ((doDelete) && (NMedIsoTausWithMu3_ != NULL)) delete NMedIsoTausWithMu3_;
  NMedIsoTausWithMu3_ = NULL;
  if ((doDelete) && (NEvents_ != NULL)) delete NEvents_;
  NEvents_ = NULL;
  if ((doDelete) && (InvMassTauHadMu3_ != NULL)) delete InvMassTauHadMu3_;
  InvMassTauHadMu3_ = NULL;
  if ((doDelete) && (PtTauHadMu3_ != NULL)) delete PtTauHadMu3_;
  PtTauHadMu3_ = NULL;
  if ((doDelete) && (InvMassUpsilonRange_ != NULL)) delete InvMassUpsilonRange_;
  InvMassUpsilonRange_ = NULL;
  if ((doDelete) && (InvMassZPeakRange_ != NULL)) delete InvMassZPeakRange_;
  InvMassZPeakRange_ = NULL;
  if ((doDelete) && (InvMassFullRange_ != NULL)) delete InvMassFullRange_;
  InvMassFullRange_ = NULL;
 if ((doDelete) && (InvMassDiMuBarrBarr_ != NULL)) delete InvMassDiMuBarrBarr_;
 InvMassDiMuBarrBarr_ = NULL;
 if ((doDelete) && (InvMassDiMuBarrEndc_ != NULL)) delete InvMassDiMuBarrEndc_;
 InvMassDiMuBarrEndc_ = NULL;
 if ((doDelete) && (InvMassDiMuEndcEndc_ != NULL)) delete InvMassDiMuEndcEndc_;
 InvMassDiMuEndcEndc_ = NULL;
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
  if ((doDelete) && (DiTauEta_ != NULL)) delete DiTauEta_;
  DiTauEta_ = NULL;
  if ((doDelete) && (DiTauPhi_ != NULL)) delete DiTauPhi_;
  DiTauPhi_ = NULL;
  if ((doDelete) && (DiMuEta_ != NULL)) delete DiMuEta_;
  DiMuEta_ = NULL;
  if ((doDelete) && (DiMudR_ != NULL)) delete DiMudR_;
  DiMudR_ = NULL;
  if ((doDelete) && (DiMuPhi_ != NULL)) delete DiMuPhi_;
  DiMuPhi_ = NULL;
  if ((doDelete) && (EtMET_ != NULL)) delete EtMET_;
  EtMET_ = NULL;
  if ((doDelete) && (DiTauDiMuMass_ != NULL)) delete DiTauDiMuMass_;
  DiTauDiMuMass_ = NULL;
  if ((doDelete) && (DiMuDiTauDeltaPhi_ != NULL)) delete DiMuDiTauDeltaPhi_;
  DiMuDiTauDeltaPhi_ = NULL;
  if ((doDelete) && (METDiTauDeltaPhi_ != NULL)) delete METDiTauDeltaPhi_;
  METDiTauDeltaPhi_ = NULL;
  if ((doDelete) && (METDiMuDeltaPhi_ != NULL)) delete METDiMuDeltaPhi_;
  METDiMuDeltaPhi_ = NULL;
  if ((doDelete) && (PtMu1vsPtMu2_ != NULL)) delete PtMu1vsPtMu2_;
  PtMu1vsPtMu2_ = NULL;
  if ((doDelete) && (PileupWeights_ != NULL)) delete PileupWeights_;
  PileupWeights_ = NULL;
  if ((doDelete) && (GenWeights_ != NULL)) delete GenWeights_;
  GenWeights_ = NULL;
  if ((doDelete) && (HighestCSVInclJet_ != NULL)) delete HighestCSVInclJet_;
  HighestCSVInclJet_ = NULL;
  if ((doDelete) && (HighestCombMVAJet_ != NULL)) delete HighestCombMVAJet_;
  HighestCombMVAJet_ = NULL;

}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiMu_ExtraPlots::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMu_ExtraPlots);

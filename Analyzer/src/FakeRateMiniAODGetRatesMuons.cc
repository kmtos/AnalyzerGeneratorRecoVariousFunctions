// -*- C++ -*-
//
// Package:    FakeRateMiniAODGetRatesMuons
// Class:      FakeRateMiniAODGetRatesMuons
// 
/**\class FakeRateMiniAODGetRatesMuons FakeRateMiniAODGetRatesMuons.cc Analyzer/src/FakeRateMiniAODGetRatesMuons.cc

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
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
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

class FakeRateMiniAODGetRatesMuons : public edm::EDAnalyzer {
   public:
      typedef reco::JetFloatAssociation::Container JetBCEnergyRatioCollection;
      explicit FakeRateMiniAODGetRatesMuons(const edm::ParameterSet&);
      ~FakeRateMiniAODGetRatesMuons();

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
      edm::EDGetTokenT<edm::View<pat::Muon> > mu1Tag_;
      edm::EDGetTokenT<edm::View<pat::Muon> > muonsTag_;
      edm::EDGetTokenT<edm::View<pat::Tau> > tauTag_;
      std::string tauIsoTag_;
      std::string decayModeFindingTag_;
      bool checkTau_;
      bool checkTauIso_;
      bool passTauIso_;
      double mu3dRMin_;
      double mu3dRMax_;
      double tauPtCut_;
      double mu3dROverlapCut_;
      bool requireRemovedMuon_;
      bool checkInvMass_;
      double checkInvMassMin_;
      double checkInvMassMax_;
      double relIsoCutVal_;
      bool passRelIso_;
      double mu12dRCut_;
      bool oppositeSign_;
      bool passdR_;
      double mu2PtCut_;
      bool isMC_;
      double xsec_;
      double lumi_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryInfo_;
      edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
      double summedWeights_;
      std::string PileupFileName_;

      //Histograms
      TH1F* NEvents_;   
      TH1F* InvMassTauMuMu1_;
      TH1F* InvMassMu1Mu2_;
      TH1F* InvMassTauMuMu2_;

      TH2F* EtavsPtMuonIso_;
      TH2F* EtavsPtMuonNoIso_;

      TH1F* MuonIsoEta_;
      TH1F* MuonNoIsoEta_;
      TGraphAsymmErrors* FinalEffMuonEta_;

      TH1F* MuonIsoPt_;
      TH1F* MuonNoIsoPt_;
      TGraphAsymmErrors* FinalEffMuonPt_;


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
FakeRateMiniAODGetRatesMuons::FakeRateMiniAODGetRatesMuons(const edm::ParameterSet& iConfig):
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  mu1Tag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("mu1Tag"))),
  muonsTag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonsTag"))),
  tauTag_(consumes<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("tauTag"))),
  tauIsoTag_(iConfig.getParameter<std::string>("tauIsoTag")),
  decayModeFindingTag_(iConfig.getParameter<std::string>("decayModeFindingTag")),
  checkTau_(iConfig.getParameter<bool>("checkTau")),
  checkTauIso_(iConfig.getParameter<bool>("checkTauIso")),
  passTauIso_(iConfig.getParameter<bool>("passTauIso")),
  mu3dRMin_(iConfig.getParameter<double>("mu3dRMin")),
  mu3dRMax_(iConfig.getParameter<double>("mu3dRMax")),
  tauPtCut_(iConfig.getParameter<double>("tauPtCut")),
  mu3dROverlapCut_(iConfig.getParameter<double>("mu3dROverlapCut")),
  requireRemovedMuon_(iConfig.getParameter<bool>("requireRemovedMuon")),
  checkInvMass_(iConfig.getParameter<bool>("checkInvMass")),
  checkInvMassMin_(iConfig.getParameter<double>("checkInvMassMin")),
  checkInvMassMax_(iConfig.getParameter<double>("checkInvMassMax")),
  relIsoCutVal_(iConfig.getParameter<double>("relIsoCutVal")),
  passRelIso_(iConfig.getParameter<bool>("passRelIso")),
  mu12dRCut_(iConfig.getParameter<double>("mu12dRCut")),
  oppositeSign_(iConfig.getParameter<bool>("oppositeSign")),
  passdR_(iConfig.existsAs<bool>("passdR")? iConfig.getParameter<bool>("passdR"):true),
  mu2PtCut_(iConfig.getParameter<double>("mu2PtCut")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  xsec_(iConfig.getParameter<double>("xsec")),
  lumi_(iConfig.getParameter<double>("lumi")),
  pileupSummaryInfo_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummaryInfo"))),
  genEventInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfoToken"))),
  summedWeights_(iConfig.getParameter<double>("summedWeights")),
  PileupFileName_(iConfig.getParameter<std::string>("PileupFileName"))

{
  reset(false);    
}//FakeRateMiniAODGetRatesMuons



FakeRateMiniAODGetRatesMuons::~FakeRateMiniAODGetRatesMuons()
{
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void FakeRateMiniAODGetRatesMuons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "\n<------------THIS IS A NEW EVENT------------>" << std::endl;
  NEvents_->Fill(0);

  //Get RECO Muons particle collection
  edm::Handle<edm::View<pat::Muon> > pMuons;
  iEvent.getByToken(muonsTag_, pMuons);

  //Get CleanJets Tau particle collection
  edm::Handle<edm::View<pat::Tau> > pTaus;
  iEvent.getByToken(tauTag_, pTaus);

  //Old Jet collection for bTagging
  edm::Handle<edm::View<pat::Muon> > pMu1;
  iEvent.getByToken(mu1Tag_, pMu1);

  pat::Muon mu1 = pMu1->at(0), mu2;

  double highestPt = -1;
  bool checkMu2 = false;
  for(edm::View<pat::Muon>::const_iterator iMuon=pMuons->begin(); iMuon!=pMuons->end();++iMuon)
  {
    if ( deltaR(*iMuon, mu1) < .001 && fabs(iMuon->pt() - mu1.pt()) < .01 )
      continue;
    if (iMuon->pt() < mu2PtCut_ )
      continue;
    reco::MuonPFIsolation iso = iMuon->pfIsolationR04();
    double dR = deltaR(mu2, mu1);
    double reliso = (iso.sumChargedHadronPt+TMath::Max(0.,iso.sumNeutralHadronEt+iso.sumPhotonEt-0.5*iso.sumPUPt) ) / iMuon->pt();
    if (((reliso < relIsoCutVal_ && passRelIso_) || (reliso > relIsoCutVal_ && !passRelIso_)) && iMuon->pt() > highestPt)
    {
      if ((dR < mu12dRCut_ && passdR_) || (dR > mu12dRCut_ && !passdR_))
      {
        if ((oppositeSign_ && mu1.pdgId() == (-1)*iMuon->pdgId() )  ||  (!oppositeSign_ && mu1.pdgId()==iMuon->pdgId() ))
        {
          mu2 = *iMuon;
          highestPt = iMuon->pt();
          checkMu2 = true;
        }//if
      }//if
    }//if iso
  }//for iMuon
  
  if (checkMu2)
    NEvents_->Fill(1);
  else
    return;

  ////////////////////////////////
  // Doing the DiMuInv Mass Check
  ////////////////////////////////
  reco::LeafCandidate::LorentzVector diMuP4 = mu1.p4() + mu2.p4();
  if (checkInvMass_ &&  (diMuP4.M() < checkInvMassMin_ || diMuP4.M() > checkInvMassMax_) )
    return;
  InvMassMu1Mu2_->Fill(diMuP4.M() );
  NEvents_->Fill(2);

  ///////////////////////////
  // iterating over the taus
  /////////////////////////// 
  if (checkTau_)
  {
    /////////////////
    // Finding a Mu3
    /////////////////
    bool checkMu3 = false;
    for (edm::View<pat::Muon>::const_iterator iMu = pMuons->begin(); iMu != pMuons->end(); ++iMu)
    { 
      double mu1dR = deltaR(mu1, *iMu), mu2dR = deltaR(mu2, *iMu);
      if (mu1dR > mu3dROverlapCut_ && mu2dR > mu3dROverlapCut_)
       checkMu3 = true;
    }//for iMu
    if (checkMu3)
      NEvents_->Fill(3);
    else
      return;

    reco::LeafCandidate::LorentzVector  diTauP4;
    pat::Muon mu3;
    bool checkDiTau = false;
    for (edm::View<pat::Tau>::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); ++iTau)
    {
      if (iTau->pt() < tauPtCut_ || fabs(iTau->eta() ) > 2.4)
        continue;
      if (iTau->tauID(decayModeFindingTag_) < .5)
        continue;
      if (checkTauIso_)
      {
        if ((iTau->tauID(tauIsoTag_)  < .5 && passTauIso_) || (iTau->tauID(tauIsoTag_)  > .5 && !passTauIso_) )
          continue;
      }//if checkTauIso
      double bestMu3dR = 10000000;
      bool checkSubMu = false;
      for (edm::View<pat::Muon>::const_iterator iMu = pMuons->begin(); iMu != pMuons->end(); ++iMu)
      {
        double currdR = deltaR(*iTau, *iMu), mu1dR = deltaR(mu1, *iMu), mu2dR = deltaR(mu2, *iMu);
        if (mu1dR < mu3dROverlapCut_ || mu2dR < mu3dROverlapCut_)
         continue;
        if (currdR < mu3dRMax_ && currdR > mu3dRMin_ && currdR < bestMu3dR)
        {
          diTauP4 = iTau->p4() + iMu->p4();
          checkSubMu = true;
          mu3 = *iMu;
        }//if
      }//for iMu
  
      if (checkSubMu || !requireRemovedMuon_)
        checkDiTau = true;
    }//for iTau
  
    if (checkDiTau)
      NEvents_->Fill(4);
    else
      return;
  
    reco::LeafCandidate::LorentzVector diMuP4_Mu1TauMu, diMuP4_Mu2TauMu;    
    diMuP4_Mu1TauMu = mu1.p4();
    diMuP4_Mu1TauMu += mu3.p4();
  
    diMuP4_Mu2TauMu = mu2.p4();
    diMuP4_Mu2TauMu += mu3.p4();
  
    InvMassTauMuMu1_->Fill(diMuP4_Mu1TauMu.M() );
    InvMassTauMuMu2_->Fill(diMuP4_Mu2TauMu.M() );
  
  }//if tau check 

  double pileupWeight = 1, genWeight = 1;  
  if (isMC_)
  {
    edm::Handle<std::vector<PileupSummaryInfo> > pPileupSummaryInfo;
    iEvent.getByToken(pileupSummaryInfo_, pPileupSummaryInfo);
  
    int nTrueVertices = 0;
    if (pPileupSummaryInfo.isValid() && pPileupSummaryInfo->size()>0)
      nTrueVertices = pPileupSummaryInfo->at(1).getTrueNumInteractions();
  
    PileupFile = new TFile(PileupFileName_.c_str());
    TH1F* Pileup_ = (TH1F*)PileupFile->Get("PileupWeights");
    pileupWeight = Pileup_->GetBinContent(nTrueVertices);
    PileupFile->Close();
  
    edm::Handle<GenEventInfoProduct> genEventInfo;
    iEvent.getByToken(genEventInfoToken_, genEventInfo);
    double eventGenWeight = genEventInfo->weight();
    genWeight = eventGenWeight * lumi_ * xsec_ / summedWeights_;
    std::cout << "genWeight=" << genWeight << "\teventGenWeight=" << eventGenWeight << "\tlumi_=" << lumi_  << "\tsummedWeights_=" << summedWeights_ << "\txsec_=" << xsec_ << "\tpileupWeight=" <<pileupWeight <<  "\tnTrueVertices=" << nTrueVertices << std::endl;

  }//if isMC
///////////////////////////////////////////////////////////////////////
// After checking the tau side of event, fill the mu fake rate histos
///////////////////////////////////////////////////////////////////////
  for(edm::View<pat::Muon>::const_iterator iMuon=pMuons->begin(); iMuon!=pMuons->end();++iMuon)
  {
    reco::MuonPFIsolation iso = iMuon->pfIsolationR04();
    double reliso = (iso.sumChargedHadronPt+TMath::Max(0.,iso.sumNeutralHadronEt+iso.sumPhotonEt-0.5*iso.sumPUPt) ) / iMuon->pt();
std::cout << "dR(iMu,mu1)= " << deltaR(*iMuon, mu1) << "\tdR(iMu,mu2)=" << deltaR(*iMuon, mu2) << std::endl;
    if ( deltaR(*iMuon, mu1) < mu3dROverlapCut_ || deltaR(*iMuon, mu2) < mu3dROverlapCut_)
      continue;
    MuonNoIsoPt_->Fill(iMuon->pt() , pileupWeight*genWeight );
    MuonNoIsoEta_->Fill(fabs( iMuon->eta()) , pileupWeight*genWeight );
    EtavsPtMuonNoIso_->Fill(iMuon->pt(), fabs( iMuon->eta()), pileupWeight*genWeight );
    if ((reliso < relIsoCutVal_ && passRelIso_) || (reliso > relIsoCutVal_ && !passRelIso_) )
    {
      MuonIsoPt_->Fill(iMuon->pt() , pileupWeight*genWeight );
      MuonIsoEta_->Fill(fabs( iMuon->eta()) , pileupWeight*genWeight );
      EtavsPtMuonIso_->Fill(iMuon->pt(), fabs( iMuon->eta()), pileupWeight*genWeight );
    }//if Iso check
  }//for

}//End FakeRateMiniAODGetRatesMuons::analyze


// ------------ method called once each job just before starting event loop  ------------
void FakeRateMiniAODGetRatesMuons::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  //Open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");


  //Book histograms
  NEvents_     = new TH1F("NEvents"    , "", 5, -.5, 4.5);
      NEvents_->GetXaxis()->SetBinLabel(1, "TotalEvents"); 
      NEvents_->GetXaxis()->SetBinLabel(2, "Mu2");
      NEvents_->GetXaxis()->SetBinLabel(3, "Mass Cut");
      NEvents_->GetXaxis()->SetBinLabel(4, "Mu3");
      NEvents_->GetXaxis()->SetBinLabel(5, "DiTau");
//      NEvents_->GetXaxis()->SetBinLabel(6, "Event with #tau_{#mu} Removed");
//      NEvents_->GetXaxis()->SetBinLabel(7, "Event with no #tau_{#mu} Removed ");
  InvMassTauMuMu1_     = new TH1F("InvMassTauMuMu1"    , "", 75, 0, 150);
  InvMassMu1Mu2_     = new TH1F("InvMassMu1Mu2"    , "", 75, 0, 150);
  InvMassTauMuMu2_     = new TH1F("InvMassTauMuMu2"    , "", 75, 0, 150);

  Float_t binsx[] = {5, 10, 25, 50, 100};
  Float_t binsy[] = {0, .8, 1.2, 2.4};
  EtavsPtMuonIso_  = new TH2F("EtavsPtMuonIso" , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx, sizeof(binsy)/sizeof(Float_t) - 1, binsy);
  EtavsPtMuonNoIso_  = new TH2F("EtavsPtMuonNoIso"     , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx, sizeof(binsy)/sizeof(Float_t) - 1, binsy);

  MuonIsoEta_  = new TH1F("MuonIsoEta"    , "", 11, -2.4, 2.4);
  MuonNoIsoEta_    = new TH1F("MuonNoIsoEta", "", 11, -2.4, 2.4);

  MuonIsoPt_  = new TH1F("MuonIsoPt"    , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx);
  MuonNoIsoPt_    = new TH1F("MuonNoIsoPt", "", sizeof(binsx)/sizeof(Float_t) - 1, binsx);

  FinalEffMuonEta_ = new TGraphAsymmErrors(11);
  FinalEffMuonPt_ = new TGraphAsymmErrors(11);

  InvMassTauMuMu1_->Sumw2();
  InvMassMu1Mu2_->Sumw2();
  InvMassTauMuMu2_->Sumw2();

  EtavsPtMuonIso_->Sumw2();
  EtavsPtMuonNoIso_->Sumw2();

  MuonIsoEta_->Sumw2();
  MuonNoIsoEta_->Sumw2();

  MuonIsoPt_->Sumw2();
  MuonNoIsoPt_->Sumw2();


}

// ------------ method called once each job just after ending the event loop  ------------
void FakeRateMiniAODGetRatesMuons::endJob()
{
  //Make the Canvases
  TCanvas NEventsCanvas("NEvents","",600,600);
  TCanvas InvMassTauMuMu1Canvas("InvMassTauMuMu1","",600,600);
  TCanvas InvMassMu1Mu2Canvas("InvMassMu1Mu2","",600,600);
  TCanvas InvMassTauMuMu2Canvas("InvMassTauMuMu2","",600,600);

  TCanvas EtavsPtMuonIsoCanvas("EtavsPtMuonIso","",600,600);
  TCanvas EtavsPtMuonNoIsoCanvas("EtavsPtMuonNoIso","",600,600);

  TCanvas MuonIsoEtaCanvas("MuonIsoEta","",600,600);
  TCanvas MuonNoIsoEtaCanvas("MuonNoIsoEta","",600,600);
  TCanvas FinalEffMuonEtaCanvas("FinalEffMuonEta","",600,600);

  TCanvas MuonIsoPtCanvas("MuonIsoPt","",600,600);
  TCanvas MuonNoIsoPtCanvas("MuonNoIsoPt","",600,600);
  TCanvas FinalEffMuonPtCanvas("FinalEffMuonPt","",600,600);

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


  VariousFunctions::formatAndDrawCanvasAndHist2D(EtavsPtMuonIsoCanvas, EtavsPtMuonIso_,
         1, 0, 0, kBlack, 7, 20, "p_{T} #mu Iso", .04, .04, 1.1, "|#eta| #mu Iso", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(EtavsPtMuonNoIsoCanvas, EtavsPtMuonNoIso_,
         1, 0, 0, kBlack, 7, 20, "p_{T} #mu Not Iso", .04, .04, 1.1, "|#eta| #mu Not Iso", .04, .04, 1.6, "", .04, .04, 1.0);

  VariousFunctions::formatAndDrawCanvasAndHist1D(MuonIsoEtaCanvas, MuonIsoEta_,
	 1, 0, 0, kBlack, 7, 20, "Eta(slimmedMuon + Iso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MuonNoIsoEtaCanvas, MuonNoIsoEta_,
	 1, 0, 0, kBlack, 7, 20, "Eta(slimmedMuon)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMuonEtaCanvas, FinalEffMuonEta_,
	 1, 1, 1, kBlack, 1, 20, "Eta(slimmedMuon + Iso / slimmedMuon)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);


  VariousFunctions::formatAndDrawCanvasAndHist1D(MuonIsoPtCanvas, MuonIsoPt_,
         1, 0, 0, kBlack, 7, 20, "Pt(slimmedMuon + Iso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MuonNoIsoPtCanvas, MuonNoIsoPt_,
         1, 0, 0, kBlack, 7, 20, "Pt(slimmedMuon)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMuonPtCanvas, FinalEffMuonPt_,
         1, 1, 1, kBlack, 1, 20, "Pt(slimmedMuon + Iso / slimmedMuon)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

std::cout << "after formatting" << std::endl;
  
////////////////////////
// For Pt Gen
////////////////////////
std::cout << "<----------------Formatted Canvases and Histos-------------->" << std::endl;

  //Write output file
  out_->cd();

  NEvents_->Write();
  InvMassTauMuMu1_->Write();
  InvMassMu1Mu2_->Write();
  InvMassTauMuMu2_->Write();
 
  EtavsPtMuonIso_->Write();
  EtavsPtMuonNoIso_->Write();

  MuonIsoEta_->Write();
  MuonNoIsoEta_->Write();
  FinalEffMuonEta_->Write();

  MuonIsoPt_->Write();
  MuonNoIsoPt_->Write();
  FinalEffMuonPt_->Write();


  out_->Write();
  out_->Close();
std::cout << "DONE" << std::endl;
}//EndJob

// ------------ method called when starting to processes a run  ------------
void FakeRateMiniAODGetRatesMuons::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void FakeRateMiniAODGetRatesMuons::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void FakeRateMiniAODGetRatesMuons::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void FakeRateMiniAODGetRatesMuons::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

//Delete Memory
void FakeRateMiniAODGetRatesMuons::reset(const bool doDelete)
{
  if ((doDelete) && (NEvents_ != NULL)) delete NEvents_;
  NEvents_ = NULL;
  if ((doDelete) && (InvMassTauMuMu1_ != NULL)) delete InvMassTauMuMu1_;
  InvMassTauMuMu1_ = NULL;
  if ((doDelete) && (InvMassMu1Mu2_ != NULL)) delete InvMassMu1Mu2_;
  InvMassMu1Mu2_ = NULL;
  if ((doDelete) && (InvMassTauMuMu2_ != NULL)) delete InvMassTauMuMu2_;
  InvMassTauMuMu2_ = NULL;

  if ((doDelete) && (EtavsPtMuonIso_ != NULL)) delete EtavsPtMuonIso_;
  EtavsPtMuonIso_ = NULL;
  if ((doDelete) && (EtavsPtMuonNoIso_ != NULL)) delete EtavsPtMuonNoIso_;
  EtavsPtMuonNoIso_ = NULL;

  if ((doDelete) && (MuonIsoEta_ != NULL)) delete MuonIsoEta_;
  MuonIsoEta_ = NULL;
  if ((doDelete) && (MuonNoIsoEta_ != NULL)) delete MuonNoIsoEta_;
  MuonNoIsoEta_ = NULL;
  if ((doDelete) && (FinalEffMuonEta_ != NULL)) delete FinalEffMuonEta_;
  FinalEffMuonEta_ = NULL;

  if ((doDelete) && (MuonIsoPt_ != NULL)) delete MuonIsoPt_;
  MuonIsoPt_ = NULL;
  if ((doDelete) && (MuonNoIsoPt_ != NULL)) delete MuonNoIsoPt_;
  MuonNoIsoPt_ = NULL;
  if ((doDelete) && (FinalEffMuonPt_ != NULL)) delete FinalEffMuonPt_;
  FinalEffMuonPt_ = NULL;


}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void FakeRateMiniAODGetRatesMuons::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FakeRateMiniAODGetRatesMuons);

// -*- C++ -*-
//
// Package:    FakeRateWithWeightsAnalyzer
// Class:      FakeRateWithWeightsAnalyzer
// 
/**\class FakeRateWithWeightsAnalyzer FakeRateWithWeightsAnalyzer.cc Analyzer/src/FakeRateWithWeightsAnalyzer.cc

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

class FakeRateWithWeightsAnalyzer : public edm::EDAnalyzer {
   public:
      typedef reco::JetFloatAssociation::Container JetBCEnergyRatioCollection;
      explicit FakeRateWithWeightsAnalyzer(const edm::ParameterSet&);
      ~FakeRateWithWeightsAnalyzer();

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
      TFile* TH2File;

      //name of output root file
      std::string outFileName_;
      edm::EDGetTokenT<vector<reco::PFJet> > akJetTag_;
      edm::EDGetTokenT<vector<reco::PFTau> > tauTag_;
      edm::EDGetTokenT<reco::PFTauDiscriminator>  looseIsoTag_;
      edm::EDGetTokenT<reco::PFTauDiscriminator>  medIsoTag_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> tightIsoTag_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> decayModeFindingTag_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> isoRawTag_;
      edm::EDGetTokenT<reco::PFJetCollection>  oldJetTag_;
      edm::EDGetTokenT<reco::JetTagCollection> csvBTag_;
      edm::EDGetTokenT<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > mu12Tag_;
      edm::EDGetTokenT<vector<reco::Muon> > muonsTag_;
      edm::EDGetTokenT<edm::ValueMap<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > > muonMapTag_;
      bool requireRemovedMuon_;
      edm::EDGetTokenT<MuonRefVector> muonSrc_;
      std::string TH2FileName_;
      bool checkInvMass_;
      double checkInvMassValue_;

      //Histograms
      TH1F* NEvents_;   
      TH1F* InvMassTauMuMu1_;
      TH1F* InvMassTauMuMu2_;
      TH1F* InvMassFakeWeight_;
      TH1F* PtMu1FakeWeight_;
      TH1F* PtMu2FakeWeight_;
      TH1F* EtaFakeWeight_;
      TH1F* DRFakeWeight_;
      TH1F* DRNoWeighting_;
      TH1F* InvMassFakeWeightZoom_;
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
FakeRateWithWeightsAnalyzer::FakeRateWithWeightsAnalyzer(const edm::ParameterSet& iConfig):
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  akJetTag_(consumes<vector<reco::PFJet> >(iConfig.getParameter<edm::InputTag>("akJetTag"))),
  tauTag_(consumes<vector<reco::PFTau> >(iConfig.getParameter<edm::InputTag>("tauTag"))),
  looseIsoTag_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("looseIsoTag"))),
  medIsoTag_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("medIsoTag"))),
  tightIsoTag_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("tightIsoTag"))),
  decayModeFindingTag_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("decayModeFindingTag"))),
  isoRawTag_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("isoRawTag"))),
  oldJetTag_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("oldJetTag"))),
  csvBTag_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("csvBTag"))),
  mu12Tag_(consumes<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > >(iConfig.getParameter<edm::InputTag>("mu12Tag"))),
  muonsTag_(consumes<vector<reco::Muon> >(iConfig.getParameter<edm::InputTag>("muonsTag"))),
  muonMapTag_(consumes<edm::ValueMap<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > >(iConfig.getParameter<edm::InputTag>("muonMapTag"))),
  requireRemovedMuon_(iConfig.getParameter<bool>("requireRemovedMuon")),
  muonSrc_(consumes<MuonRefVector>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
  TH2FileName_(iConfig.getParameter<std::string>("TH2FileName")),
  checkInvMass_(iConfig.getParameter<bool>("checkInvMass")),
  checkInvMassValue_(iConfig.getParameter<double>("checkInvMassValue"))
{
  reset(false);    
}//FakeRateWithWeightsAnalyzer



FakeRateWithWeightsAnalyzer::~FakeRateWithWeightsAnalyzer()
{
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void FakeRateWithWeightsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "\n<------------THIS IS A NEW EVENT------------>" << std::endl;
  NEvents_->Fill(0);

  //Get ak4Jets particle collection
  edm::Handle<std::vector<reco::PFJet> > pAkJets;
  iEvent.getByToken(akJetTag_, pAkJets);

  //Get CleanJets Tau particle collection
  edm::Handle<std::vector<reco::PFTau> > pTaus;
  iEvent.getByToken(tauTag_, pTaus);

  //Get Loose Iso Collection
  Handle<PFTauDiscriminator> pLooseIsoDisc; 
  iEvent.getByToken(looseIsoTag_, pLooseIsoDisc); 

  //Get Medium Iso Collection
  Handle<PFTauDiscriminator> pMedIsoDisc; 
  iEvent.getByToken(medIsoTag_, pMedIsoDisc);

  //Get Tight Iso Collection
  Handle<PFTauDiscriminator> pTightIsoDisc; 
  iEvent.getByToken(tightIsoTag_, pTightIsoDisc);

  //Get Decay Mode Finding Collection
  Handle<PFTauDiscriminator> pDMFinding; 
  iEvent.getByToken(decayModeFindingTag_, pDMFinding);

  //Get IsoRaw  Collection
  Handle<PFTauDiscriminator> pIsoRaw;
  iEvent.getByToken(isoRawTag_, pIsoRaw);

  //Old Jet collection for bTagging
  edm::Handle<reco::PFJetCollection> pOldJets;
  iEvent.getByToken(oldJetTag_, pOldJets);

  //Get combVertMVA JetTagCollection
  edm::Handle<reco::JetTagCollection> pCSV;
  iEvent.getByToken(csvBTag_, pCSV);

  //Old Jet collection for bTagging
  edm::Handle<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > pMu12;
  iEvent.getByToken(mu12Tag_, pMu12);
  reco::MuonRef mu1Ref = reco::MuonRef((*pMu12)[0] );
  reco::MuonRef mu2Ref = reco::MuonRef((*pMu12)[1] );
  reco::LeafCandidate::LorentzVector diMuP4 = mu1Ref->p4() + mu2Ref->p4();

  if (checkInvMass_ && checkInvMassValue_ > diMuP4.M() )
    return;

  //Get RECO Muons particle collection
  edm::Handle<std::vector<reco::Muon> > pMuons;
  iEvent.getByToken(muonsTag_, pMuons);
   
  //Get the Muon Refs
  Handle<MuonRefVector> muons;
  iEvent.getByToken(muonSrc_, muons);
  std::vector<unsigned int> muonRefKeys;
  for (MuonRefVector::const_iterator iMuon = muons->begin(); iMuon != muons->end(); ++iMuon)
    muonRefKeys.push_back(iMuon->key());

  //get jet-muon map
  edm::Handle<edm::ValueMap<reco::MuonRefVector> > pMuonMap;
  iEvent.getByToken(muonMapTag_, pMuonMap);


//////////////////////////////
// Begin Analyzer
//////////////////////////////
  for (std::vector<reco::PFTau>::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); ++iTau)
  {
    // Now Do some muon stuff and seeing if a muon was removed
    reco::PFTauRef PFTauRef(pTaus, iTau - pTaus->begin());
    const reco::PFJetRef& tauJetRef = (*iTau).jetRef();
    const reco::MuonRefVector& removedMuons = (*pMuonMap)[tauJetRef];
    std::vector<reco::MuonRef> removedMuonRefs;
    for (reco::MuonRefVector::const_iterator iMuon = removedMuons.begin(); iMuon != removedMuons.end(); ++iMuon) 
      removedMuonRefs.push_back(*iMuon);//for iMuon
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

    reco::MuonRef removedMuonRef;
    bool removedMu = false;
    for (unsigned int iter = 0; iter < removedMuonRefs.size(); iter++)
    { 
      double dPhi = reco::deltaPhi(removedMuonRefs[iter]->phi(), iTau->phi() );
      double dR_tauMu = sqrt( (removedMuonRefs[iter]->eta() - iTau->eta() ) * (removedMuonRefs[iter]->eta() - iTau->eta() )  +  dPhi * dPhi );
      std::cout << "\t\t\tMuRef->pt()= " << removedMuonRefs[iter]->pt() << "  \tdR_tauMu= " << dR_tauMu << std::endl;
      if (dR_tauMu < .5)
      {
        removedMu = true; 
        removedMuonRef = removedMuonRefs[iter];
        break;
      }//if
    }//for iter
        
    if ( (!removedMu && requireRemovedMuon_) || iTau->pt() < 20.0 || fabs(iTau->eta() ) > 2.4)
      continue;

    if (removedMu)
    {
      reco::LeafCandidate::LorentzVector diMuP4_1, diMuP4_2;    
      diMuP4_1 = mu1Ref->p4();
      diMuP4_1 += removedMuonRef->p4();
 
      diMuP4_2 = mu2Ref->p4();
      diMuP4_2 += removedMuonRef->p4();

      InvMassTauMuMu1_->Fill(diMuP4_1.M() );
      InvMassTauMuMu2_->Fill(diMuP4_2.M() );
    }//if removed Mu

    TH2File = new TFile(TH2FileName_.c_str());
    TCanvas* FinalFakeRateDMtoMedIsoOnlyEtavsPtCanvas = (TCanvas*)TH2File->Get("FinalFakeRateDMtoMedIsoOnlyEtavsPtCanvas");
    TH2F* FinalFakeRateDMtoMedIsoOnlyEtavsPt_ = (TH2F*)FinalFakeRateDMtoMedIsoOnlyEtavsPtCanvas->GetPrimitive("FakeRateDMFindEtavsPt");
    TAxis *xaxis = FinalFakeRateDMtoMedIsoOnlyEtavsPt_->GetXaxis();
    TAxis *yaxis = FinalFakeRateDMtoMedIsoOnlyEtavsPt_->GetYaxis();
    Int_t binx = xaxis->FindBin(iTau->pt()  );
    Int_t biny = yaxis->FindBin(iTau->eta() );
    double rate= FinalFakeRateDMtoMedIsoOnlyEtavsPt_->GetBinContent(binx, biny);
    while (rate == 0 && binx > 0)
    {
      std::cout << "\tBin(" << binx << "," << biny << ")=0. New bin is (" << binx-1 << "," << biny << ")";
      binx--;
      rate = FinalFakeRateDMtoMedIsoOnlyEtavsPt_->GetBinContent(binx, biny);
      std::cout << " with new rate=" << rate << std::endl;
    }
    double weight = rate / (1 - rate);

    InvMassFakeWeight_->Fill(diMuP4.M(),  weight);
    InvMassFakeWeightZoom_->Fill(diMuP4.M(),  weight);
    reco::LeafCandidate::LorentzVector diTauP4 =  iTau->p4() + removedMuonRef->p4();
    TauVisMass_->Fill(diTauP4.M(),  weight);
    TauVisMassZoom_->Fill(diTauP4.M(),  weight);
    PtMu1FakeWeight_->Fill(mu1Ref->pt(), weight);
    PtMu2FakeWeight_->Fill(mu2Ref->pt(), weight);
    EtaFakeWeight_->Fill(mu1Ref->eta(), weight);
    double dPhi = reco::deltaPhi(mu1Ref->phi(), mu2Ref->phi() );
    double dR_tauMu = sqrt( (mu1Ref->eta() - mu2Ref->eta() ) * (mu1Ref->eta() - mu2Ref->eta() ) +  dPhi * dPhi);
    DRFakeWeight_->Fill(dR_tauMu, weight);
    DRNoWeighting_->Fill(dR_tauMu);
  }//iTau
}//End FakeRateWithWeightsAnalyzer::analyze


// ------------ method called once each job just before starting event loop  ------------
void FakeRateWithWeightsAnalyzer::beginJob()
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
  InvMassFakeWeight_     = new TH1F("InvMassFakeWeight"    , "", 10, 0, 150);
  PtMu1FakeWeight_     = new TH1F("PtMu1FakeWeight"    , "", 10, 0, 300);
  PtMu2FakeWeight_     = new TH1F("PtMu2FakeWeight"    , "", 10, 0, 300);
  EtaFakeWeight_     = new TH1F("EtaFakeWeight"    , "", 10, -2.5, 2.5);
  DRFakeWeight_     = new TH1F("DRFakeWeight"    , "", 10, 0, 5);
  DRNoWeighting_     = new TH1F("DRNoWeighting"    , "", 10, 0, 5);
  InvMassFakeWeightZoom_     = new TH1F("InvMassFakeWeightZoom"    , "", 10, 0, 20);
  TauVisMass_     = new TH1F("TauVisMass"    , "", 10, 0, 150);
  TauVisMassZoom_     = new TH1F("TauVisMassZoom"    , "", 10, 0, 30);
}

// ------------ method called once each job just after ending the event loop  ------------
void FakeRateWithWeightsAnalyzer::endJob()
{
  //Make the Canvases
  TCanvas NEventsCanvas("NEvents","",600,600);
  TCanvas InvMassTauMuMu1Canvas("InvMassTauMuMu1","",600,600);
  TCanvas InvMassTauMuMu2Canvas("InvMassTauMuMu2","",600,600);
  TCanvas InvMassFakeWeightCanvas("InvMassFakeWeight","",600,600);
  TCanvas PtMu1FakeWeightCanvas("PtMu1FakeWeight","",600,600);
  TCanvas PtMu2FakeWeightCanvas("PtMu2FakeWeight","",600,600);
  TCanvas EtaFakeWeightCanvas("EtaFakeWeight","",600,600);
  TCanvas DRFakeWeightCanvas("DRFakeWeight","",600,600);
  TCanvas DRNoWeightingCanvas("DRNoWeighting","",600,600);
  TCanvas InvMassFakeWeightZoomCanvas("InvMassFakeWeightZoom","",600,600);
  TCanvas TauVisMassCanvas("TauVisMass","",600,600);
  TCanvas TauVisMassZoomCanvas("TauVisMassZoom","",600,600);


std::cout << "<----------------Declared Canvases-------------->" << std::endl;

  //Format the 1D plots and Draw (canvas, hist, grid, log y, log z, color, size, style, xAxisTitle, xTitleSize, xLabelSize, xTitleOffSet, yAxisTitle, yTitleSize, yLabelSize, yTitleOffset)
  VariousFunctions::formatAndDrawCanvasAndHist1D(NEventsCanvas, NEvents_,
	 1, 0, 0, kBlack, 7, 20, " Weighted", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassTauMuMu1Canvas, InvMassTauMuMu1_,
	 1, 0, 0, kBlack, 7, 20, "Mass(#tau_{#mu} #mu_{1}) Weighted", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassTauMuMu2Canvas, InvMassTauMuMu2_,
	 1, 0, 0, kBlack, 7, 20, "Mass(#tau_{#mu} #mu_{2}) Weighted", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassFakeWeightCanvas, InvMassFakeWeight_,
         1, 0, 0, kBlack, 1, 20, "Mass(#mu_{2} #mu_{1}) Weighted", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(PtMu1FakeWeightCanvas, PtMu1FakeWeight_,
         1, 0, 0, kBlack, 1, 20, "p_{T}(#mu_{1}) Weighted", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(PtMu2FakeWeightCanvas, PtMu2FakeWeight_,
         1, 0, 0, kBlack, 1, 20, "p_{T}(#mu_{2}) Weighted", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(EtaFakeWeightCanvas, EtaFakeWeight_,
         1, 0, 0, kBlack, 1, 20, "#eta(#mu_{1}) Weighted", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DRFakeWeightCanvas, DRFakeWeight_,
         1, 0, 0, kBlack, 1, 20, "#DeltaR(#mu_{1} #mu_{2}) Weighted", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DRNoWeightingCanvas, DRNoWeighting_,
         1, 0, 0, kBlack, 1, 20, "#DeltaR(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassFakeWeightZoomCanvas, InvMassFakeWeightZoom_,
         1, 0, 0, kBlack, 1, 20, "Mass(#mu_{2} #mu_{1}) Weighted", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauVisMassCanvas, TauVisMass_,
         1, 0, 0, kBlack, 1, 20, "Visible Mass(#tau_{2} #tau_{1}) Weighted", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauVisMassZoomCanvas, TauVisMassZoom_,
         1, 0, 0, kBlack, 1, 20, "Visible Mass(#tau_{2} #tau_{1}) Weighted", .04, .04, 1.1,  "", .04, .04, 1.0, false);

std::cout << "after formatting" << std::endl;
  
////////////////////////
// For Pt Gen
////////////////////////
std::cout << "<----------------Formatted Canvases and Histos-------------->" << std::endl;

  //Write output file
  out_->cd();

  NEventsCanvas.Write();
  InvMassTauMuMu1Canvas.Write();
  InvMassTauMuMu2Canvas.Write();
  InvMassFakeWeightCanvas.Write();
  PtMu1FakeWeightCanvas.Write();
  PtMu2FakeWeightCanvas.Write();
  EtaFakeWeightCanvas.Write();
  DRFakeWeightCanvas.Write();
  DRNoWeightingCanvas.Write();
  InvMassFakeWeightZoomCanvas.Write();
  TauVisMassCanvas.Write();
  TauVisMassZoomCanvas.Write();
 

  out_->Write();
  out_->Close();
std::cout << "DONE" << std::endl;
}//EndJob

// ------------ method called when starting to processes a run  ------------
void FakeRateWithWeightsAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void FakeRateWithWeightsAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void FakeRateWithWeightsAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void FakeRateWithWeightsAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

//Delete Memory
void FakeRateWithWeightsAnalyzer::reset(const bool doDelete)
{
  if ((doDelete) && (NEvents_ != NULL)) delete NEvents_;
  NEvents_ = NULL;
  if ((doDelete) && (InvMassTauMuMu1_ != NULL)) delete InvMassTauMuMu1_;
  InvMassTauMuMu1_ = NULL;
  if ((doDelete) && (InvMassTauMuMu2_ != NULL)) delete InvMassTauMuMu2_;
  InvMassTauMuMu2_ = NULL;
  if ((doDelete) && (InvMassFakeWeight_ != NULL)) delete InvMassFakeWeight_;
  InvMassFakeWeight_ = NULL;
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
  if ((doDelete) && (InvMassFakeWeightZoom_ != NULL)) delete InvMassFakeWeightZoom_;
  InvMassFakeWeightZoom_ = NULL;
  if ((doDelete) && (TauVisMass_ != NULL)) delete TauVisMass_;
  TauVisMass_ = NULL;
  if ((doDelete) && (TauVisMassZoom_ != NULL)) delete TauVisMassZoom_;
  TauVisMassZoom_ = NULL;


}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void FakeRateWithWeightsAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FakeRateWithWeightsAnalyzer);

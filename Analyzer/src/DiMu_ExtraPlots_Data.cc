// -*- C++ -*-
//
// Package:    DiMu_ExtraPlots_Data
// Class:      DiMu_ExtraPlots_Data
// 
/**\class DiMu_ExtraPlots_Data DiMu_ExtraPlots_Data.cc Analyzer/src/DiMu_ExtraPlots_Data.cc

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

class DiMu_ExtraPlots_Data : public edm::EDAnalyzer {
   public:
      typedef reco::JetFloatAssociation::Container JetBCEnergyRatioCollection;
      explicit DiMu_ExtraPlots_Data(const edm::ParameterSet&);
      ~DiMu_ExtraPlots_Data();

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
      double tauPtCut_;

      //Histograms
      TH1F* NEvents_;
      TH1F* NMedIsoTausWithMu3_;
      TH1F* InvMassMu1TauMu_;
      TH1F* InvMassMu2TauMu_;
      TH1F* InvMassTauHadMu3_;
      TH1F* PtTauHadMu3_;
      TH1F* InvMassUpsilonRange_;
      TH1F* InvMassZPeakRange_;
      TH1F* InvMassFullRange_;
      TH1F* InvMassDiMuBarrBarr_;
      TH1F* InvMassDiMuBarrEndc_;
      TH1F* InvMassDiMuEndcEndc_;
      TH1F* InvMassDiMuBarrOver_;
      TH1F* InvMassDiMuOverOver_;
      TH1F* InvMassDiMuEndcOver_;
      TH1F* Mu1Pt_;
      TH1F* Mu2Pt_;
      TH1F* DiMuPt_;
      TH1F* DiMuPtSmallerdR_;
      TH1F* Mu1Eta_;
      TH1F* Mu2Eta_;
      TH1F* DiTauEta_;
      TH1F* DiTaudR_;
      TH1F* DiTauPhi_;
      TH1F* DiMuEta_;
      TH1F* DiMudR_;
      TH1F* DiMuPhi_;
      TH1F* EtMET_;
      TH1F* DiTauDiMuMass_;
      TH1F* DiTauMassSmallerdR_;
      TH1F* DiMuDiTauDeltaPhi_;
      TH1F* METDiTauDeltaPhi_;
      TH1F* METDiMuDeltaPhi_;
      TH2F* PtMu1vsPtMu2_;
      TH1F* PileupWeights_;
      TH1F* GenWeights_;
      TH1F* HighestCSVInclJet_;
      TH1F* HighestCombMVAJet_;
      TH1F* ZMassdR_;
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
DiMu_ExtraPlots_Data::DiMu_ExtraPlots_Data(const edm::ParameterSet& iConfig):
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  mu12Tag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("mu12Tag"))),
  tauTag_(consumes<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("tauTag"))),
  mu3Tag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("mu3Tag"))),
  metTag_(consumes<edm::View<pat::MET> >(iConfig.getParameter<edm::InputTag>("metTag"))),
  jetTag_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetTag"))),
  tauPtCut_(iConfig.getParameter<double>("tauPtCut"))
{
  reset(false);    
}//DiMu_ExtraPlots_Data



DiMu_ExtraPlots_Data::~DiMu_ExtraPlots_Data()
{
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void DiMu_ExtraPlots_Data::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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


	
  pat::Muon mu1, mu2;
  double highestMuPt = -1;

  if (pMu12->size() != 2)
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

  std::cout << "mu1.pt()=" << mu1.pt() << "\tmu1.eta()=" << mu1.eta() << "\tmu1.phi()=" << mu1.phi() << std::endl;
  std::cout << "mu2.pt()=" << mu2.pt() << "\tmu2.eta()=" << mu2.eta() << "\tmu2.phi()=" << mu2.phi() << std::endl;
  reco::LeafCandidate::LorentzVector diMuP4 = mu1.p4() + mu2.p4();

  if (diMuP4.M() > 30 )
    return;

  reco::LeafCandidate::LorentzVector diTauP4, fourBody, mu1mu3, mu2mu3;
  double bestdR = 100000000000000000;
  bool checkEventDiTau = false;
  unsigned int TauRemovedMuCount = 0;
  pat::Muon mu3;
  std::cout << "pTaus->size()= " << pTaus->size() << "\tpMu3->size()= " << pMu3->size() << std::endl;
  for (edm::View<pat::Tau>::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); ++iTau)
  {
    if (iTau->pt() < tauPtCut_)
      continue;
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
        mu3 = *iMu;
        mu1mu3 = iMu->p4() + mu1.p4();
        mu2mu3 = iMu->p4() + mu2.p4();
        fourBody = iTau->p4() + iMu->p4() + mu1.p4() + mu2.p4();
        checkMu3Removed = true;
        checkEventDiTau = true;
      }//if
    }//for iMu
    if (checkMu3Removed)
      TauRemovedMuCount++;
  }//for iTau
  std::cout << "FoundDiTau=" << checkEventDiTau << std::endl;
  if (checkEventDiTau)
  {
    HighestCSVInclJet_->Fill(highestCSVInc );
    HighestCombMVAJet_->Fill(highestCMVA );
    DiTauDiMuMass_->Fill(fourBody.M() );
    if (bestdR < .4)
    {
      DiMuPtSmallerdR_->Fill(diMuP4.Pt() );
      DiTauMassSmallerdR_->Fill(diTauP4.M() );
    }//if
    DiMuDiTauDeltaPhi_->Fill(fabs(diTauP4.Phi() - diMuP4.Phi() ) );
    InvMassTauHadMu3_->Fill(diTauP4.M() );
    InvMassMu1TauMu_->Fill(mu1mu3.M() );
    if (mu1mu3.M() > 60.0 && mu1mu3.M() < 120 )
    {
      double dRZ = deltaR(mu1, mu3);
      ZMassdR_->Fill(dRZ );
    }//if 
    if (mu2mu3.M() > 60.0 && mu2mu3.M() < 120 )
    {
      double dRZ = deltaR(mu2, mu3);
      ZMassdR_->Fill(dRZ );
    }//if 
    InvMassMu2TauMu_->Fill(mu2mu3.M() );
    PtTauHadMu3_->Fill(diTauP4.Pt() );
    DiTauEta_->Fill(diTauP4.Eta() );
    DiTauPhi_->Fill(diTauP4.Phi() );
    METDiTauDeltaPhi_->Fill(fabs(diTauP4.Phi() - MET.phi() ) );
    InvMassZPeakRange_->Fill(diMuP4.M() );
    InvMassUpsilonRange_->Fill(diMuP4.M() );
    InvMassFullRange_->Fill(diMuP4.M() );
    if (mu1.eta() < .9 && mu2.eta() < .9)
      InvMassDiMuBarrBarr_->Fill(diMuP4.M() );
    else if ( (mu1.eta() < .9 && mu2.eta() < 1.2 && mu2.eta() > .9) || (mu1.eta() < 1.2 && mu1.eta() > .9 && mu2.eta() < .9) )
      InvMassDiMuBarrOver_->Fill(diMuP4.M() );
    else if ( (mu1.eta() < .9 && mu2.eta() > 1.2) || (mu1.eta() > 1.2 && mu2.eta() < .9) )
      InvMassDiMuBarrEndc_->Fill(diMuP4.M() );
    else if ( (mu1.eta() > 1.2 && mu2.eta() < 1.2 && mu2.eta() > .9) || (mu1.eta() < 1.2 && mu1.eta() > .9 && mu2.eta() > 1.2) )
      InvMassDiMuEndcOver_->Fill(diMuP4.M() );
    else if ( mu2.eta() < 1.2 && mu2.eta() > .9 && mu1.eta() < 1.2 && mu1.eta() > .9)
      InvMassDiMuOverOver_->Fill(diMuP4.M() );
    else if (mu1.eta() > 1.2 && mu2.eta() > 1.2)
      InvMassDiMuEndcEndc_->Fill(diMuP4.M() );
    Mu1Pt_->Fill(mu1.pt() );
    Mu2Pt_->Fill(mu2.pt() );
    PtMu1vsPtMu2_->Fill(mu1.pt(), mu2.pt()  );
    DiMuPt_->Fill(diMuP4.Pt() );
    Mu1Eta_->Fill(mu1.eta() );
    Mu2Eta_->Fill(mu2.eta() );
    DiMuEta_->Fill(diMuP4.Eta() );
    DiMuPhi_->Fill(diMuP4.Phi() );
    DiMudR_->Fill(reco::deltaR(mu1, mu2) );
    DiTaudR_->Fill(bestdR);
    EtMET_->Fill(MET.pt() );
    METDiMuDeltaPhi_->Fill(fabs(diMuP4.Phi() - MET.phi() ) );
    EtMET_->Fill(MET.pt() );
    METDiMuDeltaPhi_->Fill(fabs(diMuP4.Phi() - MET.phi() ) );
  }//for checkMu3Removed 
  NMedIsoTausWithMu3_->Fill(TauRemovedMuCount );


}//End InvMass::analyze


// ------------ method called once each job just before starting event loop  ------------
void DiMu_ExtraPlots_Data::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  //Open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");


  //Book histograms
  NEvents_     = new TH1F("NEvents"    , "", 5, -.5, 4.5);
      NEvents_->GetXaxis()->SetBinLabel(1, "TotalEvents");
      NEvents_->GetXaxis()->SetBinLabel(2, "M > 30");
      NEvents_->GetXaxis()->SetBinLabel(3, "Fail Di-Tau");
      NEvents_->GetXaxis()->SetBinLabel(4, "Pass Di-Tau");
//      NEvents_->GetXaxis()->SetBinLabel(5, "DiTau");
///      NEvents_->GetXaxis()->SetBinLabel(6, "Event with #tau_{#mu} Removed");
//      NEvents_->GetXaxis()->SetBinLabel(7, "Event with no #tau_{#mu} Removed ");
  NMedIsoTausWithMu3_     = new TH1F("NMedIsoTausWithMu3"    , "", 9, -.5, 8.5);
  InvMassMu1TauMu_     = new TH1F("InvMassMu1TauMu"    , "", 90, 0, 400);
  InvMassMu2TauMu_     = new TH1F("InvMassMu2TauMu"    , "", 90, 0, 400);
  InvMassTauHadMu3_     = new TH1F("InvMassTauHadMu3"    , "", 90, 0, 30);
  PtTauHadMu3_     = new TH1F("PtTauHadMu3"    , "", 50, 0, 500);
  InvMassUpsilonRange_     = new TH1F("InvMassUpsilonRange"    , "", 50, 5, 15);
  InvMassZPeakRange_     = new TH1F("InvMassZPeakRange"    , "", 500, 50, 150);
  InvMassFullRange_     = new TH1F("InvMassFullRange"    , "", 600, 0, 30);
  InvMassDiMuBarrBarr_     = new TH1F("InvMassDiMuBarrBarr"    , "", 101, 0, 101);
  InvMassDiMuBarrEndc_     = new TH1F("InvMassDiMuBarrEndc"    , "", 101, 0, 101);
  InvMassDiMuEndcEndc_     = new TH1F("InvMassDiMuEndcEndc"    , "", 101, 0, 101);
  InvMassDiMuBarrOver_     = new TH1F("InvMassDiMuBarrOver"    , "", 101, 0, 101);
  InvMassDiMuOverOver_     = new TH1F("InvMassDiMuOverOver"    , "", 101, 0, 101);
  InvMassDiMuEndcOver_     = new TH1F("InvMassDiMuEndcOver"    , "", 101, 0, 101);
  Mu1Pt_     = new TH1F("Mu1Pt"    , "", 50, 0, 500);
  Mu2Pt_     = new TH1F("Mu2Pt"    , "", 50, 0, 500);
  DiMuPt_     = new TH1F("DiMuPt"    , "", 30, 0, 500);
  DiMuPtSmallerdR_     = new TH1F("DiMuPtSmallerdR"    , "", 30, 0, 500);
  Mu1Eta_     = new TH1F("Mu1Eta"    , "", 30, -2.5, 2.5);
  Mu2Eta_     = new TH1F("Mu2Eta"    , "", 30, -2.5, 2.5);
  DiTauEta_     = new TH1F("DiTauEta"    , "", 30, -2.5, 2.5);
  DiTaudR_     = new TH1F("DiTaudR"    , "", 30, 0, 1.0);
  DiTauPhi_     = new TH1F("DiTauPhi"    , "", 30, -3.2, 3.2);
  DiMuEta_     = new TH1F("DiMuEta"    , "", 30, -2.5, 2.5);
  DiMudR_     = new TH1F("DiMudR"    , "", 30, 0, 3.5);
  DiMuPhi_     = new TH1F("DiMuPhi"    , "", 30, -3.2, 3.2);
  EtMET_     = new TH1F("EtMET"    , "", 50, 0, 500);
  DiTauDiMuMass_     = new TH1F("DiTauDiMuMass"    , "", 300, 0, 1000);
  DiTauMassSmallerdR_     = new TH1F("DiTauMassSmallerdR"    , "", 60, 0, 30);
  DiMuDiTauDeltaPhi_     = new TH1F("DiMuDiTauDeltaPhi"    , "", 30, 0, 6.5);
  METDiTauDeltaPhi_     = new TH1F("METDiTauDeltaPhi"    , "", 30, 0, 6.5);
  METDiMuDeltaPhi_     = new TH1F("METDiMuDeltaPhi"    , "", 30, 0, 6.5);
  PtMu1vsPtMu2_  = new TH2F("PtMu1vsPtMu2" , "", 50, 0, 500, 50, 0, 500);
  PileupWeights_     = new TH1F("PileupWeights"    , "", 200, 0, 2);
  GenWeights_     = new TH1F("GenWeights"    , "", 20000, -10000, 10000);
  HighestCSVInclJet_     = new TH1F("HighestCSVInclJet"    , "", 100, 0, 1);
  HighestCombMVAJet_     = new TH1F("HighestCombMVAJet"    , "", 100, 0, 1);
  ZMassdR_     = new TH1F("ZMassdR"    , "", 100, 0, 6);
}

// ------------ method called once each job just after ending the event loop  ------------
void DiMu_ExtraPlots_Data::endJob()
{
  //Make the Canvases
  TCanvas NMedIsoTausWithMu3Canvas_("NMedIsoTausWithMu3Canvas","",600,600);
  TCanvas InvMassMu1TauMuCanvas_("InvMassMu1TauMuCanvas","",600,600);
  TCanvas InvMassMu2TauMuCanvas_("InvMassMu2TauMuCanvas","",600,600);
  TCanvas InvMassTauHadMu3Canvas_("InvMassTauHadMu3Canvas","",600,600);
  TCanvas PtTauHadMu3Canvas_("PtTauHadMu3Canvas","",600,600);
  TCanvas InvMassUpsilonRangeCanvas_("InvMassUpsilonRangeCanvas","",600,600);
  TCanvas InvMassZPeakRangeCanvas_("InvMassZPeakRangeCanvas","",600,600);
  TCanvas InvMassFullRangeCanvas_("InvMassFullRangeCanvas","",600,600);
  TCanvas InvMassDiMuBarrBarrCanvas_("InvMassDiMuBarrBarrCanvas","",600,600);
  TCanvas InvMassDiMuBarrEndcCanvas_("InvMassDiMuBarrEndcCanvas","",600,600);
  TCanvas InvMassDiMuEndcEndcCanvas_("InvMassDiMuEndcEndcCanvas","",600,600);
  TCanvas InvMassDiMuBarrOverCanvas_("InvMassDiMuBarrOverCanvas","",600,600);
  TCanvas InvMassDiMuOverOverCanvas_("InvMassDiMuOverOverCanvas","",600,600);
  TCanvas InvMassDiMuEndcOverCanvas_("InvMassDiMuEndcOverCanvas","",600,600);
  TCanvas Mu1PtCanvas_("Mu1PtCanvas","",600,600);
  TCanvas Mu2PtCanvas_("Mu2PtCanvas","",600,600);
  TCanvas DiMuPtCanvas_("DiMuPtCanvas","",600,600);
  TCanvas DiMuPtSmallerdRCanvas_("DiMuPtSmallerdRCanvas","",600,600);
  TCanvas Mu1EtaCanvas_("Mu1EtaCanvas","",600,600);
  TCanvas Mu2EtaCanvas_("Mu2EtaCanvas","",600,600);
  TCanvas DiTauEtaCanvas_("DiTauEtaCanvas","",600,600);
  TCanvas DiTaudRCanvas_("DiTaudRCanvas","",600,600);
  TCanvas DiTauPhiCanvas_("DiTauPhiCanvas","",600,600);
  TCanvas DiMuEtaCanvas_("DiMuEtaCanvas","",600,600);
  TCanvas DiMudRCanvas_("DiMudRCanvas","",600,600);
  TCanvas DiMuPhiCanvas_("DiMuPhiCanvas","",600,600);
  TCanvas EtMETCanvas_("EtMETCanvas","",600,600);
  TCanvas DiTauDiMuMassCanvas_("DiTauDiMuMassCanvas","",600,600);
  TCanvas DiTauMassSmallerdRCanvas_("DiTauMassSmallerdRCanvas","",600,600);
  TCanvas DiMuDiTauDeltaPhiCanvas_("DiMuDiTauDeltaPhiCanvas","",600,600);
  TCanvas METDiTauDeltaPhiCanvas_("METDiTauDeltaPhiCanvas","",600,600);
  TCanvas METDiMuDeltaPhiCanvas_("METDiMuDeltaPhiCanvas","",600,600);
  TCanvas PtMu1vsPtMu2Canvas_("PtMu1vsPtMu2Canvas","",600,600);
  TCanvas PileupWeightsCanvas_("PileupWeightsCanvas","",600,600);
  TCanvas GenWeightsCanvas_("GenWeightsCanvas","",600,600);
  TCanvas HighestCSVInclJetCanvas_("HighestCSVInclJetCanvas","",600,600);
  TCanvas HighestCombMVAJetCanvas_("HighestCombMVAJetCanvas","",600,600);
  TCanvas ZMassdRCanvas_("ZMassdRCanvas","",600,600);

std::cout << "Setting Bin contents to zero if negative" << std::endl;

  for (int i=1; i < InvMassTauHadMu3_->GetXaxis()->GetNbins(); i++)
  {
    if (InvMassTauHadMu3_->GetBinContent(i) < 0)
      InvMassTauHadMu3_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < PtTauHadMu3_->GetXaxis()->GetNbins(); i++)
  {
    if (PtTauHadMu3_->GetBinContent(i) < 0)
      PtTauHadMu3_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < InvMassMu1TauMu_->GetXaxis()->GetNbins(); i++)
  {
    if (InvMassMu1TauMu_->GetBinContent(i) < 0)
      InvMassMu1TauMu_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < InvMassMu2TauMu_->GetXaxis()->GetNbins(); i++)
  {
    if (InvMassMu2TauMu_->GetBinContent(i) < 0)
      InvMassMu2TauMu_->SetBinContent(i, 0);
  }//for i

  for (int i=1; i < InvMassUpsilonRange_->GetXaxis()->GetNbins(); i++)
  {
    if (InvMassUpsilonRange_->GetBinContent(i) < 0)
      InvMassUpsilonRange_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < InvMassZPeakRange_->GetXaxis()->GetNbins(); i++)
  {
    if (InvMassZPeakRange_->GetBinContent(i) < 0)
      InvMassZPeakRange_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < InvMassFullRange_->GetXaxis()->GetNbins(); i++)
  {
    if (InvMassFullRange_->GetBinContent(i) < 0)
      InvMassFullRange_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < InvMassDiMuBarrBarr_->GetXaxis()->GetNbins(); i++)
  {
    if (InvMassDiMuBarrBarr_->GetBinContent(i) < 0)
      InvMassDiMuBarrBarr_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < InvMassDiMuBarrEndc_->GetXaxis()->GetNbins(); i++)
  {
    if (InvMassDiMuBarrEndc_->GetBinContent(i) < 0)
      InvMassDiMuBarrEndc_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < InvMassDiMuEndcEndc_->GetXaxis()->GetNbins(); i++)
  {
    if (InvMassDiMuEndcEndc_->GetBinContent(i) < 0)
      InvMassDiMuEndcEndc_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < InvMassDiMuBarrOver_->GetXaxis()->GetNbins(); i++)
  {
    if (InvMassDiMuBarrOver_->GetBinContent(i) < 0)
      InvMassDiMuBarrOver_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < InvMassDiMuOverOver_->GetXaxis()->GetNbins(); i++)
  {
    if (InvMassDiMuOverOver_->GetBinContent(i) < 0)
      InvMassDiMuOverOver_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < InvMassDiMuEndcOver_->GetXaxis()->GetNbins(); i++)
  {
    if (InvMassDiMuEndcOver_->GetBinContent(i) < 0)
      InvMassDiMuEndcOver_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < Mu1Pt_->GetXaxis()->GetNbins(); i++)
  {
    if (Mu1Pt_->GetBinContent(i) < 0)
      Mu1Pt_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < Mu2Pt_->GetXaxis()->GetNbins(); i++)
  {
    if (Mu2Pt_->GetBinContent(i) < 0)
      Mu2Pt_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < DiMuPt_->GetXaxis()->GetNbins(); i++)
  {
    if (DiMuPt_->GetBinContent(i) < 0)
      DiMuPt_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < DiMuPtSmallerdR_->GetXaxis()->GetNbins(); i++)
  {
    if (DiMuPtSmallerdR_->GetBinContent(i) < 0)
      DiMuPtSmallerdR_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < Mu1Eta_->GetXaxis()->GetNbins(); i++)
  {
    if (Mu1Eta_->GetBinContent(i) < 0)
      Mu1Eta_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < Mu2Eta_->GetXaxis()->GetNbins(); i++)
  {
    if (Mu2Eta_->GetBinContent(i) < 0)
      Mu2Eta_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < DiTauEta_->GetXaxis()->GetNbins(); i++)
  {
    if (DiTauEta_->GetBinContent(i) < 0)
      DiTauEta_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < DiTaudR_->GetXaxis()->GetNbins(); i++)
  {
    if (DiTaudR_->GetBinContent(i) < 0)
      DiTaudR_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < DiTauPhi_->GetXaxis()->GetNbins(); i++)
  {
    if (DiTauPhi_->GetBinContent(i) < 0)
      DiTauPhi_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < DiMuEta_->GetXaxis()->GetNbins(); i++)
  {
    if (DiMuEta_->GetBinContent(i) < 0)
      DiMuEta_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < DiMudR_->GetXaxis()->GetNbins(); i++)
  {
    if (DiMudR_->GetBinContent(i) < 0)
      DiMudR_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < DiMuPhi_->GetXaxis()->GetNbins(); i++)
  {
    if (DiMuPhi_->GetBinContent(i) < 0)
      DiMuPhi_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < EtMET_->GetXaxis()->GetNbins(); i++)
  {
    if (EtMET_->GetBinContent(i) < 0)
      EtMET_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < DiTauDiMuMass_->GetXaxis()->GetNbins(); i++)
  {
    if (DiTauDiMuMass_->GetBinContent(i) < 0)
      DiTauDiMuMass_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < DiTauMassSmallerdR_->GetXaxis()->GetNbins(); i++)
  {
    if (DiTauMassSmallerdR_->GetBinContent(i) < 0)
      DiTauMassSmallerdR_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < DiMuDiTauDeltaPhi_->GetXaxis()->GetNbins(); i++)
  {
    if (DiMuDiTauDeltaPhi_->GetBinContent(i) < 0)
      DiMuDiTauDeltaPhi_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < METDiTauDeltaPhi_->GetXaxis()->GetNbins(); i++)
  {
    if (METDiTauDeltaPhi_->GetBinContent(i) < 0)
      METDiTauDeltaPhi_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < METDiMuDeltaPhi_->GetXaxis()->GetNbins(); i++)
  {
    if (METDiMuDeltaPhi_->GetBinContent(i) < 0)
      METDiMuDeltaPhi_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < PtMu1vsPtMu2_->GetXaxis()->GetNbins(); i++)
  {
    if (PtMu1vsPtMu2_->GetBinContent(i) < 0)
      PtMu1vsPtMu2_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < PileupWeights_->GetXaxis()->GetNbins(); i++)
  {
    if (PileupWeights_->GetBinContent(i) < 0)
      PileupWeights_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < GenWeights_->GetXaxis()->GetNbins(); i++)
  {
    if (GenWeights_->GetBinContent(i) < 0)
      GenWeights_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < HighestCSVInclJet_->GetXaxis()->GetNbins(); i++)
  {
    if (HighestCSVInclJet_->GetBinContent(i) < 0)
      HighestCSVInclJet_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < HighestCombMVAJet_->GetXaxis()->GetNbins(); i++)
  {
    if (HighestCombMVAJet_->GetBinContent(i) < 0)
      HighestCombMVAJet_->SetBinContent(i, 0);
  }//for i
  for (int i=1; i < ZMassdR_->GetXaxis()->GetNbins(); i++)
  {
    if (ZMassdR_->GetBinContent(i) < 0)
      ZMassdR_->SetBinContent(i, 0);
  }//for i



std::cout << "<----------------Declared Canvases-------------->" << std::endl;

  //Format the 1D plots and Draw (canvas, hist, grid, log y, log z, color, size, style, xAxisTitle, xTitleSize, xLabelSize, xTitleOffSet, yAxisTitle, yTitleSize, yLabelSize, yTitleOffset)
  VariousFunctions::formatAndDrawCanvasAndHist1D(NMedIsoTausWithMu3Canvas_, NMedIsoTausWithMu3_,
	 1, 0, 0, kBlack, .1, 20, "# of MedIso #Taus with #tau_{#mu}", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassMu1TauMuCanvas_, InvMassMu1TauMu_,
	 1, 0, 0, kBlack, .1, 20, "InvMass(#mu_{1}, #tau_{#mu})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassMu2TauMuCanvas_, InvMassMu2TauMu_,
	 1, 0, 0, kBlack, .1, 20, "InvMass(#mu_{2}, #tau_{#mu})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassTauHadMu3Canvas_, InvMassTauHadMu3_,
	 1, 0, 0, kBlack, .1, 20, "Mass(#tau_{H} #tau_{#mu})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(PtTauHadMu3Canvas_, PtTauHadMu3_,
	 1, 0, 0, kBlack, .1, 20, "Pt(#tau_{H} #tau_{#mu})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassUpsilonRangeCanvas_, InvMassUpsilonRange_,
	 1, 0, 0, kBlack, .1, 20, "Mass(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassZPeakRangeCanvas_, InvMassZPeakRange_,
	 1, 0, 0, kBlack, .1, 20, "Mass(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassFullRangeCanvas_, InvMassFullRange_,
	 1, 0, 0, kBlack, .1, 20, "Mass(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassDiMuBarrBarrCanvas_, InvMassDiMuBarrBarr_,
         1, 0, 0, kBlack, .1, 20, "BB Mass(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassDiMuBarrEndcCanvas_, InvMassDiMuBarrEndc_,
         1, 0, 0, kBlack, .1, 20, "BE Mass(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassDiMuEndcEndcCanvas_, InvMassDiMuEndcEndc_,
         1, 0, 0, kBlack, .1, 20, "EE Mass(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassDiMuBarrOverCanvas_, InvMassDiMuBarrOver_,
         1, 0, 0, kBlack, .1, 20, "BO Mass(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassDiMuOverOverCanvas_, InvMassDiMuOverOver_,
         1, 0, 0, kBlack, .1, 20, "OO Mass(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(InvMassDiMuEndcOverCanvas_, InvMassDiMuEndcOver_,
         1, 0, 0, kBlack, .1, 20, "EO Mass(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(Mu1PtCanvas_, Mu1Pt_,
	 1, 0, 0, kBlack, .1, 20, "p_{T}(#mu_{1})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(Mu2PtCanvas_, Mu2Pt_,
	 1, 0, 0, kBlack, .1, 20, "p_{T}(#mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiMuPtCanvas_, DiMuPt_,
	 1, 0, 0, kBlack, .1, 20, "p_{T}(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiMuPtSmallerdRCanvas_, DiMuPtSmallerdR_,
	 1, 0, 0, kBlack, .1, 20, "#DeltaR(#tau_{H}#tau_{#mu}) < .4p_{T}(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(Mu1EtaCanvas_, Mu1Eta_,
	 1, 0, 0, kBlack, .1, 20, "#eta(#mu_{1})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(Mu2EtaCanvas_, Mu2Eta_,
	 1, 0, 0, kBlack, .1, 20, "#eta(#mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiTauEtaCanvas_, DiTauEta_,
	 1, 0, 0, kBlack, .1, 20, "#eta(#tau#tau)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiTaudRCanvas_, DiTaudR_,
	 1, 1, 0, kBlack, .1, 20, "#DeltaR(#tau_{H} #tau_{#mu})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiTauPhiCanvas_, DiTauPhi_,
	 1, 0, 0, kBlack, .1, 20, "#phi(#tau#tau)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiMuEtaCanvas_, DiMuEta_,
	 1, 0, 0, kBlack, .1, 20, "#eta(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiMudRCanvas_, DiMudR_,
	 1, 1, 0, kBlack, .1, 20, "#DeltaR(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiMuPhiCanvas_, DiMuPhi_,
	 1, 1, 0, kBlack, .1, 20, "#Phi(#mu_{1} #mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(EtMETCanvas_, EtMET_,
	 1, 1, 0, kBlack, .1, 20, "MET E_{T}", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiTauDiMuMassCanvas_, DiTauDiMuMass_,
         1, 0, 0, kBlack, .1, 20, "Mass(#mu#mu#tau#tau)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiTauMassSmallerdRCanvas_, DiTauMassSmallerdR_,
         1, 0, 0, kBlack, .1, 20, "Mass(#tau#tau) #DeltaR(#tau_{H}#tau_{#mu}) < .4", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DiMuDiTauDeltaPhiCanvas_, DiMuDiTauDeltaPhi_,
         1, 1, 0, kBlack, .1, 20, "#DeltaPhi(#mu#mu, #tau#tau)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(METDiTauDeltaPhiCanvas_, METDiTauDeltaPhi_,
         1, 1, 0, kBlack, .1, 20, "#DeltaPhi(MET, #tau#tau)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(METDiMuDeltaPhiCanvas_, METDiMuDeltaPhi_,
         1, 1, 0, kBlack, .1, 20, "#Delta#Phi(MET, DiMu)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist2D(PtMu1vsPtMu2Canvas_, PtMu1vsPtMu2_,
         1, 0, 0, kBlack, .1, 20, "p_{T}(#mu_{1})", .04, .04, 1.1, "p_{T}(#mu_{2})", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist1D(PileupWeightsCanvas_, PileupWeights_,
         1, 0, 0, kBlack, .1, 20, "Pileup Weight (data/MC)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(GenWeightsCanvas_, GenWeights_,
         1, 0, 0, kBlack, .1, 20, "Gen Weight (EventWeight * LumiData / LumiMC)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(HighestCSVInclJetCanvas_, HighestCSVInclJet_,
         1, 0, 0, kBlack, .1, 20, "highest CSV Inclusive value jet in Event", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(HighestCombMVAJetCanvas_, HighestCombMVAJet_,
         1, 0, 0, kBlack, .1, 20, "highest Combined MVA value jet in Event", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(ZMassdRCanvas_, ZMassdR_,
         1, 0, 0, kBlack, .1, 20, "#DeltaR #mu_{1} #tau_{#mu}", .04, .04, 1.1,  "", .04, .04, 1.0, false);


std::cout << "after formatting" << std::endl;
  
////////////////////////
// For Pt Gen
////////////////////////
std::cout << "<----------------Formatted Canvases and Histos-------------->" << std::endl;

  //Write output file
  out_->cd();

  NMedIsoTausWithMu3Canvas_.Write();
  InvMassMu1TauMuCanvas_.Write();
  InvMassMu2TauMuCanvas_.Write();
  InvMassTauHadMu3Canvas_.Write();
  PtTauHadMu3Canvas_.Write();
  InvMassUpsilonRangeCanvas_.Write();
  InvMassZPeakRangeCanvas_.Write();
  InvMassFullRangeCanvas_.Write();
  InvMassDiMuBarrBarrCanvas_.Write();
  InvMassDiMuBarrEndcCanvas_.Write();
  InvMassDiMuEndcEndcCanvas_.Write();
  InvMassDiMuBarrOverCanvas_.Write();
  InvMassDiMuOverOverCanvas_.Write();
  InvMassDiMuEndcOverCanvas_.Write();
  Mu1PtCanvas_.Write();
  Mu2PtCanvas_.Write();
  DiMuPtCanvas_.Write();
  DiMuPtSmallerdRCanvas_.Write();
  Mu1EtaCanvas_.Write();
  Mu2EtaCanvas_.Write();
  DiTauEtaCanvas_.Write();
  DiTaudRCanvas_.Write();
  DiTauPhiCanvas_.Write();
  DiMuEtaCanvas_.Write();
  DiMudRCanvas_.Write();
  DiMuPhiCanvas_.Write();
  EtMETCanvas_.Write();
  DiTauDiMuMassCanvas_.Write();
  DiTauMassSmallerdRCanvas_.Write();
  DiMuDiTauDeltaPhiCanvas_.Write();
  METDiTauDeltaPhiCanvas_.Write();
  METDiMuDeltaPhiCanvas_.Write();
  PtMu1vsPtMu2Canvas_.Write();
  PileupWeightsCanvas_.Write();
  GenWeightsCanvas_.Write();
  HighestCSVInclJetCanvas_.Write();
  HighestCombMVAJetCanvas_.Write();
  ZMassdRCanvas_.Write();

  NEvents_->Write();
  NMedIsoTausWithMu3_->Write();
  InvMassMu1TauMu_->Write();
  InvMassMu2TauMu_->Write();
  InvMassTauHadMu3_->Write();
  PtTauHadMu3_->Write();
  InvMassUpsilonRange_->Write();
  InvMassZPeakRange_->Write();
  InvMassFullRange_->Write();
  InvMassDiMuBarrBarr_->Write();
  InvMassDiMuBarrEndc_->Write();
  InvMassDiMuEndcEndc_->Write();
  InvMassDiMuBarrOver_->Write();
  InvMassDiMuOverOver_->Write();
  InvMassDiMuEndcOver_->Write();
  Mu1Pt_->Write();
  Mu2Pt_->Write();
  DiMuPt_->Write();
  DiMuPtSmallerdR_->Write();
  Mu1Eta_->Write();
  Mu2Eta_->Write();
  DiTauEta_->Write();
  DiTaudR_->Write();
  DiTauPhi_->Write();
  DiMuEta_->Write();
  DiMudR_->Write();
  DiMuPhi_->Write();
  EtMET_->Write();
  DiTauDiMuMass_->Write();
  DiTauMassSmallerdR_->Write();
  DiMuDiTauDeltaPhi_->Write();
  METDiTauDeltaPhi_->Write();
  METDiMuDeltaPhi_->Write();
  PtMu1vsPtMu2_->Write();
  PileupWeights_->Write();
  GenWeights_->Write();
  HighestCSVInclJet_->Write();
  HighestCombMVAJet_->Write();
  ZMassdR_->Write();


  out_->Write();
  out_->Close();
std::cout << "DONE" << std::endl;
}//EndJob

// ------------ method called when starting to processes a run  ------------
void DiMu_ExtraPlots_Data::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void DiMu_ExtraPlots_Data::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DiMu_ExtraPlots_Data::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DiMu_ExtraPlots_Data::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

//Delete Memory
void DiMu_ExtraPlots_Data::reset(const bool doDelete)
{
  if ((doDelete) && (NMedIsoTausWithMu3_ != NULL)) delete NMedIsoTausWithMu3_;
  NMedIsoTausWithMu3_ = NULL;
  if ((doDelete) && (InvMassMu1TauMu_ != NULL)) delete InvMassMu1TauMu_;
  InvMassMu1TauMu_ = NULL;
  if ((doDelete) && (InvMassMu2TauMu_ != NULL)) delete InvMassMu2TauMu_;
  InvMassMu2TauMu_ = NULL;
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
  if ((doDelete) && (InvMassDiMuBarrOver_ != NULL)) delete InvMassDiMuBarrOver_;
  InvMassDiMuBarrOver_ = NULL;
  if ((doDelete) && (InvMassDiMuOverOver_ != NULL)) delete InvMassDiMuOverOver_;
  InvMassDiMuOverOver_ = NULL;
  if ((doDelete) && (InvMassDiMuEndcOver_ != NULL)) delete InvMassDiMuEndcOver_;
  InvMassDiMuEndcOver_ = NULL;
  if ((doDelete) && (Mu1Pt_ != NULL)) delete Mu1Pt_;
  Mu1Pt_ = NULL;
  if ((doDelete) && (Mu2Pt_ != NULL)) delete Mu2Pt_;
  Mu2Pt_ = NULL;
  if ((doDelete) && (DiMuPt_ != NULL)) delete DiMuPt_;
  DiMuPt_ = NULL;
  if ((doDelete) && (DiMuPtSmallerdR_ != NULL)) delete DiMuPtSmallerdR_;
  DiMuPtSmallerdR_ = NULL;
  if ((doDelete) && (Mu1Eta_ != NULL)) delete Mu1Eta_;
  Mu1Eta_ = NULL;
  if ((doDelete) && (Mu2Eta_ != NULL)) delete Mu2Eta_;
  Mu2Eta_ = NULL;
  if ((doDelete) && (DiTauEta_ != NULL)) delete DiTauEta_;
  DiTauEta_ = NULL;
  if ((doDelete) && (DiTaudR_ != NULL)) delete DiTaudR_;
  DiTaudR_ = NULL;
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
  if ((doDelete) && (DiTauMassSmallerdR_ != NULL)) delete DiTauMassSmallerdR_;
  DiTauMassSmallerdR_ = NULL;
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
  if ((doDelete) && (ZMassdR_ != NULL)) delete ZMassdR_;
  ZMassdR_ = NULL;

}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiMu_ExtraPlots_Data::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMu_ExtraPlots_Data);


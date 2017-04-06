// -*- C++ -*-
//
// Package:    SkimCheck
// Class:      SkimCheck
// 
/**\class SkimCheck SkimCheck.cc Analyzer/src/SkimCheck.cc

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

class SkimCheck : public edm::EDAnalyzer {
   public:
      typedef reco::JetFloatAssociation::Container JetBCEnergyRatioCollection;
      explicit SkimCheck(const edm::ParameterSet&);
      ~SkimCheck();

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
      edm::EDGetTokenT<edm::ValueMap<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > > muonMapTag_;
      edm::EDGetTokenT<edm::ValueMap<reco::PFJetRef> >  jetValMapTag_;
      edm::EDGetTokenT<edm::RefVector<vector<reco::PFTau>,reco::PFTau,edm::refhelper::FindUsingAdvance<vector<reco::PFTau>,reco::PFTau> > > tauTag_;
      edm::EDGetTokenT<reco::PFTauDiscriminator>  medIsoTag_;
      edm::EDGetTokenT<reco::PFTauDiscriminator>  medIsoTagMVA_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> tightIsoTag_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> veryTightIsoTag_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> decayModeFindingTag_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> isoRawTag_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> isoRawTagMVA_;
      edm::EDGetTokenT<reco::PFJetCollection>  oldJetTag_;
      edm::EDGetTokenT<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > mu3Tag_;
      edm::EDGetTokenT<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > mu12Tag_;
      edm::EDGetTokenT<reco::JetTagCollection> csvBTag_;

      //Histograms
      TH1F* TauMuTauHaddR_;
      TH1F* MassDiMuRECO_;
      TH1F* NEventsCuts_;
      TH1F* MassDiTau_;
      TH1F* PtTauMu_;
      TH1F* PtOverMassMu1Mu2_;
      TH1F* PtMu1Mu2_;
      TH1F* PCSVBDisc_;
      TH1F* SumHT_;
      TH1F* IsoRaw_;
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
SkimCheck::SkimCheck(const edm::ParameterSet& iConfig):
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  akJetTag_(consumes<vector<reco::PFJet> >(iConfig.getParameter<edm::InputTag>("akJetTag"))),
  muonsTag_(consumes<vector<reco::Muon> >(iConfig.getParameter<edm::InputTag>("muonsTag"))),
  muonMapTag_(consumes<edm::ValueMap<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > >(iConfig.getParameter<edm::InputTag>("muonMapTag"))),
  jetValMapTag_(consumes<edm::ValueMap<reco::PFJetRef> >(iConfig.getParameter<edm::InputTag>("jetValMapTag"))),
  tauTag_(consumes<edm::RefVector<vector<reco::PFTau>,reco::PFTau,edm::refhelper::FindUsingAdvance<vector<reco::PFTau>,reco::PFTau> > >(iConfig.getParameter<edm::InputTag>("tauTag"))),
  medIsoTag_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("medIsoTag"))),
  medIsoTagMVA_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("medIsoTagMVA"))),
  tightIsoTag_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("tightIsoTag"))),
  veryTightIsoTag_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("veryTightIsoTag"))),
  decayModeFindingTag_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("decayModeFindingTag"))),
  isoRawTag_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("isoRawTag"))),
  isoRawTagMVA_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("isoRawTagMVA"))),
  oldJetTag_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("oldJetTag"))),
  mu3Tag_(consumes<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > >(iConfig.getParameter<edm::InputTag>("mu3Tag"))),
  mu12Tag_(consumes<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > >(iConfig.getParameter<edm::InputTag>("mu12Tag"))),
  csvBTag_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("csvBTag")))
{
  reset(false);    
}//SkimCheck



SkimCheck::~SkimCheck()
{
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void SkimCheck::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "\n<------------THIS IS A NEW EVENT------------>" << std::endl;

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

  //Get CleanJets Tau particle collection
  edm::Handle<edm::RefVector<vector<reco::PFTau>,reco::PFTau,edm::refhelper::FindUsingAdvance<vector<reco::PFTau>,reco::PFTau> > > pTaus;
  iEvent.getByToken(tauTag_, pTaus);

  //Get Medium Iso Collection
  Handle<PFTauDiscriminator> pMedIsoDisc;
  iEvent.getByToken(medIsoTag_, pMedIsoDisc);

  //Get Medium Iso Collection
  Handle<PFTauDiscriminator> pMedIsoDiscMVA;
  iEvent.getByToken(medIsoTagMVA_, pMedIsoDiscMVA);

  //Get Tight Iso Collection
  Handle<PFTauDiscriminator> pTightIsoDisc; 
  iEvent.getByToken(tightIsoTag_, pTightIsoDisc);

  //Get Very Tight Iso Collection 
  Handle<PFTauDiscriminator> pVTightIsoDisc;
  iEvent.getByToken(veryTightIsoTag_, pVTightIsoDisc);

  //Get Decay Mode Finding Collection
  Handle<PFTauDiscriminator> pDMFinding; 
  iEvent.getByToken(decayModeFindingTag_, pDMFinding);

  //Get IsoRaw  Collection
  Handle<PFTauDiscriminator> pIsoRaw;
  iEvent.getByToken(isoRawTag_, pIsoRaw);

  //Get IsoRaw  Collection
  Handle<PFTauDiscriminator> pIsoRawMVA;
  iEvent.getByToken(isoRawTagMVA_, pIsoRawMVA);

  //Old Jet collection for bTagging
  edm::Handle<reco::PFJetCollection> pOldJets;
  iEvent.getByToken(oldJetTag_, pOldJets);

  //Old Jet collection for bTagging
  edm::Handle<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > pMu12;
  iEvent.getByToken(mu12Tag_, pMu12);

  //Get combVertMVA JetTagCollection
  edm::Handle<reco::JetTagCollection> pCSV;
  iEvent.getByToken(csvBTag_, pCSV);

//////////////////////////////
// Begin Analyzer
//////////////////////////////
  //Getting Mu1 and Mu2
  reco::MuonRef mu1Ref = reco::MuonRef((*pMu12)[0] );
  reco::MuonRef mu2Ref = reco::MuonRef((*pMu12)[1] );
  reco::LeafCandidate::LorentzVector diMuP4;
  diMuP4 = mu1Ref->p4();
  diMuP4 += mu2Ref->p4();
  PtOverMassMu1Mu2_->Fill(diMuP4.Pt() / diMuP4.M() );
  MassDiMuRECO_->Fill(diMuP4.M() );
  PtMu1Mu2_->Fill(diMuP4.Pt() );  //mu1Ref->pt() + mu2Ref->pt() );

  // Getting Tau_Had
  size_t numTauHad = pTaus->size(), tauHadIndex = -1;
  double largestPt = -1;
  for ( size_t iTau=0; iTau < numTauHad; ++iTau)
  {
    if (reco::PFTauRef((*pTaus)[iTau] )->pt() > largestPt);
    {
      tauHadIndex = iTau;
      largestPt = reco::PFTauRef((*pTaus)[iTau] )->pt();
    }//if largets tauHad pt
  }//for
  reco::PFTauRef tauHadRef = reco::PFTauRef((*pTaus)[tauHadIndex] );
  const reco::PFJetRef& tauJetRef = (*tauHadRef).jetRef();
  unsigned int TightIso = (*pTightIsoDisc)[tauHadRef], VTightIso = (*pVTightIsoDisc)[tauHadRef], DM = (*pDMFinding)[tauHadRef];
  unsigned int MedIso = (*pMedIsoDisc)[tauHadRef], MedIsoMVA = (*pMedIsoDiscMVA)[tauHadRef];
  double rawIsoValue = (*pIsoRaw)[tauHadRef], rawIsoValueMVA = (*pIsoRawMVA)[tauHadRef];;
  IsoRaw_->Fill(rawIsoValueMVA );
  RelIsoRaw_->Fill(rawIsoValueMVA / tauHadRef->pt() );
  
  if (MedIso == 1)
    IsoRawCutBased_->Fill(rawIsoValue );
  std::cout << "numTauHad=" << numTauHad << "\n\tRawIso=" << rawIsoValue << "\tMedIso=" << MedIso << "\n\tRawIsoMVA=" << rawIsoValueMVA << "\tMedIsoMVA=" << MedIsoMVA << 
	       "\n\tTightIso=" << TightIso << "\tVTightIso=" << VTightIso << std::endl;

  //find the highest pT associated muon
  const reco::MuonRefVector& removedMuons = (*pMuonMap)[tauJetRef];
  std::vector<reco::MuonRef> removedMuonRefs;
  for (reco::MuonRefVector::const_iterator iMuon = removedMuons.begin(); iMuon != removedMuons.end(); ++iMuon) {removedMuonRefs.push_back(*iMuon);}//for iMuon
  std::cout << "\tremovedMuonRefs.size= " << removedMuonRefs.size() << std::endl;
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
  
  
  double dPhi = reco::deltaPhi(removedMuonRefs[0]->phi(), tauHadRef->phi() );
  double dR_diTau = sqrt( (removedMuonRefs[0]->eta() - tauHadRef->eta() ) * (removedMuonRefs[0]->eta() - tauHadRef->eta() )  +  dPhi * dPhi );
  TauMuTauHaddR_->Fill(dR_diTau);
  reco::LeafCandidate::LorentzVector diTauP4;
  diTauP4 = removedMuonRefs[0]->p4() + tauHadRef->p4();
  MassDiTau_->Fill(diTauP4.M() );
  std::vector<reco::PFCandidatePtr> JetPFCands = tauJetRef->getPFConstituents();
  PtTauMu_->Fill(removedMuonRefs[0]->pt() );
  std::cout << "removedMuonRefs[0]->pt()=" << removedMuonRefs[0]->pt() << "\tRawIso-MuPt=" << rawIsoValue - removedMuonRefs[0]->pt() << std::endl;
 
  //sumHT
  reco::LeafCandidate::LorentzVector sumHTP4;
  size_t numJets = pAkJets->size();
  bool firstHTJet = true;
  for ( size_t iJet = 0; iJet < numJets; ++iJet ) 
  {
    reco::PFJetRef jetRef(pAkJets, iJet);
    if (jetRef->pt() > 20.0 && !firstHTJet)
      sumHTP4 += jetRef->p4();
    if (jetRef->pt() > 20.0 && firstHTJet)
    {
      sumHTP4 = jetRef->p4();
      firstHTJet = false;
    }//if
  }//for pAkJets      
  if (tauHadRef->pt() > 20.0 )
    sumHTP4 -= tauHadRef->p4();
  if (mu1Ref->pt() > 20 )
    sumHTP4 -= mu1Ref->p4();
  if (mu2Ref->pt() > 20 )
    sumHTP4 -= mu2Ref->p4();
  SumHT_->Fill(sumHTP4.Pt() );
  

  double BDisc = -1;
  for (unsigned int iBTagInfo = 0; iBTagInfo != pCSV->size(); ++iBTagInfo)
  {
    if (reco::PFJetRef(pOldJets, iBTagInfo).key() == tauJetRef.key() )
    {
      const reco::JetTagCollection& bTags = *(pCSV.product());
      BDisc = bTags[iBTagInfo].second;
    }//if
  }//for iBTagInfo
  PCSVBDisc_->Fill(BDisc);


  NEventsCuts_->Fill(0);
  if (BDisc < .5) 
    NEventsCuts_->Fill(1);
  if (BDisc < .75)
    NEventsCuts_->Fill(2);
  if (BDisc < .9 )
    NEventsCuts_->Fill(3);
  if (BDisc < .9 && TightIso == 1 && DM == 1)
    NEventsCuts_->Fill(4);
  if (sumHTP4.Pt() < 300)
    NEventsCuts_->Fill(5);
  if (sumHTP4.Pt() < 400)
    NEventsCuts_->Fill(6);
  if (sumHTP4.Pt() < 300 && BDisc < .75)
    NEventsCuts_->Fill(7);
  if (sumHTP4.Pt() < 400 && BDisc < .5 )
    NEventsCuts_->Fill(8);
  if (rawIsoValue < .5)
    NEventsCuts_->Fill(9);
  if (rawIsoValue < 1)
    NEventsCuts_->Fill(10);
  if (TightIso == 1 && DM == 1)
    NEventsCuts_->Fill(11);
  if (VTightIso == 1 && DM == 1)
    NEventsCuts_->Fill(12);
 
}//End SkimCheck::analyze


// ------------ method called once each job just before starting event loop  ------------
void SkimCheck::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  //Open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //Book histograms
  TauMuTauHaddR_       = new TH1F("TauMuTauHaddR"    , "", 50, 0, 1.5);
  MassDiMuRECO_       = new TH1F("MassDiMuRECO"    , "", 60, 0, 120);
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
  Float_t bins[] = {0, 1, 2, 3, 4, 20};
  MassDiTau_        = new TH1F("MassDiTau"    , "", sizeof(bins)/sizeof(Float_t) - 1, bins);
  PtTauMu_        = new TH1F("PtTauMu"    , "", 50, 0, 100);
  PtOverMassMu1Mu2_        = new TH1F("PtOverMassMu1Mu2"    , "", 50, 0, 50);
  PtMu1Mu2_        = new TH1F("PtMu1Mu2"    , "", 50, 0, 700);
  PCSVBDisc_        = new TH1F("PCSVBDisc"    , "", 50, 0, 1);
  SumHT_        = new TH1F("SumHT"    , "", 50, 0, 800);
  IsoRaw_        = new TH1F("IsoRaw"    , "", 100, 0, 10);
  RelIsoRaw_        = new TH1F("RelIsoRaw"    , "", 100, 0, 5);
  IsoRawCutBased_        = new TH1F("IsoRawCutBased"    , "", 100, 0, 10);

}

// ------------ method called once each job just after ending the event loop  ------------
void SkimCheck::endJob()
{
  //Make the Canvases
  TCanvas TauMuTauHaddRCanvas("TauMuTauHaddR","",600,600);
  TCanvas MassDiMuRECOCanvas("MassDiMuRECO","",600,600);
  TCanvas NEventsCutsCanvas("NEventsCuts","",600,600);
  TCanvas MassDiTauCanvas("MassDiTau","",600,600);
  TCanvas PtTauMuCanvas("PtTauMu","",600,600);
  TCanvas PtOverMassMu1Mu2Canvas("PtOverMassMu1Mu2","",600,600);
  TCanvas PtMu1Mu2Canvas("PtMu1Mu2","",600,600);
  TCanvas PCSVBDiscCanvas("PCSVBDisc","",600,600);
  TCanvas SumHTCanvas("SumHT","",600,600);
  TCanvas IsoRawCanvas("IsoRaw","",600,600);
  TCanvas RelIsoRawCanvas("RelIsoRaw","",600,600);
  TCanvas IsoRawCutBasedCanvas("IsoRawCutBased","",600,600);

std::cout << "<----------------Declared Canvases-------------->" << std::endl;

  //Format the 1D plots and Draw (canvas, hist, grid, log y, log z, color, size, style, xAxisTitle, xTitleSize, xLabelSize, xTitleOffSet, yAxisTitle, yTitleSize, yLabelSize, yTitleOffset)
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauMuTauHaddRCanvas, TauMuTauHaddR_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(#tau_{mu},#tau_{H})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MassDiMuRECOCanvas, MassDiMuRECO_,
         1, 0, 0, kBlack, 7, 20, "RECO m_{#mu #mu}", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NEventsCutsCanvas, NEventsCuts_,
         1, 0, 0, kBlack, 7, 20, "", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MassDiTauCanvas, MassDiTau_,
         1, 0, 0, kBlack, 7, 20, "RECO Visible m_{#tau_{mu} #tau_{H} }", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(PtTauMuCanvas, PtTauMu_,
         1, 0, 0, kBlack, 7, 20, "p_{T} (#tau_{#mu})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(PtOverMassMu1Mu2Canvas, PtOverMassMu1Mu2_,
         1, 0, 0, kBlack, 7, 20, "p_{T}(#mu_{1}#mu_{2}) / Mass(#mu_{1}#mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(PtMu1Mu2Canvas, PtMu1Mu2_,
         1, 0, 0, kBlack, 7, 20, "p_{T}(#mu_{1}#mu_{2})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(PCSVBDiscCanvas, PCSVBDisc_,
         1, 0, 0, kBlack, 7, 20, "pCSV BDiscriminant Value", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(SumHTCanvas, SumHT_,
         1, 0, 0, kBlack, 7, 20, "HT(sum p_{T}(ak4 Jet > 20) - p_{T}(#tau_{had},#mu_{1},#mu_{2} if > 20)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoRawCanvas, IsoRaw_,
         1, 0, 0, kBlack, 7, 20, "Raw Isolation Value", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(RelIsoRawCanvas, RelIsoRaw_,
         1, 0, 0, kBlack, 7, 20, "Raw Isolation Value / p_{T}(#tau_{Had})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoRawCutBasedCanvas, IsoRawCutBased_,
         1, 0, 0, kBlack, 7, 20, "Raw Isolation Value passing Med CutBased", .04, .04, 1.1,  "", .04, .04, 1.0, false);

std::cout << "after formatting" << std::endl;
  
std::cout << "<----------------Formatted Canvases and Histos-------------->" << std::endl;

  //Write output file
  out_->cd();

  TauMuTauHaddRCanvas.Write();
  MassDiMuRECOCanvas.Write();
  NEventsCutsCanvas.Write();
  MassDiTauCanvas.Write();
  PtTauMuCanvas.Write();
  PtOverMassMu1Mu2Canvas.Write();
  PtMu1Mu2Canvas.Write();
  PCSVBDiscCanvas.Write();
  SumHTCanvas.Write();
  IsoRawCanvas.Write();
  RelIsoRawCanvas.Write();
  IsoRawCutBasedCanvas.Write();

  out_->Write();
  out_->Close();
std::cout << "DONE" << std::endl;
}//EndJob

// ------------ method called when starting to processes a run  ------------
void SkimCheck::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void SkimCheck::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void SkimCheck::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void SkimCheck::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

//Delete Memory
void SkimCheck::reset(const bool doDelete)
{
  if ((doDelete) && (TauMuTauHaddR_ != NULL)) delete TauMuTauHaddR_;
  TauMuTauHaddR_ = NULL;
  if ((doDelete) && (MassDiMuRECO_ != NULL)) delete MassDiMuRECO_;
  MassDiMuRECO_ = NULL;
  if ((doDelete) && (NEventsCuts_ != NULL)) delete NEventsCuts_;
  NEventsCuts_ = NULL;
  if ((doDelete) && (MassDiTau_ != NULL)) delete MassDiTau_;
  MassDiTau_ = NULL;
  if ((doDelete) && (PtTauMu_ != NULL)) delete PtTauMu_;
  PtTauMu_ = NULL;
  if ((doDelete) && (PtOverMassMu1Mu2_ != NULL)) delete PtOverMassMu1Mu2_;
  PtOverMassMu1Mu2_ = NULL;
  if ((doDelete) && (PtMu1Mu2_ != NULL)) delete PtMu1Mu2_;
  PtMu1Mu2_ = NULL;
  if ((doDelete) && (PCSVBDisc_ != NULL)) delete PCSVBDisc_;
  PCSVBDisc_ = NULL;
  if ((doDelete) && (SumHT_ != NULL)) delete SumHT_;
  SumHT_ = NULL;
  if ((doDelete) && (IsoRaw_ != NULL)) delete IsoRaw_;
  IsoRaw_ = NULL;
  if ((doDelete) && (RelIsoRaw_ != NULL)) delete RelIsoRaw_;
  RelIsoRaw_ = NULL;
  if ((doDelete) && (IsoRawCutBased_ != NULL)) delete IsoRawCutBased_;
  IsoRawCutBased_ = NULL;

}//void SkimCheck

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SkimCheck::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SkimCheck);

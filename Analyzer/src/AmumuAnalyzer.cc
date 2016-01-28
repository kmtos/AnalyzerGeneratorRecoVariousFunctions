// -*- C++ -*-
//
// Package:    BBAAnalyzer
// Class:      BBAAnalyzer
// 
/**\class BBAAnalyzer BBAAnalyzer.cc Analyzer/src/BBAAnalyzer.cc

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
#include <map>
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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "HLTrigger/HLTcore/interface/HLTConfigData.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTanalyzers/interface/HLTInfo.h"

#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"

#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

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

class AmumuAnalyzer : public edm::EDAnalyzer {
   public:
      typedef reco::JetFloatAssociation::Container JetBCEnergyRatioCollection;
      explicit AmumuAnalyzer(const edm::ParameterSet&);
      ~AmumuAnalyzer();

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

      //gen particle Tag
      edm::InputTag prunedGenParticleTag_;

      //PAT jet Tag
      edm::InputTag slimmedJetTag_;

      //PAT TAU Tag
      edm::InputTag slimmedTauTag_;

      //PAT Muon Tag
      edm::EDGetTokenT< pat::MuonCollection > slimmedMuonTag_;//      edm::InputTag slimmedMuonTag_;

      // DiMuTrig  Muon Tag
      std::string  mOniaTrigMatchTag_;
      
      // input for patTrigger
      edm::InputTag trigger_;

      // input for patTriggerEvent
      edm::EDGetTokenT< pat::TriggerEvent > triggerEventToken_;//      edm::InputTag triggerEventToken_;

      //Histograms
      TH1F* HighMuPt_;
      TH1F* LowMuPt_;

      TH1F* NTrigMatchMu_;
      TH2F* SlimdRtvstrigMatchdR_;
      TH2F* IsoTrigPtvsRecoPt_;
      TH2F* NonIsoTrigEtavsRecoEta_;
      TH1F* NonIsoMatchDiTauDr_;

      TH1F* NonIsoBCSVTurnOn_;
      TH1F* IsoBCSVTurnOn_;
      TH1F* NonIsoBCSVPassAll_;
      TH1F* IsoBCSVPassAll_;
      TH1F* NonIsoBCSVPassTrig_;
      TH1F* IsoBCSVPassTrig_;

      TH1F* TrigMatchdR_;
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
AmumuAnalyzer::AmumuAnalyzer(const edm::ParameterSet& iConfig):
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  prunedGenParticleTag_(iConfig.getParameter<edm::InputTag>("prunedGenParticleTag")),
  slimmedJetTag_(iConfig.getParameter<edm::InputTag>("slimmedJetTag")),
  slimmedTauTag_(iConfig.getParameter<edm::InputTag>("slimmedTauTag")),
  slimmedMuonTag_( consumes< pat::MuonCollection >( iConfig.getParameter< edm::InputTag >( "slimmedMuonTag" ) ) ),
  mOniaTrigMatchTag_( iConfig.getParameter< std::string >( "mOniaTrigMatchTag" ) ),
  trigger_( iConfig.getParameter< edm::InputTag >( "trigger" ) ),
  triggerEventToken_( consumes< pat::TriggerEvent >( iConfig.getParameter< edm::InputTag >( "triggerEvent" ) ) )
{
  reset(false);    
}//AmumuAnalyzer



AmumuAnalyzer::~AmumuAnalyzer()
{
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void AmumuAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "\n<------------THIS IS A NEW EVENT------------>" << std::endl;
///////////////////////////////////////////////////////
//  Declare Event Handles
///////////////////////////////////////////////////////

  // slimmedJet handle
  edm::Handle<std::vector<pat::Jet> > pSlimmedJetTag;
  iEvent.getByLabel(slimmedJetTag_, pSlimmedJetTag);

  // slimmedMu handle
  edm::Handle<std::vector<pat::Tau> > pSlimmedTauTag;
  iEvent.getByLabel(slimmedTauTag_, pSlimmedTauTag);

  // PAT trigger event
  edm::Handle< pat::TriggerEvent > triggerEvent;
  iEvent.getByToken( triggerEventToken_, triggerEvent );

  // PAT object collection
  edm::Handle< pat::MuonCollection > pSlimmedMuonTag;
  iEvent.getByToken( slimmedMuonTag_, pSlimmedMuonTag );

  //Get gen particle collection
  edm::Handle<reco::GenParticleCollection> pPrunGenParts;
  iEvent.getByLabel(prunedGenParticleTag_, pPrunGenParts);

///////////////////////////////////////////////////////
// Get Gen Mus
///////////////////////////////////////////////////////
  reco::GenParticleCollection::const_iterator firstGen = pPrunGenParts->begin();
  reco::GenParticleRef mu1Ref = firstGen->daughterRef(0), mu2Ref = firstGen->daughterRef(0);
  for (reco::GenParticleCollection::const_iterator iGenParticle = pPrunGenParts->begin(); iGenParticle != pPrunGenParts->end(); ++iGenParticle)
  {
    if(iGenParticle->pdgId() == 36 && iGenParticle->numberOfDaughters() == 2)
    {
      mu1Ref = iGenParticle->daughterRef(0);
      mu2Ref = iGenParticle->daughterRef(1);
      break;
    }//if pdgid == 36
  }//for iGen Particle

///////////////////////////////////////////////////////
// Slim Mus to Gen Mus 
///////////////////////////////////////////////////////
std::cout << "\t<---------Matching Mu's---------->" << std::endl;
  double dR1 = 1000, dR2 = 1000, dR1Temp = 0, dR2Temp = 0;
  double slimMu1Pt = -1, slimMu2Pt = -1, slimMu1Eta = -10, slimMu2Eta = -10, slimMu1Phi = -10, slimMu2Phi = -10;
  for (auto iMu = pSlimmedMuonTag->begin(); iMu != pSlimmedMuonTag->end(); ++iMu)
  {
    dR1Temp = sqrt( (mu1Ref->eta() - iMu->eta() ) * (mu1Ref->eta() - iMu->eta() ) + (mu1Ref->phi() - iMu->phi() ) * (mu1Ref->phi() - iMu->phi() ) );
    dR2Temp = sqrt( (mu2Ref->eta() - iMu->eta() ) * (mu2Ref->eta() - iMu->eta() ) + (mu2Ref->phi() - iMu->phi() ) * (mu2Ref->phi() - iMu->phi() ) );

    std::cout << "\t\tmu1Ref: Eta= " << mu1Ref->eta() << "  Phi= " << mu1Ref->phi() << "  mu2Ref: Eta= " << mu2Ref->eta() << "  Phi= " << mu2Ref->phi() << std::endl;
    std::cout << "\t\tiMu:    Eta= " << iMu->eta() << "  Phi= " << iMu->phi() << std::endl;
    std::cout << "\t\tdR1= " << dR1 << "  dR1Temp= " << dR1Temp << "  PtDiff1= " << abs(slimMu1Pt - iMu->pt() ) << std::endl;
    std::cout << "\t\tdR2= " << dR2 << "  dR2Temp= " << dR2Temp << "  PtDiff2= " << abs(slimMu2Pt - iMu->pt() ) << std::endl;

    //first 2 checks are basic requirements, the 2 in the "or" are making sure it doesn't matche better with the second tau
    if (dR1Temp < .4 && dR1Temp < dR1 && (dR1Temp < dR2Temp || dR2Temp > dR2) )
    {
      dR1 = dR1Temp;
      slimMu1Pt = iMu->pt();
      slimMu1Eta = iMu->eta();
      slimMu1Phi = iMu->phi();
    }//if iMu matches mu1Ref
    else if (dR2Temp < .4 && dR2Temp < dR2 && (dR2Temp < dR1Temp || dR1Temp > dR1) )
    {
      dR2 = dR2Temp;
      slimMu2Pt = iMu->pt();
      slimMu2Eta = iMu->eta();
      slimMu2Phi = iMu->phi();
    }//else
  }//for auto iMu


  std::cout << "\tSLIM Match to Gen: dR1= " << dR1 << "  dR1Temp= " << dR1Temp << "  Pt= " << slimMu1Pt << "  Eta= " << slimMu1Eta << "  Phi= " << slimMu1Phi << std::endl;
  std::cout << "\tSLIM Match to Gen: dR2= " << dR2 << "  dR2Temp= " << dR2Temp << "  Pt= " << slimMu2Pt << "  Eta= " << slimMu2Eta << "  Phi= " << slimMu2Phi << std::endl;
  double diMudR_slim = sqrt( (slimMu2Eta - slimMu1Eta)*(slimMu2Eta - slimMu1Eta)  +  (slimMu2Phi - slimMu1Phi)*(slimMu2Phi - slimMu1Phi) );
///////////////////////////////////////////////////////////////////////////////////////////////
// Loop for Non- Muon Trigger Matching 
///////////////////////////////////////////////////////////////////////////////////////////////
  std::cout << "\t!!!!!!!!!!!!! TRIGGER MATCHING !!!!!!!!!!!!!!" << std::endl;
  const pat::helper::TriggerMatchHelper matchHelper;
  double dRMatch1 = 1000, dRMatch2 = 1000, dRTEMP = 1000, muMatchPt1 = -100, muMatchPt2 = -100, muMatchPhi1 = -100, muMatchPhi2 = -100, muMatchEta1 = -100, muMatchEta2 = -100;
  int nMuMatch = 0, muMatched = 0;
  std::cout << "pSlimmedMuonTag->size()= " << pSlimmedMuonTag->size() << std::endl;
  for( size_t iMuon = 0; iMuon < pSlimmedMuonTag->size(); ++iMuon )
  {
    const pat::TriggerObjectRef trigRef( matchHelper.triggerMatchObject( pSlimmedMuonTag, iMuon, mOniaTrigMatchTag_, iEvent, *triggerEvent ) );
    if ( trigRef.isAvailable() && trigRef.isNonnull() )
    {
      cout << "<---------A Muon Trigger Match----------->" << endl; 
      dRTEMP = sqrt( (trigRef->eta() - pSlimmedMuonTag->at(iMuon).eta())*(trigRef->eta() - pSlimmedMuonTag->at(iMuon).eta()) 
		+  (trigRef->phi() - pSlimmedMuonTag->at(iMuon).phi())*(trigRef->phi() - pSlimmedMuonTag->at(iMuon).phi()) );
      std::cout << "\tSLIM Match to Gen: dTEMP1= " << dRTEMP <<  "  Pt_diff= " << fabs(trigRef->pt() - pSlimmedMuonTag->at(iMuon).pt() ) << std::endl;
      if (muMatched == 0)  //Store the matches in first muMatch elements
      {
        nMuMatch++;
        muMatchPt1 = trigRef->pt();
        muMatchEta1 = trigRef->eta();
        muMatchPhi1 = trigRef->phi();
	dRMatch1 = dRTEMP;
      }//if muMatched = 0
      else if (muMatched == 1)  //if already found 1 match, then store elements in 2nd mumatch
      {
        nMuMatch++;
        muMatchPt2 = trigRef->pt();
        muMatchEta2 = trigRef->eta();
        muMatchPhi2 = trigRef->phi();
        dRMatch2 = dRTEMP;
      }//else if muMatched = 0
      else //if > 2 muons pass the trigger, then it sees which 
      {
        nMuMatch++;
	double dR_mu13 = sqrt( (muMatchEta1 - trigRef->eta() )*(muMatchEta1 - trigRef->eta() )  +  (muMatchPhi1 - trigRef->phi() )*(muMatchPhi1 - trigRef->phi() )  );
        double dR_mu23 = sqrt( (muMatchEta2 - trigRef->eta() )*(muMatchEta2 - trigRef->eta() )  +  (muMatchPhi2 - trigRef->phi() )*(muMatchPhi2 - trigRef->phi() )  );
        double dR_mu12 = sqrt( (muMatchEta2 - muMatchEta1)*(muMatchEta2 - muMatchEta1)  +  (muMatchPhi2 - muMatchPhi1)*(muMatchPhi2 - muMatchPhi1) );

	if (dR_mu13 < dR_mu23 && dR_mu13 < dR_mu12) //if muon trig match 1 and 3 are the closest to each other
	{
          muMatchPt2 = trigRef->pt();
          muMatchEta2 = trigRef->eta();
          muMatchPhi2 = trigRef->phi();
          dRMatch2 = dRTEMP;
	}//if dR_mu13
        else if (dR_mu23 < dR_mu13 && dR_mu23 < dR_mu12) //if muon trig match 2 and 3 are the closest to each other
        {
          muMatchPt1 = trigRef->pt();
          muMatchEta1 = trigRef->eta();
          muMatchPhi1 = trigRef->phi();
          dRMatch1 = dRTEMP;
        }//else if 
      }//else	
    }//if
  }//for iMuon
  NTrigMatchMu_->Fill(nMuMatch);
  if (dRMatch1 != 1000 && dRMatch2 != 1000)
  {
    double diMudR_trig = sqrt( (muMatchEta1 - muMatchEta2)*(muMatchEta1 - muMatchEta2)  +  (muMatchPhi1 - muMatchPhi2)*(muMatchPhi1 - muMatchPhi2) );
    SlimdRtvstrigMatchdR_->Fill(diMudR_slim, diMudR_trig);
    TrigMatchdR_->Fill(diMudR_trig );
    if (muMatchPt1 > muMatchPt2)
    {
      HighMuPt_->Fill(muMatchPt1);
      LowMuPt_->Fill(muMatchPt2);
    }//fi muMatchPt1 > muMatchPt2
    else
    {  
      HighMuPt_->Fill(muMatchPt2);
      LowMuPt_->Fill(muMatchPt1);
    }//else
  }//if two Mu Trig matches
}//End AmumuAnalyzer::analyze


// ------------ method called once each job just before starting event loop  ------------
void AmumuAnalyzer::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  //Open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //Book histograms
  HighMuPt_ 	= new TH1F("HighMuPt", "", 50, 0, 100);
  LowMuPt_ 	= new TH1F("LowMuPt",    "", 50, 0, 100);

  NTrigMatchMu_      = new TH1F("NTrigMatchMu", "", 12, -.5, 11.5);
  SlimdRtvstrigMatchdR_  = new TH2F("SlimdRtvstrigMatchdR", "", 30, 0, 2, 30, 0, 2);
  IsoTrigPtvsRecoPt_  = new TH2F("IsoTrigPtvsRecoPt", "", 30, 0, 150, 30, 0, 150);
  NonIsoTrigEtavsRecoEta_  = new TH2F("NonIsoTrigEtavsRecoEta", "", 30, 0, 3, 30, 0, 3);
  NonIsoMatchDiTauDr_ = new TH1F("NonIsoMatchDiTauDr", "", 50, 0, 8);

  NonIsoBCSVTurnOn_ = new TH1F("NonIsoBCSVTurnOn", "", 33, 0, 1);
  IsoBCSVTurnOn_    = new TH1F("IsoBCSVTurnOn",    "", 33, 0, 1);
  NonIsoBCSVPassAll_ = new TH1F("NonIsoBCSVPassAll", "", 33, 0, 1);
  IsoBCSVPassAll_    = new TH1F("IsoBCSVPassAll",    "", 33, 0, 1);
  NonIsoBCSVPassTrig_ = new TH1F("NonIsoBCSVPassTrig", "", 33, 0, 1);
  IsoBCSVPassTrig_    = new TH1F("IsoBCSVPassTrig",    "", 33, 0, 1);

  TrigMatchdR_     = new TH1F("TrigMatchdR",     "", 33, 0, 2);

}

// ------------ method called once each job just after ending the event loop  ------------
void AmumuAnalyzer::endJob()
{



  //Make the Canvases
  TCanvas HighMuPtCanvas_("HighMuPtCanvas","",600,600);
  TCanvas LowMuPtCanvas_("LowMuPtCanvas","",600,600);

  TCanvas NTrigMatchMuCanvas_("NTrigMatchMuCanvas","",600,600);
  TCanvas SlimdRtvstrigMatchdRCanvas("SlimdRtvstrigMatchdR","",600,600);
  TCanvas IsoTrigPtvsRecoPtCanvas("IsoTrigPtvsRecoPt","",600,600);
  TCanvas NonIsoTrigEtavsRecoEtaCanvas("NonIsoTrigEtavsRecoEta","",600,600);

  TCanvas NonIsoMatchDiTauDrCanvas_("NonIsoMatchDiTauDrCanvas","",600,600);
  TCanvas NonIsoBCSVTurnOnCanvas_("NonIsoBCSVTurnOnCanvas","",600,600); 
  TCanvas IsoBCSVTurnOnCanvas_("IsoBCSVTurnOnCanvas","",600,600);    
  TCanvas NonIsoBCSVPassAllCanvas_("NonIsoBCSVPassAllCanvas","",600,600);
  TCanvas IsoBCSVPassAllCanvas_("IsoBCSVPassAllCanvas","",600,600);  
  TCanvas NonIsoBCSVPassTrigCanvas_("NonIsoBCSVPassTrigCanvas","",600,600);
  TCanvas IsoBCSVPassTrigCanvas_("IsoBCSVPassTrigCanvas","",600,600);  

  TCanvas TrigMatchdRCanvas_("TrigMatchdRCanvas","",600,600);

std::cout << "<----------------Declared Canvases-------------->" << std::endl;

  //Format the 1D plots and Draw (canvas, hist, grid, log y, log z, color, size, style, xAxisTitle, xTitleSize, xLabelSize, xTitleOffSet, yAxisTitle, yTitleSize, yLabelSize, yTitleOffset)
  VariousFunctions::formatAndDrawCanvasAndHist1D(HighMuPtCanvas_, HighMuPt_, 1, 0, 0, kBlack, 7, 20, "pt(#mu high)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(LowMuPtCanvas_, LowMuPt_, 1, 0, 0, kBlack, 7, 20, "pt(#mu low)", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(NTrigMatchMuCanvas_, NTrigMatchMu_, 1, 0, 0, kBlack, 7, 20, "Number of trig Mu Matches", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist2D(SlimdRtvstrigMatchdRCanvas, SlimdRtvstrigMatchdR_, 0, 0, 0, kBlack, 7, 20, "Slim #DeltaR(#mu#mu)", .04, .04, 1.1, "Trig Slim #DeltaR(#mu#mu)", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(IsoTrigPtvsRecoPtCanvas, IsoTrigPtvsRecoPt_, 0, 0, 0, kBlack, 7, 20, "Iso Trig Pt", .04, .04, 1.1, "RECO Match Trig Pt", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(NonIsoTrigEtavsRecoEtaCanvas, NonIsoTrigEtavsRecoEta_, 0, 0, 0, kBlack, 7, 20, "NonIso Trig Eta", .04, .04, 1.1, "RECO Trig Match Eta", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NonIsoMatchDiTauDrCanvas_, NonIsoMatchDiTauDr_, 1, 0, 0, kBlack, 7, 20, "Non-Iso Trig #DeltaR(#tau#tau)", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(NonIsoBCSVTurnOnCanvas_, NonIsoBCSVTurnOn_, 1, 0, 0, kBlack, 7, 20, "Non-Iso Trig CSV(bJet) TurnOn pt > 50", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoBCSVTurnOnCanvas_, IsoBCSVTurnOn_, 1, 0, 0, kBlack, 7, 20, "Iso Trig CSV(bJet) TurnOn pt > 50", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NonIsoBCSVPassAllCanvas_, NonIsoBCSVPassAll_, 1, 0, 0, kBlack, 7, 20, "Non-Iso Trig CSV(bJet) PassAll pt > 50", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoBCSVPassAllCanvas_, IsoBCSVPassAll_, 1, 0, 0, kBlack, 7, 20, "Iso Trig CSV(bJet) PassAll pt > 50", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NonIsoBCSVPassTrigCanvas_, NonIsoBCSVPassTrig_, 1, 0, 0, kBlack, 7, 20, "Non-Iso Trig CSV(bJet) PassTrig pt > 50", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(IsoBCSVPassTrigCanvas_, IsoBCSVPassTrig_, 1, 0, 0, kBlack, 7, 20, "Iso Trig CSV(bJet) PassTrig pt > 50", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(TrigMatchdRCanvas_, TrigMatchdR_, 1, 0, 0, kBlack, 7, 20, "#DeltaR(trig_b, matched_b)", .04, .04, 1.1,  "", .04, .04, 1.0, false);

std::cout << "<----------------Formatted Canvases and Histos-------------->" << std::endl;

  //Write output file
  out_->cd();

  HighMuPtCanvas_.Write();
  LowMuPtCanvas_.Write();

  NTrigMatchMuCanvas_.Write();
  SlimdRtvstrigMatchdRCanvas.Write();
  IsoTrigPtvsRecoPtCanvas.Write();
  NonIsoTrigEtavsRecoEtaCanvas.Write();
  NonIsoMatchDiTauDrCanvas_.Write();
  
  NonIsoBCSVTurnOnCanvas_.Write();
  IsoBCSVTurnOnCanvas_.Write();
  NonIsoBCSVPassAllCanvas_.Write();
  IsoBCSVPassAllCanvas_.Write();
  NonIsoBCSVPassTrigCanvas_.Write();
  IsoBCSVPassTrigCanvas_.Write();
  
  TrigMatchdRCanvas_.Write();

  out_->Write();
  out_->Close();
std::cout << "DONE" << std::endl;
}//EndJob

// ------------ method called when starting to processes a run  ------------
void AmumuAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void AmumuAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void AmumuAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void AmumuAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

//Delete Memory
void AmumuAnalyzer::reset(const bool doDelete)
{
  if ((doDelete) && (HighMuPt_ != NULL)) delete HighMuPt_;
  HighMuPt_ = NULL;
  if ((doDelete) && (LowMuPt_ != NULL)) delete LowMuPt_;
  LowMuPt_ = NULL;

  if ((doDelete) && (NTrigMatchMu_ != NULL)) delete NTrigMatchMu_;
  NTrigMatchMu_ = NULL;
  if ((doDelete) && (SlimdRtvstrigMatchdR_ != NULL)) delete SlimdRtvstrigMatchdR_;
  SlimdRtvstrigMatchdR_ = NULL;
  if ((doDelete) && (IsoTrigPtvsRecoPt_ != NULL)) delete IsoTrigPtvsRecoPt_;
  IsoTrigPtvsRecoPt_ = NULL;
  if ((doDelete) && (NonIsoTrigEtavsRecoEta_ != NULL)) delete NonIsoTrigEtavsRecoEta_;
  NonIsoTrigEtavsRecoEta_ = NULL;
  if ((doDelete) && (NonIsoMatchDiTauDr_ != NULL)) delete NonIsoMatchDiTauDr_;
  NonIsoMatchDiTauDr_ = NULL;

  if ((doDelete) && (NonIsoBCSVTurnOn_ != NULL)) delete NonIsoBCSVTurnOn_;
  NonIsoBCSVTurnOn_ = NULL;
  if ((doDelete) && (IsoBCSVTurnOn_ != NULL)) delete IsoBCSVTurnOn_;
  IsoBCSVTurnOn_ = NULL;
  if ((doDelete) && (NonIsoBCSVPassAll_ != NULL)) delete NonIsoBCSVPassAll_;
  NonIsoBCSVPassAll_ = NULL;
  if ((doDelete) && (IsoBCSVPassAll_ != NULL)) delete IsoBCSVPassAll_;
  IsoBCSVPassAll_ = NULL;
  if ((doDelete) && (NonIsoBCSVPassTrig_ != NULL)) delete NonIsoBCSVPassTrig_;
  NonIsoBCSVPassTrig_ = NULL;
  if ((doDelete) && (IsoBCSVPassTrig_ != NULL)) delete IsoBCSVPassTrig_;
  IsoBCSVPassTrig_ = NULL;

  if ((doDelete) && (TrigMatchdR_ != NULL)) delete TrigMatchdR_;
  TrigMatchdR_ = NULL;
}//void AmumuAnalyzer

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void AmumuAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AmumuAnalyzer);

// -*- C++ -*-
//
// Package:    GGHAnalyzer_IndivRECO
// Class:      GGHAnalyzer_IndivRECO
// 
/**\class GGHAnalyzer_IndivRECO GGHAnalyzer_IndivRECO.cc Analyzer/src/GGHAnalyzer_IndivRECO.cc

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

class GGHAnalyzer_IndivRECO : public edm::EDAnalyzer {
   public:
      typedef reco::JetFloatAssociation::Container JetBCEnergyRatioCollection;
      explicit GGHAnalyzer_IndivRECO(const edm::ParameterSet&);
      ~GGHAnalyzer_IndivRECO();

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
      edm::EDGetTokenT<vector<reco::GenParticle> > genParticleTag_;
      edm::EDGetTokenT<vector<reco::PFJet> > akJetTag_;
      edm::EDGetTokenT<vector<reco::Muon> > muonsTag_;
      edm::EDGetTokenT<edm::ValueMap<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > > muonMapTag_;
      edm::EDGetTokenT<edm::ValueMap<reco::PFJetRef> >  jetValMapTag_;
      edm::EDGetTokenT<vector<reco::PFTau> > tauRECOTag_;
      edm::EDGetTokenT<vector<reco::RecoTauPiZero> > pizerosTag_;      
      edm::EDGetTokenT<reco::PFTauDiscriminator>  looseIsoTagRECO_;
      edm::EDGetTokenT<reco::PFTauDiscriminator>  medIsoTagRECO_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> tightIsoTagRECO_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> decayModeFindingTagRECO_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> isoPtSumTagRECO_;
      edm::EDGetTokenT<edm::ValueMap<double> > genMatchedTauVisiblePtMapTagRECO_;
      edm::EDGetTokenT<edm::ValueMap<int> > genMatchedTauDecayModeMapTagRECO_;
      edm::EDGetTokenT<edm::ValueMap<int> > genMatchedTauMatchedMapTagRECO_;
      edm::EDGetTokenT<vector<reco::PFTau> > genMatchedRecoTausRECO_;
      edm::EDGetTokenT<reco::PFTauDiscriminator>  looseIsoTagRECOMVA_;
      edm::EDGetTokenT<reco::PFTauDiscriminator>  medIsoTagRECOMVA_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> tightIsoTagRECOMVA_;


      //Histograms
      TH1F* NEvents_;
      TH1F* MuRemovedJetPt_;
      TH1F* TauMuTauHaddR_;
      TH1F* NConstituentsRECO_;
      TH1F* NTauDecayModeGEN_;
      TH1F* NTauDecayModeRECO_;
      TH1F* GenMatchedIsoPtSumRECO_;
      TH2F* TotalMuvsSoftMu_;
      TH2F* TauHadnConstvsPt_;
      TH2F* RECOPtvsGENVisiblePtRECO_;
      TH2F* RECODecayModevsGENDecayModeRECO_;

      TH1F* MatchedLooseIsoRECOPt_;
      TH1F* MatchedMedIsoRECOPt_;
      TH1F* MatchedTightIsoRECOPt_;
      TH1F* MatchedDMFindRECOPt_;
      TH1F* MatchedRECOPt_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOPt_;
      TGraphAsymmErrors* FinalEffMedIsoRECOPt_;
      TGraphAsymmErrors* FinalEffTightIsoRECOPt_;
      TGraphAsymmErrors* FinalEffDMFindRECOPt_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOPtOverDM_;
      TGraphAsymmErrors* FinalEffMedIsoRECOPtOverDM_;
      TGraphAsymmErrors* FinalEffTightIsoRECOPtOverDM_;

      TH1F* MatchedLooseIsoRECOPtMVA_;
      TH1F* MatchedMedIsoRECOPtMVA_;
      TH1F* MatchedTightIsoRECOPtMVA_;
      TH1F* MatchedRECOPtMVA_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOPtMVA_;
      TGraphAsymmErrors* FinalEffMedIsoRECOPtMVA_;
      TGraphAsymmErrors* FinalEffTightIsoRECOPtMVA_;



      TH1F* MatchedLooseIsoRECOdR_;
      TH1F* MatchedMedIsoRECOdR_;
      TH1F* MatchedTightIsoRECOdR_;
      TH1F* MatchedDMFindRECOdR_;
      TH1F* MatchedRECOdR_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOdR_;
      TGraphAsymmErrors* FinalEffMedIsoRECOdR_;
      TGraphAsymmErrors* FinalEffTightIsoRECOdR_;
      TGraphAsymmErrors* FinalEffDMFindRECOdR_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOdROverDM_;
      TGraphAsymmErrors* FinalEffMedIsoRECOdROverDM_;
      TGraphAsymmErrors* FinalEffTightIsoRECOdROverDM_;
      
      TH1F* MatchedLooseIsoRECOdRMVA_;
      TH1F* MatchedMedIsoRECOdRMVA_;
      TH1F* MatchedTightIsoRECOdRMVA_;
      TH1F* MatchedRECOdRMVA_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOdRMVA_;
      TGraphAsymmErrors* FinalEffMedIsoRECOdRMVA_;
      TGraphAsymmErrors* FinalEffTightIsoRECOdRMVA_;



      TH1F* MatchedLooseIsoRECOPtGen_;
      TH1F* MatchedMedIsoRECOPtGen_;
      TH1F* MatchedTightIsoRECOPtGen_;
      TH1F* MatchedDMFindRECOPtGen_;
      TH1F* MatchedRECOPtGen_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOPtGen_;
      TGraphAsymmErrors* FinalEffMedIsoRECOPtGen_;
      TGraphAsymmErrors* FinalEffTightIsoRECOPtGen_;
      TGraphAsymmErrors* FinalEffDMFindRECOPtGen_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOPtGenOverDM_;
      TGraphAsymmErrors* FinalEffMedIsoRECOPtGenOverDM_;
      TGraphAsymmErrors* FinalEffTightIsoRECOPtGenOverDM_;

      TH1F* MatchedLooseIsoRECOPtGenMVA_;
      TH1F* MatchedMedIsoRECOPtGenMVA_;
      TH1F* MatchedTightIsoRECOPtGenMVA_;
      TH1F* MatchedRECOPtGenMVA_;
      TGraphAsymmErrors* FinalEffLooseIsoRECOPtGenMVA_;
      TGraphAsymmErrors* FinalEffMedIsoRECOPtGenMVA_;
      TGraphAsymmErrors* FinalEffTightIsoRECOPtGenMVA_;



      TH1F* OneProngDMRECOPt_;
      TH1F* OneProngOnePizDMRECOPt_;
      TH1F* OneProngTwoPizDMRECOPt_;
      TH1F* ThreeProngDMRECOPt_;
      TH1F* MatchedOneProngRECOPt_;
      TH1F* MatchedOneProngOnePizRECOPt_;
      TH1F* MatchedOneProngTwoPizRECOPt_;
      TH1F* MatchedThreeProngRECOPt_;
      TGraphAsymmErrors* FinalOneProngDMRECOPt_;
      TGraphAsymmErrors* FinalOneProngOnePizDMRECOPt_;
      TGraphAsymmErrors* FinalOneProngTwoPizDMRECOPt_;
      TGraphAsymmErrors* FinalThreeProngDMRECOPt_;

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
GGHAnalyzer_IndivRECO::GGHAnalyzer_IndivRECO(const edm::ParameterSet& iConfig):
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  genParticleTag_(consumes<vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleTag"))),
  akJetTag_(consumes<vector<reco::PFJet> >(iConfig.getParameter<edm::InputTag>("akJetTag"))),
  muonsTag_(consumes<vector<reco::Muon> >(iConfig.getParameter<edm::InputTag>("muonsTag"))),
  muonMapTag_(consumes<edm::ValueMap<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > >(iConfig.getParameter<edm::InputTag>("muonMapTag"))),
  jetValMapTag_(consumes<edm::ValueMap<reco::PFJetRef> >(iConfig.getParameter<edm::InputTag>("jetValMapTag"))),
  tauRECOTag_(consumes<vector<reco::PFTau> >(iConfig.getParameter<edm::InputTag>("tauRECOTag"))),
  pizerosTag_(consumes<vector<reco::RecoTauPiZero> >(iConfig.getParameter<edm::InputTag>("pizerosTag"))),
  looseIsoTagRECO_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("looseIsoTagRECO"))),
  medIsoTagRECO_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("medIsoTagRECO"))),
  tightIsoTagRECO_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("tightIsoTagRECO"))),
  decayModeFindingTagRECO_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("decayModeFindingTagRECO"))),
  isoPtSumTagRECO_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("isoPtSumTagRECO"))),
  genMatchedTauVisiblePtMapTagRECO_(consumes<edm::ValueMap<double> >(iConfig.getParameter<edm::InputTag>("genMatchedTauVisiblePtMapTagRECO"))),
  genMatchedTauDecayModeMapTagRECO_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("genMatchedTauDecayModeMapTagRECO"))),
  genMatchedTauMatchedMapTagRECO_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("genMatchedTauMatchedMapTagRECO"))),
  genMatchedRecoTausRECO_(consumes<vector<reco::PFTau> >(iConfig.getParameter<edm::InputTag>("genMatchedRecoTausRECO"))),
  looseIsoTagRECOMVA_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("looseIsoTagRECOMVA"))),
  medIsoTagRECOMVA_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("medIsoTagRECOMVA"))),
  tightIsoTagRECOMVA_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("tightIsoTagRECOMVA")))
{
  reset(false);    
}//GGHAnalyzer_IndivRECO



GGHAnalyzer_IndivRECO::~GGHAnalyzer_IndivRECO()
{
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void GGHAnalyzer_IndivRECO::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "\n<------------THIS IS A NEW EVENT------------>" << std::endl;
  NEvents_->Fill(0);

  //Get gen particle collection
  edm::Handle<reco::GenParticleCollection> pGenParts;
  iEvent.getByToken(genParticleTag_, pGenParts);

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
  edm::Handle<std::vector<reco::PFTau> > pTausRECO;
  iEvent.getByToken(tauRECOTag_, pTausRECO);

  //Get pi-Zero particle collection
  edm::Handle<std::vector<reco::RecoTauPiZero> > pPiZero;
  iEvent.getByToken(pizerosTag_, pPiZero);

  //Get Loose Iso Collection
  Handle<PFTauDiscriminator> pLooseIsoDiscRECO; 
  iEvent.getByToken(looseIsoTagRECO_, pLooseIsoDiscRECO); 

  //Get Medium Iso Collection
  Handle<PFTauDiscriminator> pMedIsoDiscRECO; 
  iEvent.getByToken(medIsoTagRECO_, pMedIsoDiscRECO);

  //Get Tight Iso Collection
  Handle<PFTauDiscriminator> pTightIsoDiscRECO; 
  iEvent.getByToken(tightIsoTagRECO_, pTightIsoDiscRECO);

  //Get Decay Mode Finding Collection
  Handle<PFTauDiscriminator> pDMFindingRECO; 
  iEvent.getByToken(decayModeFindingTagRECO_, pDMFindingRECO);

  //Get Loose Iso Collection
  Handle<PFTauDiscriminator> pLooseIsoDiscRECOMVA;
  iEvent.getByToken(looseIsoTagRECOMVA_, pLooseIsoDiscRECOMVA);

  //Get Medium Iso Collection
  Handle<PFTauDiscriminator> pMedIsoDiscRECOMVA;
  iEvent.getByToken(medIsoTagRECOMVA_, pMedIsoDiscRECOMVA);

  //Get Tight Iso Collection
  Handle<PFTauDiscriminator> pTightIsoDiscRECOMVA;
  iEvent.getByToken(tightIsoTagRECOMVA_, pTightIsoDiscRECOMVA);

  //Get Isolation Pt Sum Finding Collection
  Handle<PFTauDiscriminator> pIsoPtSumRECO;
  iEvent.getByToken(isoPtSumTagRECO_, pIsoPtSumRECO);

  //Gen Matched Visible Pt Map
  edm::Handle<edm::ValueMap<double> > pGenMatchedTauVisiblePtMapRECO;
  iEvent.getByToken(genMatchedTauVisiblePtMapTagRECO_, pGenMatchedTauVisiblePtMapRECO);

  //Gen Matched Decay Mode Map
  edm::Handle<edm::ValueMap<int> > pGenMatchedTauDecayModeMapRECO;
  iEvent.getByToken(genMatchedTauDecayModeMapTagRECO_, pGenMatchedTauDecayModeMapRECO);

  //Gen Matched Truth Map
  edm::Handle<edm::ValueMap<int> > pGenMatchedTauMatchedMapRECO;
  iEvent.getByToken(genMatchedTauMatchedMapTagRECO_, pGenMatchedTauMatchedMapRECO);

  //Get RECO Taus particle collection
  edm::Handle<std::vector<reco::PFTau> > pAccessersRECO;
  iEvent.getByToken(genMatchedRecoTausRECO_, pAccessersRECO);

////////////////////
// Find Gen diTaudR
////////////////////
  reco::GenParticleRef tau1Ref, tau2Ref;
  double GEN_diTaudR = 1000;
  for (reco::GenParticleCollection::const_iterator iGenParticle = pGenParts->begin(); iGenParticle != pGenParts->end(); ++iGenParticle)
  {
    if(iGenParticle->pdgId() == 36 && iGenParticle->numberOfDaughters() == 2 && fabs(iGenParticle->daughter(0)->pdgId() ) == 15 )
    {
      tau1Ref = iGenParticle->daughterRef(0);
      tau2Ref = iGenParticle->daughterRef(1);

      //This finds the status 2 taus by seeing if a tau is the daughter of the current tau ref
      while(VariousFunctions::findIfInDaughters(tau1Ref, 15, true) || VariousFunctions::findIfInDaughters(tau2Ref, 15, true))
      {
        if(VariousFunctions::findIfInDaughters(tau1Ref, 15, true))
          tau1Ref = VariousFunctions::findDaughterInDaughters(tau1Ref, 15, true);
        if(VariousFunctions::findIfInDaughters(tau2Ref, 15, true))
          tau2Ref = VariousFunctions::findDaughterInDaughters(tau2Ref, 15, true);
      }//while

/*
      unsigned int iDaughters1 = 0, iDaughters2 = 0;
      std::cout << "TAU1" << std::endl;
      while (iDaughters1 < tau1Ref->numberOfDaughters() )
      {
	int pdgId = abs(tau1Ref->daughter(iDaughters1)->pdgId() ) ;
	std::cout << "\ttau1Ref->daughter(" << iDaughters1 << ")->pdgId()= " << pdgId << std::endl;
	if (pdgId == 213 || pdgId == 223 || pdgId == 20213 || pdgId == 130 || pdgId == 310 || pdgId == 321)
	{
	  unsigned int iDaughters3 = 0;
	  while (iDaughters3 < tau1Ref->daughter(iDaughters1)->numberOfDaughters() )
 	  {
            std::cout << "\ttau1Ref->daughter(" << iDaughters1 << ")->daughters(" << iDaughters3 << ")->pdgId()= " << pdgId << std::endl;
	    iDaughters3++;
	  }//while
	}//if pdgId
        iDaughters1++;
      }//while iDaughters1
      std::cout << "TAU2" << std::endl;
      while (iDaughters2 < tau2Ref->numberOfDaughters() )
      {
        int pdgId = abs(tau2Ref->daughter(iDaughters2)->pdgId() );
        std::cout << "\ttau2Ref->daughter(" << iDaughters2 << ")->pdgId()= " << pdgId << std::endl;
        if (pdgId == 213 || pdgId == 223 || pdgId == 20213 || pdgId == 130 || pdgId == 310 || pdgId == 321)
        {
          unsigned int iDaughters3 = 0;
          while (iDaughters3 < tau2Ref->daughter(iDaughters2)->numberOfDaughters() )
          {
            std::cout << "\t\ttau2Ref->daughter(" << iDaughters2 << ")->daughters(" << iDaughters3 << ")->pdgId()= " << pdgId << std::endl;
            iDaughters3++;
          }//while
        }//if pdgId
        iDaughters2++;
      }//while iDaughters1
*/      

      if ( (VariousFunctions::tauDecayMode(tau1Ref)==7 && VariousFunctions::tauDecayMode(tau2Ref)<5) || (VariousFunctions::tauDecayMode(tau2Ref)==7 && VariousFunctions::tauDecayMode(tau1Ref)<5 ) )
      {
        std::cout << "A TAU_MU TAU_HAD EVENT!!!!" << std::endl;
        NEvents_->Fill(2);
        double dPhi = reco::deltaPhi(tau1Ref->phi(), tau2Ref->phi() ), dEta = tau1Ref->eta() - tau2Ref->eta();
        GEN_diTaudR = sqrt( dEta * dEta + dPhi * dPhi);
        break;
      }//if
      else if (VariousFunctions::tauDecayMode(tau1Ref) == 7 && VariousFunctions::tauDecayMode(tau2Ref) == 7 )
      {
        std::cout << "A TAU_MU TAU_MU EVENT!!!!" << std::endl;
        NEvents_->Fill(6);
        return;
      }//else if
      else
      {
        std::cout << "NOT a TAU_MU TAU_HAD EVENT!!!!" << std::endl;
        NEvents_->Fill(1);
        return;
      }//else  */
    }//it found a1 and it decayed to taus
  }//for iGen Particle

//////////////////////////////
// Looking at CleanJets
//////////////////////////////
  bool muonWasRemoved = false;
  for(size_t iTauRECO = 0; iTauRECO < pAccessersRECO->size(); ++iTauRECO)
  {
    const reco::PFTauRef accesserTauRef(pAccessersRECO, iTauRECO);//
    const reco::PFTauRef hpsTauRef(pTausRECO, iTauRECO);//
    
    const int    GenTauMatchedMapRECO   = (*pGenMatchedTauMatchedMapRECO)[accesserTauRef];
    const int    GenTauDecayModeMapRECO = (*pGenMatchedTauDecayModeMapRECO)[accesserTauRef];
    const double GenTauVisiblePtMapRECO = (*pGenMatchedTauVisiblePtMapRECO)[accesserTauRef];
    const reco::PFJetRef& tauJetRef = (*accesserTauRef).jetRef();
    std::vector<reco::PFCandidatePtr> JetPFCands = tauJetRef->getPFConstituents();
    NConstituentsRECO_->Fill(JetPFCands.size() );

    //////////////////////////////////
    // Making sure a gen match exists
    //////////////////////////////////
    if (GenTauMatchedMapRECO == 1 && fabs(hpsTauRef->eta() ) < 2.3 && hpsTauRef->pt() > 20 ) //Checking that it's gen matched
    {
      int DM = (*pDMFindingRECO)[hpsTauRef], MedIso = (*pMedIsoDiscRECO)[hpsTauRef], LooseIso = (*pLooseIsoDiscRECO)[hpsTauRef], TightIso = (*pTightIsoDiscRECO)[hpsTauRef];
      double PtSum = (*pIsoPtSumRECO)[hpsTauRef];
      GenMatchedIsoPtSumRECO_->Fill(PtSum );
      int MedIsoMVA = (*pMedIsoDiscRECOMVA)[hpsTauRef], LooseIsoMVA = (*pLooseIsoDiscRECOMVA)[hpsTauRef], TightIsoMVA = (*pTightIsoDiscRECOMVA)[hpsTauRef];

      NEvents_->Fill(3);
      MatchedRECOPt_->Fill(hpsTauRef->pt() );
      RECOPtvsGENVisiblePtRECO_->Fill(hpsTauRef->pt(), GenTauVisiblePtMapRECO );
      if (DM == 1)
        MatchedDMFindRECOPt_->Fill(hpsTauRef->pt() );
      if (TightIso == 1 && DM == 1)
        MatchedTightIsoRECOPt_->Fill(hpsTauRef->pt() );
      if (MedIso == 1 && DM == 1)
        MatchedMedIsoRECOPt_->Fill(hpsTauRef->pt() );
      if (LooseIso == 1 && DM == 1)
        MatchedLooseIsoRECOPt_->Fill(hpsTauRef->pt() );
  
      MatchedRECOPtMVA_->Fill(hpsTauRef->pt() );
      if (TightIsoMVA == 1)
        MatchedTightIsoRECOPtMVA_->Fill(hpsTauRef->pt() );
      if (MedIsoMVA == 1)
        MatchedMedIsoRECOPtMVA_->Fill(hpsTauRef->pt() );
      if (LooseIsoMVA == 1)
        MatchedLooseIsoRECOPtMVA_->Fill(hpsTauRef->pt() );
  
  
    
      MatchedRECOdR_->Fill(GEN_diTaudR );
      if (DM == 1)
        MatchedDMFindRECOdR_->Fill(GEN_diTaudR );
      if (TightIso == 1 && DM == 1)
        MatchedTightIsoRECOdR_->Fill(GEN_diTaudR );
      if (MedIso == 1 && DM == 1)
        MatchedMedIsoRECOdR_->Fill(GEN_diTaudR );
      if (LooseIso == 1 && DM == 1)
        MatchedLooseIsoRECOdR_->Fill(GEN_diTaudR );
      
      MatchedRECOdRMVA_->Fill(GEN_diTaudR );
      if (TightIsoMVA == 1)
        MatchedTightIsoRECOdRMVA_->Fill(GEN_diTaudR );
      if (MedIsoMVA == 1)
        MatchedMedIsoRECOdRMVA_->Fill(GEN_diTaudR );
      if (LooseIsoMVA == 1)
        MatchedLooseIsoRECOdRMVA_->Fill(GEN_diTaudR );
      
  
  
      MatchedRECOPtGen_->Fill(GenTauVisiblePtMapRECO );
      if (DM == 1)
        MatchedDMFindRECOPtGen_->Fill(GenTauVisiblePtMapRECO );
      if (TightIso == 1 && DM == 1)
        MatchedTightIsoRECOPtGen_->Fill(GenTauVisiblePtMapRECO );
      if (MedIso == 1 && DM == 1)
        MatchedMedIsoRECOPtGen_->Fill(GenTauVisiblePtMapRECO );
      if (LooseIso == 1 && DM == 1)
        MatchedLooseIsoRECOPtGen_->Fill(GenTauVisiblePtMapRECO );
  
      MatchedRECOPtGenMVA_->Fill(GenTauVisiblePtMapRECO );
      if (TightIsoMVA == 1)
        MatchedTightIsoRECOPtGenMVA_->Fill(GenTauVisiblePtMapRECO );
      if (MedIsoMVA == 1)
        MatchedMedIsoRECOPtGenMVA_->Fill(GenTauVisiblePtMapRECO );
      if (LooseIsoMVA == 1)
        MatchedLooseIsoRECOPtGenMVA_->Fill(GenTauVisiblePtMapRECO );
  
/*      int RECODecayMode = -1;
      if (DM == 1)
      {
         
        std::cout << "\tpTausRECO->at(iTauRECO).decayMode()= " << pTausRECO->at(iTauRECO).decayMode() << "\tGenTauDecayModeMapRECO= " << GenTauDecayModeMapRECO << std::endl;
        std::cout << "\tpTausRECO->at(iTauRECO).signalPFChargedHadrCands().size()= " << pTausRECO->at(iTauRECO).signalPFChargedHadrCands().size() << "\tpTausRECO->at(iTauRECO).signalPFGammaCands().size()= " << pTausRECO->at(iTauRECO).signalPFGammaCands().size() << std::endl;
        if (pTausRECO->at(iTauRECO).decayMode() == 0)
          RECODecayMode = 1;
        if (pTausRECO->at(iTauRECO).decayMode() == 1)
          RECODecayMode = 2;
        if (pTausRECO->at(iTauRECO).decayMode() == 2)
          RECODecayMode = 3;
        if (pTausRECO->at(iTauRECO).decayMode() == 10)
          RECODecayMode = 4;
      }//if DM == 1
*/
      NTauDecayModeRECO_->Fill(pTausRECO->at(iTauRECO).decayMode() ); 

      if (GenTauDecayModeMapRECO == 1)
      {
        RECODecayModevsGENDecayModeRECO_->Fill(pTausRECO->at(iTauRECO).decayMode(), GenTauDecayModeMapRECO); 
        NTauDecayModeGEN_->Fill(1);
        MatchedOneProngRECOPt_->Fill(GenTauVisiblePtMapRECO );
        if (DM == 1 && MedIso == 1)
  	OneProngDMRECOPt_->Fill(GenTauVisiblePtMapRECO );     
      }//if decayMode == 1
      else if (GenTauDecayModeMapRECO == 2)
      {
        RECODecayModevsGENDecayModeRECO_->Fill(pTausRECO->at(iTauRECO).decayMode(), GenTauDecayModeMapRECO); 
        NTauDecayModeGEN_->Fill(2);
        MatchedOneProngOnePizRECOPt_->Fill(GenTauVisiblePtMapRECO );
        if (DM == 1 && MedIso == 1)
          OneProngOnePizDMRECOPt_->Fill(GenTauVisiblePtMapRECO );
      }//if decayMode == 1
      else if (GenTauDecayModeMapRECO == 3)
      {
        RECODecayModevsGENDecayModeRECO_->Fill(pTausRECO->at(iTauRECO).decayMode(), GenTauDecayModeMapRECO); 
        NTauDecayModeGEN_->Fill(3);
        MatchedOneProngTwoPizRECOPt_->Fill(GenTauVisiblePtMapRECO );
        if (DM == 1 && MedIso == 1)
          OneProngTwoPizDMRECOPt_->Fill(GenTauVisiblePtMapRECO );
      }//if decayMode == 1
      else if (GenTauDecayModeMapRECO == 4)
      {
        RECODecayModevsGENDecayModeRECO_->Fill(pTausRECO->at(iTauRECO).decayMode(), GenTauDecayModeMapRECO); 
        NTauDecayModeGEN_->Fill(4);
        MatchedThreeProngRECOPt_->Fill(GenTauVisiblePtMapRECO );
        if (DM == 1 && MedIso == 1)
          ThreeProngDMRECOPt_->Fill(GenTauVisiblePtMapRECO );
      }//if decayMode == 1
    }//if GenTauMatchedMapRECO == 1 && fabs(hpsTauRef->eta() ) < 2.3
  }//if GEN tau1Ref is the had in mu+had and it is matched to CleanJets Je


  if (muonWasRemoved)
    NEvents_->Fill(5);

}//End GGHAnalyzer_IndivRECO::analyze


// ------------ method called once each job just before starting event loop  ------------
void GGHAnalyzer_IndivRECO::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  //Open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //Book histograms
  NEvents_     = new TH1F("NEvents"    , "", 9, -.5, 8.5);
      NEvents_->GetXaxis()->SetBinLabel(1, "TotalEvents");
      NEvents_->GetXaxis()->SetBinLabel(2, "NOT Gen #tau_{#mu} + #tau_{had}");
      NEvents_->GetXaxis()->SetBinLabel(3, "Gen #tau_{#mu} + #tau_{had}");
      NEvents_->GetXaxis()->SetBinLabel(4, "CJ Gen Match #tau_{had}");
      NEvents_->GetXaxis()->SetBinLabel(5, "RECO Gen Match #tau_{had}");
      NEvents_->GetXaxis()->SetBinLabel(6, "Muon was removed");
      NEvents_->GetXaxis()->SetBinLabel(7, "Gen #tau_{#mu} + #tau_{#mu}");
  MuRemovedJetPt_            = new TH1F("MuRemovedJetPt"    , "", 100, -.5, 100);
  TauMuTauHaddR_       = new TH1F("TauMuTauHaddR"    , "", 100, 0, 8);
  NConstituentsRECO_        = new TH1F("NConstituentsRECO"    , "", 50, 0, 50);
  NTauDecayModeGEN_        = new TH1F("NTauDecayModeGEN"    , "", 8, -.5, 7.5);
  NTauDecayModeRECO_        = new TH1F("NTauDecayModeRECO"    , "", 11, -.5, 10.5);
  GenMatchedIsoPtSumRECO_        = new TH1F("GenMatchedIsoPtSumRECO"    , "", 100, 0, 20);
  TotalMuvsSoftMu_  = new TH2F("TotalMuvsSoftMu" , "", 21, -.5, 20.5, 21, -.5, 20.5);
  TauHadnConstvsPt_  = new TH2F("TauHadnConstvsPt" , "", 50, 0, 50, 30, 0, 220);
  RECOPtvsGENVisiblePtRECO_  = new TH2F("RECOPtvsGENVisiblePtRECO" , "", 50, 0, 220, 50, 0, 220);
  RECODecayModevsGENDecayModeRECO_  = new TH2F("RECODecayModevsGENDecayModeRECO" , "", 11, -0.5, 10.5, 4, 0.5, 4.5);



  MatchedLooseIsoRECOPt_  = new TH1F("MatchedLooseIsoRECOPt"    , "", 11, 0, 220.0);
  MatchedMedIsoRECOPt_    = new TH1F("MatchedMedIsoRECOPt", "", 11, 0, 220);
  MatchedTightIsoRECOPt_    = new TH1F("MatchedTightIsoRECOPt", "", 11, 0, 220);
  MatchedDMFindRECOPt_    = new TH1F("MatchedDMFindRECOPt"    , "", 11, 0, 220);
  MatchedRECOPt_          = new TH1F("MatchedRECOPt"    , "", 11, 0, 220);

  MatchedLooseIsoRECOPtMVA_  = new TH1F("MatchedLooseIsoRECOPtMVA"    , "", 11, 0, 220.0);
  MatchedMedIsoRECOPtMVA_    = new TH1F("MatchedMedIsoRECOPtMVA", "", 11, 0, 220);
  MatchedTightIsoRECOPtMVA_    = new TH1F("MatchedTightIsoRECOPtMVA", "", 11, 0, 220);
  MatchedRECOPtMVA_          = new TH1F("MatchedRECOPtMVA"    , "", 11, 0, 220);



  FinalEffLooseIsoRECOPt_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOPt_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOPt_ = new TGraphAsymmErrors(30);
  FinalEffDMFindRECOPt_ = new TGraphAsymmErrors(30);
  FinalEffLooseIsoRECOPtOverDM_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOPtOverDM_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOPtOverDM_ = new TGraphAsymmErrors(30);

  FinalEffLooseIsoRECOPtMVA_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOPtMVA_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOPtMVA_ = new TGraphAsymmErrors(30);



  MatchedLooseIsoRECOdR_  = new TH1F("MatchedLooseIsoRECOdR"    , "", 8, 0, 1.5);
  MatchedMedIsoRECOdR_    = new TH1F("MatchedMedIsoRECOdR", "", 8, 0, 1.5);
  MatchedTightIsoRECOdR_    = new TH1F("MatchedTightIsoRECOdR", "", 8, 0, 1.5);
  MatchedDMFindRECOdR_    = new TH1F("MatchedDMFindRECOdR"    , "", 8, 0, 1.5);
  MatchedRECOdR_          = new TH1F("MatchedRECOdR"    , "", 8, 0, 1.5);

  MatchedLooseIsoRECOdRMVA_  = new TH1F("MatchedLooseIsoRECOdRMVA"    , "", 8, 0, 1.5);
  MatchedMedIsoRECOdRMVA_    = new TH1F("MatchedMedIsoRECOdRMVA", "", 8, 0, 1.5);
  MatchedTightIsoRECOdRMVA_    = new TH1F("MatchedTightIsoRECOdRMVA", "", 8, 0, 1.5);
  MatchedRECOdRMVA_          = new TH1F("MatchedRECOdRMVA"    , "", 8, 0, 1.5);



  FinalEffLooseIsoRECOdR_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOdR_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOdR_ = new TGraphAsymmErrors(30);
  FinalEffDMFindRECOdR_ = new TGraphAsymmErrors(30);
  FinalEffLooseIsoRECOdROverDM_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOdROverDM_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOdROverDM_ = new TGraphAsymmErrors(30);

  FinalEffLooseIsoRECOdRMVA_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOdRMVA_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOdRMVA_ = new TGraphAsymmErrors(30);



  MatchedLooseIsoRECOPtGen_  = new TH1F("MatchedLooseIsoRECOPtGen"    , "", 11, 0, 220.0);
  MatchedMedIsoRECOPtGen_    = new TH1F("MatchedMedIsoRECOPtGen", "", 11, 0, 220);
  MatchedTightIsoRECOPtGen_    = new TH1F("MatchedTightIsoRECOPtGen", "", 11, 0, 220);
  MatchedDMFindRECOPtGen_    = new TH1F("MatchedDMFindRECOPtGen"    , "", 11, 0, 220);
  MatchedRECOPtGen_          = new TH1F("MatchedRECOPtGen"    , "", 11, 0, 220);
  
  MatchedLooseIsoRECOPtGenMVA_  = new TH1F("MatchedLooseIsoRECOPtGenMVA"    , "", 11, 0, 220.0);
  MatchedMedIsoRECOPtGenMVA_    = new TH1F("MatchedMedIsoRECOPtGenMVA", "", 11, 0, 220);
  MatchedTightIsoRECOPtGenMVA_    = new TH1F("MatchedTightIsoRECOPtGenMVA", "", 11, 0, 220);
  MatchedRECOPtGenMVA_          = new TH1F("MatchedRECOPtGenMVA"    , "", 11, 0, 220);
  


  FinalEffLooseIsoRECOPtGen_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOPtGen_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOPtGen_ = new TGraphAsymmErrors(30);
  FinalEffDMFindRECOPtGen_ = new TGraphAsymmErrors(30);
  FinalEffLooseIsoRECOPtGenOverDM_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOPtGenOverDM_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOPtGenOverDM_ = new TGraphAsymmErrors(30);
  
  FinalEffLooseIsoRECOPtGenMVA_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoRECOPtGenMVA_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoRECOPtGenMVA_ = new TGraphAsymmErrors(30);


  OneProngDMRECOPt_ = new TH1F("OneProngDMRECOPt"    , "", 11, 0, 220.0);
  OneProngOnePizDMRECOPt_ = new TH1F("OneProngOnePizDMRECOPt"    , "", 11, 0, 220.0);
  OneProngTwoPizDMRECOPt_ = new TH1F("OneProngTwoPizDMRECOPt"    , "", 11, 0, 220.0);
  ThreeProngDMRECOPt_ = new TH1F("ThreeProngDMRECOPt"    , "", 11, 0, 220.0);
  MatchedOneProngRECOPt_ = new TH1F("MatchedOneProngRECOPt", "", 11, 0, 220.0);
  MatchedOneProngOnePizRECOPt_ = new TH1F("MatchedOneProngOnePizRECOPt", "", 11, 0, 220.0);
  MatchedOneProngTwoPizRECOPt_ = new TH1F("MatchedOneProngTwoPizRECOPt", "", 11, 0, 220.0);
  MatchedThreeProngRECOPt_ = new TH1F("MatchedThreeProngRECOPt", "", 11, 0, 220.0);
  FinalOneProngDMRECOPt_ = new TGraphAsymmErrors(30);
  FinalOneProngOnePizDMRECOPt_ = new TGraphAsymmErrors(30);
  FinalOneProngTwoPizDMRECOPt_ = new TGraphAsymmErrors(30);
  FinalThreeProngDMRECOPt_ = new TGraphAsymmErrors(30);
}

// ------------ method called once each job just after ending the event loop  ------------
void GGHAnalyzer_IndivRECO::endJob()
{
  //Make the Canvases
  TCanvas NEventsCanvas("NEvents","",600,600);
  TCanvas MuRemovedJetPtCanvas("MuRemovedJetPt","",600,600);
  TCanvas TauMuTauHaddRCanvas("TauMuTauHaddR","",600,600);
  TCanvas NConstituentsRECOCanvas("NConstituentsRECO","",600,600);
  TCanvas NTauDecayModeGENCanvas("NTauDecayModeGEN","",600,600);
  TCanvas NTauDecayModeRECOCanvas("NTauDecayModeRECO","",600,600);
  TCanvas GenMatchedIsoPtSumRECOCanvas("GenMatchedIsoPtSumRECO","",600,600);
  TCanvas TotalMuvsSoftMuCanvas("TotalMuvsSoftMu","",600,600);
  TCanvas TauHadnConstvsPtCanvas("TauHadnConstvsPt","",600,600);
  TCanvas RECOPtvsGENVisiblePtRECOCanvas("RECOPtvsGENVisiblePtRECO","",600,600);
  TCanvas RECODecayModevsGENDecayModeRECOCanvas("RECODecayModevsGENDecayModeRECO","",600,600);

/*  for (int x = 0; x < 4; x++)
  {
    int sum = 0;
    for (int y = 0; y < 4; y++)
      sum += RECODecayModevsGENDecayModeRECO_->GetBinContent(x,y);
    for (int y = 0; y < 4; y++)
    {
      std::cout << "Before(" << x << "," << y << ")= " << RECODecayModevsGENDecayModeRECO_->GetBinContent(x,y) << "\tsum= " << sum;
      std::cout << "\tRECODecayModevsGENDecayModeRECO_->GetBinContent(x,y) / sum = " <<  RECODecayModevsGENDecayModeRECO_->GetBinContent(x,y) / sum << std::endl;
      RECODecayModevsGENDecayModeRECO_->SetBinContent(x, y, RECODecayModevsGENDecayModeRECO_->GetBinContent(x,y) / sum );
      std::cout << "After (" << x << "," << y << ")= " << RECODecayModevsGENDecayModeRECO_->GetBinContent(x,y) << std::endl;
    }//for
  }//for x
*/  


  TCanvas MatchedLooseIsoRECOPtCanvas("MatchedLooseIsoRECOPt","",600,600);
  TCanvas MatchedMedIsoRECOPtCanvas("MatchedMedIsoRECOPt","",600,600);
  TCanvas MatchedTightIsoRECOPtCanvas("MatchedTightIsoRECOPt","",600,600);
  TCanvas MatchedDMFindRECOPtCanvas("MatchedDMFindRECOPt","",600,600);
  TCanvas MatchedRECOPtCanvas("MatchedRECOPt","",600,600);
  TCanvas FinalEffLooseIsoRECOPtCanvas("FinalEffLooseIsoRECOPt","",600,600);
  TCanvas FinalEffMedIsoRECOPtCanvas("FinalEffMedIsoRECOPt","",600,600);
  TCanvas FinalEffTightIsoRECOPtCanvas("FinalEffTightIsoRECOPt","",600,600);
  TCanvas FinalEffDMFindRECOPtCanvas("FinalEffDMFindRECOPt","",600,600);
  TCanvas FinalEffLooseIsoRECOPtOverDMCanvas("FinalEffLooseIsoRECOPtOverDM","",600,600);
  TCanvas FinalEffMedIsoRECOPtOverDMCanvas("FinalEffMedIsoRECOPtOverDM","",600,600);
  TCanvas FinalEffTightIsoRECOPtOverDMCanvas("FinalEffTightIsoRECOPtOverDM","",600,600);

  TCanvas MatchedLooseIsoRECOPtMVACanvas("MatchedLooseIsoRECOPtMVA","",600,600);
  TCanvas MatchedMedIsoRECOPtMVACanvas("MatchedMedIsoRECOPtMVA","",600,600);
  TCanvas MatchedTightIsoRECOPtMVACanvas("MatchedTightIsoRECOPtMVA","",600,600);
  TCanvas MatchedRECOPtMVACanvas("MatchedRECOPtMVA","",600,600);
  TCanvas FinalEffLooseIsoRECOPtMVACanvas("FinalEffLooseIsoRECOPtMVA","",600,600);
  TCanvas FinalEffMedIsoRECOPtMVACanvas("FinalEffMedIsoRECOPtMVA","",600,600);
  TCanvas FinalEffTightIsoRECOPtMVACanvas("FinalEffTightIsoRECOPtMVA","",600,600);



  TCanvas MatchedLooseIsoRECOdRCanvas("MatchedLooseIsoRECOdR","",600,600);
  TCanvas MatchedMedIsoRECOdRCanvas("MatchedMedIsoRECOdR","",600,600);
  TCanvas MatchedTightIsoRECOdRCanvas("MatchedTightIsoRECOdR","",600,600);
  TCanvas MatchedDMFindRECOdRCanvas("MatchedDMFindRECOdR","",600,600);
  TCanvas MatchedRECOdRCanvas("MatchedRECOdR","",600,600);
  TCanvas FinalEffLooseIsoRECOdRCanvas("FinalEffLooseIsoRECOdR","",600,600);
  TCanvas FinalEffMedIsoRECOdRCanvas("FinalEffMedIsoRECOdR","",600,600);
  TCanvas FinalEffTightIsoRECOdRCanvas("FinalEffTightIsoRECOdR","",600,600);
  TCanvas FinalEffDMFindRECOdRCanvas("FinalEffDMFindRECOdR","",600,600);
  TCanvas FinalEffLooseIsoRECOdROverDMCanvas("FinalEffLooseIsoRECOdROverDM","",600,600);
  TCanvas FinalEffMedIsoRECOdROverDMCanvas("FinalEffMedIsoRECOdROverDM","",600,600);
  TCanvas FinalEffTightIsoRECOdROverDMCanvas("FinalEffTightIsoRECOdROverDM","",600,600);
  
  TCanvas MatchedLooseIsoRECOdRMVACanvas("MatchedLooseIsoRECOdRMVA","",600,600);
  TCanvas MatchedMedIsoRECOdRMVACanvas("MatchedMedIsoRECOdRMVA","",600,600);
  TCanvas MatchedTightIsoRECOdRMVACanvas("MatchedTightIsoRECOdRMVA","",600,600);
  TCanvas MatchedRECOdRMVACanvas("MatchedRECOdRMVA","",600,600);
  TCanvas FinalEffLooseIsoRECOdRMVACanvas("FinalEffLooseIsoRECOdRMVA","",600,600);
  TCanvas FinalEffMedIsoRECOdRMVACanvas("FinalEffMedIsoRECOdRMVA","",600,600);
  TCanvas FinalEffTightIsoRECOdRMVACanvas("FinalEffTightIsoRECOdRMVA","",600,600);



  TCanvas MatchedLooseIsoRECOPtGenCanvas("MatchedLooseIsoRECOPtGen","",600,600);
  TCanvas MatchedMedIsoRECOPtGenCanvas("MatchedMedIsoRECOPtGen","",600,600);
  TCanvas MatchedTightIsoRECOPtGenCanvas("MatchedTightIsoRECOPtGen","",600,600);
  TCanvas MatchedDMFindRECOPtGenCanvas("MatchedDMFindRECOPtGen","",600,600);
  TCanvas MatchedRECOPtGenCanvas("MatchedRECOPtGen","",600,600);
  TCanvas FinalEffLooseIsoRECOPtGenCanvas("FinalEffLooseIsoRECOPtGen","",600,600);
  TCanvas FinalEffMedIsoRECOPtGenCanvas("FinalEffMedIsoRECOPtGen","",600,600);
  TCanvas FinalEffTightIsoRECOPtGenCanvas("FinalEffTightIsoRECOPtGen","",600,600);
  TCanvas FinalEffDMFindRECOPtGenCanvas("FinalEffDMFindRECOPtGen","",600,600);
  TCanvas FinalEffLooseIsoRECOPtGenOverDMCanvas("FinalEffLooseIsoRECOPtGenOverDM","",600,600);
  TCanvas FinalEffMedIsoRECOPtGenOverDMCanvas("FinalEffMedIsoRECOPtGenOverDM","",600,600);
  TCanvas FinalEffTightIsoRECOPtGenOverDMCanvas("FinalEffTightIsoRECOPtGenOverDM","",600,600);

  TCanvas MatchedLooseIsoRECOPtGenMVACanvas("MatchedLooseIsoRECOPtGenMVA","",600,600);
  TCanvas MatchedMedIsoRECOPtGenMVACanvas("MatchedMedIsoRECOPtGenMVA","",600,600);
  TCanvas MatchedTightIsoRECOPtGenMVACanvas("MatchedTightIsoRECOPtGenMVA","",600,600);
  TCanvas MatchedRECOPtGenMVACanvas("MatchedRECOPtGenMVA","",600,600);
  TCanvas FinalEffLooseIsoRECOPtGenMVACanvas("FinalEffLooseIsoRECOPtGenMVA","",600,600);
  TCanvas FinalEffMedIsoRECOPtGenMVACanvas("FinalEffMedIsoRECOPtGenMVA","",600,600);
  TCanvas FinalEffTightIsoRECOPtGenMVACanvas("FinalEffTightIsoRECOPtGenMVA","",600,600);
  
 


  TCanvas OneProngDMRECOPtCanvas("OneProngDMRECOPt","",600,600);
  TCanvas OneProngOnePizDMRECOPtCanvas("OneProngOnePizDMRECOPt","",600,600);
  TCanvas OneProngTwoPizDMRECOPtCanvas("OneProngTwoPizDMRECOPt","",600,600);
  TCanvas ThreeProngDMRECOPtCanvas("ThreeProngDMRECOPt","",600,600);
  TCanvas MatchedOneProngRECOPtCanvas("MatchedOneProngRECOPtCanvas","",600,600);
  TCanvas MatchedOneProngOnePizRECOPtCanvas("MatchedOneProngOnePizRECOPtCanvas","",600,600);
  TCanvas MatchedOneProngTwoPizRECOPtCanvas("MatchedOneProngTwoPizRECOPtCanvas","",600,600);
  TCanvas MatchedThreeProngRECOPtCanvas("MatchedThreeProngRECOPtCanvas","",600,600);
  TCanvas FinalOneProngDMRECOPtCanvas("FinalOneProngDMRECOPt","",600,600);
  TCanvas FinalOneProngOnePizDMRECOPtCanvas("FinalOneProngOnePizDMRECOPt","",600,600);
  TCanvas FinalOneProngTwoPizDMRECOPtCanvas("FinalOneProngTwoPizDMRECOPt","",600,600);
  TCanvas FinalThreeProngDMRECOPtCanvas("FinalThreeProngDMRECOPt","",600,600);


  FinalEffLooseIsoRECOPt_->Divide(MatchedLooseIsoRECOPt_, MatchedRECOPt_);
  FinalEffMedIsoRECOPt_->Divide(MatchedMedIsoRECOPt_,     MatchedRECOPt_);
  FinalEffTightIsoRECOPt_->Divide(MatchedTightIsoRECOPt_, MatchedRECOPt_);
  FinalEffDMFindRECOPt_->Divide(MatchedDMFindRECOPt_,     MatchedRECOPt_);
  FinalEffLooseIsoRECOPtOverDM_->Divide(MatchedLooseIsoRECOPt_, MatchedDMFindRECOPt_);
  FinalEffMedIsoRECOPtOverDM_->Divide(MatchedMedIsoRECOPt_,     MatchedDMFindRECOPt_);
  FinalEffTightIsoRECOPtOverDM_->Divide(MatchedTightIsoRECOPt_, MatchedDMFindRECOPt_);

  FinalEffLooseIsoRECOPtMVA_->Divide(MatchedLooseIsoRECOPtMVA_, MatchedRECOPtMVA_);
  FinalEffMedIsoRECOPtMVA_->Divide(MatchedMedIsoRECOPtMVA_,     MatchedRECOPtMVA_);
  FinalEffTightIsoRECOPtMVA_->Divide(MatchedTightIsoRECOPtMVA_, MatchedRECOPtMVA_);
  


  FinalEffLooseIsoRECOdR_->Divide(MatchedLooseIsoRECOdR_, MatchedRECOdR_);
  FinalEffMedIsoRECOdR_->Divide(MatchedMedIsoRECOdR_,     MatchedRECOdR_);
  FinalEffTightIsoRECOdR_->Divide(MatchedTightIsoRECOdR_, MatchedRECOdR_);
  FinalEffDMFindRECOdR_->Divide(MatchedDMFindRECOdR_,     MatchedRECOdR_);
  FinalEffLooseIsoRECOdROverDM_->Divide(MatchedLooseIsoRECOdR_, MatchedDMFindRECOdR_);
  FinalEffMedIsoRECOdROverDM_->Divide(MatchedMedIsoRECOdR_,     MatchedDMFindRECOdR_);
  FinalEffTightIsoRECOdROverDM_->Divide(MatchedTightIsoRECOdR_, MatchedDMFindRECOdR_);
  
  FinalEffLooseIsoRECOdRMVA_->Divide(MatchedLooseIsoRECOdRMVA_, MatchedRECOdRMVA_);
  FinalEffMedIsoRECOdRMVA_->Divide(MatchedMedIsoRECOdRMVA_,     MatchedRECOdRMVA_);
  FinalEffTightIsoRECOdRMVA_->Divide(MatchedTightIsoRECOdRMVA_, MatchedRECOdRMVA_);



  FinalEffLooseIsoRECOPtGen_->Divide(MatchedLooseIsoRECOPtGen_, MatchedRECOPtGen_);
  FinalEffMedIsoRECOPtGen_->Divide(MatchedMedIsoRECOPtGen_,     MatchedRECOPtGen_);
  FinalEffTightIsoRECOPtGen_->Divide(MatchedTightIsoRECOPtGen_, MatchedRECOPtGen_);
  FinalEffDMFindRECOPtGen_->Divide(MatchedDMFindRECOPtGen_,     MatchedRECOPtGen_);
  FinalEffLooseIsoRECOPtGenOverDM_->Divide(MatchedLooseIsoRECOPtGen_, MatchedDMFindRECOPtGen_);
  FinalEffMedIsoRECOPtGenOverDM_->Divide(MatchedMedIsoRECOPtGen_,     MatchedDMFindRECOPtGen_);
  FinalEffTightIsoRECOPtGenOverDM_->Divide(MatchedTightIsoRECOPtGen_, MatchedDMFindRECOPtGen_);

  FinalEffLooseIsoRECOPtGenMVA_->Divide(MatchedLooseIsoRECOPtGenMVA_, MatchedRECOPtGenMVA_);
  FinalEffMedIsoRECOPtGenMVA_->Divide(MatchedMedIsoRECOPtGenMVA_,     MatchedRECOPtGenMVA_);
  FinalEffTightIsoRECOPtGenMVA_->Divide(MatchedTightIsoRECOPtGenMVA_, MatchedRECOPtGenMVA_);



  FinalOneProngDMRECOPt_->Divide(OneProngDMRECOPt_, MatchedOneProngRECOPt_);
  FinalOneProngOnePizDMRECOPt_->Divide(OneProngOnePizDMRECOPt_, MatchedOneProngOnePizRECOPt_);
  FinalOneProngTwoPizDMRECOPt_->Divide(OneProngTwoPizDMRECOPt_, MatchedOneProngTwoPizRECOPt_);
  FinalThreeProngDMRECOPt_->Divide(ThreeProngDMRECOPt_, MatchedThreeProngRECOPt_);

std::cout << "<----------------Declared Canvases-------------->" << std::endl;

  //Format the 1D plots and Draw (canvas, hist, grid, log y, log z, color, size, style, xAxisTitle, xTitleSize, xLabelSize, xTitleOffSet, yAxisTitle, yTitleSize, yLabelSize, yTitleOffset)
  VariousFunctions::formatAndDrawCanvasAndHist1D(NEventsCanvas, NEvents_,
         1, 0, 0, kBlack, 7, 20, "", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MuRemovedJetPtCanvas, MuRemovedJetPt_,
         1, 0, 0, kBlack, 7, 20, "Pt_{Jet with #mu Removed}", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauMuTauHaddRCanvas, TauMuTauHaddR_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(#tau_{mu} + #tau_{H})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NConstituentsRECOCanvas, NConstituentsRECO_,
         1, 0, 0, kBlack, 7, 20, "Number of Constituents", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NTauDecayModeGENCanvas, NTauDecayModeGEN_,
         1, 0, 0, kBlack, 7, 20, "TauDecayMode GEN NoCJ", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NTauDecayModeRECOCanvas, NTauDecayModeRECO_,
         1, 0, 0, kBlack, 7, 20, "TauDecayMode RECO NoCJ", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(GenMatchedIsoPtSumRECOCanvas, GenMatchedIsoPtSumRECO_,
         1, 0, 0, kBlack, 7, 20, "dR(Gen  #tau#tau), if Mu is Removed", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist2D(TotalMuvsSoftMuCanvas, TotalMuvsSoftMu_,
         1, 0, 0, kBlack, 7, 20, "nRECO #tau's", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(TauHadnConstvsPtCanvas, TauHadnConstvsPt_,
         1, 0, 0, kBlack, 7, 20, "# Constituents(#tau_{H})", .04, .04, 1.1, "Pt(#tau_{H})", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(RECOPtvsGENVisiblePtRECOCanvas, RECOPtvsGENVisiblePtRECO_,
         1, 0, 0, kBlack, 7, 20, "RECO Pt(#tau)", .04, .04, 1.1, "GEN Visible Pt(#tau)", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(RECODecayModevsGENDecayModeRECOCanvas, RECODecayModevsGENDecayModeRECO_,
         1, 0, 0, kBlack, 7, 20, "RECO Decay Mode", .04, .04, 1.1, "GEN Decay Mode", .04, .04, 1.6, "", .04, .04, 1.0);



  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoRECOPtCanvas, MatchedLooseIsoRECOPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoRECOPtCanvas, MatchedMedIsoRECOPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoRECOPtCanvas, MatchedTightIsoRECOPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedDMFindRECOPtCanvas, MatchedDMFindRECOPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedRECOPtCanvas, MatchedRECOPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoRECOPtCanvas, FinalEffLooseIsoRECOPt_,
	 1, 0, 0, kBlack, 1, 20, "Pt(CleanJets Matched + Loose Iso + DM) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoRECOPtCanvas, FinalEffMedIsoRECOPt_,
	 1, 0, 0, kBlack, 1, 20, "Pt(CleanJets Matched + Med Iso + DM) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoRECOPtCanvas, FinalEffTightIsoRECOPt_,
	 1, 0, 0, kBlack, 1, 20, "Pt(CleanJets Matched + Tight Iso + DM) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffDMFindRECOPtCanvas, FinalEffDMFindRECOPt_,
	 1, 0, 0, kBlack, 1, 20, "Pt(CleanJets Matched + DecayModeFinding) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoRECOPtOverDMCanvas, FinalEffLooseIsoRECOPtOverDM_,
         1, 0, 0, kBlack, 1, 20, "Pt(CleanJets Matched + Loose Iso + DM) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoRECOPtOverDMCanvas, FinalEffMedIsoRECOPtOverDM_,
         1, 0, 0, kBlack, 1, 20, "Pt(CleanJets Matched + Med Iso + DM) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoRECOPtOverDMCanvas, FinalEffTightIsoRECOPtOverDM_,
         1, 0, 0, kBlack, 1, 20, "Pt(CleanJets Matched + Tight Iso + DM) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoRECOPtMVACanvas, MatchedLooseIsoRECOPtMVA_,
         1, 0, 0, kBlack, 7, 20, "PtMVA(CleanJets Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoRECOPtMVACanvas, MatchedMedIsoRECOPtMVA_,
         1, 0, 0, kBlack, 7, 20, "PtMVA(CleanJets Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoRECOPtMVACanvas, MatchedTightIsoRECOPtMVA_,
         1, 0, 0, kBlack, 7, 20, "PtMVA(CleanJets Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedRECOPtMVACanvas, MatchedRECOPtMVA_,
         1, 0, 0, kBlack, 7, 20, "PtMVA(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoRECOPtMVACanvas, FinalEffLooseIsoRECOPtMVA_,
         1, 0, 0, kBlack, 1, 20, "PtMVA(CleanJets Matched + Loose Iso + DM) / PtMVA(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoRECOPtMVACanvas, FinalEffMedIsoRECOPtMVA_,
         1, 0, 0, kBlack, 1, 20, "PtMVA(CleanJets Matched + Med Iso + DM) / PtMVA(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoRECOPtMVACanvas, FinalEffTightIsoRECOPtMVA_,
         1, 0, 0, kBlack, 1, 20, "PtMVA(CleanJets Matched + Tight Iso + DM) / PtMVA(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);




  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoRECOdRCanvas, MatchedLooseIsoRECOdR_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoRECOdRCanvas, MatchedMedIsoRECOdR_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoRECOdRCanvas, MatchedTightIsoRECOdR_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedDMFindRECOdRCanvas, MatchedDMFindRECOdR_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedRECOdRCanvas, MatchedRECOdR_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoRECOdRCanvas, FinalEffLooseIsoRECOdR_,
	 1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Loose Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoRECOdRCanvas, FinalEffMedIsoRECOdR_,
	 1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Med Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoRECOdRCanvas, FinalEffTightIsoRECOdR_,
	 1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Tight Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffDMFindRECOdRCanvas, FinalEffDMFindRECOdR_,
	 1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + DecayModeFinding) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoRECOdROverDMCanvas, FinalEffLooseIsoRECOdROverDM_,
         1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Loose Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoRECOdROverDMCanvas, FinalEffMedIsoRECOdROverDM_,
         1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Med Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoRECOdROverDMCanvas, FinalEffTightIsoRECOdROverDM_,
         1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Tight Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoRECOdRMVACanvas, MatchedLooseIsoRECOdRMVA_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoRECOdRMVACanvas, MatchedMedIsoRECOdRMVA_, 
         1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoRECOdRMVACanvas, MatchedTightIsoRECOdRMVA_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedRECOdRMVACanvas, MatchedRECOdRMVA_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoRECOdRMVACanvas, FinalEffLooseIsoRECOdRMVA_,
         1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Loose Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoRECOdRMVACanvas, FinalEffMedIsoRECOdRMVA_,
         1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Med Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoRECOdRMVACanvas, FinalEffTightIsoRECOdRMVA_,
         1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Tight Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);



  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoRECOPtGenCanvas, MatchedLooseIsoRECOPtGen_,
         1, 0, 0, kBlack, 7, 20, "PtGen(CleanJets Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoRECOPtGenCanvas, MatchedMedIsoRECOPtGen_,
         1, 0, 0, kBlack, 7, 20, "PtGen(CleanJets Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoRECOPtGenCanvas, MatchedTightIsoRECOPtGen_,
         1, 0, 0, kBlack, 7, 20, "PtGen(CleanJets Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedDMFindRECOPtGenCanvas, MatchedDMFindRECOPtGen_,
         1, 0, 0, kBlack, 7, 20, "PtGen(CleanJets Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedRECOPtGenCanvas, MatchedRECOPtGen_,
         1, 0, 0, kBlack, 7, 20, "PtGen(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoRECOPtGenCanvas, FinalEffLooseIsoRECOPtGen_,
         1, 0, 0, kBlack, 1, 20, "PtGen(CleanJets Matched + Loose Iso + DM) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoRECOPtGenCanvas, FinalEffMedIsoRECOPtGen_,
         1, 0, 0, kBlack, 1, 20, "PtGen(CleanJets Matched + Med Iso + DM) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoRECOPtGenCanvas, FinalEffTightIsoRECOPtGen_,
         1, 0, 0, kBlack, 1, 20, "PtGen(CleanJets Matched + Tight Iso + DM) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffDMFindRECOPtGenCanvas, FinalEffDMFindRECOPtGen_,
         1, 0, 0, kBlack, 1, 20, "PtGen(CleanJets Matched + DecayModeFinding) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoRECOPtGenOverDMCanvas, FinalEffLooseIsoRECOPtGenOverDM_,
         1, 0, 0, kBlack, 1, 20, "PtGen(CleanJets Matched + Loose Iso + DM) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoRECOPtGenOverDMCanvas, FinalEffMedIsoRECOPtGenOverDM_,
         1, 0, 0, kBlack, 1, 20, "PtGen(CleanJets Matched + Med Iso + DM) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoRECOPtGenOverDMCanvas, FinalEffTightIsoRECOPtGenOverDM_,
         1, 0, 0, kBlack, 1, 20, "PtGen(CleanJets Matched + Tight Iso + DM) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoRECOPtGenMVACanvas, MatchedLooseIsoRECOPtGenMVA_,
         1, 0, 0, kBlack, 7, 20, "PtGenMVA(CleanJets Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoRECOPtGenMVACanvas, MatchedMedIsoRECOPtGenMVA_,
         1, 0, 0, kBlack, 7, 20, "PtGenMVA(CleanJets Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoRECOPtGenMVACanvas, MatchedTightIsoRECOPtGenMVA_,
         1, 0, 0, kBlack, 7, 20, "PtGenMVA(CleanJets Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedRECOPtGenMVACanvas, MatchedRECOPtGenMVA_,
         1, 0, 0, kBlack, 7, 20, "PtGenMVA(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoRECOPtGenMVACanvas, FinalEffLooseIsoRECOPtGenMVA_,
         1, 0, 0, kBlack, 1, 20, "PtGenMVA(CleanJets Matched + Loose Iso + DM) / PtGenMVA(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoRECOPtGenMVACanvas, FinalEffMedIsoRECOPtGenMVA_,
         1, 0, 0, kBlack, 1, 20, "PtGenMVA(CleanJets Matched + Med Iso + DM) / PtGenMVA(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoRECOPtGenMVACanvas, FinalEffTightIsoRECOPtGenMVA_,
         1, 0, 0, kBlack, 1, 20, "PtGenMVA(CleanJets Matched + Tight Iso + DM) / PtGenMVA(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);


  
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngDMRECOPtCanvas, OneProngDMRECOPt_, 
         1, 0, 0, kBlack, 7, 20, "1 Prong p_{T}(CleanJets Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngOnePizDMRECOPtCanvas, OneProngOnePizDMRECOPt_, 
         1, 0, 0, kBlack, 7, 20, "1 Prong + 1 #pi^{0} p_{T}(CleanJets Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngTwoPizDMRECOPtCanvas, OneProngTwoPizDMRECOPt_, 
         1, 0, 0, kBlack, 7, 20, "1 Prong + 2 #pi^{0} p_{T}(CleanJets Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(ThreeProngDMRECOPtCanvas, ThreeProngDMRECOPt_, 
         1, 0, 0, kBlack, 7, 20, "3 Prong p_{T}(CleanJets Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedOneProngRECOPtCanvas, MatchedOneProngRECOPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong p_{T}(CleanJets Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedOneProngOnePizRECOPtCanvas, MatchedOneProngOnePizRECOPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 1 #pi^{0} p_{T}(CleanJets Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedOneProngTwoPizRECOPtCanvas, MatchedOneProngTwoPizRECOPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 2 #pi^{0} p_{T}(CleanJets Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedThreeProngRECOPtCanvas, MatchedThreeProngRECOPt_,
         1, 0, 0, kBlack, 7, 20, "3 Prong p_{T}(CleanJets Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngDMRECOPtCanvas, FinalOneProngDMRECOPt_, 
         1, 0, 0, kBlack, 1, 20, "1 Prong p_{T}(CleanJets Matched + DM) / 1 Prong p_{T}(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngOnePizDMRECOPtCanvas, FinalOneProngOnePizDMRECOPt_, 
         1, 0, 0, kBlack, 1, 20, "1 Prong + 1 #pi^{0} p_{T}(CleanJets Matched + DM) / 1 Prong + 1 #pi^{0} p_{T}(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngTwoPizDMRECOPtCanvas, FinalOneProngTwoPizDMRECOPt_, 
         1, 0, 0, kBlack, 1, 20, "1 Prong + 2 #pi^{0} p_{T}(CleanJets Matched + DM) / 1 Prong + 2 #pi^{0} p_{T}(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalThreeProngDMRECOPtCanvas, FinalThreeProngDMRECOPt_, 
         1, 0, 0, kBlack, 1, 20, "3 Prong p_{T}(CleanJets Matched + DM) / 3 Prong p_{T}(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);



std::cout << "after formatting" << std::endl;
std::cout << "<----------------Formatted Canvases and Histos-------------->" << std::endl;

  //Write output file
  out_->cd();

  NEventsCanvas.Write();
  MuRemovedJetPtCanvas.Write();
  TauMuTauHaddRCanvas.Write();
  NConstituentsRECOCanvas.Write();
  NTauDecayModeGENCanvas.Write();
  NTauDecayModeRECOCanvas.Write();
  GenMatchedIsoPtSumRECOCanvas.Write();
  TotalMuvsSoftMuCanvas.Write();
  TauHadnConstvsPtCanvas.Write();
  RECOPtvsGENVisiblePtRECOCanvas.Write();
  RECODecayModevsGENDecayModeRECOCanvas.Write();



  MatchedLooseIsoRECOPtCanvas.Write();
  MatchedMedIsoRECOPtCanvas.Write();
  MatchedTightIsoRECOPtCanvas.Write();
  MatchedDMFindRECOPtCanvas.Write();
  MatchedRECOPtCanvas.Write();
  FinalEffLooseIsoRECOPtCanvas.Write();
  FinalEffMedIsoRECOPtCanvas.Write();
  FinalEffTightIsoRECOPtCanvas.Write();
  FinalEffDMFindRECOPtCanvas.Write();
  FinalEffLooseIsoRECOPtOverDMCanvas.Write();
  FinalEffMedIsoRECOPtOverDMCanvas.Write();
  FinalEffTightIsoRECOPtOverDMCanvas.Write();
  
  MatchedLooseIsoRECOPtMVACanvas.Write();
  MatchedMedIsoRECOPtMVACanvas.Write();
  MatchedTightIsoRECOPtMVACanvas.Write();
  MatchedRECOPtMVACanvas.Write();
  FinalEffLooseIsoRECOPtMVACanvas.Write();
  FinalEffMedIsoRECOPtMVACanvas.Write();
  FinalEffTightIsoRECOPtMVACanvas.Write();



  MatchedLooseIsoRECOdRCanvas.Write();
  MatchedMedIsoRECOdRCanvas.Write();
  MatchedTightIsoRECOdRCanvas.Write();
  MatchedDMFindRECOdRCanvas.Write();
  MatchedRECOdRCanvas.Write();
  FinalEffLooseIsoRECOdRCanvas.Write();
  FinalEffMedIsoRECOdRCanvas.Write();
  FinalEffTightIsoRECOdRCanvas.Write();
  FinalEffDMFindRECOdRCanvas.Write();
  FinalEffLooseIsoRECOdROverDMCanvas.Write();
  FinalEffMedIsoRECOdROverDMCanvas.Write();
  FinalEffTightIsoRECOdROverDMCanvas.Write();
  
  MatchedLooseIsoRECOdRMVACanvas.Write();
  MatchedMedIsoRECOdRMVACanvas.Write();
  MatchedTightIsoRECOdRMVACanvas.Write();
  MatchedRECOdRMVACanvas.Write();
  FinalEffLooseIsoRECOdRMVACanvas.Write();
  FinalEffMedIsoRECOdRMVACanvas.Write();
  FinalEffTightIsoRECOdRMVACanvas.Write();



  MatchedLooseIsoRECOPtGenCanvas.Write();
  MatchedMedIsoRECOPtGenCanvas.Write();
  MatchedTightIsoRECOPtGenCanvas.Write();
  MatchedDMFindRECOPtGenCanvas.Write();
  MatchedRECOPtGenCanvas.Write();
  FinalEffLooseIsoRECOPtGenCanvas.Write();
  FinalEffMedIsoRECOPtGenCanvas.Write();
  FinalEffTightIsoRECOPtGenCanvas.Write();
  FinalEffDMFindRECOPtGenCanvas.Write();
  FinalEffLooseIsoRECOPtGenOverDMCanvas.Write();
  FinalEffMedIsoRECOPtGenOverDMCanvas.Write();
  FinalEffTightIsoRECOPtGenOverDMCanvas.Write();
  
  MatchedLooseIsoRECOPtGenMVACanvas.Write();
  MatchedMedIsoRECOPtGenMVACanvas.Write();
  MatchedTightIsoRECOPtGenMVACanvas.Write();
  MatchedRECOPtGenMVACanvas.Write();
  FinalEffLooseIsoRECOPtGenMVACanvas.Write();
  FinalEffMedIsoRECOPtGenMVACanvas.Write();
  FinalEffTightIsoRECOPtGenMVACanvas.Write();
  


  OneProngDMRECOPtCanvas.Write();
  OneProngOnePizDMRECOPtCanvas.Write();
  OneProngTwoPizDMRECOPtCanvas.Write();
  ThreeProngDMRECOPtCanvas.Write();
  MatchedOneProngRECOPtCanvas.Write();
  MatchedOneProngOnePizRECOPtCanvas.Write();
  MatchedOneProngTwoPizRECOPtCanvas.Write();
  MatchedThreeProngRECOPtCanvas.Write();
  FinalOneProngDMRECOPtCanvas.Write();
  FinalOneProngOnePizDMRECOPtCanvas.Write();
  FinalOneProngTwoPizDMRECOPtCanvas.Write();
  FinalThreeProngDMRECOPtCanvas.Write();

  out_->Write();
  out_->Close();
std::cout << "DONE" << std::endl;
}//EndJob

// ------------ method called when starting to processes a run  ------------
void GGHAnalyzer_IndivRECO::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void GGHAnalyzer_IndivRECO::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void GGHAnalyzer_IndivRECO::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void GGHAnalyzer_IndivRECO::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

//Delete Memory
void GGHAnalyzer_IndivRECO::reset(const bool doDelete)
{
  if ((doDelete) && (NEvents_ != NULL)) delete NEvents_;
  NEvents_ = NULL;
  if ((doDelete) && (MuRemovedJetPt_ != NULL)) delete MuRemovedJetPt_;
  MuRemovedJetPt_ = NULL;
  if ((doDelete) && (TauMuTauHaddR_ != NULL)) delete TauMuTauHaddR_;
  TauMuTauHaddR_ = NULL;
  if ((doDelete) && (NConstituentsRECO_ != NULL)) delete NConstituentsRECO_;
  NConstituentsRECO_ = NULL;
  if ((doDelete) && (NTauDecayModeGEN_ != NULL)) delete NTauDecayModeGEN_;
  NTauDecayModeGEN_ = NULL;
  if ((doDelete) && (NTauDecayModeRECO_ != NULL)) delete NTauDecayModeRECO_;
  NTauDecayModeRECO_ = NULL;
  if ((doDelete) && (GenMatchedIsoPtSumRECO_ != NULL)) delete GenMatchedIsoPtSumRECO_;
  GenMatchedIsoPtSumRECO_ = NULL;
  if ((doDelete) && (TotalMuvsSoftMu_ != NULL)) delete TotalMuvsSoftMu_;
  TotalMuvsSoftMu_ = NULL;
  if ((doDelete) && (TauHadnConstvsPt_ != NULL)) delete TauHadnConstvsPt_;
  TauHadnConstvsPt_ = NULL;
  if ((doDelete) && (RECOPtvsGENVisiblePtRECO_ != NULL)) delete RECOPtvsGENVisiblePtRECO_;
  RECOPtvsGENVisiblePtRECO_ = NULL;
  if ((doDelete) && (RECODecayModevsGENDecayModeRECO_ != NULL)) delete RECODecayModevsGENDecayModeRECO_;
  RECODecayModevsGENDecayModeRECO_ = NULL;


  if ((doDelete) && (MatchedLooseIsoRECOPt_ != NULL)) delete MatchedLooseIsoRECOPt_;
  MatchedLooseIsoRECOPt_ = NULL;
  if ((doDelete) && (MatchedMedIsoRECOPt_ != NULL)) delete MatchedMedIsoRECOPt_;
  MatchedMedIsoRECOPt_ = NULL;
  if ((doDelete) && (MatchedTightIsoRECOPt_ != NULL)) delete MatchedTightIsoRECOPt_;
  MatchedTightIsoRECOPt_ = NULL;
  if ((doDelete) && (MatchedDMFindRECOPt_ != NULL)) delete MatchedDMFindRECOPt_;
  MatchedDMFindRECOPt_ = NULL;
  if ((doDelete) && (MatchedRECOPt_ != NULL)) delete MatchedRECOPt_;
  MatchedRECOPt_ = NULL;
  if ((doDelete) && (FinalEffLooseIsoRECOPt_ != NULL)) delete FinalEffLooseIsoRECOPt_;
  FinalEffLooseIsoRECOPt_ = NULL;
  if ((doDelete) && (FinalEffMedIsoRECOPt_ != NULL)) delete FinalEffMedIsoRECOPt_;
  FinalEffMedIsoRECOPt_ = NULL;
  if ((doDelete) && (FinalEffTightIsoRECOPt_ != NULL)) delete FinalEffTightIsoRECOPt_;
  FinalEffTightIsoRECOPt_ = NULL;
  if ((doDelete) && (FinalEffDMFindRECOPt_ != NULL)) delete FinalEffDMFindRECOPt_;
  FinalEffDMFindRECOPt_ = NULL;

  if ((doDelete) && (MatchedLooseIsoRECOdR_ != NULL)) delete MatchedLooseIsoRECOdR_;
  MatchedLooseIsoRECOdR_ = NULL;
  if ((doDelete) && (MatchedMedIsoRECOdR_ != NULL)) delete MatchedMedIsoRECOdR_;
  MatchedMedIsoRECOdR_ = NULL;
  if ((doDelete) && (MatchedTightIsoRECOdR_ != NULL)) delete MatchedTightIsoRECOdR_;
  MatchedTightIsoRECOdR_ = NULL;
  if ((doDelete) && (MatchedDMFindRECOdR_ != NULL)) delete MatchedDMFindRECOdR_;
  MatchedDMFindRECOdR_ = NULL;
  if ((doDelete) && (MatchedRECOdR_ != NULL)) delete MatchedRECOdR_;
  MatchedRECOdR_ = NULL;
  if ((doDelete) && (FinalEffLooseIsoRECOdR_ != NULL)) delete FinalEffLooseIsoRECOdR_;
  FinalEffLooseIsoRECOdR_ = NULL;
  if ((doDelete) && (FinalEffMedIsoRECOdR_ != NULL)) delete FinalEffMedIsoRECOdR_;
  FinalEffMedIsoRECOdR_ = NULL;
  if ((doDelete) && (FinalEffTightIsoRECOdR_ != NULL)) delete FinalEffTightIsoRECOdR_;
  FinalEffTightIsoRECOdR_ = NULL;
  if ((doDelete) && (FinalEffDMFindRECOdR_ != NULL)) delete FinalEffDMFindRECOdR_;
  FinalEffDMFindRECOdR_ = NULL;

  if ((doDelete) && (OneProngDMRECOPt_ != NULL)) delete OneProngDMRECOPt_;
  OneProngDMRECOPt_ = NULL;
  if ((doDelete) && (OneProngOnePizDMRECOPt_ != NULL)) delete OneProngOnePizDMRECOPt_;
  OneProngOnePizDMRECOPt_ = NULL;
  if ((doDelete) && (OneProngTwoPizDMRECOPt_ != NULL)) delete OneProngTwoPizDMRECOPt_;
  OneProngOnePizDMRECOPt_ = NULL;
  if ((doDelete) && (ThreeProngDMRECOPt_ != NULL)) delete ThreeProngDMRECOPt_;
  ThreeProngDMRECOPt_ = NULL;
  if ((doDelete) && (MatchedOneProngRECOPt_ != NULL)) delete MatchedOneProngRECOPt_;
  MatchedOneProngRECOPt_ = NULL;
  if ((doDelete) && (MatchedOneProngOnePizRECOPt_ != NULL)) delete MatchedOneProngOnePizRECOPt_;
  MatchedOneProngOnePizRECOPt_ = NULL;
  if ((doDelete) && (MatchedOneProngTwoPizRECOPt_ != NULL)) delete MatchedOneProngTwoPizRECOPt_;
  MatchedOneProngTwoPizRECOPt_ = NULL;
  if ((doDelete) && (MatchedThreeProngRECOPt_ != NULL)) delete MatchedThreeProngRECOPt_;
  MatchedThreeProngRECOPt_ = NULL;
  if ((doDelete) && (FinalOneProngDMRECOPt_ != NULL)) delete FinalOneProngDMRECOPt_;
  FinalOneProngDMRECOPt_ = NULL;
  if ((doDelete) && (FinalOneProngOnePizDMRECOPt_ != NULL)) delete FinalOneProngOnePizDMRECOPt_;
  FinalOneProngOnePizDMRECOPt_ = NULL;
  if ((doDelete) && (FinalOneProngTwoPizDMRECOPt_ != NULL)) delete FinalOneProngTwoPizDMRECOPt_;
  FinalOneProngTwoPizDMRECOPt_ = NULL;
  if ((doDelete) && (FinalThreeProngDMRECOPt_ != NULL)) delete FinalThreeProngDMRECOPt_;
  FinalThreeProngDMRECOPt_ = NULL;


  if ((doDelete) && (MatchedLooseIsoRECOPtGen_ != NULL)) delete MatchedLooseIsoRECOPtGen_;
  MatchedLooseIsoRECOPtGen_ = NULL;
  if ((doDelete) && (MatchedMedIsoRECOPtGen_ != NULL)) delete MatchedMedIsoRECOPtGen_;
  MatchedMedIsoRECOPtGen_ = NULL;
  if ((doDelete) && (MatchedTightIsoRECOPtGen_ != NULL)) delete MatchedTightIsoRECOPtGen_;
  MatchedTightIsoRECOPtGen_ = NULL;
  if ((doDelete) && (MatchedDMFindRECOPtGen_ != NULL)) delete MatchedDMFindRECOPtGen_;
  MatchedDMFindRECOPtGen_ = NULL;
  if ((doDelete) && (MatchedRECOPtGen_ != NULL)) delete MatchedRECOPtGen_;
  MatchedRECOPtGen_ = NULL;
  if ((doDelete) && (FinalEffLooseIsoRECOPtGen_ != NULL)) delete FinalEffLooseIsoRECOPtGen_;
  FinalEffLooseIsoRECOPtGen_ = NULL;
  if ((doDelete) && (FinalEffMedIsoRECOPtGen_ != NULL)) delete FinalEffMedIsoRECOPtGen_;
  FinalEffMedIsoRECOPtGen_ = NULL;
  if ((doDelete) && (FinalEffTightIsoRECOPtGen_ != NULL)) delete FinalEffTightIsoRECOPtGen_;
  FinalEffTightIsoRECOPtGen_ = NULL;
  if ((doDelete) && (FinalEffDMFindRECOPtGen_ != NULL)) delete FinalEffDMFindRECOPtGen_;
  FinalEffDMFindRECOPtGen_ = NULL;


}//void GGHAnalyzer_IndivRECO

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GGHAnalyzer_IndivRECO::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GGHAnalyzer_IndivRECO);

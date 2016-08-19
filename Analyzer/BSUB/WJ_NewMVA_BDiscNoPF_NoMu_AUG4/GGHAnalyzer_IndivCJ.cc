// -*- C++ -*-
//
// Package:    GGHAnalyzer_IndivCJ
// Class:      GGHAnalyzer_IndivCJ
// 
/**\class GGHAnalyzer_IndivCJ GGHAnalyzer_IndivCJ.cc Analyzer/src/GGHAnalyzer_IndivCJ.cc

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

class GGHAnalyzer_IndivCJ : public edm::EDAnalyzer {
   public:
      typedef reco::JetFloatAssociation::Container JetBCEnergyRatioCollection;
      explicit GGHAnalyzer_IndivCJ(const edm::ParameterSet&);
      ~GGHAnalyzer_IndivCJ();

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
      edm::EDGetTokenT<vector<reco::PFTau> > tauCJTag_;
      edm::EDGetTokenT<vector<reco::RecoTauPiZero> > pizerosTag_;      
      edm::EDGetTokenT<reco::PFTauDiscriminator>  looseIsoTagCJ_;
      edm::EDGetTokenT<reco::PFTauDiscriminator>  medIsoTagCJ_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> tightIsoTagCJ_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> decayModeFindingTagCJ_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> isoPtSumTagCJ_;
      edm::EDGetTokenT<edm::ValueMap<double> > genMatchedTauVisiblePtMapTagCJ_;
      edm::EDGetTokenT<edm::ValueMap<int> > genMatchedTauDecayModeMapTagCJ_;
      edm::EDGetTokenT<edm::ValueMap<int> > genMatchedTauMatchedMapTagCJ_;
      edm::EDGetTokenT<vector<reco::PFTau> > genMatchedRecoTausCJ_;
      edm::EDGetTokenT<reco::PFTauDiscriminator>  looseIsoTagCJMVA_;
      edm::EDGetTokenT<reco::PFTauDiscriminator>  medIsoTagCJMVA_;
      edm::EDGetTokenT<reco::PFTauDiscriminator> tightIsoTagCJMVA_;
      edm::EDGetTokenT<reco::PFJetCollection>  oldJetTag_;
      edm::EDGetTokenT<reco::JetTagCollection> csvBTag_;


      //Histograms
      TH1F* NEvents_;
      TH1F* MuRemovedJetPt_;
      TH1F* TauMuTauHaddR_;
      TH1F* MatchedTauMuPt_;
      TH1F* MatchedTauHadPt_;
      TH1F* MatchedTauHadEta_;
      TH1F* MatchedHT_;
      TH1F* MatchedBDiscCSV_;
      TH1F* NEventsCuts_;
      TH1F* NConstituentsCJ_;
      TH1F* NTauDecayModeGEN_;
      TH1F* NTauDecayModeRECO_;
      TH1F* DRTauHadRemovedMu_;
      TH1F* GenDiTaudRIfRemovedMu_;
      TH1F* GenMatchedIsoPtSumCJ_;
      TH2F* TotalMuvsSoftMu_;
      TH2F* DRGenTausvsRECOTauRemovedMu_;
      TH2F* GenDiTaudRvsCJDiTaudR_;
      TH2F* TauHadnConstvsPt_;
      TH2F* RECOPtvsGENVisiblePtCJ_;
      TH2F* RECODecayModevsGENDecayModeCJ_;

      TH1F* MatchedLooseIsoCJPt_;
      TH1F* MatchedMedIsoCJPt_;
      TH1F* MatchedTightIsoCJPt_;
      TH1F* MatchedDMFindCJPt_;
      TH1F* MatchedCJPt_;
      TGraphAsymmErrors* FinalEffLooseIsoCJPt_;
      TGraphAsymmErrors* FinalEffMedIsoCJPt_;
      TGraphAsymmErrors* FinalEffTightIsoCJPt_;
      TGraphAsymmErrors* FinalEffDMFindCJPt_;
      TGraphAsymmErrors* FinalEffLooseIsoCJPtOverDM_;
      TGraphAsymmErrors* FinalEffMedIsoCJPtOverDM_;
      TGraphAsymmErrors* FinalEffTightIsoCJPtOverDM_;

      TH1F* MatchedLooseIsoCJPtMVA_;
      TH1F* MatchedMedIsoCJPtMVA_;
      TH1F* MatchedTightIsoCJPtMVA_;
      TH1F* MatchedCJPtMVA_;
      TGraphAsymmErrors* FinalEffLooseIsoCJPtMVA_;
      TGraphAsymmErrors* FinalEffMedIsoCJPtMVA_;
      TGraphAsymmErrors* FinalEffTightIsoCJPtMVA_;



      TH1F* MatchedLooseIsoCJdR_;
      TH1F* MatchedMedIsoCJdR_;
      TH1F* MatchedTightIsoCJdR_;
      TH1F* MatchedDMFindCJdR_;
      TH1F* MatchedCJdR_;
      TGraphAsymmErrors* FinalEffLooseIsoCJdR_;
      TGraphAsymmErrors* FinalEffMedIsoCJdR_;
      TGraphAsymmErrors* FinalEffTightIsoCJdR_;
      TGraphAsymmErrors* FinalEffDMFindCJdR_;
      TGraphAsymmErrors* FinalEffLooseIsoCJdROverDM_;
      TGraphAsymmErrors* FinalEffMedIsoCJdROverDM_;
      TGraphAsymmErrors* FinalEffTightIsoCJdROverDM_;
      
      TH1F* MatchedLooseIsoCJdRMVA_;
      TH1F* MatchedMedIsoCJdRMVA_;
      TH1F* MatchedTightIsoCJdRMVA_;
      TH1F* MatchedCJdRMVA_;
      TGraphAsymmErrors* FinalEffLooseIsoCJdRMVA_;
      TGraphAsymmErrors* FinalEffMedIsoCJdRMVA_;
      TGraphAsymmErrors* FinalEffTightIsoCJdRMVA_;



      TH1F* MatchedLooseIsoCJPtGen_;
      TH1F* MatchedMedIsoCJPtGen_;
      TH1F* MatchedTightIsoCJPtGen_;
      TH1F* MatchedDMFindCJPtGen_;
      TH1F* MatchedCJPtGen_;
      TGraphAsymmErrors* FinalEffLooseIsoCJPtGen_;
      TGraphAsymmErrors* FinalEffMedIsoCJPtGen_;
      TGraphAsymmErrors* FinalEffTightIsoCJPtGen_;
      TGraphAsymmErrors* FinalEffDMFindCJPtGen_;
      TGraphAsymmErrors* FinalEffLooseIsoCJPtGenOverDM_;
      TGraphAsymmErrors* FinalEffMedIsoCJPtGenOverDM_;
      TGraphAsymmErrors* FinalEffTightIsoCJPtGenOverDM_;

      TH1F* MatchedLooseIsoCJPtGenMVA_;
      TH1F* MatchedMedIsoCJPtGenMVA_;
      TH1F* MatchedTightIsoCJPtGenMVA_;
      TH1F* MatchedCJPtGenMVA_;
      TGraphAsymmErrors* FinalEffLooseIsoCJPtGenMVA_;
      TGraphAsymmErrors* FinalEffMedIsoCJPtGenMVA_;
      TGraphAsymmErrors* FinalEffTightIsoCJPtGenMVA_;



      TH1F* OneProngDMCJPt_;
      TH1F* OneProngOnePizDMCJPt_;
      TH1F* OneProngTwoPizDMCJPt_;
      TH1F* ThreeProngDMCJPt_;
      TH1F* MatchedOneProngCJPt_;
      TH1F* MatchedOneProngOnePizCJPt_;
      TH1F* MatchedOneProngTwoPizCJPt_;
      TH1F* MatchedThreeProngCJPt_;
      TGraphAsymmErrors* FinalOneProngDMCJPt_;
      TGraphAsymmErrors* FinalOneProngOnePizDMCJPt_;
      TGraphAsymmErrors* FinalOneProngTwoPizDMCJPt_;
      TGraphAsymmErrors* FinalThreeProngDMCJPt_;

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
GGHAnalyzer_IndivCJ::GGHAnalyzer_IndivCJ(const edm::ParameterSet& iConfig):
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  genParticleTag_(consumes<vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleTag"))),
  akJetTag_(consumes<vector<reco::PFJet> >(iConfig.getParameter<edm::InputTag>("akJetTag"))),
  muonsTag_(consumes<vector<reco::Muon> >(iConfig.getParameter<edm::InputTag>("muonsTag"))),
  muonMapTag_(consumes<edm::ValueMap<edm::RefVector<vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<vector<reco::Muon>,reco::Muon> > > >(iConfig.getParameter<edm::InputTag>("muonMapTag"))),
  jetValMapTag_(consumes<edm::ValueMap<reco::PFJetRef> >(iConfig.getParameter<edm::InputTag>("jetValMapTag"))),
  tauCJTag_(consumes<vector<reco::PFTau> >(iConfig.getParameter<edm::InputTag>("tauCJTag"))),
  pizerosTag_(consumes<vector<reco::RecoTauPiZero> >(iConfig.getParameter<edm::InputTag>("pizerosTag"))),
  looseIsoTagCJ_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("looseIsoTagCJ"))),
  medIsoTagCJ_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("medIsoTagCJ"))),
  tightIsoTagCJ_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("tightIsoTagCJ"))),
  decayModeFindingTagCJ_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("decayModeFindingTagCJ"))),
  isoPtSumTagCJ_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("isoPtSumTagCJ"))),
  genMatchedTauVisiblePtMapTagCJ_(consumes<edm::ValueMap<double> >(iConfig.getParameter<edm::InputTag>("genMatchedTauVisiblePtMapTagCJ"))),
  genMatchedTauDecayModeMapTagCJ_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("genMatchedTauDecayModeMapTagCJ"))),
  genMatchedTauMatchedMapTagCJ_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("genMatchedTauMatchedMapTagCJ"))),
  genMatchedRecoTausCJ_(consumes<vector<reco::PFTau> >(iConfig.getParameter<edm::InputTag>("genMatchedRecoTausCJ"))),
  looseIsoTagCJMVA_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("looseIsoTagCJMVA"))),
  medIsoTagCJMVA_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("medIsoTagCJMVA"))),
  tightIsoTagCJMVA_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("tightIsoTagCJMVA"))),
  oldJetTag_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("oldJetTag"))),
  csvBTag_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("csvBTag")))
{
  reset(false);    
}//GGHAnalyzer_IndivCJ



GGHAnalyzer_IndivCJ::~GGHAnalyzer_IndivCJ()
{
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void GGHAnalyzer_IndivCJ::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  edm::Handle<std::vector<reco::PFTau> > pTausCJ;
  iEvent.getByToken(tauCJTag_, pTausCJ);

  //Get pi-Zero particle collection
  edm::Handle<std::vector<reco::RecoTauPiZero> > pPiZero;
  iEvent.getByToken(pizerosTag_, pPiZero);

  //Get Loose Iso Collection
  Handle<PFTauDiscriminator> pLooseIsoDiscCJ; 
  iEvent.getByToken(looseIsoTagCJ_, pLooseIsoDiscCJ); 

  //Get Medium Iso Collection
  Handle<PFTauDiscriminator> pMedIsoDiscCJ; 
  iEvent.getByToken(medIsoTagCJ_, pMedIsoDiscCJ);

  //Get Tight Iso Collection
  Handle<PFTauDiscriminator> pTightIsoDiscCJ; 
  iEvent.getByToken(tightIsoTagCJ_, pTightIsoDiscCJ);

  //Get Decay Mode Finding Collection
  Handle<PFTauDiscriminator> pDMFindingCJ; 
  iEvent.getByToken(decayModeFindingTagCJ_, pDMFindingCJ);

  //Get Loose Iso Collection
  Handle<PFTauDiscriminator> pLooseIsoDiscCJMVA;
  iEvent.getByToken(looseIsoTagCJMVA_, pLooseIsoDiscCJMVA);

  //Get Medium Iso Collection
  Handle<PFTauDiscriminator> pMedIsoDiscCJMVA;
  iEvent.getByToken(medIsoTagCJMVA_, pMedIsoDiscCJMVA);

  //Get Tight Iso Collection
  Handle<PFTauDiscriminator> pTightIsoDiscCJMVA;
  iEvent.getByToken(tightIsoTagCJMVA_, pTightIsoDiscCJMVA);

  //Get Isolation Pt Sum Finding Collection
  Handle<PFTauDiscriminator> pIsoPtSumCJ;
  iEvent.getByToken(isoPtSumTagCJ_, pIsoPtSumCJ);

  //Gen Matched Visible Pt Map
  edm::Handle<edm::ValueMap<double> > pGenMatchedTauVisiblePtMapCJ;
  iEvent.getByToken(genMatchedTauVisiblePtMapTagCJ_, pGenMatchedTauVisiblePtMapCJ);

  //Gen Matched Decay Mode Map
  edm::Handle<edm::ValueMap<int> > pGenMatchedTauDecayModeMapCJ;
  iEvent.getByToken(genMatchedTauDecayModeMapTagCJ_, pGenMatchedTauDecayModeMapCJ);

  //Gen Matched Truth Map
  edm::Handle<edm::ValueMap<int> > pGenMatchedTauMatchedMapCJ;
  iEvent.getByToken(genMatchedTauMatchedMapTagCJ_, pGenMatchedTauMatchedMapCJ);

  //Get RECO Taus particle collection
  edm::Handle<std::vector<reco::PFTau> > pAccessersCJ;
  iEvent.getByToken(genMatchedRecoTausCJ_, pAccessersCJ);

  //Old Jet collection for bTagging
  edm::Handle<reco::PFJetCollection> pOldJets;
  iEvent.getByToken(oldJetTag_, pOldJets);

  //Get combVertMVA JetTagCollection
  edm::Handle<reco::JetTagCollection> pCSV;
  iEvent.getByToken(csvBTag_, pCSV);


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

/*      unsigned int iDaughters1 = 0, iDaughters2 = 0;
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
      }//else 
    }//ir found a1 and it decayed to taus
  }//for iGen Particle

//////////////////////////////
// Looking at CleanJets
//////////////////////////////
  bool muonWasRemoved = false;
  for(size_t iTauCJ = 0; iTauCJ < pAccessersCJ->size(); ++iTauCJ)
  {
    const reco::PFTauRef accesserTauRef(pAccessersCJ, iTauCJ);//
    const reco::PFTauRef hpsTauRef(pTausCJ, iTauCJ);//

    const int    GenTauMatchedMapCJ   = (*pGenMatchedTauMatchedMapCJ)[accesserTauRef];
    const int    GenTauDecayModeMapCJ = (*pGenMatchedTauDecayModeMapCJ)[accesserTauRef];
    const double GenTauVisiblePtMapCJ = (*pGenMatchedTauVisiblePtMapCJ)[accesserTauRef];
    const reco::PFJetRef& tauJetRef = (*accesserTauRef).jetRef();
    const reco::MuonRefVector& removedMuons = (*pMuonMap)[tauJetRef];
    std::vector<reco::PFCandidatePtr> JetPFCands = tauJetRef->getPFConstituents();
    NConstituentsCJ_->Fill(JetPFCands.size() );

    /////////////////////////////////////////
    // find the highest pT associated muon
    /////////////////////////////////////////
    std::vector<reco::MuonRef> removedMuonRefs;
    for (reco::MuonRefVector::const_iterator iMuonRef = removedMuons.begin(); iMuonRef != removedMuons.end(); ++iMuonRef)
      removedMuonRefs.push_back(*iMuonRef);

    if (removedMuonRefs.size() > 0)
    {
      MuRemovedJetPt_->Fill(tauJetRef->pt() );
      muonWasRemoved = true;
    }//if removedMuonRefs,.size() > 0
    TotalMuvsSoftMu_->Fill(removedMuonRefs.size(), removedMuonRefs.size() );
    for (unsigned int iter = 0; iter < removedMuonRefs.size(); iter++)
    {
      double dPhi = reco::deltaPhi(accesserTauRef->phi(), removedMuonRefs[iter]->phi() ), dEta = accesserTauRef->eta() - removedMuonRefs[iter]->eta();
      double dRTauRemovedMu = sqrt( dEta * dEta + dPhi * dPhi);
      std::cout << "\tdRTauRemovedMu= " << dRTauRemovedMu << " \tGEN_diTaudR= " << GEN_diTaudR << " \tremovedMuonRefs.size()= " << removedMuonRefs.size() <<  std::endl;
      DRTauHadRemovedMu_->Fill(dRTauRemovedMu );
      GenDiTaudRIfRemovedMu_->Fill(GEN_diTaudR );
      MuRemovedJetPt_->Fill(tauJetRef->pt() );
      GenDiTaudRvsCJDiTaudR_->Fill(GEN_diTaudR, dRTauRemovedMu);
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

    //////////////////////////////////
    // Making sure a gen match exists
    //////////////////////////////////
    if (GenTauMatchedMapCJ == 1 && fabs(hpsTauRef->eta() ) < 2.3 && hpsTauRef->pt() > 20) //Checking that it's gen matched
    {
      int DM = (*pDMFindingCJ)[hpsTauRef], MedIso = (*pMedIsoDiscCJ)[hpsTauRef], LooseIso = (*pLooseIsoDiscCJ)[hpsTauRef], TightIso = (*pTightIsoDiscCJ)[hpsTauRef];
      double PtSum = (*pIsoPtSumCJ)[hpsTauRef];
      GenMatchedIsoPtSumCJ_->Fill(PtSum );
      int MedIsoMVA = (*pMedIsoDiscCJMVA)[hpsTauRef], LooseIsoMVA = (*pLooseIsoDiscCJMVA)[hpsTauRef], TightIsoMVA = (*pTightIsoDiscCJMVA)[hpsTauRef];


      //////////////////////////////
      // Cut Option information
      //////////////////////////////
      MatchedTauHadEta_->Fill(hpsTauRef->eta() );
      MatchedTauHadPt_->Fill(hpsTauRef->pt() );
      if (removedMuonRefs.size() > 0)
        MatchedTauMuPt_->Fill(removedMuonRefs[0]->pt() );
      NEventsCuts_->Fill(0);
      double sumHT = 0;
      for (reco::PFJetCollection::const_iterator iAkJet = pAkJets->begin(); iAkJet != pAkJets->end(); ++iAkJet)
      {
        if (iAkJet->pt() > 20.0)
          sumHT += iAkJet->pt();
      }//for pAkJets      
      MatchedHT_->Fill(sumHT);

      double BDisc = -1;
      for (unsigned int iBTagInfo = 0; iBTagInfo != pCSV->size(); ++iBTagInfo)
      {
        if (reco::PFJetRef(pOldJets, iBTagInfo).key() == tauJetRef.key())
        {
          const reco::JetTagCollection& bTags = *(pCSV.product());
          BDisc = bTags[iBTagInfo].second;
          std::cout << "\tBDisc= " << BDisc << std::endl;
          MatchedBDiscCSV_->Fill(BDisc );
        }//if
      }//for iBTagInfo

 
      // Filling nEvents for passing specific cuts
      if (hpsTauRef->pt() < 200)
        NEventsCuts_->Fill(1);
      if (removedMuonRefs[0]->pt() < 100)
        NEventsCuts_->Fill(2);
      if (hpsTauRef->eta() < 1.75 );
        NEventsCuts_->Fill(3);
      if (sumHT < 300)
        NEventsCuts_->Fill(4);
      if (BDisc > .5 )
        NEventsCuts_->Fill(5);

      //////////////////////////////
      // Filling Disc Efficiency 
      //////////////////////////////
      NEvents_->Fill(3);
      MatchedCJPt_->Fill(hpsTauRef->pt() );
      RECOPtvsGENVisiblePtCJ_->Fill(hpsTauRef->pt(), GenTauVisiblePtMapCJ );
      if (DM == 1)
        MatchedDMFindCJPt_->Fill(hpsTauRef->pt() );
      if (TightIso == 1 && DM == 1)
        MatchedTightIsoCJPt_->Fill(hpsTauRef->pt() );
      if (MedIso == 1 && DM == 1)
        MatchedMedIsoCJPt_->Fill(hpsTauRef->pt() );
      if (LooseIso == 1 && DM == 1)
        MatchedLooseIsoCJPt_->Fill(hpsTauRef->pt() );
  
      MatchedCJPtMVA_->Fill(hpsTauRef->pt() );
      if (TightIsoMVA == 1)
        MatchedTightIsoCJPtMVA_->Fill(hpsTauRef->pt() );
      if (MedIsoMVA == 1)
        MatchedMedIsoCJPtMVA_->Fill(hpsTauRef->pt() );
      if (LooseIsoMVA == 1)
        MatchedLooseIsoCJPtMVA_->Fill(hpsTauRef->pt() );
  
  
    
      MatchedCJdR_->Fill(GEN_diTaudR );
      if (DM == 1)
        MatchedDMFindCJdR_->Fill(GEN_diTaudR );
      if (TightIso == 1 && DM == 1)
        MatchedTightIsoCJdR_->Fill(GEN_diTaudR );
      if (MedIso == 1 && DM == 1)
        MatchedMedIsoCJdR_->Fill(GEN_diTaudR );
      if (LooseIso == 1 && DM == 1)
        MatchedLooseIsoCJdR_->Fill(GEN_diTaudR );
      
      MatchedCJdRMVA_->Fill(GEN_diTaudR );
      if (TightIsoMVA == 1)
        MatchedTightIsoCJdRMVA_->Fill(GEN_diTaudR );
      if (MedIsoMVA == 1)
        MatchedMedIsoCJdRMVA_->Fill(GEN_diTaudR );
      if (LooseIsoMVA == 1)
        MatchedLooseIsoCJdRMVA_->Fill(GEN_diTaudR );
      
  
  
      MatchedCJPtGen_->Fill(GenTauVisiblePtMapCJ );
      if (DM == 1)
        MatchedDMFindCJPtGen_->Fill(GenTauVisiblePtMapCJ );
      if (TightIso == 1 && DM == 1)
        MatchedTightIsoCJPtGen_->Fill(GenTauVisiblePtMapCJ );
      if (MedIso == 1 && DM == 1)
        MatchedMedIsoCJPtGen_->Fill(GenTauVisiblePtMapCJ );
      if (LooseIso == 1 && DM == 1)
        MatchedLooseIsoCJPtGen_->Fill(GenTauVisiblePtMapCJ );
  
      MatchedCJPtGenMVA_->Fill(GenTauVisiblePtMapCJ );
      if (TightIsoMVA == 1)
        MatchedTightIsoCJPtGenMVA_->Fill(GenTauVisiblePtMapCJ );
      if (MedIsoMVA == 1)
        MatchedMedIsoCJPtGenMVA_->Fill(GenTauVisiblePtMapCJ );
      if (LooseIsoMVA == 1)
        MatchedLooseIsoCJPtGenMVA_->Fill(GenTauVisiblePtMapCJ );
 
/*      int RECODecayMode = -1;
      std::cout << "\tpTausCJ->at(iTauCJ).decayMode()= " << pTausCJ->at(iTauCJ).decayMode() << "\tGenTauDecayModeMapCJ= " << GenTauDecayModeMapCJ  << std::endl;
      std::cout << "\tpTausCJ->at(iTauCJ).signalPFChargedHadrCands().size()= " << pTausCJ->at(iTauCJ).signalPFChargedHadrCands().size() << "\tpTausCJ->at(iTauCJ).signalPFGammaCands().size()= " << pTausCJ->at(iTauCJ).signalPFGammaCands().size() << std::endl;
      if (pTausCJ->at(iTauCJ).decayMode() == 0)
        RECODecayMode = 1;
      if (pTausCJ->at(iTauCJ).decayMode() == 1)
        RECODecayMode = 2;
      if (pTausCJ->at(iTauCJ).decayMode() == 2)
        RECODecayMode = 3;
      if (pTausCJ->at(iTauCJ).decayMode() == 10)
        RECODecayMode = 4;
*/

      NTauDecayModeRECO_->Fill(pTausCJ->at(iTauCJ).decayMode() );
 
  
      if (GenTauDecayModeMapCJ == 1)
      {
        RECODecayModevsGENDecayModeCJ_->Fill(pTausCJ->at(iTauCJ).decayMode(), GenTauDecayModeMapCJ);
        NTauDecayModeGEN_->Fill(1);
        MatchedOneProngCJPt_->Fill(GenTauVisiblePtMapCJ );
        if (DM == 1 && MedIso == 1)
  	OneProngDMCJPt_->Fill(GenTauVisiblePtMapCJ );     
      }//if decayMode == 1
      else if (GenTauDecayModeMapCJ == 2)
      {
        RECODecayModevsGENDecayModeCJ_->Fill(pTausCJ->at(iTauCJ).decayMode(), GenTauDecayModeMapCJ);
        NTauDecayModeGEN_->Fill(2);
        MatchedOneProngOnePizCJPt_->Fill(GenTauVisiblePtMapCJ );
        if (DM == 1 && MedIso == 1)
          OneProngOnePizDMCJPt_->Fill(GenTauVisiblePtMapCJ );
      }//if decayMode == 1
      else if (GenTauDecayModeMapCJ == 3)
      {
        RECODecayModevsGENDecayModeCJ_->Fill(pTausCJ->at(iTauCJ).decayMode(), GenTauDecayModeMapCJ);
        NTauDecayModeGEN_->Fill(3);
        MatchedOneProngTwoPizCJPt_->Fill(GenTauVisiblePtMapCJ );
        if (DM == 1 && MedIso == 1)
          OneProngTwoPizDMCJPt_->Fill(GenTauVisiblePtMapCJ );
      }//if decayMode == 1
      else if (GenTauDecayModeMapCJ == 4)
      {
        RECODecayModevsGENDecayModeCJ_->Fill(pTausCJ->at(iTauCJ).decayMode(), GenTauDecayModeMapCJ);
        NTauDecayModeGEN_->Fill(4);
        MatchedThreeProngCJPt_->Fill(GenTauVisiblePtMapCJ );
        if (DM == 1 && MedIso == 1)
          ThreeProngDMCJPt_->Fill(GenTauVisiblePtMapCJ );
      }//if decayMode == 1
    }//if GenTauMatchedMapCJ == 1 && fabs(hpsTauRef->eta() ) < 2.3
  }//if GEN tau1Ref is the had in mu+had and it is matched to CleanJets Je


  if (muonWasRemoved)
    NEvents_->Fill(5);

}//End GGHAnalyzer_IndivCJ::analyze


// ------------ method called once each job just before starting event loop  ------------
void GGHAnalyzer_IndivCJ::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  //Open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //Book histograms
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
  MatchedTauMuPt_       = new TH1F("MatchedTauMuPt"    , "", 200, 0, 500);
  MatchedTauHadPt_       = new TH1F("MatchedTauHadPt"    , "", 200, 0, 500);
  MatchedTauHadEta_       = new TH1F("MatchedTauHadEta"    , "", 100, -2.5, 2.5);
  MatchedHT_       = new TH1F("MatchedHT"    , "", 100, 0, 1000);
  MatchedBDiscCSV_       = new TH1F("MatchedBDiscCSV"    , "", 100, 0, 1);
  NEventsCuts_       = new TH1F("NEventsCuts"    , "", 6, -.5, 5.5);
      NEventsCuts_->GetXaxis()->SetBinLabel(1, "TotalEvents");
      NEventsCuts_->GetXaxis()->SetBinLabel(2, "p_{T}(#tau) < 200");
      NEventsCuts_->GetXaxis()->SetBinLabel(3, "p_{T}(#mu)  < 100");
      NEventsCuts_->GetXaxis()->SetBinLabel(4, "#eta < 1.75");
      NEventsCuts_->GetXaxis()->SetBinLabel(5, "HT < 300");
      NEventsCuts_->GetXaxis()->SetBinLabel(6, "BDisc < .5");
  NConstituentsCJ_        = new TH1F("NConstituentsCJ"    , "", 50, 0, 50);
  NTauDecayModeGEN_        = new TH1F("NTauDecayModeGEN"    , "", 8, -.5, 7.5);
  NTauDecayModeRECO_        = new TH1F("NTauDecayModeRECO"    , "", 11, -.5, 10.5);
  DRTauHadRemovedMu_        = new TH1F("DRTauHadRemovedMu"    , "", 100, 0, 8);
  GenDiTaudRIfRemovedMu_        = new TH1F("GenDiTaudRIfRemovedMu"    , "", 100, 0, 8);
  GenMatchedIsoPtSumCJ_        = new TH1F("GenMatchedIsoPtSumCJ"    , "", 100, 0, 20);
  TotalMuvsSoftMu_  = new TH2F("TotalMuvsSoftMu" , "", 21, -.5, 20.5, 21, -.5, 20.5);
  DRGenTausvsRECOTauRemovedMu_  = new TH2F("DRGenTausvsRECOTauRemovedMu" , "", 21, -.5, 20.5, 21, -.5, 20.5);
  GenDiTaudRvsCJDiTaudR_  = new TH2F("GenDiTaudRvsCJDiTaudR" , "", 50, 0, 10, 50, 0, 10);
  TauHadnConstvsPt_  = new TH2F("TauHadnConstvsPt" , "", 50, 0, 50, 30, 0, 220);
  RECOPtvsGENVisiblePtCJ_  = new TH2F("RECOPtvsGENVisiblePtCJ" , "", 50, 0, 220, 50, 0, 220);
  RECODecayModevsGENDecayModeCJ_  = new TH2F("RECODecayModevsGENDecayModeCJ" , "", 11, -0.5, 10.5, 4, 0.5, 4.5);



  MatchedLooseIsoCJPt_  = new TH1F("MatchedLooseIsoCJPt"    , "", 11, 0, 220.0);
  MatchedMedIsoCJPt_    = new TH1F("MatchedMedIsoCJPt", "", 11, 0, 220);
  MatchedTightIsoCJPt_    = new TH1F("MatchedTightIsoCJPt", "", 11, 0, 220);
  MatchedDMFindCJPt_    = new TH1F("MatchedDMFindCJPt"    , "", 11, 0, 220);
  MatchedCJPt_          = new TH1F("MatchedCJPt"    , "", 11, 0, 220);

  MatchedLooseIsoCJPtMVA_  = new TH1F("MatchedLooseIsoCJPtMVA"    , "", 11, 0, 220.0);
  MatchedMedIsoCJPtMVA_    = new TH1F("MatchedMedIsoCJPtMVA", "", 11, 0, 220);
  MatchedTightIsoCJPtMVA_    = new TH1F("MatchedTightIsoCJPtMVA", "", 11, 0, 220);
  MatchedCJPtMVA_          = new TH1F("MatchedCJPtMVA"    , "", 11, 0, 220);



  FinalEffLooseIsoCJPt_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoCJPt_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoCJPt_ = new TGraphAsymmErrors(30);
  FinalEffDMFindCJPt_ = new TGraphAsymmErrors(30);
  FinalEffLooseIsoCJPtOverDM_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoCJPtOverDM_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoCJPtOverDM_ = new TGraphAsymmErrors(30);

  FinalEffLooseIsoCJPtMVA_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoCJPtMVA_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoCJPtMVA_ = new TGraphAsymmErrors(30);



  MatchedLooseIsoCJdR_  = new TH1F("MatchedLooseIsoCJdR"    , "", 8, 0, 1.5);
  MatchedMedIsoCJdR_    = new TH1F("MatchedMedIsoCJdR", "", 8, 0, 1.5);
  MatchedTightIsoCJdR_    = new TH1F("MatchedTightIsoCJdR", "", 8, 0, 1.5);
  MatchedDMFindCJdR_    = new TH1F("MatchedDMFindCJdR"    , "", 8, 0, 1.5);
  MatchedCJdR_          = new TH1F("MatchedCJdR"    , "", 8, 0, 1.5);

  MatchedLooseIsoCJdRMVA_  = new TH1F("MatchedLooseIsoCJdRMVA"    , "", 8, 0, 1.5);
  MatchedMedIsoCJdRMVA_    = new TH1F("MatchedMedIsoCJdRMVA", "", 8, 0, 1.5);
  MatchedTightIsoCJdRMVA_    = new TH1F("MatchedTightIsoCJdRMVA", "", 8, 0, 1.5);
  MatchedCJdRMVA_          = new TH1F("MatchedCJdRMVA"    , "", 8, 0, 1.5);



  FinalEffLooseIsoCJdR_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoCJdR_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoCJdR_ = new TGraphAsymmErrors(30);
  FinalEffDMFindCJdR_ = new TGraphAsymmErrors(30);
  FinalEffLooseIsoCJdROverDM_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoCJdROverDM_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoCJdROverDM_ = new TGraphAsymmErrors(30);

  FinalEffLooseIsoCJdRMVA_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoCJdRMVA_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoCJdRMVA_ = new TGraphAsymmErrors(30);



  MatchedLooseIsoCJPtGen_  = new TH1F("MatchedLooseIsoCJPtGen"    , "", 11, 0, 220.0);
  MatchedMedIsoCJPtGen_    = new TH1F("MatchedMedIsoCJPtGen", "", 11, 0, 220);
  MatchedTightIsoCJPtGen_    = new TH1F("MatchedTightIsoCJPtGen", "", 11, 0, 220);
  MatchedDMFindCJPtGen_    = new TH1F("MatchedDMFindCJPtGen"    , "", 11, 0, 220);
  MatchedCJPtGen_          = new TH1F("MatchedCJPtGen"    , "", 11, 0, 220);
  
  MatchedLooseIsoCJPtGenMVA_  = new TH1F("MatchedLooseIsoCJPtGenMVA"    , "", 11, 0, 220.0);
  MatchedMedIsoCJPtGenMVA_    = new TH1F("MatchedMedIsoCJPtGenMVA", "", 11, 0, 220);
  MatchedTightIsoCJPtGenMVA_    = new TH1F("MatchedTightIsoCJPtGenMVA", "", 11, 0, 220);
  MatchedCJPtGenMVA_          = new TH1F("MatchedCJPtGenMVA"    , "", 11, 0, 220);
  


  FinalEffLooseIsoCJPtGen_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoCJPtGen_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoCJPtGen_ = new TGraphAsymmErrors(30);
  FinalEffDMFindCJPtGen_ = new TGraphAsymmErrors(30);
  FinalEffLooseIsoCJPtGenOverDM_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoCJPtGenOverDM_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoCJPtGenOverDM_ = new TGraphAsymmErrors(30);

  FinalEffLooseIsoCJPtGenMVA_ = new TGraphAsymmErrors(30);
  FinalEffMedIsoCJPtGenMVA_ = new TGraphAsymmErrors(30);
  FinalEffTightIsoCJPtGenMVA_ = new TGraphAsymmErrors(30);


  OneProngDMCJPt_ = new TH1F("OneProngDMCJPt"    , "", 11, 0, 220.0);
  OneProngOnePizDMCJPt_ = new TH1F("OneProngOnePizDMCJPt"    , "", 11, 0, 220.0);
  OneProngTwoPizDMCJPt_ = new TH1F("OneProngTwoPizDMCJPt"    , "", 11, 0, 220.0);
  ThreeProngDMCJPt_ = new TH1F("ThreeProngDMCJPt"    , "", 11, 0, 220.0);
  MatchedOneProngCJPt_ = new TH1F("MatchedOneProngCJPt", "", 11, 0, 220.0);
  MatchedOneProngOnePizCJPt_ = new TH1F("MatchedOneProngOnePizCJPt", "", 11, 0, 220.0);
  MatchedOneProngTwoPizCJPt_ = new TH1F("MatchedOneProngTwoPizCJPt", "", 11, 0, 220.0);
  MatchedThreeProngCJPt_ = new TH1F("MatchedThreeProngCJPt", "", 11, 0, 220.0);
  FinalOneProngDMCJPt_ = new TGraphAsymmErrors(30);
  FinalOneProngOnePizDMCJPt_ = new TGraphAsymmErrors(30);
  FinalOneProngTwoPizDMCJPt_ = new TGraphAsymmErrors(30);
  FinalThreeProngDMCJPt_ = new TGraphAsymmErrors(30);

}

// ------------ method called once each job just after ending the event loop  ------------
void GGHAnalyzer_IndivCJ::endJob()
{
  //Make the Canvases
  TCanvas NEventsCanvas("NEvents","",600,600);
  TCanvas MuRemovedJetPtCanvas("MuRemovedJetPt","",600,600);
  TCanvas TauMuTauHaddRCanvas("TauMuTauHaddR","",600,600);
  TCanvas MatchedTauMuPtCanvas("MatchedTauMuPt","",600,600);
  TCanvas MatchedTauHadPtCanvas("MatchedTauHadPt","",600,600);
  TCanvas MatchedTauHadEtaCanvas("MatchedTauHadEta","",600,600);
  TCanvas MatchedHTCanvas("MatchedHT","",600,600);
  TCanvas MatchedBDiscCSVCanvas("MatchedBDiscCSV","",600,600);
  TCanvas NEventsCutsCanvas("NEventsCuts","",600,600);
  TCanvas NConstituentsCJCanvas("NConstituentsCJ","",600,600);
  TCanvas NTauDecayModeGENCanvas("NTauDecayModeGEN","",600,600);
  TCanvas NTauDecayModeRECOCanvas("NTauDecayModeRECO","",600,600);
  TCanvas DRTauHadRemovedMuCanvas("DRTauHadRemovedMu","",600,600);
  TCanvas GenDiTaudRIfRemovedMuCanvas("GenDiTaudRIfRemovedMu","",600,600);
  TCanvas GenMatchedIsoPtSumCJCanvas("GenMatchedIsoPtSumCJ","",600,600);
  TCanvas TotalMuvsSoftMuCanvas("TotalMuvsSoftMu","",600,600);
  TCanvas DRGenTausvsRECOTauRemovedMuCanvas("DRGenTausvsRECOTauRemovedMu","",600,600);
  TCanvas GenDiTaudRvsCJDiTaudRCanvas("GenDiTaudRvsCJDiTaudR","",600,600);
  TCanvas TauHadnConstvsPtCanvas("TauHadnConstvsPt","",600,600);
  TCanvas RECOPtvsGENVisiblePtCJCanvas("RECOPtvsGENVisiblePtCJ","",600,600);
  TCanvas RECODecayModevsGENDecayModeCJCanvas("RECODecayModevsGENDecayModeCJ","",600,600);

/*  for (int x = 0; x < 4; x++)
  {
    int sum = 0;
    for (int y = 0; y < 4; y++)
      sum += RECODecayModevsGENDecayModeCJ_->GetBinContent(x,y);
    for (int y = 0; y < 4; y++)
    {
      std::cout << "Before(" << x << "," << y << ")= " << RECODecayModevsGENDecayModeCJ_->GetBinContent(x,y) << "\tsum= " << sum;
      std::cout << "\tRECODecayModevsGENDecayModeCJ_->GetBinContent(x,y) / sum = " <<  RECODecayModevsGENDecayModeCJ_->GetBinContent(x,y) / sum << std::endl;
      RECODecayModevsGENDecayModeCJ_->SetBinContent(x, y, RECODecayModevsGENDecayModeCJ_->GetBinContent(x,y) / sum );
      std::cout << "After (" << x << "," << y << ")= " << RECODecayModevsGENDecayModeCJ_->GetBinContent(x,y) << std::endl;
    }//for y
  }//for x
*/


  TCanvas MatchedLooseIsoCJPtCanvas("MatchedLooseIsoCJPt","",600,600);
  TCanvas MatchedMedIsoCJPtCanvas("MatchedMedIsoCJPt","",600,600);
  TCanvas MatchedTightIsoCJPtCanvas("MatchedTightIsoCJPt","",600,600);
  TCanvas MatchedDMFindCJPtCanvas("MatchedDMFindCJPt","",600,600);
  TCanvas MatchedCJPtCanvas("MatchedCJPt","",600,600);
  TCanvas FinalEffLooseIsoCJPtCanvas("FinalEffLooseIsoCJPt","",600,600);
  TCanvas FinalEffMedIsoCJPtCanvas("FinalEffMedIsoCJPt","",600,600);
  TCanvas FinalEffTightIsoCJPtCanvas("FinalEffTightIsoCJPt","",600,600);
  TCanvas FinalEffDMFindCJPtCanvas("FinalEffDMFindCJPt","",600,600);
  TCanvas FinalEffLooseIsoCJPtOverDMCanvas("FinalEffLooseIsoCJPtOverDM","",600,600);
  TCanvas FinalEffMedIsoCJPtOverDMCanvas("FinalEffMedIsoCJPtOverDM","",600,600);
  TCanvas FinalEffTightIsoCJPtOverDMCanvas("FinalEffTightIsoCJPtOverDM","",600,600);

  TCanvas MatchedLooseIsoCJPtMVACanvas("MatchedLooseIsoCJPtMVA","",600,600);
  TCanvas MatchedMedIsoCJPtMVACanvas("MatchedMedIsoCJPtMVA","",600,600);
  TCanvas MatchedTightIsoCJPtMVACanvas("MatchedTightIsoCJPtMVA","",600,600);
  TCanvas MatchedCJPtMVACanvas("MatchedCJPtMVA","",600,600);
  TCanvas FinalEffLooseIsoCJPtMVACanvas("FinalEffLooseIsoCJPtMVA","",600,600);
  TCanvas FinalEffMedIsoCJPtMVACanvas("FinalEffMedIsoCJPtMVA","",600,600);
  TCanvas FinalEffTightIsoCJPtMVACanvas("FinalEffTightIsoCJPtMVA","",600,600);



  TCanvas MatchedLooseIsoCJdRCanvas("MatchedLooseIsoCJdR","",600,600);
  TCanvas MatchedMedIsoCJdRCanvas("MatchedMedIsoCJdR","",600,600);
  TCanvas MatchedTightIsoCJdRCanvas("MatchedTightIsoCJdR","",600,600);
  TCanvas MatchedDMFindCJdRCanvas("MatchedDMFindCJdR","",600,600);
  TCanvas MatchedCJdRCanvas("MatchedCJdR","",600,600);
  TCanvas FinalEffLooseIsoCJdRCanvas("FinalEffLooseIsoCJdR","",600,600);
  TCanvas FinalEffMedIsoCJdRCanvas("FinalEffMedIsoCJdR","",600,600);
  TCanvas FinalEffTightIsoCJdRCanvas("FinalEffTightIsoCJdR","",600,600);
  TCanvas FinalEffDMFindCJdRCanvas("FinalEffDMFindCJdR","",600,600);
  TCanvas FinalEffLooseIsoCJdROverDMCanvas("FinalEffLooseIsoCJdROverDM","",600,600);
  TCanvas FinalEffMedIsoCJdROverDMCanvas("FinalEffMedIsoCJdROverDM","",600,600);
  TCanvas FinalEffTightIsoCJdROverDMCanvas("FinalEffTightIsoCJdROverDM","",600,600);
  
  TCanvas MatchedLooseIsoCJdRMVACanvas("MatchedLooseIsoCJdRMVA","",600,600);
  TCanvas MatchedMedIsoCJdRMVACanvas("MatchedMedIsoCJdRMVA","",600,600);
  TCanvas MatchedTightIsoCJdRMVACanvas("MatchedTightIsoCJdRMVA","",600,600);
  TCanvas MatchedCJdRMVACanvas("MatchedCJdRMVA","",600,600);
  TCanvas FinalEffLooseIsoCJdRMVACanvas("FinalEffLooseIsoCJdRMVA","",600,600);
  TCanvas FinalEffMedIsoCJdRMVACanvas("FinalEffMedIsoCJdRMVA","",600,600);
  TCanvas FinalEffTightIsoCJdRMVACanvas("FinalEffTightIsoCJdRMVA","",600,600);



  TCanvas MatchedLooseIsoCJPtGenCanvas("MatchedLooseIsoCJPtGen","",600,600);
  TCanvas MatchedMedIsoCJPtGenCanvas("MatchedMedIsoCJPtGen","",600,600);
  TCanvas MatchedTightIsoCJPtGenCanvas("MatchedTightIsoCJPtGen","",600,600);
  TCanvas MatchedDMFindCJPtGenCanvas("MatchedDMFindCJPtGen","",600,600);
  TCanvas MatchedCJPtGenCanvas("MatchedCJPtGen","",600,600);
  TCanvas FinalEffLooseIsoCJPtGenCanvas("FinalEffLooseIsoCJPtGen","",600,600);
  TCanvas FinalEffMedIsoCJPtGenCanvas("FinalEffMedIsoCJPtGen","",600,600);
  TCanvas FinalEffTightIsoCJPtGenCanvas("FinalEffTightIsoCJPtGen","",600,600);
  TCanvas FinalEffDMFindCJPtGenCanvas("FinalEffDMFindCJPtGen","",600,600);
  TCanvas FinalEffLooseIsoCJPtGenOverDMCanvas("FinalEffLooseIsoCJPtGenOverDM","",600,600);
  TCanvas FinalEffMedIsoCJPtGenOverDMCanvas("FinalEffMedIsoCJPtGenOverDM","",600,600);
  TCanvas FinalEffTightIsoCJPtGenOverDMCanvas("FinalEffTightIsoCJPtGenOverDM","",600,600);

  TCanvas MatchedLooseIsoCJPtGenMVACanvas("MatchedLooseIsoCJPtGenMVA","",600,600);
  TCanvas MatchedMedIsoCJPtGenMVACanvas("MatchedMedIsoCJPtGenMVA","",600,600);
  TCanvas MatchedTightIsoCJPtGenMVACanvas("MatchedTightIsoCJPtGenMVA","",600,600);
  TCanvas MatchedCJPtGenMVACanvas("MatchedCJPtGenMVA","",600,600);
  TCanvas FinalEffLooseIsoCJPtGenMVACanvas("FinalEffLooseIsoCJPtGenMVA","",600,600);
  TCanvas FinalEffMedIsoCJPtGenMVACanvas("FinalEffMedIsoCJPtGenMVA","",600,600);
  TCanvas FinalEffTightIsoCJPtGenMVACanvas("FinalEffTightIsoCJPtGenMVA","",600,600);
  
 


  TCanvas OneProngDMCJPtCanvas("OneProngDMCJPt","",600,600);
  TCanvas OneProngOnePizDMCJPtCanvas("OneProngOnePizDMCJPt","",600,600);
  TCanvas OneProngTwoPizDMCJPtCanvas("OneProngTwoPizDMCJPt","",600,600);
  TCanvas ThreeProngDMCJPtCanvas("ThreeProngDMCJPt","",600,600);
  TCanvas MatchedOneProngCJPtCanvas("MatchedOneProngCJPtCanvas","",600,600);
  TCanvas MatchedOneProngOnePizCJPtCanvas("MatchedOneProngOnePizCJPtCanvas","",600,600);
  TCanvas MatchedOneProngTwoPizCJPtCanvas("MatchedOneProngTwoPizCJPtCanvas","",600,600);
  TCanvas MatchedThreeProngCJPtCanvas("MatchedThreeProngCJPtCanvas","",600,600);
  TCanvas FinalOneProngDMCJPtCanvas("FinalOneProngDMCJPt","",600,600);
  TCanvas FinalOneProngOnePizDMCJPtCanvas("FinalOneProngOnePizDMCJPt","",600,600);
  TCanvas FinalOneProngTwoPizDMCJPtCanvas("FinalOneProngTwoPizDMCJPt","",600,600);
  TCanvas FinalThreeProngDMCJPtCanvas("FinalThreeProngDMCJPt","",600,600);


  FinalEffLooseIsoCJPt_->Divide(MatchedLooseIsoCJPt_, MatchedCJPt_);
  FinalEffMedIsoCJPt_->Divide(MatchedMedIsoCJPt_,     MatchedCJPt_);
  FinalEffTightIsoCJPt_->Divide(MatchedTightIsoCJPt_, MatchedCJPt_);
  FinalEffDMFindCJPt_->Divide(MatchedDMFindCJPt_,     MatchedCJPt_);
  FinalEffLooseIsoCJPtOverDM_->Divide(MatchedLooseIsoCJPt_, MatchedDMFindCJPt_);
  FinalEffMedIsoCJPtOverDM_->Divide(MatchedMedIsoCJPt_,     MatchedDMFindCJPt_);
  FinalEffTightIsoCJPtOverDM_->Divide(MatchedTightIsoCJPt_, MatchedDMFindCJPt_);

  FinalEffLooseIsoCJPtMVA_->Divide(MatchedLooseIsoCJPtMVA_, MatchedCJPtMVA_);
  FinalEffMedIsoCJPtMVA_->Divide(MatchedMedIsoCJPtMVA_,     MatchedCJPtMVA_);
  FinalEffTightIsoCJPtMVA_->Divide(MatchedTightIsoCJPtMVA_, MatchedCJPtMVA_);
  


  FinalEffLooseIsoCJdR_->Divide(MatchedLooseIsoCJdR_, MatchedCJdR_);
  FinalEffMedIsoCJdR_->Divide(MatchedMedIsoCJdR_,     MatchedCJdR_);
  FinalEffTightIsoCJdR_->Divide(MatchedTightIsoCJdR_, MatchedCJdR_);
  FinalEffDMFindCJdR_->Divide(MatchedDMFindCJdR_,     MatchedCJdR_);
  FinalEffLooseIsoCJdROverDM_->Divide(MatchedLooseIsoCJdR_, MatchedDMFindCJdR_);
  FinalEffMedIsoCJdROverDM_->Divide(MatchedMedIsoCJdR_,     MatchedDMFindCJdR_);
  FinalEffTightIsoCJdROverDM_->Divide(MatchedTightIsoCJdR_, MatchedDMFindCJdR_);

  FinalEffLooseIsoCJdRMVA_->Divide(MatchedLooseIsoCJdRMVA_, MatchedCJdRMVA_);
  FinalEffMedIsoCJdRMVA_->Divide(MatchedMedIsoCJdRMVA_,     MatchedCJdRMVA_);
  FinalEffTightIsoCJdRMVA_->Divide(MatchedTightIsoCJdRMVA_, MatchedCJdRMVA_);



  FinalEffLooseIsoCJPtGen_->Divide(MatchedLooseIsoCJPtGen_, MatchedCJPtGen_);
  FinalEffMedIsoCJPtGen_->Divide(MatchedMedIsoCJPtGen_,     MatchedCJPtGen_);
  FinalEffTightIsoCJPtGen_->Divide(MatchedTightIsoCJPtGen_, MatchedCJPtGen_);
  FinalEffDMFindCJPtGen_->Divide(MatchedDMFindCJPtGen_,     MatchedCJPtGen_);
  FinalEffLooseIsoCJPtGenOverDM_->Divide(MatchedLooseIsoCJPtGen_, MatchedDMFindCJPtGen_);
  FinalEffMedIsoCJPtGenOverDM_->Divide(MatchedMedIsoCJPtGen_,     MatchedDMFindCJPtGen_);
  FinalEffTightIsoCJPtGenOverDM_->Divide(MatchedTightIsoCJPtGen_, MatchedDMFindCJPtGen_);

  FinalEffLooseIsoCJPtGenMVA_->Divide(MatchedLooseIsoCJPtGenMVA_, MatchedCJPtGenMVA_);
  FinalEffMedIsoCJPtGenMVA_->Divide(MatchedMedIsoCJPtGenMVA_,     MatchedCJPtGenMVA_);
  FinalEffTightIsoCJPtGenMVA_->Divide(MatchedTightIsoCJPtGenMVA_, MatchedCJPtGenMVA_);



  FinalOneProngDMCJPt_->Divide(OneProngDMCJPt_, MatchedOneProngCJPt_);
  FinalOneProngOnePizDMCJPt_->Divide(OneProngOnePizDMCJPt_, MatchedOneProngOnePizCJPt_);
  FinalOneProngTwoPizDMCJPt_->Divide(OneProngTwoPizDMCJPt_, MatchedOneProngTwoPizCJPt_);
  FinalThreeProngDMCJPt_->Divide(ThreeProngDMCJPt_, MatchedThreeProngCJPt_);

std::cout << "<----------------Declared Canvases-------------->" << std::endl;

  //Format the 1D plots and Draw (canvas, hist, grid, log y, log z, color, size, style, xAxisTitle, xTitleSize, xLabelSize, xTitleOffSet, yAxisTitle, yTitleSize, yLabelSize, yTitleOffset)
  VariousFunctions::formatAndDrawCanvasAndHist1D(NEventsCanvas, NEvents_,
         1, 0, 0, kBlack, 7, 20, "", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MuRemovedJetPtCanvas, MuRemovedJetPt_,
         1, 0, 0, kBlack, 7, 20, "Pt_{Jet with #mu Removed}", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauMuTauHaddRCanvas, TauMuTauHaddR_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(#tau_{mu} + #tau_{H})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTauMuPtCanvas, MatchedTauMuPt_,
         1, 0, 0, kBlack, 7, 20, "p_{T}(#tau_{mu})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTauHadPtCanvas, MatchedTauHadPt_,
         1, 0, 0, kBlack, 7, 20, "p_{T}(#tau_{H})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTauHadEtaCanvas, MatchedTauHadEta_,
         1, 0, 0, kBlack, 7, 20, "#eta(#tau_{H})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedHTCanvas, MatchedHT_,
         1, 0, 0, kBlack, 7, 20, "H_{T} of ak4PFJet with p_{T} > 20", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedBDiscCSVCanvas, MatchedBDiscCSV_,
         1, 0, 0, kBlack, 7, 20, "CSV B Discriminant Value", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NEventsCutsCanvas, NEventsCuts_,
         1, 0, 0, kBlack, 7, 20, "", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NConstituentsCJCanvas, NConstituentsCJ_,
         1, 0, 0, kBlack, 7, 20, "Number of Constituents", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NTauDecayModeGENCanvas, NTauDecayModeGEN_,
         1, 0, 0, kBlack, 7, 20, "TauDecayMode GEN CJ", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(NTauDecayModeRECOCanvas, NTauDecayModeRECO_,
         1, 0, 0, kBlack, 7, 20, "TauDecayMode RECO CJ", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(DRTauHadRemovedMuCanvas, DRTauHadRemovedMu_,
         1, 0, 0, kBlack, 7, 20, "dR(#tau_had, removed #mu)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(GenDiTaudRIfRemovedMuCanvas, GenDiTaudRIfRemovedMu_,
         1, 0, 0, kBlack, 7, 20, "Gen Matched Iso Pt Sum", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(GenMatchedIsoPtSumCJCanvas, GenMatchedIsoPtSumCJ_,
         1, 0, 0, kBlack, 7, 20, "dR(Gen  #tau#tau), if Mu is Removed", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist2D(TotalMuvsSoftMuCanvas, TotalMuvsSoftMu_,
         1, 0, 0, kBlack, 7, 20, "nRECO #tau's", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(DRGenTausvsRECOTauRemovedMuCanvas, DRGenTausvsRECOTauRemovedMu_,
         1, 0, 0, kBlack, 7, 20, "dR(Gen  #tau#tau), if Mu is Removed", .04, .04, 1.1, "dR(#tau_had, removed #mu)", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(GenDiTaudRvsCJDiTaudRCanvas, GenDiTaudRvsCJDiTaudR_,
         1, 0, 0, kBlack, 7, 20, "gen #DeltaR(#tau#tau)", .04, .04, 1.1, "CJ #DeltaR(#tau#tau)", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(TauHadnConstvsPtCanvas, TauHadnConstvsPt_,
         1, 0, 0, kBlack, 7, 20, "# Constituents(#tau_{H})", .04, .04, 1.1, "Pt(#tau_{H})", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(RECOPtvsGENVisiblePtCJCanvas, RECOPtvsGENVisiblePtCJ_,
         1, 0, 0, kBlack, 7, 20, "RECO p_{T}(#tau)", .04, .04, 1.1, "GEN Visible p_{T}(#tau)", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(RECODecayModevsGENDecayModeCJCanvas, RECODecayModevsGENDecayModeCJ_,
         1, 0, 0, kBlack, 7, 20, "RECO DecayMode No CJ", .04, .04, 1.1, "GEN Decay Mode No CJ", .04, .04, 1.6, "", .04, .04, 1.0);



  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoCJPtCanvas, MatchedLooseIsoCJPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoCJPtCanvas, MatchedMedIsoCJPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoCJPtCanvas, MatchedTightIsoCJPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedDMFindCJPtCanvas, MatchedDMFindCJPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(CleanJets Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedCJPtCanvas, MatchedCJPt_,
	 1, 0, 0, kBlack, 7, 20, "Pt(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoCJPtCanvas, FinalEffLooseIsoCJPt_,
	 1, 0, 0, kBlack, 1, 20, "Pt(CleanJets Matched + Loose Iso + DM) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoCJPtCanvas, FinalEffMedIsoCJPt_,
	 1, 0, 0, kBlack, 1, 20, "Pt(CleanJets Matched + Med Iso + DM) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoCJPtCanvas, FinalEffTightIsoCJPt_,
	 1, 0, 0, kBlack, 1, 20, "Pt(CleanJets Matched + Tight Iso + DM) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffDMFindCJPtCanvas, FinalEffDMFindCJPt_,
	 1, 0, 0, kBlack, 1, 20, "Pt(CleanJets Matched + DecayModeFinding) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoCJPtOverDMCanvas, FinalEffLooseIsoCJPtOverDM_,
         1, 0, 0, kBlack, 1, 20, "Pt(CleanJets Matched + Loose Iso + DM) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoCJPtOverDMCanvas, FinalEffMedIsoCJPtOverDM_,
         1, 0, 0, kBlack, 1, 20, "Pt(CleanJets Matched + Med Iso + DM) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoCJPtOverDMCanvas, FinalEffTightIsoCJPtOverDM_,
         1, 0, 0, kBlack, 1, 20, "Pt(CleanJets Matched + Tight Iso + DM) / Pt(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoCJPtMVACanvas, MatchedLooseIsoCJPtMVA_,
         1, 0, 0, kBlack, 7, 20, "PtMVA(CleanJets Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoCJPtMVACanvas, MatchedMedIsoCJPtMVA_,
         1, 0, 0, kBlack, 7, 20, "PtMVA(CleanJets Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoCJPtMVACanvas, MatchedTightIsoCJPtMVA_,
         1, 0, 0, kBlack, 7, 20, "PtMVA(CleanJets Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedCJPtMVACanvas, MatchedCJPtMVA_,
         1, 0, 0, kBlack, 7, 20, "PtMVA(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoCJPtMVACanvas, FinalEffLooseIsoCJPtMVA_,
         1, 0, 0, kBlack, 1, 20, "PtMVA(CleanJets Matched + Loose Iso + DM) / PtMVA(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoCJPtMVACanvas, FinalEffMedIsoCJPtMVA_,
         1, 0, 0, kBlack, 1, 20, "PtMVA(CleanJets Matched + Med Iso + DM) / PtMVA(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoCJPtMVACanvas, FinalEffTightIsoCJPtMVA_,
         1, 0, 0, kBlack, 1, 20, "PtMVA(CleanJets Matched + Tight Iso + DM) / PtMVA(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);




  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoCJdRCanvas, MatchedLooseIsoCJdR_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoCJdRCanvas, MatchedMedIsoCJdR_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoCJdRCanvas, MatchedTightIsoCJdR_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedDMFindCJdRCanvas, MatchedDMFindCJdR_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedCJdRCanvas, MatchedCJdR_,
	 1, 0, 0, kBlack, 7, 20, "#DeltaR(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoCJdRCanvas, FinalEffLooseIsoCJdR_,
	 1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Loose Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoCJdRCanvas, FinalEffMedIsoCJdR_,
	 1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Med Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoCJdRCanvas, FinalEffTightIsoCJdR_,
	 1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Tight Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffDMFindCJdRCanvas, FinalEffDMFindCJdR_,
	 1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + DecayModeFinding) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoCJdROverDMCanvas, FinalEffLooseIsoCJdROverDM_,
         1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Loose Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoCJdROverDMCanvas, FinalEffMedIsoCJdROverDM_,
         1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Med Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoCJdROverDMCanvas, FinalEffTightIsoCJdROverDM_,
         1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Tight Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoCJdRMVACanvas, MatchedLooseIsoCJdRMVA_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoCJdRMVACanvas, MatchedMedIsoCJdRMVA_, 
         1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoCJdRMVACanvas, MatchedTightIsoCJdRMVA_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(CleanJets Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedCJdRMVACanvas, MatchedCJdRMVA_,
         1, 0, 0, kBlack, 7, 20, "#DeltaR(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoCJdRMVACanvas, FinalEffLooseIsoCJdRMVA_,
         1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Loose Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoCJdRMVACanvas, FinalEffMedIsoCJdRMVA_,
         1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Med Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoCJdRMVACanvas, FinalEffTightIsoCJdRMVA_,
         1, 0, 0, kBlack, 1, 20, "#DeltaR(CleanJets Matched + Tight Iso + DM) / #DeltaR(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);



  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoCJPtGenCanvas, MatchedLooseIsoCJPtGen_,
         1, 0, 0, kBlack, 7, 20, "PtGen(CleanJets Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoCJPtGenCanvas, MatchedMedIsoCJPtGen_,
         1, 0, 0, kBlack, 7, 20, "PtGen(CleanJets Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoCJPtGenCanvas, MatchedTightIsoCJPtGen_,
         1, 0, 0, kBlack, 7, 20, "PtGen(CleanJets Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedDMFindCJPtGenCanvas, MatchedDMFindCJPtGen_,
         1, 0, 0, kBlack, 7, 20, "PtGen(CleanJets Matched + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedCJPtGenCanvas, MatchedCJPtGen_,
         1, 0, 0, kBlack, 7, 20, "PtGen(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoCJPtGenCanvas, FinalEffLooseIsoCJPtGen_,
         1, 0, 0, kBlack, 1, 20, "PtGen(CleanJets Matched + Loose Iso + DM) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoCJPtGenCanvas, FinalEffMedIsoCJPtGen_,
         1, 0, 0, kBlack, 1, 20, "PtGen(CleanJets Matched + Med Iso + DM) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoCJPtGenCanvas, FinalEffTightIsoCJPtGen_,
         1, 0, 0, kBlack, 1, 20, "PtGen(CleanJets Matched + Tight Iso + DM) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffDMFindCJPtGenCanvas, FinalEffDMFindCJPtGen_,
         1, 0, 0, kBlack, 1, 20, "PtGen(CleanJets Matched + DecayModeFinding) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoCJPtGenOverDMCanvas, FinalEffLooseIsoCJPtGenOverDM_,
         1, 0, 0, kBlack, 1, 20, "PtGen(CleanJets Matched + Loose Iso + DM) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoCJPtGenOverDMCanvas, FinalEffMedIsoCJPtGenOverDM_,
         1, 0, 0, kBlack, 1, 20, "PtGen(CleanJets Matched + Med Iso + DM) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoCJPtGenOverDMCanvas, FinalEffTightIsoCJPtGenOverDM_,
         1, 0, 0, kBlack, 1, 20, "PtGen(CleanJets Matched + Tight Iso + DM) / PtGen(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedLooseIsoCJPtGenMVACanvas, MatchedLooseIsoCJPtGenMVA_,
         1, 0, 0, kBlack, 7, 20, "PtGenMVA(CleanJets Matched + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedMedIsoCJPtGenMVACanvas, MatchedMedIsoCJPtGenMVA_,
         1, 0, 0, kBlack, 7, 20, "PtGenMVA(CleanJets Matched + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedTightIsoCJPtGenMVACanvas, MatchedTightIsoCJPtGenMVA_,
         1, 0, 0, kBlack, 7, 20, "PtGenMVA(CleanJets Matched + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedCJPtGenMVACanvas, MatchedCJPtGenMVA_,
         1, 0, 0, kBlack, 7, 20, "PtGenMVA(Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoCJPtGenMVACanvas, FinalEffLooseIsoCJPtGenMVA_,
         1, 0, 0, kBlack, 1, 20, "PtGenMVA(CleanJets Matched + Loose Iso + DM) / PtGenMVA(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoCJPtGenMVACanvas, FinalEffMedIsoCJPtGenMVA_,
         1, 0, 0, kBlack, 1, 20, "PtGenMVA(CleanJets Matched + Med Iso + DM) / PtGenMVA(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoCJPtGenMVACanvas, FinalEffTightIsoCJPtGenMVA_,
         1, 0, 0, kBlack, 1, 20, "PtGenMVA(CleanJets Matched + Tight Iso + DM) / PtGenMVA(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);


  
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngDMCJPtCanvas, OneProngDMCJPt_, 
         1, 0, 0, kBlack, 7, 20, "1 Prong p_{T}(CleanJets Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngOnePizDMCJPtCanvas, OneProngOnePizDMCJPt_, 
         1, 0, 0, kBlack, 7, 20, "1 Prong + 1 #pi^{0} p_{T}(CleanJets Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngTwoPizDMCJPtCanvas, OneProngTwoPizDMCJPt_, 
         1, 0, 0, kBlack, 7, 20, "1 Prong + 2 #pi^{0} p_{T}(CleanJets Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(ThreeProngDMCJPtCanvas, ThreeProngDMCJPt_, 
         1, 0, 0, kBlack, 7, 20, "3 Prong p_{T}(CleanJets Matched + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedOneProngCJPtCanvas, MatchedOneProngCJPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong p_{T}(CleanJets Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedOneProngOnePizCJPtCanvas, MatchedOneProngOnePizCJPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 1 #pi^{0} p_{T}(CleanJets Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedOneProngTwoPizCJPtCanvas, MatchedOneProngTwoPizCJPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 2 #pi^{0} p_{T}(CleanJets Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(MatchedThreeProngCJPtCanvas, MatchedThreeProngCJPt_,
         1, 0, 0, kBlack, 7, 20, "3 Prong p_{T}(CleanJets Matched)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngDMCJPtCanvas, FinalOneProngDMCJPt_, 
         1, 0, 0, kBlack, 1, 20, "1 Prong p_{T}(CleanJets Matched + DM) / 1 Prong p_{T}(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngOnePizDMCJPtCanvas, FinalOneProngOnePizDMCJPt_, 
         1, 0, 0, kBlack, 1, 20, "1 Prong + 1 #pi^{0} p_{T}(CleanJets Matched + DM) / 1 Prong + 1 #pi^{0} p_{T}(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngTwoPizDMCJPtCanvas, FinalOneProngTwoPizDMCJPt_, 
         1, 0, 0, kBlack, 1, 20, "1 Prong + 2 #pi^{0} p_{T}(CleanJets Matched + DM) / 1 Prong + 2 #pi^{0} p_{T}(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalThreeProngDMCJPtCanvas, FinalThreeProngDMCJPt_, 
         1, 0, 0, kBlack, 1, 20, "3 Prong p_{T}(CleanJets Matched + DM) / 3 Prong p_{T}(CleanJets Matched)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);



std::cout << "after formatting" << std::endl;
  
  FinalEffLooseIsoCJPt_->SetLineColor(kRed);
  FinalEffMedIsoCJPt_->SetLineColor(kRed);
  FinalEffTightIsoCJPt_->SetLineColor(kRed);
  FinalEffDMFindCJPt_->SetLineColor(kRed);
  FinalEffLooseIsoCJPtOverDM_->SetLineColor(kRed);
  FinalEffMedIsoCJPtOverDM_->SetLineColor(kRed);
  FinalEffTightIsoCJPtOverDM_->SetLineColor(kRed);

  FinalEffLooseIsoCJPtMVA_->SetLineColor(kRed);
  FinalEffMedIsoCJPtMVA_->SetLineColor(kRed);
  FinalEffTightIsoCJPtMVA_->SetLineColor(kRed);



  FinalEffLooseIsoCJdR_->SetLineColor(kRed);
  FinalEffMedIsoCJdR_->SetLineColor(kRed);
  FinalEffTightIsoCJdR_->SetLineColor(kRed);
  FinalEffDMFindCJdR_->SetLineColor(kRed);
  FinalEffLooseIsoCJdROverDM_->SetLineColor(kRed);
  FinalEffMedIsoCJdROverDM_->SetLineColor(kRed);
  FinalEffTightIsoCJdROverDM_->SetLineColor(kRed);
  
  FinalEffLooseIsoCJdRMVA_->SetLineColor(kRed);
  FinalEffMedIsoCJdRMVA_->SetLineColor(kRed);
  FinalEffTightIsoCJdRMVA_->SetLineColor(kRed);



  FinalEffLooseIsoCJPtGen_->SetLineColor(kRed);
  FinalEffMedIsoCJPtGen_->SetLineColor(kRed);
  FinalEffTightIsoCJPtGen_->SetLineColor(kRed);
  FinalEffDMFindCJPtGen_->SetLineColor(kRed);
  FinalEffLooseIsoCJPtGenOverDM_->SetLineColor(kRed);
  FinalEffMedIsoCJPtGenOverDM_->SetLineColor(kRed);
  FinalEffTightIsoCJPtGenOverDM_->SetLineColor(kRed);

  FinalEffLooseIsoCJPtGenMVA_->SetLineColor(kRed);
  FinalEffMedIsoCJPtGenMVA_->SetLineColor(kRed);
  FinalEffTightIsoCJPtGenMVA_->SetLineColor(kRed);



  FinalOneProngDMCJPt_->SetLineColor(kRed);
  FinalOneProngOnePizDMCJPt_->SetLineColor(kRed);
  FinalOneProngTwoPizDMCJPt_->SetLineColor(kRed);
  FinalThreeProngDMCJPt_->SetLineColor(kRed);

std::cout << "<----------------Formatted Canvases and Histos-------------->" << std::endl;

  //Write output file
  out_->cd();

  NEventsCanvas.Write();
  MuRemovedJetPtCanvas.Write();
  TauMuTauHaddRCanvas.Write();
  MatchedTauMuPtCanvas.Write();
  MatchedTauHadPtCanvas.Write();
  MatchedTauHadEtaCanvas.Write();
  MatchedHTCanvas.Write();
  MatchedBDiscCSVCanvas.Write();
  NEventsCutsCanvas.Write();
  NConstituentsCJCanvas.Write();
  NTauDecayModeGENCanvas.Write();
  NTauDecayModeRECOCanvas.Write();
  DRTauHadRemovedMuCanvas.Write();
  GenDiTaudRIfRemovedMuCanvas.Write();
  GenMatchedIsoPtSumCJCanvas.Write();
  TotalMuvsSoftMuCanvas.Write();
  DRGenTausvsRECOTauRemovedMuCanvas.Write();
  GenDiTaudRvsCJDiTaudRCanvas.Write();
  TauHadnConstvsPtCanvas.Write();
  RECOPtvsGENVisiblePtCJCanvas.Write();
  RECODecayModevsGENDecayModeCJCanvas.Write();



  MatchedLooseIsoCJPtCanvas.Write();
  MatchedMedIsoCJPtCanvas.Write();
  MatchedTightIsoCJPtCanvas.Write();
  MatchedDMFindCJPtCanvas.Write();
  MatchedCJPtCanvas.Write();
  FinalEffLooseIsoCJPtCanvas.Write();
  FinalEffMedIsoCJPtCanvas.Write();
  FinalEffTightIsoCJPtCanvas.Write();
  FinalEffDMFindCJPtCanvas.Write();
  FinalEffLooseIsoCJPtOverDMCanvas.Write();
  FinalEffMedIsoCJPtOverDMCanvas.Write();
  FinalEffTightIsoCJPtOverDMCanvas.Write();
  
  MatchedLooseIsoCJPtMVACanvas.Write();
  MatchedMedIsoCJPtMVACanvas.Write();
  MatchedTightIsoCJPtMVACanvas.Write();
  MatchedCJPtMVACanvas.Write();
  FinalEffLooseIsoCJPtMVACanvas.Write();
  FinalEffMedIsoCJPtMVACanvas.Write();
  FinalEffTightIsoCJPtMVACanvas.Write();



  MatchedLooseIsoCJdRCanvas.Write();
  MatchedMedIsoCJdRCanvas.Write();
  MatchedTightIsoCJdRCanvas.Write();
  MatchedDMFindCJdRCanvas.Write();
  MatchedCJdRCanvas.Write();
  FinalEffLooseIsoCJdRCanvas.Write();
  FinalEffMedIsoCJdRCanvas.Write();
  FinalEffTightIsoCJdRCanvas.Write();
  FinalEffDMFindCJdRCanvas.Write();
  FinalEffLooseIsoCJdROverDMCanvas.Write();
  FinalEffMedIsoCJdROverDMCanvas.Write();
  FinalEffTightIsoCJdROverDMCanvas.Write();
  
  MatchedLooseIsoCJdRMVACanvas.Write();
  MatchedMedIsoCJdRMVACanvas.Write();
  MatchedTightIsoCJdRMVACanvas.Write();
  MatchedCJdRMVACanvas.Write();
  FinalEffLooseIsoCJdRMVACanvas.Write();
  FinalEffMedIsoCJdRMVACanvas.Write();
  FinalEffTightIsoCJdRMVACanvas.Write();



  MatchedLooseIsoCJPtGenCanvas.Write();
  MatchedMedIsoCJPtGenCanvas.Write();
  MatchedTightIsoCJPtGenCanvas.Write();
  MatchedDMFindCJPtGenCanvas.Write();
  MatchedCJPtGenCanvas.Write();
  FinalEffLooseIsoCJPtGenCanvas.Write();
  FinalEffMedIsoCJPtGenCanvas.Write();
  FinalEffTightIsoCJPtGenCanvas.Write();
  FinalEffDMFindCJPtGenCanvas.Write();
  FinalEffLooseIsoCJPtGenOverDMCanvas.Write();
  FinalEffMedIsoCJPtGenOverDMCanvas.Write();
  FinalEffTightIsoCJPtGenOverDMCanvas.Write();
  
  MatchedLooseIsoCJPtGenMVACanvas.Write();
  MatchedMedIsoCJPtGenMVACanvas.Write();
  MatchedTightIsoCJPtGenMVACanvas.Write();
  MatchedCJPtGenMVACanvas.Write();
  FinalEffLooseIsoCJPtGenMVACanvas.Write();
  FinalEffMedIsoCJPtGenMVACanvas.Write();
  FinalEffTightIsoCJPtGenMVACanvas.Write();
  


  OneProngDMCJPtCanvas.Write();
  OneProngOnePizDMCJPtCanvas.Write();
  OneProngTwoPizDMCJPtCanvas.Write();
  ThreeProngDMCJPtCanvas.Write();
  MatchedOneProngCJPtCanvas.Write();
  MatchedOneProngOnePizCJPtCanvas.Write();
  MatchedOneProngTwoPizCJPtCanvas.Write();
  MatchedThreeProngCJPtCanvas.Write();
  FinalOneProngDMCJPtCanvas.Write();
  FinalOneProngOnePizDMCJPtCanvas.Write();
  FinalOneProngTwoPizDMCJPtCanvas.Write();
  FinalThreeProngDMCJPtCanvas.Write();

  out_->Write();
  out_->Close();
std::cout << "DONE" << std::endl;
}//EndJob

// ------------ method called when starting to processes a run  ------------
void GGHAnalyzer_IndivCJ::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void GGHAnalyzer_IndivCJ::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void GGHAnalyzer_IndivCJ::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void GGHAnalyzer_IndivCJ::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

//Delete Memory
void GGHAnalyzer_IndivCJ::reset(const bool doDelete)
{
  if ((doDelete) && (NEvents_ != NULL)) delete NEvents_;
  NEvents_ = NULL;
  if ((doDelete) && (MuRemovedJetPt_ != NULL)) delete MuRemovedJetPt_;
  MuRemovedJetPt_ = NULL;
  if ((doDelete) && (TauMuTauHaddR_ != NULL)) delete TauMuTauHaddR_;
  TauMuTauHaddR_ = NULL;
  if ((doDelete) && (MatchedTauMuPt_ != NULL)) delete MatchedTauMuPt_;
  MatchedTauMuPt_ = NULL;
  if ((doDelete) && (MatchedTauHadPt_ != NULL)) delete MatchedTauHadPt_;
  MatchedTauHadPt_ = NULL;
  if ((doDelete) && (MatchedTauHadEta_ != NULL)) delete MatchedTauHadEta_;
  MatchedTauHadEta_ = NULL;
  if ((doDelete) && (MatchedHT_ != NULL)) delete MatchedHT_;
  MatchedHT_ = NULL;
  if ((doDelete) && (MatchedBDiscCSV_ != NULL)) delete MatchedBDiscCSV_;
  MatchedBDiscCSV_ = NULL;
  if ((doDelete) && (NEventsCuts_ != NULL)) delete NEventsCuts_;
  NEventsCuts_ = NULL;
  if ((doDelete) && (NConstituentsCJ_ != NULL)) delete NConstituentsCJ_;
  NConstituentsCJ_ = NULL;
  if ((doDelete) && (NTauDecayModeGEN_ != NULL)) delete NTauDecayModeGEN_;
  NTauDecayModeGEN_ = NULL;
  if ((doDelete) && (NTauDecayModeRECO_ != NULL)) delete NTauDecayModeRECO_;
  NTauDecayModeRECO_ = NULL;
  if ((doDelete) && (DRTauHadRemovedMu_ != NULL)) delete DRTauHadRemovedMu_;
  DRTauHadRemovedMu_ = NULL;
  if ((doDelete) && (GenDiTaudRIfRemovedMu_ != NULL)) delete GenDiTaudRIfRemovedMu_;
  GenDiTaudRIfRemovedMu_ = NULL;
  if ((doDelete) && (GenMatchedIsoPtSumCJ_ != NULL)) delete GenMatchedIsoPtSumCJ_;
  GenMatchedIsoPtSumCJ_ = NULL;
  if ((doDelete) && (TotalMuvsSoftMu_ != NULL)) delete TotalMuvsSoftMu_;
  TotalMuvsSoftMu_ = NULL;
  if ((doDelete) && (DRGenTausvsRECOTauRemovedMu_ != NULL)) delete DRGenTausvsRECOTauRemovedMu_;
  DRGenTausvsRECOTauRemovedMu_ = NULL;
  if ((doDelete) && (GenDiTaudRvsCJDiTaudR_ != NULL)) delete GenDiTaudRvsCJDiTaudR_;
  GenDiTaudRvsCJDiTaudR_ = NULL;
  if ((doDelete) && (TauHadnConstvsPt_ != NULL)) delete TauHadnConstvsPt_;
  TauHadnConstvsPt_ = NULL;
  if ((doDelete) && (RECOPtvsGENVisiblePtCJ_ != NULL)) delete RECOPtvsGENVisiblePtCJ_;
  RECOPtvsGENVisiblePtCJ_ = NULL;
  if ((doDelete) && (RECODecayModevsGENDecayModeCJ_ != NULL)) delete RECODecayModevsGENDecayModeCJ_;
  RECODecayModevsGENDecayModeCJ_ = NULL;


  if ((doDelete) && (MatchedLooseIsoCJPt_ != NULL)) delete MatchedLooseIsoCJPt_;
  MatchedLooseIsoCJPt_ = NULL;
  if ((doDelete) && (MatchedMedIsoCJPt_ != NULL)) delete MatchedMedIsoCJPt_;
  MatchedMedIsoCJPt_ = NULL;
  if ((doDelete) && (MatchedTightIsoCJPt_ != NULL)) delete MatchedTightIsoCJPt_;
  MatchedTightIsoCJPt_ = NULL;
  if ((doDelete) && (MatchedDMFindCJPt_ != NULL)) delete MatchedDMFindCJPt_;
  MatchedDMFindCJPt_ = NULL;
  if ((doDelete) && (MatchedCJPt_ != NULL)) delete MatchedCJPt_;
  MatchedCJPt_ = NULL;
  if ((doDelete) && (FinalEffLooseIsoCJPt_ != NULL)) delete FinalEffLooseIsoCJPt_;
  FinalEffLooseIsoCJPt_ = NULL;
  if ((doDelete) && (FinalEffMedIsoCJPt_ != NULL)) delete FinalEffMedIsoCJPt_;
  FinalEffMedIsoCJPt_ = NULL;
  if ((doDelete) && (FinalEffTightIsoCJPt_ != NULL)) delete FinalEffTightIsoCJPt_;
  FinalEffTightIsoCJPt_ = NULL;
  if ((doDelete) && (FinalEffDMFindCJPt_ != NULL)) delete FinalEffDMFindCJPt_;
  FinalEffDMFindCJPt_ = NULL;

  if ((doDelete) && (MatchedLooseIsoCJdR_ != NULL)) delete MatchedLooseIsoCJdR_;
  MatchedLooseIsoCJdR_ = NULL;
  if ((doDelete) && (MatchedMedIsoCJdR_ != NULL)) delete MatchedMedIsoCJdR_;
  MatchedMedIsoCJdR_ = NULL;
  if ((doDelete) && (MatchedTightIsoCJdR_ != NULL)) delete MatchedTightIsoCJdR_;
  MatchedTightIsoCJdR_ = NULL;
  if ((doDelete) && (MatchedDMFindCJdR_ != NULL)) delete MatchedDMFindCJdR_;
  MatchedDMFindCJdR_ = NULL;
  if ((doDelete) && (MatchedCJdR_ != NULL)) delete MatchedCJdR_;
  MatchedCJdR_ = NULL;
  if ((doDelete) && (FinalEffLooseIsoCJdR_ != NULL)) delete FinalEffLooseIsoCJdR_;
  FinalEffLooseIsoCJdR_ = NULL;
  if ((doDelete) && (FinalEffMedIsoCJdR_ != NULL)) delete FinalEffMedIsoCJdR_;
  FinalEffMedIsoCJdR_ = NULL;
  if ((doDelete) && (FinalEffTightIsoCJdR_ != NULL)) delete FinalEffTightIsoCJdR_;
  FinalEffTightIsoCJdR_ = NULL;
  if ((doDelete) && (FinalEffDMFindCJdR_ != NULL)) delete FinalEffDMFindCJdR_;
  FinalEffDMFindCJdR_ = NULL;

  if ((doDelete) && (OneProngDMCJPt_ != NULL)) delete OneProngDMCJPt_;
  OneProngDMCJPt_ = NULL;
  if ((doDelete) && (OneProngOnePizDMCJPt_ != NULL)) delete OneProngOnePizDMCJPt_;
  OneProngOnePizDMCJPt_ = NULL;
  if ((doDelete) && (OneProngTwoPizDMCJPt_ != NULL)) delete OneProngTwoPizDMCJPt_;
  OneProngOnePizDMCJPt_ = NULL;
  if ((doDelete) && (ThreeProngDMCJPt_ != NULL)) delete ThreeProngDMCJPt_;
  ThreeProngDMCJPt_ = NULL;
  if ((doDelete) && (MatchedOneProngCJPt_ != NULL)) delete MatchedOneProngCJPt_;
  MatchedOneProngCJPt_ = NULL;
  if ((doDelete) && (MatchedOneProngOnePizCJPt_ != NULL)) delete MatchedOneProngOnePizCJPt_;
  MatchedOneProngOnePizCJPt_ = NULL;
  if ((doDelete) && (MatchedOneProngTwoPizCJPt_ != NULL)) delete MatchedOneProngTwoPizCJPt_;
  MatchedOneProngTwoPizCJPt_ = NULL;
  if ((doDelete) && (MatchedThreeProngCJPt_ != NULL)) delete MatchedThreeProngCJPt_;
  MatchedThreeProngCJPt_ = NULL;
  if ((doDelete) && (FinalOneProngDMCJPt_ != NULL)) delete FinalOneProngDMCJPt_;
  FinalOneProngDMCJPt_ = NULL;
  if ((doDelete) && (FinalOneProngOnePizDMCJPt_ != NULL)) delete FinalOneProngOnePizDMCJPt_;
  FinalOneProngOnePizDMCJPt_ = NULL;
  if ((doDelete) && (FinalOneProngTwoPizDMCJPt_ != NULL)) delete FinalOneProngTwoPizDMCJPt_;
  FinalOneProngTwoPizDMCJPt_ = NULL;
  if ((doDelete) && (FinalThreeProngDMCJPt_ != NULL)) delete FinalThreeProngDMCJPt_;
  FinalThreeProngDMCJPt_ = NULL;


  if ((doDelete) && (MatchedLooseIsoCJPtGen_ != NULL)) delete MatchedLooseIsoCJPtGen_;
  MatchedLooseIsoCJPtGen_ = NULL;
  if ((doDelete) && (MatchedMedIsoCJPtGen_ != NULL)) delete MatchedMedIsoCJPtGen_;
  MatchedMedIsoCJPtGen_ = NULL;
  if ((doDelete) && (MatchedTightIsoCJPtGen_ != NULL)) delete MatchedTightIsoCJPtGen_;
  MatchedTightIsoCJPtGen_ = NULL;
  if ((doDelete) && (MatchedDMFindCJPtGen_ != NULL)) delete MatchedDMFindCJPtGen_;
  MatchedDMFindCJPtGen_ = NULL;
  if ((doDelete) && (MatchedCJPtGen_ != NULL)) delete MatchedCJPtGen_;
  MatchedCJPtGen_ = NULL;
  if ((doDelete) && (FinalEffLooseIsoCJPtGen_ != NULL)) delete FinalEffLooseIsoCJPtGen_;
  FinalEffLooseIsoCJPtGen_ = NULL;
  if ((doDelete) && (FinalEffMedIsoCJPtGen_ != NULL)) delete FinalEffMedIsoCJPtGen_;
  FinalEffMedIsoCJPtGen_ = NULL;
  if ((doDelete) && (FinalEffTightIsoCJPtGen_ != NULL)) delete FinalEffTightIsoCJPtGen_;
  FinalEffTightIsoCJPtGen_ = NULL;
  if ((doDelete) && (FinalEffDMFindCJPtGen_ != NULL)) delete FinalEffDMFindCJPtGen_;
  FinalEffDMFindCJPtGen_ = NULL;


}//void GGHAnalyzer_IndivCJ

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GGHAnalyzer_IndivCJ::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GGHAnalyzer_IndivCJ);

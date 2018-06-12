// -*- C++ -*-
//
// Package:    FakeRateMiniAODGetRates
// Class:      FakeRateMiniAODGetRates
// 
/**\class FakeRateMiniAODGetRates FakeRateMiniAODGetRates.cc Analyzer/src/FakeRateMiniAODGetRates.cc

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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "AnalyzerGeneratorRecoVariousFunctions/VariousFunctions/interface/BTagCalibrationStandalone.h"

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

class FakeRateMiniAODGetRates : public edm::EDAnalyzer {
   public:
      typedef reco::JetFloatAssociation::Container JetBCEnergyRatioCollection;
      explicit FakeRateMiniAODGetRates(const edm::ParameterSet&);
      ~FakeRateMiniAODGetRates();

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
      edm::EDGetTokenT<edm::View<pat::Jet> > jetTag_;
      edm::EDGetTokenT<edm::View<pat::Muon> > muonsTag_;
      edm::EDGetTokenT<edm::View<pat::Tau> > tauTag_;
      std::string vLooseIsoTag_;
      std::string looseIsoTag_;
      std::string medIsoTag_;
      std::string tightIsoTag_;
      std::string vTightIsoTag_;
      std::string decayModeFindingTag_;
      std::string isoRawTag_;
      double mu3dRMin_;
      double mu3dRMax_;
      double tauPtCut_;
      double mu3dROverlapCut_;
//      edm::EDGetTokenT<edm::View<pat::Muon> > mu3Tag_;
      edm::EDGetTokenT<edm::View<pat::Muon> > mu12Tag_;
      bool requireRemovedMuon_;
      bool checkInvMass_;
      double checkInvMassMin_;
      double checkInvMassMax_;
      bool checkBTag_;
      string csvBTag_;
      bool isMC_;
      double xsec_;
      double lumi_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryInfo_;
      edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
      double summedWeights_;
      std::string PileupFileName_;
      std::string _fpIDs_BToF;
      std::string _fpIDs_GH;
      std::string _fpISOs_BToF;
      std::string _fpISOs_GH;
      std::string _fpTrack;
      std::string _fpTrigger_BToF;
      std::string _fpTrigger_GH;
      std::string _fp_LowPt;
      TFile *_fileIDs_BToF;
      TFile *_fileIDs_GH;
      TFile *_fileISOs_BToF;
      TFile *_fileISOs_GH;
      TFile *_fileTrack;
      TFile *_fileTrigger_BToF;
      TFile *_fileTrigger_GH;
      TFile *_file_LowPt;
      TH2F *IDsWeight_BToF;
      TH2F *IDsWeight_GH;
      TH2F *ISOsWeight_BToF;
      TH2F *ISOsWeight_GH;
      TH2F *IDsWeight_LowPt_ID;
      TH2F *ISOsWeight_LowPt_ISO;
      TGraph *TrackWeight;
      TH2F *TriggerWeight_BToF;
      TH2F *TriggerWeight_GH;
      struct TrackProperties{
        Double_t x;
        Double_t y;
        Double_t errx_up;
        Double_t errx_down;
        Double_t erry_up;
        Double_t erry_down;
      };
      std::list<TrackProperties> TrackCorr;
      struct Gauss{
        Double_t meanData;
        Double_t sigmaData;
        Double_t meanMC;
        Double_t sigmaMC;
        Double_t ptMin;
        Double_t ptMax;
      };
      struct Gauss GaussU1Corr[5];
      struct Gauss GaussU2Corr[5];
      TTree *METRecoil_;
      Double_t totalWeights_;
      Double_t f_z1_pt;
      Double_t f_u1;
      Double_t f_u2;


      //Histograms
      TH1F* NEvents_;   
      TH1F* InvMassTauMuMu1_;
      TH1F* InvMassMu1Mu2_;
      TH1F* InvMassTauMuMu2_;

      TH2F* EtavsPtTauLooseIso_;
      TH2F* EtavsPtTauMedIso_;
      TH2F* EtavsPtTauTightIso_;
      TH2F* EtavsPtTauDMFind_;
      TH2F* EtavsPtJet_;
      TH2F* EtavsPtJetSoftMuon_;
      TH2F* EtavsPtJetSoftMuon_noMu_;

      TH1F* TauLooseIsoEta_;
      TH1F* TauMedIsoEta_;
      TH1F* TauTightIsoEta_;
      TH1F* TauDMFindEta_;
      TH1F* JetEta_;
      TH1F* JetEtaWithSoftMuon_;
      TH1F* JetEtaWithSoftMuon_noMu_;
      TGraphAsymmErrors* FinalEffLooseIsoEta_;
      TGraphAsymmErrors* FinalEffMedIsoEta_;
      TGraphAsymmErrors* FinalEffTightIsoEta_;
      TGraphAsymmErrors* FinalEffDMFindEta_;

      TH1F* TauLooseIsoPt_;
      TH1F* TauMedIsoPt_;
      TH1F* TauTightIsoPt_;
      TH1F* TauDMFindPt_;
      TH1F* JetPt_;
      TH1F* JetPtWithSoftMuon_;
      TH1F* JetPtWithSoftMuon_noMu_;
      TGraphAsymmErrors* FinalEffLooseIsoPt_;
      TGraphAsymmErrors* FinalEffMedIsoPt_;
      TGraphAsymmErrors* FinalEffTightIsoPt_;
      TGraphAsymmErrors* FinalEffDMFindPt_;


     TH1F* OneProngDMEta_;
     TH1F* OneProngOnePizDMEta_;
     TH1F* OneProngTwoPizDMEta_;
     TH1F* ThreeProngDMEta_;
     TH1F* OneProngDMPt_;
     TH1F* OneProngOnePizDMPt_;
     TH1F* OneProngTwoPizDMPt_;
     TH1F* ThreeProngDMPt_;
     TGraphAsymmErrors* FinalOneProngDMEta_;
     TGraphAsymmErrors* FinalOneProngOnePizDMEta_;
     TGraphAsymmErrors* FinalOneProngTwoPizDMEta_;
     TGraphAsymmErrors* FinalThreeProngDMEta_;
     TGraphAsymmErrors* FinalOneProngDMPt_;
     TGraphAsymmErrors* FinalOneProngOnePizDMPt_;
     TGraphAsymmErrors* FinalOneProngTwoPizDMPt_;
     TGraphAsymmErrors* FinalThreeProngDMPt_;

     TH1F* OneProngMedIsoEta_;
     TH1F* OneProngOnePizMedIsoEta_;
     TH1F* OneProngTwoPizMedIsoEta_;
     TH1F* ThreeProngMedIsoEta_;
     TH1F* OneProngMedIsoPt_;
     TH1F* OneProngOnePizMedIsoPt_;
     TH1F* OneProngTwoPizMedIsoPt_;
     TH1F* ThreeProngMedIsoPt_;
     TGraphAsymmErrors* FinalOneProngMedIsoEta_;
     TGraphAsymmErrors* FinalOneProngOnePizMedIsoEta_;
     TGraphAsymmErrors* FinalOneProngTwoPizMedIsoEta_;
     TGraphAsymmErrors* FinalThreeProngMedIsoEta_;
     TGraphAsymmErrors* FinalOneProngMedIsoPt_;
     TGraphAsymmErrors* FinalOneProngOnePizMedIsoPt_;
     TGraphAsymmErrors* FinalOneProngTwoPizMedIsoPt_;
     TGraphAsymmErrors* FinalThreeProngMedIsoPt_;


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
FakeRateMiniAODGetRates::FakeRateMiniAODGetRates(const edm::ParameterSet& iConfig):
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  jetTag_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetTag"))),
  muonsTag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonsTag"))),
  tauTag_(consumes<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("tauTag"))),
  vLooseIsoTag_(iConfig.getParameter<std::string>("vLooseIsoTag")),
  looseIsoTag_(iConfig.getParameter<std::string>("looseIsoTag")),
  medIsoTag_(iConfig.getParameter<std::string>("medIsoTag")),
  tightIsoTag_(iConfig.getParameter<std::string>("tightIsoTag")),
  vTightIsoTag_(iConfig.getParameter<std::string>("vTightIsoTag")),
  decayModeFindingTag_(iConfig.getParameter<std::string>("decayModeFindingTag")),
  isoRawTag_(iConfig.getParameter<std::string>("isoRawTag")),
  mu3dRMin_(iConfig.getParameter<double>("mu3dRMin")),
  mu3dRMax_(iConfig.getParameter<double>("mu3dRMax")),
  tauPtCut_(iConfig.getParameter<double>("tauPtCut")),
  mu3dROverlapCut_(iConfig.getParameter<double>("mu3dROverlapCut")),
//  mu3Tag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("mu3Tag"))),
  mu12Tag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("mu12Tag"))),
  requireRemovedMuon_(iConfig.getParameter<bool>("requireRemovedMuon")),
  checkInvMass_(iConfig.getParameter<bool>("checkInvMass")),
  checkInvMassMin_(iConfig.getParameter<double>("checkInvMassMin")),
  checkInvMassMax_(iConfig.getParameter<double>("checkInvMassMax")),
  checkBTag_(iConfig.getParameter<bool>("checkBTag")),
  csvBTag_(iConfig.getParameter<std::string>("csvBTag")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  xsec_(iConfig.getParameter<double>("xsec")),
  lumi_(iConfig.getParameter<double>("lumi")),
  pileupSummaryInfo_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummaryInfo"))),
  genEventInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfoToken"))),
  summedWeights_(iConfig.getParameter<double>("summedWeights")),
  PileupFileName_(iConfig.getParameter<std::string>("PileupFileName")),
  _fpIDs_BToF(iConfig.getParameter<std::string>("fpIDs_BToF")),
  _fpIDs_GH(iConfig.getParameter<std::string>("fpIDs_GH")),
  _fpISOs_BToF(iConfig.getParameter<std::string>("fpISOs_BToF")),
  _fpISOs_GH(iConfig.getParameter<std::string>("fpISOs_GH")),
  _fpTrack(iConfig.getParameter<std::string>("fpTrack")),
  _fpTrigger_BToF(iConfig.getParameter<std::string>("fpTrigger_BToF")),
  _fpTrigger_GH(iConfig.getParameter<std::string>("fpTrigger_GH")),
  _fp_LowPt(iConfig.getParameter<std::string>("fp_LowPt"))

{
  _fileIDs_BToF = new TFile(_fpIDs_BToF.c_str());
  _fileIDs_GH = new TFile(_fpIDs_GH.c_str());
  _fileISOs_BToF = new TFile(_fpISOs_BToF.c_str());
  _fileISOs_GH = new TFile(_fpISOs_GH.c_str());
  _fileTrack = new TFile(_fpTrack.c_str());
  _fileTrigger_BToF = new TFile(_fpTrigger_BToF.c_str());
  _fileTrigger_GH = new TFile(_fpTrigger_GH.c_str());
  _file_LowPt = new TFile(_fp_LowPt.c_str());


  IDsWeight_BToF =  (TH2F*)_fileIDs_BToF->Get("EtavsPtMedium2016ID");
  IDsWeight_GH =  (TH2F*)_fileIDs_GH->Get("EtavsPtMedium2016ID");
  ISOsWeight_BToF = (TH2F*)_fileISOs_BToF->Get("EtavsPtLooseISO");
  ISOsWeight_GH = (TH2F*)_fileISOs_GH->Get("EtavsPtLooseISO");
  TrackWeight = (TGraph*)_fileTrack->Get("ratio_eff_aeta_dr030e030_corr");
  TriggerWeight_BToF = (TH2F*)_fileTrigger_BToF->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio");
  TriggerWeight_GH = (TH2F*)_fileTrigger_GH->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio");
  IDsWeight_LowPt_ID = (TH2F*)_file_LowPt->Get("hist_EtavsPtLooseID_DatatoMC");
  ISOsWeight_LowPt_ISO = (TH2F*)_file_LowPt->Get("hist_EtavsPtLooseISO_DatatoMC");
  std::cout << "ISOsWeight_LowPt_ISO=" << ISOsWeight_LowPt_ISO << std::endl;
  std::cout << "IDsWeight_LowPt_ID-" << IDsWeight_LowPt_ID << std::endl;
  std::cout << "IDsWeight_BToF =" << IDsWeight_BToF << std::endl;
  std::cout << "IDsWeight_GH =" << IDsWeight_GH << std::endl;
  std::cout << "ISOsWeight_BToF =" << ISOsWeight_BToF << std::endl;
  std::cout << "ISOsWeight_GH =" << ISOsWeight_GH << std::endl;
  std::cout << "TrackWeight =" << TrackWeight << std::endl;
  std::cout << "TriggerWeight_BToF =" << TriggerWeight_BToF << std::endl;

  Double_t x[TrackWeight->GetN()], y[TrackWeight->GetN()];
  for(int i=0; i<TrackWeight->GetN(); i++){
     TrackWeight->GetPoint(i, x[i], y[i]);
     cout<<"x="<<x[i]<<"; "<<"y="<<y[i]<<"."<<std::endl;
     TrackProperties val;
     val.x=x[i];
     val.y=y[i];
     val.errx_up=TrackWeight->GetErrorXhigh(i+1);
     val.errx_down=TrackWeight->GetErrorXlow(i+1);
     val.erry_up=TrackWeight->GetErrorYhigh(i+1);
     val.erry_down=TrackWeight->GetErrorYlow(i+1);
     TrackCorr.push_back(val);
  }

  GaussU1Corr[0]={1.8,19.4,1.7, 20.9,0.0,10.0};
  GaussU1Corr[1]={3.0,19.9,2.7,21.2,10.0,20.0};
  GaussU1Corr[2]={2.9, 20.6,2.4, 21.7,20.0, 30.0};
  GaussU1Corr[3]={2.5,20.5,2.0, 22.1, 30.0,50.0};
  GaussU1Corr[4]={1.8,23.3,1.1, 23.7,50.0, 1000.0};


  GaussU2Corr[0]={0.0, 19.3,0.4, 21.1, 0.0, 10.0};
  GaussU2Corr[1]={0.0, 19.6,0.0 ,21.0 , 10.0, 20.0};
  GaussU2Corr[2]={0.0, 20.0,0.0 ,21.3 ,20.0, 30.0};
  GaussU2Corr[3]={0.0, 20.5,0.0, 21.5, 30.0, 50.0};
  GaussU2Corr[4]={0.0,21.6 , 0.0, 22.1, 50.0, 1000.0};

}//FakeRateMiniAODGetRates



FakeRateMiniAODGetRates::~FakeRateMiniAODGetRates()
{
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void FakeRateMiniAODGetRates::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "\n<------------THIS IS A NEW EVENT------------>" << std::endl;
  NEvents_->Fill(0);

  //Get ak4Jets particle collection
  edm::Handle<edm::View<pat::Jet> > pJets;
  iEvent.getByToken(jetTag_, pJets);

  //Get RECO Muons particle collection
  edm::Handle<edm::View<pat::Muon> > pMuons;
  iEvent.getByToken(muonsTag_, pMuons);

  //Get CleanJets Tau particle collection
  edm::Handle<edm::View<pat::Tau> > pTaus;
  iEvent.getByToken(tauTag_, pTaus);

  //Old Jet collection for bTagging
  edm::Handle<edm::View<pat::Muon> > pMu12;
  iEvent.getByToken(mu12Tag_, pMu12);

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

//////////////////////////////
// BTaggin Scale factors
//////////////////////////////
  // BTaggin SF
  std::cout << "===> Loading the input .csv SF file..." << std::endl;

  std::string inputCSVfile = "/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/CSVv2_Moriond17_B_H.csv";
  BTagCalibration calib("csvv2", inputCSVfile);
  std::cout << "  BTagEntry::FLAV_C=" << BTagEntry::FLAV_C << "  BTagEntry::FLAV_B=" << BTagEntry::FLAV_B <<  std::endl;
  BTagCalibrationReader reader(BTagEntry::OP_MEDIUM, "central") ;
  
  reader.load(calib, BTagEntry::FLAV_B, "comb");
  reader.load(calib, BTagEntry::FLAV_C, "comb");
  reader.load(calib, BTagEntry::FLAV_UDSG, "incl");
  std::cout  << "initialized reader" << std::endl;

  double pileupWeight = 1, genWeight = 1, tauMedSF = 1.0, totalWeight = 1.0, SFsWeight = 1.0;
  if (isMC_)
  {
    double tauMedSF = 0.97;
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

  ///////////////////////////////////////
  // Do Track Iso and ID and ISO scale
  //////////////////////////////////////
    double IDs_weightBF=1.0;//every muon pass through "medium ID"
    double IDs_weightGH=1.0;
    double ISOs_weightBF=1.0; // every muon pass through 0.25 relative isolation
    double ISOs_weightGH=1.0;
    double Tracks_weight=1.0;
    float LumiFraction_GH = 16.1 / 35.9;
    float LumiFraction_BF = 19.8 / 35.9; //1.0-LumiFraction_GH;
  
    float binxIDs_BF = IDsWeight_BToF->GetXaxis()->FindBin(mu1.pt());
    float binyIDs_BF = IDsWeight_BToF->GetYaxis()->FindBin(fabs(mu1.eta()));
    float binxIDs_GH = IDsWeight_GH->GetXaxis()->FindBin(mu1.pt());
    float binyIDs_GH = IDsWeight_GH->GetYaxis()->FindBin(fabs(mu1.eta()));
    float binxISOs_BF = ISOsWeight_BToF->GetXaxis()->FindBin(mu1.pt());
    float binyISOs_BF = ISOsWeight_BToF->GetYaxis()->FindBin(fabs(mu1.eta()));
    float binxISOs_GH = ISOsWeight_GH->GetXaxis()->FindBin(mu1.pt());
    float binyISOs_GH = ISOsWeight_GH->GetYaxis()->FindBin(fabs(mu1.eta()));
    for (std::list<TrackProperties>::const_iterator it = TrackCorr.begin(); it != TrackCorr.end(); it++ )
    {
       if (fabs(mu1.eta()) >= (*it).x-(*it).errx_down && fabs(mu1.eta()) <= (*it).x+(*it).errx_up)
       {
          Tracks_weight *= (*it).y;
          break;
       }//if
    }//for it
    if (binxIDs_BF == 7)
    {
      binxIDs_BF = 6;
      binxIDs_GH = 6;
      binxISOs_BF = 6;
      binxISOs_GH = 6;
    }//if
  
    std::cout << "\n\nIDs_weightBF=" << IDs_weightBF << "\tIDsWeight_BToF->GetBinContent(binxIDs_BF, binyIDs_BF)=" << IDsWeight_BToF->GetBinContent(binxIDs_BF, binyIDs_BF) << std::endl;
    std::cout << "IDs_weightGH=" << IDs_weightGH << "\tIDsWeight_GH->GetBinContent(binxIDs_GH, binyIDs_GH)=" << IDsWeight_GH->GetBinContent(binxIDs_GH, binyIDs_GH) << std::endl;
    std::cout << "ISOs_weightBF=" << ISOs_weightBF << "\tISOsWeight_BToF->GetBinContent(binxISOs_BF, binyISOs_BF)=" << ISOsWeight_BToF->GetBinContent(binxISOs_BF, binyISOs_BF) << std::endl;
    std::cout << "ISOs_weightGH=" << ISOs_weightGH << "\tISOsWeight_GH->GetBinContent(binxISOs_GH, binyISOs_GH)=" << ISOsWeight_GH->GetBinContent(binxISOs_GH, binyISOs_GH) << std::endl;
    IDs_weightBF = IDs_weightBF * IDsWeight_BToF->GetBinContent(binxIDs_BF, binyIDs_BF);
    IDs_weightGH = IDs_weightGH * IDsWeight_BToF->GetBinContent(binxIDs_GH, binyIDs_GH);
    ISOs_weightBF = ISOs_weightBF * ISOsWeight_BToF->GetBinContent(binxISOs_BF, binyISOs_BF);
    ISOs_weightGH = ISOs_weightGH * ISOsWeight_GH->GetBinContent(binxISOs_GH, binyISOs_GH);
  
    if (mu2.pt() > 20.0)
    {
      std::cout << "Z Eff" << std::endl;
      binxIDs_BF = IDsWeight_BToF->GetXaxis()->FindBin(mu2.pt());
      binyIDs_BF = IDsWeight_BToF->GetYaxis()->FindBin(fabs(mu2.eta()));
      binxIDs_GH = IDsWeight_GH->GetXaxis()->FindBin(mu2.pt());
      binyIDs_GH = IDsWeight_GH->GetYaxis()->FindBin(fabs(mu2.eta()));
      binxISOs_BF = ISOsWeight_BToF->GetXaxis()->FindBin(mu2.pt());
      binyISOs_BF = ISOsWeight_BToF->GetYaxis()->FindBin(fabs(mu2.eta()));
      binxISOs_GH = ISOsWeight_GH->GetXaxis()->FindBin(mu2.pt());
      binyISOs_GH = ISOsWeight_GH->GetYaxis()->FindBin(fabs(mu2.eta()));
      if (binxIDs_BF == 7)
      {
        binxIDs_BF = 6;
        binxIDs_GH = 6;
        binxISOs_BF = 6;
        binxISOs_GH = 6;
      }//if
  
      std::cout << "IDs_weightBF=" << IDs_weightBF << "\tIDsWeight_BF->GetBinContent(binxIDs_BF, binyIDs_BF)=" << IDsWeight_BToF->GetBinContent(binxIDs_BF, binyIDs_BF) << std::endl;
      std::cout << "IDs_weightGH=" << IDs_weightGH << "\tIDsWeight_GH->GetBinContent(binxIDs_GH, binyIDs_GH)=" << IDsWeight_GH->GetBinContent(binxIDs_GH, binyIDs_GH) << std::endl;
      std::cout << "ISOs_weightBF=" << ISOs_weightBF << "\tISOsWeight_BF->GetBinContent(binxISOs_BF, binyISOs_BF)=" << ISOsWeight_BToF->GetBinContent(binxISOs_BF, binyISOs_BF) << std::endl;
      std::cout << "ISOs_weightGH=" << ISOs_weightGH << "\tISOsWeight_GH->GetBinContent(binxISOs_GH, binyISOs_GH)=" << ISOsWeight_GH->GetBinContent(binxISOs_GH, binyISOs_GH) << std::endl;
      IDs_weightBF = IDs_weightBF * IDsWeight_BToF->GetBinContent(binxIDs_BF, binyIDs_BF);
      IDs_weightGH = IDs_weightGH * IDsWeight_GH->GetBinContent(binxIDs_GH, binyIDs_GH);
      ISOs_weightBF = ISOs_weightBF * ISOsWeight_BToF->GetBinContent(binxISOs_BF, binyISOs_BF);
      ISOs_weightGH = ISOs_weightGH * ISOsWeight_GH->GetBinContent(binxISOs_GH, binyISOs_GH);
    }//if
    else
    {
      std::cout << "JPsi Eff" << std::endl;
  
      binyIDs_BF = IDsWeight_LowPt_ID->GetXaxis()->FindBin(mu2.pt());
      binxIDs_BF = IDsWeight_LowPt_ID->GetYaxis()->FindBin(fabs(mu2.eta()));
      binyIDs_GH = IDsWeight_LowPt_ID->GetXaxis()->FindBin(mu2.pt());
      binxIDs_GH = IDsWeight_LowPt_ID->GetYaxis()->FindBin(fabs(mu2.eta()));
      binyISOs_BF = ISOsWeight_LowPt_ISO->GetXaxis()->FindBin(mu2.pt());
      binxISOs_BF = ISOsWeight_LowPt_ISO->GetYaxis()->FindBin(fabs(mu2.eta()));
      binyISOs_GH = ISOsWeight_LowPt_ISO->GetXaxis()->FindBin(mu2.pt());
      binxISOs_GH = ISOsWeight_LowPt_ISO->GetYaxis()->FindBin(fabs(mu2.eta()));
      if (binyIDs_BF == 7)
      {
        binyIDs_BF = 6;
        binyIDs_GH = 6;
        binyISOs_BF = 6;
        binyISOs_GH = 6;
      }//if
      std::cout << "IDs_weightBF=" << IDs_weightBF << "\tIDsWeight_LowPt_ID->GetBinContent(binxIDs_BF, binyIDs_BF)=" << IDsWeight_LowPt_ID->GetBinContent(binxIDs_BF, binyIDs_BF) << std::endl;
      std::cout << "ISOs_weightGH=" << ISOs_weightGH << "\tISOsWeight_LowPt_ISO->GetBinContent(binxISOs_GH, binyISOs_GH)=" << ISOsWeight_LowPt_ISO->GetBinContent(binxISOs_GH, binyISOs_GH) << std::endl;
      IDs_weightBF = IDs_weightBF * IDsWeight_LowPt_ID->GetBinContent(binxIDs_BF, binyIDs_BF);
      IDs_weightGH = IDs_weightGH * IDsWeight_LowPt_ID->GetBinContent(binxIDs_GH, binyIDs_GH);
      ISOs_weightBF = ISOs_weightBF * ISOsWeight_LowPt_ISO->GetBinContent(binxISOs_BF, binyISOs_BF);
      ISOs_weightGH = ISOs_weightGH * ISOsWeight_LowPt_ISO->GetBinContent(binxISOs_GH, binyISOs_GH);
    }//else
  
    for (std::list<TrackProperties>::const_iterator it = TrackCorr.begin(); it != TrackCorr.end(); it++ )
    {
       if (fabs(mu2.eta()) >= (*it).x-(*it).errx_down && fabs(mu2.eta()) <= (*it).x+(*it).errx_up)
       {
          Tracks_weight *= (*it).y;
          break;
       }//if 
    }//for it
  
  
  ///////////////////////////////////////
  // Trigger efficiencies
  ///////////////////////////////////////
    double Trigger_weightBF=1.0;
    double Trigger_weightGH=1.0;
    float binxTrigger_BF = TriggerWeight_BToF->GetXaxis()->FindBin(mu1.pt());
    float binyTrigger_BF = TriggerWeight_BToF->GetYaxis()->FindBin(fabs(mu1.eta()));
    float binxTrigger_GH = TriggerWeight_GH->GetXaxis()->FindBin(mu1.pt());
    float binyTrigger_GH = TriggerWeight_GH->GetYaxis()->FindBin(fabs(mu1.eta()));
    Trigger_weightBF = Trigger_weightBF *  TriggerWeight_BToF->GetBinContent(binxTrigger_BF, binyTrigger_BF);
    Trigger_weightGH = Trigger_weightGH *  TriggerWeight_GH->GetBinContent(binxTrigger_GH, binyTrigger_GH);
  
  ///////////////////////////////////////
  // Combining weights and SF
  ///////////////////////////////////////
    std::cout << "\n\nIDs_weightBF=" << IDs_weightBF << "\tISOs_weightGH=" << ISOs_weightGH << "\tTracks_weight=" << Tracks_weight << "\tTrigger_weightGH=" << Trigger_weightGH << std::endl;
    SFsWeight = LumiFraction_BF * (IDs_weightBF * ISOs_weightBF * Trigger_weightBF) +
                       LumiFraction_GH * (IDs_weightGH * ISOs_weightGH * Trigger_weightGH);
    totalWeight = SFsWeight*pileupWeight*genWeight*tauMedSF;
    if (!isMC_)
      SFsWeight = 1.0;
    std::cout << "pileupWeight=" << pileupWeight << "\tgenWeight=" << genWeight << "\tSFsWeight=" << SFsWeight << std::endl;
  }//if isMC


  reco::LeafCandidate::LorentzVector diMuP4 = mu1.p4() + mu2.p4();
std::cout << "mu1: pt=" << mu1.pt() << "  eta=" << mu1.eta() << "  pdgID=" << mu1.pdgId() << "\nmu2: pt=" << mu2.pt() << "  eta=" << mu2.eta() << "  pdgID=" << mu2.pdgId() << std::endl;
std::cout << "DiMu: pt= " << diMuP4.Pt() << "  eta=" << diMuP4.Eta() << " Mass=" << diMuP4.M() << std::endl;

  if (checkInvMass_ && (diMuP4.M() < checkInvMassMin_ || diMuP4.M() > checkInvMassMax_) )
    return;
  InvMassMu1Mu2_->Fill(diMuP4.M() );


//////////////////////////////
// Begin Analyzer
//////////////////////////////
  ///////////////////////////
  // iterating over the taus
  /////////////////////////// 
  reco::LeafCandidate::LorentzVector  diTauP4;
  double tauDecayMode = -1, DMFind = -1, LooseIso = -1, MedIso = -1, TightIso = -1, bTagValue = -1000, smallestdRJet = 1000000000000, csvSF = 1.0;
  //double VLooseIso = -1, VTightIso = -1;
  for (edm::View<pat::Tau>::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); ++iTau)
  {
    NEvents_->Fill(1);
    double bestMu3dR = 10000000;
    bool checkSubMu = false;
    pat::Muon mu3;
    if (iTau->tauID(medIsoTag_) <= -0.5 || iTau->pt() <= tauPtCut_ || fabs(iTau->eta() ) >= 2.4 || deltaR(mu1,*iTau) <= 0.8 || deltaR(mu2,*iTau) <= 0.8)
      continue;
    for (edm::View<pat::Muon>::const_iterator iMu = pMuons->begin(); iMu != pMuons->end(); ++iMu)
    {
      double currdR = deltaR(*iTau, *iMu), mu1dR = deltaR(mu1, *iMu), mu2dR = deltaR(mu2, *iMu);
      if (mu1dR < mu3dROverlapCut_ || mu2dR < mu3dROverlapCut_)
       continue;
      if (currdR < mu3dRMax_ && currdR > mu3dRMin_ && currdR < bestMu3dR)
      {
        for (edm::View<pat::Jet>::const_iterator iJet = pJets->begin(); iJet != pJets->end(); ++iJet)
        {
          double dRJet = deltaR(*iTau, *iJet);
          if (dRJet < smallestdRJet)
          {
            bTagValue = iJet->bDiscriminator(csvBTag_);
            smallestdRJet = dRJet;
            if (isMC_)
            {
              int hadronFlavor = iJet->hadronFlavour();
              if ( hadronFlavor == 5 )
                csvSF = reader.eval(BTagEntry::FLAV_B, fabs(iJet->eta()), iJet->pt(), bTagValue);
  
              else if( hadronFlavor == 4 )
                csvSF = reader.eval(BTagEntry::FLAV_C, fabs(iJet->eta()), iJet->pt(), bTagValue);
  
              else
                csvSF = reader.eval(BTagEntry::FLAV_UDSG, fabs(iJet->eta()), iJet->pt(), bTagValue);
              std::cout << "csvSF=" << csvSF << "  hadronFlavor=" << hadronFlavor << " bTagValue=" << bTagValue << std::endl;
            }//if isMC_
          }//if
        }//for iJet
        diTauP4 = iTau->p4() + iMu->p4();
        checkSubMu = true;
        mu3 = *iMu;
        tauDecayMode = iTau->decayMode();
        DMFind = iTau->tauID(decayModeFindingTag_);
        LooseIso = iTau->tauID(looseIsoTag_);
        MedIso = iTau->tauID(medIsoTag_);
        TightIso = iTau->tauID(tightIsoTag_);
      }//if
    }//for iMu

    if (checkSubMu)
      NEvents_->Fill(2);
    else if (iTau->pt() > tauPtCut_ )
      NEvents_->Fill(3);
    else if (DMFind >= .5)
      NEvents_->Fill(4);
    else if (LooseIso >= .5)
      NEvents_->Fill(5);
    else if (MedIso >= .5)
      NEvents_->Fill(6);
    else if (TightIso >= .5)
      NEvents_->Fill(7);

    if (!checkSubMu && requireRemovedMuon_) 
      continue;

    double csvWeight = 1.0;
    if (!isMC_ && bTagValue >= 0.8484)
      continue;
    else if (isMC_ && bTagValue >= 0.8484)
      csvWeight = 1 - csvSF;
    totalWeight = totalWeight*csvWeight;

    std::cout << "csvWeight=" << csvWeight << "  totalWeight=" << totalWeight << "  genWeight=" << genWeight;
    std::cout << "  pileupWeight=" << pileupWeight << " SFsWeight=" << SFsWeight << "  tauMedSF=" <<  tauMedSF << std::endl;

    if (checkSubMu)
    {
      reco::LeafCandidate::LorentzVector diMuP4_Mu1TauMu, diMuP4_Mu2TauMu;    
      diMuP4_Mu1TauMu = mu1.p4();
      diMuP4_Mu1TauMu += mu3.p4();
 
      diMuP4_Mu2TauMu = mu2.p4();
      diMuP4_Mu2TauMu += mu3.p4();

      InvMassTauMuMu1_->Fill(diMuP4_Mu1TauMu.M() ,totalWeight );
      InvMassTauMuMu2_->Fill(diMuP4_Mu2TauMu.M() ,totalWeight );
    }//if removed Mu

    if (DMFind >= .5)
    {
      TauDMFindPt_->Fill(iTau->pt() ,totalWeight );
      TauDMFindEta_->Fill(iTau->eta() ,totalWeight );
      EtavsPtTauDMFind_->Fill(iTau->pt(), iTau->eta() ,totalWeight );
      if (tauDecayMode == 0)
      {
        OneProngDMEta_->Fill(iTau->eta() ,totalWeight );
        OneProngDMPt_->Fill(iTau->pt() ,totalWeight );  
      }//else if tauDecayMode >= .5
      else if (tauDecayMode == 1)
      {
        OneProngOnePizDMEta_->Fill(iTau->eta() ,totalWeight );
        OneProngOnePizDMPt_->Fill(iTau->pt() ,totalWeight );  
      }//else if tauDecayMode >= .5
      else if (tauDecayMode == 2)
      {
        OneProngTwoPizDMEta_->Fill(iTau->eta() ,totalWeight );
        OneProngTwoPizDMPt_->Fill(iTau->pt() ,totalWeight );  
      }//else if tauDecayMode >= .5
      else if (tauDecayMode == 10)
      {
        ThreeProngDMEta_->Fill(iTau->eta() ,totalWeight );
        ThreeProngDMPt_->Fill(iTau->pt() ,totalWeight );  
      }//else if tauDecayMode >= .5
    }//if DMFind >= .5

    if (TightIso >= .5 && DMFind >= .5)
    {
      TauTightIsoPt_->Fill(iTau->pt() ,totalWeight );
      TauTightIsoEta_->Fill(iTau->eta() ,totalWeight );
      EtavsPtTauTightIso_->Fill(iTau->pt(), iTau->eta() ,totalWeight );
    }//if TightIso >= .5 &&DMFind >= .5

    if (MedIso >= .5 && DMFind >= .5)
    {
      TauMedIsoPt_->Fill(iTau->pt() ,totalWeight );
      TauMedIsoEta_->Fill(iTau->eta() ,totalWeight );
      EtavsPtTauMedIso_->Fill(iTau->pt(), iTau->eta() ,totalWeight );
      if (tauDecayMode == 0)
      {
        OneProngMedIsoEta_->Fill(iTau->eta() ,totalWeight );
        OneProngMedIsoPt_->Fill(iTau->pt() ,totalWeight );  
      }//else if tauDecayMode >= .5
      else if (tauDecayMode == 1)
      {
        OneProngOnePizMedIsoEta_->Fill(iTau->eta() ,totalWeight );
        OneProngOnePizMedIsoPt_->Fill(iTau->pt() ,totalWeight );
      }//else if tauDecayMode >= .5
      else if (tauDecayMode == 2)
      {
        OneProngTwoPizMedIsoEta_->Fill(iTau->eta() ,totalWeight );
        OneProngTwoPizMedIsoPt_->Fill(iTau->pt() ,totalWeight );
      }//else if tauDecayMode >= .5
      else if (tauDecayMode == 10)
      {
        ThreeProngMedIsoEta_->Fill(iTau->eta() ,totalWeight );
        ThreeProngMedIsoPt_->Fill(iTau->pt() ,totalWeight );
      }//else if tauDecayMode == 1
    }//if MedIso == 1 && DMFind >= .5

    if (LooseIso >= .5 && DMFind >= .5)
    {
      TauLooseIsoPt_->Fill(iTau->pt() ,totalWeight );
      TauLooseIsoEta_->Fill(iTau->eta() ,totalWeight );
      EtavsPtTauLooseIso_->Fill(iTau->pt(), iTau->eta() ,totalWeight );
    }//if Loose DMFind == 1
  }//iTau

  for (edm::View<pat::Jet>::const_iterator iJet = pJets->begin(); iJet != pJets->end(); ++iJet)
  {
    if (iJet->pt() > tauPtCut_ && fabs(iJet->eta() ) < 2.4 )
    {
      JetEta_->Fill(iJet->eta() ,totalWeight ); 
      JetPt_->Fill(iJet->pt() ,totalWeight ); 
      EtavsPtJet_->Fill(iJet->pt(), iJet->eta() ,totalWeight );
      for (edm::View<pat::Muon>::const_iterator iMu = pMuons->begin(); iMu != pMuons->end(); ++iMu)
       {  
       double currdR = deltaR(*iJet, *iMu), mu1dR = deltaR(mu1, *iMu), mu2dR = deltaR(mu2, *iMu);
       if (mu1dR < .5 || mu2dR < .5)
         continue;
       if (currdR < mu3dRMax_ && currdR > mu3dRMin_)
        {
          JetEtaWithSoftMuon_->Fill(iJet->eta() ,totalWeight );
          JetPtWithSoftMuon_->Fill(iJet->pt() ,totalWeight );
          EtavsPtJetSoftMuon_->Fill(iJet->pt(), iJet->eta() ,totalWeight );          
          reco::LeafCandidate::LorentzVector jetP4 = iJet->p4() - iMu->p4();
          JetEtaWithSoftMuon_noMu_->Fill(jetP4.Eta() ,totalWeight );
          JetPtWithSoftMuon_noMu_->Fill(jetP4.Pt() ,totalWeight );
          EtavsPtJetSoftMuon_noMu_->Fill(jetP4.Pt(), jetP4.Eta() ,totalWeight );
        }//if iSoftMuon
      }//for i
    }//if 
  }//iJEt
}//End FakeRateMiniAODGetRates::analyze


// ------------ method called once each job just before starting event loop  ------------
void FakeRateMiniAODGetRates::beginJob()
{
  std::cout << "Begin Job" << std::endl;

  //Open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");


  //Book histograms
  NEvents_     = new TH1F("NEvents"    , "", 9, -.5, 8.5);
      NEvents_->GetXaxis()->SetBinLabel(1, "TotalEvents"); 
      NEvents_->GetXaxis()->SetBinLabel(2, "Total #tau's");
      NEvents_->GetXaxis()->SetBinLabel(3, "#tau's With #mu");
      NEvents_->GetXaxis()->SetBinLabel(4, "#tau's p_{T} > 20");
      NEvents_->GetXaxis()->SetBinLabel(5, "Pass DMFind");
      NEvents_->GetXaxis()->SetBinLabel(6, "Pass Loose Iso");
      NEvents_->GetXaxis()->SetBinLabel(7, "Pass Med Iso");
      NEvents_->GetXaxis()->SetBinLabel(8, "Pass Tight Iso");
  InvMassTauMuMu1_     = new TH1F("InvMassTauMuMu1"    , "", 75, 0, 150);
  InvMassMu1Mu2_     = new TH1F("InvMassMu1Mu2"    , "", 75, 0, 150);
  InvMassTauMuMu2_     = new TH1F("InvMassTauMuMu2"    , "", 75, 0, 150);

  Float_t binsx[] = {9.9999, 20, 30, 40, 60, 300};
  Float_t binsy[] = {0, .9, 1.5, 2.4};
  EtavsPtTauLooseIso_  = new TH2F("EtavsPtTauLooseIso" , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx, sizeof(binsy)/sizeof(Float_t) - 1, binsy);
  EtavsPtTauMedIso_  = new TH2F("EtavsPtTauMedIso"     , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx, sizeof(binsy)/sizeof(Float_t) - 1, binsy);
  EtavsPtTauTightIso_  = new TH2F("EtavsPtTauTightIso" , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx, sizeof(binsy)/sizeof(Float_t) - 1, binsy);
  EtavsPtTauDMFind_  = new TH2F("EtavsPtTauDMFind"     , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx, sizeof(binsy)/sizeof(Float_t) - 1, binsy);
  EtavsPtJet_  = new TH2F("EtavsPtJet"                 , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx, sizeof(binsy)/sizeof(Float_t) - 1, binsy);
  EtavsPtJetSoftMuon_  = new TH2F("EtavsPtJetSoftMuon" , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx, sizeof(binsy)/sizeof(Float_t) - 1, binsy);
  EtavsPtJetSoftMuon_noMu_  = new TH2F("EtavsPtJetSoftMuon_noMu" , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx, sizeof(binsy)/sizeof(Float_t) - 1, binsy);

  TauLooseIsoEta_  = new TH1F("TauLooseIsoEta"    , "", 11, -2.4, 2.4);
  TauMedIsoEta_    = new TH1F("TauMedIsoEta", "", 11, -2.4, 2.4);
  TauTightIsoEta_    = new TH1F("TauTightIsoEta", "", 11, -2.4, 2.4);
  TauDMFindEta_    = new TH1F("TauDMFindEta"    , "", 11, -2.4, 2.4);
  JetEta_          = new TH1F("JetEta"    , "", 11, -2.4, 2.4);
  JetEtaWithSoftMuon_          = new TH1F("JetEtaWithSoftMuon"    , "", 11, -2.4, 2.4);
  JetEtaWithSoftMuon_noMu_          = new TH1F("JetEtaWithSoftMuon_noMu"    , "", 11, -2.4, 2.4);

  TauLooseIsoPt_  = new TH1F("TauLooseIsoPt"    , "", 50, 20, 220.0);
  TauMedIsoPt_    = new TH1F("TauMedIsoPt", "", sizeof(binsx)/sizeof(Float_t) - 1, binsx);
  TauTightIsoPt_    = new TH1F("TauTightIsoPt", "", 50, 20, 220.0);
  TauDMFindPt_    = new TH1F("TauDMFindPt"    , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx);
  JetPt_          = new TH1F("JetPt"    , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx);
  JetPtWithSoftMuon_          = new TH1F("JetPtWithSoftMuon"    , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx);
  JetPtWithSoftMuon_noMu_          = new TH1F("JetPtWithSoftMuon_noMu"    , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx);

  FinalEffLooseIsoEta_ = new TGraphAsymmErrors(11);
  FinalEffMedIsoEta_ = new TGraphAsymmErrors(11);
  FinalEffTightIsoEta_ = new TGraphAsymmErrors(11);
  FinalEffDMFindEta_ = new TGraphAsymmErrors(11);

  FinalEffLooseIsoPt_ = new TGraphAsymmErrors(11);
  FinalEffMedIsoPt_ = new TGraphAsymmErrors(11);
  FinalEffTightIsoPt_ = new TGraphAsymmErrors(11);
  FinalEffDMFindPt_ = new TGraphAsymmErrors(11);

  OneProngDMEta_ = new TH1F("OneProngDMEta"    , "", 11, -2.4, 2.4);
  OneProngOnePizDMEta_ = new TH1F("OneProngOnePizDMEta"    , "", 11, -2.4, 2.4);
  OneProngTwoPizDMEta_ = new TH1F("OneProngTwoPizDMEta"    , "", 11, -2.4, 2.4);
  ThreeProngDMEta_ = new TH1F("ThreeProngDMEta"    , "", 11, -2.4, 2.4);
  OneProngDMPt_ = new TH1F("OneProngDMPt"    , "", 11, 20, 220.0);
  OneProngOnePizDMPt_ = new TH1F("OneProngOnePizDMPt"    , "", 11, 20, 220.0);
  OneProngTwoPizDMPt_ = new TH1F("OneProngTwoPizDMPt"    , "", 11, 20, 220.0);
  ThreeProngDMPt_ = new TH1F("ThreeProngDMPt"    , "", 11, 20, 220.0);
  FinalOneProngDMEta_ = new TGraphAsymmErrors(11);
  FinalOneProngOnePizDMEta_ = new TGraphAsymmErrors(11);
  FinalOneProngTwoPizDMEta_ = new TGraphAsymmErrors(11);
  FinalThreeProngDMEta_ = new TGraphAsymmErrors(11);
  FinalOneProngDMPt_ = new TGraphAsymmErrors(11);
  FinalOneProngOnePizDMPt_ = new TGraphAsymmErrors(11);
  FinalOneProngTwoPizDMPt_ = new TGraphAsymmErrors(11);
  FinalThreeProngDMPt_ = new TGraphAsymmErrors(11);

  OneProngMedIsoEta_ = new TH1F("OneProngMedIsoEta"    , "", 11, -2.4, 2.4);
  OneProngOnePizMedIsoEta_ = new TH1F("OneProngOnePizMedIsoEta"    , "", 11, -2.4, 2.4);
  OneProngTwoPizMedIsoEta_ = new TH1F("OneProngTwoPizMedIsoEta"    , "", 11, -2.4, 2.4);
  ThreeProngMedIsoEta_ = new TH1F("ThreeProngMedIsoEta"    , "", 11, -2.4, 2.4);
  OneProngMedIsoPt_ = new TH1F("OneProngMedIsoPt"    , "", 11, 20, 220.0);
  OneProngOnePizMedIsoPt_ = new TH1F("OneProngOnePizMedIsoPt"    , "", 11, 20, 220.0);
  OneProngTwoPizMedIsoPt_ = new TH1F("OneProngTwoPizMedIsoPt"    , "", 11, 20, 220.0);
  ThreeProngMedIsoPt_ = new TH1F("ThreeProngMedIsoPt"    , "", 11, 20, 220.0);
  FinalOneProngMedIsoEta_ = new TGraphAsymmErrors(11);
  FinalOneProngOnePizMedIsoEta_ = new TGraphAsymmErrors(11);
  FinalOneProngTwoPizMedIsoEta_ = new TGraphAsymmErrors(11);
  FinalThreeProngMedIsoEta_ = new TGraphAsymmErrors(11);
  FinalOneProngMedIsoPt_ = new TGraphAsymmErrors(11);
  FinalOneProngOnePizMedIsoPt_ = new TGraphAsymmErrors(11);
  FinalOneProngTwoPizMedIsoPt_ = new TGraphAsymmErrors(11);
  FinalThreeProngMedIsoPt_ = new TGraphAsymmErrors(11);

  InvMassTauMuMu1_->Sumw2();
  InvMassMu1Mu2_->Sumw2();
  InvMassTauMuMu2_->Sumw2();

  EtavsPtTauLooseIso_->Sumw2();
  EtavsPtTauMedIso_->Sumw2();
  EtavsPtTauTightIso_->Sumw2();
  EtavsPtTauDMFind_->Sumw2();
  EtavsPtJet_->Sumw2();
  EtavsPtJetSoftMuon_->Sumw2();
  EtavsPtJetSoftMuon_noMu_->Sumw2();

  TauLooseIsoEta_->Sumw2();
  TauMedIsoEta_->Sumw2();
  TauTightIsoEta_->Sumw2();
  TauDMFindEta_->Sumw2();
  JetEta_->Sumw2();
  JetEtaWithSoftMuon_->Sumw2();
  JetEtaWithSoftMuon_noMu_->Sumw2();

  TauLooseIsoPt_->Sumw2();
  TauMedIsoPt_->Sumw2();
  TauTightIsoPt_->Sumw2();
  TauDMFindPt_->Sumw2();
  JetPt_->Sumw2();
  JetPtWithSoftMuon_->Sumw2();
  JetPtWithSoftMuon_noMu_->Sumw2();

  OneProngDMEta_->Sumw2();
  OneProngOnePizDMEta_->Sumw2();
  OneProngTwoPizDMEta_->Sumw2();
  ThreeProngDMEta_->Sumw2();
  OneProngDMPt_->Sumw2();
  OneProngOnePizDMPt_->Sumw2();
  OneProngTwoPizDMPt_->Sumw2();
  ThreeProngDMPt_->Sumw2();

  OneProngMedIsoEta_->Sumw2();
  OneProngOnePizMedIsoEta_->Sumw2();
  OneProngTwoPizMedIsoEta_->Sumw2();
  ThreeProngMedIsoEta_->Sumw2();
  OneProngMedIsoPt_->Sumw2();
  OneProngOnePizMedIsoPt_->Sumw2();
  OneProngTwoPizMedIsoPt_->Sumw2();
  ThreeProngMedIsoPt_->Sumw2();
}

// ------------ method called once each job just after ending the event loop  ------------
void FakeRateMiniAODGetRates::endJob()
{
  //Make the Canvases
  TCanvas NEventsCanvas("NEvents","",600,600);
  TCanvas InvMassTauMuMu1Canvas("InvMassTauMuMu1","",600,600);
  TCanvas InvMassMu1Mu2Canvas("InvMassMu1Mu2","",600,600);
  TCanvas InvMassTauMuMu2Canvas("InvMassTauMuMu2","",600,600);

  TCanvas EtavsPtTauLooseIsoCanvas("EtavsPtTauLooseIso","",600,600);
  TCanvas EtavsPtTauMedIsoCanvas("EtavsPtTauMedIso","",600,600);
  TCanvas EtavsPtTauTightIsoCanvas("EtavsPtTauTightIso","",600,600);
  TCanvas EtavsPtTauDMFindCanvas("EtavsPtTauDMFind","",600,600);
  TCanvas EtavsPtJetCanvas("EtavsPtJet","",600,600);
  TCanvas EtavsPtJetSoftMuonCanvas("EtavsPtJetSoftMuon","",600,600);
  TCanvas EtavsPtJetSoftMuon_noMuCanvas("EtavsPtJetSoftMuon_noMu","",600,600);

  TCanvas TauLooseIsoEtaCanvas("TauLooseIsoEta","",600,600);
  TCanvas TauMedIsoEtaCanvas("TauMedIsoEta","",600,600);
  TCanvas TauTightIsoEtaCanvas("TauTightIsoEta","",600,600);
  TCanvas TauDMFindEtaCanvas("TauDMFindEta","",600,600);
  TCanvas JetEtaCanvas("JetEta","",600,600);
  TCanvas JetEtaWithSoftMuonCanvas("JetEtaWithSoftMuon","",600,600);
  TCanvas JetEtaWithSoftMuon_noMuCanvas("JetEtaWithSoftMuon_noMu","",600,600);
  TCanvas FinalEffLooseIsoEtaCanvas("FinalEffLooseIsoEta","",600,600);
  TCanvas FinalEffMedIsoEtaCanvas("FinalEffMedIsoEta","",600,600);
  TCanvas FinalEffTightIsoEtaCanvas("FinalEffTightIsoEta","",600,600);
  TCanvas FinalEffDMFindEtaCanvas("FinalEffDMFindEta","",600,600);

  TCanvas TauLooseIsoPtCanvas("TauLooseIsoPt","",600,600);
  TCanvas TauMedIsoPtCanvas("TauMedIsoPt","",600,600);
  TCanvas TauTightIsoPtCanvas("TauTightIsoPt","",600,600);
  TCanvas TauDMFindPtCanvas("TauDMFindPt","",600,600);
  TCanvas JetPtCanvas("JetPt","",600,600);
  TCanvas JetPtWithSoftMuonCanvas("JetPtWithSoftMuon","",600,600);
  TCanvas JetPtWithSoftMuon_noMuCanvas("JetPtWithSoftMuon_noMu","",600,600);
  TCanvas FinalEffLooseIsoPtCanvas("FinalEffLooseIsoPt","",600,600);
  TCanvas FinalEffMedIsoPtCanvas("FinalEffMedIsoPt","",600,600);
  TCanvas FinalEffTightIsoPtCanvas("FinalEffTightIsoPt","",600,600);
  TCanvas FinalEffDMFindPtCanvas("FinalEffDMFindPt","",600,600);

  TCanvas OneProngDMEtaCanvas("OneProngDMEta","",600,600);
  TCanvas OneProngOnePizDMEtaCanvas("OneProngOnePizDMEta","",600,600);
  TCanvas OneProngTwoPizDMEtaCanvas("OneProngTwoPizDMEta","",600,600);
  TCanvas ThreeProngDMEtaCanvas("ThreeProngDMEta","",600,600);
  TCanvas OneProngDMPtCanvas("OneProngDMPt","",600,600);
  TCanvas OneProngOnePizDMPtCanvas("OneProngOnePizDMPt","",600,600);
  TCanvas OneProngTwoPizDMPtCanvas("OneProngTwoPizDMPt","",600,600);
  TCanvas ThreeProngDMPtCanvas("ThreeProngDMPt","",600,600);
  TCanvas FinalOneProngDMEtaCanvas("FinalOneProngDMEta","",600,600);
  TCanvas FinalOneProngOnePizDMEtaCanvas("FinalOneProngOnePizDMEta","",600,600);
  TCanvas FinalOneProngTwoPizDMEtaCanvas("FinalOneProngTwoPizDMEta","",600,600);
  TCanvas FinalThreeProngDMEtaCanvas("FinalThreeProngDMEta","",600,600);
  TCanvas FinalOneProngDMPtCanvas("FinalOneProngDMPt","",600,600);
  TCanvas FinalOneProngOnePizDMPtCanvas("FinalOneProngOnePizDMPt","",600,600);
  TCanvas FinalOneProngTwoPizDMPtCanvas("FinalOneProngTwoPizDMPt","",600,600);
  TCanvas FinalThreeProngDMPtCanvas("FinalThreeProngDMPt","",600,600);

  TCanvas OneProngMedIsoEtaCanvas("OneProngMedIsoEta","",600,600);
  TCanvas OneProngOnePizMedIsoEtaCanvas("OneProngOnePizMedIsoEta","",600,600);
  TCanvas OneProngTwoPizMedIsoEtaCanvas("OneProngTwoPizMedIsoEta","",600,600);
  TCanvas ThreeProngMedIsoEtaCanvas("ThreeProngMedIsoEta","",600,600);
  TCanvas OneProngMedIsoPtCanvas("OneProngMedIsoPt","",600,600);
  TCanvas OneProngOnePizMedIsoPtCanvas("OneProngOnePizMedIsoPt","",600,600);
  TCanvas OneProngTwoPizMedIsoPtCanvas("OneProngTwoPizMedIsoPt","",600,600);
  TCanvas ThreeProngMedIsoPtCanvas("ThreeProngMedIsoPt","",600,600);
  TCanvas FinalOneProngMedIsoEtaCanvas("FinalOneProngMedIsoEta","",600,600);
  TCanvas FinalOneProngOnePizMedIsoEtaCanvas("FinalOneProngOnePizMedIsoEta","",600,600);
  TCanvas FinalOneProngTwoPizMedIsoEtaCanvas("FinalOneProngTwoPizMedIsoEta","",600,600);
  TCanvas FinalThreeProngMedIsoEtaCanvas("FinalThreeProngMedIsoEta","",600,600);
  TCanvas FinalOneProngMedIsoPtCanvas("FinalOneProngMedIsoPt","",600,600);
  TCanvas FinalOneProngOnePizMedIsoPtCanvas("FinalOneProngOnePizMedIsoPt","",600,600);
  TCanvas FinalOneProngTwoPizMedIsoPtCanvas("FinalOneProngTwoPizMedIsoPt","",600,600);
  TCanvas FinalThreeProngMedIsoPtCanvas("FinalThreeProngMedIsoPt","",600,600);

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


  VariousFunctions::formatAndDrawCanvasAndHist2D(EtavsPtTauLooseIsoCanvas, EtavsPtTauLooseIso_,
         1, 0, 0, kBlack, 7, 20, "p_{T} #tau Loose Iso", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(EtavsPtTauMedIsoCanvas, EtavsPtTauMedIso_,
         1, 0, 0, kBlack, 7, 20, "p_{T} #tau Med Iso", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(EtavsPtTauTightIsoCanvas, EtavsPtTauTightIso_,
         1, 0, 0, kBlack, 7, 20, "p_{T} #tau Tight Iso", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(EtavsPtTauDMFindCanvas, EtavsPtTauDMFind_,
         1, 0, 0, kBlack, 7, 20, "p_{T} #tau DMFind", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(EtavsPtJetCanvas, EtavsPtJet_,
         1, 0, 0, kBlack, 7, 20, "p_{T} Jet", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(EtavsPtJetSoftMuonCanvas, EtavsPtJetSoftMuon_,
         1, 0, 0, kBlack, 7, 20, "p_{T}(Jet with Soft Mu)", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(EtavsPtJetSoftMuon_noMuCanvas, EtavsPtJetSoftMuon_noMu_,
         1, 0, 0, kBlack, 7, 20, "p_{T}(Jet with Soft Mu Removed)", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);

  VariousFunctions::formatAndDrawCanvasAndHist1D(TauLooseIsoEtaCanvas, TauLooseIsoEta_,
	 1, 0, 0, kBlack, 7, 20, "Eta(PFTau + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauMedIsoEtaCanvas, TauMedIsoEta_,
	 1, 0, 0, kBlack, 7, 20, "Eta(PFTau + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauTightIsoEtaCanvas, TauTightIsoEta_,
	 1, 0, 0, kBlack, 7, 20, "Eta(PFTau + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauDMFindEtaCanvas, TauDMFindEta_,
	 1, 0, 0, kBlack, 7, 20, "Eta(PFTau + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(JetEtaCanvas, JetEta_,
	 1, 0, 0, kBlack, 7, 20, "Eta(ak4Jets)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(JetEtaWithSoftMuonCanvas, JetEtaWithSoftMuon_,
	 1, 0, 0, kBlack, 7, 20, "Eta(ak4Jets) with Soft Muon", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(JetEtaWithSoftMuon_noMuCanvas, JetEtaWithSoftMuon_noMu_,
	 1, 0, 0, kBlack, 7, 20, "Eta(ak4Jets noMu) with Soft Muon", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoEtaCanvas, FinalEffLooseIsoEta_,
	 1, 1, 1, kBlack, 1, 20, "Eta(PFTau + Loose Iso + DM / nJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoEtaCanvas, FinalEffMedIsoEta_,
	 1, 1, 1, kBlack, 1, 20, "Eta(PFTau + Med Iso + DM / nJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoEtaCanvas, FinalEffTightIsoEta_,
	 1, 1, 1, kBlack, 1, 20, "Eta(PFTau + Tight Iso + DM / nJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffDMFindEtaCanvas, FinalEffDMFindEta_,
	 1, 1, 1, kBlack, 1, 20, "Eta(PFTau + DecayModeFinding / nJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);


  VariousFunctions::formatAndDrawCanvasAndHist1D(TauLooseIsoPtCanvas, TauLooseIsoPt_,
         1, 0, 0, kBlack, 7, 20, "Pt(PFTau + Loose Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauMedIsoPtCanvas, TauMedIsoPt_,
         1, 0, 0, kBlack, 7, 20, "Pt(PFTau + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauTightIsoPtCanvas, TauTightIsoPt_,
         1, 0, 0, kBlack, 7, 20, "Pt(PFTau + Tight Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauDMFindPtCanvas, TauDMFindPt_,
         1, 0, 0, kBlack, 7, 20, "Pt(PFTau + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(JetPtCanvas, JetPt_,
         1, 0, 0, kBlack, 7, 20, "Pt(ak4Jets)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(JetPtWithSoftMuonCanvas, JetPtWithSoftMuon_,
         1, 0, 0, kBlack, 7, 20, "Pt(ak4Jets) with Soft Muon", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(JetPtWithSoftMuon_noMuCanvas, JetPtWithSoftMuon_noMu_,
         1, 0, 0, kBlack, 7, 20, "Pt(ak4Jets noMu) with Soft Muon", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffLooseIsoPtCanvas, FinalEffLooseIsoPt_,
         1, 1, 1, kBlack, 1, 20, "Pt(PFTau + Loose Iso + DM / nJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffMedIsoPtCanvas, FinalEffMedIsoPt_,
         1, 1, 1, kBlack, 1, 20, "Pt(PFTau + Med Iso + DM / nJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffTightIsoPtCanvas, FinalEffTightIsoPt_,
         1, 1, 1, kBlack, 1, 20, "Pt(PFTau + Tight Iso + DM / nJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalEffDMFindPtCanvas, FinalEffDMFindPt_,
         1, 1, 1, kBlack, 1, 20, "Pt(PFTau + DecayModeFinding / nJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngDMEtaCanvas, OneProngDMEta_,
         1, 0, 0, kBlack, 7, 20, "1 Prong #eta(PFTau + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngOnePizDMEtaCanvas, OneProngOnePizDMEta_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 1 #pi^{0} #eta(PFTau + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngTwoPizDMEtaCanvas, OneProngTwoPizDMEta_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 2 #pi^{0} #eta(PFTau + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(ThreeProngDMEtaCanvas, ThreeProngDMEta_,
         1, 0, 0, kBlack, 7, 20, "3 Prong #eta(PFTau + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngDMPtCanvas, OneProngDMPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong p_{T}(PFTau + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngOnePizDMPtCanvas, OneProngOnePizDMPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 1 #pi^{0} p_{T}(PFTau + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngTwoPizDMPtCanvas, OneProngTwoPizDMPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 2 #pi^{0} p_{T}(PFTau + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(ThreeProngDMPtCanvas, ThreeProngDMPt_,
         1, 0, 0, kBlack, 7, 20, "3 Prong #eta(PFTau + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngDMEtaCanvas, FinalOneProngDMEta_,
         1, 1, 1, kBlack, 1, 20, "1 Prong #eta(PFTau + DM / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngOnePizDMEtaCanvas, FinalOneProngOnePizDMEta_,
         1, 1, 1, kBlack, 1, 20, "1 Prong + 1 #pi^{0} #eta(PFTau + DM / ak4PfJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngTwoPizDMEtaCanvas, FinalOneProngTwoPizDMEta_,
         1, 1, 1, kBlack, 1, 20, "1 Prong + 2 #pi^{0} #eta(PFTau + DM / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalThreeProngDMEtaCanvas, FinalThreeProngDMEta_,
         1, 1, 1, kBlack, 1, 20, "3 Prong #eta(PFTau + DM / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngDMPtCanvas, FinalOneProngDMPt_,
         1, 1, 1, kBlack, 1, 20, "1 Prong p_{T}(PFTau + DM / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngOnePizDMPtCanvas, FinalOneProngOnePizDMPt_,
         1, 1, 1, kBlack, 1, 20, "1 Prong + 1 #pi^{0} p_{T}(PFTau + DM / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngTwoPizDMPtCanvas, FinalOneProngTwoPizDMPt_,
         1, 1, 1, kBlack, 1, 20, "1 Prong + 2 #pi^{0} p_{T}(PFTau + DM / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalThreeProngDMPtCanvas, FinalThreeProngDMPt_,
         1, 1, 1, kBlack, 1, 20, "3 Prong p_{T}(PFTau + DM / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngMedIsoEtaCanvas, OneProngMedIsoEta_,
         1, 0, 0, kBlack, 7, 20, "1 Prong #eta(PFTau + DM + MedIso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngOnePizMedIsoEtaCanvas, OneProngOnePizMedIsoEta_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 1 #pi^{0} #eta(PFTau + DM + MedIso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngTwoPizMedIsoEtaCanvas, OneProngTwoPizMedIsoEta_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 2 #pi^{0} #eta(PFTau + DM + MedIso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(ThreeProngMedIsoEtaCanvas, ThreeProngMedIsoEta_,
         1, 0, 0, kBlack, 7, 20, "3 Prong #eta(PFTau + DM + MedIso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngMedIsoPtCanvas, OneProngMedIsoPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong p_{T}(PFTau + DM + MedIso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngOnePizMedIsoPtCanvas, OneProngOnePizMedIsoPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 1 #pi^{0} p_{T}(PFTau + DM + MedIso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(OneProngTwoPizMedIsoPtCanvas, OneProngTwoPizMedIsoPt_,
         1, 0, 0, kBlack, 7, 20, "1 Prong + 2 #pi^{0} p_{T}(PFTau + DM + MedIso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(ThreeProngMedIsoPtCanvas, ThreeProngMedIsoPt_,
         1, 0, 0, kBlack, 7, 20, "3 Prong #eta(PFTau + DM + MedIso)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngMedIsoEtaCanvas, FinalOneProngMedIsoEta_,
         1, 1, 1, kBlack, 1, 20, "1 Prong #eta(PFTau + DM + MedIso / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngOnePizMedIsoEtaCanvas, FinalOneProngOnePizMedIsoEta_,
         1, 1, 1, kBlack, 1, 20, "1 Prong + 1 #pi^{0} #eta(PFTau + DM + MedIso / ak4PfJets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngTwoPizMedIsoEtaCanvas, FinalOneProngTwoPizMedIsoEta_,
         1, 1, 1, kBlack, 1, 20, "1 Prong + 2 #pi^{0} #eta(PFTau + DM + MedIso / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalThreeProngMedIsoEtaCanvas, FinalThreeProngMedIsoEta_,
         1, 1, 1, kBlack, 1, 20, "3 Prong #eta(PFTau + DM + MedIso / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngMedIsoPtCanvas, FinalOneProngMedIsoPt_,
         1, 1, 1, kBlack, 1, 20, "1 Prong p_{T}(PFTau + DM + MedIso / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngOnePizMedIsoPtCanvas, FinalOneProngOnePizMedIsoPt_,
         1, 1, 1, kBlack, 1, 20, "1 Prong + 1 #pi^{0} p_{T}(PFTau + DM + MedIso / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalOneProngTwoPizMedIsoPtCanvas, FinalOneProngTwoPizMedIsoPt_,
         1, 1, 1, kBlack, 1, 20, "1 Prong + 2 #pi^{0} p_{T}(PFTau + DM + MedIso / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndTGraphAsym(FinalThreeProngMedIsoPtCanvas, FinalThreeProngMedIsoPt_,
         1, 1, 1, kBlack, 1, 20, "3 Prong p_{T}(PFTau + DM + MedIso / ak4Jets)", .04, .04, 1.1,  "#epsilon", .04, .04, 1.0, false);

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
 
  EtavsPtTauLooseIso_->Write();
  EtavsPtTauMedIso_->Write();
  EtavsPtTauTightIso_->Write();
  EtavsPtTauDMFind_->Write();
  EtavsPtJet_->Write();
  EtavsPtJetSoftMuon_->Write();
  EtavsPtJetSoftMuon_noMu_->Write();

  TauLooseIsoEta_->Write();
  TauMedIsoEta_->Write();
  TauTightIsoEta_->Write();
  TauDMFindEta_->Write();
  JetEta_->Write();
  JetEtaWithSoftMuon_->Write();
  JetEtaWithSoftMuon_noMu_->Write();
  FinalEffLooseIsoEta_->Write();
  FinalEffMedIsoEta_->Write();
  FinalEffTightIsoEta_->Write();
  FinalEffDMFindEta_->Write();

  TauLooseIsoPt_->Write();
  TauMedIsoPt_->Write();
  TauTightIsoPt_->Write();
  TauDMFindPt_->Write();
  JetPt_->Write();
  JetPtWithSoftMuon_->Write();
  JetPtWithSoftMuon_noMu_->Write();
  FinalEffLooseIsoPt_->Write();
  FinalEffMedIsoPt_->Write();
  FinalEffTightIsoPt_->Write();
  FinalEffDMFindPt_->Write();

  OneProngDMEta_->Write();
  OneProngOnePizDMEta_->Write();
  OneProngTwoPizDMEta_->Write();
  ThreeProngDMEta_->Write();
  OneProngDMPt_->Write();
  OneProngOnePizDMPt_->Write();
  OneProngTwoPizDMPt_->Write();
  ThreeProngDMPt_->Write();
  FinalOneProngDMEta_->Write();
  FinalOneProngOnePizDMEta_->Write();
  FinalOneProngTwoPizDMEta_->Write();
  FinalThreeProngDMEta_->Write();
  FinalOneProngDMPt_->Write();
  FinalOneProngOnePizDMPt_->Write();
  FinalOneProngTwoPizDMPt_->Write();
  FinalThreeProngDMPt_->Write();

  OneProngMedIsoEta_->Write();
  OneProngOnePizMedIsoEta_->Write();
  OneProngTwoPizMedIsoEta_->Write();
  ThreeProngMedIsoEta_->Write();
  OneProngMedIsoPt_->Write();
  OneProngOnePizMedIsoPt_->Write();
  OneProngTwoPizMedIsoPt_->Write();
  ThreeProngMedIsoPt_->Write();
  FinalOneProngMedIsoEta_->Write();
  FinalOneProngOnePizMedIsoEta_->Write();
  FinalOneProngTwoPizMedIsoEta_->Write();
  FinalThreeProngMedIsoEta_->Write();
  FinalOneProngMedIsoPt_->Write();
  FinalOneProngOnePizMedIsoPt_->Write();
  FinalOneProngTwoPizMedIsoPt_->Write();
  FinalThreeProngMedIsoPt_->Write();

  out_->Write();
  out_->Close();
std::cout << "DONE" << std::endl;
}//EndJob

// ------------ method called when starting to processes a run  ------------
void FakeRateMiniAODGetRates::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void FakeRateMiniAODGetRates::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void FakeRateMiniAODGetRates::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void FakeRateMiniAODGetRates::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

//Delete Memory
void FakeRateMiniAODGetRates::reset(const bool doDelete)
{
  if ((doDelete) && (NEvents_ != NULL)) delete NEvents_;
  NEvents_ = NULL;
  if ((doDelete) && (InvMassTauMuMu1_ != NULL)) delete InvMassTauMuMu1_;
  InvMassTauMuMu1_ = NULL;
  if ((doDelete) && (InvMassMu1Mu2_ != NULL)) delete InvMassMu1Mu2_;
  InvMassMu1Mu2_ = NULL;
  if ((doDelete) && (InvMassTauMuMu2_ != NULL)) delete InvMassTauMuMu2_;
  InvMassTauMuMu2_ = NULL;

  if ((doDelete) && (EtavsPtTauLooseIso_ != NULL)) delete EtavsPtTauLooseIso_;
  EtavsPtTauLooseIso_ = NULL;
  if ((doDelete) && (EtavsPtTauMedIso_ != NULL)) delete EtavsPtTauMedIso_;
  EtavsPtTauMedIso_ = NULL;
  if ((doDelete) && (EtavsPtTauTightIso_ != NULL)) delete EtavsPtTauTightIso_;
  EtavsPtTauTightIso_ = NULL;
  if ((doDelete) && (EtavsPtTauDMFind_ != NULL)) delete EtavsPtTauDMFind_;
  EtavsPtTauDMFind_ = NULL;
  if ((doDelete) && (EtavsPtJet_ != NULL)) delete EtavsPtJet_;
  EtavsPtJet_ = NULL;
  if ((doDelete) && (EtavsPtJetSoftMuon_ != NULL)) delete EtavsPtJetSoftMuon_;
  EtavsPtJetSoftMuon_ = NULL;
  if ((doDelete) && (EtavsPtJetSoftMuon_noMu_ != NULL)) delete EtavsPtJetSoftMuon_noMu_;
  EtavsPtJetSoftMuon_noMu_ = NULL;

  if ((doDelete) && (TauLooseIsoEta_ != NULL)) delete TauLooseIsoEta_;
  TauLooseIsoEta_ = NULL;
  if ((doDelete) && (TauMedIsoEta_ != NULL)) delete TauMedIsoEta_;
  TauMedIsoEta_ = NULL;
  if ((doDelete) && (TauTightIsoEta_ != NULL)) delete TauTightIsoEta_;
  TauTightIsoEta_ = NULL;
  if ((doDelete) && (TauDMFindEta_ != NULL)) delete TauDMFindEta_;
  TauDMFindEta_ = NULL;
  if ((doDelete) && (JetEta_ != NULL)) delete JetEta_;
  JetEta_ = NULL;
  if ((doDelete) && (JetEtaWithSoftMuon_ != NULL)) delete JetEtaWithSoftMuon_;
  JetEtaWithSoftMuon_ = NULL;
  if ((doDelete) && (JetEtaWithSoftMuon_noMu_ != NULL)) delete JetEtaWithSoftMuon_noMu_;
  JetEtaWithSoftMuon_noMu_ = NULL;
  if ((doDelete) && (FinalEffLooseIsoEta_ != NULL)) delete FinalEffLooseIsoEta_;
  FinalEffLooseIsoEta_ = NULL;
  if ((doDelete) && (FinalEffMedIsoEta_ != NULL)) delete FinalEffMedIsoEta_;
  FinalEffMedIsoEta_ = NULL;
  if ((doDelete) && (FinalEffTightIsoEta_ != NULL)) delete FinalEffTightIsoEta_;
  FinalEffTightIsoEta_ = NULL;
  if ((doDelete) && (FinalEffDMFindEta_ != NULL)) delete FinalEffDMFindEta_;
  FinalEffDMFindEta_ = NULL;

  if ((doDelete) && (TauLooseIsoPt_ != NULL)) delete TauLooseIsoPt_;
  TauLooseIsoPt_ = NULL;
  if ((doDelete) && (TauMedIsoPt_ != NULL)) delete TauMedIsoPt_;
  TauMedIsoPt_ = NULL;
  if ((doDelete) && (TauTightIsoPt_ != NULL)) delete TauTightIsoPt_;
  TauTightIsoPt_ = NULL;
  if ((doDelete) && (TauDMFindPt_ != NULL)) delete TauDMFindPt_;
  TauDMFindPt_ = NULL;
  if ((doDelete) && (JetPt_ != NULL)) delete JetPt_;
  JetPt_ = NULL;
  if ((doDelete) && (JetPtWithSoftMuon_ != NULL)) delete JetPtWithSoftMuon_;
  JetPtWithSoftMuon_ = NULL;
  if ((doDelete) && (JetPtWithSoftMuon_noMu_ != NULL)) delete JetPtWithSoftMuon_noMu_;
  JetPtWithSoftMuon_noMu_ = NULL;
  if ((doDelete) && (FinalEffLooseIsoPt_ != NULL)) delete FinalEffLooseIsoPt_;
  FinalEffLooseIsoPt_ = NULL;
  if ((doDelete) && (FinalEffMedIsoPt_ != NULL)) delete FinalEffMedIsoPt_;
  FinalEffMedIsoPt_ = NULL;
  if ((doDelete) && (FinalEffTightIsoPt_ != NULL)) delete FinalEffTightIsoPt_;
  FinalEffTightIsoPt_ = NULL;
  if ((doDelete) && (FinalEffDMFindPt_ != NULL)) delete FinalEffDMFindPt_;
  FinalEffDMFindPt_ = NULL;

  if ((doDelete) && (OneProngDMEta_ != NULL)) delete OneProngDMEta_;
  OneProngDMEta_ = NULL;
  if ((doDelete) && (OneProngOnePizDMEta_ != NULL)) delete OneProngOnePizDMEta_;
  OneProngOnePizDMEta_ = NULL;
  if ((doDelete) && (OneProngTwoPizDMEta_ != NULL)) delete OneProngTwoPizDMEta_;
  OneProngOnePizDMEta_ = NULL;
  if ((doDelete) && (ThreeProngDMEta_ != NULL)) delete ThreeProngDMEta_;
  ThreeProngDMEta_ = NULL;
  if ((doDelete) && (OneProngDMPt_ != NULL)) delete OneProngDMPt_;
  OneProngDMPt_ = NULL;
  if ((doDelete) && (OneProngOnePizDMPt_ != NULL)) delete OneProngOnePizDMPt_;
  OneProngOnePizDMPt_ = NULL;
  if ((doDelete) && (OneProngTwoPizDMPt_ != NULL)) delete OneProngTwoPizDMPt_;
  OneProngTwoPizDMPt_ = NULL;
  if ((doDelete) && (ThreeProngDMPt_ != NULL)) delete ThreeProngDMPt_;
  ThreeProngDMPt_ = NULL;
  if ((doDelete) && (FinalOneProngDMEta_ != NULL)) delete FinalOneProngDMEta_;
  FinalOneProngDMEta_ = NULL;
  if ((doDelete) && (FinalOneProngOnePizDMEta_ != NULL)) delete FinalOneProngOnePizDMEta_;
  FinalOneProngOnePizDMEta_ = NULL;
  if ((doDelete) && (FinalOneProngTwoPizDMEta_ != NULL)) delete FinalOneProngTwoPizDMEta_;
  FinalOneProngTwoPizDMEta_ = NULL;
  if ((doDelete) && (FinalThreeProngDMEta_ != NULL)) delete FinalThreeProngDMEta_;
  FinalThreeProngDMEta_ = NULL;
  if ((doDelete) && (FinalOneProngDMPt_ != NULL)) delete FinalOneProngDMPt_;
  FinalOneProngDMPt_ = NULL;
  if ((doDelete) && (FinalOneProngOnePizDMPt_ != NULL)) delete FinalOneProngOnePizDMPt_;
  FinalOneProngOnePizDMPt_ = NULL;
  if ((doDelete) && (FinalOneProngTwoPizDMPt_ != NULL)) delete FinalOneProngTwoPizDMPt_;
  FinalOneProngTwoPizDMPt_ = NULL;
  if ((doDelete) && (FinalThreeProngDMPt_ != NULL)) delete FinalThreeProngDMPt_;
  FinalThreeProngDMPt_ = NULL;

  if ((doDelete) && (OneProngMedIsoEta_ != NULL)) delete OneProngMedIsoEta_;
  OneProngMedIsoEta_ = NULL;
  if ((doDelete) && (OneProngOnePizMedIsoEta_ != NULL)) delete OneProngOnePizMedIsoEta_;
  OneProngOnePizMedIsoEta_ = NULL;
  if ((doDelete) && (OneProngTwoPizMedIsoEta_ != NULL)) delete OneProngTwoPizMedIsoEta_;
  OneProngOnePizMedIsoEta_ = NULL;
  if ((doDelete) && (ThreeProngMedIsoEta_ != NULL)) delete ThreeProngMedIsoEta_;
  ThreeProngMedIsoEta_ = NULL;
  if ((doDelete) && (OneProngMedIsoPt_ != NULL)) delete OneProngMedIsoPt_;
  OneProngMedIsoPt_ = NULL;
  if ((doDelete) && (OneProngOnePizMedIsoPt_ != NULL)) delete OneProngOnePizMedIsoPt_;
  OneProngOnePizMedIsoPt_ = NULL;
  if ((doDelete) && (OneProngTwoPizMedIsoPt_ != NULL)) delete OneProngTwoPizMedIsoPt_;
  OneProngTwoPizMedIsoPt_ = NULL;
  if ((doDelete) && (ThreeProngMedIsoPt_ != NULL)) delete ThreeProngMedIsoPt_;
  ThreeProngMedIsoPt_ = NULL;
  if ((doDelete) && (FinalOneProngMedIsoEta_ != NULL)) delete FinalOneProngMedIsoEta_;
  FinalOneProngMedIsoEta_ = NULL;
  if ((doDelete) && (FinalOneProngOnePizMedIsoEta_ != NULL)) delete FinalOneProngOnePizMedIsoEta_;
  FinalOneProngOnePizMedIsoEta_ = NULL;
  if ((doDelete) && (FinalOneProngTwoPizMedIsoEta_ != NULL)) delete FinalOneProngTwoPizMedIsoEta_;
  FinalOneProngTwoPizMedIsoEta_ = NULL;
  if ((doDelete) && (FinalThreeProngMedIsoEta_ != NULL)) delete FinalThreeProngMedIsoEta_;
  FinalThreeProngMedIsoEta_ = NULL;
  if ((doDelete) && (FinalOneProngMedIsoPt_ != NULL)) delete FinalOneProngMedIsoPt_;
  FinalOneProngMedIsoPt_ = NULL;
  if ((doDelete) && (FinalOneProngOnePizMedIsoPt_ != NULL)) delete FinalOneProngOnePizMedIsoPt_;
  FinalOneProngOnePizMedIsoPt_ = NULL;
  if ((doDelete) && (FinalOneProngTwoPizMedIsoPt_ != NULL)) delete FinalOneProngTwoPizMedIsoPt_;
  FinalOneProngTwoPizMedIsoPt_ = NULL;
  if ((doDelete) && (FinalThreeProngMedIsoPt_ != NULL)) delete FinalThreeProngMedIsoPt_;
  FinalThreeProngMedIsoPt_ = NULL;

}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void FakeRateMiniAODGetRates::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FakeRateMiniAODGetRates);

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

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

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
#include "AnalyzerGeneratorRecoVariousFunctions/VariousFunctions/interface/BTagCalibrationStandalone.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

#include "RooArgSet.h"
#include "RooArgusBG.h"
#include "RooRealVar.h"
#include "RooDataSet.h"

#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooBreitWigner.h"

#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooMinuit.h"
#include "RooPlot.h"

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
      TH1F* Pileup_;
      
      //name of output root file
      std::string outFileName_;
      edm::EDGetTokenT<edm::View<pat::Muon> > mu12Tag_;
      edm::EDGetTokenT<edm::View<pat::Tau> > tauTag_;
      bool checkBTag_;
      string csvBTag_;
      string bTagSFShift_;
      edm::EDGetTokenT<edm::View<pat::Muon> > mu3Tag_;
      edm::EDGetTokenT<edm::View<pat::MET> > metTag_;
      edm::EDGetTokenT<edm::View<pat::Jet> > jetTag_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryInfo_;
      edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
      
      double tauPtCut_;
      double diMudRCut_;
      double mu3dROverlapCut_;
      double tauHadOverlapdRCut_;
      double xsec_; 
      double lumi_; 
      double summedWeights_;
      bool MC_;
      bool ApplyFR_;
      bool rooDataset_;
      std::string PileupFileName_;
      std::string _fpIDs_BToF;
      std::string _fpIDs_GH;
      std::string _fpISOs_BToF;
      std::string _fpISOs_GH;
      std::string _fpTrack;
      std::string _fpTrigger_BToF;
      std::string _fpTrigger_GH;
      std::string _fp_LowPt;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo>> PUTag_;
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

      RooRealVar *x = new RooRealVar("x","x", 0, 30);
      RooRealVar *w = new RooRealVar("w","w",-9999999, 99999999);
      RooDataSet *mumumass_dataset = new RooDataSet("mumumass_dataset", "mumumass_dataset", RooArgSet(*x,*w));

      RooRealVar *y = new RooRealVar("y","y", 0, 30);
      RooDataSet *mumutautaumass_dataset = new RooDataSet("mumutautaumass_dataset", "mumutautaumass_dataset", RooArgSet(*x,*y,*w));

      RooRealVar *y1 = new RooRealVar("y1","y1", 0, 1000);
      RooDataSet *mumufourBodymass_dataset = new RooDataSet("mumufourBodymass_dataset", "mumufourBodymass_dataset", RooArgSet(*x,*y1,*w));

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
      TH1F* Mu3Pt_;
      TH1F* TauHPt_;
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
      TH1F* SFWeights_;
      TH1F* GenWeights_;
      TH1F* HighestCSVInclJet_;
      TH1F* HighestCombMVAJet_;
      TH1F* ZMassdR_;
      TH2F* DiTauMassVSMVA_;
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
  checkBTag_(iConfig.getParameter<bool>("checkBTag")),
  csvBTag_(iConfig.getParameter<std::string>("csvBTag")),
  bTagSFShift_(iConfig.getParameter<std::string>("bTagSFShift")), 
  mu3Tag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("mu3Tag"))),
  metTag_(consumes<edm::View<pat::MET> >(iConfig.getParameter<edm::InputTag>("metTag"))),
  jetTag_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetTag"))),
  pileupSummaryInfo_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummaryInfo"))),
  genEventInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfoToken"))),
  tauPtCut_(iConfig.getParameter<double>("tauPtCut")),
  diMudRCut_(iConfig.getParameter<double>("diMudRCut")),
  mu3dROverlapCut_(iConfig.getParameter<double>("mu3dROverlapCut")),
  tauHadOverlapdRCut_(iConfig.getParameter<double>("tauHadOverlapdRCut")),
  xsec_(iConfig.getParameter<double>("xsec")),
  lumi_(iConfig.getParameter<double>("lumi")),
  summedWeights_(iConfig.getParameter<double>("summedWeights")),
  MC_(iConfig.getParameter<bool>("MC")),
  ApplyFR_(iConfig.getParameter<bool>("ApplyFR")),
  rooDataset_(iConfig.getParameter<bool>("rooDataset")),
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
std::cout << "check1" << std::endl;
  PileupFile = new TFile(PileupFileName_.c_str());
  Pileup_ = (TH1F*)PileupFile->Get("PileupWeights");

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
  NEvents_->Fill(0);
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

  edm::Handle<edm::View<pat::MET> > pMET;
  iEvent.getByToken(metTag_, pMET);
  pat::MET MET = pMET->at(0);

///////////////////////////////////////
// Get PU and Gen Weights if MC
///////////////////////////////////////
  double pileupWeight = 1.0, genWeight = 1.0;
  float nTrueVertices = -1;
  if (MC_)
  {
    edm::Handle<std::vector<PileupSummaryInfo> > pPileupSummaryInfo;
    iEvent.getByToken(pileupSummaryInfo_, pPileupSummaryInfo);
    for(vector<PileupSummaryInfo>::const_iterator cand=pPileupSummaryInfo->begin(); cand!=pPileupSummaryInfo->end();++cand)
    {
      if (cand->getBunchCrossing() == 0)
      { 
        nTrueVertices=cand->getTrueNumInteractions();
        break;
      }
    }//for cand
    if (nTrueVertices != -1)
    {
      int binx =  Pileup_->GetXaxis()->FindBin(nTrueVertices);
      pileupWeight = Pileup_->GetBinContent(binx);
    }//if nTrueVertices != -1
  
    edm::Handle<GenEventInfoProduct> genEventInfo;
    iEvent.getByToken(genEventInfoToken_, genEventInfo);
    double eventGenWeight = genEventInfo->weight();
    genWeight = eventGenWeight * lumi_ * xsec_ / summedWeights_;
  }//if MC_

///////////////////////////////////////
// Get Mu1 and Mu2 adn do InvM check
///////////////////////////////////////
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

  if (mu1.pt() < 26.0) 
    return;

  std::cout << "mu1.pt()=" << mu1.pt() << "\tmu1.eta()=" << mu1.eta() << "\tmu1.phi()=" << mu1.phi() << std::endl;
  std::cout << "mu2.pt()=" << mu2.pt() << "\tmu2.eta()=" << mu2.eta() << "\tmu2.phi()=" << mu2.phi() << std::endl;
  reco::LeafCandidate::LorentzVector diMuP4 = mu1.p4() + mu2.p4();
  std::cout << "diMuP4.M()=" << diMuP4.M() << std::endl;
  std::cout << "diMudR= " << deltaR(mu1, mu2) << std::endl;
  if (diMuP4.M() > 30.0 || deltaR(mu1, mu2) > diMudRCut_)
  {
    NEvents_->Fill(1);
    return;
  }//

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
//  std::cout << "LumiFraction_BF=" << LumiFraction_BF << "\tLumiFraction_GH=" << LumiFraction_GH << "\nIDsWeight_BF->GetBinContent(binxIDs_BF, binyIDs_BF)=" << IDsWeight_BF->GetBinContent(binxIDs_BF, binyIDs_BF) << "\tIDsWeight_GH->GetBinContent(binxIDs_GH, binyIDs_GH)=" << IDsWeight_GH->GetBinContent(binxIDs_GH, binyIDs_GH) << "\nISOsWeight_BF->GetBinContent(binxISOs_BF, binyISOs_BF)=" << ISOsWeight_BF->GetBinContent(binxISOs_BF, binyISOs_BF) << "\tISOsWeight_GH->GetBinContent(binxISOs_GH, binyISOs_GH)=" << ISOsWeight_GH->GetBinContent(binxISOs_GH, binyISOs_GH) << "\nbinxIDs_BF=" << binxIDs_BF << "\tbinyIDs_BF=" << binyIDs_BF << "\tbinxIDs_GH=" << binxIDs_GH << "\tbinyIDs_GH=" << binyIDs_GH << "\nbinxISOs_BF=" << binxISOs_BF << "\tbinyISOs_BF=" << binyISOs_BF << "\tbinxISOs_GH=" << binxISOs_GH << "\tbinyISOs_GH=" << binyISOs_GH << std::endl;
  
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
  IDs_weightGH = IDs_weightGH * IDsWeight_GH->GetBinContent(binxIDs_GH, binyIDs_GH);
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

    binyIDs_BF = IDsWeight_LowPt_ID->GetYaxis()->FindBin(mu2.pt());
    binxIDs_BF = IDsWeight_LowPt_ID->GetXaxis()->FindBin(fabs(mu2.eta()));
    binyIDs_GH = IDsWeight_LowPt_ID->GetYaxis()->FindBin(mu2.pt());
    binxIDs_GH = IDsWeight_LowPt_ID->GetXaxis()->FindBin(fabs(mu2.eta()));
    binyISOs_BF = ISOsWeight_LowPt_ISO->GetYaxis()->FindBin(mu2.pt());
    binxISOs_BF = ISOsWeight_LowPt_ISO->GetXaxis()->FindBin(fabs(mu2.eta()));
    binyISOs_GH = ISOsWeight_LowPt_ISO->GetYaxis()->FindBin(mu2.pt());
    binxISOs_GH = ISOsWeight_LowPt_ISO->GetXaxis()->FindBin(fabs(mu2.eta()));
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

 std::cout << "check6" << std::endl;

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

//////////////////////////////
// BTaggin Scale factors
//////////////////////////////
  // BTaggin SF
  std::cout << "===> Loading the input .csv SF file..." << std::endl;
  
  std::string inputCSVfile = "/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/CSVv2_Moriond17_B_H.csv";  
  BTagCalibration calib("csvv2", inputCSVfile);
  std::cout << "bTagSFShift_= " << bTagSFShift_ << "  BTagEntry::FLAV_C=" << BTagEntry::FLAV_C << "  BTagEntry::FLAV_B=" << BTagEntry::FLAV_B <<  std::endl;
  BTagCalibrationReader reader(BTagEntry::OP_MEDIUM, bTagSFShift_) ;
  
  reader.load(calib, BTagEntry::FLAV_B, "comb");
  reader.load(calib, BTagEntry::FLAV_C, "comb");
  reader.load(calib, BTagEntry::FLAV_UDSG, "incl");
  std::cout  << "initialized reader" << std::endl;

///////////////////////////////////////
// Combining weights and SF
///////////////////////////////////////
  std::cout << "\n\nIDs_weightBF=" << IDs_weightBF << "\tISOs_weightGH=" << ISOs_weightGH << "\tTracks_weight=" << Tracks_weight << "\tTrigger_weightGH=" << Trigger_weightGH << std::endl;
  double SFsWeight = LumiFraction_BF * (IDs_weightBF * ISOs_weightBF * Trigger_weightBF) + 
                     LumiFraction_GH * (IDs_weightGH * ISOs_weightGH * Trigger_weightGH);  
  if (!MC_)
    SFsWeight = 1.0;
  std::cout << "pileupWeight=" << pileupWeight << "\tgenWeight=" << genWeight << "\tSFsWeight=" << SFsWeight << std::endl;

///////////////////////////////////////
// Do DiTau Side of selection
///////////////////////////////////////
  reco::LeafCandidate::LorentzVector diTauP4, fourBody, mu1mu3, mu2mu3;
  double bestdR = 100000000000000000, bTagValue = -1, csvSF = -1;
  bool checkEventDiTau = false;
  unsigned int TauRemovedMuCount = 0;
  pat::Muon mu3; 
  pat::Tau tau;
  std::cout << "pTaus->size()= " << pTaus->size() << "\tpMu3->size()= " << pMu3->size() << std::endl;
  for (edm::View<pat::Tau>::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); ++iTau)
  {
    if (iTau->pt() < tauPtCut_ || fabs(iTau->eta() ) > 2.4 || deltaR(*iTau,mu1) < tauHadOverlapdRCut_ || 
        deltaR(*iTau, mu2) < tauHadOverlapdRCut_ || iTau->tauID("byIsolationMVArun2v1DBoldDMwLTraw") <= -0.5)
      continue;
    std::cout << "iTau=>pt=" << iTau->pt() << std::endl;
    bool checkMu3Removed = false;
    for (edm::View<pat::Muon>::const_iterator iMu = pMu3->begin(); iMu != pMu3->end(); ++iMu)
    {
      double currdR = deltaR(*iTau, *iMu);
      std::cout  << "\tdR=" << currdR << std::endl;
      if (currdR < .8 && currdR < bestdR && deltaR(*iMu, mu1) > mu3dROverlapCut_ && deltaR(*iMu, mu2) > mu3dROverlapCut_)   
      {
        bestdR = currdR;
        diTauP4 = iTau->p4() + iMu->p4();
        mu3 = *iMu; 
        tau = *iTau;
	mu1mu3 = iMu->p4() + mu1.p4();   
	mu2mu3 = iMu->p4() + mu2.p4();
        fourBody = iTau->p4() + iMu->p4() + mu1.p4() + mu2.p4();
        checkMu3Removed = true;
        checkEventDiTau = true;
        double smallestdRJet = 10000000000000000;
        for (edm::View<pat::Jet>::const_iterator iJet = pJets->begin(); iJet != pJets->end(); ++iJet)
        {
          double dRJet = deltaR(*iTau, *iJet);
          if (dRJet < smallestdRJet)
          {
            bTagValue = iJet->bDiscriminator(csvBTag_);
            smallestdRJet = dRJet;
            int hadronFlavor = iJet->hadronFlavour();
            if ( hadronFlavor == 5 )
              csvSF = reader.eval(BTagEntry::FLAV_B, fabs(iJet->eta()), iJet->pt(), bTagValue);

            else if( hadronFlavor == 4 )
              csvSF = reader.eval(BTagEntry::FLAV_C, fabs(iJet->eta()), iJet->pt(), bTagValue);

            else 
              csvSF = reader.eval(BTagEntry::FLAV_UDSG, fabs(iJet->eta()), iJet->pt(), bTagValue);
	    std::cout << "csvSF=" << csvSF << "  hadronFlavor=" << hadronFlavor << " bTagValue=" << bTagValue << std::endl;
          }//if dRJet < smallestdRJet
        }//for iJet
      }//if
    }//for iMu
    if (checkMu3Removed)
      TauRemovedMuCount++;
  }//for iTau

  double csvWeight = 1;
  if (bTagValue >= 0.8484)
    csvWeight = 1 - csvSF;

   std::cout << "csvWeight= " << csvWeight << std::endl;

  std::cout << "FoundDiTau=" << checkEventDiTau << std::endl;
  double tauMedSF = 0.97;
  if (checkEventDiTau)
  {
    NEvents_->Fill(3);
    DiTauDiMuMass_->Fill(fourBody.M(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    if (bestdR < .4)
    {
      DiMuPtSmallerdR_->Fill(diMuP4.Pt(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
      DiTauMassSmallerdR_->Fill(diTauP4.M(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    }//if
    DiMuDiTauDeltaPhi_->Fill(fabs(diTauP4.Phi() - diMuP4.Phi() ), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    InvMassTauHadMu3_->Fill(diTauP4.M(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    InvMassMu1TauMu_->Fill(mu1mu3.M(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    if (mu1mu3.M() > 60.0 && mu1mu3.M() < 120 )
    {
      double dRZ = deltaR(mu1, mu3);
      ZMassdR_->Fill(dRZ, pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    }//if 
    if (mu2mu3.M() > 60.0 && mu2mu3.M() < 120 )
    {
      double dRZ = deltaR(mu2, mu3);
      ZMassdR_->Fill(dRZ, pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    }//if 
    InvMassMu2TauMu_->Fill(mu2mu3.M() );
    PtTauHadMu3_->Fill(diTauP4.Pt(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    DiTauEta_->Fill(diTauP4.Eta(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    DiTauPhi_->Fill(diTauP4.Phi(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    METDiTauDeltaPhi_->Fill(fabs(diTauP4.Phi() - MET.phi() ), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    InvMassZPeakRange_->Fill(diMuP4.M(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    InvMassUpsilonRange_->Fill(diMuP4.M(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    InvMassFullRange_->Fill(diMuP4.M(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    GenWeights_->Fill(genWeight);
    PileupWeights_->Fill(pileupWeight);
    SFWeights_->Fill(SFsWeight*tauMedSF);
    if (mu1.eta() < .9 && mu2.eta() < .9)
      InvMassDiMuBarrBarr_->Fill(diMuP4.M(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    else if ( (mu1.eta() < .9 && mu2.eta() < 1.2 && mu2.eta() > .9) || (mu1.eta() < 1.2 && mu1.eta() > .9 && mu2.eta() < .9) )
      InvMassDiMuBarrOver_->Fill(diMuP4.M(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    else if ( (mu1.eta() < .9 && mu2.eta() > 1.2) || (mu1.eta() > 1.2 && mu2.eta() < .9) )
      InvMassDiMuBarrEndc_->Fill(diMuP4.M(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    else if ( (mu1.eta() > 1.2 && mu2.eta() < 1.2 && mu2.eta() > .9) || (mu1.eta() < 1.2 && mu1.eta() > .9 && mu2.eta() > 1.2) )
      InvMassDiMuEndcOver_->Fill(diMuP4.M(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    else if ( mu2.eta() < 1.2 && mu2.eta() > .9 && mu1.eta() < 1.2 && mu1.eta() > .9)
      InvMassDiMuOverOver_->Fill(diMuP4.M(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    else if (mu1.eta() > 1.2 && mu2.eta() > 1.2)
      InvMassDiMuEndcEndc_->Fill(diMuP4.M(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    Mu1Pt_->Fill(mu1.pt(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    Mu2Pt_->Fill(mu2.pt(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    Mu3Pt_->Fill(mu3.pt(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    TauHPt_->Fill(tau.pt(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    PtMu1vsPtMu2_->Fill(mu1.pt(), mu2.pt() , pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    DiMuPt_->Fill(diMuP4.Pt(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    Mu1Eta_->Fill(mu1.eta(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    Mu2Eta_->Fill(mu2.eta(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    DiMuEta_->Fill(diMuP4.Eta(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    DiMuPhi_->Fill(diMuP4.Phi(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    DiMudR_->Fill(reco::deltaR(mu1, mu2), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    DiTaudR_->Fill(bestdR, pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF);
    EtMET_->Fill(MET.pt(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    METDiMuDeltaPhi_->Fill(fabs(diMuP4.Phi() - MET.phi() ), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    EtMET_->Fill(MET.pt(), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    METDiMuDeltaPhi_->Fill(fabs(diMuP4.Phi() - MET.phi() ), pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF );
    DiTauMassVSMVA_->Fill(diMuP4.M(), tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw"));
  }//for checkMu3Removed 
  else
    NEvents_->Fill(2);
  NMedIsoTausWithMu3_->Fill(TauRemovedMuCount ); 

  if (rooDataset_ && checkEventDiTau)
  {
    x->setVal(diMuP4.M() );
    y->setVal(diTauP4.M() );
    y1->setVal(fourBody.M() );
    w->setVal(pileupWeight*csvWeight*genWeight*SFsWeight*tauMedSF*.001 );
    mumumass_dataset->add(RooArgSet(*x,*w));
    mumutautaumass_dataset->add(RooArgSet(*x,*y,*w));
    mumufourBodymass_dataset->add(RooArgSet(*x,*y1,*w));
  }//if
}//End InvMass::analyze


// ------------ method called once each job just before starting event loop  ------------
void DiMu_ExtraPlots::beginJob()
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
  Mu3Pt_     = new TH1F("Mu3Pt"    , "", 50, 0, 500);
  TauHPt_     = new TH1F("TauHPt"    , "", 50, 0, 500);
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
  SFWeights_     = new TH1F("SFWeights"    , "", 10000, 0, 1000);
  GenWeights_     = new TH1F("GenWeights"    , "", 20000, -10000, 10000);
  HighestCSVInclJet_     = new TH1F("HighestCSVInclJet"    , "", 100, 0, 1);
  HighestCombMVAJet_     = new TH1F("HighestCombMVAJet"    , "", 100, 0, 1);
  ZMassdR_     = new TH1F("ZMassdR"    , "", 100, 0, 6);
  DiTauMassVSMVA_     = new TH2F("DiTauMassVSMVA"    , "", 50, 0, 25, 50, -10, 1);
}

// ------------ method called once each job just after ending the event loop  ------------
void DiMu_ExtraPlots::endJob()
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
  TCanvas Mu3PtCanvas_("Mu3PtCanvas","",600,600);
  TCanvas TauHPtCanvas_("TauHPtCanvas","",600,600);
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
  TCanvas SFWeightsCanvas_("SFWeightsCanvas","",600,600);
  TCanvas GenWeightsCanvas_("GenWeightsCanvas","",600,600);
  TCanvas HighestCSVInclJetCanvas_("HighestCSVInclJetCanvas","",600,600);
  TCanvas HighestCombMVAJetCanvas_("HighestCombMVAJetCanvas","",600,600);
  TCanvas ZMassdRCanvas_("ZMassdRCanvas","",600,600);
  TCanvas DiTauMassVSMVACanvas_("DiTauMassVSMVACanvas","",600,600);

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
  VariousFunctions::formatAndDrawCanvasAndHist1D(Mu3PtCanvas_, Mu3Pt_,
	 1, 0, 0, kBlack, .1, 20, "p_{T}(#mu_{3})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauHPtCanvas_, TauHPt_,
	 1, 0, 0, kBlack, .1, 20, "p_{T}(#tau_{h})", .04, .04, 1.1,  "", .04, .04, 1.0, false);
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
  VariousFunctions::formatAndDrawCanvasAndHist1D(SFWeightsCanvas_, SFWeights_,
         1, 0, 0, kBlack, .1, 20, "Scale Factor Weight (ID,Iso,Track,Trigger)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(GenWeightsCanvas_, GenWeights_,
         1, 0, 0, kBlack, .1, 20, "Gen Weight (EventWeight * LumiData / LumiMC)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(HighestCSVInclJetCanvas_, HighestCSVInclJet_,
         1, 0, 0, kBlack, .1, 20, "highest CSV Inclusive value jet in Event", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(HighestCombMVAJetCanvas_, HighestCombMVAJet_,
         1, 0, 0, kBlack, .1, 20, "highest Combined MVA value jet in Event", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(ZMassdRCanvas_, ZMassdR_,
         1, 0, 0, kBlack, .1, 20, "#DeltaR #mu_{1} #tau_{#mu}", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist2D(DiTauMassVSMVACanvas_, DiTauMassVSMVA_,
         1, 0, 0, kBlack, 7, 20, "M(tautau)", .04, .04, 1.1, "MVA", .04, .04, 1.6, "", .04, .04, 1.0);


std::cout << "after formatting" << std::endl;
  
////////////////////////
// For Pt Gen
////////////////////////
std::cout << "<----------------Formatted Canvases and Histos-------------->" << std::endl;

  //Write output file
  out_->cd();
  if (rooDataset_)
  {
    mumumass_dataset->Write();
    mumutautaumass_dataset->Write();
    mumufourBodymass_dataset->Write();
  }//if
//  NMedIsoTausWithMu3Canvas_.Write();
//  InvMassMu1TauMuCanvas_.Write();
//  InvMassMu2TauMuCanvas_.Write();
//  InvMassTauHadMu3Canvas_.Write();
//  PtTauHadMu3Canvas_.Write();
//  InvMassUpsilonRangeCanvas_.Write();
//  InvMassZPeakRangeCanvas_.Write();
//  InvMassFullRangeCanvas_.Write();
//  InvMassDiMuBarrBarrCanvas_.Write();
//  InvMassDiMuBarrEndcCanvas_.Write();
//  InvMassDiMuEndcEndcCanvas_.Write();
//  InvMassDiMuBarrOverCanvas_.Write();
//  InvMassDiMuOverOverCanvas_.Write();
//  InvMassDiMuEndcOverCanvas_.Write();
//  Mu1PtCanvas_.Write();
//  Mu2PtCanvas_.Write();
//  DiMuPtCanvas_.Write();
//  DiMuPtSmallerdRCanvas_.Write();
//  Mu1EtaCanvas_.Write();
//  Mu2EtaCanvas_.Write();
//  DiTauEtaCanvas_.Write();
//  DiTaudRCanvas_.Write();
//  DiTauPhiCanvas_.Write();
//  DiMuEtaCanvas_.Write();
//  DiMudRCanvas_.Write();
//  DiMuPhiCanvas_.Write();
//  EtMETCanvas_.Write();
//  DiTauDiMuMassCanvas_.Write();
//  DiTauMassSmallerdRCanvas_.Write();
//  DiMuDiTauDeltaPhiCanvas_.Write();
//  METDiTauDeltaPhiCanvas_.Write();
//  METDiMuDeltaPhiCanvas_.Write();
//  PtMu1vsPtMu2Canvas_.Write();
//  PileupWeightsCanvas_.Write();
//  SFWeightsCanvas_.Write();
//  GenWeightsCanvas_.Write();
//  HighestCSVInclJetCanvas_.Write();
//  HighestCombMVAJetCanvas_.Write();
//  ZMassdRCanvas_.Write();

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
  Mu3Pt_->Write();
  TauHPt_->Write();
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
  SFWeights_->Write();
  GenWeights_->Write();
  HighestCSVInclJet_->Write();
  HighestCombMVAJet_->Write();
  ZMassdR_->Write();
  DiTauMassVSMVA_->Write();


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
  if ((doDelete) && (Mu3Pt_ != NULL)) delete Mu3Pt_;
  Mu3Pt_ = NULL;
  if ((doDelete) && (TauHPt_ != NULL)) delete TauHPt_;
  TauHPt_ = NULL;
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
  if ((doDelete) && (SFWeights_ != NULL)) delete SFWeights_;
  SFWeights_ = NULL;
  if ((doDelete) && (GenWeights_ != NULL)) delete GenWeights_;
  GenWeights_ = NULL;
  if ((doDelete) && (HighestCSVInclJet_ != NULL)) delete HighestCSVInclJet_;
  HighestCSVInclJet_ = NULL;
  if ((doDelete) && (HighestCombMVAJet_ != NULL)) delete HighestCombMVAJet_;
  HighestCombMVAJet_ = NULL;
  if ((doDelete) && (ZMassdR_ != NULL)) delete ZMassdR_;
  ZMassdR_ = NULL;
  if ((doDelete) && (DiTauMassVSMVA_ != NULL)) delete DiTauMassVSMVA_;
  DiTauMassVSMVA_ = NULL;

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

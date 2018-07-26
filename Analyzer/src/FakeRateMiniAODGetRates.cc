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
      std::string medIsoTag_;
      std::string decayModeFindingTag_;
      double mu3dRMin_;
      double mu3dRMax_;
      double tauPtCut_;
      double mu3dROverlapCut_;
      std::string minMVARaw_;
      std::string medIsoTau_;
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
      TGraph *TrackWeight_eta;
      TGraph *TrackWeight_vtx;
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
      std::list<TrackProperties> TrackCorr_eta, TrackCorr_vtx;
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
      TH2F* MVARawvsMVADisc_;

      TH2F* EtavsPtTauMedIso_;
      TH2F* EtavsPtTauDMFind_;

      TH1F* TauMedIsoEta_;
      TH1F* TauDMFindEta_;

      TH1F* TauMedIsoPt_;
      TH1F* TauDMFindPt_;
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
  medIsoTag_(iConfig.getParameter<std::string>("medIsoTag")),
  decayModeFindingTag_(iConfig.getParameter<std::string>("decayModeFindingTag")),
  mu3dRMin_(iConfig.getParameter<double>("mu3dRMin")),
  mu3dRMax_(iConfig.getParameter<double>("mu3dRMax")),
  tauPtCut_(iConfig.getParameter<double>("tauPtCut")),
  mu3dROverlapCut_(iConfig.getParameter<double>("mu3dROverlapCut")),
  minMVARaw_(iConfig.getParameter<std::string>("minMVARaw")),
  medIsoTau_(iConfig.getParameter<std::string>("medIsoTau")),
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


  IDsWeight_BToF =  (TH2F*)_fileIDs_BToF->Get("EtavsPtLooseID");
  IDsWeight_GH =  (TH2F*)_fileIDs_GH->Get("EtavsPtLooseID");
  ISOsWeight_BToF = (TH2F*)_fileISOs_BToF->Get("EtavsPtLooseISO");
  ISOsWeight_GH = (TH2F*)_fileISOs_GH->Get("EtavsPtLooseISO");
  TrackWeight_eta = (TGraph*)_fileTrack->Get("ratio_eff_eta3_dr030e030_corr");
  TrackWeight_vtx = (TGraph*)_fileTrack->Get("ratio_eff_vtx_dr030e030_corr");
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
  std::cout << "TrackWeight_eta =" << TrackWeight_eta << std::endl;
  std::cout << "TrackWeight_vtx =" << TrackWeight_vtx << std::endl;
  std::cout << "TriggerWeight_BToF =" << TriggerWeight_BToF << std::endl;

  Double_t x[TrackWeight_eta->GetN()], y[TrackWeight_eta->GetN()];
  for(int i=0; i<TrackWeight_eta->GetN(); i++){
     TrackWeight_eta->GetPoint(i, x[i], y[i]);
     cout<<"x="<<x[i]<<"; "<<"y="<<y[i]<<"."<<std::endl;
     TrackProperties val;
     val.x=x[i];
     val.y=y[i];
     val.errx_up=TrackWeight_eta->GetErrorXhigh(i);
     val.errx_down=TrackWeight_eta->GetErrorXlow(i);
     val.erry_up=TrackWeight_eta->GetErrorYhigh(i);
     val.erry_down=TrackWeight_eta->GetErrorYlow(i);
     TrackCorr_eta.push_back(val);
  }

  Double_t xval[TrackWeight_vtx->GetN()], yval[TrackWeight_vtx->GetN()];
  for(int i=0; i<TrackWeight_vtx->GetN(); i++){
     TrackWeight_vtx->GetPoint(i, xval[i], yval[i]);
     cout<<"xval="<<xval[i]<<"; "<<"yval="<<yval[i]<<"."<<std::endl;
     TrackProperties val;
     val.x=xval[i];
     val.y=yval[i];
     val.errx_up=TrackWeight_vtx->GetErrorXhigh(i);
     val.errx_down=TrackWeight_vtx->GetErrorXlow(i);
     val.erry_up=TrackWeight_vtx->GetErrorYhigh(i);
     val.erry_down=TrackWeight_vtx->GetErrorYlow(i);
     TrackCorr_vtx.push_back(val);
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

  double pileupWeight = 1, genWeight = 1, totalWeight = 1.0;
  double mu1TrackWeight=1.0, mu2TrackWeight=1.0, mu3TrackWeight=1.0, mu1IDWeight=1.0, mu2IDWeight=1.0, mu3IDWeight=1.0, mu1ISOWeight=1.0, mu2ISOWeight=1.0;
  double tauMedSF = 0.97, TriggerWeight = 1.0;
  float LumiFraction_GH = 16.1 / 35.9;
  float LumiFraction_BF = 19.8 / 35.9; //1.0-LumiFraction_GH;
 
  if (isMC_)
  {
    edm::Handle<std::vector<PileupSummaryInfo> > pPileupSummaryInfo;
    iEvent.getByToken(pileupSummaryInfo_, pPileupSummaryInfo);

    int nTrueVertices = 0;
    if (pPileupSummaryInfo.isValid() && pPileupSummaryInfo->size()>0)
      nTrueVertices = pPileupSummaryInfo->at(1).getTrueNumInteractions();

    PileupFile = new TFile(PileupFileName_.c_str());
    TH1F* Pileup_ = (TH1F*)PileupFile->Get("pileup_scale");
    pileupWeight = Pileup_->GetBinContent(nTrueVertices);
    PileupFile->Close();

    edm::Handle<GenEventInfoProduct> genEventInfo;
    iEvent.getByToken(genEventInfoToken_, genEventInfo);
    double eventGenWeight = genEventInfo->weight();
    genWeight = eventGenWeight * lumi_ * xsec_ / summedWeights_;
    std::cout << "\teventGenWeight=" << eventGenWeight << "\tlumi_=" << lumi_  << "\tsummedWeights_=" << summedWeights_ << "\txsec_=" << xsec_ << "\tnTrueVertices=" << nTrueVertices << std::endl;

  ///////////////////////////////////////
  // Do Track Iso and ID and ISO scale
  //////////////////////////////////////
    for (std::list<TrackProperties>::const_iterator it = TrackCorr_eta.begin(); it != TrackCorr_eta.end(); it++ )
    {
       std::cout << "fabs(mu1.eta()=" << fabs(mu1.eta()) << "   (*it).x-(*it).errx_down=" <<  (*it).x-(*it).errx_down << "  (*it).x+(*it).errx_up=" << (*it).x+(*it).errx_up << std::endl;
       if (fabs(mu1.eta()) >= (*it).x-(*it).errx_down && fabs(mu1.eta()) <= (*it).x+(*it).errx_up)
       {
          mu1TrackWeight *= (*it).y;
          break;
       }//if
    }//for it

    for (std::list<TrackProperties>::const_iterator it = TrackCorr_vtx.begin(); it != TrackCorr_vtx.end(); it++ )
    {
       std::cout << "nTrueVertices=" << nTrueVertices << "   (*it).x-(*it).errx_down=" <<  (*it).x-(*it).errx_down << "  (*it).x+(*it).errx_up=" << (*it).x+(*it).errx_up << std::endl;
       if (nTrueVertices >= (*it).x-(*it).errx_down && nTrueVertices <= (*it).x+(*it).errx_up)
       {
          mu1TrackWeight *= (*it).y;
          break;
       }//if
    }//for it

    for (std::list<TrackProperties>::const_iterator it = TrackCorr_eta.begin(); it != TrackCorr_eta.end(); it++ )
    {
       if (fabs(mu2.eta()) >= (*it).x-(*it).errx_down && fabs(mu2.eta()) <= (*it).x+(*it).errx_up)
       {
          mu2TrackWeight *= (*it).y;
          break;
       }//if
    }//for it

    for (std::list<TrackProperties>::const_iterator it = TrackCorr_vtx.begin(); it != TrackCorr_vtx.end(); it++ )
    {
       if (nTrueVertices >= (*it).x-(*it).errx_down && nTrueVertices <= (*it).x+(*it).errx_up)
       {
          mu2TrackWeight *= (*it).y;
          break;
       }//if
    }//for it

    float binxIDs_BF = IDsWeight_BToF->GetXaxis()->FindBin(mu1.pt());
    float binyIDs_BF = IDsWeight_BToF->GetYaxis()->FindBin(fabs(mu1.eta()));
    float binxIDs_GH = IDsWeight_GH->GetXaxis()->FindBin(mu1.pt());
    float binyIDs_GH = IDsWeight_GH->GetYaxis()->FindBin(fabs(mu1.eta()));
    float binxISOs_BF = ISOsWeight_BToF->GetXaxis()->FindBin(mu1.pt());
    float binyISOs_BF = ISOsWeight_BToF->GetYaxis()->FindBin(fabs(mu1.eta()));
    float binxISOs_GH = ISOsWeight_GH->GetXaxis()->FindBin(mu1.pt());
    float binyISOs_GH = ISOsWeight_GH->GetYaxis()->FindBin(fabs(mu1.eta()));
    if (binxIDs_BF == 7)
    {
      binxIDs_BF = 6;
      binxIDs_GH = 6;
      binxISOs_BF = 6;
      binxISOs_GH = 6;
    }//if
    mu1IDWeight  = LumiFraction_BF * IDsWeight_BToF->GetBinContent(binxIDs_BF, binyIDs_BF)    + LumiFraction_GH * IDsWeight_GH->GetBinContent(binxIDs_GH, binyIDs_GH);
    mu1ISOWeight = LumiFraction_BF * ISOsWeight_BToF->GetBinContent(binxISOs_BF, binyISOs_BF) + LumiFraction_GH * ISOsWeight_GH->GetBinContent(binxISOs_GH, binyISOs_GH);
  
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
      mu2IDWeight  = LumiFraction_BF * IDsWeight_BToF->GetBinContent(binxIDs_BF, binyIDs_BF)    + LumiFraction_GH * IDsWeight_GH->GetBinContent(binxIDs_GH, binyIDs_GH);
      mu2ISOWeight = LumiFraction_BF * ISOsWeight_BToF->GetBinContent(binxISOs_BF, binyISOs_BF) + LumiFraction_GH * ISOsWeight_GH->GetBinContent(binxISOs_GH, binyISOs_GH);
    }//if
    else
    {
      std::cout << "JPsi Eff" << std::endl;
      binyIDs_BF = IDsWeight_LowPt_ID->GetYaxis()->FindBin(mu2.pt());
      binxIDs_BF = IDsWeight_LowPt_ID->GetXaxis()->FindBin(fabs(mu2.eta()));
      binyISOs_BF = ISOsWeight_LowPt_ISO->GetYaxis()->FindBin(mu2.pt());
      binxISOs_BF = ISOsWeight_LowPt_ISO->GetXaxis()->FindBin(fabs(mu2.eta()));
      mu2IDWeight  = IDsWeight_BToF->GetBinContent(binxIDs_BF, binyIDs_BF);
      mu2ISOWeight = ISOsWeight_BToF->GetBinContent(binxISOs_BF, binyISOs_BF);
    }//else
  
  ///////////////////////////////////////
  // Trigger efficiencies
  ///////////////////////////////////////
    float binxTrigger_BF = TriggerWeight_BToF->GetXaxis()->FindBin(mu1.pt());
    float binyTrigger_BF = TriggerWeight_BToF->GetYaxis()->FindBin(fabs(mu1.eta()));
    float binxTrigger_GH = TriggerWeight_GH->GetXaxis()->FindBin(mu1.pt());
    float binyTrigger_GH = TriggerWeight_GH->GetYaxis()->FindBin(fabs(mu1.eta()));
    TriggerWeight = LumiFraction_BF *  TriggerWeight_BToF->GetBinContent(binxTrigger_BF, binyTrigger_BF) + LumiFraction_GH *  TriggerWeight_GH->GetBinContent(binxTrigger_GH, binyTrigger_GH);
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
  double DMFind = -1, MedIso = -1, bTagValue = -1000, smallestdRJet = 1000000000000, csvSF = 1.0;
  std::cout << "pTaus->size()=" << pTaus->size() << "\tpMuons->size()=" << pMuons->size() << std::endl;
  for (edm::View<pat::Tau>::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); ++iTau)
  {
    std::cout << "\tiTau->pt()=" << iTau->pt() << std::endl;
    NEvents_->Fill(1);
    double bestMu3dR = 10000000;
    bool checkSubMu = false;
    pat::Muon mu3;
    std::cout << "iTau->tauID(medIsoTag_)=" << iTau->tauID(medIsoTag_) << "  DMFind=" << iTau->tauID(decayModeFindingTag_) << "  iTau->pt()=" << iTau->pt() << "  fabs(iTau->eta() )=" << fabs(iTau->eta() ) << "  deltaR(mu1,*iTau)=" << deltaR(mu1,*iTau) << "  deltaR(mu2,*iTau)=" << deltaR(mu2,*iTau) << std::endl;
    if (iTau->tauID(minMVARaw_) <= -0.5 || iTau->pt() <= tauPtCut_ || fabs(iTau->eta() ) >= 2.4 || deltaR(mu1,*iTau) <= 0.8 || deltaR(mu2,*iTau) <= 0.8)
      continue;
    for (edm::View<pat::Muon>::const_iterator iMu = pMuons->begin(); iMu != pMuons->end(); ++iMu)
    {
      std::cout << "\tiMu->pt()=" << iMu->pt() << std::endl;
      double currdR = deltaR(*iTau, *iMu), mu1dR = deltaR(mu1, *iMu), mu2dR = deltaR(mu2, *iMu);
      if (mu1dR < mu3dROverlapCut_ || mu2dR < mu3dROverlapCut_)
       continue;
      if (currdR < mu3dRMax_ && currdR > mu3dRMin_ && currdR < bestMu3dR)
      {
        diTauP4 = iTau->p4() + iMu->p4();
        checkSubMu = true;
        mu3 = *iMu;
        DMFind = iTau->tauID(decayModeFindingTag_);
        MedIso = iTau->tauID(medIsoTau_);
        //mu Track Weight
        for (std::list<TrackProperties>::const_iterator it = TrackCorr_eta.begin(); it != TrackCorr_eta.end(); it++ )
        {
           if (fabs(mu3.eta()) >= (*it).x-(*it).errx_down && fabs(mu3.eta()) <= (*it).x+(*it).errx_up)
           {
              mu3TrackWeight *= (*it).y;
              break;
           }//if
        }//for it

        for (std::list<TrackProperties>::const_iterator it = TrackCorr_vtx.begin(); it != TrackCorr_vtx.end(); it++ )
        {
           if (fabs(mu3.eta()) >= (*it).x-(*it).errx_down && fabs(mu3.eta()) <= (*it).x+(*it).errx_up)
           {
              mu3TrackWeight *= (*it).y;
              break;
           }//if
        }//for it
        //mu3 ID weight
        if (mu3.pt() > 20.0)
        {
          std::cout << "Z Eff" << std::endl;
          double binxIDs_BF = IDsWeight_BToF->GetXaxis()->FindBin(mu3.pt());
          double binyIDs_BF = IDsWeight_BToF->GetYaxis()->FindBin(fabs(mu3.eta()));
          double binxIDs_GH = IDsWeight_GH->GetXaxis()->FindBin(mu3.pt());
          double binyIDs_GH = IDsWeight_GH->GetYaxis()->FindBin(fabs(mu3.eta()));
          if (binxIDs_BF == 7)
          {
            binxIDs_BF = 6;
            binxIDs_GH = 6;
          }//if
          mu3IDWeight  = LumiFraction_BF * IDsWeight_BToF->GetBinContent(binxIDs_BF, binyIDs_BF)    + LumiFraction_GH * IDsWeight_GH->GetBinContent(binxIDs_GH, binyIDs_GH);
        }//if
        else
        {
          std::cout << "JPsi Eff" << std::endl;
          double binyIDs_BF = IDsWeight_LowPt_ID->GetYaxis()->FindBin(mu3.pt());
          double binxIDs_BF = IDsWeight_LowPt_ID->GetXaxis()->FindBin(fabs(mu3.eta()));
          mu3IDWeight  = IDsWeight_LowPt_ID->GetBinContent(binxIDs_BF, binyIDs_BF);
        }//else
 
        //csv Weight
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
      }//if
    }//for iMu

    if (checkSubMu)
      NEvents_->Fill(2);
    else if (iTau->pt() > tauPtCut_ )
      NEvents_->Fill(3);
    else if ( bTagValue >= 0.8484)
      NEvents_->Fill(4);
    else if (DMFind >= .5)
      NEvents_->Fill(5);
    else if (MedIso >= .5)
      NEvents_->Fill(6);

    if (!checkSubMu && requireRemovedMuon_) 
      continue;

    std::cout << "Found DiTau\tbTagValue=" << bTagValue <<  std::endl;
    double csvWeight = 1.0;
    if (!isMC_ && bTagValue >= 0.8484)
      continue;
    else if (isMC_ && bTagValue >= 0.8484)
      csvWeight = 1 - csvSF;

    std::cout << "Pass CSV cut" << std::endl;
    ///////////////////////////////////////
    // Combining weights and SF
    ///////////////////////////////////////
    std::cout << "mu1ISOWeight=" << mu1ISOWeight << "  mu2ISOWeight=" << mu2ISOWeight << std::endl;
    std::cout << "mu1TrackWeight=" << mu1TrackWeight << "  mu2TrackWeight=" << mu2TrackWeight << "  mu3TrackWeight=" << mu3TrackWeight << std::endl;
    std::cout << "mu1IDWeight=" << mu1IDWeight << "  mu2IDWeight=" << mu2IDWeight << "  mu3IDWeight=" << mu3IDWeight << std::endl;
    std::cout << "pileupWeight=" << pileupWeight << "\tgenWeight=" << genWeight << "\tTriggerWeight=" << TriggerWeight << "tauMedSF=" << tauMedSF << "  csvWeight=" << csvWeight << std::endl;
    
    totalWeight = mu1ISOWeight * mu2ISOWeight * mu1TrackWeight * mu2TrackWeight * mu3TrackWeight * mu1IDWeight * mu2IDWeight * mu3IDWeight * TriggerWeight * pileupWeight * genWeight * tauMedSF * csvWeight;
    if (!isMC_)
      totalWeight = 1.0;

    std::cout << "total_weight=" << totalWeight << std::endl;
    if (checkSubMu)
    {
      reco::LeafCandidate::LorentzVector diMuP4_Mu1TauMu, diMuP4_Mu2TauMu;    
      diMuP4_Mu1TauMu = mu1.p4();
      diMuP4_Mu1TauMu += mu3.p4();
 
      diMuP4_Mu2TauMu = mu2.p4();
      diMuP4_Mu2TauMu += mu3.p4();

      InvMassTauMuMu1_->Fill(diMuP4_Mu1TauMu.M() ,totalWeight );
      InvMassTauMuMu2_->Fill(diMuP4_Mu2TauMu.M() ,totalWeight );
      MVARawvsMVADisc_->Fill(iTau->tauID(minMVARaw_), iTau->tauID(medIsoTau_), totalWeight );
    }//if removed Mu

    if (DMFind >= .5)
    {
      TauDMFindPt_->Fill(iTau->pt() ,totalWeight );
      TauDMFindEta_->Fill(iTau->eta() ,totalWeight );
      EtavsPtTauDMFind_->Fill(iTau->pt(), iTau->eta() ,totalWeight );
    }//if DMFind >= .5

    if (MedIso >= .5 && DMFind >= .5)
    {
      TauMedIsoPt_->Fill(iTau->pt() ,totalWeight );
      TauMedIsoEta_->Fill(iTau->eta() ,totalWeight );
      EtavsPtTauMedIso_->Fill(iTau->pt(), iTau->eta() ,totalWeight );
    }//if MedIso == 1 && DMFind >= .5

  }//iTau

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
      NEvents_->GetXaxis()->SetBinLabel(4, "#tau's p_{T} > 10");
      NEvents_->GetXaxis()->SetBinLabel(5, "Pass csv");
      NEvents_->GetXaxis()->SetBinLabel(6, "Pass DMFind");
      NEvents_->GetXaxis()->SetBinLabel(7, "Pass Med Iso");
      NEvents_->GetXaxis()->SetBinLabel(8, "Pass Tight Iso");
  InvMassTauMuMu1_     = new TH1F("InvMassTauMuMu1"    , "", 75, 0, 150);
  InvMassMu1Mu2_     = new TH1F("InvMassMu1Mu2"    , "", 75, 0, 150);
  InvMassTauMuMu2_     = new TH1F("InvMassTauMuMu2"    , "", 75, 0, 150);
  MVARawvsMVADisc_     = new TH2F("MVARawvsMVADisc"    , "", 75, 0, 150, 2, -0.5, 1.5);

  Float_t binsx[] = {9.9999, 20, 30, 40, 60, 300};
  Float_t binsy[] = {0, 0.9, 1.5, 2.4};
  EtavsPtTauMedIso_  = new TH2F("EtavsPtTauMedIso"     , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx, sizeof(binsy)/sizeof(Float_t) - 1, binsy);
  EtavsPtTauDMFind_  = new TH2F("EtavsPtTauDMFind"     , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx, sizeof(binsy)/sizeof(Float_t) - 1, binsy);

  TauMedIsoEta_    = new TH1F("TauMedIsoEta", "", 11, -2.4, 2.4);
  TauDMFindEta_    = new TH1F("TauDMFindEta"    , "", 11, -2.4, 2.4);

  TauMedIsoPt_    = new TH1F("TauMedIsoPt", "", sizeof(binsx)/sizeof(Float_t) - 1, binsx);
  TauDMFindPt_    = new TH1F("TauDMFindPt"    , "", sizeof(binsx)/sizeof(Float_t) - 1, binsx);

  InvMassTauMuMu1_->Sumw2();
  InvMassMu1Mu2_->Sumw2();
  InvMassTauMuMu2_->Sumw2();
  MVARawvsMVADisc_->Sumw2();

  EtavsPtTauMedIso_->Sumw2();
  EtavsPtTauDMFind_->Sumw2();

  TauMedIsoEta_->Sumw2();
  TauDMFindEta_->Sumw2();

  TauMedIsoPt_->Sumw2();
  TauDMFindPt_->Sumw2();
}

// ------------ method called once each job just after ending the event loop  ------------
void FakeRateMiniAODGetRates::endJob()
{
  //Make the Canvases
  TCanvas NEventsCanvas("NEvents","",600,600);
  TCanvas InvMassTauMuMu1Canvas("InvMassTauMuMu1","",600,600);
  TCanvas InvMassMu1Mu2Canvas("InvMassMu1Mu2","",600,600);
  TCanvas InvMassTauMuMu2Canvas("InvMassTauMuMu2","",600,600);
  TCanvas MVARawvsMVADiscCanvas("MVARawvsMVADisc","",600,600);

  TCanvas EtavsPtTauMedIsoCanvas("EtavsPtTauMedIso","",600,600);
  TCanvas EtavsPtTauDMFindCanvas("EtavsPtTauDMFind","",600,600);

  TCanvas TauMedIsoEtaCanvas("TauMedIsoEta","",600,600);
  TCanvas TauDMFindEtaCanvas("TauDMFindEta","",600,600);

  TCanvas TauMedIsoPtCanvas("TauMedIsoPt","",600,600);
  TCanvas TauDMFindPtCanvas("TauDMFindPt","",600,600);

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
  VariousFunctions::formatAndDrawCanvasAndHist2D(MVARawvsMVADiscCanvas, MVARawvsMVADisc_,
         1, 0, 0, kBlack, 7, 20, "MVA Raw Iso ", .04, .04, 1.1, "MVA Discriminator", .04, .04, 1.6, "", .04, .04, 1.0);


  VariousFunctions::formatAndDrawCanvasAndHist2D(EtavsPtTauMedIsoCanvas, EtavsPtTauMedIso_,
         1, 0, 0, kBlack, 7, 20, "p_{T} #tau Med Iso", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(EtavsPtTauDMFindCanvas, EtavsPtTauDMFind_,
         1, 0, 0, kBlack, 7, 20, "p_{T} #tau DMFind", .04, .04, 1.1, "nCleanJets #tau's", .04, .04, 1.6, "", .04, .04, 1.0);

  VariousFunctions::formatAndDrawCanvasAndHist1D(TauMedIsoEtaCanvas, TauMedIsoEta_,
	 1, 0, 0, kBlack, 7, 20, "Eta(PFTau + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauDMFindEtaCanvas, TauDMFindEta_,
	 1, 0, 0, kBlack, 7, 20, "Eta(PFTau + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);

  VariousFunctions::formatAndDrawCanvasAndHist1D(TauMedIsoPtCanvas, TauMedIsoPt_,
         1, 0, 0, kBlack, 7, 20, "Pt(PFTau + Med Iso + DM)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  VariousFunctions::formatAndDrawCanvasAndHist1D(TauDMFindPtCanvas, TauDMFindPt_,
         1, 0, 0, kBlack, 7, 20, "Pt(PFTau + DecayModeFinding)", .04, .04, 1.1,  "", .04, .04, 1.0, false);
  

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
  MVARawvsMVADisc_->Write();
 
  EtavsPtTauMedIso_->Write();
  EtavsPtTauDMFind_->Write();

  TauMedIsoEta_->Write();
  TauDMFindEta_->Write();

  TauMedIsoPt_->Write();
  TauDMFindPt_->Write();

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
  if ((doDelete) && (MVARawvsMVADisc_ != NULL)) delete MVARawvsMVADisc_;
  MVARawvsMVADisc_ = NULL;

  if ((doDelete) && (EtavsPtTauMedIso_ != NULL)) delete EtavsPtTauMedIso_;
  EtavsPtTauMedIso_ = NULL;
  if ((doDelete) && (EtavsPtTauDMFind_ != NULL)) delete EtavsPtTauDMFind_;
  EtavsPtTauDMFind_ = NULL;

  if ((doDelete) && (TauMedIsoEta_ != NULL)) delete TauMedIsoEta_;
  TauMedIsoEta_ = NULL;
  if ((doDelete) && (TauDMFindEta_ != NULL)) delete TauDMFindEta_;
  TauDMFindEta_ = NULL;

  if ((doDelete) && (TauMedIsoPt_ != NULL)) delete TauMedIsoPt_;
  TauMedIsoPt_ = NULL;
  if ((doDelete) && (TauDMFindPt_ != NULL)) delete TauDMFindPt_;
  TauDMFindPt_ = NULL;

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

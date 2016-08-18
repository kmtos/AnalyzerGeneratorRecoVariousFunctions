#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include "TCanvas.h"

void rootMacro_ZTT()
{
  gStyle->SetOptStat(kFALSE);

  TFile infileCJ("/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/ggH125a9_OLDDecayID_NewDMFind_dRchecks20_20PtCut_ForceEventMatch_JUN6/ggH125a9_OLDDecayID_NewDMFind_dRchecks20_20PtCut_ForceEventMatch_JUN6_Plots.root");
  TFile infileZTT("/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/ZTT_PU_JUN10/ZTT_PU_JUN10_Plots.root");
  TFile *outFile = new TFile("combHist_PlotFormatting_MAY4.root", "RECREATE");

cout << "Files Created" << endl;

  TCanvas *MatchedLooseIsoCJPtCanvas   = (TCanvas*)infileCJ.Get("MatchedLooseIsoCJPt");
  TCanvas *MatchedLooseIsoRECOPtCanvas = (TCanvas*)infileCJ.Get("MatchedLooseIsoRECOPt");
  TCanvas *MatchedMedIsoCJPtCanvas     = (TCanvas*)infileCJ.Get("MatchedMedIsoCJPt");
  TCanvas *MatchedMedIsoRECOPtCanvas   = (TCanvas*)infileCJ.Get("MatchedMedIsoRECOPt");
  TCanvas *MatchedTightIsoCJPtCanvas   = (TCanvas*)infileCJ.Get("MatchedTightIsoCJPt");
  TCanvas *MatchedTightIsoRECOPtCanvas = (TCanvas*)infileCJ.Get("MatchedTightIsoRECOPt");
  TCanvas *MatchedDMFindCJPtCanvas     = (TCanvas*)infileCJ.Get("MatchedDMFindCJPt");
  TCanvas *MatchedDMFindRECOPtCanvas   = (TCanvas*)infileCJ.Get("MatchedDMFindRECOPt");
  TCanvas *MatchedCJPtCanvas   = (TCanvas*)infileCJ.Get("MatchedCJPt");
  TCanvas *MatchedRECOPtCanvas = (TCanvas*)infileCJ.Get("MatchedRECOPt");

  TCanvas *MatchedLooseIsoCJdRCanvas   = (TCanvas*)infileCJ.Get("MatchedLooseIsoCJdR");
  TCanvas *MatchedLooseIsoRECOdRCanvas = (TCanvas*)infileCJ.Get("MatchedLooseIsoRECOdR");
  TCanvas *MatchedMedIsoCJdRCanvas     = (TCanvas*)infileCJ.Get("MatchedMedIsoCJdR");
  TCanvas *MatchedMedIsoRECOdRCanvas   = (TCanvas*)infileCJ.Get("MatchedMedIsoRECOdR");
  TCanvas *MatchedTightIsoCJdRCanvas   = (TCanvas*)infileCJ.Get("MatchedTightIsoCJdR");
  TCanvas *MatchedTightIsoRECOdRCanvas = (TCanvas*)infileCJ.Get("MatchedTightIsoRECOdR");
  TCanvas *MatchedDMFindCJdRCanvas     = (TCanvas*)infileCJ.Get("MatchedDMFindCJdR");
  TCanvas *MatchedDMFindRECOdRCanvas   = (TCanvas*)infileCJ.Get("MatchedDMFindRECOdR");
  TCanvas *MatchedCJdRCanvas   = (TCanvas*)infileCJ.Get("MatchedCJdR");
  TCanvas *MatchedRECOdRCanvas = (TCanvas*)infileCJ.Get("MatchedRECOdR");

  TCanvas *MatchedLooseIsoCJPtGenCanvas   = (TCanvas*)infileCJ.Get("MatchedLooseIsoCJPtGen");
  TCanvas *MatchedLooseIsoRECOPtGenCanvas = (TCanvas*)infileCJ.Get("MatchedLooseIsoRECOPtGen");
  TCanvas *MatchedMedIsoCJPtGenCanvas     = (TCanvas*)infileCJ.Get("MatchedMedIsoCJPtGen");
  TCanvas *MatchedMedIsoRECOPtGenCanvas   = (TCanvas*)infileCJ.Get("MatchedMedIsoRECOPtGen");
  TCanvas *MatchedTightIsoCJPtGenCanvas   = (TCanvas*)infileCJ.Get("MatchedTightIsoCJPtGen");
  TCanvas *MatchedTightIsoRECOPtGenCanvas = (TCanvas*)infileCJ.Get("MatchedTightIsoRECOPtGen");
  TCanvas *MatchedDMFindCJPtGenCanvas     = (TCanvas*)infileCJ.Get("MatchedDMFindCJPtGen");
  TCanvas *MatchedDMFindRECOPtGenCanvas   = (TCanvas*)infileCJ.Get("MatchedDMFindRECOPtGen");
  TCanvas *MatchedCJPtGenCanvas   = (TCanvas*)infileCJ.Get("MatchedCJPtGen");
  TCanvas *MatchedRECOPtGenCanvas = (TCanvas*)infileCJ.Get("MatchedRECOPtGen");


  TCanvas *MatchedLooseIsoRECOPtCanvasZTT   = (TCanvas*)infileZTT.Get("MatchedLooseIsoRECOPt");
  TCanvas *MatchedMedIsoRECOPtCanvasZTT     = (TCanvas*)infileZTT.Get("MatchedMedIsoRECOPt");
  TCanvas *MatchedTightIsoRECOPtCanvasZTT   = (TCanvas*)infileZTT.Get("MatchedTightIsoRECOPt");
  TCanvas *MatchedDMFindRECOPtCanvasZTT     = (TCanvas*)infileZTT.Get("MatchedDMFindRECOPt");
  TCanvas *MatchedRECOPtCanvasZTT   = (TCanvas*)infileZTT.Get("MatchedRECOPt");

  TCanvas *MatchedLooseIsoRECOdRCanvasZTT   = (TCanvas*)infileZTT.Get("MatchedLooseIsoRECOdR");
  TCanvas *MatchedMedIsoRECOdRCanvasZTT     = (TCanvas*)infileZTT.Get("MatchedMedIsoRECOdR");
  TCanvas *MatchedTightIsoRECOdRCanvasZTT   = (TCanvas*)infileZTT.Get("MatchedTightIsoRECOdR");
  TCanvas *MatchedDMFindRECOdRCanvasZTT     = (TCanvas*)infileZTT.Get("MatchedDMFindRECOdR");
  TCanvas *MatchedRECOdRCanvasZTT   = (TCanvas*)infileZTT.Get("MatchedRECOdR");
  
  TCanvas *MatchedLooseIsoRECOPtGenCanvasZTT   = (TCanvas*)infileZTT.Get("MatchedLooseIsoRECOPtGen");
  TCanvas *MatchedMedIsoRECOPtGenCanvasZTT     = (TCanvas*)infileZTT.Get("MatchedMedIsoRECOPtGen");
  TCanvas *MatchedTightIsoRECOPtGenCanvasZTT   = (TCanvas*)infileZTT.Get("MatchedTightIsoRECOPtGen");
  TCanvas *MatchedDMFindRECOPtGenCanvasZTT     = (TCanvas*)infileZTT.Get("MatchedDMFindRECOPtGen");
  TCanvas *MatchedRECOPtGenCanvasZTT   = (TCanvas*)infileZTT.Get("MatchedRECOPtGen");


cout << "Got Canvases" << endl;

  TH1F* MatchedLooseIsoCJPt_   = (TH1F*)MatchedLooseIsoCJPtCanvas->GetPrimitive("MatchedLooseIsoCJPt");
  TH1F* MatchedLooseIsoRECOPt_ = (TH1F*)MatchedLooseIsoRECOPtCanvas->GetPrimitive("MatchedLooseIsoRECOPt");
  TH1F* MatchedMedIsoCJPt_     = (TH1F*)MatchedMedIsoCJPtCanvas->GetPrimitive("MatchedMedIsoCJPt");
  TH1F* MatchedMedIsoRECOPt_   = (TH1F*)MatchedMedIsoRECOPtCanvas->GetPrimitive("MatchedMedIsoRECOPt");
  TH1F* MatchedTightIsoCJPt_   = (TH1F*)MatchedTightIsoCJPtCanvas->GetPrimitive("MatchedTightIsoCJPt");
  TH1F* MatchedTightIsoRECOPt_ = (TH1F*)MatchedTightIsoRECOPtCanvas->GetPrimitive("MatchedTightIsoRECOPt");
  TH1F* MatchedDMFindCJPt_     = (TH1F*)MatchedDMFindCJPtCanvas->GetPrimitive("MatchedDMFindCJPt");
  TH1F* MatchedDMFindRECOPt_   = (TH1F*)MatchedDMFindRECOPtCanvas->GetPrimitive("MatchedDMFindRECOPt");
  TH1F* MatchedCJPt_   = (TH1F*)MatchedCJPtCanvas->GetPrimitive("MatchedCJPt");
  TH1F* MatchedRECOPt_ = (TH1F*)MatchedRECOPtCanvas->GetPrimitive("MatchedRECOPt");

  TH1F* MatchedLooseIsoCJdR_   = (TH1F*)MatchedLooseIsoCJdRCanvas->GetPrimitive("MatchedLooseIsoCJdR");
  TH1F* MatchedLooseIsoRECOdR_ = (TH1F*)MatchedLooseIsoRECOdRCanvas->GetPrimitive("MatchedLooseIsoRECOdR");
  TH1F* MatchedMedIsoCJdR_     = (TH1F*)MatchedMedIsoCJdRCanvas->GetPrimitive("MatchedMedIsoCJdR");
  TH1F* MatchedMedIsoRECOdR_   = (TH1F*)MatchedMedIsoRECOdRCanvas->GetPrimitive("MatchedMedIsoRECOdR");
  TH1F* MatchedTightIsoCJdR_   = (TH1F*)MatchedTightIsoCJdRCanvas->GetPrimitive("MatchedTightIsoCJdR");
  TH1F* MatchedTightIsoRECOdR_ = (TH1F*)MatchedTightIsoRECOdRCanvas->GetPrimitive("MatchedTightIsoRECOdR");
  TH1F* MatchedDMFindCJdR_     = (TH1F*)MatchedDMFindCJdRCanvas->GetPrimitive("MatchedDMFindCJdR");
  TH1F* MatchedDMFindRECOdR_   = (TH1F*)MatchedDMFindRECOdRCanvas->GetPrimitive("MatchedDMFindRECOdR");
  TH1F* MatchedCJdR_   = (TH1F*)MatchedCJdRCanvas->GetPrimitive("MatchedCJdR");
  TH1F* MatchedRECOdR_ = (TH1F*)MatchedRECOdRCanvas->GetPrimitive("MatchedRECOdR");

  TH1F* MatchedLooseIsoCJPtGen_   = (TH1F*)MatchedLooseIsoCJPtGenCanvas->GetPrimitive("MatchedLooseIsoCJPtGen");
  TH1F* MatchedLooseIsoRECOPtGen_ = (TH1F*)MatchedLooseIsoRECOPtGenCanvas->GetPrimitive("MatchedLooseIsoRECOPtGen");
  TH1F* MatchedMedIsoCJPtGen_     = (TH1F*)MatchedMedIsoCJPtGenCanvas->GetPrimitive("MatchedMedIsoCJPtGen");
  TH1F* MatchedMedIsoRECOPtGen_   = (TH1F*)MatchedMedIsoRECOPtGenCanvas->GetPrimitive("MatchedMedIsoRECOPtGen");
  TH1F* MatchedTightIsoCJPtGen_   = (TH1F*)MatchedTightIsoCJPtGenCanvas->GetPrimitive("MatchedTightIsoCJPtGen");
  TH1F* MatchedTightIsoRECOPtGen_ = (TH1F*)MatchedTightIsoRECOPtGenCanvas->GetPrimitive("MatchedTightIsoRECOPtGen");
  TH1F* MatchedDMFindCJPtGen_     = (TH1F*)MatchedDMFindCJPtGenCanvas->GetPrimitive("MatchedDMFindCJPtGen");
  TH1F* MatchedDMFindRECOPtGen_   = (TH1F*)MatchedDMFindRECOPtGenCanvas->GetPrimitive("MatchedDMFindRECOPtGen");
  TH1F* MatchedCJPtGen_   = (TH1F*)MatchedCJPtGenCanvas->GetPrimitive("MatchedCJPtGen");
  TH1F* MatchedRECOPtGen_ = (TH1F*)MatchedRECOPtGenCanvas->GetPrimitive("MatchedRECOPtGen");


  TH1F* MatchedLooseIsoRECOPtZTT_ = (TH1F*)MatchedLooseIsoRECOPtCanvasZTT->GetPrimitive("MatchedLooseIsoRECOPt");
  TH1F* MatchedMedIsoRECOPtZTT_   = (TH1F*)MatchedMedIsoRECOPtCanvasZTT->GetPrimitive("MatchedMedIsoRECOPt");
  TH1F* MatchedTightIsoRECOPtZTT_ = (TH1F*)MatchedTightIsoRECOPtCanvasZTT->GetPrimitive("MatchedTightIsoRECOPt");
  TH1F* MatchedDMFindRECOPtZTT_   = (TH1F*)MatchedDMFindRECOPtCanvasZTT->GetPrimitive("MatchedDMFindRECOPt");
  TH1F* MatchedRECOPtZTT_ = (TH1F*)MatchedRECOPtCanvasZTT->GetPrimitive("MatchedRECOPt");

  TH1F* MatchedLooseIsoRECOdRZTT_ = (TH1F*)MatchedLooseIsoRECOdRCanvasZTT->GetPrimitive("MatchedLooseIsoRECOdR");
  TH1F* MatchedMedIsoRECOdRZTT_   = (TH1F*)MatchedMedIsoRECOdRCanvasZTT->GetPrimitive("MatchedMedIsoRECOdR");
  TH1F* MatchedTightIsoRECOdRZTT_ = (TH1F*)MatchedTightIsoRECOdRCanvasZTT->GetPrimitive("MatchedTightIsoRECOdR");
  TH1F* MatchedDMFindRECOdRZTT_   = (TH1F*)MatchedDMFindRECOdRCanvasZTT->GetPrimitive("MatchedDMFindRECOdR");
  TH1F* MatchedRECOdRZTT_ = (TH1F*)MatchedRECOdRCanvasZTT->GetPrimitive("MatchedRECOdR");

  TH1F* MatchedLooseIsoRECOPtGenZTT_ = (TH1F*)MatchedLooseIsoRECOPtGenCanvasZTT->GetPrimitive("MatchedLooseIsoRECOPtGen");
  TH1F* MatchedMedIsoRECOPtGenZTT_   = (TH1F*)MatchedMedIsoRECOPtGenCanvasZTT->GetPrimitive("MatchedMedIsoRECOPtGen");
  TH1F* MatchedTightIsoRECOPtGenZTT_ = (TH1F*)MatchedTightIsoRECOPtGenCanvasZTT->GetPrimitive("MatchedTightIsoRECOPtGen");
  TH1F* MatchedDMFindRECOPtGenZTT_   = (TH1F*)MatchedDMFindRECOPtGenCanvasZTT->GetPrimitive("MatchedDMFindRECOPtGen");
  TH1F* MatchedRECOPtGenZTT_ = (TH1F*)MatchedRECOPtGenCanvasZTT->GetPrimitive("MatchedRECOPtGen");


cout << "Got Histos" << endl;

cout << "Histograms assigned." << endl;

  TCanvas LoosePtCanvas("LoosePtCanvas","",600,600);
  TCanvas MedPtCanvas("MedPtCanvas","",600,600);
  TCanvas TightPtCanvas("TightPtCanvas","",600,600);
  TCanvas DMPtCanvas("DMPtCanvas","",600,600);

  TCanvas LoosedRCanvas("LoosedRCanvas","",600,600);
  TCanvas MeddRCanvas("MeddRCanvas","",600,600);
  TCanvas TightdRCanvas("TightdRCanvas","",600,600);
  TCanvas DMdRCanvas("DMdRCanvas","",600,600);

  TCanvas LoosePtGenCanvas("LoosePtGenCanvas","",600,600);
  TCanvas MedPtGenCanvas("MedPtGenCanvas","",600,600);
  TCanvas TightPtGenCanvas("TightPtGenCanvas","",600,600);
  TCanvas DMPtGenCanvas("DMPtGenCanvas","",600,600);

cout << "Canvases created" << endl;

  // To divide histograms for efficiency with new Gen Matching Old DM's
  TGraphAsymmErrors* FinalEffLooseIsoRECOPt_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoRECOPt_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoRECOPt_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindRECOPt_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoCJPt_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoCJPt_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoCJPt_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindCJPt_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoRECOdR_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoRECOdR_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoRECOdR_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindRECOdR_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoCJdR_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoCJdR_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoCJdR_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindCJdR_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoRECOPtGen_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoRECOPtGen_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoRECOPtGen_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindRECOPtGen_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoCJPtGen_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoCJPtGen_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoCJPtGen_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindCJPtGen_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoCJdRZTT_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoCJdRZTT_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoCJdRZTT_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindCJdRZTT_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoRECOPtGenZTT_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoRECOPtGenZTT_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoRECOPtGenZTT_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindRECOPtGenZTT_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoCJPtGenZTT_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoCJPtGenZTT_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoCJPtGenZTT_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindCJPtGenZTT_ = new TGraphAsymmErrors(30);


  // To divide histograms for efficiency with new Gen Matching Old DM's
  FinalEffLooseIsoRECOPt_->Divide(MatchedLooseIsoRECOPt_, MatchedRECOPt_);
  FinalEffMedIsoRECOPt_->Divide(MatchedMedIsoRECOPt_,     MatchedRECOPt_);
  FinalEffTightIsoRECOPt_->Divide(MatchedTightIsoRECOPt_, MatchedRECOPt_);
  FinalEffDMFindRECOPt_->Divide(MatchedDMFindRECOPt_,     MatchedRECOPt_);
  FinalEffLooseIsoCJPt_->Divide(MatchedLooseIsoCJPt_, MatchedCJPt_);
  FinalEffMedIsoCJPt_->Divide(MatchedMedIsoCJPt_,     MatchedCJPt_);
  FinalEffTightIsoCJPt_->Divide(MatchedTightIsoCJPt_, MatchedCJPt_);
  FinalEffDMFindCJPt_->Divide(MatchedDMFindCJPt_,     MatchedCJPt_);

  FinalEffLooseIsoRECOdR_->Divide(MatchedLooseIsoRECOdR_, MatchedRECOdR_);
  FinalEffMedIsoRECOdR_->Divide(MatchedMedIsoRECOdR_,     MatchedRECOdR_);
  FinalEffTightIsoRECOdR_->Divide(MatchedTightIsoRECOdR_, MatchedRECOdR_);
  FinalEffDMFindRECOdR_->Divide(MatchedDMFindRECOdR_,     MatchedRECOdR_);
  FinalEffLooseIsoCJdR_->Divide(MatchedLooseIsoCJdR_, MatchedCJdR_);
  FinalEffMedIsoCJdR_->Divide(MatchedMedIsoCJdR_,     MatchedCJdR_);
  FinalEffTightIsoCJdR_->Divide(MatchedTightIsoCJdR_, MatchedCJdR_);
  FinalEffDMFindCJdR_->Divide(MatchedDMFindCJdR_,     MatchedCJdR_);

  FinalEffLooseIsoRECOPtGen_->Divide(MatchedLooseIsoRECOPtGen_, MatchedRECOPtGen_);
  FinalEffMedIsoRECOPtGen_->Divide(MatchedMedIsoRECOPtGen_,     MatchedRECOPtGen_);
  FinalEffTightIsoRECOPtGen_->Divide(MatchedTightIsoRECOPtGen_, MatchedRECOPtGen_);
  FinalEffDMFindRECOPtGen_->Divide(MatchedDMFindRECOPtGen_,     MatchedRECOPtGen_);
  FinalEffLooseIsoCJPtGen_->Divide(MatchedLooseIsoCJPtGen_, MatchedCJPtGen_);
  FinalEffMedIsoCJPtGen_->Divide(MatchedMedIsoCJPtGen_,     MatchedCJPtGen_);
  FinalEffTightIsoCJPtGen_->Divide(MatchedTightIsoCJPtGen_, MatchedCJPtGen_);
  FinalEffDMFindCJPtGen_->Divide(MatchedDMFindCJPtGen_,     MatchedCJPtGen_);

  // To divide histograms for efficiency with new Gen Matching Old DM's
  FinalEffLooseIsoRECOPtZTT_->Divide(MatchedLooseIsoRECOPtZTT_, MatchedRECOPtZTT_);
  FinalEffMedIsoRECOPtZTT_->Divide(MatchedMedIsoRECOPtZTT_,     MatchedRECOPtZTT_);
  FinalEffTightIsoRECOPtZTT_->Divide(MatchedTightIsoRECOPtZTT_, MatchedRECOPtZTT_);
  FinalEffDMFindRECOPtZTT_->Divide(MatchedDMFindRECOPtZTT_,     MatchedRECOPtZTT_);

  FinalEffLooseIsoRECOdRZTT_->Divide(MatchedLooseIsoRECOdRZTT_, MatchedRECOdRZTT_);
  FinalEffMedIsoRECOdRZTT_->Divide(MatchedMedIsoRECOdRZTT_,     MatchedRECOdRZTT_);
  FinalEffTightIsoRECOdRZTT_->Divide(MatchedTightIsoRECOdRZTT_, MatchedRECOdRZTT_);
  FinalEffDMFindRECOdRZTT_->Divide(MatchedDMFindRECOdRZTT_,     MatchedRECOdRZTT_);

  FinalEffLooseIsoRECOPtGenZTT_->Divide(MatchedLooseIsoRECOPtGenZTT_, MatchedRECOPtGenZTT_);
  FinalEffMedIsoRECOPtGenZTT_->Divide(MatchedMedIsoRECOPtGenZTT_,     MatchedRECOPtGenZTT_);
  FinalEffTightIsoRECOPtGenZTT_->Divide(MatchedTightIsoRECOPtGenZTT_, MatchedRECOPtGenZTT_);
  FinalEffDMFindRECOPtGenZTT_->Divide(MatchedDMFindRECOPtGenZTT_,     MatchedRECOPtGenZTT_);

cout << "Attributes set." << endl;  

  //Set Colors for New Gen Matching Old DM's
  setGraphOptions(FinalEffLooseIsoCJPt_, kRed, 0.7, 10, "p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffLooseIsoRECOPt_, kBlack, 0.7, 10, "p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoCJPt_, kRed, 0.7, 10, "p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoRECOPt_, kBlack, 0.7, 10, "p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoCJPt_, kRed, 0.7, 10, "p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoRECOPt_, kBlack, 0.7, 10, "p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindCJPt_, kRed, 0.7, 10, "p_{T}", "#epsilon (DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindRECOPt_, kBlack, 0.7, 10, "p_{T}", "#epsilon (DMFinding + GM / GM)");

  setGraphOptions(FinalEffLooseIsoCJdR_, kRed, 0.7, 10, "#DeltaR", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffLooseIsoRECOdR_, kBlack, 0.7, 10, "#DeltaR", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoCJdR_, kRed, 0.7, 10, "#DeltaR", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoRECOdR_, kBlack, 0.7, 10, "#DeltaR", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoCJdR_, kRed, 0.7, 10, "#DeltaR", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoRECOdR_, kBlack, 0.7, 10, "#DeltaR", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindCJdR_, kRed, 0.7, 10, "#DeltaR", "#epsilon (DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindRECOdR_, kBlack, 0.7, 10, "#DeltaR", "#epsilon (DMFinding + GM / GM)");

  setGraphOptions(FinalEffLooseIsoCJPtGen_, kRed, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffLooseIsoRECOPtGen_, kBlack, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoCJPtGen_, kRed, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoRECOPtGen_, kBlack, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoCJPtGen_, kRed, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoRECOPtGen_, kBlack, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindCJPtGen_, kRed, 0.7, 10, "Gen Visible p_{T}", "#epsilon (DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindRECOPtGen_, kBlack, 0.7, 10, "Gen Visible p_{T}", "#epsilon (DMFinding + GM / GM)");

  //Set Colors for New Gen Matching Old DM's
  setGraphOptions(FinalEffLooseIsoRECOPtZTT_, kBlue, 0.7, 10, "p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoRECOPtZTT_, kBlue, 0.7, 10, "p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoRECOPtZTT_, kBlue, 0.7, 10, "p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindRECOPtZTT_, kBlue, 0.7, 10, "p_{T}", "#epsilon (DMFinding + GM / GM)");

  setGraphOptions(FinalEffLooseIsoRECOdRZTT_, kBlue, 0.7, 10, "#DeltaR", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoRECOdRZTT_, kBlue, 0.7, 10, "#DeltaR", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoRECOdRZTT_, kBlue, 0.7, 10, "#DeltaR", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindRECOdRZTT_, kBlue, 0.7, 10, "#DeltaR", "#epsilon (DMFinding + GM / GM)");

  setGraphOptions(FinalEffLooseIsoRECOPtGenZTT_, kBlue, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoRECOPtGenZTT_, kBlue, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoRECOPtGenZTT_, kBlue, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindRECOPtGenZTT_, kBlue, 0.7, 10, "Gen Visible p_{T}", "#epsilon (DMFinding + GM / GM)");

  LoosePtCanvas.SetGrid(1,1);
  MedPtCanvas.SetGrid(1,1);
  TightPtCanvas.SetGrid(1,1);
  DMPtCanvas.SetGrid(1,1);

  LoosedRCanvas.SetGrid(1,1);
  MeddRCanvas.SetGrid(1,1);
  TightdRCanvas.SetGrid(1,1);
  DMdRCanvas.SetGrid(1,1);

  LoosePtGenCanvas.SetGrid(1,1);
  MedPtGenCanvas.SetGrid(1,1);
  TightPtGenCanvas.SetGrid(1,1);
  DMPtGenCanvas.SetGrid(1,1);

cout << "Histograms Drawn" << endl;

  leg = new TLegend(0.1,0.7,0.25,0.9);
  leg->AddEntry(FinalEffDMFindRECOPtGen_, "No CJ","L");
  leg->AddEntry(FinalEffDMFindCJPtGen_, "CJ","L");
  leg->AddEntry(FinalEffLooseIsoRECOPtZTT_,, "ZTT","L");

  LoosePtCanvas.cd();
  FinalEffLooseIsoCJPt_->Draw();
  FinalEffLooseIsoRECOPt_->Draw("SAME");
  FinalEffLooseIsoCJPtZTT_->Draw("SAME");
  leg->Draw();

  MedPtCanvas.cd();
  FinalEffMedIsoCJPt_->Draw();
  FinalEffMedIsoRECOPt_->Draw("SAME");
  FinalEffMedIsoCJPtZTT_->Draw("SAME");
  leg->Draw();

  TightPtCanvas.cd();
  FinalEffTightIsoCJPt_->Draw();
  FinalEffTightIsoRECOPt_->Draw("SAME");
  FinalEffTightIsoCJPtZTT_->Draw("SAME");
  leg->Draw();

  DMPtCanvas.cd();
  FinalEffDMFindCJPt_->Draw();
  FinalEffDMFindRECOPt_->Draw("SAME");
  FinalEffDMFindCJPtZTT_->Draw("SAME");
  leg->Draw();

  outFile->cd();
  LoosePtCanvas.Write();
  MedPtCanvas.Write();
  TightPtCanvas.Write();
  DMPtCanvas.Write();

  LoosedRCanvas.Write();
  MeddRCanvas.Write();
  TightdRCanvas.Write();
  DMdRCanvas.Write();

  LoosePtGenCanvas.Write();
  MedPtGenCanvas.Write();
  TightPtGenCanvas.Write();
  DMPtGenCanvas.Write();
  outFile->Write();
  outFile->Close();

cout << "end" << endl;

}//rootMacro_BBA_combine

#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TMultiGraph.h"


void rootMacro_CleanJets_OnlyNewDMs()
{

  gStyle->SetOptStat(kFALSE);
  gStyle->SetEndErrorSize(7);
  TFile infileCJNew("/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/ggH125a9_GenTauDecayID_IndivCJ_NewDMFind_MVA_APR1/ggH125a9_GenTauDecayID_IndivCJ_NewDMFind_MVA_APR1_Plots.root");
  TFile infileRECONew("/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/ggH125a9_GenTauDecayID_IndivRECO_NewDMFind_MVA_APR1/ggH125a9_GenTauDecayID_IndivRECO_NewDMFind_MVA_APR1_Plots.root");
  TFile infileKyleNew("/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/ggH125a9_OLDDecayID_NewDMFind_MVA_APR1/ggH125a9_OLDDecayID_NewDMFind_MVA_APR1_Plots.root");

cout << "Getting Files" << std::endl;

  TFile *outFile = new TFile("combHist_CleanJets_h125a9_APR1.root", "RECREATE");

cout << "Files Created" << endl;

  // CJ/RECO Comparison histograms
  TCanvas *MatchedLooseIsoCJPtCanvasRachNew   = (TCanvas*)infileCJNew.Get("MatchedLooseIsoCJPt");
  TCanvas *MatchedLooseIsoRECOPtCanvasRachNew = (TCanvas*)infileRECONew.Get("MatchedLooseIsoRECOPt");
  TCanvas *MatchedMedIsoCJPtCanvasRachNew     = (TCanvas*)infileCJNew.Get("MatchedMedIsoCJPt");
  TCanvas *MatchedMedIsoRECOPtCanvasRachNew   = (TCanvas*)infileRECONew.Get("MatchedMedIsoRECOPt");
  TCanvas *MatchedTightIsoCJPtCanvasRachNew   = (TCanvas*)infileCJNew.Get("MatchedTightIsoCJPt");
  TCanvas *MatchedTightIsoRECOPtCanvasRachNew = (TCanvas*)infileRECONew.Get("MatchedTightIsoRECOPt");
  TCanvas *MatchedDMFindCJPtCanvasRachNew     = (TCanvas*)infileCJNew.Get("MatchedDMFindCJPt");
  TCanvas *MatchedDMFindRECOPtCanvasRachNew   = (TCanvas*)infileRECONew.Get("MatchedDMFindRECOPt");
  TCanvas *MatchedCJPtCanvasRachNew   = (TCanvas*)infileCJNew.Get("MatchedCJPt");
  TCanvas *MatchedRECOPtCanvasRachNew = (TCanvas*)infileRECONew.Get("MatchedRECOPt");
  
  TCanvas *MatchedLooseIsoCJdRCanvasRachNew   = (TCanvas*)infileCJNew.Get("MatchedLooseIsoCJdR");
  TCanvas *MatchedLooseIsoRECOdRCanvasRachNew = (TCanvas*)infileRECONew.Get("MatchedLooseIsoRECOdR");
  TCanvas *MatchedMedIsoCJdRCanvasRachNew     = (TCanvas*)infileCJNew.Get("MatchedMedIsoCJdR");
  TCanvas *MatchedMedIsoRECOdRCanvasRachNew   = (TCanvas*)infileRECONew.Get("MatchedMedIsoRECOdR");
  TCanvas *MatchedTightIsoCJdRCanvasRachNew   = (TCanvas*)infileCJNew.Get("MatchedTightIsoCJdR");
  TCanvas *MatchedTightIsoRECOdRCanvasRachNew = (TCanvas*)infileRECONew.Get("MatchedTightIsoRECOdR");
  TCanvas *MatchedDMFindCJdRCanvasRachNew     = (TCanvas*)infileCJNew.Get("MatchedDMFindCJdR");
  TCanvas *MatchedDMFindRECOdRCanvasRachNew   = (TCanvas*)infileRECONew.Get("MatchedDMFindRECOdR");
  TCanvas *MatchedCJdRCanvasRachNew   = (TCanvas*)infileCJNew.Get("MatchedCJdR");
  TCanvas *MatchedRECOdRCanvasRachNew = (TCanvas*)infileRECONew.Get("MatchedRECOdR");

  TCanvas *MatchedLooseIsoCJPtGenCanvasRachNew   = (TCanvas*)infileCJNew.Get("MatchedLooseIsoCJPtGen");
  TCanvas *MatchedLooseIsoRECOPtGenCanvasRachNew = (TCanvas*)infileRECONew.Get("MatchedLooseIsoRECOPtGen");
  TCanvas *MatchedMedIsoCJPtGenCanvasRachNew     = (TCanvas*)infileCJNew.Get("MatchedMedIsoCJPtGen");
  TCanvas *MatchedMedIsoRECOPtGenCanvasRachNew   = (TCanvas*)infileRECONew.Get("MatchedMedIsoRECOPtGen");
  TCanvas *MatchedTightIsoCJPtGenCanvasRachNew   = (TCanvas*)infileCJNew.Get("MatchedTightIsoCJPtGen");
  TCanvas *MatchedTightIsoRECOPtGenCanvasRachNew = (TCanvas*)infileRECONew.Get("MatchedTightIsoRECOPtGen");
  TCanvas *MatchedDMFindCJPtGenCanvasRachNew     = (TCanvas*)infileCJNew.Get("MatchedDMFindCJPtGen");
  TCanvas *MatchedDMFindRECOPtGenCanvasRachNew   = (TCanvas*)infileRECONew.Get("MatchedDMFindRECOPtGen");
  TCanvas *MatchedCJPtGenCanvasRachNew   = (TCanvas*)infileCJNew.Get("MatchedCJPtGen");
  TCanvas *MatchedRECOPtGenCanvasRachNew = (TCanvas*)infileRECONew.Get("MatchedRECOPtGen");

  // Kyle Gen Matching New DM 
  TCanvas *MatchedLooseIsoCJPtCanvasKyleNew   = (TCanvas*)infileKyleNew.Get("MatchedLooseIsoCJPt");
  TCanvas *MatchedLooseIsoRECOPtCanvasKyleNew = (TCanvas*)infileKyleNew.Get("MatchedLooseIsoRECOPt");
  TCanvas *MatchedMedIsoCJPtCanvasKyleNew     = (TCanvas*)infileKyleNew.Get("MatchedMedIsoCJPt");
  TCanvas *MatchedMedIsoRECOPtCanvasKyleNew   = (TCanvas*)infileKyleNew.Get("MatchedMedIsoRECOPt");
  TCanvas *MatchedTightIsoCJPtCanvasKyleNew   = (TCanvas*)infileKyleNew.Get("MatchedTightIsoCJPt");
  TCanvas *MatchedTightIsoRECOPtCanvasKyleNew = (TCanvas*)infileKyleNew.Get("MatchedTightIsoRECOPt");
  TCanvas *MatchedDMFindCJPtCanvasKyleNew     = (TCanvas*)infileKyleNew.Get("MatchedDMFindCJPt");
  TCanvas *MatchedDMFindRECOPtCanvasKyleNew   = (TCanvas*)infileKyleNew.Get("MatchedDMFindRECOPt");
  TCanvas *MatchedCJPtCanvasKyleNew   = (TCanvas*)infileKyleNew.Get("MatchedCJPt");
  TCanvas *MatchedRECOPtCanvasKyleNew = (TCanvas*)infileKyleNew.Get("MatchedRECOPt");

  TCanvas *MatchedLooseIsoCJdRCanvasKyleNew   = (TCanvas*)infileKyleNew.Get("MatchedLooseIsoCJdR");
  TCanvas *MatchedLooseIsoRECOdRCanvasKyleNew = (TCanvas*)infileKyleNew.Get("MatchedLooseIsoRECOdR");
  TCanvas *MatchedMedIsoCJdRCanvasKyleNew     = (TCanvas*)infileKyleNew.Get("MatchedMedIsoCJdR");
  TCanvas *MatchedMedIsoRECOdRCanvasKyleNew   = (TCanvas*)infileKyleNew.Get("MatchedMedIsoRECOdR");
  TCanvas *MatchedTightIsoCJdRCanvasKyleNew   = (TCanvas*)infileKyleNew.Get("MatchedTightIsoCJdR");
  TCanvas *MatchedTightIsoRECOdRCanvasKyleNew = (TCanvas*)infileKyleNew.Get("MatchedTightIsoRECOdR");
  TCanvas *MatchedDMFindCJdRCanvasKyleNew     = (TCanvas*)infileKyleNew.Get("MatchedDMFindCJdR");
  TCanvas *MatchedDMFindRECOdRCanvasKyleNew   = (TCanvas*)infileKyleNew.Get("MatchedDMFindRECOdR");
  TCanvas *MatchedCJdRCanvasKyleNew   = (TCanvas*)infileKyleNew.Get("MatchedCJdR");
  TCanvas *MatchedRECOdRCanvasKyleNew = (TCanvas*)infileKyleNew.Get("MatchedRECOdR");

  TCanvas *MatchedLooseIsoCJPtGenCanvasKyleNew   = (TCanvas*)infileKyleNew.Get("MatchedLooseIsoCJPtGen");
  TCanvas *MatchedLooseIsoRECOPtGenCanvasKyleNew = (TCanvas*)infileKyleNew.Get("MatchedLooseIsoRECOPtGen");
  TCanvas *MatchedMedIsoCJPtGenCanvasKyleNew     = (TCanvas*)infileKyleNew.Get("MatchedMedIsoCJPtGen");
  TCanvas *MatchedMedIsoRECOPtGenCanvasKyleNew   = (TCanvas*)infileKyleNew.Get("MatchedMedIsoRECOPtGen");
  TCanvas *MatchedTightIsoCJPtGenCanvasKyleNew   = (TCanvas*)infileKyleNew.Get("MatchedTightIsoCJPtGen");
  TCanvas *MatchedTightIsoRECOPtGenCanvasKyleNew = (TCanvas*)infileKyleNew.Get("MatchedTightIsoRECOPtGen");
  TCanvas *MatchedDMFindCJPtGenCanvasKyleNew     = (TCanvas*)infileKyleNew.Get("MatchedDMFindCJPtGen");
  TCanvas *MatchedDMFindRECOPtGenCanvasKyleNew   = (TCanvas*)infileKyleNew.Get("MatchedDMFindRECOPtGen");
  TCanvas *MatchedCJPtGenCanvasKyleNew   = (TCanvas*)infileKyleNew.Get("MatchedCJPtGen");
  TCanvas *MatchedRECOPtGenCanvasKyleNew = (TCanvas*)infileKyleNew.Get("MatchedRECOPtGen");

cout << "Got Canvases" << endl;

  // Getting the histograms for the New gen matching and DM New
  TH1F* MatchedLooseIsoCJPtRachNew_   = (TH1F*)MatchedLooseIsoCJPtCanvasRachNew->GetPrimitive("MatchedLooseIsoCJPt");
  TH1F* MatchedLooseIsoRECOPtRachNew_ = (TH1F*)MatchedLooseIsoRECOPtCanvasRachNew->GetPrimitive("MatchedLooseIsoRECOPt");
  TH1F* MatchedMedIsoCJPtRachNew_     = (TH1F*)MatchedMedIsoCJPtCanvasRachNew->GetPrimitive("MatchedMedIsoCJPt");
  TH1F* MatchedMedIsoRECOPtRachNew_   = (TH1F*)MatchedMedIsoRECOPtCanvasRachNew->GetPrimitive("MatchedMedIsoRECOPt");
  TH1F* MatchedTightIsoCJPtRachNew_   = (TH1F*)MatchedTightIsoCJPtCanvasRachNew->GetPrimitive("MatchedTightIsoCJPt");
  TH1F* MatchedTightIsoRECOPtRachNew_ = (TH1F*)MatchedTightIsoRECOPtCanvasRachNew->GetPrimitive("MatchedTightIsoRECOPt");
  TH1F* MatchedDMFindCJPtRachNew_     = (TH1F*)MatchedDMFindCJPtCanvasRachNew->GetPrimitive("MatchedDMFindCJPt");
  TH1F* MatchedDMFindRECOPtRachNew_   = (TH1F*)MatchedDMFindRECOPtCanvasRachNew->GetPrimitive("MatchedDMFindRECOPt");
  TH1F* MatchedCJPtRachNew_   = (TH1F*)MatchedCJPtCanvasRachNew->GetPrimitive("MatchedCJPt");
  TH1F* MatchedRECOPtRachNew_ = (TH1F*)MatchedRECOPtCanvasRachNew->GetPrimitive("MatchedRECOPt");

  TH1F* MatchedLooseIsoCJdRRachNew_   = (TH1F*)MatchedLooseIsoCJdRCanvasRachNew->GetPrimitive("MatchedLooseIsoCJdR");
  TH1F* MatchedLooseIsoRECOdRRachNew_ = (TH1F*)MatchedLooseIsoRECOdRCanvasRachNew->GetPrimitive("MatchedLooseIsoRECOdR");
  TH1F* MatchedMedIsoCJdRRachNew_     = (TH1F*)MatchedMedIsoCJdRCanvasRachNew->GetPrimitive("MatchedMedIsoCJdR");
  TH1F* MatchedMedIsoRECOdRRachNew_   = (TH1F*)MatchedMedIsoRECOdRCanvasRachNew->GetPrimitive("MatchedMedIsoRECOdR");
  TH1F* MatchedTightIsoCJdRRachNew_   = (TH1F*)MatchedTightIsoCJdRCanvasRachNew->GetPrimitive("MatchedTightIsoCJdR");
  TH1F* MatchedTightIsoRECOdRRachNew_ = (TH1F*)MatchedTightIsoRECOdRCanvasRachNew->GetPrimitive("MatchedTightIsoRECOdR");
  TH1F* MatchedDMFindCJdRRachNew_     = (TH1F*)MatchedDMFindCJdRCanvasRachNew->GetPrimitive("MatchedDMFindCJdR");
  TH1F* MatchedDMFindRECOdRRachNew_   = (TH1F*)MatchedDMFindRECOdRCanvasRachNew->GetPrimitive("MatchedDMFindRECOdR");
  TH1F* MatchedCJdRRachNew_   = (TH1F*)MatchedCJdRCanvasRachNew->GetPrimitive("MatchedCJdR");
  TH1F* MatchedRECOdRRachNew_ = (TH1F*)MatchedRECOdRCanvasRachNew->GetPrimitive("MatchedRECOdR");

  TH1F* MatchedLooseIsoCJPtGenRachNew_   = (TH1F*)MatchedLooseIsoCJPtGenCanvasRachNew->GetPrimitive("MatchedLooseIsoCJPtGen");
  TH1F* MatchedLooseIsoRECOPtGenRachNew_ = (TH1F*)MatchedLooseIsoRECOPtGenCanvasRachNew->GetPrimitive("MatchedLooseIsoRECOPtGen");
  TH1F* MatchedMedIsoCJPtGenRachNew_     = (TH1F*)MatchedMedIsoCJPtGenCanvasRachNew->GetPrimitive("MatchedMedIsoCJPtGen");
  TH1F* MatchedMedIsoRECOPtGenRachNew_   = (TH1F*)MatchedMedIsoRECOPtGenCanvasRachNew->GetPrimitive("MatchedMedIsoRECOPtGen");
  TH1F* MatchedTightIsoCJPtGenRachNew_   = (TH1F*)MatchedTightIsoCJPtGenCanvasRachNew->GetPrimitive("MatchedTightIsoCJPtGen");
  TH1F* MatchedTightIsoRECOPtGenRachNew_ = (TH1F*)MatchedTightIsoRECOPtGenCanvasRachNew->GetPrimitive("MatchedTightIsoRECOPtGen");
  TH1F* MatchedDMFindCJPtGenRachNew_     = (TH1F*)MatchedDMFindCJPtGenCanvasRachNew->GetPrimitive("MatchedDMFindCJPtGen");
  TH1F* MatchedDMFindRECOPtGenRachNew_   = (TH1F*)MatchedDMFindRECOPtGenCanvasRachNew->GetPrimitive("MatchedDMFindRECOPtGen");
  TH1F* MatchedCJPtGenRachNew_   = (TH1F*)MatchedCJPtGenCanvasRachNew->GetPrimitive("MatchedCJPtGen");
  TH1F* MatchedRECOPtGenRachNew_ = (TH1F*)MatchedRECOPtGenCanvasRachNew->GetPrimitive("MatchedRECOPtGen");
  
// Getting the histograms for the DM New
  TH1F* MatchedLooseIsoCJPtKyleNew_   = (TH1F*)MatchedLooseIsoCJPtCanvasKyleNew->GetPrimitive("MatchedLooseIsoCJPt");
  TH1F* MatchedLooseIsoRECOPtKyleNew_ = (TH1F*)MatchedLooseIsoRECOPtCanvasKyleNew->GetPrimitive("MatchedLooseIsoRECOPt");
  TH1F* MatchedMedIsoCJPtKyleNew_     = (TH1F*)MatchedMedIsoCJPtCanvasKyleNew->GetPrimitive("MatchedMedIsoCJPt");
  TH1F* MatchedMedIsoRECOPtKyleNew_   = (TH1F*)MatchedMedIsoRECOPtCanvasKyleNew->GetPrimitive("MatchedMedIsoRECOPt");
  TH1F* MatchedTightIsoCJPtKyleNew_   = (TH1F*)MatchedTightIsoCJPtCanvasKyleNew->GetPrimitive("MatchedTightIsoCJPt");
  TH1F* MatchedTightIsoRECOPtKyleNew_ = (TH1F*)MatchedTightIsoRECOPtCanvasKyleNew->GetPrimitive("MatchedTightIsoRECOPt");
  TH1F* MatchedDMFindCJPtKyleNew_     = (TH1F*)MatchedDMFindCJPtCanvasKyleNew->GetPrimitive("MatchedDMFindCJPt");
  TH1F* MatchedDMFindRECOPtKyleNew_   = (TH1F*)MatchedDMFindRECOPtCanvasKyleNew->GetPrimitive("MatchedDMFindRECOPt");
  TH1F* MatchedCJPtKyleNew_   = (TH1F*)MatchedCJPtCanvasKyleNew->GetPrimitive("MatchedCJPt");
  TH1F* MatchedRECOPtKyleNew_ = (TH1F*)MatchedRECOPtCanvasKyleNew->GetPrimitive("MatchedRECOPt");

  TH1F* MatchedLooseIsoCJdRKyleNew_   = (TH1F*)MatchedLooseIsoCJdRCanvasKyleNew->GetPrimitive("MatchedLooseIsoCJdR");
  TH1F* MatchedLooseIsoRECOdRKyleNew_ = (TH1F*)MatchedLooseIsoRECOdRCanvasKyleNew->GetPrimitive("MatchedLooseIsoRECOdR");
  TH1F* MatchedMedIsoCJdRKyleNew_     = (TH1F*)MatchedMedIsoCJdRCanvasKyleNew->GetPrimitive("MatchedMedIsoCJdR");
  TH1F* MatchedMedIsoRECOdRKyleNew_   = (TH1F*)MatchedMedIsoRECOdRCanvasKyleNew->GetPrimitive("MatchedMedIsoRECOdR");
  TH1F* MatchedTightIsoCJdRKyleNew_   = (TH1F*)MatchedTightIsoCJdRCanvasKyleNew->GetPrimitive("MatchedTightIsoCJdR");
  TH1F* MatchedTightIsoRECOdRKyleNew_ = (TH1F*)MatchedTightIsoRECOdRCanvasKyleNew->GetPrimitive("MatchedTightIsoRECOdR");
  TH1F* MatchedDMFindCJdRKyleNew_     = (TH1F*)MatchedDMFindCJdRCanvasKyleNew->GetPrimitive("MatchedDMFindCJdR");
  TH1F* MatchedDMFindRECOdRKyleNew_   = (TH1F*)MatchedDMFindRECOdRCanvasKyleNew->GetPrimitive("MatchedDMFindRECOdR");
  TH1F* MatchedCJdRKyleNew_   = (TH1F*)MatchedCJdRCanvasKyleNew->GetPrimitive("MatchedCJdR");
  TH1F* MatchedRECOdRKyleNew_ = (TH1F*)MatchedRECOdRCanvasKyleNew->GetPrimitive("MatchedRECOdR");

  TH1F* MatchedLooseIsoCJPtGenKyleNew_   = (TH1F*)MatchedLooseIsoCJPtGenCanvasKyleNew->GetPrimitive("MatchedLooseIsoCJPtGen");
  TH1F* MatchedLooseIsoRECOPtGenKyleNew_ = (TH1F*)MatchedLooseIsoRECOPtGenCanvasKyleNew->GetPrimitive("MatchedLooseIsoRECOPtGen");
  TH1F* MatchedMedIsoCJPtGenKyleNew_     = (TH1F*)MatchedMedIsoCJPtGenCanvasKyleNew->GetPrimitive("MatchedMedIsoCJPtGen");
  TH1F* MatchedMedIsoRECOPtGenKyleNew_   = (TH1F*)MatchedMedIsoRECOPtGenCanvasKyleNew->GetPrimitive("MatchedMedIsoRECOPtGen");
  TH1F* MatchedTightIsoCJPtGenKyleNew_   = (TH1F*)MatchedTightIsoCJPtGenCanvasKyleNew->GetPrimitive("MatchedTightIsoCJPtGen");
  TH1F* MatchedTightIsoRECOPtGenKyleNew_ = (TH1F*)MatchedTightIsoRECOPtGenCanvasKyleNew->GetPrimitive("MatchedTightIsoRECOPtGen");
  TH1F* MatchedDMFindCJPtGenKyleNew_     = (TH1F*)MatchedDMFindCJPtGenCanvasKyleNew->GetPrimitive("MatchedDMFindCJPtGen");
  TH1F* MatchedDMFindRECOPtGenKyleNew_   = (TH1F*)MatchedDMFindRECOPtGenCanvasKyleNew->GetPrimitive("MatchedDMFindRECOPtGen");
  TH1F* MatchedCJPtGenKyleNew_   = (TH1F*)MatchedCJPtGenCanvasKyleNew->GetPrimitive("MatchedCJPtGen");
  TH1F* MatchedRECOPtGenKyleNew_ = (TH1F*)MatchedRECOPtGenCanvasKyleNew->GetPrimitive("MatchedRECOPtGen");

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

  TCanvas LoosePtNewGenMatchCanvas("LoosePtNewGenMatchCanvas","",600,600);
  TCanvas MedPtNewGenMatchCanvas("MedPtNewGenMatchCanvas","",600,600);
  TCanvas TightPtNewGenMatchCanvas("TightPtNewGenMatchCanvas","",600,600);
  TCanvas DMPtNewGenMatchCanvas("DMPtNewGenMatchCanvas","",600,600);

  TCanvas LoosedRNewGenMatchCanvas("LoosedRNewGenMatchCanvas","",600,600);
  TCanvas MeddRNewGenMatchCanvas("MeddRNewGenMatchCanvas","",600,600);
  TCanvas TightdRNewGenMatchCanvas("TightdRNewGenMatchCanvas","",600,600);
  TCanvas DMdRNewGenMatchCanvas("DMdRNewGenMatchCanvas","",600,600);

  TCanvas LoosePtGenNewGenMatchCanvas("LoosePtGenNewGenMatchCanvas","",600,600);
  TCanvas MedPtGenNewGenMatchCanvas("MedPtGenNewGenMatchCanvas","",600,600);
  TCanvas TightPtGenNewGenMatchCanvas("TightPtGenNewGenMatchCanvas","",600,600);
  TCanvas DMPtGenNewGenMatchCanvas("DMPtGenNewGenMatchCanvas","",600,600);

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

  LoosePtNewGenMatchCanvas.SetGrid(1,1);
  MedPtNewGenMatchCanvas.SetGrid(1,1);
  TightPtNewGenMatchCanvas.SetGrid(1,1);
  DMPtNewGenMatchCanvas.SetGrid(1,1);
 
  LoosedRNewGenMatchCanvas.SetGrid(1,1);
  MeddRNewGenMatchCanvas.SetGrid(1,1);
  TightdRNewGenMatchCanvas.SetGrid(1,1);
  DMdRNewGenMatchCanvas.SetGrid(1,1);
 
  LoosePtGenNewGenMatchCanvas.SetGrid(1,1);
  MedPtGenNewGenMatchCanvas.SetGrid(1,1);
  TightPtGenNewGenMatchCanvas.SetGrid(1,1);
  DMPtGenNewGenMatchCanvas.SetGrid(1,1);


cout << "Canvases created" << endl;

  // To  divide histograms for efficiency with new Gen Matching New DM's
  TGraphAsymmErrors* FinalEffLooseIsoRECOPtRachNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoRECOPtRachNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoRECOPtRachNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindRECOPtRachNew_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoCJPtRachNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoCJPtRachNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoCJPtRachNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindCJPtRachNew_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoRECOdRRachNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoRECOdRRachNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoRECOdRRachNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindRECOdRRachNew_ = new TGraphAsymmErrors(30);
  
  TGraphAsymmErrors* FinalEffLooseIsoCJdRRachNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoCJdRRachNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoCJdRRachNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindCJdRRachNew_ = new TGraphAsymmErrors(30);
  
  TGraphAsymmErrors* FinalEffLooseIsoRECOPtGenRachNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoRECOPtGenRachNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoRECOPtGenRachNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindRECOPtGenRachNew_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoCJPtGenRachNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoCJPtGenRachNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoCJPtGenRachNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindCJPtGenRachNew_ = new TGraphAsymmErrors(30);

  // To divide histograms for efficiency with Old Gen Matching New DM's
  TGraphAsymmErrors* FinalEffLooseIsoRECOPtKyleNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoRECOPtKyleNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoRECOPtKyleNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindRECOPtKyleNew_ = new TGraphAsymmErrors(30);
  
  TGraphAsymmErrors* FinalEffLooseIsoCJPtKyleNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoCJPtKyleNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoCJPtKyleNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindCJPtKyleNew_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoRECOdRKyleNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoRECOdRKyleNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoRECOdRKyleNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindRECOdRKyleNew_ = new TGraphAsymmErrors(30);
  
  TGraphAsymmErrors* FinalEffLooseIsoCJdRKyleNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoCJdRKyleNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoCJdRKyleNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindCJdRKyleNew_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoRECOPtGenKyleNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoRECOPtGenKyleNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoRECOPtGenKyleNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindRECOPtGenKyleNew_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoCJPtGenKyleNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoCJPtGenKyleNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoCJPtGenKyleNew_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindCJPtGenKyleNew_ = new TGraphAsymmErrors(30);

  // To divide histograms for efficiency with new Gen Matching New DM's
  FinalEffLooseIsoRECOPtRachNew_->Divide(MatchedLooseIsoRECOPtRachNew_, MatchedRECOPtRachNew_);
  FinalEffMedIsoRECOPtRachNew_->Divide(MatchedMedIsoRECOPtRachNew_,     MatchedRECOPtRachNew_);
  FinalEffTightIsoRECOPtRachNew_->Divide(MatchedTightIsoRECOPtRachNew_, MatchedRECOPtRachNew_);
  FinalEffDMFindRECOPtRachNew_->Divide(MatchedDMFindRECOPtRachNew_,     MatchedRECOPtRachNew_);
  FinalEffLooseIsoCJPtRachNew_->Divide(MatchedLooseIsoCJPtRachNew_, MatchedCJPtRachNew_);
  FinalEffMedIsoCJPtRachNew_->Divide(MatchedMedIsoCJPtRachNew_,     MatchedCJPtRachNew_);
  FinalEffTightIsoCJPtRachNew_->Divide(MatchedTightIsoCJPtRachNew_, MatchedCJPtRachNew_);
  FinalEffDMFindCJPtRachNew_->Divide(MatchedDMFindCJPtRachNew_,     MatchedCJPtRachNew_);

  FinalEffLooseIsoRECOdRRachNew_->Divide(MatchedLooseIsoRECOdRRachNew_, MatchedRECOdRRachNew_);
  FinalEffMedIsoRECOdRRachNew_->Divide(MatchedMedIsoRECOdRRachNew_,     MatchedRECOdRRachNew_);
  FinalEffTightIsoRECOdRRachNew_->Divide(MatchedTightIsoRECOdRRachNew_, MatchedRECOdRRachNew_);
  FinalEffDMFindRECOdRRachNew_->Divide(MatchedDMFindRECOdRRachNew_,     MatchedRECOdRRachNew_);
  FinalEffLooseIsoCJdRRachNew_->Divide(MatchedLooseIsoCJdRRachNew_, MatchedCJdRRachNew_);
  FinalEffMedIsoCJdRRachNew_->Divide(MatchedMedIsoCJdRRachNew_,     MatchedCJdRRachNew_);
  FinalEffTightIsoCJdRRachNew_->Divide(MatchedTightIsoCJdRRachNew_, MatchedCJdRRachNew_);
  FinalEffDMFindCJdRRachNew_->Divide(MatchedDMFindCJdRRachNew_,     MatchedCJdRRachNew_);

  FinalEffLooseIsoRECOPtGenRachNew_->Divide(MatchedLooseIsoRECOPtGenRachNew_, MatchedRECOPtGenRachNew_);
  FinalEffMedIsoRECOPtGenRachNew_->Divide(MatchedMedIsoRECOPtGenRachNew_,     MatchedRECOPtGenRachNew_);
  FinalEffTightIsoRECOPtGenRachNew_->Divide(MatchedTightIsoRECOPtGenRachNew_, MatchedRECOPtGenRachNew_);
  FinalEffDMFindRECOPtGenRachNew_->Divide(MatchedDMFindRECOPtGenRachNew_,     MatchedRECOPtGenRachNew_);
  FinalEffLooseIsoCJPtGenRachNew_->Divide(MatchedLooseIsoCJPtGenRachNew_, MatchedCJPtGenRachNew_);
  FinalEffMedIsoCJPtGenRachNew_->Divide(MatchedMedIsoCJPtGenRachNew_,     MatchedCJPtGenRachNew_);
  FinalEffTightIsoCJPtGenRachNew_->Divide(MatchedTightIsoCJPtGenRachNew_, MatchedCJPtGenRachNew_);
  FinalEffDMFindCJPtGenRachNew_->Divide(MatchedDMFindCJPtGenRachNew_,     MatchedCJPtGenRachNew_);
  
  // To divide histograms for efficiency with OLD Gen Matching New DM's 
  FinalEffLooseIsoRECOPtKyleNew_->Divide(MatchedLooseIsoRECOPtKyleNew_, MatchedRECOPtKyleNew_);
  FinalEffMedIsoRECOPtKyleNew_->Divide(MatchedMedIsoRECOPtKyleNew_,     MatchedRECOPtKyleNew_);
  FinalEffTightIsoRECOPtKyleNew_->Divide(MatchedTightIsoRECOPtKyleNew_, MatchedRECOPtKyleNew_);
  FinalEffDMFindRECOPtKyleNew_->Divide(MatchedDMFindRECOPtKyleNew_,     MatchedRECOPtKyleNew_);
  FinalEffLooseIsoCJPtKyleNew_->Divide(MatchedLooseIsoCJPtKyleNew_, MatchedCJPtKyleNew_);
  FinalEffMedIsoCJPtKyleNew_->Divide(MatchedMedIsoCJPtKyleNew_,     MatchedCJPtKyleNew_);
  FinalEffTightIsoCJPtKyleNew_->Divide(MatchedTightIsoCJPtKyleNew_, MatchedCJPtKyleNew_);
  FinalEffDMFindCJPtKyleNew_->Divide(MatchedDMFindCJPtKyleNew_,     MatchedCJPtKyleNew_);

  FinalEffLooseIsoRECOdRKyleNew_->Divide(MatchedLooseIsoRECOdRKyleNew_, MatchedRECOdRKyleNew_);
  FinalEffMedIsoRECOdRKyleNew_->Divide(MatchedMedIsoRECOdRKyleNew_,     MatchedRECOdRKyleNew_);
  FinalEffTightIsoRECOdRKyleNew_->Divide(MatchedTightIsoRECOdRKyleNew_, MatchedRECOdRKyleNew_);
  FinalEffDMFindRECOdRKyleNew_->Divide(MatchedDMFindRECOdRKyleNew_,     MatchedRECOdRKyleNew_);
  FinalEffLooseIsoCJdRKyleNew_->Divide(MatchedLooseIsoCJdRKyleNew_, MatchedCJdRKyleNew_);
  FinalEffMedIsoCJdRKyleNew_->Divide(MatchedMedIsoCJdRKyleNew_,     MatchedCJdRKyleNew_);
  FinalEffTightIsoCJdRKyleNew_->Divide(MatchedTightIsoCJdRKyleNew_, MatchedCJdRKyleNew_);
  FinalEffDMFindCJdRKyleNew_->Divide(MatchedDMFindCJdRKyleNew_,     MatchedCJdRKyleNew_);

  FinalEffLooseIsoRECOPtGenKyleNew_->Divide(MatchedLooseIsoRECOPtGenKyleNew_, MatchedRECOPtGenKyleNew_);
  FinalEffMedIsoRECOPtGenKyleNew_->Divide(MatchedMedIsoRECOPtGenKyleNew_,     MatchedRECOPtGenKyleNew_);
  FinalEffTightIsoRECOPtGenKyleNew_->Divide(MatchedTightIsoRECOPtGenKyleNew_, MatchedRECOPtGenKyleNew_);
  FinalEffDMFindRECOPtGenKyleNew_->Divide(MatchedDMFindRECOPtGenKyleNew_,     MatchedRECOPtGenKyleNew_);
  FinalEffLooseIsoCJPtGenKyleNew_->Divide(MatchedLooseIsoCJPtGenKyleNew_, MatchedCJPtGenKyleNew_);
  FinalEffMedIsoCJPtGenKyleNew_->Divide(MatchedMedIsoCJPtGenKyleNew_,     MatchedCJPtGenKyleNew_);
  FinalEffTightIsoCJPtGenKyleNew_->Divide(MatchedTightIsoCJPtGenKyleNew_, MatchedCJPtGenKyleNew_);
  FinalEffDMFindCJPtGenKyleNew_->Divide(MatchedDMFindCJPtGenKyleNew_,     MatchedCJPtGenKyleNew_);

  //Set Colors for New Gen Matching New DM's
  setGraphOptions(FinalEffLooseIsoCJPtRachNew_, kBlue+1, 0.7, 10, "p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffLooseIsoRECOPtRachNew_, kGreen+1, 0.7, 10, "p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoCJPtRachNew_, kBlue+1, 0.7, 10, "p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoRECOPtRachNew_, kGreen+1, 0.7, 10, "p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoCJPtRachNew_, kBlue+1, 0.7, 10, "p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoRECOPtRachNew_, kGreen+1, 0.7, 10, "p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindCJPtRachNew_, kBlue+1, 0.7, 10, "p_{T}", "#epsilon (DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindRECOPtRachNew_, kGreen+1, 0.7, 10, "p_{T}", "#epsilon (DMFinding + GM / GM)");

  setGraphOptions(FinalEffLooseIsoCJdRRachNew_, kBlue+1, 0.7, 10, "#DeltaR", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffLooseIsoRECOdRRachNew_, kGreen+1, 0.7, 10, "#DeltaR", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoCJdRRachNew_, kBlue+1, 0.7, 10, "#DeltaR", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoRECOdRRachNew_, kGreen+1, 0.7, 10, "#DeltaR", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoCJdRRachNew_, kBlue+1, 0.7, 10, "#DeltaR", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoRECOdRRachNew_, kGreen+1, 0.7, 10, "#DeltaR", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindCJdRRachNew_, kBlue+1, 0.7, 10, "#DeltaR", "#epsilon (DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindRECOdRRachNew_, kGreen+1, 0.7, 10, "#DeltaR", "#epsilon (DMFinding + GM / GM)");
  
  setGraphOptions(FinalEffLooseIsoCJPtGenRachNew_, kBlue+1, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffLooseIsoRECOPtGenRachNew_, kGreen+1, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoCJPtGenRachNew_, kBlue+1, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoRECOPtGenRachNew_, kGreen+1, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoCJPtGenRachNew_, kBlue+1, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoRECOPtGenRachNew_, kGreen+1, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindCJPtGenRachNew_, kBlue+1, 0.7, 10, "Gen Visible p_{T}", "#epsilon (DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindRECOPtGenRachNew_, kGreen+1, 0.7, 10, "Gen Visible p_{T}", "#epsilon (DMFinding + GM / GM)");

  //Set Colors for New Gen Matching New DM's
  setGraphOptions(FinalEffLooseIsoCJPtKyleNew_, kViolet-5, 0.7, 10, "p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffLooseIsoRECOPtKyleNew_, kGray+1, 0.7, 10, "p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoCJPtKyleNew_, kViolet-5, 0.7, 10, "p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoRECOPtKyleNew_, kGray+1, 0.7, 10, "p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoCJPtKyleNew_, kViolet-5, 0.7, 10, "p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoRECOPtKyleNew_, kGray+1, 0.7, 10, "p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindCJPtKyleNew_, kViolet-5, 0.7, 10, "p_{T}", "#epsilon (DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindRECOPtKyleNew_, kGray+1, 0.7, 10, "p_{T}", "#epsilon (DMFinding + GM / GM)");
  
  setGraphOptions(FinalEffLooseIsoCJdRKyleNew_, kViolet-5, 0.7, 10, "#DeltaR", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffLooseIsoRECOdRKyleNew_, kGray+1, 0.7, 10, "#DeltaR", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoCJdRKyleNew_, kViolet-5, 0.7, 10, "#DeltaR", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoRECOdRKyleNew_, kGray+1, 0.7, 10, "#DeltaR", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoCJdRKyleNew_, kViolet-5, 0.7, 10, "#DeltaR", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoRECOdRKyleNew_, kGray+1, 0.7, 10, "#DeltaR", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindCJdRKyleNew_, kViolet-5, 0.7, 10, "#DeltaR", "#epsilon (DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindRECOdRKyleNew_, kGray+1, 0.7, 10, "#DeltaR", "#epsilon (DMFinding + GM / GM)");

  setGraphOptions(FinalEffLooseIsoCJPtGenKyleNew_, kViolet-5, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffLooseIsoRECOPtGenKyleNew_, kGray+1, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoCJPtGenKyleNew_, kViolet-5, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoRECOPtGenKyleNew_, kGray+1, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoCJPtGenKyleNew_, kViolet-5, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoRECOPtGenKyleNew_, kGray+1, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindCJPtGenKyleNew_, kViolet-5, 0.7, 10, "Gen Visible p_{T}", "#epsilon (DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindRECOPtGenKyleNew_, kGray+1, 0.7, 10, "Gen Visible p_{T}", "#epsilon (DMFinding + GM / GM)");

cout << "Attributes set." << endl;  

  leg = new TLegend(0.1,0.7,0.25,0.9);
  leg->AddEntry(FinalEffDMFindRECOPtGenRachNew_, "No  CJ,New Gen,New DM","L");
  leg->AddEntry(FinalEffDMFindRECOPtGenKyleNew_, "No  CJ,Old Gen,New DM","L");
  leg->AddEntry(FinalEffDMFindCJPtGenRachNew_,   "Yes CJ,New Gen,New DM","L");
  leg->AddEntry(FinalEffDMFindCJPtGenKyleNew_,   "Yes CJ,Old Gen,New DM","L");

  legMatch = new TLegend(0.1,0.7,0.25,0.9);
  legMatch->AddEntry(FinalEffDMFindRECOPtGenRachNew_, "No  CJ,New Gen,New DM","L");
  legMatch->AddEntry(FinalEffDMFindCJPtGenRachNew_,   "Yes CJ,New Gen,New DM","L");



  LoosePtCanvas.cd();
  FinalEffLooseIsoCJPtRachNew_->Draw();
  FinalEffLooseIsoRECOPtRachNew_->Draw("SAME");
  FinalEffLooseIsoCJPtKyleNew_->Draw("SAME");
  FinalEffLooseIsoRECOPtKyleNew_->Draw("SAME");
  leg->Draw();

  LoosePtNewGenMatchCanvas.cd();
  FinalEffLooseIsoCJPtRachNew_->Draw();
  FinalEffLooseIsoRECOPtRachNew_->Draw("SAME");
  legMatch->Draw();

  MedPtCanvas.cd();
  FinalEffMedIsoCJPtRachNew_->Draw();
  FinalEffMedIsoRECOPtRachNew_->Draw("SAME");
  FinalEffMedIsoCJPtKyleNew_->Draw("SAME");
  FinalEffMedIsoRECOPtKyleNew_->Draw("SAME");
  leg->Draw();

  MedPtNewGenMatchCanvas.cd();
  FinalEffMedIsoCJPtRachNew_->Draw();
  FinalEffMedIsoRECOPtRachNew_->Draw("SAME");
  legMatch->Draw();
  

  TightPtCanvas.cd();
  FinalEffTightIsoCJPtRachNew_->Draw();
  FinalEffTightIsoRECOPtRachNew_->Draw("SAME");
  FinalEffTightIsoCJPtKyleNew_->Draw("SAME");
  FinalEffTightIsoRECOPtKyleNew_->Draw("SAME");
  leg->Draw();

  TightPtNewGenMatchCanvas.cd();
  FinalEffTightIsoCJPtRachNew_->Draw();
  FinalEffTightIsoRECOPtRachNew_->Draw("SAME");
  legMatch->Draw();
  

  DMPtCanvas.cd();
  FinalEffDMFindCJPtRachNew_->Draw();
  FinalEffDMFindRECOPtRachNew_->Draw("SAME");
  FinalEffDMFindCJPtKyleNew_->Draw("SAME");
  FinalEffDMFindRECOPtKyleNew_->Draw("SAME");
  leg->Draw();

  DMPtNewGenMatchCanvas.cd();
  FinalEffDMFindCJPtRachNew_->Draw();
  FinalEffDMFindRECOPtRachNew_->Draw("SAME");
  legMatch->Draw();
  


  LoosedRCanvas.cd();
  FinalEffLooseIsoCJdRRachNew_->Draw();
  FinalEffLooseIsoRECOdRRachNew_->Draw("SAME");
  FinalEffLooseIsoCJdRKyleNew_->Draw("SAME");
  FinalEffLooseIsoRECOdRKyleNew_->Draw("SAME");
  leg->Draw();

  LoosedRNewGenMatchCanvas.cd();
  FinalEffLooseIsoCJdRRachNew_->Draw();
  FinalEffLooseIsoRECOdRRachNew_->Draw("SAME");
  legMatch->Draw();
  
 
  MeddRCanvas.cd();
  FinalEffMedIsoCJdRRachNew_->Draw();
  FinalEffMedIsoRECOdRRachNew_->Draw("SAME");
  FinalEffMedIsoCJdRKyleNew_->Draw("SAME");
  FinalEffMedIsoRECOdRKyleNew_->Draw("SAME");
  leg->Draw();
 
  MeddRNewGenMatchCanvas.cd();
  FinalEffMedIsoCJdRRachNew_->Draw();
  FinalEffMedIsoRECOdRRachNew_->Draw("SAME");
  legMatch->Draw();
 
  TightdRCanvas.cd();
  FinalEffTightIsoCJdRRachNew_->Draw();
  FinalEffTightIsoRECOdRRachNew_->Draw("SAME");
  FinalEffTightIsoCJdRKyleNew_->Draw("SAME");
  FinalEffTightIsoRECOdRKyleNew_->Draw("SAME");
  leg->Draw();

  TightdRNewGenMatchCanvas.cd();
  FinalEffTightIsoCJdRRachNew_->Draw();
  FinalEffTightIsoRECOdRRachNew_->Draw("SAME");
  legMatch->Draw();

  DMdRCanvas.cd();
  FinalEffDMFindCJdRRachNew_->Draw();
  FinalEffDMFindRECOdRRachNew_->Draw("SAME");
  FinalEffDMFindCJdRKyleNew_->Draw("SAME");
  FinalEffDMFindRECOdRKyleNew_->Draw("SAME");
  leg->Draw();

  DMdRNewGenMatchCanvas.cd();
  FinalEffDMFindCJdRRachNew_->Draw();
  FinalEffDMFindRECOdRRachNew_->Draw("SAME");
  legMatch->Draw();
 
  LoosePtGenCanvas.cd();
  FinalEffLooseIsoCJPtGenRachNew_->Draw();
  FinalEffLooseIsoRECOPtGenRachNew_->Draw("SAME");
  FinalEffLooseIsoCJPtGenKyleNew_->Draw("SAME");
  FinalEffLooseIsoRECOPtGenKyleNew_->Draw("SAME");
  leg->Draw();

  LoosePtGenNewGenMatchCanvas.cd();
  FinalEffLooseIsoCJPtGenRachNew_->Draw();
  FinalEffLooseIsoRECOPtGenRachNew_->Draw("SAME");
  legMatch->Draw();
 
  MedPtGenCanvas.cd();
  FinalEffMedIsoCJPtGenRachNew_->Draw();
  FinalEffMedIsoRECOPtGenRachNew_->Draw("SAME");
  FinalEffMedIsoCJPtGenKyleNew_->Draw("SAME");
  FinalEffMedIsoRECOPtGenKyleNew_->Draw("SAME");
  leg->Draw();

  MedPtGenNewGenMatchCanvas.cd();
  FinalEffMedIsoCJPtGenRachNew_->Draw();
  FinalEffMedIsoRECOPtGenRachNew_->Draw("SAME");
  legMatch->Draw();
 
  TightPtGenCanvas.cd();
  FinalEffTightIsoCJPtGenRachNew_->Draw("SAME");
  FinalEffTightIsoRECOPtGenRachNew_->Draw("SAME");
  FinalEffTightIsoCJPtGenKyleNew_->Draw("SAME");
  FinalEffTightIsoRECOPtGenKyleNew_->Draw("SAME");
  leg->Draw();

  TightPtGenNewGenMatchCanvas.cd();
  FinalEffTightIsoCJPtGenRachNew_->Draw();
  FinalEffTightIsoRECOPtGenRachNew_->Draw("SAME");
  legMatch->Draw();

  DMPtGenCanvas.cd();
  FinalEffDMFindCJPtGenRachNew_->Draw();
  FinalEffDMFindRECOPtGenRachNew_->Draw("SAME");
  FinalEffDMFindCJPtGenKyleNew_->Draw("SAME");
  FinalEffDMFindRECOPtGenKyleNew_->Draw("SAME");
  leg->Draw();

  DMPtGenNewGenMatchCanvas.cd();
  FinalEffDMFindCJPtGenRachNew_->Draw();
  FinalEffDMFindRECOPtGenRachNew_->Draw("SAME");
  legMatch->Draw();

 
cout << "Histograms Drawn" << endl;

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

  LoosePtNewGenMatchCanvas.Write();
  MedPtNewGenMatchCanvas.Write();
  TightPtNewGenMatchCanvas.Write();
  DMPtNewGenMatchCanvas.Write();

  LoosedRNewGenMatchCanvas.Write();
  MeddRNewGenMatchCanvas.Write();
  TightdRNewGenMatchCanvas.Write();
  DMdRNewGenMatchCanvas.Write();

  LoosePtGenNewGenMatchCanvas.Write();
  MedPtGenNewGenMatchCanvas.Write();
  TightPtGenNewGenMatchCanvas.Write();
  DMPtGenNewGenMatchCanvas.Write();

  MatchedCJPtCanvasRachNew->Write();
  MatchedRECOPtCanvasRachNew->Write();
  MatchedCJdRCanvasRachNew->Write();
  MatchedRECOdRCanvasRachNew->Write();
  MatchedCJPtGenCanvasRachNew->Write();
  MatchedRECOPtGenCanvasRachNew->Write();

  MatchedCJPtCanvasKyleNew->Write();
  MatchedRECOPtCanvasKyleNew->Write();
  MatchedCJdRCanvasKyleNew->Write();
  MatchedRECOdRCanvasKyleNew->Write();
  MatchedCJPtGenCanvasKyleNew->Write();
  MatchedRECOPtGenCanvasKyleNew->Write();

  outFile->Write();
  outFile->Close();

cout << "end" << endl;

}//rootMacro_BBA_combine


void setGraphOptions(TGraphAsymmErrors* graph, const Color_t color, const Size_t size, const Style_t style, const char* xAxisTitle, const char* yAxisTitle)
{
  graph->SetMarkerColor(color);
  graph->SetMarkerSize(size);
  graph->SetMarkerStyle(style);
  graph->SetLineColor(color);
  graph->SetTitle("");
  graph->GetXaxis()->SetLabelFont(42);
  graph->GetXaxis()->SetLabelOffset(0.007);
  graph->GetXaxis()->SetLabelSize(.05);
  graph->GetXaxis()->SetTitleFont(42);
  graph->GetXaxis()->SetTitleSize(.05);
  graph->GetXaxis()->SetTitleOffset(.9);
  graph->GetXaxis()->SetTitle(xAxisTitle);

  graph->GetYaxis()->SetLabelFont(42);
  graph->GetYaxis()->SetLabelOffset(0.007);
  graph->GetYaxis()->SetLabelSize(.05);
  graph->GetYaxis()->SetTitleFont(42);
  graph->GetYaxis()->SetTitleSize(.05);
  graph->GetYaxis()->SetTitleOffset(.9);
  graph->GetYaxis()->SetTitle(yAxisTitle);
  graph->GetYaxis()->SetRangeUser(0.0, 1.1);

  for (unsigned int i = 0; i < graph->GetN(); i++)
  {
    graph->SetPointEXhigh(i, 0);
    graph->SetPointEXlow(i, 0);
  }//for

}


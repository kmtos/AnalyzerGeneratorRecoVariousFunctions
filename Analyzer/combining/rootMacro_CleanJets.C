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


void rootMacro_CleanJets()
{

  gStyle->SetOptStat(kFALSE);
  gStyle->SetEndErrorSize(7);
  TFile infileCJNew("/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/ggH125a9_GenTauDecayID_IndivCJ_NewDMFind_MVA_APR1/ggH125a9_GenTauDecayID_IndivCJ_NewDMFind_MVA_APR1_Plots.root");
  TFile infileCJOld("/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/ggH125a9_GenTauDecayID_IndivCJ_OldDMFind_MVA_APR1/ggH125a9_GenTauDecayID_IndivCJ_OldDMFind_MVA_APR1_Plots.root");
  TFile infileRECOOld("/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/ggH125a9_GenTauDecayID_IndivRECO_OldMFind_MVA_APR1/ggH125a9_GenTauDecayID_IndivRECO_OldMFind_MVA_APR1_Plots.root");
  TFile infileRECONew("/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/ggH125a9_GenTauDecayID_IndivRECO_NewMFind_MVA_APR1/ggH125a9_GenTauDecayID_IndivRECO_NewMFind_MVA_APR1_Plots.root");
  TFile infileKyleNew("/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/ggH125a9_OLDDecayID_NewDMFind_MVA_APR1/ggH125a9_OLDDecayID_NewDMFind_MVA_APR1_Plots.root");
  TFile infileKyleOld("/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/ggH125a9_OLDDecayID_OldDMFind_MVA_APR1/ggH125a9_OLDDecayID_OldDMFind_MVA_APR1_Plots.root");

cout << "Getting Files" << std::endl;
//  TFile infileCJNew("/afs/cern.ch/user/k/ktos/TEMP_MOUNT_RECOVERY/ggH125a9_GenTauDecayID_IndivCJ_NewDMFind_MVA_MAR10/ggH125a9_GenTauDecayID_IndivCJ_NewDMFind_MVA_MAR10_Plots.root");
//  TFile infileCJOld("/afs/cern.ch/user/k/ktos/TEMP_MOUNT_RECOVERY/ggH125a9_GenTauDecayID_IndivCJ_OldDMFind_MVA_MAR10/ggH125a9_GenTauDecayID_IndivCJ_OldDMFind_MVA_MAR10_Plots.root");
//  TFile infileRECOOld("/afs/cern.ch/user/k/ktos/TEMP_MOUNT_RECOVERY/ggH125a9_GenTauDecayID_IndivRECO_OldDMFind_MVA_MAR10/ggH125a9_GenTauDecayID_IndivRECO_OldDMFind_MVA_MAR10_Plots.root");
//  TFile infileRECONew("/afs/cern.ch/user/k/ktos/TEMP_MOUNT_RECOVERY/ggH125a9_GenTauDecayID_IndivRECO_NewDMFind_MVA_MAR10/ggH125a9_GenTauDecayID_IndivRECO_NewDMFind_MVA_MAR10_Plots.root");
//  TFile infileKyleNew("/afs/cern.ch/user/k/ktos/TEMP_MOUNT_RECOVERY/ggH125a9_OLDDecayID_NewDMFind_MVA_MAR10/ggH125a9_OLDDecayID_NewDMFind_MVA_MAR10_Plots.root");
//  TFile infileKyleOld("/afs/cern.ch/user/k/ktos/TEMP_MOUNT_RECOVERY/ggH125a9_OLDDecayID_OldDMFind_MVA_MAR10/ggH125a9_OLDDecayID_OldDMFind_MVA_MAR10_Plots.root");


  TFile *outFile = new TFile("combHist_CleanJets_h125a9_APR1.root", "RECREATE");

cout << "Files Created" << endl;

  // CJ/RECO Comparison histograms Old DM
  TCanvas *MatchedLooseIsoCJPtCanvasRachOld   = (TCanvas*)infileCJOld.Get("MatchedLooseIsoCJPt");
  TCanvas *MatchedLooseIsoRECOPtCanvasRachOld = (TCanvas*)infileRECOOld.Get("MatchedLooseIsoRECOPt");
  TCanvas *MatchedMedIsoCJPtCanvasRachOld     = (TCanvas*)infileCJOld.Get("MatchedMedIsoCJPt");
  TCanvas *MatchedMedIsoRECOPtCanvasRachOld   = (TCanvas*)infileRECOOld.Get("MatchedMedIsoRECOPt");
  TCanvas *MatchedTightIsoCJPtCanvasRachOld   = (TCanvas*)infileCJOld.Get("MatchedTightIsoCJPt");
  TCanvas *MatchedTightIsoRECOPtCanvasRachOld = (TCanvas*)infileRECOOld.Get("MatchedTightIsoRECOPt");
  TCanvas *MatchedDMFindCJPtCanvasRachOld     = (TCanvas*)infileCJOld.Get("MatchedDMFindCJPt");
  TCanvas *MatchedDMFindRECOPtCanvasRachOld   = (TCanvas*)infileRECOOld.Get("MatchedDMFindRECOPt");
  TCanvas *MatchedCJPtCanvasRachOld   = (TCanvas*)infileCJOld.Get("MatchedCJPt");
  TCanvas *MatchedRECOPtCanvasRachOld = (TCanvas*)infileRECOOld.Get("MatchedRECOPt");

  TCanvas *MatchedLooseIsoCJdRCanvasRachOld   = (TCanvas*)infileCJOld.Get("MatchedLooseIsoCJdR");
  TCanvas *MatchedLooseIsoRECOdRCanvasRachOld = (TCanvas*)infileRECOOld.Get("MatchedLooseIsoRECOdR");
  TCanvas *MatchedMedIsoCJdRCanvasRachOld     = (TCanvas*)infileCJOld.Get("MatchedMedIsoCJdR");
  TCanvas *MatchedMedIsoRECOdRCanvasRachOld   = (TCanvas*)infileRECOOld.Get("MatchedMedIsoRECOdR");
  TCanvas *MatchedTightIsoCJdRCanvasRachOld   = (TCanvas*)infileCJOld.Get("MatchedTightIsoCJdR");
  TCanvas *MatchedTightIsoRECOdRCanvasRachOld = (TCanvas*)infileRECOOld.Get("MatchedTightIsoRECOdR");
  TCanvas *MatchedDMFindCJdRCanvasRachOld     = (TCanvas*)infileCJOld.Get("MatchedDMFindCJdR");
  TCanvas *MatchedDMFindRECOdRCanvasRachOld   = (TCanvas*)infileRECOOld.Get("MatchedDMFindRECOdR");
  TCanvas *MatchedCJdRCanvasRachOld   = (TCanvas*)infileCJOld.Get("MatchedCJdR");
  TCanvas *MatchedRECOdRCanvasRachOld = (TCanvas*)infileRECOOld.Get("MatchedRECOdR");

  TCanvas *MatchedLooseIsoCJPtGenCanvasRachOld   = (TCanvas*)infileCJOld.Get("MatchedLooseIsoCJPtGen");
  TCanvas *MatchedLooseIsoRECOPtGenCanvasRachOld = (TCanvas*)infileRECOOld.Get("MatchedLooseIsoRECOPtGen");
  TCanvas *MatchedMedIsoCJPtGenCanvasRachOld     = (TCanvas*)infileCJOld.Get("MatchedMedIsoCJPtGen");
  TCanvas *MatchedMedIsoRECOPtGenCanvasRachOld   = (TCanvas*)infileRECOOld.Get("MatchedMedIsoRECOPtGen");
  TCanvas *MatchedTightIsoCJPtGenCanvasRachOld   = (TCanvas*)infileCJOld.Get("MatchedTightIsoCJPtGen");
  TCanvas *MatchedTightIsoRECOPtGenCanvasRachOld = (TCanvas*)infileRECOOld.Get("MatchedTightIsoRECOPtGen");
  TCanvas *MatchedDMFindCJPtGenCanvasRachOld     = (TCanvas*)infileCJOld.Get("MatchedDMFindCJPtGen");
  TCanvas *MatchedDMFindRECOPtGenCanvasRachOld   = (TCanvas*)infileRECOOld.Get("MatchedDMFindRECOPtGen");
  TCanvas *MatchedCJPtGenCanvasRachOld   = (TCanvas*)infileCJOld.Get("MatchedCJPtGen");
  TCanvas *MatchedRECOPtGenCanvasRachOld = (TCanvas*)infileRECOOld.Get("MatchedRECOPtGen");

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

  // Kyle Gen Matching Old DM 
  TCanvas *MatchedLooseIsoCJPtCanvasKyleOld   = (TCanvas*)infileKyleOld.Get("MatchedLooseIsoCJPt");
  TCanvas *MatchedLooseIsoRECOPtCanvasKyleOld = (TCanvas*)infileKyleOld.Get("MatchedLooseIsoRECOPt");
  TCanvas *MatchedMedIsoCJPtCanvasKyleOld     = (TCanvas*)infileKyleOld.Get("MatchedMedIsoCJPt");
  TCanvas *MatchedMedIsoRECOPtCanvasKyleOld   = (TCanvas*)infileKyleOld.Get("MatchedMedIsoRECOPt");
  TCanvas *MatchedTightIsoCJPtCanvasKyleOld   = (TCanvas*)infileKyleOld.Get("MatchedTightIsoCJPt");
  TCanvas *MatchedTightIsoRECOPtCanvasKyleOld = (TCanvas*)infileKyleOld.Get("MatchedTightIsoRECOPt");
  TCanvas *MatchedDMFindCJPtCanvasKyleOld     = (TCanvas*)infileKyleOld.Get("MatchedDMFindCJPt");
  TCanvas *MatchedDMFindRECOPtCanvasKyleOld   = (TCanvas*)infileKyleOld.Get("MatchedDMFindRECOPt");
  TCanvas *MatchedCJPtCanvasKyleOld   = (TCanvas*)infileKyleOld.Get("MatchedCJPt");
  TCanvas *MatchedRECOPtCanvasKyleOld = (TCanvas*)infileKyleOld.Get("MatchedRECOPt");
  
  TCanvas *MatchedLooseIsoCJdRCanvasKyleOld   = (TCanvas*)infileKyleOld.Get("MatchedLooseIsoCJdR");
  TCanvas *MatchedLooseIsoRECOdRCanvasKyleOld = (TCanvas*)infileKyleOld.Get("MatchedLooseIsoRECOdR");
  TCanvas *MatchedMedIsoCJdRCanvasKyleOld     = (TCanvas*)infileKyleOld.Get("MatchedMedIsoCJdR");
  TCanvas *MatchedMedIsoRECOdRCanvasKyleOld   = (TCanvas*)infileKyleOld.Get("MatchedMedIsoRECOdR");
  TCanvas *MatchedTightIsoCJdRCanvasKyleOld   = (TCanvas*)infileKyleOld.Get("MatchedTightIsoCJdR");
  TCanvas *MatchedTightIsoRECOdRCanvasKyleOld = (TCanvas*)infileKyleOld.Get("MatchedTightIsoRECOdR");
  TCanvas *MatchedDMFindCJdRCanvasKyleOld     = (TCanvas*)infileKyleOld.Get("MatchedDMFindCJdR");
  TCanvas *MatchedDMFindRECOdRCanvasKyleOld   = (TCanvas*)infileKyleOld.Get("MatchedDMFindRECOdR");
  TCanvas *MatchedCJdRCanvasKyleOld   = (TCanvas*)infileKyleOld.Get("MatchedCJdR");
  TCanvas *MatchedRECOdRCanvasKyleOld = (TCanvas*)infileKyleOld.Get("MatchedRECOdR");

  TCanvas *MatchedLooseIsoCJPtGenCanvasKyleOld   = (TCanvas*)infileKyleOld.Get("MatchedLooseIsoCJPtGen");
  TCanvas *MatchedLooseIsoRECOPtGenCanvasKyleOld = (TCanvas*)infileKyleOld.Get("MatchedLooseIsoRECOPtGen");
  TCanvas *MatchedMedIsoCJPtGenCanvasKyleOld     = (TCanvas*)infileKyleOld.Get("MatchedMedIsoCJPtGen");
  TCanvas *MatchedMedIsoRECOPtGenCanvasKyleOld   = (TCanvas*)infileKyleOld.Get("MatchedMedIsoRECOPtGen");
  TCanvas *MatchedTightIsoCJPtGenCanvasKyleOld   = (TCanvas*)infileKyleOld.Get("MatchedTightIsoCJPtGen");
  TCanvas *MatchedTightIsoRECOPtGenCanvasKyleOld = (TCanvas*)infileKyleOld.Get("MatchedTightIsoRECOPtGen");
  TCanvas *MatchedDMFindCJPtGenCanvasKyleOld     = (TCanvas*)infileKyleOld.Get("MatchedDMFindCJPtGen");
  TCanvas *MatchedDMFindRECOPtGenCanvasKyleOld   = (TCanvas*)infileKyleOld.Get("MatchedDMFindRECOPtGen");
  TCanvas *MatchedCJPtGenCanvasKyleOld   = (TCanvas*)infileKyleOld.Get("MatchedCJPtGen");
  TCanvas *MatchedRECOPtGenCanvasKyleOld = (TCanvas*)infileKyleOld.Get("MatchedRECOPtGen");

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

  // Getting the histograms for the New Gen Matching and  DM Old
  TH1F* MatchedLooseIsoCJPtRachOld_   = (TH1F*)MatchedLooseIsoCJPtCanvasRachOld->GetPrimitive("MatchedLooseIsoCJPt");
  TH1F* MatchedLooseIsoRECOPtRachOld_ = (TH1F*)MatchedLooseIsoRECOPtCanvasRachOld->GetPrimitive("MatchedLooseIsoRECOPt");
  TH1F* MatchedMedIsoCJPtRachOld_     = (TH1F*)MatchedMedIsoCJPtCanvasRachOld->GetPrimitive("MatchedMedIsoCJPt");
  TH1F* MatchedMedIsoRECOPtRachOld_   = (TH1F*)MatchedMedIsoRECOPtCanvasRachOld->GetPrimitive("MatchedMedIsoRECOPt");
  TH1F* MatchedTightIsoCJPtRachOld_   = (TH1F*)MatchedTightIsoCJPtCanvasRachOld->GetPrimitive("MatchedTightIsoCJPt");
  TH1F* MatchedTightIsoRECOPtRachOld_ = (TH1F*)MatchedTightIsoRECOPtCanvasRachOld->GetPrimitive("MatchedTightIsoRECOPt");
  TH1F* MatchedDMFindCJPtRachOld_     = (TH1F*)MatchedDMFindCJPtCanvasRachOld->GetPrimitive("MatchedDMFindCJPt");
  TH1F* MatchedDMFindRECOPtRachOld_   = (TH1F*)MatchedDMFindRECOPtCanvasRachOld->GetPrimitive("MatchedDMFindRECOPt");
  TH1F* MatchedCJPtRachOld_   = (TH1F*)MatchedCJPtCanvasRachOld->GetPrimitive("MatchedCJPt");
  TH1F* MatchedRECOPtRachOld_ = (TH1F*)MatchedRECOPtCanvasRachOld->GetPrimitive("MatchedRECOPt");

  TH1F* MatchedLooseIsoCJdRRachOld_   = (TH1F*)MatchedLooseIsoCJdRCanvasRachOld->GetPrimitive("MatchedLooseIsoCJdR");
  TH1F* MatchedLooseIsoRECOdRRachOld_ = (TH1F*)MatchedLooseIsoRECOdRCanvasRachOld->GetPrimitive("MatchedLooseIsoRECOdR");
  TH1F* MatchedMedIsoCJdRRachOld_     = (TH1F*)MatchedMedIsoCJdRCanvasRachOld->GetPrimitive("MatchedMedIsoCJdR");
  TH1F* MatchedMedIsoRECOdRRachOld_   = (TH1F*)MatchedMedIsoRECOdRCanvasRachOld->GetPrimitive("MatchedMedIsoRECOdR");
  TH1F* MatchedTightIsoCJdRRachOld_   = (TH1F*)MatchedTightIsoCJdRCanvasRachOld->GetPrimitive("MatchedTightIsoCJdR");
  TH1F* MatchedTightIsoRECOdRRachOld_ = (TH1F*)MatchedTightIsoRECOdRCanvasRachOld->GetPrimitive("MatchedTightIsoRECOdR");
  TH1F* MatchedDMFindCJdRRachOld_     = (TH1F*)MatchedDMFindCJdRCanvasRachOld->GetPrimitive("MatchedDMFindCJdR");
  TH1F* MatchedDMFindRECOdRRachOld_   = (TH1F*)MatchedDMFindRECOdRCanvasRachOld->GetPrimitive("MatchedDMFindRECOdR");
  TH1F* MatchedCJdRRachOld_   = (TH1F*)MatchedCJdRCanvasRachOld->GetPrimitive("MatchedCJdR");
  TH1F* MatchedRECOdRRachOld_ = (TH1F*)MatchedRECOdRCanvasRachOld->GetPrimitive("MatchedRECOdR");

  TH1F* MatchedLooseIsoCJPtGenRachOld_   = (TH1F*)MatchedLooseIsoCJPtGenCanvasRachOld->GetPrimitive("MatchedLooseIsoCJPtGen");
  TH1F* MatchedLooseIsoRECOPtGenRachOld_ = (TH1F*)MatchedLooseIsoRECOPtGenCanvasRachOld->GetPrimitive("MatchedLooseIsoRECOPtGen");
  TH1F* MatchedMedIsoCJPtGenRachOld_     = (TH1F*)MatchedMedIsoCJPtGenCanvasRachOld->GetPrimitive("MatchedMedIsoCJPtGen");
  TH1F* MatchedMedIsoRECOPtGenRachOld_   = (TH1F*)MatchedMedIsoRECOPtGenCanvasRachOld->GetPrimitive("MatchedMedIsoRECOPtGen");
  TH1F* MatchedTightIsoCJPtGenRachOld_   = (TH1F*)MatchedTightIsoCJPtGenCanvasRachOld->GetPrimitive("MatchedTightIsoCJPtGen");
  TH1F* MatchedTightIsoRECOPtGenRachOld_ = (TH1F*)MatchedTightIsoRECOPtGenCanvasRachOld->GetPrimitive("MatchedTightIsoRECOPtGen");
  TH1F* MatchedDMFindCJPtGenRachOld_     = (TH1F*)MatchedDMFindCJPtGenCanvasRachOld->GetPrimitive("MatchedDMFindCJPtGen");
  TH1F* MatchedDMFindRECOPtGenRachOld_   = (TH1F*)MatchedDMFindRECOPtGenCanvasRachOld->GetPrimitive("MatchedDMFindRECOPtGen");
  TH1F* MatchedCJPtGenRachOld_   = (TH1F*)MatchedCJPtGenCanvasRachOld->GetPrimitive("MatchedCJPtGen");
  TH1F* MatchedRECOPtGenRachOld_ = (TH1F*)MatchedRECOPtGenCanvasRachOld->GetPrimitive("MatchedRECOPtGen");


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


  // Getting the histograms for the DM Old
  TH1F* MatchedLooseIsoCJPtKyleOld_   = (TH1F*)MatchedLooseIsoCJPtCanvasKyleOld->GetPrimitive("MatchedLooseIsoCJPt");
  TH1F* MatchedLooseIsoRECOPtKyleOld_ = (TH1F*)MatchedLooseIsoRECOPtCanvasKyleOld->GetPrimitive("MatchedLooseIsoRECOPt");
  TH1F* MatchedMedIsoCJPtKyleOld_     = (TH1F*)MatchedMedIsoCJPtCanvasKyleOld->GetPrimitive("MatchedMedIsoCJPt");
  TH1F* MatchedMedIsoRECOPtKyleOld_   = (TH1F*)MatchedMedIsoRECOPtCanvasKyleOld->GetPrimitive("MatchedMedIsoRECOPt");
  TH1F* MatchedTightIsoCJPtKyleOld_   = (TH1F*)MatchedTightIsoCJPtCanvasKyleOld->GetPrimitive("MatchedTightIsoCJPt");
  TH1F* MatchedTightIsoRECOPtKyleOld_ = (TH1F*)MatchedTightIsoRECOPtCanvasKyleOld->GetPrimitive("MatchedTightIsoRECOPt");
  TH1F* MatchedDMFindCJPtKyleOld_     = (TH1F*)MatchedDMFindCJPtCanvasKyleOld->GetPrimitive("MatchedDMFindCJPt");
  TH1F* MatchedDMFindRECOPtKyleOld_   = (TH1F*)MatchedDMFindRECOPtCanvasKyleOld->GetPrimitive("MatchedDMFindRECOPt");
  TH1F* MatchedCJPtKyleOld_   = (TH1F*)MatchedCJPtCanvasKyleOld->GetPrimitive("MatchedCJPt");
  TH1F* MatchedRECOPtKyleOld_ = (TH1F*)MatchedRECOPtCanvasKyleOld->GetPrimitive("MatchedRECOPt");

  TH1F* MatchedLooseIsoCJdRKyleOld_   = (TH1F*)MatchedLooseIsoCJdRCanvasKyleOld->GetPrimitive("MatchedLooseIsoCJdR");
  TH1F* MatchedLooseIsoRECOdRKyleOld_ = (TH1F*)MatchedLooseIsoRECOdRCanvasKyleOld->GetPrimitive("MatchedLooseIsoRECOdR");
  TH1F* MatchedMedIsoCJdRKyleOld_     = (TH1F*)MatchedMedIsoCJdRCanvasKyleOld->GetPrimitive("MatchedMedIsoCJdR");
  TH1F* MatchedMedIsoRECOdRKyleOld_   = (TH1F*)MatchedMedIsoRECOdRCanvasKyleOld->GetPrimitive("MatchedMedIsoRECOdR");
  TH1F* MatchedTightIsoCJdRKyleOld_   = (TH1F*)MatchedTightIsoCJdRCanvasKyleOld->GetPrimitive("MatchedTightIsoCJdR");
  TH1F* MatchedTightIsoRECOdRKyleOld_ = (TH1F*)MatchedTightIsoRECOdRCanvasKyleOld->GetPrimitive("MatchedTightIsoRECOdR");
  TH1F* MatchedDMFindCJdRKyleOld_     = (TH1F*)MatchedDMFindCJdRCanvasKyleOld->GetPrimitive("MatchedDMFindCJdR");
  TH1F* MatchedDMFindRECOdRKyleOld_   = (TH1F*)MatchedDMFindRECOdRCanvasKyleOld->GetPrimitive("MatchedDMFindRECOdR");
  TH1F* MatchedCJdRKyleOld_   = (TH1F*)MatchedCJdRCanvasKyleOld->GetPrimitive("MatchedCJdR");
  TH1F* MatchedRECOdRKyleOld_ = (TH1F*)MatchedRECOdRCanvasKyleOld->GetPrimitive("MatchedRECOdR");

  TH1F* MatchedLooseIsoCJPtGenKyleOld_   = (TH1F*)MatchedLooseIsoCJPtGenCanvasKyleOld->GetPrimitive("MatchedLooseIsoCJPtGen");
  TH1F* MatchedLooseIsoRECOPtGenKyleOld_ = (TH1F*)MatchedLooseIsoRECOPtGenCanvasKyleOld->GetPrimitive("MatchedLooseIsoRECOPtGen");
  TH1F* MatchedMedIsoCJPtGenKyleOld_     = (TH1F*)MatchedMedIsoCJPtGenCanvasKyleOld->GetPrimitive("MatchedMedIsoCJPtGen");
  TH1F* MatchedMedIsoRECOPtGenKyleOld_   = (TH1F*)MatchedMedIsoRECOPtGenCanvasKyleOld->GetPrimitive("MatchedMedIsoRECOPtGen");
  TH1F* MatchedTightIsoCJPtGenKyleOld_   = (TH1F*)MatchedTightIsoCJPtGenCanvasKyleOld->GetPrimitive("MatchedTightIsoCJPtGen");
  TH1F* MatchedTightIsoRECOPtGenKyleOld_ = (TH1F*)MatchedTightIsoRECOPtGenCanvasKyleOld->GetPrimitive("MatchedTightIsoRECOPtGen");
  TH1F* MatchedDMFindCJPtGenKyleOld_     = (TH1F*)MatchedDMFindCJPtGenCanvasKyleOld->GetPrimitive("MatchedDMFindCJPtGen");
  TH1F* MatchedDMFindRECOPtGenKyleOld_   = (TH1F*)MatchedDMFindRECOPtGenCanvasKyleOld->GetPrimitive("MatchedDMFindRECOPtGen");
  TH1F* MatchedCJPtGenKyleOld_   = (TH1F*)MatchedCJPtGenCanvasKyleOld->GetPrimitive("MatchedCJPtGen");
  TH1F* MatchedRECOPtGenKyleOld_ = (TH1F*)MatchedRECOPtGenCanvasKyleOld->GetPrimitive("MatchedRECOPtGen");

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

  // To divide histograms for efficiency with new Gen Matching Old DM's
  TGraphAsymmErrors* FinalEffLooseIsoRECOPtRachOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoRECOPtRachOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoRECOPtRachOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindRECOPtRachOld_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoCJPtRachOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoCJPtRachOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoCJPtRachOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindCJPtRachOld_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoRECOdRRachOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoRECOdRRachOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoRECOdRRachOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindRECOdRRachOld_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoCJdRRachOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoCJdRRachOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoCJdRRachOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindCJdRRachOld_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoRECOPtGenRachOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoRECOPtGenRachOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoRECOPtGenRachOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindRECOPtGenRachOld_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoCJPtGenRachOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoCJPtGenRachOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoCJPtGenRachOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindCJPtGenRachOld_ = new TGraphAsymmErrors(30);


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

  // To divide histograms for efficiency with Old Gen Matching Old DM's
  TGraphAsymmErrors* FinalEffLooseIsoRECOPtKyleOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoRECOPtKyleOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoRECOPtKyleOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindRECOPtKyleOld_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoCJPtKyleOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoCJPtKyleOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoCJPtKyleOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindCJPtKyleOld_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoRECOdRKyleOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoRECOdRKyleOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoRECOdRKyleOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindRECOdRKyleOld_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoCJdRKyleOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoCJdRKyleOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoCJdRKyleOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindCJdRKyleOld_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoRECOPtGenKyleOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoRECOPtGenKyleOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoRECOPtGenKyleOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindRECOPtGenKyleOld_ = new TGraphAsymmErrors(30);

  TGraphAsymmErrors* FinalEffLooseIsoCJPtGenKyleOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffMedIsoCJPtGenKyleOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffTightIsoCJPtGenKyleOld_ = new TGraphAsymmErrors(30);
  TGraphAsymmErrors* FinalEffDMFindCJPtGenKyleOld_ = new TGraphAsymmErrors(30);

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



  // To divide histograms for efficiency with new Gen Matching Old DM's
  FinalEffLooseIsoRECOPtRachOld_->Divide(MatchedLooseIsoRECOPtRachOld_, MatchedRECOPtRachOld_);
  FinalEffMedIsoRECOPtRachOld_->Divide(MatchedMedIsoRECOPtRachOld_,     MatchedRECOPtRachOld_);
  FinalEffTightIsoRECOPtRachOld_->Divide(MatchedTightIsoRECOPtRachOld_, MatchedRECOPtRachOld_);
  FinalEffDMFindRECOPtRachOld_->Divide(MatchedDMFindRECOPtRachOld_,     MatchedRECOPtRachOld_);
  FinalEffLooseIsoCJPtRachOld_->Divide(MatchedLooseIsoCJPtRachOld_, MatchedCJPtRachOld_);
  FinalEffMedIsoCJPtRachOld_->Divide(MatchedMedIsoCJPtRachOld_,     MatchedCJPtRachOld_);
  FinalEffTightIsoCJPtRachOld_->Divide(MatchedTightIsoCJPtRachOld_, MatchedCJPtRachOld_);
  FinalEffDMFindCJPtRachOld_->Divide(MatchedDMFindCJPtRachOld_,     MatchedCJPtRachOld_);

  FinalEffLooseIsoRECOdRRachOld_->Divide(MatchedLooseIsoRECOdRRachOld_, MatchedRECOdRRachOld_);
  FinalEffMedIsoRECOdRRachOld_->Divide(MatchedMedIsoRECOdRRachOld_,     MatchedRECOdRRachOld_);
  FinalEffTightIsoRECOdRRachOld_->Divide(MatchedTightIsoRECOdRRachOld_, MatchedRECOdRRachOld_);
  FinalEffDMFindRECOdRRachOld_->Divide(MatchedDMFindRECOdRRachOld_,     MatchedRECOdRRachOld_);
  FinalEffLooseIsoCJdRRachOld_->Divide(MatchedLooseIsoCJdRRachOld_, MatchedCJdRRachOld_);
  FinalEffMedIsoCJdRRachOld_->Divide(MatchedMedIsoCJdRRachOld_,     MatchedCJdRRachOld_);
  FinalEffTightIsoCJdRRachOld_->Divide(MatchedTightIsoCJdRRachOld_, MatchedCJdRRachOld_);
  FinalEffDMFindCJdRRachOld_->Divide(MatchedDMFindCJdRRachOld_,     MatchedCJdRRachOld_);

  FinalEffLooseIsoRECOPtGenRachOld_->Divide(MatchedLooseIsoRECOPtGenRachOld_, MatchedRECOPtGenRachOld_);
  FinalEffMedIsoRECOPtGenRachOld_->Divide(MatchedMedIsoRECOPtGenRachOld_,     MatchedRECOPtGenRachOld_);
  FinalEffTightIsoRECOPtGenRachOld_->Divide(MatchedTightIsoRECOPtGenRachOld_, MatchedRECOPtGenRachOld_);
  FinalEffDMFindRECOPtGenRachOld_->Divide(MatchedDMFindRECOPtGenRachOld_,     MatchedRECOPtGenRachOld_);
  FinalEffLooseIsoCJPtGenRachOld_->Divide(MatchedLooseIsoCJPtGenRachOld_, MatchedCJPtGenRachOld_);
  FinalEffMedIsoCJPtGenRachOld_->Divide(MatchedMedIsoCJPtGenRachOld_,     MatchedCJPtGenRachOld_);
  FinalEffTightIsoCJPtGenRachOld_->Divide(MatchedTightIsoCJPtGenRachOld_, MatchedCJPtGenRachOld_);
  FinalEffDMFindCJPtGenRachOld_->Divide(MatchedDMFindCJPtGenRachOld_,     MatchedCJPtGenRachOld_);

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

  // To divide histograms for efficiency with Old Gen Matching Old DM's 
  FinalEffLooseIsoRECOPtKyleOld_->Divide(MatchedLooseIsoRECOPtKyleOld_, MatchedRECOPtKyleOld_);
  FinalEffMedIsoRECOPtKyleOld_->Divide(MatchedMedIsoRECOPtKyleOld_,     MatchedRECOPtKyleOld_);
  FinalEffTightIsoRECOPtKyleOld_->Divide(MatchedTightIsoRECOPtKyleOld_, MatchedRECOPtKyleOld_);
  FinalEffDMFindRECOPtKyleOld_->Divide(MatchedDMFindRECOPtKyleOld_,     MatchedRECOPtKyleOld_);
  FinalEffLooseIsoCJPtKyleOld_->Divide(MatchedLooseIsoCJPtKyleOld_, MatchedCJPtKyleOld_);
  FinalEffMedIsoCJPtKyleOld_->Divide(MatchedMedIsoCJPtKyleOld_,     MatchedCJPtKyleOld_);
  FinalEffTightIsoCJPtKyleOld_->Divide(MatchedTightIsoCJPtKyleOld_, MatchedCJPtKyleOld_);
  FinalEffDMFindCJPtKyleOld_->Divide(MatchedDMFindCJPtKyleOld_,     MatchedCJPtKyleOld_);
  
  FinalEffLooseIsoRECOdRKyleOld_->Divide(MatchedLooseIsoRECOdRKyleOld_, MatchedRECOdRKyleOld_);
  FinalEffMedIsoRECOdRKyleOld_->Divide(MatchedMedIsoRECOdRKyleOld_,     MatchedRECOdRKyleOld_);
  FinalEffTightIsoRECOdRKyleOld_->Divide(MatchedTightIsoRECOdRKyleOld_, MatchedRECOdRKyleOld_);
  FinalEffDMFindRECOdRKyleOld_->Divide(MatchedDMFindRECOdRKyleOld_,     MatchedRECOdRKyleOld_);
  FinalEffLooseIsoCJdRKyleOld_->Divide(MatchedLooseIsoCJdRKyleOld_, MatchedCJdRKyleOld_);
  FinalEffMedIsoCJdRKyleOld_->Divide(MatchedMedIsoCJdRKyleOld_,     MatchedCJdRKyleOld_);
  FinalEffTightIsoCJdRKyleOld_->Divide(MatchedTightIsoCJdRKyleOld_, MatchedCJdRKyleOld_);
  FinalEffDMFindCJdRKyleOld_->Divide(MatchedDMFindCJdRKyleOld_,     MatchedCJdRKyleOld_);
  
  FinalEffLooseIsoRECOPtGenKyleOld_->Divide(MatchedLooseIsoRECOPtGenKyleOld_, MatchedRECOPtGenKyleOld_);
  FinalEffMedIsoRECOPtGenKyleOld_->Divide(MatchedMedIsoRECOPtGenKyleOld_,     MatchedRECOPtGenKyleOld_);
  FinalEffTightIsoRECOPtGenKyleOld_->Divide(MatchedTightIsoRECOPtGenKyleOld_, MatchedRECOPtGenKyleOld_);
  FinalEffDMFindRECOPtGenKyleOld_->Divide(MatchedDMFindRECOPtGenKyleOld_,     MatchedRECOPtGenKyleOld_);
  FinalEffLooseIsoCJPtGenKyleOld_->Divide(MatchedLooseIsoCJPtGenKyleOld_, MatchedCJPtGenKyleOld_);
  FinalEffMedIsoCJPtGenKyleOld_->Divide(MatchedMedIsoCJPtGenKyleOld_,     MatchedCJPtGenKyleOld_);
  FinalEffTightIsoCJPtGenKyleOld_->Divide(MatchedTightIsoCJPtGenKyleOld_, MatchedCJPtGenKyleOld_);
  FinalEffDMFindCJPtGenKyleOld_->Divide(MatchedDMFindCJPtGenKyleOld_,     MatchedCJPtGenKyleOld_);
  
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

  //Set Colors for New Gen Matching Old DM's
  setGraphOptions(FinalEffLooseIsoCJPtRachOld_, kRed, 0.7, 10, "p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffLooseIsoRECOPtRachOld_, kBlack, 0.7, 10, "p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoCJPtRachOld_, kRed, 0.7, 10, "p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoRECOPtRachOld_, kBlack, 0.7, 10, "p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoCJPtRachOld_, kRed, 0.7, 10, "p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoRECOPtRachOld_, kBlack, 0.7, 10, "p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindCJPtRachOld_, kRed, 0.7, 10, "p_{T}", "#epsilon (DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindRECOPtRachOld_, kBlack, 0.7, 10, "p_{T}", "#epsilon (DMFinding + GM / GM)");

  setGraphOptions(FinalEffLooseIsoCJdRRachOld_, kRed, 0.7, 10, "#DeltaR", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffLooseIsoRECOdRRachOld_, kBlack, 0.7, 10, "#DeltaR", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoCJdRRachOld_, kRed, 0.7, 10, "#DeltaR", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoRECOdRRachOld_, kBlack, 0.7, 10, "#DeltaR", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoCJdRRachOld_, kRed, 0.7, 10, "#DeltaR", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoRECOdRRachOld_, kBlack, 0.7, 10, "#DeltaR", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindCJdRRachOld_, kRed, 0.7, 10, "#DeltaR", "#epsilon (DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindRECOdRRachOld_, kBlack, 0.7, 10, "#DeltaR", "#epsilon (DMFinding + GM / GM)");

  setGraphOptions(FinalEffLooseIsoCJPtGenRachOld_, kRed, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffLooseIsoRECOPtGenRachOld_, kBlack, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoCJPtGenRachOld_, kRed, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoRECOPtGenRachOld_, kBlack, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoCJPtGenRachOld_, kRed, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoRECOPtGenRachOld_, kBlack, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindCJPtGenRachOld_, kRed, 0.7, 10, "Gen Visible p_{T}", "#epsilon (DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindRECOPtGenRachOld_, kBlack, 0.7, 10, "Gen Visible p_{T}", "#epsilon (DMFinding + GM / GM)");

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

  //Set Colors for New Gen Matching Old DM's
  setGraphOptions(FinalEffLooseIsoCJPtKyleOld_, kAzure+10, 0.7, 10, "p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffLooseIsoRECOPtKyleOld_, kOrange+7, 0.7, 10, "p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoCJPtKyleOld_, kAzure+10, 0.7, 10, "p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoRECOPtKyleOld_, kOrange+7, 0.7, 10, "p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoCJPtKyleOld_, kAzure+10, 0.7, 10, "p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoRECOPtKyleOld_, kOrange+7, 0.7, 10, "p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindCJPtKyleOld_, kAzure+10, 0.7, 10, "p_{T}", "#epsilon (DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindRECOPtKyleOld_, kOrange+7, 0.7, 10, "p_{T}", "#epsilon (DMFinding + GM / GM)");

  setGraphOptions(FinalEffLooseIsoCJdRKyleOld_, kAzure+10, 0.7, 10, "#DeltaR", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffLooseIsoRECOdRKyleOld_, kOrange+7, 0.7, 10, "#DeltaR", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoCJdRKyleOld_, kAzure+10, 0.7, 10, "#DeltaR", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoRECOdRKyleOld_, kOrange+7, 0.7, 10, "#DeltaR", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoCJdRKyleOld_, kAzure+10, 0.7, 10, "#DeltaR", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoRECOdRKyleOld_, kOrange+7, 0.7, 10, "#DeltaR", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindCJdRKyleOld_, kAzure+10, 0.7, 10, "#DeltaR", "#epsilon (DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindRECOdRKyleOld_, kOrange+7, 0.7, 10, "#DeltaR", "#epsilon (DMFinding + GM / GM)");

  setGraphOptions(FinalEffLooseIsoCJPtGenKyleOld_, kAzure+10, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffLooseIsoRECOPtGenKyleOld_, kOrange+7, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Loose Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoCJPtGenKyleOld_, kAzure+10, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffMedIsoRECOPtGenKyleOld_, kOrange+7, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Med Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoCJPtGenKyleOld_, kAzure+10, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffTightIsoRECOPtGenKyleOld_, kOrange+7, 0.7, 10, "Gen Visible p_{T}", "#epsilon (Tight Iso + DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindCJPtGenKyleOld_, kAzure+10, 0.7, 10, "Gen Visible p_{T}", "#epsilon (DMFinding + GM / GM)");
  setGraphOptions(FinalEffDMFindRECOPtGenKyleOld_, kOrange+7, 0.7, 10, "Gen Visible p_{T}", "#epsilon (DMFinding + GM / GM)");

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
  leg->AddEntry(FinalEffDMFindRECOPtGenRachOld_, "No  CJ,New Gen,Old DM","L");
  leg->AddEntry(FinalEffDMFindRECOPtGenRachNew_, "No  CJ,New Gen,New DM","L");
  leg->AddEntry(FinalEffDMFindRECOPtGenKyleOld_, "No  CJ,Old Gen,Old DM","L");
  leg->AddEntry(FinalEffDMFindRECOPtGenKyleNew_, "No  CJ,Old Gen,New DM","L");
  leg->AddEntry(FinalEffDMFindCJPtGenRachOld_,   "Yes CJ,New Gen,Old DM","L");
  leg->AddEntry(FinalEffDMFindCJPtGenRachNew_,   "Yes CJ,New Gen,New DM","L");
  leg->AddEntry(FinalEffDMFindCJPtGenKyleOld_,   "Yes CJ,Old Gen,Old DM","L");
  leg->AddEntry(FinalEffDMFindCJPtGenKyleNew_,   "Yes CJ,Old Gen,New DM","L");

  legMatch = new TLegend(0.1,0.7,0.25,0.9);
  legMatch->AddEntry(FinalEffDMFindRECOPtGenRachNew_, "No  CJ,New Gen,New DM","L");
  legMatch->AddEntry(FinalEffDMFindCJPtGenRachNew_,   "Yes CJ,New Gen,New DM","L");
  legMatch->AddEntry(FinalEffDMFindRECOPtGenRachOld_, "No  CJ,New Gen,Old DM","L");
  legMatch->AddEntry(FinalEffDMFindCJPtGenRachOld_,   "Yes CJ,New Gen,Old DM","L");



  LoosePtCanvas.cd();
  FinalEffLooseIsoCJPtRachOld_->Draw();
  FinalEffLooseIsoRECOPtRachOld_->Draw("SAME");
  FinalEffLooseIsoCJPtRachNew_->Draw("SAME");
  FinalEffLooseIsoRECOPtRachNew_->Draw("SAME");
  FinalEffLooseIsoCJPtKyleOld_->Draw("SAME");
  FinalEffLooseIsoRECOPtKyleOld_->Draw("SAME");
  FinalEffLooseIsoCJPtKyleNew_->Draw("SAME");
  FinalEffLooseIsoRECOPtKyleNew_->Draw("SAME");
  leg->Draw();

  LoosePtNewGenMatchCanvas.cd();
  FinalEffLooseIsoCJPtRachOld_->Draw();
  FinalEffLooseIsoRECOPtRachOld_->Draw("SAME");
  FinalEffLooseIsoCJPtRachNew_->Draw("SAME");
  FinalEffLooseIsoRECOPtRachNew_->Draw("SAME");
  legMatch->Draw();

  MedPtCanvas.cd();
  FinalEffMedIsoCJPtRachOld_->Draw();
  FinalEffMedIsoRECOPtRachOld_->Draw("SAME");
  FinalEffMedIsoCJPtRachNew_->Draw("SAME");
  FinalEffMedIsoRECOPtRachNew_->Draw("SAME");
  FinalEffMedIsoCJPtKyleOld_->Draw("SAME");
  FinalEffMedIsoRECOPtKyleOld_->Draw("SAME");
  FinalEffMedIsoCJPtKyleNew_->Draw("SAME");
  FinalEffMedIsoRECOPtKyleNew_->Draw("SAME");
  leg->Draw();

  MedPtNewGenMatchCanvas.cd();
  FinalEffMedIsoCJPtRachOld_->Draw();
  FinalEffMedIsoRECOPtRachOld_->Draw("SAME");
  FinalEffMedIsoCJPtRachNew_->Draw("SAME");
  FinalEffMedIsoRECOPtRachNew_->Draw("SAME");
  legMatch->Draw();
  

  TightPtCanvas.cd();
  FinalEffTightIsoCJPtRachOld_->Draw();
  FinalEffTightIsoRECOPtRachOld_->Draw("SAME");
  FinalEffTightIsoCJPtRachNew_->Draw("SAME");
  FinalEffTightIsoRECOPtRachNew_->Draw("SAME");
  FinalEffTightIsoCJPtKyleOld_->Draw("SAME");
  FinalEffTightIsoRECOPtKyleOld_->Draw("SAME");
  FinalEffTightIsoCJPtKyleNew_->Draw("SAME");
  FinalEffTightIsoRECOPtKyleNew_->Draw("SAME");
  leg->Draw();

  TightPtNewGenMatchCanvas.cd();
  FinalEffTightIsoCJPtRachOld_->Draw();
  FinalEffTightIsoRECOPtRachOld_->Draw("SAME");
  FinalEffTightIsoCJPtRachNew_->Draw("SAME");
  FinalEffTightIsoRECOPtRachNew_->Draw("SAME");
  legMatch->Draw();
  

  DMPtCanvas.cd();
  FinalEffDMFindCJPtRachOld_->Draw();
  FinalEffDMFindRECOPtRachOld_->Draw("SAME");
  FinalEffDMFindCJPtRachNew_->Draw("SAME");
  FinalEffDMFindRECOPtRachNew_->Draw("SAME");
  FinalEffDMFindCJPtKyleOld_->Draw("SAME");
  FinalEffDMFindRECOPtKyleOld_->Draw("SAME");
  FinalEffDMFindCJPtKyleNew_->Draw("SAME");
  FinalEffDMFindRECOPtKyleNew_->Draw("SAME");
  leg->Draw();

  DMPtNewGenMatchCanvas.cd();
  FinalEffDMFindCJPtRachOld_->Draw();
  FinalEffDMFindRECOPtRachOld_->Draw("SAME");
  FinalEffDMFindCJPtRachNew_->Draw("SAME");
  FinalEffDMFindRECOPtRachNew_->Draw("SAME");
  legMatch->Draw();
  


  LoosedRCanvas.cd();
  FinalEffLooseIsoCJdRRachOld_->Draw();
  FinalEffLooseIsoRECOdRRachOld_->Draw("SAME");
  FinalEffLooseIsoCJdRRachNew_->Draw("SAME");
  FinalEffLooseIsoRECOdRRachNew_->Draw("SAME");
  FinalEffLooseIsoCJdRKyleOld_->Draw("SAME");
  FinalEffLooseIsoRECOdRKyleOld_->Draw("SAME");
  FinalEffLooseIsoCJdRKyleNew_->Draw("SAME");
  FinalEffLooseIsoRECOdRKyleNew_->Draw("SAME");
  leg->Draw();

  LoosedRNewGenMatchCanvas.cd();
  FinalEffLooseIsoCJdRRachOld_->Draw();
  FinalEffLooseIsoRECOdRRachOld_->Draw("SAME");
  FinalEffLooseIsoCJdRRachNew_->Draw("SAME");
  FinalEffLooseIsoRECOdRRachNew_->Draw("SAME");
  legMatch->Draw();
  
 
  MeddRCanvas.cd();
  FinalEffMedIsoCJdRRachOld_->Draw();
  FinalEffMedIsoRECOdRRachOld_->Draw("SAME");
  FinalEffMedIsoCJdRRachNew_->Draw("SAME");
  FinalEffMedIsoRECOdRRachNew_->Draw("SAME");
  FinalEffMedIsoCJdRKyleOld_->Draw("SAME");
  FinalEffMedIsoRECOdRKyleOld_->Draw("SAME");
  FinalEffMedIsoCJdRKyleNew_->Draw("SAME");
  FinalEffMedIsoRECOdRKyleNew_->Draw("SAME");
  leg->Draw();
 
  MeddRNewGenMatchCanvas.cd();
  FinalEffMedIsoCJdRRachOld_->Draw();
  FinalEffMedIsoRECOdRRachOld_->Draw("SAME");
  FinalEffMedIsoCJdRRachNew_->Draw("SAME");
  FinalEffMedIsoRECOdRRachNew_->Draw("SAME");
  legMatch->Draw();
 
  TightdRCanvas.cd();
  FinalEffTightIsoCJdRRachOld_->Draw();
  FinalEffTightIsoRECOdRRachOld_->Draw("SAME");
  FinalEffTightIsoCJdRRachNew_->Draw("SAME");
  FinalEffTightIsoRECOdRRachNew_->Draw("SAME");
  FinalEffTightIsoCJdRKyleOld_->Draw("SAME");
  FinalEffTightIsoRECOdRKyleOld_->Draw("SAME");
  FinalEffTightIsoCJdRKyleNew_->Draw("SAME");
  FinalEffTightIsoRECOdRKyleNew_->Draw("SAME");
  leg->Draw();

  TightdRNewGenMatchCanvas.cd();
  FinalEffTightIsoCJdRRachOld_->Draw();
  FinalEffTightIsoRECOdRRachOld_->Draw("SAME");
  FinalEffTightIsoCJdRRachNew_->Draw("SAME");
  FinalEffTightIsoRECOdRRachNew_->Draw("SAME");
  legMatch->Draw();

  DMdRCanvas.cd();
  FinalEffDMFindCJdRRachOld_->Draw();
  FinalEffDMFindRECOdRRachOld_->Draw("SAME");
  FinalEffDMFindCJdRRachNew_->Draw("SAME");
  FinalEffDMFindRECOdRRachNew_->Draw("SAME");
  FinalEffDMFindCJdRKyleOld_->Draw("SAME");
  FinalEffDMFindRECOdRKyleOld_->Draw("SAME");
  FinalEffDMFindCJdRKyleNew_->Draw("SAME");
  FinalEffDMFindRECOdRKyleNew_->Draw("SAME");
  leg->Draw();

  DMdRNewGenMatchCanvas.cd();
  FinalEffDMFindCJdRRachOld_->Draw();
  FinalEffDMFindRECOdRRachOld_->Draw("SAME");
  FinalEffDMFindCJdRRachNew_->Draw("SAME");
  FinalEffDMFindRECOdRRachNew_->Draw("SAME");
  legMatch->Draw();
 
  LoosePtGenCanvas.cd();
  FinalEffLooseIsoCJPtGenRachOld_->Draw();
  FinalEffLooseIsoRECOPtGenRachOld_->Draw("SAME");
  FinalEffLooseIsoCJPtGenRachNew_->Draw("SAME");
  FinalEffLooseIsoRECOPtGenRachNew_->Draw("SAME");
  FinalEffLooseIsoCJPtGenKyleOld_->Draw("SAME");
  FinalEffLooseIsoRECOPtGenKyleOld_->Draw("SAME");
  FinalEffLooseIsoCJPtGenKyleNew_->Draw("SAME");
  FinalEffLooseIsoRECOPtGenKyleNew_->Draw("SAME");
  leg->Draw();

  LoosePtGenNewGenMatchCanvas.cd();
  FinalEffLooseIsoCJPtGenRachOld_->Draw();
  FinalEffLooseIsoRECOPtGenRachOld_->Draw("SAME");
  FinalEffLooseIsoCJPtGenRachNew_->Draw("SAME");
  FinalEffLooseIsoRECOPtGenRachNew_->Draw("SAME");
  legMatch->Draw();
 
  MedPtGenCanvas.cd();
  FinalEffMedIsoCJPtGenRachOld_->Draw();
  FinalEffMedIsoRECOPtGenRachOld_->Draw("SAME");
  FinalEffMedIsoCJPtGenRachNew_->Draw("SAME");
  FinalEffMedIsoRECOPtGenRachNew_->Draw("SAME");
  FinalEffMedIsoCJPtGenKyleOld_->Draw("SAME");
  FinalEffMedIsoRECOPtGenKyleOld_->Draw("SAME");
  FinalEffMedIsoCJPtGenKyleNew_->Draw("SAME");
  FinalEffMedIsoRECOPtGenKyleNew_->Draw("SAME");
  leg->Draw();

  MedPtGenNewGenMatchCanvas.cd();
  FinalEffMedIsoCJPtGenRachOld_->Draw();
  FinalEffMedIsoRECOPtGenRachOld_->Draw("SAME");
  FinalEffMedIsoCJPtGenRachNew_->Draw("SAME");
  FinalEffMedIsoRECOPtGenRachNew_->Draw("SAME");
  legMatch->Draw();
 
  TightPtGenCanvas.cd();
  FinalEffTightIsoCJPtGenRachOld_->Draw();
  FinalEffTightIsoRECOPtGenRachOld_->Draw("SAME");
  FinalEffTightIsoCJPtGenRachNew_->Draw("SAME");
  FinalEffTightIsoRECOPtGenRachNew_->Draw("SAME");
  FinalEffTightIsoCJPtGenKyleOld_->Draw("SAME");
  FinalEffTightIsoRECOPtGenKyleOld_->Draw("SAME");
  FinalEffTightIsoCJPtGenKyleNew_->Draw("SAME");
  FinalEffTightIsoRECOPtGenKyleNew_->Draw("SAME");
  leg->Draw();

  TightPtGenNewGenMatchCanvas.cd();
  FinalEffTightIsoCJPtGenRachOld_->Draw();
  FinalEffTightIsoRECOPtGenRachOld_->Draw("SAME");
  FinalEffTightIsoCJPtGenRachNew_->Draw("SAME");
  FinalEffTightIsoRECOPtGenRachNew_->Draw("SAME");
  legMatch->Draw();

  DMPtGenCanvas.cd();
  FinalEffDMFindCJPtGenRachOld_->Draw();
  FinalEffDMFindRECOPtGenRachOld_->Draw("SAME");
  FinalEffDMFindCJPtGenRachNew_->Draw("SAME");
  FinalEffDMFindRECOPtGenRachNew_->Draw("SAME");
  FinalEffDMFindCJPtGenKyleOld_->Draw("SAME");
  FinalEffDMFindRECOPtGenKyleOld_->Draw("SAME");
  FinalEffDMFindCJPtGenKyleNew_->Draw("SAME");
  FinalEffDMFindRECOPtGenKyleNew_->Draw("SAME");
  leg->Draw();

  DMPtGenNewGenMatchCanvas.cd();
  FinalEffDMFindCJPtGenRachOld_->Draw();
  FinalEffDMFindRECOPtGenRachOld_->Draw("SAME");
  FinalEffDMFindCJPtGenRachNew_->Draw("SAME");
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

  MatchedCJPtCanvasRachOld->Write();
  MatchedRECOPtCanvasRachOld->Write();
  MatchedCJdRCanvasRachOld->Write();
  MatchedRECOdRCanvasRachOld->Write();
  MatchedCJPtGenCanvasRachOld->Write();
  MatchedRECOPtGenCanvasRachOld->Write();

  MatchedCJPtCanvasRachNew->Write();
  MatchedRECOPtCanvasRachNew->Write();
  MatchedCJdRCanvasRachNew->Write();
  MatchedRECOdRCanvasRachNew->Write();
  MatchedCJPtGenCanvasRachNew->Write();
  MatchedRECOPtGenCanvasRachNew->Write();

  MatchedCJPtCanvasKyleOld->Write();
  MatchedRECOPtCanvasKyleOld->Write();
  MatchedCJdRCanvasKyleOld->Write();
  MatchedRECOdRCanvasKyleOld->Write();
  MatchedCJPtGenCanvasKyleOld->Write();
  MatchedRECOPtGenCanvasKyleOld->Write();

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


#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include "TCanvas.h"

void rootMacro_CleanJets_NTauDM()
{
  gStyle->SetOptStat(kFALSE);

  TFile infileCJ("/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/ggH125a9_GenTauDecayID_IndivCJ_NewDMFind_TauDMPlots_Scale_APR6/ggH125a9_GenTauDecayID_IndivCJ_NewDMFind_TauDMPlots_Scale_APR6_Plots.root");
  TFile infileRECO("/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/ggH125a9_GenTauDecayID_IndivRECO_NewDMFind_TauDMPlots_Scale_APR6/ggH125a9_GenTauDecayID_IndivRECO_NewDMFind_TauDMPlots_Scale_APR6_Plots.root");

  TFile *outFile = new TFile("combHist_CleanJets_h125a9_NTauDecayMode.root", "RECREATE");

cout << "Files Created" << endl;

  TCanvas *NTauDecayModeCJRecoCanvas = (TCanvas*)infileCJ.Get("NTauDecayModeRECO");
  TCanvas *NTauDecayModeRECORecoCanvas = (TCanvas*)infileRECO.Get("NTauDecayModeRECO");


  TCanvas *NTauDecayModeCJGenCanvas = (TCanvas*)infileCJ.Get("NTauDecayModeGEN");
  TCanvas *NTauDecayModeRECOGenCanvas = (TCanvas*)infileRECO.Get("NTauDecayModeGEN");

cout << "Got Canvases" << endl;

  TH1F* NTauDecayModeCJReco_ = (TH1F*)NTauDecayModeCJRecoCanvas->GetPrimitive("NTauDecayModeRECO");
  TH1F* NTauDecayModeRECOReco_ = (TH1F*)NTauDecayModeRECORecoCanvas->GetPrimitive("NTauDecayModeRECO");

  TH1F* NTauDecayModeCJGen_ = (TH1F*)NTauDecayModeCJGenCanvas->GetPrimitive("NTauDecayModeGEN");
  TH1F* NTauDecayModeRECOGen_ = (TH1F*)NTauDecayModeRECOGenCanvas->GetPrimitive("NTauDecayModeGEN");

cout << "Histograms assigned." << endl; 

  TCanvas NTauDecayModeRECO("NTauDecayModeRECO","",600,600);
  TCanvas NTauDecayModeGEN("NTauDecayModeGEN","",600,600);

cout << "Canvases created" << endl;

  NTauDecayModeCJReco_->SetLineColor(kBlack);
  NTauDecayModeRECOReco_->SetLineColor(kRed);
  
  NTauDecayModeCJGen_->SetLineColor(kBlack);
  NTauDecayModeRECOGen_->SetLineColor(kRed);
  
  NTauDecayModeCJReco_->SetMarkerColor(kBlack);
  NTauDecayModeRECOReco_->SetMarkerColor(kRed);
  
  NTauDecayModeCJGen_->SetMarkerColor(kBlack);
  NTauDecayModeRECOGen_->SetMarkerColor(kRed);
 
  NTauDecayModeCJReco_->GetXaxis()->SetTitle("Tau Decay Mode RECO");
  NTauDecayModeRECOReco_->GetXaxis()->SetTitle("Tau Decay Mode RECO");
 
  NTauDecayModeCJGen_->GetXaxis()->SetTitle("Tau Decay Mode Gen");
  NTauDecayModeRECOGen_->GetXaxis()->SetTitle("Tau Decay Mode GEN");
 
cout << "Attributes set." << endl;  

  leg = new TLegend(0.1,0.7,0.25,0.9);
  leg->AddEntry(NTauDecayModeCJReco_, "No CJ","L");
  leg->AddEntry(NTauDecayModeRECOReco_, "CJ","L");

  NTauDecayModeRECO.cd();
  NTauDecayModeCJReco_->Draw();
  NTauDecayModeRECOReco_->Draw("SAME");
  leg->Draw();

  NTauDecayModeGEN.cd();
  NTauDecayModeCJGen_->Draw();
  NTauDecayModeRECOGen_->Draw("SAME");
  leg->Draw();

cout << "Histograms Drawn" << endl;

  outFile->cd();
  NTauDecayModeRECO.Write();
  NTauDecayModeGEN.Write();
  outFile->Write();
  outFile->Close();

cout << "end" << endl;

}//rootMacro_BBA_combine

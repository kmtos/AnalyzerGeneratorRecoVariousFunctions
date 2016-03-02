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
#include "TMultiGraph.h"

void rootMacro_CleanJets()
{
  gStyle->SetOptStat(kFALSE);

  TFile infileCJ("/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_4_12_patch4/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/ggH125a9_GenTauDecayID_IndividualRCJ_FEB25/ggH125a9_GenTauDecayID_IndividualRCJ_FEB25_Plots.root");
  TFile infileRECO("/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_4_12_patch4/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/ggH125a9_GenTauDecayID_IndividualRECO_FEB25/ggH125a9_GenTauDecayID_IndividualRECO_FEB25_Plots.root");
  TFile *outFile = new TFile("combHist_CleanJets_h125a9.root", "RECREATE");

cout << "Files Created" << endl;

  TCanvas *MatchedLooseIsoCJPtCanvas = (TCanvas*)infileCJ.Get("MatchedLooseIsoCJPt");
  TCanvas *MatchedLooseIsoRECOPtCanvas = (TCanvas*)infileRECO.Get("MatchedLooseIsoRECOPt");
  TCanvas *MatchedMedIsoCJPtCanvas = (TCanvas*)infileCJ.Get("MatchedMedIsoCJPt");
  TCanvas *MatchedMedIsoRECOPtCanvas = (TCanvas*)infileRECO.Get("MatchedMedIsoRECOPt");
  TCanvas *MatchedTightIsoCJPtCanvas = (TCanvas*)infileCJ.Get("MatchedTightIsoCJPt");
  TCanvas *MatchedTightIsoRECOPtCanvas = (TCanvas*)infileRECO.Get("MatchedTightIsoRECOPt");
  TCanvas *MatchedDMFindCJPtCanvas = (TCanvas*)infileCJ.Get("MatchedDMFindCJPt");
  TCanvas *MatchedDMFindRECOPtCanvas = (TCanvas*)infileRECO.Get("MatchedDMFindRECOPt");
  TCanvas *MatchedCJPtCanvas = (TCanvas*)infileCJ.Get("MatchedCJPt");
  TCanvas *MatchedRECOPtCanvas = (TCanvas*)infileRECO.Get("MatchedRECOPt");


  TCanvas *MatchedLooseIsoCJdRCanvas = (TCanvas*)infileCJ.Get("MatchedLooseIsoCJdR");
  TCanvas *MatchedLooseIsoRECOdRCanvas = (TCanvas*)infileRECO.Get("MatchedLooseIsoRECOdR");
  TCanvas *MatchedMedIsoCJdRCanvas = (TCanvas*)infileCJ.Get("MatchedMedIsoCJdR");
  TCanvas *MatchedMedIsoRECOdRCanvas = (TCanvas*)infileRECO.Get("MatchedMedIsoRECOdR");
  TCanvas *MatchedTightIsoCJdRCanvas = (TCanvas*)infileCJ.Get("MatchedTightIsoCJdR");
  TCanvas *MatchedTightIsoRECOdRCanvas = (TCanvas*)infileRECO.Get("MatchedTightIsoRECOdR");
  TCanvas *MatchedDMFindCJdRCanvas = (TCanvas*)infileCJ.Get("MatchedDMFindCJdR");
  TCanvas *MatchedDMFindRECOdRCanvas = (TCanvas*)infileRECO.Get("MatchedDMFindRECOdR");
  TCanvas *MatchedCJdRCanvas = (TCanvas*)infileCJ.Get("MatchedCJdR");
  TCanvas *MatchedRECOdRCanvas = (TCanvas*)infileRECO.Get("MatchedRECOdR");

  TCanvas *MatchedLooseIsoCJPtGenCanvas = (TCanvas*)infileCJ.Get("MatchedLooseIsoCJPtGen");
  TCanvas *MatchedLooseIsoRECOPtGenCanvas = (TCanvas*)infileRECO.Get("MatchedLooseIsoRECOPtGen");
  TCanvas *MatchedMedIsoCJPtGenCanvas = (TCanvas*)infileCJ.Get("MatchedMedIsoCJPtGen");
  TCanvas *MatchedMedIsoRECOPtGenCanvas = (TCanvas*)infileRECO.Get("MatchedMedIsoRECOPtGen");
  TCanvas *MatchedTightIsoCJPtGenCanvas = (TCanvas*)infileCJ.Get("MatchedTightIsoCJPtGen");
  TCanvas *MatchedTightIsoRECOPtGenCanvas = (TCanvas*)infileRECO.Get("MatchedTightIsoRECOPtGen");
  TCanvas *MatchedDMFindCJPtGenCanvas = (TCanvas*)infileCJ.Get("MatchedDMFindCJPtGen");
  TCanvas *MatchedDMFindRECOPtGenCanvas = (TCanvas*)infileRECO.Get("MatchedDMFindRECOPtGen");
  TCanvas *MatchedCJPtGenCanvas = (TCanvas*)infileCJ.Get("MatchedCJPtGen");
  TCanvas *MatchedRECOPtGenCanvas = (TCanvas*)infileRECO.Get("MatchedRECOPtGen");

cout << "Got Canvases" << endl;

  TH1F* MatchedLooseIsoCJPt_ = (TH1F*)MatchedLooseIsoCJPtCanvas->GetPrimitive("MatchedLooseIsoCJPt");
  TH1F* MatchedLooseIsoRECOPt_ = (TH1F*)MatchedLooseIsoRECOPtCanvas->GetPrimitive("MatchedLooseIsoRECOPt");
  TH1F* MatchedMedIsoCJPt_ = (TH1F*)MatchedMedIsoCJPtCanvas->GetPrimitive("MatchedMedIsoCJPt");
  TH1F* MatchedMedIsoRECOPt_ = (TH1F*)MatchedMedIsoRECOPtCanvas->GetPrimitive("MatchedMedIsoRECOPt");
  TH1F* MatchedTightIsoCJPt_ = (TH1F*)MatchedTightIsoCJPtCanvas->GetPrimitive("MatchedTightIsoCJPt");
  TH1F* MatchedTightIsoRECOPt_ = (TH1F*)MatchedTightIsoRECOPtCanvas->GetPrimitive("MatchedTightIsoRECOPt");
  TH1F* MatchedDMFindCJPt_ = (TH1F*)MatchedDMFindCJPtCanvas->GetPrimitive("MatchedDMFindCJPt");
  TH1F* MatchedDMFindRECOPt_ = (TH1F*)MatchedDMFindRECOPtCanvas->GetPrimitive("MatchedDMFindRECOPt");
  TH1F* MatchedCJPt_ = (TH1F*)MatchedCJPtCanvas->GetPrimitive("MatchedCJPt");
  TH1F* MatchedRECOPt_ = (TH1F*)MatchedRECOPtCanvas->GetPrimitive("MatchedRECOPt");

  TH1F* MatchedLooseIsoCJdR_ = (TH1F*)MatchedLooseIsoCJdRCanvas->GetPrimitive("MatchedLooseIsoCJdR");
  TH1F* MatchedLooseIsoRECOdR_ = (TH1F*)MatchedLooseIsoRECOdRCanvas->GetPrimitive("MatchedLooseIsoRECOdR");
  TH1F* MatchedMedIsoCJdR_ = (TH1F*)MatchedMedIsoCJdRCanvas->GetPrimitive("MatchedMedIsoCJdR");
  TH1F* MatchedMedIsoRECOdR_ = (TH1F*)MatchedMedIsoRECOdRCanvas->GetPrimitive("MatchedMedIsoRECOdR");
  TH1F* MatchedTightIsoCJdR_ = (TH1F*)MatchedTightIsoCJdRCanvas->GetPrimitive("MatchedTightIsoCJdR");
  TH1F* MatchedTightIsoRECOdR_ = (TH1F*)MatchedTightIsoRECOdRCanvas->GetPrimitive("MatchedTightIsoRECOdR");
  TH1F* MatchedDMFindCJdR_ = (TH1F*)MatchedDMFindCJdRCanvas->GetPrimitive("MatchedDMFindCJdR");
  TH1F* MatchedDMFindRECOdR_ = (TH1F*)MatchedDMFindRECOdRCanvas->GetPrimitive("MatchedDMFindRECOdR");
  TH1F* MatchedCJdR_ = (TH1F*)MatchedCJdRCanvas->GetPrimitive("MatchedCJdR");
  TH1F* MatchedRECOdR_ = (TH1F*)MatchedRECOdRCanvas->GetPrimitive("MatchedRECOdR");

  TH1F* MatchedLooseIsoCJPtGen_ = (TH1F*)MatchedLooseIsoCJPtGenCanvas->GetPrimitive("MatchedLooseIsoCJPtGen");
  TH1F* MatchedLooseIsoRECOPtGen_ = (TH1F*)MatchedLooseIsoRECOPtGenCanvas->GetPrimitive("MatchedLooseIsoRECOPtGen");
  TH1F* MatchedMedIsoCJPtGen_ = (TH1F*)MatchedMedIsoCJPtGenCanvas->GetPrimitive("MatchedMedIsoCJPtGen");
  TH1F* MatchedMedIsoRECOPtGen_ = (TH1F*)MatchedMedIsoRECOPtGenCanvas->GetPrimitive("MatchedMedIsoRECOPtGen");
  TH1F* MatchedTightIsoCJPtGen_ = (TH1F*)MatchedTightIsoCJPtGenCanvas->GetPrimitive("MatchedTightIsoCJPtGen");
  TH1F* MatchedTightIsoRECOPtGen_ = (TH1F*)MatchedTightIsoRECOPtGenCanvas->GetPrimitive("MatchedTightIsoRECOPtGen");
  TH1F* MatchedDMFindCJPtGen_ = (TH1F*)MatchedDMFindCJPtGenCanvas->GetPrimitive("MatchedDMFindCJPtGen");
  TH1F* MatchedDMFindRECOPtGen_ = (TH1F*)MatchedDMFindRECOPtGenCanvas->GetPrimitive("MatchedDMFindRECOPtGen");
  TH1F* MatchedCJPtGen_ = (TH1F*)MatchedCJPtGenCanvas->GetPrimitive("MatchedCJPtGen");
  TH1F* MatchedRECOPtGen_ = (TH1F*)MatchedRECOPtGenCanvas->GetPrimitive("MatchedRECOPtGen");

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


  FinalEffLooseIsoCJPt_->SetLineColor(kRed);
  FinalEffLooseIsoRECOPt_->SetLineColor(kBlack);
  FinalEffMedIsoCJPt_->SetLineColor(kRed);
  FinalEffMedIsoRECOPt_->SetLineColor(kBlack);
  FinalEffTightIsoCJPt_->SetLineColor(kRed);
  FinalEffTightIsoRECOPt_->SetLineColor(kBlack);
  FinalEffDMFindCJPt_->SetLineColor(kRed);
  FinalEffDMFindRECOPt_->SetLineColor(kBlack);

  FinalEffLooseIsoCJdR_->SetLineColor(kRed);
  FinalEffLooseIsoRECOdR_->SetLineColor(kBlack);
  FinalEffMedIsoCJdR_->SetLineColor(kRed);
  FinalEffMedIsoRECOdR_->SetLineColor(kBlack);
  FinalEffTightIsoCJdR_->SetLineColor(kRed);
  FinalEffTightIsoRECOdR_->SetLineColor(kBlack);
  FinalEffDMFindCJdR_->SetLineColor(kRed);
  FinalEffDMFindRECOdR_->SetLineColor(kBlack);

  FinalEffLooseIsoCJPtGen_->SetLineColor(kRed);
  FinalEffLooseIsoRECOPtGen_->SetLineColor(kBlack);
  FinalEffMedIsoCJPtGen_->SetLineColor(kRed);
  FinalEffMedIsoRECOPtGen_->SetLineColor(kBlack);
  FinalEffTightIsoCJPtGen_->SetLineColor(kRed);
  FinalEffTightIsoRECOPtGen_->SetLineColor(kBlack);
  FinalEffDMFindCJPtGen_->SetLineColor(kRed);
  FinalEffDMFindRECOPtGen_->SetLineColor(kBlack);



  FinalEffLooseIsoCJPt_->SetMarkerColor(kRed);
  FinalEffLooseIsoRECOPt_->SetMarkerColor(kBlack);
  FinalEffMedIsoCJPt_->SetMarkerColor(kRed);
  FinalEffMedIsoRECOPt_->SetMarkerColor(kBlack);
  FinalEffTightIsoCJPt_->SetMarkerColor(kRed);
  FinalEffTightIsoRECOPt_->SetMarkerColor(kBlack);
  FinalEffDMFindCJPt_->SetMarkerColor(kRed);
  FinalEffDMFindRECOPt_->SetMarkerColor(kBlack);
  
  FinalEffLooseIsoCJdR_->SetMarkerColor(kRed);
  FinalEffLooseIsoRECOdR_->SetMarkerColor(kBlack);
  FinalEffMedIsoCJdR_->SetMarkerColor(kRed);
  FinalEffMedIsoRECOdR_->SetMarkerColor(kBlack);
  FinalEffTightIsoCJdR_->SetMarkerColor(kRed);
  FinalEffTightIsoRECOdR_->SetMarkerColor(kBlack);
  FinalEffDMFindCJdR_->SetMarkerColor(kRed);
  FinalEffDMFindRECOdR_->SetMarkerColor(kBlack);

  FinalEffLooseIsoCJPtGen_->SetMarkerColor(kRed);
  FinalEffLooseIsoRECOPtGen_->SetMarkerColor(kBlack);
  FinalEffMedIsoCJPtGen_->SetMarkerColor(kRed);
  FinalEffMedIsoRECOPtGen_->SetMarkerColor(kBlack);
  FinalEffTightIsoCJPtGen_->SetMarkerColor(kRed);
  FinalEffTightIsoRECOPtGen_->SetMarkerColor(kBlack);
  FinalEffDMFindCJPtGen_->SetMarkerColor(kRed);
  FinalEffDMFindRECOPtGen_->SetMarkerColor(kBlack);

cout << "Attributes set." << endl;  

  MatchedLooseIsoCJPt_->SetXTitle("P_{T}");
  MatchedMedIsoCJPt_->SetXTitle("P_{T}");
  MatchedTightIsoCJPt_->SetXTitle("P_{T}");
  MatchedDMFindCJPt_->SetXTitle("P_{T}");

  MatchedLooseIsoCJPt_->SetYTitle("#epsilon(Loose Iso + DMFinding + GenMatch / GenMatch)");
  MatchedMedIsoCJPt_->SetYTitle("#epsilon(Medium Iso + DMFinding + GenMatch / GenMatch)");
  MatchedTightIsoCJPt_->SetYTitle("#epsilon(Tight Iso + DMFinding + GenMatch / GenMatch)");
  MatchedDMFindCJPt_->SetYTitle("#epsilon(DMFinding + GenMatch / GenMatch)");

  MatchedLooseIsoCJdR_->SetXTitle("#DeltaR");
  MatchedMedIsoCJdR_->SetXTitle("#DeltaR");
  MatchedTightIsoCJdR_->SetXTitle("#DeltaR");
  MatchedDMFindCJdR_->SetXTitle("#DeltaR");

  MatchedLooseIsoCJdR_->SetYTitle("#epsilon(Loose Iso + DMFinding + GenMatch / GenMatch)");
  MatchedMedIsoCJdR_->SetYTitle("#epsilon(Medium Iso + DMFinding + GenMatch / GenMatch)");
  MatchedTightIsoCJdR_->SetYTitle("#epsilon(Tight Iso + DMFinding + GenMatch / GenMatch)");
  MatchedDMFindCJdR_->SetYTitle("#epsilon(DMFinding + GenMatch / GenMatch)");

  MatchedLooseIsoCJPtGen_->SetXTitle("Visible Gen P_{T}");
  MatchedMedIsoCJPtGen_->SetXTitle("Visible Gen P_{T}");
  MatchedTightIsoCJPtGen_->SetXTitle("Visible Gen P_{T}");
  MatchedDMFindCJPtGen_->SetXTitle("Visible Gen P_{T}");

  MatchedLooseIsoCJPtGen_->SetYTitle("#epsilon(Loose Iso + DMFinding + GenMatch / GenMatch)");
  MatchedMedIsoCJPtGen_->SetYTitle("#epsilon(Medium Iso + DMFinding + GenMatch / GenMatch)");
  MatchedTightIsoCJPtGen_->SetYTitle("#epsilon(Tight Iso + DMFinding + GenMatch / GenMatch)");
  MatchedDMFindCJPtGen_->SetYTitle("#epsilon(DMFinding + GenMatch / GenMatch)");

  leg = new TLegend(0.1,0.7,0.25,0.9);
  leg->AddEntry(FinalEffDMFindRECOPtGen_, "No CleanJets","L");
  leg->AddEntry(FinalEffDMFindCJPtGen_, "CleanJets","L");

  LoosePtCanvas.cd();
  FinalEffLooseIsoCJPt_->Draw();
  FinalEffLooseIsoRECOPt_->Draw("SAME");
  leg->Draw();
  MedPtCanvas.cd();
  FinalEffMedIsoCJPt_->Draw();  
  FinalEffMedIsoRECOPt_->Draw("SAME");
  leg->Draw();
  TightPtCanvas.cd();
  FinalEffTightIsoCJPt_->Draw();  
  FinalEffTightIsoRECOPt_->Draw("SAME");
  leg->Draw();
  DMPtCanvas.cd();
  FinalEffDMFindCJPt_->Draw();  
  FinalEffDMFindRECOPt_->Draw("SAME");
  leg->Draw();
  
  LoosedRCanvas.cd();
  FinalEffLooseIsoCJdR_->Draw();
  FinalEffLooseIsoRECOdR_->Draw("SAME");
  leg->Draw();
  MeddRCanvas.cd();
  FinalEffMedIsoCJdR_->Draw();
  FinalEffMedIsoRECOdR_->Draw("SAME");
  leg->Draw();
  TightdRCanvas.cd();
  FinalEffTightIsoCJdR_->Draw();
  FinalEffTightIsoRECOdR_->Draw("SAME");
  leg->Draw();
  DMdRCanvas.cd();
  FinalEffDMFindCJdR_->Draw();
  FinalEffDMFindRECOdR_->Draw("SAME");
  leg->Draw();
  
  LoosePtGenCanvas.cd();
  FinalEffLooseIsoCJPtGen_->Draw();
  FinalEffLooseIsoRECOPtGen_->Draw("SAME");
  leg->Draw();
  MedPtGenCanvas.cd();
  FinalEffMedIsoCJPtGen_->Draw();
  FinalEffMedIsoRECOPtGen_->Draw("SAME");
  leg->Draw();
  TightPtGenCanvas.cd();
  FinalEffTightIsoCJPtGen_->Draw();
  FinalEffTightIsoRECOPtGen_->Draw("SAME");
  leg->Draw();
  DMPtGenCanvas.cd();
  FinalEffDMFindCJPtGen_->Draw();
  FinalEffDMFindRECOPtGen_->Draw("SAME");
  leg->Draw();

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

  outFile->Write();
  outFile->Close();

cout << "end" << endl;

}//rootMacro_BBA_combine

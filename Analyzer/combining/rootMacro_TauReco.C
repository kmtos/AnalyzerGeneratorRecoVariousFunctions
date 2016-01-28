#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

void rootMacro_TauReco()
{
  gStyle->SetOptStat(kFALSE);

  TFile infile("/afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/BBA/Analyzer/BSUB/BBATauReco_a30_muchMore_DEC_14_v1/BBATauReco_a30_muchMore_DEC_14_v1.root");
  TFile *outFile = new TFile("combHist_BBA_a30_TauReco.root", "RECREATE");

cout << "Files Created" << endl;

//  TCanvas *CJetPt_ = (TCanvas*)infile.Get("JetPt");
  TCanvas *FinalEffLooseIsoCJPtCanvas = (TCanvas*)infile.Get("FinalEffLooseIsoCJPt");
  TCanvas *FinalEffLooseIsoRECOPtCanvas = (TCanvas*)infile.Get("FinalEffLooseIsoRECOPt");
  TCanvas *FinalEffMedIsoCJPtCanvas = (TCanvas*)infile.Get("FinalEffMedIsoCJPt");
  TCanvas *FinalEffMedIsoRECOPtCanvas = (TCanvas*)infile.Get("FinalEffMedIsoRECOPt");
  TCanvas *FinalEffTightIsoCJPtCanvas = (TCanvas*)infile.Get("FinalEffTightIsoCJPt");
  TCanvas *FinalEffTightIsoRECOPtCanvas = (TCanvas*)infile.Get("FinalEffTightIsoRECOPt");
  TCanvas *FinalEffDMFindCJPtCanvas = (TCanvas*)infile.Get("FinalEffDMFindCJPt");
  TCanvas *FinalEffDMFindRECOPtCanvas = (TCanvas*)infile.Get("FinalEffDMFindRECOPt");

  TCanvas *FinalEffLooseIsoCJdRCanvas = (TCanvas*)infile.Get("FinalEffLooseIsoCJdR");
  TCanvas *FinalEffLooseIsoRECOdRCanvas = (TCanvas*)infile.Get("FinalEffLooseIsoRECOdR");
  TCanvas *FinalEffMedIsoCJdRCanvas = (TCanvas*)infile.Get("FinalEffMedIsoCJdR");
  TCanvas *FinalEffMedIsoRECOdRCanvas = (TCanvas*)infile.Get("FinalEffMedIsoRECOdR");
  TCanvas *FinalEffTightIsoCJdRCanvas = (TCanvas*)infile.Get("FinalEffTightIsoCJdR");
  TCanvas *FinalEffTightIsoRECOdRCanvas = (TCanvas*)infile.Get("FinalEffTightIsoRECOdR");
  TCanvas *FinalEffDMFindCJdRCanvas = (TCanvas*)infile.Get("FinalEffDMFindCJdR");
  TCanvas *FinalEffDMFindRECOdRCanvas = (TCanvas*)infile.Get("FinalEffDMFindRECOdR");

  TCanvas *FinalEffLooseIsoCJPtGenCanvas = (TCanvas*)infile.Get("FinalEffLooseIsoCJPtGen");
  TCanvas *FinalEffLooseIsoRECOPtGenCanvas = (TCanvas*)infile.Get("FinalEffLooseIsoRECOPtGen");
  TCanvas *FinalEffMedIsoCJPtGenCanvas = (TCanvas*)infile.Get("FinalEffMedIsoCJPtGen");
  TCanvas *FinalEffMedIsoRECOPtGenCanvas = (TCanvas*)infile.Get("FinalEffMedIsoRECOPtGen");
  TCanvas *FinalEffTightIsoCJPtGenCanvas = (TCanvas*)infile.Get("FinalEffTightIsoCJPtGen");
  TCanvas *FinalEffTightIsoRECOPtGenCanvas = (TCanvas*)infile.Get("FinalEffTightIsoRECOPtGen");
  TCanvas *FinalEffDMFindCJPtGenCanvas = (TCanvas*)infile.Get("FinalEffDMFindCJPtGen");
  TCanvas *FinalEffDMFindRECOPtGenCanvas = (TCanvas*)infile.Get("FinalEffDMFindRECOPtGen");

cout << "Got Canvases" << endl;

  TGraphAsymmErrors* FinalEffLooseIsoCJPt_ = (TGraphAsymmErrors*)FinalEffLooseIsoCJPtCanvas->GetPrimitive("FinalEffLooseIsoCJPt");
  TGraphAsymmErrors* FinalEffLooseIsoRECOPt_ = (TGraphAsymmErrors*)FinalEffLooseIsoRECOPtCanvas->GetPrimitive("FinalEffLooseIsoRECOPt");
  TGraphAsymmErrors* FinalEffMedIsoCJPt_ = (TGraphAsymmErrors*)FinalEffMedIsoCJPtCanvas->GetPrimitive("FinalEffMedIsoCJPt");
  TGraphAsymmErrors* FinalEffMedIsoRECOPt_ = (TGraphAsymmErrors*)FinalEffMedIsoRECOPtCanvas->GetPrimitive("FinalEffMedIsoRECOPt");
  TGraphAsymmErrors* FinalEffTightIsoCJPt_ = (TGraphAsymmErrors*)FinalEffTightIsoCJPtCanvas->GetPrimitive("FinalEffTightIsoCJPt");
  TGraphAsymmErrors* FinalEffTightIsoRECOPt_ = (TGraphAsymmErrors*)FinalEffTightIsoRECOPtCanvas->GetPrimitive("FinalEffTightIsoRECOPt");
  TGraphAsymmErrors* FinalEffDMFindCJPt_ = (TGraphAsymmErrors*)FinalEffDMFindCJPtCanvas->GetPrimitive("FinalEffDMFindCJPt");
  TGraphAsymmErrors* FinalEffDMFindRECOPt_ = (TGraphAsymmErrors*)FinalEffDMFindRECOPtCanvas->GetPrimitive("FinalEffDMFindRECOPt");

  TGraphAsymmErrors* FinalEffLooseIsoCJdR_ = (TGraphAsymmErrors*)FinalEffLooseIsoCJdRCanvas->GetPrimitive("FinalEffLooseIsoCJdR");
  TGraphAsymmErrors* FinalEffLooseIsoRECOdR_ = (TGraphAsymmErrors*)FinalEffLooseIsoRECOdRCanvas->GetPrimitive("FinalEffLooseIsoRECOdR");
  TGraphAsymmErrors* FinalEffMedIsoCJdR_ = (TGraphAsymmErrors*)FinalEffMedIsoCJdRCanvas->GetPrimitive("FinalEffMedIsoCJdR");
  TGraphAsymmErrors* FinalEffMedIsoRECOdR_ = (TGraphAsymmErrors*)FinalEffMedIsoRECOdRCanvas->GetPrimitive("FinalEffMedIsoRECOdR");
  TGraphAsymmErrors* FinalEffTightIsoCJdR_ = (TGraphAsymmErrors*)FinalEffTightIsoCJdRCanvas->GetPrimitive("FinalEffTightIsoCJdR");
  TGraphAsymmErrors* FinalEffTightIsoRECOdR_ = (TGraphAsymmErrors*)FinalEffTightIsoRECOdRCanvas->GetPrimitive("FinalEffTightIsoRECOdR");
  TGraphAsymmErrors* FinalEffDMFindCJdR_ = (TGraphAsymmErrors*)FinalEffDMFindCJdRCanvas->GetPrimitive("FinalEffDMFindCJdR");
  TGraphAsymmErrors* FinalEffDMFindRECOdR_ = (TGraphAsymmErrors*)FinalEffDMFindRECOdRCanvas->GetPrimitive("FinalEffDMFindRECOdR");

  TGraphAsymmErrors* FinalEffLooseIsoCJPtGen_ = (TGraphAsymmErrors*)FinalEffLooseIsoCJPtGenCanvas->GetPrimitive("FinalEffLooseIsoCJPtGen");
  TGraphAsymmErrors* FinalEffLooseIsoRECOPtGen_ = (TGraphAsymmErrors*)FinalEffLooseIsoRECOPtGenCanvas->GetPrimitive("FinalEffLooseIsoRECOPtGen");
  TGraphAsymmErrors* FinalEffMedIsoCJPtGen_ = (TGraphAsymmErrors*)FinalEffMedIsoCJPtGenCanvas->GetPrimitive("FinalEffMedIsoCJPtGen");
  TGraphAsymmErrors* FinalEffMedIsoRECOPtGen_ = (TGraphAsymmErrors*)FinalEffMedIsoRECOPtGenCanvas->GetPrimitive("FinalEffMedIsoRECOPtGen");
  TGraphAsymmErrors* FinalEffTightIsoCJPtGen_ = (TGraphAsymmErrors*)FinalEffTightIsoCJPtGenCanvas->GetPrimitive("FinalEffTightIsoCJPtGen");
  TGraphAsymmErrors* FinalEffTightIsoRECOPtGen_ = (TGraphAsymmErrors*)FinalEffTightIsoRECOPtGenCanvas->GetPrimitive("FinalEffTightIsoRECOPtGen");
  TGraphAsymmErrors* FinalEffDMFindCJPtGen_ = (TGraphAsymmErrors*)FinalEffDMFindCJPtGenCanvas->GetPrimitive("FinalEffDMFindCJPtGen");
  TGraphAsymmErrors* FinalEffDMFindRECOPtGen_ = (TGraphAsymmErrors*)FinalEffDMFindRECOPtGenCanvas->GetPrimitive("FinalEffDMFindRECOPtGen");

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

  FinalEffLooseIsoCJPt_->SetLineColor(kBlack);
  FinalEffLooseIsoRECOPt_->SetLineColor(kRed);
  FinalEffMedIsoCJPt_->SetLineColor(kBlack);
  FinalEffMedIsoRECOPt_->SetLineColor(kRed);
  FinalEffTightIsoCJPt_->SetLineColor(kBlack);
  FinalEffTightIsoRECOPt_->SetLineColor(kRed);
  FinalEffDMFindCJPt_->SetLineColor(kBlack);
  FinalEffDMFindRECOPt_->SetLineColor(kRed);
  FinalEffLooseIsoCJdR_->SetLineColor(kBlack);
  FinalEffLooseIsoRECOdR_->SetLineColor(kRed);
  FinalEffMedIsoCJdR_->SetLineColor(kBlack);
  FinalEffMedIsoRECOdR_->SetLineColor(kRed);
  FinalEffTightIsoCJdR_->SetLineColor(kBlack);
  FinalEffTightIsoRECOdR_->SetLineColor(kRed);
  FinalEffDMFindCJdR_->SetLineColor(kBlack);
  FinalEffDMFindRECOdR_->SetLineColor(kRed);
  FinalEffLooseIsoCJPtGen_->SetLineColor(kBlack);
  FinalEffLooseIsoRECOPtGen_->SetLineColor(kRed);
  FinalEffMedIsoCJPtGen_->SetLineColor(kBlack);
  FinalEffMedIsoRECOPtGen_->SetLineColor(kRed);
  FinalEffTightIsoCJPtGen_->SetLineColor(kBlack);
  FinalEffTightIsoRECOPtGen_->SetLineColor(kRed);
  FinalEffDMFindCJPtGen_->SetLineColor(kBlack);
  FinalEffDMFindRECOPtGen_->SetLineColor(kRed);

  FinalEffLooseIsoCJPt_->SetMarkerColor(kBlack);
  FinalEffLooseIsoRECOPt_->SetMarkerColor(kRed);
  FinalEffMedIsoCJPt_->SetMarkerColor(kBlack);
  FinalEffMedIsoRECOPt_->SetMarkerColor(kRed);
  FinalEffTightIsoCJPt_->SetMarkerColor(kBlack);
  FinalEffTightIsoRECOPt_->SetMarkerColor(kRed);
  FinalEffDMFindCJPt_->SetMarkerColor(kBlack);
  FinalEffDMFindRECOPt_->SetMarkerColor(kRed);
  FinalEffLooseIsoCJdR_->SetMarkerColor(kBlack);
  FinalEffLooseIsoRECOdR_->SetMarkerColor(kRed);
  FinalEffMedIsoCJdR_->SetMarkerColor(kBlack);
  FinalEffMedIsoRECOdR_->SetMarkerColor(kRed);
  FinalEffTightIsoCJdR_->SetMarkerColor(kBlack);
  FinalEffTightIsoRECOdR_->SetMarkerColor(kRed);
  FinalEffDMFindCJdR_->SetMarkerColor(kBlack);
  FinalEffDMFindRECOdR_->SetMarkerColor(kRed);
  FinalEffLooseIsoCJPtGen_->SetMarkerColor(kBlack);
  FinalEffLooseIsoRECOPtGen_->SetMarkerColor(kRed);
  FinalEffMedIsoCJPtGen_->SetMarkerColor(kBlack);
  FinalEffMedIsoRECOPtGen_->SetMarkerColor(kRed);
  FinalEffTightIsoCJPtGen_->SetMarkerColor(kBlack);
  FinalEffTightIsoRECOPtGen_->SetMarkerColor(kRed);
  FinalEffDMFindCJPtGen_->SetMarkerColor(kBlack);
  FinalEffDMFindRECOPtGen_->SetMarkerColor(kRed);

cout << "Attributes set." << endl;  

  FinalEffLooseIsoCJPt_->SetXTitle("P_{T}(Loose Iso + DMFinding + GenMatch) / P_{T}(GenMatch)");
  FinalEffMedIsoCJPt_->SetXTitle("P_{T}(Medium Iso + DMFinding + GenMatch) / P_{T}(GenMatch)");
  FinalEffTightIsoCJPt_->SetXTitle("P_{T}(Tight Iso + DMFinding + GenMatch) / P_{T}(GenMatch)");
  FinalEffDMFindCJPt_->SetXTitle("P_{T}(DMFinding + GenMatch) / P_{T}(GenMatch)");

  FinalEffLooseIsoCJdR_->SetXTitle("#DeltaR(Loose Iso + DMFinding + GenMatch) / #DeltaR(GenMatch)");
  FinalEffMedIsoCJdR_->SetXTitle("#DeltaR(Medium Iso + DMFinding + GenMatch) / #DeltaR(GenMatch)");
  FinalEffTightIsoCJdR_->SetXTitle("#DeltaR(Tight Iso + DMFinding + GenMatch) / #DeltaR(GenMatch)");
  FinalEffDMFindIsoCJdR_->SetXTitle("#DeltaR(DMFinding + GenMatch) / #DeltaR(GenMatch)");

  FinalEffLooseIsoCJPtGen_->SetXTitle("Visible Gen P_{T}(Loose Iso + DMFinding + GenMatch) / Visible Gen P_{T}(GenMatch)");
  FinalEffMedIsoCJPtGen_->SetXTitle("Visible Gen P_{T}(Medium Iso + DMFinding + GenMatch) / Visible Gen P_{T}(GenMatch)");
  FinalEffTightIsoCJPtGen_->SetXTitle("Visible Gen P_{T}(Tight Iso + DMFinding + GenMatch) / Visible Gen P_{T}(GenMatch)");
  FinalEffDMFindCJPtGen_->SetXTitle("Visible Gen P_{T}(DMFinding + GenMatch) / Visible Gen P_{T}(GenMatch)");

  leg = new TLegend(0.1,0.7,0.25,0.9);
  leg->AddEntry(FinalEffDMFindRECOPtGen_, "No CleanJets","L");
  leg->AddEntry(FinalEffDMFindCJPtGen_, "CleanJets","L");

  LoosePtCanvas.cd()
  FinalEffLooseIsoCJPt_->Draw();
  FinalEffLooseIsoRECOPt_->Draw("SAME");
  leg->Draw();
  MedPtCanvas.cd()
  FinalEffMediumIsoCJPt_->Draw();  
  FinalEffMediumIsoRECOPt_->Draw("SAME");
  leg->Draw();
  TightPtCanvas.cd()
  FinalEffTightIsoCJPt_->Draw();  
  FinalEffTightIsoRECOPt_->Draw("SAME");
  leg->Draw();
  DMPtCanvas.cd()
  FinalEffDMFindCJPt_->Draw();  
  FinalEffDMFindRECOPt_->Draw("SAME");
  leg->Draw();
  
  LoosedRCanvas.cd()
  FinalEffLooseIsoCJdR_->Draw();
  FinalEffLooseIsoRECOdR_->Draw("SAME");
  leg->Draw();
  MeddRCanvas.cd()
  FinalEffMediumIsoCJdR_->Draw();
  FinalEffMediumIsoRECOdR_->Draw("SAME");
  leg->Draw();
  TightdRCanvas.cd()
  FinalEffTightIsoCJdR_->Draw();
  FinalEffTightIsoRECOdR_->Draw("SAME");
  leg->Draw();
  DMdRCanvas.cd()
  FinalEffDMFindCJdR_->Draw();
  FinalEffDMFindRECOdR_->Draw("SAME");
  leg->Draw();
  
  LoosePtGenCanvas.cd()
  FinalEffLooseIsoCJPtGen_->Draw();
  FinalEffLooseIsoRECOPtGen_->Draw("SAME");
  leg->Draw();
  MedPtGenCanvas.cd()
  FinalEffMediumIsoCJPtGen_->Draw();
  FinalEffMediumIsoRECOPtGen_->Draw("SAME");
  leg->Draw();
  TightPtGenCanvas.cd()
  FinalEffTightIsoCJPtGen_->Draw();
  FinalEffTightIsoRECOPtGen_->Draw("SAME");
  leg->Draw();
  DMPtGenCanvas.cd()
  FinalEffDMFindCJPtGen_->Draw();
  FinalEffDMFindRECOPtGen_->Draw("SAME");
  leg->Draw();

cout << "Histograms Drawn" << endl;

  outFile->cd();
  JetPtCSVALLCanvas.Write();
  JetPtDiffCSVALLCanvas.Write();
  JetDRCSVALLCanvas.Write();
  outFile->Write();
  outFile->Close();

cout << "end" << endl;

}//rootMacro_BBA_combine

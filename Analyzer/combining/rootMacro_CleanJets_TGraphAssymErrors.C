#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include "TCanvas.h"

void rootMacro_CleanJets()
{
  gStyle->SetOptStat(kFALSE);

  TFile infileCJ("/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_4_12_patch4/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/ggH125a9_GenTauDecayID_IndividualRCJ_FEB25/ggH125a9_GenTauDecayID_IndividualRCJ_FEB25_Plots.root");
  TFile infileRECO("/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_4_12_patch4/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/ggH125a9_GenTauDecayID_IndividualRECO_FEB25/ggH125a9_GenTauDecayID_IndividualRECO_FEB25_Plots.root");
  TFile *outFile = new TFile("combHist_CleanJets_h125a9.root", "RECREATE");

cout << "Files Created" << endl;

  TCanvas *FinalEffLooseIsoCJPtCanvas = (TCanvas*)infileCJ.Get("FinalEffLooseIsoCJPt");
  TCanvas *FinalEffLooseIsoRECOPtCanvas = (TCanvas*)infileRECO.Get("FinalEffLooseIsoRECOPt");
  TCanvas *FinalEffMedIsoCJPtCanvas = (TCanvas*)infileCJ.Get("FinalEffMedIsoCJPt");
  TCanvas *FinalEffMedIsoRECOPtCanvas = (TCanvas*)infileRECO.Get("FinalEffMedIsoRECOPt");
  TCanvas *FinalEffTightIsoCJPtCanvas = (TCanvas*)infileCJ.Get("FinalEffTightIsoCJPt");
  TCanvas *FinalEffTightIsoRECOPtCanvas = (TCanvas*)infileRECO.Get("FinalEffTightIsoRECOPt");
  TCanvas *FinalEffDMFindCJPtCanvas = (TCanvas*)infileCJ.Get("FinalEffDMFindCJPt");
  TCanvas *FinalEffDMFindRECOPtCanvas = (TCanvas*)infileRECO.Get("FinalEffDMFindRECOPt");


  TCanvas *FinalEffLooseIsoCJdRCanvas = (TCanvas*)infileCJ.Get("FinalEffLooseIsoCJdR");
  TCanvas *FinalEffLooseIsoRECOdRCanvas = (TCanvas*)infileRECO.Get("FinalEffLooseIsoRECOdR");
  TCanvas *FinalEffMedIsoCJdRCanvas = (TCanvas*)infileCJ.Get("FinalEffMedIsoCJdR");
  TCanvas *FinalEffMedIsoRECOdRCanvas = (TCanvas*)infileRECO.Get("FinalEffMedIsoRECOdR");
  TCanvas *FinalEffTightIsoCJdRCanvas = (TCanvas*)infileCJ.Get("FinalEffTightIsoCJdR");
  TCanvas *FinalEffTightIsoRECOdRCanvas = (TCanvas*)infileRECO.Get("FinalEffTightIsoRECOdR");
  TCanvas *FinalEffDMFindCJdRCanvas = (TCanvas*)infileCJ.Get("FinalEffDMFindCJdR");
  TCanvas *FinalEffDMFindRECOdRCanvas = (TCanvas*)infileRECO.Get("FinalEffDMFindRECOdR");

  TCanvas *FinalEffLooseIsoCJPtGenCanvas = (TCanvas*)infileCJ.Get("FinalEffLooseIsoCJPtGen");
  TCanvas *FinalEffLooseIsoRECOPtGenCanvas = (TCanvas*)infileRECO.Get("FinalEffLooseIsoRECOPtGen");
  TCanvas *FinalEffMedIsoCJPtGenCanvas = (TCanvas*)infileCJ.Get("FinalEffMedIsoCJPtGen");
  TCanvas *FinalEffMedIsoRECOPtGenCanvas = (TCanvas*)infileRECO.Get("FinalEffMedIsoRECOPtGen");
  TCanvas *FinalEffTightIsoCJPtGenCanvas = (TCanvas*)infileCJ.Get("FinalEffTightIsoCJPtGen");
  TCanvas *FinalEffTightIsoRECOPtGenCanvas = (TCanvas*)infileRECO.Get("FinalEffTightIsoRECOPtGen");
  TCanvas *FinalEffDMFindCJPtGenCanvas = (TCanvas*)infileCJ.Get("FinalEffDMFindCJPtGen");
  TCanvas *FinalEffDMFindRECOPtGenCanvas = (TCanvas*)infileRECO.Get("FinalEffDMFindRECOPtGen");

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

  FinalEffLooseIsoCJPt_->SetXTitle("P_{T}");
  FinalEffMedIsoCJPt_->SetXTitle("P_{T}");
  FinalEffTightIsoCJPt_->SetXTitle("P_{T}");
  FinalEffDMFindCJPt_->SetXTitle("P_{T}");

  FinalEffLooseIsoCJPt_->SetYTitle("#epsilon(Loose Iso + DMFinding + GenMatch / GenMatch)");
  FinalEffMedIsoCJPt_->SetYTitle("#epsilon(Medium Iso + DMFinding + GenMatch / GenMatch)");
  FinalEffTightIsoCJPt_->SetYTitle("#epsilon(Tight Iso + DMFinding + GenMatch / GenMatch)");
  FinalEffDMFindCJPt_->SetYTitle("#epsilon(DMFinding + GenMatch / GenMatch)");

  FinalEffLooseIsoCJdR_->SetXTitle("#DeltaR");
  FinalEffMedIsoCJdR_->SetXTitle("#DeltaR");
  FinalEffTightIsoCJdR_->SetXTitle("#DeltaR");
  FinalEffDMFindIsoCJdR_->SetXTitle("#DeltaR");

  FinalEffLooseIsoCJdR_->SetYTitle("#epsilon(Loose Iso + DMFinding + GenMatch / GenMatch)");
  FinalEffMedIsoCJdR_->SetYTitle("#epsilon(Medium Iso + DMFinding + GenMatch / GenMatch)");
  FinalEffTightIsoCJdR_->SetYTitle("#epsilon(Tight Iso + DMFinding + GenMatch / GenMatch)");
  FinalEffDMFindIsoCJdR_->SetYTitle("#epsilon(DMFinding + GenMatch / GenMatch)");

  FinalEffLooseIsoCJPtGen_->SetXTitle("Visible Gen P_{T}");
  FinalEffMedIsoCJPtGen_->SetXTitle("Visible Gen P_{T}");
  FinalEffTightIsoCJPtGen_->SetXTitle("Visible Gen P_{T}");
  FinalEffDMFindCJPtGen_->SetXTitle("Visible Gen P_{T}");

  FinalEffLooseIsoCJPtGen_->SetYTitle("#epsilon(Loose Iso + DMFinding + GenMatch / GenMatch)");
  FinalEffMedIsoCJPtGen_->SetYTitle("#epsilon(Medium Iso + DMFinding + GenMatch / GenMatch)");
  FinalEffTightIsoCJPtGen_->SetYTitle("#epsilon(Tight Iso + DMFinding + GenMatch / GenMatch)");
  FinalEffDMFindCJPtGen_->SetYTitle("#epsilon(DMFinding + GenMatch / GenMatch)");

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

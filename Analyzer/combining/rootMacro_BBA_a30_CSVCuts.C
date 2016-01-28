#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

void rootMacro_BBA_a30_CSVCuts() 
{
  gStyle->SetOptStat(kFALSE);

  TFile infile("/afs/cern.ch/user/k/ktos/BBA/CMSSW_7_4_1_patch1/src/BBA/Analyzer/BSUB/BBA_a30_AUG_26_FIX_v2/BBA_a30_AUG_26_FIX_v2.root");
  TFile *outFile = new TFile("combHist_BBA_a30_CSVCuts.root", "RECREATE");

cout << "Files Created" << endl;

//  TCanvas *CJetPt_ = (TCanvas*)infile.Get("JetPt");
  TCanvas *CJetPtCSV1_ = (TCanvas*)infile.Get("JetPtCSV1");
  TCanvas *CJetPtCSV2_ = (TCanvas*)infile.Get("JetPtCSV2");
  TCanvas *CJetPtCSV3_ = (TCanvas*)infile.Get("JetPtCSV3");
  TCanvas *CJetPtCSV4_ = (TCanvas*)infile.Get("JetPtCSV4");
  TCanvas *CJetPtCSV5_ = (TCanvas*)infile.Get("JetPtCSV5");
  TCanvas *CJetPtCSV6_ = (TCanvas*)infile.Get("JetPtCSV6");

//  TCanvas *CJetPtDiff_ = (TCanvas*)infile.Get("JetPtDiff");
  TCanvas *CJetPtDiffCSV1_ = (TCanvas*)infile.Get("JetPtDiffCSV1");
  TCanvas *CJetPtDiffCSV2_ = (TCanvas*)infile.Get("JetPtDiffCSV2");
  TCanvas *CJetPtDiffCSV3_ = (TCanvas*)infile.Get("JetPtDiffCSV3");
  TCanvas *CJetPtDiffCSV4_ = (TCanvas*)infile.Get("JetPtDiffCSV4");
  TCanvas *CJetPtDiffCSV5_ = (TCanvas*)infile.Get("JetPtDiffCSV5");
  TCanvas *CJetPtDiffCSV6_ = (TCanvas*)infile.Get("JetPtDiffCSV6");

//  TCanvas *CJetDR_ = (TCanvas*)infile.Get("JetDR");
  TCanvas *CJetDRCSV1_ = (TCanvas*)infile.Get("JetDRCSV1");
  TCanvas *CJetDRCSV2_ = (TCanvas*)infile.Get("JetDRCSV2");
  TCanvas *CJetDRCSV3_ = (TCanvas*)infile.Get("JetDRCSV3");
  TCanvas *CJetDRCSV4_ = (TCanvas*)infile.Get("JetDRCSV4");
  TCanvas *CJetDRCSV5_ = (TCanvas*)infile.Get("JetDRCSV5");
  TCanvas *CJetDRCSV6_ = (TCanvas*)infile.Get("JetDRCSV6");

cout << "Got Canvases" << endl;

//  TH1F* JetPt_ = (TH1F*)CJetPt_->GetPrimitive("JetPt");
  TH1F* JetPtCSV1_ = (TH1F*)CJetPtCSV1_->GetPrimitive("JetPtCSV1");
  TH1F* JetPtCSV2_ = (TH1F*)CJetPtCSV2_->GetPrimitive("JetPtCSV2");
  TH1F* JetPtCSV3_ = (TH1F*)CJetPtCSV3_->GetPrimitive("JetPtCSV3");
  TH1F* JetPtCSV4_ = (TH1F*)CJetPtCSV4_->GetPrimitive("JetPtCSV4");
  TH1F* JetPtCSV5_ = (TH1F*)CJetPtCSV5_->GetPrimitive("JetPtCSV5");
  TH1F* JetPtCSV6_ = (TH1F*)CJetPtCSV6_->GetPrimitive("JetPtCSV6");
  
//  TH1F* JetPtDiff_ = (TH1F*)CJetPtDiff_->GetPrimitive("JetPtDiff");
  TH1F* JetPtDiffCSV1_ = (TH1F*)CJetPtDiffCSV1_->GetPrimitive("JetPtDiffCSV1");
  TH1F* JetPtDiffCSV2_ = (TH1F*)CJetPtDiffCSV2_->GetPrimitive("JetPtDiffCSV2");
  TH1F* JetPtDiffCSV3_ = (TH1F*)CJetPtDiffCSV3_->GetPrimitive("JetPtDiffCSV3");
  TH1F* JetPtDiffCSV4_ = (TH1F*)CJetPtDiffCSV4_->GetPrimitive("JetPtDiffCSV4");
  TH1F* JetPtDiffCSV5_ = (TH1F*)CJetPtDiffCSV5_->GetPrimitive("JetPtDiffCSV5");
  TH1F* JetPtDiffCSV6_ = (TH1F*)CJetPtDiffCSV6_->GetPrimitive("JetPtDiffCSV6");
  
//  TH1F* JetDR_ = (TH1F*)CJetDR_->GetPrimitive("JetDR");
  TH1F* JetDRCSV1_ = (TH1F*)CJetDRCSV1_->GetPrimitive("JetDRCSV1");
  TH1F* JetDRCSV2_ = (TH1F*)CJetDRCSV2_->GetPrimitive("JetDRCSV2");
  TH1F* JetDRCSV3_ = (TH1F*)CJetDRCSV3_->GetPrimitive("JetDRCSV3");
  TH1F* JetDRCSV4_ = (TH1F*)CJetDRCSV4_->GetPrimitive("JetDRCSV4");
  TH1F* JetDRCSV5_ = (TH1F*)CJetDRCSV5_->GetPrimitive("JetDRCSV5");
  TH1F* JetDRCSV6_ = (TH1F*)CJetDRCSV6_->GetPrimitive("JetDRCSV6");

cout << "Histograms assigned." << endl; 

  TCanvas JetPtCSVALLCanvas("JetPtCSVALL","",600,600);
  TCanvas JetPtDiffCSVALLCanvas("JetPtDiffCSVALL","",600,600);
  TCanvas JetDRCSVALLCanvas("JetDRCSVALL","",600,600);

cout << "Canvases created" << endl;

//  JetPt_->SetLineColor(kBlack);
  JetPtCSV1_->SetLineColor(kGray+2);
  JetPtCSV2_->SetLineColor(kMagenta+2);
  JetPtCSV3_->SetLineColor(3);
  JetPtCSV4_->SetLineColor(4);
  JetPtCSV5_->SetLineColor(kOrange);
  JetPtCSV6_->SetLineColor(kMagenta-4);

//  JetPtDiff_->SetLineColor(kBlack);
  JetPtDiffCSV1_->SetLineColor(kGray+2);
  JetPtDiffCSV2_->SetLineColor(kMagenta+2);
  JetPtDiffCSV3_->SetLineColor(3);
  JetPtDiffCSV4_->SetLineColor(4);
  JetPtDiffCSV5_->SetLineColor(kOrange);
  JetPtDiffCSV6_->SetLineColor(kMagenta-4);

//  JetDR_->SetLineColor(kBlack);
  JetDRCSV1_->SetLineColor(kGray+2);
  JetDRCSV2_->SetLineColor(kMagenta+2);
  JetDRCSV3_->SetLineColor(3);
  JetDRCSV4_->SetLineColor(4);
  JetDRCSV5_->SetLineColor(kOrange);
  JetDRCSV6_->SetLineColor(kMagenta-4);

cout << "Attributes set." << endl;  

  int CSV1_hits = JetPtDiffCSV1_->GetEntries(), CSV2_hits = JetPtDiffCSV2_->GetEntries(), CSV3_hits = JetPtDiffCSV3_->GetEntries(), CSV4_hits = JetPtDiffCSV4_->GetEntries();
  int CSV5_hits = JetPtDiffCSV5_->GetEntries(), CSV6_hits = JetPtDiffCSV6_->GetEntries();//, NOCSV_hits = JetPt_->GetEntries();

  TString leg1 = Form("CSV > .1: NEntries=%d",CSV1_hits), leg2 = Form("CSV > .2: NEntries=%d",CSV2_hits), leg3 = Form("CSV > .3: NEntries=%d",CSV3_hits), leg4 = Form("CSV > .4: NEntries=%d",CSV4_hits);
  TString leg5 = Form("CSV > .5: NEntries=%d",CSV5_hits), leg6 = Form("CSV > .6: NEntries=%d",CSV6_hits);//, legNOCSV = Form("NO CSV CUT : NEntries=%d",NOCSV_hits);

  JetPtCSV2_->SetXTitle("pt(PAT b-Jet)");
  JetPtDiffCSV1_->SetXTitle("#Delta pt(PAT and matched Gen b-jet)");
  JetDRCSV1_->SetXTitle("#Delta R(PAT and matched Gen b-jet)");

  legCSV = new TLegend(0.1,0.7,0.25,0.9);
//  legCSV->AddEntry(JetPt_, legNOCSV,"L");
  legCSV->AddEntry(JetPtCSV1_,leg1,"L");
  legCSV->AddEntry(JetPtCSV2_,leg2,"L");
  legCSV->AddEntry(JetPtCSV3_,leg3,"L");
  legCSV->AddEntry(JetPtCSV4_,leg4,"L");
  legCSV->AddEntry(JetPtCSV5_,leg5,"L");
  legCSV->AddEntry(JetPtCSV6_,leg6,"L");

  JetPtCSVALLCanvas.cd();
//  JetPt_->DrawNormalized();
  JetPtCSV2_->DrawNormalized();
  JetPtCSV1_->DrawNormalized("SAME");
  JetPtCSV3_->DrawNormalized("SAME");
  JetPtCSV4_->DrawNormalized("SAME");
  JetPtCSV5_->DrawNormalized("SAME");
  JetPtCSV6_->DrawNormalized("SAME");
  legCSV->Draw();
  
  JetPtDiffCSVALLCanvas.cd();
//  JetPtDiff_->DrawNormalized();
  JetPtDiffCSV1_->DrawNormalized();
  JetPtDiffCSV2_->DrawNormalized("SAME");
  JetPtDiffCSV3_->DrawNormalized("SAME");
  JetPtDiffCSV4_->DrawNormalized("SAME");
  JetPtDiffCSV5_->DrawNormalized("SAME");
  JetPtDiffCSV6_->DrawNormalized("SAME");
  legCSV->Draw();
  
  JetDRCSVALLCanvas.cd();
//  JetDR_->DrawNormalized();
  JetDRCSV1_->DrawNormalized();
  JetDRCSV2_->DrawNormalized("SAME");
  JetDRCSV3_->DrawNormalized("SAME");
  JetDRCSV4_->DrawNormalized("SAME");
  JetDRCSV5_->DrawNormalized("SAME");
  JetDRCSV6_->DrawNormalized("SAME");
  legCSV->Draw();

cout << "Histograms Drawn" << endl;

  outFile->cd();
  JetPtCSVALLCanvas.Write();
  JetPtDiffCSVALLCanvas.Write();
  JetDRCSVALLCanvas.Write();
  outFile->Write();
  outFile->Close();

cout << "end" << endl;

}//rootMacro_BBA_combine

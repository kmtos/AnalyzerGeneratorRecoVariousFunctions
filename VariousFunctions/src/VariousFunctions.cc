#include "AnalyzerGeneratorRecoVariousFunctions/VariousFunctions/interface/VariousFunctions.h"
#include <math.h>
#include "TH1F.h"
#include "THistPainter.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TASImage.h"

//
//Format and Draw 1D Histograms
//
void VariousFunctions::formatAndDrawCanvasAndHist1D(TCanvas& canvas, TH1F* hist, const Int_t grid, const Int_t logY, const Int_t logZ, const Color_t color, const Size_t size, const Style_t style,
                                          const char* xAxisTitle, const Float_t xTitleSize, const Float_t xLabelSize, const Float_t xTitleOffset, const char* yAxisTitle, const Float_t yTitleSize,
                                          const Float_t yLabelSize, const Float_t yTitleOffset, const bool errorBars)
{
  canvas.SetFillStyle(0);
  canvas.SetFillColor(0);
  canvas.SetGrid(grid, grid);
  canvas.SetLogy(logY);
  canvas.SetLogz(logZ);
  hist->SetMarkerColor(color);
  hist->SetMarkerSize(size);
  hist->SetMarkerStyle(style);
  hist->SetLineColor(color);
  hist->SetLineWidth(2);
  
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelSize(xLabelSize);
  hist->GetXaxis()->SetTitle(xAxisTitle);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(xTitleSize);
  hist->GetXaxis()->SetTitleOffset(xTitleOffset);
  hist->GetXaxis()->SetTitle(xAxisTitle);
  
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(yLabelSize);
  hist->GetYaxis()->SetTitle(yAxisTitle);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleSize(yTitleSize);
  hist->GetYaxis()->SetTitleOffset(yTitleOffset);
  hist->GetYaxis()->SetTitle(yAxisTitle);

  canvas.cd();
  if(errorBars)
    hist->Draw("E1");
  else
    hist->Draw();
}//VariousFunctions::formatAndDrawCanvasAndHist1D

//
//Format and Draw 1D Histograms
//
void VariousFunctions::formatAndDrawCanvasAndHist1DFileService(TCanvas* canvas, TH1F* hist, const Int_t grid, const Int_t logY, const Int_t logZ, const Color_t color, const Size_t size, 
					  const Style_t style, const char* xAxisTitle, const Float_t xTitleSize, const Float_t xLabelSize, const Float_t xTitleOffset, const char* yAxisTitle, 
				 	  const Float_t yTitleSize, const Float_t yLabelSize, const Float_t yTitleOffset, const bool errorBars)
{ 
  canvas->SetFillStyle(0);
  canvas->SetFillColor(0);
  canvas->SetGrid(grid, grid);
  canvas->SetLogy(logY);
  canvas->SetLogz(logZ);
  hist->SetMarkerColor(color);
  hist->SetMarkerSize(size);
  hist->SetMarkerStyle(style);
  hist->SetLineColor(color);
  hist->SetLineWidth(2);
  
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelSize(xLabelSize);
  hist->GetXaxis()->SetTitle(xAxisTitle);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(xTitleSize);
  hist->GetXaxis()->SetTitleOffset(xTitleOffset);
  hist->GetXaxis()->SetTitle(xAxisTitle);
  
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(yLabelSize);
  hist->GetYaxis()->SetTitle(yAxisTitle);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleSize(yTitleSize);
  hist->GetYaxis()->SetTitleOffset(yTitleOffset);
  hist->GetYaxis()->SetTitle(yAxisTitle);
  
  canvas->cd();
  if(errorBars)
    hist->Draw("E1");
  else
    hist->Draw();
}//VariousFunctions::formatAndDrawCanvasAndHist1DFileService


//
//Format and Draw 2D Histograms
//
void VariousFunctions::formatAndDrawCanvasAndHist2D(TCanvas& canvas, TH2F* hist, const Int_t grid, const Int_t logY, const Int_t logZ, const Color_t color, const Size_t size, const Style_t style,
 						   const char* xAxisTitle, const Float_t xTitleSize, const Float_t xLabelSize, const Float_t xTitleOffset,
						   const char* yAxisTitle, const Float_t yTitleSize, const Float_t yLabelSize, const Float_t yTitleOffset,
                                                   const char* zAxisTitle, const Float_t zTitleSize, const Float_t zLabelSize, const Float_t zTitleOffset)
{
  canvas.SetFillStyle(0);
  canvas.SetFillColor(0);
  canvas.SetGrid(grid, grid);
  canvas.SetLogy(logY);
  canvas.SetLogz(logZ);
  canvas.cd()->SetLeftMargin(.2);
  canvas.cd()->SetTopMargin(.2);
  canvas.cd()->SetRightMargin(.2);
  canvas.cd()->SetBottomMargin(.2);

  hist->SetMarkerColor(color);
  hist->SetMarkerSize(size);
  hist->SetMarkerStyle(style);
  hist->SetLineColor(color);
  hist->SetLineWidth(1);
  hist->SetFillStyle(0);

  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelOffset(0.007);
  hist->GetXaxis()->SetLabelSize(xLabelSize);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(xTitleSize);
  hist->GetXaxis()->SetTitleOffset(xTitleOffset);
  hist->GetXaxis()->SetTitle(xAxisTitle);

  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelOffset(0.007);
  hist->GetYaxis()->SetLabelSize(yLabelSize);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleSize(yTitleSize);
  hist->GetYaxis()->SetTitleOffset(yTitleOffset);
  hist->GetYaxis()->SetTitle(yAxisTitle);

  hist->GetZaxis()->SetLabelFont(42);
  hist->GetZaxis()->SetLabelOffset(0.007);
  hist->GetZaxis()->SetLabelSize(zLabelSize);
  hist->GetZaxis()->SetTitleFont(42);
  hist->GetZaxis()->SetTitleSize(zTitleSize);
  hist->GetZaxis()->SetTitleOffset(zTitleOffset);
  hist->GetZaxis()->SetTitle(zAxisTitle);

  canvas.cd();
  hist->Draw("COLZ");

}//VariousFunctions::formatAndDrawCanvasAndHist1D

//
//Format and Draw 2D Histograms
//
void VariousFunctions::formatAndDrawCanvasAndHist2DFileService(TCanvas* canvas, TH2F* hist, const Int_t grid, const Int_t logY, const Int_t logZ, const Color_t color, const Size_t size, 
						   const Style_t style,
                                                   const char* xAxisTitle, const Float_t xTitleSize, const Float_t xLabelSize, const Float_t xTitleOffset,
                                                   const char* yAxisTitle, const Float_t yTitleSize, const Float_t yLabelSize, const Float_t yTitleOffset,
                                                   const char* zAxisTitle, const Float_t zTitleSize, const Float_t zLabelSize, const Float_t zTitleOffset)
{
  canvas->SetFillStyle(0);
  canvas->SetFillColor(0);
  canvas->SetGrid(grid, grid);
  canvas->SetLogy(logY);
  canvas->SetLogz(logZ);
  canvas->cd()->SetLeftMargin(.2);
  canvas->cd()->SetTopMargin(.2);
  canvas->cd()->SetRightMargin(.2);
  canvas->cd()->SetBottomMargin(.2);

  hist->SetMarkerColor(color);
  hist->SetMarkerSize(size);
  hist->SetMarkerStyle(style);
  hist->SetLineColor(color);
  hist->SetLineWidth(1);
  hist->SetFillStyle(0);

  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelOffset(0.007);
  hist->GetXaxis()->SetLabelSize(xLabelSize);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(xTitleSize);
  hist->GetXaxis()->SetTitleOffset(xTitleOffset);
  hist->GetXaxis()->SetTitle(xAxisTitle);

  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelOffset(0.007);
  hist->GetYaxis()->SetLabelSize(yLabelSize);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleSize(yTitleSize);
  hist->GetYaxis()->SetTitleOffset(yTitleOffset);
  hist->GetYaxis()->SetTitle(yAxisTitle);

  hist->GetZaxis()->SetLabelFont(42);
  hist->GetZaxis()->SetLabelOffset(0.007);
  hist->GetZaxis()->SetLabelSize(zLabelSize);
  hist->GetZaxis()->SetTitleFont(42);
  hist->GetZaxis()->SetTitleSize(zTitleSize);
  hist->GetZaxis()->SetTitleOffset(zTitleOffset);
  hist->GetZaxis()->SetTitle(zAxisTitle);

  canvas->cd();
  hist->Draw("COLZ");

}//VariousFunctions::formatAndDrawCanvasAndHist1D

void VariousFunctions::formatAndDrawCanvasAndTGraphAsym(TCanvas& canvas, TGraphAsymmErrors* graph, const Int_t grid, const Int_t logY, const Int_t logZ, const Color_t color, const Size_t size, const Style_t style,  
						  const char* xAxisTitle, const Float_t xTitleSize, const Float_t xLabelSize, const Float_t xTitleOffset, const char* yAxisTitle, const Float_t yTitleSize, const Float_t yLabelSize,
                                                  const Float_t yTitleOffset, const bool errorBars)
{ 
  canvas.SetFillStyle(0);
  canvas.SetFillColor(0);
  canvas.SetGrid(grid, grid);
  canvas.SetLogy(logY);
  canvas.SetLogz(logZ);
  graph->SetMarkerColor(color);
  graph->SetMarkerSize(size);
  graph->SetMarkerStyle(style);
  graph->SetLineColor(color);
  graph->SetLineWidth(2);
  
  graph->GetXaxis()->SetLabelFont(42);
  graph->GetXaxis()->SetLabelSize(xLabelSize);
  graph->GetXaxis()->SetTitle(xAxisTitle);
  graph->GetXaxis()->SetTitleFont(42);
  graph->GetXaxis()->SetTitleSize(xTitleSize);
  graph->GetXaxis()->SetTitleOffset(xTitleOffset);
  graph->GetXaxis()->SetTitle(xAxisTitle);
  
  graph->GetYaxis()->SetLabelFont(42);
  graph->GetYaxis()->SetLabelSize(yLabelSize);
  graph->GetYaxis()->SetTitle(yAxisTitle);
  graph->GetYaxis()->SetTitleFont(42);
  graph->GetYaxis()->SetTitleSize(yTitleSize);
  graph->GetYaxis()->SetTitleOffset(yTitleOffset);
  graph->GetYaxis()->SetTitle(yAxisTitle);
  
  canvas.cd();
  if(errorBars)
    graph->Draw("E1");
  else
    graph->Draw("ap*");
}//VariousFunctions::formatAndDrawCanvasAndHist1D


//
//This finds and returns the GenParticleRef a daughter with the givne pdgId of a provided mother
//
reco::GenParticleRef  VariousFunctions::findDaughterInDaughters(const reco::GenParticleRef& mother, const double pdgId, const bool abs)
{
  unsigned int iDaughters = 0;
  reco::GenParticleRef  found;
  reco::GenParticleRef childRef;
  if(abs)
  {
    while (iDaughters < mother->numberOfDaughters())
    {
      childRef = mother->daughterRef(iDaughters);
      if (fabs(childRef->pdgId()) == pdgId)
      {
	found = childRef;
        break;
      }//if fabs childRef
      iDaughters++;
    }//for iGenDaughters
  }//ifabs
 else 
  { 
    while (iDaughters < mother->numberOfDaughters())
    { 
      childRef = mother->daughterRef(iDaughters);
      if (childRef->pdgId() == pdgId)
      { 
        found = childRef;
        break;
      }//if fabs childRef
      iDaughters++;
    }//for iGenDaughters
  }//else
  return found;
}//findDaughter



bool VariousFunctions::findIfInDaughters(const reco::GenParticleRef& mother, const double pdgId, const bool abs)
{
  unsigned int iDaughters = 0;
  bool found = false;
  reco::GenParticleRef childRef;
  if(abs)
  {
    while (iDaughters < mother->numberOfDaughters())
    {
      childRef = mother->daughterRef(iDaughters);
      if (fabs(childRef->pdgId()) == pdgId)
      {
        found = true;
        break;
      }//if fabs childRef
      iDaughters++;
    }//for iGenDaughters
  }//ifabs
 else
  {
    while (iDaughters < mother->numberOfDaughters())
    {
      childRef = mother->daughterRef(iDaughters);
      if (childRef->pdgId() == pdgId)
      {
        found = true;
        break;
      }//if fabs childRef
      iDaughters++;
    }//for iGenDaughters
  }//else
  return found;
}//findDaughter

//
//Takes a tau that has had or mu or ele decay products
//Tau Decay Mode: 1= 1prong , 2= 1pr + 1pi_0 , 3= 1pr + 2pi_0 , 4=3pr , 5=other , 6=electronic , 7=muonic
//
int VariousFunctions::tauDecayMode(const reco::GenParticleRef& tauRef)
{
  if(tauRef->numberOfDaughters() == 3)
  {
  reco::GenParticleRef daughter1 = tauRef->daughterRef(0);
  reco::GenParticleRef daughter2 = tauRef->daughterRef(1);
  reco::GenParticleRef daughter3 = tauRef->daughterRef(2);

  //Leptonic checks
  if(fabs(daughter1->pdgId()) == 13 || fabs(daughter2->pdgId()) == 13 || fabs(daughter3->pdgId()) == 13) 
    return 7;
  if(fabs(daughter1->pdgId()) == 11 || fabs(daughter2->pdgId()) == 11 || fabs(daughter3->pdgId()) == 11) 
    return 6;
  }//numof Daughters == 3

  unsigned int iDaughters = 0, nProngs = 0, nPi0=0;
  while(iDaughters < tauRef->numberOfDaughters())
  {
    reco::GenParticleRef interDaughter = tauRef->daughterRef(iDaughters);
    if(interDaughter->pdgId() == 111 || interDaughter->pdgId() == 310 || interDaughter->pdgId() == 130) //checks for pi_0, KS_0, or KL_0
      nPi0++;
    if(fabs(interDaughter->pdgId()) == 211 || fabs(interDaughter->pdgId()) == 321)  //checks for +-pion or +- Kaon
      nProngs++;
    if(fabs(interDaughter->pdgId()) == 213 || fabs(interDaughter->pdgId()) == 223 || fabs(interDaughter->pdgId()) == 20213)  //checks for +-rho and omega_0
    {
      unsigned int iDaughters2 = 0;
      while(iDaughters2 < interDaughter->numberOfDaughters())
      {
	reco::GenParticleRef iGrandDaughter = interDaughter->daughterRef(iDaughters2);
        if(iGrandDaughter->pdgId() == 111 || iGrandDaughter->pdgId() == 310 || iGrandDaughter->pdgId() == 130) //checks for rho-> pi_0, KS_0, or KL_0
          nPi0++;
        if(fabs(iGrandDaughter->pdgId()) == 211 || fabs(iGrandDaughter->pdgId()) == 321)//checks for rho->+-pi or +- kaon
          nProngs++;
    	iDaughters2++;
      }//while
    }//iffabs
    iDaughters++;
  }//while

  if(nProngs == 1 && nPi0 == 0)
    return 1;
  if(nProngs == 1 && nPi0 == 1)
    return 2;
  if(nProngs == 1 && nPi0 == 2)
    return 3;
  if(nProngs == 3 && nPi0 == 0)
    return 4;
  else
    return 5;
}//tauDecayMode

//
//
//Takes a tau Ref and returns the Decay Products
//
//
reco::LeafCandidate::LorentzVector VariousFunctions::sumTauP4(const reco::GenParticleRef& particleRef, const int decayMode, const bool pi0_decay)
{
 reco::LeafCandidate::LorentzVector tauP4; 
  
  if(decayMode == 1)
  {
    reco::GenParticleRef prong1;
    if(fabs(particleRef->daughter(0)->pdgId()) == 211 || fabs(particleRef->daughter(0)->pdgId()) == 321)
      prong1 = particleRef->daughterRef(0);
    else
      prong1 = particleRef->daughterRef(1);
    tauP4 = prong1->p4(); 
  }//if decayMode ==1 

  if(decayMode == 2)
  {
    reco::GenParticleRef neutPart, prong1;
    unsigned int iDaughters = 0;
    while(iDaughters < particleRef->numberOfDaughters())
    {
      if(fabs(particleRef->daughter(iDaughters)->pdgId()) == 213 || fabs(particleRef->daughter(iDaughters)->pdgId() ) == 223 || fabs(particleRef->daughter(iDaughters)->pdgId() ) == 20213)
      {
	unsigned int iDaughters2 = 0;
	reco::GenParticleRef newRef = particleRef->daughterRef(iDaughters);
	while(iDaughters2 < newRef->numberOfDaughters())
	{
	  if(newRef->daughter(iDaughters2)->pdgId() == 111 || newRef->daughter(iDaughters2)->pdgId() == 130 || newRef->daughter(iDaughters2)->pdgId() == 310)
	    neutPart = newRef->daughterRef(iDaughters2);
	  if(fabs(newRef->daughter(iDaughters2)->pdgId()) == 211 || fabs(newRef->daughter(iDaughters2)->pdgId()) == 321)
	    prong1 = newRef->daughterRef(iDaughters2);
	  iDaughters2++;
	}//while Daughters2
      }//if rho or omega
      if(particleRef->daughter(iDaughters)->pdgId() == 111 || particleRef->daughter(iDaughters)->pdgId() == 130 || particleRef->daughter(iDaughters)->pdgId() == 310)
        neutPart = particleRef->daughterRef(iDaughters);
      if(fabs(particleRef->daughter(iDaughters)->pdgId()) == 211 || fabs(particleRef->daughter(iDaughters)->pdgId()) == 321)
        prong1 = particleRef->daughterRef(iDaughters);
      iDaughters++;
    }//while iDaughters

    if(pi0_decay && neutPart->pdgId() == 111)
      tauP4 = neutPart->daughter(0)->p4() + neutPart->daughter(1)->p4() + prong1->p4();
    else
      tauP4 = neutPart->p4() + prong1->p4();
  }//if decayMode == 2

  if(decayMode == 3)
  {
    reco::GenParticleRef neutPart1, neutPart2, prong1;
    unsigned int iDaughters = 0, check = 0; 
    while(iDaughters < particleRef->numberOfDaughters())
    {
      if(fabs(particleRef->daughter(iDaughters)->pdgId()) == 213 || fabs(particleRef->daughter(iDaughters)->pdgId() ) == 223 || fabs(particleRef->daughter(iDaughters)->pdgId() ) == 20213)
      {
        unsigned int iDaughters2 = 0;
        reco::GenParticleRef newRef = particleRef->daughterRef(iDaughters);
        while(iDaughters2 < newRef->numberOfDaughters())
        {
          if((newRef->daughter(iDaughters2)->pdgId() == 111 || newRef->daughter(iDaughters2)->pdgId() == 130 || newRef->daughter(iDaughters2)->pdgId() == 310) && check == 1 )
 	  {
            neutPart2 = newRef->daughterRef(iDaughters2);
	    check++;
	  }//if check==0
          if((newRef->daughter(iDaughters2)->pdgId() == 111 || newRef->daughter(iDaughters2)->pdgId() == 130 || newRef->daughter(iDaughters2)->pdgId() == 310) && check == 0 )
	  {
            neutPart1 = newRef->daughterRef(iDaughters2);
	    check++;
	  }//check==0
          if(fabs(newRef->daughter(iDaughters2)->pdgId()) == 211 || fabs(newRef->daughter(iDaughters2)->pdgId()) == 321)
            prong1 = newRef->daughterRef(iDaughters2);
          iDaughters2++;
        }//while Daughters2
      }//if rho or omega
      if( (particleRef->daughter(iDaughters)->pdgId() == 111 || particleRef->daughter(iDaughters)->pdgId() == 130 || particleRef->daughter(iDaughters)->pdgId() == 310 )  && check == 1 )
      {
        neutPart2 = particleRef->daughterRef(iDaughters);
	check++;
      }//check==0
      if( (particleRef->daughter(iDaughters)->pdgId() == 111 || particleRef->daughter(iDaughters)->pdgId() == 130 || particleRef->daughter(iDaughters)->pdgId() == 310 )  && check == 0 )
      { 
        neutPart1 = particleRef->daughterRef(iDaughters);
        check++;
      }//check==0
      if(fabs(particleRef->daughter(iDaughters)->pdgId()) == 211 || fabs(particleRef->daughter(iDaughters)->pdgId()) == 321)
        prong1 = particleRef->daughterRef(iDaughters);
      iDaughters++;
    }//while iDaughters

    if(pi0_decay)
    {
      if(neutPart1->pdgId() == 111 && neutPart2->pdgId() != 111)
        tauP4 = prong1->p4() + neutPart1->daughter(0)->p4() + neutPart1->daughter(1)->p4() + neutPart2->p4();
      if(neutPart2->pdgId() == 111 && neutPart1->pdgId() != 111)
	tauP4 = prong1->p4() + neutPart1->p4() + neutPart2->daughter(0)->p4() + neutPart2->daughter(1)->p4();
      if(neutPart2->pdgId() == 111 && neutPart1->pdgId() == 111)
        tauP4 = prong1->p4() + neutPart1->daughter(0)->p4() + neutPart1->daughter(1)->p4() + neutPart2->daughter(0)->p4() + neutPart2->daughter(1)->p4(); 
      else
        tauP4 = neutPart1->p4() + neutPart2->p4() + prong1->p4();
    }//if
    else
      tauP4 = neutPart1->p4() + neutPart2->p4() + prong1->p4();
  }//if decayMode ==3 


  if(decayMode == 4)
  {
    reco::GenParticleRef prong1, prong2, prong3;
    unsigned int iDaughters = 0, check = 0;
    while(iDaughters < particleRef->numberOfDaughters())
    {
      if(fabs(particleRef->daughter(iDaughters)->pdgId()) == 213 || fabs(particleRef->daughter(iDaughters)->pdgId() ) == 223  || fabs(particleRef->daughter(iDaughters)->pdgId() ) == 20213)
      {
	unsigned int iDaughters2 = 0;
	reco::GenParticleRef newRef = particleRef->daughterRef(iDaughters);
	while(iDaughters2 < newRef->numberOfDaughters())
	{  
          if( (fabs(newRef->daughter(iDaughters2)->pdgId()) == 211 || fabs(newRef->daughter(iDaughters2)->pdgId()) == 321 ) && check == 0)
          {
            prong1 = newRef->daughterRef(iDaughters2);
	    check++;
          }//if check ==0
          if( (fabs(newRef->daughter(iDaughters2)->pdgId()) == 211 || fabs(newRef->daughter(iDaughters2)->pdgId()) == 321 ) && check == 1)
          {
            prong2 = newRef->daughterRef(iDaughters2);
	    check++;
          }//check == 1
          if( (fabs(newRef->daughter(iDaughters2)->pdgId()) == 211 || fabs(newRef->daughter(iDaughters2)->pdgId()) == 321 ) && check == 2)
   	    prong3 = newRef->daughterRef(iDaughters2);
          iDaughters2++;
	}//while iDaughters2
      }//if not final state particle
      if( (fabs(particleRef->daughter(iDaughters)->pdgId()) == 211 || fabs(particleRef->daughter(iDaughters)->pdgId()) == 321 ) && check == 0)
      {
	prong1 = particleRef->daughterRef(iDaughters);
	check++;
      }// if check == 0
      if( (fabs(particleRef->daughter(iDaughters)->pdgId()) == 211 || fabs(particleRef->daughter(iDaughters)->pdgId()) == 321 ) && check == 1)
      {
	prong2 = particleRef->daughterRef(iDaughters);
	check++;
      }//check == 1
      if( (fabs(particleRef->daughter(iDaughters)->pdgId()) == 211 || fabs(particleRef->daughter(iDaughters)->pdgId()) == 321 ) && check == 2)
      {
	prong3 = particleRef->daughterRef(iDaughters);
	check++;
      }//if 2
      iDaughters++;
    }//while iDaughters

    tauP4 = prong1->p4() + prong2->p4() + prong3->p4();
  }//if decayMode == 4


  if(decayMode == 6)
  {
    reco::GenParticleRef electron;
    unsigned int nDaughters = particleRef->numberOfDaughters(), iDaughters=0;
    while (iDaughters < nDaughters)
    {
      reco::GenParticleRef iChild = particleRef->daughterRef(iDaughters);
      if(fabs(iChild->pdgId()) == 11)
	electron = iChild;
      iDaughters++;
    }//while  

    tauP4 = electron->p4();
  }//of decayMode == 6


  if(decayMode == 7)
  {
    reco::GenParticleRef muon;
    unsigned int nDaughters = particleRef->numberOfDaughters(), iDaughters=0;
    while (iDaughters < nDaughters) 
    {
      reco::GenParticleRef iChild = particleRef->daughterRef(iDaughters);
      if(fabs(iChild->pdgId()) == 13)
        muon = iChild;
      iDaughters++;
    }//while  
 
    tauP4 = muon->p4();
  }//if decayMode == 7

  return tauP4; 
}//sumDaughtersPt

//
//Calculates the dR of the diTau object
//
double VariousFunctions::getDiTauDR(const reco::GenParticleRef& tau1Ref, const reco::GenParticleRef& tau2Ref, const bool piDecay)
{
  int tau1DecayMode = VariousFunctions::tauDecayMode(tau1Ref), tau2DecayMode = VariousFunctions::tauDecayMode(tau2Ref);
  reco::LeafCandidate::LorentzVector tau1P4 = sumTauP4(tau1Ref, tau1DecayMode, piDecay), tau2P4 = VariousFunctions::sumTauP4(tau2Ref, tau2DecayMode, piDecay);
  double dPhi = reco::deltaPhi(tau1P4.Phi(), tau2P4.Phi() );
  double aDR = sqrt( (tau1P4.Eta() - tau2P4.Eta() )*(tau1P4.Eta() - tau2P4.Eta() )  +  (dPhi )*(dPhi ) );
  std::cout << "eta1= " << tau1P4.Eta() << " eta2= " << tau2P4.Eta() << " phi1= " << tau1P4.Phi() << " phi2= " << tau2P4.Phi() << " dr= " << aDR << std::endl;
  return aDR;
}//VariousFunctions::getDiTauDR

//
//Calculates the dR between a and b
//
double VariousFunctions::getABDR(const double bEta, const double bPhi, const reco::GenParticleRef& tau1Ref, const reco::GenParticleRef& tau2Ref, const bool piDecay)
{
  int tau1DecayMode = VariousFunctions::tauDecayMode(tau1Ref), tau2DecayMode = VariousFunctions::tauDecayMode(tau2Ref);
  reco::LeafCandidate::LorentzVector tau1P4 = VariousFunctions::sumTauP4(tau1Ref, tau1DecayMode, piDecay),  tau2P4 = VariousFunctions::sumTauP4(tau2Ref, tau2DecayMode, piDecay);
  reco::LeafCandidate::LorentzVector aP4 = tau1P4 + tau2P4;
  double dPhi = reco::deltaPhi(bPhi, aP4.Phi() );
  double aBDR = sqrt( (bEta - aP4.Eta() )*(bEta - aP4.Eta() )  +  (dPhi )*(dPhi ) );
  return aBDR;
}//VariousFunctions::getABDR

double* VariousFunctions::orderFour(const double o1, double const o2, const double o3, const double o4)
{
  static double a[4];
  double  o1o1 = 0, o1o2 = 0, o2o1 = 0, o2o2 = 0;
  if (o1 > o2 )
  {
    o1o1 = o1;
    o1o2 = o2;
  }//if 11 > 12
  else
  {
    o1o1 = o2;
    o1o2 = o1;
  }//else 12 > 11
  if (o3 > o4 )
  {
    o2o1 = o3;
    o2o2 = o4;
  }//if 21 > 22
  else
  {
    o2o1 = o4;
    o2o2 = o3;
  }//else 22 > 21
  if (o1o1 > o2o1)
  {
    a[0] = o1o1;
    if (o1o2 > o2o2)
    {
      a[3] = o2o2;
      if( o1o2 > o2o1 )
      {
        a[1] = o1o2;
        a[2] = o2o1;
      }//if o1o2 > o2o1
      else
      {
        a[1] = o2o1;
        a[2] = o1o2;
      }//else o2o1 > o1o2
    }//if o1o2 > o2o2
    else
    {
      a[3] = o1o2;
      a[1] = o2o1;
      a[2] = o2o2;
    }//else o2o2 > o1o2
  }// if o1o1 > o2o1
  else
  {
    a[0] = o2o1;
    if (o1o2 > o2o2)
    {
      a[3] = o2o2;
      a[1] = o1o1;
      a[2] = o1o2;
    }// o1o2 > o2o2
    else
    {
      a[3] = o1o2;
      if (o1o1 > o2o2)
      {
        a[1] = o1o1;
        a[2] = o2o2;
      }//if o1o1 >  o2o2
      else
      {
        a[1] = o2o2;
        a[2] = o1o1;
      }//else o2o2 > o1o1
    }//else o2o1 > o1o2
  }//else o1o1 < o2o1

 return a;
}//VariousFunctions::orderFour


void VariousFunctions::findAndPlotBMuons(const reco::GenParticleRef& bRef, const int treeDepth, TH1F* hist, const bool displayTree)
{
  for(unsigned int i = 0; i < bRef->numberOfDaughters(); i++)
  {
    reco::GenParticleRef childRef = bRef->daughterRef(i);
    for(int j = 0; j < treeDepth; j++)
      std::cout << "\t";
    if(displayTree)
      std::cout << "child #" << i << " has pdgId= " << childRef->pdgId() << " and has #" << childRef->numberOfDaughters() << " daughters." << std::endl;
    if(childRef->numberOfDaughters() > 1 || VariousFunctions::findIfInDaughters(bRef, 5, true))
      VariousFunctions::findAndPlotBMuons(childRef, treeDepth + 1, hist, displayTree);
    if(fabs(childRef->pdgId() ) == 15 )
    {
      std::cout << "\t\t<----BMOUN----->" << std::endl;
      hist->Fill(childRef->pt() );
    }//if pdgid == 15
  }//for i


}//VariousFunctions::FindAndPlotBMuons
#include "TStyle.h"


void VariousFunctions::setTDRStyle(bool gridOn) {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  if (gridOn)
  {
    tdrStyle->SetPadGridX(true);
    tdrStyle->SetPadGridY(true);
  }//if
  else
  {
    tdrStyle->SetPadGridX(false);
    tdrStyle->SetPadGridY(false);
  }//else
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);
  
// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  // tdrStyle->SetErrorMarker(20);
  //tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);
  
//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->cd();

} 

void VariousFunctions::SetCanvas(TCanvas* canv, int canvWidth, int canvHeight)
{
  //canv->SetCanvasSize(canvWidth, canvHeight);

  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  //canv->SetLeftMargin( 0.12/canvWidth );
  //canv->SetRightMargin( .04/canvWidth );
  //canv->SetTopMargin( .08/canvHeight );
  //canv->SetBottomMargin( 0.12/canvHeight );
  canv->SetTickx(0);
  canv->SetTicky(0);
}

void VariousFunctions::SetHist(TH1F* hist, int histLineColor, int histMarkerColor, int histFillColor, float markerSize, const char* xAxisTitle, const char* yAxisTitle)
{
  hist->GetXaxis()->SetNdivisions(6,5,0);
  hist->GetXaxis()->SetTitle(xAxisTitle);
  hist->GetYaxis()->SetNdivisions(6,5,0);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetTitle(yAxisTitle);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(markerSize);
  hist->SetMarkerColor(histMarkerColor);
  hist->SetLineColor(histLineColor);
  hist->SetFillColor(histFillColor);
}

void VariousFunctions::CMS_lumi( TPad* pad, int iPeriod, int iPosX, bool outOfFrame, bool writeExtraText, bool drawLogo)
{
  TString cmsText     = "CMS";
  float cmsTextFont   = 61;  // default is helvetic-bold
  
  TString extraText   = "Preliminary";
  float extraTextFont = 52;  // default is helvetica-italics
  
  // text sizes and text offsets with respect to the top frame
  // in unit of the top margin size
  float lumiTextSize     = 0.6;
  float lumiTextOffset   = 0.2;
  float cmsTextSize      = 0.75;
  //float cmsTextOffset    = 0.1;  // only used in outOfFrame version
  
  float relPosX    = 0.045;
  float relPosY    = 0.035;
  float relExtraDY = 1.2;
  
  // ratio of "CMS" and extra text size
  float extraOverCmsTextSize  = 0.76;
  
  TString lumi_13TeV = "20.1 fb^{-1}";
  TString lumi_8TeV  = "19.7 fb^{-1}";
  TString lumi_7TeV  = "5.1 fb^{-1}";
  TString lumi_sqrtS = "";
  
  if( iPosX/10==0 )
    outOfFrame = true;

  int alignY_=3;
  int alignX_=2;
  if( iPosX/10==0 ) alignX_=1;
  if( iPosX==0    ) alignX_=1;
  if( iPosX==0    ) alignY_=1;
  if( iPosX/10==1 ) alignX_=1;
  if( iPosX/10==2 ) alignX_=2;
  if( iPosX/10==3 ) alignX_=3;
  //if( iPosX == 0  ) relPosX = 0.12;
  int align_ = 10*alignX_ + alignY_;
  
//  float H = pad->GetWh(); //Because of Non-compiling draw_logo
//  float W = pad->GetWw(); //Because of Non=compiling draw_logo
  float l = pad->GetLeftMargin();
  float t = pad->GetTopMargin();
  float r = pad->GetRightMargin();
  float b = pad->GetBottomMargin();
  //  float e = 0.025;
 
  pad->cd();
  
  TString lumiText;
  if( iPeriod==1 )
  { 
    lumiText += lumi_7TeV;
    lumiText += " (7 TeV)";
  } 
  else if ( iPeriod==2 )
  {
    lumiText += lumi_8TeV;
    lumiText += " (8 TeV)";
  } 
  else if( iPeriod==3 )
  {
    lumiText = lumi_8TeV;
    lumiText += " (8 TeV)";
    lumiText += " + ";
    lumiText += lumi_7TeV;
    lumiText += " (7 TeV)";
  } 
  else if ( iPeriod==4 )
  {
    lumiText += lumi_13TeV;
    lumiText += " (13 TeV)";
  } 
  else if ( iPeriod==7 )
  {
    if( outOfFrame ) lumiText += "#scale[0.85]{";
    lumiText += lumi_13TeV;
    lumiText += " (13 TeV)";
    lumiText += " + ";
    lumiText += lumi_8TeV;
    lumiText += " (8 TeV)";
    lumiText += " + ";
    lumiText += lumi_7TeV;
    lumiText += " (7 TeV)";
    if( outOfFrame) lumiText += "}";
  }
  else if ( iPeriod==12 )
    lumiText += "8 TeV";
  else if ( iPeriod==0 )
    lumiText += lumi_sqrtS;

  std::cout << lumiText << std::endl;


  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);
  float extraTextSize = extraOverCmsTextSize*cmsTextSize;

  latex.SetTextFont(42);
  latex.SetTextAlign(31);
  latex.SetTextSize(lumiTextSize*t);
  latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);
std::cout << "\tlumiText=" << lumiText << "\t1-t+lumiTextOffset*t=" << 1-t+lumiTextOffset*t << std::endl;

  if( outOfFrame )
  {
    latex.SetTextFont(cmsTextFont);
    latex.SetTextAlign(11);
    latex.SetTextSize(cmsTextSize*t);
    latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);
  }

  pad->cd();

  float posX_=0;
  if( iPosX%10<=1 )
    posX_ =   l + relPosX*(1-l-r);
  else if( iPosX%10==2 )
    posX_ =  l + 0.5*(1-l-r);
  else if( iPosX%10==3 )
    posX_ =  1-r - relPosX*(1-l-r);
  float posY_ = 1-t - relPosY*(1-t-b);
  if( !outOfFrame )
  {
    if( drawLogo )
    {
/*
      posX_ =   l + 0.045*(1-l-r)*W/H;
      posY_ = 1-t - 0.045*(1-t-b);
      float xl_0 = posX_;
      float yl_0 = posY_ - 0.05;
      float xl_1 = posX_ + 0.10*H/W;
      float yl_1 = posY_ + .05;
      TASImage* CMS_logo = new TASImage("/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/doc/CMS-BW-label.png");
      TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
      pad_logo->Draw();
      pad_logo->cd();
      CMS_logo->Draw("X");
      pad_logo->Modified();
      pad->cd();
*/
    }//if drawlogo
    else
    {
      latex.SetTextFont(cmsTextFont);
      latex.SetTextSize(cmsTextSize*t);
      latex.SetTextAlign(align_);
      latex.DrawLatex(posX_, posY_, cmsText);
      if( writeExtraText )
      {
        latex.SetTextFont(extraTextFont);
        latex.SetTextAlign(align_);
        latex.SetTextSize(extraTextSize*t);
        latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText);
      }//if writeExtraText
    }//else
  }//if  !outOfFrame
  else if( writeExtraText )
  {
    if( iPosX==0)
    {
      posX_ =   l +  relPosX*(1-l-r);
      posY_ =   1-t+lumiTextOffset*t;
   }//if iPosX
   latex.SetTextFont(extraTextFont);
   latex.SetTextSize(extraTextSize*t);
   latex.SetTextAlign(align_);
   latex.DrawLatex(posX_, posY_, extraText);
  }//else
  return;
}//void CMS_lumi

#include <cstdlib>
#include <iostream> 
#include <map>
#include <string>

#include "TChain.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCut.h"

#define VERBOSE true
#define SAVE true


////////////////////////////////////////////////////////

void plotMacroTest( Int_t nBins_ = 100,
		    Float_t xMin_=-0.3,
		    Float_t xMax_=0.3,
		    Int_t nBinsPt_ = 8,
		    Float_t xMinPt_=20,
		    Float_t xMaxPt_=100,
		    Float_t MVACut_ = -0.0125,//0.03 or -0.0125
		    bool ForBarrel = true
		  )
{   
  // Open the files
  
  TFile *fSignalVBF = new TFile("nTupleVBFH115-ElecTau-powheg-PUS6_run_Open_ElecTauStream.root" ,"READ");  
  TTree *signalVBF           = (TTree*)fSignalVBF->Get("outTreePtOrd");

  
  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  
  
  TLegend* leg = new TLegend(0.55,0.42,0.85,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  leg->SetHeader(Form("#splitline{CMS Preliminary}{#splitline{ #sqrt{s}=7 TeV, VBF H#rightarrow#tau#tau m_{H}=115}{RelIsoDB<0.2 and MVA>%0f for Barrel}}",MVACut_));
  if(!ForBarrel)  leg->SetHeader(Form("#splitline{CMS Preliminary}{#splitline{ #sqrt{s}=7 TeV, VBF H#rightarrow#tau#tau m_{H}=115}{RelIsoDB<0.2 and MVA>%0f for Endcap}}",MVACut_));

  TH1F* hDanieleMVA   = new TH1F( "hDanieleMVA" ,"" , nBins_ ,xMin_ , xMax_);
  hDanieleMVA->SetXTitle("MVA discriminator");
  hDanieleMVA->SetYTitle("Events");

  TH1F* hEff   = new TH1F( "hEff" ,"" , nBins_ ,xMin_ , xMax_);
  hEff->SetXTitle("Cut on MVA EleId");
  hEff->SetYTitle("Efficiency");

  TH1F* hEffvsPt   = new TH1F( "hEffvsPt" ,"" , nBinsPt_ ,xMinPt_ , xMaxPt_);
  hEffvsPt->SetMarkerColor(kBlue);
  hEffvsPt->SetMarkerStyle(kOpenCircle);
  hEffvsPt->SetMarkerSize(1.5);
  TH1F* hEffvsPtP   = new TH1F( "hEffvsPtP" ,"" , nBinsPt_ ,xMinPt_ , xMaxPt_);
  hEffvsPtP->SetXTitle("P_{T}(e)[GeV]");
  hEffvsPtP->SetYTitle("Efficiency");
  hEffvsPtP->SetMarkerColor(kBlue);
  hEffvsPtP->SetMarkerStyle(kOpenCircle);
  hEffvsPtP->SetMarkerSize(1.5);

  TArrayF bins(4);
  bins[0] = 0;    bins[1] = 6;    bins[2] = 12;    bins[3] = 30;    
  TH1F* hEffvsnumPV   = new TH1F( "hEffvsnumPV" ,"" , 3 ,bins.GetArray());
  hEffvsnumPV->SetMarkerColor(kRed);
  hEffvsnumPV->SetMarkerStyle(kOpenCircle);
  hEffvsnumPV->SetMarkerSize(1.5);
  TH1F* hEffvsnumPVP   = new TH1F( "hEffvsnumPVP" ,"" , 3 ,bins.GetArray());
  hEffvsnumPVP->SetXTitle("numPV");
  hEffvsnumPVP->SetYTitle("Efficiency");
  hEffvsnumPVP->SetMarkerColor(kRed);
  hEffvsnumPVP->SetMarkerStyle(kOpenCircle);
  hEffvsnumPVP->SetMarkerSize(1.5);

  signalVBF->Draw("tightestDanieleMVAWP>>hDanieleMVA");
  signalVBF->Draw("ptL1>>hEffvsPt",/*"combRelIsoLeg1DBeta<0.2 &&*/" abs(etaL1)<1.5 ");
  signalVBF->Draw("ptL1>>hEffvsPtP",Form("combRelIsoLeg1DBeta<0.2 && tightestDanieleMVAWP>%0f && abs(etaL1)<1.5 ",MVACut_));
  signalVBF->Draw("numPV>>hEffvsnumPV","combRelIsoLeg1DBeta<0.2 && abs(etaL1)<1.5 ");
  signalVBF->Draw("numPV>>hEffvsnumPVP",Form("combRelIsoLeg1DBeta<0.2 && tightestDanieleMVAWP>%0f && abs(etaL1)<1.5 ",MVACut_));
  //signalVBF->Draw("numPV>>hEffvsnumPV","combRelIsoLeg1DBeta<0.2"/* && abs(etaL1)<1.5 "*/);
  //signalVBF->Draw("numPV>>hEffvsnumPVP",Form("combRelIsoLeg1DBeta<0.2 && tightestDanieleMVAWP>%0f"/* && abs(etaL1)<1.5 "*/,MVACut_));
  if (!ForBarrel){
    signalVBF->Draw("ptL1>>hEffvsPt",/*"combRelIsoLeg1DBeta<0.2 && */"abs(etaL1)>1.5");
    signalVBF->Draw("ptL1>>hEffvsPtP",Form("combRelIsoLeg1DBeta<0.2 && tightestDanieleMVAWP>%0f && abs(etaL1)>1.5 ",MVACut_));
    signalVBF->Draw("numPV>>hEffvsnumPV","combRelIsoLeg1DBeta<0.2 && abs(etaL1)>1.5");
    signalVBF->Draw("numPV>>hEffvsnumPVP",Form("combRelIsoLeg1DBeta<0.2 && tightestDanieleMVAWP>%0f && abs(etaL1)>1.5 ",MVACut_));
  }

  hEffvsPtP->Sumw2();
  hEffvsPt->Sumw2();
  hEffvsPtP->Divide(hEffvsPt);

  hEffvsnumPVP->Sumw2();
  hEffvsnumPV->Sumw2();
  hEffvsnumPVP->Divide(hEffvsnumPV);

  //hDanieleMVA->Draw();
  //leg->Draw();

  //Compute the cumulative of hDanieleMVA
  int total = 0;
  for (int i=0;i<nBins_;i++)
    {
     total += hDanieleMVA->GetBinContent(i);
     for (int k=0;k<total;k++)
       {
        hEff->AddBinContent(i);
       }
    }
  double norme = 1/hDanieleMVA->Integral();
  //if(VERBOSE)cout <<"MVA discriminator Normalization factor : "<< norme << endl;
  hEff->Scale(norme);
  //Efficiency vs MVACut
  //hEff->Draw();
  //leg->Draw();

  
  // //Efficiency vs Pt(electron)
  hEffvsPtP->Draw("P");
  leg->Draw();

  //Efficiency vs numPV
  //hEffvsnumPVP->Draw("P");
  //leg->Draw();
}

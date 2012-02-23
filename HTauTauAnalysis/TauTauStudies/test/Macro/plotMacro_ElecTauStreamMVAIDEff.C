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

void plotMacro_ElecTauStreamMVAIDEff( Int_t nBins_ = 100,
			      Float_t xMin_=-0.3,
			      Float_t xMax_=0.3,
			      Int_t nBinsPt_ = 1,
			      Float_t xMinPt_=20,
			      Float_t xMaxPt_=30,//100
			      Float_t MVACut_ = 0.01625,//0.0225 for endcap or 0.01625 for barrel
			      Float_t PtCutMin_ = 20,
			      Float_t PtCutMax_ = 30,
			      bool ForBarrel = true
			      )
{   
  // Open the files
  
  TFile *fSignal = new TFile("/data_CMS/cms/ivo/Trees/ElecTauStream_08Feb2012/OpenNtuples/nTupleDYJets-ElecTau-50-madgraph-PUS6-v2_run_Open_ElecTauStream.root" ,"READ");  
  TTree *signal           = (TTree*)fSignal->Get("outTreePtOrd");

  
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
  leg->SetHeader(Form("#splitline{CMS Preliminary}{#splitline{ #sqrt{s}=7 TeV, GGF H#rightarrow#tau#tau m_{H}=115}{RelIsoDB<0.2 and MVA>%0f for Barrel}}",MVACut_));
  if(!ForBarrel)  leg->SetHeader(Form("#splitline{CMS Preliminary}{#splitline{ #sqrt{s}=7 TeV}{MVA>%0f for Endcap}}",MVACut_));

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


  signal->Draw("tightestDanieleMVAWP>>hDanieleMVA");
  signal->Draw("ptL1>>hEffvsPt",Form("puWeight*(combRelIsoLeg1DBeta<0.2 && abs(etaL1)<1.5 && ptL1>%0f && ptL1<%0f && isElecLegMatched>0.5 && HLTx>0.5 && HLTmatch>0.5)", PtCutMin_, PtCutMax_));
  signal->Draw("ptL1>>hEffvsPtP",Form("puWeight*(combRelIsoLeg1DBeta<0.2 && tightestDanieleMVAWP>%0f && abs(etaL1)<1.5 && isElecLegMatched>0.5 && ptL1>%0f && ptL1<%0f && HLTx>0.5 && HLTmatch>0.5)",MVACut_, PtCutMin_, PtCutMax_));
  if (!ForBarrel){
    signal->Draw("ptL1>>hEffvsPt",Form("puWeight*(combRelIsoLeg1DBeta<0.2 && abs(etaL1)>1.5&& ptL1>%0f && ptL1<%0f && isElecLegMatched>0.5 && HLTx>0.5 && HLTmatch>0.5 && HLTx>0.5 && HLTmatch>0.5)", PtCutMin_, PtCutMax_));
    signal->Draw("ptL1>>hEffvsPtP",Form("puWeight*(combRelIsoLeg1DBeta<0.2 && tightestDanieleMVAWP>%0f && abs(etaL1)>1.5 && ptL1>%0f && ptL1<%0f && isElecLegMatched>0.5 && HLTx>0.5 && HLTmatch>0.5)",MVACut_, PtCutMin_, PtCutMax_));
  }

  hEffvsPtP->Sumw2();
  hEffvsPt->Sumw2();
  hEffvsPtP->Divide(hEffvsPt);

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


}

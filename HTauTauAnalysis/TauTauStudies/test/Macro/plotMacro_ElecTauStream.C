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

void plotMacro_ElecTauStream( Int_t nBins_ = 100,
			      Float_t xMin_=-0.3,
			      Float_t xMax_=0.3,
			      Int_t nBinsPt_ = 8,
			      Int_t nBinsEta_ = 5,
			      Float_t xMinPt_=20,
			      Float_t xMaxPt_=100,//100
			      Float_t xMinEta_=-2.1,
			      Float_t xMaxEta_=2.1,//100
			      Float_t MVACut_ = 0.0275,//0.03 or -0.0125
			      Float_t PtCutMin_ = 20,
			      Float_t PtCutMax_ = 30,
			      bool ForBarrel = true
			      )
{   
  // Open the files
  
  TFile *fSignal = new TFile("/data_CMS/cms/ivo/Trees/ElecTauStream_16Feb2012/OpenNtuples/nTupleDYJets-ElecTau-50-madgraph-PUS6-v2_run_Open_ElecTauStream.root" ,"READ");  
  TTree *signal           = (TTree*)fSignal->Get("outTreePtOrd");

  TFile *fQCDBCtoE2080 = new TFile("/data_CMS/cms/ivo/Trees/ElecTauStream_20Feb2012/OpenNtuples/nTupleQCD-20to80-BCtoE-ElecTau-pythia-PUS6_run_Open_ElecTauStream.root" ,"READ");  
  TTree *QCDBCtoE2080           = (TTree*)fQCDBCtoE2080->Get("outTreePtOrd");



  //Eta histos for QCD estimated from data
  TFile *fEtaQCDDaniele = new TFile("/home/llr/cms/ivo/CMSSW_4_2_8_patch7/src/Bianchi/Limits/htautau/histograms/ElecTauStream_16Feb2012/eTau_mH120_inclusive__etaL1_Daniele.root" ,"READ");  
  TH1F *hEtaQCDDaniele           = (TH1F*)fEtaQCDDaniele->Get("hQCD");
  hEtaQCDDaniele->SetXTitle("#eta(e)");
  hEtaQCDDaniele->SetFillColor(0);
  hEtaQCDDaniele->SetMarkerColor(kBlue);
  hEtaQCDDaniele->SetMarkerStyle(kOpenCircle);
  hEtaQCDDaniele->SetMarkerSize(1);
  hEtaQCDDaniele->Sumw2();

  TFile *fEtaQCDMIT = new TFile("/home/llr/cms/ivo/CMSSW_4_2_8_patch7/src/Bianchi/Limits/htautau/histograms/ElecTauStream_16Feb2012/eTau_mH120_inclusive__etaL1_MIT.root" ,"READ");  
  TH1F *hEtaQCDMIT           = (TH1F*)fEtaQCDMIT->Get("hQCD");
  hEtaQCDMIT->SetXTitle("#eta(e)");
  hEtaQCDMIT->SetFillColor(0);
  hEtaQCDMIT->SetMarkerColor(kRed);
  hEtaQCDMIT->SetMarkerStyle(kOpenCircle);
  hEtaQCDMIT->SetMarkerSize(1);
  hEtaQCDMIT->Sumw2();

  TFile *fEtaQCDWP80 = new TFile("/home/llr/cms/ivo/CMSSW_4_2_8_patch7/src/Bianchi/Limits/htautau/histograms/ElecTauStream_16Feb2012/eTau_mH120_inclusive__etaL1_CutWP80.root" ,"READ");  
  TH1F *hEtaQCDWP80           = (TH1F*)fEtaQCDWP80->Get("hQCD");
  hEtaQCDWP80->SetXTitle("#eta(e)");
  hEtaQCDWP80->SetFillColor(0);
  hEtaQCDWP80->SetMarkerColor(kGreen-9);
  hEtaQCDWP80->SetMarkerStyle(kOpenCircle);
  hEtaQCDWP80->SetMarkerSize(1);
  hEtaQCDWP80->Sumw2();




  //Pt histos for QCD estimated from data

 TFile *fQCDMVADaniele = new TFile("/home/llr/cms/ivo/CMSSW_4_2_8_patch7/src/Bianchi/Limits/htautau/histograms/ElecTauStream_16Feb2012/eTau_mH120_inclusive__tightestDanieleMVAWP_Daniele.root" ,"READ");  
  TH1F *hQCDMVADaniele           = (TH1F*)fQCDMVADaniele->Get("hQCD");
  hQCDMVADaniele->SetFillColor(0);

  TFile *fPtQCDDaniele = new TFile("/home/llr/cms/ivo/CMSSW_4_2_8_patch7/src/Bianchi/Limits/htautau/histograms/ElecTauStream_16Feb2012/eTau_mH120_inclusive__ptL1_Daniele.root" ,"READ");  
  TH1F *hPtQCDDaniele           = (TH1F*)fPtQCDDaniele->Get("hQCD");
  hPtQCDDaniele->SetXTitle("P_{T}(e)[GeV]");
  hPtQCDDaniele->SetFillColor(0);
  hPtQCDDaniele->SetMarkerColor(kBlue);
  hPtQCDDaniele->SetMarkerStyle(kOpenCircle);
  hPtQCDDaniele->SetMarkerSize(1);
  hPtQCDDaniele->Sumw2();

  TFile *fQCDMIT = new TFile("/home/llr/cms/ivo/CMSSW_4_2_8_patch7/src/Bianchi/Limits/htautau/histograms/ElecTauStream_16Feb2012/eTau_mH120_inclusive__ptL1_MIT.root" ,"READ");  
  TH1F *hPtQCDMIT           = (TH1F*)fQCDMIT->Get("hQCD");
  hPtQCDMIT->SetXTitle("P_{T}(e)[GeV]");
  hPtQCDMIT->SetFillColor(0);
  hPtQCDMIT->SetMarkerColor(kRed);
  hPtQCDMIT->SetMarkerStyle(kOpenCircle);
  hPtQCDMIT->SetMarkerSize(1);
  //hPtQCDMIT->Sumw2();

  TFile *fQCDWP80 = new TFile("/home/llr/cms/ivo/CMSSW_4_2_8_patch7/src/Bianchi/Limits/htautau/histograms/ElecTauStream_16Feb2012/eTau_mH120_inclusive__ptL1_CutWP80.root" ,"READ");  
  TH1F *hPtQCDWP80           = (TH1F*)fQCDWP80->Get("hQCD");
  hPtQCDWP80->SetXTitle("P_{T}(e)[GeV]");
  hPtQCDWP80->SetFillColor(0);
  hPtQCDWP80->SetMarkerColor(kGreen-8);
  hPtQCDWP80->SetMarkerStyle(kOpenCircle);
  hPtQCDWP80->SetMarkerSize(1);
  hPtQCDWP80->Sumw2();

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  //c1->Divide(2,1);
  
  TPad* pad1 = new TPad("pad1DEta","",0.05,0.22,0.96,0.97);
  TPad* pad2 = new TPad("pad2DEta","",0.05,0.02,0.96,0.20);

  pad1->SetFillColor(0);
  pad2->SetFillColor(0);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
//   pad1->SetLogy(logy_);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");

  TLegend* leg = new TLegend(0.55,0.42,0.85,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  leg->SetHeader("#splitline{CMS Preliminary}{ #sqrt{s}=7 TeV}");
  if(!ForBarrel)  leg->SetHeader("#splitline{CMS Preliminary}{7 TeV}}");


  TH1F* hDanieleMVA   = new TH1F( "hDanieleMVA" ,"" , nBins_ ,xMin_ , xMax_);
  hDanieleMVA->SetXTitle("MVA discriminator");
  hDanieleMVA->SetYTitle("Events");

  TH1F* hEff   = new TH1F( "hEff" ,"" , nBins_ ,xMin_ , xMax_);
  hEff->SetXTitle("Cut on MVA EleId");
  hEff->SetYTitle("Efficiency");

  //Daniele third WP
  TCut presel("(TMath::Abs(etaL1)<1.5 && sihih<0.01 && dEta<0.007 && dPhi<0.15 && HoE<0.12) || (TMath::Abs(etaL1)>1.5 && TMath::Abs(etaL1)<2.1 && sihih<0.03 && dEta<0.009 && dPhi<0.10 && HoE<0.10)");
  TCut lpt("nHits<=0 && antiConv>0.5 && (ptL1>20 && TMath::Abs(etaL1)<1.0  && tightestDanieleMVAWP>=0.013)||(ptL1>20 && TMath::Abs(etaL1)>1.0 && TMath::Abs(etaL1)<1.5  && tightestDanieleMVAWP>=0.0425)|| (ptL1>20 && TMath::Abs(etaL1)>1.5 && TMath::Abs(etaL1)<2.1  && tightestDanieleMVAWP>=0.0225)");
  TCut sbinDaniele = presel && lpt;
  //Daniele second WP
  //TCut sbinDaniele("(ptL1>20 && TMath::Abs(etaL1)>1.5 && TMath::Abs(etaL1)<2.1  && tightestDanieleMVAWP>=0.0225)||(ptL1>20 && TMath::Abs(etaL1)<1.5  && tightestDanieleMVAWP>=0.01625)");
  //Daniele first WP
//  TCut sbinDaniele("(ptL1>20 && ptL1<30 && TMath::Abs(etaL1)>1.5 && TMath::Abs(etaL1)<2.1  && tightestDanieleMVAWP>=0.0275)||(ptL1>20 && ptL1<30 && TMath::Abs(etaL1)<1.5  && tightestDanieleMVAWP>=0.035)||(ptL1>30 && TMath::Abs(etaL1)>1.5 && TMath::Abs(etaL1)<2.1 && tightestDanieleMVAWP>=0.0225)||(ptL1>30 && TMath::Abs(etaL1)<1.5  && tightestDanieleMVAWP>=0.035)");
  TCut sbinMIT("ptL1>20 && TMath::Abs(etaL1)<2.1 && tightestMVAWP>=1");
  TCut sbinWP80("ptL1>20 && TMath::Abs(etaL1)<2.1 && tightestCutBasedWP>=1");

  TH1F* hPtDanieleWP   = new TH1F( "hPtDanieleWP" ,"" , nBinsPt_ ,xMinPt_ , xMaxPt_);
  hPtDanieleWP->SetXTitle("P_{T}(e)[GeV]");
  hPtDanieleWP->SetMarkerColor(kBlue);
  hPtDanieleWP->SetMarkerStyle(kOpenCircle);
  hPtDanieleWP->SetMarkerSize(1);
  TH1F* hPtMITWP   = new TH1F( "hPtMITWP" ,"" , nBinsPt_ ,xMinPt_ , xMaxPt_);
  hPtMITWP->SetXTitle("P_{T}(e)[GeV]");
  hPtMITWP->SetMarkerColor(kRed);
  hPtMITWP->SetMarkerStyle(kOpenCircle);
  hPtMITWP->SetMarkerSize(1);
  TH1F* hPtWP80   = new TH1F( "hPtWP80" ,"" , nBinsPt_ ,xMinPt_ , xMaxPt_);
  hPtWP80->SetXTitle("P_{T}(e)[GeV]");
  hPtWP80->SetMarkerColor(kGreen-8);
  hPtWP80->SetMarkerStyle(kOpenCircle);
  hPtWP80->SetMarkerSize(1);
  TH1F* hEtaDanieleWP   = new TH1F( "hEtaDanieleWP" ,"" , nBinsEta_ ,xMinEta_ , xMaxEta_);
  hEtaDanieleWP->SetXTitle("#eta(e)");
  hEtaDanieleWP->SetMarkerColor(kBlue);
  hEtaDanieleWP->SetMarkerStyle(kOpenCircle);
  hEtaDanieleWP->SetMarkerSize(1);
  TH1F* hEtaMITWP   = new TH1F( "hEtaMITWP" ,"" , nBinsEta_ ,xMinEta_ , xMaxEta_);
  hEtaMITWP->SetXTitle("#eta(e)");
  hEtaMITWP->SetMarkerColor(kRed);
  hEtaMITWP->SetMarkerStyle(kOpenCircle);
  hEtaMITWP->SetMarkerSize(1);
  TH1F* hEtaWP80   = new TH1F( "hEtaWP80" ,"" , nBinsEta_ ,xMinEta_ , xMaxEta_);
  hEtaWP80->SetXTitle("#eta(e)");
  hEtaWP80->SetMarkerColor(kGreen-8);
  hEtaWP80->SetMarkerStyle(kOpenCircle);
  hEtaWP80->SetMarkerSize(1);

//   signal->Draw("tightestDanieleMVAWP>>hDanieleMVA");
//   signal->Draw("ptL1>>hPtDanieleWP","puWeight*(combRelIsoLeg1DBeta<0.2 && abs(etaL1)<1.5 && isElecLegMatched>0.5 && isTauLegMatched>0.5 && HLTx>0.5 && HLTmatch>0.5)"*sbinDaniele);
//   signal->Draw("ptL1>>hPtMITWP","puWeight*(combRelIsoLeg1DBeta<0.2 && abs(etaL1)<1.5 && isElecLegMatched>0.5 && isTauLegMatched>0.5 && HLTx>0.5 && HLTmatch>0.5)"*sbinMIT);
//   signal->Draw("ptL1>>hPtWP80","puWeight*(combRelIsoLeg1DBeta<0.2 && abs(etaL1)<1.5 && isElecLegMatched>0.5 && isTauLegMatched>0.5 && HLTx>0.5 && HLTmatch>0.5)"*sbinWP80);
//   signal->Draw("etaL1>>hEtaDanieleWP","puWeight*(combRelIsoLeg1DBeta<0.2 && abs(etaL1)<1.5 && isElecLegMatched>0.5 && isTauLegMatched>0.5 && HLTx>0.5 && HLTmatch>0.5)"*sbinDaniele);
//   signal->Draw("etaL1>>hEtaMITWP","puWeight*(combRelIsoLeg1DBeta<0.2 && abs(etaL1)<1.5 && isElecLegMatched>0.5 && isTauLegMatched>0.5 && HLTx>0.5 && HLTmatch>0.5)"*sbinMIT);
//   signal->Draw("etaL1>>hEtaWP80","puWeight*(combRelIsoLeg1DBeta<0.2 && abs(etaL1)<1.5 && isElecLegMatched>0.5 && isTauLegMatched>0.5 && HLTx>0.5 && HLTmatch>0.5)"*sbinWP80);
//   if (!ForBarrel){
//     signal->Draw("ptL1>>hEffvsPt",Form("puWeight*(combRelIsoLeg1DBeta<0.2 && abs(etaL1)>1.5&& ptL1>%0f && ptL1<%0f && isElecLegMatched>0.5 && isTauLegMatched>0.5 && HLTx>0.5 && HLTmatch>0.5 && HLTx>0.5 && HLTmatch>0.5)", PtCutMin_, PtCutMax_));
//     signal->Draw("ptL1>>hEffvsPtP",Form("puWeight*(combRelIsoLeg1DBeta<0.2 && tightestDanieleMVAWP>%0f && abs(etaL1)>1.5 && ptL1>%0f && ptL1<%0f && isElecLegMatched>0.5 && isTauLegMatched>0.5 && HLTx>0.5 && HLTmatch>0.5)",MVACut_, PtCutMin_, PtCutMax_));
//   }


  QCDBCtoE2080->Draw("tightestDanieleMVAWP>>hDanieleMVA");
  QCDBCtoE2080->Draw("ptL1>>hPtDanieleWP","puWeight*(combRelIsoLeg1DBeta<0.3  && HLTx>0.5 && HLTmatch>0.5)"*sbinDaniele);
  QCDBCtoE2080->Draw("ptL1>>hPtMITWP","puWeight*(combRelIsoLeg1DBeta<0.3 && HLTx>0.5 && HLTmatch>0.5)"*sbinMIT);
  QCDBCtoE2080->Draw("ptL1>>hPtWP80","puWeight*(combRelIsoLeg1DBeta<0.3  && HLTx>0.5 && HLTmatch>0.5)"*sbinWP80);
  QCDBCtoE2080->Draw("etaL1>>hEtaDanieleWP","puWeight*(combRelIsoLeg1DBeta<0.2  && HLTx>0.5 && HLTmatch>0.5)"*sbinDaniele);
  QCDBCtoE2080->Draw("etaL1>>hEtaMITWP","puWeight*(combRelIsoLeg1DBeta<0.2  && HLTx>0.5 && HLTmatch>0.5)"*sbinMIT);
  QCDBCtoE2080->Draw("etaL1>>hEtaWP80","puWeight*(combRelIsoLeg1DBeta<0.2 && HLTx>0.5 && HLTmatch>0.5)"*sbinWP80);

  cout<<"For QCD BCtoE :"<<endl;
  cout<<"After Daniele's MVA :"<<hPtDanieleWP->GetEntries()<<endl;
  cout<<"After MIT's MVA :"<<hPtMITWP->GetEntries()<<endl;
  cout<<"After CutWP80 :"<<hPtWP80->GetEntries()<<endl;


  hPtDanieleWP->Sumw2();
  hPtMITWP->Sumw2();
  hPtWP80->Sumw2();
  hEtaDanieleWP->Sumw2();
  hEtaMITWP->Sumw2();
  hEtaWP80->Sumw2();

  //hDanieleMVA->Draw();
  //hQCDMVADaniele->Draw();
  //leg->AddEntry(hDanieleMVA,"DY");
  //leg->AddEntry(hQCDMVADaniele,"QCD");
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

  


//   // Eta  spectra ot electrons passing Daniele and MIT ID for DY

//   hEtaMITWP->Draw();
//   hEtaDanieleWP->Draw("same");
//   hEtaWP80->Draw("same");
//   leg->AddEntry(hEtaMITWP,"MIT DY");
//   leg->AddEntry(hEtaDanieleWP,"Daniele DY");
//   leg->AddEntry(hEtaWP80,"WP80 DY");
//   leg->Draw();

//   pad2->cd();
//   gStyle->SetOptStat(0);
//   gStyle->SetTitleFillColor(0);
//   gStyle->SetCanvasBorderMode(0);
//   gStyle->SetCanvasColor(0);
//   gStyle->SetPadBorderMode(0);
//   gStyle->SetPadColor(0);
//   gStyle->SetTitleFillColor(0);
//   gStyle->SetTitleBorderSize(0);
//   gStyle->SetTitleH(0.07);
//   gStyle->SetTitleFontSize(0.1);
//   gStyle->SetTitleStyle(0);
//   gStyle->SetTitleOffset(1.3,"y");
  
//   TH1F* hRatio = (TH1F*)hEtaDanieleWP->Clone("hRatio");
//   //hRatio->Reset();
//   hRatio->SetXTitle("");
//   hRatio->SetYTitle("#frac{Daniele}{MIT}");
  
//   hRatio->SetMarkerStyle(kFullCircle);
//   hRatio->SetMarkerSize(0.8);
//   hRatio->SetLabelSize(0.12,"X");
//   hRatio->SetLabelSize(0.10,"Y");
//   hRatio->SetTitleSize(0.12,"Y");
//   hRatio->SetTitleOffset(0.36,"Y");

//   hRatio->Divide(hEtaMITWP);
//   hRatio->Draw("P");

//   TF1* line = new TF1("line","1",hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());
//   //line->SetLineStyle(3);
//   line->SetLineWidth(1.5);
//   line->SetLineColor(kBlack);
//   line->Draw("SAME");









  
 //  // Pt  spectra ot electrons passing Daniele and MIT ID for DY
//   hPtMITWP->Draw();
//   hPtDanieleWP->Draw("same");
//   hPtWP80->Draw("same");
//   leg->AddEntry(hPtMITWP,"MIT DY");
//   leg->AddEntry(hPtDanieleWP,"Daniele DY");
//   leg->AddEntry(hPtWP80,"WP80 DY");
//   leg->Draw();  
  
//   pad2->cd();
//   gStyle->SetOptStat(0);
//   gStyle->SetTitleFillColor(0);
//   gStyle->SetCanvasBorderMode(0);
//   gStyle->SetCanvasColor(0);
//   gStyle->SetPadBorderMode(0);
//   gStyle->SetPadColor(0);
//   gStyle->SetTitleFillColor(0);
//   gStyle->SetTitleBorderSize(0);
//   gStyle->SetTitleH(0.07);
//   gStyle->SetTitleFontSize(0.1);
//   gStyle->SetTitleStyle(0);
//   gStyle->SetTitleOffset(1.3,"y");
  
//   TH1F* hRatio = (TH1F*)hPtDanieleWP->Clone("hRatio");
//   //hRatio->Reset();
//   hRatio->SetXTitle("");
//   hRatio->SetYTitle("#frac{Daniele}{MIT}");
  
//   hRatio->SetMarkerStyle(kFullCircle);
//   hRatio->SetMarkerSize(0.8);
//   hRatio->SetLabelSize(0.12,"X");
//   hRatio->SetLabelSize(0.10,"Y");
//   hRatio->SetTitleSize(0.12,"Y");
//   hRatio->SetTitleOffset(0.36,"Y");
  
//   hRatio->Divide(hPtMITWP);
//   hRatio->Draw("P");
    
//   TF1* line = new TF1("line","1",hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());
//   //line->SetLineStyle(3);
//   line->SetLineWidth(1.5);
//   line->SetLineColor(kBlack);
//   line->Draw("SAME");








//   // Eta spectra ot electrons passing Daniele and MIT ID for QCD
//     hEtaQCDDaniele->Draw("P");
//     hEtaQCDWP80->Draw("PSAME");
//     hEtaQCDMIT->Draw("PSAME");
//     leg->AddEntry(hEtaQCDMIT,"QCD MIT", "P");
//     leg->AddEntry(hEtaQCDDaniele,"QCD Daniele","P");
//     leg->AddEntry(hEtaQCDWP80,"QCD WP80","P");
//     leg->Draw();
    
//     pad2->cd();
//     gStyle->SetOptStat(0);
//     gStyle->SetTitleFillColor(0);
//     gStyle->SetCanvasBorderMode(0);
//     gStyle->SetCanvasColor(0);
//     gStyle->SetPadBorderMode(0);
//     gStyle->SetPadColor(0);
//     gStyle->SetTitleFillColor(0);
//     gStyle->SetTitleBorderSize(0);
//     gStyle->SetTitleH(0.07);
//     gStyle->SetTitleFontSize(0.1);
//     gStyle->SetTitleStyle(0);
//     gStyle->SetTitleOffset(1.3,"y");
    
//   TH1F* hRatio = (TH1F*)hEtaQCDDaniele->Clone("hRatio");
//   //hRatio->Reset();
//   hRatio->SetXTitle("");
//   hRatio->SetYTitle("#frac{Daniele}{MIT}");
//   hRatio->SetAxisRange(0.0,2.5,"Y");

//   hRatio->SetMarkerStyle(kFullCircle);
//   hRatio->SetMarkerSize(0.8);
//   hRatio->SetLabelSize(0.12,"X");
//   hRatio->SetLabelSize(0.10,"Y");
//   hRatio->SetTitleSize(0.12,"Y");
//   hRatio->SetTitleOffset(0.36,"Y");

//   hRatio->Divide(hEtaQCDMIT);
//   hRatio->Draw("P");

//   TF1* line = new TF1("line","1",hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());
//   //line->SetLineStyle(3);
//   line->SetLineWidth(1.5);
//   line->SetLineColor(kBlack);
//   line->Draw("SAME");








//   // Pt spectra ot electrons passing Daniele and MIT ID for QCD
//     hPtQCDMIT->Draw("P");
//     hPtQCDDaniele->Draw("PSAME");
//     hPtQCDWP80->Draw("PSAME");
//     leg->AddEntry(hPtQCDMIT,"QCD MIT", "P");
//     leg->AddEntry(hPtQCDDaniele,"QCD Daniele","P");
//     leg->AddEntry(hPtQCDWP80,"QCD WP80","P");
//     leg->Draw();
    
//     pad2->cd();
//     gStyle->SetOptStat(0);
//     gStyle->SetTitleFillColor(0);
//     gStyle->SetCanvasBorderMode(0);
//     gStyle->SetCanvasColor(0);
//     gStyle->SetPadBorderMode(0);
//     gStyle->SetPadColor(0);
//     gStyle->SetTitleFillColor(0);
//     gStyle->SetTitleBorderSize(0);
//     gStyle->SetTitleH(0.07);
//     gStyle->SetTitleFontSize(0.1);
//     gStyle->SetTitleStyle(0);
//     gStyle->SetTitleOffset(1.3,"y");
    
//   TH1F* hRatio = (TH1F*)hPtQCDDaniele->Clone("hRatio");
//   //hRatio->Reset();
//   hRatio->SetXTitle("");
//   hRatio->SetYTitle("#frac{Daniele}{MIT}");
//   hRatio->SetAxisRange(0.0,2.5,"Y");

//   hRatio->SetMarkerStyle(kFullCircle);
//   hRatio->SetMarkerSize(0.8);
//   hRatio->SetLabelSize(0.12,"X");
//   hRatio->SetLabelSize(0.10,"Y");
//   hRatio->SetTitleSize(0.12,"Y");
//   hRatio->SetTitleOffset(0.36,"Y");

//   hRatio->Divide(hPtQCDMIT);
//   hRatio->Draw("P");

//   TF1* line = new TF1("line","1",hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());
//   //line->SetLineStyle(3);
//   line->SetLineWidth(1.5);
//   line->SetLineColor(kBlack);
//   line->Draw("SAME");




}

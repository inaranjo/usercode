#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TGraphAsymmErrors.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include "TCut.h"

#include <string>
#include <map>
#include <iostream>
#include <iomanip>

#define DEBUG false

void plotEfficienciesVsPt(string type = "Eff"
			  )
{

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

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

  TFile *fSignal = new TFile("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_makeWP-iter1_Pythia.root","READ"); 
//   TFile *fSignal = new TFile("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_makeWP-PtEleRelaxed-iter1_Pythia.root","READ"); 
//   TFile *fSignal = new TFile("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_makeWP_Pythia.root","READ"); 
  if(fSignal->IsZombie()){
    std::cout<<"No such file!"<<std::endl;
    return;
  }

  TTree*inTree = (TTree*)fSignal->Get("tree");

  TH1F* hMVASig   = new TH1F( "hMVASig","" ,100,-1,1);
  TH1F* hMVABkg = new TH1F( "hMVABkg","" ,100,-1,1);

  int numPVMin = 0;
  int numPVMax = 50;
  float EtaMin = 0;
  float EtaMax = 3.0;
  float PhiMin = -3;
  float PhiMax = 3;

  float nTot = 0;
  float nPassSig_MVA = 0;
  float nPassSig_MVA2_WP75 = 0;
  float nPassSig_MVA2_WP85 = 0;
  float nPassSig_MVA2_WP95 = 0;
  float nPassBkg_MVA = 0;
  float nPassBkg_MVA2_WP75 = 0;
  float nPassBkg_MVA2_WP85 = 0;
  float nPassBkg_MVA2_WP95 = 0;

  int nBinsPt = 5;
  float Pt [nBinsPt];
  Pt[0] =30;
  Pt[1] =50;
  Pt[2] =70;
  Pt[3] =90;
  Pt[4] =110;

  float EffPt_MVA [nBinsPt-1];
  float EffPt_MVA2_WP75 [nBinsPt-1];
  float EffPt_MVA2_WP85 [nBinsPt-1];
  float EffPt_MVA2_WP95 [nBinsPt-1];
  float FakePt_MVA [nBinsPt-1];
  float FakePt_MVA2_WP75 [nBinsPt-1];
  float FakePt_MVA2_WP85 [nBinsPt-1];
  float FakePt_MVA2_WP95 [nBinsPt-1];

  for(int k=0;k<nBinsPt-1;k++){

    TCut PUSelection(Form("NumPV>%i && NumPV<%i",numPVMin,numPVMax));
    TCut EtaSelection (Form("Tau_AbsEta>%0f && Tau_AbsEta<%0f",EtaMin,EtaMax));
    TCut PhiSelection (Form("Tau_Phi>%0f && Tau_Phi<%0f",PhiMin,PhiMax));
    TCut PtSelection (Form("Tau_Pt>%0f && Tau_Pt<%0f",Pt[k],Pt[k+1]));
    TCut Preselection = "Tau_AntiELoose";
    TCut Selection = PUSelection && PtSelection && EtaSelection && PhiSelection && Preselection;
        
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch"*Selection);
    nTot = hMVASig->GetEntries();  
    hMVASig->Reset();  
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA"*Selection);
    inTree->Draw("Tau_AntiEMVA>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA"*Selection);
    nPassSig_MVA = hMVASig->GetEntries();
    nPassBkg_MVA = hMVABkg->GetEntries();
    hMVASig->Reset();  
    hMVABkg->Reset(); 
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA2_WP75"*Selection);
    inTree->Draw("Tau_AntiEMVA>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA2_WP75"*Selection);
    nPassSig_MVA2_WP75 = hMVASig->GetEntries();
    nPassBkg_MVA2_WP75 = hMVABkg->GetEntries();
    hMVASig->Reset();  
    hMVABkg->Reset();
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA2_WP85"*Selection);
    inTree->Draw("Tau_AntiEMVA>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA2_WP85"*Selection);
    nPassSig_MVA2_WP85 = hMVASig->GetEntries();
    nPassBkg_MVA2_WP85 = hMVABkg->GetEntries();
    hMVASig->Reset();  
    hMVABkg->Reset();
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA2_WP95"*Selection);
    inTree->Draw("Tau_AntiEMVA>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA2_WP95"*Selection);
    nPassSig_MVA2_WP95 = hMVASig->GetEntries();
    nPassBkg_MVA2_WP95 = hMVABkg->GetEntries();
    hMVASig->Reset();  
    hMVABkg->Reset();
  
    EffPt_MVA[k] = nPassSig_MVA/nTot;
    FakePt_MVA[k] = nPassBkg_MVA/nTot;
    EffPt_MVA2_WP75[k] = nPassSig_MVA2_WP75/nTot;
    FakePt_MVA2_WP75[k] = nPassBkg_MVA2_WP75/nTot;    
    EffPt_MVA2_WP85[k] = nPassSig_MVA2_WP85/nTot;
    FakePt_MVA2_WP85[k] = nPassBkg_MVA2_WP85/nTot;
    EffPt_MVA2_WP95[k] = nPassSig_MVA2_WP95/nTot;
    FakePt_MVA2_WP95[k] = nPassBkg_MVA2_WP95/nTot;
    
    if(DEBUG){
    std::cout<<"For Pt between :"<<Pt[k]<<" and "<<Pt[k+1]<<std::endl;
//     std::cout<<" Total :"<<nTot<<std::endl;
//     std::cout<<" passing signal MVA:"<<nPassSig_MVA<<std::endl;
//     std::cout<<" passing bkg MVA:"<<nPassBkg_MVA<<std::endl;
    std::cout<<" signal eff MVA:"<<EffPt_MVA[k]<<std::endl;
    std::cout<<" bkg eff MVA:"<<FakePt_MVA[k]<<std::endl;
//     std::cout<<std::endl;
//     std::cout<<" passing signal MVA WP75:"<<nPassSig_MVA2_WP75<<std::endl;
//     std::cout<<" passing bkg MVA WP75:"<<nPassBkg_MVA2_WP75<<std::endl;
    std::cout<<" signal eff MVA WP75:"<<EffPt_MVA2_WP75[k]<<std::endl;
    std::cout<<" bkg eff MVA WP75:"<<FakePt_MVA2_WP75[k]<<std::endl;
//     std::cout<<std::endl;
//     std::cout<<" passing signal MVA WP85:"<<nPassSig_MVA2_WP85<<std::endl;
//     std::cout<<" passing bkg MVA WP85:"<<nPassBkg_MVA2_WP85<<std::endl;
    std::cout<<" signal eff MVA WP85:"<<EffPt_MVA2_WP85[k]<<std::endl;
    std::cout<<" bkg eff MVA WP85:"<<FakePt_MVA2_WP85[k]<<std::endl;
//     std::cout<<std::endl;
//     std::cout<<" passing signal MVA WP95:"<<nPassSig_MVA2_WP95<<std::endl;
//     std::cout<<" passing bkg MVA WP95:"<<nPassBkg_MVA2_WP95<<std::endl;
    std::cout<<" signal eff MVA WP95:"<<EffPt_MVA2_WP95[k]<<std::endl;
    std::cout<<" bkg eff MVA WP95:"<<FakePt_MVA2_WP95[k]<<std::endl;
//     std::cout<<std::endl; 
    }
    Pt[k] = (Pt[k]+Pt[k+1])/2;

  }
  
  TGraph* EffVsPt_MVA = new TGraph (nBinsPt-1,Pt,EffPt_MVA);
  EffVsPt_MVA->SetMarkerColor(kBlue); 
  EffVsPt_MVA->SetMarkerStyle(21); 
  EffVsPt_MVA->SetMarkerSize(1); 
  EffVsPt_MVA->SetLineColor(kBlue);
  EffVsPt_MVA->SetFillColor(0);
  EffVsPt_MVA->SetLineWidth(1);
  EffVsPt_MVA->SetTitle("");

  TGraph* EffVsPt_MVA2_WP75 = new TGraph (nBinsPt-1,Pt,EffPt_MVA2_WP75);
  EffVsPt_MVA2_WP75->SetMarkerColor(kRed); 
  EffVsPt_MVA2_WP75->SetMarkerStyle(21); 
  EffVsPt_MVA2_WP75->SetMarkerSize(1); 
  EffVsPt_MVA2_WP75->SetLineColor(kRed);
  EffVsPt_MVA2_WP75->SetFillColor(0);
  EffVsPt_MVA2_WP75->SetLineWidth(1);
  EffVsPt_MVA2_WP75->SetTitle("WP75");

  TGraph* EffVsPt_MVA2_WP85 = new TGraph (nBinsPt-1,Pt,EffPt_MVA2_WP85);
  EffVsPt_MVA2_WP85->SetMarkerColor(kGreen); 
  EffVsPt_MVA2_WP85->SetMarkerStyle(21); 
  EffVsPt_MVA2_WP85->SetMarkerSize(1); 
  EffVsPt_MVA2_WP85->SetLineColor(kGreen);
  EffVsPt_MVA2_WP85->SetFillColor(0);
  EffVsPt_MVA2_WP85->SetLineWidth(1);
  EffVsPt_MVA2_WP85->SetTitle("WP85");

  TGraph* EffVsPt_MVA2_WP95 = new TGraph (nBinsPt-1,Pt,EffPt_MVA2_WP95);
  EffVsPt_MVA2_WP95->SetMarkerColor(kMagenta); 
  EffVsPt_MVA2_WP95->SetMarkerStyle(21); 
  EffVsPt_MVA2_WP95->SetMarkerSize(1); 
  EffVsPt_MVA2_WP95->SetLineColor(kMagenta);
  EffVsPt_MVA2_WP95->SetFillColor(0);
  EffVsPt_MVA2_WP95->SetLineWidth(1);
  EffVsPt_MVA2_WP95->SetTitle("WP95");

  TMultiGraph *mg_EffVsPt = new TMultiGraph();
  mg_EffVsPt->Add(EffVsPt_MVA);
  mg_EffVsPt->Add(EffVsPt_MVA2_WP75);
  mg_EffVsPt->Add(EffVsPt_MVA2_WP85);
  mg_EffVsPt->Add(EffVsPt_MVA2_WP95);

  TLegend* leg_EffVsPt = new TLegend(0.6,0.75,0.8,0.88,NULL,"brNDC");
  leg_EffVsPt->SetFillStyle(0);
  leg_EffVsPt->SetBorderSize(0);
  leg_EffVsPt->SetFillColor(10);
  leg_EffVsPt->SetTextSize(0.03);

  leg_EffVsPt->AddEntry(EffVsPt_MVA,"old MVA");
  leg_EffVsPt->AddEntry(EffVsPt_MVA2_WP75,"new MVA WP75");
  leg_EffVsPt->AddEntry(EffVsPt_MVA2_WP85,"new MVA WP85");
  leg_EffVsPt->AddEntry(EffVsPt_MVA2_WP95,"new MVA WP95");

  ///FakeRates
  TGraph* FakeVsPt_MVA = new TGraph (nBinsPt-1,Pt,FakePt_MVA);
  FakeVsPt_MVA->SetMarkerColor(kBlue); 
  FakeVsPt_MVA->SetMarkerStyle(21); 
  FakeVsPt_MVA->SetMarkerSize(1); 
  FakeVsPt_MVA->SetLineColor(kBlue);
  FakeVsPt_MVA->SetFillColor(0);
  FakeVsPt_MVA->SetLineWidth(1);
  FakeVsPt_MVA->SetTitle("");

  TGraph* FakeVsPt_MVA2_WP75 = new TGraph (nBinsPt-1,Pt,FakePt_MVA2_WP75 );
  FakeVsPt_MVA2_WP75->SetMarkerColor(kRed); 
  FakeVsPt_MVA2_WP75->SetMarkerStyle(21); 
  FakeVsPt_MVA2_WP75->SetMarkerSize(1); 
  FakeVsPt_MVA2_WP75->SetLineColor(kRed);
  FakeVsPt_MVA2_WP75->SetFillColor(0);
  FakeVsPt_MVA2_WP75->SetLineWidth(1);
  FakeVsPt_MVA2_WP75->SetTitle("");

  TGraph* FakeVsPt_MVA2_WP85 = new TGraph (nBinsPt-1,Pt,FakePt_MVA2_WP85 );
  FakeVsPt_MVA2_WP85->SetMarkerColor(kGreen); 
  FakeVsPt_MVA2_WP85->SetMarkerStyle(21); 
  FakeVsPt_MVA2_WP85->SetMarkerSize(1); 
  FakeVsPt_MVA2_WP85->SetLineColor(kGreen);
  FakeVsPt_MVA2_WP85->SetFillColor(0);
  FakeVsPt_MVA2_WP85->SetLineWidth(1);
  FakeVsPt_MVA2_WP85->SetTitle("");

  TGraph* FakeVsPt_MVA2_WP95 = new TGraph (nBinsPt-1,Pt,FakePt_MVA2_WP95 );
  FakeVsPt_MVA2_WP95->SetMarkerColor(kMagenta); 
  FakeVsPt_MVA2_WP95->SetMarkerStyle(21); 
  FakeVsPt_MVA2_WP95->SetMarkerSize(1); 
  FakeVsPt_MVA2_WP95->SetLineColor(kMagenta);
  FakeVsPt_MVA2_WP95->SetFillColor(0);
  FakeVsPt_MVA2_WP95->SetLineWidth(1);
  FakeVsPt_MVA2_WP95->SetTitle("");

  TMultiGraph *mg_FakeVsPt = new TMultiGraph();
  mg_FakeVsPt->Add(FakeVsPt_MVA);
  mg_FakeVsPt->Add(FakeVsPt_MVA2_WP75);
  mg_FakeVsPt->Add(FakeVsPt_MVA2_WP85);
  mg_FakeVsPt->Add(FakeVsPt_MVA2_WP95);

  TLegend* leg_FakeVsPt = new TLegend(0.6,0.75,0.8,0.88,NULL,"brNDC");
  leg_FakeVsPt->SetFillStyle(0);
  leg_FakeVsPt->SetBorderSize(0);
  leg_FakeVsPt->SetFillColor(10);
  leg_FakeVsPt->SetTextSize(0.03);

  leg_FakeVsPt->AddEntry(FakeVsPt_MVA,"old MVA");
  leg_FakeVsPt->AddEntry(FakeVsPt_MVA2_WP75,"new MVA WP75");
  leg_FakeVsPt->AddEntry(FakeVsPt_MVA2_WP85,"new MVA WP85");
  leg_FakeVsPt->AddEntry(FakeVsPt_MVA2_WP95,"new MVA WP95");

  c1->Clear();

  if(type=="Eff" ){
    mg_EffVsPt->Draw("APL");
    mg_EffVsPt->GetXaxis()->SetTitle("#tau P_{T}(GeV)");
    mg_EffVsPt->GetXaxis()->SetRangeUser(0,130);
    mg_EffVsPt->GetYaxis()->SetTitle("Efficiency");
    mg_EffVsPt->SetMinimum(0.3);
    mg_EffVsPt->SetMaximum(1.);
    leg_EffVsPt->Draw("same");
  }
  if(type=="FakeRate"){
    mg_FakeVsPt->Draw("APL");
    mg_FakeVsPt->GetXaxis()->SetTitle("#tau P_{T}(GeV)");
    mg_FakeVsPt->GetXaxis()->SetRangeUser(0,130);
    mg_FakeVsPt->GetYaxis()->SetTitle("Fake Rate");
    mg_FakeVsPt->SetMinimum(0.);
    mg_FakeVsPt->SetMaximum(1.0);
    leg_FakeVsPt->Draw("same");
  }
  string outputName = Form("plots/plotEfficiencies_%sVsPt",type.data());
  c1->Print(std::string(outputName).append(".png").data());
  c1->Print(std::string(outputName).append(".pdf").data());
  
  fSignal->Close();

}



void plotEfficienciesVsEta(string type = "Eff"
			  )
{

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

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

  TFile *fSignal = new TFile("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_makeWP-iter1_Pythia.root","READ"); 
//   TFile *fSignal = new TFile("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_makeWP-PtEleRelaxed-iter1_Pythia.root","READ"); 
//   TFile *fSignal = new TFile("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_makeWP_Pythia.root","READ"); 
  if(fSignal->IsZombie()){
    std::cout<<"No such file!"<<std::endl;
    return;
  }

  TTree*inTree = (TTree*)fSignal->Get("tree");

  TH1F* hMVASig   = new TH1F( "hMVASig","" ,100,-1,1);
  TH1F* hMVABkg = new TH1F( "hMVABkg","" ,100,-1,1);

  float PtMin = 20;
  float PtMax = 100;
  int numPVMin = 0;
  int numPVMax = 50;
  float PhiMin = -3;
  float PhiMax = 3;

  float nTot = 0;
  float nPassSig_MVA = 0;
  float nPassSig_MVA2_WP75 = 0;
  float nPassSig_MVA2_WP85 = 0;
  float nPassSig_MVA2_WP95 = 0;
  float nPassBkg_MVA = 0;
  float nPassBkg_MVA2_WP75 = 0;
  float nPassBkg_MVA2_WP85 = 0;
  float nPassBkg_MVA2_WP95 = 0;

  int nBinsEta = 4;
  float Eta [nBinsEta];
  Eta[0] = 0.;
  Eta[1] = 0.8;
  Eta[2] = 1.5;
  Eta[3] = 2.3;

  float EffEta_MVA [nBinsEta-1];
  float EffEta_MVA2_WP75 [nBinsEta-1];
  float EffEta_MVA2_WP85 [nBinsEta-1];
  float EffEta_MVA2_WP95 [nBinsEta-1];
  float FakeEta_MVA [nBinsEta-1];
  float FakeEta_MVA2_WP75 [nBinsEta-1];
  float FakeEta_MVA2_WP85 [nBinsEta-1];
  float FakeEta_MVA2_WP95 [nBinsEta-1];

  for(int k=0;k<nBinsEta-1;k++){

    TCut PUSelection(Form("NumPV>%i && NumPV<%i",numPVMin,numPVMax));
    TCut EtaSelection (Form("Tau_AbsEta>%0f && Tau_AbsEta<%0f",Eta[k],Eta[k+1]));
    TCut PhiSelection (Form("Tau_Phi>%0f && Tau_Phi<%0f",PhiMin,PhiMax));
    TCut PtSelection (Form("Tau_Pt>%0f && Tau_Pt<%0f",PtMin,PtMax));
    TCut Preselection = "Tau_AntiELoose";
    TCut Selection = PUSelection && PtSelection && EtaSelection && PhiSelection && Preselection;
        
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch"*Selection);
    inTree->Draw("Tau_AntiEMVA>>hMVABkg","Tau_GenEleMatch"*Selection);
    nTot = hMVASig->GetEntries();    
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA"*Selection);
    inTree->Draw("Tau_AntiEMVA>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA"*Selection);
    nPassSig_MVA = hMVASig->GetEntries();
    nPassBkg_MVA = hMVABkg->GetEntries();
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA2_WP75"*Selection);
    inTree->Draw("Tau_AntiEMVA>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA2_WP75"*Selection);
    nPassSig_MVA2_WP75 = hMVASig->GetEntries();
    nPassBkg_MVA2_WP75 = hMVABkg->GetEntries();
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA2_WP85"*Selection);
    inTree->Draw("Tau_AntiEMVA>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA2_WP85"*Selection);
    nPassSig_MVA2_WP85 = hMVASig->GetEntries();
    nPassBkg_MVA2_WP85 = hMVABkg->GetEntries();
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA2_WP95"*Selection);
    inTree->Draw("Tau_AntiEMVA>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA2_WP95"*Selection);
    nPassSig_MVA2_WP95 = hMVASig->GetEntries();
    nPassBkg_MVA2_WP95 = hMVABkg->GetEntries();
  
    EffEta_MVA[k] = nPassSig_MVA/nTot;
    FakeEta_MVA[k] = nPassBkg_MVA/nTot;
    EffEta_MVA2_WP75[k] = nPassSig_MVA2_WP75/nTot;
    FakeEta_MVA2_WP75[k] = nPassBkg_MVA2_WP75/nTot;    
    EffEta_MVA2_WP85[k] = nPassSig_MVA2_WP85/nTot;
    FakeEta_MVA2_WP85[k] = nPassBkg_MVA2_WP85/nTot;
    EffEta_MVA2_WP95[k] = nPassSig_MVA2_WP95/nTot;
    FakeEta_MVA2_WP95[k] = nPassBkg_MVA2_WP95/nTot;
    
    if(DEBUG){
    std::cout<<"For Eta between :"<<Eta[k]<<" and "<<Eta[k+1]<<std::endl;
//     std::cout<<" Total :"<<nTot<<std::endl;
//     std::cout<<" passing signal MVA:"<<nPassSig_MVA<<std::endl;
//     std::cout<<" passing bkg MVA:"<<nPassBkg_MVA<<std::endl;
    std::cout<<" signal eff MVA:"<<EffEta_MVA[k]<<std::endl;
    std::cout<<" bkg eff MVA:"<<FakeEta_MVA[k]<<std::endl;
//     std::cout<<std::endl;
//     std::cout<<" passing signal MVA WP75:"<<nPassSig_MVA2_WP75<<std::endl;
//     std::cout<<" passing bkg MVA WP75:"<<nPassBkg_MVA2_WP75<<std::endl;
    std::cout<<" signal eff MVA WP75:"<<EffEta_MVA2_WP75[k]<<std::endl;
    std::cout<<" bkg eff MVA WP75:"<<FakeEta_MVA2_WP75[k]<<std::endl;
//     std::cout<<std::endl;
//     std::cout<<" passing signal MVA WP85:"<<nPassSig_MVA2_WP85<<std::endl;
//     std::cout<<" passing bkg MVA WP85:"<<nPassBkg_MVA2_WP85<<std::endl;
    std::cout<<" signal eff MVA WP85:"<<EffEta_MVA2_WP85[k]<<std::endl;
    std::cout<<" bkg eff MVA WP85:"<<FakeEta_MVA2_WP85[k]<<std::endl;
//     std::cout<<std::endl;
//     std::cout<<" passing signal MVA WP95:"<<nPassSig_MVA2_WP95<<std::endl;
//     std::cout<<" passing bkg MVA WP95:"<<nPassBkg_MVA2_WP95<<std::endl;
    std::cout<<" signal eff MVA WP95:"<<EffEta_MVA2_WP95[k]<<std::endl;
    std::cout<<" bkg eff MVA WP95:"<<FakeEta_MVA2_WP95[k]<<std::endl;
//     std::cout<<std::endl; 
    }
    Eta[k] = (Eta[k]+Eta[k+1])/2;
  }


  TGraph* EffVsEta_MVA = new TGraph (nBinsEta-1,Eta,EffEta_MVA);
  EffVsEta_MVA->SetMarkerColor(kBlue); 
  EffVsEta_MVA->SetMarkerStyle(21); 
  EffVsEta_MVA->SetMarkerSize(1); 
  EffVsEta_MVA->SetLineColor(kBlue);
  EffVsEta_MVA->SetFillColor(0);
  EffVsEta_MVA->SetLineWidth(1);
  EffVsEta_MVA->SetTitle("");

  TGraph* EffVsEta_MVA2_WP75 = new TGraph (nBinsEta-1,Eta,EffEta_MVA2_WP75);
  EffVsEta_MVA2_WP75->SetMarkerColor(kRed); 
  EffVsEta_MVA2_WP75->SetMarkerStyle(21); 
  EffVsEta_MVA2_WP75->SetMarkerSize(1); 
  EffVsEta_MVA2_WP75->SetLineColor(kRed);
  EffVsEta_MVA2_WP75->SetFillColor(0);
  EffVsEta_MVA2_WP75->SetLineWidth(1);
  EffVsEta_MVA2_WP75->SetTitle("WP75");

  TGraph* EffVsEta_MVA2_WP85 = new TGraph (nBinsEta-1,Eta,EffEta_MVA2_WP85);
  EffVsEta_MVA2_WP85->SetMarkerColor(kGreen); 
  EffVsEta_MVA2_WP85->SetMarkerStyle(21); 
  EffVsEta_MVA2_WP85->SetMarkerSize(1); 
  EffVsEta_MVA2_WP85->SetLineColor(kGreen);
  EffVsEta_MVA2_WP85->SetFillColor(0);
  EffVsEta_MVA2_WP85->SetLineWidth(1);
  EffVsEta_MVA2_WP85->SetTitle("WP85");

  TGraph* EffVsEta_MVA2_WP95 = new TGraph (nBinsEta-1,Eta,EffEta_MVA2_WP95);
  EffVsEta_MVA2_WP95->SetMarkerColor(kMagenta); 
  EffVsEta_MVA2_WP95->SetMarkerStyle(21); 
  EffVsEta_MVA2_WP95->SetMarkerSize(1); 
  EffVsEta_MVA2_WP95->SetLineColor(kMagenta);
  EffVsEta_MVA2_WP95->SetFillColor(0);
  EffVsEta_MVA2_WP95->SetLineWidth(1);
  EffVsEta_MVA2_WP95->SetTitle("WP95");

  TMultiGraph *mg_EffVsEta = new TMultiGraph();
  mg_EffVsEta->Add(EffVsEta_MVA);
  mg_EffVsEta->Add(EffVsEta_MVA2_WP75);
  mg_EffVsEta->Add(EffVsEta_MVA2_WP85);
  mg_EffVsEta->Add(EffVsEta_MVA2_WP95);

  TLegend* leg_EffVsEta = new TLegend(0.6,0.75,0.8,0.88,NULL,"brNDC");
  leg_EffVsEta->SetFillStyle(0);
  leg_EffVsEta->SetBorderSize(0);
  leg_EffVsEta->SetFillColor(10);
  leg_EffVsEta->SetTextSize(0.03);

  leg_EffVsEta->AddEntry(EffVsEta_MVA,"old MVA");
  leg_EffVsEta->AddEntry(EffVsEta_MVA2_WP75,"new MVA WP75");
  leg_EffVsEta->AddEntry(EffVsEta_MVA2_WP85,"new MVA WP85");
  leg_EffVsEta->AddEntry(EffVsEta_MVA2_WP95,"new MVA WP95");

  ///FakeRates
  TGraph* FakeVsEta_MVA = new TGraph (nBinsEta-1,Eta,FakeEta_MVA);
  FakeVsEta_MVA->SetMarkerColor(kBlue); 
  FakeVsEta_MVA->SetMarkerStyle(21); 
  FakeVsEta_MVA->SetMarkerSize(1); 
  FakeVsEta_MVA->SetLineColor(kBlue);
  FakeVsEta_MVA->SetFillColor(0);
  FakeVsEta_MVA->SetLineWidth(1);
  FakeVsEta_MVA->SetTitle("");

  TGraph* FakeVsEta_MVA2_WP75 = new TGraph (nBinsEta-1,Eta,FakeEta_MVA2_WP75 );
  FakeVsEta_MVA2_WP75->SetMarkerColor(kRed); 
  FakeVsEta_MVA2_WP75->SetMarkerStyle(21); 
  FakeVsEta_MVA2_WP75->SetMarkerSize(1); 
  FakeVsEta_MVA2_WP75->SetLineColor(kRed);
  FakeVsEta_MVA2_WP75->SetFillColor(0);
  FakeVsEta_MVA2_WP75->SetLineWidth(1);
  FakeVsEta_MVA2_WP75->SetTitle("");

  TGraph* FakeVsEta_MVA2_WP85 = new TGraph (nBinsEta-1,Eta,FakeEta_MVA2_WP85 );
  FakeVsEta_MVA2_WP85->SetMarkerColor(kGreen); 
  FakeVsEta_MVA2_WP85->SetMarkerStyle(21); 
  FakeVsEta_MVA2_WP85->SetMarkerSize(1); 
  FakeVsEta_MVA2_WP85->SetLineColor(kGreen);
  FakeVsEta_MVA2_WP85->SetFillColor(0);
  FakeVsEta_MVA2_WP85->SetLineWidth(1);
  FakeVsEta_MVA2_WP85->SetTitle("");

  TGraph* FakeVsEta_MVA2_WP95 = new TGraph (nBinsEta-1,Eta,FakeEta_MVA2_WP95 );
  FakeVsEta_MVA2_WP95->SetMarkerColor(kMagenta); 
  FakeVsEta_MVA2_WP95->SetMarkerStyle(21); 
  FakeVsEta_MVA2_WP95->SetMarkerSize(1); 
  FakeVsEta_MVA2_WP95->SetLineColor(kMagenta);
  FakeVsEta_MVA2_WP95->SetFillColor(0);
  FakeVsEta_MVA2_WP95->SetLineWidth(1);
  FakeVsEta_MVA2_WP95->SetTitle("");

  TMultiGraph *mg_FakeVsEta = new TMultiGraph();
  mg_FakeVsEta->Add(FakeVsEta_MVA);
  mg_FakeVsEta->Add(FakeVsEta_MVA2_WP75);
  mg_FakeVsEta->Add(FakeVsEta_MVA2_WP85);
  mg_FakeVsEta->Add(FakeVsEta_MVA2_WP95);

  TLegend* leg_FakeVsEta = new TLegend(0.6,0.75,0.8,0.88,NULL,"brNDC");
  leg_FakeVsEta->SetFillStyle(0);
  leg_FakeVsEta->SetBorderSize(0);
  leg_FakeVsEta->SetFillColor(10);
  leg_FakeVsEta->SetTextSize(0.03);

  leg_FakeVsEta->AddEntry(FakeVsEta_MVA,"old MVA");
  leg_FakeVsEta->AddEntry(FakeVsEta_MVA2_WP75,"new MVA WP75");
  leg_FakeVsEta->AddEntry(FakeVsEta_MVA2_WP85,"new MVA WP85");
  leg_FakeVsEta->AddEntry(FakeVsEta_MVA2_WP95,"new MVA WP95");

  c1->Clear();
  
  if(type=="Eff" ){
    mg_EffVsEta->Draw("APL");
    mg_EffVsEta->GetXaxis()->SetTitle("#tau |#eta|");
    mg_EffVsEta->GetXaxis()->SetRangeUser(0.,2.8);
    mg_EffVsEta->GetYaxis()->SetTitle("Efficiency");
    mg_EffVsEta->SetMinimum(0.3);
    mg_EffVsEta->SetMaximum(1.);
    leg_EffVsEta->Draw("same");
    
  }
  if(type=="FakeRate"){
    mg_FakeVsEta->Draw("APL");
    mg_FakeVsEta->GetXaxis()->SetTitle("#tau |#eta|");
    mg_FakeVsEta->GetXaxis()->SetRangeUser(0.,2.8);
    mg_FakeVsEta->GetYaxis()->SetTitle("Fake Rate");
    mg_FakeVsEta->SetMinimum(0.);
    mg_FakeVsEta->SetMaximum(1.0);
    leg_FakeVsEta->Draw("same");
  }
  
  string outputName = Form("plots/plotEfficiencies_%sVsEta",type.data());
  c1->Print(std::string(outputName).append(".png").data());
  c1->Print(std::string(outputName).append(".pdf").data());
  
  fSignal->Close();
  
}



void plotEfficienciesVsPhi(string type = "Eff"
			  )
{

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

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

  TFile *fSignal = new TFile("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_makeWP-iter1_Pythia.root","READ"); 
//   TFile *fSignal = new TFile("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_makeWP-PtEleRelaxed-iter1_Pythia.root","READ"); 
//   TFile *fSignal = new TFile("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_makeWP_Pythia.root","READ"); 
  if(fSignal->IsZombie()){
    std::cout<<"No such file!"<<std::endl;
    return;
  }

  TTree*inTree = (TTree*)fSignal->Get("tree");

  TH1F* hMVASig   = new TH1F( "hMVASig","" ,100,-1,1);
  TH1F* hMVABkg = new TH1F( "hMVABkg","" ,100,-1,1);

  float PtMin = 20;
  float PtMax = 100;
  int numPVMin = 0;
  int numPVMax = 50;
  float EtaMin = 0;
  float EtaMax = 3.0;


  float nTot = 0;
  float nPassSig_MVA = 0;
  float nPassSig_MVA2_WP75 = 0;
  float nPassSig_MVA2_WP85 = 0;
  float nPassSig_MVA2_WP95 = 0;
  float nPassBkg_MVA = 0;
  float nPassBkg_MVA2_WP75 = 0;
  float nPassBkg_MVA2_WP85 = 0;
  float nPassBkg_MVA2_WP95 = 0;

  int nBinsPhi = 7;
  float Phi [nBinsPhi];
  Phi[0] = -3.;
  Phi[1] = -2;
  Phi[2] = -1;
  Phi[3] = 0;
  Phi[4] = 1;
  Phi[5] = 2;
  Phi[6] = 3;

  float EffPhi_MVA [nBinsPhi-1];
  float EffPhi_MVA2_WP75 [nBinsPhi-1];
  float EffPhi_MVA2_WP85 [nBinsPhi-1];
  float EffPhi_MVA2_WP95 [nBinsPhi-1];
  float FakePhi_MVA [nBinsPhi-1];
  float FakePhi_MVA2_WP75 [nBinsPhi-1];
  float FakePhi_MVA2_WP85 [nBinsPhi-1];
  float FakePhi_MVA2_WP95 [nBinsPhi-1];

  for(int k=0;k<nBinsPhi-1;k++){

    TCut PUSelection(Form("NumPV>%i && NumPV<%i",numPVMin,numPVMax));
    TCut EtaSelection (Form("Tau_AbsEta>%0f && Tau_AbsEta<%0f",EtaMin,EtaMax));
    TCut PhiSelection (Form("Tau_Phi>%0f && Tau_Phi<%0f",Phi[k],Phi[k+1]));
    TCut PtSelection (Form("Tau_Pt>%0f && Tau_Pt<%0f",PtMin,PtMax));
    TCut Preselection = "Tau_AntiELoose";
    TCut Selection = PUSelection && PtSelection && EtaSelection && PhiSelection && Preselection;
        
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch"*Selection);
    inTree->Draw("Tau_AntiEMVA>>hMVABkg","Tau_GenEleMatch"*Selection);
    nTot = hMVASig->GetEntries();    
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA"*Selection);
    inTree->Draw("Tau_AntiEMVA>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA"*Selection);
    nPassSig_MVA = hMVASig->GetEntries();
    nPassBkg_MVA = hMVABkg->GetEntries();
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA2_WP75"*Selection);
    inTree->Draw("Tau_AntiEMVA>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA2_WP75"*Selection);
    nPassSig_MVA2_WP75 = hMVASig->GetEntries();
    nPassBkg_MVA2_WP75 = hMVABkg->GetEntries();
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA2_WP85"*Selection);
    inTree->Draw("Tau_AntiEMVA>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA2_WP85"*Selection);
    nPassSig_MVA2_WP85 = hMVASig->GetEntries();
    nPassBkg_MVA2_WP85 = hMVABkg->GetEntries();
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA2_WP95"*Selection);
    inTree->Draw("Tau_AntiEMVA>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA2_WP95"*Selection);
    nPassSig_MVA2_WP95 = hMVASig->GetEntries();
    nPassBkg_MVA2_WP95 = hMVABkg->GetEntries();
  
    EffPhi_MVA[k] = nPassSig_MVA/nTot;
    FakePhi_MVA[k] = nPassBkg_MVA/nTot;
    EffPhi_MVA2_WP75[k] = nPassSig_MVA2_WP75/nTot;
    FakePhi_MVA2_WP75[k] = nPassBkg_MVA2_WP75/nTot;    
    EffPhi_MVA2_WP85[k] = nPassSig_MVA2_WP85/nTot;
    FakePhi_MVA2_WP85[k] = nPassBkg_MVA2_WP85/nTot;
    EffPhi_MVA2_WP95[k] = nPassSig_MVA2_WP95/nTot;
    FakePhi_MVA2_WP95[k] = nPassBkg_MVA2_WP95/nTot;
    
    if(DEBUG){
    std::cout<<"For Phi between :"<<Phi[k]<<" and "<<Phi[k+1]<<std::endl;
//     std::cout<<" Total :"<<nTot<<std::endl;
//     std::cout<<" passing signal MVA:"<<nPassSig_MVA<<std::endl;
//     std::cout<<" passing bkg MVA:"<<nPassBkg_MVA<<std::endl;
    std::cout<<" signal eff MVA:"<<EffPhi_MVA[k]<<std::endl;
    std::cout<<" bkg eff MVA:"<<FakePhi_MVA[k]<<std::endl;
//     std::cout<<std::endl;
//     std::cout<<" passing signal MVA WP75:"<<nPassSig_MVA2_WP75<<std::endl;
//     std::cout<<" passing bkg MVA WP75:"<<nPassBkg_MVA2_WP75<<std::endl;
    std::cout<<" signal eff MVA WP75:"<<EffPhi_MVA2_WP75[k]<<std::endl;
    std::cout<<" bkg eff MVA WP75:"<<FakePhi_MVA2_WP75[k]<<std::endl;
//     std::cout<<std::endl;
//     std::cout<<" passing signal MVA WP85:"<<nPassSig_MVA2_WP85<<std::endl;
//     std::cout<<" passing bkg MVA WP85:"<<nPassBkg_MVA2_WP85<<std::endl;
    std::cout<<" signal eff MVA WP85:"<<EffPhi_MVA2_WP85[k]<<std::endl;
    std::cout<<" bkg eff MVA WP85:"<<FakePhi_MVA2_WP85[k]<<std::endl;
//     std::cout<<std::endl;
//     std::cout<<" passing signal MVA WP95:"<<nPassSig_MVA2_WP95<<std::endl;
//     std::cout<<" passing bkg MVA WP95:"<<nPassBkg_MVA2_WP95<<std::endl;
    std::cout<<" signal eff MVA WP95:"<<EffPhi_MVA2_WP95[k]<<std::endl;
    std::cout<<" bkg eff MVA WP95:"<<FakePhi_MVA2_WP95[k]<<std::endl;
//     std::cout<<std::endl; 
    }
    Phi[k] = (Phi[k]+Phi[k+1])/2;
  }


  TGraph* EffVsPhi_MVA = new TGraph (nBinsPhi-1,Phi,EffPhi_MVA);
  EffVsPhi_MVA->SetMarkerColor(kBlue); 
  EffVsPhi_MVA->SetMarkerStyle(21); 
  EffVsPhi_MVA->SetMarkerSize(1); 
  EffVsPhi_MVA->SetLineColor(kBlue);
  EffVsPhi_MVA->SetFillColor(0);
  EffVsPhi_MVA->SetLineWidth(1);
  EffVsPhi_MVA->SetTitle("");

  TGraph* EffVsPhi_MVA2_WP75 = new TGraph (nBinsPhi-1,Phi,EffPhi_MVA2_WP75);
  EffVsPhi_MVA2_WP75->SetMarkerColor(kRed); 
  EffVsPhi_MVA2_WP75->SetMarkerStyle(21); 
  EffVsPhi_MVA2_WP75->SetMarkerSize(1); 
  EffVsPhi_MVA2_WP75->SetLineColor(kRed);
  EffVsPhi_MVA2_WP75->SetFillColor(0);
  EffVsPhi_MVA2_WP75->SetLineWidth(1);
  EffVsPhi_MVA2_WP75->SetTitle("WP75");

  TGraph* EffVsPhi_MVA2_WP85 = new TGraph (nBinsPhi-1,Phi,EffPhi_MVA2_WP85);
  EffVsPhi_MVA2_WP85->SetMarkerColor(kGreen); 
  EffVsPhi_MVA2_WP85->SetMarkerStyle(21); 
  EffVsPhi_MVA2_WP85->SetMarkerSize(1); 
  EffVsPhi_MVA2_WP85->SetLineColor(kGreen);
  EffVsPhi_MVA2_WP85->SetFillColor(0);
  EffVsPhi_MVA2_WP85->SetLineWidth(1);
  EffVsPhi_MVA2_WP85->SetTitle("WP85");

  TGraph* EffVsPhi_MVA2_WP95 = new TGraph (nBinsPhi-1,Phi,EffPhi_MVA2_WP95);
  EffVsPhi_MVA2_WP95->SetMarkerColor(kMagenta); 
  EffVsPhi_MVA2_WP95->SetMarkerStyle(21); 
  EffVsPhi_MVA2_WP95->SetMarkerSize(1); 
  EffVsPhi_MVA2_WP95->SetLineColor(kMagenta);
  EffVsPhi_MVA2_WP95->SetFillColor(0);
  EffVsPhi_MVA2_WP95->SetLineWidth(1);
  EffVsPhi_MVA2_WP95->SetTitle("WP95");

  TMultiGraph *mg_EffVsPhi = new TMultiGraph();
  mg_EffVsPhi->Add(EffVsPhi_MVA);
  mg_EffVsPhi->Add(EffVsPhi_MVA2_WP75);
  mg_EffVsPhi->Add(EffVsPhi_MVA2_WP85);
  mg_EffVsPhi->Add(EffVsPhi_MVA2_WP95);

  TLegend* leg_EffVsPhi = new TLegend(0.6,0.75,0.8,0.88,NULL,"brNDC");
  leg_EffVsPhi->SetFillStyle(0);
  leg_EffVsPhi->SetBorderSize(0);
  leg_EffVsPhi->SetFillColor(10);
  leg_EffVsPhi->SetTextSize(0.03);

  leg_EffVsPhi->AddEntry(EffVsPhi_MVA,"old MVA");
  leg_EffVsPhi->AddEntry(EffVsPhi_MVA2_WP75,"new MVA WP75");
  leg_EffVsPhi->AddEntry(EffVsPhi_MVA2_WP85,"new MVA WP85");
  leg_EffVsPhi->AddEntry(EffVsPhi_MVA2_WP95,"new MVA WP95");

  ///FakeRates
  TGraph* FakeVsPhi_MVA = new TGraph (nBinsPhi-1,Phi,FakePhi_MVA);
  FakeVsPhi_MVA->SetMarkerColor(kBlue); 
  FakeVsPhi_MVA->SetMarkerStyle(21); 
  FakeVsPhi_MVA->SetMarkerSize(1); 
  FakeVsPhi_MVA->SetLineColor(kBlue);
  FakeVsPhi_MVA->SetFillColor(0);
  FakeVsPhi_MVA->SetLineWidth(1);
  FakeVsPhi_MVA->SetTitle("");

  TGraph* FakeVsPhi_MVA2_WP75 = new TGraph (nBinsPhi-1,Phi,FakePhi_MVA2_WP75 );
  FakeVsPhi_MVA2_WP75->SetMarkerColor(kRed); 
  FakeVsPhi_MVA2_WP75->SetMarkerStyle(21); 
  FakeVsPhi_MVA2_WP75->SetMarkerSize(1); 
  FakeVsPhi_MVA2_WP75->SetLineColor(kRed);
  FakeVsPhi_MVA2_WP75->SetFillColor(0);
  FakeVsPhi_MVA2_WP75->SetLineWidth(1);
  FakeVsPhi_MVA2_WP75->SetTitle("");

  TGraph* FakeVsPhi_MVA2_WP85 = new TGraph (nBinsPhi-1,Phi,FakePhi_MVA2_WP85 );
  FakeVsPhi_MVA2_WP85->SetMarkerColor(kGreen); 
  FakeVsPhi_MVA2_WP85->SetMarkerStyle(21); 
  FakeVsPhi_MVA2_WP85->SetMarkerSize(1); 
  FakeVsPhi_MVA2_WP85->SetLineColor(kGreen);
  FakeVsPhi_MVA2_WP85->SetFillColor(0);
  FakeVsPhi_MVA2_WP85->SetLineWidth(1);
  FakeVsPhi_MVA2_WP85->SetTitle("");

  TGraph* FakeVsPhi_MVA2_WP95 = new TGraph (nBinsPhi-1,Phi,FakePhi_MVA2_WP95 );
  FakeVsPhi_MVA2_WP95->SetMarkerColor(kMagenta); 
  FakeVsPhi_MVA2_WP95->SetMarkerStyle(21); 
  FakeVsPhi_MVA2_WP95->SetMarkerSize(1); 
  FakeVsPhi_MVA2_WP95->SetLineColor(kMagenta);
  FakeVsPhi_MVA2_WP95->SetFillColor(0);
  FakeVsPhi_MVA2_WP95->SetLineWidth(1);
  FakeVsPhi_MVA2_WP95->SetTitle("");

  TMultiGraph *mg_FakeVsPhi = new TMultiGraph();
  mg_FakeVsPhi->Add(FakeVsPhi_MVA);
  mg_FakeVsPhi->Add(FakeVsPhi_MVA2_WP75);
  mg_FakeVsPhi->Add(FakeVsPhi_MVA2_WP85);
  mg_FakeVsPhi->Add(FakeVsPhi_MVA2_WP95);

  TLegend* leg_FakeVsPhi = new TLegend(0.6,0.75,0.8,0.88,NULL,"brNDC");
  leg_FakeVsPhi->SetFillStyle(0);
  leg_FakeVsPhi->SetBorderSize(0);
  leg_FakeVsPhi->SetFillColor(10);
  leg_FakeVsPhi->SetTextSize(0.03);

  leg_FakeVsPhi->AddEntry(FakeVsPhi_MVA,"old MVA");
  leg_FakeVsPhi->AddEntry(FakeVsPhi_MVA2_WP75,"new MVA WP75");
  leg_FakeVsPhi->AddEntry(FakeVsPhi_MVA2_WP85,"new MVA WP85");
  leg_FakeVsPhi->AddEntry(FakeVsPhi_MVA2_WP95,"new MVA WP95");

  c1->Clear();

  if(type=="Eff" ){
    mg_EffVsPhi->Draw("APL");
    mg_EffVsPhi->GetXaxis()->SetTitle("#tau #phi");
    mg_EffVsPhi->GetXaxis()->SetRangeUser(-3,3);
    mg_EffVsPhi->GetYaxis()->SetTitle("Efficiency");
    mg_EffVsPhi->SetMinimum(0.3);
    mg_EffVsPhi->SetMaximum(1.);
    leg_EffVsPhi->Draw("same");
  }
  if(type=="FakeRate"){
    mg_FakeVsPhi->Draw("APL");
    mg_FakeVsPhi->GetXaxis()->SetTitle("#tau #phi");
    mg_FakeVsPhi->GetXaxis()->SetRangeUser(-3,3);
    mg_FakeVsPhi->GetYaxis()->SetTitle("Fake Rate");
    mg_FakeVsPhi->SetMinimum(0.);
    mg_FakeVsPhi->SetMaximum(1.0);
    leg_FakeVsPhi->Draw("same");
  }

  string outputName = Form("plots/plotEfficiencies_%sVsPhi",type.data());
  c1->Print(std::string(outputName).append(".png").data());
  c1->Print(std::string(outputName).append(".pdf").data());

  fSignal->Close();

}





void plotEfficienciesVsPV(string type = "Eff"
			  )
{

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

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

  TFile *fSignal = new TFile("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_makeWP-iter1_Pythia.root","READ"); 
//   TFile *fSignal = new TFile("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_makeWP-PtEleRelaxed-iter1_Pythia.root","READ"); 
//   TFile *fSignal = new TFile("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_makeWP_Pythia.root","READ"); 
  if(fSignal->IsZombie()){
    std::cout<<"No such file!"<<std::endl;
    return;
  }

  TTree*inTree = (TTree*)fSignal->Get("tree");

  TH1F* hMVASig   = new TH1F( "hMVASig","" ,100,-1,1);
  TH1F* hMVABkg = new TH1F( "hMVABkg","" ,100,-1,1);

  float PtMin = 20;
  float PtMax = 100;
  float PhiMin = -3.0;
  float PhiMax = 3.0;
  float EtaMin = 0;
  float EtaMax = 3.0;


  float nTot = 0;
  float nPassSig_MVA = 0;
  float nPassSig_MVA2_WP75 = 0;
  float nPassSig_MVA2_WP85 = 0;
  float nPassSig_MVA2_WP95 = 0;
  float nPassBkg_MVA = 0;
  float nPassBkg_MVA2_WP75 = 0;
  float nPassBkg_MVA2_WP85 = 0;
  float nPassBkg_MVA2_WP95 = 0;

  int nBinsPV = 6;
  float PV [nBinsPV];
  PV[0] = 0.;
  PV[1] = 5;
  PV[2] = 10;
  PV[3] = 15;
  PV[4] = 20;
  PV[5] = 30;

  float EffPV_MVA [nBinsPV-1];
  float EffPV_MVA2_WP75 [nBinsPV-1];
  float EffPV_MVA2_WP85 [nBinsPV-1];
  float EffPV_MVA2_WP95 [nBinsPV-1];
  float FakePV_MVA [nBinsPV-1];
  float FakePV_MVA2_WP75 [nBinsPV-1];
  float FakePV_MVA2_WP85 [nBinsPV-1];
  float FakePV_MVA2_WP95 [nBinsPV-1];

  for(int k=0;k<nBinsPV-1;k++){

    TCut PUSelection(Form("NumPV>%0f && NumPV<%0f",PV[k],PV[k+1]));
    TCut EtaSelection (Form("Tau_AbsEta>%0f && Tau_AbsEta<%0f",EtaMin,EtaMax));
    TCut PhiSelection (Form("Tau_Phi>%0f && Tau_Phi<%0f",PhiMin,PhiMax));
    TCut PtSelection (Form("Tau_Pt>%0f && Tau_Pt<%0f",PtMin,PtMax));
    TCut Preselection = "Tau_AntiELoose";
    TCut Selection = PUSelection && PtSelection && EtaSelection && PhiSelection && Preselection;
        
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch"*Selection);
    inTree->Draw("Tau_AntiEMVA>>hMVABkg","Tau_GenEleMatch"*Selection);
    nTot = hMVASig->GetEntries();    
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA"*Selection);
    inTree->Draw("Tau_AntiEMVA>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA"*Selection);
    nPassSig_MVA = hMVASig->GetEntries();
    nPassBkg_MVA = hMVABkg->GetEntries();
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA2_WP75"*Selection);
    inTree->Draw("Tau_AntiEMVA>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA2_WP75"*Selection);
    nPassSig_MVA2_WP75 = hMVASig->GetEntries();
    nPassBkg_MVA2_WP75 = hMVABkg->GetEntries();
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA2_WP85"*Selection);
    inTree->Draw("Tau_AntiEMVA>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA2_WP85"*Selection);
    nPassSig_MVA2_WP85 = hMVASig->GetEntries();
    nPassBkg_MVA2_WP85 = hMVABkg->GetEntries();
    inTree->Draw("Tau_AntiEMVA>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA2_WP95"*Selection);
    inTree->Draw("Tau_AntiEMVA>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA2_WP95"*Selection);
    nPassSig_MVA2_WP95 = hMVASig->GetEntries();
    nPassBkg_MVA2_WP95 = hMVABkg->GetEntries();
  
    EffPV_MVA[k] = nPassSig_MVA/nTot;
    FakePV_MVA[k] = nPassBkg_MVA/nTot;
    EffPV_MVA2_WP75[k] = nPassSig_MVA2_WP75/nTot;
    FakePV_MVA2_WP75[k] = nPassBkg_MVA2_WP75/nTot;    
    EffPV_MVA2_WP85[k] = nPassSig_MVA2_WP85/nTot;
    FakePV_MVA2_WP85[k] = nPassBkg_MVA2_WP85/nTot;
    EffPV_MVA2_WP95[k] = nPassSig_MVA2_WP95/nTot;
    FakePV_MVA2_WP95[k] = nPassBkg_MVA2_WP95/nTot;
    
    if(DEBUG){
    std::cout<<"For PV between :"<<PV[k]<<" and "<<PV[k+1]<<std::endl;
//     std::cout<<" Total :"<<nTot<<std::endl;
//     std::cout<<" passing signal MVA:"<<nPassSig_MVA<<std::endl;
//     std::cout<<" passing bkg MVA:"<<nPassBkg_MVA<<std::endl;
    std::cout<<" signal eff MVA:"<<EffPV_MVA[k]<<std::endl;
    std::cout<<" bkg eff MVA:"<<FakePV_MVA[k]<<std::endl;
//     std::cout<<std::endl;
//     std::cout<<" passing signal MVA WP75:"<<nPassSig_MVA2_WP75<<std::endl;
//     std::cout<<" passing bkg MVA WP75:"<<nPassBkg_MVA2_WP75<<std::endl;
    std::cout<<" signal eff MVA WP75:"<<EffPV_MVA2_WP75[k]<<std::endl;
    std::cout<<" bkg eff MVA WP75:"<<FakePV_MVA2_WP75[k]<<std::endl;
//     std::cout<<std::endl;
//     std::cout<<" passing signal MVA WP85:"<<nPassSig_MVA2_WP85<<std::endl;
//     std::cout<<" passing bkg MVA WP85:"<<nPassBkg_MVA2_WP85<<std::endl;
    std::cout<<" signal eff MVA WP85:"<<EffPV_MVA2_WP85[k]<<std::endl;
    std::cout<<" bkg eff MVA WP85:"<<FakePV_MVA2_WP85[k]<<std::endl;
//     std::cout<<std::endl;
//     std::cout<<" passing signal MVA WP95:"<<nPassSig_MVA2_WP95<<std::endl;
//     std::cout<<" passing bkg MVA WP95:"<<nPassBkg_MVA2_WP95<<std::endl;
    std::cout<<" signal eff MVA WP95:"<<EffPV_MVA2_WP95[k]<<std::endl;
    std::cout<<" bkg eff MVA WP95:"<<FakePV_MVA2_WP95[k]<<std::endl;
//     std::cout<<std::endl; 
    }
    PV[k] = (PV[k]+PV[k+1])/2;
  }


  TGraph* EffVsPV_MVA = new TGraph (nBinsPV-1,PV,EffPV_MVA);
  EffVsPV_MVA->SetMarkerColor(kBlue); 
  EffVsPV_MVA->SetMarkerStyle(21); 
  EffVsPV_MVA->SetMarkerSize(1); 
  EffVsPV_MVA->SetLineColor(kBlue);
  EffVsPV_MVA->SetFillColor(0);
  EffVsPV_MVA->SetLineWidth(1);
  EffVsPV_MVA->SetTitle("");

  TGraph* EffVsPV_MVA2_WP75 = new TGraph (nBinsPV-1,PV,EffPV_MVA2_WP75);
  EffVsPV_MVA2_WP75->SetMarkerColor(kRed); 
  EffVsPV_MVA2_WP75->SetMarkerStyle(21); 
  EffVsPV_MVA2_WP75->SetMarkerSize(1); 
  EffVsPV_MVA2_WP75->SetLineColor(kRed);
  EffVsPV_MVA2_WP75->SetFillColor(0);
  EffVsPV_MVA2_WP75->SetLineWidth(1);
  EffVsPV_MVA2_WP75->SetTitle("WP75");

  TGraph* EffVsPV_MVA2_WP85 = new TGraph (nBinsPV-1,PV,EffPV_MVA2_WP85);
  EffVsPV_MVA2_WP85->SetMarkerColor(kGreen); 
  EffVsPV_MVA2_WP85->SetMarkerStyle(21); 
  EffVsPV_MVA2_WP85->SetMarkerSize(1); 
  EffVsPV_MVA2_WP85->SetLineColor(kGreen);
  EffVsPV_MVA2_WP85->SetFillColor(0);
  EffVsPV_MVA2_WP85->SetLineWidth(1);
  EffVsPV_MVA2_WP85->SetTitle("WP85");

  TGraph* EffVsPV_MVA2_WP95 = new TGraph (nBinsPV-1,PV,EffPV_MVA2_WP95);
  EffVsPV_MVA2_WP95->SetMarkerColor(kMagenta); 
  EffVsPV_MVA2_WP95->SetMarkerStyle(21); 
  EffVsPV_MVA2_WP95->SetMarkerSize(1); 
  EffVsPV_MVA2_WP95->SetLineColor(kMagenta);
  EffVsPV_MVA2_WP95->SetFillColor(0);
  EffVsPV_MVA2_WP95->SetLineWidth(1);
  EffVsPV_MVA2_WP95->SetTitle("WP95");

  TMultiGraph *mg_EffVsPV = new TMultiGraph();
  mg_EffVsPV->Add(EffVsPV_MVA);
  mg_EffVsPV->Add(EffVsPV_MVA2_WP75);
  mg_EffVsPV->Add(EffVsPV_MVA2_WP85);
  mg_EffVsPV->Add(EffVsPV_MVA2_WP95);

  TLegend* leg_EffVsPV = new TLegend(0.6,0.75,0.8,0.88,NULL,"brNDC");
  leg_EffVsPV->SetFillStyle(0);
  leg_EffVsPV->SetBorderSize(0);
  leg_EffVsPV->SetFillColor(10);
  leg_EffVsPV->SetTextSize(0.03);

  leg_EffVsPV->AddEntry(EffVsPV_MVA,"old MVA");
  leg_EffVsPV->AddEntry(EffVsPV_MVA2_WP75,"new MVA WP75");
  leg_EffVsPV->AddEntry(EffVsPV_MVA2_WP85,"new MVA WP85");
  leg_EffVsPV->AddEntry(EffVsPV_MVA2_WP95,"new MVA WP95");

  ///FakeRates
  TGraph* FakeVsPV_MVA = new TGraph (nBinsPV-1,PV,FakePV_MVA);
  FakeVsPV_MVA->SetMarkerColor(kBlue); 
  FakeVsPV_MVA->SetMarkerStyle(21); 
  FakeVsPV_MVA->SetMarkerSize(1); 
  FakeVsPV_MVA->SetLineColor(kBlue);
  FakeVsPV_MVA->SetFillColor(0);
  FakeVsPV_MVA->SetLineWidth(1);
  FakeVsPV_MVA->SetTitle("");

  TGraph* FakeVsPV_MVA2_WP75 = new TGraph (nBinsPV-1,PV,FakePV_MVA2_WP75 );
  FakeVsPV_MVA2_WP75->SetMarkerColor(kRed); 
  FakeVsPV_MVA2_WP75->SetMarkerStyle(21); 
  FakeVsPV_MVA2_WP75->SetMarkerSize(1); 
  FakeVsPV_MVA2_WP75->SetLineColor(kRed);
  FakeVsPV_MVA2_WP75->SetFillColor(0);
  FakeVsPV_MVA2_WP75->SetLineWidth(1);
  FakeVsPV_MVA2_WP75->SetTitle("");

  TGraph* FakeVsPV_MVA2_WP85 = new TGraph (nBinsPV-1,PV,FakePV_MVA2_WP85 );
  FakeVsPV_MVA2_WP85->SetMarkerColor(kGreen); 
  FakeVsPV_MVA2_WP85->SetMarkerStyle(21); 
  FakeVsPV_MVA2_WP85->SetMarkerSize(1); 
  FakeVsPV_MVA2_WP85->SetLineColor(kGreen);
  FakeVsPV_MVA2_WP85->SetFillColor(0);
  FakeVsPV_MVA2_WP85->SetLineWidth(1);
  FakeVsPV_MVA2_WP85->SetTitle("");

  TGraph* FakeVsPV_MVA2_WP95 = new TGraph (nBinsPV-1,PV,FakePV_MVA2_WP95 );
  FakeVsPV_MVA2_WP95->SetMarkerColor(kMagenta); 
  FakeVsPV_MVA2_WP95->SetMarkerStyle(21); 
  FakeVsPV_MVA2_WP95->SetMarkerSize(1); 
  FakeVsPV_MVA2_WP95->SetLineColor(kMagenta);
  FakeVsPV_MVA2_WP95->SetFillColor(0);
  FakeVsPV_MVA2_WP95->SetLineWidth(1);
  FakeVsPV_MVA2_WP95->SetTitle("");

  TMultiGraph *mg_FakeVsPV = new TMultiGraph();
  mg_FakeVsPV->Add(FakeVsPV_MVA);
  mg_FakeVsPV->Add(FakeVsPV_MVA2_WP75);
  mg_FakeVsPV->Add(FakeVsPV_MVA2_WP85);
  mg_FakeVsPV->Add(FakeVsPV_MVA2_WP95);

  TLegend* leg_FakeVsPV = new TLegend(0.6,0.75,0.8,0.88,NULL,"brNDC");
  leg_FakeVsPV->SetFillStyle(0);
  leg_FakeVsPV->SetBorderSize(0);
  leg_FakeVsPV->SetFillColor(10);
  leg_FakeVsPV->SetTextSize(0.03);

  leg_FakeVsPV->AddEntry(FakeVsPV_MVA,"old MVA");
  leg_FakeVsPV->AddEntry(FakeVsPV_MVA2_WP75,"new MVA WP75");
  leg_FakeVsPV->AddEntry(FakeVsPV_MVA2_WP85,"new MVA WP85");
  leg_FakeVsPV->AddEntry(FakeVsPV_MVA2_WP95,"new MVA WP95");

  c1->Clear();

  if(type=="Eff" ){
    mg_EffVsPV->Draw("APL");
    mg_EffVsPV->GetXaxis()->SetTitle("# PV");
    mg_EffVsPV->GetXaxis()->SetRangeUser(0,50);
    mg_EffVsPV->GetYaxis()->SetTitle("Efficiency");
    mg_EffVsPV->SetMinimum(0.3);
    mg_EffVsPV->SetMaximum(1.);
    leg_EffVsPV->Draw("same");
  }
  if(type=="FakeRate"){
    mg_FakeVsPV->Draw("APL");
    mg_FakeVsPV->GetXaxis()->SetTitle("# PV");
    mg_FakeVsPV->GetXaxis()->SetRangeUser(0,50);
    mg_FakeVsPV->GetYaxis()->SetTitle("Fake Rate");
    mg_FakeVsPV->SetMinimum(0.);
    mg_FakeVsPV->SetMaximum(1.0);
    leg_FakeVsPV->Draw("same");
  }

  string outputName = Form("plots/plotEfficiencies_%sVsPV",type.data());
  c1->Print(std::string(outputName).append(".png").data());
  c1->Print(std::string(outputName).append(".pdf").data());

  fSignal->Close();

}


void plotAll(){
  plotEfficienciesVsPt("Eff");
  plotEfficienciesVsPV("Eff");
  plotEfficienciesVsPhi("Eff");
  plotEfficienciesVsEta("Eff");

  plotEfficienciesVsPt("FakeRate");
  plotEfficienciesVsPV("FakeRate");
  plotEfficienciesVsPhi("FakeRate");
  plotEfficienciesVsEta("FakeRate");
}

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
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

#define DEBUG true

void plotEfficiency(string type = "FakeRate",
		    string variable = "Pt",
		    string category = "wGwGSFwPFMVA",
		    const TString& xAxisTitle = "P_{T} #tau",
		    float xMin = 15, 
		    float xMax = 150,
		    int nBins = 1	
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

  std::string inputFileName = "/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_makeWP_DYJetsToLL.root";
  TFile* inputFile = new TFile (inputFileName.data(),"READ");
  if(inputFile->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
  TTree* inTree = (TTree*)inputFile->Get("tree");
  TH1F* hMVASig   = new TH1F( "hMVASig","" ,100,-1,1);
  TH1F* hMVABkg = new TH1F( "hMVABkg","" ,100,-1,1);


  int numPVMin = 0;
  int numPVMax = 50;
  float AbsEtaMin = 0;
  float AbsEtaMax = 3.0;
  float PhiMin = -3;
  float PhiMax = 3;
  float PtMin = 15;
  float PtMax = 150;

  float nTotSig = 0;
  float nTotBkg = 0;
  float nPassSig_Loose = 0;
  float nPassSig_Medium = 0;
  float nPassSig_Tight = 0;
  float nPassSig_MVA = 0;
  float nPassSig_MVA2_WP75 = 0;
  float nPassSig_MVA2_WP85 = 0;
  float nPassSig_MVA2_WP95 = 0;
  float nPassBkg_Loose = 0;
  float nPassBkg_Medium = 0;
  float nPassBkg_Tight = 0;
  float nPassBkg_MVA = 0;
  float nPassBkg_MVA2_WP75 = 0;
  float nPassBkg_MVA2_WP85 = 0;
  float nPassBkg_MVA2_WP95 = 0;

  float Pt [nBins+1];
  float AbsEta [nBins+1];
  float Phi [nBins+1];
  float numPV [nBins+1];

  float YMin = 1.0;

  for(int i=0;i<nBins+1;i++){
    if(variable == "Pt") Pt[i]          =  xMin + i*(xMax-xMin)/nBins;
    if(variable == "AbsEta") AbsEta[i]  =  xMin + i*(xMax-xMin)/nBins;
    if(variable == "Phi") Phi[i]        =  xMin + i*(xMax-xMin)/nBins;
    if(variable == "numPV") numPV[i]    =  xMin + i*(xMax-xMin)/nBins;
  }

  float ErrSig [nBins];
  float ErrBkg [nBins];
  float Eff_Loose [nBins];
  float Eff_Medium [nBins];
  float Eff_Tight [nBins];
  float Eff_MVA [nBins];
  float Eff_MVA2_WP75 [nBins];
  float Eff_MVA2_WP85 [nBins];
  float Eff_MVA2_WP95 [nBins];
  float Fake_Loose [nBins];
  float Fake_Medium [nBins];
  float Fake_Tight [nBins];
  float Fake_MVA [nBins];
  float Fake_MVA2_WP75 [nBins];
  float Fake_MVA2_WP85 [nBins];
  float Fake_MVA2_WP95 [nBins];

  for(int k=0;k<nBins;k++){
    Eff_Loose [0] = 0;
    Eff_Medium [0] = 0;
    Eff_Tight [0] = 0;
    Eff_MVA [k] = 0;
    Eff_MVA2_WP75 [k] = 0;
    Eff_MVA2_WP85 [k] = 0;
    Eff_MVA2_WP95 [k] = 0;
    Fake_Loose [0] = 0;
    Fake_Medium [0] = 0;
    Fake_Tight [0] = 0;
    Fake_MVA [k] = 0;
    Fake_MVA2_WP75 [k] = 0;
    Fake_MVA2_WP85 [k] = 0;
    Fake_MVA2_WP95 [k] = 0;

    TCut Preselection = "";
    TCut PUSelection(Form("NumPV>%i && NumPV<%i",numPVMin,numPVMax));
    TCut EtaSelection (Form("Tau_AbsEta>%0f && Tau_AbsEta<%0f",AbsEtaMin,AbsEtaMax));
    TCut PhiSelection (Form("Tau_Phi>%0f && Tau_Phi<%0f",PhiMin,PhiMax));
    TCut PtSelection (Form("Tau_Pt>%0f && Tau_Pt<%0f",PtMin,PtMax));
    if(variable == "Pt") PtSelection      = Form("Tau_Pt>%0f && Tau_Pt<%0f",Pt[k],Pt[k+1]);
    if(variable == "AbsEta") EtaSelection = Form("Tau_AbsEta>%0f && Tau_AbsEta<%0f",AbsEta[k],AbsEta[k+1]);
    if(variable == "Phi") PhiSelection    = Form("Tau_Phi>%0f && Tau_Phi<%0f",Phi[k],Phi[k+1]);
    if(variable == "numPV") PUSelection   = Form("NumPV>%0f && NumPV<%0f",numPV[k],numPV[k+1]);

    TCut Selection = PUSelection && PtSelection && EtaSelection && PhiSelection && Preselection ;


    inTree->Draw("Tau_Phi>>hMVASig","Tau_GenHadMatch"*Selection);
    inTree->Draw("Tau_Phi>>hMVABkg","Tau_GenEleMatch"*Selection);
    nTotSig = hMVASig->GetEntries(); 
    nTotBkg = hMVABkg->GetEntries(); 
//     cout<<"Tot sig: "<<nTotSig<<endl;
//     cout<<"Tot bkg : "<<nTotBkg<<endl;
    hMVASig->Reset();

    TCut CategorySelection = "";
    if (category == "NoEleMatch") CategorySelection   = "NoEleMatch"; 
    if (category == "woG") CategorySelection          = "woG"; 
    if (category == "wGwoGSF") CategorySelection      = "wGwoGSF";
    if (category == "wGwGSFwoPFMVA")CategorySelection = "wGwGSFwoPFMVA";
    if (category == "wGwGSFwPFMVA")CategorySelection  = "wGwGSFwPFMVA";
    Selection = PUSelection && PtSelection && EtaSelection && PhiSelection && Preselection && CategorySelection;
    inTree->Draw("Tau_Phi>>hMVASig","Tau_GenHadMatch && Tau_AntiELoose"*Selection);
    inTree->Draw("Tau_Phi>>hMVABkg","Tau_GenEleMatch && Tau_AntiELoose"*Selection);
    nPassSig_Loose = hMVASig->GetEntries();
    nPassBkg_Loose = hMVABkg->GetEntries();
    hMVASig->Reset();
    hMVABkg->Reset();
    inTree->Draw("Tau_Phi>>hMVASig","Tau_GenHadMatch && Tau_AntiEMedium"*Selection);
    inTree->Draw("Tau_Phi>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMedium"*Selection);
    nPassSig_Medium = hMVASig->GetEntries();
    nPassBkg_Medium = hMVABkg->GetEntries();
    hMVASig->Reset();
    hMVABkg->Reset();
    inTree->Draw("Tau_Phi>>hMVASig","Tau_GenHadMatch && Tau_AntiETight"*Selection);
    inTree->Draw("Tau_Phi>>hMVABkg","Tau_GenEleMatch && Tau_AntiETight"*Selection);
    nPassSig_Tight = hMVASig->GetEntries();
    nPassBkg_Tight = hMVABkg->GetEntries();
    hMVASig->Reset();
    hMVABkg->Reset();
    inTree->Draw("Tau_Phi>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA"*Selection);
    inTree->Draw("Tau_Phi>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA"*Selection);
    nPassSig_MVA = hMVASig->GetEntries();
    nPassBkg_MVA = hMVABkg->GetEntries();
    hMVASig->Reset();
    hMVABkg->Reset();
//     cout<<"Pass sig oldMVA: "<<nPassSig_MVA<<endl;
//     cout<<"Pass bkg oldMVA: "<<nPassBkg_MVA<<endl;
    inTree->Draw("Tau_Phi>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA2_WP75"*Selection);
    inTree->Draw("Tau_Phi>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA2_WP75"*Selection);
    nPassSig_MVA2_WP75 = hMVASig->GetEntries();
    nPassBkg_MVA2_WP75 = hMVABkg->GetEntries();
    hMVASig->Reset();
    hMVABkg->Reset();
//     cout<<"Pass sig 75: "<<nPassSig_MVA2_WP75<<endl;
//     cout<<"Pass bkg 75: "<<nPassBkg_MVA2_WP75<<endl; 
    inTree->Draw("Tau_Phi>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA2_WP85"*Selection);
    inTree->Draw("Tau_Phi>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA2_WP85"*Selection);
    nPassSig_MVA2_WP85 = hMVASig->GetEntries();
    nPassBkg_MVA2_WP85 = hMVABkg->GetEntries();
    hMVASig->Reset();
    hMVABkg->Reset();
//     cout<<"Pass sig 85: "<<nPassSig_MVA2_WP85<<endl;
//     cout<<"Pass bkg 85: "<<nPassBkg_MVA2_WP85<<endl;
    inTree->Draw("Tau_Phi>>hMVASig","Tau_GenHadMatch && Tau_AntiEMVA2_WP95"*Selection);
    inTree->Draw("Tau_Phi>>hMVABkg","Tau_GenEleMatch && Tau_AntiEMVA2_WP95"*Selection);
    nPassSig_MVA2_WP95 = hMVASig->GetEntries();
    nPassBkg_MVA2_WP95 = hMVABkg->GetEntries();
//     cout<<"Pass sig 95: "<<nPassSig_MVA2_WP95<<endl;
//     cout<<"Pass bkg 95: "<<nPassBkg_MVA2_WP95<<endl;
    hMVASig->Reset();
    hMVABkg->Reset();

    ErrSig[k] = 1/sqrt(nTotSig);
    ErrBkg[k] = 1/sqrt(nTotBkg);
    Eff_Loose[k] = nPassSig_Loose/nTotSig;
    Fake_Loose[k] = nPassBkg_Loose/nTotBkg;
    Eff_Medium[k] = nPassSig_Medium/nTotSig;
    Fake_Medium[k] = nPassBkg_Medium/nTotBkg;
    Eff_Tight[k] = nPassSig_Tight/nTotSig;
    Fake_Tight[k] = nPassBkg_Tight/nTotBkg;
    Eff_MVA[k] = nPassSig_MVA/nTotSig;
    Fake_MVA[k] = nPassBkg_MVA/nTotBkg;
    Eff_MVA2_WP75[k] = nPassSig_MVA2_WP75/nTotSig;
    Fake_MVA2_WP75[k] = nPassBkg_MVA2_WP75/nTotBkg;    
    Eff_MVA2_WP85[k] = nPassSig_MVA2_WP85/nTotSig;
    Fake_MVA2_WP85[k] = nPassBkg_MVA2_WP85/nTotBkg;
    Eff_MVA2_WP95[k] = nPassSig_MVA2_WP95/nTotSig;
    Fake_MVA2_WP95[k] = nPassBkg_MVA2_WP95/nTotBkg;

    if(Fake_MVA[k]<YMin) YMin = Fake_MVA[k] ;
    if(Fake_MVA2_WP75[k]<YMin) YMin = Fake_MVA2_WP75[k] ;
    if(Fake_MVA2_WP85[k]<YMin) YMin = Fake_MVA2_WP85[k] ;
    if(Fake_MVA2_WP95[k]<YMin) YMin = Fake_MVA2_WP95[k] ;

    if(DEBUG){
      std::cout<<"For category : "<<category.data()<<endl;
      if(variable == "Pt") std::cout<<" Pt between :"<<Pt[k]<<" and "<<Pt[k+1]<<std::endl;
      else std::cout<<" Pt between :"<<PtMin<<" and "<<PtMax<<std::endl;
      if(variable == "AbsEta") std::cout<<" AbsEta between :"<<AbsEta[k]<<" and "<<AbsEta[k]<<std::endl;
      else std::cout<<" AbsEta between :"<<AbsEtaMin<<" and "<<AbsEtaMax<<std::endl;
      if(variable == "Phi") std::cout<<" Phi between :"<<Phi[k]<<" and "<<Phi[k]<<std::endl;
      else std::cout<<" Phi between :"<<PhiMin<<" and "<<PhiMax<<std::endl;
      if(variable == "numPV")std::cout<<" numPV between :"<<numPV[k]<<" and "<<numPV[k]<<std::endl;
      else std::cout<<" numPV between :"<<numPVMin<<" and "<<numPVMax<<std::endl;
      std::cout<<"  signal eff Loose:"<<Eff_Loose[k]<<std::endl;
      std::cout<<"  bkg eff Loose:"<<Fake_Loose[k]<<std::endl;
      std::cout<<"  signal eff Medium:"<<Eff_Medium[k]<<std::endl;
      std::cout<<"  bkg eff Medium:"<<Fake_Medium[k]<<std::endl;
      std::cout<<"  signal eff Tight:"<<Eff_Tight[k]<<std::endl;
      std::cout<<"  bkg eff Tight:"<<Fake_Tight[k]<<std::endl;
      std::cout<<"  signal eff MVA:"<<Eff_MVA[k]<<std::endl;
      std::cout<<"  bkg eff MVA:"<<Fake_MVA[k]<<std::endl;
      std::cout<<"  signal eff MVA2 WP75:"<<Eff_MVA2_WP75[k]<<std::endl;
      std::cout<<"  bkg eff MVA2 WP75:"<<Fake_MVA2_WP75[k]<<std::endl;
      std::cout<<"  signal eff MVA2 WP85:"<<Eff_MVA2_WP85[k]<<std::endl;
      std::cout<<"  bkg eff MVA2 WP85:"<<Fake_MVA2_WP85[k]<<std::endl;
      std::cout<<"  signal eff MVA2 WP95:"<<Eff_MVA2_WP95[k]<<std::endl;
      std::cout<<"  bkg eff MVA2 WP95:"<<Fake_MVA2_WP95[k]<<std::endl;
    }

    //For plots
    if(variable == "Pt")  Pt[k] = (Pt[k]+Pt[k+1])/2;
    if(variable == "AbsEta")  AbsEta[k] = (AbsEta[k]+AbsEta[k+1])/2;
    if(variable == "Phi")  Phi[k] = (Phi[k]+Phi[k+1])/2;
    if(variable == "numPV")  numPV[k] = (numPV[k]+numPV[k+1])/2;

  }

  TGraphErrors* GraphEff_MVA = new TGraphErrors (nBins,Pt,Eff_MVA,0,ErrSig);
  if(variable == "AbsEta")GraphEff_MVA = new TGraphErrors (nBins,AbsEta,Eff_MVA,0,ErrSig);
  if(variable == "Phi")GraphEff_MVA = new TGraphErrors (nBins,Phi,Eff_MVA,0,ErrSig);
  if(variable == "numPV")GraphEff_MVA = new TGraphErrors (nBins,numPV,Eff_MVA,0,ErrSig);
  GraphEff_MVA->SetMarkerColor(kBlue); 
  GraphEff_MVA->SetMarkerStyle(21); 
  GraphEff_MVA->SetMarkerSize(1); 
  GraphEff_MVA->SetLineColor(kBlue);
  GraphEff_MVA->SetFillColor(0);
  GraphEff_MVA->SetLineWidth(1);
  GraphEff_MVA->SetTitle("");

  TGraphErrors* GraphEff_MVA2_WP75 = new TGraphErrors (nBins,Pt,Eff_MVA2_WP75,0,ErrSig);
  if(variable == "AbsEta")GraphEff_MVA2_WP75 = new TGraphErrors (nBins,AbsEta,Eff_MVA,0,ErrSig);
  if(variable == "Phi")GraphEff_MVA2_WP75 = new TGraphErrors (nBins,Phi,Eff_MVA,0,ErrSig);
  if(variable == "numPV")GraphEff_MVA2_WP75 = new TGraphErrors (nBins,numPV,Eff_MVA,0,ErrSig);
  GraphEff_MVA2_WP75->SetMarkerColor(kRed); 
  GraphEff_MVA2_WP75->SetMarkerStyle(21); 
  GraphEff_MVA2_WP75->SetMarkerSize(1); 
  GraphEff_MVA2_WP75->SetLineColor(kRed);
  GraphEff_MVA2_WP75->SetFillColor(0);
  GraphEff_MVA2_WP75->SetLineWidth(1);
  GraphEff_MVA2_WP75->SetTitle("WP75");

  TGraphErrors* GraphEff_MVA2_WP85 = new TGraphErrors (nBins,Pt,Eff_MVA2_WP85,0,ErrSig);
  if(variable == "AbsEta")GraphEff_MVA2_WP85 = new TGraphErrors (nBins,AbsEta,Eff_MVA,0,ErrSig);
  if(variable == "Phi")GraphEff_MVA2_WP85 = new TGraphErrors (nBins,Phi,Eff_MVA,0,ErrSig);
  if(variable == "numPV")GraphEff_MVA2_WP85 = new TGraphErrors (nBins,numPV,Eff_MVA,0,ErrSig);
  GraphEff_MVA2_WP85->SetMarkerColor(kGreen); 
  GraphEff_MVA2_WP85->SetMarkerStyle(21); 
  GraphEff_MVA2_WP85->SetMarkerSize(1); 
  GraphEff_MVA2_WP85->SetLineColor(kGreen);
  GraphEff_MVA2_WP85->SetFillColor(0);
  GraphEff_MVA2_WP85->SetLineWidth(1);
  GraphEff_MVA2_WP85->SetTitle("WP85");

  TGraphErrors* GraphEff_MVA2_WP95 = new TGraphErrors (nBins,Pt,Eff_MVA2_WP95,0,ErrSig);
  if(variable == "AbsEta")GraphEff_MVA2_WP95 = new TGraphErrors (nBins,AbsEta,Eff_MVA,0,ErrSig);
  if(variable == "Phi")GraphEff_MVA2_WP95 = new TGraphErrors (nBins,Phi,Eff_MVA,0,ErrSig);
  if(variable == "numPV")GraphEff_MVA2_WP95 = new TGraphErrors (nBins,numPV,Eff_MVA,0,ErrSig);
  GraphEff_MVA2_WP95->SetMarkerColor(kMagenta); 
  GraphEff_MVA2_WP95->SetMarkerStyle(21); 
  GraphEff_MVA2_WP95->SetMarkerSize(1); 
  GraphEff_MVA2_WP95->SetLineColor(kMagenta);
  GraphEff_MVA2_WP95->SetFillColor(0);
  GraphEff_MVA2_WP95->SetLineWidth(1);
  GraphEff_MVA2_WP95->SetTitle("WP95");

  TMultiGraph *mg_Eff = new TMultiGraph();
  mg_Eff->Add(GraphEff_MVA);
  mg_Eff->Add(GraphEff_MVA2_WP75);
  mg_Eff->Add(GraphEff_MVA2_WP85);
  mg_Eff->Add(GraphEff_MVA2_WP95);

  TLegend* leg_Eff = new TLegend(0.6,0.75,0.8,0.88,NULL,"brNDC");
  leg_Eff->SetFillStyle(0);
  leg_Eff->SetBorderSize(0);
  leg_Eff->SetFillColor(10);
  leg_Eff->SetTextSize(0.03);

  leg_Eff->AddEntry(GraphEff_MVA,"old MVA");
  leg_Eff->AddEntry(GraphEff_MVA2_WP75,"new MVA WP75");
  leg_Eff->AddEntry(GraphEff_MVA2_WP85,"new MVA WP85");
  leg_Eff->AddEntry(GraphEff_MVA2_WP95,"new MVA WP95");

 ///FakeRates
  TGraphErrors* GraphFake_MVA = new TGraphErrors (nBins,Pt,Fake_MVA,0,ErrBkg);
  if(variable == "AbsEta")GraphFake_MVA = new TGraphErrors (nBins,AbsEta,Eff_MVA,0,ErrBkg);
  if(variable == "Phi")GraphFake_MVA = new TGraphErrors (nBins,Phi,Eff_MVA,0,ErrBkg);
  if(variable == "numPV")GraphFake_MVA = new TGraphErrors (nBins,numPV,Eff_MVA,0,ErrBkg);
  GraphFake_MVA->SetMarkerColor(kBlue); 
  GraphFake_MVA->SetMarkerStyle(21); 
  GraphFake_MVA->SetMarkerSize(1); 
  GraphFake_MVA->SetLineColor(kBlue);
  GraphFake_MVA->SetFillColor(0);
  GraphFake_MVA->SetLineWidth(1);
  GraphFake_MVA->SetTitle("");

  TGraphErrors* GraphFake_MVA2_WP75 = new TGraphErrors (nBins,Pt,Fake_MVA2_WP75,0,ErrBkg );
  if(variable == "AbsEta")GraphFake_MVA2_WP75 = new TGraphErrors (nBins,AbsEta,Eff_MVA,0,ErrBkg);
  if(variable == "Phi")GraphFake_MVA2_WP75 = new TGraphErrors (nBins,Phi,Eff_MVA,0,ErrBkg);
  if(variable == "numPV")GraphFake_MVA2_WP75 = new TGraphErrors (nBins,numPV,Eff_MVA,0,ErrBkg);
  GraphFake_MVA2_WP75->SetMarkerColor(kRed); 
  GraphFake_MVA2_WP75->SetMarkerStyle(21); 
  GraphFake_MVA2_WP75->SetMarkerSize(1); 
  GraphFake_MVA2_WP75->SetLineColor(kRed);
  GraphFake_MVA2_WP75->SetFillColor(0);
  GraphFake_MVA2_WP75->SetLineWidth(1);
  GraphFake_MVA2_WP75->SetTitle("");

  TGraphErrors* GraphFake_MVA2_WP85 = new TGraphErrors (nBins,Pt,Fake_MVA2_WP85,0,ErrBkg );
  if(variable == "AbsEta")GraphFake_MVA2_WP85 = new TGraphErrors (nBins,AbsEta,Eff_MVA,0,ErrBkg);
  if(variable == "Phi")GraphFake_MVA2_WP85 = new TGraphErrors (nBins,Phi,Eff_MVA,0,ErrBkg);
  if(variable == "numPV")GraphFake_MVA2_WP85 = new TGraphErrors (nBins,numPV,Eff_MVA,0,ErrBkg);
  GraphFake_MVA2_WP85->SetMarkerColor(kGreen); 
  GraphFake_MVA2_WP85->SetMarkerStyle(21); 
  GraphFake_MVA2_WP85->SetMarkerSize(1); 
  GraphFake_MVA2_WP85->SetLineColor(kGreen);
  GraphFake_MVA2_WP85->SetFillColor(0);
  GraphFake_MVA2_WP85->SetLineWidth(1);
  GraphFake_MVA2_WP85->SetTitle("");

  TGraphErrors* GraphFake_MVA2_WP95 = new TGraphErrors (nBins,Pt,Fake_MVA2_WP95,0,ErrBkg );
  if(variable == "AbsEta")GraphFake_MVA2_WP95 = new TGraphErrors (nBins,AbsEta,Eff_MVA,0,ErrBkg);
  if(variable == "Phi")GraphFake_MVA2_WP95 = new TGraphErrors (nBins,Phi,Eff_MVA,0,ErrBkg);
  if(variable == "numPV")GraphFake_MVA2_WP95 = new TGraphErrors (nBins,numPV,Eff_MVA,0,ErrBkg);
  GraphFake_MVA2_WP95->SetMarkerColor(kMagenta); 
  GraphFake_MVA2_WP95->SetMarkerStyle(21); 
  GraphFake_MVA2_WP95->SetMarkerSize(1); 
  GraphFake_MVA2_WP95->SetLineColor(kMagenta);
  GraphFake_MVA2_WP95->SetFillColor(0);
  GraphFake_MVA2_WP95->SetLineWidth(1);
  GraphFake_MVA2_WP95->SetTitle("");

  TMultiGraph *mg_Fake = new TMultiGraph();
  mg_Fake->Add(GraphFake_MVA);
  mg_Fake->Add(GraphFake_MVA2_WP75);
  mg_Fake->Add(GraphFake_MVA2_WP85);
  mg_Fake->Add(GraphFake_MVA2_WP95);

  TLegend* leg_Fake = new TLegend(0.6,0.75,0.8,0.88,NULL,"brNDC");
  leg_Fake->SetFillStyle(0);
  leg_Fake->SetBorderSize(0);
  leg_Fake->SetFillColor(10);
  leg_Fake->SetTextSize(0.03);

  leg_Fake->AddEntry(GraphFake_MVA,"old MVA");
  leg_Fake->AddEntry(GraphFake_MVA2_WP75,"new MVA WP75");
  leg_Fake->AddEntry(GraphFake_MVA2_WP85,"new MVA WP85");
  leg_Fake->AddEntry(GraphFake_MVA2_WP95,"new MVA WP95");

  c1->Clear();

  if(type=="Eff" ){
    mg_Eff->Draw("AP");
    mg_Eff->GetXaxis()->SetTitle(xAxisTitle);
    mg_Eff->GetXaxis()->SetLimits(xMin,xMax);
    mg_Eff->GetYaxis()->SetTitle("Efficiency");
    mg_Eff->SetMinimum(0.);
    mg_Eff->SetMaximum(1.);
    leg_Eff->Draw("same");
  }
  if(type=="FakeRate"){
    mg_Fake->Draw("AP");
    mg_Fake->GetXaxis()->SetTitle(xAxisTitle);
    mg_Fake->GetXaxis()->SetLimits(xMin,xMax);
    mg_Fake->GetYaxis()->SetTitle("Fake Rate");
    mg_Fake->GetYaxis()->SetLimits(YMin-0.1*YMin,1);
    mg_Fake->SetMinimum(YMin-0.1*YMin);
    mg_Fake->SetMaximum(1.0);
    leg_Fake->Draw("same");
    c1->SetLogy();

  }


  string outputName = Form("plots/plotEfficiencies_%s_%sVs%s",category.data(),type.data(),variable.data());
  c1->Print(std::string(outputName).append(".png").data());
  c1->Print(std::string(outputName).append(".pdf").data());

  inputFile->Close();
}



void plotAllEfficiencies(){


  std::vector<std::string> categories;
  categories.push_back(std::string("All"));
  categories.push_back(std::string("NoEleMatch"));
  categories.push_back(std::string("woG"));
  categories.push_back(std::string("wGwoGSF"));
  categories.push_back(std::string("wGwGSFwoPFMVA"));
  categories.push_back(std::string("wGwGSFwPFMVA"));


  std::vector<std::string> variables;
  variables.push_back(std::string("Pt"));
//   variables.push_back(std::string("AbsEta"));
//   variables.push_back(std::string("Phi"));
//   variables.push_back(std::string("numPV"));



  std::map<std::string, std::string> xAxisTitles;
  xAxisTitles["Pt"]         = "P_{T}(#tau)" ;
  xAxisTitles["Eta"]        = "|#eta|(#tau)" ;
  xAxisTitles["Phi"]        = "#phi(#tau)" ;
  xAxisTitles["numPV"]      = "# Reco Vertices" ;



  std::map<std::string, float> xMinValues;
  xMinValues["Pt"]             = 15.;
  xMinValues["AbsEta"]         = 0.;
  xMinValues["Phi"]            = -3.;
  xMinValues["numPV"]          = 0.;



  std::map<std::string, float> xMaxValues;
  xMaxValues["Pt"]             = 150.;
  xMaxValues["AbsEta"]         = 2.3;
  xMaxValues["Phi"]            = 3.;
  xMaxValues["numPV"]          = 50.;




  std::map<std::string, int> nBins;
//   nBins["Pt"]              = 15;
//   nBins["AbsEta"]          = 30;
//   nBins["Phi"]             = 30;
//   nBins["numPV"]           = 50;
  nBins["Pt"]              = 1;
  nBins["AbsEta"]          = 1;
  nBins["Phi"]             = 1;
  nBins["numPV"]           = 1;


 
  for ( std::vector<std::string>::const_iterator category = categories.begin();
	category != categories.end(); ++category ) {
    for ( std::vector<std::string>::const_iterator variable = variables.begin();
	  variable!= variables.end(); ++variable ) {
//       plotEfficiency("Eff",variable->data(), category->data(), xAxisTitles[*variable].data(),xMinValues[*variable], xMaxValues[*variable], nBins[*variable]);
      plotEfficiency("FakeRate",variable->data(), category->data(), xAxisTitles[*variable].data(),xMinValues[*variable], xMaxValues[*variable], nBins[*variable]);
    }
  }
}





























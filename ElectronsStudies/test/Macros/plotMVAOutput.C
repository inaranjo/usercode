#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TGraphAsymmErrors.h>
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



void plotTestMVAOutput(string Region = "Endcap"
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

  TLegend* leg = new TLegend(0.65,0.42,0.95,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  leg->SetHeader(Form("#splitline{CMS Preliminary}{#splitline{ #sqrt{s}=7 TeV}{%s}}",Region.data()));

  std::string inFileName1 = "AntiEMVA_Pythia.root";
  cout<<"opening file : "<<inFileName1<<endl;
  TFile* inFile1 = new TFile (inFileName1.data(),"READ");
  if(inFile1->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
  TTree* t1 = (TTree*)inFile1->Get("AntiEMVAAnalyzer/tree");

  TCut RegionCut ="";
  if(Region=="Endcap") RegionCut = RegionCut && TCut("Tau_AbsEta>=1.479 || Elec_AbsEta>=1.479");
  if(Region=="Barrel") RegionCut = RegionCut && TCut("Tau_AbsEta<1.479 && Elec_AbsEta<1.479");


  TH1F* hMVA1   = new TH1F( "hMVA1","MVA Output" ,100,-1,1);
  hMVA1->SetLineColor(kBlue);
  hMVA1->SetLineWidth(2);
  hMVA1->SetTitle("");
  t1->Draw("Tau_AntiEMVA2Out>>hMVA1","Tau_GenHadMatch && Elec_GenHadMatch"*RegionCut);

  TH1F* hMVA2 = new TH1F( "hMVA2","MVA Output" ,100,-1,1);
  hMVA2->SetLineColor(kRed);
  hMVA2->SetLineWidth(2);
  hMVA2->SetTitle("");
  t1->Draw("Tau_AntiEMVA2Out>>hMVA2","Tau_GenEleMatch && Elec_GenEleMatch"*RegionCut);

  TAxis* xAxis = hMVA1->GetXaxis();
  xAxis->SetTitle("AntiEMVA2 output");
  xAxis->SetTitleOffset(1.15);
  TAxis* yAxis = hMVA1->GetYaxis();
  yAxis->SetTitle("a.u.");
  yAxis->SetTitleOffset(1.30);

  hMVA1->Scale(1./hMVA1->Integral());
  hMVA2->Scale(1./hMVA2->Integral());

  hMVA1->Draw();
  hMVA2->Draw("same");

  leg->AddEntry(hMVA1,"Signal");
  leg->AddEntry(hMVA2,"Background");
  leg->Draw("same");

  string outputName2 = Form("MVAOutput_PythiaMC_%s",Region.data());
  c1->Print(std::string(outputName2).append(".png").data());
  c1->Print(std::string(outputName2).append(".pdf").data());

}

void plotMVAOutput(string Discr_ = "",
		   string Region = "Barrel"
		   )
{
  std::vector<std::string> categories;
  categories.push_back(std::string("NoEleMatch"));
  categories.push_back(std::string("woG"));
  categories.push_back(std::string("wGwoGSF"));
  categories.push_back(std::string("wGwGSFwoPFMVA"));
  if (Discr_ == "")categories.push_back(std::string("wGwGSFwPFMVA"));

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

  TLegend* leg = new TLegend(0.65,0.42,0.95,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);

  TFile *fSignal = new TFile(Form("/data_CMS/cms/ivo/AntiEMVA/Trees/MVAOutput/MvaOutput_AntiEMVA-v5%s_Signal.root",Discr_.data()),"READ"); 
  if(fSignal->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
  TFile *fBackgrd = new TFile(Form("/data_CMS/cms/ivo/AntiEMVA/Trees/MVAOutput/MvaOutput_AntiEMVA-v5%s_Backgrd.root",Discr_.data()),"READ"); 
  if(fBackgrd->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
  TTree*inTreeSgn = (TTree*)fSignal->Get("tree");
  TTree*inTreeBkg = (TTree*)fBackgrd->Get("tree");

  TH1F* hMVASig   = new TH1F( "hMVASig","MVA Output" ,100,-1,1);
  hMVASig->SetLineColor(kBlue);
  hMVASig->SetLineWidth(2);
  hMVASig->SetTitle("");
  TAxis* xAxis = hMVASig->GetXaxis();
  xAxis->SetTitle("AntiEMVA2 output");
  xAxis->SetTitleOffset(1.15);

  TH1F* hMVABkg = new TH1F( "hMVABkg","MVA Output" ,100,-1,1);
  hMVABkg->SetLineColor(kRed);
  hMVABkg->SetLineWidth(2);
  hMVABkg->SetTitle("");

  leg->AddEntry(hMVASig,"Signal");
  leg->AddEntry(hMVABkg,"Background");
  
  for ( std::vector<std::string>::const_iterator category = categories.begin();
	category != categories.end(); ++category ) {
    inTreeSgn->Draw(Form("%s_%s>>hMVASig",category->data(),Region.data()));
    inTreeBkg->Draw(Form("%s_%s>>hMVABkg",category->data(),Region.data()));
    
    hMVASig->Scale(1./hMVASig->Integral());
    hMVABkg->Scale(1./hMVABkg->Integral());
    
    TAxis* yAxis = hMVASig->GetYaxis();
    yAxis->SetTitle("a.u.");
    yAxis->SetTitleOffset(1.30);
    yAxis->SetRangeUser(0,TMath::Max( hMVASig->GetMaximum(), hMVABkg->GetMaximum())+0.01);
    
    hMVASig->Draw();
    hMVABkg->Draw("SAME");

    leg->SetHeader(Form("#splitline{CMS Preliminary}{#splitline{#sqrt{s}=7 TeV}{%s %s}}",Region.data(),category->data()));
    leg->Draw();

    string outputName2 = Form("plots/plotMVAOutput_v5%s_%s_%s",Discr_.data(),category->data(),Region.data());
    c1->Print(std::string(outputName2).append(".png").data());
    c1->Print(std::string(outputName2).append(".pdf").data());
  }
  fSignal->Close();
  fBackgrd->Close();
}

void plotAllMVAOutput(){
  plotMVAOutput("","Barrel");
  plotMVAOutput("","Endcap");
//   plotMVAOutput("-AntiEMed","Barrel");
//   plotMVAOutput("-AntiEMed","Endcap");
}

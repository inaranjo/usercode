//--------------------------------------------------------------------------------------------------
//  
//
// PlotROCCurves
//
// Macro plotting efficiencies and ROC Curves for the AntiElectron MVA
//
// Authors: I.Naranjo
//--------------------------------------------------------------------------------------------------

#include "TTree.h"
#include "TBranch.h"
#include "TString.h"
#include "TPad.h"
#include <TFile.h>
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

#define DEBUG true


void plotCombinedROC(string discr = "",
	     string Region = "Endcap"
	     )
{
  std::vector<std::string> categories;
  categories.push_back(std::string("All"));
  categories.push_back(std::string("woG"));
  categories.push_back(std::string("wGwoGSF"));
  categories.push_back(std::string("wGwGSFwoPFMVA"));
  categories.push_back(std::string("wGwGSFwPFMVA"));

  float CatProbSig [5];
  float CatProbBkg [5];
  float nSig [5];
  float nBkg [5];

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
  //leg->SetHeader("#splitline{CMS Preliminary}{ #sqrt{s}=7 TeV}");

  std::string inFileName = Form("./tmva/tmvaRoot/TMVA%s_woG_%s.root",discr.data(),Region.data());
  cout<<"opening file : "<<inFileName<<endl;
  TFile* inFile = new TFile (inFileName.data(),"READ");
  if(inFile->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
  TH1F* hSigEff = (TH1F*)inFile->Get("Method_BDT/BDT/MVA_BDT_effS");
  hSigEff->SetLineColor(kBlue);
  hSigEff->SetLineWidth(2);
  hSigEff->SetTitle("");

  TH1F* hBkgEff = (TH1F*)inFile->Get("Method_BDT/BDT/MVA_BDT_effB");
  hBkgEff->SetLineColor(kRed);
  hBkgEff->SetLineWidth(2);
  hBkgEff->SetTitle("");



  int i = 0;
  for ( std::vector<std::string>::const_iterator category = categories.begin();
	category != categories.end(); ++category ) {
    TFile* fSig = new TFile (Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA%s_%s_Tau.root",discr.data(),category->data())) ;
    TTree* tSig = (TTree*)fSig->Get("tree");
    nSig[i] = tSig->GetEntries();
    TFile* fBkg = new TFile (Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA%s_%s_Elec.root",discr.data(),category->data())) ;
    TTree* tBkg = (TTree*)fBkg->Get("tree");
    nBkg[i] = tBkg->GetEntries();
    CatProbSig[i]= nSig[i]/nSig[0];
    CatProbBkg[i]= nBkg[i]/nBkg[0];
    cout <<"Signal Entries : "<<nSig[i]<<" Fraction of category:  "<<CatProbSig[i]<<endl;
    cout <<"Background Entries : "<<nBkg[i]<<" Fraction of category:  "<<CatProbBkg[i]<<endl;
    
    if(i==1){
      hSigEff->Scale(CatProbSig[i]);
      hBkgEff->Scale(CatProbBkg[i]);
    }
    if (i>1) {
    std::string inputFileName = Form("./tmva/tmvaRoot/TMVA%s_%s_%s.root",discr.data(),category->data(),Region.data());
    cout<<"opening file : "<<inputFileName<<endl;
    TFile* inputFile = new TFile (inputFileName.data(),"READ");
    if(inputFile->IsZombie()){
      cout << "No such file!" << endl;
      return;
    }
    
    TH1F* hSigEffTemp = (TH1F*)inputFile->Get("Method_BDT/BDT/MVA_BDT_effS");
    TH1F* hBkgEffTemp = (TH1F*)inputFile->Get("Method_BDT/BDT/MVA_BDT_effB");
    hSigEff->Add(hSigEffTemp,CatProbSig[i]);
    hBkgEff->Add(hBkgEffTemp,CatProbBkg[i]);


    }

    i++;			   
  }

  TAxis* xAxisS = hSigEff->GetXaxis();
  xAxisS->SetTitle("BDT");
  xAxisS->SetTitleOffset(1.15);
  TAxis* yAxisS = hSigEff->GetYaxis();
  yAxisS->SetTitle("Efficiency");
  yAxisS->SetTitleOffset(1.30);
  TAxis* xAxisB = hBkgEff->GetXaxis();
  xAxisB->SetTitle("BDT");
  xAxisB->SetTitleOffset(1.15);
  TAxis* yAxisB = hBkgEff->GetYaxis();
  yAxisB->SetTitle("Efficiency");
  yAxisB->SetTitleOffset(1.30);

  hSigEff->Draw();
  hBkgEff->Draw("same");

  leg->AddEntry(hSigEff,"Signal");
  leg->AddEntry(hBkgEff,"Background");
  leg->Draw();

  string outputName = Form("plots/plotROCCurves%s_%s",discr.data(),Region.data());
  c1->Print(std::string(outputName).append(".png").data());
  c1->Print(std::string(outputName).append(".pdf").data());
}


void plotAllROC(){
//   plotROC( "", "Barrel");
//   plotROC( "", "Endcap");
//   plotROC( "-AntiEMed", "Barrel");
//   plotROC( "-AntiEMed", "Endcap");
}


void plot2ROC(string category = "woG",
	     string Region = "Endcap"
	     )
{
  std::vector<std::string> categories;
  categories.push_back(std::string("All"));
  categories.push_back(std::string("woG"));
  categories.push_back(std::string("wGwoGSF"));
  categories.push_back(std::string("wGwGSFwoPFMVA"));
  categories.push_back(std::string("wGwGSFwPFMVA"));

  float CatProbSig [5];
  float CatProbBkg [5];
  float nSig [5];
  float nBkg [5];

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
  //leg->SetHeader("#splitline{CMS Preliminary}{ #sqrt{s}=7 TeV}");

  std::string inFileName1 = Form("./tmva/tmvaRoot/TMVA_%s_%s.root",category.data(),Region.data());
  cout<<"opening file : "<<inFileName1<<endl;
  TFile* inFile1 = new TFile (inFileName1.data(),"READ");
  if(inFile1->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
  TH1F* hROC1 = (TH1F*)inFile1->Get("Method_BDT/BDT/MVA_BDT_effBvsS");
  hROC1->SetLineColor(kBlue);
  hROC1->SetLineWidth(2);
  hROC1->SetTitle("");

  std::string inFileName2 = Form("./tmva/tmvaRoot/TMVA-AntiEMed_%s_%s.root",category.data(),Region.data());
  cout<<"opening file : "<<inFileName2<<endl;
  TFile* inFile2 = new TFile (inFileName2.data(),"READ");
  if(inFile2->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
  TH1F* hROC2 = (TH1F*)inFile2->Get("Method_BDT/BDT/MVA_BDT_effBvsS");
  hROC2->SetLineColor(kRed);
  hROC2->SetLineWidth(2);
  hROC2->SetTitle("");

  hROC1->Draw();
  hROC2->Draw("same");

  leg->AddEntry(hROC1,"No discr");
  leg->AddEntry(hROC2,"AntiEMed discr");
  leg->Draw();

}

void computeDiscriminator(string discr = "",
	     string Region = "Endcap"
	     )
{
  std::vector<std::string> categories;
  categories.push_back(std::string("All"));
  categories.push_back(std::string("woG"));
  categories.push_back(std::string("wGwoGSF"));
  categories.push_back(std::string("wGwGSFwoPFMVA"));
  categories.push_back(std::string("wGwGSFwPFMVA"));

  float CatProbSig [5];
  float CatProbBkg [5];
  float nSig [5];
  float nBkg [5];
  float nSigDiscr [5];
  float nBkgDiscr [5];
  float EffSigDiscr [5];
  float EffBkgDiscr [5];

  for (int k=0; k<5; k++){
    CatProbSig [k] = 1;
    CatProbBkg [k] = 1;
    nSig [k] = 1;
    nBkg [k] = 1;
    nSigDiscr [k] = 1;
    nBkgDiscr [k] = 1;
    EffSigDiscr [k] = 1;
    EffBkgDiscr [k] = 1;
  }

  //////////////Efficiency to pass AntiEMed//////////////
  if(discr =="-AntiEMed"){
    cout<<endl;
    cout<<"Calculating efficiencies to pas AntiEMedium discriminator..."<<endl; 
    int j = 0;
    for ( std::vector<std::string>::const_iterator category = categories.begin();
	category != categories.end(); ++category ) {
      std::string fSigName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA_%s_Tau.root",category->data());
      TFile* fSig = new TFile (fSigName.data(),"READ") ;
      TTree* tSig = (TTree*)fSig->Get("tree");
      nSig[j] = tSig->GetEntries();
      std::string fBkgName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA_%s_Elec.root",category->data());
      TFile* fBkg = new TFile (fBkgName.data(),"READ") ;
      TTree* tBkg = (TTree*)fBkg->Get("tree");
      nBkg[j] = tBkg->GetEntries();
      std::string fSigDiscrName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA%s_%s_Tau.root",discr.data(),category->data());
      TFile* fSigDiscr = new TFile (fSigDiscrName.data(),"READ") ;
      TTree* tSigDiscr = (TTree*)fSigDiscr->Get("tree");
      nSigDiscr[j] = tSigDiscr->GetEntries();
      std::string fBkgDiscrName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA%s_%s_Elec.root",discr.data(),category->data());
      TFile* fBkgDiscr = new TFile (fBkgDiscrName.data(),"READ") ;
      TTree* tBkgDiscr = (TTree*)fBkgDiscr->Get("tree");
      nBkgDiscr[j] = tBkgDiscr->GetEntries();
      EffSigDiscr[j]= nSigDiscr[j]/nSig[j];
      EffBkgDiscr[j]= nBkgDiscr[j]/nBkg[j];
      if (DEBUG){
      cout <<"Signal Entries : "<<nSigDiscr[j]<<Form(" Efficiency to pass AntiEMed for category: %s  ",category->data())<<EffSigDiscr[j]<<endl;
      cout <<"Background Entries : "<<nBkgDiscr[j]<<Form(" Efficiency to pass AntiEMed for category: %s  ",category->data())<<EffBkgDiscr[j]<<endl;
      }
      fSig->Close();
      fBkg->Close();    
      fSigDiscr->Close();    
      fBkgDiscr->Close(); 
      j++;			   
    }
  }


  ////////////////Probability of different categories//////////////
  int i = 0;
  for ( std::vector<std::string>::const_iterator category = categories.begin();
	category != categories.end(); ++category ) {
    std::string fSigName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA%s_%s_Tau.root",discr.data(),category->data());
    cout<<endl;
    cout<<"opening file for signal : "<<fSigName<<endl;
    TFile* fSig = new TFile (fSigName.data(),"READ") ;
    TTree* tSig = (TTree*)fSig->Get("tree");
    nSig[i] = tSig->GetEntries();
    std::string fBkgName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA%s_%s_Elec.root",discr.data(),category->data());
    cout<<"opening file for background : "<<fBkgName<<endl;
    TFile* fBkg = new TFile (fBkgName.data(),"READ") ;
    TTree* tBkg = (TTree*)fBkg->Get("tree");
    nBkg[i] = tBkg->GetEntries();
    CatProbSig[i]= nSig[i]/nSig[0];
    CatProbBkg[i]= nBkg[i]/nBkg[0];
    if(DEBUG){
    cout <<"Signal Entries : "<<nSig[i]<<Form(" Fraction of category: %s  ",category->data())<<CatProbSig[i]<<endl;
    cout <<"Background Entries : "<<nBkg[i]<<Form(" Fraction of category: %s  ",category->data())<<CatProbBkg[i]<<endl;
    }
    fSig->Close();
    fBkg->Close(); 
    i++;
  }

}


void computeAll(){
  computeDiscriminator( "", "Barrel");
  computeDiscriminator( "", "Endcap");
  computeDiscriminator( "-AntiEMed", "Barrel");
  computeDiscriminator( "-AntiEMed", "Endcap");
}

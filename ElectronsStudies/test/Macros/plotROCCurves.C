//--------------------------------------------------------------------------------------------------
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
#include <TMarker.h>
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



void plotEfficiencies(
	     )
{
  string discr = "";
  std::vector<std::string> discrs;
  discrs.push_back(std::string(""));

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

  TLegend* leg = new TLegend(0.45,0.12,0.55,0.45,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  //leg->SetHeader("#splitline{CMS Preliminary}{ #sqrt{s}=7 TeV}");

  std::string inFileName1 = Form("./tmva/tmvaRoot/TMVAOptimization%s.root",discrs[0].data());
  cout<<"opening file : "<<inFileName1<<endl;
  TFile* inFile1 = new TFile (inFileName1.data(),"READ");
  if(inFile1->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
  TH1F* hEffS = (TH1F*)inFile1->Get("Method_Cuts/Cuts/MVA_Cuts_effS");
  hEffS->SetLineColor(kBlue);
  hEffS->SetLineWidth(2);
  hEffS->SetTitle("");
  TAxis* xAxis = hEffS->GetXaxis();
  xAxis->SetTitle("Cut");
  xAxis->SetTitleOffset(1.15);
  xAxis->SetRangeUser(0.,1.0);
  TAxis* yAxis = hEffS->GetYaxis();
  yAxis->SetTitle("Efficiency");
  yAxis->SetTitleOffset(1.30);
  yAxis->SetRangeUser(0.,1.0);
  TH1F* hEffB = (TH1F*)inFile1->Get("Method_Cuts/Cuts/MVA_Cuts_effB");
  hEffB->SetLineColor(kRed);
  hEffB->SetLineWidth(2);
  hEffB->SetTitle("");

  hEffS->Draw();
  hEffB->Draw("same");

  leg->AddEntry(hEffS,"Signal");
  leg->AddEntry(hEffB,"Background");
  leg->Draw();

  string outputName = Form("plots/plotROCCurves_Efficiencies%s",discr.data());
  c1->Print(std::string(outputName).append(".png").data());
  c1->Print(std::string(outputName).append(".pdf").data());
}



void plotCombinedROC(string discr = "",
		     string Region = "Barrel"
		     )
{
  std::vector<std::string> categories;
  if (discr == ""){
    categories.push_back(std::string("All"));
    categories.push_back(std::string("woG"));
    categories.push_back(std::string("wGwoGSF"));
    categories.push_back(std::string("wGwGSFwoPFMVA"));
    categories.push_back(std::string("wGwGSFwPFMVA"));
  }
  if (discr == "-AntiEMed"){
    categories.push_back(std::string("All"));
    categories.push_back(std::string("woG"));
    categories.push_back(std::string("wGwoGSF"));
    categories.push_back(std::string("wGwGSFwoPFMVA"));
  }

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
//       if (DEBUG){
//       cout <<"Signal Entries : "<<nSigDiscr[j]<<Form(" Efficiency to pass AntiEMed for category: %s  ",category->data())<<EffSigDiscr[j]<<endl;
//       cout <<"Background Entries : "<<nBkgDiscr[j]<<Form(" Efficiency to pass AntiEMed for category: %s  ",category->data())<<EffBkgDiscr[j]<<endl;
//       }
      fSig->Close();
      fBkg->Close();    
      fSigDiscr->Close();    
      fBkgDiscr->Close(); 
      j++;			   
    }
  }
  //////////////Efficiency to pass AntiEMed//////////////


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
    if(DEBUG){
//     cout <<"Signal Entries : "<<nSig[i]<<" Fraction of category:  "<<CatProbSig[i]<<endl;
//     cout <<"Background Entries : "<<nBkg[i]<<" Fraction of category:  "<<CatProbBkg[i]<<endl;
    }
    if(i==1){
      hSigEff->Scale(CatProbSig[i]*EffSigDiscr[i]);
      hBkgEff->Scale(CatProbBkg[i]*EffBkgDiscr[i]);
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
    hSigEff->Add(hSigEffTemp,CatProbSig[i]*EffSigDiscr[i]);
    hBkgEff->Add(hBkgEffTemp,CatProbBkg[i]*EffBkgDiscr[i]);
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

  string outputName = Form("plots/plotROCCurves_Efficiency%s_%s",discr.data(),Region.data());
  c1->Print(std::string(outputName).append(".png").data());
  c1->Print(std::string(outputName).append(".pdf").data());

  //Combined ROC Curve of the two previous combined efficiency curves 
  if(DEBUG){
    cout<<"Number of bins Signal  : "<<hSigEff->GetSize()-2<<endl;
    cout<<"Number of bins Background  : "<<hBkgEff->GetSize()-2<<endl;
    cout<<"Bin content Sig : "<<hSigEff->GetBinContent(5000)<<endl;
    cout<<"Bin content Bkg : "<<hBkgEff->GetBinContent(5000)<<endl;
  }
  int nBins = hSigEff->GetSize();
  TH2F* hROC = new TH2F ("hROC","ROC Curve",nBins,0,1,nBins,0,1);
  hROC->SetLineWidth(2);
  TAxis* xAxisROC = hROC->GetXaxis();
  xAxisROC->SetTitle("Signal Efficiency");
  xAxisROC->SetTitleOffset(1.15);
  TAxis* yAxisROC = hROC->GetYaxis();
  yAxisROC->SetTitle("Bkg rejection");
  yAxisROC->SetTitleOffset(1.30);

  for(int k=0;k<nBins;k++){
    hROC->Fill(hSigEff->GetBinContent(k),1-hBkgEff->GetBinContent(k));
  }
  TLegend* leg2 = new TLegend(0.65,0.42,0.95,0.85,NULL,"brNDC");
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(10);
  leg2->SetTextSize(0.03);
  leg2->AddEntry(hROC,Region.data());
//   if(discr.data()=="-AntiEMed")leg2->AddEntry(hROC,Region.data()+" AntiEMedium applied");

  hROC->Draw();
  leg2->Draw();

  string outputName2 = Form("plots/plotROCCurves_ROCCurve%s_%s",discr.data(),Region.data());
  c1->Print(std::string(outputName2).append(".png").data());
  c1->Print(std::string(outputName2).append(".pdf").data());

}


void plot2ROC(string category = "woG",
	      string Region = "Endcap",
	      string discr = ""
	     )
{
  std::vector<std::string> categories;
  if (discr == ""){
    categories.push_back(std::string("All"));
    categories.push_back(std::string("woG"));
    categories.push_back(std::string("wGwoGSF"));
    categories.push_back(std::string("wGwGSFwoPFMVA"));
    categories.push_back(std::string("wGwGSFwPFMVA"));
  }
  if (discr == "-AntiEMed"){
    categories.push_back(std::string("All"));
    categories.push_back(std::string("woG"));
    categories.push_back(std::string("wGwoGSF"));
    categories.push_back(std::string("wGwGSFwoPFMVA"));
  }


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



void plotROCForCategories(string Region = "Barrel",
		       string discr = ""
		       )
{
  std::vector<std::string> categories;
  if (discr == ""){
    categories.push_back(std::string("NoEleMatch"));
    categories.push_back(std::string("woG"));
    categories.push_back(std::string("wGwoGSF"));
    categories.push_back(std::string("wGwGSFwoPFMVA"));
    categories.push_back(std::string("wGwGSFwPFMVA"));
  }
  if (discr == "-AntiEMed"){
    categories.push_back(std::string("NoEleMatch"));
    categories.push_back(std::string("woG"));
    categories.push_back(std::string("wGwoGSF"));
    categories.push_back(std::string("wGwGSFwoPFMVA"));
  }


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

  TLegend* leg = new TLegend(0.55,0.32,0.75,0.75,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  leg->SetHeader(Form("%s",Region.data()));

  std::string inFileName0 = Form("./tmva/tmvaRoot/TMVA_v5%s_%s_%s.root",discr.data(),categories[0].data(),Region.data());
  cout<<"opening file : "<<inFileName0<<endl;
  TFile* inFile0 = new TFile (inFileName0.data(),"READ");
  if(inFile0->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
//   TH1F* hROC0 = (TH1F*)inFile0->Get("Method_BDT/BDT/MVA_BDT_rejBvsS");
  TH1F* hROC0 = (TH1F*)inFile0->Get("Method_BDT/BDTG/MVA_BDTG_rejBvsS");
  hROC0->SetLineColor(kBlue);
  hROC0->SetLineWidth(2);
  hROC0->SetTitle("");

  std::string inFileName1 = Form("./tmva/tmvaRoot/TMVA_v5%s_%s_%s.root",discr.data(),categories[1].data(),Region.data());
  cout<<"opening file : "<<inFileName1<<endl;
  TFile* inFile1 = new TFile (inFileName1.data(),"READ");
  if(inFile1->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
//   TH1F* hROC1 = (TH1F*)inFile1->Get("Method_BDT/BDT/MVA_BDT_rejBvsS");
  TH1F* hROC1 = (TH1F*)inFile1->Get("Method_BDT/BDTG/MVA_BDTG_rejBvsS");
  hROC1->SetLineColor(kRed);
  hROC1->SetLineWidth(2);
  hROC1->SetTitle("");

  std::string inFileName2 = Form("./tmva/tmvaRoot/TMVA_v5%s_%s_%s.root",discr.data(),categories[2].data(),Region.data());
  cout<<"opening file : "<<inFileName2<<endl;
  TFile* inFile2 = new TFile (inFileName2.data(),"READ");
  if(inFile2->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
//   TH1F* hROC2 = (TH1F*)inFile2->Get("Method_BDT/BDT/MVA_BDT_rejBvsS");
  TH1F* hROC2 = (TH1F*)inFile2->Get("Method_BDT/BDTG/MVA_BDTG_rejBvsS");
  hROC2->SetLineColor(kBlack);
  hROC2->SetLineWidth(2);
  hROC2->SetTitle("");

 std::string inFileName3 = Form("./tmva/tmvaRoot/TMVA_v5%s_%s_%s.root",discr.data(),categories[3].data(),Region.data());
  cout<<"opening file : "<<inFileName3<<endl;
  TFile* inFile3 = new TFile (inFileName3.data(),"READ");
  if(inFile3->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
//   TH1F* hROC3 = (TH1F*)inFile3->Get("Method_BDT/BDT/MVA_BDT_rejBvsS");
  TH1F* hROC3 = (TH1F*)inFile3->Get("Method_BDT/BDTG/MVA_BDTG_rejBvsS");
  hROC3->SetLineColor(kMagenta);
  hROC3->SetLineWidth(2);
  hROC3->SetTitle("");

  std::string inFileName4 = Form("./tmva/tmvaRoot/TMVA_v5%s_%s_%s.root",discr.data(),categories[3].data(),Region.data());
  if (discr == "")inFileName4 = Form("./tmva/tmvaRoot/TMVA_v5%s_%s_%s.root",discr.data(),categories[4].data(),Region.data());
  cout<<"opening file : "<<inFileName3<<endl;
  TFile* inFile4 = new TFile (inFileName4.data(),"READ");
  if(inFile4->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
//   TH1F* hROC4 = (TH1F*)inFile4->Get("Method_BDT/BDT/MVA_BDT_rejBvsS");
  TH1F* hROC4 = (TH1F*)inFile4->Get("Method_BDT/BDTG/MVA_BDTG_rejBvsS");
  hROC4->SetLineColor(kGreen);
  hROC4->SetLineWidth(2);
  hROC4->SetTitle("");

  hROC0->Draw();
  hROC1->Draw("same");
  hROC2->Draw("same");
  hROC3->Draw("same");
  if(discr == "")hROC4->Draw("same");

  leg->AddEntry(hROC0,"Category 1");
  leg->AddEntry(hROC1,"Category 2");
  leg->AddEntry(hROC2,"Category 3");
  leg->AddEntry(hROC3,"Category 4");
  if(discr == "")leg->AddEntry(hROC4,"Category 5");
  leg->Draw();

  string outputName = Form("plots/plotROCCurves_AllCategories_v5%s_%s",discr.data(),Region.data());
  c1->Print(std::string(outputName).append(".png").data());
  c1->Print(std::string(outputName).append(".pdf").data());

}

void plotAllROCForCategories()
{
  plotROCForCategories("Barrel","");
  plotROCForCategories("Endcap","");
//   plotROCForCategories("Barrel","-AntiEMed");
//   plotROCForCategories("Endcap","-AntiEMed");
}



void plotFinalROC(
	     )
{
  string discr = "";
  std::vector<std::string> discrs;
  discrs.push_back(std::string(""));
  discrs.push_back(std::string("-AntiEMed"));

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

  TLegend* leg = new TLegend(0.35,0.12,0.55,0.45,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  //leg->SetHeader("#splitline{CMS Preliminary}{ #sqrt{s}=7 TeV}");
  //Put Markers for Loose, Medium and Tight AntiElectron discriminators
  float nSig;
  float nSigDiscr;
  float EffSigAntiELoose;
  float EffSigAntiEMed;
  float EffSigAntiETight;
  float EffSigAntiEMVA;
  float nBkg;
  float nBkgDiscr;
  float EffBkgAntiELoose;
  float EffBkgAntiEMed;
  float EffBkgAntiETight;
  float EffBkgAntiEMVA;
  //Compute efficiencies to pass discriminators :
  std::string fSigName = "/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA_All_Tau.root";
  TFile* fSig = new TFile (fSigName.data(),"READ") ;
  TTree* tSig = (TTree*)fSig->Get("tree");
  TH1F* hSig = new TH1F("hSig","",1,-10,10);
  tSig->Draw("Elec_Pt>>hSig");
  nSig = hSig->GetEntries();
  //AntiELoose discriminator
  TH1F* hSigDiscr = new TH1F("hSigDiscr","",1,-10,10);
  tSig->Draw("Elec_Pt>>hSigDiscr","Tau_AntiELoose>0.5");
  nSigDiscr = hSigDiscr->GetEntries();
  EffSigAntiELoose= nSigDiscr/nSig;
  //AntiEMedium discriminator
  tSig->Draw("Elec_Pt>>hSigDiscr","Tau_AntiEMedium>0.5");
  nSigDiscr = hSigDiscr->GetEntries();
  EffSigAntiEMed= nSigDiscr/nSig;
  //AntiETight discriminator
  tSig->Draw("Elec_Pt>>hSigDiscr","Tau_AntiETight>0.5");
  nSigDiscr = hSigDiscr->GetEntries();
  EffSigAntiETight= nSigDiscr/nSig;
  //AntiEMVA discriminator
  tSig->Draw("Elec_Pt>>hSigDiscr","Tau_AntiEMVA>0.5");
  nSigDiscr = hSigDiscr->GetEntries();
  EffSigAntiEMVA= nSigDiscr/nSig;

  std::string fBkgName = "/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA_All_Elec.root";
  TFile* fBkg = new TFile (fBkgName.data(),"READ") ;
  TTree* tBkg = (TTree*)fBkg->Get("tree");
  TH1F* hBkg = new TH1F("hBkg","",1,-10,10);
  tBkg->Draw("Elec_Pt>>hBkg");
  nBkg = hBkg->GetEntries();
 //AntiELoose discriminator
  TH1F* hBkgDiscr = new TH1F("hBkgDiscr","",1,-10,10);
  tBkg->Draw("Elec_Pt>>hBkgDiscr","Tau_AntiELoose>0.5");
  nBkgDiscr = hBkgDiscr->GetEntries();
  EffBkgAntiELoose= nBkgDiscr/nBkg;
  //AntiEMedium discriminator
  tBkg->Draw("Elec_Pt>>hBkgDiscr","Tau_AntiEMedium>0.5");
  nBkgDiscr = hBkgDiscr->GetEntries();
  EffBkgAntiEMed= nBkgDiscr/nBkg;
  //AntiETight discriminator
  tBkg->Draw("Elec_Pt>>hBkgDiscr","Tau_AntiETight>0.5");
  nBkgDiscr = hBkgDiscr->GetEntries();
  EffBkgAntiETight= nBkgDiscr/nBkg;
  //AntiEMVA discriminator
  tBkg->Draw("Elec_Pt>>hBkgDiscr","Tau_AntiEMVA>0.5");
  nBkgDiscr = hBkgDiscr->GetEntries();
  EffBkgAntiEMVA= nBkgDiscr/nBkg;

  TMarker* MarkAntiELoose = new TMarker (EffSigAntiELoose, 1-EffBkgAntiELoose,20);
  TMarker* MarkAntiEMedium = new TMarker (EffSigAntiEMed, 1-EffBkgAntiEMed,24);
  TMarker* MarkAntiETight = new TMarker (EffSigAntiETight, 1-EffBkgAntiETight,21);
  TMarker* MarkAntiEMVA = new TMarker (EffSigAntiEMVA, 1-EffBkgAntiEMVA,25);


  std::string inFileName1 = Form("./tmva/tmvaRoot/TMVAOptimization_v4%s.root",discrs[0].data());
  cout<<"opening file : "<<inFileName1<<endl;
  TFile* inFile1 = new TFile (inFileName1.data(),"READ");
  if(inFile1->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
  TH1F* hROC1 = (TH1F*)inFile1->Get("Method_Cuts/Cuts/MVA_Cuts_rejBvsS");
  hROC1->SetLineColor(kBlue);
  hROC1->SetLineWidth(2);
  hROC1->SetTitle("");
  TAxis* xAxis = hROC1->GetXaxis();
  xAxis->SetTitle("Signal Efficiency");
  xAxis->SetTitleOffset(1.15);
  xAxis->SetRangeUser(0.4,1.0);
  TAxis* yAxis = hROC1->GetYaxis();
  yAxis->SetTitle("Background Rejection");
  yAxis->SetTitleOffset(1.30);
  yAxis->SetRangeUser(0.3,1.0);

  /*
  std::string inFileName2 = Form("./tmva/tmvaRoot/TMVAOptimization_v4%s.root",discrs[1].data());
  cout<<"opening file : "<<inFileName2<<endl;
  TFile* inFile2 = new TFile (inFileName2.data(),"READ");
  if(inFile2->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
  //Second ROC Curve weighted by the efficiency to pass AntiEMedium discriminator
  TH1F* hSigEff = (TH1F*)inFile2->Get("Method_Cuts/Cuts/MVA_Cuts_effS");
  TH1F* hBkgEff = (TH1F*)inFile2->Get("Method_Cuts/Cuts/MVA_Cuts_effB");
  //float EffAntiESig = EffSigAntiEMed;//0.742067;
  //float EffAntiEBkg = EffBkgAntiEMed;//0.0908936;
  hSigEff->Scale(EffSigAntiEMed);
  hBkgEff->Scale(EffBkgAntiEMed);
//   hSigEff->Draw();
//   hBkgEff->Draw("same");
  int nBins = hSigEff->GetSize();
//   TH1F* hROC2 = (TH1F*)inFile2->Get("Method_Cuts/Cuts/MVA_Cuts_rejBvsS");
  TGraph* hROC2 = new TGraph (nBins-2);
  hROC2->SetLineColor(kRed);
  hROC2->SetFillColor(0);
  hROC2->SetLineWidth(2);
  hROC2->SetTitle("");

  for(int k=0;k<nBins-2;k++){
    cout<<"k "<<k<<endl;
    cout<<" X "<<hSigEff->GetBinContent(k)<<endl;
    cout<<" Y "<<1-hBkgEff->GetBinContent(k)<<endl;
    hROC2->SetPoint(k,hSigEff->GetBinContent(k),1-hBkgEff->GetBinContent(k));
  }
//   hROC2->SetPoint(nBins-2,EffSigAntiEMed,1-EffBkgAntiEMed);


  */

  hROC1->Smooth();
  hROC1->Draw();
//   hROC2->Draw("same");
  MarkAntiELoose->Draw("same");
  MarkAntiEMedium->Draw("same");
  MarkAntiETight->Draw("same");
  MarkAntiEMVA->Draw("same");

  leg->AddEntry(MarkAntiELoose,"AntiELoose");
  leg->AddEntry(MarkAntiEMedium,"AntiEMedium");
  leg->AddEntry(MarkAntiETight,"AntiETight");
  leg->AddEntry(MarkAntiEMVA,"Old AntiEMVA");
  leg->AddEntry(hROC1,"New AntiEMVA");
//   leg->AddEntry(hROC2,"AntiEMed discr");
  leg->Draw();

  string outputName = "plots/plotROCCurves_FinalROC";
  c1->Print(std::string(outputName).append(".png").data());
  c1->Print(std::string(outputName).append(".pdf").data());
}

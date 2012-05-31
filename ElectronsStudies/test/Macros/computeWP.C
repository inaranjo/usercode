//--------------------------------------------------------------------------------------------------
// computeDiscriminator
//
// Macro computing the probability to be in each category. And computing working points for given Signal/Background ratio.
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
#include <fstream>

#define DEBUG true


  ////////////////Probability of different categories//////////////
void computeCategory(Float_t out[2][10])
{
  std::vector<std::string> categories;
  categories.push_back(std::string("NoEleMatch"));
  categories.push_back(std::string("woG"));
  categories.push_back(std::string("wGwoGSF"));
  categories.push_back(std::string("wGwGSFwoPFMVA"));
  categories.push_back(std::string("wGwGSFwPFMVA"));
  
  std::vector<std::string> Regions;
  Regions.push_back(std::string("Barrel"));
  Regions.push_back(std::string("Endcap"));

  std::string fNameS = "/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA_All_Tau.root";
  TFile* fS = new TFile (fNameS.data(),"READ") ;
  TTree* tS = (TTree*)fS->Get("tree");
  int nSig = tS->GetEntries();

  std::string fNameB = "/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA_All_Elec.root";
  TFile* fB = new TFile (fNameB.data(),"READ") ;
  TTree* tB = (TTree*)fB->Get("tree");
  int nBkg = tB->GetEntries();

  fS->Close();
  fB->Close();

  int i = 0;
  for ( std::vector<std::string>::const_iterator Region = Regions.begin();
	Region != Regions.end(); ++Region ) {
    for ( std::vector<std::string>::const_iterator category = categories.begin();
	  category != categories.end(); ++category ) {

      std::string fSigName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA_%s_Tau_%s.root",category->data(),Region->data());
      TFile* fSig = new TFile (fSigName.data(),"READ") ;
      TTree* tSig = (TTree*)fSig->Get("tree");
      std::string fBkgName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA_%s_Elec_%s.root",category->data(),Region->data());
      TFile* fBkg = new TFile (fBkgName.data(),"READ") ;
      TTree* tBkg = (TTree*)fBkg->Get("tree");
      float EntriesSig = tSig->GetEntries(); 
      float EntriesBkg = tBkg->GetEntries(); 
      out[0][i] = EntriesSig/nSig;
      out[1][i] = EntriesBkg/nBkg;

      if(DEBUG){
	cout <<"Signal Entries : "<<tSig->GetEntries()<<Form(" Fraction of category: %s  ",category->data())<<out[0][i]<<endl;
	cout <<"Background Entries : "<<tBkg->GetEntries()<<Form(" Fraction of category: %s  ",category->data())<<out[1][i]<<endl;
      }
      fSig->Close();
      fBkg->Close();
      i++; 
    }
  }

  return;
}
void computeCat(){
  Float_t Prob[2][10];
  computeCategory(Prob);
}
////////////////Get cut from EffVsBDT and FRVsBDT histos//////////////
void getCut(float out [3],
	    Float_t probSig,
	    Float_t probBkg,
	    TH1* hS,
	    TH1* hB,
	    float SigOverBkg = 80
	    )
{
  float nBins = hS->GetNbinsX();
//   float nBinsB = hB->GetNbinsX();

//   cout<<"Bins S : "<<nBinsS<<endl;
//   cout<<"Bins B : "<<nBinsB<<endl;
//   cout<<"Sig prob : "<<probSig<<endl;
//   cout<<"Bkg prob : "<<probBkg<<endl;

  float WPbin = 0;
  float CutWP = -99;
  float RatioWP = 999;
  float EffWP = -99;
  float FRWP = -99;
  for ( int k = 1 ; k<=nBins ; k++){
    float Ratio = (hS->GetBinContent(k)*probSig)/(hB->GetBinContent(k)*probBkg);
    if(DEBUG)cout<<k<<"   "<<Ratio<<endl;

//Search for the WP maximizing the efficiency
    if (Ratio>SigOverBkg && (hS->GetBinContent(k)*probSig)>EffWP){
      RatioWP = Ratio ; 
      WPbin = k;
      EffWP = hS->GetBinContent(k)*probSig;
      FRWP = hB->GetBinContent(k)*probBkg;
      CutWP = hS->GetBinCenter(k);
    }
  }
  if(DEBUG){
    cout<<"Ratio chosen : "<<RatioWP<<endl;
    cout<<"Efficiency : "<<hS->GetBinContent(WPbin)<<endl;
    cout<<"Efficiency*probaCat : "<<EffWP<<endl;
    cout<<"Fake rate : "<<hB->GetBinContent(WPbin)<<endl;
    cout<<"Fake rate*probaCat : "<<FRWP<<endl;
    cout<<"bin for WP chosen: "<<WPbin<<endl;
    cout<<"Cut : "<<CutWP<<endl;
  }


  out [0] = EffWP;
  out [1] = FRWP;
  out [2] = CutWP;

}




////////////////Get WPs for all categories//////////////
double computeWP(std::ofstream& out,float SigOverBkg, double output[2])
{


  //Get probabilities to enter each one of the categories for signal and background
  Float_t Prob[2][10];
  float Eff[10];
  float FR[10];
  float Cut[10];
  for (int k=0; k<10; k++){
    Prob[0][k] = 1;
    Prob[1][k] = 1;
    Eff[k]=-99;
    FR[k]=-99;
    Cut[k]=-99;
  }
  computeCategory(Prob);

  if(DEBUG){
    for (int k=0; k<10; k++){
      cout<<Prob[0][k]<<endl;
      cout<<Prob[1][k]<<endl;
    }
  }
  //Get the EffVsBDT and FRVsBDT histograms for each category
  std::vector<std::string> categories;
  categories.push_back(std::string("NoEleMatch"));
  categories.push_back(std::string("woG"));
  categories.push_back(std::string("wGwoGSF"));
  categories.push_back(std::string("wGwGSFwoPFMVA"));
  categories.push_back(std::string("wGwGSFwPFMVA"));
  
  std::vector<std::string> Regions;
  Regions.push_back(std::string("Barrel"));
  Regions.push_back(std::string("Endcap"));

  int Iter = 0;
  for ( std::vector<std::string>::const_iterator Region = Regions.begin();
	Region != Regions.end(); ++Region ) {   
    
    for ( std::vector<std::string>::const_iterator category = categories.begin();
	  category != categories.end(); ++category ) {
      
      std::string fName = Form("./tmva/tmvaRoot/TMVA_v5_%s_%s.root",category->data(),Region->data());
      TFile* f = new TFile (fName.data(),"READ") ;
      
      TString histogramNameS ="Method_BDT/BDTG/MVA_BDTG_effS";
      TString histogramNameB ="Method_BDT/BDTG/MVA_BDTG_effB";

      TH1* histogramS = (TH1*)f->Get(histogramNameS.Data());
      TH1* histogramB = (TH1*)f->Get(histogramNameB.Data());

      float ProbaSig = Prob[0][Iter];
      float ProbaBkg = Prob[1][Iter];

      if(DEBUG)cout<<"For category : "<<category->data()<<" and Region : "<<Region->data()<<endl;
      float Result[3] ;
      getCut(Result,ProbaSig,ProbaBkg,histogramS,histogramB,SigOverBkg);
      Eff[Iter] = Result[0];
      FR[Iter] = Result[1];
      Cut[Iter] = Result[2];
      f->Close();
      Iter++;
    }
  }

  float Efficiency = 0;
  float FakeRate = 0;
  for (int k=0; k<10; k++){
    Efficiency = Efficiency + Eff[k];
    FakeRate = FakeRate + FR[k];
  }
  if(DEBUG){
    cout<<endl;
    cout<<"Efficiency : "<<Efficiency<<endl;
    cout<<"Fake Rate : "<<FakeRate<<endl;
  }
  out<<" S/B :"<<SigOverBkg<<" : Eff : "<<Efficiency<<" : FR : "<<FakeRate<<endl;
  out<<"cuts = : }"<<endl;
  for (int k=0; k<10; k++){
    out<<Cut[k]<<endl;
  }
  out<<" }"<<endl;

  output [0] = Efficiency;
  output [1] = FakeRate;

  return 0;

}



void getWP()
{
  string outName = "computeWP.txt";
  std::ofstream out(outName.c_str());
  out.precision(8);

  Double_t x[1000], Eff[1000], FR[1000];
  Int_t n = 1000;
  double Pair [2];
  for (Int_t i=1;i<=n;i++) {
    cout<<i<<endl;
    computeWP(out,i,Pair);
    x[i] = i;
    Eff[i] = Pair[0];
    FR[i] = Pair[1];
  }
  TGraph* gEff = new TGraph(n,x,Eff);
  gEff->SetLineColor(kBlue);
  TGraph* gFR = new TGraph(n,x,FR);
  gFR->SetLineColor(kRed);

  gEff->Draw("AC");
  gFR->Draw("same");

  TString outRootName = "computeWP.root";
  TFile outRootFile(outRootName,  "RECREATE");
  gEff->Write();
  gFR->Write();
  outRootFile.Close();
}

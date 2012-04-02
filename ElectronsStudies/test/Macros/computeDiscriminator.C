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

//--------------------------------------------------------------------------------------------------
// computeDiscriminator
//
// Macro computing the probability to pass AntiEMedium discriminator, and probability to be in each category
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



void computeDiscriminator1(string discr = "")
{
  std::vector<std::string> categories;
  if (discr == ""){
    categories.push_back(std::string("All"));
    categories.push_back(std::string("NoEleMatch"));
    categories.push_back(std::string("woG"));
    categories.push_back(std::string("wGwoGSF"));
    categories.push_back(std::string("wGwGSFwoPFMVA"));
    categories.push_back(std::string("wGwGSFwPFMVA"));
  }
  if (discr == "-AntiEMed"){
    categories.push_back(std::string("All"));
    categories.push_back(std::string("NoEleMatch"));
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
  cout<<endl;
  cout<<"Calculating efficiencies to pass AntiEMedium discriminator..."<<endl; 
  int j = 0;
  for ( std::vector<std::string>::const_iterator category = categories.begin();
	category != categories.end(); ++category ) {
    std::string fSigName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA_%s_Tau.root",category->data());
    TFile* fSig = new TFile (fSigName.data(),"READ") ;
    TTree* tSig = (TTree*)fSig->Get("tree");
    nSig[j] = tSig->GetEntries();

    TH1F* hSig = new TH1F("hSig","",1,-10,10);
    tSig->Draw("Elec_Pt>>hSig","Tau_AbsEta<2.3 && Tau_Pt>20");
    nSig[j] = hSig->GetEntries();

    std::string fBkgName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA_%s_Elec.root",category->data());
    TFile* fBkg = new TFile (fBkgName.data(),"READ") ;
    TTree* tBkg = (TTree*)fBkg->Get("tree");
    nBkg[j] = tBkg->GetEntries();
    std::string fSigDiscrName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA-AntiEMed_%s_Tau.root",category->data());
    TFile* fSigDiscr = new TFile (fSigDiscrName.data(),"READ") ;
    TTree* tSigDiscr = (TTree*)fSigDiscr->Get("tree");
    nSigDiscr[j] = tSigDiscr->GetEntries();

    TH1F* hSigDiscr = new TH1F("hSigDiscr","",1,-10,10);
    tSigDiscr->Draw("Elec_Pt>>hSigDiscr","Tau_AbsEta<2.3 && Tau_Pt>20");
    nSigDiscr[j] = hSigDiscr->GetEntries();

    std::string fBkgDiscrName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA-AntiEMed_%s_Elec.root",category->data());
    TFile* fBkgDiscr = new TFile (fBkgDiscrName.data(),"READ") ;
    TTree* tBkgDiscr = (TTree*)fBkgDiscr->Get("tree");
    nBkgDiscr[j] = tBkgDiscr->GetEntries();
    EffSigDiscr[j]= nSigDiscr[j]/nSig[j];
    EffBkgDiscr[j]= nBkgDiscr[j]/nBkg[j];
    if (DEBUG){
      cout <<"Signal Entries : "<<nSig[j]<<" Signal Entries after discr : "<<nSigDiscr[j]<<Form(" Efficiency to pass AntiEMed for category: %s  ",category->data())<<EffSigDiscr[j]<<endl;
      cout <<"Background Entries : "<<nBkgDiscr[j]<<Form(" Efficiency to pass AntiEMed for category: %s  ",category->data())<<EffBkgDiscr[j]<<endl;
    }
    fSig->Close();
    fBkg->Close();    
    fSigDiscr->Close();    
    fBkgDiscr->Close(); 
    j++;			   
  }


  ////////////////Probability of different categories//////////////
  cout<<endl;
  cout<<"Calculating efficiencies to be in a given category..."<<endl; 
  cout<<" No discriminator required :"<<endl; 
  int i = 0;
  for ( std::vector<std::string>::const_iterator category = categories.begin();
	category != categories.end(); ++category ) {
    std::string fSigName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA_%s_Tau.root",category->data());
    cout<<endl;
    cout<<"opening file for signal : "<<fSigName<<endl;
    TFile* fSig = new TFile (fSigName.data(),"READ") ;
    TTree* tSig = (TTree*)fSig->Get("tree");
    nSig[i] = tSig->GetEntries();
    std::string fBkgName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA_%s_Elec.root",category->data());
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
  }
  cout<<" Discriminator required :"<<endl; 
  int m = 0;
  for ( std::vector<std::string>::const_iterator category = categories.begin();
	category != categories.end(); ++category ) {
    std::string fSigName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA-AntiEMed_%s_Tau.root",category->data());
    cout<<endl;
    cout<<"opening file for signal : "<<fSigName<<endl;
    TFile* fSig = new TFile (fSigName.data(),"READ") ;
    TTree* tSig = (TTree*)fSig->Get("tree");
    nSig[i] = tSig->GetEntries();
    std::string fBkgName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA-AntiEMed_%s_Elec.root",category->data());
    cout<<"opening file for background : "<<fBkgName<<endl;
    TFile* fBkg = new TFile (fBkgName.data(),"READ") ;
    TTree* tBkg = (TTree*)fBkg->Get("tree");
    nBkg[m] = tBkg->GetEntries();
    CatProbSig[m]= nSig[m]/nSig[0];
    CatProbBkg[m]= nBkg[m]/nBkg[0];
    if(DEBUG){
      cout <<"Signal Entries : "<<nSig[i]<<Form(" Fraction of category: %s  ",category->data())<<CatProbSig[m]<<endl;
      cout <<"Background Entries : "<<nBkg[i]<<Form(" Fraction of category: %s  ",category->data())<<CatProbBkg[m]<<endl;
    }
    fSig->Close();
    fBkg->Close(); 
  }
  m++;
  return;
}
				      

void computeDiscriminator2(string discr = "")
{
  std::vector<std::string> categories;
  if (discr == ""){
    categories.push_back(std::string("All"));
    categories.push_back(std::string("NoEleMatch"));
    categories.push_back(std::string("woG"));
    categories.push_back(std::string("wGwoGSF"));
    categories.push_back(std::string("wGwGSFwoPFMVA"));
    categories.push_back(std::string("wGwGSFwPFMVA"));
  }
  if (discr == "-AntiEMed"){
    categories.push_back(std::string("All"));
    categories.push_back(std::string("NoEleMatch"));
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
  cout<<endl;
  cout<<"Calculating efficiencies to pass AntiEMedium discriminator..."<<endl; 
  int j = 0;
  for ( std::vector<std::string>::const_iterator category = categories.begin();
	category != categories.end(); ++category ) {
    std::string fSigName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA_%s_Tau.root",category->data());
    TFile* fSig = new TFile (fSigName.data(),"READ") ;
    TTree* tSig = (TTree*)fSig->Get("tree");

    TH1F* hSig = new TH1F("hSig","",1,-10,10);
    tSig->Draw("Tau_AntiEMedium>>hSig","Tau_AbsEta<2.3 && Tau_Pt>20");
    nSig[j] = hSig->GetEntries();

    TH1F* hSigDiscr = new TH1F("hSigDiscr","",1,-10,10);
    tSig->Draw("Tau_AntiEMedium>>hSigDiscr","Tau_AbsEta<2.3 && Tau_Pt>20 && Tau_AntiEMedium>0.5");
    nSigDiscr[j] = hSigDiscr->GetEntries();

    std::string fBkgName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA_%s_Elec.root",category->data());
    TFile* fBkg = new TFile (fBkgName.data(),"READ") ;
    TTree* tBkg = (TTree*)fBkg->Get("tree");

    TH1F* hBkg = new TH1F("hBkg","",1,-10,10);
    tBkg->Draw("Tau_AntiEMedium>>hBkg","Tau_AbsEta<2.3 && Tau_Pt>20");
    nBkg[j] = hBkg->GetEntries();

    TH1F* hBkgDiscr = new TH1F("hBkgDiscr","",1,-10,10);
    tBkg->Draw("Tau_AntiEMedium>>hBkgDiscr","Tau_AbsEta<2.3 && Tau_Pt>20 && Tau_AntiEMedium>0.5");
    nBkgDiscr[j] = hBkgDiscr->GetEntries();

    EffSigDiscr[j]= nSigDiscr[j]/nSig[j];
    EffBkgDiscr[j]= nBkgDiscr[j]/nBkg[j];
    if (DEBUG){
      cout <<"Signal Entries : "<<nSig[j]<<" Signal Entries after discr : "<<nSigDiscr[j]<<Form(" Efficiency to pass AntiEMed for category: %s  ",category->data())<<EffSigDiscr[j]<<endl;
      cout <<"Background Entries : "<<nBkg[j]<<" Background Entries after discr : "<<nBkgDiscr[j]<<Form(" Efficiency to pass AntiEMed for category: %s  ",category->data())<<EffBkgDiscr[j]<<endl;
    }
    fSig->Close();
    fBkg->Close();     
    j++;			   
  }


  ////////////////Probability of different categories//////////////
  cout<<endl;
  cout<<"Calculating efficiencies to be in a given category..."<<endl; 
  cout<<Form(" discriminator required : %s",discr.data())<<endl; 
  int i = 0;
  for ( std::vector<std::string>::const_iterator category = categories.begin();
	category != categories.end(); ++category ) {
    std::string fSigName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA_%s_Tau.root",category->data());
    //       cout<<"opening file for signal : "<<fSigName<<endl;
    TFile* fSig = new TFile (fSigName.data(),"READ") ;
    TTree* tSig = (TTree*)fSig->Get("tree");
    nSig[i] = tSig->GetEntries();
    std::string fBkgName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA_%s_Elec.root",category->data());
    //       cout<<"opening file for background : "<<fBkgName<<endl;
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

  return;

}


void computeAll(){
  computeDiscriminator2("");
}



void plotDiscriminators()
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

  std::vector<std::string> discriminators;
    discriminators.push_back(std::string("AntiELoose"));
    discriminators.push_back(std::string("AntiEMedium"));
    discriminators.push_back(std::string("AntiETight"));
    discriminators.push_back(std::string("AntiEMVA"));

  float nSig ;
  float nBkg ;
  float nSigDiscr ;
  float nBkgDiscr ;
  float EffSigDiscr ;
  float EffBkgDiscr ;

  //////////////Efficiency to pass AntiEMed//////////////
  cout<<endl;
  cout<<"Calculating efficiencies to pass AntiEMedium discriminator..."<<endl; 
  for ( std::vector<std::string>::const_iterator discr = discriminators.begin();
	discr != discriminators.end(); ++discr ) {
    std::string fSigName = "/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA_All_Tau.root";
    TFile* fSig = new TFile (fSigName.data(),"READ") ;
    TTree* tSig = (TTree*)fSig->Get("tree");

    TH1F* hSig = new TH1F("hSig",Form("Tau_%s",discr->data()),2,-0.5,1.5);
    tSig->Draw(Form("Tau_%s>>hSig",discr->data()),"Tau_AbsEta<2.3 && Tau_Pt>20");
    nSig = hSig->GetEntries();

    TH1F* hSigDiscr = new TH1F("hSigDiscr","",2,-0.5,1.5);
    tSig->Draw(Form("Tau_%s>>hSigDiscr",discr->data()),Form("Tau_AbsEta<2.3 && Tau_Pt>20 && Tau_%s>0.5",discr->data()));
    nSigDiscr = hSigDiscr->GetEntries();

    std::string fBkgName = "/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA_All_Elec.root";
    TFile* fBkg = new TFile (fBkgName.data(),"READ") ;
    TTree* tBkg = (TTree*)fBkg->Get("tree");

    TH1F* hBkg = new TH1F("hBkg","",2,-0.5,1.5);
    tBkg->Draw(Form("Tau_%s>>hBkg",discr->data()),"Tau_AbsEta<2.3 && Tau_Pt>20");
    nBkg = hBkg->GetEntries();

    TH1F* hBkgDiscr = new TH1F("hBkgDiscr","",2,-0.5,1.5);
    tBkg->Draw(Form("Tau_%s>>hBkgDiscr",discr->data()),Form("Tau_AbsEta<2.3 && Tau_Pt>20 && Tau_%s>0.5",discr->data()));
    nBkgDiscr = hBkgDiscr->GetEntries();

    hSig->Draw();

    EffSigDiscr= nSigDiscr/nSig;
    EffBkgDiscr= nBkgDiscr/nBkg;
    if (DEBUG){
      cout <<"Signal Entries : "<<nSig<<" Signal Entries after discr : "<<nSigDiscr<<Form(" Efficiency to pass %s for category: All  ",discr->data())<<EffSigDiscr<<endl;
      cout <<"Background Entries : "<<nBkg<<" Background Entries after discr : "<<nBkgDiscr<<Form(" Efficiency to pass %s for category: All  ",discr->data())<<EffBkgDiscr<<endl;
    }

    string outputName = Form("plots/computeDiscriminator_%s",discr->data());
    c1->Print(std::string(outputName).append(".png").data());
    c1->Print(std::string(outputName).append(".pdf").data());
  
    fSig->Close();
    fBkg->Close();     
  }
}



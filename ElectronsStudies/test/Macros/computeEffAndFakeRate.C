//--------------------------------------------------------------------------------------------------
// computeDiscriminator
//
// Macro computing the Efficiency and fake rate from the tmva optimization tesTree
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



void compute()
{
    std::string fSigName = "./tmva/tmvaRoot/V4/TMVAOptimization_v4Bis.root";
    TFile* fSig = new TFile (fSigName.data(),"READ") ;
    TTree* tSig = (TTree*)fSig->Get("TestTree");
    float ntot = tSig->GetEntries();

    TH1F* hSig = new TH1F("hSig","",1,-10,10);
    tSig->Draw("className>>hSig","classID<0.5");
    float nSig = hSig->GetEntries();

    TH1F* hBkg = new TH1F("hBkg","",1,-10,10);
    tSig->Draw("className>>hBkg","classID>0.5");
    float nBkg = hBkg->GetEntries();

    if (DEBUG){
      cout<<"Total entries :"<<ntot<<endl;
      cout<<"Signal entries :"<<nSig<<endl;
      cout<<"Bkg entries :"<<nBkg<<endl;
    }

    float NoEleMatch_Barrel_Cut = 0.0306715 ;  
    float woG_Barrel_Cut = 0.992195 ; 
    float wGwoGSF_Barrel_Cut = 0.308324 ; 
    float wGwGSFwoPFMVA_Barrel_Cut = -0.0370998 ; 
    float wGwGSFwPFMVA_Barrel_Cut = 0.864643 ; 
    float NoEleMatch_Endcap_Cut = 0.0832094 ; 
    float woG_Endcap_Cut = 0.791665 ; 
    float wGwoGSF_Endcap_Cut = 0.675537 ; 
    float wGwGSFwoPFMVA_Endcap_Cut = 0.87047 ; 
    float wGwGSFwPFMVA_Endcap_Cut = 0.233711 ; 

    TCut Cut =Form("NoEleMatch_Barrel>%f && woG_Barrel>%f && wGwoGSF_Barrel>%f && wGwGSFwoPFMVA_Barrel>%f &&wGwGSFwPFMVA_Barrel>%f &&NoEleMatch_Endcap>%f && woG_Endcap>%f && wGwoGSF_Endcap>%f && wGwGSFwoPFMVA_Endcap>%f &&wGwGSFwPFMVA_Endcap>%f",
		   NoEleMatch_Barrel_Cut,
		   woG_Barrel_Cut, 
		   wGwoGSF_Barrel_Cut, 
		   wGwGSFwoPFMVA_Barrel_Cut,
		   wGwGSFwPFMVA_Barrel_Cut,
		   NoEleMatch_Endcap_Cut, 
		   woG_Endcap_Cut, 
		   wGwoGSF_Endcap_Cut, 
		   wGwGSFwoPFMVA_Endcap_Cut,
		   wGwGSFwPFMVA_Endcap_Cut
		   );

    TH1F* hSigPass = new TH1F("hSig","",1,-10,10);
    tSig->Draw("className>>hSig","classID<0.5"*Cut);
    float nSigPass = hSigPass->GetEntries();

    TH1F* hBkgPass = new TH1F("hBkg","",1,-10,10);
    tSig->Draw("className>>hBkg","classID>0.5"*Cut);
    float nBkgPass = hBkgPass->GetEntries();

    if (DEBUG){
      cout<<"Total entries :"<<ntot<<endl;
      cout<<"Signal entries pass :"<<nSigPass<<" Efficiency -> "<<nSigPass/nSig<<endl;
      cout<<"Bkg entries pass :"<<nBkgPass<<" FakeRate -> "<<nBkgPass/nBkg<<endl;
    }
    fSig->Close();





  return;
}
				      


void compute2()
{
//   std::string fSigName = "/data_CMS/cms/ivo/AntiEMVA/Trees/MVAOutput/MvaOutput_AntiEMVA-DYJetsToLL_v4_Signal.root";
//   std::string fSigName = "/data_CMS/cms/ivo/AntiEMVA/Trees/MVAOutput/V4/MvaOutput_AntiEMVA_v4_Signal.root";
  std::string fSigName = "/data_CMS/cms/ivo/AntiEMVA/Trees/MVAOutput/MvaOutput_AntiEMVA-v5_Signal.root";
  TFile* fSig = new TFile (fSigName.data(),"READ") ;
  TTree* tSig = (TTree*)fSig->Get("tree");
  float nSig = tSig->GetEntries();
  cout<<"Signal entries :"<<nSig<<endl;
//   std::string fBkgName = "/data_CMS/cms/ivo/AntiEMVA/Trees/MVAOutput/MvaOutput_AntiEMVA-DYJetsToLL_v4_Backgrd.root";
//   std::string fBkgName = "/data_CMS/cms/ivo/AntiEMVA/Trees/MVAOutput/V4/MvaOutput_AntiEMVA_v4_Backgrd.root";
  std::string fBkgName = "/data_CMS/cms/ivo/AntiEMVA/Trees/MVAOutput/MvaOutput_AntiEMVA-v5_Backgrd.root";
  TFile* fBkg = new TFile (fBkgName.data(),"READ") ;
  TTree* tBkg = (TTree*)fBkg->Get("tree");
  float nBkg = tBkg->GetEntries();
  cout<<"Bkg entries :"<<nBkg<<endl;


  //61 % eff WP from V4_MoreStats
  float NoEleMatch_Barrel_Cut = 0.0306715 ;  
  float woG_Barrel_Cut = 0.992195 ; 
  float wGwoGSF_Barrel_Cut = 0.308324 ; 
  float wGwGSFwoPFMVA_Barrel_Cut = -0.0370998 ; 
  float wGwGSFwPFMVA_Barrel_Cut = 0.864643 ; 
  float NoEleMatch_Endcap_Cut = 0.0832094 ; 
  float woG_Endcap_Cut = 0.791665 ; 
  float wGwoGSF_Endcap_Cut = 0.675537 ; 
  float wGwGSFwoPFMVA_Endcap_Cut = 0.87047 ; 
  float wGwGSFwPFMVA_Endcap_Cut = 0.233711 ; 

  
  TCut Cut =Form("NoEleMatch_Barrel>%f && woG_Barrel>%f && wGwoGSF_Barrel>%f && wGwGSFwoPFMVA_Barrel>%f &&wGwGSFwPFMVA_Barrel>%f && NoEleMatch_Endcap>%f && woG_Endcap>%f && wGwoGSF_Endcap>%f && wGwGSFwoPFMVA_Endcap>%f &&wGwGSFwPFMVA_Endcap>%f",
		 NoEleMatch_Barrel_Cut,
		 woG_Barrel_Cut, 
		 wGwoGSF_Barrel_Cut, 
		 wGwGSFwoPFMVA_Barrel_Cut,
		 wGwGSFwPFMVA_Barrel_Cut,
		 NoEleMatch_Endcap_Cut, 
		 woG_Endcap_Cut, 
		 wGwoGSF_Endcap_Cut, 
		 wGwGSFwoPFMVA_Endcap_Cut,
		 wGwGSFwPFMVA_Endcap_Cut
		 );
  
  TH1F* hSigPass = new TH1F("hSig","",1,-10,10);
  tSig->Draw("Tau_GenEleMatch>>hSig",Cut);
  float nSigPass = hSigPass->GetEntries();
  
  TH1F* hBkgPass = new TH1F("hBkg","",1,-10,10);
  tBkg->Draw("Tau_GenEleMatch>>hBkg",Cut);
  float nBkgPass = hBkgPass->GetEntries();
  
  if (DEBUG){
    cout<<"Signal entries pass :"<<nSigPass<<" Efficiency -> "<<nSigPass/nSig<<endl;
    cout<<"Bkg entries pass :"<<nBkgPass<<" FakeRate -> "<<nBkgPass/nBkg<<endl;
  }
  fSig->Close();
  fBkg->Close();
  
  return;
}

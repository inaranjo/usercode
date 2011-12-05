#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TFrame.h"
#include "TStyle.h"
#include "TFile.h"
#include "iostream"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

TCanvas * getDefaultCanvas(float x=10,float y=30,float w=650,float h=600)
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1) ;
  TCanvas *c1 = new TCanvas("c2","canvas",x,y,w,h);
  c1->SetGridx() ;
  c1->SetGridy() ;
  c1->SetFillStyle(4000);
  c1->SetLeftMargin(0.16);

  return c1;
}

void MassPlot() 
{

  TCanvas *c1 = getDefaultCanvas();
  c1->cd() ;

  int xmax = 200;
  int bin = 60;
  THStack *hMass_MC = new THStack("hMass_MC","Invariant mass Mu Tau");
  TH1F *hMass_data = new TH1F("hMass_data","Invariant mass Mu Tau",bin,0,xmax);
  hMass_data->SetMarkerStyle(21);
  hMass_data->SetXTitle("Mass[GeV]");
  hMass_data->SetYTitle("Events");
  hMass_data->SetStats(0);
  hMass_data->GetYaxis()->SetTitleOffset(1.4);

  TH1F *hMassDY = new TH1F ("hMassDY", "hMassDY",bin,0,xmax);
  hMassDY->SetFillColor(kRed);
  hMassDY->SetXTitle("Mass[GeV]");
  TH1F *hMassWLnu = new TH1F ("hMassWLnu", "hMassWLnu",bin,0,xmax);
  hMassWLnu->SetFillColor(kBlue);
  hMassWLnu->SetXTitle("Mass[GeV]");
  TH1F *hMassQCD = new TH1F ("hMassQCD", "hMassQCD",bin,0,xmax);
  hMassQCD->SetFillColor(kGreen);
  hMassQCD->SetXTitle("Mass[GeV]");
  TH1F *hMassTTJets = new TH1F ("hMassTTJets", "hMassTTJets",bin,0,xmax);
  hMassTTJets->SetFillColor(kYellow);
  hMassTTJets->SetXTitle("Mass[GeV]");
  
  
  TFile *file0 ;
  Float_t mass;
  file0 = TFile::Open("output_data.root");
  file0->cd();
  treeTnP->SetBranchAddress("mass",&mass);
  TTree* treeTnP = (TTree*) file0->Get("treeTnP") ;
  Int_t nentries0 = (Int_t)treeTnP->GetEntries();
  for (Int_t i=0;i<nentries0;i++) {
    treeTnP->GetEntry(i);
    hMass_data->Fill(mass);
  }
  
  hMass_data->Draw("e1p");

    
  TFile *file ;
  Float_t mass_MC,weight;
  Int_t Index;
  file = TFile::Open("output_MC.root");
  file->cd();
  TTree* treeTnP_MC = (TTree*) file->Get("treeTnP") ;
  treeTnP_MC->SetBranchAddress("mass",&mass_MC);
  treeTnP_MC->SetBranchAddress("Index",&Index);
  treeTnP_MC->SetBranchAddress("weight",&weight);
  TTree* treeTnP_MC = (TTree*) file->Get("treeTnP") ;
  Int_t nentries = (Int_t)treeTnP_MC->GetEntries();
  for (Int_t i=0;i<nentries;i++) {
    if (Index == 0){
      treeTnP_MC->GetEntry(i);
      hMassDY->Fill(mass_MC,weight);
      
    }
    if (Index == 1){
      treeTnP_MC->GetEntry(i);
      hMassWLnu->Fill(mass_MC,weight);
    }
    if (Index == 2){
      treeTnP_MC->GetEntry(i);
      hMassQCD->Fill(mass_MC,weight);
    }
    if (Index == 3){
      treeTnP_MC->GetEntry(i);
      hMassTTJets->Fill(mass_MC,weight);
    }
  }
  hMass_MC->Add(hMassTTJets);
  hMass_MC->Add(hMassQCD);
  hMass_MC->Add(hMassWLnu);
  hMass_MC->Add(hMassDY);
 
  
  
  hMass_MC->Draw("same");
  hMass_data->Draw("e1p same");
  TLegend *legend = new TLegend(0.7,0.65,0.86,0.88);
  legend->AddEntry(hMass_data,"data","p");
  legend->AddEntry(hMassDY,"DY","f");
  legend->AddEntry(hMassWLnu,"WLnu","f");
  legend->AddEntry(hMassQCD,"QCD","f");
  legend->AddEntry(hMassTTJets,"TTJets","f");
  legend->Draw();
  
}

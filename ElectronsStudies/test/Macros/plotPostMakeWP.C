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

void plotTauPtRatio()
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

  TFile *fSignal = new TFile("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_makeWP_DYJetsToLL.root","READ"); 

  if(fSignal->IsZombie()){
    std::cout<<"No such file!"<<std::endl;
    return;
  }

  TTree*inTree = (TTree*)fSignal->Get("tree");

  float Tau_Pt;
  float Tau_LeadHadronPt;
  float Tau_AntiEMVA2_WP75; 

  inTree->SetBranchAddress("Tau_AntiEMVA2_WP75", &Tau_AntiEMVA2_WP75 );
  inTree->SetBranchAddress("Tau_Pt", &Tau_Pt );
  inTree->SetBranchAddress("Tau_LeadHadronPt", &Tau_LeadHadronPt );

  inTree->SetBranchStatus("Tau_AntiEMVA2_WP75", 1);
  inTree->SetBranchStatus("Tau_Pt", 1);
  inTree->SetBranchStatus("Tau_LeadHadronPt", 1);

  TH1F* hLeadHadronPtOverTauPt = new TH1F( "hLeadHadronPtOverTauPt","" ,100,-0.2,1.2);
//   hLeadHadronPtOverTauPt->Sumw2();
//   hLeadHadronPtOverTauPt->SetMarkerStyle(21);
  hLeadHadronPtOverTauPt->SetLineWidth(2);
  hLeadHadronPtOverTauPt->GetXaxis()->SetTitle("LeadHadron P_{T} / #tau P_{T} at WP75");
  int nEntries = inTree->GetEntries();

  for (int i=0;i<nEntries;i++){
    if(i%10000==0) cout<< i<<endl;

    inTree->GetEntry(i);

    if(Tau_AntiEMVA2_WP75<0.5)continue;
//     cout<<Tau_LeadHadronPt/Tau_Pt<<endl;
    hLeadHadronPtOverTauPt->Fill(Tau_LeadHadronPt/Tau_Pt);
  }


  cout<<hLeadHadronPtOverTauPt->GetEntries()<<endl;

  hLeadHadronPtOverTauPt->Draw();

  string outputName = "plots/plotPostMakeWP";
  c1->Print(std::string(outputName).append(".png").data());
  c1->Print(std::string(outputName).append(".pdf").data());
  
  fSignal->Close();

}







void plot()
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

  TFile *fSignal = new TFile("/data_CMS/cms/ivo/AntiEMVA/Trees/Trees_ForV4/AntiEMVA_AntiEMVATrees-DYJetsToLL-madgraph-PUS6.root","READ"); 

  if(fSignal->IsZombie()){
    std::cout<<"No such file!"<<std::endl;
    return;
  }

  TTree*inTree = (TTree*)fSignal->Get("AntiEMVAAnalyzer2/tree");

  TH1F* h = new TH1F( "h","" ,10,-2,5);

  cout<<inTree->GetEntries()<<endl;
  inTree->Draw("Elec_GenEleMatch>>h","Tau_Pt>20");
  cout<< h->GetEntries()<<endl;

  h->Draw();

   string outputName = "plots/test";
  c1->Print(std::string(outputName).append(".png").data());
  c1->Print(std::string(outputName).append(".pdf").data());
  fSignal->Close();

}

//--------------------------------------------------------------------------------------------------
// makeRootFilesAntiEMVA 
//
// Macro creating signal and background rootFiles from the analysis nTuples. 
// Signal is VBFH, GGFH and VH higgs MC.
// Background is DY Zee MC.
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

#define DEBUG false

void makeRoot(string type = "VBFH125"
	      )
{
  std::string inputFileName ="";
  if(type == "Zee")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/SplitDY/nTupleDYJ_EToTau_ElecTau_nominal.root";
  else if(type == "Zej")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/SplitDY/nTupleDYJ_JetToTau_ElecTau_nominal.root";
  else if(type == "VBFH110")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleVBFH110_ElecTau_nominal.root";
  else if(type == "GGFH110")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleGGFH110_ElecTau_nominal.root";
  else if(type == "VH110")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleVH110_ElecTau_nominal.root";

  else if(type == "VBFH115")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleVBFH110_ElecTau_nominal.root";
  else if(type == "GGFH115")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleGGFH110_ElecTau_nominal.root";
  else if(type == "VH115")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleVH110_ElecTau_nominal.root";
  else if(type == "VBFH120")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleVBFH110_ElecTau_nominal.root";
  else if(type == "GGFH120")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleGGFH110_ElecTau_nominal.root";
  else if(type == "VH120")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleVH110_ElecTau_nominal.root";
  else if(type == "VBFH125")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleVBFH125_ElecTau_nominal.root";
  else if(type == "GGFH125")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleGGFH125_ElecTau_nominal.root";
  else if(type == "VH125")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleVH125_ElecTau_nominal.root";
  else if(type == "VBFH110")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleVBFH110_ElecTau_nominal.root";
  else if(type == "GGFH110")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleGGFH110_ElecTau_nominal.root";
  else if(type == "VH110")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleVH110_ElecTau_nominal.root";
  else if(type == "VBFH110")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleVBFH110_ElecTau_nominal.root";
  else if(type == "GGFH110")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleGGFH110_ElecTau_nominal.root";
  else if(type == "VH110")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleVH110_ElecTau_nominal.root";
  else if(type == "VBFH110")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleVBFH110_ElecTau_nominal.root";
  else if(type == "GGFH110")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleGGFH110_ElecTau_nominal.root";
  else if(type == "VH110")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleVH110_ElecTau_nominal.root";
  else if(type == "VBFH110")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleVBFH110_ElecTau_nominal.root";
  else if(type == "GGFH110")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleGGFH110_ElecTau_nominal.root";
  else if(type == "VH110")
    inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTupleVH110_ElecTau_nominal.root";
  
  TFile* inputFile = new TFile (inputFileName.data(),"READ");
  if(inputFile->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
  TTree* inputTree = (TTree*)inputFile->Get("outTreePtOrd");

  ULong64_t run,event,lumi;
  float ptL1,ptL2,scEtaL1,etaL1,etaL2;
  float MEtMVA;
  float diTauVisMass;
  float diTauVisPt;
  float pt1,pt2;
  float visibleTauMass;
  float dPhiL1L2;
  float dPhiL1J1;
  float dPhiL1J2;
  float AntiEMVA3raw;

  inputTree->SetBranchAddress("run", &run);
  inputTree->SetBranchAddress("event", &event);
  inputTree->SetBranchAddress("lumi", &lumi);

  inputTree->SetBranchAddress("ptL1", &ptL1);
  inputTree->SetBranchAddress("etaL1", &etaL1);
  inputTree->SetBranchAddress("scEtaL1", &scEtaL1);
  inputTree->SetBranchAddress("ptL2", &ptL2);
  inputTree->SetBranchAddress("etaL2", &etaL2);
  inputTree->SetBranchAddress("MEtMVA", &MEtMVA);
  inputTree->SetBranchAddress("diTauVisMass", &diTauVisMass);
  inputTree->SetBranchAddress("diTauVisPt", &diTauVisPt);
  inputTree->SetBranchAddress("pt1", &pt1);
  inputTree->SetBranchAddress("pt2", &pt2);
  inputTree->SetBranchAddress("visibleTauMass", &visibleTauMass);
  inputTree->SetBranchAddress("dPhiL1L2", &dPhiL1L2);
  inputTree->SetBranchAddress("dPhiL1J1", &dPhiL1J1);
  inputTree->SetBranchAddress("dPhiL1J2", &dPhiL1J2);
  inputTree->SetBranchAddress("AntiEMVA3raw", &AntiEMVA3raw);

  inputTree->SetBranchStatus("run", 1);
  inputTree->SetBranchStatus("event", 1);
  inputTree->SetBranchStatus("lumi", 1);
  inputTree->SetBranchStatus("ptL1", 1);
  inputTree->SetBranchStatus("ptL2", 1);
  inputTree->SetBranchStatus("etaL1", 1);
  inputTree->SetBranchStatus("scEtaL1", 1);
  inputTree->SetBranchStatus("etaL2", 1);
  inputTree->SetBranchStatus("MEtMVA", 1);
  inputTree->SetBranchStatus("diTauVisMass", 1);
  inputTree->SetBranchStatus("pt1", 1);
  inputTree->SetBranchStatus("pt2", 1);
  inputTree->SetBranchStatus("visibleTauMass", 1);
  inputTree->SetBranchStatus("dPhiL1L2", 1);
  inputTree->SetBranchStatus("dPhiL1J1", 1);
  inputTree->SetBranchStatus("dPhiL1J2", 1);
  inputTree->SetBranchStatus("AntiEMVA3raw", 1);

  std::string outputFileName = Form("/data_CMS/cms/ivo/AntiZee/root/tree_AntiZeeMVA_%s.root",type.data());

  TFile* outputFile = new TFile (outputFileName.data(),"RECREATE");
  TTree* mytree = new TTree("tree", "tree");

  ULong64_t t_run,t_event,t_lumi;
  float t_ptL1,t_ptL2,t_etaL1,t_etaL2,t_scEtaL1;
  float t_MEtMVA;
  float t_diTauVisMass;
  float t_diTauVisPtOverPtSum;
  float t_pt1,t_pt2;
  float t_visibleTauMass;
  float t_dPhiL1L2;
  float t_dPhiL1J1;
  float t_dPhiL1J2;
  float t_AntiEMVA3raw;

  //counters
  mytree->Branch("run",&t_run,"run/l");
  mytree->Branch("event",&t_event,"event/l");
  mytree->Branch("lumi",&t_lumi,"lumi/l");

  mytree->Branch("ptL1",&t_ptL1,"ptL1/F");
  mytree->Branch("etaL1",&t_etaL1,"etaL1/F");
  mytree->Branch("scEtaL1",&t_scEtaL1,"scEtaL1/F");
  mytree->Branch("ptL2",&t_ptL2,"ptL2/F");
  mytree->Branch("etaL2",&t_etaL2,"etaL2/F");
  mytree->Branch("diTauVisPtOverPtSum",&t_diTauVisPtOverPtSum,"diTauVisPtOverPtSum/F");
  mytree->Branch("MEtMVA", &t_MEtMVA,"MEtMVA/F");
  mytree->Branch("diTauVisMass", &t_diTauVisMass,"diTauVisMass/F");
  mytree->Branch("pt1", &t_pt1,"pt1/F");
  mytree->Branch("pt2", &t_pt2,"pt2/F");
  mytree->Branch("visibleTauMass", &t_visibleTauMass,"visibleTauMass/F");
  mytree->Branch("dPhiL1L2", &t_dPhiL1L2,"dPhiL1L2/F");
  mytree->Branch("dPhiL1J1", &t_dPhiL1J1,"dPhiL1J1/F");
  mytree->Branch("dPhiL1J2", &t_dPhiL1J2,"dPhiL1J2/F");
  mytree->Branch("AntiEMVA3raw", &t_AntiEMVA3raw,"AntiEMVA3raw/F");

  int nEntries = inputTree->GetEntries();

  cout<< "Number of entries : "<<nEntries<<endl;

  for (int iEntry = 0; iEntry<nEntries ; iEntry++){
    if(iEntry%10000==0) cout << iEntry << endl;

    inputTree->GetEntry(iEntry);

    t_run = run;
    t_event = event;
    t_lumi = lumi;
    t_ptL1 = ptL1;
    t_etaL1 = etaL1;
    t_scEtaL1 = scEtaL1;
    t_ptL2 = ptL2;
    t_etaL2 = etaL2;
    t_MEtMVA = MEtMVA;
    t_diTauVisMass = diTauVisMass;
    t_diTauVisPtOverPtSum = diTauVisPt/(ptL1+ptL2);
    t_pt1 = pt1;
    t_pt2 = pt2;
    t_visibleTauMass = visibleTauMass;
    t_dPhiL1L2 = dPhiL1L2;
    t_dPhiL1J1 = dPhiL1J1;
    t_dPhiL1J2 = dPhiL1J2;
    t_AntiEMVA3raw = AntiEMVA3raw;
    t_dPhiL1L2 = dPhiL1L2;

    if(DEBUG){
      cout<<endl;
      cout<<" run : "<<t_run<<endl;
      cout<<" event : "<<t_event<<endl;
      cout<<" lumi : "<<t_lumi<<endl;
    }
    mytree->Fill();
  }
  mytree->Write();
  
  cout<<"Creating file : "<<outputFileName.data()<<endl;
  inputFile->Close();
  outputFile->Close();
  return;
}



void makeAll(){

  makeRoot("VBFH125");
//   makeRoot("GGFH125");
//   makeRoot("VH125");
//   makeRoot("Zee");
//   makeRoot("Zej");

}

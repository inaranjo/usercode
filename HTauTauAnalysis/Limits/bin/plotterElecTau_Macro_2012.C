#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>

#include "TChain.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TH1F.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCut.h"
#include "TArrayF.h"
#include "TStyle.h"

#define DEBUG false


void plot( string variable_ = "diTauVisMass",
	   string xTitle_   = " ; mass (GeV) ; Events",
	   float binWidth_  = -1,
	   float nbins_     = -1,
	   string MtCut_    = "MtLeg1<40",
	   string charge_   = "diTauCharge==0",
	   string outFile_  = "visibleMass"
	   ){
    
  int debug = 0;

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetStatColor(0);
  gStyle->SetOptStat(0000000);
  gStyle->SetOptFit(0111);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetPalette(1);

  TPad* pad1 = new TPad("pad1","",0.05,0.27,0.96,0.97);
  TPad* pad2 = new TPad("pad2","",0.05,0.02,0.96,0.26);
  pad1->SetFillColor(0);
  pad2->SetFillColor(0);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();

//   TFile fData("/data_CMS/cms/ivo/HTauTauAnalysis/Trees/ElecTauStream_06Jul2012/OpenNtuples/nTupleRun2012-ElecTau-All_run_Open_ElecTauStream_Nominal.root");
//   TFile fData("/data_CMS/cms/ivo/HTauTauAnalysis/Trees/ElecTauStream_25Jul2012/OpenNtuples/nTupleRun2012-ElecTau-All_run_Open_ElecTauStream_Nominal.root");
  TFile fData("/data_CMS/cms/ivo/HTauTauAnalysis/Trees/ElecTauStream_01Aug2012/OpenNtuples/nTupleRun2012-ElecTau-All_run_Nominal.root");
  TFile fDYJets("/data_CMS/cms/ivo/HTauTauAnalysis/Trees/ElecTauStream_01Aug2012/OpenNtuples/nTupleDYJetsToLL-ElecTau_run_Nominal.root");
  TFile fWJets("/data_CMS/cms/ivo/HTauTauAnalysis/Trees/ElecTauStream_01Aug2012/OpenNtuples/nTupleWJetsToLNu-ElecTau_run_Nominal.root");
  TFile fTTJets("/data_CMS/cms/ivo/HTauTauAnalysis/Trees/ElecTauStream_01Aug2012/OpenNtuples/nTupleTTJets-ElecTau_run_Nominal.root");

//   TFile fDYJets("/data_CMS/cms/ivo/HTauTauAnalysis/Trees/ElecTauStream_30Jul2012/OpenNtuples/nTupleDYJetsToLL-ElecTau-madgraph-tarball-iter1_run_Open_ElecTauStream_Nominal.root");
//   TFile fWJets("/data_CMS/cms/ivo/HTauTauAnalysis/Trees/ElecTauStream_30Jul2012/OpenNtuples/nTupleWJetsToLNu-ElecTau-madgraph-tarball-iter1_run_Open_ElecTauStream_Nominal.root");
//   TFile fTTJets("/data_CMS/cms/ivo/HTauTauAnalysis/Trees/ElecTauStream_30Jul2012/OpenNtuples/nTupleTTJets-ElecTau-powheg-iter1_run_Open_ElecTauStream_Nominal.root");


//   TFile fDYJets("/data_CMS/cms/ivo/HTauTauAnalysis/Trees/ElecTauStream_06Jul2012/OpenNtuples/nTupleDYJetsToLL-ElecTau-madgraph-tarball_run_Open_ElecTauStream_Nominal.root");
//   TFile fWJets("/data_CMS/cms/ivo/HTauTauAnalysis/Trees/ElecTauStream_06Jul2012/OpenNtuples/nTupleWJetsToLNu-ElecTau-madgraph-tarball_run_Open_ElecTauStream_Nominal.root");
//   TFile fTTJets("/data_CMS/cms/ivo/HTauTauAnalysis/Trees/ElecTauStream_06Jul2012/OpenNtuples/nTupleTTJets-ElecTau-powheg_run_Open_ElecTauStream_Nominal.root");

  TTree* tData   = (TTree*)fData.Get("outTreePtOrd");
  TTree* tDYJets = (TTree*)fDYJets.Get("outTreePtOrd");
  TTree* tWJets  = (TTree*)fWJets.Get("outTreePtOrd");
  TTree* tTTJets = (TTree*)fTTJets.Get("outTreePtOrd");

  TArrayF bins(nbins_+1);
  if(binWidth_<0){
  int nbins_ = 13;
  bins[0] = 0;    bins[1] = 30;    bins[2] = 40;    bins[3] = 50;    bins[4]  = 60;    bins[5] = 70;
  bins[6] = 80;   bins[7] = 90;    bins[8] = 100;   bins[9] = 110;   bins[10] = 120;  bins[11] = 140; 
  bins[12] = 160; bins[13] = 200;
  }
  if(binWidth_>0) 
    for(int k = 0 ; k <= nbins_ ; k++) bins[k] = binWidth_*k; 
  TString variable(variable_.c_str());
  TString labels(xTitle_.c_str());

  TH1F* hData    = new TH1F("hData"   ,labels, nbins_, bins.GetArray());
  TH1F* hQCD     = new TH1F("hQCD"    ,labels, nbins_, bins.GetArray());
  TH1F* hDYJets  = new TH1F("hDYJets" ,labels, nbins_, bins.GetArray());
  TH1F* hZFakes  = new TH1F("hZFakes" ,labels, nbins_, bins.GetArray());
  TH1F* hWJets   = new TH1F("hWJets"  ,labels, nbins_, bins.GetArray());
  TH1F* hWJetsSS = new TH1F("hWJetsSS",labels, nbins_, bins.GetArray());
  TH1F* hTTJets  = new TH1F("hTTJets" ,labels, nbins_, bins.GetArray());
  TH1F* hAll     = new TH1F("hAll"    ,labels, nbins_, bins.GetArray());

  hData->SetMarkerStyle(kFullCircle);
  hQCD->SetFillColor(kMagenta-9);
  hDYJets->SetFillColor(kYellow-9);
  hWJets->SetFillColor(kRed-3);
  hZFakes->SetFillColor(kBlue-6);
  hTTJets->SetFillColor(kBlue-2);

  //tightestMVAWP-> EleIDMVA
  //tightestHPSDBWP-> TauIso hps db corr
  //tightestHPSMVAWP-> TauIso hps mva 
  // full signal selection leg 1 no isolated
  TCut sbinRel(      "ptL1>24 && ptL2>20 && tightestHPSMVAWP>0 && tightestMVAWP>=1 && diTauCharge==0 && elecFlag==0 && HLTx && HLTmatch");
  // full signal selection
  TCut sbin(         "ptL1>24 && ptL2>20 && tightestHPSMVAWP>0 && tightestMVAWP>=1                   && elecFlag==0 && HLTx && HLTmatch && combRelIsoLeg1DBeta<0.10");
  // selection of opposite-sign W-enriched control region
  TCut sbinAntiW(    "ptL1>24 && ptL2>20 && tightestHPSMVAWP>0 && tightestMVAWP>=1 && diTauCharge==0 && elecFlag==0 && HLTx && HLTmatch && combRelIsoLeg1DBeta<0.10");
  // selection of same-sign enriched control-region 
  TCut sbinSSAntiW(  "ptL1>24 && ptL2>20 && tightestHPSMVAWP>0 && tightestMVAWP>=1 && diTauCharge!=0 && elecFlag==0 && HLTx && HLTmatch && combRelIsoLeg1DBeta<0.10");
  // selection of same-sign control region
  TCut sbinSS(       "ptL1>24 && ptL2>20 && tightestHPSMVAWP>0 && tightestMVAWP>=1 && diTauCharge!=0 && elecFlag==0 && HLTx && HLTmatch && combRelIsoLeg1DBeta<0.10");
  // selection of same sign loosely-isolated control region
  TCut sbinSSRelIso( "ptL1>24 && ptL2>20 && tightestHPSMVAWP>0 && tightestMVAWP>=1 && diTauCharge!=0 && elecFlag==0 && HLTx && HLTmatch && combRelIsoLeg1DBeta<0.30");


//   TCut sbinRel(      "ptL1>24 && ptL2>20 && tightestHPSMVAWP>0 && tightestMVAWP>=1 && diTauCharge==0 && HLTx && HLTmatch");
//   // full signal selection
//   TCut sbin(         "ptL1>24 && ptL2>20 && tightestHPSMVAWP>0 && tightestMVAWP>=1                   && HLTx && HLTmatch && combRelIsoLeg1DBeta<0.10");
//   // selection of opposite-sign W-enriched control region
//   TCut sbinAntiW(    "ptL1>24 && ptL2>20 && tightestHPSMVAWP>0 && tightestMVAWP>=1 && diTauCharge==0 && HLTx && HLTmatch && combRelIsoLeg1DBeta<0.10");
//   // selection of same-sign enriched control-region 
//   TCut sbinSSAntiW(  "ptL1>24 && ptL2>20 && tightestHPSMVAWP>0 && tightestMVAWP>=1 && diTauCharge!=0 && HLTx && HLTmatch && combRelIsoLeg1DBeta<0.10");
//   // selection of same-sign control region
//   TCut sbinSS(       "ptL1>24 && ptL2>20 && tightestHPSMVAWP>0 && tightestMVAWP>=1 && diTauCharge!=0 && HLTx && HLTmatch && combRelIsoLeg1DBeta<0.10");
//   // selection of same sign loosely-isolated control region
//   TCut sbinSSRelIso( "ptL1>24 && ptL2>20 && tightestHPSMVAWP>0 && tightestMVAWP>=1 && diTauCharge!=0 && HLTx && HLTmatch && combRelIsoLeg1DBeta<0.10");

//   // full signal selection leg 1 no isolated
//   TCut sbinRel(      "ptL1>20 && ptL2>20 && tightestHPSDBWP>0 && tightestMVAWP>=1 && diTauCharge==0 && elecFlag==0");
//   // full signal selection
//   TCut sbin(         "ptL1>20 && ptL2>20 && tightestHPSDBWP>0 && tightestMVAWP>=1                   && elecFlag==0 && combRelIsoLeg1DBeta<0.10");
//   // selection of opposite-sign W-enriched control region
//   TCut sbinAntiW(    "ptL1>20 && ptL2>20 && tightestHPSDBWP>0 && tightestMVAWP>=1 && diTauCharge==0 && elecFlag==0 && combRelIsoLeg1DBeta<0.10");
//   // selection of same-sign enriched control-region 
//   TCut sbinSSAntiW(  "ptL1>20 && ptL2>20 && tightestHPSDBWP>0 && tightestMVAWP>=1 && diTauCharge!=0 && elecFlag==0 && combRelIsoLeg1DBeta<0.10");
//   // selection of same-sign control region
//   TCut sbinSS(       "ptL1>20 && ptL2>20 && tightestHPSDBWP>0 && tightestMVAWP>=1 && diTauCharge!=0 && elecFlag==0 && combRelIsoLeg1DBeta<0.10");
//   // selection of same sign loosely-isolated control region
//   TCut sbinSSRelIso( "ptL1>20 && ptL2>20 && tightestHPSDBWP>0 && tightestMVAWP>=1 && diTauCharge!=0 && elecFlag==0 && combRelIsoLeg1DBeta<0.10");


  //Mt cut
  TCut highMt("MtLeg1>60");
  TCut lowMt("MtLeg1<40");
  TCut slowMt(MtCut_.c_str());
  
  //Pair Charge cut
  TCut SSCharge("diTauCharge!=0");
  TCut OSCharge("diTauCharge==0");
  TCut sCharge(charge_.c_str());

  //Tau antiElectron discr
  TCut tAntiEMVA("tightestAntiECutWP>1"); 
//   TCut tAntiEMVA("tightestAntiEMVAWP>=6"); 
//   TCut tAntiEMVA(""); 

  TCut tAntiEMVA2Tight1("tightestAntiEMVAWP>=5"); 
  TCut tAntiEMVA2Tight2("tightestAntiEMVAWP>=6"); 
  TCut tAntiEMVA2VTight1("tightestAntiEMVAWP>=7"); 
  TCut tAntiEMVA2VTight2("tightestAntiEMVAWP>=8"); 

  /////////// get approximated pu reweighting factors /////////////////
  TH1F* numPVData     = new TH1F("numPVData","",  20,0,40);
  TH1F* numPVDYJets   = new TH1F("numPVDYJets","",20,0,40);
  TH1F* numPVWJets    = new TH1F("numPVWJets","", 20,0,40);

  tData->Draw("numPV>>numPVData"    ,sbinRel);
  tDYJets->Draw("numPV>>numPVDYJets",sbinRel);
  tWJets->Draw("numPV>>numPVWJets"  ,sbinRel);

  numPVData->Scale(  1./numPVData->Integral());
  numPVDYJets->Scale(1./numPVDYJets->Integral());
  numPVWJets->Scale( 1./numPVWJets->Integral());

  float totalWeight = 0.;
  string scaleFactDYJets = "( ";
  for(int i = 1; i < numPVData->GetNbinsX(); i++){
    float binWidth = numPVDYJets->GetBinWidth(i);
    float weightBin_i = numPVDYJets->GetBinContent(i)>0 ?
      numPVData->GetBinContent(i)/numPVDYJets->GetBinContent(i) : 0.0;
    if(DEBUG){
      cout<<"binWidth : "<<binWidth<<endl;
      cout<<"numPVData->GetBinContent(i): "<<numPVData->GetBinContent(i)<<endl;
      cout<<"numPVDYJets->GetBinContent(i): "<<numPVDYJets->GetBinContent(i)<<endl;
      cout<<"weightBin_i : "<<weightBin_i<<endl;
    }

    totalWeight += weightBin_i;

    scaleFactDYJets += string( Form("(numPV>=%f && numPV<%f)*%f",binWidth*(i-1), binWidth*i, weightBin_i) );
    if(i < (numPVData->GetNbinsX() - 1)) scaleFactDYJets += " + ";
  }
  scaleFactDYJets += " )";

  scaleFactDYJets="1";//puWeight2012
//   scaleFactDYJets="puWeight2012";

  cout << scaleFactDYJets << endl;
  cout << endl;

  string scaleFactWJets = "( ";
  for(int i = 1; i < numPVData->GetNbinsX(); i++){
    float binWidth = numPVWJets->GetBinWidth(i);
    float weightBin_i = numPVWJets->GetBinContent(i)>0 ?
      numPVData->GetBinContent(i)/numPVWJets->GetBinContent(i) : 0.0;
 
    totalWeight += weightBin_i;

    scaleFactWJets += string( Form("(numPV>=%f && numPV<%f)*%f",binWidth*(i-1), binWidth*i, weightBin_i) );
    if(i < (numPVData->GetNbinsX() - 1)) scaleFactWJets += " + ";
  }
  scaleFactWJets += " )";

  scaleFactWJets="1";//puWeight2012
//   scaleFactWJets="puWeight2012";

  cout << scaleFactWJets << endl;

  
  if(debug){
    cout << totalWeight/numPVData->GetNbinsX() << endl;

    TH1F* numPVData     = new TH1F("numPVData","",  20,0,120);
    TH1F* numPVWJets    = new TH1F("numPVWJets","", 20,0,120);

    numPVData->Reset();
    numPVWJets->Reset();
    tWJets->Draw("MtLeg1>>numPVData",sbin);
    tWJets->Draw("MtLeg1>>numPVWJets",TString(scaleFactWJets.c_str())*sbin);
    
    cout << numPVWJets->Integral() << " -- " << numPVData->Integral() << endl;
    numPVWJets->SetLineColor(kRed);
    numPVWJets->Draw();
    numPVData->Draw("SAME");

    c1->SaveAs("plots/plotElecTauDataMC_ControlWJets.png");

    return;

  }


  // luminosity of data is 5082 pb
  float lumiFact = 5082/1000;

  // opposite-sign to same-sign ratio for QCD
  float OStoSSRatio = 1.11;
  if( (string(sCharge.GetTitle())).find("diTauCharge!=0")!=string::npos ){
    OStoSSRatio = 1.0;
    cout << "SS selection" << endl;
  }

  /////////////////////////////////////////////////////////////////////
  // estimation of W+jets
  TH1F* h1 = new TH1F("h1","",1,-10,10);
  tWJets->Draw("etaL1>>h1",   TString(("sampleWeight*puWeight2012*"+scaleFactWJets).c_str())*(sbin&&OSCharge&&highMt&&tAntiEMVA));//sampleWeight*puWeight2012
  float WsbinMCWEnriched  = h1->Integral()*lumiFact;
  h1->Reset();
  tWJets->Draw("etaL1>>h1",   TString(("sampleWeight*puWeight2012*"+scaleFactWJets).c_str())*(sbin&&OSCharge&&lowMt&&tAntiEMVA));
  float WsbinMC       = h1->Integral()*lumiFact;
  h1->Reset();
  tTTJets->Draw("etaL1>>h1",  TString(("sampleWeight*puWeight2012*"+scaleFactWJets).c_str())*(sbin&&OSCharge&&highMt&&tAntiEMVA));
  float TTsbinWEnriched   = h1->Integral()*lumiFact;
  h1->Reset();
  tDYJets->Draw("etaL1>>h1",  TString(("sampleWeight*puWeight2012*(  (abs(genDecay)==23*15) + (abs(genDecay)!=23*15)*( (TMath::Abs(etaL2)<1.5)*0.85 + (TMath::Abs(etaL2)>1.5)*0.65 ))*"+scaleFactDYJets).c_str())*(sbin&&OSCharge&&highMt&&tAntiEMVA));
  float ZttsbinWEnriched = h1->Integral()*lumiFact;
  h1->Reset();
  tData->Draw("etaL1>>h1",   sbin&&OSCharge&&highMt&&tAntiEMVA);
  float DatasbinWEnriched = h1->Integral();
  h1->Reset();

  float WscaleFactorOS =  1./WsbinMC;
  cout << "Wsbin (MC) = " << WsbinMC << endl;
  float Wsbin = (DatasbinWEnriched - TTsbinWEnriched - ZttsbinWEnriched)*(WsbinMC/WsbinMCWEnriched);
  cout << "Wsbin = (" << DatasbinWEnriched << " - " << TTsbinWEnriched << " - " << ZttsbinWEnriched << " )*" << WsbinMC/WsbinMCWEnriched << " = " << Wsbin << endl;
  WscaleFactorOS *= Wsbin;
  cout << " ==> scale factor (OS) = " << WscaleFactorOS << endl;

  /////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////
  // estimation of QCD
  /////////////////////////////////////////////////////////////////////

  //Wjets estimation in SS region
  h1->Reset();
  tWJets->Draw("etaL1>>h1",   TString(("sampleWeight*puWeight2012*"+scaleFactWJets).c_str())*(sbinSS&&highMt&&tAntiEMVA));
  float WsbinMCSSWEnriched  = h1->Integral()*lumiFact;
  h1->Reset();
  tWJets->Draw("etaL1>>h1",   TString(("sampleWeight*puWeight2012*"+scaleFactWJets).c_str())*(sbinSS&&lowMt&&tAntiEMVA));
  float WsbinMCSS       = h1->Integral()*lumiFact;
  h1->Reset();
  tTTJets->Draw("etaL1>>h1",  TString(("sampleWeight*puWeight2012*"+scaleFactWJets).c_str())*(sbinSS&&highMt&&tAntiEMVA));
  float TTsbinSSWEnriched = h1->Integral()*lumiFact;
  h1->Reset();
  tDYJets->Draw("etaL1>>h1",  TString(("sampleWeight*puWeight2012*(  (abs(genDecay)==23*15) + (abs(genDecay)!=23*15)*( (TMath::Abs(etaL2)<1.5)*0.85 + (TMath::Abs(etaL2)>1.5)*0.65 ))*"+scaleFactDYJets).c_str())*(sbinSS&&highMt&&tAntiEMVA));
  float ZttSSWEnriched = h1->Integral()*lumiFact;
  h1->Reset();
  tData->Draw("etaL1>>h1",   sbinSS&&highMt&&tAntiEMVA);
  float DatasbinSSWEnriched = h1->Integral();
  h1->Reset();
  tData->Draw("etaL1>>h1",   sbinSS&&lowMt&&tAntiEMVA);
  float DatasbinSS = h1->Integral();
  h1->Reset();
  
  float WscaleFactorSS =  1./WsbinMCSS;
  cout << "WsbinSS (MC) = " << WsbinMCSS << endl; 
  float WsbinSS   = (DatasbinSSWEnriched - TTsbinSSWEnriched - ZttSSWEnriched)*(WsbinMCSS/WsbinMCSSWEnriched);
  cout << "WsbinSS = (" << DatasbinSSWEnriched << " - " << TTsbinSSWEnriched << " - " << ZttSSWEnriched << " )*" << WsbinMCSS/WsbinMCSSWEnriched << " = " << WsbinSS << endl;
  WscaleFactorSS *= WsbinSS;
  cout << " ==> scale factor (SS) = " << WscaleFactorSS << endl;
  float QCDsbinSS = DatasbinSS - WsbinSS;
  float QCDsbin   = QCDsbinSS*OStoSSRatio;
  /////////////////////////////////////////////////////////////////////

  hDYJets->Sumw2();
  hWJets->Sumw2();
  hTTJets->Sumw2();
  hQCD->Sumw2();
  hAll->Sumw2();

  // Draw with cuts and weights !!!
  tData->Draw(  variable+">>hData",   "(index==0)"*(sbin&&sCharge&&slowMt&&tAntiEMVA));
  tData->Draw(  variable+">>hQCD",    "(index==0)"*(sbinSSRelIso&&slowMt&&tAntiEMVA));
  tDYJets->Draw(variable+">>hZFakes", TString(("sampleWeight*puWeight2012*(index==0)*(abs(genDecay)!=23*15)*( 1 +(TMath::Abs(etaL2)<1.5)*(0.85-1) + (TMath::Abs(etaL2)>1.5)*(0.65-1) )*"+scaleFactDYJets).c_str())*(sbin&&sCharge&&slowMt&&tAntiEMVA));
  tDYJets->Draw(variable+">>hDYJets", TString(("sampleWeight*puWeight2012*(index==0)*(abs(genDecay)==23*15)*"+scaleFactDYJets).c_str())*(sbin&&sCharge&&slowMt&&tAntiEMVA));
  tWJets->Draw( variable+">>hWJets",  TString(("sampleWeight*puWeight2012*"+scaleFactWJets).c_str())*(sbin&&sCharge&&slowMt&&tAntiEMVA));
  tWJets->Draw( variable+">>hWJetsSS",TString(("sampleWeight*puWeight2012*"+scaleFactWJets).c_str())*(sbinSS&&slowMt&&tAntiEMVA));
  tTTJets->Draw(variable+">>hTTJets", TString(("sampleWeight*puWeight2012*"+scaleFactWJets).c_str())*(sbin&&sCharge&&slowMt&&tAntiEMVA));

//   tDYJets->Draw(variable+">>hZFakes", TString("puWeight2012*sampleWeight*(index==0)*(abs(genDecay)!=23*15)*( 1 +(TMath::Abs(etaL2)<1.5)*(0.85-1) + (TMath::Abs(etaL2)>1.5)*(0.65-1) )")*(sbin&&sCharge&&slowMt));
//   tDYJets->Draw(variable+">>hDYJets", TString("puWeight2012*sampleWeight*(index==0)*(abs(genDecay)==23*15)"*(sbin&&sCharge&&slowMt)));
//   tWJets->Draw( variable+">>hWJets",  TString("puWeight2012*sampleWeight"*(sbin&&sCharge&&slowMt)));
//   tWJets->Draw( variable+">>hWJetsSS",TString("puWeight2012*sampleWeight"*(sbinSS&&slowMt)));
//   tTTJets->Draw(variable+">>hTTJets", TString("puWeight2012*sampleWeight"*(sbin&&sCharge&&slowMt)));

  //tVBFH130->Draw(variable+">>hVBFH130", "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbin);
  //tGGFH130->Draw(variable+">>hGGFH130", "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight*HqTWeight"*sbin);
  
  // Scale histograms
  hDYJets->Scale( lumiFact );
  hZFakes->Scale( lumiFact*1.2 );//20% increase due to eTau fakerate

  hWJets->Scale(  lumiFact*( int((string(sCharge.GetTitle())).find("diTauCharge==0")!=string::npos)*WscaleFactorOS +
					      int((string(sCharge.GetTitle())).find("diTauCharge!=0")!=string::npos)*WscaleFactorSS)  );
  hWJetsSS->Scale(lumiFact*WscaleFactorSS );
  hTTJets->Scale( lumiFact );
  hQCD->Add(   hWJetsSS, -1);
  hQCD->Scale( QCDsbin/hQCD->Integral());
//   hWJets->Add(hZFakes);

  hAll->Add(hDYJets);
  hAll->Add(hZFakes);
  hAll->Add(hWJets);
  hAll->Add(hTTJets);
  hAll->Add(hQCD);
  hAll->SetLineWidth(1);
  hAll->SetLineColor(kRed);
  hAll->SetFillColor(kRed);
  hAll->SetFillStyle(3003);

  for(int i = 1; i < hAll->GetNbinsX(); i++){
    float err = 0.;
    err += TMath::Power(hDYJets->GetBinError(i),2);
    err += TMath::Power(hWJets->GetBinError(i),2);
    err += TMath::Power(hTTJets->GetBinError(i),2);
    err += TMath::Power(hQCD->GetBinError(i),2);
    hAll->SetBinError(i,TMath::Sqrt(err));
  }

  // Add all together
  THStack* aStack = new THStack("aStack","");
  aStack->Add(hTTJets);
  aStack->Add(hQCD);
  aStack->Add(hWJets);
  aStack->Add(hZFakes);
  aStack->Add(hDYJets);

  hData->SetAxisRange( 0, TMath::Max(hData->GetMaximum(), hAll->GetMaximum())*1.20, "Y");

  hData->Sumw2();
  hData->Draw("P");
  aStack->Draw("HISTSAME");
  hAll->Draw("E2SAME");
  hData->Draw("PSAME");

  cout << "DYJets = " << hDYJets->Integral() << endl;
  cout << "ZFakes = " << hZFakes->Integral() << endl;
  cout << "TT = " << hTTJets->Integral() << endl;
  cout << "QCD = " << hQCD->Integral() << endl;
  cout << "WJets = " << hWJets->Integral() << endl;
  cout << "Data = " << hData->Integral() << endl;

  // Legend
  TLegend* leg = new TLegend(0.60,0.62,0.78,0.87,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);

  leg->SetHeader("#splitline{CMS Preliminary 2012}{#sqrt{s}=8 TeV, L=5.08 fb^{-1}, #tau_{e}#tau_{h}}");
  leg->AddEntry(hData,"Observed","P");
  leg->AddEntry(hDYJets,"Z#rightarrow#tau#tau","F");
  leg->AddEntry(hZFakes,"Z#rightarrow ee","F");
  leg->AddEntry(hWJets,"W+jets","F");
  leg->AddEntry(hQCD,"QCD","F");
  leg->AddEntry(hTTJets,"t#bar{t}","F");
  //leg->AddEntry(hVBFH130,"100 X H#rightarrow#tau#tau, m_{H}=130 GeV","F");
  leg->Draw();

  TH1F* hStack = (TH1F*)aStack->GetHistogram();

  pad2->cd();

  TH1F* hRatio = new TH1F("hRatio", " ; ; #frac{(DATA)}{MC}",
			  hStack->GetNbinsX(), 
			  hStack->GetXaxis()->GetXmin(), hStack->GetXaxis()->GetXmax());
  hRatio->SetMarkerStyle(kFullCircle);
  hRatio->SetMarkerSize(1.0);
  hRatio->SetLabelSize(0.12,"X");
  hRatio->SetLabelSize(0.10,"Y");
  hRatio->SetTitleSize(0.12,"Y");
  hRatio->SetTitleOffset(0.36,"Y");
  hData->Add(hAll,-1);
  hRatio->Divide( hData ,hAll,1.0,1.0);
  hRatio->SetAxisRange(-1,1,"Y");
  hRatio->Draw("P");
  TF1* line = new TF1("line","0",hRatio->GetXaxis()->GetXmin(),hStack->GetXaxis()->GetXmax());
  line->SetLineStyle(3);
  line->SetLineWidth(1.5);
  line->SetLineColor(kBlack);
  line->Draw("SAME");

  if(charge_.find("diTauCharge==0")!=string::npos){
    c1->SaveAs(("plots/plotElecTauDataMC_"+outFile_+".png").c_str());
//     c1->SaveAs(("plots/plotElecTauDataMC_"+outFile_+".pdf").c_str());
  }
  else{
    c1->SaveAs(("plots/plotElecTauDataMC_"+outFile_+"_SS.png").c_str());
//     c1->SaveAs(("plots/plotElecTauDataMC_"+outFile_+"_SS.pdf").c_str());
  }

  TCanvas *c2 = new TCanvas("c2","",800,800);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);
  c2->Divide(2,3);

  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetStatColor(0);
  gStyle->SetOptStat(0000000);
  gStyle->SetOptFit(0111);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetPalette(1);

  c2->cd(1);
  TH1F* hRatio1 = new TH1F("hRatio1", " ; ; #frac{(DATA)}{Ztt}",
			  hStack->GetNbinsX(), 
			  hStack->GetXaxis()->GetXmin(), hStack->GetXaxis()->GetXmax());
  hRatio1->Divide( hData ,hDYJets,1.0,1.0);
  hRatio1->Draw();
  //hData->Draw();
  //hDYJets->Draw();
  c2->cd(2);
  TH1F* hRatio2 = new TH1F("hRatio2", " ; ; #frac{(DATA)}{Zee}",
			  hStack->GetNbinsX(), 
			  hStack->GetXaxis()->GetXmin(), hStack->GetXaxis()->GetXmax());
  hRatio2->Divide( hData ,hZFakes,1.0,1.0);
  hRatio2->Draw();
  c2->cd(3);
  TH1F* hRatio3 = new TH1F("hRatio3", " ; ; #frac{(DATA)}{WJets}",
			  hStack->GetNbinsX(), 
			  hStack->GetXaxis()->GetXmin(), hStack->GetXaxis()->GetXmax());  
  hRatio3->Divide( hData ,hWJets,1.0,1.0);
  hRatio3->Draw();
  c2->cd(4);
  TH1F* hRatio4 = new TH1F("hRatio4", " ; ; #frac{(DATA)}{QCD}",
			  hStack->GetNbinsX(), 
			  hStack->GetXaxis()->GetXmin(), hStack->GetXaxis()->GetXmax());
  hRatio4->Divide( hData ,hQCD,1.0,1.0);
  hRatio4->Draw();
  c2->cd(5);
  TH1F* hRatio5 = new TH1F("hRatio5", " ; ; #frac{(DATA)}{TT}",
			  hStack->GetNbinsX(), 
			  hStack->GetXaxis()->GetXmin(), hStack->GetXaxis()->GetXmax());
  hRatio5->Divide( hData ,hTTJets,1.0,1.0);
  hRatio5->Draw();

//   c2->SaveAs("plots/plotElecTau_Test.png");

  return;

}


void plotAll(){



  //////// SS

  plot("diTauVisMass"    ," ; mass (GeV) ; Events"            ,10  , 30 , "MtLeg1<40" ,"diTauCharge!=0","visibleMass");
  plot("diTauSVFitMass"  ," ; full mass (GeV) ; Events"       ,18  , 13 ,"MtLeg1<40" ,"diTauCharge!=0","fullMass");
  plot("ptL1"            ," ; e p_{T} (GeV) ; Events"         ,6   , 13 ,"MtLeg1<40" ,"diTauCharge!=0","electronPt");
  plot("ptL2"            ," ; #tau p_{T} (GeV) ; Events"      ,5   , 13 ,"MtLeg1<40" ,"diTauCharge!=0","tauPt");
  plot("abs(etaL1)"      ," ; e #eta ; Events"                ,0.2 , 13 ,"MtLeg1<40" ,"diTauCharge!=0","electronEta");
  plot("abs(etaL2)"      , " ; #tau #eta ; Events"            ,0.2 , 13 ,"MtLeg1<40" ,"diTauCharge!=0","tauEta");
  plot("MtLeg1"          ," ; M_{T}(e,MET) (GeV) ; Events"    ,10  , 13 ,"MtLeg1<999","diTauCharge!=0","Mt");
  plot("MEt"             ," ; MET (GeV) ; Events"             ,8   , 13 ,"MtLeg1<999","diTauCharge!=0","MEt");
  plot("numPV"           ," ; vertices ; Events"              ,1   , 35 ,"MtLeg1<40" ,"diTauCharge!=0","PV");

//   plot("nJets30"         ," ; jet multiplicity ; Events"       ,1   ,"MtLeg1<40" ,"diTauCharge!=0","jetMult");
//   plot("pt1"             ," ; lead jet p_{T} (GeV) ; Events"   ,10  ,"MtLeg1<40" ,"diTauCharge!=0","jet1Pt");
//   plot("pt2"             ," ; sublead jet p_{T} (GeV) ; Events",10  ,"MtLeg1<40" ,"diTauCharge!=0","jet2Pt");
//   plot("abs(eta1)"       ," ; lead jet #eta ; Events"          ,0.2 ,"MtLeg1<40" ,"diTauCharge!=0","jet1Eta");
//   plot("abs(eta2)"       ," ; sublead jet #eta ; Events"       ,0.2 ,"MtLeg1<40" ,"diTauCharge!=0","jet2Eta");

  //////// OS
//   plot("diTauVisMass"    ," ; mass (GeV) ; Events"            ,10  , 30 , "MtLeg1<40" ,"diTauCharge==0","visibleMass");
//   plot("diTauSVFitMass"  ," ; full mass (GeV) ; Events"       ,18  , 13 ,"MtLeg1<40" ,"diTauCharge==0","fullMass");
//   plot("ptL1"            ," ; e p_{T} (GeV) ; Events"         ,6   , 13 ,"MtLeg1<40" ,"diTauCharge==0","electronPt");
//   plot("ptL2"            ," ; #tau p_{T} (GeV) ; Events"      ,5   , 13 ,"MtLeg1<40" ,"diTauCharge==0","tauPt");
//   plot("abs(etaL1)"      ," ; e #eta ; Events"                ,0.2 , 13 ,"MtLeg1<40" ,"diTauCharge==0","electronEta");
//   plot("abs(etaL2)"      , " ; #tau #eta ; Events"            ,0.2 , 13 ,"MtLeg1<40" ,"diTauCharge==0","tauEta");
//   plot("MtLeg1"          ," ; M_{T}(e,MET) (GeV) ; Events"    ,10  , 13 ,"MtLeg1<999","diTauCharge==0","Mt");
//   plot("MEt"             ," ; MET (GeV) ; Events"             ,8   , 13 ,"MtLeg1<999","diTauCharge==0","MEt");
//   plot("numPV"           ," ; vertices ; Events"              ,1   , 35 ,"MtLeg1<40" ,"diTauCharge==0","PV");
//   plot("nPUVertices"     ," ; Interactions ; Events"              ,3   ,"MtLeg1<40" ,"diTauCharge==0","Num interactions");

//   plot("nJets30",   " ; jet multiplicity ; Events", 1,            "MtLeg1<40", "diTauCharge==0",  "jetMult");
//   plot("pt1",   " ; lead jet p_{T} (GeV) ; Events", 10,           "MtLeg1<40", "diTauCharge==0",  "jet1Pt");
//   plot("pt2",   " ; sublead jet p_{T} (GeV) ; Events", 10,        "MtLeg1<40", "diTauCharge==0",  "jet2Pt");
//   plot("abs(eta1)",   " ; lead jet #eta ; Events", 0.2,           "MtLeg1<40", "diTauCharge==0",  "jet1Eta");


}

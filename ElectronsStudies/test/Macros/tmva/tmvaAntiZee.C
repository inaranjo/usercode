//--------------------------------------------------------------------------------------------------
// tmvaAntiZee 
//
// Macro training and testing an MVA AntiZee.
// It takes as input an analysis eTau event.
// Signal GGFHiggs
// Background DY Zee j->Tau
// Authors: I.Naranjo
//--------------------------------------------------------------------------------------------------
#include <cstdlib>
#include <iostream> 
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
#include "TFile.h"

// #include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)

#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#endif


void TMVAClassification()
{

  TMVA::Tools::Instance();

  TString outfileName( "tmvaRoot/TMVA_AntiZee_v1.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification_AntiZee_v1", outputFile, 
					      "!V:!Silent:Color:DrawProgressBar" );
 
  factory->AddVariable( "ptL1","ptL1","     ", 'F'  );   
  factory->AddVariable( "scEtaL1","scEtaL1","     ", 'F'  );   
  factory->AddVariable( "ptL2","ptL2","     ", 'F'  );   
  factory->AddVariable( "etaL2","etaL2","     ", 'F'  );   
  factory->AddVariable( "MEtMVA","MEtMVA","     ", 'F'  );   
  factory->AddVariable( "diTauVisMass","diTauVisMass","     ", 'F'  );   

  TFile *fSig = new TFile("/data_CMS/cms/htautau/PostMoriond/NTUPLES_JetIdFix/EleTau/nTupleGGFH125_ElecTau_nominal.root","READ"); 
  TFile *fBkg = new TFile("/data_CMS/cms/htautau/PostMoriond/NTUPLES_JetIdFix/EleTau/nTupleDYJetsJetToTau_ElecTau_nominal.root","READ");

  TTree *tSig = (TTree*)fSig->Get("outTreePtOrd");
  TTree *tBkg = (TTree*)fBkg->Get("outTreePtOrd");
 
  bool useMt      = true;
  string antiWcut = useMt ? "MtLeg1MVA" : "-(pZetaMVA-1.5*pZetaVisMVA)" ; 
  float antiWsgn  = useMt ? 20. :  20. ;
  float antiWsdb  = useMt ? 70. :  40. ; 
  ///// LEPT PT ///////
  TCut lpt("ptL1>24 && TMath::Abs(etaL1)<2.1");
  TCut lID("((TMath::Abs(scEtaL1)<0.80 && mvaPOGNonTrig>0.925) || (TMath::Abs(scEtaL1)<1.479 && TMath::Abs(scEtaL1)>0.80 && mvaPOGNonTrig>0.975) || (TMath::Abs(scEtaL1)>1.479 && mvaPOGNonTrig>0.985)) && nHits<0.5");
  lpt = lpt && lID;
  TCut tpt("ptL2>20 && TMath::Abs(etaL2)<2.3");
  
  ////// TAU ISO //////
  //     TCut tiso("tightestHPSMVAWP>=0 && tightestAntiECutWP>1 && (tightestAntiEMVAWP>4 || tightestAntiEMVAWP==3)"); 
  //     TCut tiso("tightestHPSDB3HWP>=0 && tightestAntiEMVA3WP>1 && tightestAntiMu2WP>0");  // pass AntiEMediumMVA3WP AntiMuLooseMVA3WP and HPSComb3Hits
  TCut tiso("tightestHPSMVAWP>=0 && tightestAntiEMVA3WP>1 && tightestAntiMuWP>0");  // pass AntiEMediumMVA3WP AntiMuLooseMVA3WP and HPSComb3Hits
  TCut ltiso("tightestHPSMVAWP>-99 && tightestAntiECutWP<1 ");
  TCut mtiso("hpsMVA>0.7");
  
  ////// E ISO ///////
  TCut liso("combRelIsoLeg1DBetav2<0.10");
  TCut laiso("combRelIsoLeg1DBetav2>0.20 && combRelIsoLeg1DBetav2<0.50");
  TCut lliso("combRelIsoLeg1DBetav2<0.30");
  
  ////// EVENT WISE //////
  TCut lveto("elecFlag==0"); //elecFlag==0
  TCut SS("diTauCharge!=0");
  TCut OS("diTauCharge==0");
  TCut pZ( Form("((%s)<%f)",antiWcut.c_str(),antiWsgn));
  TCut apZ(Form("((%s)>%f)",antiWcut.c_str(),antiWsdb));
  
  TCut apZ2(Form("((%s)>60 && (%s)<120)",antiWcut.c_str(),antiWcut.c_str()));
  TCut hltevent("vetoEvent==0 && pairIndex[2]<1 && HLTx==1");
  TCut hltmatch("HLTmatch==1");
  
  ////// CATEGORIES ///
  TCut vbf("nJets30>=2 && Mjj>500 && Deta>3.5 && isVetoInJets!=1");
//   TCut boost("nJets30>0 && pt1>30 && nJets20BTagged<1");
//   boost = boost && !vbf;

  TCut boost("nJets30>0 && pt1>30 && nJets20BTagged<1");

  TCut novbf("nJets30<1 && nJets20BTagged==0");
  
  bool removeMtCut     = false;  
  TCut MtCut       = removeMtCut     ? "(etaL1<999)" : pZ;
  TCut diTauCharge = OS; 
  
  TCut sbin;
  
  TCut sbinTmp("");
//   if(category=="inclusive") 
//     sbinTmp = "etaL1<999";
//   else if(category=="vbf")
//     sbinTmp = vbf;
//   else if(category=="novbf")
//     sbinTmp = novbf;
//   else if(category=="boost")
  sbinTmp = boost;
  
  sbin =  sbinTmp && lpt && tpt && tiso && liso && lveto && diTauCharge  && MtCut  && hltevent && hltmatch ;
  
  cout<<"Selection : "<<sbin<<endl;

  ////////////////// compute the weights
  factory->AddSignalTree( tSig );
  factory->AddBackgroundTree( tBkg );
//   factory->SetWeightExpression("puWeight");

  factory->PrepareTrainingAndTestTree( sbin,sbin,
				       "nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
  
  // Boosted Decision Trees with gradient boosting
  factory->BookMethod( TMVA::Types::kBDT, "BDTG","!H:!V:NTrees=400:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrength=50:NNodesMax=5" );
  
  factory->TrainAllMethods();
  
  factory->TestAllMethods();

  factory->EvaluateAllMethods();  

  outputFile->Close();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;      
  
  delete factory;

}

void TMVAAllClassification(){
  TMVAClassification();

}

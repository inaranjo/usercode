//--------------------------------------------------------------------------------------------------
// tmvaAntiElectronOptimization 
//
// Macro optimizing the MVA cuts for the AntiElectron MVA.
// It takes as input the MVA output for the 4*2 categories.
//
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

//#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)

#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"
#endif


void TMVAClassification(string Discr_ = "")
{

  TMVA::Tools::Instance();

//   TString outfileName("tmvaRoot/TMVAOptimization_v4"+Discr_+".root" );
//   TString outfileName("tmvaRoot/TMVAOptimization_v4-EtaAtEcal"+Discr_+".root" );
  TString outfileName("tmvaRoot/TMVAOptimization_v5-MC"+Discr_+".root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

//   TMVA::Factory *factory = new TMVA::Factory("TMVAClassificationOptimized_v4"+Discr_, outputFile, 
// 					      "!V:!Silent:Color:DrawProgressBar" );
//   TMVA::Factory *factory = new TMVA::Factory("TMVAClassificationOptimized_v4-EtaAtEcal"+Discr_, outputFile, 
// 					      "!V:!Silent:Color:DrawProgressBar" );
  TMVA::Factory *factory = new TMVA::Factory("TMVAClassificationOptimized_v5-MC"+Discr_, outputFile, 
					      "!V:!Silent:Color:DrawProgressBar" );

  factory->AddVariable( "NoEleMatch_Barrel", "mva_NoEleMatch_Barrel","",'F');
  factory->AddVariable( "woG_Barrel", "mva_woG_Barrel","",'F');
  factory->AddVariable( "wGwoGSF_Barrel", "mva_wGwoGSF_Barrel","",'F');
  factory->AddVariable( "wGwGSFwoPFMVA_Barrel", "mva_wGwGSFwoPFMVA_Barrel","",'F');
  if (Discr_ == "")factory->AddVariable( "wGwGSFwPFMVA_Barrel", "mva_wGwGSFwPFMVA_Barrel","",'F');
  factory->AddVariable( "NoEleMatch_Endcap", "mva_NoEleMatch_Endcap","",'F');
  factory->AddVariable( "woG_Endcap", "mva_woG_Endcap","",'F');
  factory->AddVariable( "wGwoGSF_Endcap", "mva_wGwoGSF_Endcap","",'F');
  factory->AddVariable( "wGwGSFwoPFMVA_Endcap", "mva_wGwGSFwoPFMVA_Endcap","",'F');
  if (Discr_ == "")factory->AddVariable( "wGwGSFwPFMVA_Endcap", "mva_wGwGSFwPFMVA_Endcap","",'F');

//   TFile *fSignal = new TFile(Form("/data_CMS/cms/ivo/AntiEMVA/Trees/MVAOutput/MvaOutput_AntiEMVA_v4%s_Signal.root",Discr_.data()),"READ"); 
//   TFile *fSignal = new TFile(Form("/data_CMS/cms/ivo/AntiEMVA/Trees/MVAOutput/MvaOutput_AntiEMVA-EtaAtEcal_v4%s_Signal.root",Discr_.data()),"READ"); 
  TFile *fSignal = new TFile(Form("/data_CMS/cms/ivo/AntiEMVA/Trees/MVAOutput/MvaOutput_AntiEMVA-v5%s_Signal.root",Discr_.data()),"READ"); 
  if(fSignal->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
//   TFile *fBackgrd = new TFile(Form("/data_CMS/cms/ivo/AntiEMVA/Trees/MVAOutput/MvaOutput_AntiEMVA_v4%s_Backgrd.root",Discr_.data()),"READ"); 
//   TFile *fBackgrd = new TFile(Form("/data_CMS/cms/ivo/AntiEMVA/Trees/MVAOutput/MvaOutput_AntiEMVA-EtaAtEcal_v4%s_Backgrd.root",Discr_.data()),"READ"); 
  TFile *fBackgrd = new TFile(Form("/data_CMS/cms/ivo/AntiEMVA/Trees/MVAOutput/MvaOutput_AntiEMVA-v5%s_Backgrd.root",Discr_.data()),"READ"); 
  if(fBackgrd->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
  TTree*inTreeSgn = (TTree*)fSignal->Get("tree");
  TTree*inTreeBkg = (TTree*)fBackgrd->Get("tree");

  TCut myCut = "";

  factory->AddSignalTree( inTreeSgn );
  factory->AddBackgroundTree(inTreeBkg );

  // v5, v5-Backup
//   factory->PrepareTrainingAndTestTree(myCut,myCut,"nTrain_Signal=100000:nTrain_Background=300000:nTest_Signal=100000:nTest_Background=300000:SplitMode=Random:NormMode=NumEvents:!V" );
//   factory->PrepareTrainingAndTestTree(myCut,myCut,"nTrain_Signal=50000:nTrain_Background=150000:nTest_Signal=50000:nTest_Background=150000:SplitMode=Random:NormMode=NumEvents:!V" );
//23Apr2012, v5-iter1   
// factory->PrepareTrainingAndTestTree(myCut,myCut,"nTrain_Signal=50000:nTrain_Background=200000:nTest_Signal=50000:nTest_Background=200000:SplitMode=Random:NormMode=NumEvents:!V" );

factory->PrepareTrainingAndTestTree(myCut,myCut,"nTrain_Signal=3000000:nTrain_Background=3000000:nTest_Signal=1000000:nTest_Background=1000000:SplitMode=Random:NormMode=NumEvents:!V" );
  
  factory->BookMethod( TMVA::Types::kCuts, "Cuts", 
		       //"!H:!V:FitMethod=GA:EffSel" );
 		       "!H:!V:FitMethod=MC:EffSel:CutRangeMin[0]=-1.:CutRangeMax[0]=1:CutRangeMin[1]=-1.:CutRangeMax[1]=1.:CutRangeMin[2]=-1:CutRangeMax[2]=1.:CutRangeMin[3]=-1:CutRangeMax[3]=1:CutRangeMin[4]=-1:CutRangeMax[4]=1:CutRangeMin[5]=-1:CutRangeMax[5]=1:CutRangeMin[6]=-1:CutRangeMax[6]=1:CutRangeMin[7]=-1:CutRangeMax[7]=1:CutRangeMin[8]=-1:CutRangeMax[8]=1:CutRangeMin[9]=-1:CutRangeMax[9]=1" );
//  		       "!H:!V:FitMethod=GA:EffSel:CutRangeMin[0]=-1.:CutRangeMax[0]=1:CutRangeMin[1]=-1.:CutRangeMax[1]=1.:CutRangeMin[2]=-1:CutRangeMax[2]=1.:CutRangeMin[3]=-1:CutRangeMax[3]=1:CutRangeMin[4]=-1:CutRangeMax[4]=1:CutRangeMin[5]=-1:CutRangeMax[5]=1");//:CutRangeMin[6]=-1:CutRangeMax[6]=1:CutRangeMin[7]=-1:CutRangeMax[7]=1" );
//  		       "!H:!V:FitMethod=GA:EffSel:CutRangeMin[0]=-1.:CutRangeMax[0]=1:CutRangeMin[1]=-1.:CutRangeMax[1]=1.:CutRangeMin[2]=-1:CutRangeMax[2]=1.:CutRangeMin[3]=-1:CutRangeMax[3]=1:CutRangeMin[4]=-1:CutRangeMax[4]=1:CutRangeMin[5]=-1:CutRangeMax[5]=1:CutRangeMin[6]=-1:CutRangeMax[6]=1:CutRangeMin[7]=-1:CutRangeMax[7]=1" );
  factory->TrainAllMethods();
  
  factory->TestAllMethods();

  factory->EvaluateAllMethods();  

  TMVA::MethodCuts* cuts = dynamic_cast<TMVA::MethodCuts*>(factory->GetMethod("Cuts"));    
  Double_t epsilon = 0.0001;
  for (Double_t eff=0.05; eff<1; eff += 0.01) {
    cuts->PrintCuts( eff+epsilon );
  } 

  outputFile->Close();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;      
  
  delete factory;

  //if (!gROOT->IsBatch()) TMVAGui( outfileName );

}

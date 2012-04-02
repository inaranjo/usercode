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
#endif


void TMVAClassification(std::string Discr_ = ""	)
{

  TMVA::Tools::Instance();

  TString outfileName("tmvaRoot/TMVAOptimization"+Discr_+".root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory("TMVAClassificationOptimized"+Discr_, outputFile, 
					      "!V:!Silent:Color:DrawProgressBar" );

  factory->AddVariable( "woG_Barrel", "mva_woG_Barrel","",'F');
  factory->AddVariable( "wGwoGSF_Barrel", "mva_wGwoGSF_Barrel","",'F');
  factory->AddVariable( "wGwGSFwoPFMVA_Barrel", "mva_wGwGSFwoPFMVA_Barrel","",'F');
  factory->AddVariable( "wGwGSFwPFMVA_Barrel", "mva_wGwGSFwPFMVA_Barrel","",'F');
  factory->AddVariable( "woG_Endcap", "mva_woG_Endcap","",'F');
  factory->AddVariable( "wGwoGSF_Endcap", "mva_wGwoGSF_Endcap","",'F');
  factory->AddVariable( "wGwGSFwoPFMVA_Endcap", "mva_wGwGSFwoPFMVA_Endcap","",'F');
  factory->AddVariable( "wGwGSFwPFMVA_Endcap", "mva_wGwGSFwPFMVA_Endcap","",'F');

  TFile *fSignal = new TFile(Form("../MvaOutput_AntiEMVA%s_Signal.root",Discr_.data()),"READ"); 
  if(fSignal->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
  TFile *fBackgrd = new TFile(Form("../MvaOutput_AntiEMVA%s_Backgrd.root",Discr_.data()),"READ"); 
  if(fBackgrd->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
  TTree*inTreeSgn = (TTree*)fSignal->Get("tree");
  TTree*inTreeBkg = (TTree*)fBackgrd->Get("tree");

  TCut myCut = "";

  factory->AddSignalTree( inTreeSgn );
  factory->AddBackgroundTree(inTreeBkg );

  factory->PrepareTrainingAndTestTree(myCut,myCut,"nTrain_Signal=20000:nTrain_Background=100000:nTest_Signal=20000:nTest_Background=100000:SplitMode=Random:NormMode=NumEvents:!V" );
  
  factory->BookMethod( TMVA::Types::kCuts, "Cuts", 
		       //"!H:!V:FitMethod=GA:EffSel" );
 		       "!H:!V:FitMethod=GA:EffSel:CutRangeMin[0]=-1.:CutRangeMax[0]=1:CutRangeMin[1]=-1.:CutRangeMax[1]=1.:CutRangeMin[2]=-1:CutRangeMax[2]=1.:CutRangeMin[3]=-1:CutRangeMax[3]=1:CutRangeMin[4]=-1:CutRangeMax[4]=1:CutRangeMin[5]=-1:CutRangeMax[5]=1:CutRangeMin[6]=-1:CutRangeMax[6]=1:CutRangeMin[7]=-1:CutRangeMax[7]=1" );
  

  factory->TrainAllMethods();
  
  factory->TestAllMethods();

  factory->EvaluateAllMethods();  

  outputFile->Close();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;      
  
  delete factory;

  //if (!gROOT->IsBatch()) TMVAGui( outfileName );

}

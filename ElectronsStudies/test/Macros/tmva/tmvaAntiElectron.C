//--------------------------------------------------------------------------------------------------
// tmvaAntiElectron 
//
// Macro training and testing an MVA AntiElectronID for Taus.
// It takes as input a PFTau-GsfElectron pair matched to a Tau Hadron (Signal), and matched with a Gen Electron
// for the background tree.
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

// #include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)

#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#endif


void TMVAClassification(std::string Cat_ = "_All",
			std::string Sel_ = "Barrel",
			std::string Discr_ = ""			
			)
{

  TMVA::Tools::Instance();

//   TString outfileName( "tmvaRoot/TMVA_v4"+Discr_+Cat_+"_"+Sel_+".root" );
  TString outfileName( "tmvaRoot/TMVA_v5"+Discr_+Cat_+"_"+Sel_+".root" );
//   TString outfileName( "tmvaRoot/TMVA_v4_EtaAtEcal"+Discr_+Cat_+"_"+Sel_+".root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

//   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification_v4"+Discr_+Cat_+"_"+Sel_, outputFile, 
// 					      "!V:!Silent:Color:DrawProgressBar" );
//   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification_v4_EtaAtEcal"+Discr_+Cat_+"_"+Sel_, outputFile, 
// 					      "!V:!Silent:Color:DrawProgressBar" );
  TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification_v5"+Discr_+Cat_+"_"+Sel_, outputFile, 
					      "!V:!Silent:Color:DrawProgressBar" );
 
 if(Cat_.find("_All")!=string::npos){    
   factory->AddVariable( "Elec_Pt","Elec_Pt","     ", 'F'  );   
   factory->AddVariable( "Elec_AbsEta","ElecAbsEta","     ", 'F'  );
   factory->AddVariable( "Elec_EtotOverPin","EtotOverPin","     ", 'F'  );
   factory->AddVariable( "Elec_EeOverPout","EeOverPout","     ", 'F'  );
   factory->AddVariable( "Elec_EgammaOverPdif","EgammaOverPdif","     ", 'F'  );
   factory->AddVariable( "Elec_EarlyBrem","EarlyBrem","     ", 'I'  );
   factory->AddVariable( "Elec_LateBrem","LateBrem","     ", 'I'  );
// //    factory->AddVariable( "Elec_Logsihih","Elec_Logsihih","     ", 'F'  );
// //    factory->AddVariable( "Elec_DeltaEta","Elec_DeltaEta","     ", 'F'  );
   factory->AddVariable( "Elec_Fbrem", "Fbrem","     ", 'F'  );
   factory->AddVariable( "Elec_Chi2KF", "Chi2KF","     ", 'F'  );
   factory->AddVariable( "Elec_Chi2GSF","Elec_Chi2GSF","     ", 'F'  ); 
   factory->AddVariable( "Elec_NumHits","NumHits","     ", 'I'  );
   factory->AddVariable( "Elec_GSFTrackResol","Elec_GSFTrackResol","     ", 'F'  );
   factory->AddVariable( "Elec_GSFTracklnPt","Elec_GSFTracklnPt","     ", 'F'  ); 
   factory->AddVariable( "Elec_GSFTrackEta","Elec_GSFTrackEta","     ", 'F'  );  
   
   factory->AddVariable( "Tau_AbsEta","PFTauEta","     ", 'F'  );
   factory->AddVariable( "Tau_Pt","TauPt","     ", 'F'  );            
   factory->AddVariable( "Tau_HasGsf","TauHasGsf","     ", 'I'  );
   factory->AddVariable( "Tau_EmFraction","TauEmFraction","     ", 'F'  );
   factory->AddVariable( "Tau_NumChargedCands","Tau_NumChargedCands","     ", 'I'  );
   factory->AddVariable( "Tau_NumGammaCands","TauNumGammaCands","     ", 'I'  );
   factory->AddVariable( "Tau_HadrHoP","TauHadrHoP","     ", 'F'  );
   factory->AddVariable( "Tau_HadrEoP","TauHadrEoP","     ", 'F'  );
   factory->AddVariable( "Tau_VisMass","TauVisMass","     ", 'F'  );
   factory->AddVariable( "Tau_GammaEtaMom","Tau_GammaEtaMom","     ", 'F'  );    
   factory->AddVariable( "Tau_GammaPhiMom","Tau_GammaPhiMom","     ", 'F'  ); 
   factory->AddVariable( "Tau_GammaEnFrac","GammaEnFrac","     ", 'F'  );
 }
 if(Cat_.find("_NoEleMatch")!=string::npos){    
//    factory->AddVariable( "Elec_Pt","Elec_Pt","     ", 'F'  );   
//    factory->AddVariable( "Elec_AbsEta","ElecAbsEta","     ", 'F'  );
//    factory->AddVariable( "Elec_EtotOverPin","EtotOverPin","     ", 'F'  );
//    factory->AddVariable( "Elec_EeOverPout","EeOverPout","     ", 'F'  );
//    factory->AddVariable( "Elec_EgammaOverPdif","EgammaOverPdif","     ", 'F'  );
//    factory->AddVariable( "Elec_EarlyBrem","EarlyBrem","     ", 'I'  );
//    factory->AddVariable( "Elec_LateBrem","LateBrem","     ", 'I'  );
// //    factory->AddVariable( "Elec_Logsihih","Elec_Logsihih","     ", 'F'  );
// //    factory->AddVariable( "Elec_DeltaEta","Elec_DeltaEta","     ", 'F'  );
//    factory->AddVariable( "Elec_Fbrem", "Fbrem","     ", 'F'  );
//    factory->AddVariable( "Elec_Chi2KF", "Chi2KF","     ", 'F'  );
//    factory->AddVariable( "Elec_Chi2GSF","Elec_Chi2GSF","     ", 'F'  ); 
//    factory->AddVariable( "Elec_NumHits","NumHits","     ", 'I'  );
//    factory->AddVariable( "Elec_GSFTrackResol","Elec_GSFTrackResol","     ", 'F'  );
//    factory->AddVariable( "Elec_GSFTracklnPt","Elec_GSFTracklnPt","     ", 'F'  ); 
//    factory->AddVariable( "Elec_GSFTrackEta","Elec_GSFTrackEta","     ", 'F'  );  
   
//for v4   factory->AddVariable( "Tau_AbsEta","PFTauEta","     ", 'F'  );
   factory->AddVariable( "Tau_EtaAtEcalEntrance","PFTauEtaAtEcalEntrance","     ", 'F'  );
//    factory->AddVariable( "Tau_Eta","PFTauEta","     ", 'F'  );
   factory->AddVariable( "Tau_Pt","TauPt","     ", 'F'  );            
//for v3     factory->AddVariable( "Tau_HasGsf","TauHasGsf","     ", 'I'  );
   factory->AddVariable( "Tau_EmFraction","TauEmFraction","     ", 'F'  );
//for v3    factory->AddVariable( "Tau_NumChargedCands","Tau_NumChargedCands","     ", 'I'  );
   factory->AddVariable( "Tau_NumGammaCands","TauNumGammaCands","     ", 'I'  );
   factory->AddVariable( "Tau_HadrHoP","TauHadrHoP","     ", 'F'  );
   factory->AddVariable( "Tau_HadrEoP","TauHadrEoP","     ", 'F'  );
   factory->AddVariable( "Tau_VisMass","TauVisMass","     ", 'F'  );
   factory->AddVariable( "Tau_GammaEtaMom","Tau_GammaEtaMom","     ", 'F'  );    
   factory->AddVariable( "Tau_GammaPhiMom","Tau_GammaPhiMom","     ", 'F'  ); 
   factory->AddVariable( "Tau_GammaEnFrac","GammaEnFrac","     ", 'F'  );
 }
 if(Cat_.find("_woG")!=string::npos){
 //for v2   factory->AddVariable( "Elec_Pt","Elec_Pt","     ", 'F'  );   
 //for v2   factory->AddVariable( "Elec_AbsEta","ElecAbsEta","     ", 'F'  );
   factory->AddVariable( "Elec_EtotOverPin","EtotOverPin","     ", 'F'  );
//    factory->AddVariable( "Elec_EeOverPout","EeOverPout","     ", 'F'  );
//    factory->AddVariable( "Elec_EgammaOverPdif","EgammaOverPdif","     ", 'F'  );
//for v2    factory->AddVariable( "Elec_EarlyBrem","EarlyBrem","     ", 'I'  );
   factory->AddVariable( "Elec_LateBrem","LateBrem","     ", 'I'  );
// //    factory->AddVariable( "Elec_Logsihih","Elec_Logsihih","     ", 'F'  );
// //    factory->AddVariable( "Elec_DeltaEta","Elec_DeltaEta","     ", 'F'  );
   factory->AddVariable( "Elec_Fbrem", "Fbrem","     ", 'F'  );
   factory->AddVariable( "Elec_Chi2KF", "Chi2KF","     ", 'F'  );
//    factory->AddVariable( "Elec_Chi2GSF","Elec_Chi2GSF","     ", 'F'  ); 
//for v2    factory->AddVariable( "Elec_NumHits","NumHits","     ", 'I'  );
   factory->AddVariable( "Elec_GSFTrackResol","Elec_GSFTrackResol","     ", 'F'  );
   factory->AddVariable( "Elec_GSFTracklnPt","Elec_GSFTracklnPt","     ", 'F'  ); 
   factory->AddVariable( "Elec_GSFTrackEta","Elec_GSFTrackEta","     ", 'F'  );  
   
//for v4   factory->AddVariable( "Tau_AbsEta","PFTauEta","     ", 'F'  );
   factory->AddVariable( "Tau_EtaAtEcalEntrance","PFTauEtaAtEcalEntrance","     ", 'F'  );
//    factory->AddVariable( "Tau_Eta","PFTauEta","     ", 'F'  );
   factory->AddVariable( "Tau_Pt","TauPt","     ", 'F'  );            
//for v3    factory->AddVariable( "Tau_HasGsf","TauHasGsf","     ", 'I'  );
   factory->AddVariable( "Tau_EmFraction","TauEmFraction","     ", 'F'  );
//    factory->AddVariable( "Tau_NumChargedCands","Tau_NumChargedCands","     ", 'I'  );
//    factory->AddVariable( "Tau_NumGammaCands","TauNumGammaCands","     ", 'I'  );
   factory->AddVariable( "Tau_HadrHoP","TauHadrHoP","     ", 'F'  );
   factory->AddVariable( "Tau_HadrEoP","TauHadrEoP","     ", 'F'  );
   factory->AddVariable( "Tau_VisMass","TauVisMass","     ", 'F'  );
//    factory->AddVariable( "Tau_GammaEtaMom","Tau_GammaEtaMom","     ", 'F'  );    
//    factory->AddVariable( "Tau_GammaPhiMom","Tau_GammaPhiMom","     ", 'F'  ); 
//    factory->AddVariable( "Tau_GammaEnFrac","GammaEnFrac","     ", 'F'  );
 }
 if(Cat_.find("_wGwoGSF")!=string::npos){
  //for v2  factory->AddVariable( "Elec_Pt","Elec_Pt","     ", 'F'  );   
  //for v2  factory->AddVariable( "Elec_AbsEta","ElecAbsEta","     ", 'F'  );
   factory->AddVariable( "Elec_EtotOverPin","EtotOverPin","     ", 'F'  );
//    factory->AddVariable( "Elec_EeOverPout","EeOverPout","     ", 'F'  );
   factory->AddVariable( "Elec_EgammaOverPdif","EgammaOverPdif","     ", 'F'  );
//for v2    factory->AddVariable( "Elec_EarlyBrem","EarlyBrem","     ", 'I'  );
   factory->AddVariable( "Elec_LateBrem","LateBrem","     ", 'I'  );
// //    factory->AddVariable( "Elec_Logsihih","Elec_Logsihih","     ", 'F'  );
// //    factory->AddVariable( "Elec_DeltaEta","Elec_DeltaEta","     ", 'F'  );
   factory->AddVariable( "Elec_Fbrem", "Fbrem","     ", 'F'  );
//    factory->AddVariable( "Elec_Chi2KF", "Chi2KF","     ", 'F'  );
   factory->AddVariable( "Elec_Chi2GSF","Elec_Chi2GSF","     ", 'F'  ); 
   factory->AddVariable( "Elec_NumHits","NumHits","     ", 'I'  );
   factory->AddVariable( "Elec_GSFTrackResol","Elec_GSFTrackResol","     ", 'F'  );
   factory->AddVariable( "Elec_GSFTracklnPt","Elec_GSFTracklnPt","     ", 'F'  ); 
   factory->AddVariable( "Elec_GSFTrackEta","Elec_GSFTrackEta","     ", 'F'  );  
   
//for v4   factory->AddVariable( "Tau_AbsEta","PFTauEta","     ", 'F'  );
   factory->AddVariable( "Tau_EtaAtEcalEntrance","PFTauEtaAtEcalEntrance","     ", 'F'  );
//    factory->AddVariable( "Tau_Eta","PFTauEta","     ", 'F'  );
   factory->AddVariable( "Tau_Pt","TauPt","     ", 'F'  );            
//    factory->AddVariable( "Tau_HasGsf","TauHasGsf","     ", 'I'  );
   factory->AddVariable( "Tau_EmFraction","TauEmFraction","     ", 'F'  );
//    factory->AddVariable( "Tau_NumChargedCands","Tau_NumChargedCands","     ", 'I'  );
   factory->AddVariable( "Tau_NumGammaCands","TauNumGammaCands","     ", 'I'  );
   factory->AddVariable( "Tau_HadrHoP","TauHadrHoP","     ", 'F'  );
   factory->AddVariable( "Tau_HadrEoP","TauHadrEoP","     ", 'F'  );
   factory->AddVariable( "Tau_VisMass","TauVisMass","     ", 'F'  );
   factory->AddVariable( "Tau_GammaEtaMom","Tau_GammaEtaMom","     ", 'F'  );    
   factory->AddVariable( "Tau_GammaPhiMom","Tau_GammaPhiMom","     ", 'F'  ); 
   factory->AddVariable( "Tau_GammaEnFrac","GammaEnFrac","     ", 'F'  );   
 }
 if(Cat_.find("_wGwGSFwoPFMVA")!=string::npos){
//for v2    factory->AddVariable( "Elec_Pt","Elec_Pt","     ", 'F'  );   
//for v2    factory->AddVariable( "Elec_AbsEta","ElecAbsEta","     ", 'F'  );
//    factory->AddVariable( "Elec_EtotOverPin","EtotOverPin","     ", 'F'  );
//    factory->AddVariable( "Elec_EeOverPout","EeOverPout","     ", 'F'  );
//    factory->AddVariable( "Elec_EgammaOverPdif","EgammaOverPdif","     ", 'F'  );
//    factory->AddVariable( "Elec_EarlyBrem","EarlyBrem","     ", 'I'  );
//    factory->AddVariable( "Elec_LateBrem","LateBrem","     ", 'I'  );
// //    factory->AddVariable( "Elec_Logsihih","Elec_Logsihih","     ", 'F'  );
// //    factory->AddVariable( "Elec_DeltaEta","Elec_DeltaEta","     ", 'F'  );
   factory->AddVariable( "Elec_Fbrem", "Fbrem","     ", 'F'  );
   factory->AddVariable( "Elec_Chi2KF", "Chi2KF","     ", 'F'  );
   factory->AddVariable( "Elec_Chi2GSF","Elec_Chi2GSF","     ", 'F'  ); 
   factory->AddVariable( "Elec_NumHits","NumHits","     ", 'I'  );
   factory->AddVariable( "Elec_GSFTrackResol","Elec_GSFTrackResol","     ", 'F'  );
   factory->AddVariable( "Elec_GSFTracklnPt","Elec_GSFTracklnPt","     ", 'F'  ); 
   factory->AddVariable( "Elec_GSFTrackEta","Elec_GSFTrackEta","     ", 'F'  );  
   
//for v4   factory->AddVariable( "Tau_AbsEta","PFTauEta","     ", 'F'  );
   factory->AddVariable( "Tau_EtaAtEcalEntrance","PFTauEtaAtEcalEntrance","     ", 'F'  );
//    factory->AddVariable( "Tau_Eta","PFTauEta","     ", 'F'  );
   factory->AddVariable( "Tau_Pt","TauPt","     ", 'F'  );            
//    factory->AddVariable( "Tau_HasGsf","TauHasGsf","     ", 'I'  );
   factory->AddVariable( "Tau_EmFraction","TauEmFraction","     ", 'F'  );
//    factory->AddVariable( "Tau_NumChargedCands","Tau_NumChargedCands","     ", 'I'  );
   factory->AddVariable( "Tau_NumGammaCands","TauNumGammaCands","     ", 'I'  );
   factory->AddVariable( "Tau_HadrHoP","TauHadrHoP","     ", 'F'  );
   factory->AddVariable( "Tau_HadrEoP","TauHadrEoP","     ", 'F'  );
   factory->AddVariable( "Tau_VisMass","TauVisMass","     ", 'F'  );
   factory->AddVariable( "Tau_GammaEtaMom","Tau_GammaEtaMom","     ", 'F'  );    
   factory->AddVariable( "Tau_GammaPhiMom","Tau_GammaPhiMom","     ", 'F'  ); 
   factory->AddVariable( "Tau_GammaEnFrac","GammaEnFrac","     ", 'F'  );
 }
 if(Cat_.find("_wGwGSFwPFMVA")!=string::npos){
//for v2    factory->AddVariable( "Elec_Pt","Elec_Pt","     ", 'F'  );   
//for v2    factory->AddVariable( "Elec_AbsEta","ElecAbsEta","     ", 'F'  );
   factory->AddVariable( "Elec_EtotOverPin","EtotOverPin","     ", 'F'  );
   factory->AddVariable( "Elec_EeOverPout","EeOverPout","     ", 'F'  );
//    factory->AddVariable( "Elec_EgammaOverPdif","EgammaOverPdif","     ", 'F'  );
//for v3     factory->AddVariable( "Elec_EarlyBrem","EarlyBrem","     ", 'I'  );
   factory->AddVariable( "Elec_LateBrem","LateBrem","     ", 'I'  );
// //    factory->AddVariable( "Elec_Logsihih","Elec_Logsihih","     ", 'F'  );
// // //    factory->AddVariable( "Elec_DeltaEta","Elec_DeltaEta","     ", 'F'  );
//    factory->AddVariable( "Elec_Fbrem", "Fbrem","     ", 'F'  );
//    factory->AddVariable( "Elec_Chi2KF", "Chi2KF","     ", 'F'  );
   factory->AddVariable( "Elec_Chi2GSF","Elec_Chi2GSF","     ", 'F'  ); 
   factory->AddVariable( "Elec_NumHits","NumHits","     ", 'I'  );
   factory->AddVariable( "Elec_GSFTrackResol","Elec_GSFTrackResol","     ", 'F'  );
   factory->AddVariable( "Elec_GSFTracklnPt","Elec_GSFTracklnPt","     ", 'F'  ); 
   factory->AddVariable( "Elec_GSFTrackEta","Elec_GSFTrackEta","     ", 'F'  );  
   
//for v4   factory->AddVariable( "Tau_AbsEta","PFTauEta","     ", 'F'  );
   factory->AddVariable( "Tau_EtaAtEcalEntrance","PFTauEtaAtEcalEntrance","     ", 'F'  );
//    factory->AddVariable( "Tau_Eta","PFTauEta","     ", 'F'  );
   factory->AddVariable( "Tau_Pt","TauPt","     ", 'F'  );            
//    factory->AddVariable( "Tau_HasGsf","TauHasGsf","     ", 'I'  );
   factory->AddVariable( "Tau_EmFraction","TauEmFraction","     ", 'F'  );
//    factory->AddVariable( "Tau_NumChargedCands","Tau_NumChargedCands","     ", 'I'  );
   factory->AddVariable( "Tau_NumGammaCands","TauNumGammaCands","     ", 'I'  );
   factory->AddVariable( "Tau_HadrHoP","TauHadrHoP","     ", 'F'  );
   factory->AddVariable( "Tau_HadrEoP","TauHadrEoP","     ", 'F'  );
   factory->AddVariable( "Tau_VisMass","TauVisMass","     ", 'F'  );
   factory->AddVariable( "Tau_GammaEtaMom","Tau_GammaEtaMom","     ", 'F'  );    
   factory->AddVariable( "Tau_GammaPhiMom","Tau_GammaPhiMom","     ", 'F'  ); 
   factory->AddVariable( "Tau_GammaEnFrac","GammaEnFrac","     ", 'F'  );
 }


  TFile *fTau = new TFile(Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA%s%s_Tau.root",Discr_.data(),Cat_.data()),"READ"); 
  TFile *fEle = new TFile(Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_AntiEMVA%s%s_Elec.root",Discr_.data(),Cat_.data()),"READ");


  TTree *tTau = (TTree*)fTau->Get("tree");
  TTree *tEle = (TTree*)fEle->Get("tree");
 
  TCut myCut = "";

  if( Sel_.find("Barrel") !=string::npos)  myCut = myCut && TCut("Elec_AbsEta<1.479 && Tau_EtaAtEcalEntrance<1.479 && Tau_EtaAtEcalEntrance>-1.479");
  if(Sel_.find("Endcap") !=string::npos)  myCut = myCut && TCut("Elec_AbsEta>1.479 && Elec_AbsEta<3.0 && (Tau_EtaAtEcalEntrance>1.479 && Tau_EtaAtEcalEntrance<2.3) || (Tau_EtaAtEcalEntrance>-2.3 && Tau_EtaAtEcalEntrance<-1.479)");
//   if( Sel_.find("Barrel") !=string::npos)  myCut = myCut && TCut("Elec_AbsEta<1.479 && Tau_Eta<1.479 && Tau_Eta>-1.479");
//   if(Sel_.find("Endcap") !=string::npos)  myCut = myCut && TCut("Elec_AbsEta>1.479 && Elec_AbsEta<3.0 && (Tau_Eta>1.479 && Tau_Eta<2.3) || (Tau_Eta>-2.3 && Tau_Eta<-1.479)");


  ////////////////// compute the weights

  factory->AddSignalTree( tTau );
  factory->AddBackgroundTree( tEle );
//   factory->SetWeightExpression("puWeight");

  factory->PrepareTrainingAndTestTree( myCut,myCut,
				       "nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
  
  
//   factory->BookMethod( TMVA::Types::kBDT, "BDT", 
// 		       "!H:!V:NTrees=400:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
  
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
  TMVAClassification("_NoEleMatch","Barrel","");
  TMVAClassification("_woG","Barrel","");
  TMVAClassification("_wGwoGSF","Barrel","");
  TMVAClassification("_wGwGSFwoPFMVA","Barrel","");
  TMVAClassification("_wGwGSFwPFMVA","Barrel","");

  TMVAClassification("_NoEleMatch","Endcap","");
  TMVAClassification("_woG","Endcap","");
  TMVAClassification("_wGwoGSF","Endcap","");
  TMVAClassification("_wGwGSFwoPFMVA","Endcap","");
  TMVAClassification("_wGwGSFwPFMVA","Endcap","");



//   TMVAClassification("_NoEleMatch","Barrel","-AntiEMed");
//   TMVAClassification("_woG","Barrel","-AntiEMed");
//   TMVAClassification("_wGwoGSF","Barrel","-AntiEMed");
//   TMVAClassification("_wGwGSFwoPFMVA","Barrel","-AntiEMed");
// //   TMVAClassification("_wGwGSFwPFMVA","Barrel","-AntiEMed");

//   TMVAClassification("_NoEleMatch","Endcap","-AntiEMed");
//   TMVAClassification("_woG","Endcap","-AntiEMed");
//   TMVAClassification("_wGwoGSF","Endcap","-AntiEMed");
//   TMVAClassification("_wGwGSFwoPFMVA","Endcap","-AntiEMed");
// //   TMVAClassification("_wGwGSFwPFMVA","Endcap","-AntiEMed");
}

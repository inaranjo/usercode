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


void TMVAClassification(std::string Cat_ = "All",
			std::string Sel_ = "Barrel"
			
			)
{

  TMVA::Tools::Instance();

  TString outfileName( "tmvaRoot/TMVA_"+Cat_+"_"+Sel_+".root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification_"+Cat_+"_"+Sel_, outputFile, 
					      "!V:!Silent:Color:DrawProgressBar" );

 
 if(Cat_.find("All")!=string::npos){
    factory->AddVariable( "Elec_Fbrem", "Fbrem",            "     "     , 'F'  );
    factory->AddVariable( "Elec_Chi2KF", "Chi2KF",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EarlyBrem","EarlyBrem",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EeOverPout","EeOverPout",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EgammaOverPdif","EgammaOverPdif",            "     "     , 'F'  );
    factory->AddVariable( "Elec_AbsEta","ElecAbsEta",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EtotOverPin","EtotOverPin",            "     "     , 'F'  );
    factory->AddVariable( "Elec_LateBrem","LateBrem",            "     "     , 'F'  );
    factory->AddVariable( "Elec_NumHits","NumHits",            "     "     , 'I'  );

    factory->AddVariable( "Tau_GammaEnFrac","GammaEnFrac",            "     "     , 'F'  );
    factory->AddVariable( "Tau_AbsEta","PFTauEta",            "     "     , 'F'  );
    factory->AddVariable( "Tau_EmFraction","TauEmFraction",            "     "     , 'F'  );
    factory->AddVariable( "Tau_HasGsf","TauHasGsf",            "     "     , 'I'  );
    factory->AddVariable( "Tau_HadrEoP","TauHadrEoP",            "     "     , 'F'  );
    factory->AddVariable( "Tau_HadrHoP","TauHadrHoP",            "     "     , 'F'  );
    factory->AddVariable( "Tau_NumGammaCands","TauNumGammaCands",            "     "     , 'I'  );
    factory->AddVariable( "Tau_VisMass","TauVisMass",            "     "     , 'F'  );
  }
  if(Cat_.find("TauNoGammas")!=string::npos){
    factory->AddVariable( "Elec_Fbrem", "Fbrem",            "     "     , 'F'  );
    factory->AddVariable( "Elec_Chi2KF", "Chi2KF",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EarlyBrem","EarlyBrem",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EeOverPout","EeOverPout",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EgammaOverPdif","EgammaOverPdif",            "     "     , 'F'  );
    factory->AddVariable( "Elec_AbsEta","ElecAbsEta",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EtotOverPin","EtotOverPin",            "     "     , 'F'  );
    factory->AddVariable( "Elec_LateBrem","LateBrem",            "     "     , 'F'  );
    factory->AddVariable( "Elec_NumHits","NumHits",            "     "     , 'I'  );

//     factory->AddVariable( "Tau_GammaEnFrac","GammaEnFrac",            "     "     , 'F'  );
    factory->AddVariable( "Tau_AbsEta","PFTauEta",            "     "     , 'F'  );
    factory->AddVariable( "Tau_EmFraction","TauEmFraction",            "     "     , 'F'  );
    factory->AddVariable( "Tau_HasGsf","TauHasGsf",            "     "     , 'I'  );
    factory->AddVariable( "Tau_HadrEoP","TauHadrEoP",            "     "     , 'F'  );
    factory->AddVariable( "Tau_HadrHoP","TauHadrHoP",            "     "     , 'F'  );
//     factory->AddVariable( "Tau_NumGammaCands","TauNumGammaCands",            "     "     , 'I'  );
    factory->AddVariable( "Tau_VisMass","TauVisMass",            "     "     , 'F'  );
  }
  if(Cat_.find("TauHasGammasNoGsfTrack")!=string::npos){
    factory->AddVariable( "Elec_Fbrem", "Fbrem",            "     "     , 'F'  );
    factory->AddVariable( "Elec_Chi2KF", "Chi2KF",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EarlyBrem","EarlyBrem",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EeOverPout","EeOverPout",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EgammaOverPdif","EgammaOverPdif",            "     "     , 'F'  );
    factory->AddVariable( "Elec_AbsEta","ElecAbsEta",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EtotOverPin","EtotOverPin",            "     "     , 'F'  );
    factory->AddVariable( "Elec_LateBrem","LateBrem",            "     "     , 'F'  );
    factory->AddVariable( "Elec_NumHits","NumHits",            "     "     , 'I'  );

    factory->AddVariable( "Tau_GammaEnFrac","GammaEnFrac",            "     "     , 'F'  );
    factory->AddVariable( "Tau_AbsEta","PFTauEta",            "     "     , 'F'  );
    factory->AddVariable( "Tau_EmFraction","TauEmFraction",            "     "     , 'F'  );
    factory->AddVariable( "Tau_HasGsf","TauHasGsf",            "     "     , 'I'  );
    factory->AddVariable( "Tau_HadrEoP","TauHadrEoP",            "     "     , 'F'  );
    factory->AddVariable( "Tau_HadrHoP","TauHadrHoP",            "     "     , 'F'  );
    factory->AddVariable( "Tau_NumGammaCands","TauNumGammaCands",            "     "     , 'I'  );
    factory->AddVariable( "Tau_VisMass","TauVisMass",            "     "     , 'F'  );
  }
  if(Cat_.find("TauHasGammasHasGsfTrackPFmvaBelow01")!=string::npos){
    factory->AddVariable( "Elec_Fbrem", "Fbrem",            "     "     , 'F'  );
    factory->AddVariable( "Elec_Chi2KF", "Chi2KF",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EarlyBrem","EarlyBrem",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EeOverPout","EeOverPout",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EgammaOverPdif","EgammaOverPdif",            "     "     , 'F'  );
    factory->AddVariable( "Elec_AbsEta","ElecAbsEta",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EtotOverPin","EtotOverPin",            "     "     , 'F'  );
    factory->AddVariable( "Elec_LateBrem","LateBrem",            "     "     , 'F'  );
    factory->AddVariable( "Elec_NumHits","NumHits",            "     "     , 'I'  );

    factory->AddVariable( "Tau_GammaEnFrac","GammaEnFrac",            "     "     , 'F'  );
    factory->AddVariable( "Tau_AbsEta","PFTauEta",            "     "     , 'F'  );
    factory->AddVariable( "Tau_EmFraction","TauEmFraction",            "     "     , 'F'  );
    factory->AddVariable( "Tau_HasGsf","TauHasGsf",            "     "     , 'I'  );
    factory->AddVariable( "Tau_HadrEoP","TauHadrEoP",            "     "     , 'F'  );
    factory->AddVariable( "Tau_HadrHoP","TauHadrHoP",            "     "     , 'F'  );
    factory->AddVariable( "Tau_NumGammaCands","TauNumGammaCands",            "     "     , 'I'  );
    factory->AddVariable( "Tau_VisMass","TauVisMass",            "     "     , 'F'  );
  }
  if(Cat_.find("TauHasGammasHasGsfTrackPFmvaOver01")!=string::npos){
    factory->AddVariable( "Elec_Fbrem", "Fbrem",            "     "     , 'F'  );
    factory->AddVariable( "Elec_Chi2KF", "Chi2KF",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EarlyBrem","EarlyBrem",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EeOverPout","EeOverPout",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EgammaOverPdif","EgammaOverPdif",            "     "     , 'F'  );
    factory->AddVariable( "Elec_AbsEta","ElecAbsEta",            "     "     , 'F'  );
    factory->AddVariable( "Elec_EtotOverPin","EtotOverPin",            "     "     , 'F'  );
    factory->AddVariable( "Elec_LateBrem","LateBrem",            "     "     , 'F'  );
    factory->AddVariable( "Elec_NumHits","NumHits",            "     "     , 'I'  );

    factory->AddVariable( "Tau_GammaEnFrac","GammaEnFrac",            "     "     , 'F'  );
    factory->AddVariable( "Tau_AbsEta","PFTauEta",            "     "     , 'F'  );
    factory->AddVariable( "Tau_EmFraction","TauEmFraction",            "     "     , 'F'  );
    factory->AddVariable( "Tau_HasGsf","TauHasGsf",            "     "     , 'I'  );
    factory->AddVariable( "Tau_HadrEoP","TauHadrEoP",            "     "     , 'F'  );
    factory->AddVariable( "Tau_HadrHoP","TauHadrHoP",            "     "     , 'F'  );
    factory->AddVariable( "Tau_NumGammaCands","TauNumGammaCands",            "     "     , 'I'  );
    factory->AddVariable( "Tau_VisMass","TauVisMass",            "     "     , 'F'  );
  }


  TFile *fTau = new TFile(Form("../root/tree_AntiEMVA_Ivo_%s_Tau.root",Cat_.data()),"READ"); 
  TFile *fEle = new TFile(Form("../root/tree_AntiEMVA_Ivo_%s_Elec.root",Cat_.data()),"READ");


  TTree *tTau = (TTree*)fTau->Get("tree");
  TTree *tEle = (TTree*)fEle->Get("tree");
 
  TCut myCut = "";

  if( Sel_.find("Barrel") !=string::npos)  myCut = myCut && TCut("Elec_AbsEta<1.479 && Tau_AbsEta<1.479");
  if(Sel_.find("Endcap") !=string::npos)  myCut = myCut && TCut("Elec_AbsEta>1.479 && Tau_AbsEta>1.479 && Elec_AbsEta<3.0 && Tau_AbsEta<3.0");



  ////////////////// compute the weights

  factory->AddSignalTree( tTau );
  factory->AddBackgroundTree( tEle );
//   factory->SetWeightExpression("puWeight");

  factory->PrepareTrainingAndTestTree( myCut,myCut,
				       "nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
  
  
  factory->BookMethod( TMVA::Types::kBDT, "BDT", 
		       "!H:!V:NTrees=400:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
  

  factory->TrainAllMethods();
  
  factory->TestAllMethods();

  factory->EvaluateAllMethods();  

  outputFile->Close();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;      
  
  delete factory;

}

void TMVAAllClassification(){
  TMVAClassification("All","Barrel");
  TMVAClassification("TauNoGammas","Barrel");
  TMVAClassification("TauHasGammasNoGsfTrack","Barrel");
  TMVAClassification("TauHasGammasHasGsfTrackPFmvaBelow01","Barrel");
  TMVAClassification("TauHasGammasHasGsfTrackPFmvaOver01","Barrel");

  TMVAClassification("All","Endcap");
  TMVAClassification("TauNoGammas","Endcap");
  TMVAClassification("TauHasGammasNoGsfTrack","Endcap");
  TMVAClassification("TauHasGammasHasGsfTrackPFmvaBelow01","Endcap");
  TMVAClassification("TauHasGammasHasGsfTrackPFmvaOver01","Endcap");
}

void TMVAClassificationApplication(double effS_ = 0.80 ) 
{   

  TMVA::Tools::Instance();

  float mva, visMass, HoP, dEtaG, dPhiG, dPtG, dPhi, dEta, sihih, emFraction;
  float signalPFGammaCands, hasGsf;

  TMVA::Reader *reader01 = new TMVA::Reader( "!Color:!Silent" );   
  reader01->AddVariable("TMath::Max(leadPFChargedHadrMva,-1.0)",&mva);
  reader01->AddVariable("visMass",&visMass);
  reader01->AddVariable("abs(gammadEta[0])",&dEtaG);
  reader01->AddVariable("TMath::Min(abs(gammadPhi[0]),0.3)",&dPhiG);
  reader01->AddVariable("gammaPt[0]/pt",&dPtG);
  reader01->BookMVA( "BDT", "weights/TMVAClassification_0_1_BDT.weights.xml" ); 

  TMVA::Reader *reader11 = new TMVA::Reader( "!Color:!Silent" );   
  reader11->AddVariable("TMath::Max(leadPFChargedHadrMva,-1.0)",&mva);
  reader11->AddVariable("signalPFGammaCands",&signalPFGammaCands);
  reader11->AddVariable("visMass",&visMass);
  reader11->AddVariable("dPhi",&dPhi);
  reader11->AddVariable("dEta",&dEta);
  reader11->AddVariable("sihih",&sihih);
  reader11->BookMVA( "BDT", "weights/TMVAClassification_1_1_BDT.weights.xml" );

  TMVA::Reader *readerX0 = new TMVA::Reader( "!Color:!Silent" );   
  readerX0->AddVariable("leadPFChargedHadrHcalEnergy/leadPFChargedHadrTrackP",&HoP);
  //readerX0->AddVariable("TMath::Max(leadPFChargedHadrMva,-1.0)",&mva);
  readerX0->AddVariable("TMath::Max(emFraction, 0.0)",&emFraction);
  //readerX0->AddVariable("hasGsf",&hasGsf);

  readerX0->BookMVA( "BDT", "weights/TMVAClassification_X_0_BDT.weights.xml" );

  TFile f11("TMVA1_1.root");
  TH1F* hEff11 = (TH1F*)f11.Get("Method_BDT/BDT/MVA_BDT_effS");
  float lik11=0;
  for(int lik=0; lik<=1000; lik++){
    if( fabs(hEff11->GetBinContent(hEff11->FindBin(lik*0.002)) - effS_)<0.0025 ) lik11=lik*0.002;
  }
  cout << lik11 << endl;
  TFile f01("TMVA0_1.root");
  TH1F* hEff01 = (TH1F*)f01.Get("Method_BDT/BDT/MVA_BDT_effS");
  float lik01=0;
  for(int lik=0; lik<=1000; lik++){
    if( fabs(hEff01->GetBinContent(hEff01->FindBin(lik*0.001)) - effS_)<0.0025 ) lik01=lik*0.001;
  }
  cout << lik01 << endl;
  TFile fX0("TMVAX_0.root");
  TH1F* hEffX0 = (TH1F*)fX0.Get("Method_BDT/BDT/MVA_BDT_effS");
  float likX0=0;
  for(int lik=0; lik<=1000; lik++){
    if( fabs(hEffX0->GetBinContent(hEffX0->FindBin(lik*0.001)) - effS_)<0.0025 ) likX0=lik*0.001;
  }
  cout << likX0 << endl;



  TFile *fTau = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/Utilities/tauNTuplizer_Tau.root","READ"); 
  TFile *fEle = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/Utilities/tauNTuplizer_Electron.root","READ");

  TString tree = "tauNTuplizer/tree";

  TTree *tTau = (TTree*)fTau->Get(tree);
  TTree *tEle = (TTree*)fEle->Get(tree);
 
  TCut myCutTau = "tightestHPSDBWP>0 && charge==1 && pt>20 && abs(eta)<2.3 && tightestAntiMWP>0 && tightestAntiEWP>1 && electronVeto==0 && signalPFChargedHadrCands==1";

  TFile* dummy = new TFile("dummy.root","RECREATE");  
  TTree* tTauSelected = (TTree*)tTau->CopyTree(myCutTau);
  TTree* tEleSelected = (TTree*)tEle->CopyTree(myCutTau);

  TH1F* h = new TH1F("h","",1,-10,10);
  tTauSelected->Draw("eta>>h","puWeight");
  float Tall = h->Integral();
  h->Reset();
  tTauSelected->Draw("eta>>h","puWeight*(tightestAntiEWP>2)");
  float Tpass = h->Integral();
  h->Reset();
  tEleSelected->Draw("eta>>h","puWeight");
  float Eall = h->Integral();
  h->Reset();
  tEleSelected->Draw("eta>>h","puWeight*(tightestAntiEWP>2)");
  float Epass = h->Integral();
  h->Reset();

  cout << "Sgn Eff of antiETight wrt antiEMedium = " << Tpass/Tall << endl;
  cout << "Cat Eff of antiETight wrt antiEMedium = " << Epass/Eall << endl;
  
  float leadPFChargedHadrMva_, visMass_, leadPFChargedHadrHcalEnergy_, leadPFChargedHadrTrackP_;
  float dPhi_,dEta_,sihih_,pt_, hasGsf_, emFraction_;
  float puWeight_;
  int signalPFGammaCands_;
  vector<float>* gammaPt_   = new vector<float>();
  vector<float>* gammadPhi_ = new vector<float>();
  vector<float>* gammadEta_ = new vector<float>();

  tTauSelected->SetBranchAddress( "leadPFChargedHadrMva", &leadPFChargedHadrMva_ );
  tTauSelected->SetBranchAddress( "visMass", &visMass_ );
  tTauSelected->SetBranchAddress( "leadPFChargedHadrHcalEnergy",&leadPFChargedHadrHcalEnergy_ );
  tTauSelected->SetBranchAddress( "leadPFChargedHadrTrackP",    &leadPFChargedHadrTrackP_ );
  tTauSelected->SetBranchAddress( "dPhi", &dPhi_ );
  tTauSelected->SetBranchAddress( "dEta", &dEta_ );
  tTauSelected->SetBranchAddress( "sihih",&sihih_ );
  tTauSelected->SetBranchAddress( "signalPFGammaCands", &signalPFGammaCands_ );
  tTauSelected->SetBranchAddress( "gammaPt",   &gammaPt_ );
  tTauSelected->SetBranchAddress( "gammadPhi", &gammadPhi_ );
  tTauSelected->SetBranchAddress( "gammadEta", &gammadEta_ );
  tTauSelected->SetBranchAddress( "puWeight",  &puWeight_ );
  tTauSelected->SetBranchAddress( "pt", &pt_ );
  tTauSelected->SetBranchAddress( "hasGsf", &hasGsf_ );
  tTauSelected->SetBranchAddress( "emFraction", &emFraction_ );
  
  float totalT = 0;
  float passT  = 0;

  float frac11T = 0;
  float frac01T = 0;
  float fracX0T = 0;

  for (Long64_t ievt=0; ievt<tTauSelected->GetEntries();ievt++) {
   
    if(ievt%20000==0)  cout << "Process event " << ievt << endl;

    tTauSelected->GetEntry(ievt);

    mva = TMath::Max(leadPFChargedHadrMva_,float(-1.0));
    visMass = visMass_;
    HoP = leadPFChargedHadrHcalEnergy_/leadPFChargedHadrTrackP_;
    dEtaG = fabs((*gammadEta_)[0]);
    dPhiG = TMath::Min(fabs((*gammadPhi_)[0]),0.1);
    dPtG  = (*gammaPt_)[0]/pt_;
    dPhi   = dPhi_;
    dEta   = dEta_;
    sihih  = sihih_;
    hasGsf = hasGsf_;
    emFraction =emFraction_;
    signalPFGammaCands = signalPFGammaCands_;

    if( signalPFGammaCands>0 && dEta>-98 ){
      if( reader11->EvaluateMVA( "BDT") > lik11) passT += puWeight_;
      frac11T += puWeight_;
    }
    else if( signalPFGammaCands>0 && dEta<-98 ){
      if( reader01->EvaluateMVA( "BDT") > lik01) passT += puWeight_;
      frac01T += puWeight_;
    }
    else if( signalPFGammaCands==0 ){
      if( readerX0->EvaluateMVA( "BDT") > likX0) passT += puWeight_;
      fracX0T += puWeight_;
    }
    totalT += puWeight_;

  }
  
  cout << "Signal Efficiency = " << passT/totalT << endl;
  cout << ">0 gammas and GSF="  << frac11T/(frac11T+frac01T+fracX0T) << endl;
  cout << ">0 gammas and !GSF=" << frac01T/(frac11T+frac01T+fracX0T) << endl;
  cout << "=0 gammas =" << fracX0T/(frac11T+frac01T+fracX0T) << endl;
  /////////////////////////////////////////////////////////////

  tEleSelected->SetBranchAddress( "leadPFChargedHadrMva", &leadPFChargedHadrMva_ );
  tEleSelected->SetBranchAddress( "visMass", &visMass_ );
  tEleSelected->SetBranchAddress( "leadPFChargedHadrHcalEnergy",&leadPFChargedHadrHcalEnergy_ );
  tEleSelected->SetBranchAddress( "leadPFChargedHadrTrackP",    &leadPFChargedHadrTrackP_ );
  tEleSelected->SetBranchAddress( "dPhi", &dPhi_ );
  tEleSelected->SetBranchAddress( "dEta", &dEta_ );
  tEleSelected->SetBranchAddress( "sihih",&sihih_ );
  tEleSelected->SetBranchAddress( "signalPFGammaCands", &signalPFGammaCands_ );
  tEleSelected->SetBranchAddress( "gammaPt",   &gammaPt_ );
  tEleSelected->SetBranchAddress( "gammadPhi", &gammadPhi_ );
  tEleSelected->SetBranchAddress( "gammadEta", &gammadEta_ );
  tEleSelected->SetBranchAddress( "puWeight",  &puWeight_ );
  tEleSelected->SetBranchAddress( "pt", &pt_ );
  tEleSelected->SetBranchAddress( "hasGsf", &hasGsf_ );
  tEleSelected->SetBranchAddress( "emFraction", &emFraction_ );
  
  float totalE = 0;
  float passE  = 0;

  float frac11E = 0;
  float frac01E = 0;
  float fracX0E = 0;

  for (Long64_t ievt=0; ievt<tEleSelected->GetEntries();ievt++) {
   
    if(ievt%20000==0)  cout << "Process event " << ievt << endl;

    tEleSelected->GetEntry(ievt);

    mva = TMath::Max(leadPFChargedHadrMva_,float(-1.0));
    visMass = visMass_;
    HoP = leadPFChargedHadrHcalEnergy_/leadPFChargedHadrTrackP_;
    dEtaG = fabs((*gammadEta_)[0]);
    dPhiG = TMath::Min(fabs((*gammadPhi_)[0]),0.1);
    dPtG  = (*gammaPt_)[0]/pt_;
    dPhi   = dPhi_;
    dEta   = dEta_;
    sihih  = sihih_;
    hasGsf = hasGsf_;
    emFraction =emFraction_;
    signalPFGammaCands = signalPFGammaCands_;

    if( signalPFGammaCands>0 && dEta>-98 ){
      if( reader11->EvaluateMVA( "BDT") > lik11) passE += puWeight_;
      frac11E += puWeight_;
    }
    else if( signalPFGammaCands>0 && dEta<-98 ){
      if( reader01->EvaluateMVA( "BDT") > lik01) passE += puWeight_;
      frac01E += puWeight_;
    }
    else if( signalPFGammaCands==0 ){
      if( readerX0->EvaluateMVA( "BDT") > likX0) passE += puWeight_;
      fracX0E += puWeight_;
    }
    totalE += puWeight_;
  }

  cout << "Cat Efficiency = " << passE/totalE << endl;
  cout << ">0 gammas and GSF="  << frac11E/(frac11E+frac01E+fracX0E) << endl;
  cout << ">0 gammas and !GSF=" << frac01E/(frac11E+frac01E+fracX0E) << endl;
  cout << "=0 gammas =" << fracX0E/(frac11E+frac01E+fracX0E) << endl;

}


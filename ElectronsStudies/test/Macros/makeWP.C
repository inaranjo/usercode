//--------------------------------------------------------------------------------------------------
// makeRootFilesAntiEMVA 
//
// Macro creating a rootfile with the branches corresponding to the different WP.
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

void makeWP(string data = "DYJetsToLL")
{
//    std::string inputFileName = "/data_CMS/cms/ivo/AntiEMVA/Trees/testAntiEMVA_ZZTo2e2tau_7TeV-powheg-pythia6-iter1.root";
   std::string inputFileName = "/data_CMS/cms/ivo/AntiEMVA/Trees/Trees_ForV4/AntiEMVA_AntiEMVATrees-DYJetsToLL-madgraph-PUS6.root";
  TFile* inputFile = new TFile (inputFileName.data(),"READ");
  if(inputFile->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }

  std::string outputFileName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/root/tree_makeWP_%s.root",data.data());

  TFile* outputFile = new TFile (outputFileName.data(),"RECREATE");
  TTree* mytree = new TTree("tree", "tree");


  ULong64_t t_run_,t_event_,t_lumi_;
  int t_NumPV_;
  int t_NumGsfEle_;
  int t_NumPFTaus_;
  int t_NumGenEle_;
  int t_NumGenHad_;
  int t_NumGenJet_;

  int t_Tau_GsfEleMatch_;
  int t_Tau_GenEleMatch_;

  int t_Tau_GenHadMatch_;
  float t_Tau_AbsEta_;
  float t_Tau_AbsEtaAtEcalEntrance_;
  float t_Tau_Pt_;
  float t_Tau_LeadHadronPt_;
  float t_Tau_Phi_;

  float t_Tau_AntiELoose_; 
  float t_Tau_AntiEMedium_; 
  float t_Tau_AntiETight_; 
  float t_Tau_AntiEMVA_; 
  float t_Tau_AntiEMVA2_; 
  float t_Tau_AntiEMVA2_WP75_; 
  float t_Tau_AntiEMVA2_WP85_; 
  float t_Tau_AntiEMVA2_WP95_; 

  float t_NoEleMatch_;
  float t_woG_;
  float t_wGwoGSF_;
  float t_wGwGSFwoPFMVA_;
  float t_wGwGSFwPFMVA_;

  //counters
  mytree->Branch("run",&t_run_,"run/l");
  mytree->Branch("event",&t_event_,"event/l");
  mytree->Branch("lumi",&t_lumi_,"lumi/l");
  mytree->Branch("NumPV",&t_NumPV_,"NumPV/I");
  mytree->Branch("NumGsfEle",&t_NumGsfEle_,"NumGsfEle/I");
  mytree->Branch("NumGenEle",&t_NumGenEle_,"NumGenEle/I");
  mytree->Branch("NumPFTaus",&t_NumPFTaus_,"NumPFTaus/I");
  mytree->Branch("NumGenHad",&t_NumGenHad_,"NumGenHad/I");
  mytree->Branch("NumGenJet",&t_NumGenJet_,"NumGenEJet/I");


  //PFTaus variables
  mytree->Branch("Tau_GsfEleMatch",&t_Tau_GsfEleMatch_,"Tau_GsfEleMatch/I");
  mytree->Branch("Tau_GenEleMatch",&t_Tau_GenEleMatch_,"Tau_GenEleMatch/I");
  mytree->Branch("Tau_GenHadMatch",&t_Tau_GenHadMatch_,"Tau_GenHadMatch/I");
  mytree->Branch("Tau_AbsEta",&t_Tau_AbsEta_,"Tau_AbsEta/F");
  mytree->Branch("Tau_AbsEta",&t_Tau_AbsEtaAtEcalEntrance_,"Tau_AbsEtaAtEcalEntrance/F");
  mytree->Branch("Tau_Pt",&t_Tau_Pt_,"Tau_Pt/F");
  mytree->Branch("Tau_LeadHadronPt",&t_Tau_LeadHadronPt_,"Tau_LeadHadronPt/F");
  mytree->Branch("Tau_Phi",&t_Tau_Phi_,"Tau_Phi/F");
  mytree->Branch("Tau_AntiELoose",&t_Tau_AntiELoose_,"Tau_AntiELoose/F");
  mytree->Branch("Tau_AntiEMedium",&t_Tau_AntiEMedium_,"Tau_AntiEMedium/F");
  mytree->Branch("Tau_AntiETight",&t_Tau_AntiETight_,"Tau_AntiETight/F");
  mytree->Branch("Tau_AntiEMVA",&t_Tau_AntiEMVA_,"Tau_AntiEMVA/F");
  mytree->Branch("Tau_AntiEMVA2",&t_Tau_AntiEMVA2_,"Tau_AntiEMVA2/F");
  mytree->Branch("Tau_AntiEMVA2_WP75",&t_Tau_AntiEMVA2_WP75_,"Tau_AntiEMVA2_WP75/F");
  mytree->Branch("Tau_AntiEMVA2_WP85",&t_Tau_AntiEMVA2_WP85_,"Tau_AntiEMVA2_WP85/F");
  mytree->Branch("Tau_AntiEMVA2_WP95",&t_Tau_AntiEMVA2_WP95_,"Tau_AntiEMVA2_WP95/F");

  //Categories
  mytree->Branch("NoEleMatch",&t_NoEleMatch_,"NoEleMatch/F");
  mytree->Branch("woG",&t_woG_,"woG/F");
  mytree->Branch("wGwoGSF",&t_wGwoGSF_,"wGwoGSF/F");
  mytree->Branch("wGwGSFwoPFMVA",&t_wGwGSFwoPFMVA_,"wGwGSFwoPFMVA/F");
  mytree->Branch("wGwGSFwPFMVA",&t_wGwGSFwPFMVA_,"wGwGSFwPFMVA/F");





  TTree* inputTree = (TTree*)inputFile->Get("AntiEMVAAnalyzer2/tree");
  int nEntries = inputTree->GetEntries();

  ULong64_t run,event,lumi;
  int NumPV;
  int NumGsfEle;
  int NumPFTaus;
  int NumGenEle;
  int NumGenHad;
  int NumGenJet;

  int Tau_GsfEleMatch;
  int Tau_GenEleMatch;
  int Tau_GenEleFromZMatch;
  int Tau_GenEleFromZTauTauMatch;
  int Tau_GenHadMatch;
  int Tau_GenJetMatch;
  float Tau_AbsEta;
  float Tau_Eta;
  float Tau_EtaAtEcalEntrance;
  float Tau_Pt;
  float Tau_LeadHadronPt;
  float Tau_Phi;
  float Tau_HasGsf; 
  float Tau_EmFraction; 
  float Tau_NumChargedCands;
  float Tau_NumGammaCands; 
  float Tau_HadrHoP; 
  float Tau_HadrEoP; 
  float Tau_VisMass; 
  float Tau_GammaEtaMom;
  float Tau_GammaPhiMom;
  float Tau_GammaEnFrac;
  float Tau_HadrMva; 

  float Tau_mvaAntiEValue; 
  float Tau_AntiELoose; 
  float Tau_AntiEMedium; 
  float Tau_AntiETight; 
  float Tau_AntiEMVA; 
  float Tau_AntiEMVA2; 

  int Elec_GenEleMatch;
  int Elec_GenEleFromZMatch;
  int Elec_GenEleFromZTauTauMatch;
//   int Elec_PFTauMatch;
  int Elec_GenHadMatch;
  int Elec_GenJetMatch;
  float Elec_AbsEta;
  float Elec_Pt;
  float Elec_PFMvaOutput;
  float Elec_Ee;
  float Elec_Egamma;
  float Elec_Pin;
  float Elec_Pout;
  float Elec_EtotOverPin;
  float Elec_EeOverPout;
  float Elec_EgammaOverPdif;
  int Elec_EarlyBrem;
  int Elec_LateBrem;
  float Elec_Logsihih;
  float Elec_DeltaEta;
  float Elec_HoHplusE;
  float Elec_Fbrem;
  float Elec_Chi2KF;
  float Elec_Chi2GSF;
  float Elec_NumHits;
  float Elec_GSFTrackResol;
  float Elec_GSFTracklnPt;
  float Elec_GSFTrackEta;

  inputTree->SetBranchAddress("run", &run);
  inputTree->SetBranchAddress("event", &event);
  inputTree->SetBranchAddress("lumi", &lumi);
  inputTree->SetBranchAddress("NumPV", &NumPV );
  inputTree->SetBranchAddress("NumGsfEle", &NumGsfEle );
  inputTree->SetBranchAddress("NumPFTaus", &NumPFTaus );
  inputTree->SetBranchAddress("NumGenEle", &NumGenEle );
  inputTree->SetBranchAddress("NumGenHad", &NumGenHad );
  inputTree->SetBranchAddress("NumGenJet", &NumGenJet );

  inputTree->SetBranchAddress("Tau_GsfEleMatch", &Tau_GsfEleMatch );
  inputTree->SetBranchAddress("Tau_GenEleMatch", &Tau_GenEleMatch );
  inputTree->SetBranchAddress("Tau_GenEleFromZMatch", &Tau_GenEleFromZMatch );
  inputTree->SetBranchAddress("Tau_GenEleFromZTauTauMatch", &Tau_GenEleFromZTauTauMatch );
  inputTree->SetBranchAddress("Tau_GenHadMatch", &Tau_GenHadMatch );
  inputTree->SetBranchAddress("Tau_GenJetMatch", &Tau_GenJetMatch );
//   inputTree->SetBranchAddress("Tau_AbsEta", &Tau_AbsEta );
  inputTree->SetBranchAddress("Tau_Eta", &Tau_Eta );
  inputTree->SetBranchAddress("Tau_EtaAtEcalEntrance", &Tau_EtaAtEcalEntrance );
  inputTree->SetBranchAddress("Tau_Pt", &Tau_Pt );
  inputTree->SetBranchAddress("Tau_LeadHadronPt", &Tau_LeadHadronPt );
  inputTree->SetBranchAddress("Tau_Phi", &Tau_Phi );
  inputTree->SetBranchAddress("Tau_HasGsf", &Tau_HasGsf ); 
  inputTree->SetBranchAddress("Tau_EmFraction", &Tau_EmFraction ); 
  inputTree->SetBranchAddress("Tau_NumChargedCands", &Tau_NumChargedCands );
  inputTree->SetBranchAddress("Tau_NumGammaCands", &Tau_NumGammaCands ); 
  inputTree->SetBranchAddress("Tau_HadrHoP", &Tau_HadrHoP ); 
  inputTree->SetBranchAddress("Tau_HadrEoP", &Tau_HadrEoP ); 
  inputTree->SetBranchAddress("Tau_VisMass", &Tau_VisMass ); 
  inputTree->SetBranchAddress("Tau_GammaEtaMom", &Tau_GammaEtaMom );
  inputTree->SetBranchAddress("Tau_GammaPhiMom", &Tau_GammaPhiMom );
  inputTree->SetBranchAddress("Tau_GammaEnFrac", &Tau_GammaEnFrac );
  inputTree->SetBranchAddress("Tau_HadrMva", &Tau_HadrMva ); 

  inputTree->SetBranchAddress("Tau_mvaAntiEValue", &Tau_mvaAntiEValue ); 
  inputTree->SetBranchAddress("Tau_AntiELoose", &Tau_AntiELoose ); 
  inputTree->SetBranchAddress("Tau_AntiEMedium", &Tau_AntiEMedium ); 
  inputTree->SetBranchAddress("Tau_AntiETight", &Tau_AntiETight ); 
  inputTree->SetBranchAddress("Tau_AntiEMVA", &Tau_AntiEMVA ); 
  inputTree->SetBranchAddress("Tau_AntiEMVA2", &Tau_AntiEMVA2 ); 

  inputTree->SetBranchAddress("Elec_GenEleMatch", &Elec_GenEleMatch );
  inputTree->SetBranchAddress("Elec_GenEleFromZMatch", &Elec_GenEleFromZMatch);
  inputTree->SetBranchAddress("Elec_GenEleFromZTauTauMatch", &Elec_GenEleFromZTauTauMatch );
//   inputTree->SetBranchAddress("Elec_PFTauMatch", &Elec_PFTauMatch );
  inputTree->SetBranchAddress("Elec_GenHadMatch", &Elec_GenHadMatch );
  inputTree->SetBranchAddress("Elec_GenJetMatch", &Elec_GenJetMatch );
  inputTree->SetBranchAddress("Elec_AbsEta", &Elec_AbsEta );
  inputTree->SetBranchAddress("Elec_Pt", &Elec_Pt );
  inputTree->SetBranchAddress("Elec_PFMvaOutput", &Elec_PFMvaOutput );
  inputTree->SetBranchAddress("Elec_Ee", &Elec_Ee );
  inputTree->SetBranchAddress("Elec_Egamma", &Elec_Egamma );
  inputTree->SetBranchAddress("Elec_Pin", &Elec_Pin );
  inputTree->SetBranchAddress("Elec_Pout", &Elec_Pout );
  inputTree->SetBranchAddress("Elec_EtotOverPin", &Elec_EtotOverPin );
  inputTree->SetBranchAddress("Elec_EeOverPout", &Elec_EeOverPout );
  inputTree->SetBranchAddress("Elec_EgammaOverPdif", &Elec_EgammaOverPdif );
  inputTree->SetBranchAddress("Elec_EarlyBrem", &Elec_EarlyBrem );
  inputTree->SetBranchAddress("Elec_LateBrem", &Elec_LateBrem );
  inputTree->SetBranchAddress("Elec_Logsihih", &Elec_Logsihih );
  inputTree->SetBranchAddress("Elec_DeltaEta", &Elec_DeltaEta );
  inputTree->SetBranchAddress("Elec_HoHplusE", &Elec_HoHplusE );
  inputTree->SetBranchAddress("Elec_Fbrem", &Elec_Fbrem );
  inputTree->SetBranchAddress("Elec_Chi2KF", &Elec_Chi2KF );
  inputTree->SetBranchAddress("Elec_Chi2GSF", &Elec_Chi2GSF );
  inputTree->SetBranchAddress("Elec_NumHits", &Elec_NumHits );
  inputTree->SetBranchAddress("Elec_GSFTrackResol", &Elec_GSFTrackResol );
  inputTree->SetBranchAddress("Elec_GSFTracklnPt", &Elec_GSFTracklnPt );
  inputTree->SetBranchAddress("Elec_GSFTrackEta", &Elec_GSFTrackEta );

  inputTree->SetBranchStatus("run", 1);
  inputTree->SetBranchStatus("event", 1);
  inputTree->SetBranchStatus("lumi", 1);
  inputTree->SetBranchStatus("NumPV", 1);
  inputTree->SetBranchStatus("NumGsfEle", 1);
  inputTree->SetBranchStatus("NumPFTaus", 1);
  inputTree->SetBranchStatus("NumGenEle", 1);
  inputTree->SetBranchStatus("NumGenHad", 1);
  inputTree->SetBranchStatus("NumGenJet", 1);

  inputTree->SetBranchStatus("Tau_GsfEleMatch", 1);
  inputTree->SetBranchStatus("Tau_GenEleMatch", 1);
  inputTree->SetBranchStatus("Tau_GenEleFromZMatch", 1);
  inputTree->SetBranchStatus("Tau_GenEleFromZTauTauMatch", 1);
  inputTree->SetBranchStatus("Tau_GenHadMatch", 1);
  inputTree->SetBranchStatus("Tau_GenJetMatch", 1);
//   inputTree->SetBranchStatus("Tau_AbsEta", 1);
  inputTree->SetBranchStatus("Tau_Eta", 1);
  inputTree->SetBranchStatus("Tau_EtaAtEcalEntrance", 1);
  inputTree->SetBranchStatus("Tau_Pt", 1);
  inputTree->SetBranchStatus("Tau_LeadHadronPt", 1);
  inputTree->SetBranchStatus("Tau_Phi", 1);
  inputTree->SetBranchStatus("Tau_HasGsf", 1); 
  inputTree->SetBranchStatus("Tau_EmFraction", 1); 
  inputTree->SetBranchStatus("Tau_NumChargedCands", 1);
  inputTree->SetBranchStatus("Tau_NumGammaCands", 1); 
  inputTree->SetBranchStatus("Tau_HadrHoP", 1); 
  inputTree->SetBranchStatus("Tau_HadrEoP", 1); 
  inputTree->SetBranchStatus("Tau_VisMass", 1); 
  inputTree->SetBranchStatus("Tau_GammaEtaMom", 1);
  inputTree->SetBranchStatus("Tau_GammaPhiMom", 1);
  inputTree->SetBranchStatus("Tau_GammaEnFrac", 1);
  inputTree->SetBranchStatus("Tau_HadrMva", 1); 

  inputTree->SetBranchStatus("Tau_mvaAntiEValue", 1); 
  inputTree->SetBranchStatus("Tau_AntiELoose", 1); 
  inputTree->SetBranchStatus("Tau_AntiEMedium", 1); 
  inputTree->SetBranchStatus("Tau_AntiETight", 1); 
  inputTree->SetBranchStatus("Tau_AntiEMVA", 1); 
  inputTree->SetBranchStatus("Tau_AntiEMVA2", 1); 

  inputTree->SetBranchStatus("Elec_GenEleMatch", 1);
  inputTree->SetBranchStatus("Elec_GenEleFromZMatch", 1);
  inputTree->SetBranchStatus("Elec_GenEleFromZTauTauMatch", 1);
//   inputTree->SetBranchStatus("Elec_PFTauMatch", 1);
  inputTree->SetBranchStatus("Elec_GenHadMatch", 1);
  inputTree->SetBranchStatus("Elec_GenJetMatch", 1);
  inputTree->SetBranchStatus("Elec_AbsEta", 1);
  inputTree->SetBranchStatus("Elec_Pt", 1);
  inputTree->SetBranchStatus("Elec_PFMvaOutput", 1);
  inputTree->SetBranchStatus("Elec_Ee", 1);
  inputTree->SetBranchStatus("Elec_Egamma", 1);
  inputTree->SetBranchStatus("Elec_Pin", 1);
  inputTree->SetBranchStatus("Elec_Pout", 1);
  inputTree->SetBranchStatus("Elec_EtotOverPin", 1);
  inputTree->SetBranchStatus("Elec_EeOverPout", 1);
  inputTree->SetBranchStatus("Elec_EgammaOverPdif", 1);
  inputTree->SetBranchStatus("Elec_EarlyBrem", 1);
  inputTree->SetBranchStatus("Elec_LateBrem", 1);
  inputTree->SetBranchStatus("Elec_Logsihih", 1);
  inputTree->SetBranchStatus("Elec_DeltaEta", 1);
  inputTree->SetBranchStatus("Elec_HoHplusE", 1);
  inputTree->SetBranchStatus("Elec_Fbrem", 1);
  inputTree->SetBranchStatus("Elec_Chi2KF", 1);
  inputTree->SetBranchStatus("Elec_Chi2GSF", 1);
  inputTree->SetBranchStatus("Elec_NumHits", 1);
  inputTree->SetBranchStatus("Elec_GSFTrackResol", 1);
  inputTree->SetBranchStatus("Elec_GSFTracklnPt", 1);
  inputTree->SetBranchStatus("Elec_GSFTrackEta", 1);

  cout<< "Number of entries : "<<nEntries<<endl;

  for (int iEntry = 0; iEntry<nEntries ; iEntry++){
    if(iEntry%10000==0) cout << iEntry << endl;

    inputTree->GetEntry(iEntry);
    
    Tau_AbsEta = TMath::Abs(Tau_Eta);

    t_NoEleMatch_ = 0;
    t_woG_ = 0;
    t_wGwoGSF_ = 0;
    t_wGwGSFwoPFMVA_ = 0;
    t_wGwGSFwPFMVA_ = 0;

    if (Elec_Pt<10) continue;

    t_Tau_AntiEMVA2_WP75_ = 1;
    t_Tau_AntiEMVA2_WP85_ = 1;
    t_Tau_AntiEMVA2_WP95_ = 1;
    if (Tau_NumChargedCands == 1){
    double mvaCut75 = 999.;
    double mvaCut85 = 999.;
    double mvaCut95 = 999.;
      if(Tau_GsfEleMatch<0.5){//NoEleMatch
	t_NoEleMatch_ = 1;
	if (Tau_AbsEta<1.5){//Barrel
	 mvaCut75 = -0.126102;
	 mvaCut85 = -0.0727222;
	 mvaCut95 = -0.101727;
// 	 mvaCut75 = -0.0351292;
// 	 mvaCut85 = -0.0835313;
// 	 mvaCut95 = -0.083524;
	}
	else {//Endcap
	 mvaCut75 = -0.0177012;
	 mvaCut85 = -0.084118;
	 mvaCut95 = -0.169389 ;
// 	 mvaCut75 = -0.0356601;
// 	 mvaCut85 = -0.128549;
// 	 mvaCut95 = -0.172681;
	}
      }
      else{//match TauEle
	if (Tau_AbsEta<1.5){//Barrel
	  if(Tau_NumGammaCands == 0){//woG
	    t_woG_ = 1;
	    mvaCut75 = -0.0458154;
	    mvaCut85 = -0.072778;
	    mvaCut95 = -0.130411;
// 	    mvaCut75 = 0.00125865;
// 	    mvaCut85 = -0.0814948;
// 	    mvaCut95 = -0.124782;
	  }
	  else if(Tau_NumGammaCands >= 1 && Tau_HasGsf<0.5){//wGwoGSF
	    t_wGwoGSF_ = 1;
	    mvaCut75 = -0.137043;
	    mvaCut85 = -0.137213;
	    mvaCut95 = -0.110745 ;
// 	    mvaCut75 = -0.206314;
// 	    mvaCut85 = -0.166072;
// 	    mvaCut95 = -0.124145;
	  }
	  else if(Tau_NumGammaCands >= 1 && Tau_HasGsf>0.5 && Elec_PFMvaOutput <-0.1){//wGwGSFwoPFMVA
	    t_wGwGSFwoPFMVA_ = 1;
	    mvaCut75 = 0.0332071;
	    mvaCut85 = -0.0948499;
	    mvaCut95 = -0.11647 ;
// 	    mvaCut75 = -0.157582;
// 	    mvaCut85 = -0.141405;
// 	    mvaCut95 = -0.179366;
	  }
	  else if(Tau_NumGammaCands >= 1 && Tau_HasGsf>0.5 && Elec_PFMvaOutput >-0.1){//wGwGSFwPFMVA
	    t_wGwGSFwPFMVA_ = 1;
	    mvaCut75 = -0.0448832;
	    mvaCut85 = -0.0600284;
	    mvaCut95 = -0.134414 ;
// 	    mvaCut75 = 0.373828;
// 	    mvaCut85 = -0.064168;
// 	    mvaCut95 = -0.0979345;
	  }
	}
	else{//Endcap
	  if(Tau_NumGammaCands == 0){//woG
	    t_woG_ = 1;
	    mvaCut75 = 0.189192;
	    mvaCut85 = 0.0648186;
	    mvaCut95 = -0.140327 ;
// 	    mvaCut75 = -0.0628364;
// 	    mvaCut85 = 0.0461575;
// 	    mvaCut95 = -0.150082;
	  }
	  else if(Tau_NumGammaCands >= 1 && Tau_HasGsf<0.5){//wGwoGSF
	    t_wGwoGSF_ = 1;
	    mvaCut75 = -0.0444424;
	    mvaCut85 = -0.0804441;
	    mvaCut95 = -0.102089 ;
// 	    mvaCut75 = 0.443144;
// 	    mvaCut85 = -0.0466969;
// 	    mvaCut95 = -0.119345;
	  }
	  else if(Tau_NumGammaCands >= 1 && Tau_HasGsf>0.5 && Elec_PFMvaOutput <-0.1){//wGwGSFwoPFMVA
	    t_wGwGSFwoPFMVA_ = 1;
	    mvaCut75 = -0.0938333;
	    mvaCut85 = -0.0227585;
	    mvaCut95 = -0.14057 ;
// 	    mvaCut75 = -0.138559;
// 	    mvaCut85 = 0.00619836;
// 	    mvaCut95 = -0.126053;
	  }
	  else if(Tau_NumGammaCands >= 1 && Tau_HasGsf>0.5 && Elec_PFMvaOutput >-0.1){//wGwGSFwPFMVA
	    t_wGwGSFwPFMVA_ = 1;
	    mvaCut75 = 0.144127;
	    mvaCut85 = -0.116097;
	    mvaCut95 = -0.0975809 ;
// 	    mvaCut75 = 0.00279379;
// 	    mvaCut85 = 0.208963;
// 	    mvaCut95 = -0.0852086;
	  }
	}
      }
    
    t_Tau_AntiEMVA2_WP75_ = (Tau_AntiEMVA2 > mvaCut75);
    t_Tau_AntiEMVA2_WP85_ = (Tau_AntiEMVA2 > mvaCut85);
    t_Tau_AntiEMVA2_WP95_ = (Tau_AntiEMVA2 > mvaCut95);
    }
 

    t_run_ = run;
    t_event_ = event;
    t_lumi_ = lumi;
    t_NumPV_ = NumPV ;
    t_NumGsfEle_ = NumGsfEle ;
    t_NumPFTaus_ = NumPFTaus ;
    t_NumGenEle_ = NumGenEle ;
    t_NumGenHad_ = NumGenHad ;
    t_NumGenJet_ = NumGenJet ;

    t_Tau_GsfEleMatch_ = Tau_GsfEleMatch ;
    t_Tau_GenEleMatch_ = Tau_GenEleMatch ;
    t_Tau_GenHadMatch_ = Tau_GenHadMatch ;
    t_Tau_AbsEta_ = TMath::Abs(Tau_Eta) ;
    t_Tau_AbsEtaAtEcalEntrance_ = TMath::Abs(Tau_EtaAtEcalEntrance) ;
    t_Tau_Pt_ = Tau_Pt ;
    t_Tau_LeadHadronPt_ = Tau_LeadHadronPt ;
    t_Tau_Phi_ = Tau_Phi ; 

    t_Tau_AntiELoose_ = Tau_AntiELoose ; 
    t_Tau_AntiEMedium_ = Tau_AntiEMedium ; 
    t_Tau_AntiETight_ = Tau_AntiETight ; 
    t_Tau_AntiEMVA_ = Tau_AntiEMVA ; 
    t_Tau_AntiEMVA2_ = Tau_AntiEMVA2 ; 
    
    if(DEBUG){
      cout<<endl;
      cout<<" run : "<<t_run_<<endl;
      cout<<" event : "<<t_event_<<endl;
      cout<<" lumi : "<<t_lumi_<<endl;
      cout<<" NumPV : "<<t_NumPV_<<endl;
      cout<<" NumGsfEle : "<<t_NumGsfEle_<<endl;
      cout<<" NumPFTaus : "<<t_NumPFTaus_<<endl;
      cout<<" NumGenEle : "<<t_NumGenEle_<<endl;
      cout<<" NumGenHad : "<<t_NumGenHad_<<endl;
      cout<<" NumGenJet : "<<t_NumGenJet_<<endl;
      
      cout<<" Tau_GsfEleMatch :"<<t_Tau_GsfEleMatch_<<endl;
      cout<<" Tau_GenEleMatch :"<<t_Tau_GenEleMatch_<<endl;
      cout<<" Tau_GenHadMatch :"<<t_Tau_GenHadMatch_<<endl;
      cout<<" Tau_AbsEta :"<<t_Tau_AbsEta_<<endl;
      cout<<" Tau_Pt : "<<t_Tau_Pt_<<endl;
 
      cout<<" Tau_AntiELoose :"<<t_Tau_AntiELoose_<<endl; 
      cout<<" Tau_AntiEMedium :"<<t_Tau_AntiEMedium_<<endl; 
      cout<<" Tau_AntiETight :"<<t_Tau_AntiETight_<<endl; 
      cout<<" Tau_AntiEMVA :"<<t_Tau_AntiEMVA_<<endl; 
    
    }

    mytree->Fill();

  }
  mytree->Write();
  
  cout<<"Creating file : "<<outputFileName.data()<<endl;
  inputFile->Close();
  outputFile->Close();
  return;
}





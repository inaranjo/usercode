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

void makeRoot(string matching = "Elec",
	      string category = "TauNoGammas",
	      string discriminator = ""
	      )
{


//   std::string inputFileName = "/data_CMS/cms/ivo/AntiEMVA/Trees/AntiEMVA_DYJetsToLL.root";
//   std::string inputFileName = "../AntiEMVA.root";
//   std::string inputFileName = "../AntiEMVA_Ivo-PFTauAntiEMed.root";
  std::string inputFileName = "../AntiEMVA_Ivo";
  if(discriminator=="PFTauAntiEMed") inputFileName.append(Form("-%s",discriminator.data()));
 inputFileName.append(".root");
  TFile* inputFile = new TFile (inputFileName.data(),"READ");
  if(inputFile->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }

//   std::string outputFileName = Form("./root/tree_AntiEMVA_Ivo-PFTauAntiEMed_%s_%s.root",category.data(),matching.data());
  std::string outputFileName = Form("./root/tree_AntiEMVA_Ivo_%s_%s.root",category.data(),matching.data());
  if(discriminator=="PFTauAntiEMed") outputFileName = Form("./root/tree_AntiEMVA_Ivo-%s_%s_%s.root",discriminator.data(),category.data(),matching.data());

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
  int t_Tau_GenEleFromZMatch_;
  int t_Tau_GenEleFromZTauTauMatch_;
  int t_Tau_GenHadMatch_;
  int t_Tau_GenJetMatch_;
  float t_Tau_AbsEta_;
  float t_Tau_Pt_;
  float t_Tau_HasGsf_; 
  float t_Tau_EmFraction_; 
  float t_Tau_NumChargedCands_;
  float t_Tau_NumGammaCands_; 
  float t_Tau_HadrHoP_; 
  float t_Tau_HadrEoP_; 
  float t_Tau_VisMass_; 
  float t_Tau_GammaEtaMom_;
  float t_Tau_GammaPhiMom_;
  float t_Tau_GammaEnFrac_;
  float t_Tau_HadrMva_; 

  int t_Elec_GenEleMatch_;
  int t_Elec_GenEleFromZMatch_;
  int t_Elec_GenEleFromZTauTauMatch_;
  int t_Elec_PFTauMatch_;
  int t_Elec_GenHadMatch_;
  int t_Elec_GenJetMatch_;
  float t_Elec_AbsEta_;
  float t_Elec_Pt_;
  float t_Elec_PFMvaOutput_;
  float t_Elec_Ee_;
  float t_Elec_Egamma_;
  float t_Elec_Pin_;
  float t_Elec_Pout_;
  float t_Elec_EtotOverPin_;
  float t_Elec_EeOverPout_;
  float t_Elec_EgammaOverPdif_;
  int t_Elec_EarlyBrem_;
  int t_Elec_LateBrem_;
  float t_Elec_Logsihih_;
  float t_Elec_DeltaEta_;
  float t_Elec_HoHplusE_;
  float t_Elec_Fbrem_;
  float t_Elec_Chi2KF_;
  float t_Elec_Chi2GSF_;
  float t_Elec_NumHits_;
  float t_Elec_GSFTrackResol_;
  float t_Elec_GSFTracklnPt_;
  float t_Elec_GSFTrackEta_;

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

  //GsfElectron variables
  mytree->Branch("Elec_GenEleMatch",&t_Elec_GenEleMatch_,"Elec_GenEleMatch/I");
  mytree->Branch("Elec_GenEleFromZMatch",&t_Elec_GenEleFromZMatch_,"Elec_GenEleFromZMatch/I");
  mytree->Branch("Elec_GenEleFromZTauTauMatch",&t_Elec_GenEleFromZTauTauMatch_,"Elec_GenEleFromZTauTauMatch/I");
  mytree->Branch("Elec_PFTauMatch",&t_Elec_PFTauMatch_,"Elec_PFTauMatch/I");
  mytree->Branch("Elec_GenHadMatch",&t_Elec_GenHadMatch_,"Elec_GenHadMatch/I");
  mytree->Branch("Elec_GenJetMatch",&t_Elec_GenJetMatch_,"Elec_GenJetMatch/I");
  mytree->Branch("Elec_AbsEta",&t_Elec_AbsEta_,"Elec_AbsEta/F");
  mytree->Branch("Elec_Pt",&t_Elec_Pt_,"Elec_Pt/F");
  mytree->Branch("Elec_PFMvaOutput",&t_Elec_PFMvaOutput_,"Elec_PFMvaOutput/F");
  mytree->Branch("Elec_Ee",&t_Elec_Ee_,"Elec_Ee/F");
  mytree->Branch("Elec_Egamma",&t_Elec_Egamma_,"Elec_Egamma/F");
  mytree->Branch("Elec_Pin",&t_Elec_Pin_,"Elec_Pin/F");
  mytree->Branch("Elec_Pout",&t_Elec_Pout_,"Elec_Pout/F");
  mytree->Branch("Elec_EtotOverPin",&t_Elec_EtotOverPin_,"Elec_EtotOverPin/F");
  mytree->Branch("Elec_EeOverPout",&t_Elec_EeOverPout_,"Elec_EeOverPout/F");
  mytree->Branch("Elec_EgammaOverPdif",&t_Elec_EgammaOverPdif_,"Elec_EgammaOverPdif/F");
  mytree->Branch("Elec_EarlyBrem",&t_Elec_EarlyBrem_,"Elec_EarlyBrem/I");
  mytree->Branch("Elec_LateBrem",&t_Elec_LateBrem_,"Elec_LateBrem/I");
  mytree->Branch("Elec_Logsihih",&t_Elec_Logsihih_,"Elec_Logsihih/F");
  mytree->Branch("Elec_DeltaEta",&t_Elec_DeltaEta_,"Elec_DeltaEta/F");
  mytree->Branch("Elec_HoHplusE",&t_Elec_HoHplusE_,"Elec_HoHplusE/F");
  mytree->Branch("Elec_Fbrem",&t_Elec_Fbrem_,"Elec_Fbrem/F");
  mytree->Branch("Elec_Chi2KF",&t_Elec_Chi2KF_,"Elec_Chi2KF/F");
  mytree->Branch("Elec_Chi2GSF",&t_Elec_Chi2GSF_,"Elec_Chi2GSF/F");
  mytree->Branch("Elec_NumHits",&t_Elec_NumHits_,"Elec_NumHits/F");
  mytree->Branch("Elec_GSFTrackResol",&t_Elec_GSFTrackResol_,"Elec_GSFTrackResol/F");
  mytree->Branch("Elec_GSFTracklnPt",&t_Elec_GSFTracklnPt_,"Elec_GSFTracklnPt/F");
  mytree->Branch("Elec_GSFTrackEta",&t_Elec_GSFTrackEta_,"Elec_GSFTrackEta/F");
  //PFTaus variables
  mytree->Branch("Tau_GsfEleMatch",&t_Tau_GsfEleMatch_,"Tau_GsfEleMatch/I");
  mytree->Branch("Tau_GenEleMatch",&t_Tau_GenEleMatch_,"Tau_GenEleMatch/I");
  mytree->Branch("Tau_GenEleFromZMatch",&t_Tau_GenEleFromZMatch_,"Tau_GenEleFromZMatch/I");
  mytree->Branch("Tau_GenEleFromZTauTauMatch",&t_Tau_GenEleFromZTauTauMatch_,"Tau_GenEleFromZTauTauMatch/I");
  mytree->Branch("Tau_GenHadMatch",&t_Tau_GenHadMatch_,"Tau_GenHadMatch/I");
  mytree->Branch("Tau_GenJetMatch",&t_Tau_GenJetMatch_,"Tau_GenJetMatch/I");
  mytree->Branch("Tau_AbsEta",&t_Tau_AbsEta_,"Tau_AbsEta/F");
  mytree->Branch("Tau_Pt",&t_Tau_Pt_,"Tau_Pt/F");
  mytree->Branch("Tau_HasGsf",&t_Tau_HasGsf_,"Tau_HasGsf/F");
  mytree->Branch("Tau_EmFraction",&t_Tau_EmFraction_,"Tau_EmFraction/F");
  mytree->Branch("Tau_NumChargedCands",&t_Tau_NumChargedCands_,"Tau_NumChargedCands/F");
  mytree->Branch("Tau_NumGammaCands",&t_Tau_NumGammaCands_,"Tau_NumGammaCands/F");
  mytree->Branch("Tau_HadrHoP",&t_Tau_HadrHoP_,"Tau_HadrHoP/F");
  mytree->Branch("Tau_HadrEoP",&t_Tau_HadrEoP_,"Tau_HadrEoP/F");
  mytree->Branch("Tau_VisMass",&t_Tau_VisMass_,"Tau_VisMass/F");
  mytree->Branch("Tau_GammaEtaMom",&t_Tau_GammaEtaMom_,"Tau_GammaEtaMom/F");
  mytree->Branch("Tau_GammaPhiMom",&t_Tau_GammaPhiMom_,"Tau_GammaPhiMom/F");
  mytree->Branch("Tau_GammaEnFrac",&t_Tau_GammaEnFrac_,"Tau_GammaEnFrac/F");
  mytree->Branch("Tau_HadrMva",&t_Tau_HadrMva_,"Tau_HadrMva/F");


  TTree* inputTree = (TTree*)inputFile->Get("AntiEMVAAnalyzer/tree");
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
  float Tau_Pt;
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

  int Elec_GenEleMatch;
  int Elec_GenEleFromZMatch;
  int Elec_GenEleFromZTauTauMatch;
  int Elec_PFTauMatch;
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
  inputTree->SetBranchAddress("Tau_AbsEta", &Tau_AbsEta );
  inputTree->SetBranchAddress("Tau_Pt", &Tau_Pt );
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

  inputTree->SetBranchAddress("Elec_GenEleMatch", &Elec_GenEleMatch );
  inputTree->SetBranchAddress("Elec_GenEleFromZMatch", &Elec_GenEleFromZMatch);
  inputTree->SetBranchAddress("Elec_GenEleFromZTauTauMatch", &Elec_GenEleFromZTauTauMatch );
  inputTree->SetBranchAddress("Elec_PFTauMatch", &Elec_PFTauMatch );
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
  inputTree->SetBranchStatus("Tau_AbsEta", 1);
  inputTree->SetBranchStatus("Tau_Pt", 1);
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

  inputTree->SetBranchStatus("Elec_GenEleMatch", 1);
  inputTree->SetBranchStatus("Elec_GenEleFromZMatch", 1);
  inputTree->SetBranchStatus("Elec_GenEleFromZTauTauMatch", 1);
  inputTree->SetBranchStatus("Elec_PFTauMatch", 1);
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
    if(iEntry%200==0) cout << iEntry << endl;

    inputTree->GetEntry(iEntry);


    if(matching == "Elec" && (Tau_GenEleMatch!=1 || Elec_GenEleMatch!=1)) continue;
    if(matching == "Tau" && (Tau_GenHadMatch!=1 || Elec_GenHadMatch!=1)) continue;

    if(category == "TauNoGammas" && Tau_NumGammaCands>0) continue;
    if(category == "TauHasGammasNoGsfTrack" && Tau_NumGammaCands<1 && Tau_HasGsf!=1) continue;
    if(category == "TauHasGammasHasGsfTrackPFmvaBelow01" && Tau_NumGammaCands<1 && Tau_HasGsf==1 && Elec_PFMvaOutput<-0.1) continue;
    if(category == "TauHasGammasHasGsfTrackPFmvaOver01" && Tau_NumGammaCands<1 && Tau_HasGsf==1 && Elec_PFMvaOutput>=-0.1) continue;

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
    t_Tau_GenEleFromZMatch_ = Tau_GenEleFromZMatch ;
    t_Tau_GenEleFromZTauTauMatch_ = Tau_GenEleFromZTauTauMatch ;
    t_Tau_GenHadMatch_ = Tau_GenHadMatch ;
    t_Tau_GenJetMatch_ = Tau_GenJetMatch ;
    t_Tau_AbsEta_ = Tau_AbsEta ;
    t_Tau_Pt_ = Tau_Pt ;
    t_Tau_HasGsf_ = Tau_HasGsf ; 
    t_Tau_EmFraction_ = Tau_EmFraction ; 
    t_Tau_NumChargedCands_ = Tau_NumChargedCands ;
    t_Tau_NumGammaCands_ = Tau_NumGammaCands ; 
    t_Tau_HadrHoP_ = Tau_HadrHoP ; 
    t_Tau_HadrEoP_ = Tau_HadrEoP ; 
    t_Tau_VisMass_ = Tau_VisMass ; 
    t_Tau_GammaEtaMom_ = Tau_GammaEtaMom ;
    t_Tau_GammaPhiMom_ = Tau_GammaPhiMom ;
    t_Tau_GammaEnFrac_ = Tau_GammaEnFrac ;
    t_Tau_HadrMva_ = Tau_HadrMva ; 
    
    t_Elec_GenEleMatch_ = Elec_GenEleMatch ;
    t_Elec_GenEleFromZMatch_ = Elec_GenEleFromZMatch;
    t_Elec_GenEleFromZTauTauMatch_ = Elec_GenEleFromZTauTauMatch ;
    t_Elec_PFTauMatch_ = Elec_PFTauMatch ;
    t_Elec_GenHadMatch_ = Elec_GenHadMatch ;
    t_Elec_GenJetMatch_ = Elec_GenJetMatch ;
    t_Elec_AbsEta_ = Elec_AbsEta ;
    t_Elec_Pt_ = Elec_Pt ;
    t_Elec_PFMvaOutput_ = Elec_PFMvaOutput ;
    t_Elec_Ee_ = Elec_Ee ;
    t_Elec_Egamma_ = Elec_Egamma ;
    t_Elec_Pin_ = Elec_Pin ;
    t_Elec_Pout_ = Elec_Pout ;
    t_Elec_EtotOverPin_ = Elec_EtotOverPin ;
    t_Elec_EeOverPout_ = Elec_EeOverPout ;
    t_Elec_EgammaOverPdif_ = Elec_EgammaOverPdif ;
    t_Elec_EarlyBrem_ = Elec_EarlyBrem ;
    t_Elec_LateBrem_ = Elec_LateBrem ;
    t_Elec_Logsihih_ = Elec_Logsihih ;
    t_Elec_DeltaEta_ = Elec_DeltaEta ;
    t_Elec_HoHplusE_ = Elec_HoHplusE ;
    t_Elec_Fbrem_ = Elec_Fbrem ;
    t_Elec_Chi2KF_ = Elec_Chi2KF ;
    t_Elec_Chi2GSF_ = Elec_Chi2GSF ;
    t_Elec_NumHits_ = Elec_NumHits ;
    t_Elec_GSFTrackResol_ = Elec_GSFTrackResol ;
    t_Elec_GSFTracklnPt_ = Elec_GSFTracklnPt ;
    t_Elec_GSFTrackEta_ = Elec_GSFTrackEta ;

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
      cout<<" Tau_GenEleFromZMatch :"<<t_Tau_GenEleFromZMatch_<<endl;
      cout<<" Tau_GenEleFromZTauTauMatch :"<<t_Tau_GenEleFromZTauTauMatch_<<endl;
      cout<<" Tau_GenHadMatch :"<<t_Tau_GenHadMatch_<<endl;
      cout<<" Tau_GenJetMatch :"<<t_Tau_GenJetMatch_<<endl;
      cout<<" Tau_AbsEta :"<<t_Tau_AbsEta_<<endl;
      cout<<" Tau_Pt : "<<t_Tau_Pt_<<endl;
      cout<<" Tau_HasGsf :"<<t_Tau_HasGsf_<<endl; 
      cout<<" Tau_EmFraction :"<<t_Tau_EmFraction_<<endl; 
      cout<<" Tau_NumChargedCands :"<<t_Tau_NumChargedCands_<<endl;
      cout<<" Tau_NumGammaCands :"<<t_Tau_NumGammaCands_<<endl; 
      cout<<" Tau_HadrHoP :"<<t_Tau_HadrHoP_<<endl; 
      cout<<" Tau_HadrEoP :"<<t_Tau_HadrEoP_<<endl; 
      cout<<" Tau_VisMass :"<<t_Tau_VisMass_<<endl; 
      cout<<" Tau_GammaEtaMom :"<<t_Tau_GammaEtaMom_<<endl;
      cout<<" Tau_GammaPhiMom :"<<t_Tau_GammaPhiMom_<<endl;
      cout<<" Tau_GammaEnFrac :"<<t_Tau_GammaEnFrac_<<endl;
      cout<<" Tau_HadrMva :"<<t_Tau_HadrMva_<<endl; 
    
      cout<<" Elec_GenEleMatch :"<<t_Elec_GenEleMatch_<<endl;
      cout<<" Elec_GenEleFromZMatch :"<<t_Elec_GenEleFromZMatch_<<endl;
      cout<<" Elec_GenEleFromZTauTauMatch :"<<t_Elec_GenEleFromZTauTauMatch_<<endl;
      cout<<" Elec_PFTauMatch :"<<t_Elec_PFTauMatch_<<endl;
      cout<<" Elec_GenHadMatch :"<<t_Elec_GenHadMatch_<<endl;
      cout<<" Elec_GenJetMatch :"<<t_Elec_GenJetMatch_<<endl;
      cout<<" Elec_AbsEta :"<<t_Elec_AbsEta_<<endl;
      cout<<" Elec_Pt :"<<t_Elec_Pt_<<endl;
      cout<<" Elec_PFMvaOutput : "<<t_Elec_PFMvaOutput_<<endl;
      cout<<" Elec_Ee :"<<t_Elec_Ee_<<endl;
      cout<<" Elec_Egamma :"<<t_Elec_Egamma_<<endl;
      cout<<" Elec_Pin :"<<t_Elec_Pin_<<endl;
      cout<<" Elec_Pout : "<<t_Elec_Pout_<<endl;
      cout<<" Elec_EtotOverPin :"<<t_Elec_EtotOverPin_<<endl;
      cout<<" Elec_EeOverPout : "<<t_Elec_EeOverPout_<<endl;
      cout<<" Elec_EgammaOverPdif :"<<t_Elec_EgammaOverPdif_<<endl;
      cout<<" Elec_EarlyBrem :"<<t_Elec_EarlyBrem_ <<endl;
      cout<<" Elec_LateBrem :"<<t_Elec_LateBrem_<<endl;
      cout<<" Elec_Logsihih :"<<t_Elec_Logsihih_<<endl;
      cout<<" Elec_DeltaEta :"<<t_Elec_DeltaEta_<<endl;
      cout<<" Elec_HoHplusE :"<<t_Elec_HoHplusE_<<endl;
      cout<<" Elec_Fbrem :"<<t_Elec_Fbrem_<<endl;
      cout<<" Elec_Chi2KF :"<<t_Elec_Chi2KF_<<endl;
      cout<<" Elec_Chi2GSF :"<<t_Elec_Chi2GSF_<<endl;
      cout<<" Elec_NumHits :"<<t_Elec_NumHits_<<endl;
      cout<<" Elec_GSFTrackResol :"<<t_Elec_GSFTrackResol_<<endl;
      cout<<" Elec_GSFTracklnPt : "<<t_Elec_GSFTracklnPt_<<endl;
      cout<<" Elec_GSFTrackEta :"<<t_Elec_GSFTrackEta_<<endl;
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

  makeRoot("Elec","All","");
  makeRoot("Tau","All","");
  makeRoot("Elec","TauNoGammas","");
  makeRoot("Tau","TauNoGammas","");
  makeRoot("Elec","TauHasGammasNoGsfTrack","");
  makeRoot("Tau","TauHasGammasNoGsfTrack","");
  makeRoot("Elec","TauHasGammasHasGsfTrackPFmvaBelow01","");
  makeRoot("Tau","TauHasGammasHasGsfTrackPFmvaBelow01","");
  makeRoot("Elec","TauHasGammasHasGsfTrackPFmvaOver01","");
  makeRoot("Tau","TauHasGammasHasGsfTrackPFmvaOver01","");

  makeRoot("Elec","All","PFTauAntiEMed");
  makeRoot("Tau","All","PFTauAntiEMed");
  makeRoot("Elec","TauNoGammas","PFTauAntiEMed");
  makeRoot("Tau","TauNoGammas","PFTauAntiEMed");
  makeRoot("Elec","TauHasGammasNoGsfTrack","PFTauAntiEMed");
  makeRoot("Tau","TauHasGammasNoGsfTrack","PFTauAntiEMed");
  makeRoot("Elec","TauHasGammasHasGsfTrackPFmvaBelow01","PFTauAntiEMed");
  makeRoot("Tau","TauHasGammasHasGsfTrackPFmvaBelow01","PFTauAntiEMed");
  makeRoot("Elec","TauHasGammasHasGsfTrackPFmvaOver01","PFTauAntiEMed");
  makeRoot("Tau","TauHasGammasHasGsfTrackPFmvaOver01","PFTauAntiEMed");

}

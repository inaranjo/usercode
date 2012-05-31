//--------------------------------------------------------------------------------------------------
// AntiElectronMVAApplication
//
// Macro applying the trained AntiEMVA to the rootTrees of the AntiEMVA Analyzer. 
// It also classifies the ouput Rootfiles in 4 categories. These trees will be used as input of the cut
// optimization (tmva/tmvaAntiElectronOptimization.C)
//
// Authors: I.Naranjo
//--------------------------------------------------------------------------------------------------

#include <TFile.h>
#include <TMath.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include <vector>

//Creates rootFile with a TTree with the ouput of MVA in it.
void MVAOutput(string discriminator = "",
	       string Type = "Signal"
	       )
{
//   std::string inputFileName ="/data_CMS/cms/ivo/AntiEMVA/Trees/AntiEMVA_Fall11DYJetsToLL-iter4.root";
//   std::string inputFileName ="/data_CMS/cms/ivo/AntiEMVA/Trees/Trees_ForV4/AntiEMVA_AntiEMVATrees-DYJetsToLL-madgraph-PUS6.root";
  std::string inputFileName ="/data_CMS/cms/ivo/AntiEMVA/Trees/Trees_ForV4/AntiEMVA_V4.root";
//   std::string inputFileName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/AntiEMVA_Fall11DYJetsToLL%s-iter3.root",discriminator.data());
  TFile* inputFile = new TFile (inputFileName.data(),"READ");
  if(inputFile->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
//   std::string outputFileName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/MVAOutput/MvaOutput_AntiEMVA-DYJetsToLL_v4%s_%s.root",discriminator.data(),Type.data());
//   std::string outputFileName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/MVAOutput/MvaOutput_AntiEMVA_v4%s_%s.root",discriminator.data(),Type.data());
  std::string outputFileName = Form("/data_CMS/cms/ivo/AntiEMVA/Trees/MVAOutput/MvaOutput_AntiEMVA-v5%s_%s.root",discriminator.data(),Type.data());
  
  TFile* outputFile = new TFile (outputFileName.data(),"RECREATE");
  TTree* mytree = new TTree("tree", "tree");

  float t_NoEleMatch_Barrel;
  float t_woG_Barrel;
  float t_wGwoGSF_Barrel;
  float t_wGwGSFwoPFMVA_Barrel;
  float t_wGwGSFwPFMVA_Barrel;
  float t_NoEleMatch_Endcap;
  float t_woG_Endcap;
  float t_wGwoGSF_Endcap;
  float t_wGwGSFwoPFMVA_Endcap;
  float t_wGwGSFwPFMVA_Endcap;
  
  int t_Elec_GenEleMatch;
  int t_Elec_GenHadMatch;
  int t_Tau_GenEleMatch;
  int t_Tau_GenHadMatch;
  
  mytree->Branch("NoEleMatch_Barrel",&t_NoEleMatch_Barrel,"NoEleMatch_Barrel/F");
  mytree->Branch("woG_Barrel",&t_woG_Barrel,"woG_Barrel/F");
  mytree->Branch("wGwoGSF_Barrel",&t_wGwoGSF_Barrel,"wGwoGSF_Barrel/F");
  mytree->Branch("wGwGSFwoPFMVA_Barrel",&t_wGwGSFwoPFMVA_Barrel,"wGwGSFwoPFMVA_Barrel/F");
  if(discriminator == ""){
    mytree->Branch("wGwGSFwPFMVA_Barrel",&t_wGwGSFwPFMVA_Barrel,"wGwGSFwPFMVA_Barrel/F");
  }
  mytree->Branch("NoEleMatch_Endcap",&t_NoEleMatch_Endcap,"NoEleMatch_Endcap/F");
  mytree->Branch("woG_Endcap",&t_woG_Endcap,"woG_Endcap/F");
  mytree->Branch("wGwoGSF_Endcap",&t_wGwoGSF_Endcap,"wGwoGSF_Endcap/F");
  mytree->Branch("wGwGSFwoPFMVA_Endcap",&t_wGwGSFwoPFMVA_Endcap,"wGwGSFwoPFMVA_Endcap/F");
  if(discriminator == ""){
    mytree->Branch("wGwGSFwPFMVA_Endcap",&t_wGwGSFwPFMVA_Endcap,"wGwGSFwPFMVA_Endcap/F");
  }
  mytree->Branch("Elec_GenEleMatch",&t_Elec_GenEleMatch,"Elec_GenEleMatch/I");
  mytree->Branch("Elec_GenHadMatch",&t_Elec_GenHadMatch,"Elec_GenHadMatch/I");
  mytree->Branch("Tau_GenEleMatch",&t_Tau_GenEleMatch,"Tau_GenEleMatch/I");
  mytree->Branch("Tau_GenHadMatch",&t_Tau_GenHadMatch,"Tau_GenHadMatch/I");



  TTree* inputTree = (TTree*)inputFile->Get("AntiEMVAAnalyzer2/tree");
//   TTree* inputTree = (TTree*)inputFile->Get("AntiEMVAAnalyzer2/tree");
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
  float Tau_AntiEMedium; 

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
  float Elec_EarlyBrem;//
  float Elec_LateBrem;//
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
  inputTree->SetBranchAddress("Tau_EtaAtEcalEntrance", &Tau_EtaAtEcalEntrance );
//   inputTree->SetBranchAddress("Tau_Eta", &Tau_Eta );
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
  inputTree->SetBranchAddress("Tau_AntiEMedium", &Tau_AntiEMedium ); 

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

  inputTree->SetBranchStatus("*", 1);



//   string WeightNoEleMatch_BL = "tmva/weights/TMVAClassification_v4_NoEleMatch_Barrel_BDT.weights.xml";
//   string WeightwoG_BL = "tmva/weights/TMVAClassification_v4_woG_Barrel_BDT.weights.xml";
//   string WeightwGwoGSF_BL = "tmva/weights/TMVAClassification_v4_wGwoGSF_Barrel_BDT.weights.xml";
//   string WeightwGwGSFwoPFMVA_BL = "tmva/weights/TMVAClassification_v4_wGwGSFwoPFMVA_Barrel_BDT.weights.xml";
//   string WeightwGwGSFwPFMVA_BL = "tmva/weights/TMVAClassification_v4_wGwGSFwPFMVA_Barrel_BDT.weights.xml";
//   string WeightNoEleMatch_EC = "tmva/weights/TMVAClassification_v4_NoEleMatch_Endcap_BDT.weights.xml";
//   string WeightwoG_EC = "tmva/weights/TMVAClassification_v4_woG_Endcap_BDT.weights.xml";
//   string WeightwGwoGSF_EC = "tmva/weights/TMVAClassification_v4_wGwoGSF_Endcap_BDT.weights.xml";
//   string WeightwGwGSFwoPFMVA_EC = "tmva/weights/TMVAClassification_v4_wGwGSFwoPFMVA_Endcap_BDT.weights.xml";
//   string WeightwGwGSFwPFMVA_EC = "tmva/weights/TMVAClassification_v4_wGwGSFwPFMVA_Endcap_BDT.weights.xml";
  
//   string WeightNoEleMatch_BL = "tmva/weights/TMVAClassification_v4_EtaAtEcal_NoEleMatch_Barrel_BDT.weights.xml";
//   string WeightwoG_BL = "tmva/weights/TMVAClassification_v4_EtaAtEcal_woG_Barrel_BDT.weights.xml";
//   string WeightwGwoGSF_BL = "tmva/weights/TMVAClassification_v4_EtaAtEcal_wGwoGSF_Barrel_BDT.weights.xml";
//   string WeightwGwGSFwoPFMVA_BL = "tmva/weights/TMVAClassification_v4_EtaAtEcal_wGwGSFwoPFMVA_Barrel_BDT.weights.xml";
//   string WeightwGwGSFwPFMVA_BL = "tmva/weights/TMVAClassification_v4_EtaAtEcal_wGwGSFwPFMVA_Barrel_BDT.weights.xml";
//   string WeightNoEleMatch_EC = "tmva/weights/TMVAClassification_v4_EtaAtEcal_NoEleMatch_Endcap_BDT.weights.xml";
//   string WeightwoG_EC = "tmva/weights/TMVAClassification_v4_EtaAtEcal_woG_Endcap_BDT.weights.xml";
//   string WeightwGwoGSF_EC = "tmva/weights/TMVAClassification_v4_EtaAtEcal_wGwoGSF_Endcap_BDT.weights.xml";
//   string WeightwGwGSFwoPFMVA_EC = "tmva/weights/TMVAClassification_v4_EtaAtEcal_wGwGSFwoPFMVA_Endcap_BDT.weights.xml";
//   string WeightwGwGSFwPFMVA_EC = "tmva/weights/TMVAClassification_v4_EtaAtEcal_wGwGSFwPFMVA_Endcap_BDT.weights.xml";
 
  string WeightNoEleMatch_BL = "tmva/weights/TMVAClassification_v5_NoEleMatch_Barrel_BDTG.weights.xml";
  string WeightwoG_BL = "tmva/weights/TMVAClassification_v5_woG_Barrel_BDTG.weights.xml";
  string WeightwGwoGSF_BL = "tmva/weights/TMVAClassification_v5_wGwoGSF_Barrel_BDTG.weights.xml";
  string WeightwGwGSFwoPFMVA_BL = "tmva/weights/TMVAClassification_v5_wGwGSFwoPFMVA_Barrel_BDTG.weights.xml";
  string WeightwGwGSFwPFMVA_BL = "tmva/weights/TMVAClassification_v5_wGwGSFwPFMVA_Barrel_BDTG.weights.xml";
  string WeightNoEleMatch_EC = "tmva/weights/TMVAClassification_v5_NoEleMatch_Endcap_BDTG.weights.xml";
  string WeightwoG_EC = "tmva/weights/TMVAClassification_v5_woG_Endcap_BDTG.weights.xml";
  string WeightwGwoGSF_EC = "tmva/weights/TMVAClassification_v5_wGwoGSF_Endcap_BDTG.weights.xml";
  string WeightwGwGSFwoPFMVA_EC = "tmva/weights/TMVAClassification_v5_wGwGSFwoPFMVA_Endcap_BDTG.weights.xml";
  string WeightwGwGSFwPFMVA_EC = "tmva/weights/TMVAClassification_v5_wGwGSFwPFMVA_Endcap_BDTG.weights.xml";
   
  if (discriminator.find("-AntiEMed")!=std::string::npos){
    WeightNoEleMatch_BL = "tmva/weights/TMVAClassification_v4-AntiEMed_NoEleMatch_Barrel_BDT.weights.xml";
    WeightwoG_BL = "tmva/weights/TMVAClassification_v4-AntiEMed_woG_Barrel_BDT.weights.xml";
    WeightwGwoGSF_BL = "tmva/weights/TMVAClassification_v4-AntiEMed_wGwoGSF_Barrel_BDT.weights.xml";
    WeightwGwGSFwoPFMVA_BL = "tmva/weights/TMVAClassification_v4-AntiEMed_wGwGSFwoPFMVA_Barrel_BDT.weights.xml";
//     WeightwGwGSFwPFMVA_BL = "tmva/weights/TMVAClassification-AntiEMed_wGwGSFwPFMVA_Barrel_BDT.weights.xml";
    WeightNoEleMatch_EC = "tmva/weights/TMVAClassification_v4-AntiEMed_NoEleMatch_Endcap_BDT.weights.xml";
    WeightwoG_EC = "tmva/weights/TMVAClassification_v4-AntiEMed_woG_Endcap_BDT.weights.xml";
    WeightwGwoGSF_EC = "tmva/weights/TMVAClassification_v4-AntiEMed_wGwoGSF_Endcap_BDT.weights.xml";
    WeightwGwGSFwoPFMVA_EC = "tmva/weights/TMVAClassification_v4-AntiEMed_wGwGSFwoPFMVA_Endcap_BDT.weights.xml";
//     WeightwGwGSFwPFMVA_EC = "tmva/weights/TMVAClassification-AntiEMed_wGwGSFwPFMVA_Endcap_BDT.weights.xml";
  }
  
  TMVA::Reader *readerNoEleMatch_BL = new TMVA::Reader( "!Color:!Silent:Error" );  
//for v4  readerNoEleMatch_BL->AddVariable("Tau_AbsEta",&Tau_AbsEta);
//   readerNoEleMatch_BL->AddVariable("Tau_Eta",&Tau_Eta);
  readerNoEleMatch_BL->AddVariable("Tau_EtaAtEcalEntrance",&Tau_EtaAtEcalEntrance);
  readerNoEleMatch_BL->AddVariable("Tau_Pt",&Tau_Pt);
  // forv3readerNoEleMatch_BL->AddVariable("Tau_HasGsf",&Tau_HasGsf);
  readerNoEleMatch_BL->AddVariable("Tau_EmFraction",&Tau_EmFraction);
   // forv3readerNoEleMatch_BL->AddVariable("Tau_NumChargedCands",&Tau_NumChargedCands);
  readerNoEleMatch_BL->AddVariable("Tau_NumGammaCands",&Tau_NumGammaCands);
  readerNoEleMatch_BL->AddVariable("Tau_HadrHoP",&Tau_HadrHoP);
  readerNoEleMatch_BL->AddVariable("Tau_HadrEoP",&Tau_HadrEoP);
  readerNoEleMatch_BL->AddVariable("Tau_VisMass",&Tau_VisMass);
  readerNoEleMatch_BL->AddVariable("Tau_GammaEtaMom",&Tau_GammaEtaMom);
  readerNoEleMatch_BL->AddVariable("Tau_GammaPhiMom",&Tau_GammaPhiMom);
  readerNoEleMatch_BL->AddVariable("Tau_GammaEnFrac",&Tau_GammaEnFrac);
  readerNoEleMatch_BL->SetVerbose(kTRUE);
  readerNoEleMatch_BL->BookMVA("BDT",WeightNoEleMatch_BL);

  TMVA::Reader *readerwoG_BL = new TMVA::Reader( "!Color:!Silent" );  
  //for v2  readerwoG_BL->AddVariable("Elec_Pt",&Elec_Pt);   
  //for v2  readerwoG_BL->AddVariable("Elec_AbsEta",&Elec_AbsEta);
  readerwoG_BL->AddVariable("Elec_EtotOverPin",&Elec_EtotOverPin);
  //for v2  readerwoG_BL->AddVariable("Elec_EarlyBrem",&Elec_EarlyBrem);
  readerwoG_BL->AddVariable("Elec_LateBrem",&Elec_LateBrem);
  readerwoG_BL->AddVariable("Elec_Fbrem",&Elec_Fbrem);
  readerwoG_BL->AddVariable("Elec_Chi2KF",&Elec_Chi2KF);
   //for v2 readerwoG_BL->AddVariable("Elec_NumHits",&Elec_NumHits);
  readerwoG_BL->AddVariable("Elec_GSFTrackResol",&Elec_GSFTrackResol);
  readerwoG_BL->AddVariable("Elec_GSFTracklnPt",&Elec_GSFTracklnPt);
  readerwoG_BL->AddVariable("Elec_GSFTrackEta",&Elec_GSFTrackEta);
//for v4  readerwoG_BL->AddVariable("Tau_AbsEta",&Tau_AbsEta);
//   readerwoG_BL->AddVariable("Tau_Eta",&Tau_Eta);
  readerwoG_BL->AddVariable("Tau_EtaAtEcalEntrance",&Tau_EtaAtEcalEntrance);
  readerwoG_BL->AddVariable("Tau_Pt",&Tau_Pt);
 //for v3   readerwoG_BL->AddVariable("Tau_HasGsf",&Tau_HasGsf);
  readerwoG_BL->AddVariable("Tau_EmFraction",&Tau_EmFraction);
 //for v3 readerwoG_BL->AddVariable("Tau_NumChargedCands",&Tau_NumChargedCands);
  readerwoG_BL->AddVariable("Tau_HadrHoP",&Tau_HadrHoP);
  readerwoG_BL->AddVariable("Tau_HadrEoP",&Tau_HadrEoP);
  readerwoG_BL->AddVariable("Tau_VisMass",&Tau_VisMass);
  readerwoG_BL->SetVerbose(kTRUE);
  readerwoG_BL->BookMVA("BDT",WeightwoG_BL  );
  
  TMVA::Reader *readerwGwoGSF_BL = new TMVA::Reader( "!Color:!Silent" ); 
   //for v2 readerwGwoGSF_BL->AddVariable("Elec_Pt",&Elec_Pt);   
   //for v2 readerwGwoGSF_BL->AddVariable("Elec_AbsEta",&Elec_AbsEta);
  readerwGwoGSF_BL->AddVariable("Elec_EtotOverPin",&Elec_EtotOverPin);
  readerwGwoGSF_BL->AddVariable("Elec_EgammaOverPdif",&Elec_EgammaOverPdif);
   //for v2 readerwGwoGSF_BL->AddVariable("Elec_EarlyBrem",&Elec_EarlyBrem);
  readerwGwoGSF_BL->AddVariable("Elec_LateBrem",&Elec_LateBrem);
  readerwGwoGSF_BL->AddVariable("Elec_Fbrem",&Elec_Fbrem);
  readerwGwoGSF_BL->AddVariable("Elec_Chi2GSF",&Elec_Chi2GSF);
  readerwGwoGSF_BL->AddVariable("Elec_NumHits",&Elec_NumHits);
  readerwGwoGSF_BL->AddVariable("Elec_GSFTrackResol",&Elec_GSFTrackResol);
  readerwGwoGSF_BL->AddVariable("Elec_GSFTracklnPt",&Elec_GSFTracklnPt);
  readerwGwoGSF_BL->AddVariable("Elec_GSFTrackEta",&Elec_GSFTrackEta);
 //for v4 readerwGwoGSF_BL->AddVariable("Tau_AbsEta",&Tau_AbsEta);
//   readerwGwoGSF_BL->AddVariable("Tau_Eta",&Tau_Eta);
  readerwGwoGSF_BL->AddVariable("Tau_EtaAtEcalEntrance",&Tau_EtaAtEcalEntrance);
  readerwGwoGSF_BL->AddVariable("Tau_Pt",&Tau_Pt);
  readerwGwoGSF_BL->AddVariable("Tau_EmFraction",&Tau_EmFraction);
  readerwGwoGSF_BL->AddVariable("Tau_NumGammaCands",&Tau_NumGammaCands);
  readerwGwoGSF_BL->AddVariable("Tau_HadrHoP",&Tau_HadrHoP);
  readerwGwoGSF_BL->AddVariable("Tau_HadrEoP",&Tau_HadrEoP);
  readerwGwoGSF_BL->AddVariable("Tau_VisMass",&Tau_VisMass);
  readerwGwoGSF_BL->AddVariable("Tau_GammaEtaMom",&Tau_GammaEtaMom);
  readerwGwoGSF_BL->AddVariable("Tau_GammaPhiMom",&Tau_GammaPhiMom);
  readerwGwoGSF_BL->AddVariable("Tau_GammaEnFrac",&Tau_GammaEnFrac);
  readerwGwoGSF_BL->SetVerbose(kTRUE);
  readerwGwoGSF_BL->BookMVA("BDT",WeightwGwoGSF_BL );   
  
  TMVA::Reader *readerwGwGSFwoPFMVA_BL = new TMVA::Reader( "!Color:!Silent" ); 
   //for v2 readerwGwGSFwoPFMVA_BL->AddVariable("Elec_Pt",&Elec_Pt);   
   //for v2 readerwGwGSFwoPFMVA_BL->AddVariable("Elec_AbsEta",&Elec_AbsEta);
  readerwGwGSFwoPFMVA_BL->AddVariable("Elec_Fbrem",&Elec_Fbrem);
  readerwGwGSFwoPFMVA_BL->AddVariable("Elec_Chi2KF",&Elec_Chi2KF);
  readerwGwGSFwoPFMVA_BL->AddVariable("Elec_Chi2GSF",&Elec_Chi2GSF);
  readerwGwGSFwoPFMVA_BL->AddVariable("Elec_NumHits",&Elec_NumHits);
  readerwGwGSFwoPFMVA_BL->AddVariable("Elec_GSFTrackResol",&Elec_GSFTrackResol);
  readerwGwGSFwoPFMVA_BL->AddVariable("Elec_GSFTracklnPt",&Elec_GSFTracklnPt);
  readerwGwGSFwoPFMVA_BL->AddVariable("Elec_GSFTrackEta",&Elec_GSFTrackEta);
 //for v4  readerwGwGSFwoPFMVA_BL->AddVariable("Tau_AbsEta",&Tau_AbsEta);
//   readerwGwGSFwoPFMVA_BL->AddVariable("Tau_Eta",&Tau_Eta);
  readerwGwGSFwoPFMVA_BL->AddVariable("Tau_EtaAtEcalEntrance",&Tau_EtaAtEcalEntrance);
  readerwGwGSFwoPFMVA_BL->AddVariable("Tau_Pt",&Tau_Pt);
  readerwGwGSFwoPFMVA_BL->AddVariable("Tau_EmFraction",&Tau_EmFraction);
  readerwGwGSFwoPFMVA_BL->AddVariable("Tau_NumGammaCands",&Tau_NumGammaCands);
  readerwGwGSFwoPFMVA_BL->AddVariable("Tau_HadrHoP",&Tau_HadrHoP);
  readerwGwGSFwoPFMVA_BL->AddVariable("Tau_HadrEoP",&Tau_HadrEoP);
  readerwGwGSFwoPFMVA_BL->AddVariable("Tau_VisMass",&Tau_VisMass);
  readerwGwGSFwoPFMVA_BL->AddVariable("Tau_GammaEtaMom",&Tau_GammaEtaMom);
  readerwGwGSFwoPFMVA_BL->AddVariable("Tau_GammaPhiMom",&Tau_GammaPhiMom);
  readerwGwGSFwoPFMVA_BL->AddVariable("Tau_GammaEnFrac",&Tau_GammaEnFrac);
  readerwGwGSFwoPFMVA_BL->SetVerbose(kTRUE);
  readerwGwGSFwoPFMVA_BL->BookMVA("BDT",WeightwGwGSFwoPFMVA_BL ); 
  
  TMVA::Reader *readerwGwGSFwPFMVA_BL = new TMVA::Reader( "!Color:!Silent" ); 
   //for v2 readerwGwGSFwPFMVA_BL->AddVariable("Elec_Pt",&Elec_Pt);   
   //for v2 readerwGwGSFwPFMVA_BL->AddVariable("Elec_AbsEta",&Elec_AbsEta);
  readerwGwGSFwPFMVA_BL->AddVariable("Elec_EtotOverPin",&Elec_EtotOverPin);
  readerwGwGSFwPFMVA_BL->AddVariable("Elec_EeOverPout",&Elec_EeOverPout);
  //for v3  readerwGwGSFwPFMVA_BL->AddVariable("Elec_EarlyBrem",&Elec_EarlyBrem);
  readerwGwGSFwPFMVA_BL->AddVariable("Elec_LateBrem",&Elec_LateBrem);
  readerwGwGSFwPFMVA_BL->AddVariable("Elec_Chi2GSF",&Elec_Chi2GSF);
  readerwGwGSFwPFMVA_BL->AddVariable("Elec_NumHits",&Elec_NumHits);
  readerwGwGSFwPFMVA_BL->AddVariable("Elec_GSFTrackResol",&Elec_GSFTrackResol);
  readerwGwGSFwPFMVA_BL->AddVariable("Elec_GSFTracklnPt",&Elec_GSFTracklnPt);
  readerwGwGSFwPFMVA_BL->AddVariable("Elec_GSFTrackEta",&Elec_GSFTrackEta);
 //for v4  readerwGwGSFwPFMVA_BL->AddVariable("Tau_AbsEta",&Tau_AbsEta);
//   readerwGwGSFwPFMVA_BL->AddVariable("Tau_Eta",&Tau_Eta);
  readerwGwGSFwPFMVA_BL->AddVariable("Tau_EtaAtEcalEntrance",&Tau_EtaAtEcalEntrance);
  readerwGwGSFwPFMVA_BL->AddVariable("Tau_Pt",&Tau_Pt);
  readerwGwGSFwPFMVA_BL->AddVariable("Tau_EmFraction",&Tau_EmFraction);
  readerwGwGSFwPFMVA_BL->AddVariable("Tau_NumGammaCands",&Tau_NumGammaCands);
  readerwGwGSFwPFMVA_BL->AddVariable("Tau_HadrHoP",&Tau_HadrHoP);
  readerwGwGSFwPFMVA_BL->AddVariable("Tau_HadrEoP",&Tau_HadrEoP);
  readerwGwGSFwPFMVA_BL->AddVariable("Tau_VisMass",&Tau_VisMass);
  readerwGwGSFwPFMVA_BL->AddVariable("Tau_GammaEtaMom",&Tau_GammaEtaMom);
  readerwGwGSFwPFMVA_BL->AddVariable("Tau_GammaPhiMom",&Tau_GammaPhiMom);
  readerwGwGSFwPFMVA_BL->AddVariable("Tau_GammaEnFrac",&Tau_GammaEnFrac);
  readerwGwGSFwPFMVA_BL->SetVerbose(kTRUE);
  if(discriminator == "")readerwGwGSFwPFMVA_BL->BookMVA("BDT",WeightwGwGSFwPFMVA_BL );  
  
  ////////////////////////

  TMVA::Reader *readerNoEleMatch_EC = new TMVA::Reader( "!Color:!Silent:Error" );  
 //for v4  readerNoEleMatch_EC->AddVariable("Tau_AbsEta",&Tau_AbsEta);
//   readerNoEleMatch_EC->AddVariable("Tau_Eta",&Tau_Eta);
  readerNoEleMatch_EC->AddVariable("Tau_EtaAtEcalEntrance",&Tau_EtaAtEcalEntrance);
  readerNoEleMatch_EC->AddVariable("Tau_Pt",&Tau_Pt);
  // forv3readerNoEleMatch_EC->AddVariable("Tau_HasGsf",&Tau_HasGsf);
  readerNoEleMatch_EC->AddVariable("Tau_EmFraction",&Tau_EmFraction);
  // forv3readerNoEleMatch_EC->AddVariable("Tau_NumChargedCands",&Tau_NumChargedCands);
  readerNoEleMatch_EC->AddVariable("Tau_NumGammaCands",&Tau_NumGammaCands);
  readerNoEleMatch_EC->AddVariable("Tau_HadrHoP",&Tau_HadrHoP);
  readerNoEleMatch_EC->AddVariable("Tau_HadrEoP",&Tau_HadrEoP);
  readerNoEleMatch_EC->AddVariable("Tau_VisMass",&Tau_VisMass);
  readerNoEleMatch_EC->AddVariable("Tau_GammaEtaMom",&Tau_GammaEtaMom);
  readerNoEleMatch_EC->AddVariable("Tau_GammaPhiMom",&Tau_GammaPhiMom);
  readerNoEleMatch_EC->AddVariable("Tau_GammaEnFrac",&Tau_GammaEnFrac);
  readerNoEleMatch_EC->SetVerbose(kTRUE);
  readerNoEleMatch_EC->BookMVA("BDT",WeightNoEleMatch_EC);

  TMVA::Reader *readerwoG_EC = new TMVA::Reader( "!Color:!Silent" ); 
   //for v2 readerwoG_EC->AddVariable("Elec_Pt",&Elec_Pt);   
   //for v2 readerwoG_EC->AddVariable("Elec_AbsEta",&Elec_AbsEta);
  readerwoG_EC->AddVariable("Elec_EtotOverPin",&Elec_EtotOverPin);
   //for v2 readerwoG_EC->AddVariable("Elec_EarlyBrem",&Elec_EarlyBrem);
  readerwoG_EC->AddVariable("Elec_LateBrem",&Elec_LateBrem);
  readerwoG_EC->AddVariable("Elec_Fbrem",&Elec_Fbrem);
  readerwoG_EC->AddVariable("Elec_Chi2KF",&Elec_Chi2KF);
   //for v2 readerwoG_EC->AddVariable("Elec_NumHits",&Elec_NumHits);
  readerwoG_EC->AddVariable("Elec_GSFTrackResol",&Elec_GSFTrackResol);
  readerwoG_EC->AddVariable("Elec_GSFTracklnPt",&Elec_GSFTracklnPt);
  readerwoG_EC->AddVariable("Elec_GSFTrackEta",&Elec_GSFTrackEta);
  //for v4 readerwoG_EC->AddVariable("Tau_AbsEta",&Tau_AbsEta);
//   readerwoG_EC->AddVariable("Tau_Eta",&Tau_Eta);
  readerwoG_EC->AddVariable("Tau_EtaAtEcalEntrance",&Tau_EtaAtEcalEntrance);
  readerwoG_EC->AddVariable("Tau_Pt",&Tau_Pt);
 //for v3   readerwoG_EC->AddVariable("Tau_HasGsf",&Tau_HasGsf);
  readerwoG_EC->AddVariable("Tau_EmFraction",&Tau_EmFraction);
 //for v3 readerwoG_EC->AddVariable("Tau_NumChargedCands",&Tau_NumChargedCands);
  readerwoG_EC->AddVariable("Tau_HadrHoP",&Tau_HadrHoP);
  readerwoG_EC->AddVariable("Tau_HadrEoP",&Tau_HadrEoP);
  readerwoG_EC->AddVariable("Tau_VisMass",&Tau_VisMass);
  readerwoG_EC->SetVerbose(kTRUE);
  readerwoG_EC->BookMVA("BDT",WeightwoG_EC  );
  
  TMVA::Reader *readerwGwoGSF_EC = new TMVA::Reader( "!Color:!Silent" );
   //for v2 readerwGwoGSF_EC->AddVariable("Elec_Pt",&Elec_Pt);   
   //for v2 readerwGwoGSF_EC->AddVariable("Elec_AbsEta",&Elec_AbsEta);
  readerwGwoGSF_EC->AddVariable("Elec_EtotOverPin",&Elec_EtotOverPin);
  readerwGwoGSF_EC->AddVariable("Elec_EgammaOverPdif",&Elec_EgammaOverPdif);
   //for v2 readerwGwoGSF_EC->AddVariable("Elec_EarlyBrem",&Elec_EarlyBrem);
  readerwGwoGSF_EC->AddVariable("Elec_LateBrem",&Elec_LateBrem);
  readerwGwoGSF_EC->AddVariable("Elec_Fbrem",&Elec_Fbrem);
  readerwGwoGSF_EC->AddVariable("Elec_Chi2GSF",&Elec_Chi2GSF);
  readerwGwoGSF_EC->AddVariable("Elec_NumHits",&Elec_NumHits);
  readerwGwoGSF_EC->AddVariable("Elec_GSFTrackResol",&Elec_GSFTrackResol);
  readerwGwoGSF_EC->AddVariable("Elec_GSFTracklnPt",&Elec_GSFTracklnPt);
  readerwGwoGSF_EC->AddVariable("Elec_GSFTrackEta",&Elec_GSFTrackEta);
  //for v4  readerwGwoGSF_EC->AddVariable("Tau_AbsEta",&Tau_AbsEta);
//   readerwGwoGSF_EC->AddVariable("Tau_Eta",&Tau_Eta);
  readerwGwoGSF_EC->AddVariable("Tau_EtaAtEcalEntrance",&Tau_EtaAtEcalEntrance);
  readerwGwoGSF_EC->AddVariable("Tau_Pt",&Tau_Pt);
  readerwGwoGSF_EC->AddVariable("Tau_EmFraction",&Tau_EmFraction);
  readerwGwoGSF_EC->AddVariable("Tau_NumGammaCands",&Tau_NumGammaCands);
  readerwGwoGSF_EC->AddVariable("Tau_HadrHoP",&Tau_HadrHoP);
  readerwGwoGSF_EC->AddVariable("Tau_HadrEoP",&Tau_HadrEoP);
  readerwGwoGSF_EC->AddVariable("Tau_VisMass",&Tau_VisMass);
  readerwGwoGSF_EC->AddVariable("Tau_GammaEtaMom",&Tau_GammaEtaMom);
  readerwGwoGSF_EC->AddVariable("Tau_GammaPhiMom",&Tau_GammaPhiMom);
  readerwGwoGSF_EC->AddVariable("Tau_GammaEnFrac",&Tau_GammaEnFrac);
  readerwGwoGSF_EC->SetVerbose(kTRUE);
  readerwGwoGSF_EC->BookMVA("BDT",WeightwGwoGSF_EC ); 

  TMVA::Reader *readerwGwGSFwoPFMVA_EC = new TMVA::Reader( "!Color:!Silent" );   
   //for v2 readerwGwGSFwoPFMVA_EC->AddVariable("Elec_Pt",&Elec_Pt);   
   //for v2 readerwGwGSFwoPFMVA_EC->AddVariable("Elec_AbsEta",&Elec_AbsEta);
  readerwGwGSFwoPFMVA_EC->AddVariable("Elec_Fbrem",&Elec_Fbrem);
  readerwGwGSFwoPFMVA_EC->AddVariable("Elec_Chi2KF",&Elec_Chi2KF);
  readerwGwGSFwoPFMVA_EC->AddVariable("Elec_Chi2GSF",&Elec_Chi2GSF);
  readerwGwGSFwoPFMVA_EC->AddVariable("Elec_NumHits",&Elec_NumHits);
  readerwGwGSFwoPFMVA_EC->AddVariable("Elec_GSFTrackResol",&Elec_GSFTrackResol);
  readerwGwGSFwoPFMVA_EC->AddVariable("Elec_GSFTracklnPt",&Elec_GSFTracklnPt);
  readerwGwGSFwoPFMVA_EC->AddVariable("Elec_GSFTrackEta",&Elec_GSFTrackEta);
  //for v4   readerwGwGSFwoPFMVA_EC->AddVariable("Tau_AbsEta",&Tau_AbsEta);
//   readerwGwGSFwoPFMVA_EC->AddVariable("Tau_Eta",&Tau_Eta);
  readerwGwGSFwoPFMVA_EC->AddVariable("Tau_EtaAtEcalEntrance",&Tau_EtaAtEcalEntrance);
  readerwGwGSFwoPFMVA_EC->AddVariable("Tau_Pt",&Tau_Pt);
  readerwGwGSFwoPFMVA_EC->AddVariable("Tau_EmFraction",&Tau_EmFraction);
  readerwGwGSFwoPFMVA_EC->AddVariable("Tau_NumGammaCands",&Tau_NumGammaCands);
  readerwGwGSFwoPFMVA_EC->AddVariable("Tau_HadrHoP",&Tau_HadrHoP);
  readerwGwGSFwoPFMVA_EC->AddVariable("Tau_HadrEoP",&Tau_HadrEoP);
  readerwGwGSFwoPFMVA_EC->AddVariable("Tau_VisMass",&Tau_VisMass);
  readerwGwGSFwoPFMVA_EC->AddVariable("Tau_GammaEtaMom",&Tau_GammaEtaMom);
  readerwGwGSFwoPFMVA_EC->AddVariable("Tau_GammaPhiMom",&Tau_GammaPhiMom);
  readerwGwGSFwoPFMVA_EC->AddVariable("Tau_GammaEnFrac",&Tau_GammaEnFrac);
  readerwGwGSFwoPFMVA_EC->SetVerbose(kTRUE);
  readerwGwGSFwoPFMVA_EC->BookMVA("BDT",WeightwGwGSFwoPFMVA_EC ); 

  TMVA::Reader *readerwGwGSFwPFMVA_EC = new TMVA::Reader( "!Color:!Silent" ); 
   //for v2 readerwGwGSFwPFMVA_EC->AddVariable("Elec_Pt",&Elec_Pt);   
   //for v2 readerwGwGSFwPFMVA_EC->AddVariable("Elec_AbsEta",&Elec_AbsEta);
  readerwGwGSFwPFMVA_EC->AddVariable("Elec_EtotOverPin",&Elec_EtotOverPin);
  readerwGwGSFwPFMVA_EC->AddVariable("Elec_EeOverPout",&Elec_EeOverPout);
//for v3    readerwGwGSFwPFMVA_EC->AddVariable("Elec_EarlyBrem",&Elec_EarlyBrem);
  readerwGwGSFwPFMVA_EC->AddVariable("Elec_LateBrem",&Elec_LateBrem);
  readerwGwGSFwPFMVA_EC->AddVariable("Elec_Chi2GSF",&Elec_Chi2GSF);
  readerwGwGSFwPFMVA_EC->AddVariable("Elec_NumHits",&Elec_NumHits);
  readerwGwGSFwPFMVA_EC->AddVariable("Elec_GSFTrackResol",&Elec_GSFTrackResol);
  readerwGwGSFwPFMVA_EC->AddVariable("Elec_GSFTracklnPt",&Elec_GSFTracklnPt);
  readerwGwGSFwPFMVA_EC->AddVariable("Elec_GSFTrackEta",&Elec_GSFTrackEta);
   //for v4  readerwGwGSFwPFMVA_EC->AddVariable("Tau_AbsEta",&Tau_AbsEta);
//   readerwGwGSFwPFMVA_EC->AddVariable("Tau_Eta",&Tau_Eta);
  readerwGwGSFwPFMVA_EC->AddVariable("Tau_EtaAtEcalEntrance",&Tau_EtaAtEcalEntrance);
  readerwGwGSFwPFMVA_EC->AddVariable("Tau_Pt",&Tau_Pt);
  readerwGwGSFwPFMVA_EC->AddVariable("Tau_EmFraction",&Tau_EmFraction);
  readerwGwGSFwPFMVA_EC->AddVariable("Tau_NumGammaCands",&Tau_NumGammaCands);
  readerwGwGSFwPFMVA_EC->AddVariable("Tau_HadrHoP",&Tau_HadrHoP);
  readerwGwGSFwPFMVA_EC->AddVariable("Tau_HadrEoP",&Tau_HadrEoP);
  readerwGwGSFwPFMVA_EC->AddVariable("Tau_VisMass",&Tau_VisMass);
  readerwGwGSFwPFMVA_EC->AddVariable("Tau_GammaEtaMom",&Tau_GammaEtaMom);
  readerwGwGSFwPFMVA_EC->AddVariable("Tau_GammaPhiMom",&Tau_GammaPhiMom);
  readerwGwGSFwPFMVA_EC->AddVariable("Tau_GammaEnFrac",&Tau_GammaEnFrac);
  readerwGwGSFwPFMVA_EC->SetVerbose(kTRUE);
  if(discriminator == "") readerwGwGSFwPFMVA_EC->BookMVA("BDT",WeightwGwGSFwPFMVA_EC ); 


  cout<< "Number of entries : "<<nEntries<<endl;

  for (int iEntry = 0; iEntry<nEntries ; iEntry++){
    if(iEntry%100000==0) cout << iEntry << endl;
    inputTree->GetEntry(iEntry);

    t_NoEleMatch_Barrel = 1;
    t_woG_Barrel = 1;
    t_wGwoGSF_Barrel = 1;
    t_wGwGSFwoPFMVA_Barrel = 1;
    t_wGwGSFwPFMVA_Barrel = 1;
    t_NoEleMatch_Endcap = 1;
    t_woG_Endcap = 1;
    t_wGwoGSF_Endcap = 1;
    t_wGwGSFwoPFMVA_Endcap = 1;
    t_wGwGSFwPFMVA_Endcap = 1;

    if(discriminator == "-AntiEMed" && Tau_AntiEMedium<0.5)continue; 
    if (Type.find("Signal")!=std::string::npos && Tau_GenHadMatch!=1)continue;
    if (Type.find("Backgrd")!=std::string::npos && Tau_GenEleMatch!=1)continue;

    t_Elec_GenEleMatch = Elec_GenEleMatch;
    t_Elec_GenHadMatch = Elec_GenHadMatch;
    t_Tau_GenEleMatch = Tau_GenEleMatch;
    t_Tau_GenHadMatch = Tau_GenHadMatch;

//     Tau_AbsEta = TMath::Abs(Tau_Eta);
    Tau_AbsEta = TMath::Abs(Tau_EtaAtEcalEntrance);
    if(Elec_AbsEta < 1.479 && Tau_AbsEta<1.479){
      if(Tau_GsfEleMatch<0.5) t_NoEleMatch_Barrel = readerNoEleMatch_BL->EvaluateMVA("BDT");
      else if(Tau_NumGammaCands==0) t_woG_Barrel = readerwoG_BL->EvaluateMVA("BDT");
      else if(Tau_NumGammaCands>0 && Tau_HasGsf<0.5)t_wGwoGSF_Barrel = readerwGwoGSF_BL->EvaluateMVA("BDT");
      else if(Tau_NumGammaCands>0 && Tau_HasGsf>0.5 && Elec_PFMvaOutput<-0.1)t_wGwGSFwoPFMVA_Barrel = readerwGwGSFwoPFMVA_BL->EvaluateMVA("BDT");
      else if(Tau_NumGammaCands>0 && Tau_HasGsf>0.5 && Elec_PFMvaOutput>-0.1&& discriminator == "")t_wGwGSFwPFMVA_Barrel = readerwGwGSFwPFMVA_BL->EvaluateMVA("BDT");
    }//For barrel
    else{
      if(Tau_GsfEleMatch<0.5) t_NoEleMatch_Endcap = readerNoEleMatch_EC->EvaluateMVA("BDT");
      else if(Tau_NumGammaCands==0) t_woG_Endcap = readerwoG_EC->EvaluateMVA("BDT");
      else if(Tau_NumGammaCands>0 && Tau_HasGsf<0.5)t_wGwoGSF_Endcap = readerwGwoGSF_EC->EvaluateMVA("BDT");
      else if(Tau_NumGammaCands>0 && Tau_HasGsf>0.5 && Elec_PFMvaOutput<-0.1)t_wGwGSFwoPFMVA_Endcap = readerwGwGSFwoPFMVA_EC->EvaluateMVA("BDT");
      else if(Tau_NumGammaCands>0 && Tau_HasGsf>0.5 && Elec_PFMvaOutput>-0.1 && discriminator == "")t_wGwGSFwPFMVA_Endcap = readerwGwGSFwPFMVA_EC->EvaluateMVA("BDT");
    }//For Endcap

    mytree->Fill();

  }

  mytree->Write();
  
  cout<<"Creating file : "<<outputFileName.data()<<endl;
  inputFile->Close();
  outputFile->Close();
  return;
}

void getAllMVAOutput(){
  MVAOutput("","Signal");
  MVAOutput("","Backgrd");
//   MVAOutput("-AntiEMed","Signal");
//   MVAOutput("-AntiEMed","Backgrd");
}

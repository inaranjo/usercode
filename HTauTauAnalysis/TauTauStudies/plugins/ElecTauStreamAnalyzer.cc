#include "Bianchi/TauTauStudies/interface/ElecTauStreamAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "TauAnalysis/CandidateTools/interface/CompositePtrCandidateT1T2MEtProducer.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include <DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h>

#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include <DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

#include <vector>
#include <utility>
#include <map>

using namespace std;
using namespace reco;

typedef std::map<double, math::XYZTLorentzVectorD ,ElecTauStreamAnalyzer::more>::iterator CImap;

ElecTauStreamAnalyzer::ElecTauStreamAnalyzer(const edm::ParameterSet & iConfig){

  diTauTag_          = iConfig.getParameter<edm::InputTag>("diTaus");
  jetsTag_           = iConfig.getParameter<edm::InputTag>("jets");
  newJetsTag_        = iConfig.getParameter<edm::InputTag>("newJets");
  metTag_            = iConfig.getParameter<edm::InputTag>("met");
  rawMetTag_         = iConfig.getParameter<edm::InputTag>("rawMet");
  electronsTag_      = iConfig.getParameter<edm::InputTag>("electrons");
  electronsRelTag_   = iConfig.getParameter<edm::InputTag>("electronsRel");
  verticesTag_       = iConfig.getParameter<edm::InputTag>("vertices");
  triggerResultsTag_ = iConfig.getParameter<edm::InputTag>("triggerResults"); 
  isMC_              = iConfig.getParameter<bool>("isMC");
  deltaRLegJet_      = iConfig.getUntrackedParameter<double>("deltaRLegJet",0.3);
  minCorrPt_         = iConfig.getUntrackedParameter<double>("minCorrPt",10.);
  minJetID_          = iConfig.getUntrackedParameter<double>("minJetID",0.5);

  //inputFileNameX0BL_ = iConfig.getParameter<edm::FileInPath>("inputFileNameX0BL");
  //inputFileName11BL_ = iConfig.getParameter<edm::FileInPath>("inputFileName11BL");
  //inputFileName01BL_ = iConfig.getParameter<edm::FileInPath>("inputFileName01BL");
  //inputFileNameX0EC_ = iConfig.getParameter<edm::FileInPath>("inputFileNameX0EC");
  //inputFileName11EC_ = iConfig.getParameter<edm::FileInPath>("inputFileName11EC");
  //inputFileName01EC_ = iConfig.getParameter<edm::FileInPath>("inputFileName01EC");


  verbose_           = iConfig.getUntrackedParameter<bool>("verbose",false);

}

void ElecTauStreamAnalyzer::beginJob(){

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","qqH tree");

  tRandom_ = new TRandom3();
 
  jetsBtagHE_  = new std::vector< double >();
  jetsBtagHP_  = new std::vector< double >();

  jetsChEfraction_  = new std::vector< float >();
  jetsChNfraction_  = new std::vector< float >();
  jetMoments_       = new std::vector< float >();

  gammadR_       = new std::vector< float >();
  gammadEta_     = new std::vector< float >();
  gammadPhi_     = new std::vector< float >();
  gammaPt_       = new std::vector< float >();

  tauXTriggers_= new std::vector< int >();
  triggerBits_ = new std::vector< int >();

  jetsP4_          = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  jetsIDP4_        = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  jetsIDUpP4_      = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  jetsIDDownP4_    = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  jetsIDL1OffsetP4_= new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genJetsIDP4_     = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();

  diTauVisP4_   = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  diTauCAP4_    = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  diTauICAP4_   = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  diTauSVfitP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();

  diTauLegsP4_    = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genDiTauLegsP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genTausP4_      = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  METP4_          = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genMETP4_       = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genVP4_         = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  
  extraElectrons_   = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  pfElectrons_      = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();

  //antiE_  = new AntiElectronIDMVA();
  //antiE_->Initialize("BDT",
  //	     inputFileNameX0BL_.fullPath().data(),
  //	     inputFileName11BL_.fullPath().data(),
  //	     inputFileName01BL_.fullPath().data(),
  //	     inputFileNameX0EC_.fullPath().data(),
  //	     inputFileName11EC_.fullPath().data(),
  //	     inputFileName01EC_.fullPath().data()
  //	     );

  std::vector< float > Summer11Lumi ;
  Double_t Summer11Lumi_f[35] = {
    1.45346E-01,
    6.42802E-02,
    6.95255E-02,
    6.96747E-02,
    6.92955E-02,
    6.84997E-02,
    6.69528E-02,
    6.45515E-02,
    6.09865E-02,
    5.63323E-02,
    5.07322E-02,
    4.44681E-02,
    3.79205E-02,
    3.15131E-02,
    2.54220E-02,
    2.00184E-02,
    1.53776E-02,
    1.15387E-02,
    8.47608E-03,
    6.08715E-03,
    4.28255E-03,
    2.97185E-03,
    2.01918E-03,
    1.34490E-03,
    8.81587E-04,
    5.69954E-04,
    3.61493E-04,
    2.28692E-04,
    1.40791E-04,
    8.44606E-05,
    5.10204E-05,
    3.07802E-05,
    1.81401E-05,
    1.00201E-05,
    5.80004E-06
  };

  std::vector< float > Data2011Lumi;
  Double_t Data2011Lumi_f[35] = {
    0.00290212,
    0.0123985, 
    0.0294783, 
    0.0504491, 
    0.0698525, 
    0.0836611, 
    0.0905799, 
    0.0914388, 
    0.0879379, 
    0.0817086, 
    0.073937, 
    0.0653785,
    0.0565162,
    0.047707, 
    0.0392591,
    0.0314457,
    0.0244864,
    0.018523, 
    0.013608, 
    0.00970977, 
    0.00673162, 
    0.00453714, 
    0.00297524, 
    0.00189981, 
    0.00118234, 
    0.000717854,
    0.00042561, 
    0.000246653,
    0.000139853,
    7.76535E-05,
    4.22607E-05,
    2.25608E-05,
    1.18236E-05,
    6.0874E-06, 
    6.04852E-06
  };

  std::vector< float > Data2011LumiExt;
  Double_t Data2011LumiExt_f[50] = {
    0.00290212,
    0.0123985,
    0.0294783,
    0.0504491,
    0.0698525,
    0.0836611,
    0.0905799,
    0.0914388,
    0.0879379,
    0.0817086,
    0.073937,
    0.0653785,
    0.0565162,
    0.047707,
    0.0392591,
    0.0314457,
    0.0244864,
    0.018523,
    0.013608,
    0.00970977,
    0.00673162,
    0.00453714,
    0.00297524,
    0.00189981,
    0.00118234,
    0.000717854,
    0.00042561,
    0.000246653,
    0.000139853,
    7.76535E-05,
    4.22607E-05,
    2.25608E-05,
    1.18236E-05,
    6.0874E-06,
    6.04852E-06,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  };

  std::vector< float > Fall11Lumi ;
  Double_t Fall11Lumi_f[50] = {
    0.003388501,
    0.010357558,
    0.024724258,
    0.042348605,
    0.058279812,
    0.068851751,
    0.072914824,
    0.071579609,
    0.066811668,
    0.060672356,
    0.054528356,
    0.04919354,
    0.044886042,
    0.041341896,
    0.0384679,
    0.035871463,
    0.03341952,
    0.030915649,
    0.028395374,
    0.025798107,
    0.023237445,
    0.020602754,
    0.0180688,
    0.015559693,
    0.013211063,
    0.010964293,
    0.008920993,
    0.007080504,
    0.005499239,
    0.004187022,
    0.003096474,
    0.002237361,
    0.001566428,
    0.001074149,
    0.000721755,
    0.000470838,
    0.00030268,
    0.000184665,
    0.000112883,
    6.74043E-05,
    3.82178E-05,
    2.22847E-05,
    1.20933E-05,
    6.96173E-06,
    3.4689E-06,
    1.96172E-06,
    8.49283E-07,
    5.02393E-07,
    2.15311E-07,
    9.56938E-08
  };

  std::vector< float > Data2011TruthLumi;
  Double_t Data2011TruthLumi_f[50] = {
    0, 
    6.48087e-05, 
    0.00122942, 
    0.0107751, 
    0.0572544, 
    0.110336, 
    0.122622, 
    0.113354, 
    0.0991184, 
    0.0907195, 
    0.0813241, 
    0.0748131, 
    0.0703376, 
    0.0625598, 
    0.0487456, 
    0.030941, 
    0.0158153, 
    0.00657034, 
    0.00232069, 
    0.000782941, 
    0.000240306, 
    6.13863e-05, 
    1.36142e-05, 
    4.99913e-07, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  };

  for( int i=0; i<35; i++) {
    Data2011Lumi.push_back(Data2011Lumi_f[i]);
    Summer11Lumi.push_back(Summer11Lumi_f[i]);
  }
  for( int i=0; i<50; i++) {
    Data2011TruthLumi.push_back(Data2011TruthLumi_f[i]);
    Data2011LumiExt.push_back(Data2011LumiExt_f[i]);
    Fall11Lumi.push_back(Fall11Lumi_f[i]);
  }
  
  //LumiWeights_ = edm::LumiReWeighting(Summer11Lumi, Data2011Lumi);
  LumiWeights_ = edm::LumiReWeighting(Fall11Lumi, Data2011LumiExt);
    
  tree_->Branch("jetsP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsP4_);
  tree_->Branch("jetsIDP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsIDP4_);
  tree_->Branch("jetsIDUpP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsIDUpP4_);
  tree_->Branch("jetsIDDownP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsIDDownP4_);
  tree_->Branch("jetsIDL1OffsetP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsIDL1OffsetP4_);
  tree_->Branch("genJetsIDP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genJetsIDP4_);
  
  tree_->Branch("jetsBtagHE","std::vector<double> ",&jetsBtagHE_);
  tree_->Branch("jetsBtagHP","std::vector<double> ",&jetsBtagHP_);

  tree_->Branch("jetMoments","std::vector<float> ",&jetMoments_);

  tree_->Branch("gammadR","std::vector<float> ",&gammadR_);
  tree_->Branch("gammadEta","std::vector<float> ",&gammadEta_);
  tree_->Branch("gammadPhi","std::vector<float> ",&gammadPhi_);
  tree_->Branch("gammaPt","std::vector<float> ",&gammaPt_);

  tree_->Branch("jetsChEfraction","std::vector<float>",&jetsChEfraction_);
  tree_->Branch("jetsChNfraction","std::vector<float>",&jetsChNfraction_);

  tree_->Branch("tauXTriggers","std::vector<int>",&tauXTriggers_);
  tree_->Branch("triggerBits","std::vector<int>",&triggerBits_);

  tree_->Branch("diTauVisP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",  &diTauVisP4_);
  tree_->Branch("diTauCAP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",   &diTauCAP4_);
  tree_->Branch("diTauICAP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",   &diTauICAP4_);
  tree_->Branch("diTauSVfitP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&diTauSVfitP4_);

  tree_->Branch("diTauNSVfitMass",       &diTauNSVfitMass_,"diTauNSVfitMass/F");
  tree_->Branch("diTauNSVfitMassErrUp",  &diTauNSVfitMassErrUp_,"diTauNSVfitMassErrUp/F");
  tree_->Branch("diTauNSVfitMassErrDown",&diTauNSVfitMassErrDown_,"diTauNSVfitMassErrDown/F");

  tree_->Branch("diTauLegsP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&diTauLegsP4_);
  tree_->Branch("genDiTauLegsP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genDiTauLegsP4_);
  tree_->Branch("genTausP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genTausP4_);

  tree_->Branch("METP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&METP4_);
  tree_->Branch("genMETP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genMETP4_);
  tree_->Branch("genVP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genVP4_);
  tree_->Branch("genDecay",&genDecay_,"genDecay/I");
  tree_->Branch("extraElectrons","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&extraElectrons_);
  tree_->Branch("pfElectrons","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&pfElectrons_);
  tree_->Branch("leg1GenDecay",&leg1GenDecay_,"leg1GenDecay/I");
  tree_->Branch("ptHighestLeg1Matched",&ptHighestLeg1Matched_,"ptHighestLeg1Matched/F");
  tree_->Branch("ratioPtMaxMatchPtGenElec",&ratioPtMaxMatchPtGenElec_,"ratioPtMaxMatchPtGenElec/F");

  tree_->Branch("sumEt",&sumEt_,"sumEt/F");
  tree_->Branch("mTauTauMin",&mTauTauMin_,"mTauTauMin/F");
  tree_->Branch("MtLeg1",&MtLeg1_,"MtLeg1/F");
  tree_->Branch("pZeta",&pZeta_,"pZeta/F");
  tree_->Branch("pZetaVis",&pZetaVis_,"pZetaVis/F");
  tree_->Branch("pZetaSig",&pZetaSig_,"pZetaSig/F");

  tree_->Branch("chIsoLeg1v1",&chIsoLeg1v1_,"chIsoLeg1v1/F");
  tree_->Branch("nhIsoLeg1v1",&nhIsoLeg1v1_,"nhIsoLeg1v1/F");
  tree_->Branch("phIsoLeg1v1",&phIsoLeg1v1_,"phIsoLeg1v1/F");
  tree_->Branch("elecIsoLeg1v1",&elecIsoLeg1v1_,"elecIsoLeg1v1/F");
  tree_->Branch("muIsoLeg1v1",&muIsoLeg1v1_,"muIsoLeg1v1/F");
  tree_->Branch("chIsoPULeg1v1",&chIsoPULeg1v1_,"chIsoPULeg1v1/F");
  tree_->Branch("nhIsoPULeg1v1",&nhIsoPULeg1v1_,"nhIsoPULeg1v1/F");
  tree_->Branch("phIsoPULeg1v1",&phIsoPULeg1v1_,"phIsoPULeg1v1/F");

  tree_->Branch("chIsoLeg1v2",&chIsoLeg1v2_,"chIsoLeg1v2/F");
  tree_->Branch("nhIsoLeg1v2",&nhIsoLeg1v2_,"nhIsoLeg1v2/F");
  tree_->Branch("phIsoLeg1v2",&phIsoLeg1v2_,"phIsoLeg1v2/F");
  tree_->Branch("elecIsoLeg1v2",&elecIsoLeg1v2_,"elecIsoLeg1v2/F");
  tree_->Branch("muIsoLeg1v2",&muIsoLeg1v2_ ,"muIsoLeg1v2/F");  
  tree_->Branch("chIsoPULeg1v2",&chIsoPULeg1v2_,"chIsoPULeg1v2/F");
  tree_->Branch("nhIsoPULeg1v2",&nhIsoPULeg1v2_,"nhIsoPULeg1v2/F");
  tree_->Branch("phIsoPULeg1v2",&phIsoPULeg1v2_,"phIsoPULeg1v2/F");

  tree_->Branch("chIsoEELeg1v2"  ,&chIsoEELeg1v2_,  "chIsoEELeg1v2/F");
  tree_->Branch("nhIsoEELeg1v2"  ,&nhIsoEELeg1v2_,  "nhIsoEELeg1v2/F");
  tree_->Branch("phIsoEELeg1v2"  ,&phIsoEELeg1v2_,  "phIsoEELeg1v2/F");
  tree_->Branch("elecIsoEELeg1v2",&elecIsoEELeg1v2_,"elecIsoEELeg1v2/F");
  tree_->Branch("muIsoEELeg1v2"  ,&muIsoEELeg1v2_ , "muIsoEELeg1v2/F");  
  tree_->Branch("nhIsoEEPULeg1v2",&nhIsoEEPULeg1v2_,"nhIsoEEPULeg1v2/F");
  tree_->Branch("phIsoEEPULeg1v2",&phIsoEEPULeg1v2_,"phIsoEEPULeg1v2/F");

  tree_->Branch("chIsoLeg2",&chIsoLeg2_,"chIsoLeg2/F");
  tree_->Branch("nhIsoLeg2",&nhIsoLeg2_,"nhIsoLeg2/F");
  tree_->Branch("phIsoLeg2",&phIsoLeg2_,"phIsoLeg2/F");
  tree_->Branch("dxy1",&dxy1_,"dxy1/F");
  tree_->Branch("dxy2",&dxy2_,"dxy2/F");
  tree_->Branch("dz1",&dz1_,"dz1/F");
  tree_->Branch("dz2",&dz2_,"dz2/F");

  // electron specific variables
  tree_->Branch("nBrehm",&nBrehm_,"nBrehm/F");
  tree_->Branch("likelihood",&likelihood_,"likelihood/F");
  tree_->Branch("nHits",&nHits_,"nHits/F");
  tree_->Branch("sihih",&sihih_,"sihih/F");
  tree_->Branch("dPhi",&dPhi_,"dPhi/F");
  tree_->Branch("dEta",&dEta_,"dEta/F");
  tree_->Branch("HoE",&HoE_,"HoE/F");
  tree_->Branch("EoP",&EoP_,"EoP/F");
  tree_->Branch("fbrem",&fbrem_,"fbrem/F");
  //tree_->Branch("isEleLikelihoodID",&isEleLikelihoodID_,"isEleLikelihoodID/I");
  //tree_->Branch("isEleCutBasedID",&isEleCutBasedID_,"isEleCutBasedID/I");
  tree_->Branch("elecFlag",&elecFlag_,"elecFlag/I");
  tree_->Branch("elecVetoRelIso",&elecVetoRelIso_,"elecVetoRelIso/F");
  tree_->Branch("hasKft",&hasKft_,"hasKft/I");
  tree_->Branch("antiConv",&antiConv_,"antiConv/F");


  tree_->Branch("run",&run_,"run/l");
  tree_->Branch("event",&event_,"event/l");
  tree_->Branch("lumi",&lumi_,"lumi/l");
  tree_->Branch("numPV",&numPV_,"numPV/F");
  tree_->Branch("numOfDiTaus",&numOfDiTaus_,"numOfDiTaus/I");
  tree_->Branch("numOfLooseIsoDiTaus",&numOfLooseIsoDiTaus_,"numOfLooseIsoDiTaus/I");
  tree_->Branch("decayMode",&decayMode_,"decayMode/I");
  tree_->Branch("tightestHPSWP",&tightestHPSWP_,"tightestHPSWP/I");
  tree_->Branch("tightestCiCWP",&tightestCiCWP_,"tightestCiCWP/I");
  tree_->Branch("tightestCutBasedWP",&tightestCutBasedWP_,"tightestCutBasedWP/I");
  tree_->Branch("tightestMVAWP",&tightestMVAWP_,"tightestMVAWP/I");
  tree_->Branch("tightestDanieleMVAWP",&tightestDanieleMVAWP_,"tightestDanieleMVAWP/F");
  tree_->Branch("tightestAntiECutsWP",&tightestAntiECutsWP_,"tightestAntiECutsWP/I");
  tree_->Branch("tightestAntiEMVAWP",&tightestAntiEMVAWP_,"tightestAntiEMVAWP/I");
  tree_->Branch("tightestHPSDBWP",&tightestHPSDBWP_,"tightestHPSDBWP/I");
  tree_->Branch("visibleTauMass",&visibleTauMass_,"visibleTauMass/F");

  tree_->Branch("leadPFChargedHadrTrackPt",&leadPFChargedHadrTrackPt_,"leadPFChargedHadrTrackPt/F");
  tree_->Branch("leadPFChargedHadrTrackP", &leadPFChargedHadrTrackP_,"leadPFChargedHadrTrackP/F");
  tree_->Branch("leadPFChargedHadrPt",&leadPFChargedHadrPt_,"leadPFChargedHadrPt/F");
  tree_->Branch("leadPFChargedHadrP", &leadPFChargedHadrP_,"leadPFChargedHadrP/F");
  tree_->Branch("leadPFChargedHadrMva",&leadPFChargedHadrMva_,"leadPFChargedHadrMva/F");
  tree_->Branch("leadPFChargedHadrHcalEnergy",&leadPFChargedHadrHcalEnergy_,"leadPFChargedHadrHcalEnergy/F");
  tree_->Branch("leadPFChargedHadrEcalEnergy",&leadPFChargedHadrEcalEnergy_,"leadPFChargedHadrEcalEnergy/F");
  tree_->Branch("leadPFCandMva",&leadPFCandMva_,"leadPFCandMva/F");
  tree_->Branch("leadPFCandHcalEnergy",&leadPFCandHcalEnergy_,"leadPFCandHcalEnergy/F");
  tree_->Branch("leadPFCandEcalEnergy",&leadPFCandEcalEnergy_,"leadPFCandEcalEnergy/F");
  tree_->Branch("leadPFCandPt",&leadPFCandPt_,"leadPFCandPt/F");
  tree_->Branch("leadPFCandP",&leadPFCandP_,"leadPFCandP/F");
  tree_->Branch("emFraction",&emFraction_,"emFraction/F");
  tree_->Branch("hasGsf",&hasGsf_,"hasGsf/F");
  tree_->Branch("signalPFChargedHadrCands",&signalPFChargedHadrCands_,"signalPFChargedHadrCands/I");
  tree_->Branch("signalPFGammaCands",&signalPFGammaCands_,"signalPFGammaCands/I");
  //tree_->Branch("mvaAntiE",&mvaAntiE_,"mvaAntiE/F");

  tree_->Branch("isTauLegMatched",&isTauLegMatched_,"isTauLegMatched/I");
  tree_->Branch("isElecLegMatched",&isElecLegMatched_,"isElecLegMatched/I");

  tree_->Branch("diTauCharge",&diTauCharge_,"diTauCharge/F");
  tree_->Branch("rhoFastJet",&rhoFastJet_,"rhoFastJet/F");
  tree_->Branch("rhoNeutralFastJet",&rhoNeutralFastJet_,"rhoNeutralFastJet/F");
  tree_->Branch("mcPUweight",&mcPUweight_,"mcPUweight/F");
  tree_->Branch("embeddingWeight",&embeddingWeight_,"embeddingWeight/F");
  tree_->Branch("nPUVertices",&nPUVertices_,"nPUVertices/I");
  tree_->Branch("nPUVerticesM1",&nPUVerticesM1_,"nPUVerticesM1/I");
  tree_->Branch("nPUVerticesP1",&nPUVerticesP1_,"nPUVerticesP1/I");
  tree_->Branch("nPUaverage",&nPUaverage_,"nPUaverage/F");


}


ElecTauStreamAnalyzer::~ElecTauStreamAnalyzer(){
  delete jetsP4_; delete jetsIDP4_; delete jetsIDL1OffsetP4_; delete jetsIDUpP4_; delete jetsIDDownP4_; 
  delete METP4_; delete diTauVisP4_; delete diTauCAP4_; delete diTauICAP4_; 
  delete diTauSVfitP4_; delete genVP4_;
  delete diTauLegsP4_; delete jetsBtagHE_; delete jetsBtagHP_; delete tauXTriggers_; delete triggerBits_;
  delete genJetsIDP4_; delete genDiTauLegsP4_; delete genMETP4_;delete extraElectrons_; delete pfElectrons_;
  delete genTausP4_;
  delete tRandom_ ; delete jetsChNfraction_; delete jetsChEfraction_; delete jetMoments_;
  delete gammadR_ ; delete gammadPhi_; delete gammadEta_; delete gammaPt_;
  //delete antiE_;
}

void ElecTauStreamAnalyzer::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup){

  jetsP4_->clear();
  jetsIDP4_->clear();
  jetsIDUpP4_->clear();
  jetsIDDownP4_->clear();
  jetsIDL1OffsetP4_->clear();
  diTauVisP4_->clear();
  diTauCAP4_->clear();
  diTauICAP4_->clear();
  diTauSVfitP4_->clear();
  diTauLegsP4_->clear();
  METP4_->clear();
  genVP4_->clear();
  extraElectrons_->clear();
  pfElectrons_->clear();
  jetsChNfraction_->clear();
  jetsChEfraction_->clear();
  jetMoments_->clear();

  gammadR_->clear();
  gammadEta_->clear();
  gammadPhi_->clear();
  gammaPt_->clear();

  genJetsIDP4_->clear();
  genDiTauLegsP4_->clear();
  genTausP4_->clear();
  genMETP4_->clear();
  jetsBtagHE_->clear();
  jetsBtagHP_->clear();
  tauXTriggers_->clear();
  triggerBits_->clear();

  
  edm::Handle<PATElecTauPairCollection> diTauHandle;
  iEvent.getByLabel(diTauTag_,diTauHandle);
  if( !diTauHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No diTau label available \n";
  const PATElecTauPairCollection* diTaus = diTauHandle.product();

  edm::Handle<pat::JetCollection> jetsHandle;
  iEvent.getByLabel(jetsTag_,jetsHandle);
  if( !jetsHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No jets label available \n";
  const pat::JetCollection* jets = jetsHandle.product();

  edm::Handle<pat::JetCollection> newJetsHandle;
  const pat::JetCollection* newJets = 0;
  if(newJetsTag_.label()!=""){
    iEvent.getByLabel(newJetsTag_,newJetsHandle);
    if( !newJetsHandle.isValid() )  
      edm::LogError("DataNotAvailable")
	<< "No newJets label available \n";
    newJets = newJetsHandle.product();
  }

  edm::Handle<reco::PFCandidateCollection> pfHandle;
  iEvent.getByLabel("particleFlow",pfHandle);
  if( !pfHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No pf particles label available \n";
  const reco::PFCandidateCollection* pfCandidates = pfHandle.product();

  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByLabel(  verticesTag_ ,pvHandle);
  if( !pvHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No PV label available \n";
  const reco::VertexCollection* vertexes = pvHandle.product();

  std::vector<float> vtxZ;
  for(unsigned int k = 0; k<vertexes->size(); k++){
    vtxZ.push_back(((*vertexes)[k].position()).z());
  }
  
  if(verbose_){
    cout << "List of vertexes " << endl;
    for(unsigned int k = 0; k<vertexes->size(); k++){
      cout << "Vtx[" << k << "] (x,y,z) = (" << ((*vertexes)[k].position()).x()
	   << "," << ((*vertexes)[k].position()).y() << "," << ((*vertexes)[k].position()).z() << ")"
	   << endl;
    }
  }
  
  numPV_ = vertexes->size();

  edm::Handle<pat::METCollection> metHandle;
  iEvent.getByLabel( metTag_ ,metHandle);
  if( !metHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No MET label available \n";
  const pat::METCollection* met = metHandle.product();

  edm::Handle<pat::METCollection> rawMetHandle;
  iEvent.getByLabel( rawMetTag_, rawMetHandle);
  if( !rawMetHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No raw MET label available \n";
  const pat::METCollection* rawMet = rawMetHandle.product();

  edm::Handle<pat::TriggerEvent> triggerHandle;
  iEvent.getByLabel(triggerResultsTag_, triggerHandle);
  if( !triggerHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No Trigger label available \n";
  const pat::TriggerEvent* trigger = triggerHandle.product();

  edm::Handle<pat::TriggerObjectStandAloneCollection > triggerObjsHandle;
  iEvent.getByLabel(edm::InputTag("patTrigger"),triggerObjsHandle);
  if( !triggerObjsHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No Trigger Objects label available \n";
  const pat::TriggerObjectStandAloneCollection* triggerObjs = triggerObjsHandle.product();
    
  genDecay_ = -99;
  edm::Handle<reco::GenParticleCollection> genHandle;
  const reco::GenParticleCollection* genParticles = 0;
  iEvent.getByLabel(edm::InputTag("genParticles"),genHandle);

  if(genHandle.isValid()){
    genParticles = genHandle.product();
    int index1 = -99;
    int index2 = -99;
    for(unsigned int k = 0; k < genParticles->size(); k ++){
      if( fabs((*genParticles)[k].pdgId()) != 15 || 
	  (*genParticles)[k].status()!=2) continue;
      if(index1<0) 
	index1 = k;
      else 
	index2 = k;
    }
    if(index1>=0 && index2>=0){
      if((*genParticles)[index1].pt()<(*genParticles)[index2].pt()){
	int indexBkp = index1;
	index1 = index2;
	index2 = indexBkp;
      }
      genTausP4_->push_back((*genParticles)[index1].p4());
      genTausP4_->push_back((*genParticles)[index2].p4());
    }
  }

  if(isMC_){
    iEvent.getByLabel(edm::InputTag("genParticles"),genHandle);
    if( !genHandle.isValid() )  
      edm::LogError("DataNotAvailable")
	<< "No gen particles label available \n";
    genParticles = genHandle.product();
    for(unsigned int k = 0; k < genParticles->size(); k ++){
      if( !( (*genParticles)[k].pdgId() == 23 || abs((*genParticles)[k].pdgId()) == 24 || (*genParticles)[k].pdgId() == 25) || (*genParticles)[k].status()!=3)
	continue;
      if(verbose_) cout << "Boson status, pt,phi " << (*genParticles)[k].status() << "," << (*genParticles)[k].pt() << "," << (*genParticles)[k].phi() << endl;
      
      genVP4_->push_back( (*genParticles)[k].p4() );
      genDecay_ = (*genParticles)[k].pdgId();
      
      int breakLoop = 0;
      for(unsigned j = 0; j< ((*genParticles)[k].daughterRefVector()).size() && breakLoop!=1; j++){
	if( abs(((*genParticles)[k].daughterRef(j))->pdgId()) == 11 ){
	  genDecay_ *= 11;
	  breakLoop = 1;
	}
	if( abs(((*genParticles)[k].daughterRef(j))->pdgId()) == 13 ){
	  genDecay_ *= 13;
	  breakLoop = 1;
	}
	if( abs(((*genParticles)[k].daughterRef(j))->pdgId()) == 15 ){
	  genDecay_ *= 15;
	  breakLoop = 1;
	}
      }
      if(verbose_) cout << "Decays to pdgId " << genDecay_/(*genParticles)[k].pdgId()  << endl;
      break;
    }
  }

  edm::Handle<reco::GenJetCollection> tauGenJetsHandle;
  edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;
  nPUVertices_       = -99;
  nPUVerticesP1_     = -99;
  nPUVerticesM1_     = -99;

  nPUtruth_          = -99;

  const reco::GenJetCollection* tauGenJets = 0;
  if(isMC_){
    iEvent.getByLabel(edm::InputTag("tauGenJetsSelectorAllHadrons"),tauGenJetsHandle);
    if( !tauGenJetsHandle.isValid() )  
      edm::LogError("DataNotAvailable")
	<< "No gen jet label available \n";
    tauGenJets = tauGenJetsHandle.product();

    // PU infos
    iEvent.getByType(puInfoH);
    if(puInfoH.isValid()){
      for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++){

	//cout << "Bunc crossing " << it->getBunchCrossing() << endl;
	if(it->getBunchCrossing() == 0 ) 
	  nPUVertices_  = it->getPU_NumInteractions();
	else if(it->getBunchCrossing() == -1)  
	  nPUVerticesM1_= it->getPU_NumInteractions();
	else if(it->getBunchCrossing() == +1)  
	  nPUVerticesP1_= it->getPU_NumInteractions();

	nPUtruth_ = it->getTrueNumInteractions();

      }
    }
  }
  if(verbose_){
    cout << "Num of PU = "          << nPUVertices_ << endl;
    cout << "Num of OOT PU = "      << nPUVerticesM1_+nPUVerticesP1_<< endl;
    cout << "Num of true interactions = " << nPUtruth_ << endl;
  }

  mcPUweight_ = LumiWeights_.weight( nPUVertices_ );

  edm::Handle<double> rhoFastJetHandle;
  iEvent.getByLabel(edm::InputTag("kt6PFJetsCentral","rho",""), rhoFastJetHandle);
  if( !rhoFastJetHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No rho label available \n";
  rhoFastJet_ = (*rhoFastJetHandle);

  edm::Handle<double> embeddingWeightHandle;
  iEvent.getByLabel(edm::InputTag("generator","weight",""), embeddingWeightHandle);
  embeddingWeight_ = embeddingWeightHandle.isValid() ? (*embeddingWeightHandle) : 1.0;


  edm::Handle<double> rhoNeutralFastJetHandle;
  iEvent.getByLabel(edm::InputTag("kt6PFJetsNeutral","rho", ""), rhoNeutralFastJetHandle);
  if( !rhoNeutralFastJetHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No rho neutral label available \n";
  rhoNeutralFastJet_ = (*rhoNeutralFastJetHandle);

  edm::Handle<pat::ElectronCollection> electronsHandle;
  iEvent.getByLabel( electronsTag_ ,electronsHandle);
  if( !electronsHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No electrons label available \n";
  const pat::ElectronCollection* electrons = electronsHandle.product();
  if(electrons->size()<1){
    cout << " No electrons !!! " << endl;
    return;
  } else if(electrons->size()>1 && verbose_){
    cout << "WARNING: "<< electrons->size() << "  electrons found in the event !!! We will select only one" << endl;
  }
  edm::Handle<pat::ElectronCollection> electronsRelHandle;
  iEvent.getByLabel(electronsRelTag_,electronsRelHandle);
  if( !electronsRelHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No electronsRel label available \n";
  const pat::ElectronCollection* electronsRel = electronsRelHandle.product();
  if(electronsRel->size()<1){
    cout << " No electronsRel !!! " << endl;
    //return;
  } else if(electronsRel->size()>1 && verbose_){
    cout << "WARNING: "<< electronsRel->size() << "  electronsRel found in the event !!! We will select only one" << endl;
  }
  
  const PATElecTauPair *theDiTau = 0;
  if(diTaus->size()<1){
    cout << " No diTau !!! " << endl;
    return;
  } else if(diTaus->size()>1 && verbose_){
    cout << "WARNING: "<< diTaus->size() << "  diTaus found in the event !!! We will select only one" << endl;
  }
  numOfDiTaus_ = diTaus->size();

  unsigned int elecIndex = 0;
  elecFlag_       = 0;
  elecVetoRelIso_ = 0;

  bool found = false;
  for(unsigned int i=0; i<electronsRel->size(); i++){
    for(unsigned int j=0; j<electrons->size(); j++){ 

      if(found) continue;

      bool mvaPreselectionNoIsoElecI = 
	(*electrons)[j].userInt("antiConv")>0.5 && (*electrons)[j].userFloat("nHits")<=0 &&
	(((*electrons)[j].isEB()
	  && (*electrons)[j].userFloat("sihih") < 0.01 
	  && TMath::Abs((*electrons)[j].userFloat("dEta")) < 0.007
	  && TMath::Abs((*electrons)[j].userFloat("dPhi")) < 0.15
	  && (*electrons)[j].userFloat("HoE") < 0.12
	  ) ||
	 ((*electrons)[j].isEE() 
	  && (*electrons)[j].userFloat("sihih") < 0.03 
	  && TMath::Abs((*electrons)[j].userFloat("dEta")) < 0.009
	  && TMath::Abs((*electrons)[j].userFloat("dPhi")) < 0.10
	  && (*electrons)[j].userFloat("HoE") < 0.10
	  ));

      if( Geom::deltaR( (*electronsRel)[i].p4(),(*electrons)[j].p4())>0.3
	  && (*electronsRel)[i].charge()*(*electrons)[j].charge()<0
	  && (*electrons)[j].userFloat("PFRelIsoDB04v3")<0.30 && mvaPreselectionNoIsoElecI
	  && (*electronsRel)[i].userFloat("PFRelIsoDB04v3")<0.20 ){
	elecFlag_       = 1;
	elecVetoRelIso_ = (*electronsRel)[i].userFloat("PFRelIsoDB04v3");
	if(verbose_) cout<< "Two electrons failing diElec veto: flag= " << elecFlag_ << endl;
	found=true;
      }
      else if( Geom::deltaR( (*electronsRel)[i].p4(),(*electrons)[j].p4())>0.3
	       && (*electronsRel)[i].charge()*(*electrons)[j].charge()>0
	       && (*electrons)[j].userFloat("PFRelIsoDB04v3")<0.30 && mvaPreselectionNoIsoElecI
	       && (*electronsRel)[i].userFloat("PFRelIsoDB04v3")<0.20 ){
	elecFlag_       = 2;
	elecVetoRelIso_ = (*electronsRel)[i].userFloat("PFRelIsoDB04v3");
	if(verbose_) cout<< "Two electrons with SS: flag= " << elecFlag_ << endl;
	found=true;
      }
    }
  }  


  // select highest-pt electron that forms at least one di-tau 
  bool found2 = false;
  for(unsigned int j=0; j<electrons->size(); j++){
    if(verbose_) cout << "Ele " << j << ": " << (*electrons)[j].pt() << ", " << (*electrons)[j].eta() << endl;
    for(unsigned int i=0; i< diTaus->size(); i++){

      if(found2) continue;

      if(verbose_) cout << "diTau " << i <<  ": Leg1 " << ((*diTaus)[i].leg1())->pt() 
			<< ", " << ((*diTaus)[i].leg1())->eta() 
			<<  ": Leg2 " << ((*diTaus)[i].leg2())->pt() 
			<< ", " << ((*diTaus)[i].leg2())->eta() 
			<< endl;
      if( Geom::deltaR(((*diTaus)[i].leg1())->p4(),((*electrons)[j]).p4())<0.01 ){
	elecIndex=j;
	found2 = true;
      }
    }
  }
  if(verbose_) cout<< "Selected electron with index " << elecIndex << endl;

  for(unsigned int i=0; i<electronsRel->size(); i++){
    if( Geom::deltaR( (*electronsRel)[i].p4(),(*electrons)[elecIndex].p4())>0.3) 
      extraElectrons_->push_back( (*electronsRel)[i].p4() );
  }

  std::vector<unsigned int> selectedDiTausFromElec;
  for(unsigned int i=0; i< diTaus->size(); i++){
    if( Geom::deltaR(((*diTaus)[i].leg1())->p4(),((*electrons)[elecIndex]).p4())<0.01 ) 
      selectedDiTausFromElec.push_back(i);
  }
  if(verbose_) cout<< "Selection on Electrons has selected " << selectedDiTausFromElec.size() << " diTaus...check tau leg now..." << endl;


  // choose the diTau with most isolated tau leg

  double sumIsoTau = 999.;
  //double highestPt = 0.;
  unsigned int index = 0;

  std::vector<unsigned int> identifiedTaus;
  std::vector<unsigned int> looseIsoTaus;
  //std::vector<unsigned int> not_identifiedTaus;

  for(unsigned int i=0; i<selectedDiTausFromElec.size(); i++){
    const pat::Tau*  tau_i = dynamic_cast<const pat::Tau*>(  ((*diTaus)[i].leg2()).get() );
    if(tau_i->tauID("decayModeFinding")<0.5) 
      continue;
    identifiedTaus.push_back(selectedDiTausFromElec[i]);
    bool isIsolatedTau = false;
    isIsolatedTau = tau_i->tauID("byLooseCombinedIsolationDeltaBetaCorr")>0.5;
    if(isIsolatedTau) looseIsoTaus.push_back(selectedDiTausFromElec[i]);
  }

  numOfLooseIsoDiTaus_ = looseIsoTaus.size(); 
  if(looseIsoTaus.size()>0 ){
    identifiedTaus.swap(looseIsoTaus);
    if(verbose_) cout << identifiedTaus.size() << "  isolated taus found..." << endl;
    for(unsigned int i=0; i<identifiedTaus.size(); i++){
      if(verbose_) cout << "Testing isolation of " << i << "th tau" << endl;
      const pat::Tau*  tau_i = dynamic_cast<const pat::Tau*>(  ((*diTaus)[ identifiedTaus[i] ].leg2()).get() );
      double sumIsoTau_i = 0.;
      sumIsoTau_i += tau_i->isolationPFChargedHadrCandsPtSum();
      sumIsoTau_i += std::max( (tau_i->isolationPFGammaCandsEtSum() -
				rhoFastJet_*TMath::Pi()*0.5*0.5), 0.0);
      //sumIsoTau_i += tau_i->isolationPFGammaCandsEtSum();
      //sumIsoTau_i += tau_i->isolationPFNeutrHadrCandsEtSum();
      if(sumIsoTau_i<sumIsoTau){
	index = identifiedTaus[i];
	sumIsoTau = sumIsoTau_i;
      } 
    }
  } 
  else if(identifiedTaus.size()>0  ) {
    index = identifiedTaus[tRandom_->Integer( identifiedTaus.size() )];
    if(verbose_) cout << "Random selection has chosen index " << index << endl;   
  }


  if(verbose_) cout << "Chosen index " << index << endl;
  identifiedTaus.clear(); /*not_identifiedTaus.clear(); */ looseIsoTaus.clear();

  theDiTau = &(*diTaus)[index];

  diTauCharge_ = theDiTau->charge();
  METP4_->push_back((*rawMet)[0].p4());
  METP4_->push_back((*met)[0].p4());
  if(isMC_) genMETP4_->push_back( (*rawMet)[0].genMET()->p4() );
  sumEt_  = (*met)[0].sumEt();
  MtLeg1_ =  theDiTau->mt1MET();
  pZeta_     =  theDiTau->pZeta();
  pZetaVis_  =  theDiTau->pZetaVis();
  pZetaSig_  =  theDiTau->pZetaSig();
  mTauTauMin_=  theDiTau->mTauTauMin();

  isElecLegMatched_  = 0;
  isTauLegMatched_   = 0;

  const pat::Electron* leg1 = dynamic_cast<const pat::Electron*>( (theDiTau->leg1()).get() );
  const pat::Tau*      leg2 = dynamic_cast<const pat::Tau*>(  (theDiTau->leg2()).get() );


  //mvaAntiE_ = antiE_->MVAValue( leg2 );
  //if(verbose_)
  //cout << "//Event " << iEvent.run() << ":" << (iEvent.eventAuxiliary()).event() << ":" << iEvent.luminosityBlock()
  // << antiE_->MVAValue( leg2 ) 
  // << endl;

  vector<string> triggerPaths;
  vector<string> XtriggerPaths;
  vector<string> HLTfiltersElec;
  vector<string> HLTfiltersTau;

  if(isMC_){

    // X-triggers
    XtriggerPaths.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v*");
    XtriggerPaths.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*");

    /*
    // for Summer11
    triggerPaths.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2");
    triggerPaths.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2");
 
    HLTfiltersElec.push_back("hltEle15CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter");
    HLTfiltersElec.push_back("hltEle20CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter");
    HLTfiltersTau.push_back("hltOverlapFilterIsoEle15IsoPFTau15");
    HLTfiltersTau.push_back("hltOverlapFilterIsoEle15IsoPFTau20");
    */

    // for Fall11
    triggerPaths.push_back("HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v1");
    triggerPaths.push_back("HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v1");

    HLTfiltersElec.push_back("hltEle18CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter");
    HLTfiltersElec.push_back("hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter");
    HLTfiltersTau.push_back("hltOverlapFilterIsoEle18MediumIsoPFTau20");
    HLTfiltersTau.push_back("hltOverlapFilterIsoEle20MediumIsoPFTau20");
  }
  else{
    
    // X-triggers
    XtriggerPaths.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v*");
    XtriggerPaths.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*");
    XtriggerPaths.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_v*");
    XtriggerPaths.push_back("HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*");

    // Single Electron triggers + X-triggers
    triggerPaths.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v1");
    triggerPaths.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2");
    triggerPaths.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v4");
    triggerPaths.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v6");
    triggerPaths.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v8");
    triggerPaths.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v9");
    triggerPaths.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_v2");
    triggerPaths.push_back("HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v1");
    triggerPaths.push_back("HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v1");
    triggerPaths.push_back("HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v5");
    triggerPaths.push_back("HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v6");
                              
    // watch out! name convention changed with time
    HLTfiltersElec.push_back("hltEle15CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter");
    HLTfiltersElec.push_back("hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter");
    // watch out! name convention changed with time
    HLTfiltersElec.push_back("hltEle18CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter");
    HLTfiltersElec.push_back("hltEle18CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter");
 
    HLTfiltersElec.push_back("hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1SingleEG18orL1SingleEG20");
 
   
    //HLTfiltersTau.push_back("hltPFTau15TrackLooseIso");
    //HLTfiltersTau.push_back("hltPFTau20TrackLooseIso");
    HLTfiltersTau.push_back("hltOverlapFilterIsoEle15IsoPFTau15");
    HLTfiltersTau.push_back("hltOverlapFilterIsoEle15IsoPFTau20");
    HLTfiltersTau.push_back("hltOverlapFilterIsoEle15TightIsoPFTau20");
    HLTfiltersTau.push_back("hltOverlapFilterIsoEle18MediumIsoPFTau20");
    HLTfiltersTau.push_back("hltOverlapFilterIsoEle20MediumIsoPFTau20");

  }

  for(unsigned int i=0;i<triggerPaths.size();i++){
    if(!trigger){
      continue;
      cout << "Invalid triggerEvent !" << endl;
    }
    const pat::TriggerPath *triggerPath =  trigger->path(triggerPaths[i]);

    if(verbose_){
      cout<<  "Testing " << triggerPaths[i] << endl;
      if(triggerPath) cout << "Is there..." << endl;
      if(triggerPath && triggerPath->wasRun()) cout << "Was run..." << endl;
      if(triggerPath && triggerPath->wasRun() && triggerPath->wasAccept()) cout << "Was accepted..." << endl;
    }
    
    if(triggerPath && triggerPath->wasRun() && 
       triggerPath->wasAccept() && 
       triggerPath->prescale()==1) triggerBits_->push_back(1);
    else if (triggerPath && triggerPath->wasRun() && 
	     triggerPath->wasAccept() && 
	     triggerPath->prescale()!=1) triggerBits_->push_back(2);
    else triggerBits_->push_back(0);
  }

  
  for(unsigned int i=0 ; i< HLTfiltersElec.size() ; i++){
    bool matched = false;
    for(pat::TriggerObjectStandAloneCollection::const_iterator it = triggerObjs->begin() ; it !=triggerObjs->end() ; it++){
      pat::TriggerObjectStandAlone *aObj = const_cast<pat::TriggerObjectStandAlone*>(&(*it));
      if(verbose_) {
	if( Geom::deltaR( aObj->triggerObject().p4(), leg1->p4() )<0.3 ){
	  for(unsigned int k =0; k < (aObj->filterLabels()).size() ; k++){
	    cout << "Object passing " << (aObj->filterLabels())[k] << " within 0.3 of electron" << endl;
	  }
	}
      }
      if( Geom::deltaR( aObj->triggerObject().p4(), leg1->p4() )<0.3  && aObj->hasFilterLabel(HLTfiltersElec[i]) ){
	matched = true;
      }
    }
    if(matched) 
      tauXTriggers_->push_back(1);
    else 
      tauXTriggers_->push_back(0);
    if(verbose_){
      if(matched) cout << "Electron matched within dR=0.3 with trigger object passing filter " << HLTfiltersElec[i] << endl;
      else cout << "!!! Electron is not trigger matched within dR=0.3 !!!" << endl;
    }
  }
  for(unsigned int i=0 ; i< HLTfiltersTau.size() ; i++){
    bool matched = false;
    for(pat::TriggerObjectStandAloneCollection::const_iterator it = triggerObjs->begin() ; it !=triggerObjs->end() ; it++){
      pat::TriggerObjectStandAlone *aObj = const_cast<pat::TriggerObjectStandAlone*>(&(*it));
      if(verbose_) {
        if( Geom::deltaR( aObj->triggerObject().p4(), leg2->p4() )<0.3 ){
          for(unsigned int k =0; k < (aObj->filterLabels()).size() ; k++){
            cout << "Object passing " << (aObj->filterLabels())[k] << " within 0.3 of tau" << endl;
          }
        }
      }
      if( Geom::deltaR( aObj->triggerObject().p4(), leg2->p4() )<0.3  && aObj->hasFilterLabel(HLTfiltersTau[i]) ){
	matched = true;
      }
    }
    if(matched) 
      tauXTriggers_->push_back(1);
    else 
      tauXTriggers_->push_back(0);
    if(verbose_){
      if(matched) cout << "Tau matched within dR=0.3 with trigger object passing filter " << HLTfiltersTau[i] << endl;
      else cout << "!!! Tau is not trigger matched within dR=0.3 !!!" << endl;
    }
  }



  // triggers Elec
  if(verbose_){
    cout << "Electron triggers" << endl;
    const pat::TriggerObjectStandAloneCollection trColl = leg1->triggerObjectMatchesByType(trigger::TriggerElectron);
    for(pat::TriggerObjectStandAloneCollection::const_iterator it = trColl.begin(); it != trColl.end(); it++){
      for(unsigned int k = 0; k < (it->pathNames(false,false)).size(); k++){
	cout << (it->pathNames(false,false))[k] << " with filter: " << endl;
	for(unsigned int m = 0; m < (it->filterLabels()).size(); m++){
	  cout << "  - " << (it->filterLabels())[m] << endl;
	}
      }
    }
  }
  // triggers Tau
  if(verbose_){
    const pat::TriggerObjectStandAloneCollection trColl = leg2->triggerObjectMatchesByType(trigger::TriggerTau);
    cout << "Tau triggers" << endl;
    for(pat::TriggerObjectStandAloneCollection::const_iterator it = trColl.begin(); it != trColl.end(); it++){
      for(unsigned int k = 0; k < (it->pathNames(false,false)).size(); k++){
	cout << (it->pathNames(false,false))[k] << " with filter: " << endl;
	for(unsigned int m = 0; m < (it->filterLabels()).size(); m++){
	  cout << "  - " << (it->filterLabels())[m] << endl;
	}
      }
    }
  }
  
  diTauLegsP4_->push_back(leg1->p4());
  diTauLegsP4_->push_back(leg2->p4());
  
  if(isMC_){
    if( (leg1->genParticleById(11,0,true)).isNonnull()  ){
      genDiTauLegsP4_->push_back( leg1->genParticleById(11,0,true)->p4() );
      isElecLegMatched_ = 1;
    }
    else{
      genDiTauLegsP4_->push_back( math::XYZTLorentzVectorD(0,0,0,0) );
      if(verbose_){
	for(unsigned int l = 0; l < leg1->genParticlesSize() ; l++){
	  if((leg1->genParticleRefs())[l]->pt() < 0.5 ) continue;
	  cout << "Elec leg matchged to particle " << (leg1->genParticleRefs())[l]->pdgId() 
	       << " with pt " << (leg1->genParticleRefs())[l]->pt()
	       << endl;
	}
      }
    }

    //Jet matching electron
    if(verbose_)cout <<"particles matching electron :"<<endl;
    leg1GenDecay_ = -99;
    ptHighestLeg1Matched_ = -99;
    ratioPtMaxMatchPtGenElec_ = -99;
    float GenElecPt = -99;
    int c=0;
    for(unsigned int k = 0; k < genParticles->size(); k ++){
      if(!(Geom::deltaR( (*genParticles)[k].p4(),leg1->p4())<0.3 && (*genParticles)[k].pt()>5))continue;
      if (verbose_){
	cout <<"    "<<c<<"th particle :"<<(*genParticles)[k].pdgId()<<endl;
	cout <<"    "<<c<<"th particle status:"<<(*genParticles)[k].status()<<endl;
	if((*genParticles)[k].mother()!=0)cout <<"    "<<c<<"th particle mother:"<<(*genParticles)[k].mother()->pdgId()<<endl;
	}
      c++;
      if(abs((*genParticles)[k].pdgId())==4 && (*genParticles)[k].status()==2){
	leg1GenDecay_ = (*genParticles)[k].pdgId();
      }
      if(abs((*genParticles)[k].pdgId())==5  && (*genParticles)[k].status()==2){
	leg1GenDecay_ = (*genParticles)[k].pdgId();
      }
      if(abs((*genParticles)[k].pdgId())==11 ){
	GenElecPt = (*genParticles)[k].pt();
	leg1GenDecay_ *= (*genParticles)[k].pdgId();
      }
      if((*genParticles)[k].pt()>ptHighestLeg1Matched_ && (*genParticles)[k].status() == 1)ptHighestLeg1Matched_ = (*genParticles)[k].pt();
    }
  

    ratioPtMaxMatchPtGenElec_ = ptHighestLeg1Matched_/GenElecPt;
    if(verbose_){
    cout<<"leg1GenDecay : "<<leg1GenDecay_<<endl;
    cout<<"Pt of  the electron : "<<GenElecPt<<endl;
    cout<<"Pt of highest matched to the electron particle : "<<ptHighestLeg1Matched_<<endl;
    cout<<"ratio of gen pt : "<<ratioPtMaxMatchPtGenElec_<<endl;
    }
  

    bool leg2IsFromElec = false;
    math::XYZTLorentzVectorD genElecP4(0,0,0,0);
    for(unsigned int k = 0; k < genParticles->size(); k ++){
      if( abs((*genParticles)[k].pdgId()) != 11  || (*genParticles)[k].status()!=3 )
        continue;
      if(Geom::deltaR( (*genParticles)[k].p4(),leg2->p4())<0.15){
	leg2IsFromElec = true;
	genElecP4 = (*genParticles)[k].p4();
      }
    }

    if( leg2->genJet() !=0 ) 
      genDiTauLegsP4_->push_back(leg2->genJet()->p4());
    else if(leg2IsFromElec)
      genDiTauLegsP4_->push_back( genElecP4 );      
    else{
      genDiTauLegsP4_->push_back( math::XYZTLorentzVectorD(0,0,0,0) );
      if(verbose_) cout << "WARNING: no genJet matched to the leg2 with eta,phi " << leg2->eta() << ", " << leg2->phi() << endl;
    }

    bool tauHadMatched = false;
    for(unsigned int k = 0; k < tauGenJets->size(); k++){
      if( Geom::deltaR( (*tauGenJets)[k].p4(),leg2->p4() ) < 0.15 ) tauHadMatched = true;
    }

    if( ( (leg2->genParticleById(15,0,true)).isNonnull() || (leg2->genParticleById(-15,0,true)).isNonnull() ) && tauHadMatched ) {
      isTauLegMatched_ = 1;
    }
    else if(verbose_){
      for(unsigned int l = 0; l < leg2->genParticlesSize() ; l++){
	if((leg2->genParticleRefs())[l]->pt() < 0.5 ) continue;
	cout << "Tau leg matchged to particle " << (leg2->genParticleRefs())[l]->pdgId() 
	     << " with pt " << (leg2->genParticleRefs())[l]->pt()
	     << endl;
      }
    }
  }

  if((leg2->signalPFChargedHadrCands()).size()==1 && (leg2->signalPFGammaCands()).size()==0) decayMode_ = 0; 
  else if((leg2->signalPFChargedHadrCands()).size()==1 && (leg2->signalPFGammaCands()).size()>0)  decayMode_ = 1; 
  else if((leg2->signalPFChargedHadrCands()).size()==3) decayMode_ = 2; 
  else  decayMode_ = -99;

  for(unsigned int k = 0 ; k < (leg2->signalPFGammaCands()).size() ; k++){
    reco::PFCandidateRef gamma = (leg2->signalPFGammaCands()).at(k);
    if( (leg2->leadPFChargedHadrCand()).isNonnull() ){
      gammadR_->push_back(   Geom::deltaR( gamma->p4(), leg2->leadPFChargedHadrCand()->p4() ) );
      gammadPhi_->push_back( Geom::deltaPhi( gamma->p4(), leg2->leadPFChargedHadrCand()->p4() ) );
      gammadEta_->push_back( gamma->eta() - leg2->leadPFChargedHadrCand()->eta() ) ;
    }
    else{
      gammadR_->push_back(   Geom::deltaR( gamma->p4(), leg2->p4() ) );
      gammadPhi_->push_back( Geom::deltaPhi( gamma->p4(), leg2->p4() ) );
      gammadEta_->push_back(  gamma->eta() - leg2->eta() ) ;
    }
    gammaPt_->push_back(  gamma->pt() );
  }


  visibleTauMass_ = leg2->mass();
 
  if((leg2->leadPFChargedHadrCand()->trackRef()).isNonnull()){
    leadPFChargedHadrTrackPt_ = leg2->leadPFChargedHadrCand()->trackRef()->pt();
    leadPFChargedHadrTrackP_  = leg2->leadPFChargedHadrCand()->trackRef()->p();
  } else if((leg2->leadPFChargedHadrCand()->gsfTrackRef()).isNonnull()){
    leadPFChargedHadrTrackPt_ = leg2->leadPFChargedHadrCand()->gsfTrackRef()->pt();
    leadPFChargedHadrTrackP_  = leg2->leadPFChargedHadrCand()->gsfTrackRef()->p();
  } else{
    leadPFChargedHadrTrackPt_ = -99;
    leadPFChargedHadrTrackP_  = -99;
  }
  leadPFChargedHadrPt_  = leg2->leadPFChargedHadrCand()->pt();  
  leadPFChargedHadrP_   = leg2->leadPFChargedHadrCand()->p();  

  leadPFChargedHadrMva_        =   leg2->electronPreIDOutput() ;	
  leadPFChargedHadrHcalEnergy_ =  (leg2->leadPFChargedHadrCand()).isNonnull() ? (leg2->leadPFChargedHadrCand())->hcalEnergy() : -99 ;
  leadPFChargedHadrEcalEnergy_ =  (leg2->leadPFChargedHadrCand()).isNonnull() ? (leg2->leadPFChargedHadrCand())->ecalEnergy() : -99;
  
  if( (leg2->leadPFCand()).isNonnull() ){
    leadPFCandMva_               =  (leg2->leadPFCand())->mva_e_pi() ;	
    leadPFCandHcalEnergy_        =  (leg2->leadPFCand())->hcalEnergy() ;
    leadPFCandEcalEnergy_        =  (leg2->leadPFCand())->ecalEnergy() ;
    leadPFCandPt_                =  (leg2->leadPFCand())->pt();
    leadPFCandP_                 =  (leg2->leadPFCand())->p();
  }
  
  signalPFChargedHadrCands_ = (leg2->signalPFChargedHadrCands()).size();
  signalPFGammaCands_       = (leg2->signalPFGammaCands()).size();
  emFraction_ = leg2->emFraction();
  hasGsf_     = ((leg2->leadPFChargedHadrCand())->gsfTrackRef()).isNonnull() ? 1. : 0.;
  

  tightestHPSWP_ = -1;
  if(leg2->tauID("byVLooseIsolation")>0.5) tightestHPSWP_++;
  if(leg2->tauID("byLooseIsolation")>0.5)  tightestHPSWP_++;
  if(leg2->tauID("byMediumIsolation")>0.5) tightestHPSWP_++;
  if(leg2->tauID("byTightIsolation")>0.5)  tightestHPSWP_++;
  tightestHPSDBWP_ = -1;
  if(leg2->tauID("byVLooseCombinedIsolationDeltaBetaCorr")>0.5) tightestHPSDBWP_++;
  if(leg2->tauID("byLooseCombinedIsolationDeltaBetaCorr")>0.5)  tightestHPSDBWP_++;
  if(leg2->tauID("byMediumCombinedIsolationDeltaBetaCorr")>0.5) tightestHPSDBWP_++;
  if(leg2->tauID("byTightCombinedIsolationDeltaBetaCorr")>0.5)  tightestHPSDBWP_++;


  dxy1_ = vertexes->size()!=0 ? leg1->gsfTrack()->dxy( (*vertexes)[0].position() ) : -99;
  dz1_ = vertexes->size()!=0 ? leg1->gsfTrack()->dz( (*vertexes)[0].position() ) : -99;
  dxy2_ = -99;
  dz2_ = -99;
  hasKft_ = 0;
  if( vertexes->size()!=0 && (leg2->leadPFChargedHadrCand()).isNonnull() ){
    if( (leg2->leadPFChargedHadrCand()->trackRef()).isNonnull() ){
      dxy2_ = leg2->leadPFChargedHadrCand()->trackRef()->dxy( (*vertexes)[0].position() );
      dz2_ = leg2->leadPFChargedHadrCand()->trackRef()->dz( (*vertexes)[0].position() );
      hasKft_ = 1;
    }
    if( (leg2->leadPFChargedHadrCand()->gsfTrackRef()).isNonnull() ){
      dxy2_ = leg2->leadPFChargedHadrCand()->gsfTrackRef()->dxy( (*vertexes)[0].position() );
      dz2_ = leg2->leadPFChargedHadrCand()->gsfTrackRef()->dz( (*vertexes)[0].position() );
    }
  }
  
  // isoDeposit definition: 2010
  isodeposit::AbsVetos vetos2010ChargedLeg1; 
  isodeposit::AbsVetos vetos2010NeutralLeg1; 
  isodeposit::AbsVetos vetos2010PhotonLeg1;
  // isoDeposit definition: 2011
  isodeposit::AbsVetos vetos2011ChargedLeg1; 
  isodeposit::AbsVetos vetos2011NeutralLeg1; 
  isodeposit::AbsVetos vetos2011PhotonLeg1;
  // isoDeposit definition: 2011 for EE
  isodeposit::AbsVetos vetos2011EEChargedLeg1; 
  isodeposit::AbsVetos vetos2011EENeutralLeg1;  
  isodeposit::AbsVetos vetos2011EEPhotonLeg1;

 
  vetos2010ChargedLeg1.push_back(new isodeposit::ConeVeto(reco::isodeposit::Direction(leg1->eta(),leg1->phi()),0.01));
  vetos2010ChargedLeg1.push_back(new isodeposit::ThresholdVeto(0.5));
  vetos2010NeutralLeg1.push_back(new isodeposit::ConeVeto(isodeposit::Direction(leg1->eta(),leg1->phi()),0.08));
  vetos2010NeutralLeg1.push_back(new isodeposit::ThresholdVeto(1.0));
  vetos2010PhotonLeg1.push_back( new isodeposit::ConeVeto(isodeposit::Direction(leg1->eta(),leg1->phi()),0.05));
  vetos2010PhotonLeg1.push_back( new isodeposit::ThresholdVeto(1.0));

  vetos2011ChargedLeg1.push_back(new isodeposit::ConeVeto(reco::isodeposit::Direction(leg1->eta(),leg1->phi()),0.01));
  vetos2011ChargedLeg1.push_back(new isodeposit::ThresholdVeto(0.0));
  vetos2011NeutralLeg1.push_back(new isodeposit::ConeVeto(isodeposit::Direction(leg1->eta(),leg1->phi()),0.01));
  vetos2011NeutralLeg1.push_back(new isodeposit::ThresholdVeto(0.5));
  vetos2011PhotonLeg1.push_back( new isodeposit::ConeVeto(isodeposit::Direction(leg1->eta(),leg1->phi()),0.01));
  vetos2011PhotonLeg1.push_back( new isodeposit::ThresholdVeto(0.5));

  vetos2011EEChargedLeg1.push_back(new isodeposit::ConeVeto(reco::isodeposit::Direction(leg1->eta(),leg1->phi()),0.015));
  vetos2011EEChargedLeg1.push_back(new isodeposit::ThresholdVeto(0.0));
  vetos2011EENeutralLeg1.push_back(new isodeposit::ConeVeto(reco::isodeposit::Direction(leg1->eta(),leg1->phi()),0.01));
  vetos2011EENeutralLeg1.push_back(new isodeposit::ThresholdVeto(0.5));
  vetos2011EEPhotonLeg1.push_back(new isodeposit::ConeVeto(reco::isodeposit::Direction(leg1->eta(),leg1->phi()),0.08));
  vetos2011EEPhotonLeg1.push_back(new isodeposit::ThresholdVeto(0.5));
  

  chIsoLeg1v1_   = 
    leg1->isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4,vetos2010ChargedLeg1).first;
  nhIsoLeg1v1_   = 
    leg1->isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4,vetos2010NeutralLeg1).first;
  phIsoLeg1v1_   = 
    leg1->isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4,vetos2010PhotonLeg1).first;
  elecIsoLeg1v1_ = 
    leg1->isoDeposit(pat::User3Iso)->depositAndCountWithin(0.4,vetos2010ChargedLeg1).first;
  muIsoLeg1v1_   = 
    leg1->isoDeposit(pat::User2Iso)->depositAndCountWithin(0.4,vetos2010ChargedLeg1).first;
  chIsoPULeg1v1_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4,vetos2010ChargedLeg1).first;
  nhIsoPULeg1v1_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4,vetos2010NeutralLeg1).first;
  phIsoPULeg1v1_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4,vetos2010PhotonLeg1).first;

  chIsoLeg1v2_   = 
    leg1->isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4,vetos2011ChargedLeg1).first;
  nhIsoLeg1v2_   = 
    leg1->isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4,vetos2011NeutralLeg1).first;
  phIsoLeg1v2_   = 
    leg1->isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4,vetos2011PhotonLeg1).first;
  elecIsoLeg1v2_ = 
    leg1->isoDeposit(pat::User3Iso)->depositAndCountWithin(0.4,vetos2011ChargedLeg1).first;
  muIsoLeg1v2_   = 
    leg1->isoDeposit(pat::User2Iso)->depositAndCountWithin(0.4,vetos2011ChargedLeg1).first;
  chIsoPULeg1v2_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4,vetos2011ChargedLeg1).first;
  nhIsoPULeg1v2_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4,vetos2011NeutralLeg1).first;
  phIsoPULeg1v2_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4,vetos2011PhotonLeg1).first;
  
  chIsoEELeg1v2_   = 
    leg1->isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, vetos2011EEChargedLeg1).first;
  nhIsoEELeg1v2_   = 
    leg1->isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, vetos2011EENeutralLeg1).first;
  phIsoEELeg1v2_   = 
    leg1->isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, vetos2011EEPhotonLeg1).first;
  nhIsoEEPULeg1v2_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2011EENeutralLeg1).first;
  phIsoEEPULeg1v2_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2011EEPhotonLeg1).first;
  elecIsoEELeg1v2_     = 
    leg1->isoDeposit(pat::User3Iso)->depositAndCountWithin(0.4,vetos2011EEChargedLeg1).first;
  muIsoEELeg1v2_       = 
    leg1->isoDeposit(pat::User2Iso)->depositAndCountWithin(0.4,vetos2011EEChargedLeg1).first;



  // loop over pfElectrons to make sure we don't have low-mass DY events
  for(unsigned int i = 0 ; i < pfCandidates->size() ; i++){
    const reco::PFCandidate cand = (*pfCandidates)[i];
    if(! (cand.particleId()== PFCandidate::e && cand.pt()>0.5) ) continue;
    float dz = 999;
    if(cand.trackRef().isNonnull()) 
      dz = (cand.trackRef())->dz(    leg1->vertex() );
    if(cand.gsfTrackRef().isNonnull())
      dz = (cand.gsfTrackRef())->dz( leg1->vertex() );
    if(dz>0.2) continue;
    pfElectrons_->push_back(cand.p4());
  }

  chIsoLeg2_ = leg2->isolationPFChargedHadrCandsPtSum();
  nhIsoLeg2_ = 0.;
  for(unsigned int i = 0; i < (leg2->isolationPFNeutrHadrCands()).size(); i++){
    nhIsoLeg2_ += (leg2->isolationPFNeutrHadrCands()).at(i)->pt();
  }
  phIsoLeg2_ = leg2->isolationPFGammaCandsEtSum();

  // cleaning
  for(unsigned int i = 0; i <vetos2010ChargedLeg1.size(); i++){
    delete vetos2010ChargedLeg1[i];
  }
  for(unsigned int i = 0; i <vetos2010NeutralLeg1.size(); i++){
    delete vetos2010NeutralLeg1[i];
    delete vetos2010PhotonLeg1[i];
  }
  for(unsigned int i = 0; i <vetos2011ChargedLeg1.size(); i++){
    delete vetos2011ChargedLeg1[i];
  }
  for(unsigned int i = 0; i <vetos2011NeutralLeg1.size(); i++){
    delete vetos2011NeutralLeg1[i];
    delete vetos2011PhotonLeg1[i];
  }
  for(unsigned int i = 0; i <vetos2011EEChargedLeg1.size(); i++){
    delete vetos2011EEChargedLeg1[i];
  }
  for(unsigned int i = 0; i <vetos2011EENeutralLeg1.size(); i++){
    delete vetos2011EENeutralLeg1[i];
    delete vetos2011EEPhotonLeg1[i];
  }
  

  //
  nBrehm_= leg1->numberOfBrems();
  likelihood_ = -99; //leg1->electronID("electronIDLH");
  nHits_ = leg1->userFloat("nHits");
  sihih_ = leg1->userFloat("sihih");
  dPhi_  = leg1->userFloat("dPhi");
  dEta_  = leg1->userFloat("dEta");
  HoE_   = leg1->hadronicOverEm();
  EoP_   = leg1->eSuperClusterOverP();
  fbrem_ = leg1->fbrem();
  //cout<<"anticonv :"<<leg1->userInt("antiConv")<<endl;
  antiConv_ = leg1->userInt("antiConv");


  tightestCutBasedWP_ = -1;
  if((leg1->userFloat("nHits")<=1
      && (  
	  (leg1->isEB() && leg1->userFloat("sihih")<0.010 && leg1->userFloat("dPhi")<0.80 &&  
	   leg1->userFloat("dEta") <0.007 && leg1->userFloat("HoE") <0.15)    
	  ||    
	  (leg1->isEE() && leg1->userFloat("sihih")<0.030 && leg1->userFloat("dPhi")<0.70 &&  
	   leg1->userFloat("dEta") <0.010 && leg1->userFloat("HoE") <0.07)    
	  )
      )
     ) tightestCutBasedWP_ ++;
  if((leg1->userFloat("nHits")==0 && leg1->userInt("antiConv")>0.5 	
      && (							
	  (leg1->pt()>=20 && (    						
			      (leg1->isEB() && leg1->userFloat("sihih")<0.010 && leg1->userFloat("dPhi")<0.06 &&    
			       leg1->userFloat("dEta")< 0.004 && leg1->userFloat("HoE") <0.04)      
			      ||						
			      (leg1->isEE() && leg1->userFloat("sihih")<0.030 && leg1->userFloat("dPhi")<0.030 &&   
			       leg1->userFloat("dEta") <0.007 && leg1->userFloat("HoE") <0.025) ))  
          || 							
	  (leg1->pt()<20 && (leg1->fbrem()>0.15 || (fabs(leg1->superClusterPosition().Eta())<1. && leg1->eSuperClusterOverP()>0.95) ) && (  
																	  (leg1->isEB() && leg1->userFloat("sihih")<0.010 && leg1->userFloat("dPhi")<0.030 &&   
																	   leg1->userFloat("dEta") <0.004 && leg1->userFloat("HoE") <0.025)     
																	  ||						
																	  (leg1->isEE() && leg1->userFloat("sihih")<0.030 && leg1->userFloat("dPhi")<0.020 &&   
																	   leg1->userFloat("dEta") <0.005 && leg1->userFloat("HoE") <0.025) ))  
	  )  						
      )     
     ) tightestCutBasedWP_ ++;
  
  tightestMVAWP_ = -1;
  bool passesMVA = 
    (leg1->pt()<=20 && fabs(leg1->superClusterPosition().Eta())>=0.0 && fabs(leg1->superClusterPosition().Eta())<1.0 && leg1->userFloat("mva")>0.133) ||
    (leg1->pt()<=20 && fabs(leg1->superClusterPosition().Eta())>=1.0 && fabs(leg1->superClusterPosition().Eta())<1.5 && leg1->userFloat("mva")>0.465) ||
    (leg1->pt()<=20 && fabs(leg1->superClusterPosition().Eta())>=1.5 && fabs(leg1->superClusterPosition().Eta())<2.5 && leg1->userFloat("mva")>0.518) ||
    (leg1->pt()>20  && fabs(leg1->superClusterPosition().Eta())>=0.0 && fabs(leg1->superClusterPosition().Eta())<1.0 && leg1->userFloat("mva")>0.942) ||
    (leg1->pt()>20  && fabs(leg1->superClusterPosition().Eta())>=1.0 && fabs(leg1->superClusterPosition().Eta())<1.5 && leg1->userFloat("mva")>0.947) ||
    (leg1->pt()>20  && fabs(leg1->superClusterPosition().Eta())>=1.5 && fabs(leg1->superClusterPosition().Eta())<2.5 && leg1->userFloat("mva")>0.878);
  bool mvaPreselectionNoIso = 
    leg1->userInt("antiConv")>0.5 && leg1->userFloat("nHits")<=0 &&
    ((leg1->isEB()
      && leg1->userFloat("sihih") < 0.01 
      && TMath::Abs(leg1->userFloat("dEta")) < 0.007
      && TMath::Abs(leg1->userFloat("dPhi")) < 0.15
      && leg1->userFloat("HoE") < 0.12
      ) ||
     (leg1->isEE() 
      && leg1->userFloat("sihih") < 0.03 
      && TMath::Abs(leg1->userFloat("dEta")) < 0.009
      && TMath::Abs(leg1->userFloat("dPhi")) < 0.10
      && leg1->userFloat("HoE") < 0.10
      ));
  bool mvaPreselection = leg1->userInt("mvaPreselection")>0.5 ;

  if( passesMVA )                        tightestMVAWP_ = 0;
  if( passesMVA && mvaPreselectionNoIso) tightestMVAWP_ = 1;
  if( passesMVA && mvaPreselection)      tightestMVAWP_ = 2;
    
  if( (int(leg1->electronID("eidCiCHZZTight"))&1)==1     )    tightestCiCWP_++;
  if( (int(leg1->electronID("eidCiCHZZSuperTight"))&1)==1)    tightestCiCWP_++;
  if( (int(leg1->electronID("eidCiCHZZHyperTight1"))&1)==1  ) tightestCiCWP_++;

  tightestDanieleMVAWP_ = leg1->userFloat("mvaDaniele");

  tightestAntiECutsWP_ = -1;
  if( leg2->tauID("againstElectronTight")>0.5 ) tightestAntiECutsWP_++;
  tightestAntiEMVAWP_ = -1;
  //std::cout<<leg2->pt()<<"   "<<leg2->eta()<<"   "<<leg2->phi()<<std::endl;
  //cout<<"againstElectronMVA  : "<<leg2->tauID("againstElectronMVA")<<endl;
  //cout<<"againstElectronCuts  : "<<leg2->tauID("againstElectronTight")<<endl;
  if( leg2->tauID("againstElectronMVA")>0.5 ) tightestAntiEMVAWP_++;


  diTauVisP4_->push_back( theDiTau->p4Vis() );
  diTauCAP4_->push_back( theDiTau->p4CollinearApprox() );
  diTauICAP4_->push_back( theDiTau->p4ImprovedCollinearApprox() );

  math::XYZTLorentzVectorD nSVfitFitP4(0,0,0,(theDiTau->p4Vis()).M() );
  if( theDiTau->hasNSVFitSolutions() && theDiTau->nSVfitSolution("psKine_MEt_logM_fit",0)!=0 /*&& theDiTau->nSVfitSolution("psKine_MEt_logM_fit",0)->isValidSolution()*/ ){
    if(verbose_) cout << "Visible mass ==> " << nSVfitFitP4.E() << endl;
    nSVfitFitP4.SetPxPyPzE( 0,0,0, theDiTau->nSVfitSolution("psKine_MEt_logM_fit",0)->mass() ) ;
    if(verbose_) cout << "SVFit fit solution ==> " << nSVfitFitP4.E() << endl;
  }
  diTauSVfitP4_->push_back( nSVfitFitP4  );

  int errFlag = 0;
  diTauNSVfitMass_        = (theDiTau->hasNSVFitSolutions() && theDiTau->nSVfitSolution("psKine_MEt_logM_int",&errFlag)!=0 && theDiTau->nSVfitSolution("psKine_MEt_logM_int",0)->isValidSolution() ) 
    ? theDiTau->nSVfitSolution("psKine_MEt_logM_int",0)->mass()        : -99; 
  diTauNSVfitMassErrUp_   = (theDiTau->hasNSVFitSolutions() && theDiTau->nSVfitSolution("psKine_MEt_logM_int",&errFlag)!=0 && theDiTau->nSVfitSolution("psKine_MEt_logM_int",0)->isValidSolution() ) 
    ? theDiTau->nSVfitSolution("psKine_MEt_logM_int",0)->massErrUp()   : -99; 
  diTauNSVfitMassErrDown_ = (theDiTau->hasNSVFitSolutions() && theDiTau->nSVfitSolution("psKine_MEt_logM_int",&errFlag)!=0 && theDiTau->nSVfitSolution("psKine_MEt_logM_int",0)->isValidSolution() ) 
    ? theDiTau->nSVfitSolution("psKine_MEt_logM_int",0)->massErrDown() : -99; 

  run_   = iEvent.run();
  event_ = (iEvent.eventAuxiliary()).event();
  lumi_  = iEvent.luminosityBlock();
  
  std::map<double, math::XYZTLorentzVectorD ,ElecTauStreamAnalyzer::more> sortedJets;
  std::map<double, math::XYZTLorentzVectorD ,ElecTauStreamAnalyzer::more> sortedJetsIDL1Offset;
  std::map<double, math::XYZTLorentzVectorD ,ElecTauStreamAnalyzer::more> sortedJetsID;
  std::map<double, math::XYZTLorentzVectorD ,ElecTauStreamAnalyzer::more> sortedJetsIDUp;
  std::map<double, math::XYZTLorentzVectorD ,ElecTauStreamAnalyzer::more> sortedJetsIDDown;
  std::map<double, math::XYZTLorentzVectorD ,ElecTauStreamAnalyzer::more> sortedGenJetsID;
  std::map<double, std::pair<float,float> ,  ElecTauStreamAnalyzer::more> bTaggers;
  std::map<double, std::pair<float,float> ,  ElecTauStreamAnalyzer::more> jetPVassociation;
  std::map<double, std::pair<float,float> ,  ElecTauStreamAnalyzer::more> jetMoments;


  for(unsigned int it = 0; it < jets->size() ; it++){

    pat::Jet* jet = const_cast<pat::Jet*>(&(*jets)[it]);

    // newJet redone using possibly different JEC ==> it may not contain infos on bTag
    // so I use it together with jet
    pat::Jet* newJet = 0;
    if(newJets!=0) newJetMatched( jet , newJets );
    if(!newJet){
      if(verbose_) cout << "No jet from the new collection can be matched ==> using old one..." << endl;
      newJet = jet;
    }
    
    math::XYZTLorentzVectorD leg2p4 = ( (leg2->pfJetRef()).isNonnull() ) ? leg2->pfJetRef()->p4() : leg2->p4();
    
    if( Geom::deltaR(jet->p4(), leg1->p4())<deltaRLegJet_ || 
	Geom::deltaR(jet->p4(), leg2p4 )<deltaRLegJet_ ){
      if(verbose_) cout << "The jet at (" <<jet->pt()<<","<<jet->eta()<<") is closer than "<<deltaRLegJet_ << " from one of the legs" << endl;  
      continue;
    }
    
    /////////////////////////////////////////////////////////////////////////
    //// use JES uncertainties
    edm::ESHandle<JetCorrectorParametersCollection> jetCorrParameters;
    // get the jet corrector parameters collection from the global tag
    iSetup.get<JetCorrectionsRecord>().get("AK5PF", jetCorrParameters);
    // get the uncertainty parameters from the collection
    JetCorrectorParameters const & param = (*jetCorrParameters)["Uncertainty"];
    // instantiate the jec uncertainty object
    JetCorrectionUncertainty* deltaJEC = new JetCorrectionUncertainty(param);
    deltaJEC->setJetEta(newJet->eta()); deltaJEC->setJetPt(newJet->pt());
    float shift  = deltaJEC->getUncertainty( true );
    /////////////////////////////////////////////////////////////////////////
    
    if(verbose_){
      cout << "Analyzing jet with pt " << (newJet->p4()).Pt() 
	   << "and eta " << (newJet->p4()).Eta() << endl;
      cout << " ==> JEC uncertainty is " << shift*100 << " %" << endl;
      for(unsigned int i = 0; i < (newJet->availableJECSets()).size() ; i++ ){
	std::cout << (newJet->availableJECSets())[i] << std::endl;
      }
      for(unsigned int i = 0; i < (newJet->availableJECLevels()).size() ; i++ ){
	std::cout << (newJet->availableJECLevels())[i] << std::endl;
      }
      std::cout << "L1FastJet ========> " << std::endl;
      std::cout << "Uncorrected " << newJet->correctedJet("Uncorrected","none","patJetCorrFactors").pt() << std::endl;
      std::cout << "L1FastJet "   << newJet->correctedJet("L1FastJet","none",  "patJetCorrFactors").pt() << std::endl;
      std::cout << "L2Relative "  << newJet->correctedJet("L2Relative","none", "patJetCorrFactors").pt() << std::endl; 
      std::cout << "L3Absolute "  << newJet->correctedJet("L3Absolute","none", "patJetCorrFactors").pt() << std::endl; 
      std::cout << "L1Offset ========> " << std::endl;
      std::cout << "Uncorrected " << newJet->jecFactor("Uncorrected","none","patJetCorrFactorsL1Offset")*newJet->pt() << std::endl;
      std::cout << "L1Offset"     << newJet->jecFactor("L1Offset","none",   "patJetCorrFactorsL1Offset")*newJet->pt() << std::endl;
      std::cout << "L2Relative "  << newJet->jecFactor("L2Relative","none", "patJetCorrFactorsL1Offset")*newJet->pt() << std::endl;
      std::cout << "L3Absolute"   << newJet->jecFactor("L3Absolute","none", "patJetCorrFactorsL1Offset")*newJet->pt() << std::endl;  
    }

    std::map<string,float> aMap;

    if( jetID( jet , &((*vertexes)[0]), vtxZ, aMap ) < minJetID_ )  continue;

    sortedJets.insert( make_pair( newJet->correctedJet("Uncorrected").p4().Pt() ,newJet->correctedJet("Uncorrected").p4() ) );
    
    if(newJet->p4().Pt() < minCorrPt_) continue;
    
    // add b-tag info
    bTaggers.insert(         make_pair( newJet->p4().Pt(), make_pair( jet->bDiscriminator("trackCountingHighEffBJetTags"),
								      jet->bDiscriminator("trackCountingHighPurBJetTags")  ) ) );
    // add pu information
    jetPVassociation.insert( make_pair( newJet->p4().Pt(), make_pair(aMap["chFracRawJetE"],
								     aMap["chFracAllChargE"]) ) );
    // add jet moments
    jetMoments.insert(       make_pair( newJet->p4().Pt(), make_pair( jet->etaetaMoment(),
								   jet->phiphiMoment()) ) );
   
    if(isMC_) 
      sortedJetsIDL1Offset.insert( make_pair( newJet->jecFactor("L3Absolute","none", "patJetCorrFactorsL1Offset")*newJet->pt() , 
					      newJet->jecFactor("L3Absolute","none", "patJetCorrFactorsL1Offset")*newJet->p4()) );   
    else 
      sortedJetsIDL1Offset.insert( make_pair( newJet->jecFactor("L2L3Residual","none", "patJetCorrFactorsL1Offset")*newJet->pt() , 
					      newJet->jecFactor("L2L3Residual","none", "patJetCorrFactorsL1Offset")*newJet->p4()) ); 

    if(verbose_) cout << "Components: "
		      << "px=" << (newJet->p4()).Px() << " (" << newJet->px() << "), "
		      << "py=" << (newJet->p4()).Py() << " (" << newJet->py() << "), "
		      << "pz=" << (newJet->p4()).Pz() << " (" << newJet->pz() << "), "
		      << "E="  << (newJet->p4()).E()  << " (" << newJet->energy()  << ")"
		      << endl;

    sortedJetsID.insert(     make_pair( newJet->p4().Pt() ,  newJet->p4() ) );
    sortedJetsIDUp.insert(   make_pair( newJet->p4().Pt() ,  newJet->p4()*(1+shift) ) );
    sortedJetsIDDown.insert( make_pair( newJet->p4().Pt() ,  newJet->p4()*(1-shift) ) );

    if(isMC_){
      if(jet->genJet() != 0) sortedGenJetsID.insert( make_pair( newJet->p4().Pt() ,jet->genJet()->p4() ) );
      else sortedGenJetsID.insert( make_pair( newJet->p4().Pt() , math::XYZTLorentzVectorD(0,0,0,0) ) );
    }
     
  }

  
  for(CImap it = sortedJets.begin(); it != sortedJets.end() ; it++){
    jetsP4_->push_back( it->second );
  }
  for(CImap it = sortedJetsID.begin(); it != sortedJetsID.end() ; it++){
    jetsIDP4_->push_back( it->second );
  }
  for(CImap it = sortedJetsIDUp.begin(); it != sortedJetsIDUp.end() ; it++){
    jetsIDUpP4_->push_back( it->second );
  }
  for(CImap it = sortedJetsIDDown.begin(); it != sortedJetsIDDown.end() ; it++){
    jetsIDDownP4_->push_back( it->second );
  }
  for(CImap it = sortedJetsIDL1Offset.begin(); it != sortedJetsIDL1Offset.end() ; it++){
    jetsIDL1OffsetP4_->push_back( it->second );
  }
  for(CImap it = sortedGenJetsID.begin(); it != sortedGenJetsID.end() ; it++){
    genJetsIDP4_->push_back( it->second );
  }
  for(std::map<double, std::pair<float,float> >::iterator it = bTaggers.begin(); it != bTaggers.end() ; it++){
    jetsBtagHE_->push_back( (it->second).first  );
    jetsBtagHP_->push_back( (it->second).second );
  }
  for(std::map<double, std::pair<float,float> >::iterator it = jetPVassociation.begin(); it != jetPVassociation.end() ; it++){
    jetsChEfraction_->push_back( (it->second).first  );
    jetsChNfraction_->push_back( (it->second).second );
  }
  for(std::map<double, std::pair<float,float> >::iterator it = jetMoments.begin(); it != jetMoments.end() ; it++){
    jetMoments_->push_back( (it->second).first  );
    jetMoments_->push_back( (it->second).second );
  }

  tree_->Fill();

}


unsigned int  ElecTauStreamAnalyzer::jetID( const pat::Jet* jet, const reco::Vertex* vtx, std::vector<float> vtxZ, std::map<std::string,float>& map_){

  if( (jet->pt())<10 ) return 99; // always pass jet ID

  std::vector<reco::PFCandidatePtr> pfCandPtrs = jet->getPFConstituents();

  int nCharged = 0;
  int nPhotons = 0;
  int nNeutral = 0;
  int nConst = 0;

  float energyCharged = 0;
  float energyPhotons = 0;
  float energyNeutral = 0;
  float energyElectrons = 0;
 
  float totalEnergyFromConst = 0;

  float chEnFractionPV   = 0;
  float chEnFractionChPV = 0;
  float chFractionPV     = 0;

  for(unsigned i=0; i<pfCandPtrs.size(); ++i) {
    const reco::PFCandidate& cand = *(pfCandPtrs[i]);

    totalEnergyFromConst +=  cand.energy();
    nConst += 1;

    switch( cand.particleId() ) {
    case reco::PFCandidate::h: 
      nCharged++;
      energyCharged += cand.energy(); 

      if((cand.trackRef()).isNonnull()){
	bool isMatched = false;
	for(reco::Vertex::trackRef_iterator it = vtx->tracks_begin() ; it!=vtx->tracks_end() && !isMatched; it++){
	  if(  (*it).id() == (cand.trackRef()).id() && (*it).key() == (cand.trackRef()).key() ){
	    isMatched = true;
	    if(verbose_) cout << (*it).id() << ", " << (*it).key() << " is matched!" << endl;
	    chEnFractionPV += cand.energy();
	    chFractionPV+=1;
	  }
	}
	if(!isMatched){
	  float dist = vtxZ.size()>0 ? fabs(vtxZ[0]-((cand.trackRef())->referencePoint()).z()) : 999.;
	  float tmpDist = 999.;
	  for(unsigned k = 1; k < vtxZ.size(); k++){
	    if( fabs(((cand.trackRef())->referencePoint()).z() - vtxZ[k] ) < tmpDist )
	      tmpDist = fabs(((cand.trackRef())->referencePoint()).z() - vtxZ[k] );
	  }
	  if(tmpDist>dist){
	    isMatched = true;
	    if(verbose_) cout << "Matched by closest vtx in z!!" << endl;
	    chEnFractionPV += cand.energy();
	    chFractionPV+=1;
	  }
	}
	if(!isMatched && verbose_) {
	  cout << "Ch with pt " << cand.pt() << " and eta " << cand.eta() << " is not matched to the PV !!!" << endl;
	  cout << "z position of PV " << (vtx->position()).z() << ", z position of track " << ((cand.trackRef())->referencePoint()).z()  << ", x position of track " << ((cand.trackRef())->referencePoint()).x()   << ", y position of track " << ((cand.trackRef())->referencePoint()).y() << endl;
	}
      }
      break;
    case reco::PFCandidate::gamma:
      nPhotons++;
      energyPhotons += cand.energy();
      break;
    case reco::PFCandidate::h0:
      nNeutral++;
      energyNeutral += cand.energy();
      break;
    case reco::PFCandidate::e: 
      energyElectrons += cand.energy(); 

      if((cand.gsfTrackRef()).isNonnull()){
	bool isMatched = false;
	float dist = vtxZ.size()>0 ? fabs(vtxZ[0]-((cand.gsfTrackRef())->referencePoint()).z()) : 999.;
	float tmpDist = 999.;
	for(unsigned k = 1; k < vtxZ.size(); k++){
	  if( fabs(((cand.gsfTrackRef())->referencePoint()).z() - vtxZ[k] ) < tmpDist )
	    tmpDist = fabs(((cand.gsfTrackRef())->referencePoint()).z() - vtxZ[k] );
	}
	if(tmpDist>dist){
	  isMatched = true;
	  if(verbose_) cout << "Matched by closest vtx in z!!" << endl;
	  chEnFractionPV += cand.energy();
	  chFractionPV+=1;
	}
	if(!isMatched && verbose_) {
	  cout << "Ele with pt " << cand.pt() << " and eta " << cand.eta() << " is not matched to the PV !!!" << endl;
	  cout << "z position of PV " << (vtx->position()).z() << ", z position of gsfTrack " << ((cand.gsfTrackRef())->referencePoint()).z()  << ", x position of gsfTrack " << ((cand.gsfTrackRef())->referencePoint()).x()   << ", y position of gsfTrack " << ((cand.gsfTrackRef())->referencePoint()).y() << endl;
	}
      }
      
      break;
    case reco::PFCandidate::h_HF: // fill neutral
      nNeutral++;
      //energyNeutral += cand.energy();
      break;
    case reco::PFCandidate::egamma_HF: // fill e/gamma
      nPhotons++;
      energyPhotons += cand.energy();
      break;
    default:
      break;
    }
  }

  chEnFractionChPV = chEnFractionPV;
  if(energyCharged>0 || energyElectrons>0)
    chEnFractionChPV /= (energyCharged+energyElectrons); 
  chEnFractionPV /= jet->correctedJet("Uncorrected").p4().E();
  chFractionPV   /= nCharged;

  //if(chFractionPV==0) cout << energyCharged << endl;

  map_["chMult"]          =  chFractionPV;
  map_["chFracRawJetE"]   =  chEnFractionPV;
  map_["chFracAllChargE"] =  chEnFractionChPV;

  bool loose=false;
  bool medium=false;
  bool tight=false;

  //loose id
  if( (TMath::Abs(jet->eta())>2.4 && 
       energyNeutral/totalEnergyFromConst<0.99 && 
       energyPhotons/totalEnergyFromConst<0.99 &&
       nConst > 1) || 
      (TMath::Abs(jet->eta())<2.4 && 
       energyNeutral/totalEnergyFromConst<0.99 && 
       energyPhotons/totalEnergyFromConst<0.99 &&
       nConst > 1 &&
       energyCharged/totalEnergyFromConst>0 &&
       nCharged>0 &&
       energyElectrons/totalEnergyFromConst<0.99
       )
      ) loose = true;
  // medium id
  if( (TMath::Abs(jet->eta())>2.4 && 
       energyNeutral/totalEnergyFromConst<0.95 && 
       energyPhotons/totalEnergyFromConst<0.95 &&
       nConst > 1) || 
      (TMath::Abs(jet->eta())<2.4 && 
       energyNeutral/totalEnergyFromConst<1 && 
       energyPhotons/totalEnergyFromConst<1 &&
       nConst > 1 &&
       energyCharged/totalEnergyFromConst>0 &&
       nCharged>0 &&
       energyElectrons/totalEnergyFromConst<1
       )
      ) medium = true;
  // tight id
  if( (TMath::Abs(jet->eta())>2.4 && 
       energyNeutral/totalEnergyFromConst<0.90 && 
       energyPhotons/totalEnergyFromConst<0.90 &&
       nConst > 1) || 
      (TMath::Abs(jet->eta())<2.4 && 
       energyNeutral/totalEnergyFromConst<1 && 
       energyPhotons/totalEnergyFromConst<1 &&
       nConst > 1 &&
       energyCharged/totalEnergyFromConst>0 &&
       nCharged>0 &&
       energyElectrons/totalEnergyFromConst<1
       )
      ) tight = true;
  
  if(loose && !medium && !tight) return 1;
  if(loose && medium && !tight)  return 2;
  if(loose && medium && tight)   return 3; 
  
  return 0;

}

pat::Jet* ElecTauStreamAnalyzer::newJetMatched( const pat::Jet* oldJet , const pat::JetCollection* newJets ){
  
  pat::Jet* matchedJet = 0;
  
  for(unsigned int it = 0; it < newJets->size() ; it++){
    if( Geom::deltaR( (*newJets)[it].p4() , oldJet->p4() )<0.01 ){
      matchedJet = const_cast<pat::Jet*>(&((*newJets)[it]));
      break;
    }
  }

  return matchedJet;

}


void ElecTauStreamAnalyzer::endJob(){}


#include "FWCore/Framework/interface/MakerMacros.h"
 
DEFINE_FWK_MODULE(ElecTauStreamAnalyzer);



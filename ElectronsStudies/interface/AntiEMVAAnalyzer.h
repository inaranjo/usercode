/* #ifndef Bianchi_TauTauStudies_AntiEMVAAnalyzer_h */
/* #define Bianchi_TauTauStudies_AntiEMVAAnalyzer_h */
#ifndef IvoNaranjo_ElectronsStudies_AntiEMVAAnalyzer_h
#define IvoNaranjo_ElectronsStudies_AntiEMVAAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"


#include "TFile.h"
#include "TTree.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateElectronExtra.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateElectronExtraFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <TMath.h>

#include <vector>
#include <string>

class AntiEMVAAnalyzer : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit AntiEMVAAnalyzer(const edm::ParameterSet&);
    
  // destructor
  ~AntiEMVAAnalyzer();
  
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();
  
  edm::InputTag srcGsfElectrons_;
  edm::InputTag srcPFTaus_;
  edm::InputTag srcGenElectrons_;
  edm::InputTag srcGenElectronsFromZ_;
  edm::InputTag srcGenElectronsFromZTauTau_;
  edm::InputTag srcGenTaus_;
  edm::InputTag srcGenJets_;
  edm::InputTag srcPrimaryVertex_;
  bool debug_;
  
  TTree* tree_;

  unsigned long run_,event_,lumi_;
  int NumPV_;
  int NumGsfEle_;
  int NumPFTaus_;
  int NumGenEle_;
  int NumGenHad_;
  int NumGenJet_;

  int Elec_GenEleMatch_[50];
  int Elec_GenEleFromZMatch_[50];
  int Elec_GenEleFromZTauTauMatch_[50];
  int Elec_PFTauMatch_[50];
  int Elec_GenHadMatch_[50];
  int Elec_GenJetMatch_[50];
  float Elec_AbsEta_[50];
  float Elec_Pt_[50];
  float Elec_PFMvaOutput_[50];
  float Elec_Ee_[50];
  float Elec_Egamma_[50];
  float Elec_Pin_[50];
  float Elec_Pout_[50];
  float Elec_EtotOverPin_[50];
  float Elec_EeOverPout_[50];
  float Elec_EgammaOverPdif_[50];
  int Elec_EarlyBrem_[50];
  int Elec_LateBrem_[50];
  float Elec_Logsihih_[50];
  float Elec_DeltaEta_[50];
  float Elec_HoHplusE_[50];
  float Elec_Fbrem_[50];
  float Elec_Chi2KF_[50];
  float Elec_Chi2GSF_[50];
  float Elec_NumHits_[50];
  float Elec_GSFTrackResol_[50];
  float Elec_GSFTracklnPt_[50];
  float Elec_GSFTrackEta_[50];

  std::vector<float>   GammasdEta;
  std::vector<float>   GammasdPhi;   
  std::vector<float>   GammasPt;
  int Tau_GsfEleMatch_[50];
  int Tau_GenEleMatch_[50];
  int Tau_GenEleFromZMatch_[50];
  int Tau_GenEleFromZTauTauMatch_[50];
  int Tau_GenHadMatch_[50];
  int Tau_GenJetMatch_[50];
  float Tau_AbsEta_[50];
  float Tau_Pt_[50];
  float Tau_HasGsf_[50]; 
  float Tau_EmFraction_[50]; 
  float Tau_NumChargedCands_[50];
  float Tau_NumGammaCands_[50]; 
  float Tau_HadrHoP_[50]; 
  float Tau_HadrEoP_[50]; 
  float Tau_VisMass_[50]; 
  float Tau_GammaEtaMom_[50];
  float Tau_GammaPhiMom_[50];
  float Tau_GammaEnFrac_[50];
  float Tau_HadrMva_[50]; 

};//AntiEMVAAnalyzer

#endif   

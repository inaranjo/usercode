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

  std::vector<float>   GammasdEta;
  std::vector<float>   GammasdPhi;
  std::vector<float>   GammasPt;
  int Tau_GsfEleMatch_Tab_[50];
  int Tau_GenEleMatch_Tab_[50];
  int Tau_GenEleFromZMatch_Tab_[50];
  int Tau_GenEleFromZTauTauMatch_Tab_[50];
  int Tau_GenHadMatch_Tab_[50];
  int Tau_GenJetMatch_Tab_[50];
  float Tau_AbsEta_Tab_[50];
  float Tau_Pt_Tab_[50];
  float Tau_HasGsf_Tab_[50]; 
  float Tau_EmFraction_Tab_[50]; 
  float Tau_NumChargedCands_Tab_[50];
  float Tau_NumGammaCands_Tab_[50]; 
  float Tau_HadrHoP_Tab_[50]; 
  float Tau_HadrEoP_Tab_[50]; 
  float Tau_VisMass_Tab_[50]; 
  float Tau_GammaEtaMom_Tab_[50];
  float Tau_GammaPhiMom_Tab_[50];
  float Tau_GammaEnFrac_Tab_[50];
  float Tau_HadrMva_Tab_[50]; 

  int Tau_GsfEleMatch_;
  int Tau_GenEleMatch_;
  int Tau_GenEleFromZMatch_;
  int Tau_GenEleFromZTauTauMatch_;
  int Tau_GenHadMatch_;
  int Tau_GenJetMatch_;
  float Tau_AbsEta_;
  float Tau_Pt_;
  float Tau_HasGsf_; 
  float Tau_EmFraction_; 
  float Tau_NumChargedCands_;
  float Tau_NumGammaCands_; 
  float Tau_HadrHoP_; 
  float Tau_HadrEoP_; 
  float Tau_VisMass_; 
  float Tau_GammaEtaMom_;
  float Tau_GammaPhiMom_;
  float Tau_GammaEnFrac_;
  float Tau_HadrMva_; 

  int Elec_GenEleMatch_;
  int Elec_GenEleFromZMatch_;
  int Elec_GenEleFromZTauTauMatch_;
  int Elec_PFTauMatch_;
  int Elec_GenHadMatch_;
  int Elec_GenJetMatch_;
  float Elec_AbsEta_;
  float Elec_Pt_;
  float Elec_HasSC_;
  float Elec_HasKF_;
  float Elec_HasGSF_;
  float Elec_PFMvaOutput_;
  float Elec_Ee_;
  float Elec_Egamma_;
  float Elec_Pin_;
  float Elec_Pout_;
  float Elec_EtotOverPin_;
  float Elec_EeOverPout_;
  float Elec_EgammaOverPdif_;
  int Elec_EarlyBrem_;
  int Elec_LateBrem_;
  float Elec_Logsihih_;
  float Elec_DeltaEta_;
  float Elec_HoHplusE_;
  float Elec_Fbrem_;
  float Elec_Chi2KF_;
  float Elec_Chi2GSF_;
  float Elec_NumHits_;
  float Elec_GSFTrackResol_;
  float Elec_GSFTracklnPt_;
  float Elec_GSFTrackEta_;


};//AntiEMVAAnalyzer

#endif   

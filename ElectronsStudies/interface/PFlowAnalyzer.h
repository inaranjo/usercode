#ifndef Bianchi_TauTauStudies_PFlowAnalyzer_h
#define Bianchi_TauTauStudies_PFlowAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "TFile.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateElectronExtra.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateElectronExtraFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include <TMath.h>

#include <vector>
#include <string>

class PFlowAnalyzer : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit PFlowAnalyzer(const edm::ParameterSet&);
    
  // destructor
  ~PFlowAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag srcPFCandidates_;
  edm::InputTag srcPFJets_;
  edm::InputTag srcPFMEt_;

  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcWeights_;

  std::string dqmDirectory_;

  TH1F* h1_;


  struct plotEntryType
  {
    plotEntryType()
    {}
    ~plotEntryType() {}

    void bookHistograms()
    {
/*       edm::Service<TFileService> fs; */
      //TFileDirectory dir = fs.mkdir(directory);
      //dir.make<TH1F>("h1", "h1", 100, 0., 100.);
    }
    
/*     void bookHistograms(DQMStore& dqmStore) */
/*     { */
/*       std::string dqmDirectory_full = dqmDirectory_; */
/*       if ( dqmSubDirectory_ != "" ) dqmDirectory_full.append("/").append(dqmSubDirectory_); */

/*       dqmStore.setCurrentFolder(dqmDirectory_full); */

/*       pfChargedHadronPt_    = dqmStore.book1D("pfChargedHadronPt",    "pfChargedHadronPt",    100, 0., 100.); */
/*       pfChargedHadronPtSum_ = dqmStore.book1D("pfChargedHadronPtSum", "pfChargedHadronPtSum", 100, 0., 100.); */
/*       pfGammaPt_            = dqmStore.book1D("pfGammaPt",            "pfGammaPt",            100, 0., 100.); */
/*       pfGammaPtSum_         = dqmStore.book1D("pfGammaPtSum",         "pfGammaPtSum",         100, 0., 100.); */
/*       pfNeutralHadronPt_    = dqmStore.book1D("pfNeutralHadronPt",    "pfNeutralHadronPt",    100, 0., 100.); */
/*       pfNeutralHadronPtSum_ = dqmStore.book1D("pfNeutralHadronPtSum", "pfNeutralHadronPtSum", 100, 0., 100.); */
/*       pfElectronPt_         = dqmStore.book1D("pfElectronPt",         "pfElectronPt",         100, 0., 100.); */
/*       pfElectronPtSum_      = dqmStore.book1D("pfElectronPtSum",      "pfElectronPtSum",      100, 0., 100.); */
/*       pfMuonPt_             = dqmStore.book1D("pfMuonPt",             "pfMuonPt",             100, 0., 100.); */
/*       pfMuonPtSum_          = dqmStore.book1D("pfMuonPtSum",          "pfMuonPtSum",          100, 0., 100.); */
      
/*       pfJetPt_              = dqmStore.book1D("pfJetPt",              "pfJetPt",              100, 0., 100.); */
/*       pfJetPtGt10PtSum_     = dqmStore.book1D("pfJetPtGt10PtSum",     "pfJetPtGt10PtSum",     100, 0., 100.); */
/*       pfJetPtGt15PtSum_     = dqmStore.book1D("pfJetPtGt15PtSum",     "pfJetPtGt15PtSum",     100, 0., 100.); */
/*       pfJetPtGt20PtSum_     = dqmStore.book1D("pfJetPtGt20PtSum",     "pfJetPtGt20PtSum",     100, 0., 100.); */
/*       pfJetPtGt25PtSum_     = dqmStore.book1D("pfJetPtGt25PtSum",     "pfJetPtGt25PtSum",     100, 0., 100.); */
/*       pfJetPtGt30PtSum_     = dqmStore.book1D("pfJetPtGt30PtSum",     "pfJetPtGt30PtSum",     100, 0., 100.); */
/*     } */
    void fillHistograms(const reco::PFCandidateCollection& pfCandidates)
    {
      double pfCandidateAbsEta = -99;
      
      for ( reco::PFCandidateCollection::const_iterator pfCandidate = pfCandidates.begin();
	    pfCandidate != pfCandidates.end(); ++pfCandidate ) {
	pfCandidateAbsEta = TMath::Abs(pfCandidate->eta());
	//h1_->Fill(pfCandidateAbsEta);
      }
    }

    std::string dqmDirectory_;
    std::string dqmSubDirectory_;

    double absEtaMin_;
    double absEtaMax_;

    MonitorElement* pfAllCandidatePt_;
    MonitorElement* pfAllCandidatePtSum_;
    MonitorElement* pfChargedHadronPt_;
    MonitorElement* pfChargedHadronPtSum_;
    MonitorElement* pfGammaPt_;
    MonitorElement* pfGammaPtSum_;
    MonitorElement* pfNeutralHadronPt_;
    MonitorElement* pfNeutralHadronPtSum_;
    MonitorElement* pfElectronPt_;
    MonitorElement* pfElectronPtSum_;
    MonitorElement* pfMuonPt_;
    MonitorElement* pfMuonPtSum_;
      
    MonitorElement* pfJetPt_;
    MonitorElement* pfAllJetPtSum_;
    MonitorElement* pfJetPtGt10PtSum_;
    MonitorElement* pfJetPtGt15PtSum_;
    MonitorElement* pfJetPtGt20PtSum_;
    MonitorElement* pfJetPtGt25PtSum_;
    MonitorElement* pfJetPtGt30PtSum_;
  };
  
  plotEntryType* plotsAllEta_;

  plotEntryType* plotsAbsEtaLt11_;
  plotEntryType* plotsAbsEta11to17_;
  plotEntryType* plotsAbsEta17to25_;
  plotEntryType* plotsAbsEtaGt25_;

  MonitorElement* dPhi_vs_pfChargedHadronPtSum_;
  MonitorElement* pfMEt_vs_pfChargedHadronPtSum_;

};

#endif   

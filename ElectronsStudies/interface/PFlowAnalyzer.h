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

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateElectronExtra.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateElectronExtraFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

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
  
  edm::InputTag srcGsfElectrons_;
  bool debug_;
  
  struct plotEntryType
  {
    plotEntryType(const std::string& directory, double absEtaMin, double absEtaMax)
      : directory_(directory),
	absEtaMin_(absEtaMin),
	absEtaMax_(absEtaMax)
    {}
    ~plotEntryType() {}
    
    void bookHistograms() 
    { 
      edm::Service<TFileService> fs;
      TFileDirectory dir = fs->mkdir(directory_);

      hEeOverPout_ = dir.make<TH1F>("hEeOverPout","hEeOverPout",100,0,4);
      hEarlyBrem_ = dir.make<TH1F>("hEarlyBbrem","hEarlyBrem",100,-2,1);
      hLateBrem_ = dir.make<TH1F>("hLateBbrem","hLateBrem",100,-2,1);
      hLogsihih_ = dir.make<TH1F>("hLogsihih","hLogsihih",100,-13,-2);
      hDeltaEta_ = dir.make<TH1F>("hDeltaEta","hDeltaEta",100,0,0.05);
      hHoE_ = dir.make<TH1F>("hHoE","hHoE",100,0,0.5);
      hHoEBc_ = dir.make<TH1F>("hHoEBc","hHoEBc",100,0,0.5);
      hFbrem_ = dir.make<TH1F>("hFbrem","hFbrem",100,-0.2,1.1);
      hChi2KF_ = dir.make<TH1F>("hChi2KF","hChi2KF",50,0,5);
      hChi2GSF_ = dir.make<TH1F>("hChi2GSF","hChi2GSF",50,0,5);
      hNHits_ = dir.make<TH1F>("hNHits","hNHits",50,0,50);
      hGSFResol_ = dir.make<TH1F>("hGSFResol","hGSFResol",100,0,1);
      hGSFlnPt_ = dir.make<TH1F>("hGSFlnPt","hGSFlnPt",100,0,15);
      hGSFEta_ = dir.make<TH1F>("hGSFEta","hGSFEta",100,-2.5,2.5);

    }
    
    void fillHistograms(const reco::GsfElectronCollection& pfElectrons,bool debug_)
    {
      float EeOverPout = -99;
      int EarlyBrem = -99;
      int LateBrem = -99;
      float Logsihih = -99;
      float DeltaEta = -99;
      float HoE = -99;
      float HoEBc = -99;
      float Fbrem = -99;
      float Chi2KF = -99;
      float Chi2GSF = -99;
      float NHits = -99;
      float GSFResol = -99;
      float GSFlnPt = -99;
      float GSFEta = -99;

      for ( reco::GsfElectronCollection::const_iterator pfElectron = pfElectrons.begin();
	    pfElectron != pfElectrons.end(); ++pfElectron ) {	      

	EeOverPout = pfElectron->eEleClusterOverPout();
	EarlyBrem = pfElectron->mvaInput().earlyBrem;
	LateBrem = pfElectron->mvaInput().lateBrem;
	Logsihih = log(pfElectron->mvaInput().sigmaEtaEta);
	DeltaEta = pfElectron->mvaInput().deltaEta;
	HoE = pfElectron->showerShape().hcalDepth2OverEcal ;
	HoEBc = pfElectron->showerShape().hcalDepth2OverEcalBc ;
	Fbrem = pfElectron->fbrem();
	if (pfElectron->closestCtfTrackRef().isNonnull()){
	  Chi2KF = pfElectron->closestCtfTrackRef()->normalizedChi2();
	  NHits = pfElectron->closestCtfTrackRef()->numberOfValidHits();
	}
	if(pfElectron->gsfTrack().isNonnull()){
	  Chi2GSF = pfElectron->gsfTrack()->normalizedChi2();
	  GSFResol = pfElectron->gsfTrack()->ptError()/pfElectron->gsfTrack()->pt();
	  GSFlnPt = log(pfElectron->gsfTrack()->pt())*TMath::Ln10();
	  GSFEta = pfElectron->gsfTrack()->eta();
	}

	if(debug_){
	std::cout<<std::endl;
	std::cout<<"E electron cluster over Pout :   "<<EeOverPout<<std::endl;
	std::cout<<"EarlyBrem :   "<<EarlyBrem<<std::endl;
	std::cout<<"LateBrem :   "<<LateBrem<<std::endl;
	std::cout<<"log(sigma EtaEta with the SC) :  "<<Logsihih<<std::endl;
	std::cout<<"PF-cluster GSF track delta-eta :  "<<DeltaEta<<std::endl;
	std::cout<<"H over E :   "<<HoE<<std::endl;
	std::cout<<"H over E Bc:   "<<HoEBc<<std::endl;
	std::cout<<"FBrem :   "<<Fbrem<<std::endl;
	std::cout<<"Normalized Chi2 KF :   "<<Chi2KF<<std::endl;
	std::cout<<"Normalized Chi2 GSF :   "<<Chi2GSF<<std::endl;
	std::cout<<"Number of valid hits KF :   "<<NHits<<std::endl;
	std::cout<<"GSF sig(Pt)/Pt :   "<<GSFResol<<std::endl;
	std::cout<<"GSF ln(Pt) :   "<<GSFlnPt<<std::endl;
	std::cout<<"GSF Eta :   "<<GSFEta<<std::endl;
	}
	hEeOverPout_->Fill(EeOverPout);
	hEarlyBrem_->Fill(EarlyBrem);
	hLateBrem_->Fill(LateBrem);
	hLogsihih_->Fill(Logsihih);
	hDeltaEta_->Fill(DeltaEta);
	hHoE_->Fill(HoE);
	hHoEBc_->Fill(HoEBc);
	hFbrem_->Fill(Fbrem);
	hChi2KF_->Fill(Chi2KF);
	hChi2GSF_->Fill(Chi2GSF);
	hNHits_->Fill(NHits);
	hGSFResol_->Fill(GSFResol);
	hGSFlnPt_->Fill(GSFlnPt);
	hGSFEta_->Fill(GSFEta);

      }// loop electrons
    }//fillHistograms
  
    std::string directory_;
    
    double absEtaMin_;
    double absEtaMax_;
    
    TH1F* hEeOverPout_;
    TH1F* hEarlyBrem_;
    TH1F* hLateBrem_;
    TH1F* hLogsihih_;
    TH1F* hDeltaEta_;   
    TH1F* hHoE_;
    TH1F* hHoEBc_;
    TH1F* hFbrem_;
    TH1F* hChi2KF_;
    TH1F* hChi2GSF_;
    TH1F* hNHits_;
    TH1F* hGSFResol_;
    TH1F* hGSFlnPt_;
    TH1F* hGSFEta_;

  };//plotEntryType
  
  plotEntryType* plotsAllEta_;
  
};//PFlowAnalyzer

#endif   

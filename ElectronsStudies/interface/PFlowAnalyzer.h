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
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

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
  edm::InputTag srcPrimaryVertex_;
  bool debug_;
  
  struct plotEntryType
  {
    plotEntryType(const std::string& directory, double absEtaMin, double absEtaMax, int numPVMin,  int numPVMax, double ptMin, double ptMax)
      : directory_(directory),
	absEtaMin_(absEtaMin),
	absEtaMax_(absEtaMax),
	numPVMin_(numPVMin),
	numPVMax_(numPVMax),
	ptMin_(ptMin),
	ptMax_(ptMax)
    {}
    ~plotEntryType() {}
    
    void bookHistograms(double absEtaMin, double absEtaMax, int numPVMin,  int numPVMax, double ptMin, double ptMax ) 
    { 
      edm::Service<TFileService> fs;
      TFileDirectory dir = fs->mkdir(directory_);

      hNumPV_ = dir.make<TH1F>("hNumPV","hNumPV",numPVMax-numPVMin,numPVMin,numPVMax);
      hElecAbsEta_ = dir.make<TH1F>("hElecAbsEta","hElecAbsEta",100,absEtaMin,absEtaMax);
      hElecPt_ = dir.make<TH1F>("hElecPt","hElecPt",100,ptMin,ptMax);
      hEtotOverPin_ = dir.make<TH1F>("hEtotOverPin","hEtotOverPin",100,0,4);
      hEeOverPout_ = dir.make<TH1F>("hEeOverPout","hEeOverPout",100,0,4);
      hEgammaOverPdif_ = dir.make<TH1F>("hEgammaOverPdif","hEgammaOverPdif",100,0,4);
      hEarlyBrem_ = dir.make<TH1F>("hEarlyBbrem","hEarlyBrem",100,-2,1);
      hLateBrem_ = dir.make<TH1F>("hLateBbrem","hLateBrem",100,-2,1);
      hLogsihih_ = dir.make<TH1F>("hLogsihih","hLogsihih",100,-13,-2);
      hDeltaEta_ = dir.make<TH1F>("hDeltaEta","hDeltaEta",100,0,0.05);
      hHoE_ = dir.make<TH1F>("hHoE","hHoE",100,0,0.5);
      hHoEBc_ = dir.make<TH1F>("hHoEBc","hHoEBc",100,0,0.5);
      hFbrem_ = dir.make<TH1F>("hFbrem","hFbrem",100,-0.2,1.1);
      hChi2KF_ = dir.make<TH1F>("hChi2KF","hChi2KF",50,0,5);
      hChi2GSF_ = dir.make<TH1F>("hChi2GSF","hChi2GSF",50,0,5);
      hNHits_ = dir.make<TH1F>("hNHits","hNHits",30,0,30);
      hGSFResol_ = dir.make<TH1F>("hGSFResol","hGSFResol",100,0,1);
      hGSFlnPt_ = dir.make<TH1F>("hGSFlnPt","hGSFlnPt",100,0,15);
      hGSFEta_ = dir.make<TH1F>("hGSFEta","hGSFEta",100,-2.5,2.5);

    }
    
    void fillHistograms(const reco::GsfElectronCollection& Electrons,bool debug_, int numPV_)
    {
      int numPV = numPV_;
      double ElecAbsEta = -99;
      double ElecPt = -99;

      double Ee = -99;
      double Egamma = -99;

      float EtotOverPin = -99;
      float EeOverPout = -99;
      float EgammaOverPdif = -99;
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

      
      if(numPV>=numPVMin_ && numPV<numPVMax_){
	
	for ( reco::GsfElectronCollection::const_iterator Electron = Electrons.begin();
	      Electron != Electrons.end(); ++Electron ) {	      
	  
	  ElecAbsEta = TMath::Abs(Electron->eta());
	  ElecPt = Electron->pt();

	  if(ElecAbsEta>absEtaMin_ && ElecAbsEta<absEtaMax_ && ElecPt>ptMin_ && ElecPt< ptMax_) {

	    reco::SuperClusterRef pfSuperCluster = Electron->pflowSuperCluster();
	    if(pfSuperCluster.isNonnull() && pfSuperCluster.isAvailable()){
	      Ee = 0.;
	      Egamma = 0.;
	      //std::cout<<"SuperCluster accessed   "<<std::endl;
	      for (reco::CaloCluster_iterator pfCluster = pfSuperCluster->clustersBegin();
		   pfCluster != pfSuperCluster->clustersEnd(); ++pfCluster ) {
		double pfClusterEn = (*pfCluster)->energy();
		if ( pfCluster == pfSuperCluster->clustersBegin() ) Ee += pfClusterEn;
		else Egamma += pfClusterEn;
	      }
	    }
	    
	    //std::cout<<"Ee :   "<<Ee<<" Egamma :  "<<Egamma<<std::endl;
	    
	    double Pin = TMath::Sqrt(Electron->trackMomentumAtVtx().Mag2());
	    double Pout = TMath::Sqrt(Electron->trackMomentumOut().Mag2()); 
	    if (debug_)std::cout<<"Pin :   "<<Pin<<" Pout :  "<<Pout<<std::endl;
	    
	    

	    EtotOverPin = (Ee+Egamma)/Pin;
	    EeOverPout = Electron->eEleClusterOverPout();
	    EgammaOverPdif = Egamma/(Pin-Pout);
	    EarlyBrem = Electron->mvaInput().earlyBrem;
	    LateBrem = Electron->mvaInput().lateBrem;
	    Logsihih = log(Electron->mvaInput().sigmaEtaEta);
	    DeltaEta = Electron->mvaInput().deltaEta;
	    HoE = Electron->showerShape().hcalDepth2OverEcal ;
	    HoEBc = Electron->showerShape().hcalDepth2OverEcalBc ;
	    Fbrem = Electron->fbrem();
	    if (Electron->closestCtfTrackRef().isNonnull()){
	      Chi2KF = Electron->closestCtfTrackRef()->normalizedChi2();
	      NHits = Electron->closestCtfTrackRef()->numberOfValidHits();
	    }
	    if(Electron->gsfTrack().isNonnull()){
	      Chi2GSF = Electron->gsfTrack()->normalizedChi2();
	      GSFResol = Electron->gsfTrack()->ptError()/Electron->gsfTrack()->pt();
	      GSFlnPt = log(Electron->gsfTrack()->pt())*TMath::Ln10();
	      GSFEta = Electron->gsfTrack()->eta();
	    }
	    
	    if(debug_){
	      std::cout<<std::endl;
	      std::cout<<"NumPV :"<<numPV<<std::endl;
	      std::cout<<"ElecAbsEta :"<<ElecAbsEta<<std::endl;
	      std::cout<<"E electron cluster plus photons over Pin :   "<<EtotOverPin<<std::endl;
	      std::cout<<"E electron cluster over Pout :   "<<EeOverPout<<std::endl;
	      std::cout<<"E photons over (Pin - Pout) :   "<<EgammaOverPdif<<std::endl;
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
	    
	    hNumPV_->Fill(numPV);
	    hElecAbsEta_->Fill(ElecAbsEta);
	    hElecPt_->Fill(ElecPt);
	    hEtotOverPin_->Fill(EtotOverPin);
	    hEeOverPout_->Fill(EeOverPout);
	    hEgammaOverPdif_->Fill(EgammaOverPdif);
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
	    
	  }//ElecAbsEta condition
	  
	}// loop electrons
	
      }//numPV condition
      
    }//fillHistograms
  
    std::string directory_;
    
    double absEtaMin_;
    double absEtaMax_;
    double ptMin_;
    double ptMax_;
    int numPVMin_;
    int numPVMax_;

    
    TH1F* hNumPV_;
    TH1F* hElecAbsEta_;
    TH1F* hElecPt_;
    TH1F* hEtotOverPin_;
    TH1F* hEeOverPout_;
    TH1F* hEgammaOverPdif_;
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
  
  plotEntryType* plotsAll_;
  plotEntryType* plotsNumPV0to10_Barrel_Pt20to30_;
  plotEntryType* plotsNumPV10to20_Barrel_Pt20to30_;
  plotEntryType* plotsNumPV20to30_Barrel_Pt20to30_;
  plotEntryType* plotsNumPV30to40_Barrel_Pt20to30_;
  plotEntryType* plotsNumPV0to10_Endcap_Pt20to30_;
  plotEntryType* plotsNumPV10to20_Endcap_Pt20to30_;
  plotEntryType* plotsNumPV20to30_Endcap_Pt20to30_;
  plotEntryType* plotsNumPV30to40_Endcap_Pt20to30_;
  plotEntryType* plotsNumPV0to10_Barrel_Pt30to40_;
  plotEntryType* plotsNumPV10to20_Barrel_Pt30to40_;
  plotEntryType* plotsNumPV20to30_Barrel_Pt30to40_;
  plotEntryType* plotsNumPV30to40_Barrel_Pt30to40_;
  plotEntryType* plotsNumPV0to10_Endcap_Pt30to40_;
  plotEntryType* plotsNumPV10to20_Endcap_Pt30to40_;
  plotEntryType* plotsNumPV20to30_Endcap_Pt30to40_;
  plotEntryType* plotsNumPV30to40_Endcap_Pt30to40_;
  plotEntryType* plotsNumPV0to10_Barrel_Pt40to50_;
  plotEntryType* plotsNumPV10to20_Barrel_Pt40to50_;
  plotEntryType* plotsNumPV20to30_Barrel_Pt40to50_;
  plotEntryType* plotsNumPV30to40_Barrel_Pt40to50_;
  plotEntryType* plotsNumPV0to10_Endcap_Pt40to50_;
  plotEntryType* plotsNumPV10to20_Endcap_Pt40to50_;
  plotEntryType* plotsNumPV20to30_Endcap_Pt40to50_;
  plotEntryType* plotsNumPV30to40_Endcap_Pt40to50_;

  int numPV_;
};//PFlowAnalyzer

#endif   

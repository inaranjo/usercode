#ifndef Bianchi_TauTauStudies_AntiEMVAVariablesAnalyzer_h
#define Bianchi_TauTauStudies_AntiEMVAVariablesAnalyzer_h

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

class AntiEMVAVariablesAnalyzer : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit AntiEMVAVariablesAnalyzer(const edm::ParameterSet&);
    
  // destructor
  ~AntiEMVAVariablesAnalyzer();
  
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();
  
  edm::InputTag srcGsfElectrons_;
  edm::InputTag srcPFTaus_;
  edm::InputTag srcGenElectrons_;
  edm::InputTag srcGenTaus_;
  edm::InputTag srcGenJets_;
  edm::InputTag srcPrimaryVertex_;
  int match_;// 0 no matching, 1 match to GenElectrons, 2 match to GenTaus, 3 match to GenJets
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
      hEarlyBrem_ = dir.make<TH1F>("hEarlyBrem","hEarlyBrem",100,-2,1);
      hLateBrem_ = dir.make<TH1F>("hLateBrem","hLateBrem",100,-2,1);
      hLogsihih_ = dir.make<TH1F>("hLogsihih","hLogsihih",50,-13,-2);
      hDeltaEta_ = dir.make<TH1F>("hDeltaEta","hDeltaEta",100,0,0.05);
      hHoHplusE_ = dir.make<TH1F>("hHoHplusE","hHoHplusE",100,0,1);
      hFbrem_ = dir.make<TH1F>("hFbrem","hFbrem",100,-0.2,1.1);
      hChi2KF_ = dir.make<TH1F>("hChi2KF","hChi2KF",50,0,5);
      hChi2GSF_ = dir.make<TH1F>("hChi2GSF","hChi2GSF",50,0,5);
      hNHits_ = dir.make<TH1F>("hNHits","hNHits",30,0,30);
      hGSFResol_ = dir.make<TH1F>("hGSFResol","hGSFResol",100,0,1);
      hGSFlnPt_ = dir.make<TH1F>("hGSFlnPt","hGSFlnPt",100,0,15);
      hGSFEta_ = dir.make<TH1F>("hGSFEta","hGSFEta",100,-2.5,2.5);

      hTauAbsEta_ = dir.make<TH1F>("hTauAbsEta","hTauAbsEta",100,absEtaMin,absEtaMax);
      hTauPt_ = dir.make<TH1F>("hTauPt","hTauPt",100,ptMin,ptMax);
      hTauHasGsf_ = dir.make<TH1F>("hTauHasGsf","hTauHasGsf",2,0,2);
      hTauEmFraction_ = dir.make<TH1F>("hTauEmFraction","hTauEmFraction",100,0,1);
      hTauSignalPFChargedCands_ = dir.make<TH1F>("hTauSignalPFChargedCands","hTauSignalPFChargedCands",5,0,5);
      hTauSignalPFGammaCands_ = dir.make<TH1F>("hTauSignalPFGammaCands","hTauSignalPFGammaCands",5,0,5);
      hTauLeadPFChargedHadrHoP_ = dir.make<TH1F>("hTauLeadPFChargedHadrHoP","hTauLeadPFChargedHadrHoP",20,0,1);
      hTauLeadPFChargedHadrEoP_ = dir.make<TH1F>("hTauLeadPFChargedHadrEoP","hTauLeadPFChargedHadrEoP",20,0,1);
      hTauVisMass_NoGsf_ = dir.make<TH1F>("hTauVisMass_NoGsf","hTauVisMass_NoGsf",34,0,1.8);
      hGammaEtaMom_NoGsf_ = dir.make<TH1F>("hGammaEtaMom_NoGsf","hGammaEtaMom_NoGsf",30,0,3);
      hGammaPhiMom_NoGsf_ = dir.make<TH1F>("hGammaPhiMom_NoGsf","hGammaPhiMom_NoGsf",50,0,5);
      hGammaEnFrac_NoGsf_ = dir.make<TH1F>("hGammaEnFrac_NoGsf","hGammaEnFrac_NoGsf",20,0,1);

      hTauVisMass_HasGsf_ = dir.make<TH1F>("hTauVisMass_HasGsf","hTauVisMass_HasGsf",34,0,1.8);
      hGammaEtaMom_HasGsf_ = dir.make<TH1F>("hGammaEtaMom_HasGsf","hGammaEtaMom_HasGsf",30,0,3);
      hGammaPhiMom_HasGsf_ = dir.make<TH1F>("hGammaPhiMom_HasGsf","hGammaPhiMom_HasGsf",50,0,5);
      hGammaEnFrac_HasGsf_ = dir.make<TH1F>("hGammaEnFrac_HasGsf","hGammaEnFrac_HasGsf",20,0,1);
      hTauLeadPFChargedHadrMva_HasGsf_ = dir.make<TH1F>("hTauLeadPFChargedHadrMva_HasGsf","hTauLeadPFChargedHadrMva_HasGsf",15,-1,0);

    }
    
    void fillHistograms(const reco::GsfElectronCollection& GsfElectrons, 
			const reco::CandidateView& GenElectrons,
			const reco::PFTauCollection& PfTaus,
			const reco::CandidateView& GenTaus, 
			const reco::CandidateView& GenJets, 
			bool debug_, 
			int numPV_,
			int match_
			)
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
      float HoHplusE = -99;
      float Fbrem = -99;
      float Chi2KF = -99;
      float Chi2GSF = -99;
      float NHits = -99;
      float GSFResol = -99;
      float GSFlnPt = -99;
      float GSFEta = -99;

      
      if(numPV>=numPVMin_ && numPV<numPVMax_){
	hNumPV_->Fill(numPV);
	if(debug_)std::cout<<"NumPV :"<<numPV<<std::endl;

	int countEle = 0;
	for ( reco::GsfElectronCollection::const_iterator GsfElectron = GsfElectrons.begin();
	      GsfElectron != GsfElectrons.end(); ++GsfElectron ) {
	  if (debug_){
	    std::cout<<std::endl;
	    std::cout<<"GsfElectron number : "<<countEle<<std::endl;
	  }
	  //////////////Matching GsfElectron with GenElectron//////////////
	  if (match_ == 1){

	    bool GsfEleGenEleMatch = false;
	    int countGenEle = 0;

	    for ( reco::CandidateView::const_iterator GenElectron = GenElectrons.begin();
		  GenElectron != GenElectrons.end(); ++GenElectron ) {
	      if(debug_){
		std::cout<<"  GenElectron number : "<<countGenEle<<std::endl;
		std::cout<<"  DeltaR GsfEle-GenEle :"<<deltaR(GsfElectron->eta(),GsfElectron->phi(),GenElectron->eta(),GenElectron->phi())<<std::endl;
	      }
	      countGenEle++;
	      if(deltaR(GsfElectron->eta(),GsfElectron->phi(),GenElectron->eta(),GenElectron->phi())<0.3)GsfEleGenEleMatch = true;
	    }
	    if(!GsfEleGenEleMatch)continue;
	  }
 	  //////////////Matching GsfElectron with GenElectron//////////////
	  
	  
 	  //////////////Matching GsfElectron with GenHadrons from Ztautau//////////////
	  if(match_ == 2){
	    bool GsfEleGenTauMatch = false;
	    int countGenTau = 0;
	    
	    for ( reco::CandidateView::const_iterator GenTau = GenTaus.begin();
		  GenTau != GenTaus.end(); ++GenTau ) {
	      if(debug_){
		std::cout<<"  GenTau number : "<<countGenTau<<std::endl;
		std::cout<<"  DeltaR GsfEle-GenTau :"<<deltaR(GsfElectron->eta(),GsfElectron->phi(),GenTau->eta(),GenTau->phi())<<std::endl;
	      }
	      countGenTau++;
	      if(deltaR(GsfElectron->eta(),GsfElectron->phi(),GenTau->eta(),GenTau->phi())<0.3)GsfEleGenTauMatch = true;
	    }
	  if(!GsfEleGenTauMatch)continue;
	  }	  
 	  //////////////Matching GsfElectron with GenHadrons from Ztautau//////////////
	  

 	  //////////////Matching GsfElectron with GenJets//////////////
	  if(match_ == 3){
	    bool GsfEleGenJetMatch = false;
	    int countGenJet = 0;
	    
	    for ( reco::CandidateView::const_iterator GenJet = GenJets.begin();
		  GenJet != GenJets.end(); ++GenJet ) {
	      if(debug_){
		std::cout<<"  GenJet number : "<<countGenJet<<std::endl;
		std::cout<<"  DeltaR GsfEle-GenJet :"<<deltaR(GsfElectron->eta(),GsfElectron->phi(),GenJet->eta(),GenJet->phi())<<std::endl;
	      }
	      countGenJet++;
	      if(deltaR(GsfElectron->eta(),GsfElectron->phi(),GenJet->eta(),GenJet->phi())<0.3)GsfEleGenJetMatch = true;
	    }
	  if(!GsfEleGenJetMatch)continue;
	  }	  
 	  //////////////Matching GsfElectron with GenJets //////////////
	  


	  ElecAbsEta = TMath::Abs(GsfElectron->eta());
	  ElecPt = GsfElectron->pt();

	  if(ElecAbsEta>absEtaMin_ && ElecAbsEta<absEtaMax_ && ElecPt>ptMin_ && ElecPt< ptMax_) {

	    reco::SuperClusterRef pfSuperCluster = GsfElectron->pflowSuperCluster();
	    if(pfSuperCluster.isNonnull() && pfSuperCluster.isAvailable()){
	      Ee = 0.;
	      Egamma = 0.;
	      if (debug_)std::cout<<"SuperCluster accessed   "<<std::endl;
	      for (reco::CaloCluster_iterator pfCluster = pfSuperCluster->clustersBegin();
		   pfCluster != pfSuperCluster->clustersEnd(); ++pfCluster ) {
		double pfClusterEn = (*pfCluster)->energy();
		if ( pfCluster == pfSuperCluster->clustersBegin() ) Ee += pfClusterEn;
		else Egamma += pfClusterEn;
	      }
	    }
	    
	    if(debug_)std::cout<<"Ee :   "<<Ee<<" Egamma :  "<<Egamma<<std::endl;
	    
	    double Pin = TMath::Sqrt(GsfElectron->trackMomentumAtVtx().Mag2());
	    double Pout = TMath::Sqrt(GsfElectron->trackMomentumOut().Mag2()); 
	    if (debug_)std::cout<<"Pin :   "<<Pin<<" Pout :  "<<Pout<<std::endl;

	    EtotOverPin = (Ee+Egamma)/Pin;
	    EeOverPout = Ee/Pout;
	    EgammaOverPdif = Egamma/(Pin-Pout);
	    EarlyBrem = GsfElectron->mvaInput().earlyBrem;
	    LateBrem = GsfElectron->mvaInput().lateBrem;
	    Logsihih = log(GsfElectron->mvaInput().sigmaEtaEta);
	    DeltaEta = GsfElectron->mvaInput().deltaEta;
	    HoHplusE = GsfElectron->mvaInput().hadEnergy/(GsfElectron->mvaInput().hadEnergy+Ee) ;
	    Fbrem = GsfElectron->fbrem();
	    if (GsfElectron->closestCtfTrackRef().isNonnull()){
	      Chi2KF = GsfElectron->closestCtfTrackRef()->normalizedChi2();
	      NHits = GsfElectron->closestCtfTrackRef()->numberOfValidHits();
	    }
	    if(GsfElectron->gsfTrack().isNonnull()){
	      Chi2GSF = GsfElectron->gsfTrack()->normalizedChi2();
	      GSFResol = GsfElectron->gsfTrack()->ptError()/GsfElectron->gsfTrack()->pt();
	      GSFlnPt = log(GsfElectron->gsfTrack()->pt())*TMath::Ln10();
	      GSFEta = GsfElectron->gsfTrack()->eta();
	    }
	    
	    if(debug_){
	      std::cout<<std::endl;
	      std::cout<<"ElecAbsEta :"<<ElecAbsEta<<std::endl;
	      std::cout<<"E electron cluster plus photons over Pin :   "<<EtotOverPin<<std::endl;
	      std::cout<<"E electron cluster over Pout :   "<<EeOverPout<<std::endl;
	      std::cout<<"E photons over (Pin - Pout) :   "<<EgammaOverPdif<<std::endl;
	      std::cout<<"EarlyBrem :   "<<EarlyBrem<<std::endl;
	      std::cout<<"LateBrem :   "<<LateBrem<<std::endl;
	      std::cout<<"log(sigma EtaEta with the SC) :  "<<Logsihih<<std::endl;
	      std::cout<<"PF-cluster GSF track delta-eta :  "<<DeltaEta<<std::endl;
	      std::cout<<"H over H plus E :   "<<HoHplusE<<std::endl;
	      std::cout<<"FBrem :   "<<Fbrem<<std::endl;
	      std::cout<<"Normalized Chi2 KF :   "<<Chi2KF<<std::endl;
	      std::cout<<"Normalized Chi2 GSF :   "<<Chi2GSF<<std::endl;
	      std::cout<<"Number of valid hits KF :   "<<NHits<<std::endl;
	      std::cout<<"GSF sig(Pt)/Pt :   "<<GSFResol<<std::endl;
	      std::cout<<"GSF ln(Pt) :   "<<GSFlnPt<<std::endl;
	      std::cout<<"GSF Eta :   "<<GSFEta<<std::endl;
	    }
	    
	    hElecAbsEta_->Fill(ElecAbsEta);
	    hElecPt_->Fill(ElecPt);
	    hEtotOverPin_->Fill(EtotOverPin);
	    hEeOverPout_->Fill(EeOverPout);
	    hEgammaOverPdif_->Fill(EgammaOverPdif);
	    hEarlyBrem_->Fill(EarlyBrem);
	    hLateBrem_->Fill(LateBrem);
	    hLogsihih_->Fill(Logsihih);
	    hDeltaEta_->Fill(DeltaEta);
	    hHoHplusE_->Fill(HoHplusE);
	    hFbrem_->Fill(Fbrem);
	    hChi2KF_->Fill(Chi2KF);
	    hChi2GSF_->Fill(Chi2GSF);
	    hNHits_->Fill(NHits);
	    hGSFResol_->Fill(GSFResol);
	    hGSFlnPt_->Fill(GSFlnPt);
	    hGSFEta_->Fill(GSFEta);
	    
	  }//ElecAbsEta condition
	 
	  countEle++;
	  
	}// loop GsfElectrons
	

	double TauAbsEta = -99;
	double TauPt = -99;
	float TauEmFraction           = -99; 
	float TauHasGsf               = -99; 
	float TauSignalPFChargedCands = -99;
	float TauSignalPFGammaCands   = -99; 

	std::vector<float>   GammasdEta;
	std::vector<float>   GammasdPhi;   
	std::vector<float>   GammasPt;     

	float TauLeadPFChargedHadrHoP = -99; 
	float TauLeadPFChargedHadrEoP = -99; 

	float TauVisMass_NoGsf              = -99; 
	float GammaEtaMom_NoGsf             = -99;
	float GammaPhiMom_NoGsf             = -99;
	float GammaEnFrac_NoGsf             = -99;
	float TauVisMass_HasGsf              = -99; 
	float GammaEtaMom_HasGsf             = -99;
	float GammaPhiMom_HasGsf             = -99;
	float GammaEnFrac_HasGsf             = -99;
	float TauLeadPFChargedHadrMva_HasGsf = -99; 


	int countPfTau = 0;
	for ( reco::PFTauCollection::const_iterator PfTau = PfTaus.begin();
	      PfTau != PfTaus.end(); ++PfTau ) {
	  if (debug_){
	    std::cout<<std::endl;
	    std::cout<<"PfTau number : "<<countPfTau<<std::endl;
	  }
	  countPfTau++;
 	  //////////////Matching PfTau with GenElectron//////////////
	  bool PfTauGenEleMatch = true;
	  if(match_ == 1){
	    PfTauGenEleMatch = false;
	    int countGenEle = 0;
	    
	    for ( reco::CandidateView::const_iterator GenElectron = GenElectrons.begin();
		  GenElectron != GenElectrons.end(); ++GenElectron ) {
	      	    if(debug_){
	      std::cout<<"  GenElectron number : "<<countGenEle<<std::endl;
	      std::cout<<"  DeltaR PfTau-GenEle :"<<deltaR(PfTau->eta(),PfTau->phi(),GenElectron->eta(),GenElectron->phi())<<std::endl;
	      	    }
	      countGenEle++;
	      if(deltaR(PfTau->eta(),PfTau->phi(),GenElectron->eta(),GenElectron->phi())<0.3)PfTauGenEleMatch = true;
	    }
	  }
	  if(!PfTauGenEleMatch)continue;
 	  //////////////Matching PfTau with GenElectron//////////////


	  //////////////Matching PfTau with GenHadrons from Ztautau//////////////
	  bool PfTauGenTauMatch = true;
	  if(match_ == 2){
	    PfTauGenTauMatch = false;
	    int countGenTau = 0;
	    
	    for ( reco::CandidateView::const_iterator GenTau = GenTaus.begin();
		  GenTau != GenTaus.end(); ++GenTau ) {
	      if (debug_){
		std::cout<<"  GenTau number : "<<countGenTau<<std::endl;
		std::cout<<"  DeltaR PfTau-GenTau :"<<deltaR(PfTau->eta(),PfTau->phi(),GenTau->eta(),GenTau->phi())<<std::endl;
	      }
	      countGenTau++;
	      if(deltaR(PfTau->eta(),PfTau->phi(),GenTau->eta(),GenTau->phi())<0.3)PfTauGenTauMatch = true;
	    }
	  }
	  if(!PfTauGenTauMatch)continue;
 	  //////////////Matching PfTau with GenHadrons from Ztautau//////////////
	  

	  //////////////Matching PfTau with GenJets//////////////
	  bool PfTauGenJetMatch = true;
	  if(match_ == 3){
	    PfTauGenJetMatch = false;
	    int countGenJet = 0;
	    
	    for ( reco::CandidateView::const_iterator GenJet = GenJets.begin();
		  GenJet != GenJets.end(); ++GenJet ) {
	      if (debug_){
		std::cout<<"  GenJet number : "<<countGenJet<<std::endl;
		std::cout<<"  DeltaR PfTau-GenJet :"<<deltaR(PfTau->eta(),PfTau->phi(),GenJet->eta(),GenJet->phi())<<std::endl;
	      }
	      countGenJet++;
	      if(deltaR(PfTau->eta(),PfTau->phi(),GenJet->eta(),GenJet->phi())<0.3)PfTauGenJetMatch = true;
	    }
	  }
	  if(!PfTauGenJetMatch)continue;
 	  //////////////Matching PfTau with GenJets//////////////
	  


	  TauAbsEta = TMath::Abs(PfTau->eta());
	  TauPt = PfTau->pt();

	  if(TauAbsEta>absEtaMin_ && TauAbsEta<absEtaMax_ && TauPt>ptMin_ && TauPt< ptMax_) {

	    for(unsigned int k = 0 ; k < (PfTau->signalPFGammaCands()).size() ; k++){
	      reco::PFCandidateRef gamma = (PfTau->signalPFGammaCands()).at(k);
	      if( (PfTau->leadPFChargedHadrCand()).isNonnull() ){
		GammasdEta.push_back( gamma->eta() - PfTau->leadPFChargedHadrCand()->eta() );
		GammasdPhi.push_back( gamma->phi() - PfTau->leadPFChargedHadrCand()->phi() );
	      }
	      else{
		GammasdEta.push_back( gamma->eta() - PfTau->eta() );
		GammasdPhi.push_back( gamma->phi() - PfTau->phi() );
	      }
	      GammasPt.push_back(  gamma->pt() );
	    }
	    
	    float sumPt  = 0;
	    float dEta   = 0;
	    float dEta2  = 0;
	    float dPhi   = 0;
	    float dPhi2  = 0;
	    float sumPt2 = 0;
	    
	    for(unsigned int k = 0 ; k < GammasPt.size() ; k++){
	      float pt_k  = GammasPt[k];
	      float phi_k = GammasdPhi[k];
	      if (GammasdPhi[k] > TMath::Pi()) phi_k = GammasdPhi[k] - 2*TMath::Pi();
	      else if(GammasdPhi[k] < -TMath::Pi()) phi_k = GammasdPhi[k] + 2*TMath::Pi();
	      float eta_k = GammasdEta[k];
	      sumPt  +=  pt_k;
	      sumPt2 += (pt_k*pt_k);
	      dEta   += (pt_k*eta_k);
	      dEta2  += (pt_k*eta_k*eta_k);
	      dPhi   += (pt_k*phi_k);
	      dPhi2  += (pt_k*phi_k*phi_k);  
	    }
	    
	    float GammadPt_ = sumPt/PfTau->pt();
	    
	  if(sumPt>0){
	    dEta  /= sumPt;
	    dPhi  /= sumPt;
	    dEta2 /= sumPt;
	    dPhi2 /= sumPt;
	  }
	  
	  TauEmFraction           = TMath::Max(PfTau->emFraction(),float(0.0));
	  TauSignalPFChargedCands = PfTau->signalPFChargedHadrCands().size();
	  TauSignalPFGammaCands   = PfTau->signalPFGammaCands().size();
	  TauHasGsf               = (PfTau->leadPFChargedHadrCand()->gsfTrackRef()).isNonnull();
	  if(TauSignalPFGammaCands==0){
	  TauLeadPFChargedHadrHoP = PfTau->leadPFChargedHadrCand()->hcalEnergy()/PfTau->leadPFChargedHadrCand()->p();
	  TauLeadPFChargedHadrEoP = PfTau->leadPFChargedHadrCand()->ecalEnergy()/PfTau->leadPFChargedHadrCand()->p();
	  }
	  if(TauSignalPFGammaCands>0 && TauHasGsf==0){
	  GammaEtaMom_NoGsf = TMath::Sqrt(dEta2)*TMath::Sqrt(GammadPt_)*PfTau->pt();
	  GammaPhiMom_NoGsf = TMath::Sqrt(dPhi2)*TMath::Sqrt(GammadPt_)*PfTau->pt();  
	  GammaEnFrac_NoGsf = GammadPt_;
	  TauVisMass_NoGsf              = PfTau->mass();
	  }
	  if(TauSignalPFGammaCands>0 && TauHasGsf==1){
	  GammaEtaMom_HasGsf = TMath::Sqrt(dEta2)*TMath::Sqrt(GammadPt_)*PfTau->pt();
	  GammaPhiMom_HasGsf = TMath::Sqrt(dPhi2)*TMath::Sqrt(GammadPt_)*PfTau->pt();  
	  GammaEnFrac_HasGsf = GammadPt_;
	  TauVisMass_HasGsf  = PfTau->mass();
	  TauLeadPFChargedHadrMva_HasGsf = TMath::Max(PfTau->electronPreIDOutput(),float(-1.0));
	  }


	    if(debug_){
	      std::cout<<std::endl;
	      std::cout<<"TauAbsEta :"<<TauAbsEta<<std::endl;
	      std::cout<<"TauPt :"<<TauPt<<std::endl;
	      std::cout<<"TauHasGsf :"<<TauHasGsf<<std::endl;
	      std::cout<<"TauEmFraction :"<<TauEmFraction<<std::endl;
	      std::cout<<"TauSignalPFChargedCands :"<<TauSignalPFChargedCands<<std::endl;
	      std::cout<<"TauSignalPFGammaCands :"<<TauSignalPFGammaCands<<std::endl;

	      std::cout<<"TauLeadPFChargedHadrHoP :"<<TauLeadPFChargedHadrHoP<<std::endl;
	      std::cout<<"TauLeadPFChargedHadrEoP :"<<TauLeadPFChargedHadrEoP<<std::endl;


	      std::cout<<"TauVisMass_NoGsf :"<<TauVisMass_NoGsf<<std::endl;
	      std::cout<<"GammaEtaMom_NoGsf :"<<GammaEtaMom_NoGsf<<std::endl;
	      std::cout<<"GammaPhiMom_NoGsf :"<<GammaPhiMom_NoGsf<<std::endl;
	      std::cout<<"GammaEnFrac_NoGsf :"<<GammaEnFrac_NoGsf<<std::endl;
	      std::cout<<"TauVisMass_HasGsf :"<<TauVisMass_HasGsf<<std::endl;
	      std::cout<<"GammaEtaMom_HasGsf :"<<GammaEtaMom_HasGsf<<std::endl;
	      std::cout<<"GammaPhiMom_HasGsf :"<<GammaPhiMom_HasGsf<<std::endl;
	      std::cout<<"GammaEnFrac_HasGsf :"<<GammaEnFrac_HasGsf<<std::endl;
	      std::cout<<"TauLeadPFChargedHadrMva_HasGsf :"<<TauLeadPFChargedHadrMva_HasGsf<<std::endl;

	    }

	    hTauAbsEta_->Fill(TauAbsEta);
	    hTauPt_->Fill(TauPt);
	    hTauHasGsf_->Fill(TauHasGsf); 
	    hTauEmFraction_->Fill(TauEmFraction); 
	    hTauSignalPFChargedCands_->Fill(TauSignalPFChargedCands); 
	    hTauSignalPFGammaCands_->Fill(TauSignalPFGammaCands);

	    hTauLeadPFChargedHadrHoP_->Fill(TauLeadPFChargedHadrHoP); 
	    hTauLeadPFChargedHadrEoP_->Fill(TauLeadPFChargedHadrEoP); 
 
	    hTauVisMass_NoGsf_->Fill(TauVisMass_NoGsf); 
	    hGammaEtaMom_NoGsf_->Fill(GammaEtaMom_NoGsf); 
	    hGammaPhiMom_NoGsf_->Fill(GammaPhiMom_NoGsf); 
	    hGammaEnFrac_NoGsf_->Fill(GammaEnFrac_NoGsf); 
	    hTauVisMass_HasGsf_->Fill(TauVisMass_HasGsf); 
	    hGammaEtaMom_HasGsf_->Fill(GammaEtaMom_HasGsf); 
	    hGammaPhiMom_HasGsf_->Fill(GammaPhiMom_HasGsf); 
	    hGammaEnFrac_HasGsf_->Fill(GammaEnFrac_HasGsf); 
	    hTauLeadPFChargedHadrMva_HasGsf_->Fill(TauLeadPFChargedHadrMva_HasGsf); 

	  }

	}// loop PfTaus
	      
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
    TH1F* hHoHplusE_;
    TH1F* hFbrem_;
    TH1F* hChi2KF_;
    TH1F* hChi2GSF_;
    TH1F* hNHits_;
    TH1F* hGSFResol_;
    TH1F* hGSFlnPt_;
    TH1F* hGSFEta_;

    TH1F* hTauAbsEta_;
    TH1F* hTauPt_;
    TH1F* hTauEmFraction_;
    TH1F* hTauSignalPFChargedCands_;
    TH1F* hTauSignalPFGammaCands_;
    TH1F* hTauLeadPFChargedHadrHoP_;
    TH1F* hTauLeadPFChargedHadrEoP_;
    TH1F* hTauHasGsf_;
    TH1F* hTauVisMass_NoGsf_;
    TH1F* hGammaEtaMom_NoGsf_;
    TH1F* hGammaPhiMom_NoGsf_;
    TH1F* hGammaEnFrac_NoGsf_;
    TH1F* hTauVisMass_HasGsf_;
    TH1F* hGammaEtaMom_HasGsf_;
    TH1F* hGammaPhiMom_HasGsf_;
    TH1F* hGammaEnFrac_HasGsf_;
    TH1F* hTauLeadPFChargedHadrMva_HasGsf_;

  };//plotEntryType
  
  plotEntryType* plotsAll_;
  plotEntryType* plotsAll_Barrel_;
  plotEntryType* plotsAll_Endcap_;
  plotEntryType* plotsNumPV0to5_Barrel_Pt20to30_;
  plotEntryType* plotsNumPV5to10_Barrel_Pt20to30_;
  plotEntryType* plotsNumPV10to20_Barrel_Pt20to30_;
  plotEntryType* plotsNumPVOver20_Barrel_Pt20to30_;
  plotEntryType* plotsNumPV0to5_Endcap_Pt20to30_;
  plotEntryType* plotsNumPV5to10_Endcap_Pt20to30_;
  plotEntryType* plotsNumPV10to20_Endcap_Pt20to30_;
  plotEntryType* plotsNumPVOver20_Endcap_Pt20to30_;
  plotEntryType* plotsNumPV0to5_Barrel_Pt30to40_;
  plotEntryType* plotsNumPV5to10_Barrel_Pt30to40_;
  plotEntryType* plotsNumPV10to20_Barrel_Pt30to40_;
  plotEntryType* plotsNumPVOver20_Barrel_Pt30to40_;
  plotEntryType* plotsNumPV0to5_Endcap_Pt30to40_;
  plotEntryType* plotsNumPV5to10_Endcap_Pt30to40_;
  plotEntryType* plotsNumPV10to20_Endcap_Pt30to40_;
  plotEntryType* plotsNumPVOver20_Endcap_Pt30to40_;
  plotEntryType* plotsNumPV0to5_Barrel_Pt40to50_;
  plotEntryType* plotsNumPV5to10_Barrel_Pt40to50_;
  plotEntryType* plotsNumPV10to20_Barrel_Pt40to50_;
  plotEntryType* plotsNumPVOver20_Barrel_Pt40to50_;
  plotEntryType* plotsNumPV0to5_Endcap_Pt40to50_;
  plotEntryType* plotsNumPV5to10_Endcap_Pt40to50_;
  plotEntryType* plotsNumPV10to20_Endcap_Pt40to50_;
  plotEntryType* plotsNumPVOver20_Endcap_Pt40to50_;

  int numPV_;
};//AntiEMVAVariablesAnalyzer

#endif   

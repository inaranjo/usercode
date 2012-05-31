#include "IvoNaranjo/ElectronsStudies/interface/AntiEMVAAnalyzer2.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
// #include "DataFormats/PatCandidates/interface/Electron.h"


AntiEMVAAnalyzer2::AntiEMVAAnalyzer2(const edm::ParameterSet& cfg)
{ 
  srcPrimaryVertex_  = cfg.getParameter<edm::InputTag>("srcPrimaryVertex");
  srcGsfElectrons_  = cfg.getParameter<edm::InputTag>("srcGsfElectrons");
  srcPFTaus_  = cfg.getParameter<edm::InputTag>("srcPFTaus");
  srcGenElectrons_  = cfg.getParameter<edm::InputTag>("srcGenElectrons");
  srcGenElectronsFromZ_  = cfg.getParameter<edm::InputTag>("srcGenElectronsFromZ");
  srcGenElectronsFromZTauTau_  = cfg.getParameter<edm::InputTag>("srcGenElectronsFromZTauTau");
  srcGenTaus_  = cfg.getParameter<edm::InputTag>("srcGenTaus");
  srcGenJets_  = cfg.getParameter<edm::InputTag>("srcGenJets");
  srcPatTaus_  = cfg.getParameter<edm::InputTag>("srcPatTaus");
  debug_  = cfg.getParameter<bool>("debug");

} 

void AntiEMVAAnalyzer2::beginJob()
{ 
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","AntiEMVA tree");
  //counters
  tree_->Branch("run",&run_,"run/l");
  tree_->Branch("event",&event_,"event/l");
  tree_->Branch("lumi",&lumi_,"lumi/l");
  tree_->Branch("NumPV",&NumPV_,"NumPV/I");
  tree_->Branch("NumGsfEle",&NumGsfEle_,"NumGsfEle/I");
  tree_->Branch("NumGenEle",&NumGenEle_,"NumGenEle/I");
  tree_->Branch("NumPFTaus",&NumPFTaus_,"NumPFTaus/I");
  tree_->Branch("NumGenHad",&NumGenHad_,"NumGenHad/I");
  tree_->Branch("NumGenJet",&NumGenJet_,"NumGenEJet/I");

  //GsfElectron variables
  tree_->Branch("Elec_GenEleMatch",&Elec_GenEleMatch_,"Elec_GenEleMatch/I");
  tree_->Branch("Elec_GenEleFromZMatch",&Elec_GenEleFromZMatch_,"Elec_GenEleFromZMatch/I");
  tree_->Branch("Elec_GenEleFromZTauTauMatch",&Elec_GenEleFromZTauTauMatch_,"Elec_GenEleFromZTauTauMatch/I");
  tree_->Branch("Elec_GenHadMatch",&Elec_GenHadMatch_,"Elec_GenHadMatch/I");
  tree_->Branch("Elec_GenJetMatch",&Elec_GenJetMatch_,"Elec_GenJetMatch/I");
  tree_->Branch("Elec_AbsEta",&Elec_AbsEta_,"Elec_AbsEta/F");
  tree_->Branch("Elec_Pt",&Elec_Pt_,"Elec_Pt/F");
  tree_->Branch("Elec_HasSC",&Elec_HasSC_,"Elec_HasSC/F");
  tree_->Branch("Elec_HasKF",&Elec_HasKF_,"Elec_HasKF/F");
  tree_->Branch("Elec_HasGSF",&Elec_HasGSF_,"Elec_HasGSF/F");
  tree_->Branch("Elec_PFMvaOutput",&Elec_PFMvaOutput_,"Elec_PFMvaOutput/F");
  tree_->Branch("Elec_Ee",&Elec_Ee_,"Elec_Ee/F");
  tree_->Branch("Elec_Egamma",&Elec_Egamma_,"Elec_Egamma/F");
  tree_->Branch("Elec_Pin",&Elec_Pin_,"Elec_Pin/F");
  tree_->Branch("Elec_Pout",&Elec_Pout_,"Elec_Pout/F");
  tree_->Branch("Elec_EtotOverPin",&Elec_EtotOverPin_,"Elec_EtotOverPin/F");
  tree_->Branch("Elec_EeOverPout",&Elec_EeOverPout_,"Elec_EeOverPout/F");
  tree_->Branch("Elec_EgammaOverPdif",&Elec_EgammaOverPdif_,"Elec_EgammaOverPdif/F");
  tree_->Branch("Elec_EarlyBrem",&Elec_EarlyBrem_,"Elec_EarlyBrem/I");
  tree_->Branch("Elec_LateBrem",&Elec_LateBrem_,"Elec_LateBrem/I");
  tree_->Branch("Elec_Logsihih",&Elec_Logsihih_,"Elec_Logsihih/F");
  tree_->Branch("Elec_DeltaEta",&Elec_DeltaEta_,"Elec_DeltaEta/F");
  tree_->Branch("Elec_HoHplusE",&Elec_HoHplusE_,"Elec_HoHplusE/F");
  tree_->Branch("Elec_Fbrem",&Elec_Fbrem_,"Elec_Fbrem/F");
  tree_->Branch("Elec_Chi2KF",&Elec_Chi2KF_,"Elec_Chi2KF/F");
  tree_->Branch("Elec_Chi2GSF",&Elec_Chi2GSF_,"Elec_Chi2GSF/F");
  tree_->Branch("Elec_NumHits",&Elec_NumHits_,"Elec_NumHits/F");
  tree_->Branch("Elec_GSFTrackResol",&Elec_GSFTrackResol_,"Elec_GSFTrackResol/F");
  tree_->Branch("Elec_GSFTracklnPt",&Elec_GSFTracklnPt_,"Elec_GSFTracklnPt/F");
  tree_->Branch("Elec_GSFTrackEta",&Elec_GSFTrackEta_,"Elec_GSFTrackEta/F");
  //PFTaus variables
  tree_->Branch("Tau_GsfEleMatch",&Tau_GsfEleMatch_,"Tau_GsfEleMatch/I");
  tree_->Branch("Tau_GenEleMatch",&Tau_GenEleMatch_,"Tau_GenEleMatch/I");
  tree_->Branch("Tau_GenEleFromZMatch",&Tau_GenEleFromZMatch_,"Tau_GenEleFromZMatch/I");
  tree_->Branch("Tau_GenEleFromZTauTauMatch",&Tau_GenEleFromZTauTauMatch_,"Tau_GenEleFromZTauTauMatch/I");
  tree_->Branch("Tau_GenHadMatch",&Tau_GenHadMatch_,"Tau_GenHadMatch/I");
  tree_->Branch("Tau_GenJetMatch",&Tau_GenJetMatch_,"Tau_GenJetMatch/I");
  tree_->Branch("Tau_Eta",&Tau_Eta_,"Tau_Eta/F");
  tree_->Branch("Tau_EtaAtEcalEntrance",&Tau_EtaAtEcalEntrance_,"Tau_EtaAtEcalEntrance/F");
  tree_->Branch("Tau_Phi",&Tau_Phi_,"Tau_Phi/F");
  tree_->Branch("Tau_Pt",&Tau_Pt_,"Tau_Pt/F");
  tree_->Branch("Tau_LeadHadronPt",&Tau_LeadHadronPt_,"Tau_LeadHadronPt/F");
  tree_->Branch("Tau_HasGsf",&Tau_HasGsf_,"Tau_HasGsf/F");
  tree_->Branch("Tau_EmFraction",&Tau_EmFraction_,"Tau_EmFraction/F");
  tree_->Branch("Tau_NumChargedCands",&Tau_NumChargedCands_,"Tau_NumChargedCands/F");
  tree_->Branch("Tau_NumGammaCands",&Tau_NumGammaCands_,"Tau_NumGammaCands/F");
  tree_->Branch("Tau_HadrHoP",&Tau_HadrHoP_,"Tau_HadrHoP/F");
  tree_->Branch("Tau_HadrEoP",&Tau_HadrEoP_,"Tau_HadrEoP/F");
  tree_->Branch("Tau_VisMass",&Tau_VisMass_,"Tau_VisMass/F");
  tree_->Branch("Tau_GammaEtaMom",&Tau_GammaEtaMom_,"Tau_GammaEtaMom/F");
  tree_->Branch("Tau_GammaPhiMom",&Tau_GammaPhiMom_,"Tau_GammaPhiMom/F");
  tree_->Branch("Tau_GammaEnFrac",&Tau_GammaEnFrac_,"Tau_GammaEnFrac/F");
  tree_->Branch("Tau_HadrMva",&Tau_HadrMva_,"Tau_HadrMva/F");
  tree_->Branch("Tau_mvaAntiEValue",&Tau_mvaAntiEValue_,"Tau_mvaAntiEValue/F");
  tree_->Branch("Tau_AntiELoose",&Tau_AntiELoose_,"Tau_AntiELoose/F");
  tree_->Branch("Tau_AntiEMedium",&Tau_AntiEMedium_,"Tau_AntiEMedium/F");
  tree_->Branch("Tau_AntiETight",&Tau_AntiETight_,"Tau_AntiETight/F");
  tree_->Branch("Tau_AntiEMVA",&Tau_AntiEMVA_,"Tau_AntiEMVA/F");
  tree_->Branch("Tau_MatchElePassVeto",&Tau_MatchElePassVeto_,"Tau_MatchElePassVeto/F");
}


void AntiEMVAAnalyzer2::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  run_   = evt.run();
  event_ = (evt.eventAuxiliary()).event();
  lumi_  = evt.luminosityBlock();
  
  edm::Handle<reco::GsfElectronCollection> GsfElectrons;
  evt.getByLabel(srcGsfElectrons_, GsfElectrons);

  edm::Handle<reco::PFTauCollection> PfTaus;
  evt.getByLabel(srcPFTaus_, PfTaus);
  const reco::PFTauCollection* Pftaus = PfTaus.product();
  if(debug_){
//     std::cout<<"num PfTaus : "<<Pftaus->size()<<std::endl;
  }

  edm::Handle<pat::TauCollection> PatTaus;
  evt.getByLabel(srcPatTaus_,PatTaus);
  const pat::TauCollection* taus = PatTaus.product();
  if(debug_){
//     std::cout<<"num PatTaus : "<<taus->size()<<std::endl;
  }
  typedef edm::View<reco::Candidate> CandidateView;

  edm::Handle<CandidateView> GenElectrons;
  evt.getByLabel(srcGenElectrons_, GenElectrons);
  if( !GenElectrons.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No GenElectrons label available \n";

  edm::Handle<CandidateView> GenElectronsFromZ;
  evt.getByLabel(srcGenElectronsFromZ_, GenElectronsFromZ);
  if( !GenElectronsFromZ.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No GenElectronsFromZ label available \n";

  edm::Handle<CandidateView> GenElectronsFromZTauTau;
  evt.getByLabel(srcGenElectronsFromZTauTau_, GenElectronsFromZTauTau);
  if( !GenElectronsFromZTauTau.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No GenElectronsFromZTauTau label available \n";

  edm::Handle<CandidateView> GenTaus;
  evt.getByLabel(srcGenTaus_, GenTaus);
  if( !GenTaus.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No GenTaus label available \n";

  edm::Handle<CandidateView> GenJets;
  evt.getByLabel(srcGenJets_, GenJets);
  if( !GenJets.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No GenJets label available \n";

  edm::Handle<reco::VertexCollection> pVertexes;
  evt.getByLabel(srcPrimaryVertex_, pVertexes);
  if( !pVertexes.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No PV label available \n";
  const reco::VertexCollection* vertexes = pVertexes.product();
  if(debug_){
//     std::cout<<"numPV : "<<vertexes->size()<<std::endl;
//     for(unsigned int k = 0; k<vertexes->size(); k++){
//       std::cout << "Vtx[" << k << "] (x,y,z) = (" << ((*vertexes)[k].position()).x()
// 		<< "," << ((*vertexes)[k].position()).y() << "," << ((*vertexes)[k].position()).z() << ")"
// 		<< std::endl;
//     }
  }
  NumPV_ = vertexes->size();

  ///////////////////////////////////////////////////////////////////////////
  ////////////////////////////////PfTaus variables///////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  NumPFTaus_ = 0;
  NumGsfEle_ = 0;

  for ( reco::PFTauCollection::const_iterator PfTau = PfTaus->begin();
	PfTau != PfTaus->end(); ++PfTau ) {

    if (NumPFTaus_>49) continue; //Accept first 50 taus


    Tau_EtaAtEcalEntrance_ = -99;
    Tau_Eta_ = -99;
    Tau_Phi_ = -99;
    Tau_Pt_ = -99;
    Tau_LeadHadronPt_ = -99;
    Tau_HasGsf_ = -99;
    Tau_EmFraction_ = -99;
    Tau_NumChargedCands_ = -99;
    Tau_NumGammaCands_ = -99;
    Tau_HadrHoP_ = -99;
    Tau_HadrEoP_ = -99;
    Tau_VisMass_ = -99;
    Tau_GammaEtaMom_ = -99;
    Tau_GammaPhiMom_ = -99;
    Tau_GammaEnFrac_ = -99;
    Tau_HadrMva_ = -99;
    Tau_mvaAntiEValue_ = -99;
    Tau_AntiELoose_ = -99;
    Tau_AntiEMedium_ = -99;
    Tau_AntiETight_ = -99;
    Tau_AntiEMVA_ = -99;
    Tau_MatchElePassVeto_ = -99;


    /////////////////////////////////////Matchings  /////////////////////////////////////
    Tau_GsfEleMatch_= 0;
    Tau_GenEleMatch_ = 0;
    Tau_GenEleFromZMatch_ = 0;
    Tau_GenEleFromZTauTauMatch_ = 0;
    Tau_GenHadMatch_ = 0;
    Tau_GenJetMatch_ = 0;
    for ( reco::CandidateView::const_iterator GenElectron = GenElectrons->begin();
	  GenElectron != GenElectrons->end(); ++GenElectron ) {

      if(deltaR(PfTau->eta(),PfTau->phi(),GenElectron->eta(),GenElectron->phi())<0.3)Tau_GenEleMatch_ = 1;
    }
    for ( reco::CandidateView::const_iterator GenElectronFromZ = GenElectronsFromZ->begin();
	  GenElectronFromZ != GenElectronsFromZ->end(); ++GenElectronFromZ ) {

      if(deltaR(PfTau->eta(),PfTau->phi(),GenElectronFromZ->eta(),GenElectronFromZ->phi())<0.3) Tau_GenEleFromZMatch_ = 1;
    }
    for ( reco::CandidateView::const_iterator GenElectronFromZTauTau = GenElectronsFromZTauTau->begin();
	  GenElectronFromZTauTau != GenElectronsFromZTauTau->end(); ++GenElectronFromZTauTau ) {

      if(deltaR(PfTau->eta(),PfTau->phi(),GenElectronFromZTauTau->eta(),GenElectronFromZTauTau->phi())<0.3) Tau_GenEleFromZTauTauMatch_ = 1;
    }
    for ( reco::CandidateView::const_iterator GenTau = GenTaus->begin();
	  GenTau != GenTaus->end(); ++GenTau ) {

      if(deltaR(PfTau->eta(),PfTau->phi(),GenTau->eta(),GenTau->phi())<0.3)Tau_GenHadMatch_ = 1;
    }
    for ( reco::CandidateView::const_iterator GenJet = GenJets->begin();
	  GenJet != GenJets->end(); ++GenJet ) {

      if(deltaR(PfTau->eta(),PfTau->phi(),GenJet->eta(),GenJet->phi())<0.3)Tau_GenJetMatch_ = 1;
    }
    /////////////////////////////////////Matchings  /////////////////////////////////////
    float sumEtaTimesEnergy = 0;
    float sumEnergy = 0;
    for(unsigned int j = 0 ; j < (PfTau->signalPFCands()).size() ; j++){
       reco::PFCandidateRef pfcandidate = (PfTau->signalPFCands()).at(j);
       sumEtaTimesEnergy += pfcandidate->positionAtECALEntrance().eta()*pfcandidate->energy();
       sumEnergy += pfcandidate->energy();
    }
    if(sumEnergy>0)Tau_EtaAtEcalEntrance_ = sumEtaTimesEnergy/sumEnergy;

    Tau_Eta_ = PfTau->eta();
    Tau_Pt_ = PfTau->pt();
    Tau_Phi_ = PfTau->phi();
    Tau_EmFraction_ = TMath::Max(PfTau->emFraction(),float(0.0));
    Tau_NumChargedCands_ = PfTau->signalPFChargedHadrCands().size();
    Tau_NumGammaCands_  = PfTau->signalPFGammaCands().size();
    Tau_LeadHadronPt_ = PfTau->leadPFChargedHadrCand()->pt();
    Tau_HasGsf_ = (PfTau->leadPFChargedHadrCand()->gsfTrackRef()).isNonnull();
    Tau_HadrHoP_ = PfTau->leadPFChargedHadrCand()->hcalEnergy()/PfTau->leadPFChargedHadrCand()->p();
    Tau_HadrEoP_ = PfTau->leadPFChargedHadrCand()->ecalEnergy()/PfTau->leadPFChargedHadrCand()->p();

    GammasdEta.clear();
    GammasdPhi.clear();
    GammasPt.clear();

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

    Tau_GammaEtaMom_ = TMath::Sqrt(dEta2)*TMath::Sqrt(GammadPt_)*PfTau->pt();
    Tau_GammaPhiMom_ = TMath::Sqrt(dPhi2)*TMath::Sqrt(GammadPt_)*PfTau->pt();  
    Tau_GammaEnFrac_ = GammadPt_;
    Tau_VisMass_ = PfTau->mass();
    Tau_HadrMva_ = TMath::Max(PfTau->electronPreIDOutput(),float(-1.0));

    if(debug_){
//       std::cout<<"PfTauLoop number : "<<NumPFTaus_<<std::endl;
//       std::cout<<"  GammaEtaMom : "<<Tau_GammaEtaMom_<<std::endl;
//       std::cout<<"  GammaPhiMom : "<<Tau_GammaPhiMom_<<std::endl;
//       std::cout<<"  GammaPt : "<<sumPt<<std::endl;
//       std::cout<<"  TauPt : "<<PfTau->pt()<<std::endl;
//       std::cout<<"  TauEta : "<<PfTau->eta()<<std::endl;
//       std::cout<<"  TauEtaAtEcalEntrance : "<<Tau_EtaAtEcalEntrance_<<std::endl;
//       std::cout<<"  sumEtaTimesEnergy : "<<sumEtaTimesEnergy<<std::endl;
//       std::cout<<"  sumEnergy : "<<sumEnergy<<std::endl;
//       std::cout<<"  GammaEnFrac : "<<Tau_GammaEnFrac_<<std::endl;
    }
    NumPatTaus_ = 0;
    int countMatch = 0;
    for (pat::TauCollection::const_iterator  PatTau = PatTaus->begin();
	  PatTau != PatTaus->end(); ++PatTau ) {
      if(deltaR(PfTau->eta(),PfTau->phi(),PatTau->eta(),PatTau->phi())<0.01){
	countMatch++;
	Tau_AntiELoose_ = PatTau->tauID("againstElectronLoose");
	Tau_AntiEMedium_ = PatTau->tauID("againstElectronMedium");
	Tau_AntiETight_ = PatTau->tauID("againstElectronTight");
	Tau_AntiEMVA_ = PatTau->tauID("againstElectronMVA");

	if(debug_){
// 	  std::cout<<"===================>PAT MATCHED!!"<<std::endl;
// 	  std::cout<<" PatTauLoop number : "<<NumPatTaus_<<std::endl;
// 	  std::cout<<"  PfTau  : "<<PfTau->pt()<<" "<<PfTau->eta()<<" "<<PfTau->phi()<<std::endl;
// 	  std::cout<<"  PatTau : "<<PatTau->pt()<<" "<<PatTau->eta()<<" "<<PatTau->phi()<<std::endl;
// 	  std::cout<<" PatTau AntiELoose : "<<PatTau->tauID("againstElectronLoose")<<std::endl;
// 	  std::cout<<" PatTau AntiEMedium : "<<PatTau->tauID("againstElectronMedium")<<std::endl;
// 	  std::cout<<" PatTau AntiETight : "<<PatTau->tauID("againstElectronTight")<<std::endl;
// 	  std::cout<<" PatTau AntiEMVA : "<<PatTau->tauID("againstElectronMVA")<<std::endl;
// 	  std::cout<<" PatTau DecayMode : "<<PatTau->tauID("decayModeFinding")<<std::endl;
// 	  std::cout<<" PatTau CombIsoDB : "<<PatTau->tauID("byLooseCombinedIsolationDeltaBetaCorr")<<std::endl;
	}
      }
      NumPatTaus_++ ;   
    }
    
    //Tau is matched to a gsfElectron passing CutsVeto for SecondElectronVeto
    bool matchElectronCutsVeto = false;
    float matchElectronCutsVetoFloat = 0;
    for ( reco::GsfElectronCollection::const_iterator GsfElectron = GsfElectrons->begin();
	  GsfElectron != GsfElectrons->end(); ++GsfElectron ) {
      const reco::Track *el_track = (const reco::Track*)((GsfElectron)->gsfTrack().get());  
      const reco::HitPattern& p_inner = el_track->trackerExpectedHitsInner(); 
      float nHits = p_inner.numberOfHits();
      float dPhi  = fabs(GsfElectron->deltaPhiSuperClusterTrackAtVtx());
      float dEta  = fabs(GsfElectron->deltaEtaSuperClusterTrackAtVtx());
      float sihih = GsfElectron->sigmaIetaIeta();
      float HoE   = GsfElectron->hadronicOverEm();
      if(debug_){
	std::cout<<"GsfElectron nHits: "<<nHits<<std::endl;
	std::cout<<"GsfElectron dPhi: "<<dPhi<<std::endl;
	std::cout<<"GsfElectron dEta: "<<dEta<<std::endl;
	std::cout<<"GsfElectron sihih: "<<sihih<<std::endl;
	std::cout<<"GsfElectron HoE: "<<HoE<<std::endl;
      }
      bool ElectronPassCutsVeto = false;
      if((nHits<=999) &&
	 ((fabs(GsfElectron->eta())<1.5) &&
	  (sihih < 0.010) &&
	  (dPhi < 0.80) &&
	  (dEta < 0.007) &&
	  (HoE < 0.15)) ||
	 ((fabs(GsfElectron->eta())>1.5) && (fabs(GsfElectron->eta())<2.3) &&
	  (sihih < 0.030) &&
	  (dPhi < 0.70) &&
	  (dEta < 0.010) &&
	  (HoE < 999))
	 ) ElectronPassCutsVeto = true ;
      if ((deltaR(PfTau->eta(),PfTau->phi(),GsfElectron->eta(),GsfElectron->phi())<0.3) && ElectronPassCutsVeto){
	matchElectronCutsVeto = true;
      } 
    }//Loop on GsfElectrons
    if(matchElectronCutsVeto)matchElectronCutsVetoFloat = 1;
    Tau_MatchElePassVeto_ = matchElectronCutsVetoFloat;
    
    
    NumPFTaus_++;

    ///////////////////////////////////////////////////////////////////////////
    ////////////////////////////////GsfElectron variables//////////////////////
    ///////////////////////////////////////////////////////////////////////////
    for ( reco::GsfElectronCollection::const_iterator GsfElectron = GsfElectrons->begin();
	  GsfElectron != GsfElectrons->end(); ++GsfElectron ) {
      if (debug_){
	//       std::cout<<std::endl;
	//       std::cout<<"GsfElectron number : "<<NumGsfEle_<<std::endl;
      }
      if(deltaR(PfTau->eta(),PfTau->phi(),GsfElectron->eta(),GsfElectron->phi())<0.3)Tau_GsfEleMatch_ = 1;


    NumGenEle_ = 0;
    NumGenHad_ = 0;
    NumGenJet_ = 0;

    if (Tau_GsfEleMatch_ != 1)continue;
    if(GsfElectron->pt()<10)continue;

    Elec_AbsEta_ = -99;
    Elec_Pt_ = -99;
    Elec_PFMvaOutput_ = -99;
    Elec_EarlyBrem_ =  -99;
    Elec_LateBrem_=  -99;
    Elec_Logsihih_ =  -99;
    Elec_DeltaEta_ = -99;
    Elec_Fbrem_ =  -99;

    //Variables related to the SC
    Elec_HasSC_ = -99;
    Elec_Ee_ = -99;
    Elec_Egamma_ = -99;
    Elec_Pin_ = -99;
    Elec_Pout_ = -99;
    Elec_EtotOverPin_ = -99;
    Elec_EeOverPout_ = -99;
    Elec_EgammaOverPdif_ = -99;
    Elec_HoHplusE_ = -99;

    Elec_HasKF_ = -99;
    Elec_Chi2KF_ = -99;
    Elec_NumHits_ = -99;

    Elec_HasGSF_ = -99;
    Elec_Chi2GSF_ = -99;
    Elec_GSFTrackResol_ = -99;
    Elec_GSFTracklnPt_ = -99;
    Elec_GSFTrackEta_ = -99;

    /////////////////////////////////////Matchings  /////////////////////////////////////
    Elec_GenEleMatch_ = 0;
    Elec_GenEleFromZMatch_ = 0;
    Elec_GenEleFromZTauTauMatch_ = 0;
    Elec_GenHadMatch_ = 0; 
    Elec_GenJetMatch_ = 0;
    for ( reco::CandidateView::const_iterator GenElectron = GenElectrons->begin();
	  GenElectron != GenElectrons->end(); ++GenElectron ) {
//       if(debug_){
// 	std::cout<<"  DeltaR GsfEle-GenEle :"<<deltaR(GsfElectron->eta(),GsfElectron->phi(),GenElectron->eta(),GenElectron->phi())<<std::endl;
//       }
      NumGenEle_++;
      if(deltaR(GsfElectron->eta(),GsfElectron->phi(),GenElectron->eta(),GenElectron->phi())<0.3) Elec_GenEleMatch_ = 1;
    }
    for ( reco::CandidateView::const_iterator GenElectronFromZ = GenElectronsFromZ->begin();
	  GenElectronFromZ != GenElectronsFromZ->end(); ++GenElectronFromZ ) {
//       if(debug_){
// 	std::cout<<"  DeltaR GsfEle-GenEleFromZ :"<<deltaR(GsfElectron->eta(),GsfElectron->phi(),GenElectronFromZ->eta(),GenElectronFromZ->phi())<<std::endl;
//       }
      if(deltaR(GsfElectron->eta(),GsfElectron->phi(),GenElectronFromZ->eta(),GenElectronFromZ->phi())<0.3) Elec_GenEleFromZMatch_ = 1;
    }
    for ( reco::CandidateView::const_iterator GenElectronFromZTauTau = GenElectronsFromZTauTau->begin();
	  GenElectronFromZTauTau != GenElectronsFromZTauTau->end(); ++GenElectronFromZTauTau ) {
//       if(debug_){
// 	std::cout<<"  DeltaR GsfEle-GenEleFromZTauTau :"<<deltaR(GsfElectron->eta(),GsfElectron->phi(),GenElectronFromZTauTau->eta(),GenElectronFromZTauTau->phi())<<std::endl;
//       }
      if(deltaR(GsfElectron->eta(),GsfElectron->phi(),GenElectronFromZTauTau->eta(),GenElectronFromZTauTau->phi())<0.3) Elec_GenEleFromZTauTauMatch_ = 1;
    }
    for ( reco::CandidateView::const_iterator GenTau = GenTaus->begin();
	  GenTau != GenTaus->end(); ++GenTau ) {
//       if(debug_){
// 	std::cout<<"  DeltaR GsfEle-GenTau :"<<deltaR(GsfElectron->eta(),GsfElectron->phi(),GenTau->eta(),GenTau->phi())<<std::endl;
//       }
      NumGenHad_++;
      if(deltaR(GsfElectron->eta(),GsfElectron->phi(),GenTau->eta(),GenTau->phi())<0.3)Elec_GenHadMatch_ = 1;
    }
    for ( reco::CandidateView::const_iterator GenJet = GenJets->begin();
	  GenJet != GenJets->end(); ++GenJet ) {
//       if(debug_){
// 	std::cout<<"  DeltaR GsfEle-GenJet :"<<deltaR(GsfElectron->eta(),GsfElectron->phi(),GenJet->eta(),GenJet->phi())<<std::endl;
//       }
      NumGenJet_++;
      if(deltaR(GsfElectron->eta(),GsfElectron->phi(),GenJet->eta(),GenJet->phi())<0.3)Elec_GenJetMatch_ = 1;
    }
    /////////////////////////////////////Matchings  /////////////////////////////////////

    Elec_AbsEta_ = TMath::Abs(GsfElectron->eta());
    Elec_Pt_ = GsfElectron->pt();
    Elec_PFMvaOutput_ = TMath::Max(GsfElectron->mvaOutput().mva,float(-1.0));
    Elec_EarlyBrem_ = GsfElectron->mvaInput().earlyBrem;
    Elec_LateBrem_= GsfElectron->mvaInput().lateBrem;
    Elec_Logsihih_ = log(GsfElectron->mvaInput().sigmaEtaEta);
    Elec_DeltaEta_ = GsfElectron->mvaInput().deltaEta;
    Elec_Fbrem_ = GsfElectron->fbrem();


    //Variables related to the SC
    Elec_HasSC_ = 0;
    Elec_Ee_ = -99;
    Elec_Egamma_ = -99;
    Elec_Pin_ = -99;
    Elec_Pout_ = -99;
    Elec_EtotOverPin_ = -99;
    Elec_EeOverPout_ = -99;
    Elec_EgammaOverPdif_ = -99;
    Elec_HoHplusE_ = -99;
    reco::SuperClusterRef pfSuperCluster = GsfElectron->pflowSuperCluster();
    if(pfSuperCluster.isNonnull() && pfSuperCluster.isAvailable()){
      Elec_HasSC_ = 1;
      Elec_Ee_ = 0.;
      Elec_Egamma_ = 0.;
//       if (debug_)std::cout<<"SuperCluster accessed   "<<std::endl;
      for (reco::CaloCluster_iterator pfCluster = pfSuperCluster->clustersBegin();
	   pfCluster != pfSuperCluster->clustersEnd(); ++pfCluster ) {
	float pfClusterEn = (*pfCluster)->energy();
	if ( pfCluster == pfSuperCluster->clustersBegin() ) Elec_Ee_ += pfClusterEn;
	else Elec_Egamma_ += pfClusterEn;
      }
      Elec_Pin_ = TMath::Sqrt(GsfElectron->trackMomentumAtVtx().Mag2());
      Elec_Pout_ = TMath::Sqrt(GsfElectron->trackMomentumOut().Mag2()); 
      Elec_EtotOverPin_ = (Elec_Ee_+Elec_Egamma_)/Elec_Pin_;
      Elec_EeOverPout_ = Elec_Ee_/Elec_Pout_;
      Elec_EgammaOverPdif_ = Elec_Egamma_/(Elec_Pin_-Elec_Pout_);
      Elec_HoHplusE_ = GsfElectron->mvaInput().hadEnergy/(GsfElectron->mvaInput().hadEnergy+Elec_Ee_) ;
      
    }

//     if(debug_)std::cout<<"Elec_Ee :   "<<Elec_Ee_<<" Elec_Egamma :  "<<Elec_Egamma_<<std::endl;
//     if (debug_)std::cout<<"Elec_Pin :   "<<Elec_Pin_<<" Elec_Pout :  "<<Elec_Pout_<<std::endl;
    
    //Variables related to the CtfTrack
    Elec_HasKF_ = 0;
    Elec_Chi2KF_ = -99;
    Elec_NumHits_ = -99;
    if (GsfElectron->closestCtfTrackRef().isNonnull()){
      Elec_HasKF_ = 1.;
      Elec_Chi2KF_ = GsfElectron->closestCtfTrackRef()->normalizedChi2();
      Elec_NumHits_ = GsfElectron->closestCtfTrackRef()->numberOfValidHits();
    }

    //Variables related to the GsfTrack
    Elec_HasGSF_ = 0;
    Elec_Chi2GSF_ = -99;
    Elec_GSFTrackResol_ = -99;
    Elec_GSFTracklnPt_ = -99;
    Elec_GSFTrackEta_ = -99;
    if(GsfElectron->gsfTrack().isNonnull()){
      Elec_HasGSF_ = 1.;
      Elec_Chi2GSF_ = GsfElectron->gsfTrack()->normalizedChi2();
      Elec_GSFTrackResol_ = GsfElectron->gsfTrack()->ptError()/GsfElectron->gsfTrack()->pt();
      Elec_GSFTracklnPt_ = log(GsfElectron->gsfTrack()->pt())*TMath::Ln10();
      Elec_GSFTrackEta_ = GsfElectron->gsfTrack()->eta();
    }
    
    if(debug_){
//       std::cout<<"Elec_AbsEta :"<<Elec_AbsEta_<<std::endl;
//       std::cout<<"Elec_Pt :"<<Elec_Pt_<<std::endl;
//       std::cout<<"Elec_HasSC :"<<Elec_HasSC_<<std::endl;
//       std::cout<<"Elec_HasKF :"<<Elec_HasKF_<<std::endl;
//       std::cout<<"Elec_HasGSF :"<<Elec_HasGSF_<<std::endl;
//       std::cout<<"Elec PFMvaOutput :"<<Elec_PFMvaOutput_<<std::endl;
//       std::cout<<"E electron cluster plus photons over Pin :   "<<Elec_EtotOverPin_<<std::endl;
//       std::cout<<"E electron cluster over Pout :   "<<Elec_EeOverPout_<<std::endl;
//       std::cout<<"E photons over (Elec_Pin - Elec_Pout) :   "<<Elec_EgammaOverPdif_<<std::endl;
//       std::cout<<"EarlyBrem :   "<<Elec_EarlyBrem_<<std::endl;
//       std::cout<<"LateBrem :   "<<Elec_LateBrem_<<std::endl;
//       std::cout<<"log(sigma EtaEta with the SC) :  "<<Elec_Logsihih_<<std::endl;
//       std::cout<<"PF-cluster GSF track delta-eta :  "<<Elec_DeltaEta_<<std::endl;
//       std::cout<<"H over H plus E :   "<<Elec_HoHplusE_<<std::endl;
//       std::cout<<"FBrem :   "<<Elec_Fbrem_<<std::endl;
//       std::cout<<"Normalized Chi2 KF :   "<<Elec_Chi2KF_<<std::endl;
//       std::cout<<"Normalized Chi2 GSF :   "<<Elec_Chi2GSF_<<std::endl;
//       std::cout<<"Number of valid hits KF :   "<<Elec_NumHits_<<std::endl;
//       std::cout<<"GSFTrack sig(Pt)/Pt :   "<<Elec_GSFTrackResol_<<std::endl;
//       std::cout<<"GSFTrack ln(Pt) :   "<<Elec_GSFTracklnPt_<<std::endl;
//       std::cout<<"GSFTrack Eta :   "<<Elec_GSFTrackEta_<<std::endl;
//       std::cout<<std::endl;
    }
    

  }//Loop on GsfElectrons
  tree_->Fill();
  
  
  }//Loop on PFTaus
    

}//analyze()

AntiEMVAAnalyzer2::~AntiEMVAAnalyzer2()
{

}
void AntiEMVAAnalyzer2::endJob()
{
// nothing to be done yet...
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(AntiEMVAAnalyzer2);




#include "IvoNaranjo/ElectronsStudies/interface/AntiEMVAAnalyzer.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


AntiEMVAAnalyzer::AntiEMVAAnalyzer(const edm::ParameterSet& cfg)
{ 
  srcPrimaryVertex_  = cfg.getParameter<edm::InputTag>("srcPrimaryVertex");
  srcGsfElectrons_  = cfg.getParameter<edm::InputTag>("srcGsfElectrons");
  srcPFTaus_  = cfg.getParameter<edm::InputTag>("srcPFTaus");
  srcGenElectrons_  = cfg.getParameter<edm::InputTag>("srcGenElectrons");
  srcGenElectronsFromZ_  = cfg.getParameter<edm::InputTag>("srcGenElectronsFromZ");
  srcGenElectronsFromZTauTau_  = cfg.getParameter<edm::InputTag>("srcGenElectronsFromZTauTau");
  srcGenTaus_  = cfg.getParameter<edm::InputTag>("srcGenTaus");
  srcGenJets_  = cfg.getParameter<edm::InputTag>("srcGenJets");
  debug_  = cfg.getParameter<bool>("debug");
} 



void AntiEMVAAnalyzer::beginJob()
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
  tree_->Branch("Elec_GenEleMatch",&Elec_GenEleMatch_,"Elec_GenEleMatch[NumGsfEle]/I");
  tree_->Branch("Elec_GenEleFromZMatch",&Elec_GenEleFromZMatch_,"Elec_GenEleFromZMatch[NumGsfEle]/I");
  tree_->Branch("Elec_GenEleFromZTauTauMatch",&Elec_GenEleFromZTauTauMatch_,"Elec_GenEleFromZTauTauMatch[NumGsfEle]/I");
  tree_->Branch("Elec_PFTauMatch",&Elec_PFTauMatch_,"Elec_PFTauMatch[NumGsfEle]/I");
  tree_->Branch("Elec_GenHadMatch",&Elec_GenHadMatch_,"Elec_GenHadMatch[NumGsfEle]/I");
  tree_->Branch("Elec_GenJetMatch",&Elec_GenJetMatch_,"Elec_GenJetMatch[NumGsfEle]/I");
  tree_->Branch("Elec_AbsEta",&Elec_AbsEta_,"Elec_AbsEta[NumGsfEle]/F");
  tree_->Branch("Elec_Pt",&Elec_Pt_,"Elec_Pt[NumGsfEle]/F");
  tree_->Branch("Elec_PFMvaOutput",&Elec_PFMvaOutput_,"Elec_PFMvaOutput[NumGsfEle]/F");
  tree_->Branch("Elec_Ee",&Elec_Ee_,"Elec_Ee[NumGsfEle]/F");
  tree_->Branch("Elec_Egamma",&Elec_Egamma_,"Elec_Egamma[NumGsfEle]/F");
  tree_->Branch("Elec_Pin",&Elec_Pin_,"Elec_Pin[NumGsfEle]/F");
  tree_->Branch("Elec_Pout",&Elec_Pout_,"Elec_Pout[NumGsfEle]/F");
  tree_->Branch("Elec_EtotOverPin",&Elec_EtotOverPin_,"Elec_EtotOverPin[NumGsfEle]/F");
  tree_->Branch("Elec_EeOverPout",&Elec_EeOverPout_,"Elec_EeOverPout[NumGsfEle]/F");
  tree_->Branch("Elec_EgammaOverPdif",&Elec_EgammaOverPdif_,"Elec_EgammaOverPdif[NumGsfEle]/F");
  tree_->Branch("Elec_EarlyBrem",&Elec_EarlyBrem_,"Elec_EarlyBrem[NumGsfEle]/I");
  tree_->Branch("Elec_LateBrem",&Elec_LateBrem_,"Elec_LateBrem[NumGsfEle]/I");
  tree_->Branch("Elec_Logsihih",&Elec_Logsihih_,"Elec_Logsihih[NumGsfEle]/F");
  tree_->Branch("Elec_DeltaEta",&Elec_DeltaEta_,"Elec_DeltaEta[NumGsfEle]/F");
  tree_->Branch("Elec_HoHplusE",&Elec_HoHplusE_,"Elec_HoHplusE[NumGsfEle]/F");
  tree_->Branch("Elec_Fbrem",&Elec_Fbrem_,"Elec_Fbrem[NumGsfEle]/F");
  tree_->Branch("Elec_Chi2KF",&Elec_Chi2KF_,"Elec_Chi2KF[NumGsfEle]/F");
  tree_->Branch("Elec_Chi2GSF",&Elec_Chi2GSF_,"Elec_Chi2GSF[NumGsfEle]/F");
  tree_->Branch("Elec_NumHits",&Elec_NumHits_,"Elec_NumHits[NumGsfEle]/F");
  tree_->Branch("Elec_GSFTrackResol",&Elec_GSFTrackResol_,"Elec_GSFTrackResol[NumGsfEle]/F");
  tree_->Branch("Elec_GSFTracklnPt",&Elec_GSFTracklnPt_,"Elec_GSFTracklnPt[NumGsfEle]/F");
  tree_->Branch("Elec_GSFTrackEta",&Elec_GSFTrackEta_,"Elec_GSFTrackEta[NumGsfEle]/F");
  //PFTaus variables
  tree_->Branch("Tau_GsfEleMatch",&Tau_GsfEleMatch_,"Tau_GsfEleMatch[NumPFTaus]/I");
  tree_->Branch("Tau_GenEleMatch",&Tau_GenEleMatch_,"Tau_GenEleMatch[NumPFTaus]/I");
  tree_->Branch("Tau_GenEleFromZMatch",&Tau_GenEleFromZMatch_,"Tau_GenEleFromZMatch[NumPFTaus]/I");
  tree_->Branch("Tau_GenEleFromZTauTauMatch",&Tau_GenEleFromZTauTauMatch_,"Tau_GenEleFromZTauTauMatch[NumPFTaus]/I");
  tree_->Branch("Tau_GenHadMatch",&Tau_GenHadMatch_,"Tau_GenHadMatch[NumPFTaus]/I");
  tree_->Branch("Tau_GenJetMatch",&Tau_GenJetMatch_,"Tau_GenJetMatch[NumPFTaus]/I");
  tree_->Branch("Tau_AbsEta",&Tau_AbsEta_,"Tau_AbsEta[NumPFTaus]/F");
  tree_->Branch("Tau_Pt",&Tau_Pt_,"Tau_Pt[NumPFTaus]/F");
  tree_->Branch("Tau_HasGsf",&Tau_HasGsf_,"Tau_HasGsf[NumPFTaus]/F");
  tree_->Branch("Tau_EmFraction",&Tau_EmFraction_,"Tau_EmFraction[NumPFTaus]/F");
  tree_->Branch("Tau_NumChargedCands",&Tau_NumChargedCands_,"Tau_NumChargedCands[NumPFTaus]/F");
  tree_->Branch("Tau_NumGammaCands",&Tau_NumGammaCands_,"Tau_NumGammaCands[NumPFTaus]/F");
  tree_->Branch("Tau_HadrHoP",&Tau_HadrHoP_,"Tau_HadrHoP[NumPFTaus]/F");
  tree_->Branch("Tau_HadrEoP",&Tau_HadrEoP_,"Tau_HadrEoP[NumPFTaus]/F");
  tree_->Branch("Tau_VisMass_",&Tau_VisMass_,"Tau_VisMass_[NumPFTaus]/F");
  tree_->Branch("Tau_GammaEtaMom_",&Tau_GammaEtaMom_,"Tau_GammaEtaMom_[NumPFTaus]/F");
  tree_->Branch("Tau_GammaPhiMom_",&Tau_GammaPhiMom_,"Tau_GammaPhiMom_[NumPFTaus]/F");
  tree_->Branch("Tau_GammaEnFrac_",&Tau_GammaEnFrac_,"Tau_GammaEnFrac_[NumPFTaus]/F");
  tree_->Branch("Tau_HadrMva_",&Tau_HadrMva_,"Tau_HadrMva_[NumPFTaus]/F");

}


void AntiEMVAAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  run_   = evt.run();
  event_ = (evt.eventAuxiliary()).event();
  lumi_  = evt.luminosityBlock();
  
  edm::Handle<reco::GsfElectronCollection> GsfElectrons;
  evt.getByLabel(srcGsfElectrons_, GsfElectrons);

  edm::Handle<reco::PFTauCollection> PfTaus;
  evt.getByLabel(srcPFTaus_, PfTaus);

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
    std::cout<<"numPV : "<<vertexes->size()<<std::endl;
    for(unsigned int k = 0; k<vertexes->size(); k++){
      std::cout << "Vtx[" << k << "] (x,y,z) = (" << ((*vertexes)[k].position()).x()
		<< "," << ((*vertexes)[k].position()).y() << "," << ((*vertexes)[k].position()).z() << ")"
		<< std::endl;
    }
  }
  NumPV_ = vertexes->size();

  ///////////////////////////////////////////////////////////////////////////
  ////////////////////////////////GsfElectron variables//////////////////////
  ///////////////////////////////////////////////////////////////////////////
  NumGsfEle_ = 0;

  for ( reco::GsfElectronCollection::const_iterator GsfElectron = GsfElectrons->begin();
	GsfElectron != GsfElectrons->end(); ++GsfElectron ) {
    if (debug_){
      std::cout<<std::endl;
      std::cout<<"GsfElectron number : "<<NumGsfEle_<<std::endl;
    }
    if (NumGsfEle_>49) continue; //Accept first 50 electrons

    NumGenEle_ = 0;
    NumGenHad_ = 0;
    NumGenJet_ = 0;


    /////////////////////////////////////Matchings  /////////////////////////////////////
    Elec_GenEleMatch_[NumGsfEle_] = 0;
    Elec_GenEleFromZMatch_[NumGsfEle_] = 0;
    Elec_GenEleFromZTauTauMatch_[NumGsfEle_] = 0;
    Elec_PFTauMatch_[NumGsfEle_] = 0;
    Elec_GenHadMatch_[NumGsfEle_] = 0;
    Elec_GenJetMatch_[NumGsfEle_] = 0;
    for ( reco::CandidateView::const_iterator GenElectron = GenElectrons->begin();
	  GenElectron != GenElectrons->end(); ++GenElectron ) {
      if(debug_){
	std::cout<<"  DeltaR GsfEle-GenEle :"<<deltaR(GsfElectron->eta(),GsfElectron->phi(),GenElectron->eta(),GenElectron->phi())<<std::endl;
      }
      NumGenEle_++;
      if(deltaR(GsfElectron->eta(),GsfElectron->phi(),GenElectron->eta(),GenElectron->phi())<0.3) Elec_GenEleMatch_[NumGsfEle_] = 1;
    }
    for ( reco::CandidateView::const_iterator GenElectronFromZ = GenElectronsFromZ->begin();
	  GenElectronFromZ != GenElectronsFromZ->end(); ++GenElectronFromZ ) {
      if(debug_){
	std::cout<<"  DeltaR GsfEle-GenEleFromZ :"<<deltaR(GsfElectron->eta(),GsfElectron->phi(),GenElectronFromZ->eta(),GenElectronFromZ->phi())<<std::endl;
      }
      if(deltaR(GsfElectron->eta(),GsfElectron->phi(),GenElectronFromZ->eta(),GenElectronFromZ->phi())<0.3) Elec_GenEleFromZMatch_[NumGsfEle_] = 1;
    }
    for ( reco::CandidateView::const_iterator GenElectronFromZTauTau = GenElectronsFromZTauTau->begin();
	  GenElectronFromZTauTau != GenElectronsFromZTauTau->end(); ++GenElectronFromZTauTau ) {
      if(debug_){
	std::cout<<"  DeltaR GsfEle-GenEleFromZTauTau :"<<deltaR(GsfElectron->eta(),GsfElectron->phi(),GenElectronFromZTauTau->eta(),GenElectronFromZTauTau->phi())<<std::endl;
      }
      if(deltaR(GsfElectron->eta(),GsfElectron->phi(),GenElectronFromZTauTau->eta(),GenElectronFromZTauTau->phi())<0.3) Elec_GenEleFromZTauTauMatch_[NumGsfEle_] = 1;
    }
    for ( reco::PFTauCollection::const_iterator PfTau = PfTaus->begin();
	  PfTau != PfTaus->end(); ++PfTau ) {
      if(debug_){
	std::cout<<"  DeltaR GsfEle-PfTau :"<<deltaR(GsfElectron->eta(),GsfElectron->phi(),PfTau->eta(),PfTau->phi())<<std::endl;
      }
      if(deltaR(GsfElectron->eta(),GsfElectron->phi(),PfTau->eta(),PfTau->phi())<0.3)Elec_PFTauMatch_[NumGsfEle_] = 1;
    }
    for ( reco::CandidateView::const_iterator GenTau = GenTaus->begin();
	  GenTau != GenTaus->end(); ++GenTau ) {
      if(debug_){
	std::cout<<"  DeltaR GsfEle-GenTau :"<<deltaR(GsfElectron->eta(),GsfElectron->phi(),GenTau->eta(),GenTau->phi())<<std::endl;
      }
      NumGenHad_++;
      if(deltaR(GsfElectron->eta(),GsfElectron->phi(),GenTau->eta(),GenTau->phi())<0.3)Elec_GenHadMatch_[NumGsfEle_] = 1;
    }
    for ( reco::CandidateView::const_iterator GenJet = GenJets->begin();
	  GenJet != GenJets->end(); ++GenJet ) {
      if(debug_){
	std::cout<<"  DeltaR GsfEle-GenJet :"<<deltaR(GsfElectron->eta(),GsfElectron->phi(),GenJet->eta(),GenJet->phi())<<std::endl;
      }
      NumGenJet_++;
      if(deltaR(GsfElectron->eta(),GsfElectron->phi(),GenJet->eta(),GenJet->phi())<0.3)Elec_GenJetMatch_[NumGsfEle_] = 1;
    }
    /////////////////////////////////////Matchings  /////////////////////////////////////

    Elec_AbsEta_[NumGsfEle_] = TMath::Abs(GsfElectron->eta());
    Elec_Pt_[NumGsfEle_] = GsfElectron->pt();
    Elec_PFMvaOutput_[NumGsfEle_] = TMath::Max(GsfElectron->mvaOutput().mva,float(-1.0));

    Elec_Ee_[NumGsfEle_] = -99;
    Elec_Egamma_[NumGsfEle_] = -99;
    reco::SuperClusterRef pfSuperCluster = GsfElectron->pflowSuperCluster();
    if(pfSuperCluster.isNonnull() && pfSuperCluster.isAvailable()){
      Elec_Ee_[NumGsfEle_] = 0.;
      Elec_Egamma_[NumGsfEle_] = 0.;
      if (debug_)std::cout<<"SuperCluster accessed   "<<std::endl;
      for (reco::CaloCluster_iterator pfCluster = pfSuperCluster->clustersBegin();
	   pfCluster != pfSuperCluster->clustersEnd(); ++pfCluster ) {
	float pfClusterEn = (*pfCluster)->energy();
	if ( pfCluster == pfSuperCluster->clustersBegin() ) Elec_Ee_[NumGsfEle_] += pfClusterEn;
	else Elec_Egamma_[NumGsfEle_] += pfClusterEn;
      }
    }
    
    if(debug_)std::cout<<"Elec_Ee :   "<<Elec_Ee_[NumGsfEle_]<<" Elec_Egamma :  "<<Elec_Egamma_[NumGsfEle_]<<std::endl;
	    
    Elec_Pin_[NumGsfEle_] = TMath::Sqrt(GsfElectron->trackMomentumAtVtx().Mag2());
    Elec_Pout_[NumGsfEle_] = TMath::Sqrt(GsfElectron->trackMomentumOut().Mag2()); 
    if (debug_)std::cout<<"Elec_Pin :   "<<Elec_Pin_<<" Elec_Pout :  "<<Elec_Pout_<<std::endl;
    
    Elec_EtotOverPin_[NumGsfEle_] = (Elec_Ee_[NumGsfEle_]+Elec_Egamma_[NumGsfEle_])/Elec_Pin_[NumGsfEle_];
    Elec_EeOverPout_[NumGsfEle_] = Elec_Ee_[NumGsfEle_]/Elec_Pout_[NumGsfEle_];
    Elec_EgammaOverPdif_[NumGsfEle_] = Elec_Egamma_[NumGsfEle_]/(Elec_Pin_[NumGsfEle_]-Elec_Pout_[NumGsfEle_]);
    Elec_EarlyBrem_[NumGsfEle_] = GsfElectron->mvaInput().earlyBrem;
    Elec_LateBrem_[NumGsfEle_]= GsfElectron->mvaInput().lateBrem;
    Elec_Logsihih_[NumGsfEle_] = log(GsfElectron->mvaInput().sigmaEtaEta);
    Elec_DeltaEta_[NumGsfEle_] = GsfElectron->mvaInput().deltaEta;
    Elec_HoHplusE_[NumGsfEle_] = GsfElectron->mvaInput().hadEnergy/(GsfElectron->mvaInput().hadEnergy+Elec_Ee_[NumGsfEle_]) ;
    Elec_Fbrem_[NumGsfEle_] = GsfElectron->fbrem();
    Elec_Chi2KF_[NumGsfEle_] = -99;
    Elec_Chi2GSF_[NumGsfEle_] = -99;
    Elec_NumHits_[NumGsfEle_] = -99;
    if (GsfElectron->closestCtfTrackRef().isNonnull()){
      Elec_Chi2KF_[NumGsfEle_] = GsfElectron->closestCtfTrackRef()->normalizedChi2();
      Elec_NumHits_[NumGsfEle_] = GsfElectron->closestCtfTrackRef()->numberOfValidHits();
    }
    Elec_Chi2GSF_[NumGsfEle_] = -99;
    Elec_GSFTrackResol_[NumGsfEle_] = -99;
    Elec_GSFTracklnPt_[NumGsfEle_] = -99;
    Elec_GSFTrackEta_[NumGsfEle_] = -99;
    if(GsfElectron->gsfTrack().isNonnull()){
      Elec_Chi2GSF_[NumGsfEle_] = GsfElectron->gsfTrack()->normalizedChi2();
      Elec_GSFTrackResol_[NumGsfEle_] = GsfElectron->gsfTrack()->ptError()/GsfElectron->gsfTrack()->pt();
      Elec_GSFTracklnPt_[NumGsfEle_] = log(GsfElectron->gsfTrack()->pt())*TMath::Ln10();
      Elec_GSFTrackEta_[NumGsfEle_] = GsfElectron->gsfTrack()->eta();
    }
    
    if(debug_){
      std::cout<<"Elec_AbsEta :"<<Elec_AbsEta_[NumGsfEle_]<<std::endl;
      std::cout<<"Elec_Pt :"<<Elec_Pt_[NumGsfEle_]<<std::endl;
      std::cout<<"Elec PFMvaOutput :"<<Elec_PFMvaOutput_[NumGsfEle_]<<std::endl;
      std::cout<<"E electron cluster plus photons over Pin :   "<<Elec_EtotOverPin_[NumGsfEle_]<<std::endl;
      std::cout<<"E electron cluster over Pout :   "<<Elec_EeOverPout_[NumGsfEle_]<<std::endl;
      std::cout<<"E photons over (Elec_Pin - Elec_Pout) :   "<<Elec_EgammaOverPdif_[NumGsfEle_]<<std::endl;
      std::cout<<"EarlyBrem :   "<<Elec_EarlyBrem_[NumGsfEle_]<<std::endl;
      std::cout<<"LateBrem :   "<<Elec_LateBrem_[NumGsfEle_]<<std::endl;
      std::cout<<"log(sigma EtaEta with the SC) :  "<<Elec_Logsihih_[NumGsfEle_]<<std::endl;
      std::cout<<"PF-cluster GSF track delta-eta :  "<<Elec_DeltaEta_[NumGsfEle_]<<std::endl;
      std::cout<<"H over H plus E :   "<<Elec_HoHplusE_[NumGsfEle_]<<std::endl;
      std::cout<<"FBrem :   "<<Elec_Fbrem_[NumGsfEle_]<<std::endl;
      std::cout<<"Normalized Chi2 KF :   "<<Elec_Chi2KF_[NumGsfEle_]<<std::endl;
      std::cout<<"Normalized Chi2 GSF :   "<<Elec_Chi2GSF_[NumGsfEle_]<<std::endl;
      std::cout<<"Number of valid hits KF :   "<<Elec_NumHits_[NumGsfEle_]<<std::endl;
      std::cout<<"GSFTrack sig(Pt)/Pt :   "<<Elec_GSFTrackResol_[NumGsfEle_]<<std::endl;
      std::cout<<"GSFTrack ln(Pt) :   "<<Elec_GSFTracklnPt_[NumGsfEle_]<<std::endl;
      std::cout<<"GSFTrack Eta :   "<<Elec_GSFTrackEta_[NumGsfEle_]<<std::endl;
      std::cout<<std::endl;
    }
    
    NumGsfEle_++;
  }//Loop on GsfElectrons


  ///////////////////////////////////////////////////////////////////////////
  ////////////////////////////////PfTaus variables///////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  NumPFTaus_ = 0;
  for ( reco::PFTauCollection::const_iterator PfTau = PfTaus->begin();
	PfTau != PfTaus->end(); ++PfTau ) {
    if (debug_){
      std::cout<<std::endl;
      std::cout<<"PfTau number : "<<NumPFTaus_<<std::endl;
    }
    if (NumPFTaus_>49) continue; //Accept first 50 taus


    /////////////////////////////////////Matchings  /////////////////////////////////////
    Tau_GsfEleMatch_[NumPFTaus_] = 0;
    Tau_GenEleMatch_[NumPFTaus_] = 0;
    Tau_GenEleFromZMatch_[NumPFTaus_] = 0;
    Tau_GenEleFromZTauTauMatch_[NumPFTaus_] = 0;
    Tau_GenHadMatch_[NumPFTaus_] = 0;
    Tau_GenJetMatch_[NumPFTaus_] = 0;
    for ( reco::GsfElectronCollection::const_iterator GsfElectron = GsfElectrons->begin();
	  GsfElectron != GsfElectrons->end(); ++GsfElectron ) {
      if(debug_){
	std::cout<<"  DeltaR PfTau-GsfEle :"<<deltaR(PfTau->eta(),PfTau->phi(),GsfElectron->eta(),GsfElectron->phi())<<std::endl;
      }
      if(deltaR(PfTau->eta(),PfTau->phi(),GsfElectron->eta(),GsfElectron->phi())<0.3)Tau_GsfEleMatch_[NumPFTaus_] = 1;
    }
    for ( reco::CandidateView::const_iterator GenElectron = GenElectrons->begin();
	  GenElectron != GenElectrons->end(); ++GenElectron ) {
      if(debug_){
	std::cout<<"  DeltaR PfTau-GenEle :"<<deltaR(PfTau->eta(),PfTau->phi(),GenElectron->eta(),GenElectron->phi())<<std::endl;
      }
      if(deltaR(PfTau->eta(),PfTau->phi(),GenElectron->eta(),GenElectron->phi())<0.3)Tau_GenEleMatch_[NumPFTaus_] = 1;
    }
    for ( reco::CandidateView::const_iterator GenElectronFromZ = GenElectronsFromZ->begin();
	  GenElectronFromZ != GenElectronsFromZ->end(); ++GenElectronFromZ ) {
      if(debug_){
	std::cout<<"  DeltaR PfTau-GenEleFromZ :"<<deltaR(PfTau->eta(),PfTau->phi(),GenElectronFromZ->eta(),GenElectronFromZ->phi())<<std::endl;
      }
      if(deltaR(PfTau->eta(),PfTau->phi(),GenElectronFromZ->eta(),GenElectronFromZ->phi())<0.3) Tau_GenEleFromZMatch_[NumPFTaus_] = 1;
    }
    for ( reco::CandidateView::const_iterator GenElectronFromZTauTau = GenElectronsFromZTauTau->begin();
	  GenElectronFromZTauTau != GenElectronsFromZTauTau->end(); ++GenElectronFromZTauTau ) {
      if(debug_){
	std::cout<<"  DeltaR PfTau-GenEleFromZTauTau :"<<deltaR(PfTau->eta(),PfTau->phi(),GenElectronFromZTauTau->eta(),GenElectronFromZTauTau->phi())<<std::endl;
      }
      if(deltaR(PfTau->eta(),PfTau->phi(),GenElectronFromZTauTau->eta(),GenElectronFromZTauTau->phi())<0.3) Tau_GenEleFromZTauTauMatch_[NumPFTaus_] = 1;
    }
    for ( reco::CandidateView::const_iterator GenTau = GenTaus->begin();
	  GenTau != GenTaus->end(); ++GenTau ) {
      if (debug_){
	std::cout<<"  DeltaR PfTau-GenTau :"<<deltaR(PfTau->eta(),PfTau->phi(),GenTau->eta(),GenTau->phi())<<std::endl;
      }
      if(deltaR(PfTau->eta(),PfTau->phi(),GenTau->eta(),GenTau->phi())<0.3)Tau_GenHadMatch_[NumPFTaus_] = 1;
    }
    for ( reco::CandidateView::const_iterator GenJet = GenJets->begin();
	  GenJet != GenJets->end(); ++GenJet ) {
      if (debug_){
	std::cout<<"  DeltaR PfTau-GenJet :"<<deltaR(PfTau->eta(),PfTau->phi(),GenJet->eta(),GenJet->phi())<<std::endl;
      }
      if(deltaR(PfTau->eta(),PfTau->phi(),GenJet->eta(),GenJet->phi())<0.3)Tau_GenJetMatch_[NumPFTaus_] = 1;
    }
    /////////////////////////////////////Matchings  /////////////////////////////////////

    Tau_AbsEta_[NumPFTaus_] = TMath::Abs(PfTau->eta());
    Tau_Pt_[NumPFTaus_] = PfTau->pt();
    Tau_EmFraction_[NumPFTaus_] = TMath::Max(PfTau->emFraction(),float(0.0));
    Tau_NumChargedCands_[NumPFTaus_] = PfTau->signalPFChargedHadrCands().size();
    Tau_NumGammaCands_[NumPFTaus_]  = PfTau->signalPFGammaCands().size();
    Tau_HasGsf_[NumPFTaus_] = (PfTau->leadPFChargedHadrCand()->gsfTrackRef()).isNonnull();
    Tau_HadrHoP_[NumPFTaus_] = PfTau->leadPFChargedHadrCand()->hcalEnergy()/PfTau->leadPFChargedHadrCand()->p();
    Tau_HadrEoP_[NumPFTaus_] = PfTau->leadPFChargedHadrCand()->ecalEnergy()/PfTau->leadPFChargedHadrCand()->p();

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

    Tau_GammaEtaMom_[NumPFTaus_] = TMath::Sqrt(dEta2)*TMath::Sqrt(GammadPt_)*PfTau->pt();
    Tau_GammaPhiMom_[NumPFTaus_] = TMath::Sqrt(dPhi2)*TMath::Sqrt(GammadPt_)*PfTau->pt();  
    Tau_GammaEnFrac_[NumPFTaus_] = GammadPt_;
    Tau_VisMass_[NumPFTaus_] = PfTau->mass();
    Tau_HadrMva_[NumPFTaus_] = TMath::Max(PfTau->electronPreIDOutput(),float(-1.0));

    if(debug_){
      std::cout<<std::endl;
      std::cout<<"TauAbsEta :"<<Tau_AbsEta_[NumPFTaus_]<<std::endl;
      std::cout<<"TauPt :"<<Tau_Pt_[NumPFTaus_]<<std::endl;
      std::cout<<"TauHasGsf :"<<Tau_HasGsf_[NumPFTaus_]<<std::endl;
      std::cout<<"TauEmFraction :"<<Tau_EmFraction_[NumPFTaus_]<<std::endl;
      std::cout<<"TauNumChargedCands :"<<Tau_NumChargedCands_[NumPFTaus_]<<std::endl;
      std::cout<<"TauNumGammaCands :"<<Tau_NumGammaCands_[NumPFTaus_]<<std::endl;      
      std::cout<<"Tau_HadrHoP :"<<Tau_HadrHoP_[NumPFTaus_]<<std::endl;
      std::cout<<"Tau_HadrEoP :"<<Tau_HadrEoP_[NumPFTaus_]<<std::endl;
      std::cout<<"Tau_VisMass_ :"<<Tau_VisMass_[NumPFTaus_]<<std::endl;
      std::cout<<"Tau_GammaEtaMom_ :"<<Tau_GammaEtaMom_[NumPFTaus_]<<std::endl;
      std::cout<<"Tau_GammaPhiMom_ :"<<Tau_GammaPhiMom_[NumPFTaus_]<<std::endl;
      std::cout<<"Tau_GammaEnFrac_ :"<<Tau_GammaEnFrac_[NumPFTaus_]<<std::endl;
      std::cout<<"Tau_HadrMva_ :"<<Tau_HadrMva_[NumPFTaus_]<<std::endl;
    }
    
    NumPFTaus_++;
  }//Loop on PFTaus

  tree_->Fill();

}//analyze()

AntiEMVAAnalyzer::~AntiEMVAAnalyzer()
{

}
void AntiEMVAAnalyzer::endJob()
{
// nothing to be done yet...
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(AntiEMVAAnalyzer);






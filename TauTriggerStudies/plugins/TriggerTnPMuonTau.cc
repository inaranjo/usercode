// My include
#include "Htautau/TriggerStudies/plugins/TriggerTnPMuonTau.h"

// C/C++ include
#include <memory>
#include <iostream>

// CMSSW include
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
//
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
// Pile UP
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
// Vertices
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
// Trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

// TrackingParticles
#include "FWCore/Framework/interface/EventSetupRecordImplementation.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include <vector>

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "CLHEP/Units/GlobalPhysicalConstants.h"

#include "TLorentzVector.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateMode.h"

// Other specific
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

// Muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

// Taus
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

// MET

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

using namespace std;
using namespace reco;
using namespace edm;
using namespace IPTools;
//using namespace math;

// ====================================================================================
TriggerTnPMuonTau::TriggerTnPMuonTau(const edm::ParameterSet& iConfig) :
VerticesTag_(iConfig.getParameter<edm::InputTag> ("VerticesTag")),
//
MuonTag_ (iConfig.getParameter<edm::InputTag> ("MuonTag")),
// Trigger Stuff
HLTTag_(iConfig.getParameter<edm::InputTag> ("HLTTag")),
triggerEventTag_(iConfig.getParameter<edm::InputTag> ("TriggerEventTag")),
//
PileupSrc_ ("addPileupInfo"),
RhoCorrection_("kt6PFJets:rho"),
SigmaRhoCorrection_("kt6PFJets:sigma"),
//
type_ (iConfig.getParameter<std::string>("type")),
aod_ (iConfig.getUntrackedParameter<bool>("AOD",true))
  
// pflow isolation


//
// ====================================================================================
{
	//now do what ever initialization is needed
	
	HLT_Paths_  = iConfig.getParameter<std::vector<std::string > >("HLTPaths");
	HLT_Filters_   = iConfig.getParameter<std::vector<edm::InputTag > >("HLTFilters");
	
	edm::Service<TFileService> fs ;
	mytree_  = fs->make <TTree>("eIDSimpleTree","eIDSimpleTree"); 
	
	// Global
	mytree_->Branch("nEvent",&nEvent,"nEvent/I");
	mytree_->Branch("nRun",&nRun,"nRun/I");
	mytree_->Branch("nLumi",&nLumi,"nLumi/I");
	
	// Pile UP
	mytree_->Branch("PU_N",&_PU_N,"PU_N/I");
	mytree_->Branch("PU_rhoCorr",&_PU_rho,"PU_rhoCorr/D");
	mytree_->Branch("PU_sigmaCorr",&_PU_sigma,"PU_sigmaCorr/D");

	// Vertices
	mytree_->Branch("vtx_N",&_vtx_N,"vtx_N/I");
	mytree_->Branch("vtx_normalizedChi2",&_vtx_normalizedChi2,"vtx_normalizedChi2[vtx_N]/D");
	mytree_->Branch("vtx_ndof",&_vtx_ndof,"vtx_ndof[vtx_N]/D");
	mytree_->Branch("vtx_nTracks",&_vtx_nTracks,"vtx_nTracks[vtx_N]/D");
	mytree_->Branch("vtx_d0",&_vtx_d0,"vtx_d0[vtx_N]/D");
	mytree_->Branch("vtx_x",&_vtx_x,"vtx_x[vtx_N]/D");
	mytree_->Branch("vtx_y",&_vtx_y,"vtx_y[vtx_N]/D");
	mytree_->Branch("vtx_z",&_vtx_z,"vtx_z[vtx_N]/D");
	
	// Trigger
	mytree_->Branch("trig_fired_names",&trig_fired_names,"trig_fired_names[50000]/C");
	//
	mytree_->Branch("trig_HLT_algoStudied", "std::vector<std::string>",&_trig_HLT_algoStudied) ;
	mytree_->Branch("trig_HLT_filterStudied", "std::vector<std::string>",&_trig_HLT_filterStudied) ;
	//
	mytree_->Branch("trig_HLT_NPath",    &_trig_HLT_NPath,     "trig_HLT_NPath/I");
	mytree_->Branch("trig_HLT_pathname", &_trig_HLT_pathname,   "trig_HLT_pathname[trig_HLT_NPath]/I");
	//
	mytree_->Branch("trig_HLT_N",      &_trig_HLT_N,     "trig_HLT_N/I");
	mytree_->Branch("trig_HLT_eta",    &_trig_HLT_eta,   "trig_HLT_eta[trig_HLT_N]/D");
	mytree_->Branch("trig_HLT_phi",    &_trig_HLT_phi,   "trig_HLT_phi[trig_HLT_N]/D");
	mytree_->Branch("trig_HLT_energy", &_trig_HLT_energy,"trig_HLT_energy[trig_HLT_N]/D");
	mytree_->Branch("trig_HLT_pt",     &_trig_HLT_pt,    "trig_HLT_pt[trig_HLT_N]/D");
	mytree_->Branch("trig_HLT_name",   &_trig_HLT_name,   "trig_HLT_name[trig_HLT_N]/I");
	mytree_->Branch("trig_HLT_id",   &_trig_HLT_id,   "trig_HLT_id[trig_HLT_N]/I");
	mytree_->Branch("trig_HLT_vids",   &_trig_HLT_vids,   "trig_HLT_vids[trig_HLT_N]/I");
	
	// Muons
	mytree_->Branch("muons_N",&_muons_N,"muons_N/I");
	m_muons = new TClonesArray ("TLorentzVector");
	mytree_->Branch("muons", "TClonesArray", &m_muons, 256000,0);
	mytree_->Branch("muons_charge",&_muons_charge,"muons_charge[muons_N]/I");
	mytree_->Branch("muons_istracker",&_muons_istracker,"muons_istracker[muons_N]/I");
	mytree_->Branch("muons_isstandalone",&_muons_isstandalone,"muons_isstandalone[muons_N]/I");
	mytree_->Branch("muons_isglobal",&_muons_isglobal,"muons_isglobal[muons_N]/I");
	//
	mytree_->Branch("muons_dxy",&_muons_dxy,"muons_dxy[muons_N]/D");
	mytree_->Branch("muons_dz",&_muons_dz,"muons_dz[muons_N]/D");
	mytree_->Branch("muons_dxyPV",&_muons_dxyPV,"muons_dxyPV[muons_N]/D");
	mytree_->Branch("muons_dzPV",&_muons_dzPV,"muons_dzPV[muons_N]/D");
	mytree_->Branch("muons_normalizedChi2",&_muons_normalizedChi2,"muons_normalizedChi2[muons_N]/D");
	mytree_->Branch("muons_NtrackerHits",&_muons_NtrackerHits,"muons_NtrackerHits[muons_N]/I");
	mytree_->Branch("muons_NpixelHits",&_muons_NpixelHits,"muons_NpixelHits[muons_N]/I");
	mytree_->Branch("muons_NmuonHits",&_muons_NmuonHits,"muons_NmuonHits[muons_N]/I");
	mytree_->Branch("muons_Nmatches",&_muons_Nmatches,"muons_Nmatches[muons_N]/I");
	//
	mytree_->Branch("muons_nTkIsoR03",&_muons_nTkIsoR03,"muons_nTkIsoR03[muons_N]/I");
	mytree_->Branch("muons_nTkIsoR05",&_muons_nTkIsoR05,"muons_nTkIsoR05[muons_N]/I");
	mytree_->Branch("muons_tkIsoR03",&_muons_tkIsoR03,"muons_tkIsoR03[muons_N]/D");
	mytree_->Branch("muons_tkIsoR05",&_muons_tkIsoR05,"muons_tkIsoR05[muons_N]/D");
	mytree_->Branch("muons_emIsoR03",&_muons_emIsoR03,"muons_emIsoR03[muons_N]/D");
	mytree_->Branch("muons_emIsoR05",&_muons_emIsoR05,"muons_emIsoR05[muons_N]/D");
	mytree_->Branch("muons_hadIsoR03",&_muons_hadIsoR03,"muons_hadIsoR03[muons_N]/D");
	mytree_->Branch("muons_hadIsoR05",&_muons_hadIsoR05,"muons_hadIsoR05[muons_N]/D");
	
	//muons TIP/LIP/IP
	mytree_->Branch("muons_Tip",&muons_Tip,"muons_Tip[muons_N]/D");
	mytree_->Branch("muons_Lip",&muons_Lip,"muons_Lip[muons_N]/D");
	mytree_->Branch("muons_STip",&muons_STip,"muons_STip[muons_N]/D");
	mytree_->Branch("muons_SLip",&muons_SLip,"muons_SLip[muons_N]/D");
	mytree_->Branch("muons_TipSignif",&muons_TipSignif,"muons_TipSignif[muons_N]/D");
	mytree_->Branch("muons_LipSignif",&muons_LipSignif,"muons_LipSignif[muons_N]/D");
	mytree_->Branch("muons_Significance3D",&muons_Significance3D,"muons_Significance3D[muons_N]/D");
	mytree_->Branch("muons_Value3D",&muons_Value3D,"muons_Value3D[muons_N]/D");
	mytree_->Branch("muons_Error3D",&muons_Error3D,"muons_Error3D[muons_N]/D");
	//muonID variables
	mytree_->Branch("muons_trkDxy",&_muons_trkDxy,"muons_trkDxy[muons_N]/D");
	mytree_->Branch("muons_trkDxyError",&_muons_trkDxyError,"muons_trkDxyError[muons_N]/D");
	mytree_->Branch("muons_trkDxyB",&_muons_trkDxyB,"muons_trkDxyB[muons_N]/D");
	mytree_->Branch("muons_trkDz",&_muons_trkDz,"muons_trkDz[muons_N]/D");
	mytree_->Branch("muons_trkDzError",&_muons_trkDzError,"muons_trkDzError[muons_N]/D");
	mytree_->Branch("muons_trkDzB",&_muons_trkDzB,"muons_trkDzB[muons_N]/D"); 
	mytree_->Branch("muons_trkChi2PerNdof",&_muons_trkChi2PerNdof,"muons_trkChi2PerNdof[muons_N]/D");
	mytree_->Branch("muons_trkCharge",&_muons_trkCharge,"muons_trkCharge[muons_N]/D");
	mytree_->Branch("muons_trkNHits",&_muons_trkNHits,"muons_trkNHits[muons_N]/D");
	mytree_->Branch("muons_trkNPixHits",&_muons_trkNPixHits,"muons_trkNPixHits[muons_N]/D");
	mytree_->Branch("muons_trkmuArbitration",&_muons_trkmuArbitration,"muons_trkmuArbitration[muons_N]/D");
	mytree_->Branch("muons_trkmu2DCompatibilityLoose",&_muons_trkmu2DCompatibilityLoose,"muons_trkmu2DCompatibilityLoose[muons_N]/D");
	mytree_->Branch("muons_trkmu2DCompatibilityTight",&_muons_trkmu2DCompatibilityTight,"muons_trkmu2DCompatibilityTight[muons_N]/D");
	mytree_->Branch("muons_trkmuOneStationLoose",&_muons_trkmuOneStationLoose,"muons_trkmuOneStationLoose[muons_N]/D");
	mytree_->Branch("muons_trkmuOneStationTight",&_muons_trkmuOneStationTight,"muons_trkmuOneStationTight[muons_N]/D");
	mytree_->Branch("muons_trkmuLastStationLoose",&_muons_trkmuLastStationLoose,"muons_trkmuLastStationLoose[muons_N]/D");
	mytree_->Branch("muons_trkmuLastStationTight",&_muons_trkmuLastStationTight,"muons_trkmuLastStationTight[muons_N]/D");
	mytree_->Branch("muons_trkmuOneStationAngLoose",&_muons_trkmuOneStationAngLoose,"muons_trkmuOneStationAngLoose[muons_N]/D");
	mytree_->Branch("muons_trkmuOneStationAngTight",&_muons_trkmuOneStationAngTight,"muons_trkmuOneStationAngTight[muons_N]/D");
	mytree_->Branch("muons_trkmuLastStationAngLoose",&_muons_trkmuLastStationAngLoose,"muons_trkmuLastStationAngLoose[muons_N]/D");
	mytree_->Branch("muons_trkmuLastStationAngTight",&_muons_trkmuLastStationAngTight,"muons_trkmuLastStationAngTight[muons_N]/D");
	mytree_->Branch("muons_trkmuLastStationOptimizedLowPtLoose",&_muons_trkmuLastStationOptimizedLowPtLoose,"muons_trkmuLastStationOptimizedLowPtLoose[muons_N]/D");
	mytree_->Branch("muons_trkmuLastStationOptimizedLowPtTight",&_muons_trkmuLastStationOptimizedLowPtTight,"muons_trkmuLastStationOptimizedLowPtTight[muons_N]/D");
	mytree_->Branch("muons_caloCompatibility",&_muons_caloCompatibility,"muons_caloCompatibility[muons_N]/D");
	mytree_->Branch("muons_segmentCompatibility",&_muons_segmentCompatibility,"muons_segmentCompatibility[muons_N]/D");
	mytree_->Branch("muons_glbmuPromptTight",&_muons_glbmuPromptTight,"muons_glbmuPromptTight[muons_N]/D");	

	// Taus	
	mytree_->Branch("hpsTau_N", &_nhpsTau, "hpsTau_N/I");
	mytree_->Branch("hpsTau_eta", &_hpsTauEta, "hpsTau_eta[hpsTau_N]/F");
	mytree_->Branch("hpsTau_phi", &_hpsTauPhi, "hpsTau_phi[hpsTau_N]/F");
	mytree_->Branch("hpsTau_pt", &_hpsTauPt, "hpsTau_pt[hpsTau_N]/F");
	mytree_->Branch("hpsTau_vx", &_hpsTauVx, "hpsTau_vx[hpsTau_N]/F");
	mytree_->Branch("hpsTau_vy", &_hpsTauVy, "hpsTau_vy[hpsTau_N]/F");
	mytree_->Branch("hpsTau_vz", &_hpsTauVz, "hpsTau_vz[hpsTau_N]/F");
	mytree_->Branch("hpsTau_jet_pt", &_hpsTauJetPt, "hpsTau_jet_pt[hpsTau_N]/F");
	mytree_->Branch("hpsTau_leadPion_pt", &_hpsTauLeadPionPt, "hpsTau_leadPion_pt[hpsTau_N]/F");
	mytree_->Branch("hpsTau_leadTrack_pt", &_hpsTauLeadTrackPt, "hpsTau_leadTrack_pt[hpsTau_N]/F");
	mytree_->Branch("hpsTau_charge", &_hpsTauCharge, "hpsTau_charge[hpsTau_N]/I");
	mytree_->Branch("hpsTau_chIso", &_hpsTauChargedIso, "hpsTau_chIso[hpsTau_N]/F");
	mytree_->Branch("hpsTau_phIso", &_hpsTauPhotonsIso, "hpsTau_phIso[hpsTau_N]/F");
	mytree_->Branch("hpsTau_decayMode", &_hpsTauDiscrByDecMode, "hpsTau_decayMode[hpsTau_N]/F");
	mytree_->Branch("hpsTau_isoL", &_hpsTauDiscrByLooseIso, "hpsTau_isoL[hpsTau_N]/F");
	mytree_->Branch("hpsTau_isoM", &_hpsTauDiscrByMediumIso, "hpsTau_isoM[hpsTau_N]/F");
	mytree_->Branch("hpsTau_antiMuL", &_hpsTauDiscrAgainstMuonLoose, "hpsTau_antiMuL[hpsTau_N]/F");
	mytree_->Branch("hpsTau_antiMuT", &_hpsTauDiscrAgainstMuonTight, "hpsTau_antiMuT[hpsTau_N]/F");
	mytree_->Branch("hpsTau_antiElL", &_hpsTauDiscrAgainstElecLoose, "hpsTau_antiElL[hpsTau_N]/F");
	mytree_->Branch("hpsTau_anitElM", &_hpsTauDiscrAgainstElecMedium, "hpsTau_antiElM[hpsTau_N]/F");
	mytree_->Branch("hpsTau_antiElT", &_hpsTauDiscrAgainstElecTight, "hpsTau_antiElT[hpsTau_N]/F");
	mytree_->Branch("hpsTau_isoL_DBSumPtCorr", &_hpsTauDiscrByLooseIsoDBSumPtCorr, "hpsTau_isoL_DBSumPtCorr[hpsTau_N]/F");
	mytree_->Branch("hpsTau_isoM_DBSumPtCorr", &_hpsTauDiscrByMediumIsoDBSumPtCorr, "hpsTau_isoM_DBSumPtCorr[hpsTau_N]/F");
	mytree_->Branch("hpsTau_isoT_DBSumPtCorr", &_hpsTauDiscrByTightIsoDBSumPtCorr, "hpsTau_isoT_DBSumPtCorr[hpsTau_N]/F");
	mytree_->Branch("hpsTau_CombinedisoL_DBSumPtCorr", &_hpsTauDiscrByLooseCombinedIsoDBSumPtCorr, "hpsTau_CombinedisoL_DBSumPtCorr[hpsTau_N]/F");
	mytree_->Branch("hpsTau_CombinedisoM_DBSumPtCorr", &_hpsTauDiscrByMediumCombinedIsoDBSumPtCorr, "hpsTau_CombinedisoM_DBSumPtCorr[hpsTau_N]/F");
	mytree_->Branch("hpsTau_CombinedisoT_DBSumPtCorr", &_hpsTauDiscrByTightCombinedIsoDBSumPtCorr, "hpsTau_CombinedisoT_DBSumPtCorr[hpsTau_N]/F");
	// MET

	
	mytree_->Branch("met_pf_et",&_met_pf_et,"met_pf_et/D");
	mytree_->Branch("met_pf_px",&_met_pf_px,"met_pf_px/D");
	mytree_->Branch("met_pf_py",&_met_pf_py,"met_pf_py/D");
	mytree_->Branch("met_pf_phi",&_met_pf_phi,"met_pf_phi/D");
	mytree_->Branch("met_pf_set",&_met_pf_set,"met_pf_set/D");
	mytree_->Branch("met_pf_sig",&_met_pf_sig,"met_pf_sig/D");
	
	_trig_HLT_algoStudied = new std::vector<std::string> ;
	for (int ipath=0;ipath< (int) HLT_Paths_.size();ipath++) 
	  _trig_HLT_algoStudied->push_back(HLT_Paths_[ipath]) ;
	
	_trig_HLT_filterStudied = new std::vector<std::string> ;
	for (int ifilter=0;ifilter< (int) HLT_Filters_.size();ifilter++) 
	  _trig_HLT_filterStudied->push_back(HLT_Filters_[ifilter].label()) ;

}

// ====================================================================================
TriggerTnPMuonTau::~TriggerTnPMuonTau()
// ====================================================================================
{
	delete m_muons;

}

// ====================================================================================
void TriggerTnPMuonTau::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	
  // Tree Maker
  //std::cout << "Init()" << std::endl;
  Init();
  FillEvent (iEvent, iSetup);

 
  //
  if (FilterTrigger (iEvent, iSetup)) {
    m_muons -> Clear();
    if (FilterMuon (iEvent, iSetup)) {
      //cout<<"MUON"<<endl;
      if (FilterTau (iEvent, iSetup)) {
	//cout<<"TAU"<<endl;
	FillMet(iEvent, iSetup);
	mytree_->Fill();
      }
    }
  }
  
} // analyze

// ====================================================================================
void TriggerTnPMuonTau::FillEvent (const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
  nEvent = iEvent.id().event();
  nRun   = iEvent.id().run();
  nLumi  = iEvent.luminosityBlock();
  
  // -----------------
  // Pile-up
  // -----------------
  if(type_ == "MC") {
    Handle<vector<PileupSummaryInfo> > PupInfo;
    iEvent.getByLabel(PileupSrc_, PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    float sum_nvtx = 0, ave_nvtx = 0 ;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      sum_nvtx += float(PVI->getPU_NumInteractions() );
    } // loop on Pile up
    if (PupInfo->size()>0) ave_nvtx = sum_nvtx/PupInfo->size() ;
    _PU_N = (int)(ave_nvtx+0.5) ;
  } // if MC

  // Rho/FastJet Correction
  Handle<double> rhoHandle, sigmaHandle;
  iEvent.getByLabel(RhoCorrection_, rhoHandle);
  iEvent.getByLabel(SigmaRhoCorrection_, sigmaHandle);
  _PU_rho   = *rhoHandle;
  _PU_sigma = *sigmaHandle;
	

  // -----------------
  // Vertices
  // -----------------
  Handle<reco::VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(VerticesTag_,recoPrimaryVertexCollection);
  
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByType(recoBeamSpotHandle);
  const reco::BeamSpot bs = *recoBeamSpotHandle;
  
  int vtx_counter=0;
  _vtx_N = recoPrimaryVertexCollection->size();
  
  // select the primary vertex as the one with higest sum of (pt)^2 of tracks                                                                               
  PrimaryVertexSorter PVSorter;
  std::vector<reco::Vertex> sortedVertices = PVSorter.sortedList( *(recoPrimaryVertexCollection.product()) );
  
  if(_vtx_N > 0) {
    GlobalPoint local_vertexPosition(sortedVertices.front().position().x(),
				     sortedVertices.front().position().y(),
				     sortedVertices.front().position().z());
    vertexPosition = local_vertexPosition;
  }
  else {
    GlobalPoint local_vertexPosition(bs.position().x(),
				     bs.position().y(),
				     bs.position().z());
    vertexPosition = local_vertexPosition;
  }
  for( std::vector<reco::Vertex>::const_iterator PV = sortedVertices.begin(); PV != sortedVertices.end(); ++PV){
    if(vtx_counter > 24 ) continue;
    
    _vtx_normalizedChi2[vtx_counter] = PV->normalizedChi2();
    _vtx_ndof[vtx_counter] = PV->ndof();
    _vtx_nTracks[vtx_counter] = PV->tracksSize();
    _vtx_d0[vtx_counter] = PV->position().Rho();
    _vtx_x[vtx_counter] = PV->x();
    _vtx_y[vtx_counter] = PV->y();
    _vtx_z[vtx_counter] = PV->z();
    
    vtx_counter++;
  } // for loop on primary vertices
  
  if(vtx_counter>24) { _vtx_N = 25; cout << "Number of primary vertices>25, vtx_N set to 25" << endl;}
	
}

// ====================================================================================
bool TriggerTnPMuonTau::FilterTrigger (const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	
  // ----------------------------------------------
  //  Get HLT info
  // ----------------------------------------------
  
  edm::Handle<edm::TriggerResults> triggerResultsHandle;
  iEvent.getByLabel (HLTTag_,triggerResultsHandle);
  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResultsHandle);
  
//   	// Get List of available Triggers in the menu
//   	for (int in=0;in<(int)triggerNames.size();in++) {
//   	  cout << " Trigger Names " << in << " = " << triggerNames.triggerName(in) << endl;
//   	} // for loop in triggernames
  
  
  // LOOP Over Trigger Results
  bool eventOk(false) ;
  strcpy(trig_fired_names,"*") ;
  int hltpathcounter(0) ;
  for (int iHLT=0; iHLT<static_cast<int>(triggerResultsHandle->size()); iHLT++) {	
    
    // triggers fired
    if (triggerResultsHandle->accept (iHLT)) {
      //cout << " Trigger Names " << iHLT << " = " << triggerNames.triggerName(iHLT) << endl;
      if ( strlen(trig_fired_names) <= 49550) {
	const char* c_str() ;
	string hlt_string = triggerNames.triggerName(iHLT) ;
	strcat(trig_fired_names,hlt_string.c_str()) ;
	strcat(trig_fired_names,"*") ;
      } 

      //  considere only event triggered matching user's requirement
      for (int ipath=0;ipath< (int) HLT_Paths_.size();ipath++)
	if (string(triggerNames.triggerName(iHLT)).find(HLT_Paths_[ipath]) != std::string::npos) {  
	  _trig_HLT_pathname[hltpathcounter] = ipath ;
	  hltpathcounter++ ;
	  eventOk = true ;
	}
    }   
  } // end LOOP Over Trigger Results
		
  _trig_HLT_NPath = hltpathcounter ;

  if (eventOk) {
	
    // ----------------------
    //  get HLT candidate
    // ----------------------
    edm::Handle<trigger::TriggerEvent> trigEvent;
    iEvent.getByLabel(triggerEventTag_, trigEvent);
	
    const Int_t N_filter(trigEvent->sizeFilters());
    std::vector<Int_t> ID_filter; 
	
//     // Print Official Filters
//     for(int ifi=0;ifi<N_filter;ifi++) {
//       cout << "filter tag " << ifi << " = " << trigEvent->filterTag(ifi) << endl;
//     } // for loop on filters
    
    int hlt_counter = 0;

    // Loop on user's Filters
    for(int itrig=0;itrig< (int) HLT_Filters_.size();itrig++) {		
		
      ID_filter.push_back(trigEvent->filterIndex(HLT_Filters_[itrig]));       
      const trigger::TriggerObjectCollection& TOC(trigEvent->getObjects());

      if( ID_filter[itrig] <  N_filter) {
	const trigger::Keys & keys (trigEvent->filterKeys(ID_filter[itrig])); 
	const trigger::Vids & VIDS (trigEvent->filterIds(ID_filter[itrig])) ;
	if (keys.size() == VIDS.size()) {
	  // Loop on HLT objects
	  int ikey = 0 ;
	  for (trigger::Keys::const_iterator keyIt = keys.begin();keyIt != keys.end();++keyIt) { 
	    const trigger::TriggerObject& obj = TOC[*keyIt];	  
// 	    cout<<itrig<<" "<<VIDS[ikey]<<"/"<<keys[ikey]<<endl ;
// 	    cout<<abs(obj.id())<<" "<<obj.eta()<<" "<<obj.phi()<<" "<<obj.et()<<" "<<obj.particle().charge()<<endl ;
	    if(hlt_counter>49) continue;
	    ikey++ ;
	    _trig_HLT_eta[hlt_counter]    = obj.eta();
	    _trig_HLT_phi[hlt_counter]    = obj.phi();
	    _trig_HLT_energy[hlt_counter] = obj.energy();
	    _trig_HLT_pt[hlt_counter]     = obj.pt();
	    _trig_HLT_id[hlt_counter]     = obj.id();
	    _trig_HLT_vids[hlt_counter]   = VIDS[ikey];
	    _trig_HLT_name[hlt_counter]   = itrig;
	    hlt_counter++;
	  }
	}
	else cout<<"BUGGGGGGG"<<endl ;

      } // if idfilter<trigevent size
    } // for loop on filters

    _trig_HLT_N = hlt_counter;
    if(hlt_counter>49) { _trig_HLT_N = 50; cout << "Number of HLT Objects>50, trig_HLT_N set to 50" << endl;}


  } // if eventOk


  return eventOk ;

} // end of FillTrigger



// ====================================================================================
bool TriggerTnPMuonTau::FilterMuon(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
  // Beam spot
  //Handle<reco::BeamSpot> beamSpotHandle;
  //iEvent.getByLabel(InputTag("offlineBeamSpot"), beamSpotHandle);
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle ;
  iEvent.getByType(recoBeamSpotHandle) ;
  const reco::BeamSpot bs = *recoBeamSpotHandle ;
  
  // Muon Retrieving
  Handle<View<reco::Muon> > MuonHandle;
  iEvent.getByLabel(MuonTag_, MuonHandle);
  
  TClonesArray &muons = *m_muons;
  int mu_counter = 0;

  
  // ----------------------------------------------
  //  Loop over Muons
  // ----------------------------------------------
  _muons_N = MuonHandle->size();
  
  
  for (edm::View<reco::Muon>::const_iterator imuons=MuonHandle->begin(); imuons!=MuonHandle->end(); ++imuons) {  
    if(mu_counter>19) continue;
    
    // 4-vector
    //edm::Ref<reco::MuonCollection> muonEdmRef(MuonHandle,i);
    setMomentum (myvector, imuons->p4());
    new (muons[mu_counter]) TLorentzVector (myvector);
    
    _muons_charge[mu_counter] = imuons->charge(); 
    
    // Provenance
    if(imuons->isTrackerMuon())    _muons_istracker[mu_counter]    = 1;
    if(imuons->isStandAloneMuon()) _muons_isstandalone[mu_counter] = 1;
    if(imuons->isGlobalMuon())     _muons_isglobal[mu_counter]     = 1;
    
    // Quality cuts
    reco::TrackRef gm = imuons->globalTrack();
    reco::TrackRef tk = imuons->innerTrack();
    
    if(imuons->isGlobalMuon()==1) {
      _muons_dxy[mu_counter]            = gm->dxy(bs.position()); //beamSpotHandle->position());
      _muons_dz[mu_counter]             = gm->dz(bs.position()); //beamSpotHandle->position());
      _muons_dxyPV[mu_counter]          = gm->dxy(math::XYZPoint(vertexPosition)); //beamSpotHandle->position());
      _muons_dzPV[mu_counter]           = gm->dz(math::XYZPoint(vertexPosition)); //beamSpotHandle->position());
      _muons_normalizedChi2[mu_counter] = gm->normalizedChi2(); 
      _muons_NmuonHits[mu_counter]      = gm->hitPattern().numberOfValidMuonHits(); // muon hit matched to global fit
    } // if Global Track
    
    if(imuons->innerTrack().isAvailable()){
      _muons_trkDxy[mu_counter]=imuons->innerTrack()->dxy();
      _muons_trkDxyError[mu_counter]=imuons->innerTrack()->dxyError();
      _muons_trkDxyB[mu_counter]=imuons->innerTrack()->dxy(bs.position()) ;
      _muons_trkDz[mu_counter]=imuons->innerTrack()->dz();
      _muons_trkDzError[mu_counter]=imuons->innerTrack()->dzError();
      _muons_trkDzB[mu_counter]=imuons->innerTrack()->dz(bs.position());
      _muons_trkChi2PerNdof[mu_counter]=imuons->innerTrack()->normalizedChi2();
      _muons_trkCharge[mu_counter]=imuons->innerTrack()->charge();
      _muons_trkNHits[mu_counter]=imuons->innerTrack()->numberOfValidHits();
      _muons_trkNPixHits[mu_counter]=imuons->innerTrack()->hitPattern().numberOfValidPixelHits();
      // Tracker muon properties
      _muons_trkmuArbitration[mu_counter]=(muon::segmentCompatibility( (*imuons),reco::Muon::SegmentAndTrackArbitration));
      _muons_trkmu2DCompatibilityLoose[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TM2DCompatibilityLoose));
      _muons_trkmu2DCompatibilityTight[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TM2DCompatibilityTight));
      _muons_trkmuOneStationLoose[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMOneStationLoose));
      _muons_trkmuOneStationTight[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMOneStationTight));
      _muons_trkmuLastStationLoose[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMLastStationLoose));
      _muons_trkmuLastStationTight[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMLastStationTight));
      _muons_trkmuOneStationAngLoose[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMOneStationAngLoose));
      _muons_trkmuOneStationAngTight[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMOneStationAngTight));
      _muons_trkmuLastStationAngLoose[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMLastStationAngLoose));
      _muons_trkmuLastStationAngTight[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMLastStationAngTight));
      _muons_trkmuLastStationOptimizedLowPtLoose[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMLastStationOptimizedLowPtLoose));
      _muons_trkmuLastStationOptimizedLowPtTight[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMLastStationOptimizedLowPtTight));
    }
    
    if(imuons->isGlobalMuon()==1 || imuons->isTrackerMuon()==1) {
      _muons_NtrackerHits[mu_counter]   = tk->hitPattern().numberOfValidTrackerHits();
      _muons_NpixelHits[mu_counter]     = tk->hitPattern().numberOfValidPixelHits();
    } // if Tracker track
    _muons_Nmatches[mu_counter]             = imuons->numberOfMatches(); // number of segments matched to muon stations
    _muons_caloCompatibility[mu_counter]    = imuons->caloCompatibility() ;
    _muons_segmentCompatibility[mu_counter] = ( muon::segmentCompatibility ( (*imuons) , reco::Muon::SegmentAndTrackArbitration) ) ;
    _muons_glbmuPromptTight[mu_counter]     = ( muon::isGoodMuon( (*imuons) , muon::GlobalMuonPromptTight) );
    
    // Isolation
    _muons_nTkIsoR03[mu_counter] = imuons->isolationR03().nTracks; 
    _muons_nTkIsoR05[mu_counter] = imuons->isolationR05().nTracks;
    _muons_tkIsoR03[mu_counter]  = imuons->isolationR03().sumPt;
    _muons_tkIsoR05[mu_counter]  = imuons->isolationR05().sumPt;
    _muons_emIsoR03[mu_counter]  = imuons->isolationR03().emEt;
    _muons_emIsoR05[mu_counter]  = imuons->isolationR05().emEt;
    _muons_hadIsoR03[mu_counter] = imuons->isolationR03().hadEt;
    _muons_hadIsoR05[mu_counter] = imuons->isolationR05().hadEt;
    
    
    mu_counter++;
  } // for loop on muons
  if(mu_counter>19) { _muons_N = 20; cout << "Number of muons>20, muons_N set to 20" << endl;}
  
  if(mu_counter>0 )return true ; 
  else return false ; 
}

// ====================================================================================
bool TriggerTnPMuonTau::FilterTau(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{

  //Handle tau objects from data
  edm::Handle<reco::PFTauCollection> hpsTaus;
  iEvent.getByLabel("hpsPFTauProducer",hpsTaus);

  edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByDecMode;
  iEvent.getByLabel("hpsPFTauDiscriminationByDecayModeFinding",hpsTauDiscrByDecMode);

  edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByLooseIsolation;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseIsolation",hpsTauDiscrByLooseIsolation);
  edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByMediumIsolation;
  iEvent.getByLabel("hpsPFTauDiscriminationByMediumIsolation",hpsTauDiscrByMediumIsolation);

  edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrAgainstElecLoose;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseElectronRejection",hpsTauDiscrAgainstElecLoose);
  edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrAgainstElecMedium;
  iEvent.getByLabel("hpsPFTauDiscriminationByMediumElectronRejection",hpsTauDiscrAgainstElecMedium);
  edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrAgainstElecTight;
  iEvent.getByLabel("hpsPFTauDiscriminationByTightElectronRejection",hpsTauDiscrAgainstElecTight);

  edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrAgainstMuonLoose;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseMuonRejection",hpsTauDiscrAgainstMuonLoose);
  edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrAgainstMuonTight;
  iEvent.getByLabel("hpsPFTauDiscriminationByTightMuonRejection",hpsTauDiscrAgainstMuonTight);

  edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByLooseIsolationDBSumPtCorr;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr",hpsTauDiscrByLooseIsolationDBSumPtCorr);
  edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByMediumIsolationDBSumPtCorr;
  iEvent.getByLabel("hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr",hpsTauDiscrByMediumIsolationDBSumPtCorr);
  edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByTightIsolationDBSumPtCorr;
  iEvent.getByLabel("hpsPFTauDiscriminationByTightIsolationDBSumPtCorr",hpsTauDiscrByTightIsolationDBSumPtCorr);
  edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByLooseCombinedIsolationDBSumPtCorr;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr",hpsTauDiscrByLooseCombinedIsolationDBSumPtCorr);
  edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByMediumCombinedIsolationDBSumPtCorr;
  iEvent.getByLabel("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr",hpsTauDiscrByMediumCombinedIsolationDBSumPtCorr);
  edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByTightCombinedIsolationDBSumPtCorr;
  iEvent.getByLabel("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr",hpsTauDiscrByTightCombinedIsolationDBSumPtCorr);

  int igoodpftau=0;

  //To got what is need for analysis
  if(hpsTaus.isValid() ) {
    reco::PFTauCollection taus = *hpsTaus;
    typedef reco::PFTauCollection::const_iterator pftauit;
    int ipftau=0;
    
    for(pftauit i=taus.begin(); i!=taus.end(); ++i, ++ipftau){
      if(igoodpftau>4) break; //accept first 5 taus
      if(i->pt()<10) continue;
      if(fabs(i->eta() )>2.3) continue;
      const PFTauRef thisTauRef(hpsTaus,ipftau);
      //check decay mode
      if(!hpsTauDiscrByDecMode.isValid() ) continue;
      else if( (*hpsTauDiscrByDecMode)[thisTauRef] < 0.5 ) continue;
      //check loose isolation (not needed)
      //if(!hpsTauDiscrByLooseIsolation.isValid() ) continue;
      //else if( (*hpsTauDiscrByLooseIsolation)[thisTauRef] < 0.5 ) continue;
      
      //fill tree
      _hpsTauEta[igoodpftau] = i->eta();
      _hpsTauPhi[igoodpftau] = i->phi();
      _hpsTauPt[igoodpftau]  = i->pt();
      _hpsTauVx[igoodpftau]  = i->vertex().x();
      _hpsTauVy[igoodpftau]  = i->vertex().y();
      _hpsTauVz[igoodpftau]  = i->vertex().z();
      _hpsTauJetPt[igoodpftau] = i->pfTauTagInfoRef()->pfjetRef()->pt();
      _hpsTauChargedIso[igoodpftau] = i->isolationPFChargedHadrCandsPtSum();
      _hpsTauPhotonsIso[igoodpftau] = i->isolationPFGammaCandsEtSum();
      _hpsTauCharge[igoodpftau] = i->charge();

      //       cout<<"vertex x  "<<_hpsTauVx[igoodpftau]<<endl;
      //       cout<<"vertex y  "<<_hpsTauVy[igoodpftau]<<endl;
      //       cout<<"vertex z  "<<_hpsTauVz[igoodpftau]<<endl;

      if( (i->leadPFNeutralCand() ).isNonnull() )
	_hpsTauLeadPionPt[igoodpftau] = i->leadPFNeutralCand()->pt();
      else         
	_hpsTauLeadPionPt[igoodpftau] = -1;
      
      if( (i->leadPFChargedHadrCand() ).isNonnull() )
	_hpsTauLeadTrackPt[igoodpftau] = i->leadPFChargedHadrCand()->pt();
      else 
	_hpsTauLeadTrackPt[igoodpftau] = -1;
      
      
      if(hpsTauDiscrByLooseIsolation.isValid() ){ 
	_hpsTauDiscrByLooseIso[igoodpftau] = (*hpsTauDiscrByLooseIsolation)[thisTauRef];}
      else{
	_hpsTauDiscrByLooseIso[igoodpftau] = -1;}
      
      if(hpsTauDiscrByMediumIsolation.isValid() ){ 
	_hpsTauDiscrByMediumIso[igoodpftau] = (*hpsTauDiscrByMediumIsolation)[thisTauRef];}
      else{
	_hpsTauDiscrByMediumIso[igoodpftau] = -1;}
      
      if(hpsTauDiscrAgainstMuonLoose.isValid() ){
	_hpsTauDiscrAgainstMuonLoose[igoodpftau] = (*hpsTauDiscrAgainstMuonLoose)[thisTauRef];}
      else{
	_hpsTauDiscrAgainstMuonLoose[igoodpftau] = -1;}

      if(hpsTauDiscrAgainstMuonTight.isValid() ){
	_hpsTauDiscrAgainstMuonTight[igoodpftau] = (*hpsTauDiscrAgainstMuonTight)[thisTauRef];}
      else{
	_hpsTauDiscrAgainstMuonTight[igoodpftau] = -1;}
      
      if(hpsTauDiscrAgainstElecLoose.isValid() ){
	_hpsTauDiscrAgainstElecLoose[igoodpftau] = (*hpsTauDiscrAgainstElecLoose)[thisTauRef];}
      else{
	_hpsTauDiscrAgainstElecLoose[igoodpftau] = -1;}

      if(hpsTauDiscrAgainstElecMedium.isValid() ){
	_hpsTauDiscrAgainstElecMedium[igoodpftau] = (*hpsTauDiscrAgainstElecMedium)[thisTauRef];}
      else{
	_hpsTauDiscrAgainstElecMedium[igoodpftau] = -1;}

      if(hpsTauDiscrAgainstElecTight.isValid() ){
	_hpsTauDiscrAgainstElecTight[igoodpftau] = (*hpsTauDiscrAgainstElecTight)[thisTauRef];}
      else{
	_hpsTauDiscrAgainstElecTight[igoodpftau] = -1;}
      
      if(hpsTauDiscrByDecMode.isValid() ){
	_hpsTauDiscrByDecMode[igoodpftau] = (*hpsTauDiscrByDecMode)[thisTauRef];}
      else{
	_hpsTauDiscrByDecMode[igoodpftau] = -1;}


      if(hpsTauDiscrByLooseIsolationDBSumPtCorr.isValid() ){ 
	//cout<<"collection"<<endl;
	_hpsTauDiscrByLooseIsoDBSumPtCorr[igoodpftau] = (*hpsTauDiscrByLooseIsolationDBSumPtCorr)[thisTauRef];}
      else{
	_hpsTauDiscrByLooseIsoDBSumPtCorr[igoodpftau] = -1;}

       if(hpsTauDiscrByMediumIsolationDBSumPtCorr.isValid() ){ 
	 _hpsTauDiscrByMediumIsoDBSumPtCorr[igoodpftau] = (*hpsTauDiscrByMediumIsolationDBSumPtCorr)[thisTauRef];}
      else{
	_hpsTauDiscrByMediumIsoDBSumPtCorr[igoodpftau] = -1;}

       if(hpsTauDiscrByTightIsolationDBSumPtCorr.isValid() ){ 
	 _hpsTauDiscrByTightIsoDBSumPtCorr[igoodpftau] = (*hpsTauDiscrByTightIsolationDBSumPtCorr)[thisTauRef];}
      else{
	_hpsTauDiscrByTightIsoDBSumPtCorr[igoodpftau] = -1;}

       if(hpsTauDiscrByLooseCombinedIsolationDBSumPtCorr.isValid() ){ 
	 _hpsTauDiscrByLooseCombinedIsoDBSumPtCorr[igoodpftau] = (*hpsTauDiscrByLooseCombinedIsolationDBSumPtCorr)[thisTauRef];}
      else{
	_hpsTauDiscrByLooseCombinedIsoDBSumPtCorr[igoodpftau] = -1;}

       if(hpsTauDiscrByMediumCombinedIsolationDBSumPtCorr.isValid() ){ 
	 _hpsTauDiscrByMediumCombinedIsoDBSumPtCorr[igoodpftau] = (*hpsTauDiscrByMediumCombinedIsolationDBSumPtCorr)[thisTauRef];}
      else{
	_hpsTauDiscrByMediumCombinedIsoDBSumPtCorr[igoodpftau] = -1;}

       if(hpsTauDiscrByTightCombinedIsolationDBSumPtCorr.isValid() ){ 
	 _hpsTauDiscrByTightCombinedIsoDBSumPtCorr[igoodpftau] = (*hpsTauDiscrByTightCombinedIsolationDBSumPtCorr)[thisTauRef];}
      else{
	_hpsTauDiscrByTightCombinedIsoDBSumPtCorr[igoodpftau] = -1;}
      
      igoodpftau++;
    }       
  } // HPS tau valid 			    
  _nhpsTau = std::min(igoodpftau,5);

  return (igoodpftau > 0 );  
}


// ====================================================================================
void TriggerTnPMuonTau::FillMet(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	
	// MET object built as the (negative) vector sum of all particles (PFCandidates) reconstructed in the event
	edm::Handle< edm::View<reco::PFMET> > pfMEThandle;
	iEvent.getByLabel("pfMet", pfMEThandle);
	
	
	// PFMET
	_met_pf_et  = (pfMEThandle->front() ).et();
	_met_pf_px  = (pfMEThandle->front() ).px();
	_met_pf_py  = (pfMEThandle->front() ).py();
	_met_pf_phi = (pfMEThandle->front() ).phi();
	_met_pf_set = (pfMEThandle->front() ).sumEt();
	_met_pf_sig = (pfMEThandle->front() ).mEtSig();


	//cout<<_met_pf_et<<endl;
	
} // end of Fill MET






// ====================================================================================
void TriggerTnPMuonTau::Init()
// ====================================================================================
{
  nEvent = 0;
  nRun = 0;
  nLumi = 0;

  //Pile-up
  _PU_N     = 0;
  _PU_rho   = 0.;
  _PU_sigma = 0.;
		
  // Vertices
  _vtx_N = 0; 
  for(int iv=0;iv<25;iv++) {
    _vtx_normalizedChi2[iv] = 0.;
    _vtx_ndof[iv] = 0.;
    _vtx_nTracks[iv] = 0.;
    _vtx_d0[iv] = 0.;
    _vtx_x[iv] = 0.;
    _vtx_y[iv] = 0.;
    _vtx_z[iv] = 0.;
  }// for loop on vertices
	
		
  // Trigger
  // HLT
  _trig_HLT_N = 0;
  for(int ihlt=0;ihlt<50;ihlt++) {
    _trig_HLT_eta[ihlt]    = 0.; 
    _trig_HLT_phi[ihlt]    = 0.; 
    _trig_HLT_energy[ihlt] = 0.; 
    _trig_HLT_pt[ihlt]     = 0.;
    _trig_HLT_name[ihlt]   = -1;
    _trig_HLT_id[ihlt]     = -1;
    _trig_HLT_vids[ihlt]   = -1;
  } // for loop on hlt
  _trig_HLT_NPath = 0 ;
  for (int ipath=0 ; ipath<128 ; ipath++) _trig_HLT_pathname[ipath] = -1 ;
	
	
  // Muons
  _muons_N = 0; 
  
  for(int im=0;im<20;im++) {
    _muons_charge[im] = 0;
    // Provenance
    _muons_istracker[im]    = 0;
    _muons_isstandalone[im] = 0;
    _muons_isglobal[im]     = 0;
    // Quality cuts
    _muons_dxy[im]            = 0.;
    _muons_dz[im]             = 0.;
    _muons_dxyPV[im]            = 0.;
    _muons_dzPV[im]             = 0.;
    _muons_normalizedChi2[im] = 0.;
    _muons_NtrackerHits[im]   = 0; 
    _muons_NpixelHits[im]     = 0; 
    _muons_NmuonHits[im]      = 0; 
    _muons_Nmatches[im]       = 0; 
    // Isolation
    _muons_nTkIsoR03[im] = 0; 
    _muons_nTkIsoR05[im] = 0; 
    _muons_tkIsoR03[im]  = 0.;
    _muons_tkIsoR05[im]  = 0.;
    _muons_emIsoR03[im]  = 0.;
    _muons_emIsoR05[im]  = 0.;
    _muons_hadIsoR03[im] = 0.;
    _muons_hadIsoR05[im] = 0.;
    
    _muons_trkDxy[im] = 0.;
    _muons_trkDxyError[im] = 0.;
    _muons_trkDxyB[im] = 0.;
    _muons_trkDz[im] = 0.;
    _muons_trkDzError[im] = 0.;
    _muons_trkDzB[im] = 0.; 
    _muons_trkChi2PerNdof[im] = 0.;
    _muons_trkCharge[im] = 0.;
    _muons_trkNHits[im] = 0.;
    _muons_trkNPixHits[im] = 0.;
    _muons_trkmuArbitration[im] = 0.;
    _muons_trkmu2DCompatibilityLoose[im] = 0.;
    _muons_trkmu2DCompatibilityTight[im] = 0.;
    _muons_trkmuOneStationLoose[im] = 0.;
    _muons_trkmuOneStationTight[im] = 0.;
    _muons_trkmuLastStationLoose[im] = 0.;
    _muons_trkmuLastStationTight[im] = 0.;
    _muons_trkmuOneStationAngLoose[im] = 0.;
    _muons_trkmuOneStationAngTight[im] = 0.;
    _muons_trkmuLastStationAngLoose[im] = 0.;
    _muons_trkmuLastStationAngTight[im] = 0.;
    _muons_trkmuLastStationOptimizedLowPtLoose[im] = 0.;
    _muons_trkmuLastStationOptimizedLowPtTight[im] = 0.;
    _muons_caloCompatibility[im] = 0.;
    _muons_segmentCompatibility[im] = 0.;
    _muons_glbmuPromptTight[im] = 0.;
    
    muons_Tip[im] = -999. ;
    muons_Lip[im] = -999. ;
    muons_STip[im] = -999. ;
    muons_SLip[im] = -999. ;
    muons_TipSignif[im] = -999. ;
    muons_LipSignif[im] = -999. ;
    muons_Significance3D[im] = -999. ;
    muons_Value3D[im] = -999. ;
    muons_Error3D[im] = -999. ;
  } // for loop on muons
  
  // Offline Taus
  _nhpsTau = 0;
  for(int itau=0; itau<5; ++itau){
    _hpsTauEta[itau] = -99;
    _hpsTauPhi[itau] = -99; 
    _hpsTauPt[itau] = -1; 
    _hpsTauJetPt[itau] = -1;
    _hpsTauLeadPionPt[itau] = -1;
    _hpsTauLeadTrackPt[itau] = -1;
    _hpsTauCharge[itau] = -99;
    _hpsTauChargedIso[itau] = -1;
    _hpsTauPhotonsIso[itau] = -1;
    _hpsTauDiscrByDecMode[itau] = -1;
    _hpsTauDiscrByLooseIso[itau] = -1;
    _hpsTauDiscrByMediumIso[itau] = -1;    
    _hpsTauDiscrAgainstMuonLoose[itau] = -1;
    _hpsTauDiscrAgainstMuonTight[itau] = -1;
    _hpsTauDiscrAgainstElecLoose[itau] = -1;
    _hpsTauDiscrAgainstElecMedium[itau] = -1;
    _hpsTauDiscrAgainstElecTight[itau] = -1;
    _hpsTauDiscrByLooseIsoDBSumPtCorr[itau] = -1;
    _hpsTauDiscrByMediumIsoDBSumPtCorr[itau] = -1;
    _hpsTauDiscrByTightIsoDBSumPtCorr[itau] = -1;
    _hpsTauDiscrByLooseCombinedIsoDBSumPtCorr[itau] = -1;
    _hpsTauDiscrByMediumCombinedIsoDBSumPtCorr[itau] = -1;
    _hpsTauDiscrByTightCombinedIsoDBSumPtCorr[itau] = -1;
  }


	
  // MET
  
  
  _met_pf_et  = 0.;
  _met_pf_px  = 0.; 
  _met_pf_py  = 0.; 
  _met_pf_phi = 0.; 
  _met_pf_set = 0.; 
  _met_pf_sig = 0.; 
  
}

// ====================================================================================
void TriggerTnPMuonTau::beginJob(const edm::ParameterSet& conf)
// ====================================================================================
{
	
}

// ====================================================================================
void TriggerTnPMuonTau::endJob() {}
// ====================================================================================

// ====================================================================================
void TriggerTnPMuonTau::setMomentum (TLorentzVector &myvector, const LorentzVector & mom)
// ====================================================================================
{
  myvector.SetPx (mom.Px());
  myvector.SetPy (mom.Py());
  myvector.SetPz (mom.Pz());
  myvector.SetE (mom.E());
}

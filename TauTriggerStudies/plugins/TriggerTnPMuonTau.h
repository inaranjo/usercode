// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "AnalysisDataFormats/Egamma/interface/ElectronID.h"
#include "AnalysisDataFormats/Egamma/interface/ElectronIDAssociation.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TVector3.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"

#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
//#include "CommonTools/RecoAlgos/src/SuperClusterToCandidate.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"

/// among the includes ///
class MultiTrajectoryStateMode ;

#include "FWCore/Framework/interface/ESHandle.h"

#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionBaseClass.h"
// for H/E
class EgammaTowerIsolation ;
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

//
// class declaration
//

class TriggerTnPMuonTau : public edm::EDAnalyzer {
public:
	explicit TriggerTnPMuonTau(const edm::ParameterSet&);
	~TriggerTnPMuonTau();
	
	typedef math::XYZTLorentzVector LorentzVector ;
	typedef edm::View<reco::Track> trackCollection ;
	
private:
	virtual void beginJob(const edm::ParameterSet& conf) ;
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	
	void Init();
	
	void FillEvent (const edm::Event&, const edm::EventSetup&);
	bool FilterTrigger (const edm::Event&, const edm::EventSetup&);
	bool FilterMuon (const edm::Event&, const edm::EventSetup&);
	bool FilterTau (const edm::Event&, const edm::EventSetup&);
	void FillMet (const edm::Event&, const edm::EventSetup&);
		
	void setMomentum (TLorentzVector &myvector, const LorentzVector & mom) ;

	// ----------member data ---------------------------
	TTree *mytree_;
	
	int nEvent, nRun, nLumi;
	
	// Vertices
	int _vtx_N;
	double _vtx_x[25], _vtx_y[25], _vtx_z[25];
	double _vtx_normalizedChi2[25], _vtx_ndof[25], _vtx_nTracks[25], _vtx_d0[25];
	GlobalPoint vertexPosition;

	//Pile-up
	int _PU_N;
	double _PU_rho, _PU_sigma;  //corrections from FastJets

	// Trigger Paths
	char trig_fired_names[50000];
	//
	std::vector<std::string> * _trig_HLT_algoStudied ;
	std::vector<std::string> * _trig_HLT_filterStudied ;
	// HLT
	int _trig_HLT_N;
	double _trig_HLT_eta[50], _trig_HLT_phi[50], _trig_HLT_energy[50], _trig_HLT_pt[50];
        int _trig_HLT_name[50], _trig_HLT_id[50], _trig_HLT_vids[50];
	int _trig_HLT_NPath ;
	int _trig_HLT_pathname[128] ;

	// Muons
	int _muons_N;
	//TClonesArray * m_muons;
	int _muons_charge[20];
	int _muons_istracker[20], _muons_isstandalone[20], _muons_isglobal[20];
	double _muons_dxy[20], _muons_dz[20], _muons_dxyPV[20], _muons_dzPV[20], _muons_normalizedChi2[20];
	int  _muons_NtrackerHits[20], _muons_NpixelHits[20], _muons_NmuonHits[20], _muons_Nmatches[20];
	int _muons_nTkIsoR03[20], _muons_nTkIsoR05[20];
	double _muons_tkIsoR03[20],_muons_tkIsoR05[20],_muons_emIsoR03[20],_muons_emIsoR05[20],_muons_hadIsoR03[20],_muons_hadIsoR05[20];
	
	double _muons_trkDxy[20], _muons_trkDxyError[20], _muons_trkDxyB[20],
	_muons_trkDz[20], _muons_trkDzError[20], _muons_trkDzB[20], _muons_trkChi2PerNdof[20], 
	_muons_trkCharge[20],_muons_trkNHits[20],_muons_trkNPixHits[20];
	// Tracker muon properties
	double _muons_trkmuArbitration[20],
	_muons_trkmu2DCompatibilityLoose[20],
	_muons_trkmu2DCompatibilityTight[20],
	_muons_trkmuOneStationLoose[20],
	_muons_trkmuOneStationTight[20],
	_muons_trkmuLastStationLoose[20],
	_muons_trkmuLastStationTight[20],
	_muons_trkmuOneStationAngLoose[20],
	_muons_trkmuOneStationAngTight[20],
	_muons_trkmuLastStationAngLoose[20],
	_muons_trkmuLastStationAngTight[20],
	_muons_trkmuLastStationOptimizedLowPtLoose[20],
	_muons_trkmuLastStationOptimizedLowPtTight[20];
	
	double _muons_caloCompatibility[20], _muons_segmentCompatibility[20], _muons_glbmuPromptTight[20] ; 
	

	// Vector for muons
	TClonesArray * m_muons ;	
	TLorentzVector myvector ;  


	//add TIP/LIP/IP variables
	double muons_Tip[20],muons_Lip[20],muons_STip[20],muons_SLip[20],muons_TipSignif[20],muons_LipSignif[20],muons_Significance3D[20],muons_Value3D[20],muons_Error3D[20] ;
	double ele_Tip[10],ele_Lip[10],ele_STip[10],ele_SLip[10],ele_TipSignif[10],ele_LipSignif[10],ele_Significance3D[10],ele_Value3D[10],ele_Error3D[10];
	

	// Taus
	int _nhpsTau;
	float _hpsTauEta[5], _hpsTauPhi[5], _hpsTauPt[5], _hpsTauJetPt[5],_hpsTauVx[5],_hpsTauVy[5],_hpsTauVz[5];
	float _hpsTauLeadPionPt[5], _hpsTauLeadTrackPt[5];
	int _hpsTauCharge[5];
	float _hpsTauChargedIso[5], _hpsTauPhotonsIso[5];
	float _hpsTauDiscrByDecMode[5];
	float _hpsTauDiscrByLooseIso[5], _hpsTauDiscrByMediumIso[5];
	float _hpsTauDiscrAgainstMuonLoose[5], _hpsTauDiscrAgainstMuonTight[5];
	float _hpsTauDiscrAgainstElecLoose[5], _hpsTauDiscrAgainstElecMedium[5], _hpsTauDiscrAgainstElecTight[5];	
	float _hpsTauDiscrByLooseIsoDBSumPtCorr[5],_hpsTauDiscrByMediumIsoDBSumPtCorr[5],_hpsTauDiscrByTightIsoDBSumPtCorr[5], _hpsTauDiscrByLooseCombinedIsoDBSumPtCorr[5],_hpsTauDiscrByMediumCombinedIsoDBSumPtCorr[5],_hpsTauDiscrByTightCombinedIsoDBSumPtCorr[5];

// MET

	double _met_pf_et,_met_pf_px, _met_pf_py, _met_pf_phi, _met_pf_set, _met_pf_sig;

	edm::InputTag VerticesTag_;

	edm::InputTag MuonTag_;


	// Trigger Stuff
	edm::InputTag HLTTag_; 
	edm::InputTag triggerEventTag_;
	std::vector<std::string > HLT_Paths_;
	std::vector<edm::InputTag > HLT_Filters_;

	//Pile-up
	edm::InputTag PileupSrc_;
	edm::InputTag RhoCorrection_, SigmaRhoCorrection_, BetaCorrection_;



	// Pflow isolation

	std::string type_;	
	bool aod_;	

	EcalClusterFunctionBaseClass* funcbase_;
};

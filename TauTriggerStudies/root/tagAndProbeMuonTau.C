#include "readJSONFile.cc"
#include "LumiReWeighting.cc"
#include "LumiReWeighting_Pascal.cc"

#include <string>
#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TRandom.h"


using namespace std;


double deltaR(double eta1, double phi1, double eta2, double phi2)
{
  double deta = eta1 - eta2 ;
  double dphi = phi1 - phi2 ;
  while (dphi > TMath::Pi() ) 
    dphi -= 2*TMath::Pi();
  while (dphi <= -TMath::Pi() ) 
    dphi += 2*TMath::Pi();
  return sqrt(deta*deta + dphi*dphi); 
}

void tagAndProbeMuonTau(string dataset     = "2011BPromptRecov1", 
						string isolation   = "Loose", 
						bool NewTauID      = true,  
						bool MassPlot      = true,  
						bool MuTau         = true,  
						bool DoMC          = false,
						int ptTresh        = 20,
						int MCtype         = 0
						)
{
  cout<<"Studying dataset : "<<dataset<<endl;

  string MassString = "";
  if (MassPlot){
    cout<<"Selection for Mass plot"<<endl;
    MassString = "Mass";
  }
  else{
    cout<<"Selection for Turn on plot"<<endl;
  }
  string Channel;
  if (MuTau){
    cout<<"Selection for MuTau channel"<<endl;
    Channel = "MuTauAna";
  }
  else{
    cout<<"Selection for ElecTau channel"<<endl;
    Channel = "ElecTauAna";
  }   
  
  bool DY(false) ;
  bool WLnu(false) ;
  bool QCD(false) ;
  bool TTJets(false) ; 
  string MCString = "" ; 
  
  cout<<"Considering offline Htautau selection"<<endl ;
  if (DoMC){ 
    cout<<"with MC data"<<endl ; 
  	if (MCtype==0) {
  		DY = true ;
 		cout<<"DY sample"<<endl ;
 		MCString = "DY";
  	}
  	if (MCtype==1) {
  		WLnu = true ;
  	  	cout<<"WLnu sample"<<endl ;
 		MCString = "WLnu";
  	}
	if (MCtype==2) {
  		QCD = true ;
  	  	cout<<"QCD sample"<<endl ;
  		MCString = "QCD";
 	}
	if (MCtype==3) {
  		TTJets = true ;
  	  	cout<<"TTJets sample"<<endl ;
 		MCString = "TTJets";
	}
  }
  
  else cout<<"with real data"<<endl ;

  //MC
  int idtriggersource = 0 ;//HLT_IsoMu17
  int idtagfilter = 0 ;//hltSingleMuIsoL3IsoFiltered17
  int idprobefilter = 3 ;//3 hltFilterIsoMu12IsoPFTau10LooseIsolation

  if(dataset == "2011AMay10ReReco"){
  idtriggersource = 0 ;//5
  idtagfilter= 1 ;//3
  idprobefilter= 7 ;//6for PromptReco 3forMC 11 for Tight tau
  if (isolation == "Loose") idprobefilter= 7 ;//6for PromptReco 3forMC 11 for Tight tau
  if (isolation == "Medium") idprobefilter= -1  ;//6for PromptReco 3forMC 11 for Tight tau
  if (isolation == "Tight") idprobefilter= -1 ;//6for PromptReco 3forMC 11 for Tight tau
  }

  if(dataset == "2011APromptRecov4" && ptTresh==20){
  idtriggersource = 0 ;//5
  idtagfilter= 1 ;//3
  idprobefilter= 8 ;//6for PromptReco 3forMC 11 for Tight tau
  if (isolation == "Loose") idprobefilter= 10 ;
  if (isolation == "Medium") idprobefilter= 12 ;
  if (isolation == "Tight") idprobefilter= 13 ;
  }

  if(dataset == "2011APromptRecov4" && ptTresh==15){
  idtriggersource = 0 ;//5
  idtagfilter= 1 ;//3
  idprobefilter= 7 ;//6for PromptReco 3forMC 11 for Tight tau
  if (isolation == "Loose") idprobefilter= 10 ;//6for PromptReco 3forMC 11 for Tight tau
  }
  if(dataset == "2011AAug05ReReco" && ptTresh==20){
  idtriggersource = 5 ;//5
  idtagfilter= 3 ;//3
  idprobefilter= 8 ;//6for PromptReco 3forMC 11 for Tight tau
  if (isolation == "Loose") idprobefilter= 11 ;//6for PromptReco 3forMC 11 for Tight tau
  if (isolation == "Medium") idprobefilter= 12 ;//6for PromptReco 3forMC 11 for Tight tau
  if (isolation == "Tight") idprobefilter= 13 ;//6for PromptReco 3forMC 11 for Tight tau
  }

  if(dataset == "2011AAug05ReReco" && ptTresh==15){
  idtriggersource = 5 ;//5
  idtagfilter= 3 ;//3
  idprobefilter= 7 ;//6for PromptReco 3forMC 11 for Tight tau
  if (isolation == "Loose") idprobefilter= 10 ;//6for PromptReco 3forMC 11 for Tight tau
  }

  if(dataset == "2011APromptRecov6" && ptTresh==20){
  idtriggersource = 5 ;//5
  idtagfilter= 3 ;//3
  idprobefilter= 8 ;//6for PromptReco 3forMC 11 for Tight tau
  if (isolation == "Loose") idprobefilter= 8 ;//6for PromptReco 3forMC 11 for Tight tau
  if (isolation == "Medium") idprobefilter= 9 ;//6for PromptReco 3forMC 11 for Tight tau
  if (isolation == "Tight") idprobefilter= 10 ;//6for PromptReco 3forMC 11 for Tight tau
  }

  if(dataset == "2011APromptRecov6" && ptTresh==15){
  idtriggersource = 5 ;//5
  idtagfilter= 3 ;//3
  idprobefilter= 7 ;//6for PromptReco 3forMC 11 for Tight tau
  if (isolation == "Loose") idprobefilter= 7 ;//6for PromptReco 3forMC 11 for Tight tau
  }

  if(dataset == "2011BPromptRecov1" && ptTresh==20){
  idtriggersource = 4 ;//HLT_IsoMu24_eta2p1
  idtagfilter= 3 ;//hltSingleMuL2QualIsoL3IsoFiltered24
  idprobefilter= 8 ;
  if (isolation == "Loose") idprobefilter= 11 ;//hltPFTau20TrackLooseIso
  if (isolation == "Medium") idprobefilter= 12 ;//hltPFTauMediumIso20TrackMediumIso
  if (isolation == "Tight") idprobefilter= 13 ;//hltPFTauTightIso20TrackTightIso
  }

  TChain * myChain = new TChain("produceNtuple/eIDSimpleTree") ;
  //std::string dir="/data_CMS/cms/paganini/Single2011Apr22ReReco" ;
  //std::string dir="/data_CMS/cms/paganini/Single2011May10ReReco" ;
  //std::string dir="/data_CMS/cms/paganini/Single2011May10ReReco3" ;
  //std::string dir="/data_CMS/cms/mbluj/Production/test/TriggerTrees/Jun7_v1" ;
  std::string dir;
  if (DoMC) {
    if(DY)dir = "/data_CMS/cms/ivo/TriggerStudies/trees/SingleMuonSummer11DY";
    if(WLnu)dir = "/data_CMS/cms/ivo/TriggerStudies/trees/SingleMuonSummer11WLnu";
    if(QCD)dir = "/data_CMS/cms/ivo/TriggerStudies/trees/SingleMuonSummer11QCD";
    if(TTJets)dir = "/data_CMS/cms/ivo/TriggerStudies/trees/SingleMuonSummer11TTJets";
  }
  else{
    if(dataset == "2011AMay10ReReco")dir = "/data_CMS/cms/ivo/TriggerStudies/trees/Trees_2011AMay10ReReco_iter1";
    if(dataset == "2011APromptRecov4")dir = "/data_CMS/cms/ivo/TriggerStudies/trees/Trees_2011APromptRecov4" ;
    if(dataset == "2011AAug05ReReco")dir = "/data_CMS/cms/ivo/TriggerStudies/trees/Trees_2011AAug05ReReco";
    if(dataset == "2011APromptRecov6")dir = "/data_CMS/cms/ivo/TriggerStudies/trees/Trees_2011APromptRecov6_iter2" ;
    if(dataset == "2011BPromptRecov1")dir = "/data_CMS/cms/ivo/TriggerStudies/trees/Trees_2011BPromptRecov1_iter4" ;
  }
    //dir = "/data_CMS/cms/ivo/TriggerStudies/trees/Trees_PromptReco_2" ;
    //dir = "/data_CMS/cms/ivo/TriggerStudies/trees/Trees_PromptReco" ;
    //dir = "/home/llr/cms/ivo/CMSSW_4_2_4/src/Htautau/TriggerStudies/root/Trees_PromptReco_v6_2/";
    //dir = "/home/llr/cms/ivo/CMSSW_4_2_4/src/Htautau/TriggerStudies/root/Trees_2011B_PromptReco_v1/";
    //dir = "/data_CMS/cms/ivo/TriggerStudies/trees/Trees_2011B_PromptReco_v1_iter1" ;
    //dir = "/data_CMS/cms/ivo/TriggerStudies/trees/Trees_2011BPromptRecov1_iter3" ;
    //dir = "/data_CMS/cms/ivo/TriggerStudies/trees/Trees_PromptReco_v6_iter1" ;
    //dir = "/data_CMS/cms/ivo/TriggerStudies/trees/Trees_2011APromptRecov6_iter2" ;
    //dir = "/data_CMS/cms/ivo/TriggerStudies/trees/Trees_PromptReco_v6" ;
    //dir = "/data_CMS/cms/ivo/TriggerStudies/trees/SingleMuon2011May10ReReco_5" ;
    //dir = "/data_CMS/cms/ivo/TriggerStudies/trees/Trees_2011AMay10ReReco_iter1" ;
    //dir = "/home/llr/cms/ivo/CMSSW_4_2_4/src/Htautau/TriggerStudies/test" ;
//     if(DY)dir = "/data_CMS/cms/ivo/TriggerStudies/trees/Trees_DY_2" ;
//     if(WLnu)dir = "/data_CMS/cms/ivo/TriggerStudies/trees/Trees_WLnu_2" ;
//     if(QCD)dir = "/data_CMS/cms/ivo/TriggerStudies/trees/Trees_QCD_2" ;
//     if(TTJets)dir = "/data_CMS/cms/ivo/TriggerStudies/trees/Trees_TTJets_2" ;


  dir += "/*.root";
  myChain->Add (dir.c_str()) ;

  // Declaration of leaf types
  Int_t           nEvent;
  Int_t           nRun;
  Int_t           nLumi;
  Int_t           PU_N;
  Double_t        PU_rhoCorr;
  Double_t        PU_sigmaCorr;
  Int_t           vtx_N;
  Double_t        vtx_normalizedChi2[25];   //[vtx_N]
  Double_t        vtx_ndof[25];   //[vtx_N]
  Double_t        vtx_nTracks[25];   //[vtx_N]
  Double_t        vtx_d0[25];   //[vtx_N]
  Double_t        vtx_x[25];   //[vtx_N]
  Double_t        vtx_y[25];   //[vtx_N]
  Double_t        vtx_z[25];   //[vtx_N]
  Char_t          trig_fired_names[50000];

  std::vector<std::string> * trig_HLT_algoStudied = 0 ;
  std::vector<std::string> * trig_HLT_filterStudied = 0 ;

  Int_t           trig_HLT_N;
  Double_t        trig_HLT_eta[50];   //[trig_HLT_N]
  Double_t        trig_HLT_phi[50];   //[trig_HLT_N]
  Double_t        trig_HLT_energy[50];   //[trig_HLT_N]
  Double_t        trig_HLT_pt[50];   //[trig_HLT_N]
  Int_t           trig_HLT_name[50];   //[trig_HLT_N]
  Int_t           trig_HLT_id[50];   //[trig_HLT_N]
  Int_t           trig_HLT_vids[50];   //[trig_HLT_N]
  int             trig_HLT_NPath ;
  int             trig_HLT_pathname[128] ;   //[trig_HLT_NPath]


  // Muons
  Int_t _muons_N;
  //TClonesArray * m_muons;
  Int_t _muons_charge[20];
  Int_t _muons_istracker[20];
  Int_t _muons_isstandalone[20];
  Int_t _muons_isglobal[20];
  Double_t _muons_dxy[20];
  Double_t _muons_dz[20];
  Double_t _muons_dxyPV[20];
  Double_t _muons_dzPV[20];
  Double_t _muons_normalizedChi2[20];
  Int_t  _muons_NtrackerHits[20];
  Int_t _muons_NpixelHits[20];
  Int_t _muons_NmuonHits[20];
  Int_t _muons_Nmatches[20];
  Int_t _muons_nTkIsoR03[20];
  Int_t _muons_nTkIsoR05[20];
  Double_t _muons_tkIsoR03[20];
  Double_t _muons_tkIsoR05[20];
  Double_t _muons_emIsoR03[20];
  Double_t _muons_emIsoR05[20];
  Double_t _muons_hadIsoR03[20];
  Double_t _muons_hadIsoR05[20];
  
  Double_t _muons_trkDxy[20];
  Double_t _muons_trkDxyError[20];
  Double_t _muons_trkDxyB[20];
  Double_t _muons_trkDz[20];
  Double_t _muons_trkDzError[20];
  Double_t _muons_trkDzB[20];
  Double_t _muons_trkChi2PerNdof[20]; 
  Double_t _muons_trkCharge[20];
  Double_t _muons_trkNHits[20];
  Double_t _muons_trkNPixHits[20];
  // Tracker muon properties
  Double_t _muons_trkmuArbitration[20];
  Double_t  _muons_trkmu2DCompatibilityLoose[20];
  Double_t  _muons_trkmu2DCompatibilityTight[20];
  Double_t  _muons_trkmuOneStationLoose[20];
  Double_t  _muons_trkmuOneStationTight[20];
  Double_t  _muons_trkmuLastStationLoose[20];
  Double_t  _muons_trkmuLastStationTight[20];
  Double_t  _muons_trkmuOneStationAngLoose[20];
  Double_t  _muons_trkmuOneStationAngTight[20];
  Double_t  _muons_trkmuLastStationAngLoose[20];
  Double_t  _muons_trkmuLastStationAngTight[20];
  Double_t  _muons_trkmuLastStationOptimizedLowPtLoose[20];
  Double_t  _muons_trkmuLastStationOptimizedLowPtTight[20];
  
  Double_t _muons_caloCompatibility[20];
  Double_t _muons_segmentCompatibility[20];
  Double_t _muons_glbmuPromptTight[20] ; 
  
  
  // Vector for muons
  TClonesArray * m_muons ;	
  TLorentzVector myvector ;  
  
  
  //add TIP/LIP/IP variables
  Double_t muons_Tip[20];
  Double_t muons_Lip[20];
  Double_t muons_STip[20];
  Double_t muons_SLip[20];
  Double_t muons_TipSignif[20];
  Double_t muons_LipSignif[20];
  Double_t muons_Significance3D[20];
  Double_t muons_Value3D[20];
  Double_t muons_Error3D[20] ;

  
  
  // Taus
  Int_t _nhpsTau;
  float _hpsTauEta[5];
  float _hpsTauPhi[5];
  float _hpsTauPt[5];

  float _hpsTauVx[5] ;
  float _hpsTauVy[5] ;
  float _hpsTauVz[5] ;

  float _hpsTauJetPt[5];
  float _hpsTauLeadPionPt[5];
  float _hpsTauLeadTrackPt[5];
  Int_t _hpsTauCharge[5];
  float _hpsTauChargedIso[5];
  float _hpsTauPhotonsIso[5];
  float _hpsTauDiscrByDecMode[5];
  float _hpsTauDiscrByLooseIso[5];
  float _hpsTauDiscrByMediumIso[5];
  float _hpsTauDiscrAgainstMuonLoose[5];
  float _hpsTauDiscrAgainstMuonTight[5];
  float _hpsTauDiscrAgainstElecLoose[5];
  float _hpsTauDiscrAgainstElecMedium[5];
  float _hpsTauDiscrAgainstElecTight[5];

  float _hpsTauDiscrByLooseIsoDBSumPtCorr[5];
  float _hpsTauDiscrByMediumIsoDBSumPtCorr[5];
  float _hpsTauDiscrByTightIsoDBSumPtCorr[5];
  float _hpsTauDiscrByLooseCombinedIsoDBSumPtCorr[5];
  float _hpsTauDiscrByMediumCombinedIsoDBSumPtCorr[5];
  float _hpsTauDiscrByTightCombinedIsoDBSumPtCorr[5];
	
// MET

  Double_t _met_pf_et;
  Double_t _met_pf_px;
  Double_t _met_pf_py;
  Double_t _met_pf_phi;
  Double_t _met_pf_set;
  Double_t _met_pf_sig;

 
  myChain->SetBranchAddress("nEvent", &nEvent);
  myChain->SetBranchAddress("nRun", &nRun);
  myChain->SetBranchAddress("nLumi", &nLumi);
  myChain->SetBranchAddress("PU_N", &PU_N);
  myChain->SetBranchAddress("PU_rhoCorr", &PU_rhoCorr);
  myChain->SetBranchAddress("PU_sigmaCorr", &PU_sigmaCorr);
  myChain->SetBranchAddress("vtx_N", &vtx_N);
  myChain->SetBranchAddress("vtx_normalizedChi2", vtx_normalizedChi2);
  myChain->SetBranchAddress("vtx_ndof", vtx_ndof);
  myChain->SetBranchAddress("vtx_nTracks", vtx_nTracks);
  myChain->SetBranchAddress("vtx_d0", vtx_d0);
  myChain->SetBranchAddress("vtx_x", vtx_x);
  myChain->SetBranchAddress("vtx_y", vtx_y);
  myChain->SetBranchAddress("vtx_z", vtx_z);

  myChain->SetBranchAddress("trig_fired_names", trig_fired_names);

  myChain->SetBranchAddress("trig_HLT_algoStudied", &trig_HLT_algoStudied);
  myChain->SetBranchAddress("trig_HLT_filterStudied", &trig_HLT_filterStudied);


  myChain->SetBranchAddress("trig_HLT_N", &trig_HLT_N);
  myChain->SetBranchAddress("trig_HLT_eta", trig_HLT_eta);
  myChain->SetBranchAddress("trig_HLT_phi", trig_HLT_phi);
  myChain->SetBranchAddress("trig_HLT_energy", trig_HLT_energy);
  myChain->SetBranchAddress("trig_HLT_pt", trig_HLT_pt);
  myChain->SetBranchAddress("trig_HLT_name", trig_HLT_name);
  myChain->SetBranchAddress("trig_HLT_id", trig_HLT_id);
  myChain->SetBranchAddress("trig_HLT_vids", trig_HLT_vids);

  myChain->SetBranchAddress("trig_HLT_NPath", &trig_HLT_NPath);
  myChain->SetBranchAddress("trig_HLT_pathname", trig_HLT_pathname);

 
  // Muons
  myChain->SetBranchAddress("muons_N",&_muons_N);
  myChain->SetBranchAddress("muons_charge",&_muons_charge);
  myChain->SetBranchAddress("muons_istracker",&_muons_istracker);
  myChain->SetBranchAddress("muons_isstandalone",&_muons_isstandalone);
  myChain->SetBranchAddress("muons_isglobal",&_muons_isglobal);
  //
  myChain->SetBranchAddress("muons_dxy",&_muons_dxy);
  myChain->SetBranchAddress("muons_dz",&_muons_dz);
  myChain->SetBranchAddress("muons_dxyPV",&_muons_dxyPV);
  myChain->SetBranchAddress("muons_dzPV",&_muons_dzPV);
  myChain->SetBranchAddress("muons_normalizedChi2",&_muons_normalizedChi2);
  myChain->SetBranchAddress("muons_NtrackerHits",&_muons_NtrackerHits);
  myChain->SetBranchAddress("muons_NpixelHits",&_muons_NpixelHits);
  myChain->SetBranchAddress("muons_NmuonHits",&_muons_NmuonHits);
  myChain->SetBranchAddress("muons_Nmatches",&_muons_Nmatches);
  //
  myChain->SetBranchAddress("muons_nTkIsoR03",&_muons_nTkIsoR03);
  myChain->SetBranchAddress("muons_nTkIsoR05",&_muons_nTkIsoR05);
  myChain->SetBranchAddress("muons_tkIsoR03",&_muons_tkIsoR03);
  myChain->SetBranchAddress("muons_tkIsoR05",&_muons_tkIsoR05);
  myChain->SetBranchAddress("muons_emIsoR03",&_muons_emIsoR03);
  myChain->SetBranchAddress("muons_emIsoR05",&_muons_emIsoR05);
  myChain->SetBranchAddress("muons_hadIsoR03",&_muons_hadIsoR03);
  myChain->SetBranchAddress("muons_hadIsoR05",&_muons_hadIsoR05);
  
  //muons TIP/LIP/IP
  myChain->SetBranchAddress("muons_Tip",&muons_Tip);
  myChain->SetBranchAddress("muons_Lip",&muons_Lip);
  myChain->SetBranchAddress("muons_STip",&muons_STip);
  myChain->SetBranchAddress("muons_SLip",&muons_SLip);
  myChain->SetBranchAddress("muons_TipSignif",&muons_TipSignif);
  myChain->SetBranchAddress("muons_LipSignif",&muons_LipSignif);
  myChain->SetBranchAddress("muons_Significance3D",&muons_Significance3D);
  myChain->SetBranchAddress("muons_Value3D",&muons_Value3D);
  myChain->SetBranchAddress("muons_Error3D",&muons_Error3D);
  //muonID variables
  myChain->SetBranchAddress("muons_trkDxy",&_muons_trkDxy);
  myChain->SetBranchAddress("muons_trkDxyError",&_muons_trkDxyError);
  myChain->SetBranchAddress("muons_trkDxyB",&_muons_trkDxyB);
  myChain->SetBranchAddress("muons_trkDz",&_muons_trkDz);
  myChain->SetBranchAddress("muons_trkDzError",&_muons_trkDzError);
  myChain->SetBranchAddress("muons_trkDzB",&_muons_trkDzB); 
  myChain->SetBranchAddress("muons_trkChi2PerNdof",&_muons_trkChi2PerNdof);
  myChain->SetBranchAddress("muons_trkCharge",&_muons_trkCharge);
  myChain->SetBranchAddress("muons_trkNHits",&_muons_trkNHits);
  myChain->SetBranchAddress("muons_trkNPixHits",&_muons_trkNPixHits);
  myChain->SetBranchAddress("muons_trkmuArbitration",&_muons_trkmuArbitration);
  myChain->SetBranchAddress("muons_trkmu2DCompatibilityLoose",&_muons_trkmu2DCompatibilityLoose);
  myChain->SetBranchAddress("muons_trkmu2DCompatibilityTight",&_muons_trkmu2DCompatibilityTight);
  myChain->SetBranchAddress("muons_trkmuOneStationLoose",&_muons_trkmuOneStationLoose);
  myChain->SetBranchAddress("muons_trkmuOneStationTight",&_muons_trkmuOneStationTight);
  myChain->SetBranchAddress("muons_trkmuLastStationLoose",&_muons_trkmuLastStationLoose);
  myChain->SetBranchAddress("muons_trkmuLastStationTight",&_muons_trkmuLastStationTight);
  myChain->SetBranchAddress("muons_trkmuOneStationAngLoose",&_muons_trkmuOneStationAngLoose);
  myChain->SetBranchAddress("muons_trkmuOneStationAngTight",&_muons_trkmuOneStationAngTight);
  myChain->SetBranchAddress("muons_trkmuLastStationAngLoose",&_muons_trkmuLastStationAngLoose);
  myChain->SetBranchAddress("muons_trkmuLastStationAngTight",&_muons_trkmuLastStationAngTight);
  myChain->SetBranchAddress("muons_trkmuLastStationOptimizedLowPtLoose",&_muons_trkmuLastStationOptimizedLowPtLoose);
  myChain->SetBranchAddress("muons_trkmuLastStationOptimizedLowPtTight",&_muons_trkmuLastStationOptimizedLowPtTight);
  myChain->SetBranchAddress("muons_caloCompatibility",&_muons_caloCompatibility);
  myChain->SetBranchAddress("muons_segmentCompatibility",&_muons_segmentCompatibility);
  myChain->SetBranchAddress("muons_glbmuPromptTight",&_muons_glbmuPromptTight);	
  
  // Taus	
  myChain->SetBranchAddress("hpsTau_N", &_nhpsTau);
  myChain->SetBranchAddress("hpsTau_eta", &_hpsTauEta);
  myChain->SetBranchAddress("hpsTau_phi", &_hpsTauPhi);
  myChain->SetBranchAddress("hpsTau_pt", &_hpsTauPt);

  myChain->SetBranchAddress("hpsTau_jet_pt", &_hpsTauJetPt);
  myChain->SetBranchAddress("hpsTau_leadPion_pt", &_hpsTauLeadPionPt);
  myChain->SetBranchAddress("hpsTau_leadTrack_pt", &_hpsTauLeadTrackPt);
  myChain->SetBranchAddress("hpsTau_charge", &_hpsTauCharge);
  myChain->SetBranchAddress("hpsTau_chIso", &_hpsTauChargedIso);
  myChain->SetBranchAddress("hpsTau_phIso", &_hpsTauPhotonsIso);
  myChain->SetBranchAddress("hpsTau_decayMode", &_hpsTauDiscrByDecMode);
  myChain->SetBranchAddress("hpsTau_isoL", &_hpsTauDiscrByLooseIso);
  myChain->SetBranchAddress("hpsTau_isoM", &_hpsTauDiscrByMediumIso);
  myChain->SetBranchAddress("hpsTau_antiMuL", &_hpsTauDiscrAgainstMuonLoose);
  myChain->SetBranchAddress("hpsTau_antiMuT", &_hpsTauDiscrAgainstMuonTight);
  myChain->SetBranchAddress("hpsTau_antiElL", &_hpsTauDiscrAgainstElecLoose);
  myChain->SetBranchAddress("hpsTau_anitElM", &_hpsTauDiscrAgainstElecMedium);
  myChain->SetBranchAddress("hpsTau_antiElT", &_hpsTauDiscrAgainstElecTight);
  //for new tauId
  if(NewTauID){
  myChain->SetBranchAddress("hpsTau_vx", &_hpsTauVx);
  myChain->SetBranchAddress("hpsTau_vy", &_hpsTauVy);
  myChain->SetBranchAddress("hpsTau_vz", &_hpsTauVz);
  myChain->SetBranchAddress("hpsTau_isoL_DBSumPtCorr", &_hpsTauDiscrByLooseIsoDBSumPtCorr);
  myChain->SetBranchAddress("hpsTau_isoM_DBSumPtCorr", &_hpsTauDiscrByMediumIsoDBSumPtCorr);
  myChain->SetBranchAddress("hpsTau_isoT_DBSumPtCorr", &_hpsTauDiscrByTightIsoDBSumPtCorr);
  myChain->SetBranchAddress("hpsTau_CombinedisoL_DBSumPtCorr", &_hpsTauDiscrByLooseCombinedIsoDBSumPtCorr);
  myChain->SetBranchAddress("hpsTau_CombinedisoM_DBSumPtCorr", &_hpsTauDiscrByMediumCombinedIsoDBSumPtCorr);
  myChain->SetBranchAddress("hpsTau_CombinedisoT_DBSumPtCorr", &_hpsTauDiscrByTightCombinedIsoDBSumPtCorr);
  }
  // MET
  myChain->SetBranchAddress("met_pf_et",&_met_pf_et);
  myChain->SetBranchAddress("met_pf_px",&_met_pf_px);
  myChain->SetBranchAddress("met_pf_py",&_met_pf_py);
  myChain->SetBranchAddress("met_pf_phi",&_met_pf_phi);
  myChain->SetBranchAddress("met_pf_set",&_met_pf_set);
  myChain->SetBranchAddress("met_pf_sig",&_met_pf_sig);

  
  myChain->SetBranchAddress("muons", &m_muons);

  m_muons = new TClonesArray ("TLorentzVector");
   
  string FileName = ""; 
if(DoMC){
  FileName = Form("./%s/%sTau%i_%s_%s_%s_%s",dataset.c_str(),isolation.c_str(),ptTresh,dataset.c_str(),Channel.c_str(),MCString.c_str(),MassString.c_str());
}
else{
  FileName = Form("./%s/%sTau%i_%s_%s_%s",dataset.c_str(),isolation.c_str(),ptTresh,dataset.c_str(),Channel.c_str(),MassString.c_str());
}
cout<<"creating file : "<<FileName<<".root"<<endl;
  TFile *f = new TFile(Form("%s.root",FileName.c_str()),"RECREATE");
  //TFile *f = new TFile("output.root","RECREATE");  
  TTree* treeTnP = new TTree("treeTnP", "treeTnP");
  float t_mass;
  float t_etraw,t_et, t_pt, t_eta,t_et_tag, t_pt_tag, t_eta_tag , t_weight;
  int t_match, t_L1match, t_Run_match, t_Run_nomatch, t_trig_HLT_N , t_NRun, t_NEvent, t_Index;
  treeTnP->Branch("mass",&t_mass, "mass/F");
  treeTnP->Branch("etraw",&t_etraw, "etraw/F");
  treeTnP->Branch("et",&t_et, "et/F");
  treeTnP->Branch("pt",&t_pt, "pt/F");
  treeTnP->Branch("eta" ,&t_eta,  "eta/F");
  treeTnP->Branch("et_tag",&t_et_tag, "et_tag/F");
  treeTnP->Branch("pt_tag",&t_pt_tag, "pt_tag/F");
  treeTnP->Branch("eta_tag" ,&t_eta_tag,  "eta_tag/F");
  treeTnP->Branch("match" ,&t_match,  "match/I");
  treeTnP->Branch("L1match" ,&t_L1match,  "L1match/I");
  treeTnP->Branch("weight" ,&t_weight,  "weight/F");
  treeTnP->Branch("Run_match" ,&t_Run_match,  "Run_match/I");
  treeTnP->Branch("Run_nomatch" ,&t_Run_nomatch,  "Run_nomatch/I");
  treeTnP->Branch("NRun" ,&t_NRun,  "NRun/I");
  treeTnP->Branch("NEvent" ,&t_NEvent,  "NEvent/I");
  treeTnP->Branch("trig_HLT_N", &t_trig_HLT_N, "trig_HLT_N/I");
  treeTnP->Branch("Index", &t_Index, "Index/I");
  treeTnP->Branch("vtx_N", &vtx_N, "vtx_N/I");

  std::cout<<"getting entries"<<std::endl;
  int numEntries = myChain->GetEntries () ;
  std::cout<<numEntries<<" entries"<<std::endl;
  
  ///////////////////////////////////////////////////////////
  
  // -------------------------------------------------------------------------------
  // JSON FILE READER
  // -------------------------------------------------------------------------------
  std::string jsonFile = "Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON.txt" ;
  //std::string jsonFile = "Cert_160404-167913_7TeV_PromptReco_Collisions11_JSON.txt" ;
  //std::string jsonFile = "Cert_160404-173692_7TeV_PromptReco_Collisions11_JSON.txt" ;
  //   std::string jsonFile = "Cert_160404-176023_7TeV_PromptReco_Collisions11_JSON.txt" ;
  //std::string jsonFile = "Cert_173236-176023_7TeV_PromptReco_Collisions11_JSON.txt" ;
  //std::string jsonFile = "Cert_175832-178078_7TeV_PromptReco_Collisions11_JSON.txt" ;
  //std::string jsonFile = "Cert_160404-177053_7TeV_PromptReco_Collisions11_JSON.txt" ;
  if (dataset == "2011AMay10ReReco") jsonFile = "Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.txt" ;
  if (dataset == "2011APromptRecov4") jsonFile = "Cert_165071-168437_7TeV_PromptReco_Collisions11_JSON.txt" ;
  if (dataset == "2011AAug05ReReco") jsonFile = "Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v2.txt";
  if (dataset == "2011APromptRecov6")jsonFile = "Cert_160404-177053_7TeV_PromptReco_Collisions11_JSON.txt" ;
  if (dataset == "2011BPromptRecov1")jsonFile = "Cert_175832-180252_7TeV_PromptReco_Collisions11_JSON.txt" ;

  std::map<int, std::vector<std::pair<int, int> > > jsonMap = readJSONFile(jsonFile); 

//Pile up reweighting
  float  *weightsPU  ;
  if (DoMC) {
    //weightsPU=LumiReWeighting("pudist_MC.root","pudist.root","h") ;
    weightsPU = LumiReWeighting() ;
    //for (int i=0 ; i<51 ; i++) cout<<weightsPU[i]<<endl ;
  }

  //MC reweighting
  float  weights = 1. ;
  if (DoMC) {
    if(DY)weights = 1734*3048.0/15800000; //1796000,842000; //LumiReWeighting("pudist_MC.root","pudist.root","h") ;
    if(WLnu)weights =1734*31314.0/15960000;   //1776000,894000; 
    if(QCD)weights =1734*84679.3/19700000;   //1580000,618000  ;
    if(TTJets)weights =1734*157.5/3760000;   //1780000,888000; 
  }
  
  f->cd() ;

  TRandom x ; 
  int nStored(0) ;

  for (int iEvent = 0 ; iEvent <numEntries ; ++iEvent ) { //iEvent <numEntries ; ++iEvent){	 

    if (iEvent%100000 == 0) cout<<iEvent<<" processed"<<endl ;
    
    myChain->GetEntry (iEvent) ;
    
    if (iEvent == 0) {
      cout<<"=============================="<<endl ;
      cout<<"Triggers envisaged: "<<endl ;
       for (int itr=0; itr<trig_HLT_algoStudied->size() ; ++itr) {
	 cout<<(*trig_HLT_algoStudied)[itr]<<" , index = "<<itr<<endl ;
       }
      cout<<"=============================="<<endl ;
      cout<<"Filters envisaged: "<<endl ;
       for (int itr=0; itr<trig_HLT_filterStudied->size() ; ++itr) {
	 cout<<(*trig_HLT_filterStudied)[itr]<<" , index = "<<itr<<endl ;
       }
      cout<<"=============================="<<endl ;

      cout<<endl ;
      cout<<"=============================="<<endl ;
      cout<<"Selecting events with trigger: "<< (*trig_HLT_algoStudied)[idtriggersource] <<endl;
      cout<<"Tag is matched to filter: "<<(*trig_HLT_filterStudied)[idtagfilter]<<endl ;
      cout<<"Probe is trying to be matched to filter: "<<(*trig_HLT_filterStudied)[idprobefilter]<<endl ;
      cout<<"=============================="<<endl ;
      cout<<endl ;
    }
 
    //if(!MC && nRun<=163261) continue ;
    //if (!MC)if(nRun<165922) continue ;//for tightTau only!!
    //if(nRun>163869) continue ;
    if (dataset == "2011APromptRecov4"&& isolation == "Tight" && !DoMC && nRun<=165922)continue ;
 
    if (dataset == "2011AMay10ReReco"){
	if(!DoMC && nRun<=163261) continue ;
    }
    if (dataset == "2011APromptRecov6"){
      if(isolation == "Medium"){
	if(!DoMC && nRun<=173212) continue ;// PromptReco v6 medium
      }
      if(isolation == "Loose" && ptTresh==15){
	if(!DoMC && nRun>173212) continue ;
      }
    }
    //Vertex numbers
//     if(vtx_N>2)continue;

//     if(vtx_N>5)continue;
//     if(vtx_N<=2)continue;
  
//     if(vtx_N>8)continue;
//     if(vtx_N<=5)continue;

//     if(vtx_N>11)continue;
//     if(vtx_N<=8)continue;

//     if(vtx_N<=11)continue;

    //cout<<"vertex :"<<vtx_N<<endl;
    bool triggersource(false) ;
    if(_muons_N < 1 ) continue;
    if(_nhpsTau < 1) continue;
    if(DoMC){
      for (int itr=0; itr<trig_HLT_NPath; itr++) {
	if (trig_HLT_pathname[itr] == idtriggersource ) triggersource = true ;     // HLT trigger source must be responsible for aquisition
      }
    }
    else {
      for (int itr=0; itr<trig_HLT_NPath; itr++) {
	if (trig_HLT_pathname[itr] == idtriggersource)triggersource = true ; 
      }
    }

    ////////////RUN SELECTION/////////////////////////////////
    bool isGoodRun = AcceptEventByRunAndLumiSection(nRun, nLumi, jsonMap);
    if(!DoMC && isGoodRun==false) continue;
    ///////////RUN SELECTION//////////////////////////////////
    int goodmuons=0 ;
    for (int mun = 0 ; mun < _muons_N ; ++mun) {

      TLorentzVector * theMuon = (TLorentzVector*) (m_muons->At (mun)) ;
      //Request not to have a second good muon
      if(theMuon->Pt()>15 && 
	 fabs(theMuon->Eta())<2.4 &&
	 _muons_isglobal[mun] && 
	 (_muons_tkIsoR03[mun] + _muons_emIsoR03[mun] + _muons_hadIsoR03[mun])/(theMuon->Pt())<0.3) goodmuons++ ; 
    }
    if(goodmuons>1) continue ;
    
    vector < pair<int,int> > TnP ;
    
    //tags loop
    for (int mun = 0 ; mun < _muons_N ; ++mun) {
      TLorentzVector * theMuon = (TLorentzVector*) (m_muons->At (mun)) ;	  
      
      //////////////
      //Tag cuts must be tighter than online 
      //////////////
      
      if(theMuon->Pt()<17) continue ; //pt higher than L1 trigger
      if(fabs(theMuon->Eta())>2.1) continue ;
      if(!_muons_istracker[mun]) continue ;//
      if(!_muons_isglobal[mun]) continue ;
      if(fabs(_muons_dxyPV[mun])>0.045) continue ;
      if(fabs(_muons_dzPV[mun])>0.2) continue ;
      if(_muons_normalizedChi2[mun]>10) continue ;
      if(_muons_NtrackerHits[mun]<10) continue ;//2
      if(_muons_NpixelHits[mun]<1) continue ;
      if(_muons_NmuonHits[mun]<1) continue ;
      if(_muons_Nmatches[mun]<2) continue ;
      if((_muons_tkIsoR03[mun] + _muons_emIsoR03[mun] + _muons_hadIsoR03[mun])/(theMuon->Pt())>0.1) continue ;//isolation 0.5

      //mT cut
      Double_t mT = sqrt(2*theMuon->Pt()*_met_pf_et*(1-cos(fabs(theMuon->Phi()-_met_pf_phi))));
      if(mT>40) continue ;//<60 
      
      // tag must be associated to trigger object
      bool tagmatch(false) ;
      for (int itr=0; itr<trig_HLT_N; itr++) {
	if (trig_HLT_name[itr] != idtagfilter) continue ; // HLT object must be from a filter belonging to the trigger responsible for aquisition
	if (trig_HLT_pt[itr] < 17)  continue ; // HLT object must be above 17 GeV
	if (deltaR(trig_HLT_eta[itr], trig_HLT_phi[itr],theMuon->Eta(),theMuon->Phi() ) < 0.3) 
	tagmatch = true ;
      }

      if (!tagmatch) continue ;
      
      // probe loop
      for (int j = 0 ; j < _nhpsTau ; ++j) {
	
	//mu-tau different charge
	if(_muons_charge[mun]*_hpsTauCharge[j]>0) continue ; 
	//if(_muons_charge[mun]*_hpsTauCharge[j]<0) continue ; //Same charge to study fakes
	
	TLorentzVector * theTau = new TLorentzVector;
	theTau->SetPtEtaPhiM(_hpsTauPt[j], _hpsTauEta[j], _hpsTauPhi[j],1.777 );
	TLorentzVector total = (*theMuon) + (*theTau);
	double PairMass = total.M();
	
	if(!MassPlot){
	if(PairMass<45) continue ;
	if(PairMass>70) continue ;
	}
	//////////////
	//offline cuts
	//////////////

	if(_hpsTauDiscrByDecMode[j]<0.5) continue ;
	if(fabs(_hpsTauEta[j])>2.3) continue ;
	//if(_hpsTauPt[j]<20) continue ;
	if(MuTau){
	  if(_hpsTauDiscrAgainstMuonTight[j]<0.5) continue ;
	  if(_hpsTauDiscrAgainstElecLoose[j]<0.5) continue ;
	}
 	if(!MuTau){
	  if(_hpsTauDiscrAgainstMuonLoose[j]<0.5) continue ;
	  if(_hpsTauDiscrAgainstElecTight[j]<0.5) continue ;
 	}
	if(!NewTauID){
	  if(_hpsTauDiscrByLooseIso[j]<0.5) continue ;
	}
	
	if(NewTauID){
	  if(_hpsTauDiscrByLooseCombinedIsoDBSumPtCorr[j]<0.5) continue ;
	  if(fabs(_hpsTauVz[j] -  vtx_z[0])>0.2) continue ;
	}
	// avoid tag and probe to be in the same  region (would articficially increase trigger efficiency)
	if (deltaR(theMuon->Eta(),theMuon->Phi(),_hpsTauEta[j] ,_hpsTauPhi[j]) <0.3)continue;// 0.52) continue ; 

	TnP.push_back(pair<int, int>(mun,j)) ;
	
      }  // loop probe
    } // loop tag 
    if (TnP.size()==0) continue ;

    ////////////////////////
    // select good pairs
    ////////////////////////

    // when several tags, choose arbitrary 1 as being source of the trigger.
    int choice = int(x.Rndm() * TnP.size()) ;
    int ref = TnP[choice].first ; // this is the tag which generated trigger
    
    for (int pair = 0 ; pair< TnP.size() ; ++pair) {
      
      int mun = TnP[pair].first ;
      int j = TnP[pair].second ;
      
      if (j==ref) continue ; // avoid bias
      
      TLorentzVector *tag = (TLorentzVector*) (m_muons->At (mun));
      TLorentzVector * probe = new TLorentzVector;
      probe->SetPtEtaPhiM(_hpsTauPt[j], _hpsTauEta[j], _hpsTauPhi[j],1.777 );

      TLorentzVector total = (*tag) + (*probe);
      double mass = total.M();

      // is probe firing trigger?
      bool probematch(false) ;
      for (int itr=0; itr<trig_HLT_N; itr++) {
	// HLT object must be from a filter belonging to the studied trigger
	if (trig_HLT_name[itr] != idprobefilter) continue ; 
	if (deltaR(trig_HLT_eta[itr], trig_HLT_phi[itr],probe->Eta(),probe->Phi()) < 0.3 && trig_HLT_pt[itr]>ptTresh) probematch = true ;
      }//loop trig_HLT_N
            
      
      ////////////
      //CHECKS
      ///////////

      t_Run_match = -1 ;
      t_Run_nomatch = -1 ;
      
      if(probematch) t_Run_match = nRun ;
      else t_Run_nomatch = nRun ;

      nStored++ ;
      t_Index = 5 ;//5 for promptReco
      if(DY)t_Index = 0 ;
      if(WLnu)t_Index = 1 ;
      if(QCD)t_Index = 2 ;
      if(TTJets)t_Index = 3 ;

      t_mass = mass ;
      t_etraw = -1;
      t_et = -1;
      t_L1match = -1;
      t_pt_tag = tag->Pt()  ;
      t_eta_tag = tag->Eta() ;
      t_pt = probe->Pt()  ;
      t_eta = probe->Eta() ;
      t_match = (int)probematch ;
      t_trig_HLT_N = trig_HLT_N ;
      t_NRun = nRun;
      t_NEvent = nEvent;
      t_weight = weights ;
      if (DoMC && PU_N<51) t_weight = weights*weightsPU[PU_N] ;//[PU_N] ;

      treeTnP->Fill() ;
      
    }//loop pair
    
  } // loop entries

  cout<<"evts enregistres "<<nStored<<endl;
  cout<<"weight  "<<t_weight<<endl;
  cout<<"==================================================================================================================================================="<<endl;
  treeTnP->Write() ;
  
  return;
}

void TnP() {

//2011AMay10ReReco
//   tagAndProbeMuonTau("2011AMay10ReReco", "Loose",true,true,true,false,10,0);
//   tagAndProbeMuonTau("2011AMay10ReReco", "Loose",true,false,true,false,10,0);
//   tagAndProbeMuonTau("2011AMay10ReReco", "Loose",true,true,false,false,10,0);
//   tagAndProbeMuonTau("2011AMay10ReReco", "Loose",true,false,false,false,10,0);

//   tagAndProbeMuonTau("2011AMay10ReReco", "Loose",true,true,true,false,15,0);
//   tagAndProbeMuonTau("2011AMay10ReReco", "Loose",true,false,true,false,15,0);
//   tagAndProbeMuonTau("2011AMay10ReReco", "Loose",true,true,false,false,15,0);
//   tagAndProbeMuonTau("2011AMay10ReReco", "Loose",true,false,false,false,15,0);

//   tagAndProbeMuonTau("2011AMay10ReReco", "Loose",true,true,true,false,20,0);
//   tagAndProbeMuonTau("2011AMay10ReReco", "Loose",true,false,true,false,20,0);
//   tagAndProbeMuonTau("2011AMay10ReReco", "Loose",true,true,false,false,20,0);
//   tagAndProbeMuonTau("2011AMay10ReReco", "Loose",true,false,false,false,20,0);


//2011APromptRecov4
//   tagAndProbeMuonTau("2011APromptRecov4", "Loose",true,false,true,false,20,0);
//   tagAndProbeMuonTau("2011APromptRecov4", "Loose",true,false,false,false,20,0);
//  tagAndProbeMuonTau("2011APromptRecov4", "Loose",true,true,true,false,20,0);
//  tagAndProbeMuonTau("2011APromptRecov4", "Loose",true,true,false,false,20,0);

//   tagAndProbeMuonTau("2011APromptRecov4", "Loose",true,false,true,false,15,0);
//   tagAndProbeMuonTau("2011APromptRecov4", "Loose",true,false,false,false,15,0);
//  tagAndProbeMuonTau("2011APromptRecov4", "Loose",true,true,true,false,15,0);
//  tagAndProbeMuonTau("2011APromptRecov4", "Loose",true,true,false,false,15);

  //Tight no savetag for 2011APromptRecov4 !!
  //   tagAndProbeMuonTau("2011APromptRecov4", "Tight",true,false,true,false,20,0);
  //   tagAndProbeMuonTau("2011APromptRecov4", "Tight",true,false,false,false,20,0);
  //   tagAndProbeMuonTau("2011APromptRecov4", "Tight",true,true,true,false,20,0);
  //   tagAndProbeMuonTau("2011APromptRecov4", "Tight",true,true,false,false,20,0);



//2011AAug05ReReco
  //tagAndProbeMuonTau("2011AAug05ReReco", "Loose",true,true,true,false,20,0);
//   tagAndProbeMuonTau("2011AAug05ReReco", "Loose",true,false,true,false,20,0);
//   tagAndProbeMuonTau("2011AAug05ReReco", "Loose",true,true,false,false,20,0);
//   tagAndProbeMuonTau("2011AAug05ReReco", "Loose",true,false,false,false,20,0);

//   tagAndProbeMuonTau("2011AAug05ReReco", "Tight",true,true,true,false,20,0);
//   tagAndProbeMuonTau("2011AAug05ReReco", "Tight",true,false,true,false,20,0);
//   tagAndProbeMuonTau("2011AAug05ReReco", "Tight",true,true,false,false,20,0);
//   tagAndProbeMuonTau("2011AAug05ReReco", "Tight",true,false,false,false,20,0);

//   tagAndProbeMuonTau("2011AAug05ReReco", "Loose",true,true,true,false,15,0);
//   tagAndProbeMuonTau("2011AAug05ReReco", "Loose",true,false,true,false,15,0);
//   tagAndProbeMuonTau("2011AAug05ReReco", "Loose",true,true,false,false,15,0);
//   tagAndProbeMuonTau("2011AAug05ReReco", "Loose",true,false,false,false,15,0);


//2011APromptRecov6
//   tagAndProbeMuonTau("2011APromptRecov6", "Loose",true,true,true,false,20,0);
//   tagAndProbeMuonTau("2011APromptRecov6", "Loose",true,false,true,false,20,0);
//   tagAndProbeMuonTau("2011APromptRecov6", "Loose",true,true,false,false,20,0);
//   tagAndProbeMuonTau("2011APromptRecov6", "Loose",true,false,false,false,20,0);

//   tagAndProbeMuonTau("2011APromptRecov6", "Medium",true,true,true,false,20,0);
//   tagAndProbeMuonTau("2011APromptRecov6", "Medium",true,false,true,false,20,0);
//   tagAndProbeMuonTau("2011APromptRecov6", "Medium",true,true,false,false,20,0);
//   tagAndProbeMuonTau("2011APromptRecov6", "Medium",true,false,false,false,20,0);

//   tagAndProbeMuonTau("2011APromptRecov6", "Tight",true,true,true,false,20,0);
//   tagAndProbeMuonTau("2011APromptRecov6", "Tight",true,false,true,false,20,0);
//   tagAndProbeMuonTau("2011APromptRecov6", "Tight",true,true,false,false,20,0);
//   tagAndProbeMuonTau("2011APromptRecov6", "Tight",true,false,false,false,20,0);

//   tagAndProbeMuonTau("2011APromptRecov6", "Loose",true,true,true,false,15,0);
//   tagAndProbeMuonTau("2011APromptRecov6", "Loose",true,false,true,false,15,0);
//   tagAndProbeMuonTau("2011APromptRecov6", "Loose",true,true,false,false,15,0);
//   tagAndProbeMuonTau("2011APromptRecov6", "Loose",true,false,false,false,15,0);


//Summer11
//MuTauAna
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,true,true,20,0);   
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,true,true,20,1);   
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,true,true,20,2);   
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,true,true,20,3);
  
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,true,true,10,0);   
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,true,true,10,1);   
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,true,true,10,2);   
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,true,true,10,3);
  
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,true,true,15,0);   
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,true,true,15,1);   
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,true,true,15,2);   
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,true,true,15,3);  

//   //ElecTauAna
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,false,true,20,0);   
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,false,true,20,1);   
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,false,true,20,2);   
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,false,true,20,3);
  
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,false,true,10,0);   
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,false,true,10,1);   
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,false,true,10,2);   
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,false,true,10,3);
  
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,false,true,15,0);   
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,false,true,15,1);   
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,false,true,15,2);   
//   tagAndProbeMuonTau("Summer11", "Loose",false,false,false,true,15,3);   

  //MassPlot
  tagAndProbeMuonTau("Summer11", "Loose",false,true,true,true,20,0);   
  tagAndProbeMuonTau("Summer11", "Loose",false,true,true,true,20,1);   
  tagAndProbeMuonTau("Summer11", "Loose",false,true,true,true,20,2);   
  tagAndProbeMuonTau("Summer11", "Loose",false,true,true,true,20,3);

   
   
   
//2011BPromptRecov1
//   tagAndProbeMuonTau("2011BPromptRecov1", "Loose",true,false,true,false,20,0);
//   tagAndProbeMuonTau("2011BPromptRecov1", "Loose",true,false,false,false,20,0);
//   tagAndProbeMuonTau("2011BPromptRecov1", "Loose",true,true,false,false,20,0);
//   tagAndProbeMuonTau("2011BPromptRecov1", "Loose",true,true,true,false,20,0);
//
//   tagAndProbeMuonTau("2011BPromptRecov1", "Medium",true,false,true,false,20,0);
//   tagAndProbeMuonTau("2011BPromptRecov1", "Medium",true,false,false,false,20,0);
//   tagAndProbeMuonTau("2011BPromptRecov1", "Medium",true,true,false,false,20,0);
//   tagAndProbeMuonTau("2011BPromptRecov1", "Medium",true,true,true,false,20,0);
//
//   tagAndProbeMuonTau("2011BPromptRecov1", "Tight",true,false,true,false,20,0);
//   tagAndProbeMuonTau("2011BPromptRecov1", "Tight",true,false,false,false,20,0);
//   tagAndProbeMuonTau("2011BPromptRecov1", "Tight",true,true,false,false,20,0);
//   tagAndProbeMuonTau("2011BPromptRecov1", "Tight",true,true,true,false,20,0);


}

#include "Bianchi/TauTauStudies/interface/TausUserEmbedded.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#define DEBUG false

TausUserEmbedded::TausUserEmbedded(const edm::ParameterSet & iConfig){

  tauTag_    = iConfig.getParameter<edm::InputTag>("tauTag");
  vertexTag_ = iConfig.getParameter<edm::InputTag>("vertexTag");
  electronTag_ = iConfig.getParameter<edm::InputTag>("electronTag");

  produces<pat::TauCollection>("");

}

TausUserEmbedded::~TausUserEmbedded(){
}

void TausUserEmbedded::produce(edm::Event & iEvent, const edm::EventSetup & iSetup){

  edm::Handle<pat::TauCollection> tausHandle;
  iEvent.getByLabel(tauTag_,tausHandle);
  const pat::TauCollection* taus = tausHandle.product();

  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByLabel(vertexTag_,vertexHandle);
  const reco::VertexCollection* vertexes = vertexHandle.product();

//   edm::Handle<pat::ElectronCollection> electronsHandle;
//   iEvent.getByLabel(electronTag_,electronsHandle);
//   const pat::ElectronCollection* electrons = electronsHandle.product();
  edm::Handle<reco::GsfElectronCollection> gsfElectronsHandle;
  iEvent.getByLabel("gsfElectrons",gsfElectronsHandle);
  const reco::GsfElectronCollection* gsfElectrons = gsfElectronsHandle.product();

  std::auto_ptr< pat::TauCollection > tausUserEmbeddedColl( new pat::TauCollection() ) ;

  for(unsigned int i = 0; i < taus->size(); i++){
    pat::Tau aTau( (*taus)[i] );

    bool matchElectronCutsVeto = false;
    if(DEBUG){
      std::cout<<std::endl;
      std::cout<<"Number of GsfElectrons: "<<gsfElectrons->size()<<std::endl;
    }

    for(unsigned int i = 0; i < gsfElectrons->size(); i++){
//       pat::Electron aElectron( (*electrons)[i] );
//       const reco::GsfElectron* aGsf = static_cast<reco::GsfElectron*>(&aElectron); 
      reco::GsfElectron aElectron( (*gsfElectrons)[i] );
      
      const reco::Track *el_track = (const reco::Track*)((aElectron).gsfTrack().get());  
      const reco::HitPattern& p_inner = el_track->trackerExpectedHitsInner(); 
      float nHits = p_inner.numberOfHits();
      float dPhi  = fabs(aElectron.deltaPhiSuperClusterTrackAtVtx());
      float dEta  = fabs(aElectron.deltaEtaSuperClusterTrackAtVtx());
      float sihih = aElectron.sigmaIetaIeta();
      float HoE   = aElectron.hadronicOverEm();
      if(DEBUG){
	std::cout<<"GsfElectron: "<<i<<std::endl;
	std::cout<<"GsfElectron nHits: "<<nHits<<std::endl;
	std::cout<<"GsfElectron dPhi: "<<dPhi<<std::endl;
	std::cout<<"GsfElectron dEta: "<<dEta<<std::endl;
	std::cout<<"GsfElectron sihih: "<<sihih<<std::endl;
	std::cout<<"GsfElectron HoE: "<<HoE<<std::endl;
      }
      bool ElectronPassCutsVeto = false;



      if((nHits<=999) &&
	 ((fabs(aElectron.eta())<1.5) &&
	  (sihih < 0.010) &&
	  (dPhi < 0.80) &&
	  (dEta < 0.007) &&
	  (HoE < 0.15)) ||
	 ((fabs(aElectron.eta())>1.5) && (fabs(aElectron.eta())<2.3) &&
	  (sihih < 0.030) &&
	  (dPhi < 0.70) &&
	  (dEta < 0.010) &&
	  (HoE < 999))
	 ) ElectronPassCutsVeto = true ;
      if(DEBUG)std::cout<<"GsfElectron Pass cuts: "<<ElectronPassCutsVeto<<std::endl;
      if (Geom::deltaR(aElectron.p4(),aTau.p4())<0.3 && ElectronPassCutsVeto){
	matchElectronCutsVeto = true;
      } 
      if(DEBUG){
	if(matchElectronCutsVeto)std::cout<<"DeltaR match: "<<Geom::deltaR(aElectron.p4(),aTau.p4())<<std::endl;
	if(!matchElectronCutsVeto)std::cout<<"DeltaR no match: "<<Geom::deltaR(aElectron.p4(),aTau.p4())<<std::endl;
      }
    }
    float matchElectronCutsVetoFloat = 0;
    if(matchElectronCutsVeto)matchElectronCutsVetoFloat = 1;
    if(DEBUG)std::cout<<"Tau matchElectronCutsVeto : "<<matchElectronCutsVetoFloat<<std::endl;
    aTau.addUserFloat("matchElectronCutsVeto", matchElectronCutsVetoFloat );
    
    float dZPV = vertexes->size()>0 ?
      fabs( aTau.vertex().z() - (*vertexes)[0].position().z() ) : -99; 
    aTau.addUserFloat("dzWrtPV", dZPV );

    float TauCat = aTau.tauID("againstElectronMVA2category");
    aTau.addUserFloat("againstElectronMVA2category", TauCat );

    tausUserEmbeddedColl->push_back(aTau);

  }


  iEvent.put( tausUserEmbeddedColl );
  return;
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TausUserEmbedded);



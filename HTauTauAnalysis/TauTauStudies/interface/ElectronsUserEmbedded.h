#ifndef Bianchi_TauTauStudies_ElectronsUserEmbedded_h
#define Bianchi_TauTauStudies_ElectronsUserEmbedded_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"

#include "HiggsAnalysis/HiggsToWW2Leptons/interface/ElectronIDMVA.h"
#include "Bianchi/Utilities/interface/ElectronMVAEstimator.h"


class ElectronsUserEmbedded : public edm::EDProducer{


 public: 

  explicit ElectronsUserEmbedded(const edm::ParameterSet&);
  virtual ~ElectronsUserEmbedded();

 private:

  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);


  edm::InputTag electronTag_;
  edm::InputTag vertexTag_;
  bool isMC_;
  bool doMVAMIT_;
  bool doMVADaniele_;
  ElectronIDMVA* fMVA_;
  ElectronMVAEstimator* fMVADaniele_;

};


#endif

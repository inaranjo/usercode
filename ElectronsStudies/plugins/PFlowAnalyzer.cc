#include "IvoNaranjo/ElectronsStudies/interface/PFlowAnalyzer.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/Common/interface/Handle.h"


PFlowAnalyzer::PFlowAnalyzer(const edm::ParameterSet& cfg)
  :  plotsAllEta_(0)
{ 
  srcGsfElectrons_  = cfg.getParameter<edm::InputTag>("srcGsfElectrons");
  debug_  = cfg.getParameter<bool>("debug");
} 

PFlowAnalyzer::~PFlowAnalyzer()
{
  delete plotsAllEta_;
}

void PFlowAnalyzer::beginJob()
{ 
  plotsAllEta_       = new plotEntryType("PFMVAInput",-0.1, 9.9);
  plotsAllEta_->bookHistograms();
}


void PFlowAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<reco::GsfElectronCollection> gsfElectrons;
  evt.getByLabel(srcGsfElectrons_, gsfElectrons);
  
  plotsAllEta_->fillHistograms(*gsfElectrons,debug_);
   
}

void PFlowAnalyzer::endJob()
{
// nothing to be done yet...
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PFlowAnalyzer);






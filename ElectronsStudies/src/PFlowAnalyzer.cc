#include "Bianchi/TauTauStudies/plugins/PFlowAnalyzer.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/Common/interface/Handle.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

PFlowAnalyzer::PFlowAnalyzer(const edm::ParameterSet& cfg)
{
  srcPFCandidates_  = cfg.getParameter<edm::InputTag>("srcPFCandidates");
} 

PFlowAnalyzer::~PFlowAnalyzer()
{
  
}

void PFlowAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;
  h1_ = fs->make<TH1F>("h1","h1",150,0,5);
  //h1_ = new plotEntryType();
  //h1_ = fs->make<TH1F>("h1","h1",150,0,100);
  }



void PFlowAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  evt.getByLabel(srcPFCandidates_, pfCandidates);

  float mvaFbrem = -99;

  for ( reco::PFCandidateCollection::const_iterator pfCandidate = pfCandidates->begin();
	pfCandidate != pfCandidates->end(); ++pfCandidate ) {
    reco::PFCandidateElectronExtraRef pfElectron = pfCandidate->electronExtraRef();
    if (pfElectron.isNonnull()){
      mvaFbrem = pfElectron->mvaVariable( reco::PFCandidateElectronExtra::MVA_Fbrem);
      std::cout<<"   "<<mvaFbrem<<std::endl;
      h1_->Fill(mvaFbrem);
    }
  }
  
//   h1_->bookHistograms();
//   h1_->fillHistograms(*pfCandidates);
   
   
}

void PFlowAnalyzer::endJob()
{
// nothing to be done yet...
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PFlowAnalyzer);






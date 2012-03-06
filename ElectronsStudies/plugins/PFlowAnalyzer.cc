#include "IvoNaranjo/ElectronsStudies/interface/PFlowAnalyzer.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


PFlowAnalyzer::PFlowAnalyzer(const edm::ParameterSet& cfg)
  :  plotsAll_(0),
     plotsNumPV0to10_Barrel_Pt20to30_(0),
     plotsNumPV10to20_Barrel_Pt20to30_(0),
     plotsNumPV20to30_Barrel_Pt20to30_(0),
     plotsNumPV30to40_Barrel_Pt20to30_(0),
     plotsNumPV0to10_Endcap_Pt20to30_(0),
     plotsNumPV10to20_Endcap_Pt20to30_(0),
     plotsNumPV20to30_Endcap_Pt20to30_(0),
     plotsNumPV30to40_Endcap_Pt20to30_(0),
     plotsNumPV0to10_Barrel_Pt30to40_(0),
     plotsNumPV10to20_Barrel_Pt30to40_(0),
     plotsNumPV20to30_Barrel_Pt30to40_(0),
     plotsNumPV30to40_Barrel_Pt30to40_(0),
     plotsNumPV0to10_Endcap_Pt30to40_(0),
     plotsNumPV10to20_Endcap_Pt30to40_(0),
     plotsNumPV20to30_Endcap_Pt30to40_(0),
     plotsNumPV30to40_Endcap_Pt30to40_(0),
     plotsNumPV0to10_Barrel_Pt40to50_(0),
     plotsNumPV10to20_Barrel_Pt40to50_(0),
     plotsNumPV20to30_Barrel_Pt40to50_(0),
     plotsNumPV30to40_Barrel_Pt40to50_(0),
     plotsNumPV0to10_Endcap_Pt40to50_(0),
     plotsNumPV10to20_Endcap_Pt40to50_(0),
     plotsNumPV20to30_Endcap_Pt40to50_(0),
     plotsNumPV30to40_Endcap_Pt40to50_(0)
{ 
  srcPrimaryVertex_  = cfg.getParameter<edm::InputTag>("srcPrimaryVertex");
  srcGsfElectrons_  = cfg.getParameter<edm::InputTag>("srcGsfElectrons");
  srcPFTaus_  = cfg.getParameter<edm::InputTag>("srcPFTaus");
  srcGenElectrons_  = cfg.getParameter<edm::InputTag>("srcGenElectrons");
  srcGenTaus_  = cfg.getParameter<edm::InputTag>("srcGenTaus");
  match_  = cfg.getParameter<int>("match");
  debug_  = cfg.getParameter<bool>("debug");
} 

PFlowAnalyzer::~PFlowAnalyzer()
{
  delete plotsAll_;
  delete plotsNumPV0to10_Barrel_Pt20to30_;
  delete plotsNumPV10to20_Barrel_Pt20to30_;
  delete plotsNumPV20to30_Barrel_Pt20to30_;
  delete plotsNumPV30to40_Barrel_Pt20to30_;
  delete plotsNumPV0to10_Endcap_Pt20to30_;
  delete plotsNumPV10to20_Endcap_Pt20to30_;
  delete plotsNumPV20to30_Endcap_Pt20to30_;
  delete plotsNumPV30to40_Endcap_Pt20to30_;
  delete plotsNumPV0to10_Barrel_Pt30to40_;
  delete plotsNumPV10to20_Barrel_Pt30to40_;
  delete plotsNumPV20to30_Barrel_Pt30to40_;
  delete plotsNumPV30to40_Barrel_Pt30to40_;
  delete plotsNumPV0to10_Endcap_Pt30to40_;
  delete plotsNumPV10to20_Endcap_Pt30to40_;
  delete plotsNumPV20to30_Endcap_Pt30to40_;
  delete plotsNumPV30to40_Endcap_Pt30to40_;
  delete plotsNumPV0to10_Barrel_Pt40to50_;
  delete plotsNumPV10to20_Barrel_Pt40to50_;
  delete plotsNumPV20to30_Barrel_Pt40to50_;
  delete plotsNumPV30to40_Barrel_Pt40to50_;
  delete plotsNumPV0to10_Endcap_Pt40to50_;
  delete plotsNumPV10to20_Endcap_Pt40to50_;
  delete plotsNumPV20to30_Endcap_Pt40to50_;
  delete plotsNumPV30to40_Endcap_Pt40to50_;
}

void PFlowAnalyzer::beginJob()
{ 
  plotsAll_       = new plotEntryType("All",0.0, 3.0, 0, 40, 0, 60);
  plotsNumPV0to10_Barrel_Pt20to30_ = new plotEntryType("NumPV0to10_Barrel_Pt20to30",0.0,1.479,0,10,0,60);
  plotsNumPV10to20_Barrel_Pt20to30_ = new plotEntryType("NumPV10to20_Barrel_Pt20to30",0.0,1.479,10,20,0,60);
  plotsNumPV20to30_Barrel_Pt20to30_ = new plotEntryType("NumPV20to30_Barrel_Pt20to30",0.0,1.479,20,30,0,60);
  plotsNumPV30to40_Barrel_Pt20to30_ = new plotEntryType("NumPV30to40_Barrel_Pt20to30",0.0,1.479,30,40,0,60);
  plotsNumPV0to10_Endcap_Pt20to30_ = new plotEntryType("NumPV0to10_Endcap_Pt20to30",1.479,2.3,0,10,0,60);
  plotsNumPV10to20_Endcap_Pt20to30_ = new plotEntryType("NumPV10to20_Endcap_Pt20to30",1.479,2.3,10,20,0,60);
  plotsNumPV20to30_Endcap_Pt20to30_ = new plotEntryType("NumPV20to30_Endcap_Pt20to30",1.479,2.3,20,30,0,60);
  plotsNumPV30to40_Endcap_Pt20to30_ = new plotEntryType("NumPV30to40_Endcap_Pt20to30",1.479,2.3,30,40,0,60);
  plotsNumPV0to10_Barrel_Pt30to40_ = new plotEntryType("NumPV0to10_Barrel_Pt30to40",0.0,1.479,0,10,0,60);
  plotsNumPV10to20_Barrel_Pt30to40_ = new plotEntryType("NumPV10to20_Barrel_Pt30to40",0.0,1.479,10,20,0,60);
  plotsNumPV20to30_Barrel_Pt30to40_ = new plotEntryType("NumPV20to30_Barrel_Pt30to40",0.0,1.479,20,30,0,60);
  plotsNumPV30to40_Barrel_Pt30to40_ = new plotEntryType("NumPV30to40_Barrel_Pt30to40",0.0,1.479,30,40,0,60);
  plotsNumPV0to10_Endcap_Pt30to40_ = new plotEntryType("NumPV0to10_Endcap_Pt30to40",1.479,2.3,0,10,0,60);
  plotsNumPV10to20_Endcap_Pt30to40_ = new plotEntryType("NumPV10to20_Endcap_Pt30to40",1.479,2.3,10,20,0,60);
  plotsNumPV20to30_Endcap_Pt30to40_ = new plotEntryType("NumPV20to30_Endcap_Pt30to40",1.479,2.3,20,30,0,60);
  plotsNumPV30to40_Endcap_Pt30to40_ = new plotEntryType("NumPV30to40_Endcap_Pt30to40",1.479,2.3,30,40,0,60);
  plotsNumPV0to10_Barrel_Pt40to50_ = new plotEntryType("NumPV0to10_Barrel_Pt40to50",0.0,1.479,0,10,0,60);
  plotsNumPV10to20_Barrel_Pt40to50_ = new plotEntryType("NumPV10to20_Barrel_Pt40to50",0.0,1.479,10,20,0,60);
  plotsNumPV20to30_Barrel_Pt40to50_ = new plotEntryType("NumPV20to30_Barrel_Pt40to50",0.0,1.479,20,30,0,60);
  plotsNumPV30to40_Barrel_Pt40to50_ = new plotEntryType("NumPV30to40_Barrel_Pt40to50",0.0,1.479,30,40,0,60);
  plotsNumPV0to10_Endcap_Pt40to50_ = new plotEntryType("NumPV0to10_Endcap_Pt40to50",1.479,2.3,0,10,0,60);
  plotsNumPV10to20_Endcap_Pt40to50_ = new plotEntryType("NumPV10to20_Endcap_Pt40to50",1.479,2.3,10,20,0,60);
  plotsNumPV20to30_Endcap_Pt40to50_ = new plotEntryType("NumPV20to30_Endcap_Pt40to50",1.479,2.3,20,30,0,60);
  plotsNumPV30to40_Endcap_Pt40to50_ = new plotEntryType("NumPV30to40_Endcap_Pt40to50",1.479,2.3,30,40,0,60);

  plotsAll_->bookHistograms(0.0, 3.0,0,40, 0, 60);
  plotsNumPV0to10_Barrel_Pt20to30_->bookHistograms(0.0,1.479,0,10, 0, 60);
  plotsNumPV10to20_Barrel_Pt20to30_->bookHistograms(0.0,1.479,10,20, 0, 60);
  plotsNumPV20to30_Barrel_Pt20to30_->bookHistograms(0.0,1.479,20,30, 0, 60);
  plotsNumPV30to40_Barrel_Pt20to30_->bookHistograms(0.0,1.479,30,40, 0, 60);
  plotsNumPV0to10_Endcap_Pt20to30_->bookHistograms(1.479,2.3,0,10, 0, 60);
  plotsNumPV10to20_Endcap_Pt20to30_->bookHistograms(1.479,2.3,10,20, 0, 60);
  plotsNumPV20to30_Endcap_Pt20to30_->bookHistograms(1.479,2.3,20,30, 0, 60);
  plotsNumPV30to40_Endcap_Pt20to30_->bookHistograms(1.479,2.3,30,40, 0, 60);
  plotsNumPV0to10_Barrel_Pt30to40_->bookHistograms(0.0,1.479,0,10, 0, 60);
  plotsNumPV10to20_Barrel_Pt30to40_->bookHistograms(0.0,1.479,10,20, 0, 60);
  plotsNumPV20to30_Barrel_Pt30to40_->bookHistograms(0.0,1.479,20,30, 0, 60);
  plotsNumPV30to40_Barrel_Pt30to40_->bookHistograms(0.0,1.479,30,40, 0, 60);
  plotsNumPV0to10_Endcap_Pt30to40_->bookHistograms(1.479,2.3,0,10, 0, 60);
  plotsNumPV10to20_Endcap_Pt30to40_->bookHistograms(1.479,2.3,10,20, 0, 60);
  plotsNumPV20to30_Endcap_Pt30to40_->bookHistograms(1.479,2.3,20,30, 0, 60);
  plotsNumPV30to40_Endcap_Pt30to40_->bookHistograms(1.479,2.3,30,40, 0, 60);
  plotsNumPV0to10_Barrel_Pt40to50_->bookHistograms(0.0,1.479,0,10, 0, 60);
  plotsNumPV10to20_Barrel_Pt40to50_->bookHistograms(0.0,1.479,10,20, 0, 60);
  plotsNumPV20to30_Barrel_Pt40to50_->bookHistograms(0.0,1.479,20,30, 0, 60);
  plotsNumPV30to40_Barrel_Pt40to50_->bookHistograms(0.0,1.479,30,40, 0, 60);
  plotsNumPV0to10_Endcap_Pt40to50_->bookHistograms(1.479,2.3,0,10, 0, 60);
  plotsNumPV10to20_Endcap_Pt40to50_->bookHistograms(1.479,2.3,10,20, 0, 60);
  plotsNumPV20to30_Endcap_Pt40to50_->bookHistograms(1.479,2.3,20,30, 0, 60);
  plotsNumPV30to40_Endcap_Pt40to50_->bookHistograms(1.479,2.3,30,40, 0, 60);
}


void PFlowAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<reco::GsfElectronCollection> gsfElectrons;
  evt.getByLabel(srcGsfElectrons_, gsfElectrons);

  edm::Handle<reco::PFTauCollection> pfTaus;
  evt.getByLabel(srcPFTaus_, pfTaus);

  typedef edm::View<reco::Candidate> CandidateView;

  edm::Handle<CandidateView> genElectrons;
  evt.getByLabel(srcGenElectrons_, genElectrons);

  edm::Handle<CandidateView> genTaus;
  evt.getByLabel(srcGenTaus_, genTaus);

  edm::Handle<reco::VertexCollection> pVertexes;
  evt.getByLabel(srcPrimaryVertex_, pVertexes);
  
  const reco::VertexCollection* vertexes = pVertexes.product();
  numPV_ = vertexes->size();
  if(debug_){
    std::cout<<"numPV : "<<vertexes->size()<<std::endl;
    for(unsigned int k = 0; k<vertexes->size(); k++){
      std::cout << "Vtx[" << k << "] (x,y,z) = (" << ((*vertexes)[k].position()).x()
		<< "," << ((*vertexes)[k].position()).y() << "," << ((*vertexes)[k].position()).z() << ")"
		<< std::endl;
    }
  }


  plotsAll_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);
  plotsNumPV0to10_Barrel_Pt20to30_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);
  plotsNumPV10to20_Barrel_Pt20to30_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);
  plotsNumPV20to30_Barrel_Pt20to30_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);   
  plotsNumPV30to40_Barrel_Pt20to30_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);
  plotsNumPV0to10_Endcap_Pt20to30_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);
  plotsNumPV10to20_Endcap_Pt20to30_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);
  plotsNumPV20to30_Endcap_Pt20to30_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);   
  plotsNumPV30to40_Endcap_Pt20to30_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);
  plotsNumPV0to10_Barrel_Pt30to40_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);
  plotsNumPV10to20_Barrel_Pt30to40_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);
  plotsNumPV20to30_Barrel_Pt30to40_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);   
  plotsNumPV30to40_Barrel_Pt30to40_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);
  plotsNumPV0to10_Endcap_Pt30to40_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);
  plotsNumPV10to20_Endcap_Pt30to40_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);
  plotsNumPV20to30_Endcap_Pt30to40_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);   
  plotsNumPV30to40_Endcap_Pt30to40_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);
  plotsNumPV0to10_Barrel_Pt40to50_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);
  plotsNumPV10to20_Barrel_Pt40to50_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);
  plotsNumPV20to30_Barrel_Pt40to50_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);   
  plotsNumPV30to40_Barrel_Pt40to50_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);
  plotsNumPV0to10_Endcap_Pt40to50_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);
  plotsNumPV10to20_Endcap_Pt40to50_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);
  plotsNumPV20to30_Endcap_Pt40to50_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);   
  plotsNumPV30to40_Endcap_Pt40to50_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, debug_, numPV_, match_);

}

void PFlowAnalyzer::endJob()
{
// nothing to be done yet...
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PFlowAnalyzer);






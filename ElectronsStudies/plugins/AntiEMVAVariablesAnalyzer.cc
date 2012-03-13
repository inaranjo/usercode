#include "IvoNaranjo/ElectronsStudies/interface/AntiEMVAVariablesAnalyzer.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


AntiEMVAVariablesAnalyzer::AntiEMVAVariablesAnalyzer(const edm::ParameterSet& cfg)
  :  plotsAll_(0),
     plotsAll_Barrel_(0),
     plotsAll_Endcap_(0),
     plotsNumPV0to5_Barrel_Pt20to30_(0),
     plotsNumPV5to10_Barrel_Pt20to30_(0),
     plotsNumPV10to20_Barrel_Pt20to30_(0),
     plotsNumPVOver20_Barrel_Pt20to30_(0),
     plotsNumPV0to5_Endcap_Pt20to30_(0),
     plotsNumPV5to10_Endcap_Pt20to30_(0),
     plotsNumPV10to20_Endcap_Pt20to30_(0),
     plotsNumPVOver20_Endcap_Pt20to30_(0),
     plotsNumPV0to5_Barrel_Pt30to40_(0),
     plotsNumPV5to10_Barrel_Pt30to40_(0),
     plotsNumPV10to20_Barrel_Pt30to40_(0),
     plotsNumPVOver20_Barrel_Pt30to40_(0),
     plotsNumPV0to5_Endcap_Pt30to40_(0),
     plotsNumPV5to10_Endcap_Pt30to40_(0),
     plotsNumPV10to20_Endcap_Pt30to40_(0),
     plotsNumPVOver20_Endcap_Pt30to40_(0),
     plotsNumPV0to5_Barrel_Pt40to50_(0),
     plotsNumPV5to10_Barrel_Pt40to50_(0),
     plotsNumPV10to20_Barrel_Pt40to50_(0),
     plotsNumPVOver20_Barrel_Pt40to50_(0),
     plotsNumPV0to5_Endcap_Pt40to50_(0),
     plotsNumPV5to10_Endcap_Pt40to50_(0),
     plotsNumPV10to20_Endcap_Pt40to50_(0),
     plotsNumPVOver20_Endcap_Pt40to50_(0)
{ 
  srcPrimaryVertex_  = cfg.getParameter<edm::InputTag>("srcPrimaryVertex");
  srcGsfElectrons_  = cfg.getParameter<edm::InputTag>("srcGsfElectrons");
  srcPFTaus_  = cfg.getParameter<edm::InputTag>("srcPFTaus");
  srcGenElectrons_  = cfg.getParameter<edm::InputTag>("srcGenElectrons");
  srcGenTaus_  = cfg.getParameter<edm::InputTag>("srcGenTaus");
  srcGenJets_  = cfg.getParameter<edm::InputTag>("srcGenJets");
  match_  = cfg.getParameter<int>("match");
  debug_  = cfg.getParameter<bool>("debug");
} 

AntiEMVAVariablesAnalyzer::~AntiEMVAVariablesAnalyzer()
{
  delete plotsAll_;
  delete plotsAll_Barrel_;
  delete plotsAll_Endcap_;
  delete plotsNumPV0to5_Barrel_Pt20to30_;
  delete plotsNumPV5to10_Barrel_Pt20to30_;
  delete plotsNumPV10to20_Barrel_Pt20to30_;
  delete plotsNumPVOver20_Barrel_Pt20to30_;
  delete plotsNumPV0to5_Endcap_Pt20to30_;
  delete plotsNumPV5to10_Endcap_Pt20to30_;
  delete plotsNumPV10to20_Endcap_Pt20to30_;
  delete plotsNumPVOver20_Endcap_Pt20to30_;
  delete plotsNumPV0to5_Barrel_Pt30to40_;
  delete plotsNumPV5to10_Barrel_Pt30to40_;
  delete plotsNumPV10to20_Barrel_Pt30to40_;
  delete plotsNumPVOver20_Barrel_Pt30to40_;
  delete plotsNumPV0to5_Endcap_Pt30to40_;
  delete plotsNumPV5to10_Endcap_Pt30to40_;
  delete plotsNumPV10to20_Endcap_Pt30to40_;
  delete plotsNumPVOver20_Endcap_Pt30to40_;
  delete plotsNumPV0to5_Barrel_Pt40to50_;
  delete plotsNumPV5to10_Barrel_Pt40to50_;
  delete plotsNumPV10to20_Barrel_Pt40to50_;
  delete plotsNumPVOver20_Barrel_Pt40to50_;
  delete plotsNumPV0to5_Endcap_Pt40to50_;
  delete plotsNumPV5to10_Endcap_Pt40to50_;
  delete plotsNumPV10to20_Endcap_Pt40to50_;
  delete plotsNumPVOver20_Endcap_Pt40to50_;
}

void AntiEMVAVariablesAnalyzer::beginJob()
{ 
  plotsAll_       = new plotEntryType("All",0.0, 3.0, 0, 50, 10, 60);
  plotsAll_Barrel_       = new plotEntryType("All_Barrel",0.0, 1.479, 0, 50, 10, 60);
  plotsAll_Endcap_       = new plotEntryType("All_Endcap",1.479, 3.0, 0, 50, 10, 60);
  plotsNumPV0to5_Barrel_Pt20to30_ = new plotEntryType("NumPV0to5_Barrel_Pt20to30",0.0,1.479,0,5,20, 30);
  plotsNumPV5to10_Barrel_Pt20to30_ = new plotEntryType("NumPV5to10_Barrel_Pt20to30",0.0,1.479,5,10,20, 30);
  plotsNumPV10to20_Barrel_Pt20to30_ = new plotEntryType("NumPV10to20_Barrel_Pt20to30",0.0,1.479,10,20,20, 30);
  plotsNumPVOver20_Barrel_Pt20to30_ = new plotEntryType("NumPVOver20_Barrel_Pt20to30",0.0,1.479,20,50,20, 30);
  plotsNumPV0to5_Endcap_Pt20to30_ = new plotEntryType("NumPV0to5_Endcap_Pt20to30",1.479,2.3,0,5,20, 30);
  plotsNumPV5to10_Endcap_Pt20to30_ = new plotEntryType("NumPV5to10_Endcap_Pt20to30",1.479,2.3,5,10,20, 30);
  plotsNumPV10to20_Endcap_Pt20to30_ = new plotEntryType("NumPV10to20_Endcap_Pt20to30",1.479,2.3,10,20,20, 30);
  plotsNumPVOver20_Endcap_Pt20to30_ = new plotEntryType("NumPVOver20_Endcap_Pt20to30",1.479,2.3,20,50,20, 30);
  plotsNumPV0to5_Barrel_Pt30to40_ = new plotEntryType("NumPV0to5_Barrel_Pt30to40",0.0,1.479,0,5,30, 40);
  plotsNumPV5to10_Barrel_Pt30to40_ = new plotEntryType("NumPV5to10_Barrel_Pt30to40",0.0,1.479,5,10,30, 40);
  plotsNumPV10to20_Barrel_Pt30to40_ = new plotEntryType("NumPV10to20_Barrel_Pt30to40",0.0,1.479,10,20,30, 40);
  plotsNumPVOver20_Barrel_Pt30to40_ = new plotEntryType("NumPVOver20_Barrel_Pt30to40",0.0,1.479,20,50,30, 40);
  plotsNumPV0to5_Endcap_Pt30to40_ = new plotEntryType("NumPV0to5_Endcap_Pt30to40",1.479,2.3,0,5,30, 40);
  plotsNumPV5to10_Endcap_Pt30to40_ = new plotEntryType("NumPV5to10_Endcap_Pt30to40",1.479,2.3,5,10,30, 40);
  plotsNumPV10to20_Endcap_Pt30to40_ = new plotEntryType("NumPV10to20_Endcap_Pt30to40",1.479,2.3,10,20,30, 40);
  plotsNumPVOver20_Endcap_Pt30to40_ = new plotEntryType("NumPVOver20_Endcap_Pt30to40",1.479,2.3,20,50,30, 40);
  plotsNumPV0to5_Barrel_Pt40to50_ = new plotEntryType("NumPV0to5_Barrel_Pt40to50",0.0,1.479,0,5,40, 50);
  plotsNumPV5to10_Barrel_Pt40to50_ = new plotEntryType("NumPV5to10_Barrel_Pt40to50",0.0,1.479,5,10,40, 50);
  plotsNumPV10to20_Barrel_Pt40to50_ = new plotEntryType("NumPV10to20_Barrel_Pt40to50",0.0,1.479,10,20,40, 50);
  plotsNumPVOver20_Barrel_Pt40to50_ = new plotEntryType("NumPVOver20_Barrel_Pt40to50",0.0,1.479,20,50,40, 50);
  plotsNumPV0to5_Endcap_Pt40to50_ = new plotEntryType("NumPV0to5_Endcap_Pt40to50",1.479,2.3,0,5,40, 50);
  plotsNumPV5to10_Endcap_Pt40to50_ = new plotEntryType("NumPV5to10_Endcap_Pt40to50",1.479,2.3,5,10,40, 50);
  plotsNumPV10to20_Endcap_Pt40to50_ = new plotEntryType("NumPV10to20_Endcap_Pt40to50",1.479,2.3,10,20,40, 50);
  plotsNumPVOver20_Endcap_Pt40to50_ = new plotEntryType("NumPVOver20_Endcap_Pt40to50",1.479,2.3,20,50,40, 50);

  plotsAll_->bookHistograms(0.0, 3.0,0,50, 10, 60);
  plotsAll_Barrel_->bookHistograms(0.0, 1.479,0,50, 10, 60);
  plotsAll_Endcap_->bookHistograms(1.479, 3.0,0,50, 10, 60);
  plotsNumPV0to5_Barrel_Pt20to30_->bookHistograms(0.0,1.479,0,5, 20, 30);
  plotsNumPV5to10_Barrel_Pt20to30_->bookHistograms(0.0,1.479,5,10, 20, 30);
  plotsNumPV10to20_Barrel_Pt20to30_->bookHistograms(0.0,1.479,10,20, 20, 30);
  plotsNumPVOver20_Barrel_Pt20to30_->bookHistograms(0.0,1.479,20,50, 20, 30);
  plotsNumPV0to5_Endcap_Pt20to30_->bookHistograms(1.479,2.3,0,5, 20, 30);
  plotsNumPV5to10_Endcap_Pt20to30_->bookHistograms(1.479,2.3,5,10, 20, 30);
  plotsNumPV10to20_Endcap_Pt20to30_->bookHistograms(1.479,2.3,10,20, 20, 30);
  plotsNumPVOver20_Endcap_Pt20to30_->bookHistograms(1.479,2.3,20,50, 20, 30);
  plotsNumPV0to5_Barrel_Pt30to40_->bookHistograms(0.0,1.479,0,5, 30, 40);
  plotsNumPV5to10_Barrel_Pt30to40_->bookHistograms(0.0,1.479,5,10, 30, 40);
  plotsNumPV10to20_Barrel_Pt30to40_->bookHistograms(0.0,1.479,10,20, 30, 40);
  plotsNumPVOver20_Barrel_Pt30to40_->bookHistograms(0.0,1.479,20,50, 30, 40);
  plotsNumPV0to5_Endcap_Pt30to40_->bookHistograms(1.479,2.3,0,5, 30, 40);
  plotsNumPV5to10_Endcap_Pt30to40_->bookHistograms(1.479,2.3,5,10, 30, 40);
  plotsNumPV10to20_Endcap_Pt30to40_->bookHistograms(1.479,2.3,10,20, 30, 40);
  plotsNumPVOver20_Endcap_Pt30to40_->bookHistograms(1.479,2.3,20,50, 30, 40);
  plotsNumPV0to5_Barrel_Pt40to50_->bookHistograms(0.0,1.479,0,5, 40, 50);
  plotsNumPV5to10_Barrel_Pt40to50_->bookHistograms(0.0,1.479,5,10, 40, 50);
  plotsNumPV10to20_Barrel_Pt40to50_->bookHistograms(0.0,1.479,10,20, 40, 50);
  plotsNumPVOver20_Barrel_Pt40to50_->bookHistograms(0.0,1.479,20,50, 40, 50);
  plotsNumPV0to5_Endcap_Pt40to50_->bookHistograms(1.479,2.3,0,5, 40, 50);
  plotsNumPV5to10_Endcap_Pt40to50_->bookHistograms(1.479,2.3,5,10, 40, 50);
  plotsNumPV10to20_Endcap_Pt40to50_->bookHistograms(1.479,2.3,10,20, 40, 50);
  plotsNumPVOver20_Endcap_Pt40to50_->bookHistograms(1.479,2.3,20,50, 40, 50);
}


void AntiEMVAVariablesAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
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

  edm::Handle<CandidateView> genJets;
  evt.getByLabel(srcGenJets_, genJets);
  const reco::CandidateView* jets = genJets.product();
  //std::cout<<"number of genjets :"<<jets->size()<<std::endl;
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


  plotsAll_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsAll_Barrel_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsAll_Endcap_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsNumPV0to5_Barrel_Pt20to30_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsNumPV5to10_Barrel_Pt20to30_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsNumPV10to20_Barrel_Pt20to30_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);   
  plotsNumPVOver20_Barrel_Pt20to30_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsNumPV0to5_Endcap_Pt20to30_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsNumPV5to10_Endcap_Pt20to30_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsNumPV10to20_Endcap_Pt20to30_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);   
  plotsNumPVOver20_Endcap_Pt20to30_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsNumPV0to5_Barrel_Pt30to40_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsNumPV5to10_Barrel_Pt30to40_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsNumPV10to20_Barrel_Pt30to40_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);   
  plotsNumPVOver20_Barrel_Pt30to40_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsNumPV0to5_Endcap_Pt30to40_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsNumPV5to10_Endcap_Pt30to40_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsNumPV10to20_Endcap_Pt30to40_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);   
  plotsNumPVOver20_Endcap_Pt30to40_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsNumPV0to5_Barrel_Pt40to50_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsNumPV5to10_Barrel_Pt40to50_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsNumPV10to20_Barrel_Pt40to50_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);   
  plotsNumPVOver20_Barrel_Pt40to50_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsNumPV0to5_Endcap_Pt40to50_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsNumPV5to10_Endcap_Pt40to50_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);
  plotsNumPV10to20_Endcap_Pt40to50_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);   
  plotsNumPVOver20_Endcap_Pt40to50_->fillHistograms(*gsfElectrons, *genElectrons, *pfTaus, *genTaus, *genJets, debug_, numPV_, match_);

}

void AntiEMVAVariablesAnalyzer::endJob()
{
// nothing to be done yet...
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(AntiEMVAVariablesAnalyzer);






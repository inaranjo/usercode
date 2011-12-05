
/**
copy from CMSSW/PhysicsTools/Utilities/src/LumiReWeighting.cc
without any CMSSW dependancies
*/

#include "TH1.h"
#include "TFile.h"
#include <string>

float * LumiReWeighting( std::string generatedFile,
		 std::string dataFile,
		 std::string GenHistName = "pileup",
		 std::string DataHistName = "pileup" )
{
  std::string generatedFileName_ = generatedFile ;
  std::string dataFileName_ = dataFile ;
  std::string GenHistName_ = GenHistName ;
  std::string DataHistName_ = DataHistName ;
  TFile *     generatedFile_;
  TFile *     dataFile_;
  TH1F *      weights_;
  
  
  generatedFile_ = new TFile(generatedFileName_.c_str()); //MC distribution
  dataFile_      = new TFile(dataFileName_.c_str());      //Data distribution
  
  weights_ = new TH1F( *(static_cast<TH1F*>(dataFile_->Get( DataHistName_.c_str() )->Clone() )));
  
  // MC * data/MC = data, so the weights are data/MC:
  
  // normalize both histograms first
  
  weights_->Scale( 1.0/ weights_->Integral() );
  weights_->SetName("lumiWeights");
  
  TH1F* den = dynamic_cast<TH1F*>(generatedFile_->Get( GenHistName_.c_str() ));
  
  den->Scale(1.0/ den->Integral());
  
  weights_->Divide( den );  // so now the average weight should be 1.0
  
  std::cout << " Lumi/Pileup Reweighting: Computed Weights per In-Time Nint " << std::endl;
  
  int NBins = weights_->GetNbinsX();
  float * weightArray = new float[NBins] ;  

  for(int ibin = 1; ibin<NBins+1; ++ibin){
    std::cout << "   " << ibin-1 << " " << weights_->GetBinContent(ibin) << std::endl;
    weightArray[ibin-1] = weights_->GetBinContent(ibin) ;
  }

  generatedFile_->Close() ;
  dataFile_->Close() ;

  return weightArray ;

}

float * LumiReWeighting()
{
	
float Summer11Lumi_f[35] = {
   1.45346E-01,
   6.42802E-02,
   6.95255E-02,
   6.96747E-02,
   6.92955E-02,
   6.84997E-02,
   6.69528E-02,
   6.45515E-02,
   6.09865E-02,
   5.63323E-02,
   5.07322E-02,
   4.44681E-02,
   3.79205E-02,
   3.15131E-02,
   2.54220E-02,
   2.00184E-02,
   1.53776E-02,
   1.15387E-02,
   8.47608E-03,
   6.08715E-03,
   4.28255E-03,
   2.97185E-03,
   2.01918E-03,
   1.34490E-03,
   8.81587E-04,
   5.69954E-04,
   3.61493E-04,
   2.28692E-04,
   1.40791E-04,
   8.44606E-05,
   5.10204E-05,
   3.07802E-05,
   1.81401E-05,
   1.00201E-05,
   5.80004E-06
 };

 
float Data2011Lumi_f[35] = {
   0.00371368,
   0.0161937, 
   0.0382116, 
   0.0644434,
   0.0872967,
   0.101477, 
   0.105734, 
   0.101861, 
   0.092792, 
   0.0811974,
   0.0689395,
   0.0571016,
   0.0462355,
   0.0365947,
   0.0282805,
   0.0213088,
   0.0156338,
   0.0111578,
   0.00774186,
   0.00522119, 
   0.00342287, 
   0.00218201, 
   0.00135334, 
   0.000817217,
   0.000480835, 
   0.000275906, 
   0.000154535,
   8.45672e-05,
   4.52579e-05,
   2.37086e-05,
   1.21683e-05,
   6.12406e-06,
   3.02473e-06,
   1.46722e-06,
   1.30168e-06
 };
   float * weightArray = new float[35] ;
  for (int i=0 ; i<35 ; i++) weightArray[i] = Data2011Lumi_f[i]/Summer11Lumi_f[i] ;
  return weightArray ;
}


// double weight(TH1F * weights_, int npv ) {
//   std::cout<<npv<<std::endl;
//   int bin = weights_->GetXaxis()->FindBin( npv );
//   std::cout<<"bin "<<bin<<std::endl;
//   return weights_->GetBinContent( bin );
// }



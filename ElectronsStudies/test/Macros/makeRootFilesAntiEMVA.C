#include "TTree.h"
#include "TBranch.h"
#include "TString.h"
#include "TPad.h"
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include "TCut.h"

#include <string>
#include <map>
#include <iostream>
#include <iomanip>

#define DEBUG false



void makeRoot(string variable = "Elec_Fbrem",
	      const TString& category = "TauNoGammas",
	      const TString& xAxisTitle = "Fbrem",
	      const TString& yAxisTitle = "a.u.",
	      float xMin = -0.2, 
	      float xMax = 1,
	      int nBins = 100, 
	      int numPVMin = 0, 
	      int numPVMax = 50,
	      float PtMin = 10, 
	      float PtMax = 60,
	      const TString& Region = "Barrel"
	      )
{


  std::string inputFileName = "/data_CMS/cms/ivo/AntiEMVA/Trees/AntiEMVA_DYJetsToLL.root";
  TFile* inputFile = new TFile (inputFileName.data(),"READ");
  if(inputFile->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
//   TTree* inputTree = (TTree*)inputFile->Get("nTupleTree");
  TTree* inputTree = (TTree*)inputFile->Get("AntiEMVAAnalyzer/tree");

  std::string outputFileName = "./root/nTuple_AntiEMVA_TauNoGammas.root";
  TFile* outputFile = new TFile (outputFileName.data(),"RECREATE");



}



void makeAll(){

}

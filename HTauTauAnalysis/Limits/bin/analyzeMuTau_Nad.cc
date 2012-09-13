#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>

#include "TChain.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TH1F.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCut.h"
#include "TArrayF.h"
#include "TStyle.h"

#include "HiggsAnalysis/CombinedLimit/interface/TH1Keys.h"

#define VERBOSE          true
#define SAVE             true
#define addVH            true
#define EMBEDDEDSAMPLES  true
#define W3JETS           true
#define LOOSEISO         true
#define kFactorSM         1.0


///////////////////////////////////////////////////////////////////////////////////////////////

void makeHistoFromDensity(TH1* hDensity, TH1* hHistogram){

  if(hDensity->GetNbinsX() != hHistogram->GetNbinsX()){
    cout << "makeHistoFromDensity: different binning" << endl;
    return;
  }

  for(int k = 1 ; k <= hDensity->GetNbinsX(); k++){
    float bink   = hDensity->GetBinContent(k);
    float widthk = hHistogram->GetBinWidth(k);
    hDensity->SetBinContent(k, bink*widthk );
  }
  hDensity->Scale(hHistogram->Integral()/hDensity->Integral());
}

///////////////////////////////////////////////////////////////////////////////////////////////

void plotMuTau( TString mH_         = "120",
		Int_t   useEmbedding_ = 0,
	        TString selection_  = "inclusiveHigh",
		TString analysis_   = "",		  
		TString variable_   = "diTauVisMass",
		TString XTitle_     = "full mass",
		TString Unities_    = "GeV",
		TString dirOut      = "/home/llr/cms/ndaci/WorkArea/HTauTau/AnalysisLorenzo/CMSSW_5_2_6/src/Bianchi/Limits/bin/nad/results/",
		TString dirIn       = "/data_CMS/cms/ndaci/ndaci_2012/HTauTau/Analysis/MuTau/ntuples/", 
		Int_t   nBins_ = 80, Float_t xMin_=0, Float_t xMax_=400,
		Float_t magnifySgn_ = 1.0,
		Float_t hltEff_     = 1.0,
		Int_t   logy_       = 0,
		Float_t maxY_       = 1.2,
		Float_t Lumi        = 5.072663, // fb-1
		TString tsLumi      = "5.07",
		Float_t lumiCorrFactor = 0.976 // correct DY th cross-section
		) 
{   

  const int nMasses = 9;
  TString hMasses[nMasses] = {"110","115","120","125","130","135","140","145","155"} ;

  //TString postfix_ = "Raw";
  TString postfix_ = "";

  ofstream out(dirOut+"/histograms/yieldsMuTau_mH"+mH_+"_"+selection_+"_"+analysis_+".txt" , ios_base::out); 
  out.precision(5);
  out << " => " << selection_ << endl;
  
  // input txt file with bins
  ifstream is;

  char* c = new char[10];
  is.open(dirOut+"/bins/bins_eTau_"+variable_+"_"+selection_+".txt"); 
  if(nBins_<0 &&  !is.good()){
    cout << "Bins file not found" << endl;
    return;
  }

  int nBinsFromFile = 0;
  while (is.good())     
    {
      is.getline(c,999,',');     
      if (is.good()){
	nBinsFromFile++;
	//cout << c << endl;
      }
    }

  // choose the number of bins
  int nBins =  nBins_>0 ? nBins_ : nBinsFromFile-1 ;
  TArrayF bins(nBins+1);

  cout << "Making histograms with " << nBins << " bins:" << endl;

  is.close();
  is.open(dirOut+"/bins/bins_eTau_"+variable_+"_"+selection_+".txt"); 
  
  nBinsFromFile = 0;

  if(nBins_>0){
    for( ; nBinsFromFile <= nBins ; nBinsFromFile++){
      bins[nBinsFromFile] =  xMin_ + nBinsFromFile*(xMax_-xMin_)/nBins_;
    }
  }
  else{
    while (is.good())  
      {
	is.getline(c,999,',');     
	if (is.good() && nBinsFromFile<=nBins) {
	  bins[nBinsFromFile] = atof(c);
	  cout << bins[nBinsFromFile] << ", " ;
	}
	nBinsFromFile++;
      }
    cout << endl;
  }

  int nInArray = bins.GetSize();
  Float_t* array = bins.GetArray();

  cout << "-------------------" << endl << "ARRAY BINNING" << endl;
  cout << "From file : " << dirOut+"/bins/bins_eTau_"+variable_+"_"+selection_+".txt" << endl << endl;
  for(int iA=0 ; iA<nInArray ; iA++)
    cout << array[iA] << "   " ;
  cout << endl << "-------------------" << endl;

  //////////////////////////////////////////////////////////////////////////

  //float TTxsectionRatio           = lumiCorrFactor*(165.8/157.5) ;
  float TTxsectionRatio           = 1.00 ;

  float OStoSSRatioQCD            = 1.07;
  float SSIsoToSSAIsoRatioQCD     = 1.00;

  float MutoTauCorrectionFactor   = 1.00;
  float JtoTauCorrectionFactor    = 1.00;

  float embeddedMEtCutEff         = 1.00;
  float madgraphMEtCutEff         = 1.00;

  // Fall11_06Dec2011
  float WcorrectionFactorOS        = 0.92;  
  float WcorrectionFactorSS        = 1.08; 
  float ExtrapolationFactorZ       = 1.0;
  float ExtrapolationFactorZDataMC = 1.0;
  float ErrorExtrapolationFactorZ  = 1.0;

  //float NoVbfExtrapolationFactorZ = 0.997;
  //float VbfExtrapolationFactorZ   = 1.37;
  //float BoostExtrapolationFactorZ = 0.98;

  float VbfExtrapolationFactorW   = 1.00;
  float BoostExtrapolationFactorW = 1.00;
  float scaleFactorTTinVBF        = 1.00;
  float scaleFactorTTOS           = 1.0;
  /////////////////  change SVfit mass here ///////////////////

  //string variableStr = "";
  //TString variable(variableStr.c_str());
  TString variable = variable_;

  //////////////////////////////////////////////////////////////

  bool useMt      = true;
  TString antiWcut = useMt ? "MtLeg1Corr" : "-(pZetaCorr-1.5*pZetaVisCorr)" ; 
  TString antiWsgn  = useMt ? "40." :  "20." ;
  TString antiWsdb  = useMt ? "70." :  "40." ; 

  bool use2Dcut   = false;
  if( use2Dcut ){
    antiWcut = "!(MtLeg1Corr<40 && (pZetaCorr-1.5*pZetaVisCorr)>-20)";
    antiWsgn = "0.5";
    antiWsdb = "0.5";
  }

  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(logy_);

  TPad* pad1 = new TPad("pad1DEta","",0.05,0.22,0.96,0.97);
  TPad* pad2 = new TPad("pad2DEta","",0.05,0.02,0.96,0.20);
 
  pad1->SetFillColor(0);
  pad2->SetFillColor(0);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  pad1->SetLogy(logy_);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");

  TLegend* leg = new TLegend(0.63,0.48,0.85,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  leg->SetHeader( "#splitline{CMS Preliminary #sqrt{s}=8 TeV}{"+tsLumi+" fb^{-1} #tau_{#mu}#tau_{had}}" );

  THStack* aStack = new THStack("aStack","");

  TH1F* hSiml    = new TH1F( "hSiml"   ,"all"               , nBins , bins.GetArray());
  TH1F* hSgn     = new TH1F( "hSgn "   ,"vbf+ggf"           , nBins , bins.GetArray());
  TH1F* hSgn1    = new TH1F( "hSgn1"   ,"vbf"               , nBins , bins.GetArray());
  TH1F* hSgn2    = new TH1F( "hSgn2"   ,"ggf"               , nBins , bins.GetArray());
  TH1F* hSgn3    = new TH1F( "hSgn3"   ,"vh"                , nBins , bins.GetArray());
  TH1F* hData    = new TH1F( "hData"   ,"        "          , nBins , bins.GetArray());
  TH1F* hDataEmb = new TH1F( "hDataEmb","Embedded"          , nBins , bins.GetArray());
  TH1F* hW       = new TH1F( "hW"      ,"W+jets"            , nBins , bins.GetArray());
  TH1F* hWLooseIso2 = new TH1F( "hWLooseIso2","W+jets(miso)", nBins , bins.GetArray());
  TH1F* hW3Jets  = new TH1F( "hW3Jets" ,"W+3jets"           , nBins , bins.GetArray());
  TH1F* hEWK     = new TH1F( "hEWK"    ,"EWK"               , nBins , bins.GetArray());
  TH1F* hZtt     = new TH1F( "hZtt"    ,"Ztautau"           , nBins , bins.GetArray());
  TH1F* hZmm     = new TH1F( "hZmm"    ,"Z+jets, mu to tau" , nBins , bins.GetArray());
  TH1F* hZmmLoose= new TH1F( "hZmmLoose","Z+jets, mu to tau (loose)", nBins , bins.GetArray());
  TH1F* hZmj     = new TH1F( "hZmj"    ,"Z+jets, jet to tau", nBins , bins.GetArray());
  TH1F* hZmjLoose= new TH1F( "hZmjLoose","Z+jets, jet to tau (loose)",nBins , bins.GetArray());
  TH1F* hZfakes  = new TH1F( "hZfakes" ,"Z+jets, jet to tau", nBins , bins.GetArray());
  TH1F* hTTb     = new TH1F( "hTTb"    ,"ttbar"             , nBins , bins.GetArray());
  TH1F* hQCD     = new TH1F( "hQCD"    ,"QCD"               , nBins , bins.GetArray());
  TH1F* hLooseIso= new TH1F( "hLooseIso","Loose Iso"        , nBins , bins.GetArray());
  TH1F* hLooseIso2= new TH1F( "hLooseIso2","Loose Iso"      , nBins , bins.GetArray());
  TH1F* hLooseIso3= new TH1F( "hLooseIso3","Loose Iso"      , nBins , bins.GetArray());

  TH1F* hAntiIso = new TH1F( "hAntiIso","Anti Iso"          , nBins , bins.GetArray());
  TH1F* hVV      = new TH1F( "hVV"     ,"Diboson"           , nBins , bins.GetArray());

  ///////////////////////// for bkg estimation //////////////////////////

  TH1F* hW3JetsLooseTauIso         = new TH1F( "hW3JetsLooseTauIso" ,  "W+3jets (loose tau-iso)" ,                      nBins , bins.GetArray());
  TH1F* hW3JetsLooseTauIsoFR       = new TH1F( "hW3JetsLooseTauIsoFR" ,"W+3jets (loose tau-iso X fake-rate)"          , nBins , bins.GetArray());

  TH1F* hDataAntiIsoLooseTauIso      = new TH1F( "hDataAntiIsoLooseTauIso"   ,"data anti-iso, loose tau-iso"            , nBins , bins.GetArray());
  TH1F* hDataAntiIsoLooseTauIsoFR    = new TH1F( "hDataAntiIsoLooseTauIsoFR" ,"data anti-iso, loose tau-iso X fake-rate", nBins , bins.GetArray());
  TH1F* hDataAntiIsoLooseTauIsoFRUp  = new TH1F( "hDataAntiIsoLooseTauIsoFRUp"   ,"data anti-iso, loose tau-iso X fake-rate (Up)", nBins , bins.GetArray());
  TH1F* hDataAntiIsoLooseTauIsoFRDown= new TH1F( "hDataAntiIsoLooseTauIsoFRDown" ,"data anti-iso, loose tau-iso X fake-rate (Down)", nBins , bins.GetArray());

  TH1F* hDataAntiIso               = new TH1F( "hDataAntiIso"   ,"data anti-iso"                                 , nBins , bins.GetArray());
  TH1F* hDataAntiIsoFR             = new TH1F( "hDataAntiIsoFR" ,"data anti-iso X fake-rate"                     , nBins , bins.GetArray());
  TH1F* hDataAntiIsoFRUp           = new TH1F( "hDataAntiIsoFRUp" , "data anti-iso X fake-rate (Up)"             , nBins , bins.GetArray());
  TH1F* hDataAntiIsoFRDown         = new TH1F( "hDataAntiIsoFRDown" ,"data anti-iso X fake-rate (Down)"          , nBins , bins.GetArray());

  // get the FR-file

  TFile FakeRate(dirOut+"/fakeRate/FakeRate.root","READ");
  if(FakeRate.IsZombie()){
    cout << "Missing FR histos... exit" << endl;
    return;
  }
 
  TF1*  frMu     = (TF1*)FakeRate.Get("fit_MuTau_Mu_ptL1_incl");
  TH1F* frMuUp   = (TH1F*)FakeRate.Get("hFakeRateErrUpMuTau_Mu_ptL1_incl");
  TH1F* frMuDown = (TH1F*)FakeRate.Get("hFakeRateErrDownMuTau_Mu_ptL1_incl");

  TF1*  frTauQCD     = (TF1*)FakeRate.Get("fitQCD_MuTau_Tau_ptL2_QCDSS02_WSS60_incl");
  TH1F* frTauQCDUp   = (TH1F*)FakeRate.Get("hFakeRateQCDErrUpMuTau_Tau_ptL2_QCDSS02_WSS60_incl");
  TH1F* frTauQCDDown = (TH1F*)FakeRate.Get("hFakeRateQCDErrDownMuTau_Tau_ptL2_QCDSS02_WSS60_incl");
  
  TF1*  frTauW     = (TF1*)FakeRate.Get("fitW_MuTau_Tau_ptL2_QCDOS02_WOS60_incl");
  TH1F* frTauWUp   = (TH1F*)FakeRate.Get("hFakeRateWErrUpMuTau_Tau_ptL2_QCDOS02_WOS60_incl");
  TH1F* frTauWDown = (TH1F*)FakeRate.Get("hFakeRateWErrDownMuTau_Tau_ptL2_QCDOS02_WOS60_incl");

  if(!frMu || !frMuUp || !frMuDown || !frTauQCD || !frTauQCDUp || !frTauQCDDown || !frTauW || !frTauWUp || !frTauWDown){
    cout << "Missing FR histos... exit" << endl;
    return;
  }

  vector<int> binsFR;
  binsFR.push_back(17);
  binsFR.push_back(20);
  binsFR.push_back(22);
  binsFR.push_back(24);
  binsFR.push_back(26);
  binsFR.push_back(28);
  binsFR.push_back(30);
  binsFR.push_back(32);
  binsFR.push_back(34);
  binsFR.push_back(36);
  binsFR.push_back(40);
  binsFR.push_back(45);
  binsFR.push_back(50);
  binsFR.push_back(60); 
  binsFR.push_back(80); 
  binsFR.push_back(100);
  binsFR.push_back(9999);

  string scaleFactMu = "( ";
  string scaleFactMuUp   = "( ";
  string scaleFactMuDown = "( ";

  string scaleFactTauQCD = "( ";
  string scaleFactTauQCDUp   = "( ";
  string scaleFactTauQCDDown = "( ";

  string scaleFactTauW = "( ";
  string scaleFactTauWUp   = "( ";
  string scaleFactTauWDown = "( ";

 for(unsigned int i = 0; i < binsFR.size()-1; i++){
    
    float min = binsFR[i];
    float max = binsFR[i+1];

    float bin = frMuUp->FindBin((max+min)/2.);
    if( bin == frMuUp->GetNbinsX() + 1) bin--;

    float weightBinMu_i     =  frMu->Eval( (max+min)/2.);
    float weightBinMu_iUp   =  frMuUp->GetBinContent( bin );
    float weightBinMu_iDown =  frMuDown->GetBinContent( bin );
    
    scaleFactMu     += string( Form("(ptL1>=%f && ptL1<%f)*%f", min , max, 1./weightBinMu_i ) );
    scaleFactMuUp   += string( Form("(ptL1>=%f && ptL1<%f)*%f", min , max, 1./weightBinMu_iUp   ) );
    scaleFactMuDown += string( Form("(ptL1>=%f && ptL1<%f)*%f", min , max, 1./weightBinMu_iDown ) );

    float weightBinTauQCD_i     =  frTauQCD->Eval( (max+min)/2.);
    float weightBinTauQCD_iUp   =  frTauQCDUp->GetBinContent( bin );
    float weightBinTauQCD_iDown =  frTauQCDDown->GetBinContent( bin );
    
    scaleFactTauQCD     += string( Form("(ptL2>=%f && ptL2<%f)*%f", min , max, weightBinTauQCD_i ) );
    scaleFactTauQCDUp   += string( Form("(ptL2>=%f && ptL2<%f)*%f", min , max, weightBinTauQCD_iUp   ) );
    scaleFactTauQCDDown += string( Form("(ptL2>=%f && ptL2<%f)*%f", min , max, weightBinTauQCD_iDown ) );

    float weightBinTauW_i     =  frTauW->Eval( (max+min)/2.);
    float weightBinTauW_iUp   =  frTauWUp->GetBinContent( bin );
    float weightBinTauW_iDown =  frTauWDown->GetBinContent( bin );
    
    scaleFactTauW     += string( Form("(ptL2>=%f && ptL2<%f)*%f", min , max, weightBinTauW_i ) );
    scaleFactTauWUp   += string( Form("(ptL2>=%f && ptL2<%f)*%f", min , max, weightBinTauW_iUp   ) );
    scaleFactTauWDown += string( Form("(ptL2>=%f && ptL2<%f)*%f", min , max, weightBinTauW_iDown ) );

    if(i < binsFR.size() - 2 ){
      scaleFactMu     += " + ";
      scaleFactMuUp   += " + ";
      scaleFactMuDown += " + ";
      
      scaleFactTauQCD     += " + ";
      scaleFactTauQCDUp   += " + ";
      scaleFactTauQCDDown += " + ";
      
      scaleFactTauW     += " + ";
      scaleFactTauWUp   += " + ";
      scaleFactTauWDown += " + ";
    }
 }
 
 scaleFactMu     += " )";
 scaleFactMuUp   += " )";
 scaleFactMuDown += " )";
 
 scaleFactTauQCD     += " )";
 scaleFactTauQCDUp   += " )";
 scaleFactTauQCDDown += " )";
 
 scaleFactTauW     += " )";
 scaleFactTauWUp   += " )";
 scaleFactTauWDown += " )";
 
 cout << scaleFactMu << endl;
 cout << scaleFactTauQCD << endl;
 cout << scaleFactTauQCDUp << endl;
 cout << scaleFactTauW << endl;

  ////////////////////////////////////////////////////////////////////////////////////////////

  TH1*  hW3JetsKeys   = 0;
  TH1*  hWKeys        = 0;
  TH1*  hLooseIsoKeys = 0;
  TH1*  hAntiIsoKeys  = 0;
  TH1*  hZmmKeys      = 0;
  TH1*  hZmjKeys      = 0;
  TH1*  hVVKeys       = 0;

  ///////////////////////////////////////////////////////////////////////////////////////////

  vector<TH1F*> hggH, hqqH, hVH;
  TH1F* hBlub;

  for(int im=0 ; im<nMasses ; im++) {
    hBlub = new TH1F( "hggH"+hMasses[im] , "ggH"+hMasses[im] , nBins , bins.GetArray() );
    hggH.push_back(hBlub);
    hBlub = 0;

    hBlub = new TH1F( "hqqH"+hMasses[im] , "qqH"+hMasses[im] , nBins , bins.GetArray() );
    hqqH.push_back(hBlub);
    hBlub = 0;
    
    hBlub = new TH1F( "hVH"+hMasses[im] , "VH"+hMasses[im] , nBins , bins.GetArray() );
    hVH.push_back(hBlub);
    hBlub = 0;
  }
  //delete hBlub;

  ///////////////////////////////////////////////////////////////////////////////////////////

  // pZeta OS, N pZ sideband OS, pZeta SS, N sideband SS, N QCD SS, OS/SS
  TH1F* hParameters   = new TH1F( "hParameters", "" ,26, 0, 26);

  // Open the files //
  TString filename = "ntuple_MuTau_" ;

  // DATA //
  TFile *fData         = new TFile(dirIn+filename+"DATA_Nominal.root", "READ");  
  TFile *fDataLooseIso = new TFile(dirIn+filename+"DATA_Nominal.root", "READ");  
  TFile *fDataEmbedded = new TFile(dirIn+filename+"DataEmb_Nominal.root", "READ");  
  //TFile *fDataQCD = new TFile(dirIn+filename+"DataQCD_Nominal.root", "READ");  

  // BACKGROUNDS //
  TFile *fBackgroundDY     = new TFile(dirIn+filename+"DYJ_LL_Nominal.root", "READ"); 
  TFile *fBackgroundWJets  = new TFile(dirIn+filename+"Wj_Ln_Nominal.root", "READ"); 
  TFile *fBackgroundW3Jets = new TFile(dirIn+filename+"Wj_Ln_Nominal.root", "READ"); 
  //TFile *fBackgroundW3Jets = new TFile(dirIn+filename+"W3j_Ln_Nominal.root", "READ"); <-- to update (when official W3J sample available)
  TFile *fBackgroundTTbar  = new TFile(dirIn+filename+"TTJ_Nominal.root", "READ"); 
  TFile *fBackgroundOthers = new TFile(dirIn+filename+"VV_Nominal.root", "READ"); 

  // SIGNAL //
  vector<TFile*> fSignalGGH_m, fSignalqqH_m, fSignalVH_m;
  TFile *fSignalVBF    = new TFile(dirIn+filename+"GG_HTT_M"+mH_+"_Nominal.root", "READ");  
  TFile *fSignalGGH    = new TFile(dirIn+filename+"VBF_HTT_M"+mH_+"_Nominal.root", "READ"); 
  TFile *fSignalVH     = new TFile(dirIn+filename+"WH_ZH_TTH_HTT_M"+mH_+"_Nominal.root", "READ");  

  TFile* currentFile;

  for(int im=0 ; im<nMasses ; im++) {

    currentFile = TFile::Open(dirIn+filename+"GG_HTT_M"+hMasses[im]+"_Nominal.root");
//     if(!currentFile)
//       cout << "Error : File "+dirIn+filename+"GG_HTT_M"+hMasses[im]+"_Nominal.root does not exist !" << endl;
      
    fSignalGGH_m.push_back( currentFile );
    currentFile = 0;

    currentFile = TFile::Open(dirIn+filename+"VBF_HTT_M"+hMasses[im]+"_Nominal.root");
//     if(!currentFile)
//       cout << "Error : File "+dirIn+filename+"VBF_HTT_M"+hMasses[im]+"_Nominal.root does not exist !" << endl;

    fSignalqqH_m.push_back( currentFile );
    currentFile = 0;

    currentFile = TFile::Open(dirIn+filename+"WH_ZH_TTH_HTT_M"+hMasses[im]+"_Nominal.root");
//     if(!currentFile)
//       cout << "Error : File "+dirIn+filename+"VH_HTT_M"+hMasses[im]+"_Nominal.root does not exist !" << endl;
    cout << "file = " << currentFile << endl;
    fSignalVH_m.push_back( currentFile );
    currentFile = 0;
  }
  //currentFile->Close();
  //delete currentFile ;
 
  // Get the trees //

  // choose the analysis: Nominal "", jet up/Down "JetUp/Down" , elec up/down "MuUp/Down" , tau up/down "TauUp/Down"
  TString tree         = "outTreePtOrd"+postfix_+analysis_ ;
  TString treeEmbedded = "outTreePtOrd"+postfix_ ;

  if(analysis_.Contains("TauUp")  ) 
    treeEmbedded = tree;
  if(analysis_.Contains("TauDown")) 
    treeEmbedded = tree;
  if(analysis_.Contains("MuUp")  ) 
    treeEmbedded = tree;
  if(analysis_.Contains("MuDown")) 
    treeEmbedded = tree;

  TTree *data         = fData                            ? (TTree*)fData->Get( "outTreePtOrd" + postfix_ )         : 0 ;
  TTree *dataLooseIso = fDataLooseIso && LOOSEISO        ? (TTree*)fDataLooseIso->Get( "outTreePtOrd" + postfix_ ) : 0 ;
  TTree *dataEmbedded = fDataEmbedded && EMBEDDEDSAMPLES ? (TTree*)fDataEmbedded->Get(treeEmbedded)                : 0 ;
  TTree *signalVBF    = fSignalVBF                       ? (TTree*)fSignalVBF->Get(tree)                           : 0 ;
  TTree *signalGGH    = fSignalGGH                       ? (TTree*)fSignalGGH->Get(tree)                           : 0 ;
  TTree *signalVH     = fSignalVH     && addVH           ? (TTree*)fSignalVH->Get(tree)                            : 0 ;

  if( !fBackgroundDY ) { cout << "fBackgroundDY not found. Exit." << endl; return; }
  TTree* backgroundDY = (TTree*)fBackgroundDY->Get(tree) ;
  if( !backgroundDY ) { cout << "backgroundDY ("+tree+")not found. Exit." << endl; return; }

  // split the DY->ll into l=e/mu and l=tau (MC level) ===> temporary, need fix !!!
  TFile* dummy1 = new TFile("dummy2.root","RECREATE");

  cout << "Now copying g/Z -> tau+ tau- " << endl;
  TTree *backgroundDYTauTau  = ((TTree*)fBackgroundDY->Get(tree))->CopyTree("abs(genDecay)==(23*15)");                 // g/Z -> tau+ tau-

  cout << "Now copying g/Z -> mu+mu- mu->tau" << endl;
  TTree *backgroundDYMutoTau = ((TTree*)fBackgroundDY->Get(tree))->CopyTree("abs(genDecay)!=(23*15) &&  leptFakeTau"); // g/Z -> mu+mu- mu->tau

  cout << "Now copying g/Z -> mu+mu- jet->tau" << endl;
  TTree *backgroundDYJtoTau  = ((TTree*)fBackgroundDY->Get(tree))->CopyTree("abs(genDecay)!=(23*15) && !leptFakeTau"); // g/Z -> mu+mu- jet->tau

  cout << "///////////////////// HERE //////////////////" << endl;

  cout << backgroundDYTauTau->GetEntries()  << " come from DY->tautau"         << endl;
  cout << backgroundDYMutoTau->GetEntries() << " come from DY->mumu, mu->tau"  << endl;
  cout << backgroundDYJtoTau->GetEntries()  << " come from DY->mumu, jet->tau" << endl;

  // Backgrounds
  TTree *backgroundTTbar     = fBackgroundTTbar ? (TTree*)fBackgroundTTbar->Get(tree) : 0 ;
  TTree *backgroundWJets     = fBackgroundWJets ? (TTree*)fBackgroundWJets->Get(tree) : 0 ;
  TTree *backgroundOthers    = fBackgroundOthers ?(TTree*)fBackgroundOthers->Get(tree): 0 ;
  TTree *backgroundW3Jets    = fBackgroundW3Jets && W3JETS ? (TTree*)fBackgroundW3Jets->Get(tree) : 0 ;

  cout << "-- Take care of Lord Signal :o" << endl;

  // Signal
  vector<TTree*> tSignalGGH_m, tSignalqqH_m, tSignalVH_m;
  TTree* curTree;

  for(int im=0 ; im<nMasses ; im++) {

    curTree = fSignalVH_m[im]  ? (TTree*)fSignalVH_m[im]->Get(tree)  : 0 ;
    tSignalVH_m.push_back(curTree);

    curTree = fSignalGGH_m[im] ? (TTree*)fSignalGGH_m[im]->Get(tree) : 0 ;
    tSignalGGH_m.push_back(curTree);

    curTree = fSignalqqH_m[im] ? (TTree*)fSignalqqH_m[im]->Get(tree) : 0 ;
    tSignalqqH_m.push_back(curTree);

  }
  
  cout << "-- Take care of Susy :)" << endl;

  // SUSY
  std::map<TString,TFile*> mapSUSYfiles;
  /*
  for(unsigned int i = 0; i < SUSYhistos.size() ; i++){
    mapSUSYfiles.insert( make_pair(SUSYhistos[i], new TFile(Form("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_PreApproval//nTuple%s-MuTau-powheg-PUS6_run_Open_MuTauStream.root",SUSYhistos[i].c_str()) ,"READ")  )  );
  }
  */

  vector<TString> SUSYhistos;
  std::map<TString,TH1F*> mapSUSYhistos;
  for(unsigned int i = 0; i < SUSYhistos.size() ; i++){
    mapSUSYhistos.insert( make_pair( SUSYhistos[i], 
				     new TH1F(SUSYhistos[i] , SUSYhistos[i] , nBins , bins.GetArray() ) ) 
			  );
  }

  std::map<TString,TTree*> mapSUSYtrees;
  for(unsigned int i = 0; i < SUSYhistos.size() ; i++){
    TTree* treeSusy = (mapSUSYfiles.find(SUSYhistos[i]))->second ? (TTree*)((mapSUSYfiles.find(SUSYhistos[i]))->second)->Get(tree) : 0;
    mapSUSYtrees.insert( make_pair( SUSYhistos[i], treeSusy )) ;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "<--- Define the cuts --->" << endl;

  //// Re-Weighting ///
  TCut MCweight("sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau") ;
  
  ///// LEPT PT ///////
  TCut lpt("ptL1>20 && isPFMuon && isTightMuon");
  TCut tpt("ptL2>20");

  if(selection_.Contains("High"))
    tpt = tpt&&TCut("ptL2>40");
  else if(selection_.Contains("Low"))
    tpt = tpt&&TCut("ptL2<40");

  ////// TAU ISO //////
  TCut tiso("tightestHPSMVAWP>=0");  //<--------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  TCut ltiso("tightestHPSMVAWP>-99");
  TCut mtiso("hpsMVA>0.0");

  ////// MU ISO ///////
  TCut liso("combRelIsoLeg1DBetav2<0.10");
  TCut laiso("combRelIsoLeg1DBetav2>0.20 && combRelIsoLeg1DBetav2<0.50");
  TCut lliso("combRelIsoLeg1DBetav2<0.30");

 
  ////// EVENT WISE //////
  TCut lveto("muFlag==0"); //<----------------------------------
  TCut SS("diTauCharge!=0");
  TCut OS("diTauCharge==0");
  TCut pZ(   antiWcut + "<" + antiWsgn );
  TCut apZ(  antiWcut + ">" + antiWsdb );
  TCut apZ2( antiWcut + ">" + antiWsdb + " && " + antiWcut + "<120");
  TCut hltevent("pairIndex<1 && HLTx==1");
  TCut hltmatch("HLTmatch==1");


  ////// CATEGORIES ///
  TCut oneJet("nJets30>=1");
  TCut twoJets("nJets30>=2");
  //TCut vbf("pt1>30 && pt2>30 && eta1*eta2<0 && Mjj>400 && Deta>4.0 && isVetoInJets!=1");     // <--- BASELINE
  //TCut vbf("pt1>30 && pt2>30 && eta1*eta2<0 && Mjj>400 && Deta>4.0 && isVetoInJets!=1 && jet1PUWP>0.5 && jet2PUWP>0.5 && (jetVetoPUWP>0.5 && jetVetoPUWP<0)");
  TCut vbf("pt1>30 && pt2>30 && (ptVeto<30 || isVetoInJets!=1) && MVAvbf>0.50");
  TCut vbfLoose("pt1>30 && pt2>30 && (ptVeto<30 || isVetoInJets!=1 && MVAvbf>-0.30)"); /// <--- -0.30

  TCut vh("pt1>30 && pt2>30 && Mjj>70 && Mjj<120 && diJetPt>150 && MVAvbf<0.80 && nJets20BTagged<1");

  TCut boost("pt1>30 && nJets20BTagged<1");
  boost = boost && !vbf /*&& !vh*/;

  TCut boost2("pt1>100 && pt1<150 && !(pt2>30 && eta1*eta2<0 && Mjj>400 && Deta>4.0 && isVetoInJets!=1)");  

  TCut bTag("nJets30<2 && nJets20BTagged>0");
  TCut nobTag("nJets30<2 && nJets20BTagged==0");

  //TCut novbf = !vbf /*&& !vh*/ && !boost && !bTag;
  TCut novbf("nJets30<1 && nJets20BTagged<1");
  
  TCut sbin; TCut sbinEmbedding; TCut sbinEmbeddingPZetaRel; TCut sbinPZetaRel; TCut sbinSS; TCut sbinPZetaRelSS; TCut sbinPZetaRev; TCut sbinPZetaRevSS; TCut sbinSSaIso; TCut sbinSSlIso; TCut sbinSSlIso2; TCut sbinSSlIso3;
  TCut sbinSSaIsoLtiso;
  TCut sbinSSltiso;
  TCut sbinLtiso;
  TCut sbinMtiso;

  TCut sbinInclusive;
  sbinInclusive            = lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch;
  TCut sbinEmbeddingInclusive;
  sbinEmbeddingInclusive   = lpt && tpt && tiso && liso && lveto && OS && pZ                         ;
  TCut sbinPZetaRelEmbeddingInclusive;
  sbinPZetaRelEmbeddingInclusive = lpt && tpt && tiso && liso && lveto && OS                         ;
  TCut sbinPZetaRelSSInclusive;
  sbinPZetaRelSSInclusive  = lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch;
  TCut sbinPZetaRelInclusive;
  sbinPZetaRelInclusive    = lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch;
  TCut sbinSSInclusive;
  sbinSSInclusive          = lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch;
  TCut sbinSSaIsoInclusive;
  sbinSSaIsoInclusive      = lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch;
  TCut sbinPZetaRelSSaIsoInclusive;
  sbinPZetaRelSSaIsoInclusive = lpt && tpt && tiso && laiso&& lveto && SS && hltevent && hltmatch;
  TCut sbinSSaIsoLtisoInclusive;
  sbinSSaIsoLtisoInclusive =  lpt && tpt && ltiso&& laiso&& lveto && SS && pZ  && hltevent && hltmatch;
  TCut sbinSSltisoInclusive;
  sbinSSltisoInclusive     =  lpt && tpt && ltiso&& liso && lveto && SS && pZ  && hltevent && hltmatch;
  TCut sbinLtisoInclusive;
  sbinLtisoInclusive       =  lpt && tpt && ltiso&& liso && lveto && OS && pZ  && hltevent && hltmatch;
  TCut sbinMtisoInclusive;
  sbinMtisoInclusive       =  lpt && tpt && mtiso&& liso && lveto && OS && pZ  && hltevent && hltmatch;
  TCut sbinPZetaRelLtisoInclusive;
  sbinPZetaRelLtisoInclusive =lpt && tpt && ltiso&& liso && lveto && OS        && hltevent && hltmatch;

  if(selection_.Contains("inclusive")){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                         ;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                               ;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch;
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch;
    sbinSSlIso2            =  lpt && tpt && mtiso&& liso && lveto && SS && pZ  && hltevent && hltmatch;
    sbinSSlIso3            =  lpt && tpt && mtiso&& lliso&& lveto && SS && pZ  && hltevent && hltmatch;
    sbinSSaIsoLtiso        =  lpt && tpt && ltiso&& laiso&& lveto && SS && pZ  && hltevent && hltmatch;
    sbinSSltiso            =  lpt && tpt && ltiso&& liso && lveto && SS && pZ  && hltevent && hltmatch;
    sbinLtiso              =  lpt && tpt && ltiso&& liso && lveto && OS && pZ  && hltevent && hltmatch;
    sbinMtiso              =  lpt && tpt && mtiso&& liso && lveto && OS && pZ  && hltevent && hltmatch;
  }
  else if(selection_.Contains("oneJet")){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch && oneJet;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                          && oneJet;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                                && oneJet;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch && oneJet;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch && oneJet;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch && oneJet;
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch && oneJet;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch && oneJet;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch && oneJet;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch && oneJet;
    sbinSSlIso2            =  lpt && tpt && mtiso&& liso && lveto && SS && pZ  && hltevent && hltmatch && oneJet;
    sbinSSlIso3            =  lpt && tpt && mtiso&& lliso&& lveto && SS && pZ  && hltevent && hltmatch && oneJet;
    sbinSSaIsoLtiso        =  lpt && tpt && ltiso&& laiso&& lveto && SS && pZ  && hltevent && hltmatch && oneJet;
    sbinSSltiso            =  lpt && tpt && ltiso&& liso && lveto && SS && pZ  && hltevent && hltmatch && oneJet;
    sbinLtiso              =  lpt && tpt && ltiso&& liso && lveto && OS && pZ  && hltevent && hltmatch && oneJet;
    sbinMtiso              =  lpt && tpt && mtiso&& liso && lveto && OS && pZ  && hltevent && hltmatch && oneJet;
  }
  else if(selection_.Contains("twoJets")){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch && twoJets;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                          && twoJets;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                                && twoJets;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch && twoJets;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch && twoJets;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch && twoJets;
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch && twoJets;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch && twoJets;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch && twoJets;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch && twoJets;
    sbinSSlIso2            =  lpt && tpt && mtiso&& liso && lveto && SS && pZ  && hltevent && hltmatch && twoJets;
    sbinSSlIso3            =  lpt && tpt && mtiso&& lliso&& lveto && SS && pZ  && hltevent && hltmatch && twoJets;
    sbinSSaIsoLtiso        =  lpt && tpt && ltiso&& laiso&& lveto && SS && pZ  && hltevent && hltmatch && twoJets;
    sbinSSltiso            =  lpt && tpt && ltiso&& liso && lveto && SS && pZ  && hltevent && hltmatch && twoJets;
    sbinLtiso              =  lpt && tpt && ltiso&& liso && lveto && OS && pZ  && hltevent && hltmatch && twoJets;
    sbinMtiso              =  lpt && tpt && mtiso&& liso && lveto && OS && pZ  && hltevent && hltmatch && twoJets;
  }
  else if(selection_.Contains("vbf") ){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch && vbf;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                          && vbf;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                                && vbf;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch && vbf;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch && vbf;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch && vbf;
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch && vbf;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch && vbf;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch && vbf;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch && vbf;
    sbinSSlIso2            =  lpt && tpt && mtiso&& liso&& lveto && SS && pZ  && hltevent && hltmatch && vbf;
    sbinSSlIso3            =  lpt && tpt && mtiso&& lliso&& lveto && SS && pZ  && hltevent && hltmatch && vbf;
    sbinSSaIsoLtiso        =  lpt && tpt && ltiso&& laiso&& lveto && SS && pZ  && hltevent && hltmatch && vbf;
    sbinSSltiso            =  lpt && tpt && ltiso&& liso && lveto && SS && pZ  && hltevent && hltmatch && vbf;
    sbinLtiso              =  lpt && tpt && ltiso&& liso && lveto && OS && pZ  && hltevent && hltmatch && vbf;
    sbinMtiso              =  lpt && tpt && mtiso&& liso && lveto && OS && pZ  && hltevent && hltmatch && vbf;
  }
  else if(selection_.Contains("vh")){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch && vh;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                          && vh;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                                && vh;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch && vh;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch && vh;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch && vh;
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch && vh;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch && vh;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch && vh;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch && vh;
    sbinSSlIso2            =  lpt && tpt && mtiso&& liso && lveto && SS && pZ  && hltevent && hltmatch && vh;
    sbinSSlIso3            =  lpt && tpt && mtiso&& lliso&& lveto && SS && pZ  && hltevent && hltmatch && vh;
    sbinSSaIsoLtiso        =  lpt && tpt && ltiso&& laiso&& lveto && SS && pZ  && hltevent && hltmatch && vh;
    sbinSSltiso            =  lpt && tpt && ltiso&& liso && lveto && SS && pZ  && hltevent && hltmatch && vh;
    sbinLtiso              =  lpt && tpt && ltiso&& liso && lveto && OS && pZ  && hltevent && hltmatch && vh;
    sbinMtiso              =  lpt && tpt && mtiso&& liso && lveto && OS && pZ  && hltevent && hltmatch && vbf;
  }
  else if(selection_.Contains("0jet")){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch && novbf;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                          && novbf;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                                && novbf;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch && novbf;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch && novbf;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch && novbf;
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch && novbf;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch && novbf;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch && novbf;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch && novbf;
    sbinSSlIso2            =  lpt && tpt && mtiso&& liso && lveto && SS && pZ  && hltevent && hltmatch && novbf;
    sbinSSlIso3            =  lpt && tpt && mtiso&& lliso&& lveto && SS && pZ  && hltevent && hltmatch && novbf;
    sbinSSaIsoLtiso        =  lpt && tpt && ltiso&& laiso&& lveto && SS && pZ  && hltevent && hltmatch && novbf;
    sbinSSltiso            =  lpt && tpt && ltiso&& liso && lveto && SS && pZ  && hltevent && hltmatch && novbf;
    sbinLtiso              =  lpt && tpt && ltiso&& liso && lveto && OS && pZ  && hltevent && hltmatch && novbf;
    sbinMtiso              =  lpt && tpt && mtiso&& liso && lveto && OS && pZ  && hltevent && hltmatch && novbf;
  }
  else if(selection_.Contains("boost") && !selection_.Contains("boost2") ){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch && boost;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                          && boost;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                                && boost;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch && boost;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch && boost;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch && boost;
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch && boost;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch && boost;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch && boost;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch && boost;
    sbinSSlIso2            =  lpt && tpt && mtiso&& liso && lveto && SS && pZ  && hltevent && hltmatch && boost;
    sbinSSlIso3            =  lpt && tpt && mtiso&& lliso&& lveto && SS && pZ  && hltevent && hltmatch && boost;
    sbinSSaIsoLtiso        =  lpt && tpt && ltiso&& laiso&& lveto && SS && pZ  && hltevent && hltmatch && boost;
    sbinSSltiso            =  lpt && tpt && ltiso&& liso && lveto && SS && pZ  && hltevent && hltmatch && boost;
    sbinLtiso              =  lpt && tpt && ltiso&& liso && lveto && OS && pZ  && hltevent && hltmatch && boost;
    sbinMtiso              =  lpt && tpt && mtiso&& liso && lveto && OS && pZ  && hltevent && hltmatch && boost;
  }
  else if(selection_.Contains("boost2") ){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch && boost2;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                          && boost2;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                                && boost2;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch && boost2;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch && boost2;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch && boost2;
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch && boost2;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch && boost2;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch && boost2;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch && boost2;
    sbinSSlIso2            =  lpt && tpt && mtiso&& liso && lveto && SS && pZ  && hltevent && hltmatch && boost2;
    sbinSSlIso3            =  lpt && tpt && mtiso&& lliso&& lveto && SS && pZ  && hltevent && hltmatch && boost2;
    sbinSSaIsoLtiso        =  lpt && tpt && ltiso&& laiso&& lveto && SS && pZ  && hltevent && hltmatch && boost2;
    sbinSSltiso            =  lpt && tpt && ltiso&& liso && lveto && SS && pZ  && hltevent && hltmatch && boost2;
    sbinLtiso              =  lpt && tpt && ltiso&& liso && lveto && OS && pZ  && hltevent && hltmatch && boost2;
    sbinMtiso              =  lpt && tpt && mtiso&& liso && lveto && OS && pZ  && hltevent && hltmatch && boost2;
  }
  else if(selection_.Contains("bTag") && !selection_.Contains("nobTag")){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch && bTag;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                          && bTag;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                                && bTag;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch && bTag;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch && bTag;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch && bTag;    
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch && bTag;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch && bTag;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch && bTag;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch && bTag;
    sbinSSlIso2            =  lpt && tpt && mtiso&& liso && lveto && SS && pZ  && hltevent && hltmatch && bTag;
    sbinSSlIso3            =  lpt && tpt && mtiso&& lliso&& lveto && SS && pZ  && hltevent && hltmatch && bTag;
    sbinSSaIsoLtiso        =  lpt && tpt && ltiso&& laiso&& lveto && SS && pZ  && hltevent && hltmatch && bTag;
    sbinSSltiso            =  lpt && tpt && ltiso&& liso && lveto && SS && pZ  && hltevent && hltmatch && bTag;
    sbinLtiso              =  lpt && tpt && ltiso&& liso && lveto && OS && pZ  && hltevent && hltmatch && bTag;
    sbinMtiso              =  lpt && tpt && mtiso&& liso && lveto && OS && pZ  && hltevent && hltmatch && bTag;
  }
  else if(selection_.Contains("nobTag")){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch && nobTag;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                          && nobTag;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                                && nobTag;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch && nobTag;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch && nobTag;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch && nobTag;
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch && nobTag;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch && nobTag;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch && nobTag;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch && nobTag;
    sbinSSlIso2            =  lpt && tpt && mtiso&& liso && lveto && SS && pZ  && hltevent && hltmatch && nobTag;
    sbinSSlIso3            =  lpt && tpt && mtiso&& lliso&& lveto && SS && pZ  && hltevent && hltmatch && nobTag;
    sbinSSaIsoLtiso        =  lpt && tpt && ltiso&& laiso&& lveto && SS && pZ  && hltevent && hltmatch && nobTag;
    sbinSSltiso            =  lpt && tpt && ltiso&& liso && lveto && SS && pZ  && hltevent && hltmatch && nobTag;
    sbinLtiso              =  lpt && tpt && ltiso&& liso && lveto && OS && pZ  && hltevent && hltmatch && nobTag;
    sbinMtiso              =  lpt && tpt && mtiso&& liso && lveto && OS && pZ  && hltevent && hltmatch && nobTag;
  }


  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  cout << "******** Extrapolation factors for Z->tautau normalization ********" << endl;
  // inclusive DY->tautau:
  TH1F* hExtrap = new TH1F("hExtrap","",nBins , bins.GetArray());
  //TH1F* hExtrap = new TH1F("hExtrap","", 100, 0, 200);
  backgroundDYTauTau->Draw(variable+">>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbinInclusive);
  //backgroundDY->Draw(variable+">>hExtrap");
  cout << backgroundDY->GetEntries()  << " come from DY->tautau"         << endl;
  cout << variable << "   " << hExtrap->Integral() << endl;
  float ExtrapDYInclusive = hExtrap->Integral()*Lumi*lumiCorrFactor*hltEff_;
  hExtrap->Reset();
  cout << "All Z->tautau = " << ExtrapDYInclusive << endl; 

  backgroundDYMutoTau->Draw(variable+">>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbinInclusive);
  float ExtrapLFakeInclusive = hExtrap->Integral()*Lumi*lumiCorrFactor*hltEff_;
  hExtrap->Reset();
  cout << "All Z->mumu, mu->tau = " << ExtrapLFakeInclusive << endl;

  backgroundDYJtoTau->Draw(variable+">>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbinInclusive);
  float ExtrapJFakeInclusive = hExtrap->Integral()*Lumi*lumiCorrFactor*hltEff_;
  hExtrap->Reset();
  cout << "All Z->mumu, j->tau = " << ExtrapJFakeInclusive << endl;

  TCut sbinInclusiveEmbeddedCut = sbinEmbeddingInclusive;
  TCut sbinEmbeddedCut          = sbinEmbedding;

  // if VBF, minimize ttbar contamination asking for 0 btag jet:
  if(selection_.Contains("vbf")){
    sbinInclusiveEmbeddedCut = sbinInclusiveEmbeddedCut && TCut("nJets20BTagged<999");   /// <--------  1
    sbinEmbeddedCut          = sbinEmbeddedCut          && TCut("nJets20BTagged<999");
  }

  dataEmbedded->Draw(variable+">>hExtrap", "(HLTTau*HLTMu*1.000000)"*sbinInclusiveEmbeddedCut);
  float ExtrapEmbedDen =  hExtrap->Integral(); 
  hExtrap->Reset();
  dataEmbedded->Draw(variable+">>hExtrap", "(HLTTau*HLTMu*1.000000)"*sbinEmbeddedCut);
  float ExtrapEmbedNum =  hExtrap->Integral();
  hExtrap->Reset();

  ExtrapolationFactorZ = ExtrapEmbedNum/ExtrapEmbedDen;

  ErrorExtrapolationFactorZ = TMath::Sqrt(ExtrapolationFactorZ*(1-ExtrapolationFactorZ)/ExtrapEmbedDen);
  cout << "Extrap. factor using embedded sample: " << ExtrapolationFactorZ << " +/- " << ErrorExtrapolationFactorZ << endl;
  backgroundDYTauTau->Draw(variable+">>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbin);
  float ExtrapolationFactorMadGraph = hExtrap->Integral()*Lumi*lumiCorrFactor*hltEff_/ExtrapDYInclusive;
  cout << "MadGraph prediction = " << ExtrapolationFactorMadGraph << endl;
  ExtrapolationFactorZDataMC  = ExtrapolationFactorZ/ExtrapolationFactorMadGraph;
  cout << " ==> data/MC = " << ExtrapolationFactorZDataMC << endl;

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  cout << "******** Extrapolation factors for QCD normalization ********" << endl;
  hExtrap->Reset();
  backgroundWJets->Draw(variable+">>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSSInclusive&&pZ));
  float ExtrapSSWinSignalRegionMC   = hExtrap->Integral();
  hExtrap->Reset();
  backgroundWJets->Draw(variable+">>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSSInclusive&&apZ));
  float ExtrapSSWinSidebandRegionMC = hExtrap->Integral();
  float ExtrapscaleFactorSS         = ExtrapSSWinSignalRegionMC>0 ? ExtrapSSWinSidebandRegionMC/ExtrapSSWinSignalRegionMC : 1.0;
  cout << " Extrapolation factor W SS (inclusive) " << ExtrapscaleFactorSS << endl;

  hExtrap->Reset();
  backgroundTTbar->Draw(variable+">>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSSInclusive&&apZ));
  float ExtrapttbarExtrSS    = hExtrap->Integral()*Lumi*hltEff_*TTxsectionRatio;

  hExtrap->Reset();
  backgroundTTbar->Draw(variable+">>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelSSInclusive&&apZ&&TCut("nJets20BTagged>1")));
  float ExtrapttbarExtrSSBtag = hExtrap->Integral()*Lumi*hltEff_*TTxsectionRatio;
  hExtrap->Reset();
  backgroundWJets->Draw(variable+">>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelSSInclusive&&apZ&&TCut("nJets20BTagged>1")));
  float ExtrapwjetsExtrSSBtag = hExtrap->Integral()*Lumi*hltEff_;
  hExtrap->Reset();
  backgroundOthers->Draw(variable+">>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelSSInclusive&&apZ&&TCut("nJets20BTagged>1")));
  float ExtrapothersExtrSSBtag = hExtrap->Integral()*Lumi*hltEff_;
  hExtrap->Reset();
  data->Draw(variable+">>hExtrap",(sbinPZetaRelSSInclusive&&apZ&&TCut("nJets20BTagged>1")));
  float dataSSBtag = hExtrap->Integral();
  float scaleFactorTTSS = ExtrapttbarExtrSSBtag>0 ? dataSSBtag/ExtrapttbarExtrSSBtag : 1.0;
  //if( (dataSSBtag-ExtrapttbarExtrSSBtag)/TMath::Sqrt(dataSSBtag)>1) scaleFactorTTSS = 1.0;
  cout << "Normalizing TTbar from sideband: " << ExtrapttbarExtrSSBtag << " events expected from TTbar, " << ExtrapwjetsExtrSSBtag
       << " expected from WJets, "  << ", expected from others " << ExtrapothersExtrSSBtag << ", observed " << dataSSBtag << endl;
  cout << "====> scale factor for SS ttbar is " << scaleFactorTTSS << endl;
  ExtrapttbarExtrSS *= scaleFactorTTSS;

  hExtrap->Reset();
  backgroundOthers->Draw(variable+">>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSSInclusive&&apZ));
  float ExtrapothersExtrSS   = hExtrap->Integral()*Lumi*hltEff_;
  hExtrap->Reset();
  backgroundDYJtoTau->Draw(variable+">>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSSInclusive&&apZ));
  float ExtrapdyjtotauExtrSS = hExtrap->Integral()*Lumi*lumiCorrFactor*hltEff_;

  hExtrap->Reset();
  data->Draw(variable+">>hExtrap", sbinPZetaRelSSInclusive&&apZ);
  float ExtrapSSWinSignalRegionDATA = hExtrap->Integral();
  cout << "Extrapolation for QCD (inclusive): total data events in sideband " << ExtrapSSWinSignalRegionDATA << endl;
  ExtrapSSWinSignalRegionDATA -= ExtrapttbarExtrSS;
  ExtrapSSWinSignalRegionDATA -= ExtrapothersExtrSS;
  ExtrapSSWinSignalRegionDATA -= ExtrapdyjtotauExtrSS;
  ExtrapSSWinSignalRegionDATA /= ExtrapscaleFactorSS;
  cout << "Extrapolation for QCD (inclusive): W+jets in SS signal region (inclusive) is estimated to be " << ExtrapSSWinSignalRegionDATA << endl;

  hExtrap->Reset();
  data->Draw(variable+">>hExtrap", sbinSSInclusive);
  float SSeventsExtrap = hExtrap->Integral();
  cout << "Extrapolation for SS events in data (inclusive) " << hExtrap->GetEntries() << endl;
  cout << "Subtracting W+jets (SS)..." << endl;
  SSeventsExtrap  -= ExtrapSSWinSignalRegionDATA;

  hExtrap->Reset();
  backgroundTTbar->Draw(variable+">>hExtrap", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbinSSInclusive);
  SSeventsExtrap  -= hExtrap->Integral()*Lumi*hltEff_*TTxsectionRatio*scaleFactorTTSS;
  hExtrap->Reset();
  backgroundDYMutoTau->Draw(variable+">>hExtrap", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbinSSInclusive);
  SSeventsExtrap  -= hExtrap->Integral()*Lumi*lumiCorrFactor*hltEff_*MutoTauCorrectionFactor;

  hExtrap->Reset();
  backgroundDYJtoTau->Draw(variable+">>hExtrap", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbinSSInclusive);
  SSeventsExtrap  -= hExtrap->Integral()*Lumi*lumiCorrFactor*hltEff_*JtoTauCorrectionFactor;
  hExtrap->Reset();

  //SSeventsExtrap *= OStoSSRatioQCD;

  dataLooseIso->Draw(variable+">>hExtrap", sbinSSaIsoInclusive);
  float SSeventsExtrapAiso = hExtrap->GetEntries();
  SSIsoToSSAIsoRatioQCD = SSeventsExtrap/SSeventsExtrapAiso ;
  cout << "The extrapolation factor Iso>0.3 / Iso<0.1 is " << SSIsoToSSAIsoRatioQCD << endl;

  cout << "************** END extrapolation *******************" << endl;
  delete hExtrap;
  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////


  // estimate the W+jets in the selection bin using pZeta extrapolation

  //TH1F* hWMt = new TH1F("hWMt","",1,-10,10);
  TH1F* hWMt = new TH1F("hWMt","",nBins , bins.GetArray());

  ///////////////////////////////////////// Doing OS first...
  hWMt->Reset();
  backgroundWJets->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRel&&pZ));
  cout << "Using  " << hWMt->GetEntries() << " entries from the W+jets OS sample" << endl;
  float OSWinSignalRegionMC   = hWMt->Integral()*Lumi*hltEff_;
  hWMt->Reset();
  backgroundWJets->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRel&&apZ));
  float OSWinSidebandRegionMC = hWMt->Integral()*Lumi*hltEff_;
  float scaleFactorOS = OSWinSignalRegionMC>0 ? OSWinSidebandRegionMC/OSWinSignalRegionMC : 1.0 ;

  if(useMt)
    cout << "Extrapolation factor for W OS : P(MtCorr>" << antiWsdb << ")/P(MtCorr<" << antiWsgn << ") ==> " << scaleFactorOS << endl;
  else
    cout << "Extrapolation factor for W OS : P(pZetaCorr<- "<< antiWsdb << ")/P(pZetaCorr>"<< antiWsgn << ") ==> " << scaleFactorOS << endl;    
 
  hWMt->Reset();
  cout << "Estimating cobtribution from Ztt, ttbar and others in OS low pZeta tail..." << endl;
  backgroundTTbar->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRel&&apZ));
  float ttbarExtrOS  = hWMt->Integral()*Lumi*hltEff_*TTxsectionRatio;
  hWMt->Reset();

  if(selection_.Contains("boost"))
    backgroundTTbar->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelInclusive&&apZ&&TCut("nJets20BTagged>1")&&TCut("pt1>30")&&(!vbf)));
  else
    backgroundTTbar->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRel&&apZ&&TCut("nJets20BTagged>1")));
  float ttbarExtrOSBtag = hWMt->Integral()*Lumi*hltEff_*TTxsectionRatio;
  hWMt->Reset();

  if(selection_.Contains("boost"))
    backgroundWJets->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelInclusive&&apZ&&TCut("nJets20BTagged>1")&&TCut("pt1>30")&&(!vbf)));
  else
    backgroundWJets->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRel&&apZ&&TCut("nJets20BTagged>1")));
  float wjetsExtrOSBtag = hWMt->Integral()*Lumi*hltEff_;
  hWMt->Reset();

  if(selection_.Contains("boost"))
    backgroundOthers->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelInclusive&&apZ&&TCut("nJets20BTagged>1")&&TCut("pt1>30")&&(!vbf)));
  else
    backgroundOthers->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRel&&apZ&&TCut("nJets20BTagged>1")));
  float othersExtrOSBtag = hWMt->Integral()*Lumi*hltEff_;
  hWMt->Reset();

  if(selection_.Contains("boost"))
    data->Draw(variable+">>hWMt",(sbinPZetaRelInclusive&&apZ&&TCut("nJets20BTagged>1")&&TCut("pt1>30")&&(!vbf)));
  else
    data->Draw(variable+">>hWMt",(sbinPZetaRel&&apZ&&TCut("nJets20BTagged>1")));

  float dataOSBtag = hWMt->Integral();
  scaleFactorTTOS = ttbarExtrOSBtag>0 ? dataOSBtag/ttbarExtrOSBtag : 1.0;
  cout << "Normalizing TTbar from sideband: " << ttbarExtrOSBtag << " events expected from TTbar, " 
       << ", from WJets" << wjetsExtrOSBtag <<  ", expected from others " << othersExtrOSBtag << ", observed " << dataOSBtag << endl;
  cout << "====> scale factor for OS ttbar is " << scaleFactorTTOS << endl;
  //if( (dataOSBtag-ttbarExtrOSBtag)/TMath::Sqrt(dataOSBtag) > 1)
  ttbarExtrOS *= scaleFactorTTOS;
  cout << "Contribution from ttbar in OS is " << ttbarExtrOS << endl;


  hWMt->Reset();
  backgroundOthers->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRel&&apZ));
  float othersExtrOS = hWMt->Integral()*Lumi*hltEff_;
  cout << "Contribution from single-t and di-boson in OS is " << othersExtrOS << endl;
  hWMt->Reset();
  backgroundDYTauTau->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRel&&apZ));
  float dytautauExtrOS = hWMt->Integral()*Lumi*lumiCorrFactor*hltEff_;
  cout << "Contribution from DY->tautau in OS is " << dytautauExtrOS << endl;
  hWMt->Reset();
  backgroundDYJtoTau->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRel&&apZ));
  float dyjtotauExtrOS = hWMt->Integral()*Lumi*lumiCorrFactor*hltEff_;
  cout << "Contribution from DY->mumu, jet->tau in OS is " << dyjtotauExtrOS << endl;
  hWMt->Reset();
  backgroundDYMutoTau->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRel&&apZ));
  float dymutotauExtrOS = hWMt->Integral()*Lumi*lumiCorrFactor*hltEff_;
  cout << "Contribution from DY->mumu, mu->tau in OS is " << dymutotauExtrOS << endl;
  hWMt->Reset();

  data->Draw(variable+">>hWMt", sbinPZetaRev);
  float OSWinSignalRegionDATA = hWMt->Integral();
  cout << "Selected events in data in low pZeta/low Mt tail " << OSWinSignalRegionDATA << endl;
  OSWinSignalRegionDATA -= ttbarExtrOS;
  OSWinSignalRegionDATA -= othersExtrOS;
  OSWinSignalRegionDATA -= dytautauExtrOS;
  OSWinSignalRegionDATA -= dyjtotauExtrOS;
  OSWinSignalRegionDATA -= dymutotauExtrOS;
  OSWinSignalRegionDATA /= scaleFactorOS;
  cout << "W+jets in signal region is estimated to be "  
       << OSWinSignalRegionDATA*scaleFactorOS << "/" << scaleFactorOS << " = " 
       << OSWinSignalRegionDATA <<  " +/- " << sqrt(OSWinSignalRegionDATA/scaleFactorOS)/scaleFactorOS << endl;
  cout << "  ===> the MC prediction was " << OSWinSignalRegionMC << endl;

  hParameters->SetBinContent(1, 1./scaleFactorOS );
  hParameters->SetBinContent(2, OSWinSignalRegionDATA*scaleFactorOS );

  ///////////////////////////////////////// Doing SS last...
  hWMt->Reset();
  backgroundWJets->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSS&&pZ));
  cout << "Using  " << hWMt->GetEntries() << " entries from the SS W+jets sample" << endl;
  float SSWinSignalRegionMC   = hWMt->Integral()*Lumi*hltEff_;
  hWMt->Reset();
  backgroundWJets->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSS&&apZ));
  float SSWinSidebandRegionMC = hWMt->Integral()*Lumi*hltEff_;
  float scaleFactorSS = SSWinSignalRegionMC>0 ? SSWinSidebandRegionMC/SSWinSignalRegionMC : 1.0;
 
  if(useMt)
    cout << "Extrapolation factor for W SS : P(MtCorr>" << antiWsdb << ")/P(MtCorr<" << antiWsgn << ") ==> " << scaleFactorSS << endl;
  else
    cout << "Extrapolation factor for W SS : P(pZetaCorr<- "<< antiWsdb << ")/P(pZetaCorr>"<< antiWsgn << ") ==> " << scaleFactorSS << endl;    

  hWMt->Reset();
  cout << "Estimating cobtribution Ztt,from ttbar and others in SS low pZeta tail..." << endl;
  backgroundTTbar->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSS&&apZ));
  float ttbarExtrSS = hWMt->Integral()*Lumi*hltEff_*TTxsectionRatio;
  hWMt->Reset();

  if(selection_.Contains("boost"))
    backgroundTTbar->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelSSInclusive&&apZ&&TCut("nJets20BTagged>1")&&TCut("pt1>30")&&(!vbf)));
  else
    backgroundTTbar->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelSS&&apZ&&TCut("nJets20BTagged>1")));
  float ttbarExtrSSBtag = hWMt->Integral()*Lumi*hltEff_*TTxsectionRatio;
  hWMt->Reset();

  if(selection_.Contains("boost"))
    backgroundWJets->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelSSInclusive&&apZ&&TCut("nJets20BTagged>1")&&TCut("pt1>30")&&(!vbf)));
  else
    backgroundWJets->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelSS&&apZ&&TCut("nJets20BTagged>1")));
  float wjetsExtrSSBtag = hWMt->Integral()*Lumi*hltEff_;
  hWMt->Reset();

  if(selection_.Contains("boost"))
    backgroundOthers->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelSSInclusive&&apZ&&TCut("nJets20BTagged>1")&&TCut("pt1>30")&&(!vbf)));
  else
    backgroundOthers->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelSS&&apZ&&TCut("nJets20BTagged>1")));
  float othersExtrSSBtag = hWMt->Integral()*Lumi*hltEff_;
  hWMt->Reset();

  if(selection_.Contains("boost"))
    data->Draw(variable+">>hWMt",(sbinPZetaRelSSInclusive&&apZ&&TCut("nJets20BTagged>1")&&TCut("pt1>30")&&(!vbf)));
  else
    data->Draw(variable+">>hWMt",(sbinPZetaRelSS&&apZ&&TCut("nJets20BTagged>1")));
  float dataSSBtag2 = hWMt->Integral();
  float scaleFactorTTSS2 = ttbarExtrSSBtag>0 ? dataSSBtag2/ttbarExtrSSBtag : 1.0;
  cout << "Normalizing TTbar from sideband: " << ttbarExtrSSBtag << " events expected from TTbar, " 
       << ", from WJets " << wjetsExtrSSBtag << ", expected from others " << othersExtrSSBtag << ", observed " << dataSSBtag2 << endl;
  cout << "====> scale factor for SS ttbar is " << scaleFactorTTSS2 << endl;
  //if( (dataSSBtag2-ttbarExtrSSBtag)/TMath::Sqrt(dataSSBtag2) > 1) 
  ttbarExtrSS *= scaleFactorTTSS2;
  cout << "Contribution from ttbar in SS is " << ttbarExtrSS << endl;


  hWMt->Reset();
  backgroundOthers->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSS&&apZ));
  float othersExtrSS = hWMt->Integral()*Lumi*hltEff_;
  cout << "Contribution from single-t and di-boson in SS is " << othersExtrSS << endl;
  hWMt->Reset();
  backgroundDYJtoTau->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSS&&apZ));
  float dyjtotauExtrSS = hWMt->Integral()*Lumi*lumiCorrFactor*hltEff_;
  cout << "Contribution from DY->mumu, jet->tau in SS is " << dyjtotauExtrSS << endl;
  hWMt->Reset();

  data->Draw(variable+">>hWMt",sbinPZetaRevSS);
  float SSWinSignalRegionDATA = hWMt->Integral();
  cout << "Selected events in data in low pZeta/low Mt tail " << SSWinSignalRegionDATA << endl;
  SSWinSignalRegionDATA -= ttbarExtrSS;
  SSWinSignalRegionDATA -= othersExtrSS;
  SSWinSignalRegionDATA -= dyjtotauExtrSS;
  SSWinSignalRegionDATA /= scaleFactorSS;
  cout << "W+jets in SS signal region is estimated to be "  
       << SSWinSignalRegionDATA*scaleFactorSS << "/" << scaleFactorSS << " = " 
       << SSWinSignalRegionDATA <<  " +/- " << sqrt(SSWinSignalRegionDATA/scaleFactorSS)/scaleFactorSS << endl;
  cout << "  ===> the MC prediction was " << SSWinSignalRegionMC << endl;

  hParameters->SetBinContent(3, 1./scaleFactorSS );
  hParameters->SetBinContent(4, SSWinSignalRegionDATA*scaleFactorSS );

  // here I choose the order in the stack
  std::vector<TString> samples;
  samples.push_back("Data");
  if(dataLooseIso){
    samples.push_back("LooseIso");
    samples.push_back("AntiIso");
  }
  samples.push_back("Others");
  samples.push_back("SS");
  samples.push_back("WJets");
  if(backgroundW3Jets)
    samples.push_back("W3Jets");
  samples.push_back("TTbar");
  samples.push_back("DYMutoTau");
  samples.push_back("DYJtoTau");
  samples.push_back("DYToTauTau");
  if(dataEmbedded)
    samples.push_back("Embedded");

  for(int i = 0 ; i < nMasses ; i++) {
    samples.push_back("ggH" + hMasses[i]);
    samples.push_back("qqH" + hMasses[i]);
    samples.push_back("VH"  + hMasses[i]);
  }

  for(u_int i = 0; i < SUSYhistos.size() ; i++){
    TTree* susyTree = (mapSUSYtrees.find( SUSYhistos[i] ))->second ;
    if( susyTree ) samples.push_back( SUSYhistos[i] );
  }


  // here I define the map between a sample name and its tree
  std::map<TString,TTree*> tMap;
  tMap["Data"]         = data;
  tMap["LooseIso"]     = dataLooseIso;
  tMap["AntiIso"]      = dataLooseIso;
  tMap["Embedded"]     = dataEmbedded;
  tMap["DYToTauTau"]   = backgroundDYTauTau;
  tMap["DYMutoTau"]    = backgroundDYMutoTau;
  tMap["DYJtoTau"]     = backgroundDYJtoTau;
  tMap["WJets"]        = backgroundWJets;
  tMap["W3Jets"]       = backgroundW3Jets;
  tMap["Others"]       = backgroundOthers;
  tMap["TTbar"]        = backgroundTTbar;
  tMap["SS"]           = data;

  for(int im=0 ; im<nMasses ; im++) {
    
    tMap["ggH"+hMasses[im]] = tSignalGGH_m[im] ;
    tMap["qqH"+hMasses[im]] = tSignalqqH_m[im] ;
    tMap["VH"+hMasses[im]]  = tSignalVH_m[im] ;    
  }

  for(unsigned int i = 0; i < SUSYhistos.size() ; i++){
    TTree* susyTree = (mapSUSYtrees.find( SUSYhistos[i] ))->second ;
    tMap[SUSYhistos[i]] = susyTree ;
  }

  //for( unsigned iter=0; iter<samples.size(); iter++){
  for( unsigned iter=0; iter<=0; iter++){

    cout << "Dealing with sample " << samples[iter] ;

    std::map<TString,TTree*>::iterator it = tMap.find(samples[iter]);
    TTree* currentTree = (it->second);
    if( !currentTree ) {
      cout << "  --> not found in files, go to next sample" << endl;
      continue;
    }
    else cout << endl;

    TString h1Name = "h1_"+it->first;
    TH1F* h1 = new TH1F( h1Name ,"" , nBins , bins.GetArray());

    if((it->first).Contains("SS")){
      
      cout << "Remove W contamination from SS data sample ... " << endl;

      float error2OnQCD = 0.0;
      
      TH1F* hHelp = (TH1F*)h1->Clone("hHelp");
      hHelp->Reset();
      currentTree->Draw(variable+">>hHelp", sbinSS);
      int SSevents = hHelp->GetEntries();
      cout << "Selected SS events in data " << hHelp->GetEntries() << endl;
      h1->Add(hHelp,1);

      hHelp->Reset();
      backgroundWJets->Draw(variable+">>hHelp", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbinSS);
      cout << "We expect " << hHelp->Integral()*Lumi*hltEff_ << " SS events from W+jets (from " << hHelp->GetEntries() << " entries)" << endl;
      float sFWSS = ( selection_.Contains("0jet")  || 
		      selection_.Contains("bTag")   || 
		      selection_.Contains("boost")  || 
		      selection_.Contains("oneJet")  || 
		      selection_.Contains("twoJets")  || 
		      selection_.Contains("inclusive")) ? 
	SSWinSignalRegionDATA/SSWinSignalRegionMC : WcorrectionFactorSS; // from the extrapolation factor DATA/MC

      if(selection_.Contains("vbf")) 
	sFWSS *= VbfExtrapolationFactorW;
      else if(selection_.Contains("boost"))
	sFWSS *= BoostExtrapolationFactorW;

      hHelp->Scale(sFWSS*Lumi*hltEff_);
      cout << "We estimate " << hHelp->Integral() << " SS events from W+jets by extrapolating" << endl;
      cout << " ==> removing W+jets from SS...." << endl;
      h1->Add(hHelp, -1 );
      if(hHelp->GetEntries()>0) error2OnQCD +=  pow( hHelp->Integral()/hHelp->GetEntries(), 2)*hHelp->GetEntries(); // error on MC W+jets SS events
      error2OnQCD +=  pow(WcorrectionFactorSS*0.06,2)*pow(hHelp->GetEntries(),2);        // error on W+jets extrapolation factor ==> 6% according to Artur
      cout << sqrt(error2OnQCD) << " <==  W" << endl;      

      hHelp->Reset();
      backgroundTTbar->Draw(variable+">>hHelp", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbinSS);
      cout << "We expect " << hHelp->Integral()*Lumi*hltEff_*TTxsectionRatio << " SS events from TTbar (from " << hHelp->GetEntries() << " entries)" << endl;
      hHelp->Scale(1.0*Lumi*hltEff_*TTxsectionRatio);
      cout << "We estimate " << hHelp->Integral() << " SS events from TTbar" << endl;
      cout << " ==> removing TTbar from SS...." << endl;
      h1->Add(hHelp, -1 );
      if(hHelp->GetEntries()>0) error2OnQCD +=  pow(hHelp->Integral()/hHelp->GetEntries(),2)*hHelp->GetEntries();   // error on MC TTbar SS events
      cout << sqrt(error2OnQCD) << " <== W + TTb" << endl;      

      hHelp->Reset();
      backgroundDYMutoTau->Draw(variable+">>hHelp", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbinSS);
      cout << "We expect " << hHelp->Integral()*Lumi*hltEff_ << " SS events from DY->mumu, mu->jet" << endl;
      hHelp->Scale(MutoTauCorrectionFactor*Lumi*lumiCorrFactor*hltEff_);
      cout << "We estimate " << hHelp->Integral() << " SS events from DY->mumu, mu->tau" << endl;
      cout << " ==> removing DY->mumu, mu->tau from SS...." << endl;
      h1->Add(hHelp, -1 );
      if(hHelp->GetEntries()>0) error2OnQCD +=  pow(hHelp->Integral()/hHelp->GetEntries(),2)*hHelp->GetEntries();   // error on MC DY->mumu, mu->tau events
      cout << sqrt(error2OnQCD) << " <== W + TTb + DY(1)" << endl;      

      hHelp->Reset();
      backgroundDYJtoTau->Draw(variable+">>hHelp", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbinSS);
      cout << "We expect " << hHelp->Integral()*Lumi*lumiCorrFactor*hltEff_ << " SS events from DY->mumu, jet->tau" << endl;
      hHelp->Scale(JtoTauCorrectionFactor*Lumi*lumiCorrFactor*hltEff_);
      cout << "We estimate " << hHelp->Integral() << " SS events from DY->mumu, jet->tau" << endl;
      cout << " ==> removing DY->mumu, mu->jet from SS...." << endl;
      h1->Add(hHelp, -1 );
      if(hHelp->GetEntries()>0) error2OnQCD +=  pow(hHelp->Integral()/hHelp->GetEntries(),2)*hHelp->GetEntries();   // error on MC DY->mumu, jet->tau events
      cout << sqrt(error2OnQCD) << " <== W + TTb + DY(1,2)" << endl;      

      //  OS/SS ratio
      h1->Scale(OStoSSRatioQCD);
      cout << "After removing the expected contribution from W+jets and rescaling by " << OStoSSRatioQCD << " we expect " 
	   << h1->Integral() << " events from QCD processes" << endl;

      hParameters->SetBinContent(5, SSevents);
      hParameters->SetBinContent(6, h1->Integral()/SSevents);

      cout << "Total unceratinty from bkg subtraction in SS region is " << sqrt(error2OnQCD) << endl;
      float totalRelErrorOnQCD = 0.02 + sqrt(error2OnQCD)/h1->Integral(); //0.02 ==> uncertainty on OS/SS ratio
      hParameters->SetBinContent(7,totalRelErrorOnQCD);

      ////////////////////////////////////////////////

      hParameters->SetBinContent(8,ExtrapolationFactorZ);
      hParameters->SetBinContent(9,ErrorExtrapolationFactorZ);
      hParameters->SetBinContent(10,ExtrapolationFactorZDataMC);
      hParameters->SetBinContent(11,SSIsoToSSAIsoRatioQCD);


    }
    else{

      if( ! (it->first).Contains("Embed") ){

	if((it->first).Contains("DYToTauTau")){
	  currentTree->Draw(variable+">>"+h1Name, "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*sbinPZetaRelInclusive);  
	  float madgraphNoMEtCut = h1->Integral();
	  h1->Reset();
	  currentTree->Draw(variable+">>"+h1Name, "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*sbinInclusive);
	  madgraphMEtCutEff = h1->Integral()/madgraphNoMEtCut;
	  cout << "Efficiency of antiW cut on madgraph " << madgraphMEtCutEff << endl;
	  h1->Reset();
	  currentTree->Draw(variable+">>"+h1Name, "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*sbin);
	}
	else if((it->first).Contains("W3Jets")){

	  currentTree->Draw(variable+">>"+h1Name,     "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  hW3JetsKeys = new TH1Keys("hW3JetsKeys","W+3jets smoothed", nBins , bins.GetArray());
	  currentTree->Draw(variable+">>hW3JetsKeys", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  cout << "Keys for W3Jets filled with integral " << hW3JetsKeys->Integral() << " and entries " << hW3JetsKeys->GetEntries() << endl;

	  cout << "Filling histo with loose taus... if vbf, apply loose vbf cut" << endl;
	  h1->Reset();
	  if((selection_.Contains("vbf")) || selection_.Contains("twoJets")){

	    if( selection_.Contains("twoJets") ){
	      vbfLoose = twoJets;
	      vbf      = twoJets;
	    }

	    currentTree->Draw(variable+">>"+h1Name,"(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinMtisoInclusive&&vbfLoose)); // <------- sbinLtisoInclusive LOOSE TAU ID?
	    hW3JetsLooseTauIso->Add(h1,1.0);

	    // NORMALIZE W3JETS ON DATA

	    //TH1F* hHelp = (TH1F*)h1->Clone("hHelp");
	    TH1F* hHelp = new TH1F( "hHelp" ,"" , nBins , bins.GetArray());

	    // SIDEBAND OS: get extrap factor
	    hHelp->Reset();
	    currentTree->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelInclusive&&apZ&&vbfLoose));
	    cout << hHelp->Integral() << endl;
	    float W3JetsVBFSdb    = hHelp->Integral();
	    float W3JetsVBFSdbInt = hHelp->GetEntries();
	    hHelp->Reset();
	    currentTree->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelInclusive&&pZ&&vbfLoose));
	    float W3JetsVBFSgn    = hHelp->Integral();
	    float W3JetsVBFSgnInt = hHelp->GetEntries();
	    cout << "W3jets MC in VBF(rel): using " << W3JetsVBFSdbInt << " events from sideband and " << W3JetsVBFSgnInt << " from signal region" << endl;
	    float ExtrapFactorW3Jets = W3JetsVBFSdb/W3JetsVBFSgn;
	    cout << " ==> ratio sdb/sgn = " <<  ExtrapFactorW3Jets << endl;
	    hHelp->Reset();
            currentTree->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelInclusive&&vbfLoose));
            cout << hHelp->Integral() << endl;
            float W3JetsVBFAll    = hHelp->Integral();
            float W3JetsVBFAlllInt = hHelp->GetEntries();

	    // SIDEBAND2 OS: get extrap factor (2)
	    hHelp->Reset();
            currentTree->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelInclusive&&apZ2&&vbfLoose));
            cout << hHelp->Integral() << endl;
            float W3JetsVBFSdb2    = hHelp->Integral();
            float W3JetsVBFSdbInt2 = hHelp->GetEntries();
            hHelp->Reset();
            currentTree->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelInclusive&&pZ&&vbfLoose));
            float W3JetsVBFSgn2    = hHelp->Integral();
            float W3JetsVBFSgnInt2 = hHelp->GetEntries();
            cout << "W3jets MC in VBF(rel): using " << W3JetsVBFSdbInt2 << " events from sideband (2) and " << W3JetsVBFSgnInt2 << " from signal region" << endl;
            float ExtrapFactorW3Jets2 = W3JetsVBFSdb2/W3JetsVBFSgn2;
            cout << " ==> ratio sdb2/sgn = " <<  ExtrapFactorW3Jets2 << endl;


	    // SIDEBAND SS: get extrap fact SS
	    hHelp->Reset();
            currentTree->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelSSInclusive&&apZ2&&vbfLoose));
	    cout << hHelp->Integral() << endl;
            float W3JetsVBFSdbSS    = hHelp->Integral();
            float W3JetsVBFSdbIntSS = hHelp->GetEntries();
            hHelp->Reset();
            currentTree->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelSSInclusive&&pZ&&vbfLoose));
	    float W3JetsVBFSgnSS    = hHelp->Integral();
            float W3JetsVBFSgnIntSS = hHelp->GetEntries();
            cout << "W3jets MC (SS) in VBF(rel): using " << W3JetsVBFSdbIntSS << " events from sideband SS (2) and " << W3JetsVBFSgnIntSS << " from signal region" << endl;
            float ExtrapFactorW3JetsSS = W3JetsVBFSdbSS/W3JetsVBFSgnSS;
	    cout << " ==> ratio sdb/sgn = " <<  ExtrapFactorW3JetsSS << endl;



	    // now, normalize to sidebands
	    cout << "#######  Start normalization OS from sidebands #######" << endl;
	    hHelp->Reset();
	    data->Draw(variable+">>hHelp",sbinPZetaRelInclusive&&apZ2&&vbf);
	    float DataVBFSdb = hHelp->Integral();
	    cout << "In VBF region, I get " << DataVBFSdb << " events in the high Mt sdb" << endl;

	    hHelp->Reset();
	    backgroundDYTauTau->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelInclusive&&apZ2&&vbf));
	    float dytotauVBFSdb = hHelp->Integral()*Lumi*hltEff_;
	    cout << "In VBF region, I expect " << dytotauVBFSdb << " events in the high Mt sdb from DYtotautau" << endl;
	    hHelp->Reset();

	    hHelp->Reset();
	    backgroundDYMutoTau->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelInclusive&&apZ2&&vbf));
	    float dymutotauVBFSdb = hHelp->Integral()*Lumi*hltEff_;
	    cout << "In VBF region, I expect " << dymutotauVBFSdb << " events in the high Mt sdb from DYmutotau" << endl;
	    hHelp->Reset();

	    hHelp->Reset();
	    backgroundDYJtoTau->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelInclusive&&apZ2&&vbf));
	    float dyjtotauVBFSdb = hHelp->Integral()*Lumi*hltEff_;
	    cout << "In VBF region, I expect " << dyjtotauVBFSdb << " events in the high Mt sdb from DYjtotau" << endl;
	    hHelp->Reset();

	    hHelp->Reset();
	    backgroundOthers->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelInclusive&&apZ2&&vbf));
	    float othersVBFSdb = hHelp->Integral()*Lumi*hltEff_;
	    cout << "In VBF region, I expect " << othersVBFSdb << " events in the high Mt sdb from others" << endl;
	    hHelp->Reset();

	    backgroundTTbar->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelInclusive&&apZ2&&vbf));
	    float TTbarVBFSdb = hHelp->Integral()*Lumi*hltEff_*TTxsectionRatio;
	    cout << "In VBF region, I expect " << TTbarVBFSdb << " events in the high Mt sdb from TTbar" << endl;
	    hHelp->Reset();

	    //////////////////  Here, normalize TTbar from sideband 

	    backgroundTTbar->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelInclusive&&apZ2&&vbf&&TCut("nJets20BTagged>0")));
	    float TTbarVBFSdbBtag = hHelp->Integral()*Lumi*hltEff_*TTxsectionRatio;
	    hHelp->Reset();
	    backgroundTTbar->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelInclusive&&apZ2&&vbf&&TCut("nJets20BTagged>-999")));
	    float TTbarVBFSdbAnyBtag = hHelp->Integral()*Lumi*hltEff_*TTxsectionRatio;
	    float ratiobTagToAny = TTbarVBFSdbBtag/TTbarVBFSdbAnyBtag;
	    cout << "The ratio between events with >0 bTagged jets to any in TTJets MC is " << ratiobTagToAny << endl;
	    hHelp->Reset();
	    data->Draw(variable+">>hHelp",sbinPZetaRelInclusive&&apZ2&&vbf&&TCut("nJets20BTagged>0"));
	    float DataVBFSdbBtag = hHelp->Integral();
	    cout << "In data, I measure " << DataVBFSdbBtag << " events b-tagged" << endl;
	    scaleFactorTTinVBF = DataVBFSdbBtag/TTbarVBFSdbBtag;
	    cout << "The SF Data/MC ratio is therefore " << scaleFactorTTinVBF << endl;
	    cout << " ==> TTbar prediction in sideband goes from " << TTbarVBFSdbAnyBtag << " to " <<  DataVBFSdbBtag/ratiobTagToAny << endl;

	    hHelp->Reset();
            backgroundTTbar->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelInclusive&&apZ2&&vbfLoose));
            cout << hHelp->Integral() << endl;
            float TTJetsVBFSdb2    = hHelp->Integral();
            float TTJetsVBFSdbInt2 = hHelp->GetEntries();
            hHelp->Reset();
            backgroundTTbar->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelInclusive&&pZ&&vbfLoose));
            float TTJetsVBFSgn2    = hHelp->Integral();
            float TTJetsVBFSgnInt2 = hHelp->GetEntries();
            cout << "TT MC in VBF(rel): using " << TTJetsVBFSdbInt2 << " events from sideband (2) and " << TTJetsVBFSgnInt2 << " from signal region" << endl;
            float ExtrapFactorTTJets2 = TTJetsVBFSdb2/TTJetsVBFSgn2;
            cout << " ==> ratio sdb2/sgn = " <<  ExtrapFactorTTJets2 << endl;

	    ///////////////////////////////////////////////////////////////////////////////////////

	    data->Draw(variable+">>hHelp",(TCut(scaleFactMu.c_str()))*(sbinPZetaRelSSaIsoInclusive&&apZ2&&vbf));
	    float QCDVBFSdb = hHelp->Integral()*OStoSSRatioQCD;
	    cout << "In VBF region, I measure " << QCDVBFSdb << " events in the high Mt sdb from anti-isolated events" << endl;
	    cout << "Subtracting the backgrounds from the sideband..." << endl;
	    //DataVBFSdb -= TTbarVBFSdb;
	    DataVBFSdb -= (DataVBFSdbBtag/ratiobTagToAny);
	    DataVBFSdb -= QCDVBFSdb;
	    DataVBFSdb -= dytotauVBFSdb;
	    DataVBFSdb -= dymutotauVBFSdb;
	    DataVBFSdb -= dyjtotauVBFSdb;
	    DataVBFSdb -= othersVBFSdb;
	    float normalization = ExtrapFactorW3Jets>0 ? DataVBFSdb/ExtrapFactorW3Jets : -99;
	    cout << "Estimation of W+jets in the VBF category is " << normalization << endl;

	    hW3JetsLooseTauIso->Scale(normalization/hW3JetsLooseTauIso->Integral());

	    cout << "#######  Start normalization SS from sidebands #######" << endl;

	    hHelp->Reset();
            data->Draw(variable+">>hHelp",sbinPZetaRelSSInclusive&&apZ2&&vbf);
            float DataVBFSdbSS = hHelp->Integral();
            cout << "In VBF SS region, I get " << DataVBFSdbSS << " events in the high Mt sdb" << endl;

	    hHelp->Reset();
	    backgroundDYMutoTau->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelSSInclusive&&apZ2&&vbf));
	    float dymutotauVBFSdbSS = hHelp->Integral()*Lumi*hltEff_;
	    cout << "In VBF region, I expect " << dymutotauVBFSdbSS << " events in the high Mt sdb SS from DYmutotau" << endl;
	    hHelp->Reset();

	    hHelp->Reset();
	    backgroundDYJtoTau->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelSSInclusive&&apZ2&&vbf));
	    float dyjtotauVBFSdbSS = hHelp->Integral()*Lumi*hltEff_;
	    cout << "In VBF region, I expect " << dyjtotauVBFSdbSS << " events in the high Mt sdb SS from DYjtotau" << endl;
	    hHelp->Reset();

	    hHelp->Reset();
	    backgroundOthers->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelSSInclusive&&apZ2&&vbf));
	    float othersVBFSdbSS = hHelp->Integral()*Lumi*hltEff_;
	    cout << "In VBF region, I expect " << othersVBFSdbSS << " events in the high Mt sdb SS from others" << endl;
	    hHelp->Reset();

            hHelp->Reset();
            backgroundTTbar->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelSSInclusive&&apZ2&&vbf));
            float TTbarVBFSdbSS = hHelp->Integral()*Lumi*hltEff_*TTxsectionRatio;
            cout << "In VBF SS region, I expect " << TTbarVBFSdbSS << " events in the high Mt sdb from TTbar" << endl;
            hHelp->Reset();

	    //////////////////  Here, normalize TTbar from sideband 

	    backgroundTTbar->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelSSInclusive&&apZ2&&vbf&&TCut("nJets20BTagged>0")));
	    float TTbarVBFSdbSSBtag = hHelp->Integral()*Lumi*hltEff_*TTxsectionRatio;
	    hHelp->Reset();
	    backgroundTTbar->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelSSInclusive&&apZ2&&vbf&&TCut("nJets20BTagged>-999")));
	    float TTbarVBFSdbSSAnyBtag = hHelp->Integral()*Lumi*hltEff_*TTxsectionRatio;
	    float ratiobTagToAnySS = TTbarVBFSdbSSBtag/TTbarVBFSdbSSAnyBtag;
	    cout << "The ratio between events with >0 bTagged jets to any in TTJets MC in SS is " << ratiobTagToAnySS << endl;
	    hHelp->Reset();
	    data->Draw(variable+">>hHelp",sbinPZetaRelSSInclusive&&apZ2&&vbf&&TCut("nJets20BTagged>0"));
	    float DataVBFSdbSSBtag = hHelp->Integral();
	    cout << "In data, I measure " << DataVBFSdbSSBtag << " SS events b-tagged" << endl;
	    cout << "The SF Data/MC ratio is therefore " << DataVBFSdbSSBtag/TTbarVBFSdbSSBtag << endl;
	    cout << " ==> TTbar prediction in sideband SS goes from " << TTbarVBFSdbSSAnyBtag << " to " <<  DataVBFSdbSSBtag/ratiobTagToAnySS << endl;

	    ///////////////////////////////////////////////////////////////////////////////////////

            data->Draw(variable+">>hHelp",(TCut(scaleFactMu.c_str()))*(sbinPZetaRelSSaIsoInclusive&&apZ2&&vbf));
            float QCDVBFSdbSS = hHelp->Integral();
            cout << "In VBF SS region, I measure " << QCDVBFSdbSS << " events in the high Mt sdb from anti-isolated events" << endl;
	    hHelp->Reset();
            data->Draw(variable+">>hHelp",sbinPZetaRelSSInclusive&&pZ&&vbf);
            float DataVBFSgnSS = hHelp->Integral();
            float QCDSgnSS = DataVBFSgnSS - (DataVBFSdbSS - DataVBFSdbSSBtag/ratiobTagToAnySS /*TTbarVBFSdbSS*/-QCDVBFSdbSS)/ExtrapFactorW3JetsSS ;
            cout << "In VBF SS region, I measure " << DataVBFSgnSS << endl;
            cout << "The bkg estimation for QCD is therefore: " << DataVBFSgnSS << " - " << (DataVBFSdbSS-TTbarVBFSdbSS-QCDVBFSdbSS)/ExtrapFactorW3JetsSS << " = "
                 << QCDSgnSS << endl;

	    hParameters->SetBinContent(13, ExtrapFactorW3Jets);
	    hParameters->SetBinContent(14, DataVBFSdb);
	    hParameters->SetBinContent(15, TTbarVBFSdb);
	    hParameters->SetBinContent(16, QCDVBFSdb);
	    hParameters->SetBinContent(17, ExtrapFactorW3JetsSS);
	    hParameters->SetBinContent(18, DataVBFSdbSS);
	    hParameters->SetBinContent(19, TTbarVBFSdbSS);
	    hParameters->SetBinContent(20, QCDVBFSdbSS);
	    hParameters->SetBinContent(21, DataVBFSgnSS);
            hParameters->SetBinContent(22, QCDSgnSS);
	    hParameters->SetBinContent(23, ExtrapFactorW3Jets2);
            hParameters->SetBinContent(24, DataVBFSdb);
            hParameters->SetBinContent(25, TTbarVBFSdb);
            hParameters->SetBinContent(26, QCDVBFSdb);


	    delete hHelp;
	  } // endif sel=vbf or 2jets
	  else{ 
	    currentTree->Draw(variable+">>"+h1Name,"(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*sbinLtiso);
	    hW3JetsLooseTauIso->Add(h1,1.0);
	  }

	  cout << "Filling histo with loose taus and FR" << endl;
	  h1->Reset();
	  if(selection_.Contains("vbf")){
	    currentTree->Draw(variable+">>"+h1Name,"(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*TCut(scaleFactTauW.c_str())*(sbinLtisoInclusive&&vbfLoose));
	    hW3JetsLooseTauIsoFR->Add(h1,1.0);

	    TH1F* hHelp = new TH1F( "hHelp" ,"" , nBins , bins.GetArray());

	    hHelp->Reset();
	    currentTree->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*TCut(scaleFactTauW.c_str())*(sbinPZetaRelLtisoInclusive&&apZ&&vbfLoose));
	    float W3JetsVBFSdb    = hHelp->Integral();
	    float W3JetsVBFSdbInt = hHelp->GetEntries();
	    hHelp->Reset();
	    currentTree->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*TCut(scaleFactTauW.c_str())*(sbinPZetaRelLtisoInclusive&&pZ&&vbfLoose));
	    float W3JetsVBFSgn    = hHelp->Integral();
	    float W3JetsVBFSgnInt = hHelp->GetEntries();
	    cout << "W3jets MC in VBF(rel): using " << W3JetsVBFSdbInt << " events from sideband and " << W3JetsVBFSgnInt << " from signal region" << endl;
	    cout << " ==> ratio sdb/sgn = " << W3JetsVBFSdb/W3JetsVBFSgn << endl;

	    hHelp->Reset();
	    data->Draw(variable+">>hHelp",sbinPZetaRelInclusive&&apZ&&vbf);
	    float DataVBFSdb = hHelp->Integral();
	    cout << "In VBF region, I get " << DataVBFSdb << " events in the high Mt sdb" << endl;
	    hHelp->Reset();
	    backgroundTTbar->Draw(variable+">>hHelp","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*(sbinPZetaRelInclusive&&apZ&&vbf));
	    float TTbarVBFSdb = hHelp->Integral()*Lumi*hltEff_*TTxsectionRatio;
	    cout << "In VBF region, I expect " << TTbarVBFSdb << " events in the high Mt sdb from TTbar" << endl;
	    hHelp->Reset();
	    currentTree->Draw(variable+">>hHelp",(TCut(scaleFactMu.c_str()))*(sbinPZetaRelSSaIsoInclusive&&apZ&&vbf));
	    float QCDVBFSdb = hHelp->Integral()*OStoSSRatioQCD;
	    cout << "In VBF region, I measure " << QCDVBFSdb << " events in the high Mt sdb from anti-isolated events" << endl;

	    cout << "Subtracting the backgrounds from the sideband..." << endl;
	    DataVBFSdb -= TTbarVBFSdb;
	    DataVBFSdb -= QCDVBFSdb;
	    float normalization = W3JetsVBFSdb/W3JetsVBFSgn>0 ? DataVBFSdb/(W3JetsVBFSdb/W3JetsVBFSgn) : -99;
	    cout << "Estimation of W+jets in the VBF category is " << normalization << endl;

	    hW3JetsLooseTauIsoFR->Scale(normalization/hW3JetsLooseTauIso->Integral());

	    delete hHelp;

	  }
	  else{
	    currentTree->Draw(variable+">>"+h1Name,"(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*TCut(scaleFactTauW.c_str())*sbinLtiso);
	    hW3JetsLooseTauIsoFR->Add(h1,1.0);
	  }
	  h1->Reset();
	  currentTree->Draw(variable+">>"+h1Name,     "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);	  
	}
	else if((it->first).Contains("WJets")){
	  currentTree->Draw(variable+">>"+h1Name,     "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  hWKeys = new TH1Keys("hWKeys","W+jets smoothed", nBins , bins.GetArray());
	  currentTree->Draw(variable+">>hWKeys", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  cout << "Keys for WJets filled with integral " << hWKeys->Integral() << " and entries " << hWKeys->GetEntries() << endl;
	}
	else if((it->first).Contains("DYMutoTau")){
	  currentTree->Draw(variable+">>"+h1Name,     "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  hZmmKeys = new TH1Keys("hZmmKeys","Z+jets, mu to tau smoothed", nBins , bins.GetArray());
	  currentTree->Draw(variable+">>hZmmKeys", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  cout << "Keys for Zmm filled with integral " << hZmmKeys->Integral() << " and entries " << hZmmKeys->GetEntries() << endl;
	  if( selection_.Contains("vbf") ){
	    TH1F* hHelp = new TH1F( "hHelp" ,"" , nBins , bins.GetArray());
	    currentTree->Draw(variable+">>hHelp",     "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*(sbinInclusive&&twoJets));
	    hHelp->Scale(ExtrapLFakeInclusive*ExtrapolationFactorZ/hHelp->Integral());
	    hZmmLoose->Add(hHelp,1.0);
	    delete hHelp;
	  }
	}
	else if((it->first).Contains("DYJtoTau")){
	  currentTree->Draw(variable+">>"+h1Name,     "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  hZmjKeys = new TH1Keys("hZmjKeys","Z+jets, jet to tau smoothed", nBins , bins.GetArray());
	  currentTree->Draw(variable+">>hZmjKeys", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  cout << "Keys for Zmj filled with integral " << hZmjKeys->Integral() << " and entries " << hZmjKeys->GetEntries() << endl;
	  if( selection_.Contains("vbf") ){
            TH1F* hHelp = new TH1F( "hHelp" ,"" , nBins , bins.GetArray());
            currentTree->Draw(variable+">>hHelp",     "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*(sbinInclusive&&twoJets));
            hHelp->Scale(ExtrapJFakeInclusive*ExtrapolationFactorZ/hHelp->Integral());
            hZmjLoose->Add(hHelp,1.0);
            delete hHelp;
          }
	}
	else if((it->first).Contains("Others")){
	  currentTree->Draw(variable+">>"+h1Name,     "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  hVVKeys = new TH1Keys("hVVKeys","Others smoothed", nBins , bins.GetArray());
	  currentTree->Draw(variable+">>hVVKeys", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  cout << "Keys for VV filled with integral " << hVVKeys->Integral() << " and entries " << hVVKeys->GetEntries() << endl;
	}

	// here, we fill histos for fake-rate
	else if((it->first).Contains("Data")){

	  h1->Reset();
	  currentTree->Draw(variable+">>"+h1Name, sbin);
	  hData->Add(h1,1.0);

	  cout << "DATA === " << hData->Integral() << endl;

	  cout << "Filling histo with anti-iso muons and loose taus" << endl;
	  if(selection_.Contains("vbf")){
	    currentTree->Draw(variable+">>"+h1Name,sbinSSaIsoLtisoInclusive&&vbfLoose);
	    hDataAntiIsoLooseTauIso->Add(h1,1.0);
	    h1->Reset();
	    currentTree->Draw(variable+">>"+h1Name,sbinSSaIsoLtiso);
	    float normalization = h1->Integral();
	    hDataAntiIsoLooseTauIso->Scale(normalization/hDataAntiIsoLooseTauIso->Integral());
	  }
	  else{
	    currentTree->Draw(variable+">>"+h1Name,sbinSSaIsoLtiso);
	    hDataAntiIsoLooseTauIso->Add(h1,1.0);
	  }

	  // template for QCD : anti-iso mu && loose-tau X fake-rate(mu) X fake-rate(tau)
	  cout << "Filling histo with anti-iso muons and loose taus and FR" << endl;
	  h1->Reset();
	  if(selection_.Contains("vbf")){
	    hDataAntiIsoLooseTauIsoFR->Add(h1,1.0);
	    h1->Reset();
	    float normalization = h1->Integral();
	  }
	  else{
	    hDataAntiIsoLooseTauIsoFR->Add(h1,1.0);
	  }

	  cout << "Filling histo with anti-iso muons and loose taus and FR (Up)" << endl;
	  h1->Reset();
	  if(selection_.Contains("vbf")){
	    hDataAntiIsoLooseTauIsoFRUp->Add(h1,1.0);
	    h1->Reset();
	    float normalization = h1->Integral()*OStoSSRatioQCD;
	  }
	  else{
	    hDataAntiIsoLooseTauIsoFRUp->Add(h1,OStoSSRatioQCD);
	  }

	  cout << "Filling histo with anti-iso muons and loose taus and FR (Down)" << endl;
	  h1->Reset();
	  if(selection_.Contains("vbf")){
	    hDataAntiIsoLooseTauIsoFRDown->Add(h1,1.0);
	    h1->Reset();
	    float normalization = h1->Integral()*OStoSSRatioQCD;
	  }
	  else{
	    hDataAntiIsoLooseTauIsoFRDown->Add(h1,OStoSSRatioQCD);
	  }

	  // template for QCD : anti-iso mu && loose-tau X fake-rate(mu) X fake-rate(tau)
	  cout << "Filling histo with anti-iso muons and tight taus" << endl;
	  h1->Reset();
	  if(selection_.Contains("vbf")){
	    hDataAntiIso->Add(h1,1.0);
	    h1->Reset();
	    float normalization = h1->Integral()*OStoSSRatioQCD;
	  }
	  else{
	    hDataAntiIso->Add(h1,OStoSSRatioQCD);
	  }

	  cout << "Filling histo with anti-iso muons and tight taus and FR" << endl;
	  h1->Reset();
	  if(selection_.Contains("vbf")){
	    hDataAntiIsoFR->Add(h1,1.0);
	    h1->Reset();
	    float normalization = h1->Integral()*OStoSSRatioQCD;
	  }
	  else{
	    hDataAntiIsoFR->Add(h1,OStoSSRatioQCD);
	  }

	  cout << "Filling histo with anti-iso muons and tight taus and FR (Up)" << endl;
	  h1->Reset();
	  if(selection_.Contains("vbf")){
	    hDataAntiIsoFRUp->Add(h1,1.0);
	    h1->Reset();
	    float normalization = h1->Integral()*OStoSSRatioQCD;
	  }
	  else{
	    hDataAntiIsoFRUp->Add(h1,OStoSSRatioQCD);
	  }

	  cout << "Filling histo with anti-iso muons and tight taus and FR (Down)" << endl;
	  h1->Reset();
	  if(selection_.Contains("vbf")){
	    hDataAntiIsoFRDown->Add(h1,1.0);
	    h1->Reset();
	    float normalization = h1->Integral()*OStoSSRatioQCD;
	  }
	  else{
	    hDataAntiIsoFRDown->Add(h1,OStoSSRatioQCD);
	  }
	}

	else if((it->first).Contains("LooseIso")){
	  currentTree->Draw(variable+">>"+h1Name,    sbinSSlIso);

	  TH1F* hHelp = new TH1F( "hHelp" ,"" , nBins , bins.GetArray());
	  currentTree->Draw(variable+">>hHelp",    sbinSSlIso2);
	  hLooseIso2->Add(hHelp,1.0);
	  hHelp->Reset();
	  currentTree->Draw(variable+">>hHelp",    sbinSSlIso3);
	  hLooseIso3->Add(hHelp,1.0);
	  delete hHelp;

	  hLooseIsoKeys = new TH1Keys("hLooseIsoKeys","Loose Iso smoothed", nBins , bins.GetArray());
	  if(  ((selection_.Contains("vbf")) || 
		selection_.Contains("boost") ||
		selection_.Contains("vh")) )
	    currentTree->Draw(variable+">>hLooseIsoKeys", sbinSSlIso);
	  cout << "Keys for LooseIso filled with integral " << hLooseIsoKeys->Integral() << " and entries " << hLooseIsoKeys->GetEntries() << endl;

	}
	else if((it->first).Contains("AntiIso")){
	  currentTree->Draw(variable+">>"+h1Name,    sbinSSaIso);
	  hAntiIsoKeys = new TH1Keys("hAntiIsoKeys","Anti Iso smoothed", nBins , bins.GetArray());
	  currentTree->Draw(variable+">>hAntiIsoKeys", sbinSSaIso);
	  cout << "Keys for AntiIso filled with integral " << hAntiIsoKeys->Integral() << " and entries " << hAntiIsoKeys->GetEntries() << endl;
	}
	else
	  currentTree->Draw(variable+">>"+h1Name, "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);

      } // if not "Embed" in Sample name

      else{ // Embed
	currentTree->Draw(variable+">>"+h1Name, "(HLTTau*HLTMu*1.000000)"*sbinPZetaRelEmbeddingInclusive);
	float embeddedNoMEtCut = h1->Integral();
	h1->Reset();
	currentTree->Draw(variable+">>"+h1Name, "(HLTTau*HLTMu*1.000000)"*sbinEmbeddingInclusive);

	embeddedMEtCutEff =  h1->Integral()/embeddedNoMEtCut;
	cout << "Efficiency of antiW cut on embedded " << embeddedMEtCutEff << endl;

	h1->Reset();
	currentTree->Draw(variable+">>"+h1Name, "(HLTTau*HLTMu*1.000000)"*sbinEmbedding);
      }

      /////////////
      // SCALING //
      /////////////

      // scale by correction factors
      if( ! ((it->first).Contains("Data") || 
	     (it->first).Contains("LooseIso") ||
	     (it->first).Contains("AntiIso")) ) 
	h1->Scale(Lumi*hltEff_);

      // scale by DATA/MC ratio in Zmumu sideband
      if((it->first).Contains("DY"))
	h1->Scale(lumiCorrFactor);

      // if W+jets, scale by extrapolation
      float sFWOS = ( selection_.Contains("0jet")  || 
		      selection_.Contains("boost")  || 
		      selection_.Contains("bTag")   || 
		      selection_.Contains("oneJet")  || 
		      selection_.Contains("twoJets")  || 
		      selection_.Contains("inclusive")) ? 
	OSWinSignalRegionDATA/OSWinSignalRegionMC : WcorrectionFactorOS;
  
      if((it->first).Contains("WJets")){

	TH1F* hHelp = new TH1F( "hHelp" ,"" , nBins , bins.GetArray());
	currentTree->Draw(variable+">>hHelp", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbinMtiso);
	hWLooseIso2->Add(hHelp,1.0);
	delete hHelp;

	if(selection_.Contains("vbf")){
	  sFWOS *= VbfExtrapolationFactorW;
	  cout << "Wjets will be rescaled by " << VbfExtrapolationFactorW << " according to the Z->mumu+j+vbf/Z->mumu+j ratio" << endl;
	}
	else if(selection_.Contains("boost")){
	  sFWOS *= BoostExtrapolationFactorW;
	  cout << "Wjets will be rescaled by " << BoostExtrapolationFactorW << " according to the Z->mumu+j+vbf/Z->mumu+j ratio" << endl;
	}	

	h1->Scale( sFWOS );
	hW->Add(h1,1.0);
	if(h1->Integral()>0) hWLooseIso2->Scale(h1->Integral()/hWLooseIso2->Integral());
      }

      // if DY->tautau, and vbf scale by ratio data/MC
      if((it->first).Contains("DYToTauTau")){
	h1->Scale( ExtrapolationFactorZDataMC );
      }

      // if DY->mumu, mu->tau, scale by fake-rate
      if((it->first).Contains("DYMutoTau")){

	float sF = MutoTauCorrectionFactor;

	sF *= ExtrapolationFactorZDataMC;

	h1->Scale(sF);
	hZmm->Add(h1,1.0);
	hZfakes->Add(h1,1.0);
      }

      // if DY->mumu, jet->tau, scale by fake-rate
      if((it->first).Contains("DYJtoTau")){

	float sF = JtoTauCorrectionFactor;

	sF *= ExtrapolationFactorZDataMC;

	h1->Scale(sF);
	hZmj->Add(h1,1.0); //<-----------------------
	hZfakes->Add(h1,1.0);
      }

    }
  
    /////////////////////////////////////////////////////////////////////////////////////

    // Build hEWK = DYmuTau+DYjetTau+Wjets+VV
    if( (it->first).Contains("DYMutoTau") ||  
	(it->first).Contains("DYJtoTau") || 
	(it->first).Contains("WJets") || 
	(it->first).Contains("Others") ) {
      hEWK->SetFillColor(kRed+2);
      hEWK->Add(h1,1.0);

      if( (it->first).Contains("Others") ){
	hVV->Add(h1,1.0);
	if(hVVKeys->Integral()>0) hVVKeys->Scale(hVV->Integral()/hVVKeys->Integral());
	hVVKeys->SetFillColor(kRed+2);
      }
    }

    if( (it->first).Contains("DYToTauTau") ) {
      hZtt->Add(h1,1.0);
      hZtt->SetFillColor(kOrange-4);
    }

    // Scale Embedded
    if( (it->first).Contains("Embedded") ) {

      hParameters->SetBinContent(12,embeddedMEtCutEff/madgraphMEtCutEff);

      h1->Scale( (ExtrapolationFactorZ*ExtrapDYInclusive)/h1->Integral());
      h1->Scale(embeddedMEtCutEff/madgraphMEtCutEff);

      hDataEmb->Add(h1,1.0);
      hDataEmb->SetFillColor(kOrange-4);
    }

    // Scale keys histograms
    if( (it->first).Contains("DYMutoTau") ) {
      if(hZmmKeys->Integral()>0) hZmmKeys->Scale(hZmm->Integral()/hZmmKeys->Integral());
      hZmmKeys->SetFillColor(kRed+2);
    }
    if( (it->first).Contains("DYJtoTau") ) {
      if(hZmjKeys->Integral()>0) hZmjKeys->Scale(hZmj->Integral()/hZmjKeys->Integral());
      hZmjKeys->SetFillColor(kRed+2);
    }
    if( (it->first).Contains("WJets") ) {
      if(hWKeys->Integral()>0) hWKeys->Scale(hW->Integral()/hWKeys->Integral());
      hWKeys->SetFillColor(kRed+2);
    }
    if( (it->first).Contains("W3Jets") ) {

      hW3Jets->Add(h1,1.0);
      if(hW3Jets->Integral()>0) hW3Jets->Scale(hW->Integral()/hW3Jets->Integral());
      hW3Jets->SetFillColor(kRed+2);
      if(hW3JetsKeys->Integral()>0) hW3JetsKeys->Scale(hW->Integral()/hW3JetsKeys->Integral());
      hW3JetsKeys->SetFillColor(kRed+2);
    }
    if( (it->first).Contains("LooseIso") ) {
      hLooseIso->Add(h1,1.0);
      if(hLooseIsoKeys->Integral()>0) hLooseIsoKeys->Scale(hLooseIso->Integral()/hLooseIsoKeys->Integral());
      hLooseIsoKeys->SetFillColor(0);
    }
    if( (it->first).Contains("AntiIso") ) {
      float sF = SSIsoToSSAIsoRatioQCD*OStoSSRatioQCD;
      cout << "Anti-isolated SS data scaled by " <<  SSIsoToSSAIsoRatioQCD << " ratio measured in inclusive sample ==>" << endl;
      cout << "   SS anti-isolated events = " << h1->GetEntries() << " ==> " << h1->Integral()*SSIsoToSSAIsoRatioQCD << " predicted in signal region" << endl;
      h1->Scale(sF);
      hAntiIso->Add(h1,1.0);
      if(hAntiIsoKeys->Integral()>0) hAntiIsoKeys->Scale(hAntiIso->Integral()/hAntiIsoKeys->Integral());
      hAntiIsoKeys->SetFillColor(0);
    }
    if( (it->first).Contains("TTbar") ) {
      hTTb->Add(h1,1.0);
      hTTb->Scale(TTxsectionRatio);
      if( selection_.Contains("vbf") )
	hTTb->Scale(scaleFactorTTinVBF);
      else
	hTTb->Scale(scaleFactorTTOS);
      hTTb->SetFillColor(kBlue-8);     
    }
    if( (it->first).Contains("SS") ) {
      hQCD->Add(h1,1.0);
      hQCD->SetFillColor(kMagenta-10);
    }
    if((it->first).Contains("qqH"+mH_)){
      hSgn1->Add(h1,1.0);
      hSgn1->Scale(magnifySgn_);
      h1->Scale(magnifySgn_);
      hSgn1->SetLineWidth(2);
      h1->SetFillColor(kBlack);
      h1->SetFillStyle(3004);
      h1->SetLineColor(kBlack);
    }
    if((it->first).Contains("ggH"+mH_)){
      hSgn2->Add(h1,1.0);
      hSgn2->Scale(magnifySgn_);
      h1->Scale(magnifySgn_);
      hSgn2->SetLineWidth(2);
      h1->SetFillColor(kBlack);
      h1->SetFillStyle(3005);
      h1->SetLineColor(kBlack);
    }
    if((it->first).Contains("VH"+mH_)){
      hSgn3->Add(h1,1.0);
      hSgn3->Scale(magnifySgn_);
      h1->Scale(magnifySgn_);
      hSgn3->SetLineWidth(2);
      h1->SetFillColor(kBlack);
      h1->SetFillStyle(3005);
      h1->SetLineColor(kBlack);
    }
    
    for(int im=0 ; im<nMasses ; im++) {

      if((it->first).Contains("ggH"+hMasses[im])){
	hggH[im]->Add(h1,1.0);
	hggH[im]->SetLineWidth(2);
      }
      if((it->first).Contains("qqH"+hMasses[im])){
	hggH[im]->Add(h1,1.0);
	hggH[im]->SetLineWidth(2);
      }
      if((it->first).Contains("VH"+hMasses[im])){
	hggH[im]->Add(h1,1.0);
	hggH[im]->SetLineWidth(2);
      }
    }

    ////////////////////////////////////////////////////////////
    if((it->first).Contains("SUSY")){
      TH1F* histoSusy =  (mapSUSYhistos.find( (it->first) ))->second;
      histoSusy->Add(h1,1.0);
      histoSusy->SetLineWidth(2);
    }
    ////////////////////////////////////////////////////////////

    if((it->first).Contains("Data")){
      hData->Sumw2();
      hData->SetMarkerStyle(20);
      hData->SetMarkerSize(1.2);
      hData->SetMarkerColor(kBlack);
      hData->SetLineColor(kBlack);
      hData->SetXTitle(XTitle_+" ("+Unities_+")");
      hData->SetYTitle(Form(" Events/(%.1f %s)", hData->GetBinWidth(1), Unities_.Data() ) );
      hData->SetTitleSize(0.04,"X");
      hData->SetTitleSize(0.05,"Y");
      hData->SetTitleOffset(0.95,"Y");
    }

    if(VERBOSE) cout<<(it->first) << " ==> "  // PROBLEM TO SOLVE
		    << h1->Integral() << " +/- " 
		    << TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries())
		    << endl;

    char* c = new char[50];
    if(h1->Integral()>=10) 
      sprintf(c,"$%.0f\\pm%.0f$",h1->Integral(),  TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries()));
    else if(h1->Integral()>=1)
      sprintf(c,"$%.1f\\pm%.1f$",h1->Integral(),  TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries()));
    else if(h1->Integral()>=0.1)
      sprintf(c,"$%.2f\\pm%.2f$",h1->Integral(),  TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries()));
    else if(h1->Integral()>=0.01)
      sprintf(c,"$%.3f\\pm%.3f$",h1->Integral(),  TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries()));
    else
      sprintf(c,"$%.5f\\pm%.5f$",h1->Integral(),  TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries()));
    out << string(c) << "  //" << (it->first) << endl;
    delete c;
  }

  out.close();

  // all signal summed together:
  hSgn->Add(hSgn1,hSgn2,1,1);
  hSgn->Add(hSgn3);
  hSgn->SetFillColor(0);
  hSgn->SetLineColor(kBlue);
  hSgn->SetLineWidth(2);
  hSgn->SetLineStyle(kDashed);

  // adding alltogether
  hSiml->Add(hQCD,1.0);
  hSiml->Add(hEWK,1.0);
  hSiml->Add(hTTb,1.0);
  if(useEmbedding_)
    hSiml->Add(hDataEmb,1.0);
  else
    hSiml->Add(hZtt,1.0);


  //Adding to the stack
  if(selection_.Contains("vbf") || selection_.Contains("twoJets") ){
    hDataAntiIsoFR->SetFillColor(kMagenta-10);
    aStack->Add(hDataAntiIsoFR);
  }
  else
    aStack->Add(hQCD);

  if(selection_.Contains("vbf") || selection_.Contains("twoJets") ){
    hZfakes->SetFillColor(kRed+2);
    hW3JetsLooseTauIso->SetFillColor(kRed+2);
    hW3JetsLooseTauIso->SetLineColor(kRed+2);
    aStack->Add(hW3JetsLooseTauIso);
    aStack->Add(hZfakes);
  }
  else
    aStack->Add(hEWK);

  aStack->Add(hTTb);
  if(useEmbedding_)
    aStack->Add(hDataEmb);
  else
    aStack->Add(hZtt);
  if(!logy_)
    aStack->Add(hSgn);

  // legend
  leg->AddEntry(hData,"Observed","P");
  leg->AddEntry(hSgn, Form("(%.0fx) H#rightarrow#tau#tau m_{H}=",magnifySgn_)+mH_,"F" );
  if(useEmbedding_)
    leg->AddEntry(hDataEmb,"Z#rightarrow#tau#tau (embedded)","F");
  else
    leg->AddEntry(hZtt,"Z#rightarrow#tau#tau","F"); 
  leg->AddEntry(hTTb,"t#bar{t}","F");
  leg->AddEntry(hEWK,"Electroweak","F");
  leg->AddEntry(hQCD,"QCD","F");
  
  hData->Draw("P");
  aStack->Draw("HISTSAME");
  hData->Draw("PSAME");
  
  TH1F* hStack = (TH1F*)aStack->GetHistogram();
  hStack->SetXTitle(XTitle_+" ("+Unities_+")");
  hStack->SetYTitle(Form(" Events/(%.0f %s)", hStack->GetBinWidth(1), Unities_.Data() ) );
  hStack->SetTitleSize(0.04,"X");
  hStack->SetTitleSize(0.05,"Y");
  hStack->SetTitleOffset(0.95,"Y");
  if(!logy_)
    hData->SetAxisRange(0.0, TMath::Max( hData->GetMaximum(), hSiml->GetMaximum() )*maxY_ ,"Y");
  else
    hData->SetAxisRange(0.1, TMath::Max( hData->GetMaximum(), hSiml->GetMaximum() )*maxY_ ,"Y");
  aStack->Draw("HISTSAME");
  hData->Draw("PSAME");
  if(logy_)
    hSgn->Draw("HISTSAME");
  

  leg->Draw();

  pad2->cd();
  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");

  //TH1F* hRatio = new TH1F("hRatio", " ; ; #frac{(DATA-MC)}{#sqrt{DATA}}",
  //		  hStack->GetNbinsX(), 
  //		  hStack->GetXaxis()->GetXmin(), hStack->GetXaxis()->GetXmax());

  TH1F* hRatio = new TH1F( "hRatio" ," ; ; #frac{(DATA-MC)}{#sqrt{DATA}}" , nBins , bins.GetArray());
  hRatio->Reset();
       hRatio->SetXTitle("");
  hRatio->SetYTitle("#frac{(DATA-MC)}{#sqrt{DATA}}");

  hRatio->SetMarkerStyle(kFullCircle);
  hRatio->SetMarkerSize(0.8);
  hRatio->SetLabelSize(0.12,"X");
  hRatio->SetLabelSize(0.10,"Y");
  hRatio->SetTitleSize(0.12,"Y");
  hRatio->SetTitleOffset(0.36,"Y");
  TH1F* hError = (TH1F*)hRatio->Clone("hError");
  hError->Reset();
  hError->SetFillColor(kRed);
  hError->SetFillStyle(3004);
  hError->SetMarkerStyle(kDot);

  float uncertZtt = 0;
  if(selection_.Contains("0jet")  || selection_.Contains("inclusive") ||
     selection_.Contains("oneJet") || selection_.Contains("twoJets")){
    uncertZtt += (0.06  * 0.06) ; // Tau-ID 
    uncertZtt += (0.035 * 0.035); // Lumi 
  }
  else if(selection_.Contains("vbf")){
    uncertZtt += (0.06  * 0.06) ; // Tau-ID 
    uncertZtt += (0.035 * 0.035); // Lumi 
    uncertZtt += (0.10  * 0.10);  // Extrap. factor 
  }
  else if(selection_.Contains("boost")){
    uncertZtt += (0.06  * 0.06) ; // Tau-ID 
    uncertZtt += (0.035 * 0.035); // Lumi 
    uncertZtt += (0.10  * 0.10);  // Extrap. factor 
  }
  uncertZtt = TMath::Sqrt(uncertZtt);
  float uncertTTb = 0;
  if(selection_.Contains("0jet") || selection_.Contains("inclusive")||
     selection_.Contains("oneJet") || selection_.Contains("twoJets")){
    uncertTTb += (0.075 * 0.075) ; // xsection 
  }
  else if(selection_.Contains("vbf")){
    uncertTTb += (0.075 * 0.075) ; // xsection 
  }
  else if(selection_.Contains("boost")){
    uncertTTb += (0.075 * 0.075) ; // xsection 
  }
  uncertTTb = TMath::Sqrt(uncertTTb);
  float uncertEWK = 0;
  if(selection_.Contains("0jet") || selection_.Contains("inclusive") ||
     selection_.Contains("oneJet") || selection_.Contains("twoJets")){
    uncertEWK += (0.07 * 0.07) ; // extrapolation 
  }
  else if(selection_.Contains("vbf")){
    uncertEWK += (0.07 * 0.07) ; // extrapolation 
    uncertEWK += (0.10 * 0.10) ; // extrapolation 
  }
  else if(selection_.Contains("boost")){
    uncertEWK += (0.07 * 0.07) ; // extrapolation 
    uncertEWK += (0.10 * 0.10) ; // extrapolation  
  }
       uncertEWK = TMath::Sqrt(uncertEWK);
  float uncertQCD = 0;
  if(selection_.Contains("0jet") || selection_.Contains("inclusive") ||
     selection_.Contains("oneJet") || selection_.Contains("twoJets")){
    uncertQCD += (0.02 * 0.02) ; // extrapolation 
  }
  else if(selection_.Contains("vbf")){
    uncertQCD += (0.02 * 0.02) ; // extrapolation 
    uncertQCD += (0.10 * 0.10) ; // extrapolation 
  }
  else if(selection_.Contains("boost")){
    uncertQCD += (0.02 * 0.02) ; // extrapolation 
    uncertQCD += (0.10 * 0.10) ; // extrapolation  
  }
  uncertQCD = TMath::Sqrt(uncertQCD);

  float maxPull = 0.;
  for(int k = 0 ; k < hRatio->GetNbinsX(); k++){
    float pull = hData->GetBinContent(k) - hSiml->GetBinContent(k);
    if(hData->GetBinContent(k)>0)
      pull /= TMath::Sqrt(hData->GetBinContent(k));
    hRatio->SetBinContent(k, pull);
    hError->SetBinContent(k, 0.0);
    float error2 = 0.0;
    error2 += pow(hZtt->GetBinContent(k) * uncertZtt, 2);
    error2 += pow(hTTb->GetBinContent(k) * uncertTTb, 2);
    error2 += pow(hEWK->GetBinContent(k) * uncertEWK, 2);
    error2 += pow(hQCD->GetBinContent(k) * uncertQCD, 2);
    if(hData->GetBinContent(k)>0)
      hError->SetBinError(k, TMath::Sqrt(error2/hData->GetBinContent(k)) );
    else
      hError->SetBinError(k,0);
    if(TMath::Abs(pull) > maxPull)
      maxPull = TMath::Abs(pull);
  }
  hRatio->SetAxisRange(-1.2*maxPull,1.2*maxPull,"Y");
  hRatio->Draw("P");
  hError->Draw("E3SAME");

  TF1* line = new TF1("line","0",hRatio->GetXaxis()->GetXmin(),hStack->GetXaxis()->GetXmax());
  line->SetLineStyle(3);
  line->SetLineWidth(1.5);
  line->SetLineColor(kBlack);
  line->Draw("SAME");
  
  //return;
       
  c1->SaveAs(dirOut+"/plots/plot_muTau_mH"+mH_+"_"+selection_+"_"+analysis_+"_"+variable_+".png");
  //c1->SaveAs(dirOut+"/plots/plot_muTau_mH"+mH_+"_"+selection_+"_"+analysis_+"_"+variable_+".pdf");
  //c1->SaveAs(dirOut+"/plots/plot_muTau_mH"+mH_+"_"+selection_+"_"+analysis_+"_"+variable_+".ps");
  //c1->SaveAs(dirOut+"/plots/plot_muTau_mH"+mH_+"_"+selection_+"_"+analysis_+"_"+variable_+".eps");

  // templates for fitting
  TFile* fout = new TFile(dirOut+"/histograms/h_muTau_mH"+mH_+"_"+selection_+"_"+analysis_+"_"+variable_+".root","RECREATE");
  fout->cd();

  cout << "Starts writing" << endl;

  hSiml->Write();
  hQCD->Write();
  hZmm->Write();
  hZmmLoose->Write();
  hZmj->Write();
  hZmjLoose->Write();
  hZfakes->Write();
  hTTb->Write();
  hZtt->Write();
  hDataEmb->Write();
  hLooseIso->Write();
  hLooseIso2->Write();
  hLooseIso3->Write();
  hAntiIso->Write();
  hW->Write();
  hWLooseIso2->Write();
  if(hWKeys) hW3Jets->Write();

  cout << "Starts writing Keys" << endl;

  TH1* hWKeysHisto = hWKeys ? ((TH1Keys*)hWKeys)->GetHisto() : 0;
  if(hWKeysHisto){
    makeHistoFromDensity(hWKeysHisto,hW);
    hWKeysHisto->SetName("hWKeys");
    hWKeysHisto->Write();
  }
  TH1* hW3JetsKeysHisto = hW3JetsKeys ? ((TH1Keys*)hW3JetsKeys)->GetHisto() : 0;
  if(hW3JetsKeysHisto){
    makeHistoFromDensity(hW3JetsKeysHisto,hW3Jets);
    hW3JetsKeysHisto->SetName("hW3JetsKeys");
    hW3JetsKeysHisto->Write();
  }
  TH1* hZmmKeysHisto = hZmmKeys ? ((TH1Keys*)hZmmKeys)->GetHisto() : 0;
  if(hZmmKeysHisto){
    makeHistoFromDensity(hZmmKeysHisto,hZmm);
    hZmmKeysHisto->SetName("hZmmKeys");
    hZmmKeysHisto->Write();
  }

  TH1* hZmjKeysHisto = hZmjKeys ? ((TH1Keys*)hZmjKeys)->GetHisto() : 0;
  if(hZmjKeysHisto){
    makeHistoFromDensity(hZmjKeysHisto,hZmj);
    hZmjKeysHisto->SetName("hZmjKeys");
    hZmjKeysHisto->Write();
  }
  TH1* hLooseIsoKeysHisto = hLooseIsoKeys ? ((TH1Keys*)hLooseIsoKeys)->GetHisto() : 0;
  if(hLooseIsoKeysHisto){
    makeHistoFromDensity(hLooseIsoKeysHisto,hLooseIso);
    hLooseIsoKeysHisto->SetName("hLooseIsoKeys");
    hLooseIsoKeysHisto->Write();
  }
  TH1* hAntiIsoKeysHisto = hAntiIsoKeys ? ((TH1Keys*)hAntiIsoKeys)->GetHisto() : 0;
  if(hAntiIsoKeysHisto){
    makeHistoFromDensity(hAntiIsoKeysHisto,hAntiIso);
    hAntiIsoKeysHisto->SetName("hAntiIsoKeys");
    hAntiIsoKeysHisto->Write();
  }
  TH1* hVVKeysHisto = hVVKeys ? ((TH1Keys*)hVVKeys)->GetHisto() : 0;
  if(hVVKeysHisto){
    makeHistoFromDensity(hVVKeysHisto, hVV);
    hVVKeysHisto->SetName("hVVKeys");
    hVVKeysHisto->Write();
  }

  cout << "Starts writing other stuff" << endl;

  if(hVV) hVV->Write();
  if(hEWK) hEWK->Write();
  if(hSgn1) hSgn1->Write();
  if(hSgn2) hSgn2->Write();
  if(hSgn3) hSgn3->Write();

  if(hW3JetsLooseTauIso) hW3JetsLooseTauIso->Write();
  if(hW3JetsLooseTauIsoFR) hW3JetsLooseTauIsoFR->Write();
  if(hDataAntiIsoLooseTauIso) hDataAntiIsoLooseTauIso->Write();
  if(hDataAntiIsoLooseTauIsoFR) hDataAntiIsoLooseTauIsoFR->Write();
  if(hDataAntiIsoLooseTauIsoFRUp) hDataAntiIsoLooseTauIsoFRUp->Write();
  if(hDataAntiIsoLooseTauIsoFRDown) hDataAntiIsoLooseTauIsoFRDown->Write();
  if(hDataAntiIso) hDataAntiIso->Write();
  if(hDataAntiIsoFR) hDataAntiIsoFR->Write();
  if(hDataAntiIsoFRUp) hDataAntiIsoFRUp->Write();
  if(hDataAntiIsoFRDown) hDataAntiIsoFRDown->Write();

  cout << "Starts writing signal masses" << endl;

  for(int im=0 ; im<nMasses ; im++) {
    if(hggH[im]) hggH[im]->Write();
    if(hqqH[im]) hqqH[im]->Write();
    if(hVH[im])  hVH[im]->Write();
  }

  cout << "Starts writing susy" << endl;

  for(unsigned int i = 0; i < SUSYhistos.size() ; i++){
    ((mapSUSYhistos.find( SUSYhistos[i] ))->second)->Write();
  }
  cout << "Starts writing other stuff" << endl;

  hData->Write();
  hParameters->Write();
  fout->Write();
  fout->Close();

  cout << "Starts deleting pointers [mwahaha]" << endl;

  delete hQCD; delete hZmm; delete hZmj; delete hZfakes; delete hTTb; delete hZtt; delete hW; delete hW3Jets; delete hLooseIso; delete hAntiIso;
  delete hZmmLoose; delete hZmjLoose;
  delete hWLooseIso2; delete hLooseIso2; delete hLooseIso3;
  delete hW3JetsLooseTauIso; delete hW3JetsLooseTauIsoFR;
  delete hDataAntiIsoLooseTauIso; delete hDataAntiIsoLooseTauIsoFR; delete hDataAntiIsoLooseTauIsoFRUp; delete hDataAntiIsoLooseTauIsoFRDown;
  delete hDataAntiIso; delete hDataAntiIsoFR; delete hDataAntiIsoFRUp; delete hDataAntiIsoFRDown;

  if(hW3JetsKeys)   delete hW3JetsKeys;
  if(hWKeys)        delete hWKeys;
  if(hZmmKeys)      delete hZmmKeys;
  if(hZmjKeys)      delete hZmjKeys;
  if(hVVKeys)       delete hVVKeys;
  if(hLooseIsoKeys) delete hLooseIsoKeys;
  if(hAntiIsoKeys)  delete hAntiIsoKeys;
  //if(hHelp)         delete hHelp;
  
  for(int im=0 ; im<nMasses ; im++) {
    if(hggH[im]) delete hggH[im];
    if(hqqH[im]) delete hqqH[im];
    if(hVH[im])  delete hVH[im];
  }

  cout << "Starts deleting Susy pointers" << endl;

  for(unsigned int i = 0; i < SUSYhistos.size() ; i++) delete mapSUSYhistos.find( SUSYhistos[i] )->second ;
  for(unsigned int i = 0; i < SUSYhistos.size() ; i++){
    (mapSUSYfiles.find( SUSYhistos[i] )->second)->Close();
    delete mapSUSYfiles.find( SUSYhistos[i] )->second ;
  }

  cout << "Starts deleting some pointers" << endl;  

  delete hVV; delete hSgn1; delete hSgn2; delete hSgn3; delete hData; 

  cout << "Delete parameters" << endl;
  //delete hParameters;

  cout << "Delete stuff" << endl;
  delete hWMt; delete aStack;  delete hEWK; delete hSiml; delete hDataEmb; delete hSgn; delete hRatio; delete line;

  cout << "delete fout" << endl;
  delete fout;

  cout << "Starts deleting signal pointers" << endl;  

  fSignalGGH->Close();       delete fSignalGGH; 
  fSignalVBF->Close();       delete fSignalVBF;
  fSignalVH->Close();        delete fSignalVH; 

  cout << "Starts deleting masses pointers" << endl;

  for(int im=0 ; im<nMasses ; im++) {
    if(fSignalGGH_m[im]) { fSignalGGH_m[im]->Close(); /*delete fSignalGGH_m[im];*/ }
    // if(fSignalqqH_m[im]) { fSignalqqH_m[im]->Close();    delete fSignalqqH_m[im]; }
    // if(fSignalVH_m[im])   { fSignalVH_m[im]->Close();    delete fSignalVH_m[im];  }
  }

  cout << "Starts deleting other pointers" << endl;

  if(fBackgroundOthers) { fBackgroundOthers->Close();delete fBackgroundOthers; }
  if(fBackgroundTTbar)  { fBackgroundTTbar->Close(); delete fBackgroundTTbar;  }
  if(fBackgroundWJets)  { fBackgroundWJets->Close(); delete fBackgroundWJets; }
  if(fData)             { fData->Close();            delete fData; }
  if(dummy1)            { dummy1->Close();           delete dummy1;}
  if(fBackgroundDY)     { fBackgroundDY->Close();    delete fBackgroundDY;}

}

///////////////////////////////////////////////////////////////////////////////////////////////

void plotMuTauAll( Int_t useEmbedded = 1, 
		   TString outputDir = "/home/llr/cms/ndaci/WorkArea/HTauTau/AnalysisLorenzo/CMSSW_5_2_6/src/Bianchi/Limits/bin/nad/results/",
		   TString dirIn     = "/data_CMS/cms/ndaci/ndaci_2012/HTauTau/Analysis/MuTau/ntuples/"
		   )
{

  const int nMasses=1;
  const int nVar=1;
  const int nCat=10;
  const int nCatMode=3;
  const int nAnalyses=3;
  const int nAnMode=2;

  TString mH[nMasses]     = {"125"};

  TString variables[nVar] = {"diTauNSVfitMass"};
  TString xTitle[nVar]    = {"mass"};
  TString unit[nVar]      = {"GeV"};

  TString categ[nCat]        = {"inclusive","oneJet","twoJets","vbf","vh","0jet","boost","boost2","bTag","nobTag"};
  TString cat_mode[nCatMode] = {"","Low","High"};

  TString analysis[nAnalyses] = {"","Tau","Jet"};
  TString an_mode[nAnMode]    = {"Up","Down"};

  TString category="";
  TString analy="";

  const int nDo = 1;
  int idx_cat[nDo] = {0};
  //int idx_cat[nDo] = {0,3,5,6,8};
  //int idx_cat[2] = {3,5};

  int nBins; float min, max;

  for(int iVar=0 ; iVar<nVar; iVar++)
    for(int iCat=0 ; iCat<nDo; iCat++)
      for(int iMode=0 ; iMode<nCatMode ; iMode++)
	for(int iAN=0 ; iAN<1 ; iAN++)
	  for(int iam=0 ; iam<nAnMode ; iam++)
	    for(int iM=0 ; iM<nMasses ; iM++) {
	      
	      category = categ[ idx_cat[iCat] ];
	      if(category!="inclusive" && category!="vbf") category += cat_mode[iMode] ;
	      
	      analy = analysis[iAN];
	      if(analy!="") analy += an_mode[iM];

	      if(categ[ idx_cat[iCat] ]=="inclusive") {
		nBins=35; min=0; max=350;
	      }
	      else {
		nBins=-1; min=0; max=100;
	      }

	      plotMuTau( mH[iM] , useEmbedded , category , analy , variables[iVar] , 
			 xTitle[iVar], unit[iVar] , outputDir, dirIn, 
			 nBins , min , max , 1.0 , 1.0 , 0 , 1.2 , 5.072663 , "5.07" , 0.976 );
	    }
  /*
void plotMuTau( TString mH_         = "120",
		Int_t   useEmbedding_ = 0,
	        TString selection_  = "inclusiveHigh",
		TString analysis_   = "",		  
		TString variable_   = "diTauVisMass",
		TString XTitle_     = "full mass",
		TString Unities_    = "GeV",
		TString dirOut      = "/home/llr/cms/ndaci/WorkArea/HTauTau/AnalysisLorenzo/CMSSW_5_2_6/src/Bianchi/Limits/bin/nad/results/",
		TString dirIn       = "/data_CMS/cms/ndaci/ndaci_2012/HTauTau/Analysis/MuTau/ntuples/", 
		Int_t   nBins_ = 80, Float_t xMin_=0, Float_t xMax_=400,
		Float_t magnifySgn_ = 1.0,
		Float_t hltEff_     = 1.0,
		Int_t   logy_       = 0,
		Float_t maxY_       = 1.2,
		Float_t Lumi        = 5.072663, // fb-1
		TString tsLumi      = "5.07",
		Float_t lumiCorrFactor = 0.944 // correct DY th cross-section
		) 
  */
     
}

int main(int argc, const char* argv[])
{

  std::cout << "plotMuTau()" << std::endl;
  gROOT->SetBatch(true);
 

  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  plotMuTauAll();

}

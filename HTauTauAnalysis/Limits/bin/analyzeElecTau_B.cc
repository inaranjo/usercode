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

void plotElecTau( Int_t mH_           = 120,
		  Int_t useEmbedding_ = 0,
		  string selection_   = "inclusive",
		  string analysis_    = "",		  
		  TString variable_   = "diTauVisMass",
		  TString XTitle_     = "full mass",
		  TString Unities_    = "GeV",
		  TString outputDir   = "./",
		  Int_t nBins_ = 80, Float_t xMin_=0, Float_t xMax_=400,
		  Float_t magnifySgn_ = 1.0,
		  Float_t hltEff_     = 1.0,
		  Int_t logy_         = 0,
		  Float_t maxY_       = 1.2
		  ) 
{   

  ofstream out(Form("/home/llr/cms/ivo/CMSSW_4_2_8_patch7/src/Bianchi/Limits/htautau/histograms/%s/yieldsElecTau_mH%d_%s_%s.txt",outputDir.Data(),mH_,selection_.c_str(), analysis_.c_str() ),ios_base::out); 
  out.precision(5);
  out<< " => " << selection_ << endl;

  // input txt file with bins
  ifstream is;

  char* c = new char[10];
  is.open(Form("/home/llr/cms/ivo/CMSSW_4_2_8_patch7/src/Bianchi/Limits/htautau/bins/bins_eTau_%s_%s.txt",variable_.Data(), selection_.c_str())); 
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
  is.open(Form("/home/llr/cms/ivo/CMSSW_4_2_8_patch7/src/Bianchi/Limits/htautau/bins/bins_eTau_%s_%s.txt",variable_.Data(), selection_.c_str())); 
  
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

  // Luminosity analyzed & parameters
  //float Lumi = (-47.4 + 215.3 + 930.7 + 410.6 + (450.6+212.7) + (735.2+254.8+778.2+682.0) )*1.00;//Lorenzo
  float Lumi = (658.886 + 413.116 + 216.294 + 930.211 + 750.188 + 284.550 + 806.800 + 300.874 + 368.488 )*1.00;//Ivo
  //float Lumi = ( 519.312 + 254.164 + 1009 + 496.835 + 2223)*1.00;//embedded
  

  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  float OStoSSRatioQCD            = 1.07;
  float SSIsoToSSAIsoRatioQCD     = 1.954;//0.32;
  float EtoTauCorrectionFactor    = 1.00;
  float JtoTauCorrectionFactor    = 1.00;
  float embeddedMEtCutEff         = 1.00;
  float madgraphMEtCutEff         = 1.00;

  /*
  // Fall11_16Nov2011
  float WcorrectionFactorOS       = 0.96;  
  float WcorrectionFactorSS       = 1.20; 
  float NoVbfExtrapolationFactorZ = 0.993;
  float VbfExtrapolationFactorZ   = 1.004;
  float BoostExtrapolationFactorZ = 0.94;
  float VbfExtrapolationFactorW   = 1.004;
  float BoostExtrapolationFactorW = 0.94;
  */

  // Fall11_06Dec2011
  float WcorrectionFactorOS       = 1.13;  
  float WcorrectionFactorSS       = 1.59; 
  float NoVbfExtrapolationFactorZ = 0.996;
  float VbfExtrapolationFactorZ   = 1.23;
  float BoostExtrapolationFactorZ = 1.11;
  float VbfExtrapolationFactorW   = 1.23;
  float BoostExtrapolationFactorW = 1.11;

  //float NoVbfExtrapolationFactorSSAntiIso = 0.918123;
  //float VbfExtrapolationFactorSSAntiIso   = 0.0013;
  //float BoostExtrapolationFactorSSAntiIso = 0.0050;

  bool useMt      = true;
  string antiWcut = useMt ? "MtLeg1Corr" : "-(pZetaCorr-1.5*pZetaVisCorr)" ;
  float antiWsgn  = useMt ? 40. :  20. ; 
  float antiWsdb  = useMt ? 60. :  40. ; 

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
  leg->SetHeader(Form("#splitline{CMS Preliminary #sqrt{s}=7 TeV}{%.1f fb^{-1} #tau_{e}#tau_{had}}", Lumi/1000. ));

  THStack* aStack = new THStack("aStack","");

  TH1F* hSiml    = new TH1F( "hSiml"   ,"all"               , nBins , bins.GetArray());
  TH1F* hSgn     = new TH1F( "hSgn "   ,"vbf+ggf+vh"        , nBins , bins.GetArray());
  TH1F* hSgn1    = new TH1F( "hSgn1"   ,"vbf"               , nBins , bins.GetArray());
  TH1F* hSgn2    = new TH1F( "hSgn2"   ,"ggf"               , nBins , bins.GetArray());
  TH1F* hSgn3    = new TH1F( "hSgn3"   ,"vh"                , nBins , bins.GetArray());
  TH1F* hData    = new TH1F( "hData"   ,"        "          , nBins , bins.GetArray());
  TH1F* hDataEmb = new TH1F( "hDataEmb","Embedded"          , nBins , bins.GetArray());
  TH1F* hW       = new TH1F( "hW"      ,"W+jets"            , nBins , bins.GetArray());
  TH1F* hW3Jets  = new TH1F( "hW3Jets" ,"W+3jets"           , nBins , bins.GetArray());
  TH1F* hEWK     = new TH1F( "hEWK"    ,"EWK"               , nBins , bins.GetArray());
  TH1F* hZtt     = new TH1F( "hZtt"    ,"Ztautau"           , nBins , bins.GetArray());
  TH1F* hZmm     = new TH1F( "hZmm"    ,"Z+jets, e to tau"  , nBins , bins.GetArray());
  TH1F* hZmj     = new TH1F( "hZmj"    ,"Z+jets, jet to tau", nBins , bins.GetArray());
  TH1F* hZfakes  = new TH1F( "hZfakes" ,"Z+jets, jet to tau", nBins , bins.GetArray());
  TH1F* hTTb     = new TH1F( "hTTb"    ,"ttbar"             , nBins , bins.GetArray());
  TH1F* hQCD     = new TH1F( "hQCD"    ,"QCD"               , nBins , bins.GetArray());
  TH1F* hLooseIso= new TH1F( "hLooseIso","Loose Iso"        , nBins , bins.GetArray());
  TH1F* hAntiIso = new TH1F( "hAntiIso","Anti Iso"          , nBins , bins.GetArray());
  TH1F* hVV      = new TH1F( "hVV"     ,"Diboson"           , nBins , bins.GetArray());

  TH1*  hW3JetsKeys   = 0;
  TH1*  hWKeys        = 0;
  TH1*  hLooseIsoKeys = 0;
  TH1*  hAntiIsoKeys  = 0;
  TH1*  hZmmKeys      = 0;
  TH1*  hZmjKeys      = 0;
  TH1*  hVVKeys       = 0;

 
  // pZeta OS, N pZ sideband OS, pZeta SS, N sideband SS, N QCD SS, OS/SS
  TH1F* hParameters   = new TH1F( "hParameters", "" , 7, 0,7);

  // Open the files

  TFile *fData              
    = new TFile("/data_CMS/cms/ivo/Trees/ElecTauStream_16Feb2012/OpenNtuples//nTupleRun2011-ElecTau-All_run_Open_ElecTauStream.root", "READ");  
  TFile *fDataLooseIso              
    = new TFile("/data_CMS/cms/ivo/Trees/ElecTauStream_16Feb2012/OpenNtuples//nTupleRun2011-ElecTau-All_run_Open_ElecTauStream.root", "READ");  
  TFile *fDataEmbedded              
    = new TFile("/data_CMS/cms/ivo/Trees/ElecTauStream_16Feb2012/OpenNtuples//nTupleRun2011-ElecTau-Embedded-All_run_Open_ElecTauStream.root", "READ");  
  TFile *fSignalVBF         
    = new TFile(Form("/data_CMS/cms/ivo/Trees/ElecTauStream_16Feb2012/OpenNtuples//nTupleVBFH%d-ElecTau-powheg-PUS6_run_Open_ElecTauStream.root",mH_) ,"READ");  
  TFile *fSignalGGH         
    = new TFile(Form("/data_CMS/cms/ivo/Trees/ElecTauStream_16Feb2012/OpenNtuples//nTupleGGFH%d-ElecTau-powheg-PUS6_run_Open_ElecTauStream.root",mH_),"READ"); 
  TFile *fSignalVH         
    = new TFile(Form("/data_CMS/cms/ivo/Trees/ElecTauStream_16Feb2012/OpenNtuples//nTupleVH%d-ElecTau-pythia-PUS6_run_Open_ElecTauStream.root",mH_),"READ");  
  TFile *fBackgroundDY
    = new TFile("/data_CMS/cms/ivo/Trees/ElecTauStream_16Feb2012/OpenNtuples//nTupleDYJets-ElecTau-50-madgraph-PUS6-v2_run_Open_ElecTauStream.root","READ"); 
  TFile *fBackgroundWJets   
    = new TFile("/data_CMS/cms/ivo/Trees/ElecTauStream_16Feb2012/OpenNtuples//nTupleWJets-ElecTau-madgraph-PUS6_run_Open_ElecTauStream.root","READ"); 
  TFile *fBackgroundW3Jets   
    = new TFile("/data_CMS/cms/ivo/Trees/ElecTauStream_16Feb2012/OpenNtuples//nTupleW3Jets-ElecTau-madgraph-PUS6_run_Open_ElecTauStream.root","READ"); 
  TFile *fBackgroundTTbar  
    = new TFile("/data_CMS/cms/ivo/Trees/ElecTauStream_16Feb2012/OpenNtuples//nTupleTTJets-ElecTau-madgraph-PUS6_run_Open_ElecTauStream.root","READ"); 
  TFile *fBackgroundOthers  
  = new TFile("/data_CMS/cms/ivo/Trees/ElecTauStream_16Feb2012/OpenNtuples//nTupleOthers-ElecTau-PUS6_run_Open_ElecTauStream.root","READ"); 
  
  // choose the analysis: Nominal "", jet up/Down "JetUp/Down" , elec up/down "ElecUp/Down" , tau up/down "TauUp/Down"
  TString tree         = "outTreePtOrd"+analysis_;
  TString treeEmbedded = "outTreePtOrd";
  if(analysis_.find("TauUp")  !=string::npos) 
    treeEmbedded = tree;
  if(analysis_.find("TauDown")!=string::npos) 
    treeEmbedded = tree;
  if(analysis_.find("ElecUp")  !=string::npos) 
    treeEmbedded = tree;
  if(analysis_.find("ElecDown")!=string::npos) 
    treeEmbedded = tree;

  TTree *data                = (TTree*)fData->Get("outTreePtOrd");
  TTree *dataLooseIso        = LOOSEISO  ? (TTree*)fDataLooseIso->Get("outTreePtOrd") : 0;
  TTree *dataEmbedded        = EMBEDDEDSAMPLES ? (TTree*)fDataEmbedded->Get(treeEmbedded) : 0;
  TTree *signalVBF           = (TTree*)fSignalVBF->Get(tree);
  TTree *signalGGH           = (TTree*)fSignalGGH->Get(tree);
  TTree *signalVH            = addVH ? (TTree*)fSignalVH->Get(tree) : 0;

  // split the DY->ll into l=e/mu and l=tau (MC level) ===> temporary, need fix !!!
  TFile* dummy1 = new TFile("dummy1.root","RECREATE");
  cout << "Now copying g/Z -> tau+ tau- " << endl;
  TTree *backgroundDYTauTau  = ((TTree*)fBackgroundDY->Get(tree))->CopyTree("abs(genDecay)==(23*15)");                 // g/Z -> tau+ tau-
  cout << "Now copying g/Z -> e+e- e->tau" << endl;
  TTree *backgroundDYEtoTau = ((TTree*)fBackgroundDY->Get(tree))->CopyTree("abs(genDecay)!=(23*15) &&  leptFakeTau"); // g/Z -> mu+mu- mu->tau
  cout << "Now copying g/Z -> e+e- jet->tau" << endl;
  TTree *backgroundDYJtoTau  = ((TTree*)fBackgroundDY->Get(tree))->CopyTree("abs(genDecay)!=(23*15) && !leptFakeTau"); // g/Z -> mu+mu- jet->tau

  cout << backgroundDYTauTau->GetEntries()  << " come from DY->tautau"         << endl;
  cout << backgroundDYEtoTau->GetEntries() << " come from DY->ee, e->tau"  << endl;
  cout << backgroundDYJtoTau->GetEntries()  << " come from DY->ee, jet->tau" << endl;

  TTree *backgroundTTbar     = (TTree*)fBackgroundTTbar->Get(tree);
  TTree *backgroundWJets     = (TTree*)fBackgroundWJets->Get(tree);
  TTree *backgroundW3Jets    = W3JETS ? (TTree*)fBackgroundW3Jets->Get(tree) : 0;
  TTree *backgroundOthers    = (TTree*)fBackgroundOthers->Get(tree);
 
  //MIT MVAWP
  //TCut lpt("ptL1>20 && TMath::Abs(etaL1)<2.1 && tightestMVAWP>=1 && antiConv>0.5");

  //   //Daniele MVA for Electron ID for working points in barrel and endcap first WP
  //   TCut lpt("(ptL1>20 && ptL1<30 && TMath::Abs(etaL1)>1.5 && TMath::Abs(etaL1)<2.1  && tightestDanieleMVAWP>=0.0275)||(ptL1>20 && ptL1<30 && TMath::Abs(etaL1)<1.5  && tightestDanieleMVAWP>=0.035)||(ptL1>30 && TMath::Abs(etaL1)>1.5 && TMath::Abs(etaL1)<2.1 && tightestDanieleMVAWP>=0.0225)||(ptL1>30 && TMath::Abs(etaL1)<1.5  && tightestDanieleMVAWP>=0.035)");
  //Daniele MVA for Electron ID for working points in barrel and endcap second WP
  //TCut lpt("antiConv>0.5 && (ptL1>20 && TMath::Abs(etaL1)>1.5 && TMath::Abs(etaL1)<2.1  && tightestDanieleMVAWP>=0.0225)||(ptL1>20 && TMath::Abs(etaL1)<1.5  && tightestDanieleMVAWP>=0.01625)");
  //Daniele MVA for Electron ID for working points in barrel and endcap third WP
  TCut presel("(TMath::Abs(etaL1)<1.5 && sihih<0.01 && dEta<0.007 && dPhi<0.15 && HoE<0.12) || (TMath::Abs(etaL1)>1.5 && TMath::Abs(etaL1)<2.1 && sihih<0.03 && dEta<0.009 && dPhi<0.10 && HoE<0.10)");
  TCut lpt("nHits<=0 && antiConv>0.5 && (ptL1>20 && TMath::Abs(etaL1)<1.0  && tightestDanieleMVAWP>=0.013)||(ptL1>20 && TMath::Abs(etaL1)>1.0 && TMath::Abs(etaL1)<1.5  && tightestDanieleMVAWP>=0.0425)|| (ptL1>20 && TMath::Abs(etaL1)>1.5 && TMath::Abs(etaL1)<2.1  && tightestDanieleMVAWP>=0.0225)");
  lpt = lpt && presel;
  
  //Cut based WP
  //TCut lpt("ptL1>20 && TMath::Abs(etaL1)<2.1 && tightestCutBasedWP>=1 && antiConv>0.5");
  
  //NO ID
  //TCut lpt("ptL1>20 && TMath::Abs(etaL1)<2.1 && antiConv>0.5");

  TCut tpt("ptL2>20 && TMath::Abs(etaL2)<2.3 && tightestAntiEMVAWP>0.5");//-1all pass or 0.5
  TCut tiso("tightestHPSDBWP>0");
  TCut laiso("combRelIsoLeg1DBeta>0.30 && combRelIsoLeg1DBeta<0.50");
  TCut lliso("combRelIsoLeg1DBeta<0.30");
  TCut liso("combRelIsoLeg1DBeta<0.10");
  TCut lveto("elecFlag==0");
  TCut SS("diTauCharge!=0");
  TCut OS("diTauCharge==0");
  TCut oneJet("nJets30>=1");
  TCut twoJets("nJets30>=2");
  TCut vbf("pt1>30 && pt2>30 && eta1*eta2<0 && Mjj>400 && Deta>4.0 && isVetoInJets!=1");
  TCut novbf("pt1<150 && pt2<30");
  TCut boost("pt1>150 && !(pt2>30 && eta1*eta2<0 && Mjj>400 && Deta>4.0 && isVetoInJets!=1)");
  TCut bTag("nJets30<=1 && nJets20BTagged>=1");
  TCut nobTag("nJets30<=1 && nJets20BTagged==0");
  TCut hltevent("HLTx==1 && (run>=163269 || run==1)");
  TCut hltmatch("HLTmatch==1");
  TCut pZ( Form("((%s)<%f)",antiWcut.c_str(),antiWsgn));
  TCut apZ(Form("((%s)>%f)",antiWcut.c_str(),antiWsdb));

  TCut sbin; TCut sbinEmbedding; TCut sbinEmbeddingPZetaRel; TCut sbinPZetaRel; TCut sbinSS; TCut sbinPZetaRelSS; TCut sbinPZetaRev; TCut sbinPZetaRevSS; TCut sbinSSaIso; TCut sbinSSlIso;

  if(selection_.find("inclusive")!=string::npos){
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
  }
  else if(selection_.find("oneJet")!=string::npos){
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
  }
  else if(selection_.find("twoJets")!=string::npos){
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
  }
  else if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos){
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
  }
  else if(selection_.find("novbf")!=string::npos){
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
  }
  else if(selection_.find("boost")!=string::npos ){
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
  }
  else if(selection_.find("bTag")!=string::npos && selection_.find("nobTag")==string::npos){
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
  }
  else if(selection_.find("nobTag")!=string::npos){
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
  }

  // estimate the W+jets in the selection bin using pZeta extrapolation

  TH1F* hWMt = new TH1F("hWMt","",1,-10,10);

  ///////////////////////////////////////// Doing OS first...
  hWMt->Reset();
  backgroundWJets->Draw("etaL1>>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRel&&pZ));
  cout << "Using  " << hWMt->GetEntries() << " entries from the W+jets OS sample" << endl;
  float OSWinSignalRegionMC   = hWMt->Integral()*Lumi*hltEff_/1000.;
  hWMt->Reset();
  backgroundWJets->Draw("etaL1>>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRel&&apZ));
  float OSWinSidebandRegionMC = hWMt->Integral()*Lumi*hltEff_/1000.;
  float scaleFactorOS = OSWinSignalRegionMC>0 ? OSWinSidebandRegionMC/OSWinSignalRegionMC : 1.0 ;

  if(useMt)
    cout << "Extrapolation factor for W OS : P(MtCorr>" << antiWsdb << ")/P(MtCorr<" << antiWsgn << ") ==> " << scaleFactorOS << endl;
  else
    cout << "Extrapolation factor for W OS : P(pZetaCorr<- "<< antiWsdb << ")/P(pZetaCorr>"<< antiWsgn << ") ==> " << scaleFactorOS << endl;    
 
  hWMt->Reset();
  cout << "Estimating cobtribution from Ztt, ttbar and others in OS low pZeta tail..." << endl;
  backgroundTTbar->Draw("etaL1>>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRel&&apZ));
  float ttbarExtrOS  = hWMt->Integral()*Lumi*hltEff_/1000.;
  cout << "Contribution from ttbar in OS is " << ttbarExtrOS << endl;
  hWMt->Reset();
  backgroundOthers->Draw("etaL1>>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRel&&apZ));
  float othersExtrOS = hWMt->Integral()*Lumi*hltEff_/1000.;
  cout << "Contribution from single-t and di-boson in OS is " << othersExtrOS << endl;
  hWMt->Reset();
  backgroundDYTauTau->Draw("etaL1>>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRel&&apZ));
  float dytautauExtrOS = hWMt->Integral()*Lumi*hltEff_/1000.;
  cout << "Contribution from DY->tautau in OS is " << dytautauExtrOS << endl;
  hWMt->Reset();
  backgroundDYJtoTau->Draw("etaL1>>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRel&&apZ));
  float dyjtotauExtrOS = hWMt->Integral()*Lumi*hltEff_/1000.;
  cout << "Contribution from DY->ee, jet->tau in OS is " << dyjtotauExtrOS << endl;
  hWMt->Reset();
  backgroundDYEtoTau->Draw("etaL1>>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRel&&apZ));
  float dyetotauExtrOS = hWMt->Integral()*Lumi*hltEff_/1000.;
  cout << "Contribution from DY->ee, e->tau in OS is " << dyetotauExtrOS << endl;

  float OSWinSignalRegionDATA = data->GetEntries(sbinPZetaRev);
  cout << "Selected events in data in low pZeta/low Mt tail " << OSWinSignalRegionDATA << endl;
  OSWinSignalRegionDATA -= ttbarExtrOS;
  OSWinSignalRegionDATA -= othersExtrOS;
  OSWinSignalRegionDATA -= dytautauExtrOS;
  OSWinSignalRegionDATA -= dyjtotauExtrOS;
  OSWinSignalRegionDATA -= dyetotauExtrOS;
  OSWinSignalRegionDATA /= scaleFactorOS;
  cout << "W+jets in signal region is estimated to be "  
       << OSWinSignalRegionDATA*scaleFactorOS << "/" << scaleFactorOS << " = " 
       << OSWinSignalRegionDATA <<  " +/- " << sqrt(OSWinSignalRegionDATA/scaleFactorOS)/scaleFactorOS << endl;
  cout << "  ===> the MC prediction was " << OSWinSignalRegionMC << endl;

  hParameters->SetBinContent(1, 1./scaleFactorOS );
  hParameters->SetBinContent(2, OSWinSignalRegionDATA*scaleFactorOS );

  ///////////////////////////////////////// Doing SS last...
  hWMt->Reset();
  backgroundWJets->Draw("etaL1>>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRelSS&&pZ));
  cout << "Using  " << hWMt->GetEntries() << " entries from the SS W+jets sample" << endl;
  float SSWinSignalRegionMC   = hWMt->Integral()*Lumi*hltEff_/1000.;
  hWMt->Reset();
  backgroundWJets->Draw("etaL1>>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRelSS&&apZ));
  float SSWinSidebandRegionMC = hWMt->Integral()*Lumi*hltEff_/1000.;
  float scaleFactorSS = SSWinSignalRegionMC>0 ? SSWinSidebandRegionMC/SSWinSignalRegionMC : 1.0;
 
  if(useMt)
    cout << "Extrapolation factor for W SS : P(MtCorr>" << antiWsdb << ")/P(MtCorr<" << antiWsgn << ") ==> " << scaleFactorSS << endl;
  else
    cout << "Extrapolation factor for W SS : P(pZetaCorr<- "<< antiWsdb << ")/P(pZetaCorr>"<< antiWsgn << ") ==> " << scaleFactorSS << endl;    

  hWMt->Reset();
  cout << "Estimating cobtribution Ztt,from ttbar and others in SS low pZeta tail..." << endl;
  backgroundTTbar->Draw("etaL1>>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRelSS&&apZ));
  float ttbarExtrSS = hWMt->Integral()*Lumi*hltEff_/1000.;
  cout << "Contribution from ttbar in SS is " << ttbarExtrSS << endl;
  hWMt->Reset();
  backgroundOthers->Draw("etaL1>>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRelSS&&apZ));
  float othersExtrSS = hWMt->Integral()*Lumi*hltEff_/1000.;
  cout << "Contribution from single-t and di-boson in SS is " << othersExtrSS << endl;
  hWMt->Reset();
  backgroundDYJtoTau->Draw("etaL1>>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRelSS&&apZ));
  float dyjtotauExtrSS = hWMt->Integral()*Lumi*hltEff_/1000.;
  cout << "Contribution from DY->ee, jet->tau in SS is " << dyjtotauExtrSS << endl;
  hWMt->Reset();
  backgroundDYEtoTau->Draw("etaL1>>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFElec*SFTau*SFEtoTau)"*(sbinPZetaRelSS&&apZ));
  float dyetotauExtrSS = hWMt->Integral()*Lumi*hltEff_/1000.;
  cout << "Contribution from DY->ee, e->tau in SS is " << dyetotauExtrSS << endl;
  
  float SSWinSignalRegionDATA = data->GetEntries(sbinPZetaRevSS);
  cout << "Selected events in data in low pZeta/low Mt tail " << SSWinSignalRegionDATA << endl;
  SSWinSignalRegionDATA -= ttbarExtrSS;
  SSWinSignalRegionDATA -= othersExtrSS;
  SSWinSignalRegionDATA -= dyjtotauExtrSS;
  SSWinSignalRegionDATA -= dyetotauExtrSS;
  SSWinSignalRegionDATA /= scaleFactorSS;
  cout << "W+jets in SS signal region is estimated to be "  
       << SSWinSignalRegionDATA*scaleFactorSS << "/" << scaleFactorSS << " = " 
       << SSWinSignalRegionDATA <<  " +/- " << sqrt(SSWinSignalRegionDATA/scaleFactorSS)/scaleFactorSS << endl;
  cout << "  ===> the MC prediction was " << SSWinSignalRegionMC << endl;

  hParameters->SetBinContent(3, 1./scaleFactorSS );
  hParameters->SetBinContent(4, SSWinSignalRegionDATA*scaleFactorSS );


  // here I choose the order in the stack
  std::vector<string> samples;
  samples.push_back("Data");
  if(dataLooseIso){
  samples.push_back("LooseIso");
  samples.push_back("AntiIso");
  }
  samples.push_back("ggH115");
  samples.push_back("qqH115");
  if(signalVH)
    samples.push_back("VH115");
  samples.push_back("Others");
  samples.push_back("TTbar");
  samples.push_back("SS");
  samples.push_back("WJets");
  if(backgroundW3Jets)
    samples.push_back("W3Jets");
  samples.push_back("DYEtoTau");
  samples.push_back("DYJtoTau");
  samples.push_back("DYToTauTau");
  if(dataEmbedded)
    samples.push_back("Embedded");
  

  // here I define the map between a sample name and its tree
  std::map<std::string,TTree*> tMap;
  tMap["Data"]         = data;
  tMap["LooseIso"]     = dataLooseIso;
  tMap["AntiIso"]      = dataLooseIso;
  tMap["Embedded"]     = dataEmbedded;
  tMap["ggH115"]       = signalGGH;
  tMap["qqH115"]       = signalVBF;
  tMap["VH115"]        = signalVH;
  tMap["DYToTauTau"]   = backgroundDYTauTau;
  tMap["DYEtoTau"]    = backgroundDYEtoTau;
  tMap["DYJtoTau"]     = backgroundDYJtoTau;
  tMap["WJets"]        = backgroundWJets;
  tMap["W3Jets"]       = backgroundW3Jets;
  tMap["Others"]       = backgroundOthers;
  tMap["TTbar"]        = backgroundTTbar;
  tMap["SS"]           = data;


  
  std::map<TString,Float_t> vMap;


  for( unsigned iter=0; iter<samples.size(); iter++){

    cout << "Dealing with sample " << samples[iter] << endl;
    
    std::map<std::string,TTree*>::iterator it = tMap.find(samples[iter]);

    TString h1Name = "h1_"+it->first;
    TH1F* h1 = new TH1F( h1Name ,"" , nBins , bins.GetArray());

    TTree* currentTree = 0;
    
    if((it->first).find("SS")!=string::npos){
      
      cout << "Remove W contamination from SS data sample ... " << endl;
      currentTree = (it->second);

      float error2OnQCD = 0.0;
      
      TH1F* hHelp = (TH1F*)h1->Clone("hHelp");
      hHelp->Reset();
      currentTree->Draw(variable_+">>hHelp", sbinSS);
      int SSevents = hHelp->GetEntries();
      cout << "Selected SS events in data " << hHelp->GetEntries() << endl;
      h1->Add(hHelp,1);

      hHelp->Reset();
      backgroundWJets->Draw(variable_+">>hHelp", "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*sbinSS);
      cout << "We expect " << hHelp->Integral()*Lumi/1000*hltEff_ << " SS events from W+jets (from " << hHelp->GetEntries() << " entries)" << endl;
      float sFWSS = ( selection_.find("novbf")!=string::npos || selection_.find("nobTag")!=string::npos || selection_.find("inclusive")!=string::npos) ? SSWinSignalRegionDATA/SSWinSignalRegionMC : WcorrectionFactorSS; // from the extrapolation factor DATA/MC
      if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos) sFWSS *= VbfExtrapolationFactorW;
      else if(selection_.find("boost")!=string::npos)  sFWSS *= BoostExtrapolationFactorW;
      hHelp->Scale(sFWSS*Lumi/1000*hltEff_);
      cout << "We estimate " << hHelp->Integral() << " SS events from W+jets by extrapolating" << endl;
      cout << " ==> removing W+jets from SS...." << endl;
      h1->Add(hHelp, -1 );
      if(hHelp->GetEntries()>0) error2OnQCD +=  pow( hHelp->Integral()/hHelp->GetEntries(), 2)*hHelp->GetEntries(); // error on MC W+jets SS events
      error2OnQCD +=  pow(WcorrectionFactorSS*0.06,2)*pow(hHelp->GetEntries(),2);        // error on W+jets extrapolation factor ==> 6% according to Artur
      cout << sqrt(error2OnQCD) << " <==  W" << endl;      

      hHelp->Reset();
      backgroundTTbar->Draw(variable_+">>hHelp", "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*sbinSS);
      cout << "We expect " << hHelp->Integral()*Lumi/1000*hltEff_ << " SS events from TTbar (from " << hHelp->GetEntries() << " entries)" << endl;
      hHelp->Scale(1.0*Lumi/1000*hltEff_);
      cout << "We estimate " << hHelp->Integral() << " SS events from TTbar" << endl;
      cout << " ==> removing TTbar from SS...." << endl;
      h1->Add(hHelp, -1 );
      if(hHelp->GetEntries()>0) error2OnQCD +=  pow(hHelp->Integral()/hHelp->GetEntries(),2)*hHelp->GetEntries();   // error on MC TTbar SS events
      cout << sqrt(error2OnQCD) << " <== W + TTb" << endl;      

      hHelp->Reset();
      backgroundDYEtoTau->Draw(variable_+">>hHelp", "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*sbinSS);
      cout << "We expect " << hHelp->Integral()*Lumi/1000*hltEff_ << " SS events from DY->ee, e->jet" << endl;
      hHelp->Scale(EtoTauCorrectionFactor*Lumi/1000*hltEff_);
      cout << "We estimate " << hHelp->Integral() << " SS events from DY->ee, e->tau" << endl;
      cout << " ==> removing DY->ee, e->tau from SS...." << endl;
      h1->Add(hHelp, -1 );
      if(hHelp->GetEntries()>0) error2OnQCD +=  pow(hHelp->Integral()/hHelp->GetEntries(),2)*hHelp->GetEntries();   // error on MC DY->mumu, mu->tau events
      cout << sqrt(error2OnQCD) << " <== W + TTb + DY(1)" << endl;      

      hHelp->Reset();
      backgroundDYJtoTau->Draw(variable_+">>hHelp", "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*sbinSS);
      cout << "We expect " << hHelp->Integral()*Lumi/1000*hltEff_ << " SS events from DY->ee, jet->tau" << endl;
      hHelp->Scale(JtoTauCorrectionFactor*Lumi/1000*hltEff_);
      cout << "We estimate " << hHelp->Integral() << " SS events from DY->ee, jet->tau" << endl;
      cout << " ==> removing DY->ee, e->jet from SS...." << endl;
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

    }
    else{

      currentTree = (it->second);

      if((it->first).find("Embed")==string::npos){

	if((it->first).find("DYToTauTau")!=string::npos){
	  currentTree->Draw(variable_+">>"+h1Name, "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFElec*SFTau)"*sbinPZetaRel);  
	  float madgraphNoMEtCut = h1->Integral();
	  h1->Reset();
	  currentTree->Draw(variable_+">>"+h1Name, "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFElec*SFTau)"*sbin);
	  madgraphMEtCutEff = h1->Integral()/madgraphNoMEtCut;
	  cout << "Efficiency of antiW cut on madgraph " << madgraphMEtCutEff << endl;
	}
	else if((it->first).find("W3Jets")!=string::npos){
	  currentTree->Draw(variable_+">>"+h1Name,     "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFElec*SFTau*HqTWeight)"*sbin);
	  hW3JetsKeys = new TH1Keys("hW3JetsKeys","W+3jets smoothed", nBins , bins.GetArray());
	  currentTree->Draw(variable_+">>hW3JetsKeys", "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFElec*SFTau*HqTWeight)"*sbin);
	  cout << "Keys for W3Jets filled with integral " << hW3JetsKeys->Integral() << " and entries " << hW3JetsKeys->GetEntries() << endl;
	}
	else if((it->first).find("WJets")!=string::npos){
	  currentTree->Draw(variable_+">>"+h1Name,     "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFElec*SFTau*HqTWeight)"*sbin);
	  hWKeys = new TH1Keys("hWKeys","W+jets smoothed", nBins , bins.GetArray());
	  currentTree->Draw(variable_+">>hWKeys", "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFElec*SFTau*HqTWeight)"*sbin);
	  cout << "Keys for WJets filled with integral " << hWKeys->Integral() << " and entries " << hWKeys->GetEntries() << endl;
	}
	else if((it->first).find("DYEtoTau")!=string::npos){
	  currentTree->Draw(variable_+">>"+h1Name,     "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFElec*SFTau*HqTWeight)"*sbin);
	  hZmmKeys = new TH1Keys("hZmmKeys","Z+jets, e to tau smoothed", nBins , bins.GetArray());
	  currentTree->Draw(variable_+">>hZmmKeys", "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFElec*SFTau*HqTWeight)"*sbin);
	  cout << "Keys for Zmm filled with integral " << hZmmKeys->Integral() << " and entries " << hZmmKeys->GetEntries() << endl;
	}
	else if((it->first).find("DYJtoTau")!=string::npos){
	  currentTree->Draw(variable_+">>"+h1Name,     "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFElec*SFTau*HqTWeight)"*sbin);
	  hZmjKeys = new TH1Keys("hZmjKeys","Z+jets, jet to tau smoothed", nBins , bins.GetArray());
	  currentTree->Draw(variable_+">>hZmjKeys", "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFElec*SFTau*HqTWeight)"*sbin);
	  cout << "Keys for Zmj filled with integral " << hZmjKeys->Integral() << " and entries " << hZmjKeys->GetEntries() << endl;
	}
	else if((it->first).find("Others")!=string::npos){
	  currentTree->Draw(variable_+">>"+h1Name,     "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFElec*SFTau*HqTWeight)"*sbin);
	  hVVKeys = new TH1Keys("hVVKeys","Others smoothed", nBins , bins.GetArray());
	  currentTree->Draw(variable_+">>hVVKeys", "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFElec*SFTau*HqTWeight)"*sbin);
	  cout << "Keys for VV filled with integral " << hVVKeys->Integral() << " and entries " << hVVKeys->GetEntries() << endl;
	}
	else if((it->first).find("LooseIso")!=string::npos){
	  currentTree->Draw(variable_+">>"+h1Name,    sbinSSlIso);
	  hLooseIsoKeys = new TH1Keys("hLooseIsoKeys","Loose Iso smoothed", nBins , bins.GetArray());
	  if(  ((selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos) || selection_.find("boost")!=string::npos) )
	    currentTree->Draw(variable_+">>hLooseIsoKeys", sbinSSlIso);
	  cout << "Keys for LooseIso filled with integral " << hLooseIsoKeys->Integral() << " and entries " << hLooseIsoKeys->GetEntries() << endl;
	}
	else if((it->first).find("AntiIso")!=string::npos){
	  currentTree->Draw(variable_+">>"+h1Name,    sbinSSaIso);
	  hAntiIsoKeys = new TH1Keys("hAntiIsoKeys","Anti Iso smoothed", nBins , bins.GetArray());
	  currentTree->Draw(variable_+">>hAntiIsoKeys", sbinSSaIso);
	  cout << "Keys for AntiIso filled with integral " << hAntiIsoKeys->Integral() << " and entries " << hAntiIsoKeys->GetEntries() << endl;
	}
	else
	  currentTree->Draw(variable_+">>"+h1Name, "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFElec*SFTau*HqTWeight)"*sbin);
      }
      else{
	currentTree->Draw(variable_+">>"+h1Name, "(HLTTau*HLTElec*embeddingWeight)"*sbinEmbeddingPZetaRel);
	float embeddedNoMEtCut = h1->Integral();
	h1->Reset();
	currentTree->Draw(variable_+">>"+h1Name, "(HLTTau*HLTElec*embeddingWeight)"*sbinEmbedding);

	embeddedMEtCutEff =  h1->Integral()/embeddedNoMEtCut;
	cout << "Efficiency of antiW cut on embedded " << embeddedMEtCutEff << endl;
      }


      // scale by correction factors
      // scale by correction factors
      if( ! ((it->first).find("Data")!=string::npos || 
	     (it->first).find("LooseIso")!=string::npos ||
	     (it->first).find("AntiIso")!=string::npos) ) 
	h1->Scale(Lumi/1000*hltEff_);

      // if W+jets, scale by extrapolation
      float sFWOS = ( selection_.find("novbf")!=string::npos || selection_.find("nobTag")!=string::npos || selection_.find("inclusive")!=string::npos) ? 
	OSWinSignalRegionDATA/OSWinSignalRegionMC : WcorrectionFactorOS;
      if((it->first).find("WJets")!=string::npos){
	if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos){
	  sFWOS *= VbfExtrapolationFactorW;
	  cout << "Wjets will be rescaled by " << VbfExtrapolationFactorW << " according to the Z->ee+j+vbf/Z->ee+j ratio" << endl;
	}
	else if(selection_.find("boost")!=string::npos){
	  sFWOS *= BoostExtrapolationFactorW;
	  cout << "Wjets will be rescaled by " << BoostExtrapolationFactorW << " according to the Z->ee+j+vbf/Z->ee+j ratio" << endl;
	}	else if(selection_.find("boost")!=string::npos){
	  sFWOS *= BoostExtrapolationFactorW;
	  cout << "Wjets will be rescaled by " << BoostExtrapolationFactorW << " according to the Z->ee+j+vbf/Z->ee+j ratio" << endl;
	}
	h1->Scale( sFWOS );
	hW->Add(h1,1.0);
      }

      // if DY->tautau, and vbf scale by ratio data/MC
      if((it->first).find("DYToTauTau")!=string::npos){

	if(selection_.find("novbf")!=string::npos){
	  cout << "DY->tautau will be rescaled by " << NoVbfExtrapolationFactorZ << " according to the Z->ee+vbf/Z->ee ratio" << endl;
	  h1->Scale( NoVbfExtrapolationFactorZ );
	}
	else if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos){
	  cout << "DY->tautau will be rescaled by " << VbfExtrapolationFactorZ << " according to the Z->ee+vbf/Z->ee ratio" << endl;
	  h1->Scale( VbfExtrapolationFactorZ );
	}
	else if(selection_.find("boost")!=string::npos){
	  cout << "DY->tautau will be rescaled by " << BoostExtrapolationFactorZ << " according to the Z->ee+vbf/Z->ee ratio" << endl;
	  h1->Scale( BoostExtrapolationFactorZ );
	}
      }

      // if DY->ee, mu->tau, scale by fake-rate
      if((it->first).find("DYEtoTau")!=string::npos){
	float sF = EtoTauCorrectionFactor;

	if(selection_.find("novbf")!=string::npos){
	  sF *= NoVbfExtrapolationFactorZ;
	  cout << "DY->tautau, e->tau will be rescaled by " << NoVbfExtrapolationFactorZ << " according to the Z->ee+vbf/Z->ee ratio" << endl;
	}
	else if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos){
	  sF *= VbfExtrapolationFactorZ;
	  cout << "DY->tautau, e->tau will be rescaled by " << VbfExtrapolationFactorZ << " according to the Z->ee+vbf/Z->ee ratio" << endl;
	}
	else if(selection_.find("boost")!=string::npos){
	  cout << "DY->tautau, e->tau will be rescaled by " << BoostExtrapolationFactorZ << " according to the Z->ee+vbf/Z->ee ratio" << endl;
	  sF *= BoostExtrapolationFactorZ;
	}
	h1->Scale(sF);
	hZmm->Add(h1,1.0);
	hZfakes->Add(h1,1.0);
      }

      // if DY->ee, jet->tau, scale by fake-rate
      if((it->first).find("DYJtoTau")!=string::npos){
	float sF = JtoTauCorrectionFactor;

	if(selection_.find("novbf")!=string::npos){
	  sF *= NoVbfExtrapolationFactorZ;
	  cout << "DY->tautau, jet->tau will be rescaled by " << NoVbfExtrapolationFactorZ << " according to the Z->ee+vbf/Z->ee ratio" << endl;
	}
	else if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos){
	  sF *= VbfExtrapolationFactorZ;
	  cout << "DY->tautau, jet->tau will be rescaled by " << VbfExtrapolationFactorZ << " according to the Z->ee+vbf/Z->ee ratio" << endl;
	}
	else if(selection_.find("boost")!=string::npos){
	  cout << "DY->tautau, jet->tau will be rescaled by " << BoostExtrapolationFactorZ << " according to the Z->ee+vbf/Z->ee ratio" << endl;
	  sF *=BoostExtrapolationFactorZ;
	}
	h1->Scale(sF);
	hZmj->Add(h1,1.0);
	hZfakes->Add(h1,1.0);
      }

      if((it->first).find("ggH115")!=string::npos && selection_.find("boost")!=string::npos){
	h1->Scale(kFactorSM);
      }


    }
  
    /////////////////////////////////////////////////////////////////////////////////////

    if( (it->first).find("DYEtoTau")!=string::npos ||  
	(it->first).find("DYJtoTau")!=string::npos || 
	(it->first).find("WJets")!=string::npos || 
	(it->first).find("Others")!=string::npos ) {
      hEWK->SetFillColor(kRed-3);
      hEWK->Add(h1,1.0);

      if( (it->first).find("Others")!=string::npos ){
	hVV->Add(h1,1.0);
	if(hVVKeys->Integral()>0) hVVKeys->Scale(hVV->Integral()/hVVKeys->Integral());
	hVVKeys->SetFillColor(kRed-3);
      }
    }

    if( (it->first).find("DYToTauTau")!=string::npos ) {
      hZtt->Add(h1,1.0);
      hZtt->SetFillColor(kYellow-9);
    }
    if( (it->first).find("Embedded")!=string::npos ) {
      if(hZtt->Integral()>0) h1->Scale(hZtt->Integral()/h1->Integral());
      h1->Scale(embeddedMEtCutEff/madgraphMEtCutEff);
      hDataEmb->Add(h1,1.0);
      hDataEmb->SetFillColor(kYellow-9);
    }
    if( (it->first).find("DYEtoTau")!=string::npos ) {
      if(hZmmKeys->Integral()>0) hZmmKeys->Scale(hZmm->Integral()/hZmmKeys->Integral());
      hZmmKeys->SetFillColor(kRed-3);
    }
    if( (it->first).find("DYJtoTau")!=string::npos ) {
      if(hZmjKeys->Integral()>0) hZmjKeys->Scale(hZmj->Integral()/hZmjKeys->Integral());
      hZmjKeys->SetFillColor(kRed-3);
    }
    if( (it->first).find("WJets")!=string::npos ) {
      if(hWKeys->Integral()>0) hWKeys->Scale(hW->Integral()/hWKeys->Integral());
      hWKeys->SetFillColor(kRed-3);
    }
    if( (it->first).find("W3Jets")!=string::npos ) {
      hW3Jets->Add(h1,1.0);
      if(hW3Jets->Integral()>0) hW3Jets->Scale(hW->Integral()/hW3Jets->Integral());
      hW3Jets->SetFillColor(kRed-3);
      if(hW3JetsKeys->Integral()>0) hW3JetsKeys->Scale(hW->Integral()/hW3JetsKeys->Integral());
      hW3JetsKeys->SetFillColor(kRed-3);
    }
    if( (it->first).find("LooseIso")!=string::npos ) {
      hLooseIso->Add(h1,1.0);
      if(hLooseIsoKeys->Integral()>0) hLooseIsoKeys->Scale(hLooseIso->Integral()/hLooseIsoKeys->Integral());
      hLooseIsoKeys->SetFillColor(0);
    }
    if( (it->first).find("AntiIso")!=string::npos ) {
      float sF = SSIsoToSSAIsoRatioQCD*OStoSSRatioQCD;
      cout << "Anti-isolated SS data scaled by " <<  SSIsoToSSAIsoRatioQCD << " ratio measured in inclusive sample ==>" << endl;
      cout << "   SS annti-isolated events = " << h1->GetEntries() << " ==> " << h1->Integral()*SSIsoToSSAIsoRatioQCD << " predicted in signal region" << endl;
      h1->Scale(sF);
      hAntiIso->Add(h1,1.0);
      if(hAntiIsoKeys->Integral()>0) hAntiIsoKeys->Scale(hAntiIso->Integral()/hAntiIsoKeys->Integral());
      hAntiIsoKeys->SetFillColor(0);
    }
    if( (it->first).find("TTbar")!=string::npos ) {
      hTTb->Add(h1,1.0);
      hTTb->SetFillColor(kBlue-2);     
    }
    if( (it->first).find("SS")!=string::npos ) {
      hQCD->Add(h1,1.0);
      hQCD->SetFillColor(kMagenta-9);
    }
    if((it->first).find("qqH115")!=string::npos){
      hSgn1->Add(h1,1.0);
      hSgn1->Scale(magnifySgn_);
      h1->Scale(magnifySgn_);
      hSgn1->SetLineWidth(2);
      h1->SetFillColor(kBlack);
      h1->SetFillStyle(3004);
      h1->SetLineColor(kBlack);
    }
    if((it->first).find("ggH115")!=string::npos){
      hSgn2->Add(h1,1.0);
      hSgn2->Scale(magnifySgn_);
      h1->Scale(magnifySgn_);
      hSgn2->SetLineWidth(2);
      h1->SetFillColor(kBlack);
      h1->SetFillStyle(3005);
      h1->SetLineColor(kBlack);
    }
    if((it->first).find("VH115")!=string::npos){
      hSgn3->Add(h1,1.0);
      hSgn3->Scale(magnifySgn_);
      h1->Scale(magnifySgn_);
      hSgn3->SetLineWidth(2);
      h1->SetFillColor(kBlack);
      h1->SetFillStyle(3005);
      h1->SetLineColor(kBlack);
    }
    if((it->first).find("Data")!=string::npos){
      hData->Add(h1,1.0);
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

    if(VERBOSE) cout<<(it->first) << " ==> " 
		    << h1->Integral() << " +/- " 
		    << TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries())
		    << endl;

    //out << (it->first) << "  " << h1->Integral() << " $\\pm$ " <<  TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries()) << endl;
   
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
  aStack->Add(hQCD);
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
  leg->AddEntry(hSgn,Form("(%.0fx) H#rightarrow#tau#tau m_{H}=%d",magnifySgn_,mH_),"F");
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
  if(logy_)
    hSgn->Draw("HISTSAME");
  
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
  if(selection_.find("novbf")!=string::npos || selection_.find("inclusive")!=string::npos){
    uncertZtt += (0.06  * 0.06) ; // Tau-ID 
    uncertZtt += (0.035 * 0.035); // Lumi 
  }
  else if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos){
    uncertZtt += (0.06  * 0.06) ; // Tau-ID 
    uncertZtt += (0.035 * 0.035); // Lumi 
    uncertZtt += (0.10  * 0.10);  // Extrap. factor 
  }
  else if(selection_.find("boost")!=string::npos){
    uncertZtt += (0.06  * 0.06) ; // Tau-ID 
    uncertZtt += (0.035 * 0.035); // Lumi 
    uncertZtt += (0.10  * 0.10);  // Extrap. factor 
  }
  uncertZtt = TMath::Sqrt(uncertZtt);
  float uncertTTb = 0;
  if(selection_.find("novbf")!=string::npos || selection_.find("inclusive")!=string::npos){
    uncertTTb += (0.075 * 0.075) ; // xsection 
  }
  else if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos){
    uncertTTb += (0.075 * 0.075) ; // xsection 
  }
  else if(selection_.find("boost")!=string::npos){
    uncertTTb += (0.075 * 0.075) ; // xsection 
  }
  uncertTTb = TMath::Sqrt(uncertTTb);
  float uncertEWK = 0;
  if(selection_.find("novbf")!=string::npos || selection_.find("inclusive")!=string::npos){
    uncertEWK += (0.07 * 0.07) ; // extrapolation 
  }
  else if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos){
    uncertEWK += (0.07 * 0.07) ; // extrapolation 
    uncertEWK += (0.10 * 0.10) ; // extrapolation 
  }
  else if(selection_.find("boost")!=string::npos){
    uncertEWK += (0.07 * 0.07) ; // extrapolation 
    uncertEWK += (0.10 * 0.10) ; // extrapolation  
  }
  uncertEWK = TMath::Sqrt(uncertEWK);
  float uncertQCD = 0;
  if(selection_.find("novbf")!=string::npos || selection_.find("inclusive")!=string::npos){
    uncertQCD += (0.02 * 0.02) ; // extrapolation 
  }
  else if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos){
    uncertQCD += (0.02 * 0.02) ; // extrapolation 
    uncertQCD += (0.10 * 0.10) ; // extrapolation 
  }
  else if(selection_.find("boost")!=string::npos){
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
    if(fabs(pull) > maxPull)
      maxPull = fabs(pull);
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

  c1->SaveAs(Form("/home/llr/cms/ivo/CMSSW_4_2_8_patch7/src/Bianchi/Limits/htautau/plots/%s/plot_eTau_mH%d_%s_%s_%s.png",outputDir.Data(), mH_,selection_.c_str(),analysis_.c_str(),variable_.Data()));
  c1->SaveAs(Form("/home/llr/cms/ivo/CMSSW_4_2_8_patch7/src/Bianchi/Limits/htautau/plots/%s/plot_eTau_mH%d_%s_%s_%s.pdf",outputDir.Data(), mH_,selection_.c_str(),analysis_.c_str(),variable_.Data()));


  // templates for fitting
  TFile* fout = new TFile(Form("/home/llr/cms/ivo/CMSSW_4_2_8_patch7/src/Bianchi/Limits/htautau/histograms/%s/eTau_mH%d_%s_%s_%s.root",outputDir.Data(), mH_,selection_.c_str(),analysis_.c_str(),variable_.Data()),"RECREATE");
  fout->cd();
  hSiml->Write();
  hQCD->Write();
  hZmm->Write();
  hZmj->Write();
  hZfakes->Write();
  hTTb->Write();
  hZtt->Write();
  hDataEmb->Write();
  hLooseIso->Write();
  hAntiIso->Write();
  hW->Write();
  hW3Jets->Write();
 
  TH1* hWKeysHisto = ((TH1Keys*)hWKeys)->GetHisto();
  if(hWKeysHisto){
    makeHistoFromDensity(hWKeysHisto,hW);
    hWKeysHisto->SetName("hWKeys");
    hWKeysHisto->Write();
  }
  TH1* hW3JetsKeysHisto = ((TH1Keys*)hW3JetsKeys)->GetHisto();
  if(hW3JetsKeysHisto){
    makeHistoFromDensity(hW3JetsKeysHisto,hW3Jets);
    hW3JetsKeysHisto->SetName("hW3JetsKeys");
    hW3JetsKeysHisto->Write();
  }
  TH1* hZmmKeysHisto = ((TH1Keys*)hZmmKeys)->GetHisto();
  if(hZmmKeysHisto){
    makeHistoFromDensity(hZmmKeysHisto,hZmm);
    hZmmKeysHisto->SetName("hZmmKeys");
    hZmmKeysHisto->Write();
  }
  TH1* hZmjKeysHisto = ((TH1Keys*)hZmjKeys)->GetHisto();
  if(hZmjKeysHisto){
    makeHistoFromDensity(hZmjKeysHisto,hZmj);
    hZmjKeysHisto->SetName("hZmjKeys");
    hZmjKeysHisto->Write();
  }
  TH1* hLooseIsoKeysHisto = ((TH1Keys*)hLooseIsoKeys)->GetHisto();
  if(hLooseIsoKeysHisto){
    makeHistoFromDensity(hLooseIsoKeysHisto,hLooseIso);
    hLooseIsoKeysHisto->SetName("hLooseIsoKeys");
    hLooseIsoKeysHisto->Write();
  }
  TH1* hAntiIsoKeysHisto = ((TH1Keys*)hAntiIsoKeys)->GetHisto();
  if(hAntiIsoKeysHisto){
    makeHistoFromDensity(hAntiIsoKeysHisto,hAntiIso);
    hAntiIsoKeysHisto->SetName("hAntiIsoKeys");
    hAntiIsoKeysHisto->Write();
  }
  TH1* hVVKeysHisto = ((TH1Keys*)hVVKeys)->GetHisto();
  if(hVVKeysHisto){
    makeHistoFromDensity(hVVKeysHisto, hVV);
    hVVKeysHisto->SetName("hVVKeys");
    hVVKeysHisto->Write();
  }
  
  hVV->Write();
  hEWK->Write();
  hSgn1->Write();
  hSgn2->Write();
  hSgn3->Write();
  hData->Write();
  hParameters->Write();
  fout->Write();
  

  fout->Close();
  


  delete hQCD; delete hZmm; delete hZmj; delete hZfakes; delete hTTb; delete hZtt; delete hW; delete hW3Jets; delete hLooseIso; delete hAntiIso;
  if(hW3JetsKeys)   delete hW3JetsKeys;
  if(hWKeys)        delete hWKeys;
  if(hZmmKeys)      delete hZmmKeys;
  if(hZmjKeys)      delete hZmjKeys;
  if(hVVKeys)       delete hVVKeys;
  if(hLooseIsoKeys) delete hLooseIsoKeys;
  if(hAntiIsoKeys)  delete hAntiIsoKeys;
  delete hVV; delete hSgn1; delete hSgn2; delete hSgn3; delete hData; delete hParameters;
  delete hWMt; delete aStack;  delete hEWK; delete hSiml; delete hDataEmb; delete hSgn; delete hRatio; delete line;
  delete fout;


  fSignalGGH->Close();       delete fSignalGGH; 
  fSignalVBF->Close();       delete fSignalVBF;
  fSignalVH->Close();        delete fSignalVH; 
  fBackgroundOthers->Close();delete fBackgroundOthers;
  fBackgroundTTbar->Close(); delete fBackgroundTTbar;
  fBackgroundWJets->Close(); delete fBackgroundWJets;
  fData->Close();            delete fData; 
  dummy1->Close();           delete dummy1;
  fBackgroundDY->Close();    delete fBackgroundDY;

}


///////////////////////////////////////////////////////////////////////////////////////////////



void plotElecTauAll( Int_t useEmbedded = 1, TString outputDir = "ElecTauStream_16Feb2012"){

  vector<string> variables;
  vector<int> mH;

  //variables.push_back("diTauVisMass");
  //variables.push_back("ptL1");
  //variables.push_back("etaL1");
  variables.push_back("tightestDanieleMVAWP");

  //mH.push_back(105);
  //mH.push_back(110);
  //mH.push_back(115);
  mH.push_back(120);
  //mH.push_back(125);
  //mH.push_back(130);
  //mH.push_back(135);
  //mH.push_back(140);
  //mH.push_back(145);
  //mH.push_back(160);

  //plotElecTau(120,0,"inclusive",""   ,"MtLeg1Corr","M_{T}","GeV" ,outputDir,40,0,160,5.0,1.0,0,1.2);

  //plotElecTau(120,0,"inclusive",""   ,"numPV","reconstructed vertexes","units" ,outputDir,30,0,30,5.0,1.0,0,1.5);

  //plotElecTau(120,1,"inclusive",""   ,"diTauVisMass","visible mass","GeV" ,outputDir,50,0,200,5.0,1.0,0,1.2);
  //plotElecTau(120,1,"inclusive",""   ,"ptL2","#tau p_{T}","GeV"           ,outputDir,30,0, 90,5.0,1.0,0,1.2);
  //plotElecTau(120,1,"inclusive",""   ,"ptL1","e p_{T}", "GeV"             ,outputDir,30,0, 90,5.0,1.0,0,1.2);

  //plotElecTau(120,1,"inclusive",""   ,"diTauSVFitMass","mass","GeV",outputDir,-1,0,100,5.0,1.0,0,1.2);
 
  //plotElecTau(120,0,"inclusive",""   ,"etaL1","e #eta", "units"           ,outputDir,25,-2.5, 2.5,5.0,1.0,0,2.);
  //plotElecTau(120,0,"inclusive",""   ,"etaL2","#tau #eta","units"         ,outputDir,25,-2.5, 2.5,5.0,1.0,0,2.);

  //plotElecTau(120,0,"inclusive",""   ,"nJets30","jet multiplicity","units",outputDir,10,0, 10,5.0,1.0,1,10);
  //plotElecTau(120,0,"inclusive",""   ,"nJets20BTagged","b-tagged jet multiplicity","units",outputDir,5,0, 5,5.0,1.0,1,10);

  //plotElecTau(120,1,"oneJet",""   ,"pt1","leading jet p_{T}","GeV"        ,outputDir,-1,0, 100,5.0,1.0,1,100);
  //plotElecTau(120,1,"oneJet",""   ,"eta1","leading jet #eta","units"      ,outputDir,21,-5, 5,5.0,1.0,0,2.);

  //plotElecTau(120,1,"twoJets",""   ,"pt1","leading jet p_{T}","GeV"       ,outputDir,-1,0, 100,5.0,1.0,1,100);
  //plotElecTau(120,1,"twoJets",""   ,"pt2","trailing jet p_{T}","GeV"      ,outputDir,-1,0, 100,5.0,1.0,1,100);
  //plotElecTau(120,1,"twoJets",""   ,"eta1","leading jet #eta","units"     ,outputDir,21,-5, 5,5.0,1.0,0,2.);
  //plotElecTau(120,1,"twoJets",""   ,"eta2","trailing jet #eta","units"    ,outputDir,21,-5, 5,5.0,1.0,0,2.);
  
  //plotElecTau(120,1,"twoJets",""   ,"Deta","|#Delta#eta|_{jj}","units"    ,outputDir,20,0, 8,   5.0,1.0,0,1.5);
  //plotElecTau(120,1,"twoJets",""   ,"Mjj","M_{jj}","GeV"                  ,outputDir,-1,0, 100,5.0,1.0,1,100);
  

  for(unsigned int i = 0 ; i < variables.size(); i++){
    for(unsigned j = 0; j < mH.size(); j++){

      //    plotElecTau(mH[j],useEmbedded,"inclusive",""   ,variables[i],"mass","GeV",outputDir,80,20,340,1.0,1.0,0,1.2);//for visMass
      //      plotElecTau(mH[j],useEmbedded,"inclusive",""   ,variables[i],"mass","GeV",outputDir,20,20,100,1.0,1.0,0,1.2);//for ptL1
      //     plotElecTau(mH[j],useEmbedded,"inclusive",""   ,variables[i],"e #eta","units",outputDir,5,-2.1,2.1,1.0,1.0,0,1.2);//for etaL1
      plotElecTau(mH[j],useEmbedded,"inclusive",""   ,variables[i],"mass","GeV",outputDir,150,-0.3,0.3,1.0,1.0,0,1.2);//for DanieleMVAdiscri

      //plotElecTau(mH[j],useEmbedded,"novbf",""       ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"novbf","TauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"novbf","TauDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"novbf","JetUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"novbf","JetDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
 
      //plotElecTau(mH[j],useEmbedded,"vbf",""         ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"vbf","TauUp"    ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"vbf","TauDown"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"vbf","JetUp"    ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"vbf","JetDown"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);

      //plotElecTau(mH[j],useEmbedded,"boost",""       ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"boost","TauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"boost","TauDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"boost","JetUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"boost","JetDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);

      //plotElecTau(mH[j],useEmbedded,"twoJets",""       ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"twoJets","TauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"twoJets","TauDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"twoJets","JetUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"twoJets","JetDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      
      //plotElecTau(mH[j],useEmbedded,"oneJet",""       ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"oneJet","TauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"oneJet","TauDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"oneJet","JetUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotElecTau(mH[j],useEmbedded,"oneJet","JetDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);      
    }
  }
  


}




int main(int argc, const char* argv[])
{

  std::cout << "plotElecTau()" << std::endl;
  gROOT->SetBatch(true);
 

  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  //plotElecTau();
  plotElecTauAll();
  //plotElecTau(120,1,"inclusive","","diTauVisMass","","","Dec2011/iter3",10,0,7000);
  //plotElecTau(120,1,"novbf","",    "diTauVisMass","","","Dec2011/iter3",10,0,7000);
  //plotElecTau(120,1,"boost","",    "diTauVisMass","","","Dec2011/iter3",10,0,7000);
  //plotElecTau(120,1,"vbf","",      "diTauVisMass","","","Dec2011/iter3",10,0,7000);

}

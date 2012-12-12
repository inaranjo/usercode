#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "./RooCMSShape.h"
#include "RooFormulaVar.h"
#include "RooConstVar.h"
#include "RooLandau.h"
#include "RooUniform.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooKeysPdf.h"
#include "RooBifurGauss.h"
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooAbsCategory.h"
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"
#include "RooTruthModel.h"
#include "RooDataHist.h"
#include "RooGenericPdf.h"
#include "RooConstVar.h"
#include "RooLognormal.h"
#include "RooProdPdf.h"

#include <vector>
#include <fstream>

#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"

#define DEBUG true

using namespace std;
using namespace RooFit;

double square(double x)
{
  return x*x;
}

              /*********/
              /* FIT   */
              /*********/

void simFit(Float_t out1[2][6],
	    const string tnp_          = "electauTnP",
	    const string categoryMC_   = "passingIsoLooseEleVetoMVA2Tight",
	    const string categoryDATA_ = "passingIsoLooseEleVetoMVA2Tight",
	    const string binLabel_     = "15to20",
	    double cutValue_           = 0.5,
	    const string bin_          = "Eta<1.5 && Eta>-1.5",
	    const float binCenter_     = 0.75,
	    const float binWidth_      = 0.75,
	    const float xLow_          = 70,
	    const float xHigh_         = 110,
	    const float nBins_         = 24,
	    bool doBinned_             = true,
	    float deltaAlpha_          = 0.0,
	    float deltaN_              = 0.0,
	    float mcLumi_              = 620,
	    bool doMC_                 = false,
	    bool PUWeight              = true,
	    bool SecondEleVeto         = true,
	    bool MVAIso                = false
	    )
{
#ifdef __CINT__
  gROOT->ProcessLineSync(".x RooCMSShape.cxx+") ;
#endif

  //Signal
  TFile fSgn("/data_CMS/cms/ivo/eTotauFakeRate/Trees/treeElecTauTnP_DYJetsToLL-8TeV-Sumer12-TnP-AntiEMVAv6.root");
  TTree *fullTreeSgn  = (TTree*)fSgn.Get((tnp_+"/fitter_tree").c_str());
  fSgn.cd("allEventsFilter");
  TH1F* totalEventsSgn = (TH1F*)gDirectory->Get("totalEvents");
//   float readEventsSgn = totalEventsSgn->GetBinContent(1);

  //data
//   TFile fdat("/data_CMS/cms/ivo/eTotauFakeRate/Trees/treeElecTauTnP_Run2012A-ElecTau-TnP-MVAIso.root");
  TFile fdat("/data_CMS/cms/ivo/eTotauFakeRate/Trees/treeElecTauTnP_Run2012B-ElecTau-TnP.root");
  TTree *fullTreeData = (TTree*)fdat.Get((tnp_+"/fitter_tree").c_str());

  //Compute MC truth efficiency
  TH1F* hS           = new TH1F("hS","",1,0,150);
  TH1F* hSP          = new TH1F("hSP","",1,0,150);
  
  if(PUWeight && !SecondEleVeto){
    fullTreeSgn->Draw("mass>>hS",Form("tag_PUMCWeight2012A*(%s && mass>%f && mass<%f && mcTrue && tag_GenDecay==23*11 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),xLow_,xHigh_));
    fullTreeSgn->Draw("mass>>hSP",Form("tag_PUMCWeight2012A*(%s && %s>=%f && mass>%f && mass<%f && mcTrue && tag_GenDecay==23*11 && DecayMode>0.5 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),categoryMC_.c_str(),cutValue_,xLow_,xHigh_));

//     if(MVAIso){
//       fullTreeSgn->Draw("mass>>hSP",Form("tag_PUMCWeight2012A*(%s && %s>=%f && mass>%f && mass<%f && mcTrue && tag_GenDecay==23*11 && HpsIsoMVALoose>0.5 && DecayMode>0.5 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),categoryMC_.c_str(),cutValue_,xLow_,xHigh_));
//     }
//     else {
//       fullTreeSgn->Draw("mass>>hSP",Form("tag_PUMCWeight2012A*(%s && %s>=%f && mass>%f && mass<%f && mcTrue && tag_GenDecay==23*11 && HpsLooseCombIsoDBCorr>0.5 && DecayMode>0.5 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),categoryMC_.c_str(),cutValue_,xLow_,xHigh_));
//     }
  }

  if(PUWeight && SecondEleVeto){
    
    fullTreeSgn->Draw("mass>>hS",Form("tag_PUMCWeight2012A*(%s && mass>%f && mass<%f && mcTrue && MatchElePassVeto<0.5 && tag_GenDecay==23*11 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),xLow_,xHigh_));
    fullTreeSgn->Draw("mass>>hSP",Form("tag_PUMCWeight2012A*(%s && %s>=%f && mass>%f && mass<%f && mcTrue && MatchElePassVeto<0.5 && tag_GenDecay==23*11 && DecayMode>0.5 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),categoryMC_.c_str(),cutValue_,xLow_,xHigh_));

//     if(MVAIso){
//       fullTreeSgn->Draw("mass>>hSP",Form("tag_PUMCWeight2012A*(%s && %s>=%f && mass>%f && mass<%f && mcTrue && MatchElePassVeto<0.5 && tag_GenDecay==23*11 && HpsIsoMVALoose>0.5 && DecayMode>0.5 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),categoryMC_.c_str(),cutValue_,xLow_,xHigh_));
//     }
//     else {
//       fullTreeSgn->Draw("mass>>hSP",Form("tag_PUMCWeight2012A*(%s && %s>=%f && mass>%f && mass<%f && mcTrue && MatchElePassVeto<0.5 && tag_GenDecay==23*11 && HpsLooseCombIsoDBCorr>0.5 && DecayMode>0.5 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),categoryMC_.c_str(),cutValue_,xLow_,xHigh_));
//     }
//     fullTreeSgn->Draw("mass>>hS",Form("tag_PUMCWeight2012A*(%s && mass>%f && mass<%f && mcTrue && tag_NumElePassVeto<2 && tag_GenDecay==23*11 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),xLow_,xHigh_));
//     fullTreeSgn->Draw("mass>>hSP",Form("tag_PUMCWeight2012A*(%s && %s>=%f && mass>%f && mass<%f && mcTrue && tag_NumElePassVeto<2 && tag_GenDecay==23*11 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),categoryMC_.c_str(),cutValue_,xLow_,xHigh_));

//     fullTreeSgn->Draw("mass>>hS",Form("tag_PUMCWeight2012A*(%s && mass>%f && mass<%f && mcTrue && passingIsoLooseSecondEleVeto<0.5 && tag_GenDecay==23*11 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),xLow_,xHigh_));
//     fullTreeSgn->Draw("mass>>hSP",Form("tag_PUMCWeight2012A*(%s && %s>=%f && mass>%f && mass<%f && mcTrue && passingIsoLooseSecondEleVeto<0.5 && tag_GenDecay==23*11 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),categoryMC_.c_str(),cutValue_,xLow_,xHigh_));

//     fullTreeSgn->Draw("mass>>hS",Form("tag_PUMCWeight2012A*(%s && mass>%f && mass<%f && mcTrue && MatchElePassVeto<0.5 && tag_GenDecay==23*11 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),xLow_,xHigh_));
//     fullTreeSgn->Draw("mass>>hSP",Form("tag_PUMCWeight2012A*(%s && %s>=%f && mass>%f && mass<%f && mcTrue && MatchElePassVeto<0.5 && HpsLooseCombIsoDBCorr>0.5 && DecayMode>0.5 && tag_GenDecay==23*11 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),categoryMC_.c_str(),cutValue_,xLow_,xHigh_));
  }

//   if(!PUWeight && !SecondEleVeto) {//No PUReweighting
//     fullTreeSgn->Draw("mass>>hS",Form("%s && mass>%f && mass<%f && mcTrue && tag_GenDecay==23*11 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match",bin_.c_str(),xLow_,xHigh_));
//     fullTreeSgn->Draw("mass>>hSP",Form("%s && %s>=%f && mass>%f && mass<%f && mcTrue && tag_GenDecay==23*11 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match",bin_.c_str(),categoryMC_.c_str(),cutValue_,xLow_,xHigh_));
//   }

//   if(!PUWeight && SecondEleVeto) {//No PUReweighting
//     fullTreeSgn->Draw("mass>>hS",Form("%s && mass>%f && mass<%f && mcTrue && tag_NumElePassVeto<2 && tag_GenDecay==23*11 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match",bin_.c_str(),xLow_,xHigh_));
//     fullTreeSgn->Draw("mass>>hSP",Form("%s && %s>=%f && mass>%f && mass<%f && mcTrue && tag_NumElePassVeto<2 && tag_GenDecay==23*11 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match",bin_.c_str(),categoryMC_.c_str(),cutValue_,xLow_,xHigh_));
//   }
  
  float SGNtrue = hS->Integral();
  float SGNtruePass = hSP->Integral();

  float McTruthEff    = SGNtruePass/SGNtrue;
  float BinomialError = TMath::Sqrt(SGNtruePass/SGNtrue*(1-SGNtruePass/SGNtrue)/SGNtrue);
  
  cout << bin_.c_str() << " ==> MCTRUTH: " << McTruthEff << " +/- " << BinomialError << endl;

  //   return;
  //   int pause ; 
  //    cin>>pause;

  delete hS; delete hSP;
 
  float reductionFactor = 1.0;

  //Trees of passing and failing for signal MC without pu reweight
  TTree* fullTreeSgnCutP = new TTree;
  TTree* fullTreeSgnCutF = new TTree;

  fullTreeSgnCutF = fullTreeSgn->CopyTree( Form("(%s< %f && %s && pair_tnpCharge==0 && event_met_pfmet<25 && mcTrue && tag_GenDecay==23*11 && tag_PfRelIso<0.1 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ,"",int(fullTreeSgn->GetEntries()*reductionFactor));
  fullTreeSgnCutP = fullTreeSgn->CopyTree( Form("(%s>=%f && %s && pair_tnpCharge==0 && event_met_pfmet<25 && mcTrue && tag_GenDecay==23*11 && DecayMode>0.5 &&tag_PfRelIso<0.1 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ,"",int(fullTreeSgn->GetEntries()*reductionFactor));

//   if(MVAIso){
//     fullTreeSgnCutP = fullTreeSgn->CopyTree( Form("(%s>=%f && %s && pair_tnpCharge==0 && event_met_pfmet<25 && mcTrue && tag_GenDecay==23*11 && HpsIsoMVALoose>0.5 && DecayMode>0.5 &&tag_PfRelIso<0.1 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ,"",int(fullTreeSgn->GetEntries()*reductionFactor));
//   }
//   else {
//     fullTreeSgnCutP = fullTreeSgn->CopyTree( Form("(%s< %f && %s && pair_tnpCharge==0 && event_met_pfmet<25 && mcTrue && tag_GenDecay==23*11 && HpsLooseCombIsoDBCorr>0.5 && DecayMode>0.5 && tag_PfRelIso<0.1 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ,"",int(fullTreeSgn->GetEntries()*reductionFactor));
//   }
  
  if(SecondEleVeto){
    fullTreeSgnCutF = fullTreeSgn->CopyTree( Form("(%s< %f && %s && pair_tnpCharge==0 && event_met_pfmet<25 && mcTrue && MatchElePassVeto<0.5 && tag_GenDecay==23*11 && tag_PfRelIso<0.1 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ,"",int(fullTreeSgn->GetEntries()*reductionFactor));
    fullTreeSgnCutP = fullTreeSgn->CopyTree( Form("(%s>=%f && %s && pair_tnpCharge==0 && event_met_pfmet<25 && mcTrue && MatchElePassVeto<0.5 && tag_GenDecay==23*11 && DecayMode>0.5 &&tag_PfRelIso<0.1 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ,"",int(fullTreeSgn->GetEntries()*reductionFactor));

//     if(MVAIso){
//       fullTreeSgnCutP = fullTreeSgn->CopyTree( Form("(%s>=%f && %s && pair_tnpCharge==0 && event_met_pfmet<25 && mcTrue && MatchElePassVeto<0.5 && tag_GenDecay==23*11 && HpsIsoMVALoose>0.5 && DecayMode>0.5 &&tag_PfRelIso<0.1 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ,"",int(fullTreeSgn->GetEntries()*reductionFactor));
//     }
//     else {
//       fullTreeSgnCutP = fullTreeSgn->CopyTree( Form("(%s< %f && %s && pair_tnpCharge==0 && event_met_pfmet<25 && mcTrue && MatchElePassVeto<0.5 && tag_GenDecay==23*11 && HpsLooseCombIsoDBCorr>0.5 && DecayMode>0.5 && tag_PfRelIso<0.1 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ,"",int(fullTreeSgn->GetEntries()*reductionFactor));
//     } 
    
  }
  /******************/
  /******************/
  /******************/
  /*ROOFIT VARIABLES*/
  /******************/
  /******************/
  /******************/

  RooRealVar mass("mass","m_{tp} (GeV/c^{2})",xLow_,xHigh_);
  mass.setBins( 10000, "fft" );
  mass.setBins( nBins_ );

              /****************/
              /* BACKGROUND   */
              /****************/

  RooRealVar meanBkgP("meanBkgP","",59,30,100);
  RooRealVar sigmaBkgP("sigmaBkgP","",11,0,50);
  //RooLandau bkgPdfP("bkgPdfP","",mass,meanBkgP,sigmaBkgP);

  RooRealVar DataCP("DataCP","",0,-10,10);//-0.5,-10,0);
  //RooExponential bkgPdfP("bkgPdfP","",mass,DataCP);

  RooRealVar alpha_CMSShapeP("alpha_CMSShapeP","alpha_CMSShapeP",65.,50.,78.);
  RooRealVar beta_CMSShapeP("beta_CMSShapeP","beta_CMSShapeP",0.01, 0.,0.1);
  RooRealVar gamma_CMSShapeP("gamma_CMSShapeP","gamma_CMSShapeP",0.06,0.,0.1);
  RooRealVar peak_CMSShapeP("peak_CMSShapeP","peak_CMSShapeP",91.2,87.,95.);
  if(categoryDATA_ == "passingIsoLooseEleVetoMedium" || categoryDATA_ == "AntiEleMedium"){
    alpha_CMSShapeP.setRange(50,100);
    alpha_CMSShapeP.setVal(75);
    beta_CMSShapeP.setRange(0.,0.1);
    beta_CMSShapeP.setVal(0.0005);
    gamma_CMSShapeP.setRange(0.,0.1);
    gamma_CMSShapeP.setVal(0.075);
    peak_CMSShapeP.setRange(87.,95);
    peak_CMSShapeP.setVal(91.2);
  }
  if(categoryDATA_ == "passingIsoLooseEleVetoMVA" || categoryDATA_ == "passingIsoMVALooseEleVetoMVA"){
    alpha_CMSShapeP.setRange(50,80);
    alpha_CMSShapeP.setVal(75);
    beta_CMSShapeP.setRange(0.,0.5);
    beta_CMSShapeP.setVal(0.07);
    gamma_CMSShapeP.setRange(0.,0.2);
    gamma_CMSShapeP.setVal(0.075);
    peak_CMSShapeP.setRange(87.,95);
    peak_CMSShapeP.setVal(91.2);
  }

  RooCMSShape bkgPdfP("bkgPdfP", "background p.d.f",mass,alpha_CMSShapeP,beta_CMSShapeP,gamma_CMSShapeP,peak_CMSShapeP);

  RooRealVar DataCF("DataCF","",0,-10,10);//-0.5,-10,0);
  //RooExponential bkgPdfF("bkgPdfF","",mass,DataCF);

  RooRealVar alpha_CMSShapeF("alpha_CMSShapeF","alpha_CMSShapeF",75.,50.,78.);
  RooRealVar beta_CMSShapeF("beta_CMSShapeF","beta_CMSShapeF",0.07,0,0.5);
  RooRealVar gamma_CMSShapeF("gamma_CMSShapeF","gamma_CMSShapeF",0.05,0,0.2);
  RooRealVar peak_CMSShapeF("peak_CMSShapeF","peak_CMSShapeF",91.2,87.,95.);
  if(categoryDATA_ == "passingIsoLooseEleVetoLoose" || categoryDATA_ == "AntiEleLoose"){
    alpha_CMSShapeF.setRange(50,80);
    alpha_CMSShapeF.setVal(75);
    beta_CMSShapeF.setRange(0.,0.1);
    beta_CMSShapeF.setVal(0.0005);
    gamma_CMSShapeF.setRange(0.,0.1);
    gamma_CMSShapeF.setVal(0.025);
    peak_CMSShapeF.setRange(87.,95.);
    peak_CMSShapeF.setVal(91.2);
  }  
  if(categoryDATA_ == "passingIsoLooseEleVetoTight" || categoryDATA_ == "AntiEleTight"){
    alpha_CMSShapeF.setRange(50,75);
    alpha_CMSShapeF.setVal(75);
    beta_CMSShapeF.setRange(0.,0.1);
    beta_CMSShapeF.setVal(0.0005);
    gamma_CMSShapeF.setRange(0.,0.1);
    gamma_CMSShapeF.setVal(0.025);
    peak_CMSShapeF.setRange(87.,95.);
    peak_CMSShapeF.setVal(91.2);
  }  

  RooCMSShape bkgPdfF("bkgPdfF", "background p.d.f",mass,alpha_CMSShapeF,beta_CMSShapeF,gamma_CMSShapeF,peak_CMSShapeF);

  mass.setBins( 50 );

              /****************/
              /* SIGNAL       */
              /****************/

  ///////////////////////////////////////////////////////////// 
  // passing:
  /////////////////////////////////////////////////////////////

  RooDataSet sgnDataSetP("sgnDataSetP","dataset for signal", RooArgSet(mass), Import( *fullTreeSgnCutP ) );
  RooDataHist sgnDataHistP("sgnDataHistP","",RooArgSet(mass),sgnDataSetP, 1.0);
  //RooHistPdf  sgnTemplatePdfP("sgnTemplatePdfP","",RooArgSet(mass),sgnDataHistP);
  //RooKeysPdf sgnTemplatePdfP("sgnTemplatePdfP","",mass,sgnDataSetP); //Signal template from MC

  // Breit-Wigner constrained for Z 
  RooRealVar meanSgnP("meanSgnP","mean",91.19/*,85,100*/);
  RooRealVar widthSgnP("widthSgnP","width",2.49/*,0.,10*/);
  RooBreitWigner bwSgnP("bwSgnP","bw",mass,meanSgnP,widthSgnP);
  if(categoryDATA_ == "passingIsoMVALooseEleVetoLoose" || categoryDATA_ == "AntiEleLoose"){
    meanSgnP.setRange(89,92);
  }  

  // Crystall Ball parameters free
  RooRealVar m1SgnP("m1SgnP","m1",0,-20,20);
  RooRealVar sigmaSgnP("sigmaSgnP","sigma",0.5,0,15);
  RooRealVar alfaSgnP("alfaSgnP","alfa", 0.5,-10,20);
  RooRealVar nSgnP("nSgnP","n", 1,1e-06,70);
  RooCBShape cbSgnP("cbSgnP","",mass,m1SgnP,sigmaSgnP,alfaSgnP,nSgnP);

  mass.setBins( 1000 , "fft");
  // BW (X) CB 
  RooFFTConvPdf bvcbSgnP("bvcbSgnP","",mass,bwSgnP, cbSgnP);

  //fit it to the signal dataset and retrieve the parameters
  RooFitResult* ResSgnFitP = bvcbSgnP.fitTo(sgnDataSetP, Minos(1), Save(1), NumCPU(4) );
  RooArgSet FitParamSgnP(ResSgnFitP->floatParsFinal());
  RooRealVar* m1SgnFitP    = (RooRealVar*)(&FitParamSgnP["m1SgnP"]);
  RooRealVar* sigmaSgnFitP = (RooRealVar*)(&FitParamSgnP["sigmaSgnP"]);
  RooRealVar* alfaSgnFitP  = (RooRealVar*)(&FitParamSgnP["alfaSgnP"]);
  RooRealVar* nSgnFitP     = (RooRealVar*)(&FitParamSgnP["nSgnP"]);
  //RooRealVar* nMeanSgnFitP = (RooRealVar*)(&FitParamSgnP["meanSgnP"]);
  //RooRealVar* nWidthSgnFitP= (RooRealVar*)(&FitParamSgnP["widthSgnP"]);

  //constrained parameters for CB
  RooRealVar m1SgnP_C("m1SgnP_C","m1",m1SgnFitP->getVal(),-5,5);
  RooRealVar sigmaSgnP_C("sigmaSgnP_C","sigma",sigmaSgnFitP->getVal(),0,15);
  RooRealVar alfaSgnP_C("alfaSgnP_C","alfa",alfaSgnFitP->getVal()*(1+deltaAlpha_),0,20);
  RooRealVar nSgnP_C("nSgnP_C","n",nSgnFitP->getVal()*(1+deltaN_),0,50);

  // choose to let BW parameters float or not
  //RooRealVar meanSgnP_C("meanSgnP_C","mean",  nMeanSgnFitP->getVal() ,80,120);
  //RooRealVar widthSgnP_C("widthSgnP_C","width",nWidthSgnFitP->getVal() /*,0.,10*/);

  //constrained CB
  RooCBShape cbSgnP_C("cbSgnP_C","",mass,m1SgnP_C,sigmaSgnP_C,alfaSgnP_C,nSgnP_C);
  
  RooLognormal alfaSgnP_CPdf("alfaSgnP_CPdf","",alfaSgnP_C,RooConst(alfaSgnFitP->getVal()),RooConst(1.5));
  RooLognormal nSgnP_CPdf("nSgnP_CPdf","",nSgnP_C,RooConst(nSgnFitP->getVal()),            RooConst(1.5));
  //RooLognormal meanSgnP_CPdf("meanSgnP_CPdf","",meanSgnP_C,RooConst(nMeanSgnFitP->getVal()),RooConst(1.5));
  //RooLognormal widthSgnP_CPdf("widthSgnP_CPdf","",widthSgnP_C,RooConst(nWidthSgnFitP->getVal()),RooConst(1.5));

  ////////////////////////////////////////////////////////////////////////////
  /////////////////// fitted BW (X) CB
  RooBreitWigner bwSgnP_C("bwSgnP_C","bw",mass,meanSgnP/*meanSgnP_C*/,widthSgnP/*widthSgnP_C*/);
  RooFFTConvPdf sgnPdfP("sgnPdfP","",mass,bwSgnP_C, cbSgnP_C);
  ///////////////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////Gaussian(X)Sgn template from MC
  //Gaussian
  RooRealVar sgnMeanResP("sgnMeanResP","",0,-5,5);
//   RooRealVar sgnMeanResP("sgnMeanResP","",0,-10,10);
  RooRealVar sgnSigmaResP("sgnSigmaResP","",0.5,0,10);
  RooGaussian resolModP("sgnResolModP","",mass,sgnMeanResP,sgnSigmaResP);  
  mass.setBins( 10000, "fft" );
  //RooFFTConvPdf sgnPdfP("sgnPdfP","",mass,sgnTemplatePdfP, resolModP);
///////////////////////////////////////////////////////////////////////////////////////////////
  mass.setBins(nBins_);
  //return;

  /////////////////////////////////////////////////////////////
  // failing:
  ///////////////////////////////////////////////////////////// 

  RooDataSet sgnDataSetF("sgnDataSetF","dataset for signal", RooArgSet(mass), Import( *fullTreeSgnCutF ) );
  RooDataHist sgnDataHistF("sgnDataHistF","",RooArgSet(mass),sgnDataSetF, 1.0);
  //RooHistPdf  sgnTemplatePdfF("sgnTemplatePdfF","",RooArgSet(mass),sgnDataHistF);
  //RooKeysPdf sgnTemplatePdfF("sgnTemplatePdfF","",mass,sgnDataSetF); //Template from MC

  // Breit-Wigner
  RooRealVar meanSgnF("meanSgnF","mean",91.19/*,85,100*/);
  RooRealVar widthSgnF("widthSgnF","width",2.49/*,0.,10*/);
  RooBreitWigner bwSgnF("bwSgnF","bw",mass,meanSgnF,widthSgnF);

  // Crystall Ball
  RooRealVar m1SgnF("m1SgnF","m1",0,-20,20);
  RooRealVar sigmaSgnF("sigmaSgnF","sigma",0.5,0,20);
  RooRealVar alfaSgnF("alfaSgnF","alfa", 0.5,-10,20);
  RooRealVar nSgnF("nSgnF","n", 1,1e-06,70);
  RooCBShape cbSgnF("cbSgnF","",mass,m1SgnF,sigmaSgnF,alfaSgnF,nSgnF);

  // BW (X) CB
  RooFFTConvPdf bvcbSgnF("bvcbSgnF","",mass,bwSgnF, cbSgnF);

  //Fit it to failing dataset
  RooFitResult* ResSgnFitF = bvcbSgnF.fitTo(sgnDataSetF, Minos(1), Save(1), NumCPU(4) );
  RooArgSet FitParamSgnF(ResSgnFitF->floatParsFinal());
  RooRealVar* m1SgnFitF    = (RooRealVar*)(&FitParamSgnF["m1SgnF"]);
  RooRealVar* sigmaSgnFitF = (RooRealVar*)(&FitParamSgnF["sigmaSgnF"]);
  RooRealVar* alfaSgnFitF  = (RooRealVar*)(&FitParamSgnF["alfaSgnF"]);
  RooRealVar* nSgnFitF     = (RooRealVar*)(&FitParamSgnF["nSgnF"]);
  //RooRealVar* nMeanSgnFitF = (RooRealVar*)(&FitParamSgnF["meanSgnF"]);
  //RooRealVar* nWidthSgnFitF= (RooRealVar*)(&FitParamSgnF["widthSgnF"]);

  RooRealVar m1SgnF_C("m1SgnF_C","m1",/*-0.7*/m1SgnFitF->getVal(),-10,10);
  RooRealVar sigmaSgnF_C("sigmaSgnF_C","sigma",/*2.0*/sigmaSgnFitF->getVal(),0,20);
  RooRealVar alfaSgnF_C("alfaSgnF_C","alfa",/*17.9*/alfaSgnFitF->getVal()*(1+deltaAlpha_),0,20);
  RooRealVar nSgnF_C("nSgnF_C","n",/*1.95*/nSgnFitF->getVal()*(1+deltaN_),0,50);

  // choose to let it float or not
  // RooRealVar meanSgnF_C("meanSgnF_C","mean",  nMeanSgnFitF->getVal() ,80,120);
  //RooRealVar widthSgnF_C("widthSgnF_C","width",nWidthSgnFitF->getVal() /*,0.,10*/);

  //Constrained CB
  RooCBShape cbSgnF_C("cbSgnF_C","",mass,m1SgnF_C,sigmaSgnF_C,alfaSgnF_C,nSgnF_C);
  
  RooLognormal alfaSgnF_CPdf("alfaSgnF_CPdf","",alfaSgnF_C,RooConst(alfaSgnFitF->getVal()),RooConst(1.5));
  RooLognormal nSgnF_CPdf("nSgnF_CPdf","",nSgnF_C,RooConst(nSgnFitF->getVal()),            RooConst(1.5));

  //RooLognormal meanSgnF_CPdf("meanSgnF_CPdf","",meanSgnF_C,RooConst(nMeanSgnFitF->getVal()),RooConst(1.5));
  //RooLognormal widthSgnF_CPdf("widthSgnF_CPdf","",widthSgnF_C,RooConst(nWidthSgnFitF->getVal()),RooConst(1.5));

  ////////////////////////////////////////////////////////////////////////////
  /////////////////// fitted BW (X) CB
  RooBreitWigner bwSgnF_C("bwSgnF_C","bw",mass,meanSgnF/*meanSgnF_C*/,widthSgnF/*widthSgnF_C*/);
  RooFFTConvPdf sgnPdfF("sgnPdfF","",mass,bwSgnF_C, cbSgnF_C);
  ///////////////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////Gaussian(X)Sgn template from MC
  RooRealVar sgnMeanResF("sgnMeanResF","",0,-10,10);
  RooRealVar sgnSigmaResF("sgnSigmaResF","",0.5,0,10);
  RooGaussian resolModF("sgnResolModF","",mass,sgnMeanResF,sgnSigmaResF);
  mass.setBins( 10000, "fft" );
  //RooFFTConvPdf sgnPdfF("sgnPdfF","",mass,sgnTemplatePdfF, resolModF);
  ///////////////////////////////////////////////////////////////////////////////////////////////

  mass.setBins(nBins_);
  //return;

  /******************/
  /******************/
  /******************/
  /*FIT             */
  /******************/
  /******************/
  /******************/

  RooCategory category("category","category") ;
  category.defineType("pass") ;
  category.defineType("fail") ;
  //Variables 
  RooRealVar DataNumBkgF("DataNumBkgF","",0,10000000);
  RooRealVar DataNumBkgP("DataNumBkgP","",0,10000000);
  RooRealVar DataNumSgn("DataNumSgn","",  0,10000000);
  RooRealVar DataEfficiency("DataEfficiency","",0.5,0,1);

  RooFormulaVar DataNumSgnP("DataNumSgnP","DataEfficiency*DataNumSgn",    RooArgSet(DataEfficiency,DataNumSgn));
  RooFormulaVar DataNumSgnF("DataNumSgnF","(1-DataEfficiency)*DataNumSgn",RooArgSet(DataEfficiency,DataNumSgn));
  
  //Pdf models for passing and failing events extracted from before
  RooAddPdf DataModelP("DataModelP","",RooArgList(/*cbSgnP_C*/ sgnPdfP /*sgnTemplatePdfP*/,bkgPdfP),RooArgList(DataNumSgnP,DataNumBkgP));
  RooAddPdf DataModelF("DataModelF","",RooArgList(/*cbSgnF_C*/ sgnPdfF /*sgnTemplatePdfF*/,bkgPdfF),RooArgList(DataNumSgnF,DataNumBkgF));
  
  TTree* fullTreeDataCutF;
  fullTreeDataCutF = fullTreeData->CopyTree( Form("(%s <%f && %s   && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_HLTEle27Match && tag_Pt>25)",categoryDATA_.c_str(),cutValue_,bin_.c_str()) );
  TTree* fullTreeDataCutP;
  if(MVAIso){
    fullTreeDataCutP = fullTreeData->CopyTree( Form("(%s>=%f && %s  && HpsIsoMVALoose>0.5 && DecayMode>0.5 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_HLTEle27Match && tag_Pt>25)",categoryDATA_.c_str(),cutValue_,bin_.c_str()) ); 
  }
  else{
    fullTreeDataCutP = fullTreeData->CopyTree( Form("(%s>=%f && %s  && HpsLooseCombIsoDBCorr>0.5 && DecayMode>0.5 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_HLTEle27Match && tag_Pt>25)",categoryDATA_.c_str(),cutValue_,bin_.c_str()) ); 
  }

  if(SecondEleVeto){
    fullTreeDataCutF = fullTreeData->CopyTree( Form("(%s <%f && %s   && MatchElePassVeto<0.5 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_HLTEle27Match && tag_Pt>25)",categoryDATA_.c_str(),cutValue_,bin_.c_str()) );
    if(MVAIso){
      fullTreeDataCutP = fullTreeData->CopyTree( Form("(%s>=%f && %s  && MatchElePassVeto<0.5 && HpsIsoMVALoose>0.5 && DecayMode>0.5 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_HLTEle27Match && tag_Pt>25)",categoryDATA_.c_str(),cutValue_,bin_.c_str()) ); 
    }
    else{
      fullTreeDataCutP = fullTreeData->CopyTree( Form("(%s>=%f && %s  && MatchElePassVeto<0.5 && HpsLooseCombIsoDBCorr>0.5 && DecayMode>0.5 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_HLTEle27Match && tag_Pt>25)",categoryDATA_.c_str(),cutValue_,bin_.c_str()) ); 
    }
  }
  //Roofit datasets of data passing and failing cut
  mass.setBins(nBins_);
  RooDataSet DataDataSetP("DataDataSetP","dataset for Data pass", RooArgSet(mass), Import( *fullTreeDataCutP ) );
  std::cout << "data dataset Pass " << DataDataSetP.numEntries() << "  " << std::endl;
  RooDataHist DataDataHistP("DataDataHistP","",RooArgSet(mass),DataDataSetP, 1.0);
  RooDataSet DataDataSetF("DataDataSetF","dataset for Data fail", RooArgSet(mass), Import( *fullTreeDataCutF ) );
  std::cout << "data dataset Fail " << DataDataSetF.numEntries() << "  " << std::endl;
  RooDataHist DataDataHistF("DataDataHistF","",RooArgSet(mass),DataDataSetF, 1.0);

  RooRealVar DataNumSgnP_("DataNumSgnP_","",0,10000);
  RooAddPdf DataModelP_("DataModelP_","",RooArgList(sgnPdfP,bkgPdfP),RooArgList(DataNumSgnP_,DataNumBkgP));
  DataModelP_.fitTo(DataDataSetP, Extended(1), Minos(1), Save(1), NumCPU(4),SumW2Error(1) /*,ExternalConstraints( RooArgSet(meanSgn_CPdf,widthSgn_CPdf) )*/);
  
  RooPlot* frame2 = mass.frame(Title("template"));
  DataDataSetP.plotOn(frame2);
  DataModelP_.plotOn(frame2, LineColor(kBlue), LineStyle(kSolid));
  DataModelP_.plotOn(frame2, Components("sgnPdfP"), LineColor(kRed), LineStyle(kSolid));
  DataModelP_.plotOn(frame2, Components("bkgPdfP"), LineColor(kGreen), LineStyle(kSolid));
  frame2->Draw();
  
  //return;
  
  // binned combined dataset
  RooDataHist DataCombData("DataCombData","combined data",mass,Index(category),Import("pass", *(DataDataSetP.createHistogram("histoDataP",mass)) ) ,Import("fail", *(DataDataSetF.createHistogram("histoDataF",mass))), Weight(0.5) ) ;
  std::cout << "data dataHist Comb " << DataCombData.sumEntries() << "  " << std::endl;
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  // unbinned combined dataset
  RooDataSet DataCombDataUnBinned("DataCombDataUnBinned","combined data",mass,Index(category),Import("pass", DataDataSetP ) ,Import("fail",DataDataSetF), Weight(0.5) ) ;
  std::cout << "data dataset Comb " << DataCombDataUnBinned.numEntries() << "  " << std::endl;
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  //return out;
  
  RooSimultaneous DataSimPdf("DataSimPdf","simultaneous pdf",category) ;
  DataSimPdf.addPdf(DataModelP,"pass") ;
  DataSimPdf.addPdf(DataModelF,"fail") ;
  
  //mass.setBins( 10000, "fft" );
  mass.setBins( nBins_ );
  RooFitResult* ResDataCombinedFit =  0;
  if(doBinned_)  ResDataCombinedFit = DataSimPdf.fitTo(DataCombData , Extended(1), Minos(1), Save(1), NumCPU(4), /*ExternalConstraints( RooArgSet(alfaSgn_CPdf,nSgn_CPdf) )*/  SumW2Error(1));
  else ResDataCombinedFit = DataSimPdf.fitTo(DataCombDataUnBinned , Extended(1), Minos(1), Save(1), NumCPU(4),  /*ExternalConstraints( RooArgSet(alfaSgn_CPdf,nSgn_CPdf) )*/ SumW2Error(1));
  
  //Get the fitted efficiency and number of signal events
  RooArgSet DataFitParam(ResDataCombinedFit->floatParsFinal());
  RooRealVar* DataEffFit      = (RooRealVar*)(&DataFitParam["DataEfficiency"]);
  //RooRealVar* DataNumSigFit   = (RooRealVar*)(&DataFitParam["DataNumSgn"]);
  
  string theSample = doMC_ ? "Simulation" : "Data" ;
  float theLumi = doMC_ ? mcLumi_ : 2676;//620
  
  RooPlot* DataFrameP = mass.frame(Bins(40),Title(Form("CMS Preliminary 2012B  #sqrt{s}=8 TeV %s  L=%.0f pb^{-1}:  passing probe",theSample.c_str(),theLumi)));
  DataCombData.plotOn(DataFrameP,Cut("category==category::pass"),Name("dataP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), LineColor(kBlue),Name("modelP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("sgnPdfP"/*"cbSgnP_C"*/), LineColor(kRed), LineStyle(kSolid),Name("signal onlyP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("bkgPdfP"), LineColor(kMagenta), LineStyle(kDashed),Name("backgroundP"));
  DataFrameP->SetTitleOffset(1,"Y");
  DataFrameP->SetTitleSize(1,"Y");

  RooPlot* DataFrameF = mass.frame(Bins(40),Title(Form("CMS Preliminary 2012B  #sqrt{s}=8 TeV %s  L=%.0f pb^{-1}:  failing probe",theSample.c_str(),theLumi)));
  DataCombData.plotOn(DataFrameF,Cut("category==category::fail"),Name("dataF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), LineColor(kBlue),Name("modelF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("sgnPdfF"/*"cbSgnF_C"*/), LineColor(kRed), LineStyle(kSolid),Name("signal onlyF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("bkgPdfF"), LineColor(kMagenta), LineStyle(kDashed),Name("backgroundF"));
  DataFrameF->SetTitleOffset(1,"Y");
  DataFrameF->SetTitleSize(1,"Y");

  TCanvas *cPass = new TCanvas("fitCanvasP","canvas",10,30,650,600);
  cPass->SetGrid(0,0);
  cPass->SetFillStyle(4000);
  cPass->SetFillColor(10);
  cPass->SetTicky();
  cPass->SetObjectStat(0);
      
  cPass->cd();
  DataFrameP->Draw();
  TLegend *leg1 = new TLegend(0.65,0.73,0.86,0.87);
  leg1->SetFillColor(kWhite);
  leg1->SetLineColor(kWhite);
  leg1->AddEntry("dataP","Data", "P");
  leg1->AddEntry("modelP","Signal + background","L");
  leg1->AddEntry("backgroundP","Background only", "L");
  leg1->AddEntry("signal onlyP","Signal only", "L");
  leg1->Draw();
  //return;

  string fileNameP = "Pass_"+tnp_+"_"+categoryDATA_;
  if(doMC_) fileNameP = fileNameP+"_MC";
  if (!doMC_) fileNameP = fileNameP+"_DATA";
  if(MVAIso)cPass->SaveAs(Form("./ElecTauTnP/plots/computeElecTauTnPMacro_Fit_MVAIso_%s_%s.pdf",fileNameP.c_str(),binLabel_.c_str()));
  else cPass->SaveAs(Form("./ElecTauTnP/plots/computeElecTauTnPMacro_Fit_%s_%s.pdf",fileNameP.c_str(),binLabel_.c_str()));

  TCanvas *cFail = new TCanvas("fitCanvasF","canvas",10,30,650,600);
  cFail->SetGrid(0,0);
  cFail->SetFillStyle(4000);
  cFail->SetFillColor(10);
  cFail->SetTicky();
  cFail->SetObjectStat(0);

  cFail->cd();
  DataFrameF->Draw();
  TLegend *leg2 = new TLegend(0.65,0.73,0.86,0.87);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);
  leg2->AddEntry("dataF","Data", "P");
  leg2->AddEntry("modelF","Signal + background","L");
  leg2->AddEntry("backgroundF","Background only", "L");
  leg2->AddEntry("signal onlyF","Signal only", "L");
  leg2->Draw();

  string fileNameF = "Fail_"+tnp_+"_"+categoryDATA_;
  if(doMC_) fileNameF = fileNameF+"_MC";
  else fileNameF = fileNameF+"_DATA";
  if(MVAIso)cFail->SaveAs(Form("./ElecTauTnP/plots/computeElecTauTnPMacro_Fit_MVAIso_%s_%s.pdf",fileNameF.c_str(), binLabel_.c_str()));
  else cFail->SaveAs(Form("./ElecTauTnP/plots/computeElecTauTnPMacro_Fit_%s_%s.pdf",fileNameF.c_str(), binLabel_.c_str()));
  ResDataCombinedFit->printArgs(std::cout);
  cout << endl;
  ResDataCombinedFit->printValue(std::cout);
  cout << endl;

  float DataErrorLo = DataEffFit->getErrorLo()<0 ? DataEffFit->getErrorLo() : (-1)*DataEffFit->getError();
  float DataErrorHi = DataEffFit->getErrorHi()>0 ? DataEffFit->getErrorHi() : DataEffFit->getError();

  cout <<"Data efficiency fit value : "<< DataEffFit->getVal() << " +/- " << DataEffFit->getError() << "  ( " << DataErrorLo << ", " << DataErrorHi << ")" <<  endl;

  out1[0][0]=(binCenter_);
  out1[0][1]=(binWidth_);
  out1[0][2]=(binWidth_);
  out1[0][3]=(McTruthEff);
  out1[0][4]=(BinomialError);
  out1[0][5]=(BinomialError);

  out1[1][0]=(binCenter_);
  out1[1][1]=(binWidth_);
  out1[1][2]=(binWidth_);
  out1[1][3]=(DataEffFit->getVal());
  out1[1][4]=((-1)*DataErrorLo);
  out1[1][5]=(DataErrorHi);

  return;

}//end simFit












/*************************/
/* FIT  FOR CLOSURE TEST */
/************************/

void simFitClosure(Float_t out1[2][6],
		   const string tnp_          = "electauTnP",
		   const string categoryMC_   = "passingIsoLooseEleVetoMVA2Tight",
		   const string categoryDATA_ = "passingIsoLooseEleVetoMVA2Tight",
		   const string binLabel_     = "15to20",
		   double cutValue_           = 0.5,
		   const string bin_          = "Eta<1.5 && Eta>-1.5",
		   const float binCenter_     = 0.75,
		   const float binWidth_      = 0.75,
		   const float xLow_          = 70,
		   const float xHigh_         = 110,
		   const float nBins_         = 24,
		   bool doBinned_             = true,
		   float deltaAlpha_          = 0.0,
		   float deltaN_              = 0.0,
		   float mcLumi_              = 620,
		   bool doMC_                 = false,
		   bool PUWeight              = true,
		   bool SecondEleVeto         = true,
		   bool MVAIso                = false
		   )
{
#ifdef __CINT__
  gROOT->ProcessLineSync(".x RooCMSShape.cxx+") ;
#endif

  //Signal
  TFile fSgn("/data_CMS/cms/ivo/eTotauFakeRate/Trees/treeElecTauTnP_DYJetsToLL-8TeV-Sumer12-TnP-AntiEMVAv6.root");
  TTree *fullTreeSgn  = (TTree*)fSgn.Get((tnp_+"/fitter_tree").c_str());
  fSgn.cd("allEventsFilter");
//   TH1F* totalEventsSgn = (TH1F*)gDirectory->Get("totalEvents");
//   float readEventsSgn = totalEventsSgn->GetBinContent(1);

  //data
  //  TFile fdat("/data_CMS/cms/ivo/eTotauFakeRate/Trees/treeElecTauTnP_Run2012A-ElecTau-TnP-MVAIso.root");
  TFile fdat("/data_CMS/cms/ivo/eTotauFakeRate/Trees/treeElecTauTnP_Run2012B-ElecTau-TnP.root");
  TTree *fullTreeData = (TTree*)fdat.Get((tnp_+"/fitter_tree").c_str());

  //Mc Bkg
  TFile fZtau("/data_CMS/cms/ivo/eTotauFakeRate/Trees/treeElecTauTnP_DYJetsToLL-8TeV-Sumer12-TnP-AntiEMVAv6.root");
  TTree *fullTreeZtt  = (TTree*)fZtau.Get((tnp_+"/fitter_tree").c_str());
  fZtau.cd("allEventsFilter");
  TH1F* totalEventsZtt = (TH1F*)gDirectory->Get("totalEvents");
  float readEventsZtt = totalEventsZtt->GetBinContent(1);

  TFile fWen("/data_CMS/cms/ivo/eTotauFakeRate/Trees/treeElecTauTnP_WJetsToLNu-8TeV-Summer12-TnP-TnP-AntiEMVAv6.root");
  TTree *fullTreeWen = (TTree*)fWen.Get((tnp_+"/fitter_tree").c_str());
  fWen.cd("allEventsFilter");
  TH1F* totalEventsWen = (TH1F*)gDirectory->Get("totalEvents");
  float readEventsWen = totalEventsWen->GetBinContent(1);

  TFile fTTb("/data_CMS/cms/ivo/eTotauFakeRate/Trees/treeElecTauTnP_TTJets-ElecTau-8TeV-Summer12-TnP-TnP-AntiEMVAv6.root");
  TTree *fullTreeTTb  = (TTree*)fTTb.Get((tnp_+"/fitter_tree").c_str());
  fTTb.cd("allEventsFilter");
  TH1F* totalEventsTTb = (TH1F*)gDirectory->Get("totalEvents");
  float readEventsTTb = totalEventsTTb->GetBinContent(1);

  //Cross sections
  float SigmaZLL = 3503.71;
  float SigmaWen = 36257.2;
  float SigmaTTb = 225.197;
  //Processed lumi
  float LumiZtt = readEventsZtt/SigmaZLL;
  float LumiWen = readEventsWen/SigmaWen;
  float LumiTTb = readEventsTTb/SigmaTTb;
  if(DEBUG){
    cout<<"Luminosity processed for Ztt :"<<LumiZtt<<" pb-1"<<endl;
    cout<<"Luminosity processed for Wen :"<<LumiWen<<" pb-1"<<endl;
    cout<<"Luminosity processed for TTb :"<<LumiTTb<<" pb-1"<<endl;
    cout<<"Luminosity used for the fits :"<<mcLumi_<<" pb-1"<<endl;
  }
  if(doMC_ && (mcLumi_>LumiZtt || mcLumi_>LumiWen || mcLumi_>LumiTTb)){
    cout<<" Warning ! Luminosity exceeds MC processed luminosity !"<<endl;
  }

  //Compute MC truth efficiency
  TH1F* hS           = new TH1F("hS","",1,0,150);
  TH1F* hSP          = new TH1F("hSP","",1,0,150);
  
  if(PUWeight && !SecondEleVeto){
    fullTreeSgn->Draw("mass>>hS",Form("tag_PUMCWeight2012A*(%s && mass>%f && mass<%f && mcTrue && tag_GenDecay==23*11 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),xLow_,xHigh_));
    if(MVAIso){
      fullTreeSgn->Draw("mass>>hSP",Form("tag_PUMCWeight2012A*(%s && %s>=%f && mass>%f && mass<%f && mcTrue && tag_GenDecay==23*11 && HpsIsoMVALoose>0.5 && DecayMode>0.5 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),categoryMC_.c_str(),cutValue_,xLow_,xHigh_));
    }
    else {
      fullTreeSgn->Draw("mass>>hSP",Form("tag_PUMCWeight2012A*(%s && %s>=%f && mass>%f && mass<%f && mcTrue && tag_GenDecay==23*11 && HpsLooseCombIsoDBCorr>0.5 && DecayMode>0.5 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),categoryMC_.c_str(),cutValue_,xLow_,xHigh_));
    }
  }
  
  if(PUWeight && SecondEleVeto){
    
    fullTreeSgn->Draw("mass>>hS",Form("tag_PUMCWeight2012A*(%s && mass>%f && mass<%f && mcTrue && MatchElePassVeto<0.5 && tag_GenDecay==23*11 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),xLow_,xHigh_));
    
    if(MVAIso){
      fullTreeSgn->Draw("mass>>hSP",Form("tag_PUMCWeight2012A*(%s && %s>=%f && mass>%f && mass<%f && mcTrue && MatchElePassVeto<0.5 && tag_GenDecay==23*11 && HpsIsoMVALoose>0.5 && DecayMode>0.5 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),categoryMC_.c_str(),cutValue_,xLow_,xHigh_));
    }
    else {
      fullTreeSgn->Draw("mass>>hSP",Form("tag_PUMCWeight2012A*(%s && %s>=%f && mass>%f && mass<%f && mcTrue && MatchElePassVeto<0.5 && tag_GenDecay==23*11 && HpsLooseCombIsoDBCorr>0.5 && DecayMode>0.5 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match)",bin_.c_str(),categoryMC_.c_str(),cutValue_,xLow_,xHigh_));
    }
  }
    
  float SGNtrue = hS->Integral();
  float SGNtruePass = hSP->Integral();

  float McTruthEff    = SGNtruePass/SGNtrue;
  float BinomialError = TMath::Sqrt(SGNtruePass/SGNtrue*(1-SGNtruePass/SGNtrue)/SGNtrue);
  
  cout << bin_.c_str() << " ==> MCTRUTH: " << McTruthEff << " +/- " << BinomialError << endl;

  //   return;
//   int pause ; 
  //    cin>>pause;

  delete hS; delete hSP;
 
  float reductionFactor = 1.0;

  //Trees of passing and failing for signal MC without pu reweight
  TTree* fullTreeSgnCutP = new TTree;
  TTree* fullTreeSgnCutF = new TTree;

  fullTreeSgnCutF = fullTreeSgn->CopyTree( Form("(%s< %f && %s && pair_tnpCharge==0 && event_met_pfmet<25 && mcTrue && tag_GenDecay==23*11 && tag_PfRelIso<0.1 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ,"",int(fullTreeSgn->GetEntries()*reductionFactor));

  if(MVAIso){
    fullTreeSgnCutP = fullTreeSgn->CopyTree( Form("(%s>=%f && %s && pair_tnpCharge==0 && event_met_pfmet<25 && mcTrue && tag_GenDecay==23*11 && HpsIsoMVALoose>0.5 && DecayMode>0.5 &&tag_PfRelIso<0.1 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ,"",int(fullTreeSgn->GetEntries()*reductionFactor));
  }
  else {
    fullTreeSgnCutP = fullTreeSgn->CopyTree( Form("(%s< %f && %s && pair_tnpCharge==0 && event_met_pfmet<25 && mcTrue && tag_GenDecay==23*11 && HpsLooseCombIsoDBCorr>0.5 && DecayMode>0.5 && tag_PfRelIso<0.1 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ,"",int(fullTreeSgn->GetEntries()*reductionFactor));
  }

  if(SecondEleVeto){
    fullTreeSgnCutF = fullTreeSgn->CopyTree( Form("(%s< %f && %s && pair_tnpCharge==0 && event_met_pfmet<25 && mcTrue && MatchElePassVeto<0.5 && tag_GenDecay==23*11 && tag_PfRelIso<0.1 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ,"",int(fullTreeSgn->GetEntries()*reductionFactor));
    if(MVAIso){
      fullTreeSgnCutP = fullTreeSgn->CopyTree( Form("(%s>=%f && %s && pair_tnpCharge==0 && event_met_pfmet<25 && mcTrue && MatchElePassVeto<0.5 && tag_GenDecay==23*11 && HpsIsoMVALoose>0.5 && DecayMode>0.5 &&tag_PfRelIso<0.1 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ,"",int(fullTreeSgn->GetEntries()*reductionFactor));
    }
    else {
      fullTreeSgnCutP = fullTreeSgn->CopyTree( Form("(%s< %f && %s && pair_tnpCharge==0 && event_met_pfmet<25 && mcTrue && MatchElePassVeto<0.5 && tag_GenDecay==23*11 && HpsLooseCombIsoDBCorr>0.5 && DecayMode>0.5 && tag_PfRelIso<0.1 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ,"",int(fullTreeSgn->GetEntries()*reductionFactor));
    }    
  }

  /******************/
  /******************/
  /******************/
  /*ROOFIT VARIABLES*/
  /******************/
  /******************/
  /******************/

  RooRealVar mass("mass","m_{tp} (GeV/c^{2})",xLow_,xHigh_);
  mass.setBins( 10000, "fft" );
  mass.setBins( nBins_ );

  TTree* fullTreeDataCutF = fullTreeData->CopyTree( Form("(%s <%f && %s   && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_HLTEle27Match && tag_Pt>25)",categoryDATA_.c_str(),cutValue_,bin_.c_str()) );

  TTree* fullTreeDataCutP;

  if(MVAIso){
    fullTreeDataCutP = fullTreeData->CopyTree( Form("(%s>=%f && %s  && HpsIsoMVALoose>0.5 && DecayMode>0.5 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_HLTEle27Match && tag_Pt>25)",categoryDATA_.c_str(),cutValue_,bin_.c_str()) ); 
  }
  else{
    fullTreeDataCutP = fullTreeData->CopyTree( Form("(%s>=%f && %s  && HpsLooseCombIsoDBCorr>0.5 && DecayMode>0.5 && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_HLTEle27Match && tag_Pt>25)",categoryDATA_.c_str(),cutValue_,bin_.c_str()) ); 
  }

  //Roofit datasets of data passing and failing cut
  mass.setBins(nBins_);
  RooDataSet DataDataSetP("DataDataSetP","dataset for Data pass", RooArgSet(mass), Import( *fullTreeDataCutP ) );
  std::cout << "data dataset Pass " << DataDataSetP.numEntries() << "  " << std::endl;
  RooDataHist DataDataHistP("DataDataHistP","",RooArgSet(mass),DataDataSetP, 1.0);
//   float nPass = DataDataHistP.sum(false);
  RooDataSet DataDataSetF("DataDataSetF","dataset for Data fail", RooArgSet(mass), Import( *fullTreeDataCutF ) );
  std::cout << "data dataset Fail " << DataDataSetF.numEntries() << "  " << std::endl;
  RooDataHist DataDataHistF("DataDataHistF","",RooArgSet(mass),DataDataSetF, 1.0);
//   float nFail = DataDataHistF.sum(false);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Construction of datasets from simulation, closure test//////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TTree* fullTreeSSDataCutP  = new TTree;
  TTree* fullTreeSSZttCutP   = new TTree;
  TTree* fullTreeSSZeeCutP   = new TTree;
  TTree* fullTreeSSWenCutP   = new TTree;
  TTree* fullTreeSSTTbCutP   = new TTree;
  TTree* fullTreeSSDataCutF  = new TTree;
  TTree* fullTreeSSZttCutF   = new TTree;
  TTree* fullTreeSSZeeCutF   = new TTree;
  TTree* fullTreeSSWenCutF   = new TTree;
  TTree* fullTreeSSTTbCutF   = new TTree;

  TTree* fullTreeZttCutP   = new TTree;
  TTree* fullTreeZeeCutP   = new TTree;
  TTree* fullTreeWenCutP   = new TTree;
  TTree* fullTreeTTbCutP   = new TTree;
  TTree* fullTreeZttCutF   = new TTree;
  TTree* fullTreeZeeCutF   = new TTree;
  TTree* fullTreeWenCutF   = new TTree;
  TTree* fullTreeTTbCutF   = new TTree;

  //SS
  fullTreeSSDataCutP = fullTreeData->CopyTree( Form("(%s>=%f && %s  && tag_PfRelIso<0.1 && (pair_tnpCharge==-2 || pair_tnpCharge==-2) && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ); 
    
  fullTreeSSZttCutP = fullTreeZtt->CopyTree( Form("(%s>=%f && %s  && tag_PfRelIso<0.1 && (pair_tnpCharge==-2 || pair_tnpCharge==-2) && event_met_pfmet<25 && tag_Pt>25 && tag_GenDecay==23*15  && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ); 

  fullTreeSSZeeCutP = fullTreeZtt->CopyTree( Form("(%s>=%f && %s  && tag_PfRelIso<0.1 && (pair_tnpCharge==-2 || pair_tnpCharge==-2) && event_met_pfmet<25 && tag_Pt>25 && tag_GenDecay==23*11  && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ); 

  fullTreeSSWenCutP = fullTreeWen->CopyTree( Form("(%s>=%f && %s  && tag_PfRelIso<0.1 && (pair_tnpCharge==-2 || pair_tnpCharge==-2) && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(), xLow_, xHigh_) ); 
  
  fullTreeSSTTbCutP = fullTreeTTb->CopyTree( Form("(%s>=%f && %s  && tag_PfRelIso<0.1 && (pair_tnpCharge==-2 || pair_tnpCharge==-2) && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(), xLow_, xHigh_) ); 
    
  fullTreeSSDataCutF = fullTreeData->CopyTree( Form("(%s< %f && %s  && tag_PfRelIso<0.1 && (pair_tnpCharge==-2 || pair_tnpCharge==-2) && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ); 
  
  fullTreeSSZttCutF = fullTreeZtt->CopyTree( Form("(%s< %f && %s  && tag_PfRelIso<0.1 && (pair_tnpCharge==-2 || pair_tnpCharge==-2) && event_met_pfmet<25 && tag_Pt>25 && tag_GenDecay==23*15  && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ); 

  fullTreeSSZeeCutF = fullTreeZtt->CopyTree( Form("(%s< %f && %s  && tag_PfRelIso<0.1 && (pair_tnpCharge==-2 || pair_tnpCharge==-2) && event_met_pfmet<25 && tag_Pt>25 && tag_GenDecay==23*11  && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ); 
  
  fullTreeSSWenCutF = fullTreeWen->CopyTree( Form("(%s< %f && %s  && tag_PfRelIso<0.1 && (pair_tnpCharge==-2 || pair_tnpCharge==-2) && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ); 
  
  fullTreeSSTTbCutF = fullTreeTTb->CopyTree( Form("(%s< %f && %s  && tag_PfRelIso<0.1 && (pair_tnpCharge==-2 || pair_tnpCharge==-2) && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ); 
  
  //OS
  fullTreeZttCutP = fullTreeZtt->CopyTree( Form("(%s>=%f && %s  && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_GenDecay==23*15  && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ); 

  fullTreeZeeCutP = fullTreeZtt->CopyTree( Form("(%s>=%f && %s  && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_GenDecay==23*11  && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ); 

  fullTreeWenCutP = fullTreeWen->CopyTree( Form("(%s>=%f && %s  && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(), xLow_, xHigh_) ); 
  
  fullTreeTTbCutP = fullTreeTTb->CopyTree( Form("(%s>=%f && %s  && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(), xLow_, xHigh_) ); 
    
  
  fullTreeZttCutF = fullTreeZtt->CopyTree( Form("(%s< %f && %s  && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_GenDecay==23*15  && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ); 
  
  fullTreeZeeCutF = fullTreeZtt->CopyTree( Form("(%s< %f && %s  && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_GenDecay==23*11  && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ); 

  fullTreeWenCutF = fullTreeWen->CopyTree( Form("(%s< %f && %s  && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ); 
  
  fullTreeTTbCutF = fullTreeTTb->CopyTree( Form("(%s< %f && %s  && tag_PfRelIso<0.1 && pair_tnpCharge==0 && event_met_pfmet<25 && tag_Pt>25 && tag_HLTEle27Match && mass>%f && mass<%f)",categoryMC_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_) ); 
  
  TH1F* hEta = new TH1F("hEta","",1,-10,10); 
  fullTreeSSDataCutP->Draw("Eta>>hEta");
//   float NumSSP = hEta->Integral();
  hEta->Reset();
  fullTreeSSZttCutP->Draw("Eta>>hEta","tag_PUMCWeight2012A");
  float NumSSZttP = hEta->Integral();
  hEta->Reset();
  fullTreeSSZeeCutP->Draw("Eta>>hEta","tag_PUMCWeight2012A");
  float NumSSZeeP = hEta->Integral();
  hEta->Reset();
  fullTreeSSWenCutP->Draw("Eta>>hEta","tag_PUMCWeight2012A");
  float NumSSWenP = hEta->Integral();
  hEta->Reset();
  fullTreeSSTTbCutP->Draw("Eta>>hEta","tag_PUMCWeight2012A");
  float NumSSTTbP = hEta->Integral();
  hEta->Reset();
  
  fullTreeSSDataCutF->Draw("Eta>>hEta");
//   float NumSSF = hEta->Integral();
  hEta->Reset();
  fullTreeSSZttCutF->Draw("Eta>>hEta","tag_PUMCWeight2012A");
  float NumSSZttF = hEta->Integral();
  hEta->Reset();
  fullTreeSSZeeCutF->Draw("Eta>>hEta","tag_PUMCWeight2012A");
  float NumSSZeeF = hEta->Integral();
  hEta->Reset();
  fullTreeSSWenCutF->Draw("Eta>>hEta","tag_PUMCWeight2012A");
  float NumSSWenF = hEta->Integral();
  hEta->Reset();
  fullTreeSSTTbCutF->Draw("Eta>>hEta","tag_PUMCWeight2012A");
  float NumSSTTbF = hEta->Integral();
  hEta->Reset();

  fullTreeDataCutP->Draw("Eta>>hEta");
//   float NumP = hEta->Integral();
  hEta->Reset();
  fullTreeZttCutP->Draw("Eta>>hEta","tag_PUMCWeight2012A");
  float NumZttP = hEta->Integral();
  hEta->Reset();
  fullTreeZeeCutP->Draw("Eta>>hEta","tag_PUMCWeight2012A");
//   float NumZeeP = hEta->Integral();
  hEta->Reset();
  fullTreeWenCutP->Draw("Eta>>hEta","tag_PUMCWeight2012A");
  float NumWenP = hEta->Integral();
  hEta->Reset();
  fullTreeTTbCutP->Draw("Eta>>hEta","tag_PUMCWeight2012A");
  float NumTTbP = hEta->Integral();
  hEta->Reset();

  fullTreeDataCutF->Draw("Eta>>hEta");
//   float NumF = hEta->Integral();
  hEta->Reset();
  fullTreeZttCutF->Draw("Eta>>hEta","tag_PUMCWeight2012A");
  float NumZttF = hEta->Integral();
  hEta->Reset();
  fullTreeZeeCutF->Draw("Eta>>hEta","tag_PUMCWeight2012A");
//   float NumZeeF = hEta->Integral();
  hEta->Reset();
  fullTreeWenCutF->Draw("Eta>>hEta","tag_PUMCWeight2012A");
  float NumWenF = hEta->Integral();
  hEta->Reset();
  fullTreeTTbCutF->Draw("Eta>>hEta","tag_PUMCWeight2012A");
  float NumTTbF = hEta->Integral();
  hEta->Reset();
  if (DEBUG)cout<<"QCD estimation :"<<endl;

  //Pass
  RooDataSet SSDataDataSetP ("SSDataDataSetP","", RooArgSet(mass), Import( *fullTreeSSDataCutP));
  RooDataHist SSDataDataHistP("SSDataDataHistP","",RooArgSet(mass),SSDataDataSetP, 1.0);
  RooDataSet SSZttDataSetP ("SSZttDataSetP","", RooArgSet(mass), Import( *fullTreeSSZttCutP));
  SSZttDataSetP.reduce(EventRange(1,(int)( NumSSZttP*mcLumi_/(readEventsZtt*reductionFactor/SigmaZLL))));
  RooDataHist SSZttDataHistP ("SSZttDataHistP","",RooArgSet(mass),SSZttDataSetP, -1.0);
  RooDataSet SSZeeDataSetP ("SSZeeDataSetP","", RooArgSet(mass), Import( *fullTreeSSZeeCutP));
  SSZeeDataSetP.reduce(EventRange(1,(int)( NumSSZeeP*mcLumi_/(readEventsZtt*reductionFactor/SigmaZLL))));
  RooDataHist SSZeeDataHistP ("SSZeeDataHistP","",RooArgSet(mass),SSZeeDataSetP, -1.0);
  RooDataSet SSWenDataSetP ("SSWenDataSetP","", RooArgSet(mass), Import( *fullTreeSSWenCutP));
  SSWenDataSetP.reduce(EventRange(1,(int)( NumSSWenP*mcLumi_/(readEventsWen*reductionFactor/SigmaWen))));
  RooDataHist SSWenDataHistP ("SSWenDataHistP","",RooArgSet(mass),SSWenDataSetP, -1.0);
  RooDataSet SSTTbDataSetP ("SSTTbDataSetP","", RooArgSet(mass), Import( *fullTreeSSTTbCutP));
  SSTTbDataSetP.reduce(EventRange(1,(int)( NumSSTTbP*mcLumi_/(readEventsTTb*reductionFactor/SigmaTTb))));
  RooDataHist SSTTbDataHistP ("SSTTbDataHistP","",RooArgSet(mass),SSTTbDataSetP, -1.0);

//   float SSDataP = SSDataDataSetP.sumEntries();
//   float SSZttP = SSZttDataSetP.sumEntries();
//   float SSWenP = SSWenDataSetP.sumEntries();
//   float SSTTbP = SSTTbDataSetP.sumEntries();

  SSDataDataHistP.add(SSZttDataHistP) ;
  SSDataDataHistP.add(SSZeeDataHistP) ;
  SSDataDataHistP.add(SSWenDataHistP) ;
  SSDataDataHistP.add(SSTTbDataHistP) ;
  RooHistPdf QCDHistPdfP("QCDHistPdfP", "", RooArgSet(mass), SSDataDataHistP, 0);

//   float QCDP = SSDataDataSetP.sumEntries();

  //Fail
  RooDataSet SSDataDataSetF ("SSDataDataSetF","", RooArgSet(mass), Import( *fullTreeSSDataCutF ) );
  RooDataHist SSDataDataHistF("SSDataDataHistF","",RooArgSet(mass),SSDataDataSetF, 1.0);
  RooDataSet SSZttDataSetF ("SSZttDataSetF","", RooArgSet(mass), Import( *fullTreeSSZttCutF ));
  SSZttDataSetF.reduce(EventRange(1,(int)( NumSSZttF*mcLumi_/(readEventsZtt*reductionFactor/SigmaZLL))));
  RooDataHist SSZttDataHistF ("SSZttDataHistF","",RooArgSet(mass),SSZttDataSetF, -1.0);
  RooDataSet SSZeeDataSetF ("SSZeeDataSetF","", RooArgSet(mass), Import( *fullTreeSSZeeCutF ));
  SSZeeDataSetF.reduce(EventRange(1,(int)( NumSSZeeF*mcLumi_/(readEventsZtt*reductionFactor/SigmaZLL))));
  RooDataHist SSZeeDataHistF ("SSZeeDataHistF","",RooArgSet(mass),SSZeeDataSetF, -1.0);
  RooDataSet SSWenDataSetF ("SSWenDataSetF","", RooArgSet(mass), Import( *fullTreeSSWenCutF ));
  SSWenDataSetF.reduce(EventRange(1,(int)( NumSSWenF*mcLumi_/(readEventsWen*reductionFactor/SigmaWen))));
  RooDataHist SSWenDataHistF ("SSWenDataHistF","",RooArgSet(mass),SSWenDataSetF, -1.0);
  RooDataSet SSTTbDataSetF ("SSTTbDataSetF","", RooArgSet(mass), Import( *fullTreeSSTTbCutF ));
  SSTTbDataSetF.reduce(EventRange(1,(int)( NumSSTTbF*mcLumi_/(readEventsTTb*reductionFactor/SigmaTTb))));
  RooDataHist SSTTbDataHistF ("SSTTbDataHistF","",RooArgSet(mass),SSTTbDataSetF, -1.0);

//   float SSDataF = SSDataDataSetF.sumEntries();
//   float SSZttF = SSZttDataSetF.sumEntries();
//   float SSWenF = SSWenDataSetF.sumEntries();
//   float SSTTbF = SSTTbDataSetF.sumEntries();

  SSDataDataHistF.add(SSZttDataHistF) ;
  SSDataDataHistF.add(SSZeeDataHistF) ;
  SSDataDataHistF.add(SSWenDataHistF) ;
  SSDataDataHistF.add(SSTTbDataHistF) ;
  RooHistPdf QCDHistPdfF("QCDHistPdfF", "", RooArgSet(mass), SSDataDataHistF, 0);

// C r e a t e   c o n s t r a i n t   p d f 
  // -----------------------------------------
  RooRealVar f("f","f",0.5,0.,1.) ;
  // Construct Gaussian constraint p.d.f on parameter f at 1.05 with resolution of 0.1
  RooGaussian fconstraint("fconstraint","fconstraint",f,RooConst(1.05),RooConst(0.1)) ;



  // M E T H O D   1   -   A d d   i n t e r n a l   c o n s t r a i n t   t o   m o d e l 
  // -------------------------------------------------------------------------------------

  // Multiply constraint term with regular p.d.f using RooProdPdf
  // Specify in fitTo() that internal constraints on parameter f should be used

  // Multiply constraint with p.d.f
  RooProdPdf QCDHistPdfF_C("QCDHistPdfF_C","QCD model with constraint",RooArgSet(QCDHistPdfF,fconstraint)) ;
  RooProdPdf QCDHistPdfP_C("QCDHistPdfP_C","QCD model with constraint",RooArgSet(QCDHistPdfP,fconstraint)) ;

//   float QCDF = SSDataDataSetF.sumEntries();

    //Unweighted MC efficiency
  float unWeightedMCEff = float(fullTreeSgnCutP->GetEntries())/(fullTreeSgnCutF->GetEntries()+fullTreeSgnCutP->GetEntries());
  //The unweighted efficiency should match the PU Weighted efficiency so we append passing signal if unWeff<Weff and failing signal if unWeff>Weff:
  int sgnPToAppend =  unWeightedMCEff<McTruthEff ? int(fullTreeSgnCutP->GetEntries()) : int(McTruthEff/(1-McTruthEff)*fullTreeSgnCutF->GetEntries()) ;
  int sgnFToAppend =  unWeightedMCEff>McTruthEff ? int(fullTreeSgnCutF->GetEntries()) : int((1-McTruthEff)/McTruthEff*fullTreeSgnCutP->GetEntries()) ; 

  cout << "MCTruthEff = " << McTruthEff << ", unWeightedMCEff = " << unWeightedMCEff << " %%%%%%%  sgnToAppend PASS " << sgnPToAppend << " -- sgnToAppend FAIL " << sgnFToAppend << endl;

  RooDataSet ZeeDataSetP("ZeeDataSetP","", RooArgSet(mass), Import( *fullTreeZeeCutP ) );
  ZeeDataSetP.reduce(EventRange(1,(int)( sgnPToAppend*mcLumi_/(readEventsZtt*reductionFactor/SigmaZLL))));
  RooDataHist ZeeDataHistP("ZeeDataHistP","",RooArgSet(mass),ZeeDataSetP, 1.0);
  RooHistPdf ZeeHistPdfP("ZeeHistPdfP", "", RooArgSet(mass), ZeeDataHistP, 0);

  RooDataSet ZttDataSetP("ZttDataSetP","", RooArgSet(mass), Import( *fullTreeZttCutP ) );
  ZttDataSetP.reduce(EventRange(1,(int)( NumZttP*mcLumi_/(readEventsZtt*reductionFactor/SigmaZLL))));
  RooDataHist ZttDataHistP("ZttDataHistP","",RooArgSet(mass),ZttDataSetP, 1.0);
  RooHistPdf ZttHistPdfP("ZttHistPdfP", "", RooArgSet(mass), ZttDataHistP, 0);

  RooDataSet WenDataSetP("WenDataSetP","", RooArgSet(mass), Import( *fullTreeWenCutP ) );
  WenDataSetP.reduce(EventRange(1,(int)( NumWenP*mcLumi_/(readEventsWen*reductionFactor/SigmaWen))));
  RooDataHist WenDataHistP("WenDataHistP","",RooArgSet(mass),WenDataSetP, 1.0);
  RooHistPdf WenHistPdfP("WenHistPdfP", "", RooArgSet(mass), WenDataHistP, 0);

  RooDataSet TTbDataSetP("TTbDataSetP","", RooArgSet(mass), Import( *fullTreeTTbCutP ) );
  TTbDataSetP.reduce(EventRange(1,(int)( NumTTbP*mcLumi_/(readEventsTTb*reductionFactor/SigmaTTb))));
  RooDataHist TTbDataHistP("TTbDataHistP","",RooArgSet(mass),TTbDataSetP, 1.0);
  RooHistPdf TTbHistPdfP("TTbHistPdfP", "", RooArgSet(mass), TTbDataHistP, 0);

  DataDataHistP.reset();
  DataDataHistP.add(((RooDataHist)SSDataDataHistP));
  DataDataHistP.add(((RooDataHist)ZttDataHistP));
  DataDataHistP.add(((RooDataHist)ZeeDataHistP));
  DataDataHistP.add(((RooDataHist)WenDataHistP));
  DataDataHistP.add(((RooDataHist)TTbDataHistP));
//   float nPassTemplate =  DataDataHistP.sum(false);

  RooDataSet ZeeDataSetF("ZeeDataSetF","", RooArgSet(mass), Import( *fullTreeZeeCutF ) );
  ZeeDataSetF.reduce(EventRange(1,(int)( sgnFToAppend*mcLumi_/(readEventsZtt*reductionFactor/SigmaZLL))));
  RooDataHist ZeeDataHistF("ZeeDataHistF","",RooArgSet(mass),ZeeDataSetF, 1.0);
  RooHistPdf ZeeHistPdfF("ZeeHistPdfF", "", RooArgSet(mass), ZeeDataHistF, 0);

  RooDataSet ZttDataSetF("ZttDataSetF","", RooArgSet(mass), Import( *fullTreeZttCutF ) );
  ZttDataSetF.reduce(EventRange(1,(int)( NumZttF*mcLumi_/(readEventsZtt*reductionFactor/SigmaZLL))));
  RooDataHist ZttDataHistF("ZttDataHistF","",RooArgSet(mass),ZttDataSetF, 1.0);
  RooHistPdf ZttHistPdfF("ZttHistPdfF", "", RooArgSet(mass), ZttDataHistF, 0);

  RooDataSet WenDataSetF("WenDataSetF","", RooArgSet(mass), Import( *fullTreeWenCutF ) );
  WenDataSetF.reduce(EventRange(1,(int)( NumWenF*mcLumi_/(readEventsWen*reductionFactor/SigmaWen))));
  RooDataHist WenDataHistF("WenDataHistF","",RooArgSet(mass),WenDataSetF, 1.0);
  RooHistPdf WenHistPdfF("WenHistPdfF", "", RooArgSet(mass), WenDataHistF, 0);

  RooDataSet TTbDataSetF("TTbDataSetF","", RooArgSet(mass), Import( *fullTreeTTbCutF ) );
  TTbDataSetF.reduce(EventRange(1,(int)( NumTTbF*mcLumi_/(readEventsTTb*reductionFactor/SigmaTTb))));
  RooDataHist TTbDataHistF("TTbDataHistF","",RooArgSet(mass),TTbDataSetF, 1.0);
  RooHistPdf TTbHistPdfF("TTbHistPdfF", "", RooArgSet(mass), TTbDataHistF, 0);

  DataDataHistF.reset();
  DataDataHistF.add(((RooDataHist)SSDataDataHistF));
  DataDataHistF.add(((RooDataHist)ZttDataHistF));
  DataDataHistF.add(((RooDataHist)ZeeDataHistF));
  DataDataHistF.add(((RooDataHist)WenDataHistF));
  DataDataHistF.add(((RooDataHist)TTbDataHistF));
//   float nFailTemplate =  DataDataHistF.sum(false);

  RooRealVar DataNumSgn("DataNumSgn","",  0,10000000);
  RooRealVar DataEfficiency("DataEfficiency","",0.5,0,1);

  RooFormulaVar DataNumSgnP("DataNumSgnP","DataEfficiency*DataNumSgn",    RooArgSet(DataEfficiency,DataNumSgn));
  RooFormulaVar DataNumSgnF("DataNumSgnF","(1-DataEfficiency)*DataNumSgn",RooArgSet(DataEfficiency,DataNumSgn));
  
  //Create pdfs for passing
  RooRealVar CoeffZeeP("CoeffZeeP","",0,10000000);
  RooRealVar CoeffZttP("CoeffZttP","",0,10000000);
  RooRealVar CoeffWenP("CoeffWenP","",0,10000000);
  RooRealVar CoeffTTbP("CoeffTTbP","",0,10000000);
  RooRealVar CoeffQCDP("CoeffQCDP","",0,10000000);

  RooAddPdf DataModelP("DataModelP", "", RooArgList(ZeeHistPdfP,ZttHistPdfP,WenHistPdfP,TTbHistPdfP,QCDHistPdfP_C), RooArgList(CoeffZeeP/*DataNumSgnP*/,CoeffZttP,CoeffWenP,CoeffTTbP,CoeffQCDP));

  //Create pdfs for failing
  RooRealVar CoeffZeeF("CoeffZeeF","",0,10000000);
  RooRealVar CoeffZttF("CoeffZttF","",0,10000000);
  RooRealVar CoeffWenF("CoeffWenF","",0,10000000);
  RooRealVar CoeffTTbF("CoeffTTbF","",0,10000000);
  RooRealVar CoeffQCDF("CoeffQCDF","",0,10000000);

  RooAddPdf DataModelF("DataModelF", "", RooArgList(ZeeHistPdfF,ZttHistPdfF,WenHistPdfF,TTbHistPdfF,QCDHistPdfF_C), RooArgList(CoeffZeeF/*DataNumSgnF*/,CoeffZttF,CoeffWenF,CoeffTTbF,CoeffQCDF));
//   //Pdf models for passing and failing events extracted from before
//   RooAddPdf DataModelP("DataModelP","",RooArgList(/*cbSgnP_C*/ sgnPdfP /*sgnTemplatePdfP*/,bkgPdfP),RooArgList(DataNumSgnP,DataNumBkgP));
//   RooAddPdf DataModelF("DataModelF","",RooArgList(/*cbSgnF_C*/ sgnPdfF /*sgnTemplatePdfF*/,bkgPdfF),RooArgList(DataNumSgnF,DataNumBkgF));
      
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Construction of datasets from simulation, closure test//////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  mass.setBins(nBins_);

  /******************/
  /******************/
  /******************/
  /*FIT             */
  /******************/
  /******************/
  /******************/

  RooCategory category("category","category") ;
  category.defineType("pass") ;
  category.defineType("fail") ;
//   //Variables 
//   RooRealVar DataNumBkgF("DataNumBkgF","",0,10000000);
//   RooRealVar DataNumBkgP("DataNumBkgP","",0,10000000);
//   RooRealVar DataNumSgn("DataNumSgn","",  0,10000000);
//   RooRealVar DataEfficiency("DataEfficiency","",0.5,0,1);

//   RooFormulaVar DataNumSgnP("DataNumSgnP","DataEfficiency*DataNumSgn",    RooArgSet(DataEfficiency,DataNumSgn));
//   RooFormulaVar DataNumSgnF("DataNumSgnF","(1-DataEfficiency)*DataNumSgn",RooArgSet(DataEfficiency,DataNumSgn));
  
//   RooRealVar DataNumSgnP_("DataNumSgnP_","",0,10000);
//   RooRealVar DataNumSgnF_("DataNumSgnF_","",0,10000);
//   RooAddPdf DataModelP_("DataModelP_","",RooArgList(sgnPdfP,bkgPdfP),RooArgList(DataNumSgnP_,DataNumBkgP));
//   DataModelP.fitTo(DataDataHistP, Extended(1), Minos(1), Save(1), NumCPU(4),SumW2Error(1) /*,ExternalConstraints( RooArgSet(meanSgn_CPdf,widthSgn_CPdf) )*/);
  
 //  RooAddPdf DataModelF_("DataModelF_","",RooArgList(sgnPdfF,bkgPdfF),RooArgList(DataNumSgnF_,DataNumBkgF));
//   DataModelF.fitTo(DataDataHistF, Extended(1), Minos(1), Save(1), NumCPU(4),SumW2Error(1) /*,ExternalConstraints( RooArgSet(meanSgn_CPdf,widthSgn_CPdf) )*/);
  
//   RooPlot* frame2 = mass.frame(Title("template"));
//   DataDataHistP.plotOn(frame2);
//   DataModelP_.plotOn(frame2, LineColor(kBlue), LineStyle(kSolid));
//   DataModelP_.plotOn(frame2, Components("sgnPdfP"), LineColor(kRed), LineStyle(kSolid));
//   DataModelP_.plotOn(frame2, Components("bkgPdfP"), LineColor(kGreen), LineStyle(kSolid));
//   frame2->Draw();
  
  //return;
  
  
  

//   // binned combined dataset
//   RooDataHist DataCombData("DataCombData","combined data",mass,Index(category),Import("pass", *(DataDataSetP.createHistogram("histoDataP",mass)) ) ,Import("fail", *(DataDataSetF.createHistogram("histoDataF",mass))), Weight(0.5) ) ;

  // binned combined dataset
  RooDataHist DataCombData("DataCombData","combined data",mass,Index(category),Import("pass", *(DataDataHistP.createHistogram("histoDataP",mass)) ) ,Import("fail", *(DataDataSetF.createHistogram("histoDataF",mass))), Weight(0.5) ) ;
  std::cout << "data dataHist Comb " << DataCombData.sumEntries() << "  " << std::endl;
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

  RooSimultaneous DataSimPdf("DataSimPdf","simultaneous pdf",category) ;
  DataSimPdf.addPdf(DataModelP,"pass") ;
  DataSimPdf.addPdf(DataModelF,"fail") ;
  
  //mass.setBins( 10000, "fft" );
  mass.setBins( nBins_ );
  RooFitResult* ResDataCombinedFit =  0;
  //Binned fit
//   ResDataCombinedFit = DataSimPdf.fitTo(DataCombData , Extended(1), Minos(1), Save(1), NumCPU(4), /*ExternalConstraints( RooArgSet(alfaSgn_CPdf,nSgn_CPdf) )*/  SumW2Error(1));
  ResDataCombinedFit = DataSimPdf.fitTo(DataCombData , Extended(1), Minos(1), Save(1), NumCPU(8),SumW2Error(1));

//   //Get the fitted efficiency and number of signal events
  RooArgSet DataFitParam(ResDataCombinedFit->floatParsFinal());
//   RooRealVar* DataEffFit      = (RooRealVar*)(&DataFitParam["DataEfficiency"]);
 
//  //RooRealVar* DataNumSigFit   = (RooRealVar*)(&DataFitParam["DataNumSgn"]);
  
  string theSample = doMC_ ? "Simulation" : "Data" ;
  float theLumi = doMC_ ? mcLumi_ : 620;
  
  RooPlot* DataFrameP = mass.frame(Bins(40),Title(Form("CMS Preliminary 2012B  #sqrt{s}=8 TeV %s  L=%.0f pb^{-1}:  passing probe",theSample.c_str(),theLumi)));
  DataCombData.plotOn(DataFrameP,Cut("category==category::pass"),Name("dataP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), LineColor(kBlack),Name("modelP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("ZeeHistPdfP"), LineColor(kBlue), LineStyle(kSolid),Name("signal onlyP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("ZttHistPdfP"), LineColor(kRed), LineStyle(kSolid),Name("ZttP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("WenHistPdfP"), LineColor(kGreen), LineStyle(kSolid),Name("WenP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("TTbHistPdfP"), LineColor(kMagenta), LineStyle(kSolid),Name("TTbP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("QCDHistPdfP"), LineColor(kOrange), LineStyle(kSolid),Name("QCDP"));
//   DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("sgnPdfP"), LineColor(kRed), LineStyle(kSolid),Name("signal onlyP"));
//   DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("bkgPdfP"), LineColor(kMagenta), LineStyle(kDashed),Name("backgroundP"));
  DataFrameP->SetTitleOffset(1,"Y");
  DataFrameP->SetTitleSize(1,"Y");

  RooPlot* DataFrameF = mass.frame(Bins(40),Title(Form("CMS Preliminary 2012B  #sqrt{s}=8 TeV %s  L=%.0f pb^{-1}:  failing probe",theSample.c_str(),theLumi)));
  DataCombData.plotOn(DataFrameF,Cut("category==category::fail"),Name("dataF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), LineColor(kBlack),Name("modelF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("ZeeHistPdfF"), LineColor(kBlue), LineStyle(kSolid),Name("signal onlyF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("ZttHistPdfF"), LineColor(kRed), LineStyle(kSolid),Name("ZttF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("WenHistPdfF"), LineColor(kGreen), LineStyle(kSolid),Name("WenF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("TTbHistPdfF"), LineColor(kMagenta), LineStyle(kSolid),Name("TTbF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("QCDHistPdfF"), LineColor(kOrange), LineStyle(kSolid),Name("QCDF"));

//   DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("sgnPdfF"), LineColor(kRed), LineStyle(kSolid),Name("signal onlyF"));
//   DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("bkgPdfF"), LineColor(kMagenta), LineStyle(kDashed),Name("backgroundF"));
  DataFrameF->SetTitleOffset(1,"Y");
  DataFrameF->SetTitleSize(1,"Y");

  TCanvas *cPass = new TCanvas("fitCanvasP","canvas",10,30,650,600);
  cPass->SetGrid(0,0);
  cPass->SetFillStyle(4000);
  cPass->SetFillColor(10);
  cPass->SetTicky();
  cPass->SetObjectStat(0);
      
  cPass->cd();
  DataFrameP->Draw();
  TLegend *leg1 = new TLegend(0.65,0.73,0.86,0.87);
  leg1->SetFillColor(kWhite);
  leg1->SetLineColor(kWhite);
  leg1->AddEntry("dataP","Data", "P");
  leg1->AddEntry("modelP","Signal + background","L");
//   leg1->AddEntry("backgroundP","Background only", "L");
  leg1->AddEntry("signal onlyP","Signal only", "L");
  leg1->AddEntry("ZttP","Ztt", "L");
  leg1->AddEntry("WenP","Wen", "L");
  leg1->AddEntry("TTbP","TTb", "L");
  leg1->AddEntry("QCDP","QCD", "L");
  leg1->Draw();
  //return;

  string fileNameP = "Pass_"+tnp_+"_"+categoryDATA_;
  if(doMC_) fileNameP = fileNameP+"_MC";
  if (!doMC_) fileNameP = fileNameP+"_DATA";
  if(MVAIso)cPass->SaveAs(Form("./ElecTauTnP/plots/computeElecTauTnPMacro_ClosureFit_MVAIso_%s_%s.pdf",fileNameP.c_str(),binLabel_.c_str()));
  else cPass->SaveAs(Form("./ElecTauTnP/plots/computeElecTauTnPMacro_ClosureFit_%s_%s.pdf",fileNameP.c_str(),binLabel_.c_str()));

  TCanvas *cFail = new TCanvas("fitCanvasF","canvas",10,30,650,600);
  cFail->SetGrid(0,0);
  cFail->SetFillStyle(4000);
  cFail->SetFillColor(10);
  cFail->SetTicky();
  cFail->SetObjectStat(0);

  cFail->cd();
  DataFrameF->Draw();
  TLegend *leg2 = new TLegend(0.65,0.73,0.86,0.87);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);
  leg2->AddEntry("dataF","Data", "P");
  leg2->AddEntry("modelF","Signal + background","L");
//   leg2->AddEntry("backgroundF","Background only", "L");
  leg2->AddEntry("signal onlyF","Signal only", "L");
  leg2->AddEntry("ZttF","Ztt", "L");
  leg2->AddEntry("WenF","Wen", "L");
  leg2->AddEntry("TTbF","TTb", "L");
  leg2->AddEntry("QCDF","QCD", "L");
  leg2->Draw();

  string fileNameF = "Fail_"+tnp_+"_"+categoryDATA_;
  if(doMC_) fileNameF = fileNameF+"_MC";
  else fileNameF = fileNameF+"_DATA";
  if(MVAIso)cFail->SaveAs(Form("./ElecTauTnP/plots/computeElecTauTnPMacro_ClosureFit_MVAIso_%s_%s.pdf",fileNameF.c_str(), binLabel_.c_str()));
  else cFail->SaveAs(Form("./ElecTauTnP/plots/computeElecTauTnPMacro_ClosureFit_%s_%s.pdf",fileNameF.c_str(), binLabel_.c_str()));
  ResDataCombinedFit->printArgs(std::cout);


//   float DataErrorLo = DataEffFit->getErrorLo()<0 ? DataEffFit->getErrorLo() : (-1)*DataEffFit->getError();
//   float DataErrorHi = DataEffFit->getErrorHi()>0 ? DataEffFit->getErrorHi() : DataEffFit->getError();

//   cout <<"Data efficiency fit value : "<< DataEffFit->getVal() << " +/- " << DataEffFit->getError() << "  ( " << DataErrorLo << ", " << DataErrorHi << ")" <<  endl;
  cout <<"Data efficiency fit value : "<< CoeffZeeP.getVal()/(CoeffZeeP.getVal()+CoeffZeeF.getVal()) << " +/- " <<   endl;

  float DataEffFit =  CoeffZeeP.getVal()/(CoeffZeeP.getVal()+CoeffZeeF.getVal());

  out1[0][0]=(binCenter_);
  out1[0][1]=(binWidth_);
  out1[0][2]=(binWidth_);
  out1[0][3]=(McTruthEff);
  out1[0][4]=(BinomialError);
  out1[0][5]=(BinomialError);

  out1[1][0]=(binCenter_);
  out1[1][1]=(binWidth_);
  out1[1][2]=(binWidth_);
//   out1[1][3]=(DataEffFit->getVal());
  out1[1][3]=(DataEffFit);
//   out1[1][4]=((-1)*DataErrorLo);
//   out1[1][5]=(DataErrorHi);

  return;


 
}//end simFitClosure





void FitTest(const string tnp_          = "electauTnP",
	     const string categoryMC_   = "passingIsoLooseEleVetoMedium",
	     const string categoryDATA_ = "passingIsoLooseEleVetoMedium",
	     const string ptRange_      = "Pt>20",
	     const string binLabel_     = "Over20",
	     double cutValue_           = 0.5,
	     const float xLow_          = 60,
	     const float xHigh_         = 120,
	     const float nBins_         = 24,
	     bool doBinned_             = true,
	     float deltaAlpha_          = 0.0,
	     float deltaN_              = 0.0,
	     float mcLumi_              = 620,
	     bool doMC_                 = false,
	     bool PUWeight              = true,
	     bool SecondEleVeto         = true,
	     bool MVAIso                = false
	     )
{ 
 
  string outName = "./ElecTauTnP/electauTnP_FitTest_"+binLabel_+"_"+categoryDATA_;
  if(MVAIso) outName = outName+"_MVAIso";
  if(doMC_) outName = outName+"_MC";
  else  outName = outName+"_DATA";
  outName = outName+".txt";
  std::ofstream out(outName.c_str());
  out.precision(4);

  Float_t bin[2][6];
  if(doMC_)simFitClosure(bin,tnp_,categoryMC_,categoryDATA_,binLabel_,cutValue_,"(Eta<2.3 && Eta>-2.3 && "+ptRange_+")", 0.75,0.75,xLow_,xHigh_,nBins_,doBinned_,deltaAlpha_,deltaN_,mcLumi_,doMC_,PUWeight,SecondEleVeto,MVAIso);
  else simFit(bin,tnp_,categoryMC_,categoryDATA_,binLabel_,cutValue_,"(Eta<2.3 && Eta>-2.3 && "+ptRange_+")", 0.75,0.75,xLow_,xHigh_,nBins_,doBinned_,deltaAlpha_,deltaN_,mcLumi_,doMC_,PUWeight,SecondEleVeto,MVAIso);

  out<<"%Tag&probe measurement :" << categoryDATA_ << ", bin: " << ptRange_ << endl;
  out<<"\\begin{tabular}[!htbp]{|c|c|c|c|}" << endl;
  out << "\\hline" << endl;
  out << "Bin & MC & Data && SF\\\\" << endl;
  out << "\\hline" << endl;
  out << "$|\\eta|<2.3$ & " << bin[0][3] << " $\\pm$ " << bin[0][4] << " & " << bin[1][3] << " $\\pm_{" << bin[1][4] << "}^{"<< bin[1][5] << "}$ "<< " & " << bin[1][3]/bin[0][3] <<"\\\\" << endl;
  out << "\\hline" << endl;
  out<<"\\end{tabular}"<<endl;
}


void FitAllTest(){
  /////////////////////////////////////////////////
  /////////DATA
  //////////////////////////////////////////////////
  
  ///MVA Isolation
//   cout<<"Barrel :"<<endl;
//   FitTest("electauTnP","AntiEleLoose","passingIsoMVALooseEleVetoLoose","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,2676,false,true,false,true);
//   FitTest("electauTnP","AntiEleMedium","passingIsoMVALooseEleVetoMedium","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,2676,false,true,false,true);
//   FitTest("electauTnP","AntiEleTight","passingIsoMVALooseEleVetoTight","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,2676,false,true,false,true);
  FitTest("electauTnP","passingIsoMVALooseEleVetoMVA","passingIsoMVALooseEleVetoMVA","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,2676,false,true,false,true);

//   cout<<"Endcap :"<<endl;
//   FitTest("electauTnP","AntiEleLoose","passingIsoMVALooseEleVetoLoose","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,2676,false,true,false,true);
//   FitTest("electauTnP","AntiEleMedium","passingIsoMVALooseEleVetoMedium","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,2676,false,true,false,true);
//   FitTest("electauTnP","AntiEleTight","passingIsoMVALooseEleVetoTight","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,2676,false,true,false,true);
  FitTest("electauTnP","passingIsoMVALooseEleVetoMVA","passingIsoMVALooseEleVetoMVA","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,2676,false,true,false,true);

  ///HPS Isolation
//   cout<<"Barrel :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoLoose","passingIsoLooseEleVetoLoose","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,2676,false,true,false,false);
//   FitTest("electauTnP","passingIsoLooseEleVetoMedium","passingIsoLooseEleVetoMedium","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,2676,false,true,false,false);
//   FitTest("electauTnP","passingIsoLooseEleVetoTight","passingIsoLooseEleVetoTight","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,2676,false,true,false,false);
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA","passingIsoLooseEleVetoMVA","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,2676,false,true,false,false);

//   cout<<"Endcap :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoLoose","passingIsoLooseEleVetoLoose","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,2676,false,true,false,false);
//   FitTest("electauTnP","passingIsoLooseEleVetoMedium","passingIsoLooseEleVetoMedium","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,2676,false,true,false,false);
//   FitTest("electauTnP","passingIsoLooseEleVetoTight","passingIsoLooseEleVetoTight","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,2676,false,true,false,false);
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA","passingIsoLooseEleVetoMVA","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,2676,false,true,false,false);

//SecondEleVeto
  FitTest("electauTnP","passingIsoMVALooseEleVetoMVA","passingIsoMVALooseEleVetoMVA","Eta<1.5 && Eta>-1.5 && Pt>20","SecondEleVetoBarrelOver20",0.5,60,120,24,true,0.0,0.0,2676,false,true,true,true);
  FitTest("electauTnP","passingIsoMVALooseEleVetoMVA","passingIsoMVALooseEleVetoMVA","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","SecondEleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,2676,false,true,true,true);
  FitTest("electauTnP","passingIsoMVALooseEleVetoMVA","passingIsoMVALooseEleVetoMVA","Eta<2.3 && Eta>-2.3 && Pt>20","SecondEleVetoAllOver20",0.5,60,120,24,true,0.0,0.0,2676,false,true,true,true);

//  cout<<"Barrel :"<<endl;
//   cout<<"passingIsoLooseEleVetoLoose :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoLoose","passingIsoLooseEleVetoLoose","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMedium :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMedium","passingIsoLooseEleVetoMedium","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoTight :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoTight","passingIsoLooseEleVetoTight","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false); 
//   cout<<"passingIsoLooseEleVetoMVA :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA","passingIsoLooseEleVetoMVA","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoMVALooseEleVetoMVA :"<<endl;
//   FitTest("electauTnP","passingIsoMVALooseEleVetoMVA","passingIsoMVALooseEleVetoMVA","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose1","passingIsoLooseEleVetoMVA2Loose1","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose2","passingIsoLooseEleVetoMVA2Loose2","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium1","passingIsoLooseEleVetoMVA2Medium1","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium2","passingIsoLooseEleVetoMVA2Medium2","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight1","passingIsoLooseEleVetoMVA2Tight1","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight2","passingIsoLooseEleVetoMVA2Tight2","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight1","passingIsoLooseEleVetoMVA2VTight1","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight2","passingIsoLooseEleVetoMVA2VTight2","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);

//  cout<<"Endcap :"<<endl;
//   cout<<"passingIsoLooseEleVetoLoose :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoLoose","passingIsoLooseEleVetoLoose","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMedium :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMedium","passingIsoLooseEleVetoMedium","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoTight :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoTight","passingIsoLooseEleVetoTight","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);  
//   cout<<"passingIsoLooseEleVetoMVA :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA","passingIsoLooseEleVetoMVA","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoMVALooseEleVetoMVA :"<<endl;
//   FitTest("electauTnP","passingIsoMVALooseEleVetoMVA","passingIsoMVALooseEleVetoMVA","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose1","passingIsoLooseEleVetoMVA2Loose1","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose2","passingIsoLooseEleVetoMVA2Loose2","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium1","passingIsoLooseEleVetoMVA2Medium1","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);  
//   cout<<"passingIsoLooseEleVetoMVA2 Medium2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium2","passingIsoLooseEleVetoMVA2Medium2","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight1","passingIsoLooseEleVetoMVA2Tight1","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight2","passingIsoLooseEleVetoMVA2Tight2","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight1","passingIsoLooseEleVetoMVA2VTight1","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight2","passingIsoLooseEleVetoMVA2VTight2","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);


//   cout<<"Barrel and Endcap :"<<endl;
//   cout<<"passingIsoLooseEleVetoLoose :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoLoose","passingIsoLooseEleVetoLoose"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMedium :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMedium","passingIsoLooseEleVetoMedium"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoTight :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoTight","passingIsoLooseEleVetoTight"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);  
//   cout<<"passingIsoLooseEleVetoMVA :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA","passingIsoLooseEleVetoMVA"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose1","passingIsoLooseEleVetoMVA2Loose1"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose2","passingIsoLooseEleVetoMVA2Loose2"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium1","passingIsoLooseEleVetoMVA2Medium1"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);  
//   cout<<"passingIsoLooseEleVetoMVA2 Medium2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium2","passingIsoLooseEleVetoMVA2Medium2"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight1","passingIsoLooseEleVetoMVA2Tight1"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight2","passingIsoLooseEleVetoMVA2Tight2"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight1","passingIsoLooseEleVetoMVA2VTight1"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight2","passingIsoLooseEleVetoMVA2VTight2"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,true,false,false);


//////////////////////////////////////////
///////MC closure test
//////////////////////////////////////////

//   FitTest("electauTnP","passingIsoLooseEleVetoLoose","passingIsoLooseEleVetoLoose","Eta<1.5 && Eta>-1.5 && Pt>20 ","BarrelOver20",0.5,60,120,24,true,0.0,0.0,2676,true,true,false,false);
//   FitTest("electauTnP","passingIsoLooseEleVetoMedium","passingIsoLooseEleVetoMedium","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,2676,true,true,false,false);
//   FitTest("electauTnP","passingIsoLooseEleVetoTight","passingIsoLooseEleVetoTight","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,2676,true,true,false,false); 
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA","passingIsoLooseEleVetoMVA","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,2676,true,true,false,false);

//   FitTest("electauTnP","passingIsoLooseEleVetoLoose","passingIsoLooseEleVetoLoose","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,2676,true,true,false,false);
//   FitTest("electauTnP","passingIsoLooseEleVetoMedium","passingIsoLooseEleVetoMedium","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,2676,true,true,false,false);
//   FitTest("electauTnP","passingIsoLooseEleVetoTight","passingIsoLooseEleVetoTight","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,2676,true,true,false,false);   
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA","passingIsoLooseEleVetoMVA","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,2676,true,true,false,false);


  /////MVA Isolation
  cout<<"Barrel :"<<endl;
//   FitTest("electauTnP","AntiEleLoose","passingIsoMVALooseEleVetoLoose","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,2676,true,true,false,true);
//   FitTest("electauTnP","AntiEleMedium","passingIsoMVALooseEleVetoMedium","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,2676,true,true,false,true);
//   FitTest("electauTnP","AntiEleTight","passingIsoMVALooseEleVetoTight","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,2676,true,true,false,true);
//   FitTest("electauTnP","passingIsoMVALooseEleVetoMVA","passingIsoMVALooseEleVetoMVA","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,2676,true,true,false,true);

  cout<<"Endcap :"<<endl;
//   FitTest("electauTnP","AntiEleLoose","passingIsoMVALooseEleVetoLoose","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,2676,true,true,false,true);
//   FitTest("electauTnP","AntiEleMedium","passingIsoMVALooseEleVetoMedium","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,2676,true,true,false,true);
//   FitTest("electauTnP","AntiEleTight","passingIsoMVALooseEleVetoTight","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,2676,true,true,false,true);
//   FitTest("electauTnP","passingIsoMVALooseEleVetoMVA","passingIsoMVALooseEleVetoMVA","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,2676,true,true,false,true);



//SECONDELEVETO
//   cout<<endl;
//   cout<<"AFTER SECOND ELE VETO : "<<endl;

//   cout<<"Barrel :"<<endl;
//   cout<<"passingIsoLooseEleVetoLoose :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoLoose","passingIsoLooseEleVetoLoose","Eta<1.5 && Eta>-1.5 && Pt>20","EleVetoBarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMedium :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMedium","passingIsoLooseEleVetoMedium","Eta<1.5 && Eta>-1.5 && Pt>20","EleVetoBarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoTight :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoTight","passingIsoLooseEleVetoTight","Eta<1.5 && Eta>-1.5 && Pt>20","EleVetoBarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);  
//   cout<<"passingIsoLooseEleVetoMVA :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA","passingIsoLooseEleVetoMVA","Eta<1.5 && Eta>-1.5 && Pt>20","EleVetoBarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose1","passingIsoLooseEleVetoMVA2Loose1","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose2","passingIsoLooseEleVetoMVA2Loose2","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium1","passingIsoLooseEleVetoMVA2Medium1","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium2","passingIsoLooseEleVetoMVA2Medium2","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight1","passingIsoLooseEleVetoMVA2Tight1","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight2","passingIsoLooseEleVetoMVA2Tight2","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight1","passingIsoLooseEleVetoMVA2VTight1","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight2","passingIsoLooseEleVetoMVA2VTight2","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);

//   cout<<"Endcap : "<<endl;
//   cout<<"passingIsoLooseEleVetoLoose :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoLoose","passingIsoLooseEleVetoLoose","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMedium :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMedium","passingIsoLooseEleVetoMedium","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoTight :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoTight","passingIsoLooseEleVetoTight","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);  
//   cout<<"passingIsoLooseEleVetoMVA :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA","passingIsoLooseEleVetoMVA","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose1","passingIsoLooseEleVetoMVA2Loose1","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose2","passingIsoLooseEleVetoMVA2Loose2","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium1","passingIsoLooseEleVetoMVA2Medium1","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium2","passingIsoLooseEleVetoMVA2Medium2","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight1","passingIsoLooseEleVetoMVA2Tight1","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight2","passingIsoLooseEleVetoMVA2Tight2","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight1","passingIsoLooseEleVetoMVA2VTight1","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight2","passingIsoLooseEleVetoMVA2VTight2","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);

//   cout<<"Barrel And Endcap : "<<endl;
//   cout<<"passingIsoLooseEleVetoLoose :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoLoose","passingIsoLooseEleVetoLoose"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMedium :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMedium","passingIsoLooseEleVetoMedium"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoTight :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoTight","passingIsoLooseEleVetoTight"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);  
//   cout<<"passingIsoLooseEleVetoMVA :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA","passingIsoLooseEleVetoMVA"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose1","passingIsoLooseEleVetoMVA2Loose1"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose2","passingIsoLooseEleVetoMVA2Loose2"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium1","passingIsoLooseEleVetoMVA2Medium1"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium2","passingIsoLooseEleVetoMVA2Medium2"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight1","passingIsoLooseEleVetoMVA2Tight1"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight2","passingIsoLooseEleVetoMVA2Tight2"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight1","passingIsoLooseEleVetoMVA2VTight1"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight2","passingIsoLooseEleVetoMVA2VTight2"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,true,true,false);




  /////////  /////////  /////////  /////////  /////////  /////////  /////////
  ///No PUWeight
  /////////  /////////  /////////  /////////  /////////  /////////

//  cout<<"Barrel :"<<endl;
//   cout<<"passingIsoLooseEleVetoLoose :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoLoose","passingIsoLooseEleVetoLoose","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMedium :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMedium","passingIsoLooseEleVetoMedium","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoTight :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoTight","passingIsoLooseEleVetoTight","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false); 
//   cout<<"passingIsoLooseEleVetoMVA :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA","passingIsoLooseEleVetoMVA","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoMVALooseEleVetoMVA :"<<endl;
//   FitTest("electauTnP","passingIsoMVALooseEleVetoMVA","passingIsoMVALooseEleVetoMVA","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose1","passingIsoLooseEleVetoMVA2Loose1","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose2","passingIsoLooseEleVetoMVA2Loose2","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium1","passingIsoLooseEleVetoMVA2Medium1","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium2","passingIsoLooseEleVetoMVA2Medium2","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight1","passingIsoLooseEleVetoMVA2Tight1","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight2","passingIsoLooseEleVetoMVA2Tight2","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight1","passingIsoLooseEleVetoMVA2VTight1","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight2","passingIsoLooseEleVetoMVA2VTight2","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);

//  cout<<"Endcap :"<<endl;
//   cout<<"passingIsoLooseEleVetoLoose :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoLoose","passingIsoLooseEleVetoLoose","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMedium :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMedium","passingIsoLooseEleVetoMedium","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoTight :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoTight","passingIsoLooseEleVetoTight","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);  
//   cout<<"passingIsoLooseEleVetoMVA :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA","passingIsoLooseEleVetoMVA","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoMVALooseEleVetoMVA :"<<endl;
//   FitTest("electauTnP","passingIsoMVALooseEleVetoMVA","passingIsoMVALooseEleVetoMVA","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose1","passingIsoLooseEleVetoMVA2Loose1","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose2","passingIsoLooseEleVetoMVA2Loose2","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium1","passingIsoLooseEleVetoMVA2Medium1","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);  
//   cout<<"passingIsoLooseEleVetoMVA2 Medium2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium2","passingIsoLooseEleVetoMVA2Medium2","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight1","passingIsoLooseEleVetoMVA2Tight1","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight2","passingIsoLooseEleVetoMVA2Tight2","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight1","passingIsoLooseEleVetoMVA2VTight1","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight2","passingIsoLooseEleVetoMVA2VTight2","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);


//   cout<<"Barrel and Endcap :"<<endl;
//   cout<<"passingIsoLooseEleVetoLoose :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoLoose","passingIsoLooseEleVetoLoose"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMedium :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMedium","passingIsoLooseEleVetoMedium"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoTight :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoTight","passingIsoLooseEleVetoTight"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);  
//   cout<<"passingIsoLooseEleVetoMVA :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA","passingIsoLooseEleVetoMVA"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose1","passingIsoLooseEleVetoMVA2Loose1"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose2","passingIsoLooseEleVetoMVA2Loose2"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium1","passingIsoLooseEleVetoMVA2Medium1"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);  
//   cout<<"passingIsoLooseEleVetoMVA2 Medium2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium2","passingIsoLooseEleVetoMVA2Medium2"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight1","passingIsoLooseEleVetoMVA2Tight1"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight2","passingIsoLooseEleVetoMVA2Tight2"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight1","passingIsoLooseEleVetoMVA2VTight1"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight2","passingIsoLooseEleVetoMVA2VTight2"," Pt>20","Over20",0.5,60,120,24,true,0.0,0.0,620,false,false,false,false);    
// //SECONDELEVETO
//   cout<<endl;
//   cout<<"AFTER SECOND ELE VETO : "<<endl;

//   cout<<"Barrel :"<<endl;
//   cout<<"passingIsoLooseEleVetoLoose :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoLoose","passingIsoLooseEleVetoLoose","Eta<1.5 && Eta>-1.5 && Pt>20","EleVetoBarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMedium :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMedium","passingIsoLooseEleVetoMedium","Eta<1.5 && Eta>-1.5 && Pt>20","EleVetoBarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoTight :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoTight","passingIsoLooseEleVetoTight","Eta<1.5 && Eta>-1.5 && Pt>20","EleVetoBarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);  
//   cout<<"passingIsoLooseEleVetoMVA :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA","passingIsoLooseEleVetoMVA","Eta<1.5 && Eta>-1.5 && Pt>20","EleVetoBarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose1","passingIsoLooseEleVetoMVA2Loose1","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose2","passingIsoLooseEleVetoMVA2Loose2","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium1","passingIsoLooseEleVetoMVA2Medium1","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium2","passingIsoLooseEleVetoMVA2Medium2","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight1","passingIsoLooseEleVetoMVA2Tight1","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight2","passingIsoLooseEleVetoMVA2Tight2","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight1","passingIsoLooseEleVetoMVA2VTight1","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight2","passingIsoLooseEleVetoMVA2VTight2","Eta<1.5 && Eta>-1.5 && Pt>20","BarrelOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);

//   cout<<"Endcap : "<<endl;
//   cout<<"passingIsoLooseEleVetoLoose :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoLoose","passingIsoLooseEleVetoLoose","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMedium :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMedium","passingIsoLooseEleVetoMedium","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoTight :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoTight","passingIsoLooseEleVetoTight","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);  
//   cout<<"passingIsoLooseEleVetoMVA :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA","passingIsoLooseEleVetoMVA","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose1","passingIsoLooseEleVetoMVA2Loose1","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose2","passingIsoLooseEleVetoMVA2Loose2","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium1","passingIsoLooseEleVetoMVA2Medium1","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium2","passingIsoLooseEleVetoMVA2Medium2","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight1","passingIsoLooseEleVetoMVA2Tight1","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight2","passingIsoLooseEleVetoMVA2Tight2","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight1","passingIsoLooseEleVetoMVA2VTight1","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight2","passingIsoLooseEleVetoMVA2VTight2","((Eta>1.5 && Eta<2.3) || (Eta<-1.5 && Eta>-2.3)) && Pt>20","EleVetoEndcapOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);

//   cout<<"Barrel And Endcap : "<<endl;
//   cout<<"passingIsoLooseEleVetoLoose :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoLoose","passingIsoLooseEleVetoLoose"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMedium :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMedium","passingIsoLooseEleVetoMedium"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoTight :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoTight","passingIsoLooseEleVetoTight"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);  
//   cout<<"passingIsoLooseEleVetoMVA :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA","passingIsoLooseEleVetoMVA"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose1","passingIsoLooseEleVetoMVA2Loose1"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Loose2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Loose2","passingIsoLooseEleVetoMVA2Loose2"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium1","passingIsoLooseEleVetoMVA2Medium1"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Medium2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Medium2","passingIsoLooseEleVetoMVA2Medium2"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight1","passingIsoLooseEleVetoMVA2Tight1"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 Tight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2Tight2","passingIsoLooseEleVetoMVA2Tight2"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight1 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight1","passingIsoLooseEleVetoMVA2VTight1"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false);
//   cout<<"passingIsoLooseEleVetoMVA2 VTight2 :"<<endl;
//   FitTest("electauTnP","passingIsoLooseEleVetoMVA2VTight2","passingIsoLooseEleVetoMVA2VTight2"," Pt>20","EleVetoOver20",0.5,60,120,24,true,0.0,0.0,620,false,false,true,false); 

}


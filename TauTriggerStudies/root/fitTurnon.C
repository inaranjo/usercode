

//////////////////////////////////////////////////////////////////////////
//
// FITTING WITH A ROOFIT-USER DEFINED CRYSTAL BALL
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TFrame.h"
#include "TStyle.h"
#include "TFile.h"
#include "iostream"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>


using namespace RooFit ;

void loadPresentationStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(0111);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYOffset(0.01);
  gStyle->SetLabelOffset(0.005, "XYZ");
  //  gStyle->SetLabelOffset(0.1, "XYZ");
  gStyle->SetTitleSize(0.07, "XYZ");
  gStyle->SetTitleFont(22,"X");
  gStyle->SetTitleFont(22,"Y");
  gStyle->SetPadBottomMargin(0.13);



  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetHistLineWidth(2);

  gROOT->ProcessLine(".L tdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");


}

void paramOn(RooPlot* frame, TString label, RooRealVar vars[4],  
	     double xmin, double ymin, double xmax, double ymax, int color)
 {
   TPaveText * text = new TPaveText (xmin,ymin,xmax,ymax,"ndc");
   text->SetTextAlign(22);
   text->SetTextColor(color);
   text->SetLineColor(color);
   text->SetTextFont(42);
   text->SetFillColor(kWhite);
   // text->SetTextSize(1.);
   text->SetShadowColor(0);
   if( label != "" ) { text->AddText( label ); }

   frame->addObject(text);
 }


TLegend* CreateLegend(  TString label, std::vector<RooRealVar> vars,
		       double xmin, double ymin, double xmax, double ymax,
		       int color, int style)
{
  TH1F* hist = new TH1F("hist","hist",10,0,1);
  hist->SetMarkerStyle(style);
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  TString* text = new TString("test\n erwer");

  TLegend * leg = new TLegend (xmin,ymin,xmax,ymax,"","ndc");

  leg->SetHeader(text->Data());
  leg->SetLineColor(color);
  leg->SetFillColor(kWhite);
  leg->SetTextColor(color);
  leg->AddEntry(hist, label, "lp");

  return leg;
}


void turnon_sigma(TH1D * h, int nSig, double val, double & errmin, double & errmax, double & frac)
{
  double nEntries = h->Integral(1,h->GetNbinsX());
  double fracLR = (1-TMath::Erf(nSig/sqrt(2)))/2. ;
  double nEntriesLR = fracLR*nEntries ;

  // left interval
  int binL, binR ;
  double nEntriesL, nEntriesR ;
  for (int bin=1 ; bin<=h->GetNbinsX() ; bin++) {
    nEntriesL = h->Integral(1,bin);
    if (nEntriesL >= nEntriesLR) {
      binL = bin ;
      break ;
    }
  }
  // right interval
  for (int bin=h->GetNbinsX() ; bin>=1 ; bin--) {
    nEntriesR = h->Integral(bin,h->GetNbinsX());
    if (nEntriesR >= nEntriesLR) {
      binR = bin ;
      break ;
    }
  }
  //check fraction
  frac = h->Integral(binL, binR)/nEntries ;

  //results
  //errmin = val-h->GetBinLowEdge(binL) ;
  //errmax = h->GetBinLowEdge(binR)-val ;
  errmin = h->GetBinLowEdge(binL) ;
  errmax = h->GetBinLowEdge(binR) ;
}



void getTurnon(RooFitResult * fit, bool EB, bool MC=false, int nbdata=1)
{
#ifdef __CINT__
  gROOT->ProcessLineSync(".x FuncCB.cxx+") ;
#endif
  cout<<"Printing fit results:"<<endl ;
  fit->Print() ;
  cout<<"results printed"<<endl ;

  RooArgList param = fit->floatParsFinal() ;
  double err[5] ;
  double mu[5] ;
  for(Int_t i = 0; i < param.getSize(); i++) 
    {
      RooRealVar* var = ( dynamic_cast<RooRealVar*>( param.at(i) ) );
      var->Print() ;
      mu[i] = var->getVal() ;
      err[i] = var->getError() ;
      //cout<<var->getVal()<<" "<<var->getErrorLo()<<" "<<var->getErrorHi()<<endl ;
    }
  
  // PARAMETRES ROOFIT CRYSTAL BALL
  double min, max ;
  
  min = mu[0]-5*err[0] ;
  if (mu[0]-5*err[0]<0) min = 0. ;
  RooRealVar alpha("alpha","#alpha",mu[0],min,mu[0]+5*err[0]);
  
  min = mu[1]-5*err[1] ;
  if (mu[1]-5*err[1]<5) min = 5. ;
  RooRealVar mean("mean","mean",mu[1],min,mu[1]+5*err[1]);
  
  min = mu[2]-5*err[2] ;
  if (mu[2]-5*err[2]<1) min = 1. ;
  RooRealVar n("n","n",mu[2],min,mu[2]+5*err[2]);
  
  min = mu[3]-5*err[3] ;
  if (mu[3]-5*err[3]<0.6) min = 0.6 ; 
  max = mu[3]+5*err[3] ;
  if (mu[3]+5*err[3]>1.) max = 1. ; 
  RooRealVar norm("norm","N",mu[3],min,max);
  
  min = mu[4]-5*err[4] ;
  if (mu[4]-5*err[4]<0.) min = 0. ; 
  RooRealVar sigma("sigma","#sigma",mu[4],min,mu[4]+5*err[4]);
  
  RooRealVar xaxis("sc_et","sc_et",0,100) ;


  bool fitok(true) ;
  for (int i=0 ; i<5 ; i++) if (mu[i]!=0 && fabs(err[i]/mu[i])>1000) fitok = false ; //error is 1000 times larger than fitted value

  RooAbsPdf* parabPdf  ;
  RooDataSet* data ;
  if (fitok) {
    cout<<"fit->createHessePdf"<<endl ;
    parabPdf = fit->createHessePdf(RooArgSet(norm,alpha,n,mean,sigma)) ;
    cout<<"parabPdf->generate"<<endl ;
    data = parabPdf->generate(RooArgSet(norm,alpha,n,mean,sigma),nbdata) ;
    TH1* tmp1 = data->createHistogram("sigma",100) ;
    //TH1* tmp1 = data->createHistogram("sigma,alpha",100,100) ;
    //tmp1->Draw() ;
  } else {
    cout<<"Fit not ok, will keep not generate values" ;
    cout<<"Correlation matrix is:"<<endl ;
    fit->correlationMatrix()->Print() ;
  }
  
  TFile * f ;
  string name("turnon") ;
  if (EB) name += "EB" ;
  else name += "EE" ;
  if (MC) name += "_MC" ;
  else name += "_data" ;
  name += ".root" ;
  if (EB) f = new TFile(name.c_str(), "RECREATE") ;
  else f = new TFile(name.c_str(), "RECREATE") ;
  
  TTree* treeFit = new TTree("treeFit", "treeFit");
  float tt_pt, tt_eff ;
  int tt_index ;
  treeFit->Branch("pt",&tt_pt, "pt/F");
  treeFit->Branch("eff" ,&tt_eff,  "eff/F");
  treeFit->Branch("index" ,&tt_index,  "index/I");

  double effplateau = 0 ;
  double et02bref = 0  ;
  double et05bref = 0  ;
  double et08bref = 0  ;
  double et095bref = 0  ;
  double etplateauref = 0 ; 

  int xbin = 3000 ;
  double xmax = 100 ;
  double binw = xmax/xbin ;
  TH1D * hcb = new TH1D("hcb","hcb",xbin,0,xmax) ;
  TH2D * h2cb = new TH2D("h2cb","h2cb",xbin,0,xmax, 1100, 0, 1.1) ;
  
  for (int i=-1 ; i<nbdata ; i++) {
    
    //if (i%100==0) cout<<"=== "<<i<<" loops"<<endl ;
    
    if (i==-1 || !fitok) {
      alpha = mu[0] ;
      mean = mu[1] ;
      n = mu[2] ;
      norm = mu[3] ;
      sigma = mu[4] ;
    }
    else {
      RooArgSet* set = data->get(i) ;
      alpha = set->getRealValue("alpha") ;
      mean = set->getRealValue("mean") ;
      n = set->getRealValue("n") ;
      norm = set->getRealValue("norm") ;
      sigma = set->getRealValue("sigma") ;
    }
    if (n.getVal()<=0) continue ; // avoid meaningless parameters
    
    FuncCB cb("cb","Crystal Ball Integree",xaxis,mean,sigma,alpha,n,norm) ;
    double epsilon = 0.00001 ;

    for (int imyet=0 ; imyet<xbin ; imyet++) {
      double myet = imyet*binw + epsilon ;
      xaxis = myet ;
      tt_index = i ;
      tt_pt = myet ;
      tt_eff = cb.getVal() ;
      treeFit->Fill() ;
    }

    if (i==-1) {
      for (int imyet=0 ; imyet<xbin ; imyet++) {
	double myet = imyet*binw + epsilon ;
	xaxis = myet ;
	hcb->Fill(myet, cb.getVal()) ;
      }
    }
    for (int imyet=0 ; imyet<xbin ; imyet++) {
      double myet = imyet*binw + epsilon ;
      xaxis = myet ; 
      h2cb->Fill(myet, cb.getVal()) ;
    }
    
    
    // find turn-on points
    double diff02b = 99. ; 
    double et02b = 0. ;
    double diff05b = 99. ; 
    double et05b = 0. ;
    double diff08b = 99. ; 
    double et08b = 0. ;
    double diff095b = 99. ; 
    double et095b = 0.  ;
    double etplateau = 0. ;
    
    double cbval ;
    
    if (i==-1) {
      xaxis = 100 ; // plateau
      effplateau = cb.getVal() ;
      for (double myet=8 ; myet<100 ; myet+=0.01) {
	xaxis = myet ;
	cbval = cb.getVal() ;
	if (0.99*effplateau-cbval<0) {
	  etplateau = myet ;
	  break ;
	}
      }
    
      for (double myet=8 ; myet<100 ; myet+=0.01) {
	xaxis = myet ;
	cbval = cb.getVal() ;
	if (cbval>=0) {
	  if ( fabs(cbval-0.2) < diff02b) {
	    diff02b = fabs(cbval-0.2) ;
	    et02b = myet ;
	  }
	  if ( fabs(cbval-0.5) < diff05b) {
	    diff05b = fabs(cbval-0.5) ;
	    et05b = myet ;
	  }
	  if ( fabs(cbval-0.8) < diff08b) {
	    diff08b = fabs(cbval-0.8) ;
	    et08b = myet ;
	  }
	  if ( fabs(cbval-0.95) < diff095b) {
	    diff095b = fabs(cbval-0.95) ;
	    et095b = myet ;
	  }
	}
      }
    
      et02bref = et02b ;
      et05bref = et05b ;
      et08bref = et08b ;
      et095bref = et095b ;
      etplateauref = etplateau ;
    }
  }
  
  // measure the 1 sigma interval
  double errmin, errmax, emost, frac ;
  cout<<"turn-on at 0.5 = "<<et05bref<<endl ;
  cout<<"Delta turn-on 0.2-0.8 = "<<et08bref-et02bref<<endl ;
  cout<<"turn-on at 0.95 = "<<et095bref<<endl ;
  cout<<"efficiency plateau (100 GeV) = "<<effplateau<<endl ;
  cout<<"0.99*plateau reached at "<<etplateauref<<endl ;
  
  //Storing
  treeFit->Write() ;
  hcb->Write() ;
  h2cb->Write() ;
  f->Close() ;
}


void getIntegralEfficiency(TTree * treeTnP,double pTmin, bool barrel, bool endcap, double massMin, double massMax)
{
  float t_mass, t_etraw,t_et, t_pt, t_eta, t_weight ;
  int t_match, t_L1match ;
  treeTnP->SetBranchAddress("mass", &t_mass);
  treeTnP->SetBranchAddress("pt", &t_pt);
  treeTnP->SetBranchAddress("eta", &t_eta);
  treeTnP->SetBranchAddress("match", &t_match);
  treeTnP->SetBranchAddress("weight", &t_weight);
  float pTmax = 1000 ; 

  TH1F * htot  = new TH1F("htot","htot",1,pTmin,pTmax) ;
  TH1F * hpass  = new TH1F("hpass","hpass",1,pTmin,pTmax) ;

  for (int nEntries=0 ; nEntries<treeTnP->GetEntries() ; nEntries++) {
    treeTnP->GetEntry (nEntries) ;

    if (t_mass<=massMin) continue ;
    if (t_mass>=massMax) continue ;
    if (barrel && endcap &&  fabs(t_eta)>3) continue ;
    if (barrel && !endcap &&  fabs(t_eta)>1.47) continue ;
    if (!barrel && endcap &&  fabs(t_eta)<1.47) continue ;
    if (t_pt<pTmin) continue ;

    htot->Fill(t_pt) ;
    if (t_match==1) hpass->Fill(t_pt) ;
  }
  TGraphAsymmErrors* averageEff = new TGraphAsymmErrors(hpass,htot) ; 
  //averageEff->Draw("ALP") ;

  cout<<"Integrated efficiency over pT>"<<pTmin<<" is "<<(averageEff->GetY())[0]<<
    " +/- "<<averageEff->GetErrorYhigh(0)<<" / "<<averageEff->GetErrorYlow(0)<<endl ;
}



RooFitResult * fitTurnon(bool MC=false, 
			 			int nbdata=1000, 
			 			double pTmin=20, 
			 			string dataset="2011AMay10ReReco",//"2011AAug05ReReco",//"2011AMay10ReReco",//"2011APromptRecov6"
			 			int PtThresh = 10
			 			)
{
	
	

	bool barrel = true;
 	bool endcap = true;
  
  int answer ;
  cout<<"Loose(0) or Medium(1) or Tight (2)?"<<endl ;
  cin >> answer ;
  bool loose(false) ;
  bool medium(false) ;
  bool tight(false) ;
  if (answer==0) loose = true ; 
  if (answer==1) medium = true ; 
  if (answer==2) tight = true ; 

  int pause ;
  gROOT->Reset();
  loadPresentationStyle();  
  gROOT->ForceStyle();
  gSystem->Load("libRooFit.so");

  //binning
  if(loose && PtThresh==20){
  const int nbins = 17;
  //Double_t bins[nbins] = {2.0, 4.0, 6.0, 8.0, 10., 12., 14., 16., 18., 20., 22.,24.,26.,28.,30., 40., 50.}; 
  Double_t bins[nbins] = {2.0, 4.0, 6.0, 8.0, 10., 12., 14., 16., 18.,20.,21., 23., 24.,28.,32., 40., 50.}; 
  }
  if(loose && PtThresh==20 && MC){
  const int nbins = 15;
  //Double_t bins[nbins] = {2.0, 4.0, 6.0, 8.0, 10., 12., 14., 16., 18., 20., 22.,24.,26.,28.,30., 40., 50.}; 
  Double_t bins[nbins] = {2.0, 4.0, 6.0, 8.0, 10., 12., 14., 16., 18.,20., 23.,28.,32., 40., 50.}; 
  }
  if(loose && PtThresh == 15){
  const int nbins = 17;
  Double_t bins[nbins] = {2.0, 4.0, 6.0, 8.0, 10., 12.,13., 14., 15., 16.,18.,21., 24.,28.,32., 40., 50.}; //Tau15
  }  
  if(loose && PtThresh == 15 && MC){
  const int nbins = 13;
  Double_t bins[nbins] = {2.0, 4.0, 6.0, 8.0, 10., 12., 15., 16.,18.,24.,30., 40., 50.}; //Tau15
  }
  if(loose && PtThresh == 10){
  const int nbins = 17;
  Double_t bins[nbins] = {2.,4.,6.,7.,8.,9.,10.,13., 14., 16.,18.,21., 24.,28.,32., 40., 50.}; //Tau10
  }
  if(loose && PtThresh == 10 && MC){
  const int nbins = 14;
  Double_t bins[nbins] = {2.,4.,6.,7.,8.,9.,10.,13., 15.,18., 24.,30., 40., 50.}; //Tau10
  }
  if(medium){
  const int nbins = 15;
  //Double_t bins[nbins] = {2.0, 4.0, 6.0, 8.0, 10., 12., 14., 16., 18., 20., 22.,24.,26.,28.,30., 40., 50.}; 
  Double_t bins[nbins] = {2.0, 4.0, 6.0, 8.0, 10., 12., 14., 16., 18., 20., 22.,25.,30., 40., 50.}; 
  }
  if(tight){
  const int nbins = 17;
  //Double_t bins[nbins] = {2.0, 4.0, 6.0, 8.0, 10., 12., 14., 16., 18., 20., 22.,24.,26.,28.,30., 40., 50.}; 
  Double_t bins[nbins] = {2.0, 4.0, 6.0, 8.0, 10., 12., 14., 16., 18.,20.,21., 22., 24.,28.,30.,40, 50.}; 
  }

  
  double fit_cuts_min = 15 ;
  if(PtThresh == 10) fit_cuts_min = 10.5 ;
  if(PtThresh == 15) fit_cuts_min = 12 ;
  if(PtThresh == 20) fit_cuts_min = 12 ;

  double fit_cuts_max = 100 ;
  
  TString name_image="efficiencyTnP" ;  
  if (barrel) name_image += "EB" ;
  else name_image += "EE" ;
  if (MC) name_image += "_MC" ;
  else name_image += "_data" ;

  int colors[5] = {kRed,kBlack,kBlue,kRed,kBlack};
  int styles[5] = {kFullCircle, kFullTriangleUp, kFullTriangleDown, kFullSquare, kOpenStar};
  
  
  
  TFile *_file0 ;
  if (MC) _file0 = TFile::Open("output_MC.root");
  else  _file0 = TFile::Open("output_data.root");
  _file0->cd();

  cout << "Fit characteristics :" << endl ;
  cout << "Fit Range [" << fit_cuts_min << "," << fit_cuts_max << "]" << endl ;
  if (barrel) cout << "in barrel" << endl ;
  if (endcap) cout << "in endcap" << endl ;
  cout << "----------------------" << endl ;
  
  
  TTree* treeTnP = (TTree*) gDirectory->Get("treeTnP") ;
  
#ifdef __CINT__
  gROOT->ProcessLineSync(".x FuncCB.cxx+") ;
#endif
  
  
  // Variables from tree
  RooRealVar xaxis("pt","P_{T} [GeV]",0,100) ;
  RooCategory cut("match","trigger efficieny") ; // discrete variables
  cut.defineType("accept",1) ;
  cut.defineType("reject",0) ;
  RooCategory cutL1("L1match","L1 trigger efficieny") ; // discrete variables
  cutL1.defineType("accept",1) ;
  cutL1.defineType("reject",0) ;
  cutL1.defineType("notset",-1) ;
  RooRealVar mass("mass","mass",45,70) ;// consider only this mass range when importing data
  //RooRealVar mass("mass","mass",0,280) ;// consider only this mass range when importing data
  RooRealVar eta("eta","eta",-3., 3) ;
  RooRealVar weight("weight","weight",-1,1000) ;
  RooDataSet * dataSet ;
  if (barrel && endcap) 
    dataSet = new RooDataSet("data","data from tree",RooArgSet(xaxis, cut, cutL1, mass, eta, weight),WeightVar(weight), Import(*treeTnP), Cut("abs(eta)<3") );
  if (barrel && !endcap) 
    dataSet = new RooDataSet("data","data from tree",RooArgSet(xaxis, cut, cutL1, mass, eta, weight),WeightVar(weight),Import(*treeTnP), Cut("abs(eta)<1.47") );
  if (!barrel && endcap) 
    //dataSet = new RooDataSet("data","data from tree",RooArgSet(xaxis, cut, cutL1, mass, eta, weight),WeightVar(weight),Import(*treeTnP), Cut("abs(eta)>1.47") );
    dataSet = new RooDataSet("data","data from tree",RooArgSet(xaxis, cut, cutL1, mass, eta, weight),Import(*treeTnP), Cut("abs(eta)>1.47") );
  
//   for(int ii=0 ; ii<dataSet->numEntries() ; ii++) {
//     RooArgSet* set = dataSet->get(ii) ;
//     xaxis = set->getRealValue("pt") ;
//     mass = set->getRealValue("mass") ;
//     eta = set->getRealValue("eta") ;
//     weight = set->getRealValue("weight") ;
//     if (xaxis.getVal()>pTmin) 
//       cout<<ii<<" "<<xaxis<<" "<<mass<<" "<<eta<<" "<<weight<<endl ;
//   }

  // PARAMETRES ROOFIT CRYSTAL BALL
  RooRealVar norm("norm","N",0.9,0.6,1);
  RooRealVar alpha("alpha","#alpha",4,0.00001,8);
  RooRealVar n("n","n",1.1,1.1,35);
  RooRealVar mean("mean","mean",PtThresh,5,30);
  RooRealVar sigma("sigma","#sigma",2.5,0.01,5);
  mean.setVal(PtThresh);
//   RooRealVar norm("norm","N",0.8 ,0.6,1);
//   RooRealVar alpha("alpha","#alpha", 0.01  ,0.00001,8);
//   RooRealVar n("n","n",3,1.1,35);
//   RooRealVar mean("mean","mean",15  ,5,30);
//   RooRealVar sigma("sigma","#sigma",0.1,0.01,5);
  
  //FONCTION CRYSTALL BALL
  FuncCB cb("cb","Crystal Ball Integree",xaxis,mean,sigma,alpha,n,norm) ;
  
  
  
  //////////////////////////////////////:: Construct efficiency p.d.f eff(cut|x)/////////////////////////
  RooEfficiency eff("eff","efficiency",cb,cut,"accept");
  
  
  
  // DATA FILLED, d a t a   e f f i c i e n c y  
  // --------------------------------------------------------
  RooBinning binning;
  binning = RooBinning(nbins - 1, bins, "binning");
  
  
  RooPlot* frame = xaxis.frame(Bins(18000),Title("Fitted efficiency")) ;
  int color=colors[0];
  int style=styles[0];
  
  dataSet->plotOn(frame,Binning(binning),Efficiency(cut),MarkerColor(1),LineColor(1), MarkerStyle(style) );
  
  
  ///////////////////////////////////////////// FITTING /////////////////////////////
  xaxis.setRange("interesting",fit_cuts_min,fit_cuts_max);
  cout<<endl<<endl<<endl<<endl<<endl ;
  cout<<"======================================================================="<<endl ;
  RooFitResult * fit = eff.fitTo(*dataSet,ConditionalObservables(xaxis),Range("interesting"),Minos(kFALSE), Save(kTRUE),SumW2Error(kTRUE)); //symmetric errors
  //RooFitResult * fit = eff.fitTo(*dataSet,ConditionalObservables(xaxis),Range("interesting"),Minos(kTRUE), Save(kTRUE),SumW2Error(kTRUE)); //asymmetric errors
  cout<<"========================================================================"<<endl ;
  cout << "Results fit" << endl ;
  getTurnon(fit, barrel, MC, nbdata) ;
  cout << " #######  fit data done" << endl;
  cout << "----------------------" << endl ;
    
  //Compute average efficienc over a pT cut
  double massMin=45, massMax=70 ;
  getIntegralEfficiency(treeTnP, pTmin, barrel, endcap, massMin, massMax) ;

  cin >> pause ;


  //draw
  cb.plotOn(frame,LineColor(1),LineWidth(2));


  
  ////////////////////////////  DRAWING PLOTS AND LEGENDS /////////////////////////////////
  // _____________________________________________________________________________________

  TCanvas* ca = new TCanvas("ca","Trigger Efficiency") ;

  ca->SetGridx();
  ca->SetGridy();
  ca->cd();
  
  //gPad->SetLogx();
  gPad->SetObjectStat(1);

  frame->GetYaxis()->SetRangeUser(0,1.05);
  frame->GetXaxis()->SetRangeUser(0.5,50.);
  frame->GetYaxis()->SetTitle("Efficiency");
  frame->Draw() ;
   
  leg = new TLegend(0.15604,0.684965,0.567114,0.851399,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.0297203);
  leg->SetLineColor(0);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  if (MC) {
    leg->AddEntry("NULL","CMS Preliminary 2011, #sqrt{s} = 7 TeV","h");
    leg->AddEntry("NULL","MC","h");
  } 
  else {
    if (dataset=="2011A" && PtThresh == 10){
      leg->AddEntry("NULL","CMS Preliminary 2011A Data","h");
      leg->AddEntry("NULL","#int L dt = 153.9 pb^{-1}","h");
    }
    if (dataset=="2011A" && PtThresh == 15){
      leg->AddEntry("NULL","CMS Preliminary 2011A Data","h");
      leg->AddEntry("NULL","#int L dt = 1.739 fb^{-1}","h");
    }
    if (dataset=="2011A" && PtThresh == 20 && loose){
      leg->AddEntry("NULL","CMS Preliminary 2011A Data","h");
      leg->AddEntry("NULL","#int L dt = 2.001 fb^{-1}","h");
    }
    if (dataset=="2011A" && PtThresh == 20 && medium){
      leg->AddEntry("NULL","CMS Preliminary 2011A Data","h");
      leg->AddEntry("NULL","#int L dt = 270.3 pb^{-1}","h");
    }
    if (dataset=="2011A" && PtThresh == 20 && tight){
      leg->AddEntry("NULL","CMS Preliminary 2011A Data","h");
      leg->AddEntry("NULL","#int L dt = 997.4 pb^{-1}","h");
    }


    if (dataset=="2011AMay10ReReco"){
      leg->AddEntry("NULL","CMS Preliminary 2011A Data","h");// 2011AMay10ReReco Data","h");
      leg->AddEntry("NULL","#int L dt = 153.9 pb^{-1}","h");
    }
    if (dataset=="2011APromptRecov4" ){
      leg->AddEntry("NULL","CMS Preliminary 2011APromptRecov4 Data","h");
      leg->AddEntry("NULL","#int L dt = 858.2 pb^{-1}","h");
    }
    if (dataset=="2011AAug05ReReco" ){
      leg->AddEntry("NULL","CMS Preliminary 2011AAug05ReReco Data","h");
      leg->AddEntry("NULL","#int L dt = 297.2 pb^{-1}","h");
    }
    if (dataset=="2011APromptRecov6" && PtThresh == 20){
      leg->AddEntry("NULL","CMS Preliminary 2011APromptRecov6 Data","h");
      if(medium)leg->AddEntry("NULL","#int L dt = 270.3 pb^{-1}","h");//972.5 pb^{-1}","h");//160 pb^{-1}","h");//716.2 pb^{-1}","h");
      else leg->AddEntry("NULL","#int L dt = 700.2 pb^{-1}","h");
      //832.67 pb^{-1}","h");
    }
    if (dataset=="2011APromptRecov6" && PtThresh == 15){
      leg->AddEntry("NULL","CMS Preliminary 2011APromptRecov6 Data","h");
      leg->AddEntry("NULL","#int L dt = 429.9 pb^{-1}","h");
    }

    if (dataset=="2011BPromptRecov1" && PtThresh == 20){
      leg->AddEntry("NULL","CMS Preliminary 2011BPromptRecov1 Data","h");
      leg->AddEntry("NULL","#int L dt = 1.73 fb^{-1}","h");
    }
    
  }
  
  leg->Draw();
  
  TPaveText *paveText = new TPaveText(0.162752,0.573357,0.380872,0.673007,"brNDC");
  paveText->SetLineColor(1);
  paveText->SetTextColor(1);
  paveText->SetTextFont(42);
  paveText->SetTextSize(0.0227273);
  paveText->SetFillColor(kWhite);
  paveText->SetShadowColor(kWhite);
  if (loose && PtThresh==20) paveText->AddText("HLT_LooseTau20 ");
  if (medium && PtThresh==20) paveText->AddText("HLT_MediumTau20 ");
  if (tight && PtThresh==20) paveText->AddText("HLT_TightTau20 ");
  if (loose && PtThresh==15) paveText->AddText("HLT_LooseTau15 ");
  if (loose && PtThresh==10) paveText->AddText("HLT_LooseTau10 ");

  paveText->Draw();


  ca->Print(name_image+".pdf");
  ca->Print(name_image+".eps");
  ca->Print(name_image+".root");

  return fit ;

}


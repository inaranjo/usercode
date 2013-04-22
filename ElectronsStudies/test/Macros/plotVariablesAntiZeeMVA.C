//--------------------------------------------------------------------------------------------------
// plotVariablesAntiEMVA 
//
// Macro plotting the input Variables for the AntiZee MVA. 
// 
//
// Authors: I.Naranjo
//--------------------------------------------------------------------------------------------------
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



void plotVariable(string variable = "diTauVisMass",
		  const TString& category = "inclusive",
		  const TString& xAxisTitle = "visible mass",
		  const TString& yAxisTitle = "a.u.",
		  float xMin = 0., 
		  float xMax = 350,
		  int nBins = 100
		  )
{
  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

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

  TLegend* leg = new TLegend(0.6,0.75,0.8,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  //leg->SetHeader("#splitline{CMS Preliminary}{ #sqrt{s}=8 TeV}");


  std::vector<TH1*> histograms;

  std::vector<std::string> Types ; 
  Types.push_back(std::string("Zee"));
  Types.push_back(std::string("Zej"));
  Types.push_back(std::string("VBFH125"));
  Types.push_back(std::string("GGFH125"));

  for ( std::vector<std::string>::const_iterator Type = Types.begin();
	Type  != Types.end(); ++Type ) {

    cout<<Type->data()<<endl;
    std::string inputFileName ="";
    if(Type->data()==std::string("Zee"))
//       inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES/EleTau/nTuple_DYJets_ElecTau.root";
      inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTuple_DYJ_EToTau_ElecTau_nominal.root";
    else if(Type->data()==std::string("Zej"))
      inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTuple_DYJ_JetToTau_ElecTau_nominal.root";
    else if(Type->data()==std::string("VBFH125"))
      inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTuple_VBFH125_ElecTau_nominal.root";
    else if(Type->data()==std::string("GGFH125"))
      inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTuple_GGFH125_ElecTau_nominal.root";
    else if(Type->data()==std::string("VH125"))
      inputFileName = "/data_CMS/cms/htautau/PostMoriond/NTUPLES_SVfitFix/EleTau/nTuple_VH125_ElecTau_nominal.root";


    TFile* inputFile = new TFile (inputFileName.data(),"READ");
    if(inputFile->IsZombie()){
      cout << "No such file!" << endl;
      return;
    }
    TTree* inputTree;
//     if(Type->data()==std::string("Zee"))
//       inputTree = ((TTree*)inputFile->Get("outTreePtOrd"))->CopyTree("abs(genDecay)!=(23*15) &&  leptFakeTau"); // g/Z -> e+e- e->tau
//     else
    inputTree = (TTree*)inputFile->Get("outTreePtOrd");
    bool useMt      = true;
    string antiWcut = useMt ? "MtLeg1MVA" : "-(pZetaMVA-1.5*pZetaVisMVA)" ; 
    float antiWsgn  = useMt ? 20. :  20. ;
    float antiWsdb  = useMt ? 70. :  40. ; 
    ///// LEPT PT ///////
    TCut lpt("ptL1>24 && TMath::Abs(etaL1)<2.1");
    TCut lID("((TMath::Abs(scEtaL1)<0.80 && mvaPOGNonTrig>0.925) || (TMath::Abs(scEtaL1)<1.479 && TMath::Abs(scEtaL1)>0.80 && mvaPOGNonTrig>0.975) || (TMath::Abs(scEtaL1)>1.479 && mvaPOGNonTrig>0.985)) && nHits<0.5");
    lpt = lpt && lID;
    TCut tpt("ptL2>20 && TMath::Abs(etaL2)<2.3");

    ////// TAU ISO //////
//     TCut tiso("tightestHPSMVAWP>=0 && tightestAntiECutWP>1 && (tightestAntiEMVAWP>4 || tightestAntiEMVAWP==3)"); 
    TCut tiso("tightestHPSMVAWP>=0 && tightestAntiEMVA3WP>1");  // pass Medium MVA3 WP
    TCut ltiso("tightestHPSMVAWP>-99 && tightestAntiECutWP<1 ");
    TCut mtiso("hpsMVA>0.7");

    ////// E ISO ///////
    TCut liso("combRelIsoLeg1DBetav2<0.10");
    TCut laiso("combRelIsoLeg1DBetav2>0.20 && combRelIsoLeg1DBetav2<0.50");
    TCut lliso("combRelIsoLeg1DBetav2<0.30");

    ////// EVENT WISE //////
    TCut lveto("elecFlag==0"); //elecFlag==0
    TCut SS("diTauCharge!=0");
    TCut OS("diTauCharge==0");
    TCut pZ( Form("((%s)<%f)",antiWcut.c_str(),antiWsgn));
    TCut apZ(Form("((%s)>%f)",antiWcut.c_str(),antiWsdb));

    TCut apZ2(Form("((%s)>60 && (%s)<120)",antiWcut.c_str(),antiWcut.c_str()));
    TCut hltevent("vetoEvent==0 && pairIndex[0]<1 && HLTx==1");
    TCut hltmatch("HLTmatch==1");

    ////// CATEGORIES ///
    TCut vbf("nJets30>=2 && Mjj>500 && Deta>3.5 && isVetoInJets!=1");
    TCut boost("nJets30>0 && pt1>30 && nJets20BTagged<1 && MEtMVA>30");
    boost = boost && !vbf ;
    TCut novbf("nJets30<1 && nJets20BTagged==0");

    bool removeMtCut     = false;
    if (std::string(variable.data()) == "MtLeg1MVA")removeMtCut = true;

    TCut MtCut       = removeMtCut     ? "(etaL1<999)" : pZ;
    TCut diTauCharge = OS; 

    TCut sbin;

    TCut sbinTmp("");
    if(category=="inclusive") 
      sbinTmp = "etaL1<999";
    else if(category=="vbf")
      sbinTmp = vbf;
    else if(category=="novbf")
      sbinTmp = novbf;
    else if(category=="boost")
      sbinTmp = boost;

    sbin =  sbinTmp && lpt && tpt && tiso && liso && lveto && diTauCharge  && MtCut  && hltevent && hltmatch ;

    cout<<"MtCut : "<<MtCut<<endl;
    cout<<"Selection : "<<category.Data()<<" "<<sbin<<endl;

    TH1F* hVariable   = new TH1F( "hVariable" ,"" , nBins ,xMin, xMax);
    hVariable->SetXTitle(Form("%s",variable.data()));

    if(Type->data()==std::string("Zee"))
      hVariable->SetLineColor(1);
    else if(Type->data()==std::string("Zej"))
      hVariable->SetLineColor(2);
    else if(Type->data()==std::string("VBFH125"))
      hVariable->SetLineColor(3);
    else if(Type->data()==std::string("GGFH125"))
      hVariable->SetLineColor(4);

    hVariable->SetLineWidth(2);

    inputTree->Draw(Form("%s>>hVariable",variable.data()));

    cout<<"Variable plotted : "<<variable<<endl;
    cout<<"  Total number of Candidates : "<<hVariable->GetEntries()<<endl;
    inputTree->Draw(Form("%s>>hVariable",variable.data()),"(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec*weightHepNup*ZeeWeight)"*sbin);
    cout<<"  Number of Cantidates after selection: "<<hVariable->GetEntries()<<endl;
    hVariable->Scale(1./hVariable->Integral());
    leg->AddEntry(hVariable,Form("%s",Type->data()));

    histograms.push_back(hVariable);
    c1->Clear();
  }

  TH1* refHistogram = histograms.front();
  refHistogram->SetStats(false);
  refHistogram->SetTitle("");

  TAxis* xAxis = refHistogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.Data());
  xAxis->SetTitleOffset(1.15);
  TAxis* yAxis = refHistogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.Data());
  yAxis->SetTitleOffset(1.30);

  int numHistograms = histograms.size();
  float YMax = 0;
  for ( int iHistogram = 0; iHistogram < numHistograms; ++iHistogram ) {
    TH1* histogram = histograms[iHistogram];
    if(histogram->GetMaximum()>YMax) YMax = histogram->GetMaximum();
  }
  for ( int iHistogram = 0; iHistogram < numHistograms; ++iHistogram ) {
    TH1* histogram = histograms[iHistogram];
    yAxis->SetRangeUser(0.,YMax+0.10*YMax);
    std::string drawOption = "hist";
    if ( iHistogram > 0 ) drawOption.append("same");
    histogram->Draw(drawOption.data());
    leg->Draw();

  }//loop matchings
  string outputName = Form("plots/plotVariablesAntiZeeMVA/%s/plotVariablesAntiZeeMVA_%s_%s",category.Data(),category.Data(),variable.data());
  c1->Print(std::string(outputName).append(".png").data());
  c1->Print(std::string(outputName).append(".pdf").data());

}



void plotAllVariables(){


  std::vector<std::string> categories;
  categories.push_back(std::string("inclusive"));
//   categories.push_back(std::string("vbf"));
//   categories.push_back(std::string("novbf"));
//   categories.push_back(std::string("boost"));

  std::vector<std::string> variables;
  variables.push_back(std::string("ptL1"));
  variables.push_back(std::string("scEtaL1"));
  variables.push_back(std::string("ptL2"));
  variables.push_back(std::string("etaL2"));
  variables.push_back(std::string("MEtMVA"));
  variables.push_back(std::string("diTauVisMass"));
  variables.push_back(std::string("pt1"));
  variables.push_back(std::string("pt2"));
  variables.push_back(std::string("visibleTauMass"));
  variables.push_back(std::string("dPhiL1L2"));
  variables.push_back(std::string("dPhiL1J1"));
  variables.push_back(std::string("dPhiL1J2"));
  variables.push_back(std::string("AntiEMVA3raw"));

//   variables.push_back(std::string("etaL1"));
//   variables.push_back(std::string("combRelIsoLeg1DBetav2"));
//   variables.push_back(std::string("MEtMVAPhi"));
//   variables.push_back(std::string("MEtCov00"));
//   variables.push_back(std::string("MEtCov01"));
//   variables.push_back(std::string("MEtCov10"));
//   variables.push_back(std::string("MEtCov11"));
//   variables.push_back(std::string("MtLeg1MVA"));
//   variables.push_back(std::string("MtLeg2MVA"));
//   variables.push_back(std::string("eta1"));
//   variables.push_back(std::string("eta2"));
//   variables.push_back(std::string("diTauNSVfitMass"));

  std::map<std::string, std::string> xAxisTitles;
  xAxisTitles["ptL1"]             = "e p_{T}" ;
  xAxisTitles["scEtaL1"]          = "e sc #eta" ;
  xAxisTitles["etaL1"]            = "e #eta" ;
  xAxisTitles["ptL2"]             = "#tau p_{T}" ;
  xAxisTitles["etaL2"]            = "#tau #eta" ;
  xAxisTitles["visibleTauMass"]   = "visible #tau_{h} mass" ;
  xAxisTitles["dPhiL1L2"]         = "#Delta #phi (e #tau)" ;
  xAxisTitles["AntiEMVA3raw"]        = "antiEMVA3 (#tau)" ;
  xAxisTitles["combRelIsoLeg1DBetav2"] = "Comb Iso DB (e)" ;
  xAxisTitles["MEtMVA"]           = "MET" ;
  xAxisTitles["MEtMVAPhi"]        = "MET #phi" ;
  xAxisTitles["MEtCov00"]         = "MET Cov 00" ;
  xAxisTitles["MEtCov01"]         = "MET Cov 01" ;
  xAxisTitles["MEtCov10"]         = "MET Cov 10" ;
  xAxisTitles["MEtCov11"]         = "MET Cov 11" ;
  xAxisTitles["MtLeg1MVA"]        = "M_{T}(e#nu)" ;
  xAxisTitles["MtLeg2MVA"]        = "M_{T}(#tau#nu)" ;
  xAxisTitles["pt1"]              = "Leading jet p_{T}" ;
  xAxisTitles["eta1"]             = "Leading jet #eta" ;
  xAxisTitles["pt2"]              = "Sub-Leading jet p_{T}" ;
  xAxisTitles["eta2"]             = "Sub-Leading jet #eta" ;
  xAxisTitles["diTauVisMass"]     = "visible mass" ;
  xAxisTitles["diTauNSVfitMass"]  = "SVFit mass" ;



  std::map<std::string, float> xMinValues;
  xMinValues["ptL1"]              = 0.;
  xMinValues["scEtaL1"]           = -2.4;
  xMinValues["etaL1"]             = -2.4;
  xMinValues["ptL2"]              = 0.;
  xMinValues["etaL2"]             = -2.4;
  xMinValues["visibleTauMass"]    = 0.;
  xMinValues["dPhiL1L2"]          = 0.;
  xMinValues["AntiEMVA3raw"]= 0.;
  xMinValues["combRelIsoLeg1DBetav2"]= 0.;
  xMinValues["MEtMVA"]            = 0.;
  xMinValues["MEtMVAPhi"]         = -3.2;
  xMinValues["MEtCov00"]          = 0.;
  xMinValues["MEtCov01"]          = -200.;
  xMinValues["MEtCov10"]          = -200.;
  xMinValues["MEtCov11"]          = 0.;
  xMinValues["MtLeg1MVA"]         = 0.;
  xMinValues["MtLeg2MVA"]         = 0.;
  xMinValues["pt1"]               = 20.;
  xMinValues["eta1"]              = -4.5;
  xMinValues["pt2"]               = 20.;
  xMinValues["eta2"]              = -4.5;
  xMinValues["diTauVisMass"]      = 0.;
  xMinValues["diTauNSVfitMass"]   = 0.;



  std::map<std::string, float> xMaxValues;
  xMaxValues["ptL1"]              = 120.;
  xMaxValues["scEtaL1"]           = 2.4;
  xMaxValues["etaL1"]             = 2.4;
  xMaxValues["ptL2"]              = 120.;
  xMaxValues["etaL2"]             = 2.4;
  xMaxValues["visibleTauMass"]    = 2.;
  xMaxValues["dPhiL1L2"]          = 3.2;
  xMaxValues["AntiEMVA3raw"]= 1.;
  xMaxValues["combRelIsoLeg1DBetav2"]= 0.1;
  xMaxValues["MEtMVA"]            = 100.;
  xMaxValues["MEtMVAPhi"]         = 3.2;
  xMaxValues["MEtCov00"]          = 600.;
  xMaxValues["MEtCov01"]          = 200.;
  xMaxValues["MEtCov10"]          = 200.;
  xMaxValues["MEtCov11"]          = 600.;
  xMaxValues["MtLeg1MVA"]         = 160.;
  xMaxValues["MtLeg2MVA"]         = 160.;
  xMaxValues["pt1"]               = 300.;
  xMaxValues["eta1"]              = 4.5;
  xMaxValues["pt2"]               = 300.;
  xMaxValues["eta2"]              = 4.5;
  xMaxValues["diTauVisMass"]      = 200;
  xMaxValues["diTauNSVfitMass"]   = 300;



  std::map<std::string, int> nBins;
  nBins["ptL1"]                   = 30.;
  nBins["scEtaL1"]                = 36.;
  nBins["etaL1"]                  = 36.;
  nBins["ptL2"]                   = 30.;
  nBins["etaL2"]                  = 36.;
  nBins["visibleTauMass"]         = 40.;
  nBins["dPhiL1L2"]               = 32;
  nBins["AntiEMVA3raw"] = 32;
  nBins["combRelIsoLeg1DBetav2"]  = 20.;
  nBins["MEtMVA"]                 = 20.;
  nBins["MEtMVAPhi"]              = 32;
  nBins["MEtCov00"]               = 60.;
  nBins["MEtCov01"]               = 40.;
  nBins["MEtCov10"]               = 40.;
  nBins["MEtCov11"]               = 60.;
  nBins["MtLeg1MVA"]              = 40.;
  nBins["MtLeg2MVA"]              = 40.;
  nBins["pt1"]                    = 27.;
  nBins["eta1"]                   = 25;
  nBins["pt2"]                    = 27.;
  nBins["eta2"]                   = 25;
  nBins["diTauVisMass"]           = 20;
  nBins["diTauNSVfitMass"]        = 30;


 
  for ( std::vector<std::string>::const_iterator category = categories.begin();
	category != categories.end(); ++category ) {
    for ( std::vector<std::string>::const_iterator variable = variables.begin();
	  variable!= variables.end(); ++variable ) {
      plotVariable(variable->data(), category->data(), xAxisTitles[*variable].data(),"a.u.",xMinValues[*variable], xMaxValues[*variable], nBins[*variable]);
    }
  }
}





//--------------------------------------------------------------------------------------------------
// plotVariablesAntiEMVA 
//
// Macro plotting the input Variables for the AntiElectron MVA. 
// It takes them from the ouput tree of the AntiEMVA analyzer.
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



void plotVariable(string variable = "Elec_Fbrem",
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
		  const TString& Region = "Endcap"
		   )
{
   string discriminator = "";
//   string discriminator = "-AntiEMed";

  float AbsEtaMin = 0; 
  float AbsEtaMax = 3.0;
  if(Region == "Barrel"){
    AbsEtaMin = 0; 
    AbsEtaMax = 1.479;
  }
  if(Region == "Endcap"){
    AbsEtaMin = 1.479; 
    AbsEtaMax = 3.0;
  }
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
  //leg->SetHeader("#splitline{CMS Preliminary}{ #sqrt{s}=7 TeV}");

//   std::string inputFileName = "/data_CMS/cms/ivo/AntiEMVA/Trees/AntiEMVA_Fall11DYJetsToLL-iter4.root";
//   std::string inputFileName = "/data_CMS/cms/ivo/AntiEMVA/Trees/Trees_ForV4/AntiEMVA_AntiEMVATrees-DYJetsToLL-madgraph-PUS6.root";
  std::string inputFileName = "/data_CMS/cms/ivo/AntiEMVA/Trees/Trees_ForV4/AntiEMVA_V4.root";
  TFile* inputFile = new TFile (inputFileName.data(),"READ");
  if(inputFile->IsZombie()){
    cout << "No such file!" << endl;
    return;
  }
  TTree* inputTree = (TTree*)inputFile->Get("AntiEMVAAnalyzer2/tree");
//   TTree* inputTree = (TTree*)inputFile->Get("AntiEMVAAnalyzer/tree");
  std::vector<TH1*> histograms;

  std::vector<std::string> matchings ; 
  matchings.push_back("GenHadMatch");
  matchings.push_back("GenEleMatch");

  for ( std::vector<std::string>::const_iterator matching = matchings.begin();
	matching  != matchings.end(); ++matching ) {


    TCut PUSelection(Form("NumPV>%i && NumPV<%i",numPVMin,numPVMax));
    TCut ElecPtSelection (Form("Elec_Pt>%0f && Elec_Pt<%0f",PtMin,PtMax));
    TCut TauPtSelection (Form("Tau_Pt>%0f && Tau_Pt<%0f",PtMin,PtMax));
    TCut ElecAbsEtaSelection (Form("Elec_AbsEta>%0f && Elec_AbsEta<%0f",AbsEtaMin,AbsEtaMax));
    TCut TauAbsEtaSelection = "";
    if(Region == "Barrel"){
      TauAbsEtaSelection = "Tau_Eta>-1.479 && Tau_Eta<1.479";
    }
    if(Region == "Endcap"){
      TauAbsEtaSelection = "(Tau_Eta>1.479 && Tau_Eta<3.0) || (Tau_Eta>-3.0 && Tau_Eta<-1.479)";
    }
    //   TCut TauAbsEtaSelection (Form("Tau_AbsEta>%0f && Tau_AbsEta<%0f",AbsEtaMin,AbsEtaMax));
    TCut ElecMatchSelection (Form("Elec_%s == 1",matching->data()));
    //   TCut ElecMatchSelection (Form("Elec_PFTauMatch && Elec_%s",matching->data()));
    TCut TauMatchSelection (Form("Tau_%s",matching->data()));
    TCut CategorySelection = "";
    if(discriminator == ""){
      if (category == "NoEleMatch") CategorySelection = "Tau_GsfEleMatch<0.5"; 
      if (category == "woG") CategorySelection = "Tau_NumGammaCands<0.5"; 
      if (category == "wGwoGSF") CategorySelection = "Tau_NumGammaCands>0.5 && Tau_HasGsf<0.5";
      if (category == "wGwGSFwoPFMVA")CategorySelection = "Tau_NumGammaCands>0.5 && Tau_HasGsf>0.5 && Elec_PFMvaOutput<-0.1";
      if (category == "wGwGSFwPFMVA")CategorySelection = "Tau_NumGammaCands>0.5 && Tau_HasGsf>0.5 && Elec_PFMvaOutput>-0.1";
    }

    if(discriminator == "-AntiEMed"){
      if (category == "NoEleMatch") CategorySelection = "Tau_GsfEleMatch<0.5"; 
      if (category == "woG") CategorySelection = "Tau_NumGammaCands<0.5"; 
      if (category == "wGwoGSF") CategorySelection = "Tau_NumGammaCands>0.5 && (Tau_HasGsf<0.5 || (Tau_HasGsf>0.5 && Elec_PFMvaOutput>-0.1))";
      if (category == "wGwGSFwoPFMVA")CategorySelection = "Tau_NumGammaCands>0.5 && Tau_HasGsf>0.5 && Elec_PFMvaOutput<-0.1";
    }

  TCut ElecSelection = CategorySelection && PUSelection && ElecPtSelection && ElecAbsEtaSelection && ElecMatchSelection ;
  TCut TauSelection = CategorySelection && PUSelection && TauPtSelection && TauAbsEtaSelection && TauMatchSelection ;
  TCut Selection;
  if (variable.find("Elec")!=std::string::npos)Selection = ElecSelection;
  if (variable.find("Tau")!=std::string::npos)Selection = TauSelection;
  

  TH1F* hVariable   = new TH1F( "hVariable" ,"" , nBins ,xMin, xMax);
  hVariable->SetXTitle(Form("%s",variable.data()));

  if (matching->find("EleMatch")!=std::string::npos){
//     hVariable->SetFillColor(kRed);
//     hVariable->SetFillStyle(3345);
    hVariable->SetLineColor(kRed);
    hVariable->SetLineWidth(2);
  }
  if (matching->find("HadMatch")!=std::string::npos){
//     hVariable->SetFillColor(kBlue);
//     hVariable->SetFillStyle(3354);
    hVariable->SetLineColor(kBlue);
    hVariable->SetLineWidth(2);
  }  
  inputTree->Draw(Form("%s>>hVariable",variable.data()));

  cout<<"Variable plotted : "<<variable<<endl;
  cout<<"Matching applied : "<<matching->data()<<endl;
  cout<<"  Total number of Candidates : "<<hVariable->GetEntries()<<endl;
  inputTree->Draw(Form("%s>>hVariable",variable.data()),Selection);
  cout<<"  Number of Cantidates after selection: "<<hVariable->GetEntries()<<endl;
  hVariable->Scale(1./hVariable->Integral());
  leg->AddEntry(hVariable,Form("%s",matching->data()));

  histograms.push_back(hVariable);
  c1->Clear();
  }
//   double yMin = +1.e+6;
//   double yMax = -1.e+6;
  TH1* refHistogram = histograms.front();
  refHistogram->SetStats(false);
  refHistogram->SetTitle("");
//   refHistogram->SetMinimum(yMin);
//   refHistogram->SetMaximum(yMax);


  if (xAxisTitle == "HoHplusE" ) {
    refHistogram->SetMaximum(1.0);
    refHistogram->SetMinimum(0.01);
    c1->SetLogy();
  }

  if(xAxisTitle == "E_{#gamma}/(P_{in}-P_{out})" ){
    refHistogram->SetMaximum(0.03);
    refHistogram->SetMinimum(0.0);
  }

  if(xAxisTitle == "HadrMva(#tau)" ){
    refHistogram->SetMaximum(0.25);
    refHistogram->SetMinimum(0.0);
  }

  TAxis* xAxis = refHistogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.Data());
  xAxis->SetTitleOffset(1.15);
  //if(variable.find("AbsEta")!=std::string::npos)xAxis->SetLimits(AbsEtaMin, AbsEtaMax);
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
  string outputName = Form("plots/plotVariablesAntiEMVA/%s/plotVariablesAntiEMVA_v4_%s_%s_%s",category.Data(),category.Data(),variable.data(),Region.Data());
  c1->Print(std::string(outputName).append(".png").data());
  c1->Print(std::string(outputName).append(".pdf").data());

}



void plotAllVariables(){


  std::vector<std::string> categories;
  categories.push_back(std::string("All"));
  categories.push_back(std::string("NoEleMatch"));
  categories.push_back(std::string("woG"));
  categories.push_back(std::string("wGwoGSF"));
  categories.push_back(std::string("wGwGSFwoPFMVA"));
  categories.push_back(std::string("wGwGSFwPFMVA"));
//   categories.push_back(std::string("All"));
//   categories.push_back(std::string("TauNoGammas"));
//   categories.push_back(std::string("TauHasGammasNoGsfTrack"));
//   categories.push_back(std::string("TauHasGammasHasGsfTrackPFmvaBelow01"));
//   categories.push_back(std::string("TauHasGammasHasGsfTrackPFmvaOver01"));

  std::vector<std::string> variables;
  variables.push_back(std::string("Elec_EeOverPout"));
  variables.push_back(std::string("Elec_EgammaOverPdif"));
  variables.push_back(std::string("Elec_EtotOverPin"));
  variables.push_back(std::string("Elec_PFMvaOutput"));
  variables.push_back(std::string("Elec_Logsihih"));
  variables.push_back(std::string("Elec_DeltaEta"));
  variables.push_back(std::string("Elec_HoHplusE"));
  variables.push_back(std::string("Elec_GSFTrackResol"));
  variables.push_back(std::string("Elec_GSFTracklnPt"));
  variables.push_back(std::string("Elec_GSFTrackEta"));
//   variables.push_back(std::string("Tau_AbsEta"));
  variables.push_back(std::string("Tau_Eta"));
  variables.push_back(std::string("Tau_EtaAtEcalEntrance"));
  variables.push_back(std::string("Tau_EmFraction"));
  variables.push_back(std::string("Tau_HasGsf"));
  variables.push_back(std::string("Tau_HadrEoP"));
  variables.push_back(std::string("Tau_HadrHoP"));
  variables.push_back(std::string("Tau_NumChargedCands"));
  variables.push_back(std::string("Tau_NumGammaCands"));
  variables.push_back(std::string("Tau_HadrMva")); 
  variables.push_back(std::string("Tau_GammaEnFrac"));
  variables.push_back(std::string("Tau_GammaEtaMom")); 
  variables.push_back(std::string("Tau_GammaPhiMom")); 

  variables.push_back(std::string("Tau_Pt"));
  variables.push_back(std::string("Tau_LeadHadronPt"));
  variables.push_back(std::string("Tau_VisMass"));
  variables.push_back(std::string("Elec_AbsEta"));
  variables.push_back(std::string("Elec_Pt"));
  variables.push_back(std::string("Elec_Fbrem"));
  variables.push_back(std::string("Elec_Chi2KF"));
  variables.push_back(std::string("Elec_Chi2GSF"));
  variables.push_back(std::string("Elec_EarlyBrem"));
  variables.push_back(std::string("Elec_LateBrem"));
  variables.push_back(std::string("Elec_NumHits"));

//   variables.push_back(std::string("Tau_GsfEleMatch"));


  std::map<std::string, std::string> xAxisTitles;
  xAxisTitles["Elec_AbsEta"]         = "|#eta|(e)" ;
  xAxisTitles["Elec_Fbrem"]          = "Fbrem";
  xAxisTitles["Elec_Chi2KF"]         = "Chi2KF";
  xAxisTitles["Elec_EarlyBrem"]      = "EarlyBrem";
  xAxisTitles["Elec_LateBrem"]       = "LateBrem";
  xAxisTitles["Elec_EeOverPout"]     = "E_{e}/P_{out}";
  xAxisTitles["Elec_EgammaOverPdif"] = "E_{#gamma}/(P_{in}-P_{out})";
  xAxisTitles["Elec_EtotOverPin"]    = "E_{tot}/P_{in}";
  xAxisTitles["Elec_NumHits"]        = "NHits";
  xAxisTitles["Tau_AbsEta"]          = "|#eta|(#tau)";
  xAxisTitles["Tau_Eta"]             = "#eta(#tau)";
  xAxisTitles["Tau_EtaAtEcalEntrance"]= "#eta(#tau at ECAL)";
  xAxisTitles["Tau_GammaEnFrac"]     = "GammaEnFrac(#tau)";
  xAxisTitles["Tau_EmFraction"]      = "EmFraction(#tau)";
  xAxisTitles["Tau_HasGsf"]          = "TauHasGsf";
  xAxisTitles["Tau_HadrEoP"]         = "TauHadrEoP";
  xAxisTitles["Tau_HadrHoP"]         = "TauHadrHoP";
  xAxisTitles["Tau_NumChargedCands"] = "TauNumChargedHadrCands";
  xAxisTitles["Tau_NumGammaCands"]   = "TauNumGammaCands";
  xAxisTitles["Tau_VisMass"]         = "VisMass(#tau)";
  xAxisTitles["Tau_GammaEtaMom"]     = "Gamma #Delta#eta RMS(#tau)";
  xAxisTitles["Tau_GammaPhiMom"]     = "Gamma #Delta#phi RMS(#tau)";

  xAxisTitles["Elec_Pt"]             = "P_{T}(e)" ;
  xAxisTitles["Elec_PFMvaOutput"]    = "PFmva(e)" ;
  xAxisTitles["Elec_Logsihih"]       = "Logsihih(e)" ;
  xAxisTitles["Elec_DeltaEta"]       = "#Delta #eta (SC-Track)" ;
  xAxisTitles["Elec_HoHplusE"]       = "H/H+E" ;
  xAxisTitles["Elec_Chi2GSF"]        = "Chi2GSF" ;
  xAxisTitles["Elec_GSFTrackResol"]  = "GSFTrack Resol" ;
  xAxisTitles["Elec_GSFTracklnPt"]   = "lnP_{T}(GSFTr)" ;
  xAxisTitles["Elec_GSFTrackEta"]    = "#eta(GSFTr)" ;
  xAxisTitles["Tau_Pt"]              = "P_{T}(#tau)";
  xAxisTitles["Tau_LeadHadronPt"]    = "P_{T}(Lead hadron #tau)";
  xAxisTitles["Tau_HadrMva"]         = "HadrMva(#tau)"; 


  std::map<std::string, float> xMinValues;
  xMinValues["Elec_AbsEta"]             = 0.;
  xMinValues["Elec_Fbrem"]              = -0.2;
  xMinValues["Elec_Chi2KF"]             = 0.;
  xMinValues["Elec_EarlyBrem"]          = -2;
  xMinValues["Elec_LateBrem"]           = -2;
  xMinValues["Elec_EeOverPout"]         = 0.;
  xMinValues["Elec_EgammaOverPdif"]     = 0.;
  xMinValues["Elec_EtotOverPin"]        = 0.;
  xMinValues["Elec_NumHits"]            = 0.;
  xMinValues["Tau_AbsEta"]              = 0.;
  xMinValues["Tau_Eta"]                 = -3.;
  xMinValues["Tau_EtaAtEcalEntrance"]   = -3.;
  xMinValues["Tau_GammaEnFrac"]         = 0.;
  xMinValues["Tau_EmFraction"]          = 0.;
  xMinValues["Tau_HasGsf"]              = 0.;
  xMinValues["Tau_HadrEoP"]             = 0.;
  xMinValues["Tau_HadrHoP"]             = 0.;
  xMinValues["Tau_NumChargedCands"]     = 0.;
  xMinValues["Tau_NumGammaCands"]       = 0.;
  xMinValues["Tau_VisMass"]             = 0.;
  xMinValues["Tau_GammaEtaMom"]         = 0.; 
  xMinValues["Tau_GammaPhiMom"]         = 0.; 

  xMinValues["Elec_Pt"]                 = 10.;
  xMinValues["Elec_PFMvaOutput"]        = -1.;
  xMinValues["Elec_Logsihih"]           = -13.;
  xMinValues["Elec_DeltaEta"]           = 0.;
  xMinValues["Elec_HoHplusE"]           = 0.;
  xMinValues["Elec_Chi2GSF"]            = 0.;
  xMinValues["Elec_GSFTrackResol"]      = 0.;
  xMinValues["Elec_GSFTracklnPt"]       = 0.;
  xMinValues["Elec_GSFTrackEta"]        = -3.;
  xMinValues["Tau_Pt"]                  = 10.;
  xMinValues["Tau_LeadHadronPt"]        = 10.;
  xMinValues["Tau_HadrMva"]             = -1.; 


  std::map<std::string, float> xMaxValues;
  xMaxValues["Elec_AbsEta"]             = 3;
  xMaxValues["Elec_Fbrem"]              = 1;
  xMaxValues["Elec_Chi2KF"]             = 5;
  xMaxValues["Elec_EarlyBrem"]          = 2;
  xMaxValues["Elec_LateBrem"]           = 2;
  xMaxValues["Elec_EeOverPout"]         = 4.0;
  xMaxValues["Elec_EgammaOverPdif"]     = 4;
  xMaxValues["Elec_EtotOverPin"]        = 4;
  xMaxValues["Elec_NumHits"]            = 30;
  xMaxValues["Tau_AbsEta"]              = 3;
  xMaxValues["Tau_Eta"]                 = 3;
  xMaxValues["Tau_EtaAtEcalEntrance"]   = 3;
  xMaxValues["Tau_GammaEnFrac"]         = 1;
  xMaxValues["Tau_EmFraction"]          = 1;
  xMaxValues["Tau_HasGsf"]              = 2;
  xMaxValues["Tau_HadrEoP"]             = 1;
  xMaxValues["Tau_HadrHoP"]             = 1;
  xMaxValues["Tau_NumChargedCands"]     = 5;
  xMaxValues["Tau_NumGammaCands"]       = 5;
  xMaxValues["Tau_VisMass"]             = 1.8;
  xMaxValues["Tau_GammaEtaMom"]         = 3.; 
  xMaxValues["Tau_GammaPhiMom"]         = 5.; 

  xMaxValues["Elec_Pt"]                 = 80.;
  xMaxValues["Elec_PFMvaOutput"]        = 1.;
  xMaxValues["Elec_Logsihih"]           = -2.;
  xMaxValues["Elec_DeltaEta"]           = 0.05;
  xMaxValues["Elec_HoHplusE"]           = 1.;
  xMaxValues["Elec_Chi2GSF"]            = 5.;
  xMaxValues["Elec_GSFTrackResol"]      = 1.;
  xMaxValues["Elec_GSFTracklnPt"]       = 15.;
  xMaxValues["Elec_GSFTrackEta"]        = 3.;
  xMaxValues["Tau_Pt"]                  = 150.;
  xMaxValues["Tau_LeadHadronPt"]        = 150.;
  xMaxValues["Tau_HadrMva"]             = 1.; 


  std::map<std::string, int> nBins;
  nBins["Elec_AbsEta"]             = 100;
  nBins["Elec_Fbrem"]              = 100;
  nBins["Elec_Chi2KF"]             = 50;
  nBins["Elec_EarlyBrem"]          = 5;
  nBins["Elec_LateBrem"]           = 5;
  nBins["Elec_EeOverPout"]         = 100;
  nBins["Elec_EgammaOverPdif"]     = 100;
  nBins["Elec_EtotOverPin"]        = 100;
  nBins["Elec_NumHits"]            = 30;
  nBins["Tau_AbsEta"]              = 100;
  nBins["Tau_Eta"]                 = 100;
  nBins["Tau_EtaAtEcalEntrance"]   = 100;
  nBins["Tau_GammaEnFrac"]         = 20;
  nBins["Tau_EmFraction"]          = 100;
  nBins["Tau_HasGsf"]              = 2;
  nBins["Tau_HadrEoP"]             = 20;
  nBins["Tau_HadrHoP"]             = 20;
  nBins["Tau_NumChargedCands"]     = 5;
  nBins["Tau_NumGammaCands"]       = 5;
  nBins["Tau_VisMass"]             = 34;
  nBins["Tau_GammaEtaMom"]         = 30.; 
  nBins["Tau_GammaPhiMom"]         = 50.; 

  nBins["Elec_Pt"]                 = 100.;
  nBins["Elec_PFMvaOutput"]        = 100.;
  nBins["Elec_Logsihih"]           = 50.;
  nBins["Elec_DeltaEta"]           = 100.;
  nBins["Elec_HoHplusE"]           = 100.;
  nBins["Elec_Chi2GSF"]            = 50.;
  nBins["Elec_GSFTrackResol"]      = 100.;
  nBins["Elec_GSFTracklnPt"]       = 100.;
  nBins["Elec_GSFTrackEta"]        = 100.;
  nBins["Tau_Pt"]                  = 100.;
  nBins["Tau_LeadHadronPt"]        = 100.;
  nBins["Tau_HadrMva"]             = 50.;

 
  for ( std::vector<std::string>::const_iterator category = categories.begin();
	category != categories.end(); ++category ) {
    for ( std::vector<std::string>::const_iterator variable = variables.begin();
	  variable!= variables.end(); ++variable ) {
      plotVariable(variable->data(), category->data(), xAxisTitles[*variable].data(),"a.u.",xMinValues[*variable], xMaxValues[*variable], nBins[*variable], 0, 50, 15, 150,"Barrel" );
      plotVariable(variable->data(), category->data(), xAxisTitles[*variable].data(),"a.u.",xMinValues[*variable], xMaxValues[*variable], nBins[*variable], 0, 50, 15, 150,"Endcap" );
    }
  }
}


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

#include <string>
#include <map>
#include <iostream>
#include <iomanip>

TH1* getHistogram(TFile* inputFile, const TString& directory, const TString& variable)
{  
  //std::cout << "inputFile = " << inputFile->GetName() << std::endl;

  TString histogramName = TString(directory);
  if ( !histogramName.EndsWith("/") ) histogramName.Append("/");
  histogramName.Append(variable);

  TH1* histogram = (TH1*)inputFile->Get(histogramName.Data());
  //std::cout << "histogramName = " << histogramName.Data() << ": histogram = " << histogram;
  //if ( histogram ) std::cout << ", integral = " << histogram->Integral();
  //std::cout << std::endl; 

  if ( histogram && !histogram->GetSumw2N() ) histogram->Sumw2();
  else if ( !histogram) 
    std::cerr << "Failed to load histogram = " << histogramName << " from file = " << inputFile->GetName() << " !!" << std::endl;

  if ( histogram->Integral() > 0. ) histogram->Scale(1./histogram->Integral());

  return histogram;
}

void showHistograms1d(std::vector<TH1*>& histograms, const std::vector<std::string>& legendEntries, 
		      const TString& xAxisTitle, const TString& yAxisTitle,
		      const std::string& outputFileName)
{
  if ( histograms.size() == 0 ) return;
  assert(legendEntries.size() == histograms.size());

  double yMin = +1.e+6;
  double yMax = -1.e+6;
  for ( std::vector<TH1*>::iterator histogram = histograms.begin();
	histogram != histograms.end(); ++histogram ) {
    int numBinsX = (*histogram)->GetNbinsX();
    for ( int iBinX = 1; iBinX <= numBinsX; ++iBinX ) {
      double y = (*histogram)->GetBinContent(iBinX);
      if ( y < yMin ) yMin = y;
      if ( y > yMax ) yMax = y;
    }
  }
  yMin -= 0.06*(yMax - yMin);
  yMax += 0.06*(yMax - yMin);

  TCanvas* canvas = new TCanvas("canvas", "canvas", 1200, 600);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  
  canvas->SetLeftMargin(0.12);
  canvas->SetRightMargin(0.40);
  canvas->SetBottomMargin(0.12);

  TH1* refHistogram = histograms.front();
  refHistogram->SetStats(false);
  refHistogram->SetTitle("");
  refHistogram->SetMinimum(yMin);
  refHistogram->SetMaximum(yMax);


  if (xAxisTitle == "HoHplusE" ) {
    refHistogram->SetMaximum(1.0);
    refHistogram->SetMinimum(0.01);
    canvas->SetLogy();
  }

  if(xAxisTitle == "EgammaOverPdif" ){
    refHistogram->SetMaximum(0.03);
    refHistogram->SetMinimum(0.0);
  }


  TAxis* xAxis = refHistogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.Data());
  xAxis->SetTitleOffset(1.15);

  TAxis* yAxis = refHistogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.Data());
  yAxis->SetTitleOffset(1.30);

  int colors[] = { 1, 2, 3, 4, 5, 6, 7, 15 };
  int numHistograms = histograms.size();
  if ( numHistograms > 8 ) {
    std::cerr << "<showHistograms1d>:" << std::endl;
    std::cerr << "Number of histograms must not exceed 8 !!" << std::endl;
    assert(0);
  }

  TLegend* legend = new TLegend(1.0, 0.65 - 0.06*histograms.size(), /*0.995*/0.6, /*0.59*/0.89, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);

  for ( int iHistogram = 0; iHistogram < numHistograms; ++iHistogram ) {
    TH1* histogram = histograms[iHistogram];
    histogram->SetLineColor(colors[iHistogram]);
    histogram->SetLineWidth(2);

    std::string drawOption = "hist";
    if ( iHistogram > 0 ) drawOption.append("same");
    histogram->Draw(drawOption.data());

    std::string legendEntry = legendEntries[iHistogram];
    legendEntry.append(Form(" (mean = %1.2f, rms = %1.2f)", histogram->GetMean(), histogram->GetRMS()));
    legend->AddEntry(histogram, legendEntry.data(), "l");
  }

  legend->Draw();

  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  
  delete legend;
  delete canvas;
}

void makeMVAVariablesPlots()
{


  std::string inputFileName = Form("../MVAElecTau.root");
  TFile* inputFile = new TFile(inputFileName.data());

  gROOT->SetBatch(true);

  std::vector<std::string> directories;
  directories.push_back(std::string("PFAnalyzer/All"));                         
  directories.push_back(std::string("PFAnalyzer/NumPV0to10_Barrel_Pt20to30"));  
  directories.push_back(std::string("PFAnalyzer/NumPV10to20_Barrel_Pt20to30"));  
  directories.push_back(std::string("PFAnalyzer/NumPV20to30_Barrel_Pt20to30"));  
  directories.push_back(std::string("PFAnalyzer/NumPV30to40_Barrel_Pt20to30"));  
//   directories.push_back(std::string("PFAnalyzer/NumPV0to10_Endcap_Pt20to30"));   
//   directories.push_back(std::string("PFAnalyzer/NumPV10to20_Endcap_Pt20to30"));  
//   directories.push_back(std::string("PFAnalyzer/NumPV20to30_Endcap_Pt20to30"));  
//   directories.push_back(std::string("PFAnalyzer/NumPV30to40_Endcap_Pt20to30"));  
//   directories.push_back(std::string("PFAnalyzer/NumPV0to10_Barrel_Pt30to40"));  
//   directories.push_back(std::string("PFAnalyzer/NumPV10to20_Barrel_Pt30to40")); 
//   directories.push_back(std::string("PFAnalyzer/NumPV20to30_Barrel_Pt30to40")); 
//   directories.push_back(std::string("PFAnalyzer/NumPV30to40_Barrel_Pt30to40")); 
//   directories.push_back(std::string("PFAnalyzer/NumPV0to10_Endcap_Pt30to40"));  
//   directories.push_back(std::string("PFAnalyzer/NumPV10to20_Endcap_Pt30to40"));  
//   directories.push_back(std::string("PFAnalyzer/NumPV20to30_Endcap_Pt30to40"));  
//   directories.push_back(std::string("PFAnalyzer/NumPV30to40_Endcap_Pt30to40"));  
//   directories.push_back(std::string("PFAnalyzer/NumPV0to10_Barrel_Pt40to50"));   
//   directories.push_back(std::string("PFAnalyzer/NumPV10to20_Barrel_Pt40to50"));  
//   directories.push_back(std::string("PFAnalyzer/NumPV20to30_Barrel_Pt40to50"));  
//   directories.push_back(std::string("PFAnalyzer/NumPV30to40_Barrel_Pt40to50"));  
//   directories.push_back(std::string("PFAnalyzer/NumPV0to10_Endcap_Pt40to50"));   
//   directories.push_back(std::string("PFAnalyzer/NumPV10to20_Endcap_Pt40to50"));  
//   directories.push_back(std::string("PFAnalyzer/NumPV20to30_Endcap_Pt40to50"));  
//   directories.push_back(std::string("PFAnalyzer/NumPV30to40_Endcap_Pt40to50"));  

 
  std::map<std::string, std::string> selectionLabels;
  selectionLabels["PFAnalyzer/All"]                           = "All";
  selectionLabels["PFAnalyzer/NumPV0to10_Barrel_Pt20to30"]    = "NumPV0to10_Barrel_Pt20to30";
  selectionLabels["PFAnalyzer/NumPV10to20_Barrel_Pt20to30"]   = "NumPV10to20_Barrel_Pt20to30";
  selectionLabels["PFAnalyzer/NumPV20to30_Barrel_Pt20to30"]   = "NumPV20to30_Barrel_Pt20to30";
  selectionLabels["PFAnalyzer/NumPV30to40_Barrel_Pt20to30"]   = "NumPV30to40_Barrel_Pt20to30";
  selectionLabels["PFAnalyzer/NumPV0to10_Endcap_Pt20to30"]    = "NumPV0to10_Endcap_Pt20to30";
  selectionLabels["PFAnalyzer/NumPV10to20_Endcap_Pt20to30"]   = "NumPV10to20_Endcap_Pt20to30";
  selectionLabels["PFAnalyzer/NumPV20to30_Endcap_Pt20to30"]   = "NumPV20to30_Endcap_Pt20to30";
  selectionLabels["PFAnalyzer/NumPV30to40_Endcap_Pt20to30"]   = "NumPV30to40_Endcap_Pt20to30";
  selectionLabels["PFAnalyzer/NumPV0to10_Barrel_Pt30to40"]    = "NumPV0to10_Barrel_Pt30to40";
  selectionLabels["PFAnalyzer/NumPV10to20_Barrel_Pt30to40"]   = "NumPV10to20_Barrel_Pt30to40";
  selectionLabels["PFAnalyzer/NumPV20to30_Barrel_Pt30to40"]   = "NumPV20to30_Barrel_Pt30to40";
  selectionLabels["PFAnalyzer/NumPV30to40_Barrel_Pt30to40"]   = "NumPV30to40_Barrel_Pt30to40";
  selectionLabels["PFAnalyzer/NumPV0to10_Endcap_Pt30to40"]    = "NumPV0to10_Endcap_Pt30to40";
  selectionLabels["PFAnalyzer/NumPV10to20_Endcap_Pt30to40"]   = "NumPV10to20_Endcap_Pt30to40";
  selectionLabels["PFAnalyzer/NumPV20to30_Endcap_Pt30to40"]   = "NumPV20to30_Endcap_Pt30to40";
  selectionLabels["PFAnalyzer/NumPV30to40_Endcap_Pt30to40"]   = "NumPV30to40_Endcap_Pt30to40";
  selectionLabels["PFAnalyzer/NumPV0to10_Barrel_Pt40to50"]    = "NumPV0to10_Barrel_Pt40to50";
  selectionLabels["PFAnalyzer/NumPV10to20_Barrel_Pt40to50"]   = "NumPV10to20_Barrel_Pt40to50";
  selectionLabels["PFAnalyzer/NumPV20to30_Barrel_Pt40to50"]   = "NumPV20to30_Barrel_Pt40to50";
  selectionLabels["PFAnalyzer/NumPV30to40_Barrel_Pt40to50"]   = "NumPV30to40_Barrel_Pt40to50";
  selectionLabels["PFAnalyzer/NumPV0to10_Endcap_Pt40to50"]    = "NumPV0to10_Endcap_Pt40to50";
  selectionLabels["PFAnalyzer/NumPV10to20_Endcap_Pt40to50"]   = "NumPV10to20_Endcap_Pt40to50";
  selectionLabels["PFAnalyzer/NumPV20to30_Endcap_Pt40to50"]   = "NumPV20to30_Endcap_Pt40to50";
  selectionLabels["PFAnalyzer/NumPV30to40_Endcap_Pt40to50"]   = "NumPV30to40_Endcap_Pt40to50";


  std::vector<std::string> histogramNames;
  histogramNames.push_back(std::string("hNumPV"));
  histogramNames.push_back(std::string("hElecAbsEta"));
  histogramNames.push_back(std::string("hElecPt"));
  histogramNames.push_back(std::string("hEtotOverPin"));
  histogramNames.push_back(std::string("hEeOverPout"));
  histogramNames.push_back(std::string("hEgammaOverPdif"));
  histogramNames.push_back(std::string("hEarlyBrem")); 
  histogramNames.push_back(std::string("hLateBrem")); 
  histogramNames.push_back(std::string("hLogsihih"));
  histogramNames.push_back(std::string("hDeltaEta"));   
  histogramNames.push_back(std::string("hHoHplusE"));
  histogramNames.push_back(std::string("hFbrem"));
  histogramNames.push_back(std::string("hChi2KF"));
  histogramNames.push_back(std::string("hChi2GSF"));
  histogramNames.push_back(std::string("hNHits"));
  histogramNames.push_back(std::string("hGSFResol"));
  histogramNames.push_back(std::string("hGSFlnPt"));
  histogramNames.push_back(std::string("hGSFEta"));
  
  histogramNames.push_back(std::string("hTauAbsEta"));
  histogramNames.push_back(std::string("hTauPt"));
  histogramNames.push_back(std::string("hTauSignalPFChargedCands"));
  histogramNames.push_back(std::string("hTauSignalPFGammaCands"));
  histogramNames.push_back(std::string("hTauLeadPFChargedHadrMva"));
  histogramNames.push_back(std::string("hTauLeadPFChargedHadrHoP"));
  histogramNames.push_back(std::string("hTauLeadPFChargedHadrEoP"));
  histogramNames.push_back(std::string("hTauHasGsf"));
  histogramNames.push_back(std::string("hTauVisMass"));
  histogramNames.push_back(std::string("hTauEmFraction"));
  histogramNames.push_back(std::string("hGammaEtaMom"));
  histogramNames.push_back(std::string("hGammaPhiMom"));
  histogramNames.push_back(std::string("hGammaEnFrac"));

  std::map<std::string, std::string> plotLabels;
  plotLabels["hNumPV"]     = "NumPV" ;
  plotLabels["hElecAbsEta"]     = "ElecAbsEta" ;
  plotLabels["hElecPt"]         = "ElecPt";
  plotLabels["hEtotOverPin"]    = "EtotOverPin";
  plotLabels["hEeOverPout"]     = "EeOverPout";
  plotLabels["hEgammaOverPdif"] = "EgammaOverPdif";
  plotLabels["hEarlyBrem"]      = "EarlyBrem";
  plotLabels["hLateBrem"]       = "LateBrem";
  plotLabels["hLogsihih"]       = "Logsihih";
  plotLabels["hDeltaEta"]       = "DeltaEta";   
  plotLabels["hHoHplusE"]            = "HoHplusE";
  plotLabels["hFbrem"]          = "Fbrem";
  plotLabels["hChi2KF"]         = "Chi2KF";
  plotLabels["hChi2GSF"]        = "Chi2GSF";
  plotLabels["hNHits"]          = "NHits";
  plotLabels["hGSFResol"]       = "GSFResol";
  plotLabels["hGSFlnPt"]        = "GSFlnPt";
  plotLabels["hGSFEta"]         = "GSFEta";
  
  plotLabels["hTauAbsEta"]               = "TauAbsEta";
  plotLabels["hTauPt"]                   = "TauPt" ;
  plotLabels["hTauSignalPFChargedCands"] = "TauSignalPFChargedCands";
  plotLabels["hTauSignalPFGammaCands"]   = "TauSignalPFGammaCands";
  plotLabels["hTauLeadPFChargedHadrMva"] = "TauLeadPFChargedHadrMva";
  plotLabels["hTauLeadPFChargedHadrHoP"] = "TauLeadPFChargedHadrHoP";
  plotLabels["hTauLeadPFChargedHadrEoP"] = "TauLeadPFChargedHadrEoP";
  plotLabels["hTauHasGsf"]               = "TauHasGsf";
  plotLabels["hTauVisMass"]              = "TauVisMass";
  plotLabels["hTauEmFraction"]           = "TauEmFraction";
  plotLabels["hGammaEtaMom"]             = "GammaEtaMom";
  plotLabels["hGammaPhiMom"]             = "GammaPhiMom";
  plotLabels["hGammaEnFrac"]             = "GammaEnFrac";

  std::map<std::string, std::string> xAxisTitles;
  xAxisTitles["hNumPV"]          = "NumPV" ;
  xAxisTitles["hElecAbsEta"]     = "GsfElectron |#eta|" ;
  xAxisTitles["hElecPt"]         = "GsfElectron P_{T}";
  xAxisTitles["hEtotOverPin"]    = "EtotOverPin";
  xAxisTitles["hEeOverPout"]     = "EeOverPout";
  xAxisTitles["hEgammaOverPdif"] = "EgammaOverPdif";
  xAxisTitles["hEarlyBrem"]      = "EarlyBrem";
  xAxisTitles["hLateBrem"]       = "LateBrem";
  xAxisTitles["hLogsihih"]       = "Logsihih";
  xAxisTitles["hDeltaEta"]       = "DeltaEta";   
  xAxisTitles["hHoHplusE"]            = "HoHplusE";
  xAxisTitles["hFbrem"]          = "Fbrem";
  xAxisTitles["hChi2KF"]         = "Chi2KF";
  xAxisTitles["hChi2GSF"]        = "Chi2GSF";
  xAxisTitles["hNHits"]          = "NHits";
  xAxisTitles["hGSFResol"]       = "GSFResol";
  xAxisTitles["hGSFlnPt"]        = "GSFlnPt";
  xAxisTitles["hGSFEta"]         = "GSFEta";
  
  xAxisTitles["hTauAbsEta"]               = "PfTau |#eta|";
  xAxisTitles["hTauPt"]                   = "PfTau P_{T}" ;
  xAxisTitles["hTauSignalPFChargedCands"] = "TauSignalPFChargedCands";
  xAxisTitles["hTauSignalPFGammaCands"]   = "TauSignalPFGammaCands";
  xAxisTitles["hTauLeadPFChargedHadrMva"] = "TauLeadPFChargedHadrMva";
  xAxisTitles["hTauLeadPFChargedHadrHoP"] = "TauLeadPFChargedHadrHoP";
  xAxisTitles["hTauLeadPFChargedHadrEoP"] = "TauLeadPFChargedHadrEoP";
  xAxisTitles["hTauHasGsf"]               = "TauHasGsf";
  xAxisTitles["hTauVisMass"]              = "TauVisMass";
  xAxisTitles["hTauEmFraction"]           = "TauEmFraction";
  xAxisTitles["hGammaEtaMom"]             = "GammaEtaMom";
  xAxisTitles["hGammaPhiMom"]             = "GammaPhiMom";
  xAxisTitles["hGammaEnFrac"]             = "GammaEnFrac";




  std::vector<std::string> matchings;
  matchings.push_back(std::string("JetMatch"));
  matchings.push_back(std::string("No_Match"));
  matchings.push_back(std::string("EleMatch"));
  matchings.push_back(std::string("HadMatch"));



  for ( std::vector<std::string>::const_iterator directory = directories.begin();
	directory != directories.end(); ++directory ) {
    for ( std::vector<std::string>::const_iterator histogramName = histogramNames.begin();
	  histogramName != histogramNames.end(); ++histogramName ) {
      std::vector<TH1*> histograms;
      std::vector<std::string> legendEntries;
      for ( std::vector<std::string>::const_iterator matching = matchings.begin();
	    matching  != matchings.end(); ++matching ) {

	std::string inputFileName = Form("/data_CMS/cms/ivo/AntiEMVA/MVAElecTau_%s.root",matching->data());
	TFile* inputFile = new TFile(inputFileName.data());

	TString directory_full = directory->data();
	TH1* histogram = getHistogram(inputFile, directory_full, histogramName->data());
	histograms.push_back(histogram);
	std::string legendEntry = Form("%s", matching->data());
	legendEntries.push_back(legendEntry);

      }
      std::string outputFileName = 
	Form("plots/makeMVAVariablesPlots_%s_%s.png", 
	     plotLabels[*histogramName].data(), selectionLabels[*directory].data());
      showHistograms1d(histograms, legendEntries, xAxisTitles[*histogramName].data(), "a.u", outputFileName);
    }
  }
    
  delete inputFile;
}




#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TMath.h>

#include <vector>
#include <iostream>
#include <iomanip>

TH1* getHistogram(TFile* inputFile, const TString& dqmDirectory, const TString& meName)
{  
  TString histogramName = dqmDirectory;
  if ( !histogramName.EndsWith("/") ) histogramName.Append("/");
  histogramName.Append(meName);
  std::cout << "histogramName = " << histogramName.Data() << std::endl; 

  TH1* histogram = (TH1*)inputFile->Get(histogramName.Data());
  std::cout << "histogram = " << histogram << std::endl; 

  if ( !histogram->GetSumw2N()     ) histogram->Sumw2();
  if (  histogram->Integral() > 0. ) histogram->Scale(1./histogram->Integral());

  return histogram;
}

TGraphAsymmErrors* makeGraph(TH1* histogramNumerator, TH1* histogramDenominator)
{
  TGraphAsymmErrors* graph = new TGraphAsymmErrors();
  graph->Divide(histogramNumerator, histogramDenominator);
  return graph;
}

struct plotEntryType
{
  plotEntryType(const TString& name)
    : name_(name),
      graphTauPt_(0),
      graphTauEta_(0),
      graphTauPhi_(0),
      graphNumVertices_(0),
      graphUnbinned_(0)
  {
    numeratorTauPt_         = bookHistogram("numeratorTauPt",         name, 20,      0.,         100.);
    numeratorTauEta_        = bookHistogram("numeratorTauEta",        name, 46,     -2.3,         +2.3);
    numeratorTauPhi_        = bookHistogram("numeratorTauPhi",        name, 36, -TMath::Pi(), +TMath::Pi());
//     numeratorNumVertices_   = bookHistogram("numeratorNumVertices",   name, 25,     -0.5,        +24.5);
    numeratorNumVertices_   = bookHistogram("numeratorNumVertices",   name, 35,     -0.5,        +34.5);
    numeratorUnbinned_      = bookHistogram("numeratorUnbinned",      name,  1,     -0.5,         +0.5);
    denominatorTauPt_       = bookHistogram("denominatorTauPt",       name, 20,      0.,         100.);
    denominatorTauEta_      = bookHistogram("denominatorTauEta",      name, 46,     -2.3,         +2.3);
    denominatorTauPhi_      = bookHistogram("denominatorTauPhi",      name, 36, -TMath::Pi(), +TMath::Pi());
//     denominatorNumVertices_ = bookHistogram("denominatorNumVertices", name, 25,     -0.5,        +24.5);
    denominatorNumVertices_ = bookHistogram("denominatorNumVertices", name, 35,     -0.5,        +34.5);
    denominatorUnbinned_    = bookHistogram("denominatorUnbinned",    name,  1,     -0.5,         +0.5);
  }
  ~plotEntryType()
  {
    delete numeratorTauPt_;
    delete numeratorTauEta_;
    delete numeratorTauPhi_;
    delete numeratorNumVertices_;
    delete numeratorUnbinned_;
    delete denominatorTauPt_;
    delete denominatorTauEta_;
    delete denominatorTauPhi_;
    delete denominatorNumVertices_;
    delete denominatorUnbinned_;
    delete graphTauPt_;
    delete graphTauEta_;
    delete graphTauPhi_;
    delete graphNumVertices_;
    delete graphUnbinned_;
  }
  TH1* bookHistogram(const TString& name1, const TString& name2, Int_t numBinsX, Float_t xMin, Float_t xMax)
  {
    TString histogramName = Form("%s_%s", name1.Data(), name2.Data());
    TH1* histogram = new TH1F(histogramName.Data(), histogramName.Data(), numBinsX, xMin, xMax);
    return histogram;
  }
  void fillHistograms(bool passesNumerator, bool passesDenominator, 
		      Float_t tauPt, Float_t tauEta, Float_t tauPhi,
		      Float_t numVertices,
		      Float_t evtWeight)
  {
    if ( passesDenominator ) {
      denominatorTauPt_->Fill(tauPt, evtWeight);
      denominatorTauEta_->Fill(tauEta, evtWeight);
      denominatorTauPhi_->Fill(tauPhi, evtWeight);
      denominatorNumVertices_->Fill(numVertices, evtWeight);
      denominatorUnbinned_->Fill(0., evtWeight);
      if ( passesNumerator ) {
	numeratorTauPt_->Fill(tauPt, evtWeight);
	numeratorTauEta_->Fill(tauEta, evtWeight);
	numeratorTauPhi_->Fill(tauPhi, evtWeight);
	numeratorNumVertices_->Fill(numVertices, evtWeight);
	numeratorUnbinned_->Fill(0., evtWeight);
      }
    }
  }
  void makeGraphs()
  {    
    graphTauPt_       = makeGraph(numeratorTauPt_,       denominatorTauPt_);
    graphTauEta_      = makeGraph(numeratorTauEta_,      denominatorTauEta_);
    graphTauPhi_      = makeGraph(numeratorTauPhi_,      denominatorTauPhi_);
    graphNumVertices_ = makeGraph(numeratorNumVertices_, denominatorNumVertices_);    
    graphUnbinned_    = makeGraph(numeratorUnbinned_,    denominatorUnbinned_); 
    std::cout << "<makeGraphs>:" << std::endl;
    Double_t x, y;
    graphUnbinned_->GetPoint(0, x, y);
    std::cout << " name = " << name_.Data() << ": numerator = " << numeratorUnbinned_->Integral() << "," 
	      << " denominator = " << denominatorUnbinned_->Integral() << " --> efficiency/fake-rate = " << y << std::endl;
  }
  TString name_;
  TH1* numeratorTauPt_;
  TH1* numeratorTauEta_;
  TH1* numeratorTauPhi_;
  TH1* numeratorNumVertices_;
  TH1* numeratorUnbinned_;
  TH1* denominatorTauPt_;
  TH1* denominatorTauEta_;
  TH1* denominatorTauPhi_;
  TH1* denominatorNumVertices_;
  TH1* denominatorUnbinned_;
  TGraphAsymmErrors* graphTauPt_;
  TGraphAsymmErrors* graphTauEta_;
  TGraphAsymmErrors* graphTauPhi_;
  TGraphAsymmErrors* graphNumVertices_;
  TGraphAsymmErrors* graphUnbinned_;
};

void showGraphs(double canvasSizeX, double canvasSizeY,
                TGraph* graph1, const std::string& legendEntry1,
                TGraph* graph2, const std::string& legendEntry2,
                TGraph* graph3, const std::string& legendEntry3,
                TGraph* graph4, const std::string& legendEntry4,
                double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, 
                std::vector<std::string>& labelTextLines, double labelTextSize,
                double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
                double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
                bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
                const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);

  canvas->SetLogy(useLogScale);

  TH1* dummyHistogram = new TH1D("dummyHistogram", "dummyHistogram", 100, xMin, xMax);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  dummyHistogram->SetMinimum(yMin);
  dummyHistogram->SetMaximum(yMax);

  TAxis* xAxis = dummyHistogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(xAxisOffset);

  TAxis* yAxis = dummyHistogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleOffset(yAxisOffset);

  dummyHistogram->Draw("axis");

  int colors[4] = { 1, 2, 3, 4 };
  int markerStyles[4] = { 20, 21, 22, 23 };

  graph1->SetLineColor(colors[0]);
  graph1->SetMarkerColor(colors[0]);
  graph1->SetMarkerStyle(markerStyles[0]);
  graph1->Draw("p");

  if ( graph2 ) {
    graph2->SetLineColor(colors[1]);
    graph2->SetMarkerColor(colors[1]);
    graph2->SetMarkerStyle(markerStyles[1]);
    graph2->Draw("p");
  }
  
  if ( graph3 ) {
    graph3->SetLineColor(colors[2]);
    graph3->SetMarkerColor(colors[2]);
    graph3->SetMarkerStyle(markerStyles[2]);
    graph3->Draw("p");
  }

  if ( graph4 ) {
    graph4->SetLineColor(colors[3]);
    graph4->SetMarkerColor(colors[3]);
    graph4->SetMarkerStyle(markerStyles[3]);
    graph4->Draw("p");
  }
  
  TLegend* legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetTextSize(legendTextSize);
  legend->AddEntry(graph1, legendEntry1.data(), "l");
  if ( graph2 ) legend->AddEntry(graph2, legendEntry2.data(), "l");
  if ( graph3 ) legend->AddEntry(graph3, legendEntry3.data(), "l");
  if ( graph4 ) legend->AddEntry(graph4, legendEntry4.data(), "l");
  legend->Draw();

  TPaveText* label = new TPaveText(labelPosX, labelPosY, labelPosX + labelSizeX, labelPosY + labelSizeY, "brNDC");
  for ( std::vector<std::string>::const_iterator labelTextLine = labelTextLines.begin();
        labelTextLine != labelTextLines.end(); ++labelTextLine ) {
    label->AddText(labelTextLine->data());
  }
  label->SetFillColor(10);
  label->SetBorderSize(0);
  label->SetTextColor(1);
  label->SetTextAlign(12);
  label->SetTextSize(labelTextSize);
  label->Draw();

  canvas->Update();
  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  
  delete dummyHistogram;
  delete label;
  delete legend;
  delete canvas;  
}

void showGraphs(double canvasSizeX, double canvasSizeY,
		plotEntryType* plots1, const std::string& legendEntry1,
		plotEntryType* plots2, const std::string& legendEntry2,
		plotEntryType* plots3, const std::string& legendEntry3,
		plotEntryType* plots4, const std::string& legendEntry4,
		double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, 
                const std::string& labelTextLine, double labelTextSize,
                double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
                bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
                const std::string& outputFileName)
{
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_base = std::string(outputFileName, 0, idx);

  std::vector<std::string> labelTextLines;
  labelTextLines.push_back(labelTextLine);

  std::string outputFileName_TauPt = std::string(outputFileName_base).append("_Pt");
  TGraph* graph1_TauPt = ( plots1 ) ? plots1->graphTauPt_ : 0;
  TGraph* graph2_TauPt = ( plots2 ) ? plots2->graphTauPt_ : 0;
  TGraph* graph3_TauPt = ( plots3 ) ? plots3->graphTauPt_ : 0;
  TGraph* graph4_TauPt = ( plots4 ) ? plots4->graphTauPt_ : 0;
  showGraphs(canvasSizeX, canvasSizeY,
	     graph1_TauPt, legendEntry1,
             graph2_TauPt, legendEntry2,
             graph3_TauPt, legendEntry3,
             graph4_TauPt, legendEntry4,
             legendTextSize, legendPosX, legendPosY, legendSizeX, legendSizeY, 
             labelTextLines, labelTextSize,
	     labelPosX, labelPosY, labelSizeX, labelSizeY,
             0., 100., "P_{T}^{#tauJet} / GeV", 1.3,
	     useLogScale, yMin, yMax, yAxisTitle, yAxisOffset,
             outputFileName_TauPt);

  std::string outputFileName_TauEta = std::string(outputFileName_base).append("_Eta");
  TGraph* graph1_TauEta = ( plots1 ) ? plots1->graphTauEta_ : 0;
  TGraph* graph2_TauEta = ( plots2 ) ? plots2->graphTauEta_ : 0;
  TGraph* graph3_TauEta = ( plots3 ) ? plots3->graphTauEta_ : 0;
  TGraph* graph4_TauEta = ( plots4 ) ? plots4->graphTauEta_ : 0;
  showGraphs(canvasSizeX, canvasSizeY,
	     graph1_TauEta, legendEntry1,
             graph2_TauEta, legendEntry2,
             graph3_TauEta, legendEntry3,
             graph4_TauEta, legendEntry4,
             legendTextSize, legendPosX, legendPosY, legendSizeX, legendSizeY, 
             labelTextLines, labelTextSize,
	     labelPosX, labelPosY, labelSizeX, labelSizeY,
             -2.3, +2.3, "#eta_{#tauJet}", 1.3,
	     useLogScale, yMin, yMax, yAxisTitle, yAxisOffset,
             outputFileName_TauEta);

  std::string outputFileName_TauPhi = std::string(outputFileName_base).append("_Phi");
  TGraph* graph1_TauPhi = ( plots1 ) ? plots1->graphTauPhi_ : 0;
  TGraph* graph2_TauPhi = ( plots2 ) ? plots2->graphTauPhi_ : 0;
  TGraph* graph3_TauPhi = ( plots3 ) ? plots3->graphTauPhi_ : 0;
  TGraph* graph4_TauPhi = ( plots4 ) ? plots4->graphTauPhi_ : 0;
  showGraphs(canvasSizeX, canvasSizeY,
	     graph1_TauPhi, legendEntry1,
             graph2_TauPhi, legendEntry2,
             graph3_TauPhi, legendEntry3,
             graph4_TauPhi, legendEntry4,
             legendTextSize, legendPosX, legendPosY, legendSizeX, legendSizeY, 
             labelTextLines, labelTextSize,
	     labelPosX, labelPosY, labelSizeX, labelSizeY,
             -TMath::Pi(), +TMath::Pi(), "#phi_{#tauJet} / Rad", 1.3,
	     useLogScale, yMin, yMax, yAxisTitle, yAxisOffset,
             outputFileName_TauPhi);

  std::string outputFileName_NumVertices = std::string(outputFileName_base).append("_NumVertices");
  TGraph* graph1_NumVertices = ( plots1 ) ? plots1->graphNumVertices_ : 0;
  TGraph* graph2_NumVertices = ( plots2 ) ? plots2->graphNumVertices_ : 0;
  TGraph* graph3_NumVertices = ( plots3 ) ? plots3->graphNumVertices_ : 0;
  TGraph* graph4_NumVertices = ( plots4 ) ? plots4->graphNumVertices_ : 0;
  showGraphs(canvasSizeX, canvasSizeY,
	     graph1_NumVertices, legendEntry1,
             graph2_NumVertices, legendEntry2,
             graph3_NumVertices, legendEntry3,
             graph4_NumVertices, legendEntry4,
             legendTextSize, legendPosX, legendPosY, legendSizeX, legendSizeY, 
             labelTextLines, labelTextSize,
	     labelPosX, labelPosY, labelSizeX, labelSizeY,
             -0.5, 34.5, "Num. reconstructed Vertices", 1.3,
	     useLogScale, yMin, yMax, yAxisTitle, yAxisOffset,
             outputFileName_NumVertices);
}

void makeTauIdEffPlots()
{
  TString inputFileName = "/data_CMS/cms/ivo/AntiEMVA/Trees/Trees_ForV4/AntiEMVA_ForTest-DYJetsToLL.root";
  TFile* inputFile = TFile::Open(inputFileName.Data());
  if ( !inputFile ) {
    std::cerr << "Failed to open input file = " << inputFileName.Data() << " !!" << std::endl;
  }

  int maxEvents = -1;

  TString treeName = "AntiEMVAAnalyzer2/tree";
  TTree* tree = dynamic_cast<TTree*>(inputFile->Get(treeName.Data()));
  if ( !tree ) {
    std::cerr << "Failed to find tree = " << treeName.Data() << " in input file = " << inputFileName.Data() << " !!" << std::endl;
  }

  ULong64_t run,event,lumi;
  Int_t NumPV;
  Int_t NumGsfEle;
  Int_t NumPFTaus;
  Int_t NumGenEle;
  Int_t NumGenHad;
  Int_t NumGenJet;

  Int_t Tau_GsfEleMatch;
  Int_t Tau_GenEleMatch;
  Int_t Tau_GenEleFromZMatch;
  Int_t Tau_GenEleFromZTauTauMatch;
  Int_t Tau_GenHadMatch;
  Int_t Tau_GenJetMatch;
  Float_t Tau_Eta;
  Float_t Tau_EtaAtEcalEntrance;
  Float_t Tau_Pt;
  Float_t Tau_LeadHadronPt;
  Float_t Tau_Phi;
  Float_t Tau_HasGsf; 
  Float_t Tau_EmFraction; 
  Float_t Tau_NumChargedCands;
  Float_t Tau_NumGammaCands; 
  Float_t Tau_HadrHoP; 
  Float_t Tau_HadrEoP; 
  Float_t Tau_VisMass; 
  Float_t Tau_GammaEtaMom;
  Float_t Tau_GammaPhiMom;
  Float_t Tau_GammaEnFrac;
  Float_t Tau_HadrMva; 

  Float_t Tau_mvaAntiEValue; 
  Float_t Tau_AntiELoose; 
  Float_t Tau_AntiEMedium; 
  Float_t Tau_AntiETight; 
  Float_t Tau_AntiEMVA; 
  Float_t Tau_AntiEMVA2; 
  Float_t Tau_AntiEMVA2_WP75; 
  Float_t Tau_AntiEMVA2_WP85; 
  Float_t Tau_AntiEMVA2_WP95; 
  Float_t Tau_AntiEMVA2_WP99; 

  Int_t Elec_GenEleMatch;
  Int_t Elec_GenEleFromZMatch;
  Int_t Elec_GenEleFromZTauTauMatch;
  Int_t Elec_GenHadMatch;
  Int_t Elec_GenJetMatch;
  Float_t Elec_AbsEta;
  Float_t Elec_Pt;
  Float_t Elec_PFMvaOutput;
  Float_t Elec_Ee;
  Float_t Elec_Egamma;
  Float_t Elec_Pin;
  Float_t Elec_Pout;
  Float_t Elec_EtotOverPin;
  Float_t Elec_EeOverPout;
  Float_t Elec_EgammaOverPdif;
  Int_t Elec_EarlyBrem;
  Int_t Elec_LateBrem;
  Float_t Elec_Logsihih;
  Float_t Elec_DeltaEta;
  Float_t Elec_HoHplusE;
  Float_t Elec_Fbrem;
  Float_t Elec_Chi2KF;
  Float_t Elec_Chi2GSF;
  Float_t Elec_NumHits;
  Float_t Elec_GSFTrackResol;
  Float_t Elec_GSFTracklnPt;
  Float_t Elec_GSFTrackEta;

  Float_t evtWeight_dummy;

  tree->SetBranchAddress("run", &run);
  tree->SetBranchAddress("event", &event);
  tree->SetBranchAddress("lumi", &lumi);
  tree->SetBranchAddress("NumPV", &NumPV );
  tree->SetBranchAddress("NumGsfEle", &NumGsfEle );
  tree->SetBranchAddress("NumPFTaus", &NumPFTaus );
  tree->SetBranchAddress("NumGenEle", &NumGenEle );
  tree->SetBranchAddress("NumGenHad", &NumGenHad );
  tree->SetBranchAddress("NumGenJet", &NumGenJet );

  tree->SetBranchAddress("Tau_GsfEleMatch", &Tau_GsfEleMatch );
  tree->SetBranchAddress("Tau_GenEleMatch", &Tau_GenEleMatch );
  tree->SetBranchAddress("Tau_GenEleFromZMatch", &Tau_GenEleFromZMatch );
  tree->SetBranchAddress("Tau_GenEleFromZTauTauMatch", &Tau_GenEleFromZTauTauMatch );
  tree->SetBranchAddress("Tau_GenHadMatch", &Tau_GenHadMatch );
  tree->SetBranchAddress("Tau_GenJetMatch", &Tau_GenJetMatch );
  tree->SetBranchAddress("Tau_Eta", &Tau_Eta );
  tree->SetBranchAddress("Tau_EtaAtEcalEntrance", &Tau_EtaAtEcalEntrance );
  tree->SetBranchAddress("Tau_Pt", &Tau_Pt );
  tree->SetBranchAddress("Tau_LeadHadronPt", &Tau_LeadHadronPt );
  tree->SetBranchAddress("Tau_Phi", &Tau_Phi );
  tree->SetBranchAddress("Tau_HasGsf", &Tau_HasGsf ); 
  tree->SetBranchAddress("Tau_EmFraction", &Tau_EmFraction ); 
  tree->SetBranchAddress("Tau_NumChargedCands", &Tau_NumChargedCands );
  tree->SetBranchAddress("Tau_NumGammaCands", &Tau_NumGammaCands ); 
  tree->SetBranchAddress("Tau_HadrHoP", &Tau_HadrHoP ); 
  tree->SetBranchAddress("Tau_HadrEoP", &Tau_HadrEoP ); 
  tree->SetBranchAddress("Tau_VisMass", &Tau_VisMass ); 
  tree->SetBranchAddress("Tau_GammaEtaMom", &Tau_GammaEtaMom );
  tree->SetBranchAddress("Tau_GammaPhiMom", &Tau_GammaPhiMom );
  tree->SetBranchAddress("Tau_GammaEnFrac", &Tau_GammaEnFrac );
  tree->SetBranchAddress("Tau_HadrMva", &Tau_HadrMva ); 

  tree->SetBranchAddress("Tau_mvaAntiEValue", &Tau_mvaAntiEValue ); 
  tree->SetBranchAddress("Tau_AntiELoose", &Tau_AntiELoose ); 
  tree->SetBranchAddress("Tau_AntiEMedium", &Tau_AntiEMedium ); 
  tree->SetBranchAddress("Tau_AntiETight", &Tau_AntiETight ); 
  tree->SetBranchAddress("Tau_AntiEMVA", &Tau_AntiEMVA ); 
  tree->SetBranchAddress("Tau_AntiEMVA2", &Tau_AntiEMVA2 ); 
  tree->SetBranchAddress("Tau_AntiEMVA2_WP75", &Tau_AntiEMVA2_WP75 ); 
  tree->SetBranchAddress("Tau_AntiEMVA2_WP85", &Tau_AntiEMVA2_WP85 ); 
  tree->SetBranchAddress("Tau_AntiEMVA2_WP95", &Tau_AntiEMVA2_WP95 ); 
  tree->SetBranchAddress("Tau_AntiEMVA2_WP99", &Tau_AntiEMVA2_WP99 );
  


  tree->SetBranchAddress("Elec_GenEleMatch", &Elec_GenEleMatch );
  tree->SetBranchAddress("Elec_GenEleFromZMatch", &Elec_GenEleFromZMatch);
  tree->SetBranchAddress("Elec_GenEleFromZTauTauMatch", &Elec_GenEleFromZTauTauMatch );
  tree->SetBranchAddress("Elec_GenHadMatch", &Elec_GenHadMatch );
  tree->SetBranchAddress("Elec_GenJetMatch", &Elec_GenJetMatch );
  tree->SetBranchAddress("Elec_AbsEta", &Elec_AbsEta );
  tree->SetBranchAddress("Elec_Pt", &Elec_Pt );
  tree->SetBranchAddress("Elec_PFMvaOutput", &Elec_PFMvaOutput );
  tree->SetBranchAddress("Elec_Ee", &Elec_Ee );
  tree->SetBranchAddress("Elec_Egamma", &Elec_Egamma );
  tree->SetBranchAddress("Elec_Pin", &Elec_Pin );
  tree->SetBranchAddress("Elec_Pout", &Elec_Pout );
  tree->SetBranchAddress("Elec_EtotOverPin", &Elec_EtotOverPin );
  tree->SetBranchAddress("Elec_EeOverPout", &Elec_EeOverPout );
  tree->SetBranchAddress("Elec_EgammaOverPdif", &Elec_EgammaOverPdif );
  tree->SetBranchAddress("Elec_EarlyBrem", &Elec_EarlyBrem );
  tree->SetBranchAddress("Elec_LateBrem", &Elec_LateBrem );
  tree->SetBranchAddress("Elec_Logsihih", &Elec_Logsihih );
  tree->SetBranchAddress("Elec_DeltaEta", &Elec_DeltaEta );
  tree->SetBranchAddress("Elec_HoHplusE", &Elec_HoHplusE );
  tree->SetBranchAddress("Elec_Fbrem", &Elec_Fbrem );
  tree->SetBranchAddress("Elec_Chi2KF", &Elec_Chi2KF );
  tree->SetBranchAddress("Elec_Chi2GSF", &Elec_Chi2GSF );
  tree->SetBranchAddress("Elec_NumHits", &Elec_NumHits );
  tree->SetBranchAddress("Elec_GSFTrackResol", &Elec_GSFTrackResol );
  tree->SetBranchAddress("Elec_GSFTracklnPt", &Elec_GSFTracklnPt );
  tree->SetBranchAddress("Elec_GSFTrackEta", &Elec_GSFTrackEta );

  plotEntryType* plotsTau_againstElectronLoose            = new plotEntryType("tau_againstElectronLoose");
  plotEntryType* plotsTau_againstElectronMedium           = new plotEntryType("tau_againstElectronMedium");
  plotEntryType* plotsTau_againstElectronTight            = new plotEntryType("tau_againstElectronTight");
  plotEntryType* plotsTau_againstElectronMVA              = new plotEntryType("tau_againstElectronMVA");
  plotEntryType* plotsTau_againstElectronVLooseMVA2      = new plotEntryType("tau_againstElectronVLooseMVA2");
  plotEntryType* plotsTau_againstElectronLooseMVA2       = new plotEntryType("tau_againstElectronLooseMVA2");
  plotEntryType* plotsTau_againstElectronMediumMVA2       = new plotEntryType("tau_againstElectronMediumMVA2");
  plotEntryType* plotsTau_againstElectronTightMVA2        = new plotEntryType("tau_againstElectronTightMVA2");
  
  plotEntryType* plotsElectron_againstElectronLoose       = new plotEntryType("e_againstElectronLoose");
  plotEntryType* plotsElectron_againstElectronMedium      = new plotEntryType("e_againstElectronMedium");
  plotEntryType* plotsElectron_againstElectronTight       = new plotEntryType("e_againstElectronTight");
  plotEntryType* plotsElectron_againstElectronMVA         = new plotEntryType("e_againstElectronMVA");
  plotEntryType* plotsElectron_againstElectronVLooseMVA2 = new plotEntryType("e_againstElectronVLooseMVA2");
  plotEntryType* plotsElectron_againstElectronLooseMVA2  = new plotEntryType("e_againstElectronLooseMVA2");
  plotEntryType* plotsElectron_againstElectronMediumMVA2  = new plotEntryType("e_againstElectronMediumMVA2");
  plotEntryType* plotsElectron_againstElectronTightMVA2   = new plotEntryType("e_againstElectronTightMVA2");

  // CV: book additional histograms
  TH1* histogramTau_leadPFChargedHadronPtDivTauPt_againstElectronTightMVA2passed = 
    new TH1F("histogramTau_leadPFChargedHadronPtDivTauPt_againstElectronVTightMVA2passed",
	     "histogramTau_leadPFChargedHadronPtDivTauPt_againstElectronVTightMVA2passed", 101, -0.005, +1.005);
//   TH1* histogramTau_leadPFCandPtDivTauPt_againstElectronVTightMVA2passed = 
//     new TH1F("histogramTau_leadPFCandPtDivTauPt_againstElectronVTightMVA2passed",
// 	     "histogramTau_leadPFCandPtDivTauPt_againstElectronVTightMVA2passed", 101, -0.005, +1.005);
  
  TH1* histogramElectron_leadPFChargedHadronPtDivTauPt_againstElectronTightMVA2passed = 
    new TH1F("histogramElectron_leadPFChargedHadronPtDivTauPt_againstElectronVTightMVA2passed",
	     "histogramElectron_leadPFChargedHadronPtDivTauPt_againstElectronVTightMVA2passed", 101, -0.005, +1.005);
//   TH1* histogramElectron_leadPFCandPtDivTauPt_againstElectronVTightMVA2passed = 
//     new TH1F("histogramElectron_leadPFCandPtDivTauPt_againstElectronVTightMVA2passed",
// 	     "histogramElectron_leadPFCandPtDivTauPt_againstElectronVTightMVA2passed", 101, -0.005, +1.005);

	  std::vector< float > Fall11Lumi ;
  Double_t Fall11Lumi_f[50] = {
    0.003388501,
    0.010357558,
    0.024724258,
    0.042348605,
    0.058279812,
    0.068851751,
    0.072914824,
    0.071579609,
    0.066811668,
    0.060672356,
    0.054528356,
    0.04919354,
    0.044886042,
    0.041341896,
    0.0384679,
    0.035871463,
    0.03341952,
    0.030915649,
    0.028395374,
    0.025798107,
    0.023237445,
    0.020602754,
    0.0180688,
    0.015559693,
    0.013211063,
    0.010964293,
    0.008920993,
    0.007080504,
    0.005499239,
    0.004187022,
    0.003096474,
    0.002237361,
    0.001566428,
    0.001074149,
    0.000721755,
    0.000470838,
    0.00030268,
    0.000184665,
    0.000112883,
    6.74043E-05,
    3.82178E-05,
    2.22847E-05,
    1.20933E-05,
    6.96173E-06,
    3.4689E-06,
    1.96172E-06,
    8.49283E-07,
    5.02393E-07,
    2.15311E-07,
    9.56938E-08
  };

  std::vector< float > Data2011Lumi;
  Double_t Data2011Lumi_f[50] = {
    0, 
    6.48087e-05, 
    0.00122942, 
    0.0107751, 
    0.0572544, 
    0.110336, 
    0.122622, 
    0.113354, 
    0.0991184, 
    0.0907195, 
    0.0813241, 
    0.0748131, 
    0.0703376, 
    0.0625598, 
    0.0487456, 
    0.030941, 
    0.0158153, 
    0.00657034, 
    0.00232069, 
    0.000782941, 
    0.000240306, 
    6.13863e-05, 
    1.36142e-05, 
    4.99913e-07, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  };

   float * weight = new float[50] ;
  for (int i=0 ; i<50 ; i++) weight[i] = Data2011Lumi_f[i]/Fall11Lumi_f[i] ;




  Int_t numEntries = tree->GetEntries();
  for ( Int_t iEntry = 0; iEntry < numEntries && (iEntry < maxEvents || maxEvents == -1); ++iEntry ) {
    if ( iEntry > 0 && (iEntry % 100000) == 0 ) {
      std::cout << "processing Event " << iEntry << std::endl;
    }

    tree->GetEntry(iEntry);

    evtWeight_dummy = 1.;
  //   if(weight[NumPV]>0)evtWeight_dummy = weight[NumPV];
//     cout<<"Weight : "<<evtWeight_dummy<<" PV : "<<NumPV<<endl;

    if ( Tau_Pt > 20. && 
	 TMath::Abs(Tau_Eta) < 2.3 ){
      if ( Tau_GenHadMatch > 0.5) {
	plotsTau_againstElectronLoose->fillHistograms(
          Tau_AntiELoose > 0.5, true, Tau_Pt, Tau_Eta, Tau_Phi, NumPV, evtWeight_dummy);
	plotsTau_againstElectronMedium->fillHistograms(
          Tau_AntiEMedium > 0.5, true, Tau_Pt, Tau_Eta, Tau_Phi, NumPV, evtWeight_dummy);
        plotsTau_againstElectronTight->fillHistograms(
          Tau_AntiETight > 0.5, true, Tau_Pt, Tau_Eta, Tau_Phi, NumPV, evtWeight_dummy);
        plotsTau_againstElectronMVA->fillHistograms(
          Tau_AntiEMVA > 0.5, true, Tau_Pt, Tau_Eta, Tau_Phi, NumPV, evtWeight_dummy);
        plotsTau_againstElectronVLooseMVA2->fillHistograms(
          Tau_AntiEMVA2_WP99 > 0.5, true, Tau_Pt, Tau_Eta, Tau_Phi, NumPV, evtWeight_dummy);
        plotsTau_againstElectronLooseMVA2->fillHistograms(
          Tau_AntiEMVA2_WP95 > 0.5, true, Tau_Pt, Tau_Eta, Tau_Phi, NumPV, evtWeight_dummy);
        plotsTau_againstElectronMediumMVA2->fillHistograms(
          Tau_AntiEMVA2_WP85 > 0.5, true, Tau_Pt, Tau_Eta, Tau_Phi, NumPV, evtWeight_dummy);
        plotsTau_againstElectronTightMVA2->fillHistograms(
          Tau_AntiEMVA2_WP75 > 0.5, true, Tau_Pt, Tau_Eta, Tau_Phi, NumPV, evtWeight_dummy);

	// fill additional histograms
	if ( Tau_AntiEMVA2_WP75 > 0.5 ) {
	  histogramTau_leadPFChargedHadronPtDivTauPt_againstElectronTightMVA2passed->Fill(
	    Tau_LeadHadronPt/Tau_Pt, evtWeight_dummy); 
//           histogramTau_leadPFCandPtDivTauPt_againstElectronVTightMVA2passed->Fill(
//             leadPFCandP4.Pt()/Tau_Pt, evtWeight_dummy); 
	}
      }
      if ( Tau_GenEleMatch > 0.5 ) {
	plotsElectron_againstElectronLoose->fillHistograms(
          Tau_AntiELoose > 0.5, true, Tau_Pt, Tau_Eta, Tau_Phi, NumPV, evtWeight_dummy);
	plotsElectron_againstElectronMedium->fillHistograms(
          Tau_AntiEMedium > 0.5, true, Tau_Pt, Tau_Eta, Tau_Phi, NumPV, evtWeight_dummy);
        plotsElectron_againstElectronTight->fillHistograms(
          Tau_AntiETight > 0.5, true, Tau_Pt, Tau_Eta, Tau_Phi, NumPV, evtWeight_dummy);
        plotsElectron_againstElectronMVA->fillHistograms(
          Tau_AntiEMVA > 0.5, true, Tau_Pt, Tau_Eta, Tau_Phi, NumPV, evtWeight_dummy);
        plotsElectron_againstElectronVLooseMVA2->fillHistograms(
          Tau_AntiEMVA2_WP99 > 0.5, true, Tau_Pt, Tau_Eta, Tau_Phi, NumPV, evtWeight_dummy);
        plotsElectron_againstElectronLooseMVA2->fillHistograms(
          Tau_AntiEMVA2_WP95 > 0.5, true, Tau_Pt, Tau_Eta, Tau_Phi, NumPV, evtWeight_dummy);
        plotsElectron_againstElectronMediumMVA2->fillHistograms(
          Tau_AntiEMVA2_WP85 > 0.5, true, Tau_Pt, Tau_Eta, Tau_Phi, NumPV, evtWeight_dummy);
        plotsElectron_againstElectronTightMVA2->fillHistograms(
          Tau_AntiEMVA2_WP75 > 0.5, true, Tau_Pt, Tau_Eta, Tau_Phi, NumPV, evtWeight_dummy);

	// fill additional histograms
	if ( Tau_AntiEMVA2_WP75 > 0.5 ) {
	  histogramElectron_leadPFChargedHadronPtDivTauPt_againstElectronTightMVA2passed->Fill(
	    Tau_LeadHadronPt/Tau_Pt, evtWeight_dummy); 
//           histogramElectron_leadPFCandPtDivTauPt_againstElectronVTightMVA2passed->Fill(
//             leadPFCandP4.Pt()/Tau_Pt, evtWeight_dummy); 
	}
	
      }
    }
  }


  plotsTau_againstElectronLoose->makeGraphs();
  plotsTau_againstElectronMedium->makeGraphs();
  plotsTau_againstElectronTight->makeGraphs();
  plotsTau_againstElectronMVA->makeGraphs();
  plotsTau_againstElectronVLooseMVA2->makeGraphs();
  plotsTau_againstElectronLooseMVA2->makeGraphs();
  plotsTau_againstElectronMediumMVA2->makeGraphs();
  plotsTau_againstElectronTightMVA2->makeGraphs(); 

  plotsElectron_againstElectronLoose->makeGraphs();
  plotsElectron_againstElectronMedium->makeGraphs();
  plotsElectron_againstElectronTight->makeGraphs();
  plotsElectron_againstElectronMVA ->makeGraphs();
  plotsElectron_againstElectronVLooseMVA2->makeGraphs();
  plotsElectron_againstElectronLooseMVA2->makeGraphs();
  plotsElectron_againstElectronMediumMVA2->makeGraphs();
  plotsElectron_againstElectronTightMVA2->makeGraphs();

  showGraphs(800, 600,
	     plotsTau_againstElectronLoose,  "againstElectronLoose",
	     plotsTau_againstElectronMedium, "againstElectronMedium",
	     plotsTau_againstElectronTight,  "againstElectronTight",
	     plotsTau_againstElectronMVA,    "againstElectronMVA",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots2_Taus_old.eps");

  showGraphs(800, 600,
	     plotsTau_againstElectronVLooseMVA2, "againstElectronVLooseMVA2",
	     plotsTau_againstElectronLooseMVA2, "againstElectronLooseMVA2",
	     plotsTau_againstElectronMediumMVA2, "againstElectronMediumMVA2",
	     plotsTau_againstElectronTightMVA2,  "againstElectronTightMVA2",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots2_Taus_new.eps");

  showGraphs(800, 600,
	     plotsElectron_againstElectronLoose,  "againstElectronLoose",
	     plotsElectron_againstElectronMedium, "againstElectronMedium",
	     plotsElectron_againstElectronTight,  "againstElectronTight",
	     plotsElectron_againstElectronMVA,    "againstElectronMVA",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow ee", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     true, 1.e-3, 1.e+1, "Fake-rate", 1.2,
	     "makeTauIdEffPlots2_Electrons_old.eps");
  showGraphs(800, 600,
	     plotsElectron_againstElectronVLooseMVA2, "againstElectronVLooseMVA2",
	     plotsElectron_againstElectronLooseMVA2, "againstElectronLooseMVA2",
	     plotsElectron_againstElectronMediumMVA2, "againstElectronMediumMVA2",
	     plotsElectron_againstElectronTightMVA2,  "againstElectronTightMVA2",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow ee", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     true, 1.e-3, 1.e+1, "Fake-rate", 1.2,
	     "makeTauIdEffPlots2_Electrons_new.eps");

  // save additional histograms
  TString outputFileName = "makeTauIdEffPlots_v4_histograms.root";
  TFile* outputFile = new TFile(outputFileName.Data(), "RECREATE");
  histogramTau_leadPFChargedHadronPtDivTauPt_againstElectronTightMVA2passed->Write();
//   histogramTau_leadPFCandPtDivTauPt_againstElectronVTightMVA2passed->Write();
  histogramElectron_leadPFChargedHadronPtDivTauPt_againstElectronTightMVA2passed->Write();
//   histogramElectron_leadPFCandPtDivTauPt_againstElectronVTightMVA2passed->Write();
  delete outputFile;

  delete inputFile;
  }

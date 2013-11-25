// #####################################################################################
// # Set of macros to quickly make some operations for the B0 --> K*0 mu+ mu- analysis #
// #####################################################################################
// # Author: Mauro Dinardo                                                             #
// #####################################################################################

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLeaf.h>
#include <TObject.h>
#include <TRandom3.h>
#include <TGraphAsymmErrors.h>
#include <TVectorD.h>
#include <TLine.h>
#include <TLatex.h>
#include <TCutG.h>
#include <TKey.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::vector;


// ####################
// # Global constants #
// ####################
#define DIREXPTCOMP "./ExperimentComparison/"

#define ordinateRange 1.0e-2

#define B0MassIntervalLeft  0.28 // [GeV/c2]
#define B0MassIntervalRight 0.28 // [GeV/c2]
#define NSigmaPsiSmall 3.0
#define NSigmaPsiBig   5.0

#define B0Mass       5.27953  // [GeV/c2]
#define JPsiMass     3.096916 // [GeV/c2]
#define PsiPrimeMass 3.68609  // [GeV/c2]

// ###########################
// # From B0 --> J/psi K* MC #
// ###########################
#define SIGMAS1JPSI 0.0272354 // [GeV/c2]
#define SIGMAS2JPSI 0.0659195 // [GeV/c2]
#define FRACJPSI    0.680522

// #############################
// # From B0 --> psi(2S) K* MC #
// #############################
#define SIGMAS1PSIP 0.026305  // [GeV/c2]
#define SIGMAS2PSIP 0.0586605 // [GeV/c2]
#define FRACPSIP    0.581335


// #######################
// # Function Definition #
// #######################
TH1D* ComputeCumulative(TH1D* hIN, string hCumulName);
void TruthMatching (string fileName, bool truthMatch);
void dBFfromGEN (string fileName);
void CompareCosMassGENRECO (string fileNameRECO, string fileNameGEN);
void ComputePileUp (string fileName);
void PlotVtxWithPileUpW (string fileNameMC, string fileNameData, unsigned int TrigCat, bool withWeights);
void PlotPileUp (string fileNameMC1, string fileNameMC2);
void PlotCutScans (string fileName, string type);
void PlotEffPlots (string fileName, unsigned int plotN, unsigned int binN);
void PlotB0vsMuMu (string fileName, bool rejectPsi);
void PlotBkgMC (string fileName, bool iFit, double scaleMCdata);
void ReduceTree (string fileNameIn, string fileNameOut);
void SampleMCforPileup (string fileNameIn, string fileNameOut);
void DivideNTuple (string fileNameIn, string fileNameOut, unsigned int n);
void SampleNTuple (string fileNameIn, string fileNameOut, double fraction);
void ComputeMCfilterEff (string fileName);
void DrawString (double Lumi);
TCutG* DrawExclusion (double Xlow, double Xhigh, double Ylow, double Yhigh, string cutName, unsigned int fillStyle, unsigned int color);
void printData (int nBins, TVectorD V1, TVectorD V2, TVectorD V3, TVectorD V4, TVectorD V5, TVectorD V6);
void offsetData (int nBins, TVectorD* V1, TVectorD* V2, TVectorD* V3, double offset);
TGraphAsymmErrors* readData (TString fileName, int dataType, int color, int markerType, bool doFill, int fillStyle, bool noHbar, double offset);
void showData (int dataType, double offset, bool noHbar);


// ###########################
// # Function Implementation #
// ###########################


// ######################################################################
// # Sub-program to compute the cumulative distribution of an histogram #
// ######################################################################
TH1D* ComputeCumulative(TH1D* hIN, string hCumulName)
{
  TH1D* hCumul = (TH1D*)hIN->Clone(hCumulName.c_str());
  for (int i = 1; i <= hIN->GetNbinsX(); i++)
    {
      hCumul->SetBinContent(i,0.0);
      hCumul->SetBinError(i,0.0);
    }      

  for (int i = 1; i <= hIN->GetNbinsX(); i++)
    for (int j = i; j <= hIN->GetNbinsX(); j++)
      hCumul->SetBinContent(j, hCumul->GetBinContent(j) + hIN->GetBinContent(i));

  cout << "Maximum of comulative: " << hCumul->GetMaximum() << endl;
  return hCumul;
}


// ###################################################################
// # Sub-program script to check the goodness of truthMatching on MC #
// ###################################################################
void TruthMatching (string fileName, bool truthMatch)
// #########################################
// # Use file with single candidate events #
// #########################################
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1112);
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetTitleOffset(1.25,"y"); 
  TGaxis::SetMaxDigits(3);


  // #################
  // # Read the tree #
  // #################
  TFile* _file0 = TFile::Open(fileName.c_str(),"READ");
  TTree* B0KstMuMuNTuple = (TTree*)_file0->Get("B0KstMuMu/B0KstMuMuNTuple");


  double minX = 4.0;
  double maxX = 6.4;
  unsigned int nBins;
  if (truthMatch == true) nBins = 300;
  else                    nBins = 40;

  TH1D* hb = new TH1D("hb","hb",nBins,minX,maxX);
  hb->SetXTitle("M(K #pi #mu^{+} #mu^{-}) (GeV)");
  hb->SetYTitle("Entries");

  TH1D* hbar = new TH1D("hbar","hbar",nBins,minX,maxX);
  hbar->SetXTitle("M(K #pi #mu^{+} #mu^{-}) (GeV)");
  hbar->SetYTitle("Entries");

  TF1 *f0, *f1;


  if (truthMatch == true)
    {
      B0KstMuMuNTuple->Draw("bMass>>hb","genSignal == 1 && truthMatchSignal == 1 && genSignHasFSR == 0");
      B0KstMuMuNTuple->Draw("bBarMass>>hbar","genSignal == 2 && truthMatchSignal == 1 && genSignHasFSR == 0");
      hb->Add(hbar);

      f0 = new TF1("f0","[3]*TMath::Gaus(x,[0],[1]) + [4]*TMath::Gaus(x,[0],[2])",minX,maxX);

      f0->SetParName(0,"G-Mean");
      f0->SetParName(1,"G-Sigma1");
      f0->SetParName(2,"G-Sigma2");
      f0->SetParName(3,"G-Ampli1");
      f0->SetParName(4,"G-Ampli2");

      f0->SetParameter(0,5.28);
      f0->SetParameter(1,0.03);
      f0->SetParameter(2,0.06);
      f0->SetParameter(3,20000.0);
      f0->SetParameter(4,10000.0);
    }
  else
    {
      B0KstMuMuNTuple->Draw("bMass>>hb","genSignal == 1 && truthMatchSignal == 0 && genSignHasFSR == 0");
      B0KstMuMuNTuple->Draw("bBarMass>>hbar","genSignal == 2 && truthMatchSignal == 0 && genSignHasFSR == 0");
      hb->Add(hbar);

      f0 = new TF1("f0","[2]*TMath::Gaus(x,[0],[1]) + [4]*TMath::Gaus(x,[0],[3]) + ([5]+[6]*x)",minX,maxX);

      f0->SetParName(0,"G-Mean1");
      f0->SetParName(1,"G-Sigma1");
      f0->SetParName(2,"G-Ampli1");

      f0->SetParName(3,"G-Sigma2");
      f0->SetParName(4,"G-Ampli2");

      f0->SetParName(5,"P-Offset");
      f0->SetParName(6,"P-Ampli");

      f0->SetParameter(0,5.28);
      f0->SetParameter(1,0.4);
      f0->SetParameter(2,90.0);

      f0->SetParameter(3,0.2);
      f0->SetParameter(4,60.0);

      f0->SetParameter(5,20.0);
      f0->SetParameter(6,-1.0);
    }


  hb->Fit("f0","0");
  hb->Draw();
  hb->GetFunction("f0")->SetLineColor(kBlue);
  hb->GetFunction("f0")->Draw("same");

  cout << "Integral of the full fit function: " << f0->Integral(minX,maxX)/((maxX - minX) / ((double)nBins)) << "\tEntries: " << hb->GetEntries() << endl;


  if (truthMatch == true)
    {
      f1 = new TF1("f1","[3]*TMath::Gaus(x,[0],[1]) + [4]*TMath::Gaus(x,[0],[2])",minX,maxX);

      f1->SetParameter(0,f0->GetParameter(0));
      f1->SetParameter(1,f0->GetParameter(1));
      f1->SetParameter(2,f0->GetParameter(2));
      f1->SetParameter(3,f0->GetParameter(3));
      f1->SetParameter(4,f0->GetParameter(4));
    }
  else
    {
      f1 = new TF1("f1","[2]*TMath::Gaus(x,[0],[1]) + [4]*TMath::Gaus(x,[0],[3])",minX,maxX);
 
      f1->SetParameter(0,f0->GetParameter(0));
      f1->SetParameter(1,f0->GetParameter(1));
      f1->SetParameter(2,f0->GetParameter(2));
      f1->SetParameter(3,f0->GetParameter(3));
      f1->SetParameter(4,f0->GetParameter(4));
    }


  cout << "Integral of the signal: " << f1->Integral(minX,maxX)/((maxX - minX) / ((double)nBins)) << endl;
}


// #########################################################
// # Sub-program to make the dBF/dq2 vs d2 plot for GEN-MC #
// #########################################################
void dBFfromGEN (string fileName)
{
  TFile* _file0 = TFile::Open(fileName.c_str(),"READ");
  TList* myList = _file0->GetListOfKeys();
  TH1D* h1 = (TH1D*)((dynamic_cast<TKey*>(myList->At(0)))->ReadObj());


  double GENMClumi = 294500.0;

  // #########
  // # J/psi #
  // #########
  double BFPsi = 7.9e-5;
  double scale = 9883.6 / GENMClumi * BFPsi / 1e-7; // Lumi ratio * BF[B0 --> J/psi or psi(2S) (mu mu) K*0 (K pi)]

  // ###########
  // # psi(2S) #
  // ###########
  // double BFPsi = 46.97e-7;
  // double scale = 9883.6 / GENMClumi * BFPsi / 1e-7; // Lumi ratio * BF[B0 --> J/psi or psi(2S) (mu mu) K*0 (K pi)]


  h1->SetBinContent(1,2888976.0 / 9122980.0 * scale / h1->GetBinWidth(1));
  h1->SetBinContent(2,5760931.0 / 9122980.0 * scale / h1->GetBinWidth(2));
  h1->SetBinContent(3,12899167.0 / 9122980.0 * scale / h1->GetBinWidth(3));
  h1->SetBinContent(4,4833302.0 / 9122980.0 * scale / h1->GetBinWidth(5));
  h1->SetBinContent(5,10130407.0 / 9122980.0 * scale / h1->GetBinWidth(5));
  h1->SetBinContent(6,4580382.0 / 9122980.0 * scale / h1->GetBinWidth(7));
  h1->SetBinContent(7,6009085.0 / 9122980.0 * scale / h1->GetBinWidth(7));
  h1->SetBinContent(8,6363197.0 / 9122980.0 * scale / h1->GetBinWidth(8));
  
  h1->SetBinError(1,h1->GetBinContent(1) * sqrt(1./2888976.0 + 1./9122980.0));
  h1->SetBinError(2,h1->GetBinContent(2) * sqrt(1./5760931.0 + 1./9122980.0));
  h1->SetBinError(3,h1->GetBinContent(3) * sqrt(1./12899167.0 + 1./9122980.0));
  h1->SetBinError(4,h1->GetBinContent(3) * sqrt(1./4833302.0 + 1./9122980.0));
  h1->SetBinError(5,h1->GetBinContent(5) * sqrt(1./10130407.0 + 1./9122980.0));
  h1->SetBinError(6,h1->GetBinContent(7) * sqrt(1./4580382.0 + 1./9122980.0));
  h1->SetBinError(7,h1->GetBinContent(7) * sqrt(1./6009085.0 + 1./9122980.0));
  h1->SetBinError(8,h1->GetBinContent(8) * sqrt(1./6363197.0 + 1./9122980.0));
  

  for (int i = 0; i < h1->GetNbinsX(); i++) cout << "--> bin #" << i+1 << "\tentry: " << h1->GetBinContent(i+1) << " +/- " << h1->GetBinError(i+1) << "\twidth: " << h1->GetBinWidth(i+1) << endl;

  TCanvas* c0 = new TCanvas("cHistoMeas","cHistoMeas",10,10,700,500);
  c0->cd();
  h1->Draw();
  c0->Update();
}


// #######################################################################################
// # Sub-program to compare cos(theta_l), cos(theta_l) and mumuMass between GEN and RECO #
// #######################################################################################
void CompareCosMassGENRECO (string fileNameRECO, string fileNameGEN)
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1112);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetTitleOffset(1.25,"y"); 
  TGaxis::SetMaxDigits(3);


  TFile* NtplFileIn1 = TFile::Open(fileNameRECO.c_str(),"READ");
  TTree* theTreeIn1 = (TTree*)NtplFileIn1->Get("B0KstMuMu/B0KstMuMuNTuple");

  TFile* NtplFileIn2 = TFile::Open(fileNameGEN.c_str(),"READ");
  TTree* theTreeIn2 = (TTree*)NtplFileIn2->Get("B0KstMuMu/B0KstMuMuNTuple");

  cout << "n entry: " << theTreeIn1->GetEntries() << endl;
  cout << "n entry: " << theTreeIn2->GetEntries() << endl;

  TCanvas* c0 = new TCanvas("c0","c0",10,10,700,500);
  TCanvas* c1 = new TCanvas("c1","c1",10,10,700,500);
  TCanvas* c2 = new TCanvas("c2","c2",10,10,700,500);

  TH1D* h0 = new TH1D("h0","h0",400,-0.5,0.5);
  TH1D* h1 = new TH1D("h1","h1",400,-0.5,0.5);
  TH1D* h2 = new TH1D("h2","h2",400,-0.5,0.5);
  h0->SetXTitle("RECO-GEN");
  h0->SetYTitle("Entries");
  h0->SetFillColor(kAzure+6);
  h1->SetXTitle("RECO-GEN");
  h1->SetYTitle("Entries");
  h1->SetFillColor(kAzure+6);
  h2->SetXTitle("RECO-GEN");
  h2->SetYTitle("Entries");
  h2->SetFillColor(kAzure+6);

  vector<double>* vec1 = new vector<double>;
  vector<double>* vec2 = new vector<double>;

  int nEntries = theTreeIn1->GetEntries();

  for (int entry = 0; entry < nEntries; entry++)
    {
      theTreeIn1->SetBranchAddress("mumuMass", &vec1);
      theTreeIn2->SetBranchAddress("mumuMass", &vec2);

      theTreeIn1->GetEntry(entry);
      theTreeIn2->GetEntry(entry);

      TLeaf* leafK1 = theTreeIn1->GetLeaf("CosThetaKArb");
      TLeaf* leafK2 = theTreeIn2->GetLeaf("CosThetaKArb");
      h0->Fill(leafK1->GetValue()-leafK2->GetValue());

      TLeaf* leafMu1 = theTreeIn1->GetLeaf("CosThetaMuArb");
      TLeaf* leafMu2 = theTreeIn2->GetLeaf("CosThetaMuArb");
      h1->Fill(leafMu1->GetValue()-leafMu2->GetValue());

      h2->Fill((*vec1)[0]-(*vec2)[0]);

      vec1->clear();
      vec2->clear();
    }

  TF1* f0 = new TF1("f0","[3]*TMath::Gaus(x,[0],[1])+[4]*TMath::Gaus(x,[0],[2])",-1,1);

  f0->SetParName(0,"Mean");
  f0->SetParName(1,"Sigma1");
  f0->SetParName(2,"Sigma2");
  f0->SetParName(3,"Ampli1");
  f0->SetParName(4,"Ampli2");

  f0->SetParameter(0,0.0);
  f0->SetParameter(1,0.02);
  f0->SetParameter(2,0.04);
  f0->SetParameter(3,600.0);
  f0->SetParameter(4,300.0);

  c0->cd();
  h0->Draw();
  h0->Fit("f0");
  c0->Update();
  cout << "Sigma cos(theta_K): " << sqrt((f0->GetParameter(3)*f0->GetParameter(1)*f0->GetParameter(1) + f0->GetParameter(4)*f0->GetParameter(2)*f0->GetParameter(2))/(f0->GetParameter(3)+f0->GetParameter(4))) << endl;

  c1->cd();
  h1->Draw();
  h1->Fit("f0");
  c1->Update();
  cout << "Sigma cos(theta_l): " << sqrt((f0->GetParameter(3)*f0->GetParameter(1)*f0->GetParameter(1) + f0->GetParameter(4)*f0->GetParameter(2)*f0->GetParameter(2))/(f0->GetParameter(3)+f0->GetParameter(4))) << endl;

  f0->SetParameter(0,0.0);
  f0->SetParameter(1,0.01);
  f0->SetParameter(2,0.05);
  f0->SetParameter(3,12000.0);
  f0->SetParameter(4,5000.0);

  c2->cd();
  h2->Draw();
  h2->Fit("f0");
  c2->Update();
  cout << "Sigma mumuMass^2: " << sqrt((f0->GetParameter(3)*f0->GetParameter(1)*f0->GetParameter(1) + f0->GetParameter(4)*f0->GetParameter(2)*f0->GetParameter(2))/(f0->GetParameter(3)+f0->GetParameter(4))) << endl;
}


// ##########################################################
// # Sub-program to compute puleup histogram from MC ntuple #
// ##########################################################
void ComputePileUp (string fileName)
{
  TFile* NtplFileIn = TFile::Open(fileName.c_str(),"READ");
  TTree* theTreeIn = (TTree*)NtplFileIn->Get("B0KstMuMu/B0KstMuMuNTuple");

  vector<double>* vec1 = new vector<double>;
  vector<double>* vec2 = new vector<double>;

  int nEntries = theTreeIn->GetEntries();
  cout << "\n@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;

  TCanvas* c0 = new TCanvas("c0","c0",10,10,700,500);
  TH1D* pileup = new TH1D("pileup","pileup",50,0,50);

  for (int entry = 0; entry < nEntries; entry++)
    {
      theTreeIn->SetBranchAddress("bunchXingMC", &vec1);
      theTreeIn->SetBranchAddress("numInteractionsMC", &vec2);
      theTreeIn->GetEntry(entry);

      for (unsigned int it = 0; it < (*vec1).size(); it++)
	if ((*vec1)[it] == 0)
	  {
	    pileup->Fill((*vec2)[it]);
	    break;
	  }

      vec1->clear();
      vec2->clear();
    }

  c0->cd();
  pileup->Draw();
  c0->Update();

  TFile fout("PileupMC.root", "RECREATE");
  fout.cd();
  pileup->Write();
  fout.Close();
  NtplFileIn->Close();
}


// ##########################################################
// # Sub-program to plot variables with pileup re-weighting #
// ##########################################################
void PlotVtxWithPileUpW (string fileNameMC, string fileNameData, unsigned int TrigCat, bool withWeights)
// ###############################################
// # if TrigCat == 0 then use the whole file     #
// # if withWeights == Ture then use the weights #
// ###############################################
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1112);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetTitleOffset(1.25,"y"); 
  TGaxis::SetMaxDigits(3);


  int nEntries;

  TCanvas* c0 = new TCanvas("c0","c0",10,10,700,500);

  TFile* NtplFileInMC = TFile::Open(fileNameMC.c_str(),"READ");
  TTree* theTreeInMC = (TTree*)NtplFileInMC->Get("B0KstMuMu/B0KstMuMuNTuple");

  TFile* NtplFileInData = TFile::Open(fileNameData.c_str(),"READ");
  TTree* theTreeInData = (TTree*)NtplFileInData->Get("B0KstMuMu/B0KstMuMuNTuple");


  // ################
  // # Read from MC #
  // ################
  nEntries = theTreeInMC->GetEntries();
  cout << "\n@@@ Total number of events in the MC tree: " << nEntries << " @@@" << endl;

  TH1D* hMC = new TH1D("hMC","hMC",50,0,50);
  hMC->Sumw2();
  hMC->SetFillColor(kAzure+6);
  hMC->SetXTitle("Vtx[#]");
  hMC->SetYTitle("a.u.");
  hMC->GetXaxis()->SetRangeUser(0.0,30.0);
  hMC->GetYaxis()->SetRangeUser(0.0,0.14);

  for (int entry = 0; entry < nEntries; entry++)
    {
      theTreeInMC->GetEntry(entry);

      TLeaf* leaf_nVtx    = theTreeInMC->GetLeaf("recoVtxN");
      TLeaf* leaf_weight  = theTreeInMC->GetLeaf("evWeight");
      TLeaf* leaf_trigCat = theTreeInMC->GetLeaf("TrigCat");

      if ((TrigCat == 0) || (leaf_trigCat->GetValue() == TrigCat))
	{
	  if (leaf_nVtx->GetValue() > 0.0)
	    {
	      if (withWeights == true) hMC->Fill(leaf_nVtx->GetValue(), leaf_weight->GetValue());
	      else hMC->Fill(leaf_nVtx->GetValue());
	    }
	}
    }


  // ##################
  // # Read from Data #
  // ##################
  nEntries = theTreeInData->GetEntries();
  cout << "\n@@@ Total number of events in the Data tree: " << nEntries << " @@@" << endl;

  TH1D* hData = new TH1D("hData","hData",50,0,50);
  hData->Sumw2();
  hData->SetMarkerStyle(20);
  hData->SetXTitle("Vtx[#]");
  hData->SetYTitle("a.u.");
  hData->GetXaxis()->SetRangeUser(0.0,30.0);
  hData->GetYaxis()->SetRangeUser(0.0,0.14);

  for (int entry = 0; entry < nEntries; entry++)
    {
      theTreeInData->GetEntry(entry);

      TLeaf* leaf_nVtx    = theTreeInData->GetLeaf("recoVtxN");
      TLeaf* leaf_weight  = theTreeInData->GetLeaf("evWeight");
      TLeaf* leaf_trigCat = theTreeInData->GetLeaf("TrigCat");

      if ((TrigCat == 0) || (leaf_trigCat->GetValue() == TrigCat))
  	if (leaf_nVtx->GetValue() > 0.0) hData->Fill(leaf_nVtx->GetValue(), leaf_weight->GetValue());
    }


  hMC->Scale(1. / hMC->Integral());
  hData->Scale(1. / hData->Integral());

  c0->cd();

  hMC->Draw("e3");
  hData->Draw("e1p same");

  TLegend* leg = new TLegend(0.7, 0.79, 0.89, 0.89, "");
  leg->AddEntry(hMC,"MC");
  leg->AddEntry(hData,"Data");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  c0->Update();
}


// ############################################################
// # Sub-program to plot the pileup of two different MC files #
// ############################################################
void PlotPileUp (string fileNameMC1, string fileNameMC2)
// #############################
// # fileNameMC1 = original MC #
// # fileNameMC2 = sampled MC  #
// #############################
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1112);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetTitleOffset(1.25,"y"); 
  TGaxis::SetMaxDigits(3);


  int nEntries;

  TCanvas* c0 = new TCanvas("c0","c0",10,10,700,500);

  TFile* fileMC1 = TFile::Open(fileNameMC1.c_str(),"READ");
  TTree* theTreeInMC1 = (TTree*)fileMC1->Get("B0KstMuMu/B0KstMuMuNTuple");

  nEntries = theTreeInMC1->GetEntries();
  cout << "\n@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;

  TH1D* hMC1 = new TH1D("hMC1","hMC1",50,0,50);
  hMC1->Sumw2();
  hMC1->SetMarkerStyle(20);
  hMC1->SetXTitle("pileup[#]");
  hMC1->SetYTitle("a.u.");
  hMC1->GetXaxis()->SetRangeUser(0.0,35.0);
  hMC1->GetYaxis()->SetRangeUser(0.0,0.11);
  theTreeInMC1->Draw("numInteractionsMC>>hMC1","evWeight*(bunchXingMC == 0)");
  hMC1->Scale(1. / hMC1->Integral());

  TH1D* hMC3 = new TH1D("hMC3","hMC3",50,0,50);
  hMC3->Sumw2();
  hMC3->SetFillColor(kGreen-7);
  hMC3->SetXTitle("pileup[#]");
  hMC3->SetYTitle("a.u."); 
  hMC3->GetXaxis()->SetRangeUser(0.0,35.0);
  hMC3->GetYaxis()->SetRangeUser(0.0,0.11);
  theTreeInMC1->Draw("numInteractionsMC>>hMC3","bunchXingMC == 0");
  hMC3->Scale(1. / hMC3->Integral());


  TFile* fileMC2 = TFile::Open(fileNameMC2.c_str(),"READ");
  TTree* theTreeInMC2 = (TTree*)fileMC2->Get("B0KstMuMu/B0KstMuMuNTuple");

  nEntries = theTreeInMC2->GetEntries();
  cout << "\n@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;

  TH1D* hMC2 = new TH1D("hMC2","hMC2",50,0,50);
  hMC2->Sumw2();
  hMC2->SetFillColor(kAzure+6);
  hMC2->SetXTitle("pileup[#]");
  hMC2->SetYTitle("a.u.");
  hMC2->GetXaxis()->SetRangeUser(0.0,35.0);
  hMC2->GetYaxis()->SetRangeUser(0.0,0.11);
  theTreeInMC2->Draw("numInteractionsMC>>hMC2","bunchXingMC == 0");
  hMC2->Scale(1. / hMC2->Integral());


  c0->cd();
  hMC3->Draw("e3");
  hMC2->Draw("e3 same");
  hMC1->Draw("e1p same");

  TLegend* leg = new TLegend(0.6, 0.7, 0.89, 0.89, "");
  leg->AddEntry(hMC3,"Original MC");
  leg->AddEntry(hMC1,"Original MC weighted");
  leg->AddEntry(hMC2,"Sampled MC");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  c0->Update();
}


// #################################
// # Sub-program to plot cut scans #
// #################################
void PlotCutScans (string fileName, string type)
{
  stringstream myString;

  TFile* _file = TFile::Open(fileName.c_str(),"READ");

  myString.clear();
  myString.str("");
  myString << "c0_" << type;
  TCanvas* c0 = (TCanvas*)_file->Get(myString.str().c_str());

  myString.clear();
  myString.str("");
  myString << "c0_" << type << "_1";
  TPad* p0 = (TPad*)c0->GetPrimitive(myString.str().c_str());
  TH1D* h0 = (TH1D*)p0->GetPrimitive("histoR1");

  myString.clear();
  myString.str("");
  myString << "c0_" << type << "_4";
  TPad* p1 = (TPad*)c0->GetPrimitive(myString.str().c_str());
  TH1D* h1 = (TH1D*)p1->GetPrimitive("histoR4");
  h1->SetYTitle("S / #sqrt{(S+B)}");


  // ##########################
  // # Set histo layout style #
  // ##########################
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1112);
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetTitleOffset(1.25,"y"); 
  TGaxis::SetMaxDigits(3);


  TCanvas* c1 = new TCanvas("c1","c1",10,10,900,500);
  c1->Divide(2,1);
  c1->cd(1);
  h0->Draw();
  c1->cd(2);
  h1->Draw();
  c1->Update();
}


// ################################################
// # Sub-program to plot the efficiency sub-plots #
// ################################################
void PlotEffPlots (string fileName, unsigned int plotN, unsigned int binN)
// ##################################
// # plotN can go from 1 to 4       #
// # binN can go from 0 to #q2 bins #
// ##################################
{
  stringstream myString;

  TFile* _file = TFile::Open(fileName.c_str(),"READ");

  TCanvas* c0 = (TCanvas*)_file->Get("cEff");

  myString.clear();
  myString.str("");
  myString << "cEff_" << plotN;
  TPad* p0 = (TPad*)c0->GetPrimitive(myString.str().c_str());

  myString.clear();
  myString.str("");
  if      (plotN == 2) myString << "vecHcosThetaK_" << binN;
  else if (plotN == 3) myString << "vecHcosThetaL_" << binN;
  else if (plotN == 4) myString << "vecHphi_" << binN;
  else exit(1);
  TH1D* h0 = (TH1D*)p0->GetPrimitive(myString.str().c_str());
  h0->GetYaxis()->SetRangeUser(0.0,ordinateRange);


  // ##########################
  // # Set histo layout style #
  // ##########################
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TGaxis::SetMaxDigits(3);


  TCanvas* c1 = new TCanvas("c1","c1",10,10,700,500);
  c1->cd();
  h0->Draw("e1");
  c1->Update();
}


// ##################################################
// # Sub-program to plot B0 vs mu-mu invariant mass #
// ##################################################
void PlotB0vsMuMu (string fileName, bool rejectPsi)
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1112);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetTitleOffset(1.25,"y"); 
  TGaxis::SetMaxDigits(3);


  int nEntries;
  double minX = 5.0;
  double maxX = 5.56;
  double minY = 0.8;
  double maxY = 5.0;

  TCanvas* c0 = new TCanvas("c0","c0",10,10,700,500);

  TFile* file0 = TFile::Open(fileName.c_str(),"READ");
  TTree* theTree = (TTree*)file0->Get("B0KstMuMu/B0KstMuMuNTuple");

  nEntries = theTree->GetEntries();
  cout << "\n@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;

  vector<double>* vec1 = new vector<double>;
  vector<double>* vec2 = new vector<double>;

  TH2D* histo = new TH2D("histo","histo",100,minX,maxX,100,minY,maxY);
  histo->SetXTitle("M(K #pi #mu^{+} #mu^{-}) (GeV)");
  histo->SetYTitle("M(#mu^{+} #mu^{-}) (GeV)");
  histo->SetZTitle("Entries / ((0.056 GeV)x(0.042 GeV))");

  for (int entry = 0; entry < nEntries; entry++)
    {
      theTree->SetBranchAddress("mumuMass", &vec1);
      theTree->SetBranchAddress("mumuMassE", &vec2);
      theTree->GetEntry(entry);
      TLeaf* theLeaf = theTree->GetLeaf("B0MassArb");
      
      if ((rejectPsi == true) && ((((*vec1)[0] > JPsiMass - NSigmaPsiBig * (*vec2)[0]) && ((*vec1)[0] < JPsiMass + NSigmaPsiSmall * (*vec2)[0])) ||
				  (((*vec1)[0] > PsiPrimeMass - NSigmaPsiSmall * (*vec2)[0]) && ((*vec1)[0] < PsiPrimeMass + NSigmaPsiSmall * (*vec2)[0])))) continue;

      histo->Fill(theLeaf->GetValue(),(*vec1)[0]);
    }
  

  c0->cd();
  histo->Draw("gcolz");
  c0->Update();
}


// #####################################################
// # Sub-program to plot B0, B+, Bs, /\b background MC #
// #####################################################
void PlotBkgMC (string fileName, bool iFit, double scaleMCdata)
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1112);
  gStyle->SetOptStat(1001100);
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetTitleOffset(1.25,"y"); 
  TGaxis::SetMaxDigits(3);


  unsigned int nBins = 20;

  TFile* fileIn = TFile::Open(fileName.c_str(),"READ");
  TTree* theTree = (TTree*)fileIn->Get("B0KstMuMu/B0KstMuMuNTuple");
  unsigned int nEntries = theTree->GetEntries();
  cout << "\n@@@ Total number of events in the input tree: " << nEntries << " @@@" << endl;

  TCanvas* c0 = new TCanvas("c0","c0",10,10,700,500);
  TCanvas* c1 = new TCanvas("c1","c1",20,20,700,500);
  TCanvas* c2 = new TCanvas("c2","c2",30,30,700,500);

  TLegend* leg0;
  TLegend* leg1;

  TH1D* h0 = new TH1D("h0","h0",nBins,B0Mass - B0MassIntervalLeft,B0Mass + B0MassIntervalRight);
  h0->Sumw2();
  h0->SetXTitle("M(K #pi #mu^{+} #mu^{-}) (GeV)");
  h0->SetYTitle("Entries / (0.028 GeV)");
  h0->SetMarkerStyle(20);

  TH1D* h1 = new TH1D("h1","h1",nBins,B0Mass - B0MassIntervalLeft,B0Mass + B0MassIntervalRight);
  h1->Sumw2();
  h1->SetXTitle("M(K #pi #mu^{+} #mu^{-}) (GeV)");
  h1->SetYTitle("Entries / (0.028 GeV)");
  h1->SetMarkerStyle(20);

  TH1D* h2 = new TH1D("h2","h2",nBins,B0Mass - B0MassIntervalLeft,B0Mass + B0MassIntervalRight);
  h2->Sumw2();
  h2->SetXTitle("M(K #pi #mu^{+} #mu^{-}) (GeV)");
  h2->SetYTitle("Entries / (0.028 GeV)");
  h2->SetMarkerStyle(20);

  
  // ################
  // # J/psi region #
  // ################
  theTree->Draw("B0MassArb>>h0","B0MassArb > 5.27953 - 0.28 && B0MassArb < 5.27953 + 0.28 && genSignal == 0 && truthMatchSignal == 1 && mumuMass > 3.096916 - 5*mumuMassE && mumuMass < 3.096916 + 3*mumuMassE");
  h0->Scale(scaleMCdata);

  TF1* f0  = new TF1("f0", "[0]*([1]*TMath::Gaus(x,[2],[3]) + (1-[1])*TMath::Gaus(x,[2],[4])) + [5]*TMath::Gaus(x,[6],[7]) + [8]*TMath::Gaus(x,[9],[10])",B0Mass - B0MassIntervalLeft,B0Mass + B0MassIntervalRight);
  TF1* f0S = new TF1("f0S","[0]*([1]*TMath::Gaus(x,[2],[3]) + (1-[1])*TMath::Gaus(x,[2],[4]))",B0Mass - B0MassIntervalLeft,B0Mass + B0MassIntervalRight);
  f0S->SetLineColor(kBlue);
  f0S->SetLineStyle(kDashed);

  // ##########
  // # Signal #
  // ##########
  f0->SetParName(0,"Ampli-S");
  f0->SetParName(1,"Fraction-S");
  f0->SetParName(2,"#mu-S");
  f0->SetParName(3,"#sigma-1S");
  f0->SetParName(4,"#sigma-2S");

  // ##############
  // # Background #
  // ##############
  f0->SetParName(5,"Ampli1-B");
  f0->SetParName(6,"#mu1-B");
  f0->SetParName(7,"#sigma1-B");
  f0->SetParName(8,"Ampli2-B");
  f0->SetParName(9,"#mu2-B");
  f0->SetParName(10,"#sigma2-B");

  f0->SetParameter(0,350.0);
  f0->FixParameter(1,FRACJPSI);
  f0->SetParameter(2,5.27);
  f0->FixParameter(3,SIGMAS1JPSI);
  f0->FixParameter(4,SIGMAS2JPSI);

  f0->SetParameter(5,1300.0);
  f0->SetParameter(6,5.0);
  f0->SetParameter(7,0.09);
  f0->SetParameter(8,120.0);
  f0->SetParameter(9,5.4);
  f0->SetParameter(10,0.09);


  // ##################
  // # psi(2S) region #
  // ##################
  theTree->Draw("B0MassArb>>h1","B0MassArb > 5.27953 - 0.28 && B0MassArb < 5.27953 + 0.28 && genSignal == 0 && truthMatchSignal == 1 && mumuMass > 3.68609 - 3*mumuMassE && mumuMass < 3.68609 + 3*mumuMassE");
  h1->Scale(scaleMCdata);

  TF1* f1  = new TF1("f1", "[0]*([1]*TMath::Gaus(x,[2],[3]) + (1-[1])*TMath::Gaus(x,[2],[4])) + [5]*exp(-(x-[6])/[7])",B0Mass - B0MassIntervalLeft,B0Mass + B0MassIntervalRight);
  TF1* f1S = new TF1("f1S","[0]*([1]*TMath::Gaus(x,[2],[3]) + (1-[1])*TMath::Gaus(x,[2],[4]))",B0Mass - B0MassIntervalLeft,B0Mass + B0MassIntervalRight);
  f1S->SetLineColor(kBlue);
  f1S->SetLineStyle(kDashed);

  // ##########
  // # Signal #
  // ##########
  f1->SetParName(0,"Ampli-S");
  f1->SetParName(1,"Fraction-S");
  f1->SetParName(2,"#mu-S");
  f1->SetParName(3,"#sigma-1S");
  f1->SetParName(4,"#sigma-2S");

  // ##############
  // # Background #
  // ##############
  f1->SetParName(5,"Ampli-B");
  f1->SetParName(6,"#mu-B");
  f1->SetParName(7,"#tau-B");

  f1->SetParameter(0,170.0);
  f1->FixParameter(1,FRACPSIP);
  f1->SetParameter(2,5.27);
  f1->FixParameter(3,SIGMAS1PSIP);
  f1->FixParameter(4,SIGMAS2PSIP);

  f1->SetParameter(5,20.0);
  f1->SetParameter(6,5.27);
  f1->SetParameter(7,0.2);


  // ###############################
  // # J/psi AND psi(2S) rejection #
  // ###############################
  theTree->Draw("B0MassArb>>h2","B0MassArb > 5.27953 - 0.28 && B0MassArb < 5.27953 + 0.28 && genSignal == 0 && truthMatchSignal == 1 && mumuMass < 3.096916 - 5*mumuMassE || mumuMass > 3.68609 + 3*mumuMassE || (mumuMass > 3.096916 + 3*mumuMassE && mumuMass < 3.68609 - 3*mumuMassE)");
  h2->Scale(scaleMCdata);


  c0->cd();
  if (iFit == true) h0->Fit("f0");
  h0->Draw("e1p");
  double intErrJpsi = f0->IntegralError(B0Mass - B0MassIntervalLeft,B0Mass + B0MassIntervalRight)/((B0MassIntervalLeft+B0MassIntervalRight) / ((double)nBins));
  f0S->SetParameter(0,f0->GetParameter(0));
  f0S->SetParameter(1,f0->GetParameter(1));
  f0S->SetParameter(2,f0->GetParameter(2));
  f0S->SetParameter(3,f0->GetParameter(3));
  f0S->SetParameter(4,f0->GetParameter(4));
  double intJpsi0  = f0->Integral(B0Mass - B0MassIntervalLeft,B0Mass + B0MassIntervalRight)/((B0MassIntervalLeft+B0MassIntervalRight) / ((double)nBins));
  double intJpsi0S = f0S->Integral(B0Mass - B0MassIntervalLeft,B0Mass + B0MassIntervalRight)/((B0MassIntervalLeft+B0MassIntervalRight) / ((double)nBins));
  if (iFit == true)
    {
      f0S->Draw("same");
      leg0 = new TLegend(0.2, 0.79, 0.4, 0.89, "");
      leg0->AddEntry(h0,"Data");
      leg0->AddEntry(f0,"Total fit");
      leg0->AddEntry(f0S,"Signal");
      leg0->SetFillColor(0);
      leg0->SetBorderSize(0);
      leg0->Draw();
    }
  c0->Update();

  c1->cd();
  if (iFit == true) h1->Fit("f1");
  h1->Draw("e1p");
  double intErrPsiP = f1->IntegralError(B0Mass - B0MassIntervalLeft,B0Mass + B0MassIntervalRight)/((B0MassIntervalLeft+B0MassIntervalRight) / ((double)nBins));
  f1S->SetParameter(0,f1->GetParameter(0));
  f1S->SetParameter(1,f1->GetParameter(1));
  f1S->SetParameter(2,f1->GetParameter(2));
  f1S->SetParameter(3,f1->GetParameter(3));
  f1S->SetParameter(4,f1->GetParameter(4));
  f1S->SetParameter(5,f1->GetParameter(5));
  double intPsiP1  = f1->Integral(B0Mass - B0MassIntervalLeft,B0Mass + B0MassIntervalRight)/((B0MassIntervalLeft+B0MassIntervalRight) / ((double)nBins));
  double intPsiP1S = f1S->Integral(B0Mass - B0MassIntervalLeft,B0Mass + B0MassIntervalRight)/((B0MassIntervalLeft+B0MassIntervalRight) / ((double)nBins));
  if (iFit == true)
    {
      f1S->Draw("same");
      leg1 = new TLegend(0.2, 0.79, 0.4, 0.89, "");
      leg1->AddEntry(h1,"Data");
      leg1->AddEntry(f1,"Total fit");
      leg1->AddEntry(f1S,"Signal");
      leg1->SetFillColor(0);
      leg1->SetBorderSize(0);
      leg1->Draw();
    }
  c1->Update();

  c2->cd();
  h2->Draw("e1p");
  c2->Update();

  cout << "\nJ/psi total integral: "   << intJpsi0  << "+/-" << intErrJpsi << endl;
  cout << "J/psi signal integral: "    << intJpsi0S << "+/-" << intErrJpsi / intJpsi0 * intJpsi0S << endl;
  cout << "\npsi(2S) total integral: " << intPsiP1  << "+/-" << intErrPsiP << endl;
  cout << "psi(2S) signal integral: "  << intPsiP1S << "+/-" << intErrPsiP / intPsiP1 * intPsiP1S << endl;
}


// ########################################################################
// # Sub-program to reduce the number of branches of the candidate ntuple #
// ########################################################################
void ReduceTree (string fileNameIn, string fileNameOut)
{
  int nEntries;

  TFile* fileIn = TFile::Open(fileNameIn.c_str(),"READ");
  TTree* theTreeIn = (TTree*)fileIn->Get("B0KstMuMu/B0KstMuMuNTuple");

  nEntries = theTreeIn->GetEntries();
  cout << "\n@@@ Total number of events in the input tree: " << nEntries << " @@@" << endl;

  theTreeIn->SetBranchStatus("*",0);

  // ######
  // # B0 #
  // ######
  theTreeIn->SetBranchStatus("bMass",1);
  theTreeIn->SetBranchStatus("bBarMass",1);
  theTreeIn->SetBranchStatus("bPx",1);
  theTreeIn->SetBranchStatus("bPy",1);
  theTreeIn->SetBranchStatus("bPz",1);
  theTreeIn->SetBranchStatus("bLBS",1);
  theTreeIn->SetBranchStatus("bLBSE",1);
  theTreeIn->SetBranchStatus("bVtxCL",1);
  theTreeIn->SetBranchStatus("bCosAlphaBS",1);

  theTreeIn->SetBranchStatus("TrigTable",1);

  // #########
  // # Mu-Mu #
  // #########
  theTreeIn->SetBranchStatus("mumuVtxCL",1);
  theTreeIn->SetBranchStatus("mumuMass",1);
  theTreeIn->SetBranchStatus("mumuMassE",1);

  // #######
  // # K*0 #
  // #######
  theTreeIn->SetBranchStatus("kstMass",1);
  theTreeIn->SetBranchStatus("kstMassE",1);
  theTreeIn->SetBranchStatus("kstBarMass",1);
  theTreeIn->SetBranchStatus("kstBarMassE",1);
  theTreeIn->SetBranchStatus("kstPx",1);
  theTreeIn->SetBranchStatus("kstPy",1);
  theTreeIn->SetBranchStatus("kstPz",1);

  // #######
  // # mu- #
  // #######
  theTreeIn->SetBranchStatus("mumHighPurity",1);
  theTreeIn->SetBranchStatus("mumNormChi2",1);
  theTreeIn->SetBranchStatus("mumCat",1);
  theTreeIn->SetBranchStatus("mumPx",1);
  theTreeIn->SetBranchStatus("mumPy",1);
  theTreeIn->SetBranchStatus("mumPz",1);
  theTreeIn->SetBranchStatus("mumDCABS",1);
  theTreeIn->SetBranchStatus("mumDCABSE",1);
  theTreeIn->SetBranchStatus("mumdxyVtx",1);
  theTreeIn->SetBranchStatus("mumdzVtx",1);
  theTreeIn->SetBranchStatus("mumTrig",1);
  theTreeIn->SetBranchStatus("mumNPixLayers",1);
  theTreeIn->SetBranchStatus("mumNTrkLayers",1);

  // #######
  // # mu+ #
  // #######
  theTreeIn->SetBranchStatus("mupHighPurity",1);
  theTreeIn->SetBranchStatus("mupNormChi2",1);
  theTreeIn->SetBranchStatus("mupCat",1);
  theTreeIn->SetBranchStatus("mupPx",1);
  theTreeIn->SetBranchStatus("mupPy",1);
  theTreeIn->SetBranchStatus("mupPz",1);
  theTreeIn->SetBranchStatus("mupDCABS",1);
  theTreeIn->SetBranchStatus("mupDCABSE",1);
  theTreeIn->SetBranchStatus("mupdxyVtx",1);
  theTreeIn->SetBranchStatus("mupdzVtx",1);
  theTreeIn->SetBranchStatus("mupTrig",1);
  theTreeIn->SetBranchStatus("mupNPixLayers",1);
  theTreeIn->SetBranchStatus("mupNTrkLayers",1);

  // ##############
  // # K*0 track- #
  // ##############
  theTreeIn->SetBranchStatus("kstTrkmHighPurity",1);
  theTreeIn->SetBranchStatus("kstTrkmMuMatch",1);
  theTreeIn->SetBranchStatus("kstTrkmDCABS",1);
  theTreeIn->SetBranchStatus("kstTrkmDCABSE",1);
  theTreeIn->SetBranchStatus("kstTrkmPx",1);
  theTreeIn->SetBranchStatus("kstTrkmPy",1);
  theTreeIn->SetBranchStatus("kstTrkmPz",1);

  // ##############
  // # K*0 track+ #
  // ##############
  theTreeIn->SetBranchStatus("kstTrkpHighPurity",1);
  theTreeIn->SetBranchStatus("kstTrkpMuMatch",1);
  theTreeIn->SetBranchStatus("kstTrkpDCABS",1);
  theTreeIn->SetBranchStatus("kstTrkpDCABSE",1);
  theTreeIn->SetBranchStatus("kstTrkpPx",1);
  theTreeIn->SetBranchStatus("kstTrkpPy",1);
  theTreeIn->SetBranchStatus("kstTrkpPz",1);

  theTreeIn->SetBranchStatus("genSignal",1);
  theTreeIn->SetBranchStatus("truthMatchSignal",1);

  TFile* fileOut = TFile::Open(fileNameOut.c_str(),"RECREATE");
  fileOut->mkdir("B0KstMuMu");
  fileOut->cd("B0KstMuMu");
  TTree* theTreeOut = theTreeIn->CloneTree();

  nEntries = theTreeOut->GetEntries();
  cout << "@@@ Total number of events in the output tree: " << nEntries << " @@@" << endl;

  theTreeOut->Print();
  fileOut->Write();
  fileIn->Close();
  fileOut->Close();
}


// #############################################################################################################
// # Sub-program to sample the MC in order to obtain a subset that has the same pileup distribution as in data #
// #############################################################################################################
void SampleMCforPileup (string fileNameIn, string fileNameOut)
{
  TRandom3* myRandom = new TRandom3();

  TFile* fileIn = TFile::Open(fileNameIn.c_str(),"READ");
  TTree* theTreeIn = (TTree*)fileIn->Get("B0KstMuMu/B0KstMuMuNTuple");

  int nEntries = theTreeIn->GetEntries();
  cout << "\n@@@ Total number of events in the input tree: " << nEntries << " @@@" << endl;


  TFile* fileOut = TFile::Open(fileNameOut.c_str(),"RECREATE");
  fileOut->mkdir("B0KstMuMu");
  fileOut->cd("B0KstMuMu");
  TTree* theTreeOut = theTreeIn->CloneTree(0);


  vector<double>* vec1 = new vector<double>;
  vector<double>* vec2 = new vector<double>;

  TH1D* pileupW = new TH1D("pileupW","pileupW",50,0,50);
  theTreeIn->Draw("numInteractionsMC>>pileupW","evWeight*(bunchXingMC == 0)","goff");
  pileupW->Scale(1. / pileupW->Integral());

  TH1D* pileupNoW = new TH1D("pileupNoW","pileupNoW",50,0,50);
  theTreeIn->Draw("numInteractionsMC>>pileupNoW","bunchXingMC == 0","goff");
  pileupNoW->Scale(1. / pileupNoW->Integral());

  TH1D* pileupRatio = (TH1D*)pileupW->Clone("pileupRatio");
  pileupRatio->Divide(pileupNoW);
  pileupRatio->Scale(1. / pileupRatio->GetBinContent(pileupRatio->GetMaximumBin()));


  for (int entry = 0; entry < nEntries; entry++)
    {
      theTreeIn->SetBranchAddress("bunchXingMC", &vec1);
      theTreeIn->SetBranchAddress("numInteractionsMC", &vec2);
      theTreeIn->GetEntry(entry);

      for (unsigned int it = 0; it < (*vec1).size(); it++)
  	if (((*vec1)[it] == 0) && (myRandom->Uniform() < pileupRatio->GetBinContent(pileupRatio->FindBin((*vec2)[it])))) theTreeOut->Fill();

      vec1->clear();
      vec2->clear();
    }

  nEntries = theTreeOut->GetEntries();
  cout << "@@@ Total number of events in the output tree: " << nEntries << " @@@" << endl;


  theTreeOut->Print();
  theTreeOut->AutoSave();
  fileIn->Close();
  fileOut->Close();
}


// ###################################################
// # Sub-program to divide the ntuple in n sub-files #
// ###################################################
void DivideNTuple (string fileNameIn, string fileNameOut, unsigned int n)
{
  vector<TTree*> theTreeInCat;
  vector<TTree*> theTreeOutCat;
  stringstream myString;
  unsigned int entryStep;
  int nEntries;
  TFile* fileOut;
  TTree* theTreeOut;
  TList* listOfTrees = new TList();

  TFile* fileIn = TFile::Open(fileNameIn.c_str(),"READ");
  TTree* theTreeIn = (TTree*)fileIn->Get("B0KstMuMu/B0KstMuMuNTuple");


  // ###################################################
  // # Divide the input tree in the HLT sub-categories #
  // ###################################################
  cout << "\n@@@ Dividing the input tree in the HLT sub-categories @@@" << endl;
  TFile* fileTmp = TFile::Open("fileTmp.root","RECREATE");
  fileTmp->mkdir("B0KstMuMu");
  fileTmp->cd("B0KstMuMu");
  theTreeInCat.push_back(theTreeIn->CopyTree("TrigCat == 1"));
  theTreeInCat.push_back(theTreeIn->CopyTree("TrigCat == 2"));
  theTreeInCat.push_back(theTreeIn->CopyTree("TrigCat == 3"));
  theTreeInCat.push_back(theTreeIn->CopyTree("TrigCat == 4"));


  fileNameOut.replace(fileNameOut.find(".root"),5,"");
  for (unsigned int i = 0; i < n; i++)
    {
      myString.clear();
      myString.str("");
      myString << fileNameOut << "_" << i << ".root";
      cout << "\n@@@ Making file n #" << i << " --> " << myString.str().c_str() << " @@@" << endl;

      fileOut = TFile::Open(myString.str().c_str(),"RECREATE");
      fileOut->mkdir("B0KstMuMu");
      fileOut->cd("B0KstMuMu");

      for (unsigned int j = 0; j <theTreeInCat.size(); j++)
	{	  
	  // ########################
	  // # Clone the input tree #
	  // ########################
	  theTreeOutCat.push_back(theTreeInCat[j]->CloneTree(0));
	  entryStep = theTreeInCat[j]->GetEntries() / n;
	  cout << "@@@ Total number of events in the input tree category #" << j+1 << " : " << entryStep << " @@@" << endl;


	  // ####################
	  // # Copy the entries #
	  // ####################
	  for (unsigned int entry = i*entryStep; entry < ((i+1)*entryStep < theTreeInCat[j]->GetEntries() ? (i+1)*entryStep : theTreeInCat[j]->GetEntries()); entry++)
	    {
	      theTreeInCat[j]->GetEntry(entry);     
	      theTreeOutCat[j]->Fill();
	    }
	  nEntries = theTreeOutCat[j]->GetEntries();
	  cout << "@@@ Total number of events in the output tree in the category #" << j+1 << " : " << nEntries << " @@@" << endl;


	  listOfTrees->Add(theTreeOutCat[j]);
	}


      theTreeOut = TTree::MergeTrees(listOfTrees);
      nEntries = theTreeOut->GetEntries();
      cout << "@@@ Total number of events in the output tree: " << nEntries << " @@@" << endl;

      theTreeOut->Write();
      fileOut->Close();
      delete fileOut;
      listOfTrees->Clear();
      theTreeOutCat.clear();
    }


  fileIn->Close();
  fileTmp->Close();
}


// #####################################################
// # Sub-program to sample the single-candidate ntuple #
// #####################################################
void SampleNTuple (string fileNameIn, string fileNameOut, double fraction)
// ############################################
// # "fraction" refers to the retained events #
// ############################################
{
  vector<TTree*> theTreeInCat;
  vector<TTree*> theTreeOutCat;
  int nEntries;


  TFile* fileIn = TFile::Open(fileNameIn.c_str(),"READ");
  TTree* theTreeIn = (TTree*)fileIn->Get("B0KstMuMu/B0KstMuMuNTuple");

  nEntries = theTreeIn->GetEntries();
  cout << "\n@@@ Total number of events in the input tree: " << nEntries << " @@@" << endl;


  TFile* fileOut = TFile::Open(fileNameOut.c_str(),"RECREATE");
  fileOut->mkdir("B0KstMuMu");
  fileOut->cd("B0KstMuMu");


  // ###################################################
  // # Divide the input tree in the HLT sub-categories #
  // ###################################################
  cout << "Dividing the input tree in the HLT sub-categories" << endl;
  theTreeInCat.push_back(theTreeIn->CopyTree("TrigCat == 1"));
  theTreeInCat.push_back(theTreeIn->CopyTree("TrigCat == 2"));
  theTreeInCat.push_back(theTreeIn->CopyTree("TrigCat == 3"));
  theTreeInCat.push_back(theTreeIn->CopyTree("TrigCat == 4"));


  // #####################################################
  // # Select from each HLT category a sub-set of events #
  // #####################################################
  for (unsigned int i = 0; i < theTreeInCat.size(); i++)
    {
      cout << "Selecting from HLT category #" << i+1 << " a sub-set of events" << endl;
      theTreeOutCat.push_back(theTreeInCat[i]->CloneTree(int(rint((double)theTreeInCat[i]->GetEntries()*fraction))));
    }


  // #####################
  // # Sum each sub-tree #
  // #####################
  TList* listOfTrees = new TList();
  for (unsigned int i = 0; i < theTreeInCat.size(); i++)
    {
      cout << "Summing sub-treee #" << i+1 << endl;
      listOfTrees->Add(theTreeOutCat[i]);
    }
  TTree* theTreeOut = TTree::MergeTrees(listOfTrees);
  nEntries = theTreeOut->GetEntries();
  cout << "@@@ Total number of events in the output tree: " << nEntries << " @@@" << endl;
  

  theTreeOut->Write();
  fileOut->Close();
  fileIn->Close();
}


// ########################################
// # Code to compute MC filter efficiency #
// ########################################
void ComputeMCfilterEff (string fileName)
{
  // #################
  // # Read the tree #
  // #################
  TFile* _file0 = TFile::Open(fileName.c_str(),"READ");
  TTree* B0KstMuMuNTuple = (TTree*)_file0->Get("B0KstMuMu/B0KstMuMuNTuple");

  int nEntries = B0KstMuMuNTuple->GetEntries();
  cout << "@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;


  unsigned int val1      = 0;
  unsigned int val2      = 0;
  unsigned long evTried  = 0;
  unsigned long evPassed = 0;

  for (int entry = 0; entry < nEntries; entry++)
    {
      B0KstMuMuNTuple->SetBranchAddress("numEventsTried",  &val1);
      B0KstMuMuNTuple->SetBranchAddress("numEventsPassed", &val2);
      B0KstMuMuNTuple->GetEntry(entry);
      
      evTried  += val1;
      evPassed += val2;
    }


  cout << "Total number of events tried: "  << evTried  << endl;
  cout << "Total number of events passed: " << evPassed << endl;
  cout << "Monte Carlo filter efficiency: " << static_cast<double>(evPassed) / static_cast<double>(evTried);
  cout << " +/- " << sqrt(static_cast<double>(evPassed) * (1. - static_cast<double>(evPassed) / static_cast<double>(evTried))) / static_cast<double>(evTried) << endl;
}


// #######################################
// # Code to plot experiment-comparisons #
// #######################################
void DrawString (double Lumi)
{
  stringstream myString;

  myString.clear();
  myString.str("");
  myString << "CMS";
  TLatex* LumiTex1 = new TLatex(0.1,0.91,myString.str().c_str());
  LumiTex1->SetTextSize(0.05);
  LumiTex1->SetTextColor(kBlack);
  LumiTex1->SetNDC(true);
  LumiTex1->DrawLatex(0.1,0.91,myString.str().c_str());

  myString.clear();
  myString.str("");
  myString << "L = " << Lumi <<  " fb#lower[0.4]{^{#font[122]{\55}1}}";
  TLatex* LumiTex2 = new TLatex(0.43,0.91,myString.str().c_str());
  LumiTex2->SetTextSize(0.05);
  LumiTex2->SetTextColor(kBlack);
  LumiTex2->SetNDC(true);
  LumiTex2->DrawLatex(0.43,0.91,myString.str().c_str());

  // ##################
  // # Custom method: #
  // ##################
  double startNDCx = 0.826;
  double startNDCy = 0.935;
  TLine* line1 = new TLine(startNDCx-0.004, startNDCy, startNDCx, startNDCy);
  line1->SetBit(TLine::kLineNDC,true);
  line1->Draw();
  TLine* line2 = new TLine(startNDCx, startNDCy, startNDCx+0.005, startNDCy-0.03);
  line2->SetBit(TLine::kLineNDC,true);
  line2->Draw();
  TLine* line3 = new TLine(startNDCx+0.005, startNDCy-0.03, startNDCx+0.010, startNDCy+0.01);
  line3->SetBit(TLine::kLineNDC,true);
  line3->Draw();
  TLine* line4 = new TLine(startNDCx+0.010, startNDCy+0.01, startNDCx+0.032, startNDCy+0.01);
  line4->SetBit(TLine::kLineNDC,true);
  line4->Draw();
  // ###################
  // # Nominal method: #
  // ###################
  // myString.clear();
  // myString.str("");
  // myString << "#sqrt{  }";
  // TLatex* LumiTex3 = new TLatex(0.82,0.9,myString.str().c_str());
  // LumiTex3->SetTextSize(0.053);
  // LumiTex3->SetTextColor(kBlack);
  // LumiTex3->SetNDC(true);
  // LumiTex3->DrawLatex(0.82,0.9,myString.str().c_str());

  myString.clear();
  myString.str("");
  myString << "s = 8 TeV";
  TLatex* LumiTex4 = new TLatex(0.84,0.91,myString.str().c_str());
  LumiTex4->SetTextSize(0.05);
  LumiTex4->SetTextColor(kBlack);
  LumiTex4->SetNDC(true);
  LumiTex4->DrawLatex(0.84,0.91,myString.str().c_str());
}


TCutG* DrawExclusion (double Xlow, double Xhigh, double Ylow, double Yhigh, string cutName, unsigned int fillStyle, unsigned int color)
{
  TCutG* ExclusionZone = new TCutG(cutName.c_str(),5);
  ExclusionZone->SetVarX("");
  ExclusionZone->SetVarY("");
  ExclusionZone->SetPoint(0,Xlow,Ylow);
  ExclusionZone->SetPoint(1,Xhigh,Ylow);
  ExclusionZone->SetPoint(2,Xhigh,Yhigh);
  ExclusionZone->SetPoint(3,Xlow,Yhigh);
  ExclusionZone->SetPoint(4,Xlow,Ylow);
  ExclusionZone->SetFillColor(color);
  ExclusionZone->SetFillStyle(fillStyle);
  ExclusionZone->Draw("F");

  return ExclusionZone;
}


void printData (int nBins, TVectorD V1, TVectorD V2, TVectorD V3, TVectorD V4, TVectorD V5, TVectorD V6)
{
  for (int i = 0; i < nBins; i++)
    cout << "Read values: " << V1[i] << "\t" << V2[i] << "\t" << V3[i] << "\t" << V4[i] << "\t" << V5[i] << "\t" << V6[i] << endl;

}


void offsetData (int nBins, TVectorD* V1, TVectorD* V2, TVectorD* V3, double offset)
{
  for (int i = 0; i < nBins; i++)
    {
      (*V1)[i] = (*V1)[i] + offset;
      if ((*V2)[i] != 0.0) (*V2)[i] = (*V2)[i] + offset;
      if ((*V3)[i] != 0.0) (*V3)[i] = (*V3)[i] - offset;
    }
}


TGraphAsymmErrors* readData (TString fileName, int dataType, int color, int markerType, bool doFill, int fillStyle, bool noHbar, double offset)
// #################################
// # dataType = 0 --> FL           #
// # dataType = 1 --> AFB          #
// # dataType = 2 --> BF           #
// # noHbar: choose horizontal bar #
// #################################
{
  int nBins = 5; // Number of q2 bins

  TString exmin = "xs-xl";
  TString exmax = "xh-xs";
  TString eymin = "sqrt(yelsta*yelsta + yelsys*yelsys)";
  TString eymax = "sqrt(yehsta*yehsta + yehsys*yehsys)";
  TString valY;
  if ((dataType == 2) && (fileName.Contains("CMS") == false) && (fileName.Contains("LHCb") == false) && (fileName.Contains("Atlas") == false))
    {
      cout << "\n\nFile name: " << fileName << endl;
      valY = "ys/(xh-xl)";
      eymin = "sqrt(yelsta*yelsta + yelsys*yelsys)/(xh-xl)";
      eymax = "sqrt(yehsta*yehsta + yehsys*yehsys)/(xh-xl)";
    }
  else valY = "ys";
  

  // #################
  // # Read the tree #
  // #################
  TTree* treeData = new TTree();
  treeData->ReadFile(fileName,"xs:xl:xh:ys:yelsta:yehsta:yelsys:yehsys");
  int nEntrie = treeData->GetEntries();
  cout << "N. entries: " << nEntrie << "\tN. bins: " << nBins << endl;
  treeData->Print();


  treeData->Draw("xs:"+exmin+":"+exmax+":"+valY,"","goff",nBins,nBins*dataType);
  cout << "Selected rows: " << treeData->GetSelectedRows() << endl;
  TVectorD V1(nBins,treeData->GetV1());
  TVectorD V2(nBins,treeData->GetV2());
  TVectorD V3(nBins,treeData->GetV3());
  TVectorD V4(nBins,treeData->GetV4());
  treeData->Draw(eymin+":"+eymax,"","goff",nBins,nBins*dataType);
  TVectorD V5(nBins,treeData->GetV1());
  TVectorD V6(nBins,treeData->GetV2());

  
  // #######################
  // # Reset x-axis errors #
  // #######################
  for (int i = 0; i < V3.GetNoElements(); i++)
    {
      if (((noHbar == true) && (fileName.Contains("Theory") == false)) || ((V5[i] == 0.0) && (V6[i] == 0.0)))
	{
	  V2[i] = 0.0;
	  V3[i] = 0.0;
	}
    }


  cout << "Data before offset:" << endl;
  printData(nBins,V1,V2,V3,V4,V5,V6);
  offsetData(nBins,&V1,&V2,&V3,offset);
  cout << "\nData after offset:" << endl;
  printData(nBins,V1,V2,V3,V4,V5,V6);

  TGraphAsymmErrors* gra = new TGraphAsymmErrors(V1,V4,V2,V3,V5,V6);
  if      (dataType == 0) gra->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}}); F_{L}");
  else if (dataType == 1) gra->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}}); A_{FB}");
  else if (dataType == 2) gra->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}}); dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}7}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");
  if (doFill == true)
    {
      gra->SetFillColor(color);
      gra->SetFillStyle(fillStyle);
    }
  else
    {
      gra->SetMarkerColor(color);
      gra->SetMarkerSize(1.2);
      gra->SetLineColor(color);
      gra->SetLineWidth(1);
      gra->SetMarkerStyle(markerType);
    }


  delete treeData;
  return gra;
}


void showData (int dataType, double offset, bool noHbar)
// #########################
// # dataType == 0 --> FL  #
// # dataType == 1 --> AFB #
// # dataType == 2 --> BF  #
// #########################
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadTopMargin(0.11);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetTitleOffset(1.1,"x");
  gStyle->SetTitleOffset(0.95,"y");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetEndErrorSize(8);
  TGaxis::SetMaxDigits(3);


  stringstream myString;

  myString.clear(); myString.str("");
  myString << DIREXPTCOMP << "BinningFile.root";
  TFile* binningFile = TFile::Open(myString.str().c_str(),"READ");
  TList* myList = binningFile->GetListOfKeys();
  TCanvas* myCanv = (TCanvas*)((dynamic_cast<TKey*>(myList->At(0)))->ReadObj());
  myList = myCanv->GetListOfPrimitives();
  TH1D* h0 = (TH1D*)(myList->At(1));
  if (dataType == 0)
    {
      h0->GetYaxis()->SetRangeUser(0.0,1.0);
      h0->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}}); F_{L}");      
    }
  else if (dataType == 1)
    {
      h0->GetYaxis()->SetRangeUser(-1.0,1.0);
      h0->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}}); A_{FB}");
    }
  else if (dataType == 2)
    {
      h0->GetYaxis()->SetRangeUser(0.0,1.2);
      h0->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}}); dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}7}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");
    }
  h0->SetLineStyle(2);


  TCanvas* cData  = new TCanvas("cData","cData",10,10,700,500);
  vector<TGraphAsymmErrors*> dVar;

  myString.clear(); myString.str("");
  myString << DIREXPTCOMP << "CMS.data";
  dVar.push_back(readData(myString.str().c_str(),dataType,1,20,false,0,noHbar,0.0*offset));

  myString.clear(); myString.str("");
  myString << DIREXPTCOMP << "LHCb_1fb.data";
  dVar.push_back(readData(myString.str().c_str(),dataType,2,21,false,0,noHbar,0.0*offset));

  // myString.clear(); myString.str("");
  // myString << DIREXPTCOMP << "Atlas.data";
  // if (dataType != 2) dVar.push_back(readData(myString.str().c_str(),dataType,6,22,false,0,noHbar,-3.0*offset));

  myString.clear(); myString.str("");
  myString << DIREXPTCOMP << "BaBar.data";
  dVar.push_back(readData(myString.str().c_str(),dataType,4,23,false,0,noHbar,2.0*offset));

  myString.clear(); myString.str("");
  myString << DIREXPTCOMP << "Belle.data";
  dVar.push_back(readData(myString.str().c_str(),dataType,8,28,false,0,noHbar,-2.0*offset));

  myString.clear(); myString.str("");
  myString << DIREXPTCOMP << "CDF.data";
  dVar.push_back(readData(myString.str().c_str(),dataType,kGray+1,29,false,0,noHbar,-1.0*offset));

  myString.clear(); myString.str("");
  myString << DIREXPTCOMP << "Theory.data";
  dVar.push_back(readData(myString.str().c_str(),dataType,kBlue,20,true,3001,noHbar,0.0*offset));


  cData->cd();
  h0->Draw();
  dVar[dVar.size()-1]->Draw("same e2");
  for (unsigned int i = 0; i < dVar.size()-1; i++) dVar[dVar.size()-i-2]->Draw("same p");


  unsigned int it = 0;
  TLegend* leg = NULL;
  leg = new TLegend(0.12, 0.6, 0.27, 0.88, "");
  leg->AddEntry(dVar[dVar.size()-1],"<SM>","F");
  leg->AddEntry(dVar[it++],"CMS","lp");
  leg->AddEntry(dVar[it++],"LHCb","lp");
  // if (dataType != 2) leg->AddEntry(dVar[it++],"Atlas","lp");
  leg->AddEntry(dVar[it++],"BaBar","lp");
  leg->AddEntry(dVar[it++],"Belle","lp");
  leg->AddEntry(dVar[it++],"CDF","lp");

  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();


  TLine* line = new TLine(2.0,0.0,19.0,0.0);
  line->SetLineStyle(kDashed);
  line->Draw();

  
  DrawString(5.2);
  DrawExclusion(08.68,10.09,-1.2,1.2,"RejectJPsi1",3001,kGray);
  DrawExclusion(12.86,14.18,-1.2,1.2,"RejectPsiP1",3001,kGray);


  cData->Update();
 }

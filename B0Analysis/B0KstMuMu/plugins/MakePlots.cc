// ########################################################################
// # Program to plot final histograms for the B0 --> K*0 mu+ mu- analysis #
// ########################################################################
// # Author: Mauro Dinardo                                                #
// ########################################################################

#ifndef __CINT__
#include <TROOT.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TCutG.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TExec.h>
#include <TGraphBentErrors.h>
#endif

#include <stdlib.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Utils.h"
#include "B0KstMuMuSingleCandTreeContent.h"

using namespace std;


// ####################
// # Global constants #
// ####################
#define ParameterFILE        "../python/ParameterFile.txt"
#define ParameterFILE_MCGEN  "../results/ParameterFile_Sig_MCGEN.txt"
#define ParameterFILE_MCRECO "../results/ParameterFile_Sig_MCRECO_sample.txt"

#define SMFL     "../../SMprediction/FLErr.dat"
#define SMAFB    "../../SMprediction/AFBErr.dat"
#define SMBF     "../../SMprediction/dBFdq2.dat"
#define SMBINFL  "../../SMprediction/BinnedFL.dat"
#define SMBINAFB "../../SMprediction/BinnedAFB.dat"
#define SMBINBF  "../../SMprediction/BinneddBFdq2.dat"

#define SingleCand_MCkstJPsi  "Candidates/MonteCarlo/SingleCand/singleCand_B0ToJPsiKst_BrMC_NTuples_Merged_pileupW.root"
#define SingleCand_MCkstPsi2S "Candidates/MonteCarlo/SingleCand/singleCand_B0ToPsi2SKst_BrMC_NTuples_Merged_pileupW.root"
#define SingleCand_Data       "Candidates/Data/singleCand_B0ToKstMuMu_DataRRPRv4v5v6v1B_NTuples_Merged.root"

#define FitSysFILE "../../Efficiency/EffSystematicsData/FitSystematics_q2Bin.txt"
#define OutHisto "Results"
#define YvalueOutsideLimits 10.0 // Value given to bins with zero error in order not to show them
#define FORPAPER true // "true" = make special layout for publication


// ####################
// # Global variables #
// ####################
Utils* Utility;


// #######################
// # Function Definition #
// #######################
void DrawString (double Lumi);
void MakeComparisonDataMC (unsigned int plotType);
TCutG* DrawExclusion (double Xlow, double Xhigh, double Ylow, double Yhigh, string cutName, unsigned int fillStyle, unsigned int color);
TGraphAsymmErrors* ReadFromASCII (string fileName, unsigned int PlotType, vector<double>* q2Bins, vector<double>* vxs, vector<double>* vys, vector<double>* vxel, vector<double>* vxeh, vector<double>* vyel, vector<double>* vyeh);
void CheckPhysicsRegion ();
void MakePhysicsPlots (unsigned int PlotType);
void EvalMultyRun (unsigned int sysType, unsigned int q2BinIndx, string fileName, double NLLinterval, double NLLlessThan);
void PlotMuMu (string fileName, bool bkgSub);
void PlotKst (string fileName, bool bkgSub);
void PlotKK (string fileName, bool bkgSub, string RECOorGEN);
void PlotMuHadMass (string fileName);
void MakeupNLLandPvalPlots (string fileName, int specBin, string PlotType);
void MakePvaluePlot (string toyMCfileName, int specBin, string PlotType);
void getHfromToy (string fileNameIn, string fileNameOut, unsigned int intVal);


// ###########################
// # Function Implementation #
// ###########################
void DrawString (double Lumi)
{
  stringstream myString;

  myString.str("");
  myString << "CMS";
  TLatex* LumiTex1 = new TLatex(0.1,0.91,myString.str().c_str());
  LumiTex1->SetTextSize(0.05);
  LumiTex1->SetTextColor(kBlack);
  LumiTex1->SetNDC(true);
  LumiTex1->DrawLatex(0.1,0.91,myString.str().c_str());

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
  // myString.str("");
  // myString << "#sqrt{  }";
  // TLatex* LumiTex3 = new TLatex(0.82,0.9,myString.str().c_str());
  // LumiTex3->SetTextSize(0.053);
  // LumiTex3->SetTextColor(kBlack);
  // LumiTex3->SetNDC(true);
  // LumiTex3->DrawLatex(0.82,0.9,myString.str().c_str());

  myString.str("");
  myString << "s = 7 TeV";
  TLatex* LumiTex4 = new TLatex(0.84,0.91,myString.str().c_str());
  LumiTex4->SetTextSize(0.05);
  LumiTex4->SetTextColor(kBlack);
  LumiTex4->SetNDC(true);
  LumiTex4->DrawLatex(0.84,0.91,myString.str().c_str());
}


void MakeComparisonDataMC (unsigned int plotType)
// #############################################
// # plotType =  0 --> B0 pT                   #
// # plotType =  1 --> B0 eta                  #

// # plotType =  2 --> mu+ pT                  #
// # plotType =  3 --> mu- pT                  #
// # plotType =  4 --> mu+ eta                 #
// # plotType =  5 --> mu- eta                 #
// # plotType =  6 --> mu+ phi, eta range      #
// # plotType =  7 --> mu+ phi, eta range      #
// # plotType =  8 --> mu+ phi, eta range      #
// # plotType =  9 --> mu+ phi, eta range      #
// # plotType = 10 --> mu- phi, eta range      #
// # plotType = 11 --> mu- phi, eta range      #
// # plotType = 12 --> mu- phi, eta range      #
// # plotType = 13 --> mu- phi, eta range      #

// # plotType = 14 --> K*0 trk+ pT             #
// # plotType = 15 --> K*0 trk- pT             #
// # plotType = 16 --> K*0 trk+ eta            #
// # plotType = 17 --> K*0 trk+ eta            #
// # plotType = 18 --> K*0 trk+ phi, eta range #
// # plotType = 19 --> K*0 trk+ phi, eta range #
// # plotType = 20 --> K*0 trk+ phi, eta range #
// # plotType = 21 --> K*0 trk+ phi, eta range #
// # plotType = 22 --> K*0 trk- phi, eta range #
// # plotType = 23 --> K*0 trk- phi, eta range #
// # plotType = 24 --> K*0 trk- phi, eta range #
// # plotType = 25 --> K*0 trk- phi, eta range #

// # plotType = 26 --> cos(theta_K)            #
// # plotType = 27 --> cos(theta_l)            #
// #############################################
{
  stringstream myString;
  const unsigned int NHisto = 2;
  TLegend* leg;
  vector<TFile*> Vfiles;
  vector<TTree*> TreeMC;
  vector<TH1D*> h1DVec;
  vector<string> queryMC;
  string weightVar = "evWeight";
  string fileName;


  // ##################
  // # Read the trees #
  // ##################
  Vfiles.push_back(TFile::Open(SingleCand_MCkstJPsi,"READ"));
  TreeMC.push_back((TTree*)Vfiles.back()->Get("B0SingleCand/B0KstMuMuSingleCandNTuple"));

  Vfiles.push_back(TFile::Open(SingleCand_MCkstPsi2S,"READ"));
  TreeMC.push_back((TTree*)Vfiles.back()->Get("B0SingleCand/B0KstMuMuSingleCandNTuple"));


  Vfiles.push_back(TFile::Open(SingleCand_Data,"READ"));
  TTree* TreeData = (TTree*)Vfiles.back()->Get("B0SingleCand/B0KstMuMuSingleCandNTuple");


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


  double minX = 0.0;
  double maxX = 40.0;
  unsigned int nBinsX = 100;
  string Xtitle = "";

  double minY = 0.0;
  double maxY = 0.0;

  double signalSigma = sqrt(Utility->GetGenericParam("FRACMASSS") * Utility->GetGenericParam("SIGMAS1") * Utility->GetGenericParam("SIGMAS1") +
			    (1. - Utility->GetGenericParam("FRACMASSS")) * Utility->GetGenericParam("SIGMAS2") * Utility->GetGenericParam("SIGMAS2"));
  cout << "\n@@@ Signal sigma: " << signalSigma << " @@@" << endl;


  // ######
  // # B0 #
  // ######
  if (plotType == 0)
    {
      Xtitle = "B^{0} p_{T} (GeV)";
      maxX = 100.0;

      fileName = "B0pT.pdf";
    }
  else if (plotType == 1)
    {
      Xtitle = "B^{0} #eta";
      minX = -3.0;
      maxX = 3.0;

      fileName = "B0eta.pdf";
    }

  // #########
  // # Muons #
  // #########
  else if (plotType == 2)
    {
      Xtitle = "#mu^{#font[122]{+}} p_{T} (GeV)";
      maxX = 40.0;

      fileName = "MuppT.pdf";
    }
  else if (plotType == 3)
    {
      Xtitle = "#mu^{#font[122]{\55}} p_{T} (GeV)";
      maxX = 40.0;

      fileName = "MumpT.pdf";
    }
  else if (plotType == 4)
    {
      Xtitle = "#mu^{#font[122]{+}} #eta";
      minX = -2.4;
      maxX = 2.4;

      fileName = "Mupeta.pdf";
    }
  else if (plotType == 5)
    {
      Xtitle = "#mu^{#font[122]{\55}} #eta";
      minX = -2.4;
      maxX = 2.4;

      fileName = "Mumeta.pdf";
    }
  else if (plotType == 6)
    {
      Xtitle = "#mu^{#font[122]{+}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = -2.4;
      maxY = -1.2;

      nBinsX = 15;

      fileName = "Mupphi_eta0.pdf";
    }
  else if (plotType == 7)
    {
      Xtitle = "#mu^{#font[122]{+}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = -1.2;
      maxY = 0;

      nBinsX = 15;

      fileName = "Mupphi_eta1.pdf";
    }
  else if (plotType == 8)
    {
      Xtitle = "#mu^{#font[122]{+}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = 0.0;
      maxY = 1.2;

      nBinsX = 15;

      fileName = "Mupphi_eta2.pdf";
    }
  else if (plotType == 9)
    {
      Xtitle = "#mu^{#font[122]{+}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = 1.2;
      maxY = 2.4;

      nBinsX = 15;

      fileName = "Mupphi_eta3.pdf";
    }
  else if (plotType == 10)
    {
      Xtitle = "#mu^{#font[122]{\55}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = -2.4;
      maxY = -1.2;

      nBinsX = 15;

      fileName = "Mumphi_eta0.pdf";
    }
  else if (plotType == 11)
    {
      Xtitle = "#mu^{#font[122]{\55}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = -1.2;
      maxY = 0;

      nBinsX = 15;

      fileName = "Mumphi_eta1.pdf";
    }
  else if (plotType == 12)
    {
      Xtitle = "#mu^{#font[122]{\55}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = 0.0;
      maxY = 1.2;

      nBinsX = 15;

      fileName = "Mumphi_eta2.pdf";
    }
  else if (plotType == 13)
    {
      Xtitle = "#mu^{#font[122]{\55}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = 1.2;
      maxY = 2.4;

      nBinsX = 15;

      fileName = "Mumphi_eta3.pdf";
    }

  // ###########
  // # Hadrons #
  // ###########
  else if (plotType == 14)
    {
      Xtitle = "#font[122]{K}^{*0} trk#font[122]{+} p_{T} (GeV)";
      maxX = 20.0;

      fileName = "KstTrkppT.pdf";
    }
  else if (plotType == 15)
    {
      Xtitle = "#font[122]{K}^{*0} trk{#font[122]{\55}} p_{T} (GeV)";
      maxX = 20.0;

      fileName = "KstTrkmpT.pdf";
    }
  else if (plotType == 16)
    {
      Xtitle = "#font[122]{K}^{*0} trk#font[122]{+} #eta";
      minX = -3.0;
      maxX = 3.0;

      fileName = "KstTrkpeta.pdf";
    }
  else if (plotType == 17)
    {
      Xtitle = "#font[122]{K}^{*0} trk{#font[122]{\55}} #eta";
      minX = -3.0;
      maxX = 3.0;

      fileName = "KstTrkmeta.pdf";
    }
  else if (plotType == 18)
    {
      Xtitle = "#font[122]{K}^{*0} trk#font[122]{+} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = -3.0;
      maxY = -1.2;

      nBinsX = 15;

      fileName = "KstTrkpphi_eta0.pdf";
    }
  else if (plotType == 19)
    {
      Xtitle = "#font[122]{K}^{*0} trk#font[122]{+} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = -1.2;
      maxY = 0.0;

      nBinsX = 15;

      fileName = "KstTrkpphi_eta1.pdf";
    }
  else if (plotType == 20)
    {
      Xtitle = "#font[122]{K}^{*0} trk#font[122]{+} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = 0.0;
      maxY = 1.2;

      nBinsX = 15;

      fileName = "KstTrkpphi_eta2.pdf";
    }
  else if (plotType == 21)
    {
      Xtitle = "#font[122]{K}^{*0} trk#font[122]{+} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = 1.2;
      maxY = 3.0;
 
      nBinsX = 15;

      fileName = "KstTrkpphi_eta3.pdf";
   }
  else if (plotType == 22)
    {
      Xtitle = "#font[122]{K}^{*0} trk{#font[122]{\55}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = -3.0;
      maxY = -1.2;
 
      nBinsX = 15;

      fileName = "KstTrkmphi_eta0.pdf";
   }
  else if (plotType == 23)
    {
      Xtitle = "#font[122]{K}^{*0} trk{#font[122]{\55}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = -1.2;
      maxY = 0.0;
 
      nBinsX = 15;

      fileName = "KstTrkmphi_eta1.pdf";
   }
  else if (plotType == 24)
    {
      Xtitle = "#font[122]{K}^{*0} trk{#font[122]{\55}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = 0.0;
      maxY = 1.2;

      nBinsX = 15;

      fileName = "KstTrkmphi_eta2.pdf";
    }
  else if (plotType == 25)
    {
      Xtitle = "#font[122]{K}^{*0} trk{#font[122]{\55}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = 1.2;
      maxY = 3.0;
 
      nBinsX = 15;

      fileName = "KstTrkmphi_eta3.pdf";
   }

  // ##########
  // # Angles #
  // ##########
  else if (plotType == 26)
    {
      Xtitle = "cos(#theta_{#font[122]{K}})";
      minX = -1.0;
      maxX = 1.0;

      nBinsX = 20;

      fileName = "CosThetaK_dataMC.pdf";
   }
  else if (plotType == 27)
    {
      Xtitle = "cos(#theta_{l})";
      minX = -1.0;
      maxX = 1.0;

      nBinsX = 20;

      fileName = "CosThetaL_dataMC.pdf";
   }


  // #################
  // # 1D histograms #
  // #################
  for (unsigned int i = 0; i < NHisto; i++)
    {
      myString.str("");
      myString << "h1D" << i;
      h1DVec.push_back(new TH1D(myString.str().c_str(),myString.str().c_str(),nBinsX,minX,maxX));
      h1DVec.back()->Sumw2();
      h1DVec.back()->SetXTitle(Xtitle.c_str());
      h1DVec.back()->SetYTitle("Norm. entries");
    }

  TH1D* hDsig1D = new TH1D("hDsig1D","hDsig1D",nBinsX,minX,maxX);
  hDsig1D->Sumw2();
  hDsig1D->SetXTitle(Xtitle.c_str());
  hDsig1D->SetYTitle("Norm. entries");
  hDsig1D->SetMarkerStyle(20);

  TH1D* hDbkg1D = new TH1D("hDbkg1D","hDbkg1D",nBinsX,minX,maxX);
  hDbkg1D->Sumw2();
  hDbkg1D->SetXTitle(Xtitle.c_str());
  hDbkg1D->SetYTitle("Norm. entries");
  hDbkg1D->SetMarkerStyle(20);


  // #################
  // # Query NTuples #
  // #################
  string sigMassQuery = "";
  string bkgMassQuery = "";
  string query        = "";
  string tmpstring    = "";
  string selection    = "";
  string aVar         = "";
  string bVar         = "";
  
  myString.clear();
  myString.str("");
  myString << "((abs(B0MassArb - " << Utility->B0Mass << ") < " << Utility->GetGenericParam("NSigmaB0S")*signalSigma << ") && ";
  myString << "((mumuMass > " << Utility->JPsiMass << "-" << Utility->GetGenericParam("NSigmaPsiBig") << "*mumuMassE && mumuMass < " << Utility->JPsiMass << "+" << Utility->GetGenericParam("NSigmaPsiSmall") << "*mumuMassE) || ";
  myString << "(abs(mumuMass - " << Utility->PsiPrimeMass << ") < " << Utility->GetGenericParam("NSigmaPsiSmall") << "*mumuMassE)))";
  sigMassQuery = myString.str();

  myString.clear();
  myString.str("");
  myString << "(((B0MassArb > " << Utility->B0Mass + Utility->GetGenericParam("NSigmaB0B")*signalSigma << " && B0MassArb < "
	   << Utility->B0Mass + (Utility->GetGenericParam("NSigmaB0B") + Utility->GetGenericParam("NSigmaB0S"))*signalSigma << ") || ";
  myString << "(B0MassArb < " << Utility->B0Mass - Utility->GetGenericParam("NSigmaB0B")*signalSigma << " && B0MassArb > "
	   << Utility->B0Mass - (Utility->GetGenericParam("NSigmaB0B") + Utility->GetGenericParam("NSigmaB0S"))*signalSigma << ")) && ";
  myString << "((mumuMass > " << Utility->JPsiMass << "-" << Utility->GetGenericParam("NSigmaPsiBig") << "*mumuMassE && mumuMass < "
	   << Utility->JPsiMass << "+" << Utility->GetGenericParam("NSigmaPsiSmall") << "*mumuMassE) || ";
  myString << "(abs(mumuMass - " << Utility->PsiPrimeMass << ") < " << Utility->GetGenericParam("NSigmaPsiSmall") << "*mumuMassE)))";
  bkgMassQuery = myString.str();
  if (plotType == 0)
    query = "B0pT"; // B0 pT
  else if (plotType == 1)
    query = "B0Eta"; // B0 eta

  else if (plotType == 2)
    query = "sqrt(mupPx*mupPx+mupPy*mupPy)"; // mu+ pT
  else if (plotType == 3)
    query = "sqrt(mumPx*mumPx+mumPy*mumPy)"; // mu- pT
  else if (plotType == 4)
    query = "0.5*log((sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) + mupPz) / (sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) - mupPz))"; // mu+ eta
  else if (plotType == 5)
    query = "0.5*log((sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) + mumPz) / (sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) - mumPz))"; // mu- eta
  else if (plotType == 6)
    {
      query = "atan(mupPy / mupPx)"; // mu+ phi in eta range
      tmpstring = "0.5*log((sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) + mupPz) / (sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) - mupPz))";
      myString.clear();
      myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // mu+ eta range
      aVar = "mupPx";
      bVar = "mupPy";
    }
  else if (plotType == 7)
    {
      query = "atan(mupPy / mupPx)"; // mu+ phi in eta range
      tmpstring = "0.5*log((sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) + mupPz) / (sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) - mupPz))";
      myString.clear();
      myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // mu+ eta range
      aVar = "mupPx";
      bVar = "mupPy";
    }
  else if (plotType == 8)
    {
      query = "atan(mupPy / mupPx)"; // mu+ phi in eta range
      tmpstring = "0.5*log((sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) + mupPz) / (sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) - mupPz))";
      myString.clear();
      myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // mu+ eta range
      aVar = "mupPx";
      bVar = "mupPy";
    }
  else if (plotType == 9)
    {
      query = "atan(mupPy / mupPx)"; // mu+ phi in eta range
      tmpstring = "0.5*log((sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) + mupPz) / (sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) - mupPz))";
      myString.clear();
      myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // mu+ eta range
      aVar = "mupPx";
      bVar = "mupPy";
    }
  else if (plotType == 10)
    {
      query = "atan(mumPy / mumPx)"; // mu- phi in eta range
      tmpstring = "0.5*log((sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) + mumPz) / (sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) - mumPz))";
      myString.clear();
      myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // mu- eta range
      aVar = "mumPx";
      bVar = "mumPy";
    }
  else if (plotType == 11)
    {
      query = "atan(mumPy / mumPx)"; // mu- phi in eta range
      tmpstring = "0.5*log((sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) + mumPz) / (sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) - mumPz))";
      myString.clear();
      myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // mu- eta range
      aVar = "mumPx";
      bVar = "mumPy";
    }
  else if (plotType == 12)
    {
      query = "atan(mumPy / mumPx)"; // mu- phi in eta range
      tmpstring = "0.5*log((sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) + mumPz) / (sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) - mumPz))";
      myString.clear();
      myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // mu- eta range
      aVar = "mumPx";
      bVar = "mumPy";
    }
  else if (plotType == 13)
    {
      query = "atan(mumPy / mumPx)"; // mu- phi in eta range
      tmpstring = "0.5*log((sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) + mumPz) / (sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) - mumPz))";
      myString.clear();
      myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // mu- eta range
      aVar = "mumPx";
      bVar = "mumPy";
    }
  
  else if (plotType == 14)
    query = "sqrt(kstTrkpPx*kstTrkpPx+kstTrkpPy*kstTrkpPy)"; // K*0 trk+ pT
  else if (plotType == 15)
    query = "sqrt(kstTrkmPx*kstTrkmPx+kstTrkmPy*kstTrkmPy)"; // K*0 trk- pT
  else if (plotType == 16)
    query = "0.5*log((sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) + kstTrkpPz) / (sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) - kstTrkpPz))"; // K*0 trk+ eta
  else if (plotType == 17)
    query = "0.5*log((sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) + kstTrkmPz) / (sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) - kstTrkmPz))"; // K*0 trk- eta
  else if (plotType == 18)
    {
      query = "atan(kstTrkpPy / kstTrkpPx)"; // K*0 trk+ phi in eta range
      tmpstring = "0.5*log((sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) + kstTrkpPz) / (sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) - kstTrkpPz))";
      myString.clear();
      myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // K*0 trk+ eta range
      aVar = "kstTrkpPx";
      bVar = "kstTrkpPy";
    }
  else if (plotType == 19)
    {
      query = "atan(kstTrkpPy / kstTrkpPx)"; // K*0 trk+ phi in eta range
      tmpstring = "0.5*log((sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) + kstTrkpPz) / (sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) - kstTrkpPz))";
      myString.clear();
      myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // K*0 trk+ eta range
      aVar = "kstTrkpPx";
      bVar = "kstTrkpPy";
    }
  else if (plotType == 20)
    {
      query = "atan(kstTrkpPy / kstTrkpPx)"; // K*0 trk+ phi in eta range
      tmpstring = "0.5*log((sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) + kstTrkpPz) / (sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) - kstTrkpPz))";
      myString.clear();
      myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // K*0 trk+ eta range
      aVar = "kstTrkpPx";
      bVar = "kstTrkpPy";
    }
  else if (plotType == 21)
    {
      query = "atan(kstTrkpPy / kstTrkpPx)"; // K*0 trk+ phi in eta range
      tmpstring = "0.5*log((sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) + kstTrkpPz) / (sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) - kstTrkpPz))";
      myString.clear();
      myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // K*0 trk+ eta range
      aVar = "kstTrkpPx";
      bVar = "kstTrkpPy";
    }
  else if (plotType == 22)
    {
      query = "atan(kstTrkmPy / kstTrkmPx)"; // K*0 trk- phi in eta range
      tmpstring = "0.5*log((sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) + kstTrkmPz) / (sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) - kstTrkmPz))";
      myString.clear();
      myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // K*0 trk- eta range
      aVar = "kstTrkmPx";
      bVar = "kstTrkmPy";
    }
  else if (plotType == 23)
    {
      query = "atan(kstTrkmPy / kstTrkmPx)"; // K*0 trk- phi in eta range
      tmpstring = "0.5*log((sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) + kstTrkmPz) / (sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) - kstTrkmPz))";
      myString.clear();
      myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // K*0 trk- eta range
      aVar = "kstTrkmPx";
      bVar = "kstTrkmPy";
    }
  else if (plotType == 24)
    {
      query = "atan(kstTrkmPy / kstTrkmPx)"; // K*0 trk- phi in eta range
      tmpstring = "0.5*log((sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) + kstTrkmPz) / (sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) - kstTrkmPz))";
      myString.clear();
      myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // K*0 trk- eta range
      aVar = "kstTrkmPx";
      bVar = "kstTrkmPy";
    }
  else if (plotType == 25)
    {
      query = "atan(kstTrkmPy / kstTrkmPx)"; // K*0 trk- phi in eta range
      tmpstring = "0.5*log((sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) + kstTrkmPz) / (sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) - kstTrkmPz))";
      myString.clear();
      myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // K*0 trk- eta range
      aVar = "kstTrkmPx";
      bVar = "kstTrkmPy";
    }
  
  else if (plotType == 26) query = "CosThetaKArb"; // cos(theta_K)
  else if (plotType == 27) query = "CosThetaMuArb"; // cos(theta_l)


  TCanvas* c0 = new TCanvas("c0","c0",10,10,700,500);


  // #######################
  // # Query the MC NTuple #
  // #######################
  TH1D* hM1D;
  cout << "\n\n@@@ Query to MC @@@" << endl;
  for (unsigned int i = 0; i < NHisto; i++)
    {
      if (((plotType >= 6) && (plotType <= 13)) || ((plotType >= 18) && (plotType <= 25)))
	{
	  TH1D* hTmp = new TH1D("hTmp","hTmp",nBinsX,minX,maxX);

	  tmpstring = weightVar + "*(" + selection + " && " + aVar + " > 0 && " + sigMassQuery + ")";
	  myString.clear();
	  myString.str("");
	  myString << query + ">>hTmp";
	  cout << "\nPlot: " << myString.str().c_str() << endl;
	  cout << "Selection: " << tmpstring.c_str() << endl;
	  TreeMC[i]->Draw(myString.str().c_str(),tmpstring.c_str());
	  h1DVec[i]->Add(hTmp);

	  tmpstring = weightVar + "*(" + selection + " && " + aVar + " < 0 && " + bVar + " < 0 && " + sigMassQuery + ")";
	  myString.clear();
	  myString.str("");
	  myString << query + "-TMath::Pi()>>hTmp";
	  cout << "\nPlot: " << myString.str().c_str() << endl;
	  cout << "Selection: " << tmpstring.c_str() << endl;
	  TreeMC[i]->Draw(myString.str().c_str(),tmpstring.c_str());
	  h1DVec[i]->Add(hTmp);

	  tmpstring = weightVar + "*(" + selection + " && " + aVar + " < 0 && " + bVar + " > 0 && " + sigMassQuery + ")";
	  myString.clear();
	  myString.str("");
	  myString << query + "+TMath::Pi()>>hTmp";
	  cout << "\nPlot: " << myString.str().c_str() << endl;
	  cout << "Selection: " << tmpstring.c_str() << endl;
	  TreeMC[i]->Draw(myString.str().c_str(),tmpstring.c_str());
	  h1DVec[i]->Add(hTmp);

	  delete hTmp;
	}
      else
	{
	  tmpstring = weightVar + "*(" + sigMassQuery + ")";
	  myString.clear();
	  myString.str("");
	  myString << "h1D" << i;
	  queryMC.push_back(query + ">>" + myString.str().c_str());
	  cout << "\nPlot: " << queryMC.back().c_str() << endl;
	  cout << "Selection: " << sigMassQuery.c_str() << endl;
	  TreeMC[i]->Draw(queryMC.back().c_str(),tmpstring.c_str());
	}
    }
  
  
  // #########################
  // # Query the Data NTuple #
  // #########################
  string queryData;
  cout << "\n\n@@@ Query to data @@@" << endl;
  if (((plotType >= 6) && (plotType <= 13)) || ((plotType >= 18) && (plotType <= 25)))
    {
      TH1D* hTmp = new TH1D("hTmp","hTmp",nBinsX,minX,maxX);
      

      tmpstring = selection + " && " + aVar + " > 0 && " + sigMassQuery;
      myString.clear();
      myString.str("");
      myString << query + ">>hTmp";
      cout << "\nPlot: " << myString.str().c_str() << endl;
      cout << "Selection: " << tmpstring.c_str() << endl;
      TreeData->Draw(myString.str().c_str(),tmpstring.c_str());
      hDsig1D->Add(hTmp);
      
      tmpstring = selection + " && " + aVar + " < 0 && " + bVar + " < 0 && " + sigMassQuery;
      myString.clear();
      myString.str("");
      myString << query + "-TMath::Pi()>>hTmp";
      cout << "\nPlot: " << myString.str().c_str() << endl;
      cout << "Selection: " << tmpstring.c_str() << endl;
      TreeData->Draw(myString.str().c_str(),tmpstring.c_str());
      hDsig1D->Add(hTmp);
      
      tmpstring = selection + " && " + aVar + " < 0 && " + bVar + " > 0 && " + sigMassQuery;
      myString.clear();
      myString.str("");
      myString << query + "+TMath::Pi()>>hTmp";
      cout << "\nPlot: " << myString.str().c_str() << endl;
      cout << "Selection: " << tmpstring.c_str() << endl;
      TreeData->Draw(myString.str().c_str(),tmpstring.c_str());
      hDsig1D->Add(hTmp);


      tmpstring = selection + " && " + aVar + " > 0 && " + bkgMassQuery;
      myString.clear();
      myString.str("");
      myString << query + ">>hTmp";
      cout << "\nPlot: " << myString.str().c_str() << endl;
      cout << "Selection: " << tmpstring.c_str() << endl;
      TreeData->Draw(myString.str().c_str(),tmpstring.c_str());
      hDbkg1D->Add(hTmp);
      
      tmpstring = selection + " && " + aVar + " < 0 && " + bVar + " < 0 && " + bkgMassQuery;
      myString.clear();
      myString.str("");
      myString << query + "-TMath::Pi()>>hTmp";
      cout << "\nPlot: " << myString.str().c_str() << endl;
      cout << "Selection: " << tmpstring.c_str() << endl;
      TreeData->Draw(myString.str().c_str(),tmpstring.c_str());
      hDbkg1D->Add(hTmp);
      
      tmpstring = selection + " && " + aVar + " < 0 && " + bVar + " > 0 && " + bkgMassQuery;
      myString.clear();
      myString.str("");
      myString << query + "+TMath::Pi()>>hTmp";
      cout << "\nPlot: " << myString.str().c_str() << endl;
      cout << "Selection: " << tmpstring.c_str() << endl;
      TreeData->Draw(myString.str().c_str(),tmpstring.c_str());
      hDbkg1D->Add(hTmp);


      hDsig1D->Add(hDbkg1D, -1.0);
      delete hTmp;
    }
  else
    {
      queryData = query + ">>hDsig1D";
      cout << "\nPlot: " << queryData.c_str() << endl;
      cout << "Selection: " << sigMassQuery.c_str() << endl;
      TreeData->Draw(queryData.c_str(),sigMassQuery.c_str());


      queryData = query + ">>hDbkg1D";
      cout << "\nPlot: " << queryData.c_str() << endl;
      cout << "Selection: " << bkgMassQuery.c_str() << endl;
      TreeData->Draw(queryData.c_str(),bkgMassQuery.c_str());


      hDsig1D->Add(hDbkg1D, -1.0);
    }


  // # B0 --> K* J/psi #
  h1DVec[0]->Scale(1./h1DVec[0]->Integral() * Utility->JPsiKpiBF);
  // # B0 --> K* psi(2S) #
  h1DVec[1]->Scale(1./h1DVec[1]->Integral() * Utility->PsiPKpiBF);


  hM1D = (TH1D*)h1DVec[0]->Clone("hM1D");
  hM1D->SetXTitle(Xtitle.c_str());
  hM1D->SetYTitle("Norm. entries");
  hM1D->SetLineColor(kBlack);
  hM1D->SetFillColor(kAzure+6);

  for (unsigned int i = 1; i < NHisto; i++) hM1D->Add(h1DVec[i]);
      
  hM1D->Scale(1./hM1D->Integral());
  hDsig1D->Scale(1./hDsig1D->Integral());

  c0->cd();
  hM1D->Draw("e3");
  hDsig1D->Draw("e1p sames");
  hM1D->GetYaxis()->SetRangeUser(0.0,(hM1D->GetBinContent(hM1D->GetMaximumBin()) > hDsig1D->GetBinContent(hDsig1D->GetMaximumBin()) ?
  				      hM1D->GetBinContent(hM1D->GetMaximumBin()) : hDsig1D->GetBinContent(hDsig1D->GetMaximumBin()))*1.1);
  hDsig1D->GetYaxis()->SetRangeUser(0.0,(hM1D->GetBinContent(hM1D->GetMaximumBin()) > hDsig1D->GetBinContent(hDsig1D->GetMaximumBin()) ?
  					 hM1D->GetBinContent(hM1D->GetMaximumBin()) : hDsig1D->GetBinContent(hDsig1D->GetMaximumBin()))*1.1);

  if (((plotType >= 6) && (plotType <= 13)) ||
      ((plotType >= 18) && (plotType <= 25)) ||
      (plotType == 26) ||
      (plotType == 27)) leg = new TLegend(0.15, 0.15, 0.34, 0.25, "");
  else if ((plotType == 1) ||
	   (plotType == 4) ||
	   (plotType == 5) ||
	   (plotType == 16) ||
	   (plotType == 17)) leg = new TLegend(0.15, 0.79, 0.34, 0.89, "");
  else leg = new TLegend(0.8, 0.15, 0.97, 0.25, "");
  leg->AddEntry(hM1D,"MC");
  leg->AddEntry(hDsig1D,"Data");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  c0->Update();

  TPaveStats* stM = (TPaveStats*)hM1D->FindObject("stats");
  stM->SetFillColor(kAzure+6);
  TPaveStats* stD = (TPaveStats*)hDsig1D->FindObject("stats");
  if ((plotType == 0) ||
      (plotType == 2) ||
      (plotType == 3) ||
      (plotType == 14) ||
      (plotType == 15) ||
      (plotType == 26))
    {
      stD->SetX1NDC(0.55);
      stD->SetX2NDC(0.75);
    }
  else if ((plotType == 1) ||
	   (plotType == 4) ||
	   (plotType == 5) ||
	   (plotType == 16) ||
	   (plotType == 17))
    {
      stD->SetY1NDC(0.4);
      stD->SetY2NDC(0.7);
    }
  else if (((plotType >= 6) && (plotType <= 13)) ||
	   ((plotType >= 18) && (plotType <= 25)) ||
	   (plotType == 27))
    {
      stD->SetX1NDC(0.4);
      stD->SetX2NDC(0.6);
      stD->SetY1NDC(0.15);
      stD->SetY2NDC(0.45);

      stM->SetX1NDC(0.65);
      stM->SetX2NDC(0.85);
      stM->SetY1NDC(0.15);
      stM->SetY2NDC(0.45);
    }


  // ######################################
  // # Estimate the percentage difference #
  // ######################################
  TH1D* hdiff = (TH1D*)hM1D->Clone("hdiff");
  hdiff->Add(hDsig1D,-1);
  for (int i = 0; i < hdiff->GetNbinsX(); i++) hdiff->SetBinContent(i+1,fabs(hdiff->GetBinContent(i+1)));
  cout << "\n\n@@@ Percentage difference: " << hdiff->Integral() / 2. << " @@@" << endl;


  // ###################
  // # data / MC ratio #
  // ###################
  TCanvas* c1 = new TCanvas("c1","c1",20,20,700,500);

  c1->cd();
  TH1D* hratio = (TH1D*)hDsig1D->Clone("hratio");
  hratio->GetYaxis()->SetRangeUser(0.6,1.4);
  hratio->SetStats(false);
  hratio->Divide(hM1D);
  hratio->Draw("pe1");
  c1->Update();


  c0->Print(fileName.c_str());
  c1->Print(fileName.replace(fileName.find(".pdf"),4,"_ratio.pdf").c_str());


  TreeMC.clear();
  h1DVec.clear();
  queryMC.clear();
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


TGraphAsymmErrors* ReadFromASCII (string fileName, unsigned int PlotType, vector<double>* q2Bins, vector<double>* vxs, vector<double>* vys, vector<double>* vxel, vector<double>* vxeh, vector<double>* vyel, vector<double>* vyeh)
{
  ifstream inputFile;

  // ############################################
  // # Variables to read values from ASCII file #
  // ############################################
  double xs  = 0.0;
  double ys  = 0.0;

  double xel = 0.0;
  double xeh = 0.0;

  double yel = 0.0;
  double yeh = 0.0;

  vxs->clear();
  vys->clear();
  vxel->clear();
  vxeh->clear();
  vyel->clear();
  vyeh->clear();

  
  inputFile.open(fileName.c_str(), ifstream::in);
  if (inputFile.good() == false)
    {
      cout << "[MakePlots::ReadFromASCII]\tError opening file : " << fileName.c_str() << endl;
      exit (1);
    }
  inputFile >> xs >> ys >> xel >> xeh >> yel >> yeh;
  while (inputFile)
    {
      if ((Utility->ValIsInPsi(q2Bins,xs) == false) && (xs >= q2Bins->operator[](0)) && (xs <= q2Bins->operator[](q2Bins->size()-1)) && (Utility->ValIsBetweenJPsiAndPsiP(q2Bins,xs) == false))
	{
	  if ((PlotType == 0) || (PlotType == 10)) // Fl
	    {
	      vxs->push_back(xs);
	      vys->push_back(ys);
	      vxel->push_back(xel);
	      vxeh->push_back(xeh);
	      vyel->push_back(yel);
	      vyeh->push_back(yeh);
	    }
	  else if ((PlotType == 1) || (PlotType == 11)) // Afb
	    {
	      vxs->push_back(xs);
	      vys->push_back(-ys);
	      vxel->push_back(xel);
	      vxeh->push_back(xeh);
	      vyel->push_back(yeh);
	      vyeh->push_back(yel);
	    }
	  else if ((PlotType == 2) || (PlotType == 12)) // Branching fraction
	    {
	      vxs->push_back(xs);
	      vys->push_back(ys/1e-7);
	      vxel->push_back(xel);
	      vxeh->push_back(xeh);
	      vyel->push_back(yel/1e-7);
	      vyeh->push_back(yeh/1e-7);
	    }
	  cout << "xs: " << xs << "\tys: " << ys << "\txeh: " << xeh << "\txel: " << xel << "\tyel: " << yel << "\tyeh: " << yeh << endl;
	}
      else
	{
	  vxs->push_back(xs);
	  vys->push_back(YvalueOutsideLimits);
	  vxel->push_back(0.0);
	  vxeh->push_back(0.0);
	  vyel->push_back(0.0);
	  vyeh->push_back(0.0);
	}
      inputFile >> xs >> ys >> xel >> xeh >> yel >> yeh;
    }
  inputFile.close();

  TGraphAsymmErrors* ge = new TGraphAsymmErrors(vxs->size(), &(*vxs)[0], &(*vys)[0], &(*vxel)[0], &(*vxeh)[0], &(*vyel)[0], &(*vyeh)[0]);

  return ge;
}


void CheckPhysicsRegion ()
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
  gStyle->SetTitleOffset(1.25,"y"); 
  TGaxis::SetMaxDigits(3);

  
  TCanvas* canv0 = new TCanvas("canv0","canv0",10,10,1200,600);
  TGraphAsymmErrors* ge;
  TGraphAsymmErrors* geTMP;
  TH1D* histo = new TH1D("histo","histo",100,-1.0,1.0);
  histo->GetYaxis()->SetRangeUser(0.0,1.0);
  histo->SetXTitle("A_{FB}");
  histo->SetYTitle("F_{L}");
  TLine* line1;
  TLine* line2;

  double LUMI = Utility->ReadLumi(ParameterFILE);
  Utility->MakeGraphVar(ParameterFILE,&ge,"Fl",true);
  ge->SetMarkerColor(kBlack);
  ge->SetMarkerStyle(20);
  Utility->MakeGraphVar(ParameterFILE,&geTMP,"Afb",true);

  for (int i = 0; i < geTMP->GetN(); i++)
    {
      ge->SetPoint(i,geTMP->GetY()[i],ge->GetY()[i]);
      ge->SetPointEXhigh(i,geTMP->GetErrorYhigh(i));
      ge->SetPointEXlow(i,geTMP->GetErrorYlow(i));
    }

  canv0->cd();
  histo->Draw();
  ge->Draw("pe1");

  line1 = new TLine(-3.0/4.0,0.0,0.0,1.0);
  line1->SetLineColor(kRed);
  line1->SetLineWidth(2);
  line1->Draw("same");

  line2 = new TLine(+3.0/4.0,0.0,0.0,1.0);
  line2->SetLineColor(kRed);
  line2->SetLineWidth(2);
  line2->Draw("same");

  DrawString(LUMI);

  canv0->Update();
}


void MakePhysicsPlots (unsigned int PlotType)
// ####################
// # PlotType:        #
// ####################
// # 0 = GEN-RECO-FL  #
// # 1 = GEN-RECO-AFB #
// # 2 = GEN-RECO-BF  #
// ####################
// # 10 = DATA-FL     #
// # 11 = DATA-AFB    #
// # 12 = DATA-BF     #
// ####################
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

  
  unsigned int DoF = 0;
  double myGlobalChi2 = 0.0;


  // ##################################
  // # Read q^2 bins from config file #
  // ##################################
  vector<double> q2Bins;
  Utility->Readq2Bins(ParameterFILE,&q2Bins);
  unsigned int JPsibin = Utility->GetJPsiBin(&q2Bins);
  unsigned int PsiPbin = Utility->GetPsiPBin(&q2Bins);
  double LUMI          = Utility->ReadLumi(ParameterFILE);
  double* q2Bins_      = Utility->MakeBinning(&q2Bins);


  TCanvas* canv0;
  if (FORPAPER == false) canv0 = new TCanvas("canv0","canv0",10,10,700,900);
  else                   canv0 = new TCanvas("canv0","canv0",10,10,700,500);
  TPad *pad1, *pad2, *pad3;
  TPaveText* paveText = NULL;
  TGraphAsymmErrors* ge0   = NULL;
  TGraphAsymmErrors* ge1  = NULL;
  TGraphAsymmErrors* ge00 = NULL;
  TGraphAsymmErrors* ge11 = NULL;
  TGraphBentErrors* geb   = NULL;
  TLine* line;


  // ############################################
  // # Variables to read values from ASCII file #
  // ############################################
  vector<double> vxs, vys;
  vector<double> vxel, vxeh;
  vector<double> vyel, vyeh;

 
  if (PlotType == 0) // Fl
    {
      Utility->MakeGraphVar(ParameterFILE_MCGEN,&ge0,"Fl",false);
      ge0->SetMarkerColor(kBlack);
      ge0->SetMarkerStyle(22);
      ge0->SetMarkerSize(1.2);
      ge0->SetFillColor(kWhite);
      ge0->SetFillStyle(0);
      ge0->SetLineColor(kBlack);
      ge0->SetLineWidth(2);
      ge0->GetXaxis()->SetRangeUser(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge0->GetYaxis()->SetRangeUser(-0.02,1.0);
      ge0->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});F_{L}");

      Utility->MakeGraphVar(ParameterFILE_MCRECO,&ge1,"Fl",false);
      ge1->SetMarkerColor(kBlue);
      ge1->SetMarkerStyle(20);
      ge1->SetMarkerSize(1.2);
      ge1->SetFillColor(kGreen);
      ge1->SetFillStyle(3001);
      ge1->SetLineColor(kBlue);
      ge1->SetLineWidth(2);
      ge1->GetXaxis()->SetRangeUser(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge1->GetYaxis()->SetRangeUser(-0.02,1.0);
      ge1->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});F_{L}");
     }
  else if (PlotType == 1) // Afb
    {
      Utility->MakeGraphVar(ParameterFILE_MCGEN,&ge0,"Afb",false);
      ge0->SetMarkerColor(kBlack);
      ge0->SetMarkerStyle(22);
      ge0->SetMarkerSize(1.2);
      ge0->SetFillColor(kWhite);
      ge0->SetFillStyle(0);
      ge0->SetLineColor(kBlack);
      ge0->SetLineWidth(2);
      ge0->GetXaxis()->SetRangeUser(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge0->GetYaxis()->SetRangeUser(-1.04,1.0);
      ge0->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});A_{FB}");
 
      Utility->MakeGraphVar(ParameterFILE_MCRECO,&ge1,"Afb",false);
      ge1->SetMarkerColor(kBlue);
      ge1->SetMarkerStyle(20);
      ge1->SetMarkerSize(1.2);
      ge1->SetFillColor(kGreen);
      ge1->SetFillStyle(3001);
      ge1->SetLineColor(kBlue);
      ge1->SetLineWidth(2);
      ge1->GetXaxis()->SetRangeUser(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge1->GetYaxis()->SetRangeUser(-1.04,1.0);
      ge1->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});A_{FB}");
    }
  else if (PlotType == 2) // Branching fraction
    {
      Utility->MakeGraphVar(ParameterFILE_MCGEN,&ge0,"dBFdq2",false);
      ge0->SetMarkerColor(kBlack);
      ge0->SetMarkerStyle(22);
      ge0->SetMarkerSize(1.2);
      ge0->SetFillColor(kWhite);
      ge0->SetFillStyle(0);
      ge0->SetLineColor(kBlack);
      ge0->SetLineWidth(2);
      ge0->GetXaxis()->SetRangeUser(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge0->GetYaxis()->SetRangeUser(0.0,1.2);
      ge0->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}7}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");
 
      Utility->MakeGraphVar(ParameterFILE_MCRECO,&ge1,"dBFdq2",false);
      ge1->SetMarkerColor(kBlue);
      ge1->SetMarkerStyle(20);
      ge1->SetMarkerSize(1.2);
      ge1->SetFillColor(kGreen);
      ge1->SetFillStyle(3001);
      ge1->SetLineColor(kBlue);
      ge1->SetLineWidth(2);
      ge1->GetXaxis()->SetRangeUser(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge1->GetYaxis()->SetRangeUser(0.0,1.2);
      ge1->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}7}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");
    }
  else if (PlotType == 10) // FL
    {
      Utility->MakeGraphVar(ParameterFILE,&ge0,"Fl",false);
      ge0->SetMarkerColor(kBlack);
      ge0->SetMarkerStyle(20);
      ge0->SetMarkerSize(1.2);
      ge0->SetLineColor(kBlack);
      ge0->SetLineWidth(2);
      ge0->GetXaxis()->SetRangeUser(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge0->GetYaxis()->SetRangeUser(-0.02,1.0);
      ge0->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^2});F_{L}");
    }
  else if (PlotType == 11) // Afb
    {
      Utility->MakeGraphVar(ParameterFILE,&ge0,"Afb",false);
      ge0->SetMarkerColor(kBlack);
      ge0->SetMarkerStyle(20);
      ge0->SetMarkerSize(1.2);
      ge0->SetLineColor(kBlack);
      ge0->SetLineWidth(2);
      ge0->GetXaxis()->SetRangeUser(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge0->GetYaxis()->SetRangeUser(-1.04,1.0);
      ge0->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});A_{FB}");
    }
  else if (PlotType == 12) // Branching fraction
    {
      Utility->MakeGraphVar(ParameterFILE,&ge0,"dBFdq2",false);
      ge0->SetMarkerColor(kBlack);
      ge0->SetMarkerStyle(20);
      ge0->SetMarkerSize(1.2);
      ge0->SetLineColor(kBlack);
      ge0->SetLineWidth(2);
      ge0->GetXaxis()->SetRangeUser(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge0->GetYaxis()->SetRangeUser(0.0,1.2);
      ge0->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}7}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");
    }
  else
    {
      cout << "[MakePlots::MakePhysicsPlots]\tWrong option number" << endl;
      exit(1);
    }


  // ##################################
  // # Read SM values from ASCII file #
  // ##################################
  TGraphAsymmErrors* geSmoothTh = NULL;
  if      ((PlotType == 0) || (PlotType == 10)) geSmoothTh = ReadFromASCII(SMFL,PlotType,&q2Bins,&vxs,&vys,&vxel,&vxeh,&vyel,&vyeh);  // Fl
  else if ((PlotType == 1) || (PlotType == 11)) geSmoothTh = ReadFromASCII(SMAFB,PlotType,&q2Bins,&vxs,&vys,&vxel,&vxeh,&vyel,&vyeh); // Afb
  else if ((PlotType == 2) || (PlotType == 12)) geSmoothTh = ReadFromASCII(SMBF,PlotType,&q2Bins,&vxs,&vys,&vxel,&vxeh,&vyel,&vyeh);  // Branching fraction
  // geSmoothTh->SetFillColor(kCyan-4);
  geSmoothTh->SetFillColor(kRed-9);
  geSmoothTh->SetFillStyle(1001);
  geSmoothTh->GetXaxis()->SetRangeUser(q2Bins[0],q2Bins[q2Bins.size()-1]);


  // ################################################
  // # Average theory over q^2 bins from ASCII file #
  // ################################################
  TGraphAsymmErrors* geStepTh = NULL;
  if      ((PlotType == 0) || (PlotType == 10)) geStepTh = ReadFromASCII(SMBINFL,PlotType,&q2Bins,&vxs,&vys,&vxel,&vxeh,&vyel,&vyeh);  // Fl
  else if ((PlotType == 1) || (PlotType == 11)) geStepTh = ReadFromASCII(SMBINAFB,PlotType,&q2Bins,&vxs,&vys,&vxel,&vxeh,&vyel,&vyeh); // Afb
  else if ((PlotType == 2) || (PlotType == 12)) geStepTh = ReadFromASCII(SMBINBF,PlotType,&q2Bins,&vxs,&vys,&vxel,&vxeh,&vyel,&vyeh);  // Branching fraction
  geStepTh->SetMarkerColor(kBlack);
  geStepTh->SetMarkerStyle(1);
  geStepTh->SetFillColor(kBlue);
  geStepTh->SetFillStyle(3001);
  geStepTh->GetXaxis()->SetRangeUser(q2Bins[0],q2Bins[q2Bins.size()-1]);


  // ############################
  // # Adding systematic errors #
  // ############################
  vector<vector<double>*> vecObs; // Vector containing the pointers to the vectors containing the fit-observable systematic errors
  Utility->ReadFitSystematics(ParameterFILE,&vecObs);
  if (PlotType == 0) // Fl
    {
      // ############################################
      // # Graph containing only statistical errors #
      // ############################################
      ge11 = new TGraphAsymmErrors(*ge1);
      ge11->SetMarkerColor(kBlue);
      ge11->SetMarkerStyle(20);
      ge11->SetMarkerSize(1.2);
      ge11->SetFillColor(kGreen);
      ge11->SetFillStyle(3001);
      ge11->SetLineColor(kBlue);
      ge11->SetLineWidth(2);


      // ######################################################
      // # Graph containing statistical and systematic errors #
      // ######################################################
      for (int i = 0; i < ge1->GetN(); i++)
      	if (Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false)
      	  {
      	    ge1->SetPointEYhigh(i,sqrt(ge1->GetErrorYhigh(i)*ge1->GetErrorYhigh(i) + vecObs[0]->operator[](i)*vecObs[0]->operator[](i)));
      	    ge1->SetPointEYlow(i, sqrt(ge1->GetErrorYlow(i)*ge1->GetErrorYlow(i)   + vecObs[1]->operator[](i)*vecObs[1]->operator[](i)));
      	  }
    }
  else if (PlotType == 1) // Afb
    {
      // ############################################
      // # Graph containing only statistical errors #
      // ############################################
      ge11 = new TGraphAsymmErrors(*ge1);
      ge11->SetMarkerColor(kBlue);
      ge11->SetMarkerStyle(20);
      ge11->SetMarkerSize(1.2);
      ge11->SetFillColor(kGreen);
      ge11->SetFillStyle(3001);
      ge11->SetLineColor(kBlue);
      ge11->SetLineWidth(2);


      // ######################################################
      // # Graph containing statistical and systematic errors #
      // ######################################################
      for (int i = 0; i < ge1->GetN(); i++)
	if (Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false)
	  {
	    ge1->SetPointEYhigh(i,sqrt(ge1->GetErrorYhigh(i)*ge1->GetErrorYhigh(i) + vecObs[2]->operator[](i)*vecObs[2]->operator[](i)));
	    ge1->SetPointEYlow(i, sqrt(ge1->GetErrorYlow(i)*ge1->GetErrorYlow(i)   + vecObs[3]->operator[](i)*vecObs[3]->operator[](i)));
	  }
    }
  else if (PlotType == 2) // Branching fraction
    {
      // ############################################
      // # Graph containing only statistical errors #
      // ############################################
      ge11 = new TGraphAsymmErrors(*ge1);
      ge11->SetMarkerColor(kBlue);
      ge11->SetMarkerStyle(20);
      ge11->SetMarkerSize(1.2);
      ge11->SetFillColor(kGreen);
      ge11->SetFillStyle(3001);
      ge11->SetLineColor(kBlue);
      ge11->SetLineWidth(2);


      // ######################################################
      // # Graph containing statistical and systematic errors #
      // ######################################################
      for (int i = 0; i < ge1->GetN(); i++)
	if (Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false)
	  {
	    ge1->SetPointEYhigh(i,sqrt(ge1->GetErrorYhigh(i)*ge1->GetErrorYhigh(i) + vecObs[8]->operator[](i)*vecObs[8]->operator[](i)));
	    ge1->SetPointEYlow(i, sqrt(ge1->GetErrorYlow(i)*ge1->GetErrorYlow(i)   + vecObs[9]->operator[](i)*vecObs[9]->operator[](i)));
	  }
    }
  else if (PlotType == 10) // Fl
    {
      geSmoothTh->GetYaxis()->SetRangeUser(-0.02,1.0);
      geSmoothTh->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});F_{L}");


      // ############################################
      // # Graph containing only statistical errors #
      // ############################################
      ge00 = new TGraphAsymmErrors(*ge0);
      ge00->SetMarkerColor(kBlack);
      ge00->SetMarkerStyle(20);
      ge00->SetMarkerSize(1.2);
      ge00->SetLineColor(kBlack);
      ge00->SetLineWidth(2);


      // ######################################################
      // # Graph containing statistical and systematic errors #
      // ######################################################
      for (int i = 0; i < ge0->GetN(); i++)
	if (Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false)
	  {
	    ge0->SetPointEXlow(i,0.0);
	    ge0->SetPointEXhigh(i,0.0);

	    ge0->SetPointEYhigh(i,sqrt(ge0->GetErrorYhigh(i)*ge0->GetErrorYhigh(i) + vecObs[0]->operator[](i)*vecObs[0]->operator[](i)));
	    ge0->SetPointEYlow(i, sqrt(ge0->GetErrorYlow(i)*ge0->GetErrorYlow(i)   + vecObs[1]->operator[](i)*vecObs[1]->operator[](i)));
	  }
    }
  else if (PlotType == 11) // Afb
    {
      geSmoothTh->GetYaxis()->SetRangeUser(-1.04,1.0);
      geSmoothTh->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});A_{FB}");


      // ############################################
      // # Graph containing only statistical errors #
      // ############################################
      ge00 = new TGraphAsymmErrors(*ge0);
      ge00->SetMarkerColor(kBlack);
      ge00->SetMarkerStyle(20);
      ge00->SetMarkerSize(1.2);
      ge00->SetLineColor(kBlack);
      ge00->SetLineWidth(2);


      // ######################################################
      // # Graph containing statistical and systematic errors #
      // ######################################################
      for (int i = 0; i < ge0->GetN(); i++)
	if (Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false)
	  {
	    ge0->SetPointEXlow(i,0.0);
	    ge0->SetPointEXhigh(i,0.0);

	    ge0->SetPointEYhigh(i,sqrt(ge0->GetErrorYhigh(i)*ge0->GetErrorYhigh(i) + vecObs[2]->operator[](i)*vecObs[2]->operator[](i)));
	    ge0->SetPointEYlow(i, sqrt(ge0->GetErrorYlow(i)*ge0->GetErrorYlow(i)   + vecObs[3]->operator[](i)*vecObs[3]->operator[](i)));
	  }
    }
  else if (PlotType == 12) // Branching fraction
    {
      geSmoothTh->GetYaxis()->SetRangeUser(0.0,1.2);
      geSmoothTh->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}7}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");


      // ############################################
      // # Graph containing only statistical errors #
      // ############################################
      ge00 = new TGraphAsymmErrors(*ge0);
      ge00->SetMarkerColor(kBlack);
      ge00->SetMarkerStyle(20);
      ge00->SetMarkerSize(1.2);
      ge00->SetLineColor(kBlack);
      ge00->SetLineWidth(2);


      // ######################################################
      // # Graph containing statistical and systematic errors #
      // ######################################################
      for (int i = 0; i < ge0->GetN(); i++)
	if (Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false)
	  {
	    ge0->SetPointEXlow(i,0.0);
	    ge0->SetPointEXhigh(i,0.0);

	    ge0->SetPointEYhigh(i,sqrt(ge0->GetErrorYhigh(i)*ge0->GetErrorYhigh(i) + vecObs[8]->operator[](i)*vecObs[8]->operator[](i)));
	    ge0->SetPointEYlow(i, sqrt(ge0->GetErrorYlow(i)*ge0->GetErrorYlow(i)   + vecObs[9]->operator[](i)*vecObs[9]->operator[](i)));
	  }


      // #################################################################
      // # Divide the theoretical branching fraction by the q2 bin width #
      // #################################################################
      for (unsigned int i = 0; i < q2Bins.size()-1; i++)
        if ((Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false) && (Utility->ValIsBetweenJPsiAndPsiP(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false))
	  {
	    geStepTh->SetPoint(i,geStepTh->GetX()[i],geStepTh->GetY()[i] /  (q2Bins[i+1] - q2Bins[i]));
	    geStepTh->SetPointEYlow(i,geStepTh->GetErrorYlow(i) / (q2Bins[i+1] - q2Bins[i]));
	    geStepTh->SetPointEYhigh(i,geStepTh->GetErrorYhigh(i) / (q2Bins[i+1] - q2Bins[i]));
	  }
    }
  

  // ################################
  // # Pad for the actual histogram #
  // ################################
  canv0->cd();
  if (FORPAPER == false)
    {
      pad1 = new TPad("pad1","pad1",0,0.5,1,1);
      pad1->SetBottomMargin(0);
    }
  else pad1 = new TPad("pad1","pad1",0,0,1,1);
  pad1->Draw();
  pad1->cd();


  TLegend* leg = new TLegend(0.8, 0.7, 0.97, 0.89, "");
  if ((PlotType == 0) || (PlotType == 1) || (PlotType == 2)) // Fl OR Afb OR Branching fraction
    {
      ge0->Draw("ape1");
      ge1->Draw("same pe1");
      ge11->Draw("same e2");

      leg->AddEntry(ge0,"GEN-MC","PL");
      leg->AddEntry(ge11,"RECO-MC","EPFL");
    }
  else if ((PlotType == 10) || (PlotType == 11) || (PlotType == 12)) // Fl OR Afb OR Branching fraction
    {
      geSmoothTh->Draw("ae2");
      geStepTh->Draw("same e2");
      
      ge00->Draw("same pe1");
      // ###################################
      // # Code to make slanted error bars #
      // ###################################
      vector<double> exld;
      vector<double> exhd;
      vector<double> eyld;
      vector<double> eyhd;
      for (int i = 0; i < ge0->GetN(); i++)
	{
	  exld.push_back(0.0);
	  exhd.push_back(0.0);
	  eyld.push_back(0.4);
	  eyhd.push_back(-0.4);
	}
      geb = new TGraphBentErrors(ge0->GetN(), ge0->GetX(), ge0->GetY(), ge0->GetEXlow(), ge0->GetEXhigh(), ge0->GetEYlow(), ge0->GetEYhigh(), &exld[0], &exhd[0], &eyld[0], &eyhd[0]);
      geb->SetMarkerColor(kBlack);
      geb->SetMarkerStyle(20);
      geb->SetMarkerSize(1.2);
      geb->SetLineColor(kBlack);
      geb->SetLineWidth(2);
      geb->SetLineStyle(kDashed);
      exld.clear();
      exhd.clear();
      eyld.clear();
      eyhd.clear();
      // geb->Draw("same pe1");
      ge0->Draw("same pez");

      leg->AddEntry(ge00,"Data","EPL");
      leg->AddEntry(geSmoothTh,"SM","F");
      leg->AddEntry(geStepTh,"<SM>","F");
    }
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
  

  // ########################
  // # Draw exclusion zones #
  // ########################
  DrawExclusion(q2Bins[JPsibin],q2Bins[JPsibin+1],-1.2,1.2,"RejectJPsi1",3001,kGray);
  DrawExclusion(q2Bins[PsiPbin],q2Bins[PsiPbin+1],-1.2,1.2,"RejectPsiP1",3001,kGray);


  // ####################
  // # Pad for the chi2 #
  // ####################
  canv0->cd();
  pad2 = new TPad("pad2","pad2",0,0.25,1,0.5);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0);
  pad2->SetLogy();
  if (FORPAPER == false) pad2->Draw();
  pad2->cd();


  TH1D* chi2Histo = new TH1D("chi2Histo","chi2Histo",q2Bins.size()-1,q2Bins_);
  TLegend* ratioLeg = new TLegend(0.8, 0.85, 0.97, 0.95, "");
  if ((PlotType == 0) || (PlotType == 1) || (PlotType == 2)) // Fl OR Afb OR Branching fraction
    {
      double tmpVar;
      for (int i = 0; i < chi2Histo->GetNbinsX(); i++)
	{
	  if (Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false)
	    {
	      tmpVar = pow(ge0->GetY()[i] - ge1->GetY()[i],2.) /
		(ge0->GetY()[i] > ge1->GetY()[i] ? pow(ge0->GetErrorYlow(i),2.) + pow(ge1->GetErrorYhigh(i),2.) : pow(ge0->GetErrorYhigh(i),2.) + pow(ge1->GetErrorYlow(i),2.));

	      myGlobalChi2 = myGlobalChi2 + tmpVar;
	      chi2Histo->SetBinContent(i+1,tmpVar);

	      DoF++;
	    }
	  else chi2Histo->SetBinContent(i+1,0.0);
	}

      ratioLeg->AddEntry(chi2Histo,"#chi^{2}(RECO, GEN)");
    }
  else if ((PlotType == 10) || (PlotType == 11) || (PlotType == 12)) // Fl OR Afb OR Branching fraction
    {      
      double tmpVar;
      for (int i = 0; i < chi2Histo->GetNbinsX(); i++)
	{
	  if ((Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false) && (Utility->ValIsBetweenJPsiAndPsiP(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false))
	    {
	      tmpVar = pow(geStepTh->GetY()[i] - ge0->GetY()[i],2.) /
		(geStepTh->GetY()[i] > ge0->GetY()[i] ? pow(geStepTh->GetErrorYlow(i),2.) + pow(ge0->GetErrorYhigh(i),2.) : pow(geStepTh->GetErrorYhigh(i),2.) + pow(ge0->GetErrorYlow(i),2.));
	      
	      myGlobalChi2 = myGlobalChi2 + tmpVar;
	      chi2Histo->SetBinContent(i+1,tmpVar);

	      DoF++;
	    }
	  else chi2Histo->SetBinContent(i+1,0.0);
	}
      
      ratioLeg->AddEntry(chi2Histo,"#chi^{2}(<SM>, Data)");
    }
  cout << "\n@@@ Global chi2 = " << myGlobalChi2 / ((double)DoF) << " (" << myGlobalChi2 << "/" << ((double)DoF) << ") @@@" << endl;
  myGlobalChi2 = myGlobalChi2 / ((double)DoF);

  chi2Histo->SetTitle(";q^{2} (GeV^{2});#chi^{2}");
  chi2Histo->GetXaxis()->SetLabelSize(0.06);
  chi2Histo->GetXaxis()->SetTitleOffset(0.8);
  chi2Histo->GetXaxis()->SetTitleSize(0.07);
  chi2Histo->GetYaxis()->SetLabelSize(0.06);
  chi2Histo->GetYaxis()->SetTitleOffset(0.7);
  chi2Histo->GetYaxis()->SetTitleSize(0.07);
  chi2Histo->SetFillColor(kAzure+6);
  chi2Histo->SetFillStyle(1001);
  chi2Histo->GetYaxis()->SetRangeUser(0.008,30.0);
  chi2Histo->Draw();

  ratioLeg->SetFillColor(0);
  ratioLeg->SetBorderSize(1);
  ratioLeg->Draw();


  // ########################
  // # Draw exclusion zones #
  // ########################
  DrawExclusion(q2Bins[JPsibin],q2Bins[JPsibin+1],0.008,30.0,"RejectJPsi2",3001,kGray);
  DrawExclusion(q2Bins[PsiPbin],q2Bins[PsiPbin+1],0.008,30.0,"RejectPsiP2",3001,kGray);


  // ###########################
  // # Pad for the probability #
  // ###########################
  canv0->cd();
  pad3 = new TPad("pad3","pad3",0,0.0,1,0.25);
  pad3->SetTopMargin(0);
  pad3->SetBottomMargin(0.15);
  pad3->SetLogy();
  if (FORPAPER == false) pad3->Draw();
  pad3->cd();

  TH1D* probHisto = new TH1D("probHisto","probHisto",q2Bins.size()-1,q2Bins_);
  TLegend* probLeg = new TLegend(0.8, 0.85, 0.97, 0.95, "");
  for (int i = 0; i < chi2Histo->GetNbinsX(); i++)
    if ((((PlotType == 10) || (PlotType == 11) || (PlotType == 12)) && ((Utility->ValIsBetweenJPsiAndPsiP(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false))) ||
	!((PlotType == 10) || (PlotType == 11) || (PlotType == 12)))
      {
	if (Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false)
	  {
	    probHisto->SetBinContent(i+1,TMath::Prob(chi2Histo->GetBinContent(i+1),1.0));
	  }
      }
    else probHisto->SetBinContent(i+1,0.0);
  probLeg->AddEntry(probHisto,"p-value");

  probHisto->SetTitle(";q^{2} (GeV^{2});p-value");
  probHisto->GetXaxis()->SetLabelSize(0.06);
  probHisto->GetXaxis()->SetTitleOffset(0.8);
  probHisto->GetXaxis()->SetTitleSize(0.07);
  probHisto->GetYaxis()->SetLabelSize(0.06);
  probHisto->GetYaxis()->SetTitleOffset(0.7);
  probHisto->GetYaxis()->SetTitleSize(0.07);
  probHisto->SetFillColor(kGreen-7);
  probHisto->SetFillStyle(1001);
  probHisto->GetYaxis()->SetRangeUser(1e-3,1.3);
  probHisto->Draw();

  probLeg->SetFillColor(0);
  probLeg->SetBorderSize(1);
  probLeg->Draw();


  // ########################
  // # Draw exclusion zones #
  // ########################
  DrawExclusion(q2Bins[JPsibin],q2Bins[JPsibin+1],1e-3,1.13,"RejectJPsi3",3001,kGray);
  DrawExclusion(q2Bins[PsiPbin],q2Bins[PsiPbin+1],1e-3,1.13,"RejectPsiP3",3001,kGray);
  DrawExclusion(q2Bins[0],q2Bins[q2Bins.size()-1],1e-3,0.05,"CL95",3001,kRed-7);


  // #################################
  // # Write global chi2 and p-value #
  // #################################
  pad1->cd();
  if (FORPAPER == false)
    {
      paveText = new TPaveText(0.8,0.05,0.97,0.15,"NDC");
      paveText->SetTextAlign(11);
      paveText->SetBorderSize(0.0);
      paveText->SetFillColor(kWhite);
      paveText->AddText(Form("%s%.2f","#chi^{2}/DoF = ",myGlobalChi2));
      paveText->AddText(Form("%s%.3f","p-value = ",TMath::Prob(myGlobalChi2*(double(DoF)),DoF)));
      paveText->SetFillStyle(0);
      paveText->SetTextSize(0.035);
      paveText->Draw();
    }
  DrawString(LUMI);


  // #################################
  // # Draw horizontal line at y = 0 #
  // #################################
  if ((PlotType == 1) ||(PlotType == 11)) // Afb
    {
      line = new TLine(q2Bins[0],0.0,q2Bins[q2Bins.size()-1],0.0);
      line->SetLineStyle(kDashed);
      line->Draw();
    }
  

  canv0->cd();
  canv0->Update();
}


void EvalMultyRun (unsigned int sysType, unsigned int q2BinIndx, string fileName, double NLLinterval, double NLLlessThan)
// #########################
// # sysType:              #
// # 0 = Fl                #
// # 1 = Afb               #
// # 2 = BF                #
// # 3 = All               #
// # 10 = Fl multy minima  #
// # 11 = Afb multy minima #
// # 12 = BF multy minima  #
// #########################
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetTitleOffset(1.25,"y");
  TGaxis::SetMaxDigits(3);


  unsigned int nBinsHisto = 50;
  ifstream inputFile;
  stringstream myString;


  TCanvas* c1 = new TCanvas("c1","c1",10,10,700,500);
  c1->cd();
  TH1D* h1 = new TH1D("h1","h1",nBinsHisto/2.,0.0,1.0);
  h1->GetXaxis()->SetTitle("F_{L}");
  h1->GetYaxis()->SetTitle("Entries");
  h1->SetFillColor(kAzure+6);

  TCanvas* c2 = new TCanvas("c2","c2",10,10,700,500);
  c2->cd();
  TH1D* h2 = new TH1D("h2","h2",nBinsHisto,-1.0,1.0);
  h2->GetXaxis()->SetTitle("A_{FB}");
  h2->GetYaxis()->SetTitle("Entries");
  h2->SetFillColor(kAzure+6);
  
  TCanvas* c3 = new TCanvas("c3","c3",10,10,700,500);
  c3->cd();
  TH1D* h3 = new TH1D("h3","h3",nBinsHisto,1.0,-1.0);
  h3->GetXaxis()->SetTitle("Signal yield OR dBF/dq^{2} (GeV^{#font[122]{\55}2})");
  h3->GetYaxis()->SetTitle("Entries");
  h3->SetFillColor(kAzure+6);

  TCanvas* c4 = new TCanvas("c4","c4",10,10,700,500);
  c4->cd();
  TH1D* h4 = new TH1D("h4","h4",nBinsHisto,1.0,-1.0);
  h4->GetXaxis()->SetTitle("NLL");
  h4->GetYaxis()->SetTitle("Entries");
  h4->SetFillColor(kAzure+6);

  TCanvas* cSc = new TCanvas("cSc","cSc",10,10,700,500);
  cSc->cd();
  TH2D* hSc = new TH2D("hSc","hSc",nBinsHisto,1.0,-1.0,nBinsHisto,1.0,-1.0);
  hSc->GetXaxis()->SetTitle("Parameter Value");
  hSc->GetYaxis()->SetTitle("NLL");
  hSc->GetZaxis()->SetTitle("Entries");
  hSc->SetFillColor(kAzure+6);


  // #########################
  // # Read values from file #
  // #########################
  fileName.erase(fileName.find(".txt"),4);
  myString.str("");
  if      (sysType == 0) myString << fileName << "_FL_"  << q2BinIndx << ".txt";
  else if (sysType == 1) myString << fileName << "_AFB_" << q2BinIndx << ".txt";
  else if (sysType == 2) myString << fileName << "_BF_"  << q2BinIndx << ".txt";
  else                   myString << fileName << ".txt";
  cout << "Opening efficiency-systematics file: " << myString.str().c_str() << endl;
  inputFile.open(myString.str().c_str(), ifstream::in);
  if (inputFile.good() == false)
    {
      cout << "[MakePlots::EvalMultyRun]\tError opening file : " << myString.str().c_str() << endl;
      exit (1);
    }


  double var0 = 0.0;
  double var1 = 0.0;
  double var2 = 0.0;
  double var3 = 0.0;
  double var4 = 0.0;
  inputFile >> var0 >> var1 >> var2 >> var3 >> var4;
  while (inputFile)
    {
      if (var1 != -2.0) h1->Fill(var1);
      if (var2 != -2.0) h2->Fill(var2);
      if (var3 != -2.0) h3->Fill(var3);
      if (var4 != -2.0) h4->Fill(var4);

      if (var4 < NLLlessThan)
	{
	  if ((sysType == 10) && (var1 != -2.0))
	    {
	      hSc->GetXaxis()->SetTitle("F_{L}");
	      hSc->Fill(var1,var4);
	    }
	  else if ((sysType == 11) && (var2 != -2.0))
	    {
	      hSc->GetXaxis()->SetTitle("A_{FB}");
	      hSc->Fill(var2,var4);
	    }
	  else if ((sysType == 12) && (var3 != -2.0))
	    {
	      hSc->GetXaxis()->SetTitle("dBF/dq^{2} (GeV^{#font[122]{\55}2})");
	      hSc->Fill(var3,var4);
	    }
	  
	  cout << "var0: " << var0 << "\tvar1: " << var1 << "\tvar2: " << var2 << "\tvar3: " << var3 << "\tvar4: " << var4 << endl;
	}
      
      inputFile >> var0 >> var1 >> var2 >> var3 >> var4;
    }


  // #####################################
  // # Rebin NLL scatter plot and refill #
  // #####################################
  if (sysType == 10)      hSc->SetBins(nBinsHisto,0.0,1.0,20,hSc->GetMean(2) - NLLinterval,hSc->GetMean(2) + NLLinterval);
  else if (sysType == 11) hSc->SetBins(nBinsHisto*2,-1.0,1.0,20,hSc->GetMean(2) - NLLinterval,hSc->GetMean(2) + NLLinterval);
  else if (sysType == 12) hSc->SetBins(nBinsHisto,1.0,-1.0,20,hSc->GetMean(2) - NLLinterval,hSc->GetMean(2) + NLLinterval);
  if ((sysType >= 10) && (sysType <= 12))
    {
      cout << "\n@@@ Rebin NLL scatter plot and refill @@@" << endl;
      hSc->Reset();

      inputFile.clear();
      inputFile.seekg(0,inputFile.beg);
      inputFile >> var0 >> var1 >> var2 >> var3 >> var4;
      while (inputFile)
	{
	  
	  if (var4 < NLLlessThan)
	    {
	      if ((sysType == 10) && (var1 != -2.0))
		{
		  hSc->GetXaxis()->SetTitle("F_{L}");
		  hSc->Fill(var1,var4);
		}
	      else if ((sysType == 11) && (var2 != -2.0))
		{
		  hSc->GetXaxis()->SetTitle("A_{FB}");
		  hSc->Fill(var2,var4);
		}
	      else if ((sysType == 12) && (var3 != -2.0))
		{
		  hSc->GetXaxis()->SetTitle("dBF/dq^{2} (GeV^{#font[122]{\55}2})");
		  hSc->Fill(var3,var4);
		}
	  
	      cout << "var0: " << var0 << "\tvar1: " << var1 << "\tvar2: " << var2 << "\tvar3: " << var3 << "\tvar4: " << var4 << endl;
	    }
      
	  inputFile >> var0 >> var1 >> var2 >> var3 >> var4;
	}
    }


  cout << "\n@@@ I'm now fitting and plotting @@@" << endl;

  c1->cd();
  h1->Draw();
  h1->Fit("gaus");
  c1->Update();

  c2->cd();
  h2->Draw();
  h2->Fit("gaus");
  c2->Update();

  c3->cd();
  h3->Draw();
  h3->Fit("gaus");
  c3->Update();

  c4->cd();
  h4->Draw();
  h4->Fit("gaus");
  c4->Update();

  cSc->cd();
  hSc->Draw("gcolz");
  cSc->Update();


  inputFile.close();
}


void PlotMuMu (string fileName, bool bkgSub)
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


  int nEntries;
  double minX = 0.8;
  double maxX = 5.0;
  unsigned int nBins = 100;
  stringstream myString;
  string sigMassQuery = "";
  string bkgMassQuery = "";

  double signalSigma = sqrt(Utility->GetGenericParam("FRACMASSS") * Utility->GetGenericParam("SIGMAS1") * Utility->GetGenericParam("SIGMAS1") +
			    (1. - Utility->GetGenericParam("FRACMASSS")) * Utility->GetGenericParam("SIGMAS2") * Utility->GetGenericParam("SIGMAS2"));
  cout << "\n@@@ Signal sigma: " << signalSigma << " @@@" << endl;


  if (bkgSub == true)
    {
      myString.clear();
      myString.str("");
      myString << "(abs(B0MassArb - " << Utility->B0Mass << ") < " << Utility->GetGenericParam("NSigmaB0S")*signalSigma << ")";
      sigMassQuery = myString.str();
      
      myString.clear();
      myString.str("");
      myString << "((B0MassArb > " << Utility->B0Mass+Utility->GetGenericParam("NSigmaB0B")*signalSigma << " && B0MassArb < "
	       << Utility->B0Mass+(Utility->GetGenericParam("NSigmaB0B")+Utility->GetGenericParam("NSigmaB0S"))*signalSigma << ") || ";
      myString << "(B0MassArb < " << Utility->B0Mass-Utility->GetGenericParam("NSigmaB0B")*signalSigma << " && B0MassArb > "
	       << Utility->B0Mass-(Utility->GetGenericParam("NSigmaB0B")+Utility->GetGenericParam("NSigmaB0S"))*signalSigma << "))";
      bkgMassQuery = myString.str(); 
    }
  else
    {
      myString.clear();
      myString.str("");
      myString << "(B0MassArb > " << Utility->B0Mass - Utility->GetGenericParam("B0MassIntervalLeft") << " && B0MassArb < " << Utility->B0Mass + Utility->GetGenericParam("B0MassIntervalRight") << ")";
      sigMassQuery = myString.str();

      myString.clear();
      myString.str("");
      myString << "(B0MassArb > " << Utility->B0Mass - Utility->GetGenericParam("B0MassIntervalLeft") << " && B0MassArb < " << Utility->B0Mass + Utility->GetGenericParam("B0MassIntervalRight") << ")";
      bkgMassQuery = myString.str(); 
    }


  TFile* file0 = TFile::Open(fileName.c_str(),"READ");
  TTree* theTree = (TTree*)file0->Get("B0SingleCand/B0KstMuMuSingleCandNTuple");

  nEntries = theTree->GetEntries();
  cout << "\n@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;


  TCanvas* c0 = new TCanvas("c0","c0",10,10,700,500);

  TH1D* hDsig = new TH1D("hDsig","hDsig",nBins,minX,maxX);
  hDsig->SetXTitle("M(#mu^{#font[122]{+}}#mu^{#font[122]{\55}}) (GeV)");
  hDsig->SetYTitle("Entries / (0.042 GeV)");
  hDsig->SetMarkerStyle(20);

  TH1D* hDbkg = new TH1D("hDbkg","hDbkg",nBins,minX,maxX);
  hDbkg->SetXTitle("M(#mu^{#font[122]{+}}#mu^{#font[122]{\55}}) (GeV)");
  hDbkg->SetYTitle("Entries / (0.042 GeV)");
  hDbkg->SetMarkerStyle(20);

  theTree->Draw("mumuMass>>hDsig",sigMassQuery.c_str(),"goff");
  theTree->Draw("mumuMass>>hDbkg",bkgMassQuery.c_str(),"goff");

  if (bkgSub == true) hDsig->Add(hDbkg, -1.0);
 
  c0->cd();
  hDsig->Draw("e1p");
  c0->Update();
}


void PlotKst (string fileName, bool bkgSub)
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


  int nEntries;
  double minX  = 0.82;
  double maxX  = 0.97;
  double extra = 0.02;
  unsigned int nBins = 100;
  stringstream myString;
  string sigMassQuery = "";
  string bkgMassQuery = "";
  string tmpstring    = "";

  double signalSigma = sqrt(Utility->GetGenericParam("FRACMASSS") * Utility->GetGenericParam("SIGMAS1") * Utility->GetGenericParam("SIGMAS1") +
			    (1. - Utility->GetGenericParam("FRACMASSS")) * Utility->GetGenericParam("SIGMAS2") * Utility->GetGenericParam("SIGMAS2"));
  cout << "\n@@@ Signal sigma: " << signalSigma << " @@@" << endl;


  if (bkgSub == true)
    {
      myString.clear();
      myString.str("");
      myString << "(abs(B0MassArb - " << Utility->B0Mass << ") < " << Utility->GetGenericParam("NSigmaB0S")*signalSigma << ")";
      sigMassQuery = myString.str();

      myString.clear();
      myString.str("");
      myString << "((B0MassArb > " << Utility->B0Mass+Utility->GetGenericParam("NSigmaB0B")*signalSigma << " && B0MassArb < "
	       << Utility->B0Mass+(Utility->GetGenericParam("NSigmaB0B")+Utility->GetGenericParam("NSigmaB0S"))*signalSigma << ") || ";
      myString << "(B0MassArb < " << Utility->B0Mass-Utility->GetGenericParam("NSigmaB0B")*signalSigma << " && B0MassArb > "
	       << Utility->B0Mass-(Utility->GetGenericParam("NSigmaB0B")+Utility->GetGenericParam("NSigmaB0S"))*signalSigma << "))";
      bkgMassQuery = myString.str(); 
    }
  else
    {
      myString.clear();
      myString.str("");
      myString << "(B0MassArb > " << Utility->B0Mass - Utility->GetGenericParam("B0MassIntervalLeft") << " && B0MassArb < " << Utility->B0Mass + Utility->GetGenericParam("B0MassIntervalRight") << ")";
      sigMassQuery = myString.str();

      myString.clear();
      myString.str("");
      myString << "(B0MassArb > " << Utility->B0Mass - Utility->GetGenericParam("B0MassIntervalLeft") << " && B0MassArb < " << Utility->B0Mass + Utility->GetGenericParam("B0MassIntervalRight") << ")";
      bkgMassQuery = myString.str();
    }


  TCanvas* c0 = new TCanvas("c0","c0",10,10,700,500);
  TCanvas* c1 = new TCanvas("c1","c1",20,20,700,500);

  TFile* file0 = TFile::Open(fileName.c_str(),"READ");
  TTree* theTree = (TTree*)file0->Get("B0SingleCand/B0KstMuMuSingleCandNTuple");

  nEntries = theTree->GetEntries();
  cout << "\n@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;

  TH1D* h1Dsig = new TH1D("h1Dsig","h1Dsig",nBins,minX - extra,maxX + extra+0.01);
  h1Dsig->SetXTitle("M(#font[122]{K}^{*0}) (GeV)");
  h1Dsig->SetYTitle("Entries / (0.002 GeV)");
  h1Dsig->SetMarkerStyle(20);

  TH1D* h2Dsig = new TH1D("h2Dsig","h2Dsig",nBins,minX - extra,maxX + extra+0.01);
  h2Dsig->SetXTitle("M(#font[122]{K}^{*0}) (GeV)");
  h2Dsig->SetYTitle("Entries / (0.002 GeV)");
  h2Dsig->SetMarkerStyle(21);
  h2Dsig->SetMarkerColor(kRed);


  TH1D* h1Dbkg = new TH1D("h1Dbkg","h1Dbkg",nBins,minX - extra,maxX + extra+0.01);
  h1Dbkg->SetXTitle("M(#font[122]{K}^{*0}) (GeV)");
  h1Dbkg->SetYTitle("Entries / (0.002 GeV)");
  h1Dbkg->SetMarkerStyle(20);

  TH1D* h2Dbkg = new TH1D("h2Dbkg","h2Dbkg",nBins,minX - extra,maxX + extra+0.01);
  h2Dbkg->SetXTitle("M(#font[122]{K}^{*0}) (GeV)");
  h2Dbkg->SetYTitle("Entries / (0.002 GeV)");
  h2Dbkg->SetMarkerStyle(21);
  h2Dbkg->SetMarkerColor(kRed);

  TF1* fsig = new TF1("fsig","[9]*([6]*TMath::Gaus(x,[0],[1]) + [7]*TMath::Gaus(x,[2],[3]) + [8]*TMath::Gaus(x,[4],[5]))",minX,maxX);
  TF1* fbkg = new TF1("fbkg","[0]",minX,maxX);
  TF1* f0 = new TF1("f0","fsig + fbkg",minX,maxX);

  fsig->SetLineColor(kBlue);
  fsig->SetFillColor(kAzure+6);
  fsig->SetFillStyle(3345);

  fbkg->SetLineColor(kRed);
  fbkg->SetFillColor(kRed);
  fbkg->SetFillStyle(3354);

  f0->SetParName(0,"#mu-1");
  f0->SetParName(1,"#sigma-1");
  f0->SetParName(2,"#mu-2");
  f0->SetParName(3,"#sigma-2");
  f0->SetParName(4,"#mu-3");
  f0->SetParName(5,"#sigma-3");
  f0->SetParName(6,"Ampli-1");
  f0->SetParName(7,"Ampli-2");
  f0->SetParName(8,"Ampli-3");

  f0->SetParName(9,"Ampli-S");
  f0->SetParName(10,"Ampli-B");


  // #######################################
  // # From fit to RECO B0 --> J/psi K* MC #
  // #######################################
  f0->FixParameter(0,8.93802e-01);
  f0->FixParameter(1,1.40485e-02);
  f0->FixParameter(2,8.96142e-01);
  f0->FixParameter(3,3.59670e-02);
  f0->FixParameter(4,9.69695e-01);
  f0->FixParameter(5,1.59102e-02);
  f0->FixParameter(6,6.93082e+02);
  f0->FixParameter(7,1.04683e+03);
  f0->FixParameter(8,1.06662e+02);

  // ####################
  // # Starting valiues #
  // ####################
  // f0->SetParameter(0,0.88);
  // f0->SetParameter(1,0.015);
  // f0->SetParameter(2,0.88);
  // f0->SetParameter(3,0.04);
  // f0->SetParameter(4,0.95);
  // f0->SetParameter(5,0.015);
  // f0->SetParameter(6,800);
  // f0->SetParameter(7,400);
  // f0->SetParameter(8,200);

  f0->SetParameter(9,1.0);
  f0->SetParameter(10,100);


  tmpstring = "(B0notB0bar == 1)";
  tmpstring += " && " + sigMassQuery;
  theTree->Draw("kstMass>>h1Dsig",tmpstring.c_str());

  tmpstring = "(B0notB0bar == 0)";
  tmpstring += " && " + sigMassQuery;
  theTree->Draw("kstBarMass>>h2Dsig",tmpstring.c_str());

  TH1D* h3Dsig = (TH1D*)h1Dsig->Clone("h3Dsig");
  h3Dsig->SetXTitle("M(#font[122]{K}^{*0}) (GeV)");
  h3Dsig->SetYTitle("Entries / (0.002 GeV)");
  h3Dsig->SetMarkerStyle(20);

  h3Dsig->Add(h2Dsig, 1.0);


  tmpstring = "(B0notB0bar == 1)";
  tmpstring += " && " + bkgMassQuery;
  theTree->Draw("kstMass>>h1Dbkg",tmpstring.c_str());

  tmpstring = "(B0notB0bar == 0)";
  tmpstring += " && " + bkgMassQuery;
  theTree->Draw("kstBarMass>>h2Dbkg",tmpstring.c_str());

  TH1D* h3Dbkg = (TH1D*)h1Dbkg->Clone("h3Dbkg");
  h3Dbkg->SetXTitle("M(#font[122]{K}^{*0}) (GeV)");
  h3Dbkg->SetYTitle("Entries / (0.002 GeV)");
  h3Dbkg->SetMarkerStyle(20);

  h3Dbkg->Add(h2Dbkg, 1.0);

  if (bkgSub == true)
    {
      h1Dsig->Add(h1Dbkg, -1.0);
      h2Dsig->Add(h2Dbkg, -1.0);
      h3Dsig->Add(h3Dbkg, -1.0);
    }


  c0->cd();
  h3Dsig->Fit("f0","R0");
  h3Dsig->Draw("e1");
  h3Dsig->GetFunction("f0")->Draw("same");

  fsig->SetParameter(0,f0->GetParameter(0));
  fsig->SetParameter(1,f0->GetParameter(1));
  fsig->SetParameter(2,f0->GetParameter(2));
  fsig->SetParameter(3,f0->GetParameter(3));
  fsig->SetParameter(4,f0->GetParameter(4));
  fsig->SetParameter(5,f0->GetParameter(5));
  fsig->SetParameter(6,f0->GetParameter(6));
  fsig->SetParameter(7,f0->GetParameter(7));
  fsig->SetParameter(8,f0->GetParameter(8));
  fsig->SetParameter(9,f0->GetParameter(9));
  fsig->Draw("same");

  fbkg->SetParameter(0,f0->GetParameter(10));
  fbkg->Draw("same");

  TLegend* leg0 = new TLegend(0.15, 0.75, 0.4, 0.89, "");
  leg0->AddEntry(h3Dsig,"Data");
  leg0->AddEntry(f0,"Total fit");
  leg0->AddEntry(fsig,"Signal");
  leg0->AddEntry(fbkg,"Background");
  leg0->SetFillColor(0);
  leg0->SetBorderSize(0);
  leg0->Draw();

  c0->Update();


  c1->cd();
  h1Dsig->Draw("e1p");
  h2Dsig->Draw("e1p same");

  TLegend* leg1 = new TLegend(0.15, 0.75, 0.4, 0.89, "");
  leg1->AddEntry(h1Dsig,"#font[122]{K}^{*0}");
  leg1->AddEntry(h2Dsig,"#font[122]{K}^{*0}_{bar}");
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->Draw();

  c1->Update();
}


void PlotKK (string fileName, bool bkgSub, string RECOorGEN)
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


  int nEntries;

  double KKminX      = 0.95;
  double KKmaxX      = 1.45;
  double KstminX     = 0.800;
  double KstmaxX     = 0.990;
  double Kstextra    = 0.015;
  double DalitzminX  = 0.65;
  double DalitzmaxX  = 0.95;
  double DalitzminY  = 14.5;
  double DalitzmaxY  = 27.0;
  unsigned int nBins = 100;
  double KKmass;
  double massPsiK;
  double massKpi;

  double signalSigma = sqrt(Utility->GetGenericParam("FRACMASSS") * Utility->GetGenericParam("SIGMAS1") * Utility->GetGenericParam("SIGMAS1") +
			    (1. - Utility->GetGenericParam("FRACMASSS")) * Utility->GetGenericParam("SIGMAS2") * Utility->GetGenericParam("SIGMAS2"));
  cout << "\n@@@ Signal sigma: " << signalSigma << " @@@" << endl;


  TFile* file0 = TFile::Open(fileName.c_str(),"READ");
  TTree* theTree = (TTree*)file0->Get("B0SingleCand/B0KstMuMuSingleCandNTuple");

  B0KstMuMuSingleCandTreeContent* NTuple = new B0KstMuMuSingleCandTreeContent();
  NTuple->Init();
  NTuple->ClearNTuple();
  NTuple->SetBranchAddresses(theTree);
  
  nEntries = theTree->GetEntries();
  cout << "\n@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;


  TCanvas* c0 = new TCanvas("c0","c0",10,10,700,500);
  TCanvas* c1 = new TCanvas("c1","c1",20,20,700,500);
  TCanvas* c2 = new TCanvas("c2","c2",30,30,700,500);

  TH1D* hKKSig = new TH1D("hKKSig","hKKSig",nBins,KKminX,KKmaxX);
  hKKSig->SetXTitle("M(#font[122]{K}^{#font[122]{+}}#font[122]{K}^{#font[122]{\55}}) (GeV)");
  hKKSig->SetYTitle("Entries / (0.005 GeV)");
  hKKSig->SetFillColor(kAzure+6);

  TH1D* hKKBkg = new TH1D("hKKBkg","hKKBkg",nBins,KKminX,KKmaxX);
  hKKBkg->SetXTitle("M(#font[122]{K}^{#font[122]{+}}#font[122]{K}^{#font[122]{\55}}) (GeV)");
  hKKBkg->SetYTitle("Entries / (0.005 GeV)");
  hKKBkg->SetFillColor(kAzure+6);

  TH1D* hKstSig = new TH1D("hKstSig","hKstSig",nBins,KstminX - Kstextra,KstmaxX + Kstextra);
  hKstSig->SetXTitle("M(K #pi) (GeV)");
  hKstSig->SetYTitle("Entries / (0.022 GeV)");
  hKstSig->SetMarkerStyle(20);

  TH1D* hKstBkg = new TH1D("hKstBkg","hKstBkg",nBins,KstminX - Kstextra,KstmaxX + Kstextra);
  hKstBkg->SetXTitle("M(K #pi) (GeV)");
  hKstBkg->SetYTitle("Entries / (0.022 GeV)");
  hKstSig->SetMarkerStyle(20);

  TH2D* hDalitzSig = new TH2D("hDalitzSig","hDalitzSig",nBins,DalitzminX,DalitzmaxX,nBins,DalitzminY,DalitzmaxY);
  hDalitzSig->SetXTitle("M^{2}(K #pi) (GeV^{2})");
  hDalitzSig->SetYTitle("M^{2}(J/#psi K) (GeV^{2})");
  hDalitzSig->SetZTitle("Entries / (0.03x0.125 (GeV^{4}))");

  TH2D* hDalitzBkg = new TH2D("hDalitzBkg","hDalitzBkg",nBins,DalitzminX,DalitzmaxX,nBins,DalitzminY,DalitzmaxY);
  hDalitzBkg->SetXTitle("M^{2}(K #pi) (GeV^{2})");
  hDalitzBkg->SetYTitle("M^{2}(J/#psi K) (GeV^{2})");
  hDalitzBkg->SetZTitle("Entries / (0.03x0.125 (GeV^{4}))");

  for (int entry = 0; entry < nEntries; entry++)
    {
      theTree->GetEntry(entry);


      if (RECOorGEN == "RECO")
	{
	  // #################
	  // # Use RECO info #
	  // #################
	  KKmass = Utility->computeInvMass(NTuple->kstTrkmPx->at(0),NTuple->kstTrkmPy->at(0),NTuple->kstTrkmPz->at(0),Utility->kaonMass,
					   NTuple->kstTrkpPx->at(0),NTuple->kstTrkpPy->at(0),NTuple->kstTrkpPz->at(0),Utility->kaonMass);
	  
	  if (NTuple->B0notB0bar == true)
	    {
	      massPsiK = Utility->computeInvMass(NTuple->mumPx->at(0) + NTuple->mupPx->at(0),NTuple->mumPy->at(0) + NTuple->mupPy->at(0),NTuple->mumPz->at(0) + NTuple->mupPz->at(0),NTuple->mumuMass->at(0),
						 NTuple->kstTrkpPx->at(0),NTuple->kstTrkpPy->at(0),NTuple->kstTrkpPz->at(0),Utility->kaonMass);
	      massKpi = Utility->computeInvMass(NTuple->kstTrkmPx->at(0),NTuple->kstTrkmPy->at(0),NTuple->kstTrkmPz->at(0),Utility->pionMass,
						NTuple->kstTrkpPx->at(0),NTuple->kstTrkpPy->at(0),NTuple->kstTrkpPz->at(0),Utility->kaonMass);
	    }
	  else
	    {
	      massPsiK = Utility->computeInvMass(NTuple->mumPx->at(0) + NTuple->mupPx->at(0),NTuple->mumPy->at(0) + NTuple->mupPy->at(0),NTuple->mumPz->at(0) + NTuple->mupPz->at(0),NTuple->mumuMass->at(0),
						 NTuple->kstTrkmPx->at(0),NTuple->kstTrkmPy->at(0),NTuple->kstTrkmPz->at(0),Utility->kaonMass);
	      massKpi = Utility->computeInvMass(NTuple->kstTrkmPx->at(0),NTuple->kstTrkmPy->at(0),NTuple->kstTrkmPz->at(0),Utility->kaonMass,
						NTuple->kstTrkpPx->at(0),NTuple->kstTrkpPy->at(0),NTuple->kstTrkpPz->at(0),Utility->pionMass);
	    }
	}
      else if (RECOorGEN == "GEN")
	{
	  // ################
	  // # Use GEN info #
	  // ################
	  KKmass = Utility->computeInvMass(NTuple->genKstTrkmPx,NTuple->genKstTrkmPy,NTuple->genKstTrkmPz,Utility->kaonMass,
					   NTuple->genKstTrkpPx,NTuple->genKstTrkpPy,NTuple->genKstTrkpPz,Utility->kaonMass);
	      
	  if (NTuple->B0notB0bar == true)
	    {
	      massPsiK = Utility->computeInvMass(NTuple->genMumPx + NTuple->genMupPx,NTuple->genMumPy + NTuple->genMupPy,NTuple->genMumPz + NTuple->genMupPz,NTuple->mumuMass->at(0),
						 NTuple->genKstTrkpPx,NTuple->genKstTrkpPy,NTuple->genKstTrkpPz,Utility->kaonMass);
	      massKpi = Utility->computeInvMass(NTuple->genKstTrkmPx,NTuple->genKstTrkmPy,NTuple->genKstTrkmPz,Utility->pionMass,
						NTuple->genKstTrkpPx,NTuple->genKstTrkpPy,NTuple->genKstTrkpPz,Utility->kaonMass);
	    }
	  else
	    {
	      massPsiK = Utility->computeInvMass(NTuple->genMumPx + NTuple->genMupPx,NTuple->genMumPy + NTuple->genMupPy,NTuple->genMumPz + NTuple->genMupPz,NTuple->mumuMass->at(0),
						 NTuple->genKstTrkmPx,NTuple->genKstTrkmPy,NTuple->genKstTrkmPz,Utility->kaonMass);
	      massKpi = Utility->computeInvMass(NTuple->genKstTrkmPx,NTuple->genKstTrkmPy,NTuple->genKstTrkmPz,Utility->kaonMass,
						NTuple->genKstTrkpPx,NTuple->genKstTrkpPy,NTuple->genKstTrkpPz,Utility->pionMass);
	    }
	}
      else
	{
	  cout << "[MakePlots::PlotKK]\tWrong parameter: " << RECOorGEN << endl;
	  exit (1);
	}


      if ((bkgSub == false) || (fabs(NTuple->B0MassArb - Utility->B0Mass) < Utility->GetGenericParam("NSigmaB0S")*signalSigma))
	{
	  // ####################
	  // # Make signal plot #
	  // ####################
	  hKKSig->Fill(KKmass);
	  hKstSig->Fill(massKpi);
	  hDalitzSig->Fill(massKpi*massKpi,massPsiK*massPsiK);
	}
      else if (((NTuple->B0MassArb > Utility->B0Mass + Utility->GetGenericParam("NSigmaB0B")*signalSigma) &&
		(NTuple->B0MassArb < Utility->B0Mass + (Utility->GetGenericParam("NSigmaB0B")+Utility->GetGenericParam("NSigmaB0S"))*signalSigma)) ||
	       ((NTuple->B0MassArb < Utility->B0Mass - Utility->GetGenericParam("NSigmaB0B")*signalSigma) &&
		(NTuple->B0MassArb > Utility->B0Mass - (Utility->GetGenericParam("NSigmaB0B")+Utility->GetGenericParam("NSigmaB0S"))*signalSigma)))
	{
	  // ########################
	  // # Make background plot #
	  // ########################
	  hKKBkg->Fill(KKmass);
	  hKstBkg->Fill(massKpi);
	  hDalitzBkg->Fill(massKpi*massKpi,massPsiK*massPsiK);
	}
    }


  if (bkgSub == true) hKKSig->Add(hKKBkg, -1.0);

  c0->cd();
  hKKSig->Draw();
  c0->Update();


  if (bkgSub == true) hKstSig->Add(hKstBkg, -1.0);

  c1->cd();
  hKstSig->Draw("e1p");
  c1->Update();


  if (bkgSub == true) hDalitzSig->Add(hDalitzBkg, -1.0);

  c2->cd();
  hDalitzSig->Draw("cont4z");
  c2->Update();
}


void PlotMuHadMass (string fileName)
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


  double minX = 0.0;
  double maxX = 5.0;
  double MuHadMass;
  int nEntries;
  unsigned int nBins = 100;


  TCanvas* c0 = new TCanvas("c0","c0",10,10,700,500);
  TCanvas* c1 = new TCanvas("c1","c1",20,20,700,500);


  TFile* file0 = TFile::Open(fileName.c_str(),"READ");
  TTree* theTree = (TTree*)file0->Get("B0SingleCand/B0KstMuMuSingleCandNTuple");

  B0KstMuMuSingleCandTreeContent* NTuple = new B0KstMuMuSingleCandTreeContent();
  NTuple->Init();
  NTuple->ClearNTuple();
  NTuple->SetBranchAddresses(theTree);
 
  nEntries = theTree->GetEntries();
  cout << "\n@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;

  TH1D* histoK = new TH1D("histoK","histoK",nBins,minX,maxX);
  histoK->SetXTitle("M(K #mu) (GeV)");
  histoK->SetYTitle("Entries / (0.056 GeV)");
  histoK->SetFillColor(kAzure+6);

  TH1D* histoPi = new TH1D("histoPi","histoPi",nBins,minX,maxX);
  histoPi->SetXTitle("M(#pi #mu) (GeV)");
  histoPi->SetYTitle("Entries / (0.056 GeV)");
  histoPi->SetFillColor(kAzure+6);

  for (int entry = 0; entry < nEntries; entry++)
    {
     theTree->GetEntry(entry);

     if ((NTuple->B0MassArb > Utility->B0Mass - Utility->GetGenericParam("B0MassIntervalLeft")) && (NTuple->B0MassArb < Utility->B0Mass + Utility->GetGenericParam("B0MassIntervalRight")))
       {
	 if (NTuple->B0notB0bar == true)
	   {
	     MuHadMass = Utility->computeInvMass(NTuple->mumPx->at(0),NTuple->mumPy->at(0),NTuple->mumPz->at(0),Utility->muonMass,
						 NTuple->kstTrkpPx->at(0),NTuple->kstTrkpPy->at(0),NTuple->kstTrkpPz->at(0),Utility->kaonMass);
	 
	     histoK->Fill(MuHadMass);
	 
	     MuHadMass = Utility->computeInvMass(NTuple->mupPx->at(0),NTuple->mupPy->at(0),NTuple->mupPz->at(0),Utility->muonMass,
						 NTuple->kstTrkmPx->at(0),NTuple->kstTrkmPy->at(0),NTuple->kstTrkmPz->at(0),Utility->pionMass);
	 
	     histoPi->Fill(MuHadMass);
	   }
	 else
	   {
	     MuHadMass = Utility->computeInvMass(NTuple->mumPx->at(0),NTuple->mumPy->at(0),NTuple->mumPz->at(0),Utility->muonMass,
						 NTuple->kstTrkpPx->at(0),NTuple->kstTrkpPy->at(0),NTuple->kstTrkpPz->at(0),Utility->pionMass);
	 
	     histoPi->Fill(MuHadMass);
	 
	     MuHadMass = Utility->computeInvMass(NTuple->mupPx->at(0),NTuple->mupPy->at(0),NTuple->mupPz->at(0),Utility->muonMass,
						 NTuple->kstTrkmPx->at(0),NTuple->kstTrkmPy->at(0),NTuple->kstTrkmPz->at(0),Utility->kaonMass);
	 
	     histoK->Fill(MuHadMass);
	   }
       }
    }

  c0->cd();
  histoK->Draw();
  c0->Update();

  c1->cd();
  histoPi->Draw();
  c1->Update();
}


void MakeupNLLandPvalPlots (string fileName, int specBin, string PlotType)
// #################################################
// # If specBin = -1  --> read difference plot     #
// # If specBin != -1 --> read NLL and PULLs plots #
// #################################################
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

  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);


  unsigned int it = 1;
  stringstream myString;
  string car;
  double val;
  TFile* fileID;
  TCanvas *c0, *c1;
  TPad* pNLL;
  TH1D* histoNLL;
  TPad* pPULL;
  TH1D* histoPULL;
  TPad* pDIFF1;
  TH1D* histoDIFF1;
  TPad* pDIFF2;
  TH1D* histoDIFF2;
  vector<vector<double>*> vecNLL;
  Utility->ReadNLL(ParameterFILE,&vecNLL);
  if      (PlotType == "FL")   { val = vecNLL[0]->operator[](specBin); car = "1"; }
  else if (PlotType == "AFB")  { val = vecNLL[1]->operator[](specBin); car = "2"; }
  else if (PlotType == "AT2")  { val = vecNLL[2]->operator[](specBin); car = "1"; }
  else if (PlotType == "ATIM") { val = vecNLL[3]->operator[](specBin); car = "1"; }
  else if (PlotType == "BF")   { val = vecNLL[4]->operator[](specBin); car = "1"; }
  else
    {
      cout << "[MakePlots::MakeupNLLandPvalPlots]\tWrong parameter: " << PlotType.c_str() << endl;
      exit (1);
    }


  if (specBin != -1)
    {
      // ################################
      // # Read NLL and pulls from file #
      // ################################

      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "@ Appended ID to histogram names: " << car.c_str() << " @" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      
      fileID = TFile::Open(fileName.c_str(),"READ");

      myString.clear(); myString.str("");
      myString << "cNLL" << car << ";" << it;
      c0 = (TCanvas*)fileID->Get(myString.str().c_str());

      myString.clear(); myString.str("");
      myString << "cNLL" << car << "_1";
      pNLL = (TPad*)c0->GetPrimitive(myString.str().c_str());

      myString.clear(); myString.str("");
      myString << "histoNLL" << car;
      histoNLL = (TH1D*)pNLL->GetPrimitive(myString.str().c_str());
      histoNLL->SetFillColor(kGreen-7);

      myString.clear(); myString.str("");
      myString << "cNLL" << car << "_2";
      pPULL = (TPad*)c0->GetPrimitive(myString.str().c_str());

      myString.clear(); myString.str("");
      myString << "histoPull" << car;
      histoPULL = (TH1D*)pPULL->GetPrimitive(myString.str().c_str());
      histoPULL->SetFillColor(kAzure+6);

      c0->Draw();
      c0->Update();

      c1 = NULL;
      it++;
      myString.clear(); myString.str("");
      myString << "cNLL" << car << ";" << it;
      c1 = (TCanvas*)fileID->Get(myString.str().c_str());
      while (c1 != NULL)
	{
	  cout << "Adding histogram #" << it << endl;

	  myString.clear(); myString.str("");
	  myString << "cNLL" << car << "_1";
	  pNLL = (TPad*)c1->GetPrimitive(myString.str().c_str());

	  myString.clear(); myString.str("");
	  myString << "histoNLL" << car;
	  histoNLL->Add((TH1D*)pNLL->GetPrimitive(myString.str().c_str()));

	  myString.clear(); myString.str("");
	  myString << "cNLL" << car << "_2";
	  pPULL = (TPad*)c1->GetPrimitive(myString.str().c_str());

	  myString.clear(); myString.str("");
	  myString << "histoPull" << car;
	  histoPULL->Add((TH1D*)pPULL->GetPrimitive(myString.str().c_str()));

	  c1 = NULL;
	  it++;
	  myString.clear(); myString.str("");
	  myString << "cNLL" << car << ";" << it;
	  c1 = (TCanvas*)fileID->Get(myString.str().c_str());
	}

      c0->cd(1);
      histoNLL->Draw();
      DrawExclusion(val,histoNLL->GetBinLowEdge(histoNLL->GetNbinsX())+histoNLL->GetBinWidth(1),0.0,histoNLL->GetMaximum()*1.1,"NLL",3001,kGray);

      c0->cd(2);
      histoPULL->Draw();
      histoPULL->Fit("gaus","0");
      histoPULL->GetFunction("gaus")->Draw("same");
    }
  else
    {
      // ##############################
      // # Read differences from file #
      // ##############################
      fileID = TFile::Open(fileName.c_str(),"READ");
      c0 = (TCanvas*)fileID->Get("cNLL3;1");
       if (PlotType != "AFB")
	 {
	   cout << "Adding histogram #" << it << endl;

	   histoDIFF1 = (TH1D*)c0->GetPrimitive("histoDiff1");
	   histoDIFF1->SetFillColor(kAzure+6);

	   c0->Draw();
	   c0->Update();

	   c1 = NULL;
	   it++;
	   myString.clear(); myString.str("");
	   myString << "cNLL3;" << it;
	   c1 = (TCanvas*)fileID->Get(myString.str().c_str());
	   while (c1 != NULL)
	     {
	       histoDIFF1->Add((TH1D*)c1->GetPrimitive("histoDiff1"));

	       c1 = NULL;
	       it++;
	       myString.clear(); myString.str("");
	       myString << "cNLL3;" << it;
	       c1 = (TCanvas*)fileID->Get(myString.str().c_str());
	     }

	   c0->cd();
	   histoDIFF1->Draw();
 	   histoDIFF1->Fit("gaus","0");
	   histoDIFF1->GetFunction("gaus")->Draw("same");
	 }
       else
	 {
	   cout << "Adding histogram #" << it << endl;

	   pDIFF1 = (TPad*)c0->GetPrimitive("cNLL3_1");
	   histoDIFF1 = (TH1D*)pDIFF1->GetPrimitive("histoDiff1");
	   histoDIFF1->SetFillColor(kAzure+6);

	   pDIFF2 = (TPad*)c0->GetPrimitive("cNLL3_2");
	   histoDIFF2 = (TH1D*)pDIFF2->GetPrimitive("histoDiff2");
	   histoDIFF2->SetFillColor(kAzure+6);

	   c0->Draw();
	   c0->Update();

	   c1 = NULL;
	   it++;
	   myString.clear(); myString.str("");
	   myString << "cNLL3;" << it;
	   c1 = (TCanvas*)fileID->Get(myString.str().c_str());
	   while (c1 != NULL)
	     {
	       pDIFF1 = (TPad*)c1->GetPrimitive("cNLL3_1");
	       histoDIFF1->Add((TH1D*)pDIFF1->GetPrimitive("histoDiff1"));

	       pDIFF2 = (TPad*)c1->GetPrimitive("cNLL3_2");
	       histoDIFF2->Add((TH1D*)pDIFF2->GetPrimitive("histoDiff2"));

	       c1 = NULL;
	       it++;
	       myString.clear(); myString.str("");
	       myString << "cNLL3;" << it;
	       c1 = (TCanvas*)fileID->Get(myString.str().c_str());
	     }

	   c0->cd(1);
	   histoDIFF1->Draw();
	   histoDIFF1->Fit("gaus","0");
	   histoDIFF1->GetFunction("gaus")->Draw("same");

	   c0->cd(2);
	   histoDIFF2->Draw();
	   histoDIFF2->Fit("gaus","0");
	   histoDIFF2->GetFunction("gaus")->Draw("same");
	 }
    }

  c0->Update();
}


void MakePvaluePlot (string toyMCfileName, int specBin, string PlotType)
// ################################################
// # If specBin == -1 then loop over all q^2 bins #
// ################################################
// # PlotType: #
// #############
// # FL        #
// # AFB       #
// # AT2       #
// # ATIM      #
// # BF        #
// #############
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
  gStyle->SetTitleOffset(1.25,"y"); 
  TGaxis::SetMaxDigits(3);


  unsigned int it = 1;
  stringstream myString;
  string car;
  double val;
  TFile* fileToy;
  TCanvas *c0, *c1;
  TPad* pNLL;
  TH1D* histoNLL;

  vector<double> q2Bins;
  Utility->Readq2Bins(ParameterFILE,&q2Bins);
  double* q2Bins_ = Utility->MakeBinning(&q2Bins);

  vector<vector<double>*> vecNLL;
  Utility->ReadNLL(ParameterFILE,&vecNLL);


  TCanvas* c2 = new TCanvas("c2","c2",10,10,700,500);
  c2->cd();
  TH1D* pval = new TH1D("pval","pval",q2Bins.size()-1,q2Bins_);
  pval->SetMarkerStyle(20);
  pval->SetXTitle("q^{2} (GeV^{2})");
  pval->SetYTitle("p-value");
  pval->SetMinimum(0.0);
  pval->SetMaximum(1.05);

  cout << "\n@@@ Computaion of the p-value from profile likelihood for NLL from paramter file: " << ParameterFILE << " @@@" << endl;
  toyMCfileName.replace(toyMCfileName.find("_NLL.root")-1,10,"");

  for (int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? (int)(q2Bins.size()-1) : specBin+1); i++)
    {
      if ((PlotType == "BF") && ((i == Utility->GetJPsiBin(&q2Bins)) || (i == Utility->GetPsiPBin(&q2Bins)))) continue;

      myString.clear();
      myString.str("");
      myString << toyMCfileName << i << "_NLL.root";
      cout << "\nReading NLL distribution from file: " << myString.str().c_str() << endl;
      fileToy = new TFile(myString.str().c_str(),"READ");
      if (fileToy->IsZombie() == false)
	{
	  if      (PlotType == "FL")   { val = vecNLL[0]->operator[](i); car = "1"; }
	  else if (PlotType == "AFB")  { val = vecNLL[1]->operator[](i); car = "2"; }
	  else if (PlotType == "AT2")  { val = vecNLL[2]->operator[](i); car = "1"; }
	  else if (PlotType == "ATIM") { val = vecNLL[3]->operator[](i); car = "1"; }
	  else if (PlotType == "BF")   { val = vecNLL[4]->operator[](i); car = "1"; }
	  else
	    {
	      cout << "[MakePlots::MakePvaluePlot]\tWrong parameter: " << PlotType.c_str() << endl;
	      exit (1);
	    }

	  myString.clear(); myString.str("");
	  myString << "cNLL" << car;
	  c0 = (TCanvas*)fileToy->Get(myString.str().c_str());

	  myString.clear(); myString.str("");
	  myString << "cNLL" << car << "_1";
	  pNLL = (TPad*)c0->GetPrimitive(myString.str().c_str());

	  myString.clear(); myString.str("");
	  myString << "histoNLL" << car;
	  histoNLL = (TH1D*)pNLL->GetPrimitive(myString.str().c_str());

	  c1 = NULL;
	  it++;
	  myString.clear(); myString.str("");
	  myString << "cNLL" << car << ";" << it;
	  c1 = (TCanvas*)fileToy->Get(myString.str().c_str());
	  while (c1 != NULL)
	    {
	      cout << "Adding histogram #" << it << endl;
	      
	      myString.clear(); myString.str("");
	      myString << "cNLL" << car << "_1";
	      pNLL = (TPad*)c1->GetPrimitive(myString.str().c_str());
	      
	      myString.clear(); myString.str("");
	      myString << "histoNLL" << car;
	      histoNLL->Add((TH1D*)pNLL->GetPrimitive(myString.str().c_str()));

	      c1 = NULL;
	      it++;
	      myString.clear(); myString.str("");
	      myString << "cNLL" << car << ";" << it;
	      c1 = (TCanvas*)fileToy->Get(myString.str().c_str());
	    }

	  histoNLL->Scale(1./histoNLL->Integral());
	  pval->SetBinContent(i+1,(isnan(histoNLL->Integral(1,histoNLL->FindBin(val))) == false ? histoNLL->Integral(histoNLL->FindBin(val),histoNLL->GetNbinsX()) : 0.0));
	  pval->SetBinError(i+1,pval->GetBinContent(i+1)*1e-3);
	  cout << "p-value for " << PlotType << " q^2 bin #" << i << " --> " << pval->GetBinContent(i+1) << endl;

	  fileToy->Close();
	}
    }
  
  c2->cd();
  pval->Draw("e1p");
  DrawExclusion(pval->GetBinLowEdge(1),pval->GetBinLowEdge(pval->GetNbinsX())+pval->GetBinWidth(pval->GetNbinsX()),0.0,0.05,"p-val",3001,kRed-7);
  c2->Update();
}


void getHfromToy (string fileNameIn, string fileNameOut, unsigned int intVal)
{
  stringstream myString;

  TFile* fIn = new TFile(fileNameIn.c_str(),"READ");
  myString.clear(); myString.str("");
  myString << "cNLL" << intVal;
  TCanvas* c0 = (TCanvas*)fIn->Get(myString.str().c_str());


  myString.clear(); myString.str("");
  myString << "cNLL" << intVal << "_1";
  TPad* p0 = (TPad*)c0->GetPrimitive(myString.str().c_str());

  myString.clear(); myString.str("");
  if      (intVal == 1) myString << "histoNLL1";
  else if (intVal == 2) myString << "histoNLL2";
  else if (intVal == 3) myString << "histoDiff1";
  else
    {
      cout << "[MakePlots::getHfromToy]\tWrong parameter: " << intVal << endl;
      exit (1);
    }
  TH1D* h0 = (TH1D*)p0->GetPrimitive(myString.str().c_str());


  myString.clear(); myString.str("");
  myString << "cNLL" << intVal << "_2";
  TPad* p1 = (TPad*)c0->GetPrimitive(myString.str().c_str());

  myString.clear(); myString.str("");
  if      (intVal == 1) myString << "histoPull1";
  else if (intVal == 2) myString << "histoPull2";
  else if (intVal == 3) myString << "histoDiff2";
  else
    {
      cout << "[MakePlots::getHfromToy]\tWrong parameter: " << intVal << endl;
      exit (1);
    }
  TH1D* h1 = (TH1D*)p1->GetPrimitive(myString.str().c_str());


  TFile* fOut = new TFile(fileNameOut.c_str(),"RECREATE");
  fOut->cd();
  h0->Write();
  h1->Write();
  fOut->Close();
  fIn->Close();

  delete fIn;
  delete fOut;

  cout << "@@@ Histogram extraction done @@@" << endl;
}


int main (int argc, char** argv)
{
  if (argc >= 2)
    {
      string option = argv[1];
      unsigned int intVal = 1;
      unsigned int q2BinIndx = 0;
      double realVal1 = 0;
      double realVal2 = 0;
      string fileName1;
      string fileName2;
      string fileName3;
      string tmpStr;
      if ((option == "EvalMultyRun") && (argc >= 4))
	{
	  intVal = atoi(argv[2]);
	  q2BinIndx = atoi(argv[3]);
	  if (argc >= 5) fileName1 = argv[4];
	  else           fileName1 = FitSysFILE;
	  if (argc >= 6) realVal1 = atof(argv[5]);
	  else           realVal1 = 0.0;
	  if (argc == 7) realVal2 = atof(argv[6]);
	  else           realVal2 = 0.0;
	}
      else if (((option == "Phy") || (option == "DataMC")) && (argc == 3)) intVal = atoi(argv[2]);
      else if (((option == "Pval") || (option == "makeupNLL")) && (argc == 5))
	{
	  fileName1 = argv[2];
	  q2BinIndx = atoi(argv[3]);
	  tmpStr    = argv[4];
	}
      else if (((option == "MuMuMass") || (option == "KstMass")) && (argc == 4))
	{
	  fileName1 = argv[2];
	  intVal    = atoi(argv[3]);
	}
      else if ((option == "KKMass") && (argc == 5))
	{
	  fileName1 = argv[2];
	  intVal    = atoi(argv[3]);
	  tmpStr    = argv[4];
	}
      else if ((option == "getHfromToy") && (argc == 5))
	{
	  fileName1 = argv[2];
	  fileName2 = argv[3];
	  intVal    = atoi(argv[4]);
	}
      else if ((option == "MuHadMass") && (argc == 3)) fileName1 = argv[2];
      else if (option != "PhyRegion")
	{
	  cout << "Parameter missing:" << endl;
	  cout << "./MakePlots [Phy EvalMultyRun DataMC PhyRegion Pval makeupNLL MuMuMass KKMass KstMass MuHadMass getHfromToy]" << endl;
	  cout << "[Phy:0-2||10-12]" << endl;
	  cout << "[EvalMultyRun: 0-3||10-12 q^2 bin index 0-7 [fileName] [NLL interval] [NLL less than]" << endl;
	  cout << "[DataMC: 0-27]" << endl;
	  cout << "[Pval or makeupNLL: ToyMCfileName q^2 bin index 0-7 PlotType]" << endl;
	  cout << "[MuMuMass or KstMass: dataFileName bkgSub]" << endl;
	  cout << "[KKMass: dataFileName bkgSub RECOorGEN]" << endl;
	  cout << "[MuHadMass: dataFileName]" << endl;
	  cout << "[getHfromToy: ToyMCfileName OutputFileName canvas# [1-3]]" << endl;

	  return 1;
	}


      cout << "\n@@@ Settings @@@" << endl;
      cout << "ParameterFILE : "       << ParameterFILE << endl;
      cout << "ParameterFILE_MCGEN: "  << ParameterFILE_MCGEN << endl;
      cout << "ParameterFILE_MCRECO: " << ParameterFILE_MCRECO << endl;

      cout << "\nSMFL: "   << SMFL << endl;
      cout << "SMAFB: "    << SMAFB << endl;
      cout << "SMBF: "     << SMBF << endl;
      cout << "SMBINFL: "  << SMBF << endl;
      cout << "SMBINAFB: " << SMBINAFB << endl;
      cout << "SMBINBF: "  << SMBINBF << endl;

      cout << "\nSingleCand_MCkstJPsi: " << SingleCand_MCkstJPsi << endl;
      cout << "SingleCand_MCkstPsi2S: "  << SingleCand_MCkstPsi2S << endl;
      cout << "SingleCand_Data: "        << SingleCand_Data << endl;

      cout << "\nFitSysFILE: "        << FitSysFILE << endl;
      cout << "OutHisto: "            << OutHisto << endl;
      cout << "YvalueOutsideLimits: " << YvalueOutsideLimits << endl;
      cout << "FORPAPER: "            << FORPAPER << endl;

      cout << "\noption: "  << option << endl;
      cout << "intVal: "    << intVal << endl;
      cout << "q2BinIndx: " << q2BinIndx << endl;
      cout << "realVal1: "  << realVal1 << endl;
      cout << "realVal2: "  << realVal2 << endl;
      cout << "fileName1: " << fileName1 << endl;
      cout << "fileName2: " << fileName2 << endl;
      cout << "fileName3: " << fileName3 << endl;
      cout << "tmpStr: "    << tmpStr << endl; 

      if (option == "getHfromToy") gROOT->SetBatch(true);
      TApplication theApp ("Applications", &argc, argv);
	

      Utility = new Utils();
      Utility->ReadGenericParam(ParameterFILE);

      if      (option == "Phy")          MakePhysicsPlots(intVal);
      else if (option == "EvalMultyRun") EvalMultyRun(intVal,q2BinIndx,fileName1,realVal1,realVal2);
      else if (option == "DataMC")       MakeComparisonDataMC(intVal);
      else if (option == "PhyRegion")    CheckPhysicsRegion();
      else if (option == "Pval")         MakePvaluePlot(fileName1,q2BinIndx,tmpStr);
      else if (option == "MuMuMass")     PlotMuMu(fileName1,intVal);
      else if (option == "KKMass")       PlotKK(fileName1,intVal,tmpStr);
      else if (option == "KstMass")      PlotKst(fileName1,intVal);
      else if (option == "MuHadMass")    PlotMuHadMass(fileName1);
      else if (option == "makeupNLL")    MakeupNLLandPvalPlots(fileName1,q2BinIndx,tmpStr);
      else if (option == "getHfromToy")  { getHfromToy(fileName1,fileName2,intVal); exit(0); }
      else
	{
	  cout << "./MakePlots [Phy EvalMultyRun DataMC PhyRegion Pval makeupNLL MuMuMass KKMass KstMass MuHadMass getHfromToy]" << endl;
	  cout << "[Phy:0-2||10-12]" << endl;
	  cout << "[EvalMultyRun: 0-3||10-12 [q^2 bin index 0-7] [fileName] [NLL interval] [NLL less than]" << endl;
	  cout << "[DataMC: 0-27]" << endl;
	  cout << "[Pval or makeupNLL: ToyMCfileName q^2 bin index 0-7 PlotType]" << endl;
	  cout << "[MuMuMass or KstMass: dataFileName bkgSub]" << endl;
	  cout << "[KKMass: dataFileName bkgSub RECOorGEN]" << endl;
	  cout << "[MuHadMass: dataFileName]" << endl;
	  cout << "[getHfromToy: ToyMCfileName OutputFileName canvas# [1-3]]" << endl;

	  cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
	  cout << "For [MuMuMass KKMass KstMass]:" << endl;
	  cout << "bkgSub = 0 (do not subtract background)" << endl;
	  cout << "bkgSub = 1 (subtract background)" << endl;

	  cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
	  cout << "For [KKMass]:" << endl;
	  cout << "RECOorGEN = RECO (use RECO information)" << endl;
	  cout << "RECOorGEN = GEN (use GEN information)" << endl;
      
	  cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
	  cout << "For [Phy]:" << endl;
	  cout << "0 = Fl GEN vs RECO-MC" << endl;
	  cout << "1 = Afb GEN vs RECO-MC" << endl;
	  cout << "2 = BF GEN vs RECO-MC" << endl;
	  cout << "10 = Fl Data vs Theory" << endl;
	  cout << "11 = Afb Data vs Theory" << endl;
	  cout << "12 = BF Data vs Theory" << endl;

	  cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
	  cout << "For [EvalMultyRun]:" << endl;
	  cout << "0 = Fl" << endl;
	  cout << "1 = Afb" << endl;
	  cout << "2 = BF" << endl;
	  cout << "3 = N.A." << endl;
	  cout << "10 = Fl multy minima" << endl;
	  cout << "11 = Afb multy minima" << endl;
	  cout << "12 = BF multy minima" << endl;

	  cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
	  cout << "For [DataMC]:" << endl;
	  cout << "0 = B0 pT" << endl;
	  cout << "1 = B0 eta" << endl;

	  cout << "2 = mu+ pT" << endl;
	  cout << "3 = mu- pT" << endl;
	  cout << "4 = mu+ eta" << endl;
	  cout << "5 = mu- eta" << endl;
	  cout << "6 = mu+ phi, eta range" << endl;
	  cout << "7 = mu+ phi, eta range" << endl;
	  cout << "8 = mu+ phi, eta range" << endl;
	  cout << "9 = mu+ phi, eta range" << endl;
	  cout << "10 = mu- phi, eta range" << endl;
	  cout << "11 = mu- phi, eta range" << endl;
	  cout << "12 = mu- phi, eta range" << endl;
	  cout << "13 = mu- phi, eta range" << endl;
	  
	  cout << "14 = K*0 trk+ pT" << endl;
	  cout << "15 = K*0 trk- pT" << endl;
	  cout << "16 = K*0 trk+ eta" << endl;
	  cout << "17 = K*0 trk- eta" << endl;
	  cout << "18 = K*0 trk+ phi, eta range" << endl;
	  cout << "19 = K*0 trk+ phi, eta range" << endl;
	  cout << "20 = K*0 trk+ phi, eta range" << endl;
	  cout << "21 = K*0 trk+ phi, eta range" << endl;
	  cout << "22 = K*0 trk- phi, eta range" << endl;
	  cout << "23 = K*0 trk- phi, eta range" << endl;
	  cout << "24 = K*0 trk- phi, eta range" << endl;
	  cout << "25 = K*0 trk- phi, eta range" << endl;

	  cout << "26 = cos(theta_K)" << endl;
	  cout << "27 = cos(theta_l)" << endl;

	  cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
	  cout << "For [Pval makeupNLL]:" << endl;
	  cout << "- FL" << endl;
	  cout << "- AFB" << endl;
	  cout << "- AT2" << endl;
	  cout << "- ATIM" << endl;
	  cout << "- BF" << endl;

	  return 1;
	}

      delete Utility;
      theApp.Run (); // Eventloop on air
      return 0;
    }
  else
    {
      cout << "./MakePlots [Phy EvalMultyRun DataMC PhyRegion Pval makeupNLL MuMuMass KKMass KstMass MuHadMass getHfromToy]" << endl;
      cout << "[Phy:0-2||10-12]" << endl;
      cout << "[EvalMultyRun: 0-3||10-12 [q^2 bin index 0-7] [fileName] [NLL interval] [NLL less than]" << endl;
      cout << "[DataMC: 0-27]" << endl;
      cout << "[Pval or makeupNLL: ToyMCfileName q^2 bin index 0-7 PlotType]" << endl;
      cout << "[MuMuMass or KstMass: dataFileName bkgSub]" << endl;
      cout << "[KKMass: dataFileName bkgSub RECOorGEN]" << endl;
      cout << "[MuHadMass: dataFileName]" << endl;
      cout << "[getHfromToy: ToyMCfileName OutputFileName canvas# [1-3]]" << endl;

      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "For [MuMuMass KKMass KstMass]:" << endl;
      cout << "bkgSub = 0 (do not subtract background)" << endl;
      cout << "bkgSub = 1 (subtract background)" << endl;

      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "For [KKMass]:" << endl;
      cout << "RECOorGEN = RECO (use RECO information)" << endl;
      cout << "RECOorGEN = GEN (use GEN information)" << endl;

      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "For [Phy]:" << endl;
      cout << "0 = Fl GEN vs RECO-MC" << endl;
      cout << "1 = Afb GEN vs RECO-MC" << endl;
      cout << "2 = BF GEN vs RECO-MC" << endl;
      cout << "10 = Fl Data vs Theory" << endl;
      cout << "11 = Afb Data vs Theory" << endl;
      cout << "12 = BF Data vs Theory" << endl;

      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "For [EvalMultyRun]:" << endl;
      cout << "0 = Fl" << endl;
      cout << "1 = Afb" << endl;
      cout << "2 = BF" << endl;
      cout << "3 = N.A." << endl;
      cout << "10 = Fl multy minima" << endl;
      cout << "11 = Afb multy minima" << endl;
      cout << "12 = BF multy minima" << endl;

      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "For [DataMC]:" << endl;
      cout << "0 = B0 pT" << endl;
      cout << "1 = B0 eta" << endl;

      cout << "2 = mu+ pT" << endl;
      cout << "3 = mu- pT" << endl;
      cout << "4 = mu+ eta" << endl;
      cout << "5 = mu- eta" << endl;
      cout << "6 = mu+ phi, eta range" << endl;
      cout << "7 = mu+ phi, eta range" << endl;
      cout << "8 = mu+ phi, eta range" << endl;
      cout << "9 = mu+ phi, eta range" << endl;
      cout << "10 = mu- phi, eta range" << endl;
      cout << "11 = mu- phi, eta range" << endl;
      cout << "12 = mu- phi, eta range" << endl;
      cout << "13 = mu- phi, eta range" << endl;

      cout << "14 = K*0 trk+ pT" << endl;
      cout << "15 = K*0 trk- pT" << endl;
      cout << "16 = K*0 trk+ eta" << endl;
      cout << "17 = K*0 trk- eta" << endl;
      cout << "18 = K*0 trk+ phi, eta range" << endl;
      cout << "19 = K*0 trk+ phi, eta range" << endl;
      cout << "20 = K*0 trk+ phi, eta range" << endl;
      cout << "21 = K*0 trk+ phi, eta range" << endl;
      cout << "22 = K*0 trk- phi, eta range" << endl;
      cout << "23 = K*0 trk- phi, eta range" << endl;
      cout << "24 = K*0 trk- phi, eta range" << endl;
      cout << "25 = K*0 trk- phi, eta range" << endl;

      cout << "26 = cos(theta_K)" << endl;
      cout << "27 = cos(theta_l)" << endl;

      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "For [Pval makeupNLL]:" << endl;
      cout << "- FL" << endl;
      cout << "- AFB" << endl;
      cout << "- AT2" << endl;
      cout << "- ATIM" << endl;
      cout << "- BF" << endl;

      return 1;
    }
  
  return 0;
}

// ############################################################################
// # Program to fit the B0 invariant mass for the B0 --> K*0 mu+ mu- analysis #
// ############################################################################
// # Measure forward-backward asymmetry (Afb) from dimuons                    #
// # Longitudinal polarization (Fl) from K*0 (K pi)                           #
// # Branching fraction with respect to B0 --> K*0 J/psi (mu+ mu-)            #
// # Branching fraction with respect to B0 --> K*0 psi(2S) (mu+ mu-)          #
// ############################################################################
// # Author: Mauro Dinardo                                                    #
// ############################################################################
// # Open tasks:                                                              #
// # - S-P wave p.d.f. is not is not implemented for option "5"               #
// # - 4D fit to B0 mass, theta_K, theta_l, phi is not implemented            #
// # Search for @TMP@ to look for temporary code options                      #
// ############################################################################

// ############ To make simultaneus fit to multiple HLT categories ############
// RooDataSet* SingleCandNTuple;
// RooDataSet* SingleCandNTuple_HLT1 = NULL;
// RooDataSet* SingleCandNTuple_HLT2 = NULL;
// HLTcategory = new RooCategory("HLTcategory","HLTcategory");
// HLTcategory->defineType("HLT1");
// HLTcategory->defineType("HLT2");
// SingleCandNTuple_HLT1 = new RooDataSet("SingleCandNTuple_HLT1","SingleCandNTuple_HLT1",Vars);
// SingleCandNTuple_HLT2 = new RooDataSet("SingleCandNTuple_HLT2","SingleCandNTuple_HLT2",Vars);
// if      (NTuple->TrigCat == 1) SingleCandNTuple_HLT1->add(Vars);
// else if (NTuple->TrigCat == 2) SingleCandNTuple_HLT2->add(Vars);
// SingleCandNTuple = new RooDataSet("SingleCandNTuple","SingleCandNTuple",Vars,Index(*HLTcategory),
// 				  Import("HLT1",*SingleCandNTuple_HLT1),
// 				  Import("HLT2",*SingleCandNTuple_HLT2),
// 				  Import("HLT3",*SingleCandNTuple_HLT3),
// 				  Import("HLT4",*SingleCandNTuple_HLT4));
// ...
// RooSimultaneous* FullHLTpdf[q2Bins->size()-1];
// InstantiateMassAngleFit(&TotalPDFq2BinsHLT1[i],useEffPDF,x,y,myString.str(),FitType,configParam,fitParam,(ID == 0 ? i : 0),effFuncs1->operator[](i),true);
// InstantiateMassAngleFit(&TotalPDFq2BinsHLT2[i],useEffPDF,x,y,myString.str(),FitType,configParam,fitParam,(ID == 0 ? i : 0),effFuncs2->operator[](i),false);
// FullHLTpdf[i] = new RooSimultaneous("FullHLTpdf","FullHLTpdf",*HLTcategory);
// FullHLTpdf[i]->addPdf(*TotalPDFq2BinsHLT1[i],"HLT1");
// FullHLTpdf[i]->addPdf(*TotalPDFq2BinsHLT2[i],"HLT2");
// TotalPDFq2Bins = FullHLTpdf[i];
// ############################################################################

#ifndef __CINT__
#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF2.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TLorentzVector.h>
#include <TMinuit.h>
#include <TLatex.h>

#include <RooRealVar.h>
#include <RooGaussian.h>
#include <RooPoisson.h>
#include <RooAddPdf.h>
#include <RooDataSet.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooArgList.h>
#include <RooFitResult.h>
#include <RooProdPdf.h>
#include <RooPolynomial.h>
#include <RooMCStudy.h>
#include <RooMinuit.h>
#include <RooWorkspace.h>
#include <Roo1DTable.h>
#include <RooConstVar.h>
#include <RooRandom.h>
#endif

#include <iostream>
#include <sstream>
#include <time.h>

#include "ReadParameters.h"
#include "Utils.h"
#include "B0KstMuMuSingleCandTreeContent.h"

using namespace std;
using namespace RooFit;


// ####################
// # Global constants #
// ####################
#define NBINS 20
#define MULTYIELD 1.0       // Multiplication factor to the number of entry in toy-MC
#define NCOEFFPOLYBKG 5     // Maximum number of coefficients (= degree+1) of the polynomial describing the background in the angular variables
#define SLEWRATECONSTR 1e-2 // Fraction of Afb allowed interval for the p.d.f. constraint to transit from 0 to max
#define NORMBIN 3           // Bin used to normalize the yields to compute dBF/dq^2
#define POLYCOEFRANGE 2.0   // Polynomial coefficients range for parameter regeneration

#define nJPSIS 70000.0
#define nJPSIB  1000.0
#define nPSIPS  5000.0
#define nPSIPB   500.0

// ##########################################
// # Internal flags to control the workflow #
// ##########################################
#define ApplyConstr   true  // Apply Gaussian constraints in the likelihood
#define SAVEPOLY      false // "true" = save bkg polynomial coefficients in new parameter file; "false" = save original values
#define RESETOBS      false // Reset polynomial coefficients for combinatorial background and angular observables for a fresh start
#define SETBATCH      false // Set batch mode
#define TESTeffFUNC   false // Check whether efficiency goes negative
#define SAVEPLOT      false
#define FUNCERRBAND   false // Show the p.d.f. error band
#define MakeMuMuPlots false
#define MAKEGRAPHSCAN false // Make graphical scan of the physics-pdf*eff or physics-pdf alone (ony valid for GEN fit type options)
#define USEMINOS      true
#define UseSPwave     true

// ##################
// # External files #
// ##################
#define PARAMETERFILEIN  "../python/ParameterFile.txt"
#define PARAMETERFILEOUT "ParameterFile.txt"
#define FitSysFILEOutput "FitSystematics_q2Bin"
#define FitSysFILEInput  "../../Efficiency/EffSystematics/FitSystematics_q2Bin" // The protocol for this file is: ID --- FL --- AFB --- dBF/dq^2 --- NLL


// ####################
// # Global variables #
// ####################
TTree* theTree;
Utils* Utility;
B0KstMuMuSingleCandTreeContent* NTuple;

ofstream fileFitResults;
ofstream fileFitSystematics;
ifstream fileFitSystematicsInput;

double PsiYield;
double PsiYerr;

double LUMI;

string ParameterFILE;
TString legNames[5];

unsigned int nToy;
unsigned int nEntryToy; // nEntryToy = nSIG + nBKGCOMB + nBKGPEAK

double* q2BinsHisto;
vector<double> q2Bins;
vector<double> cosThetaKBins;
vector<double> cosThetaLBins;
vector<double> phiBins;

vector<vector<string>*> fitParam;          // Vector containing the pointers to the vectors containing the starting values for the fit
vector<vector<unsigned int>*> configParam; // Vector containing the pointers to the vectors containing the configuration parameters for the fit
vector<TF2*> effFuncs;                     // Vector containing the analytical description of the efficiency per q^2 bin


// #################################
// # Variables holding fit results #
// #################################
double p1PolyTot[NCOEFFPOLYBKG];
double c1PolyTot[NCOEFFPOLYBKG];
double p2PolyTot[NCOEFFPOLYBKG];
double c2PolyTot[NCOEFFPOLYBKG];


// ####################################
// # Useful variables from the NTuple #
// ####################################
RooDataSet* SingleCandNTuple_JPsi;
RooDataSet* SingleCandNTuple_PsiP;
RooDataSet* SingleCandNTuple_RejectPsi;
RooDataSet* SingleCandNTuple_KeepPsi;

RooRealVar* B0MassArb;
RooRealVar* mumuMass;
RooRealVar* mumuMassE;
RooRealVar* CosThetaKArb;
RooRealVar* CosThetaMuArb;
RooRealVar* PhiKstMuMuPlaneArb;
RooRealVar* truthMatchSignal;


// #################################
// # Variables and pdf for the fit #
// #################################

// ##################
// # Signal B0 mass #
// ##################
RooRealVar* meanS;

RooRealVar* sigmaS1;
RooAbsPdf* MassS1;

RooRealVar* sigmaS2;
RooAbsPdf* MassS2;

RooRealVar* fracMassS;
RooAbsPdf* MassSignal;

// #################
// # Signal angles #
// #################
RooRealVar* FlS;
RooRealVar* AfbS;
RooRealVar* At2S;
RooRealVar* AtimS;
RooRealVar* FsS;
RooRealVar* AsS;
RooAbsPdf* AngleS;

// ####################
// # Total signal pdf #
// ####################
RooAbsPdf* Signal;

// ######################
// # Background B0 mass #
// ######################
RooRealVar* tau1;
RooRealVar* tau2;
RooAbsPdf* BkgMassExp1;
RooAbsPdf* BkgMassExp2;

RooRealVar* fracMassBExp;
RooAbsPdf* BkgMassComb;

RooRealVar* meanR1;
RooRealVar* sigmaR1;
RooRealVar* meanR2;
RooRealVar* sigmaR2;
RooAbsPdf* BkgMassRPeak1;
RooAbsPdf* BkgMassRPeak2;

RooRealVar* fracMassBRPeak;
RooAbsPdf* BkgMassRPeak;

RooRealVar* meanL1;
RooRealVar* sigmaL1;
RooRealVar* meanL2;
RooRealVar* sigmaL2;
RooAbsPdf* BkgMassLPeak1;
RooAbsPdf* BkgMassLPeak2;

RooRealVar* fracMassBLPeak;
RooAbsPdf* BkgMassLPeak;

RooRealVar* fracMassBPeak;
RooAbsPdf* BkgMassPeak;

// #####################
// # Background angles #
// #####################
RooRealVar* p1Poly[NCOEFFPOLYBKG];
RooAbsPdf* BkgAngleP1;
RooRealVar* c1Poly[NCOEFFPOLYBKG];
RooAbsPdf* BkgAngleC1;

RooRealVar* p2Poly[NCOEFFPOLYBKG];
RooAbsPdf* BkgAngleP2;
RooRealVar* c2Poly[NCOEFFPOLYBKG];
RooAbsPdf* BkgAngleC2;

RooAbsPdf* BkgAnglesC;
RooAbsPdf* BkgAnglesP;

// ########################
// # Total background pdf #
// ########################
RooAbsPdf* BkgMassAngleComb;
RooAbsPdf* BkgMassAnglePeak;

RooRealVar* nSig;
RooRealVar* nBkgComb;
RooRealVar* nBkgPeak;

// ####################################################
// # Total likelihood for B0 --> K*0 J/psi OR psi(2S) #
// ####################################################
RooAbsPdf* TotalPDFPsi;

// #########################################
// # Total likelihood for B0 --> K*0 mu mu #
// #########################################
RooAbsPdf* TotalPDFRejectPsi;

// ##############################################
// # Vector containing the Gaussian constraints #
// ##############################################
RooArgSet vecConstr;


// #######################
// # Function Definition #
// #######################
bool CheckGoodFit (RooFitResult* fitResult, TPaveText* paveText = NULL);
RooRealVar* GetVar (RooAbsPdf* pdf, string varName);
void SetValueAndErrors (RooAbsPdf* pdf, string varName, double multi, stringstream* myString, double* val, double* errLo, double* errHi);
string MakeName (RooDataSet* data, int ID);
void DrawString (double Lumi, RooPlot* myFrame = NULL);
void PrintVariables (RooArgSet* setVar, string type);
void ClearVars (RooArgSet* vecConstr);
bool IsInConstraints (RooArgSet* vecConstr, string varName);
void AddPoissonConstraint (RooArgSet* vecConstr, RooAbsPdf* TotalPDF, string varName);
void AddGaussConstraint (RooArgSet* vecConstr, RooAbsPdf* TotalPDF, string varName);
void AddPhysicsConstraint (RooArgSet* vecConstr, RooAbsPdf* TotalPDF, string varName);
void BuildMassConstraints (RooArgSet* vecConstr, RooAbsPdf* TotalPDF, string varName);
void BuildAngularConstraints (RooArgSet* vecConstr, RooAbsPdf* TotalPDF);
void BuildPhysicsConstraints (RooArgSet* vecConstr, RooAbsPdf* TotalPDF, string varName);
string MakeAngWithEffPDF (TF2* effFunc, RooRealVar* x, RooRealVar* y, RooRealVar* z, unsigned int FitType, bool useEffPDF, RooArgSet* Vars);
void DeleteFit (RooAbsPdf* TotalPDF, string DeleteType);
void ResetCombPolyParam (vector<vector<string>*>* fitParam, vector<double>* q2Bins);
void ResetAngularParam (vector<vector<string>*>* fitParam, vector<double>* q2Bins);
double StoreFitResultsInFile (RooAbsPdf** TotalPDF, RooFitResult* fitResult, RooDataSet* dataSet, RooArgSet* vecConstr);
void StorePolyResultsInFile (RooAbsPdf** TotalPDF);
vector<string>* SaveFitResults (RooAbsPdf* TotalPDF, unsigned int fitParamIndx, vector<vector<string>*>* fitParam, vector<vector<unsigned int>*>* configParam, RooArgSet* vecConstr);
unsigned int CopyFitResults (RooAbsPdf* TotalPDF, unsigned int fitParamIndx, vector<vector<string>*>* fitParam);
void GenerateParameterFile (RooAbsPdf* TotalPDF, vector<vector<string>*>* fitParam, vector<vector<unsigned int>*>* configParam, RooArgSet* vecConstr, string fileName, unsigned int fileIndx, vector<double>* q2Bins, unsigned int q2BinIndx);
void GenerateDataset (RooAbsPdf* TotalPDF, RooArgSet setVar, vector<double>* q2Bins, int specBin, vector<vector<string>*>* fitParam, string fileName);
void MakeGraphicalScan (RooAbsPdf* TotalPDF, TCanvas* Canv, RooRealVar* x, RooRealVar* y, unsigned int NBins = 25);
void FitDimuonInvMass (RooDataSet* dataSet, RooAbsPdf** TotalPDFJPsi, RooAbsPdf** TotalPDFPsiP, RooRealVar* x, TCanvas* Canv, bool justPlotMuMuMass, bool justKeepPsi, string plotName);
void MakeDataSets (B0KstMuMuSingleCandTreeContent* NTuple, unsigned int FitType);
void InstantiateMassFit (RooAbsPdf** TotalPDF, RooRealVar* x, string fitName, vector<vector<unsigned int>*>* configParam, int parIndx);
RooFitResult* MakeMassFit (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* x, TCanvas* Canv, unsigned int FitType, vector<vector<unsigned int>*>* configParam, int parIndx, RooArgSet* vecConstr, double* NLLvalue, TPaveText* extText, int ID = 0);
void IterativeMassFitq2Bins (RooDataSet* dataSet,
			     bool useEffPDF,
			     double PsiYield, double PsiYieldErr,
			     RooRealVar* x, RooRealVar* y, RooRealVar* z,
			     int specBin,
			     unsigned int FitType,
			     vector<TH1D*>* VecHistoMeas,
			     vector<double>* q2Bins,
			     vector<vector<unsigned int>*>* configParam, vector<vector<string>*>* fitParam,
			     vector<TF2*>* effFuncs,
			     RooArgSet* vecConstr,
			     unsigned int ID = 0);
void MakeMassToy (RooAbsPdf* TotalPDF, RooRealVar* x, TCanvas* Canv, unsigned int FitType, unsigned int nToy, vector<vector<unsigned int>*>* configParam, int specBin, vector<vector<string>*>* fitParam, RooArgSet* vecConstr, string fileName);
void InstantiateMassAngleFit (RooAbsPdf** TotalPDF, bool useEffPDF, RooRealVar* x, RooRealVar* y, string fitName, unsigned int FitType, vector<vector<unsigned int>*>* configParam, vector<vector<string>*>* fitParam, unsigned int parIndx, TF2* effFunc, bool isFirstDataset);
RooFitResult* MakeMassAngleFit (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* x, RooRealVar* y, TCanvas* Canv, unsigned int FitType, vector<vector<unsigned int>*>* configParam, unsigned int parIndx, RooArgSet* vecConstr, double* NLLvalue, TPaveText* extText, int ID = 0);
void IterativeMassAngleFitq2Bins (RooDataSet* dataSet,
				  bool useEffPDF,
				  double PsiYield, double PsiYieldErr,
				  RooRealVar* x, RooRealVar* y, RooRealVar* z,
				  int specBin,
				  unsigned int FitType,
				  vector<TH1D*>* VecHistoMeas,
				  vector<double>* q2Bins,
				  vector<vector<unsigned int>*>* configParam, vector<vector<string>*>* fitParam,
				  vector<TF2*>* effFuncs,
				  RooArgSet* vecConstr,
				  unsigned int ID = 0);
void MakeMassAngleToy (RooAbsPdf* TotalPDF, RooRealVar* x, RooRealVar* y, TCanvas* Canv, unsigned int FitType, unsigned int nToy, vector<vector<unsigned int>*>* configParam, int specBin, vector<vector<string>*>* fitParam, RooArgSet* vecConstr, string fileName);
void InstantiateMass2AnglesFit (RooAbsPdf** TotalPDF,
				bool useEffPDF,
				RooRealVar* x, RooRealVar* y, RooRealVar* z,
				string fitName, unsigned int FitType,
				vector<vector<unsigned int>*>* configParam,
				vector<vector<string>*>* fitParam,
				unsigned int parIndx,
				TF2* effFunc,
				bool isFirstDataset);
RooFitResult* MakeMass2AnglesFit (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* x, RooRealVar* y, RooRealVar* z, TCanvas* Canv, unsigned int FitType, vector<vector<unsigned int>*>* configParam, unsigned int parIndx, RooArgSet* vecConstr, double* NLLvalue, TPaveText* extText, int ID = 0);
void IterativeMass2AnglesFitq2Bins (RooDataSet* dataSet,
				    bool useEffPDF,
				    double PsiYield, double PsiYieldErr,
				    RooRealVar* x, RooRealVar* y, RooRealVar* z,
				    int specBin,
				    unsigned int FitType,
				    vector<TH1D*>* VecHistoMeas,
				    vector<double>* q2Bins,
				    vector<vector<unsigned int>*>* configParam, vector<vector<string>*>* fitParam,
				    vector<TF2*>* effFuncs,
				    RooArgSet* vecConstr,
				    unsigned int ID = 0);
void MakeMass2AnglesToy (RooAbsPdf* TotalPDF, RooRealVar* x, RooRealVar* y, RooRealVar* z, TCanvas* Canv, unsigned int FitType, unsigned int nToy, vector<vector<unsigned int>*>* configParam, int specBin, vector<vector<string>*>* fitParam, RooArgSet* vecConstr, string fileName);
void CloseAllAndQuit (TApplication* theApp, TFile* NtplFile);


// ###########################
// # Function Implementation #
// ###########################
bool CheckGoodFit (RooFitResult* fitResult, TPaveText* paveText)
{
  if (fitResult != NULL)
    {
      if ((fitResult->covQual() == 3) && (fitResult->status() == 0))
	{
	  if (paveText != NULL) paveText->AddText("Fit status: GOOD");
	  return true;
	}
      else
	{
	  if (paveText != NULL) paveText->AddText("Fit status: BAD");
	  return false;
	}
    }

  return false;
}


RooRealVar* GetVar (RooAbsPdf* pdf, string varName)
{
  return (RooRealVar*)(pdf->getVariables()->find(varName.c_str()));
}


void SetValueAndErrors (RooAbsPdf* pdf, string varName, double multi, stringstream* myString, double* val, double* errLo, double* errHi)
// #############################################################################
// # If the error is an empty string --> set errLo = -1 and errHi = +1         #
// # If both errLo and errHi are 0.0 --> setAsymError = -1/+1 and setError = 1 #
// #############################################################################
{
  string tmpStr;

  if (myString->str().empty() == false)
    {
      tmpStr.clear();
      (*myString) >> tmpStr;
      *val = atof(tmpStr.c_str()) * multi;
      
      tmpStr.clear();
      (*myString) >> tmpStr;
      if (tmpStr.empty() == true) *errLo = -1.0;
      else *errLo = atof(tmpStr.c_str()) * multi;
      
      tmpStr.clear();
      (*myString) >> tmpStr;
      if (tmpStr.empty() == true) *errHi = 1.0;
      else *errHi = atof(tmpStr.c_str()) * multi;
    }

  if ((pdf != NULL) && (GetVar(pdf,varName) != NULL))
    {
      pdf->getVariables()->setRealValue(varName.c_str(),*val);
      if ((*errLo == 0.0) && (*errHi == 0.0))
      	{
      	  GetVar(pdf,varName)->setAsymError(-1.0,1.0);
      	  GetVar(pdf,varName)->setError(1.0);
      	}
      else
      	{
      	  GetVar(pdf,varName)->setAsymError(-1.0*fabs(*errLo),fabs(*errHi));
      	  GetVar(pdf,varName)->setError((fabs(*errLo) + fabs(*errHi)) / 2.0);
      	}
    }
}


string MakeName (RooDataSet* data, int ID)
{
  stringstream myString;
  myString.clear(); myString.str("");
  myString << data->GetName() << "_" << ID;
  return myString.str();
}


void DrawString (double Lumi, RooPlot* myFrame)
{
  stringstream myString;

  myString.clear(); myString.str("");
  myString << "CMS";
  TLatex* LumiTex1 = new TLatex(0.1,0.91,myString.str().c_str());
  LumiTex1->SetTextSize(0.05);
  LumiTex1->SetTextColor(kBlack);
  LumiTex1->SetNDC(true);
  if (myFrame == NULL) LumiTex1->DrawLatex(0.1,0.91,myString.str().c_str());
  else
    {
      LumiTex1->Paint();
      myFrame->addObject(LumiTex1);
    }

  myString.clear(); myString.str("");
  myString << "L = " << Lumi <<  " fb#lower[0.4]{^{#font[122]{\55}1}}";
  TLatex* LumiTex2 = new TLatex(0.43,0.91,myString.str().c_str());
  LumiTex2->SetTextSize(0.05);
  LumiTex2->SetTextColor(kBlack);
  LumiTex2->SetNDC(true);
  if (myFrame == NULL) LumiTex2->DrawLatex(0.43,0.91,myString.str().c_str());
  else
    {
      LumiTex2->Paint();
      myFrame->addObject(LumiTex2);
    }

  // ##################
  // # Custom method: #
  // ##################
  double startNDCx = 0.826;
  double startNDCy = 0.935;
  TLine* line1 = new TLine(startNDCx-0.004, startNDCy, startNDCx, startNDCy);
  line1->SetBit(TLine::kLineNDC,true);
  if (myFrame == NULL) line1->Draw();
  else
    {
      line1->Paint();
      myFrame->addObject(line1);
    }
  TLine* line2 = new TLine(startNDCx, startNDCy, startNDCx+0.005, startNDCy-0.03);
  line2->SetBit(TLine::kLineNDC,true);
  if (myFrame == NULL) line2->Draw();
  else
    {
      line2->Paint();
      myFrame->addObject(line2);
    }
  TLine* line3 = new TLine(startNDCx+0.005, startNDCy-0.03, startNDCx+0.010, startNDCy+0.01);
  line3->SetBit(TLine::kLineNDC,true);
  if (myFrame == NULL) line3->Draw();
  else
    {
      line3->Paint();
      myFrame->addObject(line3);
    }
  TLine* line4 = new TLine(startNDCx+0.010, startNDCy+0.01, startNDCx+0.032, startNDCy+0.01);
  line4->SetBit(TLine::kLineNDC,true);
  if (myFrame == NULL) line4->Draw();
  else
    {
      line4->Paint();
      myFrame->addObject(line4);
    }
  // ###################
  // # Nominal method: #
  // ###################
  // myString.clear(); myString.str("");
  // myString << "#sqrt{  }";
  // TLatex* LumiTex3 = new TLatex(0.82,0.9,myString.str().c_str());
  // LumiTex3->SetTextSize(0.053);
  // LumiTex3->SetTextColor(kBlack);
  // LumiTex3->SetNDC(true);
  // if (myFrame == NULL) LumiTex3->DrawLatex(0.82,0.9,myString.str().c_str());
  // else
  //   {
  //     LumiTex3->Paint();
  //     myFrame->addObject(LumiTex3);
  //   }

  myString.clear(); myString.str("");
  myString << "s = 7 TeV";
  TLatex* LumiTex4 = new TLatex(0.84,0.91,myString.str().c_str());
  LumiTex4->SetTextSize(0.05);
  LumiTex4->SetTextColor(kBlack);
  LumiTex4->SetNDC(true);
  if (myFrame == NULL) LumiTex4->DrawLatex(0.84,0.91,myString.str().c_str());
  else
    {
      LumiTex4->Paint();
      myFrame->addObject(LumiTex4);
    }
}


void PrintVariables (RooArgSet* setVar, string type)
{
  RooRealVar* tmpVar;
  int nEleSet = setVar->getSize();

  if (type == "vars")
    {
      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout <<   "@@@ Printing variables @@@" << endl;
      cout <<   "@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      
      TIterator* it = setVar->createIterator();
      for (int i = 0; i < nEleSet; i++)
	{
	  tmpVar = (RooRealVar*)it->Next();
	  cout << "Variable: " << i;
	  cout << "\tname: "   << tmpVar->GetName();
	  cout << "\tvalue: "  << tmpVar->getVal();
	  cout << "\terr: "    << tmpVar->getError();
	  cout << "\tErrLo: "  << tmpVar->getErrorLo();
	  cout << "\tErrHi: "  << tmpVar->getErrorHi() << endl;
	}
    }
  else if (type == "cons")
    {
      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout <<   "@@@@@@@@@ Printing constraints @@@@@@@@@" << endl;
      cout <<   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

      TIterator* it = setVar->createIterator();
      for (int i = 0; i < nEleSet; i++) PrintVariables(((RooAbsPdf*)it->Next())->getVariables(),"vars");
    }
  else
    {
      cout << "[ExtractYield::PrintVariables]\tWrong parameter: " << type << endl;
      exit(1);
    }
}


void ClearVars (RooArgSet* vecConstr)
{
  int nEle = vecConstr->getSize();

  TIterator* it = vecConstr->createIterator();
  for (int i = 0; i < nEle; i++) delete it->Next();

  vecConstr->removeAll();
}


bool IsInConstraints (RooArgSet* vecConstr, string varName)
{
  string tmp;
  int nEle = vecConstr->getSize();

  TIterator* it = vecConstr->createIterator();
  for (int i = 0; i < nEle; i++)
    {
      tmp.clear();
      tmp = it->Next()->GetName();
      if (tmp.replace(tmp.find("_constr"),7,"") == varName) return true;
    }

  return false;
}


void AddPoissonConstraint (RooArgSet* vecConstr, RooAbsPdf* TotalPDF, string varName)
{
  stringstream myString;
  
  myString.clear(); myString.str("");
  myString << varName;

  RooRealVar* varConstr = GetVar(TotalPDF,myString.str().c_str());
  double mean  = TotalPDF->getVariables()->getRealValue(myString.str().c_str());
  
  myString << "_constr";
  RooPoisson* newConstr = new RooPoisson(myString.str().c_str(), myString.str().c_str(), *varConstr, RooConst(mean));
  vecConstr->add(*newConstr);
}


void AddGaussConstraint (RooArgSet* vecConstr, RooAbsPdf* TotalPDF, string varName)
{
  stringstream myString;
  
  myString.clear(); myString.str("");
  myString << varName;

  RooRealVar* varConstr = GetVar(TotalPDF,myString.str().c_str());
  double mean  = TotalPDF->getVariables()->getRealValue(myString.str().c_str());
  double sigma = (fabs(varConstr->getErrorLo()) + fabs(varConstr->getErrorHi())) / 2.0;
  
  myString << "_constr";
  RooGaussian* newConstr = new RooGaussian(myString.str().c_str(), myString.str().c_str(), *varConstr, RooConst(mean), RooConst(sigma));
  vecConstr->add(*newConstr);
}


void AddPhysicsConstraint (RooArgSet* vecConstr, RooAbsPdf* TotalPDF, string varName)
{
  stringstream myString;
  stringstream myCoeff;
  stringstream myConstr;
  string name;
  RooGenericPdf* newConstr = NULL;

  myConstr.clear(); myConstr.str("");
  if (varName == "AfbBound")
    {
      myConstr << "3/4 * (1 - " << GetVar(TotalPDF,"FlS")->getPlotLabel() << ")";
      name = "AfbS";
    }
  else if (varName == "AsBound")
    {
      myConstr << "1/2 * (" << GetVar(TotalPDF,"FsS")->getPlotLabel() << " + 3 * " << GetVar(TotalPDF,"FlS")->getPlotLabel() << " * (1-" << GetVar(TotalPDF,"FsS")->getPlotLabel() << "))";
      name = "AsS";
    }
  else
    {
      cout << "[ExtractYield::AddPhysicsConstraint]\tWrong parameter: " << varName << endl;
      exit(1);
    }

  myString.clear(); myString.str("");
  myString << varName;
  myString << "_constr";

  myCoeff.clear(); myCoeff.str("");
  myCoeff << "(1-TMath::Erf((-" << GetVar(TotalPDF,name)->getPlotLabel() << "-" << myConstr.str() << ") / ";
  myCoeff << "(" << SLEWRATECONSTR << " * 2*" << myConstr.str() << "))) * ";
  myCoeff << "(1-TMath::Erf(( " << GetVar(TotalPDF,name)->getPlotLabel() << "-" << myConstr.str() << ") / ";
  myCoeff << "(" << SLEWRATECONSTR << " * 2*" << myConstr.str() << ")))";
  cout << "\n@@@ Physical constraint on Afb/As boundaries: " << myCoeff.str().c_str() << " @@@" << endl;
  
  if (varName == "AfbBound")
    newConstr = new RooGenericPdf(myString.str().c_str(),myCoeff.str().c_str(),RooArgSet(*GetVar(TotalPDF,name),*GetVar(TotalPDF,"FlS")));
  else if (varName == "AsBound")
    newConstr = new RooGenericPdf(myString.str().c_str(),myCoeff.str().c_str(),RooArgSet(*GetVar(TotalPDF,name),*GetVar(TotalPDF,"FlS"),*GetVar(TotalPDF,"FsS")));
  vecConstr->add(*newConstr);
}


void BuildMassConstraints (RooArgSet* vecConstr, RooAbsPdf* TotalPDF, string varName)
// ############################################################
// # All:  apply constraint to all physical angular variables #
// # sign: apply constraint to signal                         #
// # peak: apply constraint to peaking background             #
// ############################################################
{
  if ((GetVar(TotalPDF,"sigmaS1")        != NULL) && ((varName == "All") || (varName == "sign"))) AddGaussConstraint(vecConstr, TotalPDF, "sigmaS1");
  if ((GetVar(TotalPDF,"sigmaS2")        != NULL) && ((varName == "All") || (varName == "sign"))) AddGaussConstraint(vecConstr, TotalPDF, "sigmaS2");
  if ((GetVar(TotalPDF,"fracMassS")      != NULL) && ((varName == "All") || (varName == "sign"))) AddGaussConstraint(vecConstr, TotalPDF, "fracMassS");
  
  if ((GetVar(TotalPDF,"meanR1")         != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, TotalPDF, "meanR1");
  if ((GetVar(TotalPDF,"sigmaR1")        != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, TotalPDF, "sigmaR1");
  if ((GetVar(TotalPDF,"meanR2")         != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, TotalPDF, "meanR2");
  if ((GetVar(TotalPDF,"sigmaR2")        != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, TotalPDF, "sigmaR2");
  if ((GetVar(TotalPDF,"fracMassBRPeak") != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, TotalPDF, "fracMassBRPeak");

  if ((GetVar(TotalPDF,"meanL1")         != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, TotalPDF, "meanL1");
  if ((GetVar(TotalPDF,"sigmaL1")        != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, TotalPDF, "sigmaL1");
  if ((GetVar(TotalPDF,"meanL2")         != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, TotalPDF, "meanL2");
  if ((GetVar(TotalPDF,"sigmaL2")        != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, TotalPDF, "sigmaL2");
  if ((GetVar(TotalPDF,"fracMassBLPeak") != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, TotalPDF, "fracMassBLPeak");

  if ((GetVar(TotalPDF,"fracMassBPeak")  != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, TotalPDF, "fracMassBPeak");

  // @TMP@
  if ((GetVar(TotalPDF,"nBkgPeak")       != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, TotalPDF, "nBkgPeak");
  // if ((GetVar(TotalPDF,"nBkgPeak")       != NULL) && ((varName == "All") || (varName == "peak"))) AddPoissonConstraint(vecConstr, TotalPDF, "nBkgPeak");
}


void BuildAngularConstraints (RooArgSet* vecConstr, RooAbsPdf* TotalPDF)
{
  stringstream myString;
  
  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
    {
      myString.clear(); myString.str("");
      myString << "p1Poly" << i;
      if (GetVar(TotalPDF,myString.str().c_str()) != NULL) AddGaussConstraint(vecConstr, TotalPDF, myString.str().c_str());
    }
  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
    {
      myString.clear(); myString.str("");
      myString << "p2Poly" << i;
      if (GetVar(TotalPDF,myString.str().c_str()) != NULL) AddGaussConstraint(vecConstr, TotalPDF, myString.str().c_str());
    }
}


void BuildPhysicsConstraints (RooArgSet* vecConstr, RooAbsPdf* TotalPDF, string varName)
// ##################################################################
// # FlS:   apply gaussian constraint to FlS                        #
// # AfbS:  apply gaussian constraint to AfbS                       #
// # At2S:  apply gaussian constraint to At2S                       #
// # AtimS: apply gaussian constraint to AtimS                      #
// ##################################################################
// # FlAfb: apply physics constraint to Afb boundaries by Fl values #
// ##################################################################
{
  if ((GetVar(TotalPDF,"FlS")   != NULL) && (varName == "FlS"))   AddGaussConstraint(vecConstr, TotalPDF, varName);
  if ((GetVar(TotalPDF,"AfbS")  != NULL) && (varName == "AfbS"))  AddGaussConstraint(vecConstr, TotalPDF, varName);
  if ((GetVar(TotalPDF,"At2S")  != NULL) && (varName == "At2S"))  AddGaussConstraint(vecConstr, TotalPDF, varName);
  if ((GetVar(TotalPDF,"AtimS") != NULL) && (varName == "AtimS")) AddGaussConstraint(vecConstr, TotalPDF, varName);

  if ((GetVar(TotalPDF,"FsS")   != NULL) && (varName == "FsS"))   AddGaussConstraint(vecConstr, TotalPDF, varName);
  if ((GetVar(TotalPDF,"AsS")   != NULL) && (varName == "AsS"))   AddGaussConstraint(vecConstr, TotalPDF, varName);


  // #########################################
  // # Special constraint for AFB boundaries #
  // #########################################
  if ((varName == "AfbBound") && (GetVar(TotalPDF,"FlS") != NULL) && (GetVar(TotalPDF,"AfbS") != NULL)) AddPhysicsConstraint(vecConstr, TotalPDF, varName);


  // ########################################
  // # Special constraint for AS boundaries #
  // ########################################
  if ((varName == "AsBound") && (GetVar(TotalPDF,"FlS") != NULL) && (GetVar(TotalPDF,"FsS") != NULL) && (GetVar(TotalPDF,"AsS") != NULL)) AddPhysicsConstraint(vecConstr, TotalPDF, varName);
}


string MakeAngWithEffPDF (TF2* effFunc, RooRealVar* x, RooRealVar* y, RooRealVar* z, unsigned int FitType, bool useEffPDF, RooArgSet* Vars)
// ###################
// # x: angle phi    #
// # y: cos(theta_l) #
// # z: cos(theta_K) #
// ###################
{
  stringstream myString;
  vector<RooRealVar*> vecParam;


  // ############################################
  // # Copy parameters from efficiency function #
  // ############################################
  for (int i = 0; i < effFunc->GetNpar(); i++)
    {
      myString.clear(); myString.str("");
      myString << "P" << i;
      vecParam.push_back(new RooRealVar(myString.str().c_str(),myString.str().c_str(),effFunc->GetParameter(i),-1.0,1.0));
      vecParam.back()->setConstant(true);
    }


  if ((FitType == 3) || (FitType == 23) || (FitType == 33) || (FitType == 43) || (FitType == 53) || (FitType == 63) || (FitType == 73) || (FitType == 83) || (FitType == 93)) // Fl-Afb-fit
    {
      // ########################################################################
      // # Make 1D signal*efficiency p.d.f.: integral over phi and cos(theta_l) #
      // ########################################################################
      FlS = new RooRealVar("FlS","F_{L}",0.0,0.0,1.0);
      FlS->setConstant(false);
      Vars->add(*FlS);

      myString.clear(); myString.str("");
      if (useEffPDF == false)
	{
	  if (UseSPwave == false)
	    {
	      // #####################
	      // # P-wave decay rate #
	      // #####################
	      myString << "(3/2 * FlS * " << x->getPlotLabel() << "*" << x->getPlotLabel();
	      myString << " + 3/4 * (1-FlS) * (1 - " << x->getPlotLabel() << "*" << x->getPlotLabel() << "))";
	    }
	  else
	    {
	      FsS = new RooRealVar("FsS","F_{s}",0.0,0.0,1.0);
	      AsS = new RooRealVar("AsS","A_{s}",0.0,-1.0,1.0);
	      FsS->setConstant(false);
	      AsS->setConstant(false);
	      Vars->add(*FsS);
	      Vars->add(*AsS);

	      // ###########################
	      // # S and P-wave decay rate #
	      // ###########################
	      myString << "(1/4 * (2*FsS + 4*AsS*" << x->getPlotLabel() << " + ";
	      myString << "3*(1 - FsS)*(1 - FlS - " << x->getPlotLabel() << "*" << x->getPlotLabel() << " + 3*FlS*" << x->getPlotLabel() << "*" << x->getPlotLabel() << ")))";
	    }
	}
      else
	{
	  if (Utility->IsThereOddDegree(effFunc) == true)
	    {
	      AfbS = new RooRealVar("AfbS","A_{FB}",0.0,-1.0,1.0);
	      AfbS->setConstant(false);
	      Vars->add(*AfbS);
	    }

	  if (UseSPwave == false)
	    {
	      // #####################
	      // # P-wave decay rate #
	      // #####################
	      myString << "(3/4*(FlS*(3*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "-1)-" << x->getPlotLabel() << "*" << x->getPlotLabel() << "+1) * ";
	      myString << "(P0+" << x->getPlotLabel() << "*(P1+" << x->getPlotLabel() << "*(P2+P3*" << x->getPlotLabel() << "))) + ";

	      myString << "3/10*(FlS*(2*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "-1)-" << x->getPlotLabel() << "*" << x->getPlotLabel() << "+1) * ";
	      myString << "(P4+" << x->getPlotLabel() << "*(P5+" << x->getPlotLabel() << "*(P6+P7*" << x->getPlotLabel() << "))) + ";

	      if (Utility->IsThereOddDegree(effFunc) == true)
		{
		  myString << "3/10*AfbS*(" << x->getPlotLabel() << "*" << x->getPlotLabel() << "-1) * ";
		  myString << "(P8+" << x->getPlotLabel() << "*(P9+" << x->getPlotLabel() << "*(P10+P11*" << x->getPlotLabel() << "))) + ";
		}

	      myString << "9/140*(FlS*(5*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "-3)-3*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "+3) * ";
	      myString << "(P12+" << x->getPlotLabel() << "*(P13+" << x->getPlotLabel() << "*(P14+P15*" << x->getPlotLabel() << "))) + ";

	      myString << "1/14*(FlS*(3*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "-2)-2*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "+2) * ";
	      myString << "(P16+" << x->getPlotLabel() << "*(P17+" << x->getPlotLabel() << "*(P18+P19*" << x->getPlotLabel() << "))))";
	    }
	  else
	    {
	      FsS = new RooRealVar("FsS","F_{s}",0.0,0.0,1.0);
	      AsS = new RooRealVar("AsS","A_{s}",0.0,-1.0,1.0);
	      FsS->setConstant(false);
	      AsS->setConstant(false);
	      Vars->add(*FsS);
	      Vars->add(*AsS);

	      // ###########################
	      // # S and P-wave decay rate #
	      // ###########################
	      myString << "(1/420 * ";

	      if (Utility->IsThereOddDegree(effFunc) == true)
		{
		  // 0-degree z
		  myString << "((-FsS)*(63*P12 + 50*P16 + 84*P4 + 126*AfbS*P8) + 3*(27*P12 + 20*P16 + 42*(P4 + AfbS*P8)) + ";

		  // 1-degree z
		  myString << "315*P1 * " << x->getPlotLabel() << " + (4*AsS*(9*P12 + 5*P16 + 21*P4) - FsS*(105*P1 + 63*P13 + 50*P17 + 84*P5 + 126*AfbS*P9) + 3*(27*P13 + 20*P17 + 42*(P5 + AfbS*P9))) * " << x->getPlotLabel() << " + ";

		  // 2-degree z
		  myString << "(-81*P12 + 81*P14 - 60*P16 + 60*P18 + 315*P2 - 126*P4 + 4*AsS*(105*P1 + 9*P13 + 5*P17 + 21*P5) + 126*P6 + FsS*(81*P12 - 63*P14 + 60*P16 - 50*P18 - 21*(5*P2 - 6*P4 + 4*P6)) - 126*AfbS*(-1 + FsS)*(P10 - P8)) * " << x->getPlotLabel() << "*" << x->getPlotLabel() << " + ";

		  // 3-degree z
		  myString << "(315*(-1 + FsS)*P1 - 81*P13 + 4*AsS*(9*P14 + 5*P18 + 21*(5*P2 + P6)) + 3*(27*P15 - 20*P17 + 20*P19 + 105*P3 - 42*P5 + 42*P7) + FsS*(81*P13 - 63*P15 + 60*P17 - 50*P19 - 21*(5*P3 - 6*P5 + 4*P7)) - 126*AfbS*(-1 + FsS)*(P11 - P9)) * " << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << " + ";

		  // 4-degree z
		  myString << "(126*AfbS*(-1 + FsS)*P10 + 81*(-1 + FsS)*P14 + 3*(-1 + FsS)*(20*P18 + 105*P2 + 42*P6) + 4*AsS*(9*P15 + 5*P19 + 21*(5*P3 + P7))) * " << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << " + ";

		  // 5-degree z
		  myString << "3*(-1 + FsS)*(42*AfbS*P11 + 27*P15 + 20*P19 + 105*P3 + 42*P7) * " << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << " + ";

		  // mix-degree z
		  myString << "105*P0*(3 - 3*FlS - FsS + 3*FlS*FsS + 4*AsS*" << x->getPlotLabel() << " - 3*(-1 + 3*FlS)*(-1 + FsS)*" << x->getPlotLabel() << "*" << x->getPlotLabel() << ") + 3*FlS*(-1 + FsS)*(20*P16 + 42*P4 + (105*P1 + 27*P13 + 20*P17 + 42*P5)*" << x->getPlotLabel() << " + 9*P12*(3 - 5*" << x->getPlotLabel() << "*" << x->getPlotLabel() << ") + " << x->getPlotLabel() << "*" << x->getPlotLabel() << "*(-30*P16 + 20*P18 + 105*P2 - 84*P4 + 42*P6 + (-315*P1 - 45*P13 + 27*P15 - 30*P17 + 20*P19 + 105*P3 - 84*P5 + 42*P7)*" << x->getPlotLabel() << " - 3*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*(10*P18 + 105*P2 + 28*P6 + (15*P15 + 10*P19 + 105*P3 + 28*P7)*" << x->getPlotLabel() << ") + 9*P14*(3 - 5*" << x->getPlotLabel() << "*" << x->getPlotLabel() << ")))))";
		}
	      else
		{
		  // 0-degree z
		  myString << "(81*P12 + 60*P16 + 126*P4 - FsS*(63*P12 + 50*P16 + 84*P4) + ";
	  
		  // 1-degree z
		  myString << "(315*P1 + 81*P13 + 60*P17 + 4*AsS*(9*P12 + 5*P16 + 21*P4) + 126*P5 - FsS*(105*P1 + 63*P13 + 50*P17 + 84*P5))*" << x->getPlotLabel() << " + ";
	  
		  // 2-degree z
		  myString << "(81*(-1 + FsS)*P12 + 9*(9 - 7*FsS)*P14 - 60*P16 + 60*P18 + 315*P2 - 126*P4 + 4*AsS*(105*P1 + 9*P13 + 5*P17 + 21*P5) + 126*P6 + FsS*(60*P16 - 50*P18 - 21*(5*P2 - 6*P4 + 4*P6)))*" << x->getPlotLabel() << "*" << x->getPlotLabel() << " + ";
	  
		  // 3-degree z
		  myString << "(315*(-1 + FsS)*P1 + 81*(-1 + FsS)*P13 + 4*AsS*(9*P14 + 5*P18 + 21*(5*P2 + P6)) + 3*(27*P15 - 20*P17 + 20*P19 + 105*P3 - 42*P5 + 42*P7) - FsS*(63*P15 - 60*P17 + 50*P19 + 21*(5*P3 - 6*P5 + 4*P7)))*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << " + ";
	  
		  // 4-degree z
		  myString << "(81*(-1 + FsS)*P14 + 3*(-1 + FsS)*(20*P18 + 105*P2 + 42*P6) + 4*AsS*(9*P15 + 5*P19 + 21*(5*P3 + P7)))*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << " + ";
	  
		  // 5-degree z
		  myString << "3*(-1 + FsS)*(27*P15 + 20*P19 + 105*P3 + 42*P7)*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << " + ";
	  
		  // mix-degree z
		  myString << "105*P0*(3 - 3*FlS - FsS + 3*FlS*FsS + 4*AsS*" << x->getPlotLabel() << " - 3*(-1 + 3*FlS)*(-1 + FsS)*" << x->getPlotLabel() << "*" << x->getPlotLabel() << ") + ";
		  myString << "3*FlS*(-1 + FsS)*(20*P16 + 42*P4 + (105*P1 + 27*P13 + 20*P17 + 42*P5)*" << x->getPlotLabel() << " + 9*P12*(3 - 5*" << x->getPlotLabel() << "*" << x->getPlotLabel() << ") + ";
		  myString << x->getPlotLabel() << "*" << x->getPlotLabel() << "*(-30*P16 + 20*P18 + 105*P2 - 84*P4 + 42*P6 + (-315*P1 - 45*P13 + 27*P15 - 30*P17 + 20*P19 + 105*P3 - 84*P5 + 42*P7)*" << x->getPlotLabel() << " - ";
		  myString << "3*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*(10*P18 + 105*P2 + 28*P6 + (15*P15 + 10*P19 + 105*P3 + 28*P7)*" << x->getPlotLabel() << ") + 9*P14*(3 - 5*" << x->getPlotLabel() << "*" << x->getPlotLabel() << ")))))";
		}
	    }
	  for (int i = 0; i < effFunc->GetNpar(); i++) Vars->add(*vecParam[i]);
	}

      Vars->add(*x);


      cout << "\n@@@ 1D signal*efficiency p.d.f. @@@" << endl;
      cout << myString.str().c_str() << endl;
    }
  else if ((FitType == 4) || (FitType == 24) || (FitType == 34) || (FitType == 44) || (FitType == 54) || (FitType == 64) || (FitType == 74) || (FitType == 84) || (FitType == 94)) // Afb-Fl-fit
    {
      // ########################################################################
      // # Make 1D signal*efficiency p.d.f.: integral over phi and cos(theta_K) #
      // ########################################################################
      FlS = new RooRealVar("FlS","F_{L}",0.0,0.0,1.0);
      AfbS = new RooRealVar("AfbS","A_{FB}",0.0,-1.0,1.0);
      FlS->setConstant(false);
      AfbS->setConstant(false);
      Vars->add(*FlS);
      Vars->add(*AfbS);

      myString.clear(); myString.str("");
      if (useEffPDF == false)
	{
	  if (UseSPwave == false)
	    {
	      // #####################
	      // # P-wave decay rate #
	      // #####################
	      myString << "(3/4 * FlS * (1 - " << x->getPlotLabel() << "*" << x->getPlotLabel() << ")";
	      myString << " + 3/8 * (1-FlS) * (1 + " << x->getPlotLabel() << "*" << x->getPlotLabel() << ")";
	      myString << " + AfbS * " << x->getPlotLabel() << ")";
	    }
	  else
	    {
	      FsS = new RooRealVar("FsS","F_{s}",0.0,0.0,1.0);
	      FsS->setConstant(false);
	      Vars->add(*FsS);

	      // ###########################
	      // # S and P-wave decay rate #
	      // ###########################
	      myString << "(1/8 * (3*(1 + FlS + FsS - FlS*FsS) + ";
	      myString << "8*AfbS*(1 - FsS)*" << x->getPlotLabel() << " + ";
	      myString << "3*(1 - 3*FlS*(1 - FsS) - 3*FsS)*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "))";
	    }
	}
      else
	{
	  if (UseSPwave == false)
	    {
	      // #####################
	      // # P-wave decay rate #
	      // #####################
	      myString << "(1/40 * ";

	      myString << "(3*(5*P0+P2+5*(P0+P2)*FlS) + 8*(5*P0+P2)*AfbS*" << x->getPlotLabel() << " + ";
	      myString << "3*(P2+P0*(5-15*FlS) - 7*P2*FlS)*" << x->getPlotLabel() << "*" << x->getPlotLabel() << ") + ";
	  
	      myString << "1/40*" << x->getPlotLabel() << "*" << x->getPlotLabel() << " * ";
	      myString << "(3*(5*P4+P6+5*(P4+P6)*FlS) + 8*(5*P4+P6)*AfbS*" << x->getPlotLabel() << " + ";
	      myString << "3*(P6+P4*(5-15*FlS) - 7*P6*FlS)*" << x->getPlotLabel() << "*" << x->getPlotLabel() << ") + ";
	  
	      myString << "1/40*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << " * ";
	      myString << "(3*(5*P8+P10+5*(P8+P10)*FlS) + 8*(5*P8+P10)*AfbS*" << x->getPlotLabel() << " + ";
	      myString << "3*(P10+P8*(5-15*FlS) - 7*P10*FlS)*" << x->getPlotLabel() << "*" << x->getPlotLabel() << ") + ";

	      myString << "1/40*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << " * ";
	      myString << "(3*(5*P12+P14+5*(P12+P14)*FlS) + 8*(5*P12+P14)*AfbS*" << x->getPlotLabel() << " + ";
	      myString << "3*(P14+P12*(5-15*FlS) - 7*P14*FlS)*" << x->getPlotLabel() << "*" << x->getPlotLabel() << ") +";

	      myString << "1/40*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << " * ";
	      myString << "(3*(5*P16+P18+5*(P16+P18)*FlS) + 8*(5*P16+P18)*AfbS*" << x->getPlotLabel() << " + ";
	      myString << "3*(P18+P16*(5-15*FlS) - 7*P18*FlS)*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "))";
	    }
	  else
	    {
	      FsS = new RooRealVar("FsS","F_{s}",0.0,0.0,1.0);
	      AsS = new RooRealVar("AsS","A_{s}",0.0,-1.0,1.0);
	      FsS->setConstant(false);
	      AsS->setConstant(false);
	      Vars->add(*FsS);
	      Vars->add(*AsS);

	      // ###########################
	      // # S and P-wave decay rate #
	      // ###########################
	      myString << "(1/420 * ";

	      // 0-degree y
	      myString << "(3*P2 + 7*FsS*P2 + 4*AsS*(5*P1 + 3*P3) - ";

	      // 1-degree y
	      myString << "8*AfbS*(-1 + FsS)*P2*" << x->getPlotLabel() << " + ";

	      // 2-degree y
	      myString << "((3 - 13*FsS)*P2 + 15*(1 + FsS)*P4 + (3 + 7*FsS)*P6 - 4*AsS*(5*P1 + 3*P3 - 5*P5 - 3*P7))*" << x->getPlotLabel() << "*" << x->getPlotLabel() << " + ";

	      // 3-degree y
	      myString << "((3 + 7*FsS)*P10 - 8*AfbS*(-1 + FsS)*(5*P4 + P6) + 15*(1 + FsS)*P8 + 4*AsS*(3*P11 + 5*P9))*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << " + ";

	      // 4-degree y
	      myString << "(15*(1 + FsS)*P12 + FsS*(7*P14 - 45*P4 - 13*P6) + 3*(P14 + 5*P4 + P6) + 4*AsS*(5*P13 + 3*P15 - 5*P5 - 3*P7) - 8*AfbS*(-1 + FsS)*(P10 + 5*P8))*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << " - ";

	      // 5-degree y
	      myString << "((-3 + 13*FsS)*P10 + 8*AfbS*(-1 + FsS)*(5*P12 + P14) + 15*(-1 + 3*FsS)*P8 + 4*AsS*(3*P11 + 5*P9))*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << " + ";

	      // 6-degree y
	      myString << "((15 - 45*FsS)*P12 + (3 - 13*FsS)*P14 + 15*(1 + FsS)*P16 + (3 + 7*FsS)*P18 - 4*AsS*(5*P13 + 3*P15 - 5*P17 - 3*P19))*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << " - ";

	      // 7-degree y
	      myString << "8*AfbS*(-1 + FsS)*(5*P16 + P18)*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << " - ";

	      // 8-degree y
	      myString << "(15*(-1 + 3*FsS)*P16 + (-3 + 13*FsS)*P18 + 4*AsS*(5*P17 + 3*P19))*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << " + ";

	      // mix-degree y
	      myString << "5*P0*(3*(1 + FlS + FsS - FlS*FsS) - 8*AfbS*(-1 + FsS)*" << x->getPlotLabel() << " + 3*(1 + 3*FlS*(-1 + FsS) - 3*FsS)*" << x->getPlotLabel() << "*" << x->getPlotLabel() << ") + ";

	      // mix-degree y
	      myString << "3*FlS*(-1 + FsS)*(P2*(-5 + 7*" << x->getPlotLabel() << "*" << x->getPlotLabel() << ") + " << x->getPlotLabel() << "*" << x->getPlotLabel() << "*(5*P4*(-1 + 3*" << x->getPlotLabel() << "*" << x->getPlotLabel() << ") - 5*(P6 + " << x->getPlotLabel() << "*(P10 + P8 + (P12 + P14)*" << x->getPlotLabel() << ")) + " << x->getPlotLabel() << "*" << x->getPlotLabel() << "*(7*P6 + " << x->getPlotLabel() << "*(7*P10 + 15*P8 + (15*P12 + 7*P14 - 5*(P16 + P18))*" << x->getPlotLabel() << " + (15*P16 + 7*P18)*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "*" << x->getPlotLabel() << "))))))";
	    }
	  for (int i = 0; i < effFunc->GetNpar(); i++) Vars->add(*vecParam[i]);
	}

      Vars->add(*x);


      cout << "\n@@@ 1D signal*efficiency p.d.f. @@@" << endl;
      cout << myString.str().c_str() << endl;
    }
  else if ((FitType == 5) || (FitType == 25) || (FitType == 35) || (FitType == 45) || (FitType == 55) || (FitType == 65) || (FitType == 75) || (FitType == 85) || (FitType == 95)) // Fl-At2-Atim-fit
    {
      // #################################################################################
      // # Make 1D signal*efficiency p.d.f.: integral over cos(theta_K) and cos(theta_l) #
      // #################################################################################
      FlS   = new RooRealVar("FlS","F_{L}",0.0,0.0,1.0);
      At2S  = new RooRealVar("At2S","A_{T}#lower[0.4]{^{2}}",0.0,-1.0,1.0);
      AtimS = new RooRealVar("AtimS","A_{T}#lower[0.4]{^{Im}}",0.0,-1.0,1.0);
      FlS->setConstant(false);
      At2S->setConstant(false);
      AtimS->setConstant(false);
      Vars->add(*FlS);
      Vars->add(*At2S);
      Vars->add(*AtimS);

      myString.clear(); myString.str("");
      if (useEffPDF == false)
	{
	  if (UseSPwave == false)
	    // #####################
	    // # P-wave decay rate #
	    // #####################
	    myString << "(1/(2*" << Utility->PI << ") * (1 + 1/2 * (1-FlS) * At2S * cos(2*" << x->getPlotLabel() << ") + AtimS * sin(2*" << x->getPlotLabel() << ")))";
	  else
	    {
	      // ###########################
	      // # S and P-wave decay rate #
	      // ###########################
	      cout << "@TMP@: not implemented yet for fit type " << FitType << " (Fl-At2-Atim-fit)" << endl;
	      exit(1);
	    }
	}
      else
	{
	  cout << "@TMP@: not implemented yet for fit type " << FitType << " (Fl-At2-Atim-fit)" << endl;
	  exit(1);
	}
      
      Vars->add(*x);
      
      
      cout << "\n@@@ 1D angular*efficiency p.d.f. @@@" << endl;
      cout << myString.str().c_str() << endl; 
    }
  else if ((FitType == 1) || (FitType == 41) || (FitType == 61) || (FitType == 81) || // Branching fraction
	   (FitType == 6) || (FitType == 26) || (FitType == 36) || (FitType == 46) || (FitType == 56) || (FitType == 66) || (FitType == 76) || (FitType == 86) || (FitType == 96)) // Fl-Afb-fit AND Afb-Fl-fit
    {
      // #######################################################
      // # Make 2D signal*efficiency p.d.f.: integral over phi #
      // #######################################################
      FlS  = new RooRealVar("FlS","F_{L}",0.0,0.0,1.0);
      AfbS = new RooRealVar("AfbS","A_{FB}",0.0,-1.0,1.0);
      FlS->setConstant(false);
      AfbS->setConstant(false);
      Vars->add(*FlS);
      Vars->add(*AfbS);

      myString.clear(); myString.str("");
      if (UseSPwave == false)
	{
	  // #####################
	  // # P-wave decay rate #
	  // #####################
	  myString << "(3/4 * (3/2 * FlS * (1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * " << z->getPlotLabel() << "*" << z->getPlotLabel() << " + ";
	  myString << "(3/8 * (1-FlS) * (1+" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") + AfbS*" << y->getPlotLabel() << ") * ";
	  myString << "(1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")))";
	}
      else
	{
	  FsS = new RooRealVar("FsS","F_{s}",0.0,0.0,1.0);
	  AsS = new RooRealVar("AsS","A_{s}",0.0,-1.0,1.0);
	  FsS->setConstant(false);
	  AsS->setConstant(false);
	  Vars->add(*FsS);
	  Vars->add(*AsS);

	  // ###########################
	  // # S and P-wave decay rate #
	  // ###########################
	  myString << "(9/16 * ((2/3*FsS + 4/3*AsS*" << z->getPlotLabel() << ") * (1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") + ";
	  myString << "(1-FsS) * ";
	  myString << "(2*FlS*" << z->getPlotLabel() << "*" << z->getPlotLabel() << " * (1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") + ";
	  myString << "1/2*(1-FlS) * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") * (1+" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") + ";
	  myString << "4/3*AfbS * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") * " << y->getPlotLabel() << ")))";
	}

      if (useEffPDF == true)
	{
	  // #############################
	  // # Make 2D efficiency p.d.f. #
	  // #############################
	  myString << " * ";
	  myString << "((P0 + P1*" << z->getPlotLabel() << " + P2*" << z->getPlotLabel() << "*" << z->getPlotLabel() << " + ";
	  myString << "P3*" << z->getPlotLabel() << "*" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") + ";
	  
	  myString << "(P4 + P5*" << z->getPlotLabel() << " + P6*" << z->getPlotLabel() << "*" << z->getPlotLabel() << " + ";
	  myString << "P7*" << z->getPlotLabel() << "*" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") * ";
	  myString << y->getPlotLabel() << "*" << y->getPlotLabel() << " + ";
	  
	  myString << "(P8 + P9*" << z->getPlotLabel() << " + P10*" << z->getPlotLabel() << "*" << z->getPlotLabel() << " + ";
	  myString << "P11*" << z->getPlotLabel() << "*" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") * ";
	  myString << y->getPlotLabel() << "*" << y->getPlotLabel() << "*" << y->getPlotLabel() << " + ";
	  
	  myString << "(P12 + P13*" << z->getPlotLabel() << " + P14*" << z->getPlotLabel() << "*" << z->getPlotLabel() << " + ";
	  myString << "P15*" << z->getPlotLabel() << "*" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") * ";
	  myString << y->getPlotLabel() << "*" << y->getPlotLabel() << "*" << y->getPlotLabel() << "*" << y->getPlotLabel() << " + ";

	  myString << "(P16 + P17*" << z->getPlotLabel() << " + P18*" << z->getPlotLabel() << "*" << z->getPlotLabel() << " + ";
	  myString << "P19*" << z->getPlotLabel() << "*" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") * ";
	  myString << y->getPlotLabel() << "*" << y->getPlotLabel() << "*" << y->getPlotLabel() << "*" << y->getPlotLabel() << "*" << y->getPlotLabel() << "*" << y->getPlotLabel() << ")";

	  for (int i = 0; i < effFunc->GetNpar(); i++) Vars->add(*vecParam[i]);
	}

      Vars->add(*y);
      Vars->add(*z);


      cout << "\n@@@ 2D angular*efficiency p.d.f. @@@" << endl;
      cout << myString.str().c_str() << endl;
    }
  

  vecParam.clear();
  return myString.str(); 
}


void DeleteFit (RooAbsPdf* TotalPDF, string DeleteType)
// ##############################################################
// # DeleteType = "JustVars" --> delete just variables          #
// # DeleteType = "All"      --> delete variables and functions #
// ##############################################################
{
  stringstream myString;

  if (TotalPDF != NULL)
    {
      if ((DeleteType == "JustVars") || (DeleteType == "All"))
	{
	  if (GetVar(TotalPDF,"meanS")          != NULL) delete GetVar(TotalPDF,"meanS");
	  if (GetVar(TotalPDF,"sigmaS1")        != NULL) delete GetVar(TotalPDF,"sigmaS1");
	  if (GetVar(TotalPDF,"sigmaS2")        != NULL) delete GetVar(TotalPDF,"sigmaS2");
	  if (GetVar(TotalPDF,"fracMassS")      != NULL) delete GetVar(TotalPDF,"fracMassS");
	  if (GetVar(TotalPDF,"tau1")           != NULL) delete GetVar(TotalPDF,"tau1");
	  if (GetVar(TotalPDF,"tau2")           != NULL) delete GetVar(TotalPDF,"tau2");
	  if (GetVar(TotalPDF,"fracMassBExp")   != NULL) delete GetVar(TotalPDF,"fracMassBExp");
	  if (GetVar(TotalPDF,"meanR1")         != NULL) delete GetVar(TotalPDF,"meanR1");
	  if (GetVar(TotalPDF,"sigmaR1")        != NULL) delete GetVar(TotalPDF,"sigmaR1");
	  if (GetVar(TotalPDF,"meanR2")         != NULL) delete GetVar(TotalPDF,"meanR2");
	  if (GetVar(TotalPDF,"sigmaR2")        != NULL) delete GetVar(TotalPDF,"sigmaR2");
	  if (GetVar(TotalPDF,"fracMassBRPeak") != NULL) delete GetVar(TotalPDF,"fracMassBRPeak");
	  if (GetVar(TotalPDF,"meanL1")         != NULL) delete GetVar(TotalPDF,"meanL1");
	  if (GetVar(TotalPDF,"sigmaL1")        != NULL) delete GetVar(TotalPDF,"sigmaL1");
	  if (GetVar(TotalPDF,"meanL2")         != NULL) delete GetVar(TotalPDF,"meanL2");
	  if (GetVar(TotalPDF,"sigmaL2")        != NULL) delete GetVar(TotalPDF,"sigmaL2");
	  if (GetVar(TotalPDF,"fracMassBLPeak") != NULL) delete GetVar(TotalPDF,"fracMassBLPeak");
	  if (GetVar(TotalPDF,"fracMassBPeak")  != NULL) delete GetVar(TotalPDF,"fracMassBPeak");
	  if (GetVar(TotalPDF,"nBkgComb")       != NULL) delete GetVar(TotalPDF,"nBkgComb");
	  if (GetVar(TotalPDF,"nBkgPeak")       != NULL) delete GetVar(TotalPDF,"nBkgPeak");
	  if (GetVar(TotalPDF,"nSig")           != NULL) delete GetVar(TotalPDF,"nSig");
	  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
	    {
	      myString.clear(); myString.str("");
	      myString << "p1Poly" << i;
	      if (GetVar(TotalPDF,myString.str().c_str()) != NULL) delete GetVar(TotalPDF,myString.str().c_str());
	    } 
	  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
	    {
	      myString.clear(); myString.str("");
	      myString << "c1Poly" << i;
	      if (GetVar(TotalPDF,myString.str().c_str()) != NULL) delete GetVar(TotalPDF,myString.str().c_str());
	    }
	  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
	    {
	      myString.clear(); myString.str("");
	      myString << "p2Poly" << i;
	      if (GetVar(TotalPDF,myString.str().c_str()) != NULL) delete GetVar(TotalPDF,myString.str().c_str());
	    } 
	  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
	    {
	      myString.clear(); myString.str("");
	      myString << "c2Poly" << i;
	      if (GetVar(TotalPDF,myString.str().c_str()) != NULL) delete GetVar(TotalPDF,myString.str().c_str());
	    }
	  if (GetVar(TotalPDF,"FlS")   != NULL) delete GetVar(TotalPDF,"FlS");
	  if (GetVar(TotalPDF,"AfbS")  != NULL) delete GetVar(TotalPDF,"AfbS");
	  if (GetVar(TotalPDF,"At2S")  != NULL) delete GetVar(TotalPDF,"At2S");
	  if (GetVar(TotalPDF,"AtimS") != NULL) delete GetVar(TotalPDF,"AtimS");
	  if (GetVar(TotalPDF,"FsS")   != NULL) delete GetVar(TotalPDF,"FsS");
	  if (GetVar(TotalPDF,"AsS")   != NULL) delete GetVar(TotalPDF,"AsS");
	}
      else if (DeleteType == "All")
	{
	  if (TotalPDF->getComponents()->find("MassS1")           != NULL) delete TotalPDF->getComponents()->find("MassS1");
	  if (TotalPDF->getComponents()->find("MassS2")           != NULL) delete TotalPDF->getComponents()->find("MassS2");
	  if (TotalPDF->getComponents()->find("MassSignal")       != NULL) delete TotalPDF->getComponents()->find("MassSignal");
	  if (TotalPDF->getComponents()->find("AngleS")           != NULL) delete TotalPDF->getComponents()->find("AngleS");
	  if (TotalPDF->getComponents()->find("Efficency")        != NULL) delete TotalPDF->getComponents()->find("Efficiency");
	  if (TotalPDF->getComponents()->find("Signal")           != NULL) delete TotalPDF->getComponents()->find("Signal");
	  if (TotalPDF->getComponents()->find("BkgMassExp1")      != NULL) delete TotalPDF->getComponents()->find("BkgMassExp1");
	  if (TotalPDF->getComponents()->find("BkgMassExp2")      != NULL) delete TotalPDF->getComponents()->find("BkgMassExp2");
	  if (TotalPDF->getComponents()->find("BkgMassRPeak1")    != NULL) delete TotalPDF->getComponents()->find("BkgMassRPeak1");
	  if (TotalPDF->getComponents()->find("BkgMassRPeak2")    != NULL) delete TotalPDF->getComponents()->find("BkgMassRPeak1");
	  if (TotalPDF->getComponents()->find("BkgMassRPeak")     != NULL) delete TotalPDF->getComponents()->find("BkgMassRPeak1");
	  if (TotalPDF->getComponents()->find("BkgMassLPeak1")    != NULL) delete TotalPDF->getComponents()->find("BkgMassLPeak1");
	  if (TotalPDF->getComponents()->find("BkgMassLPeak2")    != NULL) delete TotalPDF->getComponents()->find("BkgMassLPeak1");
	  if (TotalPDF->getComponents()->find("BkgMassLPeak")     != NULL) delete TotalPDF->getComponents()->find("BkgMassLPeak1");
	  if (TotalPDF->getComponents()->find("BkgMassComb")      != NULL) delete TotalPDF->getComponents()->find("BkgMassComb");
	  if (TotalPDF->getComponents()->find("BkgMassPeak")      != NULL) delete TotalPDF->getComponents()->find("BkgMassPeak");
	  if (TotalPDF->getComponents()->find("BkgAngle1")        != NULL) delete TotalPDF->getComponents()->find("BkgAngle1");
	  if (TotalPDF->getComponents()->find("BkgAngle2")        != NULL) delete TotalPDF->getComponents()->find("BkgAngle2");
	  if (TotalPDF->getComponents()->find("BkgAnglesC")       != NULL) delete TotalPDF->getComponents()->find("BkgAnglesC");
	  if (TotalPDF->getComponents()->find("BkgAnglesP")       != NULL) delete TotalPDF->getComponents()->find("BkgAnglesP");
	  if (TotalPDF->getComponents()->find("BkgMassAngleComb") != NULL) delete TotalPDF->getComponents()->find("BkgMassAngleComb");
	  if (TotalPDF->getComponents()->find("BkgMassAnglePeak") != NULL) delete TotalPDF->getComponents()->find("BkgMassAnglePeak");

	  if (TotalPDF != NULL) delete TotalPDF;
	}
      else
	{
	  cout << "[ExtractYield::DeleteFit]\tWrong parameter: " << DeleteType << endl;
	  exit(1);
	}
    }
}


void ResetCombPolyParam (vector<vector<string>*>* fitParam, vector<double>* q2Bins)
{
  stringstream myCoeff;

  for (unsigned int j = 0; j < fitParam->operator[](0)->size(); j++)
    {
      cout << "\n@@@ Resetting the combinatorial background polynomial coefficients for q^2 bin #" << j << " @@@" << endl;

      myCoeff.clear(); myCoeff.str("");
      myCoeff << NCOEFFPOLYBKG;
      fitParam->operator[](Utility->GetFitParamIndx("nPolyC1"))->operator[](j) = myCoeff.str();
      cout << "Degree set to: " << fitParam->operator[](Utility->GetFitParamIndx("nPolyC1"))->operator[](j).c_str() << endl;
      for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
	{
	  myCoeff.clear(); myCoeff.str("");
	  myCoeff << "c1Poly" << i;
	  fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](j) = "0.0";

	  cout << myCoeff.str().c_str() << "\t" << fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](j).c_str() << endl;
	}

      myCoeff.clear(); myCoeff.str("");
      myCoeff << NCOEFFPOLYBKG;
      fitParam->operator[](Utility->GetFitParamIndx("nPolyC2"))->operator[](j) = myCoeff.str();
      cout << "Degree set to: " << fitParam->operator[](Utility->GetFitParamIndx("nPolyC2"))->operator[](j).c_str() << endl;
      for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
	{
	  myCoeff.clear(); myCoeff.str("");
	  myCoeff << "c2Poly" << i;
	  fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](j) = "0.0";

	  cout << myCoeff.str().c_str() << "\t" << fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](j).c_str() << endl;
	}
    }
}


void ResetAngularParam (vector<vector<string>*>* fitParam, vector<double>* q2Bins)
{
  for (unsigned int j = 0; j < fitParam->operator[](0)->size(); j++)
    {
      cout << "\n@@@ Resetting the angular parameters for q^2 bin #" << j << " @@@" << endl;

      fitParam->operator[](Utility->GetFitParamIndx("FlS"))->operator[](j)   = "0.5   0.0   0.0";
      fitParam->operator[](Utility->GetFitParamIndx("AfbS"))->operator[](j)  = "0.0   0.0   0.0";
      fitParam->operator[](Utility->GetFitParamIndx("At2S"))->operator[](j)  = "0.0   0.0   0.0";
      fitParam->operator[](Utility->GetFitParamIndx("AtimS"))->operator[](j) = "0.0   0.0   0.0";
      fitParam->operator[](Utility->GetFitParamIndx("FsS"))->operator[](j)   = "0.05   0.02   0.02";
      fitParam->operator[](Utility->GetFitParamIndx("AsS"))->operator[](j)   = "0.0   0.0   0.0";

      cout << "FL: "   << "\t" << fitParam->operator[](Utility->GetFitParamIndx("FlS"))->operator[](j).c_str() << endl;
      cout << "AFB: "  << "\t" << fitParam->operator[](Utility->GetFitParamIndx("AfbS"))->operator[](j).c_str() << endl;
      cout << "AT2: "  << "\t" << fitParam->operator[](Utility->GetFitParamIndx("At2S"))->operator[](j).c_str() << endl;
      cout << "ATIM: " << "\t" << fitParam->operator[](Utility->GetFitParamIndx("AtimS"))->operator[](j).c_str() << endl;
      cout << "FS: "   << "\t" << fitParam->operator[](Utility->GetFitParamIndx("FsS"))->operator[](j).c_str() << endl;
      cout << "AS: "   << "\t" << fitParam->operator[](Utility->GetFitParamIndx("AsS"))->operator[](j).c_str() << endl;
    }
}


double StoreFitResultsInFile (RooAbsPdf** TotalPDF, RooFitResult* fitResult, RooDataSet* dataSet, RooArgSet* vecConstr)
{
  RooAbsReal* NLL;
  double signalSigma  = 0.0;
  double signalSigmaE = 0.0;
  double NLLvalue = NULL;

  if (fitResult != NULL)
    {
      fileFitResults << "Covariance quality (3=ok): " << fitResult->covQual() << endl;
      fileFitResults << "Fit status (0=ok): " << fitResult->status() << endl;
      fileFitResults << "Fit EDM: " << fitResult->edm() << endl;
      fileFitResults << "FCN(mPDG(B0),ang=0): " << (*TotalPDF)->getVal() << endl;
      if (ApplyConstr == true) NLL = (*TotalPDF)->createNLL(*dataSet,Extended(true),ExternalConstraints(*vecConstr));
      else                     NLL = (*TotalPDF)->createNLL(*dataSet,Extended(true));
      NLLvalue = NLL->getVal();
      delete NLL;
      fileFitResults << "NLL: " << NLLvalue << endl;


      if (GetVar(*TotalPDF,"fracMassS") != NULL)
	{
	  signalSigma  = sqrt((*TotalPDF)->getVariables()->getRealValue("fracMassS") * pow((*TotalPDF)->getVariables()->getRealValue("sigmaS1"),2.) +
			      (1.-(*TotalPDF)->getVariables()->getRealValue("fracMassS")) * pow((*TotalPDF)->getVariables()->getRealValue("sigmaS2"),2.));
	  signalSigmaE = 1./(2.*signalSigma) * sqrt( pow((pow((*TotalPDF)->getVariables()->getRealValue("sigmaS1"),2.)-pow((*TotalPDF)->getVariables()->getRealValue("sigmaS2"),2.)) * GetVar(*TotalPDF,"fracMassS")->getError(),2.) +
						     pow(2.*(*TotalPDF)->getVariables()->getRealValue("fracMassS") * (*TotalPDF)->getVariables()->getRealValue("sigmaS1")*GetVar(*TotalPDF,"sigmaS1")->getError(),2.) +
						     pow(2.*(1.-(*TotalPDF)->getVariables()->getRealValue("fracMassS")) * (*TotalPDF)->getVariables()->getRealValue("sigmaS2")*GetVar(*TotalPDF,"sigmaS2")->getError(),2.) );
	}
      else if (GetVar(*TotalPDF,"sigmaS1") != NULL)
	{
	  signalSigma  = (*TotalPDF)->getVariables()->getRealValue("sigmaS1");
	  signalSigmaE = GetVar(*TotalPDF,"sigmaS1")->getError();
	}
  
      if (GetVar(*TotalPDF,"meanS") != NULL)
	{
	  fileFitResults << "Mean: " << GetVar(*TotalPDF,"meanS")->getVal() << " +/- " << GetVar(*TotalPDF,"meanS")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"meanS")->getErrorHi() << "/" << GetVar(*TotalPDF,"meanS")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"sigmaS1") != NULL)
	{
	  fileFitResults << "Sigma-1: " << GetVar(*TotalPDF,"sigmaS1")->getVal() << " +/- " << GetVar(*TotalPDF,"sigmaS1")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"sigmaS1")->getErrorHi() << "/" << GetVar(*TotalPDF,"sigmaS1")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"sigmaS2") != NULL)
	{
	  fileFitResults << "Sigma-2: " << GetVar(*TotalPDF,"sigmaS2")->getVal() << " +/- " << GetVar(*TotalPDF,"sigmaS2")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"sigmaS2")->getErrorHi() << "/" << GetVar(*TotalPDF,"sigmaS2")->getErrorLo() << ")" << endl;
	  fileFitResults << "Fraction: " << GetVar(*TotalPDF,"fracMassS")->getVal() << " +/- " << GetVar(*TotalPDF,"fracMassS")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"fracMassS")->getErrorHi() << "/" << GetVar(*TotalPDF,"fracMassS")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"sigmaS1") != NULL) fileFitResults << "< Sigma >: " << signalSigma << " +/- " << signalSigmaE << endl;
      if (GetVar(*TotalPDF,"nSig") != NULL)
	{
	  fileFitResults << "Signal yield: " << GetVar(*TotalPDF,"nSig")->getVal() << " +/- " << GetVar(*TotalPDF,"nSig")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"nSig")->getErrorHi() << "/" << GetVar(*TotalPDF,"nSig")->getErrorLo() << ")" << endl;
	}
      

      if (GetVar(*TotalPDF,"tau1") != NULL)
	{
	  fileFitResults << "Background mass tau-1: " << GetVar(*TotalPDF,"tau1")->getVal() << " +/- " << tau1->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"tau1")->getErrorHi() << "/" << GetVar(*TotalPDF,"tau1")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"tau2") != NULL)
	{
	  fileFitResults << "Background mass tau-2: " << GetVar(*TotalPDF,"tau2")->getVal() << " +/- " << GetVar(*TotalPDF,"tau2")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"tau2")->getErrorHi() << "/" << GetVar(*TotalPDF,"tau2")->getErrorLo() << ")" << endl;
	  fileFitResults << "Fraction: " << GetVar(*TotalPDF,"fracMassBExp")->getVal() << " +/- " << GetVar(*TotalPDF,"fracMassBExp")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"fracMassBExp")->getErrorHi() << "/" << GetVar(*TotalPDF,"fracMassBExp")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"nBkgComb") != NULL)
	{
	  fileFitResults << "Comb. bkg. yield: " << GetVar(*TotalPDF,"nBkgComb")->getVal() << " +/- " << GetVar(*TotalPDF,"nBkgComb")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"nBkgComb")->getErrorHi() << "/" << GetVar(*TotalPDF,"nBkgComb")->getErrorLo() << ")" << endl;
	}


      if (GetVar(*TotalPDF,"meanR1") != NULL)
	{
	  fileFitResults << "Background right-peak mean-1: " << GetVar(*TotalPDF,"meanR1")->getVal() << " +/- " << GetVar(*TotalPDF,"meanR1")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"meanR1")->getErrorHi() << "/" << GetVar(*TotalPDF,"meanR1")->getErrorLo() << ")" << endl;
	  fileFitResults << "Background right-peak sigma-1: " << GetVar(*TotalPDF,"sigmaR1")->getVal() << " +/- " << GetVar(*TotalPDF,"sigmaR1")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"sigmaR1")->getErrorHi() << "/" << GetVar(*TotalPDF,"sigmaR1")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"meanR2") != NULL)
	{
	  fileFitResults << "Background right-peak mean-2: " << GetVar(*TotalPDF,"meanR2")->getVal() << " +/- " << GetVar(*TotalPDF,"meanR2")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"meanR2")->getErrorHi() << "/" << GetVar(*TotalPDF,"meanR2")->getErrorLo() << ")" << endl;
	  fileFitResults << "Background right-peak sigma-2: " << GetVar(*TotalPDF,"sigmaR2")->getVal() << " +/- " << GetVar(*TotalPDF,"sigmaR2")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"sigmaR2")->getErrorHi() << "/" << GetVar(*TotalPDF,"sigmaR2")->getErrorLo() << ")" << endl;
	  fileFitResults << "Fraction: " << GetVar(*TotalPDF,"fracMassBRPeak")->getVal() << " +/- " << GetVar(*TotalPDF,"fracMassBRPeak")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"fracMassBRPeak")->getErrorHi() << "/" << GetVar(*TotalPDF,"fracMassBRPeak")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"meanL1") != NULL)
	{
	  fileFitResults << "Background left-peak mean-1: " << GetVar(*TotalPDF,"meanL1")->getVal() << " +/- " << GetVar(*TotalPDF,"meanL1")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"meanL1")->getErrorHi() << "/" << GetVar(*TotalPDF,"meanL1")->getErrorLo() << ")" << endl;
	  fileFitResults << "Background left-peak sigma-1: " << GetVar(*TotalPDF,"sigmaL1")->getVal() << " +/- " << GetVar(*TotalPDF,"sigmaL1")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"sigmaL1")->getErrorHi() << "/" << GetVar(*TotalPDF,"sigmaL1")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"meanL2") != NULL)
	{
	  fileFitResults << "Background left-peak mean-2: " << GetVar(*TotalPDF,"meanL2")->getVal() << " +/- " << GetVar(*TotalPDF,"meanL2")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"meanL2")->getErrorHi() << "/" << GetVar(*TotalPDF,"meanL2")->getErrorLo() << ")" << endl;
	  fileFitResults << "Background left-peak sigma-2: " << GetVar(*TotalPDF,"sigmaL2")->getVal() << " +/- " << GetVar(*TotalPDF,"sigmaL2")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"sigmaL2")->getErrorHi() << "/" << GetVar(*TotalPDF,"sigmaL2")->getErrorLo() << ")" << endl;
	  fileFitResults << "Fraction: " << GetVar(*TotalPDF,"fracMassBLPeak")->getVal() << " +/- " << GetVar(*TotalPDF,"fracMassBLPeak")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"fracMassBLPeak")->getErrorHi() << "/" << GetVar(*TotalPDF,"fracMassBLPeak")->getErrorLo() << ")" << endl;
	} 
      if (GetVar(*TotalPDF,"fracMassBPeak") != NULL)
	{
	  fileFitResults << "Fraction right-left peak: " << GetVar(*TotalPDF,"fracMassBPeak")->getVal() << " +/- " << GetVar(*TotalPDF,"fracMassBPeak")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"fracMassBPeak")->getErrorHi() << "/" << GetVar(*TotalPDF,"fracMassBPeak")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"nBkgPeak") != NULL)
	{
	  fileFitResults << "Peaking bkg. yield: " << GetVar(*TotalPDF,"nBkgPeak")->getVal() << " +/- " << GetVar(*TotalPDF,"nBkgPeak")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"nBkgPeak")->getErrorHi() << "/" << GetVar(*TotalPDF,"nBkgPeak")->getErrorLo() << ")" << endl;
	}


      if (GetVar(*TotalPDF,"FlS") != NULL)
	{
	  fileFitResults << "Fl: " << GetVar(*TotalPDF,"FlS")->getVal() << " +/- " << GetVar(*TotalPDF,"FlS")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"FlS")->getErrorHi() << "/" << GetVar(*TotalPDF,"FlS")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"AfbS") != NULL)
	{
	  fileFitResults << "Afb: " << GetVar(*TotalPDF,"AfbS")->getVal() << " +/- ";
	  fileFitResults << GetVar(*TotalPDF,"AfbS")->getError() << " (" << GetVar(*TotalPDF,"AfbS")->getErrorHi() << "/" << GetVar(*TotalPDF,"AfbS")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"At2S") != NULL)
	{
	  fileFitResults << "At2: " << GetVar(*TotalPDF,"At2S")->getVal() << " +/- " << GetVar(*TotalPDF,"At2S")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"At2S")->getErrorHi() << "/" <<  GetVar(*TotalPDF,"At2S")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"AtimS") != NULL)
	{
	  fileFitResults << "Atim: " << GetVar(*TotalPDF,"AtimS")->getVal() << " +/- " << GetVar(*TotalPDF,"AtimS")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"AtimS")->getErrorHi() << "/" <<  GetVar(*TotalPDF,"AtimS")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"FsS") != NULL)
	{
	  fileFitResults << "Fs: " << GetVar(*TotalPDF,"FsS")->getVal() << " +/- " << GetVar(*TotalPDF,"FsS")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"FsS")->getErrorHi() << "/" <<  GetVar(*TotalPDF,"FsS")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"AsS") != NULL)
	{
	  fileFitResults << "As: " << GetVar(*TotalPDF,"AsS")->getVal() << " +/- " << GetVar(*TotalPDF,"AsS")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"AsS")->getErrorHi() << "/" <<  GetVar(*TotalPDF,"AsS")->getErrorLo() << ")" << endl;
	}
    }


  return NLLvalue;
}


void StorePolyResultsInFile (RooAbsPdf** TotalPDF)
{
  stringstream myString;

  if (GetVar(*TotalPDF,"nBkgComb") != NULL)
    {
      for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
	{
	  myString.clear(); myString.str("");
	  myString << "c1Poly" << i;
	  if (GetVar(*TotalPDF,myString.str().c_str()) != NULL)
	    {
	      fileFitResults << "Comb.bkg angle-1 c" << i << ": " << GetVar(*TotalPDF,myString.str().c_str())->getVal() << " +/- " << GetVar(*TotalPDF,myString.str().c_str())->getError();
	      fileFitResults << " (" << GetVar(*TotalPDF,myString.str().c_str())->getErrorHi() << "/" << GetVar(*TotalPDF,myString.str().c_str())->getErrorLo() << ")" << endl;
	    }
	}
	  
      for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
	{
	  myString.clear(); myString.str("");
	  myString << "c2Poly" << i;
	  if (GetVar(*TotalPDF,myString.str().c_str()) != NULL)
	    {
	      fileFitResults << "Comb.bkg angle-2 c" << i << ": " << GetVar(*TotalPDF,myString.str().c_str())->getVal() << " +/- " << GetVar(*TotalPDF,myString.str().c_str())->getError();
	      fileFitResults << " (" << GetVar(*TotalPDF,myString.str().c_str())->getErrorHi() << "/" << GetVar(*TotalPDF,myString.str().c_str())->getErrorLo() << ")" << endl;
	    }
	}
    }

  if (GetVar(*TotalPDF,"nBkgPeak") != NULL)
    {
      for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
	{
	  myString.clear(); myString.str("");
	  myString << "p1Poly" << i;
	  if (GetVar(*TotalPDF,myString.str().c_str()) != NULL)
	    {
	      fileFitResults << "Peak.bkg angle-1 p" << i << ": " << GetVar(*TotalPDF,myString.str().c_str())->getVal() << " +/- " << GetVar(*TotalPDF,myString.str().c_str())->getError();
	      fileFitResults << " (" << GetVar(*TotalPDF,myString.str().c_str())->getErrorHi() << "/" << GetVar(*TotalPDF,myString.str().c_str())->getErrorHi() << ")" << endl;
	    }
	}
	  
      for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
	{
	  myString.clear(); myString.str("");
	  myString << "p2Poly" << i;
	  if (GetVar(*TotalPDF,myString.str().c_str()) != NULL)
	    {
	      fileFitResults << "Peak.bkg angle-2 p" << i << ": " << GetVar(*TotalPDF,myString.str().c_str())->getVal() << " +/- " << GetVar(*TotalPDF,myString.str().c_str())->getError();
	      fileFitResults << " (" << GetVar(*TotalPDF,myString.str().c_str())->getErrorHi() << "/" << GetVar(*TotalPDF,myString.str().c_str())->getErrorHi() << ")" << endl;
	    }
	}
    }
}


vector<string>* SaveFitResults (RooAbsPdf* TotalPDF, unsigned int fitParamIndx, vector<vector<string>*>* fitParam, vector<vector<unsigned int>*>* configParam, RooArgSet* vecConstr)
{
  stringstream myString;
  stringstream myCoeff;
  vector<string>* vecParStr;
  unsigned int NCoeffPolyBKGpeak1 = 0;
  unsigned int NCoeffPolyBKGcomb1 = 0;
  unsigned int NCoeffPolyBKGpeak2 = 0;
  unsigned int NCoeffPolyBKGcomb2 = 0;

  vecParStr = new vector<string>;


  vecParStr->push_back("# Signal mean [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"meanS") != NULL)) { myString.clear(); myString.str(""); myString << TotalPDF->getVariables()->getRealValue("meanS"); vecParStr->push_back(myString.str()); }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("meanS"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# Signal sigma-1 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"sigmaS1") != NULL) && ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(string("sigmaS1") + string("_constr")).c_str()) == NULL))))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("sigmaS1") << "   " << GetVar(TotalPDF,"sigmaS1")->getErrorLo() << "   " << GetVar(TotalPDF,"sigmaS1")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# Signal sigma-2 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"sigmaS2") != NULL) && ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(string("sigmaS2") + string("_constr")).c_str()) == NULL))))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("sigmaS2") << "   " << GetVar(TotalPDF,"sigmaS2")->getErrorLo() << "   " << GetVar(TotalPDF,"sigmaS2")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("sigmaS2"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# Fraction");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"fracMassS") != NULL) && ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(string("fracMassS") + string("_constr")).c_str()) == NULL))))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("fracMassS") << "   " << GetVar(TotalPDF,"fracMassS")->getErrorLo() << "   " << GetVar(TotalPDF,"fracMassS")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("fracMassS"))->operator[](fitParamIndx).c_str());


  vecParStr->push_back("###################################################################");


  vecParStr->push_back("# Bkg tau-1 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"tau1")         != NULL)) { myString.clear(); myString.str(""); myString << TotalPDF->getVariables()->getRealValue("tau1"); vecParStr->push_back(myString.str()); }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("tau1"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# Bkg tau-2 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"tau2")         != NULL)) { myString.clear(); myString.str(""); myString << TotalPDF->getVariables()->getRealValue("tau2"); vecParStr->push_back(myString.str()); }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("tau2"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# Fraction");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"fracMassBExp") != NULL)) { myString.clear(); myString.str(""); myString << TotalPDF->getVariables()->getRealValue("fracMassBExp"); vecParStr->push_back(myString.str()); }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("fracMassBExp"))->operator[](fitParamIndx).c_str());


  vecParStr->push_back("###################################################################");


  vecParStr->push_back("# Bkg mean right-peak-1 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"meanR1") != NULL) && ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(string("meanR1") + string("_constr")).c_str()) == NULL))))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("meanR1") << "   " << GetVar(TotalPDF,"meanR1")->getErrorLo() << "   " << GetVar(TotalPDF,"meanR1")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("meanR1"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# Bkg sigma right-peak-1 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"sigmaR1") != NULL) && ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(string("sigmaR1") + string("_constr")).c_str()) == NULL))))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("sigmaR1") << "   " << GetVar(TotalPDF,"sigmaR1")->getErrorLo() << "   " << GetVar(TotalPDF,"sigmaR1")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# Bkg mean right-peak-2 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"meanR2") != NULL) && ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(string("meanR2") + string("_constr")).c_str()) == NULL))))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("meanR2") << "   " << GetVar(TotalPDF,"meanR2")->getErrorLo() << "   " << GetVar(TotalPDF,"meanR2")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("meanR2"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# Bkg sigma right-peak-2 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"sigmaR2") != NULL) && ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(string("sigmaR2") + string("_constr")).c_str()) == NULL))))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("sigmaR2") << "   " << GetVar(TotalPDF,"sigmaR2")->getErrorLo() << "   " << GetVar(TotalPDF,"sigmaR2")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# Fraction");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"fracMassBRPeak") != NULL) && ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(string("fracMassBRPeak") + string("_constr")).c_str()) == NULL))))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("fracMassBRPeak") << "   " << GetVar(TotalPDF,"fracMassBRPeak")->getErrorLo() << "   " << GetVar(TotalPDF,"fracMassBRPeak")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("fracMassBRPeak"))->operator[](fitParamIndx).c_str());


  vecParStr->push_back("# Bkg mean left-peak-1 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"meanL1") != NULL) && ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(string("meanL1") + string("_constr")).c_str()) == NULL))))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("meanL1") << "   " << GetVar(TotalPDF,"meanL1")->getErrorLo() << "   " << GetVar(TotalPDF,"meanL1")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("meanL1"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# Bkg sigma left-peak-1 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"sigmaL1") != NULL) && ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(string("sigmaL1") + string("_constr")).c_str()) == NULL))))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("sigmaL1") << "   " << GetVar(TotalPDF,"sigmaL1")->getErrorLo() << "   " << GetVar(TotalPDF,"sigmaL1")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# Bkg mean left-peak-2 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"meanL2") != NULL) && ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(string("meanL2") + string("_constr")).c_str()) == NULL))))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("meanL2") << "   " << GetVar(TotalPDF,"meanL2")->getErrorLo() << "   " << GetVar(TotalPDF,"meanL2")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("meanL2"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# Bkg sigma left-peak-2 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"sigmaL2") != NULL) && ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(string("sigmaL2") + string("_constr")).c_str()) == NULL))))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("sigmaL2") << "   " << GetVar(TotalPDF,"sigmaL2")->getErrorLo() << "   " << GetVar(TotalPDF,"sigmaL2")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# Fraction");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"fracMassBLPeak") != NULL) && ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(string("fracMassBLPeak") + string("_constr")).c_str()) == NULL))))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("fracMassBLPeak") << "   " << GetVar(TotalPDF,"fracMassBLPeak")->getErrorLo() << "   " << GetVar(TotalPDF,"fracMassBLPeak")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("fracMassBLPeak"))->operator[](fitParamIndx).c_str());

  vecParStr->push_back("# Fraction of right-peak");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"fracMassBPeak") != NULL) && ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(string("fracMassBPeak") + string("_constr")).c_str()) == NULL))))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("fracMassBPeak") << "   " << GetVar(TotalPDF,"fracMassBPeak")->getErrorLo() << "   " << GetVar(TotalPDF,"fracMassBPeak")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("fracMassBPeak"))->operator[](fitParamIndx).c_str());


  vecParStr->push_back("###################################################################");


  vecParStr->push_back("# Combinatorial bkg yield");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"nBkgComb") != NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("nBkgComb") << "   " << GetVar(TotalPDF,"nBkgComb")->getErrorLo() << "   " << GetVar(TotalPDF,"nBkgComb")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("nBkgComb"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# Peaking bkg yield");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"nBkgPeak") != NULL) && ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(string("nBkgPeak") + string("_constr")).c_str()) == NULL))))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("nBkgPeak") << "   " << GetVar(TotalPDF,"nBkgPeak")->getErrorLo() << "   " << GetVar(TotalPDF,"nBkgPeak")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("nBkgPeak"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# Signal yield");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"nSig") != NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("nSig") << "   " << GetVar(TotalPDF,"nSig")->getErrorLo() << "   " << GetVar(TotalPDF,"nSig")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("nSig"))->operator[](fitParamIndx).c_str());


  vecParStr->push_back("###################################################################");

  NCoeffPolyBKGpeak1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP1"))->operator[](fitParamIndx).c_str());
  NCoeffPolyBKGcomb1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC1"))->operator[](fitParamIndx).c_str());
  NCoeffPolyBKGpeak2 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP2"))->operator[](fitParamIndx).c_str());
  NCoeffPolyBKGcomb2 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC2"))->operator[](fitParamIndx).c_str());

  vecParStr->push_back("### Number of coefficients of the peak.bkg angle1 poly ###");
  myString.clear(); myString.str("");
  myString << NCoeffPolyBKGpeak1;
  vecParStr->push_back(myString.str());
  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
    {
      myCoeff.clear(); myCoeff.str("");
      myCoeff << "p1Poly" << i;
      myString.clear(); myString.str("");
      myString << "# Bkg angle1 poly coeff. p" << i;
      vecParStr->push_back(myString.str());
      myString.clear(); myString.str("");
      if ((SAVEPOLY == true) &&
	  (TotalPDF != NULL) &&
	  (GetVar(TotalPDF,myCoeff.str().c_str()) != NULL) &&
	  ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(myCoeff.str() + string("_constr")).c_str()) == NULL))))
	{
	  myString << TotalPDF->getVariables()->getRealValue(myCoeff.str().c_str()) << "   " << GetVar(TotalPDF,myCoeff.str().c_str())->getErrorLo() << "   " << GetVar(TotalPDF,myCoeff.str().c_str())->getErrorHi();
	  vecParStr->push_back(myString.str());
	}
      else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](fitParamIndx).c_str());
    }

  vecParStr->push_back("### Number of coefficients of the comb.bkg angle1 poly ###");
  myString.clear(); myString.str("");
  myString << NCoeffPolyBKGcomb1;
  vecParStr->push_back(myString.str());
  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
    {
      myCoeff.clear(); myCoeff.str("");
      myCoeff << "c1Poly" << i;
      myString.clear(); myString.str("");
      myString << "# Bkg angle1 poly coeff. c" << i;
      vecParStr->push_back(myString.str());
      myString.clear(); myString.str("");
      if ((SAVEPOLY == true) &&
	  (TotalPDF != NULL) &&
	  (GetVar(TotalPDF,myCoeff.str().c_str()) != NULL))
	{
	  myString << TotalPDF->getVariables()->getRealValue(myCoeff.str().c_str()) << "   " << GetVar(TotalPDF,myCoeff.str().c_str())->getErrorLo() << "   " << GetVar(TotalPDF,myCoeff.str().c_str())->getErrorHi();
	  vecParStr->push_back(myString.str());
	}
      else vecParStr->push_back("0.0");
    }
      
  vecParStr->push_back("### Number of coefficients of the peak.bkg angle2 poly ###");
  myString.clear(); myString.str("");
  myString << NCoeffPolyBKGpeak2;
  vecParStr->push_back(myString.str());
  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
    {
      myCoeff.clear(); myCoeff.str("");
      myCoeff << "p2Poly" << i;
      myString.clear(); myString.str("");
      myString << "# Bkg angle2 poly coeff. p" << i;
      vecParStr->push_back(myString.str());
      myString.clear(); myString.str("");
      if ((SAVEPOLY == true) &&
	  (TotalPDF != NULL) &&
	  (GetVar(TotalPDF,myCoeff.str().c_str()) != NULL) &&
	  ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(myCoeff.str() + string("_constr")).c_str()) == NULL))))
	{
	  myString << TotalPDF->getVariables()->getRealValue(myCoeff.str().c_str()) << "   " << GetVar(TotalPDF,myCoeff.str().c_str())->getErrorLo() << "   " << GetVar(TotalPDF,myCoeff.str().c_str())->getErrorHi();
	  vecParStr->push_back(myString.str());
	}
      else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](fitParamIndx).c_str());
    }

  vecParStr->push_back("### Number of coefficients of the comb.bkg angle2 poly ###");
  myString.clear(); myString.str("");
  myString << NCoeffPolyBKGcomb2;
  vecParStr->push_back(myString.str());
  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
    {
      myCoeff.clear(); myCoeff.str("");
      myCoeff << "c2Poly" << i;
      myString.clear(); myString.str("");
      myString << "# Bkg angle2 poly coeff. c" << i;
      vecParStr->push_back(myString.str());
      myString.clear(); myString.str("");
      if ((SAVEPOLY == true) &&
	  (TotalPDF != NULL) &&
	  (GetVar(TotalPDF,myCoeff.str().c_str()) != NULL))
	{
	  myString << TotalPDF->getVariables()->getRealValue(myCoeff.str().c_str()) << "   " << GetVar(TotalPDF,myCoeff.str().c_str())->getErrorLo() << "   " << GetVar(TotalPDF,myCoeff.str().c_str())->getErrorHi();
	  vecParStr->push_back(myString.str());
	}
      else vecParStr->push_back("0.0");
    }


  vecParStr->push_back("###################################################################");


  vecParStr->push_back("# FL +/- err");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"FlS") != NULL) && ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(string("FlS") + string("_constr")).c_str()) == NULL))))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("FlS") << "   " << GetVar(TotalPDF,"FlS")->getErrorLo() << "   " << GetVar(TotalPDF,"FlS")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("FlS"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# AFB +/- err");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"AfbS") != NULL) && ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(string("AfbS") + string("_constr")).c_str()) == NULL))))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("AfbS") << "   " << GetVar(TotalPDF,"AfbS")->getErrorLo() << "   " << GetVar(TotalPDF,"AfbS")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("AfbS"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# AT2 +/- err");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"At2S") != NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("At2S") << "   " << GetVar(TotalPDF,"At2S")->getErrorLo() << "   " << GetVar(TotalPDF,"At2S")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("At2S"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# ATIM +/- err");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"AtimS") != NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("AtimS") << "   " << GetVar(TotalPDF,"AtimS")->getErrorLo() << "   " << GetVar(TotalPDF,"AtimS")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("AtimS"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# FS +/- err");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"FsS") != NULL) && ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(string("FsS") + string("_constr")).c_str()) == NULL))))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("FsS") << "   " << GetVar(TotalPDF,"FsS")->getErrorLo() << "   " << GetVar(TotalPDF,"FsS")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("FsS"))->operator[](fitParamIndx).c_str());
  vecParStr->push_back("# AS +/- err");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"AsS") != NULL) && ((ApplyConstr == false) || ((ApplyConstr == true) && (vecConstr->find(string(string("AsS") + string("_constr")).c_str()) == NULL))))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("AsS") << "   " << GetVar(TotalPDF,"AsS")->getErrorLo() << "   " << GetVar(TotalPDF,"AsS")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("AsS"))->operator[](fitParamIndx).c_str());


  vecParStr->push_back("###################################################################");


  vecParStr->push_back("# Psi intrusion bkg: 0 = no; 1 = one gaussian; 2 = two gaussians");
  myString.clear(); myString.str("");
  myString << configParam->operator[](Utility->GetConfigParamIndx("nGaussPSIintru"))->operator[](fitParamIndx);
  vecParStr->push_back(myString.str());

  vecParStr->push_back("# Intrusion bkg: 0 = don't exist; 1 = exist 1 peak alone; 2 = exist 1 peak with rest; 3 = exist 2 peaks with rest");
  myString.clear(); myString.str("");
  myString << configParam->operator[](Utility->GetConfigParamIndx("fitPSIintru"))->operator[](fitParamIndx);
  vecParStr->push_back(myString.str());

  vecParStr->push_back("# Use two gaussians signal: 0 = no; 1 = yes");
  myString.clear(); myString.str("");
  myString << configParam->operator[](Utility->GetConfigParamIndx("SigType"))->operator[](fitParamIndx);
  vecParStr->push_back(myString.str());

  vecParStr->push_back("# Use two exponentials bkg: 0 = no; 1 = yes");
  myString.clear(); myString.str("");
  myString << configParam->operator[](Utility->GetConfigParamIndx("BkgType"))->operator[](fitParamIndx);
  vecParStr->push_back(myString.str());

  return vecParStr;
}


unsigned int CopyFitResults (RooAbsPdf* TotalPDF, unsigned int fitParamIndx, vector<vector<string>*>* fitParam)
{
  stringstream myString;
  stringstream myCoeff;
  double value, errLo, errHi;
  unsigned int NCoeffPolyBKGpeak1;
  unsigned int NCoeffPolyBKGcomb1;
  unsigned int NCoeffPolyBKGpeak2;
  unsigned int NCoeffPolyBKGcomb2;


  if (GetVar(TotalPDF,"meanS") != NULL) TotalPDF->getVariables()->setRealValue("meanS",atof(fitParam->operator[](Utility->GetFitParamIndx("meanS"))->operator[](fitParamIndx).c_str()));
  if (GetVar(TotalPDF,"sigmaS1") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,"sigmaS1",1.0,&myString,&value,&errLo,&errHi);
    }
  if (GetVar(TotalPDF,"sigmaS2") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaS2"))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,"sigmaS2",1.0,&myString,&value,&errLo,&errHi);
    }
  if (GetVar(TotalPDF,"fracMassS") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("fracMassS"))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,"fracMassS",1.0,&myString,&value,&errLo,&errHi);
    }

  if (GetVar(TotalPDF,"tau1")         != NULL) TotalPDF->getVariables()->setRealValue("tau1",         atof(fitParam->operator[](Utility->GetFitParamIndx("tau1"))->operator[](fitParamIndx).c_str()));
  if (GetVar(TotalPDF,"tau2")         != NULL) TotalPDF->getVariables()->setRealValue("tau2",         atof(fitParam->operator[](Utility->GetFitParamIndx("tau2"))->operator[](fitParamIndx).c_str()));
  if (GetVar(TotalPDF,"fracMassBExp") != NULL) TotalPDF->getVariables()->setRealValue("fracMassBExp", atof(fitParam->operator[](Utility->GetFitParamIndx("fracMassBExp"))->operator[](fitParamIndx).c_str()));

  if (GetVar(TotalPDF,"meanR1") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("meanR1"))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,"meanR1",1.0,&myString,&value,&errLo,&errHi);
    }
  if (GetVar(TotalPDF,"sigmaR1") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,"sigmaR1",1.0,&myString,&value,&errLo,&errHi);
    }
  if (GetVar(TotalPDF,"meanR2") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("meanR2"))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,"meanR2",1.0,&myString,&value,&errLo,&errHi);
    }
  if (GetVar(TotalPDF,"sigmaR2") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,"sigmaR2",1.0,&myString,&value,&errLo,&errHi);
    }
  if (GetVar(TotalPDF,"fracMassBRPeak") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("fracMassBRPeak"))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,"fracMassBRPeak",1.0,&myString,&value,&errLo,&errHi);
    }
  if (GetVar(TotalPDF,"meanL1") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("meanL1"))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,"meanL1",1.0,&myString,&value,&errLo,&errHi);
    }
  if (GetVar(TotalPDF,"sigmaL1") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,"sigmaL1",1.0,&myString,&value,&errLo,&errHi);
    }
  if (GetVar(TotalPDF,"meanL2") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("meanL2"))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,"meanL2",1.0,&myString,&value,&errLo,&errHi);
    }
  if (GetVar(TotalPDF,"sigmaL2") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,"sigmaL2",1.0,&myString,&value,&errLo,&errHi);
    }
  if (GetVar(TotalPDF,"fracMassBLPeak") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("fracMassBLPeak"))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,"fracMassBLPeak",1.0,&myString,&value,&errLo,&errHi);
    }
  if (GetVar(TotalPDF,"fracMassBPeak") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("fracMassBPeak"))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,"fracMassBPeak",1.0,&myString,&value,&errLo,&errHi);
    }

  if (GetVar(TotalPDF,"nBkgComb") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("nBkgComb"))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,"nBkgComb",MULTYIELD,&myString,&value,&errLo,&errHi);
    }
  if (GetVar(TotalPDF,"nBkgPeak") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("nBkgPeak"))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,"nBkgPeak",MULTYIELD,&myString,&value,&errLo,&errHi);
    }
  if (GetVar(TotalPDF,"nSig") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("nSig"))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,"nSig",MULTYIELD,&myString,&value,&errLo,&errHi);
    }

  NCoeffPolyBKGpeak1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP1"))->operator[](fitParamIndx).c_str());
  NCoeffPolyBKGcomb1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC1"))->operator[](fitParamIndx).c_str());
  NCoeffPolyBKGpeak2 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP2"))->operator[](fitParamIndx).c_str());
  NCoeffPolyBKGcomb2 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC2"))->operator[](fitParamIndx).c_str());
  if ((NCoeffPolyBKGpeak1 < 0) || (NCoeffPolyBKGcomb1 < 0) || (NCoeffPolyBKGpeak1 > NCOEFFPOLYBKG) || (NCoeffPolyBKGcomb1 > NCOEFFPOLYBKG) ||
      (NCoeffPolyBKGpeak2 < 0) || (NCoeffPolyBKGcomb2 < 0) || (NCoeffPolyBKGpeak2 > NCOEFFPOLYBKG) || (NCoeffPolyBKGcomb2 > NCOEFFPOLYBKG))
    {
      cout << "[ExtractYield::CopyFitResults]\tDegree of poly bkg is not within allowed limits: " << NCoeffPolyBKGpeak1 << "\t" << NCoeffPolyBKGcomb1 << "\t" << NCoeffPolyBKGpeak1 << "\t" << NCoeffPolyBKGcomb2 << endl;
      exit(1);
    }

  for (unsigned int i = 0; i < NCoeffPolyBKGpeak1; i++)
    {
      myCoeff.clear(); myCoeff.str("");
      myCoeff << "p1Poly" << i;
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,myCoeff.str(),1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0)) GetVar(TotalPDF,myCoeff.str().c_str())->setConstant(true);
    }
  for (unsigned int i = 0; i < NCoeffPolyBKGcomb1; i++)
    {
      myCoeff.clear(); myCoeff.str("");
      myCoeff << "c1Poly" << i;
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,myCoeff.str(),1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0)) GetVar(TotalPDF,myCoeff.str().c_str())->setConstant(true);
    }
  for (unsigned int i = 0; i < NCoeffPolyBKGpeak2; i++)
    {
      myCoeff.clear(); myCoeff.str("");
      myCoeff << "p2Poly" << i;
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,myCoeff.str(),1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0)) GetVar(TotalPDF,myCoeff.str().c_str())->setConstant(true);
    }
  for (unsigned int i = 0; i < NCoeffPolyBKGcomb2; i++)
    {
      myCoeff.clear(); myCoeff.str("");
      myCoeff << "c2Poly" << i;
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,myCoeff.str(),1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0)) GetVar(TotalPDF,myCoeff.str().c_str())->setConstant(true);
    }

  if (GetVar(TotalPDF,"FlS") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("FlS"))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,"FlS",1.0,&myString,&value,&errLo,&errHi);
    }
  if (GetVar(TotalPDF,"AfbS") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("AfbS"))->operator[](fitParamIndx).c_str();
      SetValueAndErrors(TotalPDF,"AfbS",1.0,&myString,&value,&errLo,&errHi);
    }
  if (GetVar(TotalPDF,"At2S") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("At2S"))->operator[](fitParamIndx);
      SetValueAndErrors(TotalPDF,"At2S",1.0,&myString,&value,&errLo,&errHi);
    }
  if (GetVar(TotalPDF,"AtimS") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("AtimS"))->operator[](fitParamIndx);
      SetValueAndErrors(TotalPDF,"AtimS",1.0,&myString,&value,&errLo,&errHi);
    }
  if (GetVar(TotalPDF,"FsS") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("FsS"))->operator[](fitParamIndx);
      SetValueAndErrors(TotalPDF,"FsS",1.0,&myString,&value,&errLo,&errHi);
    }
  if (GetVar(TotalPDF,"AsS") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("AsS"))->operator[](fitParamIndx);
      SetValueAndErrors(TotalPDF,"AsS",1.0,&myString,&value,&errLo,&errHi);
    }

  value = 0.0;
  if (GetVar(TotalPDF,"nBkgComb") != NULL) value = value + TotalPDF->getVariables()->getRealValue("nBkgComb");
  if (GetVar(TotalPDF,"nBkgPeak") != NULL) value = value + TotalPDF->getVariables()->getRealValue("nBkgPeak");
  if (GetVar(TotalPDF,"nSig")     != NULL) value = value + TotalPDF->getVariables()->getRealValue("nSig");
  return value;
}


void GenerateParameterFile (RooAbsPdf* TotalPDF, vector<vector<string>*>* fitParam, vector<vector<unsigned int>*>* configParam, RooArgSet* vecConstr, string fileName, unsigned int fileIndx, vector<double>* q2Bins, unsigned int q2BinIndx)
{
  unsigned int NCoeffPolyBKGpeak1;
  unsigned int NCoeffPolyBKGcomb1;
  unsigned int NCoeffPolyBKGpeak2;
  unsigned int NCoeffPolyBKGcomb2;

  stringstream myString;
  vector<string>* vecParStr;

  CopyFitResults(TotalPDF,q2BinIndx,fitParam);

  RooRandom::randomGenerator()->SetSeed(fileIndx*(q2Bins->size()-1) + q2BinIndx + 1);
  cout << "\n@@@ Random seed for parameter file generation set to: " << RooRandom::randomGenerator()->GetSeed() << " @@@" << endl;
  
  if (GetVar(TotalPDF,"FlS") != NULL)
    {
      TotalPDF->getVariables()->setRealValue("FlS",RooRandom::uniform()   * (GetVar(TotalPDF,"FlS")->getMax()   - GetVar(TotalPDF,"FlS")->getMin())   + GetVar(TotalPDF,"FlS")->getMin());
      cout << "Fl generation: lower bound = " << GetVar(TotalPDF,"FlS")->getMin() << "\thigher bound = " << GetVar(TotalPDF,"FlS")->getMax() << endl;
    }

  if (GetVar(TotalPDF,"AfbS") != NULL)
    {
      TotalPDF->getVariables()->setRealValue("AfbS",RooRandom::uniform()  * (GetVar(TotalPDF,"AfbS")->getMax()  - GetVar(TotalPDF,"AfbS")->getMin())  + GetVar(TotalPDF,"AfbS")->getMin());
      cout << "Afb generation: lower bound = " << GetVar(TotalPDF,"AfbS")->getMin() << "\thigher bound = " << GetVar(TotalPDF,"AfbS")->getMax() << endl;
    }

  if (GetVar(TotalPDF,"At2S") != NULL)
    {
      TotalPDF->getVariables()->setRealValue("At2S",RooRandom::uniform()  * (GetVar(TotalPDF,"At2S")->getMax()  - GetVar(TotalPDF,"At2S")->getMin())  + GetVar(TotalPDF,"At2S")->getMin());
      cout << "At2 generation: lower bound = " << GetVar(TotalPDF,"At2S")->getMin() << "\thigher bound = " << GetVar(TotalPDF,"At2S")->getMax() << endl;
    }

  if (GetVar(TotalPDF,"AtimS") != NULL)
    {
      TotalPDF->getVariables()->setRealValue("AtimS",RooRandom::uniform() * (GetVar(TotalPDF,"AtimS")->getMax() - GetVar(TotalPDF,"AtimS")->getMin()) + GetVar(TotalPDF,"AtimS")->getMin());
      cout << "Atim generation: lower bound = " << GetVar(TotalPDF,"AtimS")->getMin() << "\thigher bound = " << GetVar(TotalPDF,"AtimS")->getMax() << endl;
    }

  if (GetVar(TotalPDF,"FsS") != NULL)
    {
      TotalPDF->getVariables()->setRealValue("FsS",RooRandom::uniform()   * (GetVar(TotalPDF,"FsS")->getMax()   - GetVar(TotalPDF,"FsS")->getMin())   + GetVar(TotalPDF,"FsS")->getMin());
      cout << "Fs generation: lower bound = " << GetVar(TotalPDF,"FsS")->getMin() << "\thigher bound = " << GetVar(TotalPDF,"FsS")->getMax() << endl;
    }

  if (GetVar(TotalPDF,"AsS") != NULL)
    {
      TotalPDF->getVariables()->setRealValue("AsS",RooRandom::uniform()   * (GetVar(TotalPDF,"AsS")->getMax()   - GetVar(TotalPDF,"AsS")->getMin())   + GetVar(TotalPDF,"AsS")->getMin());
      cout << "As generation: lower bound = " << GetVar(TotalPDF,"AsS")->getMin() << "\thigher bound = " << GetVar(TotalPDF,"AsS")->getMax() << endl;
    }

  NCoeffPolyBKGpeak1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP1"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGcomb1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC1"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGpeak2 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP2"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGcomb2 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC2"))->operator[](q2BinIndx).c_str());

  // @TMP@
  // if (GetVar(TotalPDF,"nBkgPeak") != NULL)
  //   {
  //     for (unsigned int i = 0; i < NCoeffPolyBKGpeak1; i++)
  // 	{
  // 	  myString.clear(); myString.str("");
  // 	  myString << "p1Poly" << i;
  // 	  if (GetVar(TotalPDF,myString.str().c_str()) != NULL) TotalPDF->getVariables()->setRealValue(myString.str().c_str(),RooRandom::uniform() *
  // 												      (fabs(GetVar(TotalPDF,myString.str().c_str())->getErrorHi()) + fabs(GetVar(TotalPDF,myString.str().c_str())->getErrorLo())) +
  // 												      + GetVar(TotalPDF,myString.str().c_str())->getVal() - fabs(GetVar(TotalPDF,myString.str().c_str())->getErrorLo()));
  // 	}
  //     for (unsigned int i = 0; i < NCoeffPolyBKGpeak2; i++)
  // 	{
  // 	  myString.clear(); myString.str("");
  // 	  myString << "p2Poly" << i;
  // 	  if (GetVar(TotalPDF,myString.str().c_str()) != NULL) TotalPDF->getVariables()->setRealValue(myString.str().c_str(),RooRandom::uniform() *
  // 												      (fabs(GetVar(TotalPDF,myString.str().c_str())->getErrorHi()) + fabs(GetVar(TotalPDF,myString.str().c_str())->getErrorLo())) +
  // 												      + GetVar(TotalPDF,myString.str().c_str())->getVal() - fabs(GetVar(TotalPDF,myString.str().c_str())->getErrorLo()));
  // 	}
  //   }

  for (unsigned int i = 0; i < NCoeffPolyBKGcomb1; i++)
    {
      myString.clear(); myString.str("");
      myString << "c1Poly" << i;
      if (GetVar(TotalPDF,myString.str().c_str()) != NULL) TotalPDF->getVariables()->setRealValue(myString.str().c_str(),RooRandom::uniform() * POLYCOEFRANGE - POLYCOEFRANGE/2.0);
    }
  for (unsigned int i = 0; i < NCoeffPolyBKGcomb2; i++)
    {
      myString.clear(); myString.str("");
      myString << "c2Poly" << i;
      if (GetVar(TotalPDF,myString.str().c_str()) != NULL) TotalPDF->getVariables()->setRealValue(myString.str().c_str(),RooRandom::uniform() * POLYCOEFRANGE - POLYCOEFRANGE/2.0);
    }
  
  vecParStr = SaveFitResults(TotalPDF,q2BinIndx,fitParam,configParam,vecConstr);

  myString.clear(); myString.str("");
  myString << "_" << q2BinIndx << "_" << fileIndx << ".txt";
  fileName.replace(fileName.find(".root"),5,myString.str());
  Utility->SaveFitValues(fileName,vecParStr,q2BinIndx);

  delete vecParStr;
}


void GenerateDataset (RooAbsPdf* TotalPDF, RooArgSet setVar, vector<double>* q2Bins, int specBin, vector<vector<string>*>* fitParam, string fileName)
// ################################################################
// # Use "MULTYIELD" to rescale the entries of the parameter file #
// ################################################################
{
  TFile* NtplFileOut;
  TTree* theTreeOut;
  B0KstMuMuSingleCandTreeContent* NTupleOut;
  RooMCStudy* MyToy;
  RooDataSet* toySample;
  RooRealVar* var;
  unsigned int nEntryToy;


  // ######################
  // # Create output tree #
  // ######################
  NtplFileOut = new TFile(fileName.c_str(), "RECREATE");
  NtplFileOut->mkdir("B0SingleCand");
  NtplFileOut->cd("B0SingleCand");
  theTreeOut = new TTree("B0KstMuMuSingleCandNTuple","B0KstMuMuSingleCandNTuple");
  NTupleOut  = new B0KstMuMuSingleCandTreeContent();
  NTupleOut->Init();
  NTupleOut->ClearNTuple();
  NTupleOut->MakeTreeBranches(theTreeOut);


  cout << "\n@@@ Generating data-like ntuple @@@" << endl;


  // #####################
  // # Toy-MC generation #
  // #####################
  nEntryToy = CopyFitResults(TotalPDF,specBin,fitParam);
  PrintVariables(TotalPDF->getVariables(),"vars");
  MyToy = new RooMCStudy(*TotalPDF,setVar);
  MyToy->generate(1,nEntryToy,true);


  // ########################
  // # Loop over the events #
  // ########################
  toySample = (RooDataSet*)MyToy->genData(0);
  for (unsigned int entry = 0; entry < nEntryToy; entry++)
    {
      if (specBin == Utility->GetJPsiBin(q2Bins))
	{
	  NTupleOut->mumuMass->push_back(Utility->JPsiMass);
	  NTupleOut->mumuMassE->push_back(1e-3);
	}
      else if (specBin == Utility->GetPsiPBin(q2Bins))
	{
	  NTupleOut->mumuMass->push_back(Utility->PsiPrimeMass);
	  NTupleOut->mumuMassE->push_back(1e-3);
	}
      else
	{
	  NTupleOut->mumuMass->push_back(sqrt((q2Bins->operator[](specBin+1)+q2Bins->operator[](specBin)) / 2.));
	  NTupleOut->mumuMassE->push_back(0.0);
	}

      var = NULL;
      var = (RooRealVar*)(toySample->get(entry)->find("B0MassArb"));
      if (var != NULL) NTupleOut->B0MassArb = var->getVal();

      var = NULL;
      var = (RooRealVar*)(toySample->get(entry)->find("CosThetaMuArb"));
      if (var != NULL) NTupleOut->CosThetaMuArb = var->getVal();

      var = NULL;
      var = (RooRealVar*)(toySample->get(entry)->find("CosThetaKArb"));
      if (var != NULL) NTupleOut->CosThetaKArb = var->getVal();

      var = NULL;
      var = (RooRealVar*)(toySample->get(entry)->find("PhiKstMuMuPlaneArb"));
      if (var != NULL) NTupleOut->PhiKstMuMuPlaneArb = var->getVal();
      
      NTupleOut->truthMatchSignal->push_back(false);
      NTupleOut->TrigCat = 1;

      theTreeOut->Fill();
      NTupleOut->ClearNTuple();
    }


  // #############
  // # Save data #
  // #############
  NtplFileOut->cd("B0SingleCand");
  theTreeOut->Write();
  NtplFileOut->Close("R");
  NTupleOut->Destroy();
  delete NTupleOut;
  delete NtplFileOut;
  delete MyToy;
}


void MakeGraphicalScan (RooAbsPdf* TotalPDF, TCanvas* Canv, RooRealVar* x, RooRealVar* y, unsigned int NBins)
// ###################
// # x: cos(theta_l) #
// # y: cos(theta_K) #
// ###################
{
  double Flmin  = 0.0;
  double Flmax  = 1.0;
  double Afbmin = -1.0;
  double Afbmax = 1.0;
  double outputPDF;
  double negNOneg;


  // #############################
  // # Turn off all the printout #
  // #############################
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);


  Canv->cd();
  Canv->SetGrid();
  TH2D* histoVarAngles = new TH2D("histoVarAngles","histoVarAngles",NBins*2,Afbmin,Afbmax,NBins,Flmin,Flmax);
  histoVarAngles->SetXTitle("A_{FB}");
  histoVarAngles->SetYTitle("F_{L}");
  histoVarAngles->SetZTitle("FCN");
  
  for (unsigned int i = 0; i < NBins; i++)
    {
      TotalPDF->getVariables()->setRealValue("FlS",Flmin + (Flmax-Flmin) / ((double)NBins) * (double)i);
      
      for (unsigned int j = 0; j < NBins*2; j++)
	{
	  TotalPDF->getVariables()->setRealValue("AfbS",Afbmin + (Afbmax-Afbmin) / ((double)(NBins*2)) * (double)j);
	  
	  negNOneg = 1.0;
	  if ((x != NULL) && (y != NULL))
	    {
	      for (unsigned int k = 0; k <= NBins; k++)
		{
		  x->setVal(-1.0 + 2.0 / ((double)NBins) * (double)k);
		  for (unsigned int l = 0; l <= NBins; l++)
		    {
		      // Here I remove the corners
		      if (!(((k == 0) && (l == 0))     ||
			    ((k == NBins) && (l == 0)) ||
			    ((k == 0) && (l == NBins)) ||
			    ((k == NBins) && (l == NBins))))
			{
			  y->setVal(-1.0 + 2.0 / ((double)NBins) * (double)l);
			  outputPDF = TotalPDF->getVal();
			  if (outputPDF <= 0.0) negNOneg = 0.0;
			}
		      
		      if ((i*2*NBins*NBins*NBins + j*NBins*NBins + k*NBins + l) % NBins*NBins)
			{
			  cout << "\rIt #" << (i*2*NBins*NBins*NBins + j*NBins*NBins + k*NBins + l);
			  cout << "\tFL = " << TotalPDF->getVariables()->getRealValue("FlS");
			  cout << "\tAFB = " << TotalPDF->getVariables()->getRealValue("AfbS");
			  cout << "\tval = " << outputPDF;
			}
		    }
		}
	    }
	  else if (x != NULL)
	    {
	      for (unsigned int k = 0; k <= NBins; k++)
		{
		  x->setVal(-1.0 + 2.0 / ((double)NBins) * (double)k);
		  outputPDF = TotalPDF->getVal();
		  if (outputPDF <= 0.0) negNOneg = 0.0;
		  
		  if ((i*2*NBins*NBins + j*NBins + k) % NBins)
		    {
		      cout << "\rIt #" << (i*2*NBins*NBins + j*NBins + k);
		      cout << "\tFL = " << TotalPDF->getVariables()->getRealValue("FlS");
		      cout << "\tAFB = " << TotalPDF->getVariables()->getRealValue("AfbS");
		      cout << "\tval = " << outputPDF;
		    }
		}
	    }
	  else if (y != NULL)
	    {
	      for (unsigned int l = 0; l <= NBins; l++)
		{
		  y->setVal(-1.0 + 2.0 / ((double)NBins) * (double)l);
		  outputPDF = TotalPDF->getVal();
		  if (outputPDF <= 0.0) negNOneg = 0.0;

		  if ((i*2*NBins*NBins + j*NBins + l) % NBins)
		    {
		      cout << "\rIt #" << (i*2*NBins*NBins + j*NBins + l);
		      cout << "\tFL = " << TotalPDF->getVariables()->getRealValue("FlS");
		      cout << "\tAFB = " << TotalPDF->getVariables()->getRealValue("AfbS");
		      cout << "\tval = " << outputPDF;
		    }
		}
	    }

	  histoVarAngles->SetBinContent(j+1,i+1,negNOneg);
	}
    }

  histoVarAngles->Draw("gcolz");


  // ########################
  // # Turn on the printout #
  // ########################
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
}


void FitDimuonInvMass (RooDataSet* dataSet, RooAbsPdf** TotalPDFJPsi, RooAbsPdf** TotalPDFPsiP, RooRealVar* x, TCanvas* Canv, bool justPlotMuMuMass, bool justKeepPsi, string plotName)
{
  stringstream myString;
  double sigmaJPsi, sigmaJPsiE;
  double sigmaPsiP, sigmaPsiPE;
  unsigned int NBins;
  int nElements = 4;
  if (justKeepPsi == true) nElements = 3;


  if (justPlotMuMuMass == true)
    {
      myString.clear(); myString.str("");
      myString << "(mumuMass > " << 0.0 << " && mumuMass < " << 4.5 << ")";
      cout << "Cut for dimuon invariant mass: " << myString.str().c_str() << endl;
      RooDataSet* dataSetJPsi = (RooDataSet*)dataSet->reduce(myString.str().c_str());
      x->setRange("subRangeJPsi",0.0,4.5);
 
      Canv->cd();
      RooPlot* myFrameJPsi = x->frame(Range("subRangeJPsi"),Bins(25));
      dataSetJPsi->plotOn(myFrameJPsi,Name(dataSetJPsi->GetName()));
      myFrameJPsi->Draw();

      Canv->Update();
      myString.clear(); myString.str("");
      myString << plotName << ".pdf";
      if (SAVEPLOT == true) Canv->Print(myString.str().c_str());
    }
  else
    {
      Canv->Divide(1,2);


      // ### J/psi ###

      // #####################################################
      // # Define dimuon invariant mass sub-region for J/psi #
      // #####################################################
      myString.clear(); myString.str("");
      myString << "(mumuMass > " << 2.4 << " && mumuMass < " << 3.5 << ")";
      cout << "Cut for dimuon invariant mass: " << myString.str().c_str() << endl;
      RooDataSet* dataSetJPsi = (RooDataSet*)dataSet->reduce(myString.str().c_str());
      x->setRange("subRangeJPsi",2.4,3.5);
      NBins = 50;


      // ##########################################
      // # Define fit variables and pdf for J/psi #
      // ##########################################
      RooRealVar* meanJPsi1  = new RooRealVar("meanJPsi1","J/psi mean Gaussian-1",Utility->JPsiMass,"GeV");
      RooRealVar* sigmaJPsi1 = new RooRealVar("sigmaJPsi1","J/psi sigma Gaussian-1",0.02,"GeV");
      RooAbsPdf* MassJPsi1   = new RooGaussian("MassJPsi1","J/psi Gaussian-1",*x,*meanJPsi1,*sigmaJPsi1);

      RooRealVar* meanJPsi2  = new RooRealVar("meanJPsi2","J/psi mean Gaussian-2",Utility->JPsiMass,"GeV");
      RooRealVar* sigmaJPsi2 = new RooRealVar("sigmaJPsi2","J/psi sigma Gaussian-2",0.09,"GeV");
      RooAbsPdf* MassJPsi2   = new RooGaussian("MassJPsi2","J/psi Gaussian-2",*x,*meanJPsi2,*sigmaJPsi2);

      meanJPsi1->setConstant(false);
      meanJPsi2->setConstant(false);
      sigmaJPsi1->setConstant(false);
      sigmaJPsi2->setConstant(false);

      RooRealVar* fracMassJPsi = new RooRealVar("fracMassJPsi","Fraction of J/psi Gaussians",0.5,0.0,1.0);
      fracMassJPsi->setConstant(false);

      RooAbsPdf* JPsi = new RooAddPdf("JPsi","J/psi pdf",RooArgList(*MassJPsi1,*MassJPsi2),RooArgList(*fracMassJPsi));


      // #####################################################
      // # Define fit variables and pdf for Background J/psi #
      // #####################################################
      RooRealVar* a0JPsi = new RooRealVar("a0JPsi","First coefficient Bkg.poly.",0.0);
      RooPolynomial* BkgmumuMassJPsi = new RooPolynomial("BkgmumuMassJPsi","Background of dimuon mass distribution",*x,RooArgSet(*a0JPsi));
      a0JPsi->setConstant(false);


      // ###########################
      // # Define pdf coefficients #
      // ###########################
      RooRealVar* nJPsiMass = new RooRealVar("nJPsiMass","Number of J/psi events",nJPSIS,"#evts");
      RooRealVar* nBkgMassJPsi = new RooRealVar("nBkgMassJPsi","Number of background events",nJPSIB,"#evts");
      nJPsiMass->setConstant(false);
      nBkgMassJPsi->setConstant(false);

      
      myString.clear(); myString.str("");
      myString << plotName.c_str() << "_JPsi";
      if (justKeepPsi == true) *TotalPDFJPsi = new RooAddPdf(myString.str().c_str(),"Total J/psi mass extended pdf",RooArgList(*JPsi),RooArgList(*nJPsiMass));
      else                     *TotalPDFJPsi = new RooAddPdf(myString.str().c_str(),"Total J/psi mass extended pdf",RooArgList(*JPsi,*BkgmumuMassJPsi),RooArgList(*nJPsiMass,*nBkgMassJPsi));


      // ###################
      // # Make actual fit #
      // ###################
      RooFitResult* JPsiFitResult;
      JPsiFitResult = (*TotalPDFJPsi)->fitTo(*dataSetJPsi,Extended(true),Save(true),SumCoefRange("subRangeJPsi"),Range("subRangeJPsi"));


      // ##################################################
      // # Set p.d.f independent variables to known point #
      // ##################################################
      if (GetVar(*TotalPDFJPsi,x->getPlotLabel()) != NULL) (*TotalPDFJPsi)->getVariables()->setRealValue(x->getPlotLabel(),Utility->B0Mass);
      JPsiFitResult->Print("v");
      ((RooAddPdf*)(*TotalPDFJPsi))->fixCoefRange("subRangeJPsi");
      sigmaJPsi = sqrt((*TotalPDFJPsi)->getVariables()->getRealValue("fracMassJPsi") *
		       (*TotalPDFJPsi)->getVariables()->getRealValue("sigmaJPsi1")*(*TotalPDFJPsi)->getVariables()->getRealValue("sigmaJPsi1") +
		       (1. - (*TotalPDFJPsi)->getVariables()->getRealValue("fracMassJPsi")) *
		       (*TotalPDFJPsi)->getVariables()->getRealValue("sigmaJPsi2")*(*TotalPDFJPsi)->getVariables()->getRealValue("sigmaJPsi2"));
      sigmaJPsiE = 1./(2.*sigmaJPsi) * sqrt( pow((pow((*TotalPDFJPsi)->getVariables()->getRealValue("sigmaJPsi1"),2.)-pow((*TotalPDFJPsi)->getVariables()->getRealValue("sigmaJPsi2"),2.)) * GetVar(*TotalPDFJPsi,"fracMassJPsi")->getError(),2.) +
					     pow(2.*(*TotalPDFJPsi)->getVariables()->getRealValue("fracMassJPsi")  * (*TotalPDFJPsi)->getVariables()->getRealValue("sigmaJPsi1")*GetVar(*TotalPDFJPsi,"sigmaJPsi1")->getError(),2.) +
					     pow(2.*(1.-(*TotalPDFJPsi)->getVariables()->getRealValue("fracMassJPsi")) * (*TotalPDFJPsi)->getVariables()->getRealValue("sigmaJPsi2")*GetVar(*TotalPDFJPsi,"sigmaJPsi2")->getError(),2.) );
      

      // ################
      // # Plot results #
      // ################
      Canv->cd(1);
      RooPlot* myFrameJPsi = x->frame(Range("subRangeJPsi"),Bins(NBins));
      dataSetJPsi->plotOn(myFrameJPsi,Name(dataSetJPsi->GetName()));
      (*TotalPDFJPsi)->plotOn(myFrameJPsi,Name((*TotalPDFJPsi)->getPlotLabel()),LineColor(kBlack));
      (*TotalPDFJPsi)->plotOn(myFrameJPsi,Components(*JPsi),LineStyle(7),LineColor(kBlue));
      if (justKeepPsi != true) (*TotalPDFJPsi)->plotOn(myFrameJPsi,Components(*BkgmumuMassJPsi),LineStyle(4),LineColor(kRed));
      (*TotalPDFJPsi)->paramOn(myFrameJPsi,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(RooArgSet(*nJPsiMass,*nBkgMassJPsi)));

      myString.clear(); myString.str("");
      myString << (*TotalPDFJPsi)->getPlotLabel() << "_paramBox";
      TPaveText* paveTextJPsi = (TPaveText*)myFrameJPsi->findObject(myString.str().c_str());
      paveTextJPsi->AddText(Form("%s%.3f#pm%.3f","#mu_{1} = ",GetVar(*TotalPDFJPsi,"meanJPsi1")->getVal(),GetVar(*TotalPDFJPsi,"meanJPsi1")->getError()));
      paveTextJPsi->AddText(Form("%s%.3f#pm%.3f","#mu_{2} = ",GetVar(*TotalPDFJPsi,"meanJPsi2")->getVal(),GetVar(*TotalPDFJPsi,"meanJPsi2")->getError()));
      paveTextJPsi->AddText(Form("%s%.4f#pm%.4f","< #sigma > = ",sigmaJPsi,sigmaJPsiE));
      paveTextJPsi->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameJPsi->chiSquare((*TotalPDFJPsi)->getPlotLabel(),dataSetJPsi->GetName())));
      paveTextJPsi->AddText(Form("%s%.3f","p-value = ",TMath::Prob(myFrameJPsi->chiSquare((*TotalPDFJPsi)->getPlotLabel(),dataSetJPsi->GetName())*NBins,NBins)));
      CheckGoodFit(JPsiFitResult,paveTextJPsi);
      paveTextJPsi->SetBorderSize(0.0);
      paveTextJPsi->SetFillStyle(0);
      paveTextJPsi->SetTextSize(0.03);
      DrawString(LUMI,myFrameJPsi);
      myFrameJPsi->Draw();

      TLegend* legJPsi = new TLegend(0.75, 0.65, 0.97, 0.89, "");
      for (int i = 0; i < nElements; i++)
      	{
      	  TString objName = myFrameJPsi->nameOf(i);
      	  if (objName == "") continue;
      	  TObject* obj = myFrameJPsi->findObject(objName.Data());
      	  legJPsi->AddEntry(obj,legNames[i],"lp");
      	}
      legJPsi->SetFillStyle(0);
      legJPsi->SetFillColor(0);
      legJPsi->SetTextSize(0.04);
      legJPsi->SetBorderSize(0);
      legJPsi->Draw("same");


      // ### psi(2S) ###

      // #######################################################
      // # Define dimuon invariant mass sub-region for psi(2S) #
      // #######################################################
      myString.clear(); myString.str("");
      myString << "(mumuMass > " << 3.3 << " && mumuMass < " << 4.1 << ")";
      cout << "Cut for dimuon invariant mass: " << myString.str().c_str() << endl;
      RooDataSet* dataSetPsiP = (RooDataSet*)dataSet->reduce(myString.str().c_str());
      x->setRange("subRangePsiP",3.3,4.1);
      NBins = 40;


      // ############################################
      // # Define fit variables and pdf for psi(2S) #
      // ############################################
      RooRealVar* meanPsiP1  = new RooRealVar("meanPsiP1","psi(2S) mean Gaussian-1",Utility->PsiPrimeMass,"GeV");
      RooRealVar* sigmaPsiP1 = new RooRealVar("sigmaPsiP1","psi(2S) sigma Gaussian-1",0.02,"GeV");
      RooAbsPdf* MassPsiP1   = new RooGaussian("MassPsiP1","psi(2S) Gaussian-1",*x,*meanPsiP1,*sigmaPsiP1);

      RooRealVar* meanPsiP2  = new RooRealVar("meanPsiP2","psi(2S) mean Gaussian-2",Utility->PsiPrimeMass,"GeV");
      RooRealVar* sigmaPsiP2 = new RooRealVar("sigmaPsiP2","psi(2S) sigma Gaussian-2",0.09,"GeV");
      RooAbsPdf* MassPsiP2   = new RooGaussian("MassPsiP2","psi(2S) Gaussian-2",*x,*meanPsiP2,*sigmaPsiP2);

      meanPsiP1->setConstant(false);
      meanPsiP2->setConstant(false);
      sigmaPsiP1->setConstant(false);
      sigmaPsiP2->setConstant(false);

      RooRealVar* fracMassPsiP = new RooRealVar("fracMassPsiP","Fraction of psi(2S) Gaussians",0.5,0.0,1.0);
      fracMassPsiP->setConstant(false);

      RooAbsPdf* PsiP = new RooAddPdf("PsiP","psi(2S) pdf",RooArgList(*MassPsiP1,*MassPsiP2),RooArgList(*fracMassPsiP));


      // #######################################################
      // # Define fit variables and pdf for Background psi(2S) #
      // #######################################################
      RooRealVar* a0PsiP = new RooRealVar("a0PsiP","First coefficient Bkg.poly.",0.0);
      RooPolynomial* BkgmumuMassPsiP = new RooPolynomial("BkgmumuMassPsiP","Background of dimuon mass distribution",*x,RooArgSet(*a0PsiP));
      a0PsiP->setConstant(false);


      // ###########################
      // # Define pdf coefficients #
      // ###########################
      RooRealVar* nPsiPMass = new RooRealVar("nPsiPMass","Number of psi(2S) events",nPSIPS,"#evts");
      RooRealVar* nBkgMassPsiP = new RooRealVar("nBkgMassPsiP","Number of background events",nPSIPB,"#evts");
      nPsiPMass->setConstant(false);
      nBkgMassPsiP->setConstant(false);


      myString.clear(); myString.str("");
      myString << plotName.c_str() << "_PsiP";
      if (justKeepPsi == true) *TotalPDFPsiP = new RooAddPdf(myString.str().c_str(),"Total psi(2S) mass extended pdf",RooArgList(*PsiP),RooArgList(*nPsiPMass));
      else                     *TotalPDFPsiP = new RooAddPdf(myString.str().c_str(),"Total psi(2S) mass extended pdf",RooArgList(*PsiP,*BkgmumuMassPsiP),RooArgList(*nPsiPMass,*nBkgMassPsiP));


      // ###################
      // # Make actual fit #
      // ###################
      RooFitResult* PsiPFitResult;
      PsiPFitResult = (*TotalPDFPsiP)->fitTo(*dataSetPsiP,Extended(true),Save(true),SumCoefRange("subRangePsiP"),Range("subRangePsiP"));


      // ##################################################
      // # Set p.d.f independent variables to known point #
      // ##################################################
      if (GetVar(*TotalPDFPsiP,x->getPlotLabel()) != NULL) (*TotalPDFPsiP)->getVariables()->setRealValue(x->getPlotLabel(),Utility->B0Mass);
      PsiPFitResult->Print("v");
      ((RooAddPdf*)(*TotalPDFPsiP))->fixCoefRange("subRangePsiP");
      sigmaPsiP = sqrt((*TotalPDFPsiP)->getVariables()->getRealValue("fracMassPsiP") *
		       (*TotalPDFPsiP)->getVariables()->getRealValue("sigmaPsiP1")*(*TotalPDFPsiP)->getVariables()->getRealValue("sigmaPsiP1") +
		       (1. - (*TotalPDFPsiP)->getVariables()->getRealValue("fracMassPsiP")) *
		       (*TotalPDFPsiP)->getVariables()->getRealValue("sigmaPsiP2")*(*TotalPDFPsiP)->getVariables()->getRealValue("sigmaPsiP2"));
      sigmaPsiPE = 1./(2.*sigmaPsiP) * sqrt( pow((pow((*TotalPDFPsiP)->getVariables()->getRealValue("sigmaPsiP1"),2.)-pow((*TotalPDFPsiP)->getVariables()->getRealValue("sigmaPsiP2"),2.)) * GetVar(*TotalPDFPsiP,"fracMassPsiP")->getError(),2.) +
					     pow(2.*(*TotalPDFPsiP)->getVariables()->getRealValue("fracMassPsiP")  * (*TotalPDFPsiP)->getVariables()->getRealValue("sigmaPsiP1")*GetVar(*TotalPDFPsiP,"sigmaPsiP1")->getError(),2.) +
					     pow(2.*(1.-(*TotalPDFPsiP)->getVariables()->getRealValue("fracMassPsiP")) * (*TotalPDFPsiP)->getVariables()->getRealValue("sigmaPsiP2")*GetVar(*TotalPDFPsiP,"sigmaPsiP2")->getError(),2.) );


      // ################
      // # Plot results #
      // ################
      Canv->cd(2);
      RooPlot* myFramePsiP = x->frame(Range("subRangePsiP"),Bins(NBins));
      dataSetPsiP->plotOn(myFramePsiP,Name(dataSetPsiP->GetName()));
      (*TotalPDFPsiP)->plotOn(myFramePsiP,Name((*TotalPDFPsiP)->getPlotLabel()),LineColor(kBlack));
      (*TotalPDFPsiP)->plotOn(myFramePsiP,Components(*PsiP),LineStyle(7),LineColor(kBlue));
      if (justKeepPsi != true) (*TotalPDFPsiP)->plotOn(myFramePsiP,Components(*BkgmumuMassPsiP),LineStyle(4),LineColor(kRed));
      (*TotalPDFPsiP)->paramOn(myFramePsiP,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(RooArgSet(*nPsiPMass,*nBkgMassPsiP)));

      myString.clear(); myString.str("");
      myString << (*TotalPDFPsiP)->getPlotLabel() << "_paramBox";
      TPaveText* paveTextPsiP = (TPaveText*)myFramePsiP->findObject(myString.str().c_str());
      paveTextPsiP->AddText(Form("%s%.3f#pm%.3f","#mu_{1} = ",GetVar(*TotalPDFPsiP,"meanPsiP1")->getVal(),GetVar(*TotalPDFPsiP,"meanPsiP1")->getError()));
      paveTextPsiP->AddText(Form("%s%.3f#pm%.3f","#mu_{2} = ",GetVar(*TotalPDFPsiP,"meanPsiP2")->getVal(),GetVar(*TotalPDFPsiP,"meanPsiP2")->getError()));
      paveTextPsiP->AddText(Form("%s%.4f#pm%.4f","< #sigma > = ",sigmaPsiP,sigmaPsiPE));
      paveTextPsiP->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFramePsiP->chiSquare((*TotalPDFPsiP)->getPlotLabel(),dataSetPsiP->GetName())));
      paveTextPsiP->AddText(Form("%s%.3f","p-value = ",TMath::Prob(myFramePsiP->chiSquare((*TotalPDFPsiP)->getPlotLabel(),dataSetPsiP->GetName())*NBins,NBins)));
      CheckGoodFit(PsiPFitResult,paveTextPsiP);
      paveTextPsiP->SetBorderSize(0.0);
      paveTextPsiP->SetFillStyle(0);
      paveTextPsiP->SetTextSize(0.03);
      DrawString(LUMI,myFramePsiP);
      myFramePsiP->Draw();

      TLegend* legPsiP = new TLegend(0.75, 0.65, 0.97, 0.89, "");
      for (int i = 0; i < nElements; i++)
	{
	  TString objName = myFramePsiP->nameOf(i);
	  if (objName == "") continue;
	  TObject* obj = myFramePsiP->findObject(objName.Data());
	  legPsiP->AddEntry(obj,legNames[i],"lp");
	}
      legPsiP->SetFillStyle(0);
      legPsiP->SetFillColor(0);
      legPsiP->SetTextSize(0.04);
      legPsiP->SetBorderSize(0);
      legPsiP->Draw("same");


      // ##############
      // # Save plots #
      // ##############
      Canv->Update();
      myString.clear(); myString.str("");
      myString << plotName << ".pdf";
      if (SAVEPLOT == true) Canv->Print(myString.str().c_str());


      // ####################
      // # Save fit results #
      // ####################
      fileFitResults << "\n====================================================================" << endl;
      fileFitResults << "@@@@@@ " << plotName.c_str() << " @@@@@@" << endl;


      fileFitResults << "Fit to J/psi" << endl;
      fileFitResults << "Covariance quality (3=ok): " << JPsiFitResult->covQual() << endl;
      fileFitResults << "Fit status (0=ok): " << JPsiFitResult->status() << endl;
      fileFitResults << "Fit EDM: " << JPsiFitResult->edm() << endl;
      fileFitResults << "Chi2/DoF fit: " << myFrameJPsi->chiSquare((*TotalPDFJPsi)->getPlotLabel(),dataSetJPsi->GetName());
      fileFitResults << "; p-value: " << TMath::Prob(myFrameJPsi->chiSquare((*TotalPDFJPsi)->getPlotLabel(),dataSetJPsi->GetName())*NBINS,NBINS) << endl;
      fileFitResults << "Mean1: " << GetVar(*TotalPDFJPsi,"meanJPsi1")->getVal() << " +/- " << GetVar(*TotalPDFJPsi,"meanJPsi1")->getError();
      fileFitResults << " (" << GetVar(*TotalPDFJPsi,"meanJPsi1")->getErrorHi() << "/" << GetVar(*TotalPDFJPsi,"meanJPsi1")->getErrorLo() << ")" << endl;
      fileFitResults << "Mean2: " << GetVar(*TotalPDFJPsi,"meanJPsi2")->getVal() << " +/- " << GetVar(*TotalPDFJPsi,"meanJPsi2")->getError();
      fileFitResults << " (" << GetVar(*TotalPDFJPsi,"meanJPsi2")->getErrorHi() << "/" << GetVar(*TotalPDFJPsi,"meanJPsi2")->getErrorLo() << ")" << endl;
      fileFitResults << "Sigma-1: " << GetVar(*TotalPDFJPsi,"sigmaJPsi1")->getVal() << " +/- " << GetVar(*TotalPDFJPsi,"sigmaJPsi1")->getError();
      fileFitResults << " (" << GetVar(*TotalPDFJPsi,"sigmaJPsi1")->getErrorHi() << "/" << GetVar(*TotalPDFJPsi,"sigmaJPsi1")->getErrorLo() << ")" << endl;
      fileFitResults << "Sigma-2: " << GetVar(*TotalPDFJPsi,"sigmaJPsi2")->getVal() << " +/- " << GetVar(*TotalPDFJPsi,"sigmaJPsi2")->getError();
      fileFitResults << " (" << GetVar(*TotalPDFJPsi,"sigmaJPsi2")->getErrorHi() << "/" << GetVar(*TotalPDFJPsi,"sigmaJPsi2")->getErrorLo() << ")" << endl;
      fileFitResults << "< Sigma >: " << sigmaJPsi << endl;
      fileFitResults << "Fraction: " << GetVar(*TotalPDFJPsi,"fracMassJPsi")->getVal() << " +/- " << GetVar(*TotalPDFJPsi,"fracMassJPsi")->getError();
      fileFitResults << " (" << GetVar(*TotalPDFJPsi,"fracMassJPsi")->getErrorHi() << "/" << GetVar(*TotalPDFJPsi,"fracMassJPsi")->getErrorLo() << ")" << endl;
      fileFitResults << "Signal yield: " << GetVar(*TotalPDFJPsi,"nJPsiMass")->getVal() << " +/- " << GetVar(*TotalPDFJPsi,"nJPsiMass")->getError();
      fileFitResults << " (" << GetVar(*TotalPDFJPsi,"nJPsiMass")->getErrorHi() << "/" << GetVar(*TotalPDFJPsi,"nJPsiMass")->getErrorLo() << ")" << endl;
      if (justKeepPsi != true)
	{
	  fileFitResults << "Background yield: " << GetVar(*TotalPDFJPsi,"nBkgMassJPsi")->getVal() << " +/- " << GetVar(*TotalPDFJPsi,"nBkgMassJPsi")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDFJPsi,"nBkgMassJPsi")->getErrorHi() << "/" << GetVar(*TotalPDFJPsi,"nBkgMassJPsi")->getErrorLo() << ")" << endl;
	}


      fileFitResults << "\nFit to psi(2S)" << endl;
      fileFitResults << "Covariance quality (3=ok): " << PsiPFitResult->covQual() << endl;
      fileFitResults << "Fit status (0=ok): " << PsiPFitResult->status() << endl;
      fileFitResults << "Fit EDM: " << PsiPFitResult->edm() << endl;
      fileFitResults << "Chi2/DoF fit: " << myFramePsiP->chiSquare((*TotalPDFPsiP)->getPlotLabel(),dataSetPsiP->GetName());
      fileFitResults << "; p-value: " << TMath::Prob(myFramePsiP->chiSquare((*TotalPDFPsiP)->getPlotLabel(),dataSetPsiP->GetName())*NBins,NBins) << endl;
      fileFitResults << "Mean1: " << GetVar(*TotalPDFPsiP,"meanPsiP1")->getVal() << " +/- " << GetVar(*TotalPDFPsiP,"meanPsiP1")->getErrorHi();
      fileFitResults << " (" << GetVar(*TotalPDFPsiP,"meanPsiP1")->getErrorHi() << "/" << GetVar(*TotalPDFPsiP,"meanPsiP1")->getErrorLo() << ")" << endl;
      fileFitResults << "Mean2: " << GetVar(*TotalPDFPsiP,"meanPsiP2")->getVal() << " +/- " << GetVar(*TotalPDFPsiP,"meanPsiP2")->getErrorHi();
      fileFitResults << " (" << GetVar(*TotalPDFPsiP,"meanPsiP2")->getErrorHi() << "/" << GetVar(*TotalPDFPsiP,"meanPsiP2")->getErrorLo() << ")" << endl;
      fileFitResults << "Sigma-1: " << GetVar(*TotalPDFPsiP,"sigmaPsiP1")->getVal() << " +/- " << GetVar(*TotalPDFPsiP,"sigmaPsiP1")->getError();
      fileFitResults << " (" << GetVar(*TotalPDFPsiP,"sigmaPsiP1")->getErrorHi() << "/" << GetVar(*TotalPDFPsiP,"sigmaPsiP1")->getErrorLo() << ")" << endl;
      fileFitResults << "Sigma-2: " << GetVar(*TotalPDFPsiP,"sigmaPsiP2")->getVal() << " +/- " << GetVar(*TotalPDFPsiP,"sigmaPsiP2")->getError();
      fileFitResults << " (" << GetVar(*TotalPDFPsiP,"sigmaPsiP2")->getErrorHi() << "/" << GetVar(*TotalPDFPsiP,"sigmaPsiP2")->getErrorLo() << ")" << endl;
      fileFitResults << "< Sigma >: " << sigmaPsiP << endl;
      fileFitResults << "Fraction: " << GetVar(*TotalPDFPsiP,"fracMassPsiP")->getVal() << " +/- " << GetVar(*TotalPDFPsiP,"fracMassPsiP")->getError();
      fileFitResults << " (" << GetVar(*TotalPDFPsiP,"fracMassPsiP")->getErrorHi() << "/" << GetVar(*TotalPDFPsiP,"fracMassPsiP")->getErrorLo() << ")" << endl;
      fileFitResults << "Signal yield: " << GetVar(*TotalPDFPsiP,"nPsiPMass")->getVal() << " +/- " << GetVar(*TotalPDFPsiP,"nPsiPMass")->getError();
      fileFitResults << " (" << GetVar(*TotalPDFPsiP,"nPsiPMass")->getErrorHi() << "/" << GetVar(*TotalPDFPsiP,"nPsiPMass")->getErrorLo() << ")" << endl;
      if (justKeepPsi != true)	
	{
	  fileFitResults << "Background yield: " << GetVar(*TotalPDFPsiP,"nBkgMassPsiP")->getVal() << " +/- " << GetVar(*TotalPDFPsiP,"nBkgMassPsiP")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDFPsiP,"nBkgMassPsiP")->getErrorHi() << "/" << GetVar(*TotalPDFPsiP,"nBkgMassPsiP")->getErrorLo() << ")" << endl;
	}


      fileFitResults << "\nAmplitude of exclusion zone J/psi: -" << Utility->GetGenericParam("NSigmaPsiBig") << " sigmas --- +" << Utility->GetGenericParam("NSigmaPsiSmall") << " sigmas" << endl;
      fileFitResults << "Amplitude of exclusion zone psi(2S): -" << Utility->GetGenericParam("NSigmaPsiSmall") << " sigmas --- +" << Utility->GetGenericParam("NSigmaPsiSmall") << " sigmas" << endl;
      fileFitResults << "====================================================================" << endl;
    }
}


void MakeDataSets (B0KstMuMuSingleCandTreeContent* NTuple, unsigned int FitType)
{
  stringstream myString;

  // Total likelihood for dimuon invariant mass before applying rejections
  RooAbsPdf* TotalPDFJPsi;
  RooAbsPdf* TotalPDFPsiP;

  // Total likelihood for dimuon invariant mass J/psi region
  RooAbsPdf* TotalPDFJPsi_JPsi;
  RooAbsPdf* TotalPDFPsiP_JPsi;

  // Total likelihood for dimuon invariant mass after applying rejections (Psi regions rejected)
  RooAbsPdf* TotalPDFJPsi_RejPsi;
  RooAbsPdf* TotalPDFPsiP_RejPsi;

  // Total likelihood for dimuon invariant mass after applying rejections (Psi regions kept)
  RooAbsPdf* TotalPDFJPsi_KeepPsi;
  RooAbsPdf* TotalPDFPsiP_KeepPsi;

  RooDataSet* SingleCandNTuple;


  // ###########################
  // # Define useful variables #
  // ###########################
  B0MassArb          = new RooRealVar("B0MassArb","M(#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{+}}}}#kern[-0.3]{#pi}#kern[-0.3]{#lower[0.6]{^{#font[122]{\55}}}}#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#kern[-0.1]{#mu}#kern[-1.3]{#lower[0.6]{^{#font[122]{\55}}}})",Utility->B0Mass - Utility->GetGenericParam("B0MassIntervalLeft"),Utility->B0Mass + Utility->GetGenericParam("B0MassIntervalRight"),"GeV");
  mumuMass           = new RooRealVar("mumuMass","#mu#lower[0.4]{^{#font[122]{+}}}#mu#lower[0.4]{^{#font[122]{\55}}} inv. mass",0.0,6.0,"GeV");
  mumuMassE          = new RooRealVar("mumuMassE","#mu#lower[0.4]{^{#font[122]{+}}}#mu#lower[0.4]{^{#font[122]{\55}}} inv. mass error",0.0,0.5,"GeV");
  CosThetaKArb       = new RooRealVar("CosThetaKArb","cos(#theta#lower[-0.4]{_{#font[122]{K}}})",-1.0,1.0,"");
  CosThetaMuArb      = new RooRealVar("CosThetaMuArb","cos(#theta#lower[-0.4]{_{#font[12]{l}}})",-1.0,1.0,"");
  PhiKstMuMuPlaneArb = new RooRealVar("PhiKstMuMuPlaneArb","Angle (#mu#mu)--(#font[122]{K}#lower[0.4]{^{#font[122]{+}}}#pi#lower[0.4]{^{#font[122]{\55}}}) planes",-Utility->PI,Utility->PI,"rad");
  truthMatchSignal   = new RooRealVar("truthMatchSignal","Truth matching",0.0,1.0,"bool");


  if (!(((FitType >= 21) && (FitType <= 26)) || ((FitType >= 81) && (FitType <= 86)) || ((FitType >= 93) && (FitType <= 96))))
    {
      RooArgSet Vars;
      Vars.add(*B0MassArb);
      Vars.add(*mumuMass);
      Vars.add(*mumuMassE);
      Vars.add(*CosThetaKArb);
      Vars.add(*CosThetaMuArb);
      Vars.add(*PhiKstMuMuPlaneArb);
      Vars.add(*truthMatchSignal);

      SingleCandNTuple = new RooDataSet("SingleCandNTuple","SingleCandNTuple",Vars);


      // #############################
      // # Load values from the tree #
      // #############################
      NTuple->ClearNTuple();
      NTuple->SetBranchAddresses(theTree);
      int nEntries = theTree->GetEntries();
      cout << "@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;
      for (int entry = 0; entry < nEntries; entry++)
	{
	  theTree->GetEntry(entry);

	  if ((!((FitType >= 33) && (FitType <= 36)) && !((FitType >= 53) && (FitType <= 56)) && !((FitType >= 73) && (FitType <= 76)) &&
	       (NTuple->B0MassArb > Utility->B0Mass - Utility->GetGenericParam("B0MassIntervalLeft")) &&
	       (NTuple->B0MassArb < Utility->B0Mass + Utility->GetGenericParam("B0MassIntervalRight"))) ||

	      ((((FitType >= 33) && (FitType <= 36)) || ((FitType >= 53) && (FitType <= 56)) || ((FitType >= 73) && (FitType <= 76))) && 
	       (NTuple->B0pT > Utility->GetSeleCut("B0pT")) && (fabs(NTuple->B0Eta) < Utility->GetSeleCut("B0Eta"))))
	    {
	      Vars.setRealValue("B0MassArb",NTuple->B0MassArb);
	      Vars.setRealValue("mumuMass",NTuple->mumuMass->at(0));
	      Vars.setRealValue("mumuMassE",NTuple->mumuMassE->at(0));
	      Vars.setRealValue("CosThetaKArb",NTuple->CosThetaKArb);
	      Vars.setRealValue("CosThetaMuArb",NTuple->CosThetaMuArb);
	      Vars.setRealValue("PhiKstMuMuPlaneArb",NTuple->PhiKstMuMuPlaneArb);
	      Vars.setRealValue("truthMatchSignal",NTuple->truthMatchSignal->at(0));

	      SingleCandNTuple->add(Vars);
	    }
	}
      SingleCandNTuple->Print("v");


      // ###################################################################
      // # Measure the J/psi and psi(2S) widths before applying rejections #
      // ###################################################################
      if (MakeMuMuPlots == true)
	{
	  TCanvas* cmumuMass_beforeRej = new TCanvas("cmumuMass_beforeRej","cmumuMass_beforeRej",10, 10, 700, 900);
	  FitDimuonInvMass(SingleCandNTuple,&TotalPDFJPsi,&TotalPDFPsiP,mumuMass,cmumuMass_beforeRej,false,false,"TotalPDFPsi_beforeRej");
	}


      // ######################
      // # Define useful cuts #
      // ######################

      myString.clear(); myString.str("");
      // #########################################################################
      // # J/psi / psi(2S) keeping based on the event-by-event dimuon mass error #
      // #########################################################################
      if ((FitType >= 33) && (FitType <= 36))
	myString << "truthMatchSignal == 1 && (mumuMass > (" << Utility->JPsiMass << "-" << Utility->GetGenericParam("NSigmaPsiBig") << "*mumuMassE) && mumuMass < ("
		 << Utility->JPsiMass << "+" << Utility->GetGenericParam("NSigmaPsiSmall") << "*mumuMassE))";
      else if (((FitType >= 53) && (FitType <= 56)) || ((FitType >= 73) && (FitType <= 76))) myString << "truthMatchSignal == 1";
      else myString << "(mumuMass > (" << Utility->JPsiMass << "-" << Utility->GetGenericParam("NSigmaPsiBig") << "*mumuMassE) && mumuMass < ("
		    << Utility->JPsiMass << "+" << Utility->GetGenericParam("NSigmaPsiSmall") << "*mumuMassE))";
      cout << "Cut for B0 --> J/psi K*0: " << myString.str().c_str() << endl;
      SingleCandNTuple_JPsi = (RooDataSet*)SingleCandNTuple->reduce(myString.str().c_str());

      myString.clear(); myString.str("");
      if ((FitType >= 33) && (FitType <= 36))
	myString << "truthMatchSignal == 1 && (mumuMass > (" << Utility->PsiPrimeMass << "-" << Utility->GetGenericParam("NSigmaPsiSmall") << "*mumuMassE) && mumuMass < ("
		 << Utility->PsiPrimeMass << "+" << Utility->GetGenericParam("NSigmaPsiSmall") << "*mumuMassE))";
      else if (((FitType >= 53) && (FitType <= 56)) || ((FitType >= 73) && (FitType <= 76))) myString << "truthMatchSignal == 1";
      else myString << "(mumuMass > (" << Utility->PsiPrimeMass << "-" << Utility->GetGenericParam("NSigmaPsiSmall") << "*mumuMassE) && mumuMass < ("
		    << Utility->PsiPrimeMass << "+" << Utility->GetGenericParam("NSigmaPsiSmall") << "*mumuMassE))";
      cout << "Cut for B0 --> psi(2S) K*0: " << myString.str().c_str() << endl;
      SingleCandNTuple_PsiP = (RooDataSet*)SingleCandNTuple->reduce(myString.str().c_str());


      myString.clear(); myString.str("");
      // #############################################################################
      // # J/psi and psi(2S) rejection based on the event-by-event dimuon mass error #
      // #############################################################################
      myString << "(mumuMass < (" << Utility->JPsiMass << "-" << Utility->GetGenericParam("NSigmaPsiBig") << "*mumuMassE)";
      myString << " || mumuMass > (" << Utility->PsiPrimeMass << "+" << Utility->GetGenericParam("NSigmaPsiSmall") << "*mumuMassE)";
      myString << " || (mumuMass > (" << Utility->JPsiMass << "+" << Utility->GetGenericParam("NSigmaPsiSmall") << "*mumuMassE)";
      myString << " && mumuMass < (" << Utility->PsiPrimeMass << "-" << Utility->GetGenericParam("NSigmaPsiSmall") << "*mumuMassE)))";
      cout << "Cut for B0 --> mu+ mu- K*0 (outside Psi-region): " << myString.str().c_str() << endl;
      if (!((FitType >= 33) && (FitType <= 36)) && !((FitType >= 53) && (FitType <= 56)) && !((FitType >= 73) && (FitType <= 76))) SingleCandNTuple_RejectPsi = (RooDataSet*)SingleCandNTuple->reduce(myString.str().c_str());
      else SingleCandNTuple_RejectPsi = (RooDataSet*)SingleCandNTuple->reduce("mumuMass > 0");


      myString.clear(); myString.str("");
      // ###########################################################################
      // # J/psi and psi(2S) keeping based on the event-by-event dimuon mass error #
      // ###########################################################################
      myString << "((mumuMass > (" << Utility->JPsiMass << "-" << Utility->GetGenericParam("NSigmaPsiBig") << "*mumuMassE)";
      myString << " && mumuMass < (" << Utility->JPsiMass << "+" << Utility->GetGenericParam("NSigmaPsiSmall") << "*mumuMassE)) ||";
      myString << " (mumuMass > (" << Utility->PsiPrimeMass << "-" << Utility->GetGenericParam("NSigmaPsiSmall") << "*mumuMassE)";
      myString << " && mumuMass < (" << Utility->PsiPrimeMass << "+" << Utility->GetGenericParam("NSigmaPsiSmall") << "*mumuMassE)))";
      cout << "Cut for B0 --> mu+ mu- K*0 (inside Psi-region): " << myString.str().c_str() << endl;
      SingleCandNTuple_KeepPsi = (RooDataSet*)SingleCandNTuple->reduce(myString.str().c_str());


      // ##################################################################
      // # Measure the J/psi and psi(2S) widths after applying rejections #
      // ##################################################################
      if (MakeMuMuPlots == true)
	{
	  TCanvas* cmumuMass_JPsi = new TCanvas("cmumuMass_JPsi","cmumuMass_JPsi",10, 10, 700, 900);
	  FitDimuonInvMass(SingleCandNTuple_JPsi,&TotalPDFJPsi_JPsi,&TotalPDFPsiP_JPsi,mumuMass,cmumuMass_JPsi,false,true,"TotalPDFPsi_JPsi");

	  TCanvas* cmumuMass_PsiP = new TCanvas("cmumuMass_PsiP","cmumuMass_PsiP",10, 10, 700, 900);
	  FitDimuonInvMass(SingleCandNTuple_PsiP,&TotalPDFJPsi_JPsi,&TotalPDFPsiP_JPsi,mumuMass,cmumuMass_PsiP,false,true,"TotalPDFPsi_PsiP");
	  
	  TCanvas* cmumuMass_RejPsi = new TCanvas("cmumuMass_RejPsi","cmumuMass_RejPsi",10, 10, 700, 900);
	  FitDimuonInvMass(SingleCandNTuple_RejectPsi,&TotalPDFJPsi_RejPsi,&TotalPDFPsiP_RejPsi,mumuMass,cmumuMass_RejPsi,false,true,"TotalPDFPsi_RejPsi");
	  
	  TCanvas* cmumuMass_KeepPsi = new TCanvas("cmumuMass_KeepPsi","cmumuMass_KeepPsi",10, 10, 700, 900);
	  FitDimuonInvMass(SingleCandNTuple_KeepPsi,&TotalPDFJPsi_KeepPsi,&TotalPDFPsiP_KeepPsi,mumuMass,cmumuMass_KeepPsi,false,true,"TotalPDFPsi_KeepPsi");
	}
    }
}




// ==================
// ===> 1D MODEL <===
// ==================

void InstantiateMassFit (RooAbsPdf** TotalPDF, RooRealVar* x, string fitName, vector<vector<unsigned int>*>* configParam, int parIndx)
{
  // ################################
  // # Read configuration variables #
  // ################################
  unsigned int nGaussPSIintru = configParam->operator[](Utility->GetConfigParamIndx("nGaussPSIintru"))->operator[](parIndx);
  unsigned int fitPSIintru    = configParam->operator[](Utility->GetConfigParamIndx("fitPSIintru"))->operator[](parIndx);
  bool use2GaussS             = configParam->operator[](Utility->GetConfigParamIndx("SigType"))->operator[](parIndx);
  bool use2BExp               = configParam->operator[](Utility->GetConfigParamIndx("BkgType"))->operator[](parIndx);

  stringstream myString;


  // ################################################
  // # Define mass fit variables and pdf for signal #
  // ################################################
  meanS   = new RooRealVar("meanS","Signal mean",Utility->B0Mass,"GeV");

  sigmaS1 = new RooRealVar("sigmaS1","Signal sigma-1",0.0,"GeV");
  MassS1  = new RooGaussian("MassS1","Signal Gaussian-1",*x,*meanS,*sigmaS1);

  sigmaS2 = new RooRealVar("sigmaS2","Signal sigma-2",0.0,"GeV");
  MassS2  = new RooGaussian("MassS2","Signal Gaussian-2",*x,*meanS,*sigmaS2);

  meanS->setConstant(false);
  sigmaS1->setConstant(false);
  sigmaS2->setConstant(false);

  fracMassS = new RooRealVar("fracMassS","Fraction of signal Gaussian",0.0,0.0,1.0);
  fracMassS->setConstant(false);

  if (use2GaussS == true) MassSignal = new RooAddPdf("MassSignal","Signal mass pdf",RooArgList(*MassS1,*MassS2),RooArgList(*fracMassS));
  else                    MassSignal = new RooGaussian(*((RooGaussian*)MassS1),"MassSignal");


  // ##################################################################
  // # Define mass fit variables and pdf for combinatorial background #
  // ##################################################################
  tau1 = new RooRealVar("tau1","First background time constant",0.0,"GeV");
  myString.clear(); myString.str("");
  myString << "exp(-(" << x->getPlotLabel() << "-meanS)/tau1)";
  BkgMassExp1 = new RooGenericPdf("BkgMassExp1",myString.str().c_str(),RooArgSet(*x,*meanS,*tau1));

  tau2 = new RooRealVar("tau2","Second background time constant",0.0,"GeV");
  myString.clear(); myString.str("");
  myString << "exp(-(" << x->getPlotLabel() << "-meanS)/tau2)";
  BkgMassExp2 = new RooGenericPdf("BkgMassExp2",myString.str().c_str(),RooArgSet(*x,*meanS,*tau2));

  tau1->setConstant(false);
  tau2->setConstant(false);

  fracMassBExp = new RooRealVar("fracMassBExp","Fraction of background Exponential",0.0,0.0,1.0);
  fracMassBExp->setConstant(false);

  if (use2BExp == true) BkgMassComb = new RooAddPdf("BkgMassComb","Background mass comb. bkg pdf",RooArgList(*BkgMassExp1,*BkgMassExp2),RooArgList(*fracMassBExp));
  else                  BkgMassComb = new RooGenericPdf(*((RooGenericPdf*)BkgMassExp1),"BkgMassComb");


  // ############################################################
  // # Define mass fit variables and pdf for peaking background #
  // ############################################################
  meanR1         = new RooRealVar("meanR1","Bkg right peak mean-1",0.0,"GeV");
  sigmaR1        = new RooRealVar("sigmaR1","Bkg right peak sigma-1",0.0,"GeV");
  BkgMassRPeak1  = new RooGaussian("BkgMassRPeak1","Bkg right peak-1",*x,*meanR1,*sigmaR1);

  meanR2         = new RooRealVar("meanR2","Bkg right peak mean-2",0.0,"GeV");
  sigmaR2        = new RooRealVar("sigmaR2","Bkg right peak sigma-2",0.0,"GeV");
  BkgMassRPeak2  = new RooGaussian("BkgMassRPeak2","Bkg right peak-2",*x,*meanR2,*sigmaR2);

  fracMassBRPeak = new RooRealVar("fracMassBRPeak","Fraction of background right Peak",0.0,0.0,1.0);
  BkgMassRPeak   = new RooAddPdf("BkgMassRPeak","Right peaking bkg mass pdf",RooArgList(*BkgMassRPeak1,*BkgMassRPeak2),RooArgList(*fracMassBRPeak));

  meanL1         = new RooRealVar("meanL1","Bkg left peak mean-1",0.0,"GeV");
  sigmaL1        = new RooRealVar("sigmaL1","Bkg left peak sigma-1",0.0,"GeV");
  BkgMassLPeak1  = new RooGaussian("BkgMassLPeak1","Bkg left peak-1",*x,*meanL1,*sigmaL1);

  meanL2         = new RooRealVar("meanL2","Bkg left peak mean-2",0.0,"GeV");
  sigmaL2        = new RooRealVar("sigmaL2","Bkg left peak sigma-2",0.0,"GeV");
  BkgMassLPeak2  = new RooGaussian("BkgMassLPeak2","Bkg left peak-2",*x,*meanL2,*sigmaL2);

  fracMassBLPeak = new RooRealVar("fracMassBLPeak","Fraction of background left Peak",0.0,0.0,1.0);
  BkgMassLPeak   = new RooAddPdf("BkgMassLPeak","Left peaking bkg mass pdf",RooArgList(*BkgMassLPeak1,*BkgMassLPeak2),RooArgList(*fracMassBLPeak));

  fracMassBPeak  = new RooRealVar("fracMassBPeak","Fraction of background right-left Peak",0.0,0.0,1.0);

  meanR1->setConstant(false);
  sigmaR1->setConstant(false);
  meanR2->setConstant(false);
  sigmaR2->setConstant(false);
  fracMassBRPeak->setConstant(false);
  
  meanL1->setConstant(false);
  sigmaL1->setConstant(false);
  meanL2->setConstant(false);
  sigmaL2->setConstant(false);
  fracMassBLPeak->setConstant(false);
  
  fracMassBPeak->setConstant(false);
  
  if ((fitPSIintru == 1) || (fitPSIintru == 2))
    if (nGaussPSIintru != 1) BkgMassPeak = new RooAddPdf(*((RooAddPdf*)BkgMassRPeak),"BkgMassPeak");
    else                     BkgMassPeak = new RooGaussian(*((RooGaussian*)BkgMassRPeak1),"BkgMassPeak");
  else
    if (nGaussPSIintru != 1) BkgMassPeak = new RooAddPdf("BkgMassPeak","Peaking bkg mass pdf",RooArgList(*BkgMassRPeak,*BkgMassLPeak),RooArgList(*fracMassBPeak));
    else                     BkgMassPeak = new RooAddPdf("BkgMassPeak","Peaking bkg mass pdf",RooArgList(*BkgMassRPeak1,*BkgMassLPeak1),RooArgList(*fracMassBPeak));


  // ###########################
  // # Define pdf coefficients #
  // ###########################
  nSig = new RooRealVar("nSig","Number of signal events",1.0);
  nBkgComb = new RooRealVar("nBkgComb","Number of comb. background events",1.0);
  nBkgPeak = new RooRealVar("nBkgPeak","Number of peaking background events",1.0);

  nSig->setConstant(false);
  nBkgComb->setConstant(false);
  nBkgPeak->setConstant(false);

  if      (fitPSIintru == 0) *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgList(*MassSignal,*BkgMassComb),RooArgList(*nSig,*nBkgComb));
  else if (fitPSIintru == 1) *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgList(*BkgMassPeak),RooArgList(*nBkgPeak));
  else                       *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgList(*MassSignal,*BkgMassComb,*BkgMassPeak),RooArgList(*nSig,*nBkgComb,*nBkgPeak));
}


RooFitResult* MakeMassFit (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* x, TCanvas* Canv, unsigned int FitType, vector<vector<unsigned int>*>* configParam, int parIndx, RooArgSet* vecConstr, double* NLLvalue, TPaveText* extText, int ID)
{
  // ################################
  // # Read configuration variables #
  // ################################
  unsigned int fitPSIintru = configParam->operator[](Utility->GetConfigParamIndx("fitPSIintru"))->operator[](parIndx);

  RooFitResult* fitResult = NULL;
  stringstream myString;
  double signalSigma  = 0.0;
  double signalSigmaE = 0.0;
  int nElements = 4;
  if (fitPSIintru == 1) nElements = 2;
  else if (GetVar(*TotalPDF,"nBkgPeak") != NULL) nElements = 5;


  // ###################
  // # Make actual fit #
  // ###################
  if (ApplyConstr == true) fitResult = (*TotalPDF)->fitTo(*dataSet,Extended(true),ExternalConstraints(*vecConstr),Save(true),Minos(USEMINOS));
  else                     fitResult = (*TotalPDF)->fitTo(*dataSet,Extended(true),Save(true),Minos(USEMINOS));


  // ##################################################
  // # Set p.d.f independent variables to known point #
  // ##################################################
  if (GetVar(*TotalPDF,x->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(x->getPlotLabel(),Utility->B0Mass);
  if (fitResult != NULL) fitResult->Print("v");
  // fitTo parameters:
  // - NumCPU(n)
  // - Hesse(true) OR Minos(true)
  // - FitOptions(const char*):
  //   - "m" = MIGRAD only, i.e. no MINOS
  //   - "s" = Estimate step size with HESSE before starting MIGRAD
  //   - "h" = Run HESSE after MIGRAD
  //   - "e" = Perform extended MLL fit
  //   - "0" = Run MIGRAD with strategy MINUIT 0 (faster, but no corr. matrix at end)
  //           Does not apply to HESSE or MINOS, if run afterwards
  //   - "q" = Switch off verbose mode
  //   - "l" = Save log file with parameter values at each MINUIT step
  //   - "v" = Show changed parameters at each MINUIT step
  //   - "t" = Time fit
  //   - "r" = Save fit output in RooFitResult object (return value is object RFR pointer)


  if (fitPSIintru != 1)
    {
      if (GetVar(*TotalPDF,"fracMassS") != NULL)
	{
	  signalSigma  = sqrt((*TotalPDF)->getVariables()->getRealValue("fracMassS") * pow((*TotalPDF)->getVariables()->getRealValue("sigmaS1"),2.) +
			      (1.-(*TotalPDF)->getVariables()->getRealValue("fracMassS")) * pow((*TotalPDF)->getVariables()->getRealValue("sigmaS2"),2.));
	  signalSigmaE = 1./(2.*signalSigma) * sqrt( pow((pow((*TotalPDF)->getVariables()->getRealValue("sigmaS1"),2.)-pow((*TotalPDF)->getVariables()->getRealValue("sigmaS2"),2.)) * GetVar(*TotalPDF,"fracMassS")->getError(),2.) +
						     pow(2.*(*TotalPDF)->getVariables()->getRealValue("fracMassS") * (*TotalPDF)->getVariables()->getRealValue("sigmaS1")*GetVar(*TotalPDF,"sigmaS1")->getError(),2.) +
						     pow(2.*(1.-(*TotalPDF)->getVariables()->getRealValue("fracMassS")) * (*TotalPDF)->getVariables()->getRealValue("sigmaS2")*GetVar(*TotalPDF,"sigmaS2")->getError(),2.) );
	}
      else if (GetVar(*TotalPDF,"sigmaS1") != NULL)
	{
	  signalSigma  = (*TotalPDF)->getVariables()->getRealValue("sigmaS1");
	  signalSigmaE = GetVar(*TotalPDF,"sigmaS1")->getError();
	}
    }


  // ################
  // # Plot results #
  // ################
  Canv->cd();
  RooPlot* myFrameX = x->frame(NBINS);
  dataSet->plotOn(myFrameX,Name(MakeName(dataSet,ID).c_str()));
  if (FUNCERRBAND == true) (*TotalPDF)->plotOn(myFrameX,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),VisualizeError(*fitResult,1,true),VLines(),FillColor(kGreen-7));
  else                     (*TotalPDF)->plotOn(myFrameX,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack));
  if (fitPSIintru != 1)
    {
      (*TotalPDF)->plotOn(myFrameX,Components(*MassSignal),LineStyle(7),LineColor(kBlue));
      (*TotalPDF)->plotOn(myFrameX,Components(*BkgMassComb),LineStyle(4),LineColor(kRed));
    }
  if (GetVar(*TotalPDF,"nBkgPeak") != NULL) (*TotalPDF)->plotOn(myFrameX,Components(*BkgMassPeak),LineStyle(3),LineColor(kViolet));

  if      (fitPSIintru == 0) (*TotalPDF)->paramOn(myFrameX,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(RooArgSet(*nSig,*nBkgComb)));
  else if (fitPSIintru == 1) (*TotalPDF)->paramOn(myFrameX,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(RooArgSet(*nBkgPeak)));
  else                       (*TotalPDF)->paramOn(myFrameX,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(RooArgSet(*nSig,*nBkgComb,*nBkgPeak)),ShowConstants(true));
  // Format options:
  // - "N" add name
  // - "E" add error
  // - "A" show asymmetric error
  // - "U" show unit
  // - "H" hide the value


  // #########################################
  // # Compute integrals under signal region #
  // #########################################
  RooAbsReal* wholeSigIntTotalPDF     = NULL;
  RooAbsReal* wholeBkgCombIntTotalPDF = NULL;
  RooAbsReal* wholeBkgPeakIntTotalPDF = NULL;
  RooAbsReal* rangeSigIntTotalPDF     = NULL;
  RooAbsReal* rangeBkgCombIntTotalPDF = NULL;
  RooAbsReal* rangeBkgPeakIntTotalPDF = NULL;
  if (fitPSIintru != 1)
    {
      x->setRange("signal",GetVar(*TotalPDF,"meanS")->getVal() - Utility->GetGenericParam("NSigmaB0")*signalSigma,GetVar(*TotalPDF,"meanS")->getVal() + Utility->GetGenericParam("NSigmaB0")*signalSigma);
      wholeSigIntTotalPDF     = ((RooAbsPdf*)((*TotalPDF)->getComponents()->find("MassSignal")))->createIntegral(*x);
      wholeBkgCombIntTotalPDF = ((RooAbsPdf*)((*TotalPDF)->getComponents()->find("BkgMassComb")))->createIntegral(*x);
      rangeSigIntTotalPDF     = ((RooAbsPdf*)((*TotalPDF)->getComponents()->find("MassSignal")))->createIntegral(*x,Range("signal"));
      rangeBkgCombIntTotalPDF = ((RooAbsPdf*)((*TotalPDF)->getComponents()->find("BkgMassComb")))->createIntegral(*x,Range("signal"));
    }
  else x->setRange("signal",Utility->B0Mass - Utility->GetGenericParam("NSigmaB0")*signalSigma,Utility->B0Mass + Utility->GetGenericParam("NSigmaB0")*signalSigma);
  if (GetVar(*TotalPDF,"nBkgPeak") != NULL)
    {
      wholeBkgPeakIntTotalPDF = ((RooAbsPdf*)((*TotalPDF)->getComponents()->find("BkgMassPeak")))->createIntegral(*x);
      rangeBkgPeakIntTotalPDF = ((RooAbsPdf*)((*TotalPDF)->getComponents()->find("BkgMassPeak")))->createIntegral(*x,Range("signal"));
    }

  myString.clear(); myString.str("");
  myString << (*TotalPDF)->getPlotLabel() << "_paramBox";
  TPaveText* paveTextX = (TPaveText*)myFrameX->findObject(myString.str().c_str());
  if (fitPSIintru != 1)
    {
      paveTextX->AddText(Form("%s%.3f#pm%.3f","#mu = ",GetVar(*TotalPDF,"meanS")->getVal(),GetVar(*TotalPDF,"meanS")->getError()));
      paveTextX->AddText(Form("%s%.4f#pm%.4f","< #sigma > = ",signalSigma,signalSigmaE));
    }
  paveTextX->AddText(Form("%s%.5f","FCN(m#lower[0.4]{^{B0}}#lower[-0.4]{_{PDG}}) = ",(*TotalPDF)->getVal()));
  paveTextX->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameX->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));

  TLegend* legX = new TLegend(0.75, 0.65, 0.97, 0.89, "");
  for (int i = 0; i < nElements; i++)
    {
      TString objName = myFrameX->nameOf(i);
      if (objName == "") continue;
      TObject* obj = myFrameX->findObject(objName.Data());
      legX->AddEntry(obj,legNames[i],"lp");
    }
  legX->SetFillStyle(0);
  legX->SetFillColor(0);
  legX->SetTextSize(0.05);
  legX->SetBorderSize(0);


  // ####################
  // # Save fit results #
  // ####################
  fileFitResults << "====================================================================" << endl;
  fileFitResults << "@@@@@@ Make mass fit @@@@@@" << endl;

  fileFitResults << "Chi2/DoF fit: " << myFrameX->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str());
  fileFitResults << "; p-value: " << TMath::Prob(myFrameX->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())*NBINS,NBINS) << endl;

  *NLLvalue = StoreFitResultsInFile(TotalPDF,fitResult,dataSet,vecConstr);

  if (GetVar(*TotalPDF,"nBkgPeak") != NULL)
    {
      fileFitResults << "Peaking bkg. yield in sig-range (+/- " << Utility->GetGenericParam("NSigmaB0") << "*" << signalSigma << "): "
		     << GetVar(*TotalPDF,"nBkgPeak")->getVal() * rangeBkgPeakIntTotalPDF->getVal() / wholeBkgPeakIntTotalPDF->getVal() << endl;
      if ((GetVar(*TotalPDF,"nSig") != NULL) && (GetVar(*TotalPDF,"nBkgComb") != NULL))
	fileFitResults << "S/sqrt(S+B) in sig-range (+/- " << Utility->GetGenericParam("NSigmaB0") << "*" << signalSigma << "): "
		       << GetVar(*TotalPDF,"nSig")->getVal() * rangeSigIntTotalPDF->getVal() / wholeSigIntTotalPDF->getVal() /
	  sqrt(GetVar(*TotalPDF,"nSig")->getVal()*rangeSigIntTotalPDF->getVal()/wholeSigIntTotalPDF->getVal() + GetVar(*TotalPDF,"nBkgComb")->getVal() * rangeBkgCombIntTotalPDF->getVal() / wholeBkgCombIntTotalPDF->getVal() +
	       GetVar(*TotalPDF,"nBkgPeak")->getVal() * rangeBkgPeakIntTotalPDF->getVal() / wholeBkgPeakIntTotalPDF->getVal()) << endl;
    }
  else if ((GetVar(*TotalPDF,"nSig") != NULL) && (GetVar(*TotalPDF,"nBkgComb") != NULL))
    fileFitResults << "S/sqrt(S+B) in sig-range (+/- " << Utility->GetGenericParam("NSigmaB0") << "*" << signalSigma << "): " << GetVar(*TotalPDF,"nSig")->getVal() * rangeSigIntTotalPDF->getVal() / wholeSigIntTotalPDF->getVal() /
      sqrt(GetVar(*TotalPDF,"nSig")->getVal() * rangeSigIntTotalPDF->getVal() / wholeSigIntTotalPDF->getVal() + GetVar(*TotalPDF,"nBkgComb")->getVal() * rangeBkgCombIntTotalPDF->getVal() / wholeBkgCombIntTotalPDF->getVal()) << endl;

  fileFitResults << "====================================================================" << endl;


  // #######################
  // # Add NLL to the plot #
  // #######################
  paveTextX->AddText(Form("%s%.1f","NLL = ",*NLLvalue));
  CheckGoodFit(fitResult,paveTextX);
  paveTextX->SetBorderSize(0.0);
  paveTextX->SetFillStyle(0);
  paveTextX->SetTextSize(0.03);
  paveTextX->Paint();
  DrawString(LUMI,myFrameX);
  if ((extText != NULL) && (SETBATCH == false))
    {
      if (GetVar(*TotalPDF,"nSig") != NULL) extText->AddText(Form("%s%1.f%s%1.f","Signal yield: ",GetVar(*TotalPDF,"nSig")->getVal()," #pm ",GetVar(*TotalPDF,"nSig")->getError()));
      extText->Paint();
      myFrameX->addObject(extText);
    }
  myFrameX->Draw();
  legX->Draw("same");


  // ##############
  // # Save plots #
  // ##############
  Canv->Update();
  if (SAVEPLOT == true)
    {
      myString.clear(); myString.str("");
      myString << (*TotalPDF)->getPlotLabel() << ".pdf";
      Canv->Print(myString.str().c_str());
    }


  return fitResult;
}


void IterativeMassFitq2Bins (RooDataSet* dataSet,
			     bool useEffPDF,
			     double PsiYield, double PsiYieldErr,
			     RooRealVar* x, RooRealVar* y, RooRealVar* z,
			     int specBin,
			     unsigned int FitType,
			     vector<TH1D*>* VecHistoMeas,
			     vector<double>* q2Bins,
			     vector<vector<unsigned int>*>* configParam, vector<vector<string>*>* fitParam,
			     vector<TF2*>* effFuncs,
			     RooArgSet* vecConstr,
			     unsigned int ID)
{
  double value, errLo, errHi;

  vector<string>* vecParStr;
  stringstream myString;

  double effPsi  = 0.0;
  double effMuMu = 0.0;

  RooArgSet   Vars;
  TCanvas*    cq2Bins[q2Bins->size()-1];
  RooDataSet* dataSet_q2Bins[q2Bins->size()-1];
  RooAbsPdf*  TotalPDFq2Bins[q2Bins->size()-1];
  TPaveText*  extText[q2Bins->size()-1];

  RooFitResult* fitResult;
  RooAbsReal*   EffPDFintegral = NULL;
  RooAbsPdf*    EffPDFnorm     = NULL;
  RooAbsPdf*    EffPDFsign     = NULL;

  double NLLvalue;


  if ((useEffPDF == true) && (fitParam->operator[](0)->size() > 1))
    {
      // ##############################################################
      // # Compute integral of S*E for resonant channel and its error #
      // ##############################################################
 
      // ######################################
      // # Reference to normalization channel #
      // ######################################
      myString.clear(); myString.str("");
      myString << MakeAngWithEffPDF(effFuncs->operator[]((NORMBIN == 3 ? Utility->GetJPsiBin(q2Bins) : Utility->GetPsiPBin(q2Bins))),NULL,y,z,FitType,useEffPDF,&Vars);
      EffPDFnorm = new RooGenericPdf("Efficiency",myString.str().c_str(),Vars);
 
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("FlS"))->operator[]((NORMBIN == 3 ? Utility->GetJPsiBin(q2Bins) : Utility->GetPsiPBin(q2Bins)));
      SetValueAndErrors(EffPDFnorm,"FlS",1.0,&myString,&value,&errLo,&errHi);

      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("AfbS"))->operator[]((NORMBIN == 3 ? Utility->GetJPsiBin(q2Bins) : Utility->GetPsiPBin(q2Bins)));
      SetValueAndErrors(EffPDFnorm,"AfbS",1.0,&myString,&value,&errLo,&errHi);

      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("FsS"))->operator[]((NORMBIN == 3 ? Utility->GetJPsiBin(q2Bins) : Utility->GetPsiPBin(q2Bins)));
      SetValueAndErrors(EffPDFnorm,"FsS",1.0,&myString,&value,&errLo,&errHi);

      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("AsS"))->operator[]((NORMBIN == 3 ? Utility->GetJPsiBin(q2Bins) : Utility->GetPsiPBin(q2Bins)));
      SetValueAndErrors(EffPDFnorm,"AsS",1.0,&myString,&value,&errLo,&errHi);

      PrintVariables(EffPDFnorm->getVariables(),"vars");
      EffPDFintegral = EffPDFnorm->createIntegral(RooArgSet(*y,*z));
      effPsi = EffPDFintegral->getVal();
      
      PrintVariables(EffPDFnorm->getVariables(),"vars");
      cout << "\n@@@ Integral of S*E over (theta_l,theta_K) for normalization channel: " << effPsi << " @@@" << endl;
      
      delete EffPDFnorm;
    }


  fileFitResults << "@@@@@@ Differential branching fraction @@@@@@" << endl;
  for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
    {
      if (!(FitType == 41) && !(FitType == 61) && (Utility->ValIsInPsi(q2Bins,(q2Bins->operator[](i+1)+q2Bins->operator[](i))/2.) == true))
      	{
      	  // ###################################
      	  // # Save zeros into prarameter file #
      	  // ###################################
      	  vecParStr = SaveFitResults(NULL,i,fitParam,configParam,vecConstr);
      	  Utility->SaveFitValues(PARAMETERFILEOUT,vecParStr,i);
      	  delete vecParStr;

      	  continue;
      	}


      // ######################
      // # Make external text #
      // ######################
      extText[i] = new TPaveText(0.63,0.5,0.97,0.63,"NDC");
      extText[i]->AddText(Form("%s%.2f%s%.2f%s","q#lower[0.4]{^{2}}: ",q2Bins->operator[](i)," #font[122]{\55} ",q2Bins->operator[](i+1)," GeV#kern[-0.3]{#lower[0.4]{^{2}}}"));
      extText[i]->SetTextAlign(11);
      extText[i]->SetBorderSize(0.0);
      extText[i]->SetFillStyle(0);
      extText[i]->SetTextSize(0.05);


      fileFitResults << "\nBin[" << i << "]: " << q2Bins->operator[](i) << " <= q^2 < " << q2Bins->operator[](i+1) << endl;

      myString.clear(); myString.str("");
      myString << "c_" << i;
      cq2Bins[i] = new TCanvas(myString.str().c_str(), myString.str().c_str(), 20, 20, 700, 500);

      myString.clear(); myString.str("");
      myString << "(mumuMass*mumuMass) > " << q2Bins->operator[](i) << " && (mumuMass*mumuMass) <= " << q2Bins->operator[](i+1);
      cout << "\nCut string: " << myString.str() << endl;
      dataSet_q2Bins[i] = (RooDataSet*)dataSet->reduce(myString.str().c_str());

      myString.clear(); myString.str("");
      myString << "TotalPDFq2Bin_" << i;
      InstantiateMassFit(&TotalPDFq2Bins[i],x,myString.str(),configParam,i);


      // #####################
      // # Initialize p.d.f. #
      // #####################
      CopyFitResults(TotalPDFq2Bins[i],i,fitParam);


      // #####################
      // # Apply constraints #
      // #####################
      ClearVars(vecConstr);
      BuildMassConstraints(vecConstr,TotalPDFq2Bins[i],"All");


      // ###################################
      // # Print variables and constraints #
      // ###################################
      PrintVariables(TotalPDFq2Bins[i]->getVariables(),"vars");
      PrintVariables(vecConstr,"cons");


      // ###################
      // # Perform the fit #
      // ###################
      fitResult = MakeMassFit(dataSet_q2Bins[i],&TotalPDFq2Bins[i],x,cq2Bins[i],FitType,configParam,i,vecConstr,&NLLvalue,extText[i],ID);


      // ##############################################
      // # Save fit results back into prarameter file #
      // ##############################################
      vecParStr = SaveFitResults(TotalPDFq2Bins[i],i,fitParam,configParam,vecConstr);
      Utility->SaveFitValues(PARAMETERFILEOUT,vecParStr,i);
      delete vecParStr;


      if (useEffPDF == true)
	{
	  if (EffPDFsign != NULL) delete EffPDFsign;
	  // ####################################################
	  // # Compute integral of S*E for signal and its error #
	  // ####################################################

	  Vars.removeAll();
	  myString.clear(); myString.str("");
	  myString << MakeAngWithEffPDF(effFuncs->operator[](i),NULL,y,z,FitType,useEffPDF,&Vars);
	  EffPDFsign = new RooGenericPdf("Efficiency",myString.str().c_str(),Vars);

	  myString.clear(); myString.str("");
	  myString << fitParam->operator[](Utility->GetFitParamIndx("FlS"))->operator[](i);
	  SetValueAndErrors(EffPDFsign,"FlS",1.0,&myString,&value,&errLo,&errHi);

	  myString.clear(); myString.str("");
	  myString << fitParam->operator[](Utility->GetFitParamIndx("AfbS"))->operator[](i);
	  SetValueAndErrors(EffPDFsign,"AfbS",1.0,&myString,&value,&errLo,&errHi);

	  myString.clear(); myString.str("");
	  myString << fitParam->operator[](Utility->GetFitParamIndx("FsS"))->operator[](i);
	  SetValueAndErrors(EffPDFsign,"FsS",1.0,&myString,&value,&errLo,&errHi);

	  myString.clear(); myString.str("");
	  myString << fitParam->operator[](Utility->GetFitParamIndx("AsS"))->operator[](i);
	  SetValueAndErrors(EffPDFsign,"AsS",1.0,&myString,&value,&errLo,&errHi);

	  PrintVariables(EffPDFsign->getVariables(),"vars");
	  EffPDFintegral = EffPDFsign->createIntegral(RooArgSet(*y,*z));
	  effMuMu = EffPDFintegral->getVal();

	  PrintVariables(EffPDFsign->getVariables(),"vars");
	  cout << "\n@@@ Integral of S*E over (theta_l,theta_K) for signal: " << effMuMu << " @@@" << endl;
	}


      // ########################################
      // # Save outcome of the fit in histogram #
      // ########################################
      if (configParam->operator[](Utility->GetConfigParamIndx("fitPSIintru"))->operator[](i) != 1)
	{
	  VecHistoMeas->operator[](0)->SetBinContent(i+1,(GetVar(TotalPDFq2Bins[i],"nSig")->getVal() / effMuMu) / (PsiYield / effPsi) * (NORMBIN == 3 ? Utility->JPsiBF : Utility->PsiPBF) / (q2Bins->operator[](i+1) - q2Bins->operator[](i)) / 1e-7);
	  VecHistoMeas->operator[](0)->SetBinError(i+1,VecHistoMeas->operator[](0)->GetBinContent(i+1) * sqrt(pow(GetVar(TotalPDFq2Bins[i],"nSig")->getError() / GetVar(TotalPDFq2Bins[i],"nSig")->getVal(),2.) + pow(PsiYieldErr / PsiYield,2.)));
	  
	  fileFitResults << "@@@@@@ dBF/dq^2 @@@@@@" << endl;
	  fileFitResults << "dBF/dq^2: " << VecHistoMeas->operator[](0)->GetBinContent(i+1) << " -/+ " << VecHistoMeas->operator[](0)->GetBinError(i+1) << " (";
	  fileFitResults << VecHistoMeas->operator[](0)->GetBinContent(i+1) * sqrt(pow(GetVar(TotalPDFq2Bins[i],"nSig")->getErrorLo() / GetVar(TotalPDFq2Bins[i],"nSig")->getVal(),2.) + pow(PsiYieldErr / PsiYield,2.)) << "/";
	  fileFitResults << VecHistoMeas->operator[](0)->GetBinContent(i+1) * sqrt(pow(GetVar(TotalPDFq2Bins[i],"nSig")->getErrorHi() / GetVar(TotalPDFq2Bins[i],"nSig")->getVal(),2.) + pow(PsiYieldErr / PsiYield,2.)) << ")" << endl;

	  fileFitResults << "To cut and paste in config. file" << endl;
	  fileFitResults << VecHistoMeas->operator[](0)->GetBinContent(i+1) << "   ";
	  fileFitResults << VecHistoMeas->operator[](0)->GetBinContent(i+1) * sqrt(pow(GetVar(TotalPDFq2Bins[i],"nSig")->getErrorLo() / GetVar(TotalPDFq2Bins[i],"nSig")->getVal(),2.) + pow(PsiYieldErr / PsiYield,2.)) << "   ";
	  fileFitResults << VecHistoMeas->operator[](0)->GetBinContent(i+1) * sqrt(pow(GetVar(TotalPDFq2Bins[i],"nSig")->getErrorHi() / GetVar(TotalPDFq2Bins[i],"nSig")->getVal(),2.) + pow(PsiYieldErr / PsiYield,2.)) << endl;
	  fileFitResults << "====================================================================" << endl;

	  myString.clear(); myString.str("");
	  if (CheckGoodFit(fitResult) == true)
	    {
	      myString << ID << "   ";
	      myString << (((useEffPDF == true) && (GetVar(EffPDFsign,"FlS")  != NULL)) ? GetVar(EffPDFsign,"FlS")->getVal()  : -2.0) << "   ";
	      myString << (((useEffPDF == true) && (GetVar(EffPDFsign,"AfbS") != NULL)) ? GetVar(EffPDFsign,"AfbS")->getVal() : -2.0) << "   ";
	      myString << VecHistoMeas->operator[](0)->GetBinContent(i+1) << "   ";
	      myString << NLLvalue;
	    }
	  else myString << ID << "   " << -2.0 << "   " << -2.0 << "   " << -2.0 << "   " << -2.0;

	  fileFitSystematics << myString.str() << endl;
	}
    }
}


void MakeMassToy (RooAbsPdf* TotalPDF, RooRealVar* x, TCanvas* Canv, unsigned int FitType, unsigned int nToy, vector<vector<unsigned int>*>* configParam, int specBin, vector<vector<string>*>* fitParam, RooArgSet* vecConstr, string fileName)
{
  unsigned int nEntryToy;
  stringstream myString;
  unsigned int it = 1;
  RooRealVar* tmpVar;
  RooPlot* myFrame;

  
  // #####################
  // # Initialize p.d.f. #
  // #####################
  nEntryToy = CopyFitResults(TotalPDF,specBin,fitParam);


  // #####################
  // # Apply constraints #
  // #####################
  ClearVars(vecConstr);
  BuildMassConstraints(vecConstr,TotalPDF,"All");


  // ###################################
  // # Print variables and constraints #
  // ###################################
  PrintVariables(TotalPDF->getVariables(),"vars");
  PrintVariables(vecConstr,"cons");


  // #############################
  // # Toy-MC generation and fit #
  // #############################
  RooMCStudy* MyToy = new RooMCStudy(*TotalPDF,*x,Extended(true),ExternalConstraints(*vecConstr),FitOptions(Extended(true),ExternalConstraints(*vecConstr),Minos(USEMINOS))); // Possible options : "Binned()" = faster; Silence()
  MyToy->generateAndFit(nToy,nEntryToy,true);
  Canv->Divide(5,4);


  if ((GetVar(TotalPDF,"meanS") != NULL) && (GetVar(TotalPDF,"meanS")->getError() != 0.0) && (IsInConstraints(vecConstr,"meanS") == false))
    {
      RooRealVar* tmpVar = GetVar(TotalPDF,"meanS");
      RooPlot* myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("meanS"))->operator[](specBin).c_str()) -
  				       0.4*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](specBin).c_str())),
  				       atof(fitParam->operator[](Utility->GetFitParamIndx("meanS"))->operator[](specBin).c_str()) +
  				       0.4*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"sigmaS1") != NULL) && (GetVar(TotalPDF,"sigmaS1")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaS1") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaS1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](specBin).c_str()) -
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](specBin).c_str()) +
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"sigmaS2") != NULL) && (GetVar(TotalPDF,"sigmaS2")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaS2") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaS2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS2"))->operator[](specBin).c_str()) -
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS2"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS2"))->operator[](specBin).c_str()) +
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }
      
  if ((GetVar(TotalPDF,"fracMassS") != NULL) && (GetVar(TotalPDF,"fracMassS")->getError() != 0.0) && (IsInConstraints(vecConstr,"fracMassS") == false))
    {
      tmpVar = GetVar(TotalPDF,"fracMassS");
      myFrame = tmpVar->frame(-0.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"tau1") != NULL) && (GetVar(TotalPDF,"tau1")->getError() != 0.0) && (IsInConstraints(vecConstr,"tau1") == false))
    {
      tmpVar = GetVar(TotalPDF,"tau1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("tau1"))->operator[](specBin).c_str()) - 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("tau1"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("tau1"))->operator[](specBin).c_str()) + 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("tau1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"tau2") != NULL) && (GetVar(TotalPDF,"tau2")->getError() != 0.0) && (IsInConstraints(vecConstr,"tau2") == false))
    {
      tmpVar = GetVar(TotalPDF,"tau2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("tau2"))->operator[](specBin).c_str()) - 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("tau2"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("tau2"))->operator[](specBin).c_str()) + 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("tau2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"fracMassBExp") != NULL) && (GetVar(TotalPDF,"fracMassBExp")->getError() != 0.0) && (IsInConstraints(vecConstr,"fracMassBExp") == false))
    {
      tmpVar = GetVar(TotalPDF,"fracMassBExp");
      myFrame = tmpVar->frame(-0.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"meanR1") != NULL) && (GetVar(TotalPDF,"meanR1")->getError() != 0.0) && (IsInConstraints(vecConstr,"meanR1") == false))
    {
      tmpVar = GetVar(TotalPDF,"meanR1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("meanR1"))->operator[](specBin).c_str()) - 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("meanR1"))->operator[](specBin).c_str()) + 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    } 
   
  if ((GetVar(TotalPDF,"sigmaR1") != NULL) && (GetVar(TotalPDF,"sigmaR1")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaR1") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaR1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](specBin).c_str()) -
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](specBin).c_str()) +
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"meanR2") != NULL) && (GetVar(TotalPDF,"meanR2")->getError() != 0.0) && (IsInConstraints(vecConstr,"meanR2") == false))
    {
      tmpVar = GetVar(TotalPDF,"meanR2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("meanR2"))->operator[](specBin).c_str()) - 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("meanR2"))->operator[](specBin).c_str()) + 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    } 
   
  if ((GetVar(TotalPDF,"sigmaR2") != NULL) && (GetVar(TotalPDF,"sigmaR2")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaR2") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaR2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](specBin).c_str()) -
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](specBin).c_str()) +
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"fracMassBRPeak") != NULL) && (GetVar(TotalPDF,"fracMassBRPeak")->getError() != 0.0) && (IsInConstraints(vecConstr,"fracMassBRPeak") == false))
    {
      tmpVar = GetVar(TotalPDF,"fracMassBRPeak");
      myFrame = tmpVar->frame(-0.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"meanL1") != NULL) && (GetVar(TotalPDF,"meanL1")->getError() != 0.0) && (IsInConstraints(vecConstr,"meanL1") == false))
    {
      tmpVar = GetVar(TotalPDF,"meanL1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("meanL1"))->operator[](specBin).c_str()) - 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("meanL1"))->operator[](specBin).c_str()) + 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"sigmaL1") != NULL) && (GetVar(TotalPDF,"sigmaL1")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaL1") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaL1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](specBin).c_str()) -
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](specBin).c_str()) +
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"meanL2") != NULL) && (GetVar(TotalPDF,"meanL2")->getError() != 0.0) && (IsInConstraints(vecConstr,"meanL2") == false))
    {
      tmpVar = GetVar(TotalPDF,"meanL2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("meanL2"))->operator[](specBin).c_str()) - 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("meanL2"))->operator[](specBin).c_str()) + 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"sigmaL2") != NULL) && (GetVar(TotalPDF,"sigmaL2")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaL2") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaL2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](specBin).c_str()) -
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](specBin).c_str()) +
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"fracMassBLPeak") != NULL) && (GetVar(TotalPDF,"fracMassBLPeak")->getError() != 0.0) && (IsInConstraints(vecConstr,"fracMassBLPeak") == false))
    {
      tmpVar = GetVar(TotalPDF,"fracMassBLPeak");
      myFrame = tmpVar->frame(-0.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"fracMassBPeak") != NULL) && (GetVar(TotalPDF,"fracMassBPeak")->getError() != 0.0) && (IsInConstraints(vecConstr,"fracMassBPeak") == false))
    {
      tmpVar = GetVar(TotalPDF,"fracMassBPeak");
      myFrame = tmpVar->frame(-0.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }
  
  if ((GetVar(TotalPDF,"nBkgPeak") != NULL) && (GetVar(TotalPDF,"nBkgPeak")->getError() != 0.0) && (IsInConstraints(vecConstr,"nBkgPeak") == false))
    {
      tmpVar = GetVar(TotalPDF,"nBkgPeak");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgPeak"))->operator[](specBin).c_str()) -
  			      0.6*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgPeak"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgPeak"))->operator[](specBin).c_str()) +
  			      0.6*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgPeak"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"nSig") != NULL) && (GetVar(TotalPDF,"nSig")->getError() != 0.0) && (IsInConstraints(vecConstr,"nSig") == false))
    {
      tmpVar = GetVar(TotalPDF,"nSig");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("nSig"))->operator[](specBin).c_str()) - 0.6*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("nSig"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("nSig"))->operator[](specBin).c_str()) + 0.6*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("nSig"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"nBkgComb") != NULL) && (GetVar(TotalPDF,"nBkgComb")->getError() != 0.0) && (IsInConstraints(vecConstr,"nBkgComb") == false))
    {
      tmpVar = GetVar(TotalPDF,"nBkgComb");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgComb"))->operator[](specBin).c_str()) -
  			      0.6*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgComb"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgComb"))->operator[](specBin).c_str()) +
  			      0.6*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgComb"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  myFrame = MyToy->plotNLL(-2000.0,0.0);
  Canv->cd(it++);
  myFrame->Draw();
  Canv->Update();


  // #############################################
  // # Re-make all the fits and save the results #
  // #############################################
  TH1D* histoDiff1 = new TH1D("histoDiff1","histoDiff1",60,-30.0,30.0);
  histoDiff1->SetFillColor(kAzure+6);
  histoDiff1->SetXTitle("(fit #font[122]{\55} pdf)");
  histoDiff1->SetYTitle("Entries [#]");

  TH1D* histoPull1 = new TH1D("histoPull1","histoPull1",30,-5.0,5.0);
  histoPull1->SetFillColor(kAzure+6);
  histoPull1->SetXTitle("(fit #font[122]{\55} pdf) / #sigma");
  histoPull1->SetYTitle("Entries [#]");

  TH1D* histoNLL1 = new TH1D("histoNLL1","histoNLL1",30,1.0,-1.0);
  histoNLL1->SetFillColor(kGreen-7);
  histoNLL1->SetXTitle("NLL");
  histoNLL1->SetYTitle("Entries [#]");

  TCanvas* cB0Toy = new TCanvas("cB0Toy","cB0Toy", 20, 20, 700, 500);

  cout << "\n@@@ Now fit total TOY invariant mass @@@" << endl;
  RooDataSet* toySample;
  RooFitResult* fitResult;
  double NLLvalue;

  string varName = "nSig";
  for (unsigned int i = 0; i < nToy; i++)
    {
      cout << "\n@@@ Now fitting toy #" << i << " @@@\n" << endl;

      toySample = (RooDataSet*)MyToy->genData(i);
      CopyFitResults(TotalPDF,specBin,fitParam);
      fitResult = MakeMassFit(toySample,&TotalPDF,x,cB0Toy,FitType,configParam,specBin,vecConstr,&NLLvalue,NULL,i+1);


      // ######################################################
      // # To verify the outcome of the RooFit toy-MC studies #
      // ######################################################
      if ((CheckGoodFit(fitResult) == true) && (GetVar(TotalPDF,varName.c_str()) != NULL))
      	{
      	  if (TotalPDF->getVariables()->getRealValue(varName.c_str()) > atof(fitParam->operator[](Utility->GetFitParamIndx(varName.c_str()))->operator[](specBin).c_str()))
      	    histoPull1->Fill((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(fitParam->operator[](Utility->GetFitParamIndx(varName.c_str()))->operator[](specBin).c_str())) /
			     fabs(GetVar(TotalPDF,varName.c_str())->getErrorLo()));
      	  else histoPull1->Fill((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(fitParam->operator[](Utility->GetFitParamIndx(varName.c_str()))->operator[](specBin).c_str())) /
				fabs(GetVar(TotalPDF,varName.c_str())->getErrorHi()));
	  histoDiff1->Fill(TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(fitParam->operator[](Utility->GetFitParamIndx(varName.c_str()))->operator[](specBin).c_str()));
	  histoNLL1->Fill(NLLvalue);
      	}
      delete fitResult;
    }


  TCanvas* cNLL1 = new TCanvas("cNLL1","cNLL1",10,10,1000,500);
  cNLL1->Divide(3,1);
  cNLL1->cd(1);
  histoDiff1->Draw();
  cNLL1->cd(2);
  histoNLL1->Draw();
  cNLL1->cd(3);
  histoPull1->Draw();
  cNLL1->Update();


  // ##############
  // # Save plots #
  // ##############
  if (SAVEPLOT == true)
    {
      TFile* fNLL;

      Canv->Print(fileName.c_str());

      myString.clear(); myString.str("");
      myString << "_DIFF.root";
      fNLL = new TFile(fileName.replace(fileName.find(".root"),5,myString.str()).c_str(),"RECREATE");
      fNLL->cd();
      histoDiff1->Write();
      fNLL->Close();
      delete fNLL;

      myString.clear(); myString.str("");
      myString << "_PULL.root";
      fNLL = new TFile(fileName.replace(fileName.find("_DIFF.root"),10,myString.str()).c_str(),"RECREATE");
      fNLL->cd();
      histoPull1->Write();
      fNLL->Close();
      delete fNLL;

      myString.clear(); myString.str("");
      myString << "_NLL.root";
      fNLL = new TFile(fileName.replace(fileName.find("_PULL.root"),10,myString.str()).c_str(),"RECREATE");
      fNLL->cd();
      histoNLL1->Write();
      fNLL->Close();
      delete fNLL;
   }
}




// ==================
// ===> 2D MODEL <===
// ==================


void InstantiateMassAngleFit (RooAbsPdf** TotalPDF, bool useEffPDF, RooRealVar* x, RooRealVar* y, string fitName, unsigned int FitType, vector<vector<unsigned int>*>* configParam, vector<vector<string>*>* fitParam, unsigned int parIndx, TF2* effFunc, bool isFirstDataset)
// #########################################
// # x: mass                               #
// # y: angle cos(theta_l) OR cos(theta_K) #
// #########################################
{
  // ################################
  // # Read configuration variables #
  // ################################
  unsigned int nGaussPSIintru = configParam->operator[](Utility->GetConfigParamIndx("nGaussPSIintru"))->operator[](parIndx);
  unsigned int fitPSIintru    = configParam->operator[](Utility->GetConfigParamIndx("fitPSIintru"))->operator[](parIndx);
  bool use2GaussS             = configParam->operator[](Utility->GetConfigParamIndx("SigType"))->operator[](parIndx);
  bool use2BExp               = configParam->operator[](Utility->GetConfigParamIndx("BkgType"))->operator[](parIndx);

  // ###########################
  // # Read polynomial degrees #
  // ###########################
  unsigned int NCoeffPolyBKGpeak1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP1"))->operator[](parIndx).c_str());
  unsigned int NCoeffPolyBKGcomb1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC1"))->operator[](parIndx).c_str());

  stringstream myString;
  RooArgSet VarsC;
  RooArgSet VarsP;


  if (isFirstDataset == true)
    {
      // ################################################
      // # Define mass fit variables and pdf for signal #
      // ################################################
      meanS   = new RooRealVar("meanS","Signal mean",Utility->B0Mass,"GeV");

      sigmaS1 = new RooRealVar("sigmaS1","Signal sigma-1",0.0,"GeV");
      MassS1  = new RooGaussian("MassS1","Signal Gaussian-1",*x,*meanS,*sigmaS1);

      sigmaS2 = new RooRealVar("sigmaS2","Signal sigma-2",0.0,"GeV");
      MassS2  = new RooGaussian("MassS2","Signal Gaussian-2",*x,*meanS,*sigmaS2);

      meanS->setConstant(false);
      sigmaS1->setConstant(false);
      sigmaS2->setConstant(false);

      fracMassS = new RooRealVar("fracMassS","Fraction of signal Gaussian",0.0,0.0,1.0);
      fracMassS->setConstant(false);

      if (use2GaussS == true) MassSignal = new RooAddPdf("MassSignal","Signal mass pdf",RooArgList(*MassS1,*MassS2),RooArgList(*fracMassS));
      else                    MassSignal = new RooGaussian(*((RooGaussian*)MassS1),"MassSignal");


      // #################################################
      // # Define angle fit variables and pdf for signal #
      // #################################################
      RooArgSet* Vars = new RooArgSet("VarsEffPDF");
      myString.clear(); myString.str("");  
      myString << MakeAngWithEffPDF(effFunc,y,NULL,NULL,FitType,useEffPDF,Vars);
      AngleS = new RooGenericPdf("AngleS",myString.str().c_str(),*Vars);

      Signal = new RooProdPdf("Signal","Signal Mass*Angle",RooArgList(*MassSignal,*AngleS));
    }


  // ###################################################################
  // # Define angle fit variables and pdf for combinatorial background #
  // ###################################################################
  for (unsigned int i = 0; i < NCoeffPolyBKGcomb1; i++)
    {
      myString.clear(); myString.str("");
      if ((FitType == 43) || (FitType == 63) || (FitType == 45) || (FitType == 65)) myString << "c1Poly" << i;
      else if ((FitType == 4) || (FitType == 44) || (FitType == 64)) myString << "c2Poly" << i;
      else myString << "c1Poly" << i;
      c1Poly[i] = new RooRealVar(myString.str().c_str(),"Bkg.poly.coef.",0.0);
      c1Poly[i]->setConstant(false);
      VarsC.add(*c1Poly[i]);      
    }
  BkgAngleC1 = new RooPolynomial("BkgAngleC1","Comb.bkg angular shape",*y,VarsC);


  // ##################################################################
  // # Define mass fit variables and pdf for combinatorial background #
  // ##################################################################
  tau1 = new RooRealVar("tau1","First background time constant",0.0,"GeV");
  myString.clear(); myString.str("");
  myString << "exp(-(" << x->getPlotLabel() << "-meanS)/tau1)";
  BkgMassExp1 = new RooGenericPdf("BkgMassExp1",myString.str().c_str(),RooArgSet(*x,*meanS,*tau1));

  tau2 = new RooRealVar("tau2","Second background time constant",0.0,"GeV");
  myString.clear(); myString.str("");
  myString << "exp(-(" << x->getPlotLabel() << "-meanS)/tau2)";
  BkgMassExp2 = new RooGenericPdf("BkgMassExp2",myString.str().c_str(),RooArgSet(*x,*meanS,*tau2));

  tau1->setConstant(false);
  tau2->setConstant(false);

  fracMassBExp = new RooRealVar("fracMassBExp","Fraction of background Exponential",0.0,0.0,1.0);
  fracMassBExp->setConstant(false);

  if (use2BExp == true) BkgMassComb = new RooAddPdf("BkgMassComb","Background mass comb. bkg pdf",RooArgList(*BkgMassExp1,*BkgMassExp2),RooArgList(*fracMassBExp));
  else                  BkgMassComb = new RooGenericPdf(*((RooGenericPdf*)BkgMassExp1),"BkgMassComb");
  BkgMassAngleComb = new RooProdPdf("BkgMassAngleComb","Comb bkg Mass*Angle",RooArgList(*BkgMassComb,*BkgAngleC1));


  // #############################################################
  // # Define angle fit variables and pdf for peaking background #
  // #############################################################
  for (unsigned int i = 0; i < NCoeffPolyBKGpeak1; i++)
    {
      myString.clear(); myString.str("");
      if ((FitType == 43) || (FitType == 63) || (FitType == 45) || (FitType == 65)) myString << "p1Poly" << i;
      else if ((FitType == 4) || (FitType == 44) || (FitType == 64)) myString << "p2Poly" << i;
      else myString << "p1Poly" << i;
      p1Poly[i] = new RooRealVar(myString.str().c_str(),"Bkg.poly.coef.",0.0);
      p1Poly[i]->setConstant(false);
      VarsP.add(*p1Poly[i]);
    }
  BkgAngleP1 = new RooPolynomial("BkgAngleP1","Peak.bkg angular shape",*y,VarsP);


  // ############################################################
  // # Define mass fit variables and pdf for peaking background #
  // ############################################################
  meanR1         = new RooRealVar("meanR1","Bkg right peak mean-1",0.0,"GeV");
  sigmaR1        = new RooRealVar("sigmaR1","Bkg right peak sigma-1",0.0,"GeV");
  BkgMassRPeak1  = new RooGaussian("BkgMassRPeak1","Bkg right peak-1",*x,*meanR1,*sigmaR1);

  meanR2         = new RooRealVar("meanR2","Bkg right peak mean-2",0.0,"GeV");
  sigmaR2        = new RooRealVar("sigmaR2","Bkg right peak sigma-2",0.0,"GeV");
  BkgMassRPeak2  = new RooGaussian("BkgMassRPeak2","Bkg right peak-2",*x,*meanR2,*sigmaR2);

  fracMassBRPeak = new RooRealVar("fracMassBRPeak","Fraction of background right Peak",0.0,0.0,1.0);
  BkgMassRPeak   = new RooAddPdf("BkgMassRPeak","Right peaking bkg mass pdf",RooArgList(*BkgMassRPeak1,*BkgMassRPeak2),RooArgList(*fracMassBRPeak));

  meanL1         = new RooRealVar("meanL1","Bkg left peak mean-1",0.0,"GeV");
  sigmaL1        = new RooRealVar("sigmaL1","Bkg left peak sigma-1",0.0,"GeV");
  BkgMassLPeak1  = new RooGaussian("BkgMassLPeak1","Bkg left peak-1",*x,*meanL1,*sigmaL1);

  meanL2         = new RooRealVar("meanL2","Bkg left peak mean-2",0.0,"GeV");
  sigmaL2        = new RooRealVar("sigmaL2","Bkg left peak sigma-2",0.0,"GeV");
  BkgMassLPeak2  = new RooGaussian("BkgMassLPeak2","Bkg left peak-2",*x,*meanL2,*sigmaL2);

  fracMassBLPeak = new RooRealVar("fracMassBLPeak","Fraction of background left Peak",0.0,0.0,1.0);
  BkgMassLPeak   = new RooAddPdf("BkgMassLPeak","Left peaking bkg mass pdf",RooArgList(*BkgMassLPeak1,*BkgMassLPeak2),RooArgList(*fracMassBLPeak));

  fracMassBPeak  = new RooRealVar("fracMassBPeak","Fraction of background right-left Peak",0.0,0.0,1.0);

  meanR1->setConstant(false);
  sigmaR1->setConstant(false);
  meanR2->setConstant(false);
  sigmaR2->setConstant(false);
  fracMassBRPeak->setConstant(false);

  meanL1->setConstant(false);
  sigmaL1->setConstant(false);
  meanL2->setConstant(false);
  sigmaL2->setConstant(false);
  fracMassBLPeak->setConstant(false);
  
  fracMassBPeak->setConstant(false);
  
  if ((fitPSIintru == 1) || (fitPSIintru == 2))
    if (nGaussPSIintru != 1) BkgMassPeak = new RooAddPdf(*((RooAddPdf*)BkgMassRPeak),"BkgMassPeak");
    else                     BkgMassPeak = new RooGaussian(*((RooGaussian*)BkgMassRPeak1),"BkgMassPeak");
  else
    if (nGaussPSIintru != 1) BkgMassPeak = new RooAddPdf("BkgMassPeak","Peaking bkg mass pdf",RooArgList(*BkgMassRPeak,*BkgMassLPeak),RooArgList(*fracMassBPeak));
    else                     BkgMassPeak = new RooAddPdf("BkgMassPeak","Peaking bkg mass pdf",RooArgList(*BkgMassRPeak1,*BkgMassLPeak1),RooArgList(*fracMassBPeak));
  BkgMassAnglePeak = new RooProdPdf("BkgMassAnglePeak","Peaking bkg Mass*Angle",RooArgList(*BkgMassPeak,*BkgAngleP1));


  // ###########################
  // # Define pdf coefficients #
  // ###########################
  nSig = new RooRealVar("nSig","Number of signal events",1.0);
  nBkgComb = new RooRealVar("nBkgComb","Number of comb. background events",1.0);
  nBkgPeak = new RooRealVar("nBkgPeak","Number of peaking background events",1.0);

  nSig->setConstant(false);
  nBkgComb->setConstant(false);
  nBkgPeak->setConstant(false);

  if (((FitType >= 33) && (FitType <= 35)) || ((FitType >= 53) && (FitType <= 55)) || ((FitType >= 73) && (FitType <= 75)))
                             *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgList(*AngleS),RooArgList(*nSig));
  else if (fitPSIintru == 0) *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgList(*Signal,*BkgMassAngleComb),RooArgList(*nSig,*nBkgComb));
  else if (fitPSIintru == 1) *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgList(*BkgMassAnglePeak),RooArgList(*nBkgPeak));
  else                       *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgList(*Signal,*BkgMassAngleComb,*BkgMassAnglePeak),RooArgList(*nSig,*nBkgComb,*nBkgPeak));
}


RooFitResult* MakeMassAngleFit (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* x, RooRealVar* y, TCanvas* Canv, unsigned int FitType, vector<vector<unsigned int>*>* configParam, unsigned int parIndx, RooArgSet* vecConstr, double* NLLvalue, TPaveText* extText, int ID)
{
  // ################################
  // # Read configuration variables #
  // ################################
  unsigned int fitPSIintru = configParam->operator[](Utility->GetConfigParamIndx("fitPSIintru"))->operator[](parIndx);

  RooFitResult* fitResult = NULL;
  stringstream myString;
  double signalSigma  = 0.0;
  double signalSigmaE = 0.0;
  int nElements = 4;
  if (fitPSIintru == 1) nElements = 2;
  else if (GetVar(*TotalPDF,"nBkgPeak") != NULL) nElements = 5;
  TCanvas* localCanv[3];
  localCanv[0] = new TCanvas("localCanv0","localCanv0",20,20,700,500);
  localCanv[1] = new TCanvas("localCanv1","localCanv1",20,20,700,500);
  localCanv[2] = new TCanvas("localCanv2","localCanv2",20,20,700,500);


  // ###########################
  // # Divide the input canvas #
  // ###########################
  Canv->Divide(3,2);


  if (((FitType >= 33) && (FitType <= 35)) || ((FitType >= 53) && (FitType <= 55)) || ((FitType >= 73) && (FitType <= 75)))
    {
      if      ((MAKEGRAPHSCAN == true) && (FitType == 33)) MakeGraphicalScan(*TotalPDF,Canv,y,NULL);
      else if ((MAKEGRAPHSCAN == true) && (FitType == 34)) MakeGraphicalScan(*TotalPDF,Canv,NULL,y);
      else
	{
	  // ###################
	  // # Make actual fit #
	  // ###################
	  if (ApplyConstr == true) fitResult = (*TotalPDF)->fitTo(*dataSet,Extended(true),ExternalConstraints(*vecConstr),Save(true),Minos(USEMINOS));
	  else                     fitResult = (*TotalPDF)->fitTo(*dataSet,Extended(true),Save(true),Minos(USEMINOS));


	  // ##################################################
	  // # Set p.d.f independent variables to known point #
	  // ##################################################
	  if (GetVar(*TotalPDF,x->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(x->getPlotLabel(),Utility->B0Mass);
	  if (GetVar(*TotalPDF,y->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(y->getPlotLabel(),0.0);
	  if (fitResult != NULL) fitResult->Print("v");


	  // ########################
	  // # Angular plot results #
	  // ########################
	  Canv->cd(2);
	  RooPlot* myFrameY = y->frame(NBINS);
	  dataSet->plotOn(myFrameY,Name(MakeName(dataSet,ID).c_str()));
	  if (FUNCERRBAND == true) (*TotalPDF)->plotOn(myFrameY,Name((*TotalPDF)->getPlotLabel()),VisualizeError(*fitResult,1,true),VLines(),FillColor(kGreen-7));
	  else                     (*TotalPDF)->plotOn(myFrameY,Name((*TotalPDF)->getPlotLabel()));

	  (*TotalPDF)->paramOn(myFrameY,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(RooArgSet(*nSig)));

	  myString.clear(); myString.str("");
	  myString << (*TotalPDF)->getPlotLabel() << "_paramBox";
	  TPaveText* paveTextY = (TPaveText*)myFrameY->findObject(myString.str().c_str());
	  paveTextY->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));

	  TLegend* legY = new TLegend(0.75, 0.65, 0.97, 0.89, "");
	  for (int i = 0; i < 2; i++)
	    {
	      TString objName = myFrameY->nameOf(i);
	      if (objName == "") continue;
	      TObject* obj = myFrameY->findObject(objName.Data());
	      legY->AddEntry(obj,legNames[i],"lp");
	    }
	  legY->SetFillStyle(0);
	  legY->SetFillColor(0);
	  legY->SetTextSize(0.05);
	  legY->SetBorderSize(0);
 

	  // ####################
	  // # Save fit results #
	  // ####################
	  fileFitResults << "====================================================================" << endl;
	  fileFitResults << "@@@@@@ Make angle fit @@@@@@" << endl;

	  fileFitResults << "Chi2/DoF y-fit: " << myFrameY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str());
	  fileFitResults << "; p-value: " << TMath::Prob(myFrameY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())*NBINS,NBINS) << endl;

	  *NLLvalue = StoreFitResultsInFile(TotalPDF,fitResult,dataSet,vecConstr);

	  fileFitResults << "====================================================================" << endl;


	  // #######################
	  // # Add NLL to the plot #
	  // #######################
	  Canv->cd(2);
	  paveTextY->AddText(Form("%s%.1f","NLL = ",*NLLvalue));
	  CheckGoodFit(fitResult,paveTextY);
	  paveTextY->SetTextAlign(11);
	  paveTextY->SetBorderSize(0.0);
	  paveTextY->SetFillStyle(0);
	  paveTextY->SetTextSize(0.04);
	  paveTextY->Paint();
	  myFrameY->Draw();
	  legY->Draw("same");
	}
    }
  else
    {
      // #####################
      // # Make sideband fit #
      // #####################
      if (fitPSIintru != 1)
	{
	  RooAbsPdf* TmpPDF = NULL;
	  RooDataSet* sideBands = NULL;
	  RooRealVar frac("frac","Fraction",0.5,0.0,1.0);
	  RooArgSet constrSidebads;
	  ClearVars(&constrSidebads);
	  BuildAngularConstraints(&constrSidebads,*TotalPDF);


	  // ################
	  // # Save results #
	  // ################
	  fileFitResults << "====================================================================" << endl;
	  fileFitResults << "@@@@@@ B0 mass sideband fit @@@@@@" << endl;
	  fileFitResults << "Amplitude of signal region (+/- n*< Sigma >): " << Utility->GetGenericParam("NSigmaB0S") << " * " << Utility->GetB0Width() << endl;


	  // ##############
	  // # Get p.d.f. #
	  // ##############
	  if (fitPSIintru != 0) TmpPDF = new RooAddPdf("TmpPDF","Temporary p.d.f.",RooArgList(*BkgAngleC1,*BkgAngleP1),RooArgList(frac));
	  else                  TmpPDF = new RooPolynomial(*((RooPolynomial*)BkgAngleC1),"TmpPDF");


	  // #############
	  // # Sidebands #
	  // #############
	  myString.clear(); myString.str("");
	  myString << "B0MassArb < " << (*TotalPDF)->getVariables()->getRealValue("meanS") - Utility->GetGenericParam("NSigmaB0S")*Utility->GetB0Width();
	  myString << " || B0MassArb > " << (*TotalPDF)->getVariables()->getRealValue("meanS") + Utility->GetGenericParam("NSigmaB0S")*Utility->GetB0Width();
	  cout << "Cut for B0 sidebands: " << myString.str().c_str() << endl;
	  sideBands = (RooDataSet*)dataSet->reduce(myString.str().c_str());


	  // ###################
	  // # Make actual fit #
	  // ###################
	  if (ApplyConstr == true) fitResult = TmpPDF->fitTo(*sideBands,ExternalConstraints(constrSidebads),Save(true));
	  else                     fitResult = TmpPDF->fitTo(*sideBands,Save(true));
	  if (fitResult != NULL) fitResult->Print("v");


	  // ####################
	  // # Save fit results #
	  // ####################
	  StorePolyResultsInFile(TotalPDF);


	  delete TmpPDF;
	  delete sideBands;
	  ClearVars(&constrSidebads);
	  if (fitResult != NULL) delete fitResult;
	}


      // ###################
      // # Make actual fit #
      // ###################
      if (ApplyConstr == true) fitResult = (*TotalPDF)->fitTo(*dataSet,Extended(true),ExternalConstraints(*vecConstr),Save(true),Minos(USEMINOS));
      else                     fitResult = (*TotalPDF)->fitTo(*dataSet,Extended(true),Save(true),Minos(USEMINOS));


      // ##################################################
      // # Set p.d.f independent variables to known point #
      // ##################################################
      if (GetVar(*TotalPDF,x->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(x->getPlotLabel(),Utility->B0Mass);
      if (GetVar(*TotalPDF,y->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(y->getPlotLabel(),0.0);
      if (fitResult != NULL) fitResult->Print("v");
      

      if (fitPSIintru != 1)
	{
	  if (GetVar(*TotalPDF,"fracMassS") != NULL)
	    {
	      signalSigma  = sqrt((*TotalPDF)->getVariables()->getRealValue("fracMassS") * pow((*TotalPDF)->getVariables()->getRealValue("sigmaS1"),2.) +
				  (1.-(*TotalPDF)->getVariables()->getRealValue("fracMassS")) * pow((*TotalPDF)->getVariables()->getRealValue("sigmaS2"),2.));
	      signalSigmaE = 1./(2.*signalSigma) * sqrt( pow((pow((*TotalPDF)->getVariables()->getRealValue("sigmaS1"),2.)-pow((*TotalPDF)->getVariables()->getRealValue("sigmaS2"),2.)) * GetVar(*TotalPDF,"fracMassS")->getError(),2.) +
							 pow(2.*(*TotalPDF)->getVariables()->getRealValue("fracMassS") * (*TotalPDF)->getVariables()->getRealValue("sigmaS1")*GetVar(*TotalPDF,"sigmaS1")->getError(),2.) +
							 pow(2.*(1.-(*TotalPDF)->getVariables()->getRealValue("fracMassS")) * (*TotalPDF)->getVariables()->getRealValue("sigmaS2")*GetVar(*TotalPDF,"sigmaS2")->getError(),2.) );
	    }
	  else if (GetVar(*TotalPDF,"sigmaS1") != NULL)
	    {
	      signalSigma  = (*TotalPDF)->getVariables()->getRealValue("sigmaS1");
	      signalSigmaE = GetVar(*TotalPDF,"sigmaS1")->getError();
	    }
	}


      // #####################
      // # Mass plot results #
      // #####################
      Canv->cd(1);
      RooPlot* myFrameX = x->frame(NBINS);
      dataSet->plotOn(myFrameX,Name(MakeName(dataSet,ID).c_str()));
      if (FUNCERRBAND == true) (*TotalPDF)->plotOn(myFrameX,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),VisualizeError(*fitResult,1,true),VLines(),FillColor(kGreen-7));
      else                     (*TotalPDF)->plotOn(myFrameX,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack));
      if (fitPSIintru != 1)
	{
	  (*TotalPDF)->plotOn(myFrameX,Components(*Signal),LineStyle(7),LineColor(kBlue));
	  (*TotalPDF)->plotOn(myFrameX,Components(*BkgMassAngleComb),LineStyle(4),LineColor(kRed));
	}
      if (GetVar(*TotalPDF,"nBkgPeak") != NULL) (*TotalPDF)->plotOn(myFrameX,Components(*BkgMassAnglePeak),LineStyle(3),LineColor(kViolet));
      
      if      (fitPSIintru == 0) (*TotalPDF)->paramOn(myFrameX,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(RooArgSet(*nSig,*nBkgComb)));
      else if (fitPSIintru == 1) (*TotalPDF)->paramOn(myFrameX,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(RooArgSet(*nBkgPeak)));
      else                       (*TotalPDF)->paramOn(myFrameX,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(RooArgSet(*nSig,*nBkgComb,*nBkgPeak)),ShowConstants(true));

      myString.clear(); myString.str("");
      myString << (*TotalPDF)->getPlotLabel() << "_paramBox";
      TPaveText* paveTextX = (TPaveText*)myFrameX->findObject(myString.str().c_str());
      if (fitPSIintru != 1)
	{
	  paveTextX->AddText(Form("%s%.3f#pm%.3f","#mu = ",GetVar(*TotalPDF,"meanS")->getVal(),GetVar(*TotalPDF,"meanS")->getError()));
	  paveTextX->AddText(Form("%s%.4f#pm%.4f","< #sigma > = ",signalSigma,signalSigmaE));
	}
      paveTextX->AddText(Form("%s%.5f","FCN(m#lower[0.4]{^{B0}}#lower[-0.4]{_{PDG}},ang=0) = ",(*TotalPDF)->getVal()));
      paveTextX->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameX->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));

      TLegend* legX = new TLegend(0.75, 0.65, 0.97, 0.89, "");
      for (int i = 0; i < nElements; i++)
      	{
      	  TString objName = myFrameX->nameOf(i);
      	  if (objName == "") continue;
      	  TObject* obj = myFrameX->findObject(objName.Data());
	  legX->AddEntry(obj,legNames[i],"lp");
      	}
      legX->SetFillStyle(0);
      legX->SetFillColor(0);
      legX->SetTextSize(0.05);
      legX->SetBorderSize(0);


      // ########################
      // # Angular plot results #
      // ########################
      Canv->cd(2);
      RooPlot* myFrameY = y->frame(NBINS);
      dataSet->plotOn(myFrameY,Name(MakeName(dataSet,ID).c_str()));
      if (FUNCERRBAND == true) (*TotalPDF)->plotOn(myFrameY,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),VisualizeError(*fitResult,1,true),VLines(),FillColor(kGreen-7));
      else                     (*TotalPDF)->plotOn(myFrameY,Name((*TotalPDF)->getPlotLabel()));
      if (fitPSIintru != 1)
	{
	  (*TotalPDF)->plotOn(myFrameY,Components(*Signal),LineStyle(7),LineColor(kBlue));
	  (*TotalPDF)->plotOn(myFrameY,Components(*BkgMassAngleComb),LineStyle(4),LineColor(kRed));
	}
      if (GetVar(*TotalPDF,"nBkgPeak") != NULL) (*TotalPDF)->plotOn(myFrameY,Components(*BkgMassAnglePeak),LineStyle(3),LineColor(kViolet));
      
      TPaveText* paveTextY = new TPaveText(0.11,0.8,0.4,0.88,"NDC");
      paveTextY->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
      paveTextY->SetTextAlign(11);
      paveTextY->SetBorderSize(0.0);
      paveTextY->SetFillStyle(0);
      paveTextY->SetTextSize(0.04);
      paveTextY->Paint();
      myFrameY->addObject(paveTextY);
      DrawString(LUMI,myFrameY);
      if ((extText != NULL) && (SETBATCH == false))
	{
	  extText->Paint();
	  myFrameY->addObject(extText);
	}
      myFrameY->Draw();

      TLegend* legY = new TLegend(0.75, 0.65, 0.97, 0.89, "");
      for (int i = 0; i < nElements; i++)
      	{
      	  TString objName = myFrameY->nameOf(i);
      	  if (objName == "") continue;
      	  TObject* obj = myFrameY->findObject(objName.Data());
	  legY->AddEntry(obj,legNames[i],"lp");
      	}
      legY->SetFillStyle(0);
      legY->SetFillColor(0);
      legY->SetTextSize(0.05);
      legY->SetBorderSize(0);
      legY->Draw("same");


      // #####################
      // # Fill local canvas #
      // #####################
      localCanv[1]->cd();
      myFrameY->Draw();
      legY->Draw("same");


      // ################################
      // # LogLikelihood parameter scan #
      // ################################
      if (!((FitType >= 21) && (FitType <= 25)) && ((GetVar(*TotalPDF,"FlS") != NULL) || (GetVar(*TotalPDF,"AfbS") != NULL) || (GetVar(*TotalPDF,"At2S") != NULL)))
	{
	  localCanv[2]->cd();
	  RooAbsReal* NLL;
	  if (ApplyConstr == true) NLL = (*TotalPDF)->createNLL(*dataSet,Extended(true),ExternalConstraints(*vecConstr));
	  else                     NLL = (*TotalPDF)->createNLL(*dataSet,Extended(true));
	  RooPlot* myFrameNLL = NULL;
	  if (((FitType == 3) || (FitType == 43) || (FitType == 63)) && (GetVar(*TotalPDF,"FlS") != NULL)) // Fl-Afb-fit
	    {
	      myFrameNLL = GetVar(*TotalPDF,"FlS")->frame();
	      NLL->plotOn(myFrameNLL,ShiftToZero());
	      RooAbsReal* varProfile = NLL->createProfile(*(GetVar(*TotalPDF,"FlS")));
	      varProfile->plotOn(myFrameNLL,LineColor(kRed));
	      myFrameNLL->SetMinimum(0);
	      myFrameNLL->SetMaximum(5);
	    }
	  else if (((FitType == 4) || (FitType == 44) || (FitType == 64)) && (GetVar(*TotalPDF,"AfbS") != NULL)) // Afb-Fl-fit
	    {
	      myFrameNLL = GetVar(*TotalPDF,"AfbS")->frame();
	      NLL->plotOn(myFrameNLL,ShiftToZero());
	      RooAbsReal* varProfile = NLL->createProfile(*GetVar(*TotalPDF,"AfbS"));
	      varProfile->plotOn(myFrameNLL,LineColor(kRed));
	      myFrameNLL->SetMinimum(0);
	      myFrameNLL->SetMaximum(5);
	    }
	  else if (((FitType == 5) || (FitType == 45) || (FitType == 65)) && (GetVar(*TotalPDF,"At2S") != NULL)) // Fl-At2-Atim-fit
	    {
	      myFrameNLL = GetVar(*TotalPDF,"At2S")->frame();
	      NLL->plotOn(myFrameNLL,ShiftToZero());
	      RooAbsReal* varProfile = NLL->createProfile(*GetVar(*TotalPDF,"At2S"));
	      varProfile->plotOn(myFrameNLL,LineColor(kRed));
	      myFrameNLL->SetMinimum(0);
	      myFrameNLL->SetMaximum(5);
	    }
	  myFrameNLL->Draw();

	  TLegend* legNLL = new TLegend(0.15, 0.8, 0.35, 0.89, "");
	  TString legNamesNLL[] = {"NLL","NLL scan"};
	  for (int i = 0; i < 2; i++)
	    {
	      TString objName = myFrameNLL->nameOf(i);
	      if (objName == "") continue;
	      TObject* obj = myFrameNLL->findObject(objName.Data());
	      legNLL->AddEntry(obj,legNamesNLL[i],"lp");
	    }
	  legNLL->SetFillStyle(0);
	  legNLL->SetFillColor(0);
	  legNLL->SetTextSize(0.04);
	  legNLL->SetBorderSize(0);
	  legNLL->Draw("same");

	  delete NLL;
	}


      // ####################
      // # Save fit results #
      // ####################
      fileFitResults << "====================================================================" << endl;
      fileFitResults << "@@@@@@ Make mass - angle fit @@@@@@" << endl;

      fileFitResults << "Chi2/DoF x-fit: " << myFrameX->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str());
      fileFitResults << "; p-value: " << TMath::Prob(myFrameX->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())*NBINS,NBINS) << endl;
      fileFitResults << "Chi2/DoF y-fit: " << myFrameY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str());
      fileFitResults << "; p-value: " << TMath::Prob(myFrameY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())*NBINS,NBINS) << endl;

      *NLLvalue = StoreFitResultsInFile(TotalPDF,fitResult,dataSet,vecConstr);
      StorePolyResultsInFile(TotalPDF);
      fileFitResults << "====================================================================" << endl;


      // #######################
      // # Add NLL to the plot #
      // #######################
      Canv->cd(1);
      paveTextX->AddText(Form("%s%.1f","NLL = ",*NLLvalue));
      CheckGoodFit(fitResult,paveTextX);
      paveTextX->SetBorderSize(0.0);
      paveTextX->SetFillStyle(0);
      paveTextX->SetTextSize(0.03);
      paveTextX->Paint();
      DrawString(LUMI,myFrameX);
      if ((extText != NULL) && (SETBATCH == false))
	{
	  if (GetVar(*TotalPDF,"nSig") != NULL) extText->AddText(Form("%s%1.f%s%1.f","Signal yield: ",GetVar(*TotalPDF,"nSig")->getVal()," #pm ",GetVar(*TotalPDF,"nSig")->getError()));
	  extText->Paint();
	  myFrameX->addObject(extText);
	}
      myFrameX->Draw();
      legX->Draw("same");


      // #####################
      // # Fill local canvas #
      // #####################
      localCanv[0]->cd();
      myFrameX->Draw();
      legX->Draw("same");


      // #####################
      // # Plot on sidebands #
      // #####################
      if ((fitPSIintru != 1) && (ID == 0))
	{
	  // ################
	  // # Save results #
	  // ################
	  fileFitResults << "@@@@@@ B0 mass sideband fit overlaps @@@@@@" << endl;
	  fileFitResults << "Amplitude of signal region: +/- " << Utility->GetGenericParam("NSigmaB0S") << " * < Sigma >" << endl;


	  // ################
	  // # Low sideband #
	  // ################
	  x->setRange("lowSideband",(*TotalPDF)->getVariables()->getRealValue("meanS") - Utility->GetGenericParam("B0MassIntervalLeft"),
		      (*TotalPDF)->getVariables()->getRealValue("meanS") - Utility->GetGenericParam("NSigmaB0S")*signalSigma);


	  // #################
	  // # Signal region #
	  // #################
	  x->setRange("signalRegion",(*TotalPDF)->getVariables()->getRealValue("meanS") - Utility->GetGenericParam("NSigmaB0S")*signalSigma,
		      (*TotalPDF)->getVariables()->getRealValue("meanS") + Utility->GetGenericParam("NSigmaB0S")*signalSigma);


	  // #################
	  // # High sideband #
	  // #################
	  x->setRange("highSideband",(*TotalPDF)->getVariables()->getRealValue("meanS") + Utility->GetGenericParam("NSigmaB0S")*signalSigma,
		      (*TotalPDF)->getVariables()->getRealValue("meanS") + Utility->GetGenericParam("B0MassIntervalRight"));

	
	  // ###########################
	  // # Background plot results #
	  // ###########################
	  Canv->cd(3);
	  RooPlot* myFrameLowSideB = y->frame(NBINS);
	  dataSet->plotOn(myFrameLowSideB,Name(MakeName(dataSet,ID).c_str()),CutRange("lowSideband"));
	  (*TotalPDF)->plotOn(myFrameLowSideB,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),ProjectionRange("lowSideband"));
	  if (fitPSIintru != 1)
	    {
	      (*TotalPDF)->plotOn(myFrameLowSideB,Components(*Signal),LineStyle(7),LineColor(kBlue),ProjectionRange("lowSideband"));
	      (*TotalPDF)->plotOn(myFrameLowSideB,Components(*BkgMassAngleComb),LineStyle(4),LineColor(kRed),ProjectionRange("lowSideband"));
	    }
	  if (GetVar(*TotalPDF,"nBkgPeak") != NULL) (*TotalPDF)->plotOn(myFrameLowSideB,Components(*BkgMassAnglePeak),LineStyle(3),LineColor(kViolet),ProjectionRange("lowSideband"));

	  TPaveText* paveTextLowSideB = new TPaveText(0.11,0.75,0.4,0.88,"NDC");
	  paveTextLowSideB->AddText("Low sideband");
	  paveTextLowSideB->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameLowSideB->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
	  paveTextLowSideB->SetTextAlign(11);
	  paveTextLowSideB->SetBorderSize(0.0);
	  paveTextLowSideB->SetFillStyle(0);
	  paveTextLowSideB->SetTextSize(0.04);
	  paveTextLowSideB->Paint();
	  myFrameLowSideB->addObject(paveTextLowSideB);
	  DrawString(LUMI,myFrameLowSideB);
	  myFrameLowSideB->Draw();

	  TLegend* legLowSideB = new TLegend(0.75, 0.65, 0.97, 0.89, "");
	  for (int i = 0; i < (fitPSIintru != 0 ? 4 : 2); i++)
	    {
	      TString objName = myFrameLowSideB->nameOf(i);
	      if (objName == "") continue;
	      TObject* obj = myFrameLowSideB->findObject(objName.Data());
	      legLowSideB->AddEntry(obj,legNames[i],"lp");
	    }
	  legLowSideB->SetFillStyle(0);
	  legLowSideB->SetFillColor(0);
	  legLowSideB->SetTextSize(0.04);
	  legLowSideB->SetBorderSize(0);
	  legLowSideB->Draw("same");


	  // ################
	  // # Save results #
	  // ################
	  fileFitResults << "Chi2/DoF low sideband: " << myFrameLowSideB->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str());
	  fileFitResults << "; p-value: " << TMath::Prob(myFrameLowSideB->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())*NBINS,NBINS) << endl;


	  // ###########################
	  // # Background plot results #
	  // ###########################
	  Canv->cd(4);
	  RooPlot* myFrameSignalRegion = y->frame(NBINS);
	  dataSet->plotOn(myFrameSignalRegion,Name(MakeName(dataSet,ID).c_str()),CutRange("signalRegion"));
	  (*TotalPDF)->plotOn(myFrameSignalRegion,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),ProjectionRange("signalRegion"));
	  if (fitPSIintru != 1)
	    {
	      (*TotalPDF)->plotOn(myFrameSignalRegion,Components(*Signal),LineStyle(7),LineColor(kBlue),ProjectionRange("signalRegion"));
	      (*TotalPDF)->plotOn(myFrameSignalRegion,Components(*BkgMassAngleComb),LineStyle(4),LineColor(kRed),ProjectionRange("signalRegion"));
	    }
	  if (GetVar(*TotalPDF,"nBkgPeak") != NULL) (*TotalPDF)->plotOn(myFrameSignalRegion,Components(*BkgMassAnglePeak),LineStyle(3),LineColor(kViolet),ProjectionRange("signalRegion"));

	  TPaveText* paveTextSignalRegion = new TPaveText(0.11,0.75,0.4,0.88,"NDC");
	  paveTextSignalRegion->AddText(Form("%s%.1f%s","Signal region: #pm",Utility->GetGenericParam("NSigmaB0S")," < #sigma >"));
	  paveTextSignalRegion->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameSignalRegion->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
	  paveTextSignalRegion->SetTextAlign(11);
	  paveTextSignalRegion->SetBorderSize(0.0);
	  paveTextSignalRegion->SetFillStyle(0);
	  paveTextSignalRegion->SetTextSize(0.04);
	  paveTextSignalRegion->Paint();
	  myFrameSignalRegion->addObject(paveTextSignalRegion);
	  DrawString(LUMI,myFrameSignalRegion);
	  myFrameSignalRegion->Draw();

	  TLegend* legSignalRegion = new TLegend(0.75, 0.65, 0.97, 0.89, "");
	  for (int i = 0; i < (fitPSIintru != 0 ? 4 : 2); i++)
	    {
	      TString objName = myFrameSignalRegion->nameOf(i);
	      if (objName == "") continue;
	      TObject* obj = myFrameSignalRegion->findObject(objName.Data());
	      legSignalRegion->AddEntry(obj,legNames[i],"lp");
	    }
	  legSignalRegion->SetFillStyle(0);
	  legSignalRegion->SetFillColor(0);
	  legSignalRegion->SetTextSize(0.04);
	  legSignalRegion->SetBorderSize(0);
	  legSignalRegion->Draw("same");


	  // ################
	  // # Save results #
	  // ################
	  fileFitResults << "Chi2/DoF signal region: " << myFrameSignalRegion->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str());
	  fileFitResults << "; p-value: " << TMath::Prob(myFrameSignalRegion->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())*NBINS,NBINS) << endl;


	  // ###########################
	  // # Background plot results #
	  // ###########################
	  Canv->cd(5);
	  RooPlot* myFrameHighSideB = y->frame(NBINS);
	  dataSet->plotOn(myFrameHighSideB,Name(MakeName(dataSet,ID).c_str()),CutRange("highSideband"));
	  (*TotalPDF)->plotOn(myFrameHighSideB,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),ProjectionRange("highSideband"));
	  if (fitPSIintru != 1)
	    {
	      (*TotalPDF)->plotOn(myFrameHighSideB,Components(*Signal),LineStyle(7),LineColor(kBlue),ProjectionRange("highSideband"));
	      (*TotalPDF)->plotOn(myFrameHighSideB,Components(*BkgMassAngleComb),LineStyle(4),LineColor(kRed),ProjectionRange("highSideband"));
	    }
	  if (GetVar(*TotalPDF,"nBkgPeak") != NULL) (*TotalPDF)->plotOn(myFrameHighSideB,Components(*BkgMassAnglePeak),LineStyle(3),LineColor(kViolet),ProjectionRange("highSideband"));

	  TPaveText* paveTextHighSideB = new TPaveText(0.11,0.75,0.4,0.88,"NDC");
	  paveTextHighSideB->AddText("High sideband");
	  paveTextHighSideB->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameHighSideB->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
	  paveTextHighSideB->SetTextAlign(11);
	  paveTextHighSideB->SetBorderSize(0.0);
	  paveTextHighSideB->SetFillStyle(0);
	  paveTextHighSideB->SetTextSize(0.04);
	  paveTextHighSideB->Paint();
	  myFrameHighSideB->addObject(paveTextHighSideB);
	  DrawString(LUMI,myFrameHighSideB);
	  myFrameHighSideB->Draw();

	  TLegend* legHighSideB = new TLegend(0.75, 0.65, 0.97, 0.89, "");
	  for (int i = 0; i < (fitPSIintru != 0 ? 4 : 2); i++)
	    {
	      TString objName = myFrameHighSideB->nameOf(i);
	      if (objName == "") continue;
	      TObject* obj = myFrameHighSideB->findObject(objName.Data());
	      legHighSideB->AddEntry(obj,legNames[i],"lp");
	    }
	  legHighSideB->SetFillStyle(0);
	  legHighSideB->SetFillColor(0);
	  legHighSideB->SetTextSize(0.04);
	  legHighSideB->SetBorderSize(0);
	  legHighSideB->Draw("same");


	  // ################
	  // # Save results #
	  // ################
	  fileFitResults << "Chi2/DoF high sideband: " << myFrameHighSideB->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str());
	  fileFitResults << "; p-value: " << TMath::Prob(myFrameHighSideB->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())*NBINS,NBINS) << endl;


	  fileFitResults << "====================================================================" << endl;
	}
    }


  // ##############
  // # Save plots #
  // ##############
  Canv->Update();
  if (SAVEPLOT == true)
    {
      myString.clear(); myString.str("");
      myString << (*TotalPDF)->getPlotLabel() << ".pdf";
      Canv->Print(myString.str().c_str());
    }

  if (SETBATCH == true)
    {
      delete localCanv[0];
      delete localCanv[1];
      delete localCanv[2];
    }
  else
    {
      localCanv[0]->Update();
      localCanv[1]->Update();
      localCanv[2]->Update();
    }


  return fitResult;
}


void IterativeMassAngleFitq2Bins (RooDataSet* dataSet,
				  bool useEffPDF,
				  double PsiYield, double PsiYieldErr,
				  RooRealVar* x, RooRealVar* y, RooRealVar* z,
				  int specBin,
				  unsigned int FitType,
				  vector<TH1D*>* VecHistoMeas,
				  vector<double>* q2Bins,
				  vector<vector<unsigned int>*>* configParam, vector<vector<string>*>* fitParam,
				  vector<TF2*>* effFuncs,
				  RooArgSet* vecConstr,
				  unsigned int ID)
{
  double value, errLo, errHi;

  vector<string>* vecParStr;
  stringstream myString;

  double effPsi  = 0.0;
  double effMuMu = 0.0;

  RooArgSet   Vars;
  TCanvas*    cq2Bins[q2Bins->size()-1];
  RooDataSet* dataSet_q2Bins[q2Bins->size()-1];
  RooAbsPdf*  TotalPDFq2Bins[q2Bins->size()-1];
  TPaveText*  extText[q2Bins->size()-1];

  RooFitResult* fitResult;
  RooAbsReal*   EffPDFintegral = NULL;
  RooAbsPdf*    EffPDFnorm     = NULL;
  RooAbsPdf*    EffPDFsign     = NULL;

  double NLLvalue;
  
  RooRealVar* fitVar;


  // #######################
  // # Assign fit variable #
  // #######################
  if      ((FitType == 3) || (FitType == 33) || (FitType == 43) || (FitType == 53) || (FitType == 63) || (FitType == 73)) fitVar = z;
  else if ((FitType == 4) || (FitType == 34) || (FitType == 44) || (FitType == 54) || (FitType == 64) || (FitType == 74)) fitVar = y;
  else
    {
      cout << "@TMP@: not implemented yet for fit type " << FitType << " (Fl-At2-Atim-fit)" << endl;
      exit(1);
    }


  if ((useEffPDF == true) && (fitParam->operator[](0)->size() > 1))
    {
      // ##############################################################
      // # Compute integral of S*E for resonant channel and its error #
      // ##############################################################
 
      // ######################################
      // # Reference to normalization channel #
      // ######################################
      myString.clear(); myString.str("");
      myString << MakeAngWithEffPDF(effFuncs->operator[]((NORMBIN == 3 ? Utility->GetJPsiBin(q2Bins) : Utility->GetPsiPBin(q2Bins))),NULL,y,z,1,useEffPDF,&Vars);
      EffPDFnorm = new RooGenericPdf("Efficiency",myString.str().c_str(),Vars);
 
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("FlS"))->operator[]((NORMBIN == 3 ? Utility->GetJPsiBin(q2Bins) : Utility->GetPsiPBin(q2Bins)));
      SetValueAndErrors(EffPDFnorm,"FlS",1.0,&myString,&value,&errLo,&errHi);

      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("AfbS"))->operator[]((NORMBIN == 3 ? Utility->GetJPsiBin(q2Bins) : Utility->GetPsiPBin(q2Bins)));
      SetValueAndErrors(EffPDFnorm,"AfbS",1.0,&myString,&value,&errLo,&errHi);

      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("FsS"))->operator[]((NORMBIN == 3 ? Utility->GetJPsiBin(q2Bins) : Utility->GetPsiPBin(q2Bins)));
      SetValueAndErrors(EffPDFnorm,"FsS",1.0,&myString,&value,&errLo,&errHi);

      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("AsS"))->operator[]((NORMBIN == 3 ? Utility->GetJPsiBin(q2Bins) : Utility->GetPsiPBin(q2Bins)));
      SetValueAndErrors(EffPDFnorm,"AsS",1.0,&myString,&value,&errLo,&errHi);

      PrintVariables(EffPDFnorm->getVariables(),"vars");
      EffPDFintegral = EffPDFnorm->createIntegral(RooArgSet(*y,*z));
      effPsi = EffPDFintegral->getVal();
      
      PrintVariables(EffPDFnorm->getVariables(),"vars");
      cout << "\n@@@ Integral of S*E over (theta_l,theta_K) for normalization channel: " << effPsi << " @@@" << endl;

      delete EffPDFnorm;
    }


  fileFitResults << "@@@@@@ Measurements vs q^2 @@@@@@" << endl;
  for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
    {
      if (!((FitType >= 43) && (FitType <= 45)) && !((FitType >= 53) && (FitType <= 55)) &&
	  !((FitType >= 63) && (FitType <= 65)) && !((FitType >= 73) && (FitType <= 75)) &&
	  (Utility->ValIsInPsi(q2Bins,(q2Bins->operator[](i+1)+q2Bins->operator[](i))/2.) == true))
	{
	  // ###################################
	  // # Save zeros into prarameter file #
	  // ###################################
	  vecParStr = SaveFitResults(NULL,i,fitParam,configParam,vecConstr);
	  Utility->SaveFitValues(PARAMETERFILEOUT,vecParStr,i);
	  delete vecParStr;

	  continue;
	}


      // ######################
      // # Make external text #
      // ######################
      extText[i] = new TPaveText(0.63,0.5,0.97,0.63,"NDC");
      extText[i]->AddText(Form("%s%.2f%s%.2f%s","q#lower[0.4]{^{2}}: ",q2Bins->operator[](i)," #font[122]{\55} ",q2Bins->operator[](i+1)," GeV#kern[-0.3]{#lower[0.4]{^{2}}}"));
      extText[i]->SetTextAlign(11);
      extText[i]->SetBorderSize(0.0);
      extText[i]->SetFillStyle(0);
      extText[i]->SetTextSize(0.05);


      fileFitResults << "\nBin[" << i << "]: " << q2Bins->operator[](i) << " <= q^2 < " << q2Bins->operator[](i+1) << endl;
      
      myString.clear(); myString.str("");
      myString << "c_" << i;
      cq2Bins[i] = new TCanvas(myString.str().c_str(), myString.str().c_str(), 20, 20, 1800, 700);
  
      myString.clear(); myString.str("");
      myString << "(mumuMass*mumuMass) > " << q2Bins->operator[](i) << " && (mumuMass*mumuMass) <= " << q2Bins->operator[](i+1);
      cout << "Cut string: " << myString.str() << endl;
      dataSet_q2Bins[i] = (RooDataSet*)dataSet->reduce(myString.str().c_str());

      myString.clear(); myString.str("");
      myString << "TotalPDFq2Bin_" << i;
      InstantiateMassAngleFit(&TotalPDFq2Bins[i],useEffPDF,x,fitVar,myString.str(),FitType,configParam,fitParam,(ID == 0 ? i : 0),effFuncs->operator[](i),true);


      // #####################
      // # Initialize p.d.f. #
      // #####################
      CopyFitResults(TotalPDFq2Bins[i],(ID == 0 ? i : 0),fitParam);


      // #####################
      // # Apply constraints #
      // #####################
      ClearVars(vecConstr);
      if (!((FitType >= 33) && (FitType <= 35)) &&
	  !((FitType >= 53) && (FitType <= 55)) &&
	  !((FitType >= 73) && (FitType <= 75)))
	{
	  BuildMassConstraints(vecConstr,TotalPDFq2Bins[i],"All");
	  if (configParam->operator[](Utility->GetConfigParamIndx("fitPSIintru"))->operator[](i) != 1) BuildAngularConstraints(vecConstr,TotalPDFq2Bins[i]);
	  if ((FitType == 4) || (FitType == 44) || (FitType == 64) || (FitType == 5) || (FitType == 45) || (FitType == 65)) BuildPhysicsConstraints(vecConstr,TotalPDFq2Bins[i],"FlS");
	}
      if ((!((FitType >= 43) && (FitType <= 45)) && !((FitType >= 53) && (FitType <= 55)) &&
	   !((FitType >= 63) && (FitType <= 65)) && !((FitType >= 73) && (FitType <= 75))) ||
	  ((FitType >= 23) && (FitType <= 25) && (specBin != 3) && (specBin != 5)))
	{
	  BuildPhysicsConstraints(vecConstr,TotalPDFq2Bins[i],"FsS");
	  BuildPhysicsConstraints(vecConstr,TotalPDFq2Bins[i],"AsS");
	}
      BuildPhysicsConstraints(vecConstr,TotalPDFq2Bins[i],"AfbBounds"); // # Special constraint for AFB boundaries #
      BuildPhysicsConstraints(vecConstr,TotalPDFq2Bins[i],"AsBounds");  // # Special constraint for AS  boundaries #


      // ###################################
      // # Print variables and constraints #
      // ###################################
      PrintVariables(TotalPDFq2Bins[i]->getVariables(),"vars");
      PrintVariables(vecConstr,"cons");


      // ###################
      // # Perform the fit #
      // ###################
      fitResult = MakeMassAngleFit(dataSet_q2Bins[i],&TotalPDFq2Bins[i],x,fitVar,cq2Bins[i],FitType,configParam,(ID == 0 ? i : 0),vecConstr,&NLLvalue,extText[i],ID);


      // ##############################################
      // # Save fit results back into prarameter file #
      // ##############################################
      vecParStr = SaveFitResults(TotalPDFq2Bins[i],(ID == 0 ? i : 0),fitParam,configParam,vecConstr);
      Utility->SaveFitValues(PARAMETERFILEOUT,vecParStr,(ID == 0 ? i : 0));
      delete vecParStr;


      if (useEffPDF == true)
	{
	  if (EffPDFsign != NULL) delete EffPDFsign;
	  // ####################################################
	  // # Compute integral of S*E for signal and its error #
	  // ####################################################

	  Vars.removeAll();
	  myString.clear(); myString.str("");
	  myString << MakeAngWithEffPDF(effFuncs->operator[](i),NULL,y,z,1,useEffPDF,&Vars);
	  EffPDFsign = new RooGenericPdf("Efficiency",myString.str().c_str(),Vars);

	  myString.clear(); myString.str("");
	  myString << fitParam->operator[](Utility->GetFitParamIndx("FlS"))->operator[](i);
	  SetValueAndErrors(EffPDFsign,"FlS",1.0,&myString,&value,&errLo,&errHi);

	  myString.clear(); myString.str("");
	  myString << fitParam->operator[](Utility->GetFitParamIndx("AfbS"))->operator[](i);
	  SetValueAndErrors(EffPDFsign,"AfbS",1.0,&myString,&value,&errLo,&errHi);

	  myString.clear(); myString.str("");
	  myString << fitParam->operator[](Utility->GetFitParamIndx("FsS"))->operator[](i);
	  SetValueAndErrors(EffPDFsign,"FsS",1.0,&myString,&value,&errLo,&errHi);

	  myString.clear(); myString.str("");
	  myString << fitParam->operator[](Utility->GetFitParamIndx("AsS"))->operator[](i);
	  SetValueAndErrors(EffPDFsign,"AsS",1.0,&myString,&value,&errLo,&errHi);

	  PrintVariables(EffPDFsign->getVariables(),"vars");
	  EffPDFintegral = EffPDFsign->createIntegral(RooArgSet(*y,*z));
	  effMuMu = EffPDFintegral->getVal();

	  PrintVariables(EffPDFsign->getVariables(),"vars");
	  cout << "\n@@@ Integral of S*E over (theta_l,theta_K) for signal: " << effMuMu << " @@@" << endl;
	}


      // ########################################
      // # Save outcome of the fit in histogram #
      // ########################################
      if (configParam->operator[](Utility->GetConfigParamIndx("fitPSIintru"))->operator[](i) != 1)
	{
	  if ((FitType == 3) || (FitType == 33) || (FitType == 43) || (FitType == 53) || (FitType == 63) || (FitType == 73)) // Fl-Afb-fit
	    {
	      VecHistoMeas->operator[](0)->SetBinContent(i+1,GetVar(TotalPDFq2Bins[i],"FlS")->getVal());
	      VecHistoMeas->operator[](0)->SetBinError(i+1,GetVar(TotalPDFq2Bins[i],"FlS")->getError());
	      
	      VecHistoMeas->operator[](1)->SetBinContent(i+1,(GetVar(TotalPDFq2Bins[i],"AfbS") != NULL ? GetVar(TotalPDFq2Bins[i],"AfbS")->getVal() : 0.0));
	      VecHistoMeas->operator[](1)->SetBinContent(i+1,(GetVar(TotalPDFq2Bins[i],"AfbS") != NULL ? GetVar(TotalPDFq2Bins[i],"AfbS")->getError() : 0.0));
	    }
	  else if ((FitType == 4) || (FitType == 34) || (FitType == 44) || (FitType == 54) || (FitType == 64) || (FitType == 74)) // Afb-Fl-fit
	    {
	      VecHistoMeas->operator[](0)->SetBinContent(i+1,GetVar(TotalPDFq2Bins[i],"FlS")->getVal());
	      VecHistoMeas->operator[](0)->SetBinError(i+1,GetVar(TotalPDFq2Bins[i],"FlS")->getError());
	      
	      VecHistoMeas->operator[](1)->SetBinContent(i+1,GetVar(TotalPDFq2Bins[i],"AfbS")->getVal());
	      VecHistoMeas->operator[](1)->SetBinError(i+1,GetVar(TotalPDFq2Bins[i],"AfbS")->getError());
	    }
	  else if ((FitType == 5) || (FitType == 35) || (FitType == 45) || (FitType == 55) || (FitType == 65) || (FitType == 75)) // Fl-At2-Atim-fit
	    {
	      VecHistoMeas->operator[](0)->SetBinContent(i+1,GetVar(TotalPDFq2Bins[i],"FlS")->getVal());
	      VecHistoMeas->operator[](0)->SetBinError(i+1,GetVar(TotalPDFq2Bins[i],"FlS")->getError());
	      
	      VecHistoMeas->operator[](1)->SetBinContent(i+1,GetVar(TotalPDFq2Bins[i],"At2S")->getVal());
	      VecHistoMeas->operator[](1)->SetBinError(i+1,GetVar(TotalPDFq2Bins[i],"At2S")->getError());
	      
	      VecHistoMeas->operator[](2)->SetBinContent(i+1,GetVar(TotalPDFq2Bins[i],"AtimS")->getVal());
	      VecHistoMeas->operator[](2)->SetBinError(i+1,GetVar(TotalPDFq2Bins[i],"AtimS")->getError());
	    }

	  if (!((FitType >= 33) && (FitType <= 35)) && !((FitType >= 53) && (FitType <= 55)) && !((FitType >= 73) && (FitType <= 75)))
	    {
	      VecHistoMeas->operator[](3)->SetBinContent(i+1,(GetVar(TotalPDFq2Bins[i],"nSig")->getVal() / effMuMu) / (PsiYield / effPsi) * (NORMBIN == 3 ? Utility->JPsiBF : Utility->PsiPBF) / (q2Bins->operator[](i+1) - q2Bins->operator[](i)) / 1e-7);
	      VecHistoMeas->operator[](3)->SetBinError(i+1,VecHistoMeas->operator[](3)->GetBinContent(i+1) * sqrt(pow(GetVar(TotalPDFq2Bins[i],"nSig")->getError() / GetVar(TotalPDFq2Bins[i],"nSig")->getVal(),2.) + pow(PsiYieldErr / PsiYield,2.)));
	      
	      fileFitResults << "@@@@@@ dBF/dq^2 @@@@@@" << endl;
	      fileFitResults << "dBF/dq^2: " << VecHistoMeas->operator[](3)->GetBinContent(i+1) << " -/+ " << VecHistoMeas->operator[](3)->GetBinError(i+1) << " (";
	      fileFitResults << VecHistoMeas->operator[](3)->GetBinContent(i+1) * sqrt(pow(GetVar(TotalPDFq2Bins[i],"nSig")->getErrorLo() / GetVar(TotalPDFq2Bins[i],"nSig")->getVal(),2.) + pow(PsiYieldErr / PsiYield,2.)) << "/";
	      fileFitResults << VecHistoMeas->operator[](3)->GetBinContent(i+1) * sqrt(pow(GetVar(TotalPDFq2Bins[i],"nSig")->getErrorHi() / GetVar(TotalPDFq2Bins[i],"nSig")->getVal(),2.) + pow(PsiYieldErr / PsiYield,2.)) << ")" << endl;
	      
	      fileFitResults << "To cut and paste in config. file" << endl;
	      fileFitResults << VecHistoMeas->operator[](3)->GetBinContent(i+1) << "   ";
	      fileFitResults << VecHistoMeas->operator[](3)->GetBinContent(i+1) * sqrt(pow(GetVar(TotalPDFq2Bins[i],"nSig")->getErrorLo() / GetVar(TotalPDFq2Bins[i],"nSig")->getVal(),2.) + pow(PsiYieldErr / PsiYield,2.)) << "   ";
	      fileFitResults << VecHistoMeas->operator[](3)->GetBinContent(i+1) * sqrt(pow(GetVar(TotalPDFq2Bins[i],"nSig")->getErrorHi() / GetVar(TotalPDFq2Bins[i],"nSig")->getVal(),2.) + pow(PsiYieldErr / PsiYield,2.)) << endl;
	      fileFitResults << "====================================================================" << endl;
	    }

	  myString.clear(); myString.str("");
	  if (CheckGoodFit(fitResult) == true)
	    {
	      myString << ID << "   ";
	      myString << (((useEffPDF == true) && (GetVar(TotalPDFq2Bins[i],"FlS")  != NULL)) ? GetVar(TotalPDFq2Bins[i],"FlS")->getVal()  : -2.0) << "   ";
	      myString << (((useEffPDF == true) && (GetVar(TotalPDFq2Bins[i],"AfbS") != NULL)) ? GetVar(TotalPDFq2Bins[i],"AfbS")->getVal() : -2.0) << "   ";
	      myString << VecHistoMeas->operator[](3)->GetBinContent(i+1) << "   ";
	      myString << NLLvalue;
	    }
	  else myString << ID << "   " << -2.0 << "   " << -2.0 << "   " << -2.0 << "   " << -2.0;
	  
	  fileFitSystematics << myString.str() << endl;
	}
    }
}


void MakeMassAngleToy (RooAbsPdf* TotalPDF, RooRealVar* x, RooRealVar* y, TCanvas* Canv, unsigned int FitType, unsigned int nToy, vector<vector<unsigned int>*>* configParam, int specBin, vector<vector<string>*>* fitParam, RooArgSet* vecConstr, string fileName)
{
  unsigned int nEntryToy;
  stringstream myString;
  unsigned int it = 1;
  RooArgSet* constrToys;
  RooRealVar* tmpVar;
  RooPlot* myFrame;


  // #####################
  // # Initialize p.d.f. #
  // #####################
  nEntryToy = CopyFitResults(TotalPDF,specBin,fitParam);

  
  // #####################
  // # Apply constraints #
  // #####################
  ClearVars(vecConstr);
  BuildMassConstraints(vecConstr,TotalPDF,"All");
  BuildAngularConstraints(vecConstr,TotalPDF);
  if ((FitType == 24) || (FitType == 25)) BuildPhysicsConstraints(vecConstr,TotalPDF,"FlS");
  if ((specBin != 3) && (specBin != 5))
    {
      BuildPhysicsConstraints(vecConstr,TotalPDF,"FsS");
      BuildPhysicsConstraints(vecConstr,TotalPDF,"AsS");
    }
  constrToys = (RooArgSet*)vecConstr->clone("constrToys");
  BuildPhysicsConstraints(vecConstr,TotalPDF,"AfbBound"); // # Special constraint for AFB boundaries #
  BuildPhysicsConstraints(vecConstr,TotalPDF,"AsBound");  // # Special constraint for AS  boundaries #


  // ###################################
  // # Print variables and constraints #
  // ###################################
  PrintVariables(TotalPDF->getVariables(),"vars");
  PrintVariables(vecConstr,"cons");
  PrintVariables(constrToys,"cons");


  // #############################
  // # Toy-MC generation and fit #
  // #############################
  RooMCStudy* MyToy = new RooMCStudy(*TotalPDF,RooArgList(*x,*y),Extended(true),ExternalConstraints(*constrToys),FitOptions(Extended(true),ExternalConstraints(*vecConstr),Minos(USEMINOS)));
  MyToy->generateAndFit(nToy,nEntryToy,true);
  Canv->Divide(5,4);


  if ((GetVar(TotalPDF,"meanS") != NULL) && (GetVar(TotalPDF,"meanS")->getError() != 0.0) && (IsInConstraints(vecConstr,"meanS") == false))
    {
      RooRealVar* tmpVar = GetVar(TotalPDF,"meanS");
      RooPlot* myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("meanS"))->operator[](specBin).c_str()) -
				       0.4*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](specBin).c_str())),
				       atof(fitParam->operator[](Utility->GetFitParamIndx("meanS"))->operator[](specBin).c_str()) +
				       0.4*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"sigmaS1") != NULL) && (GetVar(TotalPDF,"sigmaS1")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaS1") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaS1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](specBin).c_str()) -
			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](specBin).c_str())),
			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](specBin).c_str()) +
			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"sigmaS2") != NULL) && (GetVar(TotalPDF,"sigmaS2")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaS2") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaS2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS2"))->operator[](specBin).c_str()) -
			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS2"))->operator[](specBin).c_str())),
			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS2"))->operator[](specBin).c_str()) +
			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"fracMassS") != NULL) && (GetVar(TotalPDF,"fracMassS")->getError() != 0.0) && (IsInConstraints(vecConstr,"fracMassS") == false))
    {
      tmpVar = GetVar(TotalPDF,"fracMassS");
      myFrame = tmpVar->frame(-0.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"tau1") != NULL) && (GetVar(TotalPDF,"tau1")->getError() != 0.0) && (IsInConstraints(vecConstr,"tau1") == false))
    {
      tmpVar = GetVar(TotalPDF,"tau1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("tau1"))->operator[](specBin).c_str()) - 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("tau1"))->operator[](specBin).c_str())),
			      atof(fitParam->operator[](Utility->GetFitParamIndx("tau1"))->operator[](specBin).c_str()) + 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("tau1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"tau2") != NULL) && (GetVar(TotalPDF,"tau2")->getError() != 0.0) && (IsInConstraints(vecConstr,"tau2") == false))
    {
      tmpVar = GetVar(TotalPDF,"tau2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("tau2"))->operator[](specBin).c_str()) - 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("tau2"))->operator[](specBin).c_str())),
			      atof(fitParam->operator[](Utility->GetFitParamIndx("tau2"))->operator[](specBin).c_str()) + 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("tau2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"fracMassBExp") != NULL) && (GetVar(TotalPDF,"fracMassBExp")->getError() != 0.0) && (IsInConstraints(vecConstr,"fracMassBExp") == false))
    {
      tmpVar = GetVar(TotalPDF,"fracMassBExp");
      myFrame = tmpVar->frame(-0.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"meanR1") != NULL) && (GetVar(TotalPDF,"meanR1")->getError() != 0.0) && (IsInConstraints(vecConstr,"meanR1") == false))
    {
      tmpVar = GetVar(TotalPDF,"meanR1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("meanR1"))->operator[](specBin).c_str()) - 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](specBin).c_str())),
			      atof(fitParam->operator[](Utility->GetFitParamIndx("meanR1"))->operator[](specBin).c_str()) + 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"sigmaR1") != NULL) && (GetVar(TotalPDF,"sigmaR1")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaR1") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaR1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](specBin).c_str()) -
			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](specBin).c_str())),
			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](specBin).c_str()) +
			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"meanR2") != NULL) && (GetVar(TotalPDF,"meanR2")->getError() != 0.0) && (IsInConstraints(vecConstr,"meanR2") == false))
    {
      tmpVar = GetVar(TotalPDF,"meanR2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("meanR2"))->operator[](specBin).c_str()) - 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](specBin).c_str())),
			      atof(fitParam->operator[](Utility->GetFitParamIndx("meanR2"))->operator[](specBin).c_str()) + 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    } 
   
  if ((GetVar(TotalPDF,"sigmaR2") != NULL) && (GetVar(TotalPDF,"sigmaR2")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaR2") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaR2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](specBin).c_str()) -
			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](specBin).c_str())),
			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](specBin).c_str()) +
			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }
  
  if ((GetVar(TotalPDF,"fracMassBRPeak") != NULL) && (GetVar(TotalPDF,"fracMassBRPeak")->getError() != 0.0) && (IsInConstraints(vecConstr,"fracMassBRPeak") == false))
    {
      tmpVar = GetVar(TotalPDF,"fracMassBRPeak");
      myFrame = tmpVar->frame(-0.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"meanL1") != NULL) && (GetVar(TotalPDF,"meanL1")->getError() != 0.0) && (IsInConstraints(vecConstr,"meanL1") == false))
    {
      tmpVar = GetVar(TotalPDF,"meanL1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("meanL1"))->operator[](specBin).c_str()) - 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](specBin).c_str())),
			      atof(fitParam->operator[](Utility->GetFitParamIndx("meanL1"))->operator[](specBin).c_str()) + 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"sigmaL1") != NULL) && (GetVar(TotalPDF,"sigmaL1")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaL1") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaL1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](specBin).c_str()) -
			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](specBin).c_str())),
			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](specBin).c_str()) +
			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"meanL2") != NULL) && (GetVar(TotalPDF,"meanL2")->getError() != 0.0) && (IsInConstraints(vecConstr,"meanL2") == false))
    {
      tmpVar = GetVar(TotalPDF,"meanL2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("meanL2"))->operator[](specBin).c_str()) - 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](specBin).c_str())),
			      atof(fitParam->operator[](Utility->GetFitParamIndx("meanL2"))->operator[](specBin).c_str()) + 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"sigmaL2") != NULL) && (GetVar(TotalPDF,"sigmaL2")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaL2") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaL2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](specBin).c_str()) -
			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](specBin).c_str())),
			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](specBin).c_str()) +
			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"fracMassBLPeak") != NULL) && (GetVar(TotalPDF,"fracMassBLPeak")->getError() != 0.0) && (IsInConstraints(vecConstr,"fracMassBLPeak") == false))
    {
      tmpVar = GetVar(TotalPDF,"fracMassBLPeak");
      myFrame = tmpVar->frame(-0.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"fracMassBPeak") != NULL) && (GetVar(TotalPDF,"fracMassBPeak")->getError() != 0.0) && (IsInConstraints(vecConstr,"fracMassBPeak") == false))
    {
      tmpVar = GetVar(TotalPDF,"fracMassBPeak");
      myFrame = tmpVar->frame(-0.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }
  
  if ((GetVar(TotalPDF,"nBkgPeak") != NULL) && (GetVar(TotalPDF,"nBkgPeak")->getError() != 0.0) && (IsInConstraints(vecConstr,"nBkgPeak") == false))
    {
      tmpVar = GetVar(TotalPDF,"nBkgPeak");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgPeak"))->operator[](specBin).c_str()) -
			      0.6*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgPeak"))->operator[](specBin).c_str())),
			      atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgPeak"))->operator[](specBin).c_str()) +
			      0.6*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgPeak"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"nSig") != NULL) && (GetVar(TotalPDF,"nSig")->getError() != 0.0) && (IsInConstraints(vecConstr,"nSig") == false))
    {
      tmpVar = GetVar(TotalPDF,"nSig");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("nSig"))->operator[](specBin).c_str()) - 0.6*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("nSig"))->operator[](specBin).c_str())),
			      atof(fitParam->operator[](Utility->GetFitParamIndx("nSig"))->operator[](specBin).c_str()) + 0.6*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("nSig"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"nBkgComb") != NULL) && (GetVar(TotalPDF,"nBkgComb")->getError() != 0.0) && (IsInConstraints(vecConstr,"nBkgComb") == false))
    {
      tmpVar = GetVar(TotalPDF,"nBkgComb");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgComb"))->operator[](specBin).c_str()) -
			      0.6*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgComb"))->operator[](specBin).c_str())),
			      atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgComb"))->operator[](specBin).c_str()) +
			      0.6*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgComb"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
    {
      myString.clear(); myString.str("");
      myString << "c1Poly" << i;

      if ((GetVar(TotalPDF,myString.str().c_str()) != NULL) && (GetVar(TotalPDF,myString.str().c_str())->getError() != 0.0) && (IsInConstraints(vecConstr,myString.str().c_str()) == false))
	{
	  tmpVar = GetVar(TotalPDF,myString.str().c_str());
	  myFrame = tmpVar->frame(-6.0,6.0);
	  MyToy->plotParamOn(myFrame);
	  Canv->cd(it++);
	  myFrame->Draw();
	  myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
	  Canv->cd(it++);
	  myFrame->Draw();
	}
    }
  
  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
    {
      myString.clear(); myString.str("");
      myString << "c2Poly" << i;

      if ((GetVar(TotalPDF,myString.str().c_str()) != NULL) && (GetVar(TotalPDF,myString.str().c_str())->getError() != 0.0) && (IsInConstraints(vecConstr,myString.str().c_str()) == false))
	{
	  tmpVar = GetVar(TotalPDF,myString.str().c_str());
	  myFrame = tmpVar->frame(-6.0,6.0);
	  MyToy->plotParamOn(myFrame);
	  Canv->cd(it++);
	  myFrame->Draw();
	  myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
	  Canv->cd(it++);
	  myFrame->Draw();
	}
    }

  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
    {
      myString.clear(); myString.str("");
      myString << "p1Poly" << i;
      
      if ((GetVar(TotalPDF,myString.str().c_str()) != NULL) && (GetVar(TotalPDF,myString.str().c_str())->getError() != 0.0) && (IsInConstraints(vecConstr,myString.str().c_str()) == false))
	{
	  tmpVar = GetVar(TotalPDF,myString.str().c_str());
	  myFrame = tmpVar->frame(-6.0,6.0);
	  MyToy->plotParamOn(myFrame);
	  Canv->cd(it++);
	  myFrame->Draw();
	  myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
	  Canv->cd(it++);
	  myFrame->Draw();
	}
    }
  
  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
    {
      myString.clear(); myString.str("");
      myString << "p2Poly" << i;
      
      if ((GetVar(TotalPDF,myString.str().c_str()) != NULL) && (GetVar(TotalPDF,myString.str().c_str())->getError() != 0.0) && (IsInConstraints(vecConstr,myString.str().c_str()) == false))
	{
	  tmpVar = GetVar(TotalPDF,myString.str().c_str());
	  myFrame = tmpVar->frame(-6.0,6.0);
	  MyToy->plotParamOn(myFrame);
	  Canv->cd(it++);
	  myFrame->Draw();
	  myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
	  Canv->cd(it++);
	  myFrame->Draw();
	}
    }

  if ((GetVar(TotalPDF,"FlS") != NULL) && (GetVar(TotalPDF,"FlS")->getError() != 0.0) && (IsInConstraints(vecConstr,"FlS") == false))
    {
      tmpVar = GetVar(TotalPDF,"FlS");
      myFrame = tmpVar->frame(-0.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"AfbS") != NULL) && (GetVar(TotalPDF,"AfbS")->getError() != 0.0) && (IsInConstraints(vecConstr,"AfbS") == false))
    {
      tmpVar = GetVar(TotalPDF,"AfbS");
      myFrame = tmpVar->frame(-1.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }
  
  if ((GetVar(TotalPDF,"At2S") != NULL) && (GetVar(TotalPDF,"At2S")->getError() != 0.0) && (IsInConstraints(vecConstr,"At2S") == false))
    {
      tmpVar = GetVar(TotalPDF,"At2S");
      myFrame = tmpVar->frame(-1.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"AtimS") != NULL) && (GetVar(TotalPDF,"AtimS")->getError() != 0.0) && (IsInConstraints(vecConstr,"AtimS") == false))
    {
      tmpVar = GetVar(TotalPDF,"AtimS");
      myFrame = tmpVar->frame(-1.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"FsS") != NULL) && (GetVar(TotalPDF,"FsS")->getError() != 0.0) && (IsInConstraints(vecConstr,"FsS") == false))
    {
      tmpVar = GetVar(TotalPDF,"FsS");
      myFrame = tmpVar->frame(-1.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"AsS") != NULL) && (GetVar(TotalPDF,"AsS")->getError() != 0.0) && (IsInConstraints(vecConstr,"AsS") == false))
    {
      tmpVar = GetVar(TotalPDF,"AsS");
      myFrame = tmpVar->frame(-1.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  myFrame = MyToy->plotNLL(-2000.0,0.0);
  Canv->cd(it++);
  myFrame->Draw();
  Canv->Update();


  // #############################################
  // # Re-make all the fits and save the results #
  // #############################################
  TH1D* histoDiff1 = new TH1D("histoDiff1","histoDiff1",50,-1.0,1.0);
  histoDiff1->SetFillColor(kAzure+6);
  histoDiff1->SetXTitle("(fit #font[122]{\55} pdf)");
  histoDiff1->SetYTitle("Entries [#]");

  TH1D* histoPull1 = new TH1D("histoPull1","histoPull1",30,-5.0,5.0);
  histoPull1->SetFillColor(kAzure+6);
  histoPull1->SetXTitle("(fit #font[122]{\55} pdf) / #sigma");
  histoPull1->SetYTitle("Entries [#]");

  TH1D* histoNLL1 = new TH1D("histoNLL1","histoNLL1",30,1.0,-1.0);
  histoNLL1->SetFillColor(kGreen-7);
  histoNLL1->SetXTitle("NLL");
  histoNLL1->SetYTitle("Entries [#]");

  TCanvas* cB0Toy = new TCanvas("cB0Toy","cB0Toy", 20, 20, 1800, 700);

  cout << "\n@@@ Now fit total TOY invariant mass and angle @@@" << endl;
  RooDataSet* toySample;
  RooFitResult* fitResult;
  double NLLvalue;

  string varName = "";
  if      (FitType == 23) varName = "FlS";
  else if (FitType == 24) varName = "AfbS";
  for (unsigned int i = 0; i < nToy; i++)
    {
      cout << "\n@@@ Now fitting toy #" << i << " @@@\n" << endl;

      toySample = (RooDataSet*)MyToy->genData(i);
      CopyFitResults(TotalPDF,specBin,fitParam);
      fitResult = MakeMassAngleFit(toySample,&TotalPDF,x,y,cB0Toy,FitType,configParam,specBin,vecConstr,&NLLvalue,NULL,i+1);


      // ######################################################
      // # To verify the outcome of the RooFit toy-MC studies #
      // ######################################################
      if ((CheckGoodFit(fitResult) == true) && (GetVar(TotalPDF,varName.c_str()) != NULL))
      	{
      	  myString.clear(); myString.str("");
      	  myString << fitParam->operator[](Utility->GetFitParamIndx(varName.c_str()))->operator[](specBin);
      	  if (TotalPDF->getVariables()->getRealValue(varName.c_str()) > atof(myString.str().c_str()))
      	    histoPull1->Fill((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(myString.str().c_str())) / fabs(GetVar(TotalPDF,varName.c_str())->getErrorLo()));
      	  else histoPull1->Fill((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(myString.str().c_str())) / fabs(GetVar(TotalPDF,varName.c_str())->getErrorHi()));
	  histoDiff1->Fill(TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(myString.str().c_str()));
      	  histoNLL1->Fill(NLLvalue);
      	}
      delete fitResult;
    }


  TCanvas* cNLL1 = new TCanvas("cNLL1","cNLL1",10,10,1000,500);
  cNLL1->Divide(3,1);
  cNLL1->cd(1);
  histoDiff1->Draw();
  cNLL1->cd(2);
  histoNLL1->Draw();
  cNLL1->cd(3);
  histoPull1->Draw();
  cNLL1->Update();


  // ##############
  // # Save plots #
  // ##############
  if (SAVEPLOT == true)
    {
      TFile* fNLL;

      Canv->Print(fileName.c_str());

      myString.clear(); myString.str("");
      myString << "_DIFF.root";
      fNLL = new TFile(fileName.replace(fileName.find(".root"),5,myString.str()).c_str(),"RECREATE");
      fNLL->cd();
      histoDiff1->Write();
      fNLL->Close();
      delete fNLL;

      myString.clear(); myString.str("");
      myString << "_PULL.root";
      fNLL = new TFile(fileName.replace(fileName.find("_DIFF.root"),10,myString.str()).c_str(),"RECREATE");
      fNLL->cd();
      histoPull1->Write();
      fNLL->Close();
      delete fNLL;

      myString.clear(); myString.str("");
      myString << "_NLL.root";
      fNLL = new TFile(fileName.replace(fileName.find("_PULL.root"),10,myString.str()).c_str(),"RECREATE");
      fNLL->cd();
      histoNLL1->Write();
      fNLL->Close();
      delete fNLL;
    }
}




// ==================
// ===> 3D MODEL <===
// ==================


void InstantiateMass2AnglesFit (RooAbsPdf** TotalPDF,
				bool useEffPDF,
				RooRealVar* x, RooRealVar* y, RooRealVar* z,
				string fitName, unsigned int FitType,
				vector<vector<unsigned int>*>* configParam,
				vector<vector<string>*>* fitParam,
				unsigned int parIndx,
				TF2* effFunc,
				bool isFirstDataset)
// #########################
// # x: mass               #
// # y: angle cos(theta_l) #
// # z: angle cos(theta_K) #
// #########################
{
  // ################################
  // # Read configuration variables #
  // ################################
  unsigned int nGaussPSIintru = configParam->operator[](Utility->GetConfigParamIndx("nGaussPSIintru"))->operator[](parIndx);
  unsigned int fitPSIintru    = configParam->operator[](Utility->GetConfigParamIndx("fitPSIintru"))->operator[](parIndx);
  bool use2GaussS             = configParam->operator[](Utility->GetConfigParamIndx("SigType"))->operator[](parIndx);
  bool use2BExp               = configParam->operator[](Utility->GetConfigParamIndx("BkgType"))->operator[](parIndx);

  // ###########################
  // # Read polynomial degrees #
  // ###########################
  unsigned int NCoeffPolyBKGpeak1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP1"))->operator[](parIndx).c_str());
  unsigned int NCoeffPolyBKGcomb1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC1"))->operator[](parIndx).c_str());
  unsigned int NCoeffPolyBKGpeak2 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP2"))->operator[](parIndx).c_str());
  unsigned int NCoeffPolyBKGcomb2 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC2"))->operator[](parIndx).c_str());

  stringstream myString;
  RooArgSet VarsC1, VarsC2;
  RooArgSet VarsP1, VarsP2;


  if (isFirstDataset == true)
    {
      // ################################################
      // # Define mass fit variables and pdf for signal #
      // ################################################
      meanS   = new RooRealVar("meanS","Signal mean",Utility->B0Mass,"GeV");

      sigmaS1 = new RooRealVar("sigmaS1","Signal sigma-1",0.0,"GeV");
      MassS1  = new RooGaussian("MassS1","Signal Gaussian-1",*x,*meanS,*sigmaS1);
      
      sigmaS2 = new RooRealVar("sigmaS2","Signal sigma-2",0.0,"GeV");
      MassS2  = new RooGaussian("MassS2","Signal Gaussian-2",*x,*meanS,*sigmaS2);

      meanS->setConstant(false);
      sigmaS1->setConstant(false);
      sigmaS2->setConstant(false);

      fracMassS = new RooRealVar("fracMassS","Fraction of signal Gaussian",0.0,0.0,1.0);
      fracMassS->setConstant(false);

      if (use2GaussS == true) MassSignal = new RooAddPdf("MassSignal","Signal mass pdf",RooArgList(*MassS1,*MassS2),RooArgList(*fracMassS));
      else                    MassSignal = new RooGaussian(*((RooGaussian*)MassS1),"MassSignal");


      // #################################################
      // # Define angle fit variables and pdf for signal #
      // #################################################
      FlS  = new RooRealVar("FlS","F_{L}",0.0,0.0,1.0);
      AfbS = new RooRealVar("AfbS","A_{FB}",0.0,-1.0,1.0);
      FlS->setConstant(false);
      AfbS->setConstant(false);

      RooArgSet* Vars = new RooArgSet("VarsEffPDF");
      myString.clear(); myString.str("");
      myString << MakeAngWithEffPDF(effFunc,NULL,y,z,FitType,useEffPDF,Vars);
      AngleS = new RooGenericPdf("AngleS",myString.str().c_str(),*Vars);

      Signal = new RooProdPdf("Signal","Signal Mass*Angle",RooArgList(*MassSignal,*AngleS));
    }


  // ###################################################################
  // # Define angle fit variables and pdf for combinatorial background #
  // ###################################################################
  for (unsigned int i = 0; i < NCoeffPolyBKGcomb1; i++)
    {
      myString.clear(); myString.str("");
      myString << "c1Poly" << i;
      c1Poly[i] = new RooRealVar(myString.str().c_str(),"Comb.bkg.poly.coef.",0.0);
      c1Poly[i]->setConstant(false);
      VarsC1.add(*c1Poly[i]);
    }
  for (unsigned int i = 0; i < NCoeffPolyBKGcomb2; i++)
    {
      myString.clear(); myString.str("");
      myString << "c2Poly" << i;
      c2Poly[i] = new RooRealVar(myString.str().c_str(),"Comb.bkg.poly.coef.",0.0);
      c2Poly[i]->setConstant(false);
      VarsC2.add(*c2Poly[i]);
    }
  BkgAngleC1 = new RooPolynomial("BkgAngleC1","Comb.bkg angular shape",*y,VarsC1);
  BkgAngleC2 = new RooPolynomial("BkgAngleC2","Comb.bkg angular shape",*z,VarsC2);
  BkgAnglesC = new RooProdPdf("BkgAnglesC","Background Angle1*Angle2",RooArgList(*BkgAngleC1,*BkgAngleC2));


  // ##################################################################
  // # Define mass fit variables and pdf for combinatorial background #
  // ##################################################################
  tau1 = new RooRealVar("tau1","First background time constant",0.0,"GeV");
  myString.clear(); myString.str("");
  myString << "exp(-(" << x->getPlotLabel() << "-meanS)/tau1)";
  BkgMassExp1 = new RooGenericPdf("BkgMassExp1",myString.str().c_str(),RooArgSet(*x,*meanS,*tau1));

  tau2 = new RooRealVar("tau2","Second background time constant",0.0,"GeV");
  myString.clear(); myString.str("");
  myString << "exp(-(" << x->getPlotLabel() << "-meanS)/tau2)";
  BkgMassExp2 = new RooGenericPdf("BkgMassExp2",myString.str().c_str(),RooArgSet(*x,*meanS,*tau2));

  tau1->setConstant(false);
  tau2->setConstant(false);
  
  fracMassBExp = new RooRealVar("fracMassBExp","Fraction of background Exponential",0.0,0.0,1.0);
  fracMassBExp->setConstant(false);

  if (use2BExp == true) BkgMassComb = new RooAddPdf("BkgMassComb","Background mass comb. bkg pdf",RooArgList(*BkgMassExp1,*BkgMassExp2),RooArgList(*fracMassBExp));
  else                  BkgMassComb = new RooGenericPdf(*((RooGenericPdf*)BkgMassExp1),"BkgMassComb");
  BkgMassAngleComb = new RooProdPdf("BkgMassAngleComb","Comb bkg Mass*Angle",RooArgList(*BkgMassComb,*BkgAnglesC));


  // #############################################################
  // # Define angle fit variables and pdf for peaking background #
  // #############################################################
  for (unsigned int i = 0; i < NCoeffPolyBKGpeak1; i++)
    {
      myString.clear(); myString.str("");
      myString << "p1Poly" << i;
      p1Poly[i] = new RooRealVar(myString.str().c_str(),"Peak.bkg.poly.coef.",0.0);
      p1Poly[i]->setConstant(false);
      VarsP1.add(*p1Poly[i]);
    }
  for (unsigned int i = 0; i < NCoeffPolyBKGpeak2; i++)
    {
      myString.clear(); myString.str("");
      myString << "p2Poly" << i;
      p2Poly[i] = new RooRealVar(myString.str().c_str(),"Peak.bkg.poly.coef.",0.0);
      p2Poly[i]->setConstant(false);
      VarsP2.add(*p2Poly[i]);
    }
  BkgAngleP1 = new RooPolynomial("BkgAngleP1","Peak.bkg angular shape",*y,VarsP1);
  BkgAngleP2 = new RooPolynomial("BkgAngleP2","Peak.bkg angular shape",*z,VarsP2);
  BkgAnglesP = new RooProdPdf("BkgAnglesP","Background Angle1*Angle2",RooArgList(*BkgAngleP1,*BkgAngleP2));


  // ############################################################
  // # Define mass fit variables and pdf for peaking background #
  // ############################################################
  meanR1         = new RooRealVar("meanR1","Bkg right peak mean-1",0.0,"GeV");
  sigmaR1        = new RooRealVar("sigmaR1","Bkg right peak sigma-1",0.0,"GeV");
  BkgMassRPeak1  = new RooGaussian("BkgMassRPeak1","Bkg right peak-1",*x,*meanR1,*sigmaR1);

  meanR2         = new RooRealVar("meanR2","Bkg right peak mean-2",0.0,"GeV");
  sigmaR2        = new RooRealVar("sigmaR2","Bkg right peak sigma-2",0.0,"GeV");
  BkgMassRPeak2  = new RooGaussian("BkgMassRPeak2","Bkg right peak-2",*x,*meanR2,*sigmaR2);

  fracMassBRPeak = new RooRealVar("fracMassBRPeak","Fraction of background right Peak",0.0,0.0,1.0);
  BkgMassRPeak   = new RooAddPdf("BkgMassRPeak","Right peaking bkg mass pdf",RooArgList(*BkgMassRPeak1,*BkgMassRPeak2),RooArgList(*fracMassBRPeak));

  meanL1         = new RooRealVar("meanL1","Bkg left peak mean-1",0.0,"GeV");
  sigmaL1        = new RooRealVar("sigmaL1","Bkg left peak sigma-1",0.0,"GeV");
  BkgMassLPeak1  = new RooGaussian("BkgMassLPeak1","Bkg left peak-1",*x,*meanL1,*sigmaL1);

  meanL2         = new RooRealVar("meanL2","Bkg left peak mean-2",0.0,"GeV");
  sigmaL2        = new RooRealVar("sigmaL2","Bkg left peak sigma-2",0.0,"GeV");
  BkgMassLPeak2  = new RooGaussian("BkgMassLPeak2","Bkg left peak-2",*x,*meanL2,*sigmaL2);

  fracMassBLPeak = new RooRealVar("fracMassBLPeak","Fraction of background left Peak",0.0,0.0,1.0);
  BkgMassLPeak   = new RooAddPdf("BkgMassLPeak","Left peaking bkg mass pdf",RooArgList(*BkgMassLPeak1,*BkgMassLPeak2),RooArgList(*fracMassBLPeak));

  fracMassBPeak  = new RooRealVar("fracMassBPeak","Fraction of background right-left Peak",0.0,0.0,1.0);
  
  meanR1->setConstant(false);
  sigmaR1->setConstant(false);
  meanR2->setConstant(false);
  sigmaR2->setConstant(false);
  fracMassBRPeak->setConstant(false);
  
  meanL1->setConstant(false);
  sigmaL1->setConstant(false);
  meanL2->setConstant(false);
  sigmaL2->setConstant(false);
  fracMassBLPeak->setConstant(false);
  
  fracMassBPeak->setConstant(false);

  if ((fitPSIintru == 1) || (fitPSIintru == 2))
    if (nGaussPSIintru != 1) BkgMassPeak = new RooAddPdf(*((RooAddPdf*)BkgMassRPeak),"BkgMassPeak");
    else      	             BkgMassPeak = new RooGaussian(*((RooGaussian*)BkgMassRPeak1),"BkgMassPeak");
  else
    if (nGaussPSIintru != 1) BkgMassPeak = new RooAddPdf("BkgMassPeak","Peaking bkg mass pdf",RooArgList(*BkgMassRPeak,*BkgMassLPeak),RooArgList(*fracMassBPeak));
    else                     BkgMassPeak = new RooAddPdf("BkgMassPeak","Peaking bkg mass pdf",RooArgList(*BkgMassRPeak1,*BkgMassLPeak1),RooArgList(*fracMassBPeak));
  BkgMassAnglePeak = new RooProdPdf("BkgMassAnglePeak","Peaking bkg Mass*Angle",RooArgList(*BkgMassPeak,*BkgAnglesP));


  // ###########################
  // # Define pdf coefficients #
  // ###########################
  nSig = new RooRealVar("nSig","Number of signal events",1.0);
  nBkgComb = new RooRealVar("nBkgComb","Number of comb. background events",1.0);
  nBkgPeak = new RooRealVar("nBkgPeak","Number of peaking background events",1.0);

  nSig->setConstant(false);
  nBkgComb->setConstant(false);
  nBkgPeak->setConstant(false);

  if ((FitType == 36) || (FitType == 56) || (FitType == 76))
                             *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgList(*AngleS),RooArgList(*nSig));
  else if (fitPSIintru == 0) *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgList(*Signal,*BkgMassAngleComb),RooArgList(*nSig,*nBkgComb));
  else if (fitPSIintru == 1) *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgList(*BkgMassAnglePeak),RooArgList(*nBkgPeak));
  else                       *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgList(*Signal,*BkgMassAngleComb,*BkgMassAnglePeak),RooArgList(*nSig,*nBkgComb,*nBkgPeak));
}


RooFitResult* MakeMass2AnglesFit (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* x, RooRealVar* y, RooRealVar* z, TCanvas* Canv, unsigned int FitType, vector<vector<unsigned int>*>* configParam, unsigned int parIndx, RooArgSet* vecConstr, double* NLLvalue, TPaveText* extText, int ID)
{
  // ################################
  // # Read configuration variables #
  // ################################
  unsigned int fitPSIintru = configParam->operator[](Utility->GetConfigParamIndx("fitPSIintru"))->operator[](parIndx);

  RooFitResult* fitResult = NULL;
  stringstream myString;
  double signalSigma  = 0.0;
  double signalSigmaE = 0.0;
  int nElements = 4;
  if (fitPSIintru == 1) nElements = 2;
  else if (GetVar(*TotalPDF,"nBkgPeak") != NULL) nElements = 5;
  TCanvas* localCanv[4];
  localCanv[0] = new TCanvas("localCanv0","localCanv0",20,20,700,500);
  localCanv[1] = new TCanvas("localCanv1","localCanv1",20,20,700,500);
  localCanv[2] = new TCanvas("localCanv2","localCanv2",20,20,700,500);
  localCanv[3] = new TCanvas("localCanv3","localCanv3",20,20,700,500);


  // ###########################
  // # Divide the input canvas #
  // ###########################
  Canv->Divide(3,3);


  if ((FitType == 36) || (FitType == 56) || (FitType == 76))
    {
      if (MAKEGRAPHSCAN == true) MakeGraphicalScan(*TotalPDF,Canv,y,z);
      else
	{
	  // ###################
	  // # Make actual fit #
	  // ###################
	  if (ApplyConstr == true) fitResult = (*TotalPDF)->fitTo(*dataSet,Extended(true),ExternalConstraints(*vecConstr),Save(true),Minos(USEMINOS));
	  else                     fitResult = (*TotalPDF)->fitTo(*dataSet,Extended(true),Save(true),Minos(USEMINOS));


	  // ##################################################
	  // # Set p.d.f independent variables to known point #
	  // ##################################################
	  if (GetVar(*TotalPDF,x->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(x->getPlotLabel(),Utility->B0Mass);
	  if (GetVar(*TotalPDF,y->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(y->getPlotLabel(),0.0);
	  if (GetVar(*TotalPDF,z->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(z->getPlotLabel(),0.0);
	  if (fitResult != NULL) fitResult->Print("v");


	  // ##########################
	  // # Angular-1 plot results #
	  // ##########################
	  Canv->cd(2);
	  RooPlot* myFrameY = y->frame(NBINS);
	  dataSet->plotOn(myFrameY,Name(MakeName(dataSet,ID).c_str()));
	  if (FUNCERRBAND == true) (*TotalPDF)->plotOn(myFrameY,Name((*TotalPDF)->getPlotLabel()),VisualizeError(*fitResult,1,true),VLines(),FillColor(kGreen-7));
	  else                     (*TotalPDF)->plotOn(myFrameY,Name((*TotalPDF)->getPlotLabel()));

	  (*TotalPDF)->paramOn(myFrameY,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(RooArgSet(*nSig)));
      
	  myString.clear(); myString.str("");
	  myString << (*TotalPDF)->getPlotLabel() << "_paramBox";
	  TPaveText* paveTextY = (TPaveText*)myFrameY->findObject(myString.str().c_str());
	  paveTextY->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));

	  TLegend* legY = new TLegend(0.75, 0.65, 0.97, 0.89, "");
	  for (int i = 0; i < 2; i++)
	    {
	      TString objName = myFrameY->nameOf(i);
	      if (objName == "") continue;
	      TObject* obj = myFrameY->findObject(objName.Data());
	      legY->AddEntry(obj,legNames[i],"lp");
	    }
	  legY->SetFillStyle(0);
	  legY->SetFillColor(0);
	  legY->SetTextSize(0.05);
	  legY->SetBorderSize(0);
 

	  // ##########################
	  // # Angular-2 plot results #
	  // ##########################
	  Canv->cd(3);
	  RooPlot* myFrameZ = z->frame(NBINS);
	  dataSet->plotOn(myFrameZ,Name(MakeName(dataSet,ID).c_str()));
	  if (FUNCERRBAND == true) (*TotalPDF)->plotOn(myFrameZ,Name((*TotalPDF)->getPlotLabel()),VisualizeError(*fitResult,1,true),VLines(),FillColor(kGreen-7));
	  else                     (*TotalPDF)->plotOn(myFrameZ,Name((*TotalPDF)->getPlotLabel()));

	  (*TotalPDF)->paramOn(myFrameZ,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(RooArgSet(*nSig)));
      
	  myString.clear(); myString.str("");
	  myString << (*TotalPDF)->getPlotLabel() << "_paramBox";
	  TPaveText* paveTextZ = (TPaveText*)myFrameZ->findObject(myString.str().c_str());
	  paveTextZ->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));

	  TLegend* legZ = new TLegend(0.75, 0.65, 0.97, 0.89, "");
	  for (int i = 0; i < 2; i++)
	    {
	      TString objName = myFrameZ->nameOf(i);
	      if (objName == "") continue;
	      TObject* obj = myFrameZ->findObject(objName.Data());
	      legZ->AddEntry(obj,legNames[i],"lp");
	    }
	  legZ->SetFillStyle(0);
	  legZ->SetFillColor(0);
	  legZ->SetTextSize(0.05);
	  legZ->SetBorderSize(0);


	  // ####################
	  // # Save fit results #
	  // ####################
	  fileFitResults << "====================================================================" << endl;
	  fileFitResults << "@@@@@@ Make angle fit @@@@@@" << endl;
 
	  fileFitResults << "Chi2/DoF y-fit: " << myFrameY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str());
	  fileFitResults << "; p-value: " << TMath::Prob(myFrameY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())*NBINS,NBINS) << endl;
	  fileFitResults << "Chi2/DoF z-fit: " << myFrameZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str());
	  fileFitResults << "; p-value: " << TMath::Prob(myFrameZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())*NBINS,NBINS) << endl;

	  *NLLvalue = StoreFitResultsInFile(TotalPDF,fitResult,dataSet,vecConstr);
	  
	  fileFitResults << "====================================================================" << endl;


	  // #######################
	  // # Add NLL to the plot #
	  // #######################
	  Canv->cd(2);
	  paveTextY->AddText(Form("%s%.1f","NLL = ",*NLLvalue));
	  CheckGoodFit(fitResult,paveTextY);
	  paveTextY->SetTextAlign(11);
	  paveTextY->SetBorderSize(0.0);
	  paveTextY->SetFillStyle(0);
	  paveTextY->SetTextSize(0.04);
	  paveTextY->Paint();
	  myFrameY->Draw();
	  legY->Draw("same");


	  // #######################
	  // # Add NLL to the plot #
	  // #######################
	  Canv->cd(3);
	  paveTextZ->AddText(Form("%s%.1f","NLL = ",*NLLvalue));
	  CheckGoodFit(fitResult,paveTextZ);
	  paveTextZ->SetTextAlign(11);
	  paveTextZ->SetBorderSize(0.0);
	  paveTextZ->SetFillStyle(0);
	  paveTextZ->SetTextSize(0.04);
	  paveTextZ->Paint();
	  myFrameZ->Draw();
	  legZ->Draw("same");
	}
    }
  else
    {
      // #####################
      // # Make sideband fit #
      // #####################
      if (fitPSIintru != 1)
	{
	  RooAbsPdf* TmpPDF = NULL;
	  RooDataSet* sideBands = NULL;
	  RooRealVar frac("frac","Fraction",0.5,0.0,1.0);
	  RooArgSet constrSidebads;
	  ClearVars(&constrSidebads);
	  BuildAngularConstraints(&constrSidebads,*TotalPDF);


	  // ################
	  // # Save results #
	  // ################
	  fileFitResults << "====================================================================" << endl;
	  fileFitResults << "@@@@@@ B0 mass sideband fit @@@@@@" << endl;
	  fileFitResults << "Amplitude of signal region (+/- n*< Sigma >): " << Utility->GetGenericParam("NSigmaB0S") << " * " << Utility->GetB0Width() << endl;


	  // ##############
	  // # Get p.d.f. #
	  // ##############
	  if (fitPSIintru != 0) TmpPDF = new RooAddPdf("TmpPDF","Temporary p.d.f.",RooArgList(*BkgAnglesC,*BkgAnglesP),RooArgList(frac));
	  else                  TmpPDF = new RooProdPdf(*((RooProdPdf*)BkgAnglesC),"TmpPDF");


	  // #############
	  // # Sidebands #
	  // #############
	  myString.clear(); myString.str("");
	  myString << "B0MassArb < " << (*TotalPDF)->getVariables()->getRealValue("meanS") - Utility->GetGenericParam("NSigmaB0S")*Utility->GetB0Width();
	  myString << " || B0MassArb > " << (*TotalPDF)->getVariables()->getRealValue("meanS") + Utility->GetGenericParam("NSigmaB0S")*Utility->GetB0Width();
	  cout << "Cut for B0 sidebands: " << myString.str().c_str() << endl;
	  sideBands = (RooDataSet*)dataSet->reduce(myString.str().c_str());


	  // ###################
	  // # Make actual fit #
	  // ###################
	  if (ApplyConstr == true) fitResult = TmpPDF->fitTo(*sideBands,ExternalConstraints(constrSidebads),Save(true));
	  else                     fitResult = TmpPDF->fitTo(*sideBands,Save(true));
	  if (fitResult != NULL) fitResult->Print("v");

	  
	  // ####################
	  // # Save fit results #
	  // ####################
	  StorePolyResultsInFile(TotalPDF);


	  delete TmpPDF;
	  delete sideBands;
	  ClearVars(&constrSidebads);
	  if (fitResult != NULL) delete fitResult;
	}


      // ###################
      // # Make actual fit #
      // ###################
      if (ApplyConstr == true) fitResult = (*TotalPDF)->fitTo(*dataSet,Extended(true),ExternalConstraints(*vecConstr),Save(true),Minos(USEMINOS));
      else                     fitResult = (*TotalPDF)->fitTo(*dataSet,Extended(true),Save(true),Minos(USEMINOS));


      // ##################################################
      // # Set p.d.f independent variables to known point #
      // ##################################################
      if (GetVar(*TotalPDF,x->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(x->getPlotLabel(),Utility->B0Mass);
      if (GetVar(*TotalPDF,y->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(y->getPlotLabel(),0.0);
      if (GetVar(*TotalPDF,z->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(z->getPlotLabel(),0.0);
      if (fitResult != NULL) fitResult->Print("v");


      if (fitPSIintru != 1)
	{
	  if (GetVar(*TotalPDF,"fracMassS") != NULL)
	    {
	      signalSigma  = sqrt((*TotalPDF)->getVariables()->getRealValue("fracMassS") * pow((*TotalPDF)->getVariables()->getRealValue("sigmaS1"),2.) +
				  (1.-(*TotalPDF)->getVariables()->getRealValue("fracMassS")) * pow((*TotalPDF)->getVariables()->getRealValue("sigmaS2"),2.));
	      signalSigmaE = 1./(2.*signalSigma) * sqrt( pow((pow((*TotalPDF)->getVariables()->getRealValue("sigmaS1"),2.)-pow((*TotalPDF)->getVariables()->getRealValue("sigmaS2"),2.)) * GetVar(*TotalPDF,"fracMassS")->getError(),2.) +
							 pow(2.*(*TotalPDF)->getVariables()->getRealValue("fracMassS") * (*TotalPDF)->getVariables()->getRealValue("sigmaS1")*GetVar(*TotalPDF,"sigmaS1")->getError(),2.) +
							 pow(2.*(1.-(*TotalPDF)->getVariables()->getRealValue("fracMassS")) * (*TotalPDF)->getVariables()->getRealValue("sigmaS2")*GetVar(*TotalPDF,"sigmaS2")->getError(),2.) );
	    }
	  else if (GetVar(*TotalPDF,"sigmaS1") != NULL)
	    {
	      signalSigma  = (*TotalPDF)->getVariables()->getRealValue("sigmaS1");
	      signalSigmaE = GetVar(*TotalPDF,"sigmaS1")->getError();
	    }
	}
      
      
      // #####################
      // # Mass plot results #
      // #####################
      Canv->cd(1);
      RooPlot* myFrameX = x->frame(NBINS);
      dataSet->plotOn(myFrameX,Name(MakeName(dataSet,ID).c_str()));
      if (FUNCERRBAND == true) (*TotalPDF)->plotOn(myFrameX,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),VisualizeError(*fitResult,1,true),VLines(),FillColor(kGreen-7));
      else                     (*TotalPDF)->plotOn(myFrameX,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack));
      if (fitPSIintru != 1)
	{
	  (*TotalPDF)->plotOn(myFrameX,Components(*Signal),LineStyle(7),LineColor(kBlue));
	  (*TotalPDF)->plotOn(myFrameX,Components(*BkgMassAngleComb),LineStyle(4),LineColor(kRed));
	}
      if (GetVar(*TotalPDF,"nBkgPeak") != NULL) (*TotalPDF)->plotOn(myFrameX,Components(*BkgMassAnglePeak),LineStyle(3),LineColor(kViolet));

      if      (fitPSIintru == 0) (*TotalPDF)->paramOn(myFrameX,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(RooArgSet(*nSig,*nBkgComb)));
      else if (fitPSIintru == 1) (*TotalPDF)->paramOn(myFrameX,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(RooArgSet(*nBkgPeak)));
      else                       (*TotalPDF)->paramOn(myFrameX,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(RooArgSet(*nSig,*nBkgComb,*nBkgPeak)),ShowConstants(true));
      
      myString.clear(); myString.str("");
      myString << (*TotalPDF)->getPlotLabel() << "_paramBox";
      TPaveText* paveTextX = (TPaveText*)myFrameX->findObject(myString.str().c_str());
      if (fitPSIintru != 1)
	{
	  paveTextX->AddText(Form("%s%.3f#pm%.3f","#mu = ",GetVar(*TotalPDF,"meanS")->getVal(),GetVar(*TotalPDF,"meanS")->getError()));
	  paveTextX->AddText(Form("%s%.4f#pm%.4f","< #sigma > = ",signalSigma,signalSigmaE));
	}
      paveTextX->AddText(Form("%s%.5f","FCN(m#lower[0.4]{^{B0}}#lower[-0.4]{_{PDG}},ang=0) = ",(*TotalPDF)->getVal()));
      paveTextX->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameX->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));

      TLegend* legX = new TLegend(0.75, 0.65, 0.97, 0.89, "");
      for (int i = 0; i < nElements; i++)
  	{
  	  TString objName = myFrameX->nameOf(i);
  	  if (objName == "") continue;
  	  TObject* obj = myFrameX->findObject(objName.Data());
	  legX->AddEntry(obj,legNames[i],"lp");
  	}
      legX->SetFillStyle(0);
      legX->SetFillColor(0);
      legX->SetTextSize(0.05);
      legX->SetBorderSize(0);


      // ##########################
      // # Angular-1 plot results #
      // ##########################
      Canv->cd(2);
      RooPlot* myFrameY = y->frame(NBINS);
      dataSet->plotOn(myFrameY,Name(MakeName(dataSet,ID).c_str()));
      if (FUNCERRBAND == true) (*TotalPDF)->plotOn(myFrameY,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),VisualizeError(*fitResult,1,true),VLines(),FillColor(kGreen-7));
      else                     (*TotalPDF)->plotOn(myFrameY,Name((*TotalPDF)->getPlotLabel()));
      if (fitPSIintru != 1)
	{
	  (*TotalPDF)->plotOn(myFrameY,Components(*Signal),LineStyle(7),LineColor(kBlue));
	  (*TotalPDF)->plotOn(myFrameY,Components(*BkgMassAngleComb),LineStyle(4),LineColor(kRed));
	}
      if (GetVar(*TotalPDF,"nBkgPeak") != NULL) (*TotalPDF)->plotOn(myFrameY,Components(*BkgMassAnglePeak),LineStyle(3),LineColor(kViolet));

      TPaveText* paveTextY = new TPaveText(0.11,0.8,0.4,0.86,"NDC");
      paveTextY->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
      paveTextY->SetTextAlign(11);
      paveTextY->SetBorderSize(0.0);
      paveTextY->SetFillStyle(0);
      paveTextY->SetTextSize(0.04);
      paveTextY->Paint();
      myFrameY->addObject(paveTextY);
      DrawString(LUMI,myFrameY);
      if ((extText != NULL) && (SETBATCH == false))
	{
	  extText->Paint();
	  myFrameY->addObject(extText);
	}
      myFrameY->Draw();

      TLegend* legY = new TLegend(0.75, 0.65, 0.97, 0.89, "");
      for (int i = 0; i < nElements; i++)
      	{
      	  TString objName = myFrameY->nameOf(i);
      	  if (objName == "") continue;
      	  TObject* obj = myFrameY->findObject(objName.Data());
      	  legY->AddEntry(obj,legNames[i],"lp");
      	}
      legY->SetFillStyle(0);
      legY->SetFillColor(0);
      legY->SetTextSize(0.05);
      legY->SetBorderSize(0);
      legY->Draw("same");
 

      // #####################
      // # Fill local canvas #
      // #####################
      localCanv[1]->cd();
      myFrameY->Draw();
      legY->Draw("same");


      // ##########################
      // # Angular-2 plot results #
      // ##########################
      Canv->cd(3);
      RooPlot* myFrameZ = z->frame(NBINS);
      dataSet->plotOn(myFrameZ,Name(MakeName(dataSet,ID).c_str()));
      if (FUNCERRBAND == true) (*TotalPDF)->plotOn(myFrameZ,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),VisualizeError(*fitResult,1,true),VLines(),FillColor(kGreen-7));
      else                     (*TotalPDF)->plotOn(myFrameZ,Name((*TotalPDF)->getPlotLabel()));
      if (fitPSIintru != 1)
      	{
      	  (*TotalPDF)->plotOn(myFrameZ,Components(*Signal),LineStyle(7),LineColor(kBlue));
      	  (*TotalPDF)->plotOn(myFrameZ,Components(*BkgMassAngleComb),LineStyle(4),LineColor(kRed));
      	}
      if (GetVar(*TotalPDF,"nBkgPeak") != NULL) (*TotalPDF)->plotOn(myFrameZ,Components(*BkgMassAnglePeak),LineStyle(3),LineColor(kViolet));

      TPaveText* paveTextZ = new TPaveText(0.11,0.8,0.4,0.86,"NDC");
      paveTextZ->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
      paveTextZ->SetTextAlign(11);
      paveTextZ->SetBorderSize(0.0);
      paveTextZ->SetFillStyle(0);
      paveTextZ->SetTextSize(0.04);
      paveTextZ->Paint();
      myFrameZ->addObject(paveTextZ);
      DrawString(LUMI,myFrameZ);
      if ((extText != NULL) && (SETBATCH == false))
      	{
      	  extText->Paint();
      	  myFrameZ->addObject(extText);
      	}
      myFrameZ->Draw();

      TLegend* legZ = new TLegend(0.75, 0.65, 0.97, 0.89, "");
      for (int i = 0; i < nElements; i++)
      	{
      	  TString objName = myFrameZ->nameOf(i);
      	  if (objName == "") continue;
      	  TObject* obj = myFrameZ->findObject(objName.Data());
      	  legZ->AddEntry(obj,legNames[i],"lp");
      	}
      legZ->SetFillStyle(0);
      legZ->SetFillColor(0);
      legZ->SetTextSize(0.05);
      legZ->SetBorderSize(0);
      legZ->Draw("same");


      // #####################
      // # Fill local canvas #
      // #####################
      localCanv[2]->cd();
      myFrameZ->Draw();
      legZ->Draw("same");


      // ################################
      // # LogLikelihood parameter scan #
      // ################################
      if ((FitType != 26) && (GetVar(*TotalPDF,"FlS") != NULL) && (GetVar(*TotalPDF,"AfbS") != NULL))
	{
	  localCanv[3]->cd();
	  RooAbsReal* NLL;
	  if (ApplyConstr == true) NLL = (*TotalPDF)->createNLL(*dataSet,Extended(true),ExternalConstraints(*vecConstr));
	  else                     NLL = (*TotalPDF)->createNLL(*dataSet,Extended(true));
	  RooMinuit RooMin(*NLL);
	  RooPlot* myFrameNLL = RooMin.contour(*GetVar(*TotalPDF,"AfbS"),*GetVar(*TotalPDF,"FlS"),1.0,2.0,3.0);
	  DrawString(LUMI,myFrameY);
	  myFrameNLL->Draw();
 	  delete NLL;
	}


      // ####################
      // # Save fit results #
      // ####################
      fileFitResults << "====================================================================" << endl;
      fileFitResults << "@@@@@@ Make mass - angle fit @@@@@@" << endl;

      fileFitResults << "Chi2/DoF x-fit: " << myFrameX->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str());
      fileFitResults << "; p-value: " << TMath::Prob(myFrameX->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())*NBINS,NBINS) << endl;
      fileFitResults << "Chi2/DoF y-fit: " << myFrameY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str());
      fileFitResults << "; p-value: " << TMath::Prob(myFrameY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())*NBINS,NBINS) << endl;
      fileFitResults << "Chi2/DoF z-fit: " << myFrameZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str());
      fileFitResults << "; p-value: " << TMath::Prob(myFrameZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())*NBINS,NBINS) << endl;

      *NLLvalue = StoreFitResultsInFile(TotalPDF,fitResult,dataSet,vecConstr);
      StorePolyResultsInFile(TotalPDF);
      fileFitResults << "====================================================================" << endl;


      // #######################
      // # Add NLL to the plot #
      // #######################
      Canv->cd(1);
      paveTextX->AddText(Form("%s%.1f","NLL = ",*NLLvalue));
      CheckGoodFit(fitResult,paveTextX);
      paveTextX->SetBorderSize(0.0);
      paveTextX->SetFillStyle(0);
      paveTextX->SetTextSize(0.03);
      paveTextX->Paint();
      DrawString(LUMI,myFrameX);
      if ((extText != NULL) && (SETBATCH == false))
      	{
      	  if (GetVar(*TotalPDF,"nSig") != NULL) extText->AddText(Form("%s%1.f%s%1.f","Signal yield: ",GetVar(*TotalPDF,"nSig")->getVal()," #pm ",GetVar(*TotalPDF,"nSig")->getError()));
      	  extText->Paint();
      	  myFrameX->addObject(extText);
      	}
      myFrameX->Draw();
      legX->Draw("same");


      // #####################
      // # Fill local canvas #
      // #####################
      localCanv[0]->cd();
      myFrameX->Draw();
      legX->Draw("same");


      // #####################
      // # Plot on sidebands #
      // #####################
      if ((fitPSIintru != 1) && (ID == 0))
	{
	  // ################
	  // # Save results #
	  // ################
	  fileFitResults << "@@@@@@ B0 mass sideband fit overlaps @@@@@@" << endl;
	  fileFitResults << "Amplitude of signal region: +/- " << Utility->GetGenericParam("NSigmaB0S") << " * < Sigma >" << endl;


	  // ################
	  // # Low sideband #
	  // ################
	  x->setRange("lowSideband",(*TotalPDF)->getVariables()->getRealValue("meanS") - Utility->GetGenericParam("B0MassIntervalLeft"),
		      (*TotalPDF)->getVariables()->getRealValue("meanS") - Utility->GetGenericParam("NSigmaB0S")*signalSigma);
	  y->setRange("lowSideband",y->getMin(),y->getMax());
	  z->setRange("lowSideband",z->getMin(),z->getMax());


	  // #################
	  // # Signal region #
	  // #################
	  x->setRange("signalRegion",(*TotalPDF)->getVariables()->getRealValue("meanS") - Utility->GetGenericParam("NSigmaB0S")*signalSigma,
		      (*TotalPDF)->getVariables()->getRealValue("meanS") + Utility->GetGenericParam("NSigmaB0S")*signalSigma);
	  y->setRange("signalRegion",y->getMin(),y->getMax());
	  z->setRange("signalRegion",z->getMin(),z->getMax());


	  // #################
	  // # High sideband #
	  // #################
	  x->setRange("highSideband",(*TotalPDF)->getVariables()->getRealValue("meanS") + Utility->GetGenericParam("NSigmaB0S")*signalSigma,
		      (*TotalPDF)->getVariables()->getRealValue("meanS") + Utility->GetGenericParam("B0MassIntervalRight"));
	  y->setRange("highSideband",y->getMin(),y->getMax());
	  z->setRange("highSideband",z->getMin(),z->getMax());


	  // ##############
	  // # Y-variable #
	  // ##############


	  // ###########################
	  // # Background plot results #
	  // ###########################
	  Canv->cd(4);
	  RooPlot* myFrameLowSideBY = y->frame(NBINS);
	  dataSet->plotOn(myFrameLowSideBY,Name(MakeName(dataSet,ID).c_str()),CutRange("lowSideband"));
	  (*TotalPDF)->plotOn(myFrameLowSideBY,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),ProjectionRange("lowSideband"));
	  if (fitPSIintru != 1)
	    {
	      (*TotalPDF)->plotOn(myFrameLowSideBY,Components(*Signal),LineStyle(7),LineColor(kBlue),ProjectionRange("lowSideband"));
	      (*TotalPDF)->plotOn(myFrameLowSideBY,Components(*BkgMassAngleComb),LineStyle(4),LineColor(kRed),ProjectionRange("lowSideband"));
	    }
	  if (GetVar(*TotalPDF,"nBkgPeak") != NULL) (*TotalPDF)->plotOn(myFrameLowSideBY,Components(*BkgMassAnglePeak),LineStyle(3),LineColor(kViolet),ProjectionRange("lowSideband"));

	  TPaveText* paveTextLowSideBY = new TPaveText(0.11,0.75,0.4,0.88,"NDC");
	  paveTextLowSideBY->AddText("Low sideband");
	  paveTextLowSideBY->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameLowSideBY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
	  paveTextLowSideBY->SetTextAlign(11);
	  paveTextLowSideBY->SetBorderSize(0.0);
	  paveTextLowSideBY->SetFillStyle(0);
	  paveTextLowSideBY->SetTextSize(0.04);
	  paveTextLowSideBY->Paint();
	  myFrameLowSideBY->addObject(paveTextLowSideBY);
	  DrawString(LUMI,myFrameLowSideBY);
	  myFrameLowSideBY->Draw();

	  TLegend* legLowSideBY = new TLegend(0.75, 0.65, 0.97, 0.89, "");
	  for (int i = 0; i < nElements; i++)
	    {
	      TString objName = myFrameLowSideBY->nameOf(i);
	      if (objName == "") continue;
	      TObject* obj = myFrameLowSideBY->findObject(objName.Data());
	      legLowSideBY->AddEntry(obj,legNames[i],"lp");
	    }
	  legLowSideBY->SetFillStyle(0);
	  legLowSideBY->SetFillColor(0);
	  legLowSideBY->SetTextSize(0.04);
	  legLowSideBY->SetBorderSize(0);
	  legLowSideBY->Draw("same");


	  // ################
	  // # Save results #
	  // ################
	  fileFitResults << "Chi2/DoF low sideband-Y: " << myFrameLowSideBY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str());
	  fileFitResults << "; p-value: " << TMath::Prob(myFrameLowSideBY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())*NBINS,NBINS) << endl;


	  // ###########################
	  // # Background plot results #
	  // ###########################
	  Canv->cd(5);
	  RooPlot* myFrameSignalRegionY = y->frame(NBINS);
	  dataSet->plotOn(myFrameSignalRegionY,Name(MakeName(dataSet,ID).c_str()),CutRange("signalRegion"));
	  (*TotalPDF)->plotOn(myFrameSignalRegionY,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),ProjectionRange("signalRegion"));
	  if (fitPSIintru != 1)
	    {
	      (*TotalPDF)->plotOn(myFrameSignalRegionY,Components(*Signal),LineStyle(7),LineColor(kBlue),ProjectionRange("signalRegion"));
	      (*TotalPDF)->plotOn(myFrameSignalRegionY,Components(*BkgMassAngleComb),LineStyle(4),LineColor(kRed),ProjectionRange("signalRegion"));
	    }
	  if (GetVar(*TotalPDF,"nBkgPeak") != NULL) (*TotalPDF)->plotOn(myFrameSignalRegionY,Components(*BkgMassAnglePeak),LineStyle(3),LineColor(kViolet),ProjectionRange("signalRegion"));

	  TPaveText* paveTextSignalRegionY = new TPaveText(0.11,0.75,0.4,0.88,"NDC");
	  paveTextSignalRegionY->AddText(Form("%s%.1f%s","Signal region: #pm",Utility->GetGenericParam("NSigmaB0S")," < #sigma >"));
	  paveTextSignalRegionY->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameSignalRegionY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
	  paveTextSignalRegionY->SetTextAlign(11);
	  paveTextSignalRegionY->SetBorderSize(0.0);
	  paveTextSignalRegionY->SetFillStyle(0);
	  paveTextSignalRegionY->SetTextSize(0.04);
	  paveTextSignalRegionY->Paint();
	  myFrameSignalRegionY->addObject(paveTextSignalRegionY);
	  DrawString(LUMI,myFrameSignalRegionY);
	  myFrameSignalRegionY->Draw();

	  TLegend* legSignalRegionY = new TLegend(0.75, 0.65, 0.97, 0.89, "");
	  for (int i = 0; i < nElements; i++)
	    {
	      TString objName = myFrameSignalRegionY->nameOf(i);
	      if (objName == "") continue;
	      TObject* obj = myFrameSignalRegionY->findObject(objName.Data());
	      legSignalRegionY->AddEntry(obj,legNames[i],"lp");
	    }
	  legSignalRegionY->SetFillStyle(0);
	  legSignalRegionY->SetFillColor(0);
	  legSignalRegionY->SetTextSize(0.04);
	  legSignalRegionY->SetBorderSize(0);
	  legSignalRegionY->Draw("same");


	  // ################
	  // # Save results #
	  // ################
	  fileFitResults << "Chi2/DoF signal region-Y: " << myFrameSignalRegionY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str());
	  fileFitResults << "; p-value: " << TMath::Prob(myFrameSignalRegionY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())*NBINS,NBINS) << endl;


	  // ###########################
	  // # Background plot results #
	  // ###########################
	  Canv->cd(6);
	  RooPlot* myFrameHighSideBY = y->frame(NBINS);
	  dataSet->plotOn(myFrameHighSideBY,Name(MakeName(dataSet,ID).c_str()),CutRange("highSideband"));
	  (*TotalPDF)->plotOn(myFrameHighSideBY,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),ProjectionRange("highSideband"));
	  if (fitPSIintru != 1)
	    {
	      (*TotalPDF)->plotOn(myFrameHighSideBY,Components(*Signal),LineStyle(7),LineColor(kBlue),ProjectionRange("highSideband"));
	      (*TotalPDF)->plotOn(myFrameHighSideBY,Components(*BkgMassAngleComb),LineStyle(4),LineColor(kRed),ProjectionRange("highSideband"));
	    }
	  if (GetVar(*TotalPDF,"nBkgPeak") != NULL) (*TotalPDF)->plotOn(myFrameHighSideBY,Components(*BkgMassAnglePeak),LineStyle(3),LineColor(kViolet),ProjectionRange("highSideband"));

	  TPaveText* paveTextHighSideBY = new TPaveText(0.11,0.75,0.4,0.88,"NDC");
	  paveTextHighSideBY->AddText("High sideband");
	  paveTextHighSideBY->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameHighSideBY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
	  paveTextHighSideBY->SetTextAlign(11);
	  paveTextHighSideBY->SetBorderSize(0.0);
	  paveTextHighSideBY->SetFillStyle(0);
	  paveTextHighSideBY->SetTextSize(0.04);
	  paveTextHighSideBY->Paint();
	  myFrameHighSideBY->addObject(paveTextHighSideBY);
	  DrawString(LUMI,myFrameHighSideBY);
	  myFrameHighSideBY->Draw();

	  TLegend* legHighSideBY = new TLegend(0.75, 0.65, 0.97, 0.89, "");
	  for (int i = 0; i < nElements; i++)
	    {
	      TString objName = myFrameHighSideBY->nameOf(i);
	      if (objName == "") continue;
	      TObject* obj = myFrameHighSideBY->findObject(objName.Data());
	      legHighSideBY->AddEntry(obj,legNames[i],"lp");
	    }
	  legHighSideBY->SetFillStyle(0);
	  legHighSideBY->SetFillColor(0);
	  legHighSideBY->SetTextSize(0.04);
	  legHighSideBY->SetBorderSize(0);
	  legHighSideBY->Draw("same");


	  // ################
	  // # Save results #
	  // ################
	  fileFitResults << "Chi2/DoF high sideband-Y: " << myFrameHighSideBY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str());
	  fileFitResults << "; p-value: " << TMath::Prob(myFrameHighSideBY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())*NBINS,NBINS) << endl;


	  // ##############
	  // # Z-variable #
	  // ##############


	  // ###########################
	  // # Background plot results #
	  // ###########################
	  Canv->cd(7);
	  RooPlot* myFrameLowSideBZ = z->frame(NBINS);
	  dataSet->plotOn(myFrameLowSideBZ,Name(MakeName(dataSet,ID).c_str()),CutRange("lowSideband"));
	  (*TotalPDF)->plotOn(myFrameLowSideBZ,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),ProjectionRange("lowSideband"));
	  if (fitPSIintru != 1)
	    {
	      (*TotalPDF)->plotOn(myFrameLowSideBZ,Components(*Signal),LineStyle(7),LineColor(kBlue),ProjectionRange("lowSideband"));
	      (*TotalPDF)->plotOn(myFrameLowSideBZ,Components(*BkgMassAngleComb),LineStyle(4),LineColor(kRed),ProjectionRange("lowSideband"));
	    }
	  if (GetVar(*TotalPDF,"nBkgPeak") != NULL) (*TotalPDF)->plotOn(myFrameLowSideBZ,Components(*BkgMassAnglePeak),LineStyle(3),LineColor(kViolet),ProjectionRange("lowSideband"));

	  TPaveText* paveTextLowSideBZ = new TPaveText(0.11,0.75,0.4,0.88,"NDC");
	  paveTextLowSideBZ->AddText("Low sideband");
	  paveTextLowSideBZ->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameLowSideBZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
	  paveTextLowSideBZ->SetTextAlign(11);
	  paveTextLowSideBZ->SetBorderSize(0.0);
	  paveTextLowSideBZ->SetFillStyle(0);
	  paveTextLowSideBZ->SetTextSize(0.04);
	  paveTextLowSideBZ->Paint();
	  myFrameLowSideBZ->addObject(paveTextLowSideBZ);
	  DrawString(LUMI,myFrameLowSideBZ);
	  myFrameLowSideBZ->Draw();

	  TLegend* legLowSideBZ = new TLegend(0.75, 0.65, 0.97, 0.89, "");
	  for (int i = 0; i < nElements; i++)
	    {
	      TString objName = myFrameLowSideBZ->nameOf(i);
	      if (objName == "") continue;
	      TObject* obj = myFrameLowSideBZ->findObject(objName.Data());
	      legLowSideBZ->AddEntry(obj,legNames[i],"lp");
	    }
	  legLowSideBZ->SetFillStyle(0);
	  legLowSideBZ->SetFillColor(0);
	  legLowSideBZ->SetTextSize(0.04);
	  legLowSideBZ->SetBorderSize(0);
	  legLowSideBZ->Draw("same");


	  // ################
	  // # Save results #
	  // ################
	  fileFitResults << "Chi2/DoF low sideband-Z: " << myFrameLowSideBZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str());
	  fileFitResults << "; p-value: " << TMath::Prob(myFrameLowSideBZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())*NBINS,NBINS) << endl;


	  // ###########################
	  // # Background plot results #
	  // ###########################
	  Canv->cd(8);
	  RooPlot* myFrameSignalRegionZ = z->frame(NBINS);
	  dataSet->plotOn(myFrameSignalRegionZ,Name(MakeName(dataSet,ID).c_str()),CutRange("signalRegion"));
	  (*TotalPDF)->plotOn(myFrameSignalRegionZ,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),ProjectionRange("signalRegion"));
	  if (fitPSIintru != 1)
	    {
	      (*TotalPDF)->plotOn(myFrameSignalRegionZ,Components(*Signal),LineStyle(7),LineColor(kBlue),ProjectionRange("signalRegion"));
	      (*TotalPDF)->plotOn(myFrameSignalRegionZ,Components(*BkgMassAngleComb),LineStyle(4),LineColor(kRed),ProjectionRange("signalRegion"));
	    }
	  if (GetVar(*TotalPDF,"nBkgPeak") != NULL) (*TotalPDF)->plotOn(myFrameSignalRegionZ,Components(*BkgMassAnglePeak),LineStyle(3),LineColor(kViolet),ProjectionRange("signalRegion"));

	  TPaveText* paveTextSignalRegionZ = new TPaveText(0.11,0.75,0.4,0.88,"NDC");
	  paveTextSignalRegionZ->AddText(Form("%s%.1f%s","Signal region: #pm",Utility->GetGenericParam("NSigmaB0S")," < #sigma >"));
	  paveTextSignalRegionZ->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameSignalRegionZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
	  paveTextSignalRegionZ->SetTextAlign(11);
	  paveTextSignalRegionZ->SetBorderSize(0.0);
	  paveTextSignalRegionZ->SetFillStyle(0);
	  paveTextSignalRegionZ->SetTextSize(0.04);
	  paveTextSignalRegionZ->Paint();
	  myFrameSignalRegionZ->addObject(paveTextSignalRegionZ);
	  DrawString(LUMI,myFrameSignalRegionZ);
	  myFrameSignalRegionZ->Draw();

	  TLegend* legSignalRegionZ = new TLegend(0.75, 0.65, 0.97, 0.89, "");
	  for (int i = 0; i < nElements; i++)
	    {
	      TString objName = myFrameSignalRegionZ->nameOf(i);
	      if (objName == "") continue;
	      TObject* obj = myFrameSignalRegionZ->findObject(objName.Data());
	      legSignalRegionZ->AddEntry(obj,legNames[i],"lp");
	    }
	  legSignalRegionZ->SetFillStyle(0);
	  legSignalRegionZ->SetFillColor(0);
	  legSignalRegionZ->SetTextSize(0.04);
	  legSignalRegionZ->SetBorderSize(0);
	  legSignalRegionZ->Draw("same");


	  // ################
	  // # Save results #
	  // ################
	  fileFitResults << "Chi2/DoF signal region-Z: " << myFrameSignalRegionZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str());
	  fileFitResults << "; p-value: " << TMath::Prob(myFrameSignalRegionZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())*NBINS,NBINS) << endl;


	  // ###########################
	  // # Background plot results #
	  // ###########################
	  Canv->cd(9);
	  RooPlot* myFrameHighSideBZ = z->frame(NBINS);
	  dataSet->plotOn(myFrameHighSideBZ,Name(MakeName(dataSet,ID).c_str()),CutRange("highSideband"));
	  (*TotalPDF)->plotOn(myFrameHighSideBZ,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),ProjectionRange("highSideband"));
	  if (fitPSIintru != 1)
	    {
	      (*TotalPDF)->plotOn(myFrameHighSideBZ,Components(*Signal),LineStyle(7),LineColor(kBlue),ProjectionRange("highSideband"));
	      (*TotalPDF)->plotOn(myFrameHighSideBZ,Components(*BkgMassAngleComb),LineStyle(4),LineColor(kRed),ProjectionRange("highSideband"));
	    }
	  if (GetVar(*TotalPDF,"nBkgPeak") != NULL) (*TotalPDF)->plotOn(myFrameHighSideBZ,Components(*BkgMassAnglePeak),LineStyle(3),LineColor(kViolet),ProjectionRange("highSideband"));

	  TPaveText* paveTextHighSideBZ = new TPaveText(0.11,0.75,0.4,0.88,"NDC");
	  paveTextHighSideBZ->AddText("High sideband");
	  paveTextHighSideBZ->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameHighSideBZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
	  paveTextHighSideBZ->SetTextAlign(11);
	  paveTextHighSideBZ->SetBorderSize(0.0);
	  paveTextHighSideBZ->SetFillStyle(0);
	  paveTextHighSideBZ->SetTextSize(0.04);
	  paveTextHighSideBZ->Paint();
	  myFrameHighSideBZ->addObject(paveTextHighSideBZ);
	  DrawString(LUMI,myFrameHighSideBZ);
	  myFrameHighSideBZ->Draw();

	  TLegend* legHighSideBZ = new TLegend(0.75, 0.65, 0.97, 0.89, "");
	  for (int i = 0; i < nElements; i++)
	    {
	      TString objName = myFrameHighSideBZ->nameOf(i);
	      if (objName == "") continue;
	      TObject* obj = myFrameHighSideBZ->findObject(objName.Data());
	      legHighSideBZ->AddEntry(obj,legNames[i],"lp");
	    }
	  legHighSideBZ->SetFillStyle(0);
	  legHighSideBZ->SetFillColor(0);
	  legHighSideBZ->SetTextSize(0.04);
	  legHighSideBZ->SetBorderSize(0);
	  legHighSideBZ->Draw("same");


	  // ################
	  // # Save results #
	  // ################
	  fileFitResults << "Chi2/DoF high sideband-Z: " << myFrameHighSideBZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str());
	  fileFitResults << "; p-value: " << TMath::Prob(myFrameHighSideBZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())*NBINS,NBINS) << endl;


	  fileFitResults << "====================================================================" << endl;
	}
    }
  

  // ##############
  // # Save plots #
  // ##############
  Canv->Update();
  if (SAVEPLOT == true)
    {
      myString.clear(); myString.str("");
      myString << (*TotalPDF)->getPlotLabel() << ".pdf";
      Canv->Print(myString.str().c_str());
    }

  if (SETBATCH == true)
    {
      delete localCanv[0];
      delete localCanv[1];
      delete localCanv[2];
      delete localCanv[3];
    }
  else
    {
      localCanv[0]->Update();
      localCanv[1]->Update();
      localCanv[2]->Update();
      localCanv[3]->Update();
    }


  return fitResult;
}


void IterativeMass2AnglesFitq2Bins (RooDataSet* dataSet,
				    bool useEffPDF,
				    double PsiYield, double PsiYieldErr,
				    RooRealVar* x, RooRealVar* y, RooRealVar* z,
				    int specBin,
				    unsigned int FitType,
				    vector<TH1D*>* VecHistoMeas,
				    vector<double>* q2Bins,
				    vector<vector<unsigned int>*>* configParam, vector<vector<string>*>* fitParam,
				    vector<TF2*>* effFuncs,
				    RooArgSet* vecConstr,
				    unsigned int ID)
{
  double value, errLo, errHi;

  vector<string>* vecParStr;
  stringstream myString;

  double effPsi  = 0.0;
  double effMuMu = 0.0;

  RooArgSet   Vars;
  TCanvas*    cq2Bins[q2Bins->size()-1];
  RooDataSet* dataSet_q2Bins[q2Bins->size()-1];
  RooAbsPdf*  TotalPDFq2Bins[q2Bins->size()-1];
  TPaveText*  extText[q2Bins->size()-1];

  RooFitResult* fitResult;
  RooAbsReal*   EffPDFintegral = NULL;
  RooAbsPdf*    EffPDFnorm     = NULL;
  RooAbsPdf*    EffPDFsign     = NULL;

  double NLLvalue;
  

  if ((useEffPDF == true) && (fitParam->operator[](0)->size() > 1))
    {
      // ##############################################################
      // # Compute integral of S*E for resonant channel and its error #
      // ##############################################################
 
      // ######################################
      // # Reference to normalization channel #
      // ######################################
      myString.clear(); myString.str("");
      myString << MakeAngWithEffPDF(effFuncs->operator[]((NORMBIN == 3 ? Utility->GetJPsiBin(q2Bins) : Utility->GetPsiPBin(q2Bins))),NULL,y,z,FitType,useEffPDF,&Vars);
      EffPDFnorm = new RooGenericPdf("Efficiency",myString.str().c_str(),Vars);
 
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("FlS"))->operator[]((NORMBIN == 3 ? Utility->GetJPsiBin(q2Bins) : Utility->GetPsiPBin(q2Bins)));
      SetValueAndErrors(EffPDFnorm,"FlS",1.0,&myString,&value,&errLo,&errHi);

      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("AfbS"))->operator[]((NORMBIN == 3 ? Utility->GetJPsiBin(q2Bins) : Utility->GetPsiPBin(q2Bins)));
      SetValueAndErrors(EffPDFnorm,"AfbS",1.0,&myString,&value,&errLo,&errHi);

      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("FsS"))->operator[]((NORMBIN == 3 ? Utility->GetJPsiBin(q2Bins) : Utility->GetPsiPBin(q2Bins)));
      SetValueAndErrors(EffPDFnorm,"FsS",1.0,&myString,&value,&errLo,&errHi);

      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("AsS"))->operator[]((NORMBIN == 3 ? Utility->GetJPsiBin(q2Bins) : Utility->GetPsiPBin(q2Bins)));
      SetValueAndErrors(EffPDFnorm,"AsS",1.0,&myString,&value,&errLo,&errHi);

      PrintVariables(EffPDFnorm->getVariables(),"vars");
      EffPDFintegral = EffPDFnorm->createIntegral(RooArgSet(*y,*z));
      effPsi = EffPDFintegral->getVal();
      
      PrintVariables(EffPDFnorm->getVariables(),"vars");
      cout << "\n@@@ Integral of S*E over (theta_l,theta_K) for normalization channel: " << effPsi << " @@@" << endl;

      delete EffPDFnorm;
    }


  fileFitResults << "@@@@@@ Measurements vs q^2 @@@@@@" << endl;
  for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
    {
      if ((FitType != 46) && (FitType != 56) &&
	  (FitType != 66) && (FitType != 76) &&
	  (Utility->ValIsInPsi(q2Bins,(q2Bins->operator[](i+1)+q2Bins->operator[](i))/2.) == true))
	{
	  // ###################################
	  // # Save zeros into prarameter file #
	  // ###################################
	  vecParStr = SaveFitResults(NULL,i,fitParam,configParam,vecConstr);
	  Utility->SaveFitValues(PARAMETERFILEOUT,vecParStr,i);
	  delete vecParStr;
	  
	  continue;
	}


      // ######################
      // # Make external text #
      // ######################
      extText[i] = new TPaveText(0.63,0.5,0.97,0.63,"NDC");
      extText[i]->AddText(Form("%s%.2f%s%.2f%s","q#lower[0.4]{^{2}}: ",q2Bins->operator[](i)," #font[122]{\55} ",q2Bins->operator[](i+1)," GeV#kern[-0.3]{#lower[0.4]{^{2}}}"));
      extText[i]->SetTextAlign(11);
      extText[i]->SetBorderSize(0.0);
      extText[i]->SetFillStyle(0);
      extText[i]->SetTextSize(0.05);


      fileFitResults << "\nBin[" << i << "]: " << q2Bins->operator[](i) << " <= q^2 < " << q2Bins->operator[](i+1) << endl;
      
      myString.clear(); myString.str("");
      myString << "c_" << i;
      cq2Bins[i] = new TCanvas(myString.str().c_str(), myString.str().c_str(), 20, 20, 1800, 1800);
  
      myString.clear(); myString.str("");
      myString << "(mumuMass*mumuMass) > " << q2Bins->operator[](i) << " && (mumuMass*mumuMass) <= " << q2Bins->operator[](i+1);
      cout << "Cut string: " << myString.str() << endl;
      dataSet_q2Bins[i] = (RooDataSet*)dataSet->reduce(myString.str().c_str());

      myString.clear(); myString.str("");
      myString << "TotalPDFq2Bin_" << i;
      InstantiateMass2AnglesFit(&TotalPDFq2Bins[i],useEffPDF,x,y,z,myString.str(),FitType,configParam,fitParam,(ID == 0 ? i : 0),effFuncs->operator[](i),true);
      

      // #####################
      // # Initialize p.d.f. #
      // #####################
      CopyFitResults(TotalPDFq2Bins[i],(ID == 0 ? i : 0),fitParam);


      // #####################
      // # Apply constraints #
      // #####################
      ClearVars(vecConstr);
      if ((FitType != 36) && (FitType != 56) && (FitType != 76))
	{
	  BuildMassConstraints(vecConstr,TotalPDFq2Bins[i],"All");
	  if (configParam->operator[](Utility->GetConfigParamIndx("fitPSIintru"))->operator[](i) != 1) BuildAngularConstraints(vecConstr,TotalPDFq2Bins[i]);
	}
      if (((FitType != 46) && (FitType != 56) &&
	   (FitType != 66) && (FitType != 76)) ||
	  ((FitType == 26) && (specBin != 3) && (specBin != 5)))
	{
	  BuildPhysicsConstraints(vecConstr,TotalPDFq2Bins[i],"FsS");
	  BuildPhysicsConstraints(vecConstr,TotalPDFq2Bins[i],"AsS");
	}
      BuildPhysicsConstraints(vecConstr,TotalPDFq2Bins[i],"AfbBound"); // # Special constraint for AFB boundaries #
      BuildPhysicsConstraints(vecConstr,TotalPDFq2Bins[i],"AsBound");  // # Special constraint for AS  boundaries #


      // ###################################
      // # Print variables and constraints #
      // ###################################
      PrintVariables(TotalPDFq2Bins[i]->getVariables(),"vars");
      PrintVariables(vecConstr,"cons");


      // ###################
      // # Perform the fit #
      // ###################
      fitResult = MakeMass2AnglesFit(dataSet_q2Bins[i],&TotalPDFq2Bins[i],x,y,z,cq2Bins[i],FitType,configParam,(ID == 0 ? i : 0),vecConstr,&NLLvalue,extText[i],ID);


      // ##############################################
      // # Save fit results back into prarameter file #
      // ##############################################
      vecParStr = SaveFitResults(TotalPDFq2Bins[i],(ID == 0 ? i : 0),fitParam,configParam,vecConstr);
      Utility->SaveFitValues(PARAMETERFILEOUT,vecParStr,(ID == 0 ? i : 0));
      delete vecParStr;


      if (useEffPDF == true)
	{
	  if (EffPDFsign != NULL) delete EffPDFsign;
	  // ####################################################
	  // # Compute integral of S*E for signal and its error #
	  // ####################################################

	  EffPDFsign = new RooGenericPdf(*((RooGenericPdf*)AngleS),"EffPDFsign");

	  PrintVariables(EffPDFsign->getVariables(),"vars");
	  EffPDFintegral = EffPDFsign->createIntegral(RooArgSet(*y,*z));
	  effMuMu = EffPDFintegral->getVal();

	  PrintVariables(EffPDFsign->getVariables(),"vars");
	  cout << "\n@@@ Integral of S*E over (theta_l,theta_K) for signal: " << effMuMu << " @@@" << endl;
	}

 
      // ########################################
      // # Save outcome of the fit in histogram #
      // ########################################
      if (configParam->operator[](Utility->GetConfigParamIndx("fitPSIintru"))->operator[](i) != 1)
	{
	  VecHistoMeas->operator[](0)->SetBinContent(i+1,GetVar(TotalPDFq2Bins[i],"FlS")->getVal());
	  VecHistoMeas->operator[](0)->SetBinError(i+1,GetVar(TotalPDFq2Bins[i],"FlS")->getError());
	  
	  VecHistoMeas->operator[](1)->SetBinContent(i+1,GetVar(TotalPDFq2Bins[i],"AfbS")->getVal());
	  VecHistoMeas->operator[](1)->SetBinError(i+1,GetVar(TotalPDFq2Bins[i],"AfbS")->getError());

	  if ((FitType != 36) && (FitType != 56) && (FitType != 76))
	    {
	      VecHistoMeas->operator[](2)->SetBinContent(i+1,(GetVar(TotalPDFq2Bins[i],"nSig")->getVal() / effMuMu) / (PsiYield / effPsi) * (NORMBIN == 3 ? Utility->JPsiBF : Utility->PsiPBF) / (q2Bins->operator[](i+1) - q2Bins->operator[](i)) / 1e-7);
	      VecHistoMeas->operator[](2)->SetBinError(i+1,VecHistoMeas->operator[](2)->GetBinContent(i+1) * sqrt(pow(GetVar(TotalPDFq2Bins[i],"nSig")->getError() / GetVar(TotalPDFq2Bins[i],"nSig")->getVal(),2.) + pow(PsiYieldErr / PsiYield,2.)));
	      
	      fileFitResults << "@@@@@@ dBF/dq^2 @@@@@@" << endl;
	      fileFitResults << "dBF/dq^2: " << VecHistoMeas->operator[](2)->GetBinContent(i+1) << " -/+ " << VecHistoMeas->operator[](2)->GetBinError(i+1) << " (";
	      fileFitResults << VecHistoMeas->operator[](2)->GetBinContent(i+1) * sqrt(pow(GetVar(TotalPDFq2Bins[i],"nSig")->getErrorLo() / GetVar(TotalPDFq2Bins[i],"nSig")->getVal(),2.) +pow(PsiYieldErr / PsiYield,2.)) << "/";
	      fileFitResults << VecHistoMeas->operator[](2)->GetBinContent(i+1) * sqrt(pow(GetVar(TotalPDFq2Bins[i],"nSig")->getErrorHi() / GetVar(TotalPDFq2Bins[i],"nSig")->getVal(),2.) + pow(PsiYieldErr / PsiYield,2.)) << ")" << endl;
	      
	      fileFitResults << "To cut and paste in config. file" << endl;
	      fileFitResults << VecHistoMeas->operator[](2)->GetBinContent(i+1) << "   ";
	      fileFitResults << VecHistoMeas->operator[](2)->GetBinContent(i+1) * sqrt(pow(GetVar(TotalPDFq2Bins[i],"nSig")->getErrorLo() / GetVar(TotalPDFq2Bins[i],"nSig")->getVal(),2.) +pow(PsiYieldErr / PsiYield,2.)) << "   ";
	      fileFitResults << VecHistoMeas->operator[](2)->GetBinContent(i+1) * sqrt(pow(GetVar(TotalPDFq2Bins[i],"nSig")->getErrorHi() / GetVar(TotalPDFq2Bins[i],"nSig")->getVal(),2.) + pow(PsiYieldErr / PsiYield,2.)) << endl;
	      fileFitResults << "====================================================================" << endl;
	    }

	  myString.clear(); myString.str("");
	  if (CheckGoodFit(fitResult) == true)
	    {
	      myString << ID << "   ";
	      myString << (((useEffPDF == true) && (GetVar(TotalPDFq2Bins[i],"FlS")  != NULL)) ? GetVar(TotalPDFq2Bins[i],"FlS")->getVal()  : -2.0) << "   ";
	      myString << (((useEffPDF == true) && (GetVar(TotalPDFq2Bins[i],"AfbS") != NULL)) ? GetVar(TotalPDFq2Bins[i],"AfbS")->getVal() : -2.0) << "   ";
	      myString << VecHistoMeas->operator[](2)->GetBinContent(i+1) << "   ";
	      myString << NLLvalue;
	    }
	  else myString << ID << "   " << -2.0 << "   " << -2.0 << "   " << -2.0 << "   " << -2.0;
	  
	  fileFitSystematics << myString.str() << endl;
	}
    }
}


void MakeMass2AnglesToy (RooAbsPdf* TotalPDF, RooRealVar* x, RooRealVar* y, RooRealVar* z, TCanvas* Canv, unsigned int FitType, unsigned int nToy, vector<vector<unsigned int>*>* configParam, int specBin, vector<vector<string>*>* fitParam, RooArgSet* vecConstr, string fileName)
{
  unsigned int nEntryToy;
  stringstream myString;
  unsigned int it = 1;
  RooArgSet* constrToys;
  RooRealVar* tmpVar;
  RooPlot* myFrame;

  
  // #####################
  // # Initialize p.d.f. #
  // #####################
  nEntryToy = CopyFitResults(TotalPDF,specBin,fitParam);


  // #####################
  // # Apply constraints #
  // #####################
  ClearVars(vecConstr);
  BuildMassConstraints(vecConstr,TotalPDF,"All");
  BuildAngularConstraints(vecConstr,TotalPDF);
  if ((specBin != 3) && (specBin != 5))
    {
      BuildPhysicsConstraints(vecConstr,TotalPDF,"FsS");
      BuildPhysicsConstraints(vecConstr,TotalPDF,"AsS");
    }
  constrToys = (RooArgSet*)vecConstr->clone("constrToys");
  BuildPhysicsConstraints(vecConstr,TotalPDF,"AfbBound"); // # Special constraint for AFB boundaries #
  BuildPhysicsConstraints(vecConstr,TotalPDF,"AsBound");  // # Special constraint for AS  boundaries #


  // ###################################
  // # Print variables and constraints #
  // ###################################
  PrintVariables(TotalPDF->getVariables(),"vars");
  PrintVariables(vecConstr,"cons");
  PrintVariables(constrToys,"cons");


  // #############################
  // # Toy-MC generation and fit #
  // #############################
  RooMCStudy* MyToy = new RooMCStudy(*TotalPDF,RooArgList(*x,*y,*z),Extended(true),ExternalConstraints(*constrToys),FitOptions(Extended(true),ExternalConstraints(*vecConstr),Minos(USEMINOS)));
  MyToy->generateAndFit(nToy,nEntryToy,true);
  Canv->Divide(5,4);


  if ((GetVar(TotalPDF,"meanS") != NULL) && (GetVar(TotalPDF,"meanS")->getError() != 0.0) && (IsInConstraints(vecConstr,"meanS") == false))
    {
      RooRealVar* tmpVar = GetVar(TotalPDF,"meanS");
      RooPlot* myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("meanS"))->operator[](specBin).c_str()) -
 				       0.4*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](specBin).c_str())),
 				       atof(fitParam->operator[](Utility->GetFitParamIndx("meanS"))->operator[](specBin).c_str()) +
 				       0.4*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"sigmaS1") != NULL) && (GetVar(TotalPDF,"sigmaS1")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaS1") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaS1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](specBin).c_str()) -
 			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](specBin).c_str())),
 			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](specBin).c_str()) +
 			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"sigmaS2") != NULL) && (GetVar(TotalPDF,"sigmaS2")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaS2") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaS2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS2"))->operator[](specBin).c_str()) -
 			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS2"))->operator[](specBin).c_str())),
 			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS2"))->operator[](specBin).c_str()) +
 			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaS2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"fracMassS") != NULL) && (GetVar(TotalPDF,"fracMassS")->getError() != 0.0) && (IsInConstraints(vecConstr,"fracMassS") == false))
    {
      tmpVar = GetVar(TotalPDF,"fracMassS");
      myFrame = tmpVar->frame(-0.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"tau1") != NULL) && (GetVar(TotalPDF,"tau1")->getError() != 0.0) && (IsInConstraints(vecConstr,"tau1") == false))
    {
      tmpVar = GetVar(TotalPDF,"tau1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("tau1"))->operator[](specBin).c_str()) - 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("tau1"))->operator[](specBin).c_str())),
 			      atof(fitParam->operator[](Utility->GetFitParamIndx("tau1"))->operator[](specBin).c_str()) + 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("tau1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"tau2") != NULL) && (GetVar(TotalPDF,"tau2")->getError() != 0.0) && (IsInConstraints(vecConstr,"tau2") == false))
    {
      tmpVar = GetVar(TotalPDF,"tau2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("tau2"))->operator[](specBin).c_str()) - 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("tau2"))->operator[](specBin).c_str())),
 			      atof(fitParam->operator[](Utility->GetFitParamIndx("tau2"))->operator[](specBin).c_str()) + 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("tau2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"fracMassBExp") != NULL) && (GetVar(TotalPDF,"fracMassBExp")->getError() != 0.0) && (IsInConstraints(vecConstr,"fracMassBExp") == false))
    {
      tmpVar = GetVar(TotalPDF,"fracMassBExp");
      myFrame = tmpVar->frame(-0.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"meanR1") != NULL) && (GetVar(TotalPDF,"meanR1")->getError() != 0.0) && (IsInConstraints(vecConstr,"meanR1") == false))
    {
      tmpVar = GetVar(TotalPDF,"meanR1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("meanR1"))->operator[](specBin).c_str()) - 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](specBin).c_str())),
 			      atof(fitParam->operator[](Utility->GetFitParamIndx("meanR1"))->operator[](specBin).c_str()) + 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"sigmaR1") != NULL) && (GetVar(TotalPDF,"sigmaR1")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaR1") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaR1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](specBin).c_str()) -
 			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](specBin).c_str())),
 			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](specBin).c_str()) +
 			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"meanR2") != NULL) && (GetVar(TotalPDF,"meanR2")->getError() != 0.0) && (IsInConstraints(vecConstr,"meanR2") == false))
    {
      tmpVar = GetVar(TotalPDF,"meanR2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("meanR2"))->operator[](specBin).c_str()) - 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](specBin).c_str())),
 			      atof(fitParam->operator[](Utility->GetFitParamIndx("meanR2"))->operator[](specBin).c_str()) + 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    } 
   
  if ((GetVar(TotalPDF,"sigmaR2") != NULL) && (GetVar(TotalPDF,"sigmaR2")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaR2") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaR2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](specBin).c_str()) -
 			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](specBin).c_str())),
 			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](specBin).c_str()) +
 			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }
  
  if ((GetVar(TotalPDF,"fracMassBRPeak") != NULL) && (GetVar(TotalPDF,"fracMassBRPeak")->getError() != 0.0) && (IsInConstraints(vecConstr,"fracMassBRPeak") == false))
    {
      tmpVar = GetVar(TotalPDF,"fracMassBRPeak");
      myFrame = tmpVar->frame(-0.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"meanL1") != NULL) && (GetVar(TotalPDF,"meanL1")->getError() != 0.0) && (IsInConstraints(vecConstr,"meanL1") == false))
    {
      tmpVar = GetVar(TotalPDF,"meanL1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("meanL1"))->operator[](specBin).c_str()) - 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](specBin).c_str())),
 			      atof(fitParam->operator[](Utility->GetFitParamIndx("meanL1"))->operator[](specBin).c_str()) + 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"sigmaL1") != NULL) && (GetVar(TotalPDF,"sigmaL1")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaL1") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaL1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](specBin).c_str()) -
 			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](specBin).c_str())),
 			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](specBin).c_str()) +
 			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"meanL2") != NULL) && (GetVar(TotalPDF,"meanL2")->getError() != 0.0) && (IsInConstraints(vecConstr,"meanL2") == false))
    {
      tmpVar = GetVar(TotalPDF,"meanL2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("meanL2"))->operator[](specBin).c_str()) - 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](specBin).c_str())),
 			      atof(fitParam->operator[](Utility->GetFitParamIndx("meanL2"))->operator[](specBin).c_str()) + 0.2*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"sigmaL2") != NULL) && (GetVar(TotalPDF,"sigmaL2")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaL2") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaL2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](specBin).c_str()) -
 			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](specBin).c_str())),
 			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](specBin).c_str()) +
 			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"fracMassBLPeak") != NULL) && (GetVar(TotalPDF,"fracMassBLPeak")->getError() != 0.0) && (IsInConstraints(vecConstr,"fracMassBLPeak") == false))
    {
      tmpVar = GetVar(TotalPDF,"fracMassBLPeak");
      myFrame = tmpVar->frame(-0.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"fracMassBPeak") != NULL) && (GetVar(TotalPDF,"fracMassBPeak")->getError() != 0.0) && (IsInConstraints(vecConstr,"fracMassBPeak") == false))
    {
      tmpVar = GetVar(TotalPDF,"fracMassBPeak");
      myFrame = tmpVar->frame(-0.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }
  
  if ((GetVar(TotalPDF,"nBkgPeak") != NULL) && (GetVar(TotalPDF,"nBkgPeak")->getError() != 0.0) && (IsInConstraints(vecConstr,"nBkgPeak") == false))
    {
      tmpVar = GetVar(TotalPDF,"nBkgPeak");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgPeak"))->operator[](specBin).c_str()) -
 			      0.6*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgPeak"))->operator[](specBin).c_str())),
 			      atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgPeak"))->operator[](specBin).c_str()) +
 			      0.6*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgPeak"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"nSig") != NULL) && (GetVar(TotalPDF,"nSig")->getError() != 0.0) && (IsInConstraints(vecConstr,"nSig") == false))
    {
      tmpVar = GetVar(TotalPDF,"nSig");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("nSig"))->operator[](specBin).c_str()) - 0.6*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("nSig"))->operator[](specBin).c_str())),
 			      atof(fitParam->operator[](Utility->GetFitParamIndx("nSig"))->operator[](specBin).c_str()) + 0.6*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("nSig"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"nBkgComb") != NULL) && (GetVar(TotalPDF,"nBkgComb")->getError() != 0.0) && (IsInConstraints(vecConstr,"nBkgComb") == false))
    {
      tmpVar = GetVar(TotalPDF,"nBkgComb");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgComb"))->operator[](specBin).c_str()) -
 			      0.6*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgComb"))->operator[](specBin).c_str())),
 			      atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgComb"))->operator[](specBin).c_str()) +
 			      0.6*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("nBkgComb"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
    {
      myString.clear(); myString.str("");
      myString << "c1Poly" << i;
      
      if ((GetVar(TotalPDF,myString.str().c_str()) != NULL) && (GetVar(TotalPDF,myString.str().c_str())->getError() != 0.0) && (IsInConstraints(vecConstr,myString.str().c_str()) == false))
 	{
 	  tmpVar = GetVar(TotalPDF,myString.str().c_str());
 	  myFrame = tmpVar->frame(-6.0,6.0);
 	  MyToy->plotParamOn(myFrame);
 	  Canv->cd(it++);
 	  myFrame->Draw();
 	  myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
 	  Canv->cd(it++);
 	  myFrame->Draw();
 	}
    }

  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
    {
      myString.clear(); myString.str("");
      myString << "c2Poly" << i;
      
      if ((GetVar(TotalPDF,myString.str().c_str()) != NULL) && (GetVar(TotalPDF,myString.str().c_str())->getError() != 0.0) && (IsInConstraints(vecConstr,myString.str().c_str()) == false))
 	{
 	  tmpVar = GetVar(TotalPDF,myString.str().c_str());
 	  myFrame = tmpVar->frame(-6.0,6.0);
 	  MyToy->plotParamOn(myFrame);
 	  Canv->cd(it++);
 	  myFrame->Draw();
 	  myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
 	  Canv->cd(it++);
 	  myFrame->Draw();
 	}
    }

  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
    {
      myString.clear(); myString.str("");
      myString << "p1Poly" << i;
      
      if ((GetVar(TotalPDF,myString.str().c_str()) != NULL) && (GetVar(TotalPDF,myString.str().c_str())->getError() != 0.0) && (IsInConstraints(vecConstr,myString.str().c_str()) == false))
 	{
 	  tmpVar = GetVar(TotalPDF,myString.str().c_str());
 	  myFrame = tmpVar->frame(-6.0,6.0);
 	  MyToy->plotParamOn(myFrame);
 	  Canv->cd(it++);
 	  myFrame->Draw();
 	  myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
 	  Canv->cd(it++);
 	  myFrame->Draw();
 	}
    }
  
  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
    {
      myString.clear(); myString.str("");
      myString << "p2Poly" << i;
      
      if ((GetVar(TotalPDF,myString.str().c_str()) != NULL) && (GetVar(TotalPDF,myString.str().c_str())->getError() != 0.0) && (IsInConstraints(vecConstr,myString.str().c_str()) == false))
 	{
 	  tmpVar = GetVar(TotalPDF,myString.str().c_str());
 	  myFrame = tmpVar->frame(-6.0,6.0);
 	  MyToy->plotParamOn(myFrame);
 	  Canv->cd(it++);
 	  myFrame->Draw();
 	  myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
 	  Canv->cd(it++);
 	  myFrame->Draw();
 	}
    }

 if ((GetVar(TotalPDF,"FlS") != NULL) && (GetVar(TotalPDF,"FlS")->getError() != 0.0) && (IsInConstraints(vecConstr,"FlS") == false))
    {
      tmpVar = GetVar(TotalPDF,"FlS");
      myFrame = tmpVar->frame(-0.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"AfbS") != NULL) && (GetVar(TotalPDF,"AfbS")->getError() != 0.0) && (IsInConstraints(vecConstr,"AfbS") == false))
    {
      tmpVar = GetVar(TotalPDF,"AfbS");
      myFrame = tmpVar->frame(-1.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"FsS") != NULL) && (GetVar(TotalPDF,"FsS")->getError() != 0.0) && (IsInConstraints(vecConstr,"FsS") == false))
    {
      tmpVar = GetVar(TotalPDF,"FsS");
      myFrame = tmpVar->frame(-1.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"AsS") != NULL) && (GetVar(TotalPDF,"AsS")->getError() != 0.0) && (IsInConstraints(vecConstr,"AsS") == false))
    {
      tmpVar = GetVar(TotalPDF,"AsS");
      myFrame = tmpVar->frame(-1.1,1.1);
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  myFrame = MyToy->plotNLL(-2000.0,0.0);
  Canv->cd(it++);
  myFrame->Draw();
  Canv->Update();


  // #############################################
  // # Re-make all the fits and save the results #
  // #############################################
  TH1D* histoDiff1 = new TH1D("histoDiff1","histoDiff1",50,-1.0,1.0);
  histoDiff1->SetFillColor(kAzure+6);
  histoDiff1->SetXTitle("(fit #font[122]{\55} pdf)");
  histoDiff1->SetYTitle("Entries [#]");

  TH1D* histoPull1 = new TH1D("histoPull1","histoPull1",30,-5.0,5.0);
  histoPull1->SetFillColor(kAzure+6);
  histoPull1->SetXTitle("(fit #font[122]{\55} pdf) / #sigma");
  histoPull1->SetYTitle("Entries [#]");

  TH1D* histoNLL1 = new TH1D("histoNLL1","histoNLL1",30,1.0,-1.0);
  histoNLL1->SetFillColor(kGreen-7);
  histoNLL1->SetXTitle("NLL");
  histoNLL1->SetYTitle("Entries [#]");

  TH1D* histoDiff2 = new TH1D("histoDiff2","histoDiff2",50,-1.0,1.0);
  histoDiff2->SetFillColor(kAzure+6);
  histoDiff2->SetXTitle("(fit #font[122]{\55} pdf)");
  histoDiff2->SetYTitle("Entries [#]");

  TH1D* histoPull2 = new TH1D("histoPull2","histoPull2",30,-5.0,5.0);
  histoPull2->SetFillColor(kAzure+6);
  histoPull2->SetXTitle("(fit #font[122]{\55} pdf) / #sigma");
  histoPull2->SetYTitle("Entries [#]");

  TH1D* histoNLL2 = new TH1D("histoNLL2","histoNLL2",30,1.0,-1.0);
  histoNLL2->SetFillColor(kGreen-7);
  histoNLL2->SetXTitle("NLL");
  histoNLL2->SetYTitle("Entries [#]");

  TCanvas* cB0Toy = new TCanvas("cB0Toy","cB0Toy", 20, 20, 1800, 1800);

  cout << "\n@@@ Now fit total TOY invariant mass and angles @@@" << endl;
  RooDataSet* toySample;
  RooFitResult* fitResult;
  double NLLvalue;

  string varName;
  for (unsigned int i = 0; i < nToy; i++)
    {
      cout << "\n@@@ Now fitting toy #" << i << " @@@\n" << endl;

      toySample = (RooDataSet*)MyToy->genData(i);
      CopyFitResults(TotalPDF,specBin,fitParam);
      fitResult = MakeMass2AnglesFit(toySample,&TotalPDF,x,y,z,cB0Toy,FitType,configParam,specBin,vecConstr,&NLLvalue,NULL,i+1);


      // ######################################################
      // # To verify the outcome of the RooFit toy-MC studies #
      // ######################################################
      if ((CheckGoodFit(fitResult) == true) && (GetVar(TotalPDF,"FlS") != NULL) && (GetVar(TotalPDF,"AfbS") != NULL))
	{
	  varName = "FlS";
	  myString.clear(); myString.str("");
	  myString << fitParam->operator[](Utility->GetFitParamIndx(varName.c_str()))->operator[](specBin);
	  if (TotalPDF->getVariables()->getRealValue(varName.c_str()) > atof(myString.str().c_str()))
	    histoPull1->Fill((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(myString.str().c_str())) / fabs(GetVar(TotalPDF,varName.c_str())->getErrorLo()));
	  else histoPull1->Fill((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(myString.str().c_str())) / fabs(GetVar(TotalPDF,varName.c_str())->getErrorHi()));
	  histoDiff1->Fill(TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(myString.str().c_str()));
	  histoNLL1->Fill(NLLvalue);
	  
	  varName = "AfbS";
	  myString.clear(); myString.str("");
	  myString << fitParam->operator[](Utility->GetFitParamIndx(varName.c_str()))->operator[](specBin);
	  if (TotalPDF->getVariables()->getRealValue(varName.c_str()) > atof(myString.str().c_str()))
	    histoPull2->Fill((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(myString.str().c_str())) / fabs(GetVar(TotalPDF,varName.c_str())->getErrorLo()));
	  else histoPull2->Fill((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(myString.str().c_str())) / fabs(GetVar(TotalPDF,varName.c_str())->getErrorHi()));
	  histoDiff2->Fill(TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(myString.str().c_str()));
	  histoNLL2->Fill(NLLvalue);
	}
      delete fitResult;
    }


  TCanvas* cNLL1 = new TCanvas("cNLL1","cNLL1",10,10,900,500);
  cNLL1->Divide(3,1);
  cNLL1->cd(1);
  histoDiff1->Draw();
  cNLL1->cd(2);
  histoNLL1->Draw();
  cNLL1->cd(3);
  histoPull1->Draw();
  cNLL1->Update();

  TCanvas* cNLL2 = new TCanvas("cNLL2","cNLL2",10,10,900,500);
  cNLL2->Divide(3,1);
  cNLL2->cd(1);
  histoDiff2->Draw();
  cNLL2->cd(2);
  histoNLL2->Draw();
  cNLL2->cd(3);
  histoPull2->Draw();
  cNLL2->Update();

  
  // ##############
  // # Save plots #
  // ##############
  if (SAVEPLOT == true)
    {
      TFile* fNLL;

      Canv->Print(fileName.c_str());

      myString.clear(); myString.str("");
      myString << "FL_" << specBin << "_DIFF.root";
      fNLL = new TFile(fileName.replace(fileName.find(".root")-1,6,myString.str()).c_str(),"RECREATE");
      fNLL->cd();
      histoDiff1->Write();
      fNLL->Close();
      delete fNLL;

      myString.clear(); myString.str("");
      myString << "_FL_" << specBin << "_PULL.root";
      fNLL = new TFile(fileName.replace(fileName.find("_FL_"),15,myString.str()).c_str(),"RECREATE");
      fNLL->cd();
      histoPull1->Write();
      fNLL->Close();
      delete fNLL;

      myString.clear(); myString.str("");
      myString << "_FL_" << specBin << "_NLL.root";
      fNLL = new TFile(fileName.replace(fileName.find("_FL_"),15,myString.str()).c_str(),"RECREATE");
      fNLL->cd();
      histoNLL1->Write();
      fNLL->Close();
      delete fNLL;

      myString.clear(); myString.str("");
      myString << "_AFB_" << specBin << "_DIFF.root";
      fNLL = new TFile(fileName.replace(fileName.find("_FL_"),14,myString.str()).c_str(),"RECREATE");
      fNLL->cd();
      histoDiff2->Write();
      fNLL->Close();
      delete fNLL;

      myString.clear(); myString.str("");
      myString << "_AFB_" << specBin << "_PULL.root";
      fNLL = new TFile(fileName.replace(fileName.find("_AFB_"),16,myString.str()).c_str(),"RECREATE");
      fNLL->cd();
      histoPull2->Write();
      fNLL->Close();
      delete fNLL;

      myString.clear(); myString.str("");
      myString << "_AFB_" << specBin << "_NLL.root";
      fNLL = new TFile(fileName.replace(fileName.find("_AFB_"),16,myString.str()).c_str(),"RECREATE");
      fNLL->cd();
      histoNLL2->Write();
      fNLL->Close();
      delete fNLL;
    }
}


void CloseAllAndQuit (TApplication* theApp, TFile* NtplFile)
{
  fileFitResults.close();
  fileFitSystematics.close();
  fileFitSystematicsInput.close();
  
  if (TotalPDFPsi != NULL)                DeleteFit(TotalPDFPsi,"All");
  if (TotalPDFRejectPsi != NULL)          DeleteFit(TotalPDFRejectPsi,"All");
  if (SingleCandNTuple_JPsi != NULL)      delete SingleCandNTuple_JPsi;
  if (SingleCandNTuple_PsiP != NULL)      delete SingleCandNTuple_PsiP;
  if (SingleCandNTuple_RejectPsi != NULL) delete SingleCandNTuple_RejectPsi;
  if (NtplFile != NULL)                   NtplFile->Close("R");

  gROOT->CloseFiles();
  theApp->Terminate(0);
}


int main(int argc, char** argv)
{
  if (argc >= 4)
    {
      // ##################
      // # Main variables #
      // ##################
      stringstream myString;
      vector<string>* vecParStr;
      
      string fileName           = "";
      string correct4Efficiency = "";
      string tmpFileName        = "";

      int specBin               = -1;
      unsigned int fileIndx     = 0;

      bool useEffPDF            = false;

      unsigned int FitType      = atoi(argv[1]);

      TFile* NtplFile           = NULL;

      
      if (((((FitType >= 1) && (FitType <= 6))   ||
	    ((FitType >= 33) && (FitType <= 36)) ||
	    ((FitType >= 41) && (FitType <= 46)) ||
	    ((FitType >= 61) && (FitType <= 66)) ||
	    ((FitType >= 81) && (FitType <= 86))) && (argc >= 4)) ||
	  ((FitType >= 93) && (FitType <= 96) && (argc == 6)) ||
	  ((FitType >= 21) && (FitType <= 26) && (argc == 8)) ||
	  ((((FitType >= 53) && (FitType <= 56)) || ((FitType >= 73) && (FitType <= 76))) && (argc == 4)))
	{
	  ParameterFILE = PARAMETERFILEIN;


	  // #################################
	  // # Check that FitType is correct #
	  // #################################
	  if (((FitType >= 1) && (FitType <= 6))   ||
	      ((FitType >= 21) && (FitType <= 26)) ||
	      ((FitType >= 33) && (FitType <= 36)) ||
	      ((FitType >= 41) && (FitType <= 46)) ||
	      ((FitType >= 53) && (FitType <= 56)) ||
	      ((FitType >= 61) && (FitType <= 66)) ||
	      ((FitType >= 73) && (FitType <= 76)) ||
	      ((FitType >= 81) && (FitType <= 86)) ||
	      ((FitType >= 93) && (FitType <= 96)))
	    {
	      fileName           = argv[2];
	      correct4Efficiency = argv[3];
	      if (((FitType >= 41) && (FitType <= 46)) || ((FitType >= 53) && (FitType <= 56))) specBin = 3; // # 3 = q^2 bin to which the J/psi belongs
	      if (((FitType >= 61) && (FitType <= 66)) || ((FitType >= 73) && (FitType <= 76))) specBin = 5; // # 5 = q^2 bin to which the psi(2S) belongs
	    }


	  // ###################################
	  // # Check that FitOption is correct #
	  // ###################################
	  if ((correct4Efficiency != "noEffCorr") && (correct4Efficiency != "EffCorrAnalyPDF") && (correct4Efficiency != "EffCorrGenAnalyPDF"))
	    {
	      cout << "[ExtractYield::main]\tIncorrect option parameter " << correct4Efficiency << endl;
	      exit(1);
	    }


	  // #################################
	  // # Read the q^2 bin and the rest #
	  // #################################
	  if (argc >= 5) specBin = atoi(argv[4]);
	  if ((correct4Efficiency == "EffCorrAnalyPDF") || (correct4Efficiency == "EffCorrGenAnalyPDF")) useEffPDF = true;
	  if (correct4Efficiency == "EffCorrGenAnalyPDF")
	    {
	      tmpFileName = argv[5];
	      fileIndx    = atoi(argv[6]);
	    }
	  else if ((!(((FitType >= 21) && (FitType <= 26)) || ((FitType >= 93) && FitType <= 96))) &&
		   ((correct4Efficiency == "noEffCorr") || (correct4Efficiency == "EffCorrAnalyPDF")))
	    {	
	      if (argc >= 6) tmpFileName = argv[5];
	      if (argc == 7) fileIndx    = atoi(argv[6]);
	    }
	  else if ((FitType >= 21) && (FitType <= 26))
	    {
	      nToy          = atoi(argv[5]);
	      ParameterFILE = argv[6];
	      fileIndx      = atoi(argv[7]);
	    }
	  else if ((FitType >= 93) && (FitType <= 96)) fileIndx = atoi(argv[5]);


	  cout << "@@@ Input variables from command line @@@" << endl;
	  cout << "- FitType = " <<               FitType << endl;
	  cout << "- input/outputFile.root = " << fileName.c_str() << endl;
	  cout << "- correct4Efficiency = " <<    correct4Efficiency << endl;
	  cout << "- specBin = " <<               specBin << endl;
	  cout << "- tmpFileName = " <<           tmpFileName.c_str() << endl;
	  cout << "- fileIndx = " <<              fileIndx << endl;
	  cout << "- ParameterFILE = " <<         ParameterFILE.c_str() << endl;
	  cout << "- nToy = " <<                  nToy << endl;

	  cout << "\n@@@ Internal settings @@@" << endl;
	  cout << "NBINS = " <<               NBINS << endl;
	  cout << "MULTYIELD = " <<           MULTYIELD << endl;
	  cout << "NCOEFFPOLYBKG = " <<       NCOEFFPOLYBKG << endl;
	  cout << "SLEWRATECONSTR = " <<      SLEWRATECONSTR << endl;
	  cout << "NORMBIN = " <<             NORMBIN << endl;
	  cout << "POLYCOEFRANGE = " <<       POLYCOEFRANGE << endl;

	  cout << "\nApplyConstr = " <<       ApplyConstr << endl;
	  cout << "SAVEPOLY = " <<            SAVEPOLY << endl;
	  cout << "RESETOBS = " <<            RESETOBS << endl;
	  cout << "SETBATCH  = " <<           SETBATCH << endl;
	  cout << "TESTeffFUNC = " <<         TESTeffFUNC << endl;
	  cout << "SAVEPLOT = " <<            SAVEPLOT << endl;
	  cout << "MakeMuMuPlots = " <<       MakeMuMuPlots << endl;
	  cout << "MAKEGRAPHSCAN = " <<       MAKEGRAPHSCAN << endl;
	  cout << "USEMINOS = " <<            USEMINOS << endl;
	  cout << "UseSPwave = " <<           UseSPwave << endl;

	  cout << "\nPARAMETERFILEIN = " <<   PARAMETERFILEIN << endl;
	  cout << "PARAMETERFILEOUT = " <<    PARAMETERFILEOUT << endl;
	  cout << "FitSysFILEOutput = " <<    FitSysFILEOutput << endl;
	  cout << "FitSysFILEInput = " <<     FitSysFILEInput << endl;

	  legNames[0] = "Data";
	  legNames[1] = "Total fit";
	  legNames[2] = "Signal";
	  legNames[3] = "Comb. bkg";
	  legNames[4] = "Peak. bkg";


	  if (SETBATCH == true)
	    {
	      cout << "\n@@@ Setting batch mode @@@" << endl;
	      gROOT->SetBatch(true);

	      // #############################
	      // # Turn off all the printout #
	      // #############################
	      RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
	    }
	  TApplication* theApp = new TApplication("Applications", &argc, argv);


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
	  gStyle->SetPadTopMargin(0.11);
	  gStyle->SetPadBottomMargin(0.12);
	  gStyle->SetTitleOffset(1.05,"x");
	  gStyle->SetTitleOffset(1.0,"y");
	  gStyle->SetTitleSize(0.05,"x");
	  gStyle->SetTitleSize(0.05,"y");
	  gStyle->SetLabelSize(0.05,"x");
	  gStyle->SetLabelSize(0.05,"y");
	  TGaxis::SetMaxDigits(3);
	  gStyle->SetStatY(0.9);


	  // ##############################
	  // # Initialize fit output file #
	  // ##############################
	  fileName.replace(fileName.find(".root"),5,".txt");
	  fileFitResults.open(fileName.c_str(),ios_base::out);
	  if (fileFitResults.good() == false)
	    {
	      cout << "[ExtractYield::main]\tError opening file : " << fileName.c_str() << endl;
	      CloseAllAndQuit(theApp,NtplFile);
	    }
	  fileFitResults << "@@@@@@@@@@@@@@@@@ Fit results : B0 --> K*0 mu+ mu- @@@@@@@@@@@@@@@@@" << endl;
	  fileFitResults << "Fit ID: " << fileIndx << endl;
	  fileName.replace(fileName.find(".txt"),4,".root");


 	  // ###################
	  // # Read parameters #
 	  // ###################
	  Utility = new Utils();
	  Utility->ReadBins(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins);
	  Utility->ReadGenericParam(ParameterFILE);
	  Utility->ReadSelectionCuts(ParameterFILE);
	  Utility->ReadFitStartingValues(ParameterFILE,&fitParam,&configParam,Utility->ParFileBlockN("fitValGlob"));
	  myString.clear(); myString.str("");
	  myString << fitParam[Utility->GetFitParamIndx("nSig")]->operator[]((NORMBIN == 3 ? 1 : 2)).c_str(); // Read J/psi yield from global-fit results
	  SetValueAndErrors(NULL,"",1.0,&myString,&PsiYield,&PsiYerr,&PsiYerr);
	  fileFitResults << "Normalization channel yield: " << PsiYield << " +/- " << PsiYerr << endl;
	  for (unsigned int i = 0; i < fitParam.size(); i++)
	    {
	      fitParam[i]->clear();
	      delete fitParam[i];
	    }
	  fitParam.clear();
	  for (unsigned int i = 0; i < configParam.size(); i++)
	    {
	      configParam[i]->clear();
	      delete configParam[i];
	    }
	  configParam.clear();
	  if (FitType == 2) Utility->ReadFitStartingValues(ParameterFILE,&fitParam,&configParam,Utility->ParFileBlockN("fitValGlob"));
	  else              Utility->ReadFitStartingValues((((correct4Efficiency != "EffCorrGenAnalyPDF") && tmpFileName.size() != 0) ? tmpFileName : ParameterFILE),
							   &fitParam,&configParam,
							   (((correct4Efficiency != "EffCorrGenAnalyPDF") && tmpFileName.size() != 0) ? 0 : Utility->ParFileBlockN("fitValBins")));


	  // #############################################################
	  // # Make q^2 bins for histograms and fill efficiency matrices #
	  // #############################################################
	  q2BinsHisto = new double[q2Bins.size()];
	  for (unsigned int i = 0; i < q2Bins.size(); i++) q2BinsHisto[i] = q2Bins[i];


	  // ###########################################################################
	  // # Prepare file to save fit results in case of systematic error evaluation #
	  // ###########################################################################
	  myString.clear(); myString.str("");
	  if (correct4Efficiency == "EffCorrGenAnalyPDF")
	    {
	      ResetAngularParam(&fitParam,&q2Bins);
	      myString << FitSysFILEOutput << "_" << specBin << ".txt";
	    }
	  else if (specBin != -1) myString << "FitSystematics_" << specBin << ".txt";
	  else                    myString << "FitSystematics.txt";
	  fileFitSystematics.open(myString.str().c_str(),ios_base::app);
	  if (fileFitSystematics.good() == false)
	    {
	      cout << "[ExtractYield::main]\tError opening file : " << myString.str().c_str() << endl;
	      CloseAllAndQuit(theApp,NtplFile);
	    }


	  // ###########################################################################
	  // # Reset polynomial coefficients and angular observables for a fresh start #
	  // ###########################################################################
	  if (RESETOBS == true)
	    {
	      ResetCombPolyParam(&fitParam,&q2Bins);
	      ResetAngularParam(&fitParam,&q2Bins);
	    }


	  // ##############################################
	  // # Read coefficents for analytical efficiency #
	  // ##############################################
	  if (correct4Efficiency == "EffCorrGenAnalyPDF") Utility->ReadAnalyticalEff(tmpFileName.c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&effFuncs,"effFuncs",0);
	  else                                            Utility->ReadAnalyticalEff(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,&effFuncs,"effFuncs",Utility->ParFileBlockN("analyEffAvg"));


	  // ####################################################
	  // # Check if the efficiencies are positively defined #
	  // ####################################################
	  if ((TESTeffFUNC == true) && (correct4Efficiency != "EffCorrGenAnalyPDF"))
	    for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins.size()-1 : specBin+1); i++)
	      if (Utility->EffMinValue2D(&cosThetaKBins,&cosThetaLBins,effFuncs[i]) < 0.0)
		{
		  cout << "[ExtractYield::main]\tNegative efficiency function for q^2 bin #" << i << endl;
		  CloseAllAndQuit(theApp,NtplFile);
		}


 	  // ####################################################################################################
	  // # Read other parameters : this also allow also to understand if the parameter file is well written #
 	  // ####################################################################################################
	  LUMI = Utility->ReadLumi(ParameterFILE);
	  bool IsThisData = Utility->IsThisData(ParameterFILE);
	  if (IsThisData == true) cout << "\n@@@ I recognize that this is a DATA file @@@" << endl;
	  else cout << "\n@@@ I recognize that this is a Monte Carlo file @@@" << endl;


	  // ###################################################################
	  // # Fit to Afb: copy Fl value from file to the vector of parameters #
	  // ###################################################################
	  myString.clear(); myString.str("");
	  if ((correct4Efficiency == "EffCorrGenAnalyPDF") && ((FitType == 4) || (FitType == 34) || (FitType == 44) || (FitType == 54) || (FitType == 64) || (FitType == 74)))
	    {
	      double var0;
	      double var1;
	      double var2;
	      double var3;


	      myString << FitSysFILEInput << "_FL_" << specBin << ".txt";
	      fileFitSystematicsInput.open(myString.str().c_str(),ios_base::in);
	      if (fileFitSystematicsInput.good() == false)
		{
		  cout << "[ExtractYield::main]\tError opening file : " << myString.str().c_str() << endl;
		  CloseAllAndQuit(theApp,NtplFile);
		}


	      // #################################
	      // # Assign value to fit parameter #
	      // #################################
	      var0 = 0.0;
	      var1 = 0.0;
	      var2 = 0.0;
	      var3 = 0.0;
	      for (unsigned int i = 0; i <= fileIndx; i++)
		if (fileFitSystematicsInput.eof() == true)
		  {
		    cout << "\n@@@ Reached EoF: " << myString.str().c_str() << " from efficiency systematics file " << myString.str().c_str() << " - index #" << fileIndx << " @@@" << endl;
		    myString.clear(); myString.str("");
		    myString << fileIndx << "   " << -2.0 << "   " << -2.0 << "   " << -2.0 << "   " << -2.0;
		    fileFitSystematics << myString.str() << endl;
		    CloseAllAndQuit(theApp,NtplFile);
		  }
		else fileFitSystematicsInput >> var0 >> var1 >> var2 >> var3;
	      if (var1 < 0.0)
		{
		  cout << "\n@@@ Found non converged fit parameter: " << var1 << " from efficiency systematics file " << myString.str().c_str() << " - index #" << fileIndx << " @@@" << endl;
		  myString.clear(); myString.str("");
		  myString << fileIndx << "   " << -2.0 << "   " << -2.0 << "   " << -2.0 << "   " << -2.0;
		  fileFitSystematics << myString.str() << endl;
		  CloseAllAndQuit(theApp,NtplFile);
		}
	      else
		{
		  cout << "\n@@@ I successfully read the FL starting value: " << var1 << " from efficiency systematics file " << myString.str().c_str() << " - index #" << fileIndx << " @@@" << endl;

		  myString.clear(); myString.str("");
		  myString << var1;
		  fitParam[Utility->GetFitParamIndx("FlS")]->operator[](specBin) = myString.str();
		}
	      fileFitSystematicsInput.close();
	    }


	  // ###########################################################################
	  // # Fit to BF: copy Fl and Afb values from file to the vector of parameters #
	  // ###########################################################################
	  if ((correct4Efficiency == "EffCorrGenAnalyPDF") && ((FitType == 1) || (FitType == 41) || (FitType == 61)))
	    {
	      double var0;
	      double var1;
	      double var2;
	      double var3;


	      // ############################
	      // # Copy data from q^2 bin n #
	      // ############################
	      myString.clear(); myString.str("");
	      myString << FitSysFILEInput << "_FLAFB_" << specBin << ".txt";
	      fileFitSystematicsInput.open(myString.str().c_str(),ios_base::in);
	      if (fileFitSystematicsInput.good() == false)
	    	{
	    	  cout << "[ExtractYield::main]\tError opening file : " << myString.str().c_str() << endl;
	    	  CloseAllAndQuit(theApp,NtplFile);
	    	}


	      // #################################
	      // # Assign value to fit parameter #
	      // #################################
	      var0 = 0.0;
	      var1 = 0.0;
	      var2 = 0.0;
	      var3 = 0.0;
	      for (unsigned int i = 0; i <= fileIndx; i++)
	    	if (fileFitSystematicsInput.eof() == true)
	    	  {
	    	    cout << "\n@@@ Reached EoF: " << myString.str().c_str() << " from efficiency systematics file " << myString.str().c_str() << " - index #" << fileIndx << " @@@" << endl;
		    myString.clear(); myString.str("");
		    myString << fileIndx << "   " << -2.0 << "   " << -2.0 << "   " << -2.0 << "   " << -2.0;
		    fileFitSystematics << myString.str() << endl;
	    	    CloseAllAndQuit(theApp,NtplFile);
	    	  }
	    	else fileFitSystematicsInput >> var0 >> var1 >> var2 >> var3;
	      if ((var1 < 0.0) || (var2 < -1.0))
	    	{
	    	  cout << "\n@@@ Found non converged fit parameter: " << var1 << " or " << var2 << " from efficiency systematics file " << myString.str().c_str() << " - index #" << fileIndx << " @@@" << endl;
		  myString.clear(); myString.str("");
		  myString << fileIndx << "   " << -2.0 << "   " << -2.0 << "   " << -2.0 << "   " << -2.0;
		  fileFitSystematics << myString.str() << endl;
	    	  CloseAllAndQuit(theApp,NtplFile);
	    	}
	      else
	    	{
	    	  cout << "\n@@@ I successfully read FL and AFB starting values: " << var1 << " and " << var2 << " from efficiency systematics file " << myString.str().c_str() << " - index #" << fileIndx << " @@@" << endl;

	    	  myString.clear(); myString.str("");
	    	  myString << var1;
	    	  fitParam[Utility->GetFitParamIndx("FlS")]->operator[](specBin) = myString.str();

	    	  myString.clear(); myString.str("");
	    	  myString << var2;
	    	  fitParam[Utility->GetFitParamIndx("AfbS")]->operator[](specBin) = myString.str();
	    	}
	      fileFitSystematicsInput.close();


	      // ########################################
	      // # Copy data from normalization q^2 bin #
	      // ########################################
	      myString.clear(); myString.str("");
	      myString << FitSysFILEInput << "_FLAFB_" << NORMBIN << ".txt";
	      fileFitSystematicsInput.open(myString.str().c_str(),ios_base::in);
	      if (fileFitSystematicsInput.good() == false)
	    	{
	    	  cout << "[ExtractYield::main]\tError opening file : " << myString.str().c_str() << endl;
	    	  CloseAllAndQuit(theApp,NtplFile);
	    	}


	      // #################################
	      // # Assign value to fit parameter #
	      // #################################
	      var0 = 0.0;
	      var1 = 0.0;
	      var2 = 0.0;
	      var3 = 0.0;
	      for (unsigned int i = 0; i <= fileIndx; i++)
	    	if (fileFitSystematicsInput.eof() == true)
	    	  {
	    	    cout << "\n@@@ Reached EoF: " << myString.str().c_str() << " from efficiency systematics file " << myString.str().c_str() << " - index #" << fileIndx << " @@@" << endl;
		    myString.clear(); myString.str("");
		    myString << fileIndx << "   " << -2.0 << "   " << -2.0 << "   " << -2.0 << "   " << -2.0;
		    fileFitSystematics << myString.str() << endl;
	    	    CloseAllAndQuit(theApp,NtplFile);
	    	  }
	    	else fileFitSystematicsInput >> var0 >> var1 >> var2 >> var3;
	      if ((var1 < 0.0) || (var2 < -1.0))
	    	{
	    	  cout << "\n@@@ Found non converged fit parameter: " << var1 << " or " << var2 << " from efficiency systematics file " << myString.str().c_str() << " - index #" << fileIndx << " @@@" << endl;
		  myString.clear(); myString.str("");
		  myString << fileIndx << "   " << -2.0 << "   " << -2.0 << "   " << -2.0 << "   " << -2.0;
		  fileFitSystematics << myString.str() << endl;
	    	  CloseAllAndQuit(theApp,NtplFile);
	    	}
	      else
	    	{
	    	  cout << "\n@@@ I successfully read FL and AFB starting values: " << var1 << " and " << var2 << " from efficiency systematics file " << myString.str().c_str() << " - index #" << fileIndx << " @@@" << endl;

	    	  myString.clear(); myString.str("");
	    	  myString << var1;
	    	  fitParam[Utility->GetFitParamIndx("FlS")]->operator[](NORMBIN) = myString.str();

	    	  myString.clear(); myString.str("");
	    	  myString << var2;
	    	  fitParam[Utility->GetFitParamIndx("AfbS")]->operator[](NORMBIN) = myString.str();
	    	}
	      fileFitSystematicsInput.close();
	    }


	  // ###################
	  // # Select fit type #
	  // ###################
	  if (((FitType >= 1) && (FitType <= 6))   ||
	      ((FitType >= 33) && (FitType <= 36)) ||
	      ((FitType >= 41) && (FitType <= 46)) ||
	      ((FitType >= 53) && (FitType <= 56)) ||
	      ((FitType >= 61) && (FitType <= 66)) ||
	      ((FitType >= 73) && (FitType <= 76)))
	    {
	      NtplFile = new TFile(fileName.c_str(), "READ");
	      theTree  = (TTree*) NtplFile->Get("B0SingleCand/B0KstMuMuSingleCandNTuple");
	      NTuple   = new B0KstMuMuSingleCandTreeContent();
	      NTuple->Init();
	    

	      // #################
	      // # Make datasets #
	      // #################
	      cout << "\n@@@ Making datasets @@@" << endl;
	      MakeDataSets(NTuple,FitType);

	      cout << "\n@@@ NTuple with J/psi region @@@" << endl;
	      SingleCandNTuple_JPsi->Print("v");

	      cout << "\n@@@ NTuple with psi(2S) region @@@" << endl;
	      SingleCandNTuple_PsiP->Print("v");

	      cout << "\n@@@ NTuple without J/psi and psi(2S) regions @@@" << endl;
	      SingleCandNTuple_RejectPsi->Print("v");


	      // ##############################
	      // # Select the proper fit type #
	      // ##############################
	      if ((FitType == 1) || (FitType == 2) || (FitType == 41) || (FitType == 61))
		{
		  if (FitType == 2)
		    {
		      double NLLvalue;


		      // ##########################
		      // # 1D-fit to B0 inv. mass #
		      // ##########################
		      cout << "\n@@@ Now fit total invariant mass NON resonant channel @@@" << endl;
		      TCanvas* cB0MassArbRejectPsi = new TCanvas("cB0MassArbRejectPsi","cB0MassArbRejectPsi",10, 10, 700, 500);
		      InstantiateMassFit(&TotalPDFRejectPsi,B0MassArb,"TotalPDFRejectPsi",&configParam,0);
		      CopyFitResults(TotalPDFRejectPsi,0,&fitParam);
		      ClearVars(&vecConstr);
		      BuildMassConstraints(&vecConstr,TotalPDFRejectPsi,"sign");
		      PrintVariables(TotalPDFRejectPsi->getVariables(),"vars");
		      PrintVariables(&vecConstr,"cons");
		      MakeMassFit(SingleCandNTuple_RejectPsi,&TotalPDFRejectPsi,B0MassArb,cB0MassArbRejectPsi,FitType,&configParam,0,&vecConstr,&NLLvalue,NULL,fileIndx);

		      // ##############################################
		      // # Save fit results back into prarameter file #
		      // ##############################################
		      vecParStr = SaveFitResults(TotalPDFRejectPsi,0,&fitParam,&configParam,&vecConstr);
		      Utility->SaveFitValues(PARAMETERFILEOUT,vecParStr,-1,"  signal   ");
		      delete vecParStr;


		      // ####################################
		      // # Save the workspace for debugging #
		      // ####################################
		      // @TMP@
		      // RooWorkspace w("w");
		      // w.import(*TotalPDFRejectPsi);
		      // w.importClassCode();
		      // w.import(*SingleCandNTuple_RejectPsi);
		      // w.writeToFile("Dinardo.root");


		      cout << "\n@@@ Now fit total invariant mass resonant J/psi channel @@@" << endl;
		      TCanvas* cB0MassArbJPsi = new TCanvas("cB0MassArbJPsi","cB0MassArbJPsi",10, 10, 700, 500);
		      InstantiateMassFit(&TotalPDFPsi,B0MassArb,"TotalPDFPsi",&configParam,1);
		      CopyFitResults(TotalPDFPsi,1,&fitParam);
		      ClearVars(&vecConstr);
		      BuildMassConstraints(&vecConstr,TotalPDFPsi,"sign");
		      PrintVariables(TotalPDFPsi->getVariables(),"vars");
		      PrintVariables(&vecConstr,"cons");
		      MakeMassFit(SingleCandNTuple_JPsi,&TotalPDFPsi,B0MassArb,cB0MassArbJPsi,FitType,&configParam,1,&vecConstr,&NLLvalue,NULL,fileIndx);

		      // ##############################################
		      // # Save fit results back into prarameter file #
		      // ##############################################
		      vecParStr = SaveFitResults(TotalPDFPsi,1,&fitParam,&configParam,&vecConstr);
		      Utility->SaveFitValues(PARAMETERFILEOUT,vecParStr,-1,"   J/psi   ");
		      delete vecParStr;


		      cout << "\n@@@ Now fit total invariant mass resonant psi(2S) channel @@@" << endl;
		      TCanvas* cB0MassArbPsiP = new TCanvas("cB0MassArbPsiP","cB0MassArbPsiP",10, 10, 700, 500);
		      InstantiateMassFit(&TotalPDFPsi,B0MassArb,"TotalPDFPsi",&configParam,2);
		      CopyFitResults(TotalPDFPsi,2,&fitParam);
		      ClearVars(&vecConstr);
		      BuildMassConstraints(&vecConstr,TotalPDFPsi,"sign");
		      PrintVariables(TotalPDFPsi->getVariables(),"vars");
		      PrintVariables(&vecConstr,"cons");
		      MakeMassFit(SingleCandNTuple_PsiP,&TotalPDFPsi,B0MassArb,cB0MassArbPsiP,FitType,&configParam,2,&vecConstr,&NLLvalue,NULL,fileIndx);

		      // ##############################################
		      // # Save fit results back into prarameter file #
		      // ##############################################
		      vecParStr = SaveFitResults(TotalPDFPsi,2,&fitParam,&configParam,&vecConstr);
		      Utility->SaveFitValues(PARAMETERFILEOUT,vecParStr,-1,"  psi(2S)  ");
		      delete vecParStr;
		    }
		  else
		    {
		      // ######################################
		      // # 1D-fit to B0 inv. mass per q^2 bin #
		      // ######################################
		      vector<TH1D*> VecHistoMeas;
		      VecHistoMeas.push_back(new TH1D("histoMeas0", "B#lower[0.4]{^{0}} \\rightarrow #font[122]{K}#lower[0.4]{^{*0}}(#font[122]{K}#lower[0.4]{^{#font[122]{+}}}\\pi#lower[0.4]{^{#font[122]{\55}}}) \\mu#lower[0.4]{^{#font[122]{+}}} \\mu#lower[0.4]{^{#font[122]{\55}}} / B#lower[0.4]{^{0}} \\rightarrow #font[122]{K}#lower[0.4]{^{*0}}(#font[122]{K}#lower[0.4]{^{#font[122]{+}}}\\pi#lower[0.4]{^{#font[122]{\55}}}) J/\\psi(\\mu#lower[0.4]{^{#font[122]{+}}}\\mu#lower[0.4]{^{#font[122]{\55}}}) Branching Fraction",q2Bins.size()-1, q2BinsHisto));
		      VecHistoMeas[0]->SetStats(false);
		      VecHistoMeas[0]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		      VecHistoMeas[0]->SetYTitle("dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}7}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");
		      VecHistoMeas[0]->GetYaxis()->SetRangeUser(0.0,1.2);

		      cout << "\n@@@ Now fit invariant mass per mumu q^2 bins @@@" << endl;
		      if (FitType == 1) IterativeMassFitq2Bins(SingleCandNTuple_RejectPsi,
							       useEffPDF,
							       PsiYield,PsiYerr,
							       B0MassArb,
							       CosThetaMuArb,
							       CosThetaKArb,
							       specBin,
							       FitType,
							       &VecHistoMeas,
							       &q2Bins,
							       &configParam,&fitParam,
							       &effFuncs,
							       &vecConstr,
							       fileIndx);
		      else if (FitType == 41) IterativeMassFitq2Bins(SingleCandNTuple_JPsi,
								     useEffPDF,
								     PsiYield,PsiYerr,
								     B0MassArb,
								     CosThetaMuArb,
								     CosThetaKArb,
								     specBin,
								     FitType,
								     &VecHistoMeas,
								     &q2Bins,
								     &configParam,&fitParam,
								     &effFuncs,
								     &vecConstr,
								     fileIndx);
		      else IterativeMassFitq2Bins(SingleCandNTuple_PsiP,
						  useEffPDF,
						  PsiYield,PsiYerr,
						  B0MassArb,
						  CosThetaMuArb,
						  CosThetaKArb,
						  specBin,
						  FitType,
						  &VecHistoMeas,
						  &q2Bins,
						  &configParam,&fitParam,
						  &effFuncs,
						  &vecConstr,
						  fileIndx);
		      
		      TCanvas* cHistoMeas = new TCanvas("cHistoMeas","cHistoMeas",10, 10, 900, 600);
		      cHistoMeas->Divide(2,2);
		      cHistoMeas->cd(1);
		      VecHistoMeas[0]->SetMarkerStyle(20);
		      VecHistoMeas[0]->Draw("pe1");
		      DrawString(LUMI);
		      
		      cHistoMeas->Update();
		    }
		}
	      else if ((FitType == 3) || (FitType == 33) || (FitType == 43) || (FitType == 53) || (FitType == 63) || (FitType == 73)) // Fl-Afb-fit
		{
		  // #########################
		  // # 2D-fit Fl per q^2 bin #
		  // #########################
		  vector<TH1D*> VecHistoMeas;
		  VecHistoMeas.push_back(new TH1D("histoMeas0", "B#lower[0.4]{^{0}} \\rightarrow #font[122]{K}#lower[0.4]{^{*0}}(#font[122]{K}#lower[0.4]{^{#font[122]{+}}}\\pi#lower[0.4]{^{#font[122]{\55}}}) \\mu#lower[0.4]{^{#font[122]{+}}} \\mu#lower[0.4]{^{#font[122]{\55}}} : F_{L} vs dimuon q#lower[0.4]{^{2}}", q2Bins.size()-1, q2BinsHisto));
		  VecHistoMeas.push_back(new TH1D("histoMeas1", "B#lower[0.4]{^{0}} \\rightarrow #font[122]{K}#lower[0.4]{^{*0}}(#font[122]{K}#lower[0.4]{^{#font[122]{+}}}\\pi#lower[0.4]{^{#font[122]{\55}}}) \\mu#lower[0.4]{^{#font[122]{+}}} \\mu#lower[0.4]{^{#font[122]{\55}}} : A_{FB} vs dimuon q#lower[0.4]{^{2}}", q2Bins.size()-1, q2BinsHisto));
		  VecHistoMeas.push_back(new TH1D("histoMeas2", "histoMeas2", q2Bins.size()-1, q2BinsHisto));
		  VecHistoMeas.push_back(new TH1D("histoMeas3", "B#lower[0.4]{^{0}} \\rightarrow #font[122]{K}#lower[0.4]{^{*0}}(#font[122]{K}#lower[0.4]{^{#font[122]{+}}}\\pi#lower[0.4]{^{#font[122]{\55}}}) \\mu#lower[0.4]{^{#font[122]{+}}} \\mu#lower[0.4]{^{#font[122]{\55}}} : dBF/dq#lower[0.4]{^{2}} vs dimuon q#lower[0.4]{^{2}}", q2Bins.size()-1, q2BinsHisto));
		  VecHistoMeas[0]->SetStats(false);
		  VecHistoMeas[1]->SetStats(false);
		  VecHistoMeas[2]->SetStats(false);
		  VecHistoMeas[3]->SetStats(false);
		  VecHistoMeas[0]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		  VecHistoMeas[1]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		  VecHistoMeas[2]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		  VecHistoMeas[3]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		  VecHistoMeas[0]->SetYTitle("F_{L}");
		  VecHistoMeas[1]->SetYTitle("A_{FB}");
		  VecHistoMeas[2]->SetYTitle("");
		  VecHistoMeas[3]->SetYTitle("dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}7}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");
		  VecHistoMeas[0]->GetYaxis()->SetRangeUser(-0.01,1.01);
		  VecHistoMeas[1]->GetYaxis()->SetRangeUser(-1.0 - 0.01,1.0 + 0.01);
		  VecHistoMeas[2]->GetYaxis()->SetRangeUser(0.0,1.0);
		  VecHistoMeas[3]->GetYaxis()->SetRangeUser(0.0,1.2);

		  cout << "\n@@@ Now fit invariant mass and cos(theta_K) per mumu q^2 bins @@@" << endl;
		  if ((FitType == 3) || (FitType == 33)) IterativeMassAngleFitq2Bins(SingleCandNTuple_RejectPsi,
										     useEffPDF,
										     PsiYield,PsiYerr,
										     B0MassArb,
										     CosThetaMuArb,
										     CosThetaKArb,
										     specBin,
										     FitType,
										     &VecHistoMeas,
										     &q2Bins,
										     &configParam,&fitParam,
										     &effFuncs,
										     &vecConstr,
										     fileIndx);
		  else if ((FitType == 43) || (FitType == 53)) IterativeMassAngleFitq2Bins(SingleCandNTuple_JPsi,
											   useEffPDF,
											   PsiYield,PsiYerr,
											   B0MassArb,
											   CosThetaMuArb,
											   CosThetaKArb,
											   specBin,
											   FitType,
											   &VecHistoMeas,
											   &q2Bins,
											   &configParam,&fitParam,
											   &effFuncs,
											   &vecConstr,
											   fileIndx);
		  else IterativeMassAngleFitq2Bins(SingleCandNTuple_PsiP,
						   useEffPDF,
						   PsiYield,PsiYerr,
						   B0MassArb,
						   CosThetaMuArb,
						   CosThetaKArb,
						   specBin,
						   FitType,
						   &VecHistoMeas,
						   &q2Bins,
						   &configParam,&fitParam,
						   &effFuncs,
						   &vecConstr,
						   fileIndx);
		  
		  TCanvas* cHistoMeas = new TCanvas("cHistoMeas","cHistoMeas",10, 10, 900, 600);
		  cHistoMeas->Divide(2,2);
		  cHistoMeas->cd(1);
		  VecHistoMeas[0]->SetMarkerStyle(20);
		  VecHistoMeas[0]->Draw("pe1");
		  DrawString(LUMI);
		  cHistoMeas->cd(2);
		  VecHistoMeas[1]->SetMarkerStyle(20);
		  VecHistoMeas[1]->Draw("pe1");
		  DrawString(LUMI);
		  cHistoMeas->cd(3);
		  VecHistoMeas[3]->SetMarkerStyle(20);
		  VecHistoMeas[3]->Draw("pe1");
		  DrawString(LUMI);
		  
		  cHistoMeas->Update();
		}
	      else if ((FitType == 4) || (FitType == 34) || (FitType == 44) || (FitType == 54) || (FitType == 64) || (FitType == 74)) // Afb-Fl-fit
		{
		  // #############################
		  // # 2D-fit Afb-Fl per q^2 bin #
		  // #############################
		  vector<TH1D*> VecHistoMeas;
		  VecHistoMeas.push_back(new TH1D("histoMeas0", "B#lower[0.4]{^{0}} \\rightarrow #font[122]{K}#lower[0.4]{^{*0}}(#font[122]{K}#lower[0.4]{^{#font[122]{+}}}\\pi#lower[0.4]{^{#font[122]{\55}}}) \\mu#lower[0.4]{^{#font[122]{+}}} \\mu#lower[0.4]{^{#font[122]{\55}}} : F_{L} vs dimuon q#lower[0.4]{^{2}}", q2Bins.size()-1, q2BinsHisto));
		  VecHistoMeas.push_back(new TH1D("histoMeas1", "B#lower[0.4]{^{0}} \\rightarrow #font[122]{K}#lower[0.4]{^{*0}}(#font[122]{K}#lower[0.4]{^{#font[122]{+}}}\\pi#lower[0.4]{^{#font[122]{\55}}}) \\mu#lower[0.4]{^{#font[122]{+}}} \\mu#lower[0.4]{^{#font[122]{\55}}} : A_{FB} vs dimuon q#lower[0.4]{^{2}}", q2Bins.size()-1, q2BinsHisto));
		  VecHistoMeas.push_back(new TH1D("histoMeas2", "histoMeas2", q2Bins.size()-1, q2BinsHisto));
		  VecHistoMeas.push_back(new TH1D("histoMeas3", "B#lower[0.4]{^{0}} \\rightarrow #font[122]{K}#lower[0.4]{^{*0}}(#font[122]{K}#lower[0.4]{^{#font[122]{+}}}\\pi#lower[0.4]{^{#font[122]{\55}}}) \\mu#lower[0.4]{^{#font[122]{+}}} \\mu#lower[0.4]{^{#font[122]{\55}}} : dBF/dq#lower[0.4]{^{2}} vs dimuon q#lower[0.4]{^{2}}", q2Bins.size()-1, q2BinsHisto));
		  VecHistoMeas[0]->SetStats(false);
		  VecHistoMeas[1]->SetStats(false);
		  VecHistoMeas[2]->SetStats(false);
		  VecHistoMeas[3]->SetStats(false);
		  VecHistoMeas[0]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		  VecHistoMeas[1]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		  VecHistoMeas[2]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		  VecHistoMeas[3]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		  VecHistoMeas[0]->SetYTitle("F_{L}");
		  VecHistoMeas[1]->SetYTitle("A_{FB}");
		  VecHistoMeas[2]->SetYTitle("");
		  VecHistoMeas[3]->SetYTitle("dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}7}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");
		  VecHistoMeas[0]->GetYaxis()->SetRangeUser(-0.01,1.01);
		  VecHistoMeas[1]->GetYaxis()->SetRangeUser(-1.0 - 0.01,1.0 + 0.01);
		  VecHistoMeas[2]->GetYaxis()->SetRangeUser(0.0,1.0);
		  VecHistoMeas[3]->GetYaxis()->SetRangeUser(0.0,1.2);
		  
		  cout << "\n@@@ Now fit invariant mass and cos(theta_l) per mumu q^2 bins @@@" << endl;
		  if ((FitType == 4) || (FitType == 34)) IterativeMassAngleFitq2Bins(SingleCandNTuple_RejectPsi,
										     useEffPDF,
										     PsiYield,PsiYerr,
										     B0MassArb,
										     CosThetaMuArb,
										     CosThetaKArb,
										     specBin,
										     FitType,
										     &VecHistoMeas,
										     &q2Bins,
										     &configParam,&fitParam,
										     &effFuncs,
										     &vecConstr,
										     fileIndx);
		  else if ((FitType == 44) || (FitType == 54)) IterativeMassAngleFitq2Bins(SingleCandNTuple_JPsi,
											   useEffPDF,
											   PsiYield,PsiYerr,
											   B0MassArb,
											   CosThetaMuArb,
											   CosThetaKArb,
											   specBin,
											   FitType,
											   &VecHistoMeas,
											   &q2Bins,
											   &configParam,&fitParam,
											   &effFuncs,
											   &vecConstr,
											   fileIndx);
		  else IterativeMassAngleFitq2Bins(SingleCandNTuple_PsiP,
						   useEffPDF,
						   PsiYield,PsiYerr,
						   B0MassArb,
						   CosThetaMuArb,
						   CosThetaKArb,
						   specBin,
						   FitType,
						   &VecHistoMeas,
						   &q2Bins,
						   &configParam,&fitParam,
						   &effFuncs,
						   &vecConstr,
						   fileIndx);
		  
		  TCanvas* cHistoMeas = new TCanvas("cHistoMeas","cHistoMeas",10, 10, 900, 600);
		  cHistoMeas->Divide(2,2);
		  cHistoMeas->cd(1);
		  VecHistoMeas[0]->SetMarkerStyle(20);
		  VecHistoMeas[0]->Draw("pe1");
		  DrawString(LUMI);
		  cHistoMeas->cd(2);
		  VecHistoMeas[1]->SetMarkerStyle(20);
		  VecHistoMeas[1]->Draw("pe1");
		  DrawString(LUMI);
		  cHistoMeas->cd(3);
		  VecHistoMeas[3]->SetMarkerStyle(20);
		  VecHistoMeas[3]->Draw("pe1");
		  DrawString(LUMI);
		  
		  cHistoMeas->Update();
		}
	      else if ((FitType == 5) || (FitType == 35) || (FitType == 45) || (FitType == 55) || (FitType == 65) || (FitType == 75)) // Fl-At2-Atim-fit
		{
		  // ##################################
		  // # 2D-fit Fl-At2-Atim per q^2 bin #
		  // ##################################
		  vector<TH1D*> VecHistoMeas;
		  VecHistoMeas.push_back(new TH1D("histoMeas0", "B#lower[0.4]{^{0}} \\rightarrow #font[122]{K}#lower[0.4]{^{*0}}(#font[122]{K}#lower[0.4]{^{#font[122]{+}}}\\pi#lower[0.4]{^{#font[122]{\55}}}) \\mu#lower[0.4]{^{#font[122]{+}}} \\mu#lower[0.4]{^{#font[122]{\55}}} : F_{L} vs dimuon q#lower[0.4]{^{2}}", q2Bins.size()-1, q2BinsHisto));
		  VecHistoMeas.push_back(new TH1D("histoMeas1", "B#lower[0.4]{^{0}} \\rightarrow #font[122]{K}#lower[0.4]{^{*0}}(#font[122]{K}#lower[0.4]{^{#font[122]{+}}}\\pi#lower[0.4]{^{#font[122]{\55}}}) \\mu#lower[0.4]{^{#font[122]{+}}} \\mu#lower[0.4]{^{#font[122]{\55}}} : A_{T}^{(2)} vs dimuon q#lower[0.4]{^{2}}", q2Bins.size()-1, q2BinsHisto));
		  VecHistoMeas.push_back(new TH1D("histoMeas1", "B#lower[0.4]{^{0}} \\rightarrow #font[122]{K}#lower[0.4]{^{*0}}(#font[122]{K}#lower[0.4]{^{#font[122]{+}}}\\pi#lower[0.4]{^{#font[122]{\55}}}) \\mu#lower[0.4]{^{#font[122]{+}}} \\mu#lower[0.4]{^{#font[122]{\55}}} : A_{T}^{(im)} vs dimuon q#lower[0.4]{^{2}}", q2Bins.size()-1, q2BinsHisto));
		  VecHistoMeas.push_back(new TH1D("histoMeas3", "B#lower[0.4]{^{0}} \\rightarrow #font[122]{K}#lower[0.4]{^{*0}}(#font[122]{K}#lower[0.4]{^{#font[122]{+}}}\\pi#lower[0.4]{^{#font[122]{\55}}}) \\mu#lower[0.4]{^{#font[122]{+}}} \\mu#lower[0.4]{^{#font[122]{\55}}} : dBF/dq#lower[0.4]{^{2}} vs dimuon q#lower[0.4]{^{2}}", q2Bins.size()-1, q2BinsHisto));
		  VecHistoMeas[0]->SetStats(false);
		  VecHistoMeas[1]->SetStats(false);
		  VecHistoMeas[2]->SetStats(false);
		  VecHistoMeas[3]->SetStats(false);
		  VecHistoMeas[0]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		  VecHistoMeas[1]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		  VecHistoMeas[2]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		  VecHistoMeas[3]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		  VecHistoMeas[0]->SetYTitle("F_{L}");
		  VecHistoMeas[1]->SetYTitle("A_{FB}");
		  VecHistoMeas[2]->SetYTitle("A_{T}^{(im)}");
		  VecHistoMeas[3]->SetYTitle("dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}7}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");
		  VecHistoMeas[0]->GetYaxis()->SetRangeUser(-0.01,1.01);
		  VecHistoMeas[1]->GetYaxis()->SetRangeUser(-1.01,1.01);
		  VecHistoMeas[2]->GetYaxis()->SetRangeUser(-1.01,1.01);
		  VecHistoMeas[2]->GetYaxis()->SetRangeUser(0.0,1.2);
		  
		  cout << "\n@@@ Now fit invariant mass and phi per mumu q^2 bins @@@" << endl;
		  if ((FitType == 5) || (FitType == 35)) IterativeMassAngleFitq2Bins(SingleCandNTuple_RejectPsi,
										     useEffPDF,
										     PsiYield,PsiYerr,
										     B0MassArb,
										     CosThetaMuArb,
										     PhiKstMuMuPlaneArb,
										     specBin,
										     FitType,
										     &VecHistoMeas,
										     &q2Bins,
										     &configParam,&fitParam,
										     &effFuncs,
										     &vecConstr,
										     fileIndx);
		  else if ((FitType == 45) || (FitType == 55)) IterativeMassAngleFitq2Bins(SingleCandNTuple_JPsi,
											   useEffPDF,
											   PsiYield,PsiYerr,
											   B0MassArb,
											   CosThetaMuArb,
											   PhiKstMuMuPlaneArb,
											   specBin,
											   FitType,
											   &VecHistoMeas,
											   &q2Bins,
											   &configParam,&fitParam,
											   &effFuncs,
											   &vecConstr,
											   fileIndx);
		  else IterativeMassAngleFitq2Bins(SingleCandNTuple_PsiP,
						   useEffPDF,
						   PsiYield,PsiYerr,
						   B0MassArb,
						   CosThetaMuArb,
						   PhiKstMuMuPlaneArb,
						   specBin,
						   FitType,
						   &VecHistoMeas,
						   &q2Bins,
						   &configParam,&fitParam,
						   &effFuncs,
						   &vecConstr,
						   fileIndx);

		  TCanvas* cHistoMeas = new TCanvas("cHistoMeas","cHistoMeas",10, 10, 900, 600);
		  cHistoMeas->Divide(2,2);
		  cHistoMeas->cd(1);
		  VecHistoMeas[0]->SetMarkerStyle(20);
		  VecHistoMeas[0]->Draw("pe1");
		  DrawString(LUMI);
		  cHistoMeas->cd(2);
		  VecHistoMeas[1]->SetMarkerStyle(20);
		  VecHistoMeas[1]->Draw("pe1");
		  DrawString(LUMI);
		  cHistoMeas->cd(3);
		  VecHistoMeas[2]->SetMarkerStyle(20);
		  VecHistoMeas[2]->Draw("pe1");
		  DrawString(LUMI);
		  cHistoMeas->cd(4);
		  VecHistoMeas[3]->SetMarkerStyle(20);
		  VecHistoMeas[3]->Draw("pe1");
		  DrawString(LUMI);
		  
		  cHistoMeas->Update();
		}
	      else if ((FitType == 6) || (FitType == 36) || (FitType == 46) || (FitType == 56) || (FitType == 66) || (FitType == 76)) // Fl-Afb-fit AND Afb-Fl-fit
		{
		  // #############################
		  // # 3D-fit Afb-Fl per q^2 bin #
		  // #############################
		  vector<TH1D*> VecHistoMeas;

		  VecHistoMeas.push_back(new TH1D("histoMeas0", "B#lower[0.4]{^{0}} \\rightarrow #font[122]{K}#lower[0.4]{^{*0}}(#font[122]{K}#lower[0.4]{^{#font[122]{+}}}\\pi#lower[0.4]{^{#font[122]{\55}}}) \\mu#lower[0.4]{^{#font[122]{+}}} \\mu#lower[0.4]{^{#font[122]{\55}}} : F_{L} vs dimuon q#lower[0.4]{^{2}}", q2Bins.size()-1, q2BinsHisto));
		  VecHistoMeas.push_back(new TH1D("histoMeas1", "B#lower[0.4]{^{0}} \\rightarrow #font[122]{K}#lower[0.4]{^{*0}}(#font[122]{K}#lower[0.4]{^{#font[122]{+}}}\\pi#lower[0.4]{^{#font[122]{\55}}}) \\mu#lower[0.4]{^{#font[122]{+}}} \\mu#lower[0.4]{^{#font[122]{\55}}} : A_{FB} vs dimuon q#lower[0.4]{^{2}}", q2Bins.size()-1, q2BinsHisto));
		  VecHistoMeas.push_back(new TH1D("histoMeas2", "B#lower[0.4]{^{0}} \\rightarrow #font[122]{K}#lower[0.4]{^{*0}}(#font[122]{K}#lower[0.4]{^{#font[122]{+}}}\\pi#lower[0.4]{^{#font[122]{\55}}}) \\mu#lower[0.4]{^{#font[122]{+}}} \\mu#lower[0.4]{^{#font[122]{\55}}} : dBF/dq#lower[0.4]{^{2}} vs dimuon q#lower[0.4]{^{2}}", q2Bins.size()-1, q2BinsHisto));
		  VecHistoMeas[0]->SetStats(false);
		  VecHistoMeas[1]->SetStats(false);
		  VecHistoMeas[2]->SetStats(false);
		  VecHistoMeas[0]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		  VecHistoMeas[1]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		  VecHistoMeas[2]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		  VecHistoMeas[0]->SetYTitle("F_{L}");
		  VecHistoMeas[1]->SetYTitle("A_{FB}");
		  VecHistoMeas[2]->SetYTitle("dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}7}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");
		  VecHistoMeas[0]->GetYaxis()->SetRangeUser(-0.01,1.01);
		  VecHistoMeas[1]->GetYaxis()->SetRangeUser(-1.0 - 0.01,1.0 + 0.01);
		  VecHistoMeas[2]->GetYaxis()->SetRangeUser(0.0,1.2);

		  cout << "\n@@@ Now fit invariant mass, cos(theta_K) and cos(theta_l) per mumu q^2 bins @@@" << endl;
		  if ((FitType == 6) || (FitType == 36)) IterativeMass2AnglesFitq2Bins(SingleCandNTuple_RejectPsi,
		  								       useEffPDF,
		  								       PsiYield,PsiYerr,
		  								       B0MassArb,
		  								       CosThetaMuArb,
		  								       CosThetaKArb,
		  								       specBin,
		  								       FitType,
		  								       &VecHistoMeas,
		  								       &q2Bins,
		  								       &configParam,&fitParam,
		  								       &effFuncs,
		  								       &vecConstr,
		  								       fileIndx);
		  else if ((FitType == 46) || (FitType == 56)) IterativeMass2AnglesFitq2Bins(SingleCandNTuple_JPsi,
		  									     useEffPDF,
		  									     PsiYield,PsiYerr,
		  									     B0MassArb,
		  									     CosThetaMuArb,
		  									     CosThetaKArb,
		  									     specBin,
		  									     FitType,
		  									     &VecHistoMeas,
		  									     &q2Bins,
		  									     &configParam,&fitParam,
		  									     &effFuncs,
		  									     &vecConstr,
		  									     fileIndx);
		  else IterativeMass2AnglesFitq2Bins(SingleCandNTuple_PsiP,
		  				     useEffPDF,
		  				     PsiYield,PsiYerr,
		  				     B0MassArb,
		  				     CosThetaMuArb,
		  				     CosThetaKArb,
		  				     specBin,
		  				     FitType,
		  				     &VecHistoMeas,
		  				     &q2Bins,
		  				     &configParam,&fitParam,
		  				     &effFuncs,
		  				     &vecConstr,
		  				     fileIndx);
		  
		  TCanvas* cHistoMeas = new TCanvas("cHistoMeas","cHistoMeas",10, 10, 900, 600);
		  cHistoMeas->Divide(2,2);
		  cHistoMeas->cd(1);
		  VecHistoMeas[0]->SetMarkerStyle(20);
		  VecHistoMeas[0]->Draw("pe1");
		  DrawString(LUMI);
		  cHistoMeas->cd(2);
		  VecHistoMeas[1]->SetMarkerStyle(20);
		  VecHistoMeas[1]->Draw("pe1");
		  DrawString(LUMI);
		  cHistoMeas->cd(3);
		  VecHistoMeas[2]->SetMarkerStyle(20);
		  VecHistoMeas[2]->Draw("pe1");
		  DrawString(LUMI);

		  cHistoMeas->Update();
		}
	    }
	  else if ((FitType >= 21) && (FitType <= 26))
	    {
	      // #################
	      // # Make datasets #
	      // #################
	      cout << "\n@@@ Making datasets @@@" << endl;
	      MakeDataSets(NTuple,FitType);


	      // ##########################################
	      // # Set seed for random number generator   #
	      // # in case the toy-MC studies are splited #
	      // ##########################################
	      RooRandom::randomGenerator()->SetSeed(fileIndx*(q2Bins.size()-1) + specBin + 1);
	      cout << "\n@@@ Random seed for toy-MC set to: " << RooRandom::randomGenerator()->GetSeed() << " @@@" << endl;


	      if (FitType == 21)
		{
		  // ####################################
		  // # Make toy-MC for B0 inv. mass fit #
		  // ####################################
		  cout << "\n@@@ Now make TOY-MC for fit to B0 total invariant mass @@@" << endl;
		  TCanvas* cToyMC = new TCanvas("cToyMC","cToyMC",10, 10, 1600, 900);

		  InstantiateMassFit(&TotalPDFRejectPsi,B0MassArb,"TotalPDFRejectPsi",&configParam,specBin);
		  MakeMassToy(TotalPDFRejectPsi,B0MassArb,cToyMC,FitType,nToy,&configParam,specBin,&fitParam,&vecConstr,fileName);
		}
	      else if (FitType == 23) // Fl-Afb-fit
		{
		  // ##############################################
		  // # Make toy-MC for B0 inv. mass fit and angle #
		  // ##############################################
		  cout << "\n@@@ Now make TOY-MC for fit to B0 total invariant mass and cos(theta_K) @@@" << endl;
		  TCanvas* cToyMC = new TCanvas("cToyMC","cToyMC",10, 10, 1600, 900);

		  InstantiateMassAngleFit(&TotalPDFRejectPsi,useEffPDF,B0MassArb,CosThetaKArb,"TotalPDFRejectPsi",FitType,&configParam,&fitParam,specBin,effFuncs[specBin],true);
		  MakeMassAngleToy(TotalPDFRejectPsi,B0MassArb,CosThetaKArb,cToyMC,FitType,nToy,&configParam,specBin,&fitParam,&vecConstr,fileName);
		}
	      else if (FitType == 24) // Afb-Fl-fit
		{
		  // ##############################################
		  // # Make toy-MC for B0 inv. mass fit and angle #
		  // ##############################################
		  cout << "\n@@@ Now make TOY-MC for fit to B0 total invariant mass and cos(theta_l) @@@" << endl;
		  TCanvas* cToyMC = new TCanvas("cToyMC","cToyMC",10, 10, 1600, 900);

		  InstantiateMassAngleFit(&TotalPDFRejectPsi,useEffPDF,B0MassArb,CosThetaMuArb,"TotalPDFRejectPsi",FitType,&configParam,&fitParam,specBin,effFuncs[specBin],true);
		  MakeMassAngleToy(TotalPDFRejectPsi,B0MassArb,CosThetaMuArb,cToyMC,FitType,nToy,&configParam,specBin,&fitParam,&vecConstr,fileName);
		}
	      else if (FitType == 25) // Fl-At2-Atim-fit
		{
		  // ##############################################
		  // # Make toy-MC for B0 inv. mass fit and angle #
		  // ##############################################
		  cout << "\n@@@ Now make TOY-MC for fit to B0 total invariant mass and phi @@@" << endl;
		  TCanvas* cToyMC = new TCanvas("cToyMC","cToyMC",10, 10, 1600, 900);

		  InstantiateMassAngleFit(&TotalPDFRejectPsi,useEffPDF,B0MassArb,PhiKstMuMuPlaneArb,"TotalPDFRejectPsi",FitType,&configParam,&fitParam,specBin,effFuncs[specBin],true);
		  MakeMassAngleToy(TotalPDFRejectPsi,B0MassArb,PhiKstMuMuPlaneArb,cToyMC,FitType,nToy,&configParam,specBin,&fitParam,&vecConstr,fileName);
		}
	      else if (FitType == 26) // Fl-Afb-fit AND Afb-Fl-fit
		{
		  // ###############################################
		  // # Make toy-MC for B0 inv. mass fit and angles #
		  // ###############################################
		  cout << "\n@@@ Now make TOY-MC for fit to B0 total invariant mass and cos(theta_K) and cos(theta_l) @@@" << endl;
		  TCanvas* cToyMC = new TCanvas("cToyMC","cToyMC",10, 10, 1600, 900);

		  InstantiateMass2AnglesFit(&TotalPDFRejectPsi,useEffPDF,B0MassArb,CosThetaMuArb,CosThetaKArb,"TotalPDFRejectPsi",FitType,&configParam,&fitParam,specBin,effFuncs[specBin],true);
		  MakeMass2AnglesToy(TotalPDFRejectPsi,B0MassArb,CosThetaMuArb,CosThetaKArb,cToyMC,FitType,nToy,&configParam,specBin,&fitParam,&vecConstr,fileName);
		}
	    }
	  else if ((FitType >= 81) && (FitType <= 86))
	    {
	      // #################
	      // # Make datasets #
	      // #################
	      cout << "\n@@@ Making datasets @@@" << endl;
	      MakeDataSets(NTuple,FitType);


	      if (FitType == 81)
		{
		  // ###########################################
		  // # Generate dataset for B0 inv. mass model #
		  // ###########################################
		  InstantiateMassFit(&TotalPDFRejectPsi,B0MassArb,"TotalPDFRejectPsi",&configParam,specBin);
		  GenerateDataset(TotalPDFRejectPsi,RooArgSet(*B0MassArb),&q2Bins,specBin,&fitParam,fileName);
		}
	      else if (FitType == 83) // Fl-Afb-fit
		{
		  // #####################################################
		  // # Generate dataset for B0 inv. mass and angle model #
		  // #####################################################
		  InstantiateMassAngleFit(&TotalPDFRejectPsi,useEffPDF,B0MassArb,CosThetaKArb,"TotalPDFRejectPsi",FitType,&configParam,&fitParam,specBin,effFuncs[specBin],true);
		  GenerateDataset(TotalPDFRejectPsi,RooArgSet(*B0MassArb,*CosThetaKArb),&q2Bins,specBin,&fitParam,fileName);
		}
	      else if (FitType == 84) // Afb-Fl-fit
		{
		  // #####################################################
		  // # Generate dataset for B0 inv. mass and angle model #
		  // #####################################################
		  InstantiateMassAngleFit(&TotalPDFRejectPsi,useEffPDF,B0MassArb,CosThetaMuArb,"TotalPDFRejectPsi",FitType,&configParam,&fitParam,specBin,effFuncs[specBin],true);
		  GenerateDataset(TotalPDFRejectPsi,RooArgSet(*B0MassArb,*CosThetaMuArb),&q2Bins,specBin,&fitParam,fileName);
		}
	      else if (FitType == 85) // Fl-At2-Atim-fit
		{
		  // #####################################################
		  // # Generate dataset for B0 inv. mass and angle model #
		  // #####################################################
		  InstantiateMassAngleFit(&TotalPDFRejectPsi,useEffPDF,B0MassArb,PhiKstMuMuPlaneArb,"TotalPDFRejectPsi",FitType,&configParam,&fitParam,specBin,effFuncs[specBin],true);
		  GenerateDataset(TotalPDFRejectPsi,RooArgSet(*B0MassArb,*CosThetaMuArb),&q2Bins,specBin,&fitParam,fileName);
		}
	      else if (FitType == 86) // Fl-Afb-fit AND Afb-Fl-fit
		{
		  // ######################################################
		  // # Generate dataset for B0 inv. mass and angles model #
		  // ######################################################
		  InstantiateMass2AnglesFit(&TotalPDFRejectPsi,useEffPDF,B0MassArb,CosThetaMuArb,CosThetaKArb,"TotalPDFRejectPsi",FitType,&configParam,&fitParam,specBin,effFuncs[specBin],true);
		  GenerateDataset(TotalPDFRejectPsi,RooArgSet(*B0MassArb,*CosThetaMuArb,*CosThetaKArb),&q2Bins,specBin,&fitParam,fileName);
		}
	    }
	  else if ((FitType >= 93) && (FitType <= 96))
	    {
	      // #################
	      // # Make datasets #
	      // #################
	      cout << "\n@@@ Making datasets @@@" << endl;
	      MakeDataSets(NTuple,FitType);


	      if (FitType == 93) // Fl-Afb-fit
		{
		  // ################################################################
		  // # Generate new parameter file for B0 inv. mass and angle model #
		  // ################################################################
		  InstantiateMassAngleFit(&TotalPDFRejectPsi,useEffPDF,B0MassArb,CosThetaKArb,"TotalPDFRejectPsi",FitType,&configParam,&fitParam,specBin,effFuncs[specBin],true);
		}
	      else if (FitType == 94) // Afb-Fl-fit
		{
		  // ################################################################
		  // # Generate new parameter file for B0 inv. mass and angle model #
		  // ################################################################
		  InstantiateMassAngleFit(&TotalPDFRejectPsi,useEffPDF,B0MassArb,CosThetaMuArb,"TotalPDFRejectPsi",FitType,&configParam,&fitParam,specBin,effFuncs[specBin],true);
		}
	      else if (FitType == 95) // Fl-At2-Atim-fit
		{
		  // ################################################################
		  // # Generate new parameter file for B0 inv. mass and angle model #
		  // ################################################################
		  InstantiateMassAngleFit(&TotalPDFRejectPsi,useEffPDF,B0MassArb,PhiKstMuMuPlaneArb,"TotalPDFRejectPsi",FitType,&configParam,&fitParam,specBin,effFuncs[specBin],true);
		}
	      else if (FitType == 96) // Fl-Afb-fit AND Afb-Fl-fit
		{
		  // #################################################################
		  // # Generate new parameter file for B0 inv. mass and angles model #
		  // #################################################################
		  InstantiateMass2AnglesFit(&TotalPDFRejectPsi,useEffPDF,B0MassArb,CosThetaMuArb,CosThetaKArb,"TotalPDFRejectPsi",FitType,&configParam,&fitParam,specBin,effFuncs[specBin],true);
		}
	      GenerateParameterFile(TotalPDFRejectPsi,&fitParam,&configParam,&vecConstr,fileName,fileIndx,&q2Bins,specBin);
	    }


	  fileFitResults.close();
	  fileFitSystematics.close();
	  if ((SETBATCH == true) || (correct4Efficiency == "EffCorrGenAnalyPDF") || ((FitType >= 81) && (FitType <= 86)) || ((FitType >= 93) && (FitType <= 96)))
	    {
	      cout << "Bye bye !" << endl;
	      CloseAllAndQuit(theApp,NtplFile);
	    }
	  else
	    {
	      system("say \"Let's rock and roll !\"");
 	      theApp->Run (); // Eventloop on air
	    }
	}
      else
	{
	  cout << "Wrong parameter: " << endl;
	  cout << "./ExtractYield [FitType] [input/output[if toy-MC]File.root] [noEffCorr EffCorrAnalyPDF EffCorrGenAnalyPDF]" << endl;
	  cout << "               [q^2 bin to fit (0 - ...)]" << endl;
	  cout << "               [[if EffCorrGenAnalyPDF]effFileName.txt AND indx]" << endl;
	  cout << "               [[if toy-MC]nToy AND ParameterFile.txt AND indx] " << endl;
	  cout << "               [[if 93-96]indx]" << endl;
	  cout << "               [ParameterFile.txt] [indx]" << endl;

	  cout << "\n --> noEffCorr          = no eff. correction" << endl;
	  cout << " --> EffCorrAnalyPDF    = analytical eff. correction used in the p.d.f." << endl;
	  cout << " --> EffCorrGenAnalyPDF = compute systematic error realted to analytical eff. correction used in the p.d.f." << endl;

	  return 1;
	}
    }
  else
    {
      cout << "Parameter missing: " << endl;
      cout << "./ExtractYield [FitType] [input/output[if toy-MC]File.root] [noEffCorr EffCorrAnalyPDF EffCorrGenAnalyPDF]" << endl;
      cout << "               [q^2 bin to fit (0 - ...)]" << endl;
      cout << "               [[if EffCorrGenAnalyPDF]effFileName.txt AND indx]" << endl;
      cout << "               [[if toy-MC]nToy AND ParameterFile.txt AND indx] " << endl;
      cout << "               [[if 93-96]indx]" << endl;
      cout << "               [ParameterFile.txt] [indx]" << endl;

      cout << "\n --> noEffCorr          = no eff. correction" << endl;
      cout << " --> EffCorrAnalyPDF    = analytical eff. correction used in the p.d.f." << endl;
      cout << " --> EffCorrGenAnalyPDF = compute systematic error realted to analytical eff. correction used in the p.d.f." << endl;

      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  Signa  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 1: 1D branching fraction per q^2 bin" << endl;
      cout << "FitType = 2: 1D B0 inv. mass peak" << endl;
      cout << "FitType = 3: 2D Fl (B0Mass, cos(theta_K)) per q^2 bin" << endl;
      cout << "FitType = 4: 2D Afb-Fl (B0Mass, cos(theta_l)) per q^2 bin" << endl;
      cout << "FitType = 5: 2D Fl-At2-Atim (B0Mass, phi) per q^2 bin" << endl;
      cout << "FitType = 6: 3D Afb-Fl (B0Mass, cos(theta_K), cos(theta_l)) per q^2 bin" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 21: toy-MC 1D B0 inv. mass peak per q^2 bin" << endl;
      cout << "FitType = 23: toy-MC 2D Fl (B0Mass, cos(theta_K)) per q^2 bin" << endl;
      cout << "FitType = 24: toy-MC 2D Afb-Fl (B0Mass, cos(theta_l)) per q^2 bin" << endl;
      cout << "FitType = 25: toy-MC 2D Fl-At2-Atim (B0Mass, phi) per q^2 bin" << endl;
      cout << "FitType = 26: toy-MC 3D Afb-Fl (B0Mass, cos(theta_K), cos(theta_l)) per q^2 bin" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 33: 1D Fl cos(theta_K) per q^2 bin on GEN variables" << endl;
      cout << "FitType = 34: 1D Afb-Fl cos(theta_l) per q^2 bin on GEN variables" << endl;
      cout << "FitType = 35: 1D Fl-At2-Atim phi per q^2 bin on GEN variables" << endl;
      cout << "FitType = 36: 2D Afb-Fl (cos(theta_K), cos(theta_l)) per q^2 bin on GEN variables" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  J/psi  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 41: 1D B0 inv. mass peak for B0 --> J/psi(mumu) K*0 in appropriate q^2 bin" << endl;
      cout << "FitType = 43: 2D Fl (B0Mass, cos(theta_K)) for B0 --> J/psi(mumu) K*0" << endl;
      cout << "FitType = 44: 2D Afb-Fl (B0Mass, cos(theta_l)) for B0 --> J/psi(mumu) K*0" << endl;
      cout << "FitType = 45: 2D Fl-At2-Atim (B0Mass, phi) for B0 --> J/psi(mumu) K*0" << endl;
      cout << "FitType = 46: 3D Afb-Fl (B0Mass, cos(theta_K), cos(theta_l)) for B0 --> J/psi(mumu) K*0" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 53: 1D Fl cos(theta_K) for B0 --> J/psi(mumu) K*0 on GEN variables" << endl;
      cout << "FitType = 54: 1D Afb-Fl cos(theta_l) for B0 --> J/psi(mumu) K*0 on GEN variables" << endl;
      cout << "FitType = 55: 1D Fl-At2-Atim phi for B0 --> J/psi(mumu) K*0 on GEN variables" << endl;
      cout << "FitType = 56: 2D Afb-Fl (cos(theta_K), cos(theta_l)) for B0 --> J/psi(mumu) K*0 on GEN variables" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Psi(2S) @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 61: 1D B0 inv. mass peak for B0 --> psi(2S)(mumu) K*0 in appropriate q^2 bin" << endl;
      cout << "FitType = 63: 2D Fl (B0Mass, cos(theta_K)) for B0 --> psi(2S)(mumu) K*0" << endl;
      cout << "FitType = 64: 2D Afb-Fl (B0Mass, cos(theta_l)) for B0 --> psi(2S)(mumu) K*0" << endl;
      cout << "FitType = 65: 2D Fl-At2-Atim (B0Mass, phi) for B0 --> psi(2S)(mumu) K*0" << endl;
      cout << "FitType = 66: 3D Afb-Fl (B0Mass, cos(theta_K), cos(theta_l)) for B0 --> psi(2S)(mumu) K*0" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 73: 1D Fl cos(theta_K) for B0 --> psi(2S)(mumu) K*0 on GEN variables" << endl;
      cout << "FitType = 74: 1D Afb-Fl cos(theta_l) for B0 --> psi(2S)(mumu) K*0 on GEN variables" << endl;
      cout << "FitType = 75: 1D Fl-At2-Atim phi for B0 --> psi(2S)(mumu) K*0 on GEN variables" << endl;
      cout << "FitType = 76: 2D Afb-Fl (cos(theta_K), cos(theta_l)) for B0 --> psi(2S)(mumu) K*0 on GEN variables" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Genera dataset @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 81: 1D generare dataset with B0 inv. mass distribution from pdf" << endl;
      cout << "FitType = 83: 2D generare dataset with (B0Mass, cos(theta_K)) distribution from pdf" << endl;
      cout << "FitType = 84: 2D generare dataset with (B0Mass, cos(theta_l)) distribution from pdf" << endl;
      cout << "FitType = 85: 2D generare dataset with (B0Mass, phi) distribution from pdf" << endl;
      cout << "FitType = 86: 3D generare dataset with (B0Mass, cos(theta_K), cos(theta_l)) distribution from pdf" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      
      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Genera parameter @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 93: 2D generate paramter with (B0Mass, cos(theta_K)) distribution from pdf" << endl;
      cout << "FitType = 94: 2D generate paramter with (B0Mass, cos(theta_l)) distribution from pdf" << endl;
      cout << "FitType = 95: 2D generate paramter with (B0Mass, phi) distribution from pdf" << endl;
      cout << "FitType = 96: 3D generate paramter with (B0Mass, cos(theta_K), cos(theta_l)) distribution from pdf" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      
      return 1;
    }

  return 0;
}

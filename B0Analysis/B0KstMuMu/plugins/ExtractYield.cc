// ################################################################################
// # Program to perform the full angular analysis of the decay B0 --> K*0 mu+ mu- #
// # Author: Mauro Dinardo                                                        #
// ################################################################################
// # Search for @TMP@ to look for temporary code options                          #
// ################################################################################
// # TO DO:                                                                       #
// # - when computing systematic error from MC-statistics for efficiency I do not #
// #   propagate from control channels                                            #
// # - implement 4D fits                                                          #
// ################################################################################

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
#include <TLatex.h>
#include <Math/Functor.h>

#include <RooRealVar.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooDataSet.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooArgSet.h>
#include <RooFitResult.h>
#include <RooPolynomial.h>
#include <RooMCStudy.h>
#include <RooMinuit.h>
#include <RooWorkspace.h>
#include <RooConstVar.h>
#include <RooRandom.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooFunctorBinding.h>
#include <RooStats/RooStatsUtils.h>

#include <ctime>
#include <iostream>
#include <utility>
#include <sstream>

#include "ReadParameters.h"
#include "Utils.h"
#include "B0KstMuMuSingleCandTreeContent.h"

using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::vector;
using std::ios_base;
using std::pair;
using std::make_pair;
using namespace RooFit;


// ####################
// # Global constants #
// ####################
#define NBINS         20
#define MULTYIELD     1.0 // Multiplication factor to the number of entry in toy-MC
#define NCOEFFPOLYBKG 5   // Maximum number of coefficients (= degree) of the polynomial describing the background in the angular variables
#define DEGREEINTERP  1   // Polynomial degree for efficiency histogram interpolation
#define MAXTRIALS     0   // Maximum number of trials in case of fit failure [0 = default single trial]

#define nJPSIS 230000.0
#define nJPSIB   2500.0
#define nPSIPS  15000.0
#define nPSIPB   1500.0

// ##########################################
// # Internal flags to control the workflow #
// ##########################################
#define MakeMuMuPlots false
#define USEMINOS      false
#define SETBATCH      false
#define SAVEPOLY      false // ["true" = save bkg polynomial coefficients in new parameter file; "false" = save original values]
#define SAVEPLOT      false
#define RESETsigANG   false // Reset signal angular parameters before starting the fit
#define RESETcomANG   false // Reset combinatorial background angular parameters before starting the fit
#define FUNCERRBAND   false // Show the p.d.f. error band

// ##################
// # External files #
// ##################
#define PARAMETERFILEIN  "../python/ParameterFile.txt"
#define PARAMETERFILEOUT "ParameterFile.txt"


// ############################################
// # Global variables from configuration file #
// ############################################
double PsiYieldGoodTag, PsiYieldGoodTagErr;
double PsiYieldMisTag,  PsiYieldMisTagErr;
double LUMI;

string CTRLfitWRKflow;
string ParameterFILE;

vector<double> q2Bins;
vector<double> cosThetaKBins;
vector<double> cosThetaLBins;
vector<double> phiBins;

vector<vector<string>*>             fitParam;    // Vector containing the pointers to the vectors containing the starting values for the fit
vector<vector<unsigned int>*>       configParam; // Vector containing the pointers to the vectors containing the configuration parameters for the fit
pair< vector<TF2*>*,vector<TF2*>* > effFuncs;    // Vector containing the analytical descriptions of the efficiency per q^2 bin for good and mis-tagged events
                                                 // The first vector is for good-tagged events, while the second vecor is for mis-tagged events

// ####################
// # Global variables #
// ####################
TTree* theTree;
Utils* Utility;
B0KstMuMuSingleCandTreeContent* NTuple;

ofstream fileFitResults;
ofstream fileFitSystematics;

double* q2BinsHisto;


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
RooRealVar* rightFlavorTag;


// #################################
// # Variables and pdf for the fit #
// #################################

// ##################
// # Signal B0 mass #
// ##################
RooRealVar* meanS;

RooRealVar* sigmaS1;
RooAbsPdf*  MassS1;

RooRealVar* sigmaS2;
RooAbsPdf*  MassS2;

RooRealVar* fracMassS;
RooAbsPdf*  MassSignal;

// #################
// # Signal angles #
// #################
RooRealVar* FlS;
RooRealVar* AfbS;
RooRealVar* FsS;
RooRealVar* AsS;
RooAbsPdf*  AngleS;

// ####################
// # Total signal pdf #
// ####################
RooAbsPdf* Signal;

// ####################################
// # Combinatorial background B0 mass #
// ####################################
RooRealVar* var1;
RooRealVar* var2;
RooAbsPdf*  BkgMassExp1;
RooAbsPdf*  BkgMassExp2;

RooRealVar* fracMassBExp;
RooAbsPdf*  BkgMassComb;

// #########################
// # Mistag signal B0 mass #
// #########################
RooRealVar* sigmaMisTag1;
RooAbsPdf*  MassMisTag1;
RooRealVar* sigmaMisTag2;
RooAbsPdf*  MassMisTag2;
RooRealVar* fracMisTag;
RooAbsPdf*  MassMisTag;
RooAbsPdf*  AngleMisTag;

// ##############################
// # Peaking background B0 mass #
// ##############################
RooRealVar* meanR1;
RooRealVar* sigmaR1;
RooRealVar* meanR2;
RooRealVar* sigmaR2;
RooAbsPdf*  BkgMassRPeak1;
RooAbsPdf*  BkgMassRPeak2;

RooRealVar* fracMassBRPeak;
RooAbsPdf*  BkgMassRPeak;

RooRealVar* meanL1;
RooRealVar* sigmaL1;
RooRealVar* meanL2;
RooRealVar* sigmaL2;
RooAbsPdf*  BkgMassLPeak1;
RooAbsPdf*  BkgMassLPeak2;

RooRealVar* fracMassBLPeak;
RooAbsPdf*  BkgMassLPeak;

RooRealVar* fracMassBPeak;
RooAbsPdf*  BkgMassPeak;

// #####################
// # Background angles #
// #####################
RooRealVar* p1Poly[NCOEFFPOLYBKG];
RooAbsPdf*  BkgAngleP1;
RooRealVar* c1Poly[NCOEFFPOLYBKG];
RooAbsPdf*  BkgAngleC1;

RooRealVar* p2Poly[NCOEFFPOLYBKG];
RooAbsPdf*  BkgAngleP2;
RooRealVar* c2Poly[NCOEFFPOLYBKG];
RooAbsPdf*  BkgAngleC2;

RooRealVar* p3Poly[NCOEFFPOLYBKG];
RooAbsPdf*  BkgAngleP3;
RooRealVar* c3Poly[NCOEFFPOLYBKG];
RooAbsPdf*  BkgAngleC3;

RooAbsPdf* BkgAnglesC;
RooAbsPdf* BkgAnglesP;

// ########################
// # Total background pdf #
// ########################
RooAbsPdf* BkgMassAngleComb;
RooAbsPdf* MassAngleMisTag;
RooAbsPdf* BkgMassAnglePeak;

RooRealVar* nSig;
RooRealVar* nBkgComb;
RooRealVar* nMisTagFrac;
RooRealVar* nBkgPeak;

// #############################################
// # Total pdf for B0 --> K*0 J/psi OR psi(2S) #
// #############################################
RooAbsPdf* TotalPDFPsi;

// ##################################
// # Total pdf for B0 --> K*0 mu mu #
// ##################################
RooAbsPdf* TotalPDFRejectPsi;

// ##############################################
// # Vector containing the Gaussian constraints #
// ##############################################
RooArgSet vecConstr;


// #######################
// # Function Definition #
// #######################

// ################################################
// # Structure to make the product of two p.d.f.s #
// ################################################
struct MyProdPdf
{
public:
  MyProdPdf (RooAbsPdf& pdf1, RooAbsPdf& pdf2) : _pdf1(pdf1), _pdf2(pdf2)
  { 
    const RooArgSet* allvar1 = pdf1.getVariables();
    const RooArgSet* allvar2 = pdf2.getVariables();

    _vars.add(*allvar1);
    _vars.add(*allvar2,false);

    delete allvar1;
    delete allvar2;
  }

  int ndim ()
  {
    return _vars.getSize();
  }

  const RooArgList& vars() const
  {
    return _vars;
  }
  
  double operator() (const double* v)
  {
    for (int i = 0; i < ndim(); ++i) ((RooRealVar&)_vars[i]).setVal(v[i]);
    
    return _pdf1.getVal() * _pdf2.getVal();
  }


private:
  RooAbsPdf& _pdf1; 
  RooAbsPdf& _pdf2; 
  RooArgList _vars; 
}; 

bool CheckGoodFit              (RooFitResult* fitResult, TPaveText* paveText = NULL);
RooRealVar* GetVar             (RooAbsPdf* pdf, string varName);
void SetValueAndErrors         (RooAbsPdf* pdf, string varName, double multi, stringstream* myString, double* val, double* errLo, double* errHi);
void PrintVariables            (RooArgSet* setVar, string type);
void ClearVars                 (RooArgSet* vecConstr);
void CloseAllAndQuit           (TApplication* theApp, TFile* NtplFile);

string MakeName                (RooDataSet* data, unsigned int ID);
void DrawString                (double Lumi, RooPlot* myFrame = NULL);

bool IsInConstraints           (RooArgSet* vecConstr, string varName);
void AddGaussConstraint        (RooArgSet* vecConstr, RooAbsPdf* TotalPDF, string varName);
void AddPhysicsConstraint      (RooArgSet* vecConstr, RooAbsPdf* TotalPDF, string varName);
void BuildMassConstraints      (RooArgSet* vecConstr, RooAbsPdf* TotalPDF, string varName);
void BuildAngularConstraints   (RooArgSet* vecConstr, RooAbsPdf* TotalPDF, string varName);
void BuildPhysicsConstraints   (RooArgSet* vecConstr, RooAbsPdf* TotalPDF, string varName);

RooAbsPdf* MakeAngWithEffPDF   (TF2* effFunc, RooRealVar* y, RooRealVar* z, unsigned int FitType, bool useEffPDF, RooArgSet* VarsAng, RooArgSet* VarsPoly, vector<double>* q2Bins, int q2BinIndx);
void DeleteFit                 (RooAbsPdf* TotalPDF, string DeleteType);
void ResetCombPolyParam        (vector<vector<string>*>* fitParam);
void ResetAngularParam         (vector<vector<string>*>* fitParam);
double StoreFitResultsInFile   (RooAbsPdf** TotalPDF, RooFitResult* fitResult, RooDataSet* dataSet, RooArgSet* vecConstr);
void StorePolyResultsInFile    (RooAbsPdf** TotalPDF);
vector<string>* SaveFitResults (RooAbsPdf* TotalPDF, unsigned int q2BinIndx, vector<vector<string>*>* fitParam, vector<vector<unsigned int>*>* configParam, RooArgSet* vecConstr);
unsigned int CopyFitResults    (RooAbsPdf* TotalPDF, unsigned int q2BinIndx, vector<vector<string>*>* fitParam, unsigned int countMisTag = 0, unsigned int countGoodTag = 0);

void GenerateFitParameters     (RooAbsPdf* TotalPDF, vector<vector<string>*>* fitParam, unsigned int fileIndx, vector<double>* q2Bins, unsigned int q2BinIndx, string option);
void GenerateDataset           (RooAbsPdf* TotalPDF, RooArgSet setVar, vector<double>* q2Bins, int q2BinIndx, vector<vector<string>*>* fitParam, string fileName);
string GeneratePolynomial      (RooRealVar* var, unsigned int nCoef, string sCoef);

void FitDimuonInvMass          (RooDataSet* dataSet, RooAbsPdf** TotalPDFJPsi, RooAbsPdf** TotalPDFPsiP, RooRealVar* x, TCanvas* Canv, bool justPlotMuMuMass, bool justKeepPsi, string plotName);

void MakeDataSets              (B0KstMuMuSingleCandTreeContent* NTuple, unsigned int FitType);
unsigned int GetSignalType     (unsigned int FitType, vector<double>* q2Bins, int q2BinIndx);

// ==================
// ===> 1D MODEL <===
// ==================
void InstantiateMassFit     (RooAbsPdf** TotalPDF, RooRealVar* x, string fitName, vector<vector<unsigned int>*>* configParam, unsigned int q2BinIndx);
RooFitResult* MakeMassFit   (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* x, TCanvas* Canv, RooArgSet* vecConstr, double* NLLvalue, TPaveText* extText, unsigned int ID);
void IterativeMassFitq2Bins (RooDataSet* dataSet,
			     bool useEffPDF,
			     double PsiYieldGoodTag, double PsiYieldGoodTagErr,
			     double PsiYieldMisTag,
			     RooRealVar* x,
			     int specBin,
			     unsigned int FitType,
			     vector<TH1D*>* VecHistoMeas,
			     vector<double>* q2Bins,
			     vector<vector<unsigned int>*>* configParam, vector<vector<string>*>* fitParam,
			     RooArgSet* vecConstr,
			     unsigned int ID = 0);
void MakeMassToy            (RooAbsPdf* TotalPDF, RooRealVar* x, TCanvas* Canv, unsigned int nToy, int specBin, vector<vector<string>*>* fitParam, RooArgSet* vecConstr, string fileName);

// ==================
// ===> 3D MODEL <===
// ==================
void InstantiateMass2AnglesFit    (RooAbsPdf** TotalPDF,
				   bool useEffPDF,
				   RooRealVar* x, RooRealVar* y, RooRealVar* z,
				   string fitName, unsigned int FitType,
				   vector<vector<unsigned int>*>* configParam,
				   vector<vector<string>*>* fitParam,
				   vector<double>* q2Bins,
				   unsigned int actualq2BinIndx,
				   unsigned int q2BinIndx,
				   pair<TF2*,TF2*> effFunc);
RooFitResult* MakeMass2AnglesFit   (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* x, RooRealVar* y, RooRealVar* z, TCanvas* Canv, unsigned int FitType, RooArgSet* vecConstr, double* NLLvalue, TPaveText* extText, unsigned int ID);
void IterativeMass2AnglesFitq2Bins (RooDataSet* dataSet,
				    bool useEffPDF,
				    double PsiYieldGoodTag, double PsiYieldGoodTagErr,
				    double PsiYieldMisTag,
				    RooRealVar* x, RooRealVar* y, RooRealVar* z,
				    int specBin,
				    unsigned int FitType,
				    vector<TH1D*>* VecHistoMeas,
				    vector<double>* q2Bins,
				    vector<vector<unsigned int>*>* configParam, vector<vector<string>*>* fitParam,
				    pair< vector<TF2*>*,vector<TF2*>* > effFuncs,
				    RooArgSet* vecConstr,
				    unsigned int ID = 0);
void MakeMass2AnglesToy            (RooAbsPdf* TotalPDF, RooRealVar* x, RooRealVar* y, RooRealVar* z, TCanvas* Canv, unsigned int FitType, unsigned int nToy, int specBin, vector<vector<string>*>* fitParam, RooArgSet* vecConstr, string fileName, vector<double>* q2Bins);


// ###########################
// # Function Implementation #
// ###########################
bool CheckGoodFit (RooFitResult* fitResult, TPaveText* paveText)
// ####################################################
// # Covariance matrix quality:                       #
// # -1 : "Unknown, matrix was externally provided"   #
// #  0 : "Not calculated at all"                     #
// #  1 : "Approximation only, not accurate"          #
// #  2 : "Full matrix, but forced positive-definite" #
// #  3 : "Full, accurate covariance matrix"          #
// ####################################################
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
// # If the error is an empty string --> setAsymError = -1/+1 and setError = 1 #
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
	  GetVar(pdf,varName)->setError((fabs(*errLo) + fabs(*errHi)) / 2.);
      	}
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
      exit (EXIT_FAILURE);
    }
}


void ClearVars (RooArgSet* vecConstr)
{
  int nEle = vecConstr->getSize();

  TIterator* it = vecConstr->createIterator();
  for (int i = 0; i < nEle; i++) delete it->Next();

  vecConstr->removeAll();
}


void CloseAllAndQuit (TApplication* theApp, TFile* NtplFile)
{
  fileFitResults.close();
  fileFitSystematics.close();
  
  if (TotalPDFPsi                != NULL) DeleteFit(TotalPDFPsi,"All");
  if (TotalPDFRejectPsi          != NULL) DeleteFit(TotalPDFRejectPsi,"All");
  if (SingleCandNTuple_JPsi      != NULL) delete SingleCandNTuple_JPsi;
  if (SingleCandNTuple_PsiP      != NULL) delete SingleCandNTuple_PsiP;
  if (SingleCandNTuple_RejectPsi != NULL) delete SingleCandNTuple_RejectPsi;
  if (NtplFile                   != NULL) NtplFile->Close("R");

  for (unsigned int i = 0; i < effFuncs.first->size(); i++) delete effFuncs.first->operator[](i);
  effFuncs.first->erase(effFuncs.first->begin(),effFuncs.first->end());
  delete effFuncs.first;

  for (unsigned int i = 0; i < effFuncs.second->size(); i++) delete effFuncs.second->operator[](i);
  effFuncs.second->erase(effFuncs.second->begin(),effFuncs.second->end());
  delete effFuncs.second;

  gROOT->CloseFiles();
  theApp->Terminate(0);
}


string MakeName (RooDataSet* data, unsigned int ID)
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
  TLatex* LumiTex1 = new TLatex(0.2,0.91,myString.str().c_str());
  LumiTex1->SetTextSize(0.05);
  LumiTex1->SetTextColor(kBlack);
  LumiTex1->SetNDC(true);
  if (myFrame == NULL) LumiTex1->DrawLatex(0.2,0.91,myString.str().c_str());
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

  double startNDCx = 0.826;
  double startNDCy = 0.935;
  TLine* line1 = new TLine(startNDCx-0.005, startNDCy, startNDCx, startNDCy);
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

  myString.clear(); myString.str("");
  myString << "s = 8 TeV";
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


void AddGaussConstraint (RooArgSet* vecConstr, RooAbsPdf* TotalPDF, string varName)
{
  stringstream myString;
  
  myString.clear(); myString.str("");
  myString << varName;

  RooRealVar* varConstr = GetVar(TotalPDF,myString.str().c_str());
  double mean  = TotalPDF->getVariables()->getRealValue(myString.str().c_str());
  double sigma = (fabs(varConstr->getErrorLo()) + fabs(varConstr->getErrorHi())) / 2.;
  
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
      exit (EXIT_FAILURE);
    }

  myString.clear(); myString.str("");
  myString << varName;
  myString << "_constr";

  myCoeff.clear(); myCoeff.str("");
  myCoeff << "(1-TMath::Erf((-" << GetVar(TotalPDF,name)->getPlotLabel() << "-" << myConstr.str() << ") / ";
  myCoeff << "(" << Utility->GetGenericParam("SlewRateConstr").c_str() << " * 2*" << myConstr.str() << "))) * ";
  myCoeff << "(1-TMath::Erf(( " << GetVar(TotalPDF,name)->getPlotLabel() << "-" << myConstr.str() << ") / ";
  myCoeff << "(" << Utility->GetGenericParam("SlewRateConstr").c_str() << " * 2*" << myConstr.str() << ")))";
  cout << "\n[ExtractYield::AddPhysicsConstraint]\t@@@ Physical constraint on Afb/As boundaries: " << myCoeff.str().c_str() << " @@@" << endl;
  
  if (varName == "AfbBound")
    newConstr = new RooGenericPdf(myString.str().c_str(),myCoeff.str().c_str(),RooArgSet(*GetVar(TotalPDF,name),*GetVar(TotalPDF,"FlS")));
  else if (varName == "AsBound")
    newConstr = new RooGenericPdf(myString.str().c_str(),myCoeff.str().c_str(),RooArgSet(*GetVar(TotalPDF,name),*GetVar(TotalPDF,"FlS"),*GetVar(TotalPDF,"FsS")));
  vecConstr->add(*newConstr);
}


void BuildMassConstraints (RooArgSet* vecConstr, RooAbsPdf* TotalPDF, string varName)
// ##############################################################
// # All:    apply constraint to all physical angular variables #
// # sign:   apply constraint to signal                         #
// # mistag: apply constraint to mistag                         #
// # peak:   apply constraint to peaking background             #
// ##############################################################
{
  if ((GetVar(TotalPDF,"sigmaS1")        != NULL) && ((varName == "All") || (varName == "sign"))) AddGaussConstraint(vecConstr, TotalPDF, "sigmaS1");
  if ((GetVar(TotalPDF,"sigmaS2")        != NULL) && ((varName == "All") || (varName == "sign"))) AddGaussConstraint(vecConstr, TotalPDF, "sigmaS2");
  if ((GetVar(TotalPDF,"fracMassS")      != NULL) && ((varName == "All") || (varName == "sign"))) AddGaussConstraint(vecConstr, TotalPDF, "fracMassS");
  
  if ((GetVar(TotalPDF,"sigmaMisTag1")   != NULL) && ((varName == "All") || (varName == "mistag"))) AddGaussConstraint(vecConstr, TotalPDF, "sigmaMisTag1");
  if ((GetVar(TotalPDF,"sigmaMisTag2")   != NULL) && ((varName == "All") || (varName == "mistag"))) AddGaussConstraint(vecConstr, TotalPDF, "sigmaMisTag2");
  if ((GetVar(TotalPDF,"fracMisTag")     != NULL) && ((varName == "All") || (varName == "mistag"))) AddGaussConstraint(vecConstr, TotalPDF, "fracMisTag");

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
  
  if ((atoi(Utility->GetGenericParam("CtrlMisTagWrkFlow").c_str()) == 1) && (strcmp(CTRLfitWRKflow.c_str(),"trueAll&FitFrac") != 0))
    if ((GetVar(TotalPDF,"nMisTagFrac") != NULL) && ((varName == "All") || (varName == "mistag"))) AddGaussConstraint(vecConstr, TotalPDF, "nMisTagFrac");
  if ((GetVar(TotalPDF,"nBkgPeak") != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, TotalPDF, "nBkgPeak");
}


void BuildAngularConstraints (RooArgSet* vecConstr, RooAbsPdf* TotalPDF, string varName)
// ##############################################################
// # comb: apply constraint to angular combinatorial background #
// # peak: apply constraint to angular peaking background       #
// ##############################################################
{
  stringstream myString;
  string polyType;

  if (varName == "peak") polyType = "p";
  else                   polyType = "c";

  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
    {
      myString.clear(); myString.str("");
      myString << polyType.c_str() << "1Poly" << i;
      if (GetVar(TotalPDF,myString.str().c_str()) != NULL) AddGaussConstraint(vecConstr, TotalPDF, myString.str().c_str());
    }

  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
    {
      myString.clear(); myString.str("");
      myString << polyType.c_str() << "2Poly" << i;
      if (GetVar(TotalPDF,myString.str().c_str()) != NULL) AddGaussConstraint(vecConstr, TotalPDF, myString.str().c_str());
    }

  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
    {
      myString.clear(); myString.str("");
      myString << polyType.c_str() << "3Poly" << i;
      if (GetVar(TotalPDF,myString.str().c_str()) != NULL) AddGaussConstraint(vecConstr, TotalPDF, myString.str().c_str());
    }
}


void BuildPhysicsConstraints (RooArgSet* vecConstr, RooAbsPdf* TotalPDF, string varName)
{
  if ((GetVar(TotalPDF,"FlS")   != NULL) && (varName == "FlS"))  AddGaussConstraint(vecConstr, TotalPDF, varName);
  if ((GetVar(TotalPDF,"AfbS")  != NULL) && (varName == "AfbS")) AddGaussConstraint(vecConstr, TotalPDF, varName);

  if ((GetVar(TotalPDF,"FsS")   != NULL) && (varName == "FsS"))  AddGaussConstraint(vecConstr, TotalPDF, varName);
  if ((GetVar(TotalPDF,"AsS")   != NULL) && (varName == "AsS"))  AddGaussConstraint(vecConstr, TotalPDF, varName);


  // #########################################
  // # Special constraint for AFB boundaries #
  // #########################################
  if ((varName == "AfbBound") && (GetVar(TotalPDF,"FlS") != NULL) && (GetVar(TotalPDF,"AfbS") != NULL)) AddPhysicsConstraint(vecConstr, TotalPDF, varName);


  // ########################################
  // # Special constraint for AS boundaries #
  // ########################################
  if ((varName == "AsBound") && (GetVar(TotalPDF,"FlS") != NULL) && (GetVar(TotalPDF,"FsS") != NULL) && (GetVar(TotalPDF,"AsS") != NULL)) AddPhysicsConstraint(vecConstr, TotalPDF, varName);
}


RooAbsPdf* MakeAngWithEffPDF (TF2* effFunc, RooRealVar* y, RooRealVar* z, unsigned int FitType, bool useEffPDF, RooArgSet* VarsAng, RooArgSet* VarsPoly, vector<double>* q2Bins, int q2BinIndx)
// ###################
// # y: cos(theta_l) #
// # z: cos(theta_K) #
// ###################
{
  stringstream myString;
  vector<RooRealVar*> vecParam;
  RooAbsPdf* AnglesPDF = NULL;

  if ((FitType == 1) || (FitType == 41) || (FitType == 61) || (FitType == 81) || // Branching fraction
      (FitType == 6) || (FitType == 26) || (FitType == 36) || (FitType == 46) || (FitType == 56) || (FitType == 66) || (FitType == 76) || (FitType == 86) || (FitType == 96)) // Fl-Afb-fit
    {
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
      
      // #######################################################
      // # Make 2D signal*efficiency p.d.f.: integral over phi #
      // # For correctly tagged events                         #
      // #######################################################
      // @TMP@
      FlS  = new RooRealVar("FlS","F_{L}",0.0);
      AfbS = new RooRealVar("AfbS","A_{FB}",0.0);
      // FlS  = new RooRealVar("FlS","F_{L}",0.0,0.0,1.0);
      // AfbS = new RooRealVar("AfbS","A_{FB}",0.0,-1.0,1.0);
      VarsAng->add(*FlS);
      VarsAng->add(*AfbS);
      VarsAng->add(*y);
      VarsAng->add(*z);

      myString.clear(); myString.str("");
      if (atoi(Utility->GetGenericParam("UseSPwave").c_str()) == false)
    	{
    	  // #####################
    	  // # P-wave decay rate #
    	  // #####################
	  // @TMP@
    	  myString << "(3/4 * (3/2 * (1/2+TMath::ATan(FlS)/TMath::Pi()) * (1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * " << z->getPlotLabel() << "*" << z->getPlotLabel() << " + ";
    	  myString << "(3/8 * (1-(1/2+TMath::ATan(FlS)/TMath::Pi())) * (1+" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") + (3/4*(1-(1/2+TMath::ATan(FlS)/TMath::Pi())))*2*TMath::ATan(AfbS)/TMath::Pi() *" << y->getPlotLabel() << ") * ";
    	  myString << "(1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")))";
    	  // myString << "(3/4 * (3/2 * FlS * (1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * " << z->getPlotLabel() << "*" << z->getPlotLabel() << " + ";
    	  // myString << "(3/8 * (1-FlS) * (1+" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") + AfbS*" << y->getPlotLabel() << ") * ";
    	  // myString << "(1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")))";
    	}
      else
    	{
    	  FsS = new RooRealVar("FsS","F_{s}",0.5,0.0,1.0);
    	  AsS = new RooRealVar("AsS","A_{s}",0.0,-1.0,1.0);
    	  VarsAng->add(*FsS);
    	  VarsAng->add(*AsS);

    	  // ###########################
    	  // # S and P-wave decay rate #
    	  // ###########################
   	  myString << "(9/16 * ((2/3*FsS + 4/3*AsS*" << z->getPlotLabel() << ") * (1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") + ";
    	  myString << "(1-FsS) * ";
    	  myString << "(2*FlS*" << z->getPlotLabel() << "*" << z->getPlotLabel() << " * (1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") + ";
    	  myString << "1/2*(1-FlS) * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") * (1+" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") + ";
    	  myString << "4/3*AfbS * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") * " << y->getPlotLabel() << ")))";
    	}

      VarsPoly->removeAll();
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
    	  myString << y->getPlotLabel() << " + ";
	  
    	  myString << "(P8 + P9*" << z->getPlotLabel() << " + P10*" << z->getPlotLabel() << "*" << z->getPlotLabel() << " + ";
    	  myString << "P11*" << z->getPlotLabel() << "*" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") * ";
    	  myString << y->getPlotLabel() << "*" << y->getPlotLabel() << " + ";
	  
    	  myString << "(P12 + P13*" << z->getPlotLabel() << " + P14*" << z->getPlotLabel() << "*" << z->getPlotLabel() << " + ";
    	  myString << "P15*" << z->getPlotLabel() << "*" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") * ";
    	  myString << y->getPlotLabel() << "*" << y->getPlotLabel() << "*" << y->getPlotLabel() << " + ";

    	  myString << "(P16 + P17*" << z->getPlotLabel() << " + P18*" << z->getPlotLabel() << "*" << z->getPlotLabel() << " + ";
    	  myString << "P19*" << z->getPlotLabel() << "*" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") * ";
    	  myString << y->getPlotLabel() << "*" << y->getPlotLabel() << "*" << y->getPlotLabel() << "*" << y->getPlotLabel() << " + ";

    	  myString << "(P20 + P21*" << z->getPlotLabel() << " + P22*" << z->getPlotLabel() << "*" << z->getPlotLabel() << " + ";
    	  myString << "P23*" << z->getPlotLabel() << "*" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") * ";
    	  myString << y->getPlotLabel() << "*" << y->getPlotLabel() << "*" << y->getPlotLabel() << "*" << y->getPlotLabel() << "*" << y->getPlotLabel() << ")";
	  
    	  for (int i = 0; i < effFunc->GetNpar(); i++) VarsPoly->add(*vecParam[i]);
    	}

      cout << "\n[ExtractYield::MakeAngWithEffPDF]\t@@@ 2D angular*efficiency p.d.f. @@@" << endl;
      cout << myString.str().c_str() << endl;

      AnglesPDF = new RooGenericPdf("AngleS",myString.str().c_str(),RooArgSet(*VarsAng,*VarsPoly));
    }
  else if ((FitType == 1*10) || (FitType == 41*10) || (FitType == 61*10) || (FitType == 81*10) || // Branching fraction
	   (FitType == 6*10) || (FitType == 26*10) || (FitType == 36*10) || (FitType == 46*10) || (FitType == 56*10) || (FitType == 66*10) || (FitType == 76*10) || (FitType == 86*10) || (FitType == 96*10)) // Fl-Afb-fit
    {
      // #######################################################
      // # Make 2D signal*efficiency p.d.f.: integral over phi #
      // # For incorrectly tagged events                       #
      // #######################################################

      myString.clear(); myString.str("");
      if (atoi(Utility->GetGenericParam("UseSPwave").c_str()) == false)
	{
	  // #####################
	  // # P-wave decay rate #
	  // #####################
	  // @TMP@
    	  myString << "(3/4 * (3/2 * (1/2+TMath::ATan(FlS)/TMath::Pi()) * (1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * " << z->getPlotLabel() << "*" << z->getPlotLabel() << " + ";
    	  myString << "(3/8 * (1-(1/2+TMath::ATan(FlS)/TMath::Pi())) * (1+" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") - (3/4*(1-(1/2+TMath::ATan(FlS)/TMath::Pi())))*2*TMath::ATan(AfbS)/TMath::Pi() *" << y->getPlotLabel() << ") * ";
    	  myString << "(1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")))";
    	  // myString << "(3/4 * (3/2 * FlS * (1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * " << z->getPlotLabel() << "*" << z->getPlotLabel() << " + ";
    	  // myString << "(3/8 * (1-FlS) * (1+" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") - AfbS*" << y->getPlotLabel() << ") * ";
    	  // myString << "(1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")))";
	}
      else
	{
	  // ###########################
	  // # S and P-wave decay rate #
	  // ###########################
    	  myString << "(9/16 * ((2/3*FsS - 4/3*AsS*" << z->getPlotLabel() << ") * (1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") + ";
    	  myString << "(1-FsS) * ";
    	  myString << "(2*FlS*" << z->getPlotLabel() << "*" << z->getPlotLabel() << " * (1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") + ";
    	  myString << "1/2*(1-FlS) * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") * (1+" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") - ";
    	  myString << "4/3*AfbS * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") * " << y->getPlotLabel() << ")))";
	}

      cout << "\n[ExtractYield::MakeAngWithEffPDF]\t@@@ 2D angular p.d.f. @@@" << endl;
      cout << myString.str().c_str() << endl;

      // @TMP@ : for the mis-tag I use the binned efficiency, and NOT the analytical efficiency, due to problems in fitting the binned efficiency
      VarsPoly->removeAll();
      RooHistPdf* histoEffPDF;
      if (useEffPDF == true)
	{
	  RooGenericPdf* _AnglesPDF = new RooGenericPdf("_AnglesPDF",myString.str().c_str(),RooArgSet(*VarsAng));

	  // #############################
	  // # Make 2D efficiency p.d.f. #
	  // #############################
	  RooDataHist* histoEff = new RooDataHist("histoEff","histoEff",RooArgSet(*z,*y),Import(*Utility->Get2DEffHistoq2Bin(&cosThetaKBins,&cosThetaLBins,q2BinIndx,GetSignalType(FitType,q2Bins,q2BinIndx),false,make_pair(-1.0,1.0),make_pair(-1.0,1.0)),true));
	  histoEffPDF           = new RooHistPdf("histoEffPDF","histoEffPDF",RooArgSet(*z,*y),*histoEff,DEGREEINTERP);
	  MyProdPdf* myprodpdf  = new MyProdPdf(*_AnglesPDF,*histoEffPDF);
	  ROOT::Math::Functor* prodFunctor = new ROOT::Math::Functor(*myprodpdf,myprodpdf->ndim());
	  AnglesPDF             = new RooFunctorPdfBinding("AngleMisTag","MisTag * Efficiency",*prodFunctor,myprodpdf->vars());
	}
      else AnglesPDF = new RooGenericPdf("AngleMisTag",myString.str().c_str(),RooArgSet(*VarsAng));
    }


  vecParam.clear();
  return AnglesPDF;
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

	  if (GetVar(TotalPDF,"var1")           != NULL) delete GetVar(TotalPDF,"var1");
	  if (GetVar(TotalPDF,"var2")           != NULL) delete GetVar(TotalPDF,"var2");
	  if (GetVar(TotalPDF,"fracMassBExp")   != NULL) delete GetVar(TotalPDF,"fracMassBExp");

	  if (GetVar(TotalPDF,"sigmaMisTag1")   != NULL) delete GetVar(TotalPDF,"sigmaMisTag1");
	  if (GetVar(TotalPDF,"sigmaMisTag2")   != NULL) delete GetVar(TotalPDF,"sigmaMisTag2");
	  if (GetVar(TotalPDF,"fracMisTag")     != NULL) delete GetVar(TotalPDF,"fracMisTag");

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
	  if (GetVar(TotalPDF,"nMisTagFrac")    != NULL) delete GetVar(TotalPDF,"nMisTagFrac");
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
	  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
	    {
	      myString.clear(); myString.str("");
	      myString << "p3Poly" << i;
	      if (GetVar(TotalPDF,myString.str().c_str()) != NULL) delete GetVar(TotalPDF,myString.str().c_str());
	    } 
	  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
	    {
	      myString.clear(); myString.str("");
	      myString << "c3Poly" << i;
	      if (GetVar(TotalPDF,myString.str().c_str()) != NULL) delete GetVar(TotalPDF,myString.str().c_str());
	    }
	  if (GetVar(TotalPDF,"FlS")  != NULL) delete GetVar(TotalPDF,"FlS");
	  if (GetVar(TotalPDF,"AfbS") != NULL) delete GetVar(TotalPDF,"AfbS");
	  if (GetVar(TotalPDF,"P1S")  != NULL) delete GetVar(TotalPDF,"P1S");
	  if (GetVar(TotalPDF,"P2S")  != NULL) delete GetVar(TotalPDF,"P2S");
	  if (GetVar(TotalPDF,"FsS")  != NULL) delete GetVar(TotalPDF,"FsS");
	  if (GetVar(TotalPDF,"AsS")  != NULL) delete GetVar(TotalPDF,"AsS");
	}
      else if (DeleteType == "All")
	{
	  if (TotalPDF->getComponents()->find("MassS1")             != NULL) delete TotalPDF->getComponents()->find("MassS1");
	  if (TotalPDF->getComponents()->find("MassS2")             != NULL) delete TotalPDF->getComponents()->find("MassS2");
	  if (TotalPDF->getComponents()->find("MassSignal")         != NULL) delete TotalPDF->getComponents()->find("MassSignal");
	  if (TotalPDF->getComponents()->find("AngleS")             != NULL) delete TotalPDF->getComponents()->find("AngleS");
	  if (TotalPDF->getComponents()->find("Signal")             != NULL) delete TotalPDF->getComponents()->find("Signal");

	  if (TotalPDF->getComponents()->find("EffPDFsign")         != NULL) delete TotalPDF->getComponents()->find("EffPDFsign");
	  if (TotalPDF->getComponents()->find("EffPDFnorm")         != NULL) delete TotalPDF->getComponents()->find("EffPDFnorm");

	  if (TotalPDF->getComponents()->find("BkgMassExp1")        != NULL) delete TotalPDF->getComponents()->find("BkgMassExp1");
	  if (TotalPDF->getComponents()->find("BkgMassExp2")        != NULL) delete TotalPDF->getComponents()->find("BkgMassExp2");

	  if (TotalPDF->getComponents()->find("MassMisTag1")        != NULL) delete TotalPDF->getComponents()->find("MassMisTag1");
	  if (TotalPDF->getComponents()->find("MassMisTag2")        != NULL) delete TotalPDF->getComponents()->find("MassMisTag2");
	  if (TotalPDF->getComponents()->find("MassMisTag")         != NULL) delete TotalPDF->getComponents()->find("MassMisTag");
	  if (TotalPDF->getComponents()->find("AngleMisTag")        != NULL) delete TotalPDF->getComponents()->find("AngleMisTag");

	  if (TotalPDF->getComponents()->find("BkgMassRPeak1")      != NULL) delete TotalPDF->getComponents()->find("BkgMassRPeak1");
	  if (TotalPDF->getComponents()->find("BkgMassRPeak2")      != NULL) delete TotalPDF->getComponents()->find("BkgMassRPeak1");
	  if (TotalPDF->getComponents()->find("BkgMassRPeak")       != NULL) delete TotalPDF->getComponents()->find("BkgMassRPeak1");
	  if (TotalPDF->getComponents()->find("BkgMassLPeak1")      != NULL) delete TotalPDF->getComponents()->find("BkgMassLPeak1");
	  if (TotalPDF->getComponents()->find("BkgMassLPeak2")      != NULL) delete TotalPDF->getComponents()->find("BkgMassLPeak1");
	  if (TotalPDF->getComponents()->find("BkgMassLPeak")       != NULL) delete TotalPDF->getComponents()->find("BkgMassLPeak1");
	  if (TotalPDF->getComponents()->find("BkgMassComb")        != NULL) delete TotalPDF->getComponents()->find("BkgMassComb");
	  if (TotalPDF->getComponents()->find("BkgMassPeak")        != NULL) delete TotalPDF->getComponents()->find("BkgMassPeak");
	  if (TotalPDF->getComponents()->find("BkgAnglesC")         != NULL) delete TotalPDF->getComponents()->find("BkgAnglesC");
	  if (TotalPDF->getComponents()->find("BkgAnglesP")         != NULL) delete TotalPDF->getComponents()->find("BkgAnglesP");
	  if (TotalPDF->getComponents()->find("BkgMassAngleComb")   != NULL) delete TotalPDF->getComponents()->find("BkgMassAngleComb");
	  if (TotalPDF->getComponents()->find("MassAngleMisTag")    != NULL) delete TotalPDF->getComponents()->find("MassAngleMisTag");
	  if (TotalPDF->getComponents()->find("BkgMassAnglePeak")   != NULL) delete TotalPDF->getComponents()->find("BkgMassAnglePeak");

	  if (TotalPDF != NULL) delete TotalPDF;
	}
      else
	{
	  cout << "[ExtractYield::DeleteFit]\tWrong parameter: " << DeleteType << endl;
	  exit (EXIT_FAILURE);
	}
    }
}


void ResetCombPolyParam (vector<vector<string>*>* fitParam)
{
  stringstream myString;

  for (unsigned int j = 0; j < fitParam->operator[](0)->size(); j++)
    {
      cout << "\n[ExtractYield::ResetCombPolyParam]\t@@@ Resetting the combinatorial background polynomial coefficients for q^2 bin #" << j << " @@@" << endl;

      for (int i = 0; i < atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC1"))->operator[](j).c_str()); i++)
        {
          myString.clear(); myString.str("");
          myString << "c1Poly" << i;
          fitParam->operator[](Utility->GetFitParamIndx(myString.str().c_str()))->operator[](j) = "0.0";

          cout << myString.str().c_str() << "\t" << fitParam->operator[](Utility->GetFitParamIndx(myString.str().c_str()))->operator[](j).c_str() << endl;
         }

      for (int i = 0; i < atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC2"))->operator[](j).c_str()); i++)
        {
          myString.clear(); myString.str("");
          myString << "c2Poly" << i;
          fitParam->operator[](Utility->GetFitParamIndx(myString.str().c_str()))->operator[](j) = "0.0";

          cout << myString.str().c_str() << "\t" << fitParam->operator[](Utility->GetFitParamIndx(myString.str().c_str()))->operator[](j).c_str() << endl;
        }

      for (int i = 0; i < atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC3"))->operator[](j).c_str()); i++)
        {
          myString.clear(); myString.str("");
          myString << "c3Poly" << i;
          fitParam->operator[](Utility->GetFitParamIndx(myString.str().c_str()))->operator[](j) = "0.0";

          cout << myString.str().c_str() << "\t" << fitParam->operator[](Utility->GetFitParamIndx(myString.str().c_str()))->operator[](j).c_str() << endl;
        }
    }
}


void ResetAngularParam (vector<vector<string>*>* fitParam)
{
  for (unsigned int j = 0; j < fitParam->operator[](0)->size(); j++)
    {
      cout << "\n[ExtractYield::ResetAngularParam]\t@@@ Resetting the angular parameters for q^2 bin #" << j << " @@@" << endl;

      fitParam->operator[](Utility->GetFitParamIndx("FlS"))->operator[](j)   = "0.5   0.0   0.0";
      fitParam->operator[](Utility->GetFitParamIndx("AfbS"))->operator[](j)  = "0.0   0.0   0.0";
      fitParam->operator[](Utility->GetFitParamIndx("P1S"))->operator[](j)   = "0.0   0.0   0.0";
      fitParam->operator[](Utility->GetFitParamIndx("P2S"))->operator[](j)   = "0.0   0.0   0.0";
      fitParam->operator[](Utility->GetFitParamIndx("FsS"))->operator[](j)   = "0.05   0.02   0.02";
      fitParam->operator[](Utility->GetFitParamIndx("AsS"))->operator[](j)   = "0.0   0.0   0.0";

      cout << "FL: "   << "\t" << fitParam->operator[](Utility->GetFitParamIndx("FlS"))->operator[](j).c_str() << endl;
      cout << "AFB: "  << "\t" << fitParam->operator[](Utility->GetFitParamIndx("AfbS"))->operator[](j).c_str() << endl;
      cout << "P1: "   << "\t" << fitParam->operator[](Utility->GetFitParamIndx("P1S"))->operator[](j).c_str() << endl;
      cout << "P2: "   << "\t" << fitParam->operator[](Utility->GetFitParamIndx("P2S"))->operator[](j).c_str() << endl;
      cout << "FS: "   << "\t" << fitParam->operator[](Utility->GetFitParamIndx("FsS"))->operator[](j).c_str() << endl;
      cout << "AS: "   << "\t" << fitParam->operator[](Utility->GetFitParamIndx("AsS"))->operator[](j).c_str() << endl;
    }
}


double StoreFitResultsInFile (RooAbsPdf** TotalPDF, RooFitResult* fitResult, RooDataSet* dataSet, RooArgSet* vecConstr)
{
  RooAbsReal* NLL;
  double signalSigmaGoodT  = 0.0;
  double signalSigmaGoodTE = 0.0;
  double signalSigmaMisT   = 0.0;
  double signalSigmaMisTE  = 0.0;
  double NLLvalue    = 0.0;


  if (fitResult != NULL)
    {
      fileFitResults << "Covariance quality (3=ok): " << fitResult->covQual() << endl;
      fileFitResults << "Fit status (0=ok): " << fitResult->status() << endl;
      fileFitResults << "Fit EDM: " << fitResult->edm() << endl;
      if (atoi(Utility->GetGenericParam("ApplyConstr").c_str()) == true) NLL = (*TotalPDF)->createNLL(*dataSet,Extended(true),ExternalConstraints(*vecConstr));
      else                                                               NLL = (*TotalPDF)->createNLL(*dataSet,Extended(true));
      NLLvalue = NLL->getVal();
      delete NLL;
      fileFitResults << "NLL: " << NLLvalue << endl;


      // ##########
      // # Signal #
      // ##########
      if (GetVar(*TotalPDF,"fracMassS") != NULL)
	{
	  signalSigmaGoodT  = sqrt( GetVar(*TotalPDF,"fracMassS")->getVal() * pow(GetVar(*TotalPDF,"sigmaS1")->getVal(),2.) + (1. - GetVar(*TotalPDF,"fracMassS")->getVal()) * pow(GetVar(*TotalPDF,"sigmaS2")->getVal(),2.) );
	  signalSigmaGoodTE = 1./(2.*signalSigmaGoodT) * sqrt( pow((pow(GetVar(*TotalPDF,"sigmaS1")->getVal(),2.) - pow(GetVar(*TotalPDF,"sigmaS2")->getVal(),2.)) * GetVar(*TotalPDF,"fracMassS")->getError(),2.) +
							       pow(2. * GetVar(*TotalPDF,"fracMassS")->getVal()        * GetVar(*TotalPDF,"sigmaS1")->getVal() * GetVar(*TotalPDF,"sigmaS1")->getError(),2.) +
							       pow(2. * (1. - GetVar(*TotalPDF,"fracMassS")->getVal()) * GetVar(*TotalPDF,"sigmaS2")->getVal() * GetVar(*TotalPDF,"sigmaS2")->getError(),2.) );
	}
      else if (GetVar(*TotalPDF,"sigmaS1") != NULL)
	{
	  signalSigmaGoodT  = GetVar(*TotalPDF,"sigmaS1")->getVal();
	  signalSigmaGoodTE = GetVar(*TotalPDF,"sigmaS1")->getError();
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
      if (GetVar(*TotalPDF,"sigmaS1") != NULL) fileFitResults << "< Sigma >: " << signalSigmaGoodT << " +/- " << signalSigmaGoodTE << endl;
      if (GetVar(*TotalPDF,"nSig") != NULL)
	{
	  fileFitResults << "Signal yield: " << GetVar(*TotalPDF,"nSig")->getVal() << " +/- " << GetVar(*TotalPDF,"nSig")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"nSig")->getErrorHi() << "/" << GetVar(*TotalPDF,"nSig")->getErrorLo() << ")" << endl;
	}
      

      // ###################
      // # Signal mist-tag #
      // ###################
      if (GetVar(*TotalPDF,"fracMisTag") != NULL)
	{
	  signalSigmaMisT  = sqrt( GetVar(*TotalPDF,"fracMisTag")->getVal() * pow(GetVar(*TotalPDF,"sigmaMisTag1")->getVal(),2.) + (1. - GetVar(*TotalPDF,"fracMisTag")->getVal()) * pow(GetVar(*TotalPDF,"sigmaMisTag2")->getVal(),2.) );
	  signalSigmaMisTE = 1./(2.*signalSigmaMisT) * sqrt( pow((pow(GetVar(*TotalPDF,"sigmaMisTag1")->getVal(),2.) - pow(GetVar(*TotalPDF,"sigmaMisTag2")->getVal(),2.)) * GetVar(*TotalPDF,"fracMisTag")->getError(),2.) +
							     pow(2. * GetVar(*TotalPDF,"fracMisTag")->getVal()        * GetVar(*TotalPDF,"sigmaMisTag1")->getVal() * GetVar(*TotalPDF,"sigmaMisTag1")->getError(),2.) +
							     pow(2. * (1. - GetVar(*TotalPDF,"fracMisTag")->getVal()) * GetVar(*TotalPDF,"sigmaMisTag2")->getVal() * GetVar(*TotalPDF,"sigmaMisTag2")->getError(),2.) );
	}
      else if (GetVar(*TotalPDF,"sigmaMisTag1") != NULL)
	{
	  signalSigmaMisT  = GetVar(*TotalPDF,"sigmaMisTag1")->getVal();
	  signalSigmaMisTE = GetVar(*TotalPDF,"sigmaMisTag1")->getError();
	}

      if (GetVar(*TotalPDF,"nMisTagFrac") != NULL) fileFitResults << "Background mistag mean: equal to signal mean" << endl;
      if (GetVar(*TotalPDF,"sigmaMisTag1") != NULL)
	{
	  fileFitResults << "Background mistag sigma-1: " << GetVar(*TotalPDF,"sigmaMisTag1")->getVal() << " +/- " << GetVar(*TotalPDF,"sigmaMisTag1")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"sigmaMisTag1")->getErrorHi() << "/" << GetVar(*TotalPDF,"sigmaMisTag1")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"sigmaMisTag2") != NULL)
	{
	  fileFitResults << "Background mistag sigma-2: " << GetVar(*TotalPDF,"sigmaMisTag2")->getVal() << " +/- " << GetVar(*TotalPDF,"sigmaMisTag2")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"sigmaMisTag2")->getErrorHi() << "/" << GetVar(*TotalPDF,"sigmaMisTag2")->getErrorLo() << ")" << endl;
	  fileFitResults << "Fraction: " << GetVar(*TotalPDF,"fracMisTag")->getVal() << " +/- " << GetVar(*TotalPDF,"fracMisTag")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"fracMisTag")->getErrorHi() << "/" << GetVar(*TotalPDF,"fracMisTag")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"sigmaMisTag1") != NULL) fileFitResults << "< Sigma >: " << signalSigmaMisT << " +/- " << signalSigmaMisTE << endl;
      if (GetVar(*TotalPDF,"nMisTagFrac") != NULL)
	{
	  fileFitResults << "Mistag fraction: " << GetVar(*TotalPDF,"nMisTagFrac")->getVal() << " +/- " << GetVar(*TotalPDF,"nMisTagFrac")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"nMisTagFrac")->getErrorHi() << "/" << GetVar(*TotalPDF,"nMisTagFrac")->getErrorLo() << ")" << endl;
	}


      // ################################
      // # Total signal yield and width #
      // ################################
      if ((GetVar(*TotalPDF,"nSig") != NULL) && (GetVar(*TotalPDF,"nMisTagFrac") != NULL))
	{
	  double totalYield    = GetVar(*TotalPDF,"nSig")->getVal() / (1. - GetVar(*TotalPDF,"nMisTagFrac")->getVal());
	  double totalYieldE   = totalYield * sqrt( pow(GetVar(*TotalPDF,"nSig")->getError()/GetVar(*TotalPDF,"nSig")->getVal(),2.) + pow(GetVar(*TotalPDF,"nMisTagFrac")->getError()/(1. - GetVar(*TotalPDF,"nMisTagFrac")->getVal()),2.) );

	  double fraction      = (1. - 2.*GetVar(*TotalPDF,"nMisTagFrac")->getVal()) / (1. - GetVar(*TotalPDF,"nMisTagFrac")->getVal());
	  double fractionE     = GetVar(*TotalPDF,"nMisTagFrac")->getError() / pow((1. - GetVar(*TotalPDF,"nMisTagFrac")->getVal()),2.);
	  
	  double signalSigmaT  = sqrt( fraction * pow(signalSigmaGoodT,2.) + (1. - fraction) * pow(signalSigmaMisT,2.) );
	  double signalSigmaTE = 1./(2.*signalSigmaT) * sqrt( pow((pow(signalSigmaGoodT,2.) - pow(signalSigmaMisT,2.)) * fractionE,2.) +
							      pow(2. * fraction        * signalSigmaGoodT * signalSigmaGoodTE,2.) +
							      pow(2. * (1. - fraction) * signalSigmaMisT  * signalSigmaMisTE,2.) );
	  
	  fileFitResults << "Total signal yield: " << totalYield << " +/- " << totalYieldE << endl;
	  fileFitResults << "< Total Sigma >: " << signalSigmaT << " +/- " << signalSigmaTE << endl;
	}


      // ############################
      // # Combinatorial background #
      // ############################
      if (GetVar(*TotalPDF,"var1") != NULL)
	{
	  fileFitResults << "Background mass var-1: " << GetVar(*TotalPDF,"var1")->getVal() << " +/- " << var1->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"var1")->getErrorHi() << "/" << GetVar(*TotalPDF,"var1")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"var2") != NULL)
	{
	  fileFitResults << "Background mass var-2: " << GetVar(*TotalPDF,"var2")->getVal() << " +/- " << GetVar(*TotalPDF,"var2")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"var2")->getErrorHi() << "/" << GetVar(*TotalPDF,"var2")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"fracMassBExp") != NULL)
	{
	  fileFitResults << "Fraction: " << GetVar(*TotalPDF,"fracMassBExp")->getVal() << " +/- " << GetVar(*TotalPDF,"fracMassBExp")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"fracMassBExp")->getErrorHi() << "/" << GetVar(*TotalPDF,"fracMassBExp")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"nBkgComb") != NULL)
	{
	  fileFitResults << "Comb. bkg yield: " << GetVar(*TotalPDF,"nBkgComb")->getVal() << " +/- " << GetVar(*TotalPDF,"nBkgComb")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"nBkgComb")->getErrorHi() << "/" << GetVar(*TotalPDF,"nBkgComb")->getErrorLo() << ")" << endl;
	}


      // ######################
      // # Peaking background #
      // ######################
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
	  fileFitResults << "Peaking bkg yield: " << GetVar(*TotalPDF,"nBkgPeak")->getVal() << " +/- " << GetVar(*TotalPDF,"nBkgPeak")->getError();
	  fileFitResults << " (" << GetVar(*TotalPDF,"nBkgPeak")->getErrorHi() << "/" << GetVar(*TotalPDF,"nBkgPeak")->getErrorLo() << ")" << endl;
	}


      // #######################
      // # Angular observables #
      // #######################
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
      if (GetVar(*TotalPDF,"P1S") != NULL)
	{
	  fileFitResults << "P1: " << GetVar(*TotalPDF,"P1S")->getVal() << " +/- ";
	  fileFitResults << GetVar(*TotalPDF,"P1S")->getError() << " (" << GetVar(*TotalPDF,"P1S")->getErrorHi() << "/" << GetVar(*TotalPDF,"P1S")->getErrorLo() << ")" << endl;
	}
      if (GetVar(*TotalPDF,"P2S") != NULL)
	{
	  fileFitResults << "P2: " << GetVar(*TotalPDF,"P2S")->getVal() << " +/- ";
	  fileFitResults << GetVar(*TotalPDF,"P2S")->getError() << " (" << GetVar(*TotalPDF,"P2S")->getErrorHi() << "/" << GetVar(*TotalPDF,"P2S")->getErrorLo() << ")" << endl;
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

      for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
	{
	  myString.clear(); myString.str("");
	  myString << "c3Poly" << i;
	  if (GetVar(*TotalPDF,myString.str().c_str()) != NULL)
	    {
	      fileFitResults << "Comb.bkg angle-3 c" << i << ": " << GetVar(*TotalPDF,myString.str().c_str())->getVal() << " +/- " << GetVar(*TotalPDF,myString.str().c_str())->getError();
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

      for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
	{
	  myString.clear(); myString.str("");
	  myString << "p3Poly" << i;
	  if (GetVar(*TotalPDF,myString.str().c_str()) != NULL)
	    {
	      fileFitResults << "Peak.bkg angle-3 p" << i << ": " << GetVar(*TotalPDF,myString.str().c_str())->getVal() << " +/- " << GetVar(*TotalPDF,myString.str().c_str())->getError();
	      fileFitResults << " (" << GetVar(*TotalPDF,myString.str().c_str())->getErrorHi() << "/" << GetVar(*TotalPDF,myString.str().c_str())->getErrorHi() << ")" << endl;
	    }
	}
    }
}


vector<string>* SaveFitResults (RooAbsPdf* TotalPDF, unsigned int q2BinIndx, vector<vector<string>*>* fitParam, vector<vector<unsigned int>*>* configParam, RooArgSet* vecConstr)
{
  stringstream myString;
  stringstream myCoeff;
  vector<string>* vecParStr;
  unsigned int NCoeffPolyBKGpeak1 = 0;
  unsigned int NCoeffPolyBKGcomb1 = 0;
  unsigned int NCoeffPolyBKGpeak2 = 0;
  unsigned int NCoeffPolyBKGcomb2 = 0;
  unsigned int NCoeffPolyBKGpeak3 = 0;
  unsigned int NCoeffPolyBKGcomb3 = 0;

  vecParStr = new vector<string>;


  vecParStr->push_back("# Signal mean [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"meanS") != NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("meanS") << "   " << GetVar(TotalPDF,"meanS")->getErrorLo() << "   " << GetVar(TotalPDF,"meanS")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("meanS"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# Signal sigma-1 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"sigmaS1") != NULL) && (vecConstr->find(string(string("sigmaS1") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("sigmaS1") << "   " << GetVar(TotalPDF,"sigmaS1")->getErrorLo() << "   " << GetVar(TotalPDF,"sigmaS1")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# Signal sigma-2 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"sigmaS2") != NULL) && (vecConstr->find(string(string("sigmaS2") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("sigmaS2") << "   " << GetVar(TotalPDF,"sigmaS2")->getErrorLo() << "   " << GetVar(TotalPDF,"sigmaS2")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("sigmaS2"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# Fraction");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"fracMassS") != NULL) && (vecConstr->find(string(string("fracMassS") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("fracMassS") << "   " << GetVar(TotalPDF,"fracMassS")->getErrorLo() << "   " << GetVar(TotalPDF,"fracMassS")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("fracMassS"))->operator[](q2BinIndx).c_str());


  vecParStr->push_back("###################################################################");


  vecParStr->push_back("# Bkg var-1 (tau1 or mean) [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"var1") != NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("var1") << "   " << GetVar(TotalPDF,"var1")->getErrorLo() << "   " << GetVar(TotalPDF,"var1")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("var1"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# Bkg var-2 (tau2 or width) [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"var2") != NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("var2") << "   " << GetVar(TotalPDF,"var2")->getErrorLo() << "   " << GetVar(TotalPDF,"var2")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("var2"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# Fraction");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"fracMassBExp") != NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("fracMassBExp") << "   " << GetVar(TotalPDF,"fracMassBExp")->getErrorLo() << "   " << GetVar(TotalPDF,"fracMassBExp")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("fracMassBExp"))->operator[](q2BinIndx).c_str());


  vecParStr->push_back("###################################################################");


  vecParStr->push_back("# Sigma-1 mistag [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"sigmaMisTag1") != NULL) && (vecConstr->find(string(string("sigmaMisTag1") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("sigmaMisTag1") << "   " << GetVar(TotalPDF,"sigmaMisTag1")->getErrorLo() << "   " << GetVar(TotalPDF,"sigmaMisTag1")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag1"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# Sigma-2 mistag [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"sigmaMisTag2") != NULL) && (vecConstr->find(string(string("sigmaMisTag2") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("sigmaMisTag2") << "   " << GetVar(TotalPDF,"sigmaMisTag2")->getErrorLo() << "   " << GetVar(TotalPDF,"sigmaMisTag2")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag2"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# Fraction");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"fracMisTag") != NULL) && (vecConstr->find(string(string("fracMisTag") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("fracMisTag") << "   " << GetVar(TotalPDF,"fracMisTag")->getErrorLo() << "   " << GetVar(TotalPDF,"fracMisTag")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("fracMisTag"))->operator[](q2BinIndx).c_str());


  vecParStr->push_back("###################################################################");


  vecParStr->push_back("# Bkg mean right-peak-1 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"meanR1") != NULL) && (vecConstr->find(string(string("meanR1") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("meanR1") << "   " << GetVar(TotalPDF,"meanR1")->getErrorLo() << "   " << GetVar(TotalPDF,"meanR1")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("meanR1"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# Bkg sigma right-peak-1 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"sigmaR1") != NULL) && (vecConstr->find(string(string("sigmaR1") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("sigmaR1") << "   " << GetVar(TotalPDF,"sigmaR1")->getErrorLo() << "   " << GetVar(TotalPDF,"sigmaR1")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# Bkg mean right-peak-2 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"meanR2") != NULL) && (vecConstr->find(string(string("meanR2") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("meanR2") << "   " << GetVar(TotalPDF,"meanR2")->getErrorLo() << "   " << GetVar(TotalPDF,"meanR2")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("meanR2"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# Bkg sigma right-peak-2 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"sigmaR2") != NULL) && (vecConstr->find(string(string("sigmaR2") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("sigmaR2") << "   " << GetVar(TotalPDF,"sigmaR2")->getErrorLo() << "   " << GetVar(TotalPDF,"sigmaR2")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# Fraction");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"fracMassBRPeak") != NULL) && (vecConstr->find(string(string("fracMassBRPeak") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("fracMassBRPeak") << "   " << GetVar(TotalPDF,"fracMassBRPeak")->getErrorLo() << "   " << GetVar(TotalPDF,"fracMassBRPeak")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("fracMassBRPeak"))->operator[](q2BinIndx).c_str());


  vecParStr->push_back("# Bkg mean left-peak-1 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"meanL1") != NULL) && (vecConstr->find(string(string("meanL1") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("meanL1") << "   " << GetVar(TotalPDF,"meanL1")->getErrorLo() << "   " << GetVar(TotalPDF,"meanL1")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("meanL1"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# Bkg sigma left-peak-1 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"sigmaL1") != NULL) && (vecConstr->find(string(string("sigmaL1") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("sigmaL1") << "   " << GetVar(TotalPDF,"sigmaL1")->getErrorLo() << "   " << GetVar(TotalPDF,"sigmaL1")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# Bkg mean left-peak-2 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"meanL2") != NULL) && (vecConstr->find(string(string("meanL2") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("meanL2") << "   " << GetVar(TotalPDF,"meanL2")->getErrorLo() << "   " << GetVar(TotalPDF,"meanL2")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("meanL2"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# Bkg sigma left-peak-2 [GeV/c2]");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"sigmaL2") != NULL) && (vecConstr->find(string(string("sigmaL2") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("sigmaL2") << "   " << GetVar(TotalPDF,"sigmaL2")->getErrorLo() << "   " << GetVar(TotalPDF,"sigmaL2")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# Fraction");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"fracMassBLPeak") != NULL) && (vecConstr->find(string(string("fracMassBLPeak") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("fracMassBLPeak") << "   " << GetVar(TotalPDF,"fracMassBLPeak")->getErrorLo() << "   " << GetVar(TotalPDF,"fracMassBLPeak")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("fracMassBLPeak"))->operator[](q2BinIndx).c_str());

  vecParStr->push_back("# Fraction of right-peak");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"fracMassBPeak") != NULL) && (vecConstr->find(string(string("fracMassBPeak") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("fracMassBPeak") << "   " << GetVar(TotalPDF,"fracMassBPeak")->getErrorLo() << "   " << GetVar(TotalPDF,"fracMassBPeak")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("fracMassBPeak"))->operator[](q2BinIndx).c_str());


  vecParStr->push_back("###################################################################");


  vecParStr->push_back("# Combinatorial bkg yield");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"nBkgComb") != NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("nBkgComb") << "   " << GetVar(TotalPDF,"nBkgComb")->getErrorLo() << "   " << GetVar(TotalPDF,"nBkgComb")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("nBkgComb"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# Mistag fraction");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"nMisTagFrac") != NULL) && ((atoi(Utility->GetGenericParam("SaveMisTagFrac").c_str()) == true) || (vecConstr->find(string(string("nMisTagFrac") + string("_constr")).c_str()) == NULL)))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("nMisTagFrac") << "   " << GetVar(TotalPDF,"nMisTagFrac")->getErrorLo() << "   " << GetVar(TotalPDF,"nMisTagFrac")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("nMisTagFrac"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# Peaking bkg yield");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"nBkgPeak") != NULL) && (vecConstr->find(string(string("nBkgPeak") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("nBkgPeak") << "   " << GetVar(TotalPDF,"nBkgPeak")->getErrorLo() << "   " << GetVar(TotalPDF,"nBkgPeak")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("nBkgPeak"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# Signal yield");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"nSig") != NULL) && (vecConstr->find(string(string("nSig") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("nSig") << "   " << GetVar(TotalPDF,"nSig")->getErrorLo() << "   " << GetVar(TotalPDF,"nSig")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("nSig"))->operator[](q2BinIndx).c_str());


  vecParStr->push_back("###################################################################");

  NCoeffPolyBKGpeak1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP1"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGcomb1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC1"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGpeak2 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP2"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGcomb2 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC2"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGpeak3 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP3"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGcomb3 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC3"))->operator[](q2BinIndx).c_str());

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
	  (vecConstr->find(string(myCoeff.str() + string("_constr")).c_str()) == NULL))
	{
	  myString << TotalPDF->getVariables()->getRealValue(myCoeff.str().c_str()) << "   " << GetVar(TotalPDF,myCoeff.str().c_str())->getErrorLo() << "   " << GetVar(TotalPDF,myCoeff.str().c_str())->getErrorHi();
	  vecParStr->push_back(myString.str());
	}
      else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](q2BinIndx).c_str());
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
	  (vecConstr->find(string(myCoeff.str() + string("_constr")).c_str()) == NULL))
	{
	  myString << TotalPDF->getVariables()->getRealValue(myCoeff.str().c_str()) << "   " << GetVar(TotalPDF,myCoeff.str().c_str())->getErrorLo() << "   " << GetVar(TotalPDF,myCoeff.str().c_str())->getErrorHi();
	  vecParStr->push_back(myString.str());
	}
      else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](q2BinIndx).c_str());
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

  vecParStr->push_back("### Number of coefficients of the peak.bkg angle3 poly ###");
  myString.clear(); myString.str("");
  myString << NCoeffPolyBKGpeak3;
  vecParStr->push_back(myString.str());
  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
    {
      myCoeff.clear(); myCoeff.str("");
      myCoeff << "p3Poly" << i;
      myString.clear(); myString.str("");
      myString << "# Bkg angle3 poly coeff. p" << i;
      vecParStr->push_back(myString.str());
      myString.clear(); myString.str("");
      if ((SAVEPOLY == true) &&
	  (TotalPDF != NULL) &&
	  (GetVar(TotalPDF,myCoeff.str().c_str()) != NULL) &&
	  (vecConstr->find(string(myCoeff.str() + string("_constr")).c_str()) == NULL))
	{
	  myString << TotalPDF->getVariables()->getRealValue(myCoeff.str().c_str()) << "   " << GetVar(TotalPDF,myCoeff.str().c_str())->getErrorLo() << "   " << GetVar(TotalPDF,myCoeff.str().c_str())->getErrorHi();
	  vecParStr->push_back(myString.str());
	}
      else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](q2BinIndx).c_str());
    }

  vecParStr->push_back("### Number of coefficients of the comb.bkg angle3 poly ###");
  myString.clear(); myString.str("");
  myString << NCoeffPolyBKGcomb3;
  vecParStr->push_back(myString.str());
  for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
    {
      myCoeff.clear(); myCoeff.str("");
      myCoeff << "c3Poly" << i;
      myString.clear(); myString.str("");
      myString << "# Bkg angle3 poly coeff. c" << i;
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
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"FlS") != NULL) && (vecConstr->find(string(string("FlS") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("FlS") << "   " << GetVar(TotalPDF,"FlS")->getErrorLo() << "   " << GetVar(TotalPDF,"FlS")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("FlS"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# AFB +/- err");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"AfbS") != NULL) && (vecConstr->find(string(string("AfbS") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("AfbS") << "   " << GetVar(TotalPDF,"AfbS")->getErrorLo() << "   " << GetVar(TotalPDF,"AfbS")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("AfbS"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# P1 +/- err");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"P1S") != NULL) && (vecConstr->find(string(string("P1S") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("P1S") << "   " << GetVar(TotalPDF,"P1S")->getErrorLo() << "   " << GetVar(TotalPDF,"P1S")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("P1S"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# P2 +/- err");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"P2S") != NULL) && (vecConstr->find(string(string("P2S") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("P2S") << "   " << GetVar(TotalPDF,"P2S")->getErrorLo() << "   " << GetVar(TotalPDF,"P2S")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("P2S"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# FS +/- err");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"FsS") != NULL) && (vecConstr->find(string(string("FsS") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("FsS") << "   " << GetVar(TotalPDF,"FsS")->getErrorLo() << "   " << GetVar(TotalPDF,"FsS")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("FsS"))->operator[](q2BinIndx).c_str());
  vecParStr->push_back("# AS +/- err");
  if ((TotalPDF != NULL) && (GetVar(TotalPDF,"AsS") != NULL) && (vecConstr->find(string(string("AsS") + string("_constr")).c_str()) == NULL))
    {
      myString.clear(); myString.str("");
      myString << TotalPDF->getVariables()->getRealValue("AsS") << "   " << GetVar(TotalPDF,"AsS")->getErrorLo() << "   " << GetVar(TotalPDF,"AsS")->getErrorHi();
      vecParStr->push_back(myString.str());
    }
  else vecParStr->push_back(fitParam->operator[](Utility->GetFitParamIndx("AsS"))->operator[](q2BinIndx).c_str());


  vecParStr->push_back("###################################################################");


  vecParStr->push_back("# Use signal: 0 = no; 1 = 1 gaussian; 2 = 2 gaussians");
  myString.clear(); myString.str("");
  myString << configParam->operator[](Utility->GetConfigParamIndx("SigType"))->operator[](q2BinIndx);
  vecParStr->push_back(myString.str());

  vecParStr->push_back("# Use peaking bkg: 0 = no; 1 = 1 gaussian; 2 = 2 gaussians; 11 = two 1 gaussian; 12 = two 2 gaussians");
  myString.clear(); myString.str("");
  myString << configParam->operator[](Utility->GetConfigParamIndx("PeakBkgType"))->operator[](q2BinIndx);
  vecParStr->push_back(myString.str());

  vecParStr->push_back("# Use combinatorial bkg: 0 = no; 1 = 1 exponential; 2 = 2 exponentials; 3 = erfc");
  myString.clear(); myString.str("");
  myString << configParam->operator[](Utility->GetConfigParamIndx("CombBkgType"))->operator[](q2BinIndx);
  vecParStr->push_back(myString.str());

  vecParStr->push_back("# Use mistag: 0 = no; 1 = 1 gaussian; 2 = 2 gaussians");
  myString.clear(); myString.str("");
  myString << configParam->operator[](Utility->GetConfigParamIndx("MistagType"))->operator[](q2BinIndx);
  vecParStr->push_back(myString.str());

  return vecParStr;
}


unsigned int CopyFitResults (RooAbsPdf* TotalPDF, unsigned int q2BinIndx, vector<vector<string>*>* fitParam, unsigned int countMisTag, unsigned int countGoodTag)
// ####################################################################
// # Only for polynomial coeficients:                                 #
// # If errLo and errHi are 0.0 --> set poly. coefficient as constant #
// ####################################################################
{
  stringstream myString;
  stringstream myCoeff;
  double value, errLo, errHi;
  unsigned int NCoeffPolyBKGpeak1;
  unsigned int NCoeffPolyBKGcomb1;
  unsigned int NCoeffPolyBKGpeak2;
  unsigned int NCoeffPolyBKGcomb2;
  unsigned int NCoeffPolyBKGpeak3;
  unsigned int NCoeffPolyBKGcomb3;


  if (GetVar(TotalPDF,"meanS") != NULL)
    {
      TotalPDF->getVariables()->setRealValue("meanS",atof(fitParam->operator[](Utility->GetFitParamIndx("meanS"))->operator[](q2BinIndx).c_str()));
      GetVar(TotalPDF,"meanS")->setConstant(false);
    }
  if (GetVar(TotalPDF,"sigmaS1") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"sigmaS1",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"sigmaS1")->setConstant(true);
      else                                  GetVar(TotalPDF,"sigmaS1")->setConstant(false);
    }
  if (GetVar(TotalPDF,"sigmaS2") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaS2"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"sigmaS2",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"sigmaS2")->setConstant(true);
      else                                  GetVar(TotalPDF,"sigmaS2")->setConstant(false);
    }
  if (GetVar(TotalPDF,"fracMassS") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("fracMassS"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"fracMassS",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"fracMassS")->setConstant(true);
      else                                  GetVar(TotalPDF,"fracMassS")->setConstant(false);
    }


  if (GetVar(TotalPDF,"var1") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("var1"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"var1",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"var1")->setConstant(true);
      else                                  GetVar(TotalPDF,"var1")->setConstant(false);
    }
  if (GetVar(TotalPDF,"var2") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("var2"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"var2",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"var2")->setConstant(true);
      else                                  GetVar(TotalPDF,"var2")->setConstant(false);
    }
  if (GetVar(TotalPDF,"fracMassBExp") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("fracMassBExp"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"fracMassBExp",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"fracMassBExp")->setConstant(true);
      else                                  GetVar(TotalPDF,"fracMassBExp")->setConstant(false);
    }


  if (GetVar(TotalPDF,"sigmaMisTag1") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag1"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"sigmaMisTag1",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"sigmaMisTag1")->setConstant(true);
      else                                  GetVar(TotalPDF,"sigmaMisTag1")->setConstant(false);
    }
  if (GetVar(TotalPDF,"sigmaMisTag2") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag2"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"sigmaMisTag2",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"sigmaMisTag2")->setConstant(true);
      else                                  GetVar(TotalPDF,"sigmaMisTag2")->setConstant(false);
    }
  if (GetVar(TotalPDF,"fracMisTag") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("fracMisTag"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"fracMisTag",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"fracMisTag")->setConstant(true);
      else                                  GetVar(TotalPDF,"fracMisTag")->setConstant(false);
    }


  if (GetVar(TotalPDF,"meanR1") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("meanR1"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"meanR1",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"meanR1")->setConstant(true);
      else                                  GetVar(TotalPDF,"meanR1")->setConstant(false);
    }
  if (GetVar(TotalPDF,"sigmaR1") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"sigmaR1",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"sigmaR1")->setConstant(true);
      else                                  GetVar(TotalPDF,"sigmaR1")->setConstant(false);
    }
  if (GetVar(TotalPDF,"meanR2") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("meanR2"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"meanR2",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"meandR2")->setConstant(true);
      else                                  GetVar(TotalPDF,"meandR2")->setConstant(false);
    }
  if (GetVar(TotalPDF,"sigmaR2") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"sigmaR2",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"sigmaR2")->setConstant(true);
      else                                  GetVar(TotalPDF,"sigmaR2")->setConstant(false);
    }
  if (GetVar(TotalPDF,"fracMassBRPeak") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("fracMassBRPeak"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"fracMassBRPeak",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"fracMassBRPeak")->setConstant(true);
      else                                  GetVar(TotalPDF,"fracMassBRPeak")->setConstant(false);
    }
  if (GetVar(TotalPDF,"meanL1") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("meanL1"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"meanL1",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"meanL1")->setConstant(true);
      else                                  GetVar(TotalPDF,"meanL1")->setConstant(false);
    }
  if (GetVar(TotalPDF,"sigmaL1") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"sigmaL1",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"sigmaL1")->setConstant(true);
      else                                  GetVar(TotalPDF,"sigmaL1")->setConstant(false);
    }
  if (GetVar(TotalPDF,"meanL2") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("meanL2"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"meanL2",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"meanL2")->setConstant(true);
      else                                  GetVar(TotalPDF,"meanL2")->setConstant(false);
    }
  if (GetVar(TotalPDF,"sigmaL2") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"sigmaL2",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"sigmaL2")->setConstant(true);
      else                                  GetVar(TotalPDF,"sigmaL2")->setConstant(false);
    }
  if (GetVar(TotalPDF,"fracMassBLPeak") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("fracMassBLPeak"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"fracMassBLPeak",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"fracMassBLPeak")->setConstant(true);
      else                                  GetVar(TotalPDF,"fracMassBLPeak")->setConstant(false);
    }
  if (GetVar(TotalPDF,"fracMassBPeak") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("fracMassBPeak"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"fracMassBPeak",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"fracMassBPeak")->setConstant(true);
      else                                  GetVar(TotalPDF,"fracMassBPeak")->setConstant(false);
    }


  if (GetVar(TotalPDF,"nBkgComb") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("nBkgComb"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"nBkgComb",MULTYIELD,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"nBkgComb")->setConstant(true);
      else                                  GetVar(TotalPDF,"nBkgComb")->setConstant(false);
    }
  if (GetVar(TotalPDF,"nMisTagFrac") != NULL)
    {
      if (atoi(Utility->GetGenericParam("CtrlMisTagWrkFlow").c_str()) != 3)
	{
	  myString.clear(); myString.str("");
	  myString << fitParam->operator[](Utility->GetFitParamIndx("nMisTagFrac"))->operator[](q2BinIndx).c_str();
	  SetValueAndErrors(TotalPDF,"nMisTagFrac",1.0,&myString,&value,&errLo,&errHi);
	  if ((atoi(Utility->GetGenericParam("CtrlMisTagWrkFlow").c_str()) == 0) || ((errLo == 0.0) && (errHi == 0.0))) GetVar(TotalPDF,"nMisTagFrac")->setConstant(true);
	  else                                                                                                          GetVar(TotalPDF,"nMisTagFrac")->setConstant(false);
	}
      else
	{
	  myString.clear(); myString.str("");
	  myString << static_cast<double>(countMisTag) / static_cast<double>(countMisTag + countGoodTag);
	  SetValueAndErrors(TotalPDF,"nMisTagFrac",1.0,&myString,&value,&errLo,&errHi);
	  GetVar(TotalPDF,"nMisTagFrac")->setConstant(true);
	}
    }
  if (GetVar(TotalPDF,"nBkgPeak") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("nBkgPeak"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"nBkgPeak",MULTYIELD,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"nBkgPeak")->setConstant(true);
      else                                  GetVar(TotalPDF,"nBkgPeak")->setConstant(false);
    }
  if (GetVar(TotalPDF,"nSig") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("nSig"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"nSig",MULTYIELD,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"nSig")->setConstant(true);
      else                                  GetVar(TotalPDF,"nSig")->setConstant(false);
    }


  NCoeffPolyBKGpeak1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP1"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGcomb1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC1"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGpeak2 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP2"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGcomb2 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC2"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGpeak3 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP3"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGcomb3 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC3"))->operator[](q2BinIndx).c_str());
  if ((NCoeffPolyBKGpeak1 > NCOEFFPOLYBKG) || (NCoeffPolyBKGcomb1 > NCOEFFPOLYBKG) ||
      (NCoeffPolyBKGpeak2 > NCOEFFPOLYBKG) || (NCoeffPolyBKGcomb2 > NCOEFFPOLYBKG) ||
      (NCoeffPolyBKGpeak3 > NCOEFFPOLYBKG) || (NCoeffPolyBKGcomb3 > NCOEFFPOLYBKG))
    {
      cout << "[ExtractYield::CopyFitResults]\tDegree of poly bkg is not within allowed limits: ";
      cout << NCoeffPolyBKGpeak1 << "\t" << NCoeffPolyBKGcomb1 << "\t" << NCoeffPolyBKGpeak2 << "\t" << NCoeffPolyBKGcomb2 << "\t" << NCoeffPolyBKGpeak3 << "\t" << NCoeffPolyBKGcomb3 << endl;
      exit (EXIT_FAILURE);
    }


  for (unsigned int i = 0; i < NCoeffPolyBKGpeak1; i++)
    {
      myCoeff.clear(); myCoeff.str("");
      myCoeff << "p1Poly" << i;
      if (GetVar(TotalPDF,myCoeff.str().c_str()) != NULL)
	{
	  myString.clear(); myString.str("");
	  myString << fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](q2BinIndx).c_str();
	  SetValueAndErrors(TotalPDF,myCoeff.str(),1.0,&myString,&value,&errLo,&errHi);
	  if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,myCoeff.str().c_str())->setConstant(true);
	  else                                  GetVar(TotalPDF,myCoeff.str().c_str())->setConstant(false);
	}
    }
  for (unsigned int i = 0; i < NCoeffPolyBKGcomb1; i++)
    {
      myCoeff.clear(); myCoeff.str("");
      myCoeff << "c1Poly" << i;
      if (GetVar(TotalPDF,myCoeff.str().c_str()) != NULL)
	{
	  myString.clear(); myString.str("");
	  myString << fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](q2BinIndx).c_str();
	  SetValueAndErrors(TotalPDF,myCoeff.str(),1.0,&myString,&value,&errLo,&errHi);
	  if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,myCoeff.str().c_str())->setConstant(true);
	  else                                  GetVar(TotalPDF,myCoeff.str().c_str())->setConstant(false);
	}
    }
  for (unsigned int i = 0; i < NCoeffPolyBKGpeak2; i++)
    {
      myCoeff.clear(); myCoeff.str("");
      myCoeff << "p2Poly" << i;
      if (GetVar(TotalPDF,myCoeff.str().c_str()) != NULL)
	{
	  myString.clear(); myString.str("");
	  myString << fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](q2BinIndx).c_str();
	  SetValueAndErrors(TotalPDF,myCoeff.str(),1.0,&myString,&value,&errLo,&errHi);
	  if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,myCoeff.str().c_str())->setConstant(true);
	  else                                  GetVar(TotalPDF,myCoeff.str().c_str())->setConstant(false);
	}
    }
  for (unsigned int i = 0; i < NCoeffPolyBKGcomb2; i++)
    {
      myCoeff.clear(); myCoeff.str("");
      myCoeff << "c2Poly" << i;
      if (GetVar(TotalPDF,myCoeff.str().c_str()) != NULL)
	{
	  myString.clear(); myString.str("");
	  myString << fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](q2BinIndx).c_str();
	  SetValueAndErrors(TotalPDF,myCoeff.str(),1.0,&myString,&value,&errLo,&errHi);
	  if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,myCoeff.str().c_str())->setConstant(true);
	  else                                  GetVar(TotalPDF,myCoeff.str().c_str())->setConstant(false);
	}
    }
  for (unsigned int i = 0; i < NCoeffPolyBKGpeak3; i++)
    {
      myCoeff.clear(); myCoeff.str("");
      myCoeff << "p3Poly" << i;
      if (GetVar(TotalPDF,myCoeff.str().c_str()) != NULL)
	{
	  myString.clear(); myString.str("");
	  myString << fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](q2BinIndx).c_str();
	  SetValueAndErrors(TotalPDF,myCoeff.str(),1.0,&myString,&value,&errLo,&errHi);
	  if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,myCoeff.str().c_str())->setConstant(true);
	  else                                  GetVar(TotalPDF,myCoeff.str().c_str())->setConstant(false);
	}
    }
  for (unsigned int i = 0; i < NCoeffPolyBKGcomb3; i++)
    {
      myCoeff.clear(); myCoeff.str("");
      myCoeff << "c3Poly" << i;
      if (GetVar(TotalPDF,myCoeff.str().c_str()) != NULL)
	{
	  myString.clear(); myString.str("");
	  myString << fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](q2BinIndx).c_str();
	  SetValueAndErrors(TotalPDF,myCoeff.str(),1.0,&myString,&value,&errLo,&errHi);
	  if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,myCoeff.str().c_str())->setConstant(true);
	  else                                  GetVar(TotalPDF,myCoeff.str().c_str())->setConstant(false);
	}
    }


  if (GetVar(TotalPDF,"FlS") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("FlS"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"FlS",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"FlS")->setConstant(true);
      else                                  GetVar(TotalPDF,"FlS")->setConstant(false);
    }
  if (GetVar(TotalPDF,"AfbS") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("AfbS"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"AfbS",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"AfbS")->setConstant(true);
      else                                  GetVar(TotalPDF,"AfbS")->setConstant(false);
    }
  if (GetVar(TotalPDF,"P1S") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("P1S"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"P1S",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"P1S")->setConstant(true);
      else                                  GetVar(TotalPDF,"P1S")->setConstant(false);
    }
  if (GetVar(TotalPDF,"P2S") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("P2S"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(TotalPDF,"P2S",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"P2S")->setConstant(true);
      else                                  GetVar(TotalPDF,"P2S")->setConstant(false);
    }
  if (GetVar(TotalPDF,"FsS") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("FsS"))->operator[](q2BinIndx);
      SetValueAndErrors(TotalPDF,"FsS",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"FsS")->setConstant(true);
      else                                  GetVar(TotalPDF,"FsS")->setConstant(false);
    }
  if (GetVar(TotalPDF,"AsS") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("AsS"))->operator[](q2BinIndx);
      SetValueAndErrors(TotalPDF,"AsS",1.0,&myString,&value,&errLo,&errHi);
      if ((errLo == 0.0) && (errHi == 0.0)) GetVar(TotalPDF,"AsS")->setConstant(true);
      else                                  GetVar(TotalPDF,"AsS")->setConstant(false);
    }


  value = 0.0;
  if (GetVar(TotalPDF,"nSig")        != NULL) value = value + TotalPDF->getVariables()->getRealValue("nSig");
  if (GetVar(TotalPDF,"nBkgComb")    != NULL) value = value + TotalPDF->getVariables()->getRealValue("nBkgComb");
  if (GetVar(TotalPDF,"nMisTagFrac") != NULL) value = value + TotalPDF->getVariables()->getRealValue("nSig") * TotalPDF->getVariables()->getRealValue("nMisTagFrac") / (1.0 - TotalPDF->getVariables()->getRealValue("nMisTagFrac"));
  if (GetVar(TotalPDF,"nBkgPeak")    != NULL) value = value + TotalPDF->getVariables()->getRealValue("nBkgPeak");
  return value;
}


void GenerateFitParameters (RooAbsPdf* TotalPDF, vector<vector<string>*>* fitParam, unsigned int fileIndx, vector<double>* q2Bins, unsigned int q2BinIndx, string option)
// #########################
// # option = "All"        #
// # option = "misTagFrac" #
// # option = "FlAfb"      #
// # option = "bkgAng"     #
// #########################
{
  unsigned int NCoeffPolyBKGcomb1;
  unsigned int NCoeffPolyBKGcomb2;
  unsigned int NCoeffPolyBKGcomb3;

  stringstream myString;

  CopyFitResults(TotalPDF,q2BinIndx,fitParam);

  RooRandom::randomGenerator()->SetSeed(fileIndx*(q2Bins->size()-1) + q2BinIndx + 1);
  cout << "\n[ExtractYield::GenerateFitParameters]\t@@@ Random seed for parameter file generation set to : " << RooRandom::randomGenerator()->GetSeed() << " @@@" << endl;


  if ((option == "All") || (option == "misTag"))
    {
      if ((atoi(Utility->GetGenericParam("CtrlMisTagWrkFlow").c_str()) == 0) && (GetVar(TotalPDF,"nMisTagFrac") != NULL))
	{
	  cout << "Mis-tag fraction generation: gaussian mean = " << GetVar(TotalPDF,"nMisTagFrac")->getVal() << "\tsigma = " << GetVar(TotalPDF,"nMisTagFrac")->getError() << endl;
	  TotalPDF->getVariables()->setRealValue("nMisTagFrac",RooRandom::gaussian() * GetVar(TotalPDF,"nMisTagFrac")->getError() + GetVar(TotalPDF,"nMisTagFrac")->getVal());
	  GetVar(TotalPDF,"nMisTagFrac")->setAsymError(-1.0,1.0);
	  GetVar(TotalPDF,"nMisTagFrac")->setError(1.0);
	}
    }


  if (option == "All")
    {
      if (GetVar(TotalPDF,"var1") != NULL)
	{
	  cout << "Background var1 generation: uniform lower bound = 0\thigher bound = 1" << endl;
	  TotalPDF->getVariables()->setRealValue("var1",RooRandom::uniform());
	  GetVar(TotalPDF,"var1")->setAsymError(-1.0,1.0);
	  GetVar(TotalPDF,"var1")->setError(1.0);
	}

      if (GetVar(TotalPDF,"var2") != NULL)
	{
	  cout << "Background var2 generation: uniform lower bound = 0\thigher bound = 1" << endl;
	  TotalPDF->getVariables()->setRealValue("var2",RooRandom::uniform());
	  GetVar(TotalPDF,"var2")->setAsymError(-1.0,1.0);
	  GetVar(TotalPDF,"var2")->setError(1.0);
	}

      if (GetVar(TotalPDF,"fracMassBExp") != NULL)
	{
	  cout << "Background fraction uniform generation" << endl;
	  TotalPDF->getVariables()->setRealValue("fracMassBExp",RooRandom::uniform());
	  GetVar(TotalPDF,"fracMassBExp")->setAsymError(-1.0,1.0);
	  GetVar(TotalPDF,"fracMassBExp")->setError(1.0);
	}
    }


  if ((option == "All") || (option == "FlAfb"))
    {
      if (GetVar(TotalPDF,"FlS") != NULL)
	{
	  cout << "Fl generation: uniform lower bound = " << GetVar(TotalPDF,"FlS")->getMin() << "\thigher bound = " << GetVar(TotalPDF,"FlS")->getMax() << endl;
	  TotalPDF->getVariables()->setRealValue("FlS",RooRandom::uniform() * (GetVar(TotalPDF,"FlS")->getMax() - GetVar(TotalPDF,"FlS")->getMin()) + GetVar(TotalPDF,"FlS")->getMin());
	  GetVar(TotalPDF,"FlS")->setAsymError(-1.0,1.0);
	  GetVar(TotalPDF,"FlS")->setError(1.0);
	}
      
      if (GetVar(TotalPDF,"AfbS") != NULL)
	{
	  cout << "Afb generation: uniform lower bound = " << GetVar(TotalPDF,"AfbS")->getMin() << "\thigher bound = " << GetVar(TotalPDF,"AfbS")->getMax() << endl;
	  TotalPDF->getVariables()->setRealValue("AfbS",RooRandom::uniform() * (GetVar(TotalPDF,"AfbS")->getMax() - GetVar(TotalPDF,"AfbS")->getMin()) + GetVar(TotalPDF,"AfbS")->getMin());
	  GetVar(TotalPDF,"AfbS")->setAsymError(-1.0,1.0);
	  GetVar(TotalPDF,"AfbS")->setError(1.0);
	}
    }


  if (option == "All")
    {
      if (GetVar(TotalPDF,"P1S") != NULL)
	{
	  cout << "P1 generation: uniform lower bound = " << GetVar(TotalPDF,"P1S")->getMin() << "\thigher bound = " << GetVar(TotalPDF,"P1S")->getMax() << endl;
	  TotalPDF->getVariables()->setRealValue("P1S",RooRandom::uniform() * (GetVar(TotalPDF,"P1S")->getMax() - GetVar(TotalPDF,"P1S")->getMin()) + GetVar(TotalPDF,"P1S")->getMin());
	  GetVar(TotalPDF,"P1S")->setAsymError(-1.0,1.0);
	  GetVar(TotalPDF,"P1S")->setError(1.0);
	}

      if (GetVar(TotalPDF,"P2S") != NULL)
	{
	  cout << "P2 generation: uniform lower bound = " << GetVar(TotalPDF,"P2S")->getMin() << "\thigher bound = " << GetVar(TotalPDF,"P2S")->getMax() << endl;
	  TotalPDF->getVariables()->setRealValue("P2S",RooRandom::uniform() * (GetVar(TotalPDF,"P2S")->getMax() - GetVar(TotalPDF,"P2S")->getMin()) + GetVar(TotalPDF,"P2S")->getMin());
	  GetVar(TotalPDF,"P2S")->setAsymError(-1.0,1.0);
	  GetVar(TotalPDF,"P2S")->setError(1.0);
	}

      if (GetVar(TotalPDF,"FsS") != NULL)
	{
	  cout << "Fs generation: uniform lower bound = " << GetVar(TotalPDF,"FsS")->getMin() << "\thigher bound = " << GetVar(TotalPDF,"FsS")->getMax() << endl;
	  TotalPDF->getVariables()->setRealValue("FsS",RooRandom::uniform() * (GetVar(TotalPDF,"FsS")->getMax() - GetVar(TotalPDF,"FsS")->getMin()) + GetVar(TotalPDF,"FsS")->getMin());
	  GetVar(TotalPDF,"FsS")->setAsymError(-1.0,1.0);
	  GetVar(TotalPDF,"FsS")->setError(1.0);
	}

      if (GetVar(TotalPDF,"AsS") != NULL)
	{
	  cout << "As generation: uniform lower bound = " << GetVar(TotalPDF,"AsS")->getMin() << "\thigher bound = " << GetVar(TotalPDF,"AsS")->getMax() << endl;
	  TotalPDF->getVariables()->setRealValue("AsS",RooRandom::uniform() * (GetVar(TotalPDF,"AsS")->getMax() - GetVar(TotalPDF,"AsS")->getMin()) + GetVar(TotalPDF,"AsS")->getMin());
	  GetVar(TotalPDF,"AsS")->setAsymError(-1.0,1.0);
	  GetVar(TotalPDF,"AsS")->setError(1.0);
	}
    }


  NCoeffPolyBKGcomb1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC1"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGcomb2 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC2"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGcomb3 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC3"))->operator[](q2BinIndx).c_str());

  if ((option == "All") || (option == "bkgAng"))
    {
      for (unsigned int i = 0; i < NCoeffPolyBKGcomb1; i++)
	{
	  myString.clear(); myString.str("");
	  myString << "c1Poly" << i;
	  if (GetVar(TotalPDF,myString.str().c_str()) != NULL)
	    {
	      cout << myString.str().c_str() << " generation: gaussian mean = " << GetVar(TotalPDF,myString.str().c_str())->getVal() << "\tsigma = " << GetVar(TotalPDF,myString.str().c_str())->getError() << endl;
	      TotalPDF->getVariables()->setRealValue(myString.str().c_str(),RooRandom::gaussian() * GetVar(TotalPDF,myString.str().c_str())->getError() + GetVar(TotalPDF,myString.str().c_str())->getVal());
	      GetVar(TotalPDF,myString.str().c_str())->setAsymError(-1.0,1.0);
	      GetVar(TotalPDF,myString.str().c_str())->setError(1.0);
	    }
	}
      for (unsigned int i = 0; i < NCoeffPolyBKGcomb2; i++)
	{
	  myString.clear(); myString.str("");
	  myString << "c2Poly" << i;
	  if (GetVar(TotalPDF,myString.str().c_str()) != NULL)
	    {
	      cout << myString.str().c_str() << " generation: gaussian mean = " << GetVar(TotalPDF,myString.str().c_str())->getVal() << "\tsigma = " << GetVar(TotalPDF,myString.str().c_str())->getError() << endl;
	      TotalPDF->getVariables()->setRealValue(myString.str().c_str(),RooRandom::gaussian() * GetVar(TotalPDF,myString.str().c_str())->getError() + GetVar(TotalPDF,myString.str().c_str())->getVal());
	      GetVar(TotalPDF,myString.str().c_str())->setAsymError(-1.0,1.0);
	      GetVar(TotalPDF,myString.str().c_str())->setError(1.0);
	    }
	}
      for (unsigned int i = 0; i < NCoeffPolyBKGcomb3; i++)
	{
	  myString.clear(); myString.str("");
	  myString << "c3Poly" << i;
	  if (GetVar(TotalPDF,myString.str().c_str()) != NULL)
	    {
	      cout << myString.str().c_str() << " generation: gaussian mean = " << GetVar(TotalPDF,myString.str().c_str())->getVal() << "\tsigma = " << GetVar(TotalPDF,myString.str().c_str())->getError() << endl;
	      TotalPDF->getVariables()->setRealValue(myString.str().c_str(),RooRandom::gaussian() * GetVar(TotalPDF,myString.str().c_str())->getError() + GetVar(TotalPDF,myString.str().c_str())->getVal());
	      GetVar(TotalPDF,myString.str().c_str())->setAsymError(-1.0,1.0);
	      GetVar(TotalPDF,myString.str().c_str())->setError(1.0);
	    }
	}
    }
}


void GenerateDataset (RooAbsPdf* TotalPDF, RooArgSet setVar, vector<double>* q2Bins, int q2BinIndx, vector<vector<string>*>* fitParam, string fileName)
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
  NtplFileOut->mkdir("B0KstMuMu");
  NtplFileOut->cd("B0KstMuMu");
  theTreeOut = new TTree("B0KstMuMuNTuple","B0KstMuMuNTuple");
  NTupleOut  = new B0KstMuMuSingleCandTreeContent();
  NTupleOut->Init();
  NTupleOut->ClearNTuple();
  NTupleOut->MakeTreeBranches(theTreeOut);


  cout << "\n[ExtractYield::GenerateDataset]\t@@@ Generating data-like ntuple @@@" << endl;
  RooRandom::randomGenerator()->SetSeed(q2BinIndx + 1);
  cout << "[ExtractYield::GenerateDataset]\t@@@ Random seed for parameter file generation set to : " << RooRandom::randomGenerator()->GetSeed() << " @@@" << endl;


  // #####################
  // # Toy-MC generation #
  // #####################
  nEntryToy = CopyFitResults(TotalPDF,q2BinIndx,fitParam);
  PrintVariables(TotalPDF->getVariables(),"vars");
  MyToy = new RooMCStudy(*TotalPDF,setVar);
  MyToy->generate(1,nEntryToy,true);


  // ########################
  // # Loop over the events #
  // ########################
  toySample = (RooDataSet*)MyToy->genData(0);
  for (unsigned int entry = 0; entry < nEntryToy; entry++)
    {
      if (q2BinIndx == Utility->GetJPsiBin(q2Bins))
	{
	  NTupleOut->mumuMass->push_back(Utility->JPsiMass);
	  NTupleOut->mumuMassE->push_back(1e-3);
	}
      else if (q2BinIndx == Utility->GetPsiPBin(q2Bins))
	{
	  NTupleOut->mumuMass->push_back(Utility->PsiPMass);
	  NTupleOut->mumuMassE->push_back(1e-3);
	}
      else
	{
	  NTupleOut->mumuMass->push_back(sqrt(q2Bins->operator[](q2BinIndx) + RooRandom::uniform() * (q2Bins->operator[](q2BinIndx+1) - q2Bins->operator[](q2BinIndx))));
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
      NTupleOut->rightFlavorTag = (RooRandom::uniform() < 0.5 ? false : true);
      NTupleOut->TrigCat        = 1;

      theTreeOut->Fill();
      NTupleOut->ClearNTuple();
    }


  // #############
  // # Save data #
  // #############
  NtplFileOut->cd("B0KstMuMu");
  theTreeOut->Write();
  NtplFileOut->Close("R");
  delete NTupleOut;
  delete NtplFileOut;
  delete MyToy;
}


string GeneratePolynomial (RooRealVar* var, unsigned int nCoef, string sCoef)
{
  stringstream myString;

  myString.clear(); myString.str("");
  myString << "1";
  for (unsigned int i = 0; i < nCoef; i++)
    {
      myString << " + ";
      for (unsigned int j = 0; j < i+1; j++) myString << var->getPlotLabel() << "*";
      myString << sCoef << i;
    }
  cout << "[ExtractYield::GeneratePolynomial]\tI've generated the polynomial: " << myString.str().c_str() << endl;
  
  return myString.str();
}


void FitDimuonInvMass (RooDataSet* dataSet, RooAbsPdf** TotalPDFJPsi, RooAbsPdf** TotalPDFPsiP, RooRealVar* x, TCanvas* Canv, bool justPlotMuMuMass, bool justKeepPsi, string plotName)
{
  stringstream myString;
  double sigmaJPsi, sigmaJPsiE;
  double sigmaPsiP, sigmaPsiPE;
  unsigned int NBins;
  int nElements = 4;
  if (justKeepPsi == true) nElements = 3;
  TString legNamesPsi[4];
  legNamesPsi[0] = "Data";
  legNamesPsi[1] = "Total fit";
  legNamesPsi[2] = "Signal";
  legNamesPsi[3] = "Comb. bkg";


  if (justPlotMuMuMass == true)
    {
      myString.clear(); myString.str("");
      myString << "(mumuMass > " << 0.0 << " && mumuMass < " << 4.5 << ")";
      cout << "[ExtractYield::FitDimuonInvMass]\tCut for dimuon invariant mass: " << myString.str().c_str() << endl;
      RooDataSet* dataSetJPsi = (RooDataSet*)dataSet->reduce(myString.str().c_str());
      x->setRange("subRangeJPsi",0.0,4.5);
 
      Canv->cd();
      RooPlot* myFrameJPsi = x->frame(Range("subRangeJPsi"),Bins(25));
      dataSetJPsi->plotOn(myFrameJPsi,Name(dataSetJPsi->GetName()));
      myFrameJPsi->Draw();

      Canv->Modified();
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
      cout << "[ExtractYield::FitDimuonInvMass]\tCut for dimuon invariant mass: " << myString.str().c_str() << endl;
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

      RooAbsPdf* JPsi = new RooAddPdf("JPsi","J/psi pdf",RooArgSet(*MassJPsi1,*MassJPsi2),RooArgSet(*fracMassJPsi));


      // #####################################################
      // # Define fit variables and pdf for Background J/psi #
      // #####################################################
      RooRealVar* a0JPsi = new RooRealVar("a0JPsi","First coefficient bkg poly.",0.0);
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
      if (justKeepPsi == true) *TotalPDFJPsi = new RooAddPdf(myString.str().c_str(),"Total J/psi mass extended pdf",RooArgSet(*JPsi),RooArgSet(*nJPsiMass));
      else                     *TotalPDFJPsi = new RooAddPdf(myString.str().c_str(),"Total J/psi mass extended pdf",RooArgSet(*JPsi,*BkgmumuMassJPsi),RooArgSet(*nJPsiMass,*nBkgMassJPsi));


      // ###################
      // # Make actual fit #
      // ###################
      RooFitResult* JPsiFitResult;
      JPsiFitResult = (*TotalPDFJPsi)->fitTo(*dataSetJPsi,Extended(true),Save(true),SumCoefRange("subRangeJPsi"),Range("subRangeJPsi"));


      // ##################################################
      // # Set p.d.f independent variables to known point #
      // ##################################################
      if (GetVar(*TotalPDFJPsi,x->getPlotLabel()) != NULL) (*TotalPDFJPsi)->getVariables()->setRealValue(x->getPlotLabel(),Utility->B0Mass);
      if (JPsiFitResult != NULL) JPsiFitResult->Print("v");
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
      (*TotalPDFJPsi)->plotOn(myFrameJPsi, Components(*JPsi), LineStyle(7), LineColor(kBlue));
      if (justKeepPsi != true) (*TotalPDFJPsi)->plotOn(myFrameJPsi, Components(*BkgmumuMassJPsi), LineStyle(4), LineColor(kRed));
      (*TotalPDFJPsi)->paramOn(myFrameJPsi,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(RooArgSet(*nJPsiMass,*nBkgMassJPsi)));

      myString.clear(); myString.str("");
      myString << (*TotalPDFJPsi)->getPlotLabel() << "_paramBox";
      TPaveText* paveTextJPsi = (TPaveText*)myFrameJPsi->findObject(myString.str().c_str());
      paveTextJPsi->AddText(Form("%s%.3f#pm%.3f","#mu#kern[-0.4]{#lower[-0.3]{_{1}}} = ",GetVar(*TotalPDFJPsi,"meanJPsi1")->getVal(),GetVar(*TotalPDFJPsi,"meanJPsi1")->getError()));
      paveTextJPsi->AddText(Form("%s%.3f#pm%.3f","#mu#kern[-0.4]{#lower[-0.3]{_{2}}} = ",GetVar(*TotalPDFJPsi,"meanJPsi2")->getVal(),GetVar(*TotalPDFJPsi,"meanJPsi2")->getError()));
      paveTextJPsi->AddText(Form("%s%.3f#pm%.3f","< #sigma > = ",sigmaJPsi,sigmaJPsiE));
      paveTextJPsi->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameJPsi->chiSquare((*TotalPDFJPsi)->getPlotLabel(),dataSetJPsi->GetName())));
      CheckGoodFit(JPsiFitResult,paveTextJPsi);
      paveTextJPsi->SetBorderSize(0.0);
      paveTextJPsi->SetFillStyle(0);
      paveTextJPsi->SetTextSize(0.03);
      DrawString(LUMI,myFrameJPsi);
      myFrameJPsi->Draw();

      TLegend* legJPsi = new TLegend(0.75, 0.65, 0.97, 0.88, "");
      for (int i = 0; i < nElements; i++)
      	{
      	  TString objName = myFrameJPsi->nameOf(i);
      	  if (objName == "") continue;
      	  TObject* obj = myFrameJPsi->findObject(objName.Data());
      	  legJPsi->AddEntry(obj,legNamesPsi[i],"lp");
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
      cout << "[ExtractYield::FitDimuonInvMass]\tCut for dimuon invariant mass: " << myString.str().c_str() << endl;
      RooDataSet* dataSetPsiP = (RooDataSet*)dataSet->reduce(myString.str().c_str());
      x->setRange("subRangePsiP",3.3,4.1);
      NBins = 40;


      // ############################################
      // # Define fit variables and pdf for psi(2S) #
      // ############################################
      RooRealVar* meanPsiP1  = new RooRealVar("meanPsiP1","psi(2S) mean Gaussian-1",Utility->PsiPMass,"GeV");
      RooRealVar* sigmaPsiP1 = new RooRealVar("sigmaPsiP1","psi(2S) sigma Gaussian-1",0.02,"GeV");
      RooAbsPdf* MassPsiP1   = new RooGaussian("MassPsiP1","psi(2S) Gaussian-1",*x,*meanPsiP1,*sigmaPsiP1);

      RooRealVar* meanPsiP2  = new RooRealVar("meanPsiP2","psi(2S) mean Gaussian-2",Utility->PsiPMass,"GeV");
      RooRealVar* sigmaPsiP2 = new RooRealVar("sigmaPsiP2","psi(2S) sigma Gaussian-2",0.09,"GeV");
      RooAbsPdf* MassPsiP2   = new RooGaussian("MassPsiP2","psi(2S) Gaussian-2",*x,*meanPsiP2,*sigmaPsiP2);

      meanPsiP1->setConstant(false);
      meanPsiP2->setConstant(false);
      sigmaPsiP1->setConstant(false);
      sigmaPsiP2->setConstant(false);

      RooRealVar* fracMassPsiP = new RooRealVar("fracMassPsiP","Fraction of psi(2S) Gaussians",0.5,0.0,1.0);
      fracMassPsiP->setConstant(false);

      RooAbsPdf* PsiP = new RooAddPdf("PsiP","psi(2S) pdf",RooArgSet(*MassPsiP1,*MassPsiP2),RooArgSet(*fracMassPsiP));


      // #######################################################
      // # Define fit variables and pdf for Background psi(2S) #
      // #######################################################
      RooRealVar* a0PsiP = new RooRealVar("a0PsiP","First coefficient bkg poly.",0.0);
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
      if (justKeepPsi == true) *TotalPDFPsiP = new RooAddPdf(myString.str().c_str(),"Total psi(2S) mass extended pdf",RooArgSet(*PsiP),RooArgSet(*nPsiPMass));
      else                     *TotalPDFPsiP = new RooAddPdf(myString.str().c_str(),"Total psi(2S) mass extended pdf",RooArgSet(*PsiP,*BkgmumuMassPsiP),RooArgSet(*nPsiPMass,*nBkgMassPsiP));


      // ###################
      // # Make actual fit #
      // ###################
      RooFitResult* PsiPFitResult;
      PsiPFitResult = (*TotalPDFPsiP)->fitTo(*dataSetPsiP,Extended(true),Save(true),SumCoefRange("subRangePsiP"),Range("subRangePsiP"));


      // ##################################################
      // # Set p.d.f independent variables to known point #
      // ##################################################
      if (GetVar(*TotalPDFPsiP,x->getPlotLabel()) != NULL) (*TotalPDFPsiP)->getVariables()->setRealValue(x->getPlotLabel(),Utility->B0Mass);
      if (PsiPFitResult != NULL) PsiPFitResult->Print("v");
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
      (*TotalPDFPsiP)->plotOn(myFramePsiP, Components(*PsiP), LineStyle(7), LineColor(kBlue));
      if (justKeepPsi != true) (*TotalPDFPsiP)->plotOn(myFramePsiP, Components(*BkgmumuMassPsiP), LineStyle(4), LineColor(kRed));
      (*TotalPDFPsiP)->paramOn(myFramePsiP,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(RooArgSet(*nPsiPMass,*nBkgMassPsiP)));

      myString.clear(); myString.str("");
      myString << (*TotalPDFPsiP)->getPlotLabel() << "_paramBox";
      TPaveText* paveTextPsiP = (TPaveText*)myFramePsiP->findObject(myString.str().c_str());
      paveTextPsiP->AddText(Form("%s%.3f#pm%.3f","#mu#kern[-0.4]{#lower[-0.3]{_{1}}} = ",GetVar(*TotalPDFPsiP,"meanPsiP1")->getVal(),GetVar(*TotalPDFPsiP,"meanPsiP1")->getError()));
      paveTextPsiP->AddText(Form("%s%.3f#pm%.3f","#mu#kern[-0.4]{#lower[-0.3]{_{2}}} = ",GetVar(*TotalPDFPsiP,"meanPsiP2")->getVal(),GetVar(*TotalPDFPsiP,"meanPsiP2")->getError()));
      paveTextPsiP->AddText(Form("%s%.3f#pm%.3f","< #sigma > = ",sigmaPsiP,sigmaPsiPE));
      paveTextPsiP->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFramePsiP->chiSquare((*TotalPDFPsiP)->getPlotLabel(),dataSetPsiP->GetName())));
      CheckGoodFit(PsiPFitResult,paveTextPsiP);
      paveTextPsiP->SetBorderSize(0.0);
      paveTextPsiP->SetFillStyle(0);
      paveTextPsiP->SetTextSize(0.03);
      DrawString(LUMI,myFramePsiP);
      myFramePsiP->Draw();

      TLegend* legPsiP = new TLegend(0.75, 0.65, 0.97, 0.88, "");
      for (int i = 0; i < nElements; i++)
	{
	  TString objName = myFramePsiP->nameOf(i);
	  if (objName == "") continue;
	  TObject* obj = myFramePsiP->findObject(objName.Data());
	  legPsiP->AddEntry(obj,legNamesPsi[i],"lp");
	}
      legPsiP->SetFillStyle(0);
      legPsiP->SetFillColor(0);
      legPsiP->SetTextSize(0.04);
      legPsiP->SetBorderSize(0);
      legPsiP->Draw("same");


      // ##############
      // # Save plots #
      // ##############
      Canv->Modified();
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


      fileFitResults << "\nAmplitude of exclusion zone J/psi and psi(2S): +/- " << Utility->GetGenericParam("NSigmaPsi").c_str() << " Sigmas" << endl;
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
  B0MassArb          = new RooRealVar("B0MassArb","M(#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{+}}}}#kern[-0.3]{#pi}#kern[-0.3]{#lower[0.6]{^{#font[122]{\55}}}}#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#kern[-0.1]{#mu}#kern[-1.3]{#lower[0.6]{^{#font[122]{\55}}}})",Utility->B0Mass - atof(Utility->GetGenericParam("B0MassIntervalLeft").c_str()),Utility->B0Mass + atof(Utility->GetGenericParam("B0MassIntervalRight").c_str()),"GeV");
  mumuMass           = new RooRealVar("mumuMass","#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#kern[-0.1]{#mu}#kern[-1.3]{#lower[0.6]{^{#font[122]{\55}}}} inv. mass",0.0,6.0,"GeV");
  mumuMassE          = new RooRealVar("mumuMassE","#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#kern[-0.1]{#mu}#kern[-1.3]{#lower[0.6]{^{#font[122]{\55}}}} inv. mass error",0.0,0.5,"GeV");
  CosThetaKArb       = new RooRealVar("CosThetaKArb","cos(#theta#lower[-0.4]{_{#font[122]{K}}})",-1.0,1.0,"");
  CosThetaMuArb      = new RooRealVar("CosThetaMuArb","cos(#theta#lower[-0.4]{_{#font[12]{l}}})",-1.0,1.0,"");
  PhiKstMuMuPlaneArb = new RooRealVar("PhiKstMuMuPlaneArb","Angle (#mu#mu)--(#font[122]{K}#lower[0.4]{^{#font[122]{+}}}#pi#lower[0.4]{^{#font[122]{\55}}}) planes",-Utility->PI,Utility->PI,"rad");
  truthMatchSignal   = new RooRealVar("truthMatchSignal","Truth matching",0.0,1.0,"bool");
  rightFlavorTag     = new RooRealVar("rightFlavorTag","Right flavor tag",0.0,1.0,"bool");


  if (!(((FitType >= 21) && (FitType <= 26)) || ((FitType >= 81) && (FitType <= 86)) || (FitType == 96)))
    {
      RooArgSet Vars;
      Vars.add(*B0MassArb);
      Vars.add(*mumuMass);
      Vars.add(*mumuMassE);
      Vars.add(*CosThetaKArb);
      Vars.add(*CosThetaMuArb);
      Vars.add(*PhiKstMuMuPlaneArb);
      Vars.add(*truthMatchSignal);
      Vars.add(*rightFlavorTag);

      SingleCandNTuple           = new RooDataSet("SingleCandNTuple"          ,"SingleCandNTuple"          ,Vars);
      SingleCandNTuple_JPsi      = new RooDataSet("SingleCandNTuple_JPsi"     ,"SingleCandNTuple_JPsi"     ,Vars);
      SingleCandNTuple_PsiP      = new RooDataSet("SingleCandNTuple_PsiP"     ,"SingleCandNTuple_PsiP"     ,Vars);
      SingleCandNTuple_RejectPsi = new RooDataSet("SingleCandNTuple_RejectPsi","SingleCandNTuple_RejectPsi",Vars);
      SingleCandNTuple_KeepPsi   = new RooDataSet("SingleCandNTuple_KeepPsi"  ,"SingleCandNTuple_KeepPsi"  ,Vars);


      // #############################
      // # Load values from the tree #
      // #############################
      NTuple->ClearNTuple();
      NTuple->SetBranchAddresses(theTree);
      int nEntries = theTree->GetEntries();
      cout << "[ExtractYield::MakeDataSets]\t@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;
      for (int entry = 0; entry < nEntries; entry++)
	{
	  theTree->GetEntry(entry);
	  
	  if ((!(FitType == 36) && !(FitType == 56) && !(FitType == 76) &&
	       (NTuple->B0MassArb > Utility->B0Mass - atof(Utility->GetGenericParam("B0MassIntervalLeft").c_str())) &&
	       (NTuple->B0MassArb < Utility->B0Mass + atof(Utility->GetGenericParam("B0MassIntervalRight").c_str()))) ||
	      
	      ((FitType == 36) || (FitType == 56) || (FitType == 76)))
	    {
	      Vars.setRealValue("B0MassArb",         NTuple->B0MassArb);
	      Vars.setRealValue("mumuMass",          NTuple->mumuMass->at(0));
	      Vars.setRealValue("mumuMassE",         NTuple->mumuMassE->at(0));
	      Vars.setRealValue("CosThetaKArb",      NTuple->CosThetaKArb);
	      Vars.setRealValue("CosThetaMuArb",     NTuple->CosThetaMuArb);
	      Vars.setRealValue("PhiKstMuMuPlaneArb",NTuple->PhiKstMuMuPlaneArb);
	      Vars.setRealValue("truthMatchSignal",  NTuple->truthMatchSignal->at(0));
	      Vars.setRealValue("rightFlavorTag",    NTuple->rightFlavorTag);


	      // ########################
	      // # NTuple with all data #
	      // ########################
	      SingleCandNTuple->add(Vars);


	      // ###########################################################################
	      // # J/psi and psi(2S) keeping based on the event-by-event dimuon mass error #
	      // ###########################################################################
              if (((FitType == 36) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == true) && (Utility->PsiRejection(NTuple->B0MassArb,NTuple->mumuMass->at(0),NTuple->mumuMassE->at(0),"keepJpsi",false) == true)) ||

                  (((FitType == 56) || (FitType == 76)) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == true)) ||

                  ((((strcmp(CTRLfitWRKflow.c_str(),"trueGoodtag")      == 0) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == true)) ||
                    ((strcmp(CTRLfitWRKflow.c_str(),"trueMistag")       == 0) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == false)) ||
                    (((strcmp(CTRLfitWRKflow.c_str(),"trueAll&FitFrac") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"trueAll&NoFFrac") == 0)) && (NTuple->truthMatchSignal->at(0) == true)) ||
                    (strcmp(CTRLfitWRKflow.c_str(),"allEvts")           == 0)) &&
                   (Utility->PsiRejection(NTuple->B0MassArb,NTuple->mumuMass->at(0),NTuple->mumuMassE->at(0),"keepJpsi",false) == true)))

                SingleCandNTuple_JPsi->add(Vars);

              if (((FitType == 36) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == true) && (Utility->PsiRejection(NTuple->B0MassArb,NTuple->mumuMass->at(0),NTuple->mumuMassE->at(0),"keepPsiP",false) == true)) ||

                  (((FitType == 56) || (FitType == 76)) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == true)) ||

                  ((((strcmp(CTRLfitWRKflow.c_str(),"trueGoodtag")      == 0) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == true)) ||
                    ((strcmp(CTRLfitWRKflow.c_str(),"trueMistag")       == 0) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == false)) ||
                    (((strcmp(CTRLfitWRKflow.c_str(),"trueAll&FitFrac") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"trueAll&NoFFrac") == 0)) && (NTuple->truthMatchSignal->at(0) == true)) ||
                    (strcmp(CTRLfitWRKflow.c_str(),"allEvts")           == 0)) &&
                   (Utility->PsiRejection(NTuple->B0MassArb,NTuple->mumuMass->at(0),NTuple->mumuMassE->at(0),"keepPsiP",false) == true)))
		
		SingleCandNTuple_PsiP->add(Vars);
	      
	      
	      // #############################################################################
	      // # J/psi and psi(2S) rejection based on the event-by-event dimuon mass error #
	      // #############################################################################
	      if (((!(FitType == 36) && !(FitType == 56) && !(FitType == 76)) &&
		   
		   (((strcmp(CTRLfitWRKflow.c_str(),"trueGoodtag")      == 0) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == true)) ||
		    ((strcmp(CTRLfitWRKflow.c_str(),"trueMistag")       == 0) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == false)) ||
		    (((strcmp(CTRLfitWRKflow.c_str(),"trueAll&FitFrac") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"trueAll&NoFFrac") == 0)) && (NTuple->truthMatchSignal->at(0) == true)) ||
		    (strcmp(CTRLfitWRKflow.c_str(),"allEvts")           == 0)) &&
		   (Utility->PsiRejection(NTuple->B0MassArb,NTuple->mumuMass->at(0),NTuple->mumuMassE->at(0),"rejectPsi",true) == true)) ||
		  
		  (((FitType == 36) || (FitType == 56) || (FitType == 76)) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == true)))
		
		SingleCandNTuple_RejectPsi->add(Vars);


	      // ###########################################################################
	      // # J/psi and psi(2S) keeping based on the event-by-event dimuon mass error #
	      // ###########################################################################
	      if (((!(FitType == 36) && !(FitType == 56) && !(FitType == 76)) &&
		   
		   (((strcmp(CTRLfitWRKflow.c_str(),"trueGoodtag")      == 0) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == true)) ||
		    ((strcmp(CTRLfitWRKflow.c_str(),"trueMistag")       == 0) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == false)) ||
		    (((strcmp(CTRLfitWRKflow.c_str(),"trueAll&FitFrac") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"trueAll&NoFFrac") == 0)) && (NTuple->truthMatchSignal->at(0) == true)) ||
		    (strcmp(CTRLfitWRKflow.c_str(),"allEvts")           == 0)) &&
		   (Utility->PsiRejection(NTuple->B0MassArb,NTuple->mumuMass->at(0),NTuple->mumuMassE->at(0),"keepPsi",false) == true)) ||
		  
		  (((FitType == 36) || (FitType == 56) || (FitType == 76)) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == true)))

		SingleCandNTuple_KeepPsi->add(Vars);
	    }
	}


      cout << "\n[ExtractYield::MakeDataSets]\t@@@ NTuple with all data @@@" << endl;
      SingleCandNTuple->Print("v");

      cout << "\n[ExtractYield::MakeDataSets]\t@@@ NTuple with J/psi region @@@" << endl;
      SingleCandNTuple_JPsi->Print("v");
      
      cout << "\n[ExtractYield::MakeDataSets]\t@@@ NTuple with psi(2S) region @@@" << endl;
      SingleCandNTuple_PsiP->Print("v");
      
      cout << "\n[ExtractYield::MakeDataSets]\t@@@ NTuple without J/psi and psi(2S) regions @@@" << endl;
      SingleCandNTuple_RejectPsi->Print("v");

      cout << "\n[ExtractYield::MakeDataSets]\t@@@ NTuple with just J/psi and psi(2S) regions @@@" << endl;
      SingleCandNTuple_KeepPsi->Print("v");


      // ########################################
      // # Measure the J/psi and psi(2S) widths #
      // ########################################
      if (MakeMuMuPlots == true)
	{
	  TCanvas* cmumuMass_BeforeRej = new TCanvas("cmumuMass_BeforeRej","cmumuMass_BeforeRej",10, 10, 700, 900);
	  FitDimuonInvMass(SingleCandNTuple,&TotalPDFJPsi,&TotalPDFPsiP,mumuMass,cmumuMass_BeforeRej,false,false,"TotalPDFPsi_BeforeRej");

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
  

  // ####################################################
  // # Setting initial values for independent variables #
  // ####################################################
  B0MassArb->setVal(Utility->B0Mass);
  CosThetaKArb->setVal(0.0);
  CosThetaMuArb->setVal(0.0);
  PhiKstMuMuPlaneArb->setVal(0.0);
}


unsigned int GetSignalType (unsigned int FitType, vector<double>* q2Bins, int q2BinIndx)
{
  unsigned int SignalType;

  if (q2BinIndx < 0)
    {
      cout << "[ExtractYield::GetSignalType]\tIncorrect bin index for signal type determination : " << q2BinIndx << endl;
      exit (EXIT_FAILURE);
    }

  if      (((FitType >= 01)    && (FitType < 40))    || (((FitType == 96)    || ((FitType >= 81)    && (FitType <= 86)))    && (q2BinIndx != Utility->GetJPsiBin(q2Bins)) && (q2BinIndx != Utility->GetPsiPBin(q2Bins)))) SignalType = 1;
  else if (((FitType >= 41)    && (FitType < 60))    || (((FitType == 96)    || ((FitType >= 81)    && (FitType <= 86)))    && (q2BinIndx == Utility->GetJPsiBin(q2Bins))))                                               SignalType = 3;
  else if (((FitType >= 61)    && (FitType < 80))    || (((FitType == 96)    || ((FitType >= 81)    && (FitType <= 86)))    && (q2BinIndx == Utility->GetPsiPBin(q2Bins))))                                               SignalType = 5;
  else if (((FitType >= 01*10) && (FitType < 40*10)) || (((FitType == 96*10) || ((FitType >= 81*10) && (FitType <= 86*10))) && (q2BinIndx != Utility->GetJPsiBin(q2Bins)) && (q2BinIndx != Utility->GetPsiPBin(q2Bins)))) SignalType = 1;
  else if (((FitType >= 41*10) && (FitType < 60*10)) || (((FitType == 96*10) || ((FitType >= 81*10) && (FitType <= 86*10))) && (q2BinIndx == Utility->GetJPsiBin(q2Bins))))                                               SignalType = 3;
  else if (((FitType >= 61*10) && (FitType < 80*10)) || (((FitType == 96*10) || ((FitType >= 81*10) && (FitType <= 86*10))) && (q2BinIndx == Utility->GetPsiPBin(q2Bins))))                                               SignalType = 5;
  else
    {
      cout << "[ExtractYield::GetSignalType]\tFit type not valid : " << FitType << endl;
      exit (EXIT_FAILURE);
    }

  return SignalType;
}




// ==================
// ===> 1D MODEL <===
// ==================

void InstantiateMassFit (RooAbsPdf** TotalPDF, RooRealVar* x, string fitName, vector<vector<unsigned int>*>* configParam, unsigned int q2BinIndx)
{
  // ################################
  // # Read configuration variables #
  // ################################
  unsigned int useSignal = configParam->operator[](Utility->GetConfigParamIndx("SigType"))->operator[](q2BinIndx);
  unsigned int usePeakB  = configParam->operator[](Utility->GetConfigParamIndx("PeakBkgType"))->operator[](q2BinIndx);
  unsigned int useCombB  = configParam->operator[](Utility->GetConfigParamIndx("CombBkgType"))->operator[](q2BinIndx);
  unsigned int useMisTag = configParam->operator[](Utility->GetConfigParamIndx("MistagType"))->operator[](q2BinIndx);


  // ###################
  // # Local variables #
  // ###################
  stringstream myString;


  // ################################################
  // # Define mass fit variables and pdf for signal #
  // ################################################
  meanS   = new RooRealVar("meanS","Signal mean",Utility->B0Mass,"GeV");

  sigmaS1 = new RooRealVar("sigmaS1","Signal sigma-1",0.0,"GeV");
  MassS1  = new RooGaussian("MassS1","Signal Gaussian-1",*x,*meanS,*sigmaS1);

  sigmaS2 = new RooRealVar("sigmaS2","Signal sigma-2",0.0,"GeV");
  MassS2  = new RooGaussian("MassS2","Signal Gaussian-2",*x,*meanS,*sigmaS2);

  fracMassS = new RooRealVar("fracMassS","Fraction of signal Gaussian",0.0,0.0,1.0);

  if      (useSignal == 1) MassSignal = new RooGaussian(*((RooGaussian*)MassS1),"MassSignal");
  else if (useSignal == 2) MassSignal = new RooAddPdf("MassSignal","Signal mass pdf",RooArgSet(*MassS1,*MassS2),RooArgSet(*fracMassS));
  else                     MassSignal = NULL;


  // ##################################################################
  // # Define mass fit variables and pdf for combinatorial background #
  // ##################################################################
  var1 = new RooRealVar("var1","First background variables",0.0,"GeV");
  myString.clear(); myString.str("");
  myString << "exp(-(" << x->getPlotLabel() << " - meanS) / var1)";
  BkgMassExp1 = new RooGenericPdf("BkgMassExp1",myString.str().c_str(),RooArgSet(*x,*meanS,*var1));

  var2 = new RooRealVar("var2","Second background variable",0.0,"GeV");
  myString.clear(); myString.str("");
  myString << "exp(-(" << x->getPlotLabel() << " - meanS) / var2)";
  BkgMassExp2 = new RooGenericPdf("BkgMassExp2",myString.str().c_str(),RooArgSet(*x,*meanS,*var2));

  fracMassBExp = new RooRealVar("fracMassBExp","Fraction of background Exponential",0.0,0.0,1.0);

  if      (useCombB == 1) BkgMassComb = new RooGenericPdf(*((RooGenericPdf*)BkgMassExp1),"BkgMassComb");
  else if (useCombB == 2) BkgMassComb = new RooAddPdf("BkgMassComb","Background mass comb. bkg pdf",RooArgSet(*BkgMassExp1,*BkgMassExp2),RooArgSet(*fracMassBExp));
  else if (useCombB == 3)
    {
      myString.clear(); myString.str("");
      myString << "TMath::Erfc((" << x->getPlotLabel() << " - var1) / sqrt(2*var2*var2))";
      BkgMassComb = new RooGenericPdf("BkgMass",myString.str().c_str(),RooArgSet(*x,*var1,*var2));
    }
  else BkgMassComb = NULL;


  // #######################################################
  // # Define mass fit variables and pdf for mistag signal #
  // #######################################################
  sigmaMisTag1 = new RooRealVar("sigmaMisTag1","Mistag sigma-1",0.0,"GeV");
  MassMisTag1  = new RooGaussian("MassMisTag1","Mistag-1",*x,*meanS,*sigmaMisTag1);

  sigmaMisTag2 = new RooRealVar("sigmaMisTag2","Mistag sigma-2",0.0,"GeV");
  MassMisTag2  = new RooGaussian("MassMisTag2","Mistag-2",*x,*meanS,*sigmaMisTag2);

  fracMisTag = new RooRealVar("fracMisTag","Fraction mistag Gaussian",0.0,0.0,1.0);

  if      (useMisTag == 1) MassMisTag = new RooGaussian(*((RooGaussian*)MassMisTag1),"MassMisTag");
  else if (useMisTag == 2) MassMisTag = new RooAddPdf("MassMisTag","Mistag mass pdf",RooArgSet(*MassMisTag1,*MassMisTag2),RooArgSet(*fracMisTag));
  else                     MassMisTag = NULL;


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
  BkgMassRPeak   = new RooAddPdf("BkgMassRPeak","Right peaking bkg mass pdf",RooArgSet(*BkgMassRPeak1,*BkgMassRPeak2),RooArgSet(*fracMassBRPeak));

  meanL1         = new RooRealVar("meanL1","Bkg left peak mean-1",0.0,"GeV");
  sigmaL1        = new RooRealVar("sigmaL1","Bkg left peak sigma-1",0.0,"GeV");
  BkgMassLPeak1  = new RooGaussian("BkgMassLPeak1","Bkg left peak-1",*x,*meanL1,*sigmaL1);

  meanL2         = new RooRealVar("meanL2","Bkg left peak mean-2",0.0,"GeV");
  sigmaL2        = new RooRealVar("sigmaL2","Bkg left peak sigma-2",0.0,"GeV");
  BkgMassLPeak2  = new RooGaussian("BkgMassLPeak2","Bkg left peak-2",*x,*meanL2,*sigmaL2);

  fracMassBLPeak = new RooRealVar("fracMassBLPeak","Fraction of background left Peak",0.0,0.0,1.0);
  BkgMassLPeak   = new RooAddPdf("BkgMassLPeak","Left peaking bkg mass pdf",RooArgSet(*BkgMassLPeak1,*BkgMassLPeak2),RooArgSet(*fracMassBLPeak));

  fracMassBPeak  = new RooRealVar("fracMassBPeak","Fraction of background right-left Peak",0.0,0.0,1.0);

  if      (usePeakB == 1)  BkgMassPeak = new RooGaussian(*((RooGaussian*)BkgMassRPeak1),"BkgMassPeak");
  else if (usePeakB == 2)  BkgMassPeak = new RooAddPdf(*((RooAddPdf*)BkgMassRPeak),"BkgMassPeak");
  else if (usePeakB == 11) BkgMassPeak = new RooAddPdf("BkgMassPeak","Peaking bkg mass pdf",RooArgSet(*BkgMassRPeak1,*BkgMassLPeak1),RooArgSet(*fracMassBPeak));
  else if (usePeakB == 12) BkgMassPeak = new RooAddPdf("BkgMassPeak","Peaking bkg mass pdf",RooArgSet(*BkgMassRPeak,*BkgMassLPeak),RooArgSet(*fracMassBPeak));
  else                     BkgMassPeak = NULL;


  // ###########################
  // # Define pdf coefficients #
  // ###########################
  nSig     = new RooRealVar("nSig","Number of signal events",1.0);
  nBkgComb = new RooRealVar("nBkgComb","Number of comb. background events",1.0);
  RooFormulaVar* nMisTag;
  if ((useSignal == 0) && (usePeakB == 0) && (useCombB == 0) && (useMisTag != 0))
    {
      nMisTagFrac = new RooRealVar("nMisTagFrac","Number of mistag events",1.0);
      nMisTag     = new RooFormulaVar("nMisTag","nMisTagFrac", RooArgSet(*nMisTagFrac));
    }
  else
    {
      nMisTagFrac = new RooRealVar("nMisTagFrac","Fraction of mistag",0.0,0.0,1.0);
      nMisTag     = new RooFormulaVar("nMisTag","nSig * nMisTagFrac / (1 - nMisTagFrac)", RooArgSet(*nSig,*nMisTagFrac));
    }
  nBkgPeak = new RooRealVar("nBkgPeak","Number of peaking background events",1.0);


  if      ((useSignal != 0) && (usePeakB == 0) && (useCombB == 0) && (useMisTag == 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*MassSignal),RooArgSet(*nSig));

  else if ((useSignal == 0) && (usePeakB != 0) && (useCombB == 0) && (useMisTag == 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*BkgMassPeak),RooArgSet(*nBkgPeak));

  else if ((useSignal == 0) && (usePeakB == 0) && (useCombB != 0) && (useMisTag == 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*BkgMassComb),RooArgSet(*nBkgComb));

  else if ((useSignal == 0) && (usePeakB == 0) && (useCombB == 0) && (useMisTag != 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*MassMisTag),RooArgSet(*nMisTag));

  else if ((useSignal != 0) && (usePeakB == 0) && (useCombB != 0) && (useMisTag == 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*MassSignal,*BkgMassComb),RooArgSet(*nSig,*nBkgComb));

  else if ((useSignal != 0) && (usePeakB == 0) && (useCombB == 0) && (useMisTag != 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*MassSignal,*MassMisTag),RooArgSet(*nSig,*nMisTag));

  else if ((useSignal != 0) && (usePeakB != 0) && (useCombB != 0) && (useMisTag == 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*MassSignal,*BkgMassComb,*BkgMassPeak),RooArgSet(*nSig,*nBkgComb,*nBkgPeak));

  else if ((useSignal != 0) && (usePeakB == 0) && (useCombB != 0) && (useMisTag != 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*MassSignal,*BkgMassComb,*MassMisTag),RooArgSet(*nSig,*nBkgComb,*nMisTag));

  else if ((useSignal != 0) && (usePeakB != 0) && (useCombB != 0) && (useMisTag != 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*MassSignal,*BkgMassComb,*MassMisTag,*BkgMassPeak),RooArgSet(*nSig,*nBkgComb,*nMisTag,*nBkgPeak));

  else
    {
      cout << "[ExtractYield::InstantiateMassFit]\tIncorrect configuration sequence : useSignal = " << useSignal << "\tusePeakB = " << usePeakB << "\tuseCombB = " << useCombB << "\tuseMisTag = " << useMisTag << endl;
      exit (EXIT_FAILURE);
    }
}


RooFitResult* MakeMassFit (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* x, TCanvas* Canv, RooArgSet* vecConstr, double* NLLvalue, TPaveText* extText, unsigned int ID)
{
  // ###################
  // # Local variables #
  // ###################
  RooFitResult* fitResult = NULL;
  RooArgSet VarsYield;
  stringstream myString;
  double signalSigmaGoodT  = 0.0;
  double signalSigmaGoodTE = 0.0;
  double signalSigmaMisT   = 0.0;
  double signalSigmaMisTE  = 0.0;
  double signalSigmaT      = 0.0;
  double signalSigmaTE     = 0.0;
  double totalYield        = 0.0;
  double totalYieldE       = 0.0;
  unsigned int nElements   = 0;
  TString legNames[6];


  // ###################
  // # Make actual fit #
  // ###################
  if (atoi(Utility->GetGenericParam("ApplyConstr").c_str()) == true) fitResult = (*TotalPDF)->fitTo(*dataSet,Extended(true),ExternalConstraints(*vecConstr),Save(true),Minos(USEMINOS));
  else                                                               fitResult = (*TotalPDF)->fitTo(*dataSet,Extended(true),Save(true),Minos(USEMINOS));


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


  // ################
  // # Plot results #
  // ################
  Canv->cd();
  RooPlot* myFrameX = x->frame(NBINS);

  dataSet->plotOn(myFrameX,Name(MakeName(dataSet,ID).c_str()));
  legNames[nElements++] = "Data";

  if (FUNCERRBAND == true) (*TotalPDF)->plotOn(myFrameX,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),VisualizeError(*fitResult,1,true),VLines(),FillColor(kGreen-7));
  else                     (*TotalPDF)->plotOn(myFrameX,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack));
  legNames[nElements++] = "Total p.d.f.";

  if (GetVar(*TotalPDF,"nSig") != NULL)
    {
      (*TotalPDF)->plotOn(myFrameX, Components(*MassSignal), LineStyle(7), LineColor(kBlue));
      legNames[nElements++] = "Right-tag sig";
      VarsYield.add(*nSig);
    }

  if (GetVar(*TotalPDF,"nMisTagFrac") != NULL)
    {
      (*TotalPDF)->plotOn(myFrameX, Components(*MassMisTag), LineStyle(8), LineColor(kAzure+6));
      legNames[nElements++] = "Mis-tag sig";
      VarsYield.add(*nMisTagFrac);
    }

  if (GetVar(*TotalPDF,"nBkgComb") != NULL)
    {
      (*TotalPDF)->plotOn(myFrameX, Components(*BkgMassComb), LineStyle(4), LineColor(kRed));
      legNames[nElements++] = "Comb. bkg";
      VarsYield.add(*nBkgComb);
    }

  if (GetVar(*TotalPDF,"nBkgPeak") != NULL)
    {
      (*TotalPDF)->plotOn(myFrameX, Components(*BkgMassPeak), LineStyle(3), LineColor(kViolet));
      legNames[nElements++] = "Peak. bkg";
      VarsYield.add(*nBkgPeak);
    }

  (*TotalPDF)->paramOn(myFrameX,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(VarsYield),ShowConstants(true));
  // Format options:
  // - "N" add name
  // - "E" add error
  // - "A" show asymmetric error
  // - "U" show unit
  // - "H" hide the value

  
  myString.clear(); myString.str("");
  myString << (*TotalPDF)->getPlotLabel() << "_paramBox";
  TPaveText* paveTextX = (TPaveText*)myFrameX->findObject(myString.str().c_str());

  if (GetVar(*TotalPDF,"nSig") != NULL)
    {
      if (GetVar(*TotalPDF,"fracMassS") != NULL)
	{
	  signalSigmaGoodT  = sqrt( GetVar(*TotalPDF,"fracMassS")->getVal() * pow(GetVar(*TotalPDF,"sigmaS1")->getVal(),2.) + (1. - GetVar(*TotalPDF,"fracMassS")->getVal()) * pow(GetVar(*TotalPDF,"sigmaS2")->getVal(),2.) );
	  signalSigmaGoodTE = 1./(2.*signalSigmaGoodT) * sqrt( pow((pow(GetVar(*TotalPDF,"sigmaS1")->getVal(),2.) - pow(GetVar(*TotalPDF,"sigmaS2")->getVal(),2.)) * GetVar(*TotalPDF,"fracMassS")->getError(),2.) +
							       pow(2. * GetVar(*TotalPDF,"fracMassS")->getVal()        * GetVar(*TotalPDF,"sigmaS1")->getVal() * GetVar(*TotalPDF,"sigmaS1")->getError(),2.) +
							       pow(2. * (1. - GetVar(*TotalPDF,"fracMassS")->getVal()) * GetVar(*TotalPDF,"sigmaS2")->getVal() * GetVar(*TotalPDF,"sigmaS2")->getError(),2.) );
	}
      else if (GetVar(*TotalPDF,"sigmaS1") != NULL)
	{
	  signalSigmaGoodT  = GetVar(*TotalPDF,"sigmaS1")->getVal();
	  signalSigmaGoodTE = GetVar(*TotalPDF,"sigmaS1")->getError();
	}

      signalSigmaT  = signalSigmaGoodT;
      signalSigmaTE = signalSigmaGoodTE;
    }

  if (GetVar(*TotalPDF,"nMisTagFrac") != NULL)
    {
      if (GetVar(*TotalPDF,"fracMisTag") != NULL)
	{
	  signalSigmaMisT  = sqrt( GetVar(*TotalPDF,"fracMisTag")->getVal() * pow(GetVar(*TotalPDF,"sigmaMisTag1")->getVal(),2.) + (1. - GetVar(*TotalPDF,"fracMisTag")->getVal()) * pow(GetVar(*TotalPDF,"sigmaMisTag2")->getVal(),2.) );
	  signalSigmaMisTE = 1./(2.*signalSigmaMisT) * sqrt( pow((pow(GetVar(*TotalPDF,"sigmaMisTag1")->getVal(),2.) - pow(GetVar(*TotalPDF,"sigmaMisTag2")->getVal(),2.)) * GetVar(*TotalPDF,"fracMisTag")->getError(),2.) +
							     pow(2. * GetVar(*TotalPDF,"fracMisTag")->getVal()        * GetVar(*TotalPDF,"sigmaMisTag1")->getVal() * GetVar(*TotalPDF,"sigmaMisTag1")->getError(),2.) +
							     pow(2. * (1. - GetVar(*TotalPDF,"fracMisTag")->getVal()) * GetVar(*TotalPDF,"sigmaMisTag2")->getVal() * GetVar(*TotalPDF,"sigmaMisTag2")->getError(),2.) );
	}
      else if (GetVar(*TotalPDF,"sigmaMisTag1") != NULL)
	{
	  signalSigmaMisT  = GetVar(*TotalPDF,"sigmaMisTag1")->getVal();
	  signalSigmaMisTE = GetVar(*TotalPDF,"sigmaMisTag1")->getError();
	}

      signalSigmaT  = signalSigmaMisT;
      signalSigmaTE = signalSigmaMisTE;
    }

  if ((GetVar(*TotalPDF,"nSig") != NULL) && (GetVar(*TotalPDF,"nMisTagFrac") != NULL))
    {
      totalYield  = GetVar(*TotalPDF,"nSig")->getVal() / (1. - GetVar(*TotalPDF,"nMisTagFrac")->getVal());
      totalYieldE = totalYield * sqrt( pow(GetVar(*TotalPDF,"nSig")->getError()/GetVar(*TotalPDF,"nSig")->getVal(),2.) + pow(GetVar(*TotalPDF,"nMisTagFrac")->getError()/(1. - GetVar(*TotalPDF,"nMisTagFrac")->getVal()),2.) );
      
      double fraction  = (1. - 2.*GetVar(*TotalPDF,"nMisTagFrac")->getVal()) / (1. - GetVar(*TotalPDF,"nMisTagFrac")->getVal());
      double fractionE = GetVar(*TotalPDF,"nMisTagFrac")->getError() / pow((1. - GetVar(*TotalPDF,"nMisTagFrac")->getVal()),2.);

      signalSigmaT  = sqrt( fraction * pow(signalSigmaGoodT,2.) + (1. - fraction) * pow(signalSigmaMisT,2.) );
      signalSigmaTE = 1./(2.*signalSigmaT) * sqrt( pow((pow(signalSigmaGoodT,2.) - pow(signalSigmaMisT,2.)) * fractionE,2.) +
						   pow(2. * fraction        * signalSigmaGoodT * signalSigmaGoodTE,2.) +
						   pow(2. * (1. - fraction) * signalSigmaMisT  * signalSigmaMisTE,2.) );
      
      paveTextX->AddText(Form("%s%1.f#pm%1.f","Total sig yield = ",totalYield,totalYieldE));
    }
  else if (GetVar(*TotalPDF,"nSig") != NULL)
    {
      totalYield  = GetVar(*TotalPDF,"nSig")->getVal();
      totalYieldE = GetVar(*TotalPDF,"nSig")->getError();
    }

  if ((GetVar(*TotalPDF,"nSig") != NULL) || (GetVar(*TotalPDF,"nMisTagFrac") != NULL))
    {
      paveTextX->AddText(Form("%s%.3f#pm%.3f","#mu = ",GetVar(*TotalPDF,"meanS")->getVal(),GetVar(*TotalPDF,"meanS")->getError()));
      paveTextX->AddText(Form("%s%.3f#pm%.3f","Total < #sigma > = ",signalSigmaT,signalSigmaTE));
    }
  paveTextX->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameX->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
  
  TLegend* legX = new TLegend(0.75, 0.65, 0.97, 0.88, "");
  for (unsigned int i = 0; i < nElements; i++)
    {
      TString objName = myFrameX->nameOf(i);
      if (objName == "") continue;
      TObject* obj = myFrameX->findObject(objName.Data());
      legX->AddEntry(obj,legNames[i],"lp");
    }
  legX->SetFillStyle(0);
  legX->SetFillColor(0);
  legX->SetTextSize(0.04);
  legX->SetBorderSize(0);


  // ####################
  // # Save fit results #
  // ####################
  fileFitResults << "====================================================================" << endl;
  fileFitResults << "@@@@@@ Make mass fit @@@@@@" << endl;

  fileFitResults << "Chi2/DoF fit: " << myFrameX->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str());
  fileFitResults << "; p-value: " << TMath::Prob(myFrameX->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())*NBINS,NBINS) << endl;

  *NLLvalue = StoreFitResultsInFile(TotalPDF,fitResult,dataSet,vecConstr);
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
      if (GetVar(*TotalPDF,"nSig") != NULL) extText->AddText(Form("%s%1.f#pm%1.f","Signal yield: ",totalYield,totalYieldE));
      extText->Paint();
      myFrameX->addObject(extText);
    }
  myFrameX->Draw();
  legX->Draw("same");


  // ##############
  // # Save plots #
  // ##############
  Canv->Modified();
  Canv->Update();
  if (SAVEPLOT == true)
    {
      myString.clear(); myString.str("");
      myString << (*TotalPDF)->getPlotLabel() << "_" << ID << ".pdf";
      Canv->Print(myString.str().c_str());
    }


  VarsYield.Clear();
  return fitResult;
}


void IterativeMassFitq2Bins (RooDataSet* dataSet,
			     bool useEffPDF,
			     double PsiYieldGoodTag, double PsiYieldGoodTagErr,
			     double PsiYieldMisTag,
			     RooRealVar* x,
			     int specBin,
			     unsigned int FitType,
			     vector<TH1D*>* VecHistoMeas,
			     vector<double>* q2Bins,
			     vector<vector<unsigned int>*>* configParam, vector<vector<string>*>* fitParam,
			     RooArgSet* vecConstr,
			     unsigned int ID)
{
  // ###################
  // # Local variables #
  // ###################
  vector<string>* vecParStr = NULL;
  stringstream myString;

  double effPsiGoodTag     = 1.0;
  double effPsiGoodTagErr  = 0.0;

  double effPsiMisTag      = 1.0;
  double effPsiMisTagErr   = 0.0;

  double effMuMuGoodTag    = 1.0;
  double effMuMuGoodTagErr = 0.0;

  double effMuMuMisTag     = 1.0;
  double effMuMuMisTagErr  = 0.0;

  TCanvas*    cq2Bins[q2Bins->size()-1];
  RooDataSet* dataSet_q2Bins[q2Bins->size()-1];
  RooAbsPdf*  TotalPDFq2Bins[q2Bins->size()-1];
  TPaveText*  extText[q2Bins->size()-1];

  RooFitResult* fitResult;

  double NLLvalue;


  // ####################################
  // # Read I[S*E] for resonant channel #
  // ####################################
  if (useEffPDF == true)
    {
      // ######################
      // # Good-tagged events #
      // ######################
      Utility->ReadParVsq2Bins(PARAMETERFILEIN,"I[S*E]okTag",&vecParStr);
      myString.clear(); myString.str("");
      myString << vecParStr->operator[](atoi(Utility->GetGenericParam("NormJPSInotPSIP").c_str()) == true ? Utility->GetJPsiBin(q2Bins) : Utility->GetPsiPBin(q2Bins)).c_str();
      SetValueAndErrors(NULL,"",1.0,&myString,&effPsiGoodTag,&effPsiGoodTagErr,&effPsiGoodTagErr);
      cout << "\n[ExtractYield::IterativeMassFitq2Bins]\t@@@ Integral of S*E over angular variables for normalization channel good-tagged events: " << effPsiGoodTag << " +/- " << effPsiGoodTagErr << " @@@" << endl;


      // #####################
      // # Mis-tagged events #
      // #####################
      Utility->ReadParVsq2Bins(PARAMETERFILEIN,"I[S*E]misTag",&vecParStr);
      myString.clear(); myString.str("");
      myString << vecParStr->operator[](atoi(Utility->GetGenericParam("NormJPSInotPSIP").c_str()) == true ? Utility->GetJPsiBin(q2Bins) : Utility->GetPsiPBin(q2Bins)).c_str();
      SetValueAndErrors(NULL,"",1.0,&myString,&effPsiMisTag,&effPsiMisTagErr,&effPsiMisTagErr);
      cout << "\n[ExtractYield::IterativeMassFitq2Bins]\t@@@ Integral of S*E over angular variables for normalization channel mis-tagged events: " << effPsiMisTag << " +/- " << effPsiMisTagErr << " @@@" << endl;
    }
 
 
  fileFitResults << "\n@@@@@@@@@@@@@@@@@ Differential branching fraction @@@@@@@@@@@@@@@@@@" << endl;
  for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
    {
      if (!(FitType == 41) && !(FitType == 61) && (Utility->ValIsInPsi(q2Bins,(q2Bins->operator[](i+1)+q2Bins->operator[](i))/2.) == true))
      	{
      	  // ###################################
      	  // # Save zeros into prarameter file #
      	  // ###################################
      	  vecParStr = SaveFitResults(NULL,i,fitParam,configParam,vecConstr);
      	  Utility->SaveFitValues(PARAMETERFILEOUT,vecParStr,i);
	  vecParStr->clear();
      	  delete vecParStr;

      	  continue;
      	}


      // ######################
      // # Make external text #
      // ######################
      extText[i] = new TPaveText(0.65,0.5,0.97,0.63,"NDC");
      extText[i]->AddText(Form("%s%.2f%s%.2f%s","q#lower[0.4]{^{2}}: ",q2Bins->operator[](i)," #font[122]{\55} ",q2Bins->operator[](i+1)," GeV#lower[0.4]{^{2}}"));
      extText[i]->SetTextAlign(11);
      extText[i]->SetBorderSize(0.0);
      extText[i]->SetFillStyle(0);
      extText[i]->SetTextSize(0.04);


      fileFitResults << "\nBin[" << i << "]: " << q2Bins->operator[](i) << " <= q^2 < " << q2Bins->operator[](i+1) << endl;

      myString.clear(); myString.str("");
      myString << "c_" << i;
      cq2Bins[i] = new TCanvas(myString.str().c_str(), myString.str().c_str(), 20, 20, 700, 500);

      myString.clear(); myString.str("");
      myString << "(mumuMass*mumuMass) > " << q2Bins->operator[](i) << " && (mumuMass*mumuMass) <= " << q2Bins->operator[](i+1);
      cout << "\n[ExtractYield::IterativeMassFitq2Bins]\tCut string: " << myString.str() << endl;
      dataSet_q2Bins[i] = (RooDataSet*)dataSet->reduce(myString.str().c_str());
      cout << "[ExtractYield::IterativeMassFitq2Bins]\tNumber of events: " << dataSet_q2Bins[i]->sumEntries() << endl;

      unsigned int countMisTag  = 0;
      unsigned int countGoodTag = 0;
      for (int j = 0; j < static_cast<int>(dataSet_q2Bins[i]->sumEntries()); j++)
	{
	  if (dataSet_q2Bins[i]->get(j)->getRealValue("truthMatchSignal") == true)
	    {
	      if (dataSet_q2Bins[i]->get(j)->getRealValue("rightFlavorTag") == 0.0) countMisTag++;
	      else                                                                  countGoodTag++;
	    }
	}
      cout << "[ExtractYield::IterativeMassFitq2Bins]\tDynamic mis-tag fraction: " << static_cast<double>(countMisTag) / static_cast<double>(countMisTag + countGoodTag) << " = (" << countMisTag << "/(" << countMisTag << "+" << countGoodTag << "))" << endl;

      myString.clear(); myString.str("");
      myString << "TotalPDFq2Bin_" << i;
      InstantiateMassFit(&TotalPDFq2Bins[i],x,myString.str(),configParam,i);


      // #####################
      // # Initialize p.d.f. #
      // #####################
      CopyFitResults(TotalPDFq2Bins[i],i,fitParam,countMisTag,countGoodTag);


      // #####################
      // # Apply constraints #
      // #####################
      ClearVars(vecConstr);
      BuildMassConstraints(vecConstr,TotalPDFq2Bins[i],"sign");
      BuildMassConstraints(vecConstr,TotalPDFq2Bins[i],"peak");
      if ((strcmp(CTRLfitWRKflow.c_str(),"trueAll&NoFFrac") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"trueAll&FitFrac") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"allEvts") == 0)) BuildMassConstraints(vecConstr,TotalPDFq2Bins[i],"mistag");


      // ###################################
      // # Print variables and constraints #
      // ###################################
      PrintVariables(TotalPDFq2Bins[i]->getVariables(),"vars");
      PrintVariables(vecConstr,"cons");


      // ###################
      // # Perform the fit #
      // ###################
      fitResult = MakeMassFit(dataSet_q2Bins[i],&TotalPDFq2Bins[i],x,cq2Bins[i],vecConstr,&NLLvalue,extText[i],ID);
      if (CheckGoodFit(fitResult) == true) cout << "\n[ExtractYield::IterativeMassFitq2Bins]\t@@@ Fit converged ! @@@" << endl;
      else                                 cout << "\n[ExtractYield::IterativeMassFitq2Bins]\t@@@ Fit didn't converge ! @@@" << endl;


      // ##############################################
      // # Save fit results back into prarameter file #
      // ##############################################
      vecParStr = SaveFitResults(TotalPDFq2Bins[i],i,fitParam,configParam,vecConstr);
      Utility->SaveFitValues(PARAMETERFILEOUT,vecParStr,i);
      vecParStr->clear();
      delete vecParStr;
      vecParStr = NULL;


      if (useEffPDF == true)
	{
	  // ######################
	  // # Good-tagged events #
	  // ######################
          Utility->ReadParVsq2Bins(PARAMETERFILEIN,"I[S*E]okTag",&vecParStr);
	  myString.clear(); myString.str("");
	  myString << vecParStr->operator[](i).c_str();
	  SetValueAndErrors(NULL,"",1.0,&myString,&effMuMuGoodTag,&effMuMuGoodTagErr,&effMuMuGoodTagErr);
	  cout << "\n[ExtractYield::IterativeMassFitq2Bins]\t===> Integral of S*E over angular variables for signal good-tagged events <===" << endl;
	  cout << effMuMuGoodTag << "   -" << effMuMuGoodTagErr << "   " << effMuMuGoodTagErr << endl;

	  
	  // #####################
	  // # Mis-tagged events #
	  // #####################
	  Utility->ReadParVsq2Bins(PARAMETERFILEIN,"I[S*E]misTag",&vecParStr);
	  myString.clear(); myString.str("");
	  myString << vecParStr->operator[](i).c_str();
	  SetValueAndErrors(NULL,"",1.0,&myString,&effMuMuMisTag,&effMuMuMisTagErr,&effMuMuMisTagErr);
	  cout << "\n[ExtractYield::IterativeMassFitq2Bins]\t===> Integral of S*E over angular variables for signal mis-tagged events <===" << endl;
	  cout << effMuMuMisTag << "   -" << effMuMuMisTagErr << "   " << effMuMuMisTagErr << endl;
	}


      // ########################################
      // # Save outcome of the fit in histogram #
      // ########################################
      if (GetVar(TotalPDFq2Bins[i],"nSig") != NULL)
	{
	  double nEvGoodTag    = GetVar(TotalPDFq2Bins[i],"nSig")->getVal();
	  double nEvGoodTagErr = GetVar(TotalPDFq2Bins[i],"nSig")->getError();

	  double nEvMisTag = 0.0;
	  if (GetVar(TotalPDFq2Bins[i],"nMisTagFrac") != NULL)
	    nEvMisTag = GetVar(TotalPDFq2Bins[i],"nSig")->getVal() / (1. - GetVar(TotalPDFq2Bins[i],"nMisTagFrac")->getVal()) * GetVar(TotalPDFq2Bins[i],"nMisTagFrac")->getVal();
	  
	  double num       = (nEvGoodTag / effMuMuGoodTag) + (nEvMisTag / effMuMuMisTag);
	  double den       = (PsiYieldGoodTag / effPsiGoodTag) + (PsiYieldMisTag / effPsiMisTag);
	  double dBFdq2    = num / den * (atoi(Utility->GetGenericParam("NormJPSInotPSIP").c_str()) == true ? Utility->JPsiBF : Utility->PsiPBF) / (q2Bins->operator[](i+1) - q2Bins->operator[](i)) / 1e-7;
	  double dBFdq2Err = dBFdq2 * sqrt( (pow(num * nEvGoodTagErr      / nEvGoodTag,2.)      + pow(nEvGoodTag         / (effMuMuGoodTag*effMuMuGoodTag) * effMuMuGoodTagErr,2.) + pow(nEvMisTag          / (effMuMuMisTag*effMuMuMisTag) * effMuMuMisTagErr,2.)) / (num*num) +
					    (pow(den * PsiYieldGoodTagErr / PsiYieldGoodTag,2.) + pow(PsiYieldGoodTagErr / (effPsiGoodTag*effPsiGoodTag)   * effPsiGoodTagErr,2.)  + pow(PsiYieldGoodTagErr / (effPsiMisTag*effPsiMisTag)   * effPsiMisTagErr,2.))  / (den*den) );

	  VecHistoMeas->operator[](0)->SetBinContent(i+1,dBFdq2);
	  VecHistoMeas->operator[](0)->SetBinError(i+1,dBFdq2Err);

	  fileFitResults << "@@@@@@ dBF/dq^2 @@@@@@" << endl;
	  fileFitResults << "dBF/dq^2: " << VecHistoMeas->operator[](0)->GetBinContent(i+1) << " -/+ " << VecHistoMeas->operator[](0)->GetBinError(i+1) << endl;

	  fileFitResults << "===> dBF/dq^2 to cut and paste in config. file <===" << endl;
	  fileFitResults << VecHistoMeas->operator[](0)->GetBinContent(i+1) << "   -" << dBFdq2Err << "   " << dBFdq2Err << endl;

	  fileFitResults << "====================================================================" << endl;
	}


      // #############################################
      // # Save observables in systematic error file #
      // #############################################
      myString.clear(); myString.str("");
      if (CheckGoodFit(fitResult) == true) myString << ID << "   " << -2.0 << "   " << -2.0 << "   " << VecHistoMeas->operator[](0)->GetBinContent(i+1) << "   " << -2.0 << "   " << -2.0 << "   " << NLLvalue;
      else                                 myString << ID << "   " << -2.0 << "   " << -2.0 << "   " << -2.0 << "   " << -2.0 << "   " << -2.0 << "   " << -2.0;
      fileFitSystematics << myString.str() << endl;
    }
}


void MakeMassToy (RooAbsPdf* TotalPDF, RooRealVar* x, TCanvas* Canv, unsigned int nToy, int specBin, vector<vector<string>*>* fitParam, RooArgSet* vecConstr, string fileName)
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
  BuildMassConstraints(vecConstr,TotalPDF,"sign");
  BuildMassConstraints(vecConstr,TotalPDF,"peak");
  if ((strcmp(CTRLfitWRKflow.c_str(),"trueAll&NoFFrac") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"trueAll&FitFrac") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"allEvts") == 0)) BuildMassConstraints(vecConstr,TotalPDF,"mistag");


  // ###################################
  // # Print variables and constraints #
  // ###################################
  PrintVariables(TotalPDF->getVariables(),"vars");
  PrintVariables(vecConstr,"cons");


  // #############################
  // # Toy-MC generation and fit #
  // #############################
  RooMCStudy* MyToy = new RooMCStudy(*TotalPDF,*x,Extended(true),FitOptions(Extended(true),ExternalConstraints(*vecConstr),Minos(USEMINOS))); // Possible options : "Binned()" = faster; Silence()
  MyToy->generateAndFit(nToy,nEntryToy,true);
  Canv->Divide(6,5);


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

  if ((GetVar(TotalPDF,"var1") != NULL) && (GetVar(TotalPDF,"var1")->getError() != 0.0) && (IsInConstraints(vecConstr,"var1") == false))
    {
      tmpVar = GetVar(TotalPDF,"var1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("var1"))->operator[](specBin).c_str()) - 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("var1"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("var1"))->operator[](specBin).c_str()) + 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("var1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"var2") != NULL) && (GetVar(TotalPDF,"var2")->getError() != 0.0) && (IsInConstraints(vecConstr,"var2") == false))
    {
      tmpVar = GetVar(TotalPDF,"var2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("var2"))->operator[](specBin).c_str()) - 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("var2"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("var2"))->operator[](specBin).c_str()) + 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("var2"))->operator[](specBin).c_str())));
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

  if ((GetVar(TotalPDF,"sigmaMisTag1") != NULL) && (GetVar(TotalPDF,"sigmaMisTag1")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaMisTag1") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaMisTag1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag1"))->operator[](specBin).c_str()) -
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag1"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag1"))->operator[](specBin).c_str()) +
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"sigmaMisTag2") != NULL) && (GetVar(TotalPDF,"sigmaMisTag2")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaMisTag2") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaMisTag2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag2"))->operator[](specBin).c_str()) -
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag2"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag2"))->operator[](specBin).c_str()) +
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"fracMisTag") != NULL) && (GetVar(TotalPDF,"fracMisTag")->getError() != 0.0) && (IsInConstraints(vecConstr,"fracMisTag") == false))
    {
      tmpVar = GetVar(TotalPDF,"fracMisTag");
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

  if ((GetVar(TotalPDF,"nMisTagFrac") != NULL) && (GetVar(TotalPDF,"nMisTagFrac")->getError() != 0.0) && (IsInConstraints(vecConstr,"nMisTagFrac") == false))
    {
      tmpVar = GetVar(TotalPDF,"nMisTagFrac");
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

  myFrame = MyToy->plotNLL(-2000.0,0.0);
  Canv->cd(it++);
  myFrame->Draw();
  Canv->Modified();
  Canv->Update();


  // #############################################
  // # Re-make all the fits and save the results #
  // #############################################
  TH1D* histoDiff1 = new TH1D("histoDiff1","histoDiff1",100,-30.0,30.0);
  histoDiff1->SetFillColor(kAzure+6);
  histoDiff1->SetXTitle("(fit #font[122]{\55} pdf)");
  histoDiff1->SetYTitle("Entries [#]");

  TH1D* histoPull1 = new TH1D("histoPull1","histoPull1",100,-5.0,5.0);
  histoPull1->SetFillColor(kAzure+6);
  histoPull1->SetXTitle("(fit #font[122]{\55} pdf) / #sigma");
  histoPull1->SetYTitle("Entries [#]");

  TH1D* histoChiSq1 = new TH1D("histoChiSq1","histoChiSq1",100,0.0,30.0);
  histoChiSq1->SetFillColor(kAzure+6);
  histoChiSq1->SetXTitle("(fit #font[122]{\55} pdf){^2} / #sigma{^2}");
  histoChiSq1->SetYTitle("Entries [#]");

  TH1D* histoVal1 = new TH1D("histoVal1","histoVal1",100,0.0,400.0);
  histoVal1->SetFillColor(kAzure+6);
  histoVal1->SetXTitle("fit");
  histoVal1->SetYTitle("Entries [#]");

  TH1D* histoNLL1 = new TH1D("histoNLL1","histoNLL1",100,1.0,-1.0);
  histoNLL1->SetFillColor(kGreen-7);
  histoNLL1->SetXTitle("NLL");
  histoNLL1->SetYTitle("Entries [#]");

  TCanvas* cB0Toy = new TCanvas("cB0Toy","cB0Toy", 20, 20, 700, 500);

  cout << "\n[ExtractYield::MakeMassToy]\t@@@ Now fit total TOY invariant mass @@@" << endl;
  RooDataSet* toySample;
  RooFitResult* fitResult;
  double NLLvalue;

  string varName = "nSig";
  for (unsigned int i = 0; i < nToy; i++)
    {
      cout << "\n[ExtractYield::MakeMassToy]\t@@@ Now fitting toy #" << i << " @@@\n" << endl;

      toySample = (RooDataSet*)MyToy->genData(i);
      CopyFitResults(TotalPDF,specBin,fitParam);
      fitResult = MakeMassFit(toySample,&TotalPDF,x,cB0Toy,vecConstr,&NLLvalue,NULL,i);
      if (CheckGoodFit(fitResult) == true) cout << "\n[ExtractYield::MakeMassToy]\t@@@ Fit converged ! @@@" << endl;
      else                                 cout << "\n[ExtractYield::MakeMassToy]\t@@@ Fit didn't converge ! @@@" << endl;


      // ######################################################
      // # To verify the outcome of the RooFit toy-MC studies #
      // ######################################################
      if ((CheckGoodFit(fitResult) == true) && (GetVar(TotalPDF,varName.c_str()) != NULL))
      	{
      	  if (TotalPDF->getVariables()->getRealValue(varName.c_str()) > atof(fitParam->operator[](Utility->GetFitParamIndx(varName.c_str()))->operator[](specBin).c_str()))
	    {
	      histoPull1->Fill((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(fitParam->operator[](Utility->GetFitParamIndx(varName.c_str()))->operator[](specBin).c_str())) /
			       fabs(GetVar(TotalPDF,varName.c_str())->getErrorLo()));
	      histoChiSq1->Fill(pow((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(fitParam->operator[](Utility->GetFitParamIndx(varName.c_str()))->operator[](specBin).c_str())) /
				    fabs(GetVar(TotalPDF,varName.c_str())->getErrorLo()),2.));
	    }
	  else
	    {
	      histoPull1->Fill((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(fitParam->operator[](Utility->GetFitParamIndx(varName.c_str()))->operator[](specBin).c_str())) /
			       fabs(GetVar(TotalPDF,varName.c_str())->getErrorHi()));
	      histoChiSq1->Fill(pow((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(fitParam->operator[](Utility->GetFitParamIndx(varName.c_str()))->operator[](specBin).c_str())) /
				    fabs(GetVar(TotalPDF,varName.c_str())->getErrorHi()),2.));
	    }
	  histoDiff1->Fill(TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(fitParam->operator[](Utility->GetFitParamIndx(varName.c_str()))->operator[](specBin).c_str()));
	  histoVal1->Fill(TotalPDF->getVariables()->getRealValue(varName.c_str()));
	  histoNLL1->Fill(NLLvalue);
      	}
      delete fitResult;
    }


  TCanvas* cNLL1 = new TCanvas("cNLL1","cNLL1",10,10,900,700);
  cNLL1->Divide(2,3);
  cNLL1->cd(1);
  histoDiff1->Draw();
  cNLL1->cd(2);
  histoNLL1->Draw();
  cNLL1->cd(3);
  histoPull1->Draw();
  cNLL1->cd(4);
  histoVal1->Draw();
  cNLL1->cd(5);
  histoChiSq1->Draw();
  cNLL1->Modified();
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
      myString << "_Chi2.root";
      fNLL = new TFile(fileName.replace(fileName.find("_PULL.root"),10,myString.str()).c_str(),"RECREATE");
      fNLL->cd();
      histoChiSq1->Write();
      fNLL->Close();
      delete fNLL;

      myString.clear(); myString.str("");
      myString << "_Val1.root";
      fNLL = new TFile(fileName.replace(fileName.find("_Chi2.root"),10,myString.str()).c_str(),"RECREATE");
      fNLL->cd();
      histoVal1->Write();
      fNLL->Close();
      delete fNLL;

      myString.clear(); myString.str("");
      myString << "_NLL.root";
      fNLL = new TFile(fileName.replace(fileName.find("_Val1.root"),10,myString.str()).c_str(),"RECREATE");
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
				vector<double>* q2Bins,
				unsigned int actualq2BinIndx,
				unsigned int q2BinIndx,
				pair<TF2*,TF2*> effFunc)
// #########################
// # x: mass               #
// # y: angle cos(theta_l) #
// # z: angle cos(theta_K) #
// #########################
{
  // ################################
  // # Read configuration variables #
  // ################################
  unsigned int useSignal = configParam->operator[](Utility->GetConfigParamIndx("SigType"))->operator[](q2BinIndx);
  unsigned int usePeakB  = configParam->operator[](Utility->GetConfigParamIndx("PeakBkgType"))->operator[](q2BinIndx);
  unsigned int useCombB  = configParam->operator[](Utility->GetConfigParamIndx("CombBkgType"))->operator[](q2BinIndx);
  unsigned int useMisTag = configParam->operator[](Utility->GetConfigParamIndx("MistagType"))->operator[](q2BinIndx);


  // ###########################
  // # Read polynomial degrees #
  // ###########################
  unsigned int NCoeffPolyBKGpeak1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP1"))->operator[](q2BinIndx).c_str());
  unsigned int NCoeffPolyBKGcomb1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC1"))->operator[](q2BinIndx).c_str());
  unsigned int NCoeffPolyBKGpeak2 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP2"))->operator[](q2BinIndx).c_str());
  unsigned int NCoeffPolyBKGcomb2 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC2"))->operator[](q2BinIndx).c_str());


  // ###################
  // # Local variables #
  // ###################
  stringstream myString;
  RooArgSet VarsC1, VarsC2;
  RooArgSet VarsP1, VarsP2;
  RooArgSet* VarsAng = new RooArgSet("VarsAng");


  // ################################################
  // # Define mass fit variables and pdf for signal #
  // ################################################
  meanS   = new RooRealVar("meanS","Signal mean",Utility->B0Mass,"GeV");
  
  sigmaS1 = new RooRealVar("sigmaS1","Signal sigma-1",0.0,"GeV");
  MassS1  = new RooGaussian("MassS1","Signal Gaussian-1",*x,*meanS,*sigmaS1);

  sigmaS2 = new RooRealVar("sigmaS2","Signal sigma-2",0.0,"GeV");
  MassS2  = new RooGaussian("MassS2","Signal Gaussian-2",*x,*meanS,*sigmaS2);

  fracMassS = new RooRealVar("fracMassS","Fraction of signal Gaussian",0.0,0.0,1.0);

  if      (useSignal == 1) MassSignal = new RooGaussian(*((RooGaussian*)MassS1),"MassSignal");
  else if (useSignal == 2) MassSignal = new RooAddPdf("MassSignal","Signal mass pdf",RooArgSet(*MassS1,*MassS2),RooArgSet(*fracMassS));
  else                     MassSignal = NULL;


  // ##################################################################
  // # Define angle fit variables and pdf for correctly tagged signal #
  // ##################################################################
  RooArgSet* VarsPolyGT = new RooArgSet("VarsPolyGT");
  AngleS = MakeAngWithEffPDF(effFunc.first,y,z,FitType,useEffPDF,VarsAng,VarsPolyGT,q2Bins,actualq2BinIndx);

  Signal = new RooProdPdf("Signal","Signal Mass*Angle",RooArgSet(*MassSignal,*AngleS));


  // ###################################################################
  // # Define angle fit variables and pdf for combinatorial background #
  // ###################################################################
  for (unsigned int i = 0; i < NCoeffPolyBKGcomb1; i++)
    {
      myString.clear(); myString.str("");
      myString << "c1Poly" << i;
      c1Poly[i] = new RooRealVar(myString.str().c_str(),"Comb.bkg.poly.coef.",0.0);
      VarsC1.add(*c1Poly[i]);
    }
  for (unsigned int i = 0; i < NCoeffPolyBKGcomb2; i++)
    {
      myString.clear(); myString.str("");
      myString << "c2Poly" << i;
      c2Poly[i] = new RooRealVar(myString.str().c_str(),"Comb.bkg.poly.coef.",0.0);
      VarsC2.add(*c2Poly[i]);
    }
  myString.clear(); myString.str("");
  myString << GeneratePolynomial(y,NCoeffPolyBKGcomb1,"c1Poly");
  BkgAngleC1 = new RooGenericPdf("BkgAngleC1",myString.str().c_str(),RooArgSet(*y,VarsC1));
  myString.clear(); myString.str("");
  myString << GeneratePolynomial(z,NCoeffPolyBKGcomb2,"c2Poly");
  BkgAngleC2 = new RooGenericPdf("BkgAngleC2",myString.str().c_str(),RooArgSet(*z,VarsC2));
  BkgAnglesC = new RooProdPdf("BkgAnglesC","Background Angle1*Angle2",RooArgSet(*BkgAngleC1,*BkgAngleC2));


  // ##################################################################
  // # Define mass fit variables and pdf for combinatorial background #
  // ##################################################################
  var1 = new RooRealVar("var1","First background variables",0.0,"GeV");
  myString.clear(); myString.str("");
  myString << "exp(-(" << x->getPlotLabel() << " - meanS) / var1)";
  BkgMassExp1 = new RooGenericPdf("BkgMassExp1",myString.str().c_str(),RooArgSet(*x,*meanS,*var1));

  var2 = new RooRealVar("var2","Second background variable",0.0,"GeV");
  myString.clear(); myString.str("");
  myString << "exp(-(" << x->getPlotLabel() << " - meanS) / var2)";
  BkgMassExp2 = new RooGenericPdf("BkgMassExp2",myString.str().c_str(),RooArgSet(*x,*meanS,*var2));

  fracMassBExp = new RooRealVar("fracMassBExp","Fraction of background Exponential",0.0,0.0,1.0);

  if      (useCombB == 1) BkgMassComb = new RooGenericPdf(*((RooGenericPdf*)BkgMassExp1),"BkgMassComb");
  else if (useCombB == 2) BkgMassComb = new RooAddPdf("BkgMassComb","Background mass comb. bkg pdf",RooArgSet(*BkgMassExp1,*BkgMassExp2),RooArgSet(*fracMassBExp));
  else if (useCombB == 3)
    {
      myString.clear(); myString.str("");
      myString << "TMath::Erfc((" << x->getPlotLabel() << " - var1) / sqrt(2*var2*var2))";
      BkgMassComb = new RooGenericPdf("BkgMass",myString.str().c_str(),RooArgSet(*x,*var1,*var2));
    }
  else                    BkgMassComb = NULL;

  BkgMassAngleComb = new RooProdPdf("BkgMassAngleComb","Combinatorial bkg Mass*Angle",RooArgSet(*BkgMassComb,*BkgAnglesC));

 
  // #######################################################
  // # Define mass fit variables and pdf for mistag signal #
  // #######################################################
  sigmaMisTag1 = new RooRealVar("sigmaMisTag1","Mistag sigma-1",0.0,"GeV");
  MassMisTag1  = new RooGaussian("MassMisTag1","Mistag-1",*x,*meanS,*sigmaMisTag1);

  sigmaMisTag2 = new RooRealVar("sigmaMisTag2","Mistag sigma-2",0.0,"GeV");
  MassMisTag2  = new RooGaussian("MassMisTag2","Mistag-2",*x,*meanS,*sigmaMisTag2);

  fracMisTag = new RooRealVar("fracMisTag","Fraction mistag Gaussian",0.0,0.0,1.0);

  if      (useMisTag == 1) MassMisTag = new RooGaussian(*((RooGaussian*)MassMisTag1),"MassMisTag");
  else if (useMisTag == 2) MassMisTag = new RooAddPdf("MassMisTag","Mistag mass pdf",RooArgSet(*MassMisTag1,*MassMisTag2),RooArgSet(*fracMisTag));
  else                     MassMisTag = NULL;


  // ############################################################
  // # Define angle fit variables and pdf for mis-tagged signal #
  // ############################################################
  RooArgSet* VarsPolyMT = new RooArgSet("VarsPolyMT");
  AngleMisTag = MakeAngWithEffPDF(effFunc.second,y,z,FitType*10,useEffPDF,VarsAng,VarsPolyMT,q2Bins,actualq2BinIndx);

  MassAngleMisTag = new RooProdPdf("MassAngleMisTag","Mistag bkg Mass*Angle",RooArgSet(*MassMisTag,*AngleMisTag));


  // #############################################################
  // # Define angle fit variables and pdf for peaking background #
  // #############################################################
  for (unsigned int i = 0; i < NCoeffPolyBKGpeak1; i++)
    {
      myString.clear(); myString.str("");
      myString << "p1Poly" << i;
      p1Poly[i] = new RooRealVar(myString.str().c_str(),"Peak.bkg.poly.coef.",0.0);
      VarsP1.add(*p1Poly[i]);
    }
  for (unsigned int i = 0; i < NCoeffPolyBKGpeak2; i++)
    {
      myString.clear(); myString.str("");
      myString << "p2Poly" << i;
      p2Poly[i] = new RooRealVar(myString.str().c_str(),"Peak.bkg.poly.coef.",0.0);
      VarsP2.add(*p2Poly[i]);
    }
  myString.clear(); myString.str("");
  myString << GeneratePolynomial(y,NCoeffPolyBKGpeak1,"p1Poly");
  BkgAngleP1 = new RooGenericPdf("BkgAngleP1",myString.str().c_str(),RooArgSet(*y,VarsP1));
  myString.clear(); myString.str("");
  myString << GeneratePolynomial(z,NCoeffPolyBKGpeak2,"p2Poly");
  BkgAngleP2 = new RooGenericPdf("BkgAngleP2",myString.str().c_str(),RooArgSet(*z,VarsP2));
  BkgAnglesP = new RooProdPdf("BkgAnglesP","Background Angle1*Angle2",RooArgSet(*BkgAngleP1,*BkgAngleP2));


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
  BkgMassRPeak   = new RooAddPdf("BkgMassRPeak","Right peaking bkg mass pdf",RooArgSet(*BkgMassRPeak1,*BkgMassRPeak2),RooArgSet(*fracMassBRPeak));

  meanL1         = new RooRealVar("meanL1","Bkg left peak mean-1",0.0,"GeV");
  sigmaL1        = new RooRealVar("sigmaL1","Bkg left peak sigma-1",0.0,"GeV");
  BkgMassLPeak1  = new RooGaussian("BkgMassLPeak1","Bkg left peak-1",*x,*meanL1,*sigmaL1);

  meanL2         = new RooRealVar("meanL2","Bkg left peak mean-2",0.0,"GeV");
  sigmaL2        = new RooRealVar("sigmaL2","Bkg left peak sigma-2",0.0,"GeV");
  BkgMassLPeak2  = new RooGaussian("BkgMassLPeak2","Bkg left peak-2",*x,*meanL2,*sigmaL2);

  fracMassBLPeak = new RooRealVar("fracMassBLPeak","Fraction of background left Peak",0.0,0.0,1.0);
  BkgMassLPeak   = new RooAddPdf("BkgMassLPeak","Left peaking bkg mass pdf",RooArgSet(*BkgMassLPeak1,*BkgMassLPeak2),RooArgSet(*fracMassBLPeak));

  fracMassBPeak  = new RooRealVar("fracMassBPeak","Fraction of background right-left Peak",0.0,0.0,1.0);
  
  if      (usePeakB == 1)  BkgMassPeak = new RooGaussian(*((RooGaussian*)BkgMassRPeak1),"BkgMassPeak");
  else if (usePeakB == 2)  BkgMassPeak = new RooAddPdf(*((RooAddPdf*)BkgMassRPeak),"BkgMassPeak");
  else if (usePeakB == 11) BkgMassPeak = new RooAddPdf("BkgMassPeak","Peaking bkg mass pdf",RooArgSet(*BkgMassRPeak1,*BkgMassLPeak1),RooArgSet(*fracMassBPeak));
  else if (usePeakB == 12) BkgMassPeak = new RooAddPdf("BkgMassPeak","Peaking bkg mass pdf",RooArgSet(*BkgMassRPeak,*BkgMassLPeak),RooArgSet(*fracMassBPeak));
  else                     BkgMassPeak = NULL;

  BkgMassAnglePeak = new RooProdPdf("BkgMassAnglePeak","Peaking bkg Mass*Angle",RooArgSet(*BkgMassPeak,*BkgAnglesP));


  // ###########################
  // # Define pdf coefficients #
  // ###########################
  nSig     = new RooRealVar("nSig","Number of signal events",1.0);
  nBkgComb = new RooRealVar("nBkgComb","Number of comb. background events",1.0);
  RooFormulaVar* nMisTag;
  if ((useSignal == 0) && (usePeakB == 0) && (useCombB == 0) && (useMisTag != 0))
    {
      nMisTagFrac = new RooRealVar("nMisTagFrac","Number of mistag events",1.0);
      nMisTag     = new RooFormulaVar("nMisTag","nMisTagFrac", RooArgSet(*nMisTagFrac));
    }
  else
    {
      nMisTagFrac = new RooRealVar("nMisTagFrac","Fraction of mistag",0.0,0.0,1.0);
      nMisTag     = new RooFormulaVar("nMisTag","nSig * nMisTagFrac / (1 - nMisTagFrac)", RooArgSet(*nSig,*nMisTagFrac));
    }
  nBkgPeak = new RooRealVar("nBkgPeak","Number of peaking background events",1.0);


  if ((FitType == 36) || (FitType == 56) || (FitType == 76)) *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*AngleS),RooArgSet(*nSig));

  else if ((useSignal != 0) && (usePeakB == 0) && (useCombB == 0) && (useMisTag == 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*Signal),RooArgSet(*nSig));

  else if ((useSignal == 0) && (usePeakB != 0) && (useCombB == 0) && (useMisTag == 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*BkgMassAnglePeak),RooArgSet(*nBkgPeak));

  else if ((useSignal == 0) && (usePeakB == 0) && (useCombB != 0) && (useMisTag == 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*BkgMassAngleComb),RooArgSet(*nBkgComb));

  else if ((useSignal == 0) && (usePeakB == 0) && (useCombB == 0) && (useMisTag != 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*MassAngleMisTag),RooArgSet(*nMisTag));

  else if ((useSignal != 0) && (usePeakB == 0) && (useCombB != 0) && (useMisTag == 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*Signal,*BkgMassAngleComb),RooArgSet(*nSig,*nBkgComb));

  else if ((useSignal != 0) && (usePeakB == 0) && (useCombB == 0) && (useMisTag != 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*Signal,*MassAngleMisTag),RooArgSet(*nSig,*nMisTag));
  
  else if ((useSignal != 0) && (usePeakB != 0) && (useCombB != 0) && (useMisTag == 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*Signal,*BkgMassAngleComb,*BkgMassAnglePeak),RooArgSet(*nSig,*nBkgComb,*nBkgPeak));

  else if ((useSignal != 0) && (usePeakB == 0) && (useCombB != 0) && (useMisTag != 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*Signal,*BkgMassAngleComb,*MassAngleMisTag),RooArgSet(*nSig,*nBkgComb,*nMisTag));

  else if ((useSignal != 0) && (usePeakB != 0) && (useCombB != 0) && (useMisTag != 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*Signal,*BkgMassAngleComb,*MassAngleMisTag,*BkgMassAnglePeak),RooArgSet(*nSig,*nBkgComb,*nMisTag,*nBkgPeak));

  else
    {
      cout << "[ExtractYield::InstantiateMass2AnglesFit]\tIncorrect configuration sequence : useSignal = " << useSignal << "\tusePeakB = " << usePeakB << "\tuseCombB = " << useCombB << "\tuseMisTag = " << useMisTag << endl;
      exit (EXIT_FAILURE);
    }
}


RooFitResult* MakeMass2AnglesFit (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* x, RooRealVar* y, RooRealVar* z, TCanvas* Canv, unsigned int FitType, RooArgSet* vecConstr, double* NLLvalue, TPaveText* extText, unsigned int ID)
{
  // ###################
  // # Local variables #
  // ###################
  RooFitResult* fitResult = NULL;
  RooArgSet VarsYield;
  stringstream myString;
  double signalSigmaGoodT  = 0.0;
  double signalSigmaGoodTE = 0.0;
  double signalSigmaMisT   = 0.0;
  double signalSigmaMisTE  = 0.0;
  double signalSigmaT      = 0.0;
  double signalSigmaTE     = 0.0;
  double totalYield        = 0.0;
  double totalYieldE       = 0.0;
  unsigned int nElements   = 0;
  TString legNames[6];

  const unsigned int nCanv = 6;
  TCanvas* localCanv[nCanv];
  for (unsigned int i = 0; i < nCanv; i++)
    {
      myString.clear(); myString.str("");
      myString << "localCanv" << i;
      localCanv[i] = new TCanvas(myString.str().c_str(),myString.str().c_str(),20,20,700,500);
    }


  if ((FitType == 36) || (FitType == 56) || (FitType == 76))
    {
      // ###################
      // # Make actual fit #
      // ###################
      if (atoi(Utility->GetGenericParam("ApplyConstr").c_str()) == true) fitResult = (*TotalPDF)->fitTo(*dataSet,Extended(true),ExternalConstraints(*vecConstr),Save(true),Minos(USEMINOS));
      else                                                               fitResult = (*TotalPDF)->fitTo(*dataSet,Extended(true),Save(true),Minos(USEMINOS));


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
      legNames[nElements++] = "Data";
      
      if (FUNCERRBAND == true) (*TotalPDF)->plotOn(myFrameY,Name((*TotalPDF)->getPlotLabel()),VisualizeError(*fitResult,1,true),VLines(),FillColor(kGreen-7),Project(*z));
      else                     (*TotalPDF)->plotOn(myFrameY,Name((*TotalPDF)->getPlotLabel()),Project(*z));
      legNames[nElements++] = "Total p.d.f.";

      (*TotalPDF)->paramOn(myFrameY,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(RooArgSet(*nSig)));


      myString.clear(); myString.str("");
      myString << (*TotalPDF)->getPlotLabel() << "_paramBox";
      TPaveText* paveTextY = (TPaveText*)myFrameY->findObject(myString.str().c_str());
      paveTextY->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));

      TLegend* legY = new TLegend(0.75, 0.65, 0.97, 0.88, "");
      for (unsigned int i = 0; i < nElements; i++)
	{
	  TString objName = myFrameY->nameOf(i);
	  if (objName == "") continue;
	  TObject* obj = myFrameY->findObject(objName.Data());
	  legY->AddEntry(obj,legNames[i],"lp");
	}
      legY->SetFillStyle(0);
      legY->SetFillColor(0);
      legY->SetTextSize(0.04);
      legY->SetBorderSize(0);
 

      // ##########################
      // # Angular-2 plot results #
      // ##########################
      Canv->cd(3);
      RooPlot* myFrameZ = z->frame(NBINS);

      dataSet->plotOn(myFrameZ,Name(MakeName(dataSet,ID).c_str()));

      if (FUNCERRBAND == true) (*TotalPDF)->plotOn(myFrameZ,Name((*TotalPDF)->getPlotLabel()),VisualizeError(*fitResult,1,true),VLines(),FillColor(kGreen-7),Project(*y));
      else                     (*TotalPDF)->plotOn(myFrameZ,Name((*TotalPDF)->getPlotLabel()),Project(*y));

      (*TotalPDF)->paramOn(myFrameZ,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(RooArgSet(*nSig)));
      

      myString.clear(); myString.str("");
      myString << (*TotalPDF)->getPlotLabel() << "_paramBox";
      TPaveText* paveTextZ = (TPaveText*)myFrameZ->findObject(myString.str().c_str());
      paveTextZ->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));

      TLegend* legZ = new TLegend(0.75, 0.65, 0.97, 0.88, "");
      for (unsigned int i = 0; i < nElements; i++)
	{
	  TString objName = myFrameZ->nameOf(i);
	  if (objName == "") continue;
	  TObject* obj = myFrameZ->findObject(objName.Data());
	  legZ->AddEntry(obj,legNames[i],"lp");
	}
      legZ->SetFillStyle(0);
      legZ->SetFillColor(0);
      legZ->SetTextSize(0.04);
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
      paveTextY->SetTextSize(0.03);
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
      paveTextZ->SetTextSize(0.03);
      paveTextZ->Paint();
      myFrameZ->Draw();
      legZ->Draw("same");
    }
  else
    {
      // #####################
      // # Make sideband fit #
      // #####################
      if (GetVar(*TotalPDF,"nBkgComb") != NULL)
      	{
	  cout << "[ExtractYield::MakeMass2AnglesFit]\t@@@ Making comb. angular background sideband fit @@@" << endl;

      	  RooAbsPdf* TmpPDF = NULL;
      	  RooDataSet* sideBands = NULL;
      	  RooRealVar frac("frac","Fraction",0.5,0.0,1.0);
      	  RooArgSet constrSidebads;
      	  ClearVars(&constrSidebads);
      	  BuildAngularConstraints(&constrSidebads,*TotalPDF,"peak");


      	  // ################
      	  // # Save results #
      	  // ################
      	  fileFitResults << "====================================================================" << endl;
      	  fileFitResults << "@@@@@@ B0 mass sideband fit @@@@@@" << endl;
	  fileFitResults << "Amplitude of signal region (+/- n*< Sigma >): " << Utility->GetGenericParam("NSigmaB0").c_str() << " * " << Utility->GetB0Width() << endl;


      	  // ##############
      	  // # Get p.d.f. #
      	  // ##############
      	  if (GetVar(*TotalPDF,"nBkgPeak") != NULL) TmpPDF = new RooAddPdf("TmpPDF","Temporary p.d.f.",RooArgSet(*BkgAnglesC,*BkgAnglesP),RooArgSet(frac));
      	  else                                      TmpPDF = new RooProdPdf(*((RooProdPdf*)BkgAnglesC),"TmpPDF");


      	  // #############
      	  // # Sidebands #
      	  // #############
      	  myString.clear(); myString.str("");
	  myString << "B0MassArb < " << (*TotalPDF)->getVariables()->getRealValue("meanS") - atof(Utility->GetGenericParam("NSigmaB0").c_str())*Utility->GetB0Width();
	  myString << " || B0MassArb > " << (*TotalPDF)->getVariables()->getRealValue("meanS") + atof(Utility->GetGenericParam("NSigmaB0").c_str())*Utility->GetB0Width();
      	  cout << "[ExtractYield::MakeMass2AnglesFit]\tCut for B0 sidebands: " << myString.str().c_str() << endl;
      	  sideBands = (RooDataSet*)dataSet->reduce(myString.str().c_str());


      	  // ###################
      	  // # Make actual fit #
      	  // ###################
      	  if (atoi(Utility->GetGenericParam("ApplyConstr").c_str()) == true) fitResult = TmpPDF->fitTo(*sideBands,ExternalConstraints(constrSidebads),Save(true));
      	  else                                                               fitResult = TmpPDF->fitTo(*sideBands,Save(true));
	  if (fitResult != NULL) fitResult->Print("v");


      	  // ####################
      	  // # Save fit results #
      	  // ####################
      	  StorePolyResultsInFile(TotalPDF);


      	  delete TmpPDF;
      	  delete sideBands;
      	  ClearVars(&constrSidebads);


      	  // ################################
      	  // # Fix comb. angular background #
      	  // ################################
	  cout << "[ExtractYield::MakeMass2AnglesFit]\t@@@ Fixing comb. angular background after sideband fit @@@" << endl;
	  for (unsigned int i = 0 ; i < NCOEFFPOLYBKG; i++)
            {
              myString.clear(); myString.str("");
              myString << "c1Poly" << i;
              if (GetVar(*TotalPDF,myString.str().c_str()) != NULL) GetVar(*TotalPDF,myString.str().c_str())->setConstant(true);
	      
              myString.clear(); myString.str("");
              myString << "c2Poly" << i;
              if (GetVar(*TotalPDF,myString.str().c_str()) != NULL) GetVar(*TotalPDF,myString.str().c_str())->setConstant(true);
            }


	  if (CheckGoodFit(fitResult) == true) { cout << "\n[ExtractYield::MakeMass2AnglesFit]\t@@@ Fit converged ! @@@" << endl;       delete fitResult; }
	  else                                 { cout << "\n[ExtractYield::MakeMass2AnglesFit]\t@@@ Fit didn't converge ! @@@" << endl; return fitResult; }
	}


      // ###################
      // # Make actual fit #
      // ###################
      if (atoi(Utility->GetGenericParam("ApplyConstr").c_str()) == true) fitResult = (*TotalPDF)->fitTo(*dataSet,Extended(true),ExternalConstraints(*vecConstr),Save(true),Minos(USEMINOS));
      else                                                               fitResult = (*TotalPDF)->fitTo(*dataSet,Extended(true),Save(true),Minos(USEMINOS));


      // ##################################################
      // # Set p.d.f independent variables to known point #
      // ##################################################
      if (GetVar(*TotalPDF,x->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(x->getPlotLabel(),Utility->B0Mass);
      if (GetVar(*TotalPDF,y->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(y->getPlotLabel(),0.0);
      if (GetVar(*TotalPDF,z->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(z->getPlotLabel(),0.0);
      if (fitResult != NULL) fitResult->Print("v");


      // #####################
      // # Mass plot results #
      // #####################
      Canv->cd(1);
      RooPlot* myFrameX = x->frame(NBINS);

      dataSet->plotOn(myFrameX,Name(MakeName(dataSet,ID).c_str()));
      legNames[nElements++] = "Data";

      if (FUNCERRBAND == true) (*TotalPDF)->plotOn(myFrameX,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),VisualizeError(*fitResult,1,true),VLines(),FillColor(kGreen-7),Project(RooArgSet(*y,*z)));
      else                     (*TotalPDF)->plotOn(myFrameX,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),Project(RooArgSet(*y,*z)));
      legNames[nElements++] = "Total p.d.f.";

      if (GetVar(*TotalPDF,"nSig") != NULL)
	{
	  (*TotalPDF)->plotOn(myFrameX, Components(*Signal), LineStyle(7), LineColor(kBlue), Project(RooArgSet(*y,*z)));
	  legNames[nElements++] = "Right-tag sig";
	  VarsYield.add(*nSig);
	}

      if (GetVar(*TotalPDF,"nMisTagFrac") != NULL)
	{
	  (*TotalPDF)->plotOn(myFrameX, Components(*MassAngleMisTag), LineStyle(8), LineColor(kAzure+6), Project(RooArgSet(*y,*z)));
	  legNames[nElements++] = "Mis-tag sig";
	  VarsYield.add(*nMisTagFrac);
	}

      if (GetVar(*TotalPDF,"nBkgComb") != NULL)
	{
	  (*TotalPDF)->plotOn(myFrameX, Components(*BkgMassAngleComb), LineStyle(4), LineColor(kRed), Project(RooArgSet(*y,*z)));
	  legNames[nElements++] = "Comb. bkg";
	  VarsYield.add(*nBkgComb);
	}

      if (GetVar(*TotalPDF,"nBkgPeak") != NULL)
	{
	  (*TotalPDF)->plotOn(myFrameX, Components(*BkgMassAnglePeak), LineStyle(3), LineColor(kViolet), Project(RooArgSet(*y,*z)));
	  legNames[nElements++] = "Peak. bkg";
	  VarsYield.add(*nBkgPeak);
	}

      (*TotalPDF)->paramOn(myFrameX,Format("NEU",AutoPrecision(2)),Layout(0.11,0.42,0.88),Parameters(VarsYield),ShowConstants(true));
  
      
      myString.clear(); myString.str("");
      myString << (*TotalPDF)->getPlotLabel() << "_paramBox";
      TPaveText* paveTextX = (TPaveText*)myFrameX->findObject(myString.str().c_str());

      if (GetVar(*TotalPDF,"nSig") != NULL)
	{
	  if (GetVar(*TotalPDF,"fracMassS") != NULL)
	    {
              signalSigmaGoodT  = sqrt( GetVar(*TotalPDF,"fracMassS")->getVal() * pow(GetVar(*TotalPDF,"sigmaS1")->getVal(),2.) + (1. - GetVar(*TotalPDF,"fracMassS")->getVal()) * pow(GetVar(*TotalPDF,"sigmaS2")->getVal(),2.) );
              signalSigmaGoodTE = 1./(2.*signalSigmaGoodT) * sqrt( pow((pow(GetVar(*TotalPDF,"sigmaS1")->getVal(),2.) - pow(GetVar(*TotalPDF,"sigmaS2")->getVal(),2.)) * GetVar(*TotalPDF,"fracMassS")->getError(),2.) +
                                                                   pow(2. * GetVar(*TotalPDF,"fracMassS")->getVal()        * GetVar(*TotalPDF,"sigmaS1")->getVal() * GetVar(*TotalPDF,"sigmaS1")->getError(),2.) +
                                                                   pow(2. * (1. - GetVar(*TotalPDF,"fracMassS")->getVal()) * GetVar(*TotalPDF,"sigmaS2")->getVal() * GetVar(*TotalPDF,"sigmaS2")->getError(),2.) );
	    }
	  else if (GetVar(*TotalPDF,"sigmaS1") != NULL)
	    {
	      signalSigmaGoodT  = GetVar(*TotalPDF,"sigmaS1")->getVal();
	      signalSigmaGoodTE = GetVar(*TotalPDF,"sigmaS1")->getError();
	    }
	  
	  signalSigmaT  = signalSigmaGoodT;
	  signalSigmaTE = signalSigmaGoodTE;
	}
      
      if (GetVar(*TotalPDF,"nMisTagFrac") != NULL)
	{
	  if (GetVar(*TotalPDF,"fracMisTag") != NULL)
	    {
	      signalSigmaMisT  = sqrt( GetVar(*TotalPDF,"fracMisTag")->getVal() * pow(GetVar(*TotalPDF,"sigmaMisTag1")->getVal(),2.) + (1. - GetVar(*TotalPDF,"fracMisTag")->getVal()) * pow(GetVar(*TotalPDF,"sigmaMisTag2")->getVal(),2.) );
              signalSigmaMisTE = 1./(2.*signalSigmaMisT) * sqrt( pow((pow(GetVar(*TotalPDF,"sigmaMisTag1")->getVal(),2.) - pow(GetVar(*TotalPDF,"sigmaMisTag2")->getVal(),2.)) * GetVar(*TotalPDF,"fracMisTag")->getError(),2.) +
                                                                 pow(2. * GetVar(*TotalPDF,"fracMisTag")->getVal()        * GetVar(*TotalPDF,"sigmaMisTag1")->getVal() * GetVar(*TotalPDF,"sigmaMisTag1")->getError(),2.) +
                                                                 pow(2. * (1. - GetVar(*TotalPDF,"fracMisTag")->getVal()) * GetVar(*TotalPDF,"sigmaMisTag2")->getVal() * GetVar(*TotalPDF,"sigmaMisTag2")->getError(),2.) );
	    }
	  else if (GetVar(*TotalPDF,"sigmaMisTag1") != NULL)
	    {
	      signalSigmaMisT  = GetVar(*TotalPDF,"sigmaMisTag1")->getVal();
	      signalSigmaMisTE = GetVar(*TotalPDF,"sigmaMisTag1")->getError();
	    }
	  
	  signalSigmaT  = signalSigmaMisT;
	  signalSigmaTE = signalSigmaMisTE;
	}
      
      if ((GetVar(*TotalPDF,"nSig") != NULL) && (GetVar(*TotalPDF,"nMisTagFrac") != NULL))
	{
	  totalYield  = GetVar(*TotalPDF,"nSig")->getVal() / (1. - GetVar(*TotalPDF,"nMisTagFrac")->getVal());
	  totalYieldE = totalYield * sqrt( pow(GetVar(*TotalPDF,"nSig")->getError()/GetVar(*TotalPDF,"nSig")->getVal(),2.) + pow(GetVar(*TotalPDF,"nMisTagFrac")->getError()/(1. - GetVar(*TotalPDF,"nMisTagFrac")->getVal()),2.) );
	  
	  double fraction  = (1. - 2.*GetVar(*TotalPDF,"nMisTagFrac")->getVal()) / (1. - GetVar(*TotalPDF,"nMisTagFrac")->getVal());
	  double fractionE = GetVar(*TotalPDF,"nMisTagFrac")->getError() / pow((1. - GetVar(*TotalPDF,"nMisTagFrac")->getVal()),2.);
	  
	  signalSigmaT  = sqrt( fraction * pow(signalSigmaGoodT,2.) + (1. - fraction) * pow(signalSigmaMisT,2.) );
	  signalSigmaTE = 1./(2.*signalSigmaT) * sqrt( pow((pow(signalSigmaGoodT,2.) - pow(signalSigmaMisT,2.)) * fractionE,2.) +
						       pow(2. * fraction        * signalSigmaGoodT * signalSigmaGoodTE,2.) +
						       pow(2. * (1. - fraction) * signalSigmaMisT  * signalSigmaMisTE,2.) );
	  
	  paveTextX->AddText(Form("%s%1.f#pm%1.f","Total sig yield = ",totalYield,totalYieldE));
	}
      else if (GetVar(*TotalPDF,"nSig") != NULL)
	{
	  totalYield  = GetVar(*TotalPDF,"nSig")->getVal();
	  totalYieldE = GetVar(*TotalPDF,"nSig")->getError();
	}
      
      if ((GetVar(*TotalPDF,"nSig") != NULL) || (GetVar(*TotalPDF,"nMisTagFrac") != NULL))
	{
	  paveTextX->AddText(Form("%s%.3f#pm%.3f","#mu = ",GetVar(*TotalPDF,"meanS")->getVal(),GetVar(*TotalPDF,"meanS")->getError()));
	  paveTextX->AddText(Form("%s%.3f#pm%.3f","Total < #sigma > = ",signalSigmaT,signalSigmaTE));
	}
      paveTextX->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameX->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
      
      TLegend* legX = new TLegend(0.75, 0.65, 0.97, 0.88, "");
      for (unsigned int i = 0; i < nElements; i++)
  	{
  	  TString objName = myFrameX->nameOf(i);
  	  if (objName == "") continue;
  	  TObject* obj = myFrameX->findObject(objName.Data());
	  legX->AddEntry(obj,legNames[i],"lp");
  	}
      legX->SetFillStyle(0);
      legX->SetFillColor(0);
      legX->SetTextSize(0.04);
      legX->SetBorderSize(0);


      // ##########################
      // # Angular-1 plot results #
      // ##########################
      Canv->cd(2);
      RooPlot* myFrameY = y->frame(NBINS);

      dataSet->plotOn(myFrameY,Name(MakeName(dataSet,ID).c_str()));

      if (FUNCERRBAND == true) (*TotalPDF)->plotOn(myFrameY,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),VisualizeError(*fitResult,1,true),VLines(),FillColor(kGreen-7),Project(RooArgSet(*x,*z)));
      else                     (*TotalPDF)->plotOn(myFrameY,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),Project(RooArgSet(*x,*z)));

      if (GetVar(*TotalPDF,"nSig")        != NULL) (*TotalPDF)->plotOn(myFrameY, Components(*Signal),           LineStyle(7), LineColor(kBlue),    Project(RooArgSet(*x,*z)));
      if (GetVar(*TotalPDF,"nMisTagFrac") != NULL) (*TotalPDF)->plotOn(myFrameY, Components(*MassAngleMisTag),  LineStyle(8), LineColor(kAzure+6), Project(RooArgSet(*x,*z)));
      if (GetVar(*TotalPDF,"nBkgComb")    != NULL) (*TotalPDF)->plotOn(myFrameY, Components(*BkgMassAngleComb), LineStyle(4), LineColor(kRed),     Project(RooArgSet(*x,*z)));
      if (GetVar(*TotalPDF,"nBkgPeak")    != NULL) (*TotalPDF)->plotOn(myFrameY, Components(*BkgMassAnglePeak), LineStyle(3), LineColor(kViolet),  Project(RooArgSet(*x,*z)));


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

      TLegend* legY = new TLegend(0.75, 0.65, 0.97, 0.88, "");
      for (unsigned int i = 0; i < nElements; i++)
      	{
      	  TString objName = myFrameY->nameOf(i);
      	  if (objName == "") continue;
      	  TObject* obj = myFrameY->findObject(objName.Data());
      	  legY->AddEntry(obj,legNames[i],"lp");
      	}
      legY->SetFillStyle(0);
      legY->SetFillColor(0);
      legY->SetTextSize(0.04);
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

      if (FUNCERRBAND == true) (*TotalPDF)->plotOn(myFrameZ,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),VisualizeError(*fitResult,1,true),VLines(),FillColor(kGreen-7),Project(RooArgSet(*x,*y)));
      else                     (*TotalPDF)->plotOn(myFrameZ,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),Project(RooArgSet(*x,*y)));

      if (GetVar(*TotalPDF,"nSig")        != NULL) (*TotalPDF)->plotOn(myFrameZ, Components(*Signal),           LineStyle(7), LineColor(kBlue),    Project(RooArgSet(*x,*y)));
      if (GetVar(*TotalPDF,"nMisTagFrac") != NULL) (*TotalPDF)->plotOn(myFrameZ, Components(*MassAngleMisTag),  LineStyle(8), LineColor(kAzure+6), Project(RooArgSet(*x,*y)));
      if (GetVar(*TotalPDF,"nBkgComb")    != NULL) (*TotalPDF)->plotOn(myFrameZ, Components(*BkgMassAngleComb), LineStyle(4), LineColor(kRed),     Project(RooArgSet(*x,*y)));
      if (GetVar(*TotalPDF,"nBkgPeak")    != NULL) (*TotalPDF)->plotOn(myFrameZ, Components(*BkgMassAnglePeak), LineStyle(3), LineColor(kViolet),  Project(RooArgSet(*x,*y)));


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

      TLegend* legZ = new TLegend(0.75, 0.65, 0.97, 0.88, "");
      for (unsigned int i = 0; i < nElements; i++)
      	{
      	  TString objName = myFrameZ->nameOf(i);
      	  if (objName == "") continue;
      	  TObject* obj = myFrameZ->findObject(objName.Data());
      	  legZ->AddEntry(obj,legNames[i],"lp");
      	}
      legZ->SetFillStyle(0);
      legZ->SetFillColor(0);
      legZ->SetTextSize(0.04);
      legZ->SetBorderSize(0);
      legZ->Draw("same");


      // #####################
      // # Fill local canvas #
      // #####################
      localCanv[2]->cd();
      myFrameZ->Draw();
      legZ->Draw("same");


      // ################################
      // # LogLikelihood profiling scan #
      // ################################
      // @TMP@
      // if ((FitType != 26) && (GetVar(*TotalPDF,"FlS") != NULL) && (GetVar(*TotalPDF,"AfbS") != NULL))
      	// {
      	//   // #############################
      	//   // # Turn off all the printout #
      	//   // #############################
      	//   RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);


      	//   RooAbsReal* NLL;
      	//   if (atoi(Utility->GetGenericParam("ApplyConstr").c_str()) == true) NLL = (*TotalPDF)->createNLL(*dataSet,Extended(true),ExternalConstraints(*vecConstr));
      	//   else                                                               NLL = (*TotalPDF)->createNLL(*dataSet,Extended(true));


      	//   localCanv[3]->cd();
      	//   RooMinuit RooMin(*NLL);
      	//   RooPlot* myFrameNLL = RooMin.contour(*GetVar(*TotalPDF,"AfbS"),*GetVar(*TotalPDF,"FlS"),1.0,2.0,3.0);
      	//   DrawString(LUMI,myFrameY);
      	//   myFrameNLL->Draw();


      	//   localCanv[4]->cd();
        //   RooPlot* myFrameNLLVar1 = GetVar(*TotalPDF,"FlS")->frame();
      	//   NLL->plotOn(myFrameNLLVar1,ShiftToZero());
      	//   RooAbsReal* var1Profile = NLL->createProfile(*(GetVar(*TotalPDF,"FlS")));
      	//   var1Profile->plotOn(myFrameNLLVar1,LineColor(kRed));
      	//   myFrameNLLVar1->SetMinimum(0);
      	//   myFrameNLLVar1->SetMaximum(5);
        //   myFrameNLLVar1->Draw();


      	//   localCanv[5]->cd();
        //   RooPlot* myFrameNLLVar2 = GetVar(*TotalPDF,"AfbS")->frame();
      	//   NLL->plotOn(myFrameNLLVar2,ShiftToZero());
      	//   RooAbsReal* var2Profile = NLL->createProfile(*(GetVar(*TotalPDF,"AfbS")));
      	//   var2Profile->plotOn(myFrameNLLVar2,LineColor(kRed));
      	//   myFrameNLLVar2->SetMinimum(0);
      	//   myFrameNLLVar2->SetMaximum(5);
        //   myFrameNLLVar2->Draw();


      	//   delete NLL;


      	//   // ############################
      	//   // # Turn on all the printout #
      	//   // ############################
      	//   RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
      	// }


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
      	  if (GetVar(*TotalPDF,"nSig") != NULL) extText->AddText(Form("%s%1.f#pm%1.f","Signal yield: ",totalYield,totalYieldE));
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

      // ################
      // # Save results #
      // ################
      fileFitResults << "@@@@@@ B0 mass sideband fit overlaps @@@@@@" << endl;
      fileFitResults << "Amplitude of signal region: +/- " << Utility->GetGenericParam("NSigmaB0").c_str() << " * < Sigma >" << endl;


      // ################
      // # Low sideband #
      // ################
      x->setRange("lowSideband",Utility->B0Mass - atof(Utility->GetGenericParam("B0MassIntervalLeft").c_str()),
		  (*TotalPDF)->getVariables()->getRealValue("meanS") - atof(Utility->GetGenericParam("NSigmaB0").c_str())*signalSigmaT);
      y->setRange("lowSideband",y->getMin(),y->getMax());
      z->setRange("lowSideband",z->getMin(),z->getMax());


      // #################
      // # Signal region #
      // #################
      x->setRange("signalRegion",(*TotalPDF)->getVariables()->getRealValue("meanS") - atof(Utility->GetGenericParam("NSigmaB0").c_str())*signalSigmaT,
		  (*TotalPDF)->getVariables()->getRealValue("meanS") + atof(Utility->GetGenericParam("NSigmaB0").c_str())*signalSigmaT);
      y->setRange("signalRegion",y->getMin(),y->getMax());
      z->setRange("signalRegion",z->getMin(),z->getMax());


      // #################
      // # High sideband #
      // #################
      x->setRange("highSideband",(*TotalPDF)->getVariables()->getRealValue("meanS") + atof(Utility->GetGenericParam("NSigmaB0").c_str())*signalSigmaT,
		  Utility->B0Mass + atof(Utility->GetGenericParam("B0MassIntervalRight").c_str()));
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

      (*TotalPDF)->plotOn(myFrameLowSideBY,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),ProjectionRange("lowSideband"),Project(RooArgSet(*x,*z)));

      if (GetVar(*TotalPDF,"nSig")        != NULL) (*TotalPDF)->plotOn(myFrameLowSideBY, Components(*Signal),           LineStyle(7), LineColor(kBlue),    ProjectionRange("lowSideband"), Project(RooArgSet(*x,*z)));
      if (GetVar(*TotalPDF,"nMisTagFrac") != NULL) (*TotalPDF)->plotOn(myFrameLowSideBY, Components(*MassAngleMisTag),  LineStyle(8), LineColor(kAzure+6), ProjectionRange("lowSideband"), Project(RooArgSet(*x,*z)));
      if (GetVar(*TotalPDF,"nBkgComb")    != NULL) (*TotalPDF)->plotOn(myFrameLowSideBY, Components(*BkgMassAngleComb), LineStyle(4), LineColor(kRed),     ProjectionRange("lowSideband"), Project(RooArgSet(*x,*z)));
      if (GetVar(*TotalPDF,"nBkgPeak")    != NULL) (*TotalPDF)->plotOn(myFrameLowSideBY, Components(*BkgMassAnglePeak), LineStyle(3), LineColor(kViolet),  ProjectionRange("lowSideband"), Project(RooArgSet(*x,*z)));


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

      TLegend* legLowSideBY = new TLegend(0.75, 0.65, 0.97, 0.88, "");
      for (unsigned int i = 0; i < nElements; i++)
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

      (*TotalPDF)->plotOn(myFrameSignalRegionY,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),ProjectionRange("signalRegion"),Project(RooArgSet(*x,*z)));

      if (GetVar(*TotalPDF,"nSig")        != NULL) (*TotalPDF)->plotOn(myFrameSignalRegionY, Components(*Signal),           LineStyle(7), LineColor(kBlue),    ProjectionRange("signalRegion"), Project(RooArgSet(*x,*z)));
      if (GetVar(*TotalPDF,"nMisTagFrac") != NULL) (*TotalPDF)->plotOn(myFrameSignalRegionY, Components(*MassAngleMisTag),  LineStyle(8), LineColor(kAzure+6), ProjectionRange("signalRegion"), Project(RooArgSet(*x,*z)));
      if (GetVar(*TotalPDF,"nBkgComb")    != NULL) (*TotalPDF)->plotOn(myFrameSignalRegionY, Components(*BkgMassAngleComb), LineStyle(4), LineColor(kRed),     ProjectionRange("signalRegion"), Project(RooArgSet(*x,*z)));
      if (GetVar(*TotalPDF,"nBkgPeak")    != NULL) (*TotalPDF)->plotOn(myFrameSignalRegionY, Components(*BkgMassAnglePeak), LineStyle(3), LineColor(kViolet),  ProjectionRange("signalRegion"), Project(RooArgSet(*x,*z)));


      TPaveText* paveTextSignalRegionY = new TPaveText(0.11,0.75,0.4,0.88,"NDC");
      paveTextSignalRegionY->AddText(Form("%s%.1f%s","Signal region: #pm",atof(Utility->GetGenericParam("NSigmaB0").c_str())," < #sigma >"));
      paveTextSignalRegionY->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameSignalRegionY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
      paveTextSignalRegionY->SetTextAlign(11);
      paveTextSignalRegionY->SetBorderSize(0.0);
      paveTextSignalRegionY->SetFillStyle(0);
      paveTextSignalRegionY->SetTextSize(0.04);
      paveTextSignalRegionY->Paint();
      myFrameSignalRegionY->addObject(paveTextSignalRegionY);
      DrawString(LUMI,myFrameSignalRegionY);
      myFrameSignalRegionY->Draw();

      TLegend* legSignalRegionY = new TLegend(0.75, 0.65, 0.97, 0.88, "");
      for (unsigned int i = 0; i < nElements; i++)
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

      (*TotalPDF)->plotOn(myFrameHighSideBY,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),ProjectionRange("highSideband"),Project(RooArgSet(*x,*z)));

      if (GetVar(*TotalPDF,"nSig")        != NULL) (*TotalPDF)->plotOn(myFrameHighSideBY, Components(*Signal),           LineStyle(7), LineColor(kBlue),    ProjectionRange("highSideband"), Project(RooArgSet(*x,*z)));
      if (GetVar(*TotalPDF,"nMisTagFrac") != NULL) (*TotalPDF)->plotOn(myFrameHighSideBY, Components(*MassAngleMisTag),  LineStyle(8), LineColor(kAzure+6), ProjectionRange("highSideband"), Project(RooArgSet(*x,*z)));
      if (GetVar(*TotalPDF,"nBkgComb")    != NULL) (*TotalPDF)->plotOn(myFrameHighSideBY, Components(*BkgMassAngleComb), LineStyle(4), LineColor(kRed),     ProjectionRange("highSideband"), Project(RooArgSet(*x,*z)));
      if (GetVar(*TotalPDF,"nBkgPeak")    != NULL) (*TotalPDF)->plotOn(myFrameHighSideBY, Components(*BkgMassAnglePeak), LineStyle(3), LineColor(kViolet),  ProjectionRange("highSideband"), Project(RooArgSet(*x,*z)));


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

      TLegend* legHighSideBY = new TLegend(0.75, 0.65, 0.97, 0.88, "");
      for (unsigned int i = 0; i < nElements; i++)
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

      (*TotalPDF)->plotOn(myFrameLowSideBZ,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),ProjectionRange("lowSideband"),Project(RooArgSet(*x,*y)));

      if (GetVar(*TotalPDF,"nSig")        != NULL) (*TotalPDF)->plotOn(myFrameLowSideBZ, Components(*Signal),           LineStyle(7), LineColor(kBlue),    ProjectionRange("lowSideband"), Project(RooArgSet(*x,*y)));
      if (GetVar(*TotalPDF,"nMisTagFrac") != NULL) (*TotalPDF)->plotOn(myFrameLowSideBZ, Components(*MassAngleMisTag),  LineStyle(8), LineColor(kAzure+6), ProjectionRange("lowSideband"), Project(RooArgSet(*x,*y)));
      if (GetVar(*TotalPDF,"nBkgComb")    != NULL) (*TotalPDF)->plotOn(myFrameLowSideBZ, Components(*BkgMassAngleComb), LineStyle(4), LineColor(kRed),     ProjectionRange("lowSideband"), Project(RooArgSet(*x,*y)));
      if (GetVar(*TotalPDF,"nBkgPeak")    != NULL) (*TotalPDF)->plotOn(myFrameLowSideBZ, Components(*BkgMassAnglePeak), LineStyle(3), LineColor(kViolet),  ProjectionRange("lowSideband"), Project(RooArgSet(*x,*y)));


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

      TLegend* legLowSideBZ = new TLegend(0.75, 0.65, 0.97, 0.88, "");
      for (unsigned int i = 0; i < nElements; i++)
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

      (*TotalPDF)->plotOn(myFrameSignalRegionZ,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),ProjectionRange("signalRegion"),Project(RooArgSet(*x,*y)));

      if (GetVar(*TotalPDF,"nSig")        != NULL) (*TotalPDF)->plotOn(myFrameSignalRegionZ, Components(*Signal),           LineStyle(7), LineColor(kBlue),    ProjectionRange("signalRegion"), Project(RooArgSet(*x,*y)));
      if (GetVar(*TotalPDF,"nMisTagFrac") != NULL) (*TotalPDF)->plotOn(myFrameSignalRegionZ, Components(*MassAngleMisTag),  LineStyle(8), LineColor(kAzure+6), ProjectionRange("signalRegion"), Project(RooArgSet(*x,*y)));
      if (GetVar(*TotalPDF,"nBkgComb")    != NULL) (*TotalPDF)->plotOn(myFrameSignalRegionZ, Components(*BkgMassAngleComb), LineStyle(4), LineColor(kRed),     ProjectionRange("signalRegion"), Project(RooArgSet(*x,*y)));
      if (GetVar(*TotalPDF,"nBkgPeak")    != NULL) (*TotalPDF)->plotOn(myFrameSignalRegionZ, Components(*BkgMassAnglePeak), LineStyle(3), LineColor(kViolet),  ProjectionRange("signalRegion"), Project(RooArgSet(*x,*y)));


      TPaveText* paveTextSignalRegionZ = new TPaveText(0.11,0.75,0.4,0.88,"NDC");
      paveTextSignalRegionZ->AddText(Form("%s%.1f%s","Signal region: #pm",atof(Utility->GetGenericParam("NSigmaB0").c_str())," < #sigma >"));
      paveTextSignalRegionZ->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameSignalRegionZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
      paveTextSignalRegionZ->SetTextAlign(11);
      paveTextSignalRegionZ->SetBorderSize(0.0);
      paveTextSignalRegionZ->SetFillStyle(0);
      paveTextSignalRegionZ->SetTextSize(0.04);
      paveTextSignalRegionZ->Paint();
      myFrameSignalRegionZ->addObject(paveTextSignalRegionZ);
      DrawString(LUMI,myFrameSignalRegionZ);
      myFrameSignalRegionZ->Draw();

      TLegend* legSignalRegionZ = new TLegend(0.75, 0.65, 0.97, 0.88, "");
      for (unsigned int i = 0; i < nElements; i++)
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

      (*TotalPDF)->plotOn(myFrameHighSideBZ,Name((*TotalPDF)->getPlotLabel()),LineColor(kBlack),ProjectionRange("highSideband"), Project(RooArgSet(*x,*y)));

      if (GetVar(*TotalPDF,"nSig")        != NULL) (*TotalPDF)->plotOn(myFrameHighSideBZ, Components(*Signal),           LineStyle(7), LineColor(kBlue),    ProjectionRange("highSideband"), Project(RooArgSet(*x,*y)));
      if (GetVar(*TotalPDF,"nMisTagFrac") != NULL) (*TotalPDF)->plotOn(myFrameHighSideBZ, Components(*MassAngleMisTag),  LineStyle(8), LineColor(kAzure+6), ProjectionRange("highSideband"), Project(RooArgSet(*x,*y)));
      if (GetVar(*TotalPDF,"nBkgComb")    != NULL) (*TotalPDF)->plotOn(myFrameHighSideBZ, Components(*BkgMassAngleComb), LineStyle(4), LineColor(kRed),     ProjectionRange("highSideband"), Project(RooArgSet(*x,*y)));
      if (GetVar(*TotalPDF,"nBkgPeak")    != NULL) (*TotalPDF)->plotOn(myFrameHighSideBZ, Components(*BkgMassAnglePeak), LineStyle(3), LineColor(kViolet),  ProjectionRange("highSideband"), Project(RooArgSet(*x,*y)));


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

      TLegend* legHighSideBZ = new TLegend(0.75, 0.65, 0.97, 0.88, "");
      for (unsigned int i = 0; i < nElements; i++)
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


  // ##############
  // # Save plots #
  // ##############
  Canv->Modified();
  Canv->Update();
  if (SAVEPLOT == true)
    {
      myString.clear(); myString.str("");
      myString << (*TotalPDF)->getPlotLabel() << "_" << ID << ".pdf";
      Canv->Print(myString.str().c_str());

      for (unsigned int i = 0; i < nCanv; i++)
	{
	  myString.clear(); myString.str("");
	  myString << (*TotalPDF)->getPlotLabel() << "_localCanv" << i << "_" << ID << ".pdf";
	  localCanv[i]->Print(myString.str().c_str());
	}
    }

  if (SETBATCH == true) for (unsigned int i = 0; i < nCanv; i++) delete localCanv[i];
  else
    for (unsigned int i = 0; i < nCanv; i++)
      {
	localCanv[i]->Modified();
	localCanv[i]->Update();
      }


  VarsYield.Clear();
  return fitResult;
}


void IterativeMass2AnglesFitq2Bins (RooDataSet* dataSet,
				    bool useEffPDF,
				    double PsiYieldGoodTag, double PsiYieldGoodTagErr,
				    double PsiYieldMisTag,
				    RooRealVar* x, RooRealVar* y, RooRealVar* z,
				    int specBin,
				    unsigned int FitType,
				    vector<TH1D*>* VecHistoMeas,
				    vector<double>* q2Bins,
				    vector<vector<unsigned int>*>* configParam, vector<vector<string>*>* fitParam,
				    pair< vector<TF2*>*,vector<TF2*>* > effFuncs,
				    RooArgSet* vecConstr,
				    unsigned int ID)
{
  // ###################
  // # Local variables #
  // ###################
  vector<string>* vecParStr = NULL;
  stringstream myString;

  double effPsiGoodTag     = 1.0;
  double effPsiGoodTagErr  = 0.0;

  double effPsiMisTag      = 1.0;
  double effPsiMisTagErr   = 0.0;

  double effMuMuGoodTag    = 1.0;
  double effMuMuGoodTagErr = 0.0;

  double effMuMuMisTag     = 1.0;
  double effMuMuMisTagErr  = 0.0;

  double value1 = 0.0;
  double value2 = 0.0;

  TCanvas*    cq2Bins[q2Bins->size()-1];
  RooDataSet* dataSet_q2Bins[q2Bins->size()-1];
  RooAbsPdf*  TotalPDFq2Bins[q2Bins->size()-1];
  TPaveText*  extText[q2Bins->size()-1];

  RooFitResult* fitResult;
  RooAbsReal*   EffPDFintegral = NULL;
  RooAbsPdf*    EffPDFgoodTag  = NULL;
  RooAbsPdf*    EffPDFmisTag   = NULL;

  double NLLvalue;
  

  // ####################################
  // # Read I[S*E] for resonant channel #
  // ####################################
  if (useEffPDF == true)
    {
      // ######################
      // # Good-tagged events #
      // ######################
      Utility->ReadParVsq2Bins(PARAMETERFILEIN,"I[S*E]okTag",&vecParStr);
      myString.clear(); myString.str("");
      myString << vecParStr->operator[](atoi(Utility->GetGenericParam("NormJPSInotPSIP").c_str()) == true ? Utility->GetJPsiBin(q2Bins) : Utility->GetPsiPBin(q2Bins)).c_str();
      SetValueAndErrors(NULL,"",1.0,&myString,&effPsiGoodTag,&effPsiGoodTagErr,&effPsiGoodTagErr);
      cout << "\n[ExtractYield::IterativeMass2AnglesFitq2Bins]\t@@@ Integral of S*E over angular variables for normalization channel good-tagged events: " << effPsiGoodTag << " +/- " << effPsiGoodTagErr << " @@@" << endl;


      // #####################
      // # Mis-tagged events #
      // #####################
      Utility->ReadParVsq2Bins(PARAMETERFILEIN,"I[S*E]misTag",&vecParStr);
      myString.clear(); myString.str("");
      myString << vecParStr->operator[](atoi(Utility->GetGenericParam("NormJPSInotPSIP").c_str()) == true ? Utility->GetJPsiBin(q2Bins) : Utility->GetPsiPBin(q2Bins)).c_str();
      SetValueAndErrors(NULL,"",1.0,&myString,&effPsiMisTag,&effPsiMisTagErr,&effPsiMisTagErr);
      cout << "\n[ExtractYield::IterativeMass2AnglesFitq2Bins]\t@@@ Integral of S*E over angular variables for normalization channel mis-tagged events: " << effPsiMisTag << " +/- " << effPsiMisTagErr << " @@@" << endl;
    }
  

  fileFitResults << "\n@@@@@@@@@@@@@@@@@@@@@ Angular analysis vs q^2 @@@@@@@@@@@@@@@@@@@@@@" << endl;
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
	  vecParStr->clear();
	  delete vecParStr;
	  
	  continue;
	}


      // ######################
      // # Make external text #
      // ######################
      extText[i] = new TPaveText(0.65,0.5,0.97,0.63,"NDC");
      extText[i]->AddText(Form("%s%.2f%s%.2f%s","q#lower[0.4]{^{2}}: ",q2Bins->operator[](i)," #font[122]{\55} ",q2Bins->operator[](i+1)," GeV#lower[0.4]{^{2}}"));
      extText[i]->SetTextAlign(11);
      extText[i]->SetBorderSize(0.0);
      extText[i]->SetFillStyle(0);
      extText[i]->SetTextSize(0.04);


      fileFitResults << "\nBin[" << i << "]: " << q2Bins->operator[](i) << " <= q^2 < " << q2Bins->operator[](i+1) << endl;
      
      myString.clear(); myString.str("");
      myString << "c_" << i;
      cq2Bins[i] = new TCanvas(myString.str().c_str(), myString.str().c_str(), 20, 20, 1800, 1800);
      cq2Bins[i]->Divide(3,3);

      myString.clear(); myString.str("");
      myString << "(mumuMass*mumuMass) > " << q2Bins->operator[](i) << " && (mumuMass*mumuMass) <= " << q2Bins->operator[](i+1);
      cout << "\n[ExtractYield::IterativeMass2AnglesFitq2Bins]\tCut string: " << myString.str() << endl;
      dataSet_q2Bins[i] = (RooDataSet*)dataSet->reduce(myString.str().c_str());
      cout << "[ExtractYield::IterativeMass2AnglesFitq2Bins]\tNumber of events: " << dataSet_q2Bins[i]->sumEntries() << endl;

      unsigned int countMisTag  = 0;
      unsigned int countGoodTag = 0;
      for (int j = 0; j < static_cast<int>(dataSet_q2Bins[i]->sumEntries()); j++)
	{
	  if (dataSet_q2Bins[i]->get(j)->getRealValue("truthMatchSignal") == true)
	    {
	      if (dataSet_q2Bins[i]->get(j)->getRealValue("rightFlavorTag") == 0.0) countMisTag++;
	      else                                                                  countGoodTag++;
	    }
	}
      cout << "[ExtractYield::IterativeMass2AnglesFitq2Bins]\tDynamic mis-tag fraction: " << static_cast<double>(countMisTag) / static_cast<double>(countMisTag + countGoodTag) << " = (" << countMisTag << "/(" << countMisTag << "+" << countGoodTag << "))" << endl;

      myString.clear(); myString.str("");
      myString << "TotalPDFq2Bin_" << i;
      InstantiateMass2AnglesFit(&TotalPDFq2Bins[i],useEffPDF,x,y,z,myString.str(),FitType,configParam,fitParam,q2Bins,i,i,make_pair(effFuncs.first->operator[](i),effFuncs.second->operator[](i)));


      // #####################
      // # Initialize p.d.f. #
      // #####################
      CopyFitResults(TotalPDFq2Bins[i],i,fitParam,countMisTag,countGoodTag);


      // #####################
      // # Apply constraints #
      // #####################
      ClearVars(vecConstr);
      if ((FitType != 36) && (FitType != 56) && (FitType != 76))
      	{
      	  BuildMassConstraints(vecConstr,TotalPDFq2Bins[i],"sign");
      	  BuildMassConstraints(vecConstr,TotalPDFq2Bins[i],"peak");
      	  if ((strcmp(CTRLfitWRKflow.c_str(),"trueAll&NoFFrac") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"trueAll&FitFrac") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"allEvts") == 0)) BuildMassConstraints(vecConstr,TotalPDFq2Bins[i],"mistag");
      	  if ((GetVar(TotalPDFq2Bins[i],"nSig") != NULL) || (GetVar(TotalPDFq2Bins[i],"nMisTagFrac") != NULL)) BuildAngularConstraints(vecConstr,TotalPDFq2Bins[i],"peak");
      	}
      if (((FitType != 46) && (FitType != 56) &&
      	   (FitType != 66) && (FitType != 76)) ||
      	  ((FitType == 26) && (specBin != Utility->GetJPsiBin(q2Bins)) && (specBin != Utility->GetPsiPBin(q2Bins))))
      	{
      	  BuildPhysicsConstraints(vecConstr,TotalPDFq2Bins[i],"FsS");
      	  BuildPhysicsConstraints(vecConstr,TotalPDFq2Bins[i],"AsS");
      	}
      // @TMP@
      // BuildPhysicsConstraints(vecConstr,TotalPDFq2Bins[i],"AfbBound"); // # Special constraint for AFB boundaries #
      // BuildPhysicsConstraints(vecConstr,TotalPDFq2Bins[i],"AsBound");  // # Special constraint for AS  boundaries #


      // ###################################
      // # Print variables and constraints #
      // ###################################
      PrintVariables(TotalPDFq2Bins[i]->getVariables(),"vars");
      PrintVariables(vecConstr,"cons");


      // ###################
      // # Perform the fit #
      // ###################
      fitResult = MakeMass2AnglesFit(dataSet_q2Bins[i],&TotalPDFq2Bins[i],x,y,z,cq2Bins[i],FitType,vecConstr,&NLLvalue,extText[i],ID);
      int nTrials = 0;
      while ((MAXTRIALS > 0) && (nTrials < MAXTRIALS) && (CheckGoodFit(fitResult) == false))
	{
	  delete fitResult;
	  GenerateFitParameters(TotalPDFq2Bins[i],fitParam,ID,q2Bins,i,"FlAfb");

	  cout << "\n[ExtractYield::IterativeMass2AnglesFitq2Bins]\t@@@ Failed trial #" << nTrials+1 << " @@@" << endl;
	  cout << "[ExtractYield::IterativeMass2AnglesFitq2Bins]\tNew trial with parameters: Fl = " << GetVar(TotalPDFq2Bins[i],"FlS")->getVal() << "; Afb = " << GetVar(TotalPDFq2Bins[i],"AfbS")->getVal() << "\n" << endl;

	  fitResult = MakeMass2AnglesFit(dataSet_q2Bins[i],&TotalPDFq2Bins[i],x,y,z,cq2Bins[i],FitType,vecConstr,&NLLvalue,extText[i],ID);
	  nTrials++;
	}
      if (CheckGoodFit(fitResult) == true) cout << "\n[ExtractYield::IterativeMass2AnglesFitq2Bins]\t@@@ Fit converged ! @@@" << endl;
      else                                 cout << "\n[ExtractYield::IterativeMass2AnglesFitq2Bins]\t@@@ Fit didn't converge ! @@@" << endl;


      // ##############################################
      // # Save fit results back into prarameter file #
      // ##############################################
      if (GetVar(TotalPDFq2Bins[i],"FlS")  != NULL) value1 = GetVar(TotalPDFq2Bins[i],"FlS")->getVal();
      if (GetVar(TotalPDFq2Bins[i],"AfbS") != NULL) value2 = GetVar(TotalPDFq2Bins[i],"AfbS")->getVal();

      vecParStr = SaveFitResults(TotalPDFq2Bins[i],i,fitParam,configParam,vecConstr);
      Utility->SaveFitValues(PARAMETERFILEOUT,vecParStr,i);
      vecParStr->clear();
      delete vecParStr;
      vecParStr = NULL;


      // #############################
      // # Compute I[S*E] for signal #
      // #############################
      if (useEffPDF == true)
      	{
      	  // ######################
      	  // # Good-tagged events #
      	  // ######################
	  EffPDFgoodTag = (RooAbsPdf*)AngleS->clone("EffPDFgoodTag");
      	  PrintVariables(EffPDFgoodTag->getVariables(),"vars");

      	  EffPDFintegral    = EffPDFgoodTag->createIntegral(RooArgSet(*y,*z));
      	  effMuMuGoodTag    = EffPDFintegral->getVal();
      	  effMuMuGoodTagErr = EffPDFintegral->getPropagatedError(*fitResult);

      	  cout << "\n[ExtractYield::IterativeMass2AnglesFitq2Bins]\t===> Integral of S*E over angular variables for signal good-tagged events <===" << endl;
      	  cout << effMuMuGoodTag << "   -" << effMuMuGoodTagErr << "   " << effMuMuGoodTagErr << endl;


      	  // #####################
      	  // # Mis-tagged events #
      	  // #####################
      	  if (GetVar(TotalPDFq2Bins[i],"nMisTagFrac") != NULL)
      	    {
	      EffPDFmisTag = (RooAbsPdf*)AngleMisTag->clone("EffPDFmisTag");
      	      PrintVariables(EffPDFmisTag->getVariables(),"vars");

      	      EffPDFintegral   = EffPDFmisTag->createIntegral(RooArgSet(*y,*z));
      	      effMuMuMisTag    = EffPDFintegral->getVal();
      	      effMuMuMisTagErr = EffPDFintegral->getPropagatedError(*fitResult);

      	      cout << "\n[ExtractYield::IterativeMass2AnglesFitq2Bins]\t===> Integral of S*E over angular variables for signal mis-tagged events <===" << endl;
      	      cout << effMuMuMisTag << "   -" << effMuMuMisTagErr << "   " << effMuMuMisTagErr << endl;
      	    }
      	}


      // ########################################
      // # Save outcome of the fit in histogram #
      // ########################################
      if (GetVar(TotalPDFq2Bins[i],"nSig") != NULL)
	{
	  VecHistoMeas->operator[](0)->SetBinContent(i+1,value1);
	  VecHistoMeas->operator[](0)->SetBinError(i+1,GetVar(TotalPDFq2Bins[i],"FlS")->getError());

	  VecHistoMeas->operator[](1)->SetBinContent(i+1,value2);
	  VecHistoMeas->operator[](1)->SetBinError(i+1,GetVar(TotalPDFq2Bins[i],"AfbS")->getError());

	  if ((FitType != 36) && (FitType != 56) && (FitType != 76))
	    {
	      double nEvGoodTag    = GetVar(TotalPDFq2Bins[i],"nSig")->getVal();
	      double nEvGoodTagErr = GetVar(TotalPDFq2Bins[i],"nSig")->getError();

	      double nEvMisTag = 0.0;
	      if (GetVar(TotalPDFq2Bins[i],"nMisTagFrac") != NULL)
		nEvMisTag = GetVar(TotalPDFq2Bins[i],"nSig")->getVal() / (1. - GetVar(TotalPDFq2Bins[i],"nMisTagFrac")->getVal()) * GetVar(TotalPDFq2Bins[i],"nMisTagFrac")->getVal();

	      double num       = (nEvGoodTag / effMuMuGoodTag) + (nEvMisTag / effMuMuMisTag);
	      double den       = (PsiYieldGoodTag / effPsiGoodTag) + (PsiYieldMisTag / effPsiMisTag);
	      double dBFdq2    = num / den * (atoi(Utility->GetGenericParam("NormJPSInotPSIP").c_str()) == true ? Utility->JPsiBF : Utility->PsiPBF) / (q2Bins->operator[](i+1) - q2Bins->operator[](i)) / 1e-7;
	      double dBFdq2Err = dBFdq2 * sqrt( (pow(num * nEvGoodTagErr      / nEvGoodTag,2.)      + pow(nEvGoodTag         / (effMuMuGoodTag*effMuMuGoodTag) * effMuMuGoodTagErr,2.) + pow(nEvMisTag          / (effMuMuMisTag*effMuMuMisTag) * effMuMuMisTagErr,2.)) / (num*num) +
					        (pow(den * PsiYieldGoodTagErr / PsiYieldGoodTag,2.) + pow(PsiYieldGoodTagErr / (effPsiGoodTag*effPsiGoodTag)   * effPsiGoodTagErr,2.)  + pow(PsiYieldGoodTagErr / (effPsiMisTag*effPsiMisTag)   * effPsiMisTagErr,2.))  / (den*den) );

	      VecHistoMeas->operator[](2)->SetBinContent(i+1,dBFdq2);
	      VecHistoMeas->operator[](2)->SetBinError(i+1,dBFdq2Err);
	      
	      fileFitResults << "@@@@@@ dBF/dq^2 @@@@@@" << endl;
	      fileFitResults << "dBF/dq^2: " << VecHistoMeas->operator[](2)->GetBinContent(i+1) << " -/+ " << VecHistoMeas->operator[](2)->GetBinError(i+1) << endl;

	      fileFitResults << "===> dBF/dq^2 to cut and paste in config. file <===" << endl;
	      fileFitResults << VecHistoMeas->operator[](2)->GetBinContent(i+1) << "   -" << dBFdq2Err << "   " << dBFdq2Err << endl;

	      fileFitResults << "===> I[S*E] good-tagged to cut and paste in config. file <===" << endl;
	      fileFitResults << effMuMuGoodTag << "   -" << effMuMuGoodTagErr << "   " << effMuMuGoodTagErr << endl;

	      fileFitResults << "===> I[S*E] mis-tagged to cut and paste in config. file <===" << endl;
	      fileFitResults << effMuMuMisTag << "   -" << effMuMuMisTagErr << "   " << effMuMuMisTagErr << endl;

	      fileFitResults << "====================================================================" << endl;
	    }
	}


      // #############################################
      // # Save observables in systematic error file #
      // #############################################
      myString.clear(); myString.str("");
      if (CheckGoodFit(fitResult) == true) myString << ID << "   " << value1 << "   " << value2 << "   " << VecHistoMeas->operator[](2)->GetBinContent(i+1) << "   " << effMuMuGoodTag << "   " << effMuMuMisTag << "   " << NLLvalue;
      else                                 myString << ID << "   " << -2.0 << "   " << -2.0 << "   " << -2.0 << "   " << -2.0 << "   " << -2.0 << "   " << -2.0;
      fileFitSystematics << myString.str() << endl;
    }
}


void MakeMass2AnglesToy (RooAbsPdf* TotalPDF, RooRealVar* x, RooRealVar* y, RooRealVar* z, TCanvas* Canv, unsigned int FitType, unsigned int nToy, int specBin, vector<vector<string>*>* fitParam, RooArgSet* vecConstr, string fileName, vector<double>* q2Bins)
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


  // ################################
  // # Fix comb. angular background #
  // ################################
  cout << "[ExtractYield::MakeMass2AnglesToy]\t@@@ Fixing comb. angular background after sideband fit @@@" << endl;
  for (unsigned int i = 0 ; i < NCOEFFPOLYBKG; i++)
    {
      myString.clear(); myString.str("");
      myString << "c1Poly" << i;
      if (GetVar(TotalPDF,myString.str().c_str()) != NULL) GetVar(TotalPDF,myString.str().c_str())->setConstant(true);

      myString.clear(); myString.str("");
      myString << "c2Poly" << i;
      if (GetVar(TotalPDF,myString.str().c_str()) != NULL) GetVar(TotalPDF,myString.str().c_str())->setConstant(true);
    }


  // #####################
  // # Apply constraints #
  // #####################
  ClearVars(vecConstr);
  BuildMassConstraints(vecConstr,TotalPDF,"sign");
  BuildMassConstraints(vecConstr,TotalPDF,"peak");
  if ((strcmp(CTRLfitWRKflow.c_str(),"trueAll&NoFFrac") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"trueAll&FitFrac") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"allEvts") == 0)) BuildMassConstraints(vecConstr,TotalPDF,"mistag");
  BuildAngularConstraints(vecConstr,TotalPDF,"peak");
  if ((specBin != Utility->GetJPsiBin(q2Bins)) && (specBin != Utility->GetPsiPBin(q2Bins)))
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
  RooMCStudy* MyToy = new RooMCStudy(*TotalPDF,RooArgSet(*x,*y,*z),Extended(true),ExternalConstraints(*constrToys),FitOptions(Extended(true),ExternalConstraints(*vecConstr),Minos(USEMINOS)));
  MyToy->generateAndFit(nToy,nEntryToy,true);
  Canv->Divide(6,5);


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

  if ((GetVar(TotalPDF,"var1") != NULL) && (GetVar(TotalPDF,"var1")->getError() != 0.0) && (IsInConstraints(vecConstr,"var1") == false))
    {
      tmpVar = GetVar(TotalPDF,"var1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("var1"))->operator[](specBin).c_str()) - 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("var1"))->operator[](specBin).c_str())),
 			      atof(fitParam->operator[](Utility->GetFitParamIndx("var1"))->operator[](specBin).c_str()) + 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("var1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"var2") != NULL) && (GetVar(TotalPDF,"var2")->getError() != 0.0) && (IsInConstraints(vecConstr,"var2") == false))
    {
      tmpVar = GetVar(TotalPDF,"var2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("var2"))->operator[](specBin).c_str()) - 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("var2"))->operator[](specBin).c_str())),
 			      atof(fitParam->operator[](Utility->GetFitParamIndx("var2"))->operator[](specBin).c_str()) + 0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("var2"))->operator[](specBin).c_str())));
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
  
  if ((GetVar(TotalPDF,"sigmaMisTag1") != NULL) && (GetVar(TotalPDF,"sigmaMisTag1")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaMisTag1") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaMisTag1");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag1"))->operator[](specBin).c_str()) -
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag1"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag1"))->operator[](specBin).c_str()) +
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag1"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"sigmaMisTag2") != NULL) && (GetVar(TotalPDF,"sigmaMisTag2")->getError() != 0.0) && (IsInConstraints(vecConstr,"sigmaMisTag2") == false))
    {
      tmpVar = GetVar(TotalPDF,"sigmaMisTag2");
      myFrame = tmpVar->frame(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag2"))->operator[](specBin).c_str()) -
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag2"))->operator[](specBin).c_str())),
  			      atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag2"))->operator[](specBin).c_str()) +
  			      0.8*fabs(atof(fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag2"))->operator[](specBin).c_str())));
      MyToy->plotParamOn(myFrame);
      Canv->cd(it++);
      myFrame->Draw();
      myFrame = MyToy->plotPull(*tmpVar,-5,5,30,true);
      Canv->cd(it++);
      myFrame->Draw();
    }

  if ((GetVar(TotalPDF,"fracMisTag") != NULL) && (GetVar(TotalPDF,"fracMisTag")->getError() != 0.0) && (IsInConstraints(vecConstr,"fracMisTag") == false))
    {
      tmpVar = GetVar(TotalPDF,"fracMisTag");
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

  if ((GetVar(TotalPDF,"nMisTagFrac") != NULL) && (GetVar(TotalPDF,"nMisTagFrac")->getError() != 0.0) && (IsInConstraints(vecConstr,"nMisTagFrac") == false))
    {
      tmpVar = GetVar(TotalPDF,"nMisTagFrac");
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
  Canv->Modified();
  Canv->Update();


  // #############################################
  // # Re-make all the fits and save the results #
  // #############################################
  TH1D* histoDiff1 = new TH1D("histoDiff1","histoDiff1",100,-1.0,1.0);
  histoDiff1->SetFillColor(kAzure+6);
  histoDiff1->SetXTitle("(fit #font[122]{\55} pdf)");
  histoDiff1->SetYTitle("Entries [#]");

  TH1D* histoPull1 = new TH1D("histoPull1","histoPull1",100,-5.0,5.0);
  histoPull1->SetFillColor(kAzure+6);
  histoPull1->SetXTitle("(fit #font[122]{\55} pdf) / #sigma");
  histoPull1->SetYTitle("Entries [#]");

  TH1D* histoChiSq1 = new TH1D("histoChiSq1","histoChiSq1",100,0.0,30.0);
  histoChiSq1->SetFillColor(kAzure+6);
  histoChiSq1->SetXTitle("(fit #font[122]{\55} pdf){^2} / #sigma{^2}");
  histoChiSq1->SetYTitle("Entries [#]");

  TH1D* histoVal1 = new TH1D("histoVal1","histoVal1",100,0.0,1.0);
  histoVal1->SetFillColor(kAzure+6);
  histoVal1->SetXTitle("fit");
  histoVal1->SetYTitle("Entries [#]");

  TH1D* histoNLL1 = new TH1D("histoNLL1","histoNLL1",100,1.0,-1.0);
  histoNLL1->SetFillColor(kGreen-7);
  histoNLL1->SetXTitle("NLL");
  histoNLL1->SetYTitle("Entries [#]");

  TH1D* histoDiff2 = new TH1D("histoDiff2","histoDiff2",100,-1.0,1.0);
  histoDiff2->SetFillColor(kAzure+6);
  histoDiff2->SetXTitle("(fit #font[122]{\55} pdf)");
  histoDiff2->SetYTitle("Entries [#]");

  TH1D* histoPull2 = new TH1D("histoPull2","histoPull2",100,-5.0,5.0);
  histoPull2->SetFillColor(kAzure+6);
  histoPull2->SetXTitle("(fit #font[122]{\55} pdf) / #sigma");
  histoPull2->SetYTitle("Entries [#]");

  TH1D* histoChiSq2 = new TH1D("histoChiSq2","histoChiSq2",100,0.0,30.0);
  histoChiSq2->SetFillColor(kAzure+6);
  histoChiSq2->SetXTitle("(fit #font[122]{\55} pdf){^2} / #sigma{^2}");
  histoChiSq2->SetYTitle("Entries [#]");

  TH1D* histoVal2 = new TH1D("histoVal2","histoVal2",100,-1.0,1.0);
  histoVal2->SetFillColor(kAzure+6);
  histoVal2->SetXTitle("fit");
  histoVal2->SetYTitle("Entries [#]");

  TH1D* histoNLL2 = new TH1D("histoNLL2","histoNLL2",100,1.0,-1.0);
  histoNLL2->SetFillColor(kGreen-7);
  histoNLL2->SetXTitle("NLL");
  histoNLL2->SetYTitle("Entries [#]");

  TCanvas* cB0Toy = new TCanvas("cB0Toy","cB0Toy", 20, 20, 1800, 1800);
  cB0Toy->Divide(3,3);

  cout << "\n[ExtractYield::MakeMass2AnglesToy]\t@@@ Now fit total TOY invariant mass and angles @@@" << endl;
  RooDataSet* toySample;
  RooFitResult* fitResult;
  double NLLvalue;

  string varName;
  for (unsigned int i = 0; i < nToy; i++)
    {
      cout << "\n[ExtractYield::MakeMass2AnglesToy]\t@@@ Now fitting toy #" << i << " @@@\n" << endl;

      toySample = (RooDataSet*)MyToy->genData(i);
      CopyFitResults(TotalPDF,specBin,fitParam);
      fitResult = MakeMass2AnglesFit(toySample,&TotalPDF,x,y,z,cB0Toy,FitType,vecConstr,&NLLvalue,NULL,i);
      int nTrials = 0;
      while ((MAXTRIALS > 0) && (nTrials < MAXTRIALS) && (CheckGoodFit(fitResult) == false))
	{
	  delete fitResult;
	  GenerateFitParameters(TotalPDF,fitParam,i,q2Bins,specBin,"FlAfb");

	  cout << "\n[ExtractYield::MakeMass2AnglesToy]\t@@@ Failed trial #" << nTrials+1 << " @@@" << endl;
	  cout << "[ExtractYield::MakeMass2AnglesToy]\tNew trial with parameters: Fl = " << GetVar(TotalPDF,"FlS")->getVal() << "; Afb = " << GetVar(TotalPDF,"AfbS")->getVal() << "\n" << endl;

	  fitResult = MakeMass2AnglesFit(toySample,&TotalPDF,x,y,z,cB0Toy,FitType,vecConstr,&NLLvalue,NULL,i);
	  nTrials++;
	}
      if (CheckGoodFit(fitResult) == true) cout << "\n[ExtractYield::MakeMass2AnglesToy]\t@@@ Fit converged ! @@@" << endl;
      else                                 cout << "\n[ExtractYield::MakeMass2AnglesToy]\t@@@ Fit didn't converge ! @@@" << endl;


      // ######################################################
      // # To verify the outcome of the RooFit toy-MC studies #
      // ######################################################
      if ((CheckGoodFit(fitResult) == true) && (GetVar(TotalPDF,"FlS") != NULL) && (GetVar(TotalPDF,"AfbS") != NULL))
	{
	  varName = "FlS";
	  myString.clear(); myString.str("");
	  myString << fitParam->operator[](Utility->GetFitParamIndx(varName.c_str()))->operator[](specBin);
	  if (TotalPDF->getVariables()->getRealValue(varName.c_str()) > atof(myString.str().c_str()))
	    {
	      histoPull1->Fill((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(myString.str().c_str())) / fabs(GetVar(TotalPDF,varName.c_str())->getErrorLo()));
	      histoChiSq1->Fill(pow((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(myString.str().c_str())) / fabs(GetVar(TotalPDF,varName.c_str())->getErrorLo()),2.));
	    }
	  else
	    {
	      histoPull1->Fill((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(myString.str().c_str())) / fabs(GetVar(TotalPDF,varName.c_str())->getErrorHi()));
	      histoChiSq1->Fill(pow((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(myString.str().c_str())) / fabs(GetVar(TotalPDF,varName.c_str())->getErrorHi()),2.));
	    }
	  histoDiff1->Fill(TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(myString.str().c_str()));
	  histoVal1->Fill(TotalPDF->getVariables()->getRealValue(varName.c_str()));
	  histoNLL1->Fill(NLLvalue);
	  
	  varName = "AfbS";
	  myString.clear(); myString.str("");
	  myString << fitParam->operator[](Utility->GetFitParamIndx(varName.c_str()))->operator[](specBin);
	  if (TotalPDF->getVariables()->getRealValue(varName.c_str()) > atof(myString.str().c_str()))
	    {
	      histoPull2->Fill((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(myString.str().c_str())) / fabs(GetVar(TotalPDF,varName.c_str())->getErrorLo()));
	      histoChiSq2->Fill(pow((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(myString.str().c_str())) / fabs(GetVar(TotalPDF,varName.c_str())->getErrorLo()),2.));
	    }
	  else
	    {
	      histoPull2->Fill((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(myString.str().c_str())) / fabs(GetVar(TotalPDF,varName.c_str())->getErrorHi()));
	      histoChiSq2->Fill(pow((TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(myString.str().c_str())) / fabs(GetVar(TotalPDF,varName.c_str())->getErrorHi()),2.));
	    }
	  histoDiff2->Fill(TotalPDF->getVariables()->getRealValue(varName.c_str()) - atof(myString.str().c_str()));
	  histoVal2->Fill(TotalPDF->getVariables()->getRealValue(varName.c_str()));
	  histoNLL2->Fill(NLLvalue);
	}
      delete fitResult;
    }


  TCanvas* cNLL1 = new TCanvas("cNLL1","cNLL1",10,10,900,700);
  cNLL1->Divide(2,3);
  cNLL1->cd(1);
  histoDiff1->Draw();
  cNLL1->cd(2);
  histoNLL1->Draw();
  cNLL1->cd(3);
  histoPull1->Draw();
  cNLL1->cd(4);
  histoVal1->Draw();
  cNLL1->cd(5);
  histoChiSq1->Draw();
  cNLL1->Modified();
  cNLL1->Update();

  TCanvas* cNLL2 = new TCanvas("cNLL2","cNLL2",10,10,900,700);
  cNLL2->Divide(2,3);
  cNLL2->cd(1);
  histoDiff2->Draw();
  cNLL2->cd(2);
  histoNLL2->Draw();
  cNLL2->cd(3);
  histoPull2->Draw();
  cNLL2->cd(4);
  histoVal2->Draw();
  cNLL2->cd(5);
  histoChiSq2->Draw();
  cNLL2->Modified();
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
      myString << "_FL_" << specBin << "_Chi2.root";
      fNLL = new TFile(fileName.replace(fileName.find("_FL_"),15,myString.str()).c_str(),"RECREATE");
      fNLL->cd();
      histoChiSq1->Write();
      fNLL->Close();
      delete fNLL;

      myString.clear(); myString.str("");
      myString << "_FL_" << specBin << "_Val1.root";
      fNLL = new TFile(fileName.replace(fileName.find("_FL_"),15,myString.str()).c_str(),"RECREATE");
      fNLL->cd();
      histoVal1->Write();
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
      myString << "_AFB_" << specBin << "_Chi2.root";
      fNLL = new TFile(fileName.replace(fileName.find("_AFB_"),16,myString.str()).c_str(),"RECREATE");
      fNLL->cd();
      histoChiSq2->Write();
      fNLL->Close();
      delete fNLL;

      myString.clear(); myString.str("");
      myString << "_AFB_" << specBin << "_Val2.root";
      fNLL = new TFile(fileName.replace(fileName.find("_AFB_"),16,myString.str()).c_str(),"RECREATE");
      fNLL->cd();
      histoVal2->Write();
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


int main(int argc, char** argv)
{
  if (argc >= 4)
    {
      // ##################
      // # Main variables #
      // ##################
      stringstream myString;
      
      string fileName           = "";
      string correct4Efficiency = "";
      string tmpFileName        = "";

      int specBin               = -1;
      unsigned int fileIndx     = 0;
      unsigned int nToy         = 0;
      unsigned int FitType      = atoi(argv[1]);

      bool useEffPDF            = false;

      TFile* NtplFile           = NULL;

      
      if (((((FitType >= 1)  && (FitType <= 6))   ||
	    (FitType == 36)                       ||
	    ((FitType >= 41) && (FitType <= 46))  ||
	    ((FitType >= 61) && (FitType <= 66))  ||
	    ((FitType >= 81) && (FitType <= 86))) && (argc >= 4)) ||
	  ((FitType == 96)   &&                      (argc == 6)) ||
	  ((FitType >= 21)   && (FitType <= 26)   && (argc == 8)) ||
	  (((FitType == 56)  || (FitType == 76))  && (argc == 4)))
	{
	  ParameterFILE = PARAMETERFILEIN;


 	  // ###################
	  // # Read parameters #
 	  // ###################
	  Utility = new Utils(false);
	  Utility->ReadAllBins(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins);


	  // #################################
	  // # Check that FitType is correct #
	  // #################################
	  if (((FitType >= 1)  && (FitType <= 6))  ||
	      ((FitType >= 21) && (FitType <= 26)) ||
	      (FitType == 36)                      ||
	      ((FitType >= 41) && (FitType <= 46)) ||
	      (FitType == 56)                      ||
	      ((FitType >= 61) && (FitType <= 66)) ||
	      (FitType == 76)                      ||
	      ((FitType >= 81) && (FitType <= 86)) ||
	      (FitType == 96))
	    {
	      fileName           = argv[2];
	      correct4Efficiency = argv[3];
	      if (((FitType >= 41) && (FitType <= 46)) || (FitType == 56)) specBin = Utility->GetJPsiBin(&q2Bins);
	      if (((FitType >= 61) && (FitType <= 66)) || (FitType == 76)) specBin = Utility->GetPsiPBin(&q2Bins);
	    }


	  // ###################################
	  // # Check that FitOption is correct #
	  // ###################################
	  if ((correct4Efficiency != "noEffCorr") && (correct4Efficiency != "yesEffCorr") && (correct4Efficiency != "yesEffCorrGenAnaly") && (correct4Efficiency != "yesEffCorrGenBin"))
	    {
	      cout << "[ExtractYield::main]\tIncorrect option parameter " << correct4Efficiency << endl;
	      exit (EXIT_FAILURE);
	    }


	  // #################################
	  // # Read the q^2 bin and the rest #
	  // #################################
	  if (argc >= 5) specBin = atoi(argv[4]);
	  if ((correct4Efficiency == "yesEffCorr") || (correct4Efficiency.find("yesEffCorrGen") != string::npos)) useEffPDF = true;
	  if (correct4Efficiency.find("yesEffCorrGen") != string::npos)
	    {
	      fileIndx    = atoi(argv[5]);
	      tmpFileName = argv[6];
	    }
	  else if ((!(((FitType >= 21) && (FitType <= 26)) || (FitType == 96))) &&
		   ((correct4Efficiency == "noEffCorr") || (correct4Efficiency == "yesEffCorr")))
	    {
	      if (argc >= 6)
		{
		  fileIndx = atoi(argv[5]);
		  if (argc == 7) tmpFileName = argv[6];
		}
	    }
	  else if ((FitType >= 21) && (FitType <= 26))
	    {
	      nToy          = atoi(argv[5]);
	      fileIndx      = atoi(argv[6]);
	      ParameterFILE = argv[7];
	    }
	  else if (FitType == 96) fileIndx = atoi(argv[5]);


	  cout << "\n[ExtractYield::main]\t@@@ Input variables from command line @@@" << endl;
	  cout << "- input/outputFile.root = " << fileName.c_str() << endl;
	  cout << "- correct4Efficiency = "    << correct4Efficiency << endl;
	  cout << "- tmpFileName = "           << tmpFileName.c_str() << endl;
	  cout << "- specBin = "               << specBin << endl;
	  cout << "- fileIndx = "              << fileIndx << endl;
	  cout << "- nToy = "                  << nToy << endl;
	  cout << "- FitType = "               << FitType << endl;
	  cout << "- useEffPDF = "             << useEffPDF << endl;
	  cout << "- ParameterFILE = "         << ParameterFILE.c_str() << endl;

	  cout << "\n[ExtractYield::main]\t@@@ Internal settings @@@" << endl;
	  cout << "NBINS = "         << NBINS << endl;
	  cout << "MULTYIELD = "     << MULTYIELD << endl;
	  cout << "NCOEFFPOLYBKG = " << NCOEFFPOLYBKG << endl;
	  cout << "DEGREEINTERP = "  << DEGREEINTERP << endl;
	  cout << "MAXTRIALS = "     << MAXTRIALS << endl;

	  cout << "\nMakeMuMuPlots = " << MakeMuMuPlots << endl;
	  cout << "USEMINOS = "        << USEMINOS << endl;
	  cout << "SETBATCH  = "       << SETBATCH << endl;
	  cout << "SAVEPOLY = "        << SAVEPOLY << endl;
	  cout << "SAVEPLOT = "        << SAVEPLOT << endl;
	  cout << "RESETsigANG = "     << RESETsigANG << endl;
	  cout << "RESETcomANG = "     << RESETcomANG << endl;
	  cout << "FUNCERRBAND = "     << FUNCERRBAND << endl;

	  cout << "\nPARAMETERFILEIN = " << PARAMETERFILEIN << endl;
	  cout << "PARAMETERFILEOUT = "  << PARAMETERFILEOUT << endl;


	  if (SETBATCH == true)
	    {
	      cout << "\n[ExtractYield::main]\t@@@ Setting batch mode @@@" << endl;
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
	  Utility->ReadGenericParam(ParameterFILE);
	  Utility->ReadSelectionCuts(ParameterFILE);
	  Utility->ReadFitStartingValues(ParameterFILE,&fitParam,&configParam,Utility->ParFileBlockN("fitValGlob"));

	  CTRLfitWRKflow = Utility->GetGenericParam("CtrlFitWrkFlow");

	  fileFitResults << "Slewrate for physics constraints shape: " << Utility->GetGenericParam("SlewRateConstr").c_str() << endl;
	  fileFitResults << "Normalize to J/psi and not psi(2S): "     << Utility->GetGenericParam("NormJPSInotPSIP").c_str() << " (0 = false; 1 = true)" << endl;
	  fileFitResults << "Apply constraints: "                      << Utility->GetGenericParam("ApplyConstr").c_str() << " (0 = false; 1 = true)" << endl;
	  fileFitResults << "Control fit workflow: "                   << CTRLfitWRKflow.c_str() << endl;
	  fileFitResults << "Control mis-tag fraction workflow: "      << Utility->GetGenericParam("CtrlMisTagWrkFlow").c_str() << endl;
	  fileFitResults << "Save mis-tag fraction: "                  << Utility->GetGenericParam("SaveMisTagFrac").c_str() << " (0 = false; 1 = true)" << endl;
	  

	  // ################################################
	  // # Read ctr. chn. yield from global-fit results #
	  // ################################################
	  myString.clear(); myString.str("");
	  myString << fitParam[Utility->GetFitParamIndx("nSig")]->operator[]((atoi(Utility->GetGenericParam("NormJPSInotPSIP").c_str()) == true ? 1 : 2)).c_str();
	  SetValueAndErrors(NULL,"",1.0,&myString,&PsiYieldGoodTag,&PsiYieldGoodTagErr,&PsiYieldGoodTagErr);
	  
	  double MisTagFrac, MisTagFracErr;
	  myString.clear(); myString.str("");
	  myString << fitParam[Utility->GetFitParamIndx("nMisTagFrac")]->operator[]((atoi(Utility->GetGenericParam("NormJPSInotPSIP").c_str()) == true ? 1 : 2)).c_str();
	  SetValueAndErrors(NULL,"",1.0,&myString,&MisTagFrac,&MisTagFracErr,&MisTagFracErr);

	  PsiYieldMisTag    = PsiYieldGoodTag    * MisTagFrac / (1. - MisTagFrac);
	  PsiYieldMisTagErr = PsiYieldGoodTagErr * MisTagFrac / (1. - MisTagFrac);

	  fileFitResults << "Normalization channel yield good-tagged events: " << PsiYieldGoodTag << " +/- " << PsiYieldGoodTagErr << endl;
	  fileFitResults << "Normalization channel yield mis-tagged events: "  << PsiYieldMisTag  << " +/- " << PsiYieldMisTagErr  << endl;

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
	  else              Utility->ReadFitStartingValues((((correct4Efficiency.find("yesEffCorrGen") == string::npos) && tmpFileName.size() != 0) ? tmpFileName : ParameterFILE),
							   &fitParam,&configParam,
							   (((correct4Efficiency.find("yesEffCorrGen") == string::npos) && tmpFileName.size() != 0) ? 1 : Utility->ParFileBlockN("fitValBins")));


	  // #############################################################
	  // # Make q^2 bins for histograms and fill efficiency matrices #
	  // #############################################################
	  q2BinsHisto = new double[q2Bins.size()];
	  for (unsigned int i = 0; i < q2Bins.size(); i++) q2BinsHisto[i] = q2Bins[i];


	  // ###################################
	  // # Reset signal angular parameters #
	  // ###################################
	  if (RESETsigANG == true) ResetAngularParam(&fitParam);


	  // #####################################################
	  // # Reset combinatorial background angular parameters #
	  // #####################################################
	  if (RESETcomANG == true) ResetCombPolyParam(&fitParam);


	  // ###########################################################################
	  // # Prepare file to save fit results in case of systematic error evaluation #
	  // ###########################################################################
	  myString.clear(); myString.str("");
	  if (specBin != -1) myString << "FitSystematics_" << specBin << ".txt";
	  else               myString << "FitSystematics.txt";
	  fileFitSystematics.open(myString.str().c_str(),ios_base::app);
	  if (fileFitSystematics.good() == false)
	    {
	      cout << "[ExtractYield::main]\tError opening file : " << myString.str().c_str() << endl;
	      CloseAllAndQuit(theApp,NtplFile);
	    }


	  // ##############################################
	  // # Read coefficents for analytical efficiency #
	  // ##############################################
	  effFuncs.first  = new vector<TF2*>;
	  effFuncs.second = new vector<TF2*>;
	  if (correct4Efficiency == "yesEffCorrGenAnaly")
	    {
	      Utility->ReadAnalyticalEff(tmpFileName.c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,effFuncs.first, "effFuncs",1);
	      Utility->ReadAnalyticalEff(tmpFileName.c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,effFuncs.second,"effFuncs",2);
	    }
	  else
	    {
	      if (correct4Efficiency == "yesEffCorrGenBin")
		{
		  // ####################
		  // # Remove file tail #
		  // ####################
		  tmpFileName.erase(tmpFileName.find(".txt") - (1 + static_cast<int>(log10(q2Bins.size()))),6 + static_cast<int>(log10(q2Bins.size())));

		  Utility->SetDirEfficiency("");
		  if ((FitType == 6) || (FitType == 46) || (FitType == 66)) Utility->SetHisto2DEffName(GetSignalType(FitType,&q2Bins,specBin),tmpFileName.c_str());
		}

	      Utility->ReadAnalyticalEff(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,effFuncs.first,"effFuncs",Utility->ParFileBlockN("analyEffokTag"));
	      Utility->ReadAnalyticalEff(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,effFuncs.second,"effFuncs",Utility->ParFileBlockN("analyEffmisTag"));
	    }


 	  // ###############################################################################################
	  // # Read other parameters : this also allow to understand if the parameter file is well written #
 	  // ###############################################################################################
	  LUMI = Utility->ReadLumi(ParameterFILE);
	  if (Utility->IsThisData(ParameterFILE) == true) cout << "\n[ExtractYield::main]\t@@@ I recognize that this is a DATA file @@@" << endl;
	  else                                            cout << "\n[ExtractYield::main]\t@@@ I recognize that this is a Monte Carlo file @@@" << endl;


	  // ###################
	  // # Select fit type #
	  // ###################
	  if (((FitType >= 1)  && (FitType <= 6))  ||
	      (FitType == 36)                      ||
	      ((FitType >= 41) && (FitType <= 46)) ||
	      (FitType == 56)                      ||
	      ((FitType >= 61) && (FitType <= 66)) ||
	      (FitType == 76))
	    {
	      NtplFile = new TFile(fileName.c_str(), "READ");
	      theTree  = (TTree*) NtplFile->Get("B0KstMuMu/B0KstMuMuNTuple");
	      NTuple   = new B0KstMuMuSingleCandTreeContent();
	      NTuple->Init();
	    

	      // #################
	      // # Make datasets #
	      // #################
	      cout << "\n[ExtractYield::main]\t@@@ Making datasets @@@" << endl;
	      MakeDataSets(NTuple,FitType);


	      // ##############################
	      // # Select the proper fit type #
	      // ##############################
	      if ((FitType == 1) || (FitType == 2) || (FitType == 41) || (FitType == 61))
		{
		  if (FitType == 2)
		    {
		      double NLLvalue;
		      vector<string>* vecParStr;


		      // ##########################
		      // # 1D-fit to B0 inv. mass #
		      // ##########################
		      cout << "\n[ExtractYield::main]\t@@@ Now fit total invariant mass NON resonant channel @@@" << endl;
		      TCanvas* cB0MassArbRejectPsi = new TCanvas("cB0MassArbRejectPsi","cB0MassArbRejectPsi",10, 10, 700, 500);
		      InstantiateMassFit(&TotalPDFRejectPsi,B0MassArb,"TotalPDFRejectPsi",&configParam,0);
		      CopyFitResults(TotalPDFRejectPsi,0,&fitParam);
		      ClearVars(&vecConstr);
		      BuildMassConstraints(&vecConstr,TotalPDFRejectPsi,"sign");
		      if ((strcmp(CTRLfitWRKflow.c_str(),"trueAll&NoFFrac") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"trueAll&FitFrac") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"allEvts") == 0)) BuildMassConstraints(&vecConstr,TotalPDFRejectPsi,"mistag");
		      PrintVariables(TotalPDFRejectPsi->getVariables(),"vars");
		      PrintVariables(&vecConstr,"cons");
		      MakeMassFit(SingleCandNTuple_RejectPsi,&TotalPDFRejectPsi,B0MassArb,cB0MassArbRejectPsi,&vecConstr,&NLLvalue,NULL,fileIndx);

		      // ##############################################
		      // # Save fit results back into prarameter file #
		      // ##############################################
		      vecParStr = SaveFitResults(TotalPDFRejectPsi,0,&fitParam,&configParam,&vecConstr);
		      Utility->SaveFitValues(PARAMETERFILEOUT,vecParStr,-1,"  signal   ");
		      vecParStr->clear();
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


		      cout << "\n[ExtractYield::main]\t@@@ Now fit total invariant mass resonant J/psi channel @@@" << endl;
		      TCanvas* cB0MassArbJPsi = new TCanvas("cB0MassArbJPsi","cB0MassArbJPsi",10, 10, 700, 500);
		      InstantiateMassFit(&TotalPDFPsi,B0MassArb,"TotalPDFPsi",&configParam,1);
		      CopyFitResults(TotalPDFPsi,1,&fitParam);
		      ClearVars(&vecConstr);
		      BuildMassConstraints(&vecConstr,TotalPDFPsi,"sign");
		      if ((strcmp(CTRLfitWRKflow.c_str(),"trueAll&NoFFrac") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"trueAll&FitFrac") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"allEvts") == 0)) BuildMassConstraints(&vecConstr,TotalPDFPsi,"mistag");
		      PrintVariables(TotalPDFPsi->getVariables(),"vars");
		      PrintVariables(&vecConstr,"cons");
		      MakeMassFit(SingleCandNTuple_JPsi,&TotalPDFPsi,B0MassArb,cB0MassArbJPsi,&vecConstr,&NLLvalue,NULL,fileIndx);

		      // ##############################################
		      // # Save fit results back into prarameter file #
		      // ##############################################
		      vecParStr = SaveFitResults(TotalPDFPsi,1,&fitParam,&configParam,&vecConstr);
		      Utility->SaveFitValues(PARAMETERFILEOUT,vecParStr,-1,"   J/psi   ");
		      vecParStr->clear();
		      delete vecParStr;


		      cout << "\n[ExtractYield::main]\t@@@ Now fit total invariant mass resonant psi(2S) channel @@@" << endl;
		      TCanvas* cB0MassArbPsiP = new TCanvas("cB0MassArbPsiP","cB0MassArbPsiP",10, 10, 700, 500);
		      InstantiateMassFit(&TotalPDFPsi,B0MassArb,"TotalPDFPsi",&configParam,2);
		      CopyFitResults(TotalPDFPsi,2,&fitParam);
		      ClearVars(&vecConstr);
		      BuildMassConstraints(&vecConstr,TotalPDFPsi,"sign");
		      if ((strcmp(CTRLfitWRKflow.c_str(),"trueAll&NoFFrac") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"trueAll&FitFrac") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"allEvts") == 0)) BuildMassConstraints(&vecConstr,TotalPDFPsi,"mistag");
		      PrintVariables(TotalPDFPsi->getVariables(),"vars");
		      PrintVariables(&vecConstr,"cons");
		      MakeMassFit(SingleCandNTuple_PsiP,&TotalPDFPsi,B0MassArb,cB0MassArbPsiP,&vecConstr,&NLLvalue,NULL,fileIndx);

		      // ##############################################
		      // # Save fit results back into prarameter file #
		      // ##############################################
		      vecParStr = SaveFitResults(TotalPDFPsi,2,&fitParam,&configParam,&vecConstr);
		      Utility->SaveFitValues(PARAMETERFILEOUT,vecParStr,-1,"  psi(2S)  ");
		      vecParStr->clear();
		      delete vecParStr;
		    }
		  else
		    {
		      // ######################################
		      // # 1D-fit to B0 inv. mass per q^2 bin #
		      // ######################################
		      TCanvas* cHistoMeas = new TCanvas("cHistoMeas","cHistoMeas",10, 10, 900, 600);
		      cHistoMeas->Divide(2,2);

		      vector<TH1D*> VecHistoMeas;

		      cHistoMeas->cd(1);
		      VecHistoMeas.push_back(new TH1D("histoMeas0", "B0 --> K*0(K+ pi-) mu+ mu- / B0 --> K*0(K+ pi-) J/psi(mu+ mu-) Branching Fraction", q2Bins.size()-1, q2BinsHisto));
		      VecHistoMeas[0]->SetStats(false);
		      VecHistoMeas[0]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		      VecHistoMeas[0]->SetYTitle("dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}7}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");
		      VecHistoMeas[0]->GetYaxis()->SetRangeUser(0.0,1.2);

		      cout << "\n[ExtractYield::main]\t@@@ Now fit invariant mass per mumu q^2 bins @@@" << endl;
		      if (FitType == 1) IterativeMassFitq2Bins(SingleCandNTuple_RejectPsi,
							       useEffPDF,
							       PsiYieldGoodTag,PsiYieldGoodTagErr,
							       PsiYieldMisTag,
							       B0MassArb,
							       specBin,
							       FitType,
							       &VecHistoMeas,
							       &q2Bins,
							       &configParam,&fitParam,
							       &vecConstr,
							       fileIndx);
		      else if (FitType == 41) IterativeMassFitq2Bins(SingleCandNTuple_JPsi,
								     useEffPDF,
								     PsiYieldGoodTag,PsiYieldGoodTagErr,
								     PsiYieldMisTag,
								     B0MassArb,
								     specBin,
								     FitType,
								     &VecHistoMeas,
								     &q2Bins,
								     &configParam,&fitParam,
								     &vecConstr,
								     fileIndx);
		      else IterativeMassFitq2Bins(SingleCandNTuple_PsiP,
						  useEffPDF,
						  PsiYieldGoodTag,PsiYieldGoodTagErr,
						  PsiYieldMisTag,
						  B0MassArb,
						  specBin,
						  FitType,
						  &VecHistoMeas,
						  &q2Bins,
						  &configParam,&fitParam,
						  &vecConstr,
						  fileIndx);

		      cHistoMeas->cd(1);
		      VecHistoMeas[0]->SetMarkerStyle(20);
		      VecHistoMeas[0]->Draw("pe1");
		      DrawString(LUMI);

		      cHistoMeas->Modified();
		      cHistoMeas->Update();
		    }
		}
	      else if ((FitType == 6) || (FitType == 36) || (FitType == 46) || (FitType == 56) || (FitType == 66) || (FitType == 76)) // Fl-Afb-fit
		{
		  // #############################
		  // # 3D-fit Afb-Fl per q^2 bin #
		  // #############################
		  TCanvas* cHistoMeas = new TCanvas("cHistoMeas","cHistoMeas",10, 10, 900, 600);
		  cHistoMeas->Divide(2,2);

		  vector<TH1D*> VecHistoMeas;

		  cHistoMeas->cd(1);
		  VecHistoMeas.push_back(new TH1D("histoMeas0", "B0 --> K*0(K+ pi-) mu+ mu- : FL vs dimuon q2", q2Bins.size()-1, q2BinsHisto));
		  VecHistoMeas[0]->SetStats(false);
		  VecHistoMeas[0]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		  VecHistoMeas[0]->SetYTitle("F_{L}");
		  VecHistoMeas[0]->GetYaxis()->SetRangeUser(-0.01,1.01);

		  cHistoMeas->cd(1);
		  VecHistoMeas.push_back(new TH1D("histoMeas1", "B0 --> K*0(K+ pi-) mu+ mu- : AFB vs dimuon q2", q2Bins.size()-1, q2BinsHisto));
		  VecHistoMeas[1]->SetStats(false);
		  VecHistoMeas[1]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		  VecHistoMeas[1]->SetYTitle("A_{FB}");
		  VecHistoMeas[1]->GetYaxis()->SetRangeUser(-1.0 - 0.01,1.0 + 0.01);

		  cHistoMeas->cd(2);
		  VecHistoMeas.push_back(new TH1D("histoMeas2", "B0 --> K*0(K+ pi-) mu+ mu- : dBF/dq2 vs dimuon q2", q2Bins.size()-1, q2BinsHisto));
		  VecHistoMeas[2]->SetStats(false);
		  VecHistoMeas[2]->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
		  VecHistoMeas[2]->SetYTitle("dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}7}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");
		  VecHistoMeas[2]->GetYaxis()->SetRangeUser(0.0,1.2);

		  cout << "\n[ExtractYield::main]\t@@@ Now fit invariant mass, cos(theta_K) and cos(theta_l) per mumu q^2 bins @@@" << endl;
		  if ((FitType == 6) || (FitType == 36)) IterativeMass2AnglesFitq2Bins(SingleCandNTuple_RejectPsi,
		  								       useEffPDF,
										       PsiYieldGoodTag,PsiYieldGoodTagErr,
										       PsiYieldMisTag,
		  								       B0MassArb,
		  								       CosThetaMuArb,
		  								       CosThetaKArb,
		  								       specBin,
		  								       FitType,
		  								       &VecHistoMeas,
		  								       &q2Bins,
		  								       &configParam,&fitParam,
		  								       effFuncs,
		  								       &vecConstr,
		  								       fileIndx);
		  else if ((FitType == 46) || (FitType == 56)) IterativeMass2AnglesFitq2Bins(SingleCandNTuple_JPsi,
		  									     useEffPDF,
											     PsiYieldGoodTag,PsiYieldGoodTagErr,
											     PsiYieldMisTag,
		  									     B0MassArb,
		  									     CosThetaMuArb,
		  									     CosThetaKArb,
		  									     specBin,
		  									     FitType,
		  									     &VecHistoMeas,
		  									     &q2Bins,
		  									     &configParam,&fitParam,
		  									     effFuncs,
		  									     &vecConstr,
		  									     fileIndx);
		  else IterativeMass2AnglesFitq2Bins(SingleCandNTuple_PsiP,
		  				     useEffPDF,
						     PsiYieldGoodTag,PsiYieldGoodTagErr,
						     PsiYieldMisTag,
		  				     B0MassArb,
		  				     CosThetaMuArb,
		  				     CosThetaKArb,
		  				     specBin,
		  				     FitType,
		  				     &VecHistoMeas,
		  				     &q2Bins,
		  				     &configParam,&fitParam,
		  				     effFuncs,
		  				     &vecConstr,
		  				     fileIndx);

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

		  cHistoMeas->Modified();
		  cHistoMeas->Update();
		}
	    }
	  else if ((FitType >= 21) && (FitType <= 26))
	    {
	      // #################
	      // # Make datasets #
	      // #################
	      cout << "\n[ExtractYield::main]\t@@@ Making datasets @@@" << endl;
	      MakeDataSets(NTuple,FitType);


	      // ##########################################
	      // # Set seed for random number generator   #
	      // # in case the toy-MC studies are splited #
	      // ##########################################
	      RooRandom::randomGenerator()->SetSeed(fileIndx*(q2Bins.size()-1) + specBin + 1);
	      cout << "\n[ExtractYield::main]\t@@@ Random seed for toy-MC set to : " << RooRandom::randomGenerator()->GetSeed() << " @@@" << endl;


	      if (FitType == 21)
		{
		  // ####################################
		  // # Make toy-MC for B0 inv. mass fit #
		  // ####################################
		  cout << "\n[ExtractYield::main]\t@@@ Now make TOY-MC for fit to B0 total invariant mass @@@" << endl;
		  TCanvas* cToyMC = new TCanvas("cToyMC","cToyMC",10, 10, 1200, 800);

		  InstantiateMassFit(&TotalPDFRejectPsi,B0MassArb,"TotalPDFRejectPsi",&configParam,specBin);
		  MakeMassToy(TotalPDFRejectPsi,B0MassArb,cToyMC,nToy,specBin,&fitParam,&vecConstr,fileName);
		}
	      else if (FitType == 26) // Fl-Afb-fit
		{
		  // ###############################################
		  // # Make toy-MC for B0 inv. mass fit and angles #
		  // ###############################################
		  cout << "\n[ExtractYield::main]\t@@@ Now make TOY-MC for fit to B0 total invariant mass and cos(theta_K) and cos(theta_l) @@@" << endl;
		  TCanvas* cToyMC = new TCanvas("cToyMC","cToyMC",10, 10, 1200, 800);

		  InstantiateMass2AnglesFit(&TotalPDFRejectPsi,useEffPDF,B0MassArb,CosThetaMuArb,CosThetaKArb,"TotalPDFRejectPsi",FitType,&configParam,&fitParam,&q2Bins,specBin,specBin,make_pair(effFuncs.first->operator[](specBin),effFuncs.second->operator[](specBin)));
		  MakeMass2AnglesToy(TotalPDFRejectPsi,B0MassArb,CosThetaMuArb,CosThetaKArb,cToyMC,FitType,nToy,specBin,&fitParam,&vecConstr,fileName,&q2Bins);
		}
	    }
	  else if ((FitType >= 81) && (FitType <= 86))
	    {
	      // #################
	      // # Make datasets #
	      // #################
	      cout << "\n[ExtractYield::main]\t@@@ Making datasets @@@" << endl;
	      MakeDataSets(NTuple,FitType);


	      if (FitType == 81)
                {
                  // ###########################################                                                                                                                                                                                                                                                     
                  // # Generate dataset for B0 inv. mass model #                                                                                                                                                                                                                                                     
                  // ###########################################                                                                                                                                                                                                                                                     
                  InstantiateMassFit(&TotalPDFRejectPsi,B0MassArb,"TotalPDFRejectPsi",&configParam,specBin);
                  GenerateDataset(TotalPDFRejectPsi,RooArgSet(*B0MassArb),&q2Bins,specBin,&fitParam,fileName);
                }
              else if (FitType == 86) // Fl-Afb-fit                                                                                                                                                                                                                                                                  
                {
                  // ######################################################                                                                                                                                                                                                                                          
                  // # Generate dataset for B0 inv. mass and angles model #                                                                                                                                                                                                                                          
                  // ######################################################                                                                                                                                                                                                                                          
                  InstantiateMass2AnglesFit(&TotalPDFRejectPsi,useEffPDF,B0MassArb,CosThetaMuArb,CosThetaKArb,"TotalPDFRejectPsi",FitType,&configParam,&fitParam,&q2Bins,specBin,specBin,make_pair(effFuncs.first->operator[](specBin),effFuncs.second->operator[](specBin)));
                  GenerateDataset(TotalPDFRejectPsi,RooArgSet(*B0MassArb,*CosThetaMuArb,*CosThetaKArb),&q2Bins,specBin,&fitParam,fileName);
                }
	    }
	  else if (FitType == 96)
	    {
	      // #################
	      // # Make datasets #
	      // #################
	      cout << "\n[ExtractYield::main]\t@@@ Making datasets @@@" << endl;
	      MakeDataSets(NTuple,FitType);


	      if (FitType == 96) // Fl-Afb-fit
                {
                  // #################################################################
                  // # Generate new parameter file for B0 inv. mass and angles model #
                  // #################################################################
                  InstantiateMass2AnglesFit(&TotalPDFRejectPsi,useEffPDF,B0MassArb,CosThetaMuArb,CosThetaKArb,"TotalPDFRejectPsi",FitType,&configParam,&fitParam,&q2Bins,specBin,specBin,make_pair(effFuncs.first->operator[](specBin),effFuncs.second->operator[](specBin)));
                }


	      vector<string>* vecParStr;

	      myString.clear(); myString.str("");
	      myString << "_" << specBin << "_" << fileIndx << ".txt";
	      fileName.replace(fileName.find(".root"),5,myString.str());

	      myString.clear(); myString.str("");
	      myString << (Utility->nFitParam + Utility->nConfigParam) * (q2Bins.size()-1);
	      vecParStr = new vector<string>;
	      vecParStr->push_back(myString.str().c_str());
	      Utility->SaveFitValues(fileName,vecParStr,-2,"# Number of fit parameters");

	      vecParStr->clear();
	      delete vecParStr;

	      for (unsigned int i = 0; i < q2Bins.size()-1; i++)
		{
		  if ((specBin != -1) && (i != static_cast<unsigned int>(specBin))) vecParStr = SaveFitResults(NULL,i,&fitParam,&configParam,NULL);
		  else
		    {
		      GenerateFitParameters(TotalPDFRejectPsi,&fitParam,fileIndx,&q2Bins,i,"All"); // @TMP@ : "All" "misTagFrac" "FlAfb" "bkgAng"
		      vecParStr = SaveFitResults(TotalPDFRejectPsi,i,&fitParam,&configParam,&vecConstr);
		    }

		  Utility->SaveFitValues(fileName,vecParStr,i);

		  vecParStr->clear();
		  delete vecParStr;
		}
	    }


	  fileFitResults.close();
	  fileFitSystematics.close();
	  if (SETBATCH == true)
	    {
	      cout << "Bye bye !" << endl;
	      CloseAllAndQuit(theApp,NtplFile);
	    }
	  else
	    {
	      system("say \" Let's rock and roll ! \"");
 	      theApp->Run (); // Eventloop on air
	    }
	}
      else
	{
	  cout << "Wrong parameter: " << endl;
	  cout << "./ExtractYield [FitType] [input/output[if toy-MC]File.root] [noEffCorr yesEffCorr yesEffCorrGenAnaly yesEffCorrGenBin]" << endl;
	  cout << "               [q^2 bin to fit (0 - ...)]" << endl;
	  cout << "               [[if yesEffCorrGen]indx AND effFileName.txt]" << endl;
	  cout << "               [[if toy-MC]nToy AND indx AND ParameterFile.txt] " << endl;
	  cout << "               [[if 96]indx]" << endl;
	  cout << "               [indx] [ParameterFile.txt]" << endl;

	  cout << "\n --> noEffCorr     = no eff. correction" << endl;
	  cout << " --> yesEffCorr    = use eff. correction" << endl;
	  cout << " --> yesEffCorrGen = compute systematic error associated to eff. correction" << endl;

	  return EXIT_FAILURE;
	}
    }
  else
    {
      cout << "Parameter missing: " << endl;
      cout << "./ExtractYield [FitType] [input/output[if toy-MC]File.root] [noEffCorr yesEffCorr yesEffCorrGenAnaly yesEffCorrGenBin]" << endl;
      cout << "               [q^2 bin to fit (0 - ...)]" << endl;
      cout << "               [[if yesEffCorrGen]indx AND effFileName.txt]" << endl;
      cout << "               [[if toy-MC]nToy AND indx AND ParameterFile.txt] " << endl;
      cout << "               [[if 96]indx]" << endl;
      cout << "               [indx] [ParameterFile.txt]" << endl;

      cout << "\n --> noEffCorr     = no eff. correction" << endl;
      cout << " --> yesEffCorr    = use eff. correction" << endl;
      cout << " --> yesEffCorrGen = compute systematic error associated to eff. correction" << endl;

      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  Signa  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 1: 1D branching fraction per q^2 bin" << endl;
      cout << "FitType = 2: 1D (B0Mass) peak" << endl;
      cout << "FitType = 6: 3D Afb-Fl (B0Mass, cos(theta_K), cos(theta_l)) per q^2 bin" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 21: toy-MC 1D (B0Mass) peak per q^2 bin" << endl;
      cout << "FitType = 26: toy-MC 3D Afb-Fl (B0Mass, cos(theta_K), cos(theta_l)) per q^2 bin" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 36: 2D Afb-Fl (cos(theta_K), cos(theta_l)) per q^2 bin on GEN variables" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  J/psi  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 41: 1D (B0Mass) peak for B0 --> J/psi(mumu) K*0 in appropriate q^2 bin" << endl;
      cout << "FitType = 46: 3D Afb-Fl (B0Mass, cos(theta_K), cos(theta_l)) for B0 --> J/psi(mumu) K*0" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 56: 2D Afb-Fl (cos(theta_K), cos(theta_l)) for B0 --> J/psi(mumu) K*0 on GEN variables" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Psi(2S) @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 61: 1D (B0Mass) peak for B0 --> psi(2S)(mumu) K*0 in appropriate q^2 bin" << endl;
      cout << "FitType = 66: 3D Afb-Fl (B0Mass, cos(theta_K), cos(theta_l)) for B0 --> psi(2S)(mumu) K*0" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 76: 2D Afb-Fl (cos(theta_K), cos(theta_l)) for B0 --> psi(2S)(mumu) K*0 on GEN variables" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Genera dataset @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 81: 1D generare dataset with (B0Mass) distribution from pdf" << endl;
      cout << "FitType = 86: 3D generare dataset with (B0Mass, cos(theta_K), cos(theta_l)) distribution from pdf" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      
      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Genera parameter @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 96: 3D generate paramter with (B0Mass, cos(theta_K), cos(theta_l)) distribution from pdf" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

// #########################################################################
// # Program to compute the efficiency for the B0 --> K*0 mu+ mu- analysis #
// # in bins of dimuon q^2, cos(theta_K), cos(theta_l), and phi            #
// #########################################################################
// # Author: Mauro Dinardo                                                 #
// #########################################################################

#include <TROOT.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TF12.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TText.h>
#include <TFitResult.h>

#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

#include "ReadParameters.h"
#include "Utils.h"
#include "B0KstMuMuSingleCandTreeContent.h"

using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::vector;


// ###########################################
// # How to create the analytical efficiency #
// ###########################################
// (a) Make binned efficiency with the command "Make"

// (b) Fit the theta_l variable with the command "Fit1DEff" for every q2 bin
//     - Set command option to "thetaL"

// (c) Fit the theta_K variable with the command "Fit1DEff"
//     - Set command option to "thetaK"
//     - Set "INPUT_THETAL" to the output.txt file of point (b)

// (d) Fit the phi variable with the command "Fit1DEff" for every q2 bin
//     - Set command option to "phi"

// (e) Fit the theta_l-theta_K variables with the command "FitEff2D"
//     - Set "INPUT_THETAL_THETAK" to the output.txt file of point (c)

// (f) Fit the theta_l-theta_K-phi variables with the command "FitEff3D"
//     - Set "INPUT_THETAL_THETAK" to the output.txt file of point (c)
//     - Set "INPUT_PHI"           to the output.txt file of point (d)

// (g) If you want to look at the final result run the command "Test2DEff"
//     - Set "INPUT_THETAL_THETAK" to the output.txt file of point (e)

// (h) If you want to look at the final result run the command "Test3DEff"
//     - Set "INPUT_THETAL_THETAK_PHI" to the output.txt file of point (f)


// ####################
// # Global constants #
// ####################
#define INPUT_THETAL           "ThetaL_B0ToKstMuMu.txt"
#define INPUT_PHI              "Phi_B0ToKstMuMu.txt"
#define INPUT_THETAL_THETAK    "ThetaK_B0ToKstMuMu.txt"
#define TEST_THETAL_THETAK     "ThetaKThetaL_B0ToKstMuMu.txt"
#define TEST_THETAL_THETAK_PHI "ThetaKThetaLPhi_B0ToKstMuMu.txt"

#define RIGHTtag       true
#define SavePlot       false
#define CHECKEFFatREAD false // Check if 2D or 3D efficiency go negative
#define NFILES         200
#define INPUTGenEff    "../efficiency/EffRndGenAnalyFilesSign_JPsi_Psi2S/Efficiency_RndGen.txt"
#define SETBATCH       true // Set batch mode when making binned efficiency
#define ParameterFILE  "../python/ParameterFile.txt"
#define ordinateRange  1e-2


// ####################
// # Global variables #
// ####################
Utils* Utility;

TTree* theTreeGenCandidatesNoFilter;
TTree* theTreeRecoCandidates;
TTree* theTreeSingleCand;
B0KstMuMuSingleCandTreeContent* NTupleGenCandidatesNoFilter;
B0KstMuMuSingleCandTreeContent* NTupleRecoCandidates;
B0KstMuMuSingleCandTreeContent* NTupleSingleCand;

vector<TF2*> effFuncs;
vector<TMatrixTSym<double>*> covMatrices;

vector<double> q2Bins;
vector<double> cosThetaKBins;
vector<double> cosThetaLBins;
vector<double> phiBins;

Utils::effStruct myEff;


// #######################
// # Function Definition #
// #######################
void ComputeEfficiency     (TTree* theTree, B0KstMuMuSingleCandTreeContent* NTuple, double* Vector, double* VectorErr2Pois, double* VectorErr2Weig, unsigned int type,
			    vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, int SignalType);
void MakeHistogramsAllBins (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, Utils::effStruct myEff);
void Read3DEfficiencies    (bool isSingleEff, vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
			    string fileNameInput, bool isAnalyEff, Utils::effStruct* myEff, bool CheckEffatRead, bool savePlot, int specBin = -1);
void Fit1DEfficiencies     (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
			    unsigned int SignalType, Utils::effStruct myEff, string who, unsigned int q2BinIndx, string fileNameOut);
void Fit2DEfficiencies     (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
			    unsigned int SignalType, Utils::effStruct myEff, unsigned int q2BinIndx, string fileNameOut);
void Fit3DEfficiencies     (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
			    unsigned int SignalType, Utils::effStruct myEff, unsigned int q2BinIndx, string fileNameOut);
void Test2DEfficiency      (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, Utils::effStruct myEff, unsigned int q2BinIndx, bool savePlot);
void Test3DEfficiency      (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, Utils::effStruct myEff, unsigned int q2BinIndx, bool savePlot);


// ###########################
// # Function Implementation #
// ###########################
void ComputeEfficiency (TTree* theTree, B0KstMuMuSingleCandTreeContent* NTuple, double* Vector, double* VectorErr2Pois, double* VectorErr2Weig, unsigned int type,
			vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, int SignalType)
// ##########################################################
// # Efficiency type = 1 --> total gen events before filter #
// # Efficiency type = 2 --> total gen events after filter  #
// # Efficiency type = 3 --> total reco events              #
// # Efficiency type = 4 --> total single candidate events  #
// ##########################################################
{
  // ###################
  // # Local variables #
  // ###################
  int nEntries;
  Utils::effStruct _Counter;
  double* Counter;
  double mumuq2;
  double cosThetaK;
  double cosThetaMu;
  double phiKstMuMuPlane;
  int mumuq2BinIndx;
  int cosThetaKBinIndx;
  int cosThetaMuBinIndx;
  int phiKstMuMuPlaneBinIndx;
  // ###################


  // ###################################
  // # Initialize efficiency structure #
  // ###################################
  Utility->GenEfficiency(&_Counter,q2Bins,cosThetaKBins,cosThetaLBins,phiBins);
  Utils::effValue myEffVal;
  Utility->ResetEffValue(&myEffVal,0.0);
  Utility->InitEfficiency(myEffVal,_Counter,q2Bins,cosThetaKBins,cosThetaLBins,phiBins);
  Counter = _Counter.Num1;


  NTuple->ClearNTuple();
  NTuple->SetBranchAddresses(theTree);
  nEntries = theTree->GetEntries();
  cout << "\n@@@ Computing efficiency type " << type;
  if      (type == 1) cout << " (before filter) @@@" << endl;
  else if (type == 2) cout << " (after filter) @@@" << endl;
  else if (type == 3) cout << " (reco events) @@@" << endl;
  else if (type == 4) cout << " (single candidate events) @@@" << endl;
  cout << "@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;


  for (int entry = 0; entry < nEntries; entry++)
    {
      theTree->GetEntry(entry);

      if ((NTuple->B0pT        > Utility->GetSeleCut("B0pT"))  &&
	  (fabs(NTuple->B0Eta) < Utility->GetSeleCut("B0Eta")) &&
	  
	  ((NTuple->genSignal == SignalType || NTuple->genSignal == SignalType+1)) &&
	  
	  ((type == 1 || type == 3) ||
	   
	   ((type == 2) &&
	    (sqrt(NTuple->genMumPx*NTuple->genMumPx + NTuple->genMumPy*NTuple->genMumPy)   > Utility->GetPreCut("MinMupT")) &&
	    (sqrt(NTuple->genMupPx*NTuple->genMupPx + NTuple->genMupPy*NTuple->genMupPy)   > Utility->GetPreCut("MinMupT")) &&
	    (fabs(Utility->computeEta(NTuple->genMumPx,NTuple->genMumPy,NTuple->genMumPz)) < Utility->GetPreCut("MuEta"))   &&
	    (fabs(Utility->computeEta(NTuple->genMupPx,NTuple->genMupPy,NTuple->genMupPz)) < Utility->GetPreCut("MuEta")))  ||

	   ((type == 4) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == Utility->RIGHTflavorTAG))))
	{
	  mumuq2          = NTuple->mumuMass->at(0)*NTuple->mumuMass->at(0);

	  if (NTuple->rightFlavorTag == true)
	    {
	      cosThetaK  = NTuple->CosThetaKArb;
	      cosThetaMu = NTuple->CosThetaMuArb;
	    }
	  else
	    {
	      cosThetaK  = - NTuple->CosThetaKArb;
	      cosThetaMu = - NTuple->CosThetaMuArb;
	    }
	  phiKstMuMuPlane = NTuple->PhiKstMuMuPlaneArb;

	  mumuq2BinIndx          = Utility->SearchBin(mumuq2,q2Bins);
	  cosThetaKBinIndx       = Utility->SearchBin(cosThetaK,cosThetaKBins);
	  cosThetaMuBinIndx      = Utility->SearchBin(cosThetaMu,cosThetaLBins);
	  phiKstMuMuPlaneBinIndx = Utility->SearchBin(phiKstMuMuPlane,phiBins);

	  if ((mumuq2BinIndx != -1) && (cosThetaKBinIndx != -1) && (cosThetaMuBinIndx != -1) && (phiKstMuMuPlaneBinIndx != -1))
	    {
	      if ((type == 1) || (type == 2))
		// ##################################################
		// # No event-weight for muon acceptance efficiency #
		// ##################################################
		Vector[phiKstMuMuPlaneBinIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
		       cosThetaMuBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
		       cosThetaKBinIndx*(q2Bins->size()-1) +
		       mumuq2BinIndx]++;
	      else if ((type == 3) || (type == 4))
		{
		  // ########################################
		  // # Add event-weight for reco efficiency #
		  // ########################################
		  Vector[phiKstMuMuPlaneBinIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
		  	 cosThetaMuBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
		  	 cosThetaKBinIndx*(q2Bins->size()-1) +
		  	 mumuq2BinIndx] +=
		    NTuple->evWeight;
	  
		  Counter[phiKstMuMuPlaneBinIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
			  cosThetaMuBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
			  cosThetaKBinIndx*(q2Bins->size()-1) +
			  mumuq2BinIndx]++;

		  VectorErr2Weig[phiKstMuMuPlaneBinIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaMuBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaKBinIndx*(q2Bins->size()-1) +
				 mumuq2BinIndx] +=
		    NTuple->evWeightE2;
		}
	    }
	}
    }


  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	for (unsigned int l = 0; l < phiBins->size()-1; l++)
	  {
	    if ((type == 1) || (type == 2))
	      {
		VectorErr2Pois[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] =
		  Vector[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      }

	    else if (Counter[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] <= 0.0)
	      {
		VectorErr2Pois[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = 0.0;
		VectorErr2Weig[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = 0.0;
	      }

	    else if ((type == 3) || (type == 4))
	      VectorErr2Pois[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] =
		Counter[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] *
		pow(Vector[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] /
		    Counter[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0);
	  }


  Utility->DeleteEfficiency(_Counter);
}


void MakeHistogramsAllBins (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, Utils::effStruct myEff)
{
  // ###################
  // # Local variables #
  // ###################
  stringstream myString;
  double Eff, EffErr;
  double totalEffAll;
  double totalEffSignal;
  double totalEffJPsi;
  double totalEffPsiP;
  double* q2Bins_        = Utility->MakeBinning(q2Bins);
  double* cosThetaKBins_ = Utility->MakeBinning(cosThetaKBins);
  double* cosThetaLBins_ = Utility->MakeBinning(cosThetaLBins);
  double* phiBins_       = Utility->MakeBinning(phiBins);
  // ###################
  string NumOrDen2Plot = "N1"; // It can be: "N1", "N2", "D1", "D2"
  double Yaxes = 2e5;
  double Zaxes = 1e3;
  // ###################


  // ##########################
  // # Efficiency projections #
  // ##########################
  TCanvas* cEff = new TCanvas("cEff", "cEff", 10, 10, 1200, 800);
  cEff->Divide(2,2);

  TH1D* Hq2 = new TH1D("Hq2", "Hq2", q2Bins->size()-1, q2Bins_);
  Hq2->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
  Hq2->SetYTitle("Efficiency");

  vector<TH1D*> vecHcosThetaK;
  vector<TH1D*> vecHcosThetaL;
  vector<TH1D*> vecHphi;

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    {
      myString.clear(); myString.str("");
      myString << "vecHcosThetaK_" << i;
      vecHcosThetaK.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), cosThetaKBins->size()-1, cosThetaKBins_));
      vecHcosThetaK.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
      vecHcosThetaK.back()->SetYTitle("Efficiency");

      myString.clear(); myString.str("");
      myString << "vecHcosThetaL_" << i;
      vecHcosThetaL.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), cosThetaLBins->size()-1, cosThetaLBins_));
      vecHcosThetaL.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
      vecHcosThetaL.back()->SetYTitle("Efficiency");

      myString.clear(); myString.str("");
      myString << "vecHphi_" << i;
      vecHphi.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), phiBins->size()-1, phiBins_));
      vecHphi.back()->SetXTitle("#phi");
      vecHphi.back()->SetYTitle("Efficiency");
    }


  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    {
      // #############################################
      // # Fill histogram : efficiency vs dimuon q^2 #
      // #############################################
      Utility->IntegrateEffButq2(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,i,myEff,&Eff,&EffErr);
      Hq2->SetBinContent(i+1,Eff);
      Hq2->SetBinError(i+1,EffErr);


      // ########################################################################
      // # Fill histogram : efficiency vs cos(theta_K) for different dimuon q^2 #
      // ########################################################################
      for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	{
	  Utility->IntegrateEffPhiCosThetaL(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,i,j,myEff,&Eff,&EffErr);
	  vecHcosThetaK[i]->SetBinContent(j+1,Eff);
	  vecHcosThetaK[i]->SetBinError(j+1,EffErr);
	}


      // ########################################################################
      // # Fill histogram : efficiency vs cos(theta_l) for different dimuon q^2 #
      // ########################################################################
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	{
	  Utility->IntegrateEffPhiCosThetaK(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,i,k,myEff,&Eff,&EffErr);
	  vecHcosThetaL[i]->SetBinContent(k+1,Eff);
	  vecHcosThetaL[i]->SetBinError(k+1,EffErr);
	}


      // ###############################################################
      // # Fill histogram : efficiency vs phi for different dimuon q^2 #
      // ###############################################################
      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  Utility->IntegrateEffCosThetaKCosThetaL(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,i,l,myEff,&Eff,&EffErr);
	  vecHphi[i]->SetBinContent(l+1,Eff);
	  vecHphi[i]->SetBinError(l+1,EffErr);
	}
    }


  cEff->cd(1);
  Hq2->SetMarkerStyle(20);
  Hq2->GetYaxis()->SetRangeUser(0.0,ordinateRange);
  Hq2->Draw("e1");

  cEff->cd(2);
  vecHcosThetaK[0]->SetMarkerStyle(20);
  vecHcosThetaK[0]->SetMarkerColor(1);
  vecHcosThetaK[0]->SetLineColor(1);
  vecHcosThetaK[0]->SetLineWidth(2);
  vecHcosThetaK[0]->Draw("e1");
  vecHcosThetaK[0]->GetYaxis()->SetRangeUser(0.0,ordinateRange);
  TLegend* legThetaK = new TLegend(0.88, 0.65, 0.97, 0.89, "");
  myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin 0";
  legThetaK->AddEntry(vecHcosThetaK[0],myString.str().c_str());
  for (unsigned int i = 1; i < q2Bins->size()-1; i++)
    {
      vecHcosThetaK[i]->SetMarkerStyle(20+i);
      vecHcosThetaK[i]->SetMarkerColor(1+i);
      vecHcosThetaK[i]->SetLineColor(1+i);
      vecHcosThetaK[i]->SetLineWidth(2);
      vecHcosThetaK[i]->Draw("sames e1");
      vecHcosThetaK[i]->GetYaxis()->SetRangeUser(0.0,ordinateRange);
      myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin " << i;
      legThetaK->AddEntry(vecHcosThetaK[i],myString.str().c_str());
    }
  legThetaK->SetFillColor(0);
  legThetaK->SetBorderSize(0);
  legThetaK->Draw();

  cEff->cd(3);
  vecHcosThetaL[0]->SetMarkerStyle(20);
  vecHcosThetaL[0]->SetMarkerColor(1);
  vecHcosThetaL[0]->SetLineColor(1);
  vecHcosThetaL[0]->SetLineWidth(2);
  vecHcosThetaL[0]->Draw("e1");
  vecHcosThetaL[0]->GetYaxis()->SetRangeUser(0.0,ordinateRange);
  TLegend* legThetaL = new TLegend(0.88, 0.65, 0.97, 0.89, "");
  myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin 0";
  legThetaL->AddEntry(vecHcosThetaL[0],myString.str().c_str());
  for (unsigned int i = 1; i < q2Bins->size()-1; i++)
    {
      vecHcosThetaL[i]->SetMarkerStyle(20+i);
      vecHcosThetaL[i]->SetMarkerColor(1+i);
      vecHcosThetaL[i]->SetLineColor(1+i);
      vecHcosThetaL[i]->SetLineWidth(2);
      vecHcosThetaL[i]->Draw("sames e1");
      vecHcosThetaL[i]->GetYaxis()->SetRangeUser(0.0,ordinateRange);
      myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin " << i;
      legThetaL->AddEntry(vecHcosThetaL[i],myString.str().c_str());
    }
  legThetaL->SetFillColor(0);
  legThetaL->SetBorderSize(0);
  legThetaL->Draw();

  cEff->cd(4);
  vecHphi[0]->SetMarkerStyle(20);
  vecHphi[0]->SetMarkerColor(1);
  vecHphi[0]->SetLineColor(1);
  vecHphi[0]->SetLineWidth(2);
  vecHphi[0]->Draw("e1");
  vecHphi[0]->GetYaxis()->SetRangeUser(0.0,ordinateRange);
  TLegend* legPhi = new TLegend(0.88, 0.65, 0.97, 0.89, "");
  myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin 0";
  legPhi->AddEntry(vecHphi[0],myString.str().c_str());
  for (unsigned int i = 1; i < q2Bins->size()-1; i++)
    {
      vecHphi[i]->SetMarkerStyle(20+i);
      vecHphi[i]->SetMarkerColor(1+i);
      vecHphi[i]->SetLineColor(1+i);
      vecHphi[i]->SetLineWidth(2);
      vecHphi[i]->Draw("sames e1");
      vecHphi[i]->GetYaxis()->SetRangeUser(0.0,ordinateRange);
      myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin " << i;
      legPhi->AddEntry(vecHphi[i],myString.str().c_str());
    }
  legPhi->SetFillColor(0);
  legPhi->SetBorderSize(0);
  legPhi->Draw();

  cEff->Update();


  // #############################
  // # Eff. Terms 2D projections #
  // #############################
  TCanvas* cNumCosThetaK = new TCanvas("cNumCosThetaK", "cNumCosThetaK", 10, 10, 1200, 800);
  cNumCosThetaK->Divide(2,3);

  TCanvas* cNumCosThetaL = new TCanvas("cNumCosThetaL", "cNumCosThetaL", 10, 10, 1200, 800);
  cNumCosThetaL->Divide(2,3);

  TCanvas* cNumPhi = new TCanvas("cNumPhi", "cNumPhi", 10, 10, 1200, 800);
  cNumPhi->Divide(2,3);

  vector<TH2D*> vecHq2ANDcosThetaK;
  vector<TH2D*> vecHq2ANDcosThetaL;
  vector<TH2D*> vecHq2ANDphi;

  vector<TLegend*> legHq2ANDcosThetaK;
  vector<TLegend*> legHq2ANDcosThetaL;
  vector<TLegend*> legHq2ANDphi;

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    {
      for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	{
	  myString.clear(); myString.str("");
	  myString << "vecHq2ANDcosThetaK_" << i << "_" << j;
	  vecHq2ANDcosThetaK.push_back(new TH2D(myString.str().c_str(), myString.str().c_str(), cosThetaLBins->size()-1, cosThetaLBins_, phiBins->size()-1, phiBins_));
	  vecHq2ANDcosThetaK.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
	  vecHq2ANDcosThetaK.back()->SetYTitle("#phi");
	  vecHq2ANDcosThetaK.back()->SetZTitle("Eff. Terms [#]");

	  legHq2ANDcosThetaK.push_back(new TLegend(0.88, 0.6, 0.97, 0.89, ""));
	}

      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	{
	  myString.clear(); myString.str("");
	  myString << "vecHq2ANDcosThetaL_" << i << "_" << k;
	  vecHq2ANDcosThetaL.push_back(new TH2D(myString.str().c_str(), myString.str().c_str(), cosThetaKBins->size()-1, cosThetaKBins_, phiBins->size()-1, phiBins_));
	  vecHq2ANDcosThetaL.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
	  vecHq2ANDcosThetaL.back()->SetYTitle("#phi");
	  vecHq2ANDcosThetaL.back()->SetZTitle("Eff. Terms [#]");

	  legHq2ANDcosThetaL.push_back(new TLegend(0.88, 0.6, 0.97, 0.89, ""));
	}

      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  myString.clear(); myString.str("");
	  myString << "vecHq2ANDphi_" << i << "_" << l;
	  vecHq2ANDphi.push_back(new TH2D(myString.str().c_str(), myString.str().c_str(), cosThetaLBins->size()-1, cosThetaLBins_, cosThetaKBins->size()-1, cosThetaKBins_));
	  vecHq2ANDphi.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
	  vecHq2ANDphi.back()->SetYTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
	  vecHq2ANDphi.back()->SetZTitle("Eff. Terms [#]");

	  legHq2ANDphi.push_back(new TLegend(0.88, 0.6, 0.97, 0.89, ""));
	}
    }


  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
    {
      for (unsigned int i = 0; i < q2Bins->size()-1; i++)
	{
	  // ###############################################################################################################
	  // # Fill histogram : numerator of the efficiency vs cos(theta_l) in bins of dimuon q^2 and cos(theta_K) and phi #
	  // ###############################################################################################################
	  for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	    for (unsigned int l = 0; l < phiBins->size()-1; l++)
	      {
		if      (NumOrDen2Plot == "N1") Eff = myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "N2") Eff = myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "D1") Eff = myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else                            Eff = myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

		vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->SetBinContent(k+1,l+1,Eff);
	      }
	  
	  cNumCosThetaK->cd(j+1)->SetLogz();
	  vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->SetMarkerStyle(1);
	  vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->SetMarkerColor(1+i);
	  vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->SetLineColor(1+i);
	  vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->SetLineWidth(2);
	  if (i == 0) vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->Draw("surf bb fb");
	  else        vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->Draw("sames surf bb fb");
	  vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->GetZaxis()->SetRangeUser(0.1,Zaxes);

	  myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin " << i;
	  legHq2ANDcosThetaK[j]->AddEntry(vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i],myString.str().c_str());
	}

      legHq2ANDcosThetaK[j]->SetFillColor(0);
      legHq2ANDcosThetaK[j]->SetBorderSize(0);
      legHq2ANDcosThetaK[j]->Draw();
    }
  

  for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
    {
      for (unsigned int i = 0; i < q2Bins->size()-1; i++)
	{
	  // #######################################################################################################
	  // # Fill histogram : numerator of the efficiency vs cos(theta_K) in bins of dimuon q^2 and cos(theta_l) #
	  // #######################################################################################################
	  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	    for (unsigned int l = 0; l < phiBins->size()-1; l++)
	      {
		if      (NumOrDen2Plot == "N1") Eff = myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "N2") Eff = myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "D1") Eff = myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else                            Eff = myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

		vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->SetBinContent(j+1,l+1,Eff);
	      }

	  cNumCosThetaL->cd(k+1)->SetLogz();
	  vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->SetMarkerStyle(1);
	  vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->SetMarkerColor(1+i);
	  vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->SetLineColor(1+i);
	  vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->SetLineWidth(2);
	  if (i == 0) vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->Draw("surf bb fb");
	  else vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->Draw("sames surf bb fb");
	  vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->GetZaxis()->SetRangeUser(0.1,Zaxes);

	  myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin " << i;
	  legHq2ANDcosThetaL[k]->AddEntry(vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i],myString.str().c_str());
	}

      legHq2ANDcosThetaL[k]->SetFillColor(0);
      legHq2ANDcosThetaL[k]->SetBorderSize(0);
      legHq2ANDcosThetaL[k]->Draw();
    }


  for (unsigned int l = 0; l < phiBins->size()-1; l++)
    {
      for (unsigned int i = 0; i < q2Bins->size()-1; i++)
	{
	  // ###############################################################################################################
	  // # Fill histogram : numerator of the efficiency vs phi in bins of dimuon q^2 and cos(theta_l) and cos(theta_k) #
	  // ###############################################################################################################
	  for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	      {
		if      (NumOrDen2Plot == "N1") Eff = myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "N2") Eff = myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "D1") Eff = myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else                            Eff = myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

		vecHq2ANDphi[l*(q2Bins->size()-1)+i]->SetBinContent(k+1,j+1,Eff);
	      }
	  
	  cNumPhi->cd(l+1)->SetLogz();
	  vecHq2ANDphi[l*(q2Bins->size()-1)+i]->SetMarkerStyle(1);
	  vecHq2ANDphi[l*(q2Bins->size()-1)+i]->SetMarkerColor(1+i);
	  vecHq2ANDphi[l*(q2Bins->size()-1)+i]->SetLineColor(1+i);
	  vecHq2ANDphi[l*(q2Bins->size()-1)+i]->SetLineWidth(2);
	  if (i == 0) vecHq2ANDphi[l*(q2Bins->size()-1)+i]->Draw("surf bb fb");
	  else        vecHq2ANDphi[l*(q2Bins->size()-1)+i]->Draw("sames surf bb fb");
	  vecHq2ANDphi[l*(q2Bins->size()-1)+i]->GetZaxis()->SetRangeUser(0.1,Zaxes);

	  myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin " << i;
	  legHq2ANDphi[l]->AddEntry(vecHq2ANDphi[l*(q2Bins->size()-1)+i],myString.str().c_str());
	}

      legHq2ANDphi[l]->SetFillColor(0);
      legHq2ANDphi[l]->SetBorderSize(0);
      legHq2ANDphi[l]->Draw();
    }


  cNumCosThetaL->Update();
  cNumCosThetaK->Update();
  cNumPhi->Update();


  // #############################
  // # Eff. Terms 1D projections #
  // #############################
  TCanvas* cNum = new TCanvas("cNum", "cNum", 10, 10, 1200, 800);
  cNum->Divide(2,2);

  vector<TH1D*> vecHq2ANDNumcosThetaK;
  vector<TH1D*> vecHq2ANDNumcosThetaL;
  vector<TH1D*> vecHq2ANDNumPhi;

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    {
      myString.clear(); myString.str("");
      myString << "vecHq2ANDNumcosThetaK_" << i;
      vecHq2ANDNumcosThetaK.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), cosThetaKBins->size()-1, cosThetaKBins_));
      vecHq2ANDNumcosThetaK.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
      vecHq2ANDNumcosThetaK.back()->SetYTitle("Eff. Terms [#]");

      myString.clear(); myString.str("");
      myString << "vecHq2ANDNumcosThetaL_" << i;
      vecHq2ANDNumcosThetaL.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), cosThetaLBins->size()-1, cosThetaLBins_));
      vecHq2ANDNumcosThetaL.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
      vecHq2ANDNumcosThetaL.back()->SetYTitle("Eff. Terms [#]");

      myString.clear(); myString.str("");
      myString << "vecHq2ANDNumPhi_" << i;
      vecHq2ANDNumPhi.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), phiBins->size()-1, phiBins_));
      vecHq2ANDNumPhi.back()->SetXTitle("#phi");
      vecHq2ANDNumPhi.back()->SetYTitle("Eff. Terms [#]");
    }


  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    {
      for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	{
	  Eff = 0.0;
	  
	  // ########################################################################
	  // # Fill histogram : efficiency vs cos(theta_K) for different dimuon q^2 #
	  // ########################################################################
	  for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	    for (unsigned int l = 0; l < phiBins->size()-1; l++)
	      {
		if      (NumOrDen2Plot == "N1") Eff += myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "N2") Eff += myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "D1") Eff += myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else                            Eff += myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      }

	  Eff = Eff / (cosThetaKBins->operator[](j+1) - cosThetaKBins->operator[](j));
	  vecHq2ANDNumcosThetaK[i]->SetBinContent(j+1,Eff);
	}


      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	{
	  Eff = 0.0;

	  // ########################################################################
	  // # Fill histogram : efficiency vs cos(theta_l) for different dimuon q^2 #
	  // ########################################################################
	  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	    for (unsigned int l = 0; l < phiBins->size()-1; l++)
	      {
		if      (NumOrDen2Plot == "N1") Eff += myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "N2") Eff += myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "D1") Eff += myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else                            Eff += myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      }

	  Eff = Eff / (cosThetaLBins->operator[](k+1) - cosThetaLBins->operator[](k));
	  vecHq2ANDNumcosThetaL[i]->SetBinContent(k+1,Eff);
	}


      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  Eff = 0.0;

	  // ###############################################################
	  // # Fill histogram : efficiency vs phi for different dimuon q^2 #
	  // ###############################################################
	  for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	      {
		if      (NumOrDen2Plot == "N1") Eff += myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "N2") Eff += myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "D1") Eff += myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else                            Eff += myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      }

	  Eff = Eff / (phiBins->operator[](l+1) - phiBins->operator[](l));
	  vecHq2ANDNumPhi[i]->SetBinContent(l+1,Eff);
	}
    }


  cNum->cd(2);
  vecHq2ANDNumcosThetaK[0]->SetMarkerStyle(20);
  vecHq2ANDNumcosThetaK[0]->SetMarkerColor(1);
  vecHq2ANDNumcosThetaK[0]->SetLineColor(1);
  vecHq2ANDNumcosThetaK[0]->SetLineWidth(2);
  vecHq2ANDNumcosThetaK[0]->Draw("p");
  vecHq2ANDNumcosThetaK[0]->GetYaxis()->SetRangeUser(0.1,Yaxes);
  TLegend* legNumThetaK = new TLegend(0.88, 0.65, 0.97, 0.89, "");
  myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin 0";
  legNumThetaK->AddEntry(vecHq2ANDNumcosThetaK[0],myString.str().c_str());
  for (unsigned int i = 1; i < q2Bins->size()-1; i++)
    {
      vecHq2ANDNumcosThetaK[i]->SetMarkerStyle(20+i);
      vecHq2ANDNumcosThetaK[i]->SetMarkerColor(1+i);
      vecHq2ANDNumcosThetaK[i]->SetLineColor(1+i);
      vecHq2ANDNumcosThetaK[i]->SetLineWidth(2);
      vecHq2ANDNumcosThetaK[i]->Draw("sames p");
      vecHq2ANDNumcosThetaK[i]->GetYaxis()->SetRangeUser(0.1,Yaxes);
      myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin " << i;
      legNumThetaK->AddEntry(vecHq2ANDNumcosThetaK[i],myString.str().c_str());
    }
  legNumThetaK->SetFillColor(0);
  legNumThetaK->SetBorderSize(0);
  legNumThetaK->Draw();

  cNum->cd(3);
  vecHq2ANDNumcosThetaL[0]->SetMarkerStyle(20);
  vecHq2ANDNumcosThetaL[0]->SetMarkerColor(1);
  vecHq2ANDNumcosThetaL[0]->SetLineColor(1);
  vecHq2ANDNumcosThetaL[0]->SetLineWidth(2);
  vecHq2ANDNumcosThetaL[0]->Draw("p");
  vecHq2ANDNumcosThetaL[0]->GetYaxis()->SetRangeUser(0.1,Yaxes);
  TLegend* legNumThetaL = new TLegend(0.88, 0.65, 0.97, 0.89, "");
  myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin 0";
  legNumThetaL->AddEntry(vecHq2ANDNumcosThetaL[0],myString.str().c_str());
  for (unsigned int i = 1; i < q2Bins->size()-1; i++)
    {
      vecHq2ANDNumcosThetaL[i]->SetMarkerStyle(20+i);
      vecHq2ANDNumcosThetaL[i]->SetMarkerColor(1+i);
      vecHq2ANDNumcosThetaL[i]->SetLineColor(1+i);
      vecHq2ANDNumcosThetaL[i]->SetLineWidth(2);
      vecHq2ANDNumcosThetaL[i]->Draw("sames p");
      vecHq2ANDNumcosThetaL[i]->GetYaxis()->SetRangeUser(0.1,Yaxes);
      myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin " << i;
      legNumThetaL->AddEntry(vecHq2ANDNumcosThetaL[i],myString.str().c_str());
    }
  legNumThetaL->SetFillColor(0);
  legNumThetaL->SetBorderSize(0);
  legNumThetaL->Draw();

  cNum->cd(4);
  vecHq2ANDNumPhi[0]->SetMarkerStyle(20);
  vecHq2ANDNumPhi[0]->SetMarkerColor(1);
  vecHq2ANDNumPhi[0]->SetLineColor(1);
  vecHq2ANDNumPhi[0]->SetLineWidth(2);
  vecHq2ANDNumPhi[0]->Draw("p");
  vecHq2ANDNumPhi[0]->GetYaxis()->SetRangeUser(0.1,Yaxes);
  TLegend* legNumPhi = new TLegend(0.88, 0.65, 0.97, 0.89, "");
  myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin 0";
  legNumPhi->AddEntry(vecHq2ANDNumPhi[0],myString.str().c_str());
  for (unsigned int i = 1; i < q2Bins->size()-1; i++)
    {
      vecHq2ANDNumPhi[i]->SetMarkerStyle(20+i);
      vecHq2ANDNumPhi[i]->SetMarkerColor(1+i);
      vecHq2ANDNumPhi[i]->SetLineColor(1+i);
      vecHq2ANDNumPhi[i]->SetLineWidth(2);
      vecHq2ANDNumPhi[i]->Draw("sames p");
      vecHq2ANDNumPhi[i]->GetYaxis()->SetRangeUser(0.1,Yaxes);
      myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin " << i;
      legNumPhi->AddEntry(vecHq2ANDNumPhi[i],myString.str().c_str());
    }
  legNumPhi->SetFillColor(0);
  legNumPhi->SetBorderSize(0);
  legNumPhi->Draw();

  cNum->Update();


  Utility->IntegrateEffAll(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,myEff,&totalEffAll,&EffErr);
  Utility->IntegrateEffButPsi(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,myEff,&totalEffSignal,&EffErr);
  Utility->IntegrateEffInJPsi(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,myEff,&totalEffJPsi,&EffErr);
  Utility->IntegrateEffInPsiP(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,myEff,&totalEffPsiP,&EffErr);
  cout << "\nTotal efficiency: " << totalEffAll <<  "\tTotal signal efficiency: " << totalEffSignal;
  cout << "\tTotal J/psi efficiency: " << totalEffJPsi << "\tTotal psi(2S) efficiency: " << totalEffPsiP << endl;
}


void Read3DEfficiencies (bool isSingleEff, vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
			 string fileNameInput, bool isAnalyEff, Utils::effStruct* myEff, bool CheckEffatRead, bool savePlot, int specBin)
{
  // ###################
  // # Local variables #
  // ###################
  const unsigned int MAXVAL = (isSingleEff == true ? 1 : NFILES+1);
  double totalEffAll;
  double totalEffSignal;
  double totalEffJPsi;
  double totalEffPsiP;
  double Eff, EffErr;
  double* q2Bins_        = Utility->MakeBinning(q2Bins);
  double* cosThetaKBins_ = Utility->MakeBinning(cosThetaKBins);
  double* cosThetaLBins_ = Utility->MakeBinning(cosThetaLBins);
  double* phiBins_       = Utility->MakeBinning(phiBins);
  vector<TF3*> effFuncs3D;
  string tmpString;
  stringstream myString;
  stringstream effFuncIDs;
  // ###################
  double Xaxes = 0.2e-3;
  // ###################


  TCanvas* cEff0 = new TCanvas("cEff0", "cEff0", 10, 10, 700, 500);
  TCanvas* cEff1 = new TCanvas("cEff1", "cEff1", 20, 20, 700, 500);
  TCanvas* cEff2 = new TCanvas("cEff2", "cEff2", 30, 30, 700, 500);
  TCanvas* cEff3 = new TCanvas("cEff3", "cEff3", 40, 40, 700, 500);

  tmpString = INPUTGenEff;
  tmpString.erase(tmpString.find(".txt"),4);

  vector<TH1D*> Hq2;
  for (unsigned int i = 0; i < MAXVAL; i++)
    {
      myString.clear(); myString.str("");
      myString << "Hq2_" << i;
      Hq2.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), q2Bins->size()-1, q2Bins_));
      Hq2.back()->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
      Hq2.back()->SetYTitle("Efficiency");
    }

  vector<TH1D*> HcosThetaK;
  for (unsigned int i = 0; i < MAXVAL; i++)
    {
      myString.clear(); myString.str("");
      myString << "HcosThetaK_" << i;
      HcosThetaK.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), cosThetaKBins->size()-1, cosThetaKBins_));
      HcosThetaK.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
      HcosThetaK.back()->SetYTitle("Efficiency");
    }

  vector<TH1D*> HcosThetaL;
  for (unsigned int i = 0; i < MAXVAL; i++)
    {
      myString.clear(); myString.str("");
      myString << "HcosThetaL_" << i;
      HcosThetaL.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), cosThetaLBins->size()-1, cosThetaLBins_));
      HcosThetaL.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
      HcosThetaL.back()->SetYTitle("Efficiency");
    }

  vector<TH1D*> Hphi;
  for (unsigned int i = 0; i < MAXVAL; i++)
    {
      myString.clear(); myString.str("");
      myString << "Hphi_" << i;
      Hphi.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), phiBins->size()-1, phiBins_));
      Hphi.back()->SetXTitle("#phi");
      Hphi.back()->SetYTitle("Efficiency");
    }


  // ###################################
  // # Initialize efficiency structure #
  // ###################################
  Utility->GenEfficiency(myEff,q2Bins,cosThetaKBins,cosThetaLBins,phiBins);
  Utils::effValue myEffVal;
  Utility->ResetEffValue(&myEffVal,0.0);
  Utility->InitEfficiency(myEffVal,*myEff,q2Bins,cosThetaKBins,cosThetaLBins,phiBins);


  for (unsigned int itF = 0; itF < MAXVAL; itF++)
    {
      if ((isSingleEff == true) || (itF == MAXVAL-1))
	{
	  if (isAnalyEff == false) Utility->ReadEfficiency(fileNameInput.c_str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,myEff);
	  else                     Utility->ReadAnalyticalEff(fileNameInput.c_str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,&effFuncs3D,"effFuncs3D",0);
	}
      else
	{
	  myString.clear(); myString.str("");
	  myString << tmpString << "_" << itF << ".txt";
	  effFuncIDs << "effFuncs3D_" << itF;
 	  if (isAnalyEff == false) Utility->ReadEfficiency(myString.str().c_str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,myEff);
	  else                     Utility->ReadAnalyticalEff(myString.str().c_str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,&effFuncs3D,effFuncIDs.str().c_str(),0);
	}


      // ################################################
      // # Check if analytical efficiency goes negative #
      // ################################################
      if ((isAnalyEff == true) && (CheckEffatRead == true))
      	for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
      	  if (Utility->EffMinValue3D(cosThetaKBins,cosThetaLBins,phiBins,effFuncs3D[i]) < 0.0)
	    cout << "@@@ Negative efficiency function #" << itF << " for q2 bin #" << i << " ! @@@" << endl;

      
      for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
	{
	  // #############################################
	  // # Fill histogram : efficiency vs dimuon q^2 #
	  // #############################################
	  if (isAnalyEff == false) Utility->IntegrateEffButq2(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,i,*myEff,&Eff,&EffErr);
	  else
	    {
	      Eff = 0.0;
	      Eff = Eff + effFuncs3D[i]->Integral(cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),
						  cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1),
						  phiBins->operator[](0),phiBins->operator[](phiBins->size()-1)) /
		((cosThetaKBins->operator[](cosThetaKBins->size()-1)-cosThetaKBins->operator[](0)) *
		 (cosThetaLBins->operator[](cosThetaLBins->size()-1)-cosThetaLBins->operator[](0)) *
		 (phiBins->operator[](phiBins->size()-1)-phiBins->operator[](0)));
	      EffErr = 0.0;
	    }
	  Hq2[itF]->SetBinContent(i+1,Eff);
	  Hq2[itF]->SetBinError(i+1,EffErr);
	}
      
      for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      	{
      	  // ###############################################
      	  // # Fill histogram : efficiency vs cos(theta_K) #
      	  // ###############################################
      	  if (isAnalyEff == false) Utility->IntegrateEffButCosThetaK(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,j,*myEff,&Eff,&EffErr);
      	  else
      	    {
      	      Eff = 0.0;
	      for (unsigned int indx = (specBin == -1 ? 0 : specBin); indx < (specBin == -1 ? q2Bins->size()-1 : specBin+1); indx++)
      		Eff = Eff + effFuncs3D[indx]->Integral(cosThetaKBins->operator[](j),cosThetaKBins->operator[](j+1),cosThetaLBins->operator[](0),
						       cosThetaLBins->operator[](cosThetaLBins->size()-1),
						       phiBins->operator[](0),phiBins->operator[](phiBins->size()-1)) /
      		  ((cosThetaKBins->operator[](j+1)-cosThetaKBins->operator[](j)) *
		   (cosThetaLBins->operator[](cosThetaLBins->size()-1)-cosThetaLBins->operator[](0)) *
		   (phiBins->operator[](phiBins->size()-1)-phiBins->operator[](0)));
      	      EffErr = 0.0;
      	    }
      	  HcosThetaK[itF]->SetBinContent(j+1,Eff);
      	  HcosThetaK[itF]->SetBinError(j+1,EffErr);
      	}

      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
      	{
      	  // ###############################################
      	  // # Fill histogram : efficiency vs cos(theta_l) #
      	  // ###############################################
      	  if (isAnalyEff == false) Utility->IntegrateEffButCosThetaL(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,k,*myEff,&Eff,&EffErr);
      	  else
      	    {
      	      Eff = 0.0;
	      for (unsigned int indx = (specBin == -1 ? 0 : specBin); indx < (specBin == -1 ? q2Bins->size()-1 : specBin+1); indx++)
      		Eff = Eff + effFuncs3D[indx]->Integral(cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),
						       cosThetaLBins->operator[](k),cosThetaLBins->operator[](k+1),
						       phiBins->operator[](0),phiBins->operator[](phiBins->size()-1)) /
      		  ((cosThetaKBins->operator[](cosThetaKBins->size()-1)-cosThetaKBins->operator[](0)) *
		   (cosThetaLBins->operator[](k+1)-cosThetaLBins->operator[](k)) *
		   (phiBins->operator[](phiBins->size()-1)-phiBins->operator[](0)));
      	      EffErr = 0.0;
      	    }
      	  HcosThetaL[itF]->SetBinContent(k+1,Eff);
      	  HcosThetaL[itF]->SetBinError(k+1,EffErr);
      	}

      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  // ######################################
	  // # Fill histogram : efficiency vs phi #
	  // ######################################
	  if (isAnalyEff == false) Utility->IntegrateEffButPhi(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,l,*myEff,&Eff,&EffErr);
	  else
	    {
	      Eff = 0.0;
	      for (unsigned int indx = (specBin == -1 ? 0 : specBin); indx < (specBin == -1 ? q2Bins->size()-1 : specBin+1); indx++)
		Eff = Eff + effFuncs3D[indx]->Integral(cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),
						       cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1),
						       phiBins->operator[](l),phiBins->operator[](l+1)) /
		  ((cosThetaKBins->operator[](cosThetaKBins->size()-1)-cosThetaKBins->operator[](0)) *
		   (cosThetaLBins->operator[](cosThetaLBins->size()-1)-cosThetaLBins->operator[](0)) *
		   (phiBins->operator[](l+1)-phiBins->operator[](l)));
	      EffErr = 0.0;
	    }
	  Hphi[itF]->SetBinContent(l+1,Eff);
	  Hphi[itF]->SetBinError(l+1,EffErr);
	}
      
      
      if ((isSingleEff == true) || (itF == MAXVAL-1))
	{
	  cEff0->cd();
	  Hq2[itF]->SetMarkerStyle(20);
	  Hq2[itF]->GetYaxis()->SetRangeUser(0,ordinateRange);
	  if (isSingleEff == true) Hq2[itF]->Draw("pe1");
	  else Hq2[itF]->Draw("same pe1");
	  
	  cEff1->cd();
	  HcosThetaK[itF]->SetMarkerStyle(20);
	  HcosThetaK[itF]->GetYaxis()->SetRangeUser(0,ordinateRange);
	  if (isSingleEff == true) HcosThetaK[itF]->Draw("pe1");
	  else HcosThetaK[itF]->Draw("same pe1");

	  cEff2->cd();
	  HcosThetaL[itF]->SetMarkerStyle(20);
	  HcosThetaL[itF]->GetYaxis()->SetRangeUser(0,ordinateRange);
	  if (isSingleEff == true) HcosThetaL[itF]->Draw("pe1");
	  else HcosThetaL[itF]->Draw("same pe1");

	  cEff3->cd();
	  Hphi[itF]->SetMarkerStyle(20);
	  Hphi[itF]->GetYaxis()->SetRangeUser(0,ordinateRange);
	  if (isSingleEff == true) Hphi[itF]->Draw("pe1");
	  else Hphi[itF]->Draw("same pe1");
	}
      else
	{
	  cEff0->cd();
	  Hq2[itF]->SetLineColor(kBlue);
	  Hq2[itF]->GetYaxis()->SetRangeUser(0,ordinateRange);
	  if (itF == 0) Hq2[itF]->Draw();
	  else Hq2[itF]->Draw("same");

	  cEff1->cd();
	  HcosThetaK[itF]->SetLineColor(kBlue);
	  HcosThetaK[itF]->GetYaxis()->SetRangeUser(0,ordinateRange);
	  if (itF == 0) HcosThetaK[itF]->Draw();
	  else HcosThetaK[itF]->Draw("same");

	  cEff2->cd();
	  HcosThetaL[itF]->SetLineColor(kBlue);
	  HcosThetaL[itF]->GetYaxis()->SetRangeUser(0,ordinateRange);
	  if (itF == 0) HcosThetaL[itF]->Draw();
	  else HcosThetaL[itF]->Draw("same");

	  cEff3->cd();
	  Hphi[itF]->SetLineColor(kBlue);
	  Hphi[itF]->GetYaxis()->SetRangeUser(0,ordinateRange);
	  if (itF == 0) Hphi[itF]->Draw();
	  else Hphi[itF]->Draw("same");
	}
      
      effFuncs3D.clear();
    }
 
  cEff0->Update();
  cEff1->Update();
  cEff2->Update();
  cEff3->Update();

  if (savePlot == true)
    {
      if (isAnalyEff == false) cEff0->Print("Histo1.pdf");
      cEff1->Print("Histo2.pdf");
      cEff2->Print("Histo3.pdf");
      if (isAnalyEff == false) cEff3->Print("Histo4.pdf");
    }


  TCanvas* cStatK = new TCanvas("cStatK", "cStatK", 10, 10, 1200, 800);
  cStatK->Divide(2,3);

  vector<TH1D*> HstatK;
  for (unsigned int i = 0; i < cosThetaKBins->size()-1; i++)
    {
      myString.clear(); myString.str("");
      myString << "HstatK_" << i;
      HstatK.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), 50, -Xaxes, Xaxes));
      HstatK.back()->SetXTitle("Diff. to reference");
      HstatK.back()->SetYTitle("Entries [#]");

      cStatK->cd(i+1);
      HstatK.back()->Draw();
    }


  TCanvas* cStatL = new TCanvas("cStatL", "cStatL", 10, 10, 1200, 800);
  cStatL->Divide(2,3);
  
  vector<TH1D*> HstatL;
  for (unsigned int i = 0; i < cosThetaLBins->size()-1; i++)
    {
      myString.clear(); myString.str("");
      myString << "HstatL_" << i;
      HstatL.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), 50, -Xaxes, Xaxes));
      HstatL.back()->SetXTitle("Diff. to reference");
      HstatL.back()->SetYTitle("Entries [#]");

      cStatL->cd(i+1);
      HstatL.back()->Draw();
    }


  TCanvas* cStatPhi = new TCanvas("cStatPhi", "cStatPhi", 10, 10, 1200, 800);
  cStatPhi->Divide(2,3);
  
  vector<TH1D*> HstatPhi;
  for (unsigned int i = 0; i < phiBins->size()-1; i++)
    {
      myString.clear(); myString.str("");
      myString << "HstatPhi_" << i;
      HstatPhi.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), 50, -Xaxes, Xaxes));
      HstatPhi.back()->SetXTitle("Diff. to reference");
      HstatPhi.back()->SetYTitle("Entries [#]");

      cStatPhi->cd(i+1);
      HstatPhi.back()->Draw();
    }


  if (isSingleEff == false)
    {
      for (unsigned int i = 0; i < cosThetaKBins->size()-1; i++)
  	for (unsigned int itF = 0; itF < MAXVAL-1; itF++)
  	  HstatK[i]->Fill(HcosThetaK[itF]->GetBinContent(i+1) - HcosThetaK[MAXVAL-1]->GetBinContent(i+1));

      for (unsigned int i = 0; i < HstatK.size(); i++)
  	{
  	  cStatK->cd(i+1);
  	  HstatK[i]->SetFillColor(kAzure+6);
  	  HstatK[i]->Draw();
  	}


      for (unsigned int i = 0; i < cosThetaLBins->size()-1; i++)
  	for (unsigned int itF = 0; itF < MAXVAL-1; itF++)
  	  HstatL[i]->Fill(HcosThetaL[itF]->GetBinContent(i+1) - HcosThetaL[MAXVAL-1]->GetBinContent(i+1));

      for (unsigned int i = 0; i < HstatL.size(); i++)
  	{
  	  cStatL->cd(i+1);
  	  HstatL[i]->SetFillColor(kAzure+6);
  	  HstatL[i]->Draw();
  	}
      

      for (unsigned int i = 0; i < phiBins->size()-1; i++)
  	for (unsigned int itF = 0; itF < MAXVAL-1; itF++)
  	  HstatPhi[i]->Fill(Hphi[itF]->GetBinContent(i+1) - Hphi[MAXVAL-1]->GetBinContent(i+1));

      for (unsigned int i = 0; i < HstatPhi.size(); i++)
  	{
  	  cStatL->cd(i+1);
  	  HstatPhi[i]->SetFillColor(kAzure+6);
  	  HstatPhi[i]->Draw();
  	}

      cStatK->Update();
      cStatL->Update();
      cStatPhi->Update();
    }


  Utility->IntegrateEffAll(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,*myEff,&totalEffAll,&EffErr);
  Utility->IntegrateEffButPsi(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,*myEff,&totalEffSignal,&EffErr);
  Utility->IntegrateEffInJPsi(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,*myEff,&totalEffJPsi,&EffErr);
  Utility->IntegrateEffInPsiP(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,*myEff,&totalEffPsiP,&EffErr);
  cout << "\nTotal efficiency: " << totalEffAll <<  "\tTotal signal efficiency: " << totalEffSignal;
  cout << "\tTotal J/psi efficiency: " << totalEffJPsi << "\tTotal psi(2S) efficiency: " << totalEffPsiP << endl;


  if ((isSingleEff == true) && (isAnalyEff == false)) MakeHistogramsAllBins(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,*myEff);
}


void Fit1DEfficiencies (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
			unsigned int SignalType, Utils::effStruct myEff, string who, unsigned int q2BinIndx, string fileNameOut)
{
  // ###################
  // # Local variables #
  // ###################
  TFitResultPtr fitResults;
  double Eff, EffErr;
  double* cosThetaKBins_ = Utility->MakeBinning(cosThetaKBins);
  double* cosThetaLBins_ = Utility->MakeBinning(cosThetaLBins);
  double* phiBins_       = Utility->MakeBinning(phiBins);
  vector<TH1D*> histFit;
  TF1* effFunc1D;
  stringstream myString;
  string tmpString;
  double coeff;
  ofstream fileOutput;
  ifstream fileInput;
  bool nullParameter;
  // ###################


  TCanvas* cEff = new TCanvas("cEff", "cEff", 10, 10, 1200, 800);
  cEff->Divide(3,2);

  if (who == "thetaL")
    {
      fileOutput.open(fileNameOut.replace(fileNameOut.find(".txt"),4,"L.txt").c_str(),ofstream::app);
      if (fileOutput.good() == false)
	{
	  cout << "[ComputeEfficiency::Fit1DEfficiencies]\tError opening file : " << fileNameOut.c_str() << endl;
	  exit (EXIT_FAILURE);
	}


      // ##############################
      // # Read analytical efficiency #
      // ##############################
      effFunc1D = new TF1("effFunc1D",Utility->TellMeEffFuncThetaL().c_str(),-1.0,1.0);


      for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	{
	  Utility->InitEffFuncThetaL(effFunc1D,q2BinIndx);

	  myString.clear(); myString.str("");
	  myString << "histFit_" << j;
	  histFit.push_back(new TH1D(myString.str().c_str(),myString.str().c_str(),cosThetaLBins->size()-1,cosThetaLBins_));

	  for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	    {
	      Utility->IntegrateEffPhi(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2Bins->operator[](q2BinIndx),cosThetaKBins->operator[](j),cosThetaLBins->operator[](k),myEff,&Eff,&EffErr,false);
	      
	      histFit.back()->SetBinContent(k+1,Eff);
	      histFit.back()->SetBinError(k+1,EffErr);
	    }


	  // ######################################################################
	  // # Add constraint where it is necessary to bound the function at zero #
	  // ######################################################################
	  Utility->AddConstraintThetaL(&histFit.back(),q2BinIndx,j,SignalType,q2BinIndx);


	  // #########################################################################
	  // # Perform the fit of the analytical efficiency to the binned efficiency #
	  // #########################################################################
	  histFit.back()->SetMarkerStyle(20);
	  histFit.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
	  histFit.back()->SetYTitle("Efficiency");
	  histFit.back()->GetYaxis()->SetRangeUser(-ordinateRange,ordinateRange);

	  cEff->cd(j+1);
	  fitResults = histFit.back()->Fit("effFunc1D","V");
	  histFit.back()->Draw("pe1");
	  if (fitResults != 0) exit (EXIT_FAILURE);


	  // ################
	  // # Save results #
	  // ################
	  fileOutput << (q2Bins->operator[](q2BinIndx)+q2Bins->operator[](q2BinIndx+1))/2.;
	  for (int i = 0; i < effFunc1D->GetNpar(); i++) fileOutput << "   " << effFunc1D->GetParameter(i) << "   " << effFunc1D->GetParError(i);
	  fileOutput << endl;

	  cout << "\n@@@ Fit for for bin n." << j << " @@@" << endl;
	  cout << "@@@ chi2/DoF = " << effFunc1D->GetChisquare() / effFunc1D->GetNDF() << " (" << effFunc1D->GetChisquare() << "/" << effFunc1D->GetNDF() << ")";
	  cout << "\tCL : " << TMath::Prob(effFunc1D->GetChisquare(),effFunc1D->GetNDF()) << " @@@" << endl;
	  cout << "@@@ Value at " << cosThetaLBins->operator[](0) << " : " << effFunc1D->Eval(cosThetaLBins->operator[](0)) << " @@@" << endl;
	  cout << "@@@ Value at " << cosThetaLBins->operator[](cosThetaLBins->size()-1) << " : " << effFunc1D->Eval(cosThetaLBins->operator[](cosThetaLBins->size()-1)) << " @@@\n" << endl;

	  if (Utility->EffMinValue1D(cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1),effFunc1D) < 0.0) { cout << "NEGATIVE EFFICIENCY !" << endl; exit (EXIT_FAILURE); }
	}


      fileOutput.close();
    }
  else if (who == "thetaK")
    {
      fileOutput.open(fileNameOut.replace(fileNameOut.find(".txt"),4,"K.txt").c_str(),ofstream::app);
      if (fileOutput.good() == false)
	{
	  cout << "[ComputeEfficiency::Fit1DEfficiencies]\tError opening file : " << fileNameOut.c_str() << endl;
	  exit (EXIT_FAILURE);
	}


      // ##############################
      // # Read analytical efficiency #
      // ##############################
      effFunc1D = new TF1("effFunc1D",Utility->TellMeEffFuncThetaK().c_str(),-1.0,1.0);


      for (unsigned int k = 0; k < Utility->NcoeffThetaL; k++)
	{
	  Utility->InitEffFuncThetaK(effFunc1D,q2BinIndx);

	  myString.clear(); myString.str("");
	  myString << "histFit_" << k;
	  histFit.push_back(new TH1D(myString.str().c_str(),myString.str().c_str(),cosThetaKBins->size()-1,cosThetaKBins_));


	  // ###############################
	  // # Read coefficients from file #
	  // ###############################
	  fileInput.open(INPUT_THETAL,ofstream::in);
	  if (fileInput.good() == false)
	    {
	      cout << "[ComputeEfficiency::Fit1DEfficiencies]\tError opening file : " << INPUT_THETAL << endl;
	      exit (EXIT_FAILURE);
	    }
	  for (unsigned int j = 0; j < (cosThetaKBins->size()-1)*q2BinIndx; j++) getline(fileInput,tmpString);
	  nullParameter = true;
	  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	    {
	      getline(fileInput,tmpString);
	      stringstream rawStringK(tmpString);
	      rawStringK >> coeff; // Discard q2 bin value
	      for (unsigned int i = 0; i < k*2; i++) rawStringK >> coeff;
	      
	      rawStringK >> coeff;
	      histFit.back()->SetBinContent(j+1,coeff);
	      rawStringK >> coeff;
	      histFit.back()->SetBinError(j+1,coeff);

	      cout << "Cos(theta_K) bin " << j << " --> reading coef. " << k << " for var. theta_l: " << histFit.back()->GetBinContent(j+1) << "+/-" <<  histFit.back()->GetBinError(j+1) << endl;

	      if (coeff == 0.0) nullParameter = nullParameter*true; // All parameter errors must be zero in order to ignore the parameter
	      else              nullParameter = false;
	    }


	  // #########################################################################
	  // # Perform the fit of the analytical efficiency to the binned efficiency #
	  // #########################################################################
	  histFit.back()->SetMarkerStyle(20);
	  histFit.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
	  histFit.back()->SetYTitle("Efficiency");
	  histFit.back()->GetYaxis()->SetRangeUser(-ordinateRange,ordinateRange);

	  cEff->cd(k+1);
	  fitResults = histFit.back()->Fit("effFunc1D","V");
	  histFit.back()->Draw("pe1");
	  if (fitResults != 0) exit (EXIT_FAILURE);


	  // ################
	  // # Save results #
	  // ################
	  fileOutput << (q2Bins->operator[](q2BinIndx)+q2Bins->operator[](q2BinIndx+1))/2.;
	  for (int i = 0; i < effFunc1D->GetNpar(); i++) fileOutput << "   " << effFunc1D->GetParameter(i) << "   " << (nullParameter == false ? effFunc1D->GetParError(i) : 0.0);
	  fileOutput << endl;

	  cout << "\n@@@ Fit for for bin n." << k << " @@@" << endl;
	  cout << "@@@ chi2/DoF = " << effFunc1D->GetChisquare() / effFunc1D->GetNDF() << " (" << effFunc1D->GetChisquare() << "/" << effFunc1D->GetNDF() << ")";
	  cout << "\tCL : " << TMath::Prob(effFunc1D->GetChisquare(),effFunc1D->GetNDF()) << " @@@" << endl;
	  cout << "@@@ Value at " << cosThetaKBins->operator[](0) << " : " << effFunc1D->Eval(cosThetaKBins->operator[](0)) << " @@@" << endl;
	  cout << "@@@ Value at " << cosThetaKBins->operator[](cosThetaKBins->size()-1) << " : " << effFunc1D->Eval(cosThetaKBins->operator[](cosThetaKBins->size()-1)) << " @@@\n" << endl;

	  fileInput.close();
	}


      fileOutput.close();
    }
  else if (who == "phi")
    {
      fileOutput.open(fileNameOut.replace(fileNameOut.find("Theta.txt"),9,"Phi.txt").c_str(),ofstream::app);
      if (fileOutput.good() == false)
	{
	  cout << "[ComputeEfficiency::Fit1DEfficiencies]\tError opening file : " << fileNameOut.c_str() << endl;
	  exit (EXIT_FAILURE);
	}

      
      // ##############################
      // # Read analytical efficiency #
      // ##############################
      effFunc1D = new TF1("effFunc1D",Utility->TellMeEffFuncPhi().c_str(),-Utility->PI,Utility->PI);
      Utility->InitEffFuncPhi(effFunc1D,q2BinIndx);


      myString.clear(); myString.str("");
      myString << "histFit";
      histFit.push_back(new TH1D(myString.str().c_str(),myString.str().c_str(),phiBins->size()-1,phiBins_));

      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  Utility->IntegrateEffCosThetaKCosThetaL(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2BinIndx,l,myEff,&Eff,&EffErr);

	  histFit.back()->SetBinContent(l+1,Eff);
	  histFit.back()->SetBinError(l+1,EffErr);
	}


      // ##############################
      // # Read analytical efficiency #
      // ##############################
      histFit.back()->SetMarkerStyle(20);
      histFit.back()->SetXTitle("#phi");
      histFit.back()->SetYTitle("Efficiency");
      histFit.back()->GetYaxis()->SetRangeUser(-ordinateRange,ordinateRange);

      cEff->cd(1);
      fitResults = histFit.back()->Fit("effFunc1D","V");
      histFit.back()->Draw("pe1");
      if (fitResults != 0) exit (EXIT_FAILURE);
	  

      // ################
      // # Save results #
      // ################
      fileOutput << (q2Bins->operator[](q2BinIndx)+q2Bins->operator[](q2BinIndx+1))/2.;
      for (int i = 0; i < effFunc1D->GetNpar(); i++) fileOutput << "   " << effFunc1D->GetParameter(i) << "   " << effFunc1D->GetParError(i);
      fileOutput << endl;

      cout << "@@@ Value at " << phiBins->operator[](0) << " : " << effFunc1D->Eval(phiBins->operator[](0)) << " @@@" << endl;
      cout << "@@@ Value at " << phiBins->operator[](phiBins->size()-1) << " : " << effFunc1D->Eval(phiBins->operator[](phiBins->size()-1)) << " @@@\n" << endl;

      if (Utility->EffMinValue1D(phiBins->operator[](0),phiBins->operator[](phiBins->size()-1),effFunc1D) < 0.0) { cout << "NEGATIVE EFFICIENCY !" << endl; exit (EXIT_FAILURE); }


      fileOutput.close();
    }


  cEff->Update();
  histFit.clear();
}


void Fit2DEfficiencies (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
			unsigned int SignalType, Utils::effStruct myEff, unsigned int q2BinIndx, string fileNameOut, bool savePlot)
{
  // ###################
  // # Local variables #
  // ###################
  TFitResultPtr fitResults;
  stringstream myString;
  vector<TF2*> effFuncs2D;
  // ###################


  TCanvas* cTestGlobalFit = new TCanvas("cTestGlobalFit", "cTestGlobalFit", 10, 10, 700, 500);


  // ##########################
  // # Read binned efficiency #
  // ##########################
  myString.clear(); myString.str("");
  myString << "Histo_q2Bin_" << q2BinIndx;
  TH2D* hisFunc2D = Utility->Get2DEffHitoq2Bin(myString.str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2BinIndx,myEff);


  // ##############################
  // # Read analytical efficiency #
  // ##############################
  Utility->ReadAnalyticalEff(INPUT_THETAL_THETAK,q2Bins,cosThetaKBins,cosThetaLBins,&effFuncs2D,"effFuncs2D",0);
  effFuncs2D[q2BinIndx]->SetRange(cosThetaKBins->operator[](0),
				  cosThetaLBins->operator[](0),
				  cosThetaKBins->operator[](cosThetaKBins->size()-1),
				  cosThetaLBins->operator[](cosThetaLBins->size()-1));

  for (int i = 0; i < effFuncs2D[q2BinIndx]->GetNpar(); i++)
    if (effFuncs2D[q2BinIndx]->GetParError(i) == 0.0) effFuncs2D[q2BinIndx]->FixParameter(i,effFuncs2D[q2BinIndx]->GetParameter(i));


  // ############################################################################################
  // # Add constraint along Y (= cosThetaL) where it is necessary to bound the function at zero #
  // ############################################################################################
  Utility->AddConstraintThetaKThetaL(&hisFunc2D,cosThetaKBins,q2BinIndx,SignalType,q2BinIndx);


  // #########################################################################
  // # Perform the fit of the analytical efficiency to the binned efficiency #
  // #########################################################################
  fitResults = hisFunc2D->Fit(effFuncs2D[q2BinIndx]->GetName(),"VMS");
  cout << "[ComputeEfficiency::Fit2DEfficiencies]\tFit status: " << fitResults << endl;

  if (fitResults != 0)
    {
      fitResults = hisFunc2D->Fit(effFuncs2D[q2BinIndx]->GetName(),"VS");
      cout << "[ComputeEfficiency::Fit2DEfficiencies]\tFit status: " << fitResults << endl;
    }

  TMatrixTSym<double> covMatrix(fitResults->GetCovarianceMatrix());

  if (fitResults != 0) exit (EXIT_FAILURE);
  else if (covMatrix.Determinant() <= 0.0)
    {
      cout << "[ComputeEfficiency::Fit2DEfficiencies]\tFit status: covariance matrix with negative determinant" << endl;
      exit (EXIT_FAILURE);
    }
  else if (covMatrix.IsSymmetric() == false)
    {
      cout << "[ComputeEfficiency::Fit2DEfficiencies]\tFit status: covariance matrix not symmetric" << endl;
      exit (EXIT_FAILURE);
    }

  cTestGlobalFit->cd();
  effFuncs2D[q2BinIndx]->Draw("surf1 fb");
  hisFunc2D->Draw("lego2 fb");
  cTestGlobalFit->Update();

  cout << "@@@ chi2/DoF = " << effFuncs2D[q2BinIndx]->GetChisquare() / effFuncs2D[q2BinIndx]->GetNDF() << " (" << effFuncs2D[q2BinIndx]->GetChisquare() << "/" << effFuncs2D[q2BinIndx]->GetNDF() << ")";
  cout << "\tCL : " << TMath::Prob(effFuncs2D[q2BinIndx]->GetChisquare(),effFuncs2D[q2BinIndx]->GetNDF()) << " @@@" << endl;


  // ################################################
  // # Check if analytical efficiency goes negative #
  // ################################################
  if (Utility->EffMinValue2D(cosThetaKBins,cosThetaLBins,effFuncs2D[q2BinIndx]) < 0.0) { cout << "NEGATIVE EFFICIENCY !" << endl; exit (EXIT_FAILURE); }
  Utility->SaveAnalyticalEff(fileNameOut.c_str(),effFuncs2D[q2BinIndx],(q2Bins->operator[](q2BinIndx) + q2Bins->operator[](q2BinIndx+1)) / 2.,q2Bins);
  fileNameOut.replace(fileNameOut.find(".txt"),4,"FullCovariance.txt");
  Utility->SaveAnalyticalEffFullCovariance(fileNameOut.c_str(),&covMatrix,(q2Bins->operator[](q2BinIndx) + q2Bins->operator[](q2BinIndx+1)) / 2.,q2Bins);


  // ########################################
  // # Check integrity of covariance matrix #
  // ########################################
  vector<TMatrixTSym<double>*>* covMatrices = new vector<TMatrixTSym<double>*>;
  Utility->ReadAnalyticalEffFullCovariance(fileNameOut.c_str(),covMatrices,"2D",0);


  // #############
  // # Save plot #
  // #############
  if (savePlot == true)
    {
      myString.clear(); myString.str("");
      myString << "Eff2D_q2Bin_" << q2BinIndx << ".pdf";
      cTestGlobalFit->Print(myString.str().c_str());
    }


  effFuncs2D.clear();
  covMatrix.Clear();
  covMatrices->clear();
  delete covMatrices;
}


void Fit3DEfficiencies (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
			unsigned int SignalType, Utils::effStruct myEff, unsigned int q2BinIndx, string fileNameOut, bool savePlot)
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadRightMargin(0.04);


  // ###################
  // # Local variables #
  // ###################
  ifstream fileInput;
  TFitResultPtr fitResults;
  stringstream myString;
  string tmpString;
  double coeff;
  TH1* tmpHist1D;
  TH2* tmpHist2D;
  TF1* effFunc1D;
  vector<TF2*> effFuncs2D;
  TF3* effFunc3D;
  TH3D* effHis3D;
  double Zmin, Zmax;
  // #######################
  double parErrorPhi = 1e-4;
  double Yaxes       = 2e-1;
  double Zscale      = 1.1;
  // #######################


  TCanvas* cShow2DAnaEff = new TCanvas("cShow2DAnaEff", "cShow2DAnaEff", 10, 10, 900, 500);
  cShow2DAnaEff->Divide(3,1);

  TCanvas* cTestGlobalFit2D = new TCanvas("cTestGlobalFit2D", "cTestGlobalFit2D", 10, 10, 900, 500);
  cTestGlobalFit2D->Divide(3,1);

  TCanvas* cTestGlobalFit1D = new TCanvas("cTestGlobalFit1D", "cTestGlobalFit11D", 10, 10, 900, 500);
  cTestGlobalFit1D->Divide(3,1);

  
  // ##########################
  // # Read binned efficiency #
  // ##########################
  myString.clear(); myString.str("");
  myString << "Histo_q2Bin_" << q2BinIndx;
  TH3D* hisFunc3D = Utility->Get3DEffHitoq2Bin(myString.str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2BinIndx,myEff);
  Zmin = hisFunc3D->GetMinimum();
 
  cShow2DAnaEff->cd(1);
  hisFunc3D->Draw("ISO");
  cTestGlobalFit2D->cd(1);
  hisFunc3D->Draw("ISO");
  cTestGlobalFit1D->cd(1);
  hisFunc3D->Draw("ISO");


  // #################################
  // # Read analytical efficiency    #
  // #################################
  // # cos(theta_k) and cos(theta_l) #
  // #################################
  Utility->ReadAnalyticalEff(INPUT_THETAL_THETAK,q2Bins,cosThetaKBins,cosThetaLBins,&effFuncs2D,"effFuncs2D",0);

  // #######
  // # phi #
  // #######
  effFunc1D = new TF1("effFunc1D",Utility->TellMeEffFuncPhi().c_str(),Utility->PI,Utility->PI);
  Utility->InitEffFuncPhi(effFunc1D,q2BinIndx);

  cout << "\n@@@ Reading coefficients for analytical efficiency for phi from file " << INPUT_PHI << " @@@" << endl;
  fileInput.open(INPUT_PHI,ofstream::in);
  if (fileInput.good() == false)
    {
      cout << "[ComputeEfficiency::Fit3DEfficiencies]\tError opening file : " << INPUT_PHI << endl;
      exit (EXIT_FAILURE);
    }
  for (unsigned int j = 0; j < q2BinIndx; j++) getline(fileInput,tmpString);
  getline(fileInput,tmpString);
  stringstream rawStringK(tmpString);
  rawStringK >> coeff; // Discard q2 bin value
  for (unsigned int l = 0; l < Utility->NcoeffPhi; l++)
    {
      rawStringK >> coeff;
      effFunc1D->SetParameter(l,coeff);
      rawStringK >> coeff;
      effFunc1D->SetParError(l,coeff);

      cout << "Reading coef. " << l << " for var. phi: " << effFunc1D->GetParameter(l) << "+/-" << effFunc1D->GetParError(l) << endl;
    }


  // ########################################################
  // # Putting together cos(theta_k), cos(theta_l), and phi #
  // ########################################################
  effFunc3D = new TF3("effFunc3D",Utility->TellMeEffFuncThetaKThetaLPhi().c_str(),
		      cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),
		      cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1),
		      phiBins->operator[](0)      ,phiBins->operator[](phiBins->size()-1)            );
  effFunc3D->SetMarkerStyle(20);
  effFunc3D->SetMarkerColor(kBlack);

  for (int i = 0; i < effFuncs2D[q2BinIndx]->GetNpar(); i++)
    {
      effFunc3D->SetParameter(i,effFuncs2D[q2BinIndx]->GetParameter(i));
      effFunc3D->SetParError(i,effFuncs2D[q2BinIndx]->GetParError(i));
    }

  for (int i = 0; i < effFunc1D->GetNpar(); i++)
    {
      effFunc3D->SetParameter(effFuncs2D[q2BinIndx]->GetNpar()+i,effFunc1D->GetParameter(i));
      effFunc3D->SetParError(effFuncs2D[q2BinIndx]->GetNpar()+i,effFunc1D->GetParError(i));
    }


  // ##################################################################
  // # Add constraint along X or Y or Z to bound the function at zero #
  // ##################################################################
  Utility->AddConstraintThetaKThetaLPhi(&hisFunc3D,q2BinIndx,SignalType,q2BinIndx);
  hisFunc3D->SetMarkerStyle(20);
  hisFunc3D->SetMarkerColor(kBlack);
  hisFunc3D->GetXaxis()->SetTitleOffset(1.35);
  hisFunc3D->GetYaxis()->SetTitleOffset(1.35);
  hisFunc3D->GetZaxis()->SetTitleOffset(1.35);


  // #########################################################################
  // # Perform the fit of the analytical efficiency to the binned efficiency #
  // #########################################################################
  fitResults = hisFunc3D->Fit(effFunc3D->GetName(),"VMS");
  cout << "[ComputeEfficiency::Fit3DEfficiencies]\tFit status: " << fitResults << endl;

  if (fitResults != 0)
    {
      fitResults = hisFunc3D->Fit(effFunc3D->GetName(),"VS");
      cout << "[ComputeEfficiency::Fit3DEfficiencies]\tFit status: " << fitResults << endl;
    }

  if (fitResults != 0)
    {
      for (int i = 1; i < effFunc1D->GetNpar(); i++)
	{
	  effFunc3D->SetParameter(effFuncs2D[q2BinIndx]->GetNpar()+i,0.0);
	  effFunc3D->SetParError(effFuncs2D[q2BinIndx]->GetNpar()+i,parErrorPhi);
	}

      fitResults = hisFunc3D->Fit(effFunc3D->GetName(),"VS");
      cout << "[ComputeEfficiency::Fit3DEfficiencies]\tFit status: " << fitResults << endl;
    }

  TMatrixTSym<double> covMatrix(fitResults->GetCovarianceMatrix());

  if (fitResults != 0) exit (EXIT_FAILURE);
  else if (covMatrix.Determinant() <= 0.0)
    {
      cout << "[ComputeEfficiency::Fit3DEfficiencies]\tFit status: covariance matrix with negative determinant" << endl;
      exit (EXIT_FAILURE);
    }
  else if (covMatrix.IsSymmetric() == false)
    {
      cout << "[ComputeEfficiency::Fit3DEfficiencies]\tFit status: covariance matrix not symmetric" << endl;
      exit (EXIT_FAILURE);
    }

  effHis3D = Utility->Get3DEffHitoq2Bin("effHis3D",q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2BinIndx,myEff);
  effHis3D->SetMarkerStyle(22);
  effHis3D->SetMarkerColor(kRed);
  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
    for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  coeff = effFunc3D->Eval((cosThetaKBins->operator[](j)+cosThetaKBins->operator[](j+1))/2.,
				  (cosThetaLBins->operator[](k)+cosThetaLBins->operator[](k+1))/2.,
				  (phiBins->operator[](l)+phiBins->operator[](l+1))/2.);
	  effHis3D->SetBinContent(j+1,k+1,l+1,coeff);
	  effHis3D->SetBinError(j+1,k+1,l+1,0.0);
	}

  // #######################
  // # Make 2D projections #
  // #######################
  cTestGlobalFit2D->cd(1);
  tmpHist2D = static_cast<TH2D*>(hisFunc3D->Project3D("xy"));
  Zmax = tmpHist2D->GetMaximum()*Zscale;
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("lego2 fb");

  tmpHist2D = static_cast<TH2D*>(effHis3D->Project3D("xy"));
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("same surf fb");

  cShow2DAnaEff->cd(1);
  tmpHist2D = static_cast<TH2D*>(effHis3D->Project3D("xy"));
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("surf1 fb");

  cTestGlobalFit2D->cd(2);
  tmpHist2D = static_cast<TH2D*>(hisFunc3D->Project3D("xz"));
  Zmax = tmpHist2D->GetMaximum()*Zscale;
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("lego2 fb");

  tmpHist2D = static_cast<TH2D*>(effHis3D->Project3D("xz"));
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("same surf fb");

  cShow2DAnaEff->cd(2);
  tmpHist2D = static_cast<TH2D*>(effHis3D->Project3D("xz"));
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("surf1 fb");

  cTestGlobalFit2D->cd(3);
  tmpHist2D = static_cast<TH2D*>(hisFunc3D->Project3D("yz"));
  Zmax = tmpHist2D->GetMaximum()*Zscale;
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("lego2 fb");

  tmpHist2D = static_cast<TH2D*>(effHis3D->Project3D("yz"));
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("same surf fb");

  cShow2DAnaEff->cd(3);
  tmpHist2D = static_cast<TH2D*>(effHis3D->Project3D("yz"));
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("surf1 fb");

  cTestGlobalFit2D->Update();
  cShow2DAnaEff->Update();

  // #######################
  // # Make 1D projections #
  // #######################
  cTestGlobalFit1D->cd(1);
  tmpHist1D = hisFunc3D->Project3D("x");
  tmpHist1D->GetYaxis()->SetRangeUser(0.0,Yaxes);
  tmpHist1D->Draw("pe1");
  effHis3D->Project3D("x")->Draw("same pe1");

  cTestGlobalFit1D->cd(2);
  tmpHist1D = hisFunc3D->Project3D("y");
  tmpHist1D->GetYaxis()->SetRangeUser(0.0,Yaxes);
  tmpHist1D->Draw("pe1");
  effHis3D->Project3D("y")->Draw("same pe1");

  cTestGlobalFit1D->cd(3);
  tmpHist1D = hisFunc3D->Project3D("z");
  tmpHist1D->GetYaxis()->SetRangeUser(0.0,Yaxes);
  tmpHist1D->Draw("pe1");
  effHis3D->Project3D("z")->Draw("same pe1");

  TLegend* leg = new TLegend(0.2, 0.8, 0.85, 0.95, "");
  leg->AddEntry(hisFunc3D,"Binned efficiency");
  leg->AddEntry(effHis3D,"Analytical efficiency");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  cTestGlobalFit1D->Update();

  cout << "@@@ chi2/DoF = " << fitResults->Chi2() / fitResults->Ndf() << " (" << fitResults->Chi2() << "/" << fitResults->Ndf() << ")";
  cout << "\tCL : " << TMath::Prob(fitResults->Chi2(),fitResults->Ndf()) << " @@@" << endl;


  // ################################################
  // # Check if analytical efficiency goes negative #
  // ################################################
  if (Utility->EffMinValue3D(cosThetaKBins,cosThetaLBins,phiBins,effFunc3D) < 0.0) { cout << "NEGATIVE EFFICIENCY !" << endl; exit (EXIT_FAILURE); }
  Utility->SaveAnalyticalEff(fileNameOut.c_str(),effFunc3D,(q2Bins->operator[](q2BinIndx) + q2Bins->operator[](q2BinIndx+1)) / 2.,q2Bins);
  fileNameOut.replace(fileNameOut.find(".txt"),4,"FullCovariance.txt");
  Utility->SaveAnalyticalEffFullCovariance(fileNameOut.c_str(),&covMatrix,(q2Bins->operator[](q2BinIndx) + q2Bins->operator[](q2BinIndx+1)) / 2.,q2Bins);
  

  // ########################################
  // # Check integrity of covariance matrix #
  // ########################################
  vector<TMatrixTSym<double>*>* covMatrices = new vector<TMatrixTSym<double>*>;
  Utility->ReadAnalyticalEffFullCovariance(fileNameOut.c_str(),covMatrices,"3D",0);


  // #############
  // # Save plot #
  // #############
  if (savePlot == true)
    {
      myString.clear(); myString.str("");
      myString << "Eff1D_q2Bin_" << q2BinIndx << ".pdf";
      cTestGlobalFit1D->Print(myString.str().c_str());

      myString.clear(); myString.str("");
      myString << "Eff2D_q2Bin_" << q2BinIndx << ".pdf";
      cTestGlobalFit2D->Print(myString.str().c_str());

      myString.clear(); myString.str("");
      myString << "Eff2Danaly_q2Bin_" << q2BinIndx << ".pdf";
      cShow2DAnaEff->Print(myString.str().c_str());
    }


  covMatrix.Clear();
  covMatrices->clear();
  delete covMatrices;
  delete effHis3D;
}


void Test2DEfficiency (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, Utils::effStruct myEff, unsigned int q2BinIndx, bool savePlot)
{
  // ###################
  // # Local variables #
  // ###################
  stringstream myString;
  string tmpString;
  vector<TF2*> effFuncs2D;
  TF2* effFunc2D;
  vector<TF12*> effFuncSlice;
  vector<TH1D*> histoSlice;
  // ###################


  TCanvas* cEff = new TCanvas("cEff", "cEff", 10, 10, 1200, 800);
  cEff->Divide(7,2);


  // ##########################
  // # Read binned efficiency #
  // ##########################
  myString.clear(); myString.str("");
  myString << "Histo_q2Bin_" << q2BinIndx;
  TH2D* hisFunc2D = Utility->Get2DEffHitoq2Bin(myString.str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2BinIndx,myEff);

  cEff->cd(1);
  hisFunc2D->Draw("lego2 fb");


  // ##############################
  // # Read analytical efficiency #
  // ##############################
  Utility->ReadAnalyticalEff(TEST_THETAL_THETAK,q2Bins,cosThetaKBins,cosThetaLBins,&effFuncs2D,"effFuncs2D",0);
  effFunc2D = effFuncs2D[q2BinIndx];
  cEff->cd(2);
  effFunc2D->GetZaxis()->SetTitleOffset(1.25);
  effFunc2D->Draw("surf2 fb");
  Utility->EffMinValue2D(cosThetaKBins,cosThetaLBins,effFunc2D);


  // ##############################
  // # Show slices of cos(thetaK) #
  // ##############################
  for (unsigned int binIndx = 0; binIndx < cosThetaKBins->size()-1; binIndx++)
    {
      cEff->cd(binIndx+3);

      myString.clear(); myString.str("");
      myString << "effFuncSlice_binK_" << binIndx;
      effFuncSlice.push_back(new TF12(myString.str().c_str(),effFunc2D,(cosThetaKBins->operator[](binIndx)+cosThetaKBins->operator[](binIndx+1))/2.,"y"));
      effFuncSlice.back()->GetXaxis()->SetTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
      effFuncSlice.back()->GetYaxis()->SetTitle("Efficiency");
      effFuncSlice.back()->GetYaxis()->SetRangeUser(0.0,ordinateRange);
      effFuncSlice.back()->Draw();

      myString.clear(); myString.str("");
      myString << "histoSlice_binK_" << binIndx;
      histoSlice.push_back(hisFunc2D->ProjectionY(myString.str().c_str(),binIndx+1,binIndx+1));
      histoSlice.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
      histoSlice.back()->SetMarkerStyle(20);
      histoSlice.back()->SetYTitle("Efficiency");
      histoSlice.back()->GetYaxis()->SetRangeUser(0.0,ordinateRange);
      histoSlice.back()->Draw("same pe1");
    }


  // ##############################
  // # Show slices of cos(thetaL) #
  // ##############################
  for (unsigned int binIndx = 0; binIndx < cosThetaLBins->size()-1; binIndx++)
    {
      cEff->cd(binIndx+3+cosThetaKBins->size()-1);

      myString.clear(); myString.str("");
      myString << "effFuncSlice_binL_" << binIndx;
      effFuncSlice.push_back(new TF12(myString.str().c_str(),effFunc2D,(cosThetaLBins->operator[](binIndx)+cosThetaLBins->operator[](binIndx+1))/2.,"x"));
      effFuncSlice.back()->GetXaxis()->SetTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
      effFuncSlice.back()->GetYaxis()->SetTitle("Efficiency");
      effFuncSlice.back()->GetYaxis()->SetRangeUser(0.0,ordinateRange);
      effFuncSlice.back()->Draw();

      myString.clear(); myString.str("");
      myString << "histoSlice_binL_" << binIndx;
      histoSlice.push_back(hisFunc2D->ProjectionX(myString.str().c_str(),binIndx+1,binIndx+1));
      histoSlice.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
      histoSlice.back()->SetMarkerStyle(20);
      histoSlice.back()->SetYTitle("Efficiency");
      histoSlice.back()->GetYaxis()->SetRangeUser(0.0,ordinateRange);
      histoSlice.back()->Draw("same pe1");
    }


  cEff->cd(1);
  hisFunc2D->Draw("lego2 fb");
  cEff->Update();
  if (savePlot == true)
    {
      tmpString = TEST_THETAL_THETAK;
      tmpString.erase(tmpString.find(".txt"),4);
      myString.clear(); myString.str("");
      myString << tmpString << "_" << q2BinIndx << ".pdf";
      cEff->Print(myString.str().c_str());
    }
}


void Test3DEfficiency (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, Utils::effStruct myEff, unsigned int q2BinIndx, bool savePlot)
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadBottomMargin(0.12);


  // ###################
  // # Local variables #
  // ###################
  stringstream myString;
  string tmpString;
  double coeff;
  TH1* tmpHist1D;
  TH2* tmpHist2D;
  vector<TF3*> effFuncs3D;
  TF3 *EffFunc3D;
  TH3D* effHis3D;
  double Zmin, Zmax;
  // ##################
  double Yaxes  = 2e-1;
  double Zscale =  1.1;
  // ##################


  TCanvas* cEff = new TCanvas("cEff", "cEff", 10, 10, 1000, 1200);
  cEff->Divide(3,3);


  // ##########################
  // # Read binned efficiency #
  // ##########################
  myString.clear(); myString.str("");
  myString << "Histo_q2Bin_" << q2BinIndx;
  TH3D* hisFunc3D = Utility->Get3DEffHitoq2Bin(myString.str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2BinIndx,myEff);
  Zmin = hisFunc3D->GetMinimum();
  hisFunc3D->SetMarkerStyle(20);
  hisFunc3D->SetMarkerColor(kBlack);
  hisFunc3D->GetXaxis()->SetTitleOffset(1.4);
  hisFunc3D->GetYaxis()->SetTitleOffset(1.4);
  hisFunc3D->GetZaxis()->SetTitleOffset(1.4);

 
  // ##############################
  // # Read analytical efficiency #
  // ##############################
  Utility->ReadAnalyticalEff(TEST_THETAL_THETAK_PHI,q2Bins,cosThetaKBins,cosThetaLBins,phiBins,&effFuncs3D,"effFuncs3D",0);
  EffFunc3D = effFuncs3D[q2BinIndx];
  Utility->EffMinValue3D(cosThetaKBins,cosThetaLBins,phiBins,EffFunc3D);
  
  effHis3D = Utility->Get3DEffHitoq2Bin("effHis3D",q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2BinIndx,myEff);
  effHis3D->SetMarkerStyle(22);
  effHis3D->SetMarkerColor(kRed);
  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
    for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  coeff = EffFunc3D->Eval((cosThetaKBins->operator[](j)+cosThetaKBins->operator[](j+1))/2.,
				  (cosThetaLBins->operator[](k)+cosThetaLBins->operator[](k+1))/2.,
				  (phiBins->operator[](l)+phiBins->operator[](l+1))/2.);
	  
	  effHis3D->SetBinContent(j+1,k+1,l+1,coeff);
	  effHis3D->SetBinError(j+1,k+1,l+1,0.0);
	}


  // #######################
  // # Make 2D projections #
  // #######################
  cEff->cd(1);
  tmpHist2D = static_cast<TH2D*>(hisFunc3D->Project3D("xy"));
  Zmax = tmpHist2D->GetMaximum()*Zscale;
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("lego2 fb");
  tmpHist2D = static_cast<TH2D*>(effHis3D->Project3D("xy"));
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("same surf fb");

  cEff->cd(4);
  tmpHist2D->Draw("surf2 fb");

  cEff->cd(2);
  tmpHist2D = static_cast<TH2D*>(hisFunc3D->Project3D("xz"));
  Zmax = tmpHist2D->GetMaximum()*Zscale;
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("lego2 fb");
  tmpHist2D = static_cast<TH2D*>(effHis3D->Project3D("xz"));
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("same surf fb");

  cEff->cd(5);
  tmpHist2D->Draw("surf2 fb");

  cEff->cd(3);
  tmpHist2D = static_cast<TH2D*>(hisFunc3D->Project3D("yz"));
  Zmax = tmpHist2D->GetMaximum()*Zscale;
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("lego2 fb");
  tmpHist2D = static_cast<TH2D*>(effHis3D->Project3D("yz"));
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("same surf fb");

  cEff->cd(6);
  tmpHist2D->Draw("surf2 fb");


  // #######################
  // # Make 1D projections #
  // #######################
  cEff->cd(7);
  tmpHist1D = hisFunc3D->Project3D("x");
  tmpHist1D->GetYaxis()->SetRangeUser(0.0,Yaxes);
  tmpHist1D->Draw("pe1");
  effHis3D->Project3D("x")->Draw("same pe1");

  cEff->cd(8);
  tmpHist1D = hisFunc3D->Project3D("y");
  tmpHist1D->GetYaxis()->SetRangeUser(0.0,Yaxes);
  tmpHist1D->Draw("pe1");
  effHis3D->Project3D("y")->Draw("same pe1");

  cEff->cd(9);
  tmpHist1D = hisFunc3D->Project3D("z");
  tmpHist1D->GetYaxis()->SetRangeUser(0.0,Yaxes);
  tmpHist1D->Draw("pe1");
  effHis3D->Project3D("z")->Draw("same pe1");


  cEff->Update();
  if (savePlot == true)
    {
      tmpString = TEST_THETAL_THETAK_PHI;
      tmpString.erase(tmpString.find(".txt"),4);
      myString.clear(); myString.str("");
      myString << tmpString << "_" << q2BinIndx << ".pdf";
      cEff->Print(myString.str().c_str());
    }

  delete effHis3D;
}


int main (int argc, char** argv)
{
  if (argc >= 2)
    {
      string option = argv[1];


      cout << "\n@@@ Settings @@@" << endl;
      cout << "INPUT_THETAL: "           << INPUT_THETAL << endl;
      cout << "INPUT_PHI: "              << INPUT_PHI << endl;
      cout << "INPUT_THETAL_THETAK: "    << INPUT_THETAL_THETAK << endl;
      cout << "TEST_THETAL_THETAK: "     << TEST_THETAL_THETAK << endl;
      cout << "TEST_THETAL_THETAK_PHI: " << TEST_THETAL_THETAK_PHI << endl;

      cout << "\nRIGHTtag: "     << RIGHTtag << endl;
      cout << "SavePlot: "       << SavePlot << endl;
      cout << "CHECKEFFatREAD: " << CHECKEFFatREAD << endl;
      cout << "NFILES: "         << NFILES << endl;
      cout << "INPUTGenEff: "    << INPUTGenEff << endl;
      cout << "SETBATCH: "       << SETBATCH << endl;
      cout << "ParameterFILE: "  << ParameterFILE << endl;
      cout << "ordinateRange: "  << ordinateRange << endl;


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
      

      if ((option == "Make") && (argc == 7))
	{
	  string fileNameGenCandidatesNoFilter = argv[2];
	  string fileNameRecoCandidates        = argv[3];
	  string fileNameSingleCand            = argv[4];
	  string fileNameOutput                = argv[5];
	  string SignalType                    = argv[6];

	  if (SETBATCH == true)
	    {
	      cout << "\n@@@ Setting batch mode @@@" << endl;
	      gROOT->SetBatch(true);
	    }
	  TApplication theApp ("Applications", &argc, argv);

      
	  TFile* NtplFileGenCandidatesNoFilter = new TFile(fileNameGenCandidatesNoFilter.c_str(), "READ");
	  theTreeGenCandidatesNoFilter         = (TTree*) NtplFileGenCandidatesNoFilter->Get("B0KstMuMu/B0KstMuMuNTuple");
	  NTupleGenCandidatesNoFilter          = new B0KstMuMuSingleCandTreeContent();
	  NTupleGenCandidatesNoFilter->Init();

	  TFile* NtplFileRecoCandidates = new TFile(fileNameRecoCandidates.c_str(), "READ");
	  theTreeRecoCandidates         = (TTree*) NtplFileRecoCandidates->Get("B0KstMuMu/B0KstMuMuNTuple");
	  NTupleRecoCandidates          = new B0KstMuMuSingleCandTreeContent();
	  NTupleRecoCandidates->Init();

	  TFile* NtplFileSingleCand = new TFile(fileNameSingleCand.c_str(), "READ");
	  theTreeSingleCand         = (TTree*) NtplFileSingleCand->Get("B0KstMuMu/B0KstMuMuNTuple");
	  NTupleSingleCand          = new B0KstMuMuSingleCandTreeContent();
	  NTupleSingleCand->Init();


	  Utility = new Utils(RIGHTtag);
	  if (Utility->RIGHTflavorTAG == true) Utility->ReadBins(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"goodTag");
	  else                                 Utility->ReadBins(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"misTag");
	  Utility->ReadPreselectionCut(ParameterFILE);
	  Utility->ReadSelectionCuts(ParameterFILE);


	  // ###################################
	  // # Initialize efficiency structure #
	  // ###################################
	  Utility->GenEfficiency(&myEff,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins);
	  Utils::effValue myEffVal;
	  Utility->ResetEffValue(&myEffVal,0.0);
	  Utility->InitEfficiency(myEffVal,myEff,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins);


	  ComputeEfficiency(theTreeGenCandidatesNoFilter,NTupleGenCandidatesNoFilter,myEff.Den1,myEff.Err2PoisDen1,myEff.Err2WeigDen1, 1 ,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,atoi(SignalType.c_str()));
	  ComputeEfficiency(theTreeGenCandidatesNoFilter,NTupleGenCandidatesNoFilter,myEff.Num1,myEff.Err2PoisNum1,myEff.Err2WeigNum1, 2 ,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,atoi(SignalType.c_str()));
	  ComputeEfficiency(theTreeRecoCandidates,       NTupleRecoCandidates,       myEff.Den2,myEff.Err2PoisDen2,myEff.Err2WeigDen2, 3 ,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,atoi(SignalType.c_str()));
	  ComputeEfficiency(theTreeSingleCand,           NTupleSingleCand,           myEff.Num2,myEff.Err2PoisNum2,myEff.Err2WeigNum2, 4 ,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,atoi(SignalType.c_str()));


	  Utility->SaveEfficiency(fileNameOutput.c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,myEff);
	  if (SETBATCH == false) MakeHistogramsAllBins(&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,myEff);
	  cout << "\n@@@ Efficiency computation is done @@@" << endl;


	  delete Utility;
	  if (SETBATCH == false) theApp.Run (); // Eventloop on air
	  return EXIT_SUCCESS;
	}
      else if (((option == "ReadBin") || (option == "Read3DAnaly")) && (argc >= 3))
	{
	  string fileNameInput = argv[2];
	  int specBin = -1;
	  if (argc == 4) specBin = atoi(argv[3]);

	  TApplication theApp ("Applications", &argc, argv);

    
	  Utility = new Utils(RIGHTtag);
	  if (Utility->RIGHTflavorTAG == true) Utility->ReadBins(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"goodTag");
	  else                                 Utility->ReadBins(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"misTag");
	  if (option == "ReadBin") Read3DEfficiencies(true,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,fileNameInput,false,&myEff,CHECKEFFatREAD,SavePlot,specBin);
	  else                     Read3DEfficiencies(true,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,fileNameInput,true,&myEff,CHECKEFFatREAD,SavePlot,specBin);


	  delete Utility;
	  theApp.Run (); // Eventloop on air
	  return EXIT_SUCCESS;
	}
      else if ((option == "Read3DGenAnaly") && (argc >= 3))
	{
	  string fileNameInput = argv[2];
	  int specBin = -1;
	  if (argc == 4) specBin = atoi(argv[3]);

	  TApplication theApp ("Applications", &argc, argv);


	  Utility = new Utils(RIGHTtag);
	  if (Utility->RIGHTflavorTAG == true) Utility->ReadBins(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"goodTag");
	  else                                 Utility->ReadBins(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"misTag");
	  Read3DEfficiencies(false,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,fileNameInput,true,&myEff,CHECKEFFatREAD,SavePlot,specBin);


	  delete Utility;
	  theApp.Run (); // Eventloop on air
	  return EXIT_SUCCESS;
	}
      else if (((option == "Fit1DEff") && (argc == 6)) || ((option == "Fit2DEff") && (argc == 5)) || ((option == "Fit3DEff") && (argc == 5)))
	{
	  string SignalType      = argv[2];
	  string fileNameInput   = argv[3];
	  unsigned int q2BinIndx = atoi(argv[4]);
	  string whichVar2Fit = "";
	  if (option == "Fit1DEff") whichVar2Fit = argv[5];
	  cout << "Which Var to Fit: " << whichVar2Fit << endl;

	  TApplication theApp ("Applications", &argc, argv);


	  Utility = new Utils(RIGHTtag);
	  if (Utility->RIGHTflavorTAG == true) Utility->ReadBins(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"goodTag");
	  else                                 Utility->ReadBins(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"misTag");
	  Utility->ReadEfficiency(fileNameInput.c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,&myEff);

	  if      (option == "Fit1DEff") Fit1DEfficiencies(&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,atoi(SignalType.c_str()),myEff,whichVar2Fit,q2BinIndx,"Theta.txt");
	  else if (option == "Fit2DEff") Fit2DEfficiencies(&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,atoi(SignalType.c_str()),myEff,q2BinIndx,"ThetaKThetaL.txt",SavePlot);
	  else if (option == "Fit3DEff") Fit3DEfficiencies(&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,atoi(SignalType.c_str()),myEff,q2BinIndx,"ThetaKThetaLPhi.txt",SavePlot);


	  delete Utility;
	  theApp.Run (); // Eventloop on air
	  return EXIT_SUCCESS;
	}
      else if (((option == "Test2DEff") || (option == "Test3DEff")) && (argc == 4))
	{
	  string fileNameInput   = argv[2];
	  unsigned int q2BinIndx = atoi(argv[3]);

	  TApplication theApp ("Applications", &argc, argv);


	  Utility = new Utils(RIGHTtag);
	  if (Utility->RIGHTflavorTAG == true) Utility->ReadBins(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"goodTag");
	  else                                 Utility->ReadBins(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"misTag");
	  Utility->ReadEfficiency(fileNameInput.c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,&myEff);
	  
	  if      (option == "Test2DEff") Test2DEfficiency(&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,myEff,q2BinIndx,SavePlot);
	  else if (option == "Test3DEff") Test3DEfficiency(&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,myEff,q2BinIndx,SavePlot);


	  delete Utility;
	  theApp.Run (); // Eventloop on air
	  return EXIT_SUCCESS;
	}
      else
	{
	  cout << "Efficiency = [#ev. after filter (GEN-level) / #ev. generated (GEN-level before filter)] * [#ev. in single cand. (reco-level truth-matched) / #ev. after filter (GEN-level in multi cand.)]" << endl;
	  cout << "1. process inputFileGenCandidatesNoFilter.root with AddVars2Candidates nvGen" << endl;
	  cout << "2. process inputRecoCandidates.root with AddVars2Candidates nvGen" << endl;
	  cout << "3. make inputFileSingleCand.root from inputRecoCandidates.root" << endl; 

	  cout << "\nParameter missing: " << endl;
	  cout << "./ComputeEfficiency [Make ReadBin Read3DAnaly Read3DGenAnaly Fit1DEff Fit2DEff Fit3DEff Test2DEff Test3DEff] " << endl;
	  cout << "[inputFileGenCandidatesNoFilter.root inputRecoCandidates.root inputFileSingleCand.root] " << endl;
	  cout << "[out/in]putFile.txt [SignalType] [q2 bin indx.]" << endl;
	  cout << "SignalType : if B0 --> K*0 mumu : 1; if B0 --> J/psi K*0 : 3; if B0 --> psi(2S) K*0 : 5" << endl;

	  cout << "\nMake           --> root files for efficiency computation AND outputFile.txt AND SignalType" << endl;

	  cout << "ReadBin        --> (read file with binned eff.) file with binned efficiency AND [q2 bin indx.(optional)]" << endl;
	  cout << "Read3DAnaly    --> (read file with analytical eff.) file with analytical efficiency AND [q2 bin indx.(optional)]" << endl;
	  
	  cout << "Read3DGenAnaly --> (read files generated from analytical eff.) file with analytical efficiency AND [q2 bin indx.(optional)]" << endl;

	  cout << "Fit1DEff       --> SignalType AND file with binned efficiency AND [q2 bin indx.] [thetaL thetaK phi]" << endl;
	  cout << "Fit2DEff       --> SignalType AND file with binned efficiency AND [q2 bin indx.]" << endl;
	  cout << "Fit3DEff       --> SignalType AND file with binned efficiency AND [q2 bin indx.]" << endl;

	  cout << "Test2DEff      --> file with binned efficiency AND [q2 bin indx.]" << endl;
	  cout << "Test3DEff      --> file with binned efficiency AND [q2 bin indx.]" << endl;

	  return EXIT_FAILURE;
	}
    }
  else
    {
      cout << "Efficiency = [#ev. after filter (GEN-level) / #ev. generated (GEN-level before filter)] * [#ev. in single cand. (reco-level truth-matched) / #ev. after filter (GEN-level in multi cand.)]" << endl;
      cout << "1. process inputFileGenCandidatesNoFilter.root with AddVars2Candidates nvGen" << endl;
      cout << "2. process inputRecoCandidates.root with AddVars2Candidates nvGen" << endl;
      cout << "3. make inputFileSingleCand.root from inputRecoCandidates.root" << endl; 

      cout << "\nParameter missing: " << endl;
      cout << "./ComputeEfficiency [Make ReadBin Read3DAnaly Read3DGenAnaly Fit1DEff Fit2DEff Fit3DEff Test2DEff Test3DEff] " << endl;
      cout << "[inputFileGenCandidatesNoFilter.root inputRecoCandidates.root inputFileSingleCand.root] " << endl;
      cout << "[out/in]putFile.txt [SignalType] [q2 bin indx.]" << endl;
      cout << "SignalType : if B0 --> K*0 mumu : 1; if B0 --> J/psi K*0 : 3; if B0 --> psi(2S) K*0 : 5" << endl;

      cout << "\nMake           --> root files for efficiency computation AND outputFile.txt AND SignalType" << endl;

      cout << "ReadBin        --> (read file with binned eff.) file with binned efficiency AND [q2 bin indx.(optional)]" << endl;
      cout << "Read3DAnaly    --> (read file with analytical eff.) file with analytical efficiency AND [q2 bin indx.(optional)]" << endl;
	  
      cout << "Read3DGenAnaly --> (read files generated from analytical eff.) file with analytical efficiency AND [q2 bin indx.(optional)]" << endl;

      cout << "Fit1DEff       --> SignalType AND file with binned efficiency AND [q2 bin indx.] [thetaL thetaK phi]" << endl;
      cout << "Fit2DEff       --> SignalType AND file with binned efficiency AND [q2 bin indx.]" << endl;
      cout << "Fit3DEff       --> SignalType AND file with binned efficiency AND [q2 bin indx.]" << endl;

      cout << "Test2DEff      --> file with binned efficiency AND [q2 bin indx.]" << endl;
      cout << "Test3DEff      --> file with binned efficiency AND [q2 bin indx.]" << endl;

      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

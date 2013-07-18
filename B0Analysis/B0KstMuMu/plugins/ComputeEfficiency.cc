// #########################################################################
// # Program to compute the efficiency for the B0 --> K*0 mu+ mu- analysis #
// # in bins of dimuon q^2, cos(theta_K), cos(theta_l), and phi            #
// #########################################################################
// # Author: Mauro Dinardo                                                 #
// #########################################################################

#ifndef __CINT__
#include <TROOT.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TF12.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TText.h>
#include <TFitResult.h>
#include <TVectorD.h>
#include <TDecompBK.h>
#include <TMatrixDSymEigen.h>
#endif

#include <math.h>
#include <iostream>
#include <sstream>

#include "ReadParameters.h"
#include "Utils.h"
#include "B0KstMuMuSingleCandTreeContent.h"

using namespace std;


// ###########################################
// # How to create the analytical efficiency #
// ###########################################
// (a) Make binned efficiency with the command "Make"
// (b) Fit the theta_l variable with the command "FitEff1D" for every q2 bin
//     - Set command option to "thetaL"
// (c) Fit the theta_K variable with the command "FitEff1D"
//     - Set command option to "thetaK"
//     - Set "INPUTTHETAL" to the output.txt file of point (b)
// (d) Fit the theta_l-theta_K variables with the command "FitEff2D"
//     - Set "INPUT2DEffRef" to the output.txt file of point (c)
// (e) Of you want to look at the final result run the command "TestEff"
//     - Set "INPUT2DEffRef" to the output.txt file of point (d)


// ####################
// # Global constants #
// ####################
#define INPUTTHETAL    "../../Efficiency/ThetaL_B0ToKstMuMu.txt"
#define INPUT2DEffRef  "../../Efficiency/ThetaK_B0ToKstMuMu.txt" // "ThetaK_B0ToKstMuMu.txt" OR "ThetaKThetaL_B0ToKstMuMu.txt"
#define INPUT2DEffComp "../../Efficiency/ThetaKThetaL_B0ToKstMuMu.txt"
#define SavePlot false
#define CHECKEFFatREAD false // Check if 2D efficiency goes negative
#define NFILES 200
#define INPUTGenEff "../../Efficiency/EffRndGenAnalyFilesSign_JPsi_Psi2S/Efficiency_RndGen.txt"
#define SignalType 1 // If checking MC B0 --> K*0 mumu  : 1
                     // If checking MC B0 --> J/psi K*0 : 3
                     // If checking MC B0 --> psi(2S) K*0 : 5
#define ParameterFILE "../python/ParameterFile.txt"

// ############################################################################
// # In case the final number of events with and without filter are different #
// ############################################################################
// For signal MC
#define CorrFactorNEvGenNoFilter_KstMuMu 9985833333.0 // Total number of actually generated events without filter
#define CorrFactorNEvGenFilter_KstMuMu   9943333333.0 // Total number of actually generated events with filter
// For J/psi normalization / control sample MC
#define CorrFactorNEvGenNoFilter_JPsi 1499500000.0 // Total number of actually generated events without filter
#define CorrFactorNEvGenFilter_JPsi  14995000000.0 // Total number of actually generated events with filter
// For psi(2S) normalization / control sample MC
#define CorrFactorNEvGenNoFilter_Psi2S 1499000000.0 // Total number of actually generated events without filter
#define CorrFactorNEvGenFilter_Psi2S  14985000000.0 // Total number of actually generated events with filter

// ###################
// # Fit constraints #
// ###################
#define abscissaErr   1.0e-2
#define ordinateVal   1.0e-5
#define ordinateErr   1.0e-5
#define ordinateRange 2.2e-3


// ####################
// # Global variables #
// ####################
Utils* Utility;

TTree* theTreeGenCandidatesNoFilter;
TTree* theTreeGenCandidatesFilter;
TTree* theTreeRecoCandidates;
TTree* theTreeSingleCand;
B0KstMuMuSingleCandTreeContent* NTupleGenCandidatesNoFilter;
B0KstMuMuSingleCandTreeContent* NTupleGenCandidatesFilter;
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
void ComputeEfficiency (TTree* theTree, B0KstMuMuSingleCandTreeContent* NTuple, double* Vector, double* VectorErr2Pois, double* VectorErr2Weig, unsigned int type,
			vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins);
void MakeHistogramsAllBins (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, Utils::effStruct myEff);
void ReadEfficiencies (bool isSingleEff, vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
		       string fileNameInput, bool isAnalyEff, Utils::effStruct* myEff, bool CheckEffatRead, bool saveHistos, int specBin = -1);
void Fit1DEfficiencies (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
			Utils::effStruct myEff, string who, unsigned int q2BinIndx, string fileNameOut);
void Fit2DEfficiencies (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
			Utils::effStruct myEff, unsigned int q2BinIndx, string fileNameOut);
void TestEfficiency (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, Utils::effStruct myEff, unsigned int q2BinIndx, bool saveHistos);


// ###########################
// # Function Implementation #
// ###########################
void ComputeEfficiency (TTree* theTree, B0KstMuMuSingleCandTreeContent* NTuple, double* Vector, double* VectorErr2Pois, double* VectorErr2Weig, unsigned int type,
			vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins)
// ##########################################################
// # Efficiency type = 1 --> total gen events before filter #
// # Efficiency type = 2 --> total gen events after filter  #
// # Efficiency type = 3 --> total reco events              #
// # Efficiency type = 4 --> total single candidate events  #
// ##########################################################
{
  // #################
  // Local variables #
  // #################
  int nEntries;
  Utils::effStruct _Counter;
  double* Counter;
  double CorrFactorNEvGenNoFilter = 1.0;
  double CorrFactorNEvGenFilter   = 1.0;
  double mumuq2;
  double cosThetaK;
  double cosThetaMu;
  double phiKstMuMuPlane;
  int mumuq2BinIndx;
  int cosThetaKBinIndx;
  int cosThetaMuBinIndx;
  int phiKstMuMuPlaneBinIndx;
  // #################


  if (SignalType == 1)
    {
      CorrFactorNEvGenNoFilter = CorrFactorNEvGenNoFilter_KstMuMu;
      CorrFactorNEvGenFilter   = CorrFactorNEvGenFilter_KstMuMu;
    }
  else if (SignalType == 3)
    {
      CorrFactorNEvGenNoFilter = CorrFactorNEvGenNoFilter_JPsi;
      CorrFactorNEvGenFilter   = CorrFactorNEvGenFilter_JPsi;
    }
  else if (SignalType == 5)
    {
      CorrFactorNEvGenNoFilter = CorrFactorNEvGenNoFilter_Psi2S;
      CorrFactorNEvGenFilter   = CorrFactorNEvGenFilter_Psi2S;
    }
  else
    {
      cout << "[ComputeEfficiency::ComputeEfficiency]\tNon valid signal type : " << SignalType << endl;
      exit (1);
    }
  cout << "\n@@@ Correction factor for efficiency between with/without filter: " << CorrFactorNEvGenFilter << " / " << CorrFactorNEvGenNoFilter;
  cout << " = " << CorrFactorNEvGenFilter/CorrFactorNEvGenNoFilter << " @@@" << endl;


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

      if ((NTuple->B0pT > Utility->GetSeleCut("B0pT")) && (fabs(NTuple->B0Eta) < Utility->GetSeleCut("B0Eta")) &&
	  ((NTuple->genSignal == SignalType || NTuple->genSignal == SignalType+1) &&
	   ((type == 1 || type == 2 || type == 3) || ((type == 4) && (NTuple->truthMatchSignal->at(0) == true)))))
	{
	  mumuq2          = NTuple->mumuMass->at(0)*NTuple->mumuMass->at(0);

	  cosThetaK       = NTuple->CosThetaKArb;
	  cosThetaMu      = NTuple->CosThetaMuArb;
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
		  	 mumuq2BinIndx] = 
		    Vector[phiKstMuMuPlaneBinIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
		  	   cosThetaMuBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
		  	   cosThetaKBinIndx*(q2Bins->size()-1) +
		  	   mumuq2BinIndx] +
		    NTuple->evWeight;
	  
		  Counter[phiKstMuMuPlaneBinIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
			  cosThetaMuBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
			  cosThetaKBinIndx*(q2Bins->size()-1) +
			  mumuq2BinIndx]++;

		  VectorErr2Weig[phiKstMuMuPlaneBinIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaMuBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaKBinIndx*(q2Bins->size()-1) +
				 mumuq2BinIndx] =
		    VectorErr2Weig[phiKstMuMuPlaneBinIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				   cosThetaMuBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				   cosThetaKBinIndx*(q2Bins->size()-1) +
				   mumuq2BinIndx] +
		    NTuple->evWeightE2*NTuple->evWeightE2;
		}
	    }
	}
    }

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	for (unsigned int l = 0; l < phiBins->size()-1; l++)
	  {
	    if (type == 1)
	      {
		VectorErr2Pois[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] =
		  Vector[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] *
		  pow(CorrFactorNEvGenFilter/CorrFactorNEvGenNoFilter,2.0);
		
		Vector[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] =
		  Vector[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] *
		  CorrFactorNEvGenFilter/CorrFactorNEvGenNoFilter;
	      }

	    else if (type == 2)
	      VectorErr2Pois[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] =
		Vector[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	    else if (Counter[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] <= 0.0)
	      {
		VectorErr2Weig[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = 0.0;
		VectorErr2Pois[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = 0.0;
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
  // #################
  // Local variables #
  // #################
  stringstream myString;
  double Eff, EffErr;
  double totalEffAll;
  double totalEffSignal;
  double totalEffJPsi;
  double totalEffPsiP;
  double num2;
  double* q2Bins_        = Utility->MakeBinning(q2Bins);
  double* cosThetaKBins_ = Utility->MakeBinning(cosThetaKBins);
  double* cosThetaLBins_ = Utility->MakeBinning(cosThetaLBins);
  double* phiBins_       = Utility->MakeBinning(phiBins);
  // ##################
  double Xaxes = 1.0e4;
  // ##################


  TCanvas* cEff = new TCanvas("cEff", "cEff", 10, 10, 1600, 900);
  cEff->Divide(2,2);

  TH1D* Hq2 = new TH1D("Hq2", "Hq2", q2Bins->size()-1, q2Bins_);
  Hq2->SetXTitle("q^{2} (GeV^{2})");
  Hq2->SetYTitle("Efficiency");

  vector<TH1D*> vecHcosThetaK;
  vector<TH1D*> vecHcosThetaL;
  vector<TH1D*> vecHphi;
  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    {
      myString.str("");
      myString << "vecHcosThetaK_" << i;
      vecHcosThetaK.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), cosThetaKBins->size()-1, cosThetaKBins_));
      vecHcosThetaK.back()->SetXTitle("cos(#theta_{K})");
      vecHcosThetaK.back()->SetYTitle("Efficiency");

      myString.str("");
      myString << "vecHcosThetaL_" << i;
      vecHcosThetaL.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), cosThetaLBins->size()-1, cosThetaLBins_));
      vecHcosThetaL.back()->SetXTitle("cos(#theta_{l})");
      vecHcosThetaL.back()->SetYTitle("Efficiency");

      myString.str("");
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
  myString.str(""); myString << "q^{2} bin 0";
  legThetaK->AddEntry(vecHcosThetaK[0],myString.str().c_str());
  for (unsigned int i = 1; i < q2Bins->size()-1; i++)
    {
      vecHcosThetaK[i]->SetMarkerStyle(20+i);
      vecHcosThetaK[i]->SetMarkerColor(1+i);
      vecHcosThetaK[i]->SetLineColor(1+i);
      vecHcosThetaK[i]->SetLineWidth(2);
      vecHcosThetaK[i]->Draw("sames e1");
      vecHcosThetaK[i]->GetYaxis()->SetRangeUser(0.0,ordinateRange);
      myString.str(""); myString << "q^{2} bin " << i;
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
  myString.str(""); myString << "q^{2} bin 0";
  legThetaL->AddEntry(vecHcosThetaL[0],myString.str().c_str());
  for (unsigned int i = 1; i < q2Bins->size()-1; i++)
    {
      vecHcosThetaL[i]->SetMarkerStyle(20+i);
      vecHcosThetaL[i]->SetMarkerColor(1+i);
      vecHcosThetaL[i]->SetLineColor(1+i);
      vecHcosThetaL[i]->SetLineWidth(2);
      vecHcosThetaL[i]->Draw("sames e1");
      vecHcosThetaL[i]->GetYaxis()->SetRangeUser(0.0,ordinateRange);
      myString.str(""); myString << "q^{2} bin " << i;
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
  TLegend* legphi = new TLegend(0.88, 0.65, 0.97, 0.89, "");
  myString.str(""); myString << "q^{2} bin 0";
  legphi->AddEntry(vecHphi[0],myString.str().c_str());
  for (unsigned int i = 1; i < q2Bins->size()-1; i++)
    {
      vecHphi[i]->SetMarkerStyle(20+i);
      vecHphi[i]->SetMarkerColor(1+i);
      vecHphi[i]->SetLineColor(1+i);
      vecHphi[i]->SetLineWidth(2);
      vecHphi[i]->Draw("sames e1");
      vecHphi[i]->GetYaxis()->SetRangeUser(0.0,ordinateRange);
      myString.str(""); myString << "q^{2} bin " << i;
      legphi->AddEntry(vecHphi[i],myString.str().c_str());
    }
  legphi->SetFillColor(0);
  legphi->SetBorderSize(0);
  legphi->Draw();

  cEff->Update();


  TCanvas* cNumCosThetaL = new TCanvas("cNumCosThetaL", "cNumCosThetaL", 10, 10, 1600, 900);
  cNumCosThetaL->Divide(2,3);

  TCanvas* cNumCosThetaK = new TCanvas("cNumCosThetaK", "cNumCosThetaK", 10, 10, 1600, 900);
  cNumCosThetaK->Divide(2,3);

  vector<TH1D*> vecHq2ANDcosThetaK;
  vector<TH1D*> vecHq2ANDcosThetaL;
  vector<TLegend*> legHq2ANDcosThetaK;
  vector<TLegend*> legHq2ANDcosThetaL;
  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    {
      for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	{
	  myString.str("");
	  myString << "vecHq2ANDcosThetaK_" << i << "_" << j;
	  vecHq2ANDcosThetaK.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), cosThetaLBins->size()-1, cosThetaLBins_));
	  vecHq2ANDcosThetaK.back()->SetXTitle("cos(#theta_{l})");
	  vecHq2ANDcosThetaK.back()->SetYTitle("Numerator single cand. [#]");

	  legHq2ANDcosThetaK.push_back(new TLegend(0.88, 0.6, 0.97, 0.89, ""));
	}

      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	{
	  myString.str("");
	  myString << "vecHq2ANDcosThetaL_" << i << "_" << k;
	  vecHq2ANDcosThetaL.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), cosThetaKBins->size()-1, cosThetaKBins_));
	  vecHq2ANDcosThetaL.back()->SetXTitle("cos(#theta_{K})");
	  vecHq2ANDcosThetaL.back()->SetYTitle("Numerator single cand. [#]");

	  legHq2ANDcosThetaL.push_back(new TLegend(0.88, 0.6, 0.97, 0.89, ""));
	}
    }
  
  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
    {
      for (unsigned int i = 0; i < q2Bins->size()-1; i++)
	{
	  // #######################################################################################################
	  // # Fill histogram : numerator of the efficiency vs cos(theta_l) in bins of dimuon q^2 and cos(theta_K) #
	  // #######################################################################################################
	  for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	    {
	      num2 = 0.0;
	      for (unsigned int l = 0; l < phiBins->size()-1; l++)
		num2 = num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->SetBinContent(k+1,num2);
	    }
	  
	  cNumCosThetaL->cd(j+1);
	  vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->SetMarkerStyle(20+i);
	  vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->SetMarkerColor(1+i);
	  vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->SetLineColor(1+i);
	  vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->SetLineWidth(2);
	  if (i == 0) vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->Draw("e1");
	  else vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->Draw("sames e1");
	  vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->GetYaxis()->SetRangeUser(1.0,Xaxes);

	  myString.str(""); myString << "q^{2} bin " << i;
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
	    {
	      num2 = 0.0;
	      for (unsigned int l = 0; l < phiBins->size()-1; l++)
		num2 = num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->SetBinContent(j+1,num2);
	    }
	  
	  cNumCosThetaK->cd(k+1);
	  vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->SetMarkerStyle(20+i);
	  vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->SetMarkerColor(1+i);
	  vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->SetLineColor(1+i);
	  vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->SetLineWidth(2);
	  if (i == 0) vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->Draw("e1");
	  else vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->Draw("sames e1");
	  vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->GetYaxis()->SetRangeUser(1.0,Xaxes);

	  myString.str(""); myString << "q^{2} bin " << i;
	  legHq2ANDcosThetaL[k]->AddEntry(vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i],myString.str().c_str());
	}

      legHq2ANDcosThetaL[k]->SetFillColor(0);
      legHq2ANDcosThetaL[k]->SetBorderSize(0);
      legHq2ANDcosThetaL[k]->Draw();
    }

  cNumCosThetaL->Update();
  cNumCosThetaK->Update();


  Utility->IntegrateEffAll(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,myEff,&totalEffAll,&EffErr);
  Utility->IntegrateEffButPsi(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,myEff,&totalEffSignal,&EffErr);
  Utility->IntegrateEffInJPsi(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,myEff,&totalEffJPsi,&EffErr);
  Utility->IntegrateEffInPsiP(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,myEff,&totalEffPsiP,&EffErr);
  cout << "\nTotal efficiency: " << totalEffAll <<  "\tTotal signal efficiency: " << totalEffSignal << "\tTotal J/psi efficiency: " << totalEffJPsi << "\tTotal psi(2S) efficiency: " << totalEffPsiP << endl;
}


void ReadEfficiencies (bool isSingleEff, vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
		       string fileNameInput, bool isAnalyEff, Utils::effStruct* myEff, bool CheckEffatRead, bool saveHistos, int specBin)
{
  // #################
  // Local variables #
  // #################
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
  vector<TF2*> effFuncs;
  string tmpString;
  stringstream myString;
  stringstream effFuncIDs;
  // ###################
  double Xaxes = 0.2e-3;
  double Yaxes = 1.2e-2;
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
      myString.str("");
      myString << "Hq2_" << i;
      Hq2.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), q2Bins->size()-1, q2Bins_));
      Hq2.back()->SetXTitle("q^{2} (GeV^{2})");
      Hq2.back()->SetYTitle("Efficiency");
    }

  vector<TH1D*> HcosThetaK;
  for (unsigned int i = 0; i < MAXVAL; i++)
    {
      myString.str("");
      myString << "HcosThetaK_" << i;
      HcosThetaK.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), cosThetaKBins->size()-1, cosThetaKBins_));
      HcosThetaK.back()->SetXTitle("cos(#theta_{K})");
      HcosThetaK.back()->SetYTitle("Efficiency");
    }

  vector<TH1D*> HcosThetaL;
  for (unsigned int i = 0; i < MAXVAL; i++)
    {
      myString.str("");
      myString << "HcosThetaL_" << i;
      HcosThetaL.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), cosThetaLBins->size()-1, cosThetaLBins_));
      HcosThetaL.back()->SetXTitle("cos(#theta_{l})");
      HcosThetaL.back()->SetYTitle("Efficiency");
    }

  vector<TH1D*> Hphi;
  for (unsigned int i = 0; i < MAXVAL; i++)
    {
      myString.str("");
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
	  else                     Utility->ReadAnalyticalEff(fileNameInput.c_str(),q2Bins,cosThetaKBins,cosThetaLBins,&effFuncs,"effFuncsRef",0);
	}
      else
	{
	  myString.str("");
	  myString << tmpString << "_" << itF << ".txt";
	  effFuncIDs << "effFuncsRef_" << itF;
 	  if (isAnalyEff == false) Utility->ReadEfficiency(myString.str().c_str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,myEff);
	  else                     Utility->ReadAnalyticalEff(myString.str().c_str(),q2Bins,cosThetaKBins,cosThetaLBins,&effFuncs,effFuncIDs.str().c_str(),0);
	}


      // ################################################
      // # Check if analytical efficiency goes negative #
      // ################################################
      if ((isAnalyEff == true) && (CheckEffatRead == true))
      	for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
      	  if (Utility->EffMinValue2D(cosThetaKBins,cosThetaLBins,effFuncs[i]) < 0.0) cout << "@@@ Negative efficiency function #" << itF << " for q2 bin #" << i << " ! @@@" << endl;

      
      for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
	{
	  // #############################################
	  // # Fill histogram : efficiency vs dimuon q^2 #
	  // #############################################
	  if (isAnalyEff == false) Utility->IntegrateEffButq2(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,i,*myEff,&Eff,&EffErr);
	  else
	    {
	      Eff = 0.0;
	      Eff = Eff + effFuncs[i]->Integral(cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1)) /
		((cosThetaKBins->operator[](cosThetaKBins->size()-1)-cosThetaKBins->operator[](0)) * (cosThetaLBins->operator[](cosThetaLBins->size()-1)-cosThetaLBins->operator[](0)));
	      EffErr = ordinateErr;
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
      		Eff = Eff + effFuncs[indx]->Integral(cosThetaKBins->operator[](j),cosThetaKBins->operator[](j+1),cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1)) /
      		  ((cosThetaKBins->operator[](j+1)-cosThetaKBins->operator[](j)) * (cosThetaLBins->operator[](cosThetaLBins->size()-1)-cosThetaLBins->operator[](0)));
      	      EffErr = ordinateErr;
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
      		Eff = Eff + effFuncs[indx]->Integral(cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),cosThetaLBins->operator[](k),cosThetaLBins->operator[](k+1)) /
      		  ((cosThetaKBins->operator[](cosThetaKBins->size()-1)-cosThetaKBins->operator[](0)) * (cosThetaLBins->operator[](k+1)-cosThetaLBins->operator[](k)));
      	      EffErr = ordinateErr;
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
		Eff = Eff + effFuncs[indx]->Integral(cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1)) /
		  ((cosThetaKBins->operator[](cosThetaKBins->size()-1)-cosThetaKBins->operator[](0)) * (cosThetaLBins->operator[](cosThetaLBins->size()-1)-cosThetaLBins->operator[](0)));
	      EffErr = ordinateErr;
	    }
	  Hphi[itF]->SetBinContent(l+1,Eff);
	  Hphi[itF]->SetBinError(l+1,EffErr);
	}
      
      
      if ((isSingleEff == true) || (itF == MAXVAL-1))
	{
	  cEff0->cd();
	  Hq2[itF]->SetMarkerStyle(20);
	  Hq2[itF]->GetYaxis()->SetRangeUser(0,Yaxes);
	  if (isSingleEff == true) Hq2[itF]->Draw("pe1");
	  else Hq2[itF]->Draw("same pe1");
	  
	  cEff1->cd();
	  HcosThetaK[itF]->SetMarkerStyle(20);
	  HcosThetaK[itF]->GetYaxis()->SetRangeUser(0,Yaxes);
	  if (isSingleEff == true) HcosThetaK[itF]->Draw("pe1");
	  else HcosThetaK[itF]->Draw("same pe1");

	  cEff2->cd();
	  HcosThetaL[itF]->SetMarkerStyle(20);
	  HcosThetaL[itF]->GetYaxis()->SetRangeUser(0,Yaxes);
	  if (isSingleEff == true) HcosThetaL[itF]->Draw("pe1");
	  else HcosThetaL[itF]->Draw("same pe1");

	  cEff3->cd();
	  Hphi[itF]->SetMarkerStyle(20);
	  Hphi[itF]->GetYaxis()->SetRangeUser(0,Yaxes);
	  if (isSingleEff == true) Hphi[itF]->Draw("pe1");
	  else Hphi[itF]->Draw("same pe1");
	}
      else
	{
	  cEff0->cd();
	  Hq2[itF]->SetLineColor(kBlue);
	  Hq2[itF]->GetYaxis()->SetRangeUser(0,Yaxes);
	  if (itF == 0) Hq2[itF]->Draw();
	  else Hq2[itF]->Draw("same");

	  cEff1->cd();
	  HcosThetaK[itF]->SetLineColor(kBlue);
	  HcosThetaK[itF]->GetYaxis()->SetRangeUser(0,Yaxes);
	  if (itF == 0) HcosThetaK[itF]->Draw();
	  else HcosThetaK[itF]->Draw("same");

	  cEff2->cd();
	  HcosThetaL[itF]->SetLineColor(kBlue);
	  HcosThetaL[itF]->GetYaxis()->SetRangeUser(0,Yaxes);
	  if (itF == 0) HcosThetaL[itF]->Draw();
	  else HcosThetaL[itF]->Draw("same");

	  cEff3->cd();
	  Hphi[itF]->SetLineColor(kBlue);
	  Hphi[itF]->GetYaxis()->SetRangeUser(0,Yaxes);
	  if (itF == 0) Hphi[itF]->Draw();
	  else Hphi[itF]->Draw("same");
	}
      
      effFuncs.clear();
    }
 
  cEff0->Update();
  cEff1->Update();
  cEff2->Update();
  cEff3->Update();
  
  if (saveHistos == true)
    {
      if (isAnalyEff == false) cEff0->Print("Histo1.pdf");
      cEff1->Print("Histo2.pdf");
      cEff2->Print("Histo3.pdf");
      if (isAnalyEff == false) cEff3->Print("Histo4.pdf");
    }


  TCanvas* cStatK = new TCanvas("cStatK", "cStatK", 10, 10, 1600, 900);
  cStatK->Divide(2,3);

  TCanvas* cStatL = new TCanvas("cStatL", "cStatL", 10, 10, 1600, 900);
  cStatL->Divide(2,3);

  vector<TH1D*> HstatK;
  for (unsigned int i = 0; i < cosThetaKBins->size()-1; i++)
    {
      myString.str("");
      myString << "HstatK_" << i;
      HstatK.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), 50, -Xaxes, Xaxes));
      HstatK.back()->SetXTitle("Diff. to reference");
      HstatK.back()->SetYTitle("Entries [#]");

      cStatK->cd(i+1);
      HstatK.back()->Draw();
    }

  vector<TH1D*> HstatL;
  for (unsigned int i = 0; i < cosThetaLBins->size()-1; i++)
    {
      myString.str("");
      myString << "HstatL_" << i;
      HstatL.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), 50, -Xaxes, Xaxes));
      HstatL.back()->SetXTitle("Diff. to reference");
      HstatL.back()->SetYTitle("Entries [#]");

      cStatL->cd(i+1);
      HstatL.back()->Draw();
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

      cStatK->Update();
      cStatL->Update();
    }


  Utility->IntegrateEffAll(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,*myEff,&totalEffAll,&EffErr);
  Utility->IntegrateEffButPsi(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,*myEff,&totalEffSignal,&EffErr);
  Utility->IntegrateEffInJPsi(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,*myEff,&totalEffJPsi,&EffErr);
  Utility->IntegrateEffInPsiP(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,*myEff,&totalEffPsiP,&EffErr);
  cout << "\nTotal efficiency: " << totalEffAll <<  "\tTotal signal efficiency: " << totalEffSignal << "\tTotal J/psi efficiency: " << totalEffJPsi << "\tTotal psi(2S) efficiency: " << totalEffPsiP << endl;


  if ((isSingleEff == true) && (isAnalyEff == false)) MakeHistogramsAllBins(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,*myEff);
}


void Fit1DEfficiencies (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
			Utils::effStruct myEff, string who, unsigned int q2BinIndx, string fileNameOut)
{
  // #################
  // Local variables #
  // #################
  double Eff, EffErr;
  double* cosThetaKBins_ = new double[cosThetaKBins->size()]; for (unsigned int i = 0; i < cosThetaKBins->size(); i++) cosThetaKBins_[i] = cosThetaKBins->operator[](i);
  double* cosThetaLBins_ = new double[cosThetaLBins->size()]; for (unsigned int i = 0; i < cosThetaLBins->size(); i++) cosThetaLBins_[i] = cosThetaLBins->operator[](i);
  vector<TH1D*> histFit;
  TF1* fitFun;
  stringstream myString;
  string tmpString;
  double coeff;
  ofstream fileOutput;
  ifstream fileInput;
  bool nullParameter;
  // #########################
  double YaxesThetaL = 4.0e-3;
  double YaxesThetaK = 1.0e-2;
  // #########################


  TCanvas* cEff = new TCanvas("cEff", "cEff", 10, 10, 1600, 900);
  cEff->Divide(3,2);

  if (who == "thetaL")
    {
      fileOutput.open(fileNameOut.replace(fileNameOut.find(".txt"),4,"L.txt").c_str(),ofstream::app);
      if (fileOutput.good() == false)
	{
	  cout << "[ComputeEfficiency::Fit1DEfficiencies]\tError opening file : " << fileNameOut.c_str() << endl;
	  exit (1);
	}


      // ###########################
      // # Define the fit function #
      // ###########################
      fitFun = new TF1("fitFun",Utility->TellMeEffFuncThetaL().c_str(),-1.0 - abscissaErr,1.0 + abscissaErr);


      for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	{
	  Utility->InitEffFuncThetaL(fitFun,q2BinIndx);


	  myString.str("");
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
	  Utility->AddConstraintThetaL(&histFit.back(),q2BinIndx,j,abscissaErr,ordinateVal,ordinateErr,q2BinIndx);


	  // ################
	  // # Save results #
	  // ################
	  histFit.back()->SetMarkerStyle(20);
	  histFit.back()->SetXTitle("cos(#theta_{l})");
	  histFit.back()->SetYTitle("Efficiency");
	  histFit.back()->GetYaxis()->SetRangeUser(-YaxesThetaL,YaxesThetaL);

	  cEff->cd(j+1);
	  histFit.back()->Fit("fitFun","V");
	  histFit.back()->Draw("pe1");
	  
	  fileOutput << (q2Bins->operator[](q2BinIndx)+q2Bins->operator[](q2BinIndx+1))/2.;
	  for (int i = 0; i < fitFun->GetNpar(); i++) fileOutput << "   " << fitFun->GetParameter(i) << "   " << fitFun->GetParError(i);
	  fileOutput << endl;

	  cout << "\n@@@ Fit for for bin n." << j << " @@@" << endl;
	  cout << "@@@ Value at " << cosThetaLBins->operator[](0) << " : " << fitFun->Eval(cosThetaLBins->operator[](0)) << " @@@" << endl;
	  cout << "@@@ Value at " << cosThetaLBins->operator[](cosThetaLBins->size()-1) << " : " << fitFun->Eval(cosThetaLBins->operator[](cosThetaLBins->size()-1)) << " @@@\n" << endl;

	  if (Utility->EffMinValue1D(cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1),fitFun) < 0.0) { cout << "NEGATIVE EFFICIENCY !" << endl; exit(1); }
	}

      fileOutput.close();
    }
  else if (who == "thetaK")
    {
      fileOutput.open(fileNameOut.replace(fileNameOut.find(".txt"),4,"K.txt").c_str(),ofstream::app);
      if (fileOutput.good() == false)
	{
	  cout << "[ComputeEfficiency::Fit1DEfficiencies]\tError opening file : " << fileNameOut.c_str() << endl;
	  exit (1);
	}


      // ###########################
      // # Define the fit function #
      // ###########################
      fitFun = new TF1("fitFun",Utility->TellMeEffFuncThetaK().c_str(),-1.0 - abscissaErr,1.0 + abscissaErr);


      for (unsigned int k = 0; k < Utility->NcoeffThetaL; k++)
	{
	  Utility->InitEffFuncThetaK(fitFun,q2BinIndx);


	  myString.str("");
	  myString << "histFit_" << k;
	  histFit.push_back(new TH1D(myString.str().c_str(),myString.str().c_str(),cosThetaKBins->size()-1,cosThetaKBins_));


	  // ###############################
	  // # Read coefficients from file #
	  // ###############################
	  fileInput.open(INPUTTHETAL,ofstream::in);
	  if (fileInput.good() == false)
	    {
	      cout << "[ComputeEfficiency::Fit1DEfficiencies]\tError opening file : " << INPUTTHETAL << endl;
	      exit (1);
	    }
	  for (unsigned int j = 0; j < (cosThetaKBins->size()-1)*q2BinIndx; j++) getline(fileInput,tmpString);
	  nullParameter = true;
	  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	    {
	      getline(fileInput,tmpString);
	      stringstream rawStringK(tmpString);
	      rawStringK >> coeff;
	      for (unsigned int i = 0; i < k*2; i++) rawStringK >> coeff;
	      
	      rawStringK >> coeff;
	      histFit.back()->SetBinContent(j+1,coeff);
	      rawStringK >> coeff;
	      histFit.back()->SetBinError(j+1,coeff);

	      cout << "Cos(theta_K) bin " << j << " --> reading coef. " << k << " for var. theta_l: " << histFit.back()->GetBinContent(j+1) << "+/-" <<  histFit.back()->GetBinError(j+1) << endl;

	      if (coeff == 0.0) nullParameter = nullParameter*true; // --> All parameter errors must be zero in order to ignore the parameter
	      else nullParameter = false;
	    }


	  // ################
	  // # Save results #
	  // ################
	  histFit.back()->SetMarkerStyle(20);
	  histFit.back()->SetXTitle("cos(#theta_{K})");
	  histFit.back()->SetYTitle("Efficiency");
	  histFit.back()->GetYaxis()->SetRangeUser(-YaxesThetaK,YaxesThetaK);

	  cEff->cd(k+1);
	  histFit.back()->Fit("fitFun","V");
	  histFit.back()->Draw("pe1");

	  fileOutput << (q2Bins->operator[](q2BinIndx)+q2Bins->operator[](q2BinIndx+1))/2.;
	  for (int i = 0; i < fitFun->GetNpar(); i++) fileOutput << "   " << fitFun->GetParameter(i) << "   " << (nullParameter == false ? fitFun->GetParError(i) : 0.0);
	  fileOutput << endl;

	  cout << "\n@@@ Fit for for bin n." << k << " @@@" << endl;
	  cout << "@@@ Value at " << cosThetaKBins->operator[](0) << " : " << fitFun->Eval(cosThetaKBins->operator[](0)) << " @@@" << endl;
	  cout << "@@@ Value at " << cosThetaKBins->operator[](cosThetaKBins->size()-1) << " : " << fitFun->Eval(cosThetaKBins->operator[](cosThetaKBins->size()-1)) << " @@@\n" << endl;

	  fileInput.close();
	}

      fileOutput.close();
    }
  

  cEff->Update();
  histFit.clear();
}


void Fit2DEfficiencies (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
			Utils::effStruct myEff, unsigned int q2BinIndx, string fileNameOut)
{
  // #################
  // Local variables #
  // #################
  TFitResultPtr fitResults;
  stringstream myString;
  double Eff, EffErr;
  double* cosThetaKBins_ = Utility->MakeBinning(cosThetaKBins);
  double* cosThetaLBins_ = Utility->MakeBinning(cosThetaLBins);
  vector<TF2*> effFuncsRef;
  // #################


  TCanvas* cTestGlobalFit = new TCanvas("cTestGlobalFit", "cTestGlobalFit", 10, 10, 700, 500);

  myString.str("");
  myString << "Histo_q2Bin_" << q2BinIndx;
  TH2D* Histo = new TH2D(myString.str().c_str(), myString.str().c_str(), cosThetaKBins->size()-1, cosThetaKBins_, cosThetaLBins->size()-1, cosThetaLBins_);
  Histo->SetXTitle("cos(#theta_{K})");
  Histo->GetXaxis()->SetTitleOffset(1.8);
  Histo->SetYTitle("cos(#theta_{l})");
  Histo->GetYaxis()->SetTitleOffset(1.8);
  Histo->SetZTitle("Efficiency");


  // ##########################
  // # Read binned efficiency #
  // ##########################
  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
    for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
      {
	Utility->IntegrateEffPhi(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2Bins->operator[](q2BinIndx),cosThetaKBins->operator[](j),cosThetaLBins->operator[](k),myEff,&Eff,&EffErr,false);
	Histo->SetBinContent(j+1,k+1,Eff);
	Histo->SetBinError(j+1,k+1,EffErr);
      }


  Utility->ReadAnalyticalEff(INPUT2DEffRef,q2Bins,cosThetaKBins,cosThetaLBins,&effFuncsRef,"effFuncsRef",0);
  effFuncsRef[q2BinIndx]->SetRange(cosThetaKBins->operator[](0),
				   cosThetaLBins->operator[](0) - abscissaErr,
				   cosThetaKBins->operator[](cosThetaKBins->size()-1) + abscissaErr,
				   cosThetaLBins->operator[](cosThetaLBins->size()-1) + abscissaErr);


  for (int i = 0; i < effFuncsRef[q2BinIndx]->GetNpar(); i++)
    if (effFuncsRef[q2BinIndx]->GetParError(i) == 0.0) effFuncsRef[q2BinIndx]->FixParameter(i,effFuncsRef[q2BinIndx]->GetParameter(i));


  // ############################################################################################
  // # Add constraint along Y (= cosThetaL) where it is necessary to bound the function at zero #
  // ############################################################################################
  Utility->AddConstraintThetaK(&Histo,cosThetaKBins,q2BinIndx,abscissaErr,ordinateVal,ordinateErr,q2BinIndx);


  cTestGlobalFit->cd();
  fitResults = Histo->Fit(effFuncsRef[q2BinIndx]->GetName(),"VMRS");
  TMatrixTSym<double> covMatrix(fitResults->GetCovarianceMatrix());
  effFuncsRef[q2BinIndx]->Draw("surf1");
  Histo->Draw("lego2");
  cTestGlobalFit->Update();

  cout << "@@@ chi2/DoF = " << effFuncsRef[q2BinIndx]->GetChisquare() / effFuncsRef[q2BinIndx]->GetNDF() << " (" << effFuncsRef[q2BinIndx]->GetChisquare() << "/" << effFuncsRef[q2BinIndx]->GetNDF() << ")";
  cout << "\tCL : " << TMath::Prob(effFuncsRef[q2BinIndx]->GetChisquare(),effFuncsRef[q2BinIndx]->GetNDF()) << " @@@" << endl;


  // ############################################################################################
  // # Add constraint along X (= cosThetaK) where it is necessary to bound the function at zero #
  // ############################################################################################
  if (Utility->EffMinValue2D(cosThetaKBins,cosThetaLBins,effFuncsRef[q2BinIndx]) < 0.0)
    {
      cout << "@@@ Efficiency is still negative ! @@@" << endl;

      Utility->AddConstraint2D(&Histo,abscissaErr,ordinateVal,ordinateErr,q2BinIndx,"X");
      cTestGlobalFit->cd();
      fitResults = Histo->Fit(effFuncsRef[q2BinIndx]->GetName(),"VMRS");
      TMatrixTSym<double> covMatrixConstr(fitResults->GetCovarianceMatrix());
      effFuncsRef[q2BinIndx]->Draw("surf1");
      Histo->Draw("lego2");
      cTestGlobalFit->Update();

      cout << "@@@ chi2/DoF = " << effFuncsRef[q2BinIndx]->GetChisquare() / effFuncsRef[q2BinIndx]->GetNDF() << " (" << effFuncsRef[q2BinIndx]->GetChisquare() << "/" << effFuncsRef[q2BinIndx]->GetNDF() << ")";
      cout << "\tCL : " << TMath::Prob(effFuncsRef[q2BinIndx]->GetChisquare(),effFuncsRef[q2BinIndx]->GetNDF()) << " @@@" << endl;
      Utility->EffMinValue2D(cosThetaKBins,cosThetaLBins,effFuncsRef[q2BinIndx]);

      Utility->SaveAnalyticalEff(fileNameOut.c_str(),effFuncsRef[q2BinIndx],(q2Bins->operator[](q2BinIndx) + q2Bins->operator[](q2BinIndx+1)) / 2.,q2Bins);
      fileNameOut.replace(fileNameOut.find(".txt"),4,"FullCovariance.txt");
      Utility->SaveAnalyticalEffFullCovariance(fileNameOut.c_str(),&covMatrixConstr,(q2Bins->operator[](q2BinIndx) + q2Bins->operator[](q2BinIndx+1)) / 2.,q2Bins);

      covMatrixConstr.Clear();

      if (Utility->EffMinValue2D(cosThetaKBins,cosThetaLBins,effFuncsRef[q2BinIndx]) < 0.0) { cout << "NEGATIVE EFFICIENCY !" << endl; exit(1); }
    }
  else
    {
      Utility->SaveAnalyticalEff(fileNameOut.c_str(),effFuncsRef[q2BinIndx],(q2Bins->operator[](q2BinIndx) + q2Bins->operator[](q2BinIndx+1)) / 2.,q2Bins);
      fileNameOut.replace(fileNameOut.find(".txt"),4,"FullCovariance.txt");
      Utility->SaveAnalyticalEffFullCovariance(fileNameOut.c_str(),&covMatrix,(q2Bins->operator[](q2BinIndx) + q2Bins->operator[](q2BinIndx+1)) / 2.,q2Bins);

      covMatrix.Clear();
    }


  // ########################################
  // # Check integrity of covariance matrix #
  // ########################################
  vector<TMatrixTSym<double>*>* covMatrices = new vector<TMatrixTSym<double>*>;
  Utility->ReadAnalyticalEffFullCovariance(fileNameOut.c_str(),covMatrices,0);


  effFuncsRef.clear();
  covMatrices->clear();
  delete covMatrices;
}


void TestEfficiency (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, Utils::effStruct myEff, unsigned int q2BinIndx, bool saveHistos)
{
  // #################
  // Local variables #
  // #################
  double Eff, EffErr;
  double averagePerBin;
  double integralRef;
  double integralComp;
  double integralDiff;
  double chi2point = 0.0;
  double chi2avg   = 0.0;
  double DoF       = 0.0;
  double* cosThetaKBins_ = Utility->MakeBinning(cosThetaKBins);
  double* cosThetaLBins_ = Utility->MakeBinning(cosThetaLBins);
  vector<TF2*> effFuncsRef;
  vector<TF2*> effFuncsComp;
  TF2 *EffFuncRef, *EffFuncComp;
  TF2 *absEffFuncRef, *absEffFuncComp;
  TF2* diffFunc;
  vector<TF12*> effFuncSlice;
  vector<TH1D*> histoSlice;
  string tmpString;
  stringstream myString;
  // ###################
  double Zaxes = 4.0e-3;
  // ###################


  TCanvas* cEffAnaly = new TCanvas("cEffAnaly", "cEffAnaly", 10, 10, 700, 500);
  TCanvas* cEff = new TCanvas("cEff", "cEff", 10, 10, 1600, 900);
  cEff->Divide(7,2);

  myString.str("");
  myString << "Histo_q2Bin_" << q2BinIndx;
  TH2D* Histo = new TH2D(myString.str().c_str(), myString.str().c_str(), cosThetaKBins->size()-1, cosThetaKBins_, cosThetaLBins->size()-1, cosThetaLBins_);
  Histo->SetXTitle("cos(#theta_{K})");
  Histo->GetXaxis()->SetTitleOffset(1.8);
  Histo->SetYTitle("cos(#theta_{l})");
  Histo->GetYaxis()->SetTitleOffset(1.8);
  Histo->SetZTitle("Efficiency");
  Histo->GetZaxis()->SetTitleOffset(1.25);


  // ##########################
  // # Read binned efficiency #
  // ##########################
  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
    for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
      {
	Utility->IntegrateEffPhi(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,
				 q2Bins->operator[](q2BinIndx),cosThetaKBins->operator[](j),cosThetaLBins->operator[](k),myEff,&Eff,&EffErr,false);
	Histo->SetBinContent(j+1,k+1,Eff);
	Histo->SetBinError(j+1,k+1,EffErr);
      }
  cEff->cd(1);
  Histo->Draw("lego2");


  // ###################################################################################
  // # Read analytical efficiencies (reference and comparison) and rescale if required #
  // ###################################################################################
  Utility->ReadAnalyticalEff(INPUT2DEffRef,q2Bins,cosThetaKBins,cosThetaLBins,&effFuncsRef,"effFuncsRef",0);
  Utility->ReadAnalyticalEff(INPUT2DEffComp,q2Bins,cosThetaKBins,cosThetaLBins,&effFuncsComp,"effFuncsComp",0);
  EffFuncRef  = effFuncsRef[q2BinIndx];
  EffFuncComp = effFuncsComp[q2BinIndx];
  cEff->cd(2);
  EffFuncRef->GetZaxis()->SetTitleOffset(1.25);
  EffFuncRef->Draw("surf1");


  // ##############################
  // # Count number of parameters #
  // ##############################
  for (int i = 0; i < EffFuncRef->GetNpar(); i++) if (EffFuncRef->GetParError(i) != 0.0) DoF = DoF + 1.0;


  // #######################################
  // # Compute chi2 Binned - AnalyticalRef #
  // #######################################
  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
    for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
      {
	Eff    = Histo->GetBinContent(j+1,k+1);
	EffErr = Histo->GetBinError(j+1,k+1);
	
	averagePerBin = EffFuncRef->Integral(cosThetaKBins->operator[](j),cosThetaKBins->operator[](j+1),cosThetaLBins->operator[](k),cosThetaLBins->operator[](k+1)) /
	  ((cosThetaKBins->operator[](j+1)-cosThetaKBins->operator[](j)) * (cosThetaLBins->operator[](k+1)-cosThetaLBins->operator[](k)));
	
	chi2avg   = chi2avg + pow((Eff - averagePerBin) / EffErr,2.);
	chi2point = chi2point + pow((Eff - EffFuncRef->Eval((cosThetaKBins->operator[](j)+cosThetaKBins->operator[](j+1))/2.,
							    (cosThetaLBins->operator[](k)+cosThetaLBins->operator[](k+1))/2.)) / EffErr,2.);
      }
  DoF = (cosThetaKBins->size()-1)*(cosThetaLBins->size()-1) - DoF; // DoF = number of bins - number of fit parameters
  cout << "\n@@@ chi2 test between binned and analytical efficiencies @@@" << endl;
  cout << "chi2/DoF (average over the bin) = " << chi2avg / DoF << " (" << chi2avg << "/" << DoF << ")";
  cout << "\tCL : " << TMath::Prob(chi2avg,DoF) << endl;
  cout << "chi2/DoF (by point) = " << chi2point / DoF << " (" << chi2point << "/" << DoF << ")";
  cout << "\tCL : " << TMath::Prob(chi2point,DoF) << endl;
  Utility->EffMinValue2D(cosThetaKBins,cosThetaLBins,EffFuncRef);


  // ##############################
  // # Show slices of cos(thetaK) #
  // ##############################
  for (unsigned int binIndx = 0; binIndx < cosThetaKBins->size()-1; binIndx++)
    {
      cEff->cd(binIndx+3);

      myString.str("");
      myString << "effFuncSlice_binK_" << binIndx;
      effFuncSlice.push_back(new TF12(myString.str().c_str(),EffFuncRef,(cosThetaKBins->operator[](binIndx)+cosThetaKBins->operator[](binIndx+1))/2.,"y"));
      effFuncSlice.back()->GetXaxis()->SetTitle("cos(#theta_{l})");
      effFuncSlice.back()->GetYaxis()->SetTitle("Efficiency");
      effFuncSlice.back()->GetYaxis()->SetRangeUser(0.0,Zaxes);
      effFuncSlice.back()->Draw();

      myString.str("");
      myString << "histoSlice_binK_" << binIndx;
      histoSlice.push_back(Histo->ProjectionY(myString.str().c_str(),binIndx+1,binIndx+1));
      histoSlice.back()->SetXTitle("cos(#theta_{l})");
      histoSlice.back()->SetMarkerStyle(20);
      histoSlice.back()->SetYTitle("Efficiency");
      histoSlice.back()->GetYaxis()->SetRangeUser(0.0,Zaxes);
      histoSlice.back()->Draw("same pe1");
    }


  // ##############################
  // # Show slices of cos(thetaL) #
  // ##############################
  for (unsigned int binIndx = 0; binIndx < cosThetaLBins->size()-1; binIndx++)
    {
      cEff->cd(binIndx+3+cosThetaKBins->size()-1);

      myString.str("");
      myString << "effFuncSlice_binL_" << binIndx;
      effFuncSlice.push_back(new TF12(myString.str().c_str(),EffFuncRef,(cosThetaLBins->operator[](binIndx)+cosThetaLBins->operator[](binIndx+1))/2.,"x"));
      effFuncSlice.back()->GetXaxis()->SetTitle("cos(#theta_{K})");
      effFuncSlice.back()->GetYaxis()->SetTitle("Efficiency");
      effFuncSlice.back()->GetYaxis()->SetRangeUser(0.0,Zaxes);
      effFuncSlice.back()->Draw();

      myString.str("");
      myString << "histoSlice_binL_" << binIndx;
      histoSlice.push_back(Histo->ProjectionX(myString.str().c_str(),binIndx+1,binIndx+1));
      histoSlice.back()->SetXTitle("cos(#theta_{K})");
      histoSlice.back()->SetMarkerStyle(20);
      histoSlice.back()->SetYTitle("Efficiency");
      histoSlice.back()->GetYaxis()->SetRangeUser(0.0,Zaxes);
      histoSlice.back()->Draw("same pe1");
    }


  cEff->cd(1);
  Histo->Draw("lego2");
  cEff->Update();
  if (saveHistos == true)
    {
      tmpString = INPUT2DEffRef;
      tmpString.erase(tmpString.find(".txt"),4);
      myString.str("");
      myString << tmpString << "_" << q2BinIndx << ".pdf";
      cEff->Print(myString.str().c_str());
    }

  cEffAnaly->cd();
  EffFuncRef->Draw("surf1");
  myString.str("");
  myString << "q2 bin = " << q2BinIndx;
  TText* Tag = new TText(0.65,0.85,myString.str().c_str());
  Tag->DrawText(0.65,0.85,Tag->GetTitle());
  cEffAnaly->Update();
  if (saveHistos == true)
    {
      myString.str("");
      myString << "q2Bin_" << q2BinIndx << ".pdf";
      cEffAnaly->Print(myString.str().c_str());
    }


  // ###############################################
  // # Compute area AnalyticalRef - AnalyticalComp #
  // ###############################################
  myString.str("");
  myString << "abs(" << EffFuncRef->GetName() << ")";
  absEffFuncRef = new TF2("absEffFuncRef",myString.str().c_str(),
			  cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),
			  cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1));
  integralRef = absEffFuncRef->Integral(cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1));

  
  myString.str("");
  myString << "abs(" << EffFuncComp->GetName() << ")";
  absEffFuncComp = new TF2("absEffFuncComp",myString.str().c_str(),
			   cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),
			   cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1));
  integralComp = absEffFuncComp->Integral(cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1));


  myString.str("");
  myString << "abs(" << EffFuncRef->GetName() << "*(" << integralRef+integralComp << ")/(" << 2.*integralRef << ") - ";
  myString << EffFuncComp->GetName() << "*(" << integralRef+integralComp << ")/(" << 2.*integralComp << "))";
  diffFunc = new TF2("diffFunc",myString.str().c_str(),
		     cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()),
		     cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()));
  integralDiff = diffFunc->Integral(cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1));


  cout << "\n@@@ Area difference between analytical-reference (" << INPUT2DEffRef << ") and analytical-comparison (" << INPUT2DEffComp << ") efficiencies @@@" << endl;
  cout << "Integral absolute difference: " << integralDiff << ";\tIntegral reference: " << integralRef;
  cout << ";\tPercentage of compatibility: " << (integralRef-integralDiff)/integralRef*100. << "%" << endl;
}


int main (int argc, char** argv)
{
  if (argc >= 2)
    {
      string option = argv[1];


      cout << "\n@@@ Settings @@@" << endl;
      cout << "INPUTTHETAL: "    << INPUTTHETAL << endl;
      cout << "INPUT2DEffRef: "  << INPUT2DEffRef << endl;
      cout << "INPUT2DEffComp: " << INPUT2DEffComp << endl;
      cout << "Save plot: "      << SavePlot << endl;
      cout << "CHECKEFFatREAD: " << CHECKEFFatREAD << endl;
      cout << "NFILES: "         << NFILES << endl;
      cout << "INPUTGenEff: "    << INPUTGenEff << endl;
      cout << "Signal Type: "    << SignalType << endl;
      cout << "ParameterFILE: "  <<  ParameterFILE << endl;

      cout << "\nCorrection factor N events generated without filter (signal):" << CorrFactorNEvGenNoFilter_KstMuMu << endl;
      cout << "Correction factor N events generated with filter (signal):"      << CorrFactorNEvGenFilter_KstMuMu << endl;
      cout << "Correction factor N events generated without filter (J/psi):"    << CorrFactorNEvGenNoFilter_JPsi << endl;
      cout << "Correction factor N events generated with filter (J/psi):"       << CorrFactorNEvGenFilter_JPsi << endl;
      cout << "Correction factor N events generated without filter (psi(2S)):"  << CorrFactorNEvGenNoFilter_Psi2S << endl;
      cout << "Correction factor N events generated with filter (psi(2S)):"     << CorrFactorNEvGenFilter_Psi2S << endl;

      cout << "\nabscissaErr: " << abscissaErr << endl;
      cout << "ordinateVal: "   << ordinateVal << endl;
      cout << "ordinateErr: "   << ordinateErr << endl;
      cout << "ordinateRange: " << ordinateRange << endl;


      if ((option == "Make") && (argc == 7))
	{
	  string fileNameGenCandidatesNoFilter = argv[2];
	  string fileNameGenCandidatesFilter   = argv[3];
	  string fileNameRecoCandidates        = argv[4];
	  string fileNameSingleCand            = argv[5];
	  string fileNameOutput                = argv[6];

	  TApplication theApp ("Applications", &argc, argv);


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

      
	  TFile* NtplFileGenCandidatesNoFilter = new TFile(fileNameGenCandidatesNoFilter.c_str(), "READ");
	  theTreeGenCandidatesNoFilter         = (TTree*) NtplFileGenCandidatesNoFilter->Get("B0SingleCand/B0KstMuMuSingleCandNTuple");
	  NTupleGenCandidatesNoFilter          = new B0KstMuMuSingleCandTreeContent();
	  NTupleGenCandidatesNoFilter->Init();

	  TFile* NtplFileGenCandidatesFilter = new TFile(fileNameGenCandidatesFilter.c_str(), "READ");
	  theTreeGenCandidatesFilter         = (TTree*) NtplFileGenCandidatesFilter->Get("B0SingleCand/B0KstMuMuSingleCandNTuple");
	  NTupleGenCandidatesFilter          = new B0KstMuMuSingleCandTreeContent();
	  NTupleGenCandidatesFilter->Init();

	  TFile* NtplFileRecoCandidates = new TFile(fileNameRecoCandidates.c_str(), "READ");
	  theTreeRecoCandidates         = (TTree*) NtplFileRecoCandidates->Get("B0SingleCand/B0KstMuMuSingleCandNTuple");
	  NTupleRecoCandidates          = new B0KstMuMuSingleCandTreeContent();
	  NTupleRecoCandidates->Init();

	  TFile* NtplFileSingleCand = new TFile(fileNameSingleCand.c_str(), "READ");
	  theTreeSingleCand         = (TTree*) NtplFileSingleCand->Get("B0SingleCand/B0KstMuMuSingleCandNTuple");
	  NTupleSingleCand          = new B0KstMuMuSingleCandTreeContent();
	  NTupleSingleCand->Init();


	  Utility = new Utils();
	  Utility->ReadBins(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins);
	  Utility->ReadSelectionCuts(ParameterFILE);


	  // ###################################
	  // # Initialize efficiency structure #
	  // ###################################
	  Utility->GenEfficiency(&myEff,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins);
	  Utils::effValue myEffVal;
	  Utility->ResetEffValue(&myEffVal,0.0);
	  Utility->InitEfficiency(myEffVal,myEff,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins);


	  ComputeEfficiency(theTreeGenCandidatesNoFilter,NTupleGenCandidatesNoFilter,myEff.Den1,myEff.Err2PoisDen1,myEff.Err2WeigDen1,1,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins);
	  ComputeEfficiency(theTreeGenCandidatesFilter,  NTupleGenCandidatesFilter,  myEff.Num1,myEff.Err2PoisNum1,myEff.Err2WeigNum1,2,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins);
	  ComputeEfficiency(theTreeRecoCandidates,       NTupleRecoCandidates,       myEff.Den2,myEff.Err2PoisDen2,myEff.Err2WeigDen2,3,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins);
	  ComputeEfficiency(theTreeSingleCand,           NTupleSingleCand,           myEff.Num2,myEff.Err2PoisNum2,myEff.Err2WeigNum2,4,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins);


	  Utility->SaveEfficiency(fileNameOutput.c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,myEff);
	  MakeHistogramsAllBins(&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,myEff);
	  cout << "\n@@@ Efficiency computation is done @@@" << endl;


	  delete Utility;
	  theApp.Run (); // Eventloop on air
	  return 0;
	}
      else if ((option == "ReadBin") || (option == "ReadAnaly"))
	{
	  string fileNameInput = argv[2];
	  int specBin = -1;
	  if (argc == 4) specBin = atoi(argv[3]);

	  TApplication theApp ("Applications", &argc, argv);


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

      
	  Utility = new Utils();
	  Utility->ReadBins(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins);

	  if (option == "ReadBin") ReadEfficiencies(true,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,fileNameInput,false,&myEff,CHECKEFFatREAD,SavePlot,specBin);
	  else                     ReadEfficiencies(true,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,fileNameInput,true,&myEff,CHECKEFFatREAD,SavePlot,specBin);


	  delete Utility;
	  theApp.Run (); // Eventloop on air
	  return 0;
	}
      else if (option == "ReadGenAnaly")
	{
	  string fileNameInput = argv[2];
	  int specBin = -1;
	  if (argc == 4) specBin = atoi(argv[3]);

	  TApplication theApp ("Applications", &argc, argv);


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


	  Utility = new Utils();
	  Utility->ReadBins(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins);
	  
	  ReadEfficiencies(false,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,fileNameInput,true,&myEff,CHECKEFFatREAD,SavePlot,specBin);


	  delete Utility;
	  theApp.Run (); // Eventloop on air
	  return 0;
	}
      else if (((option == "FitEff1D") && (argc == 5)) || ((option == "FitEff2D") && (argc == 4)))
	{
	  string fileNameInput   = argv[2];
	  unsigned int q2BinIndx = atoi(argv[3]);
	  string whichVar2Fit = "";
	  if (option == "FitEff1D") whichVar2Fit = argv[4];
	  cout << "Which Var to Fit: " << whichVar2Fit << endl;

	  TApplication theApp ("Applications", &argc, argv);


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


	  Utility = new Utils();
	  Utility->ReadBins(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins);
	  Utility->ReadEfficiency(fileNameInput.c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,&myEff);

	  if      (option == "FitEff1D") Fit1DEfficiencies(&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,myEff,whichVar2Fit,q2BinIndx,"Theta.txt");
	  else if (option == "FitEff2D") Fit2DEfficiencies(&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,myEff,q2BinIndx,"ThetaKThetaL.txt");


	  delete Utility;
	  theApp.Run (); // Eventloop on air
	  return 0;
	}
      else if ((option == "TestEff") && (argc == 4))
	{
	  string fileNameInput   = argv[2];
	  unsigned int q2BinIndx = atoi(argv[3]);

	  TApplication theApp ("Applications", &argc, argv);


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


	  Utility = new Utils();
	  Utility->ReadBins(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins);
	  Utility->ReadEfficiency(fileNameInput.c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,&myEff);
	  
	  TestEfficiency(&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,myEff,q2BinIndx,SavePlot);


	  delete Utility;
	  theApp.Run (); // Eventloop on air
	  return 0;
	}
      else
	{
	  cout << "Efficiency = [#ev. after filter (GEN-level) / #ev. generated (GEN-level before filter)] * [#ev. in single cand. (reco-level truth-matched) / #ev. after filter (GEN-level in multi cand.)]" << endl;
	  cout << "1. process inputFileGenCandidatesNoFilter.root with AddVars2Candidates nvGen" << endl;
	  cout << "2. process inputFileGenCandidatesFilter.root with AddVars2Candidates nvGen" << endl;
	  cout << "3. process inputRecoCandidates.root with AddVars2Candidates nvGen" << endl;
	  cout << "4. make inputFileSingleCand.root from inputRecoCandidates.root" << endl; 

	  cout << "\n@@@@@@@@@@@@ IMPORTANT for [Make] option @@@@@@@@@@@@" << endl;
	  cout << "If you have just the files:" << endl;
	  cout << "- inputFileGenCandidatesNoFilter.root (GEN-level before filter)" << endl;
	  cout << "- inputFileSingleCand.root (reco-level truth-matched)" << endl;
	  cout << "Use the synopsis:" << endl;
	  cout << "./ComputeEfficiency Make inputFileGenCandidatesNoFilter.root inputFileSingleCand.root inputFileSingleCand.root inputFileSingleCand.root outputFile.txt" << endl;
	  cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

	  cout << "\nParameter missing: " << endl;
	  cout << "./ComputeEfficiency [Make ReadBin ReadAnaly ReadGenAnaly FitEff1D FitEff2D TestEff] " << endl;
	  cout << "[inputFileGenCandidatesNoFilter.root inputFileGenCandidatesFilter.root inputRecoCandidates.root inputFileSingleCand.root] " << endl;
	  cout << "[out/in]putFile.txt [q2 bin indx.]" << endl;

	  cout << "Make         --> root files for efficiency computation AND outputFile.txt" << endl;

	  cout << "ReadBin      --> (read file with binned eff.) file with binned efficiency AND [q2 bin indx.(optional)]" << endl;
	  cout << "ReadAnaly    --> (read file with analytical eff.) file with analytical efficiency AND [q2 bin indx.(optional)]" << endl;
	  
	  cout << "ReadGenAnaly --> (read files generated from analytical eff.) file with analytical efficiency AND [q2 bin indx.(optional)]" << endl;

	  cout << "FitEff1D     --> file with binned efficiency AND [q2 bin indx.] [thetaL thetaK]" << endl;
	  cout << "FitEff2D     --> file with binned efficiency AND [q2 bin indx.]" << endl;

	  cout << "TestEff      --> file with binned efficiency AND [q2 bin indx.]" << endl;

	  return 1;
	}
    }
  else
    {
      cout << "Efficiency = [#ev. after filter (GEN-level) / #ev. generated (GEN-level before filter)] * [#ev. in single cand. (reco-level truth-matched) / #ev. after filter (GEN-level in multi cand.)]" << endl;
      cout << "1. process inputFileGenCandidatesNoFilter.root with AddVars2Candidates nvGen" << endl;
      cout << "2. process inputFileGenCandidatesFilter.root with AddVars2Candidates nvGen" << endl;
      cout << "3. process inputRecoCandidates.root with AddVars2Candidates nvGen" << endl;
      cout << "4. make inputFileSingleCand.root from inputRecoCandidates.root" << endl; 

      cout << "\n@@@@@@@@@@@@ IMPORTANT for [Make] option @@@@@@@@@@@@" << endl;
      cout << "If you have just the files:" << endl;
      cout << "- inputFileGenCandidatesNoFilter.root (GEN-level before filter)" << endl;
      cout << "- inputFileSingleCand.root (reco-level truth-matched)" << endl;
      cout << "Use the synopsis:" << endl;
      cout << "./ComputeEfficiency Make inputFileGenCandidatesNoFilter.root inputFileSingleCand.root inputFileSingleCand.root inputFileSingleCand.root outputFile.txt" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

      cout << "\nParameter missing: " << endl;
      cout << "./ComputeEfficiency [Make ReadBin ReadAnaly ReadGenAnaly FitEff1D FitEff2D TestEff] " << endl;
      cout << "[inputFileGenCandidatesNoFilter.root inputFileGenCandidatesFilter.root inputRecoCandidates.root inputFileSingleCand.root] " << endl;
      cout << "[out/in]putFile.txt [q2 bin indx.]" << endl;

      cout << "Make         --> root files for efficiency computation AND outputFile.txt" << endl;
      
      cout << "ReadBin      --> (read file with binned eff.) file with binned efficiency AND [q2 bin indx.(optional)]" << endl;
      cout << "ReadAnaly    --> (read file with analytical eff.) file with analytical efficiency AND [q2 bin indx.(optional)]" << endl;
      
      cout << "ReadGenAnaly --> (read files generated from analytical eff.) file with analytical efficiency AND [q2 bin indx.(optional)]" << endl;

      cout << "FitEff1D     --> file with binned efficiency AND [q2 bin indx.] [thetaL thetaK]" << endl;
      cout << "FitEff2D     --> file with binned efficiency AND [q2 bin indx.]" << endl;
      
      cout << "TestEff      --> file with binned efficiency AND [q2 bin indx.]" << endl;

      return 1;
    }
  
  return 0;
}

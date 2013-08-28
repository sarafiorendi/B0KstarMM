// ########################################################################
// # Program to compute event weights for the B0 --> K*0 mu+ mu- analysis #
// ########################################################################
// # Author: Mauro Dinardo                                                #
// ########################################################################

#ifndef __CINT__
#include <TROOT.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TFile.h>
#include <TLorentzVector.h>
#endif

#include <iostream>
#include <sstream>

#include "Utils.h"
#include "B0KstMuMuTreeContent.h"
#include "B0KstMuMuSingleCandTreeContent.h"

using namespace std;


// ####################
// # Global constants #
// ####################
#define PileUpMCFileName   "PileUp/PileupMCkstJPsi.root" // "PileupMCkstMuMu.root"  OR "PileupMCkstJPsi.root" OR "PileupMCkstPsi2S.root"
#define PileUpDataFileName "PileUp/PileupData_HLT"
#define B0pTFileName       "B0pTDataMC.root"
#define HadppTFileName     "HadppTDataMC.root"
#define HadmpTFileName     "HadmpTDataMC.root"
#define ThetaKFileName     "ThetaKDataMC.root"
#define ParameterFILE      "../python/ParameterFile.txt"

#define SignalType 1 // If checking MC B0 --> K*0 mumu  : 1
                     // If checking MC B0 --> J/psi K*0 : 3
                     // If checking MC B0 --> psi(2S) K*0 : 5


// ####################
// # Global variables #
// ####################
TTree* theTreeIn;
TTree* theTreeOut;
Utils* Utility;
B0KstMuMuSingleCandTreeContent* NTupleIn;
B0KstMuMuSingleCandTreeContent* NTupleOut;


// #######################
// # Function Definition #
// #######################
void AddRecoVariables (string option);
void AddGenVariables (string option);
template<class T> void AddEvWeightPileup (T* NTupleOut);
template<class T> void AddEvWeightB0pT (T* NTupleOut);
template<class T> void AddEvWeightHadpT (T* NTupleOut, string trkSign);
template<class T> void AddEvWeightThetaK (B0KstMuMuSingleCandTreeContent* NTupleOut);


// ###########################
// # Function Implementation #
// ###########################
void AddRecoVariables (string option)
{
  bool B0notB0bar = true;

  NTupleOut->ClearNTuple();
  NTupleOut->MakeTreeBranches(theTreeOut);

  NTupleIn->ClearNTuple();
  NTupleIn->SetBranchAddresses(theTreeIn);
  int nEntries = theTreeIn->GetEntries();
  cout << "\n@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;


  for (int entry = 0; entry < nEntries; entry++)
    {
      theTreeIn->GetEntry(entry);

      for (unsigned int i = 0; i < NTupleIn->bMass->size(); i++)
	{
	  if (((option == "nvTruthMatchReco") && ((NTupleIn->genSignal == SignalType) || (NTupleIn->genSignal == SignalType+1)) && (NTupleIn->truthMatchSignal->at(i) == true)) || (option == "nvAllReco"))
	    {
	      NTupleOut->CopyData(NTupleIn, i);


	      // ########################
	      // # Adding new variables #
	      // ########################
	      if (NTupleIn->genSignal == SignalType)
		{
		  B0notB0bar = true;

		  double cosThetaMup, cosThetaMupErr;
		  double cosThetaK, cosThetaKErr;
		  double phiKstMuMuPlane;
		  // ############
		  // # B0 boost #
		  // ############
		  TLorentzVector LoreVecB0;


		  TLorentzVector LoreVecMuMu;
		  LoreVecMuMu.SetXYZM(NTupleIn->mumuPx->at(i),NTupleIn->mumuPy->at(i),NTupleIn->mumuPz->at(i),NTupleIn->mumuMass->at(i));
		  TVector3 boostMuMu = LoreVecMuMu.BoostVector();
		  TLorentzVector LoreVecMup;
		  LoreVecMup.SetXYZM(NTupleIn->mupPx->at(i),NTupleIn->mupPy->at(i),NTupleIn->mupPz->at(i),Utility->muonMass);
		  // ################################
		  // # Boost mu+ to mumu ref. frame #
		  // ################################
		  LoreVecMup.Boost(-boostMuMu);
		  // ###############################
		  // # Boost B0 to mumu ref. frame #
		  // ###############################
		  LoreVecB0.SetXYZM(NTupleIn->bPx->at(i),NTupleIn->bPy->at(i),NTupleIn->bPz->at(i),NTupleIn->bMass->at(i));
		  LoreVecB0.Boost(-boostMuMu);
		  // ##################
		  // # Compute angles #
		  // ##################
		  Utility->computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
					   LoreVecMup.Px(),LoreVecMup.Py(),LoreVecMup.Pz(),
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   &cosThetaMup,&cosThetaMupErr);


		  TLorentzVector LoreVecKst;
		  LoreVecKst.SetXYZM(NTupleIn->kstPx->at(i),NTupleIn->kstPy->at(i),NTupleIn->kstPz->at(i),NTupleIn->kstMass->at(i));
		  TVector3 boostKst = LoreVecKst.BoostVector();
		  TLorentzVector LoreVecK;
		  LoreVecK.SetXYZM(NTupleIn->kstTrkpPx->at(i),NTupleIn->kstTrkpPy->at(i),NTupleIn->kstTrkpPz->at(i),Utility->kaonMass);
		  // ##############################
		  // # Boost K+ to K*0 ref. frame #
		  // ##############################
		  LoreVecK.Boost(-boostKst);
		  // ##############################
		  // # Boost B0 to K*0 ref. frame #
		  // ##############################
		  LoreVecB0.SetXYZM(NTupleIn->bPx->at(i),NTupleIn->bPy->at(i),NTupleIn->bPz->at(i),NTupleIn->bMass->at(i));
		  LoreVecB0.Boost(-boostKst);
		  // ##################
		  // # Compute angles #
		  // ##################
		  Utility->computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
					   LoreVecK.Px(),LoreVecK.Py(),LoreVecK.Pz(),
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   &cosThetaK,&cosThetaKErr);


		  // ######################################################################
		  // # Angle between [mu+ - mu-] and [K - pi] planes in the B0 ref. frame #
		  // ######################################################################
		  LoreVecB0.SetXYZM(NTupleIn->bPx->at(i),NTupleIn->bPy->at(i),NTupleIn->bPz->at(i),NTupleIn->bBarMass->at(i));
		  TVector3 boostB0 = LoreVecB0.BoostVector();
		  LoreVecMup.SetXYZM(NTupleIn->mupPx->at(i),NTupleIn->mupPy->at(i),NTupleIn->mupPz->at(i),Utility->muonMass);
		  TLorentzVector LoreVecMum;
		  LoreVecMum.SetXYZM(NTupleIn->mumPx->at(i),NTupleIn->mumPy->at(i),NTupleIn->mumPz->at(i),Utility->muonMass);
		  LoreVecK.SetXYZM(NTupleIn->kstTrkpPx->at(i),NTupleIn->kstTrkpPy->at(i),NTupleIn->kstTrkpPz->at(i),Utility->kaonMass);
		  TLorentzVector LoreVecPi;
		  LoreVecPi.SetXYZM(NTupleIn->kstTrkmPx->at(i),NTupleIn->kstTrkmPy->at(i),NTupleIn->kstTrkmPz->at(i),Utility->pionMass);

		  LoreVecMum.Boost(-boostB0);
		  LoreVecMup.Boost(-boostB0);
		  LoreVecK.Boost(-boostB0);
		  LoreVecPi.Boost(-boostB0);
		  TVector3 MuMuPlane = LoreVecMup.Vect().Cross(LoreVecMum.Vect());
		  TVector3 KstPlane  = LoreVecK.Vect().Cross(LoreVecPi.Vect());
		  if (MuMuPlane.Cross(KstPlane).Dot(-LoreVecB0.Vect()) > 0.0) phiKstMuMuPlane = MuMuPlane.Angle(KstPlane);
		  else phiKstMuMuPlane = -MuMuPlane.Angle(KstPlane);


		  NTupleOut->B0MassArb          = NTupleIn->bMass->at(i);
		  NTupleOut->CosThetaMuArb      = cosThetaMup;
		  NTupleOut->CosThetaKArb       = cosThetaK;
		  NTupleOut->PhiKstMuMuPlaneArb = phiKstMuMuPlane;
		}
	      else if (NTupleIn->genSignal == SignalType+1)
		{
		  B0notB0bar = false;
	      
		  double cosThetaMum, cosThetaMumErr;
		  double cosThetaK, cosThetaKErr;
		  double phiKstMuMuPlane;
		  // ############
		  // # B0 boost #
		  // ############
		  TLorentzVector LoreVecB0;
	

		  TLorentzVector LoreVecMuMu;
		  LoreVecMuMu.SetXYZM(NTupleIn->mumuPx->at(i),NTupleIn->mumuPy->at(i),NTupleIn->mumuPz->at(i),NTupleIn->mumuMass->at(i));
		  TVector3 boostMuMu = LoreVecMuMu.BoostVector();
		  TLorentzVector LoreVecMum;
		  LoreVecMum.SetXYZM(NTupleIn->mumPx->at(i),NTupleIn->mumPy->at(i),NTupleIn->mumPz->at(i),Utility->muonMass);
		  // ################################
		  // # Boost mu- to mumu ref. frame #
		  // ################################
		  LoreVecMum.Boost(-boostMuMu);
		  // ###############################
		  // # Boost B0 to mumu ref. frame #
		  // ###############################
		  LoreVecB0.SetXYZM(NTupleIn->bPx->at(i),NTupleIn->bPy->at(i),NTupleIn->bPz->at(i),NTupleIn->bBarMass->at(i));
		  LoreVecB0.Boost(-boostMuMu);
		  // ##################
		  // # Compute angles #
		  // ##################
		  Utility->computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
					   LoreVecMum.Px(),LoreVecMum.Py(),LoreVecMum.Pz(),
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   &cosThetaMum,&cosThetaMumErr);


		  TLorentzVector LoreVecKst;
		  LoreVecKst.SetXYZM(NTupleIn->kstPx->at(i),NTupleIn->kstPy->at(i),NTupleIn->kstPz->at(i),NTupleIn->kstBarMass->at(i));
		  TVector3 boostKst = LoreVecKst.BoostVector();
		  TLorentzVector LoreVecK;
		  LoreVecK.SetXYZM(NTupleIn->kstTrkmPx->at(i),NTupleIn->kstTrkmPy->at(i),NTupleIn->kstTrkmPz->at(i),Utility->kaonMass);
		  // ##############################
		  // # Boost K- to K*0 ref. frame #
		  // ##############################
		  LoreVecK.Boost(-boostKst);
		  // ##############################
		  // # Boost B0 to K*0 ref. frame #
		  // ##############################
		  LoreVecB0.SetXYZM(NTupleIn->bPx->at(i),NTupleIn->bPy->at(i),NTupleIn->bPz->at(i),NTupleIn->bBarMass->at(i));
		  LoreVecB0.Boost(-boostKst);
		  // ##################
		  // # Compute angles #
		  // ##################
		  Utility->computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
					   LoreVecK.Px(),LoreVecK.Py(),LoreVecK.Pz(),
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   &cosThetaK,&cosThetaKErr);


		  // ######################################################################
		  // # Angle between [mu+ - mu-] and [K - pi] planes in the B0 ref. frame #
		  // ######################################################################
		  LoreVecB0.SetXYZM(NTupleIn->bPx->at(i),NTupleIn->bPy->at(i),NTupleIn->bPz->at(i),NTupleIn->bBarMass->at(i));
		  TVector3 boostB0 = LoreVecB0.BoostVector();
		  LoreVecMum.SetXYZM(NTupleIn->mumPx->at(i),NTupleIn->mumPy->at(i),NTupleIn->mumPz->at(i),Utility->muonMass);
		  TLorentzVector LoreVecMup;
		  LoreVecMup.SetXYZM(NTupleIn->mupPx->at(i),NTupleIn->mupPy->at(i),NTupleIn->mupPz->at(i),Utility->muonMass);
		  LoreVecK.SetXYZM(NTupleIn->kstTrkmPx->at(i),NTupleIn->kstTrkmPy->at(i),NTupleIn->kstTrkmPz->at(i),Utility->kaonMass);
		  TLorentzVector LoreVecPi;
		  LoreVecPi.SetXYZM(NTupleIn->kstTrkpPx->at(i),NTupleIn->kstTrkpPy->at(i),NTupleIn->kstTrkpPz->at(i),Utility->pionMass);

		  LoreVecMum.Boost(-boostB0);
		  LoreVecMup.Boost(-boostB0);
		  LoreVecK.Boost(-boostB0);
		  LoreVecPi.Boost(-boostB0);
		  TVector3 MuMuPlane = LoreVecMum.Vect().Cross(LoreVecMup.Vect());
		  TVector3 KstPlane  = LoreVecK.Vect().Cross(LoreVecPi.Vect());
		  if (MuMuPlane.Cross(KstPlane).Dot(-LoreVecB0.Vect()) > 0.0) phiKstMuMuPlane = MuMuPlane.Angle(KstPlane);
		  else phiKstMuMuPlane = -MuMuPlane.Angle(KstPlane);


		  NTupleOut->B0MassArb          = NTupleIn->bBarMass->at(i);
		  NTupleOut->CosThetaMuArb      = cosThetaMum;
		  NTupleOut->CosThetaKArb       = cosThetaK;
		  NTupleOut->PhiKstMuMuPlaneArb = phiKstMuMuPlane;
		}

	      NTupleOut->B0notB0bar = B0notB0bar;
	      NTupleOut->B0pT       = sqrt(NTupleIn->bPx->at(i)*NTupleIn->bPx->at(i) + NTupleIn->bPy->at(i)*NTupleIn->bPy->at(i));
	      NTupleOut->B0Eta      = Utility->computeEta (NTupleIn->bPx->at(i),NTupleIn->bPy->at(i),NTupleIn->bPz->at(i));
	      NTupleOut->B0Phi      = Utility->computePhi (NTupleIn->bPx->at(i),NTupleIn->bPy->at(i),NTupleIn->bPz->at(i));

	      theTreeOut->Fill();
	      NTupleOut->ClearNTuple();

	      if (option == "nvTruthMatchReco") break;
	    }
	}
    }
}


void AddGenVariables (string option)
{
  bool B0notB0bar = true;

  NTupleOut->ClearNTuple();
  NTupleOut->MakeTreeBranches(theTreeOut);

  NTupleIn->ClearNTuple();
  NTupleIn->SetBranchAddresses(theTreeIn);
  int nEntries = theTreeIn->GetEntries();
  cout << "\n@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;


  for (int entry = 0; entry < nEntries; entry++)
    {
      theTreeIn->GetEntry(entry);

      if ((NTupleIn->genSignal == SignalType) || (NTupleIn->genSignal == SignalType+1))
	{
	  if (option == "nvGen") NTupleOut->CopyData(NTupleIn, -1);
	  else if (option == "nvGen2SingleCand")
	    {
	      NTupleOut->CopyData(NTupleIn, 0);
	      NTupleOut->B0MassArb = NTupleIn->B0MassArb;
	    }
	  else 
	    {
	      cout << "[AddVars2Candidates::AddGenVariables]\tWrong option: " << option << endl;
	      exit(1);
	    }


	  // ########################
	  // # Adding new variables #
	  // ########################
	  if (NTupleIn->genSignal == SignalType)
	    {
	      B0notB0bar = true;

	      double cosThetaMup, cosThetaMupErr;
	      double cosThetaK, cosThetaKErr;
	      double phiKstMuMuPlane;
	      // ############
	      // # B0 boost #
	      // ############
	      TLorentzVector LoreVecB0;
	

	      TLorentzVector LoreVecMuMu;
	      LoreVecMuMu.SetXYZM(NTupleIn->genMumPx+NTupleIn->genMupPx,NTupleIn->genMumPy+NTupleIn->genMupPy,NTupleIn->genMumPz+NTupleIn->genMupPz,
				  Utility->computeInvMass(NTupleIn->genMumPx,NTupleIn->genMumPy,NTupleIn->genMumPz,Utility->muonMass,
							  NTupleIn->genMupPx,NTupleIn->genMupPy,NTupleIn->genMupPz,Utility->muonMass));
	      TVector3 boostMuMu = LoreVecMuMu.BoostVector();
	      TLorentzVector LoreVecMup;
	      LoreVecMup.SetXYZM(NTupleIn->genMupPx,NTupleIn->genMupPy,NTupleIn->genMupPz,Utility->muonMass);
	      // ################################
	      // # Boost mu+ to mumu ref. frame #
	      // ################################
	      LoreVecMup.Boost(-boostMuMu);
	      // ###############################
	      // # Boost B0 to mumu ref. frame #
	      // ###############################
	      LoreVecB0.SetXYZM(NTupleIn->genB0Px,NTupleIn->genB0Py,NTupleIn->genB0Pz,NTupleIn->genB0Mass);
	      LoreVecB0.Boost(-boostMuMu);
	      // ##################
	      // # Compute angles #
	      // ##################
	      Utility->computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
				       LoreVecMup.Px(),LoreVecMup.Py(),LoreVecMup.Pz(),
				       0.0,0.0,0.0,
				       0.0,0.0,0.0,
				       0.0,0.0,0.0,
				       0.0,0.0,0.0,
				       &cosThetaMup,&cosThetaMupErr);


	      TLorentzVector LoreVecKst;
	      LoreVecKst.SetXYZM(NTupleIn->genKstPx,NTupleIn->genKstPy,NTupleIn->genKstPz,NTupleIn->genKstMass);
	      TVector3 boostKst = LoreVecKst.BoostVector();
	      TLorentzVector LoreVecK;
	      LoreVecK.SetXYZM(NTupleIn->genKstTrkpPx,NTupleIn->genKstTrkpPy,NTupleIn->genKstTrkpPz,Utility->kaonMass);
	      // ##############################
	      // # Boost K+ to K*0 ref. frame #
	      // ##############################
	      LoreVecK.Boost(-boostKst);
	      // ##############################
	      // # Boost B0 to K*0 ref. frame #
	      // ##############################
	      LoreVecB0.SetXYZM(NTupleIn->genB0Px,NTupleIn->genB0Py,NTupleIn->genB0Pz,NTupleIn->genB0Mass);
	      LoreVecB0.Boost(-boostKst);
	      // ##################
	      // # Compute angles #
	      // ##################
	      Utility->computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
				       LoreVecK.Px(),LoreVecK.Py(),LoreVecK.Pz(),
				       0.0,0.0,0.0,
				       0.0,0.0,0.0,
				       0.0,0.0,0.0,
				       0.0,0.0,0.0,
				       &cosThetaK,&cosThetaKErr);


	      // ######################################################################
	      // # Angle between [mu+ - mu-] and [K - pi] planes in the B0 ref. frame #
	      // ######################################################################
	      LoreVecB0.SetXYZM(NTupleIn->genB0Px,NTupleIn->genB0Py,NTupleIn->genB0Pz,NTupleIn->genB0Mass);
	      TVector3 boostB0 = LoreVecB0.BoostVector();
	      LoreVecMup.SetXYZM(NTupleIn->genMupPx,NTupleIn->genMupPy,NTupleIn->genMupPz,Utility->muonMass);
	      TLorentzVector LoreVecMum;
	      LoreVecMum.SetXYZM(NTupleIn->genMumPx,NTupleIn->genMumPy,NTupleIn->genMumPz,Utility->muonMass);
	      LoreVecK.SetXYZM(NTupleIn->genKstTrkpPx,NTupleIn->genKstTrkpPy,NTupleIn->genKstTrkpPz,Utility->kaonMass);
	      TLorentzVector LoreVecPi;
	      LoreVecPi.SetXYZM(NTupleIn->genKstTrkmPx,NTupleIn->genKstTrkmPy,NTupleIn->genKstTrkmPz,Utility->pionMass);

	      LoreVecMum.Boost(-boostB0);
	      LoreVecMup.Boost(-boostB0);
	      LoreVecK.Boost(-boostB0);
	      LoreVecPi.Boost(-boostB0);
	      TVector3 MuMuPlane = LoreVecMup.Vect().Cross(LoreVecMum.Vect());
	      TVector3 KstPlane  = LoreVecK.Vect().Cross(LoreVecPi.Vect());
	      if (MuMuPlane.Cross(KstPlane).Dot(-LoreVecB0.Vect()) > 0.0) phiKstMuMuPlane = MuMuPlane.Angle(KstPlane);
	      else phiKstMuMuPlane = -MuMuPlane.Angle(KstPlane);


	      NTupleOut->CosThetaMuArb      = cosThetaMup;
	      NTupleOut->CosThetaKArb       = cosThetaK;
	      NTupleOut->PhiKstMuMuPlaneArb = phiKstMuMuPlane;
	    }
	  else if (NTupleIn->genSignal == SignalType+1)
	    {
	      B0notB0bar = false;
	      
	      double cosThetaMum, cosThetaMumErr;
	      double cosThetaK, cosThetaKErr;
	      double phiKstMuMuPlane;
	      // ############
	      // # B0 boost #
	      // ############
	      TLorentzVector LoreVecB0;
	

	      TLorentzVector LoreVecMuMu;
	      LoreVecMuMu.SetXYZM(NTupleIn->genMumPx+NTupleIn->genMupPx,NTupleIn->genMumPy+NTupleIn->genMupPy,NTupleIn->genMumPz+NTupleIn->genMupPz,
				  Utility->computeInvMass(NTupleIn->genMumPx,NTupleIn->genMumPy,NTupleIn->genMumPz,Utility->muonMass,
							  NTupleIn->genMupPx,NTupleIn->genMupPy,NTupleIn->genMupPz,Utility->muonMass));
	      TVector3 boostMuMu = LoreVecMuMu.BoostVector();
	      TLorentzVector LoreVecMum;
	      LoreVecMum.SetXYZM(NTupleIn->genMumPx,NTupleIn->genMumPy,NTupleIn->genMumPz,Utility->muonMass);
	      // ################################
	      // # Boost mu- to mumu ref. frame #
	      // ################################
	      LoreVecMum.Boost(-boostMuMu);
	      // ###############################
	      // # Boost B0 to mumu ref. frame #
	      // ###############################
	      LoreVecB0.SetXYZM(NTupleIn->genB0Px,NTupleIn->genB0Py,NTupleIn->genB0Pz,NTupleIn->genB0Mass);
	      LoreVecB0.Boost(-boostMuMu);
	      // ##################
	      // # Compute angles #
	      // ##################
	      Utility->computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
				       LoreVecMum.Px(),LoreVecMum.Py(),LoreVecMum.Pz(),
				       0.0,0.0,0.0,
				       0.0,0.0,0.0,
				       0.0,0.0,0.0,
				       0.0,0.0,0.0,
				       &cosThetaMum,&cosThetaMumErr);


	      TLorentzVector LoreVecKst;
	      LoreVecKst.SetXYZM(NTupleIn->genKstPx,NTupleIn->genKstPy,NTupleIn->genKstPz,NTupleIn->genKstMass);
	      TVector3 boostKst = LoreVecKst.BoostVector();
	      TLorentzVector LoreVecK;
	      LoreVecK.SetXYZM(NTupleIn->genKstTrkmPx,NTupleIn->genKstTrkmPy,NTupleIn->genKstTrkmPz,Utility->kaonMass);
	      // ##############################
	      // # Boost K- to K*0 ref. frame #
	      // ##############################
	      LoreVecK.Boost(-boostKst);
	      // ##############################
	      // # Boost B0 to K*0 ref. frame #
	      // ##############################
	      LoreVecB0.SetXYZM(NTupleIn->genB0Px,NTupleIn->genB0Py,NTupleIn->genB0Pz,NTupleIn->genB0Mass);
	      LoreVecB0.Boost(-boostKst);
	      // ##################
	      // # Compute angles #
	      // ##################
	      Utility->computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
				       LoreVecK.Px(),LoreVecK.Py(),LoreVecK.Pz(),
				       0.0,0.0,0.0,
				       0.0,0.0,0.0,
				       0.0,0.0,0.0,
				       0.0,0.0,0.0,
				       &cosThetaK,&cosThetaKErr);


	      // ######################################################################
	      // # Angle between [mu+ - mu-] and [K - pi] planes in the B0 ref. frame #
	      // ######################################################################
	      LoreVecB0.SetXYZM(NTupleIn->genB0Px,NTupleIn->genB0Py,NTupleIn->genB0Pz,NTupleIn->genB0Mass);
	      TVector3 boostB0 = LoreVecB0.BoostVector();
	      LoreVecMum.SetXYZM(NTupleIn->genMumPx,NTupleIn->genMumPy,NTupleIn->genMumPz,Utility->muonMass);
	      TLorentzVector LoreVecMup;
	      LoreVecMup.SetXYZM(NTupleIn->genMupPx,NTupleIn->genMupPy,NTupleIn->genMupPz,Utility->muonMass);
	      LoreVecK.SetXYZM(NTupleIn->genKstTrkmPx,NTupleIn->genKstTrkmPy,NTupleIn->genKstTrkmPz,Utility->kaonMass);
	      TLorentzVector LoreVecPi;
	      LoreVecPi.SetXYZM(NTupleIn->genKstTrkpPx,NTupleIn->genKstTrkpPy,NTupleIn->genKstTrkpPz,Utility->pionMass);

	      LoreVecMum.Boost(-boostB0);
	      LoreVecMup.Boost(-boostB0);
	      LoreVecK.Boost(-boostB0);
	      LoreVecPi.Boost(-boostB0);
	      TVector3 MuMuPlane = LoreVecMum.Vect().Cross(LoreVecMup.Vect());
	      TVector3 KstPlane  = LoreVecK.Vect().Cross(LoreVecPi.Vect());
	      if (MuMuPlane.Cross(KstPlane).Dot(-LoreVecB0.Vect()) > 0.0) phiKstMuMuPlane = MuMuPlane.Angle(KstPlane);
	      else phiKstMuMuPlane = -MuMuPlane.Angle(KstPlane);


	      NTupleOut->CosThetaMuArb      = cosThetaMum;
	      NTupleOut->CosThetaKArb       = cosThetaK;
	      NTupleOut->PhiKstMuMuPlaneArb = phiKstMuMuPlane;
	    }

	  NTupleOut->B0notB0bar = B0notB0bar;
	  NTupleOut->B0pT       = sqrt(NTupleIn->genB0Px*NTupleIn->genB0Px + NTupleIn->genB0Py*NTupleIn->genB0Py);
	  NTupleOut->B0Eta      = Utility->computeEta (NTupleIn->genB0Px,NTupleIn->genB0Py,NTupleIn->genB0Pz);
	  NTupleOut->B0Phi      = Utility->computePhi (NTupleIn->genB0Px,NTupleIn->genB0Py,NTupleIn->genB0Pz);

	  NTupleOut->mumuMass->push_back(Utility->computeInvMass(NTupleIn->genMumPx,NTupleIn->genMumPy,NTupleIn->genMumPz,Utility->muonMass,
								 NTupleIn->genMupPx,NTupleIn->genMupPy,NTupleIn->genMupPz,Utility->muonMass));
	  NTupleOut->mumuMassE->push_back(0.0);

	  NTupleOut->kstTrkmPx->push_back(NTupleIn->genKstTrkmPx);
	  NTupleOut->kstTrkmPy->push_back(NTupleIn->genKstTrkmPy);
	  NTupleOut->kstTrkmPz->push_back(NTupleIn->genKstTrkmPz);

	  NTupleOut->kstTrkpPx->push_back(NTupleIn->genKstTrkpPx);
	  NTupleOut->kstTrkpPy->push_back(NTupleIn->genKstTrkpPy);
	  NTupleOut->kstTrkpPz->push_back(NTupleIn->genKstTrkpPz);

	  NTupleOut->truthMatchSignal->push_back(true);

	  NTupleOut->FillWithNull(NTupleOut->mumuMass->size());
	  theTreeOut->Fill();
	  NTupleOut->ClearNTuple();
	}
    }
}


template<class T> void AddEvWeightPileup (T* NTupleOut)
{
  stringstream myString;
  double value = 0.0;
  int HLTpathIndx;

  NTupleOut->ClearNTuple();
  NTupleOut->MakeTreeBranches(theTreeOut);

  NTupleIn->ClearNTuple();
  NTupleIn->SetBranchAddresses(theTreeIn);
  int nEntries = theTreeIn->GetEntries();
  cout << "\n@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;


  TFile* fileMC = TFile::Open(PileUpMCFileName,"READ");
  TH1D* hMC = (TH1D*)fileMC->Get("pileup");
  hMC->Scale(1. / hMC->Integral());

  vector<TFile*> fileData;
  vector<TH1D*> hData;
  vector<TH1D*> hW;
  for (unsigned int i = 0; i < Utility->GetNHLTCat(); i++)
    {
      myString.str("");
      myString << PileUpDataFileName << i+1 << ".root";
      fileData.push_back(TFile::Open(myString.str().c_str(),"READ"));
      hData.push_back((TH1D*)fileData[i]->Get("pileup"));
      hData.back()->Scale(1. / hData.back()->Integral());
  
      hW.push_back((TH1D*)hData.back()->Clone());
      hW.back()->SetName("PileupWeights");
      hW.back()->Divide(hMC);

      cout << "\n@@@ Computing weights for trigger category: " << i+1 << " (total weight = " << hW.back()->Integral() << ") @@@" << endl;
      for (int j = 0; j < hW.back()->GetNbinsX(); j++) cout << "--> bin " << j+1 << " has weight: " << hW.back()->GetBinContent(j+1) << endl;
    }


  cout << "\n@@@ Assigning the weights to the events @@@" << endl;
  for (int entry = 0; entry < nEntries; entry++)
    {
      theTreeIn->GetEntry(entry);      
      NTupleOut->CopyWholeNTuple(NTupleIn);

      if (typeid(T) == typeid(B0KstMuMuSingleCandTreeContent)) HLTpathIndx = ((B0KstMuMuSingleCandTreeContent*)NTupleOut)->TrigCat;
      else                                                     HLTpathIndx = Utility->HLTpathForEvFraction(((double)entry)/((double)nEntries));

      for (unsigned int i = 0; i < NTupleOut->bunchXingMC->size(); i++)
	{
	  value = NTupleOut->numInteractionsMC->at(i);
	  if ((NTupleOut->bunchXingMC->at(i) == 0) && (hW[HLTpathIndx-1]->FindBin(value) != 0) && (hW[HLTpathIndx-1]->FindBin(value) != hW[HLTpathIndx-1]->GetNbinsX()+1))
	    {
	      NTupleOut->evWeightE2 = NTupleOut->evWeightE2 * pow(hW[HLTpathIndx-1]->GetBinContent(hW[HLTpathIndx-1]->FindBin(value)),2.0) +
		pow(hW[HLTpathIndx-1]->GetBinError(hW[HLTpathIndx-1]->FindBin(value)) * NTupleOut->evWeight,2.0);
	      NTupleOut->evWeight = NTupleOut->evWeight * hW[HLTpathIndx-1]->GetBinContent(hW[HLTpathIndx-1]->FindBin(value));
	      break;
	    }
	}

      theTreeOut->Fill();
      NTupleOut->ClearNTuple();
    }


  for (unsigned int i = 0; i < Utility->GetNHLTCat(); i++) fileData[i]->Close();
  fileMC->Close();
}


template<class T> void AddEvWeightB0pT (T* NTupleOut)
{
  stringstream myString;
  double value = 0.0;
  int nEntries;


  NTupleOut->ClearNTuple();
  NTupleOut->MakeTreeBranches(theTreeOut);

  NTupleIn->ClearNTuple();
  NTupleIn->SetBranchAddresses(theTreeIn);
  nEntries = theTreeIn->GetEntries();
  cout << "\n@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;


  TFile* fileDataMC = TFile::Open(B0pTFileName,"READ");
  TCanvas* cTmp     = (TCanvas*)fileDataMC->Get("c0");
  TH1D* hMC         = (TH1D*)cTmp->GetPrimitive("hM1D");
  TH1D* hData       = (TH1D*)cTmp->GetPrimitive("hDsig1D");
 
  TH1D* hW = (TH1D*)hData->Clone();
  hW->SetName("B0pTWeights");
  hW->Divide(hMC);


  cout << "\n@@@ Assigning the weights to the events @@@" << endl;
  for (int entry = 0; entry < nEntries; entry++)
    {
      theTreeIn->GetEntry(entry);      
      NTupleOut->CopyWholeNTuple(NTupleIn);

      value = NTupleIn->B0pT;
      if ((hW->FindBin(value) != 0) && (hW->FindBin(value) != hW->GetNbinsX()+1))
	{
	  NTupleOut->evWeightE2 = NTupleOut->evWeightE2 * pow(hW->GetBinContent(hW->FindBin(value)),2.0) +
	    pow(hW->GetBinError(hW->FindBin(value)) * NTupleOut->evWeight,2.0);
	  NTupleOut->evWeight = NTupleOut->evWeight * hW->GetBinContent(hW->FindBin(value));
	}

      theTreeOut->Fill();
      NTupleOut->ClearNTuple();
    }
  

  fileDataMC->Close();
}


template<class T> void AddEvWeightHadpT (T* NTupleOut, string trkSign)
// ###################
// # trkSign = "pos" #
// # trkSign = "neg" #
// ###################
{
  stringstream myString;
  double value = 0.0;
  int nEntries;


  NTupleOut->ClearNTuple();
  NTupleOut->MakeTreeBranches(theTreeOut);

  NTupleIn->ClearNTuple();
  NTupleIn->SetBranchAddresses(theTreeIn);
  nEntries = theTreeIn->GetEntries();
  cout << "\n@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;


  TFile* fileDataMC;
  if      (trkSign == "pos") fileDataMC = TFile::Open(HadppTFileName,"READ");
  else if (trkSign == "neg") fileDataMC = TFile::Open(HadmpTFileName,"READ");
  else
    {	
      cout << "[AddVars2Candidates::AddEvWeightHadpT]\tWrong option: " << trkSign << endl;
      exit(1);
    }
  TCanvas* cTmp = (TCanvas*)fileDataMC->Get("c0");
  TH1D* hMC     = (TH1D*)cTmp->GetPrimitive("hM1D");
  TH1D* hData   = (TH1D*)cTmp->GetPrimitive("hDsig1D");
  
  TH1D* hW = (TH1D*)hData->Clone();
  hW->SetName("HadpTWeights");
  hW->Divide(hMC);


  cout << "\n@@@ Assigning the weights to the events @@@" << endl;
  for (int entry = 0; entry < nEntries; entry++)
    {
      theTreeIn->GetEntry(entry);      
      NTupleOut->CopyWholeNTuple(NTupleIn);
	  
      if      (trkSign == "pos") value = sqrt(NTupleOut->kstTrkpPx->at(0)*NTupleOut->kstTrkpPx->at(0) + NTupleOut->kstTrkpPy->at(0)*NTupleOut->kstTrkpPy->at(0));
      else if (trkSign == "neg") value = sqrt(NTupleOut->kstTrkmPx->at(0)*NTupleOut->kstTrkmPx->at(0) + NTupleOut->kstTrkmPy->at(0)*NTupleOut->kstTrkmPy->at(0));
      
      if ((hW->FindBin(value) != 0) && (hW->FindBin(value) != hW->GetNbinsX()+1))
	{
	  NTupleOut->evWeightE2 = NTupleOut->evWeightE2 * pow(hW->GetBinContent(hW->FindBin(value)),2.0) +
	    pow(hW->GetBinError(hW->FindBin(value)) * NTupleOut->evWeight,2.0);
	  NTupleOut->evWeight = NTupleOut->evWeight * hW->GetBinContent(hW->FindBin(value));
	}

      theTreeOut->Fill();
      NTupleOut->ClearNTuple();
    }

  
  fileDataMC->Close();
}


template<class T> void AddEvWeightThetaK (B0KstMuMuSingleCandTreeContent* NTupleOut)
{
  stringstream myString;
  double value = 0.0;
  int nEntries;


  NTupleOut->ClearNTuple();
  NTupleOut->MakeTreeBranches(theTreeOut);

  NTupleIn->ClearNTuple();
  NTupleIn->SetBranchAddresses(theTreeIn);
  nEntries = theTreeIn->GetEntries();
  cout << "\n@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;


  TFile* fileDataMC = TFile::Open(ThetaKFileName,"READ");
  TCanvas* cTmp     = (TCanvas*)fileDataMC->Get("c0");
  TH1D* hMC         = (TH1D*)cTmp->GetPrimitive("hM1D");
  TH1D* hData       = (TH1D*)cTmp->GetPrimitive("hDsig1D");
 
  TH1D* hW = (TH1D*)hData->Clone();
  hW->SetName("ThetaKWeights");
  hW->Divide(hMC);


  cout << "\n@@@ Assigning the weights to the events @@@" << endl;
  for (int entry = 0; entry < nEntries; entry++)
    {
      theTreeIn->GetEntry(entry);      
      NTupleOut->CopyWholeNTuple(NTupleIn);

      value = NTupleOut->CosThetaKArb;
      if ((hW->FindBin(value) != 0) && (hW->FindBin(value) != hW->GetNbinsX()+1))
	{
	  NTupleOut->evWeightE2 = NTupleOut->evWeightE2 * pow(hW->GetBinContent(hW->FindBin(value)),2.0) +
	    pow(hW->GetBinError(hW->FindBin(value)) * NTupleOut->evWeight,2.0);
	  NTupleOut->evWeight = NTupleOut->evWeight * hW->GetBinContent(hW->FindBin(value));
	}

      theTreeOut->Fill();
      NTupleOut->ClearNTuple();
    }
  

  fileDataMC->Close();
}


int main (int argc, char** argv)
{
  if (argc >= 4)
    {
      string option = argv[1];
      if ((((option == "pileupW") || (option == "HadpTW")) && (argc == 5)) ||
	  (((option == "B0pTW") || (option == "thetaKW") || (option == "nvAllReco") || (option == "nvTruthMatchReco") || (option == "nvGen") || (option == "nvGen2SingleCand")) && (argc == 4)))
	{
	  string fileNameIn  = argv[2];
	  string fileNameOut = argv[3];
	  string fileType = "";
	  if ((option == "pileupW") || (option == "HadpTW")) fileType = argv[4];

	  Utility = new Utils();
	  Utility->ReadTriggerPathsANDCutsANDEntries(ParameterFILE);

	  TFile* NtplFileIn = new TFile(fileNameIn.c_str(), "READ");
	  theTreeIn = (TTree*) NtplFileIn->Get("B0KstMuMu/B0KstMuMuNTuple");
	  NTupleIn = new B0KstMuMuSingleCandTreeContent();
	  NTupleIn->Init();

	  TFile* NtplFileOut = new TFile(fileNameOut.c_str(), "RECREATE");
	  NTupleOut = NULL;
	  NtplFileOut->mkdir("B0KstMuMu");
	  NtplFileOut->cd("B0KstMuMu");
	  theTreeOut = new TTree("B0KstMuMuNTuple","B0KstMuMuNTuple");
	  NTupleOut = new B0KstMuMuSingleCandTreeContent();
	  NTupleOut->Init();


	  cout << "\n@@@ Settings @@@" << endl;
	  cout << "PileUp MC file name: " << PileUpMCFileName << endl;
	  cout << "PileUp data file name: " << PileUpDataFileName << endl;
	  cout << "B0 pT data-MC file name: " << B0pTFileName << endl;
	  cout << "Positive hadron pT data-MC file name: " << HadppTFileName << endl;
	  cout << "Negative hadron pT data-MC file name: " << HadmpTFileName << endl;
	  cout << "cos(theta_K) data-MC file name: " << ThetaKFileName << endl;
	  cout << "ParameterFILE: " << ParameterFILE << endl;
	  cout << "Signal Type: " << SignalType << endl;


	  if (option == "pileupW")
	    {
	      if (fileType == "single") AddEvWeightPileup<B0KstMuMuSingleCandTreeContent>(NTupleOut);
	      else                      AddEvWeightPileup<B0KstMuMuTreeContent>(NTupleOut);
	      cout << "\n@@@ Added new event weight from pileup @@@" << endl;
	    }
	  else if (option == "B0pTW")
	    {
	      AddEvWeightB0pT<B0KstMuMuSingleCandTreeContent>(NTupleOut);
	      cout << "\n@@@ Added new event weight from B0 pT @@@" << endl;
	    }
	  else if (option == "HadpTW")
	    {
	      AddEvWeightHadpT<B0KstMuMuSingleCandTreeContent>(NTupleOut,fileType);
	      cout << "\n@@@ Added new event weight from hadron pT @@@" << endl;
	    }
	  else if (option == "thetaKW")
	    {
	      AddEvWeightThetaK<B0KstMuMuSingleCandTreeContent>(NTupleOut);
	      cout << "\n@@@ Added new event weight from cos(theta_K) @@@" << endl;
	    }
	  else if ((option == "nvAllReco") || (option == "nvTruthMatchReco") || (option == "nvGen") || (option == "nvGen2SingleCand"))
	    {
	      if ((option == "nvAllReco") || (option == "nvTruthMatchReco"))
		{
		  AddRecoVariables(option);
		  cout << "\n@@@ Added new variables to reco-candidates @@@" << endl;
		}
	      else if ((option == "nvGen") || (option == "nvGen2SingleCand"))
		{
		  AddGenVariables(option);
		  cout << "\n@@@ Added new variables to gen-events @@@" << endl;
		}
	    }


	  NtplFileOut->cd("B0KstMuMu");
	  theTreeOut->Write();
	  NtplFileOut->Close();
	  NTupleOut->Destroy();
	  delete NTupleOut;

	  NtplFileIn->Close();
	  NTupleIn->Destroy();
	  delete Utility;
	  delete NTupleIn;

	  return 0;
	}
      else
	{
	  cout << "Parameter missing: " << endl;
	  cout << "./AddVars2Candidates [pileupW B0pTW HadpTW thetaKW nvAllReco nvTruthMatchReco nvGen nvGen2SingleCand] inputFile.root outputFile.root [[if pileupW]outputFile-type(single/multi)] [[if HadpT]pos/neg]" << endl;
	  cout << "- pileupW          : change the weight to all single/multiple candidates according to pileup weight" << endl;
	  cout << "- B0pTW            : change the weight to all single candidates according to B0 pT weight" << endl;
	  cout << "- HadpTW           : change the weight to all single candidates according to hadron pT weight" << endl;
	  cout << "- thetaKW          : change the weight to all single candidates according to cos(theta_K) weight" << endl;
	  cout << "- nvAllReco        : generate new NTupleOut giveing to each candidate in NTupleIn a new TTree entry and adding the RECO single candidate variables" << endl;
	  cout << "- nvTruthMatchReco : generate new NTupleOut adding the RECO single candidate variables only for the first truth matched candidate in a signal event, discarding the other candidates" << endl;
	  cout << "- nvGen            : generate new NTupleOut adding the GEN single candidate variables to each gen-event to an NTuple computed from GEN-MC (NoFilter,Filter,Multi-candidates)" << endl;
	  cout << "- nvGen2SingleCand : generate new NTupleOut adding the GEN single candidate variables to each gen-event to an NTuple computed from a GEN-RECO-MC (Single-candidate)" << endl;

	  return 1;
	}
    }
  else
    {
      cout << "Parameter missing: " << endl;
      cout << "./AddVars2Candidates [pileupW B0pTW HadpTW thetaKW nvAllReco nvTruthMatchReco nvGen nvGen2SingleCand] inputFile.root outputFile.root [[if pileupW]outputFile-type(single/multi)] [[if HadpT]pos/neg]" << endl;
      cout << "- pileupW          : change the weight to all single/multiple candidates according to pileup weight" << endl;
      cout << "- B0pTW            : change the weight to all single candidates according to B0 pT weight" << endl;
      cout << "- HadpTW           : change the weight to all single candidates according to hadron pT weight" << endl;
      cout << "- thetaKW          : change the weight to all single candidates according to cos(theta_K) weight" << endl;
      cout << "- nvAllReco        : generate new NTupleOut giveing to each candidate in NTupleIn a new TTree entry and adding the RECO single candidate variables" << endl;
      cout << "- nvTruthMatchReco : generate new NTupleOut adding the RECO single candidate variables only for the first truth matched candidate in a signal event, discarding the other candidates" << endl;
      cout << "- nvGen            : generate new NTupleOut adding the GEN single candidate variables to each gen-event to an NTuple computed from GEN-MC (NoFilter,Filter,Multi-candidates)" << endl;
      cout << "- nvGen2SingleCand : generate new NTupleOut adding the GEN single candidate variables to each gen-event to an NTuple computed from a GEN-RECO-MC (Single-candidate)" << endl;
      
      return 1;
    }
  
  return 0;
}

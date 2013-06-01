// ##########################################################################
// # Program to select one single candidate the B0 --> K*0 mu+ mu- analysis #
// ##########################################################################
// # Author: Mauro Dinardo                                                  #
// ##########################################################################

#ifndef __CINT__
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#endif

#include <math.h>
#include <iostream>
#include <sstream>

#include "Utils.h"
#include "B0KstMuMuTreeContent.h"
#include "B0KstMuMuSingleCandTreeContent.h"

using namespace std;


// ####################
// # Global constants #
// ####################
#define SignalType 1 // If checking MC B0 --> K*0 mumu  : 1
                     // If checking MC B0 --> J/psi K*0 : 3
                     // If checking MC B0 --> psi(2S) K*0 : 5
#define DoTrigCheck 1
#define SaveHistos false
#define DoMCTruth false    // Compute the signel candidate variables from MC-truth values
#define TagFromTruth false // Assign the CP-eigenstate from MC-truth
#define ParameterFILE "../python/ParameterFile.txt"


// ####################
// # Global variables #
// ####################
Utils* Utility;
TTree* theTreeIn;
TTree* theTreeOut;
B0KstMuMuTreeContent* NTupleIn;
B0KstMuMuSingleCandTreeContent* NTupleOut;


void SelectBestCand ()
{
  unsigned int countCands;
  bool B0notB0bar;
  int BestCandIndx;
  int TrigCat;

  NTupleOut->ClearNTuple();
  NTupleOut->MakeTreeBranches(theTreeOut);

  NTupleIn->ClearNTuple();
  NTupleIn->SetBranchAddresses(theTreeIn);
  int nEntries = theTreeIn->GetEntries();
  cout << "\n@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;


  // ###########################
  // # Loop over B0 Candidates #
  // ###########################
  for (int entry = 0; entry < nEntries; entry++)
    {
      theTreeIn->GetEntry(entry);

      if (Utility->ChooseBestCand(NTupleIn, DoTrigCheck, ((double)entry)/((double)nEntries), &BestCandIndx, &B0notB0bar, &TrigCat, &countCands) == true)
	{
	  NTupleOut->CopyData(NTupleIn, BestCandIndx);
	  
	  // #############################################
	  // # Adding new variables for angular analysis #
	  // #############################################
	  if (((TagFromTruth == false) && (B0notB0bar == true)) ||
	      ((TagFromTruth == true) && (NTupleIn->genSignal == SignalType)))
	    {
	      double cosThetaMup, cosThetaMupErr;
	      double cosThetaK, cosThetaKErr;
	      double phiKstMuMuPlane;
	      if (DoMCTruth == false)
		{
		  // ############
		  // # B0 boost #
		  // ############
		  TLorentzVector LoreVecB0;


		  TLorentzVector LoreVecMuMu;
		  LoreVecMuMu.SetXYZM(NTupleIn->mumuPx->at(BestCandIndx),NTupleIn->mumuPy->at(BestCandIndx),NTupleIn->mumuPz->at(BestCandIndx),NTupleIn->mumuMass->at(BestCandIndx));
		  TVector3 boostMuMu = LoreVecMuMu.BoostVector();
		  TLorentzVector LoreVecMup;
		  LoreVecMup.SetXYZM(NTupleIn->mupPx->at(BestCandIndx),NTupleIn->mupPy->at(BestCandIndx),NTupleIn->mupPz->at(BestCandIndx),Utility->muonMass);
		  // ################################
		  // # Boost mu+ to mumu ref. frame #
		  // ################################
		  LoreVecMup.Boost(-boostMuMu);
		  // ###############################
		  // # Boost B0 to mumu ref. frame #
		  // ###############################
		  LoreVecB0.SetXYZM(NTupleIn->bPx->at(BestCandIndx),NTupleIn->bPy->at(BestCandIndx),NTupleIn->bPz->at(BestCandIndx),NTupleIn->bMass->at(BestCandIndx));
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
		  LoreVecKst.SetXYZM(NTupleIn->kstPx->at(BestCandIndx),NTupleIn->kstPy->at(BestCandIndx),NTupleIn->kstPz->at(BestCandIndx),NTupleIn->kstMass->at(BestCandIndx));
		  TVector3 boostKst = LoreVecKst.BoostVector();
		  TLorentzVector LoreVecK;
		  LoreVecK.SetXYZM(NTupleIn->kstTrkpPx->at(BestCandIndx),NTupleIn->kstTrkpPy->at(BestCandIndx),NTupleIn->kstTrkpPz->at(BestCandIndx),Utility->kaonMass);
		  // ##############################
		  // # Boost K+ to K*0 ref. frame #
		  // ##############################
		  LoreVecK.Boost(-boostKst);
		  // ##############################
		  // # Boost B0 to K*0 ref. frame #
		  // ##############################
		  LoreVecB0.SetXYZM(NTupleIn->bPx->at(BestCandIndx),NTupleIn->bPy->at(BestCandIndx),NTupleIn->bPz->at(BestCandIndx),NTupleIn->bMass->at(BestCandIndx));
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
		  LoreVecB0.SetXYZM(NTupleIn->bPx->at(BestCandIndx),NTupleIn->bPy->at(BestCandIndx),NTupleIn->bPz->at(BestCandIndx),NTupleIn->bBarMass->at(BestCandIndx));
		  TVector3 boostB0 = LoreVecB0.BoostVector();
		  LoreVecMup.SetXYZM(NTupleIn->mupPx->at(BestCandIndx),NTupleIn->mupPy->at(BestCandIndx),NTupleIn->mupPz->at(BestCandIndx),Utility->muonMass);
		  TLorentzVector LoreVecMum;
		  LoreVecMum.SetXYZM(NTupleIn->mumPx->at(BestCandIndx),NTupleIn->mumPy->at(BestCandIndx),NTupleIn->mumPz->at(BestCandIndx),Utility->muonMass);
		  LoreVecK.SetXYZM(NTupleIn->kstTrkpPx->at(BestCandIndx),NTupleIn->kstTrkpPy->at(BestCandIndx),NTupleIn->kstTrkpPz->at(BestCandIndx),Utility->kaonMass);
		  TLorentzVector LoreVecPi;
		  LoreVecPi.SetXYZM(NTupleIn->kstTrkmPx->at(BestCandIndx),NTupleIn->kstTrkmPy->at(BestCandIndx),NTupleIn->kstTrkmPz->at(BestCandIndx),Utility->pionMass);

		  LoreVecMum.Boost(-boostB0);
		  LoreVecMup.Boost(-boostB0);
		  LoreVecK.Boost(-boostB0);
		  LoreVecPi.Boost(-boostB0);
		  TVector3 MuMuPlane = LoreVecMup.Vect().Cross(LoreVecMum.Vect());
		  TVector3 KstPlane  = LoreVecK.Vect().Cross(LoreVecPi.Vect());
		  if (MuMuPlane.Cross(KstPlane).Dot(-LoreVecB0.Vect()) > 0.0) phiKstMuMuPlane = MuMuPlane.Angle(KstPlane);
		  else phiKstMuMuPlane = -MuMuPlane.Angle(KstPlane);
 		}
	      else if (DoMCTruth == true)
		{
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
		}

	      NTupleOut->B0MassArb          = NTupleIn->bMass->at(BestCandIndx);
	      NTupleOut->CosThetaMuArb      = cosThetaMup;
	      NTupleOut->CosThetaKArb       = cosThetaK;
	      NTupleOut->PhiKstMuMuPlaneArb = phiKstMuMuPlane;
	    }
	  else if (((TagFromTruth == false) && (B0notB0bar == false)) ||
		   ((TagFromTruth == true) && (NTupleIn->genSignal == SignalType+1)))
	    {
	      double cosThetaMum, cosThetaMumErr;
	      double cosThetaK, cosThetaKErr;
	      double phiKstMuMuPlane;
	      if (DoMCTruth == false)
		{
		  // ############
		  // # B0 boost #
		  // ############
		  TLorentzVector LoreVecB0;
	

		  TLorentzVector LoreVecMuMu;
		  LoreVecMuMu.SetXYZM(NTupleIn->mumuPx->at(BestCandIndx),NTupleIn->mumuPy->at(BestCandIndx),NTupleIn->mumuPz->at(BestCandIndx),NTupleIn->mumuMass->at(BestCandIndx));
		  TVector3 boostMuMu = LoreVecMuMu.BoostVector();
		  TLorentzVector LoreVecMum;
		  LoreVecMum.SetXYZM(NTupleIn->mumPx->at(BestCandIndx),NTupleIn->mumPy->at(BestCandIndx),NTupleIn->mumPz->at(BestCandIndx),Utility->muonMass);
		  // ################################
		  // # Boost mu- to mumu ref. frame #
		  // ################################
		  LoreVecMum.Boost(-boostMuMu);
		  // ###############################
		  // # Boost B0 to mumu ref. frame #
		  // ###############################
		  LoreVecB0.SetXYZM(NTupleIn->bPx->at(BestCandIndx),NTupleIn->bPy->at(BestCandIndx),NTupleIn->bPz->at(BestCandIndx),NTupleIn->bBarMass->at(BestCandIndx));
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
		  LoreVecKst.SetXYZM(NTupleIn->kstPx->at(BestCandIndx),NTupleIn->kstPy->at(BestCandIndx),NTupleIn->kstPz->at(BestCandIndx),NTupleIn->kstBarMass->at(BestCandIndx));
		  TVector3 boostKst = LoreVecKst.BoostVector();
		  TLorentzVector LoreVecK;
		  LoreVecK.SetXYZM(NTupleIn->kstTrkmPx->at(BestCandIndx),NTupleIn->kstTrkmPy->at(BestCandIndx),NTupleIn->kstTrkmPz->at(BestCandIndx),Utility->kaonMass);
		  // ##############################
		  // # Boost K- to K*0 ref. frame #
		  // ##############################
		  LoreVecK.Boost(-boostKst);
		  // ##############################
		  // # Boost B0 to K*0 ref. frame #
		  // ##############################
		  LoreVecB0.SetXYZM(NTupleIn->bPx->at(BestCandIndx),NTupleIn->bPy->at(BestCandIndx),NTupleIn->bPz->at(BestCandIndx),NTupleIn->bBarMass->at(BestCandIndx));
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
		  LoreVecB0.SetXYZM(NTupleIn->bPx->at(BestCandIndx),NTupleIn->bPy->at(BestCandIndx),NTupleIn->bPz->at(BestCandIndx),NTupleIn->bBarMass->at(BestCandIndx));
		  TVector3 boostB0 = LoreVecB0.BoostVector();
		  LoreVecMum.SetXYZM(NTupleIn->mumPx->at(BestCandIndx),NTupleIn->mumPy->at(BestCandIndx),NTupleIn->mumPz->at(BestCandIndx),Utility->muonMass);
		  TLorentzVector LoreVecMup;
		  LoreVecMup.SetXYZM(NTupleIn->mupPx->at(BestCandIndx),NTupleIn->mupPy->at(BestCandIndx),NTupleIn->mupPz->at(BestCandIndx),Utility->muonMass);
		  LoreVecK.SetXYZM(NTupleIn->kstTrkmPx->at(BestCandIndx),NTupleIn->kstTrkmPy->at(BestCandIndx),NTupleIn->kstTrkmPz->at(BestCandIndx),Utility->kaonMass);
		  TLorentzVector LoreVecPi;
		  LoreVecPi.SetXYZM(NTupleIn->kstTrkpPx->at(BestCandIndx),NTupleIn->kstTrkpPy->at(BestCandIndx),NTupleIn->kstTrkpPz->at(BestCandIndx),Utility->pionMass);

		  LoreVecMum.Boost(-boostB0);
		  LoreVecMup.Boost(-boostB0);
		  LoreVecK.Boost(-boostB0);
		  LoreVecPi.Boost(-boostB0);
		  TVector3 MuMuPlane = LoreVecMum.Vect().Cross(LoreVecMup.Vect());
		  TVector3 KstPlane  = LoreVecK.Vect().Cross(LoreVecPi.Vect());
		  if (MuMuPlane.Cross(KstPlane).Dot(-LoreVecB0.Vect()) > 0.0) phiKstMuMuPlane = MuMuPlane.Angle(KstPlane);
		  else phiKstMuMuPlane = -MuMuPlane.Angle(KstPlane);
 		}
	      else if (DoMCTruth == true)
		{
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
 		}

	      NTupleOut->B0MassArb          = NTupleIn->bBarMass->at(BestCandIndx);
	      NTupleOut->CosThetaMuArb      = cosThetaMum;
	      NTupleOut->CosThetaKArb       = cosThetaK;
	      NTupleOut->PhiKstMuMuPlaneArb = phiKstMuMuPlane;
	    }

	  NTupleOut->TrigCat    = TrigCat;

	  NTupleOut->B0notB0bar = B0notB0bar;
	  NTupleOut->B0pT       = sqrt(NTupleIn->bPx->at(BestCandIndx)*NTupleIn->bPx->at(BestCandIndx) + NTupleIn->bPy->at(BestCandIndx)*NTupleIn->bPy->at(BestCandIndx));
	  NTupleOut->B0Eta      = Utility->computeEta (NTupleIn->bPx->at(BestCandIndx),NTupleIn->bPy->at(BestCandIndx),NTupleIn->bPz->at(BestCandIndx));
	  NTupleOut->B0Phi      = Utility->computePhi (NTupleIn->bPx->at(BestCandIndx),NTupleIn->bPy->at(BestCandIndx),NTupleIn->bPz->at(BestCandIndx));

	  theTreeOut->Fill();
	  NTupleOut->ClearNTuple();
	}
    }
}


void BestCandPerformance ()
{
  // #######################################################
  // # Compute efficiency, purity and B-tagging efficiency #
  // #######################################################
  unsigned int coutGoodEv          = 0;
  unsigned int countSeleEv         = 0;
  unsigned int countGoodBTag       = 0;
  unsigned int countGoodBbarTag    = 0;
  unsigned int countGoodBSeleEv    = 0;
  unsigned int countGoodBbarSeleEv = 0;
  unsigned int countGoodSeleEv     = 0;

  unsigned int countCands;
  bool B0notB0bar;
  int BestCandIndx;
  int TrigCat;

  unsigned int counter  = 0;
  double countCandsAVG = 0.0;

  NTupleIn->ClearNTuple();
  NTupleIn->SetBranchAddresses(theTreeIn);
  int nEntries = theTreeIn->GetEntries();
  cout << "\n@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;


  // ###########################
  // # Loop over B0 Candidates #
  // ###########################
  for (int entry = 0; entry < nEntries; entry++)
    {
      theTreeIn->GetEntry(entry);


      // #####################
      // # Count good events #
      // #####################
      for (unsigned int i = 0; i < NTupleIn->bMass->size(); i++)
	if (((NTupleIn->genSignal == SignalType) || (NTupleIn->genSignal == SignalType+1)) && (NTupleIn->truthMatchSignal->at(i) == true))
	  {
	    coutGoodEv++;
	    break;
	  }


      if (Utility->ChooseBestCand(NTupleIn, DoTrigCheck, ((double)entry)/((double)nEntries), &BestCandIndx, &B0notB0bar, &TrigCat, &countCands) == true)
	{
	  countSeleEv++;

	  if ((NTupleIn->genSignal == SignalType) && (NTupleIn->truthMatchSignal->at(BestCandIndx) == true) && (B0notB0bar == true)) countGoodBTag++;
	  if ((NTupleIn->genSignal == SignalType) && (NTupleIn->truthMatchSignal->at(BestCandIndx) == true)) countGoodBSeleEv++;

	  if ((NTupleIn->genSignal == SignalType+1) && (NTupleIn->truthMatchSignal->at(BestCandIndx) == true) && (B0notB0bar == false)) countGoodBbarTag++;
	  if ((NTupleIn->genSignal == SignalType+1) && (NTupleIn->truthMatchSignal->at(BestCandIndx) == true)) countGoodBbarSeleEv++;

	  if (((NTupleIn->genSignal == SignalType) || (NTupleIn->genSignal == SignalType+1)) && (NTupleIn->truthMatchSignal->at(BestCandIndx) == true)) countGoodSeleEv++;
	}


      if (countCands != 0)
	{
	  counter++;
	  countCandsAVG += countCands;
	}
    }
  countCandsAVG = countCandsAVG / ((double)counter);

  cout << "\n@@@ Efficiency: " << ((double)countGoodSeleEv) / ((double)coutGoodEv) * 100. << " @@@" << endl;
  cout << "@@@ Purity: " << ((double)countGoodSeleEv) / ((double)countSeleEv) * 100. << " @@@" << endl;
  cout << "@@@ Total B-tagging efficiency: " << ((double)(countGoodBTag+countGoodBbarTag)) / ((double)countGoodSeleEv) * 100. << " @@@" << endl;
  cout << "@@@ B-tagging efficiency: " << ((double)countGoodBbarTag) / ((double)countGoodBSeleEv) * 100. << " @@@" << endl;
  cout << "@@@ Bbar-tagging efficiency: " << ((double)countGoodBbarTag) / ((double)countGoodBbarSeleEv) * 100. << " @@@" << endl;
  cout << "@@@ Average number of candidates after all cuts applied but B0-vtx CL: " << countCandsAVG << " @@@" << endl;
}


int main (int argc, char** argv)
{
  if (argc >= 3)
    {
      string option = argv[1];
      if (((option == "singlecand") && (argc == 4)) || (((option == "eff")) && (argc == 3)))
	{
	  string fileNameIn = argv[2];
	  string fileNameOut = "";
	  if (argc == 4) fileNameOut = argv[3];

	  Utility = new Utils();
	  Utility->ReadTriggerPathsANDCutsANDEntries(ParameterFILE);
	  Utility->ReadSelectionCuts(ParameterFILE);


	  cout << "\n@@@ Settings @@@" << endl;
	  cout << "Signal Type: "      << SignalType << endl;
	  cout << "Do Trigger Check: " << DoTrigCheck << endl;
	  cout << "SaveHistos: "       << SaveHistos << endl;
	  cout << "DoMCTruth: "        << DoMCTruth << endl;
	  cout << "TagFromTruth: "     << TagFromTruth << endl;
	  cout << "ParameterFILE: "    << ParameterFILE << endl;


	  TFile* NtplFileIn  = new TFile(fileNameIn.c_str(), "READ");
	  TFile* NtplFileOut = NULL;
	  theTreeIn          = (TTree*) NtplFileIn->Get("B0Cand/B0KstMuMuNTuple");
	  NTupleIn           = new B0KstMuMuTreeContent();
	  NTupleIn->Init();

	  if (option == "singlecand")
	    {
	      NtplFileOut = new TFile(fileNameOut.c_str(), "RECREATE");
	      NtplFileOut->mkdir("B0SingleCand");
	      NtplFileOut->cd("B0SingleCand");
	      theTreeOut  = new TTree("B0KstMuMuSingleCandNTuple","B0KstMuMuSingleCandNTuple");
	      NTupleOut   = new B0KstMuMuSingleCandTreeContent();
	      NTupleOut->Init();

	      SelectBestCand();
	    }
	  else if (option == "eff") BestCandPerformance();
	  else
	    {
	      cout << "[SingleCandB0KstMuMu::main]\tWrong option: " << option << endl;
	      exit(1);
	    }
	  cout << "\n@@@ Single candidate optimization cut done @@@" << endl;


	  // #############
	  // # Save Tree #
	  // #############
	  if (option == "singlecand")
	    {
	      NtplFileOut->cd("B0SingleCand");
	      theTreeOut->Write();
	      NtplFileOut->Close();
	      NtplFileIn->Close();
	  
	      NTupleOut->Destroy();
	      delete NTupleOut;
	    }
	  return 0;
	}
      else
	{
	  cout << "Parameter missing: " << endl;
	  cout << "./SingleCandB0KstMuMu [singlecand eff] inputFile.root [outputFile.root]" << endl;

	  return 1;
	}
    }
  else
    {
      cout << "Parameter missing: " << endl;
      cout << "./SingleCandB0KstMuMu [singlecand eff] inputFile.root [outputFile.root]" << endl;

      return 1;
    }

  return 0;
}

#include "../interface/Utils.h"
#include "../interface/ReadParameters.h"

#include <TAxis.h>
#include <TMath.h>

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>

// ####################
// # Global constants #
// ####################
#define YvalueOutsideLimits 10.0 // Value given to bins with zero error in order not to show them

Utils::Utils ()
{
  muonMass      = 0.10565837;
  pionMass      = 0.13957018;
  kaonMass      = 0.493677;
  kstMass       = 0.896;
  B0Mass        = 5.27958;
  JPsiMass      = 3.096916;
  PsiPrimeMass  = 3.686109;

  JPsiBF        =  7.95e-5; // B0 --> J/psi(mu+mu-) K*0          (1.34+/-0.06 * 5.93+/-0.06)
  JPsiKpiBF     =  5.30e-5; // B0 --> J/psi(mu+mu-) K*0(K+pi-)   (1.34+/-0.06 * 5.93+/-0.06 * 2/3)
  KstMuMuBF     =  1.06e-6; // B0 --> K*0 mu+mu-
  KstKpiMuMuBF  =  7.07e-7; // B0 --> K*0(K+pi-) mu+mu-          (1.06+/-0.10 * 2/3)
  PsiPBF        = 46.97e-7; // B0 --> psi(2S)(mu+mu-) K*0        (6.1+/-0.5 * 7.7+/-0.8)
  PsiPKpiBF     = 31.31e-7; // B0 --> psi(2S)(mu+mu-) K*0(K+pi-) (6.1+/-0.5 * 7.7+/-0.8 * 2/3)

  muonMassErr   = 0.00000035;
  pionMassErr   = 0.00000035;
  kaonMassErr   = 0.000016;
  B0MassErr     = 0.00017;
  kstSigma      = 0.05;

  nFitParam     = 51;
  nConfigParam  = 4;
  nFitObserv    = 5;

  NcoeffThetaL  = 5;
  NcoeffThetaK  = 4;

  PI = 3.141592653589793;

  ProbThreshold = 0.0; // Confidence Level for accepting the null hypothesis: the two mass hypothesis are statistically indistinguishable
  // if (F-test < ProbThreshold) --> accept the null hypothesis
  // if (F-test > ProbThreshold) --> discard the null hypothesis
  scrambleFraction = 0.0; // Fraction of events with random CP-tagging
  KstMassShape = new TF1("KstMassShape",
			 "2*sqrt(2)*[0]*[1]* sqrt([0]*[0] * ([0]*[0] + [1]*[1])) / (TMath::Pi()*sqrt([0]*[0] + sqrt([0]*[0] * ([0]*[0] + [1]*[1])))) / ((x*x - [0]*[0]) * (x*x - [0]*[0]) + [0]*[0]*[1]*[1])",
			 0.0,kstMass*2.);
  // Breit-Wigner distribution:
  // [0]: mass of the resonance
  // [1]: width of the resonance
  KstMassShape->SetParameter(0,kstMass);
  KstMassShape->SetParameter(1,kstSigma);
  KstMassShape->SetParNames("Mean","Width");


  // ################################
  // # Print out internal variables #
  // ################################
  std::cout << "\n@@@ Utils class settings @@@" << std::endl;
  std::cout << "nFitParam: "        << nFitParam << std::endl;
  std::cout << "nConfigParam: "     << nConfigParam << std::endl;
  std::cout << "nFitObserv: "       << nFitObserv << std::endl;
  std::cout << "NcoeffThetaL: "     << NcoeffThetaL << std::endl;
  std::cout << "NcoeffThetaK: "     << NcoeffThetaK << std::endl;
  std::cout << "ProbThreshold: "    << ProbThreshold << std::endl;
  std::cout << "scrambleFraction: " << scrambleFraction << std::endl;
  std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
}

double Utils::computeInvMass (double Px1,
			      double Py1,
			      double Pz1,
			      double mass1,
			      double Px2,
			      double Py2,
			      double Pz2,
			      double mass2)
{
  double Energy1 = sqrt(Px1*Px1 + Py1*Py1 + Pz1*Pz1 + mass1*mass1);
  double Energy2 = sqrt(Px2*Px2 + Py2*Py2 + Pz2*Pz2 + mass2*mass2);
  return sqrt((Energy1+Energy2) * (Energy1+Energy2) - ((Px1+Px2) * (Px1+Px2) + (Py1+Py2) * (Py1+Py2) + (Pz1+Pz2) * (Pz1+Pz2)));
}

double Utils::computeEta (double Px,
			  double Py,
			  double Pz)
{
  double P = sqrt(Px*Px + Py*Py + Pz*Pz);
  return 0.5*log((P + Pz) / (P - Pz));
}

double Utils::computePhi (double Px,
			  double Py,
			  double Pz)
{
  double phi = atan(Py / Px);
  if (Px < 0 && Py < 0) phi = phi - PI;
  if (Px < 0 && Py > 0) phi = phi + PI;
  return phi;
}

double Utils::computeEtaPhiDistance (double Px1,
				     double Py1,
				     double Pz1,
				     double Px2,
				     double Py2,
				     double Pz2)
{
  double phi1 = computePhi (Px1,Py1,Pz1);
  double eta1 = computeEta (Px1,Py1,Pz1);
  double phi2 = computePhi (Px2,Py2,Pz2);
  double eta2 = computeEta (Px2,Py2,Pz2);
  return sqrt((eta1-eta2) * (eta1-eta2) + (phi1-phi2) * (phi1-phi2));
}

void Utils::computeLS (double Vx,
		       double Vy,
		       double Vz,
		       double Wx,
		       double Wy,
		       double Wz,
		       double VxErr2,
		       double VyErr2,
		       double VzErr2,
		       double VxyCov,
		       double VxzCov,
		       double VyzCov,
		       double WxErr2,
		       double WyErr2,
		       double WzErr2,
		       double WxyCov,
		       double WxzCov,
		       double WyzCov,
		       double* deltaD,
		       double* deltaDErr)
{
  *deltaD = sqrt((Vx-Wx) * (Vx-Wx) + (Vy-Wy) * (Vy-Wy) + (Vz-Wz) * (Vz-Wz));
  if (*deltaD > 0.)
    *deltaDErr = sqrt((Vx-Wx) * (Vx-Wx) * VxErr2 +
		      (Vy-Wy) * (Vy-Wy) * VyErr2 +
		      (Vz-Wz) * (Vz-Wz) * VzErr2 +
		      
		      (Vx-Wx) * (Vy-Wy) * 2.*VxyCov +
		      (Vx-Wx) * (Vz-Wz) * 2.*VxzCov +
		      (Vy-Wy) * (Vz-Wz) * 2.*VyzCov +
		      
		      (Vx-Wx) * (Vx-Wx) * WxErr2 +
		      (Vy-Wy) * (Vy-Wy) * WyErr2 +
		      (Vz-Wz) * (Vz-Wz) * WzErr2 +
		      
		      (Vx-Wx) * (Vy-Wy) * 2.*WxyCov +
		      (Vx-Wx) * (Vz-Wz) * 2.*WxzCov +
		      (Vy-Wy) * (Vz-Wz) * 2.*WyzCov) / *deltaD;
  else *deltaDErr = 0.;
}

void Utils::computeCosAlpha (double Vx,
			     double Vy,
			     double Vz,
			     double Wx,
			     double Wy,
			     double Wz,
			     double VxErr2,
			     double VyErr2,
			     double VzErr2,
			     double VxyCov,
			     double VxzCov,
			     double VyzCov,
			     double WxErr2,
			     double WyErr2,
			     double WzErr2,
			     double WxyCov,
			     double WxzCov,
			     double WyzCov,
			     double* cosAlpha,
			     double* cosAlphaErr)
{
  double Vnorm = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
  double Wnorm = sqrt(Wx*Wx + Wy*Wy + Wz*Wz);
  double VdotW = Vx*Wx + Vy*Wy + Vz*Wz;
  
  if ((Vnorm > 0.) && (Wnorm > 0.))
    {
      *cosAlpha = VdotW / (Vnorm * Wnorm);
      *cosAlphaErr = sqrt( ((Vx*Wnorm - VdotW*Wx) * (Vx*Wnorm - VdotW*Wx) * WxErr2 +
			    (Vy*Wnorm - VdotW*Wy) * (Vy*Wnorm - VdotW*Wy) * WyErr2 +
			    (Vz*Wnorm - VdotW*Wz) * (Vz*Wnorm - VdotW*Wz) * WzErr2 +
			   
			    (Vx*Wnorm - VdotW*Wx) * (Vy*Wnorm - VdotW*Wy) * 2.*WxyCov +
			    (Vx*Wnorm - VdotW*Wx) * (Vz*Wnorm - VdotW*Wz) * 2.*WxzCov +
			    (Vy*Wnorm - VdotW*Wy) * (Vz*Wnorm - VdotW*Wz) * 2.*WyzCov) /
			   (Wnorm*Wnorm*Wnorm*Wnorm) +
			   
			   ((Wx*Vnorm - VdotW*Vx) * (Wx*Vnorm - VdotW*Vx) * VxErr2 +
			    (Wy*Vnorm - VdotW*Vy) * (Wy*Vnorm - VdotW*Vy) * VyErr2 +
			    (Wz*Vnorm - VdotW*Vz) * (Wz*Vnorm - VdotW*Vz) * VzErr2 +
			    
			    (Wx*Vnorm - VdotW*Vx) * (Wy*Vnorm - VdotW*Vy) * 2.*VxyCov +
			    (Wx*Vnorm - VdotW*Vx) * (Wz*Vnorm - VdotW*Vz) * 2.*VxzCov +
			    (Wy*Vnorm - VdotW*Vy) * (Wz*Vnorm - VdotW*Vz) * 2.*VyzCov) /
			   (Vnorm*Vnorm*Vnorm*Vnorm) ) / (Wnorm*Vnorm);
    }
  else
    {
      *cosAlpha = 0.;
      *cosAlphaErr = 0.;
    }
}

void Utils::ReadBins (std::string fileName, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins)
{
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // #################
  // # Read q^2 bins #
  // #################
  std::cout << "\n@@@ q^2 bins from file @@@" << std::endl;
  ParVector.clear();
  ParameterFile->ReadFromFile(ParFileBlockN("q2"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++)
    {
      q2Bins->push_back(atof(ParVector[i].c_str()));
      std::cout << "Bin " << i << "\t" << q2Bins->back() << std::endl;
    }


  // #######################
  // # Read cosThetaK bins #
  // #######################
  std::cout << "\n@@@ cos(theta_K) bins from file @@@" << std::endl;
  ParVector.clear();
  ParameterFile->ReadFromFile(ParFileBlockN("thetaK"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++)
    {
      cosThetaKBins->push_back(atof(ParVector[i].c_str()));
      std::cout << "Bin " << i << "\t" << cosThetaKBins->back() << std::endl;
    }


  // #######################
  // # Read cosThetaL bins #
  // #######################
  std::cout << "\n@@@ cos(theta_l) bins from file @@@" << std::endl;
  ParVector.clear();
  ParameterFile->ReadFromFile(ParFileBlockN("thetaL"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++)
    {
      cosThetaLBins->push_back(atof(ParVector[i].c_str()));
      std::cout << "Bin " << i << "\t" << cosThetaLBins->back() << std::endl;
    }


  // #################
  // # Read phi bins #
  // #################
  std::cout << "\n@@@ phi bins from file @@@" << std::endl;
  ParVector.clear();
  ParameterFile->ReadFromFile(ParFileBlockN("phi"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++)
    {
      phiBins->push_back(atof(ParVector[i].c_str()));
      std::cout << "Bin " << i << "\t" << phiBins->back() << std::endl;
    }


  ParVector.clear();
  delete ParameterFile;
}

void Utils::Readq2Bins (std::string fileName, std::vector<double>* q2Bins)
{
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // #################
  // # Read q^2 bins #
  // #################
  std::cout << "\n@@@ q^2 bins from file @@@" << std::endl;
  ParVector.clear();
  ParameterFile->ReadFromFile(ParFileBlockN("q2"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++)
    {
      q2Bins->push_back(atof(ParVector[i].c_str()));
      std::cout << "Bin " << i << "\t" << q2Bins->back() << std::endl;
    }


  ParVector.clear();
  delete ParameterFile;
}

void Utils::ReadHLTpaths (std::string fileName, std::vector<std::string>* TrigTable)
{
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ###############################
  // # Read HLT-trigger table bins #
  // ###############################
  std::cout << "\n@@@ HLT-trigger table from file @@@" << std::endl;
  ParVector.clear();
  ParameterFile->ReadFromFile(ParFileBlockN("HLTpath"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++)
    {
      TrigTable->push_back(ParVector[i]);
      std::cout << "Trigger path from config file: " << TrigTable->operator[](i).c_str() << std::endl;
    }


  ParVector.clear();
  delete ParameterFile;
}

void Utils::GenEfficiency (effStruct* myEff, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins)
{
  myEff->Num1 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
  myEff->Num2 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
  myEff->Den1 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
  myEff->Den2 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];

  myEff->Err2PoisNum1 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
  myEff->Err2PoisNum2 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
  myEff->Err2PoisDen1 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
  myEff->Err2PoisDen2 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];

  myEff->Err2WeigNum1 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
  myEff->Err2WeigNum2 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
  myEff->Err2WeigDen1 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
  myEff->Err2WeigDen2 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
}

void Utils::InitEfficiency (effValue myEffVal, effStruct myEff, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins)
{
  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	for (unsigned int l = 0; l < phiBins->size()-1; l++)
	  {
	    myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Num1;
	    myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Num2;
	    myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Den1;
	    myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Den2;

	    myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2PoisNum1;
	    myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2PoisNum2;
	    myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2PoisDen1;
	    myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2PoisDen2;

	    myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2WeigNum1;
	    myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2WeigNum2;
	    myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2WeigDen1;
	    myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2WeigDen2;
	  }
}

void Utils::SaveEfficiency (std::string fileName, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, effStruct myEff)
{
  effValue myEffVal;

  std::ofstream out;
  out.open(fileName.c_str(), std::ofstream::out);
  if (out.good() == false)
    {
      std::cout << "[Utils::SaveMatrices]\tError opening file : " << fileName << std::endl;
      exit (1);
    }

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	for (unsigned int l = 0; l < phiBins->size()-1; l++)
	  {
	    myEffVal.Num1 = myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Num2 = myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Den1 = myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Den2 = myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	    myEffVal.Err2PoisNum1 = myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2PoisNum2 = myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2PoisDen1 = myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2PoisDen2 = myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	    myEffVal.Err2WeigNum1 = myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2WeigNum2 = myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2WeigDen1 = myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2WeigDen2 = myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	    out << myEffVal.Num1 << "   " << myEffVal.Num2 << "   " << myEffVal.Den1 << "   " << myEffVal.Den2 << "   ";
	    out << myEffVal.Err2PoisNum1 << "   " << myEffVal.Err2PoisNum2 << "   " << myEffVal.Err2PoisDen1 << "   " << myEffVal.Err2PoisDen2 << "   ";
	    out << myEffVal.Err2WeigNum1 << "   " << myEffVal.Err2WeigNum2 << "   " << myEffVal.Err2WeigDen1 << "   " << myEffVal.Err2WeigDen2 << std::endl;
	  }

  out.close();
}

void Utils::ReadEfficiency (std::string fileName, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, effStruct* myEff)
{
  std::string tmpString;
  effValue myEffVal;

  std::ifstream in;
  in.open(fileName.c_str(), std::ofstream::in);
  if (in.good() == false)
    {
      std::cout << "[Utils::ReadEfficiency]\tError opening file : " << fileName << std::endl;
      exit (1);
    }
  
  GenEfficiency(myEff,q2Bins,cosThetaKBins,cosThetaLBins,phiBins);
  ResetEffValue(&myEffVal,0.0);
  InitEfficiency(myEffVal,*myEff,q2Bins,cosThetaKBins,cosThetaLBins,phiBins);

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	for (unsigned int l = 0; l < phiBins->size()-1; l++)
	  {
	    getline(in,tmpString);
	    std::stringstream rawString(tmpString);
	    rawString >> myEffVal.Num1 >> myEffVal.Num2 >> myEffVal.Den1 >> myEffVal.Den2 >> myEffVal.Err2PoisNum1 >> myEffVal.Err2PoisNum2 >> myEffVal.Err2PoisDen1 >> myEffVal.Err2PoisDen2 >> myEffVal.Err2WeigNum1 >> myEffVal.Err2WeigNum2 >> myEffVal.Err2WeigDen1 >> myEffVal.Err2WeigDen2;

	    myEff->Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Num1;
	    myEff->Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Num2;
	    myEff->Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Den1;
	    myEff->Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Den2;

	    myEff->Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2PoisNum1;
	    myEff->Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2PoisNum2;
	    myEff->Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2PoisDen1;
	    myEff->Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2PoisDen2;

	    myEff->Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2WeigNum1;
	    myEff->Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2WeigNum2;
	    myEff->Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2WeigDen1;
	    myEff->Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2WeigDen2;
	  }

  in.close();
}

void Utils::DeleteEfficiency (effStruct myEff)
{
  delete myEff.Num1;
  delete myEff.Num2;
  delete myEff.Den1;
  delete myEff.Den2;
}

int Utils::SearchBin (double val2Search, std::vector<double>* bins)
{
  unsigned int i = 0;
  for (i = 0; i < bins->size() - 1; i++) if ((val2Search >= bins->operator[](i)) && (val2Search < bins->operator[](i+1))) break;
  
  if (i != bins->size() - 1) return i;
  return -1;
}

int Utils::GetJPsiBin (std::vector<double>* q2Bins)
{
  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    if ((q2Bins->operator[](i) < (JPsiMass*JPsiMass)) && (q2Bins->operator[](i+1) > (JPsiMass*JPsiMass)))
      return i;
  
  return -1;
}

int Utils::GetPsiPBin (std::vector<double>* q2Bins)
{
  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    if ((q2Bins->operator[](i) < (PsiPrimeMass*PsiPrimeMass)) && (q2Bins->operator[](i+1) > (PsiPrimeMass*PsiPrimeMass)))
      return i;
  
  return -1;
}

bool Utils::ValIsInPsi (std::vector<double>* q2Bins, double q2Val)
{
  unsigned int JPsibin = GetJPsiBin(q2Bins);
  unsigned int PsiPbin = GetPsiPBin(q2Bins);

  if (((q2Bins->operator[](JPsibin) < q2Val) &&
       (q2Bins->operator[](JPsibin+1) > q2Val)) ||
      ((q2Bins->operator[](PsiPbin) < q2Val) &&
       (q2Bins->operator[](PsiPbin+1) > q2Val)))
    return true;
  
  return false;
}

bool Utils::ValIsBetweenJPsiAndPsiP (std::vector<double>* q2Bins, double q2Val)
{
  unsigned int JPsibin = GetJPsiBin(q2Bins);
  unsigned int PsiPbin = GetPsiPBin(q2Bins);

  if ((q2Val >= q2Bins->operator[](JPsibin+1)) && (q2Val <= q2Bins->operator[](PsiPbin))) return true;
  return false;
}

bool Utils::IsGoodq2CosThetaKCosThetaLBin (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, double mumuMass2, double cosThetaK, double cosThetaMu, bool PrintMsg)
{
  int mumuq2BinIndx;
  int cosThetaKBinIndx;
  int cosThetaLBinIndx;

  mumuq2BinIndx    = SearchBin(mumuMass2,q2Bins);
  cosThetaKBinIndx = SearchBin(cosThetaK,cosThetaKBins);
  cosThetaLBinIndx = SearchBin(cosThetaMu,cosThetaLBins);

  if (PrintMsg == true)
    {
      std::cout << "Bin assignement: ";
      std::cout << "\tq2 = " << mumuMass2 << " --> bin: " << mumuq2BinIndx;
      std::cout << "\tcos(theta_K) = " << cosThetaK << " --> bin: " << cosThetaKBinIndx;
      std::cout << "\tcos(theta_l) = " << cosThetaMu << " --> bin: " << cosThetaLBinIndx << std::endl;
    }

  if ((mumuq2BinIndx != -1) && (cosThetaKBinIndx != -1) && (cosThetaLBinIndx != -1)) return true;
  return false;
}

void Utils::IntegrateEffAll (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;

  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	for (unsigned int l = 0; l < phiBins->size()-1; l++)
	  {
	    myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	    if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	      myEffValOrg.Num1 = myEffValOrg.Num1 + pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	    if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	      myEffValOrg.Num2 = myEffValOrg.Num2 + pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	    if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	      myEffValOrg.Den1 = myEffValOrg.Den1 + pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	    if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	      myEffValOrg.Den2 = myEffValOrg.Den2 + pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	    myEffVal.Err2WeigNum1=myEffVal.Err2WeigNum1+myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2WeigNum2=myEffVal.Err2WeigNum2+myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2WeigDen1=myEffVal.Err2WeigDen1+myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2WeigDen2=myEffVal.Err2WeigDen2+myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  }
  
  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && myEffValOrg.Den1 > 0.0 && myEffValOrg.Num2 > 0.0 && myEffValOrg.Den2 > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / (myEffVal.Num1*myEffVal.Num1);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / (myEffVal.Num2*myEffVal.Num2);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / (myEffVal.Den1*myEffVal.Den1);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / (myEffVal.Den2*myEffVal.Den2);

      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((1. / myEffValOrg.Num1 - 1. / myEffValOrg.Den1 + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) + (1. / myEffValOrg.Num2 - 1. / myEffValOrg.Den2 + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffButPsi (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;

  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	for (unsigned int l = 0; l < phiBins->size()-1; l++)
	  if (!(((q2Bins->operator[](i) < (JPsiMass*JPsiMass)) && (q2Bins->operator[](i+1) > (JPsiMass*JPsiMass))) ||
		((q2Bins->operator[](i) < (PsiPrimeMass*PsiPrimeMass)) && (q2Bins->operator[](i+1) > (PsiPrimeMass*PsiPrimeMass)))))
	    {
	      myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Num1 = myEffValOrg.Num1 + pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	      if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Num2 = myEffValOrg.Num2 + pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Den1 = myEffValOrg.Den1 + pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	      if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Den2 = myEffValOrg.Den2 + pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      myEffVal.Err2WeigNum1=myEffVal.Err2WeigNum1+myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2WeigNum2=myEffVal.Err2WeigNum2+myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2WeigDen1=myEffVal.Err2WeigDen1+myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2WeigDen2=myEffVal.Err2WeigDen2+myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    }
  
  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && myEffValOrg.Den1 > 0.0 && myEffValOrg.Num2 > 0.0 && myEffValOrg.Den2 > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / (myEffVal.Num1*myEffVal.Num1);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / (myEffVal.Num2*myEffVal.Num2);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / (myEffVal.Den1*myEffVal.Den1);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / (myEffVal.Den2*myEffVal.Den2);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((1. / myEffValOrg.Num1 - 1. / myEffValOrg.Den1 + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) + (1. / myEffValOrg.Num2 - 1. / myEffValOrg.Den2 + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffInJPsi (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;

  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	for (unsigned int l = 0; l < phiBins->size()-1; l++)
	  if ((q2Bins->operator[](i) < (JPsiMass*JPsiMass)) && (q2Bins->operator[](i+1) > (JPsiMass*JPsiMass)))
	    {
	      myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Num1 = myEffValOrg.Num1 + pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	      if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Num2 = myEffValOrg.Num2 + pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Den1 = myEffValOrg.Den1 + pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	      if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Den2 = myEffValOrg.Den2 + pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 +
		myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 +
		myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 +
		myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 +
		myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    }
  
  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && myEffValOrg.Den1 > 0.0 && myEffValOrg.Num2 > 0.0 && myEffValOrg.Den2 > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / (myEffVal.Num1*myEffVal.Num1);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / (myEffVal.Num2*myEffVal.Num2);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / (myEffVal.Den1*myEffVal.Den1);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / (myEffVal.Den2*myEffVal.Den2);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((1. / myEffValOrg.Num1 - 1. / myEffValOrg.Den1 + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) + (1. / myEffValOrg.Num2 - 1. / myEffValOrg.Den2 + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffInPsiP (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;

  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	for (unsigned int l = 0; l < phiBins->size()-1; l++)
	  if ((q2Bins->operator[](i) < (PsiPrimeMass*PsiPrimeMass)) && (q2Bins->operator[](i+1) > (PsiPrimeMass*PsiPrimeMass)))
	    {
	      myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Num1 = myEffValOrg.Num1 + pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	      if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Num2 = myEffValOrg.Num2 + pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Den1 = myEffValOrg.Den1 + pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	      if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Den2 = myEffValOrg.Den2 + pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 +
		myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2WeigNum2=myEffVal.Err2WeigNum2 +
		myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 +
		myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 +
		myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    }

  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && myEffValOrg.Den1 > 0.0 && myEffValOrg.Num2 > 0.0 && myEffValOrg.Den2 > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / (myEffVal.Num1*myEffVal.Num1);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / (myEffVal.Num2*myEffVal.Num2);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / (myEffVal.Den1*myEffVal.Den1);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / (myEffVal.Den2*myEffVal.Den2);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((1. / myEffValOrg.Num1 - 1. / myEffValOrg.Den1 + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) + (1. / myEffValOrg.Num2 - 1. / myEffValOrg.Den2 + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffButq2 (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int q2Indx, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;

  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
    for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	  myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	  myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	  myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	  if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	    myEffValOrg.Num1 = myEffValOrg.Num1 + pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	      myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	    
	  if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	    myEffValOrg.Num2 = myEffValOrg.Num2 + pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	      myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	  if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	    myEffValOrg.Den1 = myEffValOrg.Den1 + pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	      myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	    
	  if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	    myEffValOrg.Den2 = myEffValOrg.Den2 + pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	      myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	  myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 +
	    myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	  myEffVal.Err2WeigNum2 =myEffVal.Err2WeigNum2 +
	    myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	  myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 +
	    myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	  myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 +
	    myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	}
  
  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && myEffValOrg.Den1 > 0.0 && myEffValOrg.Num2 > 0.0 && myEffValOrg.Den2 > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / (myEffVal.Num1*myEffVal.Num1);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / (myEffVal.Num2*myEffVal.Num2);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / (myEffVal.Den1*myEffVal.Den1);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / (myEffVal.Den2*myEffVal.Den2);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((1. / myEffValOrg.Num1 - 1. / myEffValOrg.Den1 + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) + (1. / myEffValOrg.Num2 - 1. / myEffValOrg.Den2 + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffButCosThetaK (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int cosThetaKBinIndx, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;

  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	  myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	  myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	  myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];

	  if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Num1 = myEffValOrg.Num1 +
	      pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	    
	  if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Num2 = myEffValOrg.Num2 +
	      pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];

	  if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Den1 = myEffValOrg.Den1 +
	      pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	    
	  if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Den2 = myEffValOrg.Den2 +
	      pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];

	  myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 +
	    myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	  myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 +
	    myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	  myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 +
	    myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	  myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 +
	    myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	}
  
  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && myEffValOrg.Den1 > 0.0 && myEffValOrg.Num2 > 0.0 && myEffValOrg.Den2 > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / (myEffVal.Num1*myEffVal.Num1);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / (myEffVal.Num2*myEffVal.Num2);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / (myEffVal.Den1*myEffVal.Den1);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / (myEffVal.Den2*myEffVal.Den2);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((1. / myEffValOrg.Num1 - 1. / myEffValOrg.Den1 + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) + (1. / myEffValOrg.Num2 - 1. / myEffValOrg.Den2 + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffButCosThetaL (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int cosThetaLBinIndx, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;

  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	  if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Num1 = myEffValOrg.Num1 +
	      pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	  if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Num2 = myEffValOrg.Num2 +
	      pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	  if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Den1 = myEffValOrg.Den1 +
	      pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	  if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Den2 = myEffValOrg.Den2 +
	      pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	  myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 +
	    myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Err2WeigNum2=myEffVal.Err2WeigNum2 +
	    myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 +
	    myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 +
	    myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	}
  
  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && myEffValOrg.Den1 > 0.0 && myEffValOrg.Num2 > 0.0 && myEffValOrg.Den2 > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / (myEffVal.Num1*myEffVal.Num1);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / (myEffVal.Num2*myEffVal.Num2);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / (myEffVal.Den1*myEffVal.Den1);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / (myEffVal.Den2*myEffVal.Den2);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((1. / myEffValOrg.Num1 - 1. / myEffValOrg.Den1 + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) + (1. / myEffValOrg.Num2 - 1. / myEffValOrg.Den2 + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffButPhi (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int phiIndx, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;

  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	{
	  myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	  if (myEff.Err2PoisNum1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Num1 = myEffValOrg.Num1 +
	      pow(myEff.Num1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisNum1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	  if (myEff.Err2PoisNum2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Num2 = myEffValOrg.Num2 +
	      pow(myEff.Num2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisNum2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	  if (myEff.Err2PoisDen1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Den1 = myEffValOrg.Den1 +
	      pow(myEff.Den1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisDen1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	  if (myEff.Err2PoisDen2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Den2 = myEffValOrg.Den2 +
	      pow(myEff.Den2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisDen2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	  myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 +
	    myEff.Err2WeigNum1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Err2WeigNum2=myEffVal.Err2WeigNum2 +
	    myEff.Err2WeigNum2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 +
	    myEff.Err2WeigDen1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 +
	    myEff.Err2WeigDen2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	}

  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && myEffValOrg.Den1 > 0.0 && myEffValOrg.Num2 > 0.0 && myEffValOrg.Den2 > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / (myEffVal.Num1*myEffVal.Num1);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / (myEffVal.Num2*myEffVal.Num2);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / (myEffVal.Den1*myEffVal.Den1);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / (myEffVal.Den2*myEffVal.Den2);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((1. / myEffValOrg.Num1 - 1. / myEffValOrg.Den1 + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) + (1. / myEffValOrg.Num2 - 1. / myEffValOrg.Den2 + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffPhiCosThetaL (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int q2Indx, unsigned int cosThetaKBinIndx, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;
  
  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
    for (unsigned int l = 0; l < phiBins->size()-1; l++)
      {
	myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
	myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
	myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
	myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];

	if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Num1 = myEffValOrg.Num1 +
	    pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
	    
	if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Num2 = myEffValOrg.Num2 +
	    pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];

	if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Den1 = myEffValOrg.Den1 +
	    pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
	    
	if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Den2 = myEffValOrg.Den2 +
	    pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];

	myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 +
	  myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2WeigNum2=myEffVal.Err2WeigNum2 +
	  myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 +
	  myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 +
	  myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
      }

  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && myEffValOrg.Den1 > 0.0 && myEffValOrg.Num2 > 0.0 && myEffValOrg.Den2 > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / (myEffVal.Num1*myEffVal.Num1);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / (myEffVal.Num2*myEffVal.Num2);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / (myEffVal.Den1*myEffVal.Den1);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / (myEffVal.Den2*myEffVal.Den2);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((1. / myEffValOrg.Num1 - 1. / myEffValOrg.Den1 + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) + (1. / myEffValOrg.Num2 - 1. / myEffValOrg.Den2 + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffPhiCosThetaK (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int q2Indx, unsigned int cosThetaLBinIndx, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;
  
  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);
  
  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
    for (unsigned int l = 0; l < phiBins->size()-1; l++)
      {
	myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Num1 = myEffValOrg.Num1 +
	    pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	    
	if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Num2 = myEffValOrg.Num2 +
	    pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Den1 = myEffValOrg.Den1 +
	    pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	    
	if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Den2 = myEffValOrg.Den2 +
	    pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 +
	  myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2WeigNum2=myEffVal.Err2WeigNum2 +
	  myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 +
	  myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 +
	  myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
      }
  
  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && myEffValOrg.Den1 > 0.0 && myEffValOrg.Num2 > 0.0 && myEffValOrg.Den2 > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / (myEffVal.Num1*myEffVal.Num1);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / (myEffVal.Num2*myEffVal.Num2);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / (myEffVal.Den1*myEffVal.Den1);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / (myEffVal.Den2*myEffVal.Den2);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((1. / myEffValOrg.Num1 - 1. / myEffValOrg.Den1 + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) + (1. / myEffValOrg.Num2 - 1. / myEffValOrg.Den2 + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffCosThetaKCosThetaL (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int q2Indx, unsigned int phiIndx, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;
  
  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);
  
  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
    for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
      {
	myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	if (myEff.Err2PoisNum1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Num1 = myEffValOrg.Num1 +
	    pow(myEff.Num1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisNum1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	    
	if (myEff.Err2PoisNum2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Num2 = myEffValOrg.Num2 +
	    pow(myEff.Num2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisNum2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	if (myEff.Err2PoisDen1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Den1 = myEffValOrg.Den1 +
	    pow(myEff.Den1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisDen1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	    
	if (myEff.Err2PoisDen2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Den2 = myEffValOrg.Den2 +
	    pow(myEff.Den2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisDen2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 +
	  myEff.Err2WeigNum1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2WeigNum2=myEffVal.Err2WeigNum2 +
	  myEff.Err2WeigNum2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 +
	  myEff.Err2WeigDen1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 +
	  myEff.Err2WeigDen2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
      }
  
  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && myEffValOrg.Den1 > 0.0 && myEffValOrg.Num2 > 0.0 && myEffValOrg.Den2 > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / (myEffVal.Num1*myEffVal.Num1);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / (myEffVal.Num2*myEffVal.Num2);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / (myEffVal.Den1*myEffVal.Den1);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / (myEffVal.Den2*myEffVal.Den2);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((1. / myEffValOrg.Num1 - 1. / myEffValOrg.Den1 + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) + (1. / myEffValOrg.Num2 - 1. / myEffValOrg.Den2 + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

bool Utils::IntegrateEffPhi (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, double mumuMass2, double cosThetaK, double cosThetaMu, effStruct myEff, double* Eff, double* EffErr, bool PrintMsg)
{
  int mumuq2BinIndx;
  int cosThetaKBinIndx;
  int cosThetaLBinIndx;
  
  effValue myEffVal;
  effValue myEffValOrg;
  
  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  if (IsGoodq2CosThetaKCosThetaLBin(q2Bins,cosThetaKBins,cosThetaLBins,mumuMass2,cosThetaK,cosThetaMu,PrintMsg) == true)
    {
      mumuq2BinIndx    = SearchBin(mumuMass2,q2Bins);
      cosThetaKBinIndx = SearchBin(cosThetaK,cosThetaKBins);
      cosThetaLBinIndx = SearchBin(cosThetaMu,cosThetaLBins);

      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
						     cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
						     cosThetaKBinIndx*(q2Bins->size()-1) +
						     mumuq2BinIndx];
	  myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
						     cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
						     cosThetaKBinIndx*(q2Bins->size()-1) +
						     mumuq2BinIndx];
	  myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
						     cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
						     cosThetaKBinIndx*(q2Bins->size()-1) +
						     mumuq2BinIndx];
	  myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
						     cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
						     cosThetaKBinIndx*(q2Bins->size()-1) +
						     mumuq2BinIndx];

	  if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaKBinIndx*(q2Bins->size()-1) +
				 mumuq2BinIndx] > 0.0)
	    myEffValOrg.Num1 = myEffValOrg.Num1 + pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
								 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
								 cosThetaKBinIndx*(q2Bins->size()-1) +
								 mumuq2BinIndx],2.0) /
	      myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaKBinIndx*(q2Bins->size()-1) +
				 mumuq2BinIndx];
	    
	  if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaKBinIndx*(q2Bins->size()-1) +
				 mumuq2BinIndx] > 0.0)
	    myEffValOrg.Num2 = myEffValOrg.Num2 + pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
								 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
								 cosThetaKBinIndx*(q2Bins->size()-1) +
								 mumuq2BinIndx],2.0) /
	      myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaKBinIndx*(q2Bins->size()-1) +
				 mumuq2BinIndx];

	  if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaKBinIndx*(q2Bins->size()-1) +
				 mumuq2BinIndx] > 0.0)
	    myEffValOrg.Den1 = myEffValOrg.Den1 + pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
								 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
								 cosThetaKBinIndx*(q2Bins->size()-1) +
								 mumuq2BinIndx],2.0) /
	      myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaKBinIndx*(q2Bins->size()-1) +
				 mumuq2BinIndx];
	    
	  if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaKBinIndx*(q2Bins->size()-1) +
				 mumuq2BinIndx] > 0.0)
	    myEffValOrg.Den2 = myEffValOrg.Den2 + pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
								 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
								 cosThetaKBinIndx*(q2Bins->size()-1) +
								 mumuq2BinIndx],2.0) /
	      myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaKBinIndx*(q2Bins->size()-1) +
				 mumuq2BinIndx];

	  myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 +
	    myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + mumuq2BinIndx];
	  myEffVal.Err2WeigNum2=myEffVal.Err2WeigNum2 +
	    myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + mumuq2BinIndx];
	  myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 +
	    myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + mumuq2BinIndx];
	  myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 +
	    myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + mumuq2BinIndx];
	}
      
      if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && myEffValOrg.Den1 > 0.0 && myEffValOrg.Num2 > 0.0 && myEffValOrg.Den2 > 0.0)
	{
	  myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / (myEffVal.Num1*myEffVal.Num1);
	  myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / (myEffVal.Num2*myEffVal.Num2);
	  myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / (myEffVal.Den1*myEffVal.Den1);
	  myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / (myEffVal.Den2*myEffVal.Den2);
	  
	  *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
	  *EffErr = *Eff * sqrt((1. / myEffValOrg.Num1 - 1. / myEffValOrg.Den1 + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) + (1. / myEffValOrg.Num2 - 1. / myEffValOrg.Den2 + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
	}
      else
	{
	  *Eff    = 0.0;
	  *EffErr = 0.0;
	}
      
      return true;
    }
  else return false;
}

unsigned int Utils::HLTpathForEvFraction (double evtFrac)
{
  // #####################################################################################
  // # This method is used together with the method: "ReadTriggerPathsANDCutsANDEntries" #
  // #####################################################################################
  double totalLumi = 0.0;
  double countLumi = 0.0;
  unsigned int it  = 0;

  for (unsigned int j = 0; j < VecHLTentries.size(); j++) totalLumi = totalLumi + VecHLTentries[j];
  if ((evtFrac < 0.0) || (totalLumi == 0.0)) totalLumi = 1.0;
  
  for (it = 0; it < HLTpath.size(); it++)
    {
      if (evtFrac < (countLumi + VecHLTentries[it])/totalLumi) break;
      countLumi = countLumi + VecHLTentries[it];
    }

  return (it < HLTpath.size()-1 ? it+1 : 1); // @TMP@: this is done in order to consider "HLT_Dimuon6p5_LowMass_Displaced" like "HLT_Dimuon7_LowMass_Displaced"
}

unsigned int Utils::IsInTriggerTable (B0KstMuMuTreeContent* NTupleIn, double* HLTCutVar1, double* HLTCutVar2, int index, double evtFrac)
// ##########################################################
// # if index == -1 just split the sample in HLT categories #
// # output > 0 if it's in trigger table                    #
// # else output = 0                                        #
// ##########################################################
{
  // #####################################################################################
  // # This method is used together with the method: "ReadTriggerPathsANDCutsANDEntries" #
  // #####################################################################################
  unsigned int it;
  unsigned int HLTpathIndx;

  for (unsigned int j = 0; j < HLTpath.size(); j++)
    {
      if (evtFrac >= 0.0) HLTpathIndx = HLTpathForEvFraction(evtFrac)-1;
      else HLTpathIndx = j;
      *HLTCutVar1 = VecHLTCutVar1[HLTpathIndx];
      *HLTCutVar2 = VecHLTCutVar2[HLTpathIndx];
      if (index == -1) return (j < HLTpath.size()-1 ? HLTpathIndx+1 : 1);

      for (it = 0; it < NTupleIn->TrigTable->size(); it++) if (NTupleIn->TrigTable->operator[](it).find(HLTpath[HLTpathIndx].c_str()) != std::string::npos) break;
      
      if ((it < NTupleIn->TrigTable->size()) &&
	  (NTupleIn->mumTrig->at(index).find(HLTpath[HLTpathIndx].c_str()) != std::string::npos) &&
	  (NTupleIn->mupTrig->at(index).find(HLTpath[HLTpathIndx].c_str()) != std::string::npos))
	return (j < HLTpath.size()-1 ? HLTpathIndx+1 : 1); // @TMP@: this is done in order to consider "HLT_Dimuon6p5_LowMass_Displaced" like "HLT_Dimuon7_LowMass_Displaced"

      if (evtFrac >= 0.0) break;
    }

  return 0;
}

double Utils::GetHLTentries (unsigned int HLTindx)
{
  // #####################################################################################
  // # This method is used together with the method: "ReadTriggerPathsANDCutsANDEntries" #
  // #####################################################################################
  if (HLTindx-1 < VecHLTentries.size()) return VecHLTentries[HLTindx-1];
  else return -1.0;
}

unsigned int Utils::GetNHLTCat ()
{
  // #####################################################################################
  // # This method is used together with the method: "ReadTriggerPathsANDCutsANDEntries" #
  // #####################################################################################
  return (HLTpath.size() > 0 ? HLTpath.size()-1 : 0); // @TMP@: this is done in order to consider "HLT_Dimuon6p5_LowMass_Displaced" like "HLT_Dimuon7_LowMass_Displaced"
}

void Utils::ReadTriggerPathsANDCutsANDEntries (std::string fileName)
{
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ############################################
  // # Read HLT-trigger paths and relative cuts #
  // ############################################
  ParameterFile->ReadFromFile(ParFileBlockN("HLTcuts"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i=i+4)
    {
      HLTpath.push_back(ParVector[i]);
      VecHLTCutVar1.push_back(atof(ParVector[i+1].c_str()));
      VecHLTCutVar2.push_back(atof(ParVector[i+2].c_str()));
      VecHLTentries.push_back(atof(ParVector[i+3].c_str()));

      std::cout << "\nRead trigger path from config file: " << HLTpath.back() << std::endl;
      std::cout << "Read first cut value: " << VecHLTCutVar1.back() << std::endl;
      std::cout << "Read second cut value: " << VecHLTCutVar2.back() << std::endl;
      std::cout << "Read entries in Data: " << VecHLTentries.back() << std::endl;
    }


  ParVector.clear();
  delete ParameterFile;
}

void Utils::ReadFitStartingValues (std::string fileName, std::vector<std::vector<std::string>*>* vecParam, std::vector<std::vector<unsigned int>*>* configParam, const unsigned int dataBlockN)
{
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ############################
  // # Read fit-starting values #
  // ############################
  ParameterFile->ReadFromFile(dataBlockN,&ParVector);

  for (unsigned int j = 0; j < nFitParam; j++) vecParam->push_back(new std::vector<std::string>);
  for (unsigned int j = 0; j < nConfigParam; j++) configParam->push_back(new std::vector<unsigned int>);
  
  for (unsigned int i = 0; i < ParVector.size(); i=i+nFitParam+nConfigParam)
    {

      std::cout << "\nRead set-" << (int)((double)i/(double)(nFitParam+nConfigParam)) << " of fit-parameter starting values" << std::endl;

      for (unsigned int j = 0; j < nFitParam; j++)
	{
	  std::stringstream rawString(ParVector[i+j]);
	  vecParam->operator[](j)->push_back(rawString.str());
	  std::cout << "Fit parameter-" << j << " value: " << vecParam->operator[](j)->back() << std::endl;
	}
      for (unsigned int j = 0; j < nConfigParam; j++)
	{
	  configParam->operator[](j)->push_back(atoi(ParVector[i+nFitParam+j].c_str()));
	  std::cout << "Config. parameter-" << j << " value: " << configParam->operator[](j)->back() << std::endl;
	}
    }

  
  ParVector.clear();
  delete ParameterFile;
}

void Utils::ReadFitSystematics (std::string fileName, std::vector<std::vector<double>*>* vecParam)
// vecParam[0] --> Fl +err
// vecParam[1] --> Fl -err
// vecParam[2] --> Afb +err
// vecParam[3] --> Afb -err
// vecParam[4] --> At2 +err
// vecParam[5] --> At2 -err
// vecParam[6] --> Atim +err
// vecParam[7] --> Atim -err
// vecParam[8] --> Branching-Fraction +err
// vecParam[9] --> Branching-Fraction -err
{
  double val;
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // #########################################
  // # Read fit-observable systematic errors #
  // #########################################
  ParameterFile->ReadFromFile(ParFileBlockN("fitSyst"),&ParVector);

  for (unsigned int j = 0; j < nFitObserv*2; j++) vecParam->push_back(new std::vector<double>);
  
  for (unsigned int i = 0; i < ParVector.size(); i=i+nFitObserv)
    {

      std::cout << "\nRead set-" << (int)((double)i/(double)(nFitObserv)) << " of fit-observable systematic errors" << std::endl;

      for (unsigned int j = 0; j < nFitObserv; j++)
	{
	  std::stringstream rawString(ParVector[i+j]);
	  rawString >> val;
	  vecParam->operator[](j*2)->push_back(val);
	  rawString >> val;
	  vecParam->operator[](j*2+1)->push_back(val);
	  std::cout << "Fit observable-" << j << " systematic error: +" << vecParam->operator[](j*2)->back() << "/-" << vecParam->operator[](j*2+1)->back() << std::endl;
	}
    }

  
  ParVector.clear();
  delete ParameterFile;
}

void Utils::ReadNLL (std::string fileName, std::vector<std::vector<double>*>* vecParam)
// vecParam[0] --> Fl
// vecParam[1] --> Afb
// vecParam[2] --> At2
// vecParam[3] --> Atim
// vecParam[4] --> Branching-Fraction
{
  double val;
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // #########################################
  // # Read fit-observable systematic errors #
  // #########################################
  ParameterFile->ReadFromFile(ParFileBlockN("fitNLL"),&ParVector);

  for (unsigned int j = 0; j < nFitObserv*2; j++) vecParam->push_back(new std::vector<double>);
  
  for (unsigned int i = 0; i < ParVector.size(); i=i+nFitObserv)
    {

      std::cout << "\nRead set-" << (int)((double)i/(double)(nFitObserv)) << " of fit-observable NLL" << std::endl;

      for (unsigned int j = 0; j < nFitObserv; j++)
	{
	  std::stringstream rawString(ParVector[i+j]);
	  rawString >> val;
	  vecParam->operator[](j)->push_back(val);
	  std::cout << "Fit observable-" << j << " NLL: " << vecParam->operator[](j)->back() << std::endl;
	}
    }

  
  ParVector.clear();
  delete ParameterFile;
}

double Utils::ReadLumi (std::string fileName)
{
  double val = 0.0;
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ##############################
  // # Read integrated luminosity #
  // ##############################
  ParameterFile->ReadFromFile(ParFileBlockN("lumi"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++) val = val + atof(ParVector[i].c_str());

  std::cout << "\n@@@ Recorded luminosity: " << val << " fb-1 @@@" << std::endl;

  ParVector.clear();
  delete ParameterFile;
  return val;
}

void Utils::SaveAnalyticalEff (std::string fileName, TF2* effFunc, double q2Val, std::vector<double>* q2Bins)
{
  ofstream fileOutput;
  fileOutput.open(fileName.c_str(),ofstream::app);
  if (fileOutput.good() == false)
    {
      std::cout << "[Utils::SaveAnalyticalEff]\tError opening file : " << fileName << std::endl;
      exit (1);
    }
  
  for (unsigned int k = 0; k < NcoeffThetaL; k++)
    {
      fileOutput << q2Val;
      for (unsigned int j = 0; j < NcoeffThetaK; j++) fileOutput << "   " << effFunc->GetParameter(j+k*NcoeffThetaK) << "   " << effFunc->GetParError(j+k*NcoeffThetaK);
      fileOutput << std::endl;
    }

  fileOutput.close();
}

void Utils::SaveAnalyticalEff (std::string fileName, double q2Val, std::vector<double>* q2Bins, std::vector<double>* effFuncPar, std::vector<double>* effFuncErr)
{
  ofstream fileOutput;
  fileOutput.open(fileName.c_str(),ofstream::app);
  if (fileOutput.good() == false)
    {
      std::cout << "[Utils::SaveAnalyticalEff]\tError opening file : " << fileName << std::endl;
      exit (1);
    }
  
  for (unsigned int k = 0; k < NcoeffThetaL; k++)
    {
      fileOutput << q2Val;
      for (unsigned int j = 0; j < NcoeffThetaK; j++) fileOutput << "   " << effFuncPar->operator[](j+k*NcoeffThetaK) << "   " << (effFuncErr != NULL ? effFuncErr->operator[](j+k*NcoeffThetaK) : 0.0);
      fileOutput << std::endl;
    }

  fileOutput.close();
}

void Utils::SaveAnalyticalEffFullCovariance (std::string fileName, TMatrixTSym<double>* covMatrix, double q2Val, std::vector<double>* q2Bins)
{
  ofstream fileOutput;
  fileOutput.open(fileName.c_str(),ofstream::app);
  if (fileOutput.good() == false)
    {
      std::cout << "[Utils::SaveAnalyticalEffFullCovariance]\tError opening file : " << fileName << std::endl;
      exit (1);
    }
  
  for (int i = 0; i < covMatrix->GetNrows(); i++)
    {
      fileOutput << q2Val;
      for (int j = 0; j < covMatrix->GetNrows(); j++) fileOutput << "   " << (*covMatrix)(i,j);
      fileOutput << std::endl;
    }

  fileOutput.close();
}

std::string Utils::TellMeEffFuncThetaKThetaL ()
{
  return "([0]+[1]*x+[2]*x*x+[3]*x*x*x) + ([4]+[5]*x+[6]*x*x+[7]*x*x*x)*y*y + ([8]+[9]*x+[10]*x*x+[11]*x*x*x)*y*y*y + ([12]+[13]*x+[14]*x*x+[15]*x*x*x)*y*y*y*y + ([16]+[17]*x+[18]*x*x+[19]*x*x*x)*y*y*y*y*y*y";
}

std::string Utils::TellMeEffFuncThetaK ()
{
  return "[0] + [1]*x + [2]*x*x + [3]*x*x*x";
}

std::string Utils::TellMeEffFuncThetaL ()
{
  return "[0]  + [1]*x*x + [2]*x*x*x + [3]*x*x*x*x + [4]*x*x*x*x*x*x";
}

void Utils::ReadAnalyticalEff (std::string fileNameEffParams,
			       std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins,
			       std::vector<TF2*>* effFuncs, std::string effID, const unsigned int dataBlockN)
{
  unsigned int indx;
  std::stringstream myString;
  double* coeffVec = new double[NcoeffThetaK];

  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileNameEffParams.c_str());


  // ######################################
  // # Read block of parameters from file #
  // ######################################
  ParameterFile->ReadFromFile(dataBlockN,&ParVector);

  for (unsigned int q2BinIndx = 0; q2BinIndx < q2Bins->size()-1; q2BinIndx++)
    {
      myString.str("");
      myString << effID << "_" << q2BinIndx;
      effFuncs->push_back(new TF2(myString.str().c_str(),TellMeEffFuncThetaKThetaL().c_str(),
				  cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),
				  cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1)));
      
      std::cout << "\n@@@ Reading coefficients for analytical efficiency for set-" << q2BinIndx << " from file " << fileNameEffParams.c_str() << " @@@" << std::endl;
      
      for (unsigned int k = 0; k < NcoeffThetaL; k++)
	{
	  std::stringstream rawStringK(ParVector[k+q2BinIndx*NcoeffThetaL]);
	  rawStringK >> coeffVec[0]; // Discard q2 bin value
	  indx = 0;
	  for (unsigned int j = 0; j < NcoeffThetaK*2; j = j+2)
	    {
	      rawStringK >> coeffVec[indx];
	      effFuncs->operator[](q2BinIndx)->SetParameter(indx+NcoeffThetaK*k,coeffVec[indx]);
	      rawStringK >> coeffVec[indx];
	      effFuncs->operator[](q2BinIndx)->SetParError(indx+NcoeffThetaK*k,coeffVec[indx]);

	      std::cout << "Theta_L coef. " << k << " --> reading coef. " << indx << " for var. theta_K: ";
	      std::cout << effFuncs->operator[](q2BinIndx)->GetParameter(indx+NcoeffThetaK*k) << " +/- ";
	      std::cout << effFuncs->operator[](q2BinIndx)->GetParError(indx+NcoeffThetaK*k) << std::endl;

	      indx++;
	    }
	}
        
      effFuncs->operator[](q2BinIndx)->GetXaxis()->SetTitle("cos(#theta_{K})");
      effFuncs->operator[](q2BinIndx)->GetXaxis()->SetTitleOffset(1.8);
      effFuncs->operator[](q2BinIndx)->GetYaxis()->SetTitle("cos(#theta_{l})");
      effFuncs->operator[](q2BinIndx)->GetYaxis()->SetTitleOffset(1.8);
      effFuncs->operator[](q2BinIndx)->GetZaxis()->SetTitle("Efficiency");
    }

  delete coeffVec;
  ParVector.clear();
  delete ParameterFile;
}

void Utils::ReadAnalyticalEffFullCovariance (std::string fileNameEffParams, std::vector<TMatrixTSym<double>*>* covMatrices, const unsigned int dataBlockN)
// #######################
// # x <--> cos(theta_K) #
// # y <--> cos(theta_l) #
// #######################
{
  const unsigned int Ncoeff = NcoeffThetaL*NcoeffThetaK;

  double* coeffVec = new double[Ncoeff];

  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileNameEffParams.c_str());


  // ######################################
  // # Read block of parameters from file #
  // ######################################
  ParameterFile->ReadFromFile(dataBlockN,&ParVector);

  for (int q2BinIndx = 0; q2BinIndx < int(rint(ParVector.size()/Ncoeff)); q2BinIndx++)
    {
      std::cout << "\n@@@ Reading covariance matrix for analytical efficiency for set-" << q2BinIndx << " from file " << fileNameEffParams.c_str() << " @@@" << std::endl;
      
      covMatrices->push_back(new TMatrixTSym<double>(Ncoeff));
      
      for (unsigned int j = 0; j < Ncoeff; j++)
	{
	  std::stringstream rawString(ParVector[j+q2BinIndx*Ncoeff]);
	  std::cout << "\nRow #" << j << ": " << rawString.str().c_str() << std::endl;
	  rawString >> coeffVec[0]; // Discard q2 bin value

	  for (unsigned int k = 0; k < Ncoeff; k++)
	    {
	      rawString >> coeffVec[k];
	      (*covMatrices->operator[](q2BinIndx))[j][k] = coeffVec[k];
	      std::cout << "Covariance [" << j << "][" << k << "] --> " << coeffVec[k] << std::endl;
	    }
	}
    }

  delete coeffVec;
  ParVector.clear();
  delete ParameterFile;
}

double Utils::EffMinValue1D (double minX, double maxX, TF1* effFunc)
{
  const unsigned int nsteps = 2000;
  unsigned int iMem = 0;
  double minVal = 0.0;

  // ###################################################
  // # Search for the minimal value of the 1D-function #
  // ###################################################
  for (unsigned int i = 0; i <= nsteps; i++)
    if (effFunc->Eval(maxX - (maxX - minX) / ((double)nsteps) * ((double)i)) < minVal)
      {
	minVal = effFunc->Eval(maxX - (maxX - minX) / ((double)nsteps) * ((double)i));
	iMem = i;
      }
  
  std::cout << "\n@@@ Efficiency minimal value (if less than zero): " << minVal << " at: X = " << maxX - (maxX - minX) / ((double)nsteps) * ((double)iMem);
  std::cout << " (step along X: " << (maxX - minX) / ((double)nsteps) << ") @@@" << std::endl;

  return minVal;
}

double Utils::EffMinValue2D (std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, TF2* effFunc)
{
  const unsigned int nsteps = 2000;
  unsigned int iMem = 0;
  unsigned int jMem = 0;
  double minVal = 0.0;
  
  // ###############################################################################################################################################
  // # Search for the minimal value of the efficiency scanning the domain lattice divided in nsteps per coordinate, i.e. nsteps*nsteps grid matrix #
  // ###############################################################################################################################################
  for (unsigned int i = 0; i <= nsteps; i++)
    for (unsigned int j = 0; j <= nsteps; j++)
      if (effFunc->Eval(cosThetaKBins->operator[](cosThetaKBins->size()-1) - (cosThetaKBins->operator[](cosThetaKBins->size()-1) - cosThetaKBins->operator[](0)) / ((double)nsteps) * (double)j,
			cosThetaLBins->operator[](cosThetaLBins->size()-1) - (cosThetaLBins->operator[](cosThetaLBins->size()-1) - cosThetaLBins->operator[](0)) / ((double)nsteps) * (double)i) < minVal)
	{
	  minVal = effFunc->Eval(cosThetaKBins->operator[](cosThetaKBins->size()-1) - (cosThetaKBins->operator[](cosThetaKBins->size()-1) - cosThetaKBins->operator[](0)) / ((double)nsteps) * (double)j,
				 cosThetaLBins->operator[](cosThetaLBins->size()-1) - (cosThetaLBins->operator[](cosThetaLBins->size()-1) - cosThetaLBins->operator[](0)) / ((double)nsteps) * (double)i);
	  iMem = i;
	  jMem = j;
	}
  
  std::cout << "\n@@@ Efficiency minimal value (if less than zero): " << minVal << " at: (theta_K,theta_l = ";
  std::cout << cosThetaKBins->operator[](cosThetaKBins->size()-1) - (cosThetaKBins->operator[](cosThetaKBins->size()-1) - cosThetaKBins->operator[](0)) / ((double)nsteps) * (double)jMem << ",";
  std::cout << cosThetaLBins->operator[](cosThetaLBins->size()-1) - (cosThetaLBins->operator[](cosThetaLBins->size()-1) - cosThetaLBins->operator[](0)) / ((double)nsteps) * (double)iMem << ")";
  std::cout << " (grid step along cos(theta_K): ";
  std::cout << (cosThetaKBins->operator[](cosThetaKBins->size()-1) - cosThetaKBins->operator[](0)) / ((double)nsteps);
  std::cout << "; grid step along cos(theta_l): ";
  std::cout << (cosThetaLBins->operator[](cosThetaLBins->size()-1) - cosThetaLBins->operator[](0)) / ((double)nsteps);
  std::cout << ") @@@" << std::endl;

  return minVal;
}

void Utils::MakeGraphVar (std::string parFileName, TGraphAsymmErrors** graph, std::string varName, bool allBins)
// ######################
// # varName = "dBFdq2" #
// # varName = "Fl"     #
// # varName = "Afb"    #
// # varName = "At2"    #
// # varName = "Atim"   #
// ######################
{
  double tmpVar;

  std::vector<std::vector<std::string>*> vecParam;
  std::vector<std::vector<unsigned int>*> configParam;
  std::vector<std::string> ParVector;
  std::vector<double> vxs;
  std::vector<double> vys;
  std::vector<double> vxel;
  std::vector<double> vxeh;
  std::vector<double> vyel;
  std::vector<double> vyeh;
  

  // ################################
  // # Read values from config file #
  // ################################
  std::vector<double> q2Bins;
  Readq2Bins(parFileName,&q2Bins);
  ReadFitStartingValues(parFileName,&vecParam,&configParam,ParFileBlockN("fitValBins"));
  ReadParameters* ParameterFile = new ReadParameters(parFileName.c_str());
  ParameterFile->ReadFromFile(ParFileBlockN("dBFdq2"),&ParVector);


  for (unsigned int i = 0; i < q2Bins.size()-1; i++)
    {
      std::stringstream rawString;
      if      (varName == "dBFdq2") rawString << ParVector[i];
      else if (varName == "Fl")     rawString << vecParam[GetFitParamIndx("FlS")]->operator[](i);
      else if (varName == "Afb")    rawString << vecParam[GetFitParamIndx("AfbS")]->operator[](i);
      else if (varName == "At2")    rawString << vecParam[GetFitParamIndx("At2S")]->operator[](i);
      else if (varName == "Atim")   rawString << vecParam[GetFitParamIndx("AtimS")]->operator[](i);
      else { std::cout << "[Utils::MakeGraphVar]\tVariable name unknown: " << varName << std::endl; exit(1); }

      if ((allBins == true) || (ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false))
	{
	  vxs.push_back((q2Bins[i+1] + q2Bins[i]) / 2.);
	  vxeh.push_back(q2Bins[i+1] - vxs.back());
	  vxel.push_back(vxs.back() - q2Bins[i]);
 
	  rawString >> tmpVar;
	  vys.push_back(tmpVar);
	  rawString >> tmpVar;
	  vyel.push_back(fabs(tmpVar));
	  rawString >> tmpVar;
	  vyeh.push_back(fabs(tmpVar));
	}
      else
	{
	  vxs.push_back((q2Bins[i+1] + q2Bins[i]) / 2.);
	  vxeh.push_back(0.0);
	  vxel.push_back(0.0);

	  vys.push_back(YvalueOutsideLimits);
	  vyeh.push_back(0.0);
	  vyel.push_back(0.0);
	}
    }
  *graph = new TGraphAsymmErrors(vxs.size(), &vxs[0], &vys[0], &vxel[0], &vxeh[0], &vyel[0], &vyeh[0]);


  // #################
  // # Clear vectors #
  // #################
  q2Bins.clear();
  ParVector.clear();
  for (unsigned int i = 0; i < vecParam.size(); i++) vecParam[i]->clear();
  vecParam.clear();
  vxs.clear();
  vys.clear();
  vxel.clear();
  vxeh.clear();
  vyel.clear();
  vyeh.clear();
  delete ParameterFile;
}

void Utils::InitEffFuncThetaL (TF1* fitFun, unsigned int q2BinIndx)
{
  fitFun->ReleaseParameter(0);
  fitFun->ReleaseParameter(1);
  fitFun->ReleaseParameter(2);
  fitFun->ReleaseParameter(3);
  fitFun->ReleaseParameter(4);

  fitFun->SetParameter(0,1e-3);
  fitFun->SetParameter(1,-1.0);

  if (q2BinIndx == 0)
    {
      fitFun->FixParameter(2,0.0);
      fitFun->SetParameter(3,0.0);
      fitFun->SetParameter(4,0.0);
    }
  else if (q2BinIndx == 1)
    {
      fitFun->FixParameter(2,0.0);
      fitFun->SetParameter(3,0.0);
      fitFun->FixParameter(4,0.0);
    }
  else if (q2BinIndx == 2)
    {
      fitFun->SetParameter(2,0.0);
      fitFun->SetParameter(3,0.0);
      fitFun->FixParameter(4,0.0);
    }
  else if (q2BinIndx == 3)
    {
      fitFun->SetParameter(2,0.0);
      fitFun->SetParameter(3,0.0);
      fitFun->FixParameter(4,0.0);
    }
  else if (q2BinIndx == 4)
    {
      fitFun->SetParameter(2,0.0);
      fitFun->SetParameter(3,0.0);
      fitFun->FixParameter(4,0.0);
    }
  else if (q2BinIndx == 5)
    {
      fitFun->SetParameter(2,0.0);
      fitFun->SetParameter(3,0.0);
      fitFun->FixParameter(4,0.0);
    }
  else if (q2BinIndx == 6)
    {
      fitFun->SetParameter(2,0.0);
      fitFun->FixParameter(3,0.0);
      fitFun->FixParameter(4,0.0);
    }
  else if (q2BinIndx == 7)
    {
      fitFun->SetParameter(2,0.0);
      fitFun->FixParameter(3,0.0);
      fitFun->FixParameter(4,0.0);
    }
}

void Utils::InitEffFuncThetaK (TF1* fitFun, unsigned int q2BinIndx)
{
  fitFun->ReleaseParameter(0);
  fitFun->ReleaseParameter(1);
  fitFun->ReleaseParameter(2);
  fitFun->ReleaseParameter(3);

  fitFun->SetParameter(0,0.0);
  fitFun->SetParameter(1,0.0);
  fitFun->SetParameter(2,0.0);
  fitFun->SetParameter(3,0.0);
}

void Utils::AddConstraint1D (TH1D** histo, std::string constrType, double errX, double Yval, double errY, unsigned int ID)
// ############################################
// # constrType = "low" : at lower boundary   #
// # constrType = "both" : at both boundaries #
// # constrType = "high" : at higher boundary #
// ############################################
{
  std::stringstream myString;

  TH1D* newHisto;
  TH1D* tmpHisto;

  unsigned int nNewBins;
  if ((constrType == "low") || (constrType == "high")) nNewBins = (*histo)->GetNbinsX()+2;
  else if (constrType == "both") nNewBins = (*histo)->GetNbinsX()+3;
  else { std::cout << "[ComputeEfficiency::AddConstraint1D]\tError wrong parameter" << std::endl; exit(1); }
  double* reBins;
  reBins = new double[nNewBins];


  if (constrType == "low")
    {
      reBins[0] = (*histo)->GetBinLowEdge(1) - errX;
      for (int i = 1; i <= (*histo)->GetNbinsX(); i++) reBins[i] = (*histo)->GetBinLowEdge(i);
      reBins[(*histo)->GetNbinsX()+1] = (*histo)->GetBinLowEdge((*histo)->GetNbinsX()) + (*histo)->GetBinWidth((*histo)->GetNbinsX());
    }
  else if (constrType == "high")
    {
      for (int i = 1; i <= (*histo)->GetNbinsX(); i++) reBins[i-1] = (*histo)->GetBinLowEdge(i);
      reBins[(*histo)->GetNbinsX()] = (*histo)->GetBinLowEdge((*histo)->GetNbinsX()) + (*histo)->GetBinWidth((*histo)->GetNbinsX());
      reBins[(*histo)->GetNbinsX()+1] = reBins[(*histo)->GetNbinsX()] + errX;
    }
  else
    {
      reBins[0] = (*histo)->GetBinLowEdge(1) - errX;
      for (int i = 1; i <= (*histo)->GetNbinsX(); i++) reBins[i] = (*histo)->GetBinLowEdge(i);
      reBins[(*histo)->GetNbinsX()+1] = (*histo)->GetBinLowEdge((*histo)->GetNbinsX()) + (*histo)->GetBinWidth((*histo)->GetNbinsX());
      reBins[(*histo)->GetNbinsX()+2] = reBins[(*histo)->GetNbinsX()+1] + errX;
    }


  myString.str("");
  myString << (*histo)->GetName() << "_" << ID;
  newHisto = new TH1D(myString.str().c_str(),myString.str().c_str(),nNewBins-1,reBins);


  if (constrType == "low")
    {
      newHisto->SetBinContent(1, Yval);
      newHisto->SetBinError(1,errY);
      for (int i = 1; i <= (*histo)->GetNbinsX(); i++)
	{
	  newHisto->SetBinContent(i+1,(*histo)->GetBinContent(i));
	  newHisto->SetBinError(i+1,(*histo)->GetBinError(i));
	}
    }
  else if (constrType == "high")
    {
      for (int i = 1; i <= (*histo)->GetNbinsX(); i++)
	{
	  newHisto->SetBinContent(i,(*histo)->GetBinContent(i));
	  newHisto->SetBinError(i,(*histo)->GetBinError(i));
	}
      newHisto->SetBinContent((*histo)->GetNbinsX()+1,Yval);
      newHisto->SetBinError((*histo)->GetNbinsX()+1,errY);
    }
  else
    {
      newHisto->SetBinContent(1,Yval);
      newHisto->SetBinError(1,errY);
      for (int i = 1; i <= (*histo)->GetNbinsX(); i++)
	{
	  newHisto->SetBinContent(i+1,(*histo)->GetBinContent(i));
	  newHisto->SetBinError(i+1,(*histo)->GetBinError(i));
	}
      newHisto->SetBinContent((*histo)->GetNbinsX()+2,Yval);
      newHisto->SetBinError((*histo)->GetNbinsX()+2,errY);
    }


  tmpHisto = *histo;
  *histo = newHisto;
  delete tmpHisto;
  delete reBins;
}

void Utils::AddConstraintThetaL (TH1D** histo,
				 unsigned int q2BinIndx, unsigned int cosThetaKBinIndx,
				 double constrXerr, double constrYval, double constrYerr, unsigned int ID)
{
  if (q2BinIndx == 0) AddConstraint1D(histo,"both",constrXerr,constrYval,constrYerr,ID);

  if (q2BinIndx == 1) AddConstraint1D(histo,"both",constrXerr,constrYval,constrYerr,ID);

  if ((q2BinIndx == 2) && (cosThetaKBinIndx == 0)) AddConstraint1D(histo,"high",constrXerr,constrYval,constrYerr,ID);
  if ((q2BinIndx == 2) && (cosThetaKBinIndx == 1)) AddConstraint1D(histo,"high",constrXerr,constrYval,constrYerr,ID);
  if ((q2BinIndx == 2) && (cosThetaKBinIndx == 2)) AddConstraint1D(histo,"high",constrXerr,constrYval,constrYerr,ID);
  if ((q2BinIndx == 2) && (cosThetaKBinIndx == 3)) AddConstraint1D(histo,"low",constrXerr,constrYval,constrYerr,ID);
}

void Utils::AddConstraint2D (TH2D** histo, double err, double Zval, double errZ, unsigned int ID, std::string toBeConstr, std::vector<std::string>* toBeAdded)
// ###############################################################################################################################################
// # toBeConstr = Y --> Add constraints to Y axis (= cosThetaL) to both sides, according to the toBeAdded variable (= X axis bining = cosThetaK) #
// # toBeConstr = X --> Add constraints to X axis (= cosThetaK) to the whole positive side                                                       #
// ###############################################################################################################################################
// # toBeAdded = "low" : at lower boundary                                                                                                       #
// # toBeAdded = "both" : at both boundaries                                                                                                     #
// # toBeAdded = "high" : at higher boundary                                                                                                     #
// ###############################################################################################################################################
{
  std::stringstream myString;
  
  double* reBinsX = NULL;
  double* reBinsY = NULL;

  TH2D* newHisto = NULL;
  TH2D* tmpHisto;

  TAxis* XAxis;
  TAxis* YAxis;


  if (toBeConstr == "Y")
    {
      unsigned int nNewBinsX;
      nNewBinsX = (*histo)->GetNbinsX()+1;
      reBinsX = new double[nNewBinsX];
      XAxis = (*histo)->GetXaxis();
      
      for (int i = 1; i <= (*histo)->GetNbinsX(); i++) reBinsX[i-1] = XAxis->GetBinLowEdge(i);
      reBinsX[(*histo)->GetNbinsX()] = XAxis->GetBinLowEdge((*histo)->GetNbinsX()) + XAxis->GetBinWidth((*histo)->GetNbinsX());
      
      
      unsigned int nNewBinsY;
      nNewBinsY = (*histo)->GetNbinsY()+3;
      reBinsY = new double[nNewBinsY];
      YAxis = (*histo)->GetYaxis();
      
      reBinsY[0] = YAxis->GetBinLowEdge(1) - err;
      for (int i = 1; i <= (*histo)->GetNbinsY(); i++) reBinsY[i] = YAxis->GetBinLowEdge(i);
      reBinsY[(*histo)->GetNbinsY()+1] = YAxis->GetBinLowEdge((*histo)->GetNbinsY()) + YAxis->GetBinWidth((*histo)->GetNbinsY());
      reBinsY[(*histo)->GetNbinsY()+2] = reBinsY[(*histo)->GetNbinsY()+1] + err;
      
      
      myString.str("");
      myString << (*histo)->GetName() << "_" << ID;
      newHisto = new TH2D(myString.str().c_str(),myString.str().c_str(),nNewBinsX-1,reBinsX,nNewBinsY-1,reBinsY);
      
      
      for (unsigned int i = 0; i < toBeAdded->size(); i++)
	{
	  if ((toBeAdded->operator[](i) == "low") || (toBeAdded->operator[](i) == "both"))
	    {
	      newHisto->SetBinContent(i+1,1,Zval);
	      newHisto->SetBinError(i+1,1,errZ);
	    }
	  else
	    {
	      newHisto->SetBinContent(i+1,1,0.0);
	      newHisto->SetBinError(i+1,1,0.0);
	    }
	  
	  for (int j = 1; j <= (*histo)->GetNbinsY(); j++)
	    {
	      newHisto->SetBinContent(i+1,j+1,(*histo)->GetBinContent(i+1,j));
	      newHisto->SetBinError(i+1,j+1,(*histo)->GetBinError(i+1,j));
	    }
	  
	  if ((toBeAdded->operator[](i) == "high") || (toBeAdded->operator[](i) == "both"))
	    {
	      newHisto->SetBinContent(i+1,(*histo)->GetNbinsY()+2,Zval);
	      newHisto->SetBinError(i+1,(*histo)->GetNbinsY()+2,errZ);
	    }
	  else
	    {
	      newHisto->SetBinContent(i+1,(*histo)->GetNbinsY()+2,0.0);
	      newHisto->SetBinError(i+1,(*histo)->GetNbinsY()+2,0.0);
	    }
	}
    }
  else if (toBeConstr == "X")
    {
      unsigned int nNewBinsX;
      nNewBinsX = (*histo)->GetNbinsX()+2;
      reBinsX = new double[nNewBinsX];
      XAxis = (*histo)->GetXaxis();
      
      for (int i = 1; i <= (*histo)->GetNbinsX(); i++) reBinsX[i-1] = XAxis->GetBinLowEdge(i);
      reBinsX[(*histo)->GetNbinsX()] = XAxis->GetBinLowEdge((*histo)->GetNbinsX()) + XAxis->GetBinWidth((*histo)->GetNbinsX());
      reBinsX[(*histo)->GetNbinsX()+1] = reBinsX[(*histo)->GetNbinsX()] + err;

      
      unsigned int nNewBinsY;
      nNewBinsY = (*histo)->GetNbinsY()+1;
      reBinsY = new double[nNewBinsY];
      YAxis = (*histo)->GetYaxis();
      
      for (int i = 1; i <= (*histo)->GetNbinsY(); i++) reBinsY[i-1] = YAxis->GetBinLowEdge(i);
      reBinsY[(*histo)->GetNbinsY()] = YAxis->GetBinLowEdge((*histo)->GetNbinsY()) + YAxis->GetBinWidth((*histo)->GetNbinsY());

      
      myString.str("");
      myString << (*histo)->GetName() << "_" << ID;
      newHisto = new TH2D(myString.str().c_str(),myString.str().c_str(),nNewBinsX-1,reBinsX,nNewBinsY-1,reBinsY);

      for (int i = 1; i <= (*histo)->GetNbinsX(); i++)
	for (int j = 1; j <= (*histo)->GetNbinsY(); j++)
	  {
	    newHisto->SetBinContent(i,j,(*histo)->GetBinContent(i,j));
	    newHisto->SetBinError(i,j,(*histo)->GetBinError(i,j));
	  }      
      newHisto->SetBinContent((*histo)->GetNbinsX()+1,1,Zval);
      newHisto->SetBinError((*histo)->GetNbinsX()+1,1,errZ);
      for (int j = 2; j < (*histo)->GetNbinsY(); j++)
	{
	  newHisto->SetBinContent((*histo)->GetNbinsX()+1,j,Zval);
	  newHisto->SetBinError((*histo)->GetNbinsX()+1,j,errZ);
	}
      newHisto->SetBinContent((*histo)->GetNbinsX()+1,(*histo)->GetNbinsY(),Zval);
      newHisto->SetBinError((*histo)->GetNbinsX()+1,(*histo)->GetNbinsY(),errZ);
    }


  newHisto->SetXTitle("cos(#theta_{K})");
  newHisto->GetXaxis()->SetTitleOffset(1.8);
  newHisto->SetYTitle("cos(#theta_{l})");
  newHisto->GetYaxis()->SetTitleOffset(1.8);
  newHisto->SetZTitle("Efficiency");


  tmpHisto = *histo;
  *histo = newHisto;
  delete tmpHisto;
  delete reBinsX;
  delete reBinsY;
}

void Utils::AddConstraintThetaK (TH2D** histo, std::vector<double>* cosThetaKBins,
				 unsigned int q2BinIndx, double constrXYerr, double constrZval, double constrZerr, unsigned int ID)
{
  std::vector<std::string> toBeAdded;

  for (unsigned int i = 0; i < cosThetaKBins->size()-1; i++) toBeAdded.push_back("no");

  if (q2BinIndx == 0) for (unsigned int i = 0; i < cosThetaKBins->size()-1; i++) toBeAdded[i] = "both";

  if (q2BinIndx == 1) for (unsigned int i = 0; i < cosThetaKBins->size()-1; i++) toBeAdded[i] = "both";

  if (q2BinIndx == 2) toBeAdded[0] = "high";
  if (q2BinIndx == 2) toBeAdded[1] = "high";
  if (q2BinIndx == 2) toBeAdded[2] = "high";
  if (q2BinIndx == 2) toBeAdded[3] = "low";

  AddConstraint2D(histo,constrXYerr,constrZval,constrZerr,ID,"Y",&toBeAdded);

  toBeAdded.clear();
}

bool Utils::IsThereOddDegree (TF2* effFunc)
// ##############################################################################
// # Tells if in the efficiency function there is the odd-power term vs theta_l #
// ##############################################################################
{
  if ((effFunc->GetParameter (8) == 0.0) &&
      (effFunc->GetParameter (9) == 0.0) &&
      (effFunc->GetParameter(10) == 0.0) &&
      (effFunc->GetParameter(11) == 0.0) &&

      (effFunc->GetParError (8) == 0.0) &&
      (effFunc->GetParError (9) == 0.0) &&
      (effFunc->GetParError(10) == 0.0) &&
      (effFunc->GetParError(11) == 0.0)) return false;
  
  return true;
}

bool Utils::IsThisData (std::string fileName)
{
  bool val = true;
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ###############################
  // # Read data type (Data or MC) #
  // ###############################
  ParameterFile->ReadFromFile(ParFileBlockN("dtype"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++) val = val * atoi(ParVector[i].c_str());


  ParVector.clear();
  delete ParameterFile;
  return val;
}

void Utils::SaveFitValues (std::string fileName, std::vector<std::string>* vecParStr, int indx, std::string str)
// ##############################
// # If indx == -1 then use str #
// ##############################
{
  std::stringstream myString;
  
  myString.clear(); myString.str("");
  if (indx != -1) myString << "#@@@@@@@@@@@@@@@@@@@@@@@@@@@ q^2 bin " << indx << " @@@@@@@@@@@@@@@@@@@@@@@@@@@#";
  else myString            << "#@@@@@@@@@@@@@@@@@@@@@@@@@@@" << str.c_str()   <<  "@@@@@@@@@@@@@@@@@@@@@@@@@@@#";
  vecParStr->insert(vecParStr->begin(),myString.str());

  vecParStr->insert(vecParStr->end(),"");
  vecParStr->insert(vecParStr->end(),"");

  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str(),"app");
  ParameterFile->SaveToFile(0,vecParStr);
  delete ParameterFile;
}

unsigned int Utils::ParFileBlockN (std::string blockName)
{
  // #########################
  // # Parameter file blocks #
  // #########################
  if      (blockName == "HLTpath")     return 1;
  else if (blockName == "precuts")     return 2;
  else if (blockName == "selecuts")    return 3;
  else if (blockName == "q2")          return 4;
  else if (blockName == "thetaK")      return 5;
  else if (blockName == "thetaL")      return 6;
  else if (blockName == "phi")         return 7;
  else if (blockName == "HLTcuts")     return 8;
  else if (blockName == "fitValGlob")  return 9;
  else if (blockName == "fitValBins")  return 10;
  else if (blockName == "dBFdq2")      return 11;
  else if (blockName == "fitSyst")     return 12;
  else if (blockName == "fitNLL")      return 13;

  else if (blockName == "analyEff1")   return 14;
  else if (blockName == "analyEff2")   return 15;
  else if (blockName == "analyEff3")   return 16;
  else if (blockName == "analyEff4")   return 17;
  else if (blockName == "analyEffAvg") return 18;

  else if (blockName == "genericpar")  return 19;
  else if (blockName == "lumi")        return 20;
  else if (blockName == "dtype")       return 21;

  std::cout << "[Utils::ParFileBlockN]\tError wrong index name : " << blockName << std::endl;
  exit (1);
}

unsigned int Utils::GetFitParamIndx (std::string varName)
{
  if      (varName == "meanS")          return 0;
  else if (varName == "sigmaS1")        return 1;
  else if (varName == "sigmaS2")        return 2;
  else if (varName == "fracMassS")      return 3;

  else if (varName == "tau1")           return 4;
  else if (varName == "tau2")           return 5;
  else if (varName == "fracMassBExp")   return 6;

  else if (varName == "meanR1")         return 7;
  else if (varName == "sigmaR1")        return 8;
  else if (varName == "meanR2")         return 9;
  else if (varName == "sigmaR2")        return 10;
  else if (varName == "fracMassBRPeak") return 11;

  else if (varName == "meanL1")         return 12;
  else if (varName == "sigmaL1")        return 13;
  else if (varName == "meanL2")         return 14;
  else if (varName == "sigmaL2")        return 15;
  else if (varName == "fracMassBLPeak") return 16;

  else if (varName == "fracMassBPeak")  return 17;

  else if (varName == "nBkgComb")       return 18;
  else if (varName == "nBkgPeak")       return 19;
  else if (varName == "nSig")           return 20;

  else if (varName == "nPolyP1")        return 21;
  else if (varName == "p1Poly0")        return 22;
  else if (varName == "p1Poly1")        return 23;
  else if (varName == "p1Poly2")        return 24;
  else if (varName == "p1Poly3")        return 25;
  else if (varName == "p1Poly4")        return 26;

  else if (varName == "nPolyC1")        return 27;
  else if (varName == "c1Poly0")        return 28;
  else if (varName == "c1Poly1")        return 29;
  else if (varName == "c1Poly2")        return 30;
  else if (varName == "c1Poly3")        return 31;
  else if (varName == "c1Poly4")        return 32;

  else if (varName == "nPolyP2")        return 33;
  else if (varName == "p2Poly0")        return 34;
  else if (varName == "p2Poly1")        return 35;
  else if (varName == "p2Poly2")        return 36;
  else if (varName == "p2Poly3")        return 37;
  else if (varName == "p2Poly4")        return 38;

  else if (varName == "nPolyC2")        return 39;
  else if (varName == "c2Poly0")        return 40;
  else if (varName == "c2Poly1")        return 41;
  else if (varName == "c2Poly2")        return 42;
  else if (varName == "c2Poly3")        return 43;
  else if (varName == "c2Poly4")        return 44;

  else if (varName == "FlS")            return 45;
  else if (varName == "AfbS")           return 46;
  else if (varName == "At2S")           return 47;
  else if (varName == "AtimS")          return 48;
  else if (varName == "FsS")            return 49;
  else if (varName == "AsS")            return 50;

  std::cout << "[Utils::GetFitParamIndx]\tError wrong index name : " << varName << std::endl;
  exit (1);
}

unsigned int Utils::GetConfigParamIndx (std::string varName)
{
  if      (varName == "nGaussPSIintru") return 0;
  else if (varName == "fitPSIintru")    return 1;
  else if (varName == "SigType")        return 2;
  else if (varName == "BkgType")        return 3;

  std::cout << "[Utils::GetConfigParamIndx]\tError wrong index name : " << varName << std::endl;
  exit (1);
}

bool Utils::ChooseBestCand (B0KstMuMuTreeContent* NTuple, unsigned int DoTrigCheck, double evFraction, int* BestCandIndx, bool* B0notB0bar, int* TrigCat, unsigned int* countCands)
// ##############################################################################
// # DoTrigCheck: to allow check on trigger requirements for muons:             #
// # 0: do not perform any trigger check                                        #
// # 1: perform normal trigger check                                            #
// # 2: perform trigger check and split the sample in HLT categories            #
// # 3: do NOT perform any trigger check and split the sample in HLT categories #
// ##############################################################################
{
  // #####################################################################
  // # This method is used together with the method: "ReadSelectionCuts" #
  // #####################################################################

  double CLMuMuVtx = GetSeleCut("CLMuMuVtx");
  double MinMupT   = GetSeleCut("MinMupT");
  double BestVal   = 0.0;
  double BestValTmp;

  *countCands = 0;
  *BestCandIndx = -1;
  *TrigCat = 0;
  for (unsigned int i = 0; i < NTuple->bMass->size(); i++)
    {
      // ##############################################
      // # Candidate selection through kinematic cuts #
      // ##############################################
      if (((DoTrigCheck == 0) ||
           ((DoTrigCheck == 1) && (IsInTriggerTable(NTuple, &CLMuMuVtx, &MinMupT, i) >= 1)) ||
	   ((DoTrigCheck == 2) && (IsInTriggerTable(NTuple, &CLMuMuVtx, &MinMupT, i, evFraction) >= 1)) ||
	   ((DoTrigCheck == 3) && (IsInTriggerTable(NTuple, &CLMuMuVtx, &MinMupT, -1, evFraction) >= 1))) &&
	  
	  // #####################
	  // # Choose good muons #
	  // #####################
	  (NTuple->mumNTrkHits->at(i) > int(rint(GetSeleCut("NTrkHits")))) &&
	  (NTuple->mumNPixLayers->at(i) > int(rint(GetSeleCut("NPixLayers")))) &&
	  (NTuple->mumNormChi2->at(i) < GetSeleCut("NormChi2")) &&
	  (fabs(NTuple->mumdxyVtx->at(i)) < GetSeleCut("dxyVtx")) &&
	  (fabs(NTuple->mumdzVtx->at(i)) < GetSeleCut("dzVtx")) &&
	  (NTuple->mumCat->at(i).find("TMOneStationTight") != std::string::npos) &&
	  
	  (NTuple->mupNTrkHits->at(i) > int(rint(GetSeleCut("NTrkHits")))) &&
	  (NTuple->mupNPixLayers->at(i) > int(rint(GetSeleCut("NPixLayers")))) &&
	  (NTuple->mupNormChi2->at(i) < GetSeleCut("NormChi2")) &&
	  (fabs(NTuple->mupdxyVtx->at(i)) < GetSeleCut("dxyVtx")) &&
	  (fabs(NTuple->mupdzVtx->at(i)) < GetSeleCut("dzVtx")) &&
	  (NTuple->mupCat->at(i).find("TMOneStationTight") != std::string::npos) &&
	  
	  // #########################################
	  // # Check that hadron track is NOT a muon #
	  // #########################################
	  (NTuple->kstTrkmMuMatch->at(i).find("TrackerMuonArbitrated") == std::string::npos) &&
	  (NTuple->kstTrkpMuMatch->at(i).find("TrackerMuonArbitrated") == std::string::npos) &&

	  // #####################
	  // # B0 selection cuts #
	  // #####################
	  (NTuple->bLBS->at(i)/NTuple->bLBSE->at(i) > GetSeleCut("B0LsBS")) &&
	  (NTuple->bVtxCL->at(i) > GetSeleCut("B0VtxCL")) &&
	  (NTuple->bCosAlphaBS->at(i) > GetSeleCut("B0cosAlpha")) &&
	  (sqrt(NTuple->bPx->at(i)*NTuple->bPx->at(i) + NTuple->bPy->at(i)*NTuple->bPy->at(i)) > GetSeleCut("B0pT")) &&
	  (fabs(computeEta(NTuple->bPx->at(i),NTuple->bPy->at(i),NTuple->bPz->at(i))) < GetSeleCut("B0Eta")) &&

	  // #######################
	  // # Muon selection cuts #
	  // #######################
	  (sqrt(NTuple->mumPx->at(i)*NTuple->mumPx->at(i) + NTuple->mumPy->at(i)*NTuple->mumPy->at(i)) > MinMupT) &&
	  (sqrt(NTuple->mupPx->at(i)*NTuple->mupPx->at(i) + NTuple->mupPy->at(i)*NTuple->mupPy->at(i)) > MinMupT) &&

	  // #########################
	  // # Dimuon selection cuts #
	  // #########################
	  (NTuple->mumuVtxCL->at(i) > CLMuMuVtx) &&

	  // #########################
	  // # Hadron selection cuts #
	  // #########################
	  (fabs(NTuple->kstTrkmDCABS->at(i)/NTuple->kstTrkmDCABSE->at(i)) > GetSeleCut("HadDCASBS")) &&
	  (fabs(NTuple->kstTrkpDCABS->at(i)/NTuple->kstTrkpDCABSE->at(i)) > GetSeleCut("HadDCASBS")) &&

	  (sqrt(NTuple->kstTrkmPx->at(i)*NTuple->kstTrkmPx->at(i) + NTuple->kstTrkmPy->at(i)*NTuple->kstTrkmPy->at(i)) > GetSeleCut("HadpT")) &&
	  (sqrt(NTuple->kstTrkpPx->at(i)*NTuple->kstTrkpPx->at(i) + NTuple->kstTrkpPy->at(i)*NTuple->kstTrkpPy->at(i)) > GetSeleCut("HadpT")) &&

	  // #####################
	  // # K* selection cuts #
	  // #####################
	  ((fabs(NTuple->kstMass->at(i) - kstMass) < GetSeleCut("KstMass")) || (fabs(NTuple->kstBarMass->at(i) - kstMass) < GetSeleCut("KstMass"))) &&

	  // #####################
	  // # phi rejection cut #
	  // #####################
	  (computeInvMass(NTuple->kstTrkmPx->at(i),NTuple->kstTrkmPy->at(i),NTuple->kstTrkmPz->at(i),kaonMass,
			  NTuple->kstTrkpPx->at(i),NTuple->kstTrkpPy->at(i),NTuple->kstTrkpPz->at(i),kaonMass) > GetSeleCut("KKMass")))
	{
	  (*countCands)++;

	  BestValTmp = NTuple->bVtxCL->at(i);
	  if (((fabs(NTuple->kstMass->at(i) - kstMass) < 3.0*kstSigma) || (fabs(NTuple->kstBarMass->at(i) - kstMass) < 3.0*kstSigma)) &&
	      ((fabs(NTuple->kstMass->at(i) - kstMass) > 1.0*kstSigma) || (fabs(NTuple->kstBarMass->at(i) - kstMass) > 1.0*kstSigma)) &&
	      (BestValTmp > BestVal))
	    {
	      if      (DoTrigCheck == 1) *TrigCat = IsInTriggerTable(NTuple, &CLMuMuVtx, &MinMupT, i);
              else if (DoTrigCheck == 2) *TrigCat = IsInTriggerTable(NTuple, &CLMuMuVtx, &MinMupT, i, evFraction);
              else if (DoTrigCheck == 3) *TrigCat = IsInTriggerTable(NTuple, &CLMuMuVtx, &MinMupT, -1, evFraction);

	      BestVal       = BestValTmp;
	      *BestCandIndx = i;
	    }
	}
    }


  if ((*BestCandIndx != -1) && (FlavorTagger(NTuple, *BestCandIndx, B0notB0bar) == true)) return true;
  return false;
}

bool Utils::FlavorTagger (B0KstMuMuTreeContent* NTuple, unsigned int i, bool* B0notB0bar)
{
  // #####################################################################
  // # Computation of the probability related to the K*0 mass hypothesis #
  // #####################################################################
  double distKst    = fabs(KstMassShape->Integral(kstMass,NTuple->kstMass->at(i)));
  double distKstBar = fabs(KstMassShape->Integral(kstMass,NTuple->kstBarMass->at(i)));
  if (distKst < distKstBar) *B0notB0bar = true;
  else                      *B0notB0bar = false;


  // # ###########################################
  // # Scramble the tagging of the CP-eigenstate #
  // # ###########################################
  double rndDice = ((double)rand()) / ((double)RAND_MAX);
  if (rndDice < scrambleFraction)
    {
      rndDice = ((double)rand()) / ((double)RAND_MAX);
      if (rndDice < 0.5)
	{
	  if (*B0notB0bar == false)
	    {
	      // #######################
	      // # Swap the K* - K*bar #
	      // #######################

	      distKst = NTuple->kstMass->at(i);
	      NTuple->kstMass->at(i) = NTuple->kstBarMass->at(i);
	      NTuple->kstBarMass->at(i) = distKst;

	      distKst = NTuple->kstTrkpPx->at(i);
	      NTuple->kstTrkpPx->at(i) = NTuple->kstTrkmPx->at(i);
	      NTuple->kstTrkmPx->at(i) = distKst;

	      distKst = NTuple->kstTrkpPy->at(i);
	      NTuple->kstTrkpPy->at(i) = NTuple->kstTrkmPy->at(i);
	      NTuple->kstTrkmPy->at(i) = distKst;

	      distKst = NTuple->kstTrkpPz->at(i);
	      NTuple->kstTrkpPz->at(i) = NTuple->kstTrkmPz->at(i);
	      NTuple->kstTrkmPz->at(i) = distKst;
	    }
	  *B0notB0bar = true;
	}
      else
	{
	  if (*B0notB0bar == true)
	    {
	      // #######################
	      // # Swap the K* - K*bar #
	      // #######################

	      distKst = NTuple->kstMass->at(i);
	      NTuple->kstMass->at(i) = NTuple->kstBarMass->at(i);
	      NTuple->kstBarMass->at(i) = distKst;

	      distKst = NTuple->kstTrkpPx->at(i);
	      NTuple->kstTrkpPx->at(i) = NTuple->kstTrkmPx->at(i);
	      NTuple->kstTrkmPx->at(i) = distKst;

	      distKst = NTuple->kstTrkpPy->at(i);
	      NTuple->kstTrkpPy->at(i) = NTuple->kstTrkmPy->at(i);
	      NTuple->kstTrkmPy->at(i) = distKst;

	      distKst = NTuple->kstTrkpPz->at(i);
	      NTuple->kstTrkpPz->at(i) = NTuple->kstTrkmPz->at(i);
	      NTuple->kstTrkmPz->at(i) = distKst;
	    }
	  *B0notB0bar = false;
	}
    }


  // ####################################################################
  // # Compute the Fisher-test to see whether the (mK*0 - K*0 PDG mass) #
  // # is significantly different from (mK*0bar - K*0 PDG mass)         #
  // ####################################################################
  double VarianceRatio;
  if (pow((NTuple->kstMass->at(i)-kstMass) / NTuple->kstMassE->at(i),2.) > pow((NTuple->kstBarMass->at(i)-kstMass) / NTuple->kstBarMassE->at(i),2.))
    VarianceRatio = pow((NTuple->kstMass->at(i)-kstMass) / NTuple->kstMassE->at(i),2.) / pow((NTuple->kstBarMass->at(i)-kstMass)  / NTuple->kstBarMassE->at(i),2.);
  else VarianceRatio = pow((NTuple->kstBarMass->at(i)-kstMass) / NTuple->kstBarMassE->at(i),2.) / pow((NTuple->kstMass->at(i)-kstMass) / NTuple->kstMassE->at(i),2.);

  if (TMath::FDistI(VarianceRatio,1.,1.) > ProbThreshold) return true;
  return false;
}

void Utils::ReadSelectionCuts (std::string fileName)
// #############################
// # CLMuMuVtx  = SeleCuts[0]  #
// # MinMupT    = SeleCuts[1]  #
// # B0LsBS     = SeleCuts[2]  #
// # B0VtxCL    = SeleCuts[3]  #
// # B0cosAlpha = SeleCuts[4]  #
// # B0pT       = SeleCuts[5]  #
// # B0Eta      = SeleCuts[6]  #
// # HadDCASBS  = SeleCuts[7]  #
// # HadpT      = SeleCuts[8]  #
// # KstMass    = SeleCuts[9]  #
// # KKMass     = SeleCuts[10] #
// # NTrkHits   = SeleCuts[11] #
// # NPixLayers = SeleCuts[12] #
// # NormChi2   = SeleCuts[13] #
// # dxyVtx     = SeleCuts[14] #
// # dzVtx      = SeleCuts[15] #
// #############################
{
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // #######################
  // # Read selection cuts #
  // #######################
  std::cout << "\n@@@ Selection cuts from file @@@" << std::endl;
  ParVector.clear();
  ParameterFile->ReadFromFile(ParFileBlockN("selecuts"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++)
    {
      SeleCuts.push_back(atof(ParVector[i].c_str()));
      std::cout << "Selection cut #" << i << " from config file: " << SeleCuts[i] << std::endl;
    }


  ParVector.clear();
  delete ParameterFile;
}

bool Utils::SetSeleCut (std::string cutName, double val)
{
  if      (cutName == "CLMuMuVtx")  SeleCuts[0]  = val;
  else if (cutName == "MinMupT")    SeleCuts[1]  = val;
  else if (cutName == "B0LsBS")     SeleCuts[2]  = val;
  else if (cutName == "B0VtxCL")    SeleCuts[3]  = val;
  else if (cutName == "B0cosAlpha") SeleCuts[4]  = val;
  else if (cutName == "B0pT")       SeleCuts[5]  = val;
  else if (cutName == "B0Eta")      SeleCuts[6]  = val;
  else if (cutName == "HadDCASBS")  SeleCuts[7]  = val;
  else if (cutName == "HadpT")      SeleCuts[8]  = val;
  else if (cutName == "KstMass")    SeleCuts[9]  = val;
  else if (cutName == "KKMass")     SeleCuts[10] = val;
  else if (cutName == "NTrkHits")   SeleCuts[11] = val;
  else if (cutName == "NPixLayers") SeleCuts[12] = val;
  else if (cutName == "NormChi2")   SeleCuts[13] = val;
  else if (cutName == "dxyVtx")     SeleCuts[14] = val;
  else if (cutName == "dzVtx")      SeleCuts[15] = val;
  else return false;

  return true;
}

double Utils::GetSeleCut (std::string cutName)
{
  if      (cutName == "CLMuMuVtx")  return SeleCuts[0];
  else if (cutName == "MinMupT")    return SeleCuts[1];
  else if (cutName == "B0LsBS")     return SeleCuts[2];
  else if (cutName == "B0VtxCL")    return SeleCuts[3];
  else if (cutName == "B0cosAlpha") return SeleCuts[4];
  else if (cutName == "B0pT")       return SeleCuts[5];
  else if (cutName == "B0Eta")      return SeleCuts[6];
  else if (cutName == "HadDCASBS")  return SeleCuts[7];
  else if (cutName == "HadpT")      return SeleCuts[8];
  else if (cutName == "KstMass")    return SeleCuts[9];
  else if (cutName == "KKMass")     return SeleCuts[10];
  else if (cutName == "NTrkHits")   return SeleCuts[11];
  else if (cutName == "NPixLayers") return SeleCuts[12];
  else if (cutName == "NormChi2")   return SeleCuts[13];
  else if (cutName == "dxyVtx")     return SeleCuts[14];
  else if (cutName == "dzVtx")      return SeleCuts[15];
  else
    {
      std::cout << "[Utils::GetSeleCut]\tSelection cut not valid : " << cutName << std::endl;
      exit (1);
    }
}

void Utils::ReadPreselectionCuts (std::string fileName)
// ################################
// # MuMuVtxCL      = PreCuts[0]  #
// # MuMuLsBS       = PreCuts[1]  #
// # DCAMuMu        = PreCuts[2]  #
// # DCAMuBS        = PreCuts[3]  #
// # cosAlphaMuMuBS = PreCuts[4]  #
// # MupT           = PreCuts[5]  #
// # MuEta          = PreCuts[6]  #
// # MuMupT         = PreCuts[7]  #
// # MinMuMuMass    = PreCuts[8]  #
// # MaxMuMuMass    = PreCuts[9]  #
// # MinB0Mass      = PreCuts[10] #
// # MaxB0Mass      = PreCuts[11] #
// # B0VtxCL        = PreCuts[12] #
// # KstMass        = PreCuts[13] #
// # HadDCASBS      = PreCuts[14] #
// # HadpT          = PreCuts[15] #
// ################################
{
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ###########################
  // # Read pre-selection cuts #
  // ###########################
  std::cout << "\n@@@ Pre-selection cuts from file @@@" << std::endl;
  ParVector.clear();
  ParameterFile->ReadFromFile(ParFileBlockN("precuts"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++)
    {
      PreCuts.push_back(atof(ParVector[i].c_str()));
      std::cout << "Pre-selection cut #" << i << " from config file: " << PreCuts[i] << std::endl;
    }


  ParVector.clear();
  delete ParameterFile;
}

bool Utils::SetPreCut (std::string cutName, double val)
{
  if      (cutName == "MuMuVtxCL")      PreCuts[0]  = val;
  else if (cutName == "MuMuLsBS")       PreCuts[1]  = val;
  else if (cutName == "DCAMuMu")        PreCuts[2]  = val;
  else if (cutName == "DCAMuBS")        PreCuts[3]  = val;
  else if (cutName == "cosAlphaMuMuBS") PreCuts[4]  = val;
  else if (cutName == "MupT")           PreCuts[5]  = val;
  else if (cutName == "MuEta")          PreCuts[6]  = val;
  else if (cutName == "MuMupT")         PreCuts[7]  = val;
  else if (cutName == "MinMuMuMass")    PreCuts[8]  = val;
  else if (cutName == "MaxMuMuMass")    PreCuts[9]  = val;
  else if (cutName == "MinB0Mass")      PreCuts[10] = val;
  else if (cutName == "MaxB0Mass")      PreCuts[11] = val;
  else if (cutName == "B0VtxCL")        PreCuts[12] = val;
  else if (cutName == "KstMass")        PreCuts[13] = val;
  else if (cutName == "HadDCASBS")      PreCuts[14] = val;
  else if (cutName == "HadpT")          PreCuts[15] = val;
  else return false;

  return true;
}

double Utils::GetPreCut (std::string cutName)
{
  if      (cutName == "MuMuVtxCL")      return PreCuts[0];
  else if (cutName == "MuMuLsBS")       return PreCuts[1];
  else if (cutName == "DCAMuMu")        return PreCuts[2];
  else if (cutName == "DCAMuBS")        return PreCuts[3];
  else if (cutName == "cosAlphaMuMuBS") return PreCuts[4];
  else if (cutName == "MupT")           return PreCuts[5];
  else if (cutName == "MuEta")          return PreCuts[6];
  else if (cutName == "MuMupT")         return PreCuts[7];
  else if (cutName == "MinMuMuMass")    return PreCuts[8];
  else if (cutName == "MaxMuMuMass")    return PreCuts[9];
  else if (cutName == "MinB0Mass")      return PreCuts[10];
  else if (cutName == "MaxB0Mass")      return PreCuts[11];
  else if (cutName == "B0VtxCL")        return PreCuts[12];
  else if (cutName == "KstMass")        return PreCuts[13];
  else if (cutName == "HadDCASBS")      return PreCuts[14];
  else if (cutName == "HadpT")          return PreCuts[15];
  else
    {
      std::cout << "[Utils::GetPreCut]\tPre-selection cut not valid : " << cutName << std::endl;
      exit (1);
    }
}

void Utils::ReadGenericParam (std::string fileName)
// ########################################
// # B0MassIntervalLeft  = GenericPars[0] #
// # B0MassIntervalRight = GenericPars[1] #
// # NSigmaB0            = GenericPars[2] #
// # NSigmaB0S           = GenericPars[3] #
// # NSigmaB0B           = GenericPars[4] #
// # NSigmaPsiSmall      = GenericPars[5] #
// # NSigmaPsiBig        = GenericPars[6] #
// # SIGMAS1             = GenericPars[7] #
// # SIGMAS2             = GenericPars[8] #
// # FRACMASSS           = GenericPars[9] #
// ########################################
{
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ###########################
  // # Read generic parameters #
  // ###########################
  std::cout << "\n@@@ Generic parameters from file @@@" << std::endl;
  ParVector.clear();
  ParameterFile->ReadFromFile(ParFileBlockN("genericpar"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++)
    {
      GenericPars.push_back(atof(ParVector[i].c_str()));
      std::cout << "Generic parameter #" << i << " from config file: " << GenericPars[i] << std::endl;
    }


  ParVector.clear();
  delete ParameterFile;
}

bool Utils::SetGenericParam (std::string parName, double val)
{
  if      (parName == "B0MassIntervalLeft")  GenericPars[0]  = val;
  else if (parName == "B0MassIntervalRight") GenericPars[1]  = val;
  else if (parName == "NSigmaB0")            GenericPars[2]  = val;
  else if (parName == "NSigmaB0S")           GenericPars[3]  = val;
  else if (parName == "NSigmaB0B")           GenericPars[4]  = val;
  else if (parName == "NSigmaPsiSmall")      GenericPars[5]  = val;
  else if (parName == "NSigmaPsiBig")        GenericPars[6]  = val;
  else if (parName == "SIGMAS1")             GenericPars[7]  = val;
  else if (parName == "SIGMAS2")             GenericPars[8]  = val;
  else if (parName == "FRACMASSS")           GenericPars[9]  = val;
  else return false;

  return true;
}

double Utils::GetGenericParam (std::string parName)
{
  if      (parName == "B0MassIntervalLeft")  return GenericPars[0];
  else if (parName == "B0MassIntervalRight") return GenericPars[1];
  else if (parName == "NSigmaB0")            return GenericPars[2];
  else if (parName == "NSigmaB0S")           return GenericPars[3];
  else if (parName == "NSigmaB0B")           return GenericPars[4];
  else if (parName == "NSigmaPsiSmall")      return GenericPars[5];
  else if (parName == "NSigmaPsiBig")        return GenericPars[6];
  else if (parName == "SIGMAS1")             return GenericPars[7];
  else if (parName == "SIGMAS2")             return GenericPars[8];
  else if (parName == "FRACMASSS")           return GenericPars[9];
  else
    {
      std::cout << "[Utils::GetGenericParam]\tPre-selection cut not valid : " << parName << std::endl;
      exit (1);
    }
}

double Utils::GetB0Width ()
{
  if (GetGenericParam("FRACMASSS") != 0.0) return sqrt(GetGenericParam("FRACMASSS") * pow(GetGenericParam("SIGMAS1"),2.) + (1.-GetGenericParam("FRACMASSS")) * pow(GetGenericParam("SIGMAS2"),2.));
  else                                     return GetGenericParam("SIGMAS1");
}

double* Utils::MakeBinning (std::vector<double>* STLvec)
{
  double* vec = new double[STLvec->size()];
  for (unsigned int i = 0; i < STLvec->size(); i++)
    vec[i] = STLvec->operator[](i);
  return vec;
}

void Utils::ResetEffValue (effValue* myEffVal, double value)
{
  myEffVal->Num1 = value;
  myEffVal->Num2 = value;
  myEffVal->Den1 = value;
  myEffVal->Den2 = value;
  
  // Poissonian errors
  myEffVal->Err2PoisNum1 = value;
  myEffVal->Err2PoisNum2 = value;
  myEffVal->Err2PoisDen1 = value;
  myEffVal->Err2PoisDen2 = value;
  
  // Weight errors
  myEffVal->Err2WeigNum1 = value;
  myEffVal->Err2WeigNum2 = value;
  myEffVal->Err2WeigDen1 = value;
  myEffVal->Err2WeigDen2 = value;
}

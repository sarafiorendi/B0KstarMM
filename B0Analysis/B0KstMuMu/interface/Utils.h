#ifndef UTILS_H
#define UTILS_H

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TMatrixTSym.h>
#include <TGraphAsymmErrors.h>

#include <string>
#include <vector>

#include "B0KstMuMuTreeContent.h"


class Utils
{


 public:
  
  Utils();
  ~Utils() {};
  

  // #################################
  // # Data structure for efficiency #
  // #################################
  struct _effStruct
  {
    double *Num1, *Num2;
    double *Den1, *Den2;
    
    // Poissonian errors
    double *Err2PoisNum1, *Err2PoisNum2;
    double *Err2PoisDen1, *Err2PoisDen2;
    
    // Weight errors
    double *Err2WeigNum1, *Err2WeigNum2;
    double *Err2WeigDen1, *Err2WeigDen2;

    // #########################
    // # Correspondence :      #
    // # GenFilter     <--> N1 #
    // # SingleCand    <--> N2 #
    // # GenNoFilter   <--> D1 #
    // # AllCandFilter <--> D2 #
    // #########################
  };
  typedef struct _effStruct effStruct;

  struct _effValue
  {
    double Num1, Num2;
    double Den1, Den2;
    
    // Poissonian errors
    double Err2PoisNum1, Err2PoisNum2;
    double Err2PoisDen1, Err2PoisDen2;
    
    // Weight errors
    double Err2WeigNum1, Err2WeigNum2;
    double Err2WeigDen1, Err2WeigDen2;
  };
  typedef struct _effValue effValue;


  double computeInvMass (double Px1,
                         double Py1,
                         double Pz1,
                         double mass1,
                         double Px2,
                         double Py2,
                         double Pz2,
                         double mass2);
  
  double computeEta (double Px,
                     double Py,
                     double Pz);
  
  double computePhi (double Px,
		     double Py,
		     double Pz);
  
  double computeEtaPhiDistance (double Px1,
				double Py1,
				double Pz1,
                                double Px2,
                                double Py2,
				double Pz2);
  
  void computeLS (double Vx,
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
		  double* deltaDErr);

  void computeCosAlpha (double Vx,
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
			double* cosAlphaErr);
  
  void ReadBins     (std::string fileName, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins);
  void Readq2Bins   (std::string fileName, std::vector<double>* q2Bins);
  void ReadHLTpaths (std::string fileName, std::vector<std::string>* TrigTable);

  void GenEfficiency      (effStruct* myEff, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins);
  void InitEfficiency     (effValue myEffiVal, effStruct myEff, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins);
  void SaveEfficiency     (std::string fileName, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, effStruct myEff);
  void ReadEfficiency     (std::string fileName, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, effStruct* myEff);
  void GetEffq2Bin        (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int q2Indx, unsigned int cosThetaKIndx, unsigned int cosThetaMuIndx, unsigned int phiIndx, effStruct myEff, double* Eff, double* EffErr);
  TH2D* Get2DEffHitoq2Bin (std::string histoName, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int q2Indx, effStruct myEff);
  TH3D* Get3DEffHitoq2Bin (std::string histoName, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int q2Indx, effStruct myEff);
  void DeleteEfficiency   (effStruct myEff);

  int SearchBin  (double val2Search, std::vector<double>* bins);
  int GetJPsiBin (std::vector<double>* q2Bins);
  int GetPsiPBin (std::vector<double>* q2Bins);

  bool ValIsInPsi (std::vector<double>* q2Bins, double q2Val);
  bool ValIsBetweenJPsiAndPsiP (std::vector<double>* q2Bins, double q2Val);
  bool IsGoodq2CosThetaKCosThetaLBin (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins,
				      double mumuMass2, double cosThetaK, double cosThetaMu, bool PrintMsg);

  void IntegrateEffAll                (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, effStruct myEff, double* Eff, double* EffErr);
  void IntegrateEffButPsi             (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, effStruct myEff, double* Eff, double* EffErr);
  void IntegrateEffInJPsi             (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, effStruct myEff, double* Eff, double* EffErr);
  void IntegrateEffInPsiP             (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, effStruct myEff, double* Eff, double* EffErr);
  void IntegrateEffButq2              (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int q2Indx, effStruct myEff, double* Eff, double* EffErr);
  void IntegrateEffButCosThetaK       (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int cosThetaKBinIndx, effStruct myEff, double* Eff, double* EffErr);
  void IntegrateEffButCosThetaL       (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int cosThetaMuBinIndx, effStruct myEff, double* Eff, double* EffErr);
  void IntegrateEffButPhi             (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int phiIndx, effStruct myEff, double* Eff, double* EffErr);
  void IntegrateEffPhiCosThetaL       (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int q2Indx, unsigned int cosThetaKBinIndx, effStruct myEff, double* Eff, double* EffErr);
  void IntegrateEffPhiCosThetaK       (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int q2Indx, unsigned int cosThetaLBinIndx, effStruct myEff, double* Eff, double* EffErr);
  void IntegrateEffCosThetaKCosThetaL (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int q2Indx, unsigned int phiIndx, effStruct myEff, double* Eff, double* EffErr);
  bool IntegrateEffPhi                (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, double mumuMass2, double cosThetaK, double cosThetaMu, effStruct myEff, double* Eff, double* EffErr, bool PrintMsg);

  unsigned int HLTpathForEvFraction (double evtFrac);
  unsigned int IsInTriggerTable     (B0KstMuMuTreeContent* NTupleIn, double* HLTCutVar1, double* HLTCutVar2, int index = 0, double evtFrac = -1.0);
  unsigned int GetNHLTCat           ();
  double GetHLTentries              (unsigned int HLTindx);
  double ReadLumi                   (std::string fileName);

  void ReadNLLval  (std::string fileName, std::vector<std::vector<double>*>* vecParam);
  double GetNLLval (std::vector<std::vector<double>*>* NLLvals, std::string varName, unsigned int q2BinIndx);

  void ReadTriggerPathsANDCutsANDEntries (std::string fileName);
  void ReadFitStartingValues (std::string fileName, std::vector<std::vector<std::string>*>* vecParam, std::vector<std::vector<unsigned int>*>* configParam, const unsigned int dataBlockN);
  void ReadFitSystematics (std::string fileName, std::vector<std::vector<double>*>* vecParam);

  void SaveAnalyticalEff (std::string fileName, TF2* effFunc, double q2Val, std::vector<double>* q2Bins);
  void SaveAnalyticalEff (std::string fileName, TF3* effFunc, double q2Val, std::vector<double>* q2Bins);
  void SaveAnalyticalEffFullCovariance (std::string fileName, TMatrixTSym<double>* covMatrix, double q2Val, std::vector<double>* q2Bins);

  std::string TellMeEffFuncThetaKThetaLPhi ();
  std::string TellMeEffFuncThetaKThetaL    ();
  std::string TellMeEffFuncThetaK          ();
  std::string TellMeEffFuncThetaL          ();
  std::string TellMeEffFuncPhi             ();

  void ReadAnalyticalEff (std::string fileNameEffParams,
			  std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins,
			  std::vector<TF2*>* effFuncs, std::string effID, const unsigned int dataBlockN);
  void ReadAnalyticalEff (std::string fileNameEffParams,
			  std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins,
			  std::vector<TF3*>* effFuncs, std::string effID, const unsigned int dataBlockN);
  void ReadAnalyticalEffFullCovariance (std::string fileNameEffParams, std::vector<TMatrixTSym<double>*>* covMatrices, const unsigned int dataBlockN);

  double EffMinValue1D (double minX, double maxX, TF1* effFunc);
  double EffMinValue2D (std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, TF2* effFunc);
  double EffMinValue3D (std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, TF3* effFunc);

  void MakeGraphVar (std::string parFileName, TGraphAsymmErrors** graph, std::string varName, bool allBins, double offset = 0.0);

  void InitEffFuncThetaL (TF1* fitFun, unsigned int q2BinIndx);
  void InitEffFuncThetaK (TF1* fitFun, unsigned int q2BinIndx);
  void InitEffFuncPhi    (TF1* fitFun, unsigned int q2BinIndx);

  void AddConstraint1D     (TH1D** histo, std::string constrType, double err, double Yval, double Yerr, unsigned int ID);
  void AddConstraintThetaL (TH1D** histo, unsigned int q2BinIndx, unsigned int cosThetaLBinIndx, double constrXerr, double constrYval, double constrYerr, unsigned int ID);
  void AddConstraint2D     (TH2D** histo, double err, double Zval, double Zerr, unsigned int ID, std::string toBeConstr, std::vector<std::string>* toBeAdded = NULL);
  void AddConstraintThetaK (TH2D** histo, std::vector<double>* cosThetaKBins, unsigned int q2BinIndx, double constrXYerr, double constrZval, double constrZerr, unsigned int ID);
  void AddConstraint3D     (TH3D** histo, double err, double Tval, double Terr, unsigned int ID, std::vector<int> toBeAdded[]);
  void AddConstraintThetaKThetaLPhi (TH3D** histo, unsigned int q2BinIndx, double constrXYZerr, double constrTval, double constrTerr, unsigned int ID);

  bool IsThereOddDegree (TF2* effFunc);

  bool IsThisData (std::string fileName);

  void SaveFitValues (std::string fileName, std::vector<std::string>* vecParStr, int indx, std::string str = "");

  unsigned int ParFileBlockN (std::string blockName);

  unsigned int GetFitParamIndx    (std::string varName);
  unsigned int GetConfigParamIndx (std::string varName);

  bool ChooseBestCand (B0KstMuMuTreeContent* NTuple, unsigned int DoTrigCheck, double evFraction, int* BestCandIndx, bool* B0notB0bar, int* TrigCat, unsigned int* countCands);
  bool FlavorTagger   (B0KstMuMuTreeContent* NTuple, unsigned int i, bool* B0notB0bar);

  void   ReadSelectionCuts (std::string fileName);
  bool   SetSeleCut        (std::string cutName, double val);
  double GetSeleCut        (std::string cutName);

  void   ReadPreselectionCut (std::string fileName);
  bool   SetPreCut           (std::string cutName, double val);
  double GetPreCut           (std::string cutName);

  void   ReadGenericParam (std::string fileName);
  bool   SetGenericParam  (std::string parName, double val);
  double GetGenericParam  (std::string parName);

  double GetB0Width ();

  double* MakeBinning (std::vector<double>* STLvec);

  void ResetEffValue (effValue* myEffVal, double value);


  double muonMass;
  double pionMass;
  double kaonMass;
  double kstMass;
  double B0Mass;
  double JPsiMass;
  double PsiPMass;

  double JPsiBF;
  double JPsiKpiBF;
  double KstMuMuBF;
  double KstKpiMuMuBF;
  double PsiPBF;
  double PsiPKpiBF;

  double muonMassErr;
  double pionMassErr;
  double kaonMassErr;
  double B0MassErr;
  double kstSigma;

  double PI;

  unsigned int NcoeffThetaL;
  unsigned int NcoeffThetaK;
  unsigned int NcoeffPhi;

 
 private:

  TF1* KstMassShape;

  std::vector<std::string> HLTpath;
  std::vector<double> VecHLTCutVar1;
  std::vector<double> VecHLTCutVar2;
  std::vector<double> VecHLTentries;
  std::vector<std::string> TrigTable;

  std::vector<double> PreCuts;
  std::vector<double> SeleCuts;
  std::vector<double> GenericPars;

  double ProbThreshold;
  double scrambleFraction;

  unsigned int nFitParam;
  unsigned int nConfigParam;
  unsigned int nFitObserv;

};

#endif

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
#include <TH3D.h>
#include <TF1.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TLeaf.h>
#include <TObject.h>
#include <TRandom3.h>
#include <TGraphAsymmErrors.h>
#include <TVectorD.h>
#include <TLine.h>
#include <TLatex.h>
#include <TCutG.h>
#include <TKey.h>
#include <TMath.h>
#include <TFitResult.h>

#include <RooRealVar.h>
#include <RooPlot.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>

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
using namespace RooFit;


// ####################
// # Global constants #
// ####################
#define DIRSMCOMP "Data2012B0KstMuMuResults/PredictionSM/"


// #######################
// # Function Definition #
// #######################
void SetStyle               ();
void PlotHistoEff           (string fileName, unsigned int smothDegree, string effDimension, bool RIGHTflavorTAG, double cosThetaKRange_lo = -1.0, double cosThetaKRange_hi = 1.0, double cosThetaLRange_lo = -1.0, double cosThetaLRange_hi = 1.0, double phiRange_lo = -3.15, double phiRange_hi = 3.15);
void TruthMatching          (string fileName, bool truthMatch);
void dBFfromGEN             (string fileName);
void ComputePileUp          (string fileName);
void PlotVtxWithPileUpW     (string fileNameMC, string fileNameData, unsigned int TrigCat, bool withWeights);
void PlotCutScans           (string fileName, string type);
void ReduceTree             (string fileNameIn, string fileNameOut, bool isSingleNotMultyCand);
void SampleMCforPileup      (string fileNameIn, string fileNameOut);
void SampleMCforHadpT       (string fileNameIn, string fileNameOut);
void DivideNTuple           (string fileNameIn, string fileNameOut, unsigned int n);
void SampleNTuple           (string fileNameIn, string fileNameOut, double fraction);
void ComputeMCfilterEff     (string fileName);
void ZeroCrossing           (string fileName, const double minq2 = 1.2, const double maxq2 = 7.0, const unsigned int nBins = 100);
void PlotKKMass             (string fileNameData, string fileNameMC);
// ####################
// # Plot the results #
// ####################
void DrawString             (double Lumi);
TCutG* DrawExclusion        (double Xlow, double Xhigh, double Ylow, double Yhigh, string cutName, unsigned int fillStyle, unsigned int color);
void printData              (TVectorD V1, TVectorD V2, TVectorD V3, TVectorD V4, TVectorD V5, TVectorD V6);
void offsetData             (TVectorD* V1, TVectorD* V2, TVectorD* V3, double offset);
TGraphAsymmErrors* readData (TString fileName, int dataType, int nBins, int color, int markerType, bool doFill, int fillStyle, bool noHbar, double offset);
void showData               (int dataType, double offset, bool noHbar);
void combineMeasurements    (string whichVar, int whichBin);


// ###########################
// # Function Implementation #
// ###########################


// ##########################
// # Set histo layout style #
// ##########################
void SetStyle ()
{
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetTextFont(42);

  gStyle->SetOptFit(1112);
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);

  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadTopMargin(0.11);
  gStyle->SetPadBottomMargin(0.12);

  gStyle->SetTitleFont(42,"x");
  gStyle->SetTitleFont(42,"y");
  gStyle->SetTitleOffset(1.05,"x");
  gStyle->SetTitleOffset(0.95,"y");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"y");

  gStyle->SetLabelFont(42,"x");
  gStyle->SetLabelFont(42,"y");
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");

  TGaxis::SetMaxDigits(3);
  gStyle->SetStatY(0.9);
}


// #########################################################################
// # Sub-program to plot the binned efficicency as "seen" in the final pdf #
// #########################################################################
void PlotHistoEff (string fileName, unsigned int smothDegree, string effDimension, bool RIGHTflavorTAG, double cosThetaKRange_lo, double cosThetaKRange_hi, double cosThetaLRange_lo, double cosThetaLRange_hi, double phiRange_lo, double phiRange_hi)
// #######################
// # effDimension = "2D" #
// # effDimension = "3D" #
// #######################
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  SetStyle();
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);


  ifstream inputFile;
  double xx, xw, yy, yw, zz, zw, cont, err, tmp;
  vector<double> Xbins;
  vector<double> Ybins;
  vector<double> Zbins;
  double cosThetaBound = 1.0;
  double phiBound      = TMath::Pi();


  cout << "[Macros::PlotHistoEff]\tReading binned efficiency file : " << fileName.c_str() << endl;
  inputFile.open(fileName.c_str(), ifstream::in);
  if (inputFile.good() == false)
    {
      cout << "[Macros::PlotHistoEff]\tError opening file : " << fileName.c_str() << endl;
      exit (EXIT_FAILURE);
    }

  // ##################
  // # Reading Z bins #
  // ##################
  if (effDimension == "3D")
    {
      inputFile >> xx >> xw >> yy >> yw >> zz >> zw >> cont >> err;
      tmp = yy;
      while (yy == tmp)
	{
	  if ((zz >= phiRange_lo) && (zz < phiRange_hi)) Zbins.push_back(zz);
	  inputFile >> xx >> xw >> yy >> yw >> zz >> zw >> cont >> err;
	}
      if (phiBound < phiRange_hi)          Zbins.push_back(phiBound);
      else if (Zbins.back() < phiRange_hi) Zbins.push_back(phiRange_hi);
      inputFile.clear();
      inputFile.seekg(0, ios::beg);
    }

  // ##################
  // # Reading Y bins #
  // ##################
  if (effDimension == "2D") inputFile >> xx >> xw >> yy >> yw >> cont >> err;
  else                      inputFile >> xx >> xw >> yy >> yw >> zz >> zw >> cont >> err;
  tmp = xx;
  while (xx == tmp)
    {
      if ((yy >= cosThetaLRange_lo) && (yy < cosThetaLRange_hi)) Ybins.push_back(yy);
      if (effDimension == "2D") inputFile >> xx >> xw >> yy >> yw >> cont >> err;
      else                      inputFile >> xx >> xw >> yy >> yw >> zz >> zw >> cont >> err;
    }
  if (cosThetaBound < cosThetaLRange_hi)     Ybins.push_back(cosThetaBound);
  else if (Ybins.back() < cosThetaLRange_hi) Ybins.push_back(cosThetaLRange_hi);
  inputFile.clear();
  inputFile.seekg(0, ios::beg);

  // ##################
  // # Reading X bins #
  // ##################
  if (effDimension == "2D") inputFile >> xx >> xw >> yy >> yw >> cont >> err;
  else                      inputFile >> xx >> xw >> yy >> yw >> zz >> zw >> cont >> err;
  tmp = xx;
  while (inputFile.eof() == false)
    {
      if ((xx >= cosThetaKRange_lo) && (xx < cosThetaKRange_hi)) Xbins.push_back(xx);
      while ((xx == tmp) && (inputFile.eof() == false))
	{
	  if (effDimension == "2D") inputFile >> xx >> xw >> yy >> yw >> cont >> err;
	  else                      inputFile >> xx >> xw >> yy >> yw >> zz >> zw >> cont >> err;
	}
      tmp = xx;
    }
  if (cosThetaBound < cosThetaKRange_hi)     Xbins.push_back(cosThetaBound);
  else if (Xbins.back() < cosThetaKRange_hi) Xbins.push_back(cosThetaKRange_hi);
  inputFile.clear();
  inputFile.seekg(0, ios::beg);


  cout << "[Macros::PlotHistoEff]\tNew X-axis binning" << std::endl;
  double* Xbins_ = new double[Xbins.size()];
  for (unsigned int i = 0; i < Xbins.size(); i++)
    {
      Xbins_[i] = Xbins[i];
      cout << "bin #" << i << " --> " << Xbins[i] << endl;
    }

  cout << "[Macros::PlotHistoEff]\tNew Y-axis binning" << std::endl;
  double* Ybins_ = new double[Ybins.size()];
  for (unsigned int i = 0; i < Ybins.size(); i++)
    {
      Ybins_[i] = Ybins[i];
      cout << "bin #" << i << " --> " << Ybins[i] << endl;
    }

  double* Zbins_ = NULL;
  if (effDimension == "3D")
    {
      cout << "[Macros::PlotHistoEff]\tNew Z-axis binning" << std::endl;
      Zbins_ = new double[Zbins.size()];
      for (unsigned int i = 0; i < Zbins.size(); i++)
	{
	  Zbins_[i] = Zbins[i];
	  cout << "bin #" << i << " --> " << Zbins[i] << endl;
	}
    }


  // ####################
  // # Plotting section #
  // ####################

  TPad* tmpPad;      
  TCanvas* c0 = new TCanvas("c0","c0",1200,800);
  c0->Divide(4,2);
  c0->cd(1);

  TH1D* projX;
  TH1D* projY;
  TH1D* projZ;

  TH2D* Histo2D;
  TH2D* Histo2D_clone = NULL;

  TH3D* Histo3D;
  TH3D* Histo3D_clone = NULL;

  RooRealVar thetaK("thetaK","cos(#theta#lower[-0.4]{_{#font[122]{K}}})",Xbins[0],Xbins[Xbins.size()-1],"");
  RooRealVar thetaL("thetaL","cos(#theta#lower[-0.4]{_{#font[12]{l}}})",Ybins[0],Ybins[Ybins.size()-1],"");
  RooRealVar phi;

  RooPlot* xframe = thetaK.frame(Name("thetaK"));
  RooPlot* yframe = thetaL.frame(Name("thetaL"));
  RooPlot* zframe = NULL;

  RooDataHist* histoEff;
  RooHistPdf* histoEffPDF;

  if (effDimension == "2D")
    {
      Histo2D = new TH2D("Histo2D", "Histo2D", Xbins.size()-1, Xbins_, Ybins.size()-1, Ybins_);
      Histo2D->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
      Histo2D->GetXaxis()->SetTitleOffset(1.8);
      Histo2D->SetYTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
      Histo2D->GetYaxis()->SetTitleOffset(1.8);
      Histo2D->SetZTitle("Efficiency");
      Histo2D_clone = (TH2D*)Histo2D->Clone();

      // ###################################
      // # Read binned efficiency and plot #
      // ###################################
      unsigned int j = 1;
      while (j <= Xbins.size()-1)
	{
	  unsigned int k = 1;
	  while (k <= Ybins.size()-1)
	    {
	      inputFile >> xx >> xw >> yy >> yw >> cont >> err;

	      if ((xx >= Xbins[0]) && (xx <= Xbins[Xbins.size()-1]) &&
		  (yy >= Ybins[0]) && (yy <= Ybins[Ybins.size()-1]))
		{
		  Histo2D->SetBinContent(j,k,cont);
		  Histo2D->SetBinError(j,k,err);

		  if (RIGHTflavorTAG == true)
		    {
 		      cont = Histo2D->GetBinContent(j,k) * Histo2D->GetXaxis()->GetBinWidth(j) * Histo2D->GetYaxis()->GetBinWidth(k);
		      Histo2D_clone->SetBinContent(j,k,cont);
		    }
		  else
		    {
		      cont = Histo2D->GetBinContent(j,k) * Histo2D->GetXaxis()->GetBinWidth(Histo2D->GetNbinsX()-j+1) * Histo2D->GetYaxis()->GetBinWidth(Histo2D->GetNbinsY()-k+1);
		      Histo2D_clone->SetBinContent(Histo2D->GetNbinsX()-j+1,Histo2D->GetNbinsY()-k+1,cont);
		    }
		  k++;
		}
	    }
	  if (k != 1 ) j++;
	}

      Histo2D->Draw("surf1 fb");

      // #####################
      // # Project histogram #
      // #####################
      tmpPad = static_cast<TPad*>(c0->cd(6));
      tmpPad->SetGrid();
      projX = static_cast<TH1D*>(Histo2D_clone->ProjectionX());
      projX->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
      projX->GetXaxis()->SetTitleOffset(1.0);
      projX->SetYTitle("Projected efficiency");
      projX->SetLineWidth(3);
      projX->SetLineColor(kRed);
      projX->SetMinimum(0);
      projX->Draw("hist");

      c0->cd(7);
      tmpPad = static_cast<TPad*>(c0->cd(7));
      tmpPad->SetGrid();
      projY = static_cast<TH1D*>(Histo2D_clone->ProjectionY());
      projY->SetXTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
      projY->GetXaxis()->SetTitleOffset(1.0);
      projY->SetYTitle("Projected efficiency");
      projY->SetLineWidth(3);
      projY->SetLineColor(kRed);
      projY->SetMinimum(0);
      projY->Draw("hist");

      // ########################
      // # Project data and pdf #
      // ########################
      histoEff    = new RooDataHist("histoEff","histoEff",RooArgSet(thetaK,thetaL),Import(*Histo2D,true));
      histoEffPDF = new RooHistPdf("histoEffPDF","histoEffPDF",RooArgSet(thetaK,thetaL),*histoEff,smothDegree);

      tmpPad = static_cast<TPad*>(c0->cd(2));
      tmpPad->SetGrid();
      histoEff->plotOn(xframe,DataError(RooAbsData::None));
      histoEffPDF->plotOn(xframe,Project(thetaL));
      xframe->Draw();

      tmpPad = static_cast<TPad*>(c0->cd(3));
      tmpPad->SetGrid();
      histoEff->plotOn(yframe,DataError(RooAbsData::None));
      histoEffPDF->plotOn(yframe,Project(thetaK));
      yframe->Draw();
    }
  else
    {
      Histo3D = new TH3D("Histo3D", "Histo3D", Xbins.size()-1, Xbins_, Ybins.size()-1, Ybins_, Zbins.size()-1, Zbins_);
      Histo3D->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
      Histo3D->GetXaxis()->SetTitleOffset(1.8);
      Histo3D->SetYTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
      Histo3D->GetYaxis()->SetTitleOffset(1.8);
      Histo3D->SetZTitle("#phi");
      Histo3D->GetZaxis()->SetTitleOffset(1.8);
      Histo3D_clone = (TH3D*)Histo3D->Clone();

      // ###################################
      // # Read binned efficiency and plot #
      // ###################################
      unsigned int j = 1;
      while (j <= Xbins.size()-1)
	{
	  unsigned int k = 1;
	  while (k <= Ybins.size()-1)
	    {
	      unsigned int l = 1;
	      while (l <= Zbins.size()-1)
		{
		  inputFile >> xx >> xw >> yy >> yw >> zz >> zw >> cont >> err;

		  if ((xx >= Xbins[0]) && (xx <= Xbins[Xbins.size()-1]) &&
		      (yy >= Ybins[0]) && (yy <= Ybins[Ybins.size()-1]) &&
		      (zz >= Zbins[0]) && (zz <= Zbins[Zbins.size()-1]))
		    {
		      Histo3D->SetBinContent(j,k,l,cont);
		      Histo3D->SetBinError(j,k,l,err);

		      if (RIGHTflavorTAG == true)
			{
			  cont = Histo3D->GetBinContent(j,k,l) * Histo3D->GetXaxis()->GetBinWidth(j) * Histo3D->GetYaxis()->GetBinWidth(k) * Histo3D->GetZaxis()->GetBinWidth(l);
			  Histo3D_clone->SetBinContent(j,k,l,cont);
			}
		      else
			{
			  cont = Histo3D->GetBinContent(j,k,l) * Histo3D->GetXaxis()->GetBinWidth(Histo3D->GetNbinsX()-j+1) * Histo3D->GetYaxis()->GetBinWidth(Histo3D->GetNbinsY()-k+1) * Histo3D->GetZaxis()->GetBinWidth(Histo3D->GetNbinsZ()-l+1);
			  Histo3D_clone->SetBinContent(Histo3D->GetNbinsX()-j+1,Histo3D->GetNbinsY()-k+1,Histo3D->GetNbinsZ()-l+1,cont);
			}
		      l++;
		    }
		}
	      if (l != 1 ) k++;
	    }
	  if (k != 1 ) j++;
	}

      Histo3D->Draw();

      // #####################
      // # Project histogram #
      // #####################
      tmpPad = static_cast<TPad*>(c0->cd(6));
      tmpPad->SetGrid();
      projX = static_cast<TH1D*>(Histo3D_clone->Project3D("yz"));
      projX->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
      projX->GetXaxis()->SetTitleOffset(1.0);
      projX->SetYTitle("Projected efficiency");
      projX->SetLineWidth(3);
      projX->SetLineColor(kRed);
      projX->SetMinimum(0);
      projX->Draw("hist");

      tmpPad = static_cast<TPad*>(c0->cd(7));
      tmpPad->SetGrid();
      projY = static_cast<TH1D*>(Histo3D_clone->Project3D("xz"));
      projY->SetXTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
      projY->GetXaxis()->SetTitleOffset(1.0);
      projY->SetYTitle("Projected efficiency");
      projY->SetLineWidth(3);
      projY->SetLineColor(kRed);
      projY->SetMinimum(0);
      projY->Draw("hist");

      tmpPad = static_cast<TPad*>(c0->cd(8));
      tmpPad->SetGrid();
      projZ = static_cast<TH1D*>(Histo3D_clone->Project3D("xy"));
      projZ->SetXTitle("#phi");
      projZ->GetXaxis()->SetTitleOffset(1.0);
      projZ->SetYTitle("Projected efficiency");
      projZ->SetLineWidth(3);
      projZ->SetLineColor(kRed);
      projZ->SetMinimum(0);
      projZ->Draw("hist");

      // ########################
      // # Project data and pdf #
      // ########################
      phi = RooRealVar("phi","cos(#theta#lower[-0.4]{_{#font[12]{l}}})",Zbins[0],Zbins[Zbins.size()-1],"rad");
      zframe = phi.frame(Name("phi"));

      histoEff    = new RooDataHist("histoEff","histoEff",RooArgSet(thetaK,thetaL,phi),Import(*Histo3D,true));
      histoEffPDF = new RooHistPdf("histoEffPDF","histoEffPDF",RooArgSet(thetaK,thetaL,phi),*histoEff,smothDegree);

      tmpPad = static_cast<TPad*>(c0->cd(2));
      tmpPad->SetGrid();
      histoEff->plotOn(xframe,DataError(RooAbsData::None));
      histoEffPDF->plotOn(xframe,Project(RooArgSet(thetaL,phi)));
      xframe->Draw();
      
      tmpPad = static_cast<TPad*>(c0->cd(3));
      tmpPad->SetGrid();
      histoEff->plotOn(yframe,DataError(RooAbsData::None));
      histoEffPDF->plotOn(yframe,Project(RooArgSet(thetaK,phi)));
      yframe->Draw();

      tmpPad = static_cast<TPad*>(c0->cd(4));
      tmpPad->SetGrid();
      histoEff->plotOn(zframe,DataError(RooAbsData::None));
      histoEffPDF->plotOn(zframe,Project(RooArgSet(thetaK,thetaL)));
      zframe->Draw();
    }

  c0->Modified();
  c0->Update();
}


// ###################################################################
// # Sub-program script to check the goodness of truthMatching on MC #
// ###################################################################
void TruthMatching (string fileName, bool truthMatch)
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  SetStyle();
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1001110);


  // #################
  // # Read the tree #
  // #################
  TFile* _file0 = TFile::Open(fileName.c_str(),"READ");
  TTree* B0KstMuMuNTuple = (TTree*)_file0->Get("B0KstMuMu/B0KstMuMuNTuple");


  double minX    = 4.55;
  double maxX    = 6.10;
  double intMinX = 5.00;
  double intMaxX = 5.56;
  unsigned int nBins;
  if (truthMatch == true) nBins = 600;
  else                    nBins = 100;

  TH1D* hb = new TH1D("hb","hb",nBins,minX,maxX);
  hb->SetXTitle("m(K #pi #mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}} #mu#kern[-0.9]{#lower[0.6]{^{#font[122]{\55}}}}) (GeV)");
  hb->SetYTitle("Entries");

  TH1D* hbar = new TH1D("hbar","hbar",nBins,minX,maxX);
  hbar->SetXTitle("m(K #pi #mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}} #mu#kern[-0.9]{#lower[0.6]{^{#font[122]{\55}}}}) (GeV)");
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
      f0->SetParameter(1,0.06);
      f0->SetParameter(2,0.03);
      f0->SetParameter(3,14000.0);
      f0->SetParameter(4,70000.0);
    }
  else
    {
      B0KstMuMuNTuple->Draw("bMass>>hb","genSignal == 1 && truthMatchSignal == 0 && genSignHasFSR == 0");
      B0KstMuMuNTuple->Draw("bBarMass>>hbar","genSignal == 2 && truthMatchSignal == 0 && genSignHasFSR == 0");
      hb->Add(hbar);

      f0 = new TF1("f0","[3]*TMath::Gaus(x,[0],[1]) + [4]*TMath::Gaus(x,[0],[2]) + ([5]+[6]*(x-[7])*(x-[7]))",3.8,6.8);

      f0->SetParName(0,"G-Mean1");
      f0->SetParName(1,"G-Sigma1");
      f0->SetParName(2,"G-Sigma2");

      f0->SetParName(3,"G-Ampli1");
      f0->SetParName(4,"G-Ampli2");

      f0->SetParName(5,"P-Offset");
      f0->SetParName(6,"P-Ampli");
      f0->SetParName(7,"P-Shift");

      // #############################################
      // # Measured from fit with truthMatch == true #
      // #############################################
      f0->FixParameter(0,5.27962);
      f0->FixParameter(1,0.0571392);
      f0->FixParameter(2,0.0277426);

      f0->SetParameter(3,55.0);
      f0->SetParameter(4,35.0);

      f0->SetParameter(5, 20.0);
      f0->SetParameter(6,-20.0);
      f0->SetParameter(7,5.28);
    }


  hb->Fit("f0","VMR0");
  hb->Draw();
  hb->GetFunction("f0")->SetLineColor(kBlue);
  hb->GetFunction("f0")->Draw("same");

  cout << "\nIntegrals low bound: " << intMinX << "\thigh bound: " << intMaxX << endl;
  cout << "Integral of the full fit function: " << f0->Integral(intMinX,intMaxX)/((maxX - minX) / static_cast<double>(nBins)) << "\tEntries: " << hb->GetEntries() << endl;


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


  cout << "Integral of the signal: " << f1->Integral(intMinX,intMaxX)/((maxX - minX) / static_cast<double>(nBins)) << endl;
}


// #########################################################
// # Sub-program to make the dBF/dq2 vs q2 plot for GEN-MC #
// #########################################################
void dBFfromGEN (string fileName)
{
  TFile* _file0 = TFile::Open(fileName.c_str(),"READ");
  TCanvas* c0   = (TCanvas*)_file0->Get("cHistoMeas");
  TPad* p0      = (TPad*)c0->GetPrimitive("cHistoMeas_1");
  TH1D* h0      = (TH1D*)p0->GetPrimitive("histoMeas0");


  // ##########
  // # Signal #
  // ##########
  double GENsigMCev = 5e9;
  double BFsig      = 1.06e-6; // As from MC configuration file (no error associated because it's known)
  double scaleF     = 1e7;

  // #########
  // # J/psi #
  // #########
  double GENpsiEv = 4.996e9;


  // ###########
  // # psi(2S) #
  // ###########
  // double GENpsiEv = 5e9;


  // ###################################
  // # Actual fraction of GEN-MC used  #
  // # for time constraint purposes or #
  // # files too big to be processed   #
  // ###################################
  double MCfrac = 0.5;


  // ###################################################################################################################################################
  // # The differential branching-fraction is: Y_Signal / Y_Ctr[35021242 for J/psi (35033720 for psi(2S))] * Y_GENctr / Y_GENsignal * BF[Signal] / dq2 #
  // ###################################################################################################################################################
  h0->SetBinContent(1,0833328.0 / (35021242.0 * MCfrac) * GENpsiEv / GENsigMCev * BFsig / h0->GetBinWidth(1) * scaleF);
  h0->SetBinContent(2,1661590.0 / (35021242.0 * MCfrac) * GENpsiEv / GENsigMCev * BFsig / h0->GetBinWidth(2) * scaleF);
  h0->SetBinContent(3,1303920.0 / (35021242.0 * MCfrac) * GENpsiEv / GENsigMCev * BFsig / h0->GetBinWidth(3) * scaleF);
  h0->SetBinContent(4,2415750.0 / (35021242.0 * MCfrac) * GENpsiEv / GENsigMCev * BFsig / h0->GetBinWidth(4) * scaleF);
  h0->SetBinContent(5,0.0       / (35021242.0 * MCfrac) * GENpsiEv / GENsigMCev * BFsig / h0->GetBinWidth(5) * scaleF);
  h0->SetBinContent(6,2881710.0 / (35021242.0 * MCfrac) * GENpsiEv / GENsigMCev * BFsig / h0->GetBinWidth(6) * scaleF);
  h0->SetBinContent(7,0.0       / (35021242.0 * MCfrac) * GENpsiEv / GENsigMCev * BFsig / h0->GetBinWidth(7) * scaleF);
  h0->SetBinContent(8,1734920.0 / (35021242.0 * MCfrac) * GENpsiEv / GENsigMCev * BFsig / h0->GetBinWidth(8) * scaleF);
  h0->SetBinContent(9,1836050.0 / (35021242.0 * MCfrac) * GENpsiEv / GENsigMCev * BFsig / h0->GetBinWidth(9) * scaleF);

  h0->SetBinError(1,h0->GetBinContent(1) * sqrt(1./0833328.0 + 1./35021242.0 + 1./GENpsiEv + 1./GENsigMCev));
  h0->SetBinError(2,h0->GetBinContent(2) * sqrt(1./1661590.0 + 1./35021242.0 + 1./GENpsiEv + 1./GENsigMCev));
  h0->SetBinError(3,h0->GetBinContent(3) * sqrt(1./1303920.0 + 1./35021242.0 + 1./GENpsiEv + 1./GENsigMCev));
  h0->SetBinError(4,h0->GetBinContent(4) * sqrt(1./2415750.0 + 1./35021242.0 + 1./GENpsiEv + 1./GENsigMCev));
  h0->SetBinError(5,h0->GetBinContent(5) * sqrt(0.0));
  h0->SetBinError(6,h0->GetBinContent(6) * sqrt(1./2881710.0 + 1./35021242.0 + 1./GENpsiEv + 1./GENsigMCev));
  h0->SetBinError(7,h0->GetBinContent(7) * sqrt(0.0));
  h0->SetBinError(8,h0->GetBinContent(8) * sqrt(1./1734920.0 + 1./35021242.0 + 1./GENpsiEv + 1./GENsigMCev));
  h0->SetBinError(9,h0->GetBinContent(9) * sqrt(1./1836050.0 + 1./35021242.0 + 1./GENpsiEv + 1./GENsigMCev));


  for (int i = 0; i < h0->GetNbinsX(); i++) cout << "--> bin #" << i+1 << "\tentry: " << h0->GetBinContent(i+1) << " +/- " << h0->GetBinError(i+1) << "\twidth: " << h0->GetBinWidth(i+1) << endl;

  TCanvas* c1 = new TCanvas("c1","c1",10,10,700,500);
  c1->cd();
  h0->Draw();
  c1->Modified();
  c1->Update();
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
  c0->Modified();
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
  SetStyle();
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);


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

  c0->Modified();
  c0->Update();
}


// #################################
// # Sub-program to plot cut scans #
// #################################
void PlotCutScans (string fileName, string type)
{
  stringstream myString;

  TFile* _file = TFile::Open(fileName.c_str(),"READ");

  myString.clear(); myString.str("");
  myString << "c0_" << type;
  TCanvas* c0 = (TCanvas*)_file->Get(myString.str().c_str());

  myString.clear(); myString.str("");
  myString << "c0_" << type << "_1";
  TPad* p0 = (TPad*)c0->GetPrimitive(myString.str().c_str());
  TH1D* h0 = (TH1D*)p0->GetPrimitive("histoR1");
  h0->SetLineWidth(2);

  myString.clear(); myString.str("");
  myString << "c0_" << type << "_4";
  TPad* p1 = (TPad*)c0->GetPrimitive(myString.str().c_str());
  TH1D* h1 = (TH1D*)p1->GetPrimitive("histoR4");
  h1->SetYTitle("S / sqrt(S+B)");
  h1->SetLineWidth(2);


  // ##########################
  // # Set histo layout style #
  // ##########################
  SetStyle();
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);


  TCanvas* c1 = new TCanvas("c1","c1",10,10,900,500);
  c1->Divide(2,1);
  c1->cd(1);
  h0->Draw("CE5");
  TH1D* tmp0 = (TH1D*)h0->Clone();
  tmp0->ResetAttFill();
  tmp0->Draw("C same");

  TLine* myLine0 = new TLine(h1->GetXaxis()->GetBinLowEdge(h1->GetMaximumBin()),h0->GetMaximum()*1.07,h1->GetXaxis()->GetBinLowEdge(h1->GetMaximumBin()),h0->GetMinimum()*0.86);
  myLine0->SetLineColor(kRed);
  myLine0->SetLineWidth(3);
  myLine0->SetLineStyle(2);
  myLine0->Draw("same");

  c1->cd(2);
  h1->Draw("CE5");
  TH1D* tmp1 = (TH1D*)h1->Clone();
  tmp1->ResetAttFill();
  tmp1->Draw("C same");

  TLine* myLine1 = new TLine(h1->GetXaxis()->GetBinLowEdge(h1->GetMaximumBin()),h1->GetMaximum()*1.07,h1->GetXaxis()->GetBinLowEdge(h1->GetMaximumBin()),h1->GetMinimum()*0.86);
  myLine1->SetLineColor(kRed);
  myLine1->SetLineWidth(3);
  myLine1->SetLineStyle(2);
  myLine1->Draw("same");

  
  c1->Modified();
  c1->Update();
}


// ########################################################################
// # Sub-program to reduce the number of branches of the candidate ntuple #
// ########################################################################
void ReduceTree (string fileNameIn, string fileNameOut, bool isSingleNotMultyCand)
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
  if (isSingleNotMultyCand == false)
    {
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
    }
  else
    {
      theTreeIn->SetBranchStatus("B0MassArb",1);
      theTreeIn->SetBranchStatus("kstMassArb",1);

      theTreeIn->SetBranchStatus("mumuMass",1);
      theTreeIn->SetBranchStatus("mumuMassE",1);

      theTreeIn->SetBranchStatus("truthMatchSignal",1);
      theTreeIn->SetBranchStatus("rightFlavorTag",1);

      theTreeIn->SetBranchStatus("CosThetaKArb",1);
      theTreeIn->SetBranchStatus("CosThetaMuArb",1);
      theTreeIn->SetBranchStatus("PhiKstMuMuPlaneArb",1);
    }

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


// ################################################################################################################
// # Sub-program to sample the MC in order to obtain a subset that has the same hadron pT distribution as in data #
// ################################################################################################################
void SampleMCforHadpT (string fileNameIn, string fileNameOut)
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
  vector<double>* vec3 = new vector<double>;
  vector<double>* vec4 = new vector<double>;

  
  TH1D* hadppTW = new TH1D("hadppTW","hadppTW",100,0,20);
  theTreeIn->Draw("sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy)>>hadppTW","evWeight*(((genSignal == 3) || (genSignal == 4)) && (truthMatchSignal == 1))","goff");
  hadppTW->Scale(1. / hadppTW->Integral());

  TH1D* hadppTNoW = new TH1D("hadppTNoW","hadppTNoW",100,0,20);
  theTreeIn->Draw("sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy)>>hadppTNoW","(((genSignal == 3) || (genSignal == 4)) && (truthMatchSignal == 1))","goff");
  hadppTNoW->Scale(1. / hadppTNoW->Integral());

  TH1D* hadppTRatio = (TH1D*)hadppTW->Clone("hadppTRatio");
  hadppTRatio->Divide(hadppTNoW);
  hadppTRatio->Scale(1. / hadppTRatio->GetBinContent(hadppTRatio->GetMaximumBin()));

  
  TH1D* hadmpTW = new TH1D("hadmpTW","hadmpTW",100,0,20);
  theTreeIn->Draw("sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy)>>hadmpTW","evWeight*(((genSignal == 3) || (genSignal == 4)) && (truthMatchSignal == 1))","goff");
  hadmpTW->Scale(1. / hadmpTW->Integral());

  TH1D* hadmpTNoW = new TH1D("hadmpTNoW","hadmpTNoW",100,0,20);
  theTreeIn->Draw("sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy)>>hadmpTNoW","(((genSignal == 3) || (genSignal == 4)) && (truthMatchSignal == 1))","goff");
  hadmpTNoW->Scale(1. / hadmpTNoW->Integral());

  TH1D* hadmpTRatio = (TH1D*)hadmpTW->Clone("hadmpTRatio");
  hadmpTRatio->Divide(hadmpTNoW);
  hadmpTRatio->Scale(1. / hadmpTRatio->GetBinContent(hadmpTRatio->GetMaximumBin()));


  for (int entry = 0; entry < nEntries; entry++)
    {
      theTreeIn->SetBranchAddress("kstTrkpPx", &vec1);
      theTreeIn->SetBranchAddress("kstTrkpPy", &vec2);
      
      theTreeIn->SetBranchAddress("kstTrkmPx", &vec3);
      theTreeIn->SetBranchAddress("kstTrkmPy", &vec4);
      
      theTreeIn->GetEntry(entry);

      for (unsigned int it = 0; it < (*vec1).size(); it++)
	if (myRandom->Uniform() < 0.5)
	  {
	    if  ((myRandom->Uniform() < hadppTRatio->GetBinContent(hadppTRatio->FindBin(sqrt((*vec1)[it]*(*vec1)[it] + (*vec2)[it]*(*vec2)[it]))))) theTreeOut->Fill();
	  }
	else if ((myRandom->Uniform() < hadmpTRatio->GetBinContent(hadmpTRatio->FindBin(sqrt((*vec3)[it]*(*vec3)[it] + (*vec4)[it]*(*vec4)[it]))))) theTreeOut->Fill();


      vec1->clear();
      vec2->clear();
      
      vec3->clear();      
      vec4->clear();
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
  cout << "\n@@@ Dividing the input tree in ==> FOUR <== HLT sub-categories @@@" << endl;
  TFile* fileTmp = TFile::Open("fileTmp.root","RECREATE");
  fileTmp->mkdir("B0KstMuMu");
  fileTmp->cd("B0KstMuMu");
  theTreeInCat.push_back(theTreeIn->CopyTree("TrigCat == 1 && truthMatchSignal == 1"));
  theTreeInCat.push_back(theTreeIn->CopyTree("TrigCat == 2 && truthMatchSignal == 1"));
  theTreeInCat.push_back(theTreeIn->CopyTree("TrigCat == 3 && truthMatchSignal == 1"));
  theTreeInCat.push_back(theTreeIn->CopyTree("TrigCat == 4 && truthMatchSignal == 1"));


  fileNameOut.replace(fileNameOut.find(".root"),5,"");
  for (unsigned int i = 0; i < n; i++)
    {
      myString.clear(); myString.str("");
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
  cout << "Dividing the input tree in the ==> FOUR <== HLT sub-categories" << endl;
  theTreeInCat.push_back(theTreeIn->CopyTree("TrigCat == 1 && truthMatchSignal == 1"));
  theTreeInCat.push_back(theTreeIn->CopyTree("TrigCat == 2 && truthMatchSignal == 1"));
  theTreeInCat.push_back(theTreeIn->CopyTree("TrigCat == 3 && truthMatchSignal == 1"));
  theTreeInCat.push_back(theTreeIn->CopyTree("TrigCat == 4 && truthMatchSignal == 1"));


  // #####################################################
  // # Select from each HLT category a sub-set of events #
  // #####################################################
  for (unsigned int i = 0; i < theTreeInCat.size(); i++)
    {
      theTreeOutCat.push_back(theTreeInCat[i]->CloneTree(static_cast<int>(rint(static_cast<double>(theTreeInCat[i]->GetEntries()) * fraction))));
      cout << "Selecting from HLT category #" << i+1 << " a sub-set of events : " << theTreeOutCat.back()->GetEntries() << endl;
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


  unsigned int val1 = 0;
  unsigned int val2 = 0;
  unsigned long long evTried  = 0;
  unsigned long long evPassed = 0;

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


// #############################################################################################
// # Compute zero-crossing point for AFB applying the definition of forward-backward asymmetry #
// #############################################################################################
void ZeroCrossing (string fileName, const double minq2, const double maxq2, const unsigned int nBins)
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  SetStyle();
  gStyle->SetPalette(1);


  TFile *_file0 = TFile::Open(fileName.c_str(),"READ");
  TTree* B0KstMuMuNTuple = (TTree*)_file0->Get("B0KstMuMu/B0KstMuMuNTuple");

  TH1D* hAFB = new TH1D("hAFB","hAFB",nBins,0,8.5);
  hAFB->GetXaxis()->SetTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
  hAFB->GetYaxis()->SetTitle("A_{FB}");
  hAFB->SetMarkerStyle(20);

  TH1D* h1p = new TH1D("h1p","h1p",nBins,0,8.5);
  h1p->GetXaxis()->SetTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
  h1p->GetYaxis()->SetTitle("Entries (#)");
  h1p->SetMarkerStyle(20);
  h1p->SetMarkerColor(kBlack);
  h1p->GetXaxis()->SetRangeUser(minq2,maxq2);

  TH1D* h1m = new TH1D("h1m","h1m",nBins,0,8.5);
  h1m->GetXaxis()->SetTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
  h1m->GetYaxis()->SetTitle("Entries (#)");
  h1m->SetMarkerStyle(21);
  h1m->SetMarkerColor(kRed);
  h1m->GetXaxis()->SetRangeUser(minq2,maxq2);

  TH2D* h2p = new TH2D("h2p","h2p",nBins,0,8.5,nBins,0,1.0);
  h2p->GetXaxis()->SetTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
  h2p->GetYaxis()->SetTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
  TH2D* h2m = new TH2D("h2m","h2m",nBins,0,8.5,nBins,0,-1.0);
  h2m->GetXaxis()->SetTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
  h2m->GetYaxis()->SetTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");

  TF1* ZeroCrox = new TF1("ZeroCrox","[0]*x + [1]",minq2,maxq2);
  ZeroCrox->SetParName(0,"a");
  ZeroCrox->SetParName(1,"b");
  ZeroCrox->SetLineColor(kRed);
  ZeroCrox->SetLineWidth(2);


  // ####################
  // # Quering the tree #
  // ####################
  B0KstMuMuNTuple->Draw("CosThetaMuArb:mumuMass*mumuMass>>h2p","mumuMass*mumuMass > 1.0 && mumuMass*mumuMass < 8.5 && CosThetaMuArb > 0","goff");
  B0KstMuMuNTuple->Draw("CosThetaMuArb:mumuMass*mumuMass>>h2m","mumuMass*mumuMass > 1.0 && mumuMass*mumuMass < 8.5 && CosThetaMuArb < 0","goff");


  TCanvas* c0 = new TCanvas("c0","c0",10,10,900,500);
  c0->Divide(2,1);

  c0->cd(1);
  h2p->GetXaxis()->SetRangeUser(minq2,maxq2);
  h2p->Draw("gcolz");

  c0->cd(2);
  h2m->GetXaxis()->SetRangeUser(minq2,maxq2);
  h2m->Draw("gcolz");

  c0->Modified();
  c0->Update();


  // ######################
  // # Making projections #
  // ######################
  for (unsigned int i = 0; i < nBins; i++)
    {
      double counterP = 0.0;
      double counterM = 0.0;
      
      for (unsigned int j = 0; j < nBins; j++)
	{
	  counterP += h2p->GetBinContent(i+1,j+1);
	  counterM += h2m->GetBinContent(i+1,j+1);
	}
      
      h1p->SetBinContent(i+1,counterP);
      h1m->SetBinContent(i+1,counterM);
    }


  TCanvas* c1 = new TCanvas("c1","c1",10,10,900,500);
  c1->cd();
  h1p->Draw("e1p");
  h1m->Draw("same e1p");
  c0->Modified();
  c0->Update();


  // ########################################
  // # Computing forward-backward asymmetry #
  // ########################################
  TCanvas* c2 = new TCanvas("c2","c2",10,10,900,500);
  c2->cd();


  for (unsigned int i = 0; i < nBins; i++)
    {
      if ((h1p->GetBinContent(i+1) != 0.0) || (h1m->GetBinContent(i+1) != 0.0))
	{
	  hAFB->SetBinContent(i+1,(h1p->GetBinContent(i+1) - h1m->GetBinContent(i+1)) / (h1p->GetBinContent(i+1) + h1m->GetBinContent(i+1)));
	  hAFB->SetBinError(i+1,2.*h1p->GetBinContent(i+1)*h1m->GetBinContent(i+1) /
			    (pow(h1p->GetBinContent(i+1) + h1m->GetBinContent(i+1),2.)) *
			    sqrt(pow(h1p->GetBinError(i+1) / h1p->GetBinContent(i+1),2.) + pow(h1m->GetBinError(i+1) / h1m->GetBinContent(i+1),2.)));
	}
      else
	{
	  hAFB->SetBinContent(i+1,0.0);
	  hAFB->SetBinError(i+1,0.0);
	}
    }

  
  // ###########
  // # Fitting #
  // ###########
  hAFB->GetXaxis()->SetRangeUser(minq2,maxq2);
  hAFB->Draw("e1p");
  TFitResultPtr fitResults = hAFB->Fit("ZeroCrox","S R0");
  TMatrixTSym<double> covMatrix(fitResults->GetCovarianceMatrix());
  ZeroCrox->Draw("same");


  // ################################
  // # Compute zero crossing porint #
  // ################################
  double q0  = -ZeroCrox->GetParameter(1) / ZeroCrox->GetParameter(0);
  double q0E = q0 * sqrt(pow(ZeroCrox->GetParError(0) / ZeroCrox->GetParameter(0),2.) + pow(ZeroCrox->GetParError(1) / ZeroCrox->GetParameter(1),2.) - 2./(ZeroCrox->GetParameter(0) * ZeroCrox->GetParameter(1)) * covMatrix(0,1));
  cout << "\n@@@ Zero crossing point: " << q0 << " +/- " << q0E << " @@@" << endl;
  cout << "Fit range: [" << minq2 << "-" << maxq2 << "]" << endl;
  cout << "Number of bins: " << nBins << endl;


  c2->Modified();
  c2->Update();
}


// #######################################################################
// # Overlap m(KK) plots from Data and MC and fit for the phi(1020) peak #
// #######################################################################
void PlotKKMass (string fileNameData, string fileNameMC)
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  SetStyle();
  gStyle->SetPalette(1);


  double cutVal = 1.035;
  double yRange = 0.034;
  cout << "\nCut value: " << cutVal << endl;
  cout << "Y high bound: " << yRange << "\n" << endl;


  TFile* fileData = TFile::Open(fileNameData.c_str(),"READ");
  TCanvas* cData = (TCanvas*)fileData->Get("c0");
  TH1D* hData = (TH1D*)cData->GetPrimitive("hKKSig");
  hData->Sumw2();
  hData->Scale(1. / hData->Integral());
  hData->GetYaxis()->SetRangeUser(0.0,yRange);
  hData->SetMarkerStyle(20);
  hData->SetXTitle("m(K^{+}K^{-}) (GeV)");
  hData->SetYTitle("Entries / (0.005 GeV)");

  TFile* fileMC = TFile::Open(fileNameMC.c_str(),"READ");
  TCanvas* cMC = (TCanvas*)fileMC->Get("c0");
  TH1D* hMC = (TH1D*)cMC->GetPrimitive("hKKSig");
  hMC->Sumw2();
  hMC->Scale(1. / hMC->Integral());
  hMC->GetYaxis()->SetRangeUser(0.0,yRange);
  hMC->SetFillColor(kAzure+6);
  hMC->SetXTitle("m(K^{+}K^{-}) (GeV)");
  hMC->SetYTitle("Entries / (0.005 GeV)");


  TCanvas * c0 = new TCanvas("c0","c0",10,10,900,500);
  c0->cd();
  hMC->Draw("hist");
  hData->Draw("same e1");


  TF1* func = new TF1("func","[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))",1,cutVal);
  func->SetLineColor(kViolet);
  func->SetLineWidth(2);
  func->SetParName(0,"Ampli");
  func->SetParName(1,"#mu");
  func->SetParName(2,"#sigma");

  func->SetParameter(0,0.022);
  func->SetParameter(1,1.02);
  func->SetParameter(2,0.005);

  hData->Fit("func","R");
  TPaveStats* st = (TPaveStats*)hData->FindObject("stats");
  st->Draw("same");

  TLine* myLine = new TLine(cutVal,0.0,cutVal,yRange - 0.0005);
  myLine->SetLineColor(kRed);
  myLine->SetLineWidth(3);
  myLine->SetLineStyle(2);
  myLine->Draw("same");

  c0->Modified();
  c0->Update();
}


// #######################################
// # Code to plot experiment-comparisons #
// #######################################
void DrawString (double Lumi)
{
  stringstream myString;
  double scaleRespect2CMS = 0.75;


  myString.clear(); myString.str("");
  myString << "CMS";
  TLatex* LumiTex1 = new TLatex(0.11,0.9,myString.str().c_str());
  LumiTex1->SetTextFont(61);
  LumiTex1->SetTextSize(0.05);
  LumiTex1->SetTextColor(kBlack);
  LumiTex1->SetNDC(true);
  LumiTex1->DrawLatex(0.11,0.9,myString.str().c_str());


  myString.clear(); myString.str("");
  myString << "#it{Preliminary}";
  TLatex* LumiTex2 = new TLatex(0.19,0.9,myString.str().c_str());
  LumiTex2->SetTextFont(42);
  LumiTex2->SetTextSize(0.05 * scaleRespect2CMS);
  LumiTex2->SetTextColor(kBlack);
  LumiTex2->SetNDC(true);
  LumiTex2->DrawLatex(0.19,0.9,myString.str().c_str());


  myString.clear(); myString.str("");
  myString << Lumi <<  " fb#lower[0.4]{^{#font[122]{\55}1}} (8 TeV)";
  TLatex* LumiTex3 = new TLatex(0.8,0.9,myString.str().c_str());
  LumiTex3->SetTextFont(42);
  LumiTex3->SetTextSize(0.05 * scaleRespect2CMS);
  LumiTex3->SetTextColor(kBlack);
  LumiTex3->SetNDC(true);
  LumiTex3->DrawLatex(0.8,0.9,myString.str().c_str());
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


void printData (TVectorD V1, TVectorD V2, TVectorD V3, TVectorD V4, TVectorD V5, TVectorD V6)
{
  for (int i = 0; i < V1.GetNoElements(); i++)
    cout << "Read values: " << V1[i] << "\t" << V2[i] << "\t" << V3[i] << "\t" << V4[i] << "\t" << V5[i] << "\t" << V6[i] << endl;

}


void offsetData (TVectorD* V1, TVectorD* V2, TVectorD* V3, double offset)
// offset = percentage of half of the the bin width [0,1]
{
  double shift;

  for (int i = 0; i < V1->GetNoElements(); i++)
    {
      shift = (fabs((*V2)[i]) + fabs((*V3)[i])) / 2. * offset;

      (*V1)[i] = (*V1)[i] + shift;
      if ((*V2)[i] != 0.0) (*V2)[i] = (*V2)[i] + shift;
      if ((*V3)[i] != 0.0) (*V3)[i] = (*V3)[i] - shift;
    }
}


TGraphAsymmErrors* readData (TString fileName, int dataType, int nBins, int color, int markerType, bool doFill, int fillStyle, bool noHbar, double offset)
// #################################
// # dataType = 0 --> FL           #
// # dataType = 1 --> AFB          #
// # dataType = 2 --> BF           #
// # noHbar: choose horizontal bar #
// #################################
{
  double YvalueOutsideLimits = 20.0; // Value given to bins with zero error in order not to show them

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
  cout << "\n\nReading file: " << fileName << "\tN. entries: " << nEntrie << "\tN. bins: " << nBins << endl;
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
      
      if ((V5[i] == 0.0) && (V6[i] == 0.0)) V4[i] = YvalueOutsideLimits;
    }


  cout << "\nData before offset:" << endl;
  printData(V1,V2,V3,V4,V5,V6);
  offsetData(&V1,&V2,&V3,offset);
  cout << "\nData after offset:" << endl;
  printData(V1,V2,V3,V4,V5,V6);

  TGraphAsymmErrors* gra = new TGraphAsymmErrors(V1,V4,V2,V3,V5,V6);
  if      (dataType == 0) gra->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}}); F_{L}");
  else if (dataType == 1) gra->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}}); A_{FB}");
  else if (dataType == 2) gra->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}}); dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}8}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");
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
  SetStyle();
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);


  double luminosity = 20.5; // CMS data luminosity
  stringstream myString;

  myString.clear(); myString.str("");
  myString << DIRSMCOMP << "BFvsq2Template.root";
  TFile* _file0 = TFile::Open(myString.str().c_str(),"READ");
  TCanvas* c0   = (TCanvas*)_file0->Get("cHistoMeas");
  TPad* p0      = (TPad*)c0->GetPrimitive("cHistoMeas_1");
  TH1D* h0      = (TH1D*)p0->GetPrimitive("histoMeas0");
  h0->GetXaxis()->SetRangeUser(0.0,19.0);
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
      h0->GetYaxis()->SetRangeUser(0.0,12.);
      h0->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}}); dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}8}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");
    }
  h0->SetLineStyle(2);
  cout << "I've read the templare of the q2 binning" << endl;


  TCanvas* cData  = new TCanvas("cData","cData",10,10,700,500);
  vector<TGraphAsymmErrors*> dVar;

  myString.clear(); myString.str("");
  myString << DIRSMCOMP << "CMS_7e8TeV.data";
  dVar.push_back(readData(myString.str().c_str(),dataType,h0->GetNbinsX(),1,20,false,0,noHbar,0.0*offset));

  // myString.clear(); myString.str("");
  // myString << DIRSMCOMP << "CMS_8TeV.data";
  // dVar.push_back(readData(myString.str().c_str(),dataType,h0->GetNbinsX(),1,20,false,0,noHbar,0.0*offset));

  // myString.clear(); myString.str("");
  // myString << DIRSMCOMP << "CMS_7TeV.data";
  // dVar.push_back(readData(myString.str().c_str(),dataType,h0->GetNbinsX(),1,24,false,0,noHbar,0.3*offset));

  myString.clear(); myString.str("");
  myString << DIRSMCOMP << "LHCb_1fb.data";
  dVar.push_back(readData(myString.str().c_str(),dataType,h0->GetNbinsX()-3,2,21,false,0,noHbar,-0.3*offset));

  // myString.clear(); myString.str("");
  // myString << DIRSMCOMP << "LHCb_3fb.data";
  // dVar.push_back(readData(myString.str().c_str(),dataType,h0->GetNbinsX()-1,2,21,false,0,noHbar,0.0*offset));

  // myString.clear(); myString.str("");
  // myString << DIRSMCOMP << "Atlas.data";
  // if (dataType != 2) dVar.push_back(readData(myString.str().c_str(),dataType,h0->GetNbinsX()-3,6,22,false,0,noHbar,0.8*offset));

  myString.clear(); myString.str("");
  myString << DIRSMCOMP << "BaBar.data";
  dVar.push_back(readData(myString.str().c_str(),dataType,h0->GetNbinsX()-3,4,23,false,0,noHbar,-0.4*offset));

  myString.clear(); myString.str("");
  myString << DIRSMCOMP << "Belle.data";
  dVar.push_back(readData(myString.str().c_str(),dataType,h0->GetNbinsX()-3,8,28,false,0,noHbar,-0.6*offset));

  myString.clear(); myString.str("");
  myString << DIRSMCOMP << "CDF.data";
  dVar.push_back(readData(myString.str().c_str(),dataType,h0->GetNbinsX()-3,kGray+1,29,false,0,noHbar,-0.8*offset));


  cData->cd();
  h0->Draw();
  for (unsigned int i = dVar.size(); i > 0; i--) dVar[i-1]->Draw("same p");
  dVar[0]->SetLineWidth(3);


  unsigned int it = 0;
  TLegend* leg = NULL;
  leg = new TLegend(0.12, 0.6, 0.27, 0.88, "");
  leg->AddEntry(dVar[it++],"CMS (7 TeV + 8 TeV)","lp");
  // leg->AddEntry(dVar[it++],"CMS (8 TeV)","lp");
  // leg->AddEntry(dVar[it++],"CMS (7 TeV)","lp");
  leg->AddEntry(dVar[it++],"LHCb","lp");
  // if (dataType != 2) leg->AddEntry(dVar[it++],"Atlas","lp");
  leg->AddEntry(dVar[it++],"BaBar","lp");
  leg->AddEntry(dVar[it++],"Belle","lp");
  leg->AddEntry(dVar[it++],"CDF","lp");

  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();


  TLine* line = new TLine(-1.0,0.0,19.0,0.0);
  line->SetLineStyle(kDashed);
  line->Draw();

  
  DrawString(luminosity);
  DrawExclusion(08.68,10.09,-1.2,30.,"RejectJPsi1",3001,kGray);
  DrawExclusion(12.86,14.18,-1.2,30.,"RejectPsiP1",3001,kGray);


  cData->Modified();
  cData->Update();
 }


void combineMeasurements (string whichVar, int whichBin)
// ############################
// # whichVar --> FL          #
// # whichVar --> AFB         #
// # whichVar --> BF          #
// # whichBin --> 0,1,2,4,6,7 #
// ############################
{
  // #############
  // # Variables #
  // #############
  unsigned int nMeas = 2;
  double scale       = 0.;
  double combMeas    = 0.;
  double combMeasVar = 0.;
  stringstream myString;

  TMatrixD CovHi(nMeas,nMeas);
  TMatrixD CovLo(nMeas,nMeas);
  TMatrixD Cov(nMeas,nMeas);
  vector<double> meas(nMeas);
  vector<double> weights(nMeas);


  // ##################
  // # Initialisation #
  // ##################
  if (whichVar == "FL")
    {
      if (whichBin == -1)
	{
	  // # FL specialBin #
	  meas[0]  = 0.704;
	  meas[1]  = 0.684;

	  CovHi(0,0) = 0.051*0.051 + (0.003*0.003+0.001*0.001+0.014*0.014+0.004*0.004+0.007*0.007);
	  CovHi(1,1) = 0.102*0.102 + (0.002*0.002+0.008*0.008+0.013*0.013+0.006*0.006);

	  CovLo(0,0) = 0.052*0.052 + (0.003*0.003+0.001*0.001+0.014*0.014+0.004*0.004+0.007*0.007);
	  CovLo(1,1) = 0.102*0.102 + (0.002*0.002+0.008*0.008+0.013*0.013+0.006*0.006);

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * (0.034*0.034+0.003*0.003+0.003*0.003+0.014*0.014); // cov[x,y] = rho * sigmax * sigmay
	}
      else if (whichBin == 0)
	{
	  meas[0]  = 0.639;
	  meas[1]  = 0.600;

	  CovHi(0,0) = 0.087*0.087 + (0.010*0.010+0.011*0.011+0.027*0.027+0.015*0.015+0.006*0.006);
	  CovHi(1,1) = 0.000*0.000 + (0.007*0.007+0.040*0.040+0.023*0.023+0.179*0.179);

	  CovLo(0,0) = 0.094*0.094 + (0.010*0.010+0.011*0.011+0.027*0.027+0.015*0.015+0.006*0.006);
	  CovLo(1,1) = 0.280*0.280 + (0.007*0.007+0.040*0.040+0.023*0.023+0.179*0.179);

	 CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * (0.034*0.034+0.003*0.003+0.003*0.003+0.023*0.023); // cov[x,y] = rho * sigmax * sigmay
	}
      else if (whichBin == 1)
	{
	  meas[0]  = 0.793;
	  meas[1]  = 0.650;

	  CovHi(0,0) = 0.080*0.080 + (0.003*0.003+0.004*0.004+0.003*0.003+0.004*0.004);
	  CovHi(1,1) = 0.170*0.170 + (0.007*0.007+0.003*0.003+0.019*0.019+0.003*0.003);

	  CovLo(0,0) = 0.085*0.085 + (0.003*0.003+0.004*0.004+0.003*0.003+0.004*0.004);
	  CovLo(1,1) = 0.170*0.170 + (0.007*0.007+0.003*0.003+0.019*0.019+0.003*0.003);

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * (0.034*0.034+0.004*0.004+0.003*0.003+0.024*0.024); // cov[x,y] = rho * sigmax * sigmay
	}
      else if (whichBin == 2)
	{
	  // #######################
	  // # Combine bin 2 and 3 #
	  // #######################
	  // meas[0]  = 0.613;
	  // meas[1]  = 0.486;

	  // CovHi(0,0) = 0.101*0.101 + (0.004*0.004+0.009*0.009+0.034*0.034+0.009*0.009);
	  // CovHi(1,1) = 0.061*0.061 + (0.005*0.005+0.017*0.017+0.010*0.010+0.006*0.006+0.027*0.027);

	  // CovLo(0,0) = 0.096*0.096 + (0.004*0.004+0.009*0.009+0.034*0.034+0.009*0.009);
	  // CovLo(1,1) = 0.062*0.062 + (0.005*0.005+0.017*0.017+0.010*0.010+0.006*0.006+0.027*0.027);

	  // CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 0.0;

	  meas[0]  = 0.527;
	  meas[1]  = 0.810;

	  CovHi(0,0) = 0.057*0.057;
	  CovHi(1,1) = 0.131*0.131 + (0.007*0.007+0.040*0.040+0.009*0.009+0.013*0.013+0.023*0.023);

	  CovLo(0,0) = 0.057*0.057;
	  CovLo(1,1) = 0.124*0.124 + (0.007*0.007+0.040*0.040+0.009*0.009+0.013*0.013+0.023*0.023);

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * ((0.034*0.034+0.003*0.003+0.003*0.003+0.014*0.014) + (0.034*0.034+0.003*0.003+0.003*0.003+0.021*0.021)); // cov[x,y] = rho * sigmax * sigmay
	}
      else if (whichBin == 4)
	{
	  meas[0]  = 0.392;
	  meas[1]  = 0.450;

	  CovHi(0,0) = 0.050*0.050 + (0.003*0.003+0.002*0.002+0.004*0.004);
	  CovHi(1,1) = 0.100*0.100 + (0.005*0.005+0.013*0.013+0.005*0.005+0.026*0.026+0.008*0.008);

	  CovLo(0,0) = 0.051*0.051 + (0.003*0.003+0.002*0.002+0.004*0.004);
	  CovLo(1,1) = 0.110*0.110 + (0.005*0.005+0.013*0.013+0.005*0.005+0.026*0.026+0.008*0.008);

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * (0.034*0.034+0.001*0.001+0.003*0.003+0.008*0.008); // cov[x,y] = rho * sigmax * sigmay
	}
      else if (whichBin == 6)
	{
	  meas[0]  = 0.464;
	  meas[1]  = 0.530;

	  CovHi(0,0) = 0.041*0.041 + (0.005*0.005+0.001*0.001+0.001*0.001+0.013*0.013+0.001*0.001);
	  CovHi(1,1) = 0.120*0.120 + (0.006*0.006+0.023*0.023+0.013*0.013+0.004*0.004);

	  CovLo(0,0) = 0.049*0.049 + (0.005*0.005+0.001*0.001+0.001*0.001+0.013*0.013+0.001*0.001);
	  CovLo(1,1) = 0.070*0.070 + (0.006*0.006+0.023*0.023+0.013*0.013+0.004*0.004);

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * (0.034*0.034+0.001*0.001+0.003*0.003+0.006*0.006); // cov[x,y] = rho * sigmax * sigmay
	}
      else if (whichBin == 7)
	{
	  meas[0]  = 0.379;
	  meas[1]  = 0.440;

	  CovHi(0,0) = 0.054*0.054 + (0.005*0.005+0.003*0.003+0.001*0.001+0.008*0.008);
	  CovHi(1,1) = 0.070*0.070 + (0.005*0.005+0.011*0.011+0.011*0.011+0.021*0.021);

	  CovLo(0,0) = 0.055*0.055 + (0.005*0.005+0.003*0.003+0.001*0.001+0.008*0.008);
	  CovLo(1,1) = 0.070*0.070 + (0.005*0.005+0.011*0.011+0.011*0.011+0.021*0.021);

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * (0.034*0.034+0.004*0.004+0.003*0.003+0.006*0.006); // cov[x,y] = rho * sigmax * sigmay
	}
      else { cout << "\nInvalid bin number: " << whichBin << endl; exit(EXIT_FAILURE); }
    }
  else if (whichVar == "AFB")
    {
      if (whichBin == -1)
	{
	  // # AFB specialBin #
	  meas[0]  = -0.150;
	  meas[1]  = -0.068;

	  CovHi(0,0) = 0.096*0.096 + (0.004*0.004+0.012*0.012+0.016*0.016+0.019*0.019+0.005*0.005);
	  CovHi(1,1) = 0.115*0.115 + (0.001*0.001+0.003*0.003+0.013*0.013+0.002*0.002);

	  CovLo(0,0) = 0.090*0.090 + (0.004*0.004+0.012*0.012+0.016*0.016+0.019*0.019+0.005*0.005);
	  CovLo(1,1) = 0.116*0.116 + (0.001*0.001+0.003*0.003+0.013*0.013+0.002*0.002);

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * (0.008*0.008+0.003*0.003+0.001*0.001+0.002*0.002); // cov[x,y] = rho * sigmax * sigmay
	}
      else if (whichBin == 0)
	{
	  meas[0]  = -0.270;
	  meas[1]  = -0.290;

	  CovHi(0,0) = 0.170*0.170 + (0.007*0.007+0.037*0.037+0.020*0.020+0.052*0.052+0.004*0.004);
	  CovHi(1,1) = 0.280*0.280 + (0.004*0.004+0.077*0.077+0.006*0.006+0.161*0.161);

	  CovLo(0,0) = 0.401*0.401 + (0.007*0.007+0.037*0.037+0.020*0.020+0.052*0.052+0.004*0.004);
	  CovLo(1,1) = 0.000*0.000 + (0.004*0.004+0.077*0.077+0.006*0.006+0.161*0.161);

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * (0.008*0.008+0.002*0.002+0.001*0.001+0.005*0.005); // cov[x,y] = rho * sigmax * sigmay
	}
      else if (whichBin == 1)
	{
	  meas[0]  = -0.123;
	  meas[1]  = -0.070;

	  CovHi(0,0) = 0.151*0.151 + (0.007*0.007+0.015*0.015+0.032*0.032+0.032*0.032+0.003*0.003);
	  CovHi(1,1) = 0.200*0.200 + (0.005*0.005+0.014*0.014+0.008*0.008+0.004*0.004);

	  CovLo(0,0) = 0.136*0.136 + (0.007*0.007+0.015*0.015+0.032*0.032+0.032*0.032+0.003*0.003);
	  CovLo(1,1) = 0.200*0.200 + (0.005*0.005+0.014*0.014+0.008*0.008+0.004*0.004);

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * (0.008*0.008+0.003*0.003+0.001*0.001+0.002*0.002); // cov[x,y] = rho * sigmax * sigmay
	}
      else if (whichBin == 2)
	{
	  // #######################
	  // # Combine bin 2 and 3 #
	  // #######################
	  // meas[0]  = 0.028;
	  // meas[1]  = 0.035;

	  // CovHi(0,0) = 0.152*0.152 + (0.008*0.008+0.020*0.020+0.013*0.013+0.005*0.005);
	  // CovHi(1,1) = 0.097*0.097 + (0.018*0.018+0.008*0.008+0.005*0.005);

	  // CovLo(0,0) = 0.153*0.153 + (0.008*0.008+0.020*0.020+0.013*0.013+0.005*0.005);
	  // CovLo(1,1) = 0.096*0.096 + (0.018*0.018+0.008*0.008+0.005*0.005);

	  // CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 0.0;

	  meas[0]  = 0.033;
	  meas[1]  = -0.010;

	  CovHi(0,0) = 0.083*0.083;
	  CovHi(1,1) = 0.110*0.110 + (0.005*0.005+0.017*0.017+0.014*0.014+0.002*0.002+0.022*0.022);

	  CovLo(0,0) = 0.083*0.083;
	  CovLo(1,1) = 0.110*0.110 + (0.005*0.005+0.017*0.017+0.014*0.014+0.002*0.002+0.022*0.022);

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * ((0.018*0.018+0.008*0.008+0.005*0.005) + (0.008*0.008+0.020*0.020+0.013*0.013+0.005*0.005)); // cov[x,y] = rho * sigmax * sigmay
	}
      else if (whichBin == 4)
	{
	  meas[0]  = 0.157;
	  meas[1]  = 0.400;
  
	  CovHi(0,0) = 0.059*0.059 + (0.005*0.005+0.003*0.003+0.006*0.006+0.003*0.003);
	  CovHi(1,1) = 0.080*0.080 + (0.004*0.004+0.043*0.043+0.014*0.014+0.008*0.008+0.005*0.005);

	  CovLo(0,0) = 0.060*0.060 + (0.005*0.005+0.003*0.003+0.006*0.006+0.003*0.003);
	  CovLo(1,1) = 0.080*0.080 + (0.004*0.004+0.043*0.043+0.014*0.014+0.008*0.008+0.005*0.005);

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * (0.008*0.008+0.0003*0.0003+0.001*0.001+0.002*0.002); // cov[x,y] = rho * sigmax * sigmay
	}
      else if (whichBin == 6)
	{
	  meas[0]  = 0.399;
	  meas[1]  = 0.290;
  
	  CovHi(0,0) = 0.029*0.029 + (0.005*0.005+0.002*0.002+0.006*0.006+0.003*0.003);
	  CovHi(1,1) = 0.090*0.090 + (0.003*0.003+0.014*0.014+0.011*0.011+0.038*0.038+0.017*0.017);

	  CovLo(0,0) = 0.064*0.064 + (0.005*0.005+0.002*0.002+0.006*0.006+0.003*0.003);
	  CovLo(1,1) = 0.090*0.090 + (0.003*0.003+0.014*0.014+0.011*0.011+0.038*0.038+0.017*0.017);

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * (0.008*0.008+0.002*0.002+0.001*0.001+0.001*0.001); // cov[x,y] = rho * sigmax * sigmay
	}
      else if (whichBin == 7)
	{
	  meas[0]  = 0.351;
	  meas[1]  = 0.410;

	  CovHi(0,0) = 0.065*0.065 + (0.005*0.005+0.008*0.008+0.003*0.003);
	  CovHi(1,1) = 0.050*0.050 + (0.004*0.004+0.012*0.012+0.010*0.010+0.016*0.016+0.013*0.013);

	  CovLo(0,0) = 0.066*0.066 + (0.005*0.005+0.008*0.008+0.003*0.003);
	  CovLo(1,1) = 0.050*0.050 + (0.004*0.004+0.012*0.012+0.010*0.010+0.016*0.016+0.013*0.013);

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * (0.008*0.008+0.007*0.007+0.001*0.001); // cov[x,y] = rho * sigmax * sigmay
	}
      else { cout << "\nInvalid bin number: " << whichBin << endl; exit(EXIT_FAILURE); }
    }
  else if (whichVar == "BF")
    {
      if (whichBin == -1)
	{
	  // # dBF/dq2 specialBin #
	  meas[0]  = 3.59;
	  meas[1]  = 4.39;

	  CovHi(0,0) = 0.29*0.29 + (0.7*0.7+2.3*2.3+1.0*1.0) * meas[0]*meas[0] / 1e4;
	  CovHi(1,1) = 0.58*0.58 + (0.6*0.6+1.0*1.0+5.0*5.0) * meas[1]*meas[1] / 1e4;

	  CovLo(0,0) = 0.29*0.29 + (0.7*0.7+2.3*2.3+1.0*1.0) * meas[0]*meas[0] / 1e4;
	  CovLo(1,1) = 0.55*0.55 + (0.6*0.6+1.0*1.0+5.0*5.0) * meas[1]*meas[1] / 1e4;

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * ((4.2*4.2+0.3*0.3+3.2*3.2+1.1*1.1+4.6*4.6) * meas[0]*meas[1] / 1e4) ; // cov[x,y] = rho * sigmax * sigmay
	}
      else if (whichBin == 0)
	{
	  meas[0]  = 4.70;
	  meas[1]  = 4.80;

	  CovHi(0,0) = 0.67*0.67 + (2.0*2.0+2.5*2.5+0.4*0.4) * meas[0]*meas[0] / 1e4;
	  CovHi(1,1) = 1.40*1.40 + (1.0*1.0+2.1*2.1+3.3*3.3+5.0*5.0) * meas[1]*meas[1] / 1e4;

	  CovLo(0,0) = 0.67*0.67 + (2.0*2.0+2.5*2.5+0.4*0.4) * meas[0]*meas[0] / 1e4;
	  CovLo(1,1) = 1.20*1.20 + (1.0*1.0+2.1*2.1+3.3*3.3+5.0*5.0) * meas[1]*meas[1] / 1e4;

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * ((4.5*4.5+0.2*0.2+3.2*3.2+1.5*1.5+4.6*4.6) * meas[0]*meas[1] / 1e4) ; // cov[x,y] = rho * sigmax * sigmay
	}
      else if (whichBin == 1)
	{
	  meas[0]  = 3.31;
	  meas[1]  = 3.80;

	  CovHi(0,0) = 0.44*0.44 + (1.1*1.1+1.9*1.9+0.2*0.2) * meas[0]*meas[0] / 1e4;
	  CovHi(1,1) = 0.70*0.70 + (1.0*1.0+2.7*2.7+5.0*5.0) * meas[1]*meas[1] / 1e4;

	  CovLo(0,0) = 0.44*0.44 + (1.1*1.1+1.9*1.9+0.2*0.2) * meas[0]*meas[0] / 1e4;
	  CovLo(1,1) = 0.70*0.70 + (1.0*1.0+2.7*2.7+5.0*5.0) * meas[1]*meas[1] / 1e4;

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * ((4.5*4.5+0.3*0.3+3.2*3.2+2.1*2.1+4.6*4.6) * meas[0]*meas[1] / 1e4) ; // cov[x,y] = rho * sigmax * sigmay
	}
      else if (whichBin == 2)
	{
	  // #######################
	  // # Combine bin 2 and 3 #
	  // #######################
	  // meas[0]  = 3.45;
	  // meas[1]  = 4.70;

	  // CovHi(0,0) = 0.47*0.47 + (0.8*0.8+5.5*5.5+1.2*1.2) * meas[0]*meas[0] / 1e4;
	  // CovHi(1,1) = 0.42*0.42 + (1.0*1.0+3.1*3.1+0.5*0.5+1.9*1.9) * meas[1]*meas[1] / 1e4;

	  // CovLo(0,0) = 0.47*0.47 + (0.8*0.8+5.5*5.5+1.2*1.2) * meas[0]*meas[0] / 1e4;
	  // CovLo(1,1) = 0.42*0.42 + (1.0*1.0+3.1*3.1+0.5*0.5+1.9*1.9) * meas[1]*meas[1] / 1e4;

	  // CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 0.0;

	  meas[0]  = 4.143;
	  meas[1]  = 3.70;

	  CovHi(0,0) = 0.34*0.34;
	  CovHi(1,1) = 0.70*0.70 + (1.0*1.0+1.3*1.3+2.9*2.9+3.2*3.2+5.0*5.0) * meas[1]*meas[1] / 1e4;

	  CovLo(0,0) = 0.34*0.34;
	  CovLo(1,1) = 0.70*0.70 + (1.0*1.0+1.3*1.3+2.9*2.9+3.2*3.2+5.0*5.0) * meas[1]*meas[1] / 1e4;

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * ((4.5*4.5+0.2*0.2+3.2*3.2+0.7*0.7+4.6*4.6) * meas[0]*meas[1] / 1e4 + (4.5*4.5+0.1*0.1+3.2*3.2+1.4*1.4+4.6*4.6) * meas[0]*meas[1] / 1e4); // cov[x,y] = rho * sigmax * sigmay
	}
      else if (whichBin == 4)
	{
	  meas[0]  = 6.21;
	  meas[1]  = 5.40;

	  CovHi(0,0) = 0.45*0.45 + (0.7*0.7+3.5*3.5+4.0*4.0) * meas[0]*meas[0] / 1e4;
	  CovHi(1,1) = 0.90*0.90 + (1.0*1.0+1.0*1.0+15.2*15.2+5.0*5.0) * meas[1]*meas[1] / 1e4;

	  CovLo(0,0) = 0.45*0.45 + (0.7*0.7+3.5*3.5+4.0*4.0) * meas[0]*meas[0] / 1e4;
	  CovLo(1,1) = 0.90*0.90 + (1.0*1.0+1.0*1.0+15.2*15.2+5.0*5.0) * meas[1]*meas[1] / 1e4;

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * ((4.5*4.5+3.7*3.7+3.2*3.2+0.4*0.4+4.6*4.6) * meas[0]*meas[1] / 1e4) ; // cov[x,y] = rho * sigmax * sigmay
	}
      else if (whichBin == 6)
	{
	  meas[0]  = 6.54;
	  meas[1]  = 4.60;

	  CovHi(0,0) = 0.59*0.59 + (0.7*0.7+2.0*2.0+0.1*0.1) * meas[0]*meas[0] / 1e4;
	  CovHi(1,1) = 0.90*0.90 + (1.0*1.0+2.2*2.2+1.0*1.0+7.1*7.1+5.0*5.0) * meas[1]*meas[1] / 1e4;

	  CovLo(0,0) = 0.59*0.59 + (0.7*0.7+2.0*2.0+0.1*0.1) * meas[0]*meas[0] / 1e4;
	  CovLo(1,1) = 0.80*0.80 + (1.0*1.0+2.2*2.2+1.0*1.0+7.1*7.1+5.0*5.0) * meas[1]*meas[1] / 1e4;

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * ((4.5*4.5+3.8*3.8+3.2*3.2+0.2*0.2+4.6*4.6) * meas[0]*meas[1] / 1e4) ; // cov[x,y] = rho * sigmax * sigmay
	}
      else if (whichBin == 7)
	{
	  meas[0]  = 4.19;
	  meas[1]  = 5.20;

	  CovHi(0,0) = 0.32*0.32 + (0.5*0.5+1.0*1.0+0.3*0.3) * meas[0]*meas[0] / 1e4;
	  CovHi(1,1) = 0.60*0.60 + (1.0*1.0+5.6*5.6+5.0*5.0) * meas[1]*meas[1] / 1e4;

	  CovLo(0,0) = 0.32*0.32 + (0.5*0.5+1.0*1.0+0.3*0.3) * meas[0]*meas[0] / 1e4;
	  CovLo(1,1) = 0.60*0.60 + (1.0*1.0+5.6*5.6+5.0*5.0) * meas[1]*meas[1] / 1e4;

	  CovHi(0,1) = CovHi(1,0) = CovLo(0,1) = CovLo(1,0) = 1. * ((4.5*4.5+4.1*4.1+3.2*3.2+0.2*0.2+4.6*4.6) * meas[0]*meas[1] / 1e4) ; // cov[x,y] = rho * sigmax * sigmay
	}
      else { cout << "\nInvalid bin number: " << whichBin << endl; exit(EXIT_FAILURE); }
    }
  else { cout << "\nInvalid variable: " << whichVar << endl; exit(EXIT_FAILURE); }


  if (meas[0] >  meas[1])
    {
      Cov(0,0) = CovLo(0,0);
      Cov(1,1) = CovHi(1,1);
    }
  else
    {
      Cov(0,0) = CovHi(0,0);
      Cov(1,1) = CovLo(1,1);
    }
  Cov(0,1) = Cov(1,0) = CovHi(0,1);
  

  // ###############
  // # Combination #
  // ###############
  TMatrixD invCov = Cov;
  invCov = invCov.Invert();
  
  for (unsigned int i = 0; i < nMeas; i++)
    for (unsigned int j = 0; j < nMeas; j++)
      scale += invCov(i,j);

  for (unsigned int i = 0; i < nMeas; i++)
    {
      weights[i] = 0;
      for (unsigned int j = 0; j < nMeas; j++) weights[i] += invCov(i,j);
      weights[i] = weights[i] / scale;
      
      combMeas += weights[i] * meas[i];
    }
  
  for (unsigned int i = 0; i < nMeas; i++)
    for (unsigned int j = 0; j < nMeas; j++)
      combMeasVar += weights[i] * weights[j] * Cov(i,j);


  // ##########
  // # Output #
  // ##########
  cout << "\n@@@ Printing measurement central values @@@" << endl;
  for (unsigned int i = 0; i < nMeas; i++) cout << "#" << i << "\t" << meas[i] << endl;

  cout << "\n@@@ Printing weights @@@" << endl;
  for (unsigned int i = 0; i < nMeas; i++) cout << "#" << i << "\t" << weights[i] << endl;

  cout << "\n@@@ Printing covariance @@@" << endl;
  for (unsigned int i = 0; i < nMeas; i++)
    for (unsigned int j = 0; j < nMeas; j++)
      cout << "#" << i << "," << j << "\t" << Cov(i,j) << endl;

  cout << "\n@@@ Printing covariance low @@@" << endl;
  for (unsigned int i = 0; i < nMeas; i++)
    for (unsigned int j = 0; j < nMeas; j++)
      cout << "#" << i << "," << j << "\t" << sqrt(CovLo(i,j)) << endl;

  cout << "\n@@@ Printing covariance high @@@" << endl;
  for (unsigned int i = 0; i < nMeas; i++)
    for (unsigned int j = 0; j < nMeas; j++)
      cout << "#" << i << "," << j << "\t" << sqrt(CovHi(i,j)) << endl;

  cout << "\nCombined measurement = " << combMeas << " +/- " << sqrt(combMeasVar) << endl;
}

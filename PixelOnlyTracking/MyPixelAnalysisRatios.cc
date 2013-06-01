// #############################################################
// # Program to compute efficiencies, purities and resolutions #
// # for pixel stand alone tracking and vertexing              #
// #############################################################
// # Author: Mauro Dinardo                                     #
// #############################################################

#include <TROOT.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

#include <math.h>
#include <map>
#include <stdlib.h>
#include <sstream>
#include <string>


using namespace std;


// ####################
// # Histogram TRACKS #
// ####################

// @@@@@@ Pixel tracks @@@@@@
TH1F* H_pix_trk_pt;
TH1F* H_pix_trk_eta;

TH1F* H_pix_trk_doubleCounter_pt;
TH1F* H_pix_trk_doubleCounter_eta;

TH2F* H_pix_trk_normchi2_pt_eta[2];
TH2F* H_pix_trk_eta_phi;
TH2F* H_pix_trk_eta_pt;

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @ Histograms containing ratios @
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
TH1F* H_pix_trk_doubleCounter_pt_rel;
TH1F* H_pix_trk_doubleCounter_eta_rel;

// @@@@@@ Tracker tracks @@@@@@
TH1F* H_trk_trk_pt;
TH1F* H_trk_trk_eta;

TH1F* H_trk_trk_doubleCounter_pt;
TH1F* H_trk_trk_doubleCounter_eta;

TH2F* H_trk_trk_normchi2_pt_eta[2];
TH2F* H_trk_trk_eta_phi;
TH2F* H_trk_trk_eta_pt;

TH2F* H_trk_trk_eta_phi_algo;
TH2F* H_trk_trk_eta_pt_algo;

TH1F* H_trk_toteff_pt_num;
TH1F* H_trk_toteff_pt_den;
TH1F* H_trk_toteff_eta_num;
TH1F* H_trk_toteff_eta_den;
TH1F* H_trk_toteff_phi_num[3];
TH1F* H_trk_toteff_phi_den[3];

TH1F* H_trk_eff_pt_num;
TH1F* H_trk_eff_pt_den;
TH1F* H_trk_eff_eta_num;
TH1F* H_trk_eff_eta_den;
TH1F* H_trk_eff_phi_num[3];
TH1F* H_trk_eff_phi_den[3];

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @ Histograms containing ratios @
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
TH1F* H_trk_trk_doubleCounter_pt_rel;
TH1F* H_trk_trk_doubleCounter_eta_rel;
TH1F* H_trk_toteff_pt;
TH1F* H_trk_toteff_eta;
TH1F* H_trk_toteff_phi[3];
TH1F* H_trk_eff_pt;
TH1F* H_trk_eff_eta;
TH1F* H_trk_eff_phi[3];

// @@@@@@ Common pixel&tracker tracks @@@@@@
typedef map<string, map<string, TH1F*> > myMap;
myMap MapFit_d0;      // Map to store d0 difference [pt, eta]
myMap MapFit_dz;      // Map to store dz difference [pt, eta]
myMap MapFit_pt;      // Map to store pt difference [pt, eta]
myMap MapFit_d0_pull; // Map to store d0 pulls [pt, eta]
myMap MapFit_dz_pull; // Map to store dz pulls [pt, eta]
myMap MapFit_pt_pull; // Map to store pt pulls [pt, eta]
myMap MapFit_d0err;   // Map to store generalTrack d0err [pt, eta]
myMap MapFit_dzerr;   // Map to store generalTrack dzerr [pt, eta]
myMap MapFit_pterr;   // Map to store generalTrack pterr [pt, eta]

TH2F* H_trked_eta_phi_pur;
TH2F* H_trked_eta_phi_eff;
TH2F* H_trked_eta_pt_pur;
TH2F* H_trked_eta_pt_eff;
TH2F* H_trked_normchi2_pt_pur_eta[2];
TH2F* H_trked_normchi2_pt_eff_eta[2];

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @ Histograms containing ratios @
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
TH1F* H_ptres_eta_pt[5];
TH1F* H_ptpull_eta_pt[5];

TH1F* H_dzres_eta_pt[5];
TH1F* H_dzpull_eta_pt[5];

TH1F* H_d0res_eta_pt[5];
TH1F* H_d0pull_eta_pt[5];

TH1F* H_ptres_pt_eta[2];
TH1F* H_ptpull_pt_eta[2];

TH1F* H_dzres_pt_eta[2];
TH1F* H_dzpull_pt_eta[2];

TH1F* H_d0res_pt_eta[2];
TH1F* H_d0pull_pt_eta[2];

TH1F* H_skew_pt_eta[2];

TH1F* H_purity_trked_eta;
TH1F* H_purity_trked_pt_eta[2];
TH1F* H_purity_trked_pt_whole_eta;
TH1F* H_purity_trked_phi[3];

TH1F* H_efficiency_trked_eta;
TH1F* H_efficiency_trked_pt_eta[2];
TH1F* H_efficiency_trked_pt_whole_eta;
TH1F* H_efficiency_trked_phi[3];

TH1F* H_efficiency_trked_eta_algo;
TH1F* H_efficiency_trked_pt_eta_algo[2];
TH1F* H_efficiency_trked_pt_whole_eta_algo;
TH1F* H_efficiency_trked_phi_algo[3];

TH1F* H_tk_pur_eff;

TH2F* H_purity_trked_normchi2_pt_eta[2];
TH2F* H_efficiency_trked_normchi2_pt_eta[2];


// ######################
// # Histogram VERTICES #
// ######################

// @@@@@@ Common pixel&tracker vertices @@@@@@
TH1F* H_vx_counters;

map<string, TH1F*> map_H_vx_x_diff;
map<string, TH1F*> map_H_vx_y_diff;
map<string, TH1F*> map_H_vx_z_diff;
map<string, TH1F*> map_H_vx_x_diff_pull;
map<string, TH1F*> map_H_vx_y_diff_pull;
map<string, TH1F*> map_H_vx_z_diff_pull;
map<string, TH1F*> map_H_vx_x_diff_orig;
map<string, TH1F*> map_H_vx_y_diff_orig;
map<string, TH1F*> map_H_vx_z_diff_orig;
map<string, TH1F*> map_H_vx_x_err2;
map<string, TH1F*> map_H_vx_y_err2;
map<string, TH1F*> map_H_vx_z_err2;

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @ Histograms containing ratios @
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
TH1F* H_vx_x_diff_vs_trk;
TH1F* H_vx_y_diff_vs_trk;
TH1F* H_vx_z_diff_vs_trk;
TH1F* H_vx_x_diff_pull_vs_trk;
TH1F* H_vx_y_diff_pull_vs_trk;
TH1F* H_vx_z_diff_pull_vs_trk;
TH1F* H_vx_x_diff_orig_vs_trk;
TH1F* H_vx_y_diff_orig_vs_trk;
TH1F* H_vx_z_diff_orig_vs_trk;
TH1F* H_vx_x_diff_noerr_vs_trk;
TH1F* H_vx_y_diff_noerr_vs_trk;
TH1F* H_vx_z_diff_noerr_vs_trk;
TH1F* H_vx_pur_eff;


// @@@@@@ Duplicated variables @@@@@@
double EtaThr;       // --> Must have at most 2 decimal digits
double EtaThrBin[3]; // --> Must have at most 2 decimal digits
double zRange;       // Unit: [um]
double d0Range;      // Unit: [cm]
double dzRange;      // Unit: [cm]
double ptRange;      // Unit: [GeV/c]
double d0errRange;   // Unit: [um]^2
double dzerrRange;   // Unit: [um]^2
double pterrRange;   // Unit: [MeV/c]^2
double pullsRange;
double MaxPtEquiBinning;

// @@@@@@ Non duplicated variables @@@@@@
double PtThr[6];               // --> Must have at most 2 decimal digits
unsigned int MinEntriesEffPur; // Histogram minimum number of entries for efficiency and purity

// @@@@@@ Parameters from cfg file @@@@@@
double MaxEtaTkTrk;            // Max eta that a track should have in order to be taken into account
double RangePt;                // pt bins. Unit: [GeV/c] --> Must have at most 2 decimal digits
double PtStep;                 // pt bins. Unit: [GeV/c] --> Must have at most 2 decimal digits
double RangeEta;               // eta bins --> Must have at most 2 decimal digits
double EtaStep;                // eta bins --> Must have at most 2 decimal digits
double MinPtTk;                // Mini Pt that a track must have in order to be taken into account. Unit: [GeV/c]
                               // --> Must have at most 2 decimal digits


void InitVariables ()
{
  MinEntriesEffPur = 12;

  MaxEtaTkTrk      = 2.3;
  RangePt          = 8.0;
  PtStep           = 0.25;
  RangeEta         = 2.3;
  EtaStep          = 0.1;
  MinPtTk          = 0.75;

  cout << "\n@@@@@@ CONFIGURATION PARAMETERS @@@@@@" << endl;

  cout << "@@@@@@ GENERAL PARAMETERS @@@@@@" << endl;
  cout << "MinEntriesEffPur --> " << MinEntriesEffPur << endl;

  cout << "@@@@@@ TRACK PARAMETERS @@@@@@" << endl;
  cout << "MaxEtaTkTrk --> " << MaxEtaTkTrk << endl;
  cout << "RangePt --> " << RangePt << endl;
  cout << "PtStep --> " << PtStep << endl;
  cout << "RangeEta --> " << RangeEta << endl;
  cout << "EtaStep --> " << EtaStep << endl;
  cout << "MinPtTk --> " << MinPtTk << endl;

  cout << "@@@@@@ VERTEX PARAMETERS @@@@@@" << endl;

  MaxPtEquiBinning = 3.;

  EtaThr = 1.5;
  EtaThrBin[0] = 0;
  EtaThrBin[1] = EtaThr;
  EtaThrBin[2] = MaxEtaTkTrk;

  PtThr[0] = MinPtTk;
  for (int i = 1; i < 5; i++) PtThr[i] = PtThr[i-1] + 0.5;
  PtThr[5] = RangePt;

  zRange     = 6000.;  // Unit: (um)

  d0Range    = 0.6;    // Unit: (cm)
  dzRange    = 0.6;    // Unit: (cm)
  ptRange    = 12.0;   // Unit: (GeV/c)
  d0errRange = 10000.; // Unit: (um)^2
  dzerrRange = 10000.; // Unit: (um)^2
  pterrRange = 50000.; // Unit: (MeV/c)^2

  pullsRange = 10.;
}


bool LoadHistos (string datatype, TFile* filename)
{
  if (datatype.compare("Tk") == 0)
    {
      // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      // @ Histograms containing numerators or denominators @
      // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      H_pix_trk_pt = (TH1F*)filename->Get("MyProcess/PixelTracks/Pixeltracks pt distribution");
      H_pix_trk_eta = (TH1F*)filename->Get("MyProcess/PixelTracks/Pixeltracks eta distribution");
      H_pix_trk_doubleCounter_pt = (TH1F*)filename->Get("MyProcess/PixelTracks/Doublecounter pixeltracks pt distribution");
      H_pix_trk_doubleCounter_eta = (TH1F*)filename->Get("MyProcess/PixelTracks/Doublecounter pixeltracks eta distribution");

      H_pix_trk_eta_phi = (TH2F*)filename->Get("MyProcess/PixelTracks/Pixeltracks distribution in eta and phi");
      H_pix_trk_eta_pt = (TH2F*)filename->Get("MyProcess/PixelTracks/Pixeltracks distribution in eta and pt");

      H_trk_trk_pt = (TH1F*)filename->Get("MyProcess/PixelTracks/Trackertracks pt distribution");
      H_trk_trk_eta = (TH1F*)filename->Get("MyProcess/PixelTracks/Trackertracks eta distribution");
      H_trk_trk_doubleCounter_pt = (TH1F*)filename->Get("MyProcess/PixelTracks/Doublecounter trackertracks pt distribution");
      H_trk_trk_doubleCounter_eta = (TH1F*)filename->Get("MyProcess/PixelTracks/Doublecounter trackertracks eta distribution");

      H_trk_trk_eta_phi = (TH2F*)filename->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi");
      H_trk_trk_eta_pt = (TH2F*)filename->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and pt");

      H_trk_trk_eta_phi_algo = (TH2F*)filename->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi algo");
      H_trk_trk_eta_pt_algo = (TH2F*)filename->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and pt algo");

      H_trk_eff_pt_num = (TH1F*)filename->Get("MyProcess/PixelTracks/Trackertracks efficiency vs pt numerator");
      H_trk_eff_pt_den = (TH1F*)filename->Get("MyProcess/PixelTracks/Trackertracks efficiency vs pt denominator");
      H_trk_eff_eta_num = (TH1F*)filename->Get("MyProcess/PixelTracks/Trackertracks efficiency vs eta numerator");
      H_trk_eff_eta_den = (TH1F*)filename->Get("MyProcess/PixelTracks/Trackertracks efficiency vs eta denominator");

      H_trk_toteff_pt_num = (TH1F*)filename->Get("MyProcess/PixelTracks/Total trackertracks efficiency vs pt numerator");
      H_trk_toteff_pt_den = (TH1F*)filename->Get("MyProcess/PixelTracks/Total trackertracks efficiency vs pt denominator");
      H_trk_toteff_eta_num = (TH1F*)filename->Get("MyProcess/PixelTracks/Total trackertracks efficiency vs eta numerator");
      H_trk_toteff_eta_den = (TH1F*)filename->Get("MyProcess/PixelTracks/Total trackertracks efficiency vs eta denominator");

      for (int i = 0; i < 3; i++)
      	{
      	  stringstream i_str;
      	  i_str << i + 1;
      	  TString hName;

      	  hName = "MyProcess/PixelTracks/Trackertracks efficiency vs phi numerator thr" + i_str.str(); 
      	  H_trk_eff_phi_num[i] = (TH1F*)filename->Get(hName);

      	  hName = "MyProcess/PixelTracks/Trackertracks efficiency vs phi denominator thr" + i_str.str(); 
      	  H_trk_eff_phi_den[i] = (TH1F*)filename->Get(hName);

      	  hName = "MyProcess/PixelTracks/Total trackertracks efficiency vs phi numerator thr" + i_str.str(); 
      	  H_trk_toteff_phi_num[i] = (TH1F*)filename->Get(hName);

      	  hName = "MyProcess/PixelTracks/Total trackertracks efficiency vs phi denominator thr" + i_str.str(); 
      	  H_trk_toteff_phi_den[i] = (TH1F*)filename->Get(hName);
      	}

      for (int i = 0; i < 2; i++)
	{
	  stringstream i_str;
	  i_str << i + 1;
	  TString hName;

	  hName = "MyProcess/PixelTracks/Pixeltracks vs chi2DoF and pt cuts thr" + i_str.str(); 
	  H_pix_trk_normchi2_pt_eta[i] = (TH2F*)filename->Get(hName);

	  hName = "MyProcess/PixelTracks/Trackertracks vs chi2DoF and pt cuts thr" + i_str.str(); 
	  H_trk_trk_normchi2_pt_eta[i] = (TH2F*)filename->Get(hName);
	  
	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/Matched pixeltracks vs chi2DoF and pt cuts pur thr" + i_str.str(); 
	  H_trked_normchi2_pt_pur_eta[i] = (TH2F*)filename->Get(hName);

	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/Matched pixeltracks vs chi2DoF and pt cuts eff thr" + i_str.str(); 
	  H_trked_normchi2_pt_eff_eta[i] = (TH2F*)filename->Get(hName);
	}

      H_trked_eta_phi_pur = (TH2F*)filename->Get("MyProcess/PixelTracks/CommonPxTkPlots/Matched pixeltracks distribution in eta and phi pur");
      H_trked_eta_phi_eff = (TH2F*)filename->Get("MyProcess/PixelTracks/CommonPxTkPlots/Matched pixeltracks distribution in eta and phi eff");
      H_trked_eta_pt_pur = (TH2F*)filename->Get("MyProcess/PixelTracks/CommonPxTkPlots/Matched pixeltracks distribution in eta and pt pur");
      H_trked_eta_pt_eff = (TH2F*)filename->Get("MyProcess/PixelTracks/CommonPxTkPlots/Matched pixeltracks distribution in eta and pt eff");

      stringstream PtCutStr;					
      stringstream EtaCutStr;					
      stringstream HistoName;
      PtCutStr.precision(3);
      EtaCutStr.precision(3);
      HistoName.precision(3);
      
      double PtCut;
      double EtaCut;
      for (PtCut = PtStep; rint(PtCut*100.) <= rint(RangePt*100.); PtCut += PtStep)
	{
	  PtCutStr.str("");
	  PtCutStr << PtCut;

	  for (EtaCut = EtaStep; rint(EtaCut*100.) <= rint(RangeEta*100.); EtaCut += EtaStep)
	    {
	      EtaCutStr.str("");
	      EtaCutStr << EtaCut;
   	  	      
	      HistoName.str("");
	      HistoName << "MyProcess/TrackFit/d0_pulls_pt_" << PtCut << "_eta_" << EtaCut;
	      if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
		MapFit_d0_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());

	      HistoName.str("");
	      HistoName << "MyProcess/TrackFit/d0err_pt_" << PtCut << "_eta_" << EtaCut;
	      if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
		MapFit_d0err[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());

	      HistoName.str("");
	      HistoName << "MyProcess/TrackFit/d0_pt_" << PtCut << "_eta_" << EtaCut;
	      if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
		MapFit_d0[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());

	      HistoName.str("");
	      HistoName << "MyProcess/TrackFit/dz_pulls_pt_" << PtCut << "_eta_" << EtaCut;
	      if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
		MapFit_dz_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());

	      HistoName.str("");
	      HistoName << "MyProcess/TrackFit/dzerr_pt_" << PtCut << "_eta_" << EtaCut;
	      if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
		MapFit_dzerr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());

	      HistoName.str("");
	      HistoName << "MyProcess/TrackFit/dz_pt_" << PtCut << "_eta_" << EtaCut;
	      if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
		MapFit_dz[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());

	      HistoName.str("");
	      HistoName << "MyProcess/TrackFit/pt_pulls_pt_" << PtCut << "_eta_" << EtaCut;
	      if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
		MapFit_pt_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());

	      HistoName.str("");
	      HistoName << "MyProcess/TrackFit/pterr_pt_" << PtCut << "_eta_" << EtaCut;
	      if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
		MapFit_pterr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());

	      HistoName.str("");
	      HistoName << "MyProcess/TrackFit/pt_pt_" << PtCut << "_eta_" << EtaCut;
	      if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
		MapFit_pt[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());
	    }
	}

      // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      // @ Histograms containing ratios @
      // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      H_pix_trk_doubleCounter_pt_rel = (TH1F*)filename->Get("MyProcess/PixelTracks/Doublecounter pixeltracks pt relative distribution");
      H_pix_trk_doubleCounter_eta_rel = (TH1F*)filename->Get("MyProcess/PixelTracks/Doublecounter pixeltracks eta relative distribution");
      H_pix_trk_doubleCounter_pt_rel->Reset();
      H_pix_trk_doubleCounter_eta_rel->Reset();

      H_trk_trk_doubleCounter_pt_rel = (TH1F*)filename->Get("MyProcess/PixelTracks/Doublecounter trackertracks pt relative distribution");
      H_trk_trk_doubleCounter_eta_rel = (TH1F*)filename->Get("MyProcess/PixelTracks/Doublecounter trackertracks eta relative distribution");
      H_trk_trk_doubleCounter_pt_rel->Reset();
      H_trk_trk_doubleCounter_eta_rel->Reset();

      H_trk_eff_pt = (TH1F*)filename->Get("MyProcess/PixelTracks/Trackertracks efficiency vs pt");
      H_trk_eff_eta = (TH1F*)filename->Get("MyProcess/PixelTracks/Trackertracks efficiency vs eta");
      H_trk_eff_pt->Reset();
      H_trk_eff_eta->Reset();

      H_trk_toteff_pt = (TH1F*)filename->Get("MyProcess/PixelTracks/Total trackertracks efficiency vs pt");
      H_trk_toteff_eta = (TH1F*)filename->Get("MyProcess/PixelTracks/Total trackertracks efficiency vs eta");
      H_trk_toteff_pt->Reset();
      H_trk_toteff_eta->Reset();

      for (int i = 0; i < 3; i++)
      	{
      	  stringstream i_str;
      	  i_str << i + 1;
      	  TString hName;

      	  hName = "MyProcess/PixelTracks/Trackertracks efficiency vs phi thr" + i_str.str(); 
      	  H_trk_eff_phi[i] = (TH1F*)filename->Get(hName);
      	  H_trk_eff_phi[i]->Reset();

      	  hName = "MyProcess/PixelTracks/Total trackertracks efficiency vs phi thr" + i_str.str(); 
      	  H_trk_toteff_phi[i] = (TH1F*)filename->Get(hName);
      	  H_trk_toteff_phi[i]->Reset();
      	}

      // #################################
      // # Use this for constant binning #
      // #################################
//       const int nPtBins = rint(RangePt/PtStep);
//       double* xPtBins = NULL;
//       xPtBins = new double[nPtBins + 1];
//       xPtBins[0] = 0.; 
//       for (int i = 1; i < nPtBins; i++) xPtBins[i] = xPtBins[i-1] + PtStep;
//       xPtBins[nPtBins] = RangePt;

      // #####################################
      // # Use this for NON constant binning #
      // #####################################
      const int nPtBins = rint(MaxPtEquiBinning/PtStep) + 7;
      double* xPtBins = NULL;
      xPtBins = new double[nPtBins + 1];
      xPtBins[0] = 0.;
      for (int i = 1; i <= rint(MaxPtEquiBinning/PtStep); i++) xPtBins[i] = xPtBins[i-1] + PtStep;
      xPtBins[(int)rint(MaxPtEquiBinning/PtStep) + 1] = 3.5;
      xPtBins[(int)rint(MaxPtEquiBinning/PtStep) + 2] = 4.0;
      xPtBins[(int)rint(MaxPtEquiBinning/PtStep) + 3] = 4.5;
      xPtBins[(int)rint(MaxPtEquiBinning/PtStep) + 4] = 5.25;
      xPtBins[(int)rint(MaxPtEquiBinning/PtStep) + 5] = 6.0;
      xPtBins[(int)rint(MaxPtEquiBinning/PtStep) + 6] = 7.0;
      xPtBins[(int)rint(MaxPtEquiBinning/PtStep) + 7] = RangePt;

      for (int i = 0; i < 2; i++) 
	{
	  stringstream i_str;
	  i_str << i + 1;
	  TString hName;

	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs chi2DoF and pt cuts thr" + i_str.str(); 
	  H_purity_trked_normchi2_pt_eta[i] = (TH2F*)filename->Get(hName);
	  H_purity_trked_normchi2_pt_eta[i]->Reset();

	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs pt thr" + i_str.str();    
	  H_purity_trked_pt_eta[i] = (TH1F*)filename->Get(hName);
	  H_purity_trked_pt_eta[i]->Reset();

	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs chi2DoF and pt cuts thr" + i_str.str(); 
	  H_efficiency_trked_normchi2_pt_eta[i] = (TH2F*)filename->Get(hName);
	  H_efficiency_trked_normchi2_pt_eta[i]->Reset();

	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs pt thr" + i_str.str();    
	  H_efficiency_trked_pt_eta[i] = (TH1F*)filename->Get(hName);
	  H_efficiency_trked_pt_eta[i]->Reset();

	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs pt thr" + i_str.str();    
	  H_efficiency_trked_pt_eta_algo[i] = (TH1F*)filename->Get(hName);
	  H_efficiency_trked_pt_eta_algo[i]->Reset();

	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs pt thr" + i_str.str();
	  H_ptres_pt_eta[i] = (TH1F*)filename->Get(hName);
	  H_ptres_pt_eta[i]->Reset();
	  H_ptres_pt_eta[i] = (TH1F*)H_ptres_pt_eta[i]->Rebin(nPtBins, H_ptres_pt_eta[i]->GetName(), xPtBins);

	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs pt thr" + i_str.str();
	  H_ptpull_pt_eta[i] = (TH1F*)filename->Get(hName);
	  H_ptpull_pt_eta[i]->Reset();
	  H_ptpull_pt_eta[i] = (TH1F*)H_ptpull_pt_eta[i]->Rebin(nPtBins, H_ptpull_pt_eta[i]->GetName(), xPtBins);

	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs pt thr" + i_str.str();
	  H_dzres_pt_eta[i] = (TH1F*)filename->Get(hName);
	  H_dzres_pt_eta[i]->Reset();
 
	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs pt thr" + i_str.str();
	  H_dzpull_pt_eta[i] = (TH1F*)filename->Get(hName);
	  H_dzpull_pt_eta[i]->Reset();
	  
	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs pt thr" + i_str.str();
	  H_d0res_pt_eta[i] = (TH1F*)filename->Get(hName);
	  H_d0res_pt_eta[i]->Reset();

	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/d0 pulls vs pt thr" + i_str.str();
	  H_d0pull_pt_eta[i] = (TH1F*)filename->Get(hName);
	  H_d0pull_pt_eta[i]->Reset();

	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/SkewPixTrkpt vs pt thr" + i_str.str();
	  H_skew_pt_eta[i] = (TH1F*)filename->Get(hName);
	  H_skew_pt_eta[i]->Reset();
	}

      for (int i = 0; i < 5; i++) 
	{
	  stringstream i_str;
	  i_str << i + 1;
	  TString hName;

	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs eta thr" + i_str.str(); 
	  H_ptres_eta_pt[i] = (TH1F*)filename->Get(hName);
	  H_ptres_eta_pt[i]->Reset();

	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs eta thr" + i_str.str(); 
	  H_ptpull_eta_pt[i] = (TH1F*)filename->Get(hName);
	  H_ptpull_eta_pt[i]->Reset();
 
	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs eta thr" + i_str.str(); 
	  H_dzres_eta_pt[i] = (TH1F*)filename->Get(hName);
	  H_dzres_eta_pt[i]->Reset();

	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs eta thr" + i_str.str(); 
	  H_dzpull_eta_pt[i] = (TH1F*)filename->Get(hName);
	  H_dzpull_eta_pt[i]->Reset();

	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs eta thr" + i_str.str(); 
	  H_d0res_eta_pt[i] = (TH1F*)filename->Get(hName);
	  H_d0res_eta_pt[i]->Reset();

	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/d0 pulls vs eta thr" + i_str.str(); 
	  H_d0pull_eta_pt[i] = (TH1F*)filename->Get(hName);
	  H_d0pull_eta_pt[i]->Reset();
	}

      for (int i = 0; i < 3; i++) 
	{
	  stringstream i_str;
	  i_str << i + 1;
	  TString hName;
	  
	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs phi thr" + i_str.str();
	  H_purity_trked_phi[i] = (TH1F*)filename->Get(hName);
	  H_purity_trked_phi[i]->Reset();

	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs phi thr" + i_str.str();
	  H_efficiency_trked_phi[i] = (TH1F*)filename->Get(hName);
	  H_efficiency_trked_phi[i]->Reset();

	  hName = "MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs phi thr" + i_str.str();
	  H_efficiency_trked_phi_algo[i] = (TH1F*)filename->Get(hName);
	  H_efficiency_trked_phi_algo[i]->Reset();
	}

      H_purity_trked_eta = (TH1F*)filename->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs eta");
      H_purity_trked_pt_whole_eta = (TH1F*)filename->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs pt");
      H_purity_trked_eta->Reset();
      H_purity_trked_pt_whole_eta->Reset();

      H_efficiency_trked_eta = (TH1F*)filename->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs eta");
      H_efficiency_trked_pt_whole_eta = (TH1F*)filename->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs pt");
      H_efficiency_trked_eta->Reset();
      H_efficiency_trked_pt_whole_eta->Reset();

      H_efficiency_trked_eta_algo = (TH1F*)filename->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs eta");
      H_efficiency_trked_pt_whole_eta_algo = (TH1F*)filename->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs pt");
      H_efficiency_trked_eta_algo->Reset();
      H_efficiency_trked_pt_whole_eta_algo->Reset();

      H_tk_pur_eff = (TH1F*)filename->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity Efficiency Tracks");
      H_tk_pur_eff->Reset();

      return true;
    }
  else if (datatype.compare("Vx") == 0)
    {
      // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      // @ Histograms containing numerators or denominators @
      // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      H_vx_counters = (TH1F*)filename->Get("MyProcess/PixelVertices/CommonPxTkPlots/Counters Match Purity Efficiency for Vertices");

      stringstream String;
      stringstream HistoName;
      String.precision(3);
      HistoName.precision(3);
      
      for (int i = 1; i < 100; i++)
	{
	  String << i;
	  
	  HistoName.str("");
	  HistoName << "MyProcess/VxNtrks/x_diff_Ntrk_" << String.str();
	  if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
	    map_H_vx_x_diff[String.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());

	  HistoName.str("");
	  HistoName << "MyProcess/VxNtrks/y_diff_Ntrk_" << String.str();
	  if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
	    map_H_vx_y_diff[String.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());

	  HistoName.str("");
	  HistoName << "MyProcess/VxNtrks/z_diff_Ntrk_" << String.str();
	  if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
	    map_H_vx_z_diff[String.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());
	  
	  HistoName.str("");
	  HistoName << "MyProcess/VxNtrks/x_diff_pull_Ntrk_" << String.str();
	  if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
	    map_H_vx_x_diff_pull[String.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());

	  HistoName.str("");
	  HistoName << "MyProcess/VxNtrks/y_diff_pull_Ntrk_" << String.str();
	  if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
	    map_H_vx_y_diff_pull[String.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());

	  HistoName.str("");
	  HistoName << "MyProcess/VxNtrks/z_diff_pull_Ntrk_" << String.str();
	  if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
	    map_H_vx_z_diff_pull[String.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());

	  HistoName.str("");
	  HistoName << "MyProcess/VxNtrks/x_diff_orig_Ntrk_" << String.str();
	  if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
	    map_H_vx_x_diff_orig[String.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());

	  HistoName.str("");
	  HistoName << "MyProcess/VxNtrks/y_diff_orig_Ntrk_" << String.str();
	  if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
	    map_H_vx_y_diff_orig[String.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());

	  HistoName.str("");
	  HistoName << "MyProcess/VxNtrks/z_diff_orig_Ntrk_" << String.str();
	  if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
	    map_H_vx_z_diff_orig[String.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());

	  HistoName.str("");
	  HistoName << "MyProcess/VxNtrks/x_err2_Ntrk_" << String.str();
	  if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
	    map_H_vx_x_err2[String.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());

	  HistoName.str("");
	  HistoName << "MyProcess/VxNtrks/y_err2_Ntrk_" << String.str();
	  if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
	    map_H_vx_y_err2[String.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());

	  HistoName.str("");
	  HistoName << "MyProcess/VxNtrks/z_err2_Ntrk_" << String.str();
	  if ((TH1F*)filename->Get(HistoName.str().c_str()) != NULL)
	    map_H_vx_z_err2[String.str().c_str()] = (TH1F*)filename->Get(HistoName.str().c_str());

	  String.str(""); 
	}

      // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      // @ Histograms containing ratios @
      // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      H_vx_pur_eff = (TH1F*)filename->Get("MyProcess/PixelVertices/CommonPxTkPlots/Purity Efficiency Vertices");
      H_vx_x_diff_vs_trk = (TH1F*)filename->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position difference along X vs Ntrks");
      H_vx_y_diff_vs_trk = (TH1F*)filename->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position difference along Y vs Ntrks");
      H_vx_z_diff_vs_trk = (TH1F*)filename->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position difference along Z vs Ntrks");
      H_vx_x_diff_pull_vs_trk = (TH1F*)filename->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position pull along X vs Ntrks");
      H_vx_y_diff_pull_vs_trk = (TH1F*)filename->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position pull along Y vs Ntrks");
      H_vx_z_diff_pull_vs_trk = (TH1F*)filename->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position pull along Z vs Ntrks");
      H_vx_x_diff_orig_vs_trk = (TH1F*)filename->Get("MyProcess/PixelVertices/CommonPxTkPlots/Original vertex position difference along X vs Ntrks");
      H_vx_y_diff_orig_vs_trk = (TH1F*)filename->Get("MyProcess/PixelVertices/CommonPxTkPlots/Original vertex position difference along Y vs Ntrks");
      H_vx_z_diff_orig_vs_trk = (TH1F*)filename->Get("MyProcess/PixelVertices/CommonPxTkPlots/Original vertex position difference along Z vs Ntrks");
      H_vx_x_diff_noerr_vs_trk = (TH1F*)filename->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex resolution along X vs Ntrks");
      H_vx_y_diff_noerr_vs_trk = (TH1F*)filename->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex resolution along Y vs Ntrks");
      H_vx_z_diff_noerr_vs_trk = (TH1F*)filename->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex resolution along Z vs Ntrks");
      H_vx_x_diff_vs_trk->Reset();
      H_vx_y_diff_vs_trk->Reset();
      H_vx_z_diff_vs_trk->Reset();
      H_vx_x_diff_pull_vs_trk->Reset();
      H_vx_y_diff_pull_vs_trk->Reset();
      H_vx_z_diff_pull_vs_trk->Reset();
      H_vx_x_diff_orig_vs_trk->Reset();
      H_vx_y_diff_orig_vs_trk->Reset();
      H_vx_z_diff_orig_vs_trk->Reset();
      H_vx_x_diff_noerr_vs_trk->Reset();
      H_vx_y_diff_noerr_vs_trk->Reset();
      H_vx_z_diff_noerr_vs_trk->Reset();
      H_vx_pur_eff->Reset();
      
      return true;
    }

  return false;
}


void WriteHistos (string datatype, TFile* filename)
{
  if (datatype.compare("Tk") == 0)
    {
      filename->cd("MyProcess/TrackFit");
      
      for (myMap::iterator it1 = MapFit_d0_pull.begin(); it1 != MapFit_d0_pull.end(); it1++)
	for (map<string, TH1F*>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++) it2->second->Write(it2->second->GetName(), TObject::kWriteDelete);
      for (myMap::iterator it1 = MapFit_d0err.begin(); it1 != MapFit_d0err.end(); it1++)
	for (map<string, TH1F*>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++) it2->second->Write(it2->second->GetName(), TObject::kWriteDelete);
      for (myMap::iterator it1 = MapFit_d0.begin(); it1 != MapFit_d0.end(); it1++)
	for (map<string, TH1F*>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++) it2->second->Write(it2->second->GetName(), TObject::kWriteDelete);
      
      for (myMap::iterator it1 = MapFit_dz_pull.begin(); it1 != MapFit_dz_pull.end(); it1++)
	for (map<string, TH1F*>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++) it2->second->Write(it2->second->GetName(), TObject::kWriteDelete);
      for (myMap::iterator it1 = MapFit_dzerr.begin(); it1 != MapFit_dzerr.end(); it1++)
	for (map<string, TH1F*>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++) it2->second->Write(it2->second->GetName(), TObject::kWriteDelete);
      for (myMap::iterator it1 = MapFit_dz.begin(); it1 != MapFit_dz.end(); it1++)
	for (map<string, TH1F*>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++) it2->second->Write(it2->second->GetName(), TObject::kWriteDelete);
      
      for (myMap::iterator it1 = MapFit_pt_pull.begin(); it1 != MapFit_pt_pull.end(); it1++)
	for (map<string, TH1F*>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++) it2->second->Write(it2->second->GetName(), TObject::kWriteDelete);
      for (myMap::iterator it1 = MapFit_pterr.begin(); it1 != MapFit_pterr.end(); it1++)
	for (map<string, TH1F*>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++) it2->second->Write(it2->second->GetName(), TObject::kWriteDelete);
      for (myMap::iterator it1 = MapFit_pt.begin(); it1 != MapFit_pt.end(); it1++)
	for (map<string, TH1F*>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++) it2->second->Write(it2->second->GetName(), TObject::kWriteDelete);
      
      filename->cd("MyProcess/PixelTracks");
      
      H_pix_trk_doubleCounter_eta_rel->Write(H_pix_trk_doubleCounter_eta_rel->GetName(), TObject::kWriteDelete);
      H_pix_trk_doubleCounter_pt_rel->Write(H_pix_trk_doubleCounter_pt_rel->GetName(), TObject::kWriteDelete);
      
      H_trk_trk_doubleCounter_eta_rel->Write(H_trk_trk_doubleCounter_eta_rel->GetName(), TObject::kWriteDelete);
      H_trk_trk_doubleCounter_pt_rel->Write(H_trk_trk_doubleCounter_pt_rel->GetName(), TObject::kWriteDelete);

      H_trk_eff_pt->Write(H_trk_eff_pt->GetName(), TObject::kWriteDelete);
      H_trk_eff_eta->Write(H_trk_eff_eta->GetName(), TObject::kWriteDelete);
      for (int i = 0; i < 3; i++) H_trk_eff_phi[i]->Write(H_trk_eff_phi[i]->GetName(), TObject::kWriteDelete);
      
      H_trk_toteff_pt->Write(H_trk_toteff_pt->GetName(), TObject::kWriteDelete);
      H_trk_toteff_eta->Write(H_trk_toteff_eta->GetName(), TObject::kWriteDelete);
      for (int i = 0; i < 3; i++) H_trk_toteff_phi[i]->Write(H_trk_toteff_phi[i]->GetName(), TObject::kWriteDelete);

      filename->cd("MyProcess/PixelTracks/CommonPxTkPlots");

      for (int i = 0; i < 2; i++) 
	{
	  H_purity_trked_normchi2_pt_eta[i]->Write(H_purity_trked_normchi2_pt_eta[i]->GetName(), TObject::kWriteDelete);
	  H_purity_trked_pt_eta[i]->Write(H_purity_trked_pt_eta[i]->GetName(), TObject::kWriteDelete);

	  H_efficiency_trked_normchi2_pt_eta[i]->Write(H_efficiency_trked_normchi2_pt_eta[i]->GetName(), TObject::kWriteDelete);
	  H_efficiency_trked_pt_eta[i]->Write(H_efficiency_trked_pt_eta[i]->GetName(), TObject::kWriteDelete);
	  H_efficiency_trked_pt_eta_algo[i]->Write(H_efficiency_trked_pt_eta_algo[i]->GetName(), TObject::kWriteDelete);

	  H_ptres_pt_eta[i]->Write(H_ptres_pt_eta[i]->GetName(), TObject::kWriteDelete);
	  H_ptpull_pt_eta[i]->Write(H_ptpull_pt_eta[i]->GetName(), TObject::kWriteDelete);

	  H_dzres_pt_eta[i]->Write(H_dzres_pt_eta[i]->GetName(), TObject::kWriteDelete);
	  H_dzpull_pt_eta[i]->Write(H_dzpull_pt_eta[i]->GetName(), TObject::kWriteDelete);

	  H_d0res_pt_eta[i]->Write(H_d0res_pt_eta[i]->GetName(), TObject::kWriteDelete);
	  H_d0pull_pt_eta[i]->Write(H_d0pull_pt_eta[i]->GetName(), TObject::kWriteDelete);

	  H_skew_pt_eta[i]->Write(H_skew_pt_eta[i]->GetName(), TObject::kWriteDelete);
	}

      for (int i = 0; i < 5; i++) 
	{
	  H_ptres_eta_pt[i]->Write(H_ptres_eta_pt[i]->GetName(), TObject::kWriteDelete);
	  H_ptpull_eta_pt[i]->Write(H_ptpull_eta_pt[i]->GetName(), TObject::kWriteDelete);
	  
	  H_dzres_eta_pt[i]->Write(H_dzres_eta_pt[i]->GetName(), TObject::kWriteDelete);
	  H_dzpull_eta_pt[i]->Write(H_dzpull_eta_pt[i]->GetName(), TObject::kWriteDelete);
	  
	  H_d0res_eta_pt[i]->Write(H_d0res_eta_pt[i]->GetName(), TObject::kWriteDelete);
	  H_d0pull_eta_pt[i]->Write(H_d0pull_eta_pt[i]->GetName(), TObject::kWriteDelete);
	}

      for (int i = 0; i < 3; i++) 
	{
	  H_purity_trked_phi[i]->Write(H_purity_trked_phi[i]->GetName(), TObject::kWriteDelete);

	  H_efficiency_trked_phi[i]->Write(H_efficiency_trked_phi[i]->GetName(), TObject::kWriteDelete);
	  H_efficiency_trked_phi_algo[i]->Write(H_efficiency_trked_phi_algo[i]->GetName(), TObject::kWriteDelete);
	}

      H_purity_trked_eta->Write(H_purity_trked_eta->GetName(), TObject::kWriteDelete);
      H_purity_trked_pt_whole_eta->Write(H_purity_trked_pt_whole_eta->GetName(), TObject::kWriteDelete);

      H_efficiency_trked_eta->Write(H_efficiency_trked_eta->GetName(), TObject::kWriteDelete);
      H_efficiency_trked_pt_whole_eta->Write(H_efficiency_trked_pt_whole_eta->GetName(), TObject::kWriteDelete);

      H_efficiency_trked_eta_algo->Write(H_efficiency_trked_eta_algo->GetName(), TObject::kWriteDelete);
      H_efficiency_trked_pt_whole_eta_algo->Write(H_efficiency_trked_pt_whole_eta_algo->GetName(), TObject::kWriteDelete);

      H_tk_pur_eff->Write(H_tk_pur_eff->GetName(), TObject::kWriteDelete);
    }
  else if (datatype.compare("Vx") == 0)
    {
      filename->cd("MyProcess/VxNtrks");

      for (map<string, TH1F*>::iterator it = map_H_vx_x_diff.begin(); it != map_H_vx_x_diff.end(); it++) it->second->Write(it->second->GetName(), TObject::kWriteDelete);
      for (map<string, TH1F*>::iterator it = map_H_vx_y_diff.begin(); it != map_H_vx_y_diff.end(); it++) it->second->Write(it->second->GetName(), TObject::kWriteDelete);
      for (map<string, TH1F*>::iterator it = map_H_vx_z_diff.begin(); it != map_H_vx_z_diff.end(); it++) it->second->Write(it->second->GetName(), TObject::kWriteDelete);
	  
      for (map<string, TH1F*>::iterator it = map_H_vx_x_diff_pull.begin(); it != map_H_vx_x_diff_pull.end(); it++) it->second->Write(it->second->GetName(), TObject::kWriteDelete);
      for (map<string, TH1F*>::iterator it = map_H_vx_y_diff_pull.begin(); it != map_H_vx_y_diff_pull.end(); it++) it->second->Write(it->second->GetName(), TObject::kWriteDelete);
      for (map<string, TH1F*>::iterator it = map_H_vx_z_diff_pull.begin(); it != map_H_vx_z_diff_pull.end(); it++) it->second->Write(it->second->GetName(), TObject::kWriteDelete);

      for (map<string, TH1F*>::iterator it = map_H_vx_x_diff_orig.begin(); it != map_H_vx_x_diff_orig.end(); it++) it->second->Write(it->second->GetName(), TObject::kWriteDelete);
      for (map<string, TH1F*>::iterator it = map_H_vx_y_diff_orig.begin(); it != map_H_vx_y_diff_orig.end(); it++) it->second->Write(it->second->GetName(), TObject::kWriteDelete);
      for (map<string, TH1F*>::iterator it = map_H_vx_z_diff_orig.begin(); it != map_H_vx_z_diff_orig.end(); it++) it->second->Write(it->second->GetName(), TObject::kWriteDelete);

      filename->cd("MyProcess/PixelVertices/CommonPxTkPlots");
      
      H_vx_x_diff_vs_trk->Write(H_vx_x_diff_vs_trk->GetName(), TObject::kWriteDelete);
      H_vx_y_diff_vs_trk->Write(H_vx_y_diff_vs_trk->GetName(), TObject::kWriteDelete);
      H_vx_z_diff_vs_trk->Write(H_vx_z_diff_vs_trk->GetName(), TObject::kWriteDelete);
      H_vx_x_diff_pull_vs_trk->Write(H_vx_x_diff_pull_vs_trk->GetName(), TObject::kWriteDelete);
      H_vx_y_diff_pull_vs_trk->Write(H_vx_y_diff_pull_vs_trk->GetName(), TObject::kWriteDelete);
      H_vx_z_diff_pull_vs_trk->Write(H_vx_z_diff_pull_vs_trk->GetName(), TObject::kWriteDelete);
      H_vx_x_diff_orig_vs_trk->Write(H_vx_x_diff_orig_vs_trk->GetName(), TObject::kWriteDelete);
      H_vx_y_diff_orig_vs_trk->Write(H_vx_y_diff_orig_vs_trk->GetName(), TObject::kWriteDelete);
      H_vx_z_diff_orig_vs_trk->Write(H_vx_z_diff_orig_vs_trk->GetName(), TObject::kWriteDelete);
      H_vx_x_diff_noerr_vs_trk->Write(H_vx_x_diff_noerr_vs_trk->GetName(), TObject::kWriteDelete);
      H_vx_y_diff_noerr_vs_trk->Write(H_vx_y_diff_noerr_vs_trk->GetName(), TObject::kWriteDelete);
      H_vx_z_diff_noerr_vs_trk->Write(H_vx_z_diff_noerr_vs_trk->GetName(), TObject::kWriteDelete);
      H_vx_pur_eff->Write(H_vx_pur_eff->GetName(), TObject::kWriteDelete);
    }
}


bool FitAlgorithm (TF1* FitFcn, TH1F* Histo)
{
  double constant;
  double mean;
  double sigma;

  constant = Histo->GetMaximum();
  mean     = Histo->GetMean();
  sigma    = Histo->GetRMS();
  
  FitFcn->SetParameters(constant, mean, sigma);
  Histo->Fit("Gauss", "QR");

  if ((isnan(FitFcn->GetParameter(2)) == false) && (isnan(FitFcn->GetParError(2)) == false)) return true;
  return false;
}


void FitAndFillMaps (string mapname,
		     double rangemin,
		     double rangemax)
{
  double DeconvolvErrsq;
  double ErrDeconvolvErrsq;
  stringstream PtCutStr;					
  stringstream EtaCutStr;					
  PtCutStr.precision(3);
  EtaCutStr.precision(3);
  double Val, Err;
  
  TF1* Gauss = new TF1("Gauss", "gaus", rangemin, rangemax);
 
  if (mapname.compare("d0pulls") == 0) // @@@ d0pulls MAP @@@
    {
      TH1F* d0pulls[5];
      d0pulls[0] = new TH1F("d0pulls_0", "d0pulls_0", 100, -pullsRange/2., pullsRange/2.);
      d0pulls[1] = new TH1F("d0pulls_1", "d0pulls_1", 100, -pullsRange/2., pullsRange/2.);
      d0pulls[2] = new TH1F("d0pulls_2", "d0pulls_2", 100, -pullsRange/2., pullsRange/2.);
      d0pulls[3] = new TH1F("d0pulls_3", "d0pulls_3", 100, -pullsRange/2., pullsRange/2.);
      d0pulls[4] = new TH1F("d0pulls_4", "d0pulls_4", 100, -pullsRange/2., pullsRange/2.);
      
      for (double PtCut = MinPtTk+PtStep; rint(PtCut*100.) <= rint(RangePt*100.); PtCut += PtStep)
	{
	  PtCutStr.str("");
	  PtCutStr << PtCut;
	  
	  for (double EtaCut = EtaStep; rint(EtaCut*100.) <= rint(RangeEta*100.); EtaCut += EtaStep)
	    {
	      EtaCutStr.str("");
	      EtaCutStr << EtaCut;

	      if ((MapFit_d0_pull.find(PtCutStr.str().c_str()) != MapFit_d0_pull.end()) &&
		  (MapFit_d0_pull[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) != MapFit_d0_pull[PtCutStr.str().c_str()].end()))
		{
		  for (int i = 0; i < 2; i++)
		    if ((rint(EtaCut*100.) > rint(EtaThrBin[i]*100.)) && (rint(EtaCut*100.) <= rint(EtaThrBin[i+1]*100.)))
		      d0pulls[i]->Add(MapFit_d0_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]);
		}
	    }
	  
	  for (int i = 0; i < 2; i++)
	    {
	      if (FitAlgorithm (Gauss,d0pulls[i]) == true)
		{
		  H_d0pull_pt_eta[i]->SetBinContent(H_d0pull_pt_eta[i]->GetXaxis()->FindBin(rint((PtCut-PtStep)*100.)/100.), Gauss->GetParameter(2));
		  H_d0pull_pt_eta[i]->SetBinError(H_d0pull_pt_eta[i]->GetXaxis()->FindBin(rint((PtCut-PtStep)*100.)/100.), Gauss->GetParError(2));
		}

	      d0pulls[i]->Reset();
	    }
	}
      
      for (double EtaCut = EtaStep; rint(EtaCut*100.) <= rint(RangeEta*100.); EtaCut += EtaStep)
	{
	  EtaCutStr.str("");		
	  EtaCutStr << EtaCut;
	  
	  for (double PtCut = MinPtTk+PtStep; rint(PtCut*100.) <= rint(RangePt*100.); PtCut += PtStep)
	    {
	      PtCutStr.str("");
	      PtCutStr << PtCut;

	      if ((MapFit_d0_pull.find(PtCutStr.str().c_str()) != MapFit_d0_pull.end()) &&
		  (MapFit_d0_pull[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) != MapFit_d0_pull[PtCutStr.str().c_str()].end()))
		{
		  for (int i = 0; i < 5; i++)
		    if ((rint(PtCut*100.) > rint(PtThr[i]*100.)) && (rint(PtCut*100.) <= rint(PtThr[i+1]*100.)))
		      d0pulls[i]->Add(MapFit_d0_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]);
		}
	    }

	  for (int i = 0; i < 5; i++)
	    {
	      if (FitAlgorithm (Gauss,d0pulls[i]) == true)
		{
		  H_d0pull_eta_pt[i]->SetBinContent(H_d0pull_eta_pt[i]->GetXaxis()->FindBin(rint((EtaCut-EtaStep)*100.)/100.), Gauss->GetParameter(2));
		  H_d0pull_eta_pt[i]->SetBinError(H_d0pull_eta_pt[i]->GetXaxis()->FindBin(rint((EtaCut-EtaStep)*100.)/100.), Gauss->GetParError(2));
		}

	      d0pulls[i]->Reset();
	    }
	}

      for (int i = 0; i < 5; i++) delete d0pulls[i];
    }
  else if (mapname.compare("d0") == 0) // @@@ d0 MAP @@@
    {
      TH1F* d0[5];
      d0[0] = new TH1F("d0_0", "d0_0", 100, -d0Range/2., d0Range/2.);
      d0[1] = new TH1F("d0_1", "d0_1", 100, -d0Range/2., d0Range/2.);
      d0[2] = new TH1F("d0_2", "d0_2", 100, -d0Range/2., d0Range/2.);
      d0[3] = new TH1F("d0_3", "d0_3", 100, -d0Range/2., d0Range/2.);
      d0[4] = new TH1F("d0_4", "d0_4", 100, -d0Range/2., d0Range/2.);
      
      TH1F* d0err[5];
      d0err[0] = new TH1F("d0err_0", "d0err_0", 200, 0., d0errRange);
      d0err[1] = new TH1F("d0err_1", "d0err_1", 200, 0., d0errRange);
      d0err[2] = new TH1F("d0err_2", "d0err_2", 200, 0., d0errRange);
      d0err[3] = new TH1F("d0err_3", "d0err_3", 200, 0., d0errRange);
      d0err[4] = new TH1F("d0err_4", "d0err_4", 200, 0., d0errRange);

      for (double PtCut = MinPtTk+PtStep; rint(PtCut*100.) <= rint(RangePt*100.); PtCut += PtStep)
	{
	  PtCutStr.str("");
	  PtCutStr << PtCut;

	  for (double EtaCut = EtaStep; rint(EtaCut*100.) <= rint(RangeEta*100.); EtaCut += EtaStep)
	    {
	      EtaCutStr.str("");
	      EtaCutStr << EtaCut;

	      if ((MapFit_d0.find(PtCutStr.str().c_str()) != MapFit_d0.end()) &&
		  (MapFit_d0[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) != MapFit_d0[PtCutStr.str().c_str()].end()) &&
		  (MapFit_d0err.find(PtCutStr.str().c_str()) != MapFit_d0err.end()) &&
		  (MapFit_d0err[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) != MapFit_d0err[PtCutStr.str().c_str()].end()))
		{
		  for (int i = 0; i < 2; i++)
		    if ((rint(EtaCut*100.) > rint(EtaThrBin[i]*100.)) && (rint(EtaCut*100.) <= rint(EtaThrBin[i+1]*100.)))
		      {
			d0[i]->Add(MapFit_d0[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]);
			d0err[i]->Add(MapFit_d0err[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]);
		      }
		}
	    }
	  
	  for (int i = 0; i < 2; i++)
	    {
	      DeconvolvErrsq    = d0err[i]->GetMean() / 100000000.;
	      ErrDeconvolvErrsq = d0err[i]->GetRMS() / sqrt(d0err[i]->GetEntries()) / 100000000.;
	      
	      if ((FitAlgorithm (Gauss,d0[i]) == true) && (Gauss->GetParameter(2)*Gauss->GetParameter(2) > DeconvolvErrsq))
		{
		  Val = sqrt(Gauss->GetParameter(2)*Gauss->GetParameter(2) - DeconvolvErrsq);
		  Err = 1./(2.*Val) * sqrt(powf(2.*Gauss->GetParameter(2)*Gauss->GetParError(2),2.) + powf(ErrDeconvolvErrsq,2.));
		  
		  H_d0res_pt_eta[i]->SetBinContent(H_d0res_pt_eta[i]->GetXaxis()->FindBin(rint((PtCut-PtStep)*100.)/100.), Val);
		  H_d0res_pt_eta[i]->SetBinError(H_d0res_pt_eta[i]->GetXaxis()->FindBin(rint((PtCut-PtStep)*100.)/100.), Err);
		}

	      d0[i]->Reset();
	      d0err[i]->Reset();
	    }
	}
      
      for (double EtaCut = EtaStep; rint(EtaCut*100.) <= rint(RangeEta*100.); EtaCut += EtaStep)
	{
	  EtaCutStr.str("");		
	  EtaCutStr << EtaCut;

	  for (double PtCut = MinPtTk+PtStep; rint(PtCut*100.) <= rint(RangePt*100.); PtCut += PtStep)
	    {
	      PtCutStr.str("");
	      PtCutStr << PtCut;

	      if ((MapFit_d0.find(PtCutStr.str().c_str()) != MapFit_d0.end()) &&
		  (MapFit_d0[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) != MapFit_d0[PtCutStr.str().c_str()].end()) &&
		  (MapFit_d0err.find(PtCutStr.str().c_str()) != MapFit_d0err.end()) &&
		  (MapFit_d0err[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) != MapFit_d0err[PtCutStr.str().c_str()].end()))
		{
		  for (int i = 0; i < 5; i++)
		    if ((rint(PtCut*100.) > rint(PtThr[i]*100.)) && (rint(PtCut*100.) <= rint(PtThr[i+1]*100.)))
		      {
			d0[i]->Add(MapFit_d0[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]);
			d0err[i]->Add(MapFit_d0err[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]);
		      }
		}
	    }

	  for (int i = 0; i < 5; i++)
	    {
	      DeconvolvErrsq    = d0err[i]->GetMean() / 100000000.;
	      ErrDeconvolvErrsq = d0err[i]->GetRMS() / sqrt(d0err[i]->GetEntries()) / 100000000.;

	      if ((FitAlgorithm (Gauss,d0[i]) == true) && (Gauss->GetParameter(2)*Gauss->GetParameter(2) > DeconvolvErrsq))
		{
		  Val = sqrt(Gauss->GetParameter(2)*Gauss->GetParameter(2) - DeconvolvErrsq);
		  Err = 1./(2.*Val) * sqrt(powf(2.*Gauss->GetParameter(2)*Gauss->GetParError(2),2.) + powf(ErrDeconvolvErrsq,2.));
		  
		  H_d0res_eta_pt[i]->SetBinContent(H_d0res_eta_pt[i]->GetXaxis()->FindBin(rint((EtaCut-EtaStep)*100.)/100.), Val);
		  H_d0res_eta_pt[i]->SetBinError(H_d0res_eta_pt[i]->GetXaxis()->FindBin(rint((EtaCut-EtaStep)*100.)/100.), Err);
		}

	      d0[i]->Reset();
	      d0err[i]->Reset();
	    }
	}

      for (int i = 0; i < 5; i++) { delete d0[i]; delete d0err[i]; }
    }
  else if (mapname.compare("dzpulls") == 0) // @@@ dzpulls MAP @@@
    {
      TH1F* dzpulls[5];
      dzpulls[0] = new TH1F("dzpulls_0", "dzpulls_0", 100, -pullsRange/2., pullsRange/2.);
      dzpulls[1] = new TH1F("dzpulls_1", "dzpulls_1", 100, -pullsRange/2., pullsRange/2.);
      dzpulls[2] = new TH1F("dzpulls_2", "dzpulls_2", 100, -pullsRange/2., pullsRange/2.);
      dzpulls[3] = new TH1F("dzpulls_3", "dzpulls_3", 100, -pullsRange/2., pullsRange/2.);
      dzpulls[4] = new TH1F("dzpulls_4", "dzpulls_4", 100, -pullsRange/2., pullsRange/2.);

      for (double PtCut = MinPtTk+PtStep; rint(PtCut*100.) <= rint(RangePt*100.); PtCut += PtStep)
	{
	  PtCutStr.str("");
	  PtCutStr << PtCut;
	  
	  for (double EtaCut = EtaStep; rint(EtaCut*100.) <= rint(RangeEta*100.); EtaCut += EtaStep)
	    {
	      EtaCutStr.str("");
	      EtaCutStr << EtaCut;

	      if ((MapFit_dz_pull.find(PtCutStr.str().c_str()) != MapFit_dz_pull.end()) &&
		  (MapFit_dz_pull[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) != MapFit_dz_pull[PtCutStr.str().c_str()].end()))
		{
		  for (int i = 0; i < 2; i++)
		    if ((rint(EtaCut*100.) > rint(EtaThrBin[i]*100.)) && (rint(EtaCut*100.) <= rint(EtaThrBin[i+1]*100.)))
		      dzpulls[i]->Add(MapFit_dz_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]);
		}
	    }
	  
	  for (int i = 0; i < 2; i++)
	    {
	      if (FitAlgorithm (Gauss,dzpulls[i]) == true)
		{
		  H_dzpull_pt_eta[i]->SetBinContent(H_dzpull_pt_eta[i]->GetXaxis()->FindBin(rint((PtCut-PtStep)*100.)/100.), Gauss->GetParameter(2));
		  H_dzpull_pt_eta[i]->SetBinError(H_dzpull_pt_eta[i]->GetXaxis()->FindBin(rint((PtCut-PtStep)*100.)/100.), Gauss->GetParError(2));
		}

	      dzpulls[i]->Reset();
	    }
	}

      for (double EtaCut = EtaStep; rint(EtaCut*100.) <= rint(RangeEta*100.); EtaCut += EtaStep)
	{
	  EtaCutStr.str("");		
	  EtaCutStr << EtaCut;
	  
	  for (double PtCut = MinPtTk+PtStep; rint(PtCut*100.) <= rint(RangePt*100.); PtCut += PtStep)
	    {
	      PtCutStr.str("");
	      PtCutStr << PtCut;

	      if ((MapFit_dz_pull.find(PtCutStr.str().c_str()) != MapFit_dz_pull.end()) &&
		  (MapFit_dz_pull[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) != MapFit_dz_pull[PtCutStr.str().c_str()].end()))
		{
		  for (int i = 0; i < 5; i++)
		    if ((rint(PtCut*100.) > rint(PtThr[i]*100.)) && (rint(PtCut*100.) <= rint(PtThr[i+1]*100.)))
		      dzpulls[i]->Add(MapFit_dz_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]);
		}
	    }

	  for (int i = 0; i < 5; i++)
	    {
	      if (FitAlgorithm (Gauss,dzpulls[i]) == true)
		{
		  H_dzpull_eta_pt[i]->SetBinContent(H_dzpull_eta_pt[i]->GetXaxis()->FindBin(rint((EtaCut-EtaStep)*100.)/100.), Gauss->GetParameter(2));
		  H_dzpull_eta_pt[i]->SetBinError(H_dzpull_eta_pt[i]->GetXaxis()->FindBin(rint((EtaCut-EtaStep)*100.)/100.), Gauss->GetParError(2));
		}

	      dzpulls[i]->Reset();
	    }
	}
      
      for (int i = 0; i < 5; i++) delete dzpulls[i];
    }  
  else if (mapname.compare("dz") == 0) // @@@ dz MAP @@@
    {
      TH1F* dz[5];
      dz[0] = new TH1F("dz_0", "dz_0", 100, -dzRange/2., dzRange/2.);
      dz[1] = new TH1F("dz_1", "dz_1", 100, -dzRange/2., dzRange/2.);
      dz[2] = new TH1F("dz_2", "dz_2", 100, -dzRange/2., dzRange/2.);
      dz[3] = new TH1F("dz_3", "dz_3", 100, -dzRange/2., dzRange/2.);
      dz[4] = new TH1F("dz_4", "dz_4", 100, -dzRange/2., dzRange/2.);

      TH1F* dzerr[5];
      dzerr[0] = new TH1F("dzerr_0", "dzerr_0", 200, 0., dzerrRange);
      dzerr[1] = new TH1F("dzerr_1", "dzerr_1", 200, 0., dzerrRange);
      dzerr[2] = new TH1F("dzerr_2", "dzerr_2", 200, 0., dzerrRange);
      dzerr[3] = new TH1F("dzerr_3", "dzerr_3", 200, 0., dzerrRange);
      dzerr[4] = new TH1F("dzerr_4", "dzerr_4", 200, 0., dzerrRange);

      for (double PtCut = MinPtTk+PtStep; rint(PtCut*100.) <= rint(RangePt*100.); PtCut += PtStep)
	{
	  PtCutStr.str("");
	  PtCutStr << PtCut;

	  for (double EtaCut = EtaStep; rint(EtaCut*100.) <= rint(RangeEta*100.); EtaCut += EtaStep)
	    {
	      EtaCutStr.str("");
	      EtaCutStr << EtaCut;

	      if ((MapFit_dz.find(PtCutStr.str().c_str()) != MapFit_dz.end()) &&
		  (MapFit_dz[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) != MapFit_dz[PtCutStr.str().c_str()].end()) &&
		  (MapFit_dzerr.find(PtCutStr.str().c_str()) != MapFit_dzerr.end()) &&
		  (MapFit_dzerr[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) != MapFit_dzerr[PtCutStr.str().c_str()].end()))
		{
		  for (int i = 0; i < 2; i++)
		    if ((rint(EtaCut*100.) > rint(EtaThrBin[i]*100.)) && (rint(EtaCut*100.) <= rint(EtaThrBin[i+1]*100.)))
		      {
			dz[i]->Add(MapFit_dz[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]);
			dzerr[i]->Add(MapFit_dzerr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]);
		      }
		}
	    }
	  
	  for (int i = 0; i < 2; i++)
	    {
	      DeconvolvErrsq    = dzerr[i]->GetMean() / 100000000.;
	      ErrDeconvolvErrsq = dzerr[i]->GetRMS() / sqrt(dzerr[i]->GetEntries()) / 100000000.;
	      
	      if ((FitAlgorithm (Gauss,dz[i]) == true) && (Gauss->GetParameter(2)*Gauss->GetParameter(2) > DeconvolvErrsq))
		{
		  Val = sqrt(Gauss->GetParameter(2)*Gauss->GetParameter(2) - DeconvolvErrsq);
		  Err = 1./(2.*Val) * sqrt(powf(2.*Gauss->GetParameter(2)*Gauss->GetParError(2),2.) + powf(ErrDeconvolvErrsq,2.));
		  
		  H_dzres_pt_eta[i]->SetBinContent(H_dzres_pt_eta[i]->GetXaxis()->FindBin(rint((PtCut-PtStep)*100.)/100.), Val);
		  H_dzres_pt_eta[i]->SetBinError(H_dzres_pt_eta[i]->GetXaxis()->FindBin(rint((PtCut-PtStep)*100.)/100.), Err);
		}

	      dz[i]->Reset();
	      dzerr[i]->Reset();
	    }
	}

      for (double EtaCut = EtaStep; rint(EtaCut*100.) <= rint(RangeEta*100.); EtaCut += EtaStep)
	{
	  EtaCutStr.str("");		
	  EtaCutStr << EtaCut;

	  for (double PtCut = MinPtTk+PtStep; rint(PtCut*100.) <= rint(RangePt*100.); PtCut += PtStep)
	    {
	      PtCutStr.str("");
	      PtCutStr << PtCut;

	      if ((MapFit_dz.find(PtCutStr.str().c_str()) != MapFit_dz.end()) &&
		  (MapFit_dz[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) != MapFit_dz[PtCutStr.str().c_str()].end()) &&
		  (MapFit_dzerr.find(PtCutStr.str().c_str()) != MapFit_dzerr.end()) &&
		  (MapFit_dzerr[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) != MapFit_dzerr[PtCutStr.str().c_str()].end()))
		{
		  for (int i = 0; i < 5; i++)
		    if ((rint(PtCut*100.) > rint(PtThr[i]*100.)) && (rint(PtCut*100.) <= rint(PtThr[i+1]*100.)))
		      {
			dz[i]->Add(MapFit_dz[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]);
			dzerr[i]->Add(MapFit_dzerr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]);
		      }
		}
	    }

	  for (int i = 0; i < 5; i++)
	    {
	      DeconvolvErrsq    = dzerr[i]->GetMean() / 100000000.;
	      ErrDeconvolvErrsq = dzerr[i]->GetRMS() / sqrt(dzerr[i]->GetEntries()) / 100000000.;

	      if ((FitAlgorithm (Gauss,dz[i]) == true) && (Gauss->GetParameter(2)*Gauss->GetParameter(2) > DeconvolvErrsq))
		{
		  Val = sqrt(Gauss->GetParameter(2)*Gauss->GetParameter(2) - DeconvolvErrsq);
		  Err = 1./(2.*Val) * sqrt(powf(2.*Gauss->GetParameter(2)*Gauss->GetParError(2),2.) + powf(ErrDeconvolvErrsq,2.));
		  
		  H_dzres_eta_pt[i]->SetBinContent(H_dzres_eta_pt[i]->GetXaxis()->FindBin(rint((EtaCut-EtaStep)*100.)/100.), Val);
		  H_dzres_eta_pt[i]->SetBinError(H_dzres_eta_pt[i]->GetXaxis()->FindBin(rint((EtaCut-EtaStep)*100.)/100.), Err);
		}

	      dz[i]->Reset();
	      dzerr[i]->Reset();
	    }
	}

      for (int i = 0; i < 5; i++) { delete dz[i]; delete dzerr[i]; }
    }
  else if (mapname.compare("ptpulls") == 0) // @@@ ptpulls MAP @@@
    {
      TH1F* ptpulls[5];
      ptpulls[0] = new TH1F("ptpulls_0", "ptpulls_0", 100, -pullsRange/2., pullsRange/2.);
      ptpulls[1] = new TH1F("ptpulls_1", "ptpulls_1", 100, -pullsRange/2., pullsRange/2.);
      ptpulls[2] = new TH1F("ptpulls_2", "ptpulls_2", 100, -pullsRange/2., pullsRange/2.);
      ptpulls[3] = new TH1F("ptpulls_3", "ptpulls_3", 100, -pullsRange/2., pullsRange/2.);
      ptpulls[4] = new TH1F("ptpulls_4", "ptpulls_4", 100, -pullsRange/2., pullsRange/2.);

      for (double PtCut = MinPtTk+PtStep; rint(PtCut*100.) <= rint(RangePt*100.); PtCut += PtStep)
	{
	  PtCutStr.str("");
	  PtCutStr << PtCut;	  

	  for (double EtaCut = EtaStep; rint(EtaCut*100.) <= rint(RangeEta*100.); EtaCut += EtaStep)
	    {
	      EtaCutStr.str("");
	      EtaCutStr << EtaCut;

	      if ((MapFit_pt_pull.find(PtCutStr.str().c_str()) != MapFit_pt_pull.end()) &&
		  (MapFit_pt_pull[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) != MapFit_pt_pull[PtCutStr.str().c_str()].end()))
		{
		  for (int i = 0; i < 2; i++)
		    if ((rint(EtaCut*100.) > rint(EtaThrBin[i]*100.)) && (rint(EtaCut*100.) <= rint(EtaThrBin[i+1]*100.)))
		      ptpulls[i]->Add(MapFit_pt_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]);
		}
	    }

      	  if ((rint(PtCut*100.) == rint(RangePt*100.)) || (H_ptpull_pt_eta[0]->GetXaxis()->FindBin(rint((PtCut-PtStep)*100.)/100.) != H_ptpull_pt_eta[0]->GetXaxis()->FindBin(rint(PtCut*100.)/100.)))
	    {
	      for (int i = 0; i < 2; i++)
		{
		  if (FitAlgorithm (Gauss,ptpulls[i]) == true)
		    {
		      H_ptpull_pt_eta[i]->SetBinContent(H_ptpull_pt_eta[i]->GetXaxis()->FindBin(rint((PtCut-PtStep)*100.)/100.), Gauss->GetParameter(2));
		      H_ptpull_pt_eta[i]->SetBinError(H_ptpull_pt_eta[i]->GetXaxis()->FindBin(rint((PtCut-PtStep)*100.)/100.), Gauss->GetParError(2));
		    }

		  ptpulls[i]->Reset();
		}
	    }
	}
      
      for (double EtaCut = EtaStep; rint(EtaCut*100.) <= rint(RangeEta*100.); EtaCut += EtaStep)
	{
	  EtaCutStr.str("");
	  EtaCutStr << EtaCut;
 
	  for (double PtCut = MinPtTk+PtStep; rint(PtCut*100.) <= rint(RangePt*100.); PtCut += PtStep)
	    {
	      PtCutStr.str("");
	      PtCutStr << PtCut;
 
	      if ((MapFit_pt_pull.find(PtCutStr.str().c_str()) != MapFit_pt_pull.end()) &&
		  (MapFit_pt_pull[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) != MapFit_pt_pull[PtCutStr.str().c_str()].end()))
		{
		  for (int i = 0; i < 5; i++)
		    if ((rint(PtCut*100.) > rint(PtThr[i]*100.)) && (rint(PtCut*100.) <= rint(PtThr[i+1]*100.)))
		      ptpulls[i]->Add(MapFit_pt_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]);
		}
	    }

	  for (int i = 0; i < 5; i++)
	    {
	      if (FitAlgorithm (Gauss,ptpulls[i]) == true)
		{
		  H_ptpull_eta_pt[i]->SetBinContent(H_ptpull_eta_pt[i]->GetXaxis()->FindBin(rint((EtaCut-EtaStep)*100.)/100.), Gauss->GetParameter(2));
		  H_ptpull_eta_pt[i]->SetBinError(H_ptpull_eta_pt[i]->GetXaxis()->FindBin(rint((EtaCut-EtaStep)*100.)/100.), Gauss->GetParError(2));
		}

	      ptpulls[i]->Reset();
	    }
	}
      
      for (int i = 0; i < 5; i++) delete ptpulls[i];
    }
  else if (mapname.compare("pt") == 0) // @@@ pt MAP @@@
    {
      TH1F* pt[5];
      pt[0] = new TH1F("pt_0", "pt_0", 600, -ptRange/2., ptRange/2.);
      pt[1] = new TH1F("pt_1", "pt_1", 600, -ptRange/2., ptRange/2.);
      pt[2] = new TH1F("pt_2", "pt_2", 600, -ptRange/2., ptRange/2.);
      pt[3] = new TH1F("pt_3", "pt_3", 600, -ptRange/2., ptRange/2.);
      pt[4] = new TH1F("pt_4", "pt_4", 600, -ptRange/2., ptRange/2.);

      TH1F* pterr[5];
      pterr[0] = new TH1F("pterr_0", "pterr_0", 5000, 0., pterrRange);
      pterr[1] = new TH1F("pterr_1", "pterr_1", 5000, 0., pterrRange);
      pterr[2] = new TH1F("pterr_2", "pterr_2", 5000, 0., pterrRange);
      pterr[3] = new TH1F("pterr_3", "pterr_3", 5000, 0., pterrRange);
      pterr[4] = new TH1F("pterr_4", "pterr_4", 5000, 0., pterrRange);

      double OldLowEdgeBin = MinPtTk;

      for (double PtCut = MinPtTk+PtStep; rint(PtCut*100.) <= rint(RangePt*100.); PtCut += PtStep)
	{
	  PtCutStr.str("");
	  PtCutStr << PtCut;
	  
	  for (double EtaCut = EtaStep; rint(EtaCut*100.) <= rint(RangeEta*100.); EtaCut += EtaStep)
	    {
	      EtaCutStr.str("");
	      EtaCutStr << EtaCut;

	      if ((MapFit_pt.find(PtCutStr.str().c_str()) != MapFit_pt.end()) &&
		  (MapFit_pt[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) != MapFit_pt[PtCutStr.str().c_str()].end()) &&
		  (MapFit_pterr.find(PtCutStr.str().c_str()) != MapFit_pterr.end()) &&
		  (MapFit_pterr[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) != MapFit_pterr[PtCutStr.str().c_str()].end()))
		{
		  for (int i = 0; i < 2; i++)
		    if ((rint(EtaCut*100.) > rint(EtaThrBin[i]*100.)) && (rint(EtaCut*100.) <= rint(EtaThrBin[i+1]*100.)))
		      {
			pt[i]->Add(MapFit_pt[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]);
			pterr[i]->Add(MapFit_pterr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]);
		      }
		}
	    }

	  if ((rint(PtCut*100.) == rint(RangePt*100.)) || (H_ptres_pt_eta[0]->GetXaxis()->FindBin(rint((PtCut-PtStep)*100.)/100.) != H_ptres_pt_eta[0]->GetXaxis()->FindBin(rint(PtCut*100.)/100.)))
	    {
	      for (int i = 0; i < 2; i++)
		{
		  DeconvolvErrsq    = pterr[i]->GetMean() / 1000000.;
		  ErrDeconvolvErrsq = pterr[i]->GetRMS() / sqrt(pterr[i]->GetEntries()) / 1000000.;
		  
		  if ((FitAlgorithm (Gauss,pt[i]) == true) && (Gauss->GetParameter(2)*Gauss->GetParameter(2) > DeconvolvErrsq))
		    {
		      Val = sqrt(Gauss->GetParameter(2)*Gauss->GetParameter(2) - DeconvolvErrsq);
		      Err = 1./(2.*Val) * sqrt(powf(2.*Gauss->GetParameter(2)*Gauss->GetParError(2),2.) + powf(ErrDeconvolvErrsq,2.));
		      
		      H_ptres_pt_eta[i]->SetBinContent(H_ptres_pt_eta[i]->GetXaxis()->FindBin(rint((PtCut-PtStep)*100.)/100.), Val / ((PtCut + OldLowEdgeBin)/2.));
		      H_ptres_pt_eta[i]->SetBinError(H_ptres_pt_eta[i]->GetXaxis()->FindBin(rint((PtCut-PtStep)*100.)/100.), Err / ((PtCut + OldLowEdgeBin)/2.));
		      
		      H_skew_pt_eta[i]->Fill(rint((PtCut-PtStep)*100.)/100., pt[i]->GetSkewness());
		    }

		  pt[i]->Reset();
		  pterr[i]->Reset();
		}

	      OldLowEdgeBin = PtCut;
	    }
      	}
      
      for (double EtaCut = EtaStep; rint(EtaCut*100.) <= rint(RangeEta*100.); EtaCut += EtaStep)
	{
	  EtaCutStr.str("");
	  EtaCutStr << EtaCut;
 
	  for (double PtCut = MinPtTk+PtStep; rint(PtCut*100.) <= rint(RangePt*100.); PtCut += PtStep)
	    {
	      PtCutStr.str("");
	      PtCutStr << PtCut;
 
	      if ((MapFit_pt.find(PtCutStr.str().c_str()) != MapFit_pt.end()) &&
		  (MapFit_pt[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) != MapFit_pt[PtCutStr.str().c_str()].end()) &&
		  (MapFit_pterr.find(PtCutStr.str().c_str()) != MapFit_pterr.end()) &&
		  (MapFit_pterr[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) != MapFit_pterr[PtCutStr.str().c_str()].end()))
		{
		  for (int i = 0; i < 5; i++)
		    if ((rint(PtCut*100.) > rint(PtThr[i]*100.)) && (rint(PtCut*100.) <= rint(PtThr[i+1]*100.)))
		      {
			pt[i]->Add(MapFit_pt[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]);
			pterr[i]->Add(MapFit_pterr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]);
		      }
		}
	    }

	  for (int i = 0; i < 5; i++)
	    {
	      DeconvolvErrsq    = pterr[i]->GetMean() / 1000000.;
	      ErrDeconvolvErrsq = pterr[i]->GetRMS() / sqrt(pterr[i]->GetEntries()) / 1000000.;
	      	      
	      if ((FitAlgorithm (Gauss,pt[i]) == true) && (Gauss->GetParameter(2)*Gauss->GetParameter(2) > DeconvolvErrsq))
		{
		  Val = sqrt(Gauss->GetParameter(2)*Gauss->GetParameter(2) - DeconvolvErrsq);
		  Err = 1./(2.*Val) * sqrt(powf(2.*Gauss->GetParameter(2)*Gauss->GetParError(2),2.) + powf(ErrDeconvolvErrsq,2.));

		  H_ptres_eta_pt[i]->SetBinContent(H_ptres_eta_pt[i]->GetXaxis()->FindBin(rint((EtaCut-EtaStep)*100.)/100.), Val / ((PtThr[i+1] + PtThr[i]) / 2.));
		  H_ptres_eta_pt[i]->SetBinError(H_ptres_eta_pt[i]->GetXaxis()->FindBin(rint((EtaCut-EtaStep)*100.)/100.), Err / ((PtThr[i+1] + PtThr[i]) / 2.));
		}

	      pt[i]->Reset();
	      pterr[i]->Reset();
	    }
	}
      
      for (int i = 0; i < 5; i++) { delete pt[i]; delete pterr[i]; }
    }

  delete Gauss;
}


void MakePlots(string datatype)
{
  double SumNum;
  double SumDen;

  if (datatype.compare("Vx") == 0)
    {
      // #################
      // ### VERTEXING ###
      // #################

      // ########################
      // # PURITY & EFFICIENCTY #
      // ########################
      SumNum = H_vx_counters->GetBinContent(H_vx_counters->GetXaxis()->FindBin("MatchEff"));
      SumDen = H_vx_counters->GetBinContent(H_vx_counters->GetXaxis()->FindBin("Strip"));
      H_vx_pur_eff->Fill("Efficiency", (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
      H_vx_pur_eff->SetBinError(H_vx_pur_eff->GetXaxis()->FindBin("Efficiency"), (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);

      SumNum = H_vx_counters->GetBinContent(H_vx_counters->GetXaxis()->FindBin("MatchPur"));
      SumDen = H_vx_counters->GetBinContent(H_vx_counters->GetXaxis()->FindBin("Pixel"));
      H_vx_pur_eff->Fill("Purity", (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
      H_vx_pur_eff->SetBinError(H_vx_pur_eff->GetXaxis()->FindBin("Purity"), (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);

      // #######################################
      // # Vertex Resolutions and Pulls vs Trk #
      // #######################################
      stringstream String;
      double Val,Err;

      TF1* Gauss = new TF1("Gauss", "gaus", -zRange/2., zRange/2.);

      for (map<string,TH1F*>::iterator it = map_H_vx_x_diff.begin(); it != map_H_vx_x_diff.end(); it++)
      	{
      	  String.str("");
      	  String << atoi(it->first.c_str());
	  
      	  Gauss->SetRange(-zRange/2., zRange/2.);

      	  // #############################
      	  // # Compute simple difference #
      	  // #############################
      	  if (map_H_vx_x_diff.find(String.str().c_str()) != map_H_vx_x_diff.end() && (FitAlgorithm (Gauss, map_H_vx_x_diff[String.str().c_str()]) == true))
      	    {
      	      H_vx_x_diff_vs_trk->SetBinContent(H_vx_x_diff_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Gauss->GetParameter(2));
      	      H_vx_x_diff_vs_trk->SetBinError(H_vx_x_diff_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Gauss->GetParError(2));
      	    }
      	  if (map_H_vx_y_diff.find(String.str().c_str()) != map_H_vx_y_diff.end() && (FitAlgorithm (Gauss, map_H_vx_y_diff[String.str().c_str()]) == true))
      	    {
      	      H_vx_y_diff_vs_trk->SetBinContent(H_vx_y_diff_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Gauss->GetParameter(2));
      	      H_vx_y_diff_vs_trk->SetBinError(H_vx_y_diff_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Gauss->GetParError(2));
      	    }
      	  if (map_H_vx_z_diff.find(String.str().c_str()) != map_H_vx_z_diff.end() && (FitAlgorithm (Gauss, map_H_vx_z_diff[String.str().c_str()]) == true))
      	    {
      	      H_vx_z_diff_vs_trk->SetBinContent(H_vx_z_diff_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Gauss->GetParameter(2));
      	      H_vx_z_diff_vs_trk->SetBinError(H_vx_z_diff_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Gauss->GetParError(2));
      	    }

    	  // #################################################
    	  // # Compute difference deconvolving vx resolution #
    	  // #################################################
    	  if ((map_H_vx_x_err2.find(String.str().c_str()) != map_H_vx_x_err2.end()) &&
    	      (map_H_vx_x_diff.find(String.str().c_str()) != map_H_vx_x_diff.end()) && 
    	      (FitAlgorithm (Gauss, map_H_vx_x_diff[String.str().c_str()]) == true) &&
    	      (powf(Gauss->GetParameter(2),2.) > map_H_vx_x_err2[String.str().c_str()]->GetMean()))
    	    {
    	      Val = sqrt(powf(Gauss->GetParameter(2),2.) - map_H_vx_x_err2[String.str().c_str()]->GetMean());
	      
    	      Err = 1./(2.*Val) * sqrt(powf(2.*Gauss->GetParameter(2)*Gauss->GetParError(2),2.) +
    				       powf(map_H_vx_x_err2[String.str().c_str()]->GetRMS()/sqrt(map_H_vx_x_err2[String.str().c_str()]->GetEntries()),2.));
	      
    	      H_vx_x_diff_noerr_vs_trk->SetBinContent(H_vx_x_diff_noerr_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Val);
    	      H_vx_x_diff_noerr_vs_trk->SetBinError(H_vx_x_diff_noerr_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Err);
    	    }
    	  if ((map_H_vx_y_err2.find(String.str().c_str()) != map_H_vx_y_err2.end()) &&
    	      (map_H_vx_y_diff.find(String.str().c_str()) != map_H_vx_y_diff.end()) &&
    	      (FitAlgorithm (Gauss, map_H_vx_y_diff[String.str().c_str()]) == true) &&
    	      (powf(Gauss->GetParameter(2),2.) > map_H_vx_y_err2[String.str().c_str()]->GetMean()))
    	    {
    	      Val = sqrt(powf(Gauss->GetParameter(2),2.) - map_H_vx_y_err2[String.str().c_str()]->GetMean());
	      
    	      Err = 1./(2.*Val) * sqrt(powf(2.*Gauss->GetParameter(2)*Gauss->GetParError(2),2.) +
    				       powf(map_H_vx_y_err2[String.str().c_str()]->GetRMS()/sqrt(map_H_vx_y_err2[String.str().c_str()]->GetEntries()),2.));
	      
    	      H_vx_y_diff_noerr_vs_trk->SetBinContent(H_vx_y_diff_noerr_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Val);
    	      H_vx_y_diff_noerr_vs_trk->SetBinError(H_vx_y_diff_noerr_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Err);
    	    }
    	  if ((map_H_vx_z_err2.find(String.str().c_str()) != map_H_vx_z_err2.end()) &&
    	      (map_H_vx_z_diff.find(String.str().c_str()) != map_H_vx_z_diff.end()) &&
    	      (FitAlgorithm (Gauss, map_H_vx_z_diff[String.str().c_str()]) == true) &&
    	      (powf(Gauss->GetParameter(2),2.) > map_H_vx_z_err2[String.str().c_str()]->GetMean()))
    	    {
    	      Val = sqrt(powf(Gauss->GetParameter(2),2.) - map_H_vx_z_err2[String.str().c_str()]->GetMean());
	      
    	      Err = 1./(2.*Val) * sqrt(powf(2.*Gauss->GetParameter(2)*Gauss->GetParError(2),2.) +
    				       powf(map_H_vx_z_err2[String.str().c_str()]->GetRMS()/sqrt(map_H_vx_z_err2[String.str().c_str()]->GetEntries()),2.));
	      
    	      H_vx_z_diff_noerr_vs_trk->SetBinContent(H_vx_z_diff_noerr_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Val);
    	      H_vx_z_diff_noerr_vs_trk->SetBinError(H_vx_z_diff_noerr_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Err);
    	    }

    	  Gauss->SetRange(-pullsRange/2., pullsRange/2.);

    	  // #################
    	  // # Compute pulls #
    	  // #################
    	  if (map_H_vx_x_diff_pull.find(String.str().c_str()) != map_H_vx_x_diff_pull.end() && (FitAlgorithm (Gauss, map_H_vx_x_diff_pull[String.str().c_str()]) == true))
    	    {
    	      H_vx_x_diff_pull_vs_trk->SetBinContent(H_vx_x_diff_pull_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Gauss->GetParameter(2));
    	      H_vx_x_diff_pull_vs_trk->SetBinError(H_vx_x_diff_pull_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Gauss->GetParError(2));
    	    }
    	  if (map_H_vx_y_diff_pull.find(String.str().c_str()) != map_H_vx_y_diff_pull.end() && (FitAlgorithm (Gauss, map_H_vx_y_diff_pull[String.str().c_str()]) == true))
    	    {
    	      H_vx_y_diff_pull_vs_trk->SetBinContent(H_vx_y_diff_pull_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Gauss->GetParameter(2));
    	      H_vx_y_diff_pull_vs_trk->SetBinError(H_vx_y_diff_pull_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Gauss->GetParError(2));
    	    }
    	  if (map_H_vx_z_diff_pull.find(String.str().c_str()) != map_H_vx_z_diff_pull.end() && (FitAlgorithm (Gauss, map_H_vx_z_diff_pull[String.str().c_str()]) == true))
    	    {
    	      H_vx_z_diff_pull_vs_trk->SetBinContent(H_vx_z_diff_pull_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Gauss->GetParameter(2));
    	      H_vx_z_diff_pull_vs_trk->SetBinError(H_vx_z_diff_pull_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Gauss->GetParError(2));
    	    }
    	}

      Gauss->SetRange(-zRange/2., zRange/2.);
      
      for (map<string,TH1F*>::iterator it = map_H_vx_x_diff_orig.begin(); it != map_H_vx_x_diff_orig.end(); it++)
    	{
    	  String.str("");
    	  String << atoi(it->first.c_str());

    	  // ####################################################
    	  // # Compute simple difference with original vertices #
    	  // ####################################################
    	  if (map_H_vx_x_diff_orig.find(String.str().c_str()) != map_H_vx_x_diff_orig.end() && (FitAlgorithm (Gauss, map_H_vx_x_diff_orig[String.str().c_str()]) == true))
    	    {
    	      H_vx_x_diff_orig_vs_trk->SetBinContent(H_vx_x_diff_orig_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Gauss->GetParameter(2));
    	      H_vx_x_diff_orig_vs_trk->SetBinError(H_vx_x_diff_orig_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Gauss->GetParError(2));
    	    }
    	  if (map_H_vx_y_diff_orig.find(String.str().c_str()) != map_H_vx_y_diff_orig.end() && (FitAlgorithm (Gauss, map_H_vx_y_diff_orig[String.str().c_str()]) == true))
    	    {
    	      H_vx_y_diff_orig_vs_trk->SetBinContent(H_vx_y_diff_orig_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Gauss->GetParameter(2));
    	      H_vx_y_diff_orig_vs_trk->SetBinError(H_vx_y_diff_orig_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Gauss->GetParError(2));
    	    }
    	  if (map_H_vx_z_diff_orig.find(String.str().c_str()) != map_H_vx_z_diff_orig.end() && (FitAlgorithm (Gauss, map_H_vx_z_diff_orig[String.str().c_str()]) == true))
    	    {
    	      H_vx_z_diff_orig_vs_trk->SetBinContent(H_vx_z_diff_orig_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Gauss->GetParameter(2));
    	      H_vx_z_diff_orig_vs_trk->SetBinError(H_vx_z_diff_orig_vs_trk->GetXaxis()->FindBin(atoi(it->first.c_str())), Gauss->GetParError(2));
    	    }
    	}

      delete Gauss;
    }
  else if (datatype.compare("Tk") == 0)
    {
      // ################
      // ### TRACKING ###
      // ################

      double TotalSumNum;
      double TotalSumDen;
     
      // ################################
      // # 3-PIXEL-HIT EFFICIENCY VS PT #
      // ################################
      for (int i = 1; i <= H_trk_eff_pt->GetNbinsX(); i++)
      	{
      	  H_trk_eff_pt->SetBinContent(i, (H_trk_eff_pt_den->GetBinContent(i) > H_trk_eff_pt_num->GetBinContent(i)) ?
      				      H_trk_eff_pt_num->GetBinContent(i) / H_trk_eff_pt_den->GetBinContent(i) : 0.);
      	  H_trk_eff_pt->SetBinError(i, (H_trk_eff_pt_den->GetBinContent(i) > H_trk_eff_pt_num->GetBinContent(i)) ?
      				    sqrt((H_trk_eff_pt_num->GetBinContent(i)*(H_trk_eff_pt_den->GetBinContent(i) - H_trk_eff_pt_num->GetBinContent(i))) /
      					 powf(H_trk_eff_pt_den->GetBinContent(i),3.)) : 0.);
      	}

      // #################################
      // # 3-PIXEL-HIT EFFICIENCY VS ETA #
      // #################################
      for (int i = 1; i <= H_trk_eff_eta->GetNbinsX(); i++)
      	{
      	  H_trk_eff_eta->SetBinContent(i, (H_trk_eff_eta_den->GetBinContent(i) > H_trk_eff_eta_num->GetBinContent(i)) ?
      				       H_trk_eff_eta_num->GetBinContent(i) / H_trk_eff_eta_den->GetBinContent(i) : 0.);
      	  H_trk_eff_eta->SetBinError(i, (H_trk_eff_eta_den->GetBinContent(i) > H_trk_eff_eta_num->GetBinContent(i)) ?
      				     sqrt((H_trk_eff_eta_num->GetBinContent(i)*(H_trk_eff_eta_den->GetBinContent(i) - H_trk_eff_eta_num->GetBinContent(i))) /
      					  powf(H_trk_eff_eta_den->GetBinContent(i),3.)) : 0.);
      	}

      // #################################
      // # 3-PIXEL-HIT EFFICIENCY VS PHI #
      // #################################
      for (int i = 1; i <= H_trk_eff_phi[0]->GetNbinsX(); i++)
      	{
      	  for (int j = 0; j < 3; j++)
      	    {
      	      H_trk_eff_phi[j]->SetBinContent(i, (H_trk_eff_phi_den[j]->GetBinContent(i) > H_trk_eff_phi_num[j]->GetBinContent(i)) ?
      					      H_trk_eff_phi_num[j]->GetBinContent(i) / H_trk_eff_phi_den[j]->GetBinContent(i) : 0.);
      	      H_trk_eff_phi[j]->SetBinError(i, (H_trk_eff_phi_den[j]->GetBinContent(i) > H_trk_eff_phi_num[j]->GetBinContent(i)) ?
      					    sqrt((H_trk_eff_phi_num[j]->GetBinContent(i)*(H_trk_eff_phi_den[j]->GetBinContent(i) - H_trk_eff_phi_num[j]->GetBinContent(i))) /
      						 powf(H_trk_eff_phi_den[j]->GetBinContent(i),3.)) : 0.);
      	    }
      	}

      // ########################
      // # ALL EFFICIENCY VS PT #
      // ########################
      for (int i = 1; i <= H_trk_toteff_pt->GetNbinsX(); i++)
      	{
      	  H_trk_toteff_pt->SetBinContent(i, (H_trk_toteff_pt_den->GetBinContent(i) > H_trk_toteff_pt_num->GetBinContent(i)) ?
      				      H_trk_toteff_pt_num->GetBinContent(i) / H_trk_toteff_pt_den->GetBinContent(i) : 0.);
      	  H_trk_toteff_pt->SetBinError(i, (H_trk_toteff_pt_den->GetBinContent(i) > H_trk_toteff_pt_num->GetBinContent(i)) ?
      				    sqrt((H_trk_toteff_pt_num->GetBinContent(i)*(H_trk_toteff_pt_den->GetBinContent(i) - H_trk_toteff_pt_num->GetBinContent(i))) /
      					 powf(H_trk_toteff_pt_den->GetBinContent(i),3.)) : 0.);
      	}

      // #########################
      // # ALL EFFICIENCY VS ETA #
      // #########################
      for (int i = 1; i <= H_trk_toteff_eta->GetNbinsX(); i++)
      	{
      	  H_trk_toteff_eta->SetBinContent(i, (H_trk_toteff_eta_den->GetBinContent(i) > H_trk_toteff_eta_num->GetBinContent(i)) ?
      				       H_trk_toteff_eta_num->GetBinContent(i) / H_trk_toteff_eta_den->GetBinContent(i) : 0.);
      	  H_trk_toteff_eta->SetBinError(i, (H_trk_toteff_eta_den->GetBinContent(i) > H_trk_toteff_eta_num->GetBinContent(i)) ?
      				     sqrt((H_trk_toteff_eta_num->GetBinContent(i)*(H_trk_toteff_eta_den->GetBinContent(i) - H_trk_toteff_eta_num->GetBinContent(i))) /
      					  powf(H_trk_toteff_eta_den->GetBinContent(i),3.)) : 0.);
      	}

      // #########################
      // # ALL EFFICIENCY VS PHI #
      // #########################
      for (int i = 1; i <= H_trk_toteff_phi[0]->GetNbinsX(); i++)
      	{
      	  for (int j = 0; j < 3; j++)
      	    {
      	      H_trk_toteff_phi[j]->SetBinContent(i, (H_trk_toteff_phi_den[j]->GetBinContent(i) > H_trk_toteff_phi_num[j]->GetBinContent(i)) ?
      					      H_trk_toteff_phi_num[j]->GetBinContent(i) / H_trk_toteff_phi_den[j]->GetBinContent(i) : 0.);
      	      H_trk_toteff_phi[j]->SetBinError(i, (H_trk_toteff_phi_den[j]->GetBinContent(i) > H_trk_toteff_phi_num[j]->GetBinContent(i)) ?
      					    sqrt((H_trk_toteff_phi_num[j]->GetBinContent(i)*(H_trk_toteff_phi_den[j]->GetBinContent(i) - H_trk_toteff_phi_num[j]->GetBinContent(i))) /
      						 powf(H_trk_toteff_phi_den[j]->GetBinContent(i),3.)) : 0.);
      	    }
      	}

      // ################################
      // # PIXEL-TRACKS DOUBLE COUNTING #
      // ################################
      for (int i = 1; i <= H_pix_trk_doubleCounter_eta->GetNbinsX(); i++)
	{
	  H_pix_trk_doubleCounter_eta_rel->SetBinContent(i, (H_pix_trk_eta->GetBinContent(i) > H_pix_trk_doubleCounter_eta->GetBinContent(i)) ?
							 H_pix_trk_doubleCounter_eta->GetBinContent(i) / H_pix_trk_eta->GetBinContent(i) : 0.);
	  H_pix_trk_doubleCounter_eta_rel->SetBinError(i, (H_pix_trk_eta->GetBinContent(i) > H_pix_trk_doubleCounter_eta->GetBinContent(i)) ?
						       sqrt((H_pix_trk_doubleCounter_eta->GetBinContent(i)*(H_pix_trk_eta->GetBinContent(i) - H_pix_trk_doubleCounter_eta->GetBinContent(i))) /
							    powf(H_pix_trk_eta->GetBinContent(i),3.)) : 0.);
	}
      for (int i = 1; i <= H_pix_trk_doubleCounter_pt->GetNbinsX(); i++)
	{
	  H_pix_trk_doubleCounter_pt_rel->SetBinContent(i, (H_pix_trk_pt->GetBinContent(i) > H_pix_trk_doubleCounter_pt->GetBinContent(i)) ?
							H_pix_trk_doubleCounter_pt->GetBinContent(i) / H_pix_trk_pt->GetBinContent(i) : 0.);
	  H_pix_trk_doubleCounter_pt_rel->SetBinError(i, (H_pix_trk_pt->GetBinContent(i) > H_pix_trk_doubleCounter_pt->GetBinContent(i)) ?
						      sqrt((H_pix_trk_doubleCounter_pt->GetBinContent(i)*(H_pix_trk_pt->GetBinContent(i) - H_pix_trk_doubleCounter_pt->GetBinContent(i))) /
							   powf(H_pix_trk_pt->GetBinContent(i),3.)) : 0.);
	}

      // ##################################
      // # TRACKER-TRACKS DOUBLE COUNTING #
      // ##################################
      for (int i = 1; i <= H_trk_trk_doubleCounter_eta->GetNbinsX(); i++)
      	{
      	  H_trk_trk_doubleCounter_eta_rel->SetBinContent(i, (H_trk_trk_eta->GetBinContent(i) > H_trk_trk_doubleCounter_eta->GetBinContent(i)) ?
      							 H_trk_trk_doubleCounter_eta->GetBinContent(i) / H_trk_trk_eta->GetBinContent(i) : 0.);
      	  H_trk_trk_doubleCounter_eta_rel->SetBinError(i, (H_trk_trk_eta->GetBinContent(i) > H_trk_trk_doubleCounter_eta->GetBinContent(i)) ?
      						       sqrt((H_trk_trk_doubleCounter_eta->GetBinContent(i)*(H_trk_trk_eta->GetBinContent(i) - H_trk_trk_doubleCounter_eta->GetBinContent(i))) /
      							    powf(H_trk_trk_eta->GetBinContent(i),3.)) : 0.);
      	}
      for (int i = 1; i <= H_trk_trk_doubleCounter_pt->GetNbinsX(); i++)
      	{
      	  H_trk_trk_doubleCounter_pt_rel->SetBinContent(i, (H_trk_trk_pt->GetBinContent(i) > H_trk_trk_doubleCounter_pt->GetBinContent(i)) ?
      							H_trk_trk_doubleCounter_pt->GetBinContent(i) / H_trk_trk_pt->GetBinContent(i) : 0.);
      	  H_trk_trk_doubleCounter_pt_rel->SetBinError(i, (H_trk_trk_pt->GetBinContent(i) > H_trk_trk_doubleCounter_pt->GetBinContent(i)) ?
      						      sqrt((H_trk_trk_doubleCounter_pt->GetBinContent(i)*(H_trk_trk_pt->GetBinContent(i) - H_trk_trk_doubleCounter_pt->GetBinContent(i))) /
      							   powf(H_trk_trk_pt->GetBinContent(i),3.)) : 0.);
      	}

      // #####################
      // # PURITY vs CHI2-Pt #
      // #####################
      for (int k = 0; k < 2; k++)
  	for (int i = 1; i <= H_purity_trked_normchi2_pt_eta[k]->GetNbinsX(); i++)    
  	  for (int j = 1; j <= H_purity_trked_normchi2_pt_eta[k]->GetNbinsY(); j++)
  	    if (H_pix_trk_normchi2_pt_eta[k]->GetBinContent(i,j) >= MinEntriesEffPur)
  	      H_purity_trked_normchi2_pt_eta[k]->SetBinContent(i,j, H_trked_normchi2_pt_pur_eta[k]->GetBinContent(i,j) / H_pix_trk_normchi2_pt_eta[k]->GetBinContent(i,j));
      
      // #################
      // # PURITY vs ETA #
      // #################
      TotalSumNum = 0.;
      TotalSumDen = 0.;
      for (int i = 1; i <= H_pix_trk_eta_phi->GetNbinsX(); i++)    
  	{
  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int j = 1; j <= H_pix_trk_eta_phi->GetNbinsY(); j++)
  	    {
  	      SumNum = SumNum + H_trked_eta_phi_pur->GetBinContent(i,j);
  	      SumDen = SumDen + H_pix_trk_eta_phi->GetBinContent(i,j);
  	    }

  	  TotalSumNum = TotalSumNum + SumNum;
  	  TotalSumDen = TotalSumDen + SumDen;

  	  H_purity_trked_eta->SetBinContent(i, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_purity_trked_eta->SetBinError(i, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);
  	}
      H_tk_pur_eff->Fill("Purity", (TotalSumDen >= MinEntriesEffPur) ? TotalSumNum / TotalSumDen : 0.);
      H_tk_pur_eff->SetBinError(H_tk_pur_eff->GetXaxis()->FindBin("Purity"), (TotalSumDen >= MinEntriesEffPur) ? sqrt((TotalSumNum*(TotalSumDen-TotalSumNum))/powf(TotalSumDen,3.)) : 0.);

      // ################
      // # PURITY vs Pt #
      // ################
      for (int j = H_pix_trk_eta_pt->GetYaxis()->FindBin(rint(MinPtTk*100.)/100.); j <= H_pix_trk_eta_pt->GetNbinsY(); j++)    
  	{
  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int i = H_pix_trk_eta_pt->GetXaxis()->FindBin(rint(-EtaThr*100.)/100.); i <= H_pix_trk_eta_pt->GetXaxis()->FindBin(rint(EtaThr*100.)/100.); i++)
  	    {
  	      SumNum = SumNum + H_trked_eta_pt_pur->GetBinContent(i,j);
  	      SumDen = SumDen + H_pix_trk_eta_pt->GetBinContent(i,j);
  	    }
  	  H_purity_trked_pt_eta[0]->SetBinContent(j, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_purity_trked_pt_eta[0]->SetBinError(j, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);
	  
  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int i = 1; i < H_pix_trk_eta_pt->GetXaxis()->FindBin(rint(-EtaThr*100.)/100.); i++)
  	    {
  	      SumNum = SumNum + H_trked_eta_pt_pur->GetBinContent(i,j);
  	      SumDen = SumDen + H_pix_trk_eta_pt->GetBinContent(i,j);
  	    }
  	  for (int i = H_pix_trk_eta_pt->GetXaxis()->FindBin(rint(EtaThr*100.)/100.) + 1; i <= H_pix_trk_eta_pt->GetNbinsX(); i++)
  	    {
  	      SumNum = SumNum + H_trked_eta_pt_pur->GetBinContent(i,j);
  	      SumDen = SumDen + H_pix_trk_eta_pt->GetBinContent(i,j);
  	    }
  	  H_purity_trked_pt_eta[1]->SetBinContent(j, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_purity_trked_pt_eta[1]->SetBinError(j, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);
  	}

      for (int j = H_pix_trk_eta_pt->GetYaxis()->FindBin(rint(MinPtTk*100.)/100.); j <= H_pix_trk_eta_pt->GetNbinsY(); j++)    
  	{
  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int i = H_pix_trk_eta_pt->GetXaxis()->FindBin(rint(-MaxEtaTkTrk*100.)/100.); i <= H_pix_trk_eta_pt->GetXaxis()->FindBin(rint(MaxEtaTkTrk*100.)/100.); i++)
  	    {
  	      SumNum = SumNum + H_trked_eta_pt_pur->GetBinContent(i,j);
  	      SumDen = SumDen + H_pix_trk_eta_pt->GetBinContent(i,j);
  	    }
  	  H_purity_trked_pt_whole_eta->SetBinContent(j, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_purity_trked_pt_whole_eta->SetBinError(j, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);
  	}
  
      // #################
      // # PURITY vs PHI #
      // #################
      for (int i = 1; i <= H_pix_trk_eta_phi->GetNbinsY(); i++)    
  	{
  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int j = H_pix_trk_eta_phi->GetXaxis()->FindBin(rint(-MaxEtaTkTrk*100.)/100.); j <= H_pix_trk_eta_phi->GetXaxis()->FindBin(rint(-EtaThr*100.)/100.); j++)
  	    {
  	      SumNum = SumNum + H_trked_eta_phi_pur->GetBinContent(j,i);
  	      SumDen = SumDen + H_pix_trk_eta_phi->GetBinContent(j,i);
  	    }
  	  H_purity_trked_phi[0]->SetBinContent(i, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_purity_trked_phi[0]->SetBinError(i, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);

  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int j = H_pix_trk_eta_phi->GetXaxis()->FindBin(rint(-EtaThr*100.)/100.); j <= H_pix_trk_eta_phi->GetXaxis()->FindBin(rint(EtaThr*100.)/100.); j++)
  	    {
  	      SumNum = SumNum + H_trked_eta_phi_pur->GetBinContent(j,i);
  	      SumDen = SumDen + H_pix_trk_eta_phi->GetBinContent(j,i);
  	    }
  	  H_purity_trked_phi[1]->SetBinContent(i, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_purity_trked_phi[1]->SetBinError(i, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);

  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int j = H_pix_trk_eta_phi->GetXaxis()->FindBin(rint(EtaThr*100.)/100.); j <= H_pix_trk_eta_phi->GetXaxis()->FindBin(rint(MaxEtaTkTrk*100.)/100.); j++)
  	    {
  	      SumNum = SumNum + H_trked_eta_phi_pur->GetBinContent(j,i);
  	      SumDen = SumDen + H_pix_trk_eta_phi->GetBinContent(j,i);
  	    }
  	  H_purity_trked_phi[2]->SetBinContent(i, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_purity_trked_phi[2]->SetBinError(i, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);
  	}

      // #########################
      // # EFFICIENCY vs CHI2-Pt #
      // #########################
      for (int k = 0; k < 2; k++)
  	for (int i = 1; i <= H_efficiency_trked_normchi2_pt_eta[k]->GetNbinsX(); i++)    
  	  for (int j = 1; j <= H_efficiency_trked_normchi2_pt_eta[k]->GetNbinsY(); j++)
  	    if (H_trk_trk_normchi2_pt_eta[k]->GetBinContent(i,j) >= MinEntriesEffPur)
  	      H_efficiency_trked_normchi2_pt_eta[k]->SetBinContent(i,j, H_trked_normchi2_pt_eff_eta[k]->GetBinContent(i,j) / H_trk_trk_normchi2_pt_eta[k]->GetBinContent(i,j));

      // #####################
      // # EFFICIENCY vs ETA #
      // #####################
      TotalSumNum = 0.;
      TotalSumDen = 0.;
      for (int i = 1; i <= H_trk_trk_eta_phi->GetNbinsX(); i++)    
  	{
  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int j = 1; j <= H_trk_trk_eta_phi->GetNbinsY(); j++)
  	    {
  	      SumNum = SumNum + H_trked_eta_phi_eff->GetBinContent(i,j);
  	      SumDen = SumDen + H_trk_trk_eta_phi->GetBinContent(i,j);
  	    }

  	  TotalSumNum = TotalSumNum + SumNum;
  	  TotalSumDen = TotalSumDen + SumDen;

  	  H_efficiency_trked_eta->SetBinContent(i, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_efficiency_trked_eta->SetBinError(i, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);
  	}
      H_tk_pur_eff->Fill("Efficiency", (TotalSumDen >= MinEntriesEffPur) ? TotalSumNum / TotalSumDen : 0.);      
      H_tk_pur_eff->SetBinError(H_tk_pur_eff->GetXaxis()->FindBin("Efficiency"), (TotalSumDen >= MinEntriesEffPur) ? sqrt((TotalSumNum*(TotalSumDen-TotalSumNum))/powf(TotalSumDen,3.)) : 0.);

      // ##########################
      // # ALGO-EFFICIENCY vs ETA #
      // ##########################
      TotalSumNum = 0.;
      TotalSumDen = 0.;
      for (int i = 1; i <= H_trk_trk_eta_phi_algo->GetNbinsX(); i++)    
  	{
  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int j = 1; j <= H_trk_trk_eta_phi_algo->GetNbinsY(); j++)
  	    {
  	      SumNum = SumNum + H_trked_eta_phi_eff->GetBinContent(i,j);
  	      SumDen = SumDen + H_trk_trk_eta_phi_algo->GetBinContent(i,j);
  	    }

  	  TotalSumNum = TotalSumNum + SumNum;
  	  TotalSumDen = TotalSumDen + SumDen;

  	  H_efficiency_trked_eta_algo->SetBinContent(i, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_efficiency_trked_eta_algo->SetBinError(i, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);
  	}
      H_tk_pur_eff->Fill("Algo-Efficiency", (TotalSumDen >= MinEntriesEffPur) ? TotalSumNum / TotalSumDen : 0.);
      H_tk_pur_eff->SetBinError(H_tk_pur_eff->GetXaxis()->FindBin("Algo-Efficiency"), (TotalSumDen >= MinEntriesEffPur) ? sqrt((TotalSumNum*(TotalSumDen-TotalSumNum))/powf(TotalSumDen,3.)) : 0.);

      // ####################
      // # EFFICIENCY vs Pt #
      // ####################
      for (int j = H_trk_trk_eta_pt->GetYaxis()->FindBin(rint(MinPtTk*100.)/100.); j <= H_trk_trk_eta_pt->GetNbinsY(); j++)    
  	{
  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int i = H_trk_trk_eta_pt->GetXaxis()->FindBin(rint(-EtaThr*100.)/100.); i <= H_trk_trk_eta_pt->GetXaxis()->FindBin(rint(EtaThr*100.)/100.); i++)
  	    {
  	      SumNum = SumNum + H_trked_eta_pt_eff->GetBinContent(i,j);
  	      SumDen = SumDen + H_trk_trk_eta_pt->GetBinContent(i,j);
  	    }
  	  H_efficiency_trked_pt_eta[0]->SetBinContent(j, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_efficiency_trked_pt_eta[0]->SetBinError(j, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);

  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int i = 1; i < H_trk_trk_eta_pt->GetXaxis()->FindBin(rint(-EtaThr*100.)/100.); i++)
  	    {
  	      SumNum = SumNum + H_trked_eta_pt_eff->GetBinContent(i,j);
  	      SumDen = SumDen + H_trk_trk_eta_pt->GetBinContent(i,j);
  	    }
  	  for (int i = H_trk_trk_eta_pt->GetXaxis()->FindBin(rint(EtaThr*100.)/100.) + 1; i <= H_trk_trk_eta_pt->GetNbinsX(); i++)
  	    {
  	      SumNum = SumNum + H_trked_eta_pt_eff->GetBinContent(i,j);
  	      SumDen = SumDen + H_trk_trk_eta_pt->GetBinContent(i,j);
  	    }
  	  H_efficiency_trked_pt_eta[1]->SetBinContent(j, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_efficiency_trked_pt_eta[1]->SetBinError(j, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);
  	}

      for (int j = H_trk_trk_eta_pt->GetYaxis()->FindBin(rint(MinPtTk*100.)/100.); j <= H_trk_trk_eta_pt->GetNbinsY(); j++)    
  	{
  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int i = H_trk_trk_eta_pt->GetXaxis()->FindBin(rint(-MaxEtaTkTrk*100.)/100.); i <= H_trk_trk_eta_pt->GetXaxis()->FindBin(rint(MaxEtaTkTrk*100.)/100.); i++)
  	    {
  	      SumNum = SumNum + H_trked_eta_pt_eff->GetBinContent(i,j);
  	      SumDen = SumDen + H_trk_trk_eta_pt->GetBinContent(i,j);
  	    }
  	  H_efficiency_trked_pt_whole_eta->SetBinContent(j, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_efficiency_trked_pt_whole_eta->SetBinError(j, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);
  	}

      // #########################
      // # ALGO-EFFICIENCY vs Pt #
      // #########################
      for (int j = H_trk_trk_eta_pt_algo->GetYaxis()->FindBin(rint(MinPtTk*100.)/100.); j <= H_trk_trk_eta_pt_algo->GetNbinsY(); j++)    
  	{
  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int i = H_trk_trk_eta_pt_algo->GetXaxis()->FindBin(rint(-EtaThr*100.)/100.); i <= H_trk_trk_eta_pt_algo->GetXaxis()->FindBin(rint(EtaThr*100.)/100.); i++)
  	    {
  	      SumNum = SumNum + H_trked_eta_pt_eff->GetBinContent(i,j);
  	      SumDen = SumDen + H_trk_trk_eta_pt_algo->GetBinContent(i,j);
  	    }
  	  H_efficiency_trked_pt_eta_algo[0]->SetBinContent(j, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_efficiency_trked_pt_eta_algo[0]->SetBinError(j, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);

  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int i = 1; i < H_trk_trk_eta_pt_algo->GetXaxis()->FindBin(rint(-EtaThr*100.)/100.); i++)
  	    {
  	      SumNum = SumNum + H_trked_eta_pt_eff->GetBinContent(i,j);
  	      SumDen = SumDen + H_trk_trk_eta_pt_algo->GetBinContent(i,j);
  	    }
  	  for (int i = H_trk_trk_eta_pt_algo->GetXaxis()->FindBin(rint(EtaThr*100.)/100.) + 1; i <= H_trk_trk_eta_pt_algo->GetNbinsX(); i++)
  	    {
  	      SumNum = SumNum + H_trked_eta_pt_eff->GetBinContent(i,j);
  	      SumDen = SumDen + H_trk_trk_eta_pt_algo->GetBinContent(i,j);
  	    }
  	  H_efficiency_trked_pt_eta_algo[1]->SetBinContent(j, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_efficiency_trked_pt_eta_algo[1]->SetBinError(j, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);
  	}

      for (int j = H_trk_trk_eta_pt_algo->GetYaxis()->FindBin(rint(MinPtTk*100.)/100.); j <= H_trk_trk_eta_pt_algo->GetNbinsY(); j++)    
  	{
  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int i = H_trk_trk_eta_pt_algo->GetXaxis()->FindBin(rint(-MaxEtaTkTrk*100.)/100.); i <= H_trk_trk_eta_pt_algo->GetXaxis()->FindBin(rint(MaxEtaTkTrk*100.)/100.); i++)
  	    {
  	      SumNum = SumNum + H_trked_eta_pt_eff->GetBinContent(i,j);
  	      SumDen = SumDen + H_trk_trk_eta_pt_algo->GetBinContent(i,j);
  	    }
  	  H_efficiency_trked_pt_whole_eta_algo->SetBinContent(j, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_efficiency_trked_pt_whole_eta_algo->SetBinError(j, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);
  	}

      // #####################
      // # EFFICIENCY vs PHI #
      // #####################
      for (int i = 1; i <= H_trk_trk_eta_phi->GetNbinsY(); i++)    
  	{
  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int j = H_trk_trk_eta_phi->GetXaxis()->FindBin(rint(-MaxEtaTkTrk*100.)/100.); j <= H_trk_trk_eta_phi->GetXaxis()->FindBin(rint(-EtaThr*100.)/100.); j++)
  	    {
  	      SumNum = SumNum + H_trked_eta_phi_eff->GetBinContent(j,i);
  	      SumDen = SumDen + H_trk_trk_eta_phi->GetBinContent(j,i);
  	    }
  	  H_efficiency_trked_phi[0]->SetBinContent(i, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_efficiency_trked_phi[0]->SetBinError(i, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);

  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int j = H_trk_trk_eta_phi->GetXaxis()->FindBin(rint(-EtaThr*100.)/100.); j <= H_trk_trk_eta_phi->GetXaxis()->FindBin(rint(EtaThr*100.)/100.); j++)
  	    {
  	      SumNum = SumNum + H_trked_eta_phi_eff->GetBinContent(j,i);
  	      SumDen = SumDen + H_trk_trk_eta_phi->GetBinContent(j,i);
  	    }
  	  H_efficiency_trked_phi[1]->SetBinContent(i, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_efficiency_trked_phi[1]->SetBinError(i, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);

  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int j = H_trk_trk_eta_phi->GetXaxis()->FindBin(rint(EtaThr*100.)/100.); j <= H_trk_trk_eta_phi->GetXaxis()->FindBin(rint(MaxEtaTkTrk*100.)/100.); j++)
  	    {
  	      SumNum = SumNum + H_trked_eta_phi_eff->GetBinContent(j,i);
  	      SumDen = SumDen + H_trk_trk_eta_phi->GetBinContent(j,i);
  	    }
  	  H_efficiency_trked_phi[2]->SetBinContent(i, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_efficiency_trked_phi[2]->SetBinError(i, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);
  	}

      // ##########################
      // # ALGO-EFFICIENCY vs PHI #
      // ##########################
      for (int i = 1; i <= H_trk_trk_eta_phi_algo->GetNbinsY(); i++)
  	{
  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int j = H_trk_trk_eta_phi_algo->GetXaxis()->FindBin(rint(-MaxEtaTkTrk*100.)/100.); j <= H_trk_trk_eta_phi_algo->GetXaxis()->FindBin(rint(-EtaThr*100.)/100.); j++)
  	    {
  	      SumNum = SumNum + H_trked_eta_phi_eff->GetBinContent(j,i);
  	      SumDen = SumDen + H_trk_trk_eta_phi_algo->GetBinContent(j,i);
  	    }
  	  H_efficiency_trked_phi_algo[0]->SetBinContent(i, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_efficiency_trked_phi_algo[0]->SetBinError(i, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);

  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int j = H_trk_trk_eta_phi_algo->GetXaxis()->FindBin(rint(-EtaThr*100.)/100.); j <= H_trk_trk_eta_phi_algo->GetXaxis()->FindBin(rint(EtaThr*100.)/100.); j++)
  	    {
  	      SumNum = SumNum + H_trked_eta_phi_eff->GetBinContent(j,i);
  	      SumDen = SumDen + H_trk_trk_eta_phi_algo->GetBinContent(j,i);
  	    }
  	  H_efficiency_trked_phi_algo[1]->SetBinContent(i, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_efficiency_trked_phi_algo[1]->SetBinError(i, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);

  	  SumNum = 0.;
  	  SumDen = 0.;
  	  for (int j = H_trk_trk_eta_phi_algo->GetXaxis()->FindBin(rint(EtaThr*100.)/100.); j <= H_trk_trk_eta_phi_algo->GetXaxis()->FindBin(rint(MaxEtaTkTrk*100.)/100.); j++)
  	    {
  	      SumNum = SumNum + H_trked_eta_phi_eff->GetBinContent(j,i);
  	      SumDen = SumDen + H_trk_trk_eta_phi_algo->GetBinContent(j,i);
  	    }
  	  H_efficiency_trked_phi_algo[2]->SetBinContent(i, (SumDen >= MinEntriesEffPur) ? SumNum / SumDen : 0.);
  	  H_efficiency_trked_phi_algo[2]->SetBinError(i, (SumDen >= MinEntriesEffPur) ? sqrt((SumNum*(SumDen-SumNum))/powf(SumDen,3.)) : 0.);
  	}

      FitAndFillMaps ("d0pulls", -pullsRange/2., pullsRange/2.);
      FitAndFillMaps ("d0", -d0Range/2., d0Range/2.);
      FitAndFillMaps ("dzpulls", -pullsRange/2., pullsRange/2.);
      FitAndFillMaps ("dz", -dzRange/2., dzRange/2.);
      FitAndFillMaps ("ptpulls", -pullsRange/2., pullsRange/2.);
      FitAndFillMaps ("pt", -ptRange/2., ptRange/2.);
    }
}


int main (int argc, char** argv)
{
  stringstream DataType;
  stringstream FileName;

  if (argc == 3)
    {
      FileName << argv[1];
      TFile* File = TFile::Open(FileName.str().c_str(),"UPDATE");

      DataType << argv[2];

      InitVariables();
      if (LoadHistos(DataType.str(), File) == true)
	{
	  MakePlots(DataType.str());
	  WriteHistos(DataType.str(), File);
	  File->Close();

	  return 1;
	}
      else
	{
	  cout << "Error while loading histograms from: " << FileName.str().c_str() << endl;
	  return 0;
	}
    }
  else
    {
      cout << "@@@ Synopsis @@@" << endl;
      cout << "./MyPixelAnalysisRatios filename.root datatype[Tk or Vx]" << endl;
      return 0;
    }
}

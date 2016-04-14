// -*- C++ -*-
//
// Package:    MyPixAnalyzer
// Class:      MyPixAnalyzer
// 
/**\class MyPixAnalyzer MyPixAnalyzer.h MyPixAnalysis/MyPixAnalyzer/interface/MyPixAnalyzer.h

   Description: <one line class summary>
   Study of the pixel stand alone tracking and vertexing

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Mauro Dinardo
//         Created:  Thu Dec  3 10:33:18 CET 2009
// $Id: MyPixAnalyzer.h,v 1.6 2010/08/28 15:47:10 dinardo Exp $


// System include files
#include <memory>

// User include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Provenance/interface/Provenance.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMask.h"
#include "CondFormats/L1TObjects/interface/L1GtPrescaleFactors.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMaskAlgoTrigRcd.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMaskTechTrigRcd.h"
#include "CondFormats/DataRecord/interface/L1GtPrescaleFactorsAlgoTrigRcd.h"
#include "CondFormats/DataRecord/interface/L1GtPrescaleFactorsTechTrigRcd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "SimTracker/VertexAssociation/interface/VertexAssociatorBase.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackAssociation/interface/ParametersDefinerForTP.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducerAlgorithm.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TROOT.h>
#include <TRandom.h>


#define XOR(A,B) ((A && !B) || (!A && B))
#define PI 3.141592653589793


using namespace std;
using namespace reco;


class MyPixAnalyzer : public edm::EDAnalyzer {
public:
  explicit MyPixAnalyzer (const edm::ParameterSet&);
  ~MyPixAnalyzer ();


private:
  // #########################
  // # STRUCT for TRACK DATA #
  // #########################
  struct TrackParType
  {
    unsigned int PxTrackValidHits;
    unsigned int TkTrackValidHits;
    unsigned int PxTrackLostHits;
    unsigned int TkTrackLostHits;
    double PxTrackChi2DoF;
    double TkTrackChi2DoF;
    double PxTrackPt;
    double PxTrackPtErr;
    double TkTrackPt;
    double TkTrackPtErr;
    bool BelongVx;
    double PxTrackDz;
    double PxTrackDzErr;
    double TkTrackDz;
    double TkTrackDzErr;
    double TkOriginZ; // With respec to to Beam Spot or parent vertex
    double PxTrackD0;
    double PxTrackD0Err;
    double TkTrackD0;
    double TkTrackD0Err;
    double TkOriginXY; // With respec to to Beam Spot or parent vertex
    double PxTrackTheta;
    double PxTrackThetaErr;
    double TkTrackTheta;
    double TkTrackThetaErr;
    double PxTrackPhi;
    double PxTrackPhiErr;
    double TkTrackPhi;
    double TkTrackPhiErr;
    double PxTrackEta;
    double PxTrackEtaErr;
    double TkTrackEta;
    double TkTrackEtaErr;
    bool AlgoEff;
    edm::View<Track>::const_iterator PxTkPtr;
  };


  // ##########################
  // # STRUCT for VERTEX DATA #
  // ##########################
  struct VertexParType
  {
    unsigned int NTracksMatched;
    unsigned int PxVxNTracks;
    unsigned int PxVxNTracksOrig;
    unsigned int TkVxNTracks;
    unsigned int TkVxNTracksOrig;
    double PxVxDoF;
    double PxVxDoFOrig;
    double TkVxDoF;
    double TkVxDoFOrig;
    double PxVxChi2DoF;
    double TkVxChi2DoF;
    double PxVxX;     // With respecto to Beam Spot
    double TkVxX;     // With respecto to Beam Spot
    double TkVxXOrig; // With respecto to Beam Spot
    double PxVxXerr;
    double TkVxXerr;
    double TkVxXerrOrig;
    double PxVxY;     // With respecto to Beam Spot
    double TkVxY;     // With respecto to Beam Spot
    double TkVxYOrig; // With respecto to Beam Spot
    double PxVxYerr;
    double TkVxYerr;
    double TkVxYerrOrig;
    double PxVxZ;     // With respecto to Beam Spot
    double TkVxZ;     // With respecto to Beam Spot
    double TkVxZOrig; // With respecto to Beam Spot
    double PxVxZerr;
    double TkVxZerr;
    double TkVxZerrOrig;
    const Vertex* PxVx;
    const TransientVertex* PxVxReFit;
    const Vertex* TkVxReco;
    TrackingVertexRef TkVxTP;
    const TransientVertex* TkVxReFit;
    vector<int>* MatchedTkTracksPix;
    vector<int>* MatchedTkTracksTrk;
    bool IsTrkPart;
  };


  virtual unsigned int CountSimChargedTracks(TrackingVertexRef VxTP);


  // ######################
  // # METHODS for TRACKs #
  // ######################
  virtual bool FindNpixelHits(unsigned int nHits,
			      const Track* TkTrack);
  virtual bool FindNpixelHits(unsigned int nHits,
			      TrackingParticleRef SimTrk);
  virtual void AddDataMap (string mapname,
			   double pt,
			   double eta,
			   double data,
			   int bins,
			   double beginval,
			   double endval);
  virtual void FillTrackHistos (string DataType,
				TrackParType* TrackPar);
  virtual void PixelTracksMC (const edm::Event& iEvent,
			      const edm::EventSetup& iSetup);
  virtual void PixelTracksData (const edm::Event& iEvent,
				const edm::EventSetup& iSetup);


  // ########################
  // # METHODS for VERTICES #
  // ########################
  virtual bool FillVertexHistos (string DataType,
				 VertexParType* VertexPar,
				 const edm::Event& iEvent,
				 const edm::EventSetup& iSetup);
  virtual void PixelVerticesMC (const edm::Event& iEvent,
				const edm::EventSetup& iSetup);
  virtual void PixelVerticesData (const edm::Event& iEvent,
				  const edm::EventSetup& iSetup);


  // ###############################
  // # METHODS for EVENT SELECTION #
  // ###############################
  virtual bool L1Analyzer (const edm::Event& iEvent,
			   const edm::EventSetup& iSetup,
			   DecisionWord MyAlgoMask,
			   DecisionWord MyTrigMask,
			   const bool PrintMsg,
                           const bool PrintWords);
  virtual bool HLTAnalyzer (const edm::Event& iEvent,
			    vector<string> MyHLTMask,
			    const bool PrintMsg);
  virtual bool EventSelection (const edm::Event& iEvent,
			       const edm::EventSetup& iSetup);


  // ##################
  // # GLOBAL METHODS #
  // ##################
  virtual void analyze (const edm::Event&,
			const edm::EventSetup&);
  virtual void beginJob ();
  virtual void endJob ();


  // @@@@@@ HISTOGRAMS @@@@@@
  TH1F* H_counters;


  // ####################
  // # Histogram TRACKS #
  // ####################

  // @@@@@@ Pixel tracks @@@@@@
  TH1F* H_pix_trk_normchi2[3];
  TH1F* H_pix_trk_pt[3];
  TH1F* H_pix_trk_dz[3];
  TH1F* H_pix_trk_d0[3];
  TH1F* H_pix_trk_eta[3];
  TH1F* H_pix_trk_theta[3];
  TH1F* H_pix_trk_phi[3];
  // @@@ Legenda @@@@
  // 0 = "NumE/AE"
  // 1 = "NumP"
  // 2 = "DenP"

  TH1F* H_pix_trk_doubleCounter;
  TH1F* H_pix_trk_doubleCounter_eta;
  TH1F* H_pix_trk_doubleCounter_pt;

  TH2F* H_pix_trk_normchi2_pt_eta[2];
  TH2F* H_pix_trk_eta_phi;
  TH2F* H_pix_trk_eta_pt;

  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // @ Histograms containing ratios @
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  TH1F* H_pix_trk_doubleCounter_pt_rel;
  TH1F* H_pix_trk_doubleCounter_eta_rel;

  // @@@@@@ Tracker tracks @@@@@@
  TH1F* H_trk_trk_normchi2[4];
  TH1F* H_trk_trk_pt[4];
  TH1F* H_trk_trk_dz[4];
  TH1F* H_trk_trk_d0[4];
  TH1F* H_trk_trk_eta[4];
  TH1F* H_trk_trk_theta[4];
  TH1F* H_trk_trk_phi[4];
  // @@@ Legenda @@@@
  // 0 = "NumE/AE"
  // 1 = "NumP"
  // 2 = "DenE"
  // 3 = "DenAE"

  TH1F* H_trk_trk_doubleCounter;
  TH1F* H_trk_trk_doubleCounter_eta;
  TH1F* H_trk_trk_doubleCounter_pt;

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
  myMap MapFit_d0;       // Map to store d0 difference [pt, eta]
  myMap MapFit_dz;       // Map to store dz difference [pt, eta]
  myMap MapFit_pt;       // Map to store pt difference [pt, eta]
  myMap MapFit_d0_pull;  // Map to store d0 pulls [pt, eta]
  myMap MapFit_dz_pull;  // Map to store dz pulls [pt, eta]
  myMap MapFit_pt_pull;  // Map to store pt pulls [pt, eta]
  myMap MapFit_d0err;    // Map to store generalTrack d0err [pt, eta]
  myMap MapFit_dzerr;    // Map to store generalTrack dzerr [pt, eta]
  myMap MapFit_pterr;    // Map to store generalTrack pterr [pt, eta]
  myMap MapFit_pixpterr; // Map to store generalTrack pterr [pt, eta]

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

  // @@@@@@ Pixel vertices @@@@@@
  TH1F* H_pix_vx_z;
  TH2F* H_pix_vx_xy;
  TH1F* H_pix_vx_xerr;
  TH1F* H_pix_vx_yerr;
  TH1F* H_pix_vx_zerr;
  TH1F* H_ntrk_pix_vx;
  TH1F* H_pix_vx_normchi2;
  TH1F* H_pix_vx_dof;

  // @@@@@@ Tracker vertices @@@@@@
  TH1F* H_trk_vx_z;
  TH2F* H_trk_vx_xy;
  TH1F* H_trk_vx_xerr;
  TH1F* H_trk_vx_yerr;
  TH1F* H_trk_vx_zerr;
  TH1F* H_ntrk_trk_vx;
  TH1F* H_trk_vx_normchi2;
  TH1F* H_trk_vx_dof;

  TH1F* H_trk_vx_z_nomatch;
  TH2F* H_trk_vx_xy_nomatch;
  TH1F* H_trk_vx_nomatch_trk_pt;
  TH1F* H_trk_vx_nomatch_trk_eta;
  TH1F* H_trk_vx_nomatch_normchi2;
  TH1F* H_trk_vx_nomatch_dof;

  // @@@@@@ Common pixel&tracker vertices @@@@@@
  TH1F* H_pix_vx_trk_pt;
  TH1F* H_pix_vx_trk_normchi2;
  TH1F* H_pix_vx_trk_pt_nomatch;
  TH1F* H_pix_vx_trk_eta_nomatch;
  TH1F* H_pix_vx_trk_normchi2_nomatch;

  TH1F* H_trk_vx_trk_pt;
  TH1F* H_trk_vx_trk_normchi2;
  TH1F* H_trk_vx_trk_pt_nomatch;
  TH1F* H_trk_vx_trk_eta_nomatch;
  TH1F* H_trk_vx_trk_normchi2_nomatch;
  TH1F* H_trk_vxrefit_trk_pt_nomatch;
  TH1F* H_trk_vxrefit_trk_eta_nomatch;

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
  map<string, TH1F*> map_H_vx_x_err2_orig;
  map<string, TH1F*> map_H_vx_y_err2_orig;
  map<string, TH1F*> map_H_vx_z_err2_orig;
  map<string, TH1F*> map_H_vx_x_err2_non_orig;
  map<string, TH1F*> map_H_vx_y_err2_non_orig;
  map<string, TH1F*> map_H_vx_z_err2_non_orig;

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


  // @@@@@@ Local variables @@@@@@
  DecisionWord MyAlgoMask;
  DecisionWord MyTechMask;
  vector<string> MyHLTMask;
  TString hTitle[5];
  TString hTitleEta[3];
  double EtaThr;       // --> Must have at most 2 decimal digits
  double EtaThrBin[3]; // --> Must have at most 2 decimal digits
  double d0Range;      // Unit: (cm)
  double dzRange;      // Unit: (cm)
  double ptRange;      // Unit: (GeV/c)
  double d0errRange;   // Unit: (um)^2
  double dzerrRange;   // Unit: (um)^2
  double pterrRange;   // Unit: (MeV/c)^2
  double pullsRange;
  double xRange;       // Unit: (um)
  double yRange;       // Unit: (um)
  double zRange;       // Unit: (um)
  double xerrRange;    // Unit: (um)^2
  double yerrRange;    // Unit: (um)^2
  double zerrRange;    // Unit: (um)^2
  int EvDiscardL1;
  int EvDiscardHLT;
  int EvPassL1;
  int EvPassHLT;
  int EvPassTotal;
  int EvDiscardTotal;
  int EvTotal;
  unsigned int NpixelLayers;
  TrackAssociatorBase* theTkAssociator;
  VertexAssociatorBase* theVxAssociator;
  edm::Service<TFileService> FileService;
  auto_ptr<TFileDirectory> SubDirPixTk;
  auto_ptr<TFileDirectory> SubDirPixTkCom;
  auto_ptr<TFileDirectory> SubDirPixTkFit;
  auto_ptr<TFileDirectory> SubDirPixVx;
  auto_ptr<TFileDirectory> SubDirPixVxNoM;
  auto_ptr<TFileDirectory> SubDirPixVxCom;
  auto_ptr<TFileDirectory> SubDirPixVxNtrks;


  // @@@@@@ Parameters from cfg file @@@@@@
  bool PrintMsg;
  bool IsMC;
  bool IsTrkPart;
  string AnalysisType;           // Used to chose: Vertex, Track, TrackANDVertex
  unsigned int MinTkTracks;      // Min number of tracks
  unsigned int MaxTkTracks;      // Max number of tracks
  unsigned int MinHitsMatch;     // Min number of valid hits that a pixel track must share with a tracker track
  double TkXYVxErrCorr;          // Correction factor for tracker vertex errors along XY
  double TkZVxErrCorr;           // Correction factor for tracker vertex errors along Z

  double MaxEtaTkTrk;            // Max eta that a track should have in order to be taken into account
                                 // --> Must have at most 2 decimal digits
  double MaxChi2PxTrk;           // Max chi2 for the acceptance of the pixel track
  double MaxChi2TkTrk;           // Max chi2 for the acceptance of the general track
  double RangePt;                // pt bins. Unit: (GeV/c) --> Must have at most 2 decimal digits
  double PtStep;                 // pt bins. Unit: (GeV/c) --> Must have at most 2 decimal digits
  double RangeEta;               // eta bins --> Must have at most 2 decimal digits
  double EtaStep;                // eta bins --> Must have at most 2 decimal digits
  double RangePhi;               // phi bins. Unit: (deg)
  double PhiStep;                // phi bins. Unit: (deg)
  double RangeChi2;              // chi2 bins
  double Chi2Step;               // chi2 bins
  double TrBound;                // [-r,r] bounding region for tracks in the transverse plane
  double TzBound;                // [-zB,zB] bounding region for tracks z coordinate
  unsigned int MinValidHitsPx;   // Min number of valid hits that a pixel track must have in order to be taken into account
  unsigned int MinValidHitsTk;   // Min number of valid hits that a tracker track must have in order to be taken into account
  double MinTrkVxDoF;            // Min number of tracks that a vertex must have in order to be used to compute dz and dxy
  double MinTrkVxWgt;            // Min weight of the vertex tracks
  double MinPtTk;                // Mini Pt that a track must have in order to be taken into account. Unit: (GeV/c)
                                 // --> Must have at most 2 decimal digits

  edm::InputTag VxInputTag;      // Used for vertexing
  double MaxEtaVxTrk;            // Max eta that a track should have in order to be taken into account
  double MaxChi2VxTrk;           // Max chi2 for the acceptance of the track
  double VrBound;                // [-r,r] bounding region for vertices in the transverse plane
  double VzBound;                // [-zB,zB] bounding region for vertex z coordinate
  double MinVxDoF;               // Min number of degree of freedom for vertex quality cut
  unsigned int MinVxTrkMatch;    // Min number of tracks for vertex matching
  double PxVxErrCorr;            // Correction factor for pixel vertex errors
  double MinPtVx;                // Mini Pt that a track must have in order to be taken into account. Unit: (GeV/c)

  edm::ParameterSet ParSetTrackAss;
};

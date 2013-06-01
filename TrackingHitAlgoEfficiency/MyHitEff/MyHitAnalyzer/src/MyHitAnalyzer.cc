// -*- C++ -*-
//
// Package:    MyHitAnalyzer
// Class:      MyHitAnalyzer
// 
/**\class MyHitAnalyzer MyHitAnalyzer.cc MyHitEff/MyHitAnalyzer/src/MyHitAnalyzer.cc

   Description: <one line class summary>
   Algorithm to estimate the hit collection efficiency of the Combinatorial Kalman Filter algorithm

   Comments:
   - The algorithm does not take into account the transient hits before the first hit of the track
   - Since the algorithm uses the same chi2 compatibility criteria it might have the same bias of the pattern recognition
   - In the SiPixelCluster the detId is set to zero !!!

   To do for cosmics:
   - test the effect of using the 1D covariance matrix "updator"

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Mauro Dinardo
//         Created:  Tue July 21 17:00:00 CEST 2009
// $Id: MyHitAnalyzer.cc,v 1.1 2010/10/06 18:30:37 dinardo Exp $

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

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"

#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"

#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"

#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/StripClusterParameterEstimator.h"
#include "RecoLocalTracker/SiStripRecHitConverter/interface/SiStripRecHitMatcher.h"

#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "RecoTracker/TransientTrackingRecHit/interface/TSiStripRecHit2DLocalPos.h"

#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>


using namespace std;
using namespace reco;


class MyHitAnalyzer : public edm::EDAnalyzer {
public:
  explicit MyHitAnalyzer (const edm::ParameterSet&);
  ~MyHitAnalyzer ();


private:
  virtual void beginRun (edm::Run const&, const edm::EventSetup& iSetup);
  virtual void analyze (const edm::Event&, const edm::EventSetup&);
  virtual void endJob ();

  // Histograms
  TH1F* H;
  TH1F* H_category_str;
  TH1F* H_category_pix;
  TH1F* H_category_multy_pix;
  TH1F* H_dist_pixe;
  TH1F* H_dist_rphi;
  TH1F* H_dist_ster;
  TH1F* H_length_pixe;
  TH1F* H_length_rphi;
  TH1F* H_length_ster;
  TH1F* H_chi2;  
  TH1F* H_p;
  TH1F* H_pt;
  TH1F* H_eta;
  TH1F* H_theta;
  TH1F* H_phi;
  TH1F* H_non_valid_hits;
  TH1F* H_non_sim_hits;
  TH1F* H_hitpos;
  TH1F* H_normhitpos;
  TH1F* H_rechits_on_mod_rphi;
  TH1F* H_rechits_on_mod_ster;
  TH1F* H_rechits_on_mod_pixe;
  TH1F* H_simhits_on_mod_rphi;
  TH1F* H_simhits_on_mod_ster;
  TH1F* H_simhits_on_mod_pixe;
  TH1F* H_str_module;
  TH1F* H_pix_charge;
  TH1F* H_pix_charge_doublehit_missed;
  TH1F* H_pix_charge_doublehit_trked;
  TH1F* H_pix_dist_doublehit;
  TH1F* H_str_dist_doublehit;
  TH1F* H_dist_str_samelayer;
  TH1F* H_dist_str_sameweel;
  TH1F* H_dist_pix_samelayer;
  TH1F* H_dist_pix_sameendcap;

  TH2F* H_geom;
  TH2F* H_normgeom;
  TH2F* H_chi2_length;
  TH2F* H_chi2_dist;
  TH2F* H_dist_length;
  TH2F* H_general_pix_module;
  TH2F* H_2x8_pix_module;
  TH2F* H_1x8_pix_module;
  TH2F* H_1x2_pix_module;
  TH2F* H_1x5_pix_module;
  TH2F* H_2x3_pix_module;
  TH2F* H_2x4_pix_module;
  TH2F* H_2x5_pix_module;

  TH2F* H_2DmapXY_B;
  TH2F* H_2DmapRZ_B;
  TH2F* H_2DmapXY_F;
  TH2F* H_2DmapRZ_F;
  
  // Histograms category-specific
  TH2F* H_2DmapXY_PixOver_B;
  TH2F* H_2DmapRZ_PixOver_B;
  TH2F* H_2DmapXY_PixOther_B;
  TH2F* H_2DmapRZ_PixOther_B;
  TH2F* H_2DmapXY_PixEdge_B;
  TH2F* H_2DmapRZ_PixEdge_B;
  TH2F* H_2DmapXY_PixTrans_B;
  TH2F* H_2DmapRZ_PixTrans_B;
  TH2F* H_2DmapXY_PixDoubleHit_B;
  TH2F* H_2DmapRZ_PixDoubleHit_B;

  TH2F* H_2DmapXY_PixOver_F;
  TH2F* H_2DmapRZ_PixOver_F;
  TH2F* H_2DmapXY_PixOther_F;
  TH2F* H_2DmapRZ_PixOther_F;
  TH2F* H_2DmapXY_PixEdge_F;
  TH2F* H_2DmapRZ_PixEdge_F;
  TH2F* H_2DmapXY_PixTrans_F;
  TH2F* H_2DmapRZ_PixTrans_F;
  TH2F* H_2DmapXY_PixDoubleHit_F;
  TH2F* H_2DmapRZ_PixDoubleHit_F;

  TH2F* H_2DmapXY_StrOver_B;
  TH2F* H_2DmapRZ_StrOver_B;
  TH2F* H_2DmapXY_StrLast_B;
  TH2F* H_2DmapRZ_StrLast_B;
  TH2F* H_2DmapXY_StrTrans_B;
  TH2F* H_2DmapRZ_StrTrans_B;
  TH2F* H_2DmapXY_StrOther_B;
  TH2F* H_2DmapRZ_StrOther_B;

  TH2F* H_2DmapXY_StrOver_F;
  TH2F* H_2DmapRZ_StrOver_F;
  TH2F* H_2DmapXY_StrLast_F;
  TH2F* H_2DmapRZ_StrLast_F;
  TH2F* H_2DmapXY_StrTrans_F;
  TH2F* H_2DmapRZ_StrTrans_F;
  TH2F* H_2DmapXY_StrOther_F;
  TH2F* H_2DmapRZ_StrOther_F;

  const TrackerGeometry* tracker;
  edm::ESHandle<MagneticField> theMF;
  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  edm::ESHandle<TransientTrackingRecHitBuilder> theBuilder;
  edm::ESHandle<StripClusterParameterEstimator> theCPE;
  TrackAssociatorBase* theAssociator;
  std::auto_ptr<TFileDirectory> SubDirClus;

  // Parameters from cfg file
  std::string trackRefitter_;
  bool PrintMsg;
  bool SimData;
  bool RealData;
  bool ClusterPlots;                  // It tells to the program whther it has to produce the plots of the clusters
  bool PixBackSearch;                 // Set the search of pixel hits on overlapping modules to be ouside --> in
  double MaxDistance;                 // Unit: [cm]
  double MaxDPhi;                     // Max change in phi during state propagation. Unit: [rad]. Default: 1.6
  double MaxChi2;                     // Max chi2 for the acceptance of the hit compatibility with the track
  double MinPt;                       // Mini Pt that a track must have in order to be taken into account. Unit: [GeV/c]
  double MaxEta;                      // Max eta that a track should have in order to be taken into account
  double LongitudinalDistanceOverlap; // Min distance along Z between overlapping modules. Unit: [cm]
  double RadialDistanceOverlap;       // Min distance along R between overlapping modules. Unit: [cm]
  int MinNValidHits;                  // Mini number of valid hits that a track must have in order to be taken into account
  int NHitCountOnLayers[6][18];       // Matrix containing rescaling factors: number of tracked hits per detector per layer
  unsigned int EdgeTollerance;        // Edge tollerance to consider a cluster an "edge cluster"
  unsigned int MaxTracks;
  edm::ParameterSet ParSetTrackAss;

  // Internal counters
  int NRphiModules;
  int NStereoModules;
  int NPixelModules;
  int NEventsPass;
  int PixelHit2x8Counter;
  int PixelHit1x8Counter;
  int PixelHit1x2Counter;
  int PixelHit1x5Counter;
  int PixelHit2x3Counter;
  int PixelHit2x4Counter;
  int PixelHit2x5Counter;

  int ParticleType;
  unsigned int ProcessType;
  double MaxPropLength;         // Unit: [cm]
  double MaxMomentum;           // Unit: [GeV/c]
  double MaxHitModule;          // Max value of the number of hits per module
  double BPixTollerance;        // Unit: [cm]
  map<int, const char*> TrackMap;
 
  typedef struct HitProperties_
  {
    unsigned int GeogId;
    unsigned int ParType;
    unsigned int ProcType;
    int IndexDet;
    int IndexHit;
    int LayerId;
    int TrackHitPosition;
    double PropLength;
    double Distance;
    double Chi2;
  } HitProperties;
};


MyHitAnalyzer::MyHitAnalyzer (const edm::ParameterSet& iConfig)
{
  trackRefitter_              = iConfig.getParameter<std::string>("trajectoryInput");
  PrintMsg                    = iConfig.getParameter<bool>("PrintMsg");
  SimData                     = iConfig.getParameter<bool>("SimData");
  RealData                    = iConfig.getParameter<bool>("RealData");
  ClusterPlots                = iConfig.getParameter<bool>("ClusterPlots");
  PixBackSearch               = iConfig.getParameter<bool>("PixBackSearch");
  MaxDistance                 = iConfig.getParameter<double>("MaxDistance");
  MaxDPhi                     = iConfig.getParameter<double>("MaxDPhi");
  MaxChi2                     = iConfig.getParameter<double>("MaxChi2");
  MinNValidHits               = iConfig.getParameter<int>("MinNValidHits");
  EdgeTollerance              = iConfig.getParameter<int>("EdgeTollerance");
  LongitudinalDistanceOverlap = iConfig.getParameter<double>("LongitudinalDistanceOverlap");
  RadialDistanceOverlap       = iConfig.getParameter<double>("RadialDistanceOverlap");
  MinPt                       = iConfig.getParameter<double>("MinPt");
  MaxEta                      = iConfig.getParameter<double>("MaxEta");
  MaxTracks                   = iConfig.getParameter<unsigned int>("MaxTracks");
  // ###### For simulated data ######
  if (SimData == true)
    ParSetTrackAss = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorByHitsPSet");

  // Internal counters
  NRphiModules   = 0;
  NStereoModules = 0;
  NPixelModules  = 0;
  NEventsPass    = 0;
  PixelHit2x8Counter = 0;
  PixelHit1x8Counter = 0;
  PixelHit1x2Counter = 0;
  PixelHit1x5Counter = 0;
  PixelHit2x3Counter = 0;
  PixelHit2x4Counter = 0;
  PixelHit2x5Counter = 0;

  ParticleType = 13; // Muons
  ProcessType = 2; // Primary
  MaxPropLength = 100.;
  MaxMomentum = 10.;
  MaxHitModule = 20.;
  BPixTollerance = 2.;
  TrackMap.insert(make_pair((int)PixelSubdetector::PixelBarrel,"PXB")); // PXB ID: 1
  TrackMap.insert(make_pair((int)PixelSubdetector::PixelEndcap,"PXE")); // PXE ID: 2
  TrackMap.insert(make_pair((int)SiStripDetId::TIB,"TIB")); // TIB ID: 3
  TrackMap.insert(make_pair((int)SiStripDetId::TID,"TID")); // TOB ID: 4
  TrackMap.insert(make_pair((int)SiStripDetId::TOB,"TOB")); // TID ID: 5
  TrackMap.insert(make_pair((int)SiStripDetId::TEC,"TEC")); // TEC ID: 6
}


MyHitAnalyzer::~MyHitAnalyzer ()
{
}


void MyHitAnalyzer::analyze (const edm::Event& iEvent,
			     const edm::EventSetup& iSetup)
{
  // ### Hit type: valid = 0, missing = 1, inactive = 2, bad = 3 ###
  // Internal variables
  bool belongs2Track;
  bool MomentumCounter;
  bool IsOther;
  bool IsSimHit;
  bool ToBeRecorded;
  bool DoubleHit;
  int TrackHitPosition;
  int CountHitsOnMod;
  int CountSimHitsOnMod;
  double PropLength;
  KFUpdator updator;
  vector<HitProperties> PixMissedHitsFw;   // Vector of pixel track-compatible hits during forward propagation: geographicalId, particle type, process type
  vector<HitProperties> PixMissedHitsBw;   // Vector of pixel track-compatible hits during backward propagation: geographicalId, particle type, process type
  vector<HitProperties> StrMissedHitsRphi; // Vector of rphi strip track-compatible hits during forward propagation: geographicalId, particle type, process type
  vector<HitProperties> StrMissedHitsSter; // Vector of stereo strip track-compatible hits during forward propagation: geographicalId, particle type, process type
  ReferenceCountingPointer<TransientTrackingRecHit> theTransientHit;
  TransientTrackingRecHit::RecHitPointer theTransientHitAlongB;
  TrajectoryStateOnSurface TrajStateOnSurf;
  TrajectoryStateOnSurface TrajStateOnNewSurf;
  const StripGeomDetUnit* theGeomDetStr;
  const PixelGeomDetUnit* theGeomDetPix;
  const StripTopology* RecStripTop;
  const RectangularPixelTopology* RecPixTop;

  // ###### For simulated data ######
  TrackerHitAssociator* hitAssociator = NULL;
  vector<unsigned int> TrParIDs;

  
  // ###### If RealData then check the correct run and bunch crossing ######
  //  if (RealData == true) if ((iEvent.id().run() != 124024) || (iEvent.bunchCrossing() != 51) || (iEvent.bunchCrossing() != 2724)) return;

  // ###### TRACKS ######
  edm::Handle<TrajTrackAssociationCollection> tracks;
  iEvent.getByLabel(trackRefitter_,tracks); // Refitter type

  // ###### VERTEX CONDITION ######
  if (RealData == true)
    {
      edm::Handle<reco::VertexCollection> primaryVertex;
      iEvent.getByLabel("offlinePrimaryVertices",primaryVertex);
      if (primaryVertex->size() == 0) return;
      if (primaryVertex->begin()->isFake()) return;
      if (primaryVertex->begin()->tracksSize() < 3) return;
    }

  // Hit Collections (it considers all the hits, also those not included in the track)
  edm::Handle<SiStripRecHit2DCollection> rechitsrphi;
  iEvent.getByLabel("siStripMatchedRecHits","rphiRecHit",rechitsrphi);
  edm::Handle<SiStripRecHit2DCollection> rechitsstereo;
  iEvent.getByLabel("siStripMatchedRecHits","stereoRecHit",rechitsstereo);
  edm::Handle<SiPixelRecHitCollection> rechitspixel;
  iEvent.getByLabel("siPixelRecHits",rechitspixel);

  // ###### For simulated data ######
  if (SimData == true)
    {
      edm::Handle<edm::View<reco::Track> > SomeTracks;
      iEvent.getByLabel("generalTracks",SomeTracks);
      edm::Handle<TrackingParticleCollection> TPC;
      iEvent.getByLabel("mergedtruth","MergedTrackTruth",TPC);

      hitAssociator = new TrackerHitAssociator::TrackerHitAssociator(iEvent, ParSetTrackAss);
      // Association between tracked tracks and simulated tracks
      reco::RecoToSimCollection Reco2Sim = theAssociator->associateRecoToSim(SomeTracks, TPC, &iEvent);

      // Extraction of the IDs of the tracked tracks associated with a simulated track
      // N.B.: a tracked particle that did bremsstrahlung might have two sim track IDs,
      // one for the track before the bremsstrahlung and one for the track after the bremsstrahlung
      for (unsigned int j = 0; j < SomeTracks->size(); j++) {
	edm::RefToBase<reco::Track> aTrack(SomeTracks, j);
	if (Reco2Sim.find(aTrack) != Reco2Sim.end())
	  {
	    vector<pair<TrackingParticleRef, double> > TrPar = Reco2Sim[aTrack];
	    if (TrPar.size() != 0)
	      {
		TrackingParticleRef TrParRef = TrPar.begin()->first;
		for (TrackingParticle::g4t_iterator g4T = TrParRef->g4Track_begin(); g4T != TrParRef->g4Track_end(); g4T++)
		  TrParIDs.push_back(g4T->trackId());
	      }
	  }
      }
    }

  // Loop over the tracks of the event:
  // - Look for a single track event in cosmics
  // - Look for at least one track in MC data
  if ((tracks->size() > 0) && (tracks->size() <= MaxTracks))
    {
      if (PrintMsg == true) cout << "\n\n### SELECTED EVENT --> " << iEvent.id() << " ###" << endl;
      NEventsPass++;

      // ########################
      // # LOOP OVER THE TRACKS #
      // ########################
      std::vector<TrajectoryMeasurement> Traj;
      for (TrajTrackAssociationCollection::const_iterator itTrack = tracks->begin(); itTrack != tracks->end(); itTrack++) {
	
	const Track& track = *itTrack->val;
	const Trajectory& Traj_ = *itTrack->key;
	Traj = Traj_.measurements();

	MomentumCounter = false;

	if ((track.quality(TrackBase::qualityByName("highPurity")) == true) &&
	    (track.numberOfValidHits() >= MinNValidHits) &&
	    (track.pt() >= MinPt) &&
	    (track.eta() >= -MaxEta) && (track.eta() <= MaxEta))
	  {
	    int layerId;
	    TransientTrackingRecHit::ConstRecHitPointer recHit;
	    // Loop over the hits of the track to count the number of tracked hits per detector and per layer
	    for (std::vector<TrajectoryMeasurement>::const_iterator itTraj = Traj.begin(); itTraj != Traj.end(); itTraj++) {
		  
	      // Extraction of the hits from the track
	      recHit = itTraj->recHit();

	      layerId = 0;
	      if (recHit->isValid() == true)
		{
		  if (recHit->geographicalId().subdetId() == (int)SiStripDetId::TIB) layerId = ((TIBDetId)(recHit->geographicalId())).layer();
		  if (recHit->geographicalId().subdetId() == (int)SiStripDetId::TOB) layerId = ((TOBDetId)(recHit->geographicalId())).layer();
		  if (recHit->geographicalId().subdetId() == (int)SiStripDetId::TID) layerId = ((TIDDetId)(recHit->geographicalId())).wheel() + 3*(((TIDDetId)(recHit->geographicalId())).side()-1);
		  if (recHit->geographicalId().subdetId() == (unsigned int)SiStripDetId::TEC) layerId = ((TECDetId)(recHit->geographicalId())).wheel() + 9*(((TECDetId)(recHit->geographicalId())).side()-1);
		  if (recHit->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) layerId = ((PXBDetId)(recHit->geographicalId())).layer();
		  if (recHit->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap) layerId = ((PXFDetId)(recHit->geographicalId())).disk() + 2*(((PXFDetId)(recHit->geographicalId())).side()-1);

		  NHitCountOnLayers[recHit->geographicalId().subdetId()-1][layerId-1]++;
		}
	    }

	    H->Fill("Total tracks",1.0); // Count the total number of tracks
	    H->Fill("Total valid hits",track.numberOfValidHits()); // Count the total number of valid hits in the track
	    H->Fill("Total lost hits",track.numberOfLostHits()); // Count the total number of invalid hits in the track
	    
	    if (PrintMsg == true)
	      {
		cout << "*** A track started ***" << endl;
		cout << "Number of hits: " << track.found();
		cout << "\tCharge: " << track.charge();
		cout << "\tdz: " << track.dz();
		cout << "\ttheta: " << track.theta();
		cout << "\teta: " << track.eta();
		cout << "\td0: " << track.d0();
		cout << "\tphi: " << track.phi();
		cout << "\tqoverp: " << track.qoverp();
		cout << "\tpt: " << track.pt();
		cout << "\tnumberOfLostHits: " << track.numberOfLostHits();
		cout << "\tnumberOfValidHits: " << track.numberOfValidHits() << endl;
		track.hitPattern().print();
	      }
	    

	    // ******************
	    // *** STRIP RPHI ***
	    // ******************

	    // ###############################################
	    // # SEARCH FOR TRACK-COMPATIBLE RPHI STRIP HITS #
	    // ###############################################
	    if (PrintMsg == true) cout << "\n*** SiStripTracker rphi hits ***" << endl;
	    for (SiStripRecHit2DCollection::const_iterator j = rechitsrphi->begin(); j != rechitsrphi->end(); j++) { // Loop on detectors
	      for (edmNew::DetSet<SiStripRecHit2D>::const_iterator h = j->begin(); h != j->end(); h++) { // Loop on hits on detector
		
		// The builder on the recHit add global information like: global position, magnetic field
		theTransientHit = theBuilder->build(&*h);
		
		PropLength = -1.0;
		belongs2Track = false;
		TrackHitPosition = -1;
		
		// Loop over the hits of the track
		for (std::vector<TrajectoryMeasurement>::const_iterator itTraj = Traj.begin(); itTraj != Traj.end(); itTraj++) {
		  
		  // Extraction of the hits from the track
		  recHit = itTraj->recHit();

		  // Check that the tracked hit is valid and belongs to the strip tracker
		  if ((recHit->isValid() == true) &&
		      ((recHit->geographicalId().subdetId() == (int)SiStripDetId::TIB) ||
		       (recHit->geographicalId().subdetId() == (int)SiStripDetId::TOB) ||
		       (recHit->geographicalId().subdetId() == (int)SiStripDetId::TID) ||
		       (recHit->geographicalId().subdetId() == (int)SiStripDetId::TEC)))
		    {
		      // Extract cluster information for the tracker hit
		      const SiStripRecHit2D* recHitStrip = dynamic_cast<const SiStripRecHit2D*>(recHit->hit());
		      SiStripRecHit2D::ClusterRef const& StripRecCluster = recHitStrip->cluster();
		      unsigned int FirstStripRecCluster = StripRecCluster->firstStrip();
		      unsigned int LastStripRecCluster = FirstStripRecCluster + StripRecCluster->amplitudes().size() - 1;
		      unsigned int IdStripRecCluster = recHitStrip->geographicalId()();
		      
		      // Extract cluster information for the transient hit
		      SiStripRecHit2D::ClusterRef const& StripTransCluster = h->cluster();
		      unsigned int FirstStripTransCluster = StripTransCluster->firstStrip();
		      unsigned int LastStripTransCluster = FirstStripTransCluster + StripTransCluster->amplitudes().size() - 1;
		      unsigned int IdStripTransCluster = h->geographicalId()();

		      // Does the TransiendRecHit not belong to the track?
		      if ((IdStripTransCluster != IdStripRecCluster) || (FirstStripTransCluster != FirstStripRecCluster) || (LastStripTransCluster != LastStripRecCluster))
			// 		      if (recHit->hit()->sharesInput(h,TrackingRecHit::all) == false)
			{
			  TrajectoryStateOnSurface TrajStateOnSurf = itTraj->forwardPredictedState();
			  
			  if (TrajStateOnSurf.isValid() == true)
			    {
			      // Propagation from the track hit to the theTransientRecHit surface
			      PropagatorWithMaterial thePropagator(alongMomentum, 0.1057, TrajStateOnSurf.magneticField(), MaxDPhi); // oppositeToMomentum, alongMomentum, anyDirection
			      TrajectoryStateOnSurface TrajStateOnNewSurf = thePropagator.propagate(TrajStateOnSurf, theTrackingGeometry->idToDet(h->geographicalId())->surface());
			      
			      double _PropLength = sqrt(powf((double)(TrajStateOnSurf.globalPosition().x() - theTransientHit->globalPosition().x()),2.)+
							powf((double)(TrajStateOnSurf.globalPosition().y() - theTransientHit->globalPosition().y()),2.)+
							powf((double)(TrajStateOnSurf.globalPosition().z() - theTransientHit->globalPosition().z()),2.));
			      			      
			      if (((PropLength < 0.0) || (_PropLength < PropLength)) && (TrajStateOnNewSurf.isValid() == true))
				{
				  PropLength = _PropLength;
				  TrackHitPosition = itTraj - Traj.begin();
				}
			    }
			}
		      else
			{
			  belongs2Track = true;
			  break;
			}
		    }
		}

		if ((TrackHitPosition != -1) && (belongs2Track == false))
		  {
		    KFUpdator updator;
		    TransientTrackingRecHit::ConstRecHitPointer recHit = Traj[TrackHitPosition].recHit();
		    TrajectoryStateOnSurface TrajStateOnSurf = updator.update(Traj[TrackHitPosition].forwardPredictedState(), *recHit);

		    // Propagation from the track hit to the theTransientRecHit surface
		    PropagatorWithMaterial thePropagator(alongMomentum, 0.1057, TrajStateOnSurf.magneticField(), MaxDPhi); // oppositeToMomentum, alongMomentum, anyDirection
		    TrajectoryStateOnSurface TrajStateOnNewSurf = thePropagator.propagate(TrajStateOnSurf, theTrackingGeometry->idToDet(h->geographicalId())->surface());

		    if (TrajStateOnNewSurf.isValid() == true)
		      {
			// Recalculate the hit position taking into account the direction of the magnetic field
			TransientTrackingRecHit::RecHitPointer theTransientHitAlongB = theTransientHit->clone(TrajStateOnNewSurf);

			double Distance = fabs(TrajStateOnNewSurf.localPosition().x() - theTransientHitAlongB->localPosition().x());
			
			if (Distance <= MaxDistance)
			  {			    
			    Chi2MeasurementEstimator Chi2Fw(MaxChi2);
			    MeasurementEstimator::HitReturnType Chi2FwResult = Chi2Fw.estimate(TrajStateOnNewSurf, *theTransientHitAlongB);
			    
			    if (Chi2FwResult.first == true)
			      {
				const Surface* surface = &(theTrackingGeometry->idToDet(h->geographicalId())->surface());
				GlobalPoint gp = surface->toGlobal(LocalPoint(0,0));
				int layerId = 0;
				if (h->geographicalId().subdetId() == (int)SiStripDetId::TIB) layerId = ((TIBDetId)(h->geographicalId())).layer();
				if (h->geographicalId().subdetId() == (int)SiStripDetId::TOB) layerId = ((TOBDetId)(h->geographicalId())).layer();
				if (h->geographicalId().subdetId() == (int)SiStripDetId::TID) layerId = ((TIDDetId)(h->geographicalId())).wheel() + 3*(((TIDDetId)(h->geographicalId())).side()-1);
				if (h->geographicalId().subdetId() == (int)SiStripDetId::TEC) layerId = ((TECDetId)(h->geographicalId())).wheel() + 9*(((TECDetId)(h->geographicalId())).side()-1);
				// side == 1 --> -Z; side == 2 --> +Z

				HitProperties HP;
				HP.GeogId = h->geographicalId()();
				HP.ParType = 0;
				HP.ProcType = 0;
				HP.IndexDet = j-rechitsrphi->begin();
				HP.IndexHit = h-j->begin();
				HP.LayerId = layerId;
				HP.TrackHitPosition = TrackHitPosition;
				HP.PropLength = PropLength;
				HP.Distance = Distance;
				HP.Chi2 = Chi2FwResult.second;
				
				if (SimData == true)
				  {
				    TransientTrackingRecHit::RecHitPointer theTH = theBuilder->build(&*h);				
				    // Extraction of the SimHits associated to a cluster with all their properties (e.g.
				    // particle type and process type)
				    // There might be several IDs because of the different strips/pixels in the cluster
				    std::vector<PSimHit> PSimTrackHits = hitAssociator->associateHit(*theTH);
				    if (PSimTrackHits.size() != 0)
				      {
					HP.ParType = abs(PSimTrackHits[0].particleType());
					HP.ProcType = PSimTrackHits[0].processType();
				      }
				  }
				
				unsigned int HIndex;
				for (HIndex = 0; HIndex < PixMissedHitsFw.size(); HIndex++)
				  if (StrMissedHitsRphi[HIndex].GeogId == HP.GeogId) break;
				
				// Filling histograms if it's a new hit or if it has a better chi2
				if ((HIndex == StrMissedHitsRphi.size()) || (HP.Chi2 < StrMissedHitsRphi[HIndex].Chi2))
				  {
				    if (HIndex == StrMissedHitsRphi.size()) StrMissedHitsRphi.push_back(HP);
				    else StrMissedHitsRphi[HIndex] = HP;
				  }
			      }	    
			  }
		      }
		    else if (PrintMsg == true) cout << "I was not able to re-propagate from track hit #" << TrackHitPosition << "\tPos.: " << TrajStateOnSurf.globalPosition() << " to TransRecHit #" << h-j->begin() << "\tID:" << theTransientHit->geographicalId()() << "\tPos.: " << theTransientHit->globalPosition() << endl;
		  }	
	      }
	    }

	    // ###################################################################################
	    // # FILLING HISTOGRAMS WITH DATA ASSOCIATED TO THE TRACK-COMPATIBLE STRIP RPHI HITS #
	    // ###################################################################################
	    for (unsigned int HIndex = 0; HIndex < StrMissedHitsRphi.size(); HIndex++) {

	      SiStripRecHit2DCollection::const_iterator j = rechitsrphi->begin() + StrMissedHitsRphi[HIndex].IndexDet;
	      edmNew::DetSet<SiStripRecHit2D>::const_iterator h = j->begin() + StrMissedHitsRphi[HIndex].IndexHit;
	      theTransientHit = theBuilder->build(&*h);
	      
	      recHit = Traj[StrMissedHitsRphi[HIndex].TrackHitPosition].recHit();
	      TrajStateOnSurf = updator.update(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition].forwardPredictedState(), *recHit);
	      
	      // Propagation from the track hit to the theTransientRecHit surface
	      PropagatorWithMaterial thePropagator(alongMomentum, 0.1057, TrajStateOnSurf.magneticField(), MaxDPhi); // oppositeToMomentum, alongMomentum, anyDirection
	      TrajStateOnNewSurf = thePropagator.propagate(TrajStateOnSurf, theTrackingGeometry->idToDet(h->geographicalId())->surface());
	      
	      // Recalculate the hit position taking into account the direction of the magnetic field
	      theTransientHitAlongB = theTransientHit->clone(TrajStateOnNewSurf);

	      // ######################
	      // # Hit categorization #
	      // ######################
	      // ###### For simulated data ######
	      IsSimHit = false;
	      if (SimData == true)
		{
		  // Extraction of the IDs associated to a hit (= cluster)
		  // There might be several IDs because of the different strips/pixels in the cluster
		  std::vector<SimHitIdpr> SimTrackIDs = hitAssociator->associateHitId(*theTransientHitAlongB);
		  unsigned int jj;
		  for (unsigned int ii = 0; ii < SimTrackIDs.size(); ii++) {
		    for (jj = 0; jj < TrParIDs.size(); jj++) {
		      if (SimTrackIDs[ii].first == TrParIDs[jj])
			{
			  IsSimHit = true;
			  break;
			}
		    }
		    if (jj != TrParIDs.size()) break;
		  }
		}
	      
	      IsOther = false;
	      DoubleHit = false;
	      theGeomDetStr = dynamic_cast<const StripGeomDetUnit*>(tracker->idToDet(h->geographicalId()));
	      RecStripTop = dynamic_cast<const StripTopology*>(&(theGeomDetStr->specificTopology()));
	      SiStripRecHit2D::ClusterRef const& StripTransCluster = h->cluster();
	      if ((SimData == true) && (IsSimHit == false))
		H_category_str->Fill("Non SimHit",1.0);
	      else if (h->geographicalId() == recHit->geographicalId())
		{
		  H_category_str->Fill("DoubleHit",1.0);
		  H_str_dist_doublehit->Fill(sqrt(powf((double)(recHit->localPosition().x()-h->localPosition().x()),2.) + powf((double)(recHit->localPosition().y()-h->localPosition().y()),2.)));
		  DoubleHit = true;
		}
	      else if ((h->geographicalId().subdetId() == recHit->geographicalId().subdetId()) &&
		       (((h->geographicalId().subdetId() == (int)SiStripDetId::TIB) &&
			 (((TIBDetId)(h->geographicalId())).layer() == ((TIBDetId)(recHit->geographicalId())).layer())) ||
			((h->geographicalId().subdetId() == (int)SiStripDetId::TOB) &&
			 (((TOBDetId)(h->geographicalId())).layer() == ((TOBDetId)(recHit->geographicalId())).layer())) ||
			((h->geographicalId().subdetId() == (int)SiStripDetId::TID) &&
			 (((TIDDetId)(h->geographicalId())).side() == ((TIDDetId)(recHit->geographicalId())).side()) &&
			 (((TIDDetId)(h->geographicalId())).wheel() == ((TIDDetId)(recHit->geographicalId())).wheel())) ||
			((h->geographicalId().subdetId() == (int)SiStripDetId::TEC) &&
			 (((TECDetId)(h->geographicalId())).side() == ((TECDetId)(recHit->geographicalId())).side()) &&
			 (((TECDetId)(h->geographicalId())).wheel() == ((TECDetId)(recHit->geographicalId())).wheel()))))
		{
		  if ((((h->geographicalId().subdetId() == (int)SiStripDetId::TIB) || (h->geographicalId().subdetId() == (int)SiStripDetId::TOB)) && 
		       (fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(recHit->globalPosition().x(),2.) + powf(recHit->globalPosition().y(),2.))) >= RadialDistanceOverlap)) ||
		      (((h->geographicalId().subdetId() == (int)SiStripDetId::TID) || (h->geographicalId().subdetId() == (int)SiStripDetId::TEC)) && 
		       (fabs(theTransientHitAlongB->globalPosition().z() - recHit->globalPosition().z()) >= LongitudinalDistanceOverlap)))
		    {
		      if ((((h->geographicalId().subdetId() == (int)SiStripDetId::TIB) || (h->geographicalId().subdetId() == (int)SiStripDetId::TOB)) && 
			   (fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(recHit->globalPosition().x(),2.) + powf(recHit->globalPosition().y(),2.))) >= RadialDistanceOverlap)))
			H_dist_str_samelayer->Fill(fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(recHit->globalPosition().x(),2.) + powf(recHit->globalPosition().y(),2.))));
		      else
			H_dist_str_sameweel->Fill(fabs(theTransientHitAlongB->globalPosition().z() - recHit->globalPosition().z()));

		      if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TIB) ||
			  (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TOB))
			{
			  H_2DmapXY_StrOver_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			  H_2DmapRZ_StrOver_B->Fill(theTransientHitAlongB->globalPosition().z(),
						    sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							 powf(theTransientHitAlongB->globalPosition().y(),2.)));
			}
		      else if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TID) ||
			       (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TEC))
			{
			  H_2DmapXY_StrOver_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			  H_2DmapRZ_StrOver_F->Fill(theTransientHitAlongB->globalPosition().z(),
						    sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							 powf(theTransientHitAlongB->globalPosition().y(),2.)));
			}
		      H_category_str->Fill("Overlap",1.0);
		    }
		}
	      else if ((((StrMissedHitsRphi[HIndex].TrackHitPosition - 1) >= 0) && (Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->isValid() == true)) &&
		       ((h->geographicalId().subdetId() == Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->geographicalId().subdetId()) &&
			(((h->geographicalId().subdetId() == (int)SiStripDetId::TIB) &&
			  (((TIBDetId)(h->geographicalId())).layer() == ((TIBDetId)(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->geographicalId())).layer())) ||
			 ((h->geographicalId().subdetId() == (int)SiStripDetId::TOB) &&
			  (((TOBDetId)(h->geographicalId())).layer() == ((TOBDetId)(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->geographicalId())).layer())) ||
			 ((h->geographicalId().subdetId() == (int)SiStripDetId::TID) &&
			  (((TIDDetId)(h->geographicalId())).side() == ((TIDDetId)(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->geographicalId())).side()) &&
			  (((TIDDetId)(h->geographicalId())).wheel() == ((TIDDetId)(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->geographicalId())).wheel())) ||
			 ((h->geographicalId().subdetId() == (int)SiStripDetId::TEC) &&
			  (((TECDetId)(h->geographicalId())).side() == ((TECDetId)(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->geographicalId())).side()) &&
			  (((TECDetId)(h->geographicalId())).wheel() == ((TECDetId)(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->geographicalId())).wheel())))))
		{
		  if ((((h->geographicalId().subdetId() == (int)SiStripDetId::TIB) || (h->geographicalId().subdetId() == (int)SiStripDetId::TOB)) && 
		       (fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->globalPosition().x(),2.) + powf(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->globalPosition().y(),2.))) >= RadialDistanceOverlap)) ||
		      (((h->geographicalId().subdetId() == (int)SiStripDetId::TID) || (h->geographicalId().subdetId() == (int)SiStripDetId::TEC)) && 
		       (fabs(theTransientHitAlongB->globalPosition().z() - Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->globalPosition().z()) >= LongitudinalDistanceOverlap)))
		    {
		      if ((((h->geographicalId().subdetId() == (int)SiStripDetId::TIB) || (h->geographicalId().subdetId() == (int)SiStripDetId::TOB)) && 
			   (fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->globalPosition().x(),2.) + powf(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->globalPosition().y(),2.))) >= RadialDistanceOverlap)))
			H_dist_str_samelayer->Fill(fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->globalPosition().x(),2.) + powf(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->globalPosition().y(),2.))));
		      else
			H_dist_str_sameweel->Fill(fabs(theTransientHitAlongB->globalPosition().z() - Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->globalPosition().z()));
		      
		      if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TIB) ||
			  (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TOB))
			{
			  H_2DmapXY_StrOver_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			  H_2DmapRZ_StrOver_B->Fill(theTransientHitAlongB->globalPosition().z(),
						    sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							 powf(theTransientHitAlongB->globalPosition().y(),2.)));
			}
		      else if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TID) ||
			       (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TEC))
			{
			  H_2DmapXY_StrOver_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			  H_2DmapRZ_StrOver_F->Fill(theTransientHitAlongB->globalPosition().z(),
						    sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							 powf(theTransientHitAlongB->globalPosition().y(),2.)));
			}
		      H_category_str->Fill("Overlap",1.0);
		    }
		}			    
	      else if ((((StrMissedHitsRphi[HIndex].TrackHitPosition - 1) >= 0) && (Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->isValid() == false)) &&
		       (((StrMissedHitsRphi[HIndex].TrackHitPosition - 2) >= 0) && (Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-2].recHit()->isValid() == true)) &&
		       ((h->geographicalId().subdetId() == Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-2].recHit()->geographicalId().subdetId()) &&
			(((h->geographicalId().subdetId() == (int)SiStripDetId::TIB) &&
			  (((TIBDetId)(h->geographicalId())).layer() == ((TIBDetId)(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-2].recHit()->geographicalId())).layer())) ||
			 ((h->geographicalId().subdetId() == (int)SiStripDetId::TOB) &&
			  (((TOBDetId)(h->geographicalId())).layer() == ((TOBDetId)(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-2].recHit()->geographicalId())).layer())) ||
			 ((h->geographicalId().subdetId() == (int)SiStripDetId::TID) &&
			  (((TIDDetId)(h->geographicalId())).side() == ((TIDDetId)(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-2].recHit()->geographicalId())).side()) &&
			  (((TIDDetId)(h->geographicalId())).wheel() == ((TIDDetId)(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-2].recHit()->geographicalId())).wheel())) ||
			 ((h->geographicalId().subdetId() == (int)SiStripDetId::TEC) &&
			  (((TECDetId)(h->geographicalId())).side() == ((TECDetId)(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-2].recHit()->geographicalId())).side()) &&
			  (((TECDetId)(h->geographicalId())).wheel() == ((TECDetId)(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-2].recHit()->geographicalId())).wheel())))))
		{
		  if ((((h->geographicalId().subdetId() == (int)SiStripDetId::TIB) || (h->geographicalId().subdetId() == (int)SiStripDetId::TOB)) && 
		       (fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-2].recHit()->globalPosition().x(),2.) + powf(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-2].recHit()->globalPosition().y(),2.))) >= RadialDistanceOverlap)) ||
		      (((h->geographicalId().subdetId() == (int)SiStripDetId::TID) || (h->geographicalId().subdetId() == (int)SiStripDetId::TEC)) && 
		       (fabs(theTransientHitAlongB->globalPosition().z() - Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-2].recHit()->globalPosition().z()) >= LongitudinalDistanceOverlap)))
		    {		     
		      if ((((h->geographicalId().subdetId() == (int)SiStripDetId::TIB) || (h->geographicalId().subdetId() == (int)SiStripDetId::TOB)) && 
			   (fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-2].recHit()->globalPosition().x(),2.) + powf(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-2].recHit()->globalPosition().y(),2.))) >= RadialDistanceOverlap)))
			H_dist_str_samelayer->Fill(fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-2].recHit()->globalPosition().x(),2.) + powf(Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-2].recHit()->globalPosition().y(),2.))));
		      else
			H_dist_str_sameweel->Fill(fabs(theTransientHitAlongB->globalPosition().z() - Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-2].recHit()->globalPosition().z()));
		      
		      if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TIB) ||
			  (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TOB))
			{
			  H_2DmapXY_StrOver_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			  H_2DmapRZ_StrOver_B->Fill(theTransientHitAlongB->globalPosition().z(),
						    sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							 powf(theTransientHitAlongB->globalPosition().y(),2.)));
			}
		      else if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TID) ||
			       (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TEC))
			{
			  H_2DmapXY_StrOver_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			  H_2DmapRZ_StrOver_F->Fill(theTransientHitAlongB->globalPosition().z(),
						    sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							 powf(theTransientHitAlongB->globalPosition().y(),2.)));
			}
		      H_category_str->Fill("Overlap",1.0);
		    }
		}
	      else if (StrMissedHitsRphi[HIndex].TrackHitPosition == 0)
		{
		  if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TIB) ||
		      (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TOB))
		    {
		      H_2DmapXY_StrLast_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_StrLast_B->Fill(theTransientHitAlongB->globalPosition().z(),
						sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						     powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  else if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TID) ||
			   (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TEC))
		    {
		      H_2DmapXY_StrLast_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_StrLast_F->Fill(theTransientHitAlongB->globalPosition().z(),
						sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						     powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  H_category_str->Fill("Last",1.0);
		}
	      else if ((StripTransCluster->firstStrip() <= EdgeTollerance) ||
		       ((StripTransCluster->firstStrip()+StripTransCluster->amplitudes().size()-1) >= (RecStripTop->nstrips()-1-EdgeTollerance)))
		H_category_str->Fill("Edge",1.0);
	      else if (h->geographicalId().subdetId() != recHit->geographicalId().subdetId())
		{
		  if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TIB) ||
		      (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TOB))
		    {
		      H_2DmapXY_StrTrans_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_StrTrans_B->Fill(theTransientHitAlongB->globalPosition().z(),
						 sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						      powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  else if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TID) ||
			   (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TEC))
		    {
		      H_2DmapXY_StrTrans_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_StrTrans_F->Fill(theTransientHitAlongB->globalPosition().z(),
						 sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						      powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  H_category_str->Fill("Aft.Trans.",1.0);
		}
	      else if ((((StrMissedHitsRphi[HIndex].TrackHitPosition - 1) >= 0) && (Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->isValid() == true) &&
			(h->geographicalId().subdetId() != Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->geographicalId().subdetId())) ||
		       ((((StrMissedHitsRphi[HIndex].TrackHitPosition - 1) >= 0) && (Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->isValid() == false) &&
			 (h->geographicalId().subdetId() == Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-1].recHit()->geographicalId().subdetId())) &&
			(((StrMissedHitsRphi[HIndex].TrackHitPosition - 2) >= 0) &&
			 (h->geographicalId().subdetId() != Traj[StrMissedHitsRphi[HIndex].TrackHitPosition-2].recHit()->geographicalId().subdetId()))))
		{
		  if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TIB) ||
		      (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TOB))
		    {
		      H_2DmapXY_StrTrans_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_StrTrans_B->Fill(theTransientHitAlongB->globalPosition().z(),
						 sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						      powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  else if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TID) ||
			   (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TEC))
		    {
		      H_2DmapXY_StrTrans_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_StrTrans_F->Fill(theTransientHitAlongB->globalPosition().z(),
						 sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						      powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  H_category_str->Fill("Bef.Trans.",1.0);
		}
	      else
		{
		  if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TIB) ||
		      (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TOB))
		    {
		      H_2DmapXY_StrOther_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_StrOther_B->Fill(theTransientHitAlongB->globalPosition().z(),
						 sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						      powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  else if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TID) ||
			   (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TEC))
		    {
		      H_2DmapXY_StrOther_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_StrOther_F->Fill(theTransientHitAlongB->globalPosition().z(),
						 sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						      powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  H_category_str->Fill("Others",1.0);
		  IsOther = true;
		}
				
	      if (DoubleHit == false)
		{
		  H_chi2->Fill(StrMissedHitsRphi[HIndex].Chi2);
		  H_chi2_length->Fill(StrMissedHitsRphi[HIndex].Chi2, StrMissedHitsRphi[HIndex].PropLength);
		  H_chi2_dist->Fill(StrMissedHitsRphi[HIndex].Chi2, StrMissedHitsRphi[HIndex].Distance);
		  H->Fill("R\\phi hits",1.0); // Count the total number of missed hits: not including the "DoubleHits"
		  H_dist_rphi->Fill(StrMissedHitsRphi[HIndex].Distance);
		  H_length_rphi->Fill(StrMissedHitsRphi[HIndex].PropLength);
		  H_dist_length->Fill(StrMissedHitsRphi[HIndex].Distance, StrMissedHitsRphi[HIndex].PropLength);
		  H_normgeom->Fill(TrackMap.find(h->geographicalId().subdetId())->second, StrMissedHitsRphi[HIndex].LayerId, 1.0);
		  H_hitpos->Fill(StrMissedHitsRphi[HIndex].TrackHitPosition);
		  H_normhitpos->Fill((double)StrMissedHitsRphi[HIndex].TrackHitPosition/(double)(Traj.size()-1));
		  
		  if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TIB) ||
		      (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TOB))
		    {
		      H_2DmapXY_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_B->Fill(theTransientHitAlongB->globalPosition().z(),
					sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
					 powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  else if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TID) ||
			   (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TEC))
		    {
		      H_2DmapXY_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_F->Fill(theTransientHitAlongB->globalPosition().z(),
					sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
					     powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  
		  if (ClusterPlots == true)
		    {
		      stringstream HistoName;
		      HistoName << "Rphi module ID: " << h->geographicalId()() << " #" << ++NPixelModules;
		      H_str_module = SubDirClus->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 769, -0.5, 768.5);
		      H_str_module->SetXTitle("Strips");
		      H_str_module->SetYTitle("Cluster ID [arb.units]");
		    }
		  CountHitsOnMod = 0;
		  CountSimHitsOnMod = 0;
		  for (SiStripRecHit2DCollection::const_iterator jIt = rechitsrphi->begin(); jIt != rechitsrphi->end(); jIt++) { // Loop on detectors
		    for (edmNew::DetSet<SiStripRecHit2D>::const_iterator hIt = jIt->begin(); hIt != jIt->end(); hIt++) { // Loop on hits on detector				   
		      if (h->geographicalId()() == hIt->geographicalId()())
			{
			  CountHitsOnMod++;
			  
			  // ###### For simulated data ######
			  if (SimData == true)
			    {
			      TransientTrackingRecHit::RecHitPointer theTH = theBuilder->build(&*hIt);
			      // Extraction of the IDs associated to a hit (= cluster)
			      // There might be several IDs because of the different strips/pixels in the cluster
			      std::vector<SimHitIdpr> SimTrackIDs = hitAssociator->associateHitId(*theTH);
			      unsigned int jj;
			      for (unsigned int ii = 0; ii < SimTrackIDs.size(); ii++) {
				for (jj = 0; jj < TrParIDs.size(); jj++) {
				  if (SimTrackIDs[ii].first == TrParIDs[jj])
				    {
				      CountSimHitsOnMod++;
				      
				      if (ClusterPlots == true)
					for (unsigned int k = 0; k < hIt->cluster()->amplitudes().size(); k++)
					  // H_str_module->Fill(hIt->cluster()->firstStrip()+k, hIt->cluster()->amplitudes()[k]);
					  H_str_module->Fill(hIt->cluster()->firstStrip()+k, CountSimHitsOnMod);
				      
				      break;
				    }
				}
				if (jj != TrParIDs.size()) break;
			      }
			    }
			}
		    }
		  }
		  H_rechits_on_mod_rphi->Fill(CountHitsOnMod);
		  H_simhits_on_mod_rphi->Fill(((double)CountSimHitsOnMod)/((double)CountHitsOnMod)*100.);
		  
		  MomentumCounter = true;
		  
		  if (PrintMsg == true)
		    {
		      cout << "\nI found an RPHI hit that doesn't belong to the track in the event: " <<  iEvent.id() << endl;
		      cout << "\tHit distance from the track: " << StrMissedHitsRphi[HIndex].Distance << " (<=" << MaxDistance << " cm)" << "\tPropagation length: " << StrMissedHitsRphi[HIndex].PropLength << endl;
		      cout << "\tChi2 compatibility with TransientRecHit: " << StrMissedHitsRphi[HIndex].Chi2 << " (<=" << MaxChi2 << ")" << endl;
		      cout << "\tPropagation from track hit #" << StrMissedHitsRphi[HIndex].TrackHitPosition << "\tTSOS pos.: " << TrajStateOnSurf.globalPosition() << "\tHit pos.: " << recHit->globalPosition() << endl;
		      cout << "\t\tSubDetId: " << recHit->geographicalId().subdetId() << "\tID: " << recHit->geographicalId()() << endl;
		      cout << "\tResult of propagation: " << TrajStateOnNewSurf.globalPosition() << endl;
		      cout << "\tTransRecHit #" << h-j->begin() << "\tID: " << theTransientHitAlongB->geographicalId()() << "\tPos.: " << theTransientHitAlongB->globalPosition() << endl;
		      cout << "\t\tSubDetId: " << h->geographicalId().subdetId() << "\tLayerID: " << StrMissedHitsRphi[HIndex].LayerId;
		      cout << "\tFirst strip: " << h->cluster()->firstStrip() << "\tLast strip: " << h->cluster()->firstStrip()+h->cluster()->amplitudes().size()-1 << endl;
		    }	
		}
	    }			    
	    

	    // ********************
	    // *** STRIP STEREO ***
	    // ********************

	    // #################################################
	    // # SEARCH FOR TRACK-COMPATIBLE STEREO STRIP HITS #
	    // #################################################
	    if (PrintMsg == true) cout << "\n*** SiStripTracker stereo hits ***" << endl;
	    for (SiStripRecHit2DCollection::const_iterator j = rechitsstereo->begin(); j != rechitsstereo->end(); j++) { // Loop on detectors
	      for (edmNew::DetSet<SiStripRecHit2D>::const_iterator h = j->begin(); h != j->end(); h++) { // Loop on hits on detector
		
		// The builder on the recHit add global information like: global position, magnetic field
		theTransientHit = theBuilder->build(&*h);

		PropLength = -1.0;
		belongs2Track = false;
		TrackHitPosition = -1;
		
		// Loop over the hits of the track
		for (std::vector<TrajectoryMeasurement>::const_iterator itTraj = Traj.begin(); itTraj != Traj.end(); itTraj++) {
		  
		  // Extraction of the hits from the track
		  recHit = itTraj->recHit();
 
		  // Check that the tracked hit is valid and belongs to the strip tracker
		  if ((recHit->isValid() == true) &&
		      ((recHit->geographicalId().subdetId() == (int)SiStripDetId::TIB) ||
		       (recHit->geographicalId().subdetId() == (int)SiStripDetId::TOB) ||
		       (recHit->geographicalId().subdetId() == (int)SiStripDetId::TID) ||
		       (recHit->geographicalId().subdetId() == (int)SiStripDetId::TEC)))
		    {
		      // Extract cluster information for the tracker hit
		      const SiStripRecHit2D* recHitStrip = dynamic_cast<const SiStripRecHit2D*>(recHit->hit());
		      SiStripRecHit2D::ClusterRef const& StripRecCluster = recHitStrip->cluster();
		      unsigned int FirstStripRecCluster = StripRecCluster->firstStrip();
		      unsigned int LastStripRecCluster = FirstStripRecCluster + StripRecCluster->amplitudes().size() - 1;
		      unsigned int IdStripRecCluster = recHitStrip->geographicalId()();
		      
		      // Extract cluster information for the transient hit
		      SiStripRecHit2D::ClusterRef const& StripTransCluster = h->cluster();
		      unsigned int FirstStripTransCluster = StripTransCluster->firstStrip();
		      unsigned int LastStripTransCluster = FirstStripTransCluster + StripTransCluster->amplitudes().size() - 1;
		      unsigned int IdStripTransCluster = h->geographicalId()();
		      
		      // Does the TransiendRecHit not belong to the track?
		      if ((IdStripTransCluster != IdStripRecCluster) || (FirstStripTransCluster != FirstStripRecCluster) || (LastStripTransCluster != LastStripRecCluster))
			// 		      if (recHit->hit()->sharesInput(h,TrackingRecHit::all) == false)
			{
			  TrajectoryStateOnSurface TrajStateOnSurf = itTraj->forwardPredictedState();
			  
			  if (TrajStateOnSurf.isValid() == true)
			    {
			      // Propagation from the track hit to the theTransientRecHit surface
			      PropagatorWithMaterial thePropagator(alongMomentum, 0.1057, TrajStateOnSurf.magneticField(), MaxDPhi); // oppositeToMomentum, alongMomentum, anyDirection
			      TrajectoryStateOnSurface TrajStateOnNewSurf = thePropagator.propagate(TrajStateOnSurf, theTrackingGeometry->idToDet(h->geographicalId())->surface());
			      
			      double _PropLength = sqrt(powf((double)(TrajStateOnSurf.globalPosition().x() - theTransientHit->globalPosition().x()),2.)+
							powf((double)(TrajStateOnSurf.globalPosition().y() - theTransientHit->globalPosition().y()),2.)+
							powf((double)(TrajStateOnSurf.globalPosition().z() - theTransientHit->globalPosition().z()),2.));
			      			      
			      if (((PropLength < 0.0) || (_PropLength < PropLength)) && (TrajStateOnNewSurf.isValid() == true))
				{
				  PropLength = _PropLength;
				  TrackHitPosition = itTraj - Traj.begin();
				}
			    }
			}
		      else
			{
			  belongs2Track = true;
			  break;
			}
		    }
		}
		
		if ((TrackHitPosition != -1) && (belongs2Track == false))
		  {
		    KFUpdator updator;
		    TransientTrackingRecHit::ConstRecHitPointer recHit = Traj[TrackHitPosition].recHit();
		    TrajectoryStateOnSurface TrajStateOnSurf = updator.update(Traj[TrackHitPosition].forwardPredictedState(), *recHit);

		    // Propagation from the track hit to the theTransientRecHit surface
		    PropagatorWithMaterial thePropagator(alongMomentum, 0.1057, TrajStateOnSurf.magneticField(), MaxDPhi); // oppositeToMomentum, alongMomentum, anyDirection
		    TrajectoryStateOnSurface TrajStateOnNewSurf = thePropagator.propagate(TrajStateOnSurf, theTrackingGeometry->idToDet(h->geographicalId())->surface());

		    if (TrajStateOnNewSurf.isValid() == true)
		      {
			// Recalculate the hit position taking into account the direction of the magnetic field
			TransientTrackingRecHit::RecHitPointer theTransientHitAlongB = theTransientHit->clone(TrajStateOnNewSurf);

			double Distance = fabs(TrajStateOnNewSurf.localPosition().x() - theTransientHitAlongB->localPosition().x());
			
			if (Distance <= MaxDistance)
			  {
			    Chi2MeasurementEstimator Chi2Fw(MaxChi2);
			    MeasurementEstimator::HitReturnType Chi2FwResult = Chi2Fw.estimate(TrajStateOnNewSurf, *theTransientHitAlongB);
			    
			    if (Chi2FwResult.first == true)
			      {
				const Surface* surface = &(theTrackingGeometry->idToDet(h->geographicalId())->surface());
				GlobalPoint gp = surface->toGlobal(LocalPoint(0,0));
				int layerId = 0;
				if (h->geographicalId().subdetId() == (int)SiStripDetId::TIB) layerId = ((TIBDetId)(h->geographicalId())).layer();
				if (h->geographicalId().subdetId() == (int)SiStripDetId::TOB) layerId = ((TOBDetId)(h->geographicalId())).layer();
				if (h->geographicalId().subdetId() == (int)SiStripDetId::TID) layerId = ((TIDDetId)(h->geographicalId())).wheel() + 3*(((TIDDetId)(h->geographicalId())).side()-1);
				if (h->geographicalId().subdetId() == (int)SiStripDetId::TEC) layerId = ((TECDetId)(h->geographicalId())).wheel() + 9*(((TECDetId)(h->geographicalId())).side()-1);
				// side == 1 --> -Z; side == 2 --> +Z
				
				HitProperties HP;
				HP.GeogId = h->geographicalId()();
				HP.ParType = 0;
				HP.ProcType = 0;
				HP.IndexDet = j-rechitsstereo->begin();
				HP.IndexHit = h-j->begin();
				HP.LayerId = layerId;
				HP.TrackHitPosition = TrackHitPosition;
				HP.PropLength = PropLength;
				HP.Distance = Distance;
				HP.Chi2 = Chi2FwResult.second;
				
				if (SimData == true)
				  {
				    TransientTrackingRecHit::RecHitPointer theTH = theBuilder->build(&*h);				
				    // Extraction of the SimHits associated to a cluster with all their properties (e.g.
				    // particle type and process type)
				    // There might be several IDs because of the different strips/pixels in the cluster
				    std::vector<PSimHit> PSimTrackHits = hitAssociator->associateHit(*theTH);
				    if (PSimTrackHits.size() != 0)
				      {
					HP.ParType = abs(PSimTrackHits[0].particleType());
					HP.ProcType = PSimTrackHits[0].processType();
				      }
				  }
				
				unsigned int HIndex;
				for (HIndex = 0; HIndex < PixMissedHitsFw.size(); HIndex++)
				  if (StrMissedHitsSter[HIndex].GeogId == HP.GeogId) break;
				
				// Filling histograms if it's a new hit or if it has a better chi2
				if ((HIndex == StrMissedHitsSter.size()) || (HP.Chi2 < StrMissedHitsSter[HIndex].Chi2))
				  {
				    if (HIndex == StrMissedHitsSter.size()) StrMissedHitsSter.push_back(HP);
				    else StrMissedHitsSter[HIndex] = HP;
				  }		
			      }	   
			  }
		      }
		    else if (PrintMsg == true) cout << "I was not able to re-propagate from track hit #" << TrackHitPosition << "\tPos.: " << TrajStateOnSurf.globalPosition() << " to TransRecHit #" << h-j->begin() << "\tID:" << theTransientHit->geographicalId()() << "\tPos.: " << theTransientHit->globalPosition() << endl;
		  }	
	      }
	    }
	    
	    // #####################################################################################
	    // # FILLING HISTOGRAMS WITH DATA ASSOCIATED TO THE TRACK-COMPATIBLE STRIP STEREO HITS #
	    // #####################################################################################
	    for (unsigned int HIndex = 0; HIndex < StrMissedHitsSter.size(); HIndex++) {

	      SiStripRecHit2DCollection::const_iterator j = rechitsstereo->begin() + StrMissedHitsSter[HIndex].IndexDet;
	      edmNew::DetSet<SiStripRecHit2D>::const_iterator h = j->begin() + StrMissedHitsSter[HIndex].IndexHit;
	      theTransientHit = theBuilder->build(&*h);
	      
	      recHit = Traj[StrMissedHitsSter[HIndex].TrackHitPosition].recHit();
	      TrajStateOnSurf = updator.update(Traj[StrMissedHitsSter[HIndex].TrackHitPosition].forwardPredictedState(), *recHit);
	      
	      // Propagation from the track hit to the theTransientRecHit surface
	      PropagatorWithMaterial thePropagator(alongMomentum, 0.1057, TrajStateOnSurf.magneticField(), MaxDPhi); // oppositeToMomentum, alongMomentum, anyDirection
	      TrajStateOnNewSurf = thePropagator.propagate(TrajStateOnSurf, theTrackingGeometry->idToDet(h->geographicalId())->surface());
	      
	      // Recalculate the hit position taking into account the direction of the magnetic field
	      theTransientHitAlongB = theTransientHit->clone(TrajStateOnNewSurf);

	      // ######################
	      // # Hit categorization #
	      // ######################
	      // ###### For simulated data ######
	      IsSimHit = false;
	      if (SimData == true)
		{
		  // Extraction of the IDs associated to a hit (= cluster)
		  // There might be several IDs because of the different strips/pixels in the cluster
		  std::vector<SimHitIdpr> SimTrackIDs = hitAssociator->associateHitId(*theTransientHitAlongB);
		  unsigned int jj;
		  for (unsigned int ii = 0; ii < SimTrackIDs.size(); ii++) {
		    for (jj = 0; jj < TrParIDs.size(); jj++) {
		      if (SimTrackIDs[ii].first == TrParIDs[jj])
			{
			  IsSimHit = true;
			  break;
			}
		    }
		    if (jj != TrParIDs.size()) break;
		  }
		}
	      
	      IsOther = false;
	      DoubleHit = false;
	      theGeomDetStr = dynamic_cast<const StripGeomDetUnit*>(tracker->idToDet(h->geographicalId()));
	      RecStripTop = dynamic_cast<const StripTopology*>(&(theGeomDetStr->specificTopology()));
	      SiStripRecHit2D::ClusterRef const& StripTransCluster = h->cluster();
	      if ((SimData == true) && (IsSimHit == false))
		H_category_str->Fill("Non SimHit",1.0);
	      else if (h->geographicalId() == recHit->geographicalId())
		{
		  H_category_str->Fill("DoubleHit",1.0);
		  H_str_dist_doublehit->Fill(sqrt(powf((double)(recHit->localPosition().x()-h->localPosition().x()),2.) + powf((double)(recHit->localPosition().y()-h->localPosition().y()),2.)));
		  DoubleHit = true;
		}
	      else if ((h->geographicalId().subdetId() == recHit->geographicalId().subdetId()) &&
		       (((h->geographicalId().subdetId() == (int)SiStripDetId::TIB) &&
			 (((TIBDetId)(h->geographicalId())).layer() == ((TIBDetId)(recHit->geographicalId())).layer())) ||
			((h->geographicalId().subdetId() == (int)SiStripDetId::TOB) &&
			 (((TOBDetId)(h->geographicalId())).layer() == ((TOBDetId)(recHit->geographicalId())).layer())) ||
			((h->geographicalId().subdetId() == (int)SiStripDetId::TID) &&
			 (((TIDDetId)(h->geographicalId())).side() == ((TIDDetId)(recHit->geographicalId())).side()) &&
			 (((TIDDetId)(h->geographicalId())).wheel() == ((TIDDetId)(recHit->geographicalId())).wheel())) ||
			((h->geographicalId().subdetId() == (int)SiStripDetId::TEC) &&
			 (((TECDetId)(h->geographicalId())).side() == ((TECDetId)(recHit->geographicalId())).side()) &&
			 (((TECDetId)(h->geographicalId())).wheel() == ((TECDetId)(recHit->geographicalId())).wheel()))))
		{
		  if ((((h->geographicalId().subdetId() == (int)SiStripDetId::TIB) || (h->geographicalId().subdetId() == (int)SiStripDetId::TOB)) && 
		       (fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(recHit->globalPosition().x(),2.) + powf(recHit->globalPosition().y(),2.))) >= RadialDistanceOverlap)) ||
		      (((h->geographicalId().subdetId() == (int)SiStripDetId::TID) || (h->geographicalId().subdetId() == (int)SiStripDetId::TEC)) && 
		       (fabs(theTransientHitAlongB->globalPosition().z() - recHit->globalPosition().z()) >= LongitudinalDistanceOverlap)))
		    {
		      if ((((h->geographicalId().subdetId() == (int)SiStripDetId::TIB) || (h->geographicalId().subdetId() == (int)SiStripDetId::TOB)) && 
			   (fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(recHit->globalPosition().x(),2.) + powf(recHit->globalPosition().y(),2.))) >= RadialDistanceOverlap)))
			H_dist_str_samelayer->Fill(fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(recHit->globalPosition().x(),2.) + powf(recHit->globalPosition().y(),2.))));
		      else
			H_dist_str_sameweel->Fill(fabs(theTransientHitAlongB->globalPosition().z() - recHit->globalPosition().z()));
		      
		      if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TIB) ||
			  (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TOB))
			{
			  H_2DmapXY_StrOver_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			  H_2DmapRZ_StrOver_B->Fill(theTransientHitAlongB->globalPosition().z(),
						    sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							 powf(theTransientHitAlongB->globalPosition().y(),2.)));
			}
		      else if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TID) ||
			       (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TEC))
			{
			  H_2DmapXY_StrOver_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			  H_2DmapRZ_StrOver_F->Fill(theTransientHitAlongB->globalPosition().z(),
						    sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							 powf(theTransientHitAlongB->globalPosition().y(),2.)));
			}
		      H_category_str->Fill("Overlap",1.0);
		    }
		}
	      else if ((((StrMissedHitsSter[HIndex].TrackHitPosition - 1) >= 0) && (Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->isValid() == true)) &&
		       ((h->geographicalId().subdetId() == Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->geographicalId().subdetId()) &&
			(((h->geographicalId().subdetId() == (int)SiStripDetId::TIB) &&
			  (((TIBDetId)(h->geographicalId())).layer() == ((TIBDetId)(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->geographicalId())).layer())) ||
			 ((h->geographicalId().subdetId() == (int)SiStripDetId::TOB) &&
			  (((TOBDetId)(h->geographicalId())).layer() == ((TOBDetId)(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->geographicalId())).layer())) ||
			 ((h->geographicalId().subdetId() == (int)SiStripDetId::TID) &&
			  (((TIDDetId)(h->geographicalId())).side() == ((TIDDetId)(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->geographicalId())).side()) &&
			  (((TIDDetId)(h->geographicalId())).wheel() == ((TIDDetId)(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->geographicalId())).wheel())) ||
			 ((h->geographicalId().subdetId() == (int)SiStripDetId::TEC) &&
			  (((TECDetId)(h->geographicalId())).side() == ((TECDetId)(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->geographicalId())).side()) &&
			  (((TECDetId)(h->geographicalId())).wheel() == ((TECDetId)(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->geographicalId())).wheel())))))
		{
		  if ((((h->geographicalId().subdetId() == (int)SiStripDetId::TIB) || (h->geographicalId().subdetId() == (int)SiStripDetId::TOB)) && 
		       (fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->globalPosition().x(),2.) + powf(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->globalPosition().y(),2.))) >= RadialDistanceOverlap)) ||
		      (((h->geographicalId().subdetId() == (int)SiStripDetId::TID) || (h->geographicalId().subdetId() == (int)SiStripDetId::TEC)) && 
		       (fabs(theTransientHitAlongB->globalPosition().z() - Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->globalPosition().z()) >= LongitudinalDistanceOverlap)))
		    {
		      if ((((h->geographicalId().subdetId() == (int)SiStripDetId::TIB) || (h->geographicalId().subdetId() == (int)SiStripDetId::TOB)) && 
			   (fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->globalPosition().x(),2.) + powf(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->globalPosition().y(),2.))) >= RadialDistanceOverlap)))
			H_dist_str_samelayer->Fill(fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->globalPosition().x(),2.) + powf(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->globalPosition().y(),2.))));
		      else
			H_dist_str_sameweel->Fill(fabs(theTransientHitAlongB->globalPosition().z() - Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->globalPosition().z()));

		      if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TIB) ||
			  (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TOB))
			{
			  H_2DmapXY_StrOver_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			  H_2DmapRZ_StrOver_B->Fill(theTransientHitAlongB->globalPosition().z(),
						    sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							 powf(theTransientHitAlongB->globalPosition().y(),2.)));
			}
		      else if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TID) ||
			       (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TEC))
			{
			  H_2DmapXY_StrOver_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			  H_2DmapRZ_StrOver_F->Fill(theTransientHitAlongB->globalPosition().z(),
						    sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							 powf(theTransientHitAlongB->globalPosition().y(),2.)));
			}
		      H_category_str->Fill("Overlap",1.0);
		    }
		}
	      else if ((((StrMissedHitsSter[HIndex].TrackHitPosition - 1) >= 0) && (Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->isValid() == false)) &&
		       (((StrMissedHitsSter[HIndex].TrackHitPosition - 2) >= 0) && (Traj[StrMissedHitsSter[HIndex].TrackHitPosition-2].recHit()->isValid() == true)) &&
		       ((h->geographicalId().subdetId() == Traj[StrMissedHitsSter[HIndex].TrackHitPosition-2].recHit()->geographicalId().subdetId()) &&
			(((h->geographicalId().subdetId() == (int)SiStripDetId::TIB) &&
			  (((TIBDetId)(h->geographicalId())).layer() == ((TIBDetId)(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-2].recHit()->geographicalId())).layer())) ||
			 ((h->geographicalId().subdetId() == (int)SiStripDetId::TOB) &&
			  (((TOBDetId)(h->geographicalId())).layer() == ((TOBDetId)(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-2].recHit()->geographicalId())).layer())) ||
			 ((h->geographicalId().subdetId() == (int)SiStripDetId::TID) &&
			  (((TIDDetId)(h->geographicalId())).side() == ((TIDDetId)(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-2].recHit()->geographicalId())).side()) &&
			  (((TIDDetId)(h->geographicalId())).wheel() == ((TIDDetId)(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-2].recHit()->geographicalId())).wheel())) ||
			 ((h->geographicalId().subdetId() == (int)SiStripDetId::TEC) &&
			  (((TECDetId)(h->geographicalId())).side() == ((TECDetId)(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-2].recHit()->geographicalId())).side()) &&
			  (((TECDetId)(h->geographicalId())).wheel() == ((TECDetId)(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-2].recHit()->geographicalId())).wheel())))))
		{
		  if ((((h->geographicalId().subdetId() == (int)SiStripDetId::TIB) || (h->geographicalId().subdetId() == (int)SiStripDetId::TOB)) && 
		       (fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-2].recHit()->globalPosition().x(),2.) + powf(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-2].recHit()->globalPosition().y(),2.))) >= RadialDistanceOverlap)) ||
		      (((h->geographicalId().subdetId() == (int)SiStripDetId::TID) || (h->geographicalId().subdetId() == (int)SiStripDetId::TEC)) && 
		       (fabs(theTransientHitAlongB->globalPosition().z() - Traj[StrMissedHitsSter[HIndex].TrackHitPosition-2].recHit()->globalPosition().z()) >= LongitudinalDistanceOverlap)))
		    {
		      if ((((h->geographicalId().subdetId() == (int)SiStripDetId::TIB) || (h->geographicalId().subdetId() == (int)SiStripDetId::TOB)) && 
			   (fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-2].recHit()->globalPosition().x(),2.) + powf(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-2].recHit()->globalPosition().y(),2.))) >= RadialDistanceOverlap)))
			H_dist_str_samelayer->Fill(fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-2].recHit()->globalPosition().x(),2.) + powf(Traj[StrMissedHitsSter[HIndex].TrackHitPosition-2].recHit()->globalPosition().y(),2.))));
		      else
			H_dist_str_sameweel->Fill(fabs(theTransientHitAlongB->globalPosition().z() - Traj[StrMissedHitsSter[HIndex].TrackHitPosition-2].recHit()->globalPosition().z()));
		      
		      if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TIB) ||
			  (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TOB))
			{
			  H_2DmapXY_StrOver_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			  H_2DmapRZ_StrOver_B->Fill(theTransientHitAlongB->globalPosition().z(),
						    sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							 powf(theTransientHitAlongB->globalPosition().y(),2.)));
			}
		      else if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TID) ||
			       (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TEC))
			{
			  H_2DmapXY_StrOver_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			  H_2DmapRZ_StrOver_F->Fill(theTransientHitAlongB->globalPosition().z(),
						    sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							 powf(theTransientHitAlongB->globalPosition().y(),2.)));
			}
		      H_category_str->Fill("Overlap",1.0);
		    }
		}
	      else if (StrMissedHitsSter[HIndex].TrackHitPosition == 0)
		{
		  if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TIB) ||
		      (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TOB))
		    {
		      H_2DmapXY_StrLast_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_StrLast_B->Fill(theTransientHitAlongB->globalPosition().z(),
						sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						     powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  else if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TID) ||
			   (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TEC))
		    {
		      H_2DmapXY_StrLast_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_StrLast_F->Fill(theTransientHitAlongB->globalPosition().z(),
						sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						     powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  H_category_str->Fill("Last",1.0);
		}
	      else if ((StripTransCluster->firstStrip() <= EdgeTollerance) ||
		       ((StripTransCluster->firstStrip()+StripTransCluster->amplitudes().size()-1) >= (RecStripTop->nstrips()-1-EdgeTollerance)))
		H_category_str->Fill("Edge",1.0);
	      else if (h->geographicalId().subdetId() != recHit->geographicalId().subdetId())
		{
		  if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TIB) ||
		      (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TOB))
		    {
		      H_2DmapXY_StrTrans_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_StrTrans_B->Fill(theTransientHitAlongB->globalPosition().z(),
						 sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						      powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  else if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TID) ||
			   (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TEC))
		    {
		      H_2DmapXY_StrTrans_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_StrTrans_F->Fill(theTransientHitAlongB->globalPosition().z(),
						 sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						      powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  H_category_str->Fill("Aft.Trans.",1.0);
		}
	      else if	((((StrMissedHitsSter[HIndex].TrackHitPosition - 1) >= 0) && (Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->isValid() == true) &&
			  (h->geographicalId().subdetId() != Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->geographicalId().subdetId())) ||
			 ((((StrMissedHitsSter[HIndex].TrackHitPosition - 1) >= 0) && (Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->isValid() == false) &&
			   (h->geographicalId().subdetId() == Traj[StrMissedHitsSter[HIndex].TrackHitPosition-1].recHit()->geographicalId().subdetId())) &&
			  (((StrMissedHitsSter[HIndex].TrackHitPosition - 2) >= 0) &&
			   (h->geographicalId().subdetId() != Traj[StrMissedHitsSter[HIndex].TrackHitPosition-2].recHit()->geographicalId().subdetId()))))
		{
		  if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TIB) ||
		      (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TOB))
		    {
		      H_2DmapXY_StrTrans_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_StrTrans_B->Fill(theTransientHitAlongB->globalPosition().z(),
						 sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						      powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  else if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TID) ||
			   (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TEC))
		    {
		      H_2DmapXY_StrTrans_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_StrTrans_F->Fill(theTransientHitAlongB->globalPosition().z(),
						 sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						      powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  H_category_str->Fill("Bef.Trans.",1.0);
		}
	      else
		{
		  if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TIB) ||
		      (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TOB))
		    {
		      H_2DmapXY_StrOther_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_StrOther_B->Fill(theTransientHitAlongB->globalPosition().z(),
						 sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						      powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  else if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TID) ||
			   (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TEC))
		    {
		      H_2DmapXY_StrOther_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_StrOther_F->Fill(theTransientHitAlongB->globalPosition().z(),
						 sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						      powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  H_category_str->Fill("Others",1.0);
		  IsOther = true;
		}
				
	      if (DoubleHit == false)
		{
		  H_chi2->Fill(StrMissedHitsSter[HIndex].Chi2);
		  H_chi2_length->Fill(StrMissedHitsSter[HIndex].Chi2, StrMissedHitsSter[HIndex].PropLength);
		  H_chi2_dist->Fill(StrMissedHitsSter[HIndex].Chi2, StrMissedHitsSter[HIndex].Distance);
		  H->Fill("Stereo hits",1.0); // Count the total number of missed hits: not including the "DoubleHits"
		  H_dist_ster->Fill(StrMissedHitsSter[HIndex].Distance);
		  H_length_ster->Fill(StrMissedHitsSter[HIndex].PropLength);
		  H_dist_length->Fill(StrMissedHitsSter[HIndex].Distance, StrMissedHitsSter[HIndex].PropLength);
		  H_normgeom->Fill(TrackMap.find(h->geographicalId().subdetId())->second, StrMissedHitsSter[HIndex].LayerId, 1.0);
		  H_hitpos->Fill(StrMissedHitsSter[HIndex].TrackHitPosition);
		  H_normhitpos->Fill((double)StrMissedHitsSter[HIndex].TrackHitPosition/(double)(Traj.size()-1));
		  
		  if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TIB) ||
		      (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TOB))
		    {
		      H_2DmapXY_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_B->Fill(theTransientHitAlongB->globalPosition().z(),
					sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
					     powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  else if ((theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TID) ||
			   (theTransientHitAlongB->geographicalId().subdetId() == (int)SiStripDetId::TEC))
		    {
		      H_2DmapXY_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_F->Fill(theTransientHitAlongB->globalPosition().z(),
					sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
					     powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  
		  if (ClusterPlots == true)
		    {
		      stringstream HistoName;				
		      HistoName << "Stereo module ID: " << h->geographicalId()() << " #" << ++NPixelModules;
		      H_str_module = SubDirClus->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 769, -0.5, 768.5);
		      H_str_module->SetXTitle("Strips");
		      H_str_module->SetYTitle("Charge [ADC]");
		    }
		  CountHitsOnMod = 0;
		  CountSimHitsOnMod = 0;
		  for (SiStripRecHit2DCollection::const_iterator jIt = rechitsstereo->begin(); jIt != rechitsstereo->end(); jIt++) { // Loop on detectors
		    for (edmNew::DetSet<SiStripRecHit2D>::const_iterator hIt = jIt->begin(); hIt != jIt->end(); hIt++) { // Loop on hits on detector
		      if (h->geographicalId()() == hIt->geographicalId()())
			{
			  CountHitsOnMod++;
			  
			  // ###### For simulated data ######
			  if (SimData == true)
			    {
			      TransientTrackingRecHit::RecHitPointer theTH = theBuilder->build(&*hIt);
			      // Extraction of the IDs associated to a hit (= cluster)
			      // There might be several IDs because of the different strips/pixels in the cluster
			      std::vector<SimHitIdpr> SimTrackIDs = hitAssociator->associateHitId(*theTH);
			      unsigned int jj;
			      for (unsigned int ii = 0; ii < SimTrackIDs.size(); ii++) {
				for (jj = 0; jj < TrParIDs.size(); jj++) {
				  if (SimTrackIDs[ii].first == TrParIDs[jj])
				    {
				      CountSimHitsOnMod++;
				      
				      if (ClusterPlots == true)
					for (unsigned int k = 0; k < hIt->cluster()->amplitudes().size(); k++)
					  // H_str_module->Fill(hIt->cluster()->firstStrip()+k, hIt->cluster()->amplitudes()[k]);
					  H_str_module->Fill(hIt->cluster()->firstStrip()+k, CountSimHitsOnMod);
				      
				      break;
				    }
				}
				if (jj != TrParIDs.size()) break;
			      }
			    }
			}
		    }
		  }
		  H_rechits_on_mod_ster->Fill(CountHitsOnMod);
		  H_simhits_on_mod_ster->Fill(((double)CountSimHitsOnMod)/((double)CountHitsOnMod)*100.);
		  
		  MomentumCounter = true;
		  
		  if (PrintMsg == true)
		    {
		      cout << "\nI found a STEREO hit that doesn't belong to the track in the event: " <<  iEvent.id() << endl;
		      cout << "\tHit distance from the track: " << StrMissedHitsSter[HIndex].Distance << " (<=" << MaxDistance << " cm)" << "\tPropagation length: " << StrMissedHitsSter[HIndex].PropLength << endl;
		      cout << "\tChi2 compatibility with TransientRecHit: " << StrMissedHitsSter[HIndex].Chi2 << " (<=" << MaxChi2 << ")" << endl;
		      cout << "\tPropagation from track hit #" << StrMissedHitsSter[HIndex].TrackHitPosition << "\tTSOS pos.: " << TrajStateOnSurf.globalPosition() << "\tHit pos.: " << recHit->globalPosition() << endl;
		      cout << "\t\tSubDetId: " << recHit->geographicalId().subdetId() << "\tID: " << recHit->geographicalId()() << endl;
		      cout << "\tResult of propagation: " << TrajStateOnNewSurf.globalPosition() << endl;
		      cout << "\tTransRecHit #" << h-j->begin() << "\tID: " << theTransientHitAlongB->geographicalId()() << "\tPos.: " << theTransientHitAlongB->globalPosition() << endl;
		      cout << "\t\tSubDetId: " << h->geographicalId().subdetId() << "\tLayerID: " << StrMissedHitsSter[HIndex].LayerId;
		      cout << "\tFirst strip: " << h->cluster()->firstStrip() << "\tLast strip: " << h->cluster()->firstStrip()+h->cluster()->amplitudes().size()-1 << endl;
		    }	
		}			
	    }
	    
	    
	    // **************
	    // *** PIXELS ***
	    // **************
	    
	    // ##########################################
	    // # SEARCH FOR TRACK-COMPATIBLE PIXEL HITS #
	    // ##########################################
	    if (PrintMsg == true) cout << "\n*** PixelTracker hits ***" << endl;
	    for (SiPixelRecHitCollection::const_iterator j = rechitspixel->begin(); j != rechitspixel->end(); j++) { // Loop on detectors
	      for (edmNew::DetSet<SiPixelRecHit>::const_iterator h = j->begin(); h != j->end(); h++) { // Loop on hits on detector
		  
		// The builder on the recHit add global information like: global position, magnetic field
		theTransientHit = theBuilder->build(&*h);
		
		PropLength = -1.0;
		belongs2Track = false;
		TrackHitPosition = -1;
		
		// Loop over the hits of the track
		for (std::vector<TrajectoryMeasurement>::const_iterator itTraj = Traj.begin(); itTraj != Traj.end(); itTraj++) {
		  
		  // Extraction of the hits from the track
		  recHit = itTraj->recHit();

		  // Check that the tracked hit is valid and belongs to the pixel tracker
		  if ((recHit->isValid() == true) &&
		      ((recHit->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) ||
		       (recHit->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)))
		    {
		      // Extract cluster information for the tracker hit
		      const SiPixelRecHit* recHitPixel = dynamic_cast<const SiPixelRecHit*>(recHit->hit());
		      SiPixelRecHit::ClusterRef const& PixelRecCluster = recHitPixel->cluster();
		      unsigned int FirstPixelRecClusterRow = PixelRecCluster->minPixelRow();
		      unsigned int FirstPixelRecClusterCol = PixelRecCluster->minPixelCol();
		      unsigned int LastPixelRecClusterRow = PixelRecCluster->maxPixelRow();
		      unsigned int LastPixelRecClusterCol = PixelRecCluster->maxPixelCol();
		      unsigned int IdPixelRecCluster = recHitPixel->geographicalId()();
		      
		      // Extract cluster information for the transient hit
		      SiPixelRecHit::ClusterRef const& PixelTransCluster = h->cluster();
		      unsigned int FirstPixelTransClusterRow = PixelTransCluster->minPixelRow();
		      unsigned int FirstPixelTransClusterCol = PixelTransCluster->minPixelCol();
		      unsigned int LastPixelTransClusterRow = PixelTransCluster->maxPixelRow();
		      unsigned int LastPixelTransClusterCol = PixelTransCluster->maxPixelCol();
		      unsigned int IdPixelTransCluster = h->geographicalId()();

		      // Does the TransiendRecHit not belong to the track?
		      if ((IdPixelTransCluster != IdPixelRecCluster) ||
			  (FirstPixelTransClusterRow != FirstPixelRecClusterRow) ||
			  (FirstPixelTransClusterCol != FirstPixelRecClusterCol) ||
			  (LastPixelTransClusterRow != LastPixelRecClusterRow) ||
			  (LastPixelTransClusterCol != LastPixelRecClusterCol))
			// 		      if (recHit->hit()->sharesInput(h,TrackingRecHit::all) == false)
			{
			  TrajectoryStateOnSurface TrajStateOnSurf = itTraj->forwardPredictedState();
			  
			  if (TrajStateOnSurf.isValid() == true)
			    {
			      // Propagation from the track hit to the theTransientRecHit surface
			      PropagatorWithMaterial thePropagator(alongMomentum, 0.1057, TrajStateOnSurf.magneticField(), MaxDPhi); // oppositeToMomentum, alongMomentum, anyDirection
			      TrajectoryStateOnSurface TrajStateOnNewSurf = thePropagator.propagate(TrajStateOnSurf, theTrackingGeometry->idToDet(h->geographicalId())->surface());
			      
			      double _PropLength = sqrt(powf((double)(TrajStateOnSurf.globalPosition().x() - theTransientHit->globalPosition().x()),2.)+
							powf((double)(TrajStateOnSurf.globalPosition().y() - theTransientHit->globalPosition().y()),2.)+
							powf((double)(TrajStateOnSurf.globalPosition().z() - theTransientHit->globalPosition().z()),2.));
			      			      
			      if (((PropLength < 0.0) || (_PropLength < PropLength)) && (TrajStateOnNewSurf.isValid() == true))
				{
				  PropLength = _PropLength;
				  TrackHitPosition = itTraj - Traj.begin();
				}
			    }
			}
		      else
			{
			  belongs2Track = true;
			  break;
			}
		    }
		}
		
		if ((TrackHitPosition != -1) && (belongs2Track == false))
		  {
		    KFUpdator updator;
		    TransientTrackingRecHit::ConstRecHitPointer recHit = Traj[TrackHitPosition].recHit();
		    TrajectoryStateOnSurface TrajStateOnSurf = updator.update(Traj[TrackHitPosition].forwardPredictedState(), *recHit);
 
		    // Propagation from the track hit to the theTransientRecHit surface
		    PropagatorWithMaterial thePropagator(alongMomentum, 0.1057, TrajStateOnSurf.magneticField(), MaxDPhi); // oppositeToMomentum, alongMomentum, anyDirection
		    TrajectoryStateOnSurface TrajStateOnNewSurf = thePropagator.propagate(TrajStateOnSurf, theTrackingGeometry->idToDet(h->geographicalId())->surface());

		    if (TrajStateOnNewSurf.isValid() == true)
		      {
			// Recalculate the hit position taking into account the direction of the magnetic field
			TransientTrackingRecHit::RecHitPointer theTransientHitAlongB = theTransientHit->clone(TrajStateOnNewSurf);

			double Distance = sqrt(powf(TrajStateOnNewSurf.localPosition().x() - theTransientHitAlongB->localPosition().x(),2.)+
					       powf(TrajStateOnNewSurf.localPosition().y() - theTransientHitAlongB->localPosition().y(),2.));
			
			if (Distance <= MaxDistance)
			  {			    
			    Chi2MeasurementEstimator Chi2Fw(MaxChi2);
			    MeasurementEstimator::HitReturnType Chi2FwResult = Chi2Fw.estimate(TrajStateOnNewSurf, *theTransientHitAlongB);
			    
			    if (Chi2FwResult.first == true)
			      {
				const Surface* surface = &(theTrackingGeometry->idToDet(h->geographicalId())->surface());
				GlobalPoint gp = surface->toGlobal(LocalPoint(0,0));
				int layerId = 0;
				if (h->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) layerId = ((PXBDetId)(h->geographicalId())).layer();
				if (h->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap) layerId = ((PXFDetId)(h->geographicalId())).disk() + 2*(((PXFDetId)(h->geographicalId())).side()-1);
				// side == 1 --> -Z; side == 2 --> +Z

				HitProperties HP;
				HP.GeogId = h->geographicalId()();
				HP.ParType = 0;
				HP.ProcType = 0;
				HP.IndexDet = j-rechitspixel->begin();
				HP.IndexHit = h-j->begin();
				HP.LayerId = layerId;
				HP.TrackHitPosition = TrackHitPosition;
				HP.PropLength = PropLength;
				HP.Distance = Distance;
				HP.Chi2 = Chi2FwResult.second;
				
				// ###### For simulated data ######
				if (SimData == true)
				  {
				    TransientTrackingRecHit::RecHitPointer theTH = theBuilder->build(&*h);				
				    // Extraction of the SimHits associated to a cluster with all their properties (e.g.
				    // particle type and process type)
				    // There might be several IDs because of the different strips/pixels in the cluster
				    std::vector<PSimHit> PSimTrackHits = hitAssociator->associateHit(*theTH);
				    if (PSimTrackHits.size() != 0)
				      {
					HP.ParType = abs(PSimTrackHits[0].particleType());
					HP.ProcType = PSimTrackHits[0].processType();
				      }
				  }
				
				unsigned int HIndex;
				for (HIndex = 0; HIndex < PixMissedHitsFw.size(); HIndex++)
				  if (PixMissedHitsFw[HIndex].GeogId == HP.GeogId) break;
				
				// Filling histograms if it's a new hit or if it has a better chi2
				if ((HIndex == PixMissedHitsFw.size()) || (HP.Chi2 < PixMissedHitsFw[HIndex].Chi2))
				  {
				    if (HIndex == PixMissedHitsFw.size()) PixMissedHitsFw.push_back(HP);
				    else PixMissedHitsFw[HIndex] = HP;
				  }
			      }
			  }
		      }
		    else if (PrintMsg == true) cout << "I was not able to re-propagate from track hit #" << TrackHitPosition << "\tPos.: " << TrajStateOnSurf.globalPosition() << " to TransRecHit #" << h-j->begin() << "\tID:" << theTransientHit->geographicalId()() << "\tPos.: " << theTransientHit->globalPosition() << endl;
		  }	
	      }
	    }
	    
	    // ##############################################################################
	    // # FILLING HISTOGRAMS WITH DATA ASSOCIATED TO THE TRACK-COMPATIBLE PIXEL HITS #
	    // ##############################################################################
	    unsigned int HIndex = 0;
	    while (HIndex < PixMissedHitsFw.size()) {
	      	      
	      SiPixelRecHitCollection::const_iterator j = rechitspixel->begin() + PixMissedHitsFw[HIndex].IndexDet;
	      edmNew::DetSet<SiPixelRecHit>::const_iterator h = j->begin() + PixMissedHitsFw[HIndex].IndexHit;
	      theTransientHit = theBuilder->build(&*h);

	      recHit = Traj[PixMissedHitsFw[HIndex].TrackHitPosition].recHit();
	      TrajStateOnSurf = updator.update(Traj[PixMissedHitsFw[HIndex].TrackHitPosition].forwardPredictedState(), *recHit);
	      
	      // Propagation from the track hit to the theTransientRecHit surface
	      PropagatorWithMaterial thePropagator(alongMomentum, 0.1057, TrajStateOnSurf.magneticField(), MaxDPhi); // oppositeToMomentum, alongMomentum, anyDirection
	      TrajStateOnNewSurf = thePropagator.propagate(TrajStateOnSurf, theTrackingGeometry->idToDet(h->geographicalId())->surface());
	      
	      // Recalculate the hit position taking into account the direction of the magnetic field
	      theTransientHitAlongB = theTransientHit->clone(TrajStateOnNewSurf);
	      
	      // ######################
	      // # Hit categorization #
	      // ######################
	      // ###### For simulated data ######
	      IsSimHit = false;
	      if (SimData == true)
		{
		  // Extraction of the IDs associated to a hit (= cluster)
		  // There might be several IDs because of the different strips/pixels in the cluster
		  std::vector<SimHitIdpr> SimTrackIDs = hitAssociator->associateHitId(*theTransientHitAlongB);
		  unsigned int jj;
		  for (unsigned int ii = 0; ii < SimTrackIDs.size(); ii++) {
		    for (jj = 0; jj < TrParIDs.size(); jj++) {
		      if (SimTrackIDs[ii].first == TrParIDs[jj])
			{
			  IsSimHit = true;
			  break;
			}
		    }
		    if (jj != TrParIDs.size()) break;
		  }
		}
	      
	      IsOther = false;
	      DoubleHit = false;
	      ToBeRecorded = false;
	      theGeomDetPix = dynamic_cast<const PixelGeomDetUnit*>(tracker->idToDet(h->geographicalId()));
	      RecPixTop = dynamic_cast<const RectangularPixelTopology*>(&(theGeomDetPix->specificTopology()));
	      SiPixelRecHit::ClusterRef const& PixelTransCluster = h->cluster();
	      if ((SimData == true) && (IsSimHit == false))
		{
		  H_category_pix->Fill("Non SimHit",1.0);
		  ToBeRecorded = true;
		}
	      else if (h->geographicalId() == recHit->geographicalId())
		{
		  H_category_pix->Fill("DoubleHit",1.0);
		  H_pix_dist_doublehit->Fill(sqrt(powf((double)(recHit->localPosition().x()-h->localPosition().x()),2.) + powf((double)(recHit->localPosition().y()-h->localPosition().y()),2.)));
		  if (theTransientHitAlongB->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel)
		    {
		      H_2DmapXY_PixDoubleHit_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_PixDoubleHit_B->Fill(theTransientHitAlongB->globalPosition().z(),
						     sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							  powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  else if (theTransientHitAlongB->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)
		    {
		      H_2DmapXY_PixDoubleHit_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_PixDoubleHit_F->Fill(theTransientHitAlongB->globalPosition().z(),
						     sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							  powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  H_pix_charge_doublehit_missed->Fill(h->cluster()->charge());
		  const SiPixelRecHit* recHitPix = dynamic_cast<const SiPixelRecHit*>(recHit->hit());
		  H_pix_charge_doublehit_trked->Fill(recHitPix->cluster()->charge());
		  DoubleHit = true;
		}
	      
	      // (PixBackSearch == true)
	      else if ((PixBackSearch == true) &&
		       ((h->geographicalId().subdetId() == recHit->geographicalId().subdetId()) &&
			((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) &&
			 (((PXBDetId)(h->geographicalId())).layer() == ((PXBDetId)(recHit->geographicalId())).layer()))))
		{
		  PixMissedHitsFw.erase(PixMissedHitsFw.begin() + HIndex);
		  HIndex--;
		}
	      else if ((PixBackSearch == true) &&
		       ((h->geographicalId().subdetId() == recHit->geographicalId().subdetId()) &&
			((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap) &&
			 (((PXFDetId)(h->geographicalId())).side() == ((PXFDetId)(recHit->geographicalId())).side()) &&
			 (((PXFDetId)(h->geographicalId())).disk() == ((PXFDetId)(recHit->geographicalId())).disk()))))
		{
		  PixMissedHitsFw.erase(PixMissedHitsFw.begin() + HIndex);
		  HIndex--;
		}
	      else if ((PixBackSearch == true) &&
		       ((((PixMissedHitsFw[HIndex].TrackHitPosition - 1) >= 0) && (Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->isValid() == true)) &&
			(h->geographicalId().subdetId() == Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->geographicalId().subdetId()) &&
			((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) &&
			 (((PXBDetId)(h->geographicalId())).layer() == ((PXBDetId)(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->geographicalId())).layer()))))
		{
		  PixMissedHitsFw.erase(PixMissedHitsFw.begin() + HIndex);
		  HIndex--;
		}
	      else if ((PixBackSearch == true) &&
		       ((((PixMissedHitsFw[HIndex].TrackHitPosition - 1) >= 0) && (Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->isValid() == true)) &&
			(h->geographicalId().subdetId() == Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->geographicalId().subdetId()) &&
			((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap) &&
			 (((PXFDetId)(h->geographicalId())).side() == ((PXFDetId)(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->geographicalId())).side()) &&
			 (((PXFDetId)(h->geographicalId())).disk() == ((PXFDetId)(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->geographicalId())).disk()))))
		{
		  PixMissedHitsFw.erase(PixMissedHitsFw.begin() + HIndex);
		  HIndex--;
		}
	      else if ((PixBackSearch == true) &&
		       ((((PixMissedHitsFw[HIndex].TrackHitPosition - 1) >= 0) && (Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->isValid() == false)) &&
			(((PixMissedHitsFw[HIndex].TrackHitPosition - 2) >= 0) && (Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->isValid() == true)) &&
			(h->geographicalId().subdetId() == Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->geographicalId().subdetId()) &&
			((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) &&
			 (((PXBDetId)(h->geographicalId())).layer() == ((PXBDetId)(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->geographicalId())).layer()))))
		{
		  PixMissedHitsFw.erase(PixMissedHitsFw.begin() + HIndex);
		  HIndex--;
		}
	      else if ((PixBackSearch == true) &&
		       ((((PixMissedHitsFw[HIndex].TrackHitPosition - 1) >= 0) && (Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->isValid() == false)) &&
			(((PixMissedHitsFw[HIndex].TrackHitPosition - 2) >= 0) && (Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->isValid() == true)) &&
			(h->geographicalId().subdetId() == Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->geographicalId().subdetId()) &&
			((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap) &&
			 (((PXFDetId)(h->geographicalId())).side() == ((PXFDetId)(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->geographicalId())).side()) &&
			 (((PXFDetId)(h->geographicalId())).disk() == ((PXFDetId)(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->geographicalId())).disk()))))
		{
		  PixMissedHitsFw.erase(PixMissedHitsFw.begin() + HIndex);
		  HIndex--;
		}
	      
	      // (PixBackSearch == false)
	      else if ((PixBackSearch == false) &&
		       ((h->geographicalId().subdetId() == recHit->geographicalId().subdetId()) &&
			((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) &&
			 (((PXBDetId)(h->geographicalId())).layer() == ((PXBDetId)(recHit->geographicalId())).layer()))))
		{
		  H_dist_pix_samelayer->Fill(fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(recHit->globalPosition().x(),2.) + powf(recHit->globalPosition().y(),2.))));
		  
		  if (fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(recHit->globalPosition().x(),2.) + powf(recHit->globalPosition().y(),2.))) >= RadialDistanceOverlap)
		    {
		      H_2DmapXY_PixOver_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_PixOver_B->Fill(theTransientHitAlongB->globalPosition().z(),
						sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						     powf(theTransientHitAlongB->globalPosition().y(),2.)));
		      H_category_pix->Fill("Overlap BPix",1.0);
		      ToBeRecorded = true;
		    }
		}
	      else if ((PixBackSearch == false) &&
		       ((h->geographicalId().subdetId() == recHit->geographicalId().subdetId()) &&
			((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap) &&
			 (((PXFDetId)(h->geographicalId())).side() == ((PXFDetId)(recHit->geographicalId())).side()) &&
			 (((PXFDetId)(h->geographicalId())).disk() == ((PXFDetId)(recHit->geographicalId())).disk()))))
		{
		  H_dist_pix_sameendcap->Fill(fabs(theTransientHitAlongB->globalPosition().z() - recHit->globalPosition().z()));

		  if (fabs(theTransientHitAlongB->globalPosition().z() - recHit->globalPosition().z()) >= LongitudinalDistanceOverlap)
		    {
		      H_2DmapXY_PixOver_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_PixOver_F->Fill(theTransientHitAlongB->globalPosition().z(),
						sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						     powf(theTransientHitAlongB->globalPosition().y(),2.)));
		      H_category_pix->Fill("Overlap FPix",1.0);
		      ToBeRecorded = true;
		    }
		}
	      else if ((PixBackSearch == false) &&
		       ((((PixMissedHitsFw[HIndex].TrackHitPosition - 1) >= 0) && (Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->isValid() == true)) &&
			(h->geographicalId().subdetId() == Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->geographicalId().subdetId()) &&
			((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) &&
			 (((PXBDetId)(h->geographicalId())).layer() == ((PXBDetId)(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->geographicalId())).layer()))))
		{
		  H_dist_pix_samelayer->Fill(fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->globalPosition().x(),2.) + powf(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->globalPosition().y(),2.))));

		  if (fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->globalPosition().x(),2.) + powf(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->globalPosition().y(),2.))) >= RadialDistanceOverlap)
		    {
		      H_2DmapXY_PixOver_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_PixOver_B->Fill(theTransientHitAlongB->globalPosition().z(),
						sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						     powf(theTransientHitAlongB->globalPosition().y(),2.)));
		      H_category_pix->Fill("Overlap BPix",1.0);				   
		      ToBeRecorded = true;
		    }
		}
	      else if ((PixBackSearch == false) &&
		       ((((PixMissedHitsFw[HIndex].TrackHitPosition - 1) >= 0) && (Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->isValid() == true)) &&
			(h->geographicalId().subdetId() == Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->geographicalId().subdetId()) &&
			((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap) &&
			 (((PXFDetId)(h->geographicalId())).side() == ((PXFDetId)(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->geographicalId())).side()) &&
			 (((PXFDetId)(h->geographicalId())).disk() == ((PXFDetId)(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->geographicalId())).disk()))))
		{
		  H_dist_pix_sameendcap->Fill(fabs(theTransientHitAlongB->globalPosition().z() - Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->globalPosition().z()));
		  
		  if (fabs(theTransientHitAlongB->globalPosition().z() - Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->globalPosition().z()) >= LongitudinalDistanceOverlap)
		    {
		      H_2DmapXY_PixOver_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_PixOver_F->Fill(theTransientHitAlongB->globalPosition().z(),
						sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						     powf(theTransientHitAlongB->globalPosition().y(),2.)));
		      H_category_pix->Fill("Overlap FPix",1.0);
		      ToBeRecorded = true;
		    }
		}
	      else if ((PixBackSearch == false) &&
		       ((((PixMissedHitsFw[HIndex].TrackHitPosition - 1) >= 0) && (Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->isValid() == false)) &&
			(((PixMissedHitsFw[HIndex].TrackHitPosition - 2) >= 0) && (Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->isValid() == true)) &&
			(h->geographicalId().subdetId() == Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->geographicalId().subdetId()) &&
			((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) &&
			 (((PXBDetId)(h->geographicalId())).layer() == ((PXBDetId)(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->geographicalId())).layer()))))
		{
		  H_dist_pix_samelayer->Fill(fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->globalPosition().x(),2.) + powf(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->globalPosition().y(),2.))));
		  
		  if (fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->globalPosition().x(),2.) + powf(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->globalPosition().y(),2.))) >= RadialDistanceOverlap)
		    {
		      H_2DmapXY_PixOver_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_PixOver_B->Fill(theTransientHitAlongB->globalPosition().z(),
						sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						     powf(theTransientHitAlongB->globalPosition().y(),2.)));
		      H_category_pix->Fill("Overlap BPix",1.0);				   
		      ToBeRecorded = true;
		    }
		}
	      else if ((PixBackSearch == false) &&
		       ((((PixMissedHitsFw[HIndex].TrackHitPosition - 1) >= 0) && (Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->isValid() == false)) &&
			(((PixMissedHitsFw[HIndex].TrackHitPosition - 2) >= 0) && (Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->isValid() == true)) &&
			(h->geographicalId().subdetId() == Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->geographicalId().subdetId()) &&
			((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap) &&
			 (((PXFDetId)(h->geographicalId())).side() == ((PXFDetId)(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->geographicalId())).side()) &&
			 (((PXFDetId)(h->geographicalId())).disk() == ((PXFDetId)(Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->geographicalId())).disk()))))
		{
		  H_dist_pix_sameendcap->Fill(fabs(theTransientHitAlongB->globalPosition().z() - Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->globalPosition().z()));

		  if (fabs(theTransientHitAlongB->globalPosition().z() - Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->globalPosition().z()) >= LongitudinalDistanceOverlap)
		    {
		      H_2DmapXY_PixOver_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_PixOver_F->Fill(theTransientHitAlongB->globalPosition().z(),
						sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						     powf(theTransientHitAlongB->globalPosition().y(),2.)));
		      H_category_pix->Fill("Overlap FPix",1.0);
		      ToBeRecorded = true;
		    }
		}
	      
	      else if (PixMissedHitsFw[HIndex].TrackHitPosition == 0)
		{
		  H_category_pix->Fill("Last",1.0);
		  ToBeRecorded = true;
		}
	      else if (RecPixTop->isItEdgePixel(PixelTransCluster->minPixelRow(), PixelTransCluster->minPixelCol()) ||
		       RecPixTop->isItEdgePixel(PixelTransCluster->maxPixelRow(), PixelTransCluster->maxPixelCol()) ||
		       RecPixTop->isItEdgePixel(PixelTransCluster->minPixelRow()-EdgeTollerance, PixelTransCluster->minPixelCol()) ||
		       RecPixTop->isItEdgePixel(PixelTransCluster->maxPixelRow()+EdgeTollerance, PixelTransCluster->maxPixelCol()) ||
		       RecPixTop->isItEdgePixel(PixelTransCluster->minPixelRow(), PixelTransCluster->minPixelCol()-EdgeTollerance) ||
		       RecPixTop->isItEdgePixel(PixelTransCluster->maxPixelRow(), PixelTransCluster->maxPixelCol()+EdgeTollerance) ||
		       RecPixTop->isItEdgePixel(PixelTransCluster->minPixelRow()-EdgeTollerance, PixelTransCluster->minPixelCol()-EdgeTollerance) ||
		       RecPixTop->isItEdgePixel(PixelTransCluster->maxPixelRow()+EdgeTollerance, PixelTransCluster->maxPixelCol()+EdgeTollerance))
		{
		  if (theTransientHitAlongB->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel)
		    {
		      H_2DmapXY_PixEdge_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_PixEdge_B->Fill(theTransientHitAlongB->globalPosition().z(),
						sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						     powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  else if (theTransientHitAlongB->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)
		    {
		      H_2DmapXY_PixEdge_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_PixEdge_F->Fill(theTransientHitAlongB->globalPosition().z(),
						sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						     powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  H_category_pix->Fill("Edge",1.0);
		  ToBeRecorded = true;
		}
	      else if (h->geographicalId().subdetId() != recHit->geographicalId().subdetId())
		{
		  if (theTransientHitAlongB->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel)
		    {
		      H_2DmapXY_PixTrans_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_PixTrans_B->Fill(theTransientHitAlongB->globalPosition().z(),
						 sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						      powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  else if (theTransientHitAlongB->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)
		    {
		      H_2DmapXY_PixTrans_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_PixTrans_F->Fill(theTransientHitAlongB->globalPosition().z(),
						 sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						      powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  H_category_pix->Fill("Aft.Trans.",1.0);
		  ToBeRecorded = true;
		}
	      else if ((((PixMissedHitsFw[HIndex].TrackHitPosition - 1) >= 0) && (Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->isValid() == true) &&
			(h->geographicalId().subdetId() != Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->geographicalId().subdetId())) ||
		       ((((PixMissedHitsFw[HIndex].TrackHitPosition - 1) >= 0) && (Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->isValid() == false) &&
			 (h->geographicalId().subdetId() == Traj[PixMissedHitsFw[HIndex].TrackHitPosition-1].recHit()->geographicalId().subdetId())) &&
			(((PixMissedHitsFw[HIndex].TrackHitPosition - 2) >= 0) &&
			 (h->geographicalId().subdetId() != Traj[PixMissedHitsFw[HIndex].TrackHitPosition-2].recHit()->geographicalId().subdetId()))))
		{
		  if (theTransientHitAlongB->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel)
		    {
		      H_2DmapXY_PixTrans_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_PixTrans_B->Fill(theTransientHitAlongB->globalPosition().z(),
						 sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						      powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  else if (theTransientHitAlongB->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)
		    {
		      H_2DmapXY_PixTrans_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_PixTrans_F->Fill(theTransientHitAlongB->globalPosition().z(),
						 sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						      powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  H_category_pix->Fill("Bef.Trans.",1.0);
		  ToBeRecorded = true;
		}
	      
	      else
		{
		  if (theTransientHitAlongB->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel)
		    {
		      H_2DmapXY_PixOther_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_PixOther_B->Fill(theTransientHitAlongB->globalPosition().z(),
						 sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						      powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  else if (theTransientHitAlongB->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)
		    {
		      H_2DmapXY_PixOther_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_PixOther_F->Fill(theTransientHitAlongB->globalPosition().z(),
						 sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						      powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  H_category_pix->Fill("Others",1.0);
		  ToBeRecorded = true;
		  IsOther = true; 
		}
	      
	      if ((ToBeRecorded == true) && (DoubleHit == false))
		{
		  MomentumCounter = true;
		  
		  H_chi2->Fill(PixMissedHitsFw[HIndex].Chi2);
		  H_chi2_length->Fill(PixMissedHitsFw[HIndex].Chi2, PixMissedHitsFw[HIndex].PropLength);
		  H_chi2_dist->Fill(PixMissedHitsFw[HIndex].Chi2, PixMissedHitsFw[HIndex].Distance);
		  H->Fill("Pixel hits",1.0); // Count the total number of missed hits: not including the "DoubleHits"
		  H_dist_pixe->Fill(PixMissedHitsFw[HIndex].Distance);
		  H_length_pixe->Fill(PixMissedHitsFw[HIndex].PropLength);
		  H_dist_length->Fill(PixMissedHitsFw[HIndex].Distance, PixMissedHitsFw[HIndex].PropLength);
		  H_normgeom->Fill(TrackMap.find(h->geographicalId().subdetId())->second, PixMissedHitsFw[HIndex].LayerId, 1.0);
		  H_hitpos->Fill(PixMissedHitsFw[HIndex].TrackHitPosition);
		  H_normhitpos->Fill((double)PixMissedHitsFw[HIndex].TrackHitPosition/(double)(Traj.size()-1));
		  
		  if (theTransientHitAlongB->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel)
		    {
		      H_2DmapXY_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_B->Fill(theTransientHitAlongB->globalPosition().z(),
					sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
					     powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  else if (theTransientHitAlongB->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)
		    {
		      H_2DmapXY_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
		      H_2DmapRZ_F->Fill(theTransientHitAlongB->globalPosition().z(),
					sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
					     powf(theTransientHitAlongB->globalPosition().y(),2.)));
		    }
		  
		  if (ClusterPlots == true)
		    {			       	
		      stringstream HistoName;
		      HistoName << "Pixel module ID: " << h->geographicalId()() << " #" << ++NPixelModules;
		      if (PrintMsg == true) cout << "\n" << HistoName.str().c_str() << endl;
		      H_general_pix_module = SubDirClus->make<TH2F>(HistoName.str().c_str(), HistoName.str().c_str(), 417, -0.5, 416.5, 161, -0.5, 160.5);
		      H_general_pix_module->SetXTitle("Columns");
		      H_general_pix_module->SetYTitle("Rows");
		      H_general_pix_module->SetZTitle("N.Clusters");
		    }
		  int CountHitsOnMod = 0;
		  int CountSimHitsOnMod = 0;
		  bool FirstMuon = false;
		  bool OtherMuon = true;
		  for (SiPixelRecHitCollection::const_iterator jIt = rechitspixel->begin(); jIt != rechitspixel->end(); jIt++) { // Loop on detectors
		    for (edmNew::DetSet<SiPixelRecHit>::const_iterator hIt = jIt->begin(); hIt != jIt->end(); hIt++) { // Loop on hits on detector
		      if (h->geographicalId()() == hIt->geographicalId()())
			{
			  CountHitsOnMod++;
			  
			  // ###### For simulated data ######
			  if (SimData == true)
			    {
			      TransientTrackingRecHit::RecHitPointer theTH = theBuilder->build(&*hIt);						
			      // Extraction of the SimHits associated to a cluster with all their properties (e.g.
			      // particle type and process type)
			      // There might be several IDs because of the different strips/pixels in the cluster
			      std::vector<PSimHit> PSimTrackHits = hitAssociator->associateHit(*theTH);
			      if (PrintMsg == true) cout << "New pixel cluster" << endl;
			      if ((PixMissedHitsFw[HIndex].IndexDet == jIt-rechitspixel->begin()) && (PixMissedHitsFw[HIndex].IndexHit == hIt-jIt->begin()) &&
				  (abs(PSimTrackHits[0].particleType()) == ParticleType) && (PSimTrackHits[0].processType() == ProcessType)) FirstMuon = true;
			      else if (((PixMissedHitsFw[HIndex].IndexDet != jIt-rechitspixel->begin()) || (PixMissedHitsFw[HIndex].IndexHit == hIt-jIt->begin())) &&
				       ((abs(PSimTrackHits[0].particleType()) != ParticleType) || (PSimTrackHits[0].processType() != ProcessType))) OtherMuon = false;

			      if (PrintMsg == true)
				for (unsigned int ii = 0; ii < PSimTrackHits.size(); ii++)
				  cout << "Pixel cluster\tparticle type: " << PSimTrackHits[ii].particleType() << "\tprocess type: " << PSimTrackHits[ii].processType() << endl;
			      
			      // Extraction of the IDs associated to a hit (= cluster)
			      // There might be several IDs because of the different strips/pixels in the cluster
			      std::vector<SimHitIdpr> SimTrackIDs = hitAssociator->associateHitId(*theTH);
			      unsigned int jj;
			      for (unsigned int ii = 0; ii < SimTrackIDs.size(); ii++) {
				for (jj = 0; jj < TrParIDs.size(); jj++) {
				  if (SimTrackIDs[ii].first == TrParIDs[jj])
				    {
				      CountSimHitsOnMod++;
					
				      if (ClusterPlots == true)
					for (unsigned int k = 0; k < hIt->cluster()->pixels().size(); k++)
					  // H_general_pix_module->Fill(hIt->cluster()->pixels()[k].y & 511, hIt->cluster()->pixels()[k].x, hIt->cluster()->pixels()[k].adc);
					  H_general_pix_module->Fill(hIt->cluster()->pixels()[k].y & 511, hIt->cluster()->pixels()[k].x, CountSimHitsOnMod);
				      // Only for the Y coordinate:
				      // 9 bits for Col information
				      // 1 bit for larger clusters than 9x33
				      // 6 bits for quality information
				      
				      break;
				    }
				}
				if (jj != TrParIDs.size()) break;
			      }
			    }
			}
		    }
		  }
		  H_rechits_on_mod_pixe->Fill(CountHitsOnMod);
		  H_simhits_on_mod_pixe->Fill(((double)CountSimHitsOnMod)/((double)CountHitsOnMod)*100.);
		  
		  // Take into account multiple track-compatible hits
		  if ((SimData == true) && (CountHitsOnMod > 0) && (FirstMuon == true) && (OtherMuon == true))
		    H_category_multy_pix->Fill("AllMu",1.0);
		  else if ((SimData == true) && (CountHitsOnMod > 0) && (FirstMuon == true) && (OtherMuon == false))
		    H_category_multy_pix->Fill("1Mu @least 1NoMu",1.0);
		  else if ((SimData == true) && (CountHitsOnMod > 0) && (FirstMuon == false) && (OtherMuon == false))
		    H_category_multy_pix->Fill("NoMu",1.0);
		  else if ((SimData == true) && (CountHitsOnMod > 0) && (FirstMuon == false) && (OtherMuon == true))
		    H_category_multy_pix->Fill("1NoMu o/AllMu",1.0);
		  
		  // Take into account single track-compatible hits
		  if (CountHitsOnMod > 0)
		    {
		      std::vector<PSimHit> PSimTrackHits;
		      if (SimData == true) PSimTrackHits = hitAssociator->associateHit(*theTransientHit);

		      if ((SimData == false) ||
			  ((PSimTrackHits.size() != 0) && (abs(PSimTrackHits[0].particleType()) == ParticleType) && (PSimTrackHits[0].processType() == ProcessType)))
			{
			  H_pix_charge->Fill(h->cluster()->charge());

			  if (h->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel)
			    {
			      PixelBarrelName BPixName(h->geographicalId());
			      if (BPixName.moduleType() == PixelBarrelName::v2x8)
				{
				  PixelHit2x8Counter++;
				  for (unsigned int k = 0; k < h->cluster()->pixels().size(); k++)
				    H_2x8_pix_module->Fill(h->cluster()->pixels()[k].y & 511, h->cluster()->pixels()[k].x, h->cluster()->pixels()[k].adc);
				}
			      else if (BPixName.moduleType() == PixelBarrelName::v1x8)
				{
				  PixelHit1x8Counter++;
				  for (unsigned int k = 0; k < h->cluster()->pixels().size(); k++)
				    H_1x8_pix_module->Fill(h->cluster()->pixels()[k].y & 511, h->cluster()->pixels()[k].x, h->cluster()->pixels()[k].adc);
				}
			    }
			  
			  else if (h->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)
			    {
			      PixelEndcapName FPixName(h->geographicalId());
			      if (FPixName.moduleType() == PixelEndcapName::v1x2)
				{
				  PixelHit1x2Counter++;
				  for (unsigned int k = 0; k < h->cluster()->pixels().size(); k++)
				    H_1x2_pix_module->Fill(h->cluster()->pixels()[k].y & 511, h->cluster()->pixels()[k].x, h->cluster()->pixels()[k].adc);
				}
			      else if (FPixName.moduleType() == PixelEndcapName::v1x5)
				{
				  PixelHit1x5Counter++;
				  for (unsigned int k = 0; k < h->cluster()->pixels().size(); k++)
				    H_1x5_pix_module->Fill(h->cluster()->pixels()[k].y & 511, h->cluster()->pixels()[k].x, h->cluster()->pixels()[k].adc);
				}
			      else if (FPixName.moduleType() == PixelEndcapName::v2x3)
				{
				  PixelHit2x3Counter++;
				  for (unsigned int k = 0; k < h->cluster()->pixels().size(); k++)
				    H_2x3_pix_module->Fill(h->cluster()->pixels()[k].y & 511, h->cluster()->pixels()[k].x, h->cluster()->pixels()[k].adc);
				}
			      else if (FPixName.moduleType() == PixelEndcapName::v2x4)
				{
				  PixelHit2x4Counter++;
				  for (unsigned int k = 0; k < h->cluster()->pixels().size(); k++)
				    H_2x4_pix_module->Fill(h->cluster()->pixels()[k].y & 511, h->cluster()->pixels()[k].x, h->cluster()->pixels()[k].adc);
				}
			      else if (FPixName.moduleType() == PixelEndcapName::v2x5)
				{
				  PixelHit2x5Counter++;
				  for (unsigned int k = 0; k < h->cluster()->pixels().size(); k++)
				    H_2x5_pix_module->Fill(h->cluster()->pixels()[k].y & 511, h->cluster()->pixels()[k].x, h->cluster()->pixels()[k].adc);
				}
			    }
			}
		    }
		  
		  if (PrintMsg == true)
		    {
		      cout << "\nI found a PIXEL hit that doesn't belong to the track in the event: " <<  iEvent.id() << endl;
		      cout << "\tHit distance from the track: " << PixMissedHitsFw[HIndex].Distance << " (<=" << MaxDistance << " cm)" << "\tPropagation length: " << PixMissedHitsFw[HIndex].PropLength << endl;
		      cout << "\tChi2 compatibility with TransientRecHit: " << PixMissedHitsFw[HIndex].Chi2 << " (<=" << MaxChi2 << ")" << endl;
		      cout << "\tPropagation from track hit #" << PixMissedHitsFw[HIndex].TrackHitPosition << "\tTSOS pos.: " << TrajStateOnSurf.globalPosition() << "\tHit pos.: " << recHit->globalPosition() << endl;
		      cout << "\t\tSubDetId: " << recHit->geographicalId().subdetId() << "\tID: " << recHit->geographicalId()() << endl;
		      cout << "\tResult of propagation: " << TrajStateOnNewSurf.globalPosition() << endl;
		      cout << "\tTransRecHit #" << h-j->begin() << "\tID: " << theTransientHitAlongB->geographicalId()() << "\tPos.: " << theTransientHitAlongB->globalPosition() << endl;
		      cout << "\t\tSubDetId: " << h->geographicalId().subdetId() << "\tLayerID: " << PixMissedHitsFw[HIndex].LayerId;
		      cout << "\tMin row: " << h->cluster()->minPixelRow() << "\tMin col: " << h->cluster()->minPixelCol() << "\tMax row: " << h->cluster()->maxPixelRow() << "\tMax col: " << h->cluster()->maxPixelCol() << endl;
		    }
		}
	      HIndex++;
	    }
	    
	    // ###################
	    // # For back search #
	    // ###################
	    if (PixBackSearch == true)
	      {
		// ##########################################
		// # SEARCH FOR TRACK-COMPATIBLE PIXEL HITS #
		// ##########################################
		unsigned int itHP;
		for (SiPixelRecHitCollection::const_iterator j = rechitspixel->begin(); j != rechitspixel->end(); j++) { // Loop on detectors
		  for (edmNew::DetSet<SiPixelRecHit>::const_iterator h = j->begin(); h != j->end(); h++) { // Loop on hits on detector
			   
		    for (itHP = 0; itHP < PixMissedHitsFw.size(); itHP++)
		      if ((PixMissedHitsFw[itHP].IndexDet == j-rechitspixel->begin()) && (PixMissedHitsFw[itHP].IndexHit == h-j->begin())) break;
		    if (itHP != PixMissedHitsFw.size()) continue;

		    // The builder on the recHit add global information like: global position, magnetic field
		    theTransientHit = theBuilder->build(&*h);
			   
		    PropLength = -1.0;
		    belongs2Track = false;
		    TrackHitPosition = -1;
			   
		    // Loop over the hits of the track
		    for (std::vector<TrajectoryMeasurement>::const_iterator itTraj = Traj.begin(); itTraj != Traj.end(); itTraj++) { // Loop on hits on detector
			     
		      // Extraction of the hits from the track
		      recHit = itTraj->recHit();
		      
		      // Check that the tracked hit is valid and belongs to the pixel tracker
		      if ((recHit->isValid() == true) &&
			  ((recHit->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) ||
			   (recHit->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)))
			{
			  // Extract cluster information for the tracker hit
			  const SiPixelRecHit* recHitPixel = dynamic_cast<const SiPixelRecHit*>(recHit->hit());
			  SiPixelRecHit::ClusterRef const& PixelRecCluster = recHitPixel->cluster();
			  unsigned int FirstPixelRecClusterRow = PixelRecCluster->minPixelRow();
			  unsigned int FirstPixelRecClusterCol = PixelRecCluster->minPixelCol();
			  unsigned int LastPixelRecClusterRow = PixelRecCluster->maxPixelRow();
			  unsigned int LastPixelRecClusterCol = PixelRecCluster->maxPixelCol();
			  unsigned int IdPixelRecCluster = recHitPixel->geographicalId()();
				 
			  // Extract cluster information for the transient hit
			  SiPixelRecHit::ClusterRef const& PixelTransCluster = h->cluster();
			  unsigned int FirstPixelTransClusterRow = PixelTransCluster->minPixelRow();
			  unsigned int FirstPixelTransClusterCol = PixelTransCluster->minPixelCol();
			  unsigned int LastPixelTransClusterRow = PixelTransCluster->maxPixelRow();
			  unsigned int LastPixelTransClusterCol = PixelTransCluster->maxPixelCol();
			  unsigned int IdPixelTransCluster = h->geographicalId()();
				 
			  // Does the TransiendRecHit not belong to the track?
			  if ((IdPixelTransCluster != IdPixelRecCluster) ||
			      (FirstPixelTransClusterRow != FirstPixelRecClusterRow) ||
			      (FirstPixelTransClusterCol != FirstPixelRecClusterCol) ||
			      (LastPixelTransClusterRow != LastPixelRecClusterRow) ||
			      (LastPixelTransClusterCol != LastPixelRecClusterCol))
			    // 			    if (recHit->hit()->sharesInput(h,TrackingRecHit::all) == false)
			    {
			      TrajectoryStateOnSurface TrajStateOnSurf = itTraj->backwardPredictedState();
				     
			      if (TrajStateOnSurf.isValid() == true)
				{
				  // Propagation from the track hit to the theTransientRecHit surface
				  PropagatorWithMaterial thePropagator(oppositeToMomentum, 0.1057, TrajStateOnSurf.magneticField(), MaxDPhi); // oppositeToMomentum, alongMomentum, anyDirection
				  TrajectoryStateOnSurface TrajStateOnNewSurf = thePropagator.propagate(TrajStateOnSurf, theTrackingGeometry->idToDet(h->geographicalId())->surface());
					 
				  double _PropLength = sqrt(powf((double)(TrajStateOnSurf.globalPosition().x() - theTransientHit->globalPosition().x()),2.)+
							    powf((double)(TrajStateOnSurf.globalPosition().y() - theTransientHit->globalPosition().y()),2.)+
							    powf((double)(TrajStateOnSurf.globalPosition().z() - theTransientHit->globalPosition().z()),2.));
					 					 
				  if (((PropLength < 0.0) || (_PropLength < PropLength)) && (TrajStateOnNewSurf.isValid() == true) &&
				      (((theTransientHit->globalPosition().x() > 0.) && (TrajStateOnSurf.globalPosition().x() > (0. - BPixTollerance))) ||
				       ((theTransientHit->globalPosition().x() < 0.) && (TrajStateOnSurf.globalPosition().x() < (0. + BPixTollerance)))) &&
				      (((theTransientHit->globalPosition().y() > 0.) && (TrajStateOnSurf.globalPosition().y() > (0. - BPixTollerance))) ||
				       ((theTransientHit->globalPosition().y() < 0.) && (TrajStateOnSurf.globalPosition().y() < (0. + BPixTollerance)))))
				    {
				      PropLength = _PropLength;
				      TrackHitPosition = itTraj - Traj.begin();
				    }
				}
			    }
			  else
			    {
			      belongs2Track = true;
			      break;
			    }
			}
		    }
			   
		    if ((TrackHitPosition != -1) && (belongs2Track == false))
		      {
			KFUpdator updator;
			TransientTrackingRecHit::ConstRecHitPointer recHit = Traj[TrackHitPosition].recHit();
			TrajectoryStateOnSurface TrajStateOnSurf = updator.update(Traj[TrackHitPosition].backwardPredictedState(), *recHit);
			       
			// Propagation from the track hit to the theTransientRecHit surface
			PropagatorWithMaterial thePropagator(oppositeToMomentum, 0.1057, TrajStateOnSurf.magneticField(), MaxDPhi); // oppositeToMomentum, alongMomentum, anyDirection
			TrajectoryStateOnSurface TrajStateOnNewSurf = thePropagator.propagate(TrajStateOnSurf, theTrackingGeometry->idToDet(h->geographicalId())->surface());
			       
			if (TrajStateOnNewSurf.isValid() == true)
			  {
			    // Recalculate the hit position taking into account the direction of the magnetic field
			    TransientTrackingRecHit::RecHitPointer theTransientHitAlongB = theTransientHit->clone(TrajStateOnNewSurf);
				   
			    double Distance = sqrt(powf(TrajStateOnNewSurf.localPosition().x() - theTransientHitAlongB->localPosition().x(),2.)+
						   powf(TrajStateOnNewSurf.localPosition().y() - theTransientHitAlongB->localPosition().y(),2.));

			    if (Distance <= MaxDistance)
			      {		
				Chi2MeasurementEstimator Chi2Fw(MaxChi2);
				MeasurementEstimator::HitReturnType Chi2FwResult = Chi2Fw.estimate(TrajStateOnNewSurf, *theTransientHitAlongB);
				       				       
				if (Chi2FwResult.first == true)
				  {
				    const Surface* surface = &(theTrackingGeometry->idToDet(h->geographicalId())->surface());
				    GlobalPoint gp = surface->toGlobal(LocalPoint(0,0));
				    int layerId = 0;
				    if (h->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) layerId = ((PXBDetId)(h->geographicalId())).layer();
				    if (h->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap) layerId = ((PXFDetId)(h->geographicalId())).disk() + 2*(((PXFDetId)(h->geographicalId())).side()-1);
				    // side == 1 --> -Z; side == 2 --> +Z

				    HitProperties HP;
				    HP.GeogId = h->geographicalId()();
				    HP.ParType = 0;
				    HP.ProcType = 0;
				    HP.IndexDet = j-rechitspixel->begin();
				    HP.IndexHit = h-j->begin();
				    HP.LayerId = layerId;
				    HP.TrackHitPosition = TrackHitPosition;
				    HP.PropLength = PropLength;
				    HP.Distance = Distance;
				    HP.Chi2 = Chi2FwResult.second;
				    
				    // ###### For simulated data ######
				    if (SimData == true)
				      {
					TransientTrackingRecHit::RecHitPointer theTH = theBuilder->build(&*h);				
					// Extraction of the SimHits associated to a cluster with all their properties (e.g.
					// particle type and process type)
					// There might be several IDs because of the different strips/pixels in the cluster
					std::vector<PSimHit> PSimTrackHits = hitAssociator->associateHit(*theTH);
					if (PSimTrackHits.size() != 0)
					  {
					    HP.ParType = abs(PSimTrackHits[0].particleType());
					    HP.ProcType = PSimTrackHits[0].processType();
					  }
				      }

				    unsigned int HIndex;
				    for (HIndex = 0; HIndex < PixMissedHitsBw.size(); HIndex++)
				      if (PixMissedHitsBw[HIndex].GeogId == HP.GeogId) break;
				    
				    // Filling histograms if it's a new hit or if it has a better chi2
				    if ((HIndex == PixMissedHitsBw.size()) || (HP.Chi2 < PixMissedHitsBw[HIndex].Chi2))
				      {
					if (HIndex == PixMissedHitsBw.size()) PixMissedHitsBw.push_back(HP);
					else PixMissedHitsBw[HIndex] = HP;
				      }
				  }			
			      }
			  }
			else if (PrintMsg == true) cout << "I was not able to re-propagate from track hit #" << TrackHitPosition << "\tPos.: " << TrajStateOnSurf.globalPosition() << " to TransRecHit #" << h-j->begin() << "\tID:" << theTransientHit->geographicalId()() << "\tPos.: " << theTransientHit->globalPosition() << endl;
		      }	
		  }
		}

		// ##############################################################################
		// # FILLING HISTOGRAMS WITH DATA ASSOCIATED TO THE TRACK-COMPATIBLE PIXEL HITS #
		// ##############################################################################
		unsigned int jj;
		for (unsigned int HIndex = 0; HIndex < PixMissedHitsBw.size(); HIndex++) {

		  SiPixelRecHitCollection::const_iterator j = rechitspixel->begin() + PixMissedHitsBw[HIndex].IndexDet;
		  edmNew::DetSet<SiPixelRecHit>::const_iterator h = j->begin() + PixMissedHitsBw[HIndex].IndexHit;
		  theTransientHit = theBuilder->build(&*h);
			   
		  recHit = Traj[PixMissedHitsBw[HIndex].TrackHitPosition].recHit();
		  TrajStateOnSurf = updator.update(Traj[PixMissedHitsBw[HIndex].TrackHitPosition].backwardPredictedState(), *recHit);
			   
		  // Propagation from the track hit to the theTransientRecHit surface
		  PropagatorWithMaterial thePropagator(oppositeToMomentum, 0.1057, TrajStateOnSurf.magneticField(), MaxDPhi); // oppositeToMomentum, alongMomentum, anyDirection
		  TrajectoryStateOnSurface TrajStateOnNewSurf = thePropagator.propagate(TrajStateOnSurf, theTrackingGeometry->idToDet(h->geographicalId())->surface());
			   
		  // Recalculate the hit position taking into account the direction of the magnetic field
		  theTransientHitAlongB = theTransientHit->clone(TrajStateOnNewSurf);
			
		  // ######################  
		  // # Hit categorization #
		  // ######################
		  // ###### For simulated data ######
		  IsSimHit = false;
		  if (SimData == true)
		    {
		      // Extraction of the IDs associated to a hit (= cluster)
		      // There might be several IDs because of the different strips/pixels in the cluster
		      std::vector<SimHitIdpr> SimTrackIDs = hitAssociator->associateHitId(*theTransientHitAlongB);
		      for (unsigned int ii = 0; ii < SimTrackIDs.size(); ii++) {
			for (jj = 0; jj < TrParIDs.size(); jj++) {
			  if (SimTrackIDs[ii].first == TrParIDs[jj])
			    {
			      IsSimHit = true;
			      break;
			    }
			}
			if (jj != TrParIDs.size()) break;
		      }
		    }
	   
		  IsOther = false;
		  DoubleHit = false;
		  // First check that's an "Overlap"
		  if (((h->geographicalId().subdetId() == recHit->geographicalId().subdetId()) &&
		       ((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) &&
			(((PXBDetId)(h->geographicalId())).layer() == ((PXBDetId)(recHit->geographicalId())).layer()))) ||
			       
		      ((h->geographicalId().subdetId() == recHit->geographicalId().subdetId()) &&
		       ((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap) &&
			(((PXFDetId)(h->geographicalId())).side() == ((PXFDetId)(recHit->geographicalId())).side()) &&
			(((PXFDetId)(h->geographicalId())).disk() == ((PXFDetId)(recHit->geographicalId())).disk()))) ||
			       
		      ((((PixMissedHitsBw[HIndex].TrackHitPosition + 1) < (signed)Traj.size()) && (Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->isValid() == true)) &&
		       (h->geographicalId().subdetId() == Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->geographicalId().subdetId()) &&
		       ((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) &&
			(((PXBDetId)(h->geographicalId())).layer() == ((PXBDetId)(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->geographicalId())).layer()))) ||
			       
		      ((((PixMissedHitsBw[HIndex].TrackHitPosition + 1) < (signed)Traj.size()) && (Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->isValid() == true)) &&
		       (h->geographicalId().subdetId() == Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->geographicalId().subdetId()) &&
		       ((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap) &&
			(((PXFDetId)(h->geographicalId())).side() == ((PXFDetId)(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->geographicalId())).side()) &&
			(((PXFDetId)(h->geographicalId())).disk() == ((PXFDetId)(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->geographicalId())).disk()))) ||
			       
		      ((((PixMissedHitsBw[HIndex].TrackHitPosition + 1) < (signed)Traj.size()) && (Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->isValid() == false)) &&
		       (((PixMissedHitsBw[HIndex].TrackHitPosition + 2) < (signed)Traj.size()) && (Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->isValid() == true)) &&
		       (h->geographicalId().subdetId() == Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->geographicalId().subdetId()) &&
		       ((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) &&
			(((PXBDetId)(h->geographicalId())).layer() == ((PXBDetId)(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->geographicalId())).layer()))) ||
			       
		      ((((PixMissedHitsBw[HIndex].TrackHitPosition + 1) < (signed)Traj.size()) && (Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->isValid() == false)) &&
		       (((PixMissedHitsBw[HIndex].TrackHitPosition + 2) < (signed)Traj.size()) && (Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->isValid() == true)) &&
		       (h->geographicalId().subdetId() == Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->geographicalId().subdetId()) &&
		       ((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap) &&
			(((PXFDetId)(h->geographicalId())).side() == ((PXFDetId)(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->geographicalId())).side()) &&
			(((PXFDetId)(h->geographicalId())).disk() == ((PXFDetId)(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->geographicalId())).disk()))))
		    {
		      if ((SimData == true) && (IsSimHit == false))
			H_category_pix->Fill("Non SimHit",1.0);
		      else if (h->geographicalId() == recHit->geographicalId())
			{
			  H_category_pix->Fill("DoubleHit",1.0);
			  H_pix_dist_doublehit->Fill(sqrt(powf((double)(recHit->localPosition().x()-h->localPosition().x()),2.) + powf((double)(recHit->localPosition().y()-h->localPosition().y()),2.)));
			  if (theTransientHitAlongB->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel)
			    {
			      H_2DmapXY_PixDoubleHit_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			      H_2DmapRZ_PixDoubleHit_B->Fill(theTransientHitAlongB->globalPosition().z(),
							     sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
								  powf(theTransientHitAlongB->globalPosition().y(),2.)));
			    }
			  else if (theTransientHitAlongB->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)
			    {
			      H_2DmapXY_PixDoubleHit_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			      H_2DmapRZ_PixDoubleHit_F->Fill(theTransientHitAlongB->globalPosition().z(),
							     sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
								  powf(theTransientHitAlongB->globalPosition().y(),2.)));
			    }
			  H_pix_charge_doublehit_missed->Fill(h->cluster()->charge());
			  const SiPixelRecHit* recHitPix = dynamic_cast<const SiPixelRecHit*>(recHit->hit());
			  H_pix_charge_doublehit_trked->Fill(recHitPix->cluster()->charge());
			  DoubleHit = true;
			}				    					
		      else if ((h->geographicalId().subdetId() == recHit->geographicalId().subdetId()) &&
			       ((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) &&
				(((PXBDetId)(h->geographicalId())).layer() == ((PXBDetId)(recHit->geographicalId())).layer())))
			{
			  H_dist_pix_samelayer->Fill(fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(recHit->globalPosition().x(),2.) + powf(recHit->globalPosition().y(),2.))));
			  
			  if (fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(recHit->globalPosition().x(),2.) + powf(recHit->globalPosition().y(),2.))) >= RadialDistanceOverlap)
			    {
			      H_2DmapXY_PixOver_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			      H_2DmapRZ_PixOver_B->Fill(theTransientHitAlongB->globalPosition().z(),
							sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							     powf(theTransientHitAlongB->globalPosition().y(),2.)));
			      H_category_pix->Fill("Overlap BPix",1.0);
			    }
			}
		      else if ((h->geographicalId().subdetId() == recHit->geographicalId().subdetId()) &&
			       ((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap) &&
				(((PXFDetId)(h->geographicalId())).side() == ((PXFDetId)(recHit->geographicalId())).side()) &&
				(((PXFDetId)(h->geographicalId())).disk() == ((PXFDetId)(recHit->geographicalId())).disk())))
			{
			  H_dist_pix_sameendcap->Fill(fabs(theTransientHitAlongB->globalPosition().z() - recHit->globalPosition().z()));

			  if (fabs(theTransientHitAlongB->globalPosition().z() - recHit->globalPosition().z()) >= LongitudinalDistanceOverlap)
			    {
			      H_2DmapXY_PixOver_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			      H_2DmapRZ_PixOver_F->Fill(theTransientHitAlongB->globalPosition().z(),
							sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							     powf(theTransientHitAlongB->globalPosition().y(),2.)));					
			      H_category_pix->Fill("Overlap FPix",1.0);
			    }
			}
		      else if ((((PixMissedHitsBw[HIndex].TrackHitPosition + 1) < (signed)Traj.size()) && (Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->isValid() == true)) &&
			       (h->geographicalId().subdetId() == Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->geographicalId().subdetId()) &&
			       ((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) &&
				(((PXBDetId)(h->geographicalId())).layer() == ((PXBDetId)(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->geographicalId())).layer())))
			{
			  H_dist_pix_samelayer->Fill(fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->globalPosition().x(),2.) + powf(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->globalPosition().y(),2.))));
			  
			  if (fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->globalPosition().x(),2.) + powf(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->globalPosition().y(),2.))) >= RadialDistanceOverlap)
			    {
			      H_2DmapXY_PixOver_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			      H_2DmapRZ_PixOver_B->Fill(theTransientHitAlongB->globalPosition().z(),
							sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							     powf(theTransientHitAlongB->globalPosition().y(),2.)));
			      H_category_pix->Fill("Overlap BPix",1.0);				   
			    }
			}
		      else if ((((PixMissedHitsBw[HIndex].TrackHitPosition + 1) < (signed)Traj.size()) && (Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->isValid() == true)) &&
			       (h->geographicalId().subdetId() == Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->geographicalId().subdetId()) &&
			       ((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap) &&
				(((PXFDetId)(h->geographicalId())).side() == ((PXFDetId)(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->geographicalId())).side()) &&
				(((PXFDetId)(h->geographicalId())).disk() == ((PXFDetId)(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->geographicalId())).disk())))
			{
			  H_dist_pix_sameendcap->Fill(fabs(theTransientHitAlongB->globalPosition().z() - Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->globalPosition().z()));

			  if (fabs(theTransientHitAlongB->globalPosition().z() - Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->globalPosition().z()) >= LongitudinalDistanceOverlap)
			    {
			      H_2DmapXY_PixOver_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			      H_2DmapRZ_PixOver_F->Fill(theTransientHitAlongB->globalPosition().z(),
							sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							     powf(theTransientHitAlongB->globalPosition().y(),2.)));					
			      H_category_pix->Fill("Overlap FPix",1.0);
			    }
			}
		      else if ((((PixMissedHitsBw[HIndex].TrackHitPosition + 1) < (signed)Traj.size()) && (Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->isValid() == false)) &&
			       (((PixMissedHitsBw[HIndex].TrackHitPosition + 2) < (signed)Traj.size()) && (Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->isValid() == true)) &&
			       (h->geographicalId().subdetId() == Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->geographicalId().subdetId()) &&
			       ((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) &&
				(((PXBDetId)(h->geographicalId())).layer() == ((PXBDetId)(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->geographicalId())).layer())))
			{
			  H_dist_pix_samelayer->Fill(fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->globalPosition().x(),2.) + powf(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->globalPosition().y(),2.))));
			  
			  if (fabs(sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) + powf(theTransientHitAlongB->globalPosition().y(),2.)) - sqrt(powf(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->globalPosition().x(),2.) + powf(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->globalPosition().y(),2.))) >= RadialDistanceOverlap)
			    {
			      H_2DmapXY_PixOver_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			      H_2DmapRZ_PixOver_B->Fill(theTransientHitAlongB->globalPosition().z(),
							sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							     powf(theTransientHitAlongB->globalPosition().y(),2.)));
			      H_category_pix->Fill("Overlap BPix",1.0);				   
			    }
			}
		      else if ((((PixMissedHitsBw[HIndex].TrackHitPosition + 1) < (signed)Traj.size()) && (Traj[PixMissedHitsBw[HIndex].TrackHitPosition+1].recHit()->isValid() == false)) &&
			       (((PixMissedHitsBw[HIndex].TrackHitPosition + 2) < (signed)Traj.size()) && (Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->isValid() == true)) &&
			       (h->geographicalId().subdetId() == Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->geographicalId().subdetId()) &&
			       ((h->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap) &&
				(((PXFDetId)(h->geographicalId())).side() == ((PXFDetId)(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->geographicalId())).side()) &&
				(((PXFDetId)(h->geographicalId())).disk() == ((PXFDetId)(Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->geographicalId())).disk())))
			{
			  H_dist_pix_sameendcap->Fill(fabs(theTransientHitAlongB->globalPosition().z() - Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->globalPosition().z()));

			  if (fabs(theTransientHitAlongB->globalPosition().z() - Traj[PixMissedHitsBw[HIndex].TrackHitPosition+2].recHit()->globalPosition().z()) >= LongitudinalDistanceOverlap)
			    {
			      H_2DmapXY_PixOver_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			      H_2DmapRZ_PixOver_F->Fill(theTransientHitAlongB->globalPosition().z(),
							sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
							     powf(theTransientHitAlongB->globalPosition().y(),2.)));					
			      H_category_pix->Fill("Overlap FPix",1.0);
			    }
			}
		      else IsOther = true;
		    }
		  else IsOther = true;			   

		  if ((IsOther == false) && (DoubleHit == false))
		    {
		      MomentumCounter = true;
			       
		      H_chi2->Fill(PixMissedHitsBw[HIndex].Chi2);
		      H_chi2_length->Fill(PixMissedHitsBw[HIndex].Chi2, PixMissedHitsBw[HIndex].PropLength);
		      H_chi2_dist->Fill(PixMissedHitsBw[HIndex].Chi2, PixMissedHitsBw[HIndex].Distance);
		      H->Fill("Pixel hits",1.0); // Count the total number of missed hits: not including the "DoubleHits"
		      H_dist_pixe->Fill(PixMissedHitsBw[HIndex].Distance);
		      H_length_pixe->Fill(PixMissedHitsBw[HIndex].PropLength);
		      H_dist_length->Fill(PixMissedHitsBw[HIndex].Distance, PixMissedHitsBw[HIndex].PropLength);
		      H_normgeom->Fill(TrackMap.find(h->geographicalId().subdetId())->second, PixMissedHitsBw[HIndex].LayerId, 1.0);
		      H_hitpos->Fill(PixMissedHitsBw[HIndex].TrackHitPosition);
		      H_normhitpos->Fill((double)PixMissedHitsBw[HIndex].TrackHitPosition/(double)(Traj.size()-1));
			       
		      if (theTransientHitAlongB->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel)
			{
			  H_2DmapXY_B->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			  H_2DmapRZ_B->Fill(theTransientHitAlongB->globalPosition().z(),
					    sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						 powf(theTransientHitAlongB->globalPosition().y(),2.)));
			}
		      else if (theTransientHitAlongB->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)
			{
			  H_2DmapXY_F->Fill(theTransientHitAlongB->globalPosition().x(), theTransientHitAlongB->globalPosition().y());
			  H_2DmapRZ_F->Fill(theTransientHitAlongB->globalPosition().z(),
					    sqrt(powf(theTransientHitAlongB->globalPosition().x(),2.) +
						 powf(theTransientHitAlongB->globalPosition().y(),2.)));
			}				      
			       
		      if (ClusterPlots == true)
			{
			  stringstream HistoName;					
			  HistoName << "Pixel module ID: " << h->geographicalId()() << " #" << ++NPixelModules;
			  if (PrintMsg == true) cout << "\n" << HistoName.str().c_str() << endl;					
			  H_general_pix_module = SubDirClus->make<TH2F>(HistoName.str().c_str(), HistoName.str().c_str(), 417, -0.5, 416.5, 161, -0.5, 160.5);
			  H_general_pix_module->SetXTitle("Columns");
			  H_general_pix_module->SetYTitle("Rows");
			  H_general_pix_module->SetZTitle("N.Clusters");
			}
		      int CountHitsOnMod = 0;
		      int CountSimHitsOnMod = 0;
		      bool FirstMuon = false;
		      bool OtherMuon = true;
		      for (SiPixelRecHitCollection::const_iterator jIt = rechitspixel->begin(); jIt != rechitspixel->end(); jIt++) { // Loop on detectors
			for (edmNew::DetSet<SiPixelRecHit>::const_iterator hIt = jIt->begin(); hIt != jIt->end(); hIt++) { // Loop on hits on detector
			  if (h->geographicalId()() == hIt->geographicalId()())
			    {
			      CountHitsOnMod++;
			
			      // ###### For simulated data ######
			      if (SimData == true)
				{	       
				  TransientTrackingRecHit::RecHitPointer theTH = theBuilder->build(&*hIt);
				  // Extraction of the SimHits associated to a cluster with all their properties (e.g.
				  // particle type and process type)
				  // There might be several IDs because of the different strips/pixels in the cluster
				  std::vector<PSimHit> PSimTrackHits = hitAssociator->associateHit(*theTH);
				  if (PrintMsg == true) cout << "New pixel cluster" << endl;
				  if ((PixMissedHitsBw[HIndex].IndexDet == jIt-rechitspixel->begin()) && (PixMissedHitsBw[HIndex].IndexHit == hIt-jIt->begin()) &&
				      (abs(PSimTrackHits[0].particleType()) == ParticleType) && (PSimTrackHits[0].processType() == ProcessType)) FirstMuon = true;
				  else if (((PixMissedHitsBw[HIndex].IndexDet != jIt-rechitspixel->begin()) || (PixMissedHitsBw[HIndex].IndexHit == hIt-jIt->begin())) &&
					   ((abs(PSimTrackHits[0].particleType()) != ParticleType) || (PSimTrackHits[0].processType() != ProcessType))) OtherMuon = false;

				  if (PrintMsg == true)
				    for (unsigned int ii = 0; ii < PSimTrackHits.size(); ii++)
				      cout << "Pixel cluster\tparticle type: " << PSimTrackHits[ii].particleType() << "\tprocess type: " << PSimTrackHits[ii].processType() << endl;

				  // Extraction of the IDs associated to a hit (= cluster)
				  // There might be several IDs because of the different strips/pixels in the cluster
				  std::vector<SimHitIdpr> SimTrackIDs = hitAssociator->associateHitId(*theTH);
				  unsigned int jj;
				  for (unsigned int ii = 0; ii < SimTrackIDs.size(); ii++) {
				    for (jj = 0; jj < TrParIDs.size(); jj++) {
				      if (SimTrackIDs[ii].first == TrParIDs[jj])
					{
					  CountSimHitsOnMod++;
						 
					  if (ClusterPlots == true)
					    for (unsigned int k = 0; k < hIt->cluster()->pixels().size(); k++)
					      // H_general_pix_module->Fill(hIt->cluster()->pixels()[k].y & 511, hIt->cluster()->pixels()[k].x, hIt->cluster()->pixels()[k].adc);
					      H_general_pix_module->Fill(hIt->cluster()->pixels()[k].y & 511, hIt->cluster()->pixels()[k].x, CountSimHitsOnMod);
					  // Only for the Y coordinate:
					  // 9 bits for Col information
					  // 1 bit for larger clusters than 9x33
					  // 6 bits for quality information
						 
					  break;
					}
				    }
				    if (jj != TrParIDs.size()) break;
				  }
				}
			    }
			}
		      }
		      H_rechits_on_mod_pixe->Fill(CountHitsOnMod);
		      H_simhits_on_mod_pixe->Fill(((double)CountSimHitsOnMod)/((double)CountHitsOnMod)*100.);
		      
		      // Take into account multiple track-compatible hits
		      if ((SimData == true) && (CountHitsOnMod > 0) && (FirstMuon == true) && (OtherMuon == true))
			H_category_multy_pix->Fill("AllMu",1.0);
		      else if ((SimData == true) && (CountHitsOnMod > 0) && (FirstMuon == true) && (OtherMuon == false))
			H_category_multy_pix->Fill("1Mu @least 1NoMu",1.0);
		      else if ((SimData == true) && (CountHitsOnMod > 0) && (FirstMuon == false) && (OtherMuon == false))
			H_category_multy_pix->Fill("NoMu",1.0);
		      else if ((SimData == true) && (CountHitsOnMod > 0) && (FirstMuon == false) && (OtherMuon == true))
			H_category_multy_pix->Fill("1NoMu o/AllMu",1.0);
			       
		      // Take into account single track-compatible hits
		      if (CountHitsOnMod > 0)
			{
			  std::vector<PSimHit> PSimTrackHits;
			  if (SimData == true) PSimTrackHits = hitAssociator->associateHit(*theTransientHit);

			  if ((SimData == false) ||
			      ((PSimTrackHits.size() != 0) && (abs(PSimTrackHits[0].particleType()) == ParticleType) && (PSimTrackHits[0].processType() == ProcessType)))
			    {
			      H_pix_charge->Fill(h->cluster()->charge());

			      if (h->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel)
				{
				  PixelBarrelName BPixName(h->geographicalId());
				  if (BPixName.moduleType() == PixelBarrelName::v2x8)
				    {
				      PixelHit2x8Counter++;
				      for (unsigned int k = 0; k < h->cluster()->pixels().size(); k++)
					H_2x8_pix_module->Fill(h->cluster()->pixels()[k].y & 511, h->cluster()->pixels()[k].x, h->cluster()->pixels()[k].adc);
				    }
				  else if (BPixName.moduleType() == PixelBarrelName::v1x8)
				    {
				      PixelHit1x8Counter++;
				      for (unsigned int k = 0; k < h->cluster()->pixels().size(); k++)
					H_1x8_pix_module->Fill(h->cluster()->pixels()[k].y & 511, h->cluster()->pixels()[k].x, h->cluster()->pixels()[k].adc);
				    }
				}
				       
			      else if (h->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)
				{
				  PixelEndcapName FPixName(h->geographicalId());
				  if (FPixName.moduleType() == PixelEndcapName::v1x2)
				    {
				      PixelHit1x2Counter++;
				      for (unsigned int k = 0; k < h->cluster()->pixels().size(); k++)
					H_1x2_pix_module->Fill(h->cluster()->pixels()[k].y & 511, h->cluster()->pixels()[k].x, h->cluster()->pixels()[k].adc);
				    }
				  else if (FPixName.moduleType() == PixelEndcapName::v1x5)
				    {
				      PixelHit1x5Counter++;
				      for (unsigned int k = 0; k < h->cluster()->pixels().size(); k++)
					H_1x5_pix_module->Fill(h->cluster()->pixels()[k].y & 511, h->cluster()->pixels()[k].x, h->cluster()->pixels()[k].adc);
				    }
				  else if (FPixName.moduleType() == PixelEndcapName::v2x3)
				    {
				      PixelHit2x3Counter++;
				      for (unsigned int k = 0; k < h->cluster()->pixels().size(); k++)
					H_2x3_pix_module->Fill(h->cluster()->pixels()[k].y & 511, h->cluster()->pixels()[k].x, h->cluster()->pixels()[k].adc);
				    }
				  else if (FPixName.moduleType() == PixelEndcapName::v2x4)
				    {
				      PixelHit2x4Counter++;
				      for (unsigned int k = 0; k < h->cluster()->pixels().size(); k++)
					H_2x4_pix_module->Fill(h->cluster()->pixels()[k].y & 511, h->cluster()->pixels()[k].x, h->cluster()->pixels()[k].adc);
				    }
				  else if (FPixName.moduleType() == PixelEndcapName::v2x5)
				    {
				      PixelHit2x5Counter++;
				      for (unsigned int k = 0; k < h->cluster()->pixels().size(); k++)
					H_2x5_pix_module->Fill(h->cluster()->pixels()[k].y & 511, h->cluster()->pixels()[k].x, h->cluster()->pixels()[k].adc);
				    }
				}
			    }
			}

		      if (PrintMsg == true)
			{
			  cout << "\nI found a PIXEL hit that doesn't belong to the track in the event: " <<  iEvent.id() << endl;
			  cout << "\tHit distance from the track: " << PixMissedHitsBw[HIndex].Distance << " (<=" << MaxDistance << " cm)" << "\tPropagation length: " << PixMissedHitsBw[HIndex].PropLength << endl;
			  cout << "\tChi2 compatibility with TransientRecHit: " << PixMissedHitsBw[HIndex].Chi2 << " (<=" << MaxChi2 << ")" << endl;
			  cout << "\tPropagation from track hit #" << PixMissedHitsBw[HIndex].TrackHitPosition << "\tTSOS pos.: " << TrajStateOnSurf.globalPosition() << "\tHit pos.: " << recHit->globalPosition() << endl;
			  cout << "\t\tSubDetId: " << recHit->geographicalId().subdetId() << "\tID: " << recHit->geographicalId()() << endl;
			  cout << "\tResult of propagation: " << TrajStateOnNewSurf.globalPosition() << endl;
			  cout << "\tTransRecHit #" << h-j->begin() << "\tID: " << theTransientHitAlongB->geographicalId()() << "\tPos.: " << theTransientHitAlongB->globalPosition() << endl;
			  cout << "\t\tSubDetId: " << h->geographicalId().subdetId() << "\tLayerID: " << PixMissedHitsBw[HIndex].LayerId;
			  cout << "\tMin row: " << h->cluster()->minPixelRow() << "\tMin col: " << h->cluster()->minPixelCol() << "\tMax row: " << h->cluster()->maxPixelRow() << "\tMax col: " << h->cluster()->maxPixelCol() << endl;
			}
		    }
		}
	      }

	    if (MomentumCounter == true)
	      {
		H_p->Fill(fabs(1/track.qoverp()));
		H_pt->Fill(track.pt());
		H_eta->Fill(track.eta());
		H_theta->Fill(track.theta());
		H_phi->Fill(track.phi());
		H_non_valid_hits->Fill(track.numberOfLostHits());
		
		// ###### For simulated data ######
		if (SimData == true)
		  {
		    int CountSimHits = 0;
		    for (std::vector<TrajectoryMeasurement>::const_iterator itTraj = Traj.begin(); itTraj != Traj.end(); itTraj++) { // Loop on hits on detector		
		      // Extraction of the hits from the track
		      recHit = itTraj->recHit();
		      bool IsSimHit = false;
		      // Extraction of the IDs associated to a hit (= cluster)
		      // There might be several IDs because of the different strips/pixels in the cluster
		      std::vector<SimHitIdpr> SimTrackIDs = hitAssociator->associateHitId(*recHit);
		      unsigned int jj;
		      for (unsigned int ii = 0; ii < SimTrackIDs.size(); ii++) {
			for (jj = 0; jj < TrParIDs.size(); jj++) {
			  if (SimTrackIDs[ii].first == TrParIDs[jj])
			    {
			      IsSimHit = true;
			      break;
			    }
			}
			if (jj != TrParIDs.size()) break;
		      }
		      if (IsSimHit == true) CountSimHits++;
		    }
		    H_non_sim_hits->Fill(track.numberOfValidHits() - CountSimHits);
		  }
	      }
	  }
	else if (PrintMsg == true) cout << "Track was rejected because: either number of valid hits<" << MinNValidHits << " or Pt<" << MinPt << "GeV/c or |eta|>" << MaxEta << endl;
	StrMissedHitsRphi.erase(StrMissedHitsRphi.begin(), StrMissedHitsRphi.end());
	StrMissedHitsSter.erase(StrMissedHitsSter.begin(), StrMissedHitsSter.end());
	PixMissedHitsFw.erase(PixMissedHitsFw.begin(), PixMissedHitsFw.end());
	PixMissedHitsBw.erase(PixMissedHitsBw.begin(), PixMissedHitsBw.end());
      }
    }
  
  // ###### For simulated data ######
  if (SimData == true)
    delete hitAssociator;
}


void MyHitAnalyzer::beginRun (edm::Run const&, const edm::EventSetup& iSetup)
{
  // Extraction of the geometry
  edm::ESHandle<TrackerGeometry> estracker;
  iSetup.get<TrackerDigiGeometryRecord>().get(estracker);
  tracker=&(*estracker);

  iSetup.get<IdealMagneticFieldRecord>().get(theMF);
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
  iSetup.get<TransientRecHitRecord>().get("WithTrackAngle",theBuilder);
  iSetup.get<TkStripCPERecord>().get("SimpleStripCPE",theCPE);

  // Initialization of the matrix containing rescaling factors: number of tracked hits per detector per layer
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 18; j++) {
      NHitCountOnLayers[i][j] = 0;
    }
  }

  // ###### For simulated data ######
  if (SimData == true)
    theAssociator = new TrackAssociatorByHits(ParSetTrackAss);

  edm::Service<TFileService> FileService;
  TFileDirectory SubDir2D = FileService->mkdir("2D maps");
  if (ClusterPlots == true) SubDirClus.reset(new TFileDirectory(FileService->mkdir("Clusters shape")));

  H = FileService->make<TH1F>("Number of missed hits","Number of missed hits", 6, 0., 6.);
  H->SetXTitle("");
  H->SetYTitle("Entries [#]");
  H->GetXaxis()->SetBinLabel(1,"R\\phi hits");
  H->GetXaxis()->SetBinLabel(2,"Stereo hits");
  H->GetXaxis()->SetBinLabel(3,"Pixel hits");
  H->GetXaxis()->SetBinLabel(4,"Total tracks");
  H->GetXaxis()->SetBinLabel(5,"Total valid hits");
  H->GetXaxis()->SetBinLabel(6,"Total lost hits");

  H_category_str = FileService->make<TH1F>("Categories of strip missed hits","Categories of strip missed hits", 8, 0., 8.);
  H_category_str->SetXTitle("");
  H_category_str->SetYTitle("Entries [#]");
  H_category_str->GetXaxis()->SetBinLabel(1,"Non SimHit");
  H_category_str->GetXaxis()->SetBinLabel(2,"DoubleHit");  
  H_category_str->GetXaxis()->SetBinLabel(3,"Overlap");
  H_category_str->GetXaxis()->SetBinLabel(4,"Last");
  H_category_str->GetXaxis()->SetBinLabel(5,"Edge");
  H_category_str->GetXaxis()->SetBinLabel(6,"Bef.Trans.");
  H_category_str->GetXaxis()->SetBinLabel(7,"Aft.Trans.");
  H_category_str->GetXaxis()->SetBinLabel(8,"Others");

  H_category_pix = FileService->make<TH1F>("Categories of pixel missed hits","Categories of pixel missed hits", 9, 0., 9.);
  H_category_pix->SetXTitle("");
  H_category_pix->SetYTitle("Entries [#]");
  H_category_pix->GetXaxis()->SetBinLabel(1,"Non SimHit");
  H_category_pix->GetXaxis()->SetBinLabel(2,"DoubleHit");
  H_category_pix->GetXaxis()->SetBinLabel(3,"Overlap BPix");
  H_category_pix->GetXaxis()->SetBinLabel(4,"Overlap FPix");
  H_category_pix->GetXaxis()->SetBinLabel(5,"Last");
  H_category_pix->GetXaxis()->SetBinLabel(6,"Edge");
  H_category_pix->GetXaxis()->SetBinLabel(7,"Bef.Trans.");
  H_category_pix->GetXaxis()->SetBinLabel(8,"Aft.Trans.");
  H_category_pix->GetXaxis()->SetBinLabel(9,"Others");

  H_category_multy_pix = FileService->make<TH1F>("Categories of multiple pixel missed hits","Categories of multiple pixel missed hits", 4, 0., 4.);
  H_category_multy_pix->SetXTitle("");
  H_category_multy_pix->SetYTitle("Entries [#]");
  H_category_multy_pix->GetXaxis()->SetBinLabel(1,"AllMu");
  H_category_multy_pix->GetXaxis()->SetBinLabel(2,"1Mu @least 1NoMu");
  H_category_multy_pix->GetXaxis()->SetBinLabel(3,"NoMu");
  H_category_multy_pix->GetXaxis()->SetBinLabel(4,"1NoMu o/AllMu");

  H_dist_rphi = FileService->make<TH1F>("Distance of rphi missed hits","Distance of rphi missed hits", (int)(50*MaxDistance/0.05), 0., MaxDistance);
  H_dist_rphi->SetXTitle("Distance [cm]");
  H_dist_rphi->SetYTitle("Entries [#]");

  H_dist_ster = FileService->make<TH1F>("Distance of stereo missed hits","Distance of stereo missed hits", (int)(50*MaxDistance/0.05), 0., MaxDistance);
  H_dist_ster->SetXTitle("Distance [cm]");
  H_dist_ster->SetYTitle("Entries [#]");

  H_dist_pixe = FileService->make<TH1F>("Distance of pixel missed hits","Distance of pixel missed hits", (int)(50*MaxDistance/0.05), 0., MaxDistance);
  H_dist_pixe->SetXTitle("Distance [cm]");
  H_dist_pixe->SetYTitle("Entries [#]");

  H_length_rphi = FileService->make<TH1F>("Prop length of rphi missed hits","Prop. length of rphi missed hits", 500, 0., MaxPropLength);
  H_length_rphi->SetXTitle("Length [cm]");
  H_length_rphi->SetYTitle("Entries [#]");

  H_length_ster = FileService->make<TH1F>("Prop length of stereo missed hits","Prop. length of stereo missed hits", 500, 0., MaxPropLength);
  H_length_ster->SetXTitle("Length [cm]");
  H_length_ster->SetYTitle("Entries [#]");

  H_length_pixe = FileService->make<TH1F>("Prop length of pixel missed hits","Prop. length of pixel missed hits", 500, 0., MaxPropLength);
  H_length_pixe->SetXTitle("Length [cm]");
  H_length_pixe->SetYTitle("Entries [#]");

  H_chi2 = FileService->make<TH1F>("Chi2","Chi2", 100, 0., MaxChi2);
  H_chi2->SetXTitle("\\chi^{2}");
  H_chi2->SetYTitle("Entries [#]");

  H_p = FileService->make<TH1F>("Momentum of tracks with missed hits","Momentum of tracks with missed hits", 100, 0., MaxMomentum);
  H_p->SetXTitle("Momentum [GeV/c]");
  H_p->SetYTitle("Entries [#]");

  H_pt = FileService->make<TH1F>("Pt of tracks with missed hits","Pt of tracks with missed hits", 100, 0., MaxMomentum);
  H_pt->SetXTitle("p_{t} [GeV/c]");
  H_pt->SetYTitle("Entries [#]");

  H_eta = FileService->make<TH1F>("Eta of tracks with missed hits","Eta of tracks with missed hits", 100, -2.5, 2.5);
  H_eta->SetXTitle("\\eta");
  H_eta->SetYTitle("Entries [#]");

  H_theta = FileService->make<TH1F>("Polar angle of tracks with missed hits","Polar angle of tracks with missed hits", 100, 0., 4.);
  H_theta->SetXTitle("\\theta [rad]");
  H_theta->SetYTitle("Entries [#]");

  H_phi = FileService->make<TH1F>("Azimuthal angle of tracks with missed hits","Azimuthal angle of tracks with missed hits", 100, -4., 4.);
  H_phi->SetXTitle("\\phi [rad]");
  H_phi->SetYTitle("Entries [#]");

  H_non_valid_hits = FileService->make<TH1F>("Number of non valid hits in the track with missed hits","Number of non valid hits in the track with missed hits", 10, -0.5, 9.5);
  H_non_valid_hits->SetXTitle("Number of non valid hits");
  H_non_valid_hits->SetYTitle("Entries [#]");

  H_non_sim_hits = FileService->make<TH1F>("Number of non simhits in the track with missed hits","Number of non sim-hits in the track with missed hits", 10, -0.5, 9.5);
  H_non_sim_hits->SetXTitle("Number of non sim-hits");
  H_non_sim_hits->SetYTitle("Entries [#]");

  H_hitpos = FileService->make<TH1F>("Hit position within the track","Hit position within the track", 30, -0.5, 29.5);
  H_hitpos->SetXTitle("Hit position [#]");
  H_hitpos->SetYTitle("Entries [#]");

  H_normhitpos = FileService->make<TH1F>("Norm hit position within the track","Norm. hit position within the track", 21, 0., 1.05);
  H_normhitpos->SetXTitle("Norm. hit position");
  H_normhitpos->SetYTitle("Entries [#]");

  H_rechits_on_mod_rphi = FileService->make<TH1F>("Number of hits on the module with an RPHI missed hit","Number of hits on the module with an RPHI missed hit", 20, -0.5, 19.5);
  H_rechits_on_mod_rphi->SetXTitle("Number of hits");
  H_rechits_on_mod_rphi->SetYTitle("Entries [#]");

  H_rechits_on_mod_ster = FileService->make<TH1F>("Number of hits on the module with a STEREO missed hit","Number of hits on the module with a STEREO missed hit", 20, -0.5, 19.5);
  H_rechits_on_mod_ster->SetXTitle("Number of hits");
  H_rechits_on_mod_ster->SetYTitle("Entries [#]");

  H_rechits_on_mod_pixe = FileService->make<TH1F>("Number of hits on the module with a PIXEL missed hit","Number of hits on the module with a PIXEL missed hit", 20, -0.5, 19.5);
  H_rechits_on_mod_pixe->SetXTitle("Number of hits");
  H_rechits_on_mod_pixe->SetYTitle("Entries [#]");

  H_simhits_on_mod_rphi = FileService->make<TH1F>("Percentage of SimHits on a RPHI module with a missed hit","Percentage of SimHits on a RPHI module with a missed hit", 101, -0.5, 100.5);
  H_simhits_on_mod_rphi->SetXTitle("Percentage [%]");
  H_simhits_on_mod_rphi->SetYTitle("Entries [#]");

  H_simhits_on_mod_ster = FileService->make<TH1F>("Percentage of SimHits on a STEREO module with a missed hit","Percentage of SimHits on a STEREO module with a missed hit", 101, -0.5, 100.5);
  H_simhits_on_mod_ster->SetXTitle("Percentage [%]");
  H_simhits_on_mod_ster->SetYTitle("Entries [#]");

  H_simhits_on_mod_pixe = FileService->make<TH1F>("Percentage of SimHits on a PIXEL module with a missed hit","Percentage of SimHits on a PIXEL module with a missed hit", 101, -0.5, 100.5);
  H_simhits_on_mod_pixe->SetXTitle("Percentage [%]");
  H_simhits_on_mod_pixe->SetYTitle("Entries [#]");

  H_pix_charge = FileService->make<TH1F>("Pixel missedhits cluster charge distribution","Pixel missedhits cluster charge distribution", 100, 0., 100000.);
  H_pix_charge->SetXTitle("Charge [e-]");
  H_pix_charge->SetYTitle("Entries [#]");

  H_pix_charge_doublehit_missed = FileService->make<TH1F>("Doublehit pixel missedhit charge distribution","Doublehit pixel missedhit charge distribution", 100, 0., 100000.);
  H_pix_charge_doublehit_missed->SetXTitle("Charge [e-]");
  H_pix_charge_doublehit_missed->SetYTitle("Entries [#]");

  H_pix_charge_doublehit_trked = FileService->make<TH1F>("Doublehit tracked pixelhit charge distribution","Doublehit tracked pixelhit charge distribution", 100, 0., 100000.);
  H_pix_charge_doublehit_trked->SetXTitle("Charge [e-]");
  H_pix_charge_doublehit_trked->SetYTitle("Entries [#]");

  H_pix_dist_doublehit = FileService->make<TH1F>("Distance between tracked hit and pixel missed doublehit","Distance between tracked hit and pixel missed doublehit", 200, 0., 2.);
  H_pix_dist_doublehit->SetXTitle("Distance [cm]");
  H_pix_dist_doublehit->SetYTitle("Entries [#]");

  H_str_dist_doublehit = FileService->make<TH1F>("Distance between tracked hit and strip missed doublehit","Distance between tracked hit and strip missed doublehit", 200, 0., 2.);
  H_str_dist_doublehit->SetXTitle("Distance [cm]");
  H_str_dist_doublehit->SetYTitle("Entries [#]");

  H_dist_str_samelayer = FileService->make<TH1F>("Radial distance between tracked hit and strip missed hit","Radial distance between tracked hit and strip missed hit", 200, 0., 2.);
  H_dist_str_samelayer->SetXTitle("Distance [cm]");
  H_dist_str_samelayer->SetYTitle("Entries [#]");

  H_dist_str_sameweel = FileService->make<TH1F>("Longitudinal distance between tracked hit and strip missed hit","Longitudinal distance between tracked hit and strip missed hit", 200, 0., 2.);
  H_dist_str_sameweel->SetXTitle("Distance [cm]");
  H_dist_str_sameweel->SetYTitle("Entries [#]");

  H_dist_pix_samelayer = FileService->make<TH1F>("Radial distance between tracked hit and pixel missed hit","Radial distance between tracked hit and pixel missed hit", 200, 0., 2.);
  H_dist_pix_samelayer->SetXTitle("Distance [cm]");
  H_dist_pix_samelayer->SetYTitle("Entries [#]");

  H_dist_pix_sameendcap = FileService->make<TH1F>("Longitudinal distance between tracked hit and pixel missed hit","Longitudinal distance between tracked hit and pixel missed hit", 200, 0., 2.);
  H_dist_pix_sameendcap->SetXTitle("Distance [cm]");
  H_dist_pix_sameendcap->SetYTitle("Entries [#]");

  H_normgeom = FileService->make<TH2F>("Norm missed hits per sub det","Norm. missed hits per sub det.", 6, 0., 6., 18, 0.5, 18.5);
  H_normgeom->SetXTitle("Sub.det");
  H_normgeom->SetYTitle("Layer");
  H_normgeom->SetZTitle("# missed hits / # tracked hits");
  H_normgeom->GetXaxis()->SetBinLabel(1,TrackMap[1]);
  H_normgeom->GetXaxis()->SetBinLabel(2,TrackMap[2]);
  H_normgeom->GetXaxis()->SetBinLabel(3,TrackMap[3]);
  H_normgeom->GetXaxis()->SetBinLabel(4,TrackMap[4]);
  H_normgeom->GetXaxis()->SetBinLabel(5,TrackMap[5]);
  H_normgeom->GetXaxis()->SetBinLabel(6,TrackMap[6]);

  H_geom = FileService->make<TH2F>("Valid and missed hits per sub det","Valid+missed hits per sub det.", 6, 0., 6., 18, 0.5, 18.5);
  H_geom->SetXTitle("Sub.det");
  H_geom->SetYTitle("Layer");
  H_geom->SetZTitle("# valid+missed hits");
  H_geom->GetXaxis()->SetBinLabel(1,TrackMap[1]);
  H_geom->GetXaxis()->SetBinLabel(2,TrackMap[2]);
  H_geom->GetXaxis()->SetBinLabel(3,TrackMap[3]);
  H_geom->GetXaxis()->SetBinLabel(4,TrackMap[4]);
  H_geom->GetXaxis()->SetBinLabel(5,TrackMap[5]);
  H_geom->GetXaxis()->SetBinLabel(6,TrackMap[6]);

  H_chi2_length = FileService->make<TH2F>("Chi2 vs prop length","Chi2 vs. prop. length", 100, 0., MaxChi2, 500, 0., MaxPropLength);
  H_chi2_length->SetXTitle("\\chi^{2}");
  H_chi2_length->SetYTitle("Prop. length [cm]");
  H_chi2_length->SetZTitle("Entries [#]");

  H_chi2_dist = FileService->make<TH2F>("Chi2 vs track distance from hit","Chi2 vs. track distance from hit", 100, 0., MaxChi2, (int)(50*MaxDistance/0.05), 0., MaxDistance);
  H_chi2_dist->SetXTitle("\\chi^{2}");
  H_chi2_dist->SetYTitle("Distance [cm]");
  H_chi2_dist->SetZTitle("Entries [#]");

  H_dist_length = FileService->make<TH2F>("Prop length vs prop distance from hit","Prop. length vs. prop. distance from hit", (int)(50*MaxDistance/0.05), 0., MaxDistance, 500, 0., MaxPropLength);
  H_dist_length->SetXTitle("Distance [cm]");
  H_dist_length->SetYTitle("Prop. length [cm]");
  H_dist_length->SetZTitle("Entries [#]");

  H_2x8_pix_module = FileService->make<TH2F>("Hit distribution on 2x8 modules", "Hit distribution on 2x8 modules", 416, -0.5, 415.5, 160, -0.5, 159.5);
  H_2x8_pix_module->SetXTitle("Columns");
  H_2x8_pix_module->SetYTitle("Rows");
  H_2x8_pix_module->SetZTitle("Charge [e-]");

  H_1x8_pix_module = FileService->make<TH2F>("Hit distribution on 1x8 modules", "Hit distribution on 1x8 modules", 416, -0.5, 415.5, 80, -0.5, 79.5);
  H_1x8_pix_module->SetXTitle("Columns");
  H_1x8_pix_module->SetYTitle("Rows");
  H_1x8_pix_module->SetZTitle("Charge [e-]");

  H_1x2_pix_module = FileService->make<TH2F>("Hit distribution on 1x2 modules", "Hit distribution on 1x2 modules", 104, -0.5, 103.5, 80, -0.5, 79.5);
  H_1x2_pix_module->SetXTitle("Columns");
  H_1x2_pix_module->SetYTitle("Rows");
  H_1x2_pix_module->SetZTitle("Charge [e-]");

  H_1x5_pix_module = FileService->make<TH2F>("Hit distribution on 1x5 modules", "Hit distribution on 1x5 modules", 260, -0.5, 259.5, 80, -0.5, 79.5);
  H_1x5_pix_module->SetXTitle("Columns");
  H_1x5_pix_module->SetYTitle("Rows");
  H_1x5_pix_module->SetZTitle("Charge [e-]");

  H_2x3_pix_module =  FileService->make<TH2F>("Hit distribution on 2x3 modules", "Hit distribution on 2x3 modules", 156, -0.5, 155.5, 160, -0.5, 159.5);
  H_2x3_pix_module->SetXTitle("Columns");
  H_2x3_pix_module->SetYTitle("Rows");
  H_2x3_pix_module->SetZTitle("Charge [e-]");

  H_2x4_pix_module =  FileService->make<TH2F>("Hit distribution on 2x4 modules", "Hit distribution on 2x4 modules", 208, -0.5, 207.5, 160, -0.5, 159.5);
  H_2x4_pix_module->SetXTitle("Columns");
  H_2x4_pix_module->SetYTitle("Rows");
  H_2x4_pix_module->SetZTitle("Charge [e-]");

  H_2x5_pix_module =  FileService->make<TH2F>("Hit distribution on 2x5 modules", "Hit distribution on 2x5 modules", 260, -0.5, 259.5, 160, -0.5, 159.5);
  H_2x5_pix_module->SetXTitle("Columns");
  H_2x5_pix_module->SetYTitle("Rows");
  H_2x5_pix_module->SetZTitle("Charge [e-]");

  // ****** 2D Maps ******
  H_2DmapXY_B = SubDir2D.make<TH2F>("2D map XY view b","2D map XY view b.", 480, -120., 120., 480, -120., 120.);
  H_2DmapXY_B->SetXTitle("X [cm]");
  H_2DmapXY_B->SetYTitle("Y [cm]");
  H_2DmapXY_B->SetZTitle("Entries [#]");

  H_2DmapRZ_B = SubDir2D.make<TH2F>("2D map RZ view b","2D map RZ view b.", 480, -120., 120., 240, 0., 120.);
  H_2DmapRZ_B->SetXTitle("Z [cm]");
  H_2DmapRZ_B->SetYTitle("R [cm]");
  H_2DmapRZ_B->SetZTitle("Entries [#]");

  H_2DmapXY_F = SubDir2D.make<TH2F>("2D map XY view fw","2D map XY view fw.", 480, -120., 120., 480, -120., 120.);
  H_2DmapXY_F->SetXTitle("X [cm]");
  H_2DmapXY_F->SetYTitle("Y [cm]");
  H_2DmapXY_F->SetZTitle("Entries [#]");

  H_2DmapRZ_F = SubDir2D.make<TH2F>("2D map RZ view fw","2D map RZ view fw.", 1120, -280., 280., 240, 0., 120.);
  H_2DmapRZ_F->SetXTitle("Z [cm]");
  H_2DmapRZ_F->SetYTitle("R [cm]");
  H_2DmapRZ_F->SetZTitle("Entries [#]");

  // Histograms category-specific
  H_2DmapXY_PixOver_B = SubDir2D.make<TH2F>("2D map XY view BPix overlaps","2D map XY view BPix overlaps", 60, -15., 15., 60, -15., 15.);
  H_2DmapXY_PixOver_B->SetXTitle("X [cm]");
  H_2DmapXY_PixOver_B->SetYTitle("Y [cm]");
  H_2DmapXY_PixOver_B->SetZTitle("Entries [#]");

  H_2DmapRZ_PixOver_B = SubDir2D.make<TH2F>("2D map RZ view BPix overlaps","2D map RZ view BPix overlaps", 120, -30., 30., 30, 0., 15.);
  H_2DmapRZ_PixOver_B->SetXTitle("Z [cm]");
  H_2DmapRZ_PixOver_B->SetYTitle("R [cm]");
  H_2DmapRZ_PixOver_B->SetZTitle("Entries [#]");

  H_2DmapXY_PixOther_B = SubDir2D.make<TH2F>("2D map XY view BPix others","2D map XY view BPix others", 60, -15., 15., 60, -15., 15.);
  H_2DmapXY_PixOther_B->SetXTitle("X [cm]");
  H_2DmapXY_PixOther_B->SetYTitle("Y [cm]");
  H_2DmapXY_PixOther_B->SetZTitle("Entries [#]");

  H_2DmapRZ_PixOther_B = SubDir2D.make<TH2F>("2D map RZ view BPix others","2D map RZ view BPix others", 120, -30., 30., 30, 0., 15.);
  H_2DmapRZ_PixOther_B->SetXTitle("Z [cm]");
  H_2DmapRZ_PixOther_B->SetYTitle("R [cm]");
  H_2DmapRZ_PixOther_B->SetZTitle("Entries [#]");

  H_2DmapXY_PixEdge_B = SubDir2D.make<TH2F>("2D map XY view BPix edges","2D map XY view BPix edges", 60, -15., 15., 60, -15., 15.);
  H_2DmapXY_PixEdge_B->SetXTitle("X [cm]");
  H_2DmapXY_PixEdge_B->SetYTitle("Y [cm]");
  H_2DmapXY_PixEdge_B->SetZTitle("Entries [#]");

  H_2DmapRZ_PixEdge_B = SubDir2D.make<TH2F>("2D map RZ view BPix edges","2D map RZ view BPix edges", 120, -30., 30., 30, 0., 15.);
  H_2DmapRZ_PixEdge_B->SetXTitle("Z [cm]");
  H_2DmapRZ_PixEdge_B->SetYTitle("R [cm]");
  H_2DmapRZ_PixEdge_B->SetZTitle("Entries [#]");

  H_2DmapXY_PixTrans_B = SubDir2D.make<TH2F>("2D map XY view BPix trans","2D map XY view BPix trans.", 60, -15., 15., 60, -15., 15.);
  H_2DmapXY_PixTrans_B->SetXTitle("X [cm]");
  H_2DmapXY_PixTrans_B->SetYTitle("Y [cm]");
  H_2DmapXY_PixTrans_B->SetZTitle("Entries [#]");

  H_2DmapRZ_PixTrans_B = SubDir2D.make<TH2F>("2D map RZ view BPix trans","2D map RZ view BPix trans.", 120, -30., 30., 30, 0., 15.);
  H_2DmapRZ_PixTrans_B->SetXTitle("Z [cm]");
  H_2DmapRZ_PixTrans_B->SetYTitle("R [cm]");
  H_2DmapRZ_PixTrans_B->SetZTitle("Entries [#]");

  H_2DmapXY_PixDoubleHit_B = SubDir2D.make<TH2F>("2D map XY view BPix doublehit","2D map XY view BPix doublehit", 60, -15., 15., 60, -15., 15.);
  H_2DmapXY_PixDoubleHit_B->SetXTitle("X [cm]");
  H_2DmapXY_PixDoubleHit_B->SetYTitle("Y [cm]");
  H_2DmapXY_PixDoubleHit_B->SetZTitle("Entries [#]");

  H_2DmapRZ_PixDoubleHit_B = SubDir2D.make<TH2F>("2D map RZ view BPix doublehit","2D map RZ view BPix doublehit", 120, -30., 30., 30, 0., 15.);
  H_2DmapRZ_PixDoubleHit_B->SetXTitle("Z [cm]");
  H_2DmapRZ_PixDoubleHit_B->SetYTitle("R [cm]");
  H_2DmapRZ_PixDoubleHit_B->SetZTitle("Entries [#]");

  H_2DmapXY_PixOver_F = SubDir2D.make<TH2F>("2D map XY view FPix overlaps","2D map XY view FPix overlaps", 80, -20., 20., 80, -20., 20.);
  H_2DmapXY_PixOver_F->SetXTitle("X [cm]");
  H_2DmapXY_PixOver_F->SetYTitle("Y [cm]");
  H_2DmapXY_PixOver_F->SetZTitle("Entries [#]");

  H_2DmapRZ_PixOver_F = SubDir2D.make<TH2F>("2D map RZ view FPix overlaps","2D map RZ view FPix overlaps", 220, -55., 55., 40, 0., 20.);
  H_2DmapRZ_PixOver_F->SetXTitle("Z [cm]");
  H_2DmapRZ_PixOver_F->SetYTitle("R [cm]");
  H_2DmapRZ_PixOver_F->SetZTitle("Entries [#]");

  H_2DmapXY_PixOther_F = SubDir2D.make<TH2F>("2D map XY view FPix others","2D map XY view FPix others", 80, -20., 20., 80, -20., 20.);
  H_2DmapXY_PixOther_F->SetXTitle("X [cm]");
  H_2DmapXY_PixOther_F->SetYTitle("Y [cm]");
  H_2DmapXY_PixOther_F->SetZTitle("Entries [#]");

  H_2DmapRZ_PixOther_F = SubDir2D.make<TH2F>("2D map RZ view FPix others","2D map RZ view FPix others", 220, -55., 55., 40, 0., 20.);
  H_2DmapRZ_PixOther_F->SetXTitle("Z [cm]");
  H_2DmapRZ_PixOther_F->SetYTitle("R [cm]");
  H_2DmapRZ_PixOther_F->SetZTitle("Entries [#]");

  H_2DmapXY_PixEdge_F = SubDir2D.make<TH2F>("2D map XY view FPix edges","2D map XY view FPix edges", 80, -20., 20., 80, -20., 20.);
  H_2DmapXY_PixEdge_F->SetXTitle("X [cm]");
  H_2DmapXY_PixEdge_F->SetYTitle("Y [cm]");
  H_2DmapXY_PixEdge_F->SetZTitle("Entries [#]");

  H_2DmapRZ_PixEdge_F = SubDir2D.make<TH2F>("2D map RZ view FPix edges","2D map RZ view FPix edges", 220, -55., 55., 40, 0., 20.);
  H_2DmapRZ_PixEdge_F->SetXTitle("Z [cm]");
  H_2DmapRZ_PixEdge_F->SetYTitle("R [cm]");
  H_2DmapRZ_PixEdge_F->SetZTitle("Entries [#]");

  H_2DmapXY_PixTrans_F = SubDir2D.make<TH2F>("2D map XY view FPix trans","2D map XY view FPix trans.", 80, -20., 20., 80, -20., 20.);
  H_2DmapXY_PixTrans_F->SetXTitle("X [cm]");
  H_2DmapXY_PixTrans_F->SetYTitle("Y [cm]");
  H_2DmapXY_PixTrans_F->SetZTitle("Entries [#]");

  H_2DmapRZ_PixTrans_F = SubDir2D.make<TH2F>("2D map RZ view FPix trans","2D map RZ view FPix trans.", 220, -55., 55., 40, 0., 20.);
  H_2DmapRZ_PixTrans_F->SetXTitle("Z [cm]");
  H_2DmapRZ_PixTrans_F->SetYTitle("R [cm]");
  H_2DmapRZ_PixTrans_F->SetZTitle("Entries [#]");

  H_2DmapXY_PixDoubleHit_F = SubDir2D.make<TH2F>("2D map XY view FPix doublehit","2D map XY view FPix doublehit", 80, -20., 20., 80, -20., 20.);
  H_2DmapXY_PixDoubleHit_F->SetXTitle("X [cm]");
  H_2DmapXY_PixDoubleHit_F->SetYTitle("Y [cm]");
  H_2DmapXY_PixDoubleHit_F->SetZTitle("Entries [#]");

  H_2DmapRZ_PixDoubleHit_F = SubDir2D.make<TH2F>("2D map RZ view FPix doublehit","2D map RZ view FPix doublehit", 220, -55., 55., 40, 0., 20.);
  H_2DmapRZ_PixDoubleHit_F->SetXTitle("Z [cm]");
  H_2DmapRZ_PixDoubleHit_F->SetYTitle("R [cm]");
  H_2DmapRZ_PixDoubleHit_F->SetZTitle("Entries [#]");

  H_2DmapXY_StrOver_B = SubDir2D.make<TH2F>("2D map XY view BStr overlaps","2D map XY view BStr overlaps", 480, -120., 120., 480, -120., 120.);
  H_2DmapXY_StrOver_B->SetXTitle("X [cm]");
  H_2DmapXY_StrOver_B->SetYTitle("Y [cm]");
  H_2DmapXY_StrOver_B->SetZTitle("Entries [#]");

  H_2DmapRZ_StrOver_B = SubDir2D.make<TH2F>("2D map RZ view BStr overlaps","2D map RZ view BStr overlaps", 480, -120., 120., 240, 0., 120.);
  H_2DmapRZ_StrOver_B->SetXTitle("Z [cm]");
  H_2DmapRZ_StrOver_B->SetYTitle("R [cm]");
  H_2DmapRZ_StrOver_B->SetZTitle("Entries [#]");

  H_2DmapXY_StrLast_B = SubDir2D.make<TH2F>("2D map XY view BStr last","2D map XY view BStr last", 480, -120., 120., 480, -120., 120.);
  H_2DmapXY_StrLast_B->SetXTitle("X [cm]");
  H_2DmapXY_StrLast_B->SetYTitle("Y [cm]");
  H_2DmapXY_StrLast_B->SetZTitle("Entries [#]");

  H_2DmapRZ_StrLast_B = SubDir2D.make<TH2F>("2D map RZ view BStr last","2D map RZ view BStr last", 480, -120., 120., 240, 0., 120.);
  H_2DmapRZ_StrLast_B->SetXTitle("Z [cm]");
  H_2DmapRZ_StrLast_B->SetYTitle("R [cm]");
  H_2DmapRZ_StrLast_B->SetZTitle("Entries [#]");

  H_2DmapXY_StrTrans_B = SubDir2D.make<TH2F>("2D map XY view BStr trans","2D map XY view BStr trans.", 480, -120., 120., 480, -120., 120.);
  H_2DmapXY_StrTrans_B->SetXTitle("X [cm]");
  H_2DmapXY_StrTrans_B->SetYTitle("Y [cm]");
  H_2DmapXY_StrTrans_B->SetZTitle("Entries [#]");

  H_2DmapRZ_StrTrans_B = SubDir2D.make<TH2F>("2D map RZ view BStr trans","2D map RZ view BStr trans.", 480, -120., 120., 240, 0., 120.);
  H_2DmapRZ_StrTrans_B->SetXTitle("Z [cm]");
  H_2DmapRZ_StrTrans_B->SetYTitle("R [cm]");
  H_2DmapRZ_StrTrans_B->SetZTitle("Entries [#]");

  H_2DmapXY_StrOther_B = SubDir2D.make<TH2F>("2D map XY view BStr other","2D map XY view BStr other", 480, -120., 120., 480, -120., 120.);
  H_2DmapXY_StrOther_B->SetXTitle("X [cm]");
  H_2DmapXY_StrOther_B->SetYTitle("Y [cm]");
  H_2DmapXY_StrOther_B->SetZTitle("Entries [#]");

  H_2DmapRZ_StrOther_B = SubDir2D.make<TH2F>("2D map RZ view BStr other","2D map RZ view BStr other", 480, -120., 120., 240, 0., 120.);
  H_2DmapRZ_StrOther_B->SetXTitle("Z [cm]");
  H_2DmapRZ_StrOther_B->SetYTitle("R [cm]");
  H_2DmapRZ_StrOther_B->SetZTitle("Entries [#]");

  H_2DmapXY_StrOver_F = SubDir2D.make<TH2F>("2D map XY view FStr overlaps","2D map XY view FStr overlaps", 480, -120., 120., 480, -120., 120.);
  H_2DmapXY_StrOver_F->SetXTitle("X [cm]");
  H_2DmapXY_StrOver_F->SetYTitle("Y [cm]");
  H_2DmapXY_StrOver_F->SetZTitle("Entries [#]");

  H_2DmapRZ_StrOver_F = SubDir2D.make<TH2F>("2D map RZ view FStr overlaps","2D map RZ view FStr overlaps", 1120, -280., 280., 240, 0., 120.);
  H_2DmapRZ_StrOver_F->SetXTitle("Z [cm]");
  H_2DmapRZ_StrOver_F->SetYTitle("R [cm]");
  H_2DmapRZ_StrOver_F->SetZTitle("Entries [#]");

  H_2DmapXY_StrLast_F = SubDir2D.make<TH2F>("2D map XY view FStr last","2D map XY view FStr last", 480, -120., 120., 480, -120., 120.);
  H_2DmapXY_StrLast_F->SetXTitle("X [cm]");
  H_2DmapXY_StrLast_F->SetYTitle("Y [cm]");
  H_2DmapXY_StrLast_F->SetZTitle("Entries [#]");

  H_2DmapRZ_StrLast_F = SubDir2D.make<TH2F>("2D map RZ view FStr last","2D map RZ view FStr last", 1120, -280., 280., 240, 0., 120.);
  H_2DmapRZ_StrLast_F->SetXTitle("Z [cm]");
  H_2DmapRZ_StrLast_F->SetYTitle("R [cm]");
  H_2DmapRZ_StrLast_F->SetZTitle("Entries [#]");

  H_2DmapXY_StrTrans_F = SubDir2D.make<TH2F>("2D map XY view FStr trans","2D map XY view FStr trans.", 480, -120., 120., 480, -120., 120.);
  H_2DmapXY_StrTrans_F->SetXTitle("X [cm]");
  H_2DmapXY_StrTrans_F->SetYTitle("Y [cm]");
  H_2DmapXY_StrTrans_F->SetZTitle("Entries [#]");

  H_2DmapRZ_StrTrans_F = SubDir2D.make<TH2F>("2D map RZ view FStr trans","2D map RZ view FStr trans.", 1120, -280., 280., 240, 0., 120.);
  H_2DmapRZ_StrTrans_F->SetXTitle("Z [cm]");
  H_2DmapRZ_StrTrans_F->SetYTitle("R [cm]");
  H_2DmapRZ_StrTrans_F->SetZTitle("Entries [#]");

  H_2DmapXY_StrOther_F = SubDir2D.make<TH2F>("2D map XY view FStr other","2D map XY view FStr other", 480, -120., 120., 480, -120., 120.);
  H_2DmapXY_StrOther_F->SetXTitle("X [cm]");
  H_2DmapXY_StrOther_F->SetYTitle("Y [cm]");
  H_2DmapXY_StrOther_F->SetZTitle("Entries [#]");

  H_2DmapRZ_StrOther_F = SubDir2D.make<TH2F>("2D map RZ view FStr other","2D map RZ view FStr other", 1120, -280., 280., 240, 0., 120.);
  H_2DmapRZ_StrOther_F->SetXTitle("Z [cm]");
  H_2DmapRZ_StrOther_F->SetYTitle("R [cm]");
  H_2DmapRZ_StrOther_F->SetZTitle("Entries [#]");
}


void MyHitAnalyzer::endJob ()
{ 
  if (PrintMsg == false)
    {
      cout << "\n### Counters ###" << endl;
      cout << "\t# ev. that passed the conditions: " << NEventsPass << endl; 
      cout << "\t# PixelHit2x8Counter: " << PixelHit2x8Counter << endl;
      cout << "\t# PixelHit1x8Counter: " << PixelHit1x8Counter << endl;
      cout << "\t# PixelHit1x2Counter: " << PixelHit1x2Counter << endl;
      cout << "\t# PixelHit1x5Counter: " << PixelHit1x5Counter << endl;
      cout << "\t# PixelHit2x3Counter: " << PixelHit2x3Counter << endl;
      cout << "\t# ixelHit2x4Counter: " << PixelHit2x4Counter << endl;
      cout << "\t# PixelHit2x5Counter: " << PixelHit2x5Counter << endl;
      cout << "Total:" << PixelHit2x8Counter+PixelHit1x8Counter+PixelHit1x2Counter+PixelHit1x5Counter+PixelHit2x3Counter+PixelHit2x4Counter+PixelHit2x5Counter << "\n" << endl;

      double HitCounter = 0;
      for (int i = 1; i <= 18; i++) HitCounter = HitCounter + H_normgeom->GetBinContent((int)SiStripDetId::TIB, i);
      cout << "#TIB Hits: " << HitCounter << endl;

      HitCounter = 0;
      for (int i = 1; i <= 18; i++) HitCounter = HitCounter + H_normgeom->GetBinContent((int)SiStripDetId::TOB, i);
      cout << "#TOB Hits: " << HitCounter << endl;

      HitCounter = 0;
      for (int i = 1; i <= 18; i++) HitCounter = HitCounter + H_normgeom->GetBinContent((int)SiStripDetId::TID, i);
      cout << "#TID Hits: " << HitCounter << endl;

      HitCounter = 0;
      for (int i = 1; i <= 18; i++) HitCounter = HitCounter + H_normgeom->GetBinContent((int)SiStripDetId::TEC, i);
      cout << "#TEC Hits: " << HitCounter << endl;

      HitCounter = 0;
      for (int i = 1; i <= 18; i++) HitCounter = HitCounter + H_normgeom->GetBinContent((int)PixelSubdetector::PixelBarrel, i);
      cout << "#BPix Hits: " << HitCounter << endl;

      HitCounter = 0;
      for (int i = 1; i <= 18; i++) HitCounter = HitCounter + H_normgeom->GetBinContent((int)PixelSubdetector::PixelEndcap, i);
      cout << "#FPix Hits: " << HitCounter << endl;
    }


  // ### Rescale detector entries by the number of tracked hits per layer, offline one should do H_normgeom->Divide(H_geom) ###
  for (int i = 0; i < 3; i++)
    H_geom->SetBinContent((int)PixelSubdetector::PixelBarrel,i+1,(double)NHitCountOnLayers[(int)PixelSubdetector::PixelBarrel-1][i]+H_normgeom->GetBinContent((int)PixelSubdetector::PixelBarrel,i+1));

  for (int i = 0; i < 4; i++)
    H_geom->SetBinContent((int)PixelSubdetector::PixelEndcap,i+1,(double)NHitCountOnLayers[(int)PixelSubdetector::PixelEndcap-1][i]+H_normgeom->GetBinContent((int)PixelSubdetector::PixelEndcap,i+1));

  for (int i = 0; i < 4; i++)
    H_geom->SetBinContent((int)SiStripDetId::TIB,i+1,(double)NHitCountOnLayers[(int)SiStripDetId::TIB-1][i]+H_normgeom->GetBinContent((int)SiStripDetId::TIB,i+1));

  for (int i = 0; i < 6; i++)
    H_geom->SetBinContent((int)SiStripDetId::TID,i+1,(double)NHitCountOnLayers[(int)SiStripDetId::TID-1][i]+H_normgeom->GetBinContent((int)SiStripDetId::TID,i+1));

  for (int i = 0; i < 6; i++)
    H_geom->SetBinContent((int)SiStripDetId::TOB,i+1,(double)NHitCountOnLayers[(int)SiStripDetId::TOB-1][i]+H_normgeom->GetBinContent((int)SiStripDetId::TOB,i+1));

  for (int i = 0; i < 18; i++)
    H_geom->SetBinContent((int)SiStripDetId::TEC,i+1,(double)NHitCountOnLayers[(int)SiStripDetId::TEC-1][i]+H_normgeom->GetBinContent((int)SiStripDetId::TEC,i+1));
}


//define this as a plug-in
DEFINE_FWK_MODULE(MyHitAnalyzer);

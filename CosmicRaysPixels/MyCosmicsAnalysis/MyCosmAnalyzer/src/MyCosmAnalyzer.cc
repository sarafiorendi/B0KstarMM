// -*- C++ -*-
//
// Package:    MyCosmAnalyzer
// Class:      MyCosmAnalyzer
// 
/**\class MyCosmAnalyzer MyCosmAnalyzer.cc MyCosmicsAnalysis/MyCosmAnalyzer/src/MyCosmAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Mauro Dinardo
//         Created:  Thu May 28 12:50:56 CEST 2009
// $Id: MyCosmAnalyzer.cc,v 1.1 2010/10/06 18:38:37 dinardo Exp $

// ###### TO DO ######
// Creare ntuple per poter fre tagli
// ####################

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
#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMaskAlgoTrigRcd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMask.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>


using namespace std;
using namespace reco;
using namespace trigger;


class MyCosmAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MyCosmAnalyzer (const edm::ParameterSet&);
      ~MyCosmAnalyzer ();


   private:
      virtual void beginJob ();
      virtual void analyze (const edm::Event&, const edm::EventSetup&);
      virtual void endJob ();
      virtual void MuonTiming (const edm::Event& iEvent, const bool PrintMsg);
      virtual void L1Analyzer (const edm::Event& iEvent,
			       const edm::EventSetup& iSetup,
			       const bool PrintMsg);
      virtual void HLTAnalyzer (const edm::Event& iEvent,
				const bool PrintMsg);


      // Histos 
      TH1F* histo_x;
      TH1F* histo_y;
      TH2F* histo_xy;
      TH2F* histo_alphax;
      TH2F* histo_betay;
      TH1F* histo_charge;
      TH1F* histo_t;
      TH1F* histo_terr;

      const TrackerGeometry* tracker;
      // Parameters from cfg file
      edm::InputTag trackTags_;   // Used to select what tracks to read from configuration file
      edm::InputTag inputTag_;    // Used for trigger
      std::string trackRefitter_; // Used to select the track refitter:ctfWithMaterialTracksP5
                                  // rsWithMaterialTracksP5
                                  // cosmictrackfinderP5
};


MyCosmAnalyzer::MyCosmAnalyzer (const edm::ParameterSet& iConfig)
{
  trackRefitter_ = iConfig.getParameter<std::string>("trajectoryInput");
  trackTags_     = iConfig.getUntrackedParameter<edm::InputTag>("tracks");
  inputTag_      = iConfig.getParameter<edm::InputTag>("inputTag");
  // Now do what ever initialization is needed
}


MyCosmAnalyzer::~MyCosmAnalyzer ()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


void MyCosmAnalyzer::analyze (const edm::Event& iEvent,
			      const edm::EventSetup& iSetup)
{
  // Extraction of the geometry
  edm::ESHandle<TrackerGeometry> estracker;
  iSetup.get<TrackerDigiGeometryRecord>().get(estracker);
  tracker=&(*estracker);

  const bool PrintMsg = false;

  if (PrintMsg == true)
  {
    cout << "\n\n###### AN EVENT STARTED ######" << endl;
    cout << "\n### PixelAnalyzer ###" << endl;
  }

  MuonTiming  (iEvent, PrintMsg);
  L1Analyzer  (iEvent, iSetup, PrintMsg);
  HLTAnalyzer (iEvent, PrintMsg);

  // ###### TRACKS ######
  edm::Handle<TrajTrackAssociationCollection> tracks;
  iEvent.getByLabel(trackRefitter_, tracks); // Refitter Type

  // Loop over the tracks of the event
  for(TrajTrackAssociationCollection::const_iterator itTrack = tracks->begin();
      itTrack != tracks->end();
      ++itTrack) {
    
    const Track& track      = *itTrack->val;
    const Trajectory& Traj_ = *itTrack->key;
    std::vector<TrajectoryMeasurement> Traj = Traj_.measurements();
    
    if (PrintMsg == true)
      {
	cout << "\n*** A track started ***" << endl;
	cout << "Number of hits: " << track.found();
	cout << "\tCharge: " << track.charge();
	cout << "\tdz: " << track.dz();
	cout << "\ttheta: " << track.theta();
	cout << "\td0: " << track.d0();
	cout << "\tphi: " << track.phi();
	cout << "\tqoverp: " << track.qoverp() << endl;
	cout << "*** Loop over the hits ***" << endl;
      }

    int i = 0;
    // Loop over the hit of the track
    for(std::vector<TrajectoryMeasurement>::const_iterator itTraj = Traj.begin();
	itTraj != Traj.end();
	itTraj++) {
      i++;
      // Hit extraction
      TransientTrackingRecHit::ConstRecHitPointer recHit = itTraj->recHit();
      // Is it a valid hit AND Is it on the Tracker AND Is it on BPix or FPix?
      if ((recHit->isValid()) &&
	  (recHit->geographicalId().det() == DetId::Tracker) &&
	  ((recHit->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) ||
	   (recHit->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)))
	{
	  const SiPixelRecHit* recHitPix = dynamic_cast<const SiPixelRecHit*>(recHit->hit());
	  if (PrintMsg == true)
	    {
	      cout << "Valid pixel hit found:";
	      cout << "\tHit number: " << i;
	      cout << "\tDetId: " << recHit->geographicalId().det();
	      cout << "\tSubDetId: " << recHit->geographicalId().subdetId() << endl;
	      cout << "Global X: " << recHit->globalPosition().x();
	      cout << "\tGlobal Y: " << recHit->globalPosition().y();
	      cout << "\tGlobal Z: " << recHit->globalPosition().z();
	      cout << "\tLocal X: " << recHit->localPosition().x();
	      cout << "\tLocal Y: " << recHit->localPosition().y();
	    }

	  SiPixelRecHit::ClusterRef const& cluster = recHitPix->cluster();

	  // Cluster coordinates conversion: from row-col to local x-y
	  DetId detIdObject(recHit->geographicalId());
	  const PixelGeomDetUnit* theGeomDet = dynamic_cast<const PixelGeomDetUnit*>(tracker->idToDet(detIdObject));
	  const RectangularPixelTopology* RecPixTop = dynamic_cast<const RectangularPixelTopology*>(&(theGeomDet->specificTopology()));
	  LocalPoint LP = RecPixTop->localPosition(MeasurementPoint(cluster->x(), cluster->y()));

	  // Extract local direction, alpha (angle in x-z plane), beta (angle in y-z plane)
	  TrajectoryStateCombiner TrajStateComb;
	  // Let's have N hits: 1 ... i-1 / i / i+1 ... N
	  // itTraj->forwardPredictedState() extrapolates the track calculated with all and only the hit: 1 ... i-1 to i
	  // itTraj ->backwardPredictedState() extrapolates the track calculated with all and only the hit: N ... i+1 to i
	  // Calculate for both the two extrapolations the momentum
	  // Perform a weighted sum of the two momentums
	  TrajectoryStateOnSurface TrajStateOnSurf = TrajStateComb(itTraj->forwardPredictedState(), itTraj->backwardPredictedState());
	  LocalTrajectoryParameters LocalPar = TrajStateOnSurf.localParameters();
	  // Local momentum vector in local coordinates / Modulus of the vector
	  LocalVector LocalVec = LocalPar.momentum() / LocalPar.momentum().mag();

	  if (PrintMsg == true)
	    {
	      cout << "\tCluster-LocalC X: " << LP.x();
	      cout << "\tCluster-LocalC Y: " << LP.y();
	      cout << "\tCluster X: " << cluster->x();
	      cout << "\tCluster Y: " << cluster->y();
	      cout << "\tAlpha: " << atan2(LocalVec.z(), LocalVec.x());
	      cout << "\tBeta: " << atan2(LocalVec.z(), LocalVec.y()) << endl;
	    }

	  // Fillin histograms
	  histo_charge->Fill(cluster->charge() / 1000);
	  // 	   = cluster->sizeX();
	  // 	   = cluster->sizeY();
	  // 	   = cluster->maxPixelCol();
	  // 	   = cluster->maxPixelRow();
	  // 	   = cluster->minPixelCol();
	  // 	   = cluster->minPixelRow();
	  histo_x->Fill((LP.x()-recHit->localPosition().x()) / 0.010); // Divided by the pitch on X: 100 microns
	  histo_y->Fill((LP.y()-recHit->localPosition().y()) / 0.015); // Divided by the pitch on Y: 150 microns
	  histo_xy->Fill((LP.x()-recHit->localPosition().x()) / 0.010, (LP.y()-recHit->localPosition().y()) / 0.015);
	  histo_alphax->Fill((LP.x()-recHit->localPosition().x()) / 0.010, atan2(LocalVec.z(),LocalVec.x()));
	  histo_betay->Fill((LP.y()-recHit->localPosition().y()) / 0.015, atan2(LocalVec.z(),LocalVec.y()));
	}
    }
  }
}


void MyCosmAnalyzer::HLTAnalyzer (const edm::Event& iEvent,
			      const bool PrintMsg)
{
  // ###### HLT-TRIGGER ######
  if (PrintMsg == true)
    {
      cout << "\n### HLTAnalyzer ###" << endl;
    }

  HLTConfigProvider HLTConfig;
  edm::Handle<TriggerEvent> TrigEvHandle;
  iEvent.getByLabel(edm::InputTag("hltTriggerSummaryAOD","","HLT"), TrigEvHandle);
  edm::Handle<edm::TriggerResults> TriggerResultsHandle;
  iEvent.getByLabel(edm::InputTag("TriggerResults","","HLT"), TriggerResultsHandle);
  
  if (HLTConfig.init("HLT"))
    // Check if trigger name in (new) config
    {
      const unsigned int ConfSize = HLTConfig.size();
      if (PrintMsg == true)
	{
	  HLTConfig.dump("Triggers");
	}

      for (unsigned int i=0; i != ConfSize; ++i)
	{
	  // "@" means analyze all triggers in config
	  const unsigned int TriggerIndex = HLTConfig.triggerIndex(HLTConfig.triggerName(i));
	  const unsigned int m            = HLTConfig.size(TriggerIndex);
	  const vector<string>& moduleLabels(HLTConfig.moduleLabels(TriggerIndex));

	  // Modules "Accepted" on current trigger path
	  if (TriggerResultsHandle->accept(TriggerIndex) == true)
	    {
	      if (PrintMsg == true)
		{
		  cout << "\nTrigger path: " << HLTConfig.triggerName(i) << " [" << TriggerIndex << "]";
		  cout << "\tStatus:";
		  cout << "\tWasRun = " << TriggerResultsHandle->wasrun(TriggerIndex);
		  cout << "\tAccept = " << TriggerResultsHandle->accept(TriggerIndex);
		  cout << "\tError = " << TriggerResultsHandle->error(TriggerIndex) << endl;
		}

	      const unsigned int moduleIndex = TriggerResultsHandle->index(TriggerIndex);

	      if (PrintMsg == true)
		{
		  cout << "Last active module (label/type):\t";
		  cout << moduleLabels[moduleIndex];
		  cout << "/" << HLTConfig.moduleType(moduleLabels[moduleIndex]);
		  cout << " [" << moduleIndex << " out of 0-" << (m-1) << " on this path]" << endl;
		}

	      for (unsigned int j=0; j <= moduleIndex; ++j)
		{
		  const string& moduleLabel(moduleLabels[j]);
		  const string  moduleType(HLTConfig.moduleType(moduleLabel));
		  // Check whether the module is packed up in TrigEvent product
		  const unsigned int filterIndex = TrigEvHandle->filterIndex(edm::InputTag(moduleLabel,"","HLT"));
		  if (filterIndex < TrigEvHandle->sizeFilters())
		    {
		      const Vids& VIDS (TrigEvHandle->filterIds(filterIndex));
		      const Keys& KEYS (TrigEvHandle->filterKeys(filterIndex));
		      const size_type nI(VIDS.size());
		      const size_type nK(KEYS.size());
		      const size_type n(max(nI,nK));
		      const TriggerObjectCollection& TOC(TrigEvHandle->getObjects());

		      if (PrintMsg == true)
			{
			  cout << "L3 filter in slot " << j << " (label/type):\t" << moduleLabel << "/" << moduleType << endl;
			  cout << "Trigger L3 found: " << n  << " object/s" << endl;
			  cout << "ID \t Pt \t eta \t phi \t\t mass" << endl;

			  for (size_type i=0; i != n; ++i)
			    {
			      const TriggerObject& TO(TOC[KEYS[i]]);
			      cout << TO.id() << " \t " << TO.pt() << " \t " << TO.eta() << " \t ";
			      cout << TO.phi() << " \t " << TO.mass() << endl;
			    }
			}
		    }
		}
	    }
	}
    }
}


void MyCosmAnalyzer::L1Analyzer (const edm::Event& iEvent,
				 const edm::EventSetup& iSetup,
				 const bool PrintMsg)
{
  // ###### L1-TRIGGER ######
  // Extraction of the tigger menu
  // It gets the current menu, from the edm::EventSetup, and extract the L1 algorithm names
  // It does not tell you which were masked/disabled
  edm::ESHandle<L1GtTriggerMenu> menuRcd;
  iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd);
  const L1GtTriggerMenu* menu = menuRcd.product();

  // Check which type of event triggered the acquisition
  edm::Handle<L1GlobalTriggerReadoutRecord> gtRecord;
  iEvent.getByLabel(inputTag_, gtRecord);
  DecisionWord AlgodWord = gtRecord->decisionWord();

  // Get / update the trigger mask from the EventSetup
  const L1GtTriggerMask* m_l1GtTmAlgo;
  unsigned long long m_l1GtTmAlgoCacheID = 0ULL;
  std::vector<unsigned int> m_triggerMaskAlgoTrig;
  int triggerMaskAlgoTrigBit;

  unsigned long long l1GtTmAlgoCacheID = iSetup.get<L1GtTriggerMaskAlgoTrigRcd>().cacheIdentifier();
  if (m_l1GtTmAlgoCacheID != l1GtTmAlgoCacheID)
    {
      edm::ESHandle<L1GtTriggerMask> l1GtTmAlgo;
      iSetup.get<L1GtTriggerMaskAlgoTrigRcd>().get(l1GtTmAlgo);
      m_l1GtTmAlgo = l1GtTmAlgo.product();
      m_triggerMaskAlgoTrig = m_l1GtTmAlgo->gtTriggerMask();
      m_l1GtTmAlgoCacheID = l1GtTmAlgoCacheID;
    }
 
  // Apply the algorithm trigger mask
  int iDaq = 0;
  int iBit = -1; // Bit counter
  for (std::vector<bool>::iterator itBit = AlgodWord.begin();
       itBit != AlgodWord.end();
       ++itBit) {
    iBit++;
    triggerMaskAlgoTrigBit = m_triggerMaskAlgoTrig[iBit] & (1 << iDaq);    
    if (triggerMaskAlgoTrigBit)	*itBit = false;
  }

  // Loop over the trigger algorithms
  if (PrintMsg == true)
    {
      cout << "\n### L1Analyzer ###" << endl;
    }

  for (CItAlgo algo = menu->gtAlgorithmMap().begin();
       algo != menu->gtAlgorithmMap().end();
       ++algo) {
    if ((menu->gtAlgorithmResult(algo->second.algoName(), AlgodWord) == true) && (PrintMsg == true))
      {
	cout << "Objet that passed the L1 trigger decision: " << algo->second.algoName() << endl;
      }
  }
  
  if (PrintMsg == true)
    {
      gtRecord->printGtDecision(std::cout);
      gtRecord->printTechnicalTrigger(std::cout);
    }
}


void MyCosmAnalyzer::MuonTiming (const edm::Event& iEvent,
				 const bool PrintMsg)
{
  // ###### DT MUONS: calculate the the arrival time of the muon ######
  // ###### with respect to a fixed reference, e.g. LHC-clk edge ######
  edm::Handle<MuonCollection> MuonHandle;
  iEvent.getByLabel("muons", MuonHandle);
  
  if (PrintMsg == true)
    {
      cout << "\n### MuonTiming ###" << endl;
    }

  for(MuonCollection::const_iterator it = MuonHandle->begin();
      it != MuonHandle->end();
      ++it) {
     if(it->isTimeValid() == true)
       {
	 if (PrintMsg == true)
	   {
	     cout << "Muon time: " << it->time().timeAtIpInOut << "\ttime error: " << it->time().timeAtIpInOutErr << endl;
	   }
	 histo_t->Fill(it->time().timeAtIpInOut);
	 histo_terr->Fill(it->time().timeAtIpInOutErr);
       } 
  }
}


void MyCosmAnalyzer::beginJob ()
{  
  edm::Service<TFileService> fs;
  double pi = 4.0*atan(1.0);


  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();

  histo_x      = fs->make<TH1F>("Residue-X","Residue-X",100,-1.,1.);
  histo_y      = fs->make<TH1F>("Residue-Y","Residue-Y",100,-1.,1.);
  histo_xy     = fs->make<TH2F>("Residue-XY","Residue-XY",100,-1.,1.,100,-1.,1.);
  histo_alphax = fs->make<TH2F>("Residue-X vs. Alpha","Residue-X vs. Alpha",100,-1.,1.,100,-pi,pi);
  histo_betay  = fs->make<TH2F>("Residue-Y vs. Beta","Residue-Y vs. Alpha",100,-1.,1.,100,-pi,pi);
  histo_charge = fs->make<TH1F>("Cluster Charge","Cluster Charge",100,0.,100.);
  histo_t      = fs->make<TH1F>("Muon Time","Muon Time",100,-50.,50.);
  histo_terr   = fs->make<TH1F>("Muon Time Err","Muon Time Err",100,0.,50.);

  histo_x->UseCurrentStyle();
  histo_y->UseCurrentStyle();
  histo_xy->UseCurrentStyle();
  histo_alphax->UseCurrentStyle();
  histo_betay->UseCurrentStyle();
  histo_charge->UseCurrentStyle();
  histo_t->UseCurrentStyle();
  histo_terr->UseCurrentStyle();

  histo_x->SetXTitle("X-residual [Pitch (100 \\mum)]");
  histo_x->SetYTitle("Entries [#]");
  histo_y->SetXTitle("Y-residual [Pitch (150 \\mum)]");
  histo_y->SetYTitle("Entries [#]");
  histo_xy->SetXTitle("X-residual [Pitch (100 \\mum)]");
  histo_xy->SetYTitle("Y-residual [Pitch (150 \\mum)]");
  histo_xy->SetZTitle("Entries [#]");
  histo_alphax->SetXTitle("X-residual [Pitch (100 \\mum)]");
  histo_alphax->SetYTitle("Alpha [Rad]");
  histo_alphax->SetZTitle("Entries [#]");
  histo_betay->SetXTitle("Y-residual [Pitch (150 \\mum)]");
  histo_betay->SetYTitle("Beta [Rad]");
  histo_betay->SetZTitle("Entries [#]");
  histo_charge->SetXTitle("Custer Charge [keV]");
  histo_charge->SetYTitle("Entries [#]");
  histo_t->SetXTitle("Muon Time [ns]");
  histo_t->SetYTitle("Entries [#]");
  histo_terr->SetXTitle("Muon Time Err [ns]");
  histo_terr->SetYTitle("Entries [#]");
}


void MyCosmAnalyzer::endJob ()
{
}


//define this as a plug-in
DEFINE_FWK_MODULE(MyCosmAnalyzer);

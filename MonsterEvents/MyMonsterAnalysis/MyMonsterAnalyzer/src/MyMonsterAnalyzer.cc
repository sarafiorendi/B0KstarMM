// -*- C++ -*-
//
// Package:    MyMonsterAnalyzer
// Class:      MyMonsterAnalyzer
// 
/**\class MyMonsterAnalyzer MyMonsterAnalyzer.cc MyMonsterAnalysis/MyMonsterAnalyzer/src/MyMonsterAnalyzer.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Mauro Dinardo
//         Created:  Mon Feb  1 12:50:00 CET 2010
// $Id: MyMonsterAnalyzer.cc,v 1.1 2010/10/06 18:34:53 dinardo Exp $
//
//


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
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMask.h"
#include "CondFormats/L1TObjects/interface/L1GtPrescaleFactors.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMaskAlgoTrigRcd.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMaskTechTrigRcd.h"
#include "CondFormats/DataRecord/interface/L1GtPrescaleFactorsAlgoTrigRcd.h"
#include "CondFormats/DataRecord/interface/L1GtPrescaleFactorsTechTrigRcd.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>


#define XOR(A,B) ((A && !B) || (!A && B))


using namespace std;
using namespace reco;


class MyMonsterAnalyzer : public edm::EDAnalyzer {
public:
  explicit MyMonsterAnalyzer (const edm::ParameterSet&);
  ~MyMonsterAnalyzer ();


private:
  virtual bool L1Analyzer (const edm::Event& iEvent,
			   const edm::EventSetup& iSetup,
			   DecisionWord MyAlgoMask,
			   DecisionWord MyTrigMask,
			   const bool PrintMsg,
                           const bool PrintWords);
  virtual void analyze (const edm::Event&,
			const edm::EventSetup&);
  virtual void beginJob ();
  virtual void endJob ();


  // @@@@@@ Histograms @@@@@@
  TH1F* H_counters;
  TH1F* H_monsters;
  TH1F* H_ndigi_run1;
  TH1F* H_ndigi_run2;
  TH1F* H_ndigi_run3;


  // @@@@@@ Local variables @@@@@@
  DecisionWord MyAlgoMask;
  DecisionWord MyTechMask;
  unsigned int EvTotalRun1;
  unsigned int EvTotalRun2;
  unsigned int EvTotalRun3;
  unsigned int EvPassRun1;
  unsigned int EvPassRun2;
  unsigned int EvPassRun3;

  edm::Service<TFileService> FileService;

  // @@@@@@ Parameters from cfg file @@@@@@
  edm::InputTag inputTag; // Used for trigger
  bool PrintMsg;
  unsigned int MinHits;
};


MyMonsterAnalyzer::MyMonsterAnalyzer (const edm::ParameterSet& iConfig)
{
  inputTag = iConfig.getParameter<edm::InputTag>("inputTag");
  PrintMsg = iConfig.getParameter<bool>("PrintMsg");
  MinHits  = iConfig.getParameter<unsigned int>("MinHits");

  EvTotalRun1 = 0;
  EvTotalRun2 = 0;
  EvTotalRun3 = 0;
  EvPassRun1  = 0;
  EvPassRun2  = 0;
  EvPassRun3  = 0;

  for (int i = 0; i < 128; i++)
    MyAlgoMask.push_back(false);

  for (int i = 0; i < 64; i++)
    MyTechMask.push_back(false);
}


MyMonsterAnalyzer::~MyMonsterAnalyzer()
{
}


bool MyMonsterAnalyzer::L1Analyzer (const edm::Event& iEvent,
				    const edm::EventSetup& iSetup,
				    DecisionWord MyAlgoMask,
				    DecisionWord MyTechMask,
				    const bool PrintMsg,
				    bool PrintWords)
{
  edm::ESHandle<L1GtTriggerMenu> menuRcd;
  edm::ESHandle<L1GtTriggerMask> l1GtTmAlgo;
  edm::ESHandle<L1GtTriggerMask> l1GtTmTech;

  const L1GtTriggerMenu* L1Menu;

  bool AlgoPassed = true;
  bool TechPassed = true;

  // Extract the L1 tigger menu: algo & technical
  iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd);
  L1Menu = menuRcd.product();

  // Exract the algo & technical masks
  iSetup.get<L1GtTriggerMaskAlgoTrigRcd>().get(l1GtTmAlgo);
  iSetup.get<L1GtTriggerMaskTechTrigRcd>().get(l1GtTmTech);

  // Algo & technical trigger bits
  edm::Handle<L1GlobalTriggerReadoutRecord> gtRecord;
  iEvent.getByLabel(inputTag,gtRecord);

  DecisionWord AlgodWord = gtRecord->decisionWord();
  TechnicalTriggerWord TechWord = gtRecord->technicalTriggerWord();

  if (PrintWords == true)
    {
      edm::ESHandle<L1GtPrescaleFactors> l1GtPfAlgo;
      iSetup.get<L1GtPrescaleFactorsAlgoTrigRcd>().get(l1GtPfAlgo);
      const vector<vector<int> >* _PFAlgoTrig = &(l1GtPfAlgo.product()->gtPrescaleFactors());

      edm::ESHandle<L1GtPrescaleFactors> l1GtPfTech;
      iSetup.get<L1GtPrescaleFactorsTechTrigRcd>().get(l1GtPfTech);
      const vector<vector<int> >* _PFTechTrig = &(l1GtPfTech.product()->gtPrescaleFactors());

      const std::vector<int>* PFAlgoTrig = &((*_PFAlgoTrig).at((gtRecord->gtFdlWord()).gtPrescaleFactorIndexAlgo()));
      const std::vector<int>* PFTechTrig = &((*_PFTechTrig).at((gtRecord->gtFdlWord()).gtPrescaleFactorIndexTech()));

      cout << "\n@@@ Run number: " << iEvent.id().run() << " @@@" << endl;
      cout << "@@@ Trigger Masks and Prescale Factors for Algo-bits @@@" <<endl;
      for (vector<bool>::iterator itBit = AlgodWord.begin(); itBit != AlgodWord.end(); itBit++)
	{
	  cout << "Index algo-bits #" << itBit-AlgodWord.begin() << "\tMask: " << ((l1GtTmAlgo.product()->gtTriggerMask()[itBit-AlgodWord.begin()] & (1 << 0)) ? '1' : '0');
	  cout << "\tPrescale: " << (*PFAlgoTrig)[itBit-AlgodWord.begin()] << endl;
	}

      cout << "@@@ Trigger Masks and Prescale Factors for Tech-bits @@@" <<endl;
      for (vector<bool>::const_iterator itBit = TechWord.begin(); itBit != TechWord.end(); itBit++)
	{
	  cout << "Index tech-bits #" << itBit-TechWord.begin() << "\tMask: " << ((l1GtTmTech.product()->gtTriggerMask()[itBit-TechWord.begin()] & (1 << 0)) ? '1' : '0');
	  cout << "\tPrescale: " <<  (*PFTechTrig)[itBit-TechWord.begin()] << endl;
	}
    }
  else
    {
      // Apply the mask to the algo bits
      for (vector<bool>::iterator itBit = AlgodWord.begin(); itBit != AlgodWord.end(); itBit++)
	if ((l1GtTmAlgo.product()->gtTriggerMask()[itBit-AlgodWord.begin()] & (1 << 0)) == true) *itBit = false;
      
      // Loop over the algo bits
      for (CItAlgo itAlgo = L1Menu->gtAlgorithmMap().begin(); itAlgo != L1Menu->gtAlgorithmMap().end(); itAlgo++) {
	if (L1Menu->gtAlgorithmResult(itAlgo->second.algoName(), MyAlgoMask) == true)
	  { 
	    if (L1Menu->gtAlgorithmResult(itAlgo->second.algoName(), AlgodWord) == false)
	      {
		if (PrintMsg == true) cout << "Objet didn't pass my Algo. trigger mask: " << itAlgo->second.algoName() << endl;
		AlgoPassed = false;
		break;
	      }
	    else if (L1Menu->gtAlgorithmResult(itAlgo->second.algoName(), AlgodWord) == true)
	      {
		if (PrintMsg == true) cout << "Objet passed my Algo. trigger mask: " << itAlgo->second.algoName() << endl;
	      }
	  }
      }
      
      // Apply the mask to the techinical bits
//       for (vector<bool>::iterator itBit = TechWord.begin(); itBit != TechWord.end(); itBit++)
// 	if ((l1GtTmTech.product()->gtTriggerMask()[itBit-TechWord.begin()] & (1 << 0)) == true) *itBit = false;
      
      // Loop over the technical bits
      //   for (CItAlgo itTech = L1Menu->gtTechnicalTriggerMap().begin(); itTech != L1Menu->gtTechnicalTriggerMap().end(); itTech++) {
      for (int Count = 0; Count < 64; Count++)
	{
	  if (MyTechMask[Count] == true)
	    {
	      if (TechWord[Count] == false)
		{
		  if (PrintMsg == true) cout << "Objet didn't pass my L1 Tech. trigger mask, bit: " << Count << endl;
		  TechPassed = false;
		  break;
		}
	      else if (TechWord[Count] == true)
		{
		  if (PrintMsg == true) cout << "Objet passed my L1 Tech. trigger mask, bit: " << Count << endl;
		}
	    }
	}

      // ### Temporary condition for 2009 data taking ###
//       TechPassed = false;
//       if (!(TechWord[36] || TechWord[37] || TechWord[38] || TechWord[39]) && (TechWord[40] || TechWord[41]) && !(XOR(TechWord[42],TechWord[43]))) TechPassed = true;
      // Veto Halo + Minimum Bias Selection + Veto Single Splash
//       if (!(TechWord[40] || TechWord[41]) && (TechWord[42] != TechWord[43])) TechPassed = true;
      // ### Temporary condition for 2009 data taking ###

    }
  
  return (AlgoPassed && TechPassed) ? true : false;
}


void MyMonsterAnalyzer::analyze(const edm::Event& iEvent,
				const edm::EventSetup& iSetup)
{
  edm::Handle<SiPixelRecHitCollection> rechitspixel;
  iEvent.getByLabel("siPixelRecHits",rechitspixel);

  if (((iEvent.id().run() == 124024) && (iEvent.luminosityBlock() >= 2) && (iEvent.luminosityBlock() <= 83) &&
       (iEvent.bunchCrossing() == 51)) ||
      ((iEvent.id().run() == 124120) && (iEvent.luminosityBlock() >= 1) && (iEvent.luminosityBlock() <= 59) &&
       (iEvent.bunchCrossing() == 51)) ||
      ((iEvent.id().run() == 124275) && (iEvent.luminosityBlock() >= 3) && (iEvent.luminosityBlock() <= 30) &&
       (iEvent.bunchCrossing() == 51)))
    {      
      unsigned int DigiCounter = 0;
      bool isMonster = false;

      for (SiPixelRecHitCollection::const_iterator j = rechitspixel->begin(); j != rechitspixel->end(); j++) { // Loop on detectors
	for (edmNew::DetSet<SiPixelRecHit>::const_iterator h = j->begin(); h != j->end(); h++) { // Loop on hits on detector

	  DigiCounter = DigiCounter + h->cluster()->size();	  

	  if ((DigiCounter >= MinHits) && (L1Analyzer (iEvent, iSetup, MyAlgoMask, MyTechMask, PrintMsg, false) == true)) isMonster = true;
	}
      }

      // #####################################
      // # Record the monster events per run #
      // #####################################
      if ((isMonster == true) && (iEvent.id().run() == 124024)) EvPassRun1++;
      
      if ((isMonster == true) && (iEvent.id().run() == 124120)) EvPassRun2++;
      
      if ((isMonster == true) && (iEvent.id().run() == 124275)) EvPassRun3++;

      // #############################################
      // # Record the total number of events per run #
      // #############################################
      if (iEvent.id().run() == 124024) { EvTotalRun1++; H_ndigi_run1->Fill(DigiCounter); }
      
      if (iEvent.id().run() == 124120) { EvTotalRun2++; H_ndigi_run2->Fill(DigiCounter); }
      
      if (iEvent.id().run() == 124275) { EvTotalRun3++; H_ndigi_run3->Fill(DigiCounter); }
    }
}


void MyMonsterAnalyzer::beginJob ()
{
  H_counters = FileService->make<TH1F>("Histo counters","Histo counters", 6, 0., 6.);
  H_counters->SetXTitle("");
  H_counters->SetYTitle("Entries [#]");
  H_counters->GetXaxis()->SetBinLabel(1,"TotalRun0.9TeV");
  H_counters->GetXaxis()->SetBinLabel(2,"PassRun0.9TeV");
  H_counters->GetXaxis()->SetBinLabel(3,"TotalRun2.3TeV");
  H_counters->GetXaxis()->SetBinLabel(4,"PassRun2.3TeV");
  H_counters->GetXaxis()->SetBinLabel(5,"TotalRun2.3TeV_B*");
  H_counters->GetXaxis()->SetBinLabel(6,"PassRun2.3TeV_B*");

  H_monsters = FileService->make<TH1F>("Monster events","Monster events", 3, 0., 3.);
  H_monsters->SetXTitle("");
  H_monsters->SetYTitle("Norm. events");
  H_monsters->GetXaxis()->SetBinLabel(1,"0.9TeV");
  H_monsters->GetXaxis()->SetBinLabel(2,"2.3TeV");
  H_monsters->GetXaxis()->SetBinLabel(3,"2.3TeV_B*");

  H_ndigi_run1 = FileService->make<TH1F>("N digis run1","Digi distribution per event run1", 200, 0., 30000.);
  H_ndigi_run1->SetXTitle("Digis [#]");
  H_ndigi_run1->SetYTitle("Entries [#]");

  H_ndigi_run2 = FileService->make<TH1F>("N digis run2","Digi distribution per event run2", 200, 0., 30000.);
  H_ndigi_run2->SetXTitle("Digis [#]");
  H_ndigi_run2->SetYTitle("Entries [#]");

  H_ndigi_run3 = FileService->make<TH1F>("N digis run3","Digi distribution per event run3", 200, 0., 30000.);
  H_ndigi_run3->SetXTitle("Digis [#]");
  H_ndigi_run3->SetYTitle("Entries [#]");
}


void MyMonsterAnalyzer::endJob()
{
  cout << "\n@@@ Counters @@@" << endl;
  cout << "Total events Run1: " << EvTotalRun1 << endl;
  cout << "Total passed Run1: " << EvPassRun1 << endl;
  cout << "Total events Run2: " << EvTotalRun2 << endl;
  cout << "Total passed Run2: " << EvPassRun2 << endl;
  cout << "Total events Run3: " << EvTotalRun3 << endl;
  cout << "Total passed Run3: " << EvPassRun3 << endl;
  cout << "@@@@@@@@@@@@@@@@" << endl;

  H_counters->Fill("TotalRun0.9TeV",EvTotalRun1);
  H_counters->Fill("PassRun0.9TeV",EvPassRun1);
  H_counters->Fill("TotalRun2.3TeV",EvTotalRun2);
  H_counters->Fill("PassRun2.3TeV",EvPassRun2);
  H_counters->Fill("TotalRun2.3TeV_B*",EvTotalRun3);
  H_counters->Fill("PassRun2.3TeV_B*",EvPassRun3);
  
  double NormFactorRun1 = 1.614; // Protons per Bunch: E=900GeV; Lumi 10.8E29; Beam size: 240um; LumiBlock 2-83; Bunch 4x4; Colliding: 51, 151, 2824; crossection 52.41 mb
  double NormFactorRun2 = 8.508; // Protons per Bunch: E=2.3TeV; Lumi 3.77E29; Beam size: 115um; LumiBlock 1-59; Bunch 2x2; Colliding: 51; crossection 60.53 mb
  double NormFactorRun3 = 4.508; // Protons per Bunch: E=2.3TeV; Lumi 1.45E29; Beam size: 145um; LumiBlock 3-63; Bunch 4x4; Colliding: 51; crossection 60.53 mb

  H_monsters->Fill("0.9TeV", (double)EvPassRun1 / NormFactorRun1);
  H_monsters->Fill("2.3TeV", (double)EvPassRun2 / NormFactorRun2);
  H_monsters->Fill("2.3TeV_B*", (double)EvPassRun3 / NormFactorRun3);
}


// Define this as a plug-in
DEFINE_FWK_MODULE(MyMonsterAnalyzer);

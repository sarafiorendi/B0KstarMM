// -*- C++ -*-
//
// Package:    MyPixAnalyzer
// Class:      MyPixAnalyzer
// 
/**\class MyPixAnalyzer MyPixAnalyzer.cc MyPixAnalysis/MyPixAnalyzer/src/MyPixAnalyzer.cc

   Description: <one line class summary>
   Study of the pixel stand alone tracking and vertexing

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Mauro Dinardo
//         Created:  Thu Dec  3 10:33:18 CET 2009
// $Id: MyPixAnalyzer.cc,v 1.6 2010/08/28 15:47:24 dinardo Exp $


#include "MyPixAnalysis/MyPixAnalyzer/interface/MyPixAnalyzer.h"


MyPixAnalyzer::MyPixAnalyzer (const edm::ParameterSet& iConfig)
{
  // #####################
  // # Global parameters #
  // #####################
  PrintMsg         = iConfig.getParameter<bool>("PrintMsg");
  IsMC             = iConfig.getParameter<bool>("IsMC");
  IsTrkPart        = iConfig.getParameter<bool>("IsTrkPart");
  AnalysisType     = iConfig.getParameter<string>("AnalysisType");
  MinTkTracks      = iConfig.getParameter<unsigned int>("MinTkTracks");
  MaxTkTracks      = iConfig.getParameter<unsigned int>("MaxTkTracks");
  MinHitsMatch     = iConfig.getParameter<unsigned int>("MinHitsMatch");
  TkXYVxErrCorr    = iConfig.getParameter<double>("TkXYVxErrCorr"); 
  TkZVxErrCorr     = iConfig.getParameter<double>("TkZVxErrCorr"); 

  // ####################
  // # Track parameters #
  // ####################
  MaxEtaTkTrk      = iConfig.getParameter<double>("MaxEtaTkTrk");
  MaxChi2PxTrk     = iConfig.getParameter<double>("MaxChi2PxTrk");
  MaxChi2TkTrk     = iConfig.getParameter<double>("MaxChi2TkTrk");
  RangePt          = iConfig.getParameter<double>("RangePt");
  PtStep           = iConfig.getParameter<double>("PtStep");
  RangeEta         = iConfig.getParameter<double>("RangeEta");
  EtaStep          = iConfig.getParameter<double>("EtaStep");
  RangePhi         = iConfig.getParameter<double>("RangePhi");
  PhiStep          = iConfig.getParameter<double>("PhiStep");
  RangeChi2        = iConfig.getParameter<double>("RangeChi2");
  Chi2Step         = iConfig.getParameter<double>("Chi2Step");
  TrBound          = iConfig.getParameter<double>("TrBound");
  TzBound          = iConfig.getParameter<double>("TzBound");
  MinValidHitsPx   = iConfig.getParameter<unsigned int>("MinValidHitsPx");
  MinValidHitsTk   = iConfig.getParameter<unsigned int>("MinValidHitsTk");
  MinTrkVxDoF      = iConfig.getParameter<double>("MinTrkVxDoF");
  MinTrkVxWgt      = iConfig.getParameter<double>("MinTrkVxWgt");
  MinPtTk          = iConfig.getParameter<double>("MinPtTk");

  // #####################
  // # Vertex parameters #
  // #####################
  VxInputTag       = iConfig.getParameter<edm::InputTag>("VxInputTag");
  MaxEtaVxTrk      = iConfig.getParameter<double>("MaxEtaVxTrk");
  MaxChi2VxTrk     = iConfig.getParameter<double>("MaxChi2VxTrk");
  VrBound          = iConfig.getParameter<double>("VrBound");
  VzBound          = iConfig.getParameter<double>("VzBound");
  MinVxDoF         = iConfig.getParameter<double>("MinVxDoF");
  MinVxTrkMatch    = iConfig.getParameter<unsigned int>("MinVxTrkMatch");
  PxVxErrCorr      = iConfig.getParameter<double>("PxVxErrCorr"); 
  MinPtVx          = iConfig.getParameter<double>("MinPtVx");

  if (IsTrkPart == true) ParSetTrackAss  = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorByHitsPSet");

  cout << "\n@@@@@@ CONFIGURATION PARAMETERS @@@@@@" << endl;

  cout << "@@@@@@ GENERAL PARAMETERS @@@@@@" << endl;
  cout << "PrintMsg --> " << PrintMsg << endl;
  cout << "IsMC --> " << IsMC << endl;
  cout << "IsTrkPart --> " << IsTrkPart << endl;
  cout << "AnalysisType --> " << AnalysisType << endl;
  cout << "MinTkTracks --> " << MinTkTracks << endl;
  cout << "MaxTkTracks --> " << MaxTkTracks << endl;
  cout << "MinHitsMatch --> " << MinHitsMatch << endl;
  cout << "TkXYVxErrCorr --> " << TkXYVxErrCorr << endl;
  cout << "TkZVxErrCorr --> " << TkZVxErrCorr << endl;

  cout << "@@@@@@ TRACK PARAMETERS @@@@@@" << endl;
  cout << "MaxEtaTkTrk --> " << MaxEtaTkTrk << endl;
  cout << "MaxChi2PxTrk --> " << MaxChi2PxTrk << endl;
  cout << "MaxChi2TkTrk --> " << MaxChi2TkTrk << endl;
  cout << "RangePt --> " << RangePt << endl;
  cout << "PtStep --> " << PtStep << endl;
  cout << "RangeEta --> " << RangeEta << endl;
  cout << "EtaStep --> " << EtaStep << endl;
  cout << "RangePhi --> " << RangePhi << endl;
  cout << "PhiStep --> " << PhiStep << endl;
  cout << "RangeChi2 --> " << RangeChi2 << endl;
  cout << "Chi2Step --> " << Chi2Step << endl;
  cout << "TrBound --> " << TrBound << endl;
  cout << "TzBound --> " << TzBound << endl;
  cout << "MinValidHitsPx --> " << MinValidHitsPx << endl;
  cout << "MinValidHitsTk --> " << MinValidHitsTk << endl;
  cout << "MinTrkVxDoF --> " << MinTrkVxDoF << endl;
  cout << "MinTrkVxWgt --> " << MinTrkVxWgt << endl;
  cout << "MinPtTk --> " << MinPtTk << endl;

  cout << "@@@@@@ VERTEX PARAMETERS @@@@@@" << endl;
  cout << "VxInputTag --> " << VxInputTag << endl;
  cout << "MaxEtaVxTrk --> " << MaxEtaVxTrk << endl;
  cout << "MaxChi2VxTrk --> " << MaxChi2VxTrk << endl;
  cout << "VrBound --> " << VrBound << endl;
  cout << "VzBound --> " << VzBound << endl;
  cout << "MinVxDoF --> " << MinVxDoF << endl;
  cout << "MinVxTrkMatch --> " << MinVxTrkMatch << endl;
  cout << "PxVxErrCorr --> " << PxVxErrCorr << endl;
  cout << "MinPtVx --> " << MinPtVx << endl;

  EtaThr = 1.5;
  EtaThrBin[0] = 0;
  EtaThrBin[1] = EtaThr;
  EtaThrBin[2] = MaxEtaTkTrk;

  d0Range    = 0.6;      // Unit: (cm)
  dzRange    = 0.6;      // Unit: (cm)
  ptRange    = 12.0;     // Unit: (GeV/c)
  d0errRange = 10000.;   // Unit: (um)^2
  dzerrRange = 10000.;   // Unit: (um)^2
  pterrRange = 50000.;   // Unit: (MeV/c)^2

  pullsRange = 10.;

  xRange     = 3000.;    // Unit: (um)
  yRange     = 3000.;    // Unit: (um)
  zRange     = 6000.;    // Unit: (um)
  xerrRange  = 100000.;  // Unit: (um)^2
  yerrRange  = 100000.;  // Unit: (um)^2
  zerrRange  = 100000.;  // Unit: (um)^2

  EvDiscardL1    = 0;
  EvDiscardHLT   = 0;
  EvPassL1       = 0;
  EvPassHLT      = 0;
  EvPassTotal    = 0;
  EvDiscardTotal = 0;
  EvTotal        = 0;

  NpixelLayers = 3;

  for (int i = 0; i < 128; i++)
    MyAlgoMask.push_back(false);

  for (int i = 0; i < 64; i++)
    MyTechMask.push_back(false);
  
  MyHLTMask.push_back("HLT_L1MuOpen");
}


MyPixAnalyzer::~MyPixAnalyzer()
{
}


unsigned int MyPixAnalyzer::CountSimChargedTracks(TrackingVertexRef VxTP)
{
  unsigned int ChrgTrkCounter = 0;
  
  for (TrackingParticleRefVector::iterator itTkTP = VxTP->daughterTracks_begin(); itTkTP != VxTP->daughterTracks_end(); itTkTP++)
    if (((*itTkTP)->status() == 1) && (abs((*itTkTP)->charge()) >= 1)) ChrgTrkCounter++;
  
  return ChrgTrkCounter;
}


// ##############################
// # @@@ METHODS for TRACKs @@@ #
// ##############################

bool MyPixAnalyzer::FindNpixelHits(unsigned int nHits,
				   const Track* TkTrack)
{
  unsigned int Counter = 0;
  unsigned int Layer   = 0;
  unsigned int Disk    = 0;

  for (TrackingRecHitRefVector::iterator itRecHit = TkTrack->recHitsBegin(); itRecHit != TkTrack->recHitsEnd(); itRecHit++) {
    if (((*itRecHit)->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) &&
	(((PXBDetId)((*itRecHit)->geographicalId())).layer() > Layer) &&
	((*itRecHit)->isValid() == true))
      {
	Layer = ((PXBDetId)((*itRecHit)->geographicalId())).layer();
	Counter++;
      }
    else if (((*itRecHit)->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap) &&
	     (((PXFDetId)((*itRecHit)->geographicalId())).disk() > Disk) &&
	     ((*itRecHit)->isValid() == true))
      {
	Disk = ((PXFDetId)((*itRecHit)->geographicalId())).disk();
	Counter++;
      }
  }
  
  if (Counter >= nHits) return true;
  return false;
}


bool MyPixAnalyzer::FindNpixelHits(unsigned int nHits,
				   TrackingParticleRef SimTrk)
{
  unsigned int Counter = 0;
  unsigned int Layer   = 0;
  unsigned int Disk    = 0;

  for (vector<PSimHit>::const_iterator itSimHit = SimTrk->pSimHit_begin(); itSimHit != SimTrk->pSimHit_end(); itSimHit++) {
    if ((((DetId)(itSimHit->detUnitId())).subdetId() == (int)PixelSubdetector::PixelBarrel) &&
	(((PXBDetId)(itSimHit->detUnitId())).layer() > Layer))
      {
	Layer = ((PXBDetId)(itSimHit->detUnitId())).layer();
	Counter++;
      }
    else if ((((DetId)(itSimHit->detUnitId())).subdetId() == (int)PixelSubdetector::PixelEndcap) &&
	     (((PXFDetId)(itSimHit->detUnitId())).disk() > Disk))
      {
	Disk = ((PXFDetId)(itSimHit->detUnitId())).disk();
	Counter++;
      }
  }
  
  if (Counter >= nHits) return true;
  return false;      
}


void MyPixAnalyzer::AddDataMap (string mapname,
				double pt,
				double eta,
				double data,
				int bins,
				double beginval,
				double endval)
{
  stringstream PtCutStr;					
  stringstream EtaCutStr;					
  PtCutStr.precision(3);
  EtaCutStr.precision(3);

  double PtCut;
  for (PtCut = PtStep; rint(PtCut*100.) <= rint(RangePt*100.); PtCut += PtStep)
    if (pt <= PtCut) break;
  PtCutStr << PtCut;

  double EtaCut;
  for (EtaCut = EtaStep; rint(EtaCut*100.) <= rint(RangeEta*100.); EtaCut += EtaStep)
    if (eta <= EtaCut) break;
  EtaCutStr << EtaCut;

  // ####################################################
  // # pt and eta bins are chosen in such a way that:   #
  // # PtCut - PtStep < pt <= PtCut                     #
  // # EtaCut - EtaStep < eta <= EtaCut                 #
  // # The bin RangePt+PtStep contains pt > RangePt     #
  // # The bin RangeEta+EtaStep contains eta > RangeEta #
  // ####################################################

  if (mapname.compare("d0pulls") == 0) // @@@ d0pulls MAP @@@
    {
      if ((MapFit_d0_pull.find(PtCutStr.str().c_str()) == MapFit_d0_pull.end()) || (MapFit_d0_pull[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) == MapFit_d0_pull[PtCutStr.str().c_str()].end()))
	MapFit_d0_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = NULL;
      
      if (MapFit_d0_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] == NULL)
	{
	  stringstream HistoName;					
	  HistoName << "d0_pulls_pt_" << PtCut << "_eta_" << EtaCut;
	  MapFit_d0_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = SubDirPixTkFit->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), bins, beginval, endval);
	  MapFit_d0_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetXTitle("(Pix-Vx)_{d0} pulls");
	  MapFit_d0_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetYTitle("Entries (#)");
	}
      
      MapFit_d0_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->Fill(data);
    }
  else if (mapname.compare("d0err") == 0) // @@@ d0err MAP @@@
    {
      if ((MapFit_d0err.find(PtCutStr.str().c_str()) == MapFit_d0err.end()) || (MapFit_d0err[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) == MapFit_d0err[PtCutStr.str().c_str()].end()))
	MapFit_d0err[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = NULL;
      
      if (MapFit_d0err[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] == NULL)
	{
	  stringstream HistoName;					
	  HistoName << "d0err_pt_" << PtCut << "_eta_" << EtaCut;
	  MapFit_d0err[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = SubDirPixTkFit->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), bins, beginval, endval);
	  MapFit_d0err[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetXTitle("Trk \\sigma^{2}(d_{0}) (\\mum)^{2}");
	  MapFit_d0err[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetYTitle("Entries (#)");
	}
      
      MapFit_d0err[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->Fill(data);
    }
  else if (mapname.compare("d0") == 0) // @@@ d0 MAP @@@
    {
      if ((MapFit_d0.find(PtCutStr.str().c_str()) == MapFit_d0.end()) || (MapFit_d0[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) == MapFit_d0[PtCutStr.str().c_str()].end()))
	MapFit_d0[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = NULL;

      if (MapFit_d0[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] == NULL)
	{
	  stringstream HistoName;					
	  HistoName << "d0_pt_" << PtCut << "_eta_" << EtaCut;
	  MapFit_d0[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = SubDirPixTkFit->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), bins, beginval, endval);
	  MapFit_d0[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetXTitle("(Pix-Vx)_{d0} (cm)");
	  MapFit_d0[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetYTitle("Entries (#)");
	}
      
      MapFit_d0[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->Fill(data);
    }
  else if (mapname.compare("dzpulls") == 0) // @@@ dzpulls MAP @@@
    {
      if ((MapFit_dz_pull.find(PtCutStr.str().c_str()) == MapFit_dz_pull.end()) || (MapFit_dz_pull[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) == MapFit_dz_pull[PtCutStr.str().c_str()].end()))
	MapFit_dz_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = NULL;

      if (MapFit_dz_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] == NULL)
	{
	  stringstream HistoName;					
	  HistoName << "dz_pulls_pt_" << PtCut << "_eta_" << EtaCut;
	  MapFit_dz_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = SubDirPixTkFit->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), bins, beginval, endval);
	  MapFit_dz_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetXTitle("(Pix-Vx)_{dz} pulls");
	  MapFit_dz_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetYTitle("Entries (#)");
	}
                  
      MapFit_dz_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->Fill(data);
    }
  else if (mapname.compare("dzerr") == 0) // @@@ dzerr MAP @@@
    {
      if ((MapFit_dzerr.find(PtCutStr.str().c_str()) == MapFit_dzerr.end()) || (MapFit_dzerr[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) == MapFit_dzerr[PtCutStr.str().c_str()].end()))
	MapFit_dzerr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = NULL;
      
      if (MapFit_dzerr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] == NULL)
	{
	  stringstream HistoName;					
	  HistoName << "dzerr_pt_" << PtCut << "_eta_" << EtaCut;
	  MapFit_dzerr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = SubDirPixTkFit->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), bins, beginval, endval);
	  MapFit_dzerr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetXTitle("Trk \\sigma^{2}(d_{z}) (\\mum)^{2}");
	  MapFit_dzerr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetYTitle("Entries (#)");
	}
      
      MapFit_dzerr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->Fill(data);
    }
  else if (mapname.compare("dz") == 0) // @@@ dz MAP @@@
    {
      if ((MapFit_dz.find(PtCutStr.str().c_str()) == MapFit_dz.end()) || (MapFit_dz[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) == MapFit_dz[PtCutStr.str().c_str()].end()))
	MapFit_dz[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = NULL;
      
      if (MapFit_dz[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] == NULL)
	{
	  stringstream HistoName;					
	  HistoName << "dz_pt_" << PtCut << "_eta_" << EtaCut;
	  MapFit_dz[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = SubDirPixTkFit->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), bins, beginval, endval);
	  MapFit_dz[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetXTitle("(Pix-Vx)_{dz} (cm)");
	  MapFit_dz[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetYTitle("Entries (#)");
	}
                  
      MapFit_dz[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->Fill(data);
    }
  else if (mapname.compare("ptpulls") == 0) // @@@ ptpulls MAP @@@
    {
      if ((MapFit_pt_pull.find(PtCutStr.str().c_str()) == MapFit_pt_pull.end()) || (MapFit_pt_pull[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) == MapFit_pt_pull[PtCutStr.str().c_str()].end()))
	MapFit_pt_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = NULL;

      if (MapFit_pt_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] == NULL)
	{
	  stringstream HistoName;					
	  HistoName << "pt_pulls_pt_" << PtCut << "_eta_" << EtaCut;
	  MapFit_pt_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = SubDirPixTkFit->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), bins, beginval, endval);
	  MapFit_pt_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetXTitle("(Pix-Trk)_{pt} pulls");
	  MapFit_pt_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetYTitle("Entries (#)");
	}

      MapFit_pt_pull[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->Fill(data);
    }
  else if (mapname.compare("pterr") == 0) // @@@ pterr MAP @@@
    {
      if ((MapFit_pterr.find(PtCutStr.str().c_str()) == MapFit_pterr.end()) || (MapFit_pterr[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) == MapFit_pterr[PtCutStr.str().c_str()].end()))
	MapFit_pterr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = NULL;
      
      if (MapFit_pterr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] == NULL)
	{
	  stringstream HistoName;					
	  HistoName << "pterr_pt_" << PtCut << "_eta_" << EtaCut;
	  MapFit_pterr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = SubDirPixTkFit->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), bins, beginval, endval);
	  MapFit_pterr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetXTitle("Trk \\sigma^{2}(p_{t}) (MeV/c)^{2}");
	  MapFit_pterr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetYTitle("Entries (#)");
	}
      
      MapFit_pterr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->Fill(data);
    }
  else if (mapname.compare("pixpterr") == 0) // @@@ pixpterr MAP @@@
    {
      if ((MapFit_pixpterr.find(PtCutStr.str().c_str()) == MapFit_pixpterr.end()) || (MapFit_pixpterr[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) == MapFit_pixpterr[PtCutStr.str().c_str()].end()))
	MapFit_pixpterr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = NULL;
      
      if (MapFit_pixpterr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] == NULL)
	{
	  stringstream HistoName;					
	  HistoName << "pixpterr_pt_" << PtCut << "_eta_" << EtaCut;
	  MapFit_pixpterr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = SubDirPixTkFit->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), bins, beginval, endval);
	  MapFit_pixpterr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetXTitle("Pix \\sigma(p_{t}) (MeV/c)");
	  MapFit_pixpterr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetYTitle("Entries (#)");
	}
      
      MapFit_pixpterr[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->Fill(data);
    }
  else if (mapname.compare("pt") == 0) // @@@ pt MAP @@@
    {
      if ((MapFit_pt.find(PtCutStr.str().c_str()) == MapFit_pt.end()) || (MapFit_pt[PtCutStr.str().c_str()].find(EtaCutStr.str().c_str()) == MapFit_pt[PtCutStr.str().c_str()].end()))
	MapFit_pt[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = NULL;

      if (MapFit_pt[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] == NULL)
	{
	  stringstream HistoName;					
	  HistoName << "pt_pt_" << PtCut << "_eta_" << EtaCut;
	  MapFit_pt[PtCutStr.str().c_str()][EtaCutStr.str().c_str()] = SubDirPixTkFit->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), bins, beginval, endval);
	  MapFit_pt[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetXTitle("(Pix-Trk)_{pt} (GeV/c)");
	  MapFit_pt[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->SetYTitle("Entries (#)");
	}

      MapFit_pt[PtCutStr.str().c_str()][EtaCutStr.str().c_str()]->Fill(data);
    }
}


void MyPixAnalyzer::FillTrackHistos (string DataType,
				     TrackParType* TrackPar)
{
  if ((DataType.compare("pixel") == 0) && (TrackPar->PxTrackValidHits >= MinValidHitsPx))
    {
      if ((TrackPar->PxTrackPt >= MinPtTk) && (TrackPar->PxTrackChi2DoF <= MaxChi2PxTrk))
	H_pix_trk_eta_phi->Fill(TrackPar->PxTrackEta, TrackPar->PxTrackPhi/PI*180.);
      if (TrackPar->PxTrackChi2DoF <= MaxChi2PxTrk) H_pix_trk_eta_pt->Fill(TrackPar->PxTrackEta, TrackPar->PxTrackPt);
      
      for (int k = 0; k < 2; k++)
	if ((fabs(TrackPar->PxTrackEta) <= EtaThrBin[k+1]) && (fabs(TrackPar->PxTrackEta) > EtaThrBin[k]))
	  for (int i = H_pix_trk_normchi2_pt_eta[k]->GetXaxis()->FindBin(TrackPar->PxTrackPt); i >= 1; i--)
	    for (int j = H_pix_trk_normchi2_pt_eta[k]->GetYaxis()->FindBin(TrackPar->PxTrackChi2DoF); j <= H_pix_trk_normchi2_pt_eta[k]->GetNbinsY(); j++)
	      H_pix_trk_normchi2_pt_eta[k]->SetBinContent(i,j, H_pix_trk_normchi2_pt_eta[k]->GetBinContent(i,j) + 1.);

      // ######################
      // # Denominator Purity #
      // ######################
      if ((fabs(TrackPar->PxTrackEta) <= MaxEtaTkTrk) && (TrackPar->PxTrackPt >= MinPtTk)) H_pix_trk_normchi2[2]->Fill(TrackPar->PxTrackChi2DoF);
      if ((fabs(TrackPar->PxTrackEta) <= MaxEtaTkTrk) && (TrackPar->PxTrackChi2DoF <= MaxChi2PxTrk)) H_pix_trk_pt[2]->Fill(TrackPar->PxTrackPt);
      if ((TrackPar->PxTrackPt >= MinPtTk) && (TrackPar->PxTrackChi2DoF <= MaxChi2PxTrk))
	{
	  H_pix_trk_theta[2]->Fill(TrackPar->PxTrackTheta);
	  H_pix_trk_phi[2]->Fill(TrackPar->PxTrackPhi/PI*180.);
	  H_pix_trk_eta[2]->Fill(TrackPar->PxTrackEta);     
	}
    }
  else if ((DataType.compare("strip") == 0) &&
	   (TrackPar->TkTrackValidHits >= MinValidHitsTk) &&
	   (fabs(TrackPar->TkOriginXY) <= TrBound) && (fabs(TrackPar->TkOriginZ) <= TzBound))
    {
      if ((TrackPar->TkTrackPt >= MinPtTk) && (TrackPar->TkTrackChi2DoF <= MaxChi2TkTrk))
	H_trk_trk_eta_phi->Fill(TrackPar->TkTrackEta, TrackPar->TkTrackPhi/PI*180.);
      if (TrackPar->TkTrackChi2DoF <= MaxChi2TkTrk) H_trk_trk_eta_pt->Fill(TrackPar->TkTrackEta, TrackPar->TkTrackPt);
      
      // ##########################
      // # Denominator Efficiency #
      // ##########################
      if ((fabs(TrackPar->TkTrackEta) <= MaxEtaTkTrk) && (TrackPar->TkTrackPt >= MinPtTk)) H_trk_trk_normchi2[2]->Fill(TrackPar->TkTrackChi2DoF);
      if ((fabs(TrackPar->TkTrackEta) <= MaxEtaTkTrk) && (TrackPar->TkTrackChi2DoF <= MaxChi2TkTrk)) H_trk_trk_pt[2]->Fill(TrackPar->TkTrackPt);
      if ((TrackPar->TkTrackPt >= MinPtTk) && (TrackPar->TkTrackChi2DoF <= MaxChi2TkTrk))
	{
	  H_trk_trk_theta[2]->Fill(TrackPar->TkTrackTheta);
	  H_trk_trk_phi[2]->Fill(TrackPar->TkTrackPhi/PI*180.);
	  H_trk_trk_eta[2]->Fill(TrackPar->TkTrackEta);
	}
      
      if (TrackPar->AlgoEff == true)
	{
	  if ((TrackPar->TkTrackPt >= MinPtTk) && (TrackPar->TkTrackChi2DoF <= MaxChi2TkTrk))
	    H_trk_trk_eta_phi_algo->Fill(TrackPar->TkTrackEta, TrackPar->TkTrackPhi/PI*180.);
	  if (TrackPar->TkTrackChi2DoF <= MaxChi2TkTrk) H_trk_trk_eta_pt_algo->Fill(TrackPar->TkTrackEta, TrackPar->TkTrackPt);

	  // ###############################
	  // # Denominator Algo-Efficiency #
	  // ###############################
	  if ((fabs(TrackPar->TkTrackEta) <= MaxEtaTkTrk) && (TrackPar->TkTrackPt >= MinPtTk)) H_trk_trk_normchi2[3]->Fill(TrackPar->TkTrackChi2DoF);
	  if ((fabs(TrackPar->TkTrackEta) <= MaxEtaTkTrk) && (TrackPar->TkTrackChi2DoF <= MaxChi2TkTrk)) H_trk_trk_pt[3]->Fill(TrackPar->TkTrackPt);	  
	  if ((TrackPar->TkTrackPt >= MinPtTk) && (TrackPar->TkTrackChi2DoF <= MaxChi2TkTrk))
	    {
	      H_trk_trk_theta[3]->Fill(TrackPar->TkTrackTheta);
	      H_trk_trk_phi[3]->Fill(TrackPar->TkTrackPhi/PI*180.);
	      H_trk_trk_eta[3]->Fill(TrackPar->TkTrackEta);
	    }
	}
  
      for (int k = 0; k < 2; k++)
	if ((fabs(TrackPar->TkTrackEta) <= EtaThrBin[k+1]) && (fabs(TrackPar->TkTrackEta) > EtaThrBin[k]))
	  for (int i = H_trk_trk_normchi2_pt_eta[k]->GetXaxis()->FindBin(TrackPar->TkTrackPt); i >= 1; i--)
	    for (int j = H_trk_trk_normchi2_pt_eta[k]->GetYaxis()->FindBin(TrackPar->TkTrackChi2DoF); j <= H_trk_trk_normchi2_pt_eta[k]->GetNbinsY(); j++)
	      H_trk_trk_normchi2_pt_eta[k]->SetBinContent(i,j, H_trk_trk_normchi2_pt_eta[k]->GetBinContent(i,j) + 1.);      
    }
  else if (DataType.compare("match") == 0)
    {
      if (PrintMsg == true)
	{
	  cout << "### Matched track ###" << endl;
	  cout << "\tNumber of valid hits: " << TrackPar->PxTrackValidHits << ", number of lost hits: " << TrackPar->PxTrackLostHits << endl;
	  cout << "\tChi2/DoF-pixel: " << TrackPar->PxTrackChi2DoF << "\tChi2/DoF-track: " << TrackPar->TkTrackChi2DoF << endl;
	  cout << "\tpt-pixel: " << TrackPar->PxTrackPt << "\tpt-track: " << TrackPar->TkTrackPt << endl;
	  cout << "\tdz-pixel: " << TrackPar->PxTrackDz << "\tdz-track:" << TrackPar->TkTrackDz << endl;
	  cout << "\td0-pixel: " << TrackPar->PxTrackD0 << "\td0-track: " << TrackPar->TkTrackD0 << endl;
	  cout << "\ttheta-pixel: " << TrackPar->PxTrackTheta << "\ttheta-track:" << TrackPar->TkTrackTheta << endl;
	  cout << "\tphi-pixel: " << TrackPar->PxTrackPhi << "\tphi-track: " << TrackPar->TkTrackPhi << endl;
	  cout << "\teta-pixel: " << TrackPar->PxTrackEta << "\teta-track: " << TrackPar->TkTrackEta << endl;
	}

      // ################
      // # PURITY PLOTS #
      // ################
      if (TrackPar->PxTrackValidHits >= MinValidHitsPx)
	{
	  if ((TrackPar->PxTrackPt >= MinPtTk) && (TrackPar->PxTrackChi2DoF <= MaxChi2PxTrk))
	    H_trked_eta_phi_pur->Fill(TrackPar->PxTrackEta, TrackPar->PxTrackPhi/PI*180.);
	  if (TrackPar->PxTrackChi2DoF <= MaxChi2PxTrk) H_trked_eta_pt_pur->Fill(TrackPar->PxTrackEta, TrackPar->PxTrackPt);
	  
	  for (int k = 0; k < 2; k++)
	    if ((fabs(TrackPar->PxTrackEta) <= EtaThrBin[k+1]) && (fabs(TrackPar->PxTrackEta) > EtaThrBin[k]))
	      for (int i = H_trked_normchi2_pt_pur_eta[k]->GetXaxis()->FindBin(TrackPar->PxTrackPt); i >= 1; i--)
		for (int j = H_trked_normchi2_pt_pur_eta[k]->GetYaxis()->FindBin(TrackPar->PxTrackChi2DoF); j <= H_trked_normchi2_pt_pur_eta[k]->GetNbinsY(); j++)
		  H_trked_normchi2_pt_pur_eta[k]->SetBinContent(i,j, H_trked_normchi2_pt_pur_eta[k]->GetBinContent(i,j) + 1.);	  
	}

      // ####################
      // # EFFICIENCY PLOTS #
      // ####################
      if ((TrackPar->TkTrackValidHits >= MinValidHitsTk) &&
	  (fabs(TrackPar->TkOriginXY) <= TrBound) && (fabs(TrackPar->TkOriginZ) <= TzBound))
	{
	  if ((TrackPar->TkTrackPt >= MinPtTk) && (TrackPar->TkTrackChi2DoF <= MaxChi2TkTrk))
	    H_trked_eta_phi_eff->Fill(TrackPar->TkTrackEta, TrackPar->TkTrackPhi/PI*180.);
	  if (TrackPar->TkTrackChi2DoF <= MaxChi2TkTrk) H_trked_eta_pt_eff->Fill(TrackPar->TkTrackEta, TrackPar->TkTrackPt);
	  
	  for (int k = 0; k < 2; k++)
	    if ((fabs(TrackPar->TkTrackEta) <= EtaThrBin[k+1]) && (fabs(TrackPar->TkTrackEta) > EtaThrBin[k]))
	      for (int i = H_trked_normchi2_pt_eff_eta[k]->GetXaxis()->FindBin(TrackPar->TkTrackPt); i >= 1; i--)
		for (int j = H_trked_normchi2_pt_eff_eta[k]->GetYaxis()->FindBin(TrackPar->TkTrackChi2DoF); j <= H_trked_normchi2_pt_eff_eta[k]->GetNbinsY(); j++)
		  H_trked_normchi2_pt_eff_eta[k]->SetBinContent(i,j, H_trked_normchi2_pt_eff_eta[k]->GetBinContent(i,j) + 1.);
	}

      if ((TrackPar->PxTrackValidHits >= MinValidHitsPx) && (TrackPar->PxTrackChi2DoF <= MaxChi2PxTrk) &&
	  (TrackPar->TkTrackValidHits >= MinValidHitsTk) && (TrackPar->TkTrackChi2DoF <= MaxChi2TkTrk))
	{
	  if (TrackPar->BelongVx == true)
	    {
	      AddDataMap("d0", TrackPar->TkTrackPt, fabs(TrackPar->TkTrackEta), TrackPar->PxTrackD0, 100, -d0Range/2., d0Range/2.);
	      AddDataMap("d0pulls", TrackPar->TkTrackPt, fabs(TrackPar->TkTrackEta), TrackPar->PxTrackD0/sqrt(powf(TrackPar->PxTrackD0Err,2.) + powf(TrackPar->TkTrackD0Err,2.)), 100, -pullsRange/2., pullsRange/2.);
	      AddDataMap("d0err", TrackPar->TkTrackPt, fabs(TrackPar->TkTrackEta), powf(TrackPar->TkTrackD0Err,2.)*100000000., 200, 0., d0errRange);

	      AddDataMap("dz", TrackPar->TkTrackPt, fabs(TrackPar->TkTrackEta), TrackPar->PxTrackDz, 100, -dzRange/2., dzRange/2.);
	      AddDataMap("dzpulls", TrackPar->TkTrackPt, fabs(TrackPar->TkTrackEta), TrackPar->PxTrackDz/sqrt(powf(TrackPar->PxTrackDzErr,2.) + powf(TrackPar->TkTrackDzErr,2.)), 100, -pullsRange/2., pullsRange/2.);
	      AddDataMap("dzerr", TrackPar->TkTrackPt, fabs(TrackPar->TkTrackEta), powf(TrackPar->TkTrackDzErr,2.)*100000000., 200, 0., dzerrRange);
	      
	      // #########
	      // # STRIP #
	      // #########
	      H_trk_trk_dz[0]->Fill(TrackPar->TkTrackDz);
	      H_trk_trk_d0[0]->Fill(TrackPar->TkTrackD0);
	      H_trk_trk_dz[1]->Fill(TrackPar->TkTrackDz);
	      H_trk_trk_d0[1]->Fill(TrackPar->TkTrackD0);
		      
	      // #########
	      // # PIXEL #
	      // #########
	      H_pix_trk_dz[0]->Fill(TrackPar->PxTrackDz);
	      H_pix_trk_d0[0]->Fill(TrackPar->PxTrackD0);
	      H_pix_trk_dz[1]->Fill(TrackPar->PxTrackDz);
	      H_pix_trk_d0[1]->Fill(TrackPar->PxTrackD0);
	    }

	  AddDataMap("pt", TrackPar->TkTrackPt, fabs(TrackPar->TkTrackEta), TrackPar->PxTrackPt - TrackPar->TkTrackPt, 600, -ptRange/2., ptRange/2.);
	  AddDataMap("ptpulls", TrackPar->TkTrackPt, fabs(TrackPar->TkTrackEta), (TrackPar->PxTrackPt - TrackPar->TkTrackPt)/sqrt(powf(TrackPar->PxTrackPtErr,2.) + powf(TrackPar->TkTrackPtErr,2.)), 100, -pullsRange/2., pullsRange/2.);
	  AddDataMap("pterr", TrackPar->TkTrackPt, fabs(TrackPar->TkTrackEta), powf(TrackPar->TkTrackPtErr,2.)*1000000., 5000, 0., pterrRange);

	  AddDataMap("pixpterr", TrackPar->TkTrackPt, fabs(TrackPar->TkTrackEta), TrackPar->PxTrackPtErr*1000., 5000, 0., pterrRange/10.);
	}

      // ############################################
      // # Numerator Efficiency and Algo-Efficiency #
      // ############################################
      if ((fabs(TrackPar->TkTrackEta) <= MaxEtaTkTrk) && 
	  (TrackPar->TkTrackValidHits >= MinValidHitsTk) &&
	  (fabs(TrackPar->TkOriginXY) <= TrBound) && (fabs(TrackPar->TkOriginZ) <= TzBound))
	{
	  if ((fabs(TrackPar->TkTrackEta) <= MaxEtaTkTrk) && (TrackPar->TkTrackPt >= MinPtTk))
	    {
	      // #########
	      // # STRIP #
	      // #########
	      H_trk_trk_normchi2[0]->Fill(TrackPar->TkTrackChi2DoF);

	      // #########
	      // # PIXEL #
	      // #########
	      H_pix_trk_normchi2[0]->Fill(TrackPar->PxTrackChi2DoF);
	    }

	  if ((fabs(TrackPar->TkTrackEta) <= MaxEtaTkTrk) && (TrackPar->TkTrackChi2DoF <= MaxChi2TkTrk))
	    {
	      // #########
	      // # STRIP #
	      // #########
	      H_trk_trk_pt[0]->Fill(TrackPar->TkTrackPt);
		  
	      // #########
	      // # PIXEL #
	      // #########
	      H_pix_trk_pt[0]->Fill(TrackPar->PxTrackPt);
	    }

	  if ((TrackPar->TkTrackPt >= MinPtTk) && (TrackPar->TkTrackChi2DoF <= MaxChi2TkTrk))
	    {
	      // #########
	      // # STRIP #
	      // #########
	      H_trk_trk_theta[0]->Fill(TrackPar->TkTrackTheta);
	      H_trk_trk_phi[0]->Fill(TrackPar->TkTrackPhi/PI*180.);
	      H_trk_trk_eta[0]->Fill(TrackPar->TkTrackEta);

	      // #########
	      // # PIXEL #
	      // #########
	      H_pix_trk_theta[0]->Fill(TrackPar->PxTrackTheta);
	      H_pix_trk_phi[0]->Fill(TrackPar->PxTrackPhi/PI*180.);
	      H_pix_trk_eta[0]->Fill(TrackPar->PxTrackEta);     
	    }
	}

      // ####################
      // # Numerator Purity #
      // ####################
      if (TrackPar->PxTrackValidHits >= MinValidHitsPx)
	{
	  if ((fabs(TrackPar->PxTrackEta) <= MaxEtaTkTrk) && (TrackPar->PxTrackPt >= MinPtTk))
	    {
	      // #########
	      // # STRIP #
	      // #########
	      H_trk_trk_normchi2[1]->Fill(TrackPar->TkTrackChi2DoF);

	      // #########
	      // # PIXEL #
	      // #########
	      H_pix_trk_normchi2[1]->Fill(TrackPar->PxTrackChi2DoF);
	    }

	  if ((fabs(TrackPar->PxTrackEta) <= MaxEtaTkTrk) && (TrackPar->PxTrackChi2DoF <= MaxChi2PxTrk))
	    {
	      // #########
	      // # STRIP #
	      // #########
	      H_trk_trk_pt[1]->Fill(TrackPar->TkTrackPt);
		  
	      // #########
	      // # PIXEL #
	      // #########
	      H_pix_trk_pt[1]->Fill(TrackPar->PxTrackPt);
	    }

	  if ((TrackPar->PxTrackPt >= MinPtTk) && (TrackPar->PxTrackChi2DoF <= MaxChi2PxTrk))
	    {
	      // #########
	      // # STRIP #
	      // #########
	      H_trk_trk_theta[1]->Fill(TrackPar->TkTrackTheta);
	      H_trk_trk_phi[1]->Fill(TrackPar->TkTrackPhi/PI*180.);
	      H_trk_trk_eta[1]->Fill(TrackPar->TkTrackEta);
	      
	      // #########
	      // # PIXEL #
	      // #########
	      H_pix_trk_theta[1]->Fill(TrackPar->PxTrackTheta);
	      H_pix_trk_phi[1]->Fill(TrackPar->PxTrackPhi/PI*180.);
	      H_pix_trk_eta[1]->Fill(TrackPar->PxTrackEta);     
	    }
	}	      
    }
}


void MyPixAnalyzer::PixelTracksMC (const edm::Event& iEvent,
				   const edm::EventSetup& iSetup)
{
  edm::Handle<edm::View<Track> > TkTracks;
  iEvent.getByLabel("generalTracks",TkTracks);

  edm::Handle<edm::View<Track> > PxTracks;
  iEvent.getByLabel("pixelTracks",PxTracks);

  edm::RefToBase<Track> RecoTrk;
  
  edm::Handle<TrackingParticleCollection> TPCollection;
  iEvent.getByLabel(edm::InputTag("mergedtruth","MergedTrackTruth"),TPCollection);
  const TrackingParticleCollection TPC = *(TPCollection.product());
  RecoToSimCollection PxTkReco2Sim = theTkAssociator->associateRecoToSim(PxTracks, TPCollection, &iEvent);
  SimToRecoCollection TkTkSim2Reco = theTkAssociator->associateSimToReco(TkTracks, TPCollection, &iEvent);
  vector<pair<edm::RefToBase<reco::Track>, double> > TrkReco;
  vector<pair<TrackingParticleRef, double> > TrkParticle;

  edm::ESHandle<ParametersDefinerForTP> parametersDefinerTP; 
  iSetup.get<TrackAssociatorRecord>().get("LhcParametersDefinerForTP",parametersDefinerTP);  

  TrackingParticleRef SimTrk;
  TrackingVertexRef SimVx;

  TrackParType TrackPar;

  map<unsigned int, unsigned int> PxTrkKeys_Counter;
  map<unsigned int, double> PxTrkKeys_NormChi2;
  map<unsigned int, double> PxTrkKeys_Eta;
  map<unsigned int, double> PxTrkKeys_Pt;

  double Simpt;
  double Simdz;
  double Simdxy;
  double Simtheta;
  double Simphi;
  double Simeta;
  
  ParticleBase::Vector MomentumTP;
  ParticleBase::Point VertexTP;
  ParticleBase::Point VxSimTk;

  // ##############################
  // # Loop over the pixel tracks #
  // ##############################
  for (edm::View<Track>::const_iterator itPxTrack = PxTracks->begin(); itPxTrack != PxTracks->end(); itPxTrack++)
    {       	
      TrackPar.PxTrackValidHits = itPxTrack->numberOfValidHits();
      TrackPar.PxTrackChi2DoF   = itPxTrack->normalizedChi2();
      TrackPar.PxTrackPt        = itPxTrack->pt();
      TrackPar.PxTrackPhi       = itPxTrack->phi();
      TrackPar.PxTrackEta       = itPxTrack->eta();

      FillTrackHistos("pixel", &TrackPar);

      RecoTrk = edm::RefToBase<Track>(PxTracks, itPxTrack - PxTracks->begin());

      // ##########################################
      // # Search for the corresponding sim-track #
      // ##########################################
      if (PxTkReco2Sim.find(RecoTrk) != PxTkReco2Sim.end())
	{
	  TrkParticle = PxTkReco2Sim[RecoTrk];

	  if (TrkParticle.size() != 0)
	    {
	      SimTrk = TrkParticle.begin()->first;

	      // #######################################################################
	      // # If there are more pixel tracks matched with the same tracker track, #
	      // # then select only the first one                                      #
	      // #######################################################################
	      if (abs(SimTrk->charge()) >= 1)
		{
		  if (PxTrkKeys_Counter.find(SimTrk.key()) != PxTrkKeys_Counter.end())
		    {
		      PxTrkKeys_Counter[SimTrk.key()]++;
		      if ((itPxTrack->numberOfValidHits() >= MinValidHitsPx) && (itPxTrack->normalizedChi2() <= MaxChi2PxTrk))
			{
			  MomentumTP = parametersDefinerTP->momentum(iEvent, iSetup, *(SimTrk.get()));
			  Simpt      = sqrt(MomentumTP.perp2());
			  Simeta     = MomentumTP.eta();
			  if (Simpt >= MinPtTk) H_pix_trk_doubleCounter_eta->Fill(itPxTrack->eta());
			  if (fabs(Simeta) <= MaxEtaTkTrk) H_pix_trk_doubleCounter_pt->Fill(itPxTrack->pt());
			}
		    }
		  else
		    {
		      PxTrkKeys_Counter[SimTrk.key()] = 1;

		      MomentumTP = parametersDefinerTP->momentum(iEvent, iSetup, *(SimTrk.get()));
		      VertexTP   = parametersDefinerTP->vertex(iEvent, iSetup, *(SimTrk.get())); // The SimVertex is refered to the beamspot
			  
		      Simpt    = sqrt(MomentumTP.perp2());
		      Simdz    = VertexTP.z() - (VertexTP.x()*MomentumTP.x() + VertexTP.y()*MomentumTP.y()) / Simpt * MomentumTP.z() / Simpt;
		      Simdxy   = -VertexTP.x()*sin(MomentumTP.phi()) + VertexTP.y()*cos(MomentumTP.phi());
		      Simtheta = MomentumTP.theta();
		      Simphi   = MomentumTP.phi();
		      Simeta   = MomentumTP.eta();

		      TrackPar.PxTrackValidHits = itPxTrack->numberOfValidHits();
		      TrackPar.TkTrackValidHits = SimTrk->matchedHit();
		      TrackPar.PxTrackLostHits  = itPxTrack->numberOfLostHits();
		      TrackPar.TkTrackLostHits  = 0;
		      TrackPar.PxTrackChi2DoF   = itPxTrack->normalizedChi2();
		      TrackPar.TkTrackChi2DoF   = 0.;
		      TrackPar.PxTrackPt        = itPxTrack->pt();
		      TrackPar.PxTrackPtErr     = itPxTrack->ptError();
		      TrackPar.TkTrackPt        = Simpt;
		      TrackPar.TkTrackPtErr     = 0.;

		      TrackPar.TkOriginZ        = Simdz;
		      TrackPar.TkOriginXY       = Simdxy;

		      TrackPar.BelongVx         = false;

		      if (SimTrk->status() == 1) // Track generated by Geant
			{
			  SimVx = SimTrk->parentVertex();

			  if ((SimVx->nSourceTracks() == 0) && (2.*(double)CountSimChargedTracks(SimVx)-3. >= MinTrkVxDoF))
			    {
			      VxSimTk = ParticleBase::Point(SimVx->position().x(),SimVx->position().y(),SimVx->position().z());
				  
			      TrackPar.PxTrackDz        = itPxTrack->dz(VxSimTk);
			      TrackPar.PxTrackDzErr     = itPxTrack->dzError();
			      TrackPar.TkTrackDz        = 0.;
			      TrackPar.TkTrackDzErr     = 0.;
				  
			      TrackPar.PxTrackD0        = -itPxTrack->dxy(VxSimTk);
			      TrackPar.PxTrackD0Err     = itPxTrack->dxyError();
			      TrackPar.TkTrackD0        = 0.;
			      TrackPar.TkTrackD0Err     = 0.;
				  
			      TrackPar.BelongVx = true;
			    }
			}
			  
		      TrackPar.PxTrackTheta     = itPxTrack->theta();
		      TrackPar.PxTrackThetaErr  = itPxTrack->thetaError()*itPxTrack->thetaError()*1000000.;
		      TrackPar.TkTrackTheta     = Simtheta;
		      TrackPar.TkTrackThetaErr  = 0.;
		      TrackPar.PxTrackPhi       = itPxTrack->phi();
		      TrackPar.PxTrackPhiErr    = itPxTrack->phiError()*itPxTrack->phiError()*1000000.;
		      TrackPar.TkTrackPhi       = Simphi;
		      TrackPar.TkTrackPhiErr    = 0.;
		      TrackPar.PxTrackEta       = itPxTrack->eta();
		      TrackPar.TkTrackEta       = Simeta;
	    
		      TrackPar.PxTkPtr          = itPxTrack;

		      FillTrackHistos("match", &TrackPar);
	    
		      if (PrintMsg == true) itPxTrack->hitPattern().print();
		    }
		} // End if charge >= 1
	    }
	} // End if SIM track exists
    } // End Loop Pixel Tracks
      
  for (map<unsigned int, unsigned int>::iterator it = PxTrkKeys_Counter.begin(); it != PxTrkKeys_Counter.end(); it++) H_pix_trk_doubleCounter->Fill(it->second);
  PxTrkKeys_Counter.clear();
  
  // ############################
  // # Loop over the sim-tracks #
  // ############################
  for (TrackingParticleCollection::size_type i = 0; i < TPC.size(); i++) {

    SimTrk = TrackingParticleRef(TPCollection, i);

    if (abs(SimTrk->charge()) >= 1)
      {
	MomentumTP = parametersDefinerTP->momentum(iEvent, iSetup, *(SimTrk.get()));
	VertexTP   = parametersDefinerTP->vertex(iEvent, iSetup, *(SimTrk.get()));    

	Simpt  = sqrt(MomentumTP.perp2());
	Simdz  = VertexTP.z() - (VertexTP.x()*MomentumTP.x() + VertexTP.y()*MomentumTP.y()) / Simpt * MomentumTP.z() / Simpt;
	Simdxy = -VertexTP.x()*sin(MomentumTP.phi()) + VertexTP.y()*cos(MomentumTP.phi());
	Simphi = MomentumTP.phi();
	Simeta = MomentumTP.eta();

	TrackPar.TkTrackValidHits = SimTrk->matchedHit();
	TrackPar.TkTrackChi2DoF   = 0.;
	TrackPar.TkTrackPt        = Simpt;
	TrackPar.TkOriginZ        = Simdz;
	TrackPar.TkOriginXY       = Simdxy;
	TrackPar.TkTrackPhi       = Simphi;
	TrackPar.TkTrackEta       = Simeta;
	TrackPar.AlgoEff          = false;

	if (FindNpixelHits(NpixelLayers, SimTrk) == true) TrackPar.AlgoEff = true;
	FillTrackHistos("strip", &TrackPar);
	
	if (((unsigned int)SimTrk->matchedHit() >= MinValidHitsTk) && (fabs(Simdxy) <= TrBound) && (fabs(Simdz) <= TzBound))
	  {
	    if (fabs(Simeta) <= MaxEtaTkTrk) H_trk_toteff_pt_den->Fill(Simpt);
	    if (Simpt >= MinPtTk)
	      {
		H_trk_toteff_eta_den->Fill(Simeta);
		
		if ((Simeta >= -EtaThrBin[2]) && (Simeta < -EtaThrBin[1]))
		  H_trk_toteff_phi_den[0]->Fill(Simphi/PI*180.);
		else if ((Simeta >= -EtaThrBin[1]) && (Simeta <= EtaThrBin[1]))
		  H_trk_toteff_phi_den[1]->Fill(Simphi/PI*180.);
		else if ((Simeta > EtaThrBin[1]) && (Simeta <= EtaThrBin[2]))
		  H_trk_toteff_phi_den[2]->Fill(Simphi/PI*180.);
	      }
	    
	    if (FindNpixelHits(NpixelLayers, SimTrk) == true)
	      {
		if (fabs(Simeta) <= MaxEtaTkTrk) H_trk_eff_pt_den->Fill(Simpt);
		if (Simpt >= MinPtTk)
		  {
		    H_trk_eff_eta_den->Fill(Simeta);
		    
		    if ((Simeta >= -EtaThrBin[2]) && (Simeta < -EtaThrBin[1]))
		      H_trk_eff_phi_den[0]->Fill(Simphi/PI*180.);
		    else if ((Simeta >= -EtaThrBin[1]) && (Simeta <= EtaThrBin[1]))
		      H_trk_eff_phi_den[1]->Fill(Simphi/PI*180.);
		    else if ((Simeta > EtaThrBin[1]) && (Simeta <= EtaThrBin[2]))
		      H_trk_eff_phi_den[2]->Fill(Simphi/PI*180.);
		  }
	      }
	  }

	if (TkTkSim2Reco.find(SimTrk) != TkTkSim2Reco.end())
	  {
	    TrkReco = TkTkSim2Reco[SimTrk];
	    
	    if (TrkReco.size() != 0)
	      {
		RecoTrk = TrkReco.begin()->first;

		if (((unsigned int)SimTrk->matchedHit() >= MinValidHitsTk) && (fabs(Simdxy) <= TrBound) && (fabs(Simdz) <= TzBound))
		  {
		    if (fabs(Simeta) <= MaxEtaTkTrk) H_trk_toteff_pt_num->Fill(Simpt);
		    if (Simpt >= MinPtTk)
		      {
			H_trk_toteff_eta_num->Fill(Simeta);
			  
			if ((Simeta >= -EtaThrBin[2]) && (Simeta < -EtaThrBin[1]))
			  H_trk_toteff_phi_num[0]->Fill(Simphi/PI*180.);
			else if ((Simeta >= -EtaThrBin[1]) && (Simeta <= EtaThrBin[1]))
			  H_trk_toteff_phi_num[1]->Fill(Simphi/PI*180.);
			else if ((Simeta > EtaThrBin[1]) && (Simeta <= EtaThrBin[2]))
			  H_trk_toteff_phi_num[2]->Fill(Simphi/PI*180.);
		      }
		      
		    if ((FindNpixelHits(NpixelLayers, &(*RecoTrk)) == true) &&
			(FindNpixelHits(NpixelLayers, SimTrk) == true))
		      {
			if (fabs(Simeta) <= MaxEtaTkTrk) H_trk_eff_pt_num->Fill(Simpt);
			if (Simpt >= MinPtTk)
			  {
			    H_trk_eff_eta_num->Fill(Simeta);
			      
			    if ((Simeta >= -EtaThrBin[2]) && (Simeta < -EtaThrBin[1]))
			      H_trk_eff_phi_num[0]->Fill(Simphi/PI*180.);
			    else if ((Simeta >= -EtaThrBin[1]) && (Simeta <= EtaThrBin[1]))
			      H_trk_eff_phi_num[1]->Fill(Simphi/PI*180.);
			    else if ((Simeta > EtaThrBin[1]) && (Simeta <= EtaThrBin[2]))
			      H_trk_eff_phi_num[2]->Fill(Simphi/PI*180.);
			  }
		      }
		  }
	      }
	  }
      }
  } // End Loop Pixel SimTracks
}


void MyPixAnalyzer::PixelTracksData (const edm::Event& iEvent,
				     const edm::EventSetup& iSetup)
{
  edm::Handle<BeamSpot> beamSpotH;
  iEvent.getByLabel("offlineBeamSpot",beamSpotH);
  BeamSpot beamSpot = *beamSpotH;
 
  edm::Handle<VertexCollection> offVertexCollection;
  iEvent.getByLabel("offlinePrimaryVertices",offVertexCollection);
  
  edm::Handle<edm::View<Track> > TkTracks;
  iEvent.getByLabel("generalTracks",TkTracks);

  edm::Handle<edm::View<Track> > PxTracks;
  iEvent.getByLabel("pixelTracks",PxTracks);

  TrackParType TrackPar;

  map<unsigned int, unsigned int> PxTrkKeys_Counter;
  map<unsigned int, unsigned int> TkTrkKeys_Counter;

  stringstream Name;
  bool PixelHitMatched;
  vector<string> TrkNoMatchHitLocation;
  vector<string> PixMatchHitLocation;
  vector<string> PixNoMatchHitLocation;

  // ##############################
  // # Loop over the pixel tracks #
  // ##############################
  for (edm::View<Track>::const_iterator itPxTrack = PxTracks->begin(); itPxTrack != PxTracks->end(); itPxTrack++) {	 

    TrackPar.PxTrackValidHits = itPxTrack->numberOfValidHits();
    TrackPar.PxTrackChi2DoF   = itPxTrack->normalizedChi2();
    TrackPar.PxTrackPt        = itPxTrack->pt();
    TrackPar.PxTrackPhi       = itPxTrack->phi();
    TrackPar.PxTrackEta       = itPxTrack->eta();

    FillTrackHistos("pixel", &TrackPar);
	    
    const Track* PxTrack = &(*itPxTrack);
    unsigned int FirstPixelRecPxClusterRow;
    unsigned int FirstPixelRecPxClusterCol;
    unsigned int LastPixelRecPxClusterRow;
    unsigned int LastPixelRecPxClusterCol;
    unsigned int IdPixelRecPxCluster;
    unsigned int PixTrkHitMatched;

    int TkTkIndx = -1;
    double BestTrkTrkChi2 = RangeChi2;

    // ################################
    // # Loop over the tracker tracks #
    // ################################
    for (edm::View<Track>::const_iterator itTkTrack = TkTracks->begin(); itTkTrack != TkTracks->end(); itTkTrack++) {
	      
      PixTrkHitMatched = 0;

      // #############################################
      // # Loop over the pixel hits of a pixel track #
      // #############################################
      for (trackingRecHit_iterator itPxHit = PxTrack->recHitsBegin(); itPxHit != PxTrack->recHitsEnd(); itPxHit++) {	

	if ((*itPxHit)->isValid() == true)
	  {
	    // Extract cluster information for the pixel hit
	    const SiPixelRecHit* recHitPxPixel = dynamic_cast<const SiPixelRecHit*>((*itPxHit)->clone());
	    SiPixelRecHit::ClusterRef const& PixelRecPxCluster = recHitPxPixel->cluster();
	    FirstPixelRecPxClusterRow = PixelRecPxCluster->minPixelRow();
	    FirstPixelRecPxClusterCol = PixelRecPxCluster->minPixelCol();
	    LastPixelRecPxClusterRow = PixelRecPxCluster->maxPixelRow();
	    LastPixelRecPxClusterCol = PixelRecPxCluster->maxPixelCol();
	    IdPixelRecPxCluster = recHitPxPixel->geographicalId()();
		    
	    const Track* TkTrack = &(*itTkTrack);
	    unsigned int FirstPixelRecTkClusterRow;
	    unsigned int FirstPixelRecTkClusterCol;
	    unsigned int LastPixelRecTkClusterRow;
	    unsigned int LastPixelRecTkClusterCol;
	    unsigned int IdPixelRecTkCluster;

	    PixelHitMatched = false;

	    // ###############################################
	    // # Loop over the pixel hits of a tracker track #
	    // ###############################################
	    for (trackingRecHit_iterator itTkHit = TkTrack->recHitsBegin(); itTkHit != TkTrack->recHitsEnd(); itTkHit++) {

	      if (((*itTkHit)->isValid() == true) &&
		  (((*itTkHit)->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) ||
		   ((*itTkHit)->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)))
		{
		  // Extract cluster information for the tracker hit
		  const SiPixelRecHit* recHitTkPixel = dynamic_cast<const SiPixelRecHit*>((*itTkHit)->clone());
		  SiPixelRecHit::ClusterRef const& PixelRecCluster = recHitTkPixel->cluster();
		  FirstPixelRecTkClusterRow = PixelRecCluster->minPixelRow();
		  FirstPixelRecTkClusterCol = PixelRecCluster->minPixelCol();
		  LastPixelRecTkClusterRow = PixelRecCluster->maxPixelRow();
		  LastPixelRecTkClusterCol = PixelRecCluster->maxPixelCol();
		  IdPixelRecTkCluster = recHitTkPixel->geographicalId()();
			  
		  if ((IdPixelRecTkCluster == IdPixelRecPxCluster) &&
		      (FirstPixelRecTkClusterRow == FirstPixelRecPxClusterRow) &&
		      (FirstPixelRecTkClusterCol == FirstPixelRecPxClusterCol) &&
		      (LastPixelRecTkClusterRow == LastPixelRecPxClusterRow) &&
		      (LastPixelRecTkClusterCol == LastPixelRecPxClusterCol))
		    {
		      if (PrintMsg == true)
			{
			  cout << "### Matched hit ###" << endl;
			  cout << "\tFirstPixelRecClusterRow: " << FirstPixelRecPxClusterRow << endl;
			  cout << "\tFirstPixelRecClusterCol: " << FirstPixelRecPxClusterCol << endl;
			  cout << "\tLastPixelRecClusterRow: " << LastPixelRecPxClusterRow << endl;
			  cout << "\tLastPixelRecClusterCol: " << LastPixelRecPxClusterCol << endl;
			  cout << "\tIdPixelRecCluster: " << IdPixelRecPxCluster << endl;
			}
			      
		      PixelHitMatched = true;
		      PixTrkHitMatched++;
		      delete recHitTkPixel;
		      break;
		    }
		  else
		    {
		      Name.str("");
		      if ((*itTkHit)->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel)
			{
			  Name << "L" << ((PXBDetId)((*itTkHit)->geographicalId())).layer();
			  TrkNoMatchHitLocation.push_back(Name.str());
			}
		      else if ((*itTkHit)->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)
			{
			  Name << "D" << ((PXFDetId)((*itTkHit)->geographicalId())).disk();
			  Name << "S" << ((PXFDetId)((*itTkHit)->geographicalId())).side();
			  TrkNoMatchHitLocation.push_back(Name.str());
			}
		    }
		      
		  delete recHitTkPixel;
		}
	    } // End Loop on Pixel-hits of a Tracker Track

	    if (PixelHitMatched == true)
	      {				
		Name.str("");
		if ((*itPxHit)->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel)
		  {
		    Name << "L" << ((PXBDetId)((*itPxHit)->geographicalId())).layer();
		    PixMatchHitLocation.push_back(Name.str());
		  }
		else if ((*itPxHit)->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)
		  {
		    Name << "D" << ((PXFDetId)((*itPxHit)->geographicalId())).disk();
		    Name << "S" << ((PXFDetId)((*itPxHit)->geographicalId())).side();
		    PixMatchHitLocation.push_back(Name.str());
		  }
	      }
	    else
	      {
		Name.str("");
		if ((*itPxHit)->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel)
		  {
		    Name << "L" << ((PXBDetId)((*itPxHit)->geographicalId())).layer();
		    PixNoMatchHitLocation.push_back(Name.str());
		  }
		else if ((*itPxHit)->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)
		  {
		    Name << "D" << ((PXFDetId)((*itPxHit)->geographicalId())).disk();
		    Name << "S" << ((PXFDetId)((*itPxHit)->geographicalId())).side();
		    PixNoMatchHitLocation.push_back(Name.str());
		  }				    
	      }
		
	    delete recHitPxPixel;
	  }
      } // End Loop on Pixel-hits of a Pixel Track
	      
      if (PixTrkHitMatched >= MinHitsMatch)
	{
	  // ############################################################################
	  // # If there are hits of the tracker track that are not matched and if they  #
	  // # do not belong to the same Layer or Disk of the pixel matched hits and if #
	  // # there are not matched pixel hits, then do not perform the match          #
	  // ############################################################################
	  unsigned int i;
	  for (i = 0; i < TrkNoMatchHitLocation.size(); i++)
	    if ((find(PixMatchHitLocation.begin(), PixMatchHitLocation.end(), TrkNoMatchHitLocation[i]) == PixMatchHitLocation.end()) && (PixNoMatchHitLocation.size() != 0)) break;
	  if (i == TrkNoMatchHitLocation.size())
	    {
	      // #######################################################################################################################
	      // # - If there is still NONE tracker tracks matched with the pixel track OR                                             #
	      // # - If the previous matched tracker track is one of the "already matched" tracker tracks AND                          #
	      // # the current matched tracker track is NOT one of the "already matched" tracker tracks OR                             #
	      // # - If the current matched tracker track is NOT one of the "already matched" tracker tracks AND it has a better Chi2, #
	      // # Than save the match                                                                                                 #
	      // #######################################################################################################################
	      edm::RefToBase<Track> RefTkTrack(TkTracks, itTkTrack - TkTracks->begin());
	      if (TkTkIndx != -1)
		{
		  edm::RefToBase<Track> RefOldTkTrack(TkTracks, TkTkIndx);
		      
		  if ((BestTrkTrkChi2 == RangeChi2) ||
		      ((PxTrkKeys_Counter.find(RefOldTkTrack.key()) != PxTrkKeys_Counter.end()) && (PxTrkKeys_Counter.find(RefTkTrack.key()) == PxTrkKeys_Counter.end())) ||
		      ((PxTrkKeys_Counter.find(RefTkTrack.key()) == PxTrkKeys_Counter.end()) && (itTkTrack->normalizedChi2() < BestTrkTrkChi2)))
		    {
		      BestTrkTrkChi2 = itTkTrack->normalizedChi2();
		      TkTkIndx = itTkTrack - TkTracks->begin();
		    }
		}
	      else if ((BestTrkTrkChi2 == RangeChi2) ||
		       ((PxTrkKeys_Counter.find(RefTkTrack.key()) == PxTrkKeys_Counter.end()) && (itTkTrack->normalizedChi2() < BestTrkTrkChi2)))
		{
		  BestTrkTrkChi2 = itTkTrack->normalizedChi2();
		  TkTkIndx = itTkTrack - TkTracks->begin();
		}
		  
	      edm::RefToBase<Track> RefPxTrack(PxTracks, itPxTrack - PxTracks->begin());
	      // #######################################################################
	      // # Count the multiple tracker tracks matched with the same pixel track #
	      // #######################################################################
	      if (TkTrkKeys_Counter.find(RefPxTrack.key()) != TkTrkKeys_Counter.end())
		{
		  TkTrkKeys_Counter[RefPxTrack.key()]++;
		  if ((itTkTrack->numberOfValidHits() >= MinValidHitsTk) && (itTkTrack->normalizedChi2() <= MaxChi2TkTrk))
		    {
		      if (itTkTrack->pt() >= MinPtTk) H_trk_trk_doubleCounter_eta->Fill(itTkTrack->eta());
		      if (fabs(itTkTrack->eta()) <= MaxEtaTkTrk) H_trk_trk_doubleCounter_pt->Fill(itTkTrack->pt());
		    }
		}
	      else TkTrkKeys_Counter[RefPxTrack.key()] = 1; 
	    }
	}
	  
      TrkNoMatchHitLocation.clear();
      PixMatchHitLocation.clear();
      PixNoMatchHitLocation.clear();
    } // End Loop Tracker Tracks

    for (map<unsigned int, unsigned int>::iterator it = TkTrkKeys_Counter.begin(); it != TkTrkKeys_Counter.end(); it++) H_trk_trk_doubleCounter->Fill(it->second);
    TkTrkKeys_Counter.clear();
	
    if (TkTkIndx >= 0)
      {
	edm::View<Track>::const_iterator itTkTrack = TkTracks->begin() + TkTkIndx;
	edm::RefToBase<Track> RefTkTrack(TkTracks, TkTkIndx);
	    
	// #######################################################################
	// # If there are more pixel tracks matched with the same tracker track, #
	// # then select only the first one                                      #
	// #######################################################################
	if (PxTrkKeys_Counter.find(RefTkTrack.key()) != PxTrkKeys_Counter.end())
	  {
	    PxTrkKeys_Counter[RefTkTrack.key()]++;
	    if ((itPxTrack->numberOfValidHits() >= MinValidHitsPx) && (itPxTrack->normalizedChi2() <= MaxChi2PxTrk))
	      {
		if (itTkTrack->pt() >= MinPtTk) H_pix_trk_doubleCounter_eta->Fill(itPxTrack->eta());
		if (fabs(itTkTrack->eta()) <= MaxEtaTkTrk) H_pix_trk_doubleCounter_pt->Fill(itPxTrack->pt());
	      }
	  }
	else
	  {
	    PxTrkKeys_Counter[RefTkTrack.key()] = 1;

	    TrackPar.PxTrackValidHits = itPxTrack->numberOfValidHits();
	    TrackPar.TkTrackValidHits = itTkTrack->numberOfValidHits();
	    TrackPar.PxTrackLostHits  = itPxTrack->numberOfLostHits();
	    TrackPar.TkTrackLostHits  = itTkTrack->numberOfLostHits();
	    TrackPar.PxTrackChi2DoF   = itPxTrack->normalizedChi2();
	    TrackPar.TkTrackChi2DoF   = itTkTrack->normalizedChi2();
	    TrackPar.PxTrackPt        = itPxTrack->pt();
	    TrackPar.PxTrackPtErr     = itPxTrack->ptError();
	    TrackPar.TkTrackPt        = itTkTrack->pt();
	    TrackPar.TkTrackPtErr     = itTkTrack->ptError();

	    TrackPar.TkOriginZ        = itTkTrack->dz(beamSpot.position());
	    TrackPar.TkOriginXY       = itTkTrack->dxy(beamSpot.position());

	    TrackPar.BelongVx         = false;
	    TrackPar.PxTrackDz        = TzBound;

	    vector<Vertex>::const_iterator itTkVx;
	    vector<TrackBaseRef>::const_iterator itVxTkTrack;
	    for (itTkVx = offVertexCollection->begin(); itTkVx != offVertexCollection->end(); itTkVx++)
	      if ((itTkVx->isValid() == true) && (itTkVx->isFake() == false) &&
		  (itTkVx->ndof() >= MinTrkVxDoF))
		{
		  for (itVxTkTrack = itTkVx->tracks_begin(); itVxTkTrack != itTkVx->tracks_end(); itVxTkTrack++)
		    if ((itVxTkTrack->key() == RefTkTrack.key()) && (itTkVx->trackWeight(*itVxTkTrack) >= MinTrkVxWgt)) break;
		  if (itVxTkTrack != itTkVx->tracks_end()) break;
		}

	    if (itTkVx != offVertexCollection->end())
	      {
		// #######################################################
		// # Obtain the configuration of the vertexing algorithm #
		// #######################################################
		const edm::Provenance* prov = offVertexCollection.provenance();
		edm::ParameterSetID psid = prov->psetID();
		edm::pset::Registry* psregistry = edm::pset::Registry::instance();
		edm::ParameterSet psetFromProvenance;
		psregistry->getMapped(psid, psetFromProvenance);
		    
		// ###############################################################
		// # Find the track that needs to be removed and save the others #
		// ###############################################################
		vector<TrackBaseRef> pixelLess;
		for (itVxTkTrack = itTkVx->tracks_begin(); itVxTkTrack != itTkVx->tracks_end(); itVxTkTrack++)
		  if (itVxTkTrack->key() != RefTkTrack.key()) pixelLess.push_back(*itVxTkTrack);
		    
		edm::ESHandle<TransientTrackBuilder> TTBuilder;
		iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

		vector<TransientTrack> TTrks;
		TTrks.reserve(pixelLess.size());
		for (vector<TrackBaseRef>::const_iterator it = pixelLess.begin(); it != pixelLess.end(); it++)
		  {
		    TTrks.push_back((*TTBuilder).build(*(edm::RefToBase<Track>(*it))));
		    TTrks.back().setBeamSpot(*beamSpotH);
		  }
		    
		// #####################
		// # Re-fit the vertex #
		// #####################
		auto_ptr<PrimaryVertexProducerAlgorithm> VxProdAlgo;
		VxProdAlgo.reset(new PrimaryVertexProducerAlgorithm(psetFromProvenance));
		vector<TransientVertex> MyPrimVx = VxProdAlgo->vertices(TTrks, *beamSpotH);

		if ((MyPrimVx.size() == 1) && (Vertex(MyPrimVx.front()).isValid() == true) && (Vertex(MyPrimVx.front()).isFake() == false) &&
		    (Vertex(MyPrimVx.front()).ndof() >= MinTrkVxDoF))
		  {
		    TrackPar.BelongVx     = true;

		    TrackPar.PxTrackDz    = itPxTrack->dz(Vertex(MyPrimVx.front()).position());
		    TrackPar.PxTrackDzErr = itPxTrack->dzError();
		    TrackPar.TkTrackDz    = itTkTrack->dz(Vertex(MyPrimVx.front()).position());
		    TrackPar.TkTrackDzErr = Vertex(MyPrimVx.front()).zError() * TkZVxErrCorr;

		    TrackPar.PxTrackD0    = -itPxTrack->dxy(Vertex(MyPrimVx.front()).position());
		    TrackPar.PxTrackD0Err = itPxTrack->dxyError();
		    TrackPar.TkTrackD0    = -itTkTrack->dxy(Vertex(MyPrimVx.front()).position());
		    TrackPar.TkTrackD0Err = sqrt(1/(Vertex(MyPrimVx.front()).x()*Vertex(MyPrimVx.front()).x() + Vertex(MyPrimVx.front()).y()*Vertex(MyPrimVx.front()).y()) *
						 (Vertex(MyPrimVx.front()).x()*Vertex(MyPrimVx.front()).x() * Vertex(MyPrimVx.front()).xError()*Vertex(MyPrimVx.front()).xError()
						  * TkXYVxErrCorr*TkXYVxErrCorr +
						  Vertex(MyPrimVx.front()).y()*Vertex(MyPrimVx.front()).y() * Vertex(MyPrimVx.front()).yError()*Vertex(MyPrimVx.front()).yError()
						  * TkXYVxErrCorr*TkXYVxErrCorr));
		  }

		pixelLess.clear();
		TTrks.clear();
	      }

	    TrackPar.PxTrackTheta     = itPxTrack->theta();
	    TrackPar.PxTrackThetaErr  = itPxTrack->thetaError()*itPxTrack->thetaError()*1000000.;
	    TrackPar.TkTrackTheta     = itTkTrack->theta();
	    TrackPar.TkTrackThetaErr  = itTkTrack->thetaError()*itTkTrack->thetaError()*1000000.;
	    TrackPar.PxTrackPhi       = itPxTrack->phi();
	    TrackPar.PxTrackPhiErr    = itPxTrack->phiError()*itPxTrack->phiError()*1000000.;
	    TrackPar.TkTrackPhi       = itTkTrack->phi();
	    TrackPar.TkTrackPhiErr    = itTkTrack->phiError()*itTkTrack->phiError()*1000000.;
	    TrackPar.PxTrackEta       = itPxTrack->eta();
	    TrackPar.TkTrackEta       = itTkTrack->eta();
		  
	    TrackPar.PxTkPtr          = itPxTrack;

	    FillTrackHistos("match",&TrackPar);

	    if (PrintMsg == true) itPxTrack->hitPattern().print();
	  }
      }
  } // End Loop Pixel Tracks
      
  for (map<unsigned int, unsigned int>::iterator it = PxTrkKeys_Counter.begin(); it != PxTrkKeys_Counter.end(); it++) H_pix_trk_doubleCounter->Fill(it->second);
  PxTrkKeys_Counter.clear();
      
  for (edm::View<Track>::const_iterator itTkTrack = TkTracks->begin(); itTkTrack != TkTracks->end(); itTkTrack++) {
	
    TrackPar.TkTrackValidHits = itTkTrack->numberOfValidHits();
    TrackPar.TkTrackChi2DoF   = itTkTrack->normalizedChi2();
    TrackPar.TkTrackPt        = itTkTrack->pt();
    TrackPar.TkOriginZ        = itTkTrack->dz(beamSpot.position());
    TrackPar.TkOriginXY       = itTkTrack->dxy(beamSpot.position());
    TrackPar.TkTrackPhi       = itTkTrack->phi();
    TrackPar.TkTrackEta       = itTkTrack->eta();
    TrackPar.AlgoEff          = false;

    if (FindNpixelHits(NpixelLayers, &(*itTkTrack)) == true) TrackPar.AlgoEff = true;
	
    FillTrackHistos("strip", &TrackPar);
  } // End Loop Tracker Tracks
}


// ################################
// # @@@ METHODS for VERTICES @@@ #
// ################################

bool MyPixAnalyzer::FillVertexHistos (string DataType,
				      VertexParType* VertexPar,
				      const edm::Event& iEvent,
				      const edm::EventSetup& iSetup)
{
  bool result = false;

  if ((DataType.compare("pixel") == 0) &&
      (VertexPar->PxVxDoFOrig >= MinVxDoF))
    {
      H_pix_vx_z->Fill(VertexPar->PxVxZ);
      H_pix_vx_xy->Fill(VertexPar->PxVxX*10000., VertexPar->PxVxY*10000.);
      H_pix_vx_xerr->Fill(VertexPar->PxVxXerr);
      H_pix_vx_yerr->Fill(VertexPar->PxVxYerr);
      H_pix_vx_zerr->Fill(VertexPar->PxVxZerr);
      H_pix_vx_normchi2->Fill(VertexPar->PxVxChi2DoF);
      H_pix_vx_dof->Fill(VertexPar->PxVxDoFOrig);
      H_ntrk_pix_vx->Fill(VertexPar->PxVxNTracksOrig);
      H_vx_counters->SetBinContent(H_vx_counters->GetXaxis()->FindBin("Pixel"), H_vx_counters->GetBinContent(H_vx_counters->GetXaxis()->FindBin("Pixel")) + 1.);

      result = true;
    }
  else if (DataType.compare("match") == 0)
    {
      // ##################################
      // # Pixel-Vertex tracks properties #
      // ##################################
      for (vector<TrackBaseRef>::const_iterator itPxTrack = VertexPar->PxVx->tracks_begin(); itPxTrack != VertexPar->PxVx->tracks_end(); itPxTrack++)
	if ((*itPxTrack)->normalizedChi2() <= MaxChi2VxTrk)
	  {
	    H_pix_vx_trk_pt->Fill((*itPxTrack)->pt());
	    H_pix_vx_trk_normchi2->Fill((*itPxTrack)->normalizedChi2());
	    
	    if ((find(VertexPar->MatchedTkTracksPix->begin(), VertexPar->MatchedTkTracksPix->end(), itPxTrack->key()) == VertexPar->MatchedTkTracksPix->end()))
	      {
		H_pix_vx_trk_pt_nomatch->Fill((*itPxTrack)->pt());
		H_pix_vx_trk_eta_nomatch->Fill((*itPxTrack)->eta());
		H_pix_vx_trk_normchi2_nomatch->Fill((*itPxTrack)->normalizedChi2());
	      }
	  }

      // ############################################################
      // # Efficiency & Purity and Tracker-Vertex tracks properties #
      // ############################################################
      unsigned int goodTrkCounter = 0;
      if (VertexPar->IsTrkPart == false)
	{
	  for (vector<TrackBaseRef>::const_iterator itTkTrack = VertexPar->TkVxReco->tracks_begin(); itTkTrack != VertexPar->TkVxReco->tracks_end(); itTkTrack++)
	    if (((*itTkTrack)->pt() >= MinPtVx) && (fabs((*itTkTrack)->eta()) <= MaxEtaVxTrk) && ((*itTkTrack)->normalizedChi2() <= MaxChi2VxTrk))
	      {
		goodTrkCounter++;
		
		H_trk_vx_trk_pt->Fill((*itTkTrack)->pt());
		H_trk_vx_trk_normchi2->Fill((*itTkTrack)->normalizedChi2());
		
		if ((find(VertexPar->MatchedTkTracksTrk->begin(), VertexPar->MatchedTkTracksTrk->end(), itTkTrack->key()) == VertexPar->MatchedTkTracksTrk->end()))
		  {
		    H_trk_vx_trk_pt_nomatch->Fill((*itTkTrack)->pt());
		    H_trk_vx_trk_eta_nomatch->Fill((*itTkTrack)->eta());
		    H_trk_vx_trk_normchi2_nomatch->Fill((*itTkTrack)->normalizedChi2());		 
		  }
	      }
	}
      else
	{
	  for (TrackingParticleRefVector::iterator itTkTrack = VertexPar->TkVxTP->daughterTracks_begin(); itTkTrack != VertexPar->TkVxTP->daughterTracks_end(); itTkTrack++)
	    if ((abs((*itTkTrack)->charge()) >= 1) && ((*itTkTrack)->pt() >= MinPtVx) && (fabs((*itTkTrack)->eta()) <= MaxEtaVxTrk))
	      {
		goodTrkCounter++;
		
		H_trk_vx_trk_pt->Fill((*itTkTrack)->pt());
		
		if ((find(VertexPar->MatchedTkTracksTrk->begin(), VertexPar->MatchedTkTracksTrk->end(), itTkTrack->key()) == VertexPar->MatchedTkTracksTrk->end()))
		  {
		    H_trk_vx_trk_pt_nomatch->Fill((*itTkTrack)->pt());
		    H_trk_vx_trk_eta_nomatch->Fill((*itTkTrack)->eta());
		  }
	      }
	}

      // ##############
      // # Efficiency #
      // ##############
      if ((goodTrkCounter >= MinVxTrkMatch) &&
	  (VertexPar->TkVxDoFOrig >= MinVxDoF) &&
	  (sqrt(VertexPar->TkVxXOrig*VertexPar->TkVxXOrig + VertexPar->TkVxYOrig*VertexPar->TkVxYOrig) <= VrBound) && (fabs(VertexPar->TkVxZOrig) <= VzBound))
	H_vx_counters->SetBinContent(H_vx_counters->GetXaxis()->FindBin("MatchEff"), H_vx_counters->GetBinContent(H_vx_counters->GetXaxis()->FindBin("MatchEff")) + 1.);

      // ##########
      // # Purity #
      // ##########
      if (VertexPar->PxVxDoFOrig >= MinVxDoF)
	H_vx_counters->SetBinContent(H_vx_counters->GetXaxis()->FindBin("MatchPur"), H_vx_counters->GetBinContent(H_vx_counters->GetXaxis()->FindBin("MatchPur")) + 1.);

      result = true;

      if ((VertexPar->IsTrkPart == true) || ((VertexPar->PxVxReFit != NULL) && (VertexPar->TkVxReFit != NULL)))
	{
	  // @@@@@@@@@@ MAKE CRASH @@@@@@@@@@
	  // ####################################################
	  // # Properties of the tracks of the re-fitted vertex #
	  // ####################################################
	  //       for (vector<TransientTrack>::const_iterator itTkTrack = VertexPar->TkVxReFit->originalTracks().begin(); itTkTrack != VertexPar->TkVxReFit->originalTracks().end(); itTkTrack++)
	  // 	{
	  // 	  H_trk_vxrefit_trk_pt_nomatch->Fill((*itTkTrack->trackBaseRef()).pt());
	  // 	  H_trk_vxrefit_trk_eta_nomatch->Fill((*itTkTrack->trackBaseRef()).eta());
	  // 	}
	  // @@@@@@@@@@ MAKE CRASH @@@@@@@@@@

	  // ##########################
	  // # Resolutions and Errors #
	  // ##########################
	  stringstream String;
	  String.precision(3);
	  String << rint(VertexPar->PxVxNTracks);

 	  if (map_H_vx_x_diff.find(String.str().c_str()) == map_H_vx_x_diff.end())
	    {
	      stringstream HistoName;

	      HistoName << "x_diff_Ntrk_" << String.str();
	      map_H_vx_x_diff[String.str().c_str()] = SubDirPixVxNtrks->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 300, -xRange/2., xRange/2.);
	      HistoName.str("");
	      HistoName << "y_diff_Ntrk_" << String.str();
	      map_H_vx_y_diff[String.str().c_str()] = SubDirPixVxNtrks->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 300, -yRange/2., yRange/2.);
	      HistoName.str("");
	      HistoName << "z_diff_Ntrk_" << String.str();
	      map_H_vx_z_diff[String.str().c_str()] = SubDirPixVxNtrks->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 300, -zRange/2., zRange/2.);

	      HistoName.str("");
	      HistoName << "x_diff_pull_Ntrk_" << String.str();
	      map_H_vx_x_diff_pull[String.str().c_str()] = SubDirPixVxNtrks->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 100, -pullsRange/2., pullsRange/2.);
	      HistoName.str("");
	      HistoName << "y_diff_pull_Ntrk_" << String.str();
	      map_H_vx_y_diff_pull[String.str().c_str()] = SubDirPixVxNtrks->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 100, -pullsRange/2., pullsRange/2.);
	      HistoName.str("");
	      HistoName << "z_diff_pull_Ntrk_" << String.str();
	      map_H_vx_z_diff_pull[String.str().c_str()] = SubDirPixVxNtrks->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 100, -pullsRange/2., pullsRange/2.);

	      HistoName.str("");
	      HistoName << "x_err2_Ntrk_" << String.str();
	      map_H_vx_x_err2[String.str().c_str()] = SubDirPixVxNtrks->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 1000, 0., xerrRange);
	      HistoName.str("");
	      HistoName << "y_err2_Ntrk_" << String.str();
	      map_H_vx_y_err2[String.str().c_str()] = SubDirPixVxNtrks->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 1000, 0., yerrRange);
	      HistoName.str("");
	      HistoName << "z_err2_Ntrk_" << String.str();
	      map_H_vx_z_err2[String.str().c_str()] = SubDirPixVxNtrks->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 1000, 0., zerrRange);
	    }

	  map_H_vx_x_diff[String.str().c_str()]->Fill((VertexPar->TkVxX-VertexPar->PxVxX)*10000.);
	  map_H_vx_y_diff[String.str().c_str()]->Fill((VertexPar->TkVxY-VertexPar->PxVxY)*10000.);
	  map_H_vx_z_diff[String.str().c_str()]->Fill((VertexPar->TkVxZ-VertexPar->PxVxZ)*10000.);

	  map_H_vx_x_diff_pull[String.str().c_str()]->Fill((VertexPar->TkVxX-VertexPar->PxVxX) / sqrt(VertexPar->TkVxXerr*VertexPar->TkVxXerr + VertexPar->PxVxXerr*VertexPar->PxVxXerr));
	  map_H_vx_y_diff_pull[String.str().c_str()]->Fill((VertexPar->TkVxY-VertexPar->PxVxY) / sqrt(VertexPar->TkVxYerr*VertexPar->TkVxYerr + VertexPar->PxVxYerr*VertexPar->PxVxYerr));
	  map_H_vx_z_diff_pull[String.str().c_str()]->Fill((VertexPar->TkVxZ-VertexPar->PxVxZ) / sqrt(VertexPar->TkVxZerr*VertexPar->TkVxZerr + VertexPar->PxVxZerr*VertexPar->PxVxZerr));

	  map_H_vx_x_err2[String.str().c_str()]->Fill(VertexPar->TkVxXerr*VertexPar->TkVxXerr*100000000.);
	  map_H_vx_y_err2[String.str().c_str()]->Fill(VertexPar->TkVxYerr*VertexPar->TkVxYerr*100000000.);
	  map_H_vx_z_err2[String.str().c_str()]->Fill(VertexPar->TkVxZerr*VertexPar->TkVxZerr*100000000.);

	  if (VertexPar->IsTrkPart == false)
	    {
	      String.str("");
	      String << rint(VertexPar->TkVxNTracks);

	      if (map_H_vx_x_diff_orig.find(String.str().c_str()) == map_H_vx_x_diff_orig.end())
		{
		  stringstream HistoName;

		  HistoName.str("");
		  HistoName << "x_diff_orig_Ntrk_" << String.str();
		  map_H_vx_x_diff_orig[String.str().c_str()] = SubDirPixVxNtrks->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 300, -xRange/2., xRange/2.);
		  HistoName.str("");
		  HistoName << "y_diff_orig_Ntrk_" << String.str();
		  map_H_vx_y_diff_orig[String.str().c_str()] = SubDirPixVxNtrks->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 300, -yRange/2., yRange/2.);
		  HistoName.str("");
		  HistoName << "z_diff_orig_Ntrk_" << String.str();
		  map_H_vx_z_diff_orig[String.str().c_str()] = SubDirPixVxNtrks->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 300, -zRange/2., zRange/2.);

		  HistoName.str("");
		  HistoName << "x_err2_non_orig_Ntrk_" << String.str();
		  map_H_vx_x_err2_non_orig[String.str().c_str()] = SubDirPixVxNtrks->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 1000, 0., xerrRange);
		  HistoName.str("");
		  HistoName << "y_err2_non_orig_Ntrk_" << String.str();
		  map_H_vx_y_err2_non_orig[String.str().c_str()] = SubDirPixVxNtrks->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 1000, 0., yerrRange);
		  HistoName.str("");
		  HistoName << "z_err2_non_orig_Ntrk_" << String.str();
		  map_H_vx_z_err2_non_orig[String.str().c_str()] = SubDirPixVxNtrks->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 1000, 0., zerrRange);
		}

	      map_H_vx_x_diff_orig[String.str().c_str()]->Fill((VertexPar->TkVxX-VertexPar->TkVxXOrig)*10000.);
	      map_H_vx_y_diff_orig[String.str().c_str()]->Fill((VertexPar->TkVxY-VertexPar->TkVxYOrig)*10000.);
	      map_H_vx_z_diff_orig[String.str().c_str()]->Fill((VertexPar->TkVxZ-VertexPar->TkVxZOrig)*10000.);

	      map_H_vx_x_err2_non_orig[String.str().c_str()]->Fill(VertexPar->TkVxXerr*VertexPar->TkVxXerr*100000000.);
	      map_H_vx_y_err2_non_orig[String.str().c_str()]->Fill(VertexPar->TkVxYerr*VertexPar->TkVxYerr*100000000.);
	      map_H_vx_z_err2_non_orig[String.str().c_str()]->Fill(VertexPar->TkVxZerr*VertexPar->TkVxZerr*100000000.);

	      String.str("");
	      String << rint(VertexPar->TkVxNTracksOrig);

	      if (map_H_vx_x_err2_orig.find(String.str().c_str()) == map_H_vx_x_err2_orig.end())
		{
		  stringstream HistoName;

		  HistoName.str("");
		  HistoName << "x_err2_orig_Ntrk_" << String.str();
		  map_H_vx_x_err2_orig[String.str().c_str()] = SubDirPixVxNtrks->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 1000, 0., xerrRange);
		  HistoName.str("");
		  HistoName << "y_err2_orig_Ntrk_" << String.str();
		  map_H_vx_y_err2_orig[String.str().c_str()] = SubDirPixVxNtrks->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 1000, 0., yerrRange);
		  HistoName.str("");
		  HistoName << "z_err2_orig_Ntrk_" << String.str();
		  map_H_vx_z_err2_orig[String.str().c_str()] = SubDirPixVxNtrks->make<TH1F>(HistoName.str().c_str(), HistoName.str().c_str(), 1000, 0., zerrRange);
		}

	      map_H_vx_x_err2_orig[String.str().c_str()]->Fill(VertexPar->TkVxXerrOrig*VertexPar->TkVxXerrOrig*100000000.);
	      map_H_vx_y_err2_orig[String.str().c_str()]->Fill(VertexPar->TkVxYerrOrig*VertexPar->TkVxYerrOrig*100000000.);
	      map_H_vx_z_err2_orig[String.str().c_str()]->Fill(VertexPar->TkVxZerrOrig*VertexPar->TkVxZerrOrig*100000000.);
	    }
	}
    }
  else if ((DataType.compare("strip") == 0) &&
	   (VertexPar->TkVxDoFOrig >= MinVxDoF) &&
	   (sqrt(VertexPar->TkVxXOrig*VertexPar->TkVxXOrig + VertexPar->TkVxYOrig*VertexPar->TkVxYOrig) <= VrBound) && (fabs(VertexPar->TkVxZOrig) <= VzBound))
    {
      unsigned int goodTrkCounter = 0;
      if (VertexPar->IsTrkPart == false)
	{
	  for (vector<TrackBaseRef>::const_iterator itTkTrack = VertexPar->TkVxReco->tracks_begin(); itTkTrack != VertexPar->TkVxReco->tracks_end(); itTkTrack++)
	    if (((*itTkTrack)->pt() >= MinPtVx) && (fabs((*itTkTrack)->eta()) <= MaxEtaVxTrk) && ((*itTkTrack)->normalizedChi2() <= MaxChi2VxTrk) &&
		(FindNpixelHits(NpixelLayers, &(**itTkTrack)) == true)) goodTrkCounter++;
	}
      else
	{
	  for (TrackingParticleRefVector::iterator itTkTrack = VertexPar->TkVxTP->daughterTracks_begin(); itTkTrack != VertexPar->TkVxTP->daughterTracks_end(); itTkTrack++)	    
	    if ((abs((*itTkTrack)->charge()) >= 1) && ((*itTkTrack)->pt() >= MinPtVx) && (fabs((*itTkTrack)->eta()) <= MaxEtaVxTrk) &&
		(FindNpixelHits(NpixelLayers, *itTkTrack) == true)) goodTrkCounter++;
	}
      
      if (goodTrkCounter >= MinVxTrkMatch)
	{
	  H_trk_vx_z->Fill(VertexPar->TkVxZOrig);	
	  H_trk_vx_xy->Fill(VertexPar->TkVxXOrig*10000., VertexPar->TkVxYOrig*10000.);
	  H_trk_vx_xerr->Fill(VertexPar->TkVxXerrOrig);
	  H_trk_vx_yerr->Fill(VertexPar->TkVxYerrOrig);
	  H_trk_vx_zerr->Fill(VertexPar->TkVxZerrOrig);
	  H_trk_vx_normchi2->Fill(VertexPar->TkVxChi2DoF);
	  H_trk_vx_dof->Fill(VertexPar->TkVxDoFOrig);
	  H_ntrk_trk_vx->Fill(VertexPar->TkVxNTracksOrig);
	  H_vx_counters->SetBinContent(H_vx_counters->GetXaxis()->FindBin("Strip"), H_vx_counters->GetBinContent(H_vx_counters->GetXaxis()->FindBin("Strip")) + 1.);

	  result = true;
	}
    }
  else if ((DataType.compare("nomatch") == 0) &&
	   (VertexPar->TkVxDoFOrig >= MinVxDoF) &&
	   (sqrt(VertexPar->TkVxXOrig*VertexPar->TkVxXOrig + VertexPar->TkVxYOrig*VertexPar->TkVxYOrig) <= VrBound) && (fabs(VertexPar->TkVxZOrig) <= VzBound))
    {
      H_trk_vx_z_nomatch->Fill(VertexPar->TkVxZOrig);
      H_trk_vx_xy_nomatch->Fill(VertexPar->TkVxXOrig*10000., VertexPar->TkVxYOrig*10000.);
      H_trk_vx_nomatch_normchi2->Fill(VertexPar->TkVxChi2DoF);
      H_trk_vx_nomatch_dof->Fill(VertexPar->TkVxDoFOrig);

      if (VertexPar->IsTrkPart == false)
	{
	  for (vector<TrackBaseRef>::const_iterator itTkTrack = VertexPar->TkVxReco->tracks_begin(); itTkTrack != VertexPar->TkVxReco->tracks_end(); itTkTrack++)
	    {
	      H_trk_vx_nomatch_trk_pt->Fill((*itTkTrack)->pt());
	      H_trk_vx_nomatch_trk_eta->Fill((*itTkTrack)->eta());
	    }
	}
      else
	{
	  for (TrackingParticleRefVector::iterator itTkTrack = VertexPar->TkVxTP->daughterTracks_begin(); itTkTrack != VertexPar->TkVxTP->daughterTracks_end(); itTkTrack++)	    
	    if (abs((*itTkTrack)->charge()) >= 1)
	      {
		H_trk_vx_nomatch_trk_pt->Fill((*itTkTrack)->pt());
		H_trk_vx_nomatch_trk_eta->Fill((*itTkTrack)->eta());
	      }
	}
      
      result = true;
    }
  
  return result;
}


void MyPixAnalyzer::PixelVerticesMC (const edm::Event& iEvent,
				     const edm::EventSetup& iSetup)
{
  edm::Handle<BeamSpot> beamSpotH;
  iEvent.getByLabel("offlineBeamSpot",beamSpotH);
  BeamSpot beamSpot = *beamSpotH;

  edm::Handle<edm::View<Vertex> > pixelVertexCollection;
  iEvent.getByLabel(VxInputTag,pixelVertexCollection);

  edm::Handle<edm::View<Track> > PxTracks;
  iEvent.getByLabel("pixelTracks",PxTracks);
  edm::RefToBase<Track> RecoTrk;

  edm::Handle<TrackingVertexCollection> TVCollection;
  iEvent.getByLabel(edm::InputTag("mergedtruth","MergedTrackTruth"),TVCollection);
  const TrackingVertexCollection TVC = *(TVCollection.product());
  TrackingVertexRef itTkVx;

  edm::Handle<TrackingParticleCollection> TPCollection;
  iEvent.getByLabel(edm::InputTag("mergedtruth","MergedTrackTruth"),TPCollection);
  RecoToSimCollection TkReco2Sim = theTkAssociator->associateRecoToSim(PxTracks, TPCollection, &iEvent);
  vector<pair<TrackingParticleRef, double> > TrkParticle;

  VertexParType VertexPar;

  vector<int>* MatchedTkTracksPix    = new vector<int>();
  vector<int>* MatchedTkTracksPixNew = new vector<int>();
  vector<int>* MatchedTkTracksTrk    = new vector<int>();
  vector<int>* MatchedTkTracksTrkNew = new vector<int>();
  vector<int>* MatchedTkTracksTmp;

  vector<unsigned int> MatchedTrkVx;

  int evId_event = -1;

  // ##################
  // # Pixel vertices #
  // ##################
  for (edm::View<Vertex>::const_iterator itPxVx = pixelVertexCollection->begin(); itPxVx != pixelVertexCollection->end(); itPxVx++) {

    if ((itPxVx->isValid() == true) && (itPxVx->isFake() == false))
      {
	VertexPar.PxVxNTracksOrig = itPxVx->tracksSize();
	VertexPar.PxVxDoFOrig     = itPxVx->ndof();
	VertexPar.PxVxChi2DoF     = itPxVx->normalizedChi2();
	VertexPar.PxVxX           = itPxVx->x() - beamSpot.position().x();
	VertexPar.PxVxXerr        = itPxVx->xError() * PxVxErrCorr;
	VertexPar.PxVxY           = itPxVx->y() - beamSpot.position().y();
	VertexPar.PxVxYerr        = itPxVx->yError() * PxVxErrCorr;
	VertexPar.PxVxZ           = itPxVx->z() - beamSpot.position().z();
	VertexPar.PxVxZerr        = itPxVx->zError() * PxVxErrCorr;

	FillVertexHistos("pixel", &VertexPar, iEvent, iSetup);

	int TkVxIndx = -1;
	unsigned int PixTrkMatched;
	unsigned int MaxPixTrkMatched = 0;

	// ##########################
	// # Loop over sim-vertices #
	// ##########################
	for (TrackingVertexCollection::size_type i = 0; i < TVC.size(); i++) {

	  itTkVx = TrackingVertexRef(TVCollection, i);

	  if (itTkVx->eventId().event() != evId_event)
	    {
	      // ###################################################################################################
	      // # The first event in the list is the signal, the others are pile-up                               #
	      // # The TrackingVertex and TrackingParticle have an eventID  that one can use to distinguish events #
	      // ###################################################################################################
	      // i = 0 --> eventId = 0 --> First Primary Vertex = Signal
	      // i = 1 --> eventId = 0
	      // i = 2 --> eventId = 0
	      // ...
	      // i = n --> eventId = 1 --> Second Primary Vertex = Pile Up
	      // i = n+1 --> eventId = 1
	      // i = n+2 --> eventId = 1
	      // ...
	      // i = m --> eventId = 2 --> Third Primary Vertex = Pile Up
	      // i = m+1 --> eventId = 2
	      // i = m+2 --> eventId = 2
	      // ...
	      evId_event = itTkVx->eventId().event();
	      PixTrkMatched = 0;
	      int TkTkIndx;

	      // ##############################
	      // # Loop over the pixel tracks #
	      // ##############################
	      for (vector<TrackBaseRef>::const_iterator itPxTrack = itPxVx->tracks_begin(); itPxTrack != itPxVx->tracks_end(); itPxTrack++) {

		RecoTrk = edm::RefToBase<Track>(*itPxTrack);
		TkTkIndx = -1;

		// ############################
		// # Loop over the sim-tracks #
		// ############################
		for (TrackingParticleRefVector::iterator itTkTrack = itTkVx->daughterTracks_begin(); itTkTrack != itTkVx->daughterTracks_end(); itTkTrack++) {	    	      	
    
		  // ##########################################
		  // # Search for the corresponding sim-track #
		  // ##########################################
		  if (TkReco2Sim.find(RecoTrk) != TkReco2Sim.end())
		    {
		      TrkParticle = TkReco2Sim[RecoTrk];

		      if ((TrkParticle.size() != 0) &&
			  (abs(TrkParticle.begin()->first->charge()) >= 1) &&
			  (TrkParticle.begin()->first.key() == itTkTrack->key())) { TkTkIndx = itTkTrack - itTkVx->daughterTracks_begin(); break; }
		    }
		} // End Loop Sim Tracks

		if (TkTkIndx >= 0)
		  {
		    TrackingParticleRefVector::iterator itTkTrack = itTkVx->daughterTracks_begin() + TkTkIndx;
		
		    if (find(MatchedTkTracksTrkNew->begin(), MatchedTkTracksTrkNew->end(), itTkTrack->key()) == MatchedTkTracksTrkNew->end())
		      {
			PixTrkMatched++;
			MatchedTkTracksTrkNew->push_back(itTkTrack->key());
		      }
		    MatchedTkTracksPixNew->push_back(itPxTrack->key());		
		  }
	      } // End Loop Pixel Tracks
	  
	      // ######################################################################
	      // # If there are more sim vertices matched with the same pixel vertex, #
	      // # then select only the one with the highest number of matched tracks #
	      // ######################################################################
	      if ((PixTrkMatched > MaxPixTrkMatched) && (PixTrkMatched >= MinVxTrkMatch))
		{
		  MaxPixTrkMatched = PixTrkMatched;
		  TkVxIndx = i;
	      
		  MatchedTkTracksTmp = MatchedTkTracksTrk;
		  MatchedTkTracksTrk = MatchedTkTracksTrkNew;
		  MatchedTkTracksTrkNew = MatchedTkTracksTmp;
		  MatchedTkTracksTrkNew->clear();
	      
		  MatchedTkTracksTmp = MatchedTkTracksPix;
		  MatchedTkTracksPix = MatchedTkTracksPixNew;
		  MatchedTkTracksPixNew = MatchedTkTracksTmp;
		  MatchedTkTracksPixNew->clear();
		}
	      else
		{
		  MatchedTkTracksTrkNew->clear();
		  MatchedTkTracksPixNew->clear();
		}
	    }
	} // End Loop Sim Vertices
	
	// ######################################################################
	// # If there are more pixel vertices matched with the same sim vertex, #
	// # then select only the first one                                     #
	// ######################################################################
	if ((TkVxIndx >= 0) && (find(MatchedTrkVx.begin(), MatchedTrkVx.end(), TkVxIndx) == MatchedTrkVx.end()))
	  {
	    itTkVx = TrackingVertexRef(TVCollection, TkVxIndx);

	    VertexPar.NTracksMatched     = MaxPixTrkMatched;
	    VertexPar.PxVxNTracks        = itPxVx->tracksSize();
	    VertexPar.PxVxNTracksOrig    = itPxVx->tracksSize();
	    VertexPar.TkVxNTracks        = CountSimChargedTracks(itTkVx);
	    VertexPar.TkVxNTracksOrig    = VertexPar.TkVxNTracks;
	    VertexPar.PxVxDoF            = itPxVx->ndof();
	    VertexPar.PxVxDoFOrig        = itPxVx->ndof();
	    VertexPar.TkVxDoF            = 2.*(double)VertexPar.TkVxNTracks-3.;
	    VertexPar.TkVxDoFOrig        = 2.*(double)VertexPar.TkVxNTracksOrig-3.;
	    VertexPar.PxVxChi2DoF        = itPxVx->normalizedChi2();
	    VertexPar.TkVxChi2DoF        = 0.;

	    VertexPar.PxVxX              = itPxVx->x() - beamSpot.position().x();
	    VertexPar.TkVxX              = itTkVx->position().x() - beamSpot.position().x();
	    VertexPar.TkVxXOrig          = itTkVx->position().x() - beamSpot.position().x();
	    VertexPar.PxVxXerr           = itPxVx->xError() * PxVxErrCorr;
	    VertexPar.TkVxXerr           = 0.;
	    VertexPar.TkVxXerrOrig       = 0.;

	    VertexPar.PxVxY              = itPxVx->y() - beamSpot.position().y();
	    VertexPar.TkVxY              = itTkVx->position().y() - beamSpot.position().y();
	    VertexPar.TkVxYOrig          = itTkVx->position().y() - beamSpot.position().y();
	    VertexPar.PxVxYerr           = itPxVx->yError() * PxVxErrCorr;
	    VertexPar.TkVxYerr           = 0.;
	    VertexPar.TkVxYerrOrig       = 0.;

	    VertexPar.PxVxZ              = itPxVx->z() - beamSpot.position().z();
	    VertexPar.TkVxZ              = itTkVx->position().z() - beamSpot.position().z();
	    VertexPar.TkVxZOrig          = itTkVx->position().z() - beamSpot.position().z();
	    VertexPar.PxVxZerr           = itPxVx->zError() * PxVxErrCorr;
	    VertexPar.TkVxZerr           = 0.;
	    VertexPar.TkVxZerrOrig       = 0.;

	    VertexPar.PxVx               = &(*itPxVx);
	    VertexPar.TkVxReco           = NULL;
	    VertexPar.PxVxReFit          = NULL;
	    VertexPar.TkVxReFit          = NULL;
	    VertexPar.TkVxTP             = itTkVx;
	    VertexPar.MatchedTkTracksPix = MatchedTkTracksPix;
	    VertexPar.MatchedTkTracksTrk = MatchedTkTracksTrk;
	    VertexPar.IsTrkPart          = true;

	    if (FillVertexHistos("match", &VertexPar, iEvent, iSetup) == true) MatchedTrkVx.push_back(TkVxIndx);
	  }

	MatchedTkTracksPix->clear();
	MatchedTkTracksTrk->clear();
	MatchedTkTracksPixNew->clear();
	MatchedTkTracksTrkNew->clear();
      }
  } // End Loop Pixel Vertices

  evId_event = -1;
  for (TrackingVertexCollection::size_type i = 0; i < TVC.size(); i++)
    {
      itTkVx = TrackingVertexRef(TVCollection, i);
      if (itTkVx->eventId().event() != evId_event)
	{
	  evId_event = itTkVx->eventId().event();
	  
	  VertexPar.TkVxNTracksOrig = CountSimChargedTracks(itTkVx);
	  VertexPar.TkVxDoFOrig     = 2.*(double)VertexPar.TkVxNTracksOrig-3.;
	  VertexPar.TkVxChi2DoF     = 0.;
	  VertexPar.TkVxXOrig       = itTkVx->position().x() - beamSpot.position().x();
	  VertexPar.TkVxXerrOrig    = 0.;
	  VertexPar.TkVxYOrig       = itTkVx->position().y() - beamSpot.position().y();
	  VertexPar.TkVxYerrOrig    = 0.;
	  VertexPar.TkVxZOrig       = itTkVx->position().z() - beamSpot.position().z();
	  VertexPar.TkVxZerrOrig    = 0.;
	  VertexPar.TkVxReco        = NULL;
	  VertexPar.TkVxTP          = itTkVx;
	  VertexPar.IsTrkPart       = true;

	  FillVertexHistos("strip", &VertexPar, iEvent, iSetup);
	  
	  if (find(MatchedTrkVx.begin(), MatchedTrkVx.end(), i) == MatchedTrkVx.end())
	    FillVertexHistos("nomatch", &VertexPar, iEvent, iSetup);
	}
    } // End Loop Sim Vertices
  
  delete MatchedTkTracksPix;
  delete MatchedTkTracksTrk;
  delete MatchedTkTracksPixNew;
  delete MatchedTkTracksTrkNew;
  MatchedTrkVx.clear();
}


void MyPixAnalyzer::PixelVerticesData (const edm::Event& iEvent,
				       const edm::EventSetup& iSetup)
{
  edm::Handle<BeamSpot> beamSpotH;
  iEvent.getByLabel("offlineBeamSpot",beamSpotH);
  BeamSpot beamSpot = *beamSpotH;

  edm::Handle<edm::View<Vertex> > offVertexCollection;
  iEvent.getByLabel("offlinePrimaryVertices",offVertexCollection);

  edm::Handle<edm::View<Vertex> > pixelVertexCollection;
  iEvent.getByLabel(VxInputTag,pixelVertexCollection);

  edm::Handle<edm::View<Track> > PxTracks;
  iEvent.getByLabel("pixelTracks",PxTracks);

  VertexParType VertexPar;

  stringstream Name;
  bool PixelHitMatched;
  vector<string> TrkNoMatchHitLocation;
  vector<string> PixMatchHitLocation;
  vector<string> PixNoMatchHitLocation;

  const Track* PxTrack;
  const Track* TkTrack;
  vector<int>* MatchedTkTracksPix    = new vector<int>();
  vector<int>* MatchedTkTracksPixNew = new vector<int>();
  vector<int>* MatchedTkTracksTrk    = new vector<int>();
  vector<int>* MatchedTkTracksTrkNew = new vector<int>();
  vector<int>* MatchedTkTracksTmp;

  vector<unsigned int> MatchedTrkVx;

  TRandom rnd;

  map<int, vector<int> >* MapMatchTrkPix    = new map<int, vector<int> >();
  map<int, vector<int> >* MapMatchTrkPixNew = new map<int, vector<int> >();
  map<int, vector<int> >* MapMatchTrkPixTmp;

  // ##################
  // # Pixel vertices #
  // ##################
  for (edm::View<Vertex>::const_iterator itPxVx = pixelVertexCollection->begin(); itPxVx != pixelVertexCollection->end(); itPxVx++) {
    
    if ((itPxVx->isValid() == true) && (itPxVx->isFake() == false))
      {
	VertexPar.PxVxNTracksOrig = itPxVx->tracksSize();
	VertexPar.PxVxDoFOrig     = itPxVx->ndof();
	VertexPar.PxVxChi2DoF     = itPxVx->normalizedChi2();
	VertexPar.PxVxX           = itPxVx->x() - beamSpot.position().x();
	VertexPar.PxVxXerr        = itPxVx->xError() * PxVxErrCorr;
	VertexPar.PxVxY           = itPxVx->y() - beamSpot.position().y();
	VertexPar.PxVxYerr        = itPxVx->yError() * PxVxErrCorr;
	VertexPar.PxVxZ           = itPxVx->z() - beamSpot.position().z();
	VertexPar.PxVxZerr        = itPxVx->zError() * PxVxErrCorr;

	FillVertexHistos("pixel", &VertexPar, iEvent, iSetup);

	int TkVxIndx = -1;
	unsigned int PixTrkMatched;
	unsigned int MaxPixTrkMatched = 0;
	double DoFTrkMatched;
	double MaxDoFTrkMatched = 0.;

	unsigned int FirstPixelRecPxClusterRow;
	unsigned int FirstPixelRecPxClusterCol;
	unsigned int LastPixelRecPxClusterRow;
	unsigned int LastPixelRecPxClusterCol;
	unsigned int IdPixelRecPxCluster;
	unsigned int PixTrkHitMatched;
	    
	// ####################
	// # Tracker vertices #
	// ####################
	for (edm::View<Vertex>::const_iterator itTkVx = offVertexCollection->begin(); itTkVx != offVertexCollection->end(); itTkVx++) {

	  PixTrkMatched = 0;
	  DoFTrkMatched = 0.;

	  if ((itTkVx->isValid() == true) && (itTkVx->isFake() == false))
	    {
	      int TkTkIndx;
	      double BestTrkTrkChi2;

	      // ##############################
	      // # Loop over the pixel tracks #
	      // ##############################
	      for (vector<TrackBaseRef>::const_iterator itPxTrack = itPxVx->tracks_begin(); itPxTrack != itPxVx->tracks_end(); itPxTrack++) {

		PxTrack = &(*(*itPxTrack));
		TkTkIndx = -1;
		BestTrkTrkChi2 = RangeChi2;

		// ################################
		// # Loop over the tracker tracks #
		// ################################
		for (vector<TrackBaseRef>::const_iterator itTkTrack = itTkVx->tracks_begin(); itTkTrack != itTkVx->tracks_end(); itTkTrack++) {

		  PixTrkHitMatched = 0;

		  // #############################################
		  // # Loop over the pixel hits of a pixel track #
		  // #############################################
		  for (trackingRecHit_iterator itPxHit = PxTrack->recHitsBegin(); itPxHit != PxTrack->recHitsEnd(); itPxHit++) {	
		
		    if ((*itPxHit)->isValid() == true)
		      {
			// Extract cluster information for the pixel hit
			const SiPixelRecHit* recHitPxPixel = dynamic_cast<const SiPixelRecHit*>((*itPxHit)->clone());
			SiPixelRecHit::ClusterRef const& PixelRecPxCluster = recHitPxPixel->cluster();
			FirstPixelRecPxClusterRow = PixelRecPxCluster->minPixelRow();
			FirstPixelRecPxClusterCol = PixelRecPxCluster->minPixelCol();
			LastPixelRecPxClusterRow = PixelRecPxCluster->maxPixelRow();
			LastPixelRecPxClusterCol = PixelRecPxCluster->maxPixelCol();
			IdPixelRecPxCluster = recHitPxPixel->geographicalId()();
		    
			TkTrack = &(*(*itTkTrack));
			unsigned int FirstPixelRecTkClusterRow;
			unsigned int FirstPixelRecTkClusterCol;
			unsigned int LastPixelRecTkClusterRow;
			unsigned int LastPixelRecTkClusterCol;
			unsigned int IdPixelRecTkCluster;

			PixelHitMatched = false;

			// ###############################################
			// # Loop over the pixel hits of a tracker track #
			// ###############################################
			for (trackingRecHit_iterator itTkHit = TkTrack->recHitsBegin(); itTkHit != TkTrack->recHitsEnd(); itTkHit++) {

			  if (((*itTkHit)->isValid() == true) &&
			      (((*itTkHit)->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel) ||
			       ((*itTkHit)->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)))
			    {
			      // Extract cluster information for the tracker hit
			      const SiPixelRecHit* recHitTkPixel = dynamic_cast<const SiPixelRecHit*>((*itTkHit)->clone());
			      SiPixelRecHit::ClusterRef const& PixelRecCluster = recHitTkPixel->cluster();
			      FirstPixelRecTkClusterRow = PixelRecCluster->minPixelRow();
			      FirstPixelRecTkClusterCol = PixelRecCluster->minPixelCol();
			      LastPixelRecTkClusterRow = PixelRecCluster->maxPixelRow();
			      LastPixelRecTkClusterCol = PixelRecCluster->maxPixelCol();
			      IdPixelRecTkCluster = recHitTkPixel->geographicalId()();
			  
			      if ((IdPixelRecTkCluster == IdPixelRecPxCluster) &&
				  (FirstPixelRecTkClusterRow == FirstPixelRecPxClusterRow) &&
				  (FirstPixelRecTkClusterCol == FirstPixelRecPxClusterCol) &&
				  (LastPixelRecTkClusterRow == LastPixelRecPxClusterRow) &&
				  (LastPixelRecTkClusterCol == LastPixelRecPxClusterCol))
				{
				  if (PrintMsg == true)
				    {
				      cout << "### Matched hit ###" << endl;
				      cout << "\tFirstPixelRecClusterRow: " << FirstPixelRecPxClusterRow << endl;
				      cout << "\tFirstPixelRecClusterCol: " << FirstPixelRecPxClusterCol << endl;
				      cout << "\tLastPixelRecClusterRow: " << LastPixelRecPxClusterRow << endl;
				      cout << "\tLastPixelRecClusterCol: " << LastPixelRecPxClusterCol << endl;
				      cout << "\tIdPixelRecCluster: " << IdPixelRecPxCluster << endl;
				    }
			      
				  PixelHitMatched = true;
				  PixTrkHitMatched++;
				  delete recHitTkPixel;
				  break;
				}
			      else
				{
				  Name.str("");
				  if ((*itTkHit)->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel)
				    {
				      Name << "L" << ((PXBDetId)((*itTkHit)->geographicalId())).layer();
				      TrkNoMatchHitLocation.push_back(Name.str());
				    }
				  else if ((*itTkHit)->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)
				    {
				      Name << "D" << ((PXFDetId)((*itTkHit)->geographicalId())).disk();
				      Name << "S" << ((PXFDetId)((*itTkHit)->geographicalId())).side();
				      TrkNoMatchHitLocation.push_back(Name.str());
				    }
				}

			      delete recHitTkPixel;
			    }
			} // End Loop on Pixel-hits of a Tracker Track
			    
			if (PixelHitMatched == true)
			  {				
			    Name.str("");
			    if ((*itPxHit)->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel)
			      {
				Name << "L" << ((PXBDetId)((*itPxHit)->geographicalId())).layer();
				PixMatchHitLocation.push_back(Name.str());
			      }
			    else if ((*itPxHit)->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)
			      {
				Name << "D" << ((PXFDetId)((*itPxHit)->geographicalId())).disk();
				Name << "S" << ((PXFDetId)((*itPxHit)->geographicalId())).side();
				PixMatchHitLocation.push_back(Name.str());
			      }
			  }
			else
			  {
			    Name.str("");
			    if ((*itPxHit)->geographicalId().subdetId() == (int)PixelSubdetector::PixelBarrel)
			      {
				Name << "L" << ((PXBDetId)((*itPxHit)->geographicalId())).layer();
				PixNoMatchHitLocation.push_back(Name.str());
			      }
			    else if ((*itPxHit)->geographicalId().subdetId() == (int)PixelSubdetector::PixelEndcap)
			      {
				Name << "D" << ((PXFDetId)((*itPxHit)->geographicalId())).disk();
				Name << "S" << ((PXFDetId)((*itPxHit)->geographicalId())).side();
				PixNoMatchHitLocation.push_back(Name.str());
			      }				    
			  }
			    
			delete recHitPxPixel;
		      }
		  } // End Loop on Pixel-hits of a Pixel Track
 
		  if (PixTrkHitMatched >= MinHitsMatch)
		    {
		      // ############################################################################
		      // # If there are hits of the tracker track that are not matched and if they  #
		      // # do not belong to the same Layer or Disk of the pixel matched hits and if #
		      // # there are not matched pixel hits, then do not perform the match          #
		      // ############################################################################
		      unsigned int i;
		      for (i = 0; i < TrkNoMatchHitLocation.size(); i++)
			if ((find(PixMatchHitLocation.begin(), PixMatchHitLocation.end(), TrkNoMatchHitLocation[i]) == PixMatchHitLocation.end()) && (PixNoMatchHitLocation.size() != 0)) break;
		      if (i == TrkNoMatchHitLocation.size())
			{
			  // #######################################################################
			  // # If there are more tracker tracks matched with the same pixel track, #
			  // # then select only the one with best chi2                             #
			  // #######################################################################
			  if ((*itTkTrack)->normalizedChi2() < BestTrkTrkChi2)
			    {
			      BestTrkTrkChi2 = (*itTkTrack)->normalizedChi2();
			      TkTkIndx = itTkTrack - itTkVx->tracks_begin();
			    }

			  // #############################################################################
			  // # Fill the map:                                                             #
			  // # [Track Key] --> <pixel-track key, pixel-track key, pixel-track key, ... > #
			  // # [Track Key] --> <pixel-track key, pixel-track key, pixel-track key, ... > #
			  // # .                                                                         #
			  // # .                                                                         #
			  // # .                                                                         #
			  // #############################################################################
			  if (MapMatchTrkPixNew->find(itTkTrack->key()) != MapMatchTrkPixNew->end())
			    {
			      if (find((*MapMatchTrkPixNew)[itTkTrack->key()].begin(), (*MapMatchTrkPixNew)[itTkTrack->key()].end(), itPxTrack->key()) == (*MapMatchTrkPixNew)[itTkTrack->key()].end())
				(*MapMatchTrkPixNew)[itTkTrack->key()].push_back(itPxTrack->key());
			    }
			  else (*MapMatchTrkPixNew)[itTkTrack->key()].push_back(itPxTrack->key());
			}
		    }

		  TrkNoMatchHitLocation.clear();
		  PixMatchHitLocation.clear();
		  PixNoMatchHitLocation.clear();
		} // End Loop Tracker Tracks

		if (TkTkIndx >= 0)
		  {
		    vector<TrackBaseRef>::const_iterator itTkTrack = itTkVx->tracks_begin() + TkTkIndx;

		    DoFTrkMatched = itTkVx->ndof();
		    if (find(MatchedTkTracksTrkNew->begin(), MatchedTkTracksTrkNew->end(), itTkTrack->key()) == MatchedTkTracksTrkNew->end())
		      {
			PixTrkMatched++;
			MatchedTkTracksTrkNew->push_back(itTkTrack->key());
		      }
		    MatchedTkTracksPixNew->push_back(itPxTrack->key());

		    if (PrintMsg == true)
		      {
			cout << "### Matched track ###" << endl;
			cout << "\tNumber of valid hits: " << (*itPxTrack)->numberOfValidHits() << ", number of lost hits: " << (*itPxTrack)->numberOfLostHits() << endl;
			cout << "\tChi2/DoF-pixel: " << (*itPxTrack)->normalizedChi2() << "\tChi2/DoF-track: " << (*itTkTrack)->normalizedChi2() << endl;
			cout << "\tpt-pixel: " << (*itPxTrack)->pt() << "\tpt-track: " << (*itTkTrack)->pt() << endl;
			cout << "\tdz-pixel: " << (*itPxTrack)->dz() << "\tdz-track:" << (*itTkTrack)->dz() << endl;
			cout << "\td0-pixel: " << -(*itPxTrack)->d0() << "\td0-track: " << -(*itTkTrack)->d0() << endl;
			cout << "\ttheta-pixel: " << (*itPxTrack)->theta() << "\ttheta-track:" << (*itTkTrack)->theta() << endl;
			cout << "\tphi-pixel: " << (*itPxTrack)->phi() << "\tphi-track: " << (*itTkTrack)->phi() << endl;
			cout << "\teta-pixel: " << (*itPxTrack)->eta() << "\teta-track: " << (*itTkTrack)->eta() << endl;
			(*itPxTrack)->hitPattern().print();
		      }
		  }
	      } // End Loop Pixel Tracks

	      // ###################################################################################
	      // # If there are more tracker vertices matched with the same pixel vertex,          #
	      // # then select only the one with the highest number of matched tracks and best DoF #
	      // ###################################################################################
	      if ((PixTrkMatched > MaxPixTrkMatched) && (DoFTrkMatched > MaxDoFTrkMatched) && (PixTrkMatched >= MinVxTrkMatch))
		{
		  MaxPixTrkMatched = PixTrkMatched;
		  MaxDoFTrkMatched = DoFTrkMatched;
		  TkVxIndx = itTkVx - offVertexCollection->begin();

		  MatchedTkTracksTmp = MatchedTkTracksTrk;
		  MatchedTkTracksTrk = MatchedTkTracksTrkNew;
		  MatchedTkTracksTrkNew = MatchedTkTracksTmp;
		  MatchedTkTracksTrkNew->clear();

		  MatchedTkTracksTmp = MatchedTkTracksPix;
		  MatchedTkTracksPix = MatchedTkTracksPixNew;
		  MatchedTkTracksPixNew = MatchedTkTracksTmp;
		  MatchedTkTracksPixNew->clear();

		  MapMatchTrkPixTmp = MapMatchTrkPix;
		  MapMatchTrkPix = MapMatchTrkPixNew;
		  MapMatchTrkPixNew = MapMatchTrkPixTmp;
		  for (map<int, vector<int> >::iterator itMap = MapMatchTrkPixNew->begin(); itMap != MapMatchTrkPixNew->end(); itMap++) itMap->second.clear();
		  MapMatchTrkPixNew->clear();
		}
	      else
		{
		  MatchedTkTracksTrkNew->clear();
		  MatchedTkTracksPixNew->clear();
		  for (map<int, vector<int> >::iterator itMap = MapMatchTrkPixNew->begin(); itMap != MapMatchTrkPixNew->end(); itMap++) itMap->second.clear();
		  MapMatchTrkPixNew->clear();
		}
	    }
	} // End Loop Tracker Vertices

	// ##########################################################################
	// # If there are more pixel vertices matched with the same tracker vertex, #
	// # then select only the first one                                         #
	// ##########################################################################
	if ((TkVxIndx >= 0) && (find(MatchedTrkVx.begin(), MatchedTrkVx.end(), TkVxIndx) == MatchedTrkVx.end()))
	  {
	    edm::View<Vertex>::const_iterator itTkVx = offVertexCollection->begin() + TkVxIndx;

	    // ####################################################################
	    // # Find the tracks of the full-tracking vertex that need to be kept #
	    // # Fing the tracks of the pixel-vertex that need to be removed      #
	    // ####################################################################
	    vector<TrackBaseRef> GoodFullTracks;
	    vector<int> BadPixelTracks;
	    for (vector<TrackBaseRef>::const_iterator itVxTkTrack = itTkVx->tracks_begin(); itVxTkTrack != itTkVx->tracks_end(); itVxTkTrack++)
	      if (rnd.Uniform() >= 0.5)
		{
		  GoodFullTracks.push_back(*itVxTkTrack);
		  if (MapMatchTrkPix->find(itVxTkTrack->key()) != MapMatchTrkPix->end())
		    for (unsigned int i = 0; i < (*MapMatchTrkPix)[itVxTkTrack->key()].size(); i++)
		      if (find(BadPixelTracks.begin(), BadPixelTracks.end(), (*MapMatchTrkPix)[itVxTkTrack->key()][i]) == BadPixelTracks.end())
			BadPixelTracks.push_back((*MapMatchTrkPix)[itVxTkTrack->key()][i]);
		}
		
	    // ############################################################
	    // # Find the tracks of the pixel vertex that need to be kept #
	    // ############################################################
	    vector<TrackBaseRef> GoodPixelTracks;
	    for (vector<TrackBaseRef>::const_iterator itVxPxTrack = itPxVx->tracks_begin(); itVxPxTrack != itPxVx->tracks_end(); itVxPxTrack++)
	      if (find(BadPixelTracks.begin(), BadPixelTracks.end(), itVxPxTrack->key()) == BadPixelTracks.end())
		GoodPixelTracks.push_back(*itVxPxTrack);
		
	    // ###############################
	    // # Variables for vertex re-fit #
	    // ###############################
	    const edm::Provenance* prov;
	    edm::ParameterSetID psid;
	    edm::pset::Registry* psregistry;
	    edm::ParameterSet psetFromProvenance;
	    edm::ESHandle<TransientTrackBuilder> TTBuilder;
	    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);
	    vector<TransientTrack> TTrks;
	    auto_ptr<PrimaryVertexProducerAlgorithm> VxProdAlgo;
	    vector<TransientVertex> MyPrimVx;
	    vector<TransientVertex> MyPixVx;

	    // #####################################################################
	    // # Obtain the configuration of the full-tracking vertexing algorithm #
	    // #####################################################################
	    prov       = offVertexCollection.provenance();
	    psid       = prov->psetID();
	    psregistry = edm::pset::Registry::instance();
	    psregistry->getMapped(psid, psetFromProvenance);
	    TTrks.reserve(GoodFullTracks.size());

	    for (vector<TrackBaseRef>::const_iterator it = GoodFullTracks.begin(); it != GoodFullTracks.end(); it++)
	      {
		TTrks.push_back((*TTBuilder).build(*(edm::RefToBase<Track>(*it))));
		TTrks.back().setBeamSpot(*beamSpotH);
	      }

	    // ###################################
	    // # Re-fit the full-tracking vertex #
	    // ###################################
	    VxProdAlgo.reset(new PrimaryVertexProducerAlgorithm(psetFromProvenance));
	    MyPrimVx = VxProdAlgo->vertices(TTrks, *beamSpotH);

	    TTrks.clear();

	    // #############################################################
	    // # Obtain the configuration of the pixel vertexing algorithm #
	    // #############################################################
	    prov       = pixelVertexCollection.provenance();
	    psid       = prov->psetID();
	    psregistry = edm::pset::Registry::instance();
	    psregistry->getMapped(psid, psetFromProvenance);
	    TTrks.reserve(GoodPixelTracks.size());

	    for (vector<TrackBaseRef>::const_iterator it = GoodPixelTracks.begin(); it != GoodPixelTracks.end(); it++)
	      {
		TTrks.push_back((*TTBuilder).build(*(edm::RefToBase<Track>(*it))));
		TTrks.back().setBeamSpot(*beamSpotH);
	      }

	    // ###########################
	    // # Re-fit the pixel vertex #
	    // ###########################
	    VxProdAlgo.reset(new PrimaryVertexProducerAlgorithm(psetFromProvenance));
	    MyPixVx = VxProdAlgo->vertices(TTrks, *beamSpotH);

	    VertexPar.NTracksMatched     = MaxPixTrkMatched;
	    VertexPar.PxVxNTracks        = 0;
	    VertexPar.PxVxNTracksOrig    = itPxVx->tracksSize();
	    VertexPar.TkVxNTracks        = 0;
	    VertexPar.TkVxNTracksOrig    = itTkVx->tracksSize();
	    VertexPar.PxVxDoF            = 0.;
	    VertexPar.PxVxDoFOrig        = itPxVx->ndof();
	    VertexPar.TkVxDoF            = 0.;
	    VertexPar.TkVxDoFOrig        = itTkVx->ndof();
	    VertexPar.PxVxChi2DoF        = itPxVx->normalizedChi2();
	    VertexPar.TkVxChi2DoF        = itTkVx->normalizedChi2();

	    VertexPar.PxVxX              = 0.;
	    VertexPar.TkVxX              = 0.;
	    VertexPar.TkVxXOrig          = itTkVx->x() - beamSpot.position().x();
	    VertexPar.PxVxXerr           = 0.;
	    VertexPar.TkVxXerr           = 0.;
	    VertexPar.TkVxXerrOrig       = itTkVx->xError() * TkXYVxErrCorr;

	    VertexPar.PxVxY              = 0.;
	    VertexPar.TkVxY              = 0.;
	    VertexPar.TkVxYOrig          = itTkVx->y() - beamSpot.position().y();
	    VertexPar.PxVxYerr           = 0.;
	    VertexPar.TkVxYerr           = 0.;
	    VertexPar.TkVxYerrOrig       = itTkVx->yError() * TkXYVxErrCorr;

	    VertexPar.PxVxZ              = 0.;
	    VertexPar.TkVxZ              = 0.;
	    VertexPar.TkVxZOrig          = itTkVx->z() - beamSpot.position().z();
	    VertexPar.PxVxZerr           = 0.;
	    VertexPar.TkVxZerr           = 0.;
	    VertexPar.TkVxZerrOrig       = itTkVx->zError() * TkZVxErrCorr;

	    VertexPar.PxVx               = &(*itPxVx);
	    VertexPar.TkVxReco           = &(*itTkVx);
	    VertexPar.PxVxReFit          = NULL;
	    VertexPar.TkVxReFit          = NULL;
	    VertexPar.MatchedTkTracksPix = MatchedTkTracksPix;
	    VertexPar.MatchedTkTracksTrk = MatchedTkTracksTrk;
	    VertexPar.IsTrkPart          = false;

	    if ((MyPrimVx.size() == 1) && (Vertex(MyPrimVx.front()).isValid() == true) && (Vertex(MyPrimVx.front()).isFake() == false) &&
		(MyPixVx.size() == 1) && (Vertex(MyPixVx.front()).isValid() == true) && (Vertex(MyPixVx.front()).isFake() == false))
	      {
		VertexPar.PxVxNTracks    = Vertex(MyPixVx.front()).tracksSize();
		VertexPar.TkVxNTracks    = Vertex(MyPrimVx.front()).tracksSize();
		VertexPar.PxVxDoF        = Vertex(MyPixVx.front()).ndof();
		VertexPar.TkVxDoF        = Vertex(MyPrimVx.front()).ndof();

		VertexPar.PxVxX          = Vertex(MyPixVx.front()).x() - beamSpot.position().x();
		VertexPar.PxVxXerr       = Vertex(MyPixVx.front()).xError() * PxVxErrCorr;
		VertexPar.TkVxX          = Vertex(MyPrimVx.front()).x() - beamSpot.position().x();
		VertexPar.TkVxXerr       = Vertex(MyPrimVx.front()).xError() * TkXYVxErrCorr;

		VertexPar.PxVxY          = Vertex(MyPixVx.front()).y() - beamSpot.position().y();
		VertexPar.PxVxYerr       = Vertex(MyPixVx.front()).yError() * PxVxErrCorr;
		VertexPar.TkVxY          = Vertex(MyPrimVx.front()).y() - beamSpot.position().y();
		VertexPar.TkVxYerr       = Vertex(MyPrimVx.front()).yError() * TkXYVxErrCorr;

		VertexPar.PxVxZ          = Vertex(MyPixVx.front()).z() - beamSpot.position().z();
		VertexPar.PxVxZerr       = Vertex(MyPixVx.front()).zError() * PxVxErrCorr;
		VertexPar.TkVxZ          = Vertex(MyPrimVx.front()).z() - beamSpot.position().z();
		VertexPar.TkVxZerr       = Vertex(MyPrimVx.front()).zError() * TkZVxErrCorr;

		VertexPar.PxVxReFit      = &MyPixVx.front();
		VertexPar.TkVxReFit      = &MyPrimVx.front();
	      }

	    if (FillVertexHistos("match", &VertexPar, iEvent, iSetup) == true) MatchedTrkVx.push_back(TkVxIndx);

	    GoodFullTracks.clear();
	    GoodPixelTracks.clear();
	    BadPixelTracks.clear();
	    MyPrimVx.clear();
	    MyPixVx.clear();
	    TTrks.clear();
	  }

	MatchedTkTracksPix->clear();
	MatchedTkTracksTrk->clear();
	MatchedTkTracksPixNew->clear();
	MatchedTkTracksTrkNew->clear();
	for (map<int, vector<int> >::iterator itMap = MapMatchTrkPixNew->begin(); itMap != MapMatchTrkPixNew->end(); itMap++) itMap->second.clear();
	MapMatchTrkPixNew->clear();
      }
  } // End Loop Pixel Vertices

  for (edm::View<Vertex>::const_iterator itTkVx = offVertexCollection->begin(); itTkVx != offVertexCollection->end(); itTkVx++)
    {
      if ((itTkVx->isValid() == true) && (itTkVx->isFake() == false))
	{
	  VertexPar.TkVxNTracksOrig = itTkVx->tracksSize();
	  VertexPar.TkVxDoFOrig     = itTkVx->ndof();
	  VertexPar.TkVxChi2DoF     = itTkVx->normalizedChi2();
	  VertexPar.TkVxXOrig       = itTkVx->position().x() - beamSpot.position().x();
	  VertexPar.TkVxXerrOrig    = itTkVx->xError() * TkXYVxErrCorr;
	  VertexPar.TkVxYOrig       = itTkVx->position().y() - beamSpot.position().y();
	  VertexPar.TkVxYerrOrig    = itTkVx->yError() * TkXYVxErrCorr;
	  VertexPar.TkVxZOrig       = itTkVx->position().z() - beamSpot.position().z();
	  VertexPar.TkVxZerrOrig    = itTkVx->zError() * TkZVxErrCorr;
	  VertexPar.TkVxReco        = &(*itTkVx);
	  VertexPar.IsTrkPart       = false;
	  
	  FillVertexHistos("strip", &VertexPar, iEvent, iSetup);

	  if (find(MatchedTrkVx.begin(), MatchedTrkVx.end(), itTkVx - offVertexCollection->begin()) == MatchedTrkVx.end())
	    FillVertexHistos("nomatch", &VertexPar, iEvent, iSetup);
	  
	}
    } // End Loop Tracker Vertices

  delete MatchedTkTracksPix;
  delete MatchedTkTracksTrk;
  delete MatchedTkTracksPixNew;
  delete MatchedTkTracksTrkNew;
  delete MapMatchTrkPix;
  delete MapMatchTrkPixNew;
  MatchedTrkVx.clear();
}


// #######################################
// # @@@ METHODS for EVENT SELECTION @@@ #
// #######################################

bool MyPixAnalyzer::L1Analyzer (const edm::Event& iEvent,
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

  // ################################################
  // # Extract the L1 tigger menu: algo & technical #
  // ################################################
  iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd);
  L1Menu = menuRcd.product();

  // Exract the algo & technical masks
  iSetup.get<L1GtTriggerMaskAlgoTrigRcd>().get(l1GtTmAlgo);
  iSetup.get<L1GtTriggerMaskTechTrigRcd>().get(l1GtTmTech);

  // #################################
  // # Algo & technical trigger bits #
  // #################################
  edm::Handle<L1GlobalTriggerReadoutRecord> gtRecord;
  iEvent.getByLabel("gtDigis",gtRecord);

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
      cout << "@@@ Trigger Masks and Prescale Factors for Algo-bits @@@" << endl;
      for (vector<bool>::iterator itBit = AlgodWord.begin(); itBit != AlgodWord.end(); itBit++)
	{
	  cout << "Index algo-bits #" << itBit-AlgodWord.begin() << "\tMask: " << ((l1GtTmAlgo.product()->gtTriggerMask()[itBit-AlgodWord.begin()] & (1 << 0)) ? '1' : '0');
	  cout << "\tPrescale: " << (*PFAlgoTrig)[itBit-AlgodWord.begin()] << "\tValue: " << (*itBit ? '1' : '0') << endl;
	}

      cout << "@@@ Trigger Masks and Prescale Factors for Tech-bits @@@" << endl;
      for (vector<bool>::const_iterator itBit = TechWord.begin(); itBit != TechWord.end(); itBit++)
	{
	  cout << "Index tech-bits #" << itBit-TechWord.begin() << "\tMask: " << ((l1GtTmTech.product()->gtTriggerMask()[itBit-TechWord.begin()] & (1 << 0)) ? '1' : '0');
	  cout << "\tPrescale: " <<  (*PFTechTrig)[itBit-TechWord.begin()] << "\tValue: " << (*itBit ? '1' : '0') << endl;
	}
    }
  else
    {
      // ###################################
      // # Apply the mask to the algo bits #
      // ###################################
      for (vector<bool>::iterator itBit = AlgodWord.begin(); itBit != AlgodWord.end(); itBit++)
	if ((l1GtTmAlgo.product()->gtTriggerMask()[itBit-AlgodWord.begin()] & (1 << 0)) == true) *itBit = false;
      
      // ###########################
      // # Loop over the algo bits #
      // ###########################
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
      
      // ### Temporary commented for 2009 data taking ###
      //       // Apply the mask to the techinical bits
      //       for (vector<bool>::iterator itBit = TechWord.begin(); itBit != TechWord.end(); itBit++)
      // 	if ((l1GtTmTech.product()->gtTriggerMask()[itBit-TechWord.begin()] & (1 << 0)) == true) *itBit = false;
      // ### Temporary commented for 2009 data taking ###
      
      // ################################
      // # Loop over the technical bits #
      // ################################
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
      TechPassed = false;
      if (!(TechWord[36] || TechWord[37] || TechWord[38] || TechWord[39]) && (TechWord[40] || TechWord[41]) && !(XOR(TechWord[42],TechWord[43]))) TechPassed = true;
      // Veto Halo + Minimum Bias Selection + Veto Single Splash
      // ### Temporary condition for 2009 data taking ###

      if (PrintMsg == true)
	{
	  gtRecord->printGtDecision(cout);
	  gtRecord->printTechnicalTrigger(cout);
	}
    }

  return (AlgoPassed && TechPassed) ? true : false;
}


bool MyPixAnalyzer::HLTAnalyzer (const edm::Event& iEvent,
				 vector<string> MyHLTMask,
				 const bool PrintMsg)
{
  edm::TriggerNames TrigNames;
  edm::Handle<edm::TriggerResults> HLTResults;

  unsigned int Count = 0;

  // ###########################
  // # Extract the HLT results #
  // ###########################
  iEvent.getByLabel(edm::InputTag("TriggerResults","","HLT"),HLTResults);

  if ((HLTResults.isValid() == true) && (HLTResults->size() > 0))
    {
      TrigNames = iEvent.triggerNames(*HLTResults);
      
      for (unsigned int i = 0; i < TrigNames.triggerNames().size(); i++)
	{
	  for (unsigned int i = 0; i < MyHLTMask.size(); i++)
	    {
	      if ((TrigNames.triggerName(i) == MyHLTMask[i]) &&
		  (HLTResults->wasrun(TrigNames.triggerIndex(TrigNames.triggerName(i))) == true) &&
		  (HLTResults->accept(TrigNames.triggerIndex(TrigNames.triggerName(i))) == true) &&
		  (HLTResults->error(TrigNames.triggerIndex(TrigNames.triggerName(i))) == false))
		{
		  if (PrintMsg == true) cout << __LINE__ << "\t" << __PRETTY_FUNCTION__ << "\tObjet passed my HLT trigger mask: " << TrigNames.triggerName(i) << endl;
		  Count++;
		  break;
		}
	      else if (PrintMsg == true) cout << __LINE__ << "\t" << __PRETTY_FUNCTION__ << "\tObjet didn't passed my HLT trigger mask: " << TrigNames.triggerName(i) << endl;
	    }
	}

      // ###########################
      // # AND of all the HLT bits #
      // ###########################
      //       if (Count == MyHLTMask.size()) return true;

      // ###########################
      // # OR of all the HLT bits #
      // ###########################
      if (Count != 0) return true;
    }

  return false;
}


bool MyPixAnalyzer::EventSelection (const edm::Event& iEvent,
				    const edm::EventSetup& iSetup)
{
  edm::Handle<TrackCollection> TkTracks;
  iEvent.getByLabel("generalTracks",TkTracks);

  if ((TkTracks->size() < MinTkTracks) || (TkTracks->size() > MaxTkTracks)) return false;
  if ((IsMC == true) || (L1Analyzer (iEvent, iSetup, MyAlgoMask, MyTechMask, PrintMsg, false) == true)) return true;
  
  return false;
}


// ###########################
// # @@@ GENERAL METHODS @@@ #
// ###########################

void MyPixAnalyzer::analyze(const edm::Event& iEvent,
			    const edm::EventSetup& iSetup)
{
  if (PrintMsg == true)
    {
      edm::Handle<BeamSpot> beamSpotH;
      iEvent.getByLabel("offlineBeamSpot",beamSpotH);
      BeamSpot beamSpot = *beamSpotH;
      cout << "@@@@@@ Beam Spot Values @@@@@@" << endl;
      cout << "Beam Spot X --> " << beamSpot.x0() << "\tBeam Spot Width X --> " << beamSpot.BeamWidthX() << endl;
      cout << "Beam Spot Y --> " << beamSpot.y0() << "\tBeam Spot Width Y --> " << beamSpot.BeamWidthY() << endl;
      cout << "Beam Spot Z --> " << beamSpot.z0() << "\tBeam Spot Sigma Z --> " << beamSpot.sigmaZ() << endl;
    }

  EvTotal++;

  //   if (IsMC == false)
  //     {
  //       if (HLTAnalyzer (iEvent, MyHLTMask, PrintMsg)) EvPassHLT++;
  //       else EvDiscardHLT++;
      
  //       if (L1Analyzer (iEvent, iSetup, MyAlgoMask, MyTechMask, PrintMsg, false)) EvPassL1++;
  //       else EvDiscardL1++;
  //     }
    
  //   if (EventSelection(iEvent, iSetup) == true)
  //     {
  if (PrintMsg == true)
    {
      cout << "\n\n@@@@@@ SELECTED EVENT STARTED @@@@@@" << endl;
      cout << "\tRun: " <<  iEvent.id().run();
      cout << "\tLumisection: " << iEvent.luminosityBlock();
      cout << "\tOrbit: " << iEvent.orbitNumber();
      cout << "\tBX: " << iEvent.bunchCrossing();
      cout << "\tEvent ID: " << iEvent.id();
      cout << endl;
    }

  EvPassTotal++;
	
  if ((AnalysisType.compare("TrackANDVertex") == 0) || (AnalysisType.compare("Track") == 0))
    {
      if (PrintMsg == true) cout << "\n@@@ PIXEL TRACK ANALYSIS @@@" << endl;
      if (IsTrkPart == true) PixelTracksMC(iEvent,iSetup);
      else PixelTracksData(iEvent,iSetup);
    }

  if ((AnalysisType.compare("TrackANDVertex") == 0) || (AnalysisType.compare("Vertex") == 0))
    {
      if (PrintMsg == true) cout << "\n@@@ PIXEL VERTEX ANALYSIS @@@" << endl;
      if (IsTrkPart == true) PixelVerticesMC(iEvent,iSetup);
      else PixelVerticesData(iEvent,iSetup);
    }

  if (PrintMsg == true) cout << "\n@@@@@@ EVENT FINISHED @@@@@@" << endl;
  //     }
  //   else
  //     EvDiscardTotal++;
}


void MyPixAnalyzer::beginJob ()
{
  if (IsTrkPart == true) theTkAssociator = new TrackAssociatorByHits(ParSetTrackAss);

  SubDirPixTk.reset(new TFileDirectory(FileService->mkdir("PixelTracks")));
  SubDirPixTkCom.reset(new TFileDirectory(SubDirPixTk->mkdir("CommonPxTkPlots")));
  SubDirPixTkFit.reset(new TFileDirectory(FileService->mkdir("TrackFit")));

  SubDirPixVx.reset(new TFileDirectory(FileService->mkdir("PixelVertices")));
  SubDirPixVxNoM.reset(new TFileDirectory(SubDirPixVx->mkdir("NoMatchedVx")));
  SubDirPixVxCom.reset(new TFileDirectory(SubDirPixVx->mkdir("CommonPxTkPlots")));
  SubDirPixVxNtrks.reset(new TFileDirectory(FileService->mkdir("VxNtrks")));

  H_counters = FileService->make<TH1F>("Histo counters","Histo counters", 7, 0., 7.);
  H_counters->SetXTitle("");
  H_counters->SetYTitle("Entries (#)");
  H_counters->GetXaxis()->SetBinLabel(1,"Total events");
  H_counters->GetXaxis()->SetBinLabel(2,"L1 pass");
  H_counters->GetXaxis()->SetBinLabel(3,"L1 discard");
  H_counters->GetXaxis()->SetBinLabel(4,"HLT pass");
  H_counters->GetXaxis()->SetBinLabel(5,"HLT discard");
  H_counters->GetXaxis()->SetBinLabel(6,"Total pass");
  H_counters->GetXaxis()->SetBinLabel(7,"Total discard");

  // ################
  // # Pixel Tracks #
  // ################
  hTitle[0] = "NumE/AE";
  hTitle[1] = "NumP";
  hTitle[2] = "DenP";

  for (int i = 0; i < 3; i++)
    {
      stringstream i_str;
      i_str << i + 1;
      TString hName;

      hName = "Pixeltracks pt distribution ver" + i_str.str();
      H_pix_trk_pt[i] = SubDirPixTk->make<TH1F>("Pixeltracks pt distribution","Pixel-tracks p_{t} distribution" + hTitle[i], rint(RangePt/PtStep), 0., RangePt);
      H_pix_trk_pt[i]->SetXTitle("Track p_{t} (GeV/c)");
      H_pix_trk_pt[i]->SetYTitle("Entries (#)");

      hName = "Pixeltracks chi2DoF distribution ver" + i_str.str();
      H_pix_trk_normchi2[i] = SubDirPixTk->make<TH1F>("Pixeltracks chi2DoF distribution","Pixel-tracks \\chi^{2}/DoF distribution" + hTitle[i], 100, 0., 30.);
      H_pix_trk_normchi2[i]->SetXTitle("\\chi^{2}/DoF");
      H_pix_trk_normchi2[i]->SetYTitle("Entries (#)");

      hName = "Pixeltracks dz distribution ver" + i_str.str();
      H_pix_trk_dz[i] = SubDirPixTk->make<TH1F>("Pixeltracks dz distribution","Pixel-tracks d_{z} distribution" + hTitle[i], 200, -d0Range/2., d0Range/2.);
      H_pix_trk_dz[i]->SetXTitle("Track d_{z} (cm)");
      H_pix_trk_dz[i]->SetYTitle("Entries (#)");

      hName = "Pixeltracks d0 distribution ver" + i_str.str();
      H_pix_trk_d0[i] = SubDirPixTk->make<TH1F>("Pixeltracks d0 distribution","Pixel-tracks d_{0} distribution" + hTitle[i], 200, -dzRange/2., dzRange/2.);
      H_pix_trk_d0[i]->SetXTitle("Track d_{0} (cm)");
      H_pix_trk_d0[i]->SetYTitle("Entries (#)");

      hName = "Pixeltracks eta distribution ver" + i_str.str();
      H_pix_trk_eta[i] = SubDirPixTk->make<TH1F>("Pixeltracks eta distribution","Pixel-tracks \\eta distribution" + hTitle[i],  rint(2.*RangeEta/EtaStep), -RangeEta, RangeEta);
      H_pix_trk_eta[i]->SetXTitle("Track \\eta");
      H_pix_trk_eta[i]->SetYTitle("Entries (#)");

      hName = "Pixeltracks theta distribution ver" + i_str.str();
      H_pix_trk_theta[i] = SubDirPixTk->make<TH1F>("Pixeltracks theta distribution","Pixel-tracks \\theta distribution" + hTitle[i], 100, 0., 3.5);
      H_pix_trk_theta[i]->SetXTitle("Track \\theta (rad)");
      H_pix_trk_theta[i]->SetYTitle("Entries (#)");

      hName = "Pixeltracks phi distribution ver" + i_str.str();
      H_pix_trk_phi[i] = SubDirPixTk->make<TH1F>("Pixeltracks phi distribution","Pixel-tracks \\phi distribution" + hTitle[i], 100, -190., 190.);
      H_pix_trk_phi[i]->SetXTitle("Track \\phi (deg)");
      H_pix_trk_phi[i]->SetYTitle("Entries (#)");
    }

  H_pix_trk_doubleCounter = SubDirPixTk->make<TH1F>("Doublecounter pixeltracks","Double-counter pixel-tracks", 10, 0.5, 10.5);
  H_pix_trk_doubleCounter->SetXTitle("Counts (#)");
  H_pix_trk_doubleCounter->SetYTitle("Entries (#)");

  H_pix_trk_doubleCounter_eta = SubDirPixTk->make<TH1F>("Doublecounter pixeltracks eta distribution","Doublecounter pixel-tracks \\eta distribution", rint(2.*RangeEta/EtaStep), -RangeEta, RangeEta);
  H_pix_trk_doubleCounter_eta->SetXTitle("Track \\eta");
  H_pix_trk_doubleCounter_eta->SetYTitle("Entries (#)");

  H_pix_trk_doubleCounter_pt = SubDirPixTk->make<TH1F>("Doublecounter pixeltracks pt distribution","Doublecounter pixel-tracks p_{t} distribution", rint(RangePt/PtStep), 0., RangePt);
  H_pix_trk_doubleCounter_pt->SetXTitle("Track p_{t} (GeV/c)");
  H_pix_trk_doubleCounter_pt->SetYTitle("Entries (#)");

  H_pix_trk_doubleCounter_eta_rel = SubDirPixTk->make<TH1F>("Doublecounter pixeltracks eta relative distribution","Doublecounter pixel-tracks \\eta relative distribution", rint(2.*RangeEta/EtaStep), -RangeEta, RangeEta);
  H_pix_trk_doubleCounter_eta_rel->SetXTitle("Track \\eta");
  H_pix_trk_doubleCounter_eta_rel->SetYTitle("Entries (#)");

  H_pix_trk_doubleCounter_pt_rel = SubDirPixTk->make<TH1F>("Doublecounter pixeltracks pt relative distribution","Doublecounter pixel-tracks p_{t} relative distribution", rint(RangePt/PtStep), 0., RangePt);
  H_pix_trk_doubleCounter_pt_rel->SetXTitle("Track p_{t} (GeV/c)");
  H_pix_trk_doubleCounter_pt_rel->SetYTitle("Entries (#)");

  H_pix_trk_eta_phi = SubDirPixTk->make<TH2F>("Pixeltracks distribution in eta and phi","Pixel-track distribution in \\eta and \\phi", rint(2.*RangeEta/EtaStep), -RangeEta, RangeEta, rint(2.*RangePhi/PhiStep), -RangePhi, RangePhi);
  H_pix_trk_eta_phi->SetXTitle("Track \\eta");
  H_pix_trk_eta_phi->SetYTitle("Track \\phi (deg)");
  H_pix_trk_eta_phi->SetZTitle("Entries (#)");

  H_pix_trk_eta_pt = SubDirPixTk->make<TH2F>("Pixeltracks distribution in eta and pt","Pixel-track distribution in \\eta and p_{t}", rint(2.*RangeEta/EtaStep), -RangeEta, RangeEta, rint(RangePt/PtStep), 0., RangePt);
  H_pix_trk_eta_pt->SetXTitle("Track \\eta");
  H_pix_trk_eta_pt->SetYTitle("Track p_{t} (GeV/c)");
  H_pix_trk_eta_pt->SetZTitle("Entries (#)");

  // ##################
  // # Tracker Tracks #
  // ##################
  hTitle[0] = "NumE/AE";
  hTitle[1] = "NumP";
  hTitle[2] = "DenE";
  hTitle[3] = "DenAE";

  for (int i = 0; i < 4; i++)
    {
      stringstream i_str;
      i_str << i + 1;
      TString hName;

      hName = "Trackertracks pt distribution ver" + i_str.str();
      H_trk_trk_pt[i] = SubDirPixTk->make<TH1F>("Trackertracks pt distribution","Tracker-tracks p_{t} distribution" + hTitle[i], rint(RangePt/PtStep), 0., RangePt);
      H_trk_trk_pt[i]->SetXTitle("Track p_{t} (GeV/c)");
      H_trk_trk_pt[i]->SetYTitle("Entries (#)");

      hName = "Trackertracks chi2DoF distribution ver" + i_str.str();
      H_trk_trk_normchi2[i] = SubDirPixTk->make<TH1F>("Trackertracks chi2DoF distribution","Tracker-tracks \\chi^{2}/DoF distribution" + hTitle[i], 100, 0., 30.);
      H_trk_trk_normchi2[i]->SetXTitle("\\chi^{2}/DoF");
      H_trk_trk_normchi2[i]->SetYTitle("Entries (#)");

      hName = "Trackertracks dz distribution ver" + i_str.str();
      H_trk_trk_dz[i] = SubDirPixTk->make<TH1F>("Trackertracks dz distribution","Tracker-tracks d_{z} distribution" + hTitle[i], 200, -d0Range/2., d0Range/2.);
      H_trk_trk_dz[i]->SetXTitle("Track d_{z} (cm)");
      H_trk_trk_dz[i]->SetYTitle("Entries (#)");

      hName = "Trackertracks d0 distribution ver" + i_str.str();
      H_trk_trk_d0[i] = SubDirPixTk->make<TH1F>("Trackertracks d0 distribution","Tracker-tracks d_{0} distribution" + hTitle[i], 200, -dzRange/2., dzRange/2.);
      H_trk_trk_d0[i]->SetXTitle("Track d_{0} (cm)");
      H_trk_trk_d0[i]->SetYTitle("Entries (#)");

      hName = "Trackertracks eta distribution ver" + i_str.str();
      H_trk_trk_eta[i] = SubDirPixTk->make<TH1F>("Trackertracks eta distribution","Tracker-tracks \\eta distribution" + hTitle[i], rint(2.*RangeEta/EtaStep), -RangeEta, RangeEta);
      H_trk_trk_eta[i]->SetXTitle("Track \\eta");
      H_trk_trk_eta[i]->SetYTitle("Entries (#)");

      hName = "Trackertracks theta distribution ver" + i_str.str();
      H_trk_trk_theta[i] = SubDirPixTk->make<TH1F>("Trackertracks theta distribution","Tracker-tracks \\theta distribution" + hTitle[i], 100, 0., 3.5);
      H_trk_trk_theta[i]->SetXTitle("Track \\theta (rad)");
      H_trk_trk_theta[i]->SetYTitle("Entries (#)");

      hName = "Trackertracks phi distribution ver" + i_str.str();
      H_trk_trk_phi[i] = SubDirPixTk->make<TH1F>("Trackertracks phi distribution","Tracker-tracks \\phi distribution" + hTitle[i], 100, -190., 190.);
      H_trk_trk_phi[i]->SetXTitle("Track \\phi (deg)");
      H_trk_trk_phi[i]->SetYTitle("Entries (#)");
    }

  H_trk_trk_doubleCounter = SubDirPixTk->make<TH1F>("Doublecounter trackertracks","Double-counter tracker-tracks", 10, 0.5, 10.5);
  H_trk_trk_doubleCounter->SetXTitle("Counts (#)");
  H_trk_trk_doubleCounter->SetYTitle("Entries (#)");

  H_trk_trk_doubleCounter_eta = SubDirPixTk->make<TH1F>("Doublecounter trackertracks eta distribution","Doublecounter tracker-tracks \\eta distribution", rint(2.*RangeEta/EtaStep), -RangeEta, RangeEta);
  H_trk_trk_doubleCounter_eta->SetXTitle("Track \\eta");
  H_trk_trk_doubleCounter_eta->SetYTitle("Entries (#)");

  H_trk_trk_doubleCounter_pt = SubDirPixTk->make<TH1F>("Doublecounter trackertracks pt distribution","Doublecounter tracker-tracks p_{t} distribution", rint(RangePt/PtStep), 0., RangePt);
  H_trk_trk_doubleCounter_pt->SetXTitle("Track p_{t} (GeV/c)");
  H_trk_trk_doubleCounter_pt->SetYTitle("Entries (#)");

  H_trk_trk_doubleCounter_eta_rel = SubDirPixTk->make<TH1F>("Doublecounter trackertracks eta relative distribution","Doublecounter tracker-tracks \\eta relative distribution", rint(2.*RangeEta/EtaStep), -RangeEta, RangeEta);
  H_trk_trk_doubleCounter_eta_rel->SetXTitle("Track \\eta");
  H_trk_trk_doubleCounter_eta_rel->SetYTitle("Entries (#)");

  H_trk_trk_doubleCounter_pt_rel = SubDirPixTk->make<TH1F>("Doublecounter trackertracks pt relative distribution","Doublecounter tracker-tracks p_{t} relative distribution", rint(RangePt/PtStep), 0., RangePt);
  H_trk_trk_doubleCounter_pt_rel->SetXTitle("Track p_{t} (GeV/c)");
  H_trk_trk_doubleCounter_pt_rel->SetYTitle("Entries (#)");

  H_trk_trk_eta_phi = SubDirPixTk->make<TH2F>("Trackertracks distribution in eta and phi","Tracker-track distribution in \\eta and \\phi", rint(2.*RangeEta/EtaStep), -RangeEta, RangeEta, rint(2.*RangePhi/PhiStep), -RangePhi, RangePhi);
  H_trk_trk_eta_phi->SetXTitle("Track \\eta");
  H_trk_trk_eta_phi->SetYTitle("Track \\phi (deg)");
  H_trk_trk_eta_phi->SetZTitle("Entries (#)");

  H_trk_trk_eta_pt = SubDirPixTk->make<TH2F>("Trackertracks distribution in eta and pt","Tracker-track distribution in \\eta and p_{t}", rint(2.*RangeEta/EtaStep), -RangeEta, RangeEta, rint(RangePt/PtStep), 0., RangePt);
  H_trk_trk_eta_pt->SetXTitle("Track \\eta");
  H_trk_trk_eta_pt->SetYTitle("Track p_{t} (GeV/c)");
  H_trk_trk_eta_pt->SetZTitle("Entries (#)");

  H_trk_trk_eta_phi_algo = SubDirPixTk->make<TH2F>("Trackertracks distribution in eta and phi algo","Tracker-track distribution in \\eta and \\phi for algo.eff", rint(2.*RangeEta/EtaStep), -RangeEta, RangeEta, rint(2.*RangePhi/PhiStep), -RangePhi, RangePhi);
  H_trk_trk_eta_phi_algo->SetXTitle("Track \\eta");
  H_trk_trk_eta_phi_algo->SetYTitle("Track \\phi (deg)");
  H_trk_trk_eta_phi_algo->SetZTitle("Entries (#)");

  H_trk_trk_eta_pt_algo = SubDirPixTk->make<TH2F>("Trackertracks distribution in eta and pt algo","Tracker-track distribution in \\eta and p_{t} for algo.eff", rint(2.*RangeEta/EtaStep), -RangeEta, RangeEta, rint(RangePt/PtStep), 0., RangePt);
  H_trk_trk_eta_pt_algo->SetXTitle("Track \\eta");
  H_trk_trk_eta_pt_algo->SetYTitle("Track p_{t} (GeV/c)");
  H_trk_trk_eta_pt_algo->SetZTitle("Entries (#)");

  H_trk_eff_pt_num = SubDirPixTk->make<TH1F>("Trackertracks efficiency vs pt numerator","Tracker-tracks efficiency vs p_{t} numerator", rint(RangePt/PtStep), 0., RangePt);
  H_trk_eff_pt_num->SetXTitle("Track p_{t} (GeV/c)");
  H_trk_eff_pt_num->SetYTitle("Entries (#)");

  H_trk_eff_pt_den = SubDirPixTk->make<TH1F>("Trackertracks efficiency vs pt denominator","Tracker-tracks efficiency vs p_{t} denominator", rint(RangePt/PtStep), 0., RangePt);
  H_trk_eff_pt_den->SetXTitle("Track p_{t} (GeV/c)");
  H_trk_eff_pt_den->SetYTitle("Entries (#)");

  H_trk_eff_pt = SubDirPixTk->make<TH1F>("Trackertracks efficiency vs pt","Tracker-tracks efficiency vs p_{t}", rint(RangePt/PtStep), 0., RangePt);
  H_trk_eff_pt->SetXTitle("Track p_{t} (GeV/c)");
  H_trk_eff_pt->SetYTitle("Entries (#)");

  H_trk_toteff_pt_num = SubDirPixTk->make<TH1F>("Total trackertracks efficiency vs pt numerator","Total tracker-tracks efficiency vs p_{t} numerator", rint(RangePt/PtStep), 0., RangePt);
  H_trk_toteff_pt_num->SetXTitle("Track p_{t} (GeV/c)");
  H_trk_toteff_pt_num->SetYTitle("Entries (#)");

  H_trk_toteff_pt_den = SubDirPixTk->make<TH1F>("Total trackertracks efficiency vs pt denominator","Total tracker-tracks efficiency vs p_{t} denominator", rint(RangePt/PtStep), 0., RangePt);
  H_trk_toteff_pt_den->SetXTitle("Track p_{t} (GeV/c)");
  H_trk_toteff_pt_den->SetYTitle("Entries (#)");

  H_trk_toteff_pt = SubDirPixTk->make<TH1F>("Total trackertracks efficiency vs pt","Total tracker-tracks efficiency vs p_{t}", rint(RangePt/PtStep), 0., RangePt);
  H_trk_toteff_pt->SetXTitle("Track p_{t} (GeV/c)");
  H_trk_toteff_pt->SetYTitle("Entries (#)");

  H_trk_eff_eta_num = SubDirPixTk->make<TH1F>("Trackertracks efficiency vs eta numerator","Tracker-tracks efficiency vs \\eta numerator", 100, -3.5, 3.5);
  H_trk_eff_eta_num->SetXTitle("Track \\eta");
  H_trk_eff_eta_num->SetYTitle("Entries (#)");

  H_trk_eff_eta_den = SubDirPixTk->make<TH1F>("Trackertracks efficiency vs eta denominator","Tracker-tracks efficiency vs \\eta denominator", 100, -3.5, 3.5);
  H_trk_eff_eta_den->SetXTitle("Track \\eta");
  H_trk_eff_eta_den->SetYTitle("Entries (#)");

  H_trk_eff_eta = SubDirPixTk->make<TH1F>("Trackertracks efficiency vs eta","Tracker-tracks efficiency vs \\eta", 100, -3.5, 3.5);
  H_trk_eff_eta->SetXTitle("Track \\eta");
  H_trk_eff_eta->SetYTitle("Entries (#)");

  H_trk_toteff_eta_num = SubDirPixTk->make<TH1F>("Total trackertracks efficiency vs eta numerator","Total tracker-tracks efficiency vs \\eta numerator", 100, -3.5, 3.5);
  H_trk_toteff_eta_num->SetXTitle("Track \\eta");
  H_trk_toteff_eta_num->SetYTitle("Entries (#)");

  H_trk_toteff_eta_den = SubDirPixTk->make<TH1F>("Total trackertracks efficiency vs eta denominator","Total tracker-tracks efficiency vs \\eta denominator", 100, -3.5, 3.5);
  H_trk_toteff_eta_den->SetXTitle("Track \\eta");
  H_trk_toteff_eta_den->SetYTitle("Entries (#)");

  H_trk_toteff_eta = SubDirPixTk->make<TH1F>("Total trackertracks efficiency vs eta","Total tracker-tracks efficiency vs \\eta", 100, -3.5, 3.5);
  H_trk_toteff_eta->SetXTitle("Track \\eta");
  H_trk_toteff_eta->SetYTitle("Entries (#)");

  hTitleEta[0] = "(-MaxEta<=\\eta<-thr)";
  hTitleEta[1] = "(-thr<=\\eta<=thr)";
  hTitleEta[2] = "(thr<\\eta<=MaxEta)";

  for (int i = 0; i < 3; i++)
    {
      stringstream i_str;
      i_str << i + 1;
      TString hName;

      hName = "Trackertracks efficiency vs phi numerator thr" + i_str.str();
      H_trk_eff_phi_num[i] = SubDirPixTk->make<TH1F>(hName,"Tracker-tracks efficiency vs \\phi numerator " + hTitleEta[i], rint(2.*RangePhi/PhiStep), -RangePhi, RangePhi);
      H_trk_eff_phi_num[i]->SetXTitle("Track \\phi (deg)");
      H_trk_eff_phi_num[i]->SetYTitle("Entries (#)");

      hName = "Trackertracks efficiency vs phi denominator thr" + i_str.str();
      H_trk_eff_phi_den[i] = SubDirPixTk->make<TH1F>(hName,"Tracker-tracks efficiency vs \\phi denominator " + hTitleEta[i], rint(2.*RangePhi/PhiStep), -RangePhi, RangePhi);
      H_trk_eff_phi_den[i]->SetXTitle("Track \\phi (deg)");
      H_trk_eff_phi_den[i]->SetYTitle("Entries (#)");

      hName = "Trackertracks efficiency vs phi thr" + i_str.str();
      H_trk_eff_phi[i] = SubDirPixTk->make<TH1F>(hName,"Tracker-tracks efficiency vs \\phi " + hTitleEta[i], rint(2.*RangePhi/PhiStep), -RangePhi, RangePhi);
      H_trk_eff_phi[i]->SetXTitle("Track \\phi (deg)");
      H_trk_eff_phi[i]->SetYTitle("Entries (#)");

      hName = "Total trackertracks efficiency vs phi numerator thr" + i_str.str();
      H_trk_toteff_phi_num[i] = SubDirPixTk->make<TH1F>(hName,"Total tracker-tracks efficiency vs \\phi numerator " + hTitleEta[i], rint(2.*RangePhi/PhiStep), -RangePhi, RangePhi);
      H_trk_toteff_phi_num[i]->SetXTitle("Track \\phi (deg)");
      H_trk_toteff_phi_num[i]->SetYTitle("Entries (#)");

      hName = "Total trackertracks efficiency vs phi denominator thr" + i_str.str();
      H_trk_toteff_phi_den[i] = SubDirPixTk->make<TH1F>(hName,"Total tracker-tracks efficiency vs \\phi denominator " + hTitleEta[i], rint(2.*RangePhi/PhiStep), -RangePhi, RangePhi);
      H_trk_toteff_phi_den[i]->SetXTitle("Track \\phi (deg)");
      H_trk_toteff_phi_den[i]->SetYTitle("Entries (#)");

      hName = "Total trackertracks efficiency vs phi thr" + i_str.str();
      H_trk_toteff_phi[i] = SubDirPixTk->make<TH1F>(hName,"Total tracker-tracks efficiency vs \\phi " + hTitleEta[i], rint(2.*RangePhi/PhiStep), -RangePhi, RangePhi);
      H_trk_toteff_phi[i]->SetXTitle("Track \\phi (deg)");
      H_trk_toteff_phi[i]->SetYTitle("Entries (#)");
    }

  // ######################
  // # Common track plots #
  // ######################
  hTitle[0] = "(thr0<p_{t}<=thr1 GeV/c)";
  hTitle[1] = "(thr1<p_{t}<=thr2 GeV/c)";
  hTitle[2] = "(thr2<p_{t}<=thr3 GeV/c)";
  hTitle[3] = "(thr3<p_{t}<=thr4 GeV/c)";
  hTitle[4] = "(thr4<p_{t}<=RangePt GeV/c)";

  hTitleEta[0] = "(\\eta<=thr)";
  hTitleEta[1] = "(thr<\\eta<=MaxEta)";

  for (int i = 0; i < 5; i++) 
    {
      stringstream i_str;
      i_str << i + 1;
      TString hName;

      if (i < 2)
	{
	  // ################
	  // # PIXEL-TRACKS #
	  // ################
	  hName = "Pixeltracks vs chi2DoF and pt cuts thr" + i_str.str(); 
	  H_pix_trk_normchi2_pt_eta[i] = SubDirPixTk->make<TH2F>(hName,"Pixel-tracks vs. \\chi^{2}/DoF and p_{t} cuts " + hTitleEta[i], rint(RangePt/PtStep), 0., RangePt, rint(RangeChi2/Chi2Step), 0., RangeChi2);
	  H_pix_trk_normchi2_pt_eta[i]->SetXTitle("Track p_{t} (=>cut) (GeV/c)");
	  H_pix_trk_normchi2_pt_eta[i]->SetYTitle("\\chi^{2}/DoF (<=cut)");
	  H_pix_trk_normchi2_pt_eta[i]->SetZTitle("Entries (#)");

	  // ##################
	  // # TRACKER-TRACKS #
	  // ##################
	  hName = "Trackertracks vs chi2DoF and pt cuts thr" + i_str.str();
	  H_trk_trk_normchi2_pt_eta[i] = SubDirPixTk->make<TH2F>(hName,"Tracker-tracks vs. \\chi^{2}/DoF and p_{t} cuts " + hTitleEta[0], rint(RangePt/PtStep), 0., RangePt, rint(RangeChi2/Chi2Step), 0., RangeChi2);
	  H_trk_trk_normchi2_pt_eta[i]->SetXTitle("Track p_{t} (=>cut) (GeV/c)");
	  H_trk_trk_normchi2_pt_eta[i]->SetYTitle("\\chi^{2}/DoF (<=cut)");
	  H_trk_trk_normchi2_pt_eta[i]->SetZTitle("Entries (#)");

	  // ######################
	  // # COMMON TRACK PLOTS #
	  // ######################
	  hName = "d0 resolution vs pt thr" + i_str.str();  
	  H_d0res_pt_eta[i] = SubDirPixTkCom->make<TH1F>(hName,"d_{0} resolution vs. p_{t} " + hTitleEta[i], rint(RangePt/PtStep), 0., RangePt);
	  H_d0res_pt_eta[i]->SetXTitle("Track p_{t} (GeV/c)");
	  H_d0res_pt_eta[i]->SetYTitle("\\sigma(d_{0}) (cm)");

	  hName = "d0 pulls vs pt thr" + i_str.str();  
	  H_d0pull_pt_eta[i] = SubDirPixTkCom->make<TH1F>(hName,"d_{0} pulls vs. p_{t} " + hTitleEta[i], rint(RangePt/PtStep), 0., RangePt);
	  H_d0pull_pt_eta[i]->SetXTitle("Track p_{t} (GeV/c)");
	  H_d0pull_pt_eta[i]->SetYTitle("Track d_{0} pulls");

	  hName = "dz resolution vs pt thr" + i_str.str();  
	  H_dzres_pt_eta[i] = SubDirPixTkCom->make<TH1F>(hName,"d_{z} resolution vs. p_{t} " + hTitleEta[i], rint(RangePt/PtStep), 0., RangePt);
	  H_dzres_pt_eta[i]->SetXTitle("Track p_{t} (GeV/c)");
	  H_dzres_pt_eta[i]->SetYTitle("\\sigma(d_{z}) (cm)");

	  hName = "dz pulls vs pt thr" + i_str.str();  
	  H_dzpull_pt_eta[i] = SubDirPixTkCom->make<TH1F>(hName,"d_{z} pulls vs. p_{t} " + hTitleEta[i], rint(RangePt/PtStep), 0., RangePt);
	  H_dzpull_pt_eta[i]->SetXTitle("Track p_{t} (GeV/c)");
	  H_dzpull_pt_eta[i]->SetYTitle("Track d_{z} pulls");

	  hName = "pt resolution vs pt thr" + i_str.str();
	  H_ptres_pt_eta[i] = SubDirPixTkCom->make<TH1F>(hName,"p_{t} resolution vs. p_{t} " + hTitleEta[i], rint(RangePt/PtStep), 0., RangePt);
	  H_ptres_pt_eta[i]->SetXTitle("Track p_{t} (GeV/c)");
	  H_ptres_pt_eta[i]->SetYTitle("\\sigma(p_{t})/p_{t}");

	  hName = "pt pulls vs pt thr" + i_str.str();
	  H_ptpull_pt_eta[i] = SubDirPixTkCom->make<TH1F>(hName,"p_{t} pulls vs. p_{t} " + hTitleEta[i], rint(RangePt/PtStep), 0., RangePt);
	  H_ptpull_pt_eta[i]->SetXTitle("Track p_{t} (GeV/c)");
	  H_ptpull_pt_eta[i]->SetYTitle("Track p_{t} pulls");

	  hName = "SkewPixTrkpt vs pt thr" + i_str.str();
	  H_skew_pt_eta[i] = SubDirPixTkCom->make<TH1F>(hName,"Skew(Pix-Trk)_{pt} vs. p_{t} " + hTitleEta[i], rint(RangePt/PtStep), 0., RangePt);
	  H_skew_pt_eta[i]->SetXTitle("Track p_{t} (GeV/c)");
	  H_skew_pt_eta[i]->SetYTitle("Skew(Pix-Trk)_{pt}");

	  // ##########
	  // # PURITY #
	  // ##########
	  hName = "Matched pixeltracks vs chi2DoF and pt cuts pur thr" + i_str.str();
	  H_trked_normchi2_pt_pur_eta[i] = SubDirPixTkCom->make<TH2F>(hName,"Matched pixel-tracks vs. \\chi^{2}/DoF and p_{t} cuts for purity " + hTitleEta[0], rint(RangePt/PtStep), 0., RangePt, rint(RangeChi2/Chi2Step), 0., RangeChi2);
	  H_trked_normchi2_pt_pur_eta[i]->SetXTitle("Track p_{t} (=>cut) (GeV/c)");
	  H_trked_normchi2_pt_pur_eta[i]->SetYTitle("\\chi^{2}/DoF (<=cut)");
	  H_trked_normchi2_pt_pur_eta[i]->SetZTitle("Entries (#)");

	  hName = "Purity pixeltracks vs pt thr" + i_str.str();    
	  H_purity_trked_pt_eta[i] = SubDirPixTkCom->make<TH1F>(hName,"Purity pixel-tracks vs p_{t} " + hTitleEta[i], rint(RangePt/PtStep), 0., RangePt);
	  H_purity_trked_pt_eta[i]->SetXTitle("Track p_{t} (GeV/c)");
	  H_purity_trked_pt_eta[i]->SetYTitle("Pixel-tracks purity");

	  hName = "Purity pixeltracks vs chi2DoF and pt cuts thr" + i_str.str(); 
	  H_purity_trked_normchi2_pt_eta[i] = SubDirPixTkCom->make<TH2F>(hName,"Purity pixel-tracks vs. \\chi^{2}/DoF and p_{t} cuts " +  hTitleEta[0], rint(RangePt/PtStep), 0., RangePt, rint(RangeChi2/Chi2Step), 0., RangeChi2);
	  H_purity_trked_normchi2_pt_eta[i]->SetXTitle("Track p_{t} (=>cut) (GeV/c)");
	  H_purity_trked_normchi2_pt_eta[i]->SetYTitle("\\chi^{2}/DoF (<=cut)");
	  H_purity_trked_normchi2_pt_eta[i]->SetZTitle("Pixel-tracks purity");

	  hName = "Matched pixeltracks vs chi2DoF and pt cuts eff thr" + i_str.str();     
	  H_trked_normchi2_pt_eff_eta[i] = SubDirPixTkCom->make<TH2F>(hName,"Matched pixel-tracks vs. \\chi^{2}/DoF and p_{t} cuts for efficiency " + hTitleEta[i], rint(RangePt/PtStep), 0., RangePt, rint(RangeChi2/Chi2Step), 0., RangeChi2);
	  H_trked_normchi2_pt_eff_eta[i]->SetXTitle("Track p_{t} (=>cut) (GeV/c)");
	  H_trked_normchi2_pt_eff_eta[i]->SetYTitle("\\chi^{2}/DoF (<=cut)");
	  H_trked_normchi2_pt_eff_eta[i]->SetZTitle("Entries (#)");

	  // ##############
	  // # EFFICIENCY #
	  // ##############
	  hName = "Efficiency pixeltracks vs pt thr" + i_str.str();        
	  H_efficiency_trked_pt_eta[i] = SubDirPixTkCom->make<TH1F>(hName,"Efficiency pixel-tracks vs p_{t} " + hTitleEta[i], rint(RangePt/PtStep), 0., RangePt);
	  H_efficiency_trked_pt_eta[i]->SetXTitle("Track p_{t} (GeV/c)");
	  H_efficiency_trked_pt_eta[i]->SetYTitle("Pixel-tracks efficiency");

	  hName = "Algoefficiency pixeltracks vs pt thr" + i_str.str();         
	  H_efficiency_trked_pt_eta_algo[i] = SubDirPixTkCom->make<TH1F>(hName,"Algo. efficiency pixel-tracks vs p_{t} " + hTitleEta[i], rint(RangePt/PtStep), 0., RangePt);
	  H_efficiency_trked_pt_eta_algo[i]->SetXTitle("Track p_{t} (GeV/c)");
	  H_efficiency_trked_pt_eta_algo[i]->SetYTitle("Pixel-tracks efficiency");

	  hName = "Efficiency pixeltracks vs chi2DoF and pt cuts thr" + i_str.str();     
	  H_efficiency_trked_normchi2_pt_eta[i] = SubDirPixTkCom->make<TH2F>(hName,"Efficiency pixel-tracks vs. \\chi^{2}/DoF and p_{t} cuts " + hTitleEta[i], rint(RangePt/PtStep), 0., RangePt, rint(RangeChi2/Chi2Step), 0., RangeChi2);
	  H_efficiency_trked_normchi2_pt_eta[i]->SetXTitle("Track p_{t} (=>cut) (GeV/c)");
	  H_efficiency_trked_normchi2_pt_eta[i]->SetYTitle("\\chi^{2}/DoF (<=cut)");
	  H_efficiency_trked_normchi2_pt_eta[i]->SetZTitle("Pixel-tracks efficiency");      
	} 

      hName = "d0 resolution vs eta thr" + i_str.str();
      H_d0res_eta_pt[i] = SubDirPixTkCom->make<TH1F>(hName,"d_{0} resolution vs. \\eta " + hTitle[i], rint(RangeEta/EtaStep), 0., RangeEta);
      H_d0res_eta_pt[i]->SetXTitle("Track \\eta");
      H_d0res_eta_pt[i]->SetYTitle("\\sigma(d_{0}) (cm)");

      hName = "d0 pulls vs eta thr" + i_str.str();  
      H_d0pull_eta_pt[i] = SubDirPixTkCom->make<TH1F>(hName,"d_{0} pulls vs. \\eta " + hTitle[i], rint(RangeEta/EtaStep), 0., RangeEta);
      H_d0pull_eta_pt[i]->SetXTitle("Track \\eta");
      H_d0pull_eta_pt[i]->SetYTitle("Track d_{0} pulls");

      hName = "dz resolution vs eta thr" + i_str.str();
      H_dzres_eta_pt[i] = SubDirPixTkCom->make<TH1F>(hName,"d_{z} resolution vs. \\eta " + hTitle[i], rint(RangeEta/EtaStep), 0., RangeEta);
      H_dzres_eta_pt[i]->SetXTitle("Track \\eta");
      H_dzres_eta_pt[i]->SetYTitle("\\sigma(d_{z}) (cm)");

      hName = "dz pulls vs eta thr" + i_str.str();  
      H_dzpull_eta_pt[i] = SubDirPixTkCom->make<TH1F>(hName,"d_{z} pulls vs. \\eta " + hTitle[i], rint(RangeEta/EtaStep), 0., RangeEta);
      H_dzpull_eta_pt[i]->SetXTitle("Track \\eta");
      H_dzpull_eta_pt[i]->SetYTitle("Track d_{z} pulls");

      hName = "pt resolution vs eta thr" + i_str.str();
      H_ptres_eta_pt[i] = SubDirPixTkCom->make<TH1F>(hName,"p_{t} resolution vs. \\eta " + hTitle[i], rint(RangeEta/EtaStep), 0., RangeEta);
      H_ptres_eta_pt[i]->SetXTitle("Track \\eta");
      H_ptres_eta_pt[i]->SetYTitle("\\sigma(p_{t})/p_{t}");

      hName = "pt pulls vs eta thr" + i_str.str();
      H_ptpull_eta_pt[i] = SubDirPixTkCom->make<TH1F>(hName,"p_{t} pulls vs. \\eta " + hTitle[i], rint(RangeEta/EtaStep), 0., RangeEta);
      H_ptpull_eta_pt[i]->SetXTitle("Track \\eta");
      H_ptpull_eta_pt[i]->SetYTitle("Track p_{t} pulls");
    }

  H_purity_trked_eta = SubDirPixTkCom->make<TH1F>("Purity pixeltracks vs eta","Purity pixel-tracks vs \\eta", rint(2.*RangeEta/EtaStep), -RangeEta, RangeEta);
  H_purity_trked_eta->SetXTitle("Track \\eta");
  H_purity_trked_eta->SetYTitle("Pixel-tracks purity");

  H_efficiency_trked_eta = SubDirPixTkCom->make<TH1F>("Efficiency pixeltracks vs eta","Efficiency pixel-tracks vs \\eta", rint(2.*RangeEta/EtaStep), -RangeEta, RangeEta);
  H_efficiency_trked_eta->SetXTitle("Track \\eta");
  H_efficiency_trked_eta->SetYTitle("Pixel-tracks efficiency");

  H_efficiency_trked_eta_algo = SubDirPixTkCom->make<TH1F>("Algoefficiency pixeltracks vs eta","Algo. efficiency pixel-tracks vs \\eta", rint(2.*RangeEta/EtaStep), -RangeEta, RangeEta);
  H_efficiency_trked_eta_algo->SetXTitle("Track \\eta");
  H_efficiency_trked_eta_algo->SetYTitle("Pixel-tracks algo-efficiency");

  H_efficiency_trked_pt_whole_eta = SubDirPixTkCom->make<TH1F>("Efficiency pixeltracks vs pt","Efficiency pixel-tracks vs p_{t}", rint(RangePt/PtStep), 0., RangePt);
  H_efficiency_trked_pt_whole_eta->SetXTitle("Track p_{t} (GeV/c)");
  H_efficiency_trked_pt_whole_eta->SetYTitle("Pixel-tracks efficiency");

  H_efficiency_trked_pt_whole_eta_algo = SubDirPixTkCom->make<TH1F>("Algoefficiency pixeltracks vs pt","Algo. efficiency pixel-tracks vs p_{t}", rint(RangePt/PtStep), 0., RangePt);
  H_efficiency_trked_pt_whole_eta_algo->SetXTitle("Track p_{t} (GeV/c)");
  H_efficiency_trked_pt_whole_eta_algo->SetYTitle("Pixel-tracks algo-efficiency");

  H_purity_trked_pt_whole_eta = SubDirPixTkCom->make<TH1F>("Purity pixeltracks vs pt","Purity pixel-tracks vs p_{t}", rint(RangePt/PtStep), 0., RangePt);
  H_purity_trked_pt_whole_eta->SetXTitle("Track p_{t} (GeV/c)");
  H_purity_trked_pt_whole_eta->SetYTitle("Pixel-tracks purity");

  hTitleEta[0] = "(-MaxEta<=\\eta<-thr)";
  hTitleEta[1] = "(-thr<=\\eta<=thr)";
  hTitleEta[2] = "(thr<\\eta<=MaxEta)";
 
  for (int i = 0; i < 3; i++) 
    {
      stringstream i_str;
      i_str << i + 1;
      TString hName;

      hName = "Efficiency pixeltracks vs phi thr" + i_str.str();
      H_efficiency_trked_phi[i] = SubDirPixTkCom->make<TH1F>(hName,"Efficiency pixel-tracks vs \\phi " + hTitleEta[i], rint(2.*RangePhi/PhiStep), -RangePhi, RangePhi);
      H_efficiency_trked_phi[i]->SetXTitle("Track \\phi (deg)");
      H_efficiency_trked_phi[i]->SetYTitle("Pixel-tracks efficiency");

      hName = "Algoefficiency pixeltracks vs phi thr" + i_str.str(); 
      H_efficiency_trked_phi_algo[i] = SubDirPixTkCom->make<TH1F>(hName,"Algo. efficiency pixel-tracks vs \\phi " + hTitleEta[i], rint(2.*RangePhi/PhiStep), -RangePhi, RangePhi);
      H_efficiency_trked_phi_algo[i]->SetXTitle("Track \\phi (deg)");
      H_efficiency_trked_phi_algo[i]->SetYTitle("Pixel-tracks algo-efficiency");

      hName = "Purity pixeltracks vs phi thr" + i_str.str(); 
      H_purity_trked_phi[i] = SubDirPixTkCom->make<TH1F>(hName,"Purity pixel-tracks vs \\phi " + hTitleEta[i], rint(2.*RangePhi/PhiStep), -RangePhi, RangePhi);
      H_purity_trked_phi[i]->SetXTitle("Track \\phi (deg)");
      H_purity_trked_phi[i]->SetYTitle("Pixel-tracks purity");
    }

  H_tk_pur_eff = SubDirPixTkCom->make<TH1F>("Purity Efficiency Tracks","Purity Efficiency Tracks", 3, 0., 3.);
  H_tk_pur_eff->SetXTitle("");
  H_tk_pur_eff->SetYTitle("Entries (#)");
  H_tk_pur_eff->GetXaxis()->SetBinLabel(1,"Purity");
  H_tk_pur_eff->GetXaxis()->SetBinLabel(2,"Efficiency");
  H_tk_pur_eff->GetXaxis()->SetBinLabel(3,"Algo-Efficiency");

  H_trked_eta_phi_pur = SubDirPixTkCom->make<TH2F>("Matched pixeltracks distribution in eta and phi pur","Matched pixel-tracks distribution in \\eta and \\phi for purity", rint(2.*RangeEta/EtaStep), -RangeEta, RangeEta, rint(2.*RangePhi/PhiStep), -RangePhi, RangePhi);
  H_trked_eta_phi_pur->SetXTitle("Track \\eta");
  H_trked_eta_phi_pur->SetYTitle("Track \\phi (deg)");
  H_trked_eta_phi_pur->SetZTitle("Entries (#)");

  H_trked_eta_phi_eff = SubDirPixTkCom->make<TH2F>("Matched pixeltracks distribution in eta and phi eff","Matched pixel-tracks distribution in \\eta and \\phi for efficiency", rint(2.*RangeEta/EtaStep), -RangeEta, RangeEta, rint(2.*RangePhi/PhiStep), -RangePhi, RangePhi);
  H_trked_eta_phi_eff->SetXTitle("Track \\eta");
  H_trked_eta_phi_eff->SetYTitle("Track \\phi (deg)");
  H_trked_eta_phi_eff->SetZTitle("Entries (#)");

  H_trked_eta_pt_pur = SubDirPixTkCom->make<TH2F>("Matched pixeltracks distribution in eta and pt pur","Matched pixel-tracks distribution in \\eta and p_{t} for purity", rint(2.*RangeEta/EtaStep), -RangeEta, RangeEta, rint(RangePt/PtStep), 0., RangePt);
  H_trked_eta_pt_pur->SetXTitle("Track \\eta");
  H_trked_eta_pt_pur->SetYTitle("Track p_{t} (GeV/c)");
  H_trked_eta_pt_pur->SetZTitle("Entries (#)");

  H_trked_eta_pt_eff = SubDirPixTkCom->make<TH2F>("Matched pixeltracks distribution in eta and pt eff","Matched pixel-tracks distribution in \\eta and p_{t} for efficiency", rint(2.*RangeEta/EtaStep), -RangeEta, RangeEta, rint(RangePt/PtStep), 0., RangePt);
  H_trked_eta_pt_eff->SetXTitle("Track \\eta");
  H_trked_eta_pt_eff->SetYTitle("Track p_{t} (GeV/c)");
  H_trked_eta_pt_eff->SetZTitle("Entries (#)");

  // ##################
  // # Pixel Vertices #
  // ##################
  H_pix_vx_z = SubDirPixVx->make<TH1F>("Pixelvertices distribution along Z","Pixel-vertices distribution along Z", 80, -20., 20.);
  H_pix_vx_z->SetXTitle("Vertex Z position (cm)");
  H_pix_vx_z->SetYTitle("Entries (#)");

  H_pix_vx_normchi2 = SubDirPixVx->make<TH1F>("Pixelvertices chi2DoF distribution","Pixel-vertices \\chi^{2}/DoF distribution", 50, 0., 5.);
  H_pix_vx_normchi2->SetXTitle("Vertex \\chi^{2}/DoF");
  H_pix_vx_normchi2->SetYTitle("Entries (#)");

  H_pix_vx_dof = SubDirPixVx->make<TH1F>("Pixelvertices DoF distribution","Pixel-vertices DoF distribution", 100, 0., 30.);
  H_pix_vx_dof->SetXTitle("Vertex DoF (#)");
  H_pix_vx_dof->SetYTitle("Entries (#)");

  H_pix_vx_xerr = SubDirPixVx->make<TH1F>("Pixelvertices error on X position","Pixel-vertices error on X position", 100, 0., 0.1);
  H_pix_vx_xerr->SetXTitle("Vertex X position error (cm)");
  H_pix_vx_xerr->SetYTitle("Entries (#)");

  H_pix_vx_yerr = SubDirPixVx->make<TH1F>("Pixelvertices error on Y position","Pixel-vertices error on Y position", 100, 0., 0.1);
  H_pix_vx_yerr->SetXTitle("Vertex Y position error (cm)");
  H_pix_vx_yerr->SetYTitle("Entries (#)");

  H_pix_vx_zerr = SubDirPixVx->make<TH1F>("Pixelvertices error on Z position","Pixel-vertices error on Z position", 100, 0., 0.1);
  H_pix_vx_zerr->SetXTitle("Vertex Z position error (cm)");
  H_pix_vx_zerr->SetYTitle("Entries (#)");

  H_pix_vx_xy = SubDirPixVx->make<TH2F>("Pixelvertices distribution in the transversal plane","Pixel-vertices distribution in the transversal plane",  800, -4000., 4000., 800, -4000., 4000.);
  H_pix_vx_xy->SetXTitle("Vertex X position (\\mum)");
  H_pix_vx_xy->SetYTitle("Vertex Y position (\\mum)");
  H_pix_vx_xy->SetZTitle("Entries (#)");

  H_ntrk_pix_vx = SubDirPixVx->make<TH1F>("Number of tracks per pixelvertices","Number of tracks per pixel-vertices", 50, -0.5, 49.5);
  H_ntrk_pix_vx->SetXTitle("N. tracks");
  H_ntrk_pix_vx->SetYTitle("Entries (#)");

  // ####################
  // # Tracker Vertices #
  // ####################
  H_trk_vx_z = SubDirPixVx->make<TH1F>("Trackervertices distribution along Z","Tracker-vertices distribution along Z", 80, -20., 20.);
  H_trk_vx_z->SetXTitle("Vertex Z position (cm)");
  H_trk_vx_z->SetYTitle("Entries (#)");

  H_trk_vx_normchi2 = SubDirPixVx->make<TH1F>("Trackervertices chi2DoF distribution","Tracker-vertices \\chi^{2}/DoF distribution", 50, 0., 5.);
  H_trk_vx_normchi2->SetXTitle("Vertex \\chi^{2}/DoF");
  H_trk_vx_normchi2->SetYTitle("Entries (#)");

  H_trk_vx_dof = SubDirPixVx->make<TH1F>("Trackervertices DoF distribution","Tracker-vertices DoF distribution", 100, 0., 100.);
  H_trk_vx_dof->SetXTitle("Vertex DoF (#)");
  H_trk_vx_dof->SetYTitle("Entries (#)");

  H_trk_vx_xerr = SubDirPixVx->make<TH1F>("Trackervertices error on X position","Tracker-vertices error on X position", 100, 0., 0.05);
  H_trk_vx_xerr->SetXTitle("Vertex X position error (cm)");
  H_trk_vx_xerr->SetYTitle("Entries (#)");

  H_trk_vx_yerr = SubDirPixVx->make<TH1F>("Trackervertices error on Y position","Tracker-vertices error on Y position", 100, 0., 0.05);
  H_trk_vx_yerr->SetXTitle("Vertex Y position error (cm)");
  H_trk_vx_yerr->SetYTitle("Entries (#)");

  H_trk_vx_zerr = SubDirPixVx->make<TH1F>("Trackervertices error on Z position","Tracker-vertices error on Z position", 100, 0., 0.05);
  H_trk_vx_zerr->SetXTitle("Vertex Z position error (cm)");
  H_trk_vx_zerr->SetYTitle("Entries (#)");

  H_trk_vx_xy = SubDirPixVx->make<TH2F>("Trackervertices distribution in the transversal plane","Tracker-vertices distribution in the transversal plane", 800, -4000., 4000., 800, -4000., 4000.);
  H_trk_vx_xy->SetXTitle("Vertex X position (\\mum)");
  H_trk_vx_xy->SetYTitle("Vertex Y position (\\mum)");
  H_trk_vx_xy->SetZTitle("Entries (#)");

  H_ntrk_trk_vx = SubDirPixVx->make<TH1F>("Number of tracks per trackervertices","Number of tracks per tracker-vertices", 100, -0.5, 99.5);
  H_ntrk_trk_vx->SetXTitle("N. tracks");
  H_ntrk_trk_vx->SetYTitle("Entries (#)");

  H_trk_vx_z_nomatch = SubDirPixVxNoM->make<TH1F>("Nomatched trackervertices vs Z","Nomatched tracker-vertices vs. Z", 80, -20., 20.);
  H_trk_vx_z_nomatch->SetXTitle("Vertex Z position (cm)");
  H_trk_vx_z_nomatch->SetYTitle("Entries (#)");

  H_trk_vx_xy_nomatch = SubDirPixVxNoM->make<TH2F>("Nomatched trackervertices vs XY","Nomatched tracker-vertices vs. XY", 800, -4000., 4000., 800, -4000., 4000.);
  H_trk_vx_xy_nomatch->SetXTitle("Vertex X position (\\mum)");
  H_trk_vx_xy_nomatch->SetYTitle("Vertex Y position (\\mum)");
  H_trk_vx_xy_nomatch->SetZTitle("Entries (#)");

  H_trk_vx_nomatch_trk_pt = SubDirPixVxNoM->make<TH1F>("pt distribution of tracks of nomatched trackervertices","p_{t} distribution of tracks of nomatched tracker-vertices", 100, 0., 10.);
  H_trk_vx_nomatch_trk_pt->SetXTitle("Track p_{t} (GeV/c)");
  H_trk_vx_nomatch_trk_pt->SetYTitle("Entries (#)");

  H_trk_vx_nomatch_trk_eta = SubDirPixVxNoM->make<TH1F>("eta distribution of tracks of nomatched trackervertices","\\eta distribution of tracks of nomatched tracker-vertices", 100, -3.5, 3.5);
  H_trk_vx_nomatch_trk_eta->SetXTitle("Track \\eta");
  H_trk_vx_nomatch_trk_eta->SetYTitle("Entries (#)");

  H_trk_vx_nomatch_normchi2 = SubDirPixVxNoM->make<TH1F>("Chi2DoF distribution of nomatched trackervertices","\\chi^{2}/DoF distribution of nomatched tracker-vertices", 100, 0., 30.);
  H_trk_vx_nomatch_normchi2->SetXTitle("Vertex \\chi^{2}/DoF");
  H_trk_vx_nomatch_normchi2->SetYTitle("Entries (#)");

  H_trk_vx_nomatch_dof = SubDirPixVxNoM->make<TH1F>("Nomatched trackervertices DoF distribution","Nomatched tracker-vertices DoF distribution", 100, 0., 100.);
  H_trk_vx_nomatch_dof->SetXTitle("Vertex DoF (#)");
  H_trk_vx_nomatch_dof->SetYTitle("Entries (#)");

  // #######################
  // # Common vertex plots #
  // #######################
  H_vx_x_diff_vs_trk = SubDirPixVxCom->make<TH1F>("Vertex position difference along X vs Ntrks","Vertex position difference along X vs. Ntrks", 100, -0.5, 99.5);
  H_vx_x_diff_vs_trk->SetXTitle("Number of Tracks");
  H_vx_x_diff_vs_trk->SetYTitle("\\sigma(Trk_{X}-Pix_{X}) (\\mum)");

  H_vx_y_diff_vs_trk = SubDirPixVxCom->make<TH1F>("Vertex position difference along Y vs Ntrks","Vertex position difference along Y vs. Ntrks", 100, -0.5, 99.5);
  H_vx_y_diff_vs_trk->SetXTitle("Number of Tracks");
  H_vx_y_diff_vs_trk->SetYTitle("\\sigma(Trk_{Y}-Pix_{Y}) (\\mum)");

  H_vx_z_diff_vs_trk = SubDirPixVxCom->make<TH1F>("Vertex position difference along Z vs Ntrks","Vertex position difference along Z vs. Ntrks", 100, -0.5, 99.5);
  H_vx_z_diff_vs_trk->SetXTitle("Number of Tracks");
  H_vx_z_diff_vs_trk->SetYTitle("\\sigma(Trk_{Z}-Pix_{Z}) (\\mum)");

  H_vx_x_diff_pull_vs_trk = SubDirPixVxCom->make<TH1F>("Vertex position pull along X vs Ntrks","Vertex position pull along X vs. Ntrks", 100, -0.5, 99.5);
  H_vx_x_diff_pull_vs_trk->SetXTitle("Number of Tracks");
  H_vx_x_diff_pull_vs_trk->SetYTitle("\\sigma((Trk_{X}-Pix_{X}) / (\\sigma^{2}(Pix_{X}) + \\sigma^{2}(Trk_{X}))^{1/2})");

  H_vx_y_diff_pull_vs_trk = SubDirPixVxCom->make<TH1F>("Vertex position pull along Y vs Ntrks","Vertex position pull along Y vs. Ntrks", 100, -0.5, 99.5);
  H_vx_y_diff_pull_vs_trk->SetXTitle("Number of Tracks");
  H_vx_y_diff_pull_vs_trk->SetYTitle("\\sigma((Trk_{Y}-Pix_{Y}) / (\\sigma^{2}(Pix_{X}) + \\sigma^{2}(Trk_{Y}))^{1/2})");

  H_vx_z_diff_pull_vs_trk = SubDirPixVxCom->make<TH1F>("Vertex position pull along Z vs Ntrks","Vertex position pull along Z vs. Ntrks", 100, -0.5, 99.5);
  H_vx_z_diff_pull_vs_trk->SetXTitle("Number of Tracks");
  H_vx_z_diff_pull_vs_trk->SetYTitle("\\sigma((Trk_{Z}-Pix_{Z}) / (\\sigma^{2}(Pix_{Z}) + \\sigma^{2}(Trk_{Z}))^{1/2})");

  H_vx_x_diff_orig_vs_trk = SubDirPixVxCom->make<TH1F>("Original vertex position difference along X vs Ntrks","Original vertex position difference along X vs. Ntrks", 100, -0.5, 99.5);
  H_vx_x_diff_orig_vs_trk->SetXTitle("Number of Tracks");
  H_vx_x_diff_orig_vs_trk->SetYTitle("\\sigma(Trk_{X}-TrkOrig_{X}) (\\mum)");

  H_vx_y_diff_orig_vs_trk = SubDirPixVxCom->make<TH1F>("Original vertex position difference along Y vs Ntrks","Original vertex position difference along Y vs. Ntrks", 100, -0.5, 99.5);
  H_vx_y_diff_orig_vs_trk->SetXTitle("Number of Tracks");
  H_vx_y_diff_orig_vs_trk->SetYTitle("\\sigma(Trk_{Y}-TrkOrig_{Y}) (\\mum)");

  H_vx_z_diff_orig_vs_trk = SubDirPixVxCom->make<TH1F>("Original vertex position difference along Z vs Ntrks","Original vertex position difference along Z vs. Ntrks", 100, -0.5, 99.5);
  H_vx_z_diff_orig_vs_trk->SetXTitle("Number of Tracks");
  H_vx_z_diff_orig_vs_trk->SetYTitle("\\sigma(Trk_{Z}-TrkOrig_{Z}) (\\mum)");

  H_vx_x_diff_noerr_vs_trk = SubDirPixVxCom->make<TH1F>("Vertex resolution along X vs Ntrks","Vertex resolution along X vs. Ntrks", 100, -0.5, 99.5);
  H_vx_x_diff_noerr_vs_trk->SetXTitle("Number of Tracks");
  H_vx_x_diff_noerr_vs_trk->SetYTitle("#sqrt{\\sigma^{2}(Pix_{X}-Trk_{X}) - \\mu(\\sigma^{2}_{TrkX})} (\\mum)");

  H_vx_y_diff_noerr_vs_trk = SubDirPixVxCom->make<TH1F>("Vertex resolution along Y vs Ntrks","Vertex resolution along Y vs. Ntrks", 100, -0.5, 99.5);
  H_vx_y_diff_noerr_vs_trk->SetXTitle("Number of Tracks");
  H_vx_y_diff_noerr_vs_trk->SetYTitle("#sqrt{\\sigma^{2}(Pix_{Y}-Trk_{Y}) - \\mu(\\sigma^{2}_{TrkY})} (\\mum)");

  H_vx_z_diff_noerr_vs_trk = SubDirPixVxCom->make<TH1F>("Vertex resolution along Z vs Ntrks","Vertex resolution along Z vs. Ntrks", 100, -0.5, 99.5);
  H_vx_z_diff_noerr_vs_trk->SetXTitle("Number of Tracks");
  H_vx_z_diff_noerr_vs_trk->SetYTitle("#sqrt{\\sigma^{2}(Pix_{Z}-Trk_{Z}) - \\mu(\\sigma^{2}_{TrkZ})} (\\mum)");

  H_vx_counters = SubDirPixVxCom->make<TH1F>("Counters Match Purity Efficiency for Vertices","Counters Match Purity Efficiency for Vertices", 4, 0., 4.);
  H_vx_counters->SetXTitle("");
  H_vx_counters->SetYTitle("Entries (#)");
  H_vx_counters->GetXaxis()->SetBinLabel(1,"MatchEff");
  H_vx_counters->GetXaxis()->SetBinLabel(2,"MatchPur");
  H_vx_counters->GetXaxis()->SetBinLabel(3,"Pixel");
  H_vx_counters->GetXaxis()->SetBinLabel(4,"Strip");

  H_vx_pur_eff = SubDirPixVxCom->make<TH1F>("Purity Efficiency Vertices","Purity Efficiency Vertices", 2, 0., 2.);
  H_vx_pur_eff->SetXTitle("");
  H_vx_pur_eff->SetYTitle("Entries (#)");
  H_vx_pur_eff->GetXaxis()->SetBinLabel(1,"Purity");
  H_vx_pur_eff->GetXaxis()->SetBinLabel(2,"Efficiency");

  H_pix_vx_trk_pt = SubDirPixVxCom->make<TH1F>("pt distribution of pixelvertices tracks","p_{t} distribution of pixel-vertices tracks", 100, 0., 10.);
  H_pix_vx_trk_pt->SetXTitle("Track p_{t} (GeV/c)");
  H_pix_vx_trk_pt->SetYTitle("Entries (#)");

  H_pix_vx_trk_pt_nomatch = SubDirPixVxCom->make<TH1F>("pt distribution of pixelvertices nomatch tracks","p_{t} distribution of pixel-vertices tracks not matched", 100, 0., 10.);
  H_pix_vx_trk_pt_nomatch->SetXTitle("Track p_{t} (GeV/c)");
  H_pix_vx_trk_pt_nomatch->SetYTitle("Entries (#)");

  H_pix_vx_trk_eta_nomatch = SubDirPixVxCom->make<TH1F>("eta distribution of pixelvertices nomatch tracks","\\eta distribution of pixel-vertices tracks not matched", 100, -3.5, 3.5);
  H_pix_vx_trk_eta_nomatch->SetXTitle("Track \\eta");
  H_pix_vx_trk_eta_nomatch->SetYTitle("Entries (#)");

  H_pix_vx_trk_normchi2 = SubDirPixVxCom->make<TH1F>("Chi2DoF distribution of pixelvertices tracks","\\chi^{2}/DoF distribution of pixel-vertices tracks", 100, 0., 30.);
  H_pix_vx_trk_normchi2->SetXTitle("Track \\chi^{2}/DoF");
  H_pix_vx_trk_normchi2->SetYTitle("Entries (#)");

  H_pix_vx_trk_normchi2_nomatch = SubDirPixVxCom->make<TH1F>("Chi2DoF distribution of pixelvertices nomatch tracks","\\chi^{2}/DoF distribution of pixel-vertices tracks not matched", 100, 0., 30.);
  H_pix_vx_trk_normchi2_nomatch->SetXTitle("Track \\chi^{2}/DoF");
  H_pix_vx_trk_normchi2_nomatch->SetYTitle("Entries (#)");

  H_trk_vx_trk_pt = SubDirPixVxCom->make<TH1F>("pt distribution of trackervertices tracks","p_{t} distribution of tracker-vertices tracks", 100, 0., 10.);
  H_trk_vx_trk_pt->SetXTitle("Track p_{t} (GeV/c)");
  H_trk_vx_trk_pt->SetYTitle("Entries (#)");

  H_trk_vx_trk_pt_nomatch = SubDirPixVxCom->make<TH1F>("pt distribution of trackervertices nomatch tracks","p_{t} distribution of tracker-vertices tracks not matched", 100, 0., 10.);
  H_trk_vx_trk_pt_nomatch->SetXTitle("Track p_{t} (GeV/c)");
  H_trk_vx_trk_pt_nomatch->SetYTitle("Entries (#)");

  H_trk_vxrefit_trk_pt_nomatch = SubDirPixVxCom->make<TH1F>("pt distribution of refitted trackervertices nomatch tracks","p_{t} distribution of refitted tracker-vertices tracks not matched", 100, 0., 10.);
  H_trk_vxrefit_trk_pt_nomatch->SetXTitle("Track p_{t} (GeV/c)");
  H_trk_vxrefit_trk_pt_nomatch->SetYTitle("Entries (#)");

  H_trk_vx_trk_eta_nomatch = SubDirPixVxCom->make<TH1F>("eta distribution of trackervertices nomatch tracks","\\eta distribution of tracker-vertices tracks not matched", 100, -3.5, 3.5);
  H_trk_vx_trk_eta_nomatch->SetXTitle("Track \\eta");
  H_trk_vx_trk_eta_nomatch->SetYTitle("Entries (#)");

  H_trk_vxrefit_trk_eta_nomatch = SubDirPixVxCom->make<TH1F>("eta distribution of refitted trackervertices nomatch tracks","\\eta distribution of refitted tracker-vertices tracks not matched", 100, -3.5, 3.5);
  H_trk_vxrefit_trk_eta_nomatch->SetXTitle("Track \\eta");
  H_trk_vxrefit_trk_eta_nomatch->SetYTitle("Entries (#)");

  H_trk_vx_trk_normchi2 = SubDirPixVxCom->make<TH1F>("Chi2DoF distribution of trackervertices tracks","\\chi^{2}/DoF distribution of tracker-vertices tracks", 100, 0., 30.);
  H_trk_vx_trk_normchi2->SetXTitle("Track \\chi^{2}/DoF");
  H_trk_vx_trk_normchi2->SetYTitle("Entries (#)");

  H_trk_vx_trk_normchi2_nomatch = SubDirPixVxCom->make<TH1F>("Chi2DoF distribution of trackervertices nomatch tracks","\\chi^{2}/DoF distribution of tracker-vertices tracks not matched", 100, 0., 30.);
  H_trk_vx_trk_normchi2_nomatch->SetXTitle("Track \\chi^{2}/DoF");
  H_trk_vx_trk_normchi2_nomatch->SetYTitle("Entries (#)");
}


void MyPixAnalyzer::endJob()
{
  cout << "\n@@@ Counters @@@" << endl;
  cout << "Total events: " << EvTotal << endl;
  cout << "L1 pass: " << EvPassL1 << endl;
  cout << "L1 discard: " << EvDiscardL1 << endl;
  cout << "HLT pass: " << EvPassHLT << endl;
  cout << "HLT discard: " << EvDiscardHLT << endl;
  cout << "Total pass: " << EvPassTotal << endl;
  cout << "Total discard: " << EvDiscardTotal << endl;
  cout << "@@@@@@@@@@@@@@@@" << endl;

  H_counters->Fill("Total events",EvTotal);
  H_counters->Fill("L1 pass",EvPassL1);
  H_counters->Fill("L1 discard",EvDiscardL1);
  H_counters->Fill("HLT pass",EvPassHLT);
  H_counters->Fill("HLT discard",EvDiscardHLT);
  H_counters->Fill("Total pass",EvPassTotal);
  H_counters->Fill("Total discard",EvDiscardTotal);

  if (IsTrkPart == true) delete theTkAssociator;
}


// Define this as a plug-in
DEFINE_FWK_MODULE(MyPixAnalyzer);

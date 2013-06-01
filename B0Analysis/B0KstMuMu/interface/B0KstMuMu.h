// -*- C++ -*-
//
// Package:    B0KstMuMu
// Class:      B0KstMuMu
// 
/**\class B0KstMuMu B0KstMuMu.cc B0Analysis/B0KstMuMu/src/B0KstMuMu.cc
   
   Description:
   Make rootTuple for b --> s mu+ mu- analyses
   
*/
//
// Original Author:  Mauro Dinardo
//         Created:  Mon Apr 27 09:53:19 MDT 2011
// $Id: B0KstMuMu.h,v 1.30 2012/12/01 14:15:11 dinardo Exp $


#ifndef _B0KSTMUMU
#define _B0KSTMUMU

// System include files
#include <memory>

// User include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <vector>
#include <string>
#include "TTree.h"

#include "ReadParameters.h"
#include "Utils.h"
#include "B0KstMuMuTreeContent.h"


class B0KstMuMu : public edm::EDAnalyzer
{
  
  
public:
  explicit B0KstMuMu (const edm::ParameterSet&);
  ~B0KstMuMu ();
  
  std::string getMuCat (reco::Muon const& muon);

  const reco::Candidate* skipOscillations (const reco::Candidate* Mother);

  const reco::Candidate* findDaughter (edm::Handle<reco::GenParticleCollection> genParticles,
				       std::vector<std::vector<unsigned int>* >* Dau,
				       unsigned int it);
  
  void searchForStableChargedDaughters (const reco::Candidate* Mother,
					unsigned int motherIndex,
					std::vector<std::vector<unsigned int>* >* posDau,
					std::vector<std::vector<unsigned int>* >* negDau);
  
  
private:
  virtual void beginJob ();
  virtual void analyze (const edm::Event&, const edm::EventSetup&);
  virtual void endJob ();

  // ####################
  // # HLT-trigger cuts #
  // ####################
  double CLMUMUVTX;
  double LSMUMUBS;
  double DCAMUMU;
  double DCAMUBS;
  double COSALPHAMUMUBS;
  double MUMINPT;
  double MUMAXETA;
  double MINMUMUPT;
  double MINMUMUINVMASS;
  double MAXMUMUINVMASS;

  // ######################
  // # Pre-selection cuts #
  // ######################
  double B0MASSUPLIMIT;
  double B0MASSLOWLIMIT;
  double CLB0VTX;
  double KSTMASSWINDOW;
  double HADDCASBS;
  double MINHADPT;

  std::string hltTriggerResults_;
  std::string vtxSample_;
  std::string beamSpot_;
  std::string genParticles_;
  std::string muonType_;
  std::string trackType_;
  std::string parameterFile_;
  unsigned int doGenReco_;
  bool PrintMsg;
  
  ReadParameters* ParameterFile;
  std::vector<std::string> TrigTable;
  
  TTree* theTree;
  B0KstMuMuTreeContent* NTuple;
  Utils* Utility;
};

#endif

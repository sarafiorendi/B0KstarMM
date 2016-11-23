#include "B0KstMuMuTreeContent.h"
#include <iostream>

B0KstMuMuTreeContent::B0KstMuMuTreeContent ()
{
  ClearScalars();

  // ### Trigger ###
  TrigTable     = NULL;
  TrigPrescales = NULL;

  // ### B0 Mass ###
  bMass      = NULL;
  bMassE     = NULL;
  bBarMass   = NULL;
  bBarMassE  = NULL;
  bPx        = NULL;
  bPy        = NULL;
  bPz        = NULL;

  // ### Pileup information in MC ###
  bunchXingMC           = NULL;
  numInteractionsMC     = NULL;
  trueNumInteractionsMC = NULL;

  // ### B0 Vtx ###
  bVtxCL        = NULL;
  bVtxX         = NULL;
  bVtxY         = NULL;
  bVtxZ         = NULL;
  bCosAlphaVtx  = NULL;
  bCosAlphaVtxE = NULL;
  bCosAlphaBS   = NULL;
  bCosAlphaBSE  = NULL;
  bLVtx         = NULL;
  bLVtxE        = NULL;
  bLBS          = NULL;
  bLBSE         = NULL;
  bDCAVtx       = NULL;
  bDCAVtxE      = NULL;
  bDCABS        = NULL;
  bDCABSE       = NULL;

  // ### B0 ctau ###
  bctauPVBS  = NULL;
  bctauPVBSE = NULL;

  // ### K*0 Mass ###
  kstMass     = NULL;
  kstMassE    = NULL;
  kstBarMass  = NULL;
  kstBarMassE = NULL;
  kstPx       = NULL;
  kstPy       = NULL;
  kstPz       = NULL;

  // ### K*0 Vtx ###
  kstVtxCL = NULL;
  kstVtxX  = NULL;
  kstVtxY  = NULL;
  kstVtxZ  = NULL;

  // ### mu+ mu- Mass ###
  mumuMass  = NULL;
  mumuMassE = NULL;
  mumuPx    = NULL;
  mumuPy    = NULL;
  mumuPz    = NULL;

  // ### mu+ mu- Vtx ###
  mumuVtxCL       = NULL;
  mumuVtxX        = NULL;
  mumuVtxY        = NULL;
  mumuVtxZ        = NULL;
  mumuCosAlphaBS  = NULL;
  mumuCosAlphaBSE = NULL; 
  mumuLBS         = NULL;
  mumuLBSE        = NULL;
  mumuDCA         = NULL;

  // ### mu- ###
  mumHighPurity    = NULL;
  mumCL            = NULL;
  mumNormChi2      = NULL;
  mumPx            = NULL;
  mumPy            = NULL;
  mumPz            = NULL;
  mumDCAVtx        = NULL;
  mumDCAVtxE       = NULL;
  mumDCABS         = NULL;
  mumDCABSE        = NULL;
  mumKinkChi2      = NULL;
  mumFracHits      = NULL;
  mumdxyVtx        = NULL;
  mumdzVtx         = NULL;
  mumDeltaRwithMC  = NULL;
  mumCat           = NULL;
  mumNPixHits      = NULL;
  mumNPixLayers    = NULL;
  mumNTrkHits      = NULL;
  mumNTrkLayers    = NULL;
  mumNMuonHits     = NULL;
  mumNMatchStation = NULL;
  mumTrig          = NULL;

  // ### mu+ ###
  mupHighPurity    = NULL;
  mupCL            = NULL;
  mupNormChi2      = NULL;
  mupPx            = NULL;
  mupPy            = NULL;
  mupPz            = NULL;
  mupDCAVtx        = NULL;
  mupDCAVtxE       = NULL;
  mupDCABS         = NULL;
  mupDCABSE        = NULL;
  mupKinkChi2      = NULL;
  mupFracHits      = NULL;
  mupdxyVtx        = NULL;
  mupdzVtx         = NULL;
  mupDeltaRwithMC  = NULL;
  mupCat           = NULL;
  mupNPixHits      = NULL;
  mupNPixLayers    = NULL;
  mupNTrkHits      = NULL;
  mupNTrkLayers    = NULL;
  mupNMuonHits     = NULL;
  mupNMatchStation = NULL;
  mupTrig          = NULL;

  // ### K*0 track- ###
  kstTrkmHighPurity   = NULL;
  kstTrkmCL           = NULL;
  kstTrkmNormChi2     = NULL;
  kstTrkmPx           = NULL;
  kstTrkmPy           = NULL;
  kstTrkmPz           = NULL;
  kstTrkmDCAVtx       = NULL;
  kstTrkmDCAVtxE      = NULL;
  kstTrkmDCABS        = NULL;
  kstTrkmDCABSE       = NULL;
  kstTrkmFracHits     = NULL;
  kstTrkmdxyVtx       = NULL;
  kstTrkmdzVtx        = NULL;
  kstTrkmDeltaRwithMC = NULL;
  kstTrkmNPixHits     = NULL;
  kstTrkmNPixLayers   = NULL;
  kstTrkmNTrkHits     = NULL;
  kstTrkmNTrkLayers   = NULL;
  kstTrkmMuMatch      = NULL;

  // ### K*0 track+ ###
  kstTrkpHighPurity   = NULL;
  kstTrkpCL           = NULL;
  kstTrkpNormChi2     = NULL;
  kstTrkpPx           = NULL;
  kstTrkpPy           = NULL;
  kstTrkpPz           = NULL;
  kstTrkpDCAVtx       = NULL;
  kstTrkpDCAVtxE      = NULL;
  kstTrkpDCABS        = NULL;
  kstTrkpDCABSE       = NULL;
  kstTrkpFracHits     = NULL;
  kstTrkpdxyVtx       = NULL;
  kstTrkpdzVtx        = NULL;
  kstTrkpDeltaRwithMC = NULL;
  kstTrkpNPixHits     = NULL;
  kstTrkpNPixLayers   = NULL;
  kstTrkpNTrkHits     = NULL;
  kstTrkpNTrkLayers   = NULL;
  kstTrkpMuMatch      = NULL;

  // ### Matching Between Reconstructed and Generated ###
  truthMatchSignal = NULL;
  truthMatchMum    = NULL;
  truthMatchMup    = NULL;
  truthMatchTrkm   = NULL;
  truthMatchTrkp   = NULL;
}

void B0KstMuMuTreeContent::Init ()
{
  // ### Trigger ###
  TrigTable     = new std::vector<std::string>;
  TrigPrescales = new std::vector<int>;

  // ### B0 Mass ###
  bMass      = new std::vector<double>;
  bMassE     = new std::vector<double>;
  bBarMass   = new std::vector<double>;
  bBarMassE  = new std::vector<double>;
  bPx        = new std::vector<double>;
  bPy        = new std::vector<double>;
  bPz        = new std::vector<double>;

  // ### Pileup information in MC ###
  bunchXingMC           = new std::vector<double>;
  numInteractionsMC     = new std::vector<double>;
  trueNumInteractionsMC = new std::vector<double>;

  // ### B0 Vtx ###
  bVtxCL        = new std::vector<double>;
  bVtxX         = new std::vector<double>;
  bVtxY         = new std::vector<double>;
  bVtxZ         = new std::vector<double>;
  bCosAlphaVtx  = new std::vector<double>;
  bCosAlphaVtxE = new std::vector<double>;
  bCosAlphaBS   = new std::vector<double>;
  bCosAlphaBSE  = new std::vector<double>;
  bLVtx         = new std::vector<double>;
  bLVtxE        = new std::vector<double>;
  bLBS          = new std::vector<double>;
  bLBSE         = new std::vector<double>;
  bDCAVtx       = new std::vector<double>;
  bDCAVtxE      = new std::vector<double>;
  bDCABS        = new std::vector<double>;
  bDCABSE       = new std::vector<double>;

  // ### B0 ctau ###
  bctauPVBS  = new std::vector<double>;
  bctauPVBSE = new std::vector<double>;

  // ### K*0 Mass ###
  kstMass     = new std::vector<double>;
  kstMassE    = new std::vector<double>;
  kstBarMass  = new std::vector<double>;
  kstBarMassE = new std::vector<double>;
  kstPx       = new std::vector<double>;
  kstPy       = new std::vector<double>;
  kstPz       = new std::vector<double>;

  // ### K*0 Vtx ###
  kstVtxCL = new std::vector<double>;
  kstVtxX  = new std::vector<double>;
  kstVtxY  = new std::vector<double>;
  kstVtxZ  = new std::vector<double>;

  // ### mu+ mu- Mass ###
  mumuMass  = new std::vector<double>;
  mumuMassE = new std::vector<double>;
  mumuPx    = new std::vector<double>;
  mumuPy    = new std::vector<double>;
  mumuPz    = new std::vector<double>;

  // ### mu+ mu- Vtx ###
  mumuVtxCL       = new std::vector<double>;
  mumuVtxX        = new std::vector<double>;
  mumuVtxY        = new std::vector<double>;
  mumuVtxZ        = new std::vector<double>;
  mumuCosAlphaBS  = new std::vector<double>;
  mumuCosAlphaBSE = new std::vector<double>; 
  mumuLBS         = new std::vector<double>;
  mumuLBSE        = new std::vector<double>;
  mumuDCA         = new std::vector<double>;

  // ### mu- ###
  mumHighPurity    = new std::vector<bool>;
  mumCL            = new std::vector<double>;
  mumNormChi2      = new std::vector<double>;
  mumPx            = new std::vector<double>;
  mumPy            = new std::vector<double>;
  mumPz            = new std::vector<double>;
  mumDCAVtx        = new std::vector<double>;
  mumDCAVtxE       = new std::vector<double>;
  mumDCABS         = new std::vector<double>;
  mumDCABSE        = new std::vector<double>;
  mumKinkChi2      = new std::vector<double>;
  mumFracHits      = new std::vector<double>;
  mumdxyVtx        = new std::vector<double>;
  mumdzVtx         = new std::vector<double>;
  mumDeltaRwithMC  = new std::vector<double>;
  mumCat           = new std::vector<std::string>;
  mumNPixHits      = new std::vector<int>;
  mumNPixLayers    = new std::vector<int>;
  mumNTrkHits      = new std::vector<int>;
  mumNTrkLayers    = new std::vector<int>;
  mumNMuonHits     = new std::vector<int>;
  mumNMatchStation = new std::vector<int>;
  mumTrig          = new std::vector<std::string>;

  // ### mu+ ###
  mupHighPurity    = new std::vector<bool>;
  mupCL            = new std::vector<double>; 
  mupNormChi2      = new std::vector<double>;
  mupPx            = new std::vector<double>;
  mupPy            = new std::vector<double>;
  mupPz            = new std::vector<double>;
  mupDCAVtx        = new std::vector<double>;
  mupDCAVtxE       = new std::vector<double>;
  mupDCABS         = new std::vector<double>;
  mupDCABSE        = new std::vector<double>;
  mupKinkChi2      = new std::vector<double>;
  mupFracHits      = new std::vector<double>;
  mupdxyVtx        = new std::vector<double>;
  mupdzVtx         = new std::vector<double>;
  mupDeltaRwithMC  = new std::vector<double>;
  mupCat           = new std::vector<std::string>;
  mupNPixHits      = new std::vector<int>;
  mupNPixLayers    = new std::vector<int>;
  mupNTrkHits      = new std::vector<int>;
  mupNTrkLayers    = new std::vector<int>;
  mupNMuonHits     = new std::vector<int>;
  mupNMatchStation = new std::vector<int>;
  mupTrig          = new std::vector<std::string>;

  // ### K*0 track- ###
  kstTrkmHighPurity   = new std::vector<bool>;
  kstTrkmCL           = new std::vector<double>;
  kstTrkmNormChi2     = new std::vector<double>;
  kstTrkmPx           = new std::vector<double>;
  kstTrkmPy           = new std::vector<double>;
  kstTrkmPz           = new std::vector<double>;
  kstTrkmDCAVtx       = new std::vector<double>;
  kstTrkmDCAVtxE      = new std::vector<double>;
  kstTrkmDCABS        = new std::vector<double>;
  kstTrkmDCABSE       = new std::vector<double>;
  kstTrkmFracHits     = new std::vector<double>;
  kstTrkmdxyVtx       = new std::vector<double>;
  kstTrkmdzVtx        = new std::vector<double>;
  kstTrkmDeltaRwithMC = new std::vector<double>;
  kstTrkmNPixHits     = new std::vector<int>;
  kstTrkmNPixLayers   = new std::vector<int>;
  kstTrkmNTrkHits     = new std::vector<int>;
  kstTrkmNTrkLayers   = new std::vector<int>;
  kstTrkmMuMatch      = new std::vector<std::string>;

  // ### K*0 track+ ###
  kstTrkpHighPurity   = new std::vector<bool>;
  kstTrkpCL           = new std::vector<double>;
  kstTrkpNormChi2     = new std::vector<double>;
  kstTrkpPx           = new std::vector<double>;
  kstTrkpPy           = new std::vector<double>;
  kstTrkpPz           = new std::vector<double>;
  kstTrkpDCAVtx       = new std::vector<double>;
  kstTrkpDCAVtxE      = new std::vector<double>;
  kstTrkpDCABS        = new std::vector<double>;
  kstTrkpDCABSE       = new std::vector<double>;
  kstTrkpFracHits     = new std::vector<double>;
  kstTrkpdxyVtx       = new std::vector<double>;
  kstTrkpdzVtx        = new std::vector<double>;
  kstTrkpDeltaRwithMC = new std::vector<double>;
  kstTrkpNPixHits     = new std::vector<int>;
  kstTrkpNPixLayers   = new std::vector<int>;
  kstTrkpNTrkHits     = new std::vector<int>;
  kstTrkpNTrkLayers   = new std::vector<int>;
  kstTrkpMuMatch      = new std::vector<std::string>;

  // ### Matching Between Reconstructed and Generated ###
  truthMatchSignal = new std::vector<bool>;
  truthMatchMum    = new std::vector<bool>;
  truthMatchMup    = new std::vector<bool>;
  truthMatchTrkm   = new std::vector<bool>;
  truthMatchTrkp   = new std::vector<bool>;
}

B0KstMuMuTreeContent::~B0KstMuMuTreeContent ()
{
  // ### Trigger ###
  delete TrigTable;
  delete TrigPrescales;

  // ### B0 Mass ###
  delete bMass;
  delete bMassE;
  delete bBarMass;
  delete bBarMassE;
  delete bPx;
  delete bPy;
  delete bPz;

  // ### Pileup information in MC ###
  delete bunchXingMC;
  delete numInteractionsMC;
  delete trueNumInteractionsMC;

  // ### B0 Vtx ###
  delete bVtxCL;
  delete bVtxX;
  delete bVtxY;
  delete bVtxZ;
  delete bCosAlphaVtx;
  delete bCosAlphaVtxE;
  delete bCosAlphaBS;
  delete bCosAlphaBSE;
  delete bLVtx;
  delete bLVtxE;
  delete bLBS;
  delete bLBSE;
  delete bDCAVtx;
  delete bDCAVtxE;
  delete bDCABS;
  delete bDCABSE;

  // ### B0 ctau ###
  delete bctauPVBS;
  delete bctauPVBSE;

  // ### K*0 Mass ###
  delete kstMass;
  delete kstMassE;
  delete kstBarMass;
  delete kstBarMassE;
  delete kstPx;
  delete kstPy;
  delete kstPz;

  // ### K*0 Vtx ###
  delete kstVtxCL;
  delete kstVtxX;
  delete kstVtxY;
  delete kstVtxZ;

  // ### mu+ mu- Mass ###
  delete mumuMass;
  delete mumuMassE;
  delete mumuPx;
  delete mumuPy;
  delete mumuPz;

  // ### mu+ mu- Vtx ###
  delete mumuVtxCL;
  delete mumuVtxX;
  delete mumuVtxY;
  delete mumuVtxZ;
  delete mumuCosAlphaBS;
  delete mumuCosAlphaBSE;
  delete mumuLBS;
  delete mumuLBSE;
  delete mumuDCA;

  // ### mu- ###
  delete mumHighPurity;
  delete mumCL;
  delete mumNormChi2;
  delete mumPx;
  delete mumPy;
  delete mumPz;
  delete mumDCAVtx;
  delete mumDCAVtxE;
  delete mumDCABS;
  delete mumDCABSE;
  delete mumKinkChi2;
  delete mumFracHits;
  delete mumdxyVtx;
  delete mumdzVtx;
  delete mumDeltaRwithMC;
  delete mumCat;
  delete mumNPixHits;
  delete mumNPixLayers;
  delete mumNTrkHits;
  delete mumNTrkLayers;
  delete mumNMuonHits;
  delete mumNMatchStation;
  delete mumTrig;

  // ### mu+ ###
  delete mupHighPurity;
  delete mupCL;
  delete mupNormChi2;
  delete mupPx;
  delete mupPy;
  delete mupPz;
  delete mupDCAVtx;
  delete mupDCAVtxE;
  delete mupDCABS;
  delete mupDCABSE;
  delete mupKinkChi2;
  delete mupFracHits;
  delete mupdxyVtx;
  delete mupdzVtx;
  delete mupDeltaRwithMC;
  delete mupCat;
  delete mupNPixHits;
  delete mupNPixLayers;
  delete mupNTrkHits;
  delete mupNTrkLayers;
  delete mupNMuonHits;
  delete mupNMatchStation;
  delete mupTrig;
      
  // ### K*0 track- ###
  delete kstTrkmHighPurity;
  delete kstTrkmCL;
  delete kstTrkmNormChi2;
  delete kstTrkmPx;
  delete kstTrkmPy;
  delete kstTrkmPz;
  delete kstTrkmDCAVtx;
  delete kstTrkmDCAVtxE;
  delete kstTrkmDCABS;
  delete kstTrkmDCABSE;
  delete kstTrkmFracHits;
  delete kstTrkmdxyVtx;
  delete kstTrkmdzVtx;
  delete kstTrkmDeltaRwithMC;
  delete kstTrkmNPixHits;
  delete kstTrkmNPixLayers;
  delete kstTrkmNTrkHits;
  delete kstTrkmNTrkLayers;
  delete kstTrkmMuMatch;

  // ### K*0 track+ ###
  delete kstTrkpHighPurity;
  delete kstTrkpCL;
  delete kstTrkpNormChi2;
  delete kstTrkpPx;
  delete kstTrkpPy;
  delete kstTrkpPz;
  delete kstTrkpDCAVtx;
  delete kstTrkpDCAVtxE;
  delete kstTrkpDCABS;
  delete kstTrkpDCABSE;
  delete kstTrkpFracHits;
  delete kstTrkpdxyVtx;
  delete kstTrkpdzVtx;
  delete kstTrkpDeltaRwithMC;
  delete kstTrkpNPixHits;
  delete kstTrkpNPixLayers;
  delete kstTrkpNTrkHits;
  delete kstTrkpNTrkLayers;
  delete kstTrkpMuMatch;

  // ### Matching Between Reconstructed and Generated ###
  delete truthMatchSignal;
  delete truthMatchMum;
  delete truthMatchMup;
  delete truthMatchTrkm;
  delete truthMatchTrkp;
}

void B0KstMuMuTreeContent::ClearScalars ()
{
  // ########################################################
  // # Run Number, event number, #reco vtx and event weight #
  // ########################################################
  runN            = 0;
  eventN          = 0;
  recoVtxN        = 0;
  evWeight        = 1;
  evWeightE2      = 0;
  numEventsTried  = 0;
  numEventsPassed = 0;

  nB = 0;
  
  // ### Primary Vertex and Beam Spot ###
  priVtxCL = 0;
  priVtxX  = 0;
  priVtxY  = 0;
  priVtxZ  = 0;
  bsX      = 0;
  bsY      = 0;
  
  ClearScalarsMonteCarlo();
}

void B0KstMuMuTreeContent::ClearScalarsMonteCarlo ()
{
  // ### Generated Observables ###
  genSignal       = 0;
  genMuMuBG       = 0;
  genMuMuBGnTrksm = 0;
  genMuMuBGnTrksp = 0;
  genPsiPrompt     = false;
  genSignHasFSR    = false;
  genSignKstHasFSR = false;
  genSignPsiHasFSR = false;

  // # Generated Primary Vertex #
  genPriVtxX = 0;
  genPriVtxY = 0;
  genPriVtxZ = 0;

  // ### Generated B0 Mass ###
  genB0Mass = 0;
  genB0Px   = 0;
  genB0Py   = 0;
  genB0Pz   = 0;

  // ### Generated B0 Vtx ###
  genB0VtxX = 0;
  genB0VtxY = 0;
  genB0VtxZ = 0;

  // ### Generated K*0 Mass ###
  genKstMass = 0;
  genKstPx   = 0;
  genKstPy   = 0;
  genKstPz   = 0;

  // ### Generated K*0 Vtx ###
  genKstVtxX = 0;
  genKstVtxY = 0;
  genKstVtxZ = 0;

  // ### Generated J/psi or psi(2S) Mass and Vtx ###
  genPsiMass = 0;
  genPsiVtxX = 0;
  genPsiVtxY = 0;
  genPsiVtxZ = 0;

  // ### Generated mu- ###
  genMumMother = 0;
  genMumPx     = 0;
  genMumPy     = 0;
  genMumPz     = 0;

  // ### Generated mu+ ###
  genMupMother = 0;
  genMupPx     = 0;
  genMupPy     = 0;
  genMupPz     = 0;

  // ### Generated K*0 track- ###
  genKstTrkmMother = 0;
  genKstTrkmID     = 0;
  genKstTrkmPx     = 0;
  genKstTrkmPy     = 0;
  genKstTrkmPz     = 0;

  // ### Generated K*0 track+ ###
  genKstTrkpMother = 0;
  genKstTrkpID     = 0;
  genKstTrkpPx     = 0;
  genKstTrkpPy     = 0;
  genKstTrkpPz     = 0;
}

void B0KstMuMuTreeContent::ClearVectors ()
{
  // ### Trigger ###
  TrigTable->clear();
  TrigPrescales->clear();

  // ### B0 Mass ###
  bMass->clear();
  bMassE->clear();
  bBarMass->clear();
  bBarMassE->clear();
  bPx->clear();
  bPy->clear();
  bPz->clear();

  // ### Pileup information in MC ###
  bunchXingMC->clear();
  numInteractionsMC->clear();
  trueNumInteractionsMC->clear();

  // ### B0 Vtx ###
  bVtxCL->clear();
  bVtxX->clear();
  bVtxY->clear();
  bVtxZ->clear();
  bCosAlphaVtx->clear();
  bCosAlphaVtxE->clear();
  bCosAlphaBS->clear();
  bCosAlphaBSE->clear();
  bLVtx->clear();
  bLVtxE->clear();
  bLBS->clear();
  bLBSE->clear();
  bDCAVtx->clear();
  bDCAVtxE->clear();
  bDCABS->clear();
  bDCABSE->clear();

  // ### B0 ctau ###
  bctauPVBS->clear();
  bctauPVBSE->clear();

  // ### K*0 Mass ###
  kstMass->clear();
  kstMassE->clear();
  kstBarMass->clear();
  kstBarMassE->clear();
  kstPx->clear();
  kstPy->clear();
  kstPz->clear();

  // ### K*0 Vtx ###
  kstVtxCL->clear();
  kstVtxX->clear();
  kstVtxY->clear();
  kstVtxZ->clear();

  // ### mu+ mu- Mass ###
  mumuMass->clear();
  mumuMassE->clear();
  mumuPx->clear();
  mumuPy->clear();
  mumuPz->clear();

  // ### mu+ mu- Vtx ###
  mumuVtxCL->clear();
  mumuVtxX->clear();
  mumuVtxY->clear();
  mumuVtxZ->clear();
  mumuCosAlphaBS->clear();
  mumuCosAlphaBSE->clear();
  mumuLBS->clear();
  mumuLBSE->clear();
  mumuDCA->clear();

  // ### mu- ###
  mumHighPurity->clear();
  mumCL->clear();
  mumNormChi2->clear();
  mumPx->clear();
  mumPy->clear();
  mumPz->clear();
  mumDCAVtx->clear();
  mumDCAVtxE->clear();
  mumDCABS->clear();
  mumDCABSE->clear();
  mumKinkChi2->clear();
  mumFracHits->clear();
  mumdxyVtx->clear();
  mumdzVtx->clear();
  mumDeltaRwithMC->clear();
  mumCat->clear();
  mumNPixHits->clear();
  mumNPixLayers->clear();
  mumNTrkHits->clear();
  mumNTrkLayers->clear();
  mumNMuonHits->clear();
  mumNMatchStation->clear();
  mumTrig->clear();

  // ### mu+ ###
  mupHighPurity->clear();
  mupCL->clear();
  mupNormChi2->clear();
  mupPx->clear();
  mupPy->clear();
  mupPz->clear();
  mupDCAVtx->clear();
  mupDCAVtxE->clear();
  mupDCABS->clear();
  mupDCABSE->clear();
  mupKinkChi2->clear();
  mupFracHits->clear();
  mupdxyVtx->clear();
  mupdzVtx->clear();
  mupDeltaRwithMC->clear();
  mupCat->clear();
  mupNPixHits->clear();
  mupNPixLayers->clear();
  mupNTrkHits->clear();
  mupNTrkLayers->clear();
  mupNMuonHits->clear();
  mupNMatchStation->clear();
  mupTrig->clear();

  // ### K*0 track- ###
  kstTrkmHighPurity->clear();
  kstTrkmCL->clear();
  kstTrkmNormChi2->clear();
  kstTrkmPx->clear();
  kstTrkmPy->clear();
  kstTrkmPz->clear();
  kstTrkmDCAVtx->clear();
  kstTrkmDCAVtxE->clear();
  kstTrkmDCABS->clear();
  kstTrkmDCABSE->clear();
  kstTrkmFracHits->clear();
  kstTrkmdxyVtx->clear();
  kstTrkmdzVtx->clear(); 
  kstTrkmDeltaRwithMC->clear();
  kstTrkmNPixHits->clear();
  kstTrkmNPixLayers->clear();
  kstTrkmNTrkHits->clear();
  kstTrkmNTrkLayers->clear();
  kstTrkmMuMatch->clear();

  // ### K*0 track+ ###
  kstTrkpHighPurity->clear();
  kstTrkpCL->clear();
  kstTrkpNormChi2->clear();
  kstTrkpPx->clear();
  kstTrkpPy->clear();
  kstTrkpPz->clear();
  kstTrkpDCAVtx->clear();
  kstTrkpDCAVtxE->clear();
  kstTrkpDCABS->clear();
  kstTrkpDCABSE->clear();
  kstTrkpFracHits->clear();
  kstTrkpdxyVtx->clear();
  kstTrkpdzVtx->clear();
  kstTrkpDeltaRwithMC->clear();
  kstTrkpNPixHits->clear();
  kstTrkpNPixLayers->clear();
  kstTrkpNTrkHits->clear();
  kstTrkpNTrkLayers->clear();
  kstTrkpMuMatch->clear();

  ClearVectorsMonteCarlo();
}

void B0KstMuMuTreeContent::ClearVectorsMonteCarlo ()
{
  // ### Matching Between Reconstructed and Generated ###
  truthMatchSignal->clear();
  truthMatchMum->clear();
  truthMatchMup->clear();
  truthMatchTrkm->clear();
  truthMatchTrkp->clear();
}

void B0KstMuMuTreeContent::ClearNTuple ()
{
  ClearScalars();
  ClearVectors();
}

void B0KstMuMuTreeContent::ClearMonteCarlo ()
{
  ClearScalarsMonteCarlo();
  ClearVectorsMonteCarlo();
}

void B0KstMuMuTreeContent::MakeTreeBranches (TTree* theTree)
{
  // ########################################################
  // # Run Number, event number, #reco vtx and event weight #
  // ########################################################
  theTree->Branch("runN",            &runN,            "runN/i");
  theTree->Branch("eventN",          &eventN,          "eventN/i");
  theTree->Branch("recoVtxN",        &recoVtxN,        "recoVtxN/i");
  theTree->Branch("evWeight",        &evWeight,        "evWeight/D");
  theTree->Branch("evWeightE2",      &evWeightE2,      "evWeightE2/D");
  theTree->Branch("numEventsTried",  &numEventsTried,  "numEventsTried/i");
  theTree->Branch("numEventsPassed", &numEventsPassed, "numEventsPassed/i");

  // ### Trigger ###
  theTree->Branch("TrigTable",     &TrigTable);
  theTree->Branch("TrigPrescales", &TrigPrescales);

  theTree->Branch("nB", &nB, "nB/i");

  // ### Primary Vertex and Beam Spot ###
  theTree->Branch("priVtxCL", &priVtxCL, "priVtxCL/D");
  theTree->Branch("priVtxX",  &priVtxX,  "priVtxX/D");
  theTree->Branch("priVtxY",  &priVtxY,  "priVtxY/D");
  theTree->Branch("priVtxZ",  &priVtxZ,  "priVtxZ/D");
  theTree->Branch("bsX",      &bsX,      "bsX/D");
  theTree->Branch("bsY",      &bsY,      "bsY/D");

   // ### B0 Mass ###  
  theTree->Branch("bMass",     &bMass);
  theTree->Branch("bMassE",    &bMassE);
  theTree->Branch("bBarMass",  &bBarMass);
  theTree->Branch("bBarMassE", &bBarMassE);
  theTree->Branch("bPx",       &bPx);
  theTree->Branch("bPy",       &bPy);
  theTree->Branch("bPz",       &bPz);

  // ### Pileup information in MC ###
  theTree->Branch("bunchXingMC",           &bunchXingMC);
  theTree->Branch("numInteractionsMC",     &numInteractionsMC);
  theTree->Branch("trueNumInteractionsMC", &trueNumInteractionsMC);

  // ### B0 Vtx ###
  theTree->Branch("bVtxCL",        &bVtxCL);
  theTree->Branch("bVtxX",         &bVtxX);
  theTree->Branch("bVtxY",         &bVtxY);
  theTree->Branch("bVtxZ",         &bVtxZ);
  theTree->Branch("bCosAlphaVtx",  &bCosAlphaVtx);
  theTree->Branch("bCosAlphaVtxE", &bCosAlphaVtxE);
  theTree->Branch("bCosAlphaBS",   &bCosAlphaBS);
  theTree->Branch("bCosAlphaBSE",  &bCosAlphaBSE);
  theTree->Branch("bLVtx",         &bLVtx);
  theTree->Branch("bLVtxE",        &bLVtxE);
  theTree->Branch("bLBS",          &bLBS);
  theTree->Branch("bLBSE",         &bLBSE);
  theTree->Branch("bDCAVtx",       &bDCAVtx);
  theTree->Branch("bDCAVtxE",      &bDCAVtxE);
  theTree->Branch("bDCABS",        &bDCABS);
  theTree->Branch("bDCABSE",       &bDCABSE);

  // ### B0 ctau ###
  theTree->Branch("bctauPVBS",  &bctauPVBS);
  theTree->Branch("bctauPVBSE", &bctauPVBSE);

   // ### K*0 Mass ###
  theTree->Branch("kstMass",     &kstMass);
  theTree->Branch("kstMassE",    &kstMassE);
  theTree->Branch("kstBarMass",  &kstBarMass);
  theTree->Branch("kstBarMassE", &kstBarMassE);
  theTree->Branch("kstPx",       &kstPx);
  theTree->Branch("kstPy",       &kstPy);
  theTree->Branch("kstPz",       &kstPz);

  // ### K*0 Vtx ###
  theTree->Branch("kstVtxCL", &kstVtxCL);
  theTree->Branch("kstVtxX",  &kstVtxX);
  theTree->Branch("kstVtxY",  &kstVtxY);
  theTree->Branch("kstVtxZ",  &kstVtxZ);
  
  // ### mu+ mu- Mass ###
  theTree->Branch("mumuMass",  &mumuMass);
  theTree->Branch("mumuMassE", &mumuMassE);
  theTree->Branch("mumuPx",    &mumuPx);
  theTree->Branch("mumuPy",    &mumuPy);
  theTree->Branch("mumuPz",    &mumuPz);

  // ### mu+ mu- Vtx ###
  theTree->Branch("mumuVtxCL",       &mumuVtxCL);
  theTree->Branch("mumuVtxX",        &mumuVtxX);
  theTree->Branch("mumuVtxY",        &mumuVtxY);
  theTree->Branch("mumuVtxZ",        &mumuVtxZ);
  theTree->Branch("mumuCosAlphaBS",  &mumuCosAlphaBS);
  theTree->Branch("mumuCosAlphaBSE", &mumuCosAlphaBSE);
  theTree->Branch("mumuLBS",         &mumuLBS);
  theTree->Branch("mumuLBSE",        &mumuLBSE);
  theTree->Branch("mumuDCA",         &mumuDCA);

  // ### mu- ###  
  theTree->Branch("mumHighPurity",    &mumHighPurity);
  theTree->Branch("mumCL",            &mumCL);
  theTree->Branch("mumNormChi2",      &mumNormChi2);
  theTree->Branch("mumPx",            &mumPx);
  theTree->Branch("mumPy",            &mumPy);
  theTree->Branch("mumPz",            &mumPz);
  theTree->Branch("mumDCAVtx",        &mumDCAVtx);
  theTree->Branch("mumDCAVtxE",       &mumDCAVtxE);
  theTree->Branch("mumDCABS",         &mumDCABS);
  theTree->Branch("mumDCABSE",        &mumDCABSE);
  theTree->Branch("mumKinkChi2",      &mumKinkChi2);
  theTree->Branch("mumFracHits",      &mumFracHits);
  theTree->Branch("mumdxyVtx",        &mumdxyVtx);
  theTree->Branch("mumdzVtx",         &mumdzVtx);
  theTree->Branch("mumDeltaRwithMC",  &mumDeltaRwithMC);
  theTree->Branch("mumCat",           &mumCat);
  theTree->Branch("mumNPixHits",      &mumNPixHits);
  theTree->Branch("mumNPixLayers",    &mumNPixLayers);
  theTree->Branch("mumNTrkHits",      &mumNTrkHits);
  theTree->Branch("mumNTrkLayers",    &mumNTrkLayers);
  theTree->Branch("mumNMuonHits",     &mumNMuonHits);
  theTree->Branch("mumNMatchStation", &mumNMatchStation);
  theTree->Branch("mumTrig",          &mumTrig);

  // ### mu+ ###  
  theTree->Branch("mupHighPurity",    &mupHighPurity);
  theTree->Branch("mupCL",            &mupCL);
  theTree->Branch("mupNormChi2",      &mupNormChi2);
  theTree->Branch("mupPx",            &mupPx);
  theTree->Branch("mupPy",            &mupPy);
  theTree->Branch("mupPz",            &mupPz);
  theTree->Branch("mupDCAVtx",        &mupDCAVtx);
  theTree->Branch("mupDCAVtxE",       &mupDCAVtxE);
  theTree->Branch("mupDCABS",         &mupDCABS);
  theTree->Branch("mupDCABSE",        &mupDCABSE);
  theTree->Branch("mupKinkChi2",      &mupKinkChi2);
  theTree->Branch("mupFracHits",      &mupFracHits);
  theTree->Branch("mupdxyVtx",        &mupdxyVtx);
  theTree->Branch("mupdzVtx",         &mupdzVtx);
  theTree->Branch("mupDeltaRwithMC",  &mupDeltaRwithMC);
  theTree->Branch("mupCat",           &mupCat);
  theTree->Branch("mupNPixHits",      &mupNPixHits);
  theTree->Branch("mupNPixLayers",    &mupNPixLayers);
  theTree->Branch("mupNTrkHits",      &mupNTrkHits);
  theTree->Branch("mupNTrkLayers",    &mupNTrkLayers);
  theTree->Branch("mupNMuonHits",     &mupNMuonHits);
  theTree->Branch("mupNMatchStation", &mupNMatchStation);
  theTree->Branch("mupTrig",          &mupTrig);

  // ### K*0 track- ###
  theTree->Branch("kstTrkmHighPurity",   &kstTrkmHighPurity);
  theTree->Branch("kstTrkmCL",           &kstTrkmCL);
  theTree->Branch("kstTrkmNormChi2",     &kstTrkmNormChi2);
  theTree->Branch("kstTrkmPx",           &kstTrkmPx);
  theTree->Branch("kstTrkmPy",           &kstTrkmPy);
  theTree->Branch("kstTrkmPz",           &kstTrkmPz);
  theTree->Branch("kstTrkmDCAVtx",       &kstTrkmDCAVtx);
  theTree->Branch("kstTrkmDCAVtxE",      &kstTrkmDCAVtxE);
  theTree->Branch("kstTrkmDCABS",        &kstTrkmDCABS);
  theTree->Branch("kstTrkmDCABSE",       &kstTrkmDCABSE);
  theTree->Branch("kstTrkmFracHits",     &kstTrkmFracHits);
  theTree->Branch("kstTrkmdxyVtx",       &kstTrkmdxyVtx);
  theTree->Branch("kstTrkmdzVtx",        &kstTrkmdzVtx);
  theTree->Branch("kstTrkmDeltaRwithMC", &kstTrkmDeltaRwithMC);
  theTree->Branch("kstTrkmNPixHits",     &kstTrkmNPixHits);
  theTree->Branch("kstTrkmNPixLayers",   &kstTrkmNPixLayers);
  theTree->Branch("kstTrkmNTrkHits",     &kstTrkmNTrkHits);
  theTree->Branch("kstTrkmNTrkLayers",   &kstTrkmNTrkLayers);
  theTree->Branch("kstTrkmMuMatch",      &kstTrkmMuMatch);

  // ### K*0 track+ ### 
  theTree->Branch("kstTrkpHighPurity",   &kstTrkpHighPurity);
  theTree->Branch("kstTrkpCL",           &kstTrkpCL);
  theTree->Branch("kstTrkpNormChi2",     &kstTrkpNormChi2);
  theTree->Branch("kstTrkpPx",           &kstTrkpPx);
  theTree->Branch("kstTrkpPy",           &kstTrkpPy);
  theTree->Branch("kstTrkpPz",           &kstTrkpPz);
  theTree->Branch("kstTrkpDCAVtx",       &kstTrkpDCAVtx);
  theTree->Branch("kstTrkpDCAVtxE",      &kstTrkpDCAVtxE);
  theTree->Branch("kstTrkpDCABS",        &kstTrkpDCABS);
  theTree->Branch("kstTrkpDCABSE",       &kstTrkpDCABSE);
  theTree->Branch("kstTrkpFracHits",     &kstTrkpFracHits);
  theTree->Branch("kstTrkpdxyVtx",       &kstTrkpdxyVtx);
  theTree->Branch("kstTrkpdzVtx",        &kstTrkpdzVtx);
  theTree->Branch("kstTrkpDeltaRwithMC", &kstTrkpDeltaRwithMC);
  theTree->Branch("kstTrkpNPixHits",     &kstTrkpNPixHits);
  theTree->Branch("kstTrkpNPixLayers",   &kstTrkpNPixLayers);
  theTree->Branch("kstTrkpNTrkHits",     &kstTrkpNTrkHits);
  theTree->Branch("kstTrkpNTrkLayers",   &kstTrkpNTrkLayers);
  theTree->Branch("kstTrkpMuMatch",      &kstTrkpMuMatch);
  
  // ### Generated Observables ###
  theTree->Branch("genSignal",        &genSignal,        "genSignal/I");
  theTree->Branch("genMuMuBG",        &genMuMuBG,        "genMuMuBG/I");
  theTree->Branch("genMuMuBGnTrksm",  &genMuMuBGnTrksm,  "genMuMuBGnTrksm/I");
  theTree->Branch("genMuMuBGnTrksp",  &genMuMuBGnTrksp,  "genMuMuBGnTrksp/I");
  theTree->Branch("genPsiPrompt",     &genPsiPrompt,     "genPsiPrompt/O");
  theTree->Branch("genSignHasFSR",    &genSignHasFSR,    "genSignHasFSR/O");
  theTree->Branch("genSignKstHasFSR", &genSignKstHasFSR, "genSignKstHasFSR/O");
  theTree->Branch("genSignPsiHasFSR", &genSignPsiHasFSR, "genSignPsiHasFSR/O");

  // ### Generated Primary Vertex ###
  theTree->Branch("genPriVtxX", &genPriVtxX, "genPriVtxX/D");
  theTree->Branch("genPriVtxY", &genPriVtxY, "genPriVtxY/D");
  theTree->Branch("genPriVtxZ", &genPriVtxZ, "genPriVtxZ/D");

  // ### Generated B0 Mass ###
  theTree->Branch("genB0Mass", &genB0Mass, "genB0Mass/D");
  theTree->Branch("genB0Px",   &genB0Px,   "genB0Px/D");
  theTree->Branch("genB0Py",   &genB0Py,   "genB0Py/D");
  theTree->Branch("genB0Pz",   &genB0Pz,   "genB0Pz/D");

  // ### Generated B0 Vtx ###
  theTree->Branch("genB0VtxX", &genB0VtxX, "genB0VtxX/D");
  theTree->Branch("genB0VtxY", &genB0VtxY, "genB0VtxY/D");
  theTree->Branch("genB0VtxZ", &genB0VtxZ, "genB0VtxZ/D");

  // ### Generated K*0 Mass ###
  theTree->Branch("genKstMass", &genKstMass, "genKstMass/D");
  theTree->Branch("genKstPx",   &genKstPx,   "genKstPx/D");
  theTree->Branch("genKstPy",   &genKstPy,   "genKstPy/D");
  theTree->Branch("genKstPz",   &genKstPz,   "genKstPz/D");

  // ### Generated K*0 Vtx ###
  theTree->Branch("genKstVtxX", &genKstVtxX, "genKstVtxX/D");
  theTree->Branch("genKstVtxY", &genKstVtxY, "genKstVtxY/D");
  theTree->Branch("genKstVtxZ", &genKstVtxZ, "genKstVtxZ/D");

  // ### Generated J/psi or psi(2S) Mass and Vtx ###
  theTree->Branch("genPsiMass", &genPsiMass, "genPsiMass/D");
  theTree->Branch("genPsiVtxX", &genPsiVtxX, "genPsiVtxX/D");
  theTree->Branch("genPsiVtxY", &genPsiVtxY, "genPsiVtxY/D");
  theTree->Branch("genPsiVtxZ", &genPsiVtxZ, "genPsiVtxZ/D");

  // ### Generated mu- ###
  theTree->Branch("genMumMother", &genMumMother, "genMumMother/I");
  theTree->Branch("genMumPx",     &genMumPx,     "genMumPx/D");
  theTree->Branch("genMumPy",     &genMumPy,     "genMumPy/D");
  theTree->Branch("genMumPz",     &genMumPz,     "genMumPz/D");

  // ### Generated mu+ ###
  theTree->Branch("genMupMother", &genMupMother, "genMupMother/I");
  theTree->Branch("genMupPx",     &genMupPx,     "genMupPx/D");
  theTree->Branch("genMupPy",     &genMupPy,     "genMupPy/D");
  theTree->Branch("genMupPz",     &genMupPz,     "genMupPz/D");

  // ### Generated K*0 track- ###
  theTree->Branch("genKstTrkmMother", &genKstTrkmMother, "genKstTrkmMother/I");
  theTree->Branch("genKstTrkmID",     &genKstTrkmID,     "genKstTrkmID/I");
  theTree->Branch("genKstTrkmPx",     &genKstTrkmPx,     "genKstTrkmPx/D");
  theTree->Branch("genKstTrkmPy",     &genKstTrkmPy,     "genKstTrkmPy/D");
  theTree->Branch("genKstTrkmPz",     &genKstTrkmPz,     "genKstTrkmPz/D");

  // ### Generated K*0 track+ ###
  theTree->Branch("genKstTrkpMother", &genKstTrkpMother, "genKstTrkpMother/I");
  theTree->Branch("genKstTrkpID",     &genKstTrkpID,     "genKstTrkpID/I");
  theTree->Branch("genKstTrkpPx",     &genKstTrkpPx,     "genKstTrkpPx/D");
  theTree->Branch("genKstTrkpPy",     &genKstTrkpPy,     "genKstTrkpPy/D");
  theTree->Branch("genKstTrkpPz",     &genKstTrkpPz,     "genKstTrkpPz/D");

  // ### Matching Between Reconstructed and Generated ###
  theTree->Branch("truthMatchSignal", &truthMatchSignal);
  theTree->Branch("truthMatchMum",    &truthMatchMum);
  theTree->Branch("truthMatchMup",    &truthMatchMup);
  theTree->Branch("truthMatchTrkm",   &truthMatchTrkm);
  theTree->Branch("truthMatchTrkp",   &truthMatchTrkp);
}

void B0KstMuMuTreeContent::SetBranchAddresses (TTree* theTree)
{
  // ########################################################
  // # Run Number, event number, #reco vtx and event weight #
  // ########################################################
  theTree->SetBranchAddress("runN",            &runN);
  theTree->SetBranchAddress("eventN",          &eventN);
  theTree->SetBranchAddress("recoVtxN",        &recoVtxN);
  theTree->SetBranchAddress("evWeight",        &evWeight);
  theTree->SetBranchAddress("evWeightE2",      &evWeightE2);
  theTree->SetBranchAddress("numEventsTried",  &numEventsTried);
  theTree->SetBranchAddress("numEventsPassed", &numEventsPassed);

  // ### Trigger ###
  theTree->SetBranchAddress("TrigTable",     &TrigTable);
  theTree->SetBranchAddress("TrigPrescales", &TrigPrescales);

  theTree->SetBranchAddress("nB", &nB);

  // ### Primary Vertex and Beam Spot ###
  theTree->SetBranchAddress("priVtxCL", &priVtxCL);
  theTree->SetBranchAddress("priVtxX",  &priVtxX);
  theTree->SetBranchAddress("priVtxY",  &priVtxY);
  theTree->SetBranchAddress("priVtxZ",  &priVtxZ);
  theTree->SetBranchAddress("bsX",      &bsX);
  theTree->SetBranchAddress("bsY",      &bsY);

  // ### B0 Mass ###  
  theTree->SetBranchAddress("bMass",     &bMass);
  theTree->SetBranchAddress("bMassE",    &bMassE);
  theTree->SetBranchAddress("bBarMass",  &bBarMass);
  theTree->SetBranchAddress("bBarMassE", &bBarMassE);
  theTree->SetBranchAddress("bPx",       &bPx);
  theTree->SetBranchAddress("bPy",       &bPy);
  theTree->SetBranchAddress("bPz",       &bPz);

  // ### Pileup information in MC ###
  theTree->SetBranchAddress("bunchXingMC",           &bunchXingMC);
  theTree->SetBranchAddress("numInteractionsMC",     &numInteractionsMC);
  theTree->SetBranchAddress("trueNumInteractionsMC", &trueNumInteractionsMC);

  // ### B0 Vtx ###
  theTree->SetBranchAddress("bVtxCL",        &bVtxCL);
  theTree->SetBranchAddress("bVtxX",         &bVtxX);
  theTree->SetBranchAddress("bVtxY",         &bVtxY);
  theTree->SetBranchAddress("bVtxZ",         &bVtxZ);
  theTree->SetBranchAddress("bCosAlphaVtx",  &bCosAlphaVtx);
  theTree->SetBranchAddress("bCosAlphaVtxE", &bCosAlphaVtxE);
  theTree->SetBranchAddress("bCosAlphaBS",   &bCosAlphaBS);
  theTree->SetBranchAddress("bCosAlphaBSE",  &bCosAlphaBSE);
  theTree->SetBranchAddress("bLVtx",         &bLVtx);
  theTree->SetBranchAddress("bLVtxE",        &bLVtxE);
  theTree->SetBranchAddress("bLBS",          &bLBS);
  theTree->SetBranchAddress("bLBSE",         &bLBSE);
  theTree->SetBranchAddress("bDCAVtx",       &bDCAVtx);
  theTree->SetBranchAddress("bDCAVtxE",      &bDCAVtxE);
  theTree->SetBranchAddress("bDCABS",        &bDCABS);
  theTree->SetBranchAddress("bDCABSE",       &bDCABSE);

  // ### B0 ctau ###
  theTree->SetBranchAddress("bctauPVBS",  &bctauPVBS);
  theTree->SetBranchAddress("bctauPVBSE", &bctauPVBSE);

  // ### K*0 Mass ###
  theTree->SetBranchAddress("kstMass",     &kstMass);
  theTree->SetBranchAddress("kstMassE",    &kstMassE);
  theTree->SetBranchAddress("kstBarMass",  &kstBarMass);
  theTree->SetBranchAddress("kstBarMassE", &kstBarMassE);
  theTree->SetBranchAddress("kstPx",       &kstPx);
  theTree->SetBranchAddress("kstPy",       &kstPy);
  theTree->SetBranchAddress("kstPz",       &kstPz);

  // ### K*0 Vtx ###
  theTree->SetBranchAddress("kstVtxCL", &kstVtxCL);
  theTree->SetBranchAddress("kstVtxX",  &kstVtxX);
  theTree->SetBranchAddress("kstVtxY",  &kstVtxY);
  theTree->SetBranchAddress("kstVtxZ",  &kstVtxZ);

  // ### mu+ mu- Mass ###
  theTree->SetBranchAddress("mumuMass",  &mumuMass);
  theTree->SetBranchAddress("mumuMassE", &mumuMassE);
  theTree->SetBranchAddress("mumuPx",    &mumuPx);
  theTree->SetBranchAddress("mumuPy",    &mumuPy);
  theTree->SetBranchAddress("mumuPz",    &mumuPz);

  // ### mu+ mu- Vtx ###
  theTree->SetBranchAddress("mumuVtxCL",       &mumuVtxCL);
  theTree->SetBranchAddress("mumuVtxX",        &mumuVtxX);
  theTree->SetBranchAddress("mumuVtxY",        &mumuVtxY);
  theTree->SetBranchAddress("mumuVtxZ",        &mumuVtxZ);
  theTree->SetBranchAddress("mumuCosAlphaBS",  &mumuCosAlphaBS);
  theTree->SetBranchAddress("mumuCosAlphaBSE", &mumuCosAlphaBSE);
  theTree->SetBranchAddress("mumuLBS",         &mumuLBS);
  theTree->SetBranchAddress("mumuLBSE",        &mumuLBSE);
  theTree->SetBranchAddress("mumuDCA",         &mumuDCA);

  // ### mu- ###  
  theTree->SetBranchAddress("mumHighPurity",    &mumHighPurity);
  theTree->SetBranchAddress("mumCL",            &mumCL);
  theTree->SetBranchAddress("mumNormChi2",      &mumNormChi2);
  theTree->SetBranchAddress("mumPx",            &mumPx);
  theTree->SetBranchAddress("mumPy",            &mumPy);
  theTree->SetBranchAddress("mumPz",            &mumPz);
  theTree->SetBranchAddress("mumDCAVtx",        &mumDCAVtx);
  theTree->SetBranchAddress("mumDCAVtxE",       &mumDCAVtxE);
  theTree->SetBranchAddress("mumDCABS",         &mumDCABS);
  theTree->SetBranchAddress("mumDCABSE",        &mumDCABSE);
  theTree->SetBranchAddress("mumKinkChi2",      &mumKinkChi2);
  theTree->SetBranchAddress("mumFracHits",      &mumFracHits);
  theTree->SetBranchAddress("mumdxyVtx",        &mumdxyVtx);
  theTree->SetBranchAddress("mumdzVtx",         &mumdzVtx);
  theTree->SetBranchAddress("mumDeltaRwithMC",  &mumDeltaRwithMC);
  theTree->SetBranchAddress("mumCat",           &mumCat);
  theTree->SetBranchAddress("mumNPixHits",      &mumNPixHits);
  theTree->SetBranchAddress("mumNPixLayers",    &mumNPixLayers);
  theTree->SetBranchAddress("mumNTrkHits",      &mumNTrkHits);
  theTree->SetBranchAddress("mumNTrkLayers",    &mumNTrkLayers);
  theTree->SetBranchAddress("mumNMuonHits",     &mumNMuonHits);
  theTree->SetBranchAddress("mumNMatchStation", &mumNMatchStation);
  theTree->SetBranchAddress("mumTrig",          &mumTrig);

  // ### mu+ ###  
  theTree->SetBranchAddress("mupHighPurity",    &mupHighPurity);
  theTree->SetBranchAddress("mupCL",            &mupCL);
  theTree->SetBranchAddress("mupNormChi2",      &mupNormChi2);
  theTree->SetBranchAddress("mupPx",            &mupPx);
  theTree->SetBranchAddress("mupPy",            &mupPy);
  theTree->SetBranchAddress("mupPz",            &mupPz);
  theTree->SetBranchAddress("mupDCAVtx",        &mupDCAVtx);
  theTree->SetBranchAddress("mupDCAVtxE",       &mupDCAVtxE);
  theTree->SetBranchAddress("mupDCABS",         &mupDCABS);
  theTree->SetBranchAddress("mupDCABSE",        &mupDCABSE);
  theTree->SetBranchAddress("mupKinkChi2",      &mupKinkChi2);
  theTree->SetBranchAddress("mupFracHits",      &mupFracHits);
  theTree->SetBranchAddress("mupdxyVtx",        &mupdxyVtx);
  theTree->SetBranchAddress("mupdzVtx",         &mupdzVtx);
  theTree->SetBranchAddress("mupDeltaRwithMC",  &mupDeltaRwithMC);
  theTree->SetBranchAddress("mupCat",           &mupCat);
  theTree->SetBranchAddress("mupNPixHits",      &mupNPixHits);
  theTree->SetBranchAddress("mupNPixLayers",    &mupNPixLayers);
  theTree->SetBranchAddress("mupNTrkHits",      &mupNTrkHits);
  theTree->SetBranchAddress("mupNTrkLayers",    &mupNTrkLayers);
  theTree->SetBranchAddress("mupNMuonHits",     &mupNMuonHits);
  theTree->SetBranchAddress("mupNMatchStation", &mupNMatchStation);
  theTree->SetBranchAddress("mupTrig",          &mupTrig);

  // ### K*0 track- ###
  theTree->SetBranchAddress("kstTrkmHighPurity",   &kstTrkmHighPurity);
  theTree->SetBranchAddress("kstTrkmCL",           &kstTrkmCL);
  theTree->SetBranchAddress("kstTrkmNormChi2",     &kstTrkmNormChi2);
  theTree->SetBranchAddress("kstTrkmPx",           &kstTrkmPx);
  theTree->SetBranchAddress("kstTrkmPy",           &kstTrkmPy);
  theTree->SetBranchAddress("kstTrkmPz",           &kstTrkmPz);
  theTree->SetBranchAddress("kstTrkmDCAVtx",       &kstTrkmDCAVtx);
  theTree->SetBranchAddress("kstTrkmDCAVtxE",      &kstTrkmDCAVtxE);
  theTree->SetBranchAddress("kstTrkmDCABS",        &kstTrkmDCABS);
  theTree->SetBranchAddress("kstTrkmDCABSE",       &kstTrkmDCABSE);
  theTree->SetBranchAddress("kstTrkmFracHits",     &kstTrkmFracHits);
  theTree->SetBranchAddress("kstTrkmdxyVtx",       &kstTrkmdxyVtx);
  theTree->SetBranchAddress("kstTrkmdzVtx",        &kstTrkmdzVtx);
  theTree->SetBranchAddress("kstTrkmDeltaRwithMC", &kstTrkmDeltaRwithMC);
  theTree->SetBranchAddress("kstTrkmNPixHits",     &kstTrkmNPixHits);
  theTree->SetBranchAddress("kstTrkmNPixLayers",   &kstTrkmNPixLayers);
  theTree->SetBranchAddress("kstTrkmNTrkHits",     &kstTrkmNTrkHits);
  theTree->SetBranchAddress("kstTrkmNTrkLayers",   &kstTrkmNTrkLayers);
  theTree->SetBranchAddress("kstTrkmMuMatch",      &kstTrkmMuMatch);

  // ### K*0 track+ ### 
  theTree->SetBranchAddress("kstTrkpHighPurity",   &kstTrkpHighPurity);
  theTree->SetBranchAddress("kstTrkpCL",           &kstTrkpCL);
  theTree->SetBranchAddress("kstTrkpNormChi2",     &kstTrkpNormChi2);
  theTree->SetBranchAddress("kstTrkpPx",           &kstTrkpPx);
  theTree->SetBranchAddress("kstTrkpPy",           &kstTrkpPy);
  theTree->SetBranchAddress("kstTrkpPz",           &kstTrkpPz);
  theTree->SetBranchAddress("kstTrkpDCAVtx",       &kstTrkpDCAVtx);
  theTree->SetBranchAddress("kstTrkpDCAVtxE",      &kstTrkpDCAVtxE);
  theTree->SetBranchAddress("kstTrkpDCABS",        &kstTrkpDCABS);
  theTree->SetBranchAddress("kstTrkpDCABSE",       &kstTrkpDCABSE);
  theTree->SetBranchAddress("kstTrkpFracHits",     &kstTrkpFracHits);
  theTree->SetBranchAddress("kstTrkpdxyVtx",       &kstTrkpdxyVtx);
  theTree->SetBranchAddress("kstTrkpdzVtx",        &kstTrkpdzVtx);
  theTree->SetBranchAddress("kstTrkpDeltaRwithMC", &kstTrkpDeltaRwithMC);
  theTree->SetBranchAddress("kstTrkpNPixHits",     &kstTrkpNPixHits);
  theTree->SetBranchAddress("kstTrkpNPixLayers",   &kstTrkpNPixLayers);
  theTree->SetBranchAddress("kstTrkpNTrkHits",     &kstTrkpNTrkHits);
  theTree->SetBranchAddress("kstTrkpNTrkLayers",   &kstTrkpNTrkLayers);
  theTree->SetBranchAddress("kstTrkpMuMatch",      &kstTrkpMuMatch);

  // ### Generated Observables ###
  theTree->SetBranchAddress("genSignal",        &genSignal);
  theTree->SetBranchAddress("genMuMuBG",        &genMuMuBG);
  theTree->SetBranchAddress("genMuMuBGnTrksm",  &genMuMuBGnTrksm);
  theTree->SetBranchAddress("genMuMuBGnTrksp",  &genMuMuBGnTrksp);
  theTree->SetBranchAddress("genPsiPrompt",     &genPsiPrompt);
  theTree->SetBranchAddress("genSignHasFSR",    &genSignHasFSR);
  theTree->SetBranchAddress("genSignKstHasFSR", &genSignKstHasFSR);
  theTree->SetBranchAddress("genSignPsiHasFSR", &genSignPsiHasFSR);

  // ### Generated Primary Vertex ###
  theTree->SetBranchAddress("genPriVtxX", &genPriVtxX);
  theTree->SetBranchAddress("genPriVtxY", &genPriVtxY);
  theTree->SetBranchAddress("genPriVtxZ", &genPriVtxZ);

  // ### Generated B0 Mass ###
  theTree->SetBranchAddress("genB0Mass", &genB0Mass);
  theTree->SetBranchAddress("genB0Px",   &genB0Px);
  theTree->SetBranchAddress("genB0Py",   &genB0Py);
  theTree->SetBranchAddress("genB0Pz",   &genB0Pz);

  // ### Generated B0 Vtx ###
  theTree->SetBranchAddress("genB0VtxX", &genB0VtxX);
  theTree->SetBranchAddress("genB0VtxY", &genB0VtxY);
  theTree->SetBranchAddress("genB0VtxZ", &genB0VtxZ);

  // ### Generated K*0 Mass ###
  theTree->SetBranchAddress("genKstMass", &genKstMass);
  theTree->SetBranchAddress("genKstPx",   &genKstPx);
  theTree->SetBranchAddress("genKstPy",   &genKstPy);
  theTree->SetBranchAddress("genKstPz",   &genKstPz);

  // ### Generated K*0 Vtx ###
  theTree->SetBranchAddress("genKstVtxX", &genKstVtxX);
  theTree->SetBranchAddress("genKstVtxY", &genKstVtxY);
  theTree->SetBranchAddress("genKstVtxZ", &genKstVtxZ);

  // ### Generated J/psi or psi(2S) Mass and Vtx ###
  theTree->SetBranchAddress("genPsiMass", &genPsiMass);
  theTree->SetBranchAddress("genPsiVtxX", &genPsiVtxX);
  theTree->SetBranchAddress("genPsiVtxY", &genPsiVtxY);
  theTree->SetBranchAddress("genPsiVtxZ", &genPsiVtxZ);

  // ### Generated mu- ###
  theTree->SetBranchAddress("genMumMother", &genMumMother);
  theTree->SetBranchAddress("genMumPx",     &genMumPx);
  theTree->SetBranchAddress("genMumPy",     &genMumPy);
  theTree->SetBranchAddress("genMumPz",     &genMumPz);

  // ### Generated mu+ ###
  theTree->SetBranchAddress("genMupMother", &genMupMother);
  theTree->SetBranchAddress("genMupPx",     &genMupPx);
  theTree->SetBranchAddress("genMupPy",     &genMupPy);
  theTree->SetBranchAddress("genMupPz",     &genMupPz);

  // ### Generated K*0 track- ###
  theTree->SetBranchAddress("genKstTrkmMother", &genKstTrkmMother);
  theTree->SetBranchAddress("genKstTrkmID",     &genKstTrkmID);
  theTree->SetBranchAddress("genKstTrkmPx",     &genKstTrkmPx);
  theTree->SetBranchAddress("genKstTrkmPy",     &genKstTrkmPy);
  theTree->SetBranchAddress("genKstTrkmPz",     &genKstTrkmPz);

  // ### Generated K*0 track+ ###
  theTree->SetBranchAddress("genKstTrkpMother", &genKstTrkpMother);
  theTree->SetBranchAddress("genKstTrkpID",     &genKstTrkpID);
  theTree->SetBranchAddress("genKstTrkpPx",     &genKstTrkpPx);
  theTree->SetBranchAddress("genKstTrkpPy",     &genKstTrkpPy);
  theTree->SetBranchAddress("genKstTrkpPz",     &genKstTrkpPz);

  // ### Matching Between Reconstructed and Generated ###
  theTree->SetBranchAddress("truthMatchSignal", &truthMatchSignal);
  theTree->SetBranchAddress("truthMatchMum",    &truthMatchMum);
  theTree->SetBranchAddress("truthMatchMup",    &truthMatchMup);
  theTree->SetBranchAddress("truthMatchTrkm",   &truthMatchTrkm);
  theTree->SetBranchAddress("truthMatchTrkp",   &truthMatchTrkp);
}

void B0KstMuMuTreeContent::CopyCandidate (B0KstMuMuTreeContent* NTupleIn, int index)
{
  CopyScalars(NTupleIn);
  CopyVectors(NTupleIn,index);
}

void B0KstMuMuTreeContent::CopyAllCandidates (B0KstMuMuTreeContent* NTupleIn)
{
  CopyScalars(NTupleIn);
  if (NTupleIn->bMass != NULL)
    for (unsigned int index = 0; index < NTupleIn->bMass->size(); index++) CopyVectors(NTupleIn,index);
}

void B0KstMuMuTreeContent::CopyScalars (B0KstMuMuTreeContent* NTupleIn)
{
  // ########################################################
  // # Run Number, event number, #reco vtx and event weight #
  // ########################################################
  runN            = NTupleIn->runN;
  eventN          = NTupleIn->eventN;
  recoVtxN        = NTupleIn->recoVtxN;
  evWeight        = NTupleIn->evWeight;
  evWeightE2      = NTupleIn->evWeightE2;
  numEventsTried  = NTupleIn->numEventsTried;
  numEventsPassed = NTupleIn->numEventsPassed;

  // ### Trigger ###
  for (unsigned int i = 0; i < NTupleIn->TrigTable->size(); i++)
    {
      TrigTable->push_back(NTupleIn->TrigTable->at(i));
      TrigPrescales->push_back(NTupleIn->TrigPrescales->at(i));
    }
  
  nB = NTupleIn->nB;
  
  // ### Primary Vertex and Beam Spot ###
  priVtxCL = NTupleIn->priVtxCL;
  priVtxX  = NTupleIn->priVtxX;
  priVtxY  = NTupleIn->priVtxY;
  priVtxZ  = NTupleIn->priVtxZ;
  bsX      = NTupleIn->bsX;
  bsY      = NTupleIn->bsY;

  // ### Pileup information in MC ###
  for (unsigned int i = 0; i < NTupleIn->bunchXingMC->size(); i++)
    {
      bunchXingMC->push_back(NTupleIn->bunchXingMC->at(i));
      numInteractionsMC->push_back(NTupleIn->numInteractionsMC->at(i));
      trueNumInteractionsMC->push_back(NTupleIn->trueNumInteractionsMC->at(i));
    }

  // ### Generated Observables ###
  genSignal        = NTupleIn->genSignal;
  genMuMuBG        = NTupleIn->genMuMuBG;
  genMuMuBGnTrksm  = NTupleIn->genMuMuBGnTrksm;
  genMuMuBGnTrksp  = NTupleIn->genMuMuBGnTrksp;
  genPsiPrompt     = NTupleIn->genPsiPrompt;
  genSignHasFSR    = NTupleIn->genSignHasFSR;
  genSignKstHasFSR = NTupleIn->genSignKstHasFSR;
  genSignPsiHasFSR = NTupleIn->genSignPsiHasFSR;

  // ### Generated Primary Vertex ###
  genPriVtxX = NTupleIn->genPriVtxX;
  genPriVtxY = NTupleIn->genPriVtxY;
  genPriVtxZ = NTupleIn->genPriVtxZ;

  // ### Generated B0 Mass ####
  genB0Mass = NTupleIn->genB0Mass;
  genB0Px   = NTupleIn->genB0Px;
  genB0Py   = NTupleIn->genB0Py;
  genB0Pz   = NTupleIn->genB0Pz;

  // ### Generated B0 Vtx ####
  genB0VtxX = NTupleIn->genB0VtxX;
  genB0VtxY = NTupleIn->genB0VtxY;
  genB0VtxZ = NTupleIn->genB0VtxZ;

  // ### Generated K*0 Mass ###
  genKstMass = NTupleIn->genKstMass;
  genKstPx   = NTupleIn->genKstPx;
  genKstPy   = NTupleIn->genKstPy;
  genKstPz   = NTupleIn->genKstPz;

  // ### Generated K*0 Vtx ###
  genKstVtxX = NTupleIn->genKstVtxX;
  genKstVtxY = NTupleIn->genKstVtxY;
  genKstVtxZ = NTupleIn->genKstVtxZ;

  // ### Generated J/psi or psi(2S) Mass and Vtx ###
  genPsiMass = NTupleIn->genPsiMass;
  genPsiVtxX = NTupleIn->genPsiVtxX;
  genPsiVtxY = NTupleIn->genPsiVtxY;
  genPsiVtxZ = NTupleIn->genPsiVtxZ;

  // ### Generated mu- ###
  genMumMother = NTupleIn->genMumMother;
  genMumPx     = NTupleIn->genMumPx;
  genMumPy     = NTupleIn->genMumPy;
  genMumPz     = NTupleIn->genMumPz;

  // ### Generated mu+ ###
  genMupMother = NTupleIn->genMupMother;
  genMupPx     = NTupleIn->genMupPx;
  genMupPy     = NTupleIn->genMupPy;
  genMupPz     = NTupleIn->genMupPz;

  // ### Generated K*0 track- ###
  genKstTrkmMother = NTupleIn->genKstTrkmMother;
  genKstTrkmID     = NTupleIn->genKstTrkmID;
  genKstTrkmPx     = NTupleIn->genKstTrkmPx;
  genKstTrkmPy     = NTupleIn->genKstTrkmPy;
  genKstTrkmPz     = NTupleIn->genKstTrkmPz;

  // ### Generated K*0 track+ ###
  genKstTrkpMother = NTupleIn->genKstTrkpMother;
  genKstTrkpID     = NTupleIn->genKstTrkpID;
  genKstTrkpPx     = NTupleIn->genKstTrkpPx;
  genKstTrkpPy     = NTupleIn->genKstTrkpPy;
  genKstTrkpPz     = NTupleIn->genKstTrkpPz;
}

void B0KstMuMuTreeContent::CopyVectors (B0KstMuMuTreeContent* NTupleIn, int index)
{
  if ((index >= 0) && (NTupleIn->bMass != NULL) && (NTupleIn->bMass->size() != 0))
    {
      // ### B0 Mass ###
      bMass->push_back(NTupleIn->bMass->at(index));
      bMassE->push_back(NTupleIn->bMassE->at(index));
      bBarMass->push_back(NTupleIn->bBarMass->at(index));
      bBarMassE->push_back(NTupleIn->bBarMassE->at(index));
      bPx->push_back(NTupleIn->bPx->at(index));
      bPy->push_back(NTupleIn->bPy->at(index));
      bPz->push_back(NTupleIn->bPz->at(index));

      // ### B0 Vtx ###
      bVtxCL->push_back(NTupleIn->bVtxCL->at(index));
      bVtxX->push_back(NTupleIn->bVtxX->at(index));
      bVtxY->push_back(NTupleIn->bVtxY->at(index));
      bVtxZ->push_back(NTupleIn->bVtxZ->at(index));
      bCosAlphaVtx->push_back(NTupleIn->bCosAlphaVtx->at(index)); 
      bCosAlphaVtxE->push_back(NTupleIn->bCosAlphaVtxE->at(index)); 
      bCosAlphaBS->push_back(NTupleIn->bCosAlphaBS->at(index)); 
      bCosAlphaBSE->push_back(NTupleIn->bCosAlphaBSE->at(index));
      bLVtx->push_back(NTupleIn->bLVtx->at(index));
      bLVtxE->push_back(NTupleIn->bLVtxE->at(index));
      bLBS->push_back(NTupleIn->bLBS->at(index));
      bLBSE->push_back(NTupleIn->bLBSE->at(index));
      bDCAVtx->push_back(NTupleIn->bDCAVtx->at(index));
      bDCAVtxE->push_back(NTupleIn->bDCAVtxE->at(index));
      bDCABS->push_back(NTupleIn->bDCABS->at(index));
      bDCABSE->push_back(NTupleIn->bDCABSE->at(index));

      // ### B0 ctau ###
      bctauPVBS->push_back(NTupleIn->bctauPVBS->at(index));
      bctauPVBSE->push_back(NTupleIn->bctauPVBSE->at(index));

      // ### K*0 Mass ###
      kstMass->push_back(NTupleIn->kstMass->at(index));
      kstMassE->push_back(NTupleIn->kstMassE->at(index));
      kstBarMass->push_back(NTupleIn->kstBarMass->at(index));
      kstBarMassE->push_back(NTupleIn->kstBarMassE->at(index));
      kstPx->push_back(NTupleIn->kstPx->at(index));
      kstPy->push_back(NTupleIn->kstPy->at(index));
      kstPz->push_back(NTupleIn->kstPz->at(index));
      
      // ### K*0 Vtx ###
      kstVtxCL->push_back(NTupleIn->kstVtxCL->at(index));
      kstVtxX->push_back(NTupleIn->kstVtxX->at(index));
      kstVtxY->push_back(NTupleIn->kstVtxY->at(index));
      kstVtxZ->push_back(NTupleIn->kstVtxZ->at(index));

      // ### mu+ mu- Mass ###
      mumuMass->push_back(NTupleIn->mumuMass->at(index));
      mumuMassE->push_back(NTupleIn->mumuMassE->at(index));
      mumuPx->push_back(NTupleIn->mumuPx->at(index));
      mumuPy->push_back(NTupleIn->mumuPy->at(index));
      mumuPz->push_back(NTupleIn->mumuPz->at(index));
      
      // ### mu+ mu- Vtx ###
      mumuVtxCL->push_back(NTupleIn->mumuVtxCL->at(index));
      mumuVtxX->push_back(NTupleIn->mumuVtxX->at(index));
      mumuVtxY->push_back(NTupleIn->mumuVtxY->at(index));
      mumuVtxZ->push_back(NTupleIn->mumuVtxZ->at(index));
      mumuCosAlphaBS->push_back(NTupleIn->mumuCosAlphaBS->at(index));
      mumuCosAlphaBSE->push_back(NTupleIn->mumuCosAlphaBSE->at(index));
      mumuLBS->push_back(NTupleIn->mumuLBS->at(index));
      mumuLBSE->push_back(NTupleIn->mumuLBSE->at(index));
      mumuDCA->push_back(NTupleIn->mumuDCA->at(index));

      // ### mu- ###  
      mumHighPurity->push_back(NTupleIn->mumHighPurity->at(index));
      mumCL->push_back(NTupleIn->mumCL->at(index));
      mumNormChi2->push_back(NTupleIn->mumNormChi2->at(index));
      mumPx->push_back(NTupleIn->mumPx->at(index));
      mumPy->push_back(NTupleIn->mumPy->at(index));
      mumPz->push_back(NTupleIn->mumPz->at(index));
      mumDCAVtx->push_back(NTupleIn->mumDCAVtx->at(index));
      mumDCAVtxE->push_back(NTupleIn->mumDCAVtxE->at(index));
      mumDCABS->push_back(NTupleIn->mumDCABS->at(index));
      mumDCABSE->push_back(NTupleIn->mumDCABSE->at(index));
      mumKinkChi2->push_back(NTupleIn->mumKinkChi2->at(index));
      mumFracHits->push_back(NTupleIn->mumFracHits->at(index));
      mumdxyVtx->push_back(NTupleIn->mumdxyVtx->at(index));
      mumdzVtx->push_back(NTupleIn->mumdzVtx->at(index));
      mumDeltaRwithMC->push_back(NTupleIn->mumDeltaRwithMC->at(index));
      mumCat->push_back(NTupleIn->mumCat->at(index));
      mumNPixHits->push_back(NTupleIn->mumNPixHits->at(index));
      mumNPixLayers->push_back(NTupleIn->mumNPixLayers->at(index));
      mumNTrkHits->push_back(NTupleIn->mumNTrkHits->at(index));
      mumNTrkLayers->push_back(NTupleIn->mumNTrkLayers->at(index));
      mumNMuonHits->push_back(NTupleIn->mumNMuonHits->at(index));
      mumNMatchStation->push_back(NTupleIn->mumNMatchStation->at(index));
      mumTrig->push_back(NTupleIn->mumTrig->at(index));

      // ### mu+ ###  
      mupHighPurity->push_back(NTupleIn->mupHighPurity->at(index));
      mupCL->push_back(NTupleIn->mupCL->at(index));
      mupNormChi2->push_back(NTupleIn->mupNormChi2->at(index));
      mupPx->push_back(NTupleIn->mupPx->at(index));
      mupPy->push_back(NTupleIn->mupPy->at(index));
      mupPz->push_back(NTupleIn->mupPz->at(index));
      mupDCAVtx->push_back(NTupleIn->mupDCAVtx->at(index));
      mupDCAVtxE->push_back(NTupleIn->mupDCAVtxE->at(index));
      mupDCABS->push_back(NTupleIn->mupDCABS->at(index));
      mupDCABSE->push_back(NTupleIn->mupDCABSE->at(index));
      mupKinkChi2->push_back(NTupleIn->mupKinkChi2->at(index));
      mupFracHits->push_back(NTupleIn->mupFracHits->at(index));
      mupdxyVtx->push_back(NTupleIn->mupdxyVtx->at(index));
      mupdzVtx->push_back(NTupleIn->mupdzVtx->at(index));
      mupDeltaRwithMC->push_back(NTupleIn->mupDeltaRwithMC->at(index));
      mupCat->push_back(NTupleIn->mupCat->at(index));
      mupNPixHits->push_back(NTupleIn->mupNPixHits->at(index));
      mupNPixLayers->push_back(NTupleIn->mupNPixLayers->at(index));
      mupNTrkHits->push_back(NTupleIn->mupNTrkHits->at(index));
      mupNTrkLayers->push_back(NTupleIn->mupNTrkLayers->at(index));
      mupNMuonHits->push_back(NTupleIn->mupNMuonHits->at(index));
      mupNMatchStation->push_back(NTupleIn->mupNMatchStation->at(index));
      mupTrig->push_back(NTupleIn->mupTrig->at(index));

      // ### K*0 track- ###
      kstTrkmHighPurity->push_back(NTupleIn->kstTrkmHighPurity->at(index));
      kstTrkmCL->push_back(NTupleIn->kstTrkmCL->at(index));
      kstTrkmNormChi2->push_back(NTupleIn->kstTrkmNormChi2->at(index));
      kstTrkmPx->push_back(NTupleIn->kstTrkmPx->at(index));
      kstTrkmPy->push_back(NTupleIn->kstTrkmPy->at(index));
      kstTrkmPz->push_back(NTupleIn->kstTrkmPz->at(index));
      kstTrkmDCAVtx->push_back(NTupleIn->kstTrkmDCAVtx->at(index));
      kstTrkmDCAVtxE->push_back(NTupleIn->kstTrkmDCAVtxE->at(index));
      kstTrkmDCABS->push_back(NTupleIn->kstTrkmDCABS->at(index));
      kstTrkmDCABSE->push_back(NTupleIn->kstTrkmDCABSE->at(index));
      kstTrkmFracHits->push_back(NTupleIn->kstTrkmFracHits->at(index));
      kstTrkmdxyVtx->push_back(NTupleIn->kstTrkmdxyVtx->at(index));
      kstTrkmdzVtx->push_back(NTupleIn->kstTrkmdzVtx->at(index));
      kstTrkmDeltaRwithMC->push_back(NTupleIn->kstTrkmDeltaRwithMC->at(index));
      kstTrkmNPixHits->push_back(NTupleIn->kstTrkmNPixHits->at(index));
      kstTrkmNPixLayers->push_back(NTupleIn->kstTrkmNPixLayers->at(index));
      kstTrkmNTrkHits->push_back(NTupleIn->kstTrkmNTrkHits->at(index));
      kstTrkmNTrkLayers->push_back(NTupleIn->kstTrkmNTrkLayers->at(index));
      kstTrkmMuMatch->push_back(NTupleIn->kstTrkmMuMatch->at(index));
 
      // ### K*0 track+ ###
      kstTrkpHighPurity->push_back(NTupleIn->kstTrkpHighPurity->at(index));
      kstTrkpCL->push_back(NTupleIn->kstTrkpCL->at(index));
      kstTrkpNormChi2->push_back(NTupleIn->kstTrkpNormChi2->at(index));
      kstTrkpPx->push_back(NTupleIn->kstTrkpPx->at(index));
      kstTrkpPy->push_back(NTupleIn->kstTrkpPy->at(index));
      kstTrkpPz->push_back(NTupleIn->kstTrkpPz->at(index));
      kstTrkpDCAVtx->push_back(NTupleIn->kstTrkpDCAVtx->at(index));
      kstTrkpDCAVtxE->push_back(NTupleIn->kstTrkpDCAVtxE->at(index));
      kstTrkpDCABS->push_back(NTupleIn->kstTrkpDCABS->at(index));
      kstTrkpDCABSE->push_back(NTupleIn->kstTrkpDCABSE->at(index));
      kstTrkpFracHits->push_back(NTupleIn->kstTrkpFracHits->at(index));
      kstTrkpdxyVtx->push_back(NTupleIn->kstTrkpdxyVtx->at(index));
      kstTrkpdzVtx->push_back(NTupleIn->kstTrkpdzVtx->at(index));
      kstTrkpDeltaRwithMC->push_back(NTupleIn->kstTrkpDeltaRwithMC->at(index));
      kstTrkpNPixHits->push_back(NTupleIn->kstTrkpNPixHits->at(index));
      kstTrkpNPixLayers->push_back(NTupleIn->kstTrkpNPixLayers->at(index));
      kstTrkpNTrkHits->push_back(NTupleIn->kstTrkpNTrkHits->at(index));
      kstTrkpNTrkLayers->push_back(NTupleIn->kstTrkpNTrkLayers->at(index));
      kstTrkpMuMatch->push_back(NTupleIn->kstTrkpMuMatch->at(index));
  
      // ### Matching Between Reconstructed and Generated ###
      truthMatchSignal->push_back(NTupleIn->truthMatchSignal->at(index));
      truthMatchMum->push_back(NTupleIn->truthMatchMum->at(index));
      truthMatchMup->push_back(NTupleIn->truthMatchMup->at(index));
      truthMatchTrkm->push_back(NTupleIn->truthMatchTrkm->at(index));
      truthMatchTrkp->push_back(NTupleIn->truthMatchTrkp->at(index));
    }
}

void B0KstMuMuTreeContent::FillWithNull (unsigned int upTo)
{
  // ### B0 Mass ###
  if (bMass->size() < upTo)     for (unsigned int i = bMass->size(); i < upTo; i++)     bMass->push_back(0);
  if (bMassE->size() < upTo)    for (unsigned int i = bMassE->size(); i < upTo; i++)    bMassE->push_back(0);
  if (bBarMass->size() < upTo)  for (unsigned int i = bBarMass->size(); i < upTo; i++)  bBarMass->push_back(0);
  if (bBarMassE->size() < upTo) for (unsigned int i = bBarMassE->size(); i < upTo; i++) bBarMassE->push_back(0);
  if (bPx->size() < upTo)       for (unsigned int i = bPx->size(); i < upTo; i++)       bPx->push_back(0);
  if (bPy->size() < upTo)       for (unsigned int i = bPy->size(); i < upTo; i++)       bPy->push_back(0);
  if (bPz->size() < upTo)       for (unsigned int i = bPz->size(); i < upTo; i++)       bPz->push_back(0);

  // ### B0 Vtx ###
  if (bVtxCL->size() < upTo)        for (unsigned int i = bVtxCL->size(); i < upTo; i++)        bVtxCL->push_back(0);
  if (bVtxX->size() < upTo)         for (unsigned int i = bVtxX->size(); i < upTo; i++)         bVtxX->push_back(0);
  if (bVtxY->size() < upTo)         for (unsigned int i = bVtxY->size(); i < upTo; i++)         bVtxY->push_back(0);
  if (bVtxZ->size() < upTo)         for (unsigned int i = bVtxZ->size(); i < upTo; i++)         bVtxZ->push_back(0);
  if (bCosAlphaVtx->size() < upTo)  for (unsigned int i = bCosAlphaVtx->size(); i < upTo; i++)  bCosAlphaVtx->push_back(0);
  if (bCosAlphaVtxE->size() < upTo) for (unsigned int i = bCosAlphaVtxE->size(); i < upTo; i++) bCosAlphaVtxE->push_back(0);
  if (bCosAlphaBS->size() < upTo)   for (unsigned int i = bCosAlphaBS->size(); i < upTo; i++)   bCosAlphaBS->push_back(0);
  if (bCosAlphaBSE->size() < upTo)  for (unsigned int i = bCosAlphaBSE->size(); i < upTo; i++)  bCosAlphaBSE->push_back(0);
  if (bLVtx->size() < upTo)         for (unsigned int i = bLVtx->size(); i < upTo; i++)         bLVtx->push_back(0);
  if (bLVtxE->size() < upTo)        for (unsigned int i = bLVtxE->size(); i < upTo; i++)        bLVtxE->push_back(0);
  if (bLBS->size() < upTo)          for (unsigned int i = bLBS->size(); i < upTo; i++)          bLBS->push_back(0);
  if (bLBSE->size() < upTo)         for (unsigned int i = bLBSE->size(); i < upTo; i++)         bLBSE->push_back(0);
  if (bDCAVtx->size() < upTo)       for (unsigned int i = bDCAVtx->size(); i < upTo; i++)       bDCAVtx->push_back(0);
  if (bDCAVtxE->size() < upTo)      for (unsigned int i = bDCAVtxE->size(); i < upTo; i++)      bDCAVtxE->push_back(0);
  if (bDCABS->size() < upTo)        for (unsigned int i = bDCABS->size(); i < upTo; i++)        bDCABS->push_back(0);
  if (bDCABSE->size() < upTo)       for (unsigned int i = bDCABSE->size(); i < upTo; i++)       bDCABSE->push_back(0);

  // ### B0 ctau ###
  if (bctauPVBS->size() < upTo)  for (unsigned int i = bctauPVBS->size(); i < upTo; i++)  bctauPVBS->push_back(0);
  if (bctauPVBSE->size() < upTo) for (unsigned int i = bctauPVBSE->size(); i < upTo; i++) bctauPVBSE->push_back(0);

  // ### K*0 Mass ###
  if (kstMass->size() < upTo)     for (unsigned int i = kstMass->size(); i < upTo; i++)     kstMass->push_back(0);
  if (kstMassE->size() < upTo)    for (unsigned int i = kstMassE->size(); i < upTo; i++)    kstMassE->push_back(0);
  if (kstBarMass->size() < upTo)  for (unsigned int i = kstBarMass->size(); i < upTo; i++)  kstBarMass->push_back(0);
  if (kstBarMassE->size() < upTo) for (unsigned int i = kstBarMassE->size(); i < upTo; i++) kstBarMassE->push_back(0);
  if (kstPx->size() < upTo)       for (unsigned int i = kstPx->size(); i < upTo; i++)       kstPx->push_back(0);
  if (kstPy->size() < upTo)       for (unsigned int i = kstPy->size(); i < upTo; i++)       kstPy->push_back(0);
  if (kstPz->size() < upTo)       for (unsigned int i = kstPz->size(); i < upTo; i++)       kstPz->push_back(0);

  // ### K*0 Vtx ###
  if (kstVtxCL->size() < upTo) for (unsigned int i = kstVtxCL->size(); i < upTo; i++) kstVtxCL->push_back(0);
  if (kstVtxX->size() < upTo)  for (unsigned int i = kstVtxX->size(); i < upTo; i++)  kstVtxX->push_back(0);
  if (kstVtxY->size() < upTo)  for (unsigned int i = kstVtxY->size(); i < upTo; i++)  kstVtxY->push_back(0);
  if (kstVtxZ->size() < upTo)  for (unsigned int i = kstVtxZ->size(); i < upTo; i++)  kstVtxZ->push_back(0);

  // ### mu+ mu- Mass ###
  if (mumuMass->size() < upTo)  for (unsigned int i = mumuMass->size(); i < upTo; i++)  mumuMass->push_back(0);
  if (mumuMassE->size() < upTo) for (unsigned int i = mumuMassE->size(); i < upTo; i++) mumuMassE->push_back(0);
  if (mumuPx->size() < upTo)    for (unsigned int i = mumuPx->size(); i < upTo; i++)    mumuPx->push_back(0);
  if (mumuPy->size() < upTo)    for (unsigned int i = mumuPy->size(); i < upTo; i++)    mumuPy->push_back(0);
  if (mumuPz->size() < upTo)    for (unsigned int i = mumuPz->size(); i < upTo; i++)    mumuPz->push_back(0);

  // ### mu+ mu- Vtx ###
  if (mumuVtxCL->size() < upTo)       for (unsigned int i = mumuVtxCL->size(); i < upTo; i++)       mumuVtxCL->push_back(0);
  if (mumuVtxX->size() < upTo)        for (unsigned int i = mumuVtxX->size(); i < upTo; i++)        mumuVtxX->push_back(0);
  if (mumuVtxY->size() < upTo)        for (unsigned int i = mumuVtxY->size(); i < upTo; i++)        mumuVtxY->push_back(0);
  if (mumuVtxZ->size() < upTo)        for (unsigned int i = mumuVtxZ->size(); i < upTo; i++)        mumuVtxZ->push_back(0);
  if (mumuCosAlphaBS->size() < upTo)  for (unsigned int i = mumuCosAlphaBS->size(); i < upTo; i++)  mumuCosAlphaBS->push_back(0);
  if (mumuCosAlphaBSE->size() < upTo) for (unsigned int i = mumuCosAlphaBSE->size(); i < upTo; i++) mumuCosAlphaBSE->push_back(0);
  if (mumuLBS->size() < upTo)         for (unsigned int i = mumuLBS->size(); i < upTo; i++)         mumuLBS->push_back(0);
  if (mumuLBSE->size() < upTo)        for (unsigned int i = mumuLBSE->size(); i < upTo; i++)        mumuLBSE->push_back(0);
  if (mumuDCA->size() < upTo)         for (unsigned int i = mumuDCA->size(); i < upTo; i++)         mumuDCA->push_back(0);

  // ### mu- ###  
  if (mumHighPurity->size() < upTo)      for (unsigned int i = mumHighPurity->size(); i < upTo; i++)      mumHighPurity->push_back(0);
  if (mumCL->size() < upTo)              for (unsigned int i = mumCL->size(); i < upTo; i++)              mumCL->push_back(0);
  if (mumNormChi2->size() < upTo)        for (unsigned int i = mumNormChi2->size(); i < upTo; i++)        mumNormChi2->push_back(0);
  if (mumPx->size() < upTo)              for (unsigned int i = mumPx->size(); i < upTo; i++)              mumPx->push_back(0);
  if (mumPy->size() < upTo)              for (unsigned int i = mumPy->size(); i < upTo; i++)              mumPy->push_back(0);
  if (mumPz->size() < upTo)              for (unsigned int i = mumPz->size(); i < upTo; i++)              mumPz->push_back(0);
  if (mumDCAVtx->size() < upTo)          for (unsigned int i = mumDCAVtx->size(); i < upTo; i++)          mumDCAVtx->push_back(0);
  if (mumDCAVtxE->size() < upTo)         for (unsigned int i = mumDCAVtxE->size(); i < upTo; i++)         mumDCAVtxE->push_back(0);
  if (mumDCABS->size() < upTo)           for (unsigned int i = mumDCABS->size(); i < upTo; i++)           mumDCABS->push_back(0);
  if (mumDCABSE->size() < upTo)          for (unsigned int i = mumDCABSE->size(); i < upTo; i++)          mumDCABSE->push_back(0);
  if (mumKinkChi2->size() < upTo)        for (unsigned int i = mumKinkChi2->size(); i < upTo; i++)        mumKinkChi2->push_back(0);
  if (mumFracHits->size() < upTo)        for (unsigned int i = mumFracHits->size(); i < upTo; i++)        mumFracHits->push_back(0);
  if (mumdxyVtx->size() < upTo)          for (unsigned int i = mumdxyVtx->size(); i < upTo; i++)          mumdxyVtx->push_back(0);
  if (mumdzVtx->size() < upTo)           for (unsigned int i = mumdzVtx->size(); i < upTo; i++)           mumdzVtx->push_back(0);
  if (mumDeltaRwithMC->size() < upTo)    for (unsigned int i = mumDeltaRwithMC->size(); i < upTo; i++)    mumDeltaRwithMC->push_back(0);
  if (mumCat->size() < upTo)             for (unsigned int i = mumCat->size(); i < upTo; i++)             mumCat->push_back("");
  if (mumNPixHits->size() < upTo)        for (unsigned int i = mumNPixHits->size(); i < upTo; i++)        mumNPixHits->push_back(0);
  if (mumNPixLayers->size() < upTo)      for (unsigned int i = mumNPixLayers->size(); i < upTo; i++)      mumNPixLayers->push_back(0);
  if (mumNTrkHits->size() < upTo)        for (unsigned int i = mumNTrkHits->size(); i < upTo; i++)        mumNTrkHits->push_back(0);
  if (mumNTrkLayers->size() < upTo)      for (unsigned int i = mumNTrkLayers->size(); i < upTo; i++)      mumNTrkLayers->push_back(0);
  if (mumNMuonHits->size() < upTo)       for (unsigned int i = mumNMuonHits->size(); i < upTo; i++)       mumNMuonHits->push_back(0);
  if (mumNMatchStation->size() < upTo)   for (unsigned int i = mumNMatchStation->size(); i < upTo; i++)   mumNMatchStation->push_back(0);
  if (mumTrig->size() < upTo)            for (unsigned int i = mumTrig->size(); i < upTo; i++)            mumTrig->push_back("");

  // ### mu+ ###  
  if (mupHighPurity->size() < upTo)      for (unsigned int i = mupHighPurity->size(); i < upTo; i++)      mupHighPurity->push_back(0);
  if (mupCL->size() < upTo)              for (unsigned int i = mupCL->size(); i < upTo; i++)              mupCL->push_back(0);
  if (mupNormChi2->size() < upTo)        for (unsigned int i = mupNormChi2->size(); i < upTo; i++)        mupNormChi2->push_back(0);
  if (mupPx->size() < upTo)              for (unsigned int i = mupPx->size(); i < upTo; i++)              mupPx->push_back(0);
  if (mupPy->size() < upTo)              for (unsigned int i = mupPy->size(); i < upTo; i++)              mupPy->push_back(0);
  if (mupPz->size() < upTo)              for (unsigned int i = mupPz->size(); i < upTo; i++)              mupPz->push_back(0);
  if (mupDCAVtx->size() < upTo)          for (unsigned int i = mupDCAVtx->size(); i < upTo; i++)          mupDCAVtx->push_back(0);
  if (mupDCAVtxE->size() < upTo)         for (unsigned int i = mupDCAVtxE->size(); i < upTo; i++)         mupDCAVtxE->push_back(0);
  if (mupDCABS->size() < upTo)           for (unsigned int i = mupDCABS->size(); i < upTo; i++)           mupDCABS->push_back(0);
  if (mupDCABSE->size() < upTo)          for (unsigned int i = mupDCABSE->size(); i < upTo; i++)          mupDCABSE->push_back(0);
  if (mupKinkChi2->size() < upTo)        for (unsigned int i = mupKinkChi2->size(); i < upTo; i++)        mupKinkChi2->push_back(0);
  if (mupFracHits->size() < upTo)        for (unsigned int i = mupFracHits->size(); i < upTo; i++)        mupFracHits->push_back(0);
  if (mupdxyVtx->size() < upTo)          for (unsigned int i = mupdxyVtx->size(); i < upTo; i++)          mupdxyVtx->push_back(0);
  if (mupdzVtx->size() < upTo)           for (unsigned int i = mupdzVtx->size(); i < upTo; i++)           mupdzVtx->push_back(0);
  if (mupDeltaRwithMC->size() < upTo)    for (unsigned int i = mupDeltaRwithMC->size(); i < upTo; i++)    mupDeltaRwithMC->push_back(0);
  if (mupCat->size() < upTo)             for (unsigned int i = mupCat->size(); i < upTo; i++)             mupCat->push_back("");
  if (mupNPixHits->size() < upTo)        for (unsigned int i = mupNPixHits->size(); i < upTo; i++)        mupNPixHits->push_back(0);
  if (mupNPixLayers->size() < upTo)      for (unsigned int i = mupNPixLayers->size(); i < upTo; i++)      mupNPixLayers->push_back(0);
  if (mupNTrkHits->size() < upTo)        for (unsigned int i = mupNTrkHits->size(); i < upTo; i++)        mupNTrkHits->push_back(0);
  if (mupNTrkLayers->size() < upTo)      for (unsigned int i = mupNTrkLayers->size(); i < upTo; i++)      mupNTrkLayers->push_back(0);
  if (mupNMuonHits->size() < upTo)       for (unsigned int i = mupNMuonHits->size(); i < upTo; i++)       mupNMuonHits->push_back(0);
  if (mupNMatchStation->size() < upTo)   for (unsigned int i = mupNMatchStation->size(); i < upTo; i++)   mupNMatchStation->push_back(0);
  if (mupTrig->size() < upTo)            for (unsigned int i = mupTrig->size(); i < upTo; i++)            mupTrig->push_back("");

  // ### K*0 track- ###
  if (kstTrkmHighPurity->size() < upTo)   for (unsigned int i = kstTrkmHighPurity->size(); i < upTo; i++)   kstTrkmHighPurity->push_back(0);
  if (kstTrkmCL->size() < upTo)           for (unsigned int i = kstTrkmCL->size(); i < upTo; i++)           kstTrkmCL->push_back(0);
  if (kstTrkmNormChi2->size() < upTo)     for (unsigned int i = kstTrkmNormChi2->size(); i < upTo; i++)     kstTrkmNormChi2->push_back(0);
  if (kstTrkmPx->size() < upTo)           for (unsigned int i = kstTrkmPx->size(); i < upTo; i++)           kstTrkmPx->push_back(0);
  if (kstTrkmPy->size() < upTo)           for (unsigned int i = kstTrkmPy->size(); i < upTo; i++)           kstTrkmPy->push_back(0);
  if (kstTrkmPz->size() < upTo)           for (unsigned int i = kstTrkmPz->size(); i < upTo; i++)           kstTrkmPz->push_back(0);
  if (kstTrkmDCAVtx->size() < upTo)       for (unsigned int i = kstTrkmDCAVtx->size(); i < upTo; i++)       kstTrkmDCAVtx->push_back(0);
  if (kstTrkmDCAVtxE->size() < upTo)      for (unsigned int i = kstTrkmDCAVtxE->size(); i < upTo; i++)      kstTrkmDCAVtxE->push_back(0);
  if (kstTrkmDCABS->size() < upTo)        for (unsigned int i = kstTrkmDCABS->size(); i < upTo; i++)        kstTrkmDCABS->push_back(0);
  if (kstTrkmDCABSE->size() < upTo)       for (unsigned int i = kstTrkmDCABSE->size(); i < upTo; i++)       kstTrkmDCABSE->push_back(0);
  if (kstTrkmFracHits->size() < upTo)     for (unsigned int i = kstTrkmFracHits->size(); i < upTo; i++)     kstTrkmFracHits->push_back(0);
  if (kstTrkmdxyVtx->size() < upTo)       for (unsigned int i = kstTrkmdxyVtx->size(); i < upTo; i++)       kstTrkmdxyVtx->push_back(0);
  if (kstTrkmdzVtx->size() < upTo)        for (unsigned int i = kstTrkmdzVtx->size(); i < upTo; i++)        kstTrkmdzVtx->push_back(0);
  if (kstTrkmDeltaRwithMC->size() < upTo) for (unsigned int i = kstTrkmDeltaRwithMC->size(); i < upTo; i++) kstTrkmDeltaRwithMC->push_back(0);
  if (kstTrkmNPixHits->size() < upTo)     for (unsigned int i = kstTrkmNPixHits->size(); i < upTo; i++)     kstTrkmNPixHits->push_back(0);
  if (kstTrkmNPixLayers->size() < upTo)   for (unsigned int i = kstTrkmNPixLayers->size(); i < upTo; i++)   kstTrkmNPixLayers->push_back(0);
  if (kstTrkmNTrkHits->size() < upTo)     for (unsigned int i = kstTrkmNTrkHits->size(); i < upTo; i++)     kstTrkmNTrkHits->push_back(0);
  if (kstTrkmNTrkLayers->size() < upTo)   for (unsigned int i = kstTrkmNTrkLayers->size(); i < upTo; i++)   kstTrkmNTrkLayers->push_back(0);
  if (kstTrkmMuMatch->size() < upTo)      for (unsigned int i = kstTrkmMuMatch->size(); i < upTo; i++)      kstTrkmMuMatch->push_back("");
 
  // ### K*0 track+ ###
  if (kstTrkpHighPurity->size() < upTo)   for (unsigned int i = kstTrkpHighPurity->size(); i < upTo; i++)   kstTrkpHighPurity->push_back(0);
  if (kstTrkpCL->size() < upTo)           for (unsigned int i = kstTrkpCL->size(); i < upTo; i++)           kstTrkpCL->push_back(0);
  if (kstTrkpNormChi2->size() < upTo)     for (unsigned int i = kstTrkpNormChi2->size(); i < upTo; i++)     kstTrkpNormChi2->push_back(0);
  if (kstTrkpPx->size() < upTo)           for (unsigned int i = kstTrkpPx->size(); i < upTo; i++)           kstTrkpPx->push_back(0);
  if (kstTrkpPy->size() < upTo)           for (unsigned int i = kstTrkpPy->size(); i < upTo; i++)           kstTrkpPy->push_back(0);
  if (kstTrkpPz->size() < upTo)           for (unsigned int i = kstTrkpPz->size(); i < upTo; i++)           kstTrkpPz->push_back(0);
  if (kstTrkpDCAVtx->size() < upTo)       for (unsigned int i = kstTrkpDCAVtx->size(); i < upTo; i++)       kstTrkpDCAVtx->push_back(0);
  if (kstTrkpDCAVtxE->size() < upTo)      for (unsigned int i = kstTrkpDCAVtxE->size(); i < upTo; i++)      kstTrkpDCAVtxE->push_back(0);
  if (kstTrkpDCABS->size() < upTo)        for (unsigned int i = kstTrkpDCABS->size(); i < upTo; i++)        kstTrkpDCABS->push_back(0);
  if (kstTrkpDCABSE->size() < upTo)       for (unsigned int i = kstTrkpDCABSE->size(); i < upTo; i++)       kstTrkpDCABSE->push_back(0);
  if (kstTrkpFracHits->size() < upTo)     for (unsigned int i = kstTrkpFracHits->size(); i < upTo; i++)     kstTrkpFracHits->push_back(0);
  if (kstTrkpdxyVtx->size() < upTo)       for (unsigned int i = kstTrkpdxyVtx->size(); i < upTo; i++)       kstTrkpdxyVtx->push_back(0);
  if (kstTrkpdzVtx->size() < upTo)        for (unsigned int i = kstTrkpdzVtx->size(); i < upTo; i++)        kstTrkpdzVtx->push_back(0);
  if (kstTrkpDeltaRwithMC->size() < upTo) for (unsigned int i = kstTrkpDeltaRwithMC->size(); i < upTo; i++) kstTrkpDeltaRwithMC->push_back(0);
  if (kstTrkpNPixHits->size() < upTo)     for (unsigned int i = kstTrkpNPixHits->size(); i < upTo; i++)     kstTrkpNPixHits->push_back(0);
  if (kstTrkpNPixLayers->size() < upTo)   for (unsigned int i = kstTrkpNPixLayers->size(); i < upTo; i++)   kstTrkpNPixLayers->push_back(0);
  if (kstTrkpNTrkHits->size() < upTo)     for (unsigned int i = kstTrkpNTrkHits->size(); i < upTo; i++)     kstTrkpNTrkHits->push_back(0);
  if (kstTrkpNTrkLayers->size() < upTo)   for (unsigned int i = kstTrkpNTrkLayers->size(); i < upTo; i++)   kstTrkpNTrkLayers->push_back(0);
  if (kstTrkpMuMatch->size() < upTo)      for (unsigned int i = kstTrkpMuMatch->size(); i < upTo; i++)      kstTrkpMuMatch->push_back("");
  
  // ### Matching Between Reconstructed and Generated ###
  if (truthMatchSignal->size() < upTo) for (unsigned int i = truthMatchSignal->size(); i < upTo; i++) truthMatchSignal->push_back(0);
  if (truthMatchMum->size() < upTo)    for (unsigned int i = truthMatchMum->size(); i < upTo; i++)    truthMatchMum->push_back(0);
  if (truthMatchMup->size() < upTo)    for (unsigned int i = truthMatchMup->size(); i < upTo; i++)    truthMatchMup->push_back(0);
  if (truthMatchTrkm->size() < upTo)   for (unsigned int i = truthMatchTrkm->size(); i < upTo; i++)   truthMatchTrkm->push_back(0);
  if (truthMatchTrkp->size() < upTo)   for (unsigned int i = truthMatchTrkp->size(); i < upTo; i++)   truthMatchTrkp->push_back(0);
}

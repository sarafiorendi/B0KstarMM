#ifndef _B0KSTMUMUTREE
#define _B0KSTMUMUTREE

#include <vector>
#include <string>
#include "TTree.h"


class B0KstMuMuTreeContent
{


 public:
  
  B0KstMuMuTreeContent () {};
  ~B0KstMuMuTreeContent () {};

  void Init ();
  void Destroy ();
  void ClearNTuple ();
  void ClearMonteCarlo ();
  void MakeTreeBranches (TTree* theTree);
  void SetBranchAddresses (TTree* theTree);
  void CopyData (B0KstMuMuTreeContent* NTupleIn, int index);
  void CopyWholeNTuple (B0KstMuMuTreeContent* NTupleIn);
  void FillWithNull (unsigned int upTo);


  // ########################################################
  // # Run Number, event number, #reco vtx and event weight #
  // ########################################################
  unsigned int              runN;
  unsigned int              eventN;
  unsigned int              recoVtxN;
  double                    evWeight;
  double                    evWeightE2;

  // ###########
  // # Trigger #
  // ###########
  std::vector<std::string>  *TrigTable;
  std::vector<int>          *TrigPrescales;
  
  // ###########################
  // # Number of B0 candidates #
  // ###########################
  unsigned int              nB;
  
  // ############################
  // # Pileup information in MC #
  // ############################
  std::vector<double>       *bunchXingMC, *numInteractionsMC, *trueNumInteractionsMC;

  // ################################
  // # Primary Vertex and Beam Spot #
  // ################################
  double                    priVtxCL, priVtxX, priVtxY, priVtxZ;
  double                    bsX, bsY;

  // ###########
  // # B0 Mass #
  // ###########
  std::vector<double>       *bMass, *bMassE, *bBarMass, *bBarMassE, *bPx, *bPy, *bPz;

  // ##########
  // # B0 Vtx #
  // ##########
  std::vector<double>       *bVtxCL, *bVtxX, *bVtxY, *bVtxZ;
  std::vector<double>       *bCosAlphaVtx, *bCosAlphaVtxE, *bCosAlphaBS, *bCosAlphaBSE;
  std::vector<double>       *bLVtx, *bLVtxE, *bLBS, *bLBSE;
  std::vector<double>       *bDCAVtx, *bDCAVtxE, *bDCABS, *bDCABSE;

  // ###########
  // # B0 ctau #
  // ###########
  std::vector<double>       *bctauPVBS, *bctauPVBSE;

  // ############
  // # K*0 Mass #
  // ############
  std::vector<double>       *kstMass, *kstMassE, *kstBarMass, *kstBarMassE, *kstPx, *kstPy, *kstPz;

  // ###########
  // # K*0 Vtx #
  // ###########
  std::vector<double>       *kstVtxCL, *kstVtxX, *kstVtxY, *kstVtxZ;

  // ################
  // # mu+ mu- Mass #
  // ################
  std::vector<double>       *mumuMass, *mumuMassE, *mumuPx, *mumuPy, *mumuPz;

  // ###############
  // # mu+ mu- Vtx #
  // ###############
  std::vector<double>       *mumuVtxCL, *mumuVtxX, *mumuVtxY, *mumuVtxZ;
  std::vector<double>       *mumuCosAlphaBS, *mumuCosAlphaBSE;
  std::vector<double>       *mumuLBS, *mumuLBSE;
  std::vector<double>       *mumuDCA;

  // #######
  // # mu- #
  // #######
  std::vector<bool>         *mumHighPurity;
  std::vector<double>       *mumCL, *mumNormChi2, *mumPx, *mumPy, *mumPz;
  std::vector<double>       *mumDCAVtx, *mumDCAVtxE, *mumDCABS, *mumDCABSE, *mumKinkChi2, *mumFracHits;
  std::vector<double>       *mumdxyVtx, *mumdzVtx;
  std::vector<std::string>  *mumCat;
  std::vector<int>          *mumNPixHits, *mumNPixLayers, *mumNTrkHits, *mumNTrkLayers, *mumNMuonHits, *mumNMatchStation;
  std::vector<std::string>  *mumTrig;

  // #######
  // # mu+ #
  // #######
  std::vector<bool>         *mupHighPurity;
  std::vector<double>       *mupCL, *mupNormChi2, *mupPx, *mupPy, *mupPz;
  std::vector<double>       *mupDCAVtx, *mupDCAVtxE, *mupDCABS, *mupDCABSE, *mupKinkChi2, *mupFracHits;
  std::vector<double>       *mupdxyVtx, *mupdzVtx;
  std::vector<std::string>  *mupCat;
  std::vector<int>          *mupNPixHits, *mupNPixLayers, *mupNTrkHits, *mupNTrkLayers, *mupNMuonHits, *mupNMatchStation;
  std::vector<std::string>  *mupTrig;

  // ##############
  // # K*0 track- #
  // ##############
  std::vector<bool>         *kstTrkmHighPurity;
  std::vector<double>       *kstTrkmCL, *kstTrkmNormChi2, *kstTrkmPx, *kstTrkmPy, *kstTrkmPz;
  std::vector<double>       *kstTrkmDCAVtx, *kstTrkmDCAVtxE, *kstTrkmDCABS, *kstTrkmDCABSE, *kstTrkmFracHits;
  std::vector<double>       *kstTrkmdxyVtx, *kstTrkmdzVtx;
  std::vector<int>          *kstTrkmNPixHits, *kstTrkmNPixLayers, *kstTrkmNTrkHits, *kstTrkmNTrkLayers;
  std::vector<std::string>  *kstTrkmMuMatch;

  // ##############
  // # K*0 track+ #
  // ##############
  std::vector<bool>         *kstTrkpHighPurity;
  std::vector<double>       *kstTrkpCL, *kstTrkpNormChi2, *kstTrkpPx, *kstTrkpPy, *kstTrkpPz;
  std::vector<double>       *kstTrkpDCAVtx, *kstTrkpDCAVtxE, *kstTrkpDCABS, *kstTrkpDCABSE, *kstTrkpFracHits;
  std::vector<double>       *kstTrkpdxyVtx, *kstTrkpdzVtx;
  std::vector<int>          *kstTrkpNPixHits, *kstTrkpNPixLayers, *kstTrkpNTrkHits, *kstTrkpNTrkLayers;
  std::vector<std::string>  *kstTrkpMuMatch;

  // #########################
  // # Generated Observables #
  // #########################
  int                       genSignal; // 1 = B0 --> K*0(K+pi-) mu+mu-
                                       // 2 = B0bar --> K*0bar(K-pi+) mu+mu-
                                       // 3 = B0 --> K*0(K+pi-) J/psi(mu+mu-)
                                       // 4 = B0bar --> K*0bar(K-pi+) J/psi(mu+mu-)
                                       // 5 = B0 --> K*0(K-pi+) psi(2S)(mu+mu-)
                                       // 6 = B0bar --> K*0bar(K-pi+) psi(2S)(mu+mu-)
  int                       genMuMuBG, genMuMuBGnTrksm, genMuMuBGnTrksp;
  bool                      genPsiPrompt;
  bool                      genSignHasFSR, genSignKstHasFSR, genSignPsiHasFSR;
  
  // ############################
  // # Generated Primary Vertex #
  // ############################
  double                    genPriVtxX, genPriVtxY, genPriVtxZ;

  // #####################
  // # Generated B0 Mass #
  // #####################
  double                    genB0Mass, genB0Px, genB0Py, genB0Pz;

  // ####################
  // # Generated B0 Vtx #
  // ####################
  double                    genB0VtxX, genB0VtxY, genB0VtxZ;

  // ######################
  // # Generated K*0 Mass #
  // ######################
  double                    genKstMass, genKstPx, genKstPy, genKstPz;

  // #####################
  // # Generated K*0 Vtx #
  // #####################
  double                    genKstVtxX, genKstVtxY, genKstVtxZ;

  // ###########################################
  // # Generated J/psi or psi(2S) Mass and Vtx #
  // ###########################################
  double                    genPsiMass, genPsiVtxX, genPsiVtxY, genPsiVtxZ;

  // #################
  // # Generated mu- #
  // #################
  int                       genMumMother;
  double                    genMumPx, genMumPy, genMumPz;
  bool                      trueMumTriggered, trueMumInAcceptance;

  // #################
  // # Generated mu+ #
  // #################
  int                       genMupMother;
  double                    genMupPx, genMupPy, genMupPz;
  bool                      trueMupTriggered, trueMupInAcceptance;

  // ########################
  // # Generated K*0 track- #
  // ########################
  int                       genKstTrkmMother;
  double                    genKstTrkmPx, genKstTrkmPy, genKstTrkmPz;

  // ########################
  // # Generated K*0 track+ #
  // ########################
  int                       genKstTrkpMother;
  double                    genKstTrkpPx, genKstTrkpPy, genKstTrkpPz;

  // ################################################
  // # Matching Between Reconstructed and Generated #
  // ################################################
  std::vector<bool>         *truthMatchSignal, *truthMatchMum, *truthMatchMup, *truthMatchTrkm, *truthMatchTrkp;

  
 private:

  void CopyScalars (B0KstMuMuTreeContent* NTupleIn);
  void CopyVectors (B0KstMuMuTreeContent* NTupleIn, int index);

};

#endif

#ifndef _B0KSTMUMUSINGLECANDTREE
#define _B0KSTMUMUSINGLECANDTREE

#include "TTree.h"
#include "../interface/B0KstMuMuTreeContent.h"


class B0KstMuMuSingleCandTreeContent : public B0KstMuMuTreeContent
{


 public:
  
  B0KstMuMuSingleCandTreeContent ();
  ~B0KstMuMuSingleCandTreeContent ();

  void Init ();
  void ClearNTuple ();
  void MakeTreeBranches (TTree* theTree);
  void SetBranchAddresses (TTree* theTree);
  void CopyWholeNTuple (B0KstMuMuSingleCandTreeContent* NTupleIn);


  // ################
  // # B0 variables #
  // ################
  double B0MassArb;
  double B0pT, B0Eta, B0Phi;
  bool   B0notB0bar;

  // ####################
  // # Trigger category #
  // ####################
  int TrigCat;

  // ####################
  // # Angular anaylsis #
  // ####################
  double CosThetaKArb, CosThetaMuArb, PhiKstMuMuPlaneArb;

 private:
  
};

#endif

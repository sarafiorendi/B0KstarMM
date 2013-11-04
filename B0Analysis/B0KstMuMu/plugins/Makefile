CC            = g++
CCFLAGS       = -g -Wall -m64 -O3

MYSOURCEDIR   = ../src
MYSOURCES     = $(MYSOURCEDIR)/B0KstMuMuTreeContent.cc $(MYSOURCEDIR)/B0KstMuMuSingleCandTreeContent.cc $(MYSOURCEDIR)/Utils.cc $(MYSOURCEDIR)/ReadParameters.cc

MYINCDIR      = ../interface
MYINCFLAGS    = -I $(MYINCDIR)

XERCESINFLAGS = -I $(HOME)/xerces/include
BOOSTINFLAGS  = -I $(HOME)/boost
GSLINCFLAGS   = -I $(HOME)/gsl/include

XERCESFLAGS   = -L $(HOME)/xerces/lib -lxerces-c
BOOSTFLFLAGS  = -L $(HOME)/boost/stage/lib
GSLFLAGS      = -L $(HOME)/gsl/lib -lgsl
ROOTFLAGS     = `root-config --cflags --glibs` -lMinuit
ROOFITFLAGS   = `root-config --cflags --glibs` -lRooFit -lRooFitCore -lFoam -lMinuit


SingleCandB0KstMuMu: $(MYSOURCEDIR) $(MYINCDIR) SingleCandB0KstMuMu.cc
	$(CC) $(CCFLAGS) -o SingleCandB0KstMuMu $(ROOTFLAGS) $(MYINCFLAGS) $(MYSOURCES) SingleCandB0KstMuMu.cc

AddVars2Candidates: $(MYSOURCEDIR) $(MYINCDIR) AddVars2Candidates.cc
	$(CC) $(CCFLAGS) -o AddVars2Candidates $(ROOTFLAGS) $(MYINCFLAGS) $(MYSOURCES) AddVars2Candidates.cc

ExtractYield: $(MYSOURCEDIR) $(MYINCDIR) ExtractYield.cc
	$(CC) $(CCFLAGS) -o ExtractYield $(ROOFITFLAGS) $(MYINCFLAGS) $(MYSOURCES) ExtractYield.cc

ComputeEfficiency: $(MYSOURCEDIR) $(MYINCDIR) ComputeEfficiency.cc
	$(CC) $(CCFLAGS) -o ComputeEfficiency $(ROOTFLAGS) $(MYINCFLAGS) $(MYSOURCES) ComputeEfficiency.cc

MakePlots: $(MYSOURCEDIR) $(MYINCDIR) MakePlots.cc
	$(CC) $(CCFLAGS) -o MakePlots $(ROOTFLAGS) $(MYINCFLAGS) $(MYSOURCES) MakePlots.cc

B0KstMuMuScanCuts: $(MYSOURCEDIR) $(MYINCDIR) B0KstMuMuScanCuts.cc
	$(CC) $(CCFLAGS) -o B0KstMuMuScanCuts $(ROOTFLAGS) $(MYINCFLAGS) $(MYSOURCES) B0KstMuMuScanCuts.cc

clean:
	rm -rf SingleCandB0KstMuMu AddVars2Candidates ExtractYield ComputeEfficiency MakePlots B0KstMuMuScanCuts SingleCandB0KstMuMu.dSYM AddVars2Candidates.dSYM ExtractYield.dSYM ComputeEfficiency.dSYM MakePlots.dSYM B0KstMuMuScanCuts.dSYM Macros_C.*

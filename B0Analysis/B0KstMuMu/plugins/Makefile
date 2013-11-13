CXX           := g++
CXXFLAGS      := -g -ansi -Wall -m64 -O3

SOURCEDIR     := ../src
TEST          := $(wildcard $(SOURCEDIR)/*.cc)
SOURCES       := $(filter-out $(SOURCEDIR)/B0KstMuMu.cc,$(TEST))

INCLUDEDIR    := ../interface
INCLUDES      := -I $(INCLUDEDIR)

XERCESINFLAGS := -I $(HOME)/xerces/include
BOOSTINFLAGS  := -I $(HOME)/boost
GSLINCFLAGS   := -I $(HOME)/gsl/include

XERCESFLAGS   := -L $(HOME)/xerces/lib -lxerces-c
BOOSTFLFLAGS  := -L $(HOME)/boost/stage/lib
GSLFLAGS      := -L $(HOME)/gsl/lib -lgsl

ROOTFLAGS     := `root-config --cflags --glibs` -lMinuit
ROOFITFLAGS   := $(ROOTFLAGS) -lRooFit -lRooFitCore -lFoam

LOCALCC       := $(wildcard *.cc)


SingleCandB0KstMuMu: $(SOURCEDIR) $(INCLUDEDIR) SingleCandB0KstMuMu.cc
	$(CXX) $(CXXFLAGS) -o $@ $(ROOTFLAGS) $(INCLUDES) $(SOURCES) $@.cc

AddVars2Candidates: $(SOURCEDIR) $(INCLUDEDIR) AddVars2Candidates.cc
	$(CXX) $(CXXFLAGS) -o $@ $(ROOTFLAGS) $(INCLUDES) $(SOURCES) $@.cc

B0KstMuMuScanCuts: $(SOURCEDIR) $(INCLUDEDIR) B0KstMuMuScanCuts.cc
	$(CXX) $(CXXFLAGS) -o $@ $(ROOTFLAGS) $(INCLUDES) $(SOURCES) $@.cc

ComputeEfficiency: $(SOURCEDIR) $(INCLUDEDIR) ComputeEfficiency.cc
	$(CXX) $(CXXFLAGS) -o $@ $(ROOTFLAGS) $(INCLUDES) $(SOURCES) $@.cc

MakePlots: $(SOURCEDIR) $(INCLUDEDIR) MakePlots.cc
	$(CXX) $(CXXFLAGS) -o $@ $(ROOTFLAGS) $(INCLUDES) $(SOURCES) $@.cc

ExtractYield: $(SOURCEDIR) $(INCLUDEDIR) ExtractYield.cc
	$(CXX) $(CXXFLAGS) -o $@ $(ROOFITFLAGS) $(INCLUDES) $(SOURCES) $@.cc

clean:
	$(RM)r $(LOCALCC:%.cc=%) $(LOCALCC:%.cc=%.dSYM) Macros_C.*

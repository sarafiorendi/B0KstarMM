OS            := $(shell uname)

ifeq ($(OS), Darwin)
CXX           := clang++
DEBUG         :=
endif

ifeq ($(OS), Linux)
CXX           := g++
DEBUG         := -g
endif

CXXFLAGS      := $(DEBUG) -ansi -Wall -Wextra -m64 -O3

SOURCEDIR     := ../src
SOURCES       := $(filter-out $(SOURCEDIR)/B0KstMuMu.cc,$(wildcard $(SOURCEDIR)/*.cc))

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


SingleCandB0KstMuMu: $(SOURCEDIR) $(INCLUDEDIR) SingleCandB0KstMuMu.cc
	$(CXX) $(CXXFLAGS) -o $@ $(ROOTFLAGS) $(INCLUDES) $(SOURCES) $@.cc

AddVarsToCandidates: $(SOURCEDIR) $(INCLUDEDIR) AddVarsToCandidates.cc
	$(CXX) $(CXXFLAGS) -o $@ $(ROOTFLAGS) $(INCLUDES) $(SOURCES) $@.cc

B0KstMuMuScanCuts: $(SOURCEDIR) $(INCLUDEDIR) B0KstMuMuScanCuts.cc
	$(CXX) $(CXXFLAGS) -o $@ $(ROOTFLAGS) $(INCLUDES) $(SOURCES) $@.cc

ComputeEfficiency: $(SOURCEDIR) $(INCLUDEDIR) ComputeEfficiency.cc
	$(CXX) $(CXXFLAGS) -o $@ $(ROOTFLAGS) $(INCLUDES) $(SOURCES) $@.cc

MakePlots: $(SOURCEDIR) $(INCLUDEDIR) MakePlots.cc
	$(CXX) $(CXXFLAGS) -o $@ $(ROOFITFLAGS) $(INCLUDES) $(SOURCES) -D"ROOFIT" $@.cc

ExtractYield: $(SOURCEDIR) $(INCLUDEDIR) ExtractYield.cc
	$(CXX) $(CXXFLAGS) -o $@ $(ROOFITFLAGS) $(INCLUDES) $(SOURCES) -D"ROOFIT" $@.cc

LOCALCC := $(wildcard *.cc)
clean:
	$(RM) $(LOCALCC:%.cc=%) Macros_C.*

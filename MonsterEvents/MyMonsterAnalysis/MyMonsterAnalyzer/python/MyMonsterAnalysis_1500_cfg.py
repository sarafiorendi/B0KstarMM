################################
######### CMSSW_3_3_5 ##########
### MyMonsterAnalysis_cfg.py ###
################################

import FWCore.ParameterSet.Config as cms
process = cms.Process("MonsterAnalysis")

#######################
### Configuring CMS ###
#######################

###### Geometry ######
process.load("Configuration.StandardSequences.Geometry_cff")


###### Magnetic Field ######
process.load("Configuration.StandardSequences.MagneticField_cff")


###### Services ######
process.load("Configuration.StandardSequences.Services_cff")


###### Calibration & Alignment ######
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# RERECO data taking january 29th 2010
process.GlobalTag.globaltag = "GR09_R_V6A::All"


###### Message Logger ######
process.load("FWCore.MessageService.MessageLogger_cfi")


###### Configuring Tracking ######
process.load("RecoTracker.Configuration.RecoTracker_cff")
process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")


###### My Analyzer ######
process.MyProcess = cms.EDAnalyzer("MyMonsterAnalyzer",
                                   inputTag = cms.InputTag("gtDigis"),
                                   PrintMsg = cms.bool(False),
                                   MinHits  = cms.uint32(1500))


###### Output of MyProgram ######
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("/afs/cern.ch/user/d/dinardo/scratch0/MyMonsterAnalysisOutput/MyMonsterAnalysis_1500dg_AllEvents_Jan29th.root"),
#                                   fileName = cms.string("MyMonsterAnalysis.root"),
                                   closeFileFast = cms.untracked.bool(True))


###### Which data ######
process.load("DataDec09_ReRecoMinBias_Jan29th_Skim_GoodRunsM_cff")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


###### RUN ######
process.p = cms.Path(process.siPixelRecHits*process.MyProcess)

process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.cerr.threshold = "INFO"

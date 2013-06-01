####################################
########### CMSSW_3_3_5 ############
### MyMonsterAnalysisSkim_cfg.py ###
####################################

import FWCore.ParameterSet.Config as cms
process = cms.Process("MonsterAnalysisSkim")

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

# RERECO data taking december 19th 2009
#process.GlobalTag.globaltag = "GR09_R_V5::All"

# RERECO data taking january 29th 2010
process.GlobalTag.globaltag = "GR09_R_V6A::All"


###### Message Logger ######
process.load("FWCore.MessageService.MessageLogger_cfi")


###### Output of this python ######
#MyEventContent = cms.PSet(outputCommands = cms.untracked.vstring('drop *', 'keep  *_*_*_PixAnalysis'))
process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string("MyMonsterAnalysisSkim.root"))


###### Which data ######
#process.load("DataDec09_ReRecoMinBias_Dec19th_cff")
process.load("DataDec09_ReRecoMinBias_Jan29th_cff")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


###### RUN ######
process.p = cms.Path(process.out)

process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.cerr.threshold = "INFO"

###### CMSSW_3_3_5 ######
###### MyRecoTracksCosmics_cfg.py ######

import FWCore.ParameterSet.Config as cms

process = cms.Process("CosmicAnalysis")

###### Configuring CMS ######
###### Geometry ######
process.load("Configuration.StandardSequences.Geometry_cff")

###### Magnetic Field ######
process.load("Configuration.StandardSequences.MagneticField_0T_cff")
# Magnetic field at 0.0T
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
# Magnetic field at 3.8T

###### Calibration & Alignment ######
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.connect = "frontier://FrontierProd/CMS_COND_21X_GLOBALTAG"
process.es_prefer_GlobalTag = cms.ESPrefer("PoolDBESSource","GlobalTag")
# CRAFT 0T data set
#process.GlobalTag.globaltag = "CRAFT_0T::All"
# CRAFT 3.8T data set
process.GlobalTag.globaltag = "CRAFT_ALL_V12::All"
# Il global tag serve a cercare nel data base tutti i parametri di
# calibrazione, allineamento, ecc. corrispondenti ad un certo run.
# E' quello che dice al programma di ricostruzione le condizioni del detector.

###### Message Logger ######
process.load("FWCore.MessageService.MessageLogger_cfi")
# Impostazioni di stampa dei messaggi su schermo: default

###### Configuring Tracking ######
###### Configuring Reconstruction Cosmics ######
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.ctfRefitter = process.TrackRefitterP5.clone()
process.ctfRefitter.src = "ctfWithMaterialTracksP5"
process.ctfRefitter.TrajectoryInEvent = True

###### Load Pixel Event Filters: Number of Muons and Timing ######
process.load("DPGAnalysis.SiPixelTools.FEDInRunFilter_cfi")
#process.load("DPGAnalysis.SiPixelTools.muonTOF_cfi")
#process.MuonTOFFilter_trackQuality.max_goodmuons = 2
#process.MuonTOFFilter_trackQuality.max_timeError = 15
#process.MuonTOFFilter_trackQuality.max_chi2_ndof = 15

###### My Analyzer ######
process.MyProcess = cms.EDAnalyzer("MyCosmAnalyzer",
                                   inputTag = cms.InputTag("gtDigis"),
                                   trajectoryInput = cms.string("ctfRefitter"),
                                   tracks = cms.untracked.InputTag("ctfWithMaterialTracksP5"))


process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("MyRecoTracksCosmics.root"),
                                   closeFileFast = cms.untracked.bool(True))

				   
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
# CRAFT 0T data set
#    "/store/data/Commissioning08/Cosmics/RECO/CRAFT_0T_229-v1/0000/38F99BDA-8C36-DE11-BBFC-001731AF698D.root"))
# CRAFT 3.8T data set
    "/store/data/Commissioning08/Cosmics/RAW-RECO/CRAFT_ALL_V12_229_Tosca090322_ReReco_FromSuperPointing_v1/0005/EE37F772-2137-DE11-AB10-00304867BFAA.root"))


process.maxEvents = cms.untracked.PSet(
	    input = cms.untracked.int32(-1))


###### RUN ######
process.p = cms.Path(process.fedInRunFilter*
#		process.MuonTOFFilter_trackQuality*
                process.ctfRefitter*
                process.MyProcess)


process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.cerr.threshold = "Info"

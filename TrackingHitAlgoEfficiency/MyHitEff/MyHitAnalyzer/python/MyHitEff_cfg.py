#######################
##### CMSSW_3_3_5 #####
### MyHitEff_cfg.py ###
#######################

import FWCore.ParameterSet.Config as cms
process = cms.Process("HitEffAnalysis")

#######################
### Configuring CMS ###
#######################

###### Geometry ######
process.load("Configuration.StandardSequences.Geometry_cff")


###### Magnetic Field ######
process.load("Configuration.StandardSequences.MagneticField_cff")


###### Calibration & Alignment ######
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# Montecarlo 900GeV Peak Mode (slow shaping, 1 sample)
process.GlobalTag.globaltag = "START3X_V16D::All"

# Montecarlo 10GeV Mu
#process.GlobalTag.globaltag = "MC_31X_V9::All"

# RECO data taking december 2009
#process.GlobalTag.globaltag = "GR09_P_V7::All"

# RERECO data taking december 2009
#process.GlobalTag.globaltag = "GR09_R_V4::All"


###### Message Logger ######
process.load("FWCore.MessageService.MessageLogger_cfi")
# Impostazioni base per stampa dei messaggi su schermo
#process.MessageLogger.categories = ["CkfPattern","TrackProducer"]
#process.MessageLogger.debugModules = ['*']
#process.MessageLogger.cout = cms.untracked.PSet(
#    threshold = cms.untracked.string('DEBUG'),
#    default = cms.untracked.PSet(limit = cms.untracked.int32(0)),
#    CkfPattern = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
#    TrackProducer = cms.untracked.PSet(limit = cms.untracked.int32(-1)))
# Used to enable the print out on screen of the messages tagged as "CkfPattern"


###### Configuring Tracking ######
process.load("RecoTracker.Configuration.RecoTracker_cff")
process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")
# To do only the re-fit and not the whole tracking
#process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
#process.ctfRefitter = process.TrackRefitter.clone()
#process.ctfRefitter.src = "generalTracks"
#process.ctfRefitter.TrajectoryInEvent = True

# Calculate the local position of the cluster even for those who are not in the track
process.ttrhbwr.ComputeCoarseLocalPositionFromDisk = True

# Size of the hit search window (default = 3)
process.Chi2MeasurementEstimator.nSigma = 3

# In order to generate the "CrossingFrame" we need to run the "MixingModule"
# This will allow us to have the PSimHit together with the possibility to
# use the methods "particleType" and "processType"
process.load("SimGeneral.MixingModule.mixNoPU_cfi")


###### My Analyzer ######
process.MyProcess = cms.EDAnalyzer("MyHitAnalyzer",
                                   # Use "generalTracks" and NOT "ctfRefitter" when the whole tracking is performed
                                   # Use "ctfRefitter" and NOT "generalTracks" when only the re-fit is performed
                                   trajectoryInput             = cms.string("generalTracks"),
                                   PrintMsg                    = cms.bool(False),
                                   SimData                     = cms.bool(False),
                                   RealData                    = cms.bool(False),
                                   ClusterPlots                = cms.bool(False),
                                   PixBackSearch               = cms.bool(True),
                                   MaxDistance                 = cms.double(5.00),
                                   MaxChi2                     = cms.double(30.0),
                                   MaxDPhi                     = cms.double(1.6),
                                   MinNValidHits               = cms.int32(8),
                                   EdgeTollerance              = cms.int32(1),
                                   LongitudinalDistanceOverlap = cms.double(0.05),
                                   RadialDistanceOverlap       = cms.double(0.05),
                                   MinPt                       = cms.double(1.0),
                                   MaxEta                      = cms.double(2.5),
                                   MaxTracks                   = cms.uint32(100),
                                   TrackAssociatorByHitsPSet = cms.PSet(
        Quality_SimToReco = cms.double(0.5),
        # associateRecoTracks is False if we use the "MixingModule"
        associateRecoTracks = cms.bool(True),
        UseGrouped = cms.bool(True),
        associatePixel = cms.bool(True),
        ROUList = cms.vstring("g4SimHitsTrackerHitsTIBLowTof","g4SimHitsTrackerHitsTIBHighTof",
                              "g4SimHitsTrackerHitsTIDLowTof","g4SimHitsTrackerHitsTIDHighTof",
                              "g4SimHitsTrackerHitsTOBLowTof","g4SimHitsTrackerHitsTOBHighTof",
                              "g4SimHitsTrackerHitsTECLowTof","g4SimHitsTrackerHitsTECHighTof",
                              "g4SimHitsTrackerHitsPixelBarrelLowTof","g4SimHitsTrackerHitsPixelBarrelHighTof",
                              "g4SimHitsTrackerHitsPixelEndcapLowTof","g4SimHitsTrackerHitsPixelEndcapHighTof"),
        UseSplitting = cms.bool(True),
        UsePixels = cms.bool(True),
        ThreeHitTracksAreSpecial = cms.bool(True),
        AbsoluteNumberOfHits = cms.bool(False),
        associateStrip = cms.bool(True),
        Purity_SimToReco = cms.double(0.75),
        Cut_RecoToSim = cms.double(0.75),
        SimToRecoDenominator = cms.string("sim"))) # OR "reco"


# Output of MyHitAnalyzer
process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string("/afs/cern.ch/user/d/dinardo/scratch0/MyHitEffOutput/MyHitEff.root"),
                                   fileName = cms.string("MyHitEff.root"),
                                   closeFileFast = cms.untracked.bool(True))


# Output of this python
MyEventContent = cms.PSet(outputCommands = cms.untracked.vstring('drop *', 'keep  *_*_*_HitEffAnalysis'))
process.out = cms.OutputModule("PoolOutputModule", outputCommands = MyEventContent.outputCommands,              
                               fileName = cms.untracked.string("/afs/cern.ch/user/d/dinardo/scratch0/MyHitEffOutput/MyHitEffData.root"))


###### Which data ######
#process.load("MC_3_3_5_10GeVMu_cff")
process.load("MC_MinBiasStartUp900GeV_Peak_cff")
#process.load("DataDec09_RecoMinBias_Dec14th_cff")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
#process.source.skipEvents = cms.untracked.uint32(447)

# Uncomment in order to remove outlier rejection
process.KFFittingSmootherWithOutliersRejectionAndRK.EstimateCut = -1
process.fifthFittingSmootherWithOutlierRejection.EstimateCut = -1

# Uncomment in order to remove the requirement of having the same seed hits in the rebuild
#process.GroupedCkfTrajectoryBuilderP5.requireSeedHitsInRebuild = False
#process.fifthCkfTrajectoryBuilder.requireSeedHitsInRebuild = False
#process.newTrajectoryBuilder.requireSeedHitsInRebuild = False
#process.fourthCkfTrajectoryBuilder.requireSeedHitsInRebuild = False
#process.GroupedCkfTrajectoryBuilder.requireSeedHitsInRebuild = False
#process.thCkfTrajectoryBuilder.requireSeedHitsInRebuild = False
#process.stepOneCkfTrajectoryBuilder.requireSeedHitsInRebuild = False
#process.secCkfTrajectoryBuilder.requireSeedHitsInRebuild = False

#process.GroupedCkfTrajectoryBuilderP5.bestHitOnly = False
#process.fifthCkfTrajectoryBuilder.bestHitOnly = False
#process.newTrajectoryBuilder.bestHitOnly = False
#process.fourthCkfTrajectoryBuilder.bestHitOnly = False
#process.GroupedCkfTrajectoryBuilder.bestHitOnly = False
#process.thCkfTrajectoryBuilder.bestHitOnly = False
#process.stepOneCkfTrajectoryBuilder.bestHitOnly = False
#process.secCkfTrajectoryBuilder.bestHitOnly = False

# Print the summary of the CMSSW running
#process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True)) 


###### L1 TRIGGER SELECTION ######
process.load("L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff")
#process.load("L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff")
process.load("HLTrigger.HLTfilters.hltLevel1GTSeed_cfi")
process.L1T1 = process.hltLevel1GTSeed.clone()
# Select events on Tech. trigger bits
process.L1T1.L1TechTriggerSeeding = cms.bool(True)
process.L1T1.L1SeedsLogicalExpression = cms.string("(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)")
# Select events on Algo. trigger bits
#process.L1T1.L1TechTriggerSeeding = cms.bool(False)
#process.L1T1.L1SeedsLogicalExpression = cms.string("(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)")


###### RUN ######
process.p = cms.Path(
                     ### Only with DATA ###
                     process.L1T1*

                     ### Only with MC Rel-Val ###
                     # Run "mix" in order to generate the "CrossingFrame"
#                     process.mix*

                     ### Allways ###
                     process.siPixelRecHits*
                     process.siStripMatchedRecHits*
                     # Use "ckftracks" and NOT "ctfRefitter" to do the whole tracking
                     # Use "ctfRefitter" and NOT "ckftracks" to do only the re-fit
                     process.ckftracks*
#                     process.ctfRefitter*
                     process.MyProcess)
#                     process.out)

process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.cerr.threshold = "INFO"

# Uncomment in order to print out all the python config files called
#print process.dumpPython()


### Controllare di aver impostato ###
# 1. Global Tag
# 2. TFileService: batch-run --> whole path; crab-run --> file name only
# 3. Parametri MyHitAnalyzer
# 4. True/False per associateRecoTracks
# 5. File con i dati
# 6. Numero di eventi
# 7. Fase di Run: con/senza L1T - con/senza mix

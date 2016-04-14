############################
### MyPixAnalysis_cfg.py ###
############################

import FWCore.ParameterSet.Config as cms
process = cms.Process("PixAnalysis")


#######################
### Configuring CMS ###
#######################

###### Geometry ######
process.load("Configuration.StandardSequences.Geometry_cff")


###### Magnetic Field ######
process.load("Configuration.StandardSequences.MagneticField_cff")


###### Services ######
process.load("Configuration.StandardSequences.Services_cff")


###### Message Logger ######
process.load("FWCore.MessageService.MessageLogger_cfi")


###### Configuring TrackingParticles ######
process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")
process.load("SimTracker.TrackAssociation.LhcParametersDefinerForTP_cfi")


###### Configuring LocalReco ######
process.load("RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff")
process.load("RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi")


###### Configuring Pixel-Only Tracking ######
process.load("RecoPixelVertexing.PixelTrackFitting.PixelTracks_cff")
process.load("RecoPixelVertexing.PixelVertexFinding.PixelVertexes_cfi")


###### Configuring RecoTracker TkTrackingRegions ######
process.load("RecoTracker.TkTrackingRegions.GlobalTrackingRegion_cfi")
process.load("RecoTracker.TkTrackingRegions.GlobalTrackingRegionFromBeamSpot_cfi")


###### Configuring Triplets ######
process.load("RecoPixelVertexing.PixelTriplets.PixelTripletHLTGenerator_cfi")


###### Configuring Vertexing ######
process.load("RecoVertex.Configuration.RecoVertex_cff")
process.load("RecoVertex.PrimaryVertexProducer.OfflinePixel3DPrimaryVertices_cfi")


###### Calibration & Alignment ######
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR_R_36X_V12A::All"


###################
### My Analyzer ###
###################

process.MyProcess = cms.EDAnalyzer("MyPixAnalyzer",

                                   #####################
                                   # Global parameters #
                                   #####################
                                   PrintMsg         = cms.bool(False),
                                   IsMC             = cms.bool(False),
                                   IsTrkPart        = cms.bool(False),
                                   AnalysisType     = cms.string("Track"), # "Track" OR "Vertex" OR "TrackANDVertex"
                                   MinTkTracks      = cms.uint32(10),
                                   MaxTkTracks      = cms.uint32(100),
                                   MinHitsMatch     = cms.uint32(3),
                                   TkXYVxErrCorr    = cms.double(0.89),
                                   TkZVxErrCorr     = cms.double(0.87),

                                   ####################
                                   # Track parameters #
                                   ####################
                                   MaxEtaTkTrk      = cms.double(2.3),
                                   MaxChi2PxTrk     = cms.double(15.0),
                                   MaxChi2TkTrk     = cms.double(5.0),
                                   RangePt          = cms.double(8.0), # Must be an integer number of PtStep
                                   PtStep           = cms.double(0.25),
                                   RangeEta         = cms.double(2.3), # Must be an integer number of EtaStep
                                   EtaStep          = cms.double(0.1),
                                   RangePhi         = cms.double(180.0),
                                   PhiStep          = cms.double(6.0),
                                   RangeChi2        = cms.double(30.0),
                                   Chi2Step         = cms.double(3.0),
                                   TrBound          = cms.double(0.1), # With respec to to Beam Spot or parent vertex
                                   TzBound          = cms.double(1.0), # With respec to to Beam Spot or parent vertex
                                   MinValidHitsPx   = cms.uint32(3),
                                   MinValidHitsTk   = cms.uint32(8),
                                   MinTrkVxDoF      = cms.double(40.0),
                                   MinTrkVxWgt      = cms.double(0.85),
                                   MinPtTk          = cms.double(1.0),

                                   #####################
                                   # Vertex parameters #
                                   #####################
                                   VxInputTag       = cms.InputTag("pixelVertices"),
                                   MaxEtaVxTrk      = cms.double(2.3),
                                   MaxChi2VxTrk     = cms.double(30.0),
                                   VrBound          = cms.double(0.1),  # With respecto to Beam Spot
                                   VzBound          = cms.double(10.0), # With respecto to Beam Spot
                                   MinVxDoF         = cms.double(2.0),
                                   MinVxTrkMatch    = cms.uint32(2),
                                   PxVxErrCorr      = cms.double(1.5),
                                   MinPtVx          = cms.double(0.9),

                                   TrackAssociatorByHitsPSet = cms.PSet(
        ROUList                  = cms.vstring("g4SimHitsTrackerHitsTIBLowTof","g4SimHitsTrackerHitsTIBHighTof",
                                               "g4SimHitsTrackerHitsTIDLowTof","g4SimHitsTrackerHitsTIDHighTof",
                                               "g4SimHitsTrackerHitsTOBLowTof","g4SimHitsTrackerHitsTOBHighTof",
                                               "g4SimHitsTrackerHitsTECLowTof","g4SimHitsTrackerHitsTECHighTof",
                                               "g4SimHitsTrackerHitsPixelBarrelLowTof","g4SimHitsTrackerHitsPixelBarrelHighTof",
                                               "g4SimHitsTrackerHitsPixelEndcapLowTof","g4SimHitsTrackerHitsPixelEndcapHighTof"),
        UseGrouped               = cms.bool(True),
        UseSplitting             = cms.bool(True),
        UsePixels                = cms.bool(True), # It will consider also pixel hits
        ThreeHitTracksAreSpecial = cms.bool(True), # If the track has only three hits, then all of them must be matched
        AbsoluteNumberOfHits     = cms.bool(False),
        associatePixel           = cms.bool(True),
        associateStrip           = cms.bool(True),
        associateRecoTracks      = cms.bool(True),
        Cut_RecoToSim            = cms.double(0.75),
        Quality_SimToReco        = cms.double(0.5), # Used only if SimToRecoDenominator = "sim"
        Purity_SimToReco         = cms.double(0.75),
        SimToRecoDenominator     = cms.string("reco")))


###############
### Filters ###
###############

###### Filter based on good primary vertices ######
process.primaryVertexFilter = cms.EDFilter("VertexSelector",
                                           src = cms.InputTag("offlinePrimaryVertices"),
                                           cut = cms.string("!isFake && ndof > 4 && abs(z) <= 15 && position.Rho <= 2"),
                                           filter = cms.bool(True)) # Otherwise it won't filter the events, just produce an empty vertex collection


###### Filter scraping events ######
# It cuts on the fraction of high purity tracks
# The current setup only cuts the event in which there are more than 10 high purity tracks
process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25))


###### L1 ######
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import hltLevel1GTSeed
process.L1MinBias = hltLevel1GTSeed.clone(
    L1TechTriggerSeeding = cms.bool(True),
    L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39) AND NOT ((42 AND NOT 43) OR (43 AND NOT 42))'))


###### HLT ######
process.HLTMinBias = cms.EDFilter("HLTHighLevel",
                                  TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                  HLTPaths = cms.vstring('HLT_L1_BscMinBiasOR_BptxPlusORMinus'), # Provide list of HLT paths (or patterns) you want
                                  eventSetupPathsKey = cms.string(''), # Not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
                                  andOr = cms.bool(True),              # How to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
                                  throw = cms.bool(False))             # Throw exception on unknown path names


######################
### Input / Output ###
######################

###### Which data ######
process.load("DataMay06_RecoMinBias_Run132605_cff")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


###### Output of MyProgram ######
process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string("/afs/cern.ch/user/d/dinardo/scratch0/MyPixAnalysisOutput/MyPixAnalysisData_Tk.root"),
                                   fileName = cms.string("MyPixAnalysisData_Tk.root"),
                                   closeFileFast = cms.untracked.bool(True))


###########
### RUN ###
###########

process.p = cms.Path(process.L1MinBias * process.HLTMinBias * process.noscraping * process.primaryVertexFilter * process.MyProcess)

process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.cerr.threshold = "INFO"

# Uncomment in order to print out all the python config files called
#print process.dumpPython()

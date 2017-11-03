#################
### Variables ###
#################
runDataMC          = 1 # 1 = Data; 2 = MC (Reco + Gen); 3 = MC (Gen)
useJSON            = False
printMsg           = False
run2016            = True
HLT_LMNR           = True
HLT_Jpsi           = False
HLT_PsiPrime       = False
triggerProcessName = 'HLT' # 'GEN' or 'HLT' or 'RECO' or 'TEST' or ...

print "\n@@@ CMSSW run configuration flags @@@"
print "runDataMC          : ", runDataMC
print "useJSON            : ", useJSON
print "printMsg           : ", printMsg
print "run2016            : ", run2016
print "dataset LMNR       : ", HLT_LMNR
print "dataset Jpsi       : ", HLT_Jpsi
print "dataset PsiPrime   : ", HLT_PsiPrime
print "triggerProcessName : ", triggerProcessName


#####################
### CMSSW configs ###
#####################
import FWCore.ParameterSet.Config as cms
process = cms.Process('B0KSTMUMUNTUPLIZER')

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.MessageLogger.suppressWarning = cms.untracked.vstring('B0KstMuMu')


#################
### GlobalTag ###
#################
if (runDataMC == 1):
#     if (run2012not2011 == True):
#     process.GlobalTag.globaltag = cms.string('80X_dataRun2_2016SeptRepro_v7') 
    process.GlobalTag.globaltag = cms.string('80X_dataRun2_Prompt_v16')
#     else:
#         process.GlobalTag.globaltag = cms.string('START53_LV6::All')
else:
#     if (run2012not2011 == True):
    process.GlobalTag.globaltag = cms.string('80X_mcRun2_asymptotic_2016_TrancheIV_v6')
#     else:
#         process.GlobalTag.globaltag = cms.string('FT_53_LV5_AN1::All')


#######################
### Read input file ###
#######################
import sys
if (len(sys.argv) > 2):
    readFiles = sys.argv[2]
### GEN MC ###
#    path = 'file:/mnt/hadoop/store/user/dinardo/PYTHIA6_Bd2JpsiKstar_EtaPtFilter_TuneZ2star_8TeV_GEN_NoFilter/'
#    path = 'file:/mnt/hadoop/store/user/dinardo/PYTHIA6_Bd2Psi2SKstar_EtaPtFilter_TuneZ2star_8TeV_GEN_NoFilter/'
    file = readFiles.replace(path, '')
else:
### GEN MC ###
    from B0ToKstMuMu_GEN_NoFilter_MC_cff import readFiles
#    from B0ToJPsiKst_GEN_NoFilter_MC_cff import readFiles
#    from B0ToPsi2SKst_GEN_NoFilter_MC_cff import readFiles


process.source = cms.Source('PoolSource',
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            skipEvents = cms.untracked.uint32(0),
                            fileNames = cms.untracked.vstring(
                              #readFiles
                              'file:/gwpool/users/fiorendi/p5prime/CMSSW_8_0_24/src/B0KstarMM/B0KstMuMu/plugins/0053C956-AC8E-E611-B7B8-003048FF9ABC.root'
#                               'file:/gwpool/users/fiorendi/p5prime/CMSSW_8_0_24/src/B0KstarMM/B0KstMuMu/563B1859-919C-E611-ABCC-848F69FD457A.root'
#                               'root://cms-xrd-global.cern.ch//store/data/Run2016G/Charmonium/AOD/23Sep2016-v1/110000/0E66894A-E39B-E611-A30C-848F69FD45CB.root'
                             #'/store/data/Run2016G/DoubleMuonLowMass/AOD/23Sep2016-v1/100000/0080D7F4-928D-E611-B8F8-0CC47AD99044.root'

                            ))


##################################
### Use JSON file interatively ###
##################################
if (runDataMC == 1 and useJSON == True):
    import FWCore.PythonUtilities.LumiList as LumiList
    import FWCore.ParameterSet.Types as CfgTypes
    if (run2016 == True):
        myLumis = LumiList.LumiList(filename = 'Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_MuonPhys.txt').getCMSSWString().split(',')
    else:
        myLumis = LumiList.LumiList(filename = 'Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_MuonPhys.txt').getCMSSWString().split(',')
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    process.source.lumisToProcess.extend(myLumis)


###################
### Output file ###
###################
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.maxLuminosityBlocks = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.TFileService = cms.Service('TFileService', fileName = cms.string(
    'B0ToKstMuMu.root'
### GEN MC ###
#    '/nfs/data36/cms/dinardo/B0ToKstMuMu_GEN_NoFilter_MC_NTuples/' + file
#    '/nfs/data36/cms/dinardo/B0ToJPsiKst_GEN_NoFilter_MC_NTuples/' + file
#    '/nfs/data36/cms/dinardo/B0ToPsi2SKst_GEN_NoFilter_MC_NTuples/' + file
    ), closeFileFast = cms.untracked.bool(True))


###################################
### Import PAT skeleton process ###
###################################
process.load("PhysicsTools.PatAlgos.patSequences_cff")


############################
### Add track candidates ###
############################
process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
# Require NOT to check overlap with muons and electrons
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps     = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)


#########################
### Set up PAT tracks ###
#########################
from PhysicsTools.PatAlgos.tools.trackTools import *
from PhysicsTools.PatAlgos.tools.coreTools  import *

# import pdb ; pdb.set_trace()
makeTrackCandidates(process,
                    label        = 'TrackCands',                  # output collection
                    tracks       = cms.InputTag('generalTracks'), # input track collection
                    particleType = 'pi+',                         # particle type (for assigning a mass)
                    preselection = 'pt > 0.1',                    # preselection cut on candidates
                    selection    = 'pt > 0.1',                    # selection cut on candidates
                    isolation    = {},                            # isolations to use (set to {} for None)
                    isoDeposits  = [],
                    mcAs         = None)                          # replicate MC match as the one used for Muons
# import pdb ; pdb.set_trace()


removeMCMatching(process, ['All'], outputModules = [])


#####################################
### Do trigger matching for muons ###
#####################################

if (run2016 == True):
    if (HLT_LMNR):
        process.cleanMuonTriggerMatchHLT = cms.EDProducer(
            'PATTriggerMatcherDRLessByR',
            src                   = cms.InputTag('cleanPatMuons'),
            matched               = cms.InputTag('patTrigger'),
            matchedCuts           = cms.string('path("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v*")'),
            maxDeltaR             = cms.double(0.1),
            resolveAmbiguities    = cms.bool(True),
            resolveByMatchQuality = cms.bool(True))
        process.cleanTrackTriggerMatchHLT = cms.EDProducer(
            'PATTriggerMatcherDRLessByR',
            src                   = cms.InputTag('cleanPatTrackCands'),
            matched               = cms.InputTag('patTrigger'),
            matchedCuts           = cms.string('path("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v*")'),
            maxDeltaR             = cms.double(0.1),
            resolveAmbiguities    = cms.bool(True),
            resolveByMatchQuality = cms.bool(True))
    elif (HLT_Jpsi):
        process.cleanMuonTriggerMatchHLT = cms.EDProducer(
            'PATTriggerMatcherDRLessByR',
            src                   = cms.InputTag('cleanPatMuons'),
            matched               = cms.InputTag('patTrigger'),
            matchedCuts           = cms.string('path("HLT_DoubleMu4_JpsiTrk_Displaced_v*")'),
            maxDeltaR             = cms.double(0.1),
            resolveAmbiguities    = cms.bool(True),
            resolveByMatchQuality = cms.bool(True))
        process.cleanTrackTriggerMatchHLT = cms.EDProducer(
            'PATTriggerMatcherDRLessByR',
            src                   = cms.InputTag('cleanPatTrackCands'),
            matched               = cms.InputTag('patTrigger'),
            matchedCuts           = cms.string('path("HLT_DoubleMu4_JpsiTrk_Displaced_v*")'),
            maxDeltaR             = cms.double(0.1),
            resolveAmbiguities    = cms.bool(True),
            resolveByMatchQuality = cms.bool(True))
    elif (HLT_PsiPrime):
        process.cleanMuonTriggerMatchHLT = cms.EDProducer(
            'PATTriggerMatcherDRLessByR',
            src                   = cms.InputTag('cleanPatMuons'),
            matched               = cms.InputTag('patTrigger'),
            matchedCuts           = cms.string('path("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v*")'),
            maxDeltaR             = cms.double(0.1),
            resolveAmbiguities    = cms.bool(True),
            resolveByMatchQuality = cms.bool(True))

        process.cleanTrackTriggerMatchHLT = cms.EDProducer(
            'PATTriggerMatcherDRLessByR',
            src                   = cms.InputTag('cleanPatTrackCands'),
            matched               = cms.InputTag('patTrigger'),
            matchedCuts           = cms.string('path("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v*")'),
            maxDeltaR             = cms.double(0.1),
            resolveAmbiguities    = cms.bool(True),
            resolveByMatchQuality = cms.bool(True))
else:
    process.cleanMuonTriggerMatchHLT1 = cms.EDProducer(
        'PATTriggerMatcherDRLessByR',
        src                   = cms.InputTag('cleanPatMuons'),
        matched               = cms.InputTag('patTrigger'),
        matchedCuts           = cms.string('path("HLT_Dimuon7_LowMass_Displaced_v*")'),
        maxDeltaR             = cms.double(0.1),
        resolveAmbiguities    = cms.bool(True),
        resolveByMatchQuality = cms.bool(True))
    process.cleanMuonTriggerMatchHLT2 = cms.EDProducer(
        'PATTriggerMatcherDRLessByR',
        src                   = cms.InputTag('cleanPatMuons'),
        matched               = cms.InputTag('patTrigger'),
        matchedCuts           = cms.string('path("HLT_DoubleMu4_LowMass_Displaced_v*")'),
        maxDeltaR             = cms.double(0.1),
        resolveAmbiguities    = cms.bool(True),
        resolveByMatchQuality = cms.bool(True))
    process.cleanMuonTriggerMatchHLT3 = cms.EDProducer(
        'PATTriggerMatcherDRLessByR',
        src                   = cms.InputTag('cleanPatMuons'),
        matched               = cms.InputTag('patTrigger'),
        matchedCuts           = cms.string('path("HLT_DoubleMu4p5_LowMass_Displaced_v*")'),
        maxDeltaR             = cms.double(0.1),
        resolveAmbiguities    = cms.bool(True),
        resolveByMatchQuality = cms.bool(True))
    process.cleanMuonTriggerMatchHLT4 = cms.EDProducer(
        'PATTriggerMatcherDRLessByR',
        src                   = cms.InputTag('cleanPatMuons'),
        matched               = cms.InputTag('patTrigger'),
        matchedCuts           = cms.string('path("HLT_DoubleMu5_LowMass_Displaced_v*")'),
        maxDeltaR             = cms.double(0.1),
        resolveAmbiguities    = cms.bool(True),
        resolveByMatchQuality = cms.bool(True))


##################################
### Switch on PAT trigger info ###
##################################
from PhysicsTools.PatAlgos.tools.trigTools import *
if (run2016 == True):
    switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT', 'cleanTrackTriggerMatchHLT' ], hltProcess = triggerProcessName, outputModule = '')
else:
    switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT1','cleanMuonTriggerMatchHLT2','cleanMuonTriggerMatchHLT3','cleanMuonTriggerMatchHLT4'], hltProcess = triggerProcessName, outputModule = '')


#########################
### Remove unused PAT ###
#########################
process.patDefaultSequence.remove(process.muonMatch)
process.patDefaultSequence.remove(process.electronMatch)
process.patDefaultSequence.remove(process.photonMatch)
process.patDefaultSequence.remove(process.tauMatch)
process.patDefaultSequence.remove(process.tauGenJetMatch)

process.patDefaultSequence.remove(process.tauGenJets)
process.patDefaultSequence.remove(process.tauGenJetsSelectorAllHadrons)

process.patDefaultSequence.remove(process.patJetPartonMatch)
process.patDefaultSequence.remove(process.patJetGenJetMatch)
process.patDefaultSequence.remove(process.patJetCorrFactors)
process.patDefaultSequence.remove(process.patJetCharge)
process.patDefaultSequence.remove(process.patJetPartons)
process.patDefaultSequence.remove(process.patJetPartonsLegacy)            
process.patDefaultSequence.remove(process.patJetPartonAssociationLegacy)  
process.patDefaultSequence.remove(process.patJetFlavourAssociationLegacy) 
# process.patDefaultSequence.remove(process.patJetPartonAssociation)
process.patDefaultSequence.remove(process.patJetFlavourAssociation)
process.patDefaultSequence.remove(process.patJets)

process.patDefaultSequence.remove(process.selectedPatJets)
process.patDefaultSequence.remove(process.cleanPatJets)
process.patDefaultSequence.remove(process.countPatJets)

# process.patDefaultSequence.remove(process.patCandidateSummary)
# process.patDefaultSequence.remove(process.selectedPatCandidateSummary)
# process.patDefaultSequence.remove(process.cleanPatCandidateSummary)

process.patDefaultSequence.remove(process.patMETs)

# 
# ##########################
# ### B0 --> K*0 mu+ mu- ###
# ##########################
process.B0KstMuMu = cms.EDAnalyzer('B0KstMuMu',
                                   HLTriggerResults = cms.untracked.string(triggerProcessName),
                                   VtxSample        = cms.InputTag("offlinePrimaryVertices"),
                                   BeamSpot         = cms.InputTag("offlineBeamSpot"),
                                   GenParticles     = cms.untracked.InputTag('genParticles'),
                                   MuonType         = cms.untracked.InputTag('cleanPatMuonsTriggerMatch'),
                                   TrackType        = cms.untracked.InputTag('cleanPatTrackCandsTriggerMatch'),
                                   TriggerResult    = cms.untracked.InputTag("TriggerResults::HLT"),
                                   PuInfoTag        = cms.untracked.InputTag("addPileupInfo"),
                                   GenFilterTag     = cms.untracked.InputTag("genFilterEfficiencyProducer"),
                                   ParameterFile    = cms.untracked.string('ParameterFile.txt'),
                                   doGenReco        = cms.untracked.uint32(runDataMC),
                                   TriggerNames     = cms.vstring("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v", 
                                                                  "HLT_DoubleMu4_JpsiTrk_Displaced_v",
                                                                  "HLT_DoubleMu4_PsiPrimeTrk_Displaced_v"
                                                                  ),
                                   printMsg         = cms.untracked.bool(printMsg))
# 
# 
# ###########
# ### RUN ###
# ###########
# ### Run unscheduled = create and read just what I need ### 
# sara: necessary to have all pat ingredients
process.options = cms.untracked.PSet(allowUnscheduled = cms.untracked.bool(True))

process.patPath  = cms.Path(process.patDefaultSequence)
process.ntupPath = cms.Path(process.B0KstMuMu)

if (runDataMC == 3):
    process.schedule = cms.Schedule(process.ntupPath)
else:
    process.schedule = cms.Schedule(process.patPath, process.ntupPath)

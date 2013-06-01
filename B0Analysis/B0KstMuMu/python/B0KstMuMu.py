####################
### B0KstMuMu.py ###
####################




#################
### Variables ###
#################
runDataMC = 1 # 1 = Data; 2 = MC (Reco + Gen); 3 = MC (Gen)
useJSON   = False
PrintMsg  = False
May10Code = False
triggerProcessName = 'HLT' # 'GEN' or 'HLT' or 'RECO' or 'TEST' or 'REDIGI36X'




#####################
### CMSSW configs ###
#####################
import FWCore.ParameterSet.Config as cms

process = cms.Process('NTUPLIZER')

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.MessageLogger.suppressWarning = cms.untracked.vstring('B0Cand')


#################
### GlobalTag ###
#################
if (runDataMC != 1):
    process.GlobalTag.globaltag = cms.string('START42_V17::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_R_42_V25::All')




#######################
### Read input file ###
#######################
import sys
if (len(sys.argv) > 2):
    readFiles = sys.argv[2]
### RECO Central MC ###
#    path = 'file:/bestman/storage/cms/store/mc/Fall11/B0ToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgen/GEN-SIM-RECO/HLTBPh2011_START42_V14B-v2/00000/'
#    path = 'file:/bestman/storage/cms/store/mc/Fall11/BpToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgen/GEN-SIM-RECO/HLTBPh2011_START42_V14B-v2/00000/'
#    path = 'file:/bestman/storage/cms/store/mc/Fall11/BsToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgen/GEN-SIM-RECO/HLTBPh2011_START42_V14B-v2/00000/'
#    path = 'file:/bestman/storage/cms/store/mc/Fall11/LambdaBToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgen/GEN-SIM-RECO/HLTBPh2011_START42_V14B-v2/00000/'
### RECO Private MC ###
#    path = 'file:/nfs/data36/cms/dinardo/B0ToKstMuMu_MyMCRECO001_NTuples/'
    path = 'file:/bestman/storage/cms/store/user/drell/B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_prod2_1/B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_RECO_prod2_1/d80aca2d30d1865a7fb9254b7a4518c6/'
### GEN MC ###
#    path = 'file:/bestman/storage/cms/store/user/dinardo/PYTHIA6_B0dToKstMuMuKPi_7TeV_cff_GEN_Filter/'
#    path = 'file:/bestman/storage/cms/store/user/dinardo/PYTHIA6_B0dToKstMuMuKPi_7TeV_cff_GEN_NoFilter_01/'
#    path = 'file:/bestman/storage/cms/store/user/dinardo/PYTHIA6_B0dToJPsiKstMuMuKPi_7TeV_cff_GEN_Filter/'
#    path = 'file:/bestman/storage/cms/store/user/dinardo/PYTHIA6_B0dToJPsiKstMuMuKPi_7TeV_cff_GEN_NoFilter/'
#    path = 'file:/bestman/storage/cms/store/user/dinardo/PYTHIA6_B0dToPsi2SKstMuMuKPi_7TeV_cff_GEN_Filter/'
#    path = 'file:/bestman/storage/cms/store/user/dinardo/PYTHIA6_B0dToPsi2SKstMuMuKPi_7TeV_cff_GEN_NoFilter/'
    file = readFiles.replace(path, '')
else:
### RECO Central MC ###
#    from B0ToPsiMuMu_MC_cff import readFiles
#    from BpToPsiMuMu_MC_cff import readFiles
#    from BsToPsiMuMu_MC_cff import readFiles
#    from LambdaBToPsiMuMu_MC_cff import readFiles
### RECO Private MC ###
#    from B0ToKstMuMu_MyMCRECO001_cff import readFiles
    from B0ToKstMuMu_BrMC1_cff import readFiles
### GEN MC ###
#    from B0ToKstMuMu_GEN_Filter_MC_cff import readFiles
#    from B0ToKstMuMu_GEN_NoFilter_01_MC_cff import readFiles
#    from B0ToJPsiKst_GEN_Filter_MC_cff import readFiles
#    from B0ToJPsiKst_GEN_NoFilter_MC_cff import readFiles
#    from B0ToPsi2SKst_GEN_Filter_MC_cff import readFiles
#    from B0ToPsi2SKst_GEN_NoFilter_MC_cff import readFiles
process.source = cms.Source('PoolSource',
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            skipEvents = cms.untracked.uint32(0),
                            fileNames = cms.untracked.vstring(
#    'rfio:/castor/cern.ch/user/d/dinardo/MyJPsiKstpMC/MyJPsiKstpSkim_1_1_yR4.root'
    readFiles
    ))


##################################
### Use JSON file interatively ###
##################################
if (runDataMC == 1 and useJSON == True):
    import PhysicsTools.PythonAnalysis.LumiList as LumiList
    import FWCore.ParameterSet.Types as CfgTypes
    myLumis = LumiList.LumiList(filename = 'Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt').getCMSSWString().split(',')
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    process.source.lumisToProcess.extend(myLumis)


###################
### Output file ###
###################
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.TFileService = cms.Service('TFileService', fileName = cms.string(
    'B0ToKstMuMu.root'
### RECO Central MC ###
#    '/nfs/data36/cms/dinardo/B0ToPsiMuMu_MC_NTuples/' + readFiles[len(readFiles) - 41 : len(readFiles)]
#    '/nfs/data36/cms/dinardo/BpToPsiMuMu_MC_NTuples/' + readFiles[len(readFiles) - 41 : len(readFiles)]
#    '/nfs/data36/cms/dinardo/BsToPsiMuMu_MC_NTuples/' + readFiles[len(readFiles) - 41 : len(readFiles)]
#    '/nfs/data36/cms/dinardo/LambdaBToPsiMuMu_MC_NTuples/' + readFiles[len(readFiles) - 41 : len(readFiles)]
### RECO Private MC ###
#    '/nfs/data36/cms/dinardo/B0ToKstMuMu_MyMCANALYSIS001_NTuples/' + file
#    '/nfs/data36/cms/dinardo/B0ToKstMuMu_BrMC1_NTuples/' + file
### GEN MC ###
#    '/nfs/data36/cms/dinardo/B0ToKstMuMu_GEN_Filter_MC_NTuples/' + file
#    '/nfs/data36/cms/dinardo/B0ToKstMuMu_GEN_NoFilter_01_MC_NTuples/' + file
#    '/nfs/data36/cms/dinardo/B0ToJPsiKst_GEN_Filter_MC_NTuples/' + file
#    '/nfs/data36/cms/dinardo/B0ToJPsiKst_GEN_NoFilter_MC_NTuples/' + file
#    '/nfs/data36/cms/dinardo/B0ToPsi2SKst_GEN_Filter_MC_NTuples/' + file
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
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)

from PhysicsTools.PatAlgos.tools.trackTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
if (runDataMC != 1):
    makeTrackCandidates(process,
                        label        = 'TrackCands',                  # output collection
                        tracks       = cms.InputTag('generalTracks'), # input track collection
                        particleType = 'pi+',                         # particle type (for assigning a mass)
                        preselection = 'pt > 0.1',                    # preselection cut on candidates
                        selection    = 'pt > 0.1',                    # selection on PAT Layer 1 objects
                        isolation    = {},                            # isolations to use (set to {} for None)
                        isoDeposits  = [],
                        mcAs         = 'muon'                         # replicate MC match as the one used for Muons
                        )
    
    process.patTrackCandsMCMatch.mcPdgId               = cms.vint32(211) # = pi+
    process.patTrackCandsMCMatch.mcStatus              = cms.vint32(1)
    process.patTrackCandsMCMatch.maxDeltaR             = cms.double(0.02)
    process.patTrackCandsMCMatch.resolveByMatchQuality = cms.bool(True)
else:
    makeTrackCandidates(process,
                        label        = 'TrackCands',                  # output collection
                        tracks       = cms.InputTag('generalTracks'), # input track collection
                        particleType = 'pi+',                         # particle type (for assigning a mass)
                        preselection = 'pt > 0.1',                    # preselection cut on candidates
                        selection    = 'pt > 0.1',                    # selection on PAT Layer 1 objects
                        isolation    = {},                            # isolations to use (set to {} for None)
                        isoDeposits  = [],
                        mcAs         = None                           # replicate MC match as the one used for Muons
                        )
    
    removeMCMatching(process, ['All'], outputInProcess = False)


#####################################
### Do trigger matching for muons ###
#####################################
# For May10 ReReco only
process.cleanMuonTriggerMatchHLT0 = cms.EDProducer(
    'PATTriggerMatcherDRLessByR',                          # match by DeltaR only (best match by DeltaR)
    src                   = cms.InputTag('cleanPatMuons'),
    matched               = cms.InputTag('patTrigger'),    # default producer label as defined in PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
    matchedCuts           = cms.string('path("HLT_Dimuon6p5_LowMass_Displaced_v*",0,0)'),
    maxDeltaR             = cms.double(0.1),
    resolveAmbiguities    = cms.bool(True),                # only one match per trigger object
    resolveByMatchQuality = cms.bool(False))               # take best match found per reco object (by DeltaR here, see above)
# After May10 ReReco
process.cleanMuonTriggerMatchHLT1 = cms.EDProducer(
    'PATTriggerMatcherDRLessByR',
    src                   = cms.InputTag('cleanPatMuons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('path("HLT_Dimuon7_LowMass_Displaced_v*")'),
    maxDeltaR             = cms.double(0.1),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(False))
process.cleanMuonTriggerMatchHLT2 = cms.EDProducer(
    'PATTriggerMatcherDRLessByR',
    src                   = cms.InputTag('cleanPatMuons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('path("HLT_DoubleMu4_LowMass_Displaced_v*")'),
    maxDeltaR             = cms.double(0.1),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(False))
process.cleanMuonTriggerMatchHLT3 = cms.EDProducer(
    'PATTriggerMatcherDRLessByR',
    src                   = cms.InputTag('cleanPatMuons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('path("HLT_DoubleMu4p5_LowMass_Displaced_v*")'),
    maxDeltaR             = cms.double(0.1),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(False))
process.cleanMuonTriggerMatchHLT4 = cms.EDProducer(
    'PATTriggerMatcherDRLessByR',
    src                   = cms.InputTag('cleanPatMuons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('path("HLT_DoubleMu5_LowMass_Displaced_v*")'),
    maxDeltaR             = cms.double(0.1),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(False))


##################################
### Switch on PAT trigger info ###
##################################
from PhysicsTools.PatAlgos.tools.trigTools import *
if (May10Code == True):
    # For May10 ReReco only
    switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT0'], hltProcess = triggerProcessName, outputModule = '')
else:
    # After May10 ReReco
    switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT1','cleanMuonTriggerMatchHLT2','cleanMuonTriggerMatchHLT3','cleanMuonTriggerMatchHLT4'], hltProcess = triggerProcessName, outputModule = '')


################################
### Remove not used from PAT ###
################################
process.patDefaultSequence.remove(process.patJetCorrFactors)
process.patDefaultSequence.remove(process.patJetCharge)
process.patDefaultSequence.remove(process.patJetPartonMatch)
process.patDefaultSequence.remove(process.patJetGenJetMatch)
process.patDefaultSequence.remove(process.patJetPartons)
process.patDefaultSequence.remove(process.patJetPartonAssociation)
process.patDefaultSequence.remove(process.patJetFlavourAssociation)
process.patDefaultSequence.remove(process.patJets)
process.patDefaultSequence.remove(process.metJESCorAK5CaloJet)
process.patDefaultSequence.remove(process.metJESCorAK5CaloJetMuons)
process.patDefaultSequence.remove(process.patMETs)
process.patDefaultSequence.remove(process.selectedPatJets)
process.patDefaultSequence.remove(process.cleanPatJets)
process.patDefaultSequence.remove(process.countPatJets)




##########################
### B0 --> K*0 mu+ mu- ###
##########################
process.B0Cand = cms.EDAnalyzer('B0KstMuMu',
                                HLTriggerResults = cms.untracked.string(triggerProcessName),
                                VtxSample        = cms.untracked.string('offlinePrimaryVertices'),
                                BeamSpot         = cms.untracked.string('offlineBeamSpot'),
                                GenParticles     = cms.untracked.string('genParticles'),
                                MuonType         = cms.untracked.string('cleanPatMuonsTriggerMatch'),
                                TrackType        = cms.untracked.string('cleanPatTrackCands'),
                                ParameterFile    = cms.untracked.string('ParameterFile.txt'),
                                doGenReco        = cms.untracked.uint32(runDataMC),
                                printMsg         = cms.untracked.bool(PrintMsg))


###########
### RUN ###
###########
process.patPath  = cms.Path(process.patDefaultSequence)
process.ntupPath = cms.Path(process.B0Cand)

if (runDataMC == 3):
    process.schedule = cms.Schedule(process.ntupPath)
else:
    process.schedule = cms.Schedule(process.patPath, process.ntupPath)

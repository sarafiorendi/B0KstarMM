#################
### Variables ###
#################
runDataMC = 1 # 1 = Data; 2 = MC (Reco + Gen); 3 = MC (Gen)
useJSON   = False
printMsg  = False
triggerProcessName = 'HLT' # 'GEN' or 'HLT' or 'RECO' or 'TEST' or ...


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
if (runDataMC != 1):
    process.GlobalTag.globaltag = cms.string('START53_V7G::All') # Signal MC: START53_V19F; J/psi X MC: START53_V7G
else:
    process.GlobalTag.globaltag = cms.string('FT53_V21A_AN6::All')


#######################
### Read input file ###
#######################
import sys
if (len(sys.argv) > 2):
    readFiles = sys.argv[2]
### GEN MC ###
    path = 'file:/bestman/storage/cms/store/user/dinardo/PYTHIA6_B0dToKstMuMuKPi_7TeV_cff_GEN_Filter/'
#    path = 'file:/bestman/storage/cms/store/user/dinardo/PYTHIA6_B0dToKstMuMuKPi_7TeV_cff_GEN_NoFilter_01/'

#    path = 'file:/bestman/storage/cms/store/user/dinardo/PYTHIA6_B0dToJPsiKstMuMuKPi_7TeV_cff_GEN_Filter/'
#    path = 'file:/bestman/storage/cms/store/user/dinardo/PYTHIA6_B0dToJPsiKstMuMuKPi_7TeV_cff_GEN_NoFilter/'

#    path = 'file:/bestman/storage/cms/store/user/dinardo/PYTHIA6_B0dToPsi2SKstMuMuKPi_7TeV_cff_GEN_Filter/'
#    path = 'file:/bestman/storage/cms/store/user/dinardo/PYTHIA6_B0dToPsi2SKstMuMuKPi_7TeV_cff_GEN_NoFilter/'
    file = readFiles.replace(path, '')
else:
### GEN MC ###
    from B0ToKstMuMu_GEN_Filter_MC_cff import readFiles
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
    import FWCore.PythonUtilities.LumiList as LumiList
    import FWCore.ParameterSet.Types as CfgTypes
    myLumis = LumiList.LumiList(filename = 'Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON_MuonPhys.txt').getCMSSWString().split(',')
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    process.source.lumisToProcess.extend(myLumis)


###################
### Output file ###
###################
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.TFileService = cms.Service('TFileService', fileName = cms.string(
    'B0ToKstMuMu.root'
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
# Require NOT to check overlap with muons and electrons
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps     = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)

from PhysicsTools.PatAlgos.tools.trackTools import *
from PhysicsTools.PatAlgos.tools.coreTools  import *
makeTrackCandidates(process,
                    label        = 'TrackCands',                  # output collection
                    tracks       = cms.InputTag('generalTracks'), # input track collection
                    particleType = 'pi+',                         # particle type (for assigning a mass)
                    preselection = 'pt > 0.1',                    # preselection cut on candidates
                    selection    = 'pt > 0.1',                    # selection cut on candidates
                    isolation    = {},                            # isolations to use (set to {} for None)
                    isoDeposits  = [],
                    mcAs         = None                           # replicate MC match as the one used for Muons
                )
    
removeMCMatching(process, ['All'], outputModules = [])


#####################################
### Do trigger matching for muons ###
#####################################
process.cleanMuonTriggerMatchHLT = cms.EDProducer(
    'PATTriggerMatcherDRLessByR',
    src                   = cms.InputTag('cleanPatMuons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('path("HLT_DoubleMu3p5_LowMass_Displaced_v*")'),
    maxDeltaR             = cms.double(0.1),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(True))


##################################
### Switch on PAT trigger info ###
##################################
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT'], hltProcess = triggerProcessName, outputModule = '')


###########################
### Remove not used PAT ###
###########################
process.patDefaultSequence.remove(process.patJetCorrFactors)
process.patDefaultSequence.remove(process.patJetCharge)
process.patDefaultSequence.remove(process.patJetPartonMatch)
process.patDefaultSequence.remove(process.patJetGenJetMatch)
process.patDefaultSequence.remove(process.patJetPartons)
process.patDefaultSequence.remove(process.patJetPartonAssociation)
process.patDefaultSequence.remove(process.patJetFlavourAssociation)
process.patDefaultSequence.remove(process.patJets)
process.patDefaultSequence.remove(process.patMETs)
process.patDefaultSequence.remove(process.selectedPatJets)
process.patDefaultSequence.remove(process.cleanPatJets)
process.patDefaultSequence.remove(process.countPatJets)


##########################
### B0 --> K*0 mu+ mu- ###
##########################
process.B0KstMuMu = cms.EDAnalyzer('B0KstMuMu',
                                   HLTriggerResults = cms.untracked.string(triggerProcessName),
                                   VtxSample        = cms.untracked.string('offlinePrimaryVertices'),
                                   BeamSpot         = cms.untracked.string('offlineBeamSpot'),
                                   GenParticles     = cms.untracked.string('genParticles'),
                                   MuonType         = cms.untracked.string('cleanPatMuonsTriggerMatch'),
                                   TrackType        = cms.untracked.string('cleanPatTrackCands'),
                                   ParameterFile    = cms.untracked.string('ParameterFile.txt'),
                                   doGenReco        = cms.untracked.uint32(runDataMC),
                                   printMsg         = cms.untracked.bool(printMsg))


###########
### RUN ###
###########
process.patPath  = cms.Path(process.patDefaultSequence)
process.ntupPath = cms.Path(process.B0KstMuMu)

if (runDataMC == 3):
    process.schedule = cms.Schedule(process.ntupPath)
else:
    process.schedule = cms.Schedule(process.patPath, process.ntupPath)

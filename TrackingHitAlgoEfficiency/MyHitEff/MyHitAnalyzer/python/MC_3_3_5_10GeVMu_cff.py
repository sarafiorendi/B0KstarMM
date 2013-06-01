import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
        "/store/relval/CMSSW_3_3_5/RelValSingleMuPt10/GEN-SIM-RECO/MC_31X_V9-v1/0007/20587239-CADB-DE11-B2C3-002618943920.root",
        "/store/relval/CMSSW_3_3_5/RelValSingleMuPt10/GEN-SIM-RECO/MC_31X_V9-v1/0007/9C7BAD23-C9DB-DE11-B43D-00304867902C.root",
        "/store/relval/CMSSW_3_3_5/RelValSingleMuPt10/GEN-SIM-RECO/MC_31X_V9-v1/0007/CE3E5A35-CADB-DE11-9C41-003048678FEA.root" ] )

secFiles.extend( [
        "/store/relval/CMSSW_3_3_5/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_31X_V9-v1/0007/24CDE520-C9DB-DE11-A1CC-003048678ED2.root",
        "/store/relval/CMSSW_3_3_5/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_31X_V9-v1/0007/3A768D1E-C9DB-DE11-8A8E-002618943910.root",
        "/store/relval/CMSSW_3_3_5/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_31X_V9-v1/0007/3AD3AF2F-CADB-DE11-A69C-002618943964.root",
        "/store/relval/CMSSW_3_3_5/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_31X_V9-v1/0007/46B25B34-CADB-DE11-BDFE-003048678FB8.root",
        "/store/relval/CMSSW_3_3_5/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_31X_V9-v1/0007/5EB8CF33-CADB-DE11-A68C-003048678A80.root",
        "/store/relval/CMSSW_3_3_5/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_31X_V9-v1/0007/6A326A23-C9DB-DE11-8C35-003048678BC6.root",
        "/store/relval/CMSSW_3_3_5/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_31X_V9-v1/0007/94B5D52C-CADB-DE11-8A83-003048678D52.root",
        "/store/relval/CMSSW_3_3_5/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_31X_V9-v1/0007/D0DBBC33-CADB-DE11-9C01-003048678A80.root" ] )

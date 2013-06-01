################
# Instructions #
################

# Checkout from cvs HEAD Configuration/GenProduction/python/customise_SilentMessageLogger.py
# Copy PYTHIA6_B0dToKstMuMuKPi_7TeV_cff.py in Configuration/GenProduction/python/
# If not in cvs copy Bd_MuMuKstar_mumuKpi.dec in GeneratorInterface/ExternalDecays/data/
# scram b
# source GenMyMC.sh


# Create just GEN-MC
# Remove only the mumugenfilter from PYTHIA6_B0dToKstMuMuKPi_7TeV_cff.py
#cmsDriver.py Configuration/GenProduction/python/PYTHIA6_B0dToKstMuMuKPi_7TeV_cff.py -s GEN --conditions START42_V17::All --datatier GEN-SIM --eventcontent RAWSIM --customise Configuration/GenProduction/customise_SilentMessageLogger.py -n 10000000000 --no_exec

# Create just GEN-SIM-MC for officila production tests
cmsDriver.py Configuration/GenProduction/python/PYTHIA6_B0dToPsi2SKstMuMuKPi_7TeV_cff.py -s GEN,SIM --conditions START42_V17::All --datatier GEN-SIM --eventcontent RAWSIM --customise Configuration/GenProduction/customise_SilentMessageLogger.py -n 100000 --no_exec

# Create MC up to HLT
#cmsDriver.py Configuration/GenProduction/python/PYTHIA6_B0dToKstMuMuKPi_7TeV_cff.py -s GEN,SIM,DIGI,L1,DIGI2RAW,HLT --conditions START42_V17::All --datatier GEN-SIM-RAW --eventcontent RAWSIM --customise Configuration/GenProduction/customise_SilentMessageLogger.py -n 10000000000 --no_exec

# Create RECO-MC from MC up to HLT (output RECO + SIM)
#cmsDriver.py Configuration/GenProduction/python/PYTHIA6_B0dToKstMuMuKPi_7TeV_cff.py -s RAW2DIGI,RECO --conditions START42_V17::All --datatier GEN-SIM-RECO --eventcontent RECOSIM --filein file:PYTHIA6_B0dToKstMuMuKPi_7TeV_cff_GEN_SIM_DIGI_L1_DIGI2RAW_HLT.root --customise Configuration/GenProduction/customise_SilentMessageLogger.py -n -1 --no_exec

# Create RECO-MC from MC up to HLT (output AOD + SIM)
#cmsDriver.py Configuration/GenProduction/python/PYTHIA6_B0dToKstMuMuKPi_7TeV_cff.py -s RAW2DIGI,RECO --conditions START42_V17::All --datatier GEN-SIM-RECO --eventcontent AODSIM --filein file:PYTHIA6_B0dToKstMuMuKPi_7TeV_cff_GEN_SIM_DIGI_L1_DIGI2RAW_HLT.root --customise Configuration/GenProduction/customise_SilentMessageLogger.py -n -1 --no_exec

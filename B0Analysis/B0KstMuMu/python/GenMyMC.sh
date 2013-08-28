################
# Instructions #
################

# 1. Checkout from repository HEAD Configuration/Generator/python/customise_SilentMessageLogger.py
# 2. Copy PYTHIA6_B0dToKstMuMuKPi_7TeV_cff.py in Configuration/Generator/python/
# 3. If not in repository copy Bd_MuMuKstar_mumuKpi.dec in GeneratorInterface/ExternalDecays/data/
# 4. scram b
# 5. source GenMyMC.sh


### Create just GEN-MC ###
# Remove only the mumugenfilter from PYTHIA6_B0dToKstMuMuKPi_7TeV_cff.py
#cmsDriver.py Configuration/Generator/python/PYTHIA6_B0dToKstMuMuKPi_7TeV_cff.py -s GEN --conditions START53_V19F::All --datatier GEN-SIM --eventcontent RAWSIM --customise Configuration/Generator/customise_SilentMessageLogger.py -n 10000000000 --no_exec


### Create just GEN-SIM-MC for official production tests ###
cmsDriver.py Configuration/Generator/python/PYTHIA6_B0dToPsi2SKstMuMuKPi_7TeV_cff.py -s GEN,SIM --conditions START53_V19F::All --datatier GEN-SIM --eventcontent RAWSIM --customise Configuration/Generator/customise_SilentMessageLogger.py -n 100000 --no_exec


### Create MC up to HLT ###
#cmsDriver.py Configuration/Generator/python/PYTHIA6_B0dToKstMuMuKPi_7TeV_cff.py -s GEN,SIM,DIGI,L1,DIGI2RAW,HLT --conditions START53_V19F::All --datatier GEN-SIM-RAW --eventcontent RAWSIM --customise Configuration/Generator/customise_SilentMessageLogger.py -n 10000000000 --no_exec


### Create RECO-MC from MC up to HLT (output RECO + SIM) ###
#cmsDriver.py Configuration/Generator/python/PYTHIA6_B0dToKstMuMuKPi_7TeV_cff.py -s RAW2DIGI,RECO --conditions START53_V19F::All --datatier GEN-SIM-RECO --eventcontent RECOSIM --filein file:PYTHIA6_B0dToKstMuMuKPi_7TeV_cff_GEN_SIM_DIGI_L1_DIGI2RAW_HLT.root --customise Configuration/Generator/customise_SilentMessageLogger.py -n -1 --no_exec


### Create RECO-MC from MC up to HLT (output AOD + SIM) ###
#cmsDriver.py Configuration/Generator/python/PYTHIA6_B0dToKstMuMuKPi_7TeV_cff.py -s RAW2DIGI,RECO --conditions START53_V19F::All --datatier GEN-SIM-RECO --eventcontent AODSIM --filein file:PYTHIA6_B0dToKstMuMuKPi_7TeV_cff_GEN_SIM_DIGI_L1_DIGI2RAW_HLT.root --customise Configuration/Generator/customise_SilentMessageLogger.py -n -1 --no_exec

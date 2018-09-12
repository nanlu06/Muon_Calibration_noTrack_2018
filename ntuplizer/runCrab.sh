#cmsenv
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init --voms cms

crab submit -c scrabconfig_MC_Zmumu_NoPU.py
#crab submit -c scrabconfig_MC_Zmumu_wPU.py

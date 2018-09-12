from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'Zmm_wPU_ieta27_Sep3'
config.General.workArea = 'crab_prod'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'hcalHBHEMuon_cfg.py'

config.Data.inputDataset = '/DYToMuMu_M-20_13TeV_pythia8/RunIISpring18DRPremix-100X_upgrade2018_realistic_v10-v1/GEN-SIM-RECO'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/group/dpg_hcal/comm_hcal/nlu'
config.Data.publication = True
config.Data.outputDatasetTag = 'Zmm_wPU_ieta27_Sep3'
config.Site.blacklist = ['T2_IT_Legnaro']
config.Site.storageSite = 'T2_CH_CERN'

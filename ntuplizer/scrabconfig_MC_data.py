from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'Run2018C_v3_Ieta27_Sep17'
config.General.workArea = 'crab_prod'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'hcalHBHEMuon_cfg.py'

config.Data.inputDataset = '/SingleMuon/Run2018C-HcalCalHBHEMuonFilter-PromptReco-v3/ALCARECO'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-322057_13TeV_PromptReco_Collisions18_JSON.txt'
config.Data.outLFNDirBase = '/store/group/dpg_hcal/comm_hcal/nlu'
config.Data.publication = True
config.Data.outputDatasetTag = 'Run2018C_v3_Ieta27_Sep17'
config.Site.blacklist = ['T2_IT_Legnaro']
config.Site.storageSite = 'T2_CH_CERN'

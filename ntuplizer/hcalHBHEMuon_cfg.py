import FWCore.ParameterSet.Config as cms

process = cms.Process("RaddamMuon")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")  
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("RecoJets.Configuration.CaloTowersES_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag='100X_upgrade2018_realistic_v10'
#process.GlobalTag.globaltag='101X_dataRun2_Prompt_v11'

process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")
process.load("Calibration.HcalCalibAlgos.hcalHBHEMuon_cfi")

#if 'MessageLogger' in process.__dict__:
#    process.MessageLogger.categories.append('HBHEMuon')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch///store/data/Run2018C/SingleMuon/ALCARECO/HcalCalHBHEMuonFilter-PromptReco-v3/000/319/849/00000/B4EAF9AC-488C-E811-AC9F-FA163E31F891.root'
        #'root://cms-xrd-global.cern.ch///store/data/Run2018C/DoubleMuonLowPU/RAW-RECO/LogError-PromptReco-v2/000/319/462/00000/D8D67FB5-8C86-E811-BB68-FA163E106838.root'
        #'root://cms-xrd-global.cern.ch///store/data/Run2018C/DoubleMuonLowPU/ALCARECO/MuAlZMuMu-PromptReco-v2/000/319/488/00000/A414E3F8-C187-E811-87C6-FA163E75C6E8.root'
        #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17DRStdmix/SingleMuon_13_30_HEP17/GEN-SIM-RECO/NoPU_94X_mc2017_realistic_v11-v2/60000/0C2598E5-BB32-E811-8606-FA163EBF2E80.root'
        #'root://cms-xrd-global.cern.ch//store/data/Run2017B/SingleMuon/RECO/PromptReco-v2/000/298/678/00000/C0C0C0B0-A466-E711-AE46-02163E019E8C.root',
#       'file:/afs/cern.ch/work/a/amkalsi/public/ForSunandaDa/C0C0C0B0-A466-E711-AE46-02163E019E8C.root',
#       'root://xrootd.unl.edu//store/data/Run2017B/SingleMuon/RECO/PromptReco-v2/000/298/678/00000/C0C0C0B0-A466-E711-AE46-02163E019E8C.root'
#        'root://xrootd.unl.edu//store/mc/Phys14DR/DYToMuMu_M-50_Tune4C_13TeV-pythia8/GEN-SIM-RECO/PU20bx25_tsg_castor_PHYS14_25_V1-v1/10000/184C1AC9-A775-E411-9196-002590200824.root'
        )
                            )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("Validation.root")
)

process.hcalTopologyIdeal.MergePosition = False
process.hcalHBHEMuon.useRaw = 0
process.hcalHBHEMuon.unCorrect = True
process.hcalHBHEMuon.getCharge = True
process.hcalHBHEMuon.collapseDepth = False
process.hcalHBHEMuon.isItPlan1 = True
process.hcalHBHEMuon.ignoreHECorr = False
process.hcalHBHEMuon.isItPreRecHit = True
process.hcalHBHEMuon.maxDepth = 7
#process.hcalHBHEMuon.LabelHBHERecHit = cms.InputTag("hbheprereco")
process.hcalHBHEMuon.verbosity = 0

process.p = cms.Path(process.hcalHBHEMuon)

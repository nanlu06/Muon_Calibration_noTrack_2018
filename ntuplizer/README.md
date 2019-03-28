
copy /afs/cern.ch/work/n/nlu/public/HCal/test/IsolatedParticles to the directory: CMSSW_10_1_8/src/Calibration/
Then in this directory you have HcalCalibAlgo and IsolatedParticles
put HcalHBHEMuonAnalyzer.cc file on github to CMSSW_10_1_8/src/Calibration/HcalCalibAlgos/plugins/
put hcalHBHEMuon_cfg.py to CMSSW_10_1_8/src/Calibration/HcalCalibAlgos/test/
then compile the code.

How to run:
1) to run the code locally: cmsRun hcalHBHEMuon_cfg.py
2) to submit crab job: 

1)) put runCrab.sh  scrabconfig_MC_Zmumu_NoPU.py scrabconfig_MC_Zmumu_wPU.py in CMSSW_10_1_8/src/Calibration/HcalCalibAlgos/test/
2)) create a folder: crab_prod, in the directory CMSSW_10_1_8/src/Calibration/HcalCalibAlgos/test/ 
and then ./runCrab.sh

resubmit failed jobs: resubmit.sh


March 28, 2019:

HcalHBHEMuonAnalyzer_March282019.cc

rename it to HcalHBHEMuonAnalyzer.cc if you want to use it.
The cuts are with new veto Zmm cuts to remove more Zmm bkg. The new cuts are:

Additional muon is present which is oppositely charged w.r.t. the tag muon,  pt > 20 GeV, |η| < 2.5 but not in |ieta25|, Loose ID (by default, all nanoAOD muons are loose muons), Rel Iso <0.25 and |dimuon_mass-91.18| < 4.

old cuts are:

Additional muon is present which is oppositely charged w.r.t. the tag muon,  pt > 20 GeV, |η| < 2.5 but not in |ieta25|, Medium ID, Rel Iso <0.15 and |dimuon_mass-91.18| < 2.

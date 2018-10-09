output='/afs/cern.ch/user/n/nlu/work/private/CMS/HCal/MIPCali/CMSSW_9_2_6/src/Muon_Calibration_noTrack_2018/muon_calibration_ana/bash/MuonCali_test'
#mkdir $output
cd $output

./readMC $1 $2 -1 $3 > out_${1}

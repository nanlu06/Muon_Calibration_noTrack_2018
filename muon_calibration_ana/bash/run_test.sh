output='/afs/cern.ch/work/n/nlu/private/CMS/HCal/MIPCali/CMSSW_9_2_6/src/Analyzer/MakeTree/MakeTree_NoTrack/bash/MuonCali_test'
#mkdir $output
cd $output
cp /afs/cern.ch/work/n/nlu/private/CMS/HCal/MIPCali/CMSSW_9_2_6/src/Analyzer/MakeTree/MakeTree_NoTrack/test_${1}.list .
./readMC test_${1}.list $1 -1 /eos/cms/store/group/dpg_hcal/comm_hcal/nlu/ntuples/Sept12018/HCALSample > out_${1}

source /cvmfs/cms.cern.ch/crab3/crab.sh
path=/afs/cern.ch/user/n/nlu/work/private/CMS/HCal/MIPCali_NoTrk/CMSSW_10_1_8/src/Calibration/HcalCalibAlgos/test/crab_prod
declare -a arr=(
"crab_Run2018C-HcalCalHBHEMuonFilter-PromptReco-v2_Sep1"
)
for i in "${arr[@]}"
do
   echo "$i"
   cd $i
   crab report --dir $path/$i
   #crab kill -d $path/$i
   crab status
   crab resubmit
   crab status
   cd ../ 
done

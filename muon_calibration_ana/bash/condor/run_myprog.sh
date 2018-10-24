#!/bin/bash

#set filelist=${1}
#set outfile=${2}
#set data=${3}

#I think you need to set your root environment. Best is to checkout a CMSSW package, cmsenv will give you access to root
#cd /home/spandey/t3store/public/CMSSW/CMSSW_7_5_0_pre5/src
#source /cvmfs/cms.cern.ch/cmsset_default.sh
#eval `scram runtime -sh`
#source /cvmfs/cms.cern.ch/cmsset_default.sh
#export SCRAM_ARCH=slc7_amd64_gcc630
#ulimit -c 0
#eval `scram runtime -sh`
#echo `which root`
cd /afs/cern.ch/work/i/idutta/private/HCAL_muALCA/CMSSW_10_1_6/src/Muon_Calibration_noTrack_2018/muon_calibration_ana/bash/condor/condor_output/condor_logs
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc630
ulimit -c 0
eval `scram runtime -sh`
#env 
#echo "$TMPDIR/$2"
/afs/cern.ch/work/i/idutta/private/HCAL_muALCA/CMSSW_10_1_6/src/Muon_Calibration_noTrack_2018/muon_calibration_ana/bash/condor/readMC $1 $2 -1 $3

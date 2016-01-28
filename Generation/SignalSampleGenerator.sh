#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc481
cd /afs/cern.ch/user/k/ktos/BBA/CMSSW_7_1_11_patch2/src
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh 
mv /afs/cern.ch/user/k/ktos/BBA/CMSSW_7_1_11_patch2/src/BBA/Generation/DIRNAME/SCRIPTNAME_JOBNUM.* . 
cmsRun SCRIPTNAME_JOBNUM.py
cmsStage -f DIRNAME_JOBNUM.root /store/user/ktos/DIRNAME
rm SCRIPTNAME_JOBNUM.*  DIRNAME_JOBNUM.root*
if [ -e RandomEngineState_JOBNUM.log ]
    then
    rm RandomEngineState_JOBNUM.log*
fi
if [ -e histProbFunction_JOBNUM.root ]
    then
    rm histProbFunction_JOBNUM.root*
fi

exit 0

#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc491
cd /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp  /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/BBA/RECO/BSUB/DIRNAME/GENERATOR.py .
cmsRun GENERATOR.py
cmsStage -f DIRNAME_NUM.root /store/user/ktos/DIRNAME
rm DIRNAME_NUM.root GENERATOR.py DIRNAME_DQM_NUM.root
if [ -e RandomEngineState_NUM.log ]
    then
    rm RandomEngineState_NUM.log*
fi
if [ -e histProbFunction_NUM.root ]
    then
    rm histProbFunction_NUM.root*
fi

exit 0


#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc491
cd /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp  /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/BBA/RECO/BSUB/BBATauRecoPAT/STEP_3_cfg_BBATauRecoPAT_2.py .
cmsRun STEP_3_cfg_BBATauRecoPAT_2.py
cmsStage -f BBATauRecoPAT_2.root /store/user/ktos/BBATauRecoPAT
rm BBATauRecoPAT_2.root STEP_3_cfg_BBATauRecoPAT_2.py BBATauRecoPAT_DQM_2.root
if [ -e RandomEngineState_2.log ]
    then
    rm RandomEngineState_2.log*
fi
if [ -e histProbFunction_2.root ]
    then
    rm histProbFunction_2.root*
fi

exit 0


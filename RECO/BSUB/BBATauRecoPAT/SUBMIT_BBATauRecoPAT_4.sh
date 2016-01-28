#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc491
cd /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp  /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/BBA/RECO/BSUB/BBATauRecoPAT/STEP_3_cfg_BBATauRecoPAT_4.py .
cmsRun STEP_3_cfg_BBATauRecoPAT_4.py
cmsStage -f BBATauRecoPAT_4.root /store/user/ktos/BBATauRecoPAT
rm BBATauRecoPAT_4.root STEP_3_cfg_BBATauRecoPAT_4.py BBATauRecoPAT_DQM_4.root
if [ -e RandomEngineState_4.log ]
    then
    rm RandomEngineState_4.log*
fi
if [ -e histProbFunction_4.root ]
    then
    rm histProbFunction_4.root*
fi

exit 0


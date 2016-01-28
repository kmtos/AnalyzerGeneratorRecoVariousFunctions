#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc491
cd /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp  /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/BBA/RECO/BSUB/BBATauRecoPAT/STEP_3_cfg_BBATauRecoPAT_3.py .
cmsRun STEP_3_cfg_BBATauRecoPAT_3.py
cmsStage -f BBATauRecoPAT_3.root /store/user/ktos/BBATauRecoPAT
rm BBATauRecoPAT_3.root STEP_3_cfg_BBATauRecoPAT_3.py BBATauRecoPAT_DQM_3.root
if [ -e RandomEngineState_3.log ]
    then
    rm RandomEngineState_3.log*
fi
if [ -e histProbFunction_3.root ]
    then
    rm histProbFunction_3.root*
fi

exit 0


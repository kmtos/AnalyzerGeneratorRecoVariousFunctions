#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc491
cd /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp  /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/BBA/RECO/BSUB/ttbarTauRecoPAT/STEP_3_cfg_ttbarTauRecoPAT_12.py .
cmsRun STEP_3_cfg_ttbarTauRecoPAT_12.py
cmsStage -f ttbarTauRecoPAT_12.root /store/user/ktos/ttbarTauRecoPAT
rm ttbarTauRecoPAT_12.root STEP_3_cfg_ttbarTauRecoPAT_12.py ttbarTauRecoPAT_DQM_12.root
if [ -e RandomEngineState_12.log ]
    then
    rm RandomEngineState_12.log*
fi
if [ -e histProbFunction_12.root ]
    then
    rm histProbFunction_12.root*
fi

exit 0


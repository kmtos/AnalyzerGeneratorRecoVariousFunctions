#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc491
cd /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp  /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/BBA/RECO/BSUB/ttbarTauRecoPAT/STEP_3_cfg_ttbarTauRecoPAT_6.py .
cmsRun STEP_3_cfg_ttbarTauRecoPAT_6.py
cmsStage -f ttbarTauRecoPAT_6.root /store/user/ktos/ttbarTauRecoPAT
rm ttbarTauRecoPAT_6.root STEP_3_cfg_ttbarTauRecoPAT_6.py ttbarTauRecoPAT_DQM_6.root
if [ -e RandomEngineState_6.log ]
    then
    rm RandomEngineState_6.log*
fi
if [ -e histProbFunction_6.root ]
    then
    rm histProbFunction_6.root*
fi

exit 0


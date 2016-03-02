#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc491
cd /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp  /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/BBA/RECO/BSUB/ttbarTauRecoPAT/STEP_3_cfg_ttbarTauRecoPAT_16.py .
cmsRun STEP_3_cfg_ttbarTauRecoPAT_16.py
cmsStage -f ttbarTauRecoPAT_16.root /store/user/ktos/ttbarTauRecoPAT
rm ttbarTauRecoPAT_16.root STEP_3_cfg_ttbarTauRecoPAT_16.py ttbarTauRecoPAT_DQM_16.root
if [ -e RandomEngineState_16.log ]
    then
    rm RandomEngineState_16.log*
fi
if [ -e histProbFunction_16.root ]
    then
    rm histProbFunction_16.root*
fi

exit 0


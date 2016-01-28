#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc491
cd /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp  /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/BBA/RECO/BSUB/ttbarTauRecoPAT/STEP_3_cfg_ttbarTauRecoPAT_15.py .
cmsRun STEP_3_cfg_ttbarTauRecoPAT_15.py
cmsStage -f ttbarTauRecoPAT_15.root /store/user/ktos/ttbarTauRecoPAT
rm ttbarTauRecoPAT_15.root STEP_3_cfg_ttbarTauRecoPAT_15.py ttbarTauRecoPAT_DQM_15.root
if [ -e RandomEngineState_15.log ]
    then
    rm RandomEngineState_15.log*
fi
if [ -e histProbFunction_15.root ]
    then
    rm histProbFunction_15.root*
fi

exit 0


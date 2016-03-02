#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc491
cd /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp  /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/BBA/RECO/BSUB/ttbarTauRecoPAT/STEP_3_cfg_ttbarTauRecoPAT_1.py .
cmsRun STEP_3_cfg_ttbarTauRecoPAT_1.py
cmsStage -f ttbarTauRecoPAT_1.root /store/user/ktos/ttbarTauRecoPAT
rm ttbarTauRecoPAT_1.root STEP_3_cfg_ttbarTauRecoPAT_1.py ttbarTauRecoPAT_DQM_1.root
if [ -e RandomEngineState_1.log ]
    then
    rm RandomEngineState_1.log*
fi
if [ -e histProbFunction_1.root ]
    then
    rm histProbFunction_1.root*
fi

exit 0


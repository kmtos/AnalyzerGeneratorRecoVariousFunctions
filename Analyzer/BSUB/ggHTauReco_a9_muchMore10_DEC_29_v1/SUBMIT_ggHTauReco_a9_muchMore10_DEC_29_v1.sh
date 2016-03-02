#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc481
cd /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/BBA/Analyzer/BSUB/ggHTauReco_a9_muchMore10_DEC_29_v1/ggHanalyzer_mshi_cfg_ggHTauReco_a9_muchMore10_DEC_29_v1.py .
cmsRun ggHanalyzer_mshi_cfg_ggHTauReco_a9_muchMore10_DEC_29_v1.py
rm ggHanalyzer_mshi_cfg_ggHTauReco_a9_muchMore10_DEC_29_v1.py 
exit 0
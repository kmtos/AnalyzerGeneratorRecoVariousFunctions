#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc481
cd /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/BBA/Analyzer/BSUB/ggH_a9_TauReco_fixed_DEC_10/ggHanalyzer_mshi_cfg_ggH_a9_TauReco_fixed_DEC_10.py .
cmsRun ggHanalyzer_mshi_cfg_ggH_a9_TauReco_fixed_DEC_10.py
rm ggHanalyzer_mshi_cfg_ggH_a9_TauReco_fixed_DEC_10.py 
exit 0

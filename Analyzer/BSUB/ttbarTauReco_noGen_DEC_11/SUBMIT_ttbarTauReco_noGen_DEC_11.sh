#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc481
cd /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/BBA/Analyzer/BSUB/ttbarTauReco_noGen_DEC_11/ttbarAnalyzer_cfg_ttbarTauReco_noGen_DEC_11.py .
cmsRun ttbarAnalyzer_cfg_ttbarTauReco_noGen_DEC_11.py
rm ttbarAnalyzer_cfg_ttbarTauReco_noGen_DEC_11.py 
exit 0

#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc481
cd /afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp /afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/WJ_New_BDiscNoPF_NoMu_AUG4/wJetsToLNUanalyzer_cfg_WJ_New_BDiscNoPF_NoMu_AUG4.py .
cmsRun wJetsToLNUanalyzer_cfg_WJ_New_BDiscNoPF_NoMu_AUG4.py
rm wJetsToLNUanalyzer_cfg_WJ_New_BDiscNoPF_NoMu_AUG4.py 
exit 0

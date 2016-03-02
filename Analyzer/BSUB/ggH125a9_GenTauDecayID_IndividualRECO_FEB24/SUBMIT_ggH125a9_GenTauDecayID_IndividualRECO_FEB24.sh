#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc481
cd /afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_4_12_patch4/src/
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp /afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_4_12_patch4/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/ggH125a9_GenTauDecayID_IndividualRECO_FEB24/ggHanalyzer_RECO_cfg_ggH125a9_GenTauDecayID_IndividualRECO_FEB24.py .
cmsRun ggHanalyzer_RECO_cfg_ggH125a9_GenTauDecayID_IndividualRECO_FEB24.py
rm ggHanalyzer_RECO_cfg_ggH125a9_GenTauDecayID_IndividualRECO_FEB24.py 
exit 0

#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc481
cd /afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_4_12_patch4/src/
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp /afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_4_12_patch4/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/DIRNAME/ANALYZER.py .
cmsRun ANALYZER.py
rm ANALYZER.py 
exit 0

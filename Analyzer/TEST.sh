#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc481
cd /afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_17/src
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp /afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/BSUB/DIRNAME/ANALYZER.py .
#cmsRun ANALYZER.py
echo "PENIS"
echo "aslkfjhasoidjfghcaiuyjDfgbcaksydghbcinkyshdgbfiknysHBDF mukayskgdfnicljakygdfnmiyawhsjndmiflahsjndmfjycasjhdbnmfljyhasjd mfjykcajsndfmunlkycajmnyfujansdf"
rm ANALYZER.py 
exit 0

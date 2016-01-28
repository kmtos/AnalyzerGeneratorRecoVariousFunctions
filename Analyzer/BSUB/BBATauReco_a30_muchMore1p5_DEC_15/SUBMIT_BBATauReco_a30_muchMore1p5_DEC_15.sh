#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc481
cd /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/BBA/Analyzer/BSUB/BBATauReco_a30_muchMore1p5_DEC_15/bbaTauRecoAnalyzer_cfg_BBATauReco_a30_muchMore1p5_DEC_15.py .
cmsRun bbaTauRecoAnalyzer_cfg_BBATauReco_a30_muchMore1p5_DEC_15.py
rm bbaTauRecoAnalyzer_cfg_BBATauReco_a30_muchMore1p5_DEC_15.py 
exit 0

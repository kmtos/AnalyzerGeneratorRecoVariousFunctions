#!/bin/bash
source  /afs/cern.ch/cms/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc481
if [ -r CMSSW_7_1_11_patch2/src ] ; then 
 echo release CMSSW_7_1_11_patch2 already exists
else
scram p CMSSW CMSSW_7_1_11_patch2
fi
cd CMSSW_7_1_11_patch2/src
eval `scram runtime -sh`

curl -s --insecure https://cms-pdmv.cern.ch/mcm/public/restapi/requests/get_fragment/HIG-RunIIFall14GS-00006 --retry 2 --create-dirs -o Configuration/GenProduction/python/HIG-RunIIFall14GS-00006-fragment.py 
[ -s Configuration/GenProduction/python/HIG-RunIIFall14GS-00006-fragment.py ] || exit $?;

scram b
pwd
cd ../../
pwd
cmsDriver.py Configuration/GenProduction/python/HIG-RunIIFall14GS-00006-fragment.py --fileout file:HIG-RunIIFall14GS-00006.root --mc --eventcontent RAWSIM --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1,Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM --conditions MCRUN2_71_V1::All --beamspot NominalCollision2015 --step GEN,SIM --magField 38T_PostLS1 --python_filename HIG-RunIIFall14GS-00006_1_cfg.py --no_exec -n 2141 || exit $? ; 



echo "nothing" ;cmsDriver.py Configuration/GenProduction/python/HIG-RunIIFall14GS-00006-fragment.py --fileout file:HIG-RunIIFall14GS-00006.root --mc --eventcontent DQM --datatier DQM --conditions MCRUN2_71_V1::All --beamspot NominalCollision2015 --step GEN,VALIDATION:genvalid_all --magField 38T_PostLS1  --fileout file:HIG-RunIIFall14GS-00006_genvalid.root --mc -n 45126 --python_filename HIG-RunIIFall14GS-00006_genvalid.py  --no_exec || exit $? ;

cmsDriver.py step2 --filein file:HIG-RunIIFall14GS-00006_genvalid.root --conditions MCRUN2_71_V1::All --mc -s HARVESTING:genHarvesting --harvesting AtJobEnd --python_filename HIG-RunIIFall14GS-00006_genvalid_harvesting.py --no_exec || exit $? ; 


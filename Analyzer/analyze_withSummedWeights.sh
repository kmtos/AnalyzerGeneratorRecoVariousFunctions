#!/bin/bash

#parse arguments
if [ $# -ne 6 ]
    then
    echo "Usage: ./generate.sh cfg_name script_name sample_tag queue name_addon lumi_data"
    exit 0
fi

cfg_name=$1
script_name=$2
tag=$3
queue=$4
name_addon=$5
input="/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/InputFor${tag}.txt"
lumi_data=$6

while IFS= read -r line
do
  linAr=($line)
  divisions=${linAr[0]}
  dir_name=${linAr[1]}
  file_path=${linAr[2]}
  full_dir_name="${linAr[1]}${name_addon}"
  echo ""
  echo ""
  mkdir -p BSUB/$full_dir_name
  cd BSUB/$full_dir_name
  xsec=`grep "$dir_name" ../../FILE_TESTS/CrossSections.py`
  xsec=${xsec##* }
  echo "xsec=$xsec"  
  summedWeights=`grep "$dir_name" /afs/cern.ch/user/k/ktos/NewSkimDir/CMSSW_8_0_30/src/GGHAA2Mu2TauAnalysis/SkimMuMuTauTau/test/SkimSequence/SummedWeightsFiles/SummedWeightsValues.out`
  summedWeights=${summedWeights##* }
  echo "summedWeights=$summedWeights"
  echo "Dir= $full_dir_name    count=$divisions"
  echo "Sub_dir_name = $dir_name"
  echo "File_path= $file_path"

  cp ../../src/FakeRate*  ../../src/DiMu* .
  divisions=$((divisions + 5))
  COUNT=1
  while [ $COUNT -le  $divisions ]; do
    echo "python file= ${cfg_name}_${full_dir_name}_${COUNT}.py"
    echo "Script name= ${script_name}_${full_dir_name}_${COUNT}.sh"


    sed -e "s|XSEC|$xsec|g" -e "s|LUMI_DATA|$lumi_data|g" -e "s|SUMMED_WEIGHTS|${summedWeights}|g" -e "s|FILE_PATH|${file_path}|g" -e "s|DIRNAME|${full_dir_name}|g" -e "s|NUM|${COUNT}|g" ../../${cfg_name}.py > ${cfg_name}_${full_dir_name}_${COUNT}.py
    sed -e "s|ANALYZER|${cfg_name}_${full_dir_name}_${COUNT}|g" -e "s|DIRNAME|${full_dir_name}|g"  ../../${script_name}.sh > ${script_name}_${full_dir_name}_${COUNT}.sh
    chmod u+x ${script_name}_${full_dir_name}_${COUNT}.sh
    bsub -q $queue -J ${cfg_name}_${full_dir_name}_${COUNT} < ${script_name}_${full_dir_name}_${COUNT}.sh
    echo "COUNT= $COUNT"
    let COUNT=COUNT+1
  done
  cd ../../
done <"$input"
exit 0

#!/bin/bash

#parse arguments
if [ $# -ne 5 ]
    then
    echo "Usage: ./generate.sh cfg_name script_name tag queue name_addon"
    exit 0
fi

cfg_name=$1
script_name=$2
tag=$3
queue=$4
name_addon=$5
input="/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS//InputFor${tag}.txt"
while IFS= read -r line
do
  linAr=($line)
  divisions=${linAr[0]}
  dir_name_TEMP=${linAr[1]}
  file_path=${linAr[2]}
  replaced='NoMassCut_NoMassCut'
  replacedWith="NoMassCut"
  dir_name="${dir_name_TEMP/$replaced/$replacedWith}$name_addon"
  echo ""
  echo ""
  echo "Dir= $dir_name    count=$divisions"
  echo "File_path= $file_path"
  mkdir -p BSUB/$dir_name
  cd BSUB/$dir_name
  path=$(pwd)
  echo "path= $path"
  cp  ../../src/FakeRate* .

  divisions=$((divisions + 5))
  COUNT=1
  while [ $COUNT -le  $divisions ]; do
    echo "DIRNAME    = ${dir_name}"
    echo "python file= ${cfg_name}_${dir_name}_${COUNT}.py"
    echo "Script name= ${script_name}_${dir_name}_${COUNT}.sh"
    sed -e "s%FILE_PATH%${file_path}%g" -e "s%DIRNAME%${dir_name}%g" -e "s%NUM%${COUNT}%g" ../../${cfg_name}.py > ${cfg_name}_${dir_name}_${COUNT}.py
    sed -e "s%ANALYZER%${cfg_name}_${dir_name}_${COUNT}%g" -e "s%DIRNAME%${dir_name}%g" -e "s%NUM%${COUNT}%g" ../../${script_name}.sh > ${script_name}_${dir_name}_${COUNT}.sh
    chmod u+x ${script_name}_${dir_name}_${COUNT}.sh
    bsub -q $queue -J ${cfg_name}_${dir_name}_${COUNT} < ${script_name}_${dir_name}_${COUNT}.sh
    echo "COUNT= $COUNT"
    let COUNT=COUNT+1
  done
  cd ../../
done <"$input"
exit 0

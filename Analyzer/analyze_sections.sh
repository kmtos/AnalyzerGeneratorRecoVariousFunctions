#!/bin/bash

#parse arguments
if [ $# -ne 5 ]
    then
    echo "Usage: ./generate.sh cfg_name script_name dir_name queue divisions"
    exit 0
fi

cfg_name=$1
script_name=$2
dir_name=$3
queue=$4
divisions=$5

mkdir -p BSUB/$dir_name
cd BSUB/$dir_name
path=$(pwd)
echo "path= $path"

COUNT=1
while [ $COUNT -le  $divisions ]; do
  echo "DIRNAME    = ${dir_name}"
  echo "python file= ${cfg_name}_${dir_name}_${COUNT}.py"
  cp ../../src/GGHAnalyzer.cc GGHAnalyzer_${dir_name}.cc
  echo "Script name= ${script_name}_${dir_name}_${COUNT}.sh"
  sed -e "s%DIRNAME%${dir_name}%g" -e "s%NUM%${COUNT}%g" ../../${cfg_name}.py > ${cfg_name}_${dir_name}_${COUNT}.py
  sed -e "s%ANALYZER%${cfg_name}_${dir_name}_${COUNT}%g" -e "s%DIRNAME%${dir_name}%g" -e "s%NUM%${COUNT}%g" ../../${script_name}.sh > ${script_name}_${dir_name}_${COUNT}.sh
  chmod u+x ${script_name}_${dir_name}_${COUNT}.sh
  bsub -q $queue -J ${cfg_name}_${dir_name}_${COUNT} < ${script_name}_${dir_name}_${COUNT}.sh
  echo "COUNT= $COUNT"
  let COUNT=COUNT+1
done
exit 0

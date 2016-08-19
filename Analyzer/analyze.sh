#!/bin/bash

#parse arguments
if [ $# -ne 4 ]
    then
    echo "Usage: ./generate.sh cfg_name script_name dir_name queue"
    exit 0
fi

cfg_name=$1
script_name=$2
dir_name=$3
queue=$4

mkdir -p BSUB/$dir_name
cd BSUB/$dir_name
path=$(pwd)
echo "path= $path"

echo "DIRNAME    = ${dir_name}"
echo "python file= ${cfg_name}_${dir_name}.py"
cp ../../src/GGHAnalyzer_IndivCJ.cc ../../src/GGHAnalyzer_IndivRECO.cc ../../src/GGHAnalyzer_OLD.cc ../../src/ZTTAnalyzer.cc ../../src/FakeRateAnalyzer.cc .
echo "ANALYZER   = $BBAAnalyzer_${dir_name}.cc"
echo "Script name= ${script_name}_${dir_name}.sh"
sed -e "s%DIRNAME%${dir_name}%g" ../../${cfg_name}.py > ${cfg_name}_${dir_name}.py
sed -e "s%ANALYZER%${cfg_name}_${dir_name}%g" -e "s%DIRNAME%${dir_name}%g" ../../${script_name}.sh > ${script_name}_${dir_name}.sh
chmod u+x ${script_name}_${dir_name}.sh
bsub -q $queue -J ${cfg_name}_${dir_name} < ${script_name}_${dir_name}.sh

exit 0

#!/bin/bash
#parse arguments
if [ $# -ne 5 ]
    then
    echo "Usage: ./generate.sh cfg_name script_name dir_name queue number_of_files"
    exit 0
fi

cfg_name=$1
script_name=$2
dir_name=$3
queue=$4
number_of_files=$5

mkdir -p BSUB/$dir_name
cd BSUB/$dir_name

#make directory on EOS
EOS_dir_query=`cmsLs /store/user/ktos/${dir_name}`
EOS_dir_query=`echo $EOS_dir_query | grep "No such file or directory"`
if [ "EOS_dir_query" != "" ]
    then
    cmsMkdir /store/user/ktos/${dir_name}
fi

#PUT FOR LOOP HERE REPLACING NUM WITH OTHER NUMBER
files="13 211 239 31 336 360 377 390 444 502 603 774 874 904 986 "

for file in ${files} 
  do
  sed -e "s%DIRNAME%${dir_name}%g" -e "s%NUM%${file}%g" ../../${cfg_name}.py >  ${cfg_name}_${dir_name}_${file}.py
  sed -e "s%GENERATOR%${cfg_name}_${dir_name}_${file}%g" -e "s%NUM%${file}%g" -e "s%DIRNAME%${dir_name}%g" ../../${script_name}.sh > ${script_name}_${dir_name}_${file}.sh
  echo "Number= $file"
  chmod u+x ${script_name}_${dir_name}_${file}.sh
  bsub -q $queue -J ${cfg_name}_${dir_name}_${file} < ${script_name}_${dir_name}_${file}.sh
done

cd ../../

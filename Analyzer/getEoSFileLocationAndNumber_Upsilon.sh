 #!/bin/bash

mv FILE_TESTS/InputForUpsilon.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/SingleMuon/*
do
  dirLs="$(echo $i | grep Upsilon)"
  if [ -z "$dirLs" ]; then
    echo "Skipped $i"
    continue
  fi
  echo "$dirLs" | while IFS= read -r line
  do
    numb=$(eos ls $line)
    fileName="$line/$numb/0000/"
    count="$(eos ls $fileName | grep selec | wc -l)"
    echo $fileName
    dir="$(echo ${line##*/})"
    echo "$count  $dir  $fileName" | grep Upsilon >> FILE_TESTS/InputForUpsilon.txt
    done
done


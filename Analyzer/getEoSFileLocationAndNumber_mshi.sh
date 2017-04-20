#!/bin/bash

mv FILE_TESTS/InputForRegionBWithMassCut.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/mshi/*
do 
  dirLs="$(eos ls $i | grep RegionBWithMassCut)"
  if [ -z "$dirLs" ]; then
    echo "Skipped $i"
    continue
  fi
  echo "$dirLs" | while IFS= read -r line
  do 
    var="$i/$line"
    numb=$(eos ls $var)
    echo "var/numb= $var/$numb"
    fileName="$var/$numb/0000/"
    echo "filename= $fileName"
    count="$(eos ls $fileName | grep selec | wc -l)"
    echo "$count  $line  $fileName" | grep RegionBWithMassCut >> FILE_TESTS/InputForRegionBWithMassCut.txt
    done  
done



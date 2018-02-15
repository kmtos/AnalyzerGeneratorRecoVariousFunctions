 #!/bin/bash

mv FILE_TESTS/InputForFakeeRatForB.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/SingleMuon/*
do
  dirLs="$(echo $i | grep _MedIsoMu2_TauDMA | grep FEB8)"
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
    echo "$count  $dir  $fileName" | grep _MedIsoMu2_TauDMA | grep FEB8 >> FILE_TESTS/InputForFakeRateForB.txt
    done
done

mv FILE_TESTS/InputForFakeeRatForC.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/SingleMuon/*
do
  dirLs="$(echo $i | grep AntiMedIsoMu2_TauDMM | grep FEB8)"
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
    echo "$count  $dir  $fileName" | grep AntiMedIsoMu2_TauDMM | grep FEB8 >> FILE_TESTS/InputForFakeRateForC.txt
    done
done


 #!/bin/bash

mv FILE_TESTS/InputForNoIsoDiTau.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
do 
  dirLs="$(eos ls $i | grep NoIsoDiTau_APR18)"
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
    count="$(eos ls $fileName | grep selec | wc -l)"
    echo "$count  $line  $fileName" | grep NoIsoDiTau >> FILE_TESTS/InputForNoIsoDiTau.txt
    done  
done


mv FILE_TESTS/InputForNoIsoDiMu.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
do
  dirLs="$(eos ls $i | grep NoIsoDiMu_APR18)"
  if [ -z "$dirLs" ]; then
    echo "Skipped $i"
    continue
  fi
  echo "$dirLs" | while IFS= read -r line
  do
    var="$i/$line"
    numb=$(eos ls $var)
    echo "$var/$numb"
    fileName="$var/$numb/0000/"
    count="$(eos ls $fileName | grep selec | wc -l)"
    echo "$count  $line  $fileName" | grep NoIsoDiMu >> FILE_TESTS/InputForNoIsoDiMu.txt
    done
done

mv FILE_TESTS/InputForSeparatedDiMu.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
do
  dirLs="$(eos ls $i | grep SeparatedDiMu_APR18)"
  if [ -z "$dirLs" ]; then
    echo "Skipped $i"
    continue
  fi
  echo "$dirLs" | while IFS= read -r line
  do
    var="$i/$line"
    numb=$(eos ls $var)
    echo "$var/$numb"
    fileName="$var/$numb/0000/"
    count="$(eos ls $fileName | grep selec | wc -l)"
    echo "$count  $line  $fileName" | grep SeparatedDiMu >> FILE_TESTS/InputForSeparatedDiMu.txt
    done
done




mv FILE_TESTS/InputForSignalRegion.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
do
  dirLs="$(eos ls $i | grep SignalRegion_APR18)"
  if [ -z "$dirLs" ]; then
    echo "Skipped $i"
    continue
  fi
  echo "$dirLs" | while IFS= read -r line
  do
    var="$i/$line"
    numb=$(eos ls $var)
    echo "$var/$numb"
    fileName="$var/$numb/0000/"
    count="$(eos ls $fileName | grep selec | wc -l)"
    echo "$count  $line  $fileName" | grep SignalRegion >> FILE_TESTS/InputForSignalRegion.txt
    done
done


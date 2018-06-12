 #!/bin/bash

mv FILE_TESTS/InputForData.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/Sin*
do
  dirLs="$(eos ls $i | grep MedIso)"
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
    echo "$count  $line  $fileName" | grep MedIso >> FILE_TESTS/InputForData.txt
    done
done

mv FILE_TESTS/InputForMu1Only.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
do
  dirLs="$(eos ls $i | grep _Mu1Only )"
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
    echo "$count  $line  $fileName" | grep Mu1Only >> FILE_TESTS/InputForMu1Only.txt
    done
done

mv FILE_TESTS/InputForTauDM.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
do
  dirLs="$(eos ls $i | grep _TauDM_ )"
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
    echo "$count  $line  $fileName" | grep TauDM_ >> FILE_TESTS/InputForTauDM.txt
    done
done


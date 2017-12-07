 #!/bin/bash
#
mv FILE_TESTS/InputForTauDMMedIso.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
do 
  dirLs="$(eos ls $i | grep TauDMMedIso)"
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
    echo "$count  $line  $fileName" | grep TauDMMedIso >> FILE_TESTS/InputForTauDMMedIso.txt
    done  
done

mv FILE_TESTS/InputForTauDMAntiMedIso.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
do
  dirLs="$(eos ls $i | grep TauDMAntiMedIso)"
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
    echo "$count  $line  $fileName" | grep TauDMAntiMedIso >> FILE_TESTS/InputForTauDMAntiMedIso.txt
    done
done

#!/bin/bash

mv FILE_TESTS/InputForMedIsoMu2IsoTauMedIso.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
do 
  dirLs="$(eos ls $i | grep _MedIsoMu2_TauDMMedIso | grep MAY1)"
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
    echo "$count  $line  $fileName" | grep _MedIsoMu2_TauDMMedIso | grep MAY1 >> FILE_TESTS/InputForMedIsoMu2IsoTauMedIso.txt
    done  
done

mv FILE_TESTS/InputForMedIsoMu2IsoAntiTauMedIso.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
do
  dirLs="$(eos ls $i | grep _MedIsoMu2_TauDMAntiMedIso | grep MAY1)"
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
    echo "$count  $line  $fileName" | grep _MedIsoMu2_TauDMAntiMedIso | grep MAY1 >> FILE_TESTS/InputForMedIsoMu2IsoAntiTauMedIso.txt
    done
done

mv FILE_TESTS/InputForAntiMedIsoMu2IsoAntiTauMedIso.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
do
  dirLs="$(eos ls $i | grep AntiMedIsoMu2_TauDMAntiMedIso | grep MAY1)"
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
    echo "$count  $line  $fileName" | grep AntiMedIsoMu2_TauDMAntiMedIso | grep MAY1 >> FILE_TESTS/InputForAntiMedIsoMu2IsoAntiTauMedIso.txt
    done
done

mv FILE_TESTS/InputForAntiMedIsoMu2IsoTauMedIso.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
do
  dirLs="$(eos ls $i | grep AntiMedIsoMu2_TauDMMedIso | grep MAY1)"
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
    echo "$count  $line  $fileName" | grep AntiMedIsoMu2_TauDMMedIso | grep MAY1 >> FILE_TESTS/InputForAntiMedIsoMu2IsoTauMedIso.txt
    done
done


mv FILE_TESTS/InputForMCZPeak.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
do
  dirLs="$(eos ls $i | grep _Mu1Only | grep MAY1)"
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
    echo "$count  $line  $fileName" | grep Mu1Only | grep MAY1 >> FILE_TESTS/InputForMCZPeak.txt
    done
done

mv FILE_TESTS/InputForTauDMOnly.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
do
  dirLs="$(eos ls $i | grep _TauDM_ | grep MAY1)"
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
    echo "$count  $line  $fileName" | grep _TauDM_ | grep MAY1 >> FILE_TESTS/InputForTauDMOnly.txt
    done
done


mv FILE_TESTS/InputForSIG.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
do
  dirLs="$(eos ls $i | grep SIG | grep MAY1)"
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
    echo "$count  $line  $fileName" | grep SIG | grep MAY1 >> FILE_TESTS/InputForSIG.txt
    done
done

mv FILE_TESTS/InputForMu1Only.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
do
  dirLs="$(eos ls $i | grep Mu1Only | grep MAY1)"
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
    echo "$count  $line  $fileName" | grep Mu1Only | grep MAY1 >> FILE_TESTS/InputForMu1Only.txt
    done
done

#mv FILE_TESTS/InputForSameSign.txt ~/Copies_Temp_etc/
#
#for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
#do 
#  dirLs="$(eos ls $i | grep SameSign | grep MAY1)"
#  if [ -z "$dirLs" ]; then
#    echo "Skipped $i"
#    continue
#  fi
#  echo "$dirLs" | while IFS= read -r line
#  do
#    var="$i/$line"
#    numb=$(eos ls $var)
#    echo "var/numb= $var/$numb"
#    fileName="$var/$numb/0000/"
#    count="$(eos ls $fileName | grep selec | wc -l)"
#    echo "$count  $line  $fileName" | grep SameSign | grep MAY1 >> FILE_TESTS/InputForSameSign.txt
#    done
#done
#

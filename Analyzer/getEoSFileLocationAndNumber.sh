 #!/bin/bash

mv FILE_TESTS/InputForNoMassCut_NoIsoDiMu.txt FILE_TESTS/InputForNoMassCut_NoIsoDiTau.txt FILE_TESTS/InputForNoMassCut_SameSignDiMu.txt FILE_TESTS/InputForNoMassCut_SeparatedDiMu.txt ~/Copies_Temp_etc/
for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
do 
  dirLs="$(eos ls $i | grep NoMassCut)"
  echo "$dirLs" | while IFS= read -r line
  do 
    var="$i/$line"
    numb=$(eos ls $var)
    fileName="$var/$numb/0000/"
    count="$(eos ls $fileName | grep selec | wc -l)"
    echo "$count  $line  $fileName" | grep NoIsoDiMu >> FILE_TESTS/InputForNoMassCut_NoIsoDiMu.txt
    echo "$count  $line  $fileName" | grep NoIsoDiTau >> FILE_TESTS/InputForNoMassCut_NoIsoDiTau.txt
    echo "$count  $line  $fileName" | grep SameSignDiMu >> FILE_TESTS/InputForNoMassCut_SameSignDiMu.txt
    echo "$count  $line  $fileName" | grep SeparatedDiMu >> FILE_TESTS/InputForNoMassCut_SeparatedDiMu.txt
    done  
done

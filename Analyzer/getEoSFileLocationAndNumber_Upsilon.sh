 #!/bin/bash

mv FILE_TESTS/InputForUpsilon.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/SingleMuon/*
do
  dirLs="$(echo $i | grep Upsilon)"
  echo "i= $i"
  echo "dirLs= $dirLs"
  echo "$dirLs" | while IFS= read -r line
  do
    numb=$(eos ls $i)
    echo "numb= $numb"
    fileName="$i/$numb/0000/"
    count="$(eos ls $fileName | grep selec | wc -l)"
    echo "$count  $line  $fileName" | grep Upsilon >> FILE_TESTS/InputForUpsilon.txt
    done
done


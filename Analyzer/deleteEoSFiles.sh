 #!/bin/bash

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
do 
  dirLs="$(eos ls $i | grep MAR)"
  if [ -z "$dirLs" ]; then
    echo "Skipped $i"
    continue
  fi
  echo "$dirLs" | while IFS= read -r line
  do 
    echo "Removing $i/$line"
    eos rm -r "$i/$line"
    done  
done


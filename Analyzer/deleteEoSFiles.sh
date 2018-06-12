 #!/bin/bash

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/SingleMuon/SameSign_MiniAOD_Run2016F-23Sep2016-v1/180323_172451/0000/*; 
do 
  temp=${i%.root*}
  temp2=${temp##*/} 
  num=${temp2##*_}
  if [[ $num -ge 1000 ]]
  then
    eos rm -r $i
  fi
done
#for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
#do 
#  dirLs="$(eos ls $i | grep MAR)"
#  if [ -z "$dirLs" ]; then
#    echo "Skipped $i"
#    continue
#  fi
#  echo "$dirLs" | while IFS= read -r line
#  do 
#    echo "Removing $i/$line"
#    eos rm -r "$i/$line"
#    done  
#done


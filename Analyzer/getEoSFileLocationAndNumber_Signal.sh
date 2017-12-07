 #!/bin/bash
#
mv FILE_TESTS/InputForSignalMC.txt ~/Copies_Temp_etc/

for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/SUSY*
do 
  dirLs="$(eos ls $i | grep TauDMMedIso_I)"
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
    echo "$count  $line  $fileName" | grep "TauDMMedIso_I" >> FILE_TESTS/InputForSignalMC.txt
    done  
done

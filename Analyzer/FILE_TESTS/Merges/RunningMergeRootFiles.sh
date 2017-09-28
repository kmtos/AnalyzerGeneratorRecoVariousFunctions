# .bashrc

cat /afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/LumiTree_Directories.out | while read line; 
do
  echo $line
  runNum="$(eos ls $line)"
  echo "$line/$runNum"
  sed -i "s|EOS_PATH_|${line}|g" MergeRootFiles_1000.C
  sed -i "s|RUN_NUM|${runNum}|g" MergeRootFiles_1000.C
  sed -i "s|OUTFILENAME|${line##*/}|g" MergeRootFiles_1000.C
  sed -i "s|MergeRootFiles_1000()|MergeRootFiles_1000_${line##*/}()|g" MergeRootFiles_1000.C
  #cat MergeRootFiles_1000.C
  cp MergeRootFiles_1000.C MergeRootFiles_1000_${line##*/}.C
  #root -l  MergeRootFiles_1000.C  
  sed -i "s|${line}|EOS_PATH_|g" MergeRootFiles_1000.C
  sed -i "s|${runNum}|RUN_NUM|g" MergeRootFiles_1000.C
  sed -i "s|MergeRootFiles_1000_${line##*/}()|MergeRootFiles_1000()|g" MergeRootFiles_1000.C
  sed -i "s|${line##*/}|OUTFILENAME|g" MergeRootFiles_1000.C
done


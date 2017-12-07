# .bashrc

cd BSUB/
for i in ./*;
do
  sed -i "s|TEMPLATE_DIRNAME|${i##*./}|g" ../MergeRootFiles_1000.C
  root -l .x ../MergeRootFiles_1000.C
  sed -i "s|${i##*./}|TEMPLATE_DIRNAME|g" ../MergeRootFiles_1000.C
done



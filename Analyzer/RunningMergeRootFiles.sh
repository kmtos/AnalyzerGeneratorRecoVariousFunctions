# .bashrc

cd BSUB/
for i in M*;
do
  cd $i
  hadd ../FINAL_${i}.root *.root
  cd ..
done



# .bashrc

cd BSUB/
for i in ./*; do 
  echo $i
  grep -lr "Exc" $i  | wc -l 
done

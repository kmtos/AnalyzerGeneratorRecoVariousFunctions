#!/bin/bash

#for i in BSUB/M*Single*; 
#do
#  sed -i "s|PROCESSNAME|${i##*BSUB/}|g" rootMacro_CombineShapePlots_Big.C
#  root -l rootMacro_CombineShapePlots_Big.C 
#  sed -i "s|${i##*BSUB/}|PROCESSNAME|g" rootMacro_CombineShapePlots_Big.C
#done
#
for i in BSUB/FINAL_RooDataSet_MiniAOD_SingleMu1*; 
do
  sed -i "s|SUBSCRIPT|${i##*SingleMu1_}|g" rootMacro_CombineShapePlots_Data.C
  root -l rootMacro_CombineShapePlots_Data.C 
  sed -i "s|${i##*SingleMu1_}|SUBSCRIPT|g" rootMacro_CombineShapePlots_Data.C
done


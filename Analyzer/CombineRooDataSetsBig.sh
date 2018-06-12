#!/bin/bash

for i in BSUB/M*Single*; 
do
  sed -i "s|PROCESSNAME|${i##*BSUB/}|g" rootMacro_CombineShapePlots_Big.C
  root -l rootMacro_CombineShapePlots_Big.C 
  sed -i "s|${i##*BSUB/}|PROCESSNAME|g" rootMacro_CombineShapePlots_Big.C
done

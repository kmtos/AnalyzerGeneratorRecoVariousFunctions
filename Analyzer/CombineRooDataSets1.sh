#!/bin/bash

for i in BSUB/M*300*;
do
  sed -i "s|PROCESSNAME|${i##*BSUB/}|g" rootMacro_CombineShapePlots1.C
  root -l rootMacro_CombineShapePlots1.C 
  sed -i "s|${i##*BSUB/}|PROCESSNAME|g" rootMacro_CombineShapePlots1.C
done

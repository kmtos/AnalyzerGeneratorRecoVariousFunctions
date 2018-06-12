#!/bin/bash

for i in BSUB/M*125*;
do
  sed -i "s|PROCESSNAME|${i##*BSUB/}|g" rootMacro_CombineShapePlots.C
  root -l rootMacro_CombineShapePlots.C 
  sed -i "s|${i##*BSUB/}|PROCESSNAME|g" rootMacro_CombineShapePlots.C
done

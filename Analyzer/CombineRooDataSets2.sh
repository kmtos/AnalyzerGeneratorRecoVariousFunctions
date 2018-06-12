#!/bin/bash

for i in BSUB/M*750*;
do
  sed -i "s|PROCESSNAME|${i##*BSUB/}|g" rootMacro_CombineShapePlots2.C
  root -l rootMacro_CombineShapePlots2.C 
  sed -i "s|${i##*BSUB/}|PROCESSNAME|g" rootMacro_CombineShapePlots2.C
done

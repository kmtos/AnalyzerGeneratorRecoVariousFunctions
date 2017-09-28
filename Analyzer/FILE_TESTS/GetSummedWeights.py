import ROOT

TFile = ROOT.TFile.Open("/afs/cern.ch/user/k/ktos/GroupDir/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/OUTFILENAME_GenWeights.root")

Tree = TFile.Get("LumiTree")
summedWeights = 0
for row in Tree:
    summedWeights += row.summedWeights

print '{0:f}'.format(summedWeights)


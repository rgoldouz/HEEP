import math
import gc
import sys
import ROOT
import numpy as np
import copy
import os
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1;")
ROOT.TH1.AddDirectory(ROOT.kFALSE)
ROOT.gStyle.SetOptStat(0)
from array import array
from ROOT import TColor
from ROOT import TGaxis
from ROOT import THStack
import gc
TGaxis.SetMaxDigits(2)
import math

def cutFlowTable(hists, samples, regions, ch, year,caption='2016', nsig=4):
    mcSum = list(0 for i in xrange(0,len(regions)))
    for ids, s in enumerate(samples):
        if ids==0:
            continue
        for idr, r in enumerate(regions):
            if ids<nsig:
                mcSum[idr] += hists[year][ids][ch][idr][0].Integral() 
#    table = '\\begin{sidewaystable*}' + "\n"
    table = '\\begin{table*}' + "\n"
    table += '\\centering' + "\n"
    table += '\\caption{Number of expected events from \ttbar, tW, and from the remaining backgrounds (other), total background contribution and observed events in data, collected during three years (2016, 2017, and 2018), after all selections  for signal (1 b-tagged) and control ($>$ 1 b-tagged) regions. The expected signal yields for single top quark production and top decays via the scalar, vector and tensor CLFV interactions, assuming $C_x/\Lambda^2 = 1$ TeV$^{-2}$ are also shown.  The uncertainties correspond to the statistical contribution only.  }' + "\n"
#    table += '\\resizebox{\\textwidth}{!}{ \n'
    table += '\\begin{tabular}{lllllll}' + "\n"
    table += '\\hline' + "\n"
    table += 'Channel & ' + ' & '.join(regions) + '\\\\' + "\n"
    table += '\\hline' + "\n"
    for ids, s in enumerate(samples):
        if ids==0:
            continue
        if ids>=nsig:
            continue
        table += s 
        for idr, r in enumerate(regions):
            if ids<nsig:
                table += ' & ' + str(int(hists[year][ids][ch][idr][0].Integral())) + '$\pm$' + str(int(math.sqrt(hists[year][ids][ch][idr][0].GetSumw2().GetSum() )))
#            else:
#                table += (' & ' + str(round(hists[year][ids][ch][idr][0].Integral(),2))) 
        table += '\\\\' + "\n"    
    table += '\\hline' + "\n"
    table += 'Prediction '
    for idr, r in enumerate(mcSum):
        table += (' & ' + str(int(r)))
    table += '\\\\' + "\n"
    table += 'Data '
    for idr, r in enumerate(regions):
        table += (' & ' + str(int(hists[year][0][ch][idr][0].Integral())))
    table += '\\\\' + "\n"
#    table += '\\hline' + "\n"
#    table += 'Data$/$Pred. '
#    for idr, r in enumerate(mcSum):
#        table += (' & ' + str(round(hists[year][0][ch][idr][0].Integral()/r,2)))
#    table += '\\\\' + "\n"
    table += '\\hline' + "\n"
    for ids, s in enumerate(samples):
        if ids<nsig:
            continue
        table += s
        for idr, r in enumerate(regions):
            table += ' & ' + str(round(hists[year][ids][ch][idr][0].Integral(),1))+ '$\pm$' + str(round(math.sqrt(hists[year][ids][ch][idr][0].GetSumw2().GetSum() ),1))
        table += '\\\\' + "\n"
    table += '\\hline' + "\n"
    table += '\\end{tabular}' + "\n"
    table += '\\end{table*}' + "\n"
#    table += '\\end{sidewaystable*}' + "\n"
    print table

#year=['2016','2017','2018','All']
year=['All']
#regions=["ll","llOffZ","llB1", "llBg1"]
regions=["llB1", "llBg1"]
regionsName=["1 b-tagged", "$>$ 1 b-tagged"]
#channels=["ee","emu","mumu"]
channels=["emu"]
#variables=["lep1Pt","lep1Eta","lep1Phi","lep2Pt","lep2Eta","lep2Phi","llM","llPt","llDr","llDphi","jet1Pt","jet1Eta","jet1Phi","njet","nbjet","Met","MetPhi","nVtx","llMZw"]
variables=["lep1Phi"]
variablesName=["p_{T}(leading lepton)","#eta(leading lepton)","#Phi(leading lepton)","p_{T}(sub-leading lepton)","#eta(sub-leading lepton)","#Phi(sub-leading lepton)","M(ll)","p_{T}(ll)","#Delta R(ll)","#Delta #Phi(ll)","p_{T}(leading jet)","#eta(leading jet)","#Phi(leading jet)","Number of jets","Number of b-tagged jets","MET","#Phi(MET)","Number of vertices", "M(ll) [z window]"]



HistAddress = '/user/rgoldouz/NewAnalysis2020/Analysis/hists/'

Samples = ['data.root','TTTo2L2Nu.root', 'ST_tW.root','DYwjetOther.root', 'LFVTtVecU.root', 'LFVStVecU.root', 'LFVTtVecC.root','LFVStVecC.root','LFVTtTensorU.root','LFVStTensorU.root','LFVTtTensorC.root','LFVStTensorC.root','LFVTtScalarU.root','LFVStScalarU.root','LFVTtScalarC.root','LFVStScalarC.root']
SamplesName = ['Data','t#bar{t}', 'tW','Others',
 'CLFV top decay - vector - $e\mu tu$','CLFV single top - vector - $e\mu tu$',
 'CLFV top decay - vector - $e\mu tc$','CLFV single top - vector - $e\mu tc$',
 'CLFV top decay - tensor - $e\mu tu$','CLFV single top - tensor - $e\mu tu$',
 'CLFV top decay - tensor - $e\mu tc$','CLFV single top - tensor - $e\mu tc$',
 'CLFV top decay - scalar - $e\mu tu$','CLFV single top - scalar - $e\mu tu$',
 'CLFV top decay - scalar - $e\mu tc$','CLFV single top - scalar - $e\mu tc$']
SamplesNameLatex = ['Data','t#bar{t}', 'tW','Others',
 'CLFV top decay - vector - $e\mu tu$','CLFV single top - vector - $e\mu tu$',
 'CLFV top decay - vector - $e\mu tc$','CLFV single top - vector - $e\mu tc$',
 'CLFV top decay - tensor - $e\mu tu$','CLFV single top - tensor - $e\mu tu$',
 'CLFV top decay - tensor - $e\mu tc$','CLFV single top - tensor - $e\mu tc$',
 'CLFV top decay - scalar - $e\mu tu$','CLFV single top - scalar - $e\mu tu$',
 'CLFV top decay - scalar - $e\mu tc$','CLFV single top - scalar - $e\mu tc$']

colors =  [ROOT.kBlack,ROOT.kYellow,ROOT.kGreen,ROOT.kBlue-3,ROOT.kRed-4,ROOT.kOrange-3, ROOT.kOrange-6, ROOT.kCyan-6,ROOT.kGreen,ROOT.kBlue-3,ROOT.kRed-4,ROOT.kOrange-3, ROOT.kOrange-6, ROOT.kCyan-6,ROOT.kBlue-3,ROOT.kRed-4,ROOT.kOrange-3, ROOT.kOrange-6, ROOT.kCyan-6,ROOT.kGreen,ROOT.kBlue-3,ROOT.kRed-4,ROOT.kOrange-3, ROOT.kOrange-6, ROOT.kCyan-6]

Hists = []
for numyear, nameyear in enumerate(year):
    l0=[]
    Files = []
    for f in range(len(Samples)):
        l1=[]
        Files.append(ROOT.TFile.Open(HistAddress + nameyear+ '_' + Samples[f]))
        for numch, namech in enumerate(channels):
            l2=[]
            for numreg, namereg in enumerate(regions):
                l3=[]
                for numvar, namevar in enumerate(variables):
                    h= Files[f].Get(namech + '_' + namereg + '_' + namevar)
                    h.SetFillColor(colors[f])
                    h.SetLineColor(colors[f])
                    l3.append(h)
                l2.append(l3)
            l1.append(l2)
        l0.append(l1)
    Hists.append(l0)       

for numyear, nameyear in enumerate(year):
    for numch, namech in enumerate(channels):
        for numreg, namereg in enumerate(regions):
            for numvar, namevar in enumerate(variables):
                HH=[]
                HHsignal=[]
                for f in range(len(Samples)):
                    if 'LFV' in Samples[f]:
                        HHsignal.append(Hists[numyear][f][numch][numreg][numvar])
                    else:
                        HH.append(Hists[numyear][f][numch][numreg][numvar])

#                stackPlots(HH, HHsignal, SamplesName, namech, namereg, nameyear,namevar,variablesName[numvar])

le = '\\documentclass{article}' + "\n"
le += '\\usepackage{rotating}' + "\n"
le += '\\usepackage{rotating}' + "\n"
le += '\\begin{document}' + "\n"

print le
for numyear, nameyear in enumerate(year):
    for numch, namech in enumerate(channels):
        cutFlowTable(Hists, SamplesNameLatex, regionsName, numch, numyear, nameyear + ' ' + namech, 4 )
print '\\end{document}' + "\n"



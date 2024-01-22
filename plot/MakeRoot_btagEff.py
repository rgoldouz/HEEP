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
#TGaxis.SetMaxDigits(2)

def draw2dHist(A, name):
    can_name = 'can'
    canvas = ROOT.TCanvas(can_name,can_name,10,10,1100,628)
    canvas.cd()
    A.SetTitle(name)
    A.GetXaxis().SetTitleSize(0.05)
    A.GetYaxis().SetTitleSize(0.05)
    A.GetYaxis().SetTitleOffset(0.7)
#    A.GetXaxis().SetTitle("eta");
#    A.GetYaxis().SetTitle("pt")
#    A.GetZaxis().SetRangeUser(0, 6)
    A.Draw("COLZ")
    canvas.Print(name + ".png")
    del canvas
    gc.collect()

num=['h2_BTaggingEff_Num_b', 'h2_BTaggingEff_Num_c', 'h2_BTaggingEff_Num_udsg']
denom=['h2_BTaggingEff_Denom_b', 'h2_BTaggingEff_Denom_c', 'h2_BTaggingEff_Denom_udsg']
Name = ['h2_BTaggingEff_b', 'h2_BTaggingEff_c', 'h2_BTaggingEff_udsg']
HistAddress = '/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/hists/'

Samples = ['mc_2016preVFP.root','mc_2016postVFP.root','mc_2017.root','mc_2018.root']
year = ['2016preVFP', '2016postVFP', '2017', '2018']
Hists = []

for numyear, nameyear in enumerate(Samples):
    f = ROOT.TFile.Open(HistAddress + nameyear)
    for numvar, namevar in enumerate(denom):
        h1= f.Get(num[numvar])
        h2= f.Get(denom[numvar])
        h1.Divide(h2)
        h1.SetName(year[numyear] + '_' + Name[numvar])
        Hists.append(h1)
        draw2dHist(h1, 'DeepJet_MediumWP_' + year[numyear] + '_' + Name[numvar])
    f.Close()

myfile = ROOT.TFile( 'btagEff.root', 'RECREATE')
for H in Hists:
    H.Write()

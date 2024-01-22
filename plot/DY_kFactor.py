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
from ROOT import TH2F
import gc
#TGaxis.SetMaxDigits(2)

leptonPTbins = array( 'd',[0,25,50,75,100,125,150,175, 200, 250, 300, 350, 400, 600, 800, 1100, 3000] )

HistAddress = '/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/hists/2017_DY.root'
A=[]
F = ROOT.TFile.Open(HistAddress)
hDY50 = F.Get("massDY50")
hDY50Bins = F.Get("massDY50Bins")
hDY50 = hDY50.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
hDY50Bins = hDY50Bins.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
hDY50Bins.Divide(hDY50)
canvas = ROOT.TCanvas('R+I','R+I',10,10,1100,628)
canvas.cd()
hDY50Bins.Draw()
canvas.Print("DYkFactormass.png")
del canvas
lPtDY50 = F.Get("lPtDY50")
lPtDY50Bins = F.Get("lPtDY50Bins")
lPtDY50 = lPtDY50.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
lPtDY50Bins = lPtDY50Bins.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
lPtDY50Bins.Divide(lPtDY50)
canvas = ROOT.TCanvas('R+I','R+I',10,10,1100,628)
canvas.cd()
lPtDY50Bins.Draw()
canvas.Print("DYkFactorlPt.png")

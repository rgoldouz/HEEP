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


HistAddress = '/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/hists/'
dir_list = os.listdir(HistAddress)
for key in dir_list:
    if 'BNV' not in key:
        continue
    print key
    FR = ROOT.TFile.Open(HistAddress +key)
    HEFT = FR.Get("crossSection")
    HEFT.GetSumFit().save('Coup/' + key.split(".")[0]+'tex')

mt=172.5
pi=3.14
s=1
t=1
Lam = 1000
A = 4* s**2
B = 4 * t**2
C= 2*s*t
topWidth=1.33
ttLOxs = 2*445
gammaBNV = (mt**5 * (A+B+C))/(192*16* pi**3 * Lam**4 * topWidth)
print str(gammaBNV)
#gammaBNV = (mt**5 * 505)/(192*16* pi**3 * Lam**4 * topWidth)
#gammaBNV = (mt**5)/(192*16* pi**3 * Lam**4)
#print str(4*gammaBNV) +' Ct^2 + ' + str(2*gammaBNV) +' CtCs + ' + str(4*gammaBNV) +' Cs^2'

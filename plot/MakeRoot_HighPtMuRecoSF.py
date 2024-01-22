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

Name=['h2_HighPtMuRecoSF']
#HistAddress = '/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/'
HistAddress = '/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/plot/'
Samples = ['2016preVFP.root','2016postVFP.root','2017.root','2018.root']
year = {
'2016preVFP' :[[[1,0],[0.9936,0.0009],[0.993,0.001],[0.993,0.002],[0.990,0.004],[0.990,0.003],[0.989,0.004],[0.8,0.3]], 
              [[1,0],[0.9930,0.0010],[0.991,0.001],[0.985,0.001],[0.981,0.002],[0.979,0.004],[0.978,0.005],[0.9,0.2]]],
'2016postVFP':[[[1,0],[0.9936,0.0009],[0.993,0.001],[0.993,0.002],[0.990,0.004],[0.990,0.003],[0.989,0.004],[0.8,0.3]],
              [[1,0],[0.9930,0.0010],[0.991,0.001],[0.985,0.001],[0.981,0.002],[0.979,0.004],[0.978,0.005],[0.9,0.2]]],
'2017'       :[[[1,0],[0.9950,0.0007],[0.996,0.001],[0.996,0.001],[0.994,0.001],[1.003,0.006],[0.987,0.003],[0.9,0.1]],
              [[1,0],[0.9930,0.0010],[0.989,0.001],[0.986,0.001],[0.989,0.001],[0.983,0.003],[0.986,0.006],[1.01,0.01]]],     
'2018'       :[[[1,0],[0.9948,0.0007],[0.995,0.009],[0.994,0.001],[0.9914,0.009],[0.993,0.002],[0.991,0.004],[1,0.1]],
              [[1,0],[0.9930,0.0010],[0.990,0.001],[0.988,0.001],[0.981,0.002],[0.983,0.003],[0.978,0.006],[0.98,0.03]]]
}
Hists = []
binsP = array( 'd',[50,100,150,200,300,400,600,1500,3500] )
binsEta = array( 'd',[0,1.6,2.4] )

for key, value in year.items():
    myfile = ROOT.TFile( HistAddress + 'HighPtMuRecoSF_'+ key +'.root', 'RECREATE')
    hpxpy  = TH2F( 'h2_HighPtMuRecoSF_pVsAbsEta', 'h2_HighPtMuRecoSF_pVsAbsEta', len(binsP)-1, binsP, len(binsEta)-1, binsEta)
    for p in range(1,len(binsP)):
        for e in range(1,len(binsEta)):
            hpxpy.SetBinContent(p,e,value[e-1][p-1][0])
            hpxpy.SetBinError(p,e,value[e-1][p-1][1])
    hpxpy.Write()
    myfile.Close()


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

bins = array( 'd',[-1,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1] )
BDTmin=-1
BDTmax=1

def compare2Hist(A, B, textA="A", textB="B", label_name="sample", can_name="can", axis_name="eta"):
    a=b=d=''
    c=0
    cc="All"
    canvas = ROOT.TCanvas(can_name,can_name,50,50,865,780)
    canvas.cd()

    pad1=ROOT.TPad("pad1", "pad1", 0, 0.315, 1, 0.99 , 0)#used for the hist plot
    pad2=ROOT.TPad("pad2", "pad2", 0, 0.0, 1, 0.305 , 0)#used for the ratio plot
    pad1.Draw()
    pad2.Draw()
    pad2.SetGridy()
    pad2.SetTickx()
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.14)
    pad1.SetRightMargin(0.05)
    pad2.SetTopMargin(0.1)
    pad2.SetBottomMargin(0.4)
    pad2.SetLeftMargin(0.14)
    pad2.SetRightMargin(0.05)
    pad2.SetFillStyle(0)
    pad1.SetFillStyle(0)
    pad1.cd()
    pad1.SetLogx(ROOT.kFALSE)
    pad2.SetLogx(ROOT.kFALSE)
    pad1.SetLogy(ROOT.kFALSE)

    A.SetLineColor( 2 )
    B.SetLineColor( 4 )
    A.SetTitle("")
    A.GetXaxis().SetTitle(axis_name)
    A.GetXaxis().CenterTitle()
    if c==0:
        A.GetYaxis().SetTitle('TP rate')
    else:
        A.GetYaxis().SetTitle('Stub rate '+d)
    A.GetXaxis().SetTitleSize(0.05)
    A.GetYaxis().SetTitleSize(0.05)
    A.GetXaxis().SetLabelSize(0)
    if b == "nstub":
        A.GetYaxis().SetTitle('Number of module')
    if b == "type":
        A.GetXaxis().SetBinLabel(1,"All")
        A.GetXaxis().SetBinLabel(2,"Genuine")
        A.GetXaxis().SetBinLabel(3,"Combinatoric")
        A.GetXaxis().SetBinLabel(4,"Unknown")
    A.SetMaximum(1.4*max(A.GetMaximum(),B.GetMaximum()));
    A.SetMinimum(0.8*min(A.GetMinimum(),B.GetMinimum()));
    A.Draw()
    B.Draw('esame')
    A.Draw("AXISSAMEY+")
    A.Draw("AXISSAMEX+")

    legend = ROOT.TLegend(0.67,0.67,0.9,0.85)
    legend.AddEntry(A ,textA,'l')
    legend.AddEntry(B ,textB,'l')
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    legend.Draw("same")
    label = ROOT.TLatex()
    label.SetTextAlign(12)
    label.SetTextFont(42)
    label.SetTextSize(0.06)
    label.SetNDC(ROOT.kTRUE)
    label.DrawLatex(0.25,0.95,"CMS Phase 2 Simulation Preliminary")
    if a == 'Barrel':
        label.DrawLatex(0.2,0.85, a + ", Layer " + str(cc))
    if a == 'Endcap':
        label.DrawLatex(0.2,0.85, a + ", Disk " + str(cc))

    pad2.cd()
    ratio = A.Clone()
    ratio.Divide(B)
    ratio.SetTitle("")
    ratio.SetMaximum(1)
    ratio.SetMinimum(0)
    ratio.GetXaxis().SetTitle(axis_name)
    ratio.GetYaxis().CenterTitle()
    ratio.GetXaxis().SetMoreLogLabels()
    ratio.GetXaxis().SetNoExponent()
    ratio.GetXaxis().SetTitleSize(0.04/0.3)
    ratio.GetYaxis().SetTitleSize(0.04/0.3)
    ratio.GetXaxis().SetTitleFont(42)
    ratio.GetYaxis().SetTitleFont(42)
    ratio.GetXaxis().SetTickLength(0.05)
    ratio.GetYaxis().SetTickLength(0.05)
    ratio.GetXaxis().SetLabelSize(0.115)
    ratio.GetYaxis().SetLabelSize(0.089)
    ratio.GetXaxis().SetLabelOffset(0.02)
    ratio.GetYaxis().SetLabelOffset(0.01)
    ratio.GetYaxis().SetTitleOffset(0.42)
    ratio.GetXaxis().SetTitleOffset(1.1)
    ratio.GetYaxis().SetNdivisions(504)
    ratio.SetStats(ROOT.kFALSE)
    ratio.GetYaxis().SetTitle('Ratio')

    ratio.Draw("e")
    ratio.Draw("AXISSAMEY+")
    ratio.Draw("AXISSAMEX+")
    print "2H_" + can_name +"_"+ label_name + " = " + str(A.Integral()) + ' , ' +str(B.Integral())
    if A.Integral()>0:
        canvas.Print("2H_" + can_name +"_"+ label_name + ".png")
    del canvas
    gc.collect()



#year=['2016','2017','2018','All']
year=['2017']
#year=['2016preVFP', '2016postVFP', '2017','2018']
regions=["ll","llOffZ","llB1"]
#regions=["llB1"]
channels=["ll"];
variables=["lep1Phi"]
#variables=["lep1Pt"]

HistAddress = '/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/hists/'


Samples = [
#'ttbar.root',
#'UL18_TTTo2L2Nu.root', 
'STBNV_TBCE.root',
'STBNV_TBUE.root',
'STBNV_TDCE.root',
'STBNV_TDUE.root',
'STBNV_TSCE.root',
'STBNV_TSUE.root',
'TTBNV_TBCE.root',
'TTBNV_TBUE.root',
'TTBNV_TDCE.root',
'TTBNV_TDUE.root',
'TTBNV_TSCE.root',
'TTBNV_TSUE.root',
]

SamplesMu = [
#'ttbar.root',
#'UL18_TTTo2L2Nu.root',
'STBNV_TBCMu.root',
'STBNV_TBUMu.root',
'STBNV_TDCMu.root',
'STBNV_TDUMu.root',
'STBNV_TSCMu.root',
'STBNV_TSUMu.root',
'TTBNV_TBCMu.root',
'TTBNV_TBUMu.root',
'TTBNV_TDCMu.root',
'TTBNV_TDUMu.root',
'TTBNV_TSCMu.root',
'TTBNV_TSUMu.root',
]
#'LFVVecU.root', 'LFVScalarU.root','LFVTensorU.root', 'LFVVecC.root', 'LFVScalarC.root','LFVTensorC.root']
#SamplesName = ['t#bar{t}', 'LFV-Vector [e#mutu]',  'LFV-Scalar [e#mutu]',  'LFV-Tensor [e#mutu]', 'LFV-Vector [e#mutc]',  'LFV-Scalar [e#mutc]',  'LFV-Tensor [e#mutc]']

#colors =  [ROOT.kRed-4,ROOT.kOrange-6, ROOT.kCyan-6,ROOT.kOrange-6,ROOT.kCyan-6,ROOT.kOrange-6, ROOT.kCyan-6]
#Style =[1,1,1,7,7,3,3]
colors =  [ROOT.kRed-4,ROOT.kBlue, ROOT.kBlue-7, ROOT.kViolet, ROOT.kViolet-5,ROOT.kAzure-9, ROOT.kAzure+10,
ROOT.kGreen, ROOT.kGreen+2, ROOT.kTeal-9, ROOT.kTeal+10, ROOT.kSpring-7, ROOT.kSpring+9]

Style =[1,1,1,1,1,1,1,2,2,2,2,2,2,2,2]
wc1 = ROOT.WCPoint("EFTrwgt1_cS_0_cT_1")

effPlots=[]
testFile = ROOT.TFile.Open("/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/ANoutput.root") 
compare2Hist(testFile.Get("eleEffNum"), testFile.Get("eleEffDen"))

CS = []
Hists = []
for numyear, nameyear in enumerate(year):
    l0=[]
    CS0=[]
    Files = []
    for f in range(len(Samples)):
        l1=[]
        Files.append(ROOT.TFile.Open(HistAddress + nameyear+ '_' + Samples[f]))
        h = Files[f].Get("crossSection")
        h.Scale(wc1)
        CS0.append(h)
        for numch, namech in enumerate(channels):
            l2=[]
            for numreg, namereg in enumerate(regions):
                l3=[]
                for numvar, namevar in enumerate(variables):
                    h= Files[f].Get(namech + '_' + namereg + '_' + namevar)
                    h.SetFillColor(colors[f])
                    h.SetLineColor(colors[f])
#                    h.SetLineStyle(Style[f])
                    h.SetLineWidth(2)
                    if 'BNV' in Samples[f]:
                        h.Scale(wc1)
                    l3.append(h)
                l2.append(l3)
            l1.append(l2)
        l0.append(l1)
    Hists.append(l0)    
    CS.append(CS0)   

CSMu = []
HistsMu = []
for numyear, nameyear in enumerate(year):
    l0=[]
    CS0=[]
    Files = []
    for f in range(len(SamplesMu)):
        l1=[]
        Files.append(ROOT.TFile.Open(HistAddress + nameyear+ '_' + SamplesMu[f]))
        h = Files[f].Get("crossSection")
        h.Scale(wc1)
        CS0.append(h)
        for numch, namech in enumerate(channels):
            l2=[]
            for numreg, namereg in enumerate(regions):
                l3=[]
                for numvar, namevar in enumerate(variables):
                    h= Files[f].Get(namech + '_' + namereg + '_' + namevar)
                    h.SetFillColor(colors[f])
                    h.SetLineColor(colors[f])
#                    h.SetLineStyle(Style[f])
                    h.SetLineWidth(2)
                    if 'BNV' in SamplesMu[f]:
                        h.Scale(wc1)
                    l3.append(h)
                l2.append(l3)
            l1.append(l2)
        l0.append(l1)
    HistsMu.append(l0)
    CSMu.append(CS0)

for numyear, nameyear in enumerate(year):
    for numch, namech in enumerate(channels):
        for numreg, namereg in enumerate(regions):
            for numvar, namevar in enumerate(variables):
                for f in range(len(Samples)):
#                    print Samples[f]+'_'+nameyear +'_'+namech +'_' + namereg
#                    print str(Hists[numyear][f][numch][numreg][numvar].Integral())
#                    print str(41530*CS[numyear][f].Integral())
#                    print str((100* Hists[numyear][f][numch][numreg][numvar].Integral())/(41530*CS[numyear][f].Integral()))
#                    print SamplesMu[f]+'_'+nameyear +'_'+namech +'_' + namereg
#                    print str(HistsMu[numyear][f][numch][numreg][numvar].Integral())
#                    print str(41530*CSMu[numyear][f].Integral())
#                    print str((100* HistsMu[numyear][f][numch][numreg][numvar].Integral())/(41530*CSMu[numyear][f].Integral()))
                    print "Combined eff:" + str((100* (Hists[numyear][f][numch][numreg][numvar].Integral() + HistsMu[numyear][f][numch][numreg][numvar].Integral()))/(41530*(CSMu[numyear][f].Integral()+CS[numyear][f].Integral())))

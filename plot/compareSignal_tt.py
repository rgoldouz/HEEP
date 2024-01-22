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

def EFTtoNormal(H, wc):
    hpx    = ROOT.TH1F( H.GetName(), H.GetName(), H.GetXaxis().GetNbins(), H.GetXaxis().GetXmin(),H.GetXaxis().GetXmax() )
    r=1
    for b in range(hpx.GetNbinsX()):
        if H.GetBinContent(b+1,ROOT.WCPoint("NONE"))>0:
            r = H.GetBinError(b+1)/H.GetBinContent(b+1,ROOT.WCPoint("NONE"))
        hpx.SetBinContent(b+1, H.GetBinContent(b+1,wc))
        hpx.SetBinError(b+1, r*H.GetBinContent(b+1,wc))
        hpx.SetLineColor(H.GetLineColor())
        hpx.SetLineStyle(H.GetLineStyle())
    if hpx.Integral()>0:
        hpx.Scale(1/hpx.Integral())
    return hpx

def compareHists(hists,Fnames, ch = "channel", reg = "region", var="sample", varname="v", WC="EFTrwgt1_cS_1_cT_1"):
    for num in range(len(hists)):
        if (hists[num].Integral() <= 0):
            return  
    Fol = 'compareHists'
    if not os.path.exists(Fol):
       os.makedirs(Fol)
    if not os.path.exists(Fol + '/' + ch):
       os.makedirs(Fol + '/' + ch)
    if not os.path.exists(Fol + '/' + ch +'/'+reg):
       os.makedirs(Fol + '/' + ch +'/'+reg)

    canvas = ROOT.TCanvas(ch+reg+var,ch+reg+var,50,50,865,780)
    canvas.SetGrid();
    canvas.SetBottomMargin(0.17)
    canvas.cd()

    legend = ROOT.TLegend(0.1,0.1,0.9,0.9)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.11)
    legend.SetNColumns(2);

    pad1=ROOT.TPad("pad1", "pad1", 0.05, 0.05, 1, 0.8 , 0)#used for the hist plot
    pad1.Draw()
    pad1.SetLogx(ROOT.kFALSE)
    pad1.SetLogy(ROOT.kFALSE)
    pad2=ROOT.TPad("pad2", "pad2", 0, 0.8, 1,1 , 0)#used for the hist plot
    pad2.Draw()
    pad2.SetLogx(ROOT.kFALSE)
    pad2.SetLogy(ROOT.kFALSE)

    pad1.cd()
    y_min=0
    y_max=0
    for num in range(0,len(hists)):
        if y_max<hists[num].GetMaximum():
            y_max=hists[num].GetMaximum()
    hists[0].SetTitle("")
    hists[0].GetYaxis().SetTitle('Fraction')
    hists[0].GetXaxis().SetLabelSize(0.03)
    hists[0].GetYaxis().SetTitleOffset(1)
    hists[0].GetYaxis().SetTitleSize(0.05)
    hists[0].GetYaxis().SetLabelSize(0.04)
    hists[0].GetYaxis().SetRangeUser(y_min,1.1*y_max)
    hists[0].GetXaxis().SetTitle(varname)
    hists[0].Draw("Hist")
    hists[0].SetLineWidth(2)
    hists[0].SetFillColor(0)
    for H in range(1,len(hists)):
        hists[H].SetLineWidth(2)
        hists[H].SetFillColor(0)
        if 'BDT' in varname:
            hists[H].GetXaxis().SetRangeUser(BDTmin, BDTmax)
            hists[H].GetXaxis().SetRangeUser(BDTmin, BDTmax)
        hists[H].Draw("HISTSAME")
    hists[0].Draw("AXISSAMEY+")
    hists[0].Draw("AXISSAMEX+")
    pad1.Update()
    pad2.cd()
    for num in range(0,len(hists)):
        legend.AddEntry(hists[num],Fnames[num],'L')
    legend.Draw("same")

    pad2.Update()
    canvas.Print(Fol + '/' + ch +'/'+reg+'/'+var+WC + ".png")
    del canvas
    gc.collect()

#year=['2016','2017','2018','All']
year=['2017']
year=['2016preVFP', '2016postVFP', '2017','2018']
regions=["llB1"]
channels=["emu"];
variables=["lep1Pt","lep1Eta","lep1Phi","lep2Pt","lep2Eta","lep2Phi","llM","llPt","llDr","llDphi","jet1Pt","jet1Eta","jet1Phi","njet","nbjet","Met","MetPhi","nVtx","llMZw", "topMass","topL1Dphi","topL1Dr","topL1DptOsumPt","topPt", "BDT"]
#variables=["lep1Pt"]
variablesName=["p_{T}(leading lepton)","#eta(leading lepton)","#Phi(leading lepton)","p_{T}(sub-leading lepton)","#eta(sub-leading lepton)","#Phi(sub-leading lepton)","M(ll)","p_{T}(ll)","#Delta R(ll)","#Delta #Phi(ll)","p_{T}(leading jet)","#eta(leading jet)","#Phi(leading jet)","Number of jets","Number of b-tagged jets","MET","#Phi(MET)","Number of vertices", "M(ll) [z window]", "top mass", "#Delta #Phi(ll, top)", "#Delta R(ll, top)", "|pt_top - pt_l1|/(pt_top + pt_l1)", "p_{T}(top)", "BDT"]

HistAddress = '/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/hists/'


Samples = [
'ttbar.root',
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
#'LFVVecU.root', 'LFVScalarU.root','LFVTensorU.root', 'LFVVecC.root', 'LFVScalarC.root','LFVTensorC.root']
#SamplesName = ['t#bar{t}', 'LFV-Vector [e#mutu]',  'LFV-Scalar [e#mutu]',  'LFV-Tensor [e#mutu]', 'LFV-Vector [e#mutc]',  'LFV-Scalar [e#mutc]',  'LFV-Tensor [e#mutc]']

#colors =  [ROOT.kRed-4,ROOT.kOrange-6, ROOT.kCyan-6,ROOT.kOrange-6,ROOT.kCyan-6,ROOT.kOrange-6, ROOT.kCyan-6]
#Style =[1,1,1,7,7,3,3]
colors =  [ROOT.kRed-4,ROOT.kBlue, ROOT.kBlue-7, ROOT.kViolet, ROOT.kViolet-5,ROOT.kAzure-9, ROOT.kAzure+10,
ROOT.kGreen, ROOT.kGreen+2, ROOT.kTeal-9, ROOT.kTeal+10, ROOT.kSpring-7, ROOT.kSpring+9]

Style =[1,1,1,1,1,1,1,2,2,2,2,2,2,2,2]
wc1 = ROOT.WCPoint("EFTrwgt1_cS_5_cT_5")
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
#                    h.SetLineStyle(Style[f])
                    h.SetLineWidth(2)
                    if 'BNV' in Samples[f]:
                        h.Scale(wc1)
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
                HHname=[]
                for f in range(len(Samples)):
                    if 'ttbar' in Samples[f]:
                        wc1 = ROOT.WCPoint("EFTrwgt1_cS_1_cT_1")
                        HH.append(EFTtoNormal(Hists[numyear][f][numch][numreg][numvar],wc1))
                        HHname.append(Samples[f][:-5])
                    else:
                        if namevar=='lep1Phi':
                            print Samples[f]+'_'+nameyear +'_'+namech +'_' + namereg+'_' +str(Hists[numyear][f][numch][numreg][numvar].Integral())
                        wc1 = ROOT.WCPoint("EFTrwgt1_cS_1_cT_1")
                        HH.append(EFTtoNormal(Hists[numyear][f][numch][numreg][numvar],wc1))
                        HHname.append(Samples[f][:-5] + "(S=1, T=1)")
                compareHists(HH,HHname, nameyear +'_'+ namech +'_' + namereg, namereg,namevar,variablesName[numvar],"EFTrwgt1_cS_1_cT_1")

for numyear, nameyear in enumerate(year):
    for numch, namech in enumerate(channels):
        for numreg, namereg in enumerate(regions):
            for numvar, namevar in enumerate(variables):
                HH=[]
                HHname=[]
                for f in range(len(Samples)):
                    if 'tbar' in Samples[f]:
                        wc1 = ROOT.WCPoint("EFTrwgt1_cS_1_cT_1")
                        HH.append(EFTtoNormal(Hists[numyear][f][numch][numreg][numvar],wc1))
                        HHname.append(Samples[f][:-5]) 
                    if 'TDUE' in Samples[f]:
                        wc1 = ROOT.WCPoint("EFTrwgt1_cS_1_cT_1")
                        AH = Hists[numyear][f][numch][numreg][numvar]
                        AH.SetLineStyle(1)
                        HH.append(EFTtoNormal(AH,wc1))
                        HHname.append(Samples[f][:-5] + "(S=1, T=1)")
                        wc1 = ROOT.WCPoint("EFTrwgt1_cS_1_cT_0")
                        AH = Hists[numyear][f][numch][numreg][numvar]
                        AH.SetLineStyle(2)
                        HH.append(EFTtoNormal(AH,wc1))
                        HHname.append(Samples[f][:-5] + "(S=1, T=0)")
                        wc1 = ROOT.WCPoint("EFTrwgt1_cS_0_cT_1")
                        AH = Hists[numyear][f][numch][numreg][numvar]
                        AH.SetLineStyle(3)
                        HH.append(EFTtoNormal(AH,wc1))
                        HHname.append(Samples[f][:-5] + "(S=0, T=1)")
                compareHists(HH,HHname, nameyear +'_'+ namech +'_' + namereg, namereg,namevar,variablesName[numvar],"STvar")


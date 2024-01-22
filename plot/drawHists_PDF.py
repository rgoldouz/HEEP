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
from operator import truediv
import copy
TGaxis.SetMaxDigits(2)



def draw(hists, sys, ch = "channel", reg = "region", year='2016', var="sample", varname="v"):
    canvas = ROOT.TCanvas(year+ch+reg+var,year+ch+reg+var,50,50,865,780)
    canvas.SetGrid();
    canvas.cd()
    hists.Draw()
    canvas.Print('sys/'+ year + '/' + ch +'/'+reg+'/'+ sys +var + ".png")
#    del legend
#    del mg
    del canvas
    gc.collect()

def compareError(histsup, ch = "channel", reg = "region", year='2016', var="sample", varname="v", prefix = 'Theory'):
    if not os.path.exists('PDF/'+year):
       os.makedirs('PDF/'+ year)
    if not os.path.exists('PDF/'+year + '/' + ch):
       os.makedirs('PDF/'+year + '/' + ch)
    if not os.path.exists('PDF/'+year + '/' + ch +'/'+reg):
       os.makedirs('PDF/'+year + '/' + ch +'/'+reg)

    canvas = ROOT.TCanvas(year+ch+reg+var,year+ch+reg+var,50,50,865,780)
    canvas.SetGrid();
    canvas.cd()

    legend = ROOT.TLegend(0.35,0.7,0.9,0.88)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.03)
    legend.SetNColumns(3);

    pad2=ROOT.TPad("pad2", "pad2", 0.0, 0.0, 1, 1 , 0)#used for the ratio plot
    pad2.Draw()
#    pad2.SetGridy()
#    pad2.SetGridx()
    pad2.SetTickx()
    pad2.SetBottomMargin(0.1)
    pad2.SetLeftMargin(0.11)
    pad2.SetRightMargin(0.1)
    pad2.SetFillStyle(0)
    pad2.SetLogx(ROOT.kFALSE)

    pad2.cd()
    maxi=0
    for n,G in enumerate(histsup):
        histsup[n].SetLineColor(n+1)
        histsup[n].SetLineWidth(2)
        histsup[n].SetFillColor(0)
        legend.AddEntry(histsup[n],str(n),'L')
        if(histsup[n].GetMaximum()>maxi):
            maxi=G.GetMaximum()
        if n==4:
            histsup[n].SetLineColor(ROOT.kOrange)
        if n==7:
            histsup[n].SetLineColor(ROOT.kGreen-1)
        if n==8:
            histsup[n].SetLineColor(28)
        if n==9:
            histsup[n].SetLineColor(46)
        if n==10:
            histsup[n].SetLineColor(30)
        if n==11:
            histsup[n].SetLineColor(38)
        if n==12:
            histsup[n].SetLineColor(17)
        if 'BDT' in varname:
            histsup[n].GetXaxis().SetRangeUser(-0.6, 0.8)
    histsup[0].SetTitle( '' )
    histsup[0].GetYaxis().SetTitle( 'Uncertainty (%)' )
    histsup[0].GetXaxis().SetTitle(varname)
    histsup[0].GetXaxis().SetLabelSize(0.04)
    histsup[0].GetYaxis().SetLabelSize(0.03)
    histsup[0].GetXaxis().SetTitleSize(0.04)
    histsup[0].GetYaxis().SetTitleSize(0.04)
    histsup[0].GetXaxis().SetTitleOffset(0.95)
    histsup[0].GetYaxis().SetTitleOffset(1)
    histsup[0].GetYaxis().SetNdivisions(804)
    histsup[0].GetXaxis().SetNdivisions(808)   
#    histsup[0].GetYaxis().SetRangeUser(-1.4*maxi,2*maxi)
    histsup[0].GetYaxis().SetRangeUser(0.75,1.25)
    histsup[0].Draw('hist')
    for n,G in enumerate(histsup):
        histsup[n].Draw('samehist')
    histsup[0].Draw('samehist')
    histsup[0].Draw("AXISSAMEY+")
    histsup[0].Draw("AXISSAMEX+")
    Lumi = '137.19'
    if (year == '2016'):
        Lumi = '35.92'
    if (year == '2017'):
        Lumi = '41.53'
    if (year == '2018'):
        Lumi = '59.74'
    label_cms="CMS Simulation Preliminary"
    Label_cms = ROOT.TLatex(0.22,0.92,label_cms)
    Label_cms.SetTextSize(0.035)
    Label_cms.SetNDC()
    Label_cms.SetTextFont(61)
    Label_cms.Draw()
    Label_lumi = ROOT.TLatex(0.65,0.92,Lumi+" fb^{-1} (13 TeV)")
    Label_lumi.SetTextSize(0.035)
    Label_lumi.SetNDC()
    Label_lumi.SetTextFont(42)
    Label_lumi.Draw("same")
    Label_channel = ROOT.TLatex(0.15,0.8,year)
    Label_channel.SetNDC()
    Label_channel.SetTextFont(42)
    Label_channel.Draw("same")

    Label_channel2 = ROOT.TLatex(0.15,0.75,ch+" ("+reg+")")
    Label_channel2.SetNDC()
    Label_channel2.SetTextFont(42)
    Label_channel2.Draw("same")

#    legend.Draw("same")
    canvas.Print('PDF/'+ year + '/' + ch +'/'+reg+'/PDF'+ prefix +'_'+var + ".png")
    del canvas
    gc.collect()


year=['2016','2017','2018','All']
#year=['2018']
LumiErr = [0.025, 0.023, 0.025, 0.018]
regions=["llB1", "llBg1"]
channels=["ee", "emu", "mumu"];
#channels=["emu"];
variables=["lep1Pt","lep1Eta","lep1Phi","lep2Pt","lep2Eta","lep2Phi","llM","llPt","llDr","llDphi","jet1Pt","jet1Eta","jet1Phi","njet","nbjet","Met","MetPhi","nVtx","llMZw","BDT","muPt","elePt","muEta","eleEta"]
#variables=["jet1Pt"]
variablesName=["p_{T}(leading lepton)","#eta(leading lepton)","#Phi(leading lepton)","p_{T}(sub-leading lepton)","#eta(sub-leading lepton)","#Phi(sub-leading lepton)","M(ll)","p_{T}(ll)","#Delta R(ll)","#Delta #Phi(ll)","p_{T}(leading jet)","#eta(leading jet)","#Phi(leading jet)","Number of jets","Number of b-tagged jets","MET","#Phi(MET)","Number of vertices", "M(ll) [z window]", "BDT output", "p_{T}(muon)","p_{T}(electron)","#eta(muon)","#eta(electron)"]
sys = ["eleRecoSf", "eleIDSf", "muIdSf", "muIsoSf", "bcTagSF", "udsgTagSF","pu", "prefiring", "trigSF", "jes", "jer", "unclusMET","muonScale","electronScale","muonRes"]
#sys = ["trigSF","eleRecoSf", "eleIDSf", "muIdSf", "muIsoSf", "unclusMET","muonScale","electronScale","muonRes"]

HistAddress = '/user/rgoldouz/NewAnalysis2020/Analysis/hists/'

Samples = ['data.root','WJetsToLNu.root','others.root', 'DY.root', 'TTTo2L2Nu.root', 'ST_tW.root', 'LFVVecC.root', 'LFVVecU.root']
SamplesName = ['Data','Jets','Others', 'DY', 't#bar{t}', 'tW' , 'LFV-vec [c_{e#mutc}] #times 100', 'LFV-vec [c_{e#mutu}] #times 10']
SamplesNameLatex = ['Data','Jets','Others', 'DY', 'tt', 'tW',  'LFV-vector(emutc)', 'LFV-vector(emutu)']
NormalizationErr = [0, 0.5, 0.5, 0.3, 0.05, 0.1, 0,0]

colors =  [ROOT.kBlack,ROOT.kYellow,ROOT.kGreen,ROOT.kBlue-3,ROOT.kRed-4,ROOT.kOrange-3, ROOT.kOrange-6, ROOT.kCyan-6]
#SignalSamples = ['LFVVecC', 'LFVVecU', 'LFVScalarC', 'LFVScalarU', 'LFVTensorC', 'LFVTensorU']
SignalSamples = ['LFVStVecU']
pdfHists = []

for f in range(len(SignalSamples)):
    for numyear, nameyear in enumerate(year):
        sysfile = ROOT.TFile.Open(HistAddress + nameyear+ '_'+ SignalSamples[f]+'.root')
        for numreg, namereg in enumerate(regions):
            for numvar, namevar in enumerate(variables):
                hNom=sysfile.Get('reweightingSys/emu' + '_' + namereg + '_' + namevar + '_QscalePDF_50')
                pdfHists=[]
                for numsys in range(50,145):
                    h=sysfile.Get('reweightingSys/emu' +'_' + namereg +  '_' + namevar + '_QscalePDF_'+str(numsys))
                    h.SetBinContent(h.GetXaxis().GetNbins(), h.GetBinContent(h.GetXaxis().GetNbins()) + h.GetBinContent(h.GetXaxis().GetNbins()+1))
                    h.Divide(hNom)
                    pdfHists.append(h)
                compareError(pdfHists, 'emu', namereg, nameyear,namevar,namevar, namevar+SignalSamples[f]) 

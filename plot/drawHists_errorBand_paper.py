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
from ROOT import TFile
from array import array
from ROOT import TColor
from ROOT import TGaxis
from ROOT import THStack
import gc
from operator import truediv
import copy
TGaxis.SetMaxDigits(4)
import random

#bins = array( 'd',[-1,-0.6,-0.4,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.35,0.6,0.8,1] )
#bins = array( 'd',[-1,-0.6,-0.5,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.20,0.25,0.4,1] )
#bins = array( 'd',[-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1] )
#bins = array( 'd',[-1,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.75,0.9,1] )
bins = array( 'd',[-1,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1] )
leptonPTbins = array( 'd',[0,25,50,75,100,125,150,175, 200, 250, 300, 350, 400, 600, 800, 1100, 1500] )
#bins = array( 'd',[-0.4,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.35,0.6] )
BDTmin=-1
BDTmax=1
Blinded=False
col=[1,632,416,600,400,616,432,800,30,
1,632,416,600,400,616,432,800,30
]
style=[1,1,1,1,1,1,1,1,1,
2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
]
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

def TH1EFTtoTH1(H, wc):
    hpx    = ROOT.TH1F( H.GetName(), H.GetName(), H.GetXaxis().GetNbins(), H.GetXaxis().GetXmin(),H.GetXaxis().GetXmax() )
    r=1
    for b in range(hpx.GetNbinsX()):
        if H.GetBinContent(b+1,ROOT.WCPoint("NONE"))>0:
            r = H.GetBinError(b+1)/H.GetBinContent(b+1,ROOT.WCPoint("NONE"))
        hpx.SetBinContent(b+1, H.GetBinContent(b+1,wc))
        hpx.SetBinError(b+1, r*H.GetBinContent(b+1,wc))
    hpx.SetLineColor(H.GetLineColor())
    hpx.SetLineStyle(H.GetLineStyle())
    return hpx


def SumofWeight(addlist):
    genEventSumw = 0
    genEventSumwScale = [0]*9
    genEventSumwPdf = [0]*100
    for add in addlist:
        for root, dirs, files in os.walk(add):
            if len(files) == 0:
                continue
            for f in files:
                filename = root + '/' + f
                if 'fail' in f:
                    continue
                fi = TFile.Open(filename)
                tree_meta = fi.Get('Runs')
                for i in range( tree_meta.GetEntries() ):
                    tree_meta.GetEntry(i)
                    genEventSumw += tree_meta.genEventSumw
                    for pdf in range(100):
                        genEventSumwPdf[pdf] += tree_meta.LHEPdfSumw[pdf]*tree_meta.genEventSumw
                    for Q in range(len(tree_meta.LHEScaleSumw)):
                        genEventSumwScale[Q] += tree_meta.LHEScaleSumw[Q]*tree_meta.genEventSumw
                tree_meta.Reset()
                tree_meta.Delete()
                fi.Close()
    if genEventSumwScale[8]==0:
        del genEventSumwScale[8]
    return [genEventSumw/x for x in genEventSumwScale] , [genEventSumw/x for x in genEventSumwPdf]


def cutFlowTable(hists, samples, regions, ch, year,caption='2016', nsig=6):
    mcSum = list(0 for i in xrange(0,len(regions))) 
#    table = '\\begin{sidewaystable*}' + "\n"
    table = '\\begin{table*}' + "\n"
    table += '\\centering' + "\n"
    table += '\\caption{' + caption +"}\n"
    table += '\\resizebox{\\textwidth}{!}{ \n'
    table += '\\begin{tabular}{|l|l|l|l|l|l|l|l|l|l|l|}' + "\n"
    table += '\\hline' + "\n"
    table += 'Samples & ' + ' & '.join(regions) + '\\\\' + "\n"
    table += '\\hline' + "\n"
    for ids, s in enumerate(samples):
        if ids==0:
            continue
        table += s 
        for idr, r in enumerate(regions):
            table += (' & ' + str(round(hists[year][ids][ch][idr][2].Integral(),2)))
            if ids<nsig:
                mcSum[idr] += hists[year][ids][ch][idr][2].Integral()
        table += '\\\\' + "\n"    
    table += '\\hline' + "\n"
    table += 'Prediction '
    for idr, r in enumerate(mcSum):
        table += (' & ' + str(round(r,2)))
    table += '\\\\' + "\n"
    table += '\\hline' + "\n"
    table += 'Data '
    for idr, r in enumerate(regions):
        table += (' & ' + str(hists[year][0][ch][idr][2].Integral()))
    table += '\\\\' + "\n"
    table += '\\hline' + "\n"
    table += 'Data$/$Pred. '
    for idr, r in enumerate(mcSum):
        table += (' & ' + str(round(hists[year][0][ch][idr][2].Integral()/r,2)))
    table += '\\\\' + "\n"
    table += '\\hline' + "\n"
    table += '\\end{tabular}}' + "\n"
    table += '\\end{table*}' + "\n"
#    table += '\\end{sidewaystable*}' + "\n"
    print table

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

def compareError(histsup,histsdown, sys, ch = "channel", reg = "region", year='2016', var="sample", varname="v", prefix = 'Theory'):
    if not os.path.exists('sys/'+year):
       os.makedirs('sys/'+ year)
    if not os.path.exists('sys/'+year + '/' + ch):
       os.makedirs('sys/'+year + '/' + ch)
    if not os.path.exists('sys/'+year + '/' + ch +'/'+reg):
       os.makedirs('sys/'+year + '/' + ch +'/'+reg)

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
        histsup[n].SetLineColor(col[n])
        histsup[n].SetLineStyle(style[n])
        histsup[n].SetLineWidth(2)
        histsup[n].SetFillColor(0)
        legend.AddEntry(histsup[n],sys[n],'L')
        if(histsup[n].GetMaximum()>maxi):
            maxi=G.GetMaximum()
        if (histsup[n].GetMaximum()>100):
            print sys[n] + 'has large error:' + str(histsup[n].GetMaximum())+  year+ch+reg+var
        histsdown[n].SetLineColor(col[n])
        histsdown[n].SetLineStyle(style[n])
        histsdown[n].SetFillColor(0)
        histsdown[n].SetLineWidth(2)
        if 'BDT' in varname:
            histsup[n].GetXaxis().SetRangeUser(BDTmin, BDTmax)
            histsdown[n].GetXaxis().SetRangeUser(BDTmin, BDTmax)
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
    histsup[0].GetYaxis().SetRangeUser(-1.4*maxi,2*maxi)
    histsup[0].Draw('hist')
    for n,G in enumerate(histsup):
        histsup[n].Draw('samehist')
        histsdown[n].Draw('samehist')
    histsup[0].Draw('samehist')
    histsdown[0].Draw('samehist')
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

    legend.Draw("same")
    canvas.Print('sys/'+ year + '/' + ch +'/'+reg+'/sys'+ prefix +'_'+var + ".png")
#    canvas.Print('sys/'+ year + '/' + ch +'/'+reg+'/sys'+ prefix +'_'+var + ".pdf")
    del canvas
    gc.collect()

def stackPlotsError(hists, SignalHists,error, errorRatio, Fnames, ch = "channel", reg = "region", year='2016', var="sample", varname="v"):
    setlog=False
    if 'BDT' in varname:
        setlog=True
    if not os.path.exists('sys/'+year):
       os.makedirs('sys/'+ year)
    if not os.path.exists('sys/'+year + '/' + ch):
       os.makedirs('sys/'+year + '/' + ch)
    if not os.path.exists('sys/'+year + '/' + ch +'/'+reg):
       os.makedirs('sys/'+year + '/' + ch +'/'+reg)
    for n,G in enumerate(hists):
        if 'BDT' in varname:
            hists[n].GetXaxis().SetRangeUser(BDTmin, BDTmax)
    for n,G in enumerate(SignalHists):
#        if 'tc' in Fnames[len(hists)+n]:
#            SignalHists[n].Scale(10)
#        if 'tu' in Fnames[len(hists)+n]:
#            SignalHists[n].Scale(1)
        if 'BDT' in varname:
            SignalHists[n].GetXaxis().SetRangeUser(BDTmin, BDTmax)

    hs = ROOT.THStack("hs","")
    for num in range(1,len(hists)):
        hs.Add(hists[num])

    dummy = hists[0].Clone()
    
    canvas = ROOT.TCanvas(year+ch+reg+var,year+ch+reg+var,50,50,865,780)
    canvas.SetGrid();
    canvas.SetBottomMargin(0.17)
    canvas.cd()

    legend = ROOT.TLegend(0.4,0.67,0.9,0.87)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.05)
    legend.SetNColumns(2);

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
    if setlog:
        pad1.SetLogy(ROOT.kTRUE)

    y_min=10
    y_max=1.5*dummy.GetMaximum()
    if setlog:
        y_max=50*dummy.GetMaximum()
    dummy.SetMarkerStyle(20)
    dummy.SetMarkerSize(1.1)
#Blinding strategy
#    if ('BDT' in var or "lep1Pt" in var or "lep2Pt" in var or "muPt" in var or "elePt" in var) and reg == 'llB1':
#        dummy.SetLineColor(0)
#        dummy.SetMarkerSize(0)
    dummy.SetTitle("")
    dummy.GetYaxis().SetTitle('Events')
    dummy.GetXaxis().SetLabelSize(0)
    dummy.GetYaxis().SetTitleOffset(0.8)
    dummy.GetYaxis().SetTitleSize(0.07)
    dummy.GetYaxis().SetLabelSize(0.05)
    dummy.GetYaxis().SetRangeUser(y_min,y_max)
    if 'BDT' in varname:
        dummy.GetXaxis().SetRangeUser(BDTmin, BDTmax)
    if 'BDT' in var and reg == 'llB1' and Blinded:
        dummy.SetLineColor(0)
        dummy.SetMarkerSize(0)
    dummy.Draw("ex0")
    hs.Draw("histSAME")
#    dummy.Draw("ex0SAME")
    for h in range(len(SignalHists)):
        SignalHists[h].SetLineWidth(2)
        SignalHists[h].SetFillColor(0)
        SignalHists[h].SetLineStyle(h+1)
        SignalHists[h].Draw("histSAME")
    if not ('BDT' in var  and reg == 'llB1' and Blinded):
        dummy.Draw("ex0SAME")

    error.SetFillColor(1)
    error.SetLineColor(1)
    error.SetFillStyle(3004)
    error.Draw("2same")
    dummy.Draw("AXISSAMEY+")
    dummy.Draw("AXISSAMEX+")  
    Lumi = '138'
    if (year == '2016preVFP'):
        Lumi = '19.52'
    if (year == '2016postVFP'):
        Lumi = '16.81'
    if (year == '2017'):
        Lumi = '41.48'
    if (year == '2018'):
        Lumi = '59.83' 
    label_cms="CMS"
    Label_cms = ROOT.TLatex(0.25,0.92,label_cms)
    Label_cms.SetNDC()
    Label_cms.SetTextSize(0.08)
    Label_cms.Draw()
    Label_cmsprelim = ROOT.TLatex(0.36,0.92,"Preliminary")
    Label_cmsprelim.SetNDC()
    Label_cmsprelim.SetTextSize(0.066)
    Label_cmsprelim.SetTextFont(51)
    Label_cmsprelim.Draw()
    Label_lumi = ROOT.TLatex(0.71,0.92,Lumi+" fb^{-1} (13 TeV)")
    Label_lumi.SetNDC()
    Label_lumi.SetTextFont(42)
    Label_lumi.SetTextSize(0.06)
    Label_lumi.Draw("same") 
    if reg=="llB1":
        Label_channel = ROOT.TLatex(0.18,0.8,ch)
        Label_channel.SetNDC()
        Label_channel.SetTextSize(0.06)
        Label_channel.SetTextFont(42)
        Label_channel.Draw("same")
        Label_channelb = ROOT.TLatex(0.18,0.74,'1 b-tagged')
        Label_channelb.SetNDC()
        Label_channelb.SetTextSize(0.06)
        Label_channelb.SetTextFont(42)
        Label_channelb.Draw("same")
    elif reg=="llBg1":
        Label_channel = ROOT.TLatex(0.18,0.8,ch)
        Label_channel.SetNDC()
        Label_channel.SetTextSize(0.06)
        Label_channel.SetTextFont(42)
        Label_channel.Draw("same")
        Label_channelb = ROOT.TLatex(0.18,0.74,'>1 b-tagged')
        Label_channelb.SetNDC()
        Label_channelb.SetTextSize(0.06)
        Label_channelb.SetTextFont(42)
        Label_channelb.Draw("same")
    else:
        Label_channel = ROOT.TLatex(0.18,0.8,ch)
        Label_channel.SetNDC()
        Label_channel.SetTextSize(0.06)
        Label_channel.SetTextFont(42)
        Label_channel.Draw("same")

    Label_channel2 = ROOT.TLatex(0.2,0.8,"#it{e#mu}")
    Label_channel2.SetNDC()
    Label_channel2.SetTextFont(42)
    Label_channel2.SetTextSize(0.06)
#    Label_channel2.Draw("same")

    legend.AddEntry(dummy,Fnames[0],'ep')
    for num in range(1,len(hists)):
#        if 'Jets' in Fnames[num] or 'DY' in Fnames[num]:
#            continue
        legend.AddEntry(hists[num],Fnames[num],'F')
#    for H in range(len(SignalHists)):
#        legend.AddEntry(SignalHists[H], Fnames[len(hists)+H],'L')
#    legend.AddEntry(None,'','')
    legend.AddEntry(error,'Stat. #oplus syst. ','F')
    legend.AddEntry(SignalHists[0], Fnames[len(hists)+0],'L')
    legend.AddEntry(None,'','')
    legend.AddEntry(SignalHists[1], Fnames[len(hists)+1],'L')
    legend.Draw("same")

    if (hs.GetStack().Last().Integral()>0):
        Label_DM = ROOT.TLatex(0.17,0.68,"Data/MC = " + str(round(hists[0].Integral()/hs.GetStack().Last().Integral(),2)))
        Label_DM.SetNDC()
        Label_DM.SetTextFont(42)
        Label_DM.Draw("same")
    pad1.Update()

    pad2.cd()
    SumofMC = hs.GetStack().Last()
    dummy_ratio = dummy.Clone()
    dummy_ratio.SetTitle("")
    dummy_ratio.GetXaxis().SetTitle(varname)
#    dummy_ratio.GetXaxis().CenterTitle()
    dummy_ratio.GetYaxis().CenterTitle()
    dummy_ratio.GetXaxis().SetMoreLogLabels()
    dummy_ratio.GetXaxis().SetNoExponent()  
    dummy_ratio.GetXaxis().SetTitleSize(0.04/0.3)
    dummy_ratio.GetYaxis().SetTitleSize(0.04/0.3)
    dummy_ratio.GetXaxis().SetTitleFont(42)
    dummy_ratio.GetYaxis().SetTitleFont(42)
    dummy_ratio.GetXaxis().SetTickLength(0.05)
    dummy_ratio.GetYaxis().SetTickLength(0.05)
    dummy_ratio.GetXaxis().SetLabelSize(0.12)
    dummy_ratio.GetYaxis().SetLabelSize(0.115)
    dummy_ratio.GetXaxis().SetLabelOffset(0.02)
    dummy_ratio.GetYaxis().SetLabelOffset(0.01)
    dummy_ratio.GetYaxis().SetTitleOffset(0.42)
    dummy_ratio.GetXaxis().SetTitleOffset(1.1)
    dummy_ratio.GetYaxis().SetNdivisions(504)    
    dummy_ratio.GetYaxis().SetRangeUser(0.6,1.4)
    if 'BDT' in varname and Blinded:
        dummy_ratio.GetXaxis().SetRangeUser(BDTmin, BDTmax)
        for b in range(dummy_ratio.GetNbinsX()):
            dummy_ratio.SetBinContent(b+1,100)
    dummy_ratio.Divide(SumofMC)
    dummy_ratio.SetStats(ROOT.kFALSE)
    dummy_ratio.GetYaxis().SetTitle('Data/Pred.')
    dummy_ratio.Draw('ex0')
    dummy_ratio.Draw("AXISSAMEY+")
    dummy_ratio.Draw("AXISSAMEX+")
    errorRatio.SetFillColor(1)
    errorRatio.SetLineColor(1)
    errorRatio.SetFillStyle(3004)
    errorRatio.Draw("2same")
    canvas.Print('sys/'+ year + '/' + ch +'/'+reg+'/'+var + ".png")
#    canvas.Print('sys/'+ year + '/' + ch +'/'+reg+'/'+var + ".pdf")

    del canvas
    gc.collect()

#year=['2016','2017','2018','All']
year=['2016preVFP', '2016postVFP', '2017','2018', 'All']
#year=['2017']
LumiErr = [0.038, 0.038, 0.038, 0.038, 0.038]
#regions=["ll","llOffZ","llB1", "llBg1"]
regions=["ll","llOffZ","llB1", "llBg1"]
regionsName=["1 b-tag", "$>$ 1 b-tag"]
channels=["ee", "emu", "mumu", "ll"];
variables=["lep1Pt","lep1Eta","lep1Phi","lep2Pt","lep2Eta","lep2Phi","llM","llPt","llDr","llDphi","jet1Pt","jet1Eta","jet1Phi","njet","nbjet","Met","MetPhi","nVtx","llMZw", "topMass","topL1Dphi","topL1Dr","topL1DptOsumPt","topPt", "BDT"]
variables=["BDT"]
variablesName=["p_{T}(leading lepton)","#eta(leading lepton)","#Phi(leading lepton)","p_{T}(sub-leading lepton)","#eta(sub-leading lepton)","#Phi(sub-leading lepton)","M(ll)","p_{T}(ll)","#Delta R(ll)","#Delta #Phi(ll)","p_{T}(leading jet)","#eta(leading jet)","#Phi(leading jet)","Number of jets","Number of b-tagged jets","MET","#Phi(MET)","Number of vertices", "M(ll) [z window]", "top mass", "#Delta #Phi(ll, top)", "#Delta R(ll, top)", "|pt_top - pt_l1|/(pt_top + pt_l1)", "p_{T}(top)", "BDT"]
variablesName=["BDT"]
sys = ["eleRecoSf", "eleIDSf", "muIdIsoSf", "bcTagSf", "LTagSf","pu", "prefiring", "trigSF","jes", "jer","muonScale","electronScale","muonRes", "unclusMET", "bcTagSfUnCorr", "LTagSfUnCorr","JetPuID","topPt"]
#sys = ["muonScale","electronScale","muonRes", "unclusMET"]

#sys = ["jes", "jer","muonScale","electronScale","muonRes", "unclusMET","JetPuID"]
#sys = ["eleRecoSf", "eleIDSf", "muIdSf", "muIsoSf", "bcTagSF", "udsgTagSF","pu", "prefiring", "trigSF", "jes", "jer","electronScale"]
#sys = ["trigSF","eleRecoSf", "eleIDSf", "muIdSf", "muIsoSf", "unclusMET","muonScale","electronScale","muonRes"]
HistAddress = '/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/hists/'

#Samples = ['data.root','WJetsToLNu.root','others.root', 'DY.root', 'TTTo2L2Nu.root', 'ST_tW.root', 'LFVVecC.root', 'LFVVecU.root']
#SamplesName = ['Data','Jets','Others', 'DY', 't#bar{t}', 'tW' , 'LFV-vec [c_{e#mutc}] #times 100', 'LFV-vec [c_{e#mutu}] #times 10']
Samples = ['data.root','other.root', 'DY.root', 'tW.root', 'ttbar.root', 'STBNV_TDUE.root', 'STBNV_TDUMu.root']
SamplesName = ['Data','Other', 'DY', 'tW','t#bar{t}', 'BNV tdue', 'BNV tdu#mu']
SamplesNameLatex = ['Data','Others', 'DY', 'tt', 'tW',  'LFV-vector(emutc)', 'LFV-vector(emutu)']
NormalizationErr = [0, 0.5, 0.3, 0.1, 0.05, 0,0]

colors =  [ROOT.kBlack,ROOT.kGreen,ROOT.kBlue-3,ROOT.kOrange-3,ROOT.kRed-4, ROOT.kGray+1,ROOT.kGray+3,
]

Hists = []
HistsSysUp = []
HistsSysDown = []
Hists_copy =[]
wc1 = ROOT.WCPoint("EFTrwgt1_cS_1_cT_1")
for numyear, nameyear in enumerate(year):
    l0=[]
    copyl0=[]
    SysUpl0=[]
    SysDownl0=[]
    Files = []
    for f in range(len(Samples)):
        l1=[]
        copyl1=[]
        SysUpl1=[]
        SysDownl1=[]
        Files.append(ROOT.TFile.Open(HistAddress + nameyear+ '_' + Samples[f]))
        for numch, namech in enumerate(channels):
            l2=[]
            copyl2=[]
            SysUpl2=[]
            SysDownl2=[]
            for numreg, namereg in enumerate(regions):
                l3=[]
                copyl3=[]
                SysUpl3=[]
                SysDownl3=[]
                for numvar, namevar in enumerate(variables):
                    SysUpl4=[]
                    SysDownl4=[]
                    h= Files[f].Get(namech + '_' + namereg + '_' + namevar)
                    if 'BNV' in Samples[f]:
                        h.Scale(wc1)
                    h.SetFillColor(colors[f])
                    h.SetLineColor(colors[f])
                    h.SetBinContent(h.GetXaxis().GetNbins(), h.GetBinContent(h.GetXaxis().GetNbins(),wc1) + h.GetBinContent(h.GetXaxis().GetNbins()+1,wc1))
                    if 'BDT' in namevar:
                        h=h.Rebin(len(bins)-1,"",bins)
                    if 'lep1Pt' in namevar:
                        h=h.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
                    l3.append(h)
                    copyl3.append(h.Clone())
                    for numsys, namesys in enumerate(sys):
                        if 'data' in Samples[f]:
                            continue
                        h= Files[f].Get('sys' + namech + '/' + namech + '_' + namereg + '_' + namevar+ '_' + namesys+ '_Up')
                        if 'BNV' in Samples[f]:
                            h.Scale(wc1)
                        h.SetFillColor(colors[f])
                        h.SetLineColor(colors[f])
                        h.SetBinContent(h.GetXaxis().GetNbins(), h.GetBinContent(h.GetXaxis().GetNbins(),wc1) + h.GetBinContent(h.GetXaxis().GetNbins()+1,wc1))
                        if 'BDT' in namevar:
                            h=h.Rebin(len(bins)-1,"",bins)
                        if 'lep1Pt' in namevar:
                            h=h.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
                        SysUpl4.append(h)
                        h= Files[f].Get('sys' + namech + '/' + namech + '_' + namereg + '_' + namevar+ '_' + namesys+ '_Down')
                        if 'BNV' in Samples[f]:
                            h.Scale(wc1)
                        h.SetFillColor(colors[f])
                        h.SetLineColor(colors[f])
                        h.SetBinContent(h.GetXaxis().GetNbins(), h.GetBinContent(h.GetXaxis().GetNbins(),wc1) + h.GetBinContent(h.GetXaxis().GetNbins()+1,wc1))
                        if 'BDT' in namevar:
                            h=h.Rebin(len(bins)-1,"",bins)
                        if 'lep1Pt' in namevar:
                            h=h.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
                        SysDownl4.append(h)
                    SysUpl3.append(SysUpl4)
                    SysDownl3.append(SysDownl4)
                l2.append(l3)
                copyl2.append(copyl3)
                SysUpl2.append(SysUpl3)
                SysDownl2.append(SysDownl3)
            l1.append(l2)
            copyl1.append(copyl2)
            SysUpl1.append(SysUpl2)
            SysDownl1.append(SysDownl2)
        l0.append(l1)
        copyl0.append(copyl1)
        SysUpl0.append(SysUpl1)
        SysDownl0.append(SysDownl1)
    Hists.append(l0)
    Hists_copy.append(copyl0)
    HistsSysUp.append(SysUpl0)       
    HistsSysDown.append(SysDownl0)


pdfGraph=[]
qscaleGraph=[]
ISRGraph=[]
FSRGraph=[]
CRGraph=[]
TuneGraph=[]
hdampGraph=[]


year2 = ['UL16preVFP','UL16postVFP','UL17','UL18']
#year2 = ['UL17', 'UL18']
for numyear, nameyear in enumerate(year):
    t1Pdf=[]
    t1Qscale=[]
    t1ISR=[]
    t1FSR=[]
    t1CR=[]
    t1Tune=[]
    t1hdamp=[]
    sysfile = ROOT.TFile.Open(HistAddress + nameyear+ '_ttbar.root')
    CR1file = ROOT.TFile.Open(HistAddress + nameyear+ '_TTTo2L2Nu_sys_CR1.root')
    CR2file = ROOT.TFile.Open(HistAddress + nameyear+ '_TTTo2L2Nu_sys_CR2.root')
    erdfile = ROOT.TFile.Open(HistAddress + nameyear+ '_TTTo2L2Nu_sys_erdON.root')
    TuneCP5upfile = ROOT.TFile.Open(HistAddress +nameyear+ '_TTTo2L2Nu_sys_TuneCP5up.root')
    TuneCP5downfile = ROOT.TFile.Open(HistAddress + nameyear+ '_TTTo2L2Nu_sys_TuneCP5down.root')
    hdampupfile = ROOT.TFile.Open(HistAddress +  nameyear+ '_TTTo2L2Nu_sys_hdampUP.root')
    hdampdownfile = ROOT.TFile.Open(HistAddress + nameyear+ '_TTTo2L2Nu_sys_hdampDOWN.root')
    if 'All' in nameyear: 
        SWscale, SWpdf =  SumofWeight(['/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL16preVFP/v2/UL16preVFP_TTTo2L2Nu', '/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL16postVFP/v2/UL16postVFP_TTTo2L2Nu', '/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL17/v2/UL17_TTTo2L2Nu', '/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL18/v2/UL18_TTTo2L2Nu'] )
    else:
        SWscale, SWpdf =  SumofWeight(['/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL' + nameyear[2:]+ '/v2/UL' + nameyear[2:]+ '_TTTo2L2Nu'])
    for numch, namech in enumerate(channels):
        tChPdf=[]
        tChQscale=[]
        tChISR=[]
        tChFSR=[]
        tChCR=[]
        tChTune=[]
        tChhdamp=[]
        for numreg, namereg in enumerate(regions):
            t2Pdf=[]
            t2Qscale=[]
            t2ISR=[]
            t2FSR=[]
            t2CR=[]
            t2Tune=[]
            t2hdamp=[]
            for numvar, namevar in enumerate(variables):
                if 'BDT' not in namevar or numreg<2:
                    continue
                pdfHists=[]
                QscaleHists=[]
                for numsys in range(9):
                    hEFT=sysfile.Get('reweightingSys/' + namech + '_' + namereg + '_' + namevar+ '_Qscale_'+str(numsys))
                    h=TH1EFTtoTH1(hEFT,wc1)
                    h.Scale(SWscale[numsys])
                    if 'BDT' in namevar:
                        h=h.Rebin(len(bins)-1,"",bins)
                    if 'lep1Pt' in namevar:
                        h=h.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
                    QscaleHists.append(h)
                for numsys in range(100):
                    hEFT=sysfile.Get('reweightingSys/' + namech +'_' + namereg + '_' + namevar+ '_PDF_'+str(numsys))
                    h=TH1EFTtoTH1(hEFT,wc1)
                    if 'BDT' in namevar:
                        h=h.Rebin(len(bins)-1,"",bins)
                    if 'lep1Pt' in namevar:
                        h=h.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
                    h.Scale(SWpdf[numsys])
                    pdfHists.append(h)
                hISRupEFT = sysfile.Get('reweightingSys/' + namech + '_' + namereg + '_' + namevar+ '_PS_2')
                hISRdownEFT = sysfile.Get('reweightingSys/' + namech + '_' + namereg + '_' + namevar+ '_PS_0')
                hFSRupEFT = sysfile.Get('reweightingSys/' + namech + '_' + namereg + '_' + namevar+ '_PS_3')
                hFSRdownEFT = sysfile.Get('reweightingSys/' + namech + '_' + namereg + '_' + namevar+ '_PS_1')
                hCR1EFT = CR1file.Get(namech + '_' + namereg + '_' + namevar)
                hCR2EFT = CR2file.Get(namech + '_' + namereg + '_' + namevar)
                herdEFT = erdfile.Get(namech + '_' + namereg + '_' + namevar)
                hTuneCP5upEFT = TuneCP5upfile.Get(namech + '_' + namereg + '_' + namevar)
                hTuneCP5downEFT = TuneCP5downfile.Get( namech + '_' + namereg + '_' + namevar)
                hhdampupEFT = hdampupfile.Get(namech + '_' + namereg + '_' + namevar)
                hhdampdownEFT = hdampdownfile.Get( namech + '_' + namereg + '_' + namevar)

                hISRup=TH1EFTtoTH1(hISRupEFT,wc1)
                hISRdown=TH1EFTtoTH1(hISRdownEFT,wc1)
                hFSRup=TH1EFTtoTH1(hFSRupEFT,wc1)
                hFSRdown=TH1EFTtoTH1(hFSRdownEFT,wc1)
                hCR1=TH1EFTtoTH1(hCR1EFT,wc1)
                hCR2=TH1EFTtoTH1(hCR2EFT,wc1)
                herd=TH1EFTtoTH1(herdEFT,wc1)
                hTuneCP5up=TH1EFTtoTH1(hTuneCP5upEFT,wc1)
                hTuneCP5down=TH1EFTtoTH1(hTuneCP5downEFT,wc1)
                hhdampup=TH1EFTtoTH1(hhdampupEFT,wc1)
                hhdampdown=TH1EFTtoTH1(hhdampdownEFT,wc1)

                if 'BDT' in namevar:
                    hISRup= hISRup.Rebin(len(bins)-1,"",bins)
                    hISRdown= hISRdown.Rebin(len(bins)-1,"",bins)
                    hFSRup=hFSRup.Rebin(len(bins)-1,"",bins)
                    hFSRdown=hFSRdown.Rebin(len(bins)-1,"",bins)
                    hCR1=hCR1.Rebin(len(bins)-1,"",bins)
                    hCR2=hCR2.Rebin(len(bins)-1,"",bins)
                    herd=herd.Rebin(len(bins)-1,"",bins)
                    hTuneCP5up=hTuneCP5up.Rebin(len(bins)-1,"",bins)
                    hTuneCP5down=hTuneCP5down.Rebin(len(bins)-1,"",bins)
                    hhdampup=hhdampup.Rebin(len(bins)-1,"",bins)
                    hhdampdown=hhdampdown.Rebin(len(bins)-1,"",bins)
                if 'lep1Pt' in namevar:
                    hISRup= hISRup.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
                    hISRdown= hISRdown.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
                    hFSRup=hFSRup.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
                    hFSRdown=hFSRdown.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
                    hCR1=hCR1.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
                    hCR2=hCR2.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
                    herd=herd.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
                    hTuneCP5up=hTuneCP5up.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
                    hTuneCP5down=hTuneCP5down.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
                    hhdampup=hhdampup.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
                    hhdampdown=hhdampdown.Rebin(len(leptonPTbins)-1,"",leptonPTbins)
                binwidth= array( 'd' )
                bincenter= array( 'd' )
                yvalue= array( 'd' )
                yerrupQscale= array( 'd' )
                yerrdownQscale= array( 'd' )
                yerrupPDF= array( 'd' )
                yerrdownPDF= array( 'd' )
                yerrupISR = array( 'd' )
                yerrdownISR = array( 'd' )
                yerrupFSR = array( 'd' )
                yerrdownFSR = array( 'd' )
                yerrupCR = array( 'd' )
                yerrdownCR = array( 'd' )
                yerrupTune= array( 'd' )
                yerrdownTune= array( 'd' )
                yerruphdamp= array( 'd' )
                yerrdownhdamp= array( 'd' )
                for b in range(QscaleHists[0].GetNbinsX()):
                    QS=np.zeros(9)
                    PDF=0
                    binwidth.append(QscaleHists[0].GetBinWidth(b+1)/2)
                    bincenter.append(QscaleHists[0].GetBinCenter(b+1))
                    yvalue.append(0)
                    nomRatio = 1
    #                if QscaleHists[0].GetBinContent(b+1) > 0:
    #                    nomRatio = 100/QscaleHists[0].GetBinContent(b+1)
                    for numsys in range(9):
                        if numsys==2 or numsys==6: 
                            continue
                        QS[numsys] = QscaleHists[numsys].GetBinContent(b+1) - QscaleHists[4].GetBinContent(b+1)
#                        print str(QscaleHists[numsys].GetBinContent(b+1)) + "***" + str(QscaleHists[0].GetBinContent(b+1)) + "***" + str(Hists[numyear][4][numch][numreg][numvar].GetBinContent(b+1,wc1))
                    yerrupQscale.append((abs(max(QS)))*nomRatio) 
                    yerrdownQscale.append((abs(min(QS)))*nomRatio)
                    for numsys in range(100):
                        PDF = PDF + (pdfHists[numsys].GetBinContent(b+1) - QscaleHists[4].GetBinContent(b+1))**2
                    yerrupPDF.append((math.sqrt(PDF))*nomRatio)
                    yerrdownPDF.append((math.sqrt(PDF))*nomRatio)
                    yerrupISR.append((abs(max(hISRup.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1) , hISRdown.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1),0)))*nomRatio)
                    yerrdownISR.append((abs(min(hISRup.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1) , hISRdown.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1),0)))*nomRatio)
                    yerrupFSR.append((abs(max(hFSRup.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1) , hFSRdown.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1,0),0)))*nomRatio)
                    yerrdownFSR.append((abs(min(hFSRup.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1) , hFSRdown.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1),0)))*nomRatio)        
                    yerrupCR.append((abs(max(herd.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1) , hCR1.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1), hCR2.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1),0)))*nomRatio)
                    yerrdownCR.append((abs(min(herd.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1) , hCR1.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1), hCR2.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1),0)))*nomRatio)
                    yerrupTune.append((abs(max(hTuneCP5up.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1) , hTuneCP5down.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1),0)))*nomRatio)
                    yerrdownTune.append((abs(min(hTuneCP5up.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1) , hTuneCP5down.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1),0)))*nomRatio)
                    yerruphdamp.append((abs(max(hhdampup.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1) , hhdampdown.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1),0)))*nomRatio)
                    yerrdownhdamp.append((abs(min(hhdampup.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1) , hhdampdown.GetBinContent(b+1)- QscaleHists[4].GetBinContent(b+1),0)))*nomRatio)
    #                del QS
    #                del PDF
    #                del nomRatio
    #                gc.collect()
                t2Qscale.append(ROOT.TGraphAsymmErrors(len(bincenter),bincenter,yvalue,binwidth,binwidth,yerrdownQscale,yerrupQscale))
                t2Pdf.append(ROOT.TGraphAsymmErrors(len(bincenter),bincenter,yvalue,binwidth,binwidth,yerrdownPDF,yerrupPDF))
                t2ISR.append(ROOT.TGraphAsymmErrors(len(bincenter),bincenter,yvalue,binwidth,binwidth,yerrdownISR,yerrupISR))
                t2FSR.append(ROOT.TGraphAsymmErrors(len(bincenter),bincenter,yvalue,binwidth,binwidth,yerrdownFSR,yerrupFSR))
                t2CR.append(ROOT.TGraphAsymmErrors(len(bincenter),bincenter,yvalue,binwidth,binwidth,yerrdownCR,yerrupCR))
                t2Tune.append(ROOT.TGraphAsymmErrors(len(bincenter),bincenter,yvalue,binwidth,binwidth,yerrdownTune,yerrupTune))
                t2hdamp.append(ROOT.TGraphAsymmErrors(len(bincenter),bincenter,yvalue,binwidth,binwidth,yerrdownhdamp,yerruphdamp))
                del binwidth
                del bincenter
                del yvalue
                del yerrupQscale
                del yerrdownQscale
                del yerrupPDF
                del yerrdownPDF
                del yerrupISR
                del yerrdownISR
                del yerrupFSR
                del yerrdownFSR
                del pdfHists
                del QscaleHists
                gc.collect()
            tChQscale.append(t2Qscale)
            tChPdf.append(t2Pdf)
            tChISR.append(t2ISR)
            tChFSR.append(t2FSR)
            tChCR.append(t2CR)
            tChTune.append(t2Tune) 
            tChhdamp.append(t2hdamp)
        t1Qscale.append(tChQscale)
        t1Pdf.append(tChPdf)
        t1ISR.append(tChISR)
        t1FSR.append(tChFSR)
        t1CR.append(tChCR)
        t1Tune.append(tChTune)
        t1hdamp.append(tChhdamp)
    pdfGraph.append(t1Pdf)
    qscaleGraph.append(t1Qscale)
    ISRGraph.append(t1ISR)
    FSRGraph.append(t1FSR)
    CRGraph.append(t1CR)
    TuneGraph.append(t1Tune)
    hdampGraph.append(t1hdamp)
    sysfile.Close()
    CR1file.Close()
    erdfile.Close()
    TuneCP5upfile.Close()
    TuneCP5downfile.Close()
    hdampupfile.Close()
    hdampdownfile.Close()
    del sysfile
    del CR1file
    del erdfile
    del TuneCP5upfile
    del TuneCP5downfile
    del hdampupfile
    del hdampdownfile
    gc.collect()

Gttsys = []
Gttsys.append(pdfGraph)
Gttsys.append(qscaleGraph)
Gttsys.append(ISRGraph)
Gttsys.append(FSRGraph)
Gttsys.append(CRGraph)
Gttsys.append(TuneGraph)
Gttsys.append(hdampGraph)
ttSys = ['pdf','QS','ISR','FSR','CR','Tune','hdamp']

for numyear, nameyear in enumerate(year):
    for numch, namech in enumerate(channels):
        for numreg, namereg in enumerate(regions):
            for numvar, namevar in enumerate(variables):
                if 'BDT' not in namevar or numreg<2:
                    continue
                glistup = []
                glistdown = []
                for g in range(len(Gttsys)):
                    hup = Hists[numyear][4][numch][numreg][numvar].Clone()
                    hdown = Hists[numyear][4][numch][numreg][numvar].Clone()
#                    print str(hup.GetNbinsX())
                    for b in range(hup.GetNbinsX()):
                        rb = 0
                        constant = 1
#                        print ttSys[g] + str(hup.GetBinContent(b+1,wc1)) +': ' + str(Gttsys[g][numyear][numch][numreg][0].GetErrorYhigh(b)) 
                        if hup.GetBinContent(b+1,wc1)>0:
                            rb = 100/hup.GetBinContent(b+1,wc1)
                        if hup.GetBinContent(b+1,wc1)<10:
                            constant = 0
                        hup.SetBinContent(b+1, 0 + Gttsys[g][numyear][numch][numreg][0].GetErrorYhigh(b)*rb*constant)
                        hdown.SetBinContent(b+1,0 - Gttsys[g][numyear][numch][numreg][0].GetErrorYlow(b)*rb*constant)
#                        print str(Gttsys[g][numyear][numch][numreg][0].GetErrorYhigh(b)*rb*constant)
                    glistup.append(hup)
                    glistdown.append(hdown)
                compareError(glistup,glistdown, ttSys, namech, namereg, nameyear,namevar,variablesName[numvar], 'ttTheory')        
#                del glist
                gc.collect()



tgraph_nominal = []
tgraph_ratio = []
errup = 0
errdown =0
for numyear, nameyear in enumerate(year):
    t1nominal = []
    t1ratio = []
    for numch, namech in enumerate(channels): 
        t2nominal = []
        t2ratio = []
        for numreg, namereg in enumerate(regions):
            t3nominal = []
            t3ratio = []
            for numvar, namevar in enumerate(variables):
                for f in range(1,len(Samples)-3):
                    Hists_copy[numyear][f+1][numch][numreg][numvar].Add(Hists_copy[numyear][f][numch][numreg][numvar])
                for numsys, namesys in enumerate(sys):
                    for f in range(1,len(Samples)-3):
                        HistsSysUp[numyear][f+1][numch][numreg][numvar][numsys].Add(HistsSysUp[numyear][f][numch][numreg][numvar][numsys]) 
                        HistsSysDown[numyear][f+1][numch][numreg][numvar][numsys].Add(HistsSysDown[numyear][f][numch][numreg][numvar][numsys])
                binwidth= array( 'f' )
                bincenter= array( 'f' )
                yvalue= array( 'f' )
                yerrup= array( 'f')
                yerrdown= array( 'f' )
                yvalueRatio= array( 'f' )
                yerrupRatio= array( 'f' )
                yerrdownRatio= array( 'f' )
                content=0
                for b in range(Hists_copy[numyear][len(Samples)-3][numch][numreg][numvar].GetNbinsX()):
                    errup = 0
                    errdown =0
                    binwidth.append(Hists_copy[numyear][len(Samples)-3][numch][numreg][numvar].GetBinWidth(b+1)/2)
                    bincenter.append(Hists_copy[numyear][len(Samples)-3][numch][numreg][numvar].GetBinCenter(b+1))
                    if Hists_copy[numyear][len(Samples)-3][numch][numreg][numvar].GetBinContent(b+1,wc1)>0:
                        content = Hists_copy[numyear][len(Samples)-3][numch][numreg][numvar].GetBinContent(b+1,wc1)
                    else:
                        content =0.0000001
                    yvalue.append(content)
                    yvalueRatio.append(content/content)
                    for numsys2, namesys2 in enumerate(sys):
                        if HistsSysUp[numyear][len(Samples)-3][numch][numreg][numvar][numsys2].Integral()==0 or HistsSysDown[numyear][len(Samples)-3][numch][numreg][numvar][numsys2].Integral()==0 or Hists_copy[numyear][len(Samples)-3][numch][numreg][numvar].GetBinContent(b+1,wc1)<=0:
                            continue
                        if HistsSysUp[numyear][len(Samples)-3][numch][numreg][numvar][numsys2].GetBinContent(b+1,wc1) - Hists_copy[numyear][len(Samples)-3][numch][numreg][numvar].GetBinContent(b+1,wc1)  > 0:
                            errup = errup + (HistsSysUp[numyear][len(Samples)-3][numch][numreg][numvar][numsys2].GetBinContent(b+1,wc1) - Hists_copy[numyear][len(Samples)-3][numch][numreg][numvar].GetBinContent(b+1,wc1))**2
                        else:
                            errdown = errdown + (HistsSysUp[numyear][len(Samples)-3][numch][numreg][numvar][numsys2].GetBinContent(b+1,wc1) - Hists_copy[numyear][len(Samples)-3][numch][numreg][numvar].GetBinContent(b+1,wc1))**2
                        if HistsSysDown[numyear][len(Samples)-3][numch][numreg][numvar][numsys2].GetBinContent(b+1,wc1) - Hists_copy[numyear][len(Samples)-3][numch][numreg][numvar].GetBinContent(b+1,wc1)  > 0:
                            errup = errup + (HistsSysDown[numyear][len(Samples)-3][numch][numreg][numvar][numsys2].GetBinContent(b+1,wc1) - Hists_copy[numyear][len(Samples)-3][numch][numreg][numvar].GetBinContent(b+1,wc1))**2
                        else:
                            errdown = errdown + (HistsSysDown[numyear][len(Samples)-3][numch][numreg][numvar][numsys2].GetBinContent(b+1,wc1) - Hists_copy[numyear][len(Samples)-3][numch][numreg][numvar].GetBinContent(b+1,wc1))**2
#statistical error
                    errup = errup + Hists_copy[numyear][len(Samples)-3][numch][numreg][numvar].GetBinError(b+1)**2
                    errdown = errdown + Hists_copy[numyear][len(Samples)-3][numch][numreg][numvar].GetBinError(b+1)**2
#Add lumi error
                    errup = errup + (LumiErr[numyear]*Hists_copy[numyear][len(Samples)-3][numch][numreg][numvar].GetBinContent(b+1,wc1))**2
                    errdown = errdown + (LumiErr[numyear]*Hists_copy[numyear][len(Samples)-3][numch][numreg][numvar].GetBinContent(b+1,wc1))**2
#                    if 'lep1Pt' in namevar:
#                        print str(math.sqrt(errup)/content) + '-' +str(math.sqrt(errdown)/content)
                  
                    for f in range(len(Samples)):
                        if 'DY' in Samples[f] and numreg==0:
                            errup = errup + (0.1*Hists[numyear][f][numch][numreg][numvar].GetBinContent(b+1,wc1))**2                    
                            errdown = errdown + (0.1*Hists[numyear][f][numch][numreg][numvar].GetBinContent(b+1,wc1))**2
                        else:
                            errup = errup + (NormalizationErr[f]*Hists[numyear][f][numch][numreg][numvar].GetBinContent(b+1,wc1))**2
                            errdown = errdown + (NormalizationErr[f]*Hists[numyear][f][numch][numreg][numvar].GetBinContent(b+1,wc1))**2
#add ttbar theory errors
                    if 'BDT' in namevar and numreg>1:
                        errup = errup + (pdfGraph[numyear][numch][numreg][0].GetErrorYhigh(b))**2
                        errup = errup + (qscaleGraph[numyear][numch][numreg][0].GetErrorYhigh(b))**2
                        errup = errup + (ISRGraph[numyear][numch][numreg][0].GetErrorYhigh(b))**2
                        errup = errup + (FSRGraph[numyear][numch][numreg][0].GetErrorYhigh(b))**2
                        errup = errup + (CRGraph[numyear][numch][numreg][0].GetErrorYhigh(b))**2
                        errup = errup + (TuneGraph[numyear][numch][numreg][0].GetErrorYhigh(b))**2
                        errup = errup + (hdampGraph[numyear][numch][numreg][0].GetErrorYhigh(b))**2
                        errdown = errdown + (pdfGraph[numyear][numch][numreg][0].GetErrorYlow(b))**2
                        errdown = errdown + (qscaleGraph[numyear][numch][numreg][0].GetErrorYlow(b))**2
                        errdown = errdown + (ISRGraph[numyear][numch][numreg][0].GetErrorYlow(b))**2
                        errdown = errdown + (FSRGraph[numyear][numch][numreg][0].GetErrorYlow(b))**2
                        errdown = errdown + (CRGraph[numyear][numch][numreg][0].GetErrorYlow(b))**2
                        errdown = errdown + (TuneGraph[numyear][numch][numreg][0].GetErrorYlow(b))**2
                        errdown = errdown + (hdampGraph[numyear][numch][numreg][0].GetErrorYlow(b))**2
                    yerrup.append(math.sqrt(errup))
                    yerrdown.append(math.sqrt(errdown))
                    yerrupRatio.append(math.sqrt(errup)/content)
                    yerrdownRatio.append(math.sqrt(errdown)/content)
                t3nominal.append(ROOT.TGraphAsymmErrors(len(bincenter),bincenter,yvalue,binwidth,binwidth,yerrdown,yerrup))
                t3ratio.append(ROOT.TGraphAsymmErrors(len(bincenter),bincenter,yvalueRatio,binwidth,binwidth,yerrdownRatio,yerrupRatio))
            t2nominal.append(t3nominal)
            t2ratio.append(t3ratio)
        t1nominal.append(t2nominal)
        t1ratio.append(t2ratio)
    tgraph_nominal.append(t1nominal)
    tgraph_ratio.append(t1ratio)


for numyear, nameyear in enumerate(year):
    for numch, namech in enumerate(channels):
        for numreg, namereg in enumerate(regions):
            for numvar, namevar in enumerate(variables):
                HH=[]
                HHsignal=[]
                for f in range(len(Samples)):
                    if 'BNV' in Samples[f]:
                        HHsignal.append(Hists[numyear][f][numch][numreg][numvar])
                    else:
                        HH.append(Hists[numyear][f][numch][numreg][numvar])
                stackPlotsError(HH, HHsignal,tgraph_nominal[numyear][numch][numreg][numvar], tgraph_ratio[numyear][numch][numreg][numvar],SamplesName, namech, namereg, nameyear,namevar,variablesName[numvar])
#                print nameyear+'_'+namech+'_'+namereg+'_'+namevar+'= '+ str(Hists[numyear][2][numch][numreg][numvar].Integral())

for numyear, nameyear in enumerate(year):
    for numch, namech in enumerate(channels):
        for numreg, namereg in enumerate(regions):
            for numvar, namevar in enumerate(variables):
                glistup = []
                glistdown = []
                for numsys2, namesys2 in enumerate(sys):
                    hup = HistsSysUp[numyear][len(Samples)-3][numch][numreg][numvar][numsys2].Clone()
                    hdown = HistsSysDown[numyear][len(Samples)-3][numch][numreg][numvar][numsys2].Clone()
                    if hup.Integral()>0 or hdown.Integral()>0:
                        for b in range(hup.GetNbinsX()):
                            cv = Hists_copy[numyear][len(Samples)-3][numch][numreg][numvar].GetBinContent(b+1,wc1)
                            rb = 0
                            if cv>0:
                                rb = 100/cv
                            hup.SetBinContent(b+1, 0 + abs(max((HistsSysUp[numyear][len(Samples)-3][numch][numreg][numvar][numsys2].GetBinContent(b+1,wc1)-cv)*rb, (HistsSysDown[numyear][len(Samples)-3][numch][numreg][numvar][numsys2].GetBinContent(b+1,wc1)-cv)*rb,0)))
                            hdown.SetBinContent(b+1, 0 - abs(min((HistsSysUp[numyear][len(Samples)-3][numch][numreg][numvar][numsys2].GetBinContent(b+1,wc1)-cv)*rb, (HistsSysDown[numyear][len(Samples)-3][numch][numreg][numvar][numsys2].GetBinContent(b+1,wc1)-cv)*rb,0)))
                    glistup.append(hup)
                    glistdown.append(hdown)
                compareError(glistup,glistdown, sys, namech, namereg, nameyear,namevar,variablesName[numvar], 'Exp')





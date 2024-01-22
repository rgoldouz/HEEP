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
TGaxis.SetMaxDigits(1)

#bins = array( 'd',[-1,-0.6,-0.4,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.35,0.6,0.8,1] )
#bins = array( 'd',[-1,-0.6,-0.4,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.3,0.6,1] )
bins = array( 'd',[-0.6,-0.4,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.20,0.26,0.6,0.8] )

def stackPlotsError(histspre, SignalHists, histspost, SignalHistspost,Fnames, ch = "channel", reg = "region", year='2016', var="sample", varname="v"):
    if not os.path.exists('sys/'+year):
       os.makedirs('sys/'+ year)
    if not os.path.exists('sys/'+year + '/' + ch):
       os.makedirs('sys/'+year + '/' + ch)
    if not os.path.exists('sys/'+year + '/' + ch +'/'+reg):
       os.makedirs('sys/'+year + '/' + ch +'/'+reg)
    for n,G in enumerate(histspre):
        if 'BDT' in varname:
            histspre[n].GetXaxis().SetRangeUser(-0.4, 0.6)
    for n,G in enumerate(SignalHists):
        if 'tc' in Fnames[len(histspre)+n]:
            SignalHists[n].Scale(10)
        if 'tu' in Fnames[len(histspre)+n]:
            SignalHists[n].Scale(1)
        if 'BDT' in varname:
            SignalHists[n].GetXaxis().SetRangeUser(-0.4, 0.6)

    hs = ROOT.THStack("hs","")
    for num in range(1,len(histspre)):
        hs.Add(histspre[num])

    dummy = histspre[0].Clone()
    
    canvas = ROOT.TCanvas(year+ch+reg+var,year+ch+reg+var,50,50,900,1000)
    canvas.SetGrid();
    canvas.SetBottomMargin(0.17)
    canvas.cd()

    legend = ROOT.TLegend(0.6,0.45,0.88,0.87)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.05)
#    legend.SetNColumns(2);

    pad1=ROOT.TPad("pad1", "pad1", 0, 0.365, 1, 0.99 , 0)#used for the hist plot
    pad2=ROOT.TPad("pad2", "pad2", 0, 0.220, 1, 0.360 , 0)#used for the ratio plot
    pad3=ROOT.TPad("pad3", "pad3", 0, 0.0, 1, 0.215 , 0)#used for the ratio plot
    pad1.Draw()
    pad2.Draw() 
    pad2.SetGridy()
    pad2.SetTickx()
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.14)
    pad1.SetRightMargin(0.05)
    pad2.SetTopMargin(0.1)
    pad2.SetBottomMargin(0.02)
    pad2.SetLeftMargin(0.14)
    pad2.SetRightMargin(0.05)
    pad2.SetFillStyle(0)
    pad3.Draw()
    pad3.SetGridy()
    pad3.SetTickx()
    pad3.SetTopMargin(0.1)
    pad3.SetBottomMargin(0.4)
    pad3.SetLeftMargin(0.14)
    pad3.SetRightMargin(0.05)
    pad3.SetFillStyle(0)
    pad1.SetFillStyle(0)
    pad1.cd()
    pad1.SetLogx(ROOT.kFALSE)
    pad2.SetLogx(ROOT.kFALSE)
    pad1.SetLogy(ROOT.kTRUE)


    y_min=10
    y_max=10000000
#    y_max=1.7*dummy.GetMaximum()
    dummy.SetMarkerStyle(20)
    dummy.SetMarkerSize(1.1)
#Blinding strategy
#    if ('BDT' in var or "lep1Pt" in var or "lep2Pt" in var or "muPt" in var or "elePt" in var) and reg == 'llB1':
#        dummy.SetLineColor(0)
#        dummy.SetMarkerSize(0)
    frame = pad1.DrawFrame(-0.4, y_min, 0.6, y_max)
    frame.SetTitle("")
    frame.GetYaxis().SetTitle('Events')
    frame.GetYaxis().SetTitleOffset(0.8)
    frame.GetYaxis().SetTitleSize(0.07)
    frame.GetYaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetRangeUser(y_min,y_max)
    frame.GetXaxis().SetLabelSize(0)
    frame.GetXaxis().SetRangeUser(-0.4,0.6)
    pad1.Update()
    dummy.Draw("ex0same")

    hs.Draw("histSAME")
    for h in range(len(SignalHists)):
        SignalHists[h].SetLineWidth(2)
        SignalHists[h].SetFillColor(0)
        SignalHists[h].SetLineStyle(h+1)
        SignalHists[h].Draw("histSAME")
#    if not (('BDT' in var or "lep1Pt" in var or "lep2Pt" in var or "muPt" in var or "elePt" in var) and reg == 'llB1'):
#        dummy.Draw("ex0SAME")

#    error.SetFillColor(13)
#    error.SetLineColor(13)
#    error.SetFillStyle(3004)
#    error.Draw("2same")
    dummy.Draw("AXISSAMEY+")
    dummy.Draw("AXISSAMEX+")    
    dummy.Draw("AXISSAMEX")
    Lumi = '137'
    if (year == '2016'):
        Lumi = '35.92'
    if (year == '2017'):
        Lumi = '41.53'
    if (year == '2018'):
        Lumi = '59.74'
    label_cms="CMS"
    Label_cms = ROOT.TLatex(0.15,0.92,label_cms)
    Label_cms.SetNDC()
    Label_cms.SetTextSize(0.08)
    Label_cms.Draw()
    Label_cmsprelim = ROOT.TLatex(0.26,0.92,"Preliminary")
    Label_cmsprelim.SetNDC()
    Label_cmsprelim.SetTextSize(0.066)
    Label_cmsprelim.SetTextFont(51)
#    Label_cmsprelim.Draw()
    Label_lumi = ROOT.TLatex(0.64,0.92,Lumi+" fb^{-1} (13 TeV)")
    Label_lumi.SetNDC()
    Label_lumi.SetTextFont(42)
    Label_lumi.SetTextSize(0.06)
    Label_lumi.Draw("same")
    if reg=="llB1":
        Label_channel = ROOT.TLatex(0.18,0.8,'e#mu')
        Label_channel.SetNDC()
        Label_channel.SetTextSize(0.06)
        Label_channel.SetTextFont(42)
        Label_channel.Draw("same")
        Label_channelb = ROOT.TLatex(0.18,0.74,'1 b-tagged')
        Label_channelb.SetNDC()
        Label_channelb.SetTextSize(0.06)
        Label_channelb.SetTextFont(42)
        Label_channelb.Draw("same")
    if reg=="llBg1":
        Label_channel = ROOT.TLatex(0.18,0.8,'e#mu')
        Label_channel.SetNDC()
        Label_channel.SetTextSize(0.06)
        Label_channel.SetTextFont(42)
        Label_channel.Draw("same")
        Label_channelb = ROOT.TLatex(0.18,0.74,'>1 b-tagged')
        Label_channelb.SetNDC()
        Label_channelb.SetTextSize(0.06)
        Label_channelb.SetTextFont(42)
        Label_channelb.Draw("same")

    Label_channel2 = ROOT.TLatex(0.2,0.8,"#it{e#mu}")
    Label_channel2.SetNDC()
    Label_channel2.SetTextFont(42)
    Label_channel2.SetTextSize(0.06)
#    Label_channel2.Draw("same")

    legend.AddEntry(dummy,Fnames[0],'ep')
    for num in range(1,len(histspre)):
        if 'Jets' in Fnames[num] or 'DY' in Fnames[num]:
            continue
        legend.AddEntry(histspre[num],Fnames[num],'F')
#    for H in range(len(SignalHists)):
#        legend.AddEntry(SignalHists[H], Fnames[len(histspre)+H],'L')
#    legend.AddEntry(None,'','')

    if (hs.GetStack().Last().Integral()>0):
        Label_DM = ROOT.TLatex(0.17,0.7,"Data/MC = " + str(round(histspre[0].Integral()/hs.GetStack().Last().Integral(),2)))
        Label_DM.SetNDC()
        Label_DM.SetTextFont(42)
#        Label_DM.Draw("same")

    for num in range(1,len(histspre)):
        for iii in range (1,histspre[num].GetSize()-2):
            if histspre[num].GetBinError(iii) > histspre[num].GetBinContent(iii):
                histspre[num].SetBinError(iii,0)

    SumofMC = histspre[1].Clone()
    for num in range(2,len(histspre)):
        SumofMC.Add(histspre[num])
    SumofMC.SetFillColor(1)
    SumofMC.SetLineColor(1)
    SumofMC.SetFillStyle(3004)
    SumofMC.Draw("e2same")
    dummy.Draw("ex0same")

    legend.AddEntry(SumofMC,'Stat. #oplus syst. ','F')
    legend.AddEntry(SignalHists[0], Fnames[len(histspre)+0],'L')
#    legend.AddEntry(None,'','')
    legend.AddEntry(SignalHists[1], Fnames[len(histspre)+1],'L')
    legend.Draw("same")
    
    pad1.Update()

    pad2.cd()
    frame2 =pad2.DrawFrame(-0.4, 0.75, 0.6, 1.25)
    pad2.Update()
    dummy_ratio = dummy.Clone()
    SumofMCE0=SumofMC.Clone()
    errorRatio=SumofMC.Clone()
    errorRatio.Divide(SumofMCE0)
    for iii in range (1,SumofMCE0.GetSize()-2):
         SumofMCE0.SetBinError(iii,0)
    frame2.SetTitle("")
    frame2.GetXaxis().SetTitle(varname)
    frame2.GetYaxis().CenterTitle()
    frame2.GetXaxis().SetMoreLogLabels()
    frame2.GetXaxis().SetNoExponent()  
    frame2.GetXaxis().SetTitleSize(0)
    frame2.GetYaxis().SetTitleSize(0.23)
    frame2.GetXaxis().SetTitleFont(42)
    frame2.GetYaxis().SetTitleFont(42)
    frame2.GetXaxis().SetTickLength(0.05)
    frame2.GetYaxis().SetTickLength(0.05)
    frame2.GetXaxis().SetLabelSize(0)
    frame2.GetYaxis().SetLabelSize(0.19)
    frame2.GetXaxis().SetLabelOffset(0.02)
    frame2.GetYaxis().SetLabelOffset(0.01)
    frame2.GetYaxis().SetTitleOffset(0.25)
    frame2.GetXaxis().SetTitleOffset(1.1)
    frame2.GetYaxis().SetNdivisions(504)    
    frame2.GetYaxis().SetRangeUser(0.85,1.15)
    if 'BDT' in varname:
        dummy_ratio.GetXaxis().SetRangeUser(-0.4, 0.6)
    dummy_ratio.Divide(SumofMCE0)
    frame2.SetStats(ROOT.kFALSE)
    frame2.GetYaxis().SetTitle('#splitline{ Data/Pred.}{ Pre-fit}')
    dummy_ratio.Draw('ex0same')
    dummy_ratio.Draw("AXISSAMEY+")
    dummy_ratio.Draw("AXISSAMEX+")
#    errorRatio.SetFillColor(13)
#    errorRatio.SetLineColor(13)
#    errorRatio.SetFillStyle(3004)
    errorRatio.Draw("e2same")
#    canvas.Print('sys/'+ year + '/' + ch +'/'+reg+'/'+var + ".png")

    for num in range(1,len(histspost)):
        for iii in range (1,histspost[num].GetSize()-2):
            if histspost[num].GetBinError(iii) > histspost[num].GetBinContent(iii):
                histspost[num].SetBinError(iii,0)
    SumofMC2 = histspost[1].Clone()
    for num in range(2,len(histspost)):
        SumofMC2.Add(histspost[num])
    SumofMC2.SetFillColor(1)
    SumofMC2.SetLineColor(1)
    SumofMC2.SetFillStyle(3004)
    pad3.cd()
    frame3 =pad3.DrawFrame(-0.4, 0.75, 0.6, 1.25)
    pad3.Update()
    dummy_ratio2 = dummy.Clone()
    SumofMCE02=SumofMC2.Clone()
    errorRatio2=SumofMC2.Clone()
    errorRatio2.Divide(SumofMCE02)
    for iii in range (1,SumofMCE02.GetSize()-2):
         SumofMCE02.SetBinError(iii,0)
    frame3.SetTitle("")
    frame3.GetXaxis().SetTitle(varname)
    frame3.GetYaxis().CenterTitle()
    frame3.GetXaxis().SetMoreLogLabels()
    frame3.GetXaxis().SetNoExponent()
    frame3.GetXaxis().SetTitleSize(0.2)
    frame3.GetYaxis().SetTitleSize(0.14)
    frame3.GetXaxis().SetTitleFont(42)
    frame3.GetYaxis().SetTitleFont(42)
    frame3.GetXaxis().SetTickLength(0.05)
    frame3.GetYaxis().SetTickLength(0.05)
    frame3.GetXaxis().SetLabelSize(0.15)
    frame3.GetYaxis().SetLabelSize(0.11)
    frame3.GetXaxis().SetLabelOffset(0.02)
    frame3.GetYaxis().SetLabelOffset(0.01)
    frame3.GetYaxis().SetTitleOffset(0.41)
    frame3.GetXaxis().SetTitleOffset(0.8)
    frame3.GetYaxis().SetNdivisions(504)
    frame3.GetYaxis().SetRangeUser(0.85,1.15)
    if 'BDT' in varname:
        dummy_ratio2.GetXaxis().SetRangeUser(-0.4, 0.6)
    dummy_ratio2.Divide(SumofMCE02)
    frame3.SetStats(ROOT.kFALSE)
    frame3.GetYaxis().SetTitle('#splitline{Data/Pred.}{Post-fit}')
    dummy_ratio2.Draw('ex0same')
    dummy_ratio2.Draw("AXISSAMEY+")
    dummy_ratio2.Draw("AXISSAMEX+")
#    errorRatio.SetFillColor(13)
#    errorRatio.SetLineColor(13)
#    errorRatio.SetFillStyle(3004)
    errorRatio2.Draw("e2same")
    canvas.Print('PostFit/'+ year + '_' + ch +'_'+reg+'_'+var + ".pdf")

    del canvas
    gc.collect()

os.system('mkdir PostFit')
year=['All']
regions=["llB1", "llBg1"]
channels=["emu"];
#variables=["lep1Pt","lep1Eta","lep1Phi","lep2Pt","lep2Eta","lep2Phi","llM","llPt","llDr","llDphi","jet1Pt","jet1Eta","jet1Phi","njet","nbjet","Met","MetPhi","nVtx","llMZw","BDT","muPt","elePt","muEta","eleEta"]
variables=["BDT"]
variablesName=["BDT discriminant"]

HistAddress = '/user/rgoldouz/NewAnalysis2020/Analysis/hists/'
#'combine/CombinedFiles_postFit/postfit_shapes.root'

#Samples = ['data.root','WJetsToLNu.root','others.root', 'DY.root', 'TTTo2L2Nu.root', 'ST_tW.root', 'LFVVecC.root', 'LFVVecU.root']
#SamplesName = ['Data','Jets','Others', 'DY', 't#bar{t}', 'tW' , 'LFV-vec [c_{e#mutc}] #times 100', 'LFV-vec [c_{e#mutu}] #times 10']
Samples = ['data_obs','Others', 'DY', 'tW','tt', 'LFVVecC.root', 'LFVVecU.root']
SamplesName = ['Data','Other', 'DY', 'tW','t#bar{t}', 'e#mutc-Vector #times 10', 'e#mutu-Vector']
SamplesNameLatex = ['Data','Other', 'DY', 'tt', 'tW',  'e#mutc-Vector #times 10', 'e#mutu-Vector']
postfitR = {}
postfitR["llB1"] = ['ch1','ch2','ch3']
postfitR["llBg1"] = ['ch4','ch5','ch6']
FitType = 'prefit'
FitType2 = 'postfit'

colors =  [ROOT.kBlack,ROOT.kGreen,ROOT.kGreen,ROOT.kOrange-3,ROOT.kRed-4, ROOT.kBlack, ROOT.kBlue]

Histsprefit = []
for numyear, nameyear in enumerate(year):
    l0=[]
    Files = []
    for f in range(len(Samples)):
        l1=[]
        if f>4: 
            Files.append(ROOT.TFile.Open(HistAddress + nameyear+ '_' + Samples[f]))
        else :
            Files.append(ROOT.TFile.Open('/user/rgoldouz/NewAnalysis2020/Analysis/combine/PostFit/LFVVecU_Combined.root'))
        for numch, namech in enumerate(channels):
            l2=[]
            for numreg, namereg in enumerate(regions):
                l3=[]
                for numvar, namevar in enumerate(variables):
                    h = ROOT.TH1F()
                    if f>4:
                        h= Files[f].Get(namech + '_' + namereg + '_' + namevar)
                        h=h.Rebin(len(bins)-1,"",bins)
                    else:
                        print postfitR[namereg][1] + '_' + FitType + '/' + Samples[f]
                       
                        h = Files[f].Get(postfitR[namereg][0] + '_' + FitType + '/' + Samples[f])
                        h1 = Files[f].Get(postfitR[namereg][1] + '_' + FitType + '/' + Samples[f])
                        h2 = Files[f].Get(postfitR[namereg][2] + '_' + FitType + '/' + Samples[f])
                        print (type(h1))
                        if type(h1) is ROOT.TH1F:
                            h.Add(h1)
                        if type(h2) is ROOT.TH1F:
                            h.Add(h2)
                    h.SetFillColor(colors[f])
                    h.SetLineColor(colors[f])
                    h.SetBinContent(h.GetXaxis().GetNbins(), h.GetBinContent(h.GetXaxis().GetNbins()) + h.GetBinContent(h.GetXaxis().GetNbins()+1))
                    l3.append(h)
                l2.append(l3)
            l1.append(l2)
        l0.append(l1)
    Histsprefit.append(l0)



Histspostfit = []
for numyear, nameyear in enumerate(year):
    l0=[]
    Files = []
    for f in range(len(Samples)):
        l1=[]
        if f>4:
            Files.append(ROOT.TFile.Open(HistAddress + nameyear+ '_' + Samples[f]))
        else :
            Files.append(ROOT.TFile.Open('/user/rgoldouz/NewAnalysis2020/Analysis/combine/PostFit/LFVVecU_Combined.root'))
        for numch, namech in enumerate(channels):
            l2=[]
            for numreg, namereg in enumerate(regions):
                l3=[]
                for numvar, namevar in enumerate(variables):
                    h = ROOT.TH1F()
                    if f>4:
                        h= Files[f].Get(namech + '_' + namereg + '_' + namevar)
                        h=h.Rebin(len(bins)-1,"",bins)
                    else:
                        print postfitR[namereg][1] + '_' + FitType2 + '/' + Samples[f]

                        h = Files[f].Get(postfitR[namereg][0] + '_' + FitType2 + '/' + Samples[f])
                        h1 = Files[f].Get(postfitR[namereg][1] + '_' + FitType2 + '/' + Samples[f])
                        h2 = Files[f].Get(postfitR[namereg][2] + '_' + FitType2 + '/' + Samples[f])
                        print (type(h1))
                        if type(h1) is ROOT.TH1F:
                            h.Add(h1)
                        if type(h2) is ROOT.TH1F:
                            h.Add(h2)
                    h.SetFillColor(colors[f])
                    h.SetLineColor(colors[f])
                    h.SetBinContent(h.GetXaxis().GetNbins(), h.GetBinContent(h.GetXaxis().GetNbins()) + h.GetBinContent(h.GetXaxis().GetNbins()+1))
                    l3.append(h)
                l2.append(l3)
            l1.append(l2)
        l0.append(l1)
    Histspostfit.append(l0)


Samples = ['data_obs','Others', 'DY', 'tW','tt', 'LFVVecC.root', 'LFVVecU.root']


for numyear, nameyear in enumerate(year):
    for numch, namech in enumerate(channels):
        for numreg, namereg in enumerate(regions):
            for numvar, namevar in enumerate(variables):
                HHpre=[]
                HHpresignal=[]
                HHpost=[]
                HHpostsignal=[]
                for f in range(len(Samples)):
                    error=ROOT.Double(0)
                    print nameyear+namech+namereg+" "+Samples[f]
                    print str(Histsprefit[numyear][f][numch][numreg][numvar].IntegralAndError(1,Histsprefit[numyear][f][numch][numreg][numvar].GetNbinsX(),error))+' +-/ '+str(error)
                    if 'LFV' in Samples[f]:
                        HHpresignal.append(Histsprefit[numyear][f][numch][numreg][numvar])
                        HHpostsignal.append(Histspostfit[numyear][f][numch][numreg][numvar])
                    else:
                        HHpre.append(Histsprefit[numyear][f][numch][numreg][numvar])
                        HHpost.append(Histspostfit[numyear][f][numch][numreg][numvar])
                stackPlotsError(HHpre, HHpresignal,HHpost, HHpostsignal,SamplesName, FitType, namereg, nameyear,namevar+FitType,variablesName[numvar])

os.system('cp PostFit/All_prefit_llB1_BDTprefit.pdf paper-fig/Figure_004-a.pdf')
os.system('cp PostFit/All_prefit_llBg1_BDTprefit.pdf paper-fig/Figure_004-b.pdf')
os.system('cp PostFit/All_postfit_llB1_BDTpostfit.pdf paper-fig/Figure_004-c.pdf')
os.system('cp PostFit/All_postfit_llBg1_BDTpostfit.pdf paper-fig/Figure_004-d.pdf')




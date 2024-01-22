#In order to run this code you need to have the postfit plots first. So run the following commands first
# combineCards.py cT_TDUE_ee_2016preVFP_llB1.txt cT_TDUE_emu_2016preVFP_llB1.txt cT_TDUE_mumu_2016preVFP_llB1.txt cT_TDUE_ee_2016postVFP_llB1.txt cT_TDUE_emu_2016postVFP_llB1.txt cT_TDUE_mumu_2016postVFP_llB1.txt cT_TDUE_ee_2017_llB1.txt cT_TDUE_emu_2017_llB1.txt cT_TDUE_mumu_2017_llB1.txt cT_TDUE_ee_2018_llB1.txt cT_TDUE_emu_2018_llB1.txt cT_TDUE_mumu_2018_llB1.txt  > cT_TDUE_2016preVFP_2016postVFP_2017_2018_com.txt
# text2workspace.py  cT_TDUE_2016preVFP_2016postVFP_2017_2018_com.txt -m 125
# combine -M FitDiagnostics  cT_TDUE_2016preVFP_2016postVFP_2017_2018_com.txt --robustFit 1  --rMin -0.15 --rMax 0.07
# PostFitShapesFromWorkspace -o postfit_shapes.root -f fitDiagnosticsTest.root:fit_s --postfit --sampling --print -d cT_TDUE_2016preVFP_2016postVFP_2017_2018_com.txt -w cT_TDUE_2016preVFP_2016postVFP_2017_2018_com.root

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
#bins = array( 'd',[-0.6,-0.4,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.20,0.26,0.6,0.8] )
bins = array( 'd',[-1,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.75,0.9,1] )

def EFTtoNormal(H, wc):
    hpx    = ROOT.TH1F( H.GetName(), H.GetName(), H.GetXaxis().GetNbins(), H.GetXaxis().GetXmin(),H.GetXaxis().GetXmax() )
    r=1
    for b in range(hpx.GetNbinsX()):
#        if H.GetBinContent(b+1,ROOT.WCPoint("NONE"))>0:
#            r = H.GetBinError(b+1)/H.GetBinContent(b+1,ROOT.WCPoint("NONE"))
        hpx.SetBinContent(b+1, H.GetBinContent(b+1,wc))
        hpx.SetBinError(b+1, H.GetBinError(b+1))
    hpx.SetLineColor(H.GetLineColor())
    hpx.SetLineStyle(H.GetLineStyle())
    if hpx.Integral()>0:
        hpx.Scale(1/hpx.Integral())
    return hpx

def TH1EFTtoTH1(H, wc):
    hpx    = ROOT.TH1F( H.GetName(), H.GetName(), H.GetXaxis().GetNbins(), H.GetXaxis().GetXmin(),H.GetXaxis().GetXmax() )
    r=1
    for b in range(hpx.GetNbinsX()):
 #       if H.GetBinContent(b+1,ROOT.WCPoint(wc))>0:
#            r = H.GetBinError(b+1)/H.GetBinContent(b+1,ROOT.WCPoint("NONE"))
        hpx.SetBinContent(b+1, H.GetBinContent(b+1,wc))
        hpx.SetBinError(b+1, H.GetBinError(b+1))
    hpx.SetLineColor(H.GetLineColor())
    hpx.SetLineStyle(H.GetLineStyle())
    return hpx


def stackPlotsError(hists, SignalHists, Fnames, ch = "channel", reg = "region", year='2016', var="sample", varname="v"):
    if not os.path.exists('sys/'+year):
       os.makedirs('sys/'+ year)
    if not os.path.exists('sys/'+year + '/' + ch):
       os.makedirs('sys/'+year + '/' + ch)
    if not os.path.exists('sys/'+year + '/' + ch +'/'+reg):
       os.makedirs('sys/'+year + '/' + ch +'/'+reg)
    for n,G in enumerate(hists):
        if 'BDT' in varname:
            hists[n].GetXaxis().SetRangeUser(-1, 1)
    for n,G in enumerate(SignalHists):
        if 'BDT' in varname:
            SignalHists[n].GetXaxis().SetRangeUser(-1, 1)

    hs = ROOT.THStack("hs","")
    for num in range(1,len(hists)):
        hs.Add(hists[num])

    dummy = hists[0].Clone()
    
#    canvas = ROOT.TCanvas(year+ch+reg+var,year+ch+reg+var,50,50,865,780)
    canvas = ROOT.TCanvas(year+ch+reg+var,year+ch+reg+var,50,50,865,600)
    canvas.SetGrid();
    canvas.SetBottomMargin(0.17)
    canvas.cd()

    legend = ROOT.TLegend(0.4,0.7,0.9,0.88)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
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
    pad1.SetLogy(ROOT.kTRUE)


    y_min=10
    y_max=100000
    if year=='All':
        y_max=200000
#    y_max=1.7*dummy.GetMaximum()
    dummy.SetMarkerStyle(20)
    dummy.SetMarkerSize(1.1)
#Blinding strategy
#    if ('BDT' in var or "lep1Pt" in var or "lep2Pt" in var or "muPt" in var or "elePt" in var) and reg == 'llB1':
#        dummy.SetLineColor(0)
#        dummy.SetMarkerSize(0)
    frame = pad1.DrawFrame(-1, y_min, 1, y_max)
    frame.SetTitle("")
    frame.GetYaxis().SetTitle('Events/bin')
    frame.GetYaxis().SetTitleOffset(0.8)
    frame.GetYaxis().SetTitleSize(0.07)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetRangeUser(y_min,y_max)
    frame.GetXaxis().SetLabelSize(0)
    frame.GetXaxis().SetRangeUser(-1,1)
    pad1.Update()
    dummy.Draw("ex0same")

    hs.Draw("histSAME")
    for H in SignalHists:
        H.SetLineWidth(2)
        H.SetFillColor(0)
        H.SetLineStyle(9)
        H.Draw("histSAME")
#    if not (('BDT' in var or "lep1Pt" in var or "lep2Pt" in var or "muPt" in var or "elePt" in var) and reg == 'llB1'):
#        dummy.Draw("ex0SAME")

#    error.SetFillColor(13)
#    error.SetLineColor(13)
#    error.SetFillStyle(3004)
#    error.Draw("2same")
    dummy.Draw("AXISSAMEY+")
    dummy.Draw("AXISSAMEX+")    
    dummy.Draw("AXISSAMEX")
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
    Label_cms = ROOT.TLatex(0.14,0.92,label_cms)
    Label_cms.SetNDC()
    Label_cms.SetTextSize(0.08)
    Label_cms.Draw()
    Label_cmsprelim = ROOT.TLatex(0.25,0.92,"Preliminary")
    Label_cmsprelim.SetNDC()
    Label_cmsprelim.SetTextSize(0.066)
    Label_cmsprelim.SetTextFont(51)
    Label_cmsprelim.Draw()
    Label_lumi = ROOT.TLatex(0.68,0.92,Lumi+" fb^{-1} (13 TeV)")
    Label_lumi.SetNDC()
    Label_lumi.SetTextFont(42)
    Label_lumi.SetTextSize(0.06)
    Label_lumi.Draw("same")
    chname='ee'
    figname='Figure_002-a'
    if ch=='emu':
        chname='e#mu'
        figname='Figure_002-c'
    if ch=='mumu':
        chname='#mu#mu'
        figname='Figure_002-b'
    if reg=="llB1":
        Label_channel = ROOT.TLatex(0.52,0.64,chname+',')
        Label_channel.SetNDC()
        Label_channel.SetTextSize(0.06)
        Label_channel.SetTextFont(42)
        Label_channel.Draw("same")
        Label_channelb = ROOT.TLatex(0.52,0.64,'     1 b-tagged')
        Label_channelb.SetNDC()
        Label_channelb.SetTextSize(0.06)
        Label_channelb.SetTextFont(42)
        Label_channelb.Draw("same")
    if reg=="llBg1":
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


    Label_channel2 = ROOT.TLatex(0.2,0.8,"#it{e#mu}")
    Label_channel2.SetNDC()
    Label_channel2.SetTextFont(42)
    Label_channel2.SetTextSize(0.06)
#    Label_channel2.Draw("same")

    legend.AddEntry(dummy,Fnames[0],'ep')
    for num in range(1,len(hists)):
        if 'Jets' in Fnames[num]:
            continue
        legend.AddEntry(hists[num],Fnames[num],'F')
#    for H in range(len(SignalHists)):
#        legend.AddEntry(SignalHists[H], Fnames[len(hists)+H],'L')
#    legend.AddEntry(None,'','')

    if (hs.GetStack().Last().Integral()>0):
        Label_DM = ROOT.TLatex(0.17,0.7,"Data/MC = " + str(round(hists[0].Integral()/hs.GetStack().Last().Integral(),2)))
        Label_DM.SetNDC()
        Label_DM.SetTextFont(42)
        if year !='All':
            Label_DM.Draw("same")

    for num in range(1,len(hists)):
        for iii in range (1,hists[num].GetSize()-2):
            if hists[num].GetBinError(iii) > hists[num].GetBinContent(iii):
                hists[num].SetBinError(iii,0)

    SumofMC = hists[1].Clone()
    for num in range(2,len(hists)):
        SumofMC.Add(hists[num])
    SumofMC.SetFillColor(13)
    SumofMC.SetLineColor(13)
    SumofMC.SetFillStyle(3244)
    SumofMC.Draw("e2same")
    dummy.Draw("ex0same")

    legend.AddEntry(SumofMC,'Stat. #oplus syst. ','F')
    legend.AddEntry(SignalHists[0], Fnames[len(hists)+0],'L')
    legend.AddEntry(SignalHists[1], Fnames[len(hists)+1],'L')
    legend.Draw("same")
    
    error=ROOT.Double(0.0)
    integral = SumofMC.IntegralAndError(0, 12, error)
    print year + '/' + ch +'/'+reg +  ' ;; ' + str(SumofMC.Integral()) + ": total MC =" + str(integral) + " +- " + str(error) + ", data=" + str(dummy.Integral()) 

    pad1.Update()

    pad2.cd()
    frame2 =pad2.DrawFrame(-1, 0.75, 1, 1.25)
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
    frame2.GetXaxis().SetTitleSize(0.04/0.3)
    frame2.GetYaxis().SetTitleSize(0.04/0.3)
    frame2.GetXaxis().SetTitleFont(42)
    frame2.GetYaxis().SetTitleFont(42)
    frame2.GetXaxis().SetTickLength(0.05)
    frame2.GetYaxis().SetTickLength(0.05)
    frame2.GetXaxis().SetLabelSize(0.115)
    frame2.GetYaxis().SetLabelSize(0.089)
    frame2.GetXaxis().SetLabelOffset(0.02)
    frame2.GetYaxis().SetLabelOffset(0.01)
    frame2.GetYaxis().SetTitleOffset(0.42)
    frame2.GetXaxis().SetTitleOffset(1.1)
    frame2.GetYaxis().SetNdivisions(504)    
    frame2.GetYaxis().SetRangeUser(0.75,1.25)
    if 'BDT' in varname:
        dummy_ratio.GetXaxis().SetRangeUser(-1, 1)
    dummy_ratio.Divide(SumofMCE0)
    frame2.SetStats(ROOT.kFALSE)
    frame2.GetYaxis().SetTitle('Data/Pred.')
    dummy_ratio.Draw('ex0same')
    dummy_ratio.Draw("AXISSAMEY+")
    dummy_ratio.Draw("AXISSAMEX+")
#    errorRatio.SetFillColor(1)
#    errorRatio.SetLineColor(13)
#    errorRatio.SetFillStyle(3004)
    errorRatio.Draw("e2same")
    dummy_ratio.Draw('ex0same')
#    canvas.Print('sys/'+ year + '/' + ch +'/'+reg+'/'+var + ".png")
    canvas.Print('PostFit/'+ year + '_' + ch +'_'+reg+'_'+var + ".pdf")
    if year=='All':
        canvas.Print('PostFit/'+figname+ ".pdf")

    del canvas
    gc.collect()

os.system('mkdir PostFit')
#year=['All']
year=['2016preVFP', '2016postVFP', '2017','2018', 'All']
regions=["llB1"]
channels=["ee", "emu", "mumu"];
#channels=["emu"];
#variables=["lep1Pt","lep1Eta","lep1Phi","lep2Pt","lep2Eta","lep2Phi","llM","llPt","llDr","llDphi","jet1Pt","jet1Eta","jet1Phi","njet","nbjet","Met","MetPhi","nVtx","llMZw","BDT","muPt","elePt","muEta","eleEta"]
variables=["BDT"]
variablesName=["BDT discriminant"]

HistAddress = '/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/hists/'
#HistAddress = '/user/rgoldouz/NewAnalysis2020/Analysis/hists/'
#'combine/CombinedFiles_postFit/postfit_shapes.root'

#Samples = ['data.root','WJetsToLNu.root','others.root', 'DY.root', 'TTTo2L2Nu.root', 'ST_tW.root', 'LFVVecC.root', 'LFVVecU.root']
#SamplesName = ['Data','Jets','Others', 'DY', 't#bar{t}', 'tW' , 'LFV-vec [c_{e#mutc}] #times 100', 'LFV-vec [c_{e#mutu}] #times 10']
Samples = ['data_obs','Others', 'DY', 'tW', 'tt', 'STBNV_TDUE.root', 'STBNV_TDUMu.root']
SamplesName = ['Data','Other', 'DY', 'tW','t#bar{t}', 'BNV tdue', 'BNV tdu#mu']
SamplesNameLatex = ['Data','Others', 'DY', 'tt', 'tW',  'LFV-vector(emutc)', 'LFV-vector(emutu)']

postfitR = {}
postfitR["ee"] = ['ch1','ch4','ch7','ch10']
postfitR["emu"] = ['ch2','ch5','ch8','ch11']
postfitR["mumu"] = ['ch3','ch6','ch9','ch12']
#FitType = 'prefit'
FitType = 'postfit'

colors =  [ROOT.kBlack,ROOT.kGreen,ROOT.kBlue-3,ROOT.kOrange-3,ROOT.kRed-4, ROOT.kGray+1,ROOT.kGray+3,]
#colors =  [ROOT.kBlack,ROOT.kGreen,ROOT.kGreen,ROOT.kOrange-3,ROOT.kRed-4, ROOT.kOrange-6, ROOT.kCyan-6]
wc1 = ROOT.WCPoint("EFTrwgt1_cT_1_cS_1")
Hists = []
for numyear, nameyear in enumerate(year):
    l0=[]
    Files = []
    for f in range(len(Samples)):
        l1=[]
        if f>4: 
            Files.append(ROOT.TFile.Open(HistAddress + nameyear+ '_' + Samples[f]))
        else :
            Files.append(ROOT.TFile.Open('/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/plot/postfit_shapes.root'))
        for numch, namech in enumerate(channels):
            l2=[]
            for numreg, namereg in enumerate(regions):
                l3=[]
                for numvar, namevar in enumerate(variables):
                    h = ROOT.TH1F()
                    if f>4:
                        hEFT= Files[f].Get(namech + '_' + namereg + '_' + namevar)
                        hEFT.Scale(wc1)
                        hNormal = TH1EFTtoTH1(hEFT,wc1)
                        h=hNormal.Rebin(len(bins)-1,"",bins)
                    else:
                        if nameyear=='All':
                            h = Files[f].Get(postfitR[namech][0] + '_' + FitType + '/' + Samples[f])
                            h1 = Files[f].Get(postfitR[namech][1] + '_' + FitType + '/' + Samples[f])
                            h2 = Files[f].Get(postfitR[namech][2] + '_' + FitType + '/' + Samples[f])
                            h3 = Files[f].Get(postfitR[namech][3] + '_' + FitType + '/' + Samples[f])
                            if type(h1) is ROOT.TH1F:
                                h.Add(h1)
                            if type(h2) is ROOT.TH1F:
                                h.Add(h2)
                            if type(h3) is ROOT.TH1F:
                                h.Add(h3)
                        else:
                            h = Files[f].Get(postfitR[namech][numyear] + '_' + FitType + '/' + Samples[f])
                    if type(h) is not ROOT.TH1F:
                        h = Files[f].Get(postfitR[namech][numyear] + '_' + FitType + '/' + Samples[4])
                        h.Scale(0.00001)
                    h.SetFillColor(colors[f])
                    h.SetLineColor(colors[f])
                    h.SetBinContent(h.GetXaxis().GetNbins(), h.GetBinContent(h.GetXaxis().GetNbins()) + h.GetBinContent(h.GetXaxis().GetNbins()+1))
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
                    if 'BNV' in Samples[f]:
                        HHsignal.append(Hists[numyear][f][numch][numreg][numvar])
                    else:
                        HH.append(Hists[numyear][f][numch][numreg][numvar])
                stackPlotsError(HH, HHsignal,SamplesName, namech,namereg, nameyear,namevar+FitType,variablesName[numvar])


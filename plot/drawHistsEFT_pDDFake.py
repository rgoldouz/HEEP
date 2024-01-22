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
TGaxis.SetMaxDigits(2)

def cutFlowTable(hists, samples, regions, ch, year,caption='2016', nsig=6):
    mcSum = list(0 for i in xrange(0,len(regions)))
    for ids, s in enumerate(samples):
        if ids==0:
            continue
        for idr, r in enumerate(regions):
            if ids<nsig:
                mcSum[idr] += hists[year][ids][ch][idr][2].Integral() 
#    table = '\\begin{sidewaystable*}' + "\n"
    table = '\\begin{table*}' + "\n"
    table += '\\centering' + "\n"
    table += '\\caption{' + caption +": Number of expected signal and background events, compared to the event yields in the data, after various selection steps. Percentage event fractions of the MC predictions are given in brackets.}\n"
#    table += '\\resizebox{\\textwidth}{!}{ \n'
    table += '\\begin{tabular}{|l|l|l|l|l|l|l|l|l|l|l|}' + "\n"
    table += '\\hline' + "\n"
    table += 'Samples & ' + ' & '.join(regions) + '\\\\' + "\n"
    table += '\\hline' + "\n"
    for ids, s in enumerate(samples):
        if ids==0:
            continue
        table += s 
        for idr, r in enumerate(regions):
            if ids<nsig:
                table += (' & ' + str(round(hists[year][ids][ch][idr][2].Integral(),2)) + '[' + str(round((100*hists[year][ids][ch][idr][2].Integral())/mcSum[idr],2)) +'\%]')
            else:
                table += (' & ' + str(round(hists[year][ids][ch][idr][2].Integral(),2)))
#            if hists[year][ids][ch][idr][2].Integral()>0:
#                print s+' ***********stat Error:' +str(math.sqrt(hists[year][ids][ch][idr][2].GetSumw2().GetSum())/hists[year][ids][ch][idr][2].Integral()) 
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
    table += '\\end{tabular}' + "\n"
    table += '\\end{table*}' + "\n"
#    table += '\\end{sidewaystable*}' + "\n"
    print table

def stackPlots(hists, SignalHists, Fnames,FnamesS, ch = "channel", reg = "region", year='2016', var="sample", varname="v"):
    if not os.path.exists(year):
       os.makedirs(year)
    if not os.path.exists(year + '/' + ch):
       os.makedirs(year + '/' + ch)
    if not os.path.exists(year + '/' + ch +'/'+reg):
       os.makedirs(year + '/' + ch +'/'+reg)

    for n,G in enumerate(hists):
        if 'BDT' in varname:
            hists[n].GetXaxis().SetRangeUser(-0.5, 0.7)
    for n,G in enumerate(SignalHists):
        if 'BDT' in varname:
            SignalHists[n].GetXaxis().SetRangeUser(-0.5, 0.7)

    hs = ROOT.THStack("hs","")
#    for num in range(len(hists)):
#        hists[num].SetBinContent(hists[num].GetXaxis().GetNbins(), hists[num].GetBinContent(hists[num].GetXaxis().GetNbins()) + hists[num].GetBinContent(hists[num].GetXaxis().GetNbins()+1))
#    for num in range(len(SignalHists)):
 #       SignalHists[num].SetBinContent(SignalHists[num].GetXaxis().GetNbins(),SignalHists[num].GetBinContent(SignalHists[num].GetXaxis().GetNbins()) + SignalHists[num].GetBinContent(SignalHists[num].GetXaxis().GetNbins()+1))
    for num in range(1,len(hists)):
        hs.Add(hists[num])

    dummy = hists[0].Clone()

    
    canvas = ROOT.TCanvas(year+ch+reg+var,year+ch+reg+var,50,50,865,780)
    canvas.SetGrid();
    canvas.SetBottomMargin(0.17)
    canvas.cd()

    legend = ROOT.TLegend(0.7,0.45,0.9,0.88)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.035)

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

    y_min=0
    y_max=1.6*dummy.GetMaximum()
    dummy.SetMarkerStyle(20)
    dummy.SetMarkerSize(1.2)
    dummy.SetTitle("")
    dummy.GetYaxis().SetTitle('Events')
    dummy.GetXaxis().SetLabelSize(0)
    dummy.GetYaxis().SetTitleOffset(0.8)
    dummy.GetYaxis().SetTitleSize(0.07)
    dummy.GetYaxis().SetLabelSize(0.04)
    dummy.GetYaxis().SetRangeUser(y_min,y_max)
    if 'BDT' in varname:
        dummy.GetXaxis().SetRangeUser(-0.5, 0.7)
    if 'BDT' in var and reg == 'llB1':
        dummy.SetLineColor(0)
        dummy.SetMarkerSize(0)
    dummy.Draw("e")
    hs.Draw("histSAME")
    for h in range(len(SignalHists)):
        SignalHists[h].SetLineWidth(2)
        SignalHists[h].SetFillColor(0)
        SignalHists[h].SetLineStyle(h+1)
        SignalHists[h].Draw("histSAME")
    if 'BDT' not in varname:
        dummy.Draw("eSAME")
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
    label_cms="CMS Preliminary"
    Label_cms = ROOT.TLatex(0.2,0.92,label_cms)
    Label_cms.SetNDC()
    Label_cms.SetTextFont(61)
    Label_cms.Draw()
    Label_lumi = ROOT.TLatex(0.71,0.92,Lumi+" fb^{-1} (13 TeV)")
    Label_lumi.SetNDC()
    Label_lumi.SetTextFont(42)
    Label_lumi.Draw("same")
    Label_channel = ROOT.TLatex(0.2,0.8,year +" / "+ch+" ("+reg+")")
    Label_channel.SetNDC()
    Label_channel.SetTextFont(42)
    Label_channel.Draw("same")


    legend.AddEntry(dummy,Fnames[0],'ep')
    for num in range(1,len(hists)):
        legend.AddEntry(hists[num],Fnames[num],'F')
    for H in range(len(SignalHists)):
        legend.AddEntry(SignalHists[H], FnamesS[H],'L')
    legend.Draw("same")

    if (hs.GetStack().Last().Integral()>0):
        Label_DM = ROOT.TLatex(0.2,0.75,"Data/MC = " + str(round(hists[0].Integral()/hs.GetStack().Last().Integral(),2)))
        Label_DM.SetNDC()
        Label_DM.SetTextFont(42)
        Label_DM.Draw("same")

    pad1.Update()

    pad2.cd()
    SumofMC = hs.GetStack().Last()
    dummy_ratio = hists[0].Clone()
    dummy_ratio.SetTitle("")
    dummy_ratio.SetMarkerStyle(20)
    dummy_ratio.SetMarkerSize(1.2)
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
    dummy_ratio.GetXaxis().SetLabelSize(0.115)
    dummy_ratio.GetYaxis().SetLabelSize(0.089)
    dummy_ratio.GetXaxis().SetLabelOffset(0.02)
    dummy_ratio.GetYaxis().SetLabelOffset(0.01)
    dummy_ratio.GetYaxis().SetTitleOffset(0.42)
    dummy_ratio.GetXaxis().SetTitleOffset(1.1)
    dummy_ratio.GetYaxis().SetNdivisions(504)    
    dummy_ratio.GetYaxis().SetRangeUser(0.8,1.2)
    dummy_ratio.Divide(SumofMC)
    dummy_ratio.SetStats(ROOT.kFALSE)
    dummy_ratio.GetYaxis().SetTitle('Data/Pred.')
    if 'BDT' in varname:
        dummy_ratio.GetXaxis().SetRangeUser(-0.5, 0.7)
    dummy_ratio.Draw("AXISSAMEY+")
    dummy_ratio.Draw("AXISSAMEX+")
    if 'BDT' not in varname:
        dummy_ratio.Draw()
    canvas.Print(year + '/' + ch +'/'+reg+'/'+var + ".png")
    del canvas
    gc.collect()


#year=['2016','2017','2018','All']
#year=['2016preVFP', '2016postVFP', '2017','2018','All']
year=['2017']
regions=["ll","llB1", "llBg1"]
#regions=["ll","llOffZ"]
#regions=["ll","llB1", "llBg1"]
regionsName=["ll","llOffZ","llB1", "llBg1"]
channels=[
"OSee", "OSemu", "OSmumu","SSee", "SSemu", "SSmumu",
"SSefe","SSfefe","SSfemu"
#"EpEm", "MUpMUm", "EpmMUmp","LLpp","LLmm", "3LonZ", "3LoffZp", "3LoffZm","4L","LLOSpp","LLOSmm",
#                               "LFpp", "FFpp", "LFmm", "FFmm",
#                                "LLFonZ", "LFFonZ","FFFonZ",
#                                "LLFoffZp", "LFFoffZp","FFFoffZp",
#                                "LLFoffZm", "LFFoffZm","FFFoffZm"
];
channelsFake=["SSefe","SSfefe","SSfemu"]
variables=["lep1Pt","lep1Eta","lep1Phi","lep2Pt","lep2Eta","lep2Phi","llM","llPt","llDr","llDphi","jet1Pt","jet1Eta","jet1Phi","njet","nbjet","Met","MetPhi","nVtx","llMZw"]
#variables=["lep1Pt"]
variablesName=["p_{T}(leading lepton)","#eta(leading lepton)","#Phi(leading lepton)","p_{T}(sub-leading lepton)","#eta(sub-leading lepton)","#Phi(sub-leading lepton)","M(ll)","p_{T}(ll)","#Delta R(ll)","#Delta #Phi(ll)","p_{T}(leading jet)","#eta(leading jet)","#Phi(leading jet)","Number of jets","Number of b-tagged jets","MET","#Phi(MET)","Number of vertices", "M(ll) [z window]"]

HistAddress = '/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/hists/'


#Samples = ['data.root','Triboson.root', 'Diboson.root', 'ttbar.root', 'ST.root','DY.root', 'Conv.root','TTX.root','WJets.root','FCNCProduction.root','FCNCDecay.root']
#SamplesName = ['Data','Triboson', 'Diboson', 't#bar{t}', 'Single top','DY', 'Conv','TTX+TX', 'WJets','FCNC-Production','FCNC-Decay']# , 'BNV_ST_TBCE', 'BNV_ST_TBUE', 'BNV_ST_TDCE',  'BNV_ST_TDUE',  'BNV_ST_TSCE',  'BNV_ST_TSUE']
Samples = ['data.root','WJets.root','other.root', 'DY.root', 'ttbar.root', 'tW.root']
SamplesName = ['Data','Jets','Others', 'DY', 't#bar{t}', 'tW','BNV tdue', 'BNV tdu#mu']# , 'BNV_ST_TBCE', 'BNV_ST_TBUE', 'BNV_ST_TDCE',  'BNV_ST_TDUE',  'BNV_ST_TSCE',  'BNV_ST_TSUE']
SamplesNameLatex = ['Data','Jets','Others', 'DY', 'tt', 'tW', 'BNV tdue', 'BNV tdu#mu']#,  'LFV-vector(emutc)', 'LFV-vector(emutu)']

colors =  [ROOT.kBlack,ROOT.kYellow,ROOT.kGreen,ROOT.kBlue-3,ROOT.kRed-4,ROOT.kOrange-3, ROOT.kPink, ROOT.kViolet-1,ROOT.kViolet+8,ROOT.kRed-5,ROOT.kSpring-1,ROOT.kGray+3]


#Samples = ['data.root','Triboson.root', 'WZTo3LNu.root', 'ZZTo4L.root', 'WWTo2L2Nu.root', 'ttbar.root', 'ST.root','DY.root', 'Conv.root','TTX.root','WJets.root']
#SamplesName = ['Data','Triboson', 'WZ','ZZ','WW', 't#bar{t}', 'Single top','DY', 'Conv','TTX+TX', 'WJets']# , 'BNV_ST_TBCE', 'BNV_ST_TBUE', 'BNV_ST_TDCE',  'BNV_ST_TDUE',  'BNV_ST_TSCE',  'BNV_ST_TSUE']
#colors =  [ROOT.kBlack,ROOT.kGreen,ROOT.kBlue+8,ROOT.kBlue+3,ROOT.kBlue-3,ROOT.kRed-4,ROOT.kOrange-3, ROOT.kCyan+1, ROOT.kYellow, ROOT.kMagenta-4,ROOT.kYellow+2]



wc1 = ROOT.WCPoint("EFTrwgt1_cS_1_cT_1")

Hists = []
HistsFake = []
for numyear, nameyear in enumerate(year):
    l0=[]
    Files = []
    for f in range(len(Samples)):
        l1=[]
        Files.append(ROOT.TFile.Open(HistAddress + nameyear+ '_' + Samples[f]))
        print HistAddress + nameyear+ '_' + Samples[f]
        for numch, namech in enumerate(channels):
            l2=[]
            for numreg, namereg in enumerate(regions):
                l3=[]
                for numvar, namevar in enumerate(variables):
#                    print namech + '_' + namereg + '_' + namevar
                    h= Files[f].Get(namech + '_' + namereg + '_' + namevar)
#                    print namevar + ":" + str(h.Integral())
#                    if 'njet' in namevar:
#                        print namech + '_' + namereg + '_' + namevar + str(h.GetBinFit(3).getDim())
#                        print namech + '_' + namereg + '_' + namevar + str(h.GetBinContent(3,wc1))
                    h.SetFillColor(colors[f])
                    h.SetLineColor(colors[f])
                    l3.append(h)
                l2.append(l3)
            l1.append(l2)
        l0.append(l1)
    Hists.append(l0)       

for numyear, nameyear in enumerate(year):
    l0=[]
    Files = []
    for f in range(1):
        l1=[]
        Files.append(ROOT.TFile.Open(HistAddress + nameyear+ '_' + Samples[f]))
#        print HistAddress + nameyear+ '_' + Samples[f]
        for numch, namech in enumerate(channelsFake):
            l2=[]
            for numreg, namereg in enumerate(regions):
                l3=[]
                for numvar, namevar in enumerate(variables):
                    h= Files[f].Get(namech + '_' + namereg + '_' + namevar)
#                    print namevar + ":" + str(h.Integral())
#                    if 'njet' in namevar:
#                        print namech + '_' + namereg + '_' + namevar + str(h.GetBinFit(3).getDim())
#                        print namech + '_' + namereg + '_' + namevar + str(h.GetBinContent(3,wc1))
                    h.SetFillColor(ROOT.kYellow+2)
                    h.SetLineColor(ROOT.kYellow+2)
                    l3.append(h)
                l2.append(l3)
            l1.append(l2)
        l0.append(l1)
    HistsFake.append(l0)

for numyear, nameyear in enumerate(year):
    for numch, namech in enumerate(channels):
        for numreg, namereg in enumerate(regions):
            for numvar, namevar in enumerate(variables):
                HH=[]
                HHsignal=[]
                SN=[]
                SNsignal=[]
                for f in range(len(Samples)):
                    Hists[numyear][f][numch][numreg][numvar]
                    HH.append(Hists[numyear][f][numch][numreg][numvar])
                    SN.append(SamplesName[f])
                stackPlots(HH, HHsignal, SN, SNsignal, namech, namereg, nameyear,namevar,variablesName[numvar])
#                stackPlots(HH, HHsignal, SamplesName, namech, namereg, nameyear,namevar,variablesName[numvar])
    os.system('tar -cvf '+ nameyear +'DD.tar ' + nameyear)

#"OSee", "OSemu", "OSmumu","SSee", "SSemu", "SSmumu","SSefe","SSfefe","SSfemu"
#for numyear, nameyear in enumerate(year):
#    for numreg, namereg in enumerate(regions):
#        for numvar, namevar in enumerate(variables):
#            HistsFake[numyear][0][channelsFake.index('SSefe')][numreg][numvar].Add(HistsFake[numyear][0][channelsFake.index('')][numreg][numvar],-1)
#            HistsFake[numyear][0][channelsFake.index('LFmm')][numreg][numvar].Add(HistsFake[numyear][0][channelsFake.index('FFmm')][numreg][numvar],-1)
#
#fakeMap={
#"LLpp":"LFpp",
#"LLmm":"LFmm", 
#"3LonZ":"LLFonZ", 
#"3LoffZp":"LLFoffZp", 
#"3LoffZm":"LLFoffZm"
#}
#
#
#
#for numyear, nameyear in enumerate(year):
#    for numch, namech in enumerate(channels):
#        for numreg, namereg in enumerate(regions):
#            for numvar, namevar in enumerate(variables):
#                HH=[]
#                HHsignal=[]
#                SN=[]
#                SNsignal=[]
#                for f in range(len(Samples)):
##                    if namech in fakeMap and Samples[f] in fakeMC:
##                        continue
#                    if '3L' in namech and 'DY' in SamplesName[f]:
#                        continue
#                    if 'FCNC' in Samples[f]:
#                        text = 'EFTrwgt4_cpQM_1.0_cpt_1.0_ctA_1.0_ctZ_0.5_ctG_0.1_cQlM_1.0_cQe_1.0_ctl_1.0_cte_1.0_ctlS_1.0_ctlT_0.05_ctp_1.0'
#                        wc1 = ROOT.WCPoint(text)
#                        Hists[numyear][f][numch][numreg][numvar].Scale(wc1)
#                        HHsignal.append(Hists[numyear][f][numch][numreg][numvar])
#                        SNsignal.append(SamplesName[f])
#                    else:
#                        HH.append(Hists[numyear][f][numch][numreg][numvar])
#                        SN.append(SamplesName[f])
#                if namech in fakeMap:
#                    HH.append(HistsFake[numyear][0][channelsFake.index(fakeMap[namech])][numreg][numvar])
#                    SN.append("Fake")
#                if namech=="LLpp":
#                    Hists[numyear][0][channels.index("LLOSpp")][numreg][numvar].SetFillColor(ROOT.kSpring+1)
#                    Hists[numyear][0][channels.index("LLOSpp")][numreg][numvar].SetLineColor(ROOT.kSpring+1)
#                    HH.append(Hists[numyear][0][channels.index("LLOSpp")][numreg][numvar])
#                    SN.append("ChargeFlip")
#                if namech=="LLmm":
#                    Hists[numyear][0][channels.index("LLOSmm")][numreg][numvar].SetFillColor(ROOT.kSpring+1)
#                    Hists[numyear][0][channels.index("LLOSmm")][numreg][numvar].SetLineColor(ROOT.kSpring+1)
#                    HH.append(Hists[numyear][0][channels.index("LLOSmm")][numreg][numvar])
#                    SN.append("ChargeFlip")
#                stackPlots(HH, HHsignal, SN, SNsignal, namech, namereg, nameyear,namevar,variablesName[numvar])
#
#le = '\\documentclass{article}' + "\n"
#le += '\\usepackage{rotating}' + "\n"
#le += '\\usepackage{rotating}' + "\n"
#le += '\\begin{document}' + "\n"
#
#print le
##for numyear, nameyear in enumerate(year):
##    for numch, namech in enumerate(channels):
##        cutFlowTable(Hists, SamplesNameLatex, regionsName, numch, numyear, nameyear + ' ' + namech, 6 )
#print '\\end{document}' + "\n"



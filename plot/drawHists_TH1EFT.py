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

def EFTtoNormal(H, wc):
    hpx    = ROOT.TH1F( H.GetName(), H.GetName(), H.GetXaxis().GetNbins(), H.GetXaxis().GetXmin(),H.GetXaxis().GetXmax() )
    r=1
    for b in range(hpx.GetNbinsX()):               
        if H.GetBinContent(b+1,ROOT.WCPoint("NONE"))>0:
            r = H.GetBinError(b+1)/H.GetBinContent(b+1,ROOT.WCPoint("NONE")) 
        hpx.SetBinContent(b+1, H.GetBinContent(b+1,wc))
        hpx.SetBinError(b+1, r*H.GetBinContent(b+1,wc)) 
    if hpx.Integral()>0:
        hpx.Scale(1/hpx.Integral())
    return hpx

def stackPlots(hists, SignalHists, Fnames, ch = "channel", reg = "region", year='2016', var="sample", varname="v"):
    if not os.path.exists(year):
       os.makedirs(year)
    if not os.path.exists(year + '/' + ch):
       os.makedirs(year + '/' + ch)
    if not os.path.exists(year + '/' + ch +'/'+reg):
       os.makedirs(year + '/' + ch +'/'+reg)
    hs = ROOT.THStack("hs","")
    for num in range(len(hists)):
        hists[num].SetBinContent(hists[num].GetXaxis().GetNbins(), hists[num].GetBinContent(hists[num].GetXaxis().GetNbins()) + hists[num].GetBinContent(hists[num].GetXaxis().GetNbins()+1))
    for num in range(len(SignalHists)):
        SignalHists[num].SetBinContent(SignalHists[num].GetXaxis().GetNbins(),SignalHists[num].GetBinContent(SignalHists[num].GetXaxis().GetNbins()) + SignalHists[num].GetBinContent(SignalHists[num].GetXaxis().GetNbins()+1))
    for num in range(1,len(hists)):
        hs.Add(hists[num])

    dummy = hists[0].Clone()

    
    canvas = ROOT.TCanvas(year+ch+reg+var,year+ch+reg+var,50,50,865,780)
    canvas.SetGrid();
    canvas.SetBottomMargin(0.17)
    canvas.cd()

    legend = ROOT.TLegend(0.7,0.55,0.9,0.88)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)

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
    dummy.Draw("e")
    hs.Draw("histSAME")
    for H in SignalHists:
        H.SetLineWidth(2)
        H.SetFillColor(0)
        H.SetLineStyle(9)
        H.Draw("histSAME")
    dummy.Draw("eSAME")
    dummy.Draw("AXISSAMEY+")
    dummy.Draw("AXISSAMEX+")

    Lumi = '137.42'
    if (year == '2016'):
        Lumi = '35.92'
    if (year == '2017'):
        Lumi = '41.53'
    if (year == '2018'):
        Lumi = '59.97'
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
        legend.AddEntry(SignalHists[H], Fnames[len(hists)+H],'L')
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
    dummy_ratio.GetYaxis().SetRangeUser(0,2)
    dummy_ratio.Divide(SumofMC)
    dummy_ratio.SetStats(ROOT.kFALSE)
    dummy_ratio.GetYaxis().SetTitle('Data/Pred.')
    dummy_ratio.Draw()
    dummy_ratio.Draw("AXISSAMEY+")
    dummy_ratio.Draw("AXISSAMEX+")
    canvas.Print(year + '/' + ch +'/'+reg+'/'+var + ".png")
    del canvas
    gc.collect()

def drawHistallW(A,AN, R, I, Xlabel):
    canvas = ROOT.TCanvas(R+I,R+I,10,10,1100,628)
    canvas.cd()

    pad_name = "pad"
    pad1=ROOT.TPad(pad_name, pad_name, 0, 0, 1, 1 , 0)
    pad1.Draw()
    pad1.cd()
#    pad1.SetLogy()
    A[0].SetTitle('')
    A[0].GetXaxis().SetTitle(R)
    A[0].GetYaxis().SetTitle(Xlabel)

    maxi=A[0].GetMaximum()
    mini=A[0].GetMinimum()
    for i in range(len(A)):
        A[i].SetLineColor(i+1)
        if (i+1)==5:
            A[i].SetLineColor(28)
        if A[i].GetMaximum()>maxi:
            maxi=A[i].GetMaximum()
        if A[i].GetMinimum()<mini:
            mini=A[i].GetMinimum()
        if i==9:
            A[i].SetMarkerColor(12);

    A[0].SetMinimum(mini);
    A[0].SetMaximum(maxi * 1.2);

    leg = ROOT.TLegend(0.65,0.70,0.88,0.88)

    for i in range(len(A)):
        A[i].Draw("HIST SAME")
        leg.AddEntry(A[i], AN[i]                           , "L");

    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    leg.Draw("same")

    canvas.Print("SW_"+R+"_"+I +".png")
    del canvas
    gc.collect()

#variables=["All","l1Pt","l2Pt","l1Eta","l2Eta", "Mll","Drll","Dphill","Njet","jet1pt","jet1eta", "MET"]
#variablesName=["All","leading lepton pt","sub-leading lepton pt","leading lepton eta","sub-leading lepton eta", "M(ll)","Dr(ll)","Dphi(ll)","Njet","leading jet pt","leading jet eta", "MET"]
##HistAddress = '/afs/crc.nd.edu/user/r/rgoldouz/BNV/analysis/hists/'

HistAddress = '/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/hists/'

#couplingsName = [
#'TDUE',
##'TDCE',
#'TSUE',
#'TSCE',
##'TBUE',
##'TBCE'
#]
#
#processname = ['BNV_ST_', 'BNV_TT_']
#
#colors =  [ROOT.kBlack,ROOT.kYellow,ROOT.kGreen,ROOT.kBlue-3,ROOT.kRed-4,ROOT.kOrange-3, ROOT.kOrange-6, ROOT.kCyan-6,ROOT.kViolet, ROOT.kBlue-9, ROOT.kYellow-2]
#
#Hists = []
#Files=[]
#for f in range(len(couplingsName)):
#    l0=[]
#    for i in range(len(processname)):
#        l1=[]
#        Files = ROOT.TFile.Open(HistAddress +processname[i] + couplingsName[f] + '.root')
#        for numvar, namevar in enumerate(variables):
#            h= Files.Get(namevar)
#            h.SetFillColor(0)
#            h.SetLineColor(colors[f+i])
#            l1.append(h)
#        l0.append(l1)
#    Hists.append(l0)       
variables=["lep1Pt","lep1Eta","lep1Phi","lep2Pt","lep2Eta","lep2Phi","llM","llPt","llDr","llDphi","jet1Pt","jet1Eta","jet1Phi","njet","nbjet","Met","MetPhi","nVtx","llMZw", "topMass","topL1Dphi","topL1Dr","topL1DptOsumPt","topPt"]
FR = ROOT.TFile.Open(HistAddress + '2017_BNV_TT_TBCE.root')
for H in variables:
    HEFT = FR.Get('emu_llB1_'+H)
    text = 'EFTrwgt1_cS_1_cT_10'
    wc1 = ROOT.WCPoint(text)
    print str(HEFT.GetBinContent(0,wc1))
    HEFT.Scale(wc1)
    canvas = ROOT.TCanvas('a','a',50,50,865,780)
    HEFT.Draw()
    canvas.Print(H+"BNV.png")


#FR = ROOT.TFile.Open(HistAddress + 'ANoutput.root')
#print HistAddress + 'ANoutput.root'
#HEFT = FR.Get("crossSection")
#HEFT.GetSumFit().save('Coup.tex')
#text = 'EFTrwgt1_cS_1_cT_1'
#wc1 = ROOT.WCPoint(text)
#HEFT.Scale(wc1)
#print str(HEFT.Integral())
#canvas = ROOT.TCanvas('a','a',50,50,865,780)
#HEFT.Draw()
#canvas.Print("BNV.png")




#canvas = ROOT.TCanvas('a','a',50,50,865,780)
#A=FR.Get("l1Pt")
##A.Draw()
#B=A.Clone()
#wc1=ROOT.WCPoint("smpt")
#B.Scale(wc1)
#B.SetLineColor(2)
##B.Draw("same")
#B.ScaleFits(1/B.Integral())
#C=A.Clone()
#wc2=ROOT.WCPoint("EFTrwgt1_ctlT2_10_ctlT1_10")
#C.Scale(wc2)
#C.ScaleFits(1/C.Integral())
#C.SetLineColor(1)
#C.Divide(B)
#C.SetMinimum(0);
#C.SetMaximum(2);
#C.Draw("same")
#D=A.Clone()
#wc3=ROOT.WCPoint("EFTrwgt1_cQl32_10_cQl31_10")
#D.Scale(wc3)
#D.SetLineColor(6)
##D.Draw("same")
#canvas.Print("W.png")
#
#for numvar, namevar in enumerate(variables):
#    HH=[]
#    HHName=[]
#    for f in range(len(couplingsName)):
#        for i in range(len(processname)):
#            h2 = Hists[f][i][numvar].Clone()
##            h1 = FR.Get(namevar).Clone()
##            h1.GetSumFit().save("FIT.txt"+namevar)
#            text = 'EFTrwgt1_cS_1_cT_1' 
##            print couplingsName[f] + ' , ' + str(h1.GetSumFit().evalPoint(ROOT.WCPoint(text))/41530.0)
#            wc1 = ROOT.WCPoint(text)
##            h2.Scale(wc1)
##            h2.ScaleFits(1/h2.Integral())
##            print str(h2.Integral())
#            HH.append(EFTtoNormal(h2,wc1))
#            HHName.append(processname[i]+couplingsName[f])
##            HH.append(h2)
##            HHName.append('nominal ' + couplingsName[f] + '=1')
#    drawHistallW(HH,HHName, namevar , "","Normalized entries")
#
#for numvar, namevar in enumerate(variables):
#    for f in range(len(couplingsName)):
#        for i in range(len(processname)):
#            HH=[]
#            HHName=[]
#            h2 = Hists[f][i][numvar].Clone()
#            text = 'EFTrwgt1_cS_1_cT_1'
#            wc1 = ROOT.WCPoint(text)
#            h2.Scale(wc1)
#            HH.append(h2)
#            HHName.append(processname[i]+couplingsName[f]+"(S=1,T=1)")
#            h2 = Hists[f][i][numvar].Clone()
#            text = 'EFTrwgt1_cS_1_cT_0'
#            wc1 = ROOT.WCPoint(text)
#            h2.Scale(wc1)
#            HH.append(h2)
#            HHName.append(processname[i]+couplingsName[f]+"(S=1,T=0)")
#            h2 = Hists[f][i][numvar].Clone()
#            text = 'EFTrwgt1_cS_0_cT_1'
#            wc1 = ROOT.WCPoint(text)
#            h2.Scale(wc1)
#            HH.append(h2)
#            HHName.append(processname[i]+couplingsName[f]+"(S=0,T=1)")
#            drawHistallW(HH,HHName, namevar , processname[i]+couplingsName[f],"N events (2017)")

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

def compare3Hist(A, B, C, textA="A", textB="B", textC="C",label_name="sample", can_name="can"):

    canvas = ROOT.TCanvas(can_name,can_name,10,10,1100,628)
    canvas.SetRightMargin(0.15)
    canvas.cd()

    pad_name = "pad"
    pad1=ROOT.TPad(pad_name, pad_name, 0.05, 0.3, 1, 0.99 , 0)
    pad1.Draw()
    pad1.SetLogy()
    pad2=ROOT.TPad(pad_name, pad_name, 0.05, 0.05, 1, 0.3 , 0)
    pad2.SetGridy();
    pad2.Draw()
    pad1.cd()

    A.SetLineColor( 1 )
    B.SetLineColor( 2 )
    C.SetLineColor( 4 )

    A.SetTitle("")
    A.GetXaxis().SetTitle('BDT output')
    A.GetYaxis().SetTitle('Event ')
    A.GetXaxis().SetTitleSize(0.05)
    A.GetYaxis().SetTitleSize(0.05)
    A.SetMaximum(1.2*max(A.GetMaximum(),B.GetMaximum(),C.GetMaximum()));
    A.SetMinimum(0.1);
    B.SetMinimum(0.1);
    C.SetMinimum(0.1);
    A.GetYaxis().SetTitleOffset(0.7)
    A.Draw()
    B.Draw('esame')
    C.Draw('esame')

    legend = ROOT.TLegend(0.7,0.75,1,1)
    legend.AddEntry(A ,textA,'l')
    legend.AddEntry(B ,textB,'l')
    legend.AddEntry(C ,textC,'l')
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.05)
    legend.Draw("same")

    Label_channel2 = ROOT.TLatex(0.15,0.65,'#color[2]{'+can_name.split('_')[-1]+'}')
    Label_channel2.SetNDC()
    Label_channel2.SetTextFont(42)
    Label_channel2.SetTextSize(0.085)
    Label_channel2.Draw("same")

    pad2.cd()
    ratioB = A.Clone()
    ratioB.Divide(B)
    ratioB.SetLineColor( 2 )
#    r = ratioB.Clone()
    fontScale = 2
    nbin = ratioB.GetNbinsX()
    x_min= ratioB.GetBinLowEdge(1)
    x_max= ratioB.GetBinLowEdge(nbin)+ratioB.GetBinWidth(nbin)
#    ratio_y_min=0.95*r.GetBinContent(r.FindFirstBinAbove(0))
#    ratio_y_max=1.05*r.GetBinContent(r.GetMaximumBin())
    dummy_ratio = ROOT.TH2D("dummy_ratio","",nbin,x_min,x_max,1,0.8,1.2)
    dummy_ratio.SetStats(ROOT.kFALSE)
    dummy_ratio.GetYaxis().SetTitle('Ratio')
    dummy_ratio.GetXaxis().SetTitle("")
    dummy_ratio.GetXaxis().SetTitleSize(0.05*fontScale)
    dummy_ratio.GetXaxis().SetLabelSize(0.05*fontScale)
    dummy_ratio.GetXaxis().SetMoreLogLabels()
    dummy_ratio.GetXaxis().SetNoExponent()
    dummy_ratio.GetYaxis().SetNdivisions(505)
    dummy_ratio.GetYaxis().SetTitleSize(0.07*fontScale)
    dummy_ratio.GetYaxis().SetLabelSize(0.05 *fontScale)
    dummy_ratio.GetYaxis().SetTitleOffset(0.3)
    dummy_ratio.Draw()
#    ratioB.Draw("esame")

    ratioC = C.Clone()
    ratioC.Divide(B)
    ratioC.SetLineColor( 4 )
    ratioC.Draw("esame")

    Label_channel = ROOT.TLatex(0.15,0.8,textB + '/'+textC)
    Label_channel.SetNDC()
    Label_channel.SetTextFont(42)
    Label_channel.Draw("same")

    canvas.Print("3H_" + can_name + ".png")
    del canvas
    gc.collect()


def puDraw(hists,Fnames,year,can):
    ratio=[]
    canvas = ROOT.TCanvas(can,can,50,50,865,780)
    canvas.SetGrid();
    canvas.SetBottomMargin(0.17)
    canvas.cd()

    legend = ROOT.TLegend(0.7,0.65,0.9,0.88)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.03)

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
    y_max=1.2*hists[0].GetMaximum()
    hists[0].SetTitle("")
    hists[0].GetYaxis().SetTitle('Probability')
    hists[0].GetXaxis().SetLabelSize(0)
    hists[0].GetYaxis().SetTitleOffset(0.8)
    hists[0].GetYaxis().SetTitleSize(0.07)
    hists[0].GetYaxis().SetLabelSize(0.04)
    hists[0].GetYaxis().SetRangeUser(y_min,y_max)
    hists[0].Draw("hist")
    nomhist = hists[0].Clone()
    for H in range(0,len(hists)):
        hists[H].SetLineWidth(2)
        if 'nominal' in Fnames[H]:
            hists[H].SetLineColor(1)
            legend.AddEntry(hists[H],'Data','L')
            nomhist = hists[H].Clone()
        if 'up' in Fnames[H]:
            hists[H].SetLineColor(16)
            hists[H].SetLineStyle(9)
            legend.AddEntry(hists[H],'Data (+1 #sigma)','L')
            ratio.append(hists[H].Clone())
        if 'down' in Fnames[H]:
            hists[H].SetLineColor(16)
            hists[H].SetLineStyle(6)
            legend.AddEntry(hists[H],'Data (-1 #sigma)','L')
            ratio.append(hists[H].Clone())
        if 'MC' in Fnames[H]:
            hists[H].SetLineColor(2)
            legend.AddEntry(hists[H],'MC','L')
            ratio.append(hists[H].Clone())
        hists[H].Draw("histSAME")
    hists[0].Draw("AXISSAMEY+")
    hists[0].Draw("AXISSAMEX+")

    legend.Draw("same")
    Lumi = '137.42'
    
    if (year == 'UL2016preVFP'):
        Lumi = '19.52'
    if (year == 'UL2016postVFP'):
        Lumi = '16.81'
    if (year == 'UL2017'):
        Lumi = '41.48'
    if (year == 'UL2018'):
        Lumi = '59.83'
    label_cms="CMS Preliminary - " + year
    Label_cms = ROOT.TLatex(0.2,0.92,label_cms)
    Label_cms.SetNDC()
    Label_cms.SetTextFont(61)
    Label_cms.Draw()
    Label_lumi = ROOT.TLatex(0.71,0.92,Lumi+" fb^{-1} (13 TeV)")
    Label_lumi.SetNDC()
    Label_lumi.SetTextFont(42)
    Label_lumi.Draw("same")
    pad1.Update()

    pad2.cd()
    print ratio
    for rH in ratio: 
        rH.SetTitle("")
        rH.SetMarkerStyle(20)
        rH.SetMarkerSize(1.2)
        rH.GetXaxis().SetTitle('True Number of Interactions')
        rH.GetYaxis().CenterTitle()
        rH.GetXaxis().SetMoreLogLabels()
        rH.GetXaxis().SetNoExponent()
        rH.GetXaxis().SetTitleSize(0.04/0.3)
        rH.GetYaxis().SetTitleSize(0.04/0.3)
        rH.GetXaxis().SetTitleFont(42)
        rH.GetYaxis().SetTitleFont(42)
        rH.GetXaxis().SetTickLength(0.05)
        rH.GetYaxis().SetTickLength(0.05)
        rH.GetXaxis().SetLabelSize(0.115)
        rH.GetYaxis().SetLabelSize(0.089)
        rH.GetXaxis().SetLabelOffset(0.02)
        rH.GetYaxis().SetLabelOffset(0.01)
        rH.GetYaxis().SetTitleOffset(0.42)
        rH.GetXaxis().SetTitleOffset(1.1)
        rH.GetYaxis().SetNdivisions(504)
        rH.GetYaxis().SetRangeUser(0.5,2)
        rH.Divide(nomhist)
        for b in range(rH.GetNbinsX()):
            if rH.GetBinContent(b+1)>0:
                rH.SetBinContent(b+1,1/rH.GetBinContent(b+1))
        rH.SetStats(ROOT.kFALSE)
        rH.GetYaxis().SetTitle('Data/x')
        rH.Draw('histSAME')
        rH.Draw("AXISSAMEY+")
        rH.Draw("AXISSAMEX+")


    canvas.Print(can + ".png")
    del canvas
    gc.collect()



samples = {}
#/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/
samples["UL2017_down"]=[    "PileupHistogram-goldenJSON-13tev-2017-66000ub-99bins.root"]
samples["UL2017_nominal"]=[    "PileupHistogram-goldenJSON-13tev-2017-69200ub-99bins.root"]
samples["UL2017_up"]=[    "PileupHistogram-goldenJSON-13tev-2017-72400ub-99bins.root"]

samples["UL2016preVFP_down"]=[    "PileupHistogram-goldenJSON-13tev-2016-preVFP-66000ub-99bins.root"]
samples["UL2016preVFP_nominal"]=[    "PileupHistogram-goldenJSON-13tev-2016-preVFP-69200ub-99bins.root"]
samples["UL2016preVFP_up"]=[    "PileupHistogram-goldenJSON-13tev-2016-preVFP-72400ub-99bins.root"]

samples["UL2016postVFP_down"]=[    "PileupHistogram-goldenJSON-13tev-2016-postVFP-66000ub-99bins.root"]
samples["UL2016postVFP_nominal"]=[    "PileupHistogram-goldenJSON-13tev-2016-postVFP-69200ub-99bins.root"]
samples["UL2016postVFP_up"]=[    "PileupHistogram-goldenJSON-13tev-2016-postVFP-72400ub-99bins.root"]

samples["UL2018_down"]=[    "PileupHistogram-goldenJSON-13tev-2018-66000ub-99bins.root"]
samples["UL2018_nominal"]=[    "PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root"]
samples["UL2018_up"]=[    "PileupHistogram-goldenJSON-13tev-2018-72400ub-99bins.root"]

MC_bins = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 
27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 
39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 
51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 
63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 
75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 
87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98
]

#https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2016_25ns_UltraLegacy_PoissonOOTPU_cfi.py
MC2016_value=[
    1.00402360149e-05, 5.76498797172e-05, 7.37891400294e-05, 0.000110932895295, 0.000158857714773,
    0.000368637432599, 0.000893114107873, 0.00189700774575, 0.00358880167437, 0.00636052573486,
    0.0104173961179, 0.0158122597405, 0.0223785660712, 0.0299186888073, 0.0380275944896,
    0.0454313901624, 0.0511181088317, 0.0547434577348, 0.0567906239028, 0.0577145461461,
    0.0578176902735, 0.0571251566494, 0.0555456541498, 0.053134383488, 0.0501519041462,
    0.0466815838899, 0.0429244592524, 0.0389566776898, 0.0348507152776, 0.0307356862528,
    0.0267712092206, 0.0229720184534, 0.0193388653099, 0.0159602510813, 0.0129310510552,
    0.0102888654183, 0.00798782770975, 0.00606651703058, 0.00447820948367, 0.00321589786478,
    0.0022450422045, 0.00151447388514, 0.000981183695515, 0.000609670479759, 0.000362193408119,
    0.000211572646801, 0.000119152364744, 6.49133515399e-05, 3.57795801581e-05, 1.99043569043e-05,
    1.13639319832e-05, 6.49624103579e-06, 3.96626216416e-06, 2.37910222874e-06, 1.50997403362e-06,
    1.09816650247e-06, 7.31298519122e-07, 6.10398791529e-07, 3.74845774388e-07, 2.65177281359e-07,
    2.01923536742e-07, 1.39347583555e-07, 8.32600052913e-08, 6.04932421298e-08, 6.52536630583e-08,
    5.90574603808e-08, 2.29162474068e-08, 1.97294602668e-08, 1.7731096903e-08, 3.57547932012e-09,
    1.35039815662e-09, 8.50071242076e-09, 5.0279187473e-09, 4.93736669066e-10, 8.13919708923e-10,
    5.62778926097e-09, 5.15140589469e-10, 8.21676746568e-10, 0.0, 1.49166873577e-09,
    8.43517992503e-09, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0]


#https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2017_25ns_UltraLegacy_PoissonOOTPU_cfi.py
MC2017_value =[
    1.1840841518e-05, 3.46661037703e-05, 8.98772521472e-05, 7.47400487733e-05, 0.000123005176624,
    0.000156501700614, 0.000154660478659, 0.000177496185603, 0.000324149805611, 0.000737524009713,
    0.00140432980253, 0.00244424508696, 0.00380027898037, 0.00541093042612, 0.00768803501793,
    0.010828224552, 0.0146608623707, 0.01887739113, 0.0228418813823, 0.0264817796874,
    0.0294637401336, 0.0317960986171, 0.0336645950831, 0.0352638818387, 0.036869429333,
    0.0382797316998, 0.039386705577, 0.0398389681346, 0.039646211131, 0.0388392805703,
    0.0374195678161, 0.0355377892706, 0.0333383902828, 0.0308286549265, 0.0282914440969,
    0.0257860718304, 0.02341635055, 0.0213126338243, 0.0195035612803, 0.0181079838989,
    0.0171991315458, 0.0166377598339, 0.0166445341361, 0.0171943735369, 0.0181980997278,
    0.0191339792146, 0.0198518804356, 0.0199714909193, 0.0194616474094, 0.0178626975229,
    0.0153296785464, 0.0126789254325, 0.0100766041988, 0.00773867100481, 0.00592386091874,
    0.00434706240169, 0.00310217013427, 0.00213213401899, 0.0013996000761, 0.000879148859271,
    0.000540866009427, 0.000326115560156, 0.000193965828516, 0.000114607606623, 6.74262828734e-05,
    3.97805301078e-05, 2.19948704638e-05, 9.72007976207e-06, 4.26179259146e-06, 2.80015581327e-06,
    1.14675436465e-06, 2.52452411995e-07, 9.08394910044e-08, 1.14291987912e-08, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0]
#https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2018_25ns_UltraLegacy_PoissonOOTPU_cfi.py
MC2018_value =[
    8.89374611122e-07, 1.1777062868e-05, 3.99725585118e-05, 0.000129888015252, 0.000265224848687,
    0.000313088635109, 0.000353781668514, 0.000508787237162, 0.000873670065767, 0.00147166880932,
    0.00228230649018, 0.00330375581273, 0.00466047608406, 0.00624959203029, 0.00810375867901,
    0.010306521821, 0.0129512453978, 0.0160303925502, 0.0192913204592, 0.0223108613632,
    0.0249798930986, 0.0273973789867, 0.0294402350483, 0.031029854302, 0.0324583524255,
    0.0338264469857, 0.0351267479019, 0.0360320204259, 0.0367489568401, 0.0374133183052,
    0.0380352633799, 0.0386200967002, 0.039124376968, 0.0394201612616, 0.0394673457109,
    0.0391705388069, 0.0384758587461, 0.0372984548399, 0.0356497876549, 0.0334655175178,
    0.030823567063, 0.0278340752408, 0.0246009685048, 0.0212676009273, 0.0180250593982,
    0.0149129830776, 0.0120582333486, 0.00953400069415, 0.00738546929512, 0.00563442079939,
    0.00422052915668, 0.00312446316347, 0.00228717533955, 0.00164064894334, 0.00118425084792,
    0.000847785826565, 0.000603466454784, 0.000419347268964, 0.000291768785963, 0.000199761337863,
    0.000136624574661, 9.46855200945e-05, 6.80243180179e-05, 4.94806013765e-05, 3.53122628249e-05,
    2.556765786e-05, 1.75845711623e-05, 1.23828210848e-05, 9.31669724108e-06, 6.0713272037e-06,
    3.95387384933e-06, 2.02760874107e-06, 1.22535149516e-06, 9.79612472109e-07, 7.61730246474e-07,
    4.2748847738e-07, 2.41170461205e-07, 1.38701083552e-07, 3.37678010922e-08, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0]


gr2016preVFP=[]
gr2016postVFP=[]
gr2017=[]
gr2018=[]
Ngr2016preVFP=[]
Ngr2016postVFP=[]
Ngr2017=[]
Ngr2018=[]

f = ROOT.TFile.Open("2018_DY.root")
h_2018_DY = f.Get('truePU')
h_2018_DY.Scale(1/h_2018_DY.Integral())

f = ROOT.TFile.Open("2017_DY.root")
h_2017_DY = f.Get('truePU')
h_2017_DY.Scale(1/h_2017_DY.Integral())

f = ROOT.TFile.Open("2016preVFP_DY.root")
h_2016preVFP_DY = f.Get('truePU')
h_2016preVFP_DY.Scale(1/h_2016preVFP_DY.Integral())

f = ROOT.TFile.Open("2016postVFP_DY.root")
h_2016postVFP_DY = f.Get('truePU')
h_2016postVFP_DY.Scale(1/h_2016postVFP_DY.Integral())

weights={}
for key, value in samples.items():
    if '2016preVFP' in key:
        W = []
        f = ROOT.TFile.Open(value[0])
        h = f.Get('pileup')
        h.Scale(1/h.Integral()) 
        gr2016preVFP.append(h)
        Ngr2016preVFP.append(key)
        for b in range(len(MC_bins)):
            if h_2016postVFP_DY.GetBinContent(b+1)>0:
#                W.append(h.GetBinContent(b+1)/MC2016_value[b])
                W.append(h.GetBinContent(b+1)/h_2016postVFP_DY.GetBinContent(b+1))
        weights[key] = W
    if '2016postVFP' in key:
        W = []
        f = ROOT.TFile.Open(value[0])
        h = f.Get('pileup')
        h.Scale(1/h.Integral())
        gr2016postVFP.append(h)
        Ngr2016postVFP.append(key)
        for b in range(len(MC_bins)):
            if h_2016postVFP_DY.GetBinContent(b+1)>0:
#                W.append(h.GetBinContent(b+1)/MC2016_value[b])
                W.append(h.GetBinContent(b+1)/h_2016postVFP_DY.GetBinContent(b+1))
        weights[key] = W
    if '2017' in key:
        W = []
        f = ROOT.TFile.Open(value[0])
        h = f.Get('pileup')
        h.Scale(1/h.Integral())
        gr2017.append(h)
        Ngr2017.append(key)
        for b in range(len(MC_bins)):
            if h_2017_DY.GetBinContent(b+1)>0:
#                W.append(h.GetBinContent(b+1)/MC2017_value[b])
                W.append(h.GetBinContent(b+1)/h_2017_DY.GetBinContent(b+1))
        weights[key] = W
    if '2018' in key:
        W = []
        f = ROOT.TFile.Open(value[0])
        h = f.Get('pileup')
        h.Scale(1/h.Integral())
        gr2018.append(h)
        Ngr2018.append(key)
        for b in range(len(MC_bins)):
            if h_2018_DY.GetBinContent(b+1)>0:
#                W.append(h.GetBinContent(b+1)/MC2018_value[b])
                W.append(h.GetBinContent(b+1)/h_2018_DY.GetBinContent(b+1))
        weights[key] = W

hcopy = gr2016preVFP[0].Clone()
for b in range(hcopy.GetNbinsX()):
    hcopy.SetBinContent(b+1,MC2016_value[b])
gr2016preVFP.append(hcopy)
Ngr2016preVFP.append('MC_2016preVFP')

hcopy = gr2016postVFP[0].Clone()
for b in range(hcopy.GetNbinsX()):
    hcopy.SetBinContent(b+1,MC2016_value[b])
gr2016postVFP.append(hcopy)
Ngr2016postVFP.append('MC_2016preVFP')

hcopy = gr2017[0].Clone()
for b in range(hcopy.GetNbinsX()):
    hcopy.SetBinContent(b+1,MC2017_value[b])
gr2017.append(hcopy)
Ngr2017.append('MC_2017')

hcopy = gr2018[0].Clone()
for b in range(hcopy.GetNbinsX()):
    hcopy.SetBinContent(b+1,MC2018_value[b])
gr2018.append(hcopy)
Ngr2018.append('MC_2018')

for key, value in samples.items():
    print 'double pu' + key + '['+str(len(weights[key])) +']={'
    print weights[key]
print sum(MC2016_value)
print sum(MC2016_value)
print sum(MC2017_value)
print sum(MC2018_value)

puDraw(gr2016preVFP,Ngr2016preVFP,'UL2016preVFP','pu2016preVFP')
puDraw(gr2016postVFP,Ngr2016postVFP,'UL2016postVFP','pu2016postVFP')
puDraw(gr2017,Ngr2017,'UL2017','pu2017')
puDraw(gr2018,Ngr2018,'UL2018','pu2018')

#f = ROOT.TFile.Open("2018_DY.root")
#h = f.Get('truePU')
#h.Scale(1/h.Integral())
#hcopy = gr2018[0].Clone()
#for b in range(hcopy.GetNbinsX()):
#    if h.GetBinContent(b+1)>0:
#        hcopy.SetBinContent(b+1,h.GetBinContent(b+1))
#    else:
#        hcopy.SetBinContent(b+1,00000001)
#MCpu2018=[hcopy]
#MCpu2018Name=["DYMC"]
#for i in range(len(gr2018)):
#    if 'up' not in Ngr2018[i] and 'down' not in Ngr2018[i]:
#        MCpu2018.append(gr2018[i])
#        MCpu2018Name.append(Ngr2018[i])
#print MCpu2018Name
#compare3Hist(MCpu2018[1],MCpu2018[0],MCpu2018[2],MCpu2018Name[1],MCpu2018Name[0],MCpu2018Name[2], 'PU2018' ,'pu2018')
#
#f = ROOT.TFile.Open("2017_DY.root")
#h = f.Get('truePU')
#h.Scale(1/h.Integral())
#hcopy = gr2017[0].Clone()
#for b in range(hcopy.GetNbinsX()):
#    if h.GetBinContent(b+1)>0:
#        hcopy.SetBinContent(b+1,h.GetBinContent(b+1))
#    else:
#        hcopy.SetBinContent(b+1,00000001)
#MCpu2017=[hcopy]
#MCpu2017Name=["DYMC"]
#for i in range(len(gr2017)):
#    if 'up' not in Ngr2017[i] and 'down' not in Ngr2017[i]:
#        MCpu2017.append(gr2017[i])
#        MCpu2017Name.append(Ngr2017[i])
#compare3Hist(MCpu2017[1],MCpu2017[0],MCpu2017[2],MCpu2017Name[1],MCpu2017Name[0],MCpu2017Name[2], 'PU2017' ,'pu2017')
#
#f = ROOT.TFile.Open("2016postVFP_DY.root")
#h = f.Get('truePU')
#h.Scale(1/h.Integral())
#hcopy = gr2016postVFP[0].Clone()
#for b in range(hcopy.GetNbinsX()):
#    if h.GetBinContent(b+1)>0:
#        hcopy.SetBinContent(b+1,h.GetBinContent(b+1))
#    else:
#        hcopy.SetBinContent(b+1,00000001)
#MCpu2016postVFP=[hcopy]
#MCpu2016postVFPName=["DYMC"]
#for i in range(len(gr2016postVFP)):
#    if 'up' not in Ngr2016postVFP[i] and 'down' not in Ngr2016postVFP[i]:
#        MCpu2016postVFP.append(gr2016postVFP[i])
#        MCpu2016postVFPName.append(Ngr2016postVFP[i])
#compare3Hist(MCpu2016postVFP[1],MCpu2016postVFP[0],MCpu2016postVFP[2],MCpu2016postVFPName[1],MCpu2016postVFPName[0],MCpu2016postVFPName[2], 'PU2016postVFP' ,'pu2016postVFP')
#
#f = ROOT.TFile.Open("2016preVFP_DY.root")
#h = f.Get('truePU')
#h.Scale(1/h.Integral())
#hcopy = gr2016preVFP[0].Clone()
#for b in range(hcopy.GetNbinsX()):
#    if h.GetBinContent(b+1)>0:
#        hcopy.SetBinContent(b+1,h.GetBinContent(b+1))
#    else:
#        hcopy.SetBinContent(b+1,00000001)
#MCpu2016preVFP=[hcopy]
#MCpu2016preVFPName=["DYMC"]
#for i in range(len(gr2016preVFP)):
#    if 'up' not in Ngr2016preVFP[i] and 'down' not in Ngr2016preVFP[i]:
#        MCpu2016preVFP.append(gr2016preVFP[i])
#        MCpu2016preVFPName.append(Ngr2016preVFP[i])
#compare3Hist(MCpu2016preVFP[1],MCpu2016preVFP[0],MCpu2016preVFP[2],MCpu2016preVFPName[1],MCpu2016preVFPName[0],MCpu2016preVFPName[2], 'PU2016preVFP' ,'pu2016preVFP')



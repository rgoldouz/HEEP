import math
import gc
import os,sys
import sys
import ROOT
import numpy as npi
import copy
from array import array
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1;")
ROOT.TH1.AddDirectory(ROOT.kFALSE)
ROOT.gStyle.SetOptStat(0)


import gc

################################## MY SIGNAL AND SM BG ################################
def flatten(hist,x1,x2):
    hist.Fit("pol0","","",x1,x2)
    func = hist.GetFunction("pol0")
    norm = func.GetParameter(0)
    new_hist = ROOT.TH1F("fit","fit",1,x1,x2)
    new_hist.SetBinContent(1, norm)
    return new_hist

def compare6Hist(A,T,label_name="sample", can_name="can"):
    MarkerStyle=[25,20,21,22,23,4,24,26]
    MarkerColor=[6,3,4,2,5,7,8,9]
    a,b,c,d="h_pt_Barrel_totalEt25".split("_")

    if len(label_name.split("_"))==4:
        a,b,c,d = label_name.split("_")
    if len(label_name.split("_"))==3:
        a,b,d = label_name.split("_")
        c=""
    if 'pt' in b:
        b='pt'
    canvas = ROOT.TCanvas(can_name,can_name,10,10,1100,628)
    canvas.cd()
    canvas.SetGrid()
    pad_name = "pad"
    pad1=ROOT.TPad(pad_name, pad_name, 0, 0.02, 1, 1 , 0)
    pad1.SetGrid()
    pad1.Draw()
    pad1.cd()
    leg = ROOT.TLegend(0.6,0.8,0.9,0.94)
#    pad1.SetLogy()
    A[0].SetTitle("")
    A[0].GetXaxis().SetTitle(b)
    A[0].GetYaxis().SetTitle('Efficiency')
    A[0].SetTitle(d)
    A[0].GetXaxis().SetTitleSize(0.05)
    A[0].GetYaxis().SetTitleSize(0.05)
    A[0].GetYaxis().SetTitleOffset(0.7)
    A[0].GetXaxis().SetTitleOffset(0.85)
    if 'Barrel' in label_name+can_name:
        A[0].SetTitle("Barrel")
    elif 'Endcap' in label_name+can_name:
        A[0].SetTitle("Endcaps")
    else:
        A[0].SetTitle(d)
    for i in range(len(A)):
        A[i].SetMarkerStyle(MarkerStyle[i])
        A[i].SetMarkerColor(MarkerColor[i])
        A[i].SetLineColor(MarkerColor[i])
        A[i].SetMarkerSize(1);
        A[i].GetYaxis().SetRangeUser(0.8*A[0].GetMinimum(),1.05*A[0].GetMaximum())
        if 'pt' in b:
            A[0].GetXaxis().SetTitle("E_{T}^{SC} (GeV)")
        if 'pt' in b and 'Et' in d:
            A[i].GetXaxis().SetRangeUser(35,500)
            A[i].GetYaxis().SetRangeUser(0.8,1.05)
        if 'Et' not in d:
            A[i].GetXaxis().SetRangeUser(35,400)
            A[i].GetYaxis().SetRangeUser(0.8,1.05)
        A[i].Draw("e HIST SAME")
        leg.AddEntry(A[i], T[i]                           , "p");
    A[0].Draw("AXISSAMEY+")
    A[0].Draw("AXISSAMEX+")
#    Fitf = ROOT.TF1( "pol0" ,"pol0" ,A[0].GetBinLowEdge(1), A[0].GetBinLowEdge(A[0].GetNbinsX() + 1))
#    A[len(A)-1].Fit(Fitf,"E")
#    A[len(A)-1].Draw("SAME E")
#    if A[len(A)-1].Integral()>0:
#        fitHist=flatten(A[len(A)-1],35, A[0].GetBinLowEdge(A[0].GetNbinsX() + 1))
#        fitHist.SetLineStyle(5)
#        fitHist.SetLineColor(1)
#        fitHist.Draw("SAME")
#        leg.AddEntry(fitHist, 'Fit '+T[len(A)-1]                           , "L");   

    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.Draw("same")

    label = ROOT.TLatex()
    label.SetTextAlign(12)
    label.SetTextFont(42)
    label.SetTextSize(0.06)
    label.SetNDC(ROOT.kTRUE)
    label.DrawLatex(0.1,0.95, can_name)

    if not os.path.exists("plots_6D_stubeff"):
        os.mkdir( "plots_6D_stubeff", 0755 );

    canvas.Print("plots_6D_stubeff/" + "6D_" + can_name +"_"+ label_name + ".png")
    del canvas
    gc.collect()

def draw1dHist(A,textA="A", label_name="sample", can_name="can"):
    a,b,c,d="h_pt_Barrel_totalEt25".split("_")

    if len(label_name.split("_"))==4:
        a,b,c,d = label_name.split("_")
    if len(label_name.split("_"))==3:
        a,b,d = label_name.split("_")
        c=""

    canvas = ROOT.TCanvas(can_name,can_name,10,10,1100,628)
    canvas.cd()

    pad_name = "pad"
    pad1=ROOT.TPad(pad_name, pad_name, 0.05, 0.05, 1, 0.99 , 0)
    pad1.Draw()
    pad1.SetGridy()
    pad1.SetGridx()
    pad1.cd()

    A.GetXaxis().SetTitle(b)
    A.GetYaxis().SetTitle('Efficiency')
    A.SetLineColor( 1 )
    A.SetLineWidth( 2 )
    A.SetTitle(c + ',' + d)
    A.GetXaxis().SetTitleSize(0.05)
    A.GetYaxis().SetTitleSize(0.05)
    A.GetYaxis().SetTitleOffset(0.7)
    A.GetYaxis().SetRangeUser(0.8*A.GetMinimum(),1.2*A.GetMaximum())
    if 'pt' in b:
        b='pt'
    if 'pt' in b and 'Et' in d:
        A.GetXaxis().SetRangeUser(20,500)
    A.SetMarkerColor(2);
    A.SetMarkerSize(1.5);
    A.SetMarkerStyle(22);
    A.Draw("p HIST SAME")
    legend = ROOT.TLegend(0.6,0.9,1,1)
    legend.AddEntry(A ,textA,'p')
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
#    legend.Draw("same")
    canvas.Print("1D_" + can_name +"_"+ label_name + ".png")
    del canvas
    gc.collect()


samples = [
'2016preVFP_data.root',
'2016postVFP_data.root',
'2017_data.root',
'2018_data.root',
#'data_2017_Bv3_SingleElectron.root',
#'data_2017_B_DoubleEG.root'
#'2018_DY50.root',
#'data_2017_B_SingleElectron.root',
#'data_2017_Bv2_SingleElectron.root'
]

samplename = [
'2016preVFP',
'2016postVFP',
'2017',
'2018',
#'2017_BSigleEle',
#'2017_BDoubleEG',
#'2018_Bv2',
#'2018_DY50'
]

effType=['Et','WPTightEt35','CaloIdLMWPMS2']
effType=['HE', 'ClusterShape','PixelMatch', 'CaloIdLMWPMS2','HESeeded', 'ClusterShapeSeeded','PixelMatchSeeded', 'CaloIdLMWPMS2Seeded']
region = ['Barrel', 'Endcaps']
variable1 = ["pt", "pteta1", "pteta2", "pteta3"]
variable1 = ["pt"]
variable1Name={
'Barrelpteta1': '|#eta|<0.79',
'Barrelpteta2': '0.79<|#eta|<1.10',
'Barrelpteta3': '1.10<|#eta|<1.44',
'Endcapspteta1': '1.56<|#eta|<1.70',
'Endcapspteta2': '1.70<|#eta|<2.10',
'Endcapspteta3': '2.10<|#eta|<2.50',
}
variable2 = ["eta","mass"]
variable2 =["pt"]
variable3 = ["lep1Pt","lep1Eta","lep1Phi","lep2Pt","lep2Eta","lep2Phi","llM","llPt","llDr","llDphi"]
directory = '/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/hists/'

Hists1 = []
for nums, s in enumerate(samples):
    file1 = ROOT.TFile.Open(directory + s)
    l0=[]
    for numr, r in enumerate(region):
        l1=[]
        for nume, e in enumerate(effType):
            l2=[]
            for numv, v in enumerate(variable1):
                print 'h_'+v+'_'+r+'_'+e
                histA = file1.Get('h_'+v+'_'+r+'_pass'+e)
                histB = file1.Get('h_'+v+'_'+r+'_total'+e)
                histA.Divide(histB)
                l2.append(histA)                
            l1.append(l2)
        l0.append(l1)
    Hists1.append(l0)
    file1.Close()

Hists2 = []
##for nums, s in enumerate(samples):
##    file1 = ROOT.TFile.Open(directory + s)
##    l0=[]
##    for nume, e in enumerate(effType):
##        l1=[]
##        for numv, v in enumerate(variable2):
##            print 'h_'+v+'_'+e
##            histA = file1.Get('h_'+v+'_pass'+e)
##            histB = file1.Get('h_'+v+'_total'+e)
##            histA.Divide(histB)
##            l1.append(histA)
##        l0.append(l1)
##    Hists2.append(l0)
##    file1.Close()

Hists3 = []
for nums, s in enumerate(samples):
    file1 = ROOT.TFile.Open(directory + s)
    l0=[]
    for numv, v in enumerate(variable3):
        print 'h_'+v+'_'+e
        histA = file1.Get('OSee_llOffZ_'+v)
        histB = file1.Get('OSee_ll_'+v)
        print v + ':' +  str(histB.Integral()) +":"+ str(histA.Integral()/histB.Integral())
        histA.Divide(histB)
        l0.append(histA)
    Hists3.append(l0)
    file1.Close()

for nums, s in enumerate(samples):
    for numr, r in enumerate(region):
        for nume, e in enumerate(effType):
            for numv, v in enumerate(variable1):
                canName='h_'+v+'_'+r+'_'+e
                draw1dHist(Hists1[nums][numr][nume][numv],"A", canName+samplename[nums], samplename[nums])


for nums, s in enumerate(samples):
    for nume, e in enumerate(effType):
        for numv, v in enumerate(variable2):
            canName='h_'+v+'_'+e
#            draw1dHist(Hists2[nums][nume][numv],"A", canName, samplename[nums])

for nums, s in enumerate(samples):
    for numr, r in enumerate(region):
        for numv, v in enumerate(variable1):
            Ah=[]
            Th=[]
            AhSeeded=[]
            ThSeeded=[]
            for nume, e in enumerate(effType):
                if 'Seeded' in e:
                    AhSeeded.append(Hists1[nums][numr][nume][numv])
                    ThSeeded.append(e)
                if 'Seeded' not in e:
                    Ah.append(Hists1[nums][numr][nume][numv])
                    Th.append(e)
            if 'postVFP' in s:
                Ah=[Hists1[nums][numr][effType.index("CaloIdLMWPMS2")][numv]]
                Th=["CaloIdLGsfTrkIdVLMWPMS2"]
            if 'preVFP' in s:
                Ah=[Hists1[nums][numr][effType.index("CaloIdLMWPMS2")][numv]]
                Th=["CaloIdL(GsfTrkIdVL/MWPMS2)"]
            if 'preVFP' in s:
                AhSeeded=[Hists1[nums][numr][effType.index("CaloIdLMWPMS2")][numv]]
                ThSeeded=["MWwrtPM"]
            canName='h_'+v+'_'+r+'_'+s
            compare6Hist(Ah,Th,canName, samplename[nums])
            canName='hSeeded_'+v+'_'+r+'_'+s
            compare6Hist(AhSeeded,ThSeeded,canName, samplename[nums])

for nums, s in enumerate(samples):
    for numr, r in enumerate(region):
        for numv, v in enumerate(variable1):
            Ah=[]
            Th=[]
            for nume, e in enumerate(effType):
                if "CaloIdLMWPMS2" not in e:
                    continue
                if 'Seeded' in e:
                    Ah.append(Hists1[nums][numr][nume][numv])
                    Th.append(e)
                if 'Seeded' not in e:
                    Ah.append(Hists1[nums][numr][nume][numv])
                    Th.append(e+'UnSeeded')
            canName='hSUS_'+v+'_'+r+'_'+s
            compare6Hist(Ah,Th,canName, samplename[nums])


for nums, s in enumerate(samples):
    for numr, r in enumerate(region):
        for numv, v in enumerate(variable1):
            Ah=[]
            Th=[]
            H1=Hists1[nums][numr][effType.index("CaloIdLMWPMS2")][numv].Clone()
            H1.Divide(Hists1[nums][numr][effType.index("PixelMatch")][numv])
            Ah.append(H1)
            Th.append('CaloIdLMWPMS2/PixelMatch')
            canName='hRatio_'+v+'_'+r+'_'+s
            compare6Hist(Ah,Th,canName, samplename[nums])

hfile = ROOT.TFile( 'DiEleCaloIdLMWPMS2_HEEPeff.root', 'RECREATE', 'DiEle*CaloIdLMWPMS2 trigger efficiency' )
for nums, s in enumerate(samples):
    for numr, r in enumerate(region):
        Ah=Hists1[nums][numr][effType.index("CaloIdLMWPMS2")][variable1.index("pt")].Clone()
        Ah.SetLineColor(1)
        Ah.SetName('UL'+samplename[nums]+'_'+r+'_Et')
        Ah.Write()
hfile.Write()
hfile.Close()

for nums, s in enumerate(samples):
    for numr, r in enumerate(region):
        for nume, e in enumerate(effType):
            Ah=[]
            Th=[]
            for numv, v in enumerate(variable1):
                if 'pteta'not in v:
                    continue
                canName='h_'+v+'_'+r+'_'+e
                Ah.append(Hists1[nums][numr][nume][numv])
                Th.append(variable1Name[r+v])
#            compare6Hist(Ah,Th,canName, samplename[nums])
            Ah=[]
            Th=[]
for nums, s in enumerate(samples):
    for nume, e in enumerate(effType):
        for numv, v in enumerate(variable1):
            if v!='pt':
                continue
            Ah=[]
            Th=[]
            for numr, r in enumerate(region):
                canName='h_'+v+'_'+r+'_'+e
                Ah.append(Hists1[nums][numr][nume][numv])
                Th.append(r)
#            compare6Hist(Ah,Th,canName, samplename[nums])

for numr, r in enumerate(region):
    for nume, e in enumerate(effType):
        for numv, v in enumerate(variable1):
            if v!='pt':
                continue
            Ah=[]
            Th=[]
            for nums, s in enumerate(samples):
                canName='h_'+v+'_'+r+'_'+e
                Ah.append(Hists1[nums][numr][nume][numv])
                Th.append(samplename[nums])
#            compare6Hist(Ah,Th,canName, 'DM'+samplename[nums])

for numv, v in enumerate(variable3):
    Ah=[]
    Th=[]
    for nums, s in enumerate(samples):
        canName='h_'+v+'_r_dieleToDiph'
        Ah.append(Hists3[nums][numv])
        Th.append(samplename[nums])
#    compare6Hist(Ah,Th,canName, 'phTrigg'+samplename[nums])

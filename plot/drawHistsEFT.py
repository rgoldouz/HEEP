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

bins = array( 'd',[-1,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.55,0.70,0.85,1] )
bins = array( 'd',[-1,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.75,0.9,1] )
binsmll= array( 'd',[0,20,40,60,80,100,140,200,300,400,500,800,1200,2000] )
binsptl = array( 'd',[0,40,80,100,150,200,300,400,500,1000,1500] )
colors =  [ROOT.kBlack,ROOT.kYellow,ROOT.kGreen,ROOT.kBlue-3,ROOT.kRed-4,ROOT.kOrange-3, ROOT.kPink, ROOT.kViolet-1,ROOT.kViolet+8,ROOT.kRed-5,ROOT.kSpring-1,ROOT.kGray+3]

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

def stackPlots(hists, SignalHists, Fnames, ch = "channel", reg = "region", year='2016', var="sample", varname="v"):
    setlog=False
    if 'BDT' in varname:
        setlog=True
    if not os.path.exists(year):
       os.makedirs(year)
    if not os.path.exists(year + '/' + ch):
       os.makedirs(year + '/' + ch)
    if not os.path.exists(year + '/' + ch +'/'+reg):
       os.makedirs(year + '/' + ch +'/'+reg)

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
#    pad1.SetLogy(ROOT.kFALSE)
    pad1.SetLogy(ROOT.kTRUE)
    if setlog:
        pad1.SetLogy(ROOT.kTRUE)

    y_min=0.1
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
    dummy_ratio.GetYaxis().SetRangeUser(0.5,2)
    if 'BDT' in varname:
        dummy_ratio.GetYaxis().SetRangeUser(0.8,1.2)
    dummy_ratio.Divide(SumofMC)
#    dummy_ratio.SetStats(ROOT.kFALSE)
    dummy_ratio.GetYaxis().SetTitle('Data/Pred.')
    dummy_ratio.Draw("AXISSAMEY+")
    dummy_ratio.Draw("AXISSAMEX+")
    dummy_ratio.Draw()
    SumofMC = hists[1].Clone()
    for num in range(2,len(hists)):
        SumofMC.Add(hists[num])
    SumofMC.SetFillColor(13)
    SumofMC.SetLineColor(13)
    SumofMC.SetFillStyle(3244)
    errorRatio=SumofMC.Clone()
    errorRatio.Divide(errorRatio)
    errorRatio.SetFillColor(1)
    errorRatio.SetLineColor(13)
    errorRatio.SetFillStyle(3004)    
    errorRatio.Draw("e2same")
    dummy_ratio.Draw("same")    
    canvas.Print(year + '/' + ch +'/'+reg+'/'+var + ".png")
    del canvas
    gc.collect()


#year=['2016','2017','2018','All']
#year=['2016preVFP', '2016postVFP', '2017','2018','All']
year=['2017']
regions=["ll","llOffZ","llB1", "llBg1"]
#regions=["ll","llOffZ"]
#regions=["ll","llOffZ","llB1", "llBg1","llB1InZ","llB1-BDTm1to0", "llB1-BDT0to0p8", "llB1-BDT0p8to1","llptl1G300", "llB1ptl1G300",   "ll-BB",    "ll-BE",    "ll-EE", "ll-leadMuP", "ll-leadMuN"]
regionsName=["ll","llOffZ","llB1", "llBg1","llB1InZ","llB1-BDTm1to0", "llB1-BDT0to0p8", "llB1-BDT0p8to1","llptl1G300", "llB1ptl1G300",  "ll-BB",    "ll-BE",    "ll-EE", "ll-leadMuP", "ll-leadMuN"]
channels=["ee", "emu", "mumu","ll"];

variables=["lep1Pt","lep1Eta","lep1Phi","lep2Pt","lep2Eta","lep2Phi","llM","llPt","llDr","llDphi","jet1Pt","jet1Eta","jet1Phi","njet","nbjet","Met","MetPhi","nVtx","llMZw", "topMass","topL1Dphi","topL1Dr","topL1DptOsumPt","topPt", "BDT"]
variables=["BDT"]
#"leadElePt","leadEleEta","leadElePhi","leadMuPt",  "leadMuEta", "leadMuPhi","mu1bDr",       "mu2bDr",       "mu1jDrMin",    "mu2jDrMin",    "PtRelMu1jet",  "PtRelMu2jet",  "ele1bDr",      "ele2bDr",      "ele1jDrMin",   "ele2jDrMin",   "PtRelEle1jet", "PtRelEle2jet"]#,"llMZwPt0to100","llMZwPt100to200","llMZwPt200to350","llMZwPt350toInf"]#
variablesName=["p_{T}(leading lepton)","#eta(leading lepton)","#Phi(leading lepton)","p_{T}(sub-leading lepton)","#eta(sub-leading lepton)","#Phi(sub-leading lepton)","M(ll)","p_{T}(ll)","#Delta R(ll)","#Delta #Phi(ll)","p_{T}(leading jet)","#eta(leading jet)","#Phi(leading jet)","Number of jets","Number of b-tagged jets","MET","#Phi(MET)","Number of vertices", "M(ll) [z window]", "top mass", "#Delta #Phi(ll, top)", "#Delta R(ll, top)", "|pt_top - pt_l1|/(pt_top + pt_l1)", "p_{T}(top)", "BDT","muTkIso","leadElePt","leadEleEta","leadElePhi","leadMuPt",  "leadMuEta", "leadMuPhi",
"mu1bDr",       "mu2bDr",       "mu1jDrMin",    "mu2jDrMin",    "PtRelMu1jet",  "PtRelMu2jet",  "ele1bDr",      "ele2bDr",      "ele1jDrMin",   "ele2jDrMin",   "PtRelEle1jet", "PtRelEle2jet"]#,"llMZwPt0to100","llMZwPt100to200","llMZwPt200to350","llMZwPt350toInf"]
HistAddress = '/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/hists/'


#{'2018ee': 1.138521042076288, '2017emu': 1.3261253697472475, '2016preVFPee': 1.1669105550937087, '2017mumu': 1.3442823492649438, '2016postVFPemu': 1.462483170835793, 'Allmumu': 1.255328224389967, 'Allemu': 1.2366488112328748, '2016preVFPmumu': 1.1121532788371327, 'Allee': 1.2182473496657449, '2016postVFPmumu': 1.5249467869094302, '2016postVFPee': 1.4025781380297742, '2016preVFPemu': 1.139202966971701, '2018mumu': 1.1809180192474082, '2017ee': 1.3082136332808982, '2018emu': 1.159525771115168}
SFDY={'2018ee': 1.1365510987619387, '2017emu': 1.3261806003138998, '2016preVFPee': 1.1633451628199767, '2017mumu': 1.3465242028877409, '2016postVFPemu': 1.462616384646022, 'Allmumu': 1.2578307798620003, 'Allemu': 1.2367413449370244, '2016preVFPmumu': 1.1158537954764174, 'Allee': 1.216005506276805, '2016postVFPmumu': 1.528904852150822, '2016postVFPee': 1.3992019749466855, '2016preVFPemu': 1.1393520594538817, '2018mumu': 1.183191225338832, '2017ee': 1.3061443536455783, '2018emu': 1.15963670483662}

Samples = ['data.root','WJets.root','other.root', 'DY.root', 'ttbar.root', 'tW.root', 'BNV_ST_TBCE.root', 'BNV_ST_TBUE.root' , 'BNV_ST_TDCE.root',  'BNV_ST_TDUE.root',  'BNV_ST_TSCE.root',  'BNV_ST_TSUE.root']
Samples = ['data.root','WJets.root','other.root', 'DY.root', 'ttbar.root', 'tW.root', 'STBNV_TDUE.root', 'STBNV_TDUMu.root']
SamplesName = ['Data','Jets','Others', 'DY', 't#bar{t}', 'tW','BNV tdue', 'BNV tdu#mu']# , 'BNV_ST_TBCE', 'BNV_ST_TBUE', 'BNV_ST_TDCE',  'BNV_ST_TDUE',  'BNV_ST_TSCE',  'BNV_ST_TSUE']
SamplesNameLatex = ['Data','Jets','Others', 'DY', 'tt', 'tW', 'BNV tdue', 'BNV tdu#mu']#,  'LFV-vector(emutc)', 'LFV-vector(emutu)']

colors =  [ROOT.kBlack,ROOT.kYellow,ROOT.kGreen,ROOT.kBlue-3,ROOT.kRed-4,ROOT.kOrange-3, ROOT.kPink, ROOT.kViolet-1,ROOT.kViolet+8,ROOT.kRed-5,ROOT.kSpring-1,ROOT.kGray+3]
wc1 = ROOT.WCPoint("EFTrwgt1_cS_1_cT_1")

Hists = []
for numyear, nameyear in enumerate(year):
    l0=[]
    Files = []
    for f in range(len(Samples)):
        l1=[]
        Files.append(ROOT.TFile.Open(HistAddress + nameyear+ '_' + Samples[f]))
#        print HistAddress + nameyear+ '_' + Samples[f]
        for numch, namech in enumerate(channels):
            l2=[]
            for numreg, namereg in enumerate(regions):
                l3=[]
                for numvar, namevar in enumerate(variables):
                    h= Files[f].Get(namech + '_' + namereg + '_' + namevar)
                    h = TH1EFTtoTH1(h,wc1)
                    if 'BDT' in namevar:
                        h=h.Rebin(len(bins)-1,"",bins)
                    if 'llM'== namevar:
                        h=h.Rebin(len(binsmll)-1,"",binsmll)
                    if 'lep1Pt'== namevar or 'lep2Pt'== namevar:
                        h=h.Rebin(len(binsptl)-1,"",binsptl)
                    if "DY" in Samples[f] and 'B1' in namereg and 'll' not in namech and '-' not in namereg:
                        h.Scale(SFDY[nameyear+namech])
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
    for numch, namech in enumerate(channels):
        for numreg, namereg in enumerate(regions):
            for numvar, namevar in enumerate(variables):
                HH=[]
                HHsignal=[]
                for f in range(len(Samples)):
                    if 'BNV' in Samples[f]:
                        text = 'EFTrwgt1_cS_20_cT_20'
                        wc1 = ROOT.WCPoint(text)
#                        Hists[numyear][f][numch][numreg][numvar].Scale(wc1)
                        HHsignal.append(Hists[numyear][f][numch][numreg][numvar])
                    else:
                        Hists[numyear][f][numch][numreg][numvar]
                        HH.append(Hists[numyear][f][numch][numreg][numvar])
#                stackPlots(HH, HHsignal, SamplesName, namech, namereg, nameyear,namevar,variablesName[numvar])
                print nameyear+'_'+namech+'_'+namereg+'_'+namevar+'= '+ str(Hists[numyear][3][numch][numreg][numvar].Integral())
    os.system('tar -cvf '+ nameyear +'.tar ' + nameyear)

#A=Hists[0][7][2][regions.index('llB1-BDT0p8to1')][variables.index('BDT')].Integral()
#B=Hists[0][7][2][regions.index('llB1-BDT0p8to1')][variables.index('l2jDrMin')].Integral()
#print str(B/(A+B))

le = '\\documentclass{article}' + "\n"
le += '\\usepackage{rotating}' + "\n"
le += '\\usepackage{rotating}' + "\n"
le += '\\begin{document}' + "\n"

print le
for numyear, nameyear in enumerate(year):
    for numch, namech in enumerate(channels):
        cutFlowTable(Hists, SamplesNameLatex, regionsName, numch, numyear, nameyear + ' ' + namech, 6 )
print '\\end{document}' + "\n"



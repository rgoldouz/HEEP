#define MyAnalysis_cxx
#include "MyAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "triggerEffAnalysis.h"
#include "PU_reWeighting.h"
#include "sumOfWeights.h"
#include "sumOfWeightsSignal.h"
#include "lepton_candidate.h"
#include "jet_candidate.h"
#include "XYMETCorrection_withUL17andUL18andUL16.h"
#include "lester_mt2_bisect.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TDirectory.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <time.h>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <vector>
#include "RoccoR.h"
#include "BTagCalibrationStandalone.h"
#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "CondFormats/Serialization/interface/Archive.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
//#include "Archive.h"
//#include "JetCorrectorParameters.h"
//#include "JetCorrectionUncertainty.h"
#include "GEScaleSyst.h"
#include "Utils.h"
#include "correction.h"
#include "WCPoint.h"
#include "WCFit.h"
#include "TH1EFT.h"
#include "PrescaleProvider.h"
#include <map>
#include "sys/types.h"
#include "sys/sysinfo.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "TMath.h"
#include <math.h>     
#endif

using namespace correction;
using namespace std;

void resetVec(std::vector<int> &K){
  for (int i=0;i<K.size();++i){  
    K[i]=-1;
  }
}

bool useRegions(std::vector<int> K){
  bool ifuse=false;
  for (int i=0;i<K.size();++i){
    if(K[i]>=0) ifuse=true;
  }
  return ifuse;
}

float DYMassPowheg(float M){
  return (1.0678 - 0.000120666*M + 0.0000000322*pow(M,2) - 0.0000000000039*pow(M,3));
}

float topPtPowhegData(float pt){
//Data top pt reweighting
  if (pt<500) return TMath::Exp(0.0615- (0.0005*pt));
  else return TMath::Exp(0.0615- (0.0005*500));
//MC top pt reweighting
//  return (0.973 - (0.000134 * pt) + (0.103 * exp(pt * (-0.0118))));
}

float topPtPowhegMC(float pt){
//MC top pt reweighting
  return (0.973 - (0.000134 * pt) + (0.103 * exp(pt * (-0.0118))));
}

float topPtMGLO(float x){
  return (0.688 -  0.0000174*x + 0.824*exp(-0.0000253*x)/(pow(x,0.2185)));
}

int vInd(std::map<TString, std::vector<float>> V, TString name){
  return V.find(name)->second.at(0);
}

bool ifSysNeeded(std::vector<lepton_candidate*> *lep, float cut){
  bool pass=false;
  if (lep->size()==2){
    if(((*lep)[0]->p4_ + (*lep)[1]->p4_).M()> cut) pass=true;
  }
  return pass;
}     

int getVecPos(std::vector<TString> vec, string element){
    int i;
    for(i = 0; i < vec.size(); i++){
        if(vec[i] == element){
            break;
        }
    }
    if(i == vec.size()){
        std::cout<<"No such element as "<<element<<" found. Please enter again: ";
        std::cin>>element;
        i = getVecPos(vec, element);
    }

    return i;
}

void MyAnalysis::Loop(TString fname, TString data, TString dataset ,string year, TString RUN, float xs, float lumi, float Nevent, int iseft, int nRuns, MyAnalysis *Evt){
  triggerEffAnalysis triggerEffA(Evt,year);
  bool ifSys=false;
  bool ifCR=false;
  double memoryInit=getValue();
  TH1EFT  *crossSection = new TH1EFT("crossSection","crossSection",1,0,1);
  TH1F  *eleEffNum = new TH1F("eleEffNum","eleEffNum",10,0,1000);
  TH1F  *eleEffDen = new TH1F("eleEffDen","eleEffDen",10,0,1000);
  TH1F  *truePU = new TH1F("truePU","truePU",100,0,100);

  TH2F  btagEff_b_H;
  TH2F  btagEff_c_H;
  TH2F  btagEff_udsg_H;
  TH2F  sf_triggeree_H;
  TH2F  sf_triggeremu_H;
  TH2F  sf_triggermumu_H;
  TH2F  jetVetoMaps_H;
  TH2F  highPtMuRecoSF_pVsAbsEta_H;
  TH2F  eff_triggermumu_mc_H;
  TH2F  eff_triggermumu_data_H;
  TH2F  eff_triggermumu_dataUp_H;
  TH2F  eff_triggermumu_dataDown_H;
  std::string rochesterFile;
  std::string triggerDataFile;
  std::string btagFile;
  GEScaleSyst *GE = new GEScaleSyst();
  PU wPU;
  RoccoR  rc;
  if(year == "2016preVFP")    rochesterFile = "/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/RoccoR2016aUL.txt";
  if(year == "2016postVFP")    rochesterFile = "/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/RoccoR2016bUL.txt";
  if(year == "2017")    rochesterFile = "/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/RoccoR2017UL.txt";
  if(year == "2018")    rochesterFile = "/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/RoccoR2018UL.txt";
  rc.init(rochesterFile);

  if(data == "mc"){
    TFile *f_btagEff_Map = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/btagEff.root");
    btagEff_b_H = *(TH2F*)f_btagEff_Map->Get((year + "_h2_BTaggingEff_b").c_str());
    btagEff_c_H = *(TH2F*)f_btagEff_Map->Get((year +"_h2_BTaggingEff_c").c_str());
    btagEff_udsg_H = *(TH2F*)f_btagEff_Map->Get((year +"_h2_BTaggingEff_udsg").c_str());

    TFile *f_trigger = new TFile(("/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/TriggerSF_" + year + "_UL.root").c_str());
    sf_triggeree_H = *(TH2F*)f_trigger->Get("h2D_SF_ee_lepABpt_FullError");
    sf_triggeremu_H = *(TH2F*)f_trigger->Get("h2D_SF_emu_lepABpt_FullError");
    sf_triggermumu_H = *(TH2F*)f_trigger->Get("h2D_SF_mumu_lepABpt_FullError");

    TFile *f_HighPtMuRecoSF = new TFile(("/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/HighPtMuRecoSF_" + year + ".root").c_str());
    highPtMuRecoSF_pVsAbsEta_H = *(TH2F*)f_HighPtMuRecoSF->Get("h2_HighPtMuRecoSF_pVsAbsEta");

    f_btagEff_Map->Close();
    f_trigger->Close();
    f_HighPtMuRecoSF->Close();
    delete f_btagEff_Map;
    delete f_trigger;
    delete f_HighPtMuRecoSF;
  }

  string eleSF="";
  string muSF="";
  string bSF="";
  string jetSF="";
  if(year == "2016preVFP"){
    eleSF="/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/data/POG/EGM/2016preVFP_UL/electron.json.gz";
    muSF= "/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/data/POG/MUO/2016preVFP_UL/muon_Z.json.gz";
    bSF="/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/data/POG/BTV/2016preVFP_UL/btagging.json.gz";
    jetSF="/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/data/POG/JME/2016preVFP_UL/UL16preVFP_jmar.json.gz";
    TFile *Map2016preVFP = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/hotjets-UL16.root");
    jetVetoMaps_H = *(TH2F*)Map2016preVFP->Get("h2hot_ul16_plus_hbm2_hbp12_qie11");
    Map2016preVFP->Close();
    delete Map2016preVFP;
    TFile *f_triggerEffMuMu = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root");
    eff_triggermumu_mc_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2016_DY_var");
    eff_triggermumu_data_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2016_Data_var");
    eff_triggermumu_dataUp_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2016_Data_errorUpper");
    eff_triggermumu_dataDown_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2016_Data_errorLower");
    f_triggerEffMuMu->Close();
    delete f_triggerEffMuMu;
    triggerDataFile="/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/triggerData/triggerData2016/triggerData2016";
  }
  if(year == "2016postVFP"){
    eleSF="/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/data/POG/EGM/2016postVFP_UL/electron.json.gz";
    muSF= "/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/data/POG/MUO/2016postVFP_UL/muon_Z.json.gz";
    bSF="/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/data/POG/BTV/2016postVFP_UL/btagging.json.gz";
    jetSF="/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/data/POG/JME/2016postVFP_UL/UL16postVFP_jmar.json.gz";
    TFile *Map2016postVFP = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/hotjets-UL16.root");
    jetVetoMaps_H = *(TH2F*)Map2016postVFP->Get("h2hot_ul16_plus_hbm2_hbp12_qie11");
    Map2016postVFP->Close();
    delete Map2016postVFP;
    TFile *f_triggerEffMuMu = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root");
    eff_triggermumu_mc_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2016_DY_var");
    eff_triggermumu_data_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2016_Data_var");
    eff_triggermumu_dataUp_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2016_Data_errorUpper");
    eff_triggermumu_dataDown_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2016_Data_errorLower");
    f_triggerEffMuMu->Close();
    delete f_triggerEffMuMu;
    triggerDataFile="/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/triggerData/triggerData2016/triggerData2016";
  }
  if(year == "2017"){
    eleSF="/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/data/POG/EGM/2017_UL/electron.json.gz";
    muSF= "/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/data/POG/MUO/2017_UL/muon_Z.json.gz";
    bSF="/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/data/POG/BTV/2017_UL/btagging.json.gz";
    jetSF="/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/data/POG/JME/2017_UL/UL17_jmar.json.gz";
    TFile *Map2017 = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/hotjets-UL17_v2.root");
    jetVetoMaps_H = *(TH2F*)Map2017->Get("h2hot_ul17_plus_hep17_plus_hbpw89");
    Map2017->Close();
    delete Map2017;
    TFile *f_triggerEffMuMu = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root");
    eff_triggermumu_mc_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2017_DY_var");
    eff_triggermumu_data_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2017_Data_var");
    eff_triggermumu_dataUp_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2017_Data_errorUpper");
    eff_triggermumu_dataDown_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2017_Data_errorLower");
    f_triggerEffMuMu->Close();
    delete f_triggerEffMuMu;
    triggerDataFile="/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/triggerData/triggerData2017/triggerData2017";
  }
  if(year == "2018"){
    eleSF="/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/data/POG/EGM/2018_UL/electron.json.gz";
    muSF= "/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/data/POG/MUO/2018_UL/muon_Z.json.gz";
    bSF="/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/data/POG/BTV/2018_UL/btagging.json.gz";
    jetSF="/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/data/POG/JME/2018_UL/UL18_jmar.json.gz";
    TFile *Map2018 = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/hotjets-UL18.root");
    jetVetoMaps_H = *(TH2F*)Map2018->Get("h2hot_ul18_plus_hem1516_and_hbp2m1");
    Map2018->Close();
    delete Map2018;
    TFile *f_triggerEffMuMu = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root");
    eff_triggermumu_mc_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2018_DY_var");
    eff_triggermumu_data_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2018_Data_var");
    eff_triggermumu_dataUp_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2018_Data_errorUpper");
    eff_triggermumu_dataDown_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2018_Data_errorLower");
    f_triggerEffMuMu->Close();
    delete f_triggerEffMuMu;
    triggerDataFile="/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/triggerData/triggerData2018/triggerData2018";
  }


  PrescaleProvider psProv(triggerDataFile);
  map<string,float> HEEPSF_B{ {"2016preVFP",0.985}, { "2016postVFP",0.985},{"2017",0.979},{"2018",0.973} };
  map<string,float> HEEPSF_E{ {"2016preVFP",0.990}, { "2016postVFP",0.990},{"2017",0.987},{"2018",0.980} };

  auto csetFileEleSF = CorrectionSet::from_file(eleSF);
  auto csetEleIdReco = csetFileEleSF->at("UL-Electron-ID-SF");

  auto csetFileMuSF = CorrectionSet::from_file(muSF);
  auto csetMuTightId = csetFileMuSF->at("NUM_HighPtID_DEN_TrackerMuons");
  auto csetMuTightRelIso = csetFileMuSF->at("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut");

  auto csetFilebSF = CorrectionSet::from_file(bSF);
  auto csetLightJetSF = csetFilebSF->at("deepJet_incl");
//  auto csetBcJetSF = csetFilebSF->at("deepCSV_mujets");
  auto csetBcJetSF = csetFilebSF->at("deepJet_comb");

  auto csetFileJetSF = CorrectionSet::from_file(jetSF);
  auto csetJetPuID = csetFileJetSF->at("PUJetID_eff");
  TRandom3 Tr;

  Double_t ptBins[11] = {30., 40., 60., 80., 100., 150., 200., 300., 400., 500., 1000.};
  Double_t etaBins [4]= {0., 0.6, 1.2, 2.4};
  TH2D *h2_BTaggingEff_Denom_b    = new TH2D("h2_BTaggingEff_Denom_b"   , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
  TH2D *h2_BTaggingEff_Denom_c    = new TH2D("h2_BTaggingEff_Denom_c"   , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
  TH2D *h2_BTaggingEff_Denom_udsg = new TH2D("h2_BTaggingEff_Denom_udsg", ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
  TH2D *h2_BTaggingEff_Num_b      = new TH2D("h2_BTaggingEff_Num_b"     , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
  TH2D *h2_BTaggingEff_Num_c      = new TH2D("h2_BTaggingEff_Num_c"     , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
  TH2D *h2_BTaggingEff_Num_udsg   = new TH2D("h2_BTaggingEff_Num_udsg"  , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);


  typedef vector<std::shared_ptr<TH1EFT>> Dim1;
  typedef vector<Dim1> Dim2;
  typedef vector<Dim2> Dim3;
  typedef vector<Dim3> Dim4;

  std::vector<TString> regions{"ll","llOffZ","llB1", "llBg1"};
  if (ifCR && !ifSys ) {
    regions.push_back("llB1InZ");
    regions.push_back("llB1-BDTm1to0");
    regions.push_back("llB1-BDT0to0p8");
    regions.push_back("llB1-BDT0p8to1");
    regions.push_back("llptl1G300");
    regions.push_back("llB1ptl1G300");
    regions.push_back("ll-BB");
    regions.push_back("ll-BE");
    regions.push_back("ll-EE");
    regions.push_back("ll-leadMuP");
    regions.push_back("ll-leadMuN");
  }
  std::vector<TString> channels{"OSee", "OSemu", "OSmumu","SSee", "SSemu", "SSmumu","SSefe","SSfefe","SSfemu"};

  const std::map<TString, std::vector<float>> vars =
  {
    {"lep1Pt",                         {0,      60,   0,  1500}},
    {"lep1Eta",                        {1,      20,   -3, 3   }},
    {"lep1Phi",                        {2,      25,   -4, 4   }},
    {"lep2Pt",                         {3,      25,   0,  1000}},
    {"lep2Eta",                        {4,      20,   -3, 3   }},
    {"lep2Phi",                        {5,      25,   -4, 4   }},
    {"llM",                            {6,      100,    0, 2000 }},
    {"llPt",                           {7,      20,    0, 200 }},
    {"llDr",                           {8,      25,    0, 7   }},
    {"llDphi",                         {9,      15,    0, 4   }},
    {"jet1Pt",                         {10,     20,    0, 300 }},
    {"jet1Eta",                        {11,     20,    -3, 3  }},
    {"jet1Phi",                        {12,     25,    -4, 4  }},
    {"njet",                           {13,     10,    0, 10  }},
    {"nbjet",                          {14,     6,     0, 6   }},
    {"Met",                            {15,     30,    0, 210 }},
    {"MetPhi",                         {16,     20,    -4, 4  }},
    {"nVtx",                           {17,     70,    0, 70  }},
    {"llMZw",                          {18,     80,    70, 110}},
    {"BDT",                            {19,     100,    -1, 1  }},
    {"topMass",                        {20,     25,    0, 500 }},
    {"topL1Dphi",                      {21,     15,    0, 4   }},
    {"topL1Dr",                        {22,     25,    0, 7   }},
    {"topL1DptOsumPt",                 {23,     20,    0, 1   }},
    {"topPt",                          {24,     25,    0, 500 }},
//   {"muTkIso",                        {25,     20,    0, 0.2 }},
//   {"leadElePt",                      {26,     60,    0, 1500 }},
//    {"leadEleEta",                     {27,     20,   -3, 3   }},
//    {"leadElePhi",                     {28,     25,   -4, 4   }},
//    {"leadMuPt",                       {29,     60,    0, 1500 }},
//    {"leadMuEta",                      {30,     20,   -3, 3   }},
//    {"leadMuPhi",                      {31,     25,   -4, 4   }},
//    {"mu1bDr",                          {32,      25,    0, 7   }},
//    {"mu2bDr",                          {33,      25,    0, 7   }},
//    {"mu1jDrMin",                       {34,      25,    0, 1   }},
//    {"mu2jDrMin",                       {35,      25,    0, 1   }},
//    {"PtRelMu1jet",                     {36,     30,    0, 300 }},
//    {"PtRelMu2jet",                     {37,     30,    0, 300 }},
//    {"ele1bDr",                         {38,      25,    0, 7   }},
//    {"ele2bDr",                         {39,      25,    0, 7   }},
//    {"ele1jDrMin",                      {40,      25,    0, 1   }},
//    {"ele2jDrMin",                      {41,      25,    0, 1   }},
//    {"PtRelEle1jet",                    {42,     30,    0, 300 }},
//    {"PtRelEle2jet",                    {43,     30,    0, 300 }},
  };

//  D3HistsContainer Hists;
  Hists.resize(channels.size());
  for (int i=0;i<channels.size();++i){
    Hists[i].resize(regions.size());
    for (int k=0;k<regions.size();++k){
      Hists[i][k].resize(vars.size());
    }
  }

  std::stringstream name;
  TH1EFT *h_test;

  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
      for( auto it = vars.cbegin() ; it != vars.cend() ; ++it ){
        name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first;
        h_test = new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3));
        h_test->StatOverflows(kTRUE);
        h_test->Sumw2(kTRUE);
        Hists[i][k][it->second.at(0)] = h_test;
        name.str("");
      }
    }
  }

  std::vector<string> wc_names_lst={};
  std::vector<string> wc_names_lst_HEEP={"cT", "cS"};
  std::vector<string> wc_names_lst_FCNC={"ctlS", "cte", "ctl", "ctlT", "ctZ", "cpt", "cpQM", "ctA", "cQe", "ctG", "cQlM"};
  if (fname.Contains("HEEP")) wc_names_lst = wc_names_lst_HEEP;
  if (fname.Contains("FCNC")) wc_names_lst = wc_names_lst_FCNC;

  std::vector<TString> sys{"eleRecoSf", "eleIDSf", "JetPuID", "muIdIsoSf", "bcTagSf", "LTagSf","pu", "prefiring", "trigSF","jes", "jer","muonScale","electronScale","muonRes", "unclusMET", "bcTagSfUnCorr", "LTagSfUnCorr","topPt"};

  if(data == "mc" && !fname.Contains("sys")  && ifSys){
    HistsSysUp.resize(channels.size());
    for (int i=0;i<channels.size();++i){
      HistsSysUp[i].resize(regions.size());
      for (int k=0;k<regions.size();++k){
        HistsSysUp[i][k].resize(vars.size());
        for (int n=0;n<vars.size();++n){
          HistsSysUp[i][k][n].resize(sys.size());
        }
      }
    }
  
    HistsSysDown.resize(channels.size());
    for (int i=0;i<channels.size();++i){
      HistsSysDown[i].resize(regions.size());
      for (int k=0;k<regions.size();++k){
        HistsSysDown[i][k].resize(vars.size());
        for (int n=0;n<vars.size();++n){
          HistsSysDown[i][k][n].resize(sys.size());
        }
      }
    }

    for (int i=0;i<channels.size();++i){
      for (int k=0;k<regions.size();++k){
        for( auto it = vars.cbegin() ; it != vars.cend() ; ++it ){
          for (int n=0;n<sys.size();++n){
            name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first<<"_"<<sys[n]<<"_Up";
            h_test = new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3));
            h_test->StatOverflows(kTRUE);
            h_test->Sumw2(kTRUE);
            HistsSysUp[i][k][it->second.at(0)][n] = h_test;
            name.str("");
            name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first<<"_"<<sys[n]<<"_Down";
            h_test = new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3));
            h_test->StatOverflows(kTRUE);
            h_test->Sumw2(kTRUE);
            HistsSysDown[i][k][it->second.at(0)][n] = h_test;
            name.str("");
          }
        }
      }
    }
  }

 
  std::string JECFile;
  if(year == "2016preVFP")    JECFile = "/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/Summer19UL16APV_V7_MC/Summer19UL16APV_V7_MC_UncertaintySources_AK4PFchs.txt";
  if(year == "2016postVFP")   JECFile = "/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/Summer19UL16_V7_MC/Summer19UL16_V7_MC_UncertaintySources_AK4PFchs.txt";
  if(year == "2017")          JECFile = "/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/Summer19UL17_V5_MC/Summer19UL17_V5_MC_UncertaintySources_AK4PFchs.txt";
  if(year == "2018")          JECFile = "/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/Summer19UL18_V5_MC/Summer19UL18_V5_MC_UncertaintySources_AK4PFchs.txt";

  std::vector<TString> sysJecNames{"AbsoluteMPFBias","AbsoluteScale","AbsoluteStat","FlavorQCD","Fragmentation","PileUpDataMC","PileUpPtBB","PileUpPtEC1","PileUpPtEC2","PileUpPtHF","PileUpPtRef","RelativeFSR","RelativePtBB","RelativePtEC1","RelativePtEC2","RelativePtHF","RelativeBal","RelativeSample","RelativeStatEC","RelativeStatFSR","RelativeStatHF","SinglePionECAL","SinglePionHCAL","TimePtEta"};
  const int nsrc = 24;
  const char* srcnames[nsrc] = {"AbsoluteMPFBias","AbsoluteScale","AbsoluteStat","FlavorQCD","Fragmentation","PileUpDataMC","PileUpPtBB","PileUpPtEC1","PileUpPtEC2","PileUpPtHF","PileUpPtRef","RelativeFSR","RelativePtBB","RelativePtEC1","RelativePtEC2","RelativePtHF","RelativeBal","RelativeSample","RelativeStatEC","RelativeStatFSR","RelativeStatHF","SinglePionECAL","SinglePionHCAL","TimePtEta"};
  std::vector<JetCorrectionUncertainty*> vsrc(nsrc);
  for (int isrc = 0; isrc < nsrc; isrc++) {
    JetCorrectorParameters *p = new JetCorrectorParameters(JECFile, srcnames[isrc]);
    JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
    vsrc[isrc] = unc;
  }
 
  JetCorrectorParameters *pJes = new JetCorrectorParameters(JECFile, "Total");
  JetCorrectionUncertainty *uncJes = new JetCorrectionUncertainty(*pJes);

  JetCorrectorParameters *pJer = new JetCorrectorParameters(JECFile, "RelativeJEREC1");
  JetCorrectionUncertainty *uncJer = new JetCorrectionUncertainty(*pJer);

  if(data == "mc" && !fname.Contains("sys") && ifSys){
    HistsJecUp.resize(channels.size());
    for (int i=0;i<channels.size();++i){
      HistsJecUp[i].resize(regions.size()-2);
      for (int k=0;k<regions.size()-2;++k){
        HistsJecUp[i][k].resize(1);
        for (int n=0;n<1;++n){
          HistsJecUp[i][k][n].resize(sysJecNames.size());
        }
      }
    }
  
    HistsJecDown.resize(channels.size());
    for (int i=0;i<channels.size();++i){
      HistsJecDown[i].resize(regions.size()-2);
      for (int k=0;k<regions.size()-2;++k){
        HistsJecDown[i][k].resize(1);
        for (int n=0;n<1;++n){
          HistsJecDown[i][k][n].resize(sysJecNames.size());
        }
      }
    }
    for (int i=0;i<channels.size();++i){
      for (int k=0;k<regions.size()-2;++k){
        for( auto it = vars.cbegin() ; it != vars.cend() ; ++it ){
          if(it->first !="BDT") continue;
          for (int n=0;n<sysJecNames.size();++n){
            name<<channels[i]<<"_"<<regions[k+2]<<"_"<<it->first<<"_"<<sysJecNames[n]<<"_Up";
            h_test = new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3));
            h_test->StatOverflows(kTRUE);
            h_test->Sumw2(kTRUE);
            HistsJecUp[i][k][0][n] = h_test;
            name.str("");
            name<<channels[i]<<"_"<<regions[k+2]<<"_"<<it->first<<"_"<<sysJecNames[n]<<"_Down";
            h_test = new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3));
            h_test->StatOverflows(kTRUE);
            h_test->Sumw2(kTRUE);
            HistsJecDown[i][k][0][n] = h_test;
            name.str("");
          }
        }
      }
    }
  }
// [0] is renscfact=0.5d0 facscfact=0.5d0 ; [1] is renscfact=0.5d0 facscfact=1d0 ; [2] is renscfact=0.5d0 facscfact=2d0 ; [3] is renscfact=1d0 facscfact=0.5d0 ; [4] is renscfact=1d0 facscfact=1d0 ; [5] is renscfact=1d0 facscfact=2d0 ; [6] is renscfact=2d0 facscfact=0.5d0 ; [7] is renscfact=2d0 facscfact=1d0 ; [8] is renscfact=2d0 facscfact=2d0 *
  int nScale = 9;
// LHA IDs NNPDF31_nnlo_hessian_pdfas 306000 - 306102*
  int nPdf = 100;
//[0] is ISR=2 FSR=1; [1] is ISR=1 FSR=2[2] is ISR=0.5 FSR=1; [3] is ISR=1 FSR=0.5;*
  int nPS = 4;

  if(data == "mc" && ifSys){
    if ((fname.Contains("TTTo2L2Nu") && !fname.Contains("sys")) || fname.Contains("HEEP")){
    HistsSysReweightsQscale.resize(channels.size());
    for (int i=0;i<channels.size();++i){
      HistsSysReweightsQscale[i].resize(2);
      for (int k=0;k<2;++k){
        HistsSysReweightsQscale[i][k].resize(1);
        for (int n=0;n<1;++n){
          HistsSysReweightsQscale[i][k][n].resize(nScale);
        }
      }
    }
  
    HistsSysReweightsPDF.resize(channels.size());
    for (int i=0;i<channels.size();++i){
      HistsSysReweightsPDF[i].resize(2);
      for (int k=0;k<2;++k){
        HistsSysReweightsPDF[i][k].resize(1);
        for (int n=0;n<1;++n){
          HistsSysReweightsPDF[i][k][n].resize(nPdf);
        }
      }
    }
  
    HistsSysReweightsPS.resize(channels.size());
    for (int i=0;i<channels.size();++i){
      HistsSysReweightsPS[i].resize(2);
      for (int k=0;k<2;++k){
        HistsSysReweightsPS[i][k].resize(1);
        for (int n=0;n<1;++n){
          HistsSysReweightsPS[i][k][n].resize(nPS);
        }
      }
    }
      for (int i=0;i<channels.size();++i){
        for (int k=2;k<regions.size();++k){
        for( auto it = vars.cbegin() ; it != vars.cend() ; ++it ){
            if(it->first !="BDT") continue;
            for (int n=0;n<nScale;++n){
              name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first<<"_Qscale_"<<n;
              h_test = new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3));
              h_test->StatOverflows(kTRUE);
              h_test->Sumw2(kTRUE);
              HistsSysReweightsQscale[i][k-2][0][n] = h_test;
              name.str("");
            }
            for (int n=0;n<nPdf;++n){
              name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first<<"_PDF_"<<n;
              h_test = new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3));
              h_test->StatOverflows(kTRUE);
              h_test->Sumw2(kTRUE);
              HistsSysReweightsPDF[i][k-2][0][n] = h_test;
              name.str("");
            }
            for (int n=0;n<nPS;++n){
              name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first<<"_PS_"<<n;
              h_test = new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3));
              h_test->StatOverflows(kTRUE);
              h_test->Sumw2(kTRUE);
              HistsSysReweightsPS[i][k-2][0][n] = h_test;
              name.str("");
            }
          }
        }
      }
    }
  }
//  TFile file_out ("ANoutput.root","RECREATE");
  float lep1Pt_;
  float lep1Eta_;
  float lep2Pt_;
  float lep2Eta_;
  float llM_;
  float llPt_;
  float llDr_;
  float llDphi_;
  float jet1Pt_;
  float jet1Eta_;
  float topMass_;
  float topL1Dphi_;
  float topL1Dr_;
  float topL1DptOsumPt_;
  float topPt_;
  float weight_;
  float BDToutput_;

  TTree tree_out("HEEP","Top HEEP analysis") ;
  tree_out.Branch("lep1Pt"      , &lep1Pt_ , "lep1Pt/F" ) ;
  tree_out.Branch("lep1Eta"      , &lep1Eta_ , "lep1Eta/F" ) ;
  tree_out.Branch("lep2Pt"      , &lep2Pt_ , "lep2Pt/F" ) ;
  tree_out.Branch("lep2Eta"      , &lep2Eta_ , "lep2Eta/F" ) ;
  tree_out.Branch("llM"      , &llM_ , "llM/F" ) ;
  tree_out.Branch("llPt"      , &llPt_ , "llPt/F" ) ;
  tree_out.Branch("llDr"      , &llDr_ , "llDr/F" ) ;
  tree_out.Branch("llDphi"      , &llDphi_ , "llDphi/F" ) ;
  tree_out.Branch("jet1Pt"      , &jet1Pt_ , "jet1Pt/F" ) ;
  tree_out.Branch("jet1Eta"      , &jet1Eta_ , "jet1Eta/F"  ) ;
  tree_out.Branch("topMass"      , &topMass_ , "topMass/F" ) ;
  tree_out.Branch("topL1Dphi"      , &topL1Dphi_ , "topL1Dphi/F" ) ;
  tree_out.Branch("topL1Dr"      , &topL1Dr_ , "topL1Dr/F" ) ;
  tree_out.Branch("topL1DptOsumPt"      , &topL1DptOsumPt_ , "topL1DptOsumPt/F" ) ;
  tree_out.Branch("topPt"      , &topPt_ , "topPt/F" ) ;
  tree_out.Branch("weight"      , &weight_ , "weight/F" ) ;
  tree_out.Branch("BDToutput"      , &BDToutput_, "BDToutput/F" ) ;

  TMVA::Tools::Instance();
  TMVA::Reader *readerMVA = new TMVA::Reader( "!Color:!Silent" );
  Float_t MVA_lep1Pt;
  Float_t MVA_lep2Pt;
  Float_t MVA_llM;
  Float_t MVA_llPt;
  Float_t MVA_llDr;
  Float_t MVA_llDphi;
  Float_t MVA_topL1Dphi;
  Float_t MVA_topL1Dr;
  Float_t MVA_topL1DptOsumPt;
  Float_t MVA_topPt;
  readerMVA->AddVariable ("lep1Pt"      , &MVA_lep1Pt) ;
  readerMVA->AddVariable ("lep2Pt"      , &MVA_lep2Pt) ;
  readerMVA->AddVariable ("llM"      , &MVA_llM) ;
  readerMVA->AddVariable ("llPt"      , &MVA_llPt) ;
  readerMVA->AddVariable ("llDr"      , &MVA_llDr) ;
  readerMVA->AddVariable ("llDphi"      , &MVA_llDphi) ;
  readerMVA->AddVariable ("topL1Dphi"      , &MVA_topL1Dphi) ;
  readerMVA->AddVariable ("topL1Dr"      , &MVA_topL1Dr) ;
  readerMVA->AddVariable ("topL1DptOsumPt"      , &MVA_topL1Dr) ;
  readerMVA->AddVariable ("topPt"      , &MVA_topPt) ;
  readerMVA->BookMVA( "BDTG", "/afs/crc.nd.edu/user/r/rgoldouz/HEEP/NanoAnalysis/input/TMVAClassification_BDTG.weights.xml");

  std::vector<lepton_candidate*> *selectedLeptons;
  std::vector<lepton_candidate*> *selectedLeptonsMuScaleUp;
  std::vector<lepton_candidate*> *selectedLeptonsMuScaleDown;
  std::vector<lepton_candidate*> *selectedLeptonsEleScaleUp;
  std::vector<lepton_candidate*> *selectedLeptonsEleScaleDown;
  std::vector<lepton_candidate*> *selectedLeptonsMuResUp;
  std::vector<lepton_candidate*> *selectedLeptonsMuResDown;
  std::vector<jet_candidate*> *selectedJets;
  std::vector<jet_candidate*> *selectedJetsJerUp;
  std::vector<jet_candidate*> *selectedJetsJerDown;
  std::vector<jet_candidate*> *selectedJetsJesUp;
  std::vector<jet_candidate*> *selectedJetsJesDown;
  std::vector<jet_candidate*> *JECJetsUp;
  std::vector<jet_candidate*> *JECJetsDown;
  std::vector<std::vector<jet_candidate*>> *JECsysUp;
  std::vector<std::vector<jet_candidate*>> *JECsysDown;
  std::vector<int> *JECsysNbtagUp;
  std::vector<int> *JECsysNbtagDown;
  std::vector<float> *JECsysMETUp;
  std::vector<float> *JECsysMETDown;
  std::vector<float> *JECsysMETPhiUp;
  std::vector<float> *JECsysMETPhiDown;
  float MetJetsPtJesUp;
  float MetJetsPtJesDown;
  float MetJetsPtJerUp;
  float MetJetsPtJerDown;
  float MetJetsPhiJesUp;
  float MetJetsPhiJesDown;
  float MetJetsPhiJerUp;
  float MetJetsPhiJerDown;
  float MetXJetsJesUp;
  float MetXJetsJesDown;
  float MetXJetsJerUp;
  float MetXJetsJerDown;
  float MetYJetsJesUp;
  float MetYJetsJesDown;
  float MetYJetsJerUp;
  float MetYJetsJerDown;
  WCFit *eft_fit;
  std::pair<double,double> METXYCorr;

  TLorentzVector wp, wm, b, ab, top, atop, lep, alep;
  std::vector<float> nominalWeights;
  TLorentzVector recoTop, recoBjet, recoW, recoNu, recoL1, recoL2, highPtMu;
  nominalWeights.assign(sys.size(), 1);
  std::vector<float> sysUpWeights;
  sysUpWeights.assign(sys.size(), 1);
  std::vector<float> sysDownWeights;
  sysDownWeights.assign(sys.size(), 1);
  bool leptonPass;
  bool triggerPass;
  bool DyPass;
  bool DyInZPass;
  bool MetPass;
  bool triggerPassEE;
  bool triggerPassEMu;
  bool triggerPassMuMu;
  bool triggerPassTrigEff;
  bool metFilterPass;
  bool ifTopPt=false;
  float ttKFactor;
  int ch;
  float MetCut=60;
  float sf_Ele_Reco;
  float sf_Ele_ID;
  float sf_Mu_ID;
  float sf_Mu_ISO;
  float sf_Mu_RECO;
  float sf_Trigger;
  float sf_JetPuId;
  float weight_PU;
  float weight_Lumi;
  float weight_lep;
  float weight_lepB;
  float weight_EFT;
  float weight_prefiring;
  float weight_topPtPowhegData;
  float weight_topPtMGLO;
  float MVAoutputJerUp;
  float MVAoutputJerDown;
  float MVAoutputJesUp;
  float MVAoutputJesDown;
  float MVAoutputMuScaleUp;
  float MVAoutputMuScaleDown;
  float MVAoutputEleScaleUp;
  float MVAoutputEleScaleDown;
  float MVAoutputMuResUp;
  float MVAoutputMuResDown;
  float MVAoutputUnclusMETUp;
  float MVAoutputUnclusMETDown;
  float metUnclusMETUp;
  float metUnclusMETDown;
  float metUnclusMETPhiUp;
  float metUnclusMETPhiDown;
  float MVAoutput;
  float mT2output;
  double P_bjet_data;
  double P_bjet_mc;
  int nAccept=0;
  float sumPuWeight=0;
  float sumPreFireWeight=0;
  float sumWeightMuID=0;
  float sumWeightMuIso=0;
  float sumWeighttTrigger=0;
  float sumWeightEleID=0;
  float sumWeightEleReco=0;
  float sumWeightBtag=0;
  float correctMuonPt=0;
  float correctMuonPtUp=0;
  float correctMuonPtDown=0;
  int nbjetGen;
  int nbjet;
  int lbjet;
  int nbjetJesUp;
  int nbjetJesDown;
  int nbjetJerUp;
  int nbjetJerDown;
  float JECMETUpx;
  float JECMETUpy;
  float JECMETDownx;
  float JECMETDowny;
  float pt_res;
  double muPtSFRochester;
  int R;
  double sup = 0;
  double sdw = 0;
  bool jetlepfail;
  float BJetSF;
  float CJetSF;
  float LJetSF;
  float BJetSF_UpCorr;
  float CJetSF_UpCorr;
  float LJetSF_UpCorr;
  float BJetSF_UpUnCorr;
  float CJetSF_UpUnCorr;
  float LJetSF_UpUnCorr;
  float BJetSF_DownCorr;
  float CJetSF_DownCorr;
  float LJetSF_DownCorr;
  float BJetSF_DownUnCorr;
  float CJetSF_DownUnCorr;
  float LJetSF_DownUnCorr;
  float BJetEff;
  float CJetEff;
  float LJetEff;
  float SmearedMuonPt;
  float myMETpt;
  float myMETphi;
  int nLHEl;
  float mllLHE;
  int weightSign;
  float effTmc1;   
  float effTmc2;   
  float effTdata1; 
  float effTdata2; 
  float effTdata1Up;
  float effTdata2Up;
  float effTdata1Down;
  float effTdata2Down;
  std::vector<int> reg(regions.size());
  std::vector<float> wgt(regions.size());
  std::vector<WCFit> wcfit(regions.size());
  int JECnbUp;
  int JECnbDown;
  float fakeRate;
  float fi;

  if (fname.Contains("TTTo2L2Nu") || fname.Contains("sys") || fname.Contains("TTHEEP")) ifTopPt=true;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  Long64_t ntr = fChain->GetEntries ();

//Loop over events
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
//  for (Long64_t jentry=0; jentry<50000;jentry++) {
std::cout <<"hlt ps "<<psProv.hltPrescale("HLT_Photon120_v",276501,800)<<std::endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    displayProgress(jentry, ntr) ;
    triggerEffA.fillTriggerHists();
    nLHEl=0;
    if(data == "mc" && !fname.Contains("pythia")){
      for (int l=0;l<nLHEPart;l++){
        if(abs(LHEPart_pdgId[l]) ==11 || abs(LHEPart_pdgId[l]) ==13 ){
          if(LHEPart_pt[l]>20 && abs(LHEPart_eta[l])<2.4) nLHEl++;
        }
      }
      for (int l=0;l<nGenPart;l++){
        if(abs(GenPart_pdgId[l])==11 || GenPart_pdgId[l]==13){
          if(abs(GenPart_pdgId[GenPart_genPartIdxMother[l]])==24 && GenPart_pt[l]>20 &&  abs(GenPart_eta[l])<2.4) nLHEl++;
        }
      }
    }

    if(nLHEl >1) nAccept++;

    triggerPassEE = false;
    triggerPassEMu = false;
    triggerPassMuMu = false;
    metFilterPass = false;
    leptonPass = false;
    triggerPass = false;
    DyPass = false;
    DyInZPass  = false;
    MetPass = false;
    ch =99;
    ttKFactor=1;
    sf_Ele_Reco =1;
    sf_Ele_ID =1;
    sf_Mu_ID =1;
    sf_Mu_ISO =1;
    sf_Mu_RECO =1;
    sf_Trigger =1;
    sf_JetPuId =1;
    muPtSFRochester=1;
    weight_PU =1;
    weight_Lumi =1;
    weight_lep =1;
    weight_lepB =1;
    weight_EFT =1;
    weight_prefiring =1;
    weight_topPtPowhegData =1;
    weight_topPtMGLO =1;
    P_bjet_data =1;
    P_bjet_mc =1;
    MVAoutput=0;
    mT2output=0;
    nbjetGen=0;
    nbjet=0;
    lbjet=0;
    nbjetJerUp=0;
    nbjetJerDown=0;
    nbjetJesUp=0;
    nbjetJesDown=0;
    metUnclusMETUp=0;
    metUnclusMETDown=0;
    metUnclusMETPhiUp=0;
    metUnclusMETPhiDown=0;
    if(year == "2016preVFP") METXYCorr = METXYCorr_Met_MetPhi(MET_pt, MET_phi, run, "2016APV", data == "mc", PV_npvs, true, false);
    else if(year == "2016postVFP") METXYCorr = METXYCorr_Met_MetPhi(MET_pt, MET_phi, run, "2016nonAPV", data == "mc", PV_npvs, true, false);
    else METXYCorr = METXYCorr_Met_MetPhi(MET_pt, MET_phi, run, year, data == "mc", PV_npvs, true, false);
    myMETpt= METXYCorr.first;
    myMETphi= METXYCorr.second;

    BJetSF=1;
    CJetSF=1;
    LJetSF=1;
    BJetSF_UpCorr=1;
    CJetSF_UpCorr=1;
    LJetSF_UpCorr=1;
    BJetSF_UpUnCorr=1;
    CJetSF_UpUnCorr=1;
    LJetSF_UpUnCorr=1;
    BJetSF_DownCorr=1;
    CJetSF_DownCorr=1;
    LJetSF_DownCorr=1;
    BJetSF_DownUnCorr=1;
    CJetSF_DownUnCorr=1;
    LJetSF_DownUnCorr=1;
    BJetEff=1;
    CJetEff=1;
    LJetEff=1;
    mllLHE=0;
    weightSign=1;
    effTmc1=1;
    effTmc2=1;
    effTdata1=1;
    effTdata2=1;
    fakeRate=1;
    fi=1;
   
    if (data == "mc" && fname.Contains("DY50")){
      for (int l=0;l<nLHEPart;l++){
        if(LHEPart_pdgId[l]==11 || LHEPart_pdgId[l]==13 || LHEPart_pdgId[l]==15) lep.SetPtEtaPhiM(LHEPart_pt[l], LHEPart_eta[l], LHEPart_phi[l], LHEPart_mass[l]) ;
        if(LHEPart_pdgId[l]==-11 || LHEPart_pdgId[l]==-13 || LHEPart_pdgId[l]==-15) alep.SetPtEtaPhiM(LHEPart_pt[l], LHEPart_eta[l], LHEPart_phi[l], LHEPart_mass[l]) ;
      }
      mllLHE=(lep+alep).M();
      if(year == "2016preVFP" && mllLHE>200) continue;
      if(year != "2016preVFP" && year != "2016postVFP" && mllLHE>100) continue;
    }
 
    for (int n=0;n<sys.size();++n){
      nominalWeights[n] =1;
      sysUpWeights[n] =1;
      sysDownWeights[n] =1;
    }

//MET filters

    if (iseft) {
      eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, 1.0/nRuns);
//      for (UInt_t i=0;i<nWCnames;++i){
//         char ch[4];
////         for(int j = 0; j < 4; j++) ch[j] = (WCnames[i] >> (4-1-j)*8) & 0xFF;
//         cout<< " - "<<WCnames[i] <<endl;
////         for(int n=0 ; n<4 ; ++n) cout << ch[n];
////         cout<<" : "<<endl;;
//      }
      
    }
    else {
      eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0); 
    }
    crossSection->Fill(0.5, 1,*eft_fit);

    delete eft_fit; 
//You should add Flag_BadPFMuonDzFilter to the MET filter list but since it is not available in v8, lets remove it for now.
    if(year == "2017" || year == "2018"){
      if ( Flag_goodVertices  &&  Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter &&  Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_eeBadScFilter && Flag_ecalBadCalibFilter && Flag_BadPFMuonDzFilter) metFilterPass = true;
    }
    else{
      if ( Flag_goodVertices  &&  Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter &&  Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_eeBadScFilter && Flag_BadPFMuonDzFilter) metFilterPass = true;
    }

//trigger MC
    if(data == "mc" && year == "2016preVFP"){
      if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele27_WPTight_Gsf || HLT_Photon175 ) triggerPassEE =true;
      if(HLT_Mu50 || HLT_TkMu50) triggerPassEMu =true;
      if(HLT_Mu50 || HLT_TkMu50) triggerPassMuMu =true;
    }

    if(data == "mc" && year == "2016postVFP"){
      if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele27_WPTight_Gsf || HLT_Photon175 ) triggerPassEE =true;
      if(HLT_Mu50 || HLT_TkMu50) triggerPassEMu =true;
      if(HLT_Mu50 || HLT_TkMu50) triggerPassMuMu =true;
    }
    if(data == "mc" && year == "2017"){
      if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele35_WPTight_Gsf || HLT_Photon200) triggerPassEE =true;
      if(HLT_Mu50 || HLT_TkMu100 || HLT_OldMu100) triggerPassEMu =true;
      if(HLT_Mu50 || HLT_TkMu100 || HLT_OldMu100) triggerPassMuMu =true;
    }

    if(data == "mc" && year == "2018"){
      if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele32_WPTight_Gsf || HLT_Photon200) triggerPassEE =true;
      if(HLT_Mu50 || HLT_TkMu100 || HLT_Mu100) triggerPassEMu =true;
      if(HLT_Mu50 || HLT_TkMu100 || HLT_Mu100) triggerPassMuMu =true;
    }

//trigger DATA
    if(data == "data"){
      if(year == "2016preVFP"){
        if(dataset=="SingleElectron"){
          if(!(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && (HLT_Ele27_WPTight_Gsf || HLT_Photon175)) triggerPassEE =true;
        }
        if(dataset=="SingleMuon"){
          if(HLT_Mu50 || HLT_TkMu50) triggerPassMuMu =true;
          if(HLT_Mu50 || HLT_TkMu50) triggerPassEMu =true;
        }
        if(dataset=="DoubleEG"){
          if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) triggerPassEE =true;
        }
      }
      if(year == "2016postVFP"){
        if(dataset=="SingleElectron"){
          if(!(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && (HLT_Ele27_WPTight_Gsf||HLT_Photon175)) triggerPassEE =true;
        }
        if(dataset=="SingleMuon"){
          if(HLT_Mu50 || HLT_TkMu50) triggerPassMuMu =true;
          if(HLT_Mu50 || HLT_TkMu50) triggerPassEMu =true;
        }
        if(dataset=="DoubleEG"){
          if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) triggerPassEE =true;
        }
      }
      if(year == "2017"){
        if(dataset=="SingleElectron"){
          if(!HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL && (HLT_Ele35_WPTight_Gsf || HLT_Photon200)) triggerPassEE =true;
        }
        if(dataset=="SingleMuon"){
          if(HLT_Mu50 || HLT_TkMu100 || HLT_OldMu100) triggerPassMuMu =true;
          if(HLT_Mu50 || HLT_TkMu100 || HLT_OldMu100) triggerPassEMu =true;
        }
        if(dataset=="DoubleEG"){
          if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL) triggerPassEE =true;
        }
      }
      if(year == "2018"){
        if(dataset=="EGamma"){
          if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele32_WPTight_Gsf || HLT_Photon200) triggerPassEE =true;
        }
        if(dataset=="SingleMuon"){
          if(HLT_Mu50 || HLT_TkMu100 || HLT_Mu100) triggerPassMuMu =true;
          if(HLT_Mu50 || HLT_TkMu100 || HLT_Mu100) triggerPassEMu =true;
        }
      }
    }

    if(year == "2017" || HLT_Ele35_WPTight_Gsf) triggerPassTrigEff=true;
    if(year != "2017" || HLT_Ele32_WPTight_Gsf) triggerPassTrigEff=true;

//    if(!(triggerPassEE || triggerPassEMu || triggerPassMuMu)) continue;
    if(!metFilterPass) continue;

    selectedLeptons = new std::vector<lepton_candidate*>();
// electron
    for (int l=0;l<nElectron;l++){
      if(abs(Electron_eta[l]) > 2.4 || (abs(Electron_eta[l])> 1.4442 && (abs(Electron_eta[l])< 1.566))) continue;
      if(Electron_pt[l] > 35) eleEffDen->Fill(Electron_pt[l]);
//      if(Electron_cutBased[l] < 4) continue;
      if (year == "2018" && Electron_eta[l] < -1.3 && Electron_phi[l] < -0.87 && Electron_phi[l] > -1.57) continue;      
      if(!Electron_cutBased_HEEP[l]) continue;
//      if(!Electron_convVeto[l]) continue;
      //remove HEM affected region
      if(Electron_pt[l] <35) continue;
      eleEffNum->Fill(Electron_pt[l]);
      selectedLeptons->push_back(new lepton_candidate(Electron_pt[l],Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
    }
    sort(selectedLeptons->begin(), selectedLeptons->end(), ComparePtLep);
//Check if analysis cuts are passed and then categorize dilepton channels
    if(selectedLeptons->size() ==2){ 
      if ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 2 && abs((*selectedLeptons)[0]->eta_)<1.444  && abs((*selectedLeptons)[1]->eta_)<1.444 && HLT_DoublePhoton70) ch = getVecPos(channels,"OSee");
//Fill histograms
      if (data == "mc"){
        if(!fname.Contains("pythia")) weightSign = signnum_typical(LHEWeight_originalXWGTUP);
        weight_lep  = weight_lep * weight_Lumi * weightSign * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
        weight_lepB = weight_lep * (P_bjet_data/P_bjet_mc);
        weight_EFT = lumi * (1000.0/nRuns)  * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      }
      if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
      else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
      resetVec(reg);
      if(ch<30){
        reg[0]=0;
        wgt[0]=weight_lep;
        wcfit[0]=*eft_fit;
      }
      if(ch<30 && HLT_DoubleEle33_CaloIdL_MW){
        reg[1]=1;
        wgt[1]=weight_lep;
        wcfit[1]=*eft_fit;
      }
      if(leptonPass && ch<30 && nbjet==1){
        reg[2]=2;
        wgt[2]=weight_lepB;
        wcfit[2]=*eft_fit;
      }
  
      if(leptonPass && ch<30 && nbjet>1){
        reg[3]=3;
        wgt[3]=weight_lepB;
        wcfit[3]=*eft_fit;
      }
      FillD3Hists(Hists, ch, reg, vInd(vars,"lep1Pt"), (*selectedLeptons)[0]->pt_ ,wgt, wcfit);
      FillD3Hists(Hists, ch, reg, vInd(vars,"lep1Eta"), (*selectedLeptons)[0]->eta_ ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"lep1Phi"), (*selectedLeptons)[0]->phi_ ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"lep2Pt"), (*selectedLeptons)[1]->pt_ ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"lep2Eta"), (*selectedLeptons)[1]->eta_ ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"lep2Phi"), (*selectedLeptons)[1]->phi_ ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"llM"), ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"llPt"), ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt() ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"llDr"), deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_) ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"llDphi"), abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_)) ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"llMZw"), ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit);
      delete eft_fit;
    }

    for (int l=0;l<selectedLeptons->size();l++){
      delete (*selectedLeptons)[l];
    }
    selectedLeptons->clear();
    selectedLeptons->shrink_to_fit();
    delete selectedLeptons;
  }
  cout<<"Loop is completed"<<endl;
  cout<<"from "<<ntr<<" events, "<<nAccept<<" events are accepted"<<endl;
  cout<<"Total Virtual Memory just after the loop: "<<(getValue()-memoryInit)/1000.0<<" MB"<<endl;
//  cout<<"sumPuWeight "<<sumPuWeight<<endl;
//  cout<<"sumPreFireWeight "<<sumPreFireWeight<<endl;
//  cout<<"sumWeightMuID "<<sumWeightMuID<<endl;
//  cout<<"sumWeightMuIso "<<sumWeightMuIso<<endl;
//  cout<<"sumWeighttTrigger "<<sumWeighttTrigger<<endl;
//  cout<<"sumWeightEleReco "<<sumWeightEleReco<<endl;
//  cout<<"sumWeightEleID "<<sumWeightEleID<<endl;  
//  cout<<"sumWeightBtag "<<sumWeightBtag<<endl;
//cout<<"Integral of the region "<<regions[getVecPos(regions,"llB1-BmuExc")]<<":"<<Hists[3][getVecPos(regions,"llB1-BmuExc")][vInd(vars,"lep1Pt")]->Integral()<<endl;
//cout<<"Integral of the region "<<regions[getVecPos(regions,"llB1")]<<":"<<Hists[3][getVecPos(regions,"llB1")][vInd(vars,"lep1Pt")]->Integral()<<endl;

  TFile file_out ("ANoutput.root","RECREATE");
  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
      for (int l=0;l<vars.size();++l){
        Hists[i][k][l]  ->Write("",TObject::kOverwrite);
      }
    }
  }

  if(data=="mc" && !fname.Contains("sys")  && ifSys){
    for (int i=0;i<channels.size();++i){
      file_out.mkdir("sys"+channels[i]);
      file_out.cd("sys"+channels[i]+"/");
      for (int k=0;k<regions.size();++k){
        for (int l=0;l<vars.size();++l){
          for (int n=0;n<sys.size();++n){
            HistsSysUp[i][k][l][n]->Write("",TObject::kOverwrite);
            HistsSysDown[i][k][l][n]->Write("",TObject::kOverwrite);
          }
        }
      }
      file_out.cd("");
    }
  }
  file_out.cd("");
  h2_BTaggingEff_Denom_b   ->Write("",TObject::kOverwrite);
  h2_BTaggingEff_Denom_c   ->Write("",TObject::kOverwrite);
  h2_BTaggingEff_Denom_udsg->Write("",TObject::kOverwrite);
  h2_BTaggingEff_Num_b     ->Write("",TObject::kOverwrite);
  h2_BTaggingEff_Num_c     ->Write("",TObject::kOverwrite);
  h2_BTaggingEff_Num_udsg  ->Write("",TObject::kOverwrite);
  crossSection             ->Write("",TObject::kOverwrite);
  eleEffNum                ->Write("",TObject::kOverwrite);
  eleEffDen                ->Write("",TObject::kOverwrite);
  truePU                   ->Write("",TObject::kOverwrite);
  if(data=="mc" && !fname.Contains("sys")  && ifSys){
    file_out.mkdir("JECSys");
    file_out.cd("JECSys/");
    for (int i=0;i<channels.size();++i){
      for (int k=2;k<regions.size();++k){
        for (int n=0;n<sysJecNames.size();++n){
          HistsJecUp[i][k-2][0][n]->Write("",TObject::kOverwrite);
          HistsJecDown[i][k-2][0][n]->Write("",TObject::kOverwrite);
        }
      }
    }
  }
  file_out.cd("");
  if(data == "mc" && ifSys){
    if ((fname.Contains("TTTo2L2Nu")&& !fname.Contains("sys")) || fname.Contains("HEEP")){
      file_out.mkdir("reweightingSys");
      file_out.cd("reweightingSys/");
      for (int i=0;i<channels.size();++i){
        for (int k=2;k<regions.size();++k){
          for (int n=0;n<nScale;++n){
            HistsSysReweightsQscale[i][k-2][0][n]->Write("",TObject::kOverwrite);
          }
          for (int n=0;n<nPdf;++n){
            HistsSysReweightsPDF[i][k-2][0][n]->Write("",TObject::kOverwrite);
          }
          for (int n=0;n<nPS;++n){
            HistsSysReweightsPS[i][k-2][0][n]->Write("",TObject::kOverwrite);
          }
        }
      }
    }
    cout<<"Cleaning the memory"<<endl;
    for (int i=0;i<channels.size();++i){
      for (int k=0;k<regions.size();++k){
        for (int l=0;l<vars.size();++l){
          delete Hists[i][k][l];
          for (int n=0;n<sys.size();++n){
            delete HistsSysUp[i][k][l][n];
            delete HistsSysDown[i][k][l][n];
          }
        }
      }
    }
    if(data=="mc" && !fname.Contains("sys")){
      for (int i=0;i<channels.size();++i){
        for (int k=2;k<regions.size();++k){
          for (int n=0;n<sysJecNames.size();++n){
            delete HistsJecUp[i][k-2][0][n];
            delete HistsJecDown[i][k-2][0][n];
          }
        }
      }
    }
    if ((fname.Contains("TTTo2L2Nu")&& !fname.Contains("sys")) || fname.Contains("HEEP")){
      for (int i=0;i<channels.size();++i){
        for (int k=2;k<regions.size();++k){
          for (int n=0;n<nScale;++n){
            delete HistsSysReweightsQscale[i][k-2][0][n];
          }
          for (int n=0;n<nPdf;++n){
            delete HistsSysReweightsPDF[i][k-2][0][n];
          }
          for (int n=0;n<nPS;++n){
            delete HistsSysReweightsPS[i][k-2][0][n];
          }
        }
      }
    }
  }
  file_out.cd("");
  cout<<"All Analysis histograms are written on the output"<<endl;
  triggerEffA.writeTriggerHists();
  triggerEffA.cleanTriggerHists();
  tree_out.Write() ;
  file_out.Close() ;
  Hists.clear();
  Hists.shrink_to_fit();
  cout<<"Hists cleaned"<<endl;
  HistsSysUp.clear();
  HistsSysUp.shrink_to_fit();
  cout<<"HistsSysUp cleaned"<<endl;
  HistsSysDown.clear();
  HistsSysDown.shrink_to_fit();
  cout<<"HistsSysDown cleaned"<<endl;
  HistsJecUp.clear();
  HistsJecUp.shrink_to_fit();
  cout<<"HistsJecUp cleaned"<<endl;
  HistsJecDown.clear();
  HistsJecDown.shrink_to_fit();
  cout<<"HistsJecDown cleaned"<<endl;
  HistsSysReweightsPDF.clear();
  HistsSysReweightsPDF.shrink_to_fit();
  cout<<"HistsSysReweightPDF cleaned"<<endl;
  HistsSysReweightsQscale.clear();
  HistsSysReweightsQscale.shrink_to_fit();
  cout<<"HistsSysReweightsQscale cleaned"<<endl;
  HistsSysReweightsPS.clear();
  HistsSysReweightsPS.shrink_to_fit();
  cout<<"HistsSysReweightPS cleaned"<<endl;
  cout<<"Job is finished"<<endl;

cout<<"Total Virtual Memory: "<<(getValue()-memoryInit)/1000.0<<" MB"<<endl;

}


void MyAnalysis::FillD3Hists(D3HistsContainer H3, int v1, std::vector<int> v2, int v3, float value, std::vector<float> weight, std::vector<WCFit> wcfit){
  for (int i = 0; i < v2.size(); ++i) {
    if(v2[i]>=0) H3[v1][v2[i]][v3]->Fill(value, weight[i], wcfit[i]);
  }
}

void MyAnalysis::FillD4Hists(D4HistsContainer H4, int v1, std::vector<int> v2, int v3, int v4, float value, std::vector<float> weight, std::vector<WCFit> wcfit){
  for (int i = 0; i < v2.size(); ++i) {
    if(v2[i]>=0) H4[v1][v2[i]][v3][v4]->Fill(value, weight[i], wcfit[i]);
  }
}

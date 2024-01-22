#include "triggerEffAnalysis.h"
#include "Utils.h"
#include<array> 

//TrigObj_filterBits : Int_t extra bits of associated information:   *
// 1 = CaloIdL_TrackIdL_IsoVL, 2 = 1e (WPTight), 4 = 1e (WPLoose), 8 = OverlapFilter PFTau, 16 = 2e, 32 = 1e-1mu, 64 = 1e-1tau, 128 = 3e, 256 = 2e-1mu, 512 = 1e-2mu, 1024 = 1e (32_L1DoubleEG_AND_L1SingleEGOr), 2048 = 1e (CaloIdVT_GsfTrkIdT), 4096 = 1e (PFJet), 8192 = 1e (Photon175_OR_Photon200), 16384 = CaloIdLMWPMS2, 32768 = CaloIdLMWPMS2Unseeded, 65536 = CaloIdVTGsfTrkIdTGsfDphi for Electron (PixelMatched e/gamma);
//             bit 0               "filter('*CaloIdLTrackIdLIsoVL*TrackIso*Filter') + " \
//             bit 1               "2*filter('hltEle*WPTight*TrackIsoFilter*') + " \
//             bit 2               "4*filter('hltEle*WPLoose*TrackIsoFilter') + " \
//             bit 3               "8*filter('hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter') + " \
//             bit 4               "16*filter('hltEG*HEFilter') + " \
//             bit 5               "32*filter('hltEG*CaloIdLClusterShapeFilter')  + " \
//             bit 6               "64*filter('hltEle*CaloIdLPixelMatchFilter') + " \
//             bit 7               "128*filter('hltDiEG*HEUnseededFilter') + " \
//             bit 8               "256*filter('hltDiEG*CaloIdLClusterShapeUnseededFilter') + " \
//             bit 9               "512*filter('hltDiEle*CaloIdLPixelMatchUnseededFilter') + " \
//             bit 10              "1024*filter('hltDiEle*CaloIdLGsfTrkIdVLMWPMS2UnseededFilter') + " \
//             bit 11              "2048*filter('hltDiEle*CaloIdLGsfTrkIdVLDPhiUnseededFilter') + " \
//             bit 12              "4096*path('HLT_Ele*PFJet*') + " \
//             bit 13              "8192*max(filter('hltEG175HEFilter'),filter('hltEG200HEFilter')) + " \
//             bit 14              "16384     * filter('hltEle*CaloIdLMWPMS2Filter') + " \
//             bit 15              "32768     * filter('hltDiEle*CaloIdLMWPMS2UnseededFilter') + " \
//             bit 16              "65536     * filter('hltEle*CaloIdVTGsfTrkIdTGsfDphiFilter')"

triggerEffAnalysis::triggerEffAnalysis(MyAnalysis *evt, string Y){
  year=Y;
  eI = evt;
  TriggerHists.resize(vars.size());
  std::stringstream name;
  TH1F *h_test;
  for( auto it = vars.cbegin() ; it != vars.cend() ; ++it ){
    name<<it->first;
    h_test = new TH1F((name.str()).c_str(),(name.str()).c_str(), it->second.first.second-1, it->second.second);
    h_test->StatOverflows(kTRUE);
    h_test->Sumw2(kTRUE);
    TriggerHists[it->second.first.first]=h_test;
    name.str("");
  }

  if(year != "2018") etProbe = 33;
  if(year == "2018") etProbe = 25;
}

int triggerEffAnalysis::WInd(std::map<TString, std::pair<std::pair<int, int>, float*>> V, TString name){
  return V.find(name)->second.first.first;
}

int triggerEffAnalysis::etaRegion(float eta){
    int i=0;
    if(abs(eta)<0.79) i=1;
    if(abs(eta)>0.79 && abs(eta)<1.10) i=2;
    if(abs(eta)>1.10 && abs(eta)<1.44) i=3;
    if(abs(eta)>1.56 && abs(eta)<1.70) i=-1;
    if(abs(eta)>1.70 && abs(eta)<2.10) i=-2; 
    if(abs(eta)>2.10 && abs(eta)<2.5) i=-3;
    return i;
}

bool triggerEffAnalysis::isTag(int l){
  if (!eI->Electron_cutBased_HEEP[l]) return false;
  if (abs(eI->Electron_eta_SC[l]) > 1.4442) return false;
  int ProbeMatchedTO=-1;
  for(int j=0;j<eI->nTrigObj;j++){
    if (eI->TrigObj_id[j]  != 11 ) continue;
    float tag_objDR=dR(eI->Electron_eta_SC[l],  eI->Electron_phi[l], eI->TrigObj_eta[j], eI->TrigObj_phi[j]);
    if(tag_objDR<0.1) ProbeMatchedTO=j ;
  }
  if(ProbeMatchedTO<0) return false;
  if(eI->TrigObj_l1pt[ProbeMatchedTO]<40) return false;
  if(year == "2016postVFP" || year == "2016preVFP"){
    if(eI->TrigObj_pt[ProbeMatchedTO]>35 && abs(eI->TrigObj_id[ProbeMatchedTO])==11  && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 1) & 1)) return true;
  }
  else{
    if(eI->TrigObj_pt[ProbeMatchedTO]>35 && abs(eI->TrigObj_id[ProbeMatchedTO])==11  && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 1) & 1) && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 14) & 1)) return true;
  }
  return false;
}

bool triggerEffAnalysis::isProbe(int l){
  if (!eI->Electron_cutBased_HEEP[l]) return false;
  return true;
}

bool triggerEffAnalysis::isProbeSeeded(int l){
  if(ProbeMatchedTO<0) return false;
  if((eI->TrigObj_filterBits[ProbeMatchedTO] >> 3) & 1) return true;
  return false;
}

bool triggerEffAnalysis::isProbeMatchedWPTight(int l){
  if(ProbeMatchedTO<0) return false;
  if((eI->TrigObj_filterBits[ProbeMatchedTO] >> 1) & 1) return true;
  return false;
}

bool triggerEffAnalysis::isProbeMatchedEt(int l, float Et){
  if(ProbeMatchedTO<0) return false;
  if(eI->TrigObj_pt[ProbeMatchedTO]>Et) return true;
  return false;
}

bool triggerEffAnalysis::isProbeMatchedHE(int l, float Et){
  if(ProbeMatchedTO<0) return false;
  if(eI->TrigObj_pt[ProbeMatchedTO]>Et && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 7) & 1)) return true;
  return false;
}

bool triggerEffAnalysis::isProbeMatchedClusterShape(int l, float Et){
  if(ProbeMatchedTO<0) return false;
  if(eI->TrigObj_pt[ProbeMatchedTO]>Et && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 7) & 1) && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 8) & 1)) return true;
  return false;
}

bool triggerEffAnalysis::isProbeMatchedPixelMatch(int l, float Et){
  if(ProbeMatchedTO<0) return false;
  if((year=="2017" || year=="2018") && eI->TrigObj_pt[ProbeMatchedTO]>Et && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 7) & 1) && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 8) & 1) && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 9) & 1)) return true;
  return false;
}

bool triggerEffAnalysis::isProbeMatchedCaloIdLMWPMS2(int l, float Et){
  if(ProbeMatchedTO<0) return false;
  if((year=="2017" || year=="2018") && eI->TrigObj_pt[ProbeMatchedTO]>Et && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 15) & 1)) return true;
  if((year=="2016postVFP") && eI->TrigObj_pt[ProbeMatchedTO]>Et && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 10) & 1)) return true;
  if((year=="2016preVFP" && eI->run>276453) && eI->TrigObj_pt[ProbeMatchedTO]>Et && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 11) & 1)) return true;
  if((year=="2016preVFP" && eI->run<276453) && eI->TrigObj_pt[ProbeMatchedTO]>Et && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 9) & 1)) return true;
  return false;
}

bool triggerEffAnalysis::isProbeMatchedHESeeded(int l, float Et){
  if(ProbeMatchedTO<0) return false;
  if(eI->TrigObj_pt[ProbeMatchedTO]>Et && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 4) & 1)) return true;
  return false;
}

bool triggerEffAnalysis::isProbeMatchedClusterShapeSeeded(int l, float Et){
  if(ProbeMatchedTO<0) return false;
  if(eI->TrigObj_pt[ProbeMatchedTO]>Et && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 4) & 1) && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 5) & 1)) return true;
  return false;
}

bool triggerEffAnalysis::isProbeMatchedPixelMatchSeeded(int l, float Et){
  if(ProbeMatchedTO<0) return false;
  if((year=="2017" || year=="2018") && eI->TrigObj_pt[ProbeMatchedTO]>Et && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 4) & 1) && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 5) & 1) && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 6) & 1)) return true;
  if(year=="2016preVFP" && eI->run>276453 && eI->TrigObj_pt[ProbeMatchedTO]>Et) return true;
  if((year=="2016preVFP" && eI->run<276453)&& eI->TrigObj_pt[ProbeMatchedTO]>Et && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 9) & 1)) return true;
  return false;
}

bool triggerEffAnalysis::isProbeMatchedCaloIdLMWPMS2Seeded(int l, float Et){
  if(ProbeMatchedTO<0)return false;
  if((year=="2017" || year=="2018") && eI->TrigObj_pt[ProbeMatchedTO]>Et && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 14) & 1)) return true;
  if(year=="2016preVFP" && eI->run>276453 && eI->TrigObj_pt[ProbeMatchedTO]>Et) return true;
  if((year=="2016preVFP" && eI->run<276453) && eI->TrigObj_pt[ProbeMatchedTO]>Et && ((eI->TrigObj_filterBits[ProbeMatchedTO] >> 15) & 1)) return true;
  return false;
}



void triggerEffAnalysis::fillTriggerHists(){
  TLorentzVector t, p;
  if(year=="2016preVFP" || year=="2016postVFP")   tagTriggerPass=eI->HLT_Ele27_WPTight_Gsf;
  else tagTriggerPass=eI->HLT_Ele35_WPTight_Gsf;
  std::vector<int> tagVec;
  std::vector<int> probVec;
  ProbeMatchedTO=-1;
  tagMatchTO=-1;
  int reg=0;
  bool passWPTight;
  bool passEt;
  bool passMW;
  if(!tagTriggerPass) return;
  for(int l=0;l<eI->nElectron;l++){
    if(etaRegion(eI->Electron_eta_SC[l])==0) continue;
    if(isTag(l)) tagVec.push_back(l);
    if(isProbe(l)) probVec.push_back(l);
  }
cout<<probVec.size()<<":"<<tagVec.size()<<endl;
  if(probVec.size()!=2 || tagVec.size()==0) return;
//one tag one probe
  int pl=0;
  int tl=0;
  float etP=0;
  float etT=0;
  if(tagVec.size()==1){
    for(int l=0;l<probVec.size();l++){
      if(probVec[l]==tagVec[0]) tl=l;   
    }
    if(tl==pl) pl=1;
  }
  if(tagVec.size()==2){
    if(eI->event % 2 == 0 && eI->Electron_phi[probVec[1]] < eI->Electron_phi[probVec[0]]) pl=1;
    if(eI->event % 2 != 0 && eI->Electron_phi[probVec[1]] > eI->Electron_phi[probVec[0]]) pl=1;
    if (pl==0) tl=1;
  }
//find the matched trigger objects
  for(int j=0;j<eI->nTrigObj;j++){
    if (eI->TrigObj_id[j]  != 11 ) continue;
    if(dR(eI->Electron_eta_SC[probVec[pl]],  eI->Electron_phi[probVec[pl]], eI->TrigObj_eta[j], eI->TrigObj_phi[j])<0.1) ProbeMatchedTO=j;
    if(dR(eI->Electron_eta_SC[probVec[tl]],  eI->Electron_phi[probVec[tl]], eI->TrigObj_eta[j], eI->TrigObj_phi[j])<0.1) tagMatchTO=j;    
  }
  reg = etaRegion(eI->Electron_eta_SC[probVec[pl]]);
  passWPTight=isProbeMatchedWPTight(probVec[pl]);
  passMW=isProbeMatchedCaloIdLMWPMS2(probVec[pl],etProbe);
  etP=eI->Electron_pt[probVec[pl]]*(eI->Electron_scEtOverPt[probVec[pl]]+1); 
  etT=eI->Electron_pt[probVec[tl]]*(eI->Electron_scEtOverPt[probVec[tl]]+1);
  t.SetPtEtaPhiM(etT, eI->Electron_eta_SC[probVec[tl]], eI->Electron_phi[probVec[tl]], 0.000511) ;
  p.SetPtEtaPhiM(etP, eI->Electron_eta_SC[probVec[pl]], eI->Electron_phi[probVec[pl]], 0.000511) ;

  TriggerHists[WInd(vars,"h_eta_totalEt")]->Fill(eI->Electron_eta_SC[probVec[pl]]);
  TriggerHists[WInd(vars,"h_mass_totalEt")]->Fill((t+p).M());

  TriggerHists[WInd(vars,"h_eta_totalWPTightEt35")]->Fill(eI->Electron_eta_SC[probVec[pl]]);
  TriggerHists[WInd(vars,"h_mass_totalWPTightEt35")]->Fill((t+p).M());
  if(reg>0){
    TriggerHists[WInd(vars,"h_pt_Barrel_totalEt")]->Fill(etP);    
    if(reg==1) TriggerHists[WInd(vars,"h_pteta1_Barrel_totalEt")]->Fill(etP);
    if(reg==2) TriggerHists[WInd(vars,"h_pteta2_Barrel_totalEt")]->Fill(etP);
    if(reg==3) TriggerHists[WInd(vars,"h_pteta3_Barrel_totalEt")]->Fill(etP);

    TriggerHists[WInd(vars,"h_pt_Barrel_totalWPTightEt35")]->Fill(etP);
    if(reg==1) TriggerHists[WInd(vars,"h_pteta1_Barrel_totalWPTightEt35")]->Fill(etP);
    if(reg==2) TriggerHists[WInd(vars,"h_pteta2_Barrel_totalWPTightEt35")]->Fill(etP);
    if(reg==3) TriggerHists[WInd(vars,"h_pteta3_Barrel_totalWPTightEt35")]->Fill(etP);
  }
  if(reg<0){
    TriggerHists[WInd(vars,"h_pt_Endcaps_totalEt")]->Fill(etP);
    if(reg==-1) TriggerHists[WInd(vars,"h_pteta1_Endcaps_totalEt")]->Fill(etP);
    if(reg==-2) TriggerHists[WInd(vars,"h_pteta2_Endcaps_totalEt")]->Fill(etP);
    if(reg==-3) TriggerHists[WInd(vars,"h_pteta3_Endcaps_totalEt")]->Fill(etP);

    TriggerHists[WInd(vars,"h_pt_Endcaps_totalWPTightEt35")]->Fill(etP);
    if(reg==-1) TriggerHists[WInd(vars,"h_pteta1_Endcaps_totalWPTightEt35")]->Fill(etP);
    if(reg==-2) TriggerHists[WInd(vars,"h_pteta2_Endcaps_totalWPTightEt35")]->Fill(etP);
    if(reg==-3) TriggerHists[WInd(vars,"h_pteta3_Endcaps_totalWPTightEt35")]->Fill(etP);
  }
  //Et turn on
  if(isProbeMatchedEt(probVec[pl],  etProbe)){
    TriggerHists[WInd(vars,"h_eta_totalCaloIdLMWPMS2")]->Fill(eI->Electron_eta_SC[probVec[pl]]);
    TriggerHists[WInd(vars,"h_mass_totalCaloIdLMWPMS2")]->Fill((t+p).M());
    TriggerHists[WInd(vars,"h_eta_passEt")]->Fill(eI->Electron_eta_SC[probVec[pl]]);
    TriggerHists[WInd(vars,"h_mass_passEt")]->Fill((t+p).M());
    if(reg>0){
      TriggerHists[WInd(vars,"h_pt_Barrel_totalCaloIdLMWPMS2")]->Fill(etP);
      TriggerHists[WInd(vars,"h_mass_Barrel_totalCaloIdLMWPMS2")]->Fill((t+p).M());
      if(reg==1) TriggerHists[WInd(vars,"h_pteta1_Barrel_totalCaloIdLMWPMS2")]->Fill(etP);
      if(reg==2) TriggerHists[WInd(vars,"h_pteta2_Barrel_totalCaloIdLMWPMS2")]->Fill(etP);
      if(reg==3) TriggerHists[WInd(vars,"h_pteta3_Barrel_totalCaloIdLMWPMS2")]->Fill(etP);

      TriggerHists[WInd(vars,"h_pt_Barrel_passEt")]->Fill(etP);
      if(reg==1) TriggerHists[WInd(vars,"h_pteta1_Barrel_passEt")]->Fill(etP);
      if(reg==2) TriggerHists[WInd(vars,"h_pteta2_Barrel_passEt")]->Fill(etP);
      if(reg==3) TriggerHists[WInd(vars,"h_pteta3_Barrel_passEt")]->Fill(etP);
    }
    if(reg<0){
      TriggerHists[WInd(vars,"h_pt_Endcaps_totalCaloIdLMWPMS2")]->Fill(etP);
      TriggerHists[WInd(vars,"h_mass_Endcaps_totalCaloIdLMWPMS2")]->Fill((t+p).M());
      if(reg==-1) TriggerHists[WInd(vars,"h_pteta1_Endcaps_totalCaloIdLMWPMS2")]->Fill(etP);
      if(reg==-2) TriggerHists[WInd(vars,"h_pteta2_Endcaps_totalCaloIdLMWPMS2")]->Fill(etP);
      if(reg==-3) TriggerHists[WInd(vars,"h_pteta3_Endcaps_totalCaloIdLMWPMS2")]->Fill(etP);

      TriggerHists[WInd(vars,"h_pt_Endcaps_passEt")]->Fill(etP);
      if(reg==-1) TriggerHists[WInd(vars,"h_pteta1_Endcaps_passEt")]->Fill(etP);
      if(reg==-2) TriggerHists[WInd(vars,"h_pteta2_Endcaps_passEt")]->Fill(etP);
      if(reg==-3) TriggerHists[WInd(vars,"h_pteta3_Endcaps_passEt")]->Fill(etP);
    }
  }
  //WP tight + ET
  if(isProbeMatchedEt(probVec[pl],35) && passWPTight){
    TriggerHists[WInd(vars,"h_eta_passWPTightEt35")]->Fill(eI->Electron_eta_SC[probVec[pl]]);
    TriggerHists[WInd(vars,"h_mass_passWPTightEt35")]->Fill((t+p).M());
    if(reg>0){
      TriggerHists[WInd(vars,"h_pt_Barrel_passWPTightEt35")]->Fill(etP);
      if(reg==1) TriggerHists[WInd(vars,"h_pteta1_Barrel_passWPTightEt35")]->Fill(etP);
      if(reg==2) TriggerHists[WInd(vars,"h_pteta2_Barrel_passWPTightEt35")]->Fill(etP);
      if(reg==3) TriggerHists[WInd(vars,"h_pteta3_Barrel_passWPTightEt35")]->Fill(etP);
    }
    if(reg<0){
      TriggerHists[WInd(vars,"h_pt_Endcaps_passWPTightEt35")]->Fill(etP);
      if(reg==-1) TriggerHists[WInd(vars,"h_pteta1_Endcaps_passWPTightEt35")]->Fill(etP);
      if(reg==-2) TriggerHists[WInd(vars,"h_pteta2_Endcaps_passWPTightEt35")]->Fill(etP);
      if(reg==-3) TriggerHists[WInd(vars,"h_pteta3_Endcaps_passWPTightEt35")]->Fill(etP);
    }
  }
  //CaloIdLMWPMS2
  if(isProbeMatchedEt(probVec[pl],etProbe)){
    if(reg>0){
      TriggerHists[WInd(vars,"h_pt_Barrel_totalHE")]->Fill(etP);
      TriggerHists[WInd(vars,"h_pt_Barrel_totalClusterShape")]->Fill(etP);
      TriggerHists[WInd(vars,"h_pt_Barrel_totalPixelMatch")]->Fill(etP);
    }
    if(reg<0){
      TriggerHists[WInd(vars,"h_pt_Endcaps_totalHE")]->Fill(etP);
      TriggerHists[WInd(vars,"h_pt_Endcaps_totalClusterShape")]->Fill(etP);
      TriggerHists[WInd(vars,"h_pt_Endcaps_totalPixelMatch")]->Fill(etP);
    }
  }
  if(isProbeMatchedHE(probVec[pl],etProbe)){
    if(reg>0){
      TriggerHists[WInd(vars,"h_pt_Barrel_passHE")]->Fill(etP);
    }
    if(reg<0){
      TriggerHists[WInd(vars,"h_pt_Endcaps_passHE")]->Fill(etP);
    }
  }
  if(isProbeMatchedClusterShape(probVec[pl],etProbe)){
    if(reg>0){
      TriggerHists[WInd(vars,"h_pt_Barrel_passClusterShape")]->Fill(etP);
    }
    if(reg<0){
      TriggerHists[WInd(vars,"h_pt_Endcaps_passClusterShape")]->Fill(etP);
    }
  }
  if(isProbeMatchedPixelMatch(probVec[pl],etProbe)){
    if(reg>0){
      TriggerHists[WInd(vars,"h_pt_Barrel_passPixelMatch")]->Fill(etP);
    }
    if(reg<0){
      TriggerHists[WInd(vars,"h_pt_Endcaps_passPixelMatch")]->Fill(etP);
    }
  }

  if(isProbeMatchedEt(probVec[pl],etProbe) && passMW){
    TriggerHists[WInd(vars,"h_eta_passCaloIdLMWPMS2")]->Fill(eI->Electron_eta_SC[probVec[pl]]);
    TriggerHists[WInd(vars,"h_mass_passCaloIdLMWPMS2")]->Fill((t+p).M());
    if(reg>0){
      TriggerHists[WInd(vars,"h_pt_Barrel_passCaloIdLMWPMS2")]->Fill(etP);
      TriggerHists[WInd(vars,"h_mass_Barrel_passCaloIdLMWPMS2")]->Fill((t+p).M());
      if(reg==1) TriggerHists[WInd(vars,"h_pteta1_Barrel_passCaloIdLMWPMS2")]->Fill(etP);
      if(reg==2) TriggerHists[WInd(vars,"h_pteta2_Barrel_passCaloIdLMWPMS2")]->Fill(etP);
      if(reg==3) TriggerHists[WInd(vars,"h_pteta3_Barrel_passCaloIdLMWPMS2")]->Fill(etP);
    }
    if(reg<0){
      TriggerHists[WInd(vars,"h_pt_Endcaps_passCaloIdLMWPMS2")]->Fill(etP);
      TriggerHists[WInd(vars,"h_mass_Endcaps_passCaloIdLMWPMS2")]->Fill((t+p).M());
      if(reg==-1) TriggerHists[WInd(vars,"h_pteta1_Endcaps_passCaloIdLMWPMS2")]->Fill(etP);
      if(reg==-2) TriggerHists[WInd(vars,"h_pteta2_Endcaps_passCaloIdLMWPMS2")]->Fill(etP);
      if(reg==-3) TriggerHists[WInd(vars,"h_pteta3_Endcaps_passCaloIdLMWPMS2")]->Fill(etP);
    }
  }

//seeded eff
  if(isProbeSeeded(probVec[pl])){
    if(isProbeMatchedEt(probVec[pl],etProbe) ){
      if(reg>0){
        TriggerHists[WInd(vars,"h_pt_Barrel_totalHESeeded")]->Fill(etP);
        TriggerHists[WInd(vars,"h_pt_Barrel_totalClusterShapeSeeded")]->Fill(etP);
        TriggerHists[WInd(vars,"h_pt_Barrel_totalPixelMatchSeeded")]->Fill(etP);
        TriggerHists[WInd(vars,"h_pt_Barrel_totalCaloIdLMWPMS2Seeded")]->Fill(etP);
      }
      if(reg<0){
        TriggerHists[WInd(vars,"h_pt_Endcaps_totalHESeeded")]->Fill(etP);
        TriggerHists[WInd(vars,"h_pt_Endcaps_totalClusterShapeSeeded")]->Fill(etP);
        TriggerHists[WInd(vars,"h_pt_Endcaps_totalPixelMatchSeeded")]->Fill(etP);
        TriggerHists[WInd(vars,"h_pt_Endcaps_totalCaloIdLMWPMS2Seeded")]->Fill(etP);

      }
    }
    if(isProbeMatchedHESeeded(probVec[pl],etProbe)){
      if(reg>0){
        TriggerHists[WInd(vars,"h_pt_Barrel_passHESeeded")]->Fill(etP);
      }
      if(reg<0){
        TriggerHists[WInd(vars,"h_pt_Endcaps_passHESeeded")]->Fill(etP);
      }
    }
    if(isProbeMatchedClusterShapeSeeded(probVec[pl],etProbe)){
      if(reg>0){
        TriggerHists[WInd(vars,"h_pt_Barrel_passClusterShapeSeeded")]->Fill(etP);
      }
      if(reg<0){
        TriggerHists[WInd(vars,"h_pt_Endcaps_passClusterShapeSeeded")]->Fill(etP);
      }
    }
    if(isProbeMatchedPixelMatchSeeded(probVec[pl],etProbe)){
      if(reg>0){
        TriggerHists[WInd(vars,"h_pt_Barrel_passPixelMatchSeeded")]->Fill(etP);
      }
      if(reg<0){
        TriggerHists[WInd(vars,"h_pt_Endcaps_passPixelMatchSeeded")]->Fill(etP);
      }
    }
    if(isProbeMatchedCaloIdLMWPMS2Seeded(probVec[pl],etProbe)){
      if(reg>0){
        TriggerHists[WInd(vars,"h_pt_Barrel_passCaloIdLMWPMS2Seeded")]->Fill(etP);
      }
      if(reg<0){
        TriggerHists[WInd(vars,"h_pt_Endcaps_passCaloIdLMWPMS2Seeded")]->Fill(etP);
      }
    }
  }
}

void triggerEffAnalysis::writeTriggerHists(){
  cout<<"The following Trigger level distributions are written on the output file"<<endl;
  for (int i=0;i<TriggerHists.size();++i){
    cout<<i<<":"<<TriggerHists[i]->GetName()<<":"<<TriggerHists[i]->Integral()<<endl;
    TriggerHists[i]->Write("",TObject::kOverwrite);
  }
}

void triggerEffAnalysis::cleanTriggerHists(){
  for (int i=0;i<TriggerHists.size();++i){
    delete TriggerHists[i];
  }
  TriggerHists.clear();
  TriggerHists.shrink_to_fit();
} 

triggerEffAnalysis::~triggerEffAnalysis(){}

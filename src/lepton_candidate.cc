#include "lepton_candidate.h"

lepton_candidate::lepton_candidate(float pt_in, float eta_in, float phi_in, int charge_in, int ind_in, int lep_in ){
  pt_ = pt_in;
  eta_ = eta_in;
  phi_ = phi_in;
  charge_ = charge_in;
  lep_ = lep_in;
  if(lep_in == 1)  p4_.SetPtEtaPhiM(pt_, eta_, phi_, 0.000511) ;
  if(lep_in == 10)  p4_.SetPtEtaPhiM(pt_, eta_, phi_, 0.10566) ;
  indice_ = ind_in;
}

float lepton_candidate::eleFakeRate(TString year){
  float fr;
  if(abs(eta_)<1.444){
    if(35<=pt_< 131.6) fr= 0.14 - 0.0029*pt_ + 2.56*10e-5*pow(pt_,2) - 8.48*10e-8*pow(pt_,3);
    if(131.6 <=pt_< 359.3) fr= 0.02 - 0.00013*pt_ + 3.5*10e-7*pow(pt_,2) - 2.9*10e-10*pow(pt_,3);
    if(pt_>=359.3) fr= 0.00514 + 4.73*10e-7*pt_;
  }
  else if(1.444<abs(eta_)<2.0){
    if(35<= pt_< 125) fr= 0.1012 - 0.00094*pt_ + 3.37*10e-6*pow(pt_,2);
    if(125<=pt_< 226.3) fr= 0.0488 - 11.37*10e-5*pt_;
    if(pt_>=226.3) fr= 0.0241 - 1.24*10e-6*pt_;
  }
  else{
    if(35<=pt_<152) fr= 0.0622 - 0.00012*pt_;
    if(pt_>= 152) fr= 0.0387;
  }
  return fr;
}


lepton_candidate::~lepton_candidate(){}



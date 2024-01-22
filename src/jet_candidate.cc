#include "jet_candidate.h"

jet_candidate::jet_candidate(float pt_in, float eta_in, float phi_in, float M_in, float btag_in, TString year, int ind_in){
  pt_ = pt_in;
  eta_ = eta_in;
  phi_ = phi_in;
  btag_ = isb(btag_in,year);
  p4_.SetPtEtaPhiM(pt_, eta_, phi_, M_in) ;
  flavor_ = ind_in;
  mass_ = M_in;
}


int jet_candidate::isb(float btag_in , TString year){
  int R = 0;
  if (year == "2016preVFP" && btag_in > 0.25980) R=1;
  if (year == "2016postVFP" && btag_in > 0.2489) R=1;
  if (year == "2017" && btag_in > 0.3040) R=1;
  if (year == "2018" && btag_in > 0.2783) R=1;
  return R;
}
  
jet_candidate::~jet_candidate(){}



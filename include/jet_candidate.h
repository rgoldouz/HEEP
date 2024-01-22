#ifndef MY_jet_candidate
#define MY_jet_candidate

#include<cmath>
#include<string>
#include<iostream>
#include<vector>
#include<complex>
#include <TLorentzVector.h>

using namespace std;
//using namespace math;
class jet_candidate {
  
public:
  jet_candidate(float, float, float, float, float, TString, int );
  virtual ~jet_candidate();
  float pt_;
  float eta_;
  float phi_;
  float mass_;
  int flavor_;
  int btag_;
  int isb(float, TString);
  TLorentzVector p4_;


private:
  
};

#endif


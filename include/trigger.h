#ifndef MY_trigger
#define MY_trigger

#include<cmath>
#include<string>
#include<iostream>
#include<vector>
#include<complex>
#include <TLorentzVector.h>

using namespace std;
//using namespace math;
class trigger {
  
public:
  trigger(float, float, float, int, int );
  ~trigger();
  void evalIsHeep(float, float, float, float, float, float, float, float, float, float, float, float, float);
  void evalIsHeep_tkiso_removed(float, float, float, float, float, float, float, float, float, float, float, float, float);
  float et_;
  float eta_;
  float phi_;
  int charge_;
  int region_;
  int indice_;

  int truthMatched_;
  int passTrigger_;
  int isHeep_;
  int isHeep_tkiso_removed_;
  TLorentzVector p4_;

private:
  const float mEle = 0.000511;
};

#endif


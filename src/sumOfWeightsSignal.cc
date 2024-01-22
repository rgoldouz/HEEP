#include "sumOfWeightsSignal.h"

sumOfWeightsSignal::sumOfWeightsSignal(){};
sumOfWeightsSignal::~sumOfWeightsSignal(){};

std::vector<float> sumOfWeightsSignal::LHEWeightSignal(TString Name){
  if (Name == "2016_LFVStVecC") return sumOfLHEWeight_2016_LFVStVecC;
  if (Name == "2016_LFVStVecU") return sumOfLHEWeight_2016_LFVStVecU;
  if (Name == "2016_LFVTtVecC") return sumOfLHEWeight_2016_LFVTtVecC;
  if (Name == "2016_LFVTtVecU") return sumOfLHEWeight_2016_LFVTtVecU;
  if (Name == "2016_LFVStScalarC") return sumOfLHEWeight_2016_LFVStScalarC;
  if (Name == "2016_LFVStScalarU") return sumOfLHEWeight_2016_LFVStScalarU;
  if (Name == "2016_LFVTtScalarC") return sumOfLHEWeight_2016_LFVTtScalarC;
  if (Name == "2016_LFVTtScalarU") return sumOfLHEWeight_2016_LFVTtScalarU;
  if (Name == "2016_LFVStTensorC") return sumOfLHEWeight_2016_LFVStTensorC;
  if (Name == "2016_LFVStTensorU") return sumOfLHEWeight_2016_LFVStTensorU;
  if (Name == "2016_LFVTtTensorC") return sumOfLHEWeight_2016_LFVTtTensorC;
  if (Name == "2016_LFVTtTensorU") return sumOfLHEWeight_2016_LFVTtTensorU;

  if (Name == "2017_LFVStVecC") return sumOfLHEWeight_2017_LFVStVecC;
  if (Name == "2017_LFVStVecU") return sumOfLHEWeight_2017_LFVStVecU;
  if (Name == "2017_LFVTtVecC") return sumOfLHEWeight_2017_LFVTtVecC;
  if (Name == "2017_LFVTtVecU") return sumOfLHEWeight_2017_LFVTtVecU;
  if (Name == "2017_LFVStScalarC") return sumOfLHEWeight_2017_LFVStScalarC;
  if (Name == "2017_LFVStScalarU") return sumOfLHEWeight_2017_LFVStScalarU;
  if (Name == "2017_LFVTtScalarC") return sumOfLHEWeight_2017_LFVTtScalarC;
  if (Name == "2017_LFVTtScalarU") return sumOfLHEWeight_2017_LFVTtScalarU;
  if (Name == "2017_LFVStTensorC") return sumOfLHEWeight_2017_LFVStTensorC;
  if (Name == "2017_LFVStTensorU") return sumOfLHEWeight_2017_LFVStTensorU;
  if (Name == "2017_LFVTtTensorC") return sumOfLHEWeight_2017_LFVTtTensorC;
  if (Name == "2017_LFVTtTensorU") return sumOfLHEWeight_2017_LFVTtTensorU;

  if (Name == "2018_LFVStVecC") return sumOfLHEWeight_2018_LFVStVecC;
  if (Name == "2018_LFVStVecU") return sumOfLHEWeight_2018_LFVStVecU;
  if (Name == "2018_LFVTtVecC") return sumOfLHEWeight_2018_LFVTtVecC;
  if (Name == "2018_LFVTtVecU") return sumOfLHEWeight_2018_LFVTtVecU;
  if (Name == "2018_LFVStScalarC") return sumOfLHEWeight_2018_LFVStScalarC;
  if (Name == "2018_LFVStScalarU") return sumOfLHEWeight_2018_LFVStScalarU;
  if (Name == "2018_LFVTtScalarC") return sumOfLHEWeight_2018_LFVTtScalarC;
  if (Name == "2018_LFVTtScalarU") return sumOfLHEWeight_2018_LFVTtScalarU;
  if (Name == "2018_LFVStTensorC") return sumOfLHEWeight_2018_LFVStTensorC;
  if (Name == "2018_LFVStTensorU") return sumOfLHEWeight_2018_LFVStTensorU;
  if (Name == "2018_LFVTtTensorC") return sumOfLHEWeight_2018_LFVTtTensorC;
  if (Name == "2018_LFVTtTensorU") return sumOfLHEWeight_2018_LFVTtTensorU;
}

std::vector<float> sumOfWeightsSignal::GenWeightSignal(TString Name){
  if (Name == "2016_LFVStVecC") return sumOfMcWeight_2016_LFVStVecC;
  if (Name == "2016_LFVStVecU") return sumOfMcWeight_2016_LFVStVecU;
  if (Name == "2016_LFVTtVecC") return sumOfMcWeight_2016_LFVTtVecC;
  if (Name == "2016_LFVTtVecU") return sumOfMcWeight_2016_LFVTtVecU;
  if (Name == "2016_LFVStScalarC") return sumOfMcWeight_2016_LFVStScalarC;
  if (Name == "2016_LFVStScalarU") return sumOfMcWeight_2016_LFVStScalarU;
  if (Name == "2016_LFVTtScalarC") return sumOfMcWeight_2016_LFVTtScalarC;
  if (Name == "2016_LFVTtScalarU") return sumOfMcWeight_2016_LFVTtScalarU;
  if (Name == "2016_LFVStTensorC") return sumOfMcWeight_2016_LFVStTensorC;
  if (Name == "2016_LFVStTensorU") return sumOfMcWeight_2016_LFVStTensorU;
  if (Name == "2016_LFVTtTensorC") return sumOfMcWeight_2016_LFVTtTensorC;
  if (Name == "2016_LFVTtTensorU") return sumOfMcWeight_2016_LFVTtTensorU;

  if (Name == "2017_LFVStVecC") return sumOfMcWeight_2017_LFVStVecC;
  if (Name == "2017_LFVStVecU") return sumOfMcWeight_2017_LFVStVecU;
  if (Name == "2017_LFVTtVecC") return sumOfMcWeight_2017_LFVTtVecC;
  if (Name == "2017_LFVTtVecU") return sumOfMcWeight_2017_LFVTtVecU;
  if (Name == "2017_LFVStScalarC") return sumOfMcWeight_2017_LFVStScalarC;
  if (Name == "2017_LFVStScalarU") return sumOfMcWeight_2017_LFVStScalarU;
  if (Name == "2017_LFVTtScalarC") return sumOfMcWeight_2017_LFVTtScalarC;
  if (Name == "2017_LFVTtScalarU") return sumOfMcWeight_2017_LFVTtScalarU;
  if (Name == "2017_LFVStTensorC") return sumOfMcWeight_2017_LFVStTensorC;
  if (Name == "2017_LFVStTensorU") return sumOfMcWeight_2017_LFVStTensorU;
  if (Name == "2017_LFVTtTensorC") return sumOfMcWeight_2017_LFVTtTensorC;
  if (Name == "2017_LFVTtTensorU") return sumOfMcWeight_2017_LFVTtTensorU;

  if (Name == "2018_LFVStVecC") return sumOfMcWeight_2018_LFVStVecC;
  if (Name == "2018_LFVStVecU") return sumOfMcWeight_2018_LFVStVecU;
  if (Name == "2018_LFVTtVecC") return sumOfMcWeight_2018_LFVTtVecC;
  if (Name == "2018_LFVTtVecU") return sumOfMcWeight_2018_LFVTtVecU;
  if (Name == "2018_LFVStScalarC") return sumOfMcWeight_2018_LFVStScalarC;
  if (Name == "2018_LFVStScalarU") return sumOfMcWeight_2018_LFVStScalarU;
  if (Name == "2018_LFVTtScalarC") return sumOfMcWeight_2018_LFVTtScalarC;
  if (Name == "2018_LFVTtScalarU") return sumOfMcWeight_2018_LFVTtScalarU;
  if (Name == "2018_LFVStTensorC") return sumOfMcWeight_2018_LFVStTensorC;
  if (Name == "2018_LFVStTensorU") return sumOfMcWeight_2018_LFVStTensorU;
  if (Name == "2018_LFVTtTensorC") return sumOfMcWeight_2018_LFVTtTensorC;
  if (Name == "2018_LFVTtTensorU") return sumOfMcWeight_2018_LFVTtTensorU;
}



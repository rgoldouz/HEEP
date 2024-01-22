#include "../include/MyAnalysis.h"
#include <chrono>
using namespace std::chrono;
int main(){
auto start = high_resolution_clock::now();
    TChain* ch    = new TChain("Events") ;
//    ch ->Add("/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL17/v2/UL17_DY50/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_UL17_DY50/221008_092819/0000/tree_1.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/reNanoAODUL18/UL18/data_UL18_D_EGamma/EGamma/crab_data_UL18_D_EGamma/230430_194443/0000/NANO_NANO_114.root");
//    ch ->Add("/afs/crc.nd.edu/user/r/rgoldouz/HEEP/UL_v9/CMSSW_10_6_27/src/nanoTest/NANO_NANO.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/reNanoAODUL18/UL18/data_UL18_D_EGamma/EGamma/crab_data_UL18_D_EGamma/230502_103045/0000/NANO_NANO_3-119.root");
//    ch ->Add("/afs/crc.nd.edu/user/r/rgoldouz/HEEP/UL_v9/CMSSW_10_6_27/src/nanoTest/NANO_NANO.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/reNanoAODUL18/UL18/data_UL18_A_EGamma/EGamma/crab_data_UL18_A_EGamma/230502_102927/0000/NANO_NANO_1-134.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/reNanoAODUL18/UL18/data_UL18_A_EGamma/EGamma/crab_data_UL18_A_EGamma/230502_102927/0000/NANO_NANO_1-161.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/reNanoAODUL18/UL18/data_UL18_A_EGamma/EGamma/crab_data_UL18_A_EGamma/230502_102927/0001/NANO_NANO_1134.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/reNanoAODUL18/UL18/data_UL18_A_EGamma/EGamma/crab_data_UL18_A_EGamma/230502_102927/0000/NANO_NANO_148.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/reNanoAODUL18/UL18/data_UL18_A_EGamma/EGamma/crab_data_UL18_A_EGamma/230502_102927/0000/NANO_NANO_195.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/reNanoAODUL18/UL18/data_UL18_A_EGamma/EGamma/crab_data_UL18_A_EGamma/230502_102927/0000/NANO_NANO_2-58.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/reNanoAODUL18/UL18/data_UL18_A_EGamma/EGamma/crab_data_UL18_A_EGamma/230502_102927/0000/NANO_NANO_2-96.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/reNanoAODUL18/UL18/data_UL18_A_EGamma/EGamma/crab_data_UL18_A_EGamma/230502_102927/0000/NANO_NANO_391.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/reNanoAODUL18/UL18/data_UL18_A_EGamma/EGamma/crab_data_UL18_A_EGamma/230502_102927/0000/NANO_NANO_676.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/reNanoAODUL18/UL18/data_UL18_A_EGamma/EGamma/crab_data_UL18_A_EGamma/230502_102927/0000/NANO_NANO_722.root");
//   ch ->Add("/hadoop/store/user/rgoldouz/reNanoAODUL17/UL17/data_UL17_B_SingleElectron/SingleElectron/crab_data_UL17_B_SingleElectron/231003_145205/0000/NANO_NANO_46.root");
//   ch ->Add("/hadoop/store/user/rgoldouz/reNanoAODUL17/UL17/testLOB/DoubleEG/crab_testLOB/231010_103216/0000/tree*");
   ch ->Add("/hadoop/store/user/rgoldouz/reNanoAODUL_HEEP/UL16preVFP/data_UL16preVFP_C_SingleElectron/SingleElectron/crab_data_UL16preVFP_C_SingleElectron/231015_105958/0000/tree_385.root");
    MyAnalysis * t1 = new MyAnalysis(ch);
//    t1->Loop("UL17_DY2000to3000", "mc","none", "2017", "none", 1, 41.48, 494000, 0, 38,t1);
    t1->Loop("data_UL17_B_SingleElectron", "data", "SingleElectron", "2016preVFP", "B", 87.31, 19.52, 32368794.5651, 0, 1,t1);
    delete t1;
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout << "time for running the code in second:"<<duration.count() << endl;
}

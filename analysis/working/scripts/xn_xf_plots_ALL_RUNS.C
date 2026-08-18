#include "xn_xf_plots_per_run.C"
#include "TROOT.h"

const int nDet = 24;

void xn_xf_plots_ALL_RUNS(){
    gROOT->SetBatch(kTRUE);
    for (int i = 0; i < nDet; ++i){
        xn_xf_plots_per_run(i);
    }
}

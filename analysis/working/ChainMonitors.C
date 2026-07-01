#include "Monitors.C+" // the plus sign mean compilation

void ChainMonitors(int RUNNUM = -1, int RUNNUM2 = -1, bool saveCanvas = false, bool isTraceON = false)
{
    ///default saveCanvas = false, no save Cavas
    ///                   = true, save Canvas
    TChain* chain = new TChain("gen_tree");

    if ( RUNNUM == -1){
        /// this list only for manual Chain sort
        ///********** start Marker for AutoCalibration.

        // Medium target 114 - runs 9-13
        chain->Add("../root_data/gen_run009.root");
        chain->Add("../root_data/gen_run01[0-3].root");

        // Thin target 88 - runs 14-17, [fix DAQ], 20-23
        chain->Add("../root_data/gen_run01[4-7].root");
        chain->Add("../root_data/gen_run02[0-3].root");

        // Junk runs 18-19

        // "Thick" target 119 24-25
        chain->Add("../root_data/gen_run02[4-5].root");

        // Background run 26

        ///********** end Marker for AutoCalibration.


    }
    else{
        // No manual chain sorting
        TString fileName;
        int endRUNNUM = RUNNUM2;
        if( RUNNUM2 == -1) endRUNNUM = RUNNUM;

        for( int i = RUNNUM ; i <= endRUNNUM ; i++){
            if( isTraceON == false )fileName.Form("../root_data/gen_run%03d.root", i);
            if( isTraceON == true ) fileName.Form("../root_data/trace_run%03d.root", i);
            chain->Add(fileName);
        }
    }

    // Print useful information regarding the sorting
    printf("\033[0;31m========================================== Number of Files : %2d\n",chain->GetListOfFiles()->GetEntries());
    chain->GetListOfFiles()->Print();
    printf("========================================== Number of Files : %2d\033[0m\n",chain->GetListOfFiles()->GetEntries());
    printf(" number of entries : %llu \n", chain->GetEntries());

    //Simple call TSelector
    //chain->Process("Monitors.C+");

    //Some input to TSelector
    Monitors* selector = new Monitors();
    if( saveCanvas ) selector->printControl(0); //quit after terminated
    chain->Process(selector, "");
    return;
}

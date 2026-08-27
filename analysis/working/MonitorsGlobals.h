#ifndef MONITORSGLOBALS_H
#define MONITORSGLOBALS_H

class TCanvas;
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TLatex.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TTree.h"

#include <vector>

#define NROW 4 // number of side of array
#define NCOL  (NARRAY / NROW) // 12 detectors, 3 rows
namespace monitorsglobals{
    inline const int rawEnergyRange[2] = {     100,    3000};       /// share with e, ring, xf, xn
    inline const int    energyRange[2] = {     0,     14};       /// in the E-Z plot
    inline const int     rdtDERange[2] = {     0,    1500};
    inline const int      rdtERange[2] = {     0,    1500};
    inline const int    apolloRange[2] = {     0,    1000};
    inline const int      crdtRange[2] = {     0,    8000};
    inline const int      elumRange[2] = {     0,    16000};
    inline const int       TACRange[3] = { 300,   2000,   6000};  /// #bin, min, max
    inline const int      TAC2Range[3] = { 100,    400,    500};
    inline const int   thetaCMRange[2] = {0, 80};

    inline const double     exRange[3] = {  40,    -2,     10};  /// bin [keV], low[MeV], high[MeV]

    inline const int  coinTimeRange[2] = { -200, 200};
    inline const int  timeRangeUser[2] = {0, 99999999}; /// min, use when cannot find time, this set the min and max

    inline const int  icRange [3] = {1000, 1000, 500}; /// max of IC0,1,2

    inline const bool isUseArrayTrace = false;
    inline const bool isUseRDTTrace = true;

    //---Gate
    inline const bool isTimeGateOn     = true;
    inline const int timeGate[2]       = {-30, 20};             /// min, max, 1 ch = 10 ns
    inline const double eCalCut[2]     = {0.5, 50};             /// lower & higher limit for eCal
    inline const bool  isTACGate       = false;
    inline const int tacGate[2]        = {-8000, -2000};
    inline const int dEgate[2]         = {  500,  1500};
    inline const int Eresgate[2]       = { 1000,  4000};
    inline const double thetaCMGate    = 10;                    /// deg
    inline const double xGate          = 0.8;                  ///cut out the edge
    inline const std::vector<int> skipDetID = {2, 11, 20, 21};

    inline TString rdtCutFile1 = "rdtCuts.root";//"test_gate_17N.root";
    inline TString rdtCutFile2 = "";//"test2_gate_17N.root";
    inline const TString ezCutFile   = "";//"ezCut.root";

    //TODO switches for histograms on/off
    //############################################ end of user setting

    inline ULong64_t NumEntries = 0;
    inline ULong64_t ProcessedEntries = 0;
    inline Float_t Frac = 0.1; ///Progress bar
    inline TStopwatch StpWatch;

    //======= Canvas
    inline TCanvas *cCanvas;
    inline TString canvasTitle;
    inline int lastRunID;
    inline bool contFlag;
    inline double runTime=0;

    //======= Recoil Cut
    inline TCutG* cutG; //!  //general temeprary pointer to cut

    inline TObjArray * cutList1;
    inline Bool_t isCutFileOpen1;
    inline int numCut1;

    inline TObjArray * cutList2;
    inline Bool_t isCutFileOpen2;
    inline int numCut2;

    //======= Other Cuts
    inline TCutG* EZCut;
    inline Bool_t isEZCutFileOpen;

    /******************************************************************
     *   variable and histogram naming rules                           *
     *   name are case sensitive, so as any C/C++ code                 *
     *                                                                 *
     *   ID is dettector ID                                            *
     *                                                                 *
     *   raw data from gen_tree are e, xf, xn, ring.                   *
     *   the x from raw data is x                                      *
     *                                                                 *
     *   xf + xn = xs, s for sum                                       *
     *                                                                 *
     *   calibrated data are  eCal, xfCal, xnCal, ringCal.             *
     *   the x from cal data is xCal                                   *
     *                                                                 *
     *   xfCal + xnCal = xsCal                                         *
     *                                                                 *
     *   since the z is always from xCal, so it calls z.               *
     *                                                                 *
     *   Excitation energy calls Ex                                    *
     *                                                                 *
     *                                                                 *
     *   TH2D is always using "V" to seperate 2 variables, like eVx    *
     *                                                                 *
     *  histogram with TCutG, add suffix "GC" for Graphical-Cut.       *
     *                                                                 *
     *******************************************************************/
    //======== raw data
    inline TH1F* he[NARRAY];
    inline TH1F* hring[NARRAY];
    inline TH1F* hxf[NARRAY];
    inline TH1F* hxn[NARRAY];

    inline TH2F* heVring[NARRAY];
    inline TH2F* hxfVxn[NARRAY];
    inline TH2F* heVxs[NARRAY];

    //TH1F* hMultiHit; //TODO

    inline TH2F* heVID;
    inline TH2F* heVIDG; //gated
    inline TH2F* hringVID;
    inline TH2F* hxfVID;
    inline TH2F* hxnVID;

    inline TH2F* heVx[NARRAY]; // e vs (xf-xn)/e
    inline TH2F* hringVx[NARRAY]; // ring vs (xf-xn)/e

    //====== cal data
    inline TH2F* heVxsCal[NARRAY]; // raw e vs xf

    inline TH2F* heCalVxCal[NARRAY]; // eCal vs xCal
    inline TH2F* heCalVxCalG[NARRAY]; // eCal vs xCal
    inline TH1F* heCal[NARRAY];
    inline TH2F* heCalID; // e vs detID
    inline TH2F* hxfCalVxnCal[NARRAY];

    inline TH2F* heCalVz;
    inline TH2F* heCalVzGC;
    inline TH2F* hecalVzRow[NROW];
    inline TH2F* hecalVzRowG[NROW];

    //====== Ex data
    inline TH1F* hEx;
    inline TH1F* hExi[NARRAY];
    inline TH2F* hExVxCal[NARRAY];
    inline TH1F* hExc[NARRAY/NROW];

    inline TH2F* hExThetaCM;

    inline TH1F* hExCut1;
    inline TH1F* hExCut2;

    //====== TAC
    inline TH1F* htac;   // by TAC
    inline TH1F* htac2;  // by timestamp
    inline TH2F* htacEx;
    inline TH2F* htac2Ex;

    inline TH1I* htacArray[NARRAY];

    inline TH2F* htacTdiff;
    inline TH2F* htacTdiffg;

    inline TH2F* htacRecoil[8];
    inline TH2F* htacRecoilsum[4];

    //======= APOLLO
    inline TH1F* hApollo[20];

    //======= Recoil
    inline TH2F* hrdtID;
    inline TH1F* hrdt[NRDT]; // single recoil
    inline TH1F* hrdtg[NRDT];

    inline TH2F* hrdt2D[4];
    inline TH2F* hrdt2Dsum[4];
    inline TH2F* hrdt2Dg[4];

    inline TH2F* hrdtMatrix; // coincident between rdt

    inline TH1F* hrdtRate1;
    inline TH1F* hrdtRate2;

    //======= Circular Recoil
    inline TH2F* hcrdtID;
    inline TH1F* hcrdt[16];
    inline TH2F* hcrdtPolar;

    //======= ELUM
    inline TH1F* helum[NELUM];
    inline TH2F* helumID;

    inline TH1F* helum4D; // elum rate for (d,d)
    inline TH1F* helum4C; // elum rate for (12C, 12C)

    inline TH1F* hBIC; // BIC, beam integrated current

    //======= EZero, or IonChamber when recoil also use
    inline TH1F* hic0; //ionChamber ch0
    inline TH1F* hic1;
    inline TH1F* hic2;

    inline TH2F* hic01; //ionChamber ch0-ch1
    inline TH2F* hic02;
    inline TH2F* hic12;

    //======= multi-Hit
    inline TH2I *hmult;
    inline TH1I *hmultEZ;
    inline TH2I *hArrayRDTMatrix;
    inline TH2I *hArrayRDTMatrixG;

    //======= ARRAY-RDT time diff
    inline TH1I *htdiff;
    inline TH1I *htdiffg;

    /***************************
     ***************************/
    //==== global variables
    inline Float_t x[NARRAY],z[NARRAY];
    inline Float_t xcal[NARRAY],xfcal[NARRAY],xncal[NARRAY],eCal[NARRAY];
    inline Int_t tacA[NARRAY];

    //==== correction parameters
    inline Float_t xnCorr[NARRAY];
    inline Float_t xfxneCorr[NARRAY][2];
    inline Float_t xScale[NARRAY];
    inline Float_t eCorr[NARRAY][2];
    inline Float_t rdtCorr[8][2];
    //==== parameters for Ex and thetaCM calcualtion

    inline double zRange[2] = {-1000, 0}; // zMin, zMax

    inline double Ex, thetaCM;
    inline double q, alpha, Et, betRel, gamm, G, massB, mass; //variables for Ex calculation
    inline bool isReaction;

    inline int padID = 0;

    inline TLatex text;



}

#endif

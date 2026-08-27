//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr  4 14:39:10 2017 by ROOT version 5.32/03
// from TTree psd_tree/PSD Tree
// found on file: psd.root
//////////////////////////////////////////////////////////
#ifndef Monitors_h
#define Monitors_h

class TBranch;
class TCanvas;
class TChain;
class TCutG;

#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TLatex.h"
#include "TList.h"
#include "TObject.h"
#include "TSelector.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TTree.h"

#include "GeneralSortMapping.h"
#include "../Armory/AnalysisLibrary.h"

// Header file for the classes stored in the TTree if any.
// Fixed size dimensions of array or collections stored in the TTree if any.

class Monitors : public TSelector {
    public :
        TTree          *fChain;   //!pointer to the analyzed TTree or TChain

        // Declaration of leaf types
        UShort_t        runID;
        Float_t         e[NARRAY];
        ULong64_t       e_t[NARRAY];
        Float_t         xf[NARRAY];
        ULong64_t       xf_t[NARRAY];
        Float_t         xn[NARRAY];
        ULong64_t       xn_t[NARRAY];
        Float_t         ring[NARRAY];
        ULong64_t       ring_t[NARRAY];
        Float_t         rdt[NRDT];
        ULong64_t       rdt_t[NRDT];
        Float_t         tac[NTAC];
        ULong64_t       tac_t[NTAC];
        Float_t         elum[NELUM];
        ULong64_t       elum_t[NELUM];
        Float_t         ezero[NEZERO];
        ULong64_t       ezero_t[NEZERO];
        Float_t         crdt  [NCRDT];
        ULong64_t       crdt_t[NCRDT];
        Float_t         apollo  [NAPOLLO];
        ULong64_t       apollo_t[NAPOLLO];

        // List of branches
        TBranch        *b_runID;   //!
        TBranch        *b_Energy;   //!
        TBranch        *b_EnergyTimestamp;   //!
        TBranch        *b_XF;   //!
        TBranch        *b_XFTimestamp;   //!
        TBranch        *b_XN;   //!
        TBranch        *b_XNTimestamp;   //!
        TBranch        *b_Ring;   //!
        TBranch        *b_RingTimestamp;   //!
        TBranch        *b_RDT;   //!
        TBranch        *b_RDTTimestamp;   //!
        TBranch        *b_TAC;   //!
        TBranch        *b_TACTimestamp;   //!
        TBranch        *b_ELUM;   //!
        TBranch        *b_ELUMTimestamp;   //!
        TBranch        *b_EZERO;   //!
        TBranch        *b_EZEROTimestamp;   //!
        TBranch        *b_CRDT;   //!
        TBranch        *b_CRDTTimestamp;   //!
        TBranch        *b_APOLLO;   //!
        TBranch        *b_APOLLOTimestamp;   //!

        // trace analysis data
        Float_t         te[NARRAY];
        Float_t         te_r[NARRAY];
        Float_t         te_t[NARRAY];
        Float_t         trdt[NRDT];
        Float_t         trdt_t[NRDT];
        Float_t         trdt_r[NRDT];
        Float_t         tcrdt[NCRDT];
        Float_t         tcrdt_t[NCRDT];
        Float_t         tcrdt_r[NCRDT];

        TBranch        *b_Trace_Energy;   //!
        TBranch        *b_Trace_Energy_RiseTime;   //!
        TBranch        *b_Trace_Energy_Time;   //!
        TBranch        *b_Trace_RDT;   //!
        TBranch        *b_Trace_RDT_Time;   //!
        TBranch        *b_Trace_RDT_RiseTime;   //!
        TBranch        *b_Trace_CRDT;   //!
        TBranch        *b_Trace_CRDT_Time;   //!
        TBranch        *b_Trace_CRDT_RiseTime;   //!

        int numCol; // number of detector in same side
        int numRow; // number of sides

        bool isArrayTraceExist;
        bool isRDTTraceExist;

        bool isRDTExist;
        bool isCRDTExist;
        bool isTACExist;
        bool isELUMExist;
        bool isEZEROExist;
        bool isAPOLLOExist;

        int printControlID;
        void printControl(int n){this->printControlID = n;}

        DetGeo detGeo;
        ReactionConfig reactionConfig;

        void LoadDetGeoAndReactionConfigFile();
        void LoadReactionPars();
        void LoadXFXNCorr();
        void LoadXFXN2ECorr();
        void LoadXScaleCorr();
        void LoadECorr();
        void LoadRDTCorr();

        double z[100];

        ULong64_t startTime ;
        ULong64_t endTime ;

        Monitors(TTree * /*tree*/ =0) : fChain(0) { }
        virtual ~Monitors() { }
        virtual Int_t   Version() const { return 2; }
        virtual void    Begin(TTree *tree);
        virtual void    SlaveBegin(TTree *tree);
        virtual void    Init(TTree *tree);
        virtual Bool_t  Notify();
        virtual Bool_t  Process(Long64_t entry);
        virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
        virtual void    SetOption(const char *option) { fOption = option; }
        virtual void    SetObject(TObject *obj) { fObject = obj; }
        virtual void    SetInputList(TList *input) { fInput = input; }
        virtual TList  *GetOutputList() const { return fOutput; }
        virtual void    SlaveTerminate();
        virtual void    Terminate();


        TString FindStartEndTime(TTree * tree,TString BranchName);
        void Draw2DHist(TH2F * hist, int colzThreshold = 3000);

        void PlotEZ(bool isRaw);
        void PlotTDiff(bool isGated, bool isLog);
        void PlotRDT(int id, bool isRaw);
        void PlotCRDTPolar();

        ClassDef(Monitors,0);
};

#define NROW 4 // number of side of array
#define NCOL  (NARRAY / NROW) // 12 detectors, 3 rows
//---histogram setting
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

#endif

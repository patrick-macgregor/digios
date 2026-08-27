//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr  4 14:39:10 2017 by ROOT version 5.32/03
// from TTree psd_tree/PSD Tree
// found on file: psd.root
//////////////////////////////////////////////////////////
#ifndef Monitors_h
#define Monitors_h

class TBranch;
class TChain;
class TCutG;
class TH2F;

#include "TList.h"
#include "TObject.h"
#include "TSelector.h"
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


#endif

#include "TH2F.h"
#include "TFile.h"
#include "TChain.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TCutG.h"
#include "TString.h"
#include "TObjArray.h"

#include <fstream>

namespace rdtcutcreator{
    Float_t rdtCorr2[8][2];

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void LoadRDTCorr2(){
        //========================================= e correction
        printf(" loading rdt correction.");
        std::ifstream file;
        file.open("correction_rdt.dat");

        if( file.is_open() ){
            double a, b;
            int i = 0;
            while( file >> a >> b){
                if( i >= 8) break;
                rdtCorr2[i][0] = a; /// a1
                rdtCorr2[i][1] = b; /// a0, e' = e * a1 + a0
                i++;
            }
            printf("............ done.\n");
            file.close();
        }
        else{
            printf("............ fail.\n");
            for( int i = 0; i < 8 ; i++){
                rdtCorr2[i][0] = 1.;
                rdtCorr2[i][1] = 0.;
            }
        }

        for( int i = 0; i < 8; i++){
            printf("%d | %f, %f \n", i, rdtCorr2[i][0], rdtCorr2[i][1]);
        }
        printf("=================================\n");

        return;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    TCanvas* InitialiseCanvas(bool isLogz){
        TCanvas * cCutCreator = new TCanvas("cCutCreator", "RDT Cut Creator", 100, 100, 800, 800);
        if( !cCutCreator->GetShowToolBar() ) cCutCreator->ToggleToolBar();
        cCutCreator->Update();
        if( isLogz ) cCutCreator->cd()->SetLogz();
        gStyle->SetOptStat("neiou");
        return cCutCreator;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    int CreateCuts(TChain* chain, TCanvas*& cCutCreator, TObjArray*& cutList, const TString& gate, bool useTrace, const int eRange, const int deRange){
        TString varX, varY, tag;
        TCutG * cut = nullptr;
        TString expression[10];
        TH2F * h[4];
        int count = 0;

        for (Int_t i = 0; i < 4; i++){
            printf("======== make a graphic cut on the plot (double click to stop), %d-th cut: ", i );

            if (useTrace){
                varX.Form("trdt[%d]",2*i); varY.Form("trdt[%d]",2*i+1);
            }
            else{
                ///varX.Form("rdt[%d]",i+4); varY.Form("rdt[%d]",i); // dE grouped
                varX.Form("rdt[%d]",2*i); varY.Form("rdt[%d]",2*i+1);
            }
            h[i] = new TH2F(Form("h%d", i), Form("%s - %s", varY.Data(), varX.Data()), 500, 0, eRange, 500, 0, deRange);

            expression[i].Form("%s*%f +%f:%s*%f+%f>>h%d",
                    varY.Data(), rdtcutcreator::rdtCorr2[2*i+1][0], rdtcutcreator::rdtCorr2[2*i+1][1],
                    varX.Data(), rdtcutcreator::rdtCorr2[2*i][0],   rdtcutcreator::rdtCorr2[2*i][1],
                    i);

            ///printf("%s \n", expression[i].Data());

            chain->Draw(expression[i], gate, "col");

            if (h[i]->Integral() < 1000){
                h[i]->SetMarkerStyle(20);
                h[i]->SetMarkerSize(0.4);
                h[i]->Draw("");
            }

            cCutCreator->Modified(); cCutCreator->Update();

            gPad->WaitPrimitive("CUTG", "");

            cut = (TCutG*) gROOT->FindObject("CUTG");

            if( !cut ){
                printf(" stopped by user. no file saved or changed. \n");
                break;
            }

            TString name; name.Form("cut%d", i);
            cut->SetName(name);
            cut->SetVarX(varX.Data());
            cut->SetVarY(varY.Data());
            cut->SetTitle(tag);
            cut->SetLineColor(i == 2 ? 8 : i+1 );
            cutList->Add(cut);

            printf(" cut-%d \n", i);

            count++;
        }
        return count;
    }

} // namespace rdtcutcreator

void RDTCutCreator(TString dataList,
                   TString saveFileName = "rdtCuts_.root",
                   int eRange=6000,
                   int deRange=2000,
                   bool isLogz = false,
                   TString gate = "",
                   TString treeName = "gen_tree",
                   bool useTrace = false
    ){

    printf("================ Graphic Cut Creator for RDT ============== \n");

    rdtcutcreator::LoadRDTCorr2();
    TChain * chain = new TChain(treeName);
    chain->Add(dataList);
    chain->GetListOfFiles()->Print();
    TCanvas* cCutCreator = rdtcutcreator::InitialiseCanvas(isLogz);

    TObjArray * cutList = new TObjArray();
    int count = rdtcutcreator::CreateCuts(chain, cCutCreator, cutList, gate, useTrace, eRange, deRange);

    TFile* cutFile = new TFile(saveFileName, "recreate");
    cutList->Write("cutList", TObject::kSingleKey);

    printf("====> saved %d cuts into %s\n", count, saveFileName.Data());
    gROOT->ProcessLine(".q");
    return;
}


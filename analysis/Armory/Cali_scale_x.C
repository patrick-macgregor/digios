#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TSpectrum.h"
#include "TGraph.h"

#include <array>
#include <fstream>

namespace calibrationscaleX{
    const int rowDet = 4;
    const int colDet = 6;
    const int nDet = rowDet*colDet;

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    TCanvas* CreateCanvas(){
        gStyle->SetOptStat(11111111);
        Int_t Div[2] = {colDet, rowDet}; //x,y
        Int_t size[2] = {230,230}; //x,y
        TCanvas* cCali_scaleX = new TCanvas("cCali_scaleX", "cCali_scaleX", 0, 0, size[0]*Div[0], size[1]*Div[1]);
        cCali_scaleX->Divide(Div[0],Div[1]);

        for ( int i = 1; i <= Div[0]*Div[1] ; i++){
            cCali_scaleX->cd(i)->SetGrid();
        }

        gStyle->SetOptStat(0);
        gStyle->SetStatY(1.0);
        gStyle->SetStatX(0.99);
        gStyle->SetStatW(0.2);
        gStyle->SetStatH(0.1);
        gStyle->SetLineScalePS(1);

        if (cCali_scaleX->GetShowEditor())cCali_scaleX->ToggleEditor();
        if (cCali_scaleX->GetShowToolBar())cCali_scaleX->ToggleToolBar();
        return cCali_scaleX;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    bool LoadXFXNCorr(std::array<double,nDet>& xnCorr){
        std::ifstream file;
        file.open("correction_xf_xn.dat");

        if (!file.is_open()){
            printf("... fail.\n");
            return false;
        }

        double a;
        int i = 0;
        while(file >> a){
            if ( i >= nDet) break;
            xnCorr[i] = a;
            i = i + 1;
        }

        printf("... done.\n");
        file.close();
        return true;

    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    bool LoadXFXNECorr(std::array<std::array<double,2>,nDet>& xfxneCorr){
        std::ifstream file;
        file.open("correction_xfxn_e.dat");

        if (!file.is_open()){
            printf("... fail.\n");
            return false;
            // for (int i = 0; i < nDet; i++){
            //     xfxneCorr[i][0] = 0;
            //     xfxneCorr[i][1] = 1;
            // }
            // return true;
        }

        double a, b;
        int i = 0;
        while( file >> a >> b){
            if ( i >= nDet) break;
            xfxneCorr[i][0] = a;
            xfxneCorr[i][1] = b;
            i = i + 1;
        }

        printf("... done.\n");
        file.close();
        return true;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    bool LoadECorr(std::array<std::array<double,2>,nDet>& eCorr){
        std::ifstream file;
        file.open("correction_e.dat");

        if (!file.is_open()){
            printf("... fail.\n");
            return false;
            // for ( int i = 0; i < 24 ; i++){
            //     eCorr[i][0] = 1.;
            //     eCorr[i][1] = 0.;
            // }
            // return true;
        }

        double a, b;
        int i = 0;
        while( file >> a >> b){
            if ( i >= nDet) break;
            eCorr[i][0] = a;  // 1/a1
            eCorr[i][1] = b;  //  a0 , e' = e * a1 + a0
                              //printf("\n%2d, e0: %9.4f, e1: %9.4f", i, eCorr[i][0], eCorr[i][1]);
            i = i + 1;
        }
        printf("... done.\n");

        file.close();
        return true;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    std::array<TH1F*,nDet> CreateHistogramArray(){
        std::array<TH1F*,nDet> d;
        for ( int i = 0; i < nDet; i++){
            TString name;
            name.Form("d%d", i);
            d[i] = new TH1F(name, name, 240, -1.2, 1.2 );
            d[i]->SetXTitle("x");
            d[i]->SetYTitle("count / 0.01 unit");
        }
        return std::move(d);
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void CalculateXForEachEntryInTheTree(TTree*& tree,
                                         std::array<double,nDet>&               xnCorr,
                                         std::array<std::array<double,2>,nDet>& xfxneCorr,
                                         std::array<std::array<double,2>,nDet>& eCorr,
                                         std::array<TH1F*,nDet>&                d
    ){
        Float_t e[100];
        Float_t xf[100];
        Float_t xn[100];

        TBranch *b_Energy;//!
        TBranch *b_XF;//!
        TBranch *b_XN;//!

        tree->SetBranchAddress("e", e, &b_Energy);
        tree->SetBranchAddress("xf", xf, &b_XF);
        tree->SetBranchAddress("xn", xn, &b_XN);

        for ( Long64_t entry = 0; entry < tree->GetEntries(); entry ++){
            Long64_t local_entry = tree->LoadTree(entry);
            b_Energy->GetEntry(local_entry,0);
            b_XF->GetEntry(local_entry,0);
            b_XN->GetEntry(local_entry,0);

            for (int idet = 0 ; idet < nDet; idet++){
                double  eC = TMath::QuietNaN() ;
                double xfC = TMath::QuietNaN() ;
                double xnC = TMath::QuietNaN() ;
                double   x = TMath::QuietNaN() ;

                if (!TMath::IsNaN(e[idet]) || e[idet] > 0){
                    eC = e[idet]/eCorr[idet][0] + eCorr[idet][1];
                }
                else{
                    continue; // when e is invalid, nothing need to do
                }

                if (eC < 1.5) continue;

                if ( !TMath::IsNaN(xf[idet]) || xf[idet] > 0) xfC = xf[idet]*xfxneCorr[idet][1] + xfxneCorr[idet][0] ;
                if ( !TMath::IsNaN(xn[idet]) || xn[idet] > 0) xnC = xn[idet]*xnCorr[idet]*xfxneCorr[idet][1] + xfxneCorr[idet][0];

                //========= calculate x
                if (xf[idet] > 0 && xn[idet] > 0 ){
                    x = (xfC-xnC)/(xfC+xnC);
                }
                else if ((xf[idet] == 0 || TMath::IsNaN(xf[idet])) && (xn[idet] > 0 || !TMath::IsNaN(xn[idet]))){
                    x = (1-2*xnC/e[idet]);
                }
                else if ((xf[idet] > 0 || !TMath::IsNaN(xf[idet])) && (xn[idet] == 0 || TMath::IsNaN(xn[idet]))){
                    x = (2*xfC/e[idet]-1);
                }
                else{
                    continue;
                }

                d[idet]->Fill(x);
            } //detector
        } //entries
        return;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    std::array<double,nDet> FindPeaks(TCanvas*& cCali_scaleX, std::array<TH1F*,nDet>& d){
        TSpectrum* spec = new TSpectrum(50);
        double* x = nullptr;
        TF1* fit = new TF1("fit", "gaus");
        fit->SetNpx(1000);
        fit->SetLineColor(2);

        std::array<double,nDet> scaleX;

        for (int i = 0; i < nDet ; i ++){
            cCali_scaleX->cd(i+1);

            spec->Search(d[i], 2, "", 0.5);
            x = spec->GetPositionX();

            d[i]->Fit("fit", "Rq", "", x[0]-0.05, x[0]+0.05);

            scaleX[i] = TMath::Abs(fit->GetParameter(1));

            printf("%2d | x : %7.4f (%6.4f) \n", i, scaleX[i], fit->GetParError(1));

            cCali_scaleX->Update();
            gSystem->ProcessEvents();
        }
        return std::move(scaleX);
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void SaveScaleParameters(std::array<double,nDet>& scaleX){
        int dummy = 0;
        printf("0 for end, 1 for save x correction: ");
        int temp = scanf("%d", &dummy);
        if (dummy == 0) return;

        FILE* paraOut;
        paraOut = fopen ("correction_scaleX_.dat", "w+");

        for ( int i = 0; i < nDet; i++){
            if (i == 11) scaleX[i] = 1.0;
            fprintf(paraOut, "%9.6f\n", scaleX[i]);
        }

        fflush(paraOut);
        fclose(paraOut);

        printf("=========== save x-correction parameters to %s \n", "correction_scaleX.dat");
        return;

    }
}// anon. namespace

void Cali_scale_x(TTree *tree){
    printf("============================================================= \n");
    printf("==================== Cali_scale_x.C ========================= \n");
    printf("============================================================= \n");
    printf(" scale x to full range from (-1,1)\n");
    printf("------------------------------------------------- \n");
    printf("=========== Total #Entry: %10lld \n", tree->GetEntries());

    //======================================================== Browser or Canvas
    TCanvas* cCali_scaleX = calibrationscaleX::CreateCanvas();

    //========================================= xf = xn correction
    printf("----- loading xf-xn correction.");
    std::array<double,calibrationscaleX::nDet> xnCorr;
    if(!calibrationscaleX::LoadXFXNCorr(xnCorr)){return;}

        //========================================= e = xf + xn correction
    std::array<std::array<double,2>,calibrationscaleX::nDet> xfxneCorr;
    printf("----- loading xf/xn-e correction.");
    if(!calibrationscaleX::LoadXFXNECorr(xfxneCorr)){return;}

    //========================================= e correction
    std::array<std::array<double,2>,calibrationscaleX::nDet> eCorr;
    printf("----- loading e correction.");
    if(!calibrationscaleX::LoadECorr(eCorr)){return;}

    /**///========================================================= Analysis
    printf("creating x histogram for each detector.... please wait.\n");
    std::array<TH1F*,calibrationscaleX::nDet> d = calibrationscaleX::CreateHistogramArray();
    calibrationscaleX::CalculateXForEachEntryInTheTree(tree,xnCorr,xfxneCorr,eCorr,d);

    for ( int i = 0; i < calibrationscaleX::nDet; i ++){
        cCali_scaleX->cd(i+1);
        d[i]->Draw();
        cCali_scaleX->Update();
        gSystem->ProcessEvents();
    }

    //=================================================== Find the peak
    printf("======= Finding Peak and fit with a single Gauss \n");
    std::array<double,calibrationscaleX::nDet> scaleX = calibrationscaleX::FindPeaks(cCali_scaleX,d);

    //===== save correction parameter
    calibrationscaleX::SaveScaleParameters(scaleX);

    // Print canvas
    TString name = "plots/scaleX_calibration.pdf";
    cCali_scaleX->SaveAs(name);
    printf("Saved %s\n", name.Data());

}


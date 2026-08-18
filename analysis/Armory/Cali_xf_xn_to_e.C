#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TSpectrum.h>
#include <TGraph.h>
#include <fstream>
#include "../Armory/AnalysisLibrary.h"

namespace{
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void LoadDetectorGeometry(int& colDet, int& rowDet, int& nDet){
        std::string detGeoFileName = "detectorGeo.txt";
        printf("loading detector geometery : %s.", detGeoFileName.c_str());

        DetGeo detGeo;

        TMacro* haha = new TMacro();
        if( haha->ReadFile(detGeoFileName.c_str()) > 0 ) {
            detGeo = LoadDetectorGeo(haha);
            PrintDetGeo(detGeo);
            printf("... done.\n");
        }
        else{
            printf("... fail\n");
            return;
        }

        double length = detGeo.detLength;
        std::vector<double> pos = detGeo.detPos;

        colDet = detGeo.nDet;
        rowDet = detGeo.mDet;
        nDet = colDet * rowDet;
        delete haha;

    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    double* LoadXNXFCorrection(const int nDet){
        double* xnCorr = new double[nDet];
        std::ifstream file;
        file.open("correction_xf_xn.dat");
        if( file.is_open() ){
            double a;
            int i = 0;
            while( file >> a ){
                if( i >= nDet) break;
                xnCorr[i] = a;
                i = i + 1;
            }

            printf("... done.\n");
        }else{
            printf("... fail.\n");
            for( int i = 0 ; i < nDet ; i++){
                xnCorr[i] = 1;
            }
        }
        file.close();
        return xnCorr;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    TCanvas* CreateCorrectionCanvas(const int nDet, const int rowDet, const int colDet){
        Int_t Div[2] = {colDet,rowDet};  //x,y
        Int_t size[2] = {230,230}; //x,y
        TCanvas * cCali_xf_xn_e = new TCanvas("cCali_xf_xn_e", "cCali_xf_xn_e", 0, 0, size[0]*Div[0], size[1]*Div[1]);
        cCali_xf_xn_e->Divide(Div[0],Div[1]);

        for( int i = 1; i <= Div[0]*Div[1] ; i++){
            cCali_xf_xn_e->cd(i)->SetGrid();
        }

        gStyle->SetOptStat(0);
        gStyle->SetStatY(1.0);
        gStyle->SetStatX(0.99);
        gStyle->SetStatW(0.2);
        gStyle->SetStatH(0.1);

        if(cCali_xf_xn_e->GetShowEditor())cCali_xf_xn_e->ToggleEditor();
        if(cCali_xf_xn_e->GetShowToolBar())cCali_xf_xn_e->ToggleToolBar();

        return cCali_xf_xn_e;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    TH2F** Create2DCorrectionHistograms(TCanvas* cCali_xf_xn_e, const int nDet, double* xnCorr, TTree* tree){
        TH2F ** d = new TH2F*[nDet];
        for( int i = 0; i < nDet; i ++){
            TString name;
            name.Form("d%d", i);
            d[i] = new TH2F(name, name, 200, 0, 5000, 200, 0, 5000);
            d[i]->SetXTitle("xf+xn");
            d[i]->SetYTitle("e");

            TString expression;
            expression.Form("e[%d]:(xf[%d]+xn[%d] * %f)>> %s", i,i,i, xnCorr[i], name.Data());
            TString gate;
            // gate.Form("e[%d]>0 && xf[%d]>0 && xn[%d]>0 && ring[%d] < 50",i, i, i, i);
            gate.Form("e[%d]>0 && xf[%d]>0 && xn[%d]>0",i, i, i);
            //gate.Form("e[%d]>0 && xf[%d]>0 && xn[%d]>0 && TMath::Abs(xf[%d] - xn[%d]*%f)< 100",i, i, i, i, i, xnCorr[i]);

            cCali_xf_xn_e->cd(i+1);
            tree->Draw(expression, gate, "");
            cCali_xf_xn_e->Update();
            gSystem->ProcessEvents();
        }

        return d;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void FitProfileToCorrectionHistograms(Double_t* slope, Double_t* intep, TCanvas* cCali_xf_xn_e, const int nDet, TH2F** d){
        TF1 *fit = new TF1("fit", "pol1", 0, 3500);

        for( int i = 0; i < nDet; i ++){

            cCali_xf_xn_e->cd(i+1);

            fit->SetParameter(0, 0);
            fit->SetParameter(1, 1);

            d[i]->Fit("fit", "qR");

            slope[i] = fit->GetParameter(1);
            intep[i] = fit->GetParameter(0);

            printf("%2d, %9.6f, %9.6f \n", i, intep[i], slope[i]);

            cCali_xf_xn_e->Update();
            gSystem->ProcessEvents();

        }

        cCali_xf_xn_e->Update();

    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void SaveCorrectionParameters(const int nDet, Double_t* intep, Double_t* slope){
        int dummy = 0;
        printf("0 for end, 1 for save e-xf+xn correction: ");
        int temp = scanf("%d", &dummy);
        if( dummy == 0 ) return;
        if( dummy == 1 ){

            FILE * paraOut;
            paraOut = fopen ("correction_xfxn_e.dat", "w+");

            for( int i = 0; i < nDet; i++){

                fprintf(paraOut, "%9.6f  %9.6f\n", intep[i], slope[i]);
            }

            fflush(paraOut);
            fclose(paraOut);
        }
    }
} // anon. namespace








void Cali_xf_xn_to_e(TTree *tree){
    /**///======================================================== initial input

    int energyRange[2] = {500, 5000};

    gStyle->SetOptStat(11111111);

    printf("============================================================= \n");
    printf("=====================  Cali_xf_xn_e.C  ====================== \n");
    printf("============================================================= \n");
    printf("==== making xf + xn = e \n");
    printf("------------------------------------------------------------- \n");
    printf("=========== Total #Entry: %10lld \n", tree->GetEntries());

    /**///========================================================= load parameters
    int colDet, rowDet, nDet;
    printf("======================= loading parameters files .... \n");
    LoadDetectorGeometry(colDet, rowDet, nDet);
    printf("----- loading xf-xn correction.");
    double* xnCorr = LoadXNXFCorrection(nDet);

    /**///======================================================== Browser or Canvas
    TCanvas* cCali_xf_xn_e = CreateCorrectionCanvas(nDet, rowDet, colDet);

    /**///========================================================= Analysis

    //======== create histogram for each detector
    printf("creating xf-xn histogram for each detector.... please wait.\n");
    TH2F ** d = Create2DCorrectionHistograms(cCali_xf_xn_e, nDet, xnCorr, tree);

    //======== profileX
    printf("fitting slope.\n");
    Double_t* slope = new Double_t[nDet];
    Double_t* intep = new Double_t[nDet];
    FitProfileToCorrectionHistograms(slope, intep, cCali_xf_xn_e, nDet, d);

    //===== plot particular detector to examine further <PTM>
    TCanvas* cPTM = new TCanvas("cPTM", "cPTM", 1400, 900);
    cPTM->cd();
    const int CHOSEN_DETECTOR = 13;
    TH2F* h2 = new TH2F("CHOSEN", "CHOSEN", 200, 0, 5000, 200, 0, 5000);
    h2->SetXTitle("xf+xn");
    h2->SetYTitle("e");
    TString expression = Form("e[%d]:(xf[%d]+xn[%d] * %f)>>CHOSEN", CHOSEN_DETECTOR, CHOSEN_DETECTOR, CHOSEN_DETECTOR, xnCorr[CHOSEN_DETECTOR]);
    TString gate = Form("e[%d]>0 && xf[%d]>0 && xn[%d]>0",CHOSEN_DETECTOR, CHOSEN_DETECTOR, CHOSEN_DETECTOR);
    tree->Draw(expression, gate, "");
    cPTM->Modified(); cPTM->Update();
    gSystem->ProcessEvents();
    // </PTM>

    //===== save correction parameter
    SaveCorrectionParameters(nDet, intep, slope);
    printf("=========== save xfxn-e-correction parameters to %s \n", "correction_xfxn_e.dat");

        //cCali_xf_xn_e->SaveAs("xfxn_e_correction.pdf");
        //printf("=========== save canvas to %s \n", "xfxn_e.pdf");

}
    /*
    // make correction
    TH2F ** dc = new TH2F[nDet];
    for( int i = 0; i < nDet; i ++){
    TString name;
    name.Form("dc%d", i);
    dc[i] = new TH2F(name, name , 300, -100 , 3500 , 300, -100 , 3500);
    dc[i]->SetXTitle("xf+xn");
    dc[i]->SetYTitle("e");

    TString expression;
    expression.Form("e[%d]:%f*(xf[%d]+xn[%d])+%f>> dc%d" , i, slope[i] , i,i, intep[i], i);
    TString gate;
    gate.Form("e[%d]!=0 && xf[%d]!=0 && xn[%d]!=0",i,i,i);
    //gate.Form("xf[%d]!=0 && xn[%d]!=0",i,i);
    tree->Draw(expression, gate , "");

    //printf("========  i = %d \n", i);
    }
    */



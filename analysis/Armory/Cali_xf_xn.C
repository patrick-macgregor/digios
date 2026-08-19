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
#include "TGraph.h"
#include "TLine.h"
#include "TSpectrum.h"

#include "../Armory/AnalysisLibrary.h"
#include "../Armory/AutoFit.C"
#include "../working/GeneralSortMapping.h"
#include "../working/scripts/fitparameters.h"
#include "RtypesCore.h"

#include <array>
#include <vector>

namespace{
    // Constants
    const int rowDet = 4;
    const int colDet = NARRAY/rowDet;
    const int numDet = rowDet*colDet;

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void Kill(){
        gROOT->ProcessLine(".q");
        return;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void PrintFunctionIntro(const Long64_t entries){
        printf("============================================================= \n");
        printf("====================== Cali_xf_xn.C ========================= \n");
        printf("============================================================= \n");
        printf("==   calibration for PSD detectors using alpha souce. \n");
        printf("------------------------------------------------------------- \n");
        printf("1, calibration energy using charateristic energy of alpha souce. \n");
        printf("2, calibration xf-xn with energy-gate. \n");
        printf("------------------------------------------------------------- \n");
        printf("=========== Total #Entry: %10lld \n", entries);
        return;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void GetEnergyRangeFromUser(int* energyRange){
        int temp = 0;
        printf(" Min Raw Energy [ch] (default = %d ch) : ", energyRange[1]);
        scanf("%d", &temp);
        energyRange[1] = temp;

        printf(" Max Raw Energy [ch] (default = %d ch) : ", energyRange[2]);
        scanf("%d", &temp);
        energyRange[2] = temp;

        printf("     Raw Energy is now %d ch to %d ch\n", energyRange[1], energyRange[2]);
        return;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    TCanvas* InitialiseCanvas(Int_t* size, Int_t* Div){
        TCanvas* cAlpha = new TCanvas("cAlpha", "cAlpha", 0, 0, size[0]*Div[0], size[1]*Div[1]);
        cAlpha->Divide(Div[0],Div[1]);

        for( int i = 1; i <= Div[0]*Div[1] ; i++){
            cAlpha->cd(i)->SetGrid();
        }

        gStyle->SetOptStat(0);
        gStyle->SetStatY(1.0);
        gStyle->SetStatX(0.99);
        gStyle->SetStatW(0.2);
        gStyle->SetStatH(0.1);
        gStyle->SetLineScalePS(1);

        if (cAlpha->GetShowEditor())cAlpha->ToggleEditor();
        if (cAlpha->GetShowToolBar())cAlpha->ToggleToolBar();
        return cAlpha;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void InitialiseAlphaSpectraHistograms(int* energyRange, TTree* tree, TCanvas* cAlpha, TH1F**& q, TString*& gate){
        for( int i = 0; i < numDet; i++){
            TString name;
            name.Form("q%d", i);
            q[i] = new TH1F(name, name, energyRange[0], energyRange[1], energyRange[2]);
            q[i]->SetLineWidth(1);
            q[i]->SetXTitle(name);

            TString expression;
            expression.Form("e[%d] >> q%d",i, i);
            //gate[i].Form("ring[%d]==0 && !TMath::IsNaN(xf[%d]) && !TMath::IsNaN(xn[%d])", i, i, i);
            //gate[i].Form("!TMath::IsNaN(xf[%d]) && !TMath::IsNaN(xn[%d])", i, i);
            //gate[i].Form("e[%d] > 0", i);

            cAlpha->cd(i+1);
            tree->Draw(expression, gate[i], "");
            cAlpha->Update();
            gSystem->ProcessEvents();
        }
        // For some reason these lines are needed
        cAlpha->Update();
        gSystem->ProcessEvents();
        return;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    int GetMethodToFindEdgeFromUser(){
        //   printf("1) maximum peak edge.\n");
        printf("2) TSpectrum search\n");
        printf("3) TSpectrum fit.\n");
        printf("9) Exit.\n");
        printf("====== choice : ");
        int method = 0;
        scanf("%d", &method);
        return method;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void FindPeaksUsingTSpectrumClass(double& threshold, TH1F**& q, std::array<std::vector<double>,numDet>& energy, TCanvas*& cAlpha){
        printf("------ finding peak using TSpectrum Class...\n");

        printf(" peak threshold (default = %.3f) : ", threshold);
        scanf("%lf", &threshold);
        printf("     threshold is now %.3f\n", threshold);

        for( int i = 0; i < numDet; i++){

            TSpectrum* spec = new TSpectrum();
            int nPeaks = spec->Search(q[i], 3, "", threshold);
            printf("%2d | found %d peaks | ", i,  nPeaks);

            double* xpos = spec->GetPositionX();
            double* ypos = spec->GetPositionY();

            std::vector<double> height;

            int* inX = new int[nPeaks];
            TMath::Sort(nPeaks, xpos, inX, 0 );
            for( int j = 0; j < nPeaks; j++){
                energy[i].push_back(xpos[inX[j]]);
                height.push_back(ypos[inX[j]]);
            }

            for( int j = 0; j < nPeaks; j++){
                printf("%7.2f, ", energy[i][j]);
            }
            printf("\n");
        }

        for( int i = 0; i < numDet; i++){
            cAlpha->cd(i+1);
            q[i]->Draw();
            cAlpha->Update();
            gSystem->ProcessEvents();
        }
        return;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void FindPeaksUsingAutoFit(TCanvas* cAlpha, std::array<std::vector<double>,numDet>& energy, TH1F**& q){
        printf("------ find and fit peak with Gaussian using AutoFit.C \n");
        for( int i = 0; i < numDet; i++){
            cAlpha->cd(i+1);
            energy[i] = fitAuto(q[i], -1, 0.3, 20, 4, "");
            int nPeaks = energy[i].size();
            printf("%2d | found %d peaks | ", i,  nPeaks);
            for( int j = 0; j < nPeaks; j++){
                printf("%7.2f, ", energy[i][j]);
            }
            printf("\n");
        }
        return;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    int GetDetectorToBeTheReferenceFromUser(){
        int refID = 0;
        printf("========== which detector to be the reference?\n");
        printf(" X =  det-X reference\n");
        printf("-1 =  manual reference\n");
        printf("-2 =  use 228Th, first 7 strongest peaks \n");
        printf("-3 =  use 148Gd + 244Cm, (3.1828, 5.8048 MeV)\n");
        //printf("-3 =  use 241Am, 5.481 MeV \n");
        printf("-9 =  stop \n");
        printf("your choice = ");
        scanf("%d", &refID);
        return refID;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void UseDetectorAsReference(std::array<std::vector<double>,numDet>& energy, std::vector<double>& refEnergy, const int refID){
        int n = energy[refID].size();
        for (int k = 0; k < n; k++) refEnergy.push_back(energy[refID][k]);
        return;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void UseManualAsReference(std::vector<double>& refEnergy){
        int n = 0;
        float eng = -1;
        do{
            printf("%2d-th peak energy (< 0 to stop):", n);
            int temp = scanf("%f", &eng);
            printf("             input: %f \n", eng);
            if ( eng >= 0 ) refEnergy.push_back(eng);
            n++;
        } while(eng >= 0);
        return;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void Use228ThAsReference(std::vector<double>& refEnergy){
        refEnergy.clear();
        refEnergy.push_back(5.34);
        refEnergy.push_back(5.423);
        refEnergy.push_back(5.685);
        refEnergy.push_back(6.050);
        refEnergy.push_back(6.288);
        refEnergy.push_back(6.778);
        refEnergy.push_back(8.785);
        return;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void Use148Gd244CmAsReference(std::vector<double>& refEnergy){
        refEnergy.clear();
        refEnergy.push_back(3.1828);
        refEnergy.push_back(5.8048);
        return;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void DoCalibration(std::array<std::vector<double>,numDet>& energy, std::vector<double>& refEnergy, const int refID, TCanvas*& cAlpha, std::array<FitParameters,numDet>& fp){
        for( int i = 0; i < numDet; i ++){
            printf("------ refID - %d, nPeaks: %lu \n", i, energy[i].size());

            if ( refID >= 0 && refID == i ){
                fp[i] = FitParameters();
                fp[i].intercept = 0;
                fp[i].slope = 1;
                printf("skipped - itself\n");
                continue;
            }

            if ( energy[i].size() == 0) {
                fp[i] = FitParameters();
                fp[i].intercept = 0;
                fp[i].slope = 1;
                printf("skipped\n");
                continue;
            }

            printf("   Energy : ");
            for( int k = 0; k < energy[i].size(); k++){ printf("%.1f, ", energy[i][k]);};printf("\n");
            std::vector<std::vector<double>> output =  FindMatchingPair(energy[i], refEnergy);

            std::vector<double> haha1 = output[0];
            std::vector<double> haha2 = output[1];

            TGraph* graph = new TGraph(haha1.size(), &haha1[0], &haha2[0] );
            cAlpha->cd(i+1);
            graph->Draw("A*");

            TF1* fit = new TF1("fit", "pol1" );
            TFitResultPtr r = graph->Fit("fit", "qS");
            if (r.Get()){
                fp[i] = FitParameters(fit, r);
            }
            else{
                fp[i] = FitParameters();
            }

            printf("%2d | a0: %6.3f, a1: %6.3f (%14.8f) \n", i, fp[i].intercept, fp[i].slope, 1./fp[i].slope);

        }
        return;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void PlotAdjustedSpectrum(const int* const energyRange, const std::vector<double>& refEnergy, const std::array<FitParameters, numDet>& fp, TString* gate, TTree* tree, const int method){
        TCanvas* cAux = new TCanvas ("cAux", "cAux", 600, 400);
        TH2F** p = new TH2F*[numDet];
        double yMax = 0;

        for ( int i = 0; i < numDet; i ++){
            TString name;
            name.Form("p%d", i);
            p[i] = new TH2F(name, name, 24, 0, 24, energyRange[0], refEnergy[0]*0.9, refEnergy.back()*1.1);

            TString expression;
            expression.Form("e*%.8f + %.8f : Iteration$ >> p%d", fp[i].slope, fp[i].intercept, i);
            gate[i].Form("Iteration$ == %d", i);
            tree->Draw(expression, gate[i] , "colz");

            cAux->Update();
            gSystem->ProcessEvents();
        }

        cAux->cd(1);
        TString hhhName;
        if ( method == 2 ) hhhName="Calibrated energy using TSpectrum peak";
        if ( method == 3 ) hhhName="Calibrated energy using Gaussian fitting";
        TH2F* hhh = new TH2F("hhh", hhhName, 24, 0, 24, energyRange[0], refEnergy[0]*0.9, refEnergy.back()*1.1);
        for( int  i = 0; i < numDet; i++){
            hhh->Add(p[i]);
        }
        hhh->Draw("colz");
        cAux->Update();
        gSystem->ProcessEvents();

        return;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void SaveEAlphaFitParameters(const std::array<FitParameters, numDet>& fp){
        printf("0 for no-save, any other integer for save e-alpha calibration: ");
        int dummy = 0;
        scanf("%d", &dummy);
        if ( dummy == 0 ) return;

        FILE* paraOut;
        TString filename;
        filename.Form("correction_e_alpha_.dat");
        paraOut = fopen (filename.Data(), "w+");
        printf("=========== save e-correction parameters to %s \n", filename.Data());
        for( int i = 0; i < numDet; i++){
            fprintf(paraOut, "%14.8f\t%9.6f\n", 1./fp[i].slope, fp[i].intercept);
        }
        fflush(paraOut);
        fclose(paraOut);
        return;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    double GetEGate(const std::vector<double>& refEnergy){
        double eGate = 0;
        gSystem->ProcessEvents();
        int peakID = 0;
        printf("------ pick the i-th peak (0, 1, ... , %d, -1 to stop): ", (int) refEnergy.size() - 1);
        int temp = scanf("%d", &peakID);
        if ( peakID < 0 ) {
            Kill();
        }
        eGate = refEnergy[peakID];
        printf("------ using the peak at : %f \n", eGate);
        return eGate;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void PlotXNXFHistograms(TH2F**& h,TCanvas*& cAlpha, const int* const energyRange, TString*& gate, TTree* tree, const double eGate, const std::array<FitParameters,numDet>& fp){
        for( int i = 0; i < numDet; i ++){
            TString name = Form("h%d", i);
            h[i] = new TH2F(name, name,  energyRange[0], 0, energyRange[2], energyRange[0], 0, energyRange[2]);
            name.Form("xf[%d]", i); h[i]->SetYTitle(name);
            name.Form("xn[%d]", i); h[i]->SetXTitle(name);

            TString expression = Form("xf[%d]:xn[%d]>> h%d", i, i, i);
            gate[i].Form("xf[%d]>0 && xn[%d]>0 && TMath::Abs(e[%d]*%f + %f - %f)< %f", i, i, i, fp[i].slope, fp[i].intercept, eGate, eGate*0.05);

            cAlpha->cd(i+1);
            tree->Draw(expression, gate[i], "colz");

            cAlpha->Update();
            gSystem->ProcessEvents();

        }
        // For some reason these lines are needed
        cAlpha->Update();
        gSystem->ProcessEvents();
        return;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void FitXNXFProfileFunction(TCanvas*& cAlpha, std::array<std::vector<double>, numDet>& energy, TH2F**& h, std::array<FitParameters,numDet>& fpxnxf){
        TF1* fit = new TF1("fit", "pol1");
        for( int i = 0; i < numDet; i++){
            cAlpha->cd(i+1);

            if ( energy[i].size() == 0) {
                fpxnxf[i] = FitParameters();
                printf("skipped - detID-%d\n", i);
                continue;
            }

            TFitResultPtr r = h[i]->ProfileX()->Fit("fit", "QS");
            if (!r.Get()){
                fpxnxf[i] = FitParameters();
            }
            else{
                fpxnxf[i] = FitParameters(fit,r);
            }
        }
        return;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void CorrectXNXFFunctions(const int* const energyRange, const std::array<FitParameters,numDet> fpxnxf, TCanvas*& cAlpha, TTree* tree, TString* gate){
        TH2F** k = new TH2F*[numDet];
        TLine line(0,0,0,0);
        line.SetLineColor(4);
        for( int i = 0; i < numDet; i ++){
            TString name;
            name.Form("k%d", i);
            k[i] = new TH2F(name, name,  energyRange[0], 0, energyRange[2],  energyRange[0], 0, energyRange[2]);
            name.Form("xf[%d]", i); k[i]->SetYTitle(name);
            name.Form("xn[%d]", i); k[i]->SetXTitle(name);

            TString expression;
            expression.Form("xf[%d]:xn[%d]*%f>> k%d" ,i ,i, -fpxnxf[i].slope, i);
            //gate[i].Form("xf[%d]!=0 && xn[%d]!=0", i, i);

            cAlpha->cd(i+1);

            tree->Draw(expression, gate[i] , "colz");
            line.SetX2(fpxnxf[i].intercept);
            line.SetY1(fpxnxf[i].intercept);
            line.Draw("same");

            cAlpha->Update();
            gSystem->ProcessEvents();
        }
        // For some reason these lines are needed
        cAlpha->Update();
        gSystem->ProcessEvents();

        return;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void SaveXNXFCorrectionParameters(const std::array<FitParameters,numDet>& fpxnxf){
        printf("0 for no-save, any other integer for save xf-xn-correction : ");
        int dummy = 0;
        scanf("%d", &dummy);
        if ( dummy == 0 ) return;

        FILE* paraOut;
        TString filename;
        filename.Form("correction_xf_xn_.dat");
        paraOut = fopen (filename.Data(), "w+");
        printf("=========== save xf_xn-correction parameters to %s \n", filename.Data());
        for( int i = 0; i < numDet; i++){
            fprintf(paraOut, "%9.6f\n", -fpxnxf[i].slope);
        }
        fflush(paraOut);
        fclose(paraOut);
        //cAlpha->SaveAs("alpha_xf_xn_corrected.pdf");
        return;
    }
} // anon. namespace



void Cali_xf_xn(TTree* tree){
    // Initialise variables
    int energyRange[3] = {400, 1000, 2600}; // bin, min, max - these are DEFAULTS
    double threshold = 0.2;
    Int_t Div[2] = {colDet,rowDet};  //x,y
    Int_t size[2] = {230,230}; //x,y

    // Print header
    PrintFunctionIntro(tree->GetEntries());

    // Initialise canvas and energy range
    GetEnergyRangeFromUser(energyRange);
    TCanvas* cAlpha = InitialiseCanvas(size, Div);

    printf("############## e correction \n");
    TH1F** q = new TH1F*[numDet];
    TString* gate = new TString[numDet];
    InitialiseAlphaSpectraHistograms(energyRange, tree, cAlpha, q, gate);

    //----------- 2, find the edge of the energy
    printf("============== method to find edge:\n");
    int method = GetMethodToFindEdgeFromUser();
    if (method == 9)Kill();

    std::array<std::vector<double>,numDet> energy;
    std::vector<double> refEnergy;

    if ( method == 2 )     FindPeaksUsingTSpectrumClass(threshold, q, energy, cAlpha);
    else if ( method == 3) FindPeaksUsingAutoFit(cAlpha, energy, q);
    TString name = "plots/alpha_calibration.pdf";
    cAlpha->SaveAs(name);
    printf("Saved %s\n", name.Data());

    //------------ 3, correction
    int refID = GetDetectorToBeTheReferenceFromUser();
    if (refID == -9) Kill();
    else if (refID >=  0) UseDetectorAsReference(energy, refEnergy, refID);
    else if (refID == -1) UseManualAsReference(refEnergy);
    else if (refID == -2) Use228ThAsReference(refEnergy);
    else if (refID == -3) Use148Gd244CmAsReference(refEnergy);

    printf("------ adjusting the energy to det-%d......\n", refID);
    for( int k = 0; k < refEnergy.size(); k++) printf("%2d-th peak : %f \n", k,  refEnergy[k]);

    const std::vector<double> refEnergy0 = refEnergy;
    std::array<FitParameters, numDet> fp;

    DoCalibration(energy, refEnergy, refID, cAlpha, fp);
    PlotAdjustedSpectrum(energyRange, refEnergy, fp, gate, tree, method);
    SaveEAlphaFitParameters(fp);
    gSystem->ProcessEvents();

    //############################################################  for xf-xn correction
    printf("############## xf - xn correction \n");
    double eGate = GetEGate(refEnergy);

    printf("------ plotting xf vs xn with energy gate near the peak...\n");
    TH2F** h = new TH2F*[numDet];
    PlotXNXFHistograms(h,cAlpha,energyRange,gate,tree,eGate,fp);

    printf("------ profile and obtain the fit function...\n");
    std::array<FitParameters,numDet> fpxnxf;
    FitXNXFProfileFunction(cAlpha, energy, h, fpxnxf);

    printf("------ correcting...\n");
    CorrectXNXFFunctions(energyRange, fpxnxf, cAlpha, tree, gate);
    //
    //--------- 4, pause for saving correction parameter
    SaveXNXFCorrectionParameters(fpxnxf);
    name = "plots/xnxf_calibration.pdf";
    cAlpha->SaveAs(name);
    printf("Saved %s\n", name.Data());
    return;
}



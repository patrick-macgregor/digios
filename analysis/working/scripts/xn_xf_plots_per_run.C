#include <fstream>
#include <iostream>
#include <vector>

#include "RtypesCore.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "../../Armory/AnalysisLibrary.h"

const std::vector<int> kRunNumbers = {
     9, 10, 11, 12, 13,
    14, 15, 16, 17,
    20, 21, 22, 23,
    24, 25
};

struct FitParameters{
    double intercept;
    double slope;
    double error_intercept;
    double error_slope;
    double covariance_intercept_slope;
    double chi2;
    double ndf;
    double prob;

    FitParameters(TF1* fitfunc, TFitResultPtr r){
        intercept = fitfunc->GetParameter(0);
        slope = fitfunc->GetParameter(1);
        error_intercept = fitfunc->GetParError(0);
        error_slope = fitfunc->GetParError(1);
        covariance_intercept_slope = r->CovMatrix(0,1);
        chi2 = r->Chi2();
        ndf = r->Ndf();
        prob = r->Prob();
    }

    FitParameters(){
        intercept = 0.0;
        error_intercept = 0.0;
        slope = 0.0;
        error_slope = 0.0;
        covariance_intercept_slope = 0.0;
        chi2 = 0.0;
        ndf = 0.0;
        prob = 0.0;
    }
};


namespace{
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void LoadDetectorGeometry(int& colDet, int& rowDet, int& nDet){
        std::string detGeoFileName = "../detectorGeo.txt";
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
        file.open("../correction_xf_xn.dat");
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
    void FitProfileToCorrectionHistograms(TCanvas* c, const int i, TH2F* d, FitParameters*& fp){
        TF1 *fit = new TF1("fit", "pol1", 0, 3500);

        c->cd(i+1);

        fit->SetParameter(0, 0);
        fit->SetParameter(1, 1);

        TFitResultPtr p = d->Fit("fit", "qRS");
        TString warn = "";

        if(!p.Get()){
            fp = new FitParameters();
            warn = "[DUMMY]";
        }
        else{
            fp = new FitParameters(fit,p);
        }

        printf("%2d, %9.6f, %9.6f %s\n", kRunNumbers[i], fp->intercept, fp->slope, warn.Data());

        c->Update();
        gSystem->ProcessEvents();

        c->Update();
        return;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: //
    void SaveFitToFile(const int N, const int CHOSEN_DETECTOR, FitParameters** fp){
        TString outfilename = Form("xnxf_detector_%d_fit_pars.txt", CHOSEN_DETECTOR);
        std::ofstream outfile;
        outfile.open(outfilename.Data());

        if (!outfile.is_open()){
            printf("CANNOT SAVE FIT TO FILE\n");
            return;
        }

        outfile << "y = Mx + C\n";
        outfile << "Run\t C\tE[C]\t M\tE[M]\t COV[M,C]\t Chi2\t Ndf\t Prob\n";

        for (int i = 0; i < N; ++i){
            outfile << kRunNumbers[i] << "\t" <<
                fp[i]->intercept << "\t" <<
                fp[i]->error_intercept << "\t" <<
                fp[i]->slope << "\t" <<
                fp[i]->error_slope << "\t" <<
                fp[i]->covariance_intercept_slope << "\t" <<
                fp[i]->chi2 << "\t" <<
                fp[i]->ndf << "\t" <<
                fp[i]->prob << "\n";
        }

        outfile.close();
        printf("file %s has been created\n", outfilename.Data());
        return;
    }
} // anon. namespace








void xn_xf_plots_per_run(const int CHOSEN_DETECTOR){
    // Get detector geometry and XN-XF correction
    int colDet, rowDet, nDet;
    printf("======================= loading parameters files .... \n");
    LoadDetectorGeometry(colDet, rowDet, nDet);
    printf("----- loading xf-xn correction.");
    double* xnCorr = LoadXNXFCorrection(nDet);

    // Create canvas
    gStyle->SetOptStat(11111111);
    TCanvas* c = new TCanvas(Form("DETECTOR %d", CHOSEN_DETECTOR), Form("DETECTOR %d", CHOSEN_DETECTOR), 2800, 1800);
    c->Divide(5,3); // HACK: this is hardcoded and may need to be changed!
    FitParameters** fp = new FitParameters*[kRunNumbers.size()];

    for (std::size_t i = 0; i < kRunNumbers.size(); ++i){
        const int run = kRunNumbers[i];
        // Generate filename
        std::string filename = "../../root_data/gen_run0";
        if (run < 10){
            filename += "0";
        }
        filename += std::to_string(run) + ".root";

        // Open the file
        TFile* f = new TFile(filename.data());
        TTree* t = (TTree*)f->Get("gen_tree");

        TString name = Form("r%d.d%d", run, CHOSEN_DETECTOR);
        TH2F* d = new TH2F(name, name, 200, 0, 5000, 200, 0, 5000);
        TString expr = Form("e[%d]:xn[%d]+xf[%d]>>%s", CHOSEN_DETECTOR, CHOSEN_DETECTOR, CHOSEN_DETECTOR, name.Data());
        TString cut = Form("!TMath::IsNaN(xn[%d]) && !TMath::IsNaN(xf[%d])", CHOSEN_DETECTOR, CHOSEN_DETECTOR);

        // Plot the data
        c->cd(i+1);
        t->Draw(expr, cut);
        d->SetDirectory(0);
        FitProfileToCorrectionHistograms(c, i, d, fp[i]);

        // Close the file
        if (f->IsOpen())f->Close();
    }

    c->Modified();
    c->Update();

    // Plot fit curves
    SaveFitToFile(kRunNumbers.size(), CHOSEN_DETECTOR, fp);
    TString outfilename = Form("xnxf_detector_%d.pdf", CHOSEN_DETECTOR);
    c->SaveAs(outfilename);

    // Delete fit parameters
    for (std::size_t i = 0; i < kRunNumbers.size(); ++i){
        delete fp[i];
    }
    delete[] fp;

}

#include "Monitors.h"

#include "TBox.h"
#include "TCanvas.h"
#include "TCutG.h"
#include "TDatime.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMacro.h"
#include "TMath.h"
#include "TMD5.h"
#include "TMultiGraph.h"
#include "TObjArray.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"


#include <fstream>
#include <vector>

#include "../Cleopatra/Isotope.h"
#include "GeneralSortMapping.h"
#include "MonitorsGlobals.h"

//############################################ User setting

void Monitors::Init(TTree *tree){

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runID", &runID, &b_runID);
   fChain->SetBranchAddress("e", e, &b_Energy);
   fChain->SetBranchAddress("e_t", e_t, &b_EnergyTimestamp);
   fChain->SetBranchAddress("xf", xf, &b_XF);
   fChain->SetBranchAddress("xf_t", xf_t, &b_XFTimestamp);
   fChain->SetBranchAddress("xn", xn, &b_XN);
   fChain->SetBranchAddress("xn_t", xn_t, &b_XNTimestamp);
   fChain->SetBranchAddress("ring", ring, &b_Ring);
   fChain->SetBranchAddress("ring_t", ring_t, &b_RingTimestamp);

   TBranch * br = (TBranch *) fChain->GetListOfBranches()->FindObject("rdt");
   if( br == NULL ){
      printf(" ++++++++ no Recoil.\n");
      isRDTExist = false;
   }else{
      isRDTExist = true;
      fChain->SetBranchAddress("rdt"  , rdt,   &b_RDT);
      fChain->SetBranchAddress("rdt_t", rdt_t, &b_RDTTimestamp);
   }
   
   br = (TBranch *) fChain->GetListOfBranches()->FindObject("crdt");
   if( br == NULL ){
      printf(" ++++++++ no Circular Recoil.\n");
      isCRDTExist = false;
   }else{
      isCRDTExist = true;
      fChain->SetBranchAddress("crdt"  , crdt,   &b_CRDT);
      fChain->SetBranchAddress("crdt_t", crdt_t, &b_CRDTTimestamp);
   }

   br = (TBranch *) fChain->GetListOfBranches()->FindObject("apollo");
   if( br == NULL ){
      printf(" ++++++++ no APOLLO.\n");
      isAPOLLOExist = false;
 }else{
      isAPOLLOExist = true;
      fChain->SetBranchAddress("apollo"  , apollo,   &b_APOLLO);
      fChain->SetBranchAddress("apollo_t", apollo_t, &b_APOLLOTimestamp);
   }
   
   br = (TBranch *) fChain->GetListOfBranches()->FindObject("tac");
   if( br == NULL ){
      printf(" ++++++++ no TAC.\n");
      isTACExist = false;
   }else{
      isTACExist = true;
      fChain->SetBranchAddress("tac", tac, &b_TAC);
      fChain->SetBranchAddress("tac_t", tac_t, &b_TACTimestamp);
   }

   br = (TBranch *) fChain->GetListOfBranches()->FindObject("elum");
   if( br == NULL ){
      printf(" ++++++++ no ELUM.\n");
      isELUMExist = false;
   }else{
      isELUMExist = true;
      fChain->SetBranchAddress("elum", elum, &b_ELUM);
      fChain->SetBranchAddress("elum_t", elum_t, &b_ELUMTimestamp);
   }

   br = (TBranch *) fChain->GetListOfBranches()->FindObject("ezero");
   if( br == NULL ){
      printf(" ++++++++ no EZERO.\n");
      isEZEROExist = false;
   }else{
      isEZEROExist = true;
      fChain->SetBranchAddress("ezero", ezero, &b_EZERO);
      fChain->SetBranchAddress("ezero_t", ezero_t, &b_EZEROTimestamp);
   }
   
   br = (TBranch *) fChain->GetListOfBranches()->FindObject("te");
   if( br == NULL ){
      printf(" ++++++++ no Array trace.\n");
      isArrayTraceExist = false;
   }else{
      isArrayTraceExist = true;
      if( monitorsglobals::isUseArrayTrace ){
         fChain->SetBranchAddress("te", e, &b_Trace_Energy); // replace  e with te
         printf("************ using Trace in array \n");
      }else{
         fChain->SetBranchAddress("e", e, &b_Energy);
      }
      fChain->SetBranchAddress("te_r", te_r, &b_Trace_Energy_RiseTime);
      fChain->SetBranchAddress("te_t", te_t, &b_Trace_Energy_Time);
   }
   
   br = (TBranch *) fChain->GetListOfBranches()->FindObject("trdt");
   if( br == NULL ){
      printf(" ++++++++ no Recoil trace.\n");
      isRDTTraceExist = false;
   }else{
      isRDTTraceExist = true;
      if( monitorsglobals::isUseRDTTrace ) {
         fChain->SetBranchAddress("trdt", rdt, &b_Trace_RDT); // replace  rdt with trdt
         printf("************ using Trace in recoil \n");
      }else{
         fChain->SetBranchAddress("trdt", rdt, &b_RDT);
      }
      fChain->SetBranchAddress("trdt_t", trdt_t, &b_Trace_RDT_Time);
      fChain->SetBranchAddress("trdt_r", trdt_r, &b_Trace_RDT_RiseTime);
   }
   
   startTime = 0;
   endTime = 0;
   
   printf("=================================== End of Branch Pointer Inititization. \n");
}

Bool_t Monitors::Notify(){
   return kTRUE;
}

TString Monitors::FindStartEndTime(TTree * tree, TString BranchName){
   
   ULong64_t timeStamp[100];
   
   int totEnetries = tree->GetEntries();
   
   printf(">>>>>>> Number of Entries : %d, Check Branch : %s\n", totEnetries, BranchName.Data());
   
   tree->SetBranchAddress(BranchName, timeStamp);
   for( int ev = 0; ev < totEnetries; ev ++ ){
      tree->GetEntry(ev);
      bool breakFlag = false;
      for( int id = 0; id < NARRAY; id ++){
         if( timeStamp[id] > 0 ) {
            startTime = timeStamp[id];
            breakFlag = true;
            break;
         }
      }
      if( breakFlag ) break;
   }
   
   for( int ev =  totEnetries - 1 ; ev > 0 ;  ev -- ){
      tree->GetEntry(ev);
      bool breakFlag = false;
      for( int id = 0; id < NARRAY; id ++){
         if( timeStamp[id] > 0 ) {
            endTime = timeStamp[id];
            breakFlag = true;
            break;
         }
      }
      if( breakFlag ) break;
   }
   
   return BranchName;
}


void DrawLine(TH1 * hist, double pos, int colorCode = 2){
   
   double yMax = hist->GetMaximum();
   TLine * line = new TLine(pos, 0, pos, yMax);
   line->SetLineColor(colorCode);
   line->Draw("");
   
}

void DrawBox(TH1* hist, double x1, double x2, Color_t color, float alpha){

   double yMax = hist->GetMaximum();
   TBox * box = new TBox (x1, 0, x2, yMax);
   box->SetFillColorAlpha(color, alpha);
   box->Draw();

}

void Monitors::Draw2DHist(TH2F * hist, int colzThreshold){
   
   if( hist->Integral() < colzThreshold ){
      hist->SetMarkerStyle(20);
      hist->SetMarkerSize(0.3);
      hist->Draw("box");
   }else{
      hist->Draw("colz"); 
   }
}

void Monitors::SlaveBegin(TTree * /*tree*/){
   /// not use, if use, place in Monitors.C
   TString option = GetOption();
}


void Monitors::SlaveTerminate(){
  /// not use, if use, place in Monitors.C
}


double solid_angle( double th ) {
   
   double ang = 2. * TMath::Pi() * TMath::Sin( th*TMath::DegToRad() );
   return 1.0 / ( ang * TMath::RadToDeg() );
   
}

void Monitors::LoadDetGeoAndReactionConfigFile(){

  string detGeoFileName = "detectorGeo.txt";
  printf(" loading detector geometery : %s.", detGeoFileName.c_str());
  TMacro * haha = new TMacro();
  if( haha->ReadFile(detGeoFileName.c_str()) > 0 ) {
    detGeo = LoadDetectorGeo(haha);
    printf("... done.\n");
    PrintDetGeo(detGeo);
    monitorsglobals::zRange[0] = detGeo.zMin - 50 ;
    monitorsglobals::zRange[1] = detGeo.zMax + 50 ;

    numCol = detGeo.nDet;
    numRow = detGeo.mDet;
  }else{
    printf("... fail\n");
  }

  string reactionConfigFileName = "reactionConfig.txt";
  printf(" loading reaction config : %s.", reactionConfigFileName.c_str());
  TMacro * kaka = new TMacro();
  if( kaka->ReadFile(reactionConfigFileName.c_str()) > 0 ) {
    reactionConfig  = LoadReactionConfig(kaka);
    printf("..... done.\n");
    PrintReactionConfig(reactionConfig);
  }else{
    printf("..... fail\n");
  }

}

void Monitors::LoadXFXNCorr(){
   //========================================= xf = xn correction
   printf(" loading xf-xn correction.");
   ifstream file;
   file.open("correction_xf_xn.dat");
   if( file.is_open() ){
      double a;
      int i = 0;
      while( file >> a ){
         if( i >= NARRAY) break;
         monitorsglobals::xnCorr[i] = a;
         i = i + 1;
      }
      
      printf(".......... done.\n");
   }else{
      printf(".......... fail.\n");
      
      for(int i = 0; i < NARRAY; i++){
         monitorsglobals::xnCorr[i] = 1;
      }
   }
   file.close();
   
}

void Monitors::LoadXFXN2ECorr(){
   //========================================= e = xf + xn correction
   printf(" loading xf/xn-e correction.");
   ifstream file;
   file.open("correction_xfxn_e.dat");
   if( file.is_open() ){
      double a, b;
      int i = 0;
      while( file >> a >> b){
         if( i >= NARRAY) break;
         monitorsglobals::xfxneCorr[i][0] = a;
         monitorsglobals::xfxneCorr[i][1] = b;
         i = i + 1;
      }
      printf("........ done.\n");
   }else{
      printf("........ fail.\n");
      for(int i = 0; i < NARRAY; i++){
         monitorsglobals::xfxneCorr[i][0] = 0;
         monitorsglobals::xfxneCorr[i][1] = 1;
      }
   }
   file.close();
} 

void Monitors::LoadXScaleCorr(){
   //========================================= X-Scale correction
   printf(" loading x-Scale correction.");
   ifstream file;
   file.open("correction_scaleX.dat");
   if( file.is_open() ){
      double a, b;
      int i = 0;
      while( file >> a ){
         if( i >= NARRAY) break;
         monitorsglobals::xScale[i] = a;  
         i = i + 1;
      }
      printf("........ done.\n");
      
   }else{
      printf("........ fail.\n");
      for( int i = 0; i < NARRAY ; i++){
         monitorsglobals::xScale[i] = 1.;
      }
   }
   file.close();
   
}


void Monitors::LoadECorr(){
   //========================================= e correction
   printf(" loading e correction.");
   ifstream file;
   file.open("correction_e.dat");
   if( file.is_open() ){
      double a, b;
      int i = 0;
      while( file >> a >> b){
         if( i >= NARRAY) break;
         monitorsglobals::eCorr[i][0] = a;  // 1/a1
         monitorsglobals::eCorr[i][1] = b;  //  a0 , e' = e * a1 + a0
         //printf("\n%2d, e0: %9.4f, e1: %9.4f", i, eCorr[i][0], eCorr[i][1]);
         i = i + 1;
      }
      printf(".............. done.\n");
      
   }else{
      printf(".............. fail.\n");
      for( int i = 0; i < NARRAY ; i++){
         monitorsglobals::eCorr[i][0] = 1.;
         monitorsglobals::eCorr[i][1] = 0.;
      }
   }
   file.close();

}

void Monitors::LoadRDTCorr(){
   //========================================= e correction
   printf(" loading rdt correction.");
   ifstream file;
   file.open("correction_rdt.dat");
   if( file.is_open() ){
      double a, b;
      int i = 0;
      while( file >> a >> b){
         if( i >= NRDT ) break;
         monitorsglobals::rdtCorr[i][0] = a;  //  a1
         monitorsglobals::rdtCorr[i][1] = b;  //  a0 , e' = e * a1 + a0
         i = i + 1;
      }
      printf("............ done.\n");
      
   }else{
      printf("............ fail.\n");
      for( int i = 0; i < NRDT ; i++){
         monitorsglobals::rdtCorr[i][0] = 1.;
         monitorsglobals::rdtCorr[i][1] = 0.;
      }
   }
   file.close();
   
}

void Monitors::LoadReactionPars(){
   
   //========================================= reaction parameters
   //check is the transfer.root is using the latest reactionConfig.txt   
   //sicne reaction.dat is generated as a by-product of transfer.root
   TFile * transfer = new TFile("transfer.root");
   TString aaa1 = "";
   TString aaa2 = "";
   if( transfer->IsOpen() ){
      TMacro * reactionConfig = (TMacro *) transfer->FindObjectAny("reactionConfig");
      TMacro presentReactionConfig ("reactionConfig.txt");
      aaa1 = ((TMD5*) reactionConfig->Checksum())->AsString();
      aaa2 = ((TMD5*) presentReactionConfig.Checksum())->AsString();
   }
   //printf("%s\n", aaa1.Data());
   //printf("%s\n", aaa2.Data());

   //if( aaa1 != aaa2 ) {
   //  printf("########################## recalculate transfer.root \n");
   //  system("../Cleopatra/Transfer");
   //  printf("########################## transfer.root updated\n");
   //}
   printf(" loading reaction parameters");
   ifstream file;
   file.open("reaction.dat");
   monitorsglobals::isReaction = false;
   if( file.is_open() ){
      string x;
      int i = 0;
      while( file >> x ){
         if( x.substr(0,2) == "//" )  continue;
         if( i == 0 ) monitorsglobals::mass = atof(x.c_str());
         if( i == 1 ) monitorsglobals::q    = atof(x.c_str());
         if( i == 2 ) monitorsglobals::betRel = atof(x.c_str()); 
         if( i == 3 ) monitorsglobals::Et   = atof(x.c_str()); 
         if( i == 4 ) monitorsglobals::massB = atof(x.c_str()); 
         i = i + 1;
      }
      printf("........ done.\n");

      monitorsglobals::isReaction = true;
      monitorsglobals::alpha = 299.792458 * abs(detGeo.Bfield) * monitorsglobals::q / TMath::TwoPi()/1000.; //MeV/mm
      monitorsglobals::gamm = 1./TMath::Sqrt(1-monitorsglobals::betRel*monitorsglobals::betRel);
      monitorsglobals::G = monitorsglobals::alpha * monitorsglobals::gamm * monitorsglobals::betRel * detGeo.detPerpDist ;
      printf("\tmass-b    : %f MeV/c2 \n", monitorsglobals::mass);
      printf("\tcharge-b  : %f \n", monitorsglobals::q);
      printf("\tE-total   : %f MeV \n", monitorsglobals::Et);
      printf("\tmass-B    : %f MeV/c2 \n", monitorsglobals::massB);
      printf("\tbeta      : %f \n", monitorsglobals::betRel);
      printf("\tB-field   : %f T \n", detGeo.Bfield);
      printf("\tslope     : %f MeV/mm \n", monitorsglobals::alpha * monitorsglobals::betRel);
      printf("\tdet radius: %f mm \n", detGeo.detPerpDist);
      printf("\tG-coeff   : %f MeV \n", monitorsglobals::G);
      printf("=================================\n");

   }else{
      printf("........ fail.\n");
      monitorsglobals::isReaction = false;
   }
   file.close();
   
}


/*###########################################################
 * Plotting Function
###########################################################*/

void Monitors::PlotEZ(bool isRaw){
   monitorsglobals::padID++; monitorsglobals::cCanvas->cd(monitorsglobals::padID);

   if( isRaw ) {
      Draw2DHist(monitorsglobals::heCalVz);
      monitorsglobals::heCalVz->SetTitle("E vs Z | " + monitorsglobals::canvasTitle + " | " + monitorsglobals::rdtCutFile1);
      if( monitorsglobals::skipDetID.size() > 0 ){
         monitorsglobals::text.DrawLatex(0.15, 0.3, "skipped Detector:");
         for( int i = 0; i < (int) monitorsglobals::skipDetID.size(); i++){
            monitorsglobals::text.DrawLatex(0.15 + 0.1*i, 0.25, Form("%d", monitorsglobals::skipDetID[i]));
         }
      }

      monitorsglobals::text.DrawLatex(0.15, 0.8, Form("%.1f < eCal < %.1f MeV", monitorsglobals::eCalCut[0], monitorsglobals::eCalCut[1]));
      if( monitorsglobals::xGate < 1 ) monitorsglobals::text.DrawLatex(0.15, 0.75, Form("with |x-0.5|<%.4f", monitorsglobals::xGate/2.));

   }else{
      Draw2DHist(monitorsglobals::heCalVzGC);
   
      if( monitorsglobals::isCutFileOpen1 ) monitorsglobals::text.DrawLatex(0.15, 0.8, "with Recoil gate");
      if( monitorsglobals::xGate < 1 )   monitorsglobals::text.DrawLatex(0.15, 0.75, Form("with |x-0.5|<%.4f", monitorsglobals::xGate/2.));
      if( monitorsglobals::isTACGate )   monitorsglobals::text.DrawLatex(0.15, 0.7, Form("%d < TAC < %d", monitorsglobals::tacGate[0], monitorsglobals::tacGate[1]));
      if( monitorsglobals::isTimeGateOn) monitorsglobals::text.DrawLatex(0.15, 0.7, Form("%d < coinTime < %d", monitorsglobals::timeGate[0], monitorsglobals::timeGate[1]));

   }
   
   TFile * transfer = new TFile("transfer.root");
   TObjArray * gList  = NULL ;
   TObjArray * fxList = NULL ;
   int nGList = 0;
   int nFxList = 0;
   if( transfer->IsOpen() ) {
     gList  = (TObjArray *) transfer->FindObjectAny("gList");
     nGList = gList->GetLast() + 1;
     fxList = (TObjArray *) transfer->FindObjectAny("fxList");
     nFxList = fxList->GetLast() +1 ;
   }
   
   ///the constant thetaCM line
   if( transfer->IsOpen() ) gList->At(0)->Draw("same");
   ///the e-z line for excitation 
   if( transfer->IsOpen() ){
      for( int i = 0 ; i < nFxList ; i++){
         ((TF1*)fxList->At(i))->SetLineColor(6);
         fxList->At(i)->Draw("same");
      }
   }

}

void Monitors::PlotTDiff(bool isGated, bool isLog){
   monitorsglobals::padID++; monitorsglobals::cCanvas->cd(monitorsglobals::padID); 
   if( isLog ) monitorsglobals::cCanvas->cd(monitorsglobals::padID)->SetLogy(1);
   double yMax = 0;
   if( isGated ){
      yMax = monitorsglobals::htdiff->GetMaximum()*1.2;
      if( isLog ){
         monitorsglobals::htdiff->GetYaxis()->SetRangeUser(1, yMax);
      }else{
         monitorsglobals::htdiff->GetYaxis()->SetRangeUser(0, yMax);
      }
   }
   monitorsglobals::htdiff->Draw();
   if( isGated ){
      monitorsglobals::htdiffg->SetLineColor(2);
      monitorsglobals::htdiffg->Draw("same");
   }
   
   if( monitorsglobals::isCutFileOpen1 ) monitorsglobals::text.DrawLatex(0.15, 0.8, "with Recoil gate");
   if(monitorsglobals::isTimeGateOn)monitorsglobals::text.DrawLatex(0.15, 0.7, Form("%d < coinTime < %d", monitorsglobals::timeGate[0], monitorsglobals::timeGate[1])); 
   DrawBox(monitorsglobals::htdiff, monitorsglobals::timeGate[0], monitorsglobals::timeGate[1], kGreen, 0.0);
}

void Monitors::PlotRDT(int id, bool isRaw){
   monitorsglobals::padID++; monitorsglobals::cCanvas->cd(monitorsglobals::padID);

   if( isRaw ){
      Draw2DHist(monitorsglobals::hrdt2D[id]);
   }else{
      Draw2DHist(monitorsglobals::hrdt2Dg[id]);
   }
   if(monitorsglobals::isTimeGateOn)monitorsglobals::text.DrawLatex(0.15, 0.8, Form("%d < coinTime < %d", monitorsglobals::timeGate[0], monitorsglobals::timeGate[1])); 
   if( monitorsglobals::isTACGate ) monitorsglobals::text.DrawLatex(0.15, 0.7, Form("%d < TAC < %d", monitorsglobals::tacGate[0], monitorsglobals::tacGate[1]));
   if( monitorsglobals::isCutFileOpen1 && monitorsglobals::numCut1 > id ) {monitorsglobals::cutG = (TCutG *)monitorsglobals::cutList1->At(id) ; monitorsglobals::cutG->Draw("same");}
   if( monitorsglobals::isCutFileOpen2 && monitorsglobals::numCut2 > id ) {monitorsglobals::cutG = (TCutG *)monitorsglobals::cutList2->At(id) ; monitorsglobals::cutG->Draw("same");}

}

void Monitors::PlotCRDTPolar(){
  monitorsglobals::padID++; monitorsglobals::cCanvas->cd(monitorsglobals::padID);
  monitorsglobals::cCanvas->cd(monitorsglobals::padID)->DrawFrame(-50, -50, 50, 50);
  monitorsglobals::hcrdtPolar->Draw("same colz pol");
}

/***************************
 ***************************/
void Monitors::Begin(TTree *tree)
{
   TString option = GetOption();
   monitorsglobals::NumEntries = tree->GetEntries();
   
   monitorsglobals::canvasTitle.Form("Runs: ");
   monitorsglobals::lastRunID = -1;
   monitorsglobals::contFlag = false;
   
   printf("####################################################\n");
   printf("##########           Monitors.C            #########\n");
   printf("####################################################\n");
   
   printf("printControlID = %d \n", printControlID);
   
   //===================================================== loading parameter
   printf("################## loading parameter files\n"); 
   
   LoadDetGeoAndReactionConfigFile();
   LoadXFXNCorr();
   LoadXFXN2ECorr();
   LoadXScaleCorr();
   LoadECorr();
   LoadRDTCorr();
   LoadReactionPars();

   //================  Get Recoil cuts;
   monitorsglobals::cutG = new TCutG();

   TFile * fCut1 = new TFile(monitorsglobals::rdtCutFile1);
   monitorsglobals::isCutFileOpen1 = fCut1->IsOpen();
   if(!monitorsglobals::isCutFileOpen1) {
      printf( "Failed to open rdt-cutfile 1 : %s\n" , monitorsglobals::rdtCutFile1.Data());
      monitorsglobals::rdtCutFile1 = "";
   }
   monitorsglobals::numCut1 = 0 ;

   if( monitorsglobals::isCutFileOpen1 ){
      monitorsglobals::cutList1 = (TObjArray *) fCut1->FindObjectAny("cutList");
      monitorsglobals::numCut1 = monitorsglobals::cutList1->GetEntries();
      printf("=========== found %d cutG in %s \n", monitorsglobals::numCut1, fCut1->GetName());

      for(int i = 0; i < monitorsglobals::numCut1 ; i++){
         printf("cut name : %s , VarX: %s, VarY: %s, numPoints: %d \n",
            monitorsglobals::cutList1->At(i)->GetName(),
            ((TCutG*)monitorsglobals::cutList1->At(i))->GetVarX(),
            ((TCutG*)monitorsglobals::cutList1->At(i))->GetVarY(),
            ((TCutG*)monitorsglobals::cutList1->At(i))->GetN());
      }
   }
   
   TFile * fCut2 = new TFile(monitorsglobals::rdtCutFile2);
   monitorsglobals::isCutFileOpen2 = fCut2->IsOpen();
   if(!monitorsglobals::isCutFileOpen2) {
      printf( "Failed to open rdt-cutfile 2 : %s\n" , monitorsglobals::rdtCutFile2.Data());
      monitorsglobals::rdtCutFile2 = "";
   }
   monitorsglobals::numCut2 = 0 ;
   if( monitorsglobals::isCutFileOpen2 ){
      monitorsglobals::cutList2 = (TObjArray *) fCut2->FindObjectAny("cutList");
      monitorsglobals::numCut2 = monitorsglobals::cutList2->GetEntries();
      printf("=========== found %d cutG in %s \n", monitorsglobals::numCut2, fCut2->GetName());

      for(int i = 0; i < monitorsglobals::numCut2 ; i++){
         printf("cut name : %s , VarX: %s, VarY: %s, numPoints: %d \n",
            monitorsglobals::cutList2->At(i)->GetName(),
            ((TCutG*)monitorsglobals::cutList2->At(i))->GetVarX(),
            ((TCutG*)monitorsglobals::cutList2->At(i))->GetVarY(),
            ((TCutG*)monitorsglobals::cutList2->At(i))->GetN());
      }
   }
   
   //================  Get EZ cuts;
   TFile * fCutez = new TFile(monitorsglobals::ezCutFile);
   monitorsglobals::isEZCutFileOpen = fCutez->IsOpen(); 
   if( !monitorsglobals::isEZCutFileOpen)  printf( "Failed to open E-Z cutfile : %s\n" , monitorsglobals::ezCutFile.Data());
   if( monitorsglobals::isEZCutFileOpen ) {
      
      monitorsglobals::EZCut = (TCutG *) fCutez->FindObjectAny("cutEZ");
      if( monitorsglobals::EZCut != NULL ) {
         printf("Found EZ cut| name : %s, VarX: %s, VarY: %s, numPoints: %d \n", 
                  monitorsglobals::EZCut->GetName(),
                  monitorsglobals::EZCut->GetVarX(),
                  monitorsglobals::EZCut->GetVarY(),
                  monitorsglobals::EZCut->GetN()
            );
      }
   }
   
   //=============== Get Time stamp range
   //TODO use the info macro
   TString branch = FindStartEndTime(tree, "e_t");
   if( (startTime == 0 && endTime == 0 ) || (startTime == endTime)  ) branch = FindStartEndTime(tree, "rdt_t");
   if( (startTime == 0 && endTime == 0 ) || (startTime == endTime)  ) branch = FindStartEndTime(tree, "ezero_t");
   if( (startTime == 0 && endTime == 0 ) || (startTime == endTime)  ) branch = FindStartEndTime(tree, "elum_t");
   
   monitorsglobals::runTime = (endTime-startTime)/1e8/60.;
   
   printf("================ Time Range in %s timestamp branch, total Entries : %lld \n", branch.Data(),  tree->GetEntries());
   printf("start time : %llu = %.1f sec = %.0f min\n", startTime, startTime/1e8, startTime/1e8/60.);
   printf("  end time : %llu = %.1f sec = %.0f min\n",  endTime, endTime/1e8, endTime/1e8/60.);
   printf("  Duration : %llu = %.1f sec = %.0f min\n",  endTime - startTime, (endTime-startTime)/1e8, monitorsglobals::runTime);
   
   int timeRange[2];
   timeRange[0] = TMath::Max((int) TMath::Floor(startTime/1e8/60.), monitorsglobals::timeRangeUser[0]);
   timeRange[1] = TMath::Min((int) TMath::Ceil(endTime/1e8/60.)   , monitorsglobals::timeRangeUser[1]);
   
   int dTime = timeRange[1] - timeRange[0];
   
   if( dTime > 0 ){
      timeRange[0] = timeRange[0] - dTime*0.1;
      timeRange[1] = timeRange[1] + dTime*0.1;  
   }else{
      dTime = 1000;
      timeRange[0] = 0;
      timeRange[1] = 1000;
   }
   
   //========================= Generate all of the histograms needed for drawing later on
   printf("======================================== Histograms declaration\n");

   gROOT->cd();
   
   monitorsglobals::heVID    = new TH2F("heVID",    "Raw e vs channel", NARRAY, 0, NARRAY, 500, monitorsglobals::rawEnergyRange[0], monitorsglobals::rawEnergyRange[1]);
   monitorsglobals::heVIDG   = new TH2F("heVIDG",   "Raw e vs channel (gated)", NARRAY, 0, NARRAY, 500, monitorsglobals::rawEnergyRange[0], monitorsglobals::rawEnergyRange[1]);
   monitorsglobals::hringVID = new TH2F("hringVID", "Raw Ring vs channel", NARRAY, 0, NARRAY, 500, monitorsglobals::rawEnergyRange[0], monitorsglobals::rawEnergyRange[1]);
   monitorsglobals::hxfVID   = new TH2F("hxfVID",   "Raw xf vs channel", NARRAY, 0, NARRAY, 500, monitorsglobals::rawEnergyRange[0], monitorsglobals::rawEnergyRange[1]);
   monitorsglobals::hxnVID   = new TH2F("hxnVID",   "Raw xn vs channel", NARRAY, 0, NARRAY, 500, monitorsglobals::rawEnergyRange[0], monitorsglobals::rawEnergyRange[1]);

   for (Int_t i = 0; i < NARRAY; i++) {//array loop
       monitorsglobals::he[i]      = new TH1F(Form("he%d", i),     Form("Raw e (ch=%d); e (channel); count", i),                    2000, monitorsglobals::rawEnergyRange[0], monitorsglobals::rawEnergyRange[1]);
       monitorsglobals::hring[i]   = new TH1F(Form("hring%d", i),  Form("Raw ring (ch=%d); ring (channel); count", i),              200, monitorsglobals::rawEnergyRange[0], monitorsglobals::rawEnergyRange[1]);
       monitorsglobals::hxf[i]     = new TH1F(Form("hxf%d", i),    Form("Raw xf (ch=%d); xf (channel); count", i),                  200, monitorsglobals::rawEnergyRange[0], monitorsglobals::rawEnergyRange[1]);
       monitorsglobals::hxn[i]     = new TH1F(Form("hxn%d", i),    Form("Raw xn (ch=%d); xn (channel); count", i),                  200, monitorsglobals::rawEnergyRange[0], monitorsglobals::rawEnergyRange[1]);
       monitorsglobals::heVring[i] = new TH2F(Form("heVring%d",i), Form("Raw e vs. ring (ch=%d);ring (channel);e (channel)",i),     500, monitorsglobals::rawEnergyRange[0], monitorsglobals::rawEnergyRange[1], 500, monitorsglobals::rawEnergyRange[0], monitorsglobals::rawEnergyRange[1]);
       monitorsglobals::hxfVxn[i]  = new TH2F(Form("hxfVxn%d",i),  Form("Raw xf vs. xn (ch=%d);xf (channel);xn (channel)",i),       500,                 0, monitorsglobals::rawEnergyRange[1], 500,                 0, monitorsglobals::rawEnergyRange[1]);
       monitorsglobals::heVxs[i]   = new TH2F(Form("heVxs%d", i),  Form("Raw e vs xf+xn (ch=%d); xf+xn (channel); e (channel)", i), 500, monitorsglobals::rawEnergyRange[0], monitorsglobals::rawEnergyRange[1], 500, monitorsglobals::rawEnergyRange[0], monitorsglobals::rawEnergyRange[1]);
       monitorsglobals::heCal[i]        = new TH1F(Form("heCal%d", i),       Form("Corrected e (ch=%d); e (MeV); count", i),                                         2000,   monitorsglobals::energyRange[0], monitorsglobals::energyRange[1]);
       monitorsglobals::hxfCalVxnCal[i] = new TH2F(Form("hxfCalVxnCal%d",i), Form("Corrected XF vs. XN (ch=%d);XF (channel);XN (channel)",i),                         500,                0, monitorsglobals::rawEnergyRange[1], 500,                 0, monitorsglobals::rawEnergyRange[1]);      
       monitorsglobals::heVxsCal[i]     = new TH2F(Form("heVxsCal%d", i),    Form("Raw e vs Corrected xf+xn (ch=%d); corrected xf+xn (channel); Raw e (channel)", i), 500,monitorsglobals::rawEnergyRange[0], monitorsglobals::rawEnergyRange[1], 500, monitorsglobals::rawEnergyRange[0], monitorsglobals::rawEnergyRange[1]);           
       monitorsglobals::heVx[i]         = new TH2F(Form("heVx%d",i),         Form("Raw PSD E vs. X (ch=%d);X (channel);E (channel)",i),                               500,             -0.1,               1.1, 500, monitorsglobals::rawEnergyRange[0], monitorsglobals::rawEnergyRange[1]);
       monitorsglobals::hringVx[i]      = new TH2F(Form("hringVx%d",i),      Form("Ring vs, X (ch=%d);X (channel);Ring (channel)",i),                                 500,             -0.1,               1.1, 500, monitorsglobals::rawEnergyRange[0], monitorsglobals::rawEnergyRange[1]);
       monitorsglobals::heCalVxCal[i]   = new TH2F(Form("heCalVxCal%d",i),   Form("Cal PSD E vs. X (ch=%d);X (cm);E (MeV)",i),                                        500,             -2.5,  detGeo.detLength + 2.5, 500,    monitorsglobals::energyRange[0],    monitorsglobals::energyRange[1]);
       monitorsglobals::heCalVxCalG[i]  = new TH2F(Form("heCalVxCalG%d",i),  Form("Cal PSD E vs. X (ch=%d);X (cm);E (MeV)",i),                                        500,             -2.5,  detGeo.detLength + 2.5, 500,    monitorsglobals::energyRange[0],    monitorsglobals::energyRange[1]);
       monitorsglobals::hExVxCal[i]     = new TH2F(Form("hExVxCal%d",i),     Form("Ex vs X (ch=%d); X (cm); Ex (MeV)", i),                                            500,             -0.1,               1.1, (int) (monitorsglobals::exRange[2]-monitorsglobals::exRange[1])/monitorsglobals::exRange[0]*1000, monitorsglobals::exRange[1], monitorsglobals::exRange[2]);
   }
   
   monitorsglobals::heCalID = new TH2F("heCalID", "Corrected E vs detID; detID; E / 10 keV", NARRAY, 0, NARRAY, 2000, monitorsglobals::energyRange[0], monitorsglobals::energyRange[1]);
   
   //====================== E-Z plot
   monitorsglobals::heCalVz   = new TH2F("heCalVz",  "E vs. Z;Z (mm);E (MeV)"      , 400, monitorsglobals::zRange[0], monitorsglobals::zRange[1], 400, monitorsglobals::energyRange[0], monitorsglobals::energyRange[1]);
   monitorsglobals::heCalVzGC = new TH2F("heCalVzGC","E vs. Z gated;Z (mm);E (MeV)", 400, monitorsglobals::zRange[0], monitorsglobals::zRange[1], 400, 0, monitorsglobals::energyRange[1]);
   
   for( int i = 0; i < numRow; i++){
      monitorsglobals::hecalVzRow[i] = new TH2F(Form("heCalVzRow%d", i), Form("E vs. Z (ch=%d-%d); Z (cm); E (MeV)", numCol*i, numCol*(i+1)-1), 500, monitorsglobals::zRange[0], monitorsglobals::zRange[1], 500, monitorsglobals::energyRange[0], monitorsglobals::energyRange[1]);
      monitorsglobals::hecalVzRowG[i] = new TH2F(Form("heCalVzRow%dG", i), Form("E vs. Z (ch=%d-%d); Z (cm); E (MeV)", numCol*i, numCol*(i+1)-1), 500, monitorsglobals::zRange[0], monitorsglobals::zRange[1], 500, monitorsglobals::energyRange[0], monitorsglobals::energyRange[1]);
   }

   //===================== Recoils
   for (Int_t i = 0; i < NRDT; i++) {
      if( i % 2 == 0 ) monitorsglobals::hrdt[i]  = new TH1F(Form("hrdt%d",i), Form("Raw Recoil E(ch=%d); E (channel)",i),         500, monitorsglobals::rdtERange[0], monitorsglobals::rdtERange[1]);
      if( i % 2 == 1 ) monitorsglobals::hrdt[i]  = new TH1F(Form("hrdt%d",i), Form("Raw Recoil DE(ch=%d); DE (channel)",i),       500,monitorsglobals::rdtDERange[0],monitorsglobals::rdtDERange[1]);

      if( i % 2 == 0 ) monitorsglobals::hrdtg[i] = new TH1F(Form("hrdt%dg",i),Form("Raw Recoil E(ch=%d) gated; E (channel)",i),   500, monitorsglobals::rdtERange[0], monitorsglobals::rdtERange[1]);
      if( i % 2 == 1 ) monitorsglobals::hrdtg[i] = new TH1F(Form("hrdt%dg",i),Form("Raw Recoil DE(ch=%d) gated; DE (channel)",i), 500,monitorsglobals::rdtDERange[0],monitorsglobals::rdtDERange[1]);
      
      ///dE vs E      
      if( i % 2 == 0 ) {
         int tempID = i / 2;
         monitorsglobals::hrdt2D[tempID]    = new TH2F(Form("hrdt2D%d",tempID),    Form("Raw Recoil DE vs Eres (dE=%d, E=%d); Eres (channel); DE (channel)", i+1, i),       500,monitorsglobals::rdtERange[0],monitorsglobals::rdtERange[1],500,monitorsglobals::rdtDERange[0],monitorsglobals::rdtDERange[1]);
         monitorsglobals::hrdt2Dsum[tempID] = new TH2F(Form("hrdt2Dsum%d",tempID), Form("Raw Recoil DE vs Eres+DE (dE=%d, E=%d); Eres+DE (channel); DE (channel)", i+1, i), 500,monitorsglobals::rdtERange[0],monitorsglobals::rdtERange[1]+monitorsglobals::rdtDERange[1],500,monitorsglobals::rdtDERange[0],monitorsglobals::rdtDERange[1]);
         monitorsglobals::hrdt2Dg[tempID]   = new TH2F(Form("hrdt2Dg%d",tempID),   Form("Gated Raw Recoil DE vs Eres (dE=%d, E=%d); Eres (channel); DE (channel)",i+1, i),  500,monitorsglobals::rdtERange[0],monitorsglobals::rdtERange[1],500,monitorsglobals::rdtDERange[0], monitorsglobals::rdtDERange[1]);
      }
   }
   monitorsglobals::hrdtID = new TH2F("hrdtID", "RDT vs ID; ID; energy [ch]", 8, 0, 8, 500, TMath::Min(monitorsglobals::rdtERange[0], monitorsglobals::rdtDERange[0]), TMath::Max(monitorsglobals::rdtERange[1], monitorsglobals::rdtDERange[1])); 
   
   monitorsglobals::hrdtMatrix = new TH2F("hrdtMatrix", "RDT ID vs RDT ID", 16 , 0, 8, 16, 0, 8);
   
   monitorsglobals::hrdtRate1 = new TH1F("hrdtRate1", "recoil rate 1 / min; min; count / 1 min", timeRange[1] - timeRange[0], timeRange[0], timeRange[1]);
   monitorsglobals::hrdtRate1->SetLineColor(2);
   
   monitorsglobals::hrdtRate2 = new TH1F("hrdtRate2", "recoil rate 2 / min; min; count / 1 min", timeRange[1] - timeRange[0], timeRange[0], timeRange[1]);
   monitorsglobals::hrdtRate2->SetLineColor(4);

   //===================== APOLLO
   for (Int_t i = 0 ; i < NAPOLLO; i++) {
      monitorsglobals::hApollo[i] = new TH1F(Form("hApollo%d",i),Form("Raw Apollo E(ch=%d); E (channel)",i), 250, monitorsglobals::apolloRange[0],monitorsglobals::apolloRange[1]);
   }
   
   //===================== Circular Recoil
   if( NCRDT > 0 ){
      monitorsglobals::hcrdtID = new TH2F("hcrdtID", "Circular Recoil ID; Angular ID; Radial ID;", NCRDT/2, 0, NCRDT/2, NCRDT/2, 0, NCRDT/2);
   }
   monitorsglobals::hcrdtPolar = new TH2F("hcrdtPolar", "Polar ID", 8, -TMath::Pi(), TMath::Pi(),8, 10, 50);

   for( int i = 0; i < NCRDT; i++){
       monitorsglobals::hcrdt[i] = new TH1F(Form("hcrdt%d", i), Form("Raw Circular Recoil-%d", i), 500, monitorsglobals::crdtRange[0], monitorsglobals::crdtRange[1] );
   }

   //===================== multiplicity
   monitorsglobals::hmult             = new TH2I("hmult","Array Multiplicity vs Recoil Multiplicity; Array ; Recoil",10,0,10,10,0,10);
   monitorsglobals::hmultEZ           = new TH1I("hmultEZ","Filled EZ with coinTime and recoil",10,0,10);
   monitorsglobals::hArrayRDTMatrix   = new TH2I("hArrayRDTMatrix", "Array ID vs Recoil ID; Array ID; Recoil ID",NARRAY,0,NARRAY,NRDT,0,NRDT);
   monitorsglobals::hArrayRDTMatrixG  = new TH2I("hArrayRDTMatrixG","Array ID vs Recoil ID / g; Array ID; Recoil ID",NARRAY,0,NARRAY,NRDT,0,NRDT);

   //===================== coincident time 
   monitorsglobals::htdiff  = new TH1I("htdiff" ,"Coincident time (recoil-dE - array); time [ch = 10ns]; count", monitorsglobals::coinTimeRange[1] - monitorsglobals::coinTimeRange[0], monitorsglobals::coinTimeRange[0], monitorsglobals::coinTimeRange[1]);   
   monitorsglobals::htdiffg = new TH1I("htdiffg","Coincident time (recoil-dE - array) w/ recoil gated; time [ch = 10ns]; count", monitorsglobals::coinTimeRange[1] - monitorsglobals::coinTimeRange[0], monitorsglobals::coinTimeRange[0], monitorsglobals::coinTimeRange[1]);
 
   //===================== TAC
    monitorsglobals::htac  = new TH1F("htac","Array-RF TAC; kind of time diff [a.u.]; Counts", monitorsglobals::TACRange[0], monitorsglobals::TACRange[1], monitorsglobals::TACRange[2]);
   monitorsglobals::htac2 = new TH1F("htac2","Array-RF TAC2; time diff [10 ns]; Counts", monitorsglobals::TAC2Range[0], monitorsglobals::TAC2Range[1], monitorsglobals::TAC2Range[2]);

   for (Int_t i = 0 ; i < NARRAY; i++) {
      monitorsglobals::htacArray[i] = new TH1I(Form("htacArray%d",i), Form("Array-RDT TAC for ch%d",i), 200, -100,100);
   }
   
   monitorsglobals::htacEx  = new TH2F("htacEx", "Ex - TAC ; TAC [a.u.]; Ex [MeV]", monitorsglobals::TACRange[0], monitorsglobals::TACRange[1], monitorsglobals::TACRange[2], (int) (monitorsglobals::exRange[2]-monitorsglobals::exRange[1])/monitorsglobals::exRange[0]*1000, monitorsglobals::exRange[1], monitorsglobals::exRange[2]);
   monitorsglobals::htac2Ex = new TH2F("htac2Ex", "Ex - TAC2 ; TAC [10ns]; Ex [MeV]", monitorsglobals::TAC2Range[0], monitorsglobals::TAC2Range[1], monitorsglobals::TAC2Range[2], (int) (monitorsglobals::exRange[2]-monitorsglobals::exRange[1])/monitorsglobals::exRange[0]*1000, monitorsglobals::exRange[1], monitorsglobals::exRange[2]);
   
   monitorsglobals::htacTdiff  = new TH2F("htacTdiff", "TDiff vs TAC; TAC [a.u.]; tDiff [ch=10ns]", monitorsglobals::TACRange[0], monitorsglobals::TACRange[1], monitorsglobals::TACRange[2], monitorsglobals::coinTimeRange[1] - monitorsglobals::coinTimeRange[0], monitorsglobals::coinTimeRange[0], monitorsglobals::coinTimeRange[1]);
   monitorsglobals::htacTdiffg = new TH2F("htacTdiffg", "TDiff vs TAC (recoil gate); TAC [a.u.]; tDiff [ch=10ns]", monitorsglobals::TACRange[0], monitorsglobals::TACRange[1], monitorsglobals::TACRange[2], monitorsglobals::coinTimeRange[1] - monitorsglobals::coinTimeRange[0], monitorsglobals::coinTimeRange[0], monitorsglobals::coinTimeRange[1]);

   for (Int_t i = 0; i < NRDT; i++){
      monitorsglobals::htacRecoil[i] = new TH2F(Form("htacRecoil%d", i), Form("RDT-%d - TAC; TAC ; RDT ", i), monitorsglobals::TACRange[0], monitorsglobals::TACRange[1] , monitorsglobals::TACRange[2], 200, 0, 4000);
      
      if( i % 2 == 0 ) {
         int tempID2 = i / 2;
         monitorsglobals::htacRecoilsum[tempID2] = new TH2F(Form("htacRecoilsum%d", tempID2), Form("(RDT-%d+RDT-%d) - TAC; TAC ; Eres+dE ", i, i+1), monitorsglobals::TACRange[0], monitorsglobals::TACRange[1], monitorsglobals::TACRange[2], 200, 0, 8000);
      }
   }
   
   //===================== energy spectrum
   monitorsglobals::hEx    = new TH1F("hEx",Form("excitation spectrum w/ goodFlag; Ex [MeV] ; Count / %4.0f keV", monitorsglobals::exRange[0]), (int) (monitorsglobals::exRange[2]-monitorsglobals::exRange[1])/monitorsglobals::exRange[0]*1000, monitorsglobals::exRange[1], monitorsglobals::exRange[2]);
   
   monitorsglobals::hExCut1  = new TH1F("hExCut1",Form("excitation spectrum w/ goodFlag; Ex [MeV] ; Count / %4.0f keV", monitorsglobals::exRange[0]), (int) (monitorsglobals::exRange[2]-monitorsglobals::exRange[1])/monitorsglobals::exRange[0]*1000, monitorsglobals::exRange[1], monitorsglobals::exRange[2]);
   monitorsglobals::hExCut1->SetLineColor(9);
   
   monitorsglobals::hExCut2  = new TH1F("hExCut2",Form("excitation spectrum w/ goodFlag; Ex [MeV] ; Count / %4.0f keV", monitorsglobals::exRange[0]), (int) (monitorsglobals::exRange[2]-monitorsglobals::exRange[1])/monitorsglobals::exRange[0]*1000, monitorsglobals::exRange[1], monitorsglobals::exRange[2]);
   monitorsglobals::hExCut2->SetLineColor(8);
   
   for(int i = 0 ; i < NARRAY; i++){
      monitorsglobals::hExi[i] = new TH1F(Form("hExi%02d", i), Form("Ex (det=%i) w/goodFlag; Ex [MeV]; Count / %4.0f keV",i, monitorsglobals::exRange[0]), (int) (monitorsglobals::exRange[2]-monitorsglobals::exRange[1])/monitorsglobals::exRange[0]*1000, monitorsglobals::exRange[1], monitorsglobals::exRange[2]);
   }
   
   for(int i = 0 ; i < numCol; i++){
      monitorsglobals::hExc[i] = new TH1F(Form("hExc%d", i), Form("Ex (col=%d) w/goodFlag; Ex [MeV]; Count / %4.0f keV", i, monitorsglobals::exRange[0]), (int) (monitorsglobals::exRange[2]-monitorsglobals::exRange[1])/monitorsglobals::exRange[0]*1000, monitorsglobals::exRange[1], monitorsglobals::exRange[2]);
   }

   monitorsglobals::hExThetaCM = new TH2F("hExThetaCM", "Ex vs ThetaCM; ThetaCM [deg]; Ex [MeV]", 200, monitorsglobals::thetaCMRange[0], monitorsglobals::thetaCMRange[1],  (int) (monitorsglobals::exRange[2]-monitorsglobals::exRange[1])/monitorsglobals::exRange[0]*1000, monitorsglobals::exRange[1], monitorsglobals::exRange[2]);

   //===================== ELUM
   for( int i = 0; i < NELUM; i++){
      monitorsglobals::helum[i] = new TH1F(Form("helum%d", i), Form("Elum-%d", i), 500, monitorsglobals::elumRange[0], monitorsglobals::elumRange[1]);
   }
   if( NELUM > 0 ) monitorsglobals::helumID = new TH2F("helumID", "Elum vs ID", NELUM, 0 , NELUM, 500, monitorsglobals::elumRange[0], monitorsglobals::elumRange[1]);

   monitorsglobals::helum4D = new TH1F("helum4d", "Elum rate for Z = 1; time [min]; count / min", timeRange[1]-timeRange[0], timeRange[0], timeRange[1]); // elum rate for (d,d)
   monitorsglobals::helum4C = new TH1F("helum4C", "Elum rate for carbon; time [min]; count / min", timeRange[1]-timeRange[0], timeRange[0], timeRange[1]); // elum rate for (12C, 12C)
   monitorsglobals::hBIC = new TH1F("hBIC", "BIC rate ; time [min]; count / min", timeRange[1]-timeRange[0], timeRange[0], timeRange[1]); // elum rate for (d,d)
   
   //===================== EZERO
   monitorsglobals::hic0 = new TH1F("hic0", "IC0; IC-0 [ch]; count", 500,  0, monitorsglobals::icRange[0]);
   monitorsglobals::hic1 = new TH1F("hic1", "IC1; IC-1 [ch]; count", 500,  0, monitorsglobals::icRange[1]);
   monitorsglobals::hic2 = new TH1F("hic2", "IC2; IC-2 [ch]; count", 500,  0, monitorsglobals::icRange[2]);

   monitorsglobals::hic01 = new TH2F("hic01", "IC0 - IC1; IC-1 [ch]; IC-0[ch]", 500,  0, monitorsglobals::icRange[1], 500, 0, monitorsglobals::icRange[0]);
   monitorsglobals::hic02 = new TH2F("hic02", "IC0 vs IC0+IC1; IC-2 [ch]; IC-0[ch]", 500, 0, monitorsglobals::icRange[1]+monitorsglobals::icRange[0], 500, 0, monitorsglobals::icRange[0]);
   monitorsglobals::hic12 = new TH2F("hic12", "IC1 - IC2; IC-2 [ch]; IC-1[ch]", 500, 0, monitorsglobals::icRange[1], 500, 0, monitorsglobals::icRange[2]);
   
   printf("======================================== End of histograms Declaration\n");
   monitorsglobals::StpWatch.Start();

}

/*###########################################################
 * Process
###########################################################*/
Bool_t Monitors::Process(Long64_t entry){

   monitorsglobals::ProcessedEntries++;
   
   /*********** Progress Bar ******************************************/ 
   if (monitorsglobals::ProcessedEntries>monitorsglobals::NumEntries*monitorsglobals::Frac-1) {
      TString msg; msg.Form("%llu", monitorsglobals::NumEntries/1000);
      int len = msg.Sizeof();
      printf(" %3.0f%% (%*llu/%llu k) processed in %6.1f sec | expect %6.1f sec",
         monitorsglobals::Frac*100, len, monitorsglobals::ProcessedEntries/1000,monitorsglobals::NumEntries/1000,monitorsglobals::StpWatch.RealTime(), monitorsglobals::StpWatch.RealTime()/monitorsglobals::Frac);
         monitorsglobals::StpWatch.Start(kFALSE);
      printf("\n\033[A\r");
      monitorsglobals::Frac+=0.1;
   }

   b_runID->GetEntry(entry);
   b_Energy->GetEntry(entry);
   b_XF->GetEntry(entry);
   b_XN->GetEntry(entry);
   b_Ring->GetEntry(entry);
   b_EnergyTimestamp->GetEntry(entry);
   
   if( isRDTExist ){
      b_RDT->GetEntry(entry);
      b_RDTTimestamp->GetEntry(entry);
   }
   
   if( isCRDTExist ){
      b_CRDT->GetEntry(entry);
      b_CRDTTimestamp->GetEntry(entry);
   }
   
   if( isAPOLLOExist ){
      b_APOLLO->GetEntry(entry);
      b_APOLLOTimestamp->GetEntry(entry);
   }
   
   if( isTACExist ){
      b_TAC->GetEntry(entry);
      b_TACTimestamp->GetEntry(entry);
   }
   
   if( isELUMExist ){
      b_ELUM->GetEntry(entry);
      b_ELUMTimestamp->GetEntry(entry);
   }
   
   if( isEZEROExist ) {
      b_EZERO->GetEntry(entry);
      b_EZEROTimestamp->GetEntry(entry);
   }   
   
   if( isArrayTraceExist ) {
      ///b_Trace_Energy->GetEntry(entry);
      b_Trace_Energy_RiseTime->GetEntry(entry);
      b_Trace_Energy_Time->GetEntry(entry);
   }
   
   if( isRDTTraceExist ){
      ///b_Trace_RDT->GetEntry(entry); 
      b_Trace_RDT_Time->GetEntry(entry);
      b_Trace_RDT_RiseTime->GetEntry(entry);
   }
   
   /*********** forming canvas Title **********************************/ 
    if( entry == 0 ) {
       if( runID == monitorsglobals::lastRunID + 1 ) {
          int len = monitorsglobals::canvasTitle.Sizeof();
          if( monitorsglobals::contFlag == false) {
             monitorsglobals::canvasTitle.Remove(len-3);
             monitorsglobals::canvasTitle += " - ";
          }
          if( monitorsglobals::contFlag == true){
             monitorsglobals::canvasTitle.Remove(len-6);
          }
          monitorsglobals::contFlag = true;
       }
       if( runID > monitorsglobals::lastRunID + 1 ) monitorsglobals::contFlag = false;
       monitorsglobals::canvasTitle += Form("%03d, ", runID );
       monitorsglobals::lastRunID = runID;
    }
    
    /*********** initization ******************************************/
    for( int i = 0 ; i < NARRAY; i++){
       monitorsglobals::z[i] = TMath::QuietNaN();
       monitorsglobals::x[i] = TMath::QuietNaN();
       monitorsglobals::xcal[i] = TMath::QuietNaN();
       monitorsglobals::eCal[i] = TMath::QuietNaN();
    }
    
    double rdtot[4] = {TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN()};
    
    /*********** TAC **************************************************/ 
    monitorsglobals::htac->Fill(tac[2]);
   
    //if( TMath::IsNaN(tac[0]) ) return kTRUE;
    //if( !(tacGate[0] < tac[0] &&  tac[0] < tacGate[1]) ) {isTACGate=true; return kTRUE;}
    
    /*********** ELUM *************************************************/
    for( int i = 0; i < NELUM; i++){
       monitorsglobals::helum[i]->Fill(elum[i]);
       monitorsglobals::helumID->Fill(i, elum[i]);
    }
      
    if( 800 < elum[0]  && elum[0] < 1200 ) monitorsglobals::helum4D->Fill(elum_t[0]/1e8/60.); 
    
    if( NELUM > 1 && !TMath::IsNaN(elum[1]) ) monitorsglobals::hBIC->Fill(elum_t[1]/1e8/60.);
    
    int tac2 = tac_t[1]-elum_t[0];        
    monitorsglobals::htac2->Fill(tac2);
    
    ///if( 4000 < elum[0]  && elum[0] < 6000 ) helum4C->Fill(elum_t[0]/1e8/60.); 
    
    /*********** Apply Recoil correction here *************************/
    
    for( int i = 0 ; i < NRDT; i++){
       rdt[i] = rdt[i]*monitorsglobals::rdtCorr[i][0] + monitorsglobals::rdtCorr[i][1];
    }
    
    /*********** Array ************************************************/ 
    //Do calculations and fill histograms
    Int_t recoilMulti = 0;
    Int_t arrayMulti = 0;
    Int_t multiEZ = 0;
    bool rdtgate1 = false;
    bool rdtgate2 = false;
    bool coinFlag = false;
    bool ezGate = false;
    bool isGoodEventFlag = false;

    for (Int_t detID = 0; detID < NARRAY; detID++) {
      
      //================== Filling raw data
      monitorsglobals::he[detID]->Fill(e[detID]);
      monitorsglobals::hring[detID]->Fill(ring[detID]);
      monitorsglobals::hxf[detID]->Fill(xf[detID]);
      monitorsglobals::hxn[detID]->Fill(xn[detID]);
      monitorsglobals::heVring[detID]->Fill(ring[detID],e[detID]);
      monitorsglobals::hxfVxn[detID]->Fill(xf[detID],xn[detID]);
      monitorsglobals::heVxs[detID]->Fill(xf[detID]+xn[detID], e[detID]);

      monitorsglobals::heVID->Fill(detID, e[detID]);
      monitorsglobals::hringVID->Fill(detID, ring[detID]);
      monitorsglobals::hxfVID->Fill(detID, xf[detID]);
      monitorsglobals::hxnVID->Fill(detID, xn[detID]);
      
            
      //==================== Basic gate
      if( TMath::IsNaN(e[detID]) ) continue ; 
      ///if( ring[detID] < -100 || ring[detID] > 100 ) continue; 
      ///if( ring[detID] > 300 ) continue; 
      if( TMath::IsNaN(xn[detID]) &&  TMath::IsNaN(xf[detID]) ) continue ; 
      
      //==================== Skip detector 
      bool skipFlag = false;
      for( unsigned int kk = 0; kk < monitorsglobals::skipDetID.size() ; kk++){
         if( detID == monitorsglobals::skipDetID[kk] ) {
            skipFlag = true;
            break;
         }
      }
      if (skipFlag ) continue;

      //==================== Calibrations go here
      monitorsglobals::xfcal[detID] = xf[detID] * monitorsglobals::xfxneCorr[detID][1] + monitorsglobals::xfxneCorr[detID][0];
      monitorsglobals::xncal[detID] = xn[detID] * monitorsglobals::xnCorr[detID] * monitorsglobals::xfxneCorr[detID][1] + monitorsglobals::xfxneCorr[detID][0];
      monitorsglobals::eCal[detID] = e[detID] / monitorsglobals::eCorr[detID][0] + monitorsglobals::eCorr[detID][1];

      if( monitorsglobals::eCal[detID] < monitorsglobals::eCalCut[0] ) continue;
      if( monitorsglobals::eCal[detID] > monitorsglobals::eCalCut[1] ) continue;
      
      //===================== fill Calibrated  data
      monitorsglobals::heCal[detID]->Fill(monitorsglobals::eCal[detID]);
      monitorsglobals::heCalID->Fill(detID, monitorsglobals::eCal[detID]);
      monitorsglobals::hxfCalVxnCal[detID]->Fill(monitorsglobals::xfcal[detID], monitorsglobals::xncal[detID]);
      monitorsglobals::heVxsCal[detID]->Fill(monitorsglobals::xncal[detID] + monitorsglobals::xfcal[detID], e[detID]);
      
      //===================== calculate X
      if( (xf[detID] > 0 || !TMath::IsNaN(xf[detID])) && ( xn[detID]>0 || !TMath::IsNaN(xn[detID])) ) {
        ///x[detID] = 0.5*((xf[detID]-xn[detID]) / (xf[detID]+xn[detID]))+0.5;
        monitorsglobals::x[detID] = 0.5*((xf[detID]-xn[detID]) / e[detID])+0.5;
      }
      
      /// range of x is (0, 1)
      if  ( !TMath::IsNaN(xf[detID]) && !TMath::IsNaN(xn[detID]) ) monitorsglobals::xcal[detID] = 0.5 + 0.5 * (monitorsglobals::xfcal[detID] - monitorsglobals::xncal[detID] ) / e[detID];
      if  ( !TMath::IsNaN(xf[detID]) &&  TMath::IsNaN(xn[detID]) ) monitorsglobals::xcal[detID] = monitorsglobals::xfcal[detID]/ e[detID];
      if  (  TMath::IsNaN(xf[detID]) && !TMath::IsNaN(xn[detID]) ) monitorsglobals::xcal[detID] = 1.0 - monitorsglobals::xncal[detID]/ e[detID];
      
      //======= Scale xcal from (0,1)      
      monitorsglobals::xcal[detID] = (monitorsglobals::xcal[detID]-0.5)/monitorsglobals::xScale[detID] + 0.5; /// if include this scale, need to also inclused in Cali_littleTree
      
      if( abs(monitorsglobals::xcal[detID] - 0.5) > monitorsglobals::xGate/2. ) continue; 
      
      //==================== calculate Z
      if( detGeo.firstPos > 0 ) {
        z[detID] = detGeo.detLength*(1.0-monitorsglobals::xcal[detID]) + detGeo.detPos[detID%numCol];
      }else{
        z[detID] = detGeo.detLength*(monitorsglobals::xcal[detID]-1.0) + detGeo.detPos[detID%numCol];
      }

      //===================== multiplicity
      arrayMulti++; /// multi-hit when both e, xf, xn are not NaN

      //=================== Array fill
      monitorsglobals::heVx[detID]->Fill(monitorsglobals::x[detID],e[detID]);
      monitorsglobals::hringVx[detID]->Fill(monitorsglobals::x[detID],ring[detID]);
     
      monitorsglobals::heCalVxCal[detID]->Fill(monitorsglobals::xcal[detID]*detGeo.detLength,monitorsglobals::eCal[detID]);
      monitorsglobals::heCalVz->Fill(z[detID],monitorsglobals::eCal[detID]);

      //=================== Recoil Gate
      if( isRDTExist && (monitorsglobals::isCutFileOpen1 || monitorsglobals::isCutFileOpen2)){
        for(int i = 0 ; i < monitorsglobals::numCut1 ; i++ ){
            monitorsglobals::cutG = (TCutG *)monitorsglobals::cutList1->At(i) ;
            if(monitorsglobals::cutG->IsInside(rdt[2*i],rdt[2*i+1])) {
            // if(cutG->IsInside(rdt[2*i] + rdt[2*i+1],rdt[2*i+1])) {

            rdtgate1= true;
            break; /// only one is enough
          }
        }
        
        for(int i = 0 ; i < monitorsglobals::numCut2 ; i++ ){
          monitorsglobals::cutG = (TCutG *)monitorsglobals::cutList2->At(i) ;
          if(monitorsglobals::cutG->IsInside(rdt[2*i],rdt[2*i+1])) {
          //if(cutG->IsInside(rdt[2*i]+ rdt[2*i+1],rdt[2*i+1])) {

            rdtgate2= true;
            break; /// only one is enough
          }
        }
        
      }else{
        rdtgate1 = true;
        rdtgate2 = true;
      }
      
      //================ coincident with Recoil when z is calculated.
      if( !TMath::IsNaN(z[detID]) ) { 
         for( int j = 0; j < NRDT ; j++){
            if( TMath::IsNaN(rdt[j]) ) continue; 
   
            int tdiff = rdt_t[j] - e_t[detID];
   
            if( j%2 == 0) { // it is plotting dE now, was 1 and E
               monitorsglobals::htdiff->Fill(tdiff);
               monitorsglobals::htacTdiff->Fill( tac[0], tdiff);
               if((rdtgate1 || rdtgate2) && (monitorsglobals::eCalCut[1] > monitorsglobals::eCal[detID] && monitorsglobals::eCal[detID]>monitorsglobals::eCalCut[0])) {
                  monitorsglobals::htdiffg->Fill(tdiff);
                  monitorsglobals::htacTdiffg->Fill( tac[0], tdiff);
               }
            }

            monitorsglobals::hArrayRDTMatrix->Fill(detID, j); 
      
            if( monitorsglobals::isTimeGateOn && monitorsglobals::timeGate[0] < tdiff && tdiff < monitorsglobals::timeGate[1] ) {
               if (monitorsglobals::isTACGate && !(monitorsglobals::tacGate[0] < tac[0] &&  tac[0] < monitorsglobals::tacGate[1])) continue;
               if(j % 2 == 0 ) monitorsglobals::hrdt2Dg[j/2]->Fill(rdt[j],rdt[j+1]); /// x=E, y=dE
               ///if(j % 2 == 0 ) hrdt2Dg[j/2]->Fill(rdt[j+1],rdt[j]); /// x=dE, y=E
               monitorsglobals::hArrayRDTMatrixG->Fill(detID, j); 
               ///if( rdtgate1) hArrayRDTMatrixG->Fill(detID, j); 
               
               monitorsglobals::hrdtg[j]->Fill(rdt[j]);
               coinFlag = true;
               
            }
         }
      }
      
      if( !monitorsglobals::isTimeGateOn ) coinFlag = true;
      
      //================ E-Z gate
      if( monitorsglobals::isEZCutFileOpen ) {
         
         if( monitorsglobals::EZCut->IsInside(z[detID], monitorsglobals::eCal[detID])  ){
            ezGate = true;
         }
         
      }else{
         ezGate = true;
      }
      
      if( coinFlag && (rdtgate1 || rdtgate2) && ezGate){ 
         monitorsglobals::heCalVzGC->Fill( z[detID] , monitorsglobals::eCal[detID] );
        
         monitorsglobals::heCalVxCalG[detID]->Fill(monitorsglobals::xcal[detID]*detGeo.detLength,monitorsglobals::eCal[detID]);
         monitorsglobals::heVIDG->Fill(detID, e[detID]);

         for( int i = 0; i < numRow; i++){
            for(int j = 0; j < numCol; j++){
               int k = numCol*i+j;
               monitorsglobals::hecalVzRowG[i] -> Fill( z[k], monitorsglobals::eCal[k]);
            }
         }
      
         multiEZ ++;
         isGoodEventFlag = true;

      }
      
   }//end of array loop
   
   if( !monitorsglobals::isEZCutFileOpen ) ezGate = true;
   
   //=========== fill eCal Vs z for each row
   for( int i = 0; i < numRow; i++){
      for(int j = 0; j < numCol; j++){
         int k = numCol*i+j;
         monitorsglobals::hecalVzRow[i] -> Fill( z[k], monitorsglobals::eCal[k]);
      }
   }
   
   /*********** RECOILS ***********************************************/    
   for( int i = 0; i < NRDT ; i++){
      monitorsglobals::hrdtID->Fill(i, rdt[i]);
      monitorsglobals::hrdt[i]->Fill(rdt[i]);
      
      for( int j = 0; j < NRDT ; j++){
         if( rdt[i] > 0 && rdt[j] > 0 )  monitorsglobals::hrdtMatrix->Fill(i, j);
      }
      
      if( i % 2 == 0  ){        
        if ( monitorsglobals::isTACGate && !(monitorsglobals::tacGate[0] < tac[0] &&  tac[0] < monitorsglobals::tacGate[1]) ) continue;        
         recoilMulti++; // when both dE and E are hit
         rdtot[i/2] = rdt[i]+rdt[i+1];
         monitorsglobals::htacRecoilsum[i/2]->Fill(tac[0],rdtot[i/2]);
         monitorsglobals::hrdt2D[i/2]->Fill(rdt[i],rdt[i+1]); //E-dE
         monitorsglobals::hrdt2Dsum[i/2]->Fill(rdtot[i/2],rdt[i+1]);//dE-(dE+E)

         monitorsglobals::htacRecoil[i]->Fill(tac[0],rdt[i]);
         monitorsglobals::htacRecoil[i+1]->Fill(tac[0],rdt[i+1]);
      }
   }
   
   /*********** Apollo ***********************************************/    
   
   for( int i = 0; i < NAPOLLO ; i++){
      monitorsglobals::hApollo[i]->Fill(apollo[i]);
   }
   
   ///if( rdt_t[4] > 0 ){
   ///   if( abs(rdt[4] - 1658) < 40) hrdtRate1->Fill(rdt_t[4]/1e8/60.);
   ///   if( abs(rdt[4] - 1783) < 40) hrdtRate2->Fill(rdt_t[4]/1e8/60.);
   ///}

   /******************* Circular Recoil *******************************/
   ///======= 0 -  7 is angular 
   ///======= 8 - 15 is radial
   
   for( int i = 0; i < NCRDT/2 ; i++){
     if( TMath::IsNaN(crdt[i]) ) continue;
     monitorsglobals::hcrdt[i]->Fill(crdt[i]);
     
     for( int j = NCRDT/2; j < NCRDT; j++){
      monitorsglobals::hcrdtID->Fill(i, j);

      double theta = -TMath::Pi() + 2*TMath::Pi()/8.*(i+0.5);
      double rho   = 10.+40./8.*(j+0.5);

      monitorsglobals::hcrdtPolar->Fill( theta, rho );
    }
   }
   
   /******************* Multi-hit *************************************/
   monitorsglobals::hmultEZ->Fill(multiEZ);
   monitorsglobals::hmult->Fill(recoilMulti,arrayMulti);


   /*********** EZERO *************************************************/ 
   //if( ezGate ) {
   // hic0->Fill(ezero[0]);
   // hic1->Fill(ezero[1]);
   // hic2->Fill(ezero[2]);
   //
   // hic01->Fill(ezero[1], ezero[0]);
   // hic02->Fill(ezero[1]+ezero[0], ezero[0]);
   // hic12->Fill(ezero[2], ezero[1]);
   //   
   //}
   
   /*********** Good event Gate ***************************************/ 
   if( !isGoodEventFlag ) return kTRUE;

   /*********** Ex and thetaCM ****************************************/ 
   for(Int_t detID = 0; detID < NARRAY ; detID++){
     	
     if( TMath::IsNaN(e[detID]) ) continue ; 
     if( TMath::IsNaN(z[detID]) ) continue ;
     if( monitorsglobals::eCal[detID] <  monitorsglobals::eCalCut[0] ) continue ;
     if( monitorsglobals::eCal[detID] >  monitorsglobals::eCalCut[1] ) continue ;

     if( monitorsglobals::isReaction ){
       ///======== Ex calculation by Ryan 
       double y = monitorsglobals::eCal[detID] + monitorsglobals::mass; // to give the KE + mass of proton;
       double Z = monitorsglobals::alpha * monitorsglobals::gamm * monitorsglobals::betRel * z[detID];
       double H = TMath::Sqrt(TMath::Power(monitorsglobals::gamm * monitorsglobals::betRel,2) * (y*y - monitorsglobals::mass * monitorsglobals::mass) ) ;
 
       if( TMath::Abs(Z) < H ) {
         ///using Newton's method to solve 0 ==	H * sin(phi) - G * tan(phi) - Z = f(phi) 
         double tolerrence = 0.001;
         double phi = 0;  ///initial phi = 0 -> ensure the solution has f'(phi) > 0
         double nPhi = 0; /// new phi

         int iter = 0;
         do{
           phi = nPhi;
           nPhi = phi - (H * TMath::Sin(phi) - monitorsglobals::G * TMath::Tan(phi) - Z) / (H * TMath::Cos(phi) - monitorsglobals::G /TMath::Power( TMath::Cos(phi), 2));					 
           iter ++;
           if( iter > 10 || TMath::Abs(nPhi) > TMath::PiOver2()) break;
         }while( TMath::Abs(phi - nPhi ) > tolerrence);
         phi = nPhi;

         /// check f'(phi) > 0
         double Df = H * TMath::Cos(phi) - monitorsglobals::G / TMath::Power( TMath::Cos(phi),2);
         if( Df > 0 && TMath::Abs(phi) < TMath::PiOver2()  ){
           double K = H * TMath::Sin(phi);
           double x = TMath::ACos( monitorsglobals::mass / ( y * monitorsglobals::gamm - K));
           double momt = monitorsglobals::mass * TMath::Tan(x); /// momentum of particel b or B in CM frame
           double EB = TMath::Sqrt(monitorsglobals::mass*monitorsglobals::mass + monitorsglobals::Et*monitorsglobals::Et - 2*monitorsglobals::Et*TMath::Sqrt(momt*momt + monitorsglobals::mass * monitorsglobals::mass));
           monitorsglobals::Ex = EB - monitorsglobals::massB;

           double hahaha1 = monitorsglobals::gamm* TMath::Sqrt(monitorsglobals::mass * monitorsglobals::mass + momt * momt) - y;
           double hahaha2 = monitorsglobals::gamm* monitorsglobals::betRel * momt;
           monitorsglobals::thetaCM = TMath::ACos(hahaha1/hahaha2) * TMath::RadToDeg();

         }else{
           monitorsglobals::Ex = TMath::QuietNaN();
           monitorsglobals::thetaCM = TMath::QuietNaN();
         }
       }else{
         monitorsglobals::Ex = TMath::QuietNaN();
         monitorsglobals::thetaCM = TMath::QuietNaN();
       }
     }else{
       monitorsglobals::Ex = TMath::QuietNaN();
       monitorsglobals::thetaCM = TMath::QuietNaN();
     }
     
     monitorsglobals::htacEx->Fill(tac[2], monitorsglobals::Ex);
     monitorsglobals::htac2Ex->Fill(tac_t[1]-e_t[detID], monitorsglobals::Ex);
     
     if( monitorsglobals::thetaCM > monitorsglobals::thetaCMGate ) {

         monitorsglobals::hEx->Fill(monitorsglobals::Ex);

         monitorsglobals::hExThetaCM->Fill(monitorsglobals::thetaCM, monitorsglobals::Ex);
         
         if( rdtgate1 ) {
            monitorsglobals::hExCut1->Fill(monitorsglobals::Ex);
            monitorsglobals::hExThetaCM->Fill(monitorsglobals::thetaCM, monitorsglobals::Ex);
         }
         if( rdtgate2 ) {
            monitorsglobals::hExCut2->Fill(monitorsglobals::Ex);
            monitorsglobals::hExThetaCM->Fill(monitorsglobals::thetaCM, monitorsglobals::Ex);
         }
         
         monitorsglobals::hExi[detID]->Fill(monitorsglobals::Ex);
         monitorsglobals::hExVxCal[detID]->Fill(monitorsglobals::xcal[detID], monitorsglobals::Ex);
         monitorsglobals::hExc[detID%numCol]->Fill(monitorsglobals::Ex);
         
      }
   }
  
   return kTRUE;
}


void Monitors::Terminate()
{
   printf("\n");
   printf("============================== finishing.\n");

   gROOT->cd();
   
   int strLen = monitorsglobals::canvasTitle.Sizeof();
   monitorsglobals::canvasTitle.Remove(strLen-3);
   
   TString runTimeStr = "";
   if( monitorsglobals::runTime > 0. ) {
      runTimeStr = Form("%.0f min", monitorsglobals::runTime);
      monitorsglobals::canvasTitle += " | " + runTimeStr;
   }

   //############################################ User is free to edit this section
   //--- Canvas Size
   int canvasDiv[2] = {3,3};
   int padSize = 550;
   int canvasXY[2] = { canvasDiv[0] * padSize, canvasDiv[1] * padSize };
   monitorsglobals::cCanvas  = new TCanvas("monitorsglobals::cCanvas",monitorsglobals::canvasTitle + " | " + monitorsglobals::rdtCutFile1 + "," + monitorsglobals::rdtCutFile2,canvasXY[0],canvasXY[1]);
   monitorsglobals::cCanvas->Modified(); monitorsglobals::cCanvas->Update();
   monitorsglobals::cCanvas->cd(); monitorsglobals::cCanvas->Divide(canvasDiv[0],canvasDiv[1]);

   gStyle->SetOptStat("neiou");
   
   monitorsglobals::text.SetNDC();
   monitorsglobals::text.SetTextFont(82);
   monitorsglobals::text.SetTextSize(0.04);
   monitorsglobals::text.SetTextColor(2);

   double yMax = 0;

   Isotope hRecoil(reactionConfig.recoilHeavyA, reactionConfig.recoilHeavyZ);
   double Sn = hRecoil.CalSp(0,1);
   double Sp = hRecoil.CalSp(1,0);
   double Sa = hRecoil.CalSp2(4,2);
   
   //TODO, Module each plot
   ///----------------------------------- Canvas - 1
   PlotEZ(1); /// raw EZ

   ///----------------------------------- Canvas - 2
   PlotEZ(0); ///gated EZ

   ///----------------------------------- Canvas - 3
   PlotTDiff(1, 1); ///with Gated Tdiff, isLog

   ///----------------------------------- Canvas - 4
   monitorsglobals::padID++; monitorsglobals::cCanvas->cd(monitorsglobals::padID); 
   
   if( monitorsglobals::hExCut1->GetMaximum() > monitorsglobals::hExCut2->GetMaximum() ){
     monitorsglobals::hExCut1->Draw("");
     monitorsglobals::hExCut2->Draw("same");
   }else{
     monitorsglobals::hExCut2->Draw("");
     monitorsglobals::hExCut1->Draw("same");
   }
   DrawLine(monitorsglobals::hEx, Sn);
   DrawLine(monitorsglobals::hEx, Sa);
   DrawLine(monitorsglobals::hEx, 4.03, 1);
   
   if(monitorsglobals::isTimeGateOn)monitorsglobals::text.DrawLatex(0.15, 0.8, Form("%d < coinTime < %d", monitorsglobals::timeGate[0], monitorsglobals::timeGate[1])); 
   if( monitorsglobals::xGate < 1 ) monitorsglobals::text.DrawLatex(0.15, 0.75, Form("with |x-0.5|<%.4f", monitorsglobals::xGate/2.));
   if( monitorsglobals::isCutFileOpen1 ) monitorsglobals::text.DrawLatex(0.15, 0.7, "with recoil gated"); 

   ///----------------------------------- Canvas - 5
   PlotRDT(0, 1);

   ///----------------------------------- Canvas - 6
   PlotRDT(1, 1);

   ///----------------------------------- Canvas - 7
   monitorsglobals::padID++; monitorsglobals::cCanvas->cd(monitorsglobals::padID);

   monitorsglobals::hEx->Draw();
   DrawLine(monitorsglobals::hEx, Sn);
   DrawLine(monitorsglobals::hEx, Sa);
   DrawLine(monitorsglobals::hEx, 4.03, 8);

   if(monitorsglobals::isTimeGateOn)monitorsglobals::text.DrawLatex(0.15, 0.8, Form("%d < coinTime < %d", monitorsglobals::timeGate[0], monitorsglobals::timeGate[1])); 
   if( monitorsglobals::xGate < 1 ) monitorsglobals::text.DrawLatex(0.15, 0.75, Form("with |x-0.5|<%.4f", monitorsglobals::xGate/2.));

   ///----------------------------------- Canvas - 8
   PlotRDT(3, 1);
   
   ///----------------------------------- Canvas - 9
   PlotRDT(2, 1);
   
   /************************************/
   gStyle->GetAttDate()->SetTextSize(0.02);
   gStyle->SetOptDate(1);
   gStyle->SetDateX(0);
   gStyle->SetDateY(0);

   /************************************/
   monitorsglobals::StpWatch.Start(kFALSE);
   
   gROOT->ProcessLine(".L ../Armory/Monitors_Util.C+");
   //gROOT->ProcessLine(Form("FindBesCanvasDivision(%d)", numDet));
   printf("=============== loaded Monitors_Utils.C\n");
   gROOT->ProcessLine(".L ../Armory/AutoFit.C");
   printf("=============== loaded Armory/AutoFit.C\n");
   gROOT->ProcessLine(".L ../Armory/RDTCutCreator.C");
   printf("=============== loaded Armory/RDTCutCreator.C\n");
   gROOT->ProcessLine(".L ../Armory/Check_rdtGate.C");
   printf("=============== loaded Armory/Check_rdtGate.C\n");
   // gROOT->ProcessLine(".L ../Armory/readTrace.C");
   // printf("=============== loaded Armory/readTrace.C\n");
   // gROOT->ProcessLine(".L ../Armory/readRawTrace.C");
   // printf("=============== loaded Armory/readRawTrace.C\n");
   gROOT->ProcessLine("listDraws()");
   
   /************************* Save histograms to root file*/
   
   gROOT->cd();
   //TString outFileNameTemp = monitorsglobals::canvasTitle;
   //outFileNameTemp.ReplaceAll(" - ", "-").ReplaceAll(" ", "_").ReplaceAll(":", "").ReplaceAll(",","");
   //gROOT->GetList()->SaveAs(outFileName + ".root");
   
   if( printControlID == 0 ) {
      TDatime dateTime;
      TString outFileName;
      outFileName.Form("Canvas_%d%02d%02d_%06d.png", dateTime.GetYear(), dateTime.GetMonth(), dateTime.GetDay(), dateTime.GetTime());
      
      monitorsglobals::cCanvas->SaveAs("Canvas/"+outFileName);
      
      /// run puhs_to_websrv.sh
      TString cmd;
      printf("|%s|\n", monitorsglobals::canvasTitle.Data());
      
      cmd.Form(".! ../Armory/push_to_websrv.sh %s %s", outFileName.Data(), monitorsglobals::canvasTitle.ReplaceAll(" ", "").ReplaceAll("|","_").Data()); 
      gROOT->ProcessLine(cmd);
       
      ///======================== exit after 
      gROOT->ProcessLine(".q");
   }
   
   
   /************************************/
   //gROOT->ProcessLine("recoils()");


}

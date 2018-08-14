// Standard C++ headers
#include "od_analysis.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// ROOT Headers
#include "TROOT.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TMath.h"
#include "TTimeStamp.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TBranch.h"
#include "TF1.h"

using namespace std;
void Loggify(TAxis* axis) {
  int bins = axis->GetNbins();

  Axis_t from = TMath::Log10(axis->GetXmin());
  Axis_t to = TMath::Log10(axis->GetXmax());
  Axis_t width = (to - from) / bins;
  Axis_t *new_bins = new Axis_t[bins + 1];
  
  for (int i = 0; i <= bins; i++) new_bins[i] = TMath::Power(10, from + i * width);
  axis->Set(bins, new_bins); 
  delete[] new_bins; 
}

void LoggifyX(TH1* h)  { Loggify(h->GetXaxis()); }
void LoggifyXY(TH2* h) { Loggify(h->GetXaxis()); Loggify(h->GetYaxis()); }
void LoggifyX(TH2* h)  { Loggify(h->GetXaxis()); }
void LoggifyY(TH2* h)  { Loggify(h->GetYaxis()); }




void od_analysis(){
  gROOT->SetBatch(kTRUE);
  TStopwatch* clock = new TStopwatch();
  clock->Start();


  TString sample="background";
  //TString sample="NA22";
  //TString sample="Kr83";
  //  TString sample="AmLi";
  TString outname  = sample+".root";
  TString txtFileList = "RQfile.list."+sample;
  
  outfile = new TFile(outname, "recreate");
  cout << "Writing output to: "<<outfile->GetName()<<endl;

  TChain* chain = new TChain("Events", "Events");
  load_chain(txtFileList, chain);
  MyEvent* evt = new MyEvent();
  evt->LoadBranches(chain);


  //define histos
  TH1D* h_singlescatter_s1Area_phd = new TH1D("singlescatter_s1Area_phd", "singlescatter_s1Area_phd", 10000, 0, 10000);
  TH1D* h_multiplescatter_s1Area_phd = new TH1D("multiplescatter_s1Area_phd", "multiplescatter_s1Area_phd", 10000, 0, 10000);
//  TH1D* h_doublescatter_s1Area_phd = new TH1D("doublescatter_s1Area_phd", "doublescatter_s1Area_phd", 10000, 0, 10000);
//  TH1D* h_triplescatter_s1Area_phd = new TH1D("triplescatter_s1Area_phd", "triplescatter_s1Area_phd", 10000, 0, 10000);
//  TH1D* h_quadscatter_s1Area_phd = new TH1D("quadscatter_s1Area_phd", "quadscatter_s1Area_phd", 10000, 0, 10000);
//  TH1D* h_quintscatter_s1Area_phd = new TH1D("quintscatter_s1Area_phd", "quintscatter_s1Area_phd", 10000, 0, 10000);
 
   TH1D* h_s1pulsewidth_TPCHG = new TH1D("s1pulsewidth_TPCHG","s1pulsewidth_TPCHG",10000,0,10000);
   TH1D* h_s2pulsewidth_TPCHG = new TH1D("s2pulsewidth_TPCHG","s2pulsewidth_TPCHG",10000,0,10000);

   TH1D* h_pulsearea_nos2_TPCHG = new TH1D("pulsearea_nos2_TPCHG","pulsearea_nos2_TPCHG",10000,0,10000);
   TH1D* h_pulsearea_yess2_TPCHG = new TH1D("pulsearea_yess2_TPCHG","pulsearea_yess2_TPCHG",10000,0,10000);

  //------------------------------------------------
  // Main event loop
  //------------------------------------------------
  const Int_t nevents = evt->chain->GetEntries();
  
  evt->GetEntry(0);

  int processed_events=0;
  //  for (Int_t n=0; n<nevents; ++n) {
  for (Int_t n=0; n<1000; ++n) {
    if (n%1000 == 0) cout << "Processing "<< n << "/"<<nevents<<endl;
    processed_events++;
    evt->GetEntry(n);

   for (int lop = 0; lop < evt->singlescatter_s1Area_phd.size(); lop++){
    if(evt->nSingleScatters>0){
      if (evt->nSingleScatters>1) cout<<"WARNING, nSingleScatters>1"<<endl;
      h_singlescatter_s1Area_phd->Fill(evt->singlescatter_s1Area_phd[lop]);
    }
   }
   for (int lop = 0; lop < evt->multiplescatter_s1Area_phd.size(); lop++){
   if(evt->nMultipleScatters>0){
      if (evt->nMultipleScatters>1) cout<<"It's more than 1 in the nMultipleScatters"<<endl;
      h_multiplescatter_s1Area_phd->Fill(evt->multiplescatter_s1Area_phd[lop]);
    }//int nevents
    }

  if ((evt->pulseStartTime_ns_TPCHG.size()) >= 1) {
//    int start = evt->pulseStartTime_ns_TPCHG[0];
//    int end = evt->pulseEndTime_ns_TPCHG[0];
//    int diff = end - start;
//   cout<< evt->singlescatter_s1Area_phd.size()<<endl;
//   cout<< evt->multiplescatter_s1Area_phd.size()<<endl;
//   cout<< evt->pulseStartTime_ns_TPCHG.size()<<endl;
//    cout << start << endl;
//    cout << end << endl;
//    cout<< evt->pulseArea_phd_TPCHG[0]<<endl;

//    cout << "" << endl;
//    cout << evt->s1Probability_TPCHG[0] << endl;
//    cout << evt->s2Probability_TPCHG[0] << endl;
   for (int lop = 0; lop < (evt->pulseStartTime_ns_TPCHG.size()-1); lop++) {
    int start = evt->pulseStartTime_ns_TPCHG[lop];
    int end = evt->pulseEndTime_ns_TPCHG[lop];
    int diff = end - start;


    if (evt->s1Probability_TPCHG[0] > 0) {
	h_s1pulsewidth_TPCHG->Fill(diff);
	}
    else {
	if (evt->s2Probability_TPCHG[0] > 0) {
		h_s2pulsewidth_TPCHG->Fill(diff);
		h_pulsearea_yess2_TPCHG->Fill(evt->pulseArea_phd_TPCHG[lop]);
		}
	else {
		h_pulsearea_nos2_TPCHG->Fill(evt->pulseArea_phd_TPCHG[lop]);
	}
}
   }


    }


  }  

  //write and close output file
  outfile->Write();
  outfile->Close();

  delete chain;
  delete clock;
  
  }

void load_chain(TString txtFileList, TChain* chain){
  cout << "Loading file names from "<<txtFileList << " into "<<chain->GetName()<<endl;
  ifstream fileList(txtFileList);
  string file;
  if (fileList.is_open()) {
    while ( getline(fileList, file) ) {
      chain->AddFile(file.c_str());
    }
    fileList.close();
  }
  else{
    cout<<"The file "<< txtFileList <<" doesn't exist. Exiting !!"<<endl;
    exit(-1);
  }
}


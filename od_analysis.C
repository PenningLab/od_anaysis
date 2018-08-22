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

TH1D* pulseArea[8];		// 0:total, 1:s1 only, 2:s2 only, 3:other only, 4:max S1, 5:max S2, 6:sub S1, 7:sub s2
TH1D* areaFrac5[8];
TH1D* areaFrac95[8];
TH1D* index[8];
TH1D* s1prob[8];
TH1D* s2prob[8];
TH1D* peakAmp[8];
TH2D* area_vs_length[8];

for  (int i=0; i<8; i++) {
char hname[50];
sprintf(hname,"pulseArea_%u",i);
pulseArea[i] = new TH1D(hname,"",10000,0,10000);
sprintf(hname,"areaFrac5_%u",i);
areaFrac5[i] = new TH1D(hname,"",300,0,3000);
sprintf(hname,"areaFrac95_%u",i);
areaFrac95[i] = new TH1D(hname,"",300,0,3000);
sprintf(hname,"index_%u",i);
index[i] = new TH1D(hname,"",500,0,500);
sprintf(hname,"s1prob_%u",i);
s1prob[i] = new TH1D(hname,"",2,0,2);
sprintf(hname,"s2prob_%u",i);
s2prob[i] = new TH1D(hname,"",2,0,2);
sprintf(hname,"peakAmp_%u",i);
peakAmp[i] = new TH1D(hname,"",600,0,3);
sprintf(hname,"area_vs_length_%u",i);
area_vs_length[i] = new TH2D(hname,"",10000,0,10000,300,0,3000);
}

  //define histos
  TH1D* h_singlescatter_s1Area_phd = new TH1D("singlescatter_s1Area_phd", "singlescatter_s1Area_phd", 10000, 0, 10000);
  TH1D* h_multiplescatter_s1Area_phd = new TH1D("multiplescatter_s1Area_phd", "multiplescatter_s1Area_phd", 10000, 0, 10000);
//  TH1D* h_doublescatter_s1Area_phd = new TH1D("doublescatter_s1Area_phd", "doublescatter_s1Area_phd", 10000, 0, 10000);
//  TH1D* h_triplescatter_s1Area_phd = new TH1D("triplescatter_s1Area_phd", "triplescatter_s1Area_phd", 10000, 0, 10000);
//  TH1D* h_quadscatter_s1Area_phd = new TH1D("quadscatter_s1Area_phd", "quadscatter_s1Area_phd", 10000, 0, 10000);
//  TH1D* h_quintscatter_s1Area_phd = new TH1D("quintscatter_s1Area_phd", "quintscatter_s1Area_phd", 10000, 0, 10000);
 
   TH1D* h_s1pulsewidth_TPCHG = new TH1D("s1pulsewidth_TPCHG","s1pulsewidth_TPCHG",10000,0,10000);
   TH1D* h_s2pulsewidth_TPCHG = new TH1D("s2pulsewidth_TPCHG","s2pulsewidth_TPCHG",10000,0,10000);

//   TH1D* h_pulsearea_nos2_TPCHG = new TH1D("pulsearea_nos2_TPCHG","pulsearea_nos2_TPCHG",10000,0,10000);
//   TH1D* h_pulsearea_yess2_TPCHG = new TH1D("pulsearea_yess2_TPCHG","pulsearea_yess2_TPCHG",10000,0,10000);

//TH1D* h_s1prob_total = new TH1D("s1prob_total","s1prob_total",2,0,2);
//TH1D* h_s2prob_total = new TH1D("s2prob_total","s2prob_total",2,0,2);
//TH1D* h_pulsearea_total = new TH1D("pulsearea_total","pulsearea_total",2000,0,2000);
//TH1D* h_singlePEprob_total = new TH1D("singlePEprob_total","singlePEprob_total",2,0,2);
//TH1D* h_otherprob_total = new TH1D("otherprob_total","otherprob_total",2,0,2);
//TH1D* h_others2prob_total = new TH1D("others2prob_total","others2prob_total",2,0,2);
//TH1D* h_areafrac5_total = new TH1D("areafrac5_total","areafrac5_total",300,0,3000);
//TH1D* h_areafrac95_total = new TH1D("areafrac95_total","areafrac95_total",300,0,3000);
//TH1D* h_peakamp_total = new TH1D("peakamp_total","peakamp_total",1000,0,50);

//TH1D* h_s1prob_s1 = new TH1D("s1prob_s1","s1prob_s1",2,0,2);
//TH1D* h_s2prob_s1 = new TH1D("s2prob_s1","s2prob_s1",2,0,2);
//TH1D* h_pulsearea_s1 = new TH1D("pulsearea_s1","pulsearea_s1",2000,0,2000);
//TH1D* h_pulsearea_maxs1 = new TH1D("pulsearea_maxs1","pulsearea_maxs1",2000,0,2000);
//TH1D* h_pulsearea_subs1 = new TH1D("pulsearea_subs1","pulsearea_subs1",2000,0,2000);
//TH1D* h_singlePEprob_s1 = new TH1D("singlePEprob_s1","singlePEprob_s1",2,0,2);
//TH1D* h_otherprob_s1 = new TH1D("otherprob_s1","otherprob_s1",2,0,2);
//TH1D* h_others2prob_s1 = new TH1D("others2prob_s1","others2prob_s1",2,0,2);
//TH1D* h_areafrac5_s1 = new TH1D("areafrac5_s1","areafrac5_s1",300,0,3000);
//TH1D* h_areafrac95_s1 = new TH1D("areafrac95_s1","areafrac95_s1",300,0,3000);
//TH1D* h_peakamp_maxs1 = new TH1D("peakamp_maxs1","peakamp_maxs1",1000,0,50);
//TH1D* h_index_maxs1 = new TH1D("index_maxs1","index_maxs1",100,0,1000);
//TH1D* h_index_subs1 = new TH1D("index_subs1","index_subs1",100,0,1000);

//TH1D* h_s1prob_s2 = new TH1D("s1prob_s2","s1prob_s2",2,0,2);
//TH1D* h_s2prob_s2 = new TH1D("s2prob_s2","s2prob_s2",2,0,2);
//TH1D* h_pulsearea_s2 = new TH1D("pulsearea_s2","pulsearea_s2",2000,0,2000);
//TH1D* h_pulsearea_maxs2 = new TH1D("pulsearea_maxs2","pulsearea_maxs2",2000,0,2000);
//TH1D* h_pulsearea_subs2 = new TH1D("pulsearea_subs2","pulsearea_subs2",2000,0,2000);
//TH1D* h_singlePEprob_s2 = new TH1D("singlePEprob_s2","singlePEprob_s2",2,0,2);
//TH1D* h_otherprob_s2 = new TH1D("otherprob_s2","otherprob_s2",2,0,2);
//TH1D* h_others2prob_s2 = new TH1D("others2prob_s2","others2prob_s2",2,0,2);
//TH1D* h_areafrac5_s2 = new TH1D("areafrac5_s2","areafrac5_s2",300,0,3000);
//TH1D* h_areafrac95_s2 = new TH1D("areafrac95_s2","areafrac95_s2",600,0,6000);
//TH1D* h_peakamp_maxs2 = new TH1D("peakamp_maxs2","peakamp_maxs2",1000,0,50);
//TH1D* h_index_maxs2 = new TH1D("index_maxs2","index_maxs2",100,0,1000);
//TH1D* h_index_subs2 = new TH1D("index_subs2","index_subs2",100,0,1000);

//TH1D* h_s1prob_other = new TH1D("s1prob_other","s1prob_other",2,0,2);
//TH1D* h_s2prob_other = new TH1D("s2prob_other","s2prob_other",2,0,2);
//TH1D* h_pulsearea_other = new TH1D("pulsearea_other","pulsearea_other",2000,0,2000);
//TH1D* h_singlePEprob_other = new TH1D("singlePEprob_other","singlePEprob_other",2,0,2);
//TH1D* h_otherprob_other = new TH1D("otherprob_other","otherprob_other",2,0,2);
//TH1D* h_others2prob_other = new TH1D("others2prob_other","others2prob_other",2,0,2);
//TH1D* h_areafrac5_other = new TH1D("areafrac5_other","areafrac5_other",300,0,3000);
//TH1D* h_areafrac95_other = new TH1D("areafrac95_other","areafrac95_other",300,0,3000);
//TH1D* h_peakamp_other = new TH1D("peakamp_other","peakamp_other",1000,0,50);

//TH2D* h_area_versus_time_total = new TH2D("area_vs_time_total","area_vs_time_total",2000,0,2000,600,0,6000);
//TH2D* h_area_versus_time_s1 = new TH2D("area_vs_time_s1","area_vs_time_s1",2000,0,2000,600,0,6000);
//TH2D* h_area_versus_time_s2 = new TH2D("area_vs_time_s2","area_vs_time_total",2000,0,2000,600,0,6000);
//TH2D* h_area_versus_time_other = new TH2D("area_vs_time_other","area_vs_time_total",2000,0,2000,600,0,6000);

//TH2D* h_area_vs_length_maxs1 = new TH2D("area_vs_length_maxs1","area_vs_length_maxs1",2000,0,2000,600,0,6000);
//TH2D* h_area_vs_length_maxs2 = new TH2D("area_vs_length_maxs2","area_vs_length_maxs2",2000,0,2000,600,0,6000);
//TH2D* h_area_vs_length_subs1 = new TH2D("area_vs_length_subs1","area_vs_length_subs1",2000,0,2000,600,0,6000);
//TH2D* h_area_vs_length_subs2 = new TH2D("area_vs_length_subs2","area_vs_length_subs2",2000,0,2000,600,0,6000);
//TH2D* h_area_vs_length_other = new TH2D("area_vs_length_other","area_vs_length_other",2000,0,2000,600,0,6000);

TH2D* h_logxx_vs_logyy_total = new TH2D("logxx_vs_logyy_total","logxx_vs_logyy_total",500,-3.5,0.5,500,-5.5,0.5);
TH2D* h_logxx_vs_logyy_single = new TH2D("logxx_vs_logyy_single","logxx_vs_logyy_single",500,-3.5,0.5,500,-5.5,0.5);
TH2D* h_logxx_vs_logyy_multiple = new TH2D("logxx_vs_logyy_multiple","logxx_vs_logyy_multiple",500,-3.5,0.5,500,-5.5,0.5);

TH2D* h_maxs1_vs_maxs2_area = new TH2D("maxs1_vs_maxs2_area","maxs1_vs_maxs2_area",2000,0,2000,2000,0,2000);
TH2D* h_subs1_vs_subs2_area = new TH2D("h_subs1_vs_subs2_area","h_subs1_vs_subs2_area",2000,0,2000,2000,0,2000);

  //------------------------------------------------
  // Main event loop
  //------------------------------------------------
  const Int_t nevents = evt->chain->GetEntries();
  
  evt->GetEntry(0);
//  cout << nevents << endl;
  int processed_events=0;
   for (Int_t n=0; n<nevents; ++n) {
//	for (Int_t n=0; n<1000; ++n) {
	if (n%1000 == 0) cout << n << "/" << nevents << "\r" << flush;
	evt->GetEntry(n);

	int s1pulseID = -1;
	int s2pulseID = -1;
	
	int maxS1pulseID = -1;
	int maxS2pulseID = -1;
	int subS2pulseID = -1;
	int subS1pulseID = -1;
	float maxS1area = 0;
	float maxS2area = 0;
	float subS2area = 0;
	float subS1area = 0;


//	cout <<"npulses" << evt->nPulses_TPCHG << endl;
	for (int p=0; p<evt->nPulses_TPCHG; ++p) {
//		cout <<"s1prob" << evt->s1Probability_TPCHG[p] << endl;
		s1prob[0]->Fill(evt->s1Probability_TPCHG[p]);
//		cout <<"s2prob" << evt->s2Probability_TPCHG[p] << endl;
		s2prob[0]->Fill(evt->s2Probability_TPCHG[p]);
//		cout <<"singlepeprob" << evt->singlePEprobability_TPCHG[p] << endl;
//		h_singlePEprob_total->Fill(evt->singlePEprobability_TPCHG[p]);
//		cout <<"others2prob" << evt->otherS2Probability_TPCHG[p] << endl;
//		h_others2prob_total->Fill(evt->otherS2Probability_TPCHG[p]);
//		cout <<"otherprob" << evt->otherProbability_TPCHG[p] << endl;
//		h_otherprob_total->Fill(evt->otherProbability_TPCHG[p]);
//		cout <<"area" << evt->pulseArea_phd_TPCHG[p] << endl;
		pulseArea[0]->Fill(evt->pulseArea_phd_TPCHG[p]);
//		cout <<"area frac 5 " << evt->areaFractionTime5_ns_TPCHG[p] << endl;
		areaFrac5[0]->Fill(evt->areaFractionTime5_ns_TPCHG[p]);
//		cout <<"area frac 95" << evt->areaFractionTime95_ns_TPCHG[p] << endl;
		areaFrac95[0]->Fill(evt->areaFractionTime95_ns_TPCHG[p]);
//		cout << "peakamp" << evt->peakAmp_TPCHG[p] << endl;
		peakAmp[0]->Fill(evt->peakAmp_TPCHG[p]);
		area_vs_length[0]->Fill(evt->pulseArea_phd_TPCHG[p],evt->areaFractionTime95_ns_TPCHG[p]);
		index[0]->Fill(p);


		if (evt->s1Probability_TPCHG[p] == 1) {

//			cout <<"s1prob" << evt->s1Probability_TPCHG[p] << endl;
			s1prob[1]->Fill(evt->s1Probability_TPCHG[p]);
//			cout <<"s2prob" << evt->s2Probability_TPCHG[p] << endl;
			s2prob[1]->Fill(evt->s2Probability_TPCHG[p]);
//			cout <<"singlepeprob" << evt->singlePEprobability_TPCHG[p] << endl;
//			h_singlePEprob_s1->Fill(evt->singlePEprobability_TPCHG[p]);
//			cout <<"others2prob" << evt->otherS2Probability_TPCHG[p] << endl;
//			h_others2prob_s1->Fill(evt->otherS2Probability_TPCHG[p]);
//			cout <<"otherprob" << evt->otherProbability_TPCHG[p] << endl;
//			h_otherprob_s1->Fill(evt->otherProbability_TPCHG[p]);
//			cout <<"area" << evt->pulseArea_phd_TPCHG[p] << endl;
			pulseArea[1]->Fill(evt->pulseArea_phd_TPCHG[p]);
//			cout <<"area frac 5 " << evt->areaFractionTime5_ns_TPCHG[p] << endl;
			areaFrac5[1]->Fill(evt->areaFractionTime5_ns_TPCHG[p]);
//			cout <<"area frac 95" << evt->areaFractionTime95_ns_TPCHG[p] << endl;
			areaFrac95[1]->Fill(evt->areaFractionTime95_ns_TPCHG[p]);
//			cout << "peakamp" << evt->peakAmp_TPCHG[p] << endl;
			peakAmp[1]->Fill(evt->peakAmp_TPCHG[p]);
			area_vs_length[1]->Fill(evt->pulseArea_phd_TPCHG[p],evt->areaFractionTime95_ns_TPCHG[p]);
//			h_area_vs_length[1]->Fill(evt->pulseArea_phd_TPCHG[p],(evt->pulseEndTime_ns_TPCHG[p]) - (evt->pulseStartTime_ns_TPCHG[p]));
			index[1]->Fill(p);
		}
		if (evt->s2Probability_TPCHG[p] == 1) {

//			cout <<"s1prob" << evt->s1Probability_TPCHG[p] << endl;
			s1prob[2]->Fill(evt->s1Probability_TPCHG[p]);
//			cout <<"s2prob" << evt->s2Probability_TPCHG[p] << endl;
			s2prob[2]->Fill(evt->s2Probability_TPCHG[p]);
//			cout <<"singlepeprob" << evt->singlePEprobability_TPCHG[p] << endl;
//			h_singlePEprob_s2->Fill(evt->singlePEprobability_TPCHG[p]);
//			cout <<"others2prob" << evt->otherS2Probability_TPCHG[p] << endl;
//			h_others2prob_s2->Fill(evt->otherS2Probability_TPCHG[p]);
//			cout <<"otherprob" << evt->otherProbability_TPCHG[p] << endl;
//			h_otherprob_s2->Fill(evt->otherProbability_TPCHG[p]);
//			cout <<"area" << evt->pulseArea_phd_TPCHG[p] << endl;
			pulseArea[2]->Fill(evt->pulseArea_phd_TPCHG[p]);
//			cout <<"area frac 5 " << evt->areaFractionTime5_ns_TPCHG[p] << endl;
			areaFrac5[2]->Fill(evt->areaFractionTime5_ns_TPCHG[p]);
//			cout <<"area frac 95" << evt->areaFractionTime95_ns_TPCHG[p] << endl;
			areaFrac95[2]->Fill(evt->areaFractionTime95_ns_TPCHG[p]);
//			cout << "peakamp" << evt->peakAmp_TPCHG[p] << endl;
			peakAmp[2]->Fill(evt->peakAmp_TPCHG[p]);
			area_vs_length[2]->Fill(evt->pulseArea_phd_TPCHG[p],evt->areaFractionTime95_ns_TPCHG[p]);
//			h_area_vs_length_s2->Fill(evt->pulseArea_phd_TPCHG[p],(evt->pulseEndTime_ns_TPCHG[p]) - (evt->pulseStartTime_ns_TPCHG[p]));
			index[2]->Fill(p);
		}
		if (evt->s1Probability_TPCHG[p] == 0 && evt->s2Probability_TPCHG[p] == 0) {

//			cout <<"s1prob" << evt->s1Probability_TPCHG[p] << endl;
			s1prob[3]->Fill(evt->s1Probability_TPCHG[p]);
//			cout <<"s2prob" << evt->s2Probability_TPCHG[p] << endl;
			s2prob[3]->Fill(evt->s2Probability_TPCHG[p]);
//			cout <<"singlepeprob" << evt->singlePEprobability_TPCHG[p] << endl;
//			h_singlePEprob_other->Fill(evt->singlePEprobability_TPCHG[p]);
//			cout <<"others2prob" << evt->otherS2Probability_TPCHG[p] << endl;
//			h_others2prob_other->Fill(evt->otherS2Probability_TPCHG[p]);
//			cout <<"otherprob" << evt->otherProbability_TPCHG[p] << endl;
//			h_otherprob_other->Fill(evt->otherProbability_TPCHG[p]);
//			cout <<"area" << evt->pulseArea_phd_TPCHG[p] << endl;
			pulseArea[3]->Fill(evt->pulseArea_phd_TPCHG[p]);
//			cout <<"area frac 5 " << evt->areaFractionTime5_ns_TPCHG[p] << endl;
			areaFrac5[3]->Fill(evt->areaFractionTime5_ns_TPCHG[p]);
//			cout <<"area frac 95" << evt->areaFractionTime95_ns_TPCHG[p] << endl;
			areaFrac95[3]->Fill(evt->areaFractionTime95_ns_TPCHG[p]);
//			cout << "peakamp" << evt->peakAmp_TPCHG[p] << endl;
			peakAmp[3]->Fill(evt->peakAmp_TPCHG[p]);
			area_vs_length[3]->Fill(evt->pulseArea_phd_TPCHG[p],evt->areaFractionTime95_ns_TPCHG[p]);
//			h_area_vs_length_other->Fill(evt->pulseArea_phd_TPCHG[p],(evt->pulseEndTime_ns_TPCHG[p]) - (evt->pulseStartTime_ns_TPCHG[p]));
			index[3]->Fill(p);
		}





//		cout << "" << endl;

		if (evt->s1Probability_TPCHG[p] == 1 && evt->pulseArea_phd_TPCHG[p] > maxS1area) {
			maxS1pulseID = p;
//			cout << "changed maxs1" << endl;
			maxS1area = evt->pulseArea_phd_TPCHG[p];
		}
		if (evt->s2Probability_TPCHG[p] == 1 && evt->pulseArea_phd_TPCHG[p] > maxS2area) {
			maxS2pulseID = p;
//			cout << "changed maxs2" << endl;
			maxS2area = evt->pulseArea_phd_TPCHG[p];
		}
	}
	
	
	if (maxS1pulseID != -1 && maxS2pulseID != -1) {
		h_maxs1_vs_maxs2_area->Fill(maxS1area,maxS2area);
	}
//	h_peakamp_maxs1->Fill(evt->peakAmp_TPCHG[maxS1pulseID]);
//	h_peakamp_maxs2->Fill(evt->peakAmp_TPCHG[maxS2pulseID]);

	if (maxS1pulseID != -1) {
                pulseArea[4]->Fill(maxS1area);
                index[4]->Fill(maxS1pulseID);
		areaFrac5[4]->Fill(evt->areaFractionTime5_ns_TPCHG[maxS1pulseID]);
		areaFrac95[4]->Fill(evt->areaFractionTime95_ns_TPCHG[maxS1pulseID]);
		const float end = evt->areaFractionTime95_ns_TPCHG[maxS1pulseID];
		const float start = evt->areaFractionTime5_ns_TPCHG[maxS1pulseID];
		const float diff = end - start;
		area_vs_length[4]->Fill(maxS1area,diff);
		peakAmp[4]->Fill(evt->peakAmp_TPCHG[maxS1pulseID]);
		s1prob[4]->Fill(evt->s1Probability_TPCHG[maxS1pulseID]);
		s2prob[4]->Fill(evt->s2Probability_TPCHG[maxS1pulseID]);
	}
	if (maxS2pulseID != -1) {
                pulseArea[5]->Fill(maxS2area);
                index[5]->Fill(maxS2pulseID);
		const float end = evt->areaFractionTime95_ns_TPCHG[maxS2pulseID];
		const float start = evt->areaFractionTime5_ns_TPCHG[maxS2pulseID];
		const float diff = end - start;
		area_vs_length[5]->Fill(maxS2area,diff);
		peakAmp[5]->Fill(evt->peakAmp_TPCHG[maxS2pulseID]);
		areaFrac5[5]->Fill(evt->areaFractionTime5_ns_TPCHG[maxS2pulseID]);
		areaFrac95[5]->Fill(evt->areaFractionTime95_ns_TPCHG[maxS2pulseID]);
		s1prob[5]->Fill(evt->s1Probability_TPCHG[maxS2pulseID]);
		s2prob[5]->Fill(evt->s2Probability_TPCHG[maxS2pulseID]);
	}

	for (int p=0; p<evt->nPulses_TPCHG; ++p) {

		if(evt->s2Probability_TPCHG[p] == 1 && evt->pulseArea_phd_TPCHG[p] > subS2area && evt->pulseArea_phd_TPCHG[p] < maxS2area) {
			subS2pulseID = p;
			subS2area = evt->pulseArea_phd_TPCHG[p];
		}
                if(evt->s1Probability_TPCHG[p] == 1 && evt->pulseArea_phd_TPCHG[p] > subS1area && evt->pulseArea_phd_TPCHG[p] < maxS1area) {
                        subS1pulseID = p;
                        subS1area = evt->pulseArea_phd_TPCHG[p];
                }
	}
	if (maxS1pulseID == -1 || maxS2pulseID == -1) {
	continue;
	}
        if (subS1pulseID != -1) {
	        pulseArea[6]->Fill(subS1area);
	        index[6]->Fill(subS1pulseID);
                const float end = evt->areaFractionTime95_ns_TPCHG[subS1pulseID];
                const float start = evt->areaFractionTime5_ns_TPCHG[subS1pulseID];
                const float diff = end - start;
                area_vs_length[6]->Fill(subS1area,diff);
		peakAmp[6]->Fill(evt->peakAmp_TPCHG[subS1pulseID]);
		s1prob[6]->Fill(evt->s1Probability_TPCHG[subS1pulseID]);
		s2prob[6]->Fill(evt->s2Probability_TPCHG[subS1pulseID]);
		areaFrac5[6]->Fill(evt->areaFractionTime5_ns_TPCHG[subS1pulseID]);
		areaFrac95[6]->Fill(evt->areaFractionTime95_ns_TPCHG[subS1pulseID]);
        }
        if (subS2pulseID != -1) {
	        pulseArea[7]->Fill(subS2area);
	        index[7]->Fill(subS2pulseID);
                const float end = evt->areaFractionTime95_ns_TPCHG[subS2pulseID];
                const float start = evt->areaFractionTime5_ns_TPCHG[subS2pulseID];
                const float diff = end - start;
                area_vs_length[7]->Fill(subS2area,diff);
		peakAmp[7]->Fill(evt->peakAmp_TPCHG[subS2pulseID]);
		s1prob[7]->Fill(evt->s1Probability_TPCHG[subS2pulseID]);
                s2prob[7]->Fill(evt->s2Probability_TPCHG[subS2pulseID]);
                areaFrac5[7]->Fill(evt->areaFractionTime5_ns_TPCHG[subS2pulseID]);
                areaFrac95[7]->Fill(evt->areaFractionTime95_ns_TPCHG[subS2pulseID]);


        }
	if (subS1pulseID != -1 && subS2pulseID != -1) {
		h_subs1_vs_subs2_area->Fill(subS1area,subS2area);
	}


	for (int p=0; p<maxS1pulseID; ++p) {

		if (evt->s2Probability_TPCHG[p] == 1 ||
		    evt->singlePEprobability_TPCHG[p] == 1 ||
		    evt->otherS2Probability_TPCHG[p] == 1 ||
		    evt->otherProbability_TPCHG[p] == 1) {
			continue;
		}

	}

	const float x1 = 1e-1;
	const float x2 = 3e-2;
	const float y1 = 5e-4;
	const float y2 = 4e-2;
	const float logx1 = TMath::Log10(x1);
	const float logx2 = TMath::Log10(x2);
	const float logy1 = TMath::Log10(y1);
	const float logy2 = TMath::Log10(y2);
	const float m = (logy2-logy1)/(logx2-logx1);
	const float b = logy1-m*logx1;
	const float xx = subS2area/maxS2area;
	const float length = evt->areaFractionTime95_ns_TPCHG[subS2pulseID] - evt->areaFractionTime5_ns_TPCHG[subS2pulseID];
	const float height = evt->peakAmp_TPCHG[subS2pulseID];
	const float yy = height/length;

//	cout << TMath::Log10(xx) << "	" << TMath::Log10(yy) << endl;
	h_logxx_vs_logyy_total->Fill(TMath::Log10(xx),TMath::Log10(yy));

	if (maxS2pulseID != -1 && subS2pulseID == -1) {
		s1pulseID = maxS1pulseID;
		s2pulseID = maxS2pulseID;
//		fill single scatters here
		h_singlescatter_s1Area_phd->Fill(evt->pulseArea_phd_TPCHG[maxS1pulseID]);
		h_logxx_vs_logyy_single->Fill(TMath::Log10(xx),TMath::Log10(yy));
		continue;
	}


	if (subS2pulseID == -1) {
		continue;
	}

	if (TMath::Log10(yy) > m * TMath::Log10(xx) + b) {
		h_multiplescatter_s1Area_phd->Fill(evt->pulseArea_phd_TPCHG[maxS1pulseID]);
		h_logxx_vs_logyy_multiple->Fill(TMath::Log10(xx),TMath::Log10(yy));
		continue;
	}

	s1pulseID = maxS1pulseID;
	s2pulseID = maxS2pulseID;
	//its a single scatter if it hasn't brokent the loop yet
	h_singlescatter_s1Area_phd->Fill(evt->pulseArea_phd_TPCHG[maxS1pulseID]);
	h_logxx_vs_logyy_single->Fill(TMath::Log10(xx),TMath::Log10(yy));


//   for (int lop = 0; lop < evt->singlescatter_s1Area_phd.size(); lop++){
//    if(evt->nSingleScatters>0){
//      if (evt->nSingleScatters>1) cout<<"WARNING, nSingleScatters>1"<<endl;
//      h_singlescatter_s1Area_phd->Fill(evt->singlescatter_s1Area_phd[lop]);
//    }
//   }
//   for (int lop = 0; lop < evt->multiplescatter_s1Area_phd.size(); lop++){
//   if(evt->nMultipleScatters>0){
//      if (evt->nMultipleScatters>1) cout<<"It's more than 1 in the nMultipleScatters"<<endl;
//      h_multiplescatter_s1Area_phd->Fill(evt->multiplescatter_s1Area_phd[lop]);
//    }//int nevents
//    }

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
//		h_pulsearea_yess2_TPCHG->Fill(evt->pulseArea_phd_TPCHG[lop]);
		}
	else {
//		h_pulsearea_nos2_TPCHG->Fill(evt->pulseArea_phd_TPCHG[lop]);
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


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


  //TString sample="background";
  //TString sample="NA22";
  //TString sample="Kr83";
  TString sample="AmLi";
  TString outname  = sample+".root";
  TString txtFileList = "RQfile.list."+sample;
  
  outfile = new TFile(outname, "recreate");
  cout << "Writing output to: "<<outfile->GetName()<<endl;

  TChain* chain = new TChain("Events", "Events");
  load_chain(txtFileList, chain);
  MyEvent* evt = new MyEvent();
  evt->LoadBranches(chain);


  std::vector<std::vector<unsigned long long>> livetime;
  int firstPulseStart = 0;
  int lastPulseEnd = 0;
  int duration = 0;
  
  TH1D* h_triggerTimeStamp_ns = new TH1D("triggerTimeStamp_ns", "triggerTimeStamp_ns", 100000, 0, 2.5e+08);
  TH1D* h_triggerTimeStamp_s = new TH1D("triggerTimeStamp_s", "triggerTimeStamp_s", 10000, 0, 10000);
  //OD
  TH1D* h_triggerTimeStamp_s_ODLG = new TH1D("triggerTimeStamp_s_ODLG", "triggerTimeStamp_s_ODLG", 10000, 0, 10000);
  TH1F* h_pulse_duration_ms_ODLG = new TH1F("pulse_duration_ms_ODLG", "pulse_duration_ms_ODLG", 300, 0, 10);
  TH1F* h_pulseArea_phd_ODLG = new TH1F("pulseArea_phd_ODLG", "pulseArea_phd_ODLG", 200, 0, 1500);
  TH1F* h_pulseArea_phd_ODLG_sum = new TH1F("pulseArea_phd_ODLG_sum", "pulseArea_phd_ODLG_sum", 100, 0, 1500);
  TH1F* h_positiveArea_phd_ODLG = new TH1F("positiveArea_phd_ODLG", "positiveArea_phd_ODLG", 300, 0, 300);
  TH1F* h_negativeArea_phd_ODLG = new TH1F("negativeArea_phd_ODLG", "negativeArea_phd_ODLG", 300, -100, 10);
  TH1F* h_nPulses_ODLG = new TH1F("nPulses_ODLG", "nPulses_ODLG", 100, 0, 100);
  TH1F* h_chPulseArea_phd_ODLG = new TH1F("chPulseArea_phd_ODLG", "chPulseArea_phd_ODLG", 20, 0, 20);
  TH1F* h_chPeakTime_ns_ODLG = new TH1F("chPeakTime_ns_ODLG", "chPeakTime_ns_ODLG", 700, 0, 700);
  TH1F* h_channelID_ODLG = new TH1F("channelID_ODLG", "channelID_ODLG", 1000, 0, 1000);
  TH1F* h_photonCount_ODLG = new TH1F("photonCount_ODLG", "photonCount_ODLG", 200, 0, 1500);
  TH1F* h_pulseStartTime_ns_ODLG = new TH1F("pulseStartTime_ns_ODLG", "pulseStartTime_ns_ODLG", 200, 0, 1500);
  TH1F* h_pulseEndTime_ns_ODLG = new TH1F("pulseEndTime_ns_ODLG", "pulseEndTime_ns_ODLG", 200, 0, 2000);
  TH1F* h_peakAmp_ODLG = new TH1F("peakAmp_ODLG", "peakAmp_ODLG", 20, 0, 20);
  TH1F* h_peakTime_ns_ODLG = new TH1F("peakTime_ns_ODLG", "peakTime_ns_ODLG", 200, 0, 1500);
  TH1F* h_areaFractionTime5_ns_ODLG = new TH1F("areaFractionTime5_ns_ODLG", "areaFractionTime5_ns_ODLG", 200, 0, 1500);
  TH1F* h_areaFractionTime25_ns_ODLG = new TH1F("areaFractionTime25_ns_ODLG", "areaFractionTime25_ns_ODLG", 200, 0, 1500);
  TH1F* h_areaFractionTime50_ns_ODLG = new TH1F("areaFractionTime50_ns_ODLG", "areaFractionTime50_ns_ODLG", 200, 0, 1500);
  TH1F* h_areaFractionTime75_ns_ODLG = new TH1F("areaFractionTime75_ns_ODLG", "areaFractionTime75_ns_ODLG", 200, 0, 1500);
  TH1F* h_areaFractionTime95_ns_ODLG = new TH1F("areaFractionTime95_ns_ODLG", "areaFractionTime95_ns_ODLG", 200, 0, 1500);
  TH1F* h_promptFraction50ns_ODLG = new TH1F("promptFraction50ns_ODLG", "promptFraction50ns_ODLG", 20, 0, 20);
  TH1F* h_rmsWidth_ns_ODLG = new TH1F("rmsWidth_ns_ODLG", "rmsWidth_ns_ODLG", 200, 0, 2000);
  TH1F* h_coincidence_ODLG = new TH1F("coincidence_ODLG", "coincidence_ODLG", 20, 0, 20);
  TH1F* h_s1Probability_ODLG = new TH1F("s1Probability_ODLG", "s1Probability_ODLG", 10, 0, 1);
  TH1F* h_s2Probability_ODLG = new TH1F("s2Probability_ODLG", "s2Probability_ODLG", 10, 0, 1);
  TH1F* h_singleElectronProbability_ODLG = new TH1F("singleElectronProbability_ODLG", "singleElectronProbability_ODLG", 10, 0, 1);
  TH1F* h_singlePEprobability_ODLG = new TH1F("singlePEprobability_ODLG", "singlePEprobability_ODLG", 10, 0, 1);
  TH1F* h_otherProbability_ODLG = new TH1F("otherProbability_ODLG", "otherProbability_ODLG", 10, 0, 1);
  TH1F* h_otherS2Probability_ODLG = new TH1F("otherS2Probability_ODLG", "otherS2Probability_ODLG", 10, 0, 1);


  TH2F* h_pulse_area_width_ODLG = new TH2F("pulse_area_width_ODLG", "pulse_area_width_ODLG", 200, 0, 1500, 500, 0, 3000);
  TH2F* h_pulse_start_area_ODLG = new TH2F("pulse_start_area_ODLG", "pulsestart_area_ODLG", 200, 0, 4000, 200, 0, 1500);

  //TPC
  TH1D* h_triggerTimeStamp_s_TPCLG = new TH1D("triggerTimeStamp_s_TPCLG", "triggerTimeStamp_s_TPCLG", 10000, 0, 10000);
  TH1F* h_pulse_duration_ms_TPCLG = new TH1F("pulse_duration_ms_TPCLG", "pulse_duration_ms_TPCLG", 300, 0, 10);
  TH1F* h_pulseArea_phd_TPCLG = new TH1F("pulseArea_phd_TPCLG", "pulseArea_phd_TPCLG", 200, 0, 1500);
  TH1F* h_pulseArea_phd_TPCLG_sum = new TH1F("pulseArea_phd_TPCLG_sum", "pulseArea_phd_TPCLG_sum", 100, 0, 1500);
  TH1F* h_positiveArea_phd_TPCLG = new TH1F("positiveArea_phd_TPCLG", "positiveArea_phd_TPCLG", 300, 0, 300);
  TH1F* h_negativeArea_phd_TPCLG = new TH1F("negativeArea_phd_TPCLG", "negativeArea_phd_TPCLG", 300, -100, 10);
  TH1F* h_nPulses_TPCLG = new TH1F("nPulses_TPCLG", "nPulses_TPCLG", 100, 0, 100);
  TH1F* h_chPulseArea_phd_TPCLG = new TH1F("chPulseArea_phd_TPCLG", "chPulseArea_phd_TPCLG", 20, 0, 20);
  TH1F* h_chPeakTime_ns_TPCLG = new TH1F("chPeakTime_ns_TPCLG", "chPeakTime_ns_TPCLG", 700, 0, 700);
  TH1F* h_channelID_TPCLG = new TH1F("channelID_TPCLG", "channelID_TPCLG", 1000, 0, 1000);
  TH1F* h_photonCount_TPCLG = new TH1F("photonCount_TPCLG", "photonCount_TPCLG", 200, 0, 1500);
  TH1F* h_pulseStartTime_ns_TPCLG = new TH1F("pulseStartTime_ns_TPCLG", "pulseStartTime_ns_TPCLG", 200, 0, 1500);
  TH1F* h_pulseEndTime_ns_TPCLG = new TH1F("pulseEndTime_ns_TPCLG", "pulseEndTime_ns_TPCLG", 200, 0, 2000);
  TH1F* h_peakAmp_TPCLG = new TH1F("peakAmp_TPCLG", "peakAmp_TPCLG", 20, 0, 20);
  TH1F* h_peakTime_ns_TPCLG = new TH1F("peakTime_ns_TPCLG", "peakTime_ns_TPCLG", 200, 0, 1500);
  TH1F* h_areaFractionTime5_ns_TPCLG = new TH1F("areaFractionTime5_ns_TPCLG", "areaFractionTime5_ns_TPCLG", 200, 0, 1500);
  TH1F* h_areaFractionTime25_ns_TPCLG = new TH1F("areaFractionTime25_ns_TPCLG", "areaFractionTime25_ns_TPCLG", 200, 0, 1500);
  TH1F* h_areaFractionTime50_ns_TPCLG = new TH1F("areaFractionTime50_ns_TPCLG", "areaFractionTime50_ns_TPCLG", 200, 0, 1500);
  TH1F* h_areaFractionTime75_ns_TPCLG = new TH1F("areaFractionTime75_ns_TPCLG", "areaFractionTime75_ns_TPCLG", 200, 0, 1500);
  TH1F* h_areaFractionTime95_ns_TPCLG = new TH1F("areaFractionTime95_ns_TPCLG", "areaFractionTime95_ns_TPCLG", 200, 0, 1500);
  TH1F* h_promptFraction50ns_TPCLG = new TH1F("promptFraction50ns_TPCLG", "promptFraction50ns_TPCLG", 20, 0, 20);
  TH1F* h_rmsWidth_ns_TPCLG = new TH1F("rmsWidth_ns_TPCLG", "rmsWidth_ns_TPCLG", 200, 0, 2000);
  TH1F* h_coincidence_TPCLG = new TH1F("coincidence_TPCLG", "coincidence_TPCLG", 20, 0, 20);
  TH1F* h_s1Probability_TPCLG = new TH1F("s1Probability_TPCLG", "s1Probability_TPCLG", 10, 0, 1);
  TH1F* h_s2Probability_TPCLG = new TH1F("s2Probability_TPCLG", "s2Probability_TPCLG", 10, 0, 1);
  TH1F* h_singleElectronProbability_TPCLG = new TH1F("singleElectronProbability_TPCLG", "singleElectronProbability_TPCLG", 10, 0, 1);
  TH1F* h_singlePEprobability_TPCLG = new TH1F("singlePEprobability_TPCLG", "singlePEprobability_TPCLG", 10, 0, 1);
  TH1F* h_otherProbability_TPCLG = new TH1F("otherProbability_TPCLG", "otherProbability_TPCLG", 10, 0, 1);
  TH1F* h_otherS2Probability_TPCLG = new TH1F("otherS2Probability_TPCLG", "otherS2Probability_TPCLG", 10, 0, 1);



  TH2F* h_pulse_area_width_TPCLG = new TH2F("pulse_area_width_TPCLG", "pulse_area_width_TPCLG", 200, 0, 1500, 500, 0, 3000);
  TH2F* h_pulse_start_area_TPCLG = new TH2F("pulse_start_area_TPCLG", "pulse_start_area_TPCLG", 200, 0, 4000, 200, 0, 1500);


  //for position reconstruction
  bool map_done=0;
  TH2F* h_map=new TH2F("map","map",20, 0.5,20.5, 6,0.5,6.5);
  // TH2F* h_z_phd=new TH2F("z_phd","z_phid",6,0.5,6.5, 200,0,1000);
  // TH2F* h_i_phd=new TH2F("i_phd","i_phid",20, 0.5,20.5, 200,0,1000);
  // TH2F* h_i_z=new TH2F("i_z","i_z", 20, -10, 10, 6, -3,3);

  TH2F* h_npulses_phd=new TH2F("npulses_phd","npulses_phd",100, 9, 130, 50,0,100);
  TH2F* h_chPktime_phd=new TH2F("chPktime_phd","chPktime_phd",200, 0, 1000, 200, 0, 100);
  TH2F* h_dz_phd=new TH2F("z_phd","z_phd",6, -225-35, 225-35, 200,0,100);
  TH2F* h_dphi_phd=new TH2F("i_phd","i_phd", 20, -180, 180, 200,0,100);
  TH2F* h_dphi_dz=new TH2F("dphi_dz","dphi_dz", 20, -180, 180, 6, -225-35, 225-35);

  std::vector <TH1F*> hvec_chPulseArea_phd_ODLG;
  for (int i=0; i<120; ++i){
    TString name="chPulseArea_phd_ODLG_"; name+=i;
    TH1F* h1 = new TH1F(name,name,40,0,40);
    hvec_chPulseArea_phd_ODLG.push_back(h1);
  }
  

  //------------------------------------------------
  // Main event loop
  //------------------------------------------------
  const Int_t nevents = evt->chain->GetEntries();
  
  evt->GetEntry(0);
  double first_ns=evt->triggerTimeStamp_ns;
  double first_s=evt->triggerTimeStamp_s;
  evt->GetEntry(nevents-1);
  double last_ns=evt->triggerTimeStamp_ns;
  double last_s=evt->triggerTimeStamp_s;
  cout<<fixed<<"duration of run: "<<(last_s-first_s)<< "sec "<<endl;

  float min_time_ns=-99;
  float max_time_ns=-99;
  float min_time_s=-99;
  float max_time_s=-99;
  
  min_time_ns=evt->triggerTimeStamp_ns;
  max_time_ns=evt->triggerTimeStamp_ns;
  min_time_s=evt->triggerTimeStamp_s;
  max_time_s=evt->triggerTimeStamp_s;
  int processed_events=0;
  double sum_pulse_time=0;
  int number_of_pulses=0;
  for (Int_t n=0; n<nevents; ++n) {
    //for (Int_t n=0; n<1000; ++n) {
    if (n%1000 == 0) cout << "Processing "<< n << "/"<<nevents<<endl;
    processed_events++;
    
    if(evt->triggerTimeStamp_s<first_s) first_s=evt->triggerTimeStamp_s;
    if(evt->triggerTimeStamp_s>last_s) last_s=evt->triggerTimeStamp_s;


    evt->GetEntry(n);
    h_triggerTimeStamp_ns->Fill(evt->triggerTimeStamp_ns);
    h_triggerTimeStamp_s->Fill(evt->triggerTimeStamp_s);

    //check timing
    firstPulseStart = -99;
    lastPulseEnd = -99;
    if(evt->nPulses_ODLG!=0){
      if(evt->pulseStartTime_ns_ODLG[0] < firstPulseStart || firstPulseStart==-99) firstPulseStart = evt->pulseStartTime_ns_ODLG[0];
      if(evt->pulseEndTime_ns_ODLG[(evt->nPulses_ODLG)-1] > lastPulseEnd || lastPulseEnd>-99) lastPulseEnd = evt->pulseEndTime_ns_ODLG[(evt->nPulses_ODLG)-1];
      h_triggerTimeStamp_s_ODLG->Fill(evt->triggerTimeStamp_s);
    }
    duration = lastPulseEnd - firstPulseStart;
    sum_pulse_time+=duration;
    h_pulse_duration_ms_ODLG->Fill(duration/1E6); //fill in ms

    //*****************
    // OD Low gain
    //*****************
    bool started_filling=0;
    number_of_pulses+=evt->nPulses_ODLG;
    h_nPulses_ODLG->Fill(evt->nPulses_ODLG);
    float all_chPulseArea_phd=0;
    float all_chPhotonCount_phd=0;
    
    //BP play with scatter
    // if( evt->nSingleScatters !=0 ){
    //   cout<<"nMultipleScatters "<<evt->nMultipleScatters <<" nSingleScatters "<< evt->nSingleScatters<<endl;
    //   cout<<"number of pulses in OD: "<<evt->nODPromptPulses[0]<<endl;
    //   //      cout<<"number of delayed pulses in OD: "<<evt->nODDelayedPulses[0]<<endl;
    // }
    
    // if( evt->nSingleScatters !=1 )        
    //   continue;

    for(int i=0; i<evt->nPulses_ODLG; ++i){
      duration = 0;
      h_pulseArea_phd_ODLG->Fill(evt->pulseArea_phd_ODLG[i]);
      h_photonCount_ODLG->Fill(evt->photonCount_ODLG[i]);
      h_pulseStartTime_ns_ODLG->Fill(evt->pulseStartTime_ns_ODLG[i]);
      h_pulseEndTime_ns_ODLG->Fill(evt->pulseEndTime_ns_ODLG[i]);
      h_peakAmp_ODLG->Fill(evt->peakAmp_ODLG[i]);
      h_peakTime_ns_ODLG->Fill(evt->peakTime_ns_ODLG[i]);
      h_areaFractionTime5_ns_ODLG->Fill(evt->areaFractionTime5_ns_ODLG[i]);
      h_areaFractionTime25_ns_ODLG->Fill(evt->areaFractionTime25_ns_ODLG[i]);
      h_areaFractionTime50_ns_ODLG->Fill(evt->areaFractionTime50_ns_ODLG[i]);
      h_areaFractionTime75_ns_ODLG->Fill(evt->areaFractionTime75_ns_ODLG[i]);
      h_areaFractionTime95_ns_ODLG->Fill(evt->areaFractionTime95_ns_ODLG[i]);
      h_promptFraction50ns_ODLG->Fill(evt->promptFraction50ns_ODLG[i]);
      h_rmsWidth_ns_ODLG->Fill(evt->rmsWidth_ns_ODLG[i]);
      h_coincidence_ODLG->Fill(evt->coincidence_ODLG[i]);
      h_s1Probability_ODLG->Fill(evt->s1Probability_ODLG[i]);
      h_s2Probability_ODLG->Fill(evt->s2Probability_ODLG[i]);
      h_singlePEprobability_ODLG->Fill(evt->singlePEprobability_ODLG[i]);
      h_singleElectronProbability_ODLG->Fill(evt->singlePEprobability_ODLG[i]);
      h_otherProbability_ODLG->Fill(evt->otherProbability_ODLG[i]);
      h_otherS2Probability_ODLG->Fill(evt->otherS2Probability_ODLG[i]);


      for (unsigned int ii=0; ii<evt->chPulseArea_phd_ODLG[i].size(); ii++) {
	h_channelID_ODLG->Fill(evt->channelID_ODLG[i][ii]-1000);
	int pmtid =(evt->channelID_ODLG[i][ii])-1801;
	hvec_chPulseArea_phd_ODLG.at(pmtid)->Fill(evt->chPulseArea_phd_ODLG[i][ii]);

	if (!map_done && evt->eventID==8) {
	  h_map->SetBinContent(pmtid/6+1, pmtid%6+1, ( h_map->GetBinContent(pmtid/6+1, pmtid%6+1) + evt->chPulseArea_phd_ODLG[i][ii]));
	  started_filling=true;
	}
	//	if(evt->chPulseArea_phd_ODLG[i][ii]>3 && evt->multiplePEprobability_TPCLG[i]==1){
	if(evt->chPulseArea_phd_ODLG[i][ii]>3){
	  all_chPulseArea_phd+=evt->chPulseArea_phd_ODLG[i][ii];
	  h_chPulseArea_phd_ODLG->Fill(evt->chPulseArea_phd_ODLG[i][ii]);
	  h_chPeakTime_ns_ODLG->Fill(evt->chPeakTime_ns_ODLG[i][ii]);
	}
      }//run over phd channels

      //now try to find where the light actually is //bp
      float max_area_phd=0;
      float max_so_far=0;
      int max_area_phd_pmt=-1;
      if (evt->chPulseArea_phd_ODLG[i].size()<5) continue;
      for (unsigned int ii=0; ii<evt->chPulseArea_phd_ODLG[i].size(); ii++) {
	int pmtid =(evt->channelID_ODLG[i][ii])-1801;
	//	if(evt->chPulseArea_phd_ODLG[i][ii]>3){
	  max_so_far=evt->chPulseArea_phd_ODLG[i][ii];
	  //	}
	if(max_area_phd<=max_so_far) {
	  max_area_phd=max_so_far;
	  max_area_phd_pmt=pmtid;
	}
	//      cout<<"ii "<<ii<<", phd: "<<evt->chPulseArea_phd_ODLG[i][ii] << ", ID: "<< (evt->channelID_ODLG[i][ii])-1801<<", z: "<<pmtid%6+1 << ", l: "<<pmtid/6+1<<endl;
      }
      //      cout<<" max_area_phd_pmt :"<<max_area_phd_pmt<<", max_area_phd: "<<max_area_phd<<endl;
      
      //if (evt->chPulseArea_phd_ODLG[i].size()<10) continue;
      for (unsigned int ii=0; ii<evt->chPulseArea_phd_ODLG[i].size(); ii++) {
	int pmtid =(evt->channelID_ODLG[i][ii])-1801;
	//if (pmtid==max_area_phd_pmt) continue;
	// float pmtid_z=pmtid%6+1;
	// float pmtid_i=pmtid/6+1;
	// float max_pmt_z=max_area_phd_pmt%6+1;
	// float max_pmt_i=max_area_phd_pmt/6+1;
	float pmtid_z=(pmtid%6+1)*77+25;
	float pmtid_i=(pmtid/6+1)*18+9;
	float max_pmt_z=(max_area_phd_pmt%6+1)*77+25;
	float max_pmt_i=(max_area_phd_pmt/6+1)*18+9;

	//cout<<" pmtid: "<<pmtid<<", max pmt id: "<<max_area_phd_pmt<<", pmtid_z:"<<pmtid_z<<", max_pmt_z: "<<max_pmt_z<<"/"<<pmtid%6+1<<", diff:"<<max_pmt_z-pmtid_z<<"/"<<max_area_phd_pmt%6+1<<endl;

	h_npulses_phd->Fill(evt->chPulseArea_phd_ODLG[i].size(), evt->chPulseArea_phd_ODLG[i][ii]);
	h_chPktime_phd->Fill(evt->chPeakTime_ns_ODLG[i][ii], evt->chPulseArea_phd_ODLG[i][ii]);
	h_dz_phd->Fill(max_pmt_z-pmtid_z, evt->chPulseArea_phd_ODLG[i][ii]);
	h_dphi_phd->Fill(max_pmt_i-pmtid_i, evt->chPulseArea_phd_ODLG[i][ii]);
	h_dphi_dz->Fill(max_pmt_i-pmtid_i, max_pmt_z-pmtid_z);
	//	if ( evt->chPulseArea_phd_ODLG[i][ii]>3 )
	  //  cout<<"event ID "<< evt->eventID << " ii "<<ii<<", phd: "<<evt->chPulseArea_phd_ODLG[i][ii] << ", ID: "<< (evt->channelID_ODLG[i][ii])-1801<<endl; //BP
	//	cout<<"ii "<<ii<<", phd: "<<evt->chPulseArea_phd_ODLG[i][ii] << ", ID: "<< (evt->channelID_ODLG[i][ii])-1801<<", dz: "<<max_pmt_z-pmtid_z << ", dl: "<<max_pmt_i-pmtid_i<<endl;
      }
      

      //end position reco

      duration+=evt->pulseEndTime_ns_ODLG[i]-evt->pulseStartTime_ns_ODLG[i];
      h_pulse_area_width_ODLG->Fill(evt->pulseArea_phd_ODLG[i], evt->pulseEndTime_ns_ODLG[i]-evt->pulseStartTime_ns_ODLG[i]);
      h_pulse_start_area_ODLG->Fill(evt->pulseStartTime_ns_ODLG[i], evt->pulseArea_phd_ODLG[i]);
    }//Pulses
    if ( started_filling ) map_done=1;
    h_pulseArea_phd_ODLG_sum->Fill(all_chPulseArea_phd);
    //    h_pulse_area_width_ODLG->Fill(all_chPulseArea_phd, duration);

    //*****************
    // TPC Low gain
    //*****************
    number_of_pulses+=evt->nPulses_TPCLG;
    h_nPulses_TPCLG->Fill(evt->nPulses_TPCLG);
    all_chPulseArea_phd=0;
    all_chPhotonCount_phd=0;
    for(int i=0; i<evt->nPulses_TPCLG; ++i){
      duration = 0;
      h_pulseArea_phd_TPCLG->Fill(evt->pulseArea_phd_TPCLG[i]);
      h_photonCount_TPCLG->Fill(evt->photonCount_TPCLG[i]);
      h_pulseStartTime_ns_TPCLG->Fill(evt->pulseStartTime_ns_TPCLG[i]);
      h_pulseEndTime_ns_TPCLG->Fill(evt->pulseEndTime_ns_TPCLG[i]);
      h_peakAmp_TPCLG->Fill(evt->peakAmp_TPCLG[i]);
      h_peakTime_ns_TPCLG->Fill(evt->peakTime_ns_TPCLG[i]);
      h_areaFractionTime5_ns_TPCLG->Fill(evt->areaFractionTime5_ns_TPCLG[i]);
      h_areaFractionTime25_ns_TPCLG->Fill(evt->areaFractionTime25_ns_TPCLG[i]);
      h_areaFractionTime50_ns_TPCLG->Fill(evt->areaFractionTime50_ns_TPCLG[i]);
      h_areaFractionTime75_ns_TPCLG->Fill(evt->areaFractionTime75_ns_TPCLG[i]);
      h_areaFractionTime95_ns_TPCLG->Fill(evt->areaFractionTime95_ns_TPCLG[i]);
      h_promptFraction50ns_TPCLG->Fill(evt->promptFraction50ns_TPCLG[i]);
      h_rmsWidth_ns_TPCLG->Fill(evt->rmsWidth_ns_TPCLG[i]);
      h_coincidence_TPCLG->Fill(evt->coincidence_TPCLG[i]);
      h_s1Probability_TPCLG->Fill(evt->s1Probability_TPCLG[i]);
      h_s2Probability_TPCLG->Fill(evt->s2Probability_TPCLG[i]);
      h_singlePEprobability_TPCLG->Fill(evt->singlePEprobability_TPCLG[i]);
      h_singleElectronProbability_TPCLG->Fill(evt->singlePEprobability_TPCLG[i]);
      h_otherProbability_TPCLG->Fill(evt->otherProbability_TPCLG[i]);
      h_otherS2Probability_TPCLG->Fill(evt->otherS2Probability_TPCLG[i]);

      for (unsigned int ii=0; ii<evt->chPulseArea_phd_TPCLG[i].size(); ii++) {
	h_channelID_TPCLG->Fill(evt->channelID_TPCLG[i][ii]-1000);
	// int pmtid =(evt->channelID_TPCLG[i][ii])-100;
	// hvec_chPulseArea_phd_TPCLG.at(pmtid)->Fill(evt->chPulseArea_phd_TPCLG[i][ii]);
	if(evt->chPulseArea_phd_TPCLG[i][ii]>3){
	  all_chPulseArea_phd+=evt->chPulseArea_phd_TPCLG[i][ii];
	  h_chPulseArea_phd_TPCLG->Fill(evt->chPulseArea_phd_TPCLG[i][ii]);
	  h_chPeakTime_ns_TPCLG->Fill(evt->chPeakTime_ns_TPCLG[i][ii]);
	}
      }//run over phd channels
      duration+=evt->pulseEndTime_ns_TPCLG[i]-evt->pulseStartTime_ns_TPCLG[i];
      h_pulse_area_width_TPCLG->Fill(evt->pulseArea_phd_TPCLG[i], evt->pulseEndTime_ns_TPCLG[i]-evt->pulseStartTime_ns_TPCLG[i]);
      h_pulse_start_area_TPCLG->Fill(evt->pulseStartTime_ns_TPCLG[i], evt->pulseArea_phd_TPCLG[i]);
    }//Pulses TPC
    h_pulseArea_phd_TPCLG_sum->Fill(all_chPulseArea_phd);
  }//int nevents
    
  //--------------------------------------------------
  // Convert to rate
  //--------------------------------------------------
  int count = 0;
  unsigned long long time = 0;
  unsigned long long time_day = 0;
  double rate = 0;
    //------------------------------------------------
  // end event loop
  //------------------------------------------------


  //write and close output file
  outfile->Write();
  outfile->Close();

  std::cout<<"times_ns: "<<min_time_ns<<" "<<max_time_ns<<" "<<processed_events<<endl;
  std::cout<<"summed time of all pulses "<<sum_pulse_time/1E9<<" sec"<<endl;
  std::cout<<"Overall time "<<last_s-first_s <<" sec, rate: "<<processed_events/(last_s-first_s)<<" evts/sec"<<", pulse rate "<<number_of_pulses/(last_s-first_s) <<endl;
  cout << "Done!"<<" "<<clock->RealTime()<<" s."<<endl;



  //check time
  evt->GetEntry(0);
  double begin_s=evt->triggerTimeStamp_s;
  evt->GetEntry(evt->chain->GetEntries()-1);
  double end_s=evt->triggerTimeStamp_s;
  std::cout<<" rate: "<<evt->chain->GetEntries()/(end_s-begin_s)<<" evt/sec"<<std::endl;


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


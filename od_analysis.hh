/**
 * # Macro to determine per pmt spe rate (hour and day)
 * Based on example macro from MDC2 kick-off
 * This macro reads all the files listed in the file `RQfile.list` and generate an output rootfile `speOD.root`.
 * ## Usage: 
 * Creation of the Rqlib, need to do only once (or as RQEvent format changes). 
 * ```bash
 * root -b -q configure.C
 * ```
 * Run the analysis macro:
 * ```bash
 *root -b -q load.C macro.C+ 
 * ```
 *
 *
 * The output file contains all the histograms ceated during the processing.
 * When you use the comman ```root -b -q load.C macro.C+```, root will automaticaly execute the function `macro`.
 *
 *---
 *# Description of the macro
 * 
 * In the following section, each part of the macro will be described. The major blocks are: 
 * 1. [Headers] (#Headers)
 * 2. [Event structure] (#event_struct)
 * 3. [Cut definition] (#cuts)
 * 4. [Main executable] (#main)
 * 5. [Load chain] (#load)
 * 6. [Event loop] (#loop)
 *
 */


/**
 * ## Headers <a name="Headers"></a>
 * Include all the headers
 *
 */
// Standard C++ headers
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

// STD namespace
using namespace std;

// Functions prototype
void load_chain(TString txtFileList, TChain* chain);
void eventLoop(TChain* events);
void init_histos();

/**
 * Data type: Background or AmLi
 */
TString option = "Bkg";

TFile* outfile;


/**
 * ## Lite event structure <a name="event_struct"></a>
 *
 * Description:
 *   Keep in this event structure only the required RQs for this analysis. 
 *
 * ### Methods:
 * - `GetEntry(const int n)`: Get the nth event of the chain
 */

class MyEvent {
public:
  MyEvent(){}
  virtual ~MyEvent(){}
  TChain* chain;
  Int_t currentTree;
  TBranch *b_eventID, *b_triggerTimeStamp_s, *b_triggerTimeStamp_ns;
  

  //MultipleScatter
  TBranch *b_nMultipleScatters;
  TBranch *b_multiplescatter_s1Area_phd;

  //SingleScatter
  TBranch *b_nSingleScatters, *b_singlescatter_odPromptArea, *b_singlescatter_nODPromptPulses, *b_singlescatter_odPromptPulseIDs, *b_singlescatter_odDelayedAreas, *b_singlescatter_nODDelayedPulses, *b_singlescatter_odDelayedPulseIDs,
    *b_singlescatter_energyER_keV, *b_singlescatter_energyNR_keV, *b_singlescatter_s1Area_phd, *b_singlescatter_s2Area_phd, *b_singlescatter_correctedS1Area_phd, *b_singlescatter_correctedS2Area_phd, *b_singlescatter_x_cm, *b_singlescatter_y_cm;

  //OD
  TBranch *b_nPulsesODLG, *b_singlePEprobability_ODLG, *b_multiplePEprobability_ODLG, *b_channelID_ODLG, *b_pulseStartTime_ns_ODLG, *b_pulseEndTime_ns_ODLG, 
    *b_pulseArea_phd_ODLG, *b_positiveArea_phd_ODLG, *b_negativeArea_phd_ODLG, *b_chPulseArea_phd_ODLG, *b_chPeakTime_ns_ODLG, *b_chPhotonCount_ODLG, *b_photonCount_ODLG, *b_peakAmp_ODLG, *b_peakTime_ns_ODLG,
    *b_areaFractionTime5_ns_ODLG, *b_areaFractionTime25_ns_ODLG, *b_areaFractionTime50_ns_ODLG, *b_areaFractionTime75_ns_ODLG, *b_areaFractionTime95_ns_ODLG, *b_promptFraction50ns_ODLG, *b_rmsWidth_ns_ODLG, *b_coincidence_ODLG,
    *b_s1Probability_ODLG, *b_s2Probability_ODLG, *b_singleElectronProbability_ODLG, *b_otherProbability_ODLG, *b_otherS2Probability_ODLG;

  //TPC LG
  TBranch *b_nPulsesTPCLG, *b_singlePEprobability_TPCLG, *b_multiplePEprobability_TPCLG, *b_channelID_TPCLG, *b_pulseStartTime_ns_TPCLG, *b_pulseEndTime_ns_TPCLG, 
    *b_pulseArea_phd_TPCLG, *b_positiveArea_phd_TPCLG, *b_negativeArea_phd_TPCLG, *b_chPulseArea_phd_TPCLG, *b_chPeakTime_ns_TPCLG, *b_chPhotonCount_TPCLG, *b_photonCount_TPCLG, *b_peakAmp_TPCLG, *b_peakTime_ns_TPCLG,
    *b_areaFractionTime5_ns_TPCLG, *b_areaFractionTime25_ns_TPCLG, *b_areaFractionTime50_ns_TPCLG, *b_areaFractionTime75_ns_TPCLG, *b_areaFractionTime95_ns_TPCLG, *b_promptFraction50ns_TPCLG, *b_rmsWidth_ns_TPCLG, *b_coincidence_TPCLG,
    *b_s1Probability_TPCLG, *b_s2Probability_TPCLG, *b_singleElectronProbability_TPCLG, *b_otherProbability_TPCLG, *b_otherS2Probability_TPCLG;

  //TPC HG
  TBranch *b_nPulsesTPCHG, *b_singlePEprobability_TPCHG, *b_multiplePEprobability_TPCHG, *b_channelID_TPCHG, *b_pulseStartTime_ns_TPCHG, *b_pulseEndTime_ns_TPCHG, 
    *b_pulseArea_phd_TPCHG, *b_positiveArea_phd_TPCHG, *b_negativeArea_phd_TPCHG, *b_chPulseArea_phd_TPCHG, *b_chPeakTime_ns_TPCHG, *b_chPhotonCount_TPCHG, *b_photonCount_TPCHG, *b_peakAmp_TPCHG, *b_peakTime_ns_TPCHG,
    *b_areaFractionTime5_ns_TPCHG, *b_areaFractionTime25_ns_TPCHG, *b_areaFractionTime50_ns_TPCHG, *b_areaFractionTime75_ns_TPCHG, *b_areaFractionTime95_ns_TPCHG, *b_promptFraction50ns_TPCHG, *b_rmsWidth_ns_TPCHG, *b_coincidence_TPCHG,
    *b_s1Probability_TPCHG, *b_s2Probability_TPCHG, *b_singleElectronProbability_TPCHG, *b_otherProbability_TPCHG, *b_otherS2Probability_TPCHG;

  //------------------------------------------------
  // Variable for extraction
  //------------------------------------------------

  //Multiple Scatter
  int nMultipleScatters;
  vector<float> multiplescatter_s1Area_phd;
  
  //Single Scatter
  int nSingleScatters;
  vector<float> singlescatter_odPromptArea;
  vector<int> singlescatter_nODPromptPulses;
  vector<vector<int>> singlescatter_odPromptPulseIDs, singlescatter_nODDelayedPulses; 
  //  vector<vector<float>> singlescatter_odDelayedAreas, singlescatter_odDelayedPulseIDs, singlescatter_energyER_keV, singlescatter_energyNR_keV, singlescatter_s1Area_phd, singlescatter_s2Area_phd, singlescatter_correctedS1Area_phd, singlescatter_correctedS2Area_phd;
  vector<float> singlescatter_odDelayedAreas, singlescatter_odDelayedPulseIDs, singlescatter_energyER_keV, singlescatter_energyNR_keV, singlescatter_s1Area_phd, singlescatter_s2Area_phd, singlescatter_correctedS1Area_phd, singlescatter_correctedS2Area_phd, singlescatter_x_cm, singlescatter_y_cm;

  //Event Header
  int eventID, nPulses_ODLG;
  unsigned long triggerTimeStamp_ns, triggerTimeStamp_s;
  
  //OD
  vector<int> pulseStartTime_ns_ODLG, pulseEndTime_ns_ODLG, photonCount_ODLG;
  vector<float> pulseArea_phd_ODLG, positiveArea_phd_ODLG, negativeArea_phd_ODLG, peakAmp_ODLG, promptFraction50ns_ODLG, s1Probability_ODLG, s2Probability_ODLG, 
    singleElectronProbability_ODLG, otherProbability_ODLG, otherS2Probability_ODLG, singlePEprobability_ODLG, multiplePEprobability_ODLG;
  vector<int>  peakTime_ns_ODLG, areaFractionTime5_ns_ODLG, areaFractionTime25_ns_ODLG, areaFractionTime50_ns_ODLG, areaFractionTime75_ns_ODLG, areaFractionTime95_ns_ODLG,
    rmsWidth_ns_ODLG, coincidence_ODLG;
  vector<vector<int>>   channelID_ODLG, chPeakTime_ns_ODLG, chPhotonCount_ODLG;
  vector<vector<float>> chPulseArea_phd_ODLG;

  //TPC LG
  int nPulses_TPCLG;
  vector<int> pulseStartTime_ns_TPCLG, pulseEndTime_ns_TPCLG, photonCount_TPCLG;
  vector<float> pulseArea_phd_TPCLG, positiveArea_phd_TPCLG, negativeArea_phd_TPCLG, peakAmp_TPCLG, promptFraction50ns_TPCLG, s1Probability_TPCLG, s2Probability_TPCLG, 
    singleElectronProbability_TPCLG, otherProbability_TPCLG, otherS2Probability_TPCLG, singlePEprobability_TPCLG, multiplePEprobability_TPCLG;
  vector<int>  peakTime_ns_TPCLG, areaFractionTime5_ns_TPCLG, areaFractionTime25_ns_TPCLG, areaFractionTime50_ns_TPCLG, areaFractionTime75_ns_TPCLG, areaFractionTime95_ns_TPCLG,
    rmsWidth_ns_TPCLG, coincidence_TPCLG;
  vector<vector<int>>   channelID_TPCLG, chPeakTime_ns_TPCLG, chPhotonCount_TPCLG;
  vector<vector<float>> chPulseArea_phd_TPCLG;

  //TPC HG
  int nPulses_TPCHG;
  vector<int> pulseStartTime_ns_TPCHG, pulseEndTime_ns_TPCHG, photonCount_TPCHG;
  vector<float> pulseArea_phd_TPCHG, positiveArea_phd_TPCHG, negativeArea_phd_TPCHG, peakAmp_TPCHG, promptFraction50ns_TPCHG, s1Probability_TPCHG, s2Probability_TPCHG, 
    singleElectronProbability_TPCHG, otherProbability_TPCHG, otherS2Probability_TPCHG, singlePEprobability_TPCHG, multiplePEprobability_TPCHG;
  vector<int>  peakTime_ns_TPCHG, areaFractionTime5_ns_TPCHG, areaFractionTime25_ns_TPCHG, areaFractionTime50_ns_TPCHG, areaFractionTime75_ns_TPCHG, areaFractionTime95_ns_TPCHG,
    rmsWidth_ns_TPCHG, coincidence_TPCHG;
  vector<vector<int>>   channelID_TPCHG, chPeakTime_ns_TPCHG, chPhotonCount_TPCHG;
  vector<vector<float>> chPulseArea_phd_TPCHG;


  // //define histos
  // TH1D* h_triggerTimeStamp_ns, h_triggerTimeStamp_s, h_triggerTimeStamp_s_ODLG;
  // TH1F* h_pulse_duration_ms_ODLG, h_pulseArea_phd_ODLG, h_pulseArea_phd_ODLG_sum, h_positiveArea_phd_ODLG, h_negativeArea_phd_ODLG, h_nPulses_ODLG, h_chPulseArea_phd_ODLG, h_chPeakTime_ns_ODLG, h_channelID_ODLG, h_photonCount_ODLG,
  // h_pulseStartTime_ns_ODLG, h_pulseEndTime_ns_ODLG, h_peakAmp_ODLG, h_peakTime_ns_ODLG, h_areaFractionTime5_ns_ODLG, h_areaFractionTime25_ns_ODLG, h_areaFractionTime50_ns_ODLG, h_areaFractionTime75_ns_ODLG, h_areaFractionTime95_ns_ODLG, h_promptFraction50ns_ODLG, h_rmsWidth_ns_ODLG, h_coincidence_ODLG, h_s1Probability_ODLG, h_s2Probability_ODLG, h_singleElectronProbability_ODLG, h_otherProbability_ODLG;

  // TH1F* h_pulse_duration_ms_OD_lt, h_pulseArea_phd_OD_lt, h_pulseArea_phd_OD_lt_sum, h_positiveArea_phd_OD_lt, h_negativeArea_phd_OD_lt, h_nPulses_OD_lt, h_chPulseArea_phd_OD_lt, h_chPeakTime_ns_OD_lt, h_channelID_OD_lt, h_photonCount_OD_lt,  h_pulseStartTime_ns_OD_lt, h_pulseEndTime_ns_OD_lt, h_peakAmp_OD_lt, h_peakTime_ns_OD_lt, h_areaFractionTime5_ns_OD_lt, h_areaFractionTime25_ns_OD_lt, h_areaFractionTime50_ns_OD_lt, h_areaFractionTime75_ns_OD_lt, h_areaFractionTime95_ns_OD_lt, h_promptFraction50ns_OD_lt, h_rmsWidth_ns_OD_lt, h_coincidence_OD_lt, h_s1Probability_OD_lt, h_s2Probability_OD_lt, h_singleElectronProbability_OD_lt, h_otherProbability_OD_lt;

  // TH1F* h_pulse_duration_ms_OD_peak, h_pulseArea_phd_OD_peak, h_pulseArea_phd_OD_peak_sum, h_positiveArea_phd_OD_peak, h_negativeArea_phd_OD_peak, h_nPulses_OD_peak, h_chPulseArea_phd_OD_peak, h_chPeakTime_ns_OD_peak, h_channelID_OD_peakg, h_photonCount_OD_peak,  h_pulseStartTime_ns_OD_peak, h_pulseEndTime_ns_OD_peak, h_peakAmp_OD_peak, h_peakTime_ns_OD_peak, h_areaFractionTime5_ns_OD_peak, h_areaFractionTime25_ns_OD_peak, h_areaFractionTime50_ns_OD_peak, h_areaFractionTime75_ns_OD_peak,eh_areaFractionTime95_ns_OD_peak, h_promptFraction50ns_OD_peak, h_rmsWidth_ns_OD_peak, h_coincidence_OD_peak, h_s1Probability_OD_peak, h_s2Probability_OD_peak, h_singleElectronProbability_OD_peak, h_otherProbability_OD_ht;


  void LoadBranches(TChain* tree);
  void GetEntry(const int n);
  void DefHistos();
};


void MyEvent::LoadBranches(TChain* tree)
{
  chain = tree;
  chain->SetMakeClass(1);

  chain->SetBranchAddress("eventHeader.eventID",                  &eventID,                          &b_eventID);
  chain->SetBranchAddress("eventHeader.triggerTimeStamp_s",       &triggerTimeStamp_s,               &b_triggerTimeStamp_s);

  chain->SetBranchAddress("multipleScatters.nMultipleScatters",   &nMultipleScatters,                &b_nMultipleScatters);
  chain->SetBranchAddress("multipleScatters.s1Area_phd",          &multiplescatter_s1Area_phd,       &b_multiplescatter_s1Area_phd);

  //Single Scatter
  chain->SetBranchAddress("singleScatters.nSingleScatters",       &nSingleScatters,                  &b_nSingleScatters);
  chain->SetBranchAddress("singleScatters.odPromptArea",          &singlescatter_odPromptArea,                     &b_singlescatter_odPromptArea);
  chain->SetBranchAddress("singleScatters.nODPromptPulses",       &singlescatter_nODPromptPulses,                  &b_singlescatter_nODPromptPulses);
  chain->SetBranchAddress("singleScatters.odPromptPulseIDs",      &singlescatter_odPromptPulseIDs,                 &b_singlescatter_odPromptPulseIDs);
  chain->SetBranchAddress("singleScatters.odDelayedAreas",        &singlescatter_odDelayedAreas,                   &b_singlescatter_odDelayedAreas);
  chain->SetBranchAddress("singleScatters.nODDelayedPulses",      &singlescatter_nODDelayedPulses,                 &b_singlescatter_nODDelayedPulses);
  chain->SetBranchAddress("singleScatters.odDelayedPulseIDs",     &singlescatter_odDelayedPulseIDs,                &b_singlescatter_odDelayedPulseIDs);
  chain->SetBranchAddress("singleScatters.energyER_keV",          &singlescatter_energyER_keV,                     &b_singlescatter_energyER_keV);
  chain->SetBranchAddress("singleScatters.energyNR_keV",          &singlescatter_energyNR_keV,                     &b_singlescatter_energyNR_keV);
  chain->SetBranchAddress("singleScatters.s1Area_phd",            &singlescatter_s1Area_phd,                       &b_singlescatter_s1Area_phd);
  chain->SetBranchAddress("singleScatters.s2Area_phd",            &singlescatter_s2Area_phd,                       &b_singlescatter_s2Area_phd);
  chain->SetBranchAddress("singleScatters.correctedS1Area_phd",   &singlescatter_correctedS1Area_phd,              &b_singlescatter_correctedS1Area_phd);
  chain->SetBranchAddress("singleScatters.correctedS2Area_phd",   &singlescatter_correctedS2Area_phd,              &b_singlescatter_correctedS2Area_phd);
  chain->SetBranchAddress("singleScatters.x_cm",                  &singlescatter_x_cm,                             &b_singlescatter_x_cm);
  chain->SetBranchAddress("singleScatters.y_cm",                  &singlescatter_y_cm,                             &b_singlescatter_y_cm);
  

  //OD
  chain->SetBranchAddress("eventHeader.triggerTimeStamp_ns",      &triggerTimeStamp_ns,              &b_triggerTimeStamp_ns);
  chain->SetBranchAddress("pulsesODLG.nPulses",                   &nPulses_ODLG,                     &b_nPulsesODLG);
  chain->SetBranchAddress("pulsesODLG.singlePEprobability",       &singlePEprobability_ODLG,         &b_singlePEprobability_ODLG);
  //  chain->SetBranchAddress("pulsesODLG.multiplePEprobability",     &multiplePEprobability_ODLG,       &b_multiplePEprobability_ODLG);
  chain->SetBranchAddress("pulsesODLG.chID",                      &channelID_ODLG,                   &b_channelID_ODLG);
  chain->SetBranchAddress("pulsesODLG.pulseStartTime_ns",         &pulseStartTime_ns_ODLG,           &b_pulseStartTime_ns_ODLG);
  chain->SetBranchAddress("pulsesODLG.pulseEndTime_ns",           &pulseEndTime_ns_ODLG,             &b_pulseEndTime_ns_ODLG);
  chain->SetBranchAddress("pulsesODLG.pulseArea_phd",             &pulseArea_phd_ODLG,               &b_pulseArea_phd_ODLG);
  chain->SetBranchAddress("pulsesODLG.positiveArea_phd",          &positiveArea_phd_ODLG,            &b_positiveArea_phd_ODLG);
  chain->SetBranchAddress("pulsesODLG.negativeArea_phd",          &negativeArea_phd_ODLG,            &b_negativeArea_phd_ODLG);
  chain->SetBranchAddress("pulsesODLG.chPulseArea_phd",           &chPulseArea_phd_ODLG,             &b_chPulseArea_phd_ODLG);
  chain->SetBranchAddress("pulsesODLG.chPeakTime_ns",             &chPeakTime_ns_ODLG,               &b_chPeakTime_ns_ODLG);
  chain->SetBranchAddress("pulsesODLG.chPhotonCount",             &chPhotonCount_ODLG,               &b_chPhotonCount_ODLG);
  chain->SetBranchAddress("pulsesODLG.photonCount",               &photonCount_ODLG,                 &b_photonCount_ODLG);
  chain->SetBranchAddress("pulsesODLG.peakAmp",                   &peakAmp_ODLG,                     &b_peakAmp_ODLG);
  chain->SetBranchAddress("pulsesODLG.peakTime_ns",               &peakTime_ns_ODLG,                 &b_peakTime_ns_ODLG);
  chain->SetBranchAddress("pulsesODLG.areaFractionTime5_ns",      &areaFractionTime5_ns_ODLG,        &b_areaFractionTime5_ns_ODLG);
  chain->SetBranchAddress("pulsesODLG.areaFractionTime25_ns",     &areaFractionTime25_ns_ODLG,       &b_areaFractionTime25_ns_ODLG);
  chain->SetBranchAddress("pulsesODLG.areaFractionTime50_ns",     &areaFractionTime50_ns_ODLG,       &b_areaFractionTime50_ns_ODLG);
  chain->SetBranchAddress("pulsesODLG.areaFractionTime75_ns",     &areaFractionTime75_ns_ODLG,       &b_areaFractionTime75_ns_ODLG);
  chain->SetBranchAddress("pulsesODLG.areaFractionTime95_ns",     &areaFractionTime95_ns_ODLG,       &b_areaFractionTime95_ns_ODLG);
  chain->SetBranchAddress("pulsesODLG.promptFraction50ns",        &promptFraction50ns_ODLG,          &b_promptFraction50ns_ODLG);
  chain->SetBranchAddress("pulsesODLG.rmsWidth_ns",               &rmsWidth_ns_ODLG,                 &b_rmsWidth_ns_ODLG);
  chain->SetBranchAddress("pulsesODLG.coincidence",               &coincidence_ODLG,                 &b_coincidence_ODLG);
  chain->SetBranchAddress("pulsesODLG.s1Probability",             &s1Probability_ODLG,               &b_s1Probability_ODLG);
  chain->SetBranchAddress("pulsesODLG.s2Probability",             &s2Probability_ODLG,               &b_s2Probability_ODLG);
  chain->SetBranchAddress("pulsesODLG.singleElectronProbability", &singleElectronProbability_ODLG,   &b_singleElectronProbability_ODLG);
  chain->SetBranchAddress("pulsesODLG.otherProbability",          &otherProbability_ODLG,            &b_otherProbability_ODLG);
  chain->SetBranchAddress("pulsesODLG.otherS2Probability",        &otherS2Probability_ODLG,          &b_otherS2Probability_ODLG);

  //TPC LG
  chain->SetBranchAddress("pulsesTPCLG.nPulses",                   &nPulses_TPCLG,                     &b_nPulsesTPCLG);
  chain->SetBranchAddress("pulsesTPCLG.singlePEprobability",       &singlePEprobability_TPCLG,         &b_singlePEprobability_TPCLG);
  //  chain->SetBranchAddress("pulsesTPCLG.multiplePEprobability",     &multiplePEprobability_TPCLG,       &b_multiplePEprobability_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.chID",                      &channelID_TPCLG,                   &b_channelID_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.pulseStartTime_ns",         &pulseStartTime_ns_TPCLG,           &b_pulseStartTime_ns_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.pulseEndTime_ns",           &pulseEndTime_ns_TPCLG,             &b_pulseEndTime_ns_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.pulseArea_phd",             &pulseArea_phd_TPCLG,               &b_pulseArea_phd_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.positiveArea_phd",          &positiveArea_phd_TPCLG,            &b_positiveArea_phd_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.negativeArea_phd",          &negativeArea_phd_TPCLG,            &b_negativeArea_phd_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.chPulseArea_phd",           &chPulseArea_phd_TPCLG,             &b_chPulseArea_phd_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.chPeakTime_ns",             &chPeakTime_ns_TPCLG,               &b_chPeakTime_ns_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.chPhotonCount",             &chPhotonCount_TPCLG,               &b_chPhotonCount_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.photonCount",               &photonCount_TPCLG,                 &b_photonCount_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.peakAmp",                   &peakAmp_TPCLG,                     &b_peakAmp_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.peakTime_ns",               &peakTime_ns_TPCLG,                 &b_peakTime_ns_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.areaFractionTime5_ns",      &areaFractionTime5_ns_TPCLG,        &b_areaFractionTime5_ns_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.areaFractionTime25_ns",     &areaFractionTime25_ns_TPCLG,       &b_areaFractionTime25_ns_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.areaFractionTime50_ns",     &areaFractionTime50_ns_TPCLG,       &b_areaFractionTime50_ns_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.areaFractionTime75_ns",     &areaFractionTime75_ns_TPCLG,       &b_areaFractionTime75_ns_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.areaFractionTime95_ns",     &areaFractionTime95_ns_TPCLG,       &b_areaFractionTime95_ns_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.promptFraction50ns",        &promptFraction50ns_TPCLG,          &b_promptFraction50ns_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.rmsWidth_ns",               &rmsWidth_ns_TPCLG,                 &b_rmsWidth_ns_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.coincidence",               &coincidence_TPCLG,                 &b_coincidence_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.s1Probability",             &s1Probability_TPCLG,               &b_s1Probability_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.s2Probability",             &s2Probability_TPCLG,               &b_s2Probability_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.singleElectronProbability", &singleElectronProbability_TPCLG,   &b_singleElectronProbability_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.otherProbability",          &otherProbability_TPCLG,            &b_otherProbability_TPCLG);
  chain->SetBranchAddress("pulsesTPCLG.otherS2Probability",        &otherS2Probability_TPCLG,          &b_otherS2Probability_TPCLG);

  //TPC HG
  chain->SetBranchAddress("pulsesTPCHG.nPulses",                   &nPulses_TPCHG,                     &b_nPulsesTPCHG);
  chain->SetBranchAddress("pulsesTPCHG.singlePEprobability",       &singlePEprobability_TPCHG,         &b_singlePEprobability_TPCHG);
  //  chain->SetBranchAddress("pulsesTPCHG.multiplePEprobability",     &multiplePEprobability_TPCHG,       &b_multiplePEprobability_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.chID",                      &channelID_TPCHG,                   &b_channelID_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.pulseStartTime_ns",         &pulseStartTime_ns_TPCHG,           &b_pulseStartTime_ns_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.pulseEndTime_ns",           &pulseEndTime_ns_TPCHG,             &b_pulseEndTime_ns_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.pulseArea_phd",             &pulseArea_phd_TPCHG,               &b_pulseArea_phd_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.positiveArea_phd",          &positiveArea_phd_TPCHG,            &b_positiveArea_phd_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.negativeArea_phd",          &negativeArea_phd_TPCHG,            &b_negativeArea_phd_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.chPulseArea_phd",           &chPulseArea_phd_TPCHG,             &b_chPulseArea_phd_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.chPeakTime_ns",             &chPeakTime_ns_TPCHG,               &b_chPeakTime_ns_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.chPhotonCount",             &chPhotonCount_TPCHG,               &b_chPhotonCount_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.photonCount",               &photonCount_TPCHG,                 &b_photonCount_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.peakAmp",                   &peakAmp_TPCHG,                     &b_peakAmp_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.peakTime_ns",               &peakTime_ns_TPCHG,                 &b_peakTime_ns_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.areaFractionTime5_ns",      &areaFractionTime5_ns_TPCHG,        &b_areaFractionTime5_ns_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.areaFractionTime25_ns",     &areaFractionTime25_ns_TPCHG,       &b_areaFractionTime25_ns_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.areaFractionTime50_ns",     &areaFractionTime50_ns_TPCHG,       &b_areaFractionTime50_ns_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.areaFractionTime75_ns",     &areaFractionTime75_ns_TPCHG,       &b_areaFractionTime75_ns_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.areaFractionTime95_ns",     &areaFractionTime95_ns_TPCHG,       &b_areaFractionTime95_ns_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.promptFraction50ns",        &promptFraction50ns_TPCHG,          &b_promptFraction50ns_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.rmsWidth_ns",               &rmsWidth_ns_TPCHG,                 &b_rmsWidth_ns_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.coincidence",               &coincidence_TPCHG,                 &b_coincidence_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.s1Probability",             &s1Probability_TPCHG,               &b_s1Probability_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.s2Probability",             &s2Probability_TPCHG,               &b_s2Probability_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.singleElectronProbability", &singleElectronProbability_TPCHG,   &b_singleElectronProbability_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.otherProbability",          &otherProbability_TPCHG,            &b_otherProbability_TPCHG);
  chain->SetBranchAddress("pulsesTPCHG.otherS2Probability",        &otherS2Probability_TPCHG,          &b_otherS2Probability_TPCHG);


}


void MyEvent::GetEntry(const int n)
{
  Long64_t ientry = chain->LoadTree(n);
  if (chain->GetTreeNumber() != currentTree) {
    currentTree = chain->GetTreeNumber();
  }
  b_eventID                           ->GetEntry(ientry);
  b_triggerTimeStamp_s                ->GetEntry(ientry);
  b_triggerTimeStamp_ns               ->GetEntry(ientry);

  //Mult Scatter
  b_nMultipleScatters                 ->GetEntry(ientry);
  b_multiplescatter_s1Area_phd        ->GetEntry(ientry);

  //SingleScatter
  b_nSingleScatters                   ->GetEntry(ientry);
  b_singlescatter_odPromptArea                      ->GetEntry(ientry);
  b_singlescatter_nODPromptPulses                   ->GetEntry(ientry);
  b_singlescatter_odPromptPulseIDs                  ->GetEntry(ientry);
  b_singlescatter_odDelayedAreas                    ->GetEntry(ientry);
  b_singlescatter_nODDelayedPulses                  ->GetEntry(ientry);
  b_singlescatter_odDelayedPulseIDs                 ->GetEntry(ientry);
  b_singlescatter_energyER_keV                      ->GetEntry(ientry);
  b_singlescatter_energyNR_keV                      ->GetEntry(ientry);
  b_singlescatter_s1Area_phd                        ->GetEntry(ientry);
  b_singlescatter_s2Area_phd                        ->GetEntry(ientry);
  b_singlescatter_correctedS1Area_phd               ->GetEntry(ientry);
  b_singlescatter_correctedS2Area_phd               ->GetEntry(ientry);
  b_singlescatter_x_cm                              ->GetEntry(ientry);
  b_singlescatter_y_cm                              ->GetEntry(ientry);



  //OD
  b_nPulsesODLG                      ->GetEntry(ientry);
  b_singlePEprobability_ODLG          ->GetEntry(ientry);
  //  b_multiplePEprobability_ODLG        ->GetEntry(ientry);
  b_channelID_ODLG                    ->GetEntry(ientry);
  b_pulseStartTime_ns_ODLG            ->GetEntry(ientry);
  b_pulseEndTime_ns_ODLG              ->GetEntry(ientry);
  b_pulseArea_phd_ODLG                ->GetEntry(ientry);
  b_positiveArea_phd_ODLG             ->GetEntry(ientry);
  b_negativeArea_phd_ODLG             ->GetEntry(ientry);
  b_chPulseArea_phd_ODLG              ->GetEntry(ientry);
  b_chPeakTime_ns_ODLG                ->GetEntry(ientry);
  b_chPhotonCount_ODLG                ->GetEntry(ientry);
  b_photonCount_ODLG                  ->GetEntry(ientry);
  b_peakAmp_ODLG                      ->GetEntry(ientry);
  b_peakTime_ns_ODLG                  ->GetEntry(ientry);
  b_areaFractionTime5_ns_ODLG         ->GetEntry(ientry);
  b_areaFractionTime25_ns_ODLG        ->GetEntry(ientry);
  b_areaFractionTime50_ns_ODLG        ->GetEntry(ientry);
  b_areaFractionTime75_ns_ODLG        ->GetEntry(ientry);
  b_areaFractionTime95_ns_ODLG        ->GetEntry(ientry);
  b_promptFraction50ns_ODLG           ->GetEntry(ientry);
  b_rmsWidth_ns_ODLG                  ->GetEntry(ientry);
  b_coincidence_ODLG                  ->GetEntry(ientry);
  b_s1Probability_ODLG                ->GetEntry(ientry); 
  b_s2Probability_ODLG                 ->GetEntry(ientry);
  b_singleElectronProbability_ODLG    ->GetEntry(ientry);
  b_otherProbability_ODLG             ->GetEntry(ientry);
  b_otherS2Probability_ODLG           ->GetEntry(ientry);

  //TPC  LG
  b_nPulsesTPCLG                        ->GetEntry(ientry);
  b_singlePEprobability_TPCLG          ->GetEntry(ientry);
  //  b_multiplePEprobability_TPCLG        ->GetEntry(ientry);
  b_channelID_TPCLG                    ->GetEntry(ientry);
  b_pulseStartTime_ns_TPCLG            ->GetEntry(ientry);
  b_pulseEndTime_ns_TPCLG              ->GetEntry(ientry);
  b_pulseArea_phd_TPCLG                ->GetEntry(ientry);
  b_positiveArea_phd_TPCLG             ->GetEntry(ientry);
  b_negativeArea_phd_TPCLG             ->GetEntry(ientry);
  b_chPulseArea_phd_TPCLG              ->GetEntry(ientry);
  b_chPeakTime_ns_TPCLG                ->GetEntry(ientry);
  b_chPhotonCount_TPCLG                ->GetEntry(ientry);
  b_photonCount_TPCLG                  ->GetEntry(ientry);
  b_peakAmp_TPCLG                      ->GetEntry(ientry);
  b_peakTime_ns_TPCLG                  ->GetEntry(ientry);
  b_areaFractionTime5_ns_TPCLG         ->GetEntry(ientry);
  b_areaFractionTime25_ns_TPCLG        ->GetEntry(ientry);
  b_areaFractionTime50_ns_TPCLG        ->GetEntry(ientry);
  b_areaFractionTime75_ns_TPCLG        ->GetEntry(ientry);
  b_areaFractionTime95_ns_TPCLG        ->GetEntry(ientry);
  b_promptFraction50ns_TPCLG           ->GetEntry(ientry);
  b_rmsWidth_ns_TPCLG                  ->GetEntry(ientry);
  b_coincidence_TPCLG                  ->GetEntry(ientry);
  b_s1Probability_TPCLG                ->GetEntry(ientry); 
  b_s2Probability_TPCLG                 ->GetEntry(ientry);
  b_singleElectronProbability_TPCLG    ->GetEntry(ientry);
  b_otherProbability_TPCLG             ->GetEntry(ientry);
  b_otherS2Probability_TPCLG           ->GetEntry(ientry);

  //TPC  HG
  b_nPulsesTPCHG                        ->GetEntry(ientry);
  b_singlePEprobability_TPCHG          ->GetEntry(ientry);
  //  b_multiplePEprobability_TPCHG        ->GetEntry(ientry);
  b_channelID_TPCHG                    ->GetEntry(ientry);
  b_pulseStartTime_ns_TPCHG            ->GetEntry(ientry);
  b_pulseEndTime_ns_TPCHG              ->GetEntry(ientry);
  b_pulseArea_phd_TPCHG                ->GetEntry(ientry);
  b_positiveArea_phd_TPCHG             ->GetEntry(ientry);
  b_negativeArea_phd_TPCHG             ->GetEntry(ientry);
  b_chPulseArea_phd_TPCHG              ->GetEntry(ientry);
  b_chPeakTime_ns_TPCHG                ->GetEntry(ientry);
  b_chPhotonCount_TPCHG                ->GetEntry(ientry);
  b_photonCount_TPCHG                  ->GetEntry(ientry);
  b_peakAmp_TPCHG                      ->GetEntry(ientry);
  b_peakTime_ns_TPCHG                  ->GetEntry(ientry);
  b_areaFractionTime5_ns_TPCHG         ->GetEntry(ientry);
  b_areaFractionTime25_ns_TPCHG        ->GetEntry(ientry);
  b_areaFractionTime50_ns_TPCHG        ->GetEntry(ientry);
  b_areaFractionTime75_ns_TPCHG        ->GetEntry(ientry);
  b_areaFractionTime95_ns_TPCHG        ->GetEntry(ientry);
  b_promptFraction50ns_TPCHG           ->GetEntry(ientry);
  b_rmsWidth_ns_TPCHG                  ->GetEntry(ientry);
  b_coincidence_TPCHG                  ->GetEntry(ientry);
  b_s1Probability_TPCHG                ->GetEntry(ientry); 
  b_s2Probability_TPCHG                 ->GetEntry(ientry);
  b_singleElectronProbability_TPCHG    ->GetEntry(ientry);
  b_otherProbability_TPCHG             ->GetEntry(ientry);
  b_otherS2Probability_TPCHG           ->GetEntry(ientry);
}


// void MyEvent::DefHistos(){
//   cout<<"bla"<<endl;
//   *h_triggerTimeStamp_ns = new TH1D("triggerTimeStamp_ns", "triggerTimeStamp_ns", 100000, 0, 2.5e+08);
// }


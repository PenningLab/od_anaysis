// h_pulseArea_phd_ODLG_lt
// h_photonCount_ODLG_lt
// h_pulseStartTime_ns_ODLG_lt
// h_pulseEndTime_ns_ODLG_lt
// h_peakAmp_ODLG_lt
// h_peakTime_ns_ODLG_lt
// h_areaFractionTime5_ns_ODLG_lt
// h_areaFractionTime25_ns_ODLG_lt
// h_areaFractionTime50_ns_ODLG_lt
// h_areaFractionTime75_ns_ODLG_lt
// h_areaFractionTime95_ns_ODLG_lt
// h_promptFraction50ns_ODLG_lt
// h_rmsWidth_ns_ODLG_lt
// h_coincidence_ODLG_lt


void compare_plots(TString histname, TFile* _file0){
  cout<<histname<<endl;
  gStyle->SetOptStat(000000);  
  TH1F* h_OD = (TH1F*) _file0->Get(histname);
  h_OD->Print();
  histname.ReplaceAll("OD", "TPC");
  TH1F* h_TPC = (TH1F*) _file0->Get(histname);

  TCanvas* c2=new TCanvas(histname);
  c2->SetLogy();

  h_OD->SetLineColor(46); h_OD->SetLineWidth(2); h_OD->Scale(1/h_OD->Integral()); 
  h_TPC->SetLineColor(38); h_TPC->SetLineWidth(2); h_TPC->Scale(1/h_TPC->Integral()); 

  h_OD->SetTitle("");
  h_OD->SetXTitle(histname);
  h_OD->Draw("histo");
  h_TPC->Draw("histosame");
  
  
  TLegend* leg = new TLegend(0.7,0.7,0.89,0.89);
  leg->SetBorderSize(0);
  leg->AddEntry(h_OD,"OD LG","l");
  leg->AddEntry(h_TPC,"TPC LG","l");
  leg->Draw();
  c2->SaveAs(histname+"_comp.png");

}

void plot2D(TString histname, TString xtitle, TString ytitle, TFile* _file0){
  cout<<histname<<endl;
  gStyle->SetOptStat(000000);  
  TH2F* h2 = (TH2F*) _file0->Get(histname);
  TCanvas* c1=new TCanvas(histname);
  h2->SetTitle("");
  h2->SetXTitle(xtitle);
  h2->SetYTitle(ytitle);
  h2->Draw("COLZ");
  c1->SaveAs(histname+".png");

}

void plot(){

  TFile* file0=new TFile("AmLi.root", "OPEN");
   plot2D("pulse_area_width_ODLG", "pulse area [phd]", "pulse width [ns]", file0);
   plot2D("pulse_area_width_TPCLG", "pulse area [phd]", "pulse width [ns]", file0); 
   plot2D("pulse_start_area_ODLG", "start time [ns]", "pulse area [phd]", file0); 
   plot2D("pulse_start_area_TPCLG", "start time [ns]", "pulse area [phd]", file0); 

//    compare_plots("triggerTimeStamp_s", file0);
//    compare_plots("pulse_duration_ms", file0);
  compare_plots("pulseArea_phd_ODLG", file0);
  compare_plots("pulseArea_phd_ODLG_sum", file0);
  compare_plots("positiveArea_phd_ODLG", file0);
  compare_plots("negativeArea_phd_ODLG", file0);
  compare_plots("nPulses_ODLG", file0);
  compare_plots("chPulseArea_phd_ODLG", file0);
  compare_plots("chPeakTime_ns_ODLG", file0);
  compare_plots("channelID_ODLG", file0);
  compare_plots("photonCount_ODLG", file0);
  compare_plots("pulseStartTime_ns_ODLG", file0);
  compare_plots("pulseEndTime_ns_ODLG", file0);
  compare_plots("peakAmp_ODLG", file0);
  compare_plots("peakTime_ns_ODLG", file0);
  compare_plots("areaFractionTime5_ns_ODLG", file0);
  compare_plots("areaFractionTime25_ns_ODLG", file0);
  compare_plots("areaFractionTime50_ns_ODLG", file0);
  compare_plots("areaFractionTime75_ns_ODLG", file0);
  compare_plots("areaFractionTime95_ns_ODLG", file0);
  compare_plots("promptFraction50ns_ODLG", file0);
  compare_plots("rmsWidth_ns_ODLG", file0);
  compare_plots("coincidence_ODLG", file0);
  compare_plots("s1Probability_ODLG", file0);
  compare_plots("s2Probability_ODLG", file0);
  compare_plots("singleElectronProbability_ODLG", file0);
  compare_plots("singlePEprobability_ODLG", file0);
  compare_plots("otherProbability_ODLG", file0);
  compare_plots("otherS2Probability_ODLG", file0);

//   TH1F* h_photonCount = (TH1F*) _file0->Get("photonCount_ODLG");
//   TH1F* h_pulseArea_phd_ODLG = (TH1F*) _file0->Get("pulseArea_phd_ODLG");
//   TCanvas* c2=new TCanvas("photonCount");
//   pulseArea_phd_ODLG->SetLineColor(2);
//   c2->SetLogy();
//   h_photonCount->SetLineWidth(2);
//   h_pulseArea_phd_ODLG->SetLineWidth(2);

//   h_photonCount->Rebin(2);
//   h_pulseArea_phd_ODLG->Rebin(2);

//   h_photonCount->SetTitle("");
//   h_photonCount->SetXTitle("photon count");
//   h_photonCount->Draw("histo");
//   h_pulseArea_phd_ODLG->Draw("same");
  
//   cout<<"bla "<<h_pulseArea_phd_ODLG->Integral()/h_photonCount->Integral()<<endl;
  
//   leg = new TLegend(0.7,0.7,0.89,0.89);
//   leg->SetBorderSize(0);
//   leg->AddEntry(h_photonCount,"photonCount","l");
//   leg->AddEntry(h_pulseArea_phd_ODLG,"pulseArea","l");
//   leg->Draw();
//   c2->SaveAs("h_photonCount_ODLG.pdf");

//   compare_plots("pulseArea_phd_ODLG", file0);
//   compare_plots("photonCount_ODLG", file0);
//   compare_plots("pulseStartTime_ns_ODLG", file0);
//   compare_plots("pulseEndTime_ns_ODLG", file0);
//   compare_plots("peakAmp_ODLG", file0);
//   compare_plots("peakTime_ns_ODLG", file0);
//   compare_plots("areaFractionTime5_ns_ODLG", file0);
//   compare_plots("areaFractionTime25_ns_ODLG", file0);
//   compare_plots("areaFractionTime50_ns_ODLG", file0);
//   compare_plots("areaFractionTime75_ns_ODLG", file0);
//   compare_plots("areaFractionTime95_ns_ODLG", file0);
//   compare_plots("promptFraction50ns_ODLG", file0);
//   compare_plots("rmsWidth_ns_ODLG", file0);
//   compare_plots("coincidence_ODLG", file0);

}

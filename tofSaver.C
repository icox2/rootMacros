{
  double banGate(double x){return 200+4*(0.95*939.6*107*107/(2*x*x)-8.0*(1-exp(-1*0.1*pow(939.6*107*107/(2*x*x),0.9))));}
  TString fname = (TString)_file0->GetName();
  Int_t st = fname.First("_");
  Int_t end = fname.Length();
  fname.Remove(st,end-st);
  cout<<fname<<endl;
  //OutputTree->Draw("vandle_qdc:vandle_corTof>>qc(2000,-30,970,1000,0,20000)","dT>0 && vandle_qdc<banGate(vandle_corTof) && cloverAAB_E>1150 && cloverAAB_E<1165","colz");
  //OutputTree->Draw("vandle_qdc:vandle_corTof>>qcb(2000,-30,970,1000,0,20000)","dT<0 && vandle_qdc<banGate(vandle_corTof) && cloverAAB_E>1150 && cloverAAB_E<1165","colz");
  OutputTree->Draw("vandle_qdc:vandle_corTof>>qc( 2000,-30,970,1000,0,20000)","abs(dT)<500 && dT>0 && dr<0.1 ","colz");
  OutputTree->Draw("vandle_qdc:vandle_corTof>>qcb(2000,-30,970,1000,0,20000)","abs(dT)<500 && dT<0 && dr<0.1","colz");
  qc->Sumw2();
  qcb->Sumw2();
  TH2D* qca = (TH2D*)qc->Clone();
  qc->Add(qcb,-1);
  qc->Draw("colz");
  qc->ProjectionX();
  //qc->GetZaxis()->SetRangeUser(1,500);
  qc_px->Draw();
  int rebin = 1;
  qc_px->Rebin(rebin);
  qc_px->GetXaxis()->SetRangeUser(25,300);
  TString histName = fname+" Banana Gated Neutron Spectrum;Time of Flight (ns);Counts/2ns";
  qc_px->SetTitle(histName);
  //TF1 *bg = new TF1("bg","banGate(x)",20,970);
  //bg->Draw("same");
  TH2D* hh_qdc_tof_WalkC_BarC = (TH2D*)qc->Clone();
  hh_qdc_tof_WalkC_BarC->SetTitle("hh_qdc_tof_WalkC_BarC");
  hh_qdc_tof_WalkC_BarC->SetName("hh_qdc_tof_WalkC_BarC");
  TString fileName = "vandleHistograms/"+fname+"_vandle_hists.root";
  TFile *_file2 = TFile::Open(fileName,"UPDATE");
  hh_qdc_tof_WalkC_BarC->Write();
  _file2->Close();
}

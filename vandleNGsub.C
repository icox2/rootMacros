{
  double banGate(double x){return 200+4*(0.95*939.6*107*107/(2*x*x)-8.0*(1-exp(-1*0.1*pow(939.6*107*107/(2*x*x),0.9))));}
  TString fname = (TString)_file0->GetName();
  Int_t st = fname.First("_");
  Int_t end = fname.Length();
  fname.Remove(st,end-st);
  cout<<fname<<endl;
  //OutputTree->Draw("vandle_qdc:vandle_corTof>>qc(2000,-30,970,1000,0,20000)","dT>0 && vandle_qdc<banGate(vandle_corTof) && cloverAAB_E>1150 && cloverAAB_E<1165","colz");
  //OutputTree->Draw("vandle_qdc:vandle_corTof>>qcb(2000,-30,970,1000,0,20000)","dT<0 && vandle_qdc<banGate(vandle_corTof) && cloverAAB_E>1150 && cloverAAB_E<1165","colz");
  //OutputTree->Draw("vandle_ng_corTof>>qc(2000,-30,970)", "dT>0 && dr<0.1 && clover_ng_E>1330 && clover_ng_E<1344 && vandle_ng_qdc>50","colz");
  //OutputTree->Draw("vandle_ng_corTof>>qcb(2000,-30,970)","dT<0 && dr<0.1 && clover_ng_E>1330 && clover_ng_E<1344 && vandle_ng_qdc>50","colz");
  OutputTree->Draw("vandle_ng_corTof>>qc(2000,-30,970)", "dT>0 && dr<0.1 && clover_ng_E>794 && clover_ng_E<804 && vandle_ng_qdc>50","colz");
  OutputTree->Draw("vandle_ng_corTof>>qcb(2000,-30,970)","dT<0 && dr<0.1 && clover_ng_E>794 && clover_ng_E<804 && vandle_ng_qdc>50","colz");
  //OutputTree->Draw("vandle_ng_corTof>>qc(2000,-30,970)","dT>0 && dr<0.1 && clover_ng_E>538 && clover_ng_E<548 && vandle_ng_qdc>50","colz");
  //OutputTree->Draw("vandle_ng_corTof>>qcb(2000,-30,970)","dT<0 && dr<0.1 && clover_ng_E>538 && clover_ng_E<548 && vandle_ng_qdc>50","colz");
  qc->Sumw2();
  qcb->Sumw2();
  //qc->Add(qcb,-1);
  qc->Draw();
  int rebin = 4;
  //qc->Rebin(rebin);
  qc->GetXaxis()->SetRangeUser(25,300);
  //TString histName = fname+" Banana Gated Neutron Spectrum;Time of Flight (ns);Counts/"+0.5*rebin+"ns";
  //qc->SetTitle(histName);
  TH1D* proj = (TH1D*)qc->Clone();
  proj->SetTitle("proj");
  proj->SetName("proj");
  TString fileName = "vandleHistograms/"+fname+"_vandle_hists.root";
  //TFile *_file2 = TFile::Open(fileName,"UPDATE");
  //proj->Write();
  //qc->Write();
  //qcb->Write();
  //_file2->Close();
}

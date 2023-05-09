#include "Rtypes.h"
{
  double banGate(double x){return 200+4*(0.95*939.6*107*107/(2*x*x)-8.0*(1-exp(-1*0.1*pow(939.6*107*107/(2*x*x),0.9))));}
  TString fname = (TString)_file0->GetName();
  Int_t st = fname.First("_");
  Int_t end = fname.Length();
  fname.Remove(st,end-st);
  cout<<fname<<endl;
  int rebin = 1;
  Color_t a = 632;
  OutputTree->Draw("vandle_qdc:vandle_corTof>>qc(2000,-30,970,1000,0,20000)","dT>0 && vandle_qdc<banGate(vandle_corTof) && dr<0.1 && vandle_bar%2==0","colz");
  OutputTree->Draw("vandle_qdc:vandle_corTof>>qcb(2000,-30,970,1000,0,20000)","dT<0 && vandle_qdc<banGate(vandle_corTof)&& dr<0.1 && vandle_bar%2==0","colz");
  //OutputTree->Draw("vandle_qdc:vandle_corTof>>qc(2000,-30,970,1000,0,20000)","dT>0 && dr<0.1 ","colz");
  //OutputTree->Draw("vandle_qdc:vandle_corTof>>qcb(2000,-30,970,1000,0,20000)","dT<0 && dr<0.1","colz");
  qc->Sumw2();
  qcb->Sumw2();
  qc->Add(qcb,-1);
  qc->Draw("colz");
  qc->ProjectionX();
  qc_px->Rebin(rebin);
  qc_px->SetLineColor(1);
  qc_px->SetFillColor(4);
  qc_px->SetFillStyle(3001);
  qc_px->GetXaxis()->SetRangeUser(25,300);
  TString histName = fname+" Banana Gated Neutron Spectrum;Time of Flight (ns);Counts/"+to_string(0.5*rebin)+"ns";
  qc_px->SetTitle(histName);
  qc_px->Draw("hist");

  cout<<"Finished first background subtraction"<<endl;

  OutputTree->Draw("vandle_qdc:vandle_corTof>>qv(2000,-30,970,1000,0,20000)","dT>0 && vandle_qdc<banGate(vandle_corTof) && dr<0.1 && vandle_bar%2==1","colz");
  OutputTree->Draw("vandle_qdc:vandle_corTof>>qvb(2000,-30,970,1000,0,20000)","dT<0 && vandle_qdc<banGate(vandle_corTof)&& dr<0.1 && vandle_bar%2==1","colz");
  qv->Sumw2();
  qvb->Sumw2();
  qv->Add(qvb,-1);
  qv->Draw("colz");
  qv->ProjectionX();
  qv_px->Rebin(rebin);
  qv_px->SetLineColor(a+2);
  qv_px->SetFillColor(a);
  qv_px->SetFillStyle(3001);

  cout<<"Finished second background subtraction"<<endl;

  qc_px->Draw("hist");
  qv_px->Draw("same hist");
  //qv_px->Draw("same");
  //qc_px->Draw("same");
  qc_px->Draw("same AXIS");

  auto legend = new TLegend(0.6,0.7,0.9,0.9);
  legend->AddEntry(qc_px,"Front Layer","f");
  legend->AddEntry(qv_px,"Back Layer","f");
  legend->Draw();

}

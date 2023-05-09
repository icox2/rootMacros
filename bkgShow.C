#include "Rtypes.h"
#include <iostream>
{
  double banGate(double x){return 200+4*(0.95*939.6*107*107/(2*x*x)-8.0*(1-exp(-1*0.1*pow(939.6*107*107/(2*x*x),0.9))));}
  //gStyle->SetPalette();
  TString fname = (TString)_file0->GetName();
  Int_t st = fname.First("_");
  Int_t end = fname.Length();
  fname.Remove(st,end-st);
  cout<<fname<<endl;
  Color_t a = 632;
  //OutputTree->Draw("vandle_qdc:vandle_corTof>>qc(2000,-30,970,1000,0,20000)","dT>0 && vandle_qdc<banGate(vandle_corTof) && vandle_qdc>1000","colz");
  //OutputTree->Draw("vandle_qdc:vandle_corTof>>qcb(2000,-30,970,1000,0,20000)","dT<0 && vandle_qdc<banGate(vandle_corTof) && vandle_qdc>1000","colz");
  OutputTree->Draw("vandle_qdc:vandle_corTof>>qc(2000,-30,970,1000,0,20000)","dT>0 && vandle_qdc<banGate(vandle_corTof) ","colz");
  OutputTree->Draw("vandle_qdc:vandle_corTof>>qcb(2000,-30,970,1000,0,20000)","dT<0 && vandle_qdc<banGate(vandle_corTof)","colz");
  qc->Sumw2();
  qcb->Sumw2();
  qc->ProjectionX();
  qcb->ProjectionX();
  qc_px->Rebin(4);
  qc_px->SetLineColor(1);
  qc_px->SetFillColor(4);
  qc_px->SetFillStyle(3001);
  qc_px->Draw("hist");
  //TH1D *h1 = (TH1D*)qc_px->Clone();
  qcb_px->Rebin(4);
  qcb_px->SetLineColor(a+2);
  qcb_px->SetFillColor(a);
  qcb_px->SetFillStyle(3001);
  //TH1D *h2 = (TH1D*)qcb_px->Clone();
  qc_px->GetXaxis()->SetRangeUser(25,300);
  TString histName = fname+" Neutron Spectrum;Time of Flight (ns);Counts/2ns";
  qc_px->SetTitle(histName);
  qcb_px->Draw("same hist");
  qcb_px->Draw("same");
  qc_px->Draw("same");
  //qc_px->GetXaxis()->Draw();
  //qc_px->GetYaxis()->Draw();
  qc_px->Draw("same AXIS");
}

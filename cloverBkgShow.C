#include "Rtypes.h"
#include <iostream>
{
  //gStyle->SetPalette();
  TString fname = (TString)_file0->GetName();
  Int_t st = fname.First("_");
  Int_t end = fname.Length();
  fname.Remove(st,end-st);
  cout<<fname<<endl;
  Color_t a = 632;
  OutputTree->Draw("cloverAAB_E>>ce(3900,10,4000)","dT>0","colz");
  OutputTree->Draw("cloverAAB_E>>ceb(3900,10,4000)","dT<0","colz");
  ce->Sumw2();
  ceb->Sumw2();
  int rebin=1;
  ce->Rebin(rebin);
  ce->SetLineColor(4);
  ce->SetFillColor(4);
  ce->SetFillStyle(3001);
  ce->Draw("hist");
  //TH1D *h1 = (TH1D*)qc_px->Clone();
  ceb->Rebin(rebin);
  ceb->SetLineColor(a+2);
  ceb->SetFillColor(a);
  ceb->SetFillStyle(3001);
  //TH1D *h2 = (TH1D*)qcb_px->Clone();
  TString histName = fname+" AddBack #gamma Spectrum;Energy (keV);Counts/"+to_string(rebin)+"keV";
  ce->SetTitle(histName);
  ceb->Draw("same hist");
  ceb->Draw("same");
  //ce->Draw("same");
  //qc_px->GetXaxis()->Draw();
  //qc_px->GetYaxis()->Draw();
  ce->Draw("same AXIS");
}

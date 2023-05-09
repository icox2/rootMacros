#include "Rtypes.h"
#include <iostream>
{
  //gStyle->SetPalette();
  TString fName = (TString)_file0->GetName();
  Int_t start = fName.First("_");
  Int_t ending = fName.Length();
  fName.Remove(start,ending-start);
  cout<<fName<<endl;
  OutputTree->Draw("cloverAAB_E>>ces(3900,10,4000)","dT>0","colz");
  OutputTree->Draw("cloverAAB_E>>cesb(3900,10,4000)","dT<0","colz");
  ces->Sumw2();
  cesb->Sumw2();
  ces->Add(cesb,-1);
  ces->Draw("hist");
  bool reb=false;
  int reBins = 1;
  //ce->SetLineColor(4);
  //ce->SetFillColor(4);
  //ce->SetFillStyle(3001);
  //ce->Draw("hist");
  ////TH1D *h1 = (TH1D*)qc_px->Clone();
  if(reb==true){
    reBins=4;
    ces->Rebin(reBins);
  }
  //ceb->SetLineColor(a+2);
  //ceb->SetFillColor(a);
  //ceb->SetFillStyle(3001);
  ////TH1D *h2 = (TH1D*)qcb_px->Clone();
  TString histoName = fName+" Bkg Subtracted AddBack #gamma Spectrum;Energy (keV);Counts/"+to_string(reBins)+"keV";
  ces->SetTitle(histoName);
  //ceb->Draw("same hist");
  //ceb->Draw("same");
  ////ce->Draw("same");
  ////qc_px->GetXaxis()->Draw();
  ////qc_px->GetYaxis()->Draw();
  //ce->Draw("same AXIS");
  OutputTree->Draw("clover_E>>cens(3900,10,4000)","dT>0","colz");
  OutputTree->Draw("clover_E>>censb(3900,10,4000)","dT<0","colz");
  cens->Sumw2();
  censb->Sumw2();
  cens->Add(censb,-1);
  ces->SetLineColor(4);
  cens->SetLineColor(2);
  ces->Draw("hist");
  cens->Draw("same hist");
}

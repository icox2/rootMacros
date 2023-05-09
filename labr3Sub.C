#include "Rtypes.h"
#include <iostream>
{
  //gStyle->SetPalette();
  TString fName = (TString)_file0->GetName();
  Int_t start = fName.First("_");
  Int_t ending = fName.Length();
  fName.Remove(start,ending-start);
  cout<<fName<<endl;
  OutputTree->Draw("gammascint_E>>ces(3900,10,4000)","dT>0 && dr<0.1 ","colz");
  OutputTree->Draw("gammascint_E>>cesb(3900,10,4000)","dT<0 && dr<0.1","colz");
  ces->Sumw2();
  cesb->Sumw2();
  ces->Add(cesb,-1);
  ces->Draw("hist");
  int reBins = 1;
  //ce->SetLineColor(4);
  //ce->SetFillColor(4);
  //ce->SetFillStyle(3001);
  //ce->Draw("hist");
  ////TH1D *h1 = (TH1D*)qc_px->Clone();
  ces->Rebin(reBins);
  //ceb->SetLineColor(a+2);
  //ceb->SetFillColor(a);
  //ceb->SetFillStyle(3001);
  ////TH1D *h2 = (TH1D*)qcb_px->Clone();
  TString histoName = fName+" Bkg Subtracted LaBr_{3} #gamma Spectrum;Energy (keV);Counts/"+to_string(reBins)+"keV";
  ces->SetTitle(histoName);
  //ceb->Draw("same hist");
  //ceb->Draw("same");
  ////ce->Draw("same");
  ////qc_px->GetXaxis()->Draw();
  ////qc_px->GetYaxis()->Draw();
  //ce->Draw("same AXIS");
  //TString fileName = "vandleHistograms/"+fName+"_vandle_hists.root";
  //TFile *_file2 = TFile::Open(fileName,"UPDATE");
  //ces->Write();
  //_file2->Close();
}

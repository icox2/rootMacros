#include "Rtypes.h"
#include <iostream>
{
  //gStyle->SetPalette();
  TString fName = (TString)_file0->GetName();
  Int_t start = fName.First("_");
  Int_t ending = fName.Length();
  fName.Remove(start,ending-start);
  cout<<fName<<endl;
  OutputTree->Draw("gg1_E:gg2_E>>ggs(2000,10,4010,2000,10,4010)","dT>0 && dr<0.1 && abs(dT)<1800 && abs(gg1_T)<100 && abs(gg2_T)<100","colz");
  OutputTree->Draw("gg1_E:gg2_E>>ggb(2000,10,4010,2000,10,4010)","dT<0 && dr<0.1 && abs(dT)<1800 && abs(gg1_T)<100 && abs(gg2_T)<100","colz");

  ggs->Sumw2();
  ggb->Sumw2();
  ggs->Add(ggb,-1);
  ggs->Draw("colz hist");
  double Energy = 543;
  double tlow = Energy*0.994;
  double thigh = Energy*1.006;
  int nBins = ggs->GetNbinsX();
  double xMax = ggs->GetXaxis()->GetXmax();
  double xMin = ggs->GetXaxis()->GetXmin();
  double binning = (xMax-xMin)/((double)nBins);
  int firstxbin = (int)((tlow-xMin)/binning);
  int lastxbin = (int)((thigh-xMin)/binning);
  TH1D* ggs_py = ggs->ProjectionY("_py",firstxbin,lastxbin);
  TString histoName = fName+" #gamma="+to_string((int)Energy)+" keV Gated AddBack #gamma Spectrum;Energy (keV);Counts/keV";
  ggs_py->SetTitle(histoName);
  ggs_py->Draw("hist");
  
  // 0 - full matrix, 1 - Prompt Gamma peak, -1 - No saving
  int ind = -1;
  if(ind==0){
    TH2F *hh_ggmatrix = (TH2F*)ggs->Clone();
    hh_ggmatrix->SetTitle("hh_ggmatrix");
    hh_ggmatrix->SetName("hh_ggmatrix");
    TString fileName = "vandleHistograms/"+fName+"_vandle_hists.root";
    TFile *_file2 = TFile::Open(fileName,"UPDATE");
    hh_ggmatrix->Write();
    _file2->Close();
  }
  else if(ind==1){
    TH2F *hh_ggmatrixPrompt = (TH2F*)ggs->Clone();
    hh_ggmatrixPrompt->SetTitle("hh_ggmatrixPrompt");
    hh_ggmatrixPrompt->SetName( "hh_ggmatrixPrompt");
    TString fileName = "vandleHistograms/"+fName+"_vandle_hists.root";
    TFile *_file2 = TFile::Open(fileName,"UPDATE");
    hh_ggmatrixPrompt->Write();
    _file2->Close();
  }
  //TString fileName = "vandleHistograms/"+fName+"_vandle_hists.root";
  //TFile *_file2 = TFile::Open(fileName,"UPDATE");
  //ces->Write();
  //_file2->Close();
}

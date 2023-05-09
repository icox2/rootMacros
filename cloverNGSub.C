#include "Rtypes.h"
#include <iostream>
{
  //gStyle->SetPalette();
  TString fName = (TString)_file0->GetName();
  Int_t start = fName.First("_");
  Int_t ending = fName.Length();
  fName.Remove(start,ending-start);
  cout<<fName<<endl;
  OutputTree->Draw("clover_ng_E:vandle_ng_corTof>>ngs(2000,-30,970,3900,10,4000)"," abs(dT)<500 && dT>0 && dr<0.1 && vandle_ng_qdc>50","colz");
  OutputTree->Draw("clover_ng_E:vandle_ng_corTof>>ngsb(2000,-30,970,3900,10,4000)","abs(dT)<500 && dT<0 && dr<0.1 && vandle_ng_qdc>50","colz");
  ngs->Sumw2();
  ngsb->Sumw2();
  ngs->Add(ngsb,-1);
  ngs->Draw("colx");
  //int reBins = 1;
  //ce->SetLineColor(4);
  //ce->SetFillColor(4);
  //ce->SetFillStyle(3001);
  //ce->Draw("hist");
  ////TH1D *h1 = (TH1D*)qc_px->Clone();
  ///ngs->Rebin(reBins);
  //ceb->SetLineColor(a+2);
  //ceb->SetFillColor(a);
  //ceb->SetFillStyle(3001);
  ////TH1D *h2 = (TH1D*)qcb_px->Clone();
  //TString histoName = fName+" Bkg Subtracted AddBack #gamma Spectrum;Energy (keV);Counts/"+to_string(reBins)+"keV";
  //ngs->SetTitle(histoName);
  //ngs->SetName("hh_gABnmatrix_BetaGammaPrompt");
  //ngs->SetTitle("hh_gABnmatrix_BetaGammaPrompt");
  TH2D* hh_gABnmatrix_BetaGammaPrompt = (TH2D*)ngs->Clone();
  hh_gABnmatrix_BetaGammaPrompt->SetName("hh_gABnmatrix_BetaGammaPrompt");
  hh_gABnmatrix_BetaGammaPrompt->SetTitle("hh_gABnmatrix_BetaGammaPrompt");
  //ceb->Draw("same hist");
  //ceb->Draw("same");
  ////ce->Draw("same");
  ////qc_px->GetXaxis()->Draw();
  ////qc_px->GetYaxis()->Draw();
  //ce->Draw("same AXIS");
  TString fileName = "vandleHistograms/"+fName+"_vandle_hists.root";
  TFile *_file2 = TFile::Open(fileName,"UPDATE");
  hh_gABnmatrix_BetaGammaPrompt->Write();
  _file2->Close();
}

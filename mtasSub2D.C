#include "Rtypes.h"
#include <iostream>
{
  //gStyle->SetPalette();
  TCanvas *c1 = new TCanvas();
  TString fName = (TString)_file0->GetName();
  Int_t start = fName.First("_");
  Int_t ending = fName.Length();
  fName.Remove(start,ending-start);
  cout<<fName<<endl;
  OutputTree->Draw("mtas_E:mtas_sum>>ggs(1500,10,15010,1500,10,15010)","dT>0 && mtas_ringnum!=3 && mtas_ringnum!=4","colz");
  OutputTree->Draw("mtas_E:mtas_sum>>ggb(1500,10,15010,1500,10,15010)","dT<0 && mtas_ringnum!=3 && mtas_ringnum!=4","colz");
  ggs->Sumw2();
  ggb->Sumw2();
  ggs->Add(ggb,-1);
  ggs->Draw("colz hist");
  ggs->SetMinimum(0);
  c1->SetLogz(1);

  //TString histoName = fName+" Bkg Subtracted AddBack #gamma Spectrum;Energy (keV);Counts/"+to_string(reBins)+"keV";
  //ces->SetTitle(histoName);

}

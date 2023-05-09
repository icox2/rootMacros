#include "Rtypes.h"
#include <iostream>
{
  TCanvas *c1 = new TCanvas();
  TString fname = (TString)_file0->GetName();
  Int_t st = fname.First("_");
  Int_t end = fname.Length();
  fname.Remove(st,end-st);
  cout<<fname<<endl;
  Color_t a = 632;

  OutputTree->Draw("mtas_sum>>es(5000,10,10010)","dT>0","");
  OutputTree->Draw("mtas_sum>>eb(5000,10,10010)","dT<0","");
  //es->Sumw2();
  //eb->Sumw2();
  es->SetLineColor(4);
  es->SetFillColor(4);
  es->SetFillStyle(3001);
  es->Draw("hist");

  eb->SetLineColor(a+2);
  eb->SetFillColor(a);
  eb->SetFillStyle(3001);

  TString histName = fname+" MTAS Sum Energy Spectrum;Energy (keV);Counts/2keV";
  es->SetTitle(histName);

  eb->Draw("same hist");

  es->Draw("same AXIS");
}

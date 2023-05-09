#include "Rtypes.h"
#include <iostream>
#include "/home/icox/git/rootMacros/gfunc.h"
void mtasSubSum(){
  //TCanvas *c1 = new TCanvas();
  TString fname = (TString)_file0->GetName();
  Int_t st = fname.First("_");
  Int_t end = fname.Length();
  fname.Remove(st,end-st);
  cout<<fname<<endl;
  TTree *OutputTree = (TTree*)_file0->Get("OutputTree");

  bool ifFit = 0;
  double Sn = 5168.9; // Sn of 45Cl in keV
  double bkgRatio = 1.;

  TH1D *es = new TH1D("es","es",1500,10,15010);
  TH1D *eb = new TH1D("eb","eb",1500,10,15010);
  OutputTree->Draw("mtas_sum>>es","dT>0 && abs(dT)<1500","");
  OutputTree->Draw("mtas_sum>>eb","dT<0 && abs(dT)<1500","");
  //OutputTree->Draw("mtas_sum>>es","betaEnergySum<7000 && dT>0 && abs(mtas_t)<75 && abs(dT)<1500","");
  //OutputTree->Draw("mtas_sum>>eb","betaEnergySum<7000 && dT<0 && abs(mtas_t)<75 && abs(dT)<1500","");
  es->Sumw2();
  eb->Sumw2();
  es->Add(eb,-1/bkgRatio);
  int rebin = 1;
  es->Rebin(rebin);
  TString histName = fname+" MTAS Sum Energy Spectrum;Energy (keV);Counts/"+to_string(2*rebin)+"keV";
  es->SetTitle(histName);
  es->Draw("hist");
  //es->GetYaxis()->SetRangeUser(-5,400);

  if(ifFit){  
    const int numPeaks = 10;
    double Peaks[numPeaks] = {542,850,1340,1417,1736,1771,2757,3295,3950,4325};
    TF1 *g = new TF1("g",gfunc,300,5200,2*numPeaks);
    for(int i=0;i<numPeaks;i++) {
      if(i!=1 && i!=2 && i!=3 && i!=9)
        g->FixParameter(2*i+1,Peaks[i]);
      else
        g->SetParameter(2*i+1,Peaks[i]);
      g->SetParLimits(2*i,0,10000);
      g->SetParameter(2*i,100);
    }
    es->Fit(g,"RL","same");
    g->Draw("same");
    for(int i=0;i<numPeaks;i++){
      double counts = 0;
      double fwhmRatio = 0.1;
      double con = g->GetParameter(2*i);
      double sig = g->GetParameter(2*i+1)*fwhmRatio/2.3548;
      counts = con*sqrt(TMath::Pi()*2*sig*sig)/(2*rebin);
      cout<<"Number of counts at "<<g->GetParameter(2*i+1)<<" keV = "<<counts<<endl;
    }
  }
  cout<<"Total counts in hist: "<<es->Integral(0,5000)<<endl;
}

#include "Rtypes.h"
#include <iostream>
#include "/home/icox/git/rootMacros/gfunc.h"
void mtasSumCompENSDF(){
  TCanvas *c1 = new TCanvas();
  TLegend *leg = new TLegend(0.6,0.7,0.85,0.85);
  TString fname = (TString)_file0->GetName();
  Int_t st = fname.First("_");
  Int_t end = fname.Length();
  fname.Remove(st,end-st);
  cout<<fname<<endl;
  TTree *OutputTree = (TTree*)_file0->Get("OutputTree");

  bool ifFit = 0;
  double Sn = 5168.9; // Sn of 45Cl in keV
  double Pn = 0.27;
  double numBetas = 11128;
  double bkgRatio = 4;
  double maxE = 15010;
  double minE = 10;
  double nBins = 1500;
  int binsPerkeV = (int)(maxE-minE)/nBins;

  TH1D *es = new TH1D("es","es",nBins,minE,maxE);
  TH1D *eb = new TH1D("eb","eb",nBins,minE,maxE);
  OutputTree->Draw("mtas_sum>>es","dT>0","");
  OutputTree->Draw("mtas_sum>>eb","dT<0","");
  //OutputTree->Draw("mtas_sum>>es","betaEnergySum<7000 && dT>0 && abs(mtas_t)<75","");
  //OutputTree->Draw("mtas_sum>>eb","betaEnergySum<7000 && dT<0 && abs(mtas_t)<75","");
  es->Sumw2();
  eb->Sumw2();
  es->Add(eb,-1/bkgRatio);
  int rebin = 10;
  es->Rebin(rebin);
  TString histName = fname+" MTAS Sum Energy Spectrum;Energy (keV);Counts/"+to_string(binsPerkeV*rebin)+"keV";
  es->SetTitle(histName);
  es->SetLineColor(1);
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
  cout<<"Total counts in data hist: "<<es->Integral(1,nBins/rebin)<<endl;

  TFile *_file1 = TFile::Open("/SCRATCH/DScratch3/e21069/simResults/Cl45Test.root");
  //TFile *_file1 = TFile::Open("/SCRATCH/DScratch3/e21069/simResults/Cl45notoplevels.root");
  bool Pnscaling = 1;
  TTree *EventInfo = (TTree*)_file1->Get("EventInfo");
  TH1D *t = new TH1D("t","t",nBins/rebin,minE,maxE);
  EventInfo->Draw("Total>>t","","");
  double numSims = t->Integral(1,nBins/rebin);
  double scale = numBetas/numSims*(Pnscaling ? (1-Pn) : 1);
  t->Scale(scale);
  t->SetLineColor(2);
  TString hName = fname+" MTAS Simulation Comp;Energy (keV);Counts/"+to_string(binsPerkeV*rebin)+"keV";
  t->SetTitle(hName);
  cout<<"Total counts in simulated hist: "<<t->Integral(1,nBins/rebin)<<endl;


  t->Draw("hist");
  t->GetYaxis()->SetRangeUser(-100,1100);
  leg->AddEntry(t,"Simulated ENSDF","l");
  _file0->cd();
  es->Draw("hist same");
  leg->AddEntry(es,"MTAS Sum Data","l");

  leg->Draw();
}

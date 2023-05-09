// pulserTiming.cpp
#include "TROOT.h"
#include "TGraph2D.h"
#include "TCanvas.h"

{
  const int numFiles = 6;
  TCanvas *c = new TCanvas();
  c->Divide(3,2);
  char* files[numFiles] = {"pulser015014R.dat","pulser015012R.dat","pulser015000R.dat","pulser015112R.dat","pulser015115R.dat","pulser115114.dat"};
  TGraph2D *graphs[numFiles];
  for(int i=0;i<numFiles;i++){
    graphs[i] = new TGraph2D(files[i]);
    c->cd(i+1);
    graphs[i]->Draw("colz");
  }
  c->cd();
}
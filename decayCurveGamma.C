{
  double energy = 542.5;
  double range = 3.5;
  double tMax = 1800;
  int nBins = 300;
  // Histograms
  TH1F *d = new TH1F("d","d",nBins,-1*tMax,tMax);
  TH1F *neg = new TH1F("neg","neg",nBins,-1*tMax,tMax);
  TH1F *pos = new TH1F("pos","pos",nBins,-1*tMax,tMax);
  //Gates
  TString uGate = "dr<0.1 && ";
  TString dGate = "cloverAAB_E>"+to_string((int)(energy-range))+" && cloverAAB_E<"+to_string((int)(energy+range));
  TString nGate = "cloverAAB_E>"+to_string((int)(energy-3*range))+" && cloverAAB_E<"+to_string((int)(energy-range));
  TString pGate = "cloverAAB_E>"+to_string((int)(energy+range))+" && cloverAAB_E<"+to_string((int)(energy+3*range));
  //Plotting
  OutputTree->Draw("dT>>d",uGate+dGate,"");
  OutputTree->Draw("dT>>neg",uGate+nGate,"");
  OutputTree->Draw("dT>>pos",uGate+pGate,"");

  d->Sumw2();
  neg->Sumw2();
  pos->Sumw2();
  d->Add(neg,-0.5);
  d->Add(pos,-0.5);

  d->Draw();
}

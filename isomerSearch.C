// isomerSearch.C
// Macro for analyzing pspmt tree for isomerSearching

double cloverWC(double cloverTime, double dlTime, double cloverEnergy){
  return cloverTime-dlTime-1.83824e+03*pow(cloverEnergy,-4.97513e-01)+1.17211e+02;
}

void isomerSearch(TTree *tree, double ratio=0){
  std::cout<<"Starting"<<std::endl;
  TCanvas *c1 = new TCanvas;
  // Histograms for Xu's analysis, tree==OutputTree
  // TH2F *pidRaw = new TH2F("pidRaw","pidRaw",500,2.3,3.5,500,5,16);
  // TH2F *pidGate = new TH2F("pidGate","pidGate",500,2.3,3.5,500,5,16);
  // tree->Draw("Zed:AoQ>>pidRaw","","colz");
  // tree->Draw("Zed:AoQ>>pidGate","clover_Twc>1000 && clover_Twc<45000","colz");
  // Histograms for trace analyzer data, tree=pspmt (or sipm?)
  TH2F *pidRaw = new TH2F("pidRaw","pidRaw",500,-260,-210,500,2000,6000);
  TH2F *pidGate = new TH2F("pidGate","pidGate",500,-260,-210,500,2000,6000);
  tree->Draw("pid_vec_.cross_pin_1_energy:pid_vec_.tof2>>pidRaw","","colz");
  tree->Draw("pid_vec_.cross_pin_1_energy:pid_vec_.tof2>>pidGate","cloverWC(clover_vec_.time,low_gain_.time_,clover_vec_.energy)>400 && cloverWC(clover_vec_.time,low_gain_.time_,clover_vec_.energy)<4500 && low_gain_.time_>0","colz");
  double rawMax = pidRaw->GetMaximum();
  // double rawMax = pidRaw->Integral();
  std::cout<<"Bin loc: "<<pidRaw->GetMaximumBin()<<" Bin cent: "<<pidRaw->GetBinContent(pidRaw->GetMaximumBin())<<" Hist max: "<<rawMax<<std::endl;
  // double gateMax = pidGate->GetMaximum();
  double gateMax = pidGate->GetBinContent(pidRaw->GetMaximumBin());
  // double gateMax = pidGate->Integral();
  std::cout<<"Scaling factor is: "<<gateMax/rawMax<<std::endl;
  std::cout<<"Bin loc: "<<pidGate->GetMaximumBin()<<" Bin cent: "<<pidGate->GetBinContent(pidGate->GetMaximumBin())<<" Hist max: "<<gateMax<<std::endl;
  TH2F *pidSub = (TH2F*)pidGate->Clone("pidSub");
  if(ratio==0)
    ratio = gateMax/rawMax;
  pidSub->Add(pidRaw,-1*ratio);
  pidSub->Draw("colz");
  return;
}

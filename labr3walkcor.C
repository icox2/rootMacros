{
  TCanvas *c1 = new TCanvas();

  TString fname = (TString)_file0->GetName();
  Int_t st = fname.First("_");
  Int_t end = fname.Length();
  fname.Remove(0,st+1);
  st = fname.First("_");
  fname.Remove(0,st+1);
  st = fname.First("_");
  end = fname.Length();
  fname.Remove(st,end-st);
  int runnum = atoi(fname.Data());
  cout<<"Run number: "<<fname<<" "<<runnum<<endl;

  TString outName = "/home/icox/git/rootMacros/labr3tofShifts.dat";
  cout<<"Writing values to: "<<outName<<endl;
  ofstream fout(outName.Data(),std::ofstream::app);
 
  vector<double> means, sigmas;
  vector<TH1F *> hv;
  for(int i=0;i<15;i++){
    cout<<"det: "<<i<<endl;
    TString gate = TString("output_vec_.detNum=="+to_string(i)); 
    TString var = TString("output_vec_.time-input_.input_.high_gain_.hrtime_>>t"+to_string(i));
    TH1F *temp = new TH1F(TString("t"+to_string(i)),TString("t"+to_string(i)),5000,-900,-400);
    twoCrateAfterTree->Draw(var,gate,"");
    hv.push_back(temp);

    int maxBin = temp->GetMaximumBin();
    int nBins = temp->GetNbinsX();
    double xMax = temp->GetXaxis()->GetXmax();
    double xMin = temp->GetXaxis()->GetXmin();
    double binning = (xMax-xMin)/((double)nBins);
    double hMax = maxBin*binning+xMin;
    cout<<"Detector "<<i<<" hist max at: "<<hMax<<endl;

    TF1 *gfit = new TF1("gfit","gaus",0,1);
    double fitWidth = 1.0;
    temp->Fit(gfit,"R","",hMax-fitWidth,hMax+fitWidth);
    means.push_back(gfit->GetParameter(1));
    sigmas.push_back(gfit->GetParameter(2));

    delete gfit;
    delete temp;
  }
  

  fout<<runnum<<" ";
  cout<<runnum<<" ";
  for(int i=0;i<means.size();i++){
    fout<<means.at(i)<<" ";
    cout<<means.at(i)<<" ";
  }
  fout<<endl;
  cout<<endl;

  return;
}

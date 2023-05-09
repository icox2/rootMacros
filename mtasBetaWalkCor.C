{
  TH2D *bh = new TH2D("bh","bh",300,-600,600,5000,0,5000);
  OutputTree->Draw("mtas_E:mtas_t>>bh","","colz");
  TF1 *func = new TF1("func","[0]*pow(abs(x),[1])+[2]");
  vector<double> par0, par1, par2;
  for(int i=0;i<25;i++){
    TString gate = "mtas_segID=="+to_string(i);
    OutputTree->Draw("mtas_E:mtas_t>>bh",gate,"colz");
    TH1D *prof = bh->ProfileY();
    prof->Draw("hist");
    func->SetParameter(0,600);
    func->SetParameter(1,-0.3);
    func->SetParameter(2,-30.);
    prof->Fit(func,"R","",5,3000);
    par0.push_back(func->GetParameter(0));
    par1.push_back(func->GetParameter(1));
    par2.push_back(func->GetParameter(2));
  }
  cout<<"double twcPar0["<<par0.size()<<"] = {";
  for(int i=0;i<par0.size();i++){
    if(i!=0) cout<<",";
    cout<<par0.at(i);
  }
  cout<<"};"<<endl;
  cout<<"double twcPar1["<<par0.size()<<"] = {";
  for(int i=0;i<par1.size();i++){
    if(i!=0) cout<<",";
    cout<<par1.at(i);
  }
  cout<<"};"<<endl;
  cout<<"double twcPar2["<<par2.size()<<"] = {";
  for(int i=0;i<par2.size();i++){
    if(i!=0) cout<<",";
    cout<<par2.at(i);
  }
  cout<<"};"<<endl;
}

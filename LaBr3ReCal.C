{
  ifstream fin("FDSiLaBr3ReCalData.dat");
  vector<double> calE;
  vector<vector< double> > inputE;
  const int numHagrids = 15;
  const int numDataPoints = 8;
  string junk;
  int cnt = 0;
  
  // Get the data points from the file
  while(fin.good()){
    if(cnt==0){
      double n=0;
      fin>>junk>>junk;
      for(int i=0;i<numHagrids;i++)
        fin>>n;
      cnt++;
    }
    else{
      double E,oldE;
      vector<double> inp;
      fin>>E;
      calE.push_back(E);
      for(int i=0;i<numHagrids;i++){
        fin>>oldE;
        inp.push_back(oldE);
      }
      inputE.push_back(inp);
      cnt++;
    }
  }

  // Put the data into TGraphs for fitting
  vector<TGraph*> hagrids;
  TMultiGraph *graph = new TMultiGraph();
  vector<TF1*> fits;
  for(int i=0;i<numHagrids;i++){
    TGraph *tmp = new TGraph(numDataPoints);
    for(int j=0;j<numDataPoints;j++)
      tmp->SetPoint(j,inputE.at(j).at(i),calE.at(j));
    tmp->SetMarkerStyle(8);
    tmp->SetMarkerColor(i+1);
    TF1 *pol3 = new TF1("pol3","pol3",0,1500);
    tmp->Fit(pol3,"Q");
    fits.push_back(pol3);
    graph->Add(tmp);
    hagrids.push_back(tmp);
  }
  cout<<"Recalibration Parameters using a pol3"<<endl;
  for(int i=0;i<fits.at(0)->GetNpar();i++){
    cout<<"double par"<<i<<"["<<fits.size()<<"] = {";
    for(int j=0;j<fits.size();j++){
      if(j!=0) cout<<",";
      cout<<fits.at(j)->GetParameter(i);
    }
    cout<<"};"<<endl;
  }

  graph->Draw("ap");
}

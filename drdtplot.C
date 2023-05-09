{
  OutputTree->Draw("2*dr:dT>>rt(600,-300,300,1000,0,0.2)","","colz");
  rt->SetTitle(";Decay Time (ms); Correlation Radius (a.u.)");
}

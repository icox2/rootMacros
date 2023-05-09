{ 
  TString fName = (TString)_file0->GetName();
  Int_t start = fName.First("_");
  Int_t ending = fName.Length();
  fName.Remove(start,ending-start);
  cout<<fName<<endl;
  OutputTree->Draw("cloverAAB_E:cloverAAB_Twc>>ces(4000,-2000,2000,5000,10,5010)","dT>0 && dr<0.1","colz");  
  OutputTree->Draw("cloverAAB_E:cloverAAB_Twc>>ceb(4000,-2000,2000,5000,10,5010)","dT<0 && dr<0.1","colz");
  ces->Sumw2();
  ceb->Sumw2();
  ces->Add(ceb,-1);
  TString histoName = fName+" Bkg Subtracted AddBack #gamma Spectrum;Clover-#beta Time (ns);Clover Energy (keV)";
  ces->SetTitle(histoName);
  ces->Draw("colz");

  // Saving The 2D Histogram
  //TH2F *hh_gammaAAB_energytime = (TH2F*)ces->Clone();
  //hh_gammaAAB_energytime->SetTitle("hh_gammaAAB_energytime");
  //hh_gammaAAB_energytime->SetName("hh_gammaAAB_energytime");
  //TString fileName = "vandleHistograms/"+fName+"_vandle_hists.root";
  //TFile *_file2 = TFile::Open(fileName,"UPDATE");
  //hh_gammaAAB_energytime->Write();
  //_file2->Close();
}

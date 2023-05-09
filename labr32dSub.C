// labr32DSub.C

{
  TCanvas *c1 = new TCanvas();

  OutputTree->Draw("gammascint_E:gammascint_Twc>>gt(2000,-10,10,2000,0,10000)","dT>0","colz");
  OutputTree->Draw("gammascint_E:gammascint_Twc>>gb(2000,-10,10,2000,0,10000)","dT<0","colz");

  gt->Sumw2();
  gb->Sumw2();

  gt->Add(gb,-1);
  gt->Draw("colz");
  gt->GetZaxis()->SetRangeUser(1,50);


  return;
}

{
    mergedBeta->Draw("2*sqrt(pow(output_vec_.low_gain_.pos_x_-input_.high_gain_.pos_x_,2)+pow(output_vec_.low_gain_.pos_y_-input_.high_gain_.pos_y_,2))>>rs2(1000,0,0.2)","input_.external_ts_high_-output_vec_.external_ts_high_>0","");
    mergedBeta->Draw("2*sqrt(pow(output_vec_.low_gain_.pos_x_-input_.high_gain_.pos_x_,2)+pow(output_vec_.low_gain_.pos_y_-input_.high_gain_.pos_y_,2))>>rb(1000,0,0.2)","input_.external_ts_high_-output_vec_.external_ts_high_<0","");
    rs2->Add(rb,-1);
    rs2->Draw();
    rs2->SetTitle("^{37}Al Correlation Radius;Radius (a.u.);Counts");
    TFile *_file2 = TFile::Open("~/histograms/imCorr_hists.root","UPDATE");
    rs2->Write();
    _file2->Close();
}

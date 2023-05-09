#include "ReadVandleSetup.C"

double FitGaussian(TH1D* hist, double fit_low, double fit_high){
   double result=0;
   TF1 *f1 = new TF1("f1","gaus",fit_low,fit_high);
   hist->Fit(f1,"RQ");
   result = f1->GetParameter(1);
   return result;
}

void timeCal(TTree *tree, TString name){
   //TFile *file = new TFile(name,"read");
   ReadVandleSetup("/home/icox/rootMacros/medium.dat");
   //TTree *tree = (TTree*)file->Get("PixTree");
   TH2F *htof2D = new TH2F("htof2D","htof2D",1000,-100,900,90,0,90);
   TH2F *htdiff2D = new TH2F("htdiff2D","htdiff2D",100,-50,50,90,0,90);
   tree->Draw("vandle_vec_.barNum:vandle_vec_.tof>>htof2D","","colz");
   tree->Draw("vandle_vec_.barNum:vandle_vec_.tDiff>>htdiff2D","","colz");
   ofstream myfile;
   char myfilename[50];
   //sprintf(myfilename,"TimeCal_medium_%s.xml","test");
   myfile.open("/SCRATCH/DScratch7/e21069rootfile/utkscan/crate1/tofConfigs/"+name+".xml");
   myfile<<"<medium>"<<endl;
   for(int i=0; i<NumOfBars; i++){
      char tdiffname[200];
      char tofname[200];
      sprintf(tdiffname, "tdiff_%d", i);
      sprintf(tofname, "tof_%d", i);
      TH1D *tdiff_bar = htdiff2D->ProjectionX(tdiffname, i+1, i+1);
      TH1D *tof_bar = htof2D->ProjectionX(tofname, i+1, i+1);
      //tof_bar->GetXaxis()->SetRangeUser(500,1050);
      double tdiff_max = tdiff_bar->GetMaximum();
      double tdiff_maxPos = tdiff_bar->GetXaxis()->GetBinCenter(tdiff_bar->GetMaximumBin());
      double tof_max = tof_bar->GetMaximum();
      double tof_maxPos = tof_bar->GetXaxis()->GetBinCenter(tof_bar->GetMaximumBin());
      //cout<<"For bar #"<<i<<": maximum tof count "<<tof_max<<" at "<<tof_maxPos<<endl;
      //cout<<"         "<<i<<": maximum tdiff count "<<tdiff_max<<" at "<<tdiff_maxPos<<endl;
      double fit_low=tdiff_maxPos-30;
      double fit_high=tdiff_maxPos+30;
      double meanTdiff = FitGaussian(tdiff_bar,fit_low,fit_high);
      //cout<<"For bar #"<<i<<": max tdiff at "<<tdiff_maxPos<<", mean tdiff at "<<meanTdiff<<endl;
      fit_low=tof_maxPos-2;
      fit_high=tof_maxPos+2;
      double gflash = FitGaussian(tof_bar,fit_low,fit_high);
      //cout<<"        #"<<i<<": max tof at "<<tof_maxPos<<", gamma flash at "<<gflash<<endl;
      double tdiffOffset = -1*meanTdiff;
      double gflashBin = VANDLEZ0[i]/speedOfLight; // Expected gamma flash location
      double tofOffset = (gflashBin-gflash);
      cout<<"For bar #"<<i<<": tDiff offset = "<<tdiffOffset<<endl;
      cout<<"        #"<<i<<": gamma-flash offset = "<<tofOffset<<endl;
      delete tdiff_bar;
      delete tof_bar;
      myfile<<"   <Bar number=\""<<i<<"\" lroffset=\""<<tdiffOffset<<"\" z0=\""<<VANDLEZ0[i]<<"\" xoffset=\""<<VANDLEXOffset[i]<<"\" zoffset=\"0\">"<<endl;
      myfile<<"      <TofOffset location=\"0\" offset=\""<<tofOffset<<"\"/>"<<endl;
      myfile<<"   </Bar>"<<endl;
   }
   myfile<<"</medium>"<<endl;
   myfile.close();
}

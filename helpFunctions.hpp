#include <fstream>

#ifndef HELPFUNCTIONS
#define HELPFUNCTIONS

void histToTxt(TH1 *hist, TString filename){
  ofstream fout;
  fout.open(filename.Data());
  double x,y;
  for(int i=0;i<hist->GetNbinsX();++i){
    x = hist->GetBinCenter(i);
    y = hist->GetBinContent(i);
    fout<<x<<" "<<y<<endl;
  }
  return;
}

#endif

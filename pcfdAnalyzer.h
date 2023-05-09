#include <cmath>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TROOT.h"
#include <vector>

TF1 *fit = new TF1("fit","[0]*exp(([1]-x)/[2])*(1-exp(-1*pow(x-[1],4)/[3]))+[4]",0,300);
TF1 *fpol1 = new TF1("fName1","pol1",0,300);
TF1 *fpol3 = new TF1("fName3","pol3",0,1);

pair<double,double> pcfdAnalyzer(vector<unsigned int> traceU, double baseline){
  int maxPos=0;
  vector<double> trace;
  for(int i=0;i<traceU.size();i++) trace.emplace_back(traceU.at(i));
  vector<double>::iterator it;
  TGraph *fTraces = new TGraph();
  double fmax=0.0;
  double pmax=0.0;
  //double baseline=0.0;
  double thresh=0.0;
  double phase=0.0;
  std::pair <double,double> points(0,0);

  if(trace.size()>10){
    it = max_element(trace.begin(),trace.end());
    maxPos = distance(trace.begin(),it);
  }
  else{
    return make_pair(-9999.,-9999.);
  }
  for (int iv = 0;iv<trace.size();iv++){
    fTraces->SetPoint(iv, iv, trace.at(iv));
  }
  
  std::pair <double, double> range((maxPos-2),(maxPos+1));
  fpol3->SetRange(range.first,range.second);
  fTraces->Fit(fpol3,"RNQSW");
  if(!(fpol3->IsValid()))return make_pair(-5555.,-5555.);
  fmax = fpol3->GetMaximum(range.first,range.second);
  pmax = trace.at(maxPos);
  //baseline = baselineCalc(trace).first;
  //if(pmax<3*baselineCalc(trace).second)cout<<"peak too small "<<baselineCalc(trace).second<<endl;
  //if(maxPos>125)return make_pair(-7777.,-7777.);
  thresh = (fmax-baseline)*0.4+baseline; //here 0.4 is the Constant Fraction
  for(int i=maxPos;i>0;i--){
    if((trace.at(i)>=thresh) && (trace.at(i-1)<thresh)){
      points.first = i-1;
      points.second = i; //setting points of trace around threshold 
    }
  }
  fpol1->SetRange(points.first,points.second);
  fTraces->Fit(fpol1,"RNQSW");
  phase = fpol1->GetX(thresh,points.first,points.second);
  return make_pair(phase,fmax);
}

//definition for cfd analysis
#include "cfdTiming.hpp"

pair <double, double> baselineCalc(vector<unsigned int> trace);

double cfdTiming(vector<unsigned int> trace){
  int maxPos=0;
  vector<unsigned int>::iterator it;
  TGraph *fTraces = new TGraph();
  TF1 *fpol3 = new TF1("fName3","pol3",0,1);
  TF1 *fpol1 = new TF1("fName1","pol1",0,1);
  double fmax=0.0;
  double pmax=0.0;
  double baseline=0.0;
  double thresh=0.0;
  double phase=0.0;
  std::pair <double,double> points(0,0);

  if(trace.size()>10){
    it = max_element(trace.begin(),trace.end());
    maxPos = distance(trace.begin(),it);
  }
  else{
    return -9999.;
  }
  for (int iv = 0;iv<trace.size();iv++){
    fTraces->SetPoint(iv, iv, trace.at(iv));
  }

  std::pair <double, double> range((maxPos-1),(maxPos+2));
  fpol3->SetRange(range.first,range.second);
  fTraces->Fit(fpol3,"RNQSW");
  fmax = fpol3->GetMaximum(range.first,range.second);
  pmax = trace.at(maxPos);
  baseline = baselineCalc(trace).first;
  thresh = (fmax-baseline)*0.5+baseline; //here 0.5 is the Constant Fraction
  for(int i=maxPos;i>0;i--){
    if((trace.at(i)>=thresh) && (trace.at(i-1)<thresh)){
      points.first = i-1;
      points.second = i; //setting points of trace around threshold 
    }
  }
  fpol1->SetRange(points.first,points.second);
  fTraces->Fit(fpol1,"RNQSW");
  phase = fpol1->GetX(thresh,points.first,points.second);
  return phase;
}

pair <double, double> baselineCalc(vector<unsigned int> trace){
    //calculating the baseline
    double baseSum = 0.;
    for(int j=0;j<20;j++){
        baseSum += trace.at(j);
    }
    double  baseline = baseSum/20.;
    
    //calculating the standard dev
    double stddev = 0.0;
    for(int j=0;j<20;j++){
        stddev += pow(trace[j]-baseline,2);
    }
    stddev = sqrt(stddev/20);
    
    return std::make_pair (baseline, stddev);
}
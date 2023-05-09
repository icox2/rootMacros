//functions.hpp
//Just a bunch of random functions for trace analysis
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
using namespace std;

#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

pair<double, double> baselineCalc(vector<unsigned int> trace){
    //calculating the baseline
    int tSize = trace.size();
    double baseSum = 0.;
    for(int j=0;j<35;j++){
        baseSum += trace[j];
    }
    double  baseline = baseSum/35.;
    
    //calculating the standard dev
    double stddev = 0.0;
    for(int j=0;j<35;j++){
        stddev += pow(trace[j]-baseline,2);
    }
    stddev = sqrt(stddev/35);
    
    return std::make_pair (baseline, stddev);
}
pair<double, double> baselineCalc(vector<double> trace){
    //calculating the baseline
    int tSize = trace.size();
    double baseSum = 0.;
    for(int j=0;j<35;j++){
        baseSum += trace[j];
    }
    double  baseline = baseSum/35.;
    
    //calculating the standard dev
    double stddev = 0.0;
    for(int j=0;j<35;j++){
        stddev += pow(trace[j]-baseline,2);
    }
    stddev = sqrt(stddev/35);
    
    return std::make_pair (baseline, stddev);
}

vector<double> normalization(vector<double> tr){
  double base = baselineCalc(tr).first;
  for(int i=0;i<tr.size();i++)
    tr.at(i)-= base;
  double max = *max_element(tr.begin(), tr.end());
  vector<double> results;
  for(int i=0;i<tr.size();i++)
    results.push_back(tr.at(i)/max);
  return results;
}

double QDCcalculator(vector<unsigned int> trace, unsigned int lowBound, unsigned int highBound){
    double  baseline = baselineCalc(trace).first;
   //calculate integral after baseline
    double integral = 0.0;
    for (int i=lowBound; i<highBound;i++){
        //integral += 0.5 * (abs(double(trace[i-1]) - baseline) + abs(double(trace[i]) - baseline));
        integral += abs(double(trace[i])-baseline);
    }
    return integral;
}
double QDCcalculator(vector<double> trace){
    double qdc=0.0;
    unsigned int highbound, lowbound;
    double baseline = baselineCalc(trace).first;
    int maxPos=0;
    // auto maxPos = max_element(trace.begin(),trace.end());
    vector<double>::iterator it;

    if(trace.size()>10){
      //it = max_element(trace.begin(),trace.end());
      //maxPos = distance(trace.begin(),it);
      maxPos = max_element(trace.begin(),trace.end()) - trace.begin();
    }
    // else if(maxPos<180){
    //   return -4444.;
    // }
    else{
      return -9999.;
    }

    highbound = maxPos+40;
    lowbound = maxPos-20;

    for(int i=lowbound;i<highbound;i++)
      qdc += 0.5 * double(abs(trace[i-1] - baseline) + abs(trace[i] - baseline));

    return qdc;
}


vector<double> der(vector<double> tr){
  vector<double> derivative;
  derivative.push_back(0);
  for(int i=0;i<tr.size()-1;i++)
    derivative.push_back(tr.at(i+1)-tr.at(i));
  return derivative;
}
#endif
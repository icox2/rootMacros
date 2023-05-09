#include "TSystem.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TChain.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <ProcessorRootStruc.hpp>
#include <cmath>
#include "position.h"

pair<double,double> pcfdAnalyzer(vector<double> trace,double thresh);
double regCfd(vector<double> &trace);
double derAnalyzer(vector<double> &trace);
pair <double, double> baselineCalc(vector<double> trace);
double traceAnalyzer(vector<double> trace);
double fittingAnalyzer(vector<double> trace);
double gFit(vector<double> trace);
vector<double> reducer(vector<double> trace);
double centroid(vector<double> trace);
TF1 *fit = new TF1("fit","[0]*exp(([1]-x)/[2])*(1-exp(-1*pow(x-[1],[5])/[3]))+[4]",0,300);
TF1 *fpol1 = new TF1("fName1","pol1",0,300);
TF1 *fpol2 = new TF1("fName2","pol2",0,1);
TF1 *fpol3 = new TF1("fName3","pol3",0,1);
TF1 *gauss = new TF1("guass","gaus",0,1);

TFile *newFile = new TFile("gaggLaBr3TimingLong10Amp.root","UPDATE");
TTree *timing = new TTree("timing","tree filled with traces and energies etc.");

void longAnalysis(){
  TChain ch("PixTree");
  string prefix = "gagg60CoLaBr3TimingLG10Amp_001-", ending = "_DD.root", fileName;
  int numfiles = 8;
  vector<TFile*> _file0;
  _file0.clear();
  
  TH1D *plot1=new TH1D("plot1","plot1",10000.,-2000.,2000.);
  TH2D *energy=new TH2D("energy","energy",5000.,0.,10000.,5000.,0.,10000.);
  TH2D *phase=new TH2D("phase","phase",1000.,80.,105.,1000.,80.,105.);

    //variables
    int eventNum = 0;
    double energyone=0.0, energytwo=0.0, energythree=0.0;
    vector<double> traceone, tracetwo, xatrace, xbtrace, yatrace, ybtrace;
    vector<unsigned int> tone, ttwo, xat, xbt, yat, ybt;
    int ione=0, itwo=0, ithree=0, twocount=0;
    double hrTimeone=0.0, hrTimetwo=0.0, hrTimethree=0.0;
    double cfdTimeone=0.0, cfdTimetwo=0.0, cfdTimethree=0.0;
    double phaseone=0.0, phasetwo=0.0, tmaxone=0., tmaxtwo=0.;
    double timeone=0.0, timetwo=0.0, timethree=0.0;
    double qdcone=0.0, qdctwo=0.0;
    double tdiff=0.0, baselineone=0.0, baselinetwo=0.0;
    double ecaltwo=0.0;
    int maxvalone=0, maxvaltwo=0;
    double xa=0., xb=0., ya=0., yb=0.;
    double xaqdc=0., xbqdc=0., yaqdc=0., ybqdc=0.;
    double xatm=0., xbtm=0., yatm=0., ybtm=0.;
    double xposq=0., yposq=0.;
    double xpost=0., ypost=0.;

  timing->Branch("eventNum", &eventNum);
  /*timing->Branch("energyone", &energyone);
  timing->Branch("energytwo", &energytwo);
  timing->Branch("energythree", &energythree);
  timing->Branch("energycal", &ecaltwo);*/
  //timing->Branch("traceone", &traceone);
  //timing->Branch("tracetwo", &tracetwo);
  /*timing->Branch("ione", &ione);
  timing->Branch("itwo", &itwo);
  timing->Branch("hrTimeone", &hrTimeone);
  timing->Branch("hrTimetwo", &hrTimetwo);
  timing->Branch("cfdTimeone", &cfdTimeone);
  timing->Branch("cfdTimetwo", &cfdTimetwo);*/
  timing->Branch("phaseone", &phaseone);
  timing->Branch("phasetwo", &phasetwo);
  /*timing->Branch("timeone", &timeone);
  timing->Branch("timetwo", &timetwo);
  timing->Branch("tmaxone", &tmaxone);
  timing->Branch("tmaxtwo", &tmaxtwo);*/
  timing->Branch("qdcone", &qdcone);
  timing->Branch("qdctwo", &qdctwo);
  timing->Branch("tdiff", &tdiff);
  timing->Branch("xposq", &xposq);
  timing->Branch("yposq", &yposq);
  timing->Branch("xpost", &xpost);
  timing->Branch("ypost", &ypost);
  //timing->Branch("bone", &baselineone);
  //timing->Branch("btwo", &baselinetwo);

  for(int i=0;i<=numfiles;i++){
    /*if(i<10)
      fileName = prefix+"00"+i+ending;
    else if(i<100)
      fileName = prefix+"0"+i+ending;
    else*/
      fileName = prefix+i+ending;
    _file0.emplace_back(TFile::Open(fileName.c_str()));
    cout<<endl<<fileName<<endl;

    TTree *GS = (TTree*)_file0.at(i)->Get("PixTree");
    TTreeReader singe;
    singe.SetTree( GS );
    TTreeReaderArray<processor_struct::ROOTDEV> rd = {singe, "root_dev_vec_"};  //gives vector of stucture 
    cout<<"readermade"<<endl;

    newFile->cd();

    fit->SetParLimits(0,500,50000);
    fit->SetParLimits(2,5,70);
    fit->SetParLimits(3,40,200);
    fit->SetParLimits(5,1.0,1.8); //power for fit
    std::vector<unsigned int> *trace;
    while(singe.Next()){

    if(eventNum%1000==0) cout<<"\r" << eventNum<< flush;
    tone.clear(); ttwo.clear();
    traceone.clear();
    tracetwo.clear();
    xat.clear(); xbt.clear(); yat.clear(); ybt.clear();
    xatrace.clear(); xbtrace.clear();
    yatrace.clear(); ybtrace.clear();
    ione=0; itwo=0; ithree=0;
    hrTimeone=0.0; hrTimetwo=0.0; hrTimethree=0.0;
    energyone=-1000.0; energytwo=-1000.0; energythree=-1000.0;
    cfdTimeone=0.0; cfdTimetwo=0.0; cfdTimethree=0.;
    phaseone=0.0; phasetwo=0.0;
    timeone=0.0; timetwo=0.0, timethree=0.0;
    tmaxone=0.; tmaxtwo=0.; tdiff=0.0;
    ecaltwo=0.0; baselineone=0.0; baselinetwo=0.0;
    qdcone=0.0;qdctwo=0.0;
    maxvalone=0;maxvaltwo=0;
    xa=0; xb=0.; ya=0.; yb=0.;
    xposq=0.; yposq=0.;
    xpost=0.; ypost=0.;
    xaqdc=0.; xbqdc=0.; yaqdc=0.; ybqdc=0.;
    xatm=0.; xbtm=0.; yatm=0.; ybtm=0.;

    for(auto itC = rd.begin(); itC!=rd.end();itC++){
        trace = &(itC->trace);
        double energy = itC->energy;
        double time = itC->timeSansCfd;
        double highResTime = itC->highResTime;
        std::string dtype = itC->subtype.Data();
        std::string dgroup = itC->group.Data();
        int channel = itC->chanNum; 
        int module = itC->modNum;
        int maxValue = itC->maxVal;
        if(eventNum==0) cout<<trace->size()<<" "<<channel<<" "<<eventNum<<endl;

        if (channel==4){ 
            ione++;
            tone = *trace;
            double hold =0.0;
            for(int i=0;i<tone.size();i++){
              hold = (double) tone[i];
              traceone.push_back(hold);
              //if(traceone[tmaxone]<=traceone[i]) tmaxone=i;
            }
            //for(int i=0;i<traceone.size()-2;i++)traceone[i] = (traceone[i]+traceone[i+1])/2;
            energyone = energy;
            hrTimeone = highResTime;
            cfdTimeone = time;
            timeone = time;
            maxvalone = maxValue;
            //phaseone = cdfAnalyzer(traceone);
            //cfdTimeone += 2*phaseone;
            qdcone = traceAnalyzer(traceone);
            baselineone = baselineCalc(traceone).first;
        }

        else if(channel==5){ //hagrid
            itwo++;
            ttwo = *trace;
            double hold =0;
            for(int i=0;i<ttwo.size();i++){
              hold = (double) ttwo[i];
              tracetwo.push_back(hold);
              if(tracetwo[tmaxtwo]<=tracetwo[i]) tmaxtwo=i;
            }
            //for(int i=0;i<tracetwo.size()-2;i++)tracetwo[i] = (tracetwo[i]+tracetwo[i+1])/2;
            energytwo = energy;
            hrTimetwo = highResTime;
            cfdTimetwo = time;
            timetwo = time;
            maxvaltwo=maxValue;
            //phasetwo = cdfAnalyzer(tracetwo);
            //cfdTimetwo += 2*phasetwo;
            qdctwo = traceAnalyzer(tracetwo);
            baselinetwo = baselineCalc(tracetwo).first;
        }
        else if(channel==6){ //hagrid
            ithree++;
            energythree = energy;
            hrTimethree = highResTime;
            //cfdTimethree = time;
        }

        else if(channel==0){ 
          xat = *trace;
          double hold =0;
          for(int i=0;i<xat.size();i++){
              hold = (double) xat[i];
              if(hold>xatm) xatm = hold;
              xatrace.push_back(hold);
            }
          xa = energy;
          xaqdc = traceAnalyzer(xatrace);
        }
        else if(channel==1){ 
          xbt = *trace;
          double hold =0;
          for(int i=0;i<xbt.size();i++){
              hold = (double) xbt[i];
              if(hold>xbtm) xbtm = hold;
              xbtrace.push_back(hold);
            }
          xb = energy;
          xbqdc = traceAnalyzer(xbtrace);
        }
        else if(channel==2){ 
          yat = *trace;
          double hold =0;
          for(int i=0;i<yat.size();i++){
              hold = (double) yat[i];
              if(hold>yatm) yatm = hold;
              yatrace.push_back(hold);
            }
          ya = energy;
          yaqdc = traceAnalyzer(yatrace);
        }
        else if(channel==3){ 
          ybt = *trace;
          double hold =0;
          for(int i=0;i<ybt.size();i++){
              hold = (double) ybt[i];
              if(hold>ybtm) ybtm = hold;
              ybtrace.push_back(hold);
            }
          yb = energy;
          ybqdc = traceAnalyzer(ybtrace);
        }
    }
    if(ione>0 && itwo>0){
      if(1){
        phaseone = 2*pcfdAnalyzer(traceone,0.45).first;
        tmaxone = pcfdAnalyzer(traceone,0.45).second;
        cfdTimeone += phaseone;
        phasetwo = 2*pcfdAnalyzer(tracetwo,0.45).first;
        tmaxtwo = pcfdAnalyzer(tracetwo,0.45).second;
        cfdTimetwo += phasetwo;
        tdiff = cfdTimeone-phaseone-cfdTimetwo+phaseone;
        ecaltwo = 0.09369787313*energytwo-1.639424969;
        xposq = position(xaqdc,xbqdc,yaqdc,ybqdc).first;
        yposq = position(xaqdc,xbqdc,yaqdc,ybqdc).second;
        xpost = position(xatm,xbtm,yatm,ybtm).first;
        ypost = position(xatm,xbtm,yatm,ybtm).second;
      }

      //plot1->Fill(cfdTimeone-cfdTimetwo);
      //energy->Fill(energyone,energytwo);
      //phase->Fill(phaseone,phasetwo);
      newFile->cd();
      timing->Fill();
      //cout<<"Phase calculations "<<eventNum<<endl;
    }
    eventNum++;
  }
  newFile->cd();
  timing->Write("", TObject::kOverwrite);
  //plot1->Write();
  //phase->Write();
  //energy->Write();
  _file0.at(i-1)->Close();
  }
}

pair<double,double> pcfdAnalyzer(vector<double> trace,double thresh){
  int maxPos=0;
  vector<double>::iterator it;
  TGraph *fTraces = new TGraph();
  double fmax=0.0;
  double pmax=0.0;
  double baseline=0.0;
  //double thresh=0.0;
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
  
  std::pair <double, double> range((maxPos-1),(maxPos+1));
  fpol2->SetRange(range.first,range.second);
  fTraces->Fit(fpol2,"RNQSW");
  if(!(fpol2->IsValid()))return make_pair(-5555.,-5555.);
  fmax = fpol2->GetMaximum(range.first,range.second);
  pmax = trace.at(maxPos);
  baseline = baselineCalc(trace).first;
  if(pmax<3*baselineCalc(trace).second)cout<<"peak too small "<<baselineCalc(trace).second<<endl;
  //if(maxPos>125)return -7777.;
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
  return make_pair(phase,fmax-baseline);
}

double regCfd(vector<double> &trace){
  double phase=0.0;
  vector<double> invTrace(trace.size(),0.0);
  double baseline=0.0;
  baseline = baselineCalc(trace).first;
  int highbound=0, lowbound=1000;
  int maxPos=0;
  vector<double>::iterator it;
  vector<double>::iterator itMin;
  vector<double>::iterator itMax;
  TGraph *fTraces = new TGraph();

  if(trace.size()>10){
    it = max_element(trace.begin(),trace.end());
    maxPos = distance(trace.begin(),it);
  }
  else
    return -9999.; // This is about the range of the waveform

  double totMin=0;
  double totMax=0;
  double zeroPos=0;

  //if(trace[maxPos]-baseline<100)
    //return 600;

  for(int i=0;i<trace.size();i++){
    trace[i] -= baseline;
    invTrace[i]=(-0.7*trace[i+3]);
    trace[i]+=invTrace[i];
  }
  for(int i=maxPos;i>maxPos-3;i--){
    if(trace[i]<0){
      zeroPos=i;
      break;
    }
  }

  for (int iv = 0;iv<trace.size();iv++){
    fTraces->SetPoint(iv, iv, trace.at(iv));
  }
  
  double m=0.0, b=0.0;
  m = (trace.at(zeroPos+1)-trace.at(zeroPos))/1;
  b = trace.at(zeroPos+1)-m*(zeroPos+1);
  phase = -1*b/m; 

  if(phase<0.0) return -1000;

  /*fpol1->SetRange(zeroPos,zeroPos+1);
  fTraces->Fit(fpol1,"RNQSW");
  phase = fpol1->GetX(0.0,zeroPos,zeroPos+1);*/
  return phase;
  
}

double derAnalyzer(vector<double> &trace){
  int maxPos=0;
  vector<double>::iterator it;
  TGraph *fTraces = new TGraph();
  TGraph *fTraces1 = new TGraph();
  TF1 *fpol2 = new TF1("fName2","pol2",0,1);
  double fmax=0.0;
  double phase=0.0;
  double par[2]={0.0,0.0};
  vector<double> der;
  double hold=0.0;

  if(trace.size()<=10){
    return -9999.;
  }

  for (int iv = 0;iv<trace.size();iv++){
    fTraces1->SetPoint(iv, iv, trace.at(iv));
  }
  der.push_back(0.0);
  for(int i=1;i<trace.size()-1;i++){
    /*fpol1->SetRange(i-1,i+1);
    fTraces1->Fit(fpol1,"RNQSW");
    fpol1->GetParameters(&par[0]);
    der.push_back(par[1]);*/
    hold = (trace[i+1]-trace[i-1])/2.0;
    der.push_back(hold);
  }
  der.push_back(0.0);
  /*for(int i=0;i<trace.size();i++){
    //if(der[i]>der[maxPos]) maxPos=i;
    trace[i] = der[i];
  }*/
  it = max_element(der.begin(),der.end());
  maxPos = distance(der.begin(),it);
  for (int iv = 0;iv<der.size();iv++){
    fTraces->SetPoint(iv, iv, der.at(iv));
  }

  pair <double, double> range((maxPos-1),(maxPos+1));
  fpol2->SetRange(range.first,range.second);
  fTraces->Fit(fpol2,"RNQSW");
  fmax = fpol2->GetMaximum(range.first,range.second);
  phase = fpol2->GetX(fmax,range.first,range.second);
  return phase;

}

pair <double, double> baselineCalc(vector<double> trace){
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

double traceAnalyzer(vector<double> trace){
    double qdc=0.0;
    unsigned int highbound, lowbound;
    double baseline = baselineCalc(trace).first;
    int maxPos=0;
    vector<double>::iterator it;

    if(trace.size()>10){
    it = max_element(trace.begin(),trace.end());
    maxPos = distance(trace.begin(),it);
    }
    else if(maxPos<180){
      return -4444.;
    }
    else{
      return -9999.;
    }

    highbound = maxPos+40;
    lowbound = maxPos-20;

    for(int i=lowbound;i<highbound;i++)
      qdc += 0.5 * double(abs(trace[i-1] - baseline) + abs(trace[i] - baseline));

    return qdc;
}

double fittingAnalyzer(vector<double> trace){
  int maxLoc=0,dMax=0;
  double phase=0.0;
  double baseline=baselineCalc(trace).first;
  TGraph *fTraces = new TGraph();
  double fmax=0.0;
  double der=0.0;
  
  double thresh=0.4;
  std::pair <double,double> points(0,0);

  if(trace.size()<10) return -9999.;

  //need to determine fitting range..want to start around beginning of rise
  for(int i=0;i<trace.size();i++)
    if(trace[i]>trace[maxLoc])maxLoc=i;

  for(int i=0;i<trace.size()-1;i++){
    //this is taking the derivative to determine the max risetime, which should be on the initial kink up
    if(trace[i+1]-trace[i]>der && i<=maxLoc)dMax=i; 
    der = trace[dMax+1]-trace[dMax];
    fTraces->SetPoint(i, i, trace.at(i));
  }
  fTraces->SetPoint(trace.size()-1, trace.size()-1, trace.at(trace.size()-1));
  
  //general ranges for parameters, to help the fit actually work.
  fit->SetParameter(1,dMax);
  fit->SetParLimits(1,dMax-5,dMax+10);
  fit->FixParameter(4,baseline);

  std::pair <double, double> range((dMax-2),dMax+50);
  fit->SetRange(range.first,range.second);
  fTraces->Fit(fit,"RNSQW");
  fmax = fit->GetMaximum(range.first,range.second);
  //phase = fit->GetX(fmax,range.first,range.second);
  //phase = fit->GetParameter(1);
  phase = fit->GetParameter(1);
  //return phase;
  thresh = (fmax-baseline)*0.5+baseline; //here 0.5 is the constant fraction
  for(int i=maxLoc;i>0;i--){
    if((trace.at(i)>=thresh) && (trace.at(i-1)<thresh)){
      points.first = i-1;
      points.second = i;
    }
  }
  fpol1->SetRange(points.first,points.second);
  fTraces->Fit(fpol1,"RNSQW");
  phase = fpol1->GetX(thresh,points.first,points.second);

  return phase;
}

double gFit(vector<double> trace){
  int maxPos=0;
  vector<double>::iterator it;
  TGraph *fTraces = new TGraph();
  double fmax=0.0;
  double pmax=0.0;
  double baseline=0.0;
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
  
  std::pair <double, double> range((maxPos-2),(maxPos+2));
  gauss->SetRange(range.first,range.second);
  fTraces->Fit(gauss,"RNQSW");
  if(!(gauss->IsValid()))return -5555.;
  fmax = gauss->GetMaximum(range.first,range.second);
  phase = gauss->GetX(fmax,range.first,range.second);
  return phase;
}

vector<double> reducer(vector<double> trace){
  int n;
  n = trace.size();
  vector<double> t2;
  for(int i=0;i<n;i++){
    if(i%2==0) t2.push_back(0.5*(trace[i]+trace[i+1]));
  }
  return t2;
}

double centroid(vector<double> trace){

  double cent=0.0, ampSum=0.0, time=0.0;
  double baseline = baselineCalc(trace).first;
  double stddev = baselineCalc(trace).second;
  int maxLoc=0;
  pair <int, int> wave = std::make_pair(0,0);

  for(int i=0;i<trace.size();i++){
    trace[i]-=baseline;
    if(trace[i]>trace[maxLoc]) maxLoc=i;
  }
  for(int i=maxLoc-4;i<trace.size();i++){
    if(trace[i]*i>cent)cent=trace[i]*i;
    ampSum+=trace[i];
  }

  return cent/ampSum;
}

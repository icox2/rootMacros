#include "TSystem.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
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

double cfdAnalyzer(vector<double> trace);
double regCfd(vector<double> &trace);
double derAnalyzer(vector<double> &trace);
pair <double, double> baselineCalc(vector<double> trace);
double traceAnalyzer(vector<double> trace);
double fittingAnalyzer(vector<double> trace, bool pow);
TF1 *fit = new TF1("fit","[0]*exp(([1]-x)/[2])*(1-exp(-1*pow(x-[1],2)/[3]))+[4]",0,300);
TF1 *fpol1 = new TF1("fName1","pol1",0,300);
TF1 *fpol3 = new TF1("fName3","pol3",0,1);

void timingAnalyzer(){
  TFile *_file0 = TFile::Open("yso_hag2_1400V_DD.root");
  TTree *GS = (TTree*)_file0->Get("PixTree");
  TTreeReader singe;
  cout<<"readermade"<<endl;
  singe.SetTree( GS );
  TTreeReaderArray<processor_struct::ROOTDEV> rd = {singe, "root_dev_vec_"};  //gives vector of stucture 

    TH1D *plot1=new TH1D("plot1","plot1",10000.,-2000.,2000.);
    TH2D *energy=new TH2D("energy","energy",5000.,0.,10000.,5000.,0.,10000.);
    TH2D *phase=new TH2D("phase","phase",1000.,80.,105.,1000.,80.,105.);

    //variables
    int eventNum = 0;
    double energyone=0.0, energytwo=0.0;
    vector<double> traceone, tracetwo;
    vector<unsigned int> tone, ttwo;
    int ione=0, itwo=0, tiny=0;
    double hrTimeone=0.0, hrTimetwo=0.0;
    double cfdTimeone=0.0, cfdTimetwo=0.0;
    double phaseone=0.0, phasetwo=0.0;
    double timeone=0.0, timetwo=0.0;
    double qdcone=0.0, qdctwo=0.0;
    double tdiff=0.0;
    double ecaltwo=0.0;
    int maxvalone=0, maxvaltwo=0;
    double xa=0., xb=0., ya=0., yb=0.;
    double xpos=0., ypos=0.;

  TFile *newFile = new TFile("yso_hag2_1400V.root","RECREATE");
  TTree *timing = new TTree("timing","tree filled with traces and energies etc.");

  timing->Branch("eventNum", &eventNum);
  timing->Branch("energyone", &energyone);
  timing->Branch("energytwo", &energytwo);
  timing->Branch("energycal", &ecaltwo);
  timing->Branch("traceone", &traceone);
  timing->Branch("tracetwo", &tracetwo);
  timing->Branch("ione", &ione);
  timing->Branch("itwo", &itwo);
  timing->Branch("hrTimeone", &hrTimeone);
  timing->Branch("hrTimetwo", &hrTimetwo);
  timing->Branch("cfdTimeone", &cfdTimeone);
  timing->Branch("cfdTimetwo", &cfdTimetwo);
  timing->Branch("phaseone", &phaseone);
  timing->Branch("phasetwo", &phasetwo);
  timing->Branch("timeone", &timeone);
  timing->Branch("timetwo", &timetwo);
  timing->Branch("tiny", &tiny);
  timing->Branch("qdcone", &qdcone);
  timing->Branch("qdctwo", &qdctwo);
  timing->Branch("tdiff", &tdiff);
  timing->Branch("xpos", &xpos);
  timing->Branch("ypos", &ypos);

  fit->SetParLimits(0,500,50000);
  fit->SetParLimits(2,5,60);
  fit->SetParLimits(3,1,60);

  std::vector<unsigned int> *trace;
  while(singe.Next()){

    if(eventNum%1000==0) cout<<"\r" << eventNum<< flush;
    tone.clear(); ttwo.clear();
    traceone.clear();
    tracetwo.clear();
    ione=0; itwo=0;
    hrTimeone=0.0; hrTimetwo=0.0;
    energyone=-1000.0; energytwo=-1000.0;
    cfdTimeone=0.0; cfdTimetwo=0.0;
    phaseone=0.0; phasetwo=0.0;
    timeone=0.0; timetwo=0.0;
    tiny=0; tdiff=0.0;
    ecaltwo=0.0;
    qdcone=0.0;qdctwo=0.0;
    maxvalone=0;maxvaltwo=0;
    xa=0; xb=0.; ya=0.; yb=0.;
    xpos=0.; ypos=0.;

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
        if(eventNum==0) cout<<trace->size()<<" "<<channel<<endl;

        if (channel==0){ //12b250
            ione++;
            tone = *trace;
            double hold =0.0;
            for(int i=0;i<tone.size();i++){
              hold = (double) tone[i];
              traceone.push_back(hold);
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
        }

        else if(channel==11){ //14b500
            itwo++;
            ttwo = *trace;
            double hold =0;
            for(int i=0;i<ttwo.size();i++){
              hold = (double) ttwo[i];
              tracetwo.push_back(hold);
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
        }

        else if(channel==1) xa = energy;
        else if(channel==2) xb = energy;
        else if(channel==3) ya = energy;
        else if(channel==4) yb = energy;
    }
    if(ione>0 && itwo>0){
      if(1){
        phaseone = cfdAnalyzer(traceone);
        cfdTimeone += 2*phaseone;
        phasetwo = cfdAnalyzer(tracetwo);
        cfdTimetwo += 2*phasetwo;
        tdiff = cfdTimeone-phaseone-cfdTimetwo+phaseone;
        ecaltwo = 0.09369787313*energytwo-1.639424969;
        xpos = position(xa,xb,ya,yb).first;
        ypos = position(xa,xb,ya,yb).second;
      }

      plot1->Fill(cfdTimeone-cfdTimetwo);
      energy->Fill(energyone,energytwo);
      phase->Fill(phaseone,phasetwo);
      timing->Fill();
    }
    eventNum++;
  }
  timing->Write();
  plot1->Write();
  phase->Write();
  energy->Write();
}

double cfdAnalyzer(vector<double> trace){
  int maxPos=0;
  vector<double>::iterator it;
  TGraph *fTraces = new TGraph();
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
  
  std::pair <double, double> range((maxPos-2),(maxPos+1));
  fpol3->SetRange(range.first,range.second);
  fTraces->Fit(fpol3,"RNQSW");
  if(!(fpol3->IsValid()))return -5555.;
  fmax = fpol3->GetMaximum(range.first,range.second);
  pmax = trace.at(maxPos);
  baseline = baselineCalc(trace).first;
  if(pmax<3*baselineCalc(trace).second)cout<<"peak too small "<<baselineCalc(trace).second<<endl;
  if(maxPos>125)return -7777.;
  thresh = (fmax-baseline)*0.3+baseline; //here 0.3 is the Constant Fraction
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

double regCfd(vector<double> &trace){
  double phase=0.0;
  vector<double> invTrace;
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
    return -9999.;
  highbound = maxPos+60;
  lowbound = maxPos-20; // This is about the range of the waveform

  double totMin=0;
  double totMax=0;
  double zeroPos=0;

  if(trace[maxPos]-baseline<100)
    return 600;

  for(int i=0;i<trace.size();i++){
    trace[i] -= baseline;

    if(i>=lowbound+2 && i<highbound+2){
      invTrace.push_back(-0.5*trace[i+2]);
    }
    else 
      invTrace.push_back(0.0);

    trace[i]+=invTrace[i];

    if(trace[i]<trace[totMin]) totMin=i;
    if(trace[i]>trace[totMax]) totMax=i;
  }
  for(int i=totMin;i<totMax;i++){
    if(abs(trace[i])<abs(trace[zeroPos])) zeroPos=i;
  }

  for (int iv = 0;iv<trace.size();iv++){
    fTraces->SetPoint(iv, iv, trace.at(iv));
  }
  
  /*double m=0.0, b=0.0;
  cout<<zeroPos<<endl;
  m = (trace.at(zeroPos+1)-trace.at(zeroPos-1))/2;
  b = trace.at(zeroPos+1)-m*(zeroPos+1);
  phase = -1*b/m; */

  fpol1->SetRange(zeroPos-1,zeroPos+1);
  fTraces->Fit(fpol1,"RNQSW");
  phase = fpol1->GetX(0,zeroPos-1,zeroPos+1);
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

/*double calculateP2(const short &x0, unsigned short *y, double *p){
	double x1[3], x2[3];
	for(size_t i = 0; i < 3; i++){
		x1[i] = (x0+i);
		x2[i] = std::pow(x0+i, 2);
	}

	double denom = (x1[1]*x2[2]-x2[1]*x1[2]) - x1[0]*(x2[2]-x2[1]*1) + x2[0]*(x1[2]-x1[1]*1);

	p[0] = (y[0]*(x1[1]*x2[2]-x2[1]*x1[2]) - x1[0]*(y[1]*x2[2]-x2[1]*y[2]) + x2[0]*(y[1]*x1[2]-x1[1]*y[2]))/denom;
	p[1] = ((y[1]*x2[2]-x2[1]*y[2]) - y[0]*(x2[2]-x2[1]*1) + x2[0]*(y[2]-y[1]*1))/denom;
	p[2] = ((x1[1]*y[2]-y[1]*x1[2]) - x1[0]*(y[2]-y[1]*1) + y[0]*(x1[2]-x1[1]*1))/denom;
	
	// Calculate the maximum of the polynomial.
	return (p[0] - p[1]*p[1]/(4*p[2]));
}
*/

double fittingAnalyzer(vector<double> trace, bool pow){
  int maxLoc=0,dMax=0;
  double phase=0.0;
  double baseline=baselineCalc(trace).first;
  TGraph *fTraces = new TGraph();
  double fmax=0.0;
  double der=0.0;
  
  double thresh=0.0;
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
  fit->SetParLimits(1,dMax-8,dMax+20);
  fit->SetParameter(4,baseline);

  std::pair <double, double> range((dMax-4),trace.size());
  fit->SetRange(range.first,range.second);
  fTraces->Fit(fit,"RNSQW");
  fmax = fit->GetMaximum(range.first,range.second);
  //phase = fit->GetX(fmax,range.first,range.second);
  //phase = fit->GetParameter(1);

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
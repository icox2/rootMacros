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
//#include <ProcessorRootStruc.hpp>
#include <PaassRootStruct.hpp>
#include <cmath>
#include <tuple>
#include "position.h"

void createPileup(vector<double> lasttrace, vector<double> &trace, int phase);
void normalizeTrace(vector<double> &trace);
std::tuple<double,double,double> pcfdAnalyzer(vector<double> trace, double frac);
double regCfd(vector<double> &trace);
double derAnalyzer(vector<double> &trace);
pair <double, double> baselineCalc(vector<double> trace);
double traceAnalyzer(vector<double> trace);
double fittingAnalyzer(vector<double> trace);
double gFit(vector<double> trace);
vector<double> reducer(vector<double> trace);
double centroid(vector<double> trace);
TF1 *fit = new TF1("fit","(x>[1])*[0]*exp(([1]-x)/[2])*(1-exp(-1*pow(x-[1],[5])/[3]))+[4]",0,300);
TF1 *fpol1 = new TF1("fName1","pol1",0,300);
TF1 *fpol3 = new TF1("fName3","pol3",0,1);
TF1 *fpol2 = new TF1("fName2","pol2",0,1);
TF1 *gauss = new TF1("gauss","gaus",0,1);
TRandom2 *randGenerator = new TRandom2(); 

void getTraces(){
  TFile *_file0 = TFile::Open("ysoSingle_DD.root");
  TTree *GS = (TTree*)_file0->Get("PixTree");
  TTreeReader singe;
  cout<<"readermade"<<endl;
  singe.SetTree( GS );
  //TTreeReaderArray<processor_struct::ROOTDEV> rd = {singe, "root_dev_vec_"};  //gives vector of stucture 
  TTreeReaderArray<processor_struct::ROOTDEV> rd = {singe, "rootdev_vec_"};  //gives vector of stucture 

  double percentPileup = 0.4;
  double phaseHigh = 100.;
  double phaseLow = 2.;

  //variables
  int eventNum = 0;

  double energyone=0.0;
  vector<double> traceone, lasttrace;
  vector<unsigned int> tone;
  int ione=0;
  double hrTimeone=0.0;
  double cfdTimeone=0.0;
  double phaseone=0.0, tmaxone=0.;
  double timeone=0.0;
  double qdcone=0.;
  double baselineone=0.0;
  int maxvalone=0;
  double phase=0.;
  bool pileup=0;

  TFile *newFile = new TFile("ysoTracesNoPileup.root","RECREATE");
  TTree *timing = new TTree("timing","tree filled with traces and energies etc.");

//  timing->Branch("eventNum", &eventNum);
  timing->Branch("trace", &traceone);
  timing->Branch("phase", &phase);
  timing->Branch("pileup", &pileup);

  std::vector<unsigned int> *trace;
  while(singe.Next()){

    if(eventNum%1000==0) cout<<"\r" << eventNum<< flush;
    tone.clear();
    traceone.clear();
//    lasttrace.clear();
    ione=0;
    hrTimeone=0.0;
    energyone= -1000.0;
    cfdTimeone=0.0;
    phaseone=0.0;
    timeone=0.0;
    tmaxone=0.0;
    baselineone=0.0;
    qdcone=0.0;
    maxvalone=0;
    phase = -99.;
    pileup = 0;

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

        if (channel==0 && module==1){ 
            ione++;
            tone = *trace;
            double hold =0.0;
            for(int i=0;i<tone.size();i++){
              hold = (double) tone[i];
              traceone.push_back(hold);
              if(traceone[tmaxone]<=traceone[i]) tmaxone=i;
            }
            //for(int i=0;i<traceone.size()-2;i++)traceone[i] = (traceone[i]+traceone[i+1])/2;
            energyone = energy;
            hrTimeone = highResTime;
            cfdTimeone = time;
            timeone = time;
            maxvalone = traceone[tmaxone];
            //phaseone = cdfAnalyzer(traceone);
            //cfdTimeone += 2*phaseone;
            qdcone = traceAnalyzer(traceone);
            auto baselineVals = baselineCalc(traceone);
            baselineone = baselineVals.first;
            //sigmaone = baselineVals.second;
        }
    }
    if(maxvalone>4000){
      normalizeTrace(traceone);
      //vector<double> beforepile = traceone;
      //if(randGenerator->Rndm()<percentPileup && lasttrace.size()>0){
      //  pileup = 1;
      //  phase = randGenerator->Rndm()*(phaseHigh-phaseLow)+phaseLow;
      //  createPileup(lasttrace,traceone,phase);
      //}
      //lasttrace = beforepile;
      //normalizeTrace(traceone);
      timing->Fill();
    }
    eventNum++;
  }
  timing->Write();
  cout<<endl;
}

////////// Function Definitions ////////////
void createPileup(vector<double> lasttrace, vector<double> &trace, int phase){
  double scale = randGenerator->Rndm()+0.5; // Allows the lasttrace to be scaled by 0.5x up to 1.5x
  for(int i=phase;i<trace.size();i++){
    trace.at(i) += lasttrace.at(i-phase);
  }
  return;
}

void normalizeTrace(vector<double> &trace){
  double baseline = baselineCalc(trace).first;
  auto it = max_element(std::begin(trace), std::end(trace));
  double maxVal = *it;
  for(int i=0;i<trace.size();i++){
    trace.at(i) = (trace.at(i)-baseline)/maxVal;
  }
  return;
}


std::tuple<double,double,double> pcfdAnalyzer(vector<double> trace, double frac){
  int maxPos=0;
  vector<double>::iterator it;
  TGraph *fTraces = new TGraph();
  double fmax=0.0;
  double pmax=0.0;
  double baseline=834.;
  double thresh=0.0;
  double phase=0.0;
  std::pair <double,double> points(0.,0.);

  if(trace.size()>10){
    it = max_element(trace.begin(),trace.end());
    maxPos = distance(trace.begin(),it);
  }
  else{
    return make_tuple(-9999.,-9999.,-9999.);
  }
  for (int iv = 0;iv<trace.size();iv++){
    fTraces->SetPoint(iv, iv, trace.at(iv));
  }
  
  std::pair <double, double> range((maxPos-2),(maxPos+1));
  fpol3->SetRange(range.first,range.second);
  fTraces->Fit(fpol3,"RNQSW");
  if(!(fpol3->IsValid()))return make_tuple(-5555.,-5555.,-5555.);
  fmax = fpol3->GetMaximum(range.first,range.second);
  pmax = trace.at(maxPos);
  //if(frac>0.11)
  baseline = baselineCalc(trace).first;
  if(pmax<3*baselineCalc(trace).second)cout<<"peak too small "<<baselineCalc(trace).second<<endl;
  //if(maxPos>125)return -7777.;
  thresh = (fmax-baseline)*frac+baseline; //here frac is the Constant Fraction
  for(int i=maxPos;i>0;i--){
    if((trace.at(i)>=thresh) && (trace.at(i-1)<thresh)){
      points.first = i-1;
      points.second = i; //setting points of trace around threshold 
    }
  }
  fpol1->SetRange(points.first,points.second);
  fTraces->Fit(fpol1,"RNQSW");
  phase = fpol1->GetX(thresh,points.first,points.second);
  return make_tuple(phase,fmax-baseline,fpol1->GetParameter(1));
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
  //TF1 *fpol2 = new TF1("fName2","pol2",0,1);
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
    const int numBins = 20;
    double baseSum = 0.;
    for(int j=0;j<numBins;j++){
        baseSum += trace.at(j);
    }
    double  baseline = baseSum/(double(numBins));
    
    //calculating the standard dev
    double stddev = 0.0;
    for(int j=0;j<numBins;j++){
        stddev += pow(trace.at(j)-baseline,2);
    }
    stddev = sqrt(stddev/(double(numBins)));
    
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
  if(trace.size()<10) return -9999.;
  double baseline=baselineCalc(trace).first;
  //TGraph *fTraces = new TGraph();
  TH1D *hist = new TH1D("hist","trace",250,0.,250.);
  double fmax=0.0;
  double der=0.0;
  
  double thresh=0.4;
  std::pair <double,double> points(0,0);


  //need to determine fitting range..want to start around beginning of rise
  for(int i=0;i<trace.size();i++){
    if(trace.at(i)>trace.at(maxLoc)) maxLoc=i;
    //fTraces->SetPoint(i, i, trace.at(i));
    if(trace[i+1]-trace[i]>der && i<=maxLoc && i<trace.size()-1){
      dMax=i; 
      der = trace[dMax+1]-trace[dMax];
    }
    for(int k=0;k<(int)trace.at(i);k++) hist->Fill(double(i));
  }

  /*for(int i=0;i<trace.size()-1;i++){
    //this is taking the derivative to determine the max risetime, which should be on the initial kink up
    if(trace[i+1]-trace[i]>der && i<=maxLoc)dMax=i; 
    der = trace[dMax+1]-trace[dMax];
  }*/
  for(int j=0;j<trace.size()-1;j++){
    double hold=0., temp=0.;
    temp = trace.at(j+1)-trace.at(j);
    if(temp>0) hold = temp;
    else  hold = -0.2*temp;
    //fTraces->SetPoint(j, 2*j, hold);
    double weight = hold/der+0.001;
    //fTraces->SetPoint(j, 2*j,weight);
    for(int k=0;k<(int)trace.at(j);k++) hist->Fill(double(j));
    hist->SetBinError(j,1/weight);
  }
  
  //general ranges for parameters, to help the fit actually work.
  //fit->SetParameter(1,dMax);
  fit->SetParLimits(0,500,50000);
  fit->SetParLimits(1,dMax-12,dMax);
  fit->SetParLimits(2,1,100);
  fit->SetParLimits(3,1,2000);
  fit->FixParameter(4,baseline);
  fit->FixParameter(5,3);

  std::pair <double, double> range(dMax-20,maxLoc+10);
  fit->SetRange(range.first,range.second);
  hist->Fit(fit,"RNSQW");
  /*fmax = fit->GetMaximum(range.first,range.second);
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
*/
  delete hist;
  return fit->GetParameter(1);
  //return phase;
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

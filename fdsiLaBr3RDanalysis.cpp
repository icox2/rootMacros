#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TCanvas.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <PaassRootStruct.hpp>
#include <vector>

std::tuple<double,double,double> pcfdAnalyzer(vector<double> trace, double frac);
pair <double, double> baselineCalc(vector<double> trace);
TF1 *fpol1 = new TF1("fName1","pol1",0,300);
TF1 *fpol3 = new TF1("fName3","pol3",0,1);

void fdsiLaBr3RDanalysis(){
  const int numHag = 15;
  vector<double> hagEnergy(numHag,0);
  vector<double> hagTime(numHag,0);
  double xPosHG = 0., yPosHG = 0.;
  int eventNum=0;
  double hgDynEnergy = 0.;
  double hgDynTime = 0., hgDynHRTime=0., hgDynPcfd=0.;
  vector<double> traceone;
  double dynBase = 0., dynTMax=0.;
  double dynStd = 0.;

  TFile *newFile = new TFile("e21069co60RDtest.root","RECREATE");
  TTree *newTree = new TTree("newTree","tree filled with traces and energies etc.");

  newTree->Branch("eventNum", &eventNum);
  newTree->Branch("high_gain_xPos", &xPosHG);
  newTree->Branch("high_gain_yPos", &yPosHG);
  newTree->Branch("high_gain_dynE",&hgDynEnergy);
  newTree->Branch("high_gain_dynT",&hgDynTime);
  newTree->Branch("high_gain_dynB",&dynBase);
  newTree->Branch("high_gain_dynS",&dynStd);
  newTree->Branch("high_gain_dynM",&dynTMax);
  newTree->Branch("high_gain_dynPcfdT",&hgDynPcfd);
  newTree->Branch("high_gain_dynHighResT",&hgDynHRTime);
	newTree->Branch("hagE",&hagEnergy);
	newTree->Branch("hagT",&hagTime);

  for(int iter=0;iter<4;iter++){
    TString fname = TString("/SCRATCH/DScratch7/e21069rootfile/utkscan/crate2/crate2_RD_run_0021-0"+to_string(iter)+"_DD.root");
    TFile *_file0 = TFile::Open(fname);
    TTree *GS = (TTree*)_file0->Get("PixTree");
    TTreeReader singe;
    cout<<endl<<"readermade for file "<<iter<<endl;
    singe.SetTree( GS );
    TTreeReaderArray<processor_struct::ROOTDEV> rd = {singe, "rootdev_vec_"};  //gives vector of stucture 

    eventNum=0;
    newFile->cd(); 

    std::vector<unsigned int> *trace;
    while(singe.Next()){	
	  	if(eventNum%10000==0) cout<<"\r" << eventNum<< flush;
	  	eventNum++;
	  	//clearing the previous values from last event
	  	xPosHG = 0.; yPosHG=0.;
	  	for(int i=0;i<numHag;i++){
	  		hagEnergy[i]=-990.;
	  		hagTime[i]=-990.;
      }
      hgDynEnergy = -990;
      hgDynTime = -990;
      hgDynHRTime = -990;
      dynBase=0.; dynStd=0.;
      int aHag = 0;
      //trace->clear();
      //traceone.clear();

	  	for(auto itC = rd.begin(); itC!=rd.end();itC++){
	  		//trace = &(itC->trace);
	  		double energy = itC->energy;
	  		double time = itC->time;
	  		double highResTime = itC->highResTime;
	  		std::string dtype = itC->subtype.Data();
	  		int detNum = itC->detNum;
        bool dsat = itC->saturation;
        double base = itC->baseline;
        double max = itC->maxVal;
        double stdbase = itC->stdBaseline;
        // Saturated events are needed in order to properly reconstruct the position 
        //if(dsat) continue;

	  		int channel = itC->chanNum;  //used to try and gate on channel number rather than group/subtype
	  		if(dtype=="dynode_high"){
          hgDynEnergy = energy;
          hgDynTime = time;
          hgDynHRTime = highResTime;
          dynBase = base;
          dynStd = stdbase;
          dynTMax = max;
          //vector<unsigned int> tone = *trace;
          //for(int i=0;i<tone.size();i++){
          //  traceone.push_back((double)tone[i]);
          //}
          //std::tuple<double,double,double> result;
          //if(!traceone.empty()){
          //  result = pcfdAnalyzer(traceone,0.45);
          //  hgDynPcfd = time + 4*(std::get<0>(result));
          //}
	  		}
        else if (dtype == "smallhag"){
          hagEnergy[detNum] = energy;
          //hagTime[detNum] = highResTime;
          hagTime[detNum] = time;
          aHag++;
	  	  }
	  	}
      if(aHag>0)
	      newTree->Fill();
	  	eventNum++;
    } 
    _file0->TFile::Close();
  }
  newTree->Write();
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
  if(range.first==range.second){
    cout<<"Error in range: "<<range.first<<" "<<range.second<<endl;
    return make_tuple(-880.,-880.,-880.);
  }
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
  if(points.first==points.second){
    //cout<<"Error in points: "<<points.first<<" "<<points.second<<endl;
    return make_tuple(-770.,-770.,-770.);
  }
  fTraces->Fit(fpol1,"RNQSW");
  phase = fpol1->GetX(thresh,points.first,points.second);
  return make_tuple(phase,fmax-baseline,fpol1->GetParameter(1));
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

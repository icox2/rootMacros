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
//#include <ProcessorRootStruc.hpp>
#include <PaassRootStruct.hpp>
#include "position.h"
#include "pcfdAnalyzer.h"
//#include "pileupAnalysis.h"

double QDCcalculator(vector<unsigned int> trace, unsigned int lowBound, unsigned int highBound);
pair <int, int> boundsCalc(vector<unsigned int> trace, int maxPos);
pair <double, double> baselineCalc(vector<unsigned int> trace);
pair<double, double> position(double xa, double xb, double ya, double yb);

void pspmtRDanalysis(/*const char* filename=""*/){
  //if(sizeof(filename)==0){ 
  //  TString fname = "impTestMarchRD5_DD.root";
  //  filename = fname.Data();
  //}
  TFile *_file0 = TFile::Open("impTestMarchRD5_DD.root");
  TTree *GS = (TTree*)_file0->Get("PixTree");
  TTreeReader singe;
  cout<<"readermade"<<endl;
  singe.SetTree( GS );
  TTreeReaderArray<processor_struct::ROOTDEV> rd = {singe, "rootdev_vec_"};  //gives vector of stucture 

 
  double xal=0.0, xbl=0.0, yal=0.0, ybl=0.0;
  double xposl=0.0, yposl=0.0;  //low gain
  double xah=0.0, xbh=0., yah=0., ybh=0.;
  double xposh=0., yposh=0.;  //high gain
  double dynTime=0.;
  double hagTimean = 0.0, hagTimed = 0.0;
  double dynodeElow = 0., dynodeEhigh=0., hagEdy = 0., hagEan;
  double qdcShortAl = 0.0, qdcLongAl = 0.0;
  double qdcShortAh = 0.0, qdcLongAh = 0.0;
  double qdcShortH = 0.0, qdcLongH = 0.0;
  double qdcShortL = 0.0, qdcLongL = 0.0;
  double bD =0.0, bxa =0.0, bxb=0.0, bya=0.0, byb=0.0;
  double anodeRatio[4];
  double dphase=0.0, xaphase=0.0, xbphase=0.0, yaphase=0.0, ybphase=0.0;
  double xatime=0.0, xbtime=0.0, yatime=0.0, ybtime=0.0;
  int maxLoc = 0, iDyn = 0, ixa=0,ixb=0,iya=0,iyb=0, minimumLoc=0;
  pair<int,int> Bound;
  vector<unsigned int> dynodeTraceh, dynodeTracel, hagridTrace, xaTracel, xbTracel, yaTracel, ybTracel, anodeSum;
  vector<unsigned int> xaTraceh, xbTraceh, yaTraceh, ybTraceh;
  vector<unsigned int> xaqdcSumsVec,xbqdcSumsVec,yaqdcSumsVec,ybqdcSumsVec;
  int eventNum=0, ihaga=0, ihagd=0, undershoot=0, noCorr=0, numPile;    
  double traceMax =0.0, stddev = 0.0, traceMin = 0.0;
  double qdcS[8];//first 4 are low gain, second 4 are high gain
  double qdcL[8];
  unsigned int tmax[8];
  double qdcSum[8];
  double scale = 1., baseline = 0.0;
  double hrTimeha = 0.0, hrTimehd = 0.0, hrTimedh = 0.0;
  vector<double> derxa, derxb, derya, deryb, derdh, derdl, der2dl, der2dh;
	
  TFile *newFile = new TFile("impTestMarchRD5.root","RECREATE");
  TTree *newTree = new TTree("newTree","tree filled with traces and energies etc.");
  
  newTree->Branch("eventNum", &eventNum);
  newTree->Branch("dynodeTraceh", &dynodeTraceh);
  newTree->Branch("dynodeTracel", &dynodeTracel);
  newTree->Branch("dynodeEnergyh", &dynodeEhigh);
  newTree->Branch("dynodeEnergyl", &dynodeElow);
  newTree->Branch("hagEnergydyn", &hagEdy);
  newTree->Branch("hagEnergyan", &hagEan);
  newTree->Branch("hagTrace", &hagridTrace);
  newTree->Branch("shortQDCl", &qdcShortAl);
  newTree->Branch("longQDCl",&qdcLongAl);
  newTree->Branch("shortQDCh", &qdcShortAh);
  newTree->Branch("longQDCh",&qdcLongAh);
  newTree->Branch("shortQDCdynh", &qdcShortH);
  newTree->Branch("longQDCdynh",&qdcLongH);
  newTree->Branch("shortQDCdynl", &qdcShortL);
  newTree->Branch("longQDCdynl",&qdcLongL);
  newTree->Branch("xposl", &xposl);
  newTree->Branch("yposl", &yposl);
  newTree->Branch("xposh", &xposh);
  newTree->Branch("yposh", &yposh);
  newTree->Branch("iDyn", &iDyn);
  newTree->Branch("ixa", &ixa);
  newTree->Branch("ixb", &ixb);
  newTree->Branch("iya", &iya);
  newTree->Branch("iyb", &iyb);
  newTree->Branch("ihaga", &ihaga);
  newTree->Branch("ihagd", &ihagd);
  newTree->Branch("hagTimean", &hagTimean);
  newTree->Branch("hagTimed", &hagTimed);
  newTree->Branch("dynTime", &dynTime);
  newTree->Branch("tmax", &traceMax);
  newTree->Branch("tmin", &traceMin);
  newTree->Branch("stddev", &stddev);
  newTree->Branch("baseline", &baseline);
  newTree->Branch("anodeRatio", anodeRatio);
  newTree->Branch("qdcxa", &qdcSum[4]);
  newTree->Branch("qdcxb", &qdcSum[5]);
  newTree->Branch("qdcya", &qdcSum[6]);
  newTree->Branch("qdcyb", &qdcSum[7]);
  //newTree->Branch("qdcSum", &qdcSum[4]+&qdcSum[5]+&qdcSum[6]+&qdcSum[7])
  newTree->Branch("undershootcount", &undershoot);
  newTree->Branch("xaTrace", &xaTraceh);
  newTree->Branch("xbTrace", &xbTraceh);
  newTree->Branch("yaTrace", &yaTraceh);
  newTree->Branch("ybTrace", &ybTraceh);
  newTree->Branch("xaqdcsumsvec", &xaqdcSumsVec);
  newTree->Branch("xbqdcsumsvec", &xbqdcSumsVec);
  newTree->Branch("yaqdcsumsvec", &yaqdcSumsVec);
  newTree->Branch("ybqdcsumsvec", &ybqdcSumsVec);
  newTree->Branch("xatime", &xaphase);
  newTree->Branch("xbtime", &xbphase);
  newTree->Branch("yatime", &yaphase);
  newTree->Branch("ybtime", &ybphase);
  newTree->Branch("dPhaseh", &dphase);
  newTree->Branch("xaEnergy", &xah);
  newTree->Branch("xbEnergy", &xbh);
  newTree->Branch("yaEnergy", &yah);
  newTree->Branch("ybEnergy", &ybh);
  newTree->Branch("highResTimeha", &hrTimeha);
  newTree->Branch("highResTimedh", &hrTimedh);
  newTree->Branch("highResTimedy", &hrTimehd);
  newTree->Branch("bDyn", &bD);
  newTree->Branch("bxa", &bxa);
  newTree->Branch("bxb", &bxb);
  newTree->Branch("bya", &bya);
  newTree->Branch("byb", &byb);
  newTree->Branch("derh", &derdh);
  newTree->Branch("derl", &derdl);
  newTree->Branch("der2l", &der2dl);
  newTree->Branch("der2h", &der2dh);
  newTree->Branch("pile", &numPile);

  std::vector<unsigned> *trace;
  while(singe.Next()){
	//clearing the previous values from last event
	iDyn = 0;
	ixa = 0;
	ixb = 0;
	iya = 0;
	iyb = 0;
  dphase=0.0; xaphase=0.0; xbphase=0.0; yaphase=0.0; ybphase=0.0;
  xatime=0.0; xbtime=0.0; yatime=0.0; ybtime=0.0;
	ihagd = 0; ihaga = 0; undershoot=0;
	qdcShortAl = -999.; qdcLongAl = -999.;
	qdcShortH = -999.; qdcLongH = -999.;
	qdcShortL = -999.; qdcLongL = -999.;
	hagEdy = -999.; hagEan = -999.;
	dynodeEhigh = -999.; dynodeElow= -999.;
  bD =0.0; bxa=0.0; bxb=0.0; bya=0.0; byb=0.0;
  baseline = -999.; numPile=0;
	xposl = 0.0;
	yposl = 0.0;
  xal=0.0; xbl=0.0; yal=0.0; ybl=0.0;
  xah=0.0; xbh=0.0; yah=0.0; ybh=0.0;
  stddev = 0.0; traceMax = 0.0; traceMin = 0.0;
	for(int i=0;i<4;i++){
		anodeRatio[i]=-999.;
    qdcSum[i]=-999.;
    qdcSum[i+4]=-999.;
		scale=1.;
		qdcS[2*i] = -999.;
		qdcL[2*i] = -999.;
		qdcS[2*i+1] = -999.;
		qdcL[2*i+1] = -999.;
		tmax[2*i]=0;
		tmax[2*i+1]=0;
	}
	
	if (eventNum%5000==0) std::cout<<"\r"<<eventNum<<flush;
   	
	dynodeTraceh.clear();
	dynodeTracel.clear();
	hagridTrace.clear();
  xaTracel.clear();
  xbTracel.clear();
  yaTracel.clear();
  ybTracel.clear();
  xaTraceh.clear();
  xbTraceh.clear();
  yaTraceh.clear();
  ybTraceh.clear();
  anodeSum.clear();
  xaqdcSumsVec.clear();
  xbqdcSumsVec.clear();
  yaqdcSumsVec.clear();
  ybqdcSumsVec.clear();
  derxa.clear();
  derxb.clear();
  derya.clear();
  deryb.clear();
  derdh.clear();
  derdl.clear();
  der2dl.clear();
  der2dh.clear();
  hagTimean=0.0;
  hagTimed=0.0;
  dynTime=0.0;
  hrTimeha = 0.0; hrTimedh = 0.0;hrTimehd = 0.0;

  for(auto itC = rd.begin(); itC!=rd.end();itC++){
	  trace = &(itC->trace);
    vector<unsigned int> onboardQDCs = itC->qdcSums;
	  double energy = itC->energy;
	  double time = itC->time;
    double highResTime = itC->highResTime;
	  std::string dtype = itC->subtype.Data();
	  std::string dgroup = itC->group.Data();
	  int channel = itC->chanNum;  //used to try and gate on channel number rather than group/subtype
	  int module = itC->modNum;  //used to try and gate on channel number rather than group/subtype


    if (channel == 1 && module==0){
        iDyn++;
	      dynodeTracel = *trace;
	      dynTime = time;
        dynodeElow = energy;
        if(dynodeTracel.size()==0) continue;
        maxLoc = std::distance(dynodeTracel.begin(),std::max_element(dynodeTracel.begin(),dynodeTracel.end()));
	      minimumLoc = std::distance(dynodeTracel.begin(),std::min_element(dynodeTracel.begin(),dynodeTracel.end()));
        Bound = boundsCalc(dynodeTracel, maxLoc);
	      traceMax = dynodeTracel[maxLoc];
	      //if(traceMax>16380) maxLoc += 5;
	      traceMin = dynodeTracel[minimumLoc];
        stddev = baselineCalc(dynodeTracel).second;
        baseline = baselineCalc(dynodeTracel).first;
        unsigned int lowBoundShort = Bound.first;
        unsigned int highBoundShort = Bound.second;
        unsigned int lowBoundLong = 1 + highBoundShort;
        unsigned int highBoundLong = 300;
        qdcShortL = QDCcalculator(dynodeTracel, maxLoc-7, maxLoc+13);
	      qdcLongL = QDCcalculator(dynodeTracel, lowBoundLong, highBoundLong);
        for(int i=0;i<dynodeTracel.size()-1;i++){
          derdl.emplace_back(double(dynodeTracel.at(i+1))-double(dynodeTracel.at(i)));
        }
        for(int j=0;j<derdl.size()-1;j++){
          der2dl.emplace_back(double(derdl.at(j+1))-double(derdl.at(j)));
        }
    }
//    if (dtype == "hagrid_one"){
//        ihagd++;
//        hagTimed = time;
//        hagEdy = energy;
//        hrTimehd = highResTime;
//      }
//
//    else if(dtype == "hagrid_two"){
//	      ihaga++;
//	      hagridTrace = *trace;
//        hagTimean = time;
//        hagEan = energy;
//        hrTimeha = highResTime;
//      }
//
////    else if (dtype == "anode_low" && dgroup == "xa"){
//    else if(channel==4 && module==1){
//      ixa++;
//      xaTracel = *trace;
//      xal = energy;
//      maxLoc = std::distance(xaTracel.begin(),std::max_element(xaTracel.begin(),xaTracel.end()));
//      /*tmax[0] = *std::max_element(xaTracel.begin(),xaTracel.end());
//      Bound = boundsCalc(xaTracel, maxLoc);
//      if(maxLoc>450) continue;*/
//      if(xaTracel.size()==0) continue;
//      qdcSum[0] = QDCcalculator(xaTracel,maxLoc-5,maxLoc+35);
//      //bxa = baselineCalc(xaTracel).first;
//      }
////    else if (dtype == "anode_low" && dgroup == "xb"){
//    else if(channel==5 && module==1){
//      ixb++;
//      xbTracel = *trace;
//      xbl = energy;
//      maxLoc = std::distance(xbTracel.begin(),std::max_element(xbTracel.begin(),xbTracel.end()));
//      /*tmax[1] = *std::max_element(xbTracel.begin(),xbTracel.end());
//      Bound = boundsCalc(xbTracel, maxLoc);
//      if(maxLoc>450) continue;*/
//      if(xbTracel.size()==0) continue;
//      qdcSum[1] = QDCcalculator(xbTracel,maxLoc-5,maxLoc+35);
//      //bxb = baselineCalc(xbTracel).first;
//      }
////    else if (dtype == "anode_low" && dgroup == "ya"){
//    else if(channel==6 && module==1){
//	    iya++;
//	    yaTracel = *trace;
//	    yal = energy;
//	    maxLoc = std::distance(yaTracel.begin(),std::max_element(yaTracel.begin(),yaTracel.end()));
//      /*tmax[2] = *std::max_element(yaTracel.begin(),yaTracel.end());
//      Bound = boundsCalc(yaTracel, maxLoc);
//      if(maxLoc>450) continue;*/
//      if(yaTracel.size()==0) continue;
//      qdcSum[2] = QDCcalculator(yaTracel,maxLoc-5,maxLoc+35);
//      //bya = baselineCalc(yaTracel).first;
//      }
////    else if (dtype == "anode_low" && dgroup == "yb"){
//    else if(channel==7 && module==1){
//	    iyb++;
//	    ybTracel = *trace;
//	    ybl = energy;
//      maxLoc = std::distance(ybTracel.begin(),std::max_element(ybTracel.begin(),ybTracel.end()));
//      /*tmax[3] = *std::max_element(ybTracel.begin(),ybTracel.end());
//      Bound = boundsCalc(ybTracel, maxLoc);
//      if(maxLoc>450) continue;*/
//      if(ybTracel.size()==0) continue;
//      qdcSum[3] = QDCcalculator(ybTracel,maxLoc-5,maxLoc+35);
//      //byb = baselineCalc(ybTracel).first;
//    }
//    else if(dtype == "dynode_high" || (module==0 && channel==0)){
//      //iDyn++;
//      dynodeEhigh = energy;
//      hrTimedh = time;
//      dynodeTraceh = *trace;
//      maxLoc = std::distance(dynodeTraceh.begin(),std::max_element(dynodeTraceh.begin(),dynodeTraceh.end()));
//      minimumLoc = std::distance(dynodeTraceh.begin(),std::min_element(dynodeTraceh.begin(),dynodeTraceh.end()));
//      Bound = boundsCalc(dynodeTraceh, maxLoc);
//      traceMax = dynodeTraceh[maxLoc];
//      //if(traceMax>16380) maxLoc += 5;
//      traceMin = dynodeTraceh[minimumLoc];
//      //stddev = baselineCalc(dynodeTraceh).second;
//      //baseline = baselineCalc(dynodeTraceh).first;
//      if(traceMin<baselineCalc(dynodeTraceh).first-3*stddev){undershoot++;}
//      unsigned int lowBoundShort = Bound.first;
//      unsigned int highBoundShort = Bound.second;
//      unsigned int lowBoundLong = highBoundShort;
//      unsigned int highBoundLong = 15+highBoundShort;
//      if(maxLoc>250){/*cout<<"Broken Bounds"<<endl;*/ continue;}
//      //qdcShortH = QDCcalculator(dynodeTraceh, maxLoc-7, maxLoc+13);
//      if(dynodeTraceh.size()==0) continue;
//      qdcShortH = QDCcalculator(dynodeTraceh, maxLoc-15, maxLoc+5);
//      qdcLongH = QDCcalculator(dynodeTraceh, lowBoundLong, highBoundLong);
//      bD = baselineCalc(dynodeTraceh).first;
//      //dphase = 2*pcfdAnalyzer(dynodeTraceh, bD).first;
//      //hrTimedh += 2*dphase;
//      for(int i=0;i<dynodeTraceh.size()-1;i++){
//        derdh.emplace_back(double(dynodeTraceh.at(i+1))-double(dynodeTraceh.at(i)));
//      }
//      for(int j=0;j<derdh.size()-1;j++){
//        der2dh.emplace_back(double(derdh.at(j+1))-double(derdh.at(j)));
//      }
//      //numPile = pileupAnalysis(dynodeTraceh);
//    }
//    else if (dtype == "anode_high" && dgroup == "xa"){
    else if(channel==0 && module==1){
	    //ixa++;
	    xaTraceh = *trace;      
	    xah = energy;
      //xatime = time;
      if(xaTraceh.size()>0){
        maxLoc = std::distance(xaTraceh.begin(),std::max_element(xaTraceh.begin(),xaTraceh.end()));
        //tmax[4] = *std::max_element(xaTraceh.begin(),xaTraceh.end());
        Bound = boundsCalc(xaTraceh, maxLoc);
        if(maxLoc>450){qdcSum[0]=-2000.; continue;}
        //qdcSum[4] = QDCcalculator(xaTraceh, maxLoc-15,maxLoc+40);
        qdcSum[4] = QDCcalculator(xaTraceh, maxLoc-15,maxLoc+100);
        //qdcS[4] = QDCcalculator(xaTraceh, Bound.first, Bound.second);
        //qdcL[4] = QDCcalculator(xaTraceh, Bound.second+1, 500);
        bxa = baselineCalc(xaTraceh).first;
        //xaphase = 2*pcfdAnalyzer(xaTraceh, bxa).first;
        xatime = xaphase;
        //tmax[4] = pcfdAnalyzer(xaTraceh, bxa).second;
      }
      else if(onboardQDCs.size()>3){
        xaqdcSumsVec = onboardQDCs;
        qdcSum[4] = onboardQDCs.at(0)-onboardQDCs.at(2);
      }
      //for(int i=0;i<xaTraceh.size()-1;i++){
      //  derxa.emplace_back(double(xaTraceh.at(i+1))-double(xaTraceh.at(i)));
      //}
      //tmax[4] = pcfdAnalyzer(derxa, bxa).second;
    }
//    else if (dtype == "anode_high" && dgroup == "xb"){
    else if(channel==1 && module==1){
	    //ixb++;
      xbTraceh = *trace;
	    xbh = energy;
      //xbtime = time;
      if(xbTraceh.size()>0){
        maxLoc = std::distance(xbTraceh.begin(),std::max_element(xbTraceh.begin(),xbTraceh.end()));
        //tmax[5] = *std::max_element(xbTraceh.begin(),xbTraceh.end());
        Bound = boundsCalc(xbTraceh, maxLoc);
        if(maxLoc>450){qdcSum[1]=-2000.; continue;}
        //qdcSum[5] = QDCcalculator(xbTraceh, maxLoc-15,maxLoc+40);
        qdcSum[5] = QDCcalculator(xbTraceh, maxLoc-15,maxLoc+100);
        //qdcS[5] = QDCcalculator(xbTraceh, Bound.first, Bound.second);
        //qdcL[5] = QDCcalculator(xbTraceh, Bound.second+1, 500);
        bxb = baselineCalc(xbTraceh).first;
        //xbphase = 2*pcfdAnalyzer(xbTraceh, bxb).first;
        xbtime = xbphase;
        //tmax[5] = pcfdAnalyzer(xbTraceh, bxb).second;
      }
      else if(onboardQDCs.size()>3){
        xbqdcSumsVec = onboardQDCs;
        qdcSum[5] = onboardQDCs.at(0)-onboardQDCs.at(2);
      }
      //for(int i=0;i<xbTraceh.size()-1;i++){
      //  derxb.emplace_back(double(xbTraceh.at(i+1))-double(xbTraceh.at(i)));
      //}
      //tmax[5] = pcfdAnalyzer(derxb, bxb).second;
    }
//    else if (dtype == "anode_high" && dgroup == "ya"){
    else if(channel==2 && module==1){
	    //iya++;
	    yaTraceh = *trace;      
	    yah = energy;
      //yatime = time;
      if(yaTraceh.size()>0){
        maxLoc = std::distance(yaTraceh.begin(),std::max_element(yaTraceh.begin(),yaTraceh.end()));
        //tmax[6] = *std::max_element(yaTraceh.begin(),yaTraceh.end());
        Bound = boundsCalc(yaTraceh, maxLoc);
        if(maxLoc>450){qdcSum[2]=-2000.; continue;}
        //qdcSum[6] = QDCcalculator(yaTraceh, maxLoc-5,maxLoc+35);
        qdcSum[6] = QDCcalculator(yaTraceh,maxLoc-15,maxLoc+100);
        //qdcS[6] = QDCcalculator(yaTraceh, Bound.first, Bound.second);
        //qdcL[6] = QDCcalculator(yaTraceh, Bound.second+1, 500);

        bya = baselineCalc(yaTraceh).first;
        //yaphase = 2*pcfdAnalyzer(yaTraceh, bya).first;
        yatime = yaphase;
        //tmax[6] = pcfdAnalyzer(yaTraceh, bya).second;
      }
      else if(onboardQDCs.size()>3){
        yaqdcSumsVec = onboardQDCs;
        qdcSum[6] = onboardQDCs.at(0)-onboardQDCs.at(2);
      }
      //for(int i=0;i<yaTraceh.size()-1;i++){
      //  derya.emplace_back(double(yaTraceh.at(i+1))-double(yaTraceh.at(i)));
      //}
      //tmax[6] = pcfdAnalyzer(derya, bya).second;
    }
//    else if (dtype == "anode_high" && dgroup == "yb"){
    else if(channel==3 && module==1){
	    //iyb++;
     	ybTraceh = *trace;
	    ybh = energy;
      //ybtime = time;
      if(ybTraceh.size()>0){
        maxLoc = std::distance(ybTraceh.begin(),std::max_element(ybTraceh.begin(),ybTraceh.end()));
        tmax[7] = *std::max_element(ybTraceh.begin(),ybTraceh.end());
        Bound = boundsCalc(ybTraceh, maxLoc);
        if(maxLoc>450){qdcSum[3]=-2000.; continue;}
        //qdcSum[7] = QDCcalculator(ybTraceh, maxLoc-5,maxLoc+35);
        qdcSum[7] = QDCcalculator(ybTraceh, maxLoc-15,maxLoc+100);
        //qdcS[7] = QDCcalculator(ybTraceh, Bound.first, Bound.second);
        //qdcL[7] = QDCcalculator(ybTraceh, Bound.second+1, 500);
        byb = baselineCalc(ybTraceh).first;
        //ybphase = 2*pcfdAnalyzer(ybTraceh, byb).first;
        //ybtime = ybphase;
        //tmax[7] = pcfdAnalyzer(ybTraceh, byb).second;
      }
      else if(onboardQDCs.size()>3){
        ybqdcSumsVec = onboardQDCs;
        qdcSum[7] = onboardQDCs.at(0)-onboardQDCs.at(2);
      }
      //for(int i=0;i<ybTraceh.size()-1;i++){
      //  deryb.emplace_back(double(ybTraceh.at(i+1))-double(ybTraceh.at(i)));
      //}
      //ybphase = 2*pcfdAnalyzer(deryb, byb).first;
      //tmax[7] = pcfdAnalyzer(deryb, byb).second;
    }

    }

      //Section for filling plots
    if(ixa>0 && ixb>0 && iya>0 && iyb>0){
      qdcShortAl = 0.0;
	    qdcLongAl = 0.0;
      qdcShortAh = 0.0;
	    qdcLongAh = 0.0;
	    /*for(int i=0;i<4;i++){
	      scale = qdcS[0]/qdcS[i];
	      qdcSum[i] = qdcS[i] + qdcL[i];
	      qdcShortAl += scale*qdcS[i];
	      scale = qdcL[0]/qdcL[i];
	      qdcLongAl += scale*qdcL[i];
	      anodeRatio[i] = qdcS[i]/qdcSum[i];
	      scale = qdcS[4]/qdcS[i+4];
	      qdcShortAh += scale*qdcS[i+4];
	      scale = qdcL[4]/qdcL[i+4];
	      qdcLongAh += scale*qdcL[i+4];
	    }*/
      int maxQ=0;
      //for(int i=0;i<4;i++){
      //  //qdcSum[i] = qdcS[i] + qdcL[i];
      //  if(qdcSum[i]>qdcSum[maxQ]) maxQ=i;
      //}
      //for(int i=0;i<4;i++) qdcSum[i] = qdcSum[i]/qdcSum[maxQ];

    }

  /*  if(qdcSum[4]>0 && qdcSum[5]>0 && qdcSum[6]>0 && qdcSum[7]>0){
      xposh = position(qdcSum[4],qdcSum[5],qdcSum[6],qdcSum[7]).first;
      yposh = position(qdcSum[4],qdcSum[5],qdcSum[6],qdcSum[7]).second;
      for(int i=0;i<300;i++){
        anodeSum.push_back(xaTraceh[i]+xbTraceh[i]+yaTraceh[i]+ybTraceh[i]);
      }
    }*/
    if(qdcSum[4]>0 && qdcSum[5]>0 && qdcSum[6]>0 && qdcSum[7]>0){
      xposh = position(qdcSum[4],qdcSum[5],qdcSum[6],qdcSum[7]).first;
      yposh = position(qdcSum[4],qdcSum[5],qdcSum[6],qdcSum[7]).second;
    }
    if(xah>0 && xbh>0 && yah>0 && ybh>0){
      xposl = position(xah, xbh, yah, ybh).first;
      yposl = position(xah, xbh, yah, ybh).second;
    }

   /*if(tmax[4]>0 && tmax[5]>0 && tmax[6]>0 && tmax[7]>0){
      xposl = position(tmax[4],tmax[5],tmax[6],tmax[7]).first;
      yposl = position(tmax[4],tmax[5],tmax[6],tmax[7]).second;
    }*/
/*    if(tmax[4]>0 && tmax[5]>0 && tmax[6]>0 && tmax[7]>0){
      xposl = position(tmax[4], tmax[5], tmax[6], tmax[7]).first;
      yposl = position(tmax[4], tmax[5], tmax[6], tmax[7]).second;
    }*/
    /*if(tmax[0]>0 && tmax[1]>0 && tmax[2]>0 && tmax[3]>0){
      xposl = position(tmax[0], tmax[1], tmax[2], tmax[3]).first;
      yposl = position(tmax[0], tmax[1], tmax[2], tmax[3]).second;
    }*/

  if((qdcSum[4]>0 && qdcSum[5]>0 && qdcSum[6]>0 && qdcSum[7]>0) ||(tmax[0]>0 && tmax[1]>0 && tmax[2]>0 && tmax[3]>0) || dynodeElow>0 || dynodeEhigh>0 || 1){ 
    newTree->Fill();
  }
  else
    noCorr++;

	eventNum++;

  }
  
  newTree->Write();
  _file0->TFile::Close();

  cout<<endl<<noCorr<<" Events without anode correlation"<<endl;
  }

//Function for calculating the QDC

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

//Function for calculating upper and lower bounds for short due to pulse height
pair <int, int> boundsCalc(vector<unsigned int> trace, int maxPos){
    double maxVal = trace[maxPos];
    return std::make_pair (maxPos-10,maxPos+5);
}

pair <double, double> baselineCalc(vector<unsigned int> trace){
    //calculating the baseline
    int tSize = trace.size();
    if(tSize==0) return std::make_pair(-999.,-999.);
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

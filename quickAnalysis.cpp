#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <PaassRootStruct.hpp>
#include <cmath>
#include <string>
#include "position.h"
#include "pileupAnalysis.h"
#include "functions.hpp"

void quickAnalysis(){
  char* fileName = "crate2highgainRDrun301_DD.root";
  TFile *_file0 = TFile::Open(fileName);
  TTree *GS = (TTree*)_file0->Get("PixTree");
  TTreeReader singe;
  cout<<"readermade"<<endl;
  singe.SetTree( GS );
  TTreeReaderArray<processor_struct::ROOTDEV> rd = {singe, "rootdev_vec_"};  //gives vector of stucture 

    TH1D *plot1=new TH1D("plot1","plot1",10000.,-2000.,2000.);
    TH2D *energy=new TH2D("energy","energy",5000.,0.,10000.,5000.,0.,10000.);
    TH2D *phase=new TH2D("phase","phase",1000.,80.,105.,1000.,80.,105.);

    //variables
    int eventNum = 0;
    double energyone=0.0, energytwo=0.0, energythree=0.0;
    vector<double> traceone, tracetwo, xatrace, xbtrace, yatrace, ybtrace;
    vector<unsigned int> tone, ttwo, xat, xbt, yat, ybt;
    int ione=0, itwo=0, ithree=0, twocount=0;
    double cfdTimeone=0.0, cfdTimetwo=0.0, cfdTimethree=0.0;
    double phaseone=0.0, phasetwo=0., tmaxone=0., tmaxtwo=0.;
    double timeone=0.0, timetwo=0.0, timethree=0.0;
    double qdcone=0.0, qdctwo=0.0;
    double tdiff=0.0, baselineone=0.0, baselinetwo=0.0;
    double ecaltwo=0.0;
    int maxvalone=0, maxvaltwo=0, pileNum=0;
    double xa=0., xb=0., ya=0., yb=0.;
    double xaqdc=0., xbqdc=0., yaqdc=0., ybqdc=0.;
    double xatm=0., xbtm=0., yatm=0., ybtm=0.;
    double xposq=0., yposq=0.;
    double xpost=0., ypost=0.;
    bool valid = false;

  TFile *newFile = new TFile("crate2highgainRDrun301.root","RECREATE");
  TTree *timing = new TTree("timing","tree filled with traces and energies etc.");

  timing->Branch("eventNum", &eventNum);
  timing->Branch("energyone", &energyone);
  timing->Branch("energytwo", &energytwo);
  timing->Branch("energythree", &energythree);
  timing->Branch("energycal", &ecaltwo);
  timing->Branch("traceone", &traceone);
  timing->Branch("tracetwo", &tracetwo);
  timing->Branch("xaTrace", &xatrace);
  timing->Branch("xbTrace", &xbtrace);
  timing->Branch("yaTrace", &yatrace);
  timing->Branch("ybTrace", &ybtrace);
  timing->Branch("ione", &ione);
  timing->Branch("itwo", &itwo);
  timing->Branch("cfdTimeone", &cfdTimeone);
  timing->Branch("cfdTimetwo", &cfdTimetwo);
  timing->Branch("phaseone", &phaseone);
  timing->Branch("phasetwo", &phasetwo);
  timing->Branch("timeone", &timeone);
  timing->Branch("timetwo", &timetwo);
  timing->Branch("tmaxone", &tmaxone);
  timing->Branch("tmaxtwo", &tmaxtwo);
  timing->Branch("qdcone", &qdcone);
  timing->Branch("qdctwo", &qdctwo);
  timing->Branch("tdiff", &tdiff);
  timing->Branch("xposq", &xposq);
  timing->Branch("yposq", &yposq);
  timing->Branch("xpost", &xpost);
  timing->Branch("ypost", &ypost);
  timing->Branch("bone", &baselineone);
  timing->Branch("btwo", &baselinetwo);
  timing->Branch("pile", &pileNum);
  timing->Branch("valid", &valid);

  std::vector<unsigned int> *trace;
  while(singe.Next()){

    if(eventNum%10000==0) cout<<"\r" << eventNum<< flush;
    tone.clear(); ttwo.clear();
    traceone.clear();
    tracetwo.clear();
    xat.clear(); xbt.clear(); yat.clear(); ybt.clear();
    xatrace.clear(); xbtrace.clear();
    yatrace.clear(); ybtrace.clear();
    ione=0; itwo=0; ithree=0;
    energyone= -1000.0; energytwo= -1000.0; energythree= -1000.0;
    cfdTimeone=0.0; cfdTimetwo=0.0; cfdTimethree=0.;
    phaseone=0.0; phasetwo=0.0;
    timeone=0.0; timetwo=0.0, timethree=0.0;
    tmaxone=0.0; tmaxtwo=0.0; tdiff=0.0;
    ecaltwo=0.0; baselineone=0.0; baselinetwo=0.0;
    qdcone=0.0;qdctwo=0.0;
    maxvalone=0;maxvaltwo=0;pileNum=0;
    xa=0; xb=0.; ya=0.; yb=0.;
    xposq=0.; yposq=0.;
    xpost=0.; ypost=0.;
    xaqdc=0.; xbqdc=0.; yaqdc=0.; ybqdc=0.;
    xatm=0.; xbtm=0.; yatm=0.; ybtm=0.;
    valid = false;

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

        if (module==0 && channel==0){ 
            ione++;
            tone = *trace;
            double hold =0.0;
            int maxLoc = 0;
            for(int i=0;i<tone.size();i++){
              hold = (double) tone[i];
              traceone.push_back(hold);
              if(traceone[maxLoc]<=traceone[i]) maxLoc=i;
            }
            //for(int i=0;i<traceone.size()-2;i++)traceone[i] = (traceone[i]+traceone[i+1])/2;
            energyone = energy;
            cfdTimeone = time;
            tmaxone = traceone.at(maxLoc);
            timeone = time;
            maxvalone = maxValue;
            //phaseone = cdfAnalyzer(traceone);
            //cfdTimeone += 2*phaseone;
            qdcone = QDCcalculator(traceone);
            pileNum = pileupAnalysis(traceone);
        }

        else if(channel==5 && false){ //external trig
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
            cfdTimetwo = time;
            timetwo = time;
            maxvaltwo=maxValue;
            //phasetwo = cdfAnalyzer(tracetwo);
            //cfdTimetwo += 2*phasetwo;
            qdctwo = QDCcalculator(tracetwo);
            baselinetwo = baselineCalc(tracetwo).first;
        }
        else if(channel==6 && false){ //hagrid
            ithree++;
            energythree = energy;
            //cfdTimethree = time;
        }

        else if(channel==4){ 
          xat = *trace;
          double hold =0;
          for(int i=0;i<xat.size();i++){
              hold = (double) xat[i];
              if(hold>xatm) xatm = hold;
              xatrace.push_back(hold);
            }
          xa = energy;
          xaqdc = QDCcalculator(xatrace);
        }
        else if(channel==5){ 
          xbt = *trace;
          double hold =0;
          for(int i=0;i<xbt.size();i++){
              hold = (double) xbt[i];
              if(hold>xbtm) xbtm = hold;
              xbtrace.push_back(hold);
            }
          xb = energy;
          xbqdc = QDCcalculator(xbtrace);
        }
        else if(channel==6){ 
          yat = *trace;
          double hold =0;
          for(int i=0;i<yat.size();i++){
              hold = (double) yat[i];
              if(hold>yatm) yatm = hold;
              yatrace.push_back(hold);
            }
          ya = energy;
          yaqdc = QDCcalculator(yatrace);
        }
        else if(channel==7){ 
          ybt = *trace;
          double hold =0;
          for(int i=0;i<ybt.size();i++){
              hold = (double) ybt[i];
              if(hold>ybtm) ybtm = hold;
              ybtrace.push_back(hold);
            }
          yb = energy;
          ybqdc = QDCcalculator(ybtrace);
        }
    }
    if(xaqdc>0 && xbqdc>0 && yaqdc>0 && ybqdc>0){
      xposq = position(xaqdc,xbqdc,yaqdc,ybqdc).first;
      yposq = position(xaqdc,xbqdc,yaqdc,ybqdc).second;
      xpost = position(xatm,xbtm,yatm,ybtm).first;
      ypost = position(xatm,xbtm,yatm,ybtm).second;
      valid = true;
    }
    else if(xa>0 && xb>0 && ya>0 && yb>0){
      xpost = position(xa,xb,ya,yb).first;
      ypost = position(xa,xb,ya,yb).second;
      valid = true;
    }
    timing->Fill();
    
    eventNum++;
  }
  timing->Write();
  //plot1->Write();
  //phase->Write();
  //energy->Write();
  cout<<endl;
}

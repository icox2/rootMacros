#include "TSystem.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TChain.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <ProcessorRootStruc.hpp>

void xpeakFitting(){
   cout<<"in loop"<<endl;
    int numPeak=0;
    std::cout<<"Enter number of peaks to fit (max 5)."<<std::endl;
    std::cin>>numPeak;
    if(numPeak>5) numPeak=5;
    std::vector<double> lb (numPeak,0);
    std::vector<double> ub (numPeak,0);
    int lowbound=100000, highbound=0;
    for(int i=0;i<numPeak;i++){
        std::cout<<"Enter lower bound for peak "<<i<<std::endl;
        std::cin>>lb[i];
        std::cout<<"Enter upper bound for peak "<<i<<std::endl;
        std::cin>>ub[i];
        if(lb[i]>=ub[i]){
            cout<<"Bad Bounds.";
            return;
        }
        if(lowbound>lb[i]) lowbound=lb[i];
        if(highbound<ub[i]) highbound=ub[i];
    }

    TF1 *g1 = new TF1("m1","gaus",lb[0],ub[0]);
    TF1 *g2 = new TF1("m2","gaus",lb[1],ub[1]);
    TF1 *g3 = new TF1("m3","gaus",lb[2],ub[2]);
    TF1 *g4 = new TF1("m4","gaus",lb[3],ub[3]);
    TF1 *g5 = new TF1("m5","gaus",lb[4],ub[4]);
    TF1 *total = new TF1("total","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)",lowbound,highbound);

    //std::cout<<"Enter histogram name"<<std::endl;
    //std::cin>>histName;

    if (numPeak>0) e1->Fit(g1,"R");
    if (numPeak>1) e1->Fit(g2,"R+");
    if (numPeak>2) e1->Fit(g3,"R+");
    if (numPeak>3) e1->Fit(g4,"R+");
    if (numPeak>4) e1->Fit(g5,"R+");
    
    double par[15] = {0};
    g1->GetParameters(&par[0]);
    g2->GetParameters(&par[3]);
    g3->GetParameters(&par[6]);
    g4->GetParameters(&par[9]);
    g5->GetParameters(&par[12]);

    bool good=0;
    std::cout<<"Total plot? (1-yes, 0-no)"<<std::endl;
    std::cin>>good;

    if(good){
        total->SetParameters(par);
        total->SetLineColor(30);
        e1->Fit("total","R+","Q");
    }
}
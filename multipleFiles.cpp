//taking multiple files and creating histograms to plot in the same TGraph

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

void multipleFiles(){

    TFile *_file0 = TFile::Open("hagrid2_wpspmt_241am_2_DD.root");
    TTree *GS = (TTree*)_file0->Get("PixTree");
    TTreeReader singe;
    cout<<"readermade0"<<endl;
    singe.SetTree( GS );
    TTreeReaderArray<processor_struct::ROOTDEV> rd = {singe, "root_dev_vec_"};  //gives vector of stucture 

    TFile *_file1 = TFile::Open("hagrid2_wpspmt_241am_DD.root");
    TTree *GS1 = (TTree*)_file1->Get("PixTree");
    TTreeReader singe1;
    cout<<"readermade1"<<endl;
    singe1.SetTree( GS1 );
    TTreeReaderArray<processor_struct::ROOTDEV> rd1 = {singe1, "root_dev_vec_"};  //gives vector of stucture 

    TH1D *plot1=new TH1D("plot1","plot1",5000.,0.,5000.);
    TH1D *plot2=new TH1D("plot2","plot2",5000.,0.,5000.);

    while(singe.Next()){
        for(auto itC = rd.begin(); itC!=rd.end();itC++){
	        double energy = itC->energy;
	        std::string dtype = itC->subtype.Data();
	        std::string dgroup = itC->group.Data();
    	    int channel = itC->chanNum;  //used to try and gate on channel number rather than group/subtype

            if(dtype=="hagrid_an") plot1->Fill(energy);
        }
    }
    while(singe1.Next()){
        for(auto itC = rd1.begin(); itC!=rd1.end();itC++){
	        double energy = itC->energy;
	        std::string dtype = itC->subtype.Data();
	        std::string dgroup = itC->group.Data();
    	    int channel = itC->chanNum;  //used to try and gate on channel number rather than group/subtype

            if(dtype=="hagrid_an") plot2->Fill(energy);
        }
    }
    plot1->Draw();
    plot2->SetLineColor(30);
    plot2->Draw("same");
}
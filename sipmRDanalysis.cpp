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

pair<double, double> posCalcRecon(vector<double> anodeEnergy);
pair<int, int> posCalcRaw(int pos);
pair<int, int> multCounter(vector<double> anodeEnergy);

void sipmRDanalysis(){
  TFile *_file0 = TFile::Open("run172RD_DD.root");
  TTree *GS = (TTree*)_file0->Get("PixTree");
  TTreeReader singe;
  cout<<"readermade"<<endl;
  singe.SetTree( GS );
  TTreeReaderArray<processor_struct::ROOTDEV> rd = {singe, "rootdev_vec_"};  //gives vector of stucture 

  vector<double> anodeEnergyHG(64,0);
  double xReconPosHG = 0., yReconPosHG = 0.;
  int xRawPosHG=0, yRawPosHG=0;
  int eventNum=0, multiplicityHG=0;
  int xMultHG=0, yMultHG=0;
	double energySumHG=0.;
  vector<double> anodeEnergyLG(64,0);
  double xReconPosLG = 0., yReconPosLG = 0.;
  int xRawPosLG=0, yRawPosLG=0;
  int multiplicityLG=0;
  int xMultLG=0, yMultLG=0;
	double energySumLG=0.;
	
  TFile *newFile = new TFile("run172RD.root","RECREATE");
  TTree *newTree = new TTree("newTree","tree filled with traces and energies etc.");
  
  newTree->Branch("eventNum", &eventNum);
  newTree->Branch("high_gain_xReconPos", &xReconPosHG);
  newTree->Branch("high_gain_yReconPos", &yReconPosHG);
  newTree->Branch("high_gain_xRawPos", &xRawPosHG);
  newTree->Branch("high_gain_yRawPos", &yRawPosHG);
  newTree->Branch("high_gain_multiplicity", &multiplicityHG);
  newTree->Branch("high_gain_xMult", &xMultHG);
  newTree->Branch("high_gain_yMult", &yMultHG);
	newTree->Branch("high_gain_anodeE",&anodeEnergyHG);
	newTree->Branch("high_gain_energySum", &energySumHG);
  newTree->Branch("low_gain_xReconPos", &xReconPosLG);
  newTree->Branch("low_gain_yReconPos", &yReconPosLG);
  newTree->Branch("low_gain_xRawPos", &xRawPosLG);
  newTree->Branch("low_gain_yRawPos", &yRawPosLG);
  newTree->Branch("low_gain_multiplicity", &multiplicityLG);
  newTree->Branch("low_gain_xMult", &xMultLG);
  newTree->Branch("low_gain_yMult", &yMultLG);
	newTree->Branch("low_gain_anodeE",&anodeEnergyLG);
	newTree->Branch("low_gain_energySum", &energySumLG);

  std::vector<unsigned> *trace;
  while(singe.Next()){	
		if(eventNum%10000==0) cout<<"\r" << eventNum<< flush;
		eventNum++;
		//clearing the previous values from last event
		xRawPosHG = 0; yRawPosHG=0;
		xReconPosHG = 0.; yReconPosHG=0.;
		for(int i=0;i<64;i++)
			anodeEnergyHG[i]=0.;
		int maxEnergyHG=0;
		multiplicityHG = 0; energySumHG=0.;
		xMultHG=0; yMultHG=0;
		xRawPosLG = 0; yRawPosLG=0;
		xReconPosLG = 0.; yReconPosLG=0.;
		for(int i=0;i<64;i++)
			anodeEnergyLG[i]=0.;
		int maxEnergyLG=0;
		xMultLG=0; yMultLG=0;
		multiplicityLG = 0;energySumLG=0.;

		for(auto itC = rd.begin(); itC!=rd.end();itC++){
			//trace = &(itC->trace);
			double energy = itC->energy;
			//double time = itC->time;
			//double highResTime = itC->highResTime;
			std::string dtype = itC->subtype.Data();
			int dgroup = itC->detNum;
      bool dsat = itC->saturation;
      // Saturated events are needed in order to properly reconstruct the position 
      //if(dsat) continue;
			// int groupNum = stoi(dgroup);
			int groupNum = dgroup;
			int channel = itC->chanNum;  //used to try and gate on channel number rather than group/subtype
			if(groupNum>=0 && groupNum<=64 && dtype=="anode_h"){
				multiplicityHG++;
				energySumHG+=energy;
				anodeEnergyHG[groupNum] = energy;
				if(energy>maxEnergyHG){
					maxEnergyHG = energy;
					auto tmpPair = posCalcRaw(groupNum);
					xRawPosHG = tmpPair.first;
					yRawPosHG = tmpPair.second;
				}
			}
			else if(groupNum>=0 && groupNum<=64 && dtype=="anode_l"){
				multiplicityLG++;
        //if(anodeEnergyLG[groupNum]>0) cout<<"!!groupNum already filled!!\n"<<groupNum<<endl;
				anodeEnergyLG[groupNum] = energy;
				energySumLG+=energy;
				if(energy>maxEnergyLG){
					maxEnergyLG = energy;
					auto tmpPair = posCalcRaw(groupNum);
					xRawPosLG = tmpPair.first;
					yRawPosLG = tmpPair.second;
				}
			}
			/*if (dtype == "hagrid_one"){
					ihagd++;
					hagTimed = time;
					hagEdy = energy;
					hrTimehd = highResTime;
				}*/
		}
		pair<double, double> reconPos = posCalcRecon(anodeEnergyHG);
		xReconPosHG = reconPos.first;
		yReconPosHG = reconPos.second;
		pair<int, int> mults = multCounter(anodeEnergyHG);
		xMultHG = mults.first;
		yMultHG = mults.second;
		reconPos = posCalcRecon(anodeEnergyLG);
		xReconPosLG = reconPos.first;
		yReconPosLG = reconPos.second;
		mults = multCounter(anodeEnergyLG);
		xMultLG = mults.first;
		yMultLG = mults.second;
		if((xMultLG>1 && yMultLG>1) || (xMultHG>1 && yMultHG>1))
			newTree->Fill();
		eventNum++;
  } 
  newTree->Write();
  _file0->TFile::Close();
}

pair<double, double> posCalcRecon(vector<double> anodeEnergy){
	double xpos=0, ypos=0;
	double total=0;
	for(int i=0;i<64;i++){
		if(anodeEnergy[i]>100000 || anodeEnergy[i]<50){
        //if(anodeEnergy[i]!=0)cout<<"anodeE: "<<anodeEnergy[i]<<endl;
			continue;}
		double xW = i%8;
		double yW = 7-i/8;
		total+= 7*anodeEnergy[i];
		xpos+= xW*anodeEnergy[i];
		ypos+= yW*anodeEnergy[i];
	}
	return make_pair(xpos/total,ypos/total);
}
pair<int, int> posCalcRaw(int pos){
	int ypos = 8-pos/8;
	int xpos = pos%8+1;	
	return make_pair(xpos,ypos);
}

pair<int, int> multCounter(vector<double> anodeEnergy){
	int xmult[8]={0}, ymult[8]={0};
	int xcnt=0, ycnt=0;
	for(int i=0;i<anodeEnergy.size();i++){
		if(anodeEnergy[i]>100000 || anodeEnergy[i]<1)
			continue;
		int ypos = 7-i/8;
		int xpos = i%8;	
		xmult[xpos]++;
		ymult[ypos]++;
	}
	for(int i=0;i<8;i++){
		if(xmult[i]>0) xcnt++;	
		if(ymult[i]>0) ycnt++;
	}
	return make_pair(xcnt,ycnt);
}

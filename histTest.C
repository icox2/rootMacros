#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

void histTest(){
	TFile *file = new TFile("crate1_run_0120-00_DD.root","read");

	TTree *tree = (TTree*)file->Get("PixTree");

	TH2F *tofHist = new TH2F("tofHist","tofHist",1000,-100,900,90,0,90);
	tree->Draw("vandle_vec_.barNum:vandle_vec_.tof>>tofHist","","colz");
	//ttof->Draw("colz");
}

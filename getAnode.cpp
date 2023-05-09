// getAnode.cpp
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <ProcessorRootStruc.hpp>
#include <cmath>
#include <string>
#include "position.h"
#include "pileupAnalysis.h"
#include "functions.hpp"

void getAnode(){
  char* fileName = "distortion_test_25g.root";
  TFile *_file0 = TFile::Open(fileName);
  TTree *tree = (TTree*)_file0->Get("data");

  double xrecon=0., yrecon=0., xphoton=0., yphoton=0.;
  int eventNum=0;

  tree->SetBranchAddress("reconComX",&xrecon);
  tree->SetBranchAddress("reconComY",&yrecon);
  tree->SetBranchAddress("photonComX",&xphoton);
  tree->SetBranchAddress("photonComY",&yphoton);

  for(int iEntry=0;tree->LoadTree(iEntry)>=0;++iEntry){
    //xrecon=0., yrecon=0., xphoton=0., yphoton=0.;
    tree->GetEntry(iEntry);
    std::cout<<xrecon<<std::endl;
    if(iEntry>5)
      break;
  }
}

#include "TGraph.h"
#include "TProfile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"

void walkCor(TH2F *his){
    TGraph *w = new TGraph();

    for(int i=1;i<1001;i++){
        int bm = his->ProjectionX("",i,i)->GetMaximumBin();
        double cent = bm *20./1000.-50;
        w->SetPoint(i-1,i*600000./1000.,cent);
    }
    TCanvas *c2 = new TCanvas;
    w->Draw();
}
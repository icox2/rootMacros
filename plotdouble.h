#include "TGraph.h"
#include "TProfile.h"
#include "doublepeak.h"
#include "TF1.h"
// double res1 = p[0]*exp(-1*(x[0]-p[1])/p[2])*(1-exp(-1*pow(x[0]-p[1],3)/p[3]));

// double res2 = p[5]*exp(-1*(x[0]-p[6])/p[7])*(1-exp(-1*pow(x[0]-p[6],3)/p[8]));

void plotdouble(TGraph *gn){
	gn->GetListOfFunctions()->Clear();
	gn->GetXaxis()->SetRangeUser(100,300);
	double plow[9] = {1000,165,0.1,1.0,-20,1000,210,0.1,1.0}; 
	double phi[9] = {20000,180,100,3000,20,20000,230,100,3000}; 
	TF1 *f1 = new TF1("f1",doublepeak,0,600,9);
	f1->SetNpx(1000);
	for(int i=0;i<9;i++) f1->SetParLimits(i,plow[i],phi[i]);
	//f1->SetParLimits(1,170,180);
	//f1->SetParLimits(4,500,1000);
	//f1->SetParLimits(6,175,200);
	gn->Fit(f1,"N","",150,250);
	//gn->Fit(f1,"N","",160,200);
	gn->GetXaxis()->SetRangeUser(150,350);
	gn->Draw();
	f1->Draw("same");
	return;
}

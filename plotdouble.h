#include "TGraph.h"
#include "TProfile.h"
#include "doublepeak.h"
#include "singlepeak.h"
#include "TF1.h"
#include <cmath>
#include <iostream>
#include <iomanip>
// double res1 = p[0]*exp(-1*(x[0]-p[1])/p[2])*(1-exp(-1*pow(x[0]-p[1],3)/p[3]));

// double res2 = p[5]*exp(-1*(x[0]-p[6])/p[7])*(1-exp(-1*pow(x[0]-p[6],3)/p[8]));

void plotdouble(TGraph *gn){
	gn->GetListOfFunctions()->Clear();
	gn->GetXaxis()->SetRangeUser(100,400);
	//double plow[12] = {0.1,165,0.1,1.0,-10,1.1,0.1,190,0.1,1.0,-10,1.1}; 
	//double phi[12] = {10,190,100,100,0,1.8,10,210,100,100,0,1.8}; 
	double plow[12] = {0.1,128,0.1,10,-0.2,1.01,0.1,129,0.1,1.0,-0.2,1.1}; 
	double phi[12] = {3,140,100,70,0,1.8,3,270,100,70,0,1.8}; 
	TF1 *f1 = new TF1("f1",doublepeak,0,600,12);
	//TF1 *f2 = new TF1("f2","(x>[1])*([0]*exp(-1*(x-[1])/[2])*(1-exp(-1*pow(x-[1],[5])/[3]))+[4])",0,600);
	//TF1 *f3 = new TF1("f3","(x>[1])*([0]*exp(-1*(x-[1])/[2])*(1-exp(-1*pow(x-[1],[5])/[3]))+[4])",0,600);
	TF1 *f2 = new TF1("f2",singlepeak,0,600,12);
	TF1 *f3 = new TF1("f3",singlepeak,0,600,12);
	f1->SetNpx(1000);
	for(int i=0;i<12;i++) f1->SetParLimits(i,plow[i],phi[i]);
	//f1->FixParameter(0,3.64);
	gn->Fit(f1,"N","",150,450);
	gn->GetXaxis()->SetRangeUser(150,350);
	gn->Draw();
	f1->Draw("same");
	for(int i=0;i<6;i++){
		f2->SetParameter(i,f1->GetParameter(i));
		f3->SetParameter(i,f1->GetParameter(i+6));
		cout<<i<<" "<<setw(10)<<f1->GetParameter(i)<<" "<<i+6<<" "<<f1->GetParameter(i+6)<<endl;
	}
	f2->SetLineColor(3);
	f3->SetLineColor(4);
	f2->Draw("same");
	f3->Draw("same");
	return;
}

#include "TGraph.h"
#include "TProfile.h"
#include "tailfit.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH1D.h"

void plottail(TH1F *his){
    //his->GetListOfFunctions()->Clear();
	his->GetXaxis()->SetRangeUser(100,300);
	double plow[4] = {10.,-55.,-10.,0.6};
	double phi[4] = {1000.,-30,-0.1,7.0};
	TF1 *f1 = new TF1("f1",tailfit,20,35,4);
	cout<<"in"<<endl;
	f1->SetNpx(1000);
	for(int i=0;i<5;i++) f1->SetParLimits(i,plow[i],phi[i]);
	//f1->SetParameter(3,0.49);
	his->Fit(f1,"N","",20,-35);
	his->GetXaxis()->SetRangeUser(20,35);
	his->Draw();
	f1->Draw("same");
	return;
}

#include "TGraph.h"
#include "TProfile.h"
#include "tailfit.h"
#include "TF1.h"
#include "TH1F.h"

void plottail(TH1F *his){
    //his->GetListOfFunctions()->Clear();
	his->GetXaxis()->SetRangeUser(100,300);
	double plow[5] = {100.,-39.,0.5,10.,0.2};
	double phi[5] = {500.,-35.,3.0,100.,1.0};
	TF1 *f1 = new TF1("f1",tailfit,-55,-25,5);
	cout<<"in"<<endl;
	f1->SetNpx(1000);
	//for(int i=0;i<5;i++) f1->SetParLimits(i,plow[i],phi[i]);
	f1->SetParameter(3,0.49);
	his->Fit(f1,"N","",-55,-25);
	his->GetXaxis()->SetRangeUser(-55,-25);
	his->Draw();
	f1->Draw("same");
	return;
}

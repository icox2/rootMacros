#include <cmath>
double doublepeak(double *x, double *p){

	double res1 = p[0]*exp(-1*(x[0]-p[1])/p[2])*(1-exp(-1*pow(x[0]-p[1],p[5])/p[3]));
	double res2 = p[6]*exp(-1*(x[0]-p[7])/p[8])*(1-exp(-1*pow(x[0]-p[7],p[11])/p[9]));

	if(x[0]>=p[1] && x[0]<p[7]) return res1 + p[4];
	else if(x[0]>=p[7]) return res1 + res2 + p[4] + p[10];
	else if(x[0]>=p[1] && x[0]<p[6] && res1<0) return 0;
	else return 0;
}

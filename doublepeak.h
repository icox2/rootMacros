#include <cmath>
double doublepeak(double *x, double *p){

	double res1 = p[0]*exp(-1*(x[0]-p[1])/p[2])*(1-exp(-1*pow(x[0]-p[1],4)/p[3]));
		
	double res2 = p[5]*exp(-1*(x[0]-p[6])/p[7])*(1-exp(-1*pow(x[0]-p[6],4)/p[8]));

	if(x[0]>=p[1] && x[0]<p[6]) return res1 + p[4];
	else if(x[0]>=p[6]) return res1 + res2 + p[4];
	else return p[4];
}

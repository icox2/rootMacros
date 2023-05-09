//singlepeak.h
#include <cmath>
double singlepeak(double *x, double *p){

	double res1 = p[0]*exp(-1*(x[0]-p[1])/p[2])*(1-exp(-1*pow(x[0]-p[1],p[5])/p[3]));

	if(x[0]>=p[1]) return res1 + p[4];
	else return 0;
}

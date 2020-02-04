#include <cmath>

double tailfit(double *x, double *p){
	
/*	double gauss = p[0]*exp(-0.5*(pow((x[0]-p[1])/p[2],2)));
	double expo = exp(p[3]+p[4]*x[0]);

	if(x[0]<-38.66) return expo;
	else return gauss
*/
	double comb = p[0]*exp(-0.5*(pow((x[0]-p[1])/p[2],2)))*exp(-p[3]*(p[1]-p[4]));

	return comb;
}

#include <cmath>
#include "TMath.h"
using namespace TMath;

double tailfit(double *x, double *p){
	
/*	double gauss = p[0]*exp(-0.5*(pow((x[0]-p[1])/p[2],2)));
	double expo = exp(p[3]+p[4]*x[0]);

	if(x[0]<-38.66) return expo;
	else return gauss
*/
	//double comb = p[0]*exp(-0.5*(pow((x[0]-p[1])/p[2],2)))*exp(-p[3]*(p[1]-p[4]));
	//h = p[0], amp of gaussian, mu = p[1], sigma = p[2], tau = p[3]/log(2), for half-life
	//for larger tau this function might need to be used, including the error function
	//double func = p[0]*p[2]/(p[3]/log(2))*sqrt(Pi()/2)*exp(0.5*pow(p[2]/(p[3]/log(2)),2)-(x[0]-p[1])/(p[3]/log(2)))*Erfc(1/sqrt(2)*(p[2]/(p[3]/log(2))-(x[0]-p[1])/p[2]));
	//p[0] is amplitude, lambda = p[2], mean = p[1], sigma = p[3]
	double func = [0]*exp([2]/2*(2*[1]+[2]*pow([3],2)-2*x))*erfc(([1]+[2]*pow([3],2)-x)/(sqrt(2)*[3]));
	return func;
}

#ifndef GFUNC
#define GFUNC

double gfunc(double *x, double *p){
  double sum = 0;
  int arrLen = sizeof(p)/sizeof(double);
  arrLen = 20;
  double fwhmRatio = 0.1;
  for(int i=0;i<arrLen/2;i++){
    sum += p[2*i]*exp(-0.5*pow((x[0]-p[2*i+1])/(fwhmRatio*p[2*i+1]/2.3548),2));
  }
  return sum;
//  double g1 = p[0]*exp(-0.5*pow((x[0]-p[1])/(0.1*p[1]/2.3548),2));
//  double g2 = p[2]*exp(-0.5*pow((x[0]-p[3])/(0.1*p[3]/2.3548),2));
//  double g3 = p[4]*exp(-0.5*pow((x[0]-p[5])/(0.1*p[5]/2.3548),2));
//  double g4 = p[6]*exp(-0.5*pow((x[0]-p[7])/(0.1*p[7]/2.3548),2));
//  double g5 = p[8]*exp(-0.5*pow((x[0]-p[9])/(0.1*p[9]/2.3548),2));
//  return g1+g2+g3+g4+g5;
}

#endif

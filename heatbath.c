#include <stdio.h>
#include <math.h>
#include "functions.h"
int main() {

const int N = 1000;
const double gamma = 1.2;
double omegaMin=pow(N,-0.7988), omegaMax=pow(N,1.0688);
double omega[N];
double oscMass;
double M;
double invM[N];
double p[N+1];
double q[N+1];
double k[N];


computeOmega(omega,omegaMin,omegaMax,N);
//printf("%f", omega[999]);
//computeMasses(invM,oscMass,omega,omegaMin,gamma);
printArray_(omega);







}

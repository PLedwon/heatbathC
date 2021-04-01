//#include <cstdio>
#include <iostream>
#include <cmath>
#include "functions.h"
int main() {

const int N = 10;
const double gamma = 1.2;
double oscMass = 1.0;
double M = 1.0;
double omegaMin=pow(N,-0.7988), omegaMax=pow(N,1.0688);
double omega[N];
double masses[N];
double k[N];
double invM[N+1];
double q0[N+1];
double p0[N+1];
double p[N+1];
double q[N+1];

// setting the bathparameters
computeOmega(omega,omegaMin,omegaMax,N);
computeMasses(masses,oscMass,omega,omegaMin,gamma,N);
computeSpringConstants(k, masses, omega, N);
invertMasses(invM, M, masses, N);

//set the initial conditions for the distinguished particle
q0[0]=1;
//p0[0]=0;

//printf("%f", H(q,p,k,invM,N));
//printArray_(q0,N+1);







}

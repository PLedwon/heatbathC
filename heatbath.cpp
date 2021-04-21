//#include <cstdio>
#include <iostream>
#include <cmath>
#include "functions.h"
#include <random>
int main() {

const int N = 1000;
const double gamma = 1.2;
const double beta = 1.0;
const double tspan[2] = {0, pow(10,1)};
double oscMass = 1.0;
double M = 1.0;
double omegaMin=pow(N,-0.7988), omegaMax=omegaMin*pow(N,1.0688);
double omega[N];
double masses[N];
double k[N+1];
double invM[N+1];
double q0[N+1];
double p0[N+1];
double q1[N+1];
double p1[N+1];

// setting the bathparameters
computeOmega(omega,omegaMin,omegaMax,N);
computeMasses(masses,oscMass,omega,omegaMin,gamma,N);
computeSpringConstants(k, masses, omega, N);
invertMasses(invM, M, masses, N);
generateInitialConditions(q0,p0,M,masses,k,beta,N);


// solve EOM and save Q at times tsave

//save Q to file

printArray_(masses,N);









}

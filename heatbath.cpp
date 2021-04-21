//#include <cstdio>
#include <iostream>
#include <cmath>
#include "functions.h"
int main() {

const int N = 100; //# of harmonic oscillators in our heatbath
const double GAMMA = 1.2; // expected superdiffusion exponent
const double BETA = 1.0; //kB*T
const double TSPAN[2] = {0, pow(10,1)};
const double DT = pow(10,-2);
double oscMass = 1.0; //mass of heaviest bath oscillator
double M = 1.0; // mass of distinguished particle
double omegaMin=pow(N,-0.7988), omegaMax=omegaMin*pow(N,1.0688); //highest and lowest eigenfrequency of the bath
double omega[N];
double masses[N];
double k[N+1];
double invM[N+1];
double q0[N+1];
double p0[N+1];
double q1[N+1];
double p1[N+1];


// setting the bathparameters
setEigenfrequencies(omega,omegaMin,omegaMax,N);
computeMasses(masses,oscMass,omega,omegaMin,GAMMA,N);
computeSpringConstants(k, masses, omega, N);
invertMasses(invM, M, masses, N);

// set initial conditions and solve the eom
generateInitialConditions(q0,p0,M,masses,k,BETA,N);

solveEOM(q1,p1,q0,p0,invM,k,TSPAN,DT,N);


// solve EOM and save Q at times tsave

//save Q to file

//printArray_(masses,N);









}

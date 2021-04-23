//#include <cstdio>
#include <iostream>
#include <cmath>
#include <array>
#include "functions.h"

int main() {

const int N = 5; //# of harmonic oscillators in our heatbath
const int NTOTAL = N + 1;
const double GAMMA = 1.2; // expected superdiffusion exponent
const double BETA = 1.0; //kB*T
const double TSPAN[2] = {0, pow(10,1)};
const double DT = pow(10,-2);
double oscMass = 1.0; //mass of heaviest bath oscillator
double M = 1.0; // mass of distinguished particle
double omegaMin=pow(N,-0.7988), omegaMax=omegaMin*pow(N,1.0688); //highest and lowest eigenfrequency of the bath

//double omega[N];
std::array<double, N> omega;
std::array<double, N> masses;
std::array<double, NTOTAL> k;
std::array<double, NTOTAL> invM;
std::array<double, NTOTAL> q0;
std::array<double, NTOTAL> p0;
std::array<double, NTOTAL> q1;
std::array<double, NTOTAL> p1;


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

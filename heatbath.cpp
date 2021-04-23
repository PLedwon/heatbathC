//#include <cstdio>
#include <iostream>
#include <cmath>
#include <vector>
#include "functions.h"

int main() {

const int N = 500; //# of harmonic oscillators in our heatbath
const int NTOTAL = N + 1;
const double GAMMA = 1.2; // expected superdiffusion exponent
const double BETA = 1.0; //kB*T
const double TSPAN[2] = {0, pow(10,1)};
const double DT = pow(10,-2);
double oscMass = 1.0; //mass of heaviest bath oscillator
double M = 1.0; // mass of distinguished particle
double omegaMin=pow(N,-0.7988), omegaMax=omegaMin*pow(N,1.0688); //highest and lowest eigenfrequency of the bath
//double omega[N];
std::vector<double> omega(N,0);
std::vector<double> masses(NTOTAL,0);
std::vector<double> k(NTOTAL,0);
std::vector<double> invM(NTOTAL,0);
std::vector<double> q0(NTOTAL,0);
std::vector<double> p0(NTOTAL,0);
std::vector<double> q1(NTOTAL,0);
std::vector<double> p1(NTOTAL,0);


// setting the bathparameters
setEigenfrequencies(omega,omegaMin,omegaMax);
computeMasses(masses,oscMass,M,omega,omegaMin,GAMMA);
computeSpringConstants(k, masses, omega);
invertMasses(invM, M, masses);
// set initial conditions and solve the eom
generateInitialConditions(q0,p0,M,masses,k,BETA);

//solveEOM(q1,p1,q0,p0,invM,k,TSPAN,DT,N);

// solve EOM and save Q at times tsave
//save Q to file

//printArray_(masses,N);









}

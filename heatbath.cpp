//#include <cstdio>
#include <iostream>
#include <cmath>
#include <vector>
#include "functions.h"
using std::vector;

int main() {

const int N = 5; //# of harmonic oscillators in our heatbath
const int NTOTAL = N + 1; // adding the distinguished particle
const double GAMMA = 1.2; // expected superdiffusion exponent
const double BETA = 1.0; //kB*T
const double TSPAN[2] = {0, pow(10,0)};
const double DT = pow(10,-0);
double oscMass = 1.0; //mass of heaviest bath oscillator
double M = 1.0; // mass of distinguished particle
double omegaMin=pow(N,-0.7988), omegaMax=omegaMin*pow(N,1.0688); //highest and lowest eigenfrequency of the bath


// setting the bathparameters
const vector<double> omega = setEigenfrequencies(omegaMin,omegaMax,N);
const vector<double> masses = computeMasses(oscMass,M,omega,omegaMin,GAMMA);
const vector<double> k = computeSpringConstants(masses, omega);
const vector<double> invM = invertMasses(masses);
printVector_(invM);
struct heatbath bath;
// set initial conditions and solve the eom


try {
bath = generateInitialConditions(M, masses, k, BETA);
} catch (const char* msg) {
    std::cerr << msg << std::endl;
}
updatePositions(bath,k,DT,N);

printVector_(invM);



/*
printf("Energy before solving heatbath: %E \n", H(bath.q,bath.p,k,invM));
solveEOM(bath,invM,k,TSPAN,DT,N);
printf("Energy after solving heatbath: %E \n", H(bath.q,bath.p,k,invM));
printf("momentum after solving heatbath: %E \n", sum(bath.p)) ;
*/

// solve EOM and save Q at times tsave
//save Q to file
}
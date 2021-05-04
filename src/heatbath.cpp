//#include <cstdio>
#include <iostream>
#include <cmath>
#include <vector>
#include "functions.h"

const int N = 10000; //# of harmonic oscillators in our heatbath
//const int Heatbath::size = N + 1; // adding the distinguished particle

int main() {

const int NTOTAL = N + 1; // adding the distinguished particle
const double GAMMA = 1.2; // expected superdiffusion exponent
const double BETA = 1.0; //kB*T
const double TSPAN[2] = {0, pow(10,2)};
const double DT = pow(10,-4);
double oscMass = pow(10,1); //mass of heaviest bath oscillator
double M = pow(10,-3); // mass of distinguished particle
double omegaMin=pow(N,-0.7988), omegaMax=omegaMin*pow(N,1.0688); //highest and lowest eigenfrequency of the bath


// setting the bathparameters
double omega[N];
double masses[NTOTAL];
double k[NTOTAL];
double invM[NTOTAL];
setEigenfrequencies(omega,omegaMin,omegaMax,N);
computeMasses(masses,oscMass,M,omega,omegaMin,GAMMA,NTOTAL);
computeSpringConstants(k, masses, omega, NTOTAL);
invertMasses(invM,masses, NTOTAL);

// set initial conditions and solve the eom

double q0[NTOTAL];
double p0[NTOTAL];


try {
generateInitialConditions(q0, p0, M, masses, k, BETA, NTOTAL);
} catch (const char* msg) {
    std::cerr << msg << std::endl;
}
double initialEnergy = 0;
double initialMomentum = 0;
Heatbath bath(invM, k, q0, p0, initialEnergy, initialMomentum, NTOTAL);
printf("%e", bath.initialEnergy);
printArray_(bath.p,NTOTAL);

//solveEOM(bath,invM,k,TSPAN,DT,N);
//printf("absolute momentum error: %e \n", momentumError(bath));
//printf("rel. energy error = %e \n", energyError(bath,k,invM));
//
//write_csv("./data/trajectory.csv","trajectory", bath.trajectory);

}
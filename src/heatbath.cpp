#include <iostream>
#include <cmath>
#include <array>
#include "functions.h"
using std::array;
//const int Heatbath::size = N + 1; // adding the distinguished particle

/*
struct heatbath
{
//    static array<double, NTOTAL> q; // store phasespace coordinates of recent timestep

    static array<double, NTOTAL> q; // store phasespace coordinates of recent timestep
    static array<double, NTOTAL> p; // store phasespace coordinates of recent timestep
    static array<double, NTOTAL> invM;
    static array<double, NTOTAL> k;
    static array<double, NTIMESTEPS> trajectory; // save distinguished particle trajectory
    double initialEnergy;
    double initialMomentum;
    heatbath(array<double, NTOTAL> q, array<double, NTOTAL> p, array<double, NTOTAL> invM, array<double, NTOTAL> k,  array<double, NTIMESTEPS> trajectory, double initialEnergy, double initialMomentum);
};

heatbath::heatbath(array<double, NTOTAL> q, array<double, NTOTAL> p, array<double, NTOTAL> invM,
                   array<double, NTOTAL> k, array<double, NTIMESTEPS> trajectory, double initialEnergy,
                   double initialMomentum) {
    q = q;
    p = p;
    invM = invM;
    k = k;
    printArray_(k);
    trajectory = trajectory;
    initialEnergy = initialEnergy;
    initialMomentum = initialMomentum;
}
*/

int main() {

const int N = 100; //# of harmonic oscillators in our heatbath
const int NTOTAL = N + 1; // adding the distinguished particle
const double TSPAN[2] = {0, pow(10,3)};
const double DT = pow(10,-5);
const int NTIMESTEPS = ceil((TSPAN[1]-TSPAN[0])/DT);
const double GAMMA = 1.2; // expected superdiffusion exponent
const double BETA = 1.0; //kB*T
const int NSAVE = (int) fmin(pow(10,6),NTIMESTEPS); // max outfile size capped at about 10 MB

double oscMass = pow(10,3); //mass of heaviest bath oscillator
double M = pow(10,-5); // mass of distinguished particle
double omegaMin=pow(N,-0.7988), omegaMax=omegaMin*pow(N,1.0688); //highest and lowest eigenfrequency of the bath
// setting the bathparameters
static array<double, N> omega;
static array<double, NTOTAL> masses;
static array<double, NTOTAL> kTemp;
static array<double, NTOTAL> invMTemp;
static array<double, NTOTAL> q;
static array<double, NTOTAL> p;
static array <double, NSAVE> trajectory;

setEigenfrequencies(omega,omegaMin,omegaMax);
computeMasses(masses,oscMass,M,omega,omegaMin,GAMMA);
computeSpringConstants(kTemp, masses, omega);
invertMasses(invMTemp,masses);

const static array<double, NTOTAL> k = kTemp;
const static array<double, NTOTAL> invM = invMTemp;

try {
generateInitialConditions(q, p, M, masses, k, BETA);
} catch (const char* msg) {
    std::cerr << msg << std::endl;
}

const double initialEnergy = H(q,p,invM,k);
const double initialMomentum = sum(p);

    /*
Heatbath bath(invM, k, q0, p0, initialEnergy, initialMomentum, NTOTAL);
printf("%e", bath.initialEnergy);
printArray_(bath.p,NTOTAL);
*/

solveEOM(q,p,invM,k,trajectory,DT,NTIMESTEPS);
printf("absolute momentum error: %e \n", momentumError(p, initialMomentum ));
printf("rel. energy error = %e \n", energyError(q,p,invM,k,initialEnergy));
//


int name = rand() % (int) pow(10,6); // generate random name for output file

write_csv("./data/trajectory" + std::to_string(name) + ".csv","trajectory" , trajectory);

}
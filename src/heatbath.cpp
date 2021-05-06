#include <iostream>
#include <cmath>
#include <array>
#include "functions.h"
using std::array;
//const int Heatbath::size = N + 1; // adding the distinguished particle

const int N = 10; //# of harmonic oscillators in our heatbath
const int NTOTAL = N + 1; // adding the distinguished particle
const double GAMMA = 1.2; // expected superdiffusion exponent
const double BETA = 1.0; //kB*T
const double TSPAN[2] = {0, pow(10,2)};
const double DT = pow(10,-4);
const int NTIMESTEPS = (TSPAN[1]-TSPAN[0])/DT;

struct heatbath
{
//    static array<double, NTOTAL> q; // store phasespace coordinates of recent timestep
    array<double, NTOTAL> q; // store phasespace coordinates of recent timestep
    array<double, NTOTAL> p; // store phasespace coordinates of recent timestep
    array<double, NTOTAL> invM;
    array<double, NTOTAL> k;
    array<double, NTIMESTEPS> trajectory; // save distinguished particle trajectory
    double initialEnergy;
    double initialMomentum;
    heatbath(array<double, NTOTAL> q, array<double, NTOTAL> p, array<double, NTOTAL> invM, array<double, NTOTAL> k,  array<double, NTIMESTEPS> trajectory, double initialEnergy, double initialMomentum) {
        q = q;
        p = p;
        invM = invM;
        k = k;
        trajectory = trajectory;
        initialEnergy = initialEnergy;
        initialMomentum = initialMomentum;

    }
};

int main() {

double oscMass = pow(10,1); //mass of heaviest bath oscillator
double M = pow(10,-3); // mass of distinguished particle
double omegaMin=pow(N,-0.7988), omegaMax=omegaMin*pow(N,1.0688); //highest and lowest eigenfrequency of the bath
// setting the bathparameters
array<double, N> omega;
static array<double, NTOTAL> masses;
static array<double, NTOTAL> k;
static array<double, NTOTAL> invM;
static array<double, NTOTAL> q;
static array<double, NTOTAL> p;
static array <double, NTIMESTEPS> trajectory;
double initialEnergy = 0;
double initialMomentum = 0;

setEigenfrequencies(omega,omegaMin,omegaMax);
computeMasses(masses,oscMass,M,omega,omegaMin,GAMMA);
computeSpringConstants(k, masses, omega);
    printArray_(k);
invertMasses(invM,masses);
    printf("%e ",invM[0]);
// set initial conditions and solve the eom
heatbath bath(q, p, invM, k, trajectory, initialEnergy, initialMomentum);
printArray_(bath.k);
/*


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
*/
//solveEOM(bath,invM,k,TSPAN,DT,N);
//printf("absolute momentum error: %e \n", momentumError(bath));
//printf("rel. energy error = %e \n", energyError(bath,k,invM));
//
//write_csv("./data/trajectory.csv","trajectory", bath.trajectory);

}
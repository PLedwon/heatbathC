#include <iostream>
#include <cmath>
#include <array>
#include <random>
#include <ctime>
#include <fstream>
#include "functions.h"
using std::array;


int main() {

const int N = 20000; //# of harmonic oscillators in our heatbath
const int NTOTAL = N + 1; // adding the distinguished particle
const double TSPAN[2] = {0, pow(10,-1)};
const double DT = pow(10,-6);
const int NTIMESTEPS = ceil((TSPAN[1]-TSPAN[0])/DT);
const double GAMMA = 1.2; // expected superdiffusion exponent
const double BETA = 1.0; //kB*T
const int NSAVE = (int) fmin(pow(10,6),NTIMESTEPS); // max outfile size capped at about 10 MB

double oscMass = pow(10,2); //mass of heaviest bath oscillator
double M = pow(10,-3); // mass of distinguished particle
double omegaMin=pow(N,-0.8323), omegaMax=omegaMin*pow(N,1.05764); //highest and lowest eigenfrequency of the bath

// setting the bathparameters
static array<double, N> omega;
static array<double, NTOTAL> masses;
static array<double, NTOTAL> kTemp;
static array<double, NTOTAL> invMTemp;
static array<double, NTOTAL> q;
static array<double, NTOTAL> p;
static array <double, NSAVE> trajectory;
static array <double, NSAVE> energyErrArray;
static array <double, NSAVE> momentumErrorArray;

setEigenfrequencies(omega,omegaMin,omegaMax);
computeMasses(masses,oscMass,M,omega,omegaMin,GAMMA);
computeSpringConstants(kTemp, masses, omega);
invertMasses(invMTemp,masses);

const static array<double, NTOTAL> k = kTemp;
const static array<double, NTOTAL> invM = invMTemp;

Heatbath<static array, NTOTAL, NTIMESTEPS >(NTOTAL, k, invM);


try {
generateInitialConditions(q, p, M, masses, k, BETA);
} catch (const char* msg) {
    std::cerr << msg << std::endl;
}

const double initialEnergy = H(q,p,invM,k);
const double initialMomentum = sum(p);

time_t begin,end; // save runtime

time(&begin);

solveEOM(q,p,invM,k,trajectory,energyErrArray,momentumErrorArray,DT,NTIMESTEPS);
double mError= momentumError(p, initialMomentum );
double eError= energyError(q,p,invM,k,initialEnergy);

time(&end);
double difference = difftime(end,begin)/3600.0;


//generate a random name for outputfiles, write to logfile
std::random_device rd;
std::uniform_int_distribution<int> dist(0, 999999);
int name = dist(rd);
write_csv("./data/trajectory" + std::to_string(name) + ".csv","trajectory" , trajectory);

std::string filename("./data/log/logfile.txt");
std::fstream file;
file.open(filename, std::ios_base::app | std::ios_base::in);
if (file.is_open())
	file << std::to_string(name) << ", " << std::to_string(eError) << ", " << std::to_string(mError) << ", " << std::to_string(difference) << std::endl;
file.close();

return 0;
}

#include <iostream>
#include <cmath>
#include <array>
#include <random>
#include <ctime>
#include <cstring>
#include <fstream>
#include <algorithm>
#include "functions.h"
using std::array;

constexpr int N = 4000; //# of harmonic oscillators in our heatbath
constexpr int NTOTAL = N + 1; // adding the distinguished particle
constexpr double TSPAN[2] = {0, pow(10,3)};
constexpr double DT =1*pow(10,-6);
const long long NTIMESTEPS = ceil((TSPAN[1]-TSPAN[0])/DT);
const double GAMMA = 1.2; // expected superdiffusion exponent
const double BETA = 1.0; //kB*T
constexpr int NSAVE = (int)fmin(pow(10,4),NTIMESTEPS); // max outfile size capped at about 10 MB

//initialize static heatbath members
double Heatbath::k[NTOTAL] = {0};
double Heatbath::invM[NTOTAL] = {0};
double Heatbath::q[NTOTAL] = {0};
double Heatbath::p[NTOTAL] = {0};
double Heatbath::trajectory[NSAVE] = {0};
double Heatbath::energyErr[NSAVE] = {0};
double Heatbath::momentumErr[NSAVE] = {0};
int Heatbath::size = NTOTAL;
int Heatbath::nSave = NSAVE;
double Heatbath::initialEnergy;
double Heatbath::initialMomentum;


int main() {
double oscMass = pow(10,2); //mass of heaviest bath oscillator
double M = pow(10,-3); // mass of distinguished particle
double omegaMin=pow(N,-0.7988), omegaMax=omegaMin*pow(N,1.0688); //highest and lowest eigenfrequency of the bath

double omega[N];
double masses[NTOTAL];

//setting bath parameters
setEigenfrequencies(omega,omegaMin,omegaMax);
computeMasses(masses,oscMass,M,omega,omegaMin,GAMMA);
computeSpringConstants(Heatbath::k, masses, omega);
invertMasses(Heatbath::invM,masses);

Heatbath bath;

try {
generateInitialConditions(bath, M, masses, BETA);
} catch (const char* msg) {
    std::cerr << msg << std::endl;
}

time_t begin,end; // save runtime
time(&begin);

try {
    solveEOM(bath, DT, NTIMESTEPS);
} catch (const char* msg) {
    std::cerr << msg << std::endl;
}
time(&end);

double maxEnergyErr = *std::max_element(bath.energyErr,bath.energyErr+bath.size);
double maxMomentumErr = *std::max_element(bath.momentumErr,bath.momentumErr+bath.size);
double difference = difftime(end,begin)/3600.0;

//generate a random name for outputfiles, write to logfile
std::random_device rd;
std::uniform_int_distribution<int> dist(0, 999999);
int name = dist(rd);

//save trajectory of distinguished particle to file
write_csv("./data/trajectory" + std::to_string(name) + ".csv","trajectory" , bath);


//save errors and runtime to logfile
std::string filename("./data/log/logfile.txt");
std::fstream file;
file.open(filename, std::ios_base::app | std::ios_base::in);
if (file.is_open())
    file.precision(8);
	file << std::scientific << std::to_string(name) << ", " << maxEnergyErr << ", " << maxMomentumErr << ", " << difference << std::endl;
file.close();

return 0;
}

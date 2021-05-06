#ifndef functions
#define functions
#include <cmath>
#include <iostream>
#include <random>
#include <fstream>
#include <array>
using std::array;

template<typename T,size_t n>
void printArray_(array<T,n> a) {
    int i;
    for (i = a.front(); i < a.size(); i++) {
        //std::cout << a[i] << ' ';
        printf("%e ", a[i]);
    }
    std::cout << '\n';
}

template<typename T,size_t n>
double avg(array<T,n> a) {
    double sum=0;
    for (int i = a.front(); i < a.size(); ++i) {
       sum += a[i];
    }
    return (double) sum/a.size();
}

template<typename T,size_t n>
void computeMasses(array<T,n> (&masses), double oscMass, double M, array<T, n-1> omega, double omegaMin, const double GAMMA){
  masses[0]=M;
  for (int i = 1; i < masses.size() ; i++) {
    masses[i]=oscMass*pow((omega[i-1]/omegaMin),(GAMMA-3))*exp(-omega[i-1]);
  }

}

template<typename T,size_t n>
void computeSpringConstants(array<T,n> (&k), array<T,n> masses, array<T,n-1> omega) {
  k[0] = 0;
  for (int i = 1; i < k.size(); i++) {
    k[i]=masses[i]*pow(omega[i-1],2);
  }
}
/*
double H(heatbath bath, const vector<double> k, const vector<double> invM) { // compute total energy of the system
  double E = bath.p[0] * bath.p[0] * invM[0];
  for (int i = 1; i < bath.q.size() ; ++i) {
    E += bath.p[i] * bath.p[i] * invM[i] + k[i] * pow(bath.q[i] - bath.q[0], 2);
    //printf("accumulated energy %e \n", E);
    }
  E *= 0.5;
  return  E;
}
*/
template<typename T,size_t n>
double sum(array<T,n> p) {
    double mom = 0;
    for (int i = p.front(); i < p.size(); ++i) {
        mom += p[i];
    }
    return mom;
}

template<typename T, size_t n>
void setEigenfrequencies(array<T,n> (&omega), double omegaMin, double omegaMax) {
    double c;
    c = (omegaMax - omegaMin)/(omega.size()-1);
    for(int i = 0; i < omega.size() ; ++i) // equidistant distribution of eigenfrequencies of the harmonic oscillators
        omega[i] = omegaMin + i*c;
    //omega.back() = omegaMax;
}

template<typename T, size_t n>
void invertMasses(array<T,n> (&invM), array<T,n> masses) {
    for (int i = 0; i < invM.size(); ++i) {
       invM[i] = 1/masses[i];
    }
}
/*
void generateInitialConditions(double (&q0)[], double (&p0)[], double M, double masses[], double k[], const double BETA, const int NTOTAL) {

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{0,1};
    double pref=pow(BETA,-0.5);

    //set the initial conditions for the distinguished particle
    q0[0] = 0;
    p0[0] = pref * pow(M,0.5) * d(gen);

    for (int i = 1; i < NTOTAL; ++i) {
       q0[i] = q0[0] + pref*pow(k[i],-0.5)*d(gen);
       p0[i] = pref*pow(masses[i],0.5)* d(gen);
    }
    //initialize heatbath with vanishing center of mass velocity
    double avgMomentum = avg(p0);
    double p0sum = sum(p0,NTOTAL);
    for (int i = 0; i < NTOTAL; ++i) {
        p0[i] -= avgMomentum;
    }
    p0sum = sum(p0,NTOTAL);

    //check that total momentum is close to zero
    if (std::abs(p0sum) > pow(10,-13)) {
        throw "Error: CoM velocity is not 0 while initializing heatbath";
    }
}
*/

 /*
void updateMomenta(heatbath &bath, vector<double> k,const double DT,const int N) {
    double s;
    for (int i = 1; i < N+1; ++i) {
        s = k[i] * (bath.q[0]-bath.q[i])*DT;
        //printf("s=%e \n", s);
        bath.p[0] -= s;
       // printf("updated P = %e \n", bath.p[0]);
        bath.p[i] += s;
    }
}// have to be of same length

void updatePositions(heatbath &bath, const vector<double> invM, const double DT, const int N) {
    for (int i = 0; i < N+1 ; ++i) {
       bath.q[i] += bath.p[i]*invM[i]*DT;
    }
}

void makeTimestep(heatbath &bath, vector<double> k, vector<double> invM,const double DT,const int N) {
    updateMomenta(bath,k,DT,N); //update momenta first for a symplectic Euler algorithm
    updatePositions(bath,invM,DT,N);

}

void solveEOM(heatbath &bath, vector<double> invM, vector<double> k, const double TSPAN[], const double DT, const int N) {
    int nTimesteps = ceil((TSPAN[1]-TSPAN[0])/DT);
    double initialEnergy = H(bath,k,invM);
    bath.initialEnergy = initialEnergy;
    for (int i = 0; i < nTimesteps ; ++i) {
        makeTimestep(bath,k,invM,DT,N);
        bath.trajectory.push_back(bath.q[0]); //  save most recent position

    }

}

double energyError(heatbath bath, vector<double> k, vector<double> invM){
 return (H(bath,k,invM)-bath.initialEnergy)/bath.initialEnergy;
}

double momentumError(Heatbath Bath) {
    return Bath.initialMomentum-sum(Bath.p,Bath.size);
}
*/
//////////////////////////////////////////////////////////
void write_csv(std::string filename, std::string colname, std::vector<double> vals){
    std::ofstream myFile(filename);
    myFile << colname << "\n";

    // Send data to the stream
    for(int i = 0; i < vals.size(); ++i)
    {
        myFile << vals.at(i) << "\n";
    }

    myFile.close();
}


#endif

#ifndef functions
#define functions
#include <cmath>
#include <random>
#include <array>

void printArray_(double a[], int len) {
    int i;
    for (i = 0; i < len; i++) std::cout << a[i] << ' ';
    printf("\n")
}

double avg(double a[], int N) {
    double avg=0;
    for (int i = 0; i < N + 1; ++i) {
       avg+=a[i]/ (double) (N+1) ;
    }
//avg /= (double) N+1;
    return avg;
}




void computeMasses(double masses[], double oscMass, double  omega[], double omegaMin,const double GAMMA, const int N){
//  int n = nElems(masses);
  for (int i = 0; i < N ; i++) {
    masses[i]=oscMass*pow((omega[i]/omegaMin),(GAMMA-3))*exp(-omega[i]);
  }
}

void computeSpringConstants(double k[] ,double masses[], double omega[],const int N) {
  //int n = nElems(k);
  k[0]=0.0;
  for (int i = 1; i < N+1; i++) {
    k[i]=masses[i]*pow(omega[i],2);
  }
}

double H(double q[] , double p[], double k[], double invM[], const int N) { // compute total energy of the system
  double E = p[0] * p[0] * invM[0];
  for (int i = 1; i < N+1 ; ++i) {
    E += p[i] * p[i] * invM[i] + k[i] * pow(q[i] - q[0], 2);
    }
  E *= 0.5;
  return  E;
}

double sum(double p[], int N) {
    double mom = 0;
    for (int i = 0; i < N+1 ; ++i) {
        mom = mom  + p[i];
    }
    return mom;
}

void setEigenfrequencies(double omega[], double omegaMin, double omegaMax, const int NTOTAL) {
    double c;
    c = (omegaMax - omegaMin)/(NTOTAL-1);
    for(int i = 0; i < NTOTAL - 1; ++i) // equidistant distribution of eigenfrequencies of the harmonic oscillators
        omega[i] = omegaMin + i*c;
    omega[NTOTAL] = omegaMax;
    printArray_(omega);
}

void invertMasses(double invM[], double M, double masses[], const int N) {
    invM[0]=1/M;
    for (int i = 0; i < N ; ++i) {
       invM[i+1] = 1/masses[i];
    }
}

void generateInitialConditions(double q0[], double p0[], double M, double masses[], double k[],const double BETA ,const int N) {

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{0,1};
    double pref=pow(BETA,-0.5);

    //set the initial conditions for the distinguished particle
    q0[0] = 0;
    p0[0] = pref * pow(M,0.5) * d(gen);

    for (int i = 1; i < N + 1; ++i) {
       q0[i] = q0[0] + pref*pow(k[i],-0.5)*d(gen);
       p0[i] = pref*pow(masses[i-1],0.5)* d(gen);
    }
    //initialize heatbath with vanishing center of mass velocity
    double avgMomentum = avg(p0,N);
    double p0sum = sum(p0,N);
    printf("total momentum: %E \n", p0sum);
    printArray_(p0,N);
    for (int i = 0; i < N + 1; ++i) {
        p0[i] = p0[i] - avgMomentum;
    }
    p0sum = sum(p0,N+1);
   // printf("avg momentum: %E \n", avg(p0,N));
   // printf("total momentum: %E \n", p0sum);
    //printArray_(p0,N);
}

void updateMomenta(double p1[], double q0[], double p0[], double k[],const double DT,const int N) {
    double s;
    p1[0] = p0[0];
    for (int i = 1; i < N+1; ++i) {
        s = k[i] * (q0[0]-q0[i])*DT;
        p1[0] -= s;
        p1[i] = s + p0[i];
    }
}// have to be of same length

void updatePositions(double q1[], double p1[], double q0[], double invM[], const double DT, const int N) {
    for (int i = 0; i < N+1 ; ++i) {
       q1[i] = q0[i] + p1[i]*invM[i]*DT;
    }
}

void makeTimestep(double q1[], double p1[], double q0[],double p0[], double k[], double invM[],const double DT,const int N) {
    updateMomenta(p1,q0,p0,k,DT,N); //update momenta first for a symplectic Euler algorithm
    updatePositions(q1,p1,q0,invM,DT,N);
}

void solveEOM(double q1[], double p1[], double q0[], double p0[], double invM[], double k[], const double TSPAN[], const double DT, const int N) {
    int nTimesteps = ceil((TSPAN[1]-TSPAN[0])/DT);
    double initialEnergy;
    for (int i = 0; i < nTimesteps ; ++i) {
        makeTimestep(q1,p1,q0,p0,k,invM,DT,N);
        q0=q1; //updated coordinates replace old ones for next timestep
        p0=p1;
    }

}

//////////////////////////////////////////////////////////


#endif

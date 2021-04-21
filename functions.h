#ifndef functions
#define functions
#include <cmath>
#include <random>

#define nElems(x) (sizeof(x) / sizeof(x[0]))


void computeMasses(double masses[], double oscMass, double  omega[], double omegaMin, double gamma, int N){
//  int n = nElems(masses);
  for (int i = 0; i < N ; i++) {
    masses[i]=oscMass*pow((omega[i]/omegaMin),(gamma-3))*exp(-omega[i]);
  }
}

void computeSpringConstants(double k[] ,double masses[], double omega[], int N) {
  //int n = nElems(k);
  k[0]=0.0;
  for (int i = 1; i < N+1; i++) {
    k[i]=masses[i]*pow(omega[i],2);
  }
}

double H(double q[] , double p[], double k[], double invM[], int N) { // compute total energy of the system
  double E = p[0] * p[0] * invM[0];
  for (int i = 1; i < N+1 ; ++i) {
    E += p[i] * p[i] * invM[i] + k[i] * pow(q[i] - q[0], 2);
    }
  E *= 0.5;
  return  E;
}

double totalMomentum(double p[], int N) {
    double mom = 0;
    for (int i = 1; i < N+1 ; ++i) {
        mom += p[i];
    }
}

void computeOmega(double omega[], double omegaMin, double omegaMax, int N) {
    double c;
    c = (omegaMax - omegaMin)/(N - 1);
    for(int i = 0; i < N - 1; ++i) // equidistant distribution of eigenfrequencies of the harmonic oscillators
        omega[i] = omegaMin + i*c;
    omega[N - 1] = omegaMax;
}

void invertMasses(double invM[], double M, double masses[], int N) {
    invM[0]=1/M;
    for (int i = 0; i < N ; ++i) {
       invM[i+1] = 1/masses[i];
    }
}

void generateInitialConditions(double q0[], double p0[], double M, double masses[],double k[], double beta ,int N) {

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{0,1};
    double pref=pow(beta,-0.5);

    //set the initial conditions for the distinguished particle
    q0[0]=0;
    p0[0]=pref*pow(M,0.5)*d(gen);

    for (int i = 1; i < N+1; ++i) {
       q0[i+1] = q0[0] + pref*pow(k[i],-0.5)*d(gen);
       p0[i+1] = pref*pow(masses[i-1],0.5)* d(gen);
    }
}

void updateMomenta(double p1[], double q0[], double p0[], double k[], double dt, int N) {
    double s;
    p1[0] = p0[0];
    for (int i = 1; i < N+1; ++i) {
        s = k[i] * (q0[0]-q0[i])*dt;
        p1[0] -= s;
        p1[i] = s + p0[i];
    }
}// have to be of same length

void updatePositions(double q1[], double p1[], double q0[], double invM[], double dt, int N) {
    for (int i = 0; i < N+1 ; ++i) {
       q1[i] = q0[i] + p1[i]*invM[i]*dt;
    }
}



void makeTimestep(double q1[], double p1[], double q0[],double p0[], double k[], double invM[], double dt, int N) {
    updateMomenta(p1,q0,p0,k,dt,N); //update momenta first for a symplectic Euler algorithm
    updatePositions()
}

void solveEOM(double q1[], double p1[], double q0[], double p0[], double invM[], double k[], double tspan[], double dt, int N) {
    int nTimesteps = ceil((tspan[1]-tspan[0])/dt);
    for (int i = 0; i < nTimesteps ; ++i) {
        makeTimestep(q1,p1,q0,p0,k,invM,dt,N);
        q0=q1;
        p0=p1;
    }

}

//////////////////////////////////////////////////////////





void printArray_(double a[], int len) {
  int i;
    for (i = 0; i < len; i++) std::cout << a[i] << ' ';
}








#endif

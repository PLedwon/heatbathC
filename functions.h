#ifndef functions
#define functions
#include <math.h>

//#define nElems(x) (sizeof(x) / sizeof(x[0]))


void computeMasses(double masses[], double oscMass, double  omega[], double omegaMin, double gamma, int N){
//  int n = nElems(masses);
  for (int i = 0; i < N ; i++) {
    masses[i]=oscMass*pow((omega[i]/omegaMin),(gamma-3))*exp(-omega[i]);
  }
}

void computeSpringConstants(double k[] ,double masses[], double omega[], int N) {
  //int n = nElems(k);
  for (int i = 0; i < N; i++) {
    k[i]=masses[i]*pow(omega[i],2);
  }
}

double H(double q[] , double p[], double k[], double invM[], int N) {
  double E = p[0] * p[0] * invM[0];
  for (int i = 1; i < N+1 ; ++i) {
    E += p[i] * p[i] * invM[i] + k[i - 1] * pow(q[i] - q[0], 2);
    }

  E *= 0.5;
  return  E;
}

void computeOmega(double omega[], double omegaMin, double omegaMax, int N) {
    double c;
    c = (omegaMax - omegaMin)/(N - 1);
    for(int i = 0; i < N - 1; ++i)
        omega[i] = omegaMin + i*c;
    omega[N - 1] = omegaMax;
}

void invertMasses(double invM[], double M, double masses[], int N) {
    invM[0]=1/M;
    for (int i = 0; i < N ; ++i) {
       invM[i+1] = 1/masses[i];
    }
}

//void generateInitialConditions(double q0[], double p0[], )

//////////////////////////////////////////////////////////


/*
void printArray_(double a[], int len) {
  int i;
    for (i = 0; i < len; i++) std::cout << a[i];
}
*/







#endif

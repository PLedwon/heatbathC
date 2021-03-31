#ifndef functions
#define functions
#include <math.h>

#define nElems(x) (sizeof(x) / sizeof(x[0]))


void computeMasses(double masses[], double oscMass, double  omega[], double omegaMin, double gamma){
  int n = nElems(masses);
  for (int i = 0; i < n ; i++) {
    masses[i]=oscMass*pow((omega[i]/omegaMin),(gamma-3))*exp(-omega[i]);
  }
}

void computeSpringConstants(double k[] ,double masses[], double omega[]) {
  int n = nElems(k);
  for (int i = 0; i < n; i++) {
    k[i]=masses[i]*pow(omega[i],2);
  }
}

double H(double q[] ,double p[],double k[],double mInv[]) {
  return q[0];
}

void computeOmega(double omega[], double omegaMin, double omegaMax, int N) {
    double c;
    c = (omegaMax - omegaMin)/(N - 1);
    for(int i = 0; i < N - 1; ++i)
        omega[i] = omegaMin + i*c;
    omega[N - 1] = omegaMax;
}





//////////////////////////////////////////////////////////

double* linspace(double a, double b, int n, double u[])
{
    double c;
    int i;

    if(n < 2 || u == 0)
        return (void*)0;

    c = (b - a)/(n - 1);

    for(i = 0; i < n - 1; ++i)
        u[i] = a + i*c;

    u[n - 1] = b;

    return u;
}

void printArray_(double a[]) {
  int i;
    for (i = 0; i < nElems(a); i++) printf("%f ", a[4]);
}








#endif

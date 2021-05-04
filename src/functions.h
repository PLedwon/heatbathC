#ifndef functions
#define functions
#include <cmath>
#include <random>
#include <vector>
#include <fstream>
using std::vector;

/*
struct heatbath
{
    vector<double> q; // store phasespace coordinates of recent timestep
    vector<double> p;
    vector<double> trajectory = {0}; // position of distinguished particle in time
    double initialEnergy;
    double initialMomentum;

};
*/





















void printArray_(double a[], int n) {
    int i;
    //for (i = 0; i < a.size(); i++) std::cout << a[i] << ' ';
    for (i = 0; i < n; i++) printf("%e ", a[i]);
    printf("\n");
}

double avg(double a[], int n) {
    double sum=0;
    for (int i = 0; i < n; ++i) {
       sum += a[i];
    }
    return (double) sum/n;
}

void computeMasses(double (&masses)[], double oscMass, double M, double omega[], double omegaMin, const double GAMMA, int NTOTAL){
  masses[0]=M;
  for (int i = 1; i < NTOTAL ; i++) {
    masses[i]=oscMass*pow((omega[i-1]/omegaMin),(GAMMA-3))*exp(-omega[i-1]);
  }

}

void computeSpringConstants(double (&k)[], double masses[], double omega[], int NTOTAL) {
  k[0] = 0;
  for (int i = 1; i < NTOTAL; i++) {
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
double sum(double p[], int n) {
    double mom = 0;
    for (int i = 0; i < n; ++i) {
        mom += p[i];
    }
    return mom;
}

void setEigenfrequencies(double (&omega)[], double omegaMin, double omegaMax, int N) {
    double c;
    c = (omegaMax - omegaMin)/(N-1);
    for(int i = 0; i < N ; ++i) // equidistant distribution of eigenfrequencies of the harmonic oscillators
        omega[i] = omegaMin + i*c;
    //omega.back() = omegaMax;
}

void invertMasses(double (&invM)[], double masses[], int NTOTAL) {
    for (int i = 0; i < NTOTAL; ++i) {
       invM[i] = 1/masses[i];
    }
}

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
    double avgMomentum = avg(p0,NTOTAL);
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

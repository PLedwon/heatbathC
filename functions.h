#ifndef functions
#define functions
#include <cmath>
#include <random>
#include <vector>

void printVector_(std::vector<double> a) {
    int i;
    for (i = 0; i < a.size(); i++) std::cout << a[i] << ' ';
    printf("\n");
}

double avg(std::vector<double> a) {
    double average=std::accumulate(a.begin(),a.end(),0.0);
    return (double) average/a.size();
}




void computeMasses(std::vector<double> masses, double oscMass,double M, std::vector<double>  omega, double omegaMin,const double GAMMA){
  masses.front()=M;
  for (int i = 1; i < masses.size() ; i++) {
    masses[i]=oscMass*pow((omega[i]/omegaMin),(GAMMA-3))*exp(-omega[i]);
  }
}

void computeSpringConstants(std::vector<double> k ,std::vector<double> masses, std::vector<double> omega) {
  //int n = nElems(k);
  k[0]=0.0;
  for (int i = 1; i < k.size(); i++) {
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

double sum(std::vector<double> p) {
    double mom = std::accumulate(p.begin(),p.end(),0.0);
    return mom;
}

void setEigenfrequencies(std::vector<double> omega, double omegaMin, double omegaMax) {
    double c;
    c = (omegaMax - omegaMin)/(omega.size()-1);
    for(int i = 0; i < omega.size() - 1; ++i) // equidistant distribution of eigenfrequencies of the harmonic oscillators
        omega[i] = omegaMin + i*c;
    omega.back() = omegaMax;
    printVector_(omega);
}

void invertMasses(std::vector<double> invM, double M, std::vector<double> masses) {
    for (int i = 0; i < invM.size() ; ++i) {
       invM[i] = 1/masses[i];
    }
}

void generateInitialConditions(std::vector<double> q0, std::vector<double> p0, double M, std::vector<double> masses, std::vector<double> k,const double BETA) {

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{0,1};
    double pref=pow(BETA,-0.5);

    //set the initial conditions for the distinguished particle
    q0[0] = 0;
    p0[0] = pref * pow(M,0.5) * d(gen);

    for (int i = 1; i < q0.size(); ++i) {
       q0[i] = q0[0] + pref*pow(k[i],-0.5)*d(gen);
       p0[i] = pref*pow(masses[i],0.5)* d(gen);
       printf("rnd momentum: %E \n", p0[i]);
    }
    //initialize heatbath with vanishing center of mass velocity
    double avgMomentum = avg(p0);
    double p0sum = sum(p0);
    printf("total momentum before: %E \n", p0sum);
    printVector_(p0);
    for (int i = 0; i < p0.size(); ++i) {
        p0[i] -= avgMomentum;
    }
    p0sum = sum(p0);
    printf("avg momentum before shift: %E \n", avgMomentum);
    printf("total momentum after shift: %E \n", p0sum);
    printVector_(p0);
}

void updateMomenta(std::vector<double> p1, std::vector<double> q0, std::vector<double> p0, std::vector<double> k,const double DT,const int N) {
    double s;
    p1[0] = p0[0];
    for (int i = 1; i < N+1; ++i) {
        s = k[i] * (q0[0]-q0[i])*DT;
        p1[0] -= s;
        p1[i] = s + p0[i];
    }
}// have to be of same length

void updatePositions(std::vector<double> q1, std::vector<double> p1, std::vector<double> q0, std::vector<double> invM, const double DT, const int N) {
    for (int i = 0; i < N+1 ; ++i) {
       q1[i] = q0[i] + p1[i]*invM[i]*DT;
    }
}

void makeTimestep(std::vector<double> q1, std::vector<double> p1, std::vector<double> q0,std::vector<double> p0, std::vector<double> k, std::vector<double> invM,const double DT,const int N) {
    updateMomenta(p1,q0,p0,k,DT,N); //update momenta first for a symplectic Euler algorithm
    updatePositions(q1,p1,q0,invM,DT,N);
}

void solveEOM(std::vector<double> q1, std::vector<double> p1, std::vector<double> q0, std::vector<double> p0, std::vector<double> invM, std::vector<double> k, const double TSPAN[], const double DT, const int N) {
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

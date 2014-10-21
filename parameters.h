#ifndef PARAMETERS_H
#define PARAMETERS_H

#include<string>
using namespace std;

const double hbarC=0.197327053;  //GeV*fm

const int Maxparticle=400;            //size of array for storage of the particles
const int Maxdecaychannel=13;
const int Maxdecaypart=5;

const int eta_s_npts = 40;
const double eta_s_i = 0.0;
const double eta_s_f = 4.0;
const int qnpts = 20;
const double delta_q = 0.01;
const double init_q = delta_q/2.0;

const double n_localp_T = 5;
const double localp_T_min = 0.6;
const double localp_T_max = 1.0;
const double n_localp_phi = 5;
const double localp_phi_min = 0.0;
const double localp_phi_max = 2*M_PI;


const string path = "results";

const double tol = 1e-15;  //tolarence
const int flagneg = 1;     //neglect all points that are negative

const int MCint_calls = 5000;  //# of calls for monte carlo integration

const size_t fit_max_iterations = 1000;  // stop at this point if not converged 
const double fit_tolarence = 1e-6;

#endif

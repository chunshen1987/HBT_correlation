#ifndef READINDATA_H
#define READINDATA_H

#include "parameters.h"
#include<fstream>
using namespace std;

typedef struct
{
   int IEOS;            //EOS selector
   double b;            //Impact parameter
   double tau_0;        //thermalization time
   double e_dec;        //kinetic decoupling energy density
   double T_dec;
   double eta_s;        //specific shear viscosity
   double zeta_s;       //specific bulk viscosity
   int NXD;             //number of skipped points along x direction
   int NYD;             //number of skipped points along y direction
   int NTauD;           //number of skipped points along tau direction
   double dtau;         //lattice spacing in tau direction
   double dx;           //lattice spacing in x direction
   double dy;           //lattice spacing in y direction
   int NLS;             //number of lattice size (assume square lattice NX=NY)
}hydropara;

typedef struct
{
   double tau, xpt, ypt;
   double da0, da1, da2;
   double vx, vy;
   double Edec, Tdec, Pdec;
   double Bn, muB, muS;
   double pi00, pi01, pi02, pi11, pi12, pi22, pi33;
   double particle_mu[Maxparticle];
}FO_surf;

typedef struct
{
  int monval;			// Montecarlo number according PDG
  string name;
  double mass;
  double width;
  int gspin;			// spin degeneracy
  int baryon;
  int strange;
  int charm;
  int bottom;
  int gisospin;			// isospin degeneracy
  int charge;
  int decays;			// amount of decays listed for this resonance
  int stable;			// defines whether this particle is considered as stable
  int decays_Npart[Maxdecaychannel];
  double decays_branchratio[Maxdecaychannel];
  int decays_part[Maxdecaychannel][Maxdecaypart];
  double mu;
  int sing;       //Bose-Einstein or Dirac-Fermi statistics
}particle_info;

void read_hydropar(hydropara* hp);
int get_filelength(string filepath);
void read_decdat(int length, FO_surf* surf_ptr);
void read_surfdat(int length, FO_surf* surf_ptr);
void read_decdat_mu(int FO_length, int N_stable, double** particle_mu);
int read_resonance(particle_info* particle);
void calculate_particle_mu(int IEOS, int Nparticle, FO_surf* FOsurf_ptr, int FO_length, particle_info* particle, double** particle_mu);

#endif

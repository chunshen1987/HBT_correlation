// ver 1.1
#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<stdlib.h>

#include "parameters.h"
#include "readindata.h"
#include "arsenal.h"
#include "ParameterReader.h"
#include "Table.h"

using namespace std;

read_FOdata::read_FOdata(ParameterReader* paraRdr_in, string path_in)
{
   paraRdr = paraRdr_in;
   path = path_in;
   mode = paraRdr->getVal("hydro_mode");
   turn_on_bulk = paraRdr->getVal("turn_on_bulk");
   n_eta_skip = 0;
}

read_FOdata::~read_FOdata()
{

}

int read_FOdata::get_number_of_freezeout_cells()
{
   int number_of_cells = 0;
   if (mode == 0)  // outputs from VISH2+1
   {
      ostringstream decdatfile;
      decdatfile << path << "/decdat2.dat";
      Table block_file(decdatfile.str().c_str());
      number_of_cells = block_file.getNumberOfRows();
   }
   else if (mode == 1)  // outputs from MUSIC boost-invariant
   {
      ostringstream surface_file;
      surface_file << path << "/surface.dat";
      Table block_file(surface_file.str().c_str());

      // determine number of the eta slides that are output
      double eta_target = block_file.get(4, 1);
      int i_eta = 1;
      int num_temp = 0;
      for(int i = 0; i < block_file.getNumberOfRows(); i++)
      {
         if(block_file.get(4, i+1) == eta_target)
            num_temp++;
      }
      number_of_cells = num_temp;
      n_eta_skip = block_file.getNumberOfRows()/number_of_cells;
   }
   else if (mode == 2)  // outputs from MUSIC full (3+1)-d
   {
      ostringstream surface_file;
      surface_file << path << "/surface.dat";
      Table block_file(surface_file.str().c_str());
      number_of_cells = block_file.getNumberOfRows();
   }
 
   return(number_of_cells);
}

void read_FOdata::read_in_freeze_out_data(int length, FO_surf* surf_ptr)
{
   if(mode == 0)     // VISH2+1 outputs
      read_FOsurfdat_VISH2p1(length, surf_ptr);
   else if (mode == 1)   // MUSIC boost invariant outputs
      read_FOsurfdat_MUSIC_boost_invariant(length, surf_ptr);
   else if (mode == 2)   // MUSIC full (3+1)-d outputs
      read_FOsurfdat_MUSIC(length, surf_ptr);

   return;
}

int read_FOdata::read_in_chemical_potentials(string path, int FO_length, FO_surf* surf_ptr, particle_info* particle_ptr)
{
   int Nparticle = 0;
   int N_stableparticle;
   Table mu_table;
   if(mode == 0)      // VISH2+1 output
   {
       ifstream particletable("EOS/EOS_particletable.dat");
       particletable >> N_stableparticle;
       particletable.close();
   }
   else if (mode == 1 || mode == 2)   // music output
   {
       // determine the type of the EOS 
       ostringstream config_file;
       config_file << path << "/input";
       ifstream configuration(config_file.str().c_str(), ios::in);
       string temp1;
       string temp_name;
       while(!configuration.eof())
       {
           getline(configuration, temp1);
           stringstream ss(temp1);
           ss >> temp_name;
           if(temp_name == "EOS_to_use")
           {
              ss >> IEOS_music;
              break;
           }
       }
       ifstream particletable;
       if(IEOS_music == 2)
           N_stableparticle = 0;
       else if(IEOS_music == 3)         // s95p-v1-PCE150
           particletable.open("EOS/EOS_tables/s95p-v1-PCE150/EOS_particletable.dat");
       else if(IEOS_music == 4)   // s95p-v1-PCE155
           particletable.open("EOS/EOS_tables/s95p-v1-PCE155/EOS_particletable.dat");
       else if(IEOS_music == 5)   // s95p-v1-PCE160
           particletable.open("EOS/EOS_tables/s95p-v1-PCE160/EOS_particletable.dat");
       else if(IEOS_music == 6)   // s95p-v0-PCE165
           particletable.open("EOS/EOS_tables/s95p-v0-PCE165/EOS_particletable.dat");
       else
       {
           cout << "invalid IEOS_music: " << IEOS_music << endl;
           exit(-1);
       }
       particletable >> N_stableparticle;
       particletable.close();
   }
   
   //read particle resonance decay table
   for (int i = 0; i < Maxparticle; i++) 
       particle_ptr[i].decays = 0; // to avoid infinite loop
   Nparticle = read_resonances_list(particle_ptr);
   cout << "total number of particle species: " << Nparticle << endl;

   if(N_stableparticle > 0)
   {
       cout << " -- EOS is partially chemical equilibrium " << endl;
       double** particle_mu = new double* [N_stableparticle];
       for(int i = 0; i < N_stableparticle; i++)
           particle_mu[i] = new double [FO_length];

       if(mode == 0)
           read_decdat_mu(FO_length, N_stableparticle, particle_mu);
       else if (mode == 1 || mode == 2)
           read_chemical_potentials_music(FO_length, surf_ptr, N_stableparticle, particle_mu);

       calculate_particle_mu(Nparticle, surf_ptr, FO_length, particle_ptr, particle_mu);
       for(int i = 0; i < N_stableparticle; i++)
           delete [] particle_mu[i];
       delete [] particle_mu;
   }
   else
   {
      cout << " -- EOS is chemical equilibrium. " << endl;
      for(int i = 0; i < Nparticle; i++)
        for(int j = 0; j < FO_length; j++)
           surf_ptr[j].particle_mu[i] = 0.0e0;
   }

   return(Nparticle);
}

void read_FOdata::read_decdat(int length, FO_surf* surf_ptr)
{
  double temp, temp_vx, temp_vy;
  cout<<" -- Read in information on freeze out surface...";
  ostringstream decdat_stream;
  decdat_stream << path << "/decdat2.dat";
  ifstream decdat(decdat_stream.str().c_str(), ios::in);
  for(int i=0; i<length; i++)
  {
     decdat >> surf_ptr[i].tau;
     surf_ptr[i].eta = 0.0;

     decdat >> surf_ptr[i].da0;
     decdat >> surf_ptr[i].da1;
     decdat >> surf_ptr[i].da2;
     surf_ptr[i].da3 = 0.0;

     decdat >> temp_vx;
     decdat >> temp_vy;
     surf_ptr[i].u0 = 1./sqrt(1. - temp_vx*temp_vx - temp_vy*temp_vy);
     surf_ptr[i].u1 = surf_ptr[i].u0*temp_vx;
     surf_ptr[i].u2 = surf_ptr[i].u0*temp_vy;
     surf_ptr[i].u3 = 0.0;

     decdat >> surf_ptr[i].Edec;
     decdat >> surf_ptr[i].Bn;
     decdat >> surf_ptr[i].Tdec;
     decdat >> surf_ptr[i].muB;
     decdat >> surf_ptr[i].muS;
     decdat >> surf_ptr[i].Pdec;

     decdat >> surf_ptr[i].pi33;
     decdat >> surf_ptr[i].pi00;
     decdat >> surf_ptr[i].pi01;
     decdat >> surf_ptr[i].pi02;
     surf_ptr[i].pi03 = 0.0;
     decdat >> surf_ptr[i].pi11;
     decdat >> surf_ptr[i].pi12;
     surf_ptr[i].pi13 = 0.0;
     decdat >> surf_ptr[i].pi22;
     surf_ptr[i].pi23 = 0.0;

     decdat >> temp;
     if(turn_on_bulk == 1)
         surf_ptr[i].bulkPi = temp;
     else
         surf_ptr[i].bulkPi = 0.0;
  }
  decdat.close();
  cout<<"done"<<endl;
  return;
}

void read_FOdata::read_surfdat(int length, FO_surf* surf_ptr)
{
  cout<<" -- Read spatial positions of freeze out surface...";
  ostringstream surfdat_stream;
  double dummy;
  char rest_dummy[512];
  surfdat_stream << path << "/surface.dat";
  ifstream surfdat(surfdat_stream.str().c_str(), ios::in);
  for(int i=0; i<length; i++)
  {
     surfdat >> dummy >> dummy;
     surfdat >> surf_ptr[i].xpt;
     surfdat >> surf_ptr[i].ypt;
     surfdat.getline(rest_dummy, 512);
  }
  surfdat.close();
  cout<<"done"<<endl;
  return;
}

void read_FOdata::read_FOsurfdat_VISH2p1(int length, FO_surf* surf_ptr)
{
  cout << " -- Loading the decoupling data from VISH2+1 ...." << endl;
  //read the data arrays for the decoupling information
  read_decdat(length, surf_ptr);
  //read the positions of the freeze out surface
  read_surfdat(length, surf_ptr);
  return;
}

void read_FOdata::read_FOsurfdat_MUSIC_boost_invariant(int length, FO_surf* surf_ptr)
{
  cout<<" -- Read spatial positions of freeze out surface from MUSIC (boost-invariant) ...";
  ostringstream surfdat_stream;
  double dummy;
  double temp;
  string input;
  double temp_tau, temp_xpt, temp_ypt, temp_eta;
  double eta_target;
  double deta = 100.0;
  int idx = 0;
  char rest_dummy[512];
  surfdat_stream << path << "/surface.dat";
  ifstream surfdat(surfdat_stream.str().c_str(), ios::in);
  for(int i = 0; i < length*n_eta_skip; i++)
  {
     getline(surfdat, input, '\n' );
     stringstream ss(input);
     ss >> temp_tau >> temp_xpt >> temp_ypt >> temp_eta;
     if(i == 0) 
         eta_target = temp_eta;
     double eta_diff = fabs(temp_eta - eta_target);
     if(eta_diff < 1e-15)
     {
         // freeze out position
         surf_ptr[idx].tau = temp_tau;
         surf_ptr[idx].xpt = temp_xpt;
         surf_ptr[idx].ypt = temp_ypt;
         surf_ptr[idx].eta = temp_eta;

         // freeze out normal vectors
         ss >> surf_ptr[idx].da0;
         ss >> surf_ptr[idx].da1;
         ss >> surf_ptr[idx].da2;
         ss >> surf_ptr[idx].da3;

         // flow velocity
         ss >> surf_ptr[idx].u0;
         ss >> surf_ptr[idx].u1;
         ss >> surf_ptr[idx].u2;
         ss >> surf_ptr[idx].u3;

         // thermodynamic quantities at freeze out
         ss >> dummy;
         surf_ptr[idx].Edec = dummy*hbarC;   
         ss >> dummy;
         surf_ptr[idx].Tdec = dummy*hbarC;
         ss >> dummy;
         surf_ptr[idx].muB = dummy*hbarC;
         ss >> dummy;              // (e+P)/T
         surf_ptr[idx].Pdec = dummy*surf_ptr[idx].Tdec - surf_ptr[idx].Edec;
         surf_ptr[idx].Bn = 0.0;
         surf_ptr[idx].muS = 0.0;

         // dissipative quantities at freeze out
         ss >> dummy;                       // 1/fm^4
         surf_ptr[idx].pi00 = dummy*hbarC;  // GeV/fm^3
         ss >> dummy;
         surf_ptr[idx].pi01 = dummy*hbarC;
         ss >> dummy;
         surf_ptr[idx].pi02 = dummy*hbarC;
         ss >> dummy;
         surf_ptr[idx].pi03 = dummy*hbarC;
         ss >> dummy;
         surf_ptr[idx].pi11 = dummy*hbarC;
         ss >> dummy;
         surf_ptr[idx].pi12 = dummy*hbarC;
         ss >> dummy;
         surf_ptr[idx].pi13 = dummy*hbarC;
         ss >> dummy;
         surf_ptr[idx].pi22 = dummy*hbarC;
         ss >> dummy;
         surf_ptr[idx].pi23 = dummy*hbarC;
         ss >> dummy;
         surf_ptr[idx].pi33 = dummy*hbarC;
         if(turn_on_bulk == 1)
         {
             ss >> dummy;
             surf_ptr[idx].bulkPi = dummy*hbarC;
         }
         else
             surf_ptr[idx].bulkPi = 0.0;
         idx++;
     }
     else
     {
         if(eta_diff < deta) deta = eta_diff;
     }
  }
  surfdat.close();

  // take out deta weight
  if(n_eta_skip > 1)
  {
     for(int i = 0; i < length; i++)
     {
         surf_ptr[i].da0 = surf_ptr[i].da0/deta;
         surf_ptr[i].da1 = surf_ptr[i].da1/deta;
         surf_ptr[i].da2 = surf_ptr[i].da2/deta;
         surf_ptr[i].da3 = surf_ptr[i].da3/deta;
     }
  }
  else
  {
       // determine deta
       ostringstream config_file;
       config_file << path << "/input";
       double eta_size;
       int n_eta;
       ifstream configuration(config_file.str().c_str(), ios::in);
       string temp1;
       string temp_name;
       while(!configuration.eof())
       {
           getline(configuration, temp1);
           stringstream ss(temp1);
           ss >> temp_name;
           if(temp_name == "Eta_grid_size")
           {
              ss >> eta_size;
              break;
           }
       }
       configuration.close();
       configuration.open(config_file.str().c_str(), ios::in);
       while(!configuration.eof())
       {
           getline(configuration, temp1);
           stringstream ss(temp1);
           ss >> temp_name;
           if(temp_name == "Grid_size_in_eta")
           {
              ss >> n_eta;
              break;
           }
       }
       configuration.close();
       deta = eta_size/n_eta;
       for(int i = 0; i < length; i++)
       {
          surf_ptr[i].da0 = surf_ptr[i].da0/deta;
          surf_ptr[i].da1 = surf_ptr[i].da1/deta;
          surf_ptr[i].da2 = surf_ptr[i].da2/deta;
          surf_ptr[i].da3 = surf_ptr[i].da3/deta;
       }
  }
  cout << "done" << endl;
  return;
}

void read_FOdata::read_FOsurfdat_MUSIC(int length, FO_surf* surf_ptr)
{
  cout<<" -- Read spatial positions of freeze out surface from MUSIC...";
  ostringstream surfdat_stream;
  double dummy;
  char rest_dummy[512];
  surfdat_stream << path << "/surface.dat";
  ifstream surfdat(surfdat_stream.str().c_str(), ios::in);
  for(int i=0; i<length; i++)
  {
     // freeze out position
     surfdat >> surf_ptr[i].tau;
     surfdat >> surf_ptr[i].xpt;
     surfdat >> surf_ptr[i].ypt;
     surfdat >> surf_ptr[i].eta;

     // freeze out normal vectors
     surfdat >> surf_ptr[i].da0;
     surfdat >> surf_ptr[i].da1;
     surfdat >> surf_ptr[i].da2;
     surfdat >> surf_ptr[i].da3;

     // flow velocity
     surfdat >> surf_ptr[i].u0;
     surfdat >> surf_ptr[i].u1;
     surfdat >> surf_ptr[i].u2;
     surfdat >> surf_ptr[i].u3;

     // thermodynamic quantities at freeze out
     surfdat >> dummy;
     surf_ptr[i].Edec = dummy*hbarC;
     surfdat >> dummy;
     surf_ptr[i].Tdec = dummy*hbarC;
     surfdat >> dummy;
     surf_ptr[i].muB = dummy*hbarC;
     surfdat >> dummy;                    //(e+p)/T
     surf_ptr[i].Pdec = dummy*surf_ptr[i].Tdec - surf_ptr[i].Edec;
     surf_ptr[i].Bn = 0.0;
     surf_ptr[i].muS = 0.0;

     // dissipative quantities at freeze out
     surfdat >> dummy;
     surf_ptr[i].pi00 = dummy*hbarC;
     surfdat >> dummy;
     surf_ptr[i].pi01 = dummy*hbarC;
     surfdat >> dummy;
     surf_ptr[i].pi02 = dummy*hbarC;
     surfdat >> dummy;
     surf_ptr[i].pi03 = dummy*hbarC;
     surfdat >> dummy;
     surf_ptr[i].pi11 = dummy*hbarC;
     surfdat >> dummy;
     surf_ptr[i].pi12 = dummy*hbarC;
     surfdat >> dummy;
     surf_ptr[i].pi13 = dummy*hbarC;
     surfdat >> dummy;
     surf_ptr[i].pi22 = dummy*hbarC;
     surfdat >> dummy;
     surf_ptr[i].pi23 = dummy*hbarC;
     surfdat >> dummy;
     surf_ptr[i].pi33 = dummy*hbarC;
     if(turn_on_bulk == 1)
     {
         surfdat >> dummy;
         surf_ptr[i].bulkPi = dummy*hbarC;
     }
     else
         surf_ptr[i].bulkPi = 0.0;
  }
  surfdat.close();
  cout << "done" << endl;
  return;
}

void read_FOdata::read_decdat_mu(int FO_length, int N_stable, double** particle_mu)
{
  cout<<" -- Read chemical potential for stable particles...";
  ostringstream decdat_mu_stream;
  double dummy;
  decdat_mu_stream << path << "/decdat_mu.dat";
  ifstream decdat_mu(decdat_mu_stream.str().c_str(), ios::in);

  //For backward compatibility: decdat_mu.dat can be one line or FO_length lines
  for(int j=0; j<FO_length; j++)
  {
    decdat_mu >> dummy;  //not used in the code plz ignore it

    if(decdat_mu.eof())
    {
      for(int k=j; k<FO_length; k++)
        for(int i=0; i<N_stable; i++)
           particle_mu[i][k]=particle_mu[i][j-1];
      break;
    }

    for(int i=0; i<N_stable; i++)
    {
       decdat_mu >> particle_mu[i][j];
    }
  }

  cout<<"done" << endl;
  return;
}

void read_FOdata::read_chemical_potentials_music(int FO_length, FO_surf* FOsurf_ptr, int N_stable, double** particle_mu)
{
  cout << " -- Interpolating chemical potentials for stable particles (MUSIC IEOS = " << IEOS_music << ") ...";

  Table mu_table;
  if(IEOS_music == 3)
      mu_table.loadTableFromFile("EOS/EOS_tables/s95p-v1-PCE150/EOS_Mu.dat");
  else if (IEOS_music == 4)
      mu_table.loadTableFromFile("EOS/EOS_tables/s95p-v1-PCE155/EOS_Mu.dat");
  else if (IEOS_music == 5)
      mu_table.loadTableFromFile("EOS/EOS_tables/s95p-v1-PCE160/EOS_Mu.dat");
  else if (IEOS_music == 6)
      mu_table.loadTableFromFile("EOS/EOS_tables/s95p-v1-PCE165/EOS_Mu.dat");

  double edec_pre = 0.0e0;
  for(int j = 0; j < FO_length; j++)
  {
    double edec = FOsurf_ptr[j].Edec;
    if(fabs(edec - edec_pre) > 1e-15)
    {
       edec_pre = edec;
       for(int i = 0; i < N_stable; i++)
          particle_mu[i][j] = mu_table.interp(1, i+2, edec);
    }
    else
    {
       for(int i = 0; i < N_stable; i++)
          particle_mu[i][j] = particle_mu[i][j-1];
    }
  }

  cout<<"done" << endl;
  return;
}

int read_FOdata::read_resonances_list(particle_info* particle)
{
   double eps = 1e-15;
   int Nparticle=0;
   cout << " -- Read in particle resonance decay table...";
   ifstream resofile("EOS/pdg.dat", ios::in);
   int local_i = 0;
   int dummy_int;
   while (!resofile.eof())
   {
      resofile >> particle[local_i].monval;
      resofile >> particle[local_i].name;
      resofile >> particle[local_i].mass;
      resofile >> particle[local_i].width;
      resofile >> particle[local_i].gspin;	      //spin degeneracy
      resofile >> particle[local_i].baryon;
      resofile >> particle[local_i].strange;
      resofile >> particle[local_i].charm;
      resofile >> particle[local_i].bottom;
      resofile >> particle[local_i].gisospin;     //isospin degeneracy
      resofile >> particle[local_i].charge;
      resofile >> particle[local_i].decays;
      for (int j = 0; j < particle[local_i].decays; j++)
      {
         resofile >> dummy_int;
         resofile >> particle[local_i].decays_Npart[j];
         resofile >> particle[local_i].decays_branchratio[j];
         resofile >> particle[local_i].decays_part[j][0];
         resofile >> particle[local_i].decays_part[j][1];
         resofile >> particle[local_i].decays_part[j][2];
         resofile >> particle[local_i].decays_part[j][3];
         resofile >> particle[local_i].decays_part[j][4];
      }

      //decide whether particle is stable under strong interactions
      if(particle[local_i].decays_Npart[0] == 1)
         particle[local_i].stable = 1;
      else
         particle[local_i].stable = 0;

      //add anti-particle entry
      if(particle[local_i].baryon == 1)
      {
         local_i++;
         particle[local_i].monval = -particle[local_i-1].monval;
         ostringstream antiname;
         antiname << "Anti-" << particle[local_i-1].name;
         particle[local_i].name = antiname.str();
         particle[local_i].mass = particle[local_i-1].mass;
         particle[local_i].width = particle[local_i-1].width;
         particle[local_i].gspin = particle[local_i-1].gspin;
         particle[local_i].baryon = -particle[local_i-1].baryon;
         particle[local_i].strange = -particle[local_i-1].strange;
         particle[local_i].charm = -particle[local_i-1].charm;
         particle[local_i].bottom = -particle[local_i-1].bottom;
         particle[local_i].gisospin = particle[local_i-1].gisospin;
         particle[local_i].charge = -particle[local_i-1].charge;
         particle[local_i].decays = particle[local_i-1].decays;
         particle[local_i].stable = particle[local_i-1].stable;
         for (int j = 0; j < particle[local_i].decays; j++)
         {
            particle[local_i].decays_Npart[j]=particle[local_i-1].decays_Npart[j];
            particle[local_i].decays_branchratio[j]=particle[local_i-1].decays_branchratio[j];
            for (int k=0; k< Maxdecaypart; k++)
            {
               if(particle[local_i-1].decays_part[j][k] == 0)
                  particle[local_i].decays_part[j][k]= particle[local_i-1].decays_part[j][k];
               else
               {
                  int idx; 
                  for(idx = 0; idx < local_i; idx++) // find the index for decay particle
                     if(particle[idx].monval == particle[local_i-1].decays_part[j][k])
                        break;
                  if(idx == local_i && particle[local_i-1].stable == 0 && particle[local_i-1].decays_branchratio[j] > eps)  // check
                  {
                     cout << "Error: can not find decay particle index for anti-baryon!" << endl;
                     cout << "particle monval : " << particle[local_i-1].decays_part[j][k] << endl;
                     exit(1);
                  }
                  if(particle[idx].baryon == 0 && particle[idx].charge == 0 && particle[idx].strange == 0)
                     particle[local_i].decays_part[j][k]= particle[local_i-1].decays_part[j][k];
                  else
                     particle[local_i].decays_part[j][k]= -particle[local_i-1].decays_part[j][k];
               }
            }
         }
       }
       local_i++;	// Add one to the counting variable "i" for the meson/baryon
   }
   resofile.close();
   Nparticle=local_i-1; //take account the final fake one
   for(int i=0; i < Nparticle; i++)
   {
      if(particle[i].baryon==0)
         particle[i].sign=-1;
      else
         particle[i].sign=1;
   }
   return(Nparticle);
}

void read_FOdata::calculate_particle_mu(int Nparticle, FO_surf* FOsurf_ptr, int FO_length, particle_info* particle, double** particle_mu)
{
   int Nstable_particle;
   int Idummy;
   char cdummy[256];
   cout << " -- Read particle table and calculating chemical potential for particles..." << endl;
   ifstream particletable;
   if(mode == 0)
   {
      particletable.open("EOS/EOS_particletable.dat", ios::in);
   }
   else if (mode == 1 || mode == 2)
   {
      if(IEOS_music == 3)
         particletable.open("EOS/EOS_tables/s95p-v1-PCE150/EOS_particletable.dat", ios::in);
      else if(IEOS_music == 4)
         particletable.open("EOS/EOS_tables/s95p-v1-PCE155/EOS_particletable.dat", ios::in);
      else if(IEOS_music == 5)
         particletable.open("EOS/EOS_tables/s95p-v1-PCE160/EOS_particletable.dat", ios::in);
      else if(IEOS_music == 6)
         particletable.open("EOS/EOS_tables/s95p-v1-PCE165/EOS_particletable.dat", ios::in);
      else
      {
         cout << "invalid EOS option for MUSIC: " << IEOS_music << endl;
         exit(-1);
      }
   }
   else
   {
      cout << "invalid hydro mode: " << mode << endl;
      exit(-1);
   }

   particletable >> Nstable_particle;
   double *stable_particle_monval = new double [Nstable_particle];
   for(int i = 0; i < Nstable_particle; i++)
   {
       particletable >> Idummy >> stable_particle_monval[i];
       particletable.getline(cdummy, 256);
   }
   particletable.close();

   // assign chemical potentials for stable particles
   for(int i = 0; i < Nstable_particle; i++)
      for(int j=0; j<Nparticle; j++)
         if(particle[j].monval == stable_particle_monval[i])
         {
            particle[j].stable = 1;
            for(int k=0; k<FO_length; k++)
                FOsurf_ptr[k].particle_mu[j] = particle_mu[i][k];
            break;
         }

   // calculating chemical potentials for unstable resonances
   print_progressbar(-1);
   for(int i=0; i < Nparticle ; i++)
   {
      if(particle[i].stable==0)
      {
         for(int j=0; j < particle[i].decays; j++)
         {
            for(int k=0; k < abs(particle[i].decays_Npart[j]); k++)
            {
               for(int l=0; l < Nparticle; l++)
               {
                  if(particle[i].decays_part[j][k] == particle[l].monval)
                  {
                     for(int m=0; m<FO_length; m++)
                       FOsurf_ptr[m].particle_mu[i] += particle[i].decays_branchratio[j]*FOsurf_ptr[m].particle_mu[l];
                     break;
                  }
                  if(l==Nparticle-1)
                     cout<<"warning: can not find particle" <<  particle[i].name << endl;
               }
            }
         }
      }
      print_progressbar((double)(i)/Nparticle);
   }
   print_progressbar(1);

   return;
}


//===============================================================================
//  calculate the HBT radii from VISH2+1
//
//
//  Programmer: Chun Shen
//       Email: shen.201@asc.ohio-state.edu
//
//        Date: 11/22/11
//  The HBT part of the code is reorganized in the class style for easy 
//  maintainance and furture extension.  --- 04/15/2012
//
//===============================================================================


#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include<math.h>
#include<sys/time.h>


#include<gsl/gsl_sf_bessel.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#include "Stopwatch.h"
#include "parameters.h"
#include "readindata.h"
#include "HBT.h"

using namespace std;

int main()
{
   Stopwatch sw;
   Stopwatch sw_total;
   sw_total.tic();
   sw.tic();
   hydropara* hydropara_ptr=new hydropara;
   //load hydro parameters
   read_hydropar(hydropara_ptr);
   
   //load freeze out information
   int FO_length = 0;
   ostringstream decdatfile;
   cout << "Loading the decoupling data...." << endl;
   decdatfile << path << "/decdat2.dat";
   FO_length=get_filelength(decdatfile.str().c_str());
   cout <<"Total number of freeze out fluid cell: " <<  FO_length << endl;

   //read the data arrays for the decoupling information
   FO_surf* FOsurf_ptr = new FO_surf[FO_length];
   read_decdat(FO_length, FOsurf_ptr);
   
   //read the positions of the freeze out surface
   read_surfdat(FO_length, FOsurf_ptr);
   
   //read the chemical potential on the freeze out surface
   int N_stableparticle;
   ifstream particletable("EOS/EOS_particletable.dat");
   particletable >> N_stableparticle;
   double** particle_mu = new double* [N_stableparticle];
   for(int i=0; i<N_stableparticle; i++)
      particle_mu[i] = new double [FO_length];
   for(int i=0; i<N_stableparticle; i++)
      for(int j=0; j<FO_length; j++)
         particle_mu[i][j] = 0.0;
   if(N_stableparticle >0)
   {
      if(hydropara_ptr->IEOS==7)       //for s95p_PCE
         read_decdat_mu(FO_length, N_stableparticle, particle_mu);
   }

   //read particle resonance decay table
   particle_info *particle = new particle_info [Maxparticle];
   int Nparticle=read_resonance(particle);
   cout <<"read in total " << Nparticle << " particles!" << endl;
   if(N_stableparticle >0)
   {
      cout << " EOS is partically chemical equilibrium " << endl;
      calculate_particle_mu(hydropara_ptr->IEOS, Nparticle, FOsurf_ptr, FO_length, particle, particle_mu);
   }
   else
   {
      cout << " EOS is chemical equilibrium. " << endl;
      for(int i=0; i<FO_length; i++)
      for(int j=0; j<Nparticle; j++)
         FOsurf_ptr[i].particle_mu[j] = 0.0e0;
   }
   //for(int i=0;i<Nparticle;i++)
     //cout << particle[i].mu << "  " << particle[i].sing << endl;
   sw.toc();
   cout << "read in data finished!" << endl;
   cout << "Used " << sw.takeTime() << " sec." << endl;

   //HBT calculations begins ...
   int particle_idx=1;  //for pion+
   cout << "Calculating "<< particle[particle_idx].name << endl;
   
   ostringstream OPHBTradii;
   OPHBTradii << path << "/HBT_radii.dat" ;
   ofstream OPfile(OPHBTradii.str().c_str());
   
   double localy = 0.0e0;
   double dp_T = (localp_T_max - localp_T_min)/(n_localp_T-1+1e-100);
   double dp_phi = (localp_phi_max - localp_phi_min)/(n_localp_phi-1+1e-100);
   for(int i = 0; i<n_localp_T; i++)
   {
      double localp_T = localp_T_min + i*dp_T;
      for(int j = 0; j<n_localp_phi; j++)
      {
         double localp_phi = localp_phi_min + j*dp_phi;
         cout << "K_T = " << localp_T << "  K_phi = " << localp_phi << endl;
         sw.tic();
         if(fabs(localy) > 1e-16)
         {
            cout << "not support y not equals 0 yet! Bye bye!" << endl;
            return 0;
         }
         
         HBT HBT_hadron(&particle[particle_idx], localp_T, localp_phi, localy, FO_length, particle_idx);
         HBT_hadron.SetEmissionData(FOsurf_ptr);
         HBT_hadron.Cal_HBTRadii_fromEmissionfunction();
         HBT_hadron.Cal_correlationfunction_1D_MC();
         HBT_hadron.Fit_Correlationfunction1D();
         HBT_hadron.Output_Correlationfunction_1D();
         exit(0);
         HBT_hadron.Cal_correlationfunction_3D_MC();
         HBT_hadron.Fit_Correlationfunction3D();

         OPfile << scientific << setprecision(7) << setw(15)
                << localp_T << "  " << localp_phi << "  "  
                << HBT_hadron.get_lambda_Correl() << "  " << HBT_hadron.get_lambda_Correl_err() << "  "
                << HBT_hadron.get_Rout_Correl() << "  " << HBT_hadron.get_Rout_Correl_err() << "  "
                << HBT_hadron.get_Rside_Correl() << "  " << HBT_hadron.get_Rside_Correl_err() << "  "
                << HBT_hadron.get_Rlong_Correl() << "  " << HBT_hadron.get_Rlong_Correl_err() << "  "
                << HBT_hadron.get_Ros_Correl() << "  " << HBT_hadron.get_Ros_Correl_err()
                << endl;
         sw.toc();
         cout << "Finished in " << sw.takeTime() << " sec." << endl;
      }
   }
   OPfile.close();
   sw_total.toc();
   cout << "Program totally finished in " << sw_total.takeTime() << " sec." << endl;
   return 0;
}

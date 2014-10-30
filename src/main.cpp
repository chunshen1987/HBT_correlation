//===============================================================================
//  calculate the HBT radii from VISH2+1
//
//
//  Programmer: Chun Shen
//       Email: shen.201@asc.ohio-state.edu
//
//        Date: 11/22/11
//  The HBT part of the code is reorganized in the class style for easy 
//  maintenance and future extension.  --- 04/15/2012
//
//=============================================================================

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
#include "arsenal.h"
#include "ParameterReader.h"

using namespace std;

int main(int argc, char *argv[])
{
   cout << endl
        << "                  iHBT                   " << endl
        << endl
        << "  Ver 1.2   ----- Chun Shen, 10/2014   " << endl;
   cout << endl << "**********************************************************" << endl;
   display_logo(2); // Hail to the king~
   cout << endl << "**********************************************************" << endl << endl;
   
   // Read-in parameters
   ParameterReader *paraRdr = new ParameterReader;
   paraRdr->readFromFile("parameters.dat");
   paraRdr->readFromArguments(argc, argv);
   paraRdr->echo();
   
   string path="results";

   Stopwatch sw;
   Stopwatch sw_total;
   sw_total.tic();
   sw.tic();

   //load freeze out information
   read_FOdata freeze_out_data(paraRdr, path);

   int FO_length = 0;
   FO_length = freeze_out_data.get_number_of_freezeout_cells();
   cout <<"total number of cells: " <<  FO_length << endl;

   FO_surf* FOsurf_ptr = new FO_surf[FO_length];
   for(int i=0; i<FO_length; i++)
     for(int j=0; j<Maxparticle; j++)
         FOsurf_ptr[i].particle_mu[j] = 0.0e0;

   freeze_out_data.read_in_freeze_out_data(FO_length, FOsurf_ptr);

   //read the chemical potential on the freeze out surface
   particle_info *particle = new particle_info [Maxparticle];
   int Nparticle = freeze_out_data.read_in_chemical_potentials(path, FO_length, FOsurf_ptr, particle);
   
   cout << endl << " -- Read in data finished!" << endl << endl;
   sw.toc();
   cout << "Used " << sw.takeTime() << " sec." << endl;

   //HBT calculations begins ...
   int particle_idx=1;  //for pion+
   HBT HBT_hadron(path, paraRdr, particle, particle_idx, FOsurf_ptr, FO_length);
   HBT_hadron.calculate_azimuthal_averaged_HBT_radii(0.0);
   /*
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
   */
   sw_total.toc();
   cout << "Program totally finished in " << sw_total.takeTime() << " sec." << endl;
   return 0;
}

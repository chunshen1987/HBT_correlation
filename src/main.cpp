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
   HBT_hadron.calculation_HBT_correlation(0.0);

   sw_total.toc();
   cout << "Program totally finished in " << sw_total.takeTime() << " sec." << endl;
   return 0;
}

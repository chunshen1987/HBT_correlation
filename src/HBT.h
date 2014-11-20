#ifndef HBT_H
#define HBT_H

#include<fstream>

#include "readindata.h"
#include "parameters.h"
#include "arsenal.h"
#include "ParameterReader.h"

using namespace std;

struct Emissionfunction_data
{
   double x;
   double y;
   double z;
   double t;
   double **data;
};

class HBT
{
   private:
      string path;
      ParameterReader* paraRdr;
      FO_surf* FOsurf_ptr;
      int FO_length;        // number of the freeze out surface cells
      particle_info* particle_ptr; // particle information
      int particle_id;             // particle id
     
      // array for eta_s
      int eta_s_npts;
      double *eta_s, *eta_s_weight;

      int INCLUDE_SHEAR_DELTAF, INCLUDE_BULK_DELTAF;
      int bulk_deltaf_type;

      int azimuthal_flag;
      int kT_differenitial_flag;

      double Psi_ev;
      int n_KT, n_Kphi;
      double *KT_array, *KT_weight;
      double *Kphi, *Kphi_weight;

      Emissionfunction_data* emission_S_K;

      // Emission function
      int Emissionfunction_length;  // length of the emission function array

      int flag_neg;

      int qnpts;
      double *q_out, *q_side, *q_long;
 
      int flag_1D_projection;
      int ndir;      // number of directions for 1D correlation function
      //store correlation functions
      double **Correl_1D_num, **Correl_1D_denorm;
      double ***Correl_3D_num, ***Correl_3D_denorm;

      double ***Correl_1D_phidiff_num, ***Correl_1D_phidiff_denorm;
      double ****Correl_3D_phidiff_num, ****Correl_3D_phidiff_denorm;

      int flag_MC;
      int MC_samples;
      double q_max;
      double *q_out_MC, *q_side_MC, *q_long_MC;
      double *Correl_MC_num, *Correl_MC_denorm;
      double **Correl_MC_phidiff_num, **Correl_MC_phidiff_denorm;
      
      //HBT radii calculated from emission functions
      int flag_HBT_from_source_variance;
      double R_out_EM;
      double R_side_EM;
      double R_long_EM;

   public:
      HBT(string path_in, ParameterReader* paraRdr, particle_info* particle_in, int particle_idx, FO_surf* FOsurf_ptr_in, int FOarray_length, double event_plane);
      ~HBT();

      void SetEmissionData(FO_surf* FO_surface, double K_rap);

      double Emissionfunction(double p0, double px, double py, double pz, FO_surf* surf);

      void Cal_HBTRadii_fromEmissionfunction(double K_y);
      void calculation_HBT_correlation(double y);
      void calculate_azimuthal_dependent_HBT_radii(double y);
      void calculate_azimuthal_averaged_HBT_radii(double y);
      void calculate_azimuthal_averaged_KT_integrated_HBT_radii(double y);
      void calculate_azimuthal_dependent_KT_integrated_HBT_radii(double y);

      void Cal_azimuthal_averaged_correlationfunction_1D(int iKT, double K_y);
      void Cal_azimuthal_averaged_correlationfunction_3D(int iKT, double K_y);
      void Cal_azimuthal_averaged_correlationfunction_MC(int iKT, double K_y);
      void Cal_azimuthal_averaged_KT_inte_correlationfunction_1D(double K_y);
      void Cal_azimuthal_averaged_KT_inte_correlationfunction_3D(double K_y);
      void Cal_azimuthal_averaged_KT_inte_correlationfunction_MC(double K_y);
      void Cal_azimuthal_dependent_correlationfunction_1D(int iKT, double K_y);
      void Cal_azimuthal_dependent_correlationfunction_3D(int iKT, double K_y);
      void Cal_azimuthal_dependent_correlationfunction_MC(int iKT, double K_y);
      void Cal_azimuthal_dependent_KT_inte_correlationfunction_1D(double K_y);
      void Cal_azimuthal_dependent_KT_inte_correlationfunction_3D(double K_y);
      void Cal_azimuthal_dependent_KT_inte_correlationfunction_MC(double K_y);

      void Output_Correlationfunction_1D(int iKT);
      void Output_Correlationfunction_3D(int iKT);
      void Output_Correlationfunction_MC(int iKT);
      void Output_Correlationfunction_azimuthal_dependent_1D(int iKT);
      void Output_Correlationfunction_azimuthal_dependent_3D(int iKT);
      void Output_Correlationfunction_azimuthal_dependent_MC(int iKT);

};

#endif

#ifndef HBT_H
#define HBT_H

#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>

#include<gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>            // gsl random number generators
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting


#include "readindata.h"
#include "parameters.h"
#include "Arsenal.h"

using namespace std;

struct Correlationfunction3D_data
{
  size_t data_length;
  double *q_o;
  double *q_s;
  double *q_l;
  double *y;
  double *sigma;
};

struct Correlationfunction1D_data
{
  size_t data_length;
  double *q;
  double *y;
  double *sigma;
};
      
int Fittarget_correlfun1D_f (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr);
int Fittarget_correlfun1D_df (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr);
int Fittarget_correlfun1D_fdf (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr);
int Fittarget_correlfun3D_f (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr);
int Fittarget_correlfun3D_df (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr);
int Fittarget_correlfun3D_fdf (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr);
int Fittarget_correlfun3D_f_withlambda (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr);
int Fittarget_correlfun3D_df_withlambda (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr);
int Fittarget_correlfun3D_fdf_withlambda (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr);

class HBT
{
   private:
      FO_surf* FOsurf_ptr;
      particle_info* particle_ptr;

      //particle information 
      string particle_name;
      double particle_mass;
      int particle_id;     //particle id
      double particle_sign;   //+/- 1 for Fermi/Bose statistic for baryon/meson
      double particle_gspin;  //particle degeneracy 
     
      //pair momentum
      double K_T, K_phi, K_y;
     
      //Emission function
      double* Emissionfunction_Data;
      double* Emissionfunction_t;
      double* Emissionfunction_x;
      double* Emissionfunction_y;
      double* Emissionfunction_z;
      double* Emissionfunction_Data_CDF;
      int Emissionfunction_length;
      int FO_length;

      double* eta_s;
      double* eta_s_weight;

      double spectra;

      double* q_out;
      double* q_side;
      double* q_long;

      //store correlation functions
      double* Correl_1D_out;
      double* Correl_1D_out_err;
      double* Correl_1D_side;
      double* Correl_1D_side_err;
      double* Correl_1D_long;
      double* Correl_1D_long_err;
      double*** Correl_3D;
      double*** Correl_3D_err;

      
      //HBT radii calculated from emission functions
      double R_out_EM;
      double R_side_EM;
      double R_long_EM;

      //HBT radii calculated from fitting correlaction functions
      double lambda_Correl;
      double R_out_Correl;
      double R_side_Correl;
      double R_long_Correl;
      double R_os_Correl;
      double lambda_Correl_err;
      double R_out_Correl_err;
      double R_side_Correl_err;
      double R_long_Correl_err;
      double R_os_Correl_err;

   public:
      HBT(particle_info* particle_in, int particle_idx, FO_surf* FOsurf_ptr_in, int FOarray_length);
      ~HBT();

      double get_lambda_Correl() {return(lambda_Correl);};
      double get_lambda_Correl_err() {return(lambda_Correl_err);};
      double get_Rout_Correl() {return(R_out_Correl);};
      double get_Rout_Correl_err() {return(R_out_Correl_err);};
      double get_Rside_Correl() {return(R_side_Correl);};
      double get_Rside_Correl_err() {return(R_side_Correl_err);};
      double get_Rlong_Correl() {return(R_long_Correl);};
      double get_Rlong_Correl_err() {return(R_long_Correl_err);};
      double get_Ros_Correl() {return(R_os_Correl);};
      double get_Ros_Correl_err() {return(R_os_Correl_err);};

      void SetEmissionData(FO_surf* FOsurf_ptr);

      double Emissionfunction(double p0, double px, double py, double pz, FO_surf* surf);

      void Cal_HBTRadii_fromEmissionfunction();
      void calculate_azimuthal_dependent_HBT_radii(double p_T, double p_phi, double y);

      void Cal_correlationfunction_1D();
      void Cal_correlationfunction_3D();
      void Cal_correlationfunction_1D_MC();
      void Cal_correlationfunction_3D_MC();
      int binary_search(double* dataset, int data_length, double value);

      void Output_Correlationfunction_1D();
      void Output_Correlationfunction_3D();

      void Fit_Correlationfunction1D();
      void Fit_Correlationfunction3D();
      void Fit_Correlationfunction3D_withlambda();
      int print_fit_state_1D (size_t iteration, gsl_multifit_fdfsolver * solver_ptr);
      int print_fit_state_3D (size_t iteration, gsl_multifit_fdfsolver * solver_ptr);
      int print_fit_state_3D_withlambda (size_t iteration, gsl_multifit_fdfsolver * solver_ptr);
      inline double get_fit_results(int i, gsl_multifit_fdfsolver * solver_ptr);
      inline double get_fit_err (int i, gsl_matrix * covariance_ptr);

};

#endif

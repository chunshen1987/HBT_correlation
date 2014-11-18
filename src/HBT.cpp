#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>

#include "HBT.h"

using namespace std;

HBT::HBT(string path_in, ParameterReader* paraRdr_in, particle_info* particle_in, int particle_idx, FO_surf* FOsurf_ptr_in, int FOarray_length, double event_plane)
{
   path = path_in;
   paraRdr =  paraRdr_in;
   particle_ptr = particle_in;
   particle_id = particle_idx;

   FOsurf_ptr = FOsurf_ptr_in;
   FO_length = FOarray_length;

   Psi_ev = event_plane;

   // initialize eta_s array
   eta_s_npts = paraRdr->getVal("eta_s_npts");
   double eta_s_f = paraRdr->getVal("eta_s_f");
   eta_s = new double [eta_s_npts];
   eta_s_weight = new double [eta_s_npts];
   gauss_quadrature(eta_s_npts, 1, 0.0, 0.0, 0.0, eta_s_f, eta_s, eta_s_weight);

   azimuthal_flag = paraRdr->getVal("azimuthal_flag");
   kT_differenitial_flag = paraRdr->getVal("kT_differenitial_flag");
   flag_1D_projection = paraRdr->getVal("flag_1D_projection");
   flag_HBT_from_source_variance = paraRdr->getVal("flag_HBT_from_source_variance");
   ndir = paraRdr->getVal("ndir");

   // initialize KT and Kphi array
   n_KT = paraRdr->getVal("n_KT");
   KT_array = new double [n_KT];
   KT_weight = new double [n_KT];
   double KT_min = paraRdr->getVal("KT_min");
   double KT_max = paraRdr->getVal("KT_max");
   if(kT_differenitial_flag == 0)
       gauss_quadrature(n_KT, 1, 0.0, 0.0, KT_min, KT_max, KT_array, KT_weight);
   else
   {
       double dKT = (KT_max - KT_min)/(n_KT - 1);
       for(int i = 0; i < n_KT; i++)
           KT_array[i] = KT_min + i*dKT;
   }
   n_Kphi = paraRdr->getVal("n_Kphi");
   Kphi = new double [n_Kphi];
   Kphi_weight = new double [n_Kphi];
   gauss_quadrature(n_Kphi, 1, 0.0, 0.0, 0.0+Psi_ev, 2*M_PI+Psi_ev, Kphi, Kphi_weight);

   // initialize emission function
   Emissionfunction_length = FO_length*eta_s_npts;
   emission_S_K = new Emissionfunction_data [Emissionfunction_length];
   for(int i=0; i<Emissionfunction_length; i++)
   {
      emission_S_K[i].t = 0.0;
      emission_S_K[i].x = 0.0;
      emission_S_K[i].y = 0.0;
      emission_S_K[i].z = 0.0;
      emission_S_K[i].data = new double* [n_KT];
      for(int j = 0; j < n_KT; j++)
      {
         emission_S_K[i].data[j] = new double [n_Kphi];
         for(int k = 0; k < n_Kphi; k++)
            emission_S_K[i].data[j][k] = 0.0;
      }
   }

   INCLUDE_SHEAR_DELTAF = paraRdr->getVal("turn_on_shear");
   INCLUDE_BULK_DELTAF = paraRdr->getVal("turn_on_bulk");
   flag_neg = paraRdr->getVal("flag_neg");

   // initialize correlation function
   qnpts = paraRdr->getVal("qnpts");
   double init_q = paraRdr->getVal("init_q");
   double delta_q = paraRdr->getVal("delta_q");
   q_out = new double [qnpts];
   q_side = new double [qnpts];
   q_long = new double [qnpts];
   for(int i=0; i<qnpts; i++)
   {
      q_out[i] = init_q + (double)i * delta_q;
      q_side[i] = init_q + (double)i * delta_q;
      q_long[i] = init_q + (double)i * delta_q;
   }

   flag_MC = paraRdr->getVal("flag_MC");
   if(flag_MC == 1)
   {
       MC_samples = paraRdr->getVal("MC_samples");
       q_max = paraRdr->getVal("q_max");
       q_out_MC = new double [MC_samples];
       q_side_MC = new double [MC_samples];
       q_long_MC = new double [MC_samples];
       if(azimuthal_flag == 0)
       {
          Correl_MC_num = new double [MC_samples];
          Correl_MC_denorm = new double [MC_samples];
       }
       else
       {
          Correl_MC_phidiff_num = new double* [n_Kphi];
          Correl_MC_phidiff_denorm = new double* [n_Kphi];
          for(int i = 0; i < n_Kphi; i++)
          {
              Correl_MC_phidiff_num[i] = new double [MC_samples];
              Correl_MC_phidiff_denorm[i] = new double [MC_samples];
          }
       }
   }
   else
   {
      if(azimuthal_flag == 0)
      {
         if(flag_1D_projection == 1)
         {
            Correl_1D_num = new double* [ndir];
            Correl_1D_denorm = new double* [ndir];
            for(int i = 0; i < ndir; i++)
            {
               Correl_1D_num[i] = new double [qnpts];
               Correl_1D_denorm[i] = new double [qnpts];
               for(int j = 0; j < qnpts; j++)
               {
                  Correl_1D_num[i][j] = 0.0;
                  Correl_1D_denorm[i][j] = 0.0;
               }
            }
         }
         else
         {
            Correl_3D_num = new double** [qnpts];
            Correl_3D_denorm = new double** [qnpts];
            for(int i=0; i<qnpts; i++)
            {
               Correl_3D_num[i] = new double* [qnpts];
               Correl_3D_denorm[i] = new double* [qnpts];
               for(int j=0; j<qnpts; j++)
               {
                  Correl_3D_num[i][j] = new double [qnpts];
                  Correl_3D_denorm[i][j] = new double [qnpts];
                  for(int k=0; k<qnpts; k++)
                  {
                     Correl_3D_num[i][j][k] = 0.0;
                     Correl_3D_denorm[i][j][k] = 0.0;
                  }
               }
            }
         }
      }
      else
      {
         if(flag_1D_projection == 1)
         {
            Correl_1D_phidiff_num = new double** [n_Kphi];
            Correl_1D_phidiff_denorm = new double** [n_Kphi];
            for(int i = 0; i < n_Kphi; i++)
            {
                Correl_1D_phidiff_num[i] = new double* [ndir];
                Correl_1D_phidiff_denorm[i] = new double* [ndir];
                for(int j = 0; j < ndir; j++)
                {
                    Correl_1D_phidiff_num[i][j] = new double [qnpts];
                    Correl_1D_phidiff_denorm[i][j] = new double [qnpts];
                    for(int k = 0; k < qnpts; k++)
                    {
                        Correl_1D_phidiff_num[i][j][k] = 0.0;
                        Correl_1D_phidiff_denorm[i][j][k] = 0.0;
                    }
                }
            }
         }
         else
         {
            Correl_3D_phidiff_num = new double*** [n_Kphi];
            Correl_3D_phidiff_denorm = new double*** [n_Kphi];
            for(int l = 0; l < n_Kphi; l++)
            {
               Correl_3D_phidiff_num[l] = new double** [qnpts];
               Correl_3D_phidiff_denorm[l] = new double** [qnpts];
               for(int i=0; i<qnpts; i++)
               {
                  Correl_3D_phidiff_num[l][i] = new double* [qnpts];
                  Correl_3D_phidiff_denorm[l][i] = new double* [qnpts];
                  for(int j=0; j<qnpts; j++)
                  {
                     Correl_3D_phidiff_num[l][i][j] = new double [qnpts];
                     Correl_3D_phidiff_denorm[l][i][j] = new double [qnpts];
                     for(int k=0; k<qnpts; k++)
                     {
                        Correl_3D_phidiff_num[l][i][j][k] = 0.0;
                        Correl_3D_phidiff_denorm[l][i][j][k] = 0.0;
                     }
                  }
               }
            }
         }
      }
   }
   return;
}

HBT::~HBT()
{

   delete [] KT_array;
   delete [] KT_weight;
   delete [] Kphi;
   delete [] Kphi_weight;

   for(int i = 0; i < Emissionfunction_length; i++)
   {
       for(int j = 0; j < n_KT; j++)
           delete [] emission_S_K[i].data[j];
       delete [] emission_S_K[i].data;
   }
   delete [] emission_S_K;

   delete [] q_out;
   delete [] q_side;
   delete [] q_long;

   if(flag_MC == 1)
   {
       delete [] q_out_MC;
       delete [] q_side_MC;
       delete [] q_long_MC;
       if(azimuthal_flag == 0)
       {
           delete [] Correl_MC_num;
           delete [] Correl_MC_denorm;
       }
       else
       {
           for(int i = 0; i < n_Kphi; i++)
           {
               delete [] Correl_MC_phidiff_num[i];
               delete [] Correl_MC_phidiff_denorm[i];
           }
           delete [] Correl_MC_phidiff_num;
           delete [] Correl_MC_phidiff_denorm;
       }
   }
   else
   {
      if(azimuthal_flag == 0)
      {
         if(flag_1D_projection == 1)
         {
            for(int i = 0; i < ndir; i++)
            {
                delete [] Correl_1D_num[i];
                delete [] Correl_1D_denorm[i];
            }
            delete [] Correl_1D_num;
            delete [] Correl_1D_denorm;
         }
         else
         {
            for(int i=0; i<qnpts; i++)
            {
               for(int j=0; j< qnpts; j++)
               {
                   delete [] Correl_3D_num[i][j];
                   delete [] Correl_3D_denorm[i][j];
               }
               delete [] Correl_3D_num[i];
               delete [] Correl_3D_denorm[i];
            }
            delete [] Correl_3D_num;
            delete [] Correl_3D_denorm;
         }
      }
      else
      {
         if(flag_1D_projection == 1)
         {
            for(int i = 0; i < n_Kphi; i++)
            {
               for(int j = 0; j < ndir; j++)
               {
                   delete [] Correl_1D_phidiff_num[i][j];
                   delete [] Correl_1D_phidiff_denorm[i][j];
               }
               delete [] Correl_1D_phidiff_num[i];
               delete [] Correl_1D_phidiff_denorm[i];
            }
            delete [] Correl_1D_phidiff_num;
            delete [] Correl_1D_phidiff_denorm;
         }
         else
         {
            for(int i = 0; i < n_Kphi; i++)
            {
               for(int j = 0; j < qnpts; j++)
               {
                  for(int k = 0; k < qnpts; k++)
                  {
                     delete [] Correl_3D_phidiff_num[i][j][k];
                     delete [] Correl_3D_phidiff_denorm[i][j][k];
                  }
                  delete [] Correl_3D_phidiff_num[i][j];
                  delete [] Correl_3D_phidiff_denorm[i][j];
               }
               delete [] Correl_3D_phidiff_num[i];
               delete [] Correl_3D_phidiff_denorm[i];
            }
            delete [] Correl_3D_phidiff_num;
            delete [] Correl_3D_phidiff_denorm;
         }
      }
   }
   return;
}

void HBT::calculation_HBT_correlation(double y)
{
   if(kT_differenitial_flag == 1)
   {
       if(azimuthal_flag == 1)
           calculate_azimuthal_dependent_HBT_radii(y);
       else
           calculate_azimuthal_averaged_HBT_radii(y);
   }
   else
   {
       if(azimuthal_flag == 1)
           calculate_azimuthal_dependent_KT_integrated_HBT_radii(y);
       else
           calculate_azimuthal_averaged_KT_integrated_HBT_radii(y);
   }
}

void HBT::calculate_azimuthal_dependent_HBT_radii(double y)
{
   cout << "Calculating "<< particle_ptr[particle_id].name << endl;

   SetEmissionData(FOsurf_ptr, y);

   if(flag_HBT_from_source_variance == 1)
       Cal_HBTRadii_fromEmissionfunction(y);

   for(int iKT = 0; iKT < n_KT; iKT++)
   {
       if(flag_MC == 1)
       {
          Cal_azimuthal_dependent_correlationfunction_MC(iKT, y);
          Output_Correlationfunction_azimuthal_dependent_MC(iKT);
       }
       else
       {
          if(flag_1D_projection == 1)
          {
             Cal_azimuthal_dependent_correlationfunction_1D(iKT, y);
             Output_Correlationfunction_azimuthal_dependent_1D(iKT);
          }
          else
          {
             Cal_azimuthal_dependent_correlationfunction_3D(iKT, y);
             Output_Correlationfunction_azimuthal_dependent_3D(iKT);
          }
       }
   }
}

void HBT::calculate_azimuthal_averaged_HBT_radii(double y)
{
   cout << "Calculating "<< particle_ptr[particle_id].name << endl;

   SetEmissionData(FOsurf_ptr, y);

   for(int iKT = 0; iKT < n_KT; iKT++)
   {
       if(flag_MC == 1)
       {
           Cal_azimuthal_averaged_correlationfunction_MC(iKT, y);
           Output_Correlationfunction_MC(iKT);
       }
       else
       {
           if(flag_1D_projection == 1)
           {
               Cal_azimuthal_averaged_correlationfunction_1D(iKT, y);
               Output_Correlationfunction_1D(iKT);
           }
           else
           {
               Cal_azimuthal_averaged_correlationfunction_3D(iKT, y);
               Output_Correlationfunction_3D(iKT);
           }
       }
   }
}

void HBT::calculate_azimuthal_averaged_KT_integrated_HBT_radii(double y)
{
   cout << "Calculating "<< particle_ptr[particle_id].name << endl;
   
   SetEmissionData(FOsurf_ptr, y);

   if(flag_MC == 1)
   {
       Cal_azimuthal_averaged_KT_inte_correlationfunction_MC(y);
       Output_Correlationfunction_MC(-1);
   }
   else
   {
       if(flag_1D_projection == 1)
       {
           Cal_azimuthal_averaged_KT_inte_correlationfunction_1D(y);
           Output_Correlationfunction_1D(-1);
       }
       else
       {
           Cal_azimuthal_averaged_KT_inte_correlationfunction_3D(y);
           Output_Correlationfunction_3D(-1);
       }
   }
}

void HBT::calculate_azimuthal_dependent_KT_integrated_HBT_radii(double y)
{
   cout << "Calculating "<< particle_ptr[particle_id].name << endl;
   
   SetEmissionData(FOsurf_ptr, y);

   if(flag_MC == 1)
   {
       Cal_azimuthal_dependent_KT_inte_correlationfunction_MC(y);
       Output_Correlationfunction_azimuthal_dependent_MC(-1);
   }
   else
   {
       if(flag_1D_projection == 1)
       {
           Cal_azimuthal_dependent_KT_inte_correlationfunction_1D(y);
           Output_Correlationfunction_azimuthal_dependent_1D(-1);
       }
       else
       {
           Cal_azimuthal_dependent_KT_inte_correlationfunction_3D(y);
           Output_Correlationfunction_azimuthal_dependent_3D(-1);
       }
   }

}

void HBT::SetEmissionData(FO_surf* FO_surface, double K_rap)
// compute emission function at a given pair momentum
{
  double tol = 1e-15;
  double mass = particle_ptr[particle_id].mass;

  double *mT = new double [n_KT];
  double **K_x = new double* [n_KT];
  double **K_y = new double* [n_KT];
  for(int i = 0; i < n_KT; i++)
  {
      double kT_local = KT_array[i];
      mT[i] = sqrt(mass*mass + kT_local*kT_local);
      K_x[i] = new double [n_Kphi];
      K_y[i] = new double [n_Kphi];
      for(int j = 0; j < n_Kphi; j++)
      {
          K_x[i][j] = kT_local*cos(Kphi[j]);
          K_y[i][j] = kT_local*sin(Kphi[j]);
      }
  }

  int idx = 0;
  for(int i=0; i<eta_s_npts; i++)
  {
      double local_eta_s = eta_s[i];
      double ch_localetas = cosh(local_eta_s);
      double sh_localetas = sinh(local_eta_s);

      double ch_y_minus_etas = cosh(K_rap - local_eta_s);
      double sh_y_minus_etas = sinh(K_rap - local_eta_s);
      
      for (int j = 0; j < FO_length; j++)
	{
          for(int ikT = 0; ikT < n_KT; ikT++)
          {
              double K_0 = mT[ikT]*ch_y_minus_etas;
              double K_z = mT[ikT]*sh_y_minus_etas;
              for(int iphi = 0; iphi < n_Kphi; iphi++)
              {
                  double K_x_local = K_x[ikT][iphi];
                  double K_y_local = K_y[ikT][iphi];
                  double S_p = Emissionfunction(K_0, K_x_local, K_y_local, K_z, &FO_surface[j]);
                  if (flag_neg == 1 && S_p < tol)
                  {
                     S_p = 0.0e0;
                  }
	            else
                  {
                     double S_p_withweight = S_p*FO_surface[j].tau*eta_s_weight[i];
                     emission_S_K[idx].data[ikT][iphi] = S_p_withweight;
                  }
              }
          }
          emission_S_K[idx].t = FO_surface[j].tau*ch_localetas;
          emission_S_K[idx].x = FO_surface[j].xpt;
          emission_S_K[idx].y = FO_surface[j].ypt;
          emission_S_K[idx].z = FO_surface[j].tau*sh_localetas;
          idx++;
      }
  }
  Emissionfunction_length = idx;

  // clean up
  delete [] mT;
  for(int i = 0; i < n_KT; i++)
  {
      delete [] K_x[i];
      delete [] K_y[i];
  }
  delete [] K_x;
  delete [] K_y;

  return;
}

double HBT::Emissionfunction(double p0, double px, double py, double pz, FO_surf* surf)
{
   double mu = surf->particle_mu[particle_id];
   double sign = particle_ptr[particle_id].sign;
   double degen = particle_ptr[particle_id].gspin;
   double mass = particle_ptr[particle_id].mass;

   double gammaT = surf->u0;
   double ux = surf->u1;
   double uy = surf->u2;
   double Tdec = surf->Tdec;
   double Pdec = surf->Pdec;
   double Edec = surf->Edec;
   double da0 = surf->da0;
   double da1 = surf->da1;
   double da2 = surf->da2;
   double pi00 = surf->pi00;
   double pi01 = surf->pi01;
   double pi02 = surf->pi02;
   double pi11 = surf->pi11;
   double pi12 = surf->pi12;
   double pi22 = surf->pi22;
   double pi33 = surf->pi33;

   double E_over_T = (p0*gammaT - px*ux - py*uy)/Tdec;
   double expon = E_over_T - mu/Tdec;
   //double expon = ((p0*gammaT - px*ux - py*uy) - mu) / Tdec;
   double f0 = 1./(exp(expon)+sign);       //thermal equilibrium distributions

   //p^mu d^3sigma_mu: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
   double pdsigma = p0*da0 + px*da1 + py*da2;

   //viscous corrections
   double delta_f_shear = 0.0;
   double delta_f_bulk = 0.0;
   if(INCLUDE_SHEAR_DELTAF)
   {
       double Wfactor = p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33;
       delta_f_shear = (1. - sign*f0)*Wfactor/(2.0*Tdec*Tdec*(Edec+Pdec));
   }
   if(INCLUDE_BULK_DELTAF)
   {
       double bulkPi = surf->bulkPi/hbarC;   // convert it to fm^-4
       // parameterization from JF
       double C_bulk, e2;

       // A Polynomial fit to each coefficient -- X is the temperature in fm^-1
       // Both fits are reliable between T=100 -- 180 MeV , do not trust it beyond
       double Tfm = Tdec/hbarC;  // convert it to fm^-1
       double T_power[11];
       T_power[0] = 1.0;
       for(int ipow = 1; ipow < 11; ipow++)
           T_power[ipow] = T_power[ipow-1]*Tfm;

       C_bulk = (  642096.624265727 - 8163329.49562861*T_power[1] 
                 + 47162768.4292073*T_power[2] - 162590040.002683*T_power[3] 
                 + 369637951.096896*T_power[4] - 578181331.809836*T_power[5] 
                 + 629434830.225675*T_power[6] - 470493661.096657*T_power[7] 
                 + 230936465.421*T_power[8] - 67175218.4629078*T_power[9] 
                 + 8789472.32652964*T_power[10]);

       e2 = (  1.18171174036192 - 17.6740645873717*T_power[1] 
             + 136.298469057177*T_power[2] - 635.999435106846*T_power[3] 
             + 1918.77100633321*T_power[4] - 3836.32258307711*T_power[5]
             + 5136.35746882372*T_power[6] - 4566.22991441914*T_power[7]
             + 2593.45375240886*T_power[8] - 853.908199724349*T_power[9]
             + 124.260460450113*T_power[10]);

      // bulk delta f is
      delta_f_bulk = -1.0*(1.-sign*f0)/E_over_T*C_bulk*(mass*mass/Tdec/Tdec/3. - e2*E_over_T*E_over_T)*bulkPi;
   }

   double dN_dyd2pTdphi;
   //if ((deltaf + delta_f_bulk) < -1.0)  // delta f correction is too large
   //   dN_dyd2pTdphi = 0.0;
   //else
   dN_dyd2pTdphi = 1.0*degen/(8.0*(M_PI*M_PI*M_PI))*pdsigma*f0*(1. + delta_f_shear + delta_f_bulk);
   //out << "Spectral funct = " << dN_dyd2pTdphi << endl;

   return (dN_dyd2pTdphi);
}

void HBT::Cal_HBTRadii_fromEmissionfunction(double K_y)
{
  for(int ikT = 0; ikT < n_KT; ikT++)
  {
      double KT_local = KT_array[ikT];
      for(int iphi = 0; iphi < n_Kphi; iphi++)
      {
          double K_phi_local = Kphi[iphi];
          double* resultsX = new double[15];
          for(int i = 0; i < 15; i++)
             resultsX[i] = 0.0e0;

          for(int i=0; i < Emissionfunction_length; i++)
          {
             double S_p = emission_S_K[i].data[ikT][iphi];
             double tpt = emission_S_K[i].t;
             double xpt = emission_S_K[i].x;
             double ypt = emission_S_K[i].y;
             double zpt = emission_S_K[i].z;
 
             for(int ii=0; ii<2; ii++) // assuming reflection symmetry in eta_s
             {
                zpt = zpt*(-1);
                resultsX[0]  += S_p;             //single particle spectra <1>
                resultsX[1]  += S_p*xpt;         //<x>
                resultsX[2]  += S_p*ypt;         //<y>
                resultsX[3]  += S_p*zpt;         //<z>
                resultsX[4]  += S_p*xpt*ypt;     //<xy>
                resultsX[5]  += S_p*xpt*xpt;     //<xx>
                resultsX[6]  += S_p*ypt*ypt;     //<yy>
                resultsX[7]  += S_p*tpt;         //<t>
                resultsX[8]  += S_p*tpt*xpt;     //<tx>
                resultsX[9]  += S_p*tpt*ypt;     //<ty>
                resultsX[10] += S_p*tpt*zpt;     //<tz>
                resultsX[11] += S_p*zpt*zpt;     //<zz>
                resultsX[12] += S_p*tpt*tpt;     //<tt>
                resultsX[13] += S_p*xpt*zpt;     //<xz>
                resultsX[14] += S_p*ypt*zpt;     //<yz>
             }
          }
                                                  	
          for(int i=0; i<15; i++)     //change to correct unit
             resultsX[i] = resultsX[i]/hbarC/hbarC/hbarC;
          
          double spectra = resultsX[0];
          double meanx  = resultsX[1];
          double meany  = resultsX[2];
          double meanz  = resultsX[3];
          double meanxy = resultsX[4];
          double meanxx = resultsX[5];
          double meanyy = resultsX[6];
          double meant  = resultsX[7];
          double meanxt = resultsX[8];
          double meanyt = resultsX[9];
          double meanzt = resultsX[10];
          double meanzz = resultsX[11];
          double meantt = resultsX[12];
          double meanxz = resultsX[13];
          double meanyz = resultsX[14];
          
          //calculate the components of S^{\mu\nu} tensor
          double S00 = meantt/spectra - meant/spectra * meant/spectra;
          double S01 = meanxt/spectra - meanx/spectra * meant/spectra;
          double S02 = meanyt/spectra - meany/spectra * meant/spectra;
          double S03 = meanzt/spectra - meanz/spectra * meant/spectra;
          double S11 = meanxx/spectra - meanx/spectra * meanx/spectra;
          double S12 = meanxy/spectra - meanx/spectra * meany/spectra;
          double S13 = meanxz/spectra - meanx/spectra * meanz/spectra;
          double S22 = meanyy/spectra - meany/spectra * meany/spectra;
          double S23 = meanyz/spectra - meany/spectra * meanz/spectra;
          double S33 = meanzz/spectra - meanz/spectra * meanz/spectra;
          
          //calculate HBT radii from single particle emission function
          double mass = particle_ptr[particle_id].mass;
          double m_T = sqrt(mass*mass + KT_local*KT_local);
          double beta_T = KT_local/m_T;
          double beta_L = 0.5*log((1+K_y)/(1-K_y));

          double R_out2 = S11*cos(K_phi_local)*cos(K_phi_local) + S22*sin(K_phi_local)*sin(K_phi_local) 
                          + S12*sin(2*K_phi_local) - 2*beta_T*(S01*cos(K_phi_local) 
                          + S02*sin(K_phi_local))+ beta_T*beta_T*S00;          //R_out^2
          double R_side2 = S11*sin(K_phi_local)*sin(K_phi_local) + S22*cos(K_phi_local)*cos(K_phi_local) 
                           - S12*sin(2*K_phi_local);                           //R_side^2
          double R_long2 = S33 - 2*beta_L*S03 + beta_L*beta_L*S00;       //R_long^2

          R_out_EM = sqrt(R_out2);
          R_side_EM = sqrt(R_side2);
          R_long_EM = sqrt(R_long2);

          cout << "dN/(dypTdpTdphi) = " << scientific << setw(12) << setprecision(6) << spectra << "  1/(GeV^2)." << endl;
          cout << "R_out = " << R_out_EM << "  fm." << endl;
          cout << "R_side = " << R_side_EM << "  fm." << endl;
          cout << "R_long = " << R_long_EM << "  fm." << endl;

          delete [] resultsX;
      }
  }
  return;
}                                         	

void HBT::Cal_azimuthal_averaged_correlationfunction_1D(int iKT, double K_y)
{
   double hbarC_inv = 1./hbarC;
   if(fabs(K_y) > 1e-16)
   {
       cout<<"HBT:: not support for y is not equal to 0 yet!" << endl;
       return;
   }
   
   double K_T = KT_array[iKT];
   cout << "generating the 1d slices of the correlation function along q_out, q_side, and q_long direction for K_T = " << K_T << " GeV..." << endl;
   double mass = particle_ptr[particle_id].mass;

   double *cosK_phi = new double [n_Kphi];
   double *sinK_phi = new double [n_Kphi];
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
      cosK_phi[iphi] = cos(Kphi[iphi]);
      sinK_phi[iphi] = sin(Kphi[iphi]);
   }

   double spectra = 0.0;

   for(int k = 0; k < Emissionfunction_length; k++)
   {
      for(int iphi = 0; iphi < n_Kphi; iphi++)
      {
         double ss = emission_S_K[k].data[iKT][iphi]*Kphi_weight[iphi];
         spectra += ss*2;
      }
   }

   double local_q_out, local_q_side, local_q_long;
   for(int i = 0; i < qnpts; i++)
   {
      cout << "calculating q_mu = " << q_out[i] << " GeV..." << endl;
      for (int l = 0; l < ndir; l++)
      {
         switch (l)
         {
            case 0:
            {
               local_q_out  = q_out[i];
               local_q_side = 0.0e0;
               local_q_long = 0.0e0;
               break;
            }
   	      case 1:
            {
               local_q_out  = 0.0e0;
               local_q_side = q_side[i];
               local_q_long = 0.0e0;
               break;
            }
            case 2:
            {
               local_q_out  = 0.0e0;
               local_q_side = 0.0e0;
               local_q_long = q_long[i];
               break;
            }
            case 3:
            {
               local_q_out  = q_out[i];
               local_q_side = q_side[i];
               local_q_long = 0.0;
               break;
            }
            case 4:
            {
               local_q_out  = q_out[i];
               local_q_side = -q_side[i];
               local_q_long = 0.0;
               break;
            }
            default:
            {
               cout << "error in assigning q values! "<< endl;
               break;
            }
         }

     	   double xsi  = K_T*K_T + mass*mass + (local_q_out*local_q_out + local_q_side*local_q_side + local_q_long*local_q_long)/4.0;  //Set Xsi
         double E1sq = xsi + K_T*local_q_out;
         double E2sq = xsi - K_T*local_q_out;
         double qt = sqrt(E1sq) - sqrt(E2sq);
         double qz = local_q_long;

         double integ1 = 0.0;  // numerator cosine part
         double integ2 = 0.0;  // numerator sine part

         for(int k = 0; k < Emissionfunction_length; k++)
         {
            double tpt = emission_S_K[k].t;
            double xpt = emission_S_K[k].x;
            double ypt = emission_S_K[k].y;
            double zpt = emission_S_K[k].z;

            for(int iphi = 0; iphi < n_Kphi; iphi++)
            {
                double ss = emission_S_K[k].data[iKT][iphi]*Kphi_weight[iphi];
                double qx = local_q_out*cosK_phi[iphi] - local_q_side*sinK_phi[iphi];
                double qy = local_q_side*cosK_phi[iphi] + local_q_out*sinK_phi[iphi];
                
                double temp_arg = tpt*qt - qx*xpt - qy*ypt;
                for(int ii=0; ii<2; ii++)
                {
                   zpt = zpt*(-1);   //using the symmetry along z axis
                   double arg = (temp_arg - qz*zpt)*hbarC_inv;
                   integ1 += cos(arg)*ss;
                   integ2 += sin(arg)*ss;
                }
             }
         }
         double localvalue = integ1*integ1+integ2*integ2;
         Correl_1D_num[l][i]  = localvalue;
         Correl_1D_denorm[l][i]  = spectra*spectra;
      }
   }

   delete [] cosK_phi;
   delete [] sinK_phi;

   return;
} 

void HBT::Cal_azimuthal_averaged_KT_inte_correlationfunction_1D(double K_y)
{
   if(fabs(K_y) > 1e-16)
   {
       cout<<"HBT:: not support for y is not equal to 0 yet!" << endl;
       return;
   }
   
   cout << "generating the 1d slices of the correlation function along q_out, q_side, and q_long direction ..." << endl;
   double mass = particle_ptr[particle_id].mass;

   double *cosK_phi = new double [n_Kphi];
   double *sinK_phi = new double [n_Kphi];
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
      cosK_phi[iphi] = cos(Kphi[iphi]);
      sinK_phi[iphi] = sin(Kphi[iphi]);
   }

   // calculate particle spectra
   double spectra = 0.0;
   for(int k = 0; k < Emissionfunction_length; k++)
   {
      for(int iKT = 0; iKT < n_KT; iKT++)
      {
         double KT_local = KT_array[iKT];
         double KT_weight_local = KT_weight[iKT];
         for(int iphi = 0; iphi < n_Kphi; iphi++)
         {
            double ss = emission_S_K[k].data[iKT][iphi]*Kphi_weight[iphi]*KT_local*KT_weight_local;
            spectra += ss*2;
         }
      }
   }
 
   // calculate correlation function
   double local_q_out, local_q_side, local_q_long;
   for(int i = 0; i < qnpts; i++)
   {
      cout << "calculating q_mu = " << q_out[i] << " GeV..." << endl;
      for (int l = 0; l < ndir; l++)
      {
         switch (l)
         {
            case 0:
            {
               local_q_out  = q_out[i];
               local_q_side = 0.0e0;
               local_q_long = 0.0e0;
               break;
            }
   	      case 1:
            {
               local_q_out  = 0.0e0;
               local_q_side = q_side[i];
               local_q_long = 0.0e0;
               break;
            }
            case 2:
            {
               local_q_out  = 0.0e0;
               local_q_side = 0.0e0;
               local_q_long = q_long[i];
               break;
            }
            case 3:
            {
               local_q_out  = q_out[i];
               local_q_side = q_side[i];
               local_q_long = 0.0;
               break;
            }
            case 4:
            {
               local_q_out  = q_out[i];
               local_q_side = -q_side[i];
               local_q_long = 0.0;
               break;
            }
            default:
            {
               cout << "error in assigning q values! "<< endl;
               break;
            }
         }

         double integ1 = 0.0;  // numerator cosine part
         double integ2 = 0.0;  // numerator sine part

         for(int iKT = 0; iKT < n_KT; iKT++)
         {
             double K_T = KT_array[iKT];
             double KT_weight_local = KT_weight[iKT];
     	       double xsi  = K_T*K_T + mass*mass + (local_q_out*local_q_out + local_q_side*local_q_side + local_q_long*local_q_long)/4.0;  //Set Xsi
             double E1sq = xsi + K_T*local_q_out;
             double E2sq = xsi - K_T*local_q_out;
             double qt = sqrt(E1sq) - sqrt(E2sq);
             double qz = local_q_long;

             for(int iphi = 0; iphi < n_Kphi; iphi++)
             {
                 double qx = local_q_out*cosK_phi[iphi] - local_q_side*sinK_phi[iphi];
                 double qy = local_q_side*cosK_phi[iphi] + local_q_out*sinK_phi[iphi];

                 for(int k = 0; k < Emissionfunction_length; k++)
                 {
                    double ss  = emission_S_K[k].data[iKT][iphi]*Kphi_weight[iphi]*K_T*KT_weight_local;
                    double tpt = emission_S_K[k].t;
                    double xpt = emission_S_K[k].x;
                    double ypt = emission_S_K[k].y;
                    double zpt = emission_S_K[k].z;
                    
                    for(int ii=0; ii<2; ii++)
                    {
                       zpt = zpt*(-1);   //using the symmetry along z axis
                       double arg = (tpt*qt - (qx*xpt + qy*ypt + qz*zpt))/hbarC;
                       integ1 += cos(arg)*ss;
                       integ2 += sin(arg)*ss;
                    }
                 }
             }
         }
         double localvalue = integ1*integ1+integ2*integ2;
         Correl_1D_num[l][i]  = localvalue;
         Correl_1D_denorm[l][i]  = spectra*spectra;
      }
   }

   delete [] cosK_phi;
   delete [] sinK_phi;

   return;
} 

void HBT::Cal_azimuthal_dependent_correlationfunction_1D(int iKT, double K_y)
{
   if(fabs(K_y) > 1e-16)
   {
       cout<<"HBT:: not support for y is not equal to 0 yet!" << endl;
       return;
   }
   
   double K_T = KT_array[iKT];
   cout << "generating the 1d slices of the correlation function along q_out, q_side, and q_long direction for K_T = " << K_T << " GeV..." << endl;
   double mass = particle_ptr[particle_id].mass;

   double *cosK_phi = new double [n_Kphi];
   double *sinK_phi = new double [n_Kphi];
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
      cosK_phi[iphi] = cos(Kphi[iphi]);
      sinK_phi[iphi] = sin(Kphi[iphi]);
   }

   double *spectra = new double [n_Kphi];
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
      double ss = 0.0;
      for(int k = 0; k < Emissionfunction_length; k++)
         ss += emission_S_K[k].data[iKT][iphi];
      spectra[iphi] = ss*2;
   }

   for(int i = 0; i < qnpts; i++)
   {
      cout << "calculating q_mu = " << q_out[i] << " GeV..." << endl;
      for (int l = 0; l < ndir; l++)
      {
         double local_q_out=0.0, local_q_side=0.0, local_q_long=0.0;
         switch (l)
         {
            case 0:
            {
               local_q_out  = q_out[i];
               local_q_side = 0.0e0;
               local_q_long = 0.0e0;
               break;
            }
   	      case 1:
            {
               local_q_out  = 0.0e0;
               local_q_side = q_side[i];
               local_q_long = 0.0e0;
               break;
            }
            case 2:
            {
               local_q_out  = 0.0e0;
               local_q_side = 0.0e0;
               local_q_long = q_long[i];
               break;
            }
            case 3:
            {
               local_q_out  = q_out[i];
               local_q_side = q_side[i];
               local_q_long = 0.0e0;
               break;
            }
            case 4:
            {
               local_q_out  = q_out[i];
               local_q_side = -q_side[i];
               local_q_long = 0.0e0;
               break;
            }
            default:
            {
               cout << "error in assigning q values! "<< endl;
               break;
            }
         }

     	   double xsi  = K_T*K_T + mass*mass + (local_q_out*local_q_out + local_q_side*local_q_side + local_q_long*local_q_long)/4.0;  //Set Xsi
         double E1sq = xsi + K_T*local_q_out;
         double E2sq = xsi - K_T*local_q_out;
         double qt = sqrt(E1sq) - sqrt(E2sq);
         double qz = local_q_long;


         for(int iphi = 0; iphi < n_Kphi; iphi++)
         {
             double qx = local_q_out*cosK_phi[iphi] - local_q_side*sinK_phi[iphi];
             double qy = local_q_side*cosK_phi[iphi] + local_q_out*sinK_phi[iphi];
             
             double integ1 = 0.0;  // numerator cosine part
             double integ2 = 0.0;  // numerator sine part

             for(int k = 0; k < Emissionfunction_length; k++)
             {
                double ss  = emission_S_K[k].data[iKT][iphi];
                double tpt = emission_S_K[k].t;
                double xpt = emission_S_K[k].x;
                double ypt = emission_S_K[k].y;
                double zpt = emission_S_K[k].z;
                
                for(int ii=0; ii<2; ii++)
                {
                   zpt = zpt*(-1);   //using the symmetry along z axis
                   double arg = (tpt*qt - (qx*xpt + qy*ypt + qz*zpt))/hbarC;
                   integ1 += cos(arg)*ss;
                   integ2 += sin(arg)*ss;
                }
             }
             double localvalue = integ1*integ1+integ2*integ2;
             Correl_1D_phidiff_num[iphi][l][i] = localvalue;
             Correl_1D_phidiff_denorm[iphi][l][i] = spectra[iphi]*spectra[iphi];
         }
      }
   }

   delete [] cosK_phi;
   delete [] sinK_phi;
   delete [] spectra;

   return;
} 

void HBT::Cal_azimuthal_dependent_KT_inte_correlationfunction_1D(double K_y)
{
   if(fabs(K_y) > 1e-16)
   {
       cout<<"HBT:: not support for y is not equal to 0 yet!" << endl;
       return;
   }
   
   cout << "generating the 1d slices of the correlation function along q_out, q_side, and q_long direction ..." << endl;
   double mass = particle_ptr[particle_id].mass;

   double *cosK_phi = new double [n_Kphi];
   double *sinK_phi = new double [n_Kphi];
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
      cosK_phi[iphi] = cos(Kphi[iphi]);
      sinK_phi[iphi] = sin(Kphi[iphi]);
   }

   double *spectra = new double [n_Kphi];
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
      double ss = 0.0;
      for(int iKT = 0; iKT < n_KT; iKT++)
      {
          double KT_local = KT_array[iKT];
          double KT_weight_local = KT_weight[iKT];
          for(int k = 0; k < Emissionfunction_length; k++)
              ss += emission_S_K[k].data[iKT][iphi]*KT_local*KT_weight_local;
      }
      spectra[iphi] = ss*2;
   }

   for(int i = 0; i < qnpts; i++)
   {
      cout << "calculating q_mu = " << q_out[i] << " GeV..." << endl;
      for (int l = 0; l < ndir; l++)
      {
         double local_q_out=0.0, local_q_side=0.0, local_q_long=0.0;
         switch (l)
         {
            case 0:
            {
               local_q_out  = q_out[i];
               local_q_side = 0.0e0;
               local_q_long = 0.0e0;
               break;
            }
   	      case 1:
            {
               local_q_out  = 0.0e0;
               local_q_side = q_side[i];
               local_q_long = 0.0e0;
               break;
            }
            case 2:
            {
               local_q_out  = 0.0e0;
               local_q_side = 0.0e0;
               local_q_long = q_long[i];
               break;
            }
            case 3:
            {
               local_q_out  = q_out[i];
               local_q_side = q_side[i];
               local_q_long = 0.0e0;
               break;
            }
            case 4:
            {
               local_q_out  = q_out[i];
               local_q_side = -q_side[i];
               local_q_long = 0.0e0;
               break;
            }
            default:
            {
               cout << "error in assigning q values! "<< endl;
               break;
            }
         }

         for(int iphi = 0; iphi < n_Kphi; iphi++)
         {
             double qx = local_q_out*cosK_phi[iphi] - local_q_side*sinK_phi[iphi];
             double qy = local_q_side*cosK_phi[iphi] + local_q_out*sinK_phi[iphi];
             
             double integ1 = 0.0;  // numerator cosine part
             double integ2 = 0.0;  // numerator sine part

             for(int iKT = 0; iKT < n_KT; iKT++)
             {
                 double K_T = KT_array[iKT];
                 double KT_weight_local = KT_weight[iKT];
     	           double xsi  = K_T*K_T + mass*mass + (local_q_out*local_q_out + local_q_side*local_q_side + local_q_long*local_q_long)/4.0;  //Set Xsi
                 double E1sq = xsi + K_T*local_q_out;
                 double E2sq = xsi - K_T*local_q_out;
                 double qt = sqrt(E1sq) - sqrt(E2sq);
                 double qz = local_q_long;

                 for(int k = 0; k < Emissionfunction_length; k++)
                 {
                    double ss  = emission_S_K[k].data[iKT][iphi]*K_T*KT_weight_local;
                    double tpt = emission_S_K[k].t;
                    double xpt = emission_S_K[k].x;
                    double ypt = emission_S_K[k].y;
                    double zpt = emission_S_K[k].z;
                    
                    for(int ii=0; ii<2; ii++)
                    {
                       zpt = zpt*(-1);   //using the symmetry along z axis
                       double arg = (tpt*qt - (qx*xpt + qy*ypt + qz*zpt))/hbarC;
                       integ1 += cos(arg)*ss;
                       integ2 += sin(arg)*ss;
                    }
                 }
             }
            double localvalue = integ1*integ1+integ2*integ2;
            Correl_1D_phidiff_num[iphi][l][i] = localvalue;
            Correl_1D_phidiff_denorm[iphi][l][i] = spectra[iphi]*spectra[iphi];
         }
      }
   }

   delete [] cosK_phi;
   delete [] sinK_phi;
   delete [] spectra;

   return;
} 

void HBT::Cal_azimuthal_averaged_correlationfunction_3D(int iKT, double K_y)
{
   if(fabs(K_y) > 1e-16)
   {
       cout<<"not support for y not equals 0 yet!" << endl;
       return;
   }

   double K_T = KT_array[iKT];
   
   cout << "generating correlation function in 3D for K_T = " 
        << K_T << " GeV ... " << endl;

   double mass = particle_ptr[particle_id].mass;
   
   double *cosK_phi = new double [n_Kphi];
   double *sinK_phi = new double [n_Kphi];
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
      cosK_phi[iphi] = cos(Kphi[iphi]);
      sinK_phi[iphi] = sin(Kphi[iphi]);
   }
   
   double spectra = 0.0e0;
   for(int k = 0; k < Emissionfunction_length; k++)
   {
      for(int iphi = 0; iphi < n_Kphi; iphi++)
      {
         double ss  = emission_S_K[k].data[iKT][iphi]*Kphi_weight[iphi];
         spectra += ss*2;
      }
   }

   for(int i = 0; i < qnpts; i++)  // q_out loop
   {
      double local_q_out = q_out[i];
      for(int j = 0; j < qnpts; j++)  // q_side loop
      {
         double local_q_side = q_side[j];
         for(int k = 0; k < qnpts; k++)  // q_long loop
         {
            double local_q_long = q_long[k];
            cout << "q_out = " << local_q_out << " GeV, "
                 << "q_side = " << local_q_side << " GeV, "
                 << "q_long = " << local_q_long << " GeV... " << endl;

            double integ1 = 0.0;                         
            double integ2 = 0.0;
            double sum = 0.0;

     	      double xsi = K_T*K_T + mass*mass + (local_q_out*local_q_out + local_q_side*local_q_side + local_q_long*local_q_long)/4.0;  //Set Xsi
            double E1sq = xsi + K_T*local_q_out;
            double E2sq = xsi - K_T*local_q_out;
            double qt = sqrt(E1sq) - sqrt(E2sq);
            double qz = local_q_long;

            for(int iphi = 0; iphi < n_Kphi; iphi++)
            {
               double qx = local_q_out*cosK_phi[iphi] - local_q_side*sinK_phi[iphi];
               double qy = local_q_side*cosK_phi[iphi] + local_q_out*sinK_phi[iphi];
               
               for(int m = 0; m < Emissionfunction_length; m++)
               {
                  double ss = emission_S_K[m].data[iKT][iphi]*Kphi_weight[iphi];
                  double tpt = emission_S_K[m].t;
                  double xpt = emission_S_K[m].x;
                  double ypt = emission_S_K[m].y;
                  double zpt = emission_S_K[m].z;

                  for(int ii=0; ii<2; ii++)
                  {
                     zpt = zpt*(-1);
                     double arg = (tpt*qt - (qx*xpt + qy*ypt + qz*zpt))/hbarC;
                     integ1 += cos(arg)*ss;
                     integ2 += sin(arg)*ss;
                  }
               }
            }
            sum = integ1*integ1+integ2*integ2;
            Correl_3D_num[i][j][k] = sum;
            Correl_3D_denorm[i][j][k] = spectra*spectra;
         }
      }
   }

   delete [] cosK_phi;
   delete [] sinK_phi;

   return;
}

void HBT::Cal_azimuthal_averaged_KT_inte_correlationfunction_3D(double K_y)
{
   if(fabs(K_y) > 1e-16)
   {
       cout<<"not support for y not equals 0 yet!" << endl;
       return;
   }

   cout << "generating correlation function in 3D  ... " << endl;

   double mass = particle_ptr[particle_id].mass;
   
   double *cosK_phi = new double [n_Kphi];
   double *sinK_phi = new double [n_Kphi];
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
      cosK_phi[iphi] = cos(Kphi[iphi]);
      sinK_phi[iphi] = sin(Kphi[iphi]);
   }
   
   double spectra = 0.0e0;
   for(int k = 0; k < Emissionfunction_length; k++)
   {
      for(int iKT = 0; iKT < n_KT; iKT++)
      {
         for(int iphi = 0; iphi < n_Kphi; iphi++)
         {
            double ss  = emission_S_K[k].data[iKT][iphi]*Kphi_weight[iphi]*KT_array[iKT]*KT_weight[iKT];
            spectra += ss*2;
         }
      }
   }

   for(int i = 0; i < qnpts; i++)  // q_out loop
   {
      double local_q_out = q_out[i];
      for(int j = 0; j < qnpts; j++)  // q_side loop
      {
         double local_q_side = q_side[j];
         for(int k = 0; k < qnpts; k++)  // q_long loop
         {
            double local_q_long = q_long[k];
            cout << "q_out = " << local_q_out << " GeV, "
                 << "q_side = " << local_q_side << " GeV, "
                 << "q_long = " << local_q_long << " GeV... " << endl;

            double integ1 = 0.0;                         
            double integ2 = 0.0;
            double sum = 0.0;

            for(int iKT = 0; iKT < n_KT; iKT++)
            {
                double K_T = KT_array[iKT];
                double KT_weight_local = KT_weight[iKT];
     	          double xsi = K_T*K_T + mass*mass + (local_q_out*local_q_out + local_q_side*local_q_side + local_q_long*local_q_long)/4.0;  //Set Xsi
                double E1sq = xsi + K_T*local_q_out;
                double E2sq = xsi - K_T*local_q_out;
                double qt = sqrt(E1sq) - sqrt(E2sq);
                double qz = local_q_long;

                for(int iphi = 0; iphi < n_Kphi; iphi++)
                {
                   double qx = local_q_out*cosK_phi[iphi] - local_q_side*sinK_phi[iphi];
                   double qy = local_q_side*cosK_phi[iphi] + local_q_out*sinK_phi[iphi];
                   
                   for(int m = 0; m < Emissionfunction_length; m++)
                   {
                      double ss = emission_S_K[m].data[iKT][iphi]*Kphi_weight[iphi]*K_T*KT_weight_local;
                      double tpt = emission_S_K[m].t;
                      double xpt = emission_S_K[m].x;
                      double ypt = emission_S_K[m].y;
                      double zpt = emission_S_K[m].z;

                      for(int ii=0; ii<2; ii++)
                      {
                         zpt = zpt*(-1);
                         double arg = (tpt*qt - (qx*xpt + qy*ypt + qz*zpt))/hbarC;
                         integ1 += cos(arg)*ss;
                         integ2 += sin(arg)*ss;
                      }
                   }
                }
            }
            sum = integ1*integ1+integ2*integ2;
            Correl_3D_num[i][j][k] = sum;
            Correl_3D_denorm[i][j][k] = spectra*spectra;
         }
      }
   }

   delete [] cosK_phi;
   delete [] sinK_phi;

   return;
}

void HBT::Cal_azimuthal_averaged_correlationfunction_MC(int iKT, double K_y)
{
   if(fabs(K_y) > 1e-16)
   {
       cout<<"not support for y not equals 0 yet!" << endl;
       return;
   }

   double K_T = KT_array[iKT];
   
   cout << "generating correlation function in 3D for K_T = " 
        << K_T << " GeV ... " << endl;

   double mass = particle_ptr[particle_id].mass;
   
   double *cosK_phi = new double [n_Kphi];
   double *sinK_phi = new double [n_Kphi];
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
      cosK_phi[iphi] = cos(Kphi[iphi]);
      sinK_phi[iphi] = sin(Kphi[iphi]);
   }
   
   double spectra = 0.0e0;
   for(int k = 0; k < Emissionfunction_length; k++)
   {
      for(int iphi = 0; iphi < n_Kphi; iphi++)
      {
         double ss  = emission_S_K[k].data[iKT][iphi]*Kphi_weight[iphi];
         spectra += ss*2;
      }
   }

   double local_q_out, local_q_side, local_q_long;
   for(int i = 0; i < MC_samples; i++)
   {
       if(i < MC_samples*2/3)
       {
           local_q_out = drand48()*q_max;
           local_q_side = drand48()*q_max;
           local_q_long = 0.0;
       }
       else
       {
           local_q_out = 0.0;
           local_q_side = 0.0;
           local_q_long = drand48()*q_max;
       }

       q_out_MC[i] = local_q_out;
       q_side_MC[i] = local_q_side;
       q_long_MC[i] = local_q_long;

       cout << "q_out = " << local_q_out << " GeV, "
            << "q_side = " << local_q_side << " GeV, "
            << "q_long = " << local_q_long << " GeV... " << endl;

       double integ1 = 0.0;                         
       double integ2 = 0.0;
       double sum = 0.0;

     	 double xsi = K_T*K_T + mass*mass + (local_q_out*local_q_out + local_q_side*local_q_side + local_q_long*local_q_long)/4.0;  //Set Xsi
       double E1sq = xsi + K_T*local_q_out;
       double E2sq = xsi - K_T*local_q_out;
       double qt = sqrt(E1sq) - sqrt(E2sq);
       double qz = local_q_long;

       for(int iphi = 0; iphi < n_Kphi; iphi++)
       {
          double qx = local_q_out*cosK_phi[iphi] - local_q_side*sinK_phi[iphi];
          double qy = local_q_side*cosK_phi[iphi] + local_q_out*sinK_phi[iphi];
          
          for(int m = 0; m < Emissionfunction_length; m++)
          {
             double ss = emission_S_K[m].data[iKT][iphi]*Kphi_weight[iphi];
             double tpt = emission_S_K[m].t;
             double xpt = emission_S_K[m].x;
             double ypt = emission_S_K[m].y;
             double zpt = emission_S_K[m].z;

             for(int ii=0; ii<2; ii++)
             {
                zpt = zpt*(-1);
                double arg = (tpt*qt - (qx*xpt + qy*ypt + qz*zpt))/hbarC;
                integ1 += cos(arg)*ss;
                integ2 += sin(arg)*ss;
             }
          }
       }
       Correl_MC_num[i] = integ1*integ1+integ2*integ2;
       Correl_MC_denorm[i] = spectra*spectra;
   }

   delete [] cosK_phi;
   delete [] sinK_phi;

   return;
}

void HBT::Cal_azimuthal_averaged_KT_inte_correlationfunction_MC(double K_y)
{
   if(fabs(K_y) > 1e-16)
   {
       cout<<"not support for y not equals 0 yet!" << endl;
       return;
   }

   
   cout << "generating correlation function in MC ... " << endl;

   double mass = particle_ptr[particle_id].mass;
   
   double *cosK_phi = new double [n_Kphi];
   double *sinK_phi = new double [n_Kphi];
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
      cosK_phi[iphi] = cos(Kphi[iphi]);
      sinK_phi[iphi] = sin(Kphi[iphi]);
   }
   
   double spectra = 0.0e0;
   for(int k = 0; k < Emissionfunction_length; k++)
   {
      for(int iKT = 0; iKT < n_KT; iKT++)
      {
         for(int iphi = 0; iphi < n_Kphi; iphi++)
         {
            double ss  = emission_S_K[k].data[iKT][iphi]*Kphi_weight[iphi]*KT_array[iKT]*KT_weight[iKT];
            spectra += ss*2;
         }
      }
   }

   double local_q_out, local_q_side, local_q_long;
   for(int i = 0; i < MC_samples; i++)
   {
       if(i < MC_samples*2/3)
       {
           local_q_out = drand48()*q_max;
           local_q_side = drand48()*q_max;
           local_q_long = 0.0;
       }
       else
       {
           local_q_out = 0.0;
           local_q_side = 0.0;
           local_q_long = drand48()*q_max;
       }

       q_out_MC[i] = local_q_out;
       q_side_MC[i] = local_q_side;
       q_long_MC[i] = local_q_long;

       cout << "q_out = " << local_q_out << " GeV, "
            << "q_side = " << local_q_side << " GeV, "
            << "q_long = " << local_q_long << " GeV... " << endl;

       double integ1 = 0.0;                         
       double integ2 = 0.0;
       double sum = 0.0;

       for(int iKT = 0; iKT < n_KT; iKT++)
       {
           double K_T = KT_array[iKT];
           double KT_weight_local = KT_weight[iKT];
     	     double xsi = K_T*K_T + mass*mass + (local_q_out*local_q_out + local_q_side*local_q_side + local_q_long*local_q_long)/4.0;  //Set Xsi
           double E1sq = xsi + K_T*local_q_out;
           double E2sq = xsi - K_T*local_q_out;
           double qt = sqrt(E1sq) - sqrt(E2sq);
           double qz = local_q_long;

           for(int iphi = 0; iphi < n_Kphi; iphi++)
           {
              double qx = local_q_out*cosK_phi[iphi] - local_q_side*sinK_phi[iphi];
              double qy = local_q_side*cosK_phi[iphi] + local_q_out*sinK_phi[iphi];
              
              for(int m = 0; m < Emissionfunction_length; m++)
              {
                 double ss = emission_S_K[m].data[iKT][iphi]*Kphi_weight[iphi]*K_T*KT_weight_local;
                 double tpt = emission_S_K[m].t;
                 double xpt = emission_S_K[m].x;
                 double ypt = emission_S_K[m].y;
                 double zpt = emission_S_K[m].z;

                 for(int ii=0; ii<2; ii++)
                 {
                    zpt = zpt*(-1);
                    double arg = (tpt*qt - (qx*xpt + qy*ypt + qz*zpt))/hbarC;
                    integ1 += cos(arg)*ss;
                    integ2 += sin(arg)*ss;
                 }
              }
           }
       }
       Correl_MC_num[i] = integ1*integ1+integ2*integ2;
       Correl_MC_denorm[i] = spectra*spectra;
   }

   delete [] cosK_phi;
   delete [] sinK_phi;

   return;
}

void HBT::Cal_azimuthal_dependent_correlationfunction_3D(int iKT, double K_y)
{
   if(fabs(K_y) > 1e-16)
   {
       cout<<"not support for y not equals 0 yet!" << endl;
       return;
   }
   
   double K_T = KT_array[iKT];
   cout << "generating correlation function in 3D for K_T = " 
        << K_T << " GeV ... " << endl;

   double mass = particle_ptr[particle_id].mass;
   
   double *cosK_phi = new double [n_Kphi];
   double *sinK_phi = new double [n_Kphi];
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
      cosK_phi[iphi] = cos(Kphi[iphi]);
      sinK_phi[iphi] = sin(Kphi[iphi]);
   }
   
   double *spectra = new double [n_Kphi];
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
      double ss = 0.0;
      for(int k = 0; k < Emissionfunction_length; k++)
         ss += emission_S_K[k].data[iKT][iphi];
      spectra[iphi] = ss*2;
   }

   for(int i = 0; i < qnpts; i++)  // q_out loop
   {
      double local_q_out = q_out[i];
      for(int j = 0; j < qnpts; j++)  // q_side loop
      {
         double local_q_side = q_side[j];
         for(int k = 0; k < qnpts; k++)  // q_long loop
         {
            double local_q_long = q_long[k];
            cout << "q_out = " << local_q_out << " GeV, "
                 << "q_side = " << local_q_side << " GeV, "
                 << "q_long = " << local_q_long << " GeV... " << endl;

     	      double xsi = K_T*K_T + mass*mass + (local_q_out*local_q_out + local_q_side*local_q_side + local_q_long*local_q_long)/4.0;  //Set Xsi
            double E1sq = xsi + K_T*local_q_out;
            double E2sq = xsi - K_T*local_q_out;
            double qt = sqrt(E1sq) - sqrt(E2sq);
            double qz = local_q_long;

            for(int iphi = 0; iphi < n_Kphi; iphi++)
            {
               double integ1 = 0.0;                         
               double integ2 = 0.0;

               double qx = local_q_out*cosK_phi[iphi] - local_q_side*sinK_phi[iphi];
               double qy = local_q_side*cosK_phi[iphi] + local_q_out*sinK_phi[iphi];
               
               for(int m = 0; m < Emissionfunction_length; m++)
               {
                  double ss = emission_S_K[m].data[iKT][iphi];
                  double tpt = emission_S_K[m].t;
                  double xpt = emission_S_K[m].x;
                  double ypt = emission_S_K[m].y;
                  double zpt = emission_S_K[m].z;

                  for(int ii=0; ii<2; ii++)
                  {
                     zpt = zpt*(-1);
                     double arg = (tpt*qt - (qx*xpt + qy*ypt + qz*zpt))/hbarC;
                     integ1 += cos(arg)*ss;
                     integ2 += sin(arg)*ss;
                  }
               }
               Correl_3D_phidiff_num[iphi][i][j][k] = integ1*integ1+integ2*integ2;
               Correl_3D_phidiff_denorm[iphi][i][j][k] = spectra[iphi]*spectra[iphi];
            }
         }
      }
   }

   delete [] cosK_phi;
   delete [] sinK_phi;
   delete [] spectra;

   return;
}

void HBT::Cal_azimuthal_dependent_KT_inte_correlationfunction_3D(double K_y)
{
   if(fabs(K_y) > 1e-16)
   {
       cout<<"not support for y not equals 0 yet!" << endl;
       return;
   }
   
   cout << "generating correlation function in 3D ... " << endl;

   double mass = particle_ptr[particle_id].mass;
   
   double *cosK_phi = new double [n_Kphi];
   double *sinK_phi = new double [n_Kphi];
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
      cosK_phi[iphi] = cos(Kphi[iphi]);
      sinK_phi[iphi] = sin(Kphi[iphi]);
   }
   
   double *spectra = new double [n_Kphi];
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
      double ss = 0.0;
      for(int iKT = 0; iKT < n_KT; iKT++)
      {
          double KT_local = KT_array[iKT];
          double KT_weight_local = KT_weight[iKT];
          for(int k = 0; k < Emissionfunction_length; k++)
              ss += emission_S_K[k].data[iKT][iphi]*KT_local*KT_weight_local;
      }
      spectra[iphi] = ss*2;
   }

   for(int i = 0; i < qnpts; i++)  // q_out loop
   {
      double local_q_out = q_out[i];
      for(int j = 0; j < qnpts; j++)  // q_side loop
      {
         double local_q_side = q_side[j];
         for(int k = 0; k < qnpts; k++)  // q_long loop
         {
            double local_q_long = q_long[k];
            cout << "q_out = " << local_q_out << " GeV, "
                 << "q_side = " << local_q_side << " GeV, "
                 << "q_long = " << local_q_long << " GeV... " << endl;

            for(int iphi = 0; iphi < n_Kphi; iphi++)
            {
               double integ1 = 0.0;                         
               double integ2 = 0.0;

               double qx = local_q_out*cosK_phi[iphi] - local_q_side*sinK_phi[iphi];
               double qy = local_q_side*cosK_phi[iphi] + local_q_out*sinK_phi[iphi];
     	      
            
               for(int iKT = 0; iKT < n_KT; iKT++)
               {
                   double K_T = KT_array[iKT];
                   double KT_weight_local = KT_weight[iKT];
                   double xsi = K_T*K_T + mass*mass + (local_q_out*local_q_out + local_q_side*local_q_side + local_q_long*local_q_long)/4.0;  //Set Xsi
                   double E1sq = xsi + K_T*local_q_out;
                   double E2sq = xsi - K_T*local_q_out;
                   double qt = sqrt(E1sq) - sqrt(E2sq);
                   double qz = local_q_long;

               
                   for(int m = 0; m < Emissionfunction_length; m++)
                   {
                      double ss = emission_S_K[m].data[iKT][iphi]*K_T*KT_weight_local;
                      double tpt = emission_S_K[m].t;
                      double xpt = emission_S_K[m].x;
                      double ypt = emission_S_K[m].y;
                      double zpt = emission_S_K[m].z;

                      for(int ii=0; ii<2; ii++)
                      {
                         zpt = zpt*(-1);
                         double arg = (tpt*qt - (qx*xpt + qy*ypt + qz*zpt))/hbarC;
                         integ1 += cos(arg)*ss;
                         integ2 += sin(arg)*ss;
                      }
                   }
               }
               Correl_3D_phidiff_num[iphi][i][j][k] = integ1*integ1+integ2*integ2;
               Correl_3D_phidiff_denorm[iphi][i][j][k] = spectra[iphi]*spectra[iphi];
            }
         }
      }
   }

   delete [] cosK_phi;
   delete [] sinK_phi;
   delete [] spectra;

   return;
}

void HBT::Cal_azimuthal_dependent_correlationfunction_MC(int iKT, double K_y)
{
   if(fabs(K_y) > 1e-16)
   {
       cout<<"not support for y not equals 0 yet!" << endl;
       return;
   }
   
   double K_T = KT_array[iKT];
   cout << "generating correlation function in 3D for K_T = " 
        << K_T << " GeV ... " << endl;

   double mass = particle_ptr[particle_id].mass;
   
   double *cosK_phi = new double [n_Kphi];
   double *sinK_phi = new double [n_Kphi];
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
      cosK_phi[iphi] = cos(Kphi[iphi]);
      sinK_phi[iphi] = sin(Kphi[iphi]);
   }
   
   double *spectra = new double [n_Kphi];
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
      double ss = 0.0;
      for(int k = 0; k < Emissionfunction_length; k++)
         ss += emission_S_K[k].data[iKT][iphi];
      spectra[iphi] = ss*2;
   }

   double local_q_out, local_q_side, local_q_long;
   for(int i = 0; i < MC_samples; i++)
   {
       if(i < MC_samples*2/3)
       {
           local_q_out = drand48()*q_max;
           local_q_side = drand48()*q_max;
           local_q_long = 0.0;
       }
       else
       {
           local_q_out = 0.0;
           local_q_side = 0.0;
           local_q_long = drand48()*q_max;
       }

       q_out_MC[i] = local_q_out;
       q_side_MC[i] = local_q_side;
       q_long_MC[i] = local_q_long;
       cout << "q_out = " << local_q_out << " GeV, "
            << "q_side = " << local_q_side << " GeV, "
            << "q_long = " << local_q_long << " GeV... " << endl;

     	 double xsi = K_T*K_T + mass*mass + (local_q_out*local_q_out + local_q_side*local_q_side + local_q_long*local_q_long)/4.0;  //Set Xsi
       double E1sq = xsi + K_T*local_q_out;
       double E2sq = xsi - K_T*local_q_out;
       double qt = sqrt(E1sq) - sqrt(E2sq);
       double qz = local_q_long;

       for(int iphi = 0; iphi < n_Kphi; iphi++)
       {
          double integ1 = 0.0;                         
          double integ2 = 0.0;

          double qx = local_q_out*cosK_phi[iphi] - local_q_side*sinK_phi[iphi];
          double qy = local_q_side*cosK_phi[iphi] + local_q_out*sinK_phi[iphi];
          
          for(int m = 0; m < Emissionfunction_length; m++)
          {
             double ss = emission_S_K[m].data[iKT][iphi];
             double tpt = emission_S_K[m].t;
             double xpt = emission_S_K[m].x;
             double ypt = emission_S_K[m].y;
             double zpt = emission_S_K[m].z;

             for(int ii=0; ii<2; ii++)
             {
                zpt = zpt*(-1);
                double arg = (tpt*qt - (qx*xpt + qy*ypt + qz*zpt))/hbarC;
                integ1 += cos(arg)*ss;
                integ2 += sin(arg)*ss;
             }
          }
          Correl_MC_phidiff_num[iphi][i] = integ1*integ1+integ2*integ2;
          Correl_MC_phidiff_denorm[iphi][i] = spectra[iphi]*spectra[iphi];
       }
   }

   delete [] cosK_phi;
   delete [] sinK_phi;
   delete [] spectra;

   return;
}

void HBT::Cal_azimuthal_dependent_KT_inte_correlationfunction_MC(double K_y)
{
   if(fabs(K_y) > 1e-16)
   {
       cout<<"not support for y not equals 0 yet!" << endl;
       return;
   }
   
   cout << "generating correlation function in MC ... " << endl;

   double mass = particle_ptr[particle_id].mass;
   
   double *cosK_phi = new double [n_Kphi];
   double *sinK_phi = new double [n_Kphi];
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
      cosK_phi[iphi] = cos(Kphi[iphi]);
      sinK_phi[iphi] = sin(Kphi[iphi]);
   }
   
   double *spectra = new double [n_Kphi];
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
      double ss = 0.0;
      for(int iKT = 0; iKT < n_KT; iKT++)
      {
          double KT_local = KT_array[iKT];
          double KT_weight_local = KT_weight[iKT];
          for(int k = 0; k < Emissionfunction_length; k++)
              ss += emission_S_K[k].data[iKT][iphi]*KT_local*KT_weight_local;
      }
      spectra[iphi] = ss*2;
   }

   double local_q_out, local_q_side, local_q_long;
   for(int i = 0; i < MC_samples; i++)
   {
       if(i < MC_samples*2/3)
       {
           local_q_out = drand48()*q_max;
           local_q_side = drand48()*q_max;
           local_q_long = 0.0;
       }
       else
       {
           local_q_out = 0.0;
           local_q_side = 0.0;
           local_q_long = drand48()*q_max;
       }

       q_out_MC[i] = local_q_out;
       q_side_MC[i] = local_q_side;
       q_long_MC[i] = local_q_long;
       cout << "q_out = " << local_q_out << " GeV, "
            << "q_side = " << local_q_side << " GeV, "
            << "q_long = " << local_q_long << " GeV... " << endl;

       for(int iphi = 0; iphi < n_Kphi; iphi++)
       {
          double integ1 = 0.0;                         
          double integ2 = 0.0;

          double qx = local_q_out*cosK_phi[iphi] - local_q_side*sinK_phi[iphi];
          double qy = local_q_side*cosK_phi[iphi] + local_q_out*sinK_phi[iphi];

          for(int iKT = 0; iKT < n_KT; iKT++)
          {
              double K_T = KT_array[iKT];
              double KT_weight_local = KT_weight[iKT];

     	        double xsi = K_T*K_T + mass*mass + (local_q_out*local_q_out + local_q_side*local_q_side + local_q_long*local_q_long)/4.0;  //Set Xsi
              double E1sq = xsi + K_T*local_q_out;
              double E2sq = xsi - K_T*local_q_out;
              double qt = sqrt(E1sq) - sqrt(E2sq);
              double qz = local_q_long;

          
              for(int m = 0; m < Emissionfunction_length; m++)
              {
                 double ss  = emission_S_K[m].data[iKT][iphi]*K_T*KT_weight_local;
                 double tpt = emission_S_K[m].t;
                 double xpt = emission_S_K[m].x;
                 double ypt = emission_S_K[m].y;
                 double zpt = emission_S_K[m].z;

                 for(int ii=0; ii<2; ii++)
                 {
                    zpt = zpt*(-1);
                    double arg = (tpt*qt - (qx*xpt + qy*ypt + qz*zpt))/hbarC;
                    integ1 += cos(arg)*ss;
                    integ2 += sin(arg)*ss;
                 }
              }
          }
          Correl_MC_phidiff_num[iphi][i] = integ1*integ1+integ2*integ2;
          Correl_MC_phidiff_denorm[iphi][i] = spectra[iphi]*spectra[iphi];
       }
   }

   delete [] cosK_phi;
   delete [] sinK_phi;
   delete [] spectra;

   return;
}

void HBT::Output_Correlationfunction_1D(int iKT)
{
   double error = 1e-4;   // "fake" error
   double local_q_out, local_q_side, local_q_long;

   ostringstream oCorrelfun_1D_stream;
   if(iKT < 0)
   {
       double KT_min = paraRdr->getVal("KT_min");
       double KT_max = paraRdr->getVal("KT_max");
       oCorrelfun_1D_stream << path << "/correlfunct" << "_" << particle_ptr[particle_id].name << "_kt_" << KT_min << "_" << KT_max << ".dat";
   }
   else
       oCorrelfun_1D_stream << path << "/correlfunct" << "_" << particle_ptr[particle_id].name << "_kt_" << KT_array[iKT] << ".dat";
   ofstream oCorrelfun_1D;
   oCorrelfun_1D.open(oCorrelfun_1D_stream.str().c_str());
   for(int j = 0; j < ndir; j++)
   {
       for(int k = 0; k < qnpts; k++)
       {
           switch (j)
           {
              case 0:
              {
                 local_q_out  = q_out[k];
                 local_q_side = 0.0e0;
                 local_q_long = 0.0e0;
                 break;
              }
              case 1:
              {
                 local_q_out  = 0.0e0;
                 local_q_side = q_side[k];
                 local_q_long = 0.0e0;
                 break;
              }
              case 2:
              {
                 local_q_out  = 0.0e0;
                 local_q_side = 0.0e0;
                 local_q_long = q_long[k];
                 break;
              }
              case 3:
              {
                 local_q_out  = q_out[k];
                 local_q_side = q_side[k];
                 local_q_long = 0.0e0;
                 break;
              }
              case 4:
              {
                 local_q_out  = q_out[k];
                 local_q_side = -q_side[k];
                 local_q_long = 0.0e0;
                 break;
              }
              default:
              {
                 cout << "error in assigning q values! "<< endl;
                 break;
              }
           }
           oCorrelfun_1D << scientific << setprecision(7) << setw(15)
                         << local_q_out << "  " << local_q_side << "  " 
                         << local_q_long << "  " 
                         << Correl_1D_num[j][k] << "  " 
                         << Correl_1D_denorm[j][k] << "  " 
                         << Correl_1D_num[j][k]/Correl_1D_denorm[j][k] 
                         << "  " << error << endl;
       }
   }
   oCorrelfun_1D.close();
   return;
}

void HBT::Output_Correlationfunction_azimuthal_dependent_1D(int iKT)
{
   double error = 1e-4;   // "fake" error
   double local_q_out, local_q_side, local_q_long;
   for(int i = 0; i < n_Kphi; i++)
   {
       ostringstream oCorrelfun_1D_stream;
       if(iKT < 0)
       {
           double KT_min = paraRdr->getVal("KT_min");
           double KT_max = paraRdr->getVal("KT_max");
           oCorrelfun_1D_stream << path << "/correlfunct" << "_" << particle_ptr[particle_id].name << "_kt_" << KT_min << "_" << KT_max << "_Kphi_" << Kphi[i]-Psi_ev << ".dat";
       }
       else
           oCorrelfun_1D_stream << path << "/correlfunct" << "_" << particle_ptr[particle_id].name << "_kt_" << KT_array[iKT] << "_Kphi_" << Kphi[i]-Psi_ev << ".dat";

       ofstream oCorrelfun_1D;
       oCorrelfun_1D.open(oCorrelfun_1D_stream.str().c_str());
       for(int j = 0; j < ndir; j++)
       {
           for(int k = 0; k < qnpts; k++)
           {
               switch (j)
               {
                  case 0:
                  {
                     local_q_out  = q_out[k];
                     local_q_side = 0.0e0;
                     local_q_long = 0.0e0;
                     break;
                  }
   	            case 1:
                  {
                     local_q_out  = 0.0e0;
                     local_q_side = q_side[k];
                     local_q_long = 0.0e0;
                     break;
                  }
                  case 2:
                  {
                     local_q_out  = 0.0e0;
                     local_q_side = 0.0e0;
                     local_q_long = q_long[k];
                     break;
                  }
                  case 3:
                  {
                     local_q_out  = q_out[k];
                     local_q_side = q_side[k];
                     local_q_long = 0.0e0;
                     break;
                  }
                  case 4:
                  {
                     local_q_out  = q_out[k];
                     local_q_side = -q_side[k];
                     local_q_long = 0.0e0;
                     break;
                  }
                  default:
                  {
                     cout << "error in assigning q values! "<< endl;
                     break;
                  }
               }
               oCorrelfun_1D << scientific << setprecision(7) << setw(15)
                             << local_q_out << "  " << local_q_side << "  " 
                             << local_q_long << "  " 
                             << Correl_1D_phidiff_num[i][j][k] << "  " 
                             << Correl_1D_phidiff_denorm[i][j][k] << "  " 
                             << Correl_1D_phidiff_num[i][j][k]/Correl_1D_phidiff_denorm[i][j][k] 
                             << "  " << error << endl;
            }
       }
       oCorrelfun_1D.close();
   }
   return;
}


void HBT::Output_Correlationfunction_3D(int iKT)
{
   double error = 1e-4;
   ostringstream oCorrelfun_3D_stream;

   if(iKT < 0)
   {
       double KT_min = paraRdr->getVal("KT_min");
       double KT_max = paraRdr->getVal("KT_max");
       oCorrelfun_3D_stream << path << "/correlfunct" << "_" << particle_ptr[particle_id].name << "_kt_" << KT_min << "_" << KT_max << ".dat";
   }
   else
       oCorrelfun_3D_stream << path << "/correlfunct" << "_" << particle_ptr[particle_id].name << "_kt_" << KT_array[iKT] << ".dat";

   ofstream oCorrelfun_3D;
   oCorrelfun_3D.open(oCorrelfun_3D_stream.str().c_str());
   for(int i=0; i < qnpts; i++)
      for(int j=0; j < qnpts; j++)
         for(int k=0; k < qnpts; k++)
              oCorrelfun_3D << scientific << setprecision(7) << setw(15)
                            << q_out[i] << "  " << q_side[j] << "  " 
                            << q_long[k] << "  " << Correl_3D_num[i][j][k] << "  " 
                            << Correl_3D_denorm[i][j][k] << "  "
                            << Correl_3D_num[i][j][k]/Correl_3D_denorm[i][j][k] 
                            << "  " << error << endl;
   oCorrelfun_3D.close();
   return;
}

void HBT::Output_Correlationfunction_azimuthal_dependent_3D(int iKT)
{
   double error = 1e-4;   // "fake" error
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
       ostringstream oCorrelfun_3D_stream;
       if(iKT < 0)
       {
           double KT_min = paraRdr->getVal("KT_min");
           double KT_max = paraRdr->getVal("KT_max");
           oCorrelfun_3D_stream << path << "/correlfunct" << "_" << particle_ptr[particle_id].name << "_kt_" << KT_min << "_" << KT_max << "_Kphi_" << Kphi[iphi]-Psi_ev << ".dat";
       }
       else
           oCorrelfun_3D_stream << path << "/correlfunct" << "_" << particle_ptr[particle_id].name << "_kt_" << KT_array[iKT] << "_Kphi_" << Kphi[iphi]-Psi_ev << ".dat";

       ofstream oCorrelfun_3D;
       oCorrelfun_3D.open(oCorrelfun_3D_stream.str().c_str());
       for(int i=0; i < qnpts; i++)
          for(int j=0; j < qnpts; j++)
             for(int k=0; k < qnpts; k++)
                  oCorrelfun_3D << scientific << setprecision(7) << setw(15)
                                << q_out[i] << "  " << q_side[j] << "  " 
                                << q_long[k] << "  " 
                                << Correl_3D_phidiff_num[iphi][i][j][k] << "  " 
                                << Correl_3D_phidiff_denorm[iphi][i][j][k] << "  "
                                << Correl_3D_phidiff_num[iphi][i][j][k]/Correl_3D_phidiff_denorm[iphi][i][j][k] 
                                << "  " << error << endl;
       oCorrelfun_3D.close();
   }
   return;
}

void HBT::Output_Correlationfunction_MC(int iKT)
{
   double error = 1e-4;
   ostringstream oCorrelfun_MC_stream;
   if(iKT < 0)
   {
       double KT_min = paraRdr->getVal("KT_min");
       double KT_max = paraRdr->getVal("KT_max");
       oCorrelfun_MC_stream << path << "/correlfunct" << "_" << particle_ptr[particle_id].name << "_kt_" << KT_min << "_" << KT_max << ".dat";
   }
   else
       oCorrelfun_MC_stream << path << "/correlfunct" << "_" << particle_ptr[particle_id].name << "_kt_" << KT_array[iKT] << ".dat";

   ofstream oCorrelfun_MC;
   oCorrelfun_MC.open(oCorrelfun_MC_stream.str().c_str());
   for(int i = 0; i < MC_samples; i++)
       oCorrelfun_MC << scientific << setprecision(7) << setw(15)
                     << q_out_MC[i] << "  " << q_side_MC[i] << "  " 
                     << q_long_MC[i] << "  " << Correl_MC_num[i] << "  " 
                     << Correl_MC_denorm[i] << "  "
                     << Correl_MC_num[i]/Correl_MC_denorm[i]
                     << "  " << error << endl;
   oCorrelfun_MC.close();
   return;
}

void HBT::Output_Correlationfunction_azimuthal_dependent_MC(int iKT)
{
   double error = 1e-4;   // "fake" error
   for(int iphi = 0; iphi < n_Kphi; iphi++)
   {
       ostringstream oCorrelfun_MC_stream;
       if(iKT < 0)
       {
           double KT_min = paraRdr->getVal("KT_min");
           double KT_max = paraRdr->getVal("KT_max");
           oCorrelfun_MC_stream << path << "/correlfunct" << "_" << particle_ptr[particle_id].name << "_kt_" << KT_min << "_" << KT_max << "_Kphi_" << Kphi[iphi]-Psi_ev << ".dat";
       }
       else
           oCorrelfun_MC_stream << path << "/correlfunct" << "_" << particle_ptr[particle_id].name << "_kt_" << KT_array[iKT] << "_Kphi_" << Kphi[iphi]-Psi_ev << ".dat";

       ofstream oCorrelfun_MC;
       oCorrelfun_MC.open(oCorrelfun_MC_stream.str().c_str());
       for(int i = 0; i < MC_samples; i++)
           oCorrelfun_MC << scientific << setprecision(7) << setw(15)
                         << q_out_MC[i] << "  " << q_side_MC[i] << "  " 
                         << q_long_MC[i] << "  " 
                         << Correl_MC_phidiff_num[iphi][i] << "  " 
                         << Correl_MC_phidiff_denorm[iphi][i] << "  "
                         << Correl_MC_phidiff_num[iphi][i]/Correl_MC_phidiff_denorm[iphi][i]
                         << "  " << error << endl;
       oCorrelfun_MC.close();
   }
   return;
}

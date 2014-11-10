#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>

#include "HBT.h"

using namespace std;

HBT::HBT(string path_in, ParameterReader* paraRdr_in, particle_info* particle_in, int particle_idx, FO_surf* FOsurf_ptr_in, int FOarray_length)
{
   path = path_in;
   paraRdr =  paraRdr_in;
   particle_ptr = particle_in;
   particle_id = particle_idx;

   FOsurf_ptr = FOsurf_ptr_in;
   FO_length = FOarray_length;

   // initialize eta_s array
   eta_s_npts = paraRdr->getVal("eta_s_npts");
   double eta_s_f = paraRdr->getVal("eta_s_f");
   eta_s = new double [eta_s_npts];
   eta_s_weight = new double [eta_s_npts];
   gauss_quadrature(eta_s_npts, 1, 0.0, 0.0, 0.0, eta_s_f, eta_s, eta_s_weight);

   azimuthal_flag = paraRdr->getVal("azimuthal_flag");

   // initialize Kphi array
   n_Kphi = paraRdr->getVal("n_Kphi");
   Kphi = new double [n_Kphi];
   Kphi_weight = new double [n_Kphi];
   gauss_quadrature(n_Kphi, 1, 0.0, 0.0, 0.0, 2*M_PI, Kphi, Kphi_weight);

   // initialize emission function
   Emissionfunction_length = FO_length*eta_s_npts;
   emission_S_K = new Emissionfunction_data [Emissionfunction_length];
   for(int i=0; i<Emissionfunction_length; i++)
   {
      emission_S_K[i].t = 0.0;
      emission_S_K[i].x = 0.0;
      emission_S_K[i].y = 0.0;
      emission_S_K[i].z = 0.0;
      emission_S_K[i].data = new double [n_Kphi];
      for(int j = 0; j < n_Kphi; j++)
          emission_S_K[i].data[j] = 0.0;
   }

   INCLUDE_SHEAR_DELTAF = paraRdr->getVal("turn_on_shear");
   INCLUDE_BULK_DELTAF = paraRdr->getVal("turn_on_bulk");
   flag_neg = paraRdr->getVal("flag_neg");

   // initialize correlation function
   qnpts = paraRdr->getVal("qnpts");
   MCint_calls = paraRdr->getVal("MCint_calls");

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

   Correl_1D_out_num = new double [qnpts];
   Correl_1D_side_num = new double [qnpts];
   Correl_1D_long_num = new double [qnpts];
   Correl_1D_out_denorm = new double [qnpts];
   Correl_1D_side_denorm = new double [qnpts];
   Correl_1D_long_denorm = new double [qnpts];
   for(int i=0; i<qnpts; i++)
   {
      Correl_1D_out_num[i] = 0.0;
      Correl_1D_side_num[i] = 0.0;
      Correl_1D_long_num[i] = 0.0;
      Correl_1D_out_denorm[i] = 0.0;
      Correl_1D_side_denorm[i] = 0.0;
      Correl_1D_long_denorm[i] = 0.0;
   }
   
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

   return;
}

HBT::~HBT()
{

   delete [] Kphi;
   delete [] Kphi_weight;

   for(int i = 0; i < Emissionfunction_length; i++)
       delete [] emission_S_K[i].data;
   delete [] emission_S_K;

   delete [] q_out;
   delete [] q_side;
   delete [] q_long;

   delete [] Correl_1D_out_num;
   delete [] Correl_1D_side_num;
   delete [] Correl_1D_long_num;
   delete [] Correl_1D_out_denorm;
   delete [] Correl_1D_side_denorm;
   delete [] Correl_1D_long_denorm;

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

   return;
}

void HBT::calculate_azimuthal_dependent_HBT_radii(double p_T, double y)
{
   cout << "Calculating "<< particle_ptr[particle_id].name << endl;

   SetEmissionData(FOsurf_ptr, y, p_T);
   Cal_HBTRadii_fromEmissionfunction(p_T, y);
   
   // under construction...
}

void HBT::calculate_azimuthal_averaged_HBT_radii(double y)
{
   cout << "Calculating "<< particle_ptr[particle_id].name << endl;

   double KT_min = paraRdr->getVal("KT_min");
   double KT_max = paraRdr->getVal("KT_max");
   double n_KT = paraRdr->getVal("n_KT");
   double dKT = (KT_max - KT_min)/(n_KT - 1);

   for(int i = 0; i < n_KT; i++)
   {
       double KT_local = KT_min + i*dKT;
       SetEmissionData(FOsurf_ptr, y, KT_local);
       Cal_azimuthal_averaged_correlationfunction_1D(KT_local, y);
       Output_Correlationfunction_1D(KT_local);
       //Cal_azimuthal_averaged_correlationfunction_3D(KT_local, y);
       //Output_Correlationfunction_1D(KT_local);
   }
}

void HBT::calculate_azimuthal_averaged_KT_integrated_HBT_radii(double y)
{
   cout << "Calculating "<< particle_ptr[particle_id].name << endl;
   
   // under construction...
}

void HBT::SetEmissionData(FO_surf* FO_surface, double K_rap, double K_T)
// compute emission function at a given pair momentum
{
  double tol = 1e-15;
  double mass = particle_ptr[particle_id].mass;
  double mT = sqrt(mass*mass + K_T*K_T);

  double *K_x = new double [n_Kphi];
  double *K_y = new double [n_Kphi];
  for(int i = 0; i < n_Kphi; i++)
  {
      K_x[i] = K_T*cos(Kphi[i]);
      K_y[i] = K_T*sin(Kphi[i]);
  }

  int idx = 0;
  for(int i=0; i<eta_s_npts; i++)
  {
      double local_eta_s = eta_s[i];
      double ch_localetas = cosh(local_eta_s);
      double sh_localetas = sinh(local_eta_s);

      double K_0 = mT*cosh(K_rap - local_eta_s);
      double K_z = mT*sinh(K_rap - local_eta_s);
      
      for (int j = 0; j < FO_length; j++)
	{
          for(int iphi = 0; iphi < n_Kphi; iphi++)
          {
              double K_x_local = K_x[iphi];
              double K_y_local = K_y[iphi];
              double S_p = Emissionfunction(K_0, K_x_local, K_y_local, K_z, &FO_surface[j]);
              if (flag_neg == 1 && S_p < tol)
              {
                 S_p = 0.0e0;
              }
	        else
              {
                 double S_p_withweight = S_p*FO_surface[j].tau*eta_s_weight[i];
                 emission_S_K[idx].data[iphi] = S_p_withweight;
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

void HBT::Cal_HBTRadii_fromEmissionfunction(double K_T, double K_y)
{
  for(int iphi = 0; iphi < n_Kphi; iphi++)
  {
      double K_phi_local = Kphi[iphi];
      double* resultsX = new double[15];
      for(int i = 0; i < 15; i++)
         resultsX[i] = 0.0e0;

      for(int i=0; i < Emissionfunction_length; i++)
      {
         double S_p = emission_S_K[i].data[iphi];
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
      double m_T = sqrt(mass*mass + K_T*K_T);
      double beta_T = K_T/m_T;
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
  return;
}                                         	

void HBT::Cal_azimuthal_averaged_correlationfunction_1D(double K_T, double K_y)
{
   if(fabs(K_y) > 1e-16)
   {
       cout<<"HBT:: not support for y is not equal to 0 yet!" << endl;
       return;
   }
   
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
         double ss  = emission_S_K[k].data[iphi]*Kphi_weight[iphi];
         spectra += ss*2;
      }
   }

   for(int i = 0; i < qnpts; i++)
   {
      cout << "calculating q_mu = " << q_out[i] << " GeV..." << endl;
      double values_num[3];
      for (int ops = 0; ops < 3; ops++)
         values_num[ops] = 0.0e0;
      for (int l = 0; l < 3; l++)
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

         for(int iphi = 0; iphi < n_Kphi; iphi++)
         {
             double qx = local_q_out*cosK_phi[iphi] - local_q_side*sinK_phi[iphi];
             double qy = local_q_side*cosK_phi[iphi] + local_q_out*sinK_phi[iphi];

             for(int k = 0; k < Emissionfunction_length; k++)
             {
                double ss  = emission_S_K[k].data[iphi]*Kphi_weight[iphi];
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
         values_num[l] = localvalue;
      }
      Correl_1D_out_num[i]  = values_num[0];
      Correl_1D_side_num[i] = values_num[1];
      Correl_1D_long_num[i] = values_num[2];
      Correl_1D_out_denorm[i]  = spectra*spectra;
      Correl_1D_side_denorm[i] = spectra*spectra;
      Correl_1D_long_denorm[i] = spectra*spectra;
   }

   delete [] cosK_phi;
   delete [] sinK_phi;

   return;
} 

void HBT::Cal_azimuthal_averaged_correlationfunction_3D(double K_T, double K_y)
{
   if(fabs(K_y) > 1e-16)
   {
       cout<<"not support for y not equals 0 yet!" << endl;
       return;
   }
   
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
         double ss  = emission_S_K[k].data[iphi]*Kphi_weight[iphi];
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
                  double ss = emission_S_K[m].data[iphi]*Kphi_weight[iphi];
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

int HBT::binary_search(double* dataset, int data_length, double value)
{
  int lowbin = 0;
  int midbin = 0;
  int highbin = data_length;
  int stop = 0;
  int dbin = highbin - lowbin;
  if (dbin == 0) 
  {
     cout << "binary search::You screwed up in the table lookup dummy. Fix it." << endl;
     exit(1);
  }

  while(stop == 0) 
  {
      dbin = highbin - lowbin;
      midbin = (int)(lowbin + dbin/2);
      if(dbin == 1)
      {
	  stop = 1;
      }
      else if(value > dataset[midbin])
      {
	  lowbin = midbin;
      }
      else
      {
	  highbin = midbin;
      }
  }
  return(lowbin);
}


void HBT::Output_Correlationfunction_1D(double K_T)
{
   double error = 1e-4;   // "fake" error
   ostringstream oCorrelfun_1D_stream;
   oCorrelfun_1D_stream << path << "/correlfunct" << "_" << particle_ptr[particle_id].name << "_kt_" << K_T << ".dat";
   ofstream oCorrelfun_1D;
   oCorrelfun_1D.open(oCorrelfun_1D_stream.str().c_str());
   // out direction
   for(int i=0; i < qnpts; i++)
       oCorrelfun_1D << scientific << setprecision(7) << setw(15)
                     << q_out[i] << "  " << 0.0e0 << "  " 
                     << 0.0e0 << "  " << Correl_1D_out_num[i] << "  " 
                     << Correl_1D_out_denorm[i] << "  " 
                     << Correl_1D_out_num[i]/Correl_1D_out_denorm[i] << "  "
                     << error << endl;
   // side direction
   for(int i=0; i < qnpts; i++)
       oCorrelfun_1D << scientific << setprecision(7) << setw(15)
                     << 0.0e0 << "  " << q_side[i] << "  " 
                     << 0.0e0 << "  " << Correl_1D_side_num[i] << "  " 
                     << Correl_1D_side_denorm[i] << "  "
                     << Correl_1D_side_num[i]/Correl_1D_side_denorm[i] << "  "
                     << error << endl;
   // long direction
   for(int i=0; i < qnpts; i++)
       oCorrelfun_1D << scientific << setprecision(7) << setw(15)
                     << 0.0e0 << "  " << 0.0e0 << "  " 
                     << q_long[i] << "  " << Correl_1D_long_num[i] << "  "
                     << Correl_1D_long_denorm[i] << "  "
                     << Correl_1D_long_num[i]/Correl_1D_long_denorm[i] << "  "
                     << error << endl;
   return;
}

void HBT::Output_Correlationfunction_3D(double K_T)
{
   double error = 1e-4;
   ostringstream oCorrelfun_3D_stream;
   oCorrelfun_3D_stream << path << "/correlfunct" << "_" << particle_ptr[particle_id].name << "_kt_" << K_T << ".dat";
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
                            << error << endl;
   return;
}

#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>

#include<gsl/gsl_sf_bessel.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#include "Arsenal.h"
using namespace std;

unsigned long int random_seed()
{
  unsigned long int seed = 0;
  unsigned long int *seed_ptr = &seed;

  ifstream dev_urandom ("/dev/urandom", ios::in | ios::binary);
  
  dev_urandom.read((char *) seed_ptr, sizeof (long int));

  dev_urandom.close();
  return(*seed_ptr);
}


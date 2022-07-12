/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __fluid__
#define __fluid__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
const int kb = 1.;
int n_props, iv, ik, it, ie, iw, ip, ig;
double vtail, ptail, bin_size, sd;
const int nbins = 100;
double bins[nbins];
double walker[m_props];

// averages
double blk_av[m_props], bin_av[nbins], blk_norm, accepted, attempted;
double glob_av[m_props], glob_av2[m_props], glob_bin[nbins], glob_bin2[nbins];
double stima_pot, stima_pres, stima_kin, stima_etot, stima_temp, stima_g[nbins];
double err_pot, err_pres, err_kin, err_etot, err_temp, err_gdir, err_g[nbins];

//configuration
const int m_part=108;
double x[m_part],    y[m_part],    z[m_part];
double xold[m_part], yold[m_part], zold[m_part];
double vx[m_part],  vy[m_part],   vz[m_part];

// thermodynamical state
int npart;
double m_beta,temp,energy,vol,rho,box,rcut;

// simulation
int iNVET, nstep, nstep_eq, nblk, restart;
double delta;

//pigreco
const double pi=3.1415927;

//functions
void Input(char *);
void Reset(int);
void Accumulate(void);
void Averages(int, char**);
void Move(void);
void ConfFinal(char**);
void ConfXYZ(int);
void Measure(void);
double Boltzmann(double, double, double, int);
double Pbc(double);
double Error(double,double,int);
double Force(int, int);

void PrePrint(int, double*);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

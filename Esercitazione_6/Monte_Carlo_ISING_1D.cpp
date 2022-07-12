/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <string>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(){
  J = 1;
  int npoint = 10;
  double T [npoint  + 1];

  // Caso h = 0
  h = 0;
  //Misuro le osservabili per 10 valori di temperatura nell'intervallo [0.5, 2.0]
  double dt = 1.5/npoint;
  for(int i = 0; i <= npoint; i++){
    T[i] = 0.5 + i*dt;
  }

  // i = 0 --> Gibbs
  // i = 1 --> Metropolis
  for(int i = 0; i < 2; i++){
    metro = i;
    //Per ogni valore di temperatura, misuro le grandezze richieste, dopo aver equilibrato la simulazione
    for(auto t:T){
      Input(h, t);
      // Equilibrazione di 10000 passi
      nblk = 10000;
      for(int i = 0; i < nblk; i++){
        Move(metro);
      }

      //Eseguo l'algoritmo per 20 blocchi di 1000 step
      nblk = 20;
      nstep = 40000;
      for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
      {
        Reset(iblk);
        for(int istep=1; istep <= nstep; ++istep)
        {
          Move(metro);
          Measure();
          Accumulate();
        }
        Averages(iblk);
      }
      Averages(t, h);
      
    }
  }
  
  cout << endl;

  // Caso h = 0.02
  h = 0.02;
  //Misuro le osservabili per 10 valori di temperatura nell'intervallo [0.5, 2.0]

  //Per ogni valore di temperatura, misuro le grandezze richieste, dopo aver equilibrato la simulazione
  for(int i = 0; i < 2; i++){
    metro = i;
    for(auto t:T){
      Input(h, t);
      // Equilibrazione di 10000 passi
      nblk = 10000;
      for(int i = 0; i < nblk; i++){
        Move(metro);
      }

      //Eseguo l'algoritmo per 20 blocchi di 1000 step
      nblk = 20;
      nstep = 40000;
      for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
      {
        Reset(iblk);
        for(int istep=1; istep <= nstep; ++istep)
        {
          Move(metro);
          Measure();
          Accumulate();
        }
        Averages(iblk);
      }
      Averages(t, h);
      
    }
  }

  return 0;
}


void Input(double m_h, double m_t)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  cout << "Temperature = " << m_t << endl;
  temp = m_t;
  beta = 1.0/temp;

  nspin = 50;
  cout << "Number of spins = " << nspin << endl;

  J = 1;
  cout << "Exchange interaction = " << J << endl;

  h = m_h;
  cout << "External field = " << m_h << endl << endl;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


  //Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}


void Move(int metro)
{
  int o;
  // double p, energy_old, energy_new, sm;
  // double energy_up, energy_down;

  double delta_E;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //FUNZIONAMENTO DEL METROPOLIS
      //  - Ho selezionato una particella a caso
      //  - genero un nuovo stato scegliendo un'opportuna T(new|old). Nel nostro caso inverto lo spin della singola
      //    particella.
      //  - Accetto il cambio di stato appena eseguito con una probabilità A(new|old) = min[1, p(new)/p(old)]. 
      //    Nel caso di un sistema di N particelle, p(stato) è data dal peso di Bolzmann, perciò il rapporto
      //    tra le probabilità è exp[-beta(E_new - E_old)], con E_new - E_old =  2*J*s_k^old \sum s_i^old.
      //    Nel caso in cui p(new)/p(old) < 1 (E_new - E_old > 0), come valutiamo la probabilità A di accettare
      //    la mossa? Scegliamo un numero casuale r tra 0 e 1. Se r < A, allora accettiamo, altrimenti no.
      //    La mossa scartata non è buttata: semplicemente si lascia il sistema com'è.
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // delta_E = 2 * J * s[o] * (s[Pbc(o-1)] + s[Pbc(o+1)]);
      delta_E = Boltzmann(- s[o], o) - Boltzmann(s[o], o);
      double A = fmin(1, exp(-beta*delta_E) );
      double r = rnd.Rannyu();
      if (r < A){
        s[o] = -s[o];
        accepted++;
      } 
      attempted++;
    }
    else //Gibbs sampling
    {
      delta_E = Boltzmann(- 1, o) - Boltzmann(1, o); 
      // Probabilità che s[o] sia posto a +1
      double p = 1./ (1 + exp(-beta*delta_E));
      double r = rnd.Rannyu();
      if (r < p){
        s[o] = +1;
      } else {
        s[o] = -1;
      }
    }
  }
}

void Measure()
{
  // int bin;
  double u = 0.0, m = 0.0;

  //cycle over spins
  for (int i=0; i<nspin; ++i) {
    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
    m += s[i];
  }

  walker[iu] = u;
  walker[im] = m;
  walker[ic] = u*u;
  walker[ix] = m*m;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
    // ofstream Ene, Heat, Mag, Chi;
    // const int wd=14;
    
    cout << "\rBlock number " << iblk << " - Acceptance rate " << accepted/attempted << " " << flush;
    
    // Ene.open("output.ene.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    // Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    // Ene.close();

    // Mag.open("output.mag.0",ios::app);
    stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetizzazione
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    // Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    // Mag.close();

    stima_c = pow(beta,2)*(blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm, 2))/(double)nspin; // Calore specifico
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);

    stima_x = beta*blk_av[ix]/blk_norm/(double)nspin; // Suscettività
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);

    // cout << "----------------------------" << endl << endl;
}

void Averages(double t, double m_h) //Print results for current block
{
    
  ofstream Ene, Heat, Mag, Chi;
  const int wd=14;
    
    // cout << "Block number " << iblk << endl;
    // cout << "Acceptance rate " << accepted/attempted << endl << endl;

  string str_metro = metro?"metro":"gibbs";
  string str_h = to_string(m_h);
    
  Ene.open("output.ene.T.h" + str_h + "." + str_metro,ios::app);
  Ene << setw(wd) << t << setw(wd) << glob_av[iu]/(double)nblk << setw(wd) << err_u << endl;
  Ene.close();

  Mag.open("output.mag.T.h" + str_h + "." + str_metro,ios::app);
  Mag << setw(wd) << t << setw(wd) << glob_av[im]/(double)nblk << setw(wd) << err_m << endl;
  Mag.close();

  Heat.open("output.cal.T.h" + str_h + "." + str_metro,ios::app);
  Heat << setw(wd) << t << setw(wd) << glob_av[ic]/(double)nblk << setw(wd) << err_c << endl;
  Heat.close();

  Chi.open("output.chi.T.h" + str_h + "." + str_metro,ios::app);
  Chi << setw(wd) << t << setw(wd) << glob_av[ix]/(double)nblk << setw(wd) << err_x << endl;
  Chi.close();

    // cout << "----------------------------" << endl << endl;
}



void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

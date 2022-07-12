#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"

#define DELTA 2

int accepted;
int throws;
double DELTA2;
double mu, sigma, muold, sigmaold;

using namespace std;

double p_met(Random & rnd, double x); // Funzione che campiona il modulo quadro di psi (è una densità di probabilità).
double Psi2(double x);
double Psi(double x);
double V(double x);
double Eloc(double x);
double Integrale(Random &, double &); //Funzione che calcola l'integrale per il valor medio dell'eneergia usando integrazione montecarlo

double H_T(Random & rnd, int, int, bool);

double Update_par(Random &rnd, double beta);

void initialize_rnd(Random &);


int main (int argc, char *argv[]){

  // Simulated Annealing.
  //   

  Random rnd;
  initialize_rnd(rnd);

  // double en = 0;
  // double y = 1;
  // int tot = 100000;
  // for (int i = 0; i < tot; i++)
  // {
  //   en += Integrale(rnd, y, 0.8, 0.63);
  // }

  // cout << en/tot << endl;

  // Stima di H
  mu = 0.5;
  sigma = 0.5;
  double H;
  
  int L = 10; // Numero di punti estratti in ogni blocco (a temperatura costamte)
  // int N = 100;
  // int M = N*L;

  vector <double> betas;
  double bmin = 1;
  double dbeta = 2.5;

  ofstream plot; plot.open("3dplot.dat");
  ofstream out; out.open("annealing.dat");

  double prec = 0.0015; // Precisione a cui si arresta l'annealing
  double average = 0, average2 = 0, error = prec, blk_err = prec,  beta;
  accepted = throws =  0;

  for(int k = 0; blk_err >= prec; k++){
    double est_H = 0, est_H2 = 0; // Stima di <H>_T

    betas.push_back(bmin + k*dbeta); // Aggiorno la temperatura
    beta = betas.at(k);
    DELTA2 = 1./sqrt(beta);
  //Ciclo nel blocco a temperatura costante. Ad ogni passo cerco una nuova configurazione di parametri con Metropolis
    for(int j = 0; j < L; j++){
      H = Update_par(rnd, beta);
      est_H += H;
      est_H2 += H*H;

      // if(!(j%5)) cout << "Step " << j << "/" << L << " \r" << flush;
    }
  // Deviazione standard sulla stima di H calcolata a temperatura costante.
    blk_err = sqrt(est_H2/L - est_H*est_H/L/L);
  // Calcolo di errore e di media progressiva
    est_H = est_H/L;

  // Formato: # blocco, media progressiva, errore progressivo, media del blocco, errore del singolo blocco
    out << k + 1 << "," << average / (k + 1) << "," << error << "," << est_H << "," << blk_err << endl;

    cout << endl << "BLOCCO " << k << " - t = " << 1./beta << " - passo = " << DELTA2 << endl;
    cout << " Mu = " << mu << ", sigma = " << sigma << endl;
    cout << " Errore su <H> = " << blk_err << endl;
    cout << " Tasso accettazione metropolis: " << (double)accepted/throws << endl << endl;
  }

  cout << endl << "mu = " << mu << endl << "sigma = " << sigma << endl;
  out.close();

  cout << endl << "Proseguimento della ricerca a passo costante per 10 blocchi da 100 passi..." << endl;
  L = 100;
  int N = 10;
  out.open("annealing_final.dat");
  bmin = betas.at(betas.size()-1);
  betas.resize(0);
  accepted = throws = 0;
  average = average2 = 0;
  double avmu = 0, avsigma = 0;

  vector <double> mus, sigmas;

  for(int k = 0; k < N; k++){
    betas.push_back(bmin + k*dbeta); // Aggiorno la temperatura
    beta = betas.at(k);
    DELTA2 = 1./beta/2;
    double est_H = 0, est_H2 = 0; // Stima di <H>_T
    avmu = avsigma = 0;
  //Ciclo nel blocco a temperatura costante. Ad ogni passo cerco una nuova configurazione di parametri con Metropolis
    for(int j = 0; j < L; j++){
      H = Update_par(rnd, beta);
      avmu += mu;
      avsigma += sigma;
      est_H += H;

      // if(!(j%5)) cout << "Step " << j << "/" << L << " \r" << flush;
    }
  // Calcolo di errore e di media progressiva
    est_H = est_H/L;
    average += est_H;
    average2 += est_H*est_H;
    error = k ? sqrt( ( average2/(k+1) - pow(average/(k+1),2) )/k ) : prec;

  // Stima di mu e sigma
    mus.push_back(avmu/L);
    sigmas.push_back(avsigma/L);

  // Formato: # blocco, media progressiva, errore progressivo, media del blocco, errore del singolo blocco, mu, sigma
    out << k + 1 << "," << average / (k + 1) << "," << error << "," << est_H << "," << blk_err << "," << mus.at(k) << "," << sigmas.at(k) << endl;

    cout << endl << "BLOCCO " << k << "/" << N << " - t = " << 1./beta << " - passo = " << DELTA2 << endl;
    cout << " Mu = " << mus.at(k) << ", sigma = " << sigmas.at(k) << endl;
    cout << " Errore su <H> = " << error << endl;
    cout << " Tasso accettazione metropolis: " << (double)accepted/throws << endl << endl;
  }

  cout << endl << "mu = " << mus.at(mus.size()-1) << endl << "sigma = " << sigmas.at(sigmas.size() - 1) << endl;

// Codice per mostrare la variazione della stima di <H> in funzione dei blocchi
  double Htmp = H_T(rnd, -1, -1, true);

// Codice per l'istogramma col campionamento del modulo quadro
  // vector<double> integral;
  ofstream histo; histo.open("histo.dat");
  mu = mus.at(mus.size()-1);
  sigma = sigmas.at(sigmas.size() - 1);
  int max = 50000;
  double x = 1;
  for(int i = 0; i < max; i++){
    // integral.push_back(Integrale(rnd, x));
    x = p_met(rnd, x);
    histo << p_met(rnd, x) << endl;
  }

  histo.close();

    
    
  cout << endl;
  rnd.SaveSeed();
  return 0;
}

void initialize_rnd(Random &rnd){
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
}

double Integrale(Random & rnd, double &x){
    double result = 0;
    x = p_met(rnd, x);
    result = Eloc(x);
    
    return result;
}

double H_T(Random &rnd, int M, int L, bool print){
  int N;
  ofstream out; 

  if(print){
    out.open("averageH.dat");
    N = 100;
    L = 5000;
  }
  else{
    N = M/L;
  }
  
  int EQ = 100;
  double tmp;
  double x = 1, est_I; // Punto del campionamento con metropolis

  // Equilibrazione
  for(int i = 0; i < EQ; i++){
    tmp = Integrale(rnd, x);
  }

  double average = 0, average2 = 0, error;
  for(int k = 0; k < N; k++){
    est_I = 0;
    //Ciclo nel blocco
    for(int j = 0; j < L; j++){
      est_I += Integrale(rnd, x);
    }
    est_I = est_I / L;
    average += est_I;

    if(print){
      average2 += est_I*est_I;
      error = k ? sqrt( ( average2/(k+1) - pow(average/(k+1),2) )/k ) : 0;
      out << k + 1 << "," << average / (k + 1) << "," << error << endl;
      cout << "Progresso: blocco " << k << " di " << N << " - " << int( 100*(double)k/N ) << "% \r" << flush;
    }
    
  }

  if(print) out.close();

  return average / N;
}

double p_met(Random &rnd, double x){
  // Campionamento della densità di probabilità tramite Metropolis
  // throws++;
  double xnew = x + rnd.Rannyu(-DELTA, DELTA);
  double p = fmin(1., Psi2(xnew)/Psi2(x) ); // Probabilità con cui accettare la mossa del Metropolis
  double thrw = rnd.Rannyu();
  if(thrw < p){
    // accepted++;
    return xnew;
  }
  else return x;
}

double Psi2(double x){
  // Funzione che calcola il quadrato della funzione d'onda di test
  double a1 = (x - mu)/M_SQRT2/sigma;
  double a2 = (x + mu)/M_SQRT2/sigma;
  double psi = exp(-a1*a1) + exp(-a2*a2);
  return psi*psi;
}

double Psi(double x){
  // Funzione che calcola il quadrato della funzione d'onda di test
  double a1 = (x - mu)/M_SQRT2/sigma;
  double a2 = (x + mu)/M_SQRT2/sigma;
  double psi = exp(-a1*a1) + exp(-a2*a2);
  return psi;
}

double Eloc(double x){
  // Funzione che restituisce E_loc, ossia H*Psi/Psi
  double alpha = pow((x+mu)/sigma,2);
  double beta = pow((x-mu)/sigma,2);
  double derivative2 = -0.5*1./sigma/sigma * ( alpha*exp(-alpha/2) + beta*exp(-beta/2) - exp(-alpha/2) - exp(-beta/2) );

  return (derivative2 + V(x)*Psi(x))/Psi(x);
}

double V(double x){
  return pow(x, 4) - 5./2. * x*x;
}

double Update_par(Random &rnd, double beta){
  throws ++;
  muold = mu;
  sigmaold = sigma;
  double H_old = H_T(rnd, 70000, 1, false);

  mu = mu + rnd.Rannyu(-DELTA2, DELTA2);
  sigma = sigma + rnd.Rannyu(-DELTA2, DELTA2);  
  double H_new = H_T(rnd, 70000, 1, false);
  double P = (H_new - H_old > 0) ? exp(-beta * (H_new - H_old)) : 1;

  double thrw = rnd.Rannyu();

  if(thrw < P){
    accepted ++;
    return H_new;
  }
  else {
    mu = muold;
    sigma = sigmaold;
    return H_old;
  }
}
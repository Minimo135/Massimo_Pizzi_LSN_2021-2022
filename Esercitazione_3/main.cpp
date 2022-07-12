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
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <numeric>
#include <cmath>
#include "random.h"

//#define DEBUG

using namespace std;

void initialize_rnd(Random &);
double black_scholes(string option, double S, double t, double T, double K, double r, double sigma);
double N(double d);
double GBM(double tmax, double nstep, double S0, double mu, double sigma, Random & rnd);
 
int main (int argc, char *argv[]){

   Random rnd;
   initialize_rnd(rnd);

   /*
   Esercizio 3:
   Calcolo dell'European call-option price e European put-option price

   1. Campiono direttamente S(T) un numero M di volte.
   2. Divido [0,T] in nstep e campiono il moto browniano geometrico (sempre M volte)
   */

   double S0 = 100.;
   double T = 1.;
   double K = 100.;
   double r = 0.1;
   double sigma = 0.25;

   int M = 50000;
   int N = 100;
   int L = M/N;

   double call_par, put_par, t, S;
   double call_par_100, put_par_100, S100;

   vector<double> call_mean, call_mean_100, call_mean2, call_mean_1002;
   vector<double> put_mean, put_mean_100, put_mean2, put_mean_1002;

   vector<double> call_progressivemean, call_progressivemean100;
   vector<double> put_progressivemean, put_progressivemean100;
   vector<double> err_call, err_call100;
   vector<double> err_put, err_put100;

   ofstream out, out100;
   out.open("ES_03_1.dat");
   out100.open("ES_03_100.dat");

   //Ciclo sui blocchi
   for(int k = 0; k < N; k++){
      //Ciclo nel blocco
      call_par = put_par = call_par_100 = put_par_100 = 0;
      for(int j = 0; j < L; j++){
         t = T;
         S = GBM(t, 1, S0, r, sigma, rnd);
         S100 = GBM(t, 100, S0, r, sigma, rnd);

         call_par += exp(-r*T)*max(0., S - K);
         put_par += exp(-r*T)*max(0., K - S);
         call_par_100 += exp(-r*T)*max(0., S100 - K);
         put_par_100 += exp(-r*T)*max(0., K - S100);
      }
      call_mean.push_back(call_par/L);
      put_mean.push_back(put_par/L);
      call_mean_100.push_back(call_par_100/L);
      put_mean_100.push_back(put_par_100/L);

      call_mean2.push_back( pow(call_par/L,2) );
      call_mean_1002.push_back( pow(call_par_100/L,2) );
      put_mean2.push_back( pow(put_par/L,2) );
      put_mean_1002.push_back( pow(put_par_100/L,2) );

      #ifdef DEBUG
      cout << "call_media = " << call_mean.at(0) << endl;
      #endif

      call_progressivemean.push_back( accumulate(call_mean.begin(), call_mean.end(), 0.)/call_mean.size() );
      put_progressivemean.push_back( accumulate(put_mean.begin(), put_mean.end(), 0.)/put_mean.size() );
      call_progressivemean100.push_back( accumulate(call_mean_100.begin(), call_mean_100.end(), 0.)/call_mean_100.size() );
      put_progressivemean100.push_back( accumulate(put_mean_100.begin(), put_mean_100.end(), 0.)/put_mean_100.size() );

      //Errore...
      err_call.push_back( sqrt( accumulate(call_mean2.begin(), call_mean2.end(), 0.)/call_mean2.size() - pow( call_progressivemean.at(k) ,2) )/sqrt(err_call.size()) );
      err_call100.push_back( sqrt( accumulate(call_mean_1002.begin(), call_mean_1002.end(), 0.)/call_mean_1002.size() - pow( call_progressivemean100.at(k) ,2) )/sqrt(err_call100.size()) );
      err_put.push_back( sqrt( accumulate(put_mean2.begin(), put_mean2.end(), 0.)/put_mean2.size() - pow( put_progressivemean.at(k) ,2) )/sqrt(err_put.size()) );
      err_put100.push_back( sqrt( accumulate(put_mean_1002.begin(), put_mean_1002.end(), 0.)/put_mean_1002.size() - pow( put_progressivemean100.at(k) ,2) )/sqrt(err_put100.size()) );

      out << k << "," << call_progressivemean.at(k) << "," << err_call.at(k) << "," << put_progressivemean.at(k) << "," << err_put.at(k) << endl;
      out100 << k << "," << call_progressivemean100.at(k) << "," << err_call100.at(k) << "," << put_progressivemean100.at(k) << "," << err_put100.at(k) << endl;
   }

   out.close();
   out100.close();

   cout << "Call (single step) = " << call_progressivemean.at(call_progressivemean.size()-1) << endl;
   cout << "Put (single step) = " << put_progressivemean.at(put_progressivemean.size()-1) << endl;
   cout << "Call (100 step) = " << call_progressivemean100.at(call_progressivemean100.size()-1) << endl;
   cout << "Put (100 step) = " << put_progressivemean100.at(put_progressivemean100.size()-1) << endl;

   rnd.SaveSeed();
   return 0;
}

void initialize_rnd(Random &rnd){
   int seed[4];
   int p1, p2;
   int tmp;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> tmp >> tmp  >> tmp >> p1 >> p2 ;
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

double black_scholes(string option, double S, double t, double T, double K, double r, double sigma){
   double d_1 = 1./(sigma*sqrt(T-t)) * ( log(S/K) + (r + sigma*sigma/2.*(T-t) ) ) ;
   double d_2 = d_1 - sigma*sqrt(T-t);
   double output = 0;
   if(option == "call"){
      output = S*N(d_1) - K*exp(-r*(T-t)) * N(d_2);

      #ifdef DEBUG
      cout << "d1 = " << d_1 << endl;
      #endif

      return output;
   }
   else if(option == "put"){
      output = S*( N(d_1) - 1 ) - K*exp(-r*(T-t) ) * (N(d_2) - 1);

      #ifdef DEBUG
      cout << "d2 = " << d_2 << endl;
      #endif

      return output;
   }
   else {
      cerr << "Inserita opzione non valida in black_scholes" << endl;
      return -1;
   }
}

double N(double d){
   return 0.5*(1 + erf(d/M_SQRT2));
}

double GBM(double tmax, double nstep, double S0, double mu, double sigma, Random & rnd){
   double S = S0;
   double step = tmax/nstep;

   if(tmax > 0){
      for(int i = 1; i <= nstep; i++){
         //cout << i*step << endl;
         S = S * exp( (mu - 0.5*sigma*sigma )*( step ) + sigma*rnd.Gauss(0,1)*sqrt( step ) );
         
         #ifdef DEBUG
         cout << "price at t = " << i*step << ", " << S << endl;
         #endif
      }
   }

   #ifdef DEBUG
   cout << "final price = " << S << endl;
   #endif

   return S;
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

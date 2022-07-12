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
#include <cmath>
#include <vector>
#include <numeric>
#include "random.h"

using namespace std;

void rnd_initialize(Random &, int *);
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   rnd_initialize(rnd, seed);
    
   ofstream out;
    
   int M = 10000; //Numero totale di "lanci"
   int N = 100; //Numero di blocchi
      
    //ES 01.2 - 1
    //Esponenziale
    out.open("out4");
    out.seekp(0);
    double N_tot[] = {1,2,10,100};
    vector<double> S_N;
    int dime = 4;
    double partial;
    double lambda = 1;
    for (int i = 0; i < M; i++){
        S_N.resize(0);
        out << i+1;
        //Calcolo 4  versioni diverse per S_N, N = 1,2,10,100
        for(int j = 0; j < dime; j++){
            partial = 0;
            for(int k = 0; k < N_tot[j]; k++){
                partial += rnd.Exp(lambda);
            }
            S_N.push_back(partial / (N_tot[j]) );
            out << "," << S_N.at(j);
        }
        out << endl;
    }
    
    out.close();
    
    //Lorentziana
    out.open("out5");
    out.seekp(0);
    double mean = 0;
    double sigma = 1;
    for (int i = 0; i < M; i++){
        S_N.resize(0);
        out << i+1;
        //Calcolo 4  versioni diverse per S_N, N = 1,2,10,100
        for(int j = 0; j < dime; j++){
            partial = 0;
            for(int k = 0; k < N_tot[j]; k++){
                partial += rnd.Lorentz(sigma,mean);
            }
            S_N.push_back(partial / (N_tot[j]) );
            out << "," << S_N.at(j);
        }
        out << endl;
    }
    
    out.close();
    
    //Uniforme
    out.open("out6");
    out.seekp(0);
    for (int i = 0; i < M; i++){
        S_N.resize(0);
        out << i+1;
        //Calcolo 4  versioni diverse per S_N, N = 1,2,10,100
        for(int j = 0; j < dime; j++){
            partial = 0;
            for(int k = 0; k < N_tot[j]; k++){
                partial += rnd.Rannyu();
            }
            S_N.push_back(partial / (N_tot[j]) );
            out << "," << S_N.at(j);
        }
        out << endl;
    }
    
   out.close();
  

   rnd.SaveSeed();
   return 0;
}

void rnd_initialize(Random &rnd, int *seed){
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


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

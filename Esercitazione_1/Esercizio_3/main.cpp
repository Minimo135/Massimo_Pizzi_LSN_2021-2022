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
bool thrw(double, double, Random &);
 
int main (int argc, char *argv[]){

    Random rnd;
    int seed[4];
    rnd_initialize(rnd, seed);
    
    ofstream out;
    out.open("out7");
    
    //ES 01.3
    
    int M = 10000; //Numero totale di "lanci"
    int N = 100; //Numero di blocchi
    int L = double(M)/N;
    vector<double> PI;
    vector<double> PI2;
    vector<double> PI_final;
    vector<double> tmp2;
    vector<double> error;
    double len = 1;
    double d = 1.5;
    bool hit = false;
    int count = 0;
    double prob;
    for(int i = 0; i < N; i++){
        count = 0;
        for(int j = 0; j < L; j++){
            if(thrw(len, d, rnd)) count++;
        }
        prob = double(count)/L;
        PI.push_back(2*len/prob/d);
        PI2.push_back(pow(PI.at(i),2));
        
        PI_final.push_back( accumulate(PI.begin(), PI.end(), 0.)/PI.size() );
        tmp2.push_back( accumulate(PI2.begin(), PI2.end(), 0.)/PI2.size() );
        error.push_back( sqrt(( tmp2.at(i) - pow(PI_final.at(i), 2) )/error.size()) );
        out << (i+1)*L << "," << PI_final.at(i) << "," << error.at(i) << endl;
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

bool thrw(double len, double d, Random &rnd){
    bool hit = false;
    //Determino posizione estremo inferiore (uniformemente distribuito)
    double p = rnd.Rannyu(0,d);
    
    //Determino un angolo theta casuale
    double theta, x, y;
    bool is_inside = false;
    while(!is_inside){
        x = rnd.Rannyu();
        y = rnd.Rannyu();
        is_inside = x*x + y*y < 1;
    }
    theta = atan2(y, x);
    
    hit = p > d - len*sin(theta);

    return hit;
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

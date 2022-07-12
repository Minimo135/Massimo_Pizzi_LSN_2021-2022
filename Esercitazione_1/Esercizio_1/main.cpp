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
   int L = int(M/N); //Step per blocco
   
   vector<double> A_i(L, 0.);
   vector<double> error;
   vector<double> tmp2;
   vector<double> A;
   vector<double> A2;
   vector<double> A_final;
    
   //ES. 01.1 - 1
   out.open("out1");
   for(int i=0; i<N; i++){
       for(int j = 0; j < L; j++){
           A_i.at(j) = rnd.Rannyu();
       }
       A.push_back( accumulate(A_i.begin(), A_i.end(), 0.)/A_i.size() );
       A2.push_back( pow(A.at(i),2) );
       
       A_final.push_back( accumulate(A.begin(), A.end(), 0.)/A.size() );
       tmp2.push_back( accumulate(A2.begin(), A2.end(), 0.)/A2.size() );
       error.push_back( sqrt(( tmp2.at(i) - pow(A_final.at(i), 2) )/error.size()) );
       out << (i+1)*L << "," << A_final.at(i) << "," << error.at(i) << endl;
   }
   out.close();
    
   A.resize(0); A2.resize(0); A_final.resize(0); tmp2.resize(0); error.resize(0);
    
   //ES 01.1 - 2
   out.open("out2");
   out.seekp(0);
   for(int i=0; i<N; i++){
       for(int j = 0; j < L; j++){
           A_i.at(j) = pow(rnd.Rannyu() - 0.5,2);
       }
       A.push_back( accumulate(A_i.begin(), A_i.end(), 0.)/A_i.size() );
       A2.push_back( pow(A.at(i),2) );
       
       A_final.push_back( accumulate(A.begin(), A.end(), 0.)/A.size() );
       tmp2.push_back( accumulate(A2.begin(), A2.end(), 0.)/A2.size() );
       error.push_back( sqrt(( tmp2.at(i) - pow(A_final.at(i), 2) )/error.size()) );
       out << (i+1)*L << "," << A_final.at(i) << "," << error.at(i) << endl;
   }
   out.close();
    
   A.resize(0); A2.resize(0); A_final.resize(0); tmp2.resize(0); error.resize(0);
    
   //ES 01.1 - 3
   double ran;
   vector<double> random;
   vector<double> Chi_i;
   double Chi;
   int count = 0;
   out.open("out3");
   out.seekp(0);
   //Divido [0,1] in N sottointervalli
   for(int i = 0; i < N; i++){
       count = 0;
       //genero M numeri casuali e valuto X_i per ogni sottointervallo
       for(int j = 0; j < M; j++){
           //random.push_back(rnd.Rannyu());
           ran = rnd.Rannyu();
           if (ran > i*1./N && ran < (i+1)*1./N){
               count ++;
           }
       }
       Chi_i.push_back( pow( count - double(M)/N ,2)/(double(M)/N) );
       out << i << "," << Chi_i.at(i) << endl;
   }
   Chi = accumulate(Chi_i.begin(), Chi_i.end(), 0.);
   cout << Chi << endl;
    
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

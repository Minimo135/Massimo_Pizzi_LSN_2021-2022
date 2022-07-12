#include <iostream>
#include <fstream>
#include "SingletonRand.h"

using namespace std;


SingletonRand::SingletonRand(){
    m_rnd = new Random;
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
                m_rnd->SetRandom(seed,p1,p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;
}

SingletonRand *  SingletonRand::get_instance(){
    if (singleton_rnd == NULL) {
        singleton_rnd = new SingletonRand;
    }
    return singleton_rnd;    
}

SingletonRand::~SingletonRand(){
    delete(m_rnd);
}
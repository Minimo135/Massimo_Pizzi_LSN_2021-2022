#ifndef _GEN_
#define _GEN_

#include <vector>
#include <cstdio>
#include "global.h"

using namespace std;

//#define DEBUG_GEN

////////////////////////////////////////////////////////////////////////////////////////////////////////////
class chromosome
{
private:
  int m_dim;  //Numero totale di città nel cromosoma
  vector<int> m_cities; //Vettore con le etichette delle città, ordinate secondo il verso di percorrenza
  //vector<city> m_cities_list; //Vettore con le città generate
  double m_fit;  //Parametro che indica l'adattamento - basso è meglio.
  // double m_lenght;

  chromosome & Swap();
  chromosome & Shift();
  chromosome & Perm();
  chromosome & Inverse();
public:

    bool generated_correctly;

    chromosome(int ncity /*, double len*/ );
    chromosome();
    ~chromosome();

    bool check();
    void fprint();
    void print();
    chromosome & mutate();
    //void estimate_fit();

    chromosome & set_fit(double fit) {m_fit = fit; return *this;};
    double get_fit() { return m_fit; };
    int get_dim() {return m_dim; };
    void set_city(int pos, int j) {m_cities.at(pos) = j;};
    int get_city(int i) {return m_cities.at(i); };

    void operator = (const chromosome &D ) { 
         m_dim = D.m_dim;
         m_fit = D.m_fit;
         m_cities = D.m_cities;
      }
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Population
{
private:
  int m_dim_pop; //Il numero di cromosomi totali
  vector<chromosome> m_population; //array dei cromosomi (la popolazione)
  vector<city> m_cities_list; //Lsta delle città

public:
  Population(int dim_pop);
  ~Population();

  chromosome get_chromosome(int i) {return m_population.at(i);};
  vector<chromosome> & get_chromosomes() {return m_population; };
  chromosome * get_best();
  vector<city> get_cities(){return m_cities_list;}
  //chromosome select();
  vector<chromosome> select(bool );

  Population & set_population(vector<chromosome> new_pop) {
    #ifdef DEBUG
    printf("Sostituzione pop. Dim vecchia = %lu, dim nuova = %lu", m_population.size(), new_pop.size());
    #endif
    m_population = new_pop;
    return *this;
  };
  Population & fprint_cities();
  Population & print_chromosomes();
  Population & estimate_fit();
  double estimate_half_av_fit();

  void operator = (const Population &P ) { 
         m_cities_list = P.m_cities_list;
         m_dim_pop = P.m_dim_pop;
         m_population = P.m_population;
      }

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool crossingover(vector<chromosome> &parents);
bool Compare(chromosome a, chromosome b);
double estimate_fit(int * get_city, int size, vector<city> m_cities_list);

#endif
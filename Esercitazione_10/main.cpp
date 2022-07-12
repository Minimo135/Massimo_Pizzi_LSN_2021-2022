#include "mpi.h"

#include <vector>
#include <cmath>
#include <climits>
#include "global.h"
#include "Genetics.h"

// #define DEBUG

using namespace std;

Random rnd;

int main (int argc, char *argv[]){

  int size, rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status stat[size];
  MPI_Request req[size];

  // if (size == 1){
  //   cout << "Sono l'unico nodo che hai utilizzato!" << endl;
  // }
  // else {
  //   cout << "Sono il nodo " << rank << " dei " << size << " che hai utilizzato!" << endl;
  // }

  int to_nodes[size];
  int from_nodes[size];
  

  // cout << "Nodo " << rank << ", mando a " << to_nodes[rank] << ", ricevo da " << from_nodes[rank] << endl;


  initialize_rnd(rnd, rank);

  int dim_pop = 300;
  int nstep = 500;
  int N_migr = 20;
  vector<double> best_fit;

  chromosome parent1, parent2;
  vector<chromosome> parents;
  chromosome child1, child2;
  vector<chromosome> new_population;

  ///////////////////////////////////////////
  // Città Americane //
  ///////////////////////////////////////////
  //Genero la popolazione di "cromosomi" iniziale
  Population population(dim_pop);
  chromosome *best = population.get_best();
  ofstream best_out;
  //Estraggo il migliore all'inizio (solo nodo 0)
  if(rank == 0){
    best_out.open("Best_pre.dat");
    for(int i = 0; i < best->get_dim(); i++){
      best_out << population.get_cities().at( best->get_city(i) ).x << "," << population.get_cities().at( best->get_city(i) ).y << endl;
    }
    best_out.close();
  }

  population.estimate_fit();

  //Inizio l'algoritmo
  for(int i = 0; i < nstep; i++){
    #ifdef DEBUG
    cout << "\n Inizio step " << i << "..." << endl;
    #endif
    //Finché non ho generato un'intera nuova popolazione:
    //1) Scelgo due cromosomi, con una probabilità maggiore se meglio adattati
    //2) Muto i geni scelti
    //3) Eseguo il crossing-over
    //4) Sostituisco i due cromosomi di partenza con i "figli"
    new_population.clear();
    for(int i = 0; i < dim_pop/2; i++){
      #ifdef DEBUG
      cout << "\nRicombino la popolazione, passo " << i << " di " << dim_pop/2 - 1 << endl;
      #endif
      parents = population.select(1);
      parents.at(0).mutate();
      parents.at(1).mutate();
      crossingover(parents);
      if (!parents.at(0).check() or !parents.at(1).check()){
        cerr << "Errore in mutazione o crossingover" << endl;
        parents.at(0).print(); parents.at(1).print();
        return -1;
      }
      new_population.push_back(parents.at(0)); new_population.push_back(parents.at(1));
    }
    #ifdef DEBUG
    cout << "Generata nuova popolazione, sostituisco la precedente...\n";
    #endif
    population.set_population(new_population);
    population.estimate_fit();

    //Scambio popolazione migliore tra nodi ogni N_migr
    if(!(i%N_migr)){
      //Genero coppie per lo scambio
      Select_pairs(to_nodes, from_nodes, rank, size);
      //Estraggo il cromosoma migliore (il suo percorso in paricolare)
      best = population.get_best();
      
      // cout << "Step " << i << ", nodo " << rank << " - fit del best "  << best->get_fit() << endl; 

      int lenght = best->get_dim();
      int message[lenght];
      int recvmsg[lenght];
      for(int k = 0; k < lenght; k++){
        message[k] = best->get_city(k);
      }

      //Effettuo gli scambi
      MPI_Send(&message[0], lenght, MPI_INTEGER, to_nodes[rank], to_nodes[rank], MPI_COMM_WORLD/* , &req[rank] */);
      MPI_Recv(&recvmsg[0], lenght, MPI_INTEGER, from_nodes[rank], rank, MPI_COMM_WORLD, &stat[rank]);

      // cout << "Step " << i << ", nodo " << rank << " - Inviato a nodo " << to_nodes[rank] << ", ricevuto da nodo " << from_nodes[rank] << "! " << endl; 

      //Sostituisco il migliore con quello ricevuto
      for(int k = 0; k < lenght; k++){
        best->set_city(k, recvmsg[k]);
      }

      population.get_best()->set_fit(0);
      population.estimate_fit();

      // cout << "Step " << i << ", nodo " << rank << " - fit del best "  << best->get_fit() << endl; 
    }

    if(not (i%10) ){
      cout << "Elaboro ... " << (int)( (double)i/nstep*100 ) << "% \r" << flush;
    }

    best_fit.push_back(population.get_best()->get_fit());

    #ifdef DEBUG
    cout << "\nLa nuova popolazione allo step " << i <<":" << endl;
    for(int i = 0; i < dim_pop; i++){
    cout << i << ": [";
    for(int j = 0; j < 50; j++){
      cout << population.get_chromosome(i).get_city(j) << " ";
    }
    cout << "]" << ", fit = " << population.get_chromosome(i).get_fit() << endl;
    }
    #endif

  }

  population.estimate_fit();

  //Se c'è più di un nodo, allora il nodo 0 trova il migliore tra i percorsi rovati
  if(size > 1){
    // Il nodo 0 raccoglie i migliori cromosomi da tutti i nodi e prende quello con fit maggiore.
    best = population.get_best();
    int lenght = best->get_dim();
    int message[lenght];
    int allmsg[(size - 1)*lenght];
    int *bests[size-1];
    for(int k = 0; k < lenght; k++){
          message[k] = best->get_city(k);
    }
    MPI_Gather(&message[0], lenght, MPI_INTEGER, &allmsg[0], lenght, MPI_INTEGER, 0, MPI_COMM_WORLD);
    if(rank == 0){
      // cout << endl << "[";
      // for(int i = 1; i <= (size - 1)*lenght; i++){
      //   cout << allmsg[i-1] << ", ";
      //   if ( !(i%lenght) ) cout  << "]" << endl << "[";
      // }
      int k = 0;
      for(int i = 1; i <= (size - 1)*lenght; i++){
        if ( !(i%lenght) ){
          bests[k] = allmsg + k*lenght;
          k ++;
        }
      }

      // for(int i = 0; i < size - 1; i++){
      //   cout << endl << "[";
      //   for(int j = 0; j < lenght; j++){
      //     cout << bests[i][j] << ", ";
      //   }
      //   cout << "]" << endl;
      // }

      //Valuto il percorso migliore
      double fit = 0;
      double min = INT_MAX;
      int index = -1;
      for(int i = 0; i < size - 1; i++){
        fit = estimate_fit(bests[i], population.get_cities().size(), population.get_cities());
        // cout << "fit = " << fit << ", passo "cleclear << i << endl;
        if(fit < min){
          min = fit;
          index = i;
        }
      }

      for(int i = 0; i < population.get_cities().size(); i++){
        best->set_city(i, bests[index][i]);
      }
      best ->set_fit(fit);
      // Stampo su schermo il percorso migliore:
      cout << "Il percorso migliore trovato: "<<endl;
      best->print();

    }
  }
  else {
    best = population.get_best();
    // Stampo su schermo il percorso migliore:
    cout << "Il percorso migliore trovato: "<<endl;
    best->print();
  }

  cout << endl;

  


  //Prendo il miglior cromosoma e ne stampo le città (solo nodo 0)
  if(rank == 0){
    best_out.open("Best.dat");
    for(int i = 0; i < best->get_dim(); i++){
      best_out << population.get_cities().at( best->get_city(i) ).x << "," << population.get_cities().at( best->get_city(i) ).y << endl;
    }
    best_out.close();
  }

  //Genero il file per il grafico della variazione del fit migliore ad ogni iterazione
  ofstream best_fit_out;
  string name("best_fit_out");
  name.append(to_string(rank));
  name.append(".dat");
  best_fit_out.open(name);
  for(uint i = 0; i < best_fit.size(); i++){
    best_fit_out << i << "," << best_fit.at(i) << endl;
  }
  best_fit_out.close();

  rnd.SaveSeed();

  MPI_Finalize();

  return 0;

}

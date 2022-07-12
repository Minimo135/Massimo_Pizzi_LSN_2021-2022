#include "global.h"

using namespace std;

int mod(int n, int m){
   //Funzione modulo personalizzata: ignora la prima posizione della lista, che per costruzione deve essere fissata
  if(n >= m) return n%m + 1;
  return n;
}


void initialize_rnd(Random &rnd, int ii){
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      for(int i = 0; i < ii; i++){
         Primes >> p1 >> p2 ;
      }
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

void Select_pairs(int *to_nodes, int* from_nodes, int rank, int size){
   //Il nodo 0 decide in maniera casuale come assegnare le coppie di nodi per lo scambio del percorso migliore
  if(rank == 0){
    //Inserisco in to_nodes l'elenco dei nodi in ordine crescente
    for(int i = 0; i < size; i++){
      to_nodes[i] = i;
    }
    //Eseguo uno shuffle dell'array
    //Rimescolo casualmente le città
    int tmp = -1;
    for(int i = 0; i < size; i++){
      int new_pos = int(rnd.Rannyu()*size);
      //Scambio l'elemento i con il new_pos-esimo
      tmp = to_nodes[new_pos];
      to_nodes[new_pos] = to_nodes[i];
      to_nodes[i] = tmp;
    }
    // Ora l'array to_nodes specifica le coppie: il nodo in posizione i da ilproprio percorso a to_nodes[i]
    // Per semplicità costruisco anche il vettore from_nodes, che rappresenta il nodo da cui ricevere la città
    for(int i = 0; i < size; i++){
      from_nodes[to_nodes[i]] = i;
    }
  }
  else{
    for(int i = 0; i < size; i++){
      from_nodes[i] = to_nodes[i] = 0;
    }
  }
  //Il nodo 0 invia in broadcast a tutti gli altri nodi questi due array
  MPI_Bcast(to_nodes, size, MPI_INTEGER, 0, MPI_COMM_WORLD);
  MPI_Bcast(from_nodes, size, MPI_INTEGER, 0, MPI_COMM_WORLD);
}
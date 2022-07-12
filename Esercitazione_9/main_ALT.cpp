#include <vector>
#include <cmath>
#include "global.h"
#include "Genetics.h"

// Decommentare sotto per usare il debug
// #define DEBUG

using namespace std;

Random rnd;

int main (int argc, char *argv[]){

  initialize_rnd(rnd);
  //auto rnd = SingletonRand::get_instance()->m_rnd;

  double r = 5;
  int ncity = 34;
  int dim_pop = 1500;
  int nstep = 1000;

  vector<double> best_fit;
  vector<double> best_av_fit;

  vector<chromosome> parents;
  vector<chromosome> new_population;

  ///////////////////////////////////////////
  // 1) Caso di città su una circonferenza //
  ///////////////////////////////////////////
  bool on_circ = true; //Flag per generare punti su una circonferenza

  //Genero la popolazione di "cromosomi" iniziale
  Population population(on_circ, ncity, dim_pop, r);
  
  //Estraggo il miglior percorso all'inizio
  chromosome best = population.estimate_fit().get_best();
  ofstream best_out;
  //Scrivo la configurazione migliore di partenza su file
  best_out.open("Best_pre.dat");
  for(int i = 0; i < best.get_dim(); i++){
    best_out << population.get_cities().at( best.get_city(i) ).x << "," << population.get_cities().at( best.get_city(i) ).y << endl;
  }
  best_out.close();

  cout << "Caso 1: punti su una circonferenza" << endl;
  // Stampo su terminale i primi 3 e gli ultimi 3 cromosomi generati
  cout << "Sono stati generati " << dim_pop << " individui:" << endl;
  for(int i = 0; i < 3; i++){
    cout << i << ": ";
    population.get_chromosome(i).print();
    cout << endl;
  }
  cout << "..." << endl;
  for(int i = dim_pop - 3; i < dim_pop; i++){
    cout << i << ": ";
    population.get_chromosome(i).print();
    cout << endl;
  }

  //Inizio l'algoritmo
  for(int k = 0; k < nstep; k++){
    #ifdef DEBUG
    cout << "\n Inizio step " << k << "..." << endl;
    #endif
    //Finché non ho generato un'intera nuova popolazione:
    //1) Scelgo due cromosomi
    //2) Eseguo il crossing-over
    //3) Muto i geni dei "figli"
    //4) Sostituisco i due cromosomi di partenza con i "figli"
    new_population.clear(); //Azzero la popolazione (qui per sicurezza)
    //Ciclo di generazione della nuova popolazione - metodo (0)
    for(int i = 0; i < dim_pop/2; i++){

      #ifdef DEBUG
      //cout << "\nRicombino la popolazione, passo " << i << " di " << dim_pop/2 - 1 << endl;
      #endif

      //Seleziono la coppia di genitori. La funzione prende in ingresso la flag tournament
      parents = population.select(0);
      //Eseguo il crossingover
      crossingover(parents);
      //Muto i cromosomi e contemporaneamente verifico che siano stati mutati correttamente
      if (!parents.at(0).mutate().check() or !parents.at(1).mutate().check()){
        cerr << "Errore in mutazione o crossingover" << endl;
        parents.at(0).print(); parents.at(1).print();
        return -1;
      }
      // Inserisco i figli appena generati all'interno della nuova popolazione
      new_population.push_back(parents.at(0)); new_population.push_back(parents.at(1));
    }
    #ifdef DEBUG
    cout << "Generata nuova popolazione, sostituisco la precedente...\n";
    #endif

    //Ho terminato la generazione di un'intera nuova popolazione, allora la sovrascrivo a quella attuale e stimo di nuovo il fit
    population.set_population(new_population).estimate_fit();

    // Inserisco nei vettori per gli output su file la lunghezza del miglior percorso trovato in questo step
    // e la media delle lunghezze eseguita su metà della popolazione (la metà migliore)
    best_fit.push_back(population.get_best().get_fit());
    best_av_fit.push_back(population.estimate_half_av_fit());

    // Mostro su terminale la percentuale di completamento, assieme alla lunghezza migliore trovata
    if( !(k%10)){
      cout << "\r" << int(100 * (double)k/nstep) << "% ";
      cout << "- lunghezza migliore (metodo 0): " << population.get_best().get_fit() << " "  << flush;
    }

    #ifdef DEBUG
    cout << "\nLa nuova popolazione allo step " << k <<":" << endl;
    for(int i = 0; i < dim_pop; i++){
    cout << i << ": [";
    for(int j = 0; j < ncity; j++){
      cout << population.get_chromosome(i).get_city(j) << " ";
    }
    cout << "]" << ", fit = " << population.get_chromosome(i).get_fit() << endl;
    }
    #endif

  }

  // Stimo il fit di ogni cromosoma sull'ultima popolazione generata
  population.estimate_fit();

  //Mostro su terminale i primi tre e gli ultimi tre cromosomi della popolazione finale
  cout << "\nPopolazione finale (metodo 0):" << endl;
  for(int i = 0; i < 3; i++){
    cout << i << ": ";
    population.get_chromosome(i).print();
    cout << endl;
  }
  cout << "..." << endl;
  for(int i = dim_pop - 3; i < dim_pop; i++){
    cout << i << ": ";
    population.get_chromosome(i).print();
    cout << endl;
  }

  //Prendo il miglior cromosoma e ne stampo le città (con le coordinate) su file
  best = population.get_best();
  best_out.open("Best.dat");
  for(int i = 0; i < best.get_dim(); i++){
    best_out << population.get_cities().at( best.get_city(i) ).x << "," << population.get_cities().at( best.get_city(i) ).y << endl;
  }
  best_out.close();

  //Genero il file per il grafico della variazione della lunghezza migliore ad ogni iterazione
  ofstream best_fit_out;
  best_fit_out.open("best_fit_out.dat");
  for(uint i = 0; i < best_fit.size(); i++){
    best_fit_out << i << "," << best_fit.at(i) << endl;
  }
  best_fit_out.close();

  //Genero il file per il grafico della variazione della lunghezza media sulal metà migliore della popolazione ad ogni iterazione
  ofstream best_av_fit_out;
  best_av_fit_out.open("best_av_fit_out.dat");
  for(uint i = 0; i < best_fit.size(); i++){
    best_av_fit_out << i << "," << best_av_fit.at(i) << endl;
  }
  best_av_fit_out.close();


  /////////////////////////////////////////////////
  // 2) Caso di città all'interno di un quadrato //
  /////////////////////////////////////////////////
  on_circ = false;
  best_fit.resize(0);
  best_av_fit.resize(0);

  //Genero la popolazione di "cromosomi" iniziale, partendo dalla lista di città generate dal programma principale
  Population population2("city_list.dat", dim_pop);
  // population2.fprint_cities();
  
  //Estraggo il miglior percorso all'inizio e scrivo su file le posizioni delle città per farne il grafico
  best = population2.estimate_fit().get_best();
  best_out.open("Best_pre_square.dat");
  for(int i = 0; i < best.get_dim(); i++){
    best_out << population2.get_cities().at( best.get_city(i) ).x << "," << population2.get_cities().at( best.get_city(i) ).y << endl;
  }
  best_out.close();

  cout << "\nCaso 2: punti all'interno di un quadrato" << endl;
  
  //Mostro su terminale i primi tre e gli ultimi tre cromosomi della popolazione casuale iniziale
  cout << "Sono stati generati i segeuenti individui:" << endl;
  for(int i = 0; i < 3; i++){
    cout << i << ": ";
    population2.get_chromosome(i).print();
    cout << endl;
  }
  cout << "..." << endl;
  for(int i = dim_pop-3; i < dim_pop; i++){
    cout << i << ": ";
    population2.get_chromosome(i).print();
    cout << endl;
  }

  //Inizio l'algoritmo
  for(int i = 0; i < nstep; i++){
    #ifdef DEBUG
    cout << "\n Inizio step " << i << "..." << endl;
    #endif

    //Finché non ho generato un'intera nuova popolazione:
    //1) Scelgo due cromosomi
    //2) Eseguo il crossing-over
    //3) Muto i geni scelti
    //4) Sostituisco i due cromosomi di partenza con i "figli"

    //Azzero la nuova popolazione 
    new_population.clear();

    //Ciclo di generazione nuova popolazione - metodo 0
    for(int i = 0; i < dim_pop/2; i++){
      #ifdef DEBUG
      cout << "\nRicombino la popolazione, passo " << i << " di " << dim_pop/2 - 1 << endl;
      #endif

      //Seleziono i due genutori
      parents = population2.select(0);
      // Eseguo il crossingover
      crossingover(parents);
      //Muto i cromosomi e verifico contemporaneamente che siano stati mutati correttamente
      if (!parents.at(0).mutate().check() or !parents.at(1).mutate().check()){
        cerr << "Errore in mutazione o crossingover" << endl;
        parents.at(0).print(); parents.at(1).print();
        return -1;
      }

      // Inserisco i due nuovi figli nella nuova popolazione
      new_population.push_back(parents.at(0)); new_population.push_back(parents.at(1));
    }

    #ifdef DEBUG
    cout << "Generata nuova popolazione, sostituisco la precedente...\n";
    #endif

    //Generata un'intera nuova popolazione, la sostituisco a quella precedente
    population2.set_population(new_population).estimate_fit();

    // Inserisco nei vettori per gli output su file la lunghezza migliore trovata in questa iterazione
    // E la lunghezza media della metà migliore della popolazione
    best_fit.push_back(population2.get_best().get_fit());
    best_av_fit.push_back(population2.estimate_half_av_fit());

    //Stampo a video il progresso e la lunghezza dell'attuale percorso migliore
    if( !(i%10)){
      cout << "\r" << int(100 * (double)i/nstep) << "% ";
      cout << "- lunghezza migliore (metodo 0 - 1): " << population2.get_best().get_fit() << " " << flush;
    }


    #ifdef DEBUG
    cout << "\nLa nuova popolazione allo step " << i <<":" << endl;
    for(int i = 0; i < dim_pop; i++){
    cout << i << ": [";
    for(int j = 0; j < ncity; j++){
      cout << population2.get_chromosome(i).get_city(j) << " ";
    }
    cout << "]" << ", fit = " << population2.get_chromosome(i).get_fit() << endl;
    }
    #endif

  }

  // Calcolo il fit di ogni cromosoma della popolazione finale
  population2.estimate_fit();

  //Stampo su terminale i primi tre e gli ultimi tre cromosomi della popolazione finale
  cout << "\nPopolazione finale (metodo 0):" << endl;
  for(int i = 0; i < 3; i++){
    cout << i << ": ";
    population2.get_chromosome(i).print();
    cout << endl;
  }
  cout << "..." << endl;
  for(int i = dim_pop-3; i < dim_pop; i++){
    cout << i << ": ";
    population2.get_chromosome(i).print();
    cout << endl;
  }

  //Prendo il miglior cromosoma e ne stampo le città su file
  best = population2.get_best();
  best_out.open("Best_square.dat");
  for(int i = 0; i < best.get_dim(); i++){
    best_out << population2.get_cities().at( best.get_city(i) ).x << "," << population2.get_cities().at( best.get_city(i) ).y << endl;
  }
  best_out.close();

  // Genero il file per il grafico della variazione della lunghezza migliore ad ogni iterazione
  // e della variazione della lunghezza medio sulla metà migliore della popolazione ad ogni iterazione.
  best_fit_out.open("best_fit_out_square.dat");
  for(uint i = 0; i < best_fit.size(); i++){
    best_fit_out << i << "," << best_fit.at(i) << endl;
  }
  best_fit_out.close();

  best_av_fit_out.open("best_av_fit_out_square.dat");
  for(uint i = 0; i < best_av_fit.size(); i++){
    best_av_fit_out << i << "," << best_av_fit.at(i) << endl;
  }
  best_av_fit_out.close();



  rnd.SaveSeed();
  return 0;

}

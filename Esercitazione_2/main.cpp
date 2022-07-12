#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"

using namespace std;

double Integrale(Random &, bool);

struct posizione{
  double x;
  double y;
  double z;
};

void initialize_rnd(Random &);
void incremento(vector<posizione> & traj, Random & rnd);
void incremento_cont(vector<posizione> & traj, Random & rnd);
double norma2(posizione punto);
double media(int passo, int end, vector<vector<double>> vec);
double errore(int passo, int N, vector<vector<double>> medie_blocchi, vector<vector<double>> medie_blocchi2);

posizione random_sfera(Random & rnd);

int main (int argc, char *argv[]){

  Random rnd;
  initialize_rnd(rnd);
  int N = 100;
  int M = 50000;
  int L = M/N;
    
  // Esercizio 1: Integrale Montecarlo
  cout << "Esercizio 1: stima dell'integrale con importance sampling." << endl;
  // Si vuole calcolare I = \int_0^1 PI/2 * cos(PI*x/2)dx [ = 1]
  //Ciclo sui blocchi
  double  average = 0, averageIMP = 0, average2 = 0, averageIMP2 = 0, err, errIMP;
    
  ofstream Int, Int_IMP;
  Int.open("Integral.dat");
  Int_IMP.open("Integral_IMP.dat");
    
  for(int k = 0; k < N; k++){
    double est_I = 0;
    double est_I_IMP = 0;
    //Ciclo nel blocco
    for(int j = 0; j < L; j++){
      est_I += Integrale(rnd, 0);
      est_I_IMP += Integrale(rnd, 1);
    }
    est_I = est_I/L;
    est_I_IMP = est_I_IMP/L;
    average += est_I;
    average2 += est_I*est_I;
    averageIMP += est_I_IMP;
    averageIMP2 += est_I_IMP*est_I_IMP;
    
    err = k?sqrt( (average2 / (k+1) - pow(average / (k+1), 2) )/k ) : 0;
    errIMP = k?sqrt( (averageIMP2 / (k+1) - pow(averageIMP / (k+1), 2))/k ) : 0;
      
    Int << k + 1 <<"," << average / (k+1) << "," << err << endl;
    Int_IMP << k + 1 << "," << averageIMP / (k+1) << "," << errIMP << endl;
      
    cout << "Progresso: blocco " << k << " di " << N << " - " << int( 100*(double)k/N ) << "% \r" << flush;
    
  }
    
    cout << endl;
    
    Int.close();
    Int_IMP.close();

    
  //Esercizio 2: Random Walk
  cout << "Esercizio 2: random walk." << endl;

  ofstream out;
  out.open("ES2_rndw.dat");
    
  int nstep = 100;

  vector<vector<vector<posizione>>> rndw;
  rndw.resize(N);
  for(int k = 0; k<N; k++){
    rndw.at(k).resize(L);
  }

  vector<vector<double>> medie_blocchi;
  medie_blocchi.resize(N);
  for(int k = 0; k < N; k++){
    medie_blocchi.at(k).resize(nstep);
  }
   
  vector<vector<double>> medie_blocchi2;
  medie_blocchi2.resize(N);
  for(int k = 0; k < N; k++){
    medie_blocchi2.at(k).resize(nstep);
  }

  vector<double> media_prog;
  media_prog.resize(nstep);
   
  vector<double> err_prog;
  err_prog.resize(nstep);

  double media_par;

  /*
  Parte a)
  Random walk su array cubico
  */

  //Ciclo sugli nstep
  for(int i = 0; i < nstep; i++){
    //Ciclo sui blocchi
    for(int k = 0; k < N; k++){
      media_par = 0;
      //Ciclo nel blocco
      for(int j = 0; j < L; j++){
        if(i == 0){
          posizione init{0,0,0};
          rndw.at(k).at(j).push_back(init);
        }
        incremento(rndw.at(k).at(j), rnd);
        media_par += norma2(rndw.at(k).at(j).at(i));
      }
      medie_blocchi.at(k).at(i) = media_par/L;
      medie_blocchi2.at(k).at(i) = pow(media_par/L, 2);
      cout << "Progresso: blocco " << k << " di " << N << " - ";
      cout << "step " << i << " di " << nstep << " - " << int( 100*(double)k/N ) << "% \r" << flush;
    }
    media_prog.at(i) = sqrt(media(i,N,medie_blocchi));
    err_prog.at(i) = errore(i, N, medie_blocchi, medie_blocchi2)/( 2*sqrt(media_prog.at(i)) );
    out << i << "," << media_prog.at(i) << "," << err_prog.at(i) << endl;
  }
    
    cout << endl;

  //Cancello i random walk e apro un nuovo file...
  for(int k = 0; k < N; k++){
    for(int j = 0; j < L; j++){
      rndw.at(k).at(j).resize(0);
    }
  }
  out.close();
  out.open("ES2_rndwc.dat");
  
  /*
  Parte b)
  Random walk nel continuo
  */
  //Ciclo sugli nstep
  for(int i = 0; i < nstep; i++){
    //Ciclo sui blocchi
    for(int k = 0; k < N; k++){
      media_par = 0;
      //Ciclo nel blocco
      for(int j = 0; j < L; j++){
        if(i == 0){
          posizione init{0,0,0};
          rndw.at(k).at(j).push_back(init);
        }
        incremento_cont(rndw.at(k).at(j), rnd);
        media_par += norma2(rndw.at(k).at(j).at(i));
      }
      medie_blocchi.at(k).at(i) = media_par/L;
      medie_blocchi2.at(k).at(i) = pow(media_par/L, 2);
    }
    media_prog.at(i) = sqrt(media(i,N,medie_blocchi));
    err_prog.at(i) = errore(i, N, medie_blocchi, medie_blocchi2)/( 2*sqrt(media_prog.at(i)) );
    out << i << "," << media_prog.at(i) << "," << err_prog.at(i) << endl;
  }
  out.close();

  //Salvo su un file qualche rndw
  cout << endl << "Scrivo su file le posizioni di " << L << " random walk dal bloccho #" << 4 << " ..." << endl;
    
  ofstream WriteXYZ;
  for(int i = 0; i < L; i++){
    WriteXYZ.open("./rndws/rndw_" + to_string(i) + ".xyz");
    for(int j = 0; j < nstep; j++){
      WriteXYZ << rndw.at(4).at(i).at(j).x << "," << rndw.at(4).at(i).at(j).y << "," << rndw.at(4).at(i).at(j).z << endl;
    }
    if (!(i%20)) cout << int(100*(double)i/L) << "% \r" << flush;
    WriteXYZ.close();
  }

  out.close();
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

double Integrale(Random & rnd, bool iIMPORTANCE){
    double result = 0;
    double thrw, d, p;
    if(iIMPORTANCE){
        // Importance sampling
        // Probabilità che si vuole usare: p(x) = 2 - 2x
        // Inversa della cumulativa: y(x) = 1 - sqrt(1 - x)
        thrw = rnd.Rannyu();
        d = 1 - sqrt(1 - thrw);
        p = 2 - 2*d;
        result = M_PI/2*cos(M_PI*d/2.)/p;
    }
    else{
        // Metodo della media (importance sampling con distribuzione uniforme)
        thrw = rnd.Rannyu();
        result = M_PI/2*cos(M_PI*thrw/2.);
    }
    return result;
}

void incremento(vector<posizione> & traj, Random & rnd){
  //La funzione aggiunge una posizione ad un array di posizioni, calcolata sulla base di quella precedente.
  //Una delle tre direzioni x,y,z, a caso, viene incrementata di +1 o -1, scelti sempre a caso
  int dice = rnd.Rannyu()*3;
  int versod = rnd.Rannyu()*2;
  int verso;
  posizione punto = traj.at(traj.size()-1);

  //cout << versod << endl;

  if(versod == 0) {verso = -1;}
  else if(versod == 1) {verso = 1;}
  if(dice == 0){
    punto.x += verso;
  }
  else if(dice == 1){
    punto.y += verso;
  }else if(dice == 2){
    punto.z += verso;
  }
  traj.push_back(punto);
}

void incremento_cont(vector<posizione> & traj, Random & rnd){
  //La funzione aggiunge una posizione ad un array di posizioni, calcolata sulla base di quella precedente.
  //La nuova posizione si trova a distanza l = 1 dlla precedente, ad una posizione angolare casuale.
  posizione next = random_sfera(rnd);
  posizione punto = traj.at(traj.size()-1);

  punto.x += next.x;
  punto.y += next.y;
  punto.z += next.z;
  
  traj.push_back(punto);
}

double norma2(posizione punto){
  return punto.x*punto.x + punto.y*punto.y + punto.z*punto.z;
}

double media(int passo, int end, vector<vector<double>> vec){
  double parziale = 0.;
  for(int k = 0; k < end; k++){
    parziale += vec.at(k).at(passo);
  }
  return parziale / double(end);
}

double errore(int passo, int N, vector<vector<double>> medie_blocchi, vector<vector<double>> medie_blocchi2){
	double parziale = 0.;
	double parziale2 = 0.;
	for(int k = 0; k < N; k++){
		parziale += medie_blocchi.at(k).at(passo);
		parziale2 += medie_blocchi2.at(k).at(passo);
	}
	
	return sqrt( (parziale2/double(N) - pow(parziale/double(N), 2) ) /double(N-1) );
}

posizione random_sfera(Random & rnd){
  double x = rnd.Rannyu()*2-1;
  double y = rnd.Rannyu()*2-1;
  double z = rnd.Rannyu()*2-1;
  while(x*x + y*y + z*z > 1){
    x = rnd.Rannyu()*2-1;
    y = rnd.Rannyu()*2-1;
    z = rnd.Rannyu()*2-1;
  }

  int dir = rnd.Rannyu()*2;
  x = (dir?-1:1) * sqrt( x*x/ (x*x + y*y + z*z) );
  dir = rnd.Rannyu()*2;
  y = (dir?-1:1) * sqrt( y*y/ (x*x + y*y + z*z) );
  dir = rnd.Rannyu()*2;
  z = (dir?-1:1) * sqrt( z*z/ (x*x + y*y + z*z) );

  posizione punto{x,y,z};
  return punto;
}
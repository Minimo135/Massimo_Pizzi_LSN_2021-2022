//#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iterator>
#include "global.h"
#include "Genetics.h"


using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////
chromosome::chromosome(int ncity /*, double len*/)
{
  
  //m_lenght = len;
  m_dim = ncity;
  //La prima città da visitare è sempre la numero 1
  //Genero un array di città in ordine
  for(int i = 0; i < m_dim; i++){
    m_cities.push_back(i);
    //cout << "controllo: " << i << endl;
  }

  //Rimescolo casualmente le città (lasciando la prima inalterata)
  for(int i = 1; i < m_dim; i++){
    int new_pos = 1 + int(rnd.Rannyu()*(m_dim-1));
    //Scambio l'elemento i con il new_pos-esimo
    int tmp = m_cities.at(new_pos);
    m_cities.at(new_pos) = m_cities.at(i);
    m_cities.at(i) = tmp;
  }
  generated_correctly = check();
  #ifdef DEBUG_GEN
  cout << flush << "Ho generato:"; this->print();
  #endif
}

chromosome::chromosome()
{
  
}

bool chromosome::check(){
  //Verifico che il cromosoma sia ben costruito
  //Verifico che la prima città sia la numero 0
  if (m_cities.at(0)) return false;

  //Per ogni città, verifico che sia unica
  for(int i = 1; i < m_dim; i++){
    for(int j = i + 1; j < m_dim; j++){
      if (m_cities.at(i) == m_cities.at(j)) return false;
    }
  }

  return true;
}

/*void chromosome::estimate_fit(){
  
}*/

void chromosome::fprint(){
  ofstream out;
  out.open("chromosome_list.dat", ios::app);
  
  out << m_cities.at(0);
  for(int i = 1; i < m_dim; i++){
    out << "," << m_cities.at(i);
  }
  out << endl;

  out.close();
}

void chromosome::print(){

  cout << "[ " << m_cities.at(0);
  for(int i = 1; i < m_dim; i++){
    cout << " " << m_cities.at(i);
  }
  cout << " ], fit = " << m_fit;
}

chromosome & chromosome::mutate(){

  //Decido di mutare i cromosomi con una probabilità proporzionale al loro livello di adattamento:
  //I meno adatti mutano con più probabilità
  int np = 4;
  double perm[np] {0.015, 0.001, 0.01, 0.015};

  #ifdef DEBUG_GEN
  cout << "\nTentativo di mutazione:" << endl;
  cout << "\tScambio - probabilità: " << perm[0] << endl;
  cout << "\tScorrimento - probabilità: " << perm[1] << endl;
  cout << "\tPermutazione - probabilità: " << perm[2] << endl;
  cout << "\t\Inversione - probabilità: " << perm[3] << endl;
  #endif
  
  bool chosenp [np] {false, false, false, false}; //0: swap, 1: shift, 2: perm, 3: inverse
  //Funzione di selezione della permutazione
  double threshold = 0;

  for(int k = 0; k < np; k++){
    double thrw = rnd.Rannyu();
    for(int i = 0; i < np; i++ ){
      threshold += perm[i];
      if(thrw < threshold){
        chosenp[i] = true;
        break;
      }
    }
  }

  //Mutazione di scambios
  if(chosenp[0]) {
    this->Swap();
  }

  //Mutazione di tipo 2: shift di M posizioni di un gruppo di N città verso sinistra
  if(chosenp[1]){
    this->Shift();
  }

  //Mutazione di tipo 3: permutazione di M città contigue con altre M città contigue
  if(chosenp[2]){
    this->Perm();
  }

  if(chosenp[3]){
    this->Inverse();
  }

  bool fail[] = {false, false, false, false};
  if (chosenp == fail){
    #ifdef DEBUG_GEN
    cout << "Tentativo di mutare fallito!\n";
    #endif
    
    return *this;
  }

  return *this;

}

chromosome & chromosome::Swap(){
    #ifdef DEBUG_GEN
    cout << "\tMuto il cromosoma (scambio)." << endl;
    cout << "Prima: "; this->print();
    #endif

    int site1, site2;
    site1 = int( rnd.Rannyu()*(m_dim-1) ) + 1;
    do{
      site2 = int( rnd.Rannyu()*(m_dim-1) ) + 1;
    }while(site2 == site1);

    double tmp = m_cities.at(site1);
    m_cities.at(site1) = m_cities.at(site2);
    m_cities.at(site2) = tmp;

    #ifdef DEBUG_GEN
    cout << endl << "Dopo:  "; this->print();
    #endif

    generated_correctly = check();

    return *this;
}

chromosome & chromosome::Shift(){
int M = int( rnd.Rannyu()* (m_dim - 3) ) + 1; //Lunghezza di scorrimento
    int N = int( rnd.Rannyu()* (m_dim - 2) ) + 2; //# di città da spostare
    int start = int( rnd.Rannyu()* (m_dim - 1) ) + 1; //Inizio del blocco da spostare

    #ifdef DEBUG_GEN
    cout << "\tMuto il cromosoma (scorrimento)." << endl;
    cout << "Sposto " << N << " città di " << M << " posti, partendo dal sito " << start << endl;
    cout << "Prima: "; this->print();
    #endif

    int new_pos;
    vector<int> new_cities(m_dim, -1); new_cities.resize(m_dim);
    new_cities.at(0) = 0;
    //Copio le città nella nuova posizione traslata
    for(int i = start; i < start + N; i++){
      new_pos = mod(mod(i, m_dim) + M, m_dim);
      //new_pos = (i%m_dim + 1 + M)%m_dim + 1;
      //cout << "copio la città (" << m_cities.at(mod(i, m_dim)) << ") in posizione " << new_pos << endl;
      new_cities.at(new_pos) = m_cities.at(mod(i, m_dim));
    }
    //Cancello dal cromosoma le città copiate
    for(auto city:new_cities){
      if(city != -1 and city != 0){
        //cout << "Cerco la città (" << city << ") tra quelle vecchie e la cancello" << endl; 
        m_cities.erase( find(m_cities.begin(), m_cities.end(), city) );
      } 
    }

    //Sposto le città rimanenti nello stesso ordine in cui si trovavano
    for(uint i = 1; i < m_cities.size(); i++){
      //cout << "Tra le restanti, prendo la città (" << m_cities.at(i) << ")..." << endl;
      for(uint j = 0; j < new_cities.size(); j++){
        //cout << "...controllo il sito " << j << " del nuovo array con la città " << new_cities.at(j) << "..." << endl;
        if(new_cities.at(j) == -1){
          new_cities.at(j) = m_cities.at(i);
          break;
        }
      }
    }
    m_cities = new_cities;

    #ifdef DEBUG_GEN
    cout << endl << "Dopo:  "; this->print();
    #endif

    generated_correctly = check();
    return *this;
}

chromosome & chromosome::Perm(){
    int M = int( rnd.Rannyu()* (m_dim/2 - 1) ) + 1;
    int start1 = int( rnd.Rannyu()* (m_dim - 1) ) + 1;
    int start2;
    do{
      start2 = int( rnd.Rannyu()* (m_dim - 1) ) + 1;
    }while(start2 < mod(start1 + M, m_dim) );

    #ifdef DEBUG_GEN
    cout << "\tMuto il cromosoma (permutazione)." << endl;
    cout << "Permuto " << M << " città contigue con altre " << M << ". Punti di partenza: " << start1 << " e " << start2 << endl;
    cout << "Prima: "; this->print();
    #endif

    vector<int> tmp;
    for(int i = 0; i < M; i++){
      tmp.push_back( m_cities.at( mod(start1 + i, m_dim) ));
      m_cities.at( mod(start1 + i, m_dim) ) = m_cities.at( mod(start2 + i, m_dim) );
      m_cities.at( mod(start2 + i, m_dim) ) = tmp.at(i);
    }

    #ifdef DEBUG_GEN
    cout << endl << "Dopo:  "; this->print();
    #endif

    generated_correctly = check();
    return *this;
}

chromosome & chromosome::Inverse(){
    int start = rnd.Rannyu() * (m_dim - 1) + 1;
    int M = rnd.Rannyu() * (m_dim - 3) + 2;
    vector<int> tmp;
    
    #ifdef DEBUG_GEN
    cout << "\tMuto il cromosoma (inversione)." << endl;
    cout << "Inveeto l'ordine di " << M << " città contigue, partendo dalla posizione " << start << endl;
    cout << "Prima: "; this->print();
    #endif

    for(int i = start; i < start + M; i++){ 
      #ifdef DEBUG_GEN
      cout << "Copio la città (" << m_cities.at( mod(i, m_dim) ) << ") che sta in posizione " << mod(i, m_dim) << endl;
      #endif
      tmp.push_back(m_cities.at( mod(i,m_dim) ));
    }
    for(int i = 0; i < M; i++){
      #ifdef DEBUG_GEN
      cout << "Reinserisco la città (" << tmp.at(tmp.size() - i - 1) << ") nella posizione " << mod(start + i, m_dim) << endl;
      #endif
      m_cities.at(mod(start + i, m_dim)) = tmp.at(tmp.size() - i - 1);
    }

    #ifdef DEBUG_GEN
    cout << endl << "Dopo:  "; this->print();
    #endif

    generated_correctly = check();
    return *this;
}

chromosome Population::get_best(){
  double min = INT32_MAX;
  double min_pos = 0;
  for(int i = 0; i < m_dim_pop; i++){
    if(m_population.at(i).get_fit() < min){
      min = m_population.at(i).get_fit();
      min_pos = i;
    }
  }
  return m_population.at(min_pos);
}

chromosome::~chromosome()
{
  m_cities.clear();
}

////////////////////////////////////////////////////////////////////////////////////////////////

Population::Population(bool on_circ, int ncity, int dim_pop, double lenght)
{
  //auto rnd = SingletonRand::get_instance().m_rnd;
  #ifdef DEBUG_GEN
  cout << "\nGenero Popolazione..." << endl;
  #endif

  m_dim_pop = dim_pop;
  //Genero le città distribuite casualmente
  if(on_circ){
    //Genero il vettore di città distribuite casualmente su una circonferenza
    for(int i = 0; i < ncity; i++ ){
      double phi = rnd.Rannyu()*M_PI*2;
      double x = lenght*cos(phi);
      double y = lenght*sin(phi);
      city single_city{x,y,i};
      m_cities_list.push_back(single_city);
    }
  }
  else {
    //Genero il vettore di città distribuite casualmente su una circonferenza
    for(int i = 0; i < ncity; i++ ){
      double x = lenght*rnd.Rannyu();
      double y = lenght*rnd.Rannyu();
      city single_city{x,y,i};
      m_cities_list.push_back(single_city);
    }
  }

  //Genero la popolazione di cromosomi iniziale
  for(int i = 0; i < dim_pop; i++){
    #ifdef DEBUG_GEN
    cout << endl << "Genero cromosoma numero " << i << endl;
    #endif
    chromosome indiv(ncity /*, lenght*/ );
    m_population.push_back(indiv);
  }

  #ifdef DEBUG_GEN
  cout << endl << "Ho generato la popolazione" << endl;
  #endif
}

Population::Population(string inname, int dim_pop)
{
  //auto rnd = SingletonRand::get_instance().m_rnd;
  #ifdef DEBUG_GEN
  cout << "\nGenero Popolazione..." << endl;
  #endif

  ifstream input; input.open(inname);
  int tot = 0;
  city single_city;
  string tmp;
  do
  {
    getline(input, tmp, ',');
    cout << "letto " << tmp <<endl;
    single_city.label = stoi(tmp);
    getline(input, tmp, ',');
    cout << "letto " << tmp <<endl;
    single_city.x = stod(tmp);
    getline(input, tmp);
    cout << "letto " << tmp <<endl;
    single_city.y = stod(tmp);
    m_cities_list.push_back(single_city);
    tot ++;
  }while (!input.eof());
  input.close();
  

  m_dim_pop = dim_pop;

  //Genero la popolazione di cromosomi iniziale
  for(int i = 0; i < m_dim_pop; i++){
    #ifdef DEBUG_GEN
    cout << endl << "Genero cromosoma numero " << i << endl;
    #endif
    chromosome indiv(tot);
    m_population.push_back(indiv);
  }

  #ifdef DEBUG_GEN
  cout << endl << "Ho generato la popolazione" << endl;
  #endif
}

Population::~Population(){}

Population & Population::fprint_cities(){
  ofstream cities;
  cities.open("city_list.dat");
  for(auto city:m_cities_list){
    cities << endl << city.label << "," << city.x << "," << city.y;
  }
  cities.close();

  return *this;
}

Population & Population::print_chromosomes(){
  for(int i = 0; i < m_dim_pop; i++){
    //cout << "Stampo su file il cromosoma " << i  << endl;
    cout << i << ": ";
    m_population.at(i).print();
    cout << endl;
  }

  return *this;
}

Population & Population::estimate_fit(){
  #ifdef DEBUG_GEN
  cout << "Inizio stima fit..." << endl ;
  #endif

  for(auto &chrom:m_population){
    double fit = 0;
    for(int i = 1; i < chrom.get_dim(); i++){
      int label = chrom.get_city(i);
      int label2 = chrom.get_city(i-1);
      fit += pow(m_cities_list.at(label).x - m_cities_list.at(label2).x, 2) + pow(m_cities_list.at(label).y - m_cities_list.at(label2).y, 2);
    }
    int label = chrom.get_city(chrom.get_dim()-1);
    fit += pow(m_cities_list.at(0).x - m_cities_list.at(label).x, 2) + pow(m_cities_list.at(0).y - m_cities_list.at(label).y, 2);
    chrom.set_fit(fit);
  }

  #ifdef DEBUG_GEN
  cout << flush << "Fine stima fit" << endl ;
  #endif

  return *this;
}

double Population::estimate_half_av_fit(){
  #ifdef DEBUG_GEN
  cout << "Inizio stima fit..." << endl ;
  #endif

  //Seleziono la metà dei cromosomi migliori
  //Ordino i cromosomi dal migliore al peggiore
  sort(m_population.begin(), m_population.end(), Compare);

  double avfit = 0;
  for(int i = 0; i < m_dim_pop/2; i++){
    avfit += m_population.at(i).get_fit();
    //cout << m_population.at(i).get_fit() << endl;
  }

  #ifdef DEBUG_GEN
  cout << flush << "Fine stima fit" << endl ;
  #endif

  return avfit*2/m_dim_pop;
}

vector<chromosome> Population::select(bool torunament){
  //auto rnd = SingletonRand::get_instance().m_rnd;

  int first_chosen = -1 /* m_dim_pop-1 */, second_chosen = -1 /* m_dim_pop-1 */;
  vector<chromosome> parents;

  //bool torunament = true;

  if(torunament){

    //Vengono selezionati due sottogruppi casuale di cromosomi da cui verranno estratti i due migliori (da ciascuno) come genitori
    //La dimensione del sottogruppo è casualmente posta tra il 5% e il 10% della dimenisione della popolazione
    int group_dim1 = ceil(m_dim_pop*rnd.Rannyu(0, 0.1));
    int group_dim2 = ceil(m_dim_pop*rnd.Rannyu(0, 0.1));
    //Eseguo uno shuffle dei cromosomi, che potrebbero essere stati ordinati dalla funzione che calcola il fit medio
    random_shuffle(m_population.begin(), m_population.end());
    //Genero un numero casuale compreso tra 0 e m_dim_pop - gorup_dim
    int stp1 = int(rnd.Rannyu()*(m_dim_pop - group_dim1));
    int stp2 = int(rnd.Rannyu()*(m_dim_pop - group_dim2));
    //Inserisco i cromosomi dalla posizione stp nel gruppo
    vector<chromosome> group1;
    vector<chromosome> group2;
    for(int i = stp1; i < m_dim_pop; i++){
      group1.push_back(m_population.at(i));
    }
    for(int i = stp2; i < m_dim_pop; i++){
      group2.push_back(m_population.at(i));
    }
    //Ordino il gruppo ponendo all'inizio i migliori
    sort(group1.begin(), group1.end(), Compare);
    sort(group2.begin(), group2.end(), Compare);
    parents.push_back(group1.at(0)); parents.push_back(group2.at(0));

  } else{
    #ifdef DEBUG_GEN
    cout << "\nInizio operazione di selezione...\n" ;
    #endif
    
    double total_fit = 0;
    for(auto i:m_population){
      total_fit += 1./i.get_fit();
      // total_fit += i.get_fit();
    }

    //Viene scelto un cromosoma con probabilità (1/fit) / total_fit.
    bool chosen = false;
    double curr_threshold = 0;

    double thrw = rnd.Rannyu();
    #ifdef DEBUG_GEN
    cout << "\tThrow per primo = " << thrw << endl;
    #endif
    for(int k = 0; k < m_dim_pop and !chosen; k++){
        curr_threshold += 1./m_population.at(k).get_fit()/total_fit;
        // curr_threshold += 1./(m_dim_pop-1)*(1-m_population.at(k).get_fit()/total_fit);
        if(thrw < curr_threshold) {
          chosen = true;
          first_chosen = k;
        } 
    }
    if (!chosen) first_chosen = m_dim_pop - 1;

    //cout << "first chosen = " << first_chosen << endl;
    do{
      curr_threshold = 0;
      chosen = false;
      thrw = rnd.Rannyu();
      #ifdef DEBUG_GEN
      cout << "\tThrow per secondo = " << thrw << endl;
      #endif
      for(int k = 0; k<m_dim_pop and !chosen; k++){
        curr_threshold += 1./m_population.at(k).get_fit()/total_fit;
        // curr_threshold += 1./(m_dim_pop-1)*(1-m_population.at(k).get_fit()/total_fit);
        if(thrw < curr_threshold) {
          chosen = true;
          second_chosen = k;
        } 
      }
    }while(first_chosen == second_chosen);
    //cout << "2nd chosen = " << second_chosen << endl;
    
    #ifdef DEBUG_GEN
    cout << "\tScelti i cromosomi:\n";
    cout << "\t" << first_chosen << ": "; m_population.at(first_chosen).print();
    cout << "\t" << second_chosen << ": "; m_population.at(second_chosen).print();
    #endif

    parents.push_back(m_population.at(first_chosen)); parents.push_back(m_population.at(second_chosen));
  }

  return parents;
  
}

//////////////////////////////////////////////////////////////////////////////////////

bool crossingover(vector<chromosome> &parents){
  //auto rnd = SingletonRand::get_instance().m_rnd;

  //La funzione esegue il crossingover con una probabilità > 50%
  double p_cross = 0.7;
  double thrw = rnd.Rannyu();
  if (thrw < p_cross){
    #ifdef DEBUG_GEN
    cout << endl << "Tentativo di crossingover fallito!\n";
    #endif
    return true;
  }

  //Effettua il crossingover tra il cromosoma a e b, modificandoli direttamente
  //Taglio in una posizione casuale
  chromosome a = parents.at(0);
  chromosome b = parents.at(1);
  int cut = int( rnd.Rannyu()*(a.get_dim()-2))  + 1; //evito di tagliare prima del primo posto e dopo l'ultimo
  #ifdef DEBUG_GEN
  cout << "\nInizio il crossingover tra i cromosomi scelti" << endl;
  cout << "\ta: "; a.print();
  cout << "\tb: "; b.print();
  cout << "\tTaglio dal " << cut << " elemento" << endl;
  #endif

  //Sostituisco le città dopo il taglio con quelle del cromosoma opposto, in ordine di come si presentano
  vector<int> missing_a, missing_b; //città che saranno cambiate dopo il taglio
  for(int i = cut; i < a.get_dim(); i++){
    missing_a.push_back(a.get_city(i));
    missing_b.push_back(b.get_city(i));
  }

  for(int i = cut; i < a.get_dim(); i++){
    for(int j = 0; j < a.get_dim(); j++){
      //Verifico se b.get_city(j) è presente in missing_a
      //  Se non è presente, j++
      //  Se è presente: parents.at(0).set_city(i, b.get_city(j));
      auto e = find(missing_a.begin(), missing_a.end(), b.get_city(j));
      if ( e != missing_a.end() ){
         parents.at(0).set_city(i, b.get_city(j));
         missing_a.erase(e);
         break;
      }
    }
  }

  for(int i = cut; i < a.get_dim(); i++){
    for(int j = 0; j < a.get_dim(); j++){
      //Verifico se a.get_city(j) è presente in missing_b
      //  Se non è presente, j++
      //  Se è presente: parents.at(1).set_city(i, a.get_city(j));
      auto e = find(missing_b.begin(), missing_b.end(), a.get_city(j));
      if ( e != missing_b.end() ){
         parents.at(1).set_city(i, a.get_city(j));
         missing_b.erase(e);
         break;
      }
    }
  }

  #ifdef DEBUG_GEN
  cout << "\tI figli:" << endl;
  cout << "\ta: [";
    for(int j = 0; j < a.get_dim(); j++){
      cout << parents.at(0).get_city(j) << " ";
    }
  cout << "]" << endl;
  cout << "\tb: [";
    for(int j = 0; j < a.get_dim(); j++){
      cout << parents.at(1).get_city(j) << " ";
    }
  cout << "]" << endl;
  #endif

  return parents.at(0).check() and parents.at(1).check();
  
}

bool Compare(chromosome a, chromosome b){
  return a.get_fit() < b.get_fit();
}